const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const stageA = require('../public/stage-A-master-table.json');

const tsContent = fs.readFileSync(path.join(__dirname, '..', 'src', 'data', 'sparc-datasets.ts'), 'utf8');

function parseRCData(ts) {
  const rcMap = {};
  const re = /"([^"]+)":\s*\{[^[]*data:\s*\[([\s\S]*?)\]/g;
  let m;
  while ((m = re.exec(ts)) !== null) {
    const name = m[1];
    const dataStr = m[2];
    const points = [];
    const ptRe = /r:\s*([\d.]+)\s*,\s*v:\s*([\d.]+)/g;
    let pm;
    while ((pm = ptRe.exec(dataStr)) !== null) {
      points.push({ r: parseFloat(pm[1]), v: parseFloat(pm[2]) });
    }
    if (points.length >= 3) rcMap[name] = points;
  }
  return rcMap;
}

const rcMapRaw = parseRCData(tsContent);
console.log(`Parsed RC data for ${Object.keys(rcMapRaw).length} galaxies`);

function normalize(n) {
  return n.replace(/\s+/g,'').replace(/^NGC0*/,'NGC').replace(/^DDO0*/,'DDO').replace(/^IC0*/,'IC').replace(/^UGC0*/,'UGC').replace(/^PGC0*/,'PGC').toUpperCase();
}

const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) {
  rcMap[normalize(name)] = pts;
}

function pearsonR(x, y) {
  const n = x.length;
  if (n < 4) return NaN;
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) {
    const dx = x[i] - mx, dy = y[i] - my;
    num += dx * dy; dx2 += dx * dx; dy2 += dy * dy;
  }
  return dx2 > 0 && dy2 > 0 ? num / Math.sqrt(dx2 * dy2) : 0;
}

function multiR2(X, y) {
  const n = y.length; const nv = X[0].length;
  const my = y.reduce((a, b) => a + b, 0) / n;
  const mx = Array(nv).fill(0);
  for (let j = 0; j < nv; j++) { for (let i = 0; i < n; i++) mx[j] += X[i][j]; mx[j] /= n; }
  const XTX = Array.from({ length: nv }, () => Array(nv).fill(0));
  const XTy = Array(nv).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < nv; j++) {
      XTy[j] += (X[i][j] - mx[j]) * (y[i] - my);
      for (let k = 0; k < nv; k++) XTX[j][k] += (X[i][j] - mx[j]) * (X[i][k] - mx[k]);
    }
  }
  const aug = XTX.map((row, i) => [...row, XTy[i]]);
  for (let col = 0; col < nv; col++) {
    let maxRow = col;
    for (let row = col + 1; row < nv; row++) if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row;
    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
    if (Math.abs(aug[col][col]) < 1e-12) continue;
    for (let row = col + 1; row < nv; row++) {
      const f = aug[row][col] / aug[col][col];
      for (let j = col; j <= nv; j++) aug[row][j] -= f * aug[col][j];
    }
  }
  const beta = Array(nv).fill(0);
  for (let i = nv - 1; i >= 0; i--) {
    beta[i] = aug[i][nv];
    for (let j = i + 1; j < nv; j++) beta[i] -= aug[i][j] * beta[j];
    beta[i] /= aug[i][i] || 1;
  }
  let sse = 0, sst = 0;
  const residuals = [];
  for (let i = 0; i < n; i++) {
    let pred = my;
    for (let j = 0; j < nv; j++) pred += beta[j] * (X[i][j] - mx[j]);
    residuals.push(y[i] - pred);
    sse += (y[i] - pred) ** 2;
    sst += (y[i] - my) ** 2;
  }
  return { R2: sst > 0 ? 1 - sse / sst : 0, beta, residuals };
}

function partialR(x, y, z) {
  if (x.length < 5) return NaN;
  const n = x.length;
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  const mz = z.reduce((a, b) => a + b, 0) / n;
  let sxz = 0, szz = 0, syz = 0;
  for (let i = 0; i < n; i++) { sxz += (x[i]-mx)*(z[i]-mz); szz += (z[i]-mz)**2; syz += (y[i]-my)*(z[i]-mz); }
  const bxz = szz > 0 ? sxz/szz : 0;
  const byz = szz > 0 ? syz/szz : 0;
  const rx = x.map((v,i) => v - mx - bxz*(z[i]-mz));
  const ry = y.map((v,i) => v - my - byz*(z[i]-mz));
  return pearsonR(rx, ry);
}

function partialR_multi(x, y, zArr) {
  const n = x.length;
  if (n < 5) return NaN;
  const resX = multiR2(Array.from({length:n},(_,i) => zArr.map(z=>z[i])), x).residuals;
  const resY = multiR2(Array.from({length:n},(_,i) => zArr.map(z=>z[i])), y).residuals;
  return pearsonR(resX, resY);
}

function interpolateV(points, targetR) {
  if (targetR <= points[0].r) return points[0].v;
  if (targetR >= points[points.length - 1].r) return points[points.length - 1].v;
  for (let i = 0; i < points.length - 1; i++) {
    if (points[i].r <= targetR && points[i + 1].r >= targetR) {
      const frac = (targetR - points[i].r) / (points[i + 1].r - points[i].r);
      return points[i].v + frac * (points[i + 1].v - points[i].v);
    }
  }
  return points[points.length - 1].v;
}


const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const stageAMap = {};
stageA.galaxies.forEach(g => { stageAMap[g.name] = g; stageAMap[normalize(g.name)] = g; });

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  const sa = stageAMap[g.name] || stageAMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;

  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 4) continue;

  const Rdisk = sp.Rdisk;
  const Vflat = sp.Vflat;

  const V_at_Rd = interpolateV(rc, Rdisk);
  const V_at_2Rd = interpolateV(rc, 2.2 * Rdisk);

  const nOuter = Math.min(4, Math.floor(rc.length / 3));
  const outerPts = rc.slice(-nOuter);
  let outerSlope = 0;
  if (outerPts.length >= 2) {
    const dr = outerPts[outerPts.length-1].r - outerPts[0].r;
    const dv = outerPts[outerPts.length-1].v - outerPts[0].v;
    outerSlope = dr > 0 ? dv / dr : 0;
  }

  let R70 = rc[rc.length - 1].r;
  for (let i = 0; i < rc.length; i++) {
    if (rc[i].v >= 0.7 * Vflat) { R70 = rc[i].r; break; }
  }

  let R90 = rc[rc.length - 1].r;
  for (let i = 0; i < rc.length; i++) {
    if (rc[i].v >= 0.9 * Vflat) { R90 = rc[i].r; break; }
  }

  const innerGrad = rc.length >= 3 && rc[1].r > rc[0].r ?
    (rc[1].v - rc[0].v) / (rc[1].r - rc[0].r) : NaN;

  const concIdx = V_at_2Rd / (rc[rc.length - 1].v || 1);

  const logVflat = Math.log10(Vflat);
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(Rdisk, 0.01));
  const logMbar = Math.log10(Math.max(sp.L36 * 0.5 + sp.MHI * 1.33, 0.001) * 1e9);

  gals.push({
    name: g.name,
    logA0: g.logA0,
    Vflat, logVflat, logL36, logRdisk, logMbar,
    logMHI: g.logMHI,
    morphT: sp.T,
    logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    inc: sp.inc, Q: sp.Q,
    envCode: g.envCode,
    rcWig: g.rcWiggliness,
    logSig0: g.logSigma0 || (sa ? sa.logSigma0 : 0),
    logMR: g.logMeanRun || (sa ? sa.logMeanRun : 0),
    Rdisk,

    V_Rd_norm: V_at_Rd / Vflat,
    V_2Rd_norm: V_at_2Rd / Vflat,
    outerSlope: outerSlope,
    outerSlopeNorm: outerSlope / (Vflat / (rc[rc.length-1].r || 1)),
    R70_norm: R70 / Rdisk,
    R90_norm: R90 / Rdisk,
    innerGrad: innerGrad,
    innerGradNorm: isFinite(innerGrad) ? innerGrad / (Vflat / Rdisk) : NaN,
    concIdx: concIdx,

    nRC: rc.length,
    Rmax: rc[rc.length - 1].r,
    Rmax_norm: rc[rc.length - 1].r / Rdisk,
  });
}

const vfModel = multiR2(
  gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]),
  gals.map(g => g.logVflat)
);
for (let i = 0; i < gals.length; i++) gals[i].VfResid = vfModel.residuals[i];

const baseStructVars = [g=>g.logMbar, g=>g.logL36, g=>g.logRdisk, g=>g.morphT, g=>g.logMHI, g=>g.logSBdisk];
const a0Model = multiR2(
  gals.map(g => baseStructVars.map(f => f(g))),
  gals.map(g => g.logA0)
);
for (let i = 0; i < gals.length; i++) gals[i].a0Resid = a0Model.residuals[i];

const valid = gals.filter(g =>
  isFinite(g.VfResid) && isFinite(g.a0Resid) && isFinite(g.V_Rd_norm) &&
  isFinite(g.R70_norm) && isFinite(g.outerSlope) && isFinite(g.innerGradNorm)
);

const rawR = pearsonR(valid.map(g => g.VfResid), valid.map(g => g.logA0));
const residResidR = pearsonR(valid.map(g => g.VfResid), valid.map(g => g.a0Resid));

console.log('═'.repeat(70));
console.log('PHASE 406: RC-SHAPE VARIABLE HUNT');
console.log('Does any rotation curve shape parameter explain the hidden channel?');
console.log('═'.repeat(70));
console.log(`\nWorking sample: N=${valid.length} (with RC data + all controls)`);
console.log(`Raw r(VfResid, logA0) = ${rawR.toFixed(3)}`);
console.log(`Residual r(VfResid, a₀_resid) = ${residResidR.toFixed(3)}`);


const rcVars = [
  { name: 'V_Rd_norm', extract: g => g.V_Rd_norm, desc: 'V(Rdisk)/Vflat — inner RC shape' },
  { name: 'V_2Rd_norm', extract: g => g.V_2Rd_norm, desc: 'V(2.2Rd)/Vflat — rise completeness' },
  { name: 'innerGradNorm', extract: g => g.innerGradNorm, desc: 'Normalized inner gradient' },
  { name: 'concIdx', extract: g => g.concIdx, desc: 'V(2.2Rd)/V(Rlast) — RC concentration' },
  { name: 'R70_norm', extract: g => g.R70_norm, desc: 'R(0.7Vflat)/Rdisk — rise radius' },
  { name: 'R90_norm', extract: g => g.R90_norm, desc: 'R(0.9Vflat)/Rdisk — flat onset' },
  { name: 'outerSlope', extract: g => g.outerSlope, desc: 'Outer RC slope (km/s/kpc)' },
  { name: 'outerSlopeNorm', extract: g => g.outerSlopeNorm, desc: 'Outer slope normalized' },
  { name: 'Rmax_norm', extract: g => g.Rmax_norm, desc: 'RC extent / Rdisk' },
];


console.log('\n' + '▓'.repeat(70));
console.log('406A: SINGLE RC-SHAPE VARIABLES VS THE CHANNEL');
console.log('▓'.repeat(70));

console.log('\n--- Correlations with VfResid, a₀, and a₀_resid ---\n');
console.log('  Variable          r(.,VfR)  r(.,a₀)  r(.,a0res)  r(.,Vfl)  taut?  desc');
console.log('  ' + '─'.repeat(85));

const singleResults = [];
for (const rv of rcVars) {
  const vf = valid.filter(g => isFinite(rv.extract(g)));
  if (vf.length < 10) continue;
  const x = vf.map(g => rv.extract(g));
  const rVfR = pearsonR(x, vf.map(g => g.VfResid));
  const rA0 = pearsonR(x, vf.map(g => g.logA0));
  const rA0res = pearsonR(x, vf.map(g => g.a0Resid));
  const rVfl = pearsonR(x, vf.map(g => g.logVflat));

  const taut = Math.abs(rVfl) > 0.7 ? 'HIGH' : (Math.abs(rVfl) > 0.4 ? 'MED' : 'LOW');

  console.log(`  ${rv.name.padEnd(18)} ${(rVfR>=0?'+':'')+rVfR.toFixed(3)}   ${(rA0>=0?'+':'')+rA0.toFixed(3)}   ${(rA0res>=0?'+':'')+rA0res.toFixed(3)}     ${(rVfl>=0?'+':'')+rVfl.toFixed(3)}   ${taut.padEnd(4)}   ${rv.desc}`);
  singleResults.push({ name: rv.name, rVfResid: rVfR, rA0, rA0resid: rA0res, rVflat: rVfl, tautRisk: taut });
}


console.log('\n\n' + '▓'.repeat(70));
console.log('406B: ABSORPTION TEST — DOES ANY RC-SHAPE ABSORB THE CHANNEL?');
console.log('▓'.repeat(70));

console.log('\n--- Partial r(VfResid, a₀ | RC-shape variable) ---\n');

console.log('  Variable          raw r    partial r   change   absorbs?');
console.log('  ' + '─'.repeat(60));

const absorptionResults = [];
for (const rv of rcVars) {
  const vf = valid.filter(g => isFinite(rv.extract(g)));
  if (vf.length < 10) continue;
  const pr = partialR(vf.map(g => g.VfResid), vf.map(g => g.logA0), vf.map(g => rv.extract(g)));
  const change = pr - rawR;
  const absorbs = pr < rawR - 0.10;
  console.log(`  ${rv.name.padEnd(18)} ${rawR.toFixed(3)}    ${(pr>=0?'+':'')+pr.toFixed(3)}      ${(change>=0?'+':'')+change.toFixed(3)}    ${absorbs ? 'YES' : 'no'}`);
  absorptionResults.push({ name: rv.name, rawR, partialR: pr, change, absorbs });
}

console.log('\n--- Partial r(VfResid, a₀_RESID | RC-shape variable) ---');
console.log('    This is the CLEAN test: both sides purged of structure\n');

console.log('  Variable          resid r  partial r   change   absorbs?');
console.log('  ' + '─'.repeat(60));

const cleanAbsorption = [];
for (const rv of rcVars) {
  const vf = valid.filter(g => isFinite(rv.extract(g)));
  if (vf.length < 10) continue;
  const pr = partialR(vf.map(g => g.VfResid), vf.map(g => g.a0Resid), vf.map(g => rv.extract(g)));
  const change = pr - residResidR;
  const absorbs = pr < residResidR - 0.10;
  console.log(`  ${rv.name.padEnd(18)} ${residResidR.toFixed(3)}    ${(pr>=0?'+':'')+pr.toFixed(3)}      ${(change>=0?'+':'')+change.toFixed(3)}    ${absorbs ? 'YES' : 'no'}`);
  cleanAbsorption.push({ name: rv.name, residR: residResidR, partialR: pr, change, absorbs });
}


console.log('\n\n' + '▓'.repeat(70));
console.log('406C: MULTI-VARIABLE RC-SHAPE MODEL');
console.log('Can ALL RC-shape vars together absorb the channel?');
console.log('▓'.repeat(70));

const lowTautVars = rcVars.filter(rv => {
  const rVfl = pearsonR(valid.map(g => rv.extract(g)), valid.map(g => g.logVflat));
  return Math.abs(rVfl) < 0.5 && valid.every(g => isFinite(rv.extract(g)));
});

console.log(`\n  Low-tautology RC vars (|r(.,Vflat)| < 0.5): ${lowTautVars.map(v=>v.name).join(', ')}`);

if (lowTautVars.length >= 2) {
  const pr_lowTaut = partialR_multi(
    valid.map(g => g.VfResid), valid.map(g => g.logA0),
    lowTautVars.map(rv => valid.map(g => rv.extract(g)))
  );
  console.log(`  partial r(VfResid, a₀ | low-taut RC vars) = ${pr_lowTaut.toFixed(3)} (raw ${rawR.toFixed(3)})`);

  const pr_lowTaut_clean = partialR_multi(
    valid.map(g => g.VfResid), valid.map(g => g.a0Resid),
    lowTautVars.map(rv => valid.map(g => rv.extract(g)))
  );
  console.log(`  partial r(VfResid, a₀_resid | low-taut RC vars) = ${pr_lowTaut_clean.toFixed(3)} (residual ${residResidR.toFixed(3)})`);
}

const allRCValid = rcVars.filter(rv => valid.every(g => isFinite(rv.extract(g))));
if (allRCValid.length >= 3) {
  const pr_allRC = partialR_multi(
    valid.map(g => g.VfResid), valid.map(g => g.logA0),
    allRCValid.map(rv => valid.map(g => rv.extract(g)))
  );
  console.log(`\n  partial r(VfResid, a₀ | ALL ${allRCValid.length} RC vars) = ${pr_allRC.toFixed(3)} (raw ${rawR.toFixed(3)})`);

  const pr_allRC_clean = partialR_multi(
    valid.map(g => g.VfResid), valid.map(g => g.a0Resid),
    allRCValid.map(rv => valid.map(g => rv.extract(g)))
  );
  console.log(`  partial r(VfResid, a₀_resid | ALL ${allRCValid.length} RC vars) = ${pr_allRC_clean.toFixed(3)} (residual ${residResidR.toFixed(3)})`);
}

console.log('\n--- RC-shape + structural controls combined ---');
const structAndRC = [
  ...baseStructVars.map(f => valid.map((g,i) => f(g))),
  ...allRCValid.map(rv => valid.map(g => rv.extract(g)))
];
const pr_everything = partialR_multi(
  valid.map(g => g.VfResid), valid.map(g => g.logA0),
  structAndRC
);
console.log(`  partial r(VfResid, a₀ | 6 struct + ${allRCValid.length} RC) = ${pr_everything.toFixed(3)}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('406D: WHAT DO RC-SHAPE VARS EXPLAIN ABOUT VfResid?');
console.log('▓'.repeat(70));

const rcModel = multiR2(
  valid.map(g => allRCValid.map(rv => rv.extract(g))),
  valid.map(g => g.VfResid)
);
console.log(`\n  R²(VfResid ~ all RC-shape vars) = ${rcModel.R2.toFixed(3)}`);

const structModel = multiR2(
  valid.map(g => baseStructVars.map(f => f(g))),
  valid.map(g => g.VfResid)
);
console.log(`  R²(VfResid ~ 6 structural vars) = ${structModel.R2.toFixed(3)}`);

const bothModel = multiR2(
  valid.map(g => [...baseStructVars.map(f => f(g)), ...allRCValid.map(rv => rv.extract(g))]),
  valid.map(g => g.VfResid)
);
console.log(`  R²(VfResid ~ struct + RC-shape) = ${bothModel.R2.toFixed(3)}`);
console.log(`  RC-shape adds: +${(bothModel.R2 - structModel.R2).toFixed(3)} above structure`);


console.log('\n\n' + '▓'.repeat(70));
console.log('406E: TAUTOLOGY CHECK — ARE RC-SHAPE VARS JUST Vflat IN DISGUISE?');
console.log('▓'.repeat(70));

console.log('\n  r(RC-var, logVflat) and r(RC-var, logRdisk):\n');
for (const rv of rcVars) {
  const rV = pearsonR(valid.map(g => rv.extract(g)), valid.map(g => g.logVflat));
  const rR = pearsonR(valid.map(g => rv.extract(g)), valid.map(g => g.logRdisk));
  const rM = pearsonR(valid.map(g => rv.extract(g)), valid.map(g => g.logMbar));
  console.log(`  ${rv.name.padEnd(18)} r(Vfl)=${(rV>=0?'+':'')+rV.toFixed(3)}  r(Rd)=${(rR>=0?'+':'')+rR.toFixed(3)}  r(Mbar)=${(rM>=0?'+':'')+rM.toFixed(3)}`);
}

const rcModelVflat = multiR2(
  valid.map(g => [g.logVflat, g.logRdisk]),
  valid.map(g => g.V_Rd_norm)
);
console.log(`\n  R²(V_Rd_norm ~ logVflat + logRdisk) = ${rcModelVflat.R2.toFixed(3)}`);
console.log(`  → V_Rd_norm is ${rcModelVflat.R2 > 0.9 ? 'TAUTOLOGICAL' : (rcModelVflat.R2 > 0.5 ? 'PARTIALLY redundant' : 'INDEPENDENT')} of VfResid building blocks`);


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 406 GRAND SYNTHESIS');
console.log('═'.repeat(70));

const bestAbsorber = absorptionResults.reduce((a, b) => a.change < b.change ? a : b);
const bestClean = cleanAbsorption.reduce((a, b) => a.change < b.change ? a : b);

console.log(`
Key findings:

1. Best single absorber (raw): ${bestAbsorber.name}
   partial r = ${bestAbsorber.partialR.toFixed(3)} (change ${bestAbsorber.change.toFixed(3)})
   ${bestAbsorber.absorbs ? '→ ABSORBS part of the channel' : '→ Does NOT absorb the channel'}

2. Best single absorber (clean residual): ${bestClean.name}
   partial r = ${bestClean.partialR.toFixed(3)} (change ${bestClean.change.toFixed(3)})
   ${bestClean.absorbs ? '→ ABSORBS part of the channel' : '→ Does NOT absorb the channel'}

3. RC-shape vars explain R²=${rcModel.R2.toFixed(3)} of VfResid variance
   Adding RC-shape to structure: R²=${bothModel.R2.toFixed(3)} (+${(bothModel.R2-structModel.R2).toFixed(3)})

4. After ALL RC-shape + structure controls:
   partial r = ${pr_everything.toFixed(3)} (raw ${rawR.toFixed(3)})
`);

const anyAbsorbs = absorptionResults.some(r => r.absorbs) || cleanAbsorption.some(r => r.absorbs);

if (!anyAbsorbs) {
  console.log('VERDICT: NO RC-shape variable absorbs the VfResid–a₀ channel.');
  console.log('The hidden component is NOT encoded in the rotation curve shape.');
  console.log('It operates at a level BELOW what RC morphology can capture.');
} else {
  const absorbers = [...new Set([
    ...absorptionResults.filter(r => r.absorbs).map(r => r.name),
    ...cleanAbsorption.filter(r => r.absorbs).map(r => r.name),
  ])];
  console.log(`VERDICT: ${absorbers.join(', ')} partially absorb the channel.`);
  console.log('These RC-shape variables may carry part of the hidden information.');
}

const outPath = path.join(__dirname, '..', 'public', 'phase406-rc-shape-hunt.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '406',
  title: 'RC-Shape Variable Hunt',
  timestamp: new Date().toISOString(),
  N: valid.length,
  rawR, residResidR,
  singleCorrelations: singleResults,
  absorptionRaw: absorptionResults,
  absorptionClean: cleanAbsorption,
  rcModelR2: rcModel.R2,
  structModelR2: structModel.R2,
  bothModelR2: bothModel.R2,
  partialR_everything: pr_everything,
}, null, 2));

console.log(`\nSaved: ${outPath}`);
