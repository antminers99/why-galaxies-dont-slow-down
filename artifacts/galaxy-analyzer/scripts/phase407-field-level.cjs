const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');
const stageA = require('../public/stage-A-master-table.json');

const G = 4.3009e-6;

const tsContent = fs.readFileSync(path.join(__dirname, '..', 'src', 'data', 'sparc-datasets.ts'), 'utf8');
function parseRCData(ts) {
  const rcMap = {};
  const re = /"([^"]+)":\s*\{[^[]*data:\s*\[([\s\S]*?)\]/g;
  let m;
  while ((m = re.exec(ts)) !== null) {
    const points = [];
    const ptRe = /r:\s*([\d.]+)\s*,\s*v:\s*([\d.]+)/g;
    let pm;
    while ((pm = ptRe.exec(m[2])) !== null) {
      points.push({ r: parseFloat(pm[1]), v: parseFloat(pm[2]) });
    }
    if (points.length >= 3) rcMap[m[1]] = points;
  }
  return rcMap;
}

function normalize(n) {
  return n.replace(/\s+/g,'').replace(/^NGC0*/,'NGC').replace(/^DDO0*/,'DDO').replace(/^IC0*/,'IC').replace(/^UGC0*/,'UGC').replace(/^PGC0*/,'PGC').toUpperCase();
}

const rcMapRaw = parseRCData(tsContent);
const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) rcMap[normalize(name)] = pts;

const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[g.name] = g; resultsMap[normalize(g.name)] = g; });

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

const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const stageAMap = {};
stageA.galaxies.forEach(g => { stageAMap[g.name] = g; stageAMap[normalize(g.name)] = g; });
const baseStructVars = [g=>g.logMbar, g=>g.logL36, g=>g.logRdisk, g=>g.morphT, g=>g.logMHI, g=>g.logSBdisk];

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  const sa = stageAMap[g.name] || stageAMap[normalize(g.name)];
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0 || !sr) continue;

  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 5) continue;

  const Vflat = sp.Vflat;
  const Rdisk = sp.Rdisk;
  const logVflat = Math.log10(Vflat);
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(Rdisk, 0.01));
  const logMbar = Math.log10(Math.max(sp.L36 * 0.5 + sp.MHI * 1.33, 0.001) * 1e9);

  const M_newt = sr.models.newtonian.M;
  const M_halo = sr.models.dark_halo_linear.M;
  const k_halo = sr.models.dark_halo_linear.k;

  const profile = rc.map(pt => {
    const V_obs = pt.v;
    const V_bar = Math.sqrt(G * M_newt / pt.r);
    const V_halo_model = Math.sqrt(Math.max(0, G * M_halo / pt.r + k_halo * pt.r));

    const g_obs = V_obs * V_obs / pt.r;
    const g_bar = V_bar * V_bar / pt.r;

    const accRatio = g_bar > 0 ? g_obs / g_bar : 0;
    const V_excess = V_obs - V_bar;
    const V2_excess = V_obs * V_obs - V_bar * V_bar;

    return {
      r: pt.r, V_obs, V_bar, V_halo_model,
      g_obs, g_bar, accRatio,
      V_excess, V2_excess,
      r_norm: pt.r / Rdisk,
    };
  });

  const accRatios = profile.map(p => p.accRatio);
  const meanAccRatio = accRatios.reduce((a,b) => a+b, 0) / accRatios.length;
  const logMeanAccRatio = Math.log10(Math.max(meanAccRatio, 0.01));

  const nHalf = Math.floor(profile.length / 2);
  const innerAccRatio = profile.slice(0, nHalf).reduce((a,p) => a + p.accRatio, 0) / nHalf;
  const outerAccRatio = profile.slice(nHalf).reduce((a,p) => a + p.accRatio, 0) / (profile.length - nHalf);
  const accGradient = outerAccRatio - innerAccRatio;

  let cumMismatch = 0;
  for (let i = 1; i < profile.length; i++) {
    const dr = profile[i].r - profile[i-1].r;
    const avgExcess = (Math.abs(profile[i].V2_excess) + Math.abs(profile[i-1].V2_excess)) / 2;
    cumMismatch += avgExcess * dr;
  }
  const totalExtent = profile[profile.length-1].r - profile[0].r;
  const normMismatch = totalExtent > 0 ? cumMismatch / (Vflat * Vflat * totalExtent) : 0;

  let R_trans = profile[profile.length-1].r;
  for (const p of profile) {
    if (p.accRatio > 1.5) { R_trans = p.r; break; }
  }
  const R_trans_norm = R_trans / Rdisk;

  let transStart = -1, transEnd = -1;
  for (let i = 0; i < profile.length; i++) {
    if (profile[i].accRatio > 1.2 && transStart < 0) transStart = profile[i].r;
    if (profile[i].accRatio > 3.0 && transEnd < 0) transEnd = profile[i].r;
  }
  const transWidth = (transStart > 0 && transEnd > 0) ? (transEnd - transStart) / Rdisk : NaN;

  const residProfile = profile.map(p => p.V_excess / Vflat);
  let coherenceSum = 0, coherenceCount = 0;
  for (let i = 1; i < residProfile.length; i++) {
    coherenceSum += residProfile[i] * residProfile[i-1];
    coherenceCount++;
  }
  const radialCoherence = coherenceCount > 0 ? coherenceSum / coherenceCount : 0;

  const logAccProfile = profile.map(p => Math.log10(Math.max(p.g_obs, 1e-15)));
  const logBarProfile = profile.map(p => Math.log10(Math.max(p.g_bar, 1e-15)));
  const rProfile = profile.map(p => p.r);

  let slopeAccMismatch = NaN;
  if (profile.length >= 4) {
    const logRatios = profile.map(p => Math.log10(Math.max(p.accRatio, 0.01)));
    const logR = profile.map(p => Math.log10(Math.max(p.r, 0.01)));
    const n = logR.length;
    const mX = logR.reduce((a,b)=>a+b,0)/n;
    const mY = logRatios.reduce((a,b)=>a+b,0)/n;
    let num=0, den=0;
    for(let i=0;i<n;i++){num+=(logR[i]-mX)*(logRatios[i]-mY);den+=(logR[i]-mX)**2;}
    slopeAccMismatch = den > 0 ? num/den : 0;
  }

  const darkDominance = k_halo > 0 && M_halo > 0 ? Math.sqrt(G * M_halo / k_halo) : NaN;
  const darkDom_norm = isFinite(darkDominance) ? darkDominance / Rdisk : NaN;

  const haloStrength = k_halo > 0 ? Math.log10(k_halo) : NaN;

  const barFrac_Rd = profile.length > 0 ?
    (() => {
      let closest = profile[0];
      for (const p of profile) if (Math.abs(p.r - Rdisk) < Math.abs(closest.r - Rdisk)) closest = p;
      return closest.g_bar > 0 ? closest.g_bar / closest.g_obs : NaN;
    })() : NaN;

  const barFrac_2Rd = profile.length > 0 ?
    (() => {
      let closest = profile[0];
      for (const p of profile) if (Math.abs(p.r - 2*Rdisk) < Math.abs(closest.r - 2*Rdisk)) closest = p;
      return closest.g_bar > 0 ? closest.g_bar / closest.g_obs : NaN;
    })() : NaN;

  gals.push({
    name: g.name, logA0: g.logA0,
    Vflat, logVflat, logL36, logRdisk, logMbar,
    logMHI: g.logMHI, morphT: sp.T,
    logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),

    logMeanAccRatio,
    accGradient,
    normMismatch,
    R_trans_norm,
    transWidth,
    radialCoherence,
    slopeAccMismatch,
    darkDom_norm,
    haloStrength,
    barFrac_Rd,
    barFrac_2Rd,
    innerAccRatio: Math.log10(Math.max(innerAccRatio, 0.01)),
    outerAccRatio: Math.log10(Math.max(outerAccRatio, 0.01)),
  });
}

const vfModel = multiR2(
  gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]),
  gals.map(g => g.logVflat)
);
for (let i = 0; i < gals.length; i++) gals[i].VfResid = vfModel.residuals[i];

const a0Model = multiR2(
  gals.map(g => baseStructVars.map(f => f(g))),
  gals.map(g => g.logA0)
);
for (let i = 0; i < gals.length; i++) gals[i].a0Resid = a0Model.residuals[i];

const valid = gals.filter(g =>
  isFinite(g.VfResid) && isFinite(g.a0Resid) && isFinite(g.logMeanAccRatio) &&
  isFinite(g.normMismatch) && isFinite(g.slopeAccMismatch)
);

const rawR = pearsonR(valid.map(g => g.VfResid), valid.map(g => g.logA0));
const residResidR = pearsonR(valid.map(g => g.VfResid), valid.map(g => g.a0Resid));

console.log('═'.repeat(70));
console.log('PHASE 407: FIELD-LEVEL / RADIAL COUPLING VARIABLES');
console.log('Does the radial baryon-dynamics mismatch explain the channel?');
console.log('═'.repeat(70));
console.log(`\nWorking sample: N=${valid.length}`);
console.log(`Raw r(VfResid, logA0) = ${rawR.toFixed(3)}`);
console.log(`Residual r(VfResid, a₀_resid) = ${residResidR.toFixed(3)}`);

const fieldVars = [
  { name: 'logMeanAccRatio', extract: g => g.logMeanAccRatio, desc: 'log mean g_obs/g_bar across RC' },
  { name: 'accGradient', extract: g => g.accGradient, desc: 'outer - inner acceleration ratio' },
  { name: 'normMismatch', extract: g => g.normMismatch, desc: 'integrated V² mismatch / (Vfl²·Rext)' },
  { name: 'R_trans_norm', extract: g => g.R_trans_norm, desc: 'transition R (accRatio>1.5) / Rdisk' },
  { name: 'transWidth', extract: g => g.transWidth, desc: 'transition zone width / Rdisk' },
  { name: 'radialCoherence', extract: g => g.radialCoherence, desc: 'autocorr of V_excess profile' },
  { name: 'slopeAccMismatch', extract: g => g.slopeAccMismatch, desc: 'd(logAccRatio)/d(logR) slope' },
  { name: 'darkDom_norm', extract: g => g.darkDom_norm, desc: 'dark dominance radius / Rdisk' },
  { name: 'haloStrength', extract: g => g.haloStrength, desc: 'log(k_halo) from linear model' },
  { name: 'barFrac_Rd', extract: g => g.barFrac_Rd, desc: 'g_bar/g_obs at R=Rdisk' },
  { name: 'barFrac_2Rd', extract: g => g.barFrac_2Rd, desc: 'g_bar/g_obs at R=2Rdisk' },
  { name: 'innerAccRatio', extract: g => g.innerAccRatio, desc: 'log mean accRatio inner half' },
  { name: 'outerAccRatio', extract: g => g.outerAccRatio, desc: 'log mean accRatio outer half' },
];


console.log('\n' + '▓'.repeat(70));
console.log('407A: FIELD-LEVEL CORRELATIONS');
console.log('▓'.repeat(70));

console.log('\n  Variable          r(.,VfR)  r(.,a₀)  r(.,a0res)  r(.,Vfl)  taut?');
console.log('  ' + '─'.repeat(75));

const singleResults = [];
for (const fv of fieldVars) {
  const vf = valid.filter(g => isFinite(fv.extract(g)));
  if (vf.length < 10) { console.log(`  ${fv.name.padEnd(18)} SKIP (N=${vf.length})`); continue; }
  const x = vf.map(g => fv.extract(g));
  const rVfR = pearsonR(x, vf.map(g => g.VfResid));
  const rA0 = pearsonR(x, vf.map(g => g.logA0));
  const rA0res = pearsonR(x, vf.map(g => g.a0Resid));
  const rVfl = pearsonR(x, vf.map(g => g.logVflat));
  const taut = Math.abs(rVfl) > 0.7 ? 'HIGH' : (Math.abs(rVfl) > 0.4 ? 'MED' : 'LOW');
  console.log(`  ${fv.name.padEnd(18)} ${(rVfR>=0?'+':'')+rVfR.toFixed(3)}   ${(rA0>=0?'+':'')+rA0.toFixed(3)}   ${(rA0res>=0?'+':'')+rA0res.toFixed(3)}     ${(rVfl>=0?'+':'')+rVfl.toFixed(3)}   ${taut}`);
  singleResults.push({ name: fv.name, rVfResid: rVfR, rA0, rA0resid: rA0res, rVflat: rVfl, tautRisk: taut, desc: fv.desc });
}


console.log('\n\n' + '▓'.repeat(70));
console.log('407B: ABSORPTION TEST — DOES ANY FIELD VAR ABSORB THE CHANNEL?');
console.log('▓'.repeat(70));

console.log('\n--- Raw: partial r(VfResid, a₀ | field var) ---\n');
console.log('  Variable          raw r    partial r   change   absorbs?');
console.log('  ' + '─'.repeat(60));

const absorptionRaw = [];
for (const fv of fieldVars) {
  const vf = valid.filter(g => isFinite(fv.extract(g)));
  if (vf.length < 10) continue;
  const pr = partialR(vf.map(g=>g.VfResid), vf.map(g=>g.logA0), vf.map(g=>fv.extract(g)));
  const change = pr - rawR;
  const absorbs = pr < rawR - 0.10;
  console.log(`  ${fv.name.padEnd(18)} ${rawR.toFixed(3)}    ${(pr>=0?'+':'')+pr.toFixed(3)}      ${(change>=0?'+':'')+change.toFixed(3)}    ${absorbs?'YES':'no'}`);
  absorptionRaw.push({ name: fv.name, rawR, partialR: pr, change, absorbs });
}

console.log('\n--- Clean: partial r(VfResid, a₀_RESID | field var) ---\n');
console.log('  Variable          resid r  partial r   change   absorbs?');
console.log('  ' + '─'.repeat(60));

const absorptionClean = [];
for (const fv of fieldVars) {
  const vf = valid.filter(g => isFinite(fv.extract(g)));
  if (vf.length < 10) continue;
  const pr = partialR(vf.map(g=>g.VfResid), vf.map(g=>g.a0Resid), vf.map(g=>fv.extract(g)));
  const change = pr - residResidR;
  const absorbs = pr < residResidR - 0.10;
  console.log(`  ${fv.name.padEnd(18)} ${residResidR.toFixed(3)}    ${(pr>=0?'+':'')+pr.toFixed(3)}      ${(change>=0?'+':'')+change.toFixed(3)}    ${absorbs?'YES':'no'}`);
  absorptionClean.push({ name: fv.name, residR: residResidR, partialR: pr, change, absorbs });
}


console.log('\n\n' + '▓'.repeat(70));
console.log('407C: MULTI-VARIABLE FIELD-LEVEL MODEL');
console.log('▓'.repeat(70));

const lowTautField = fieldVars.filter(fv => {
  const vf = valid.filter(g => isFinite(fv.extract(g)));
  if (vf.length < 10) return false;
  const rVfl = pearsonR(vf.map(g=>fv.extract(g)), vf.map(g=>g.logVflat));
  return Math.abs(rVfl) < 0.5 && valid.every(g => isFinite(fv.extract(g)));
});

console.log(`\n  Low-taut field vars: ${lowTautField.map(v=>v.name).join(', ')}`);

if (lowTautField.length >= 2) {
  const pr = partialR_multi(
    valid.map(g=>g.VfResid), valid.map(g=>g.logA0),
    lowTautField.map(fv=>valid.map(g=>fv.extract(g)))
  );
  console.log(`  partial r(VfResid, a₀ | low-taut field) = ${pr.toFixed(3)} (raw ${rawR.toFixed(3)})`);

  const pr2 = partialR_multi(
    valid.map(g=>g.VfResid), valid.map(g=>g.a0Resid),
    lowTautField.map(fv=>valid.map(g=>fv.extract(g)))
  );
  console.log(`  partial r(VfResid, a₀_resid | low-taut field) = ${pr2.toFixed(3)} (residual ${residResidR.toFixed(3)})`);
}

const allFieldValid = fieldVars.filter(fv => valid.every(g => isFinite(fv.extract(g))));
if (allFieldValid.length >= 3) {
  const pr_all = partialR_multi(
    valid.map(g=>g.VfResid), valid.map(g=>g.logA0),
    allFieldValid.map(fv=>valid.map(g=>fv.extract(g)))
  );
  console.log(`\n  partial r(VfResid, a₀ | ALL ${allFieldValid.length} field vars) = ${pr_all.toFixed(3)}`);

  const pr_all2 = partialR_multi(
    valid.map(g=>g.VfResid), valid.map(g=>g.a0Resid),
    allFieldValid.map(fv=>valid.map(g=>fv.extract(g)))
  );
  console.log(`  partial r(VfResid, a₀_resid | ALL field vars) = ${pr_all2.toFixed(3)}`);
}


console.log('\n\n' + '▓'.repeat(70));
console.log('407D: FIELD VARS vs VfResid EXPLANATORY POWER');
console.log('▓'.repeat(70));

if (allFieldValid.length >= 2) {
  const fieldModel = multiR2(
    valid.map(g => allFieldValid.map(fv=>fv.extract(g))),
    valid.map(g => g.VfResid)
  );
  console.log(`\n  R²(VfResid ~ all field vars) = ${fieldModel.R2.toFixed(3)}`);

  const structModel = multiR2(
    valid.map(g => baseStructVars.map(f=>f(g))),
    valid.map(g => g.VfResid)
  );
  console.log(`  R²(VfResid ~ 6 struct vars) = ${structModel.R2.toFixed(3)}`);

  const bothModel = multiR2(
    valid.map(g => [...baseStructVars.map(f=>f(g)), ...allFieldValid.map(fv=>fv.extract(g))]),
    valid.map(g => g.VfResid)
  );
  console.log(`  R²(VfResid ~ struct + field) = ${bothModel.R2.toFixed(3)}`);
  console.log(`  Field adds: +${(bothModel.R2 - structModel.R2).toFixed(3)} above structure`);
}


console.log('\n\n' + '▓'.repeat(70));
console.log('407E: TAUTOLOGY CHECK');
console.log('▓'.repeat(70));
console.log('\n  r(field-var, logVflat), r(field-var, logRdisk), r(field-var, logA0):');
for (const fv of fieldVars) {
  const vf = valid.filter(g => isFinite(fv.extract(g)));
  if (vf.length < 10) continue;
  const x = vf.map(g => fv.extract(g));
  const rV = pearsonR(x, vf.map(g=>g.logVflat));
  const rR = pearsonR(x, vf.map(g=>g.logRdisk));
  const rA = pearsonR(x, vf.map(g=>g.logA0));
  console.log(`  ${fv.name.padEnd(18)} Vfl=${(rV>=0?'+':'')+rV.toFixed(3)}  Rd=${(rR>=0?'+':'')+rR.toFixed(3)}  a₀=${(rA>=0?'+':'')+rA.toFixed(3)}`);
}


console.log('\n\n' + '▓'.repeat(70));
console.log('407F: COMBINED STRUCT + FIELD + a₀ — FINAL PARTIAL r');
console.log('▓'.repeat(70));

const structAndField = [
  ...baseStructVars.map(f => valid.map(g => f(g))),
  ...allFieldValid.map(fv => valid.map(g => fv.extract(g)))
];
const pr_everything = partialR_multi(
  valid.map(g=>g.VfResid), valid.map(g=>g.logA0),
  structAndField
);
console.log(`\n  partial r(VfResid, a₀ | 6 struct + ${allFieldValid.length} field) = ${pr_everything.toFixed(3)}`);
console.log(`  (raw = ${rawR.toFixed(3)})`);


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 407 GRAND SYNTHESIS');
console.log('═'.repeat(70));

const bestRaw = absorptionRaw.length ? absorptionRaw.reduce((a,b) => a.change<b.change?a:b) : null;
const bestClean = absorptionClean.length ? absorptionClean.reduce((a,b) => a.change<b.change?a:b) : null;
const anyAbsorbs = (absorptionRaw.some(r=>r.absorbs) || absorptionClean.some(r=>r.absorbs));

if (bestRaw) console.log(`\n  Best raw absorber: ${bestRaw.name} (Δ=${bestRaw.change.toFixed(3)})`);
if (bestClean) console.log(`  Best clean absorber: ${bestClean.name} (Δ=${bestClean.change.toFixed(3)})`);

if (!anyAbsorbs) {
  console.log('\n  VERDICT: NO field-level variable absorbs the channel.');
  console.log('  The hidden coupling is IRREDUCIBLE to baryon-dynamics mismatch profiles.');
} else {
  const absorbers = [...new Set([
    ...absorptionRaw.filter(r=>r.absorbs).map(r=>r.name),
    ...absorptionClean.filter(r=>r.absorbs).map(r=>r.name),
  ])];
  console.log(`\n  VERDICT: ${absorbers.join(', ')} partially absorb the channel.`);
}

const outPath = path.join(__dirname, '..', 'public', 'phase407-field-level.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '407',
  title: 'Field-Level / Radial Coupling Variables',
  timestamp: new Date().toISOString(),
  N: valid.length,
  rawR, residResidR,
  singleCorrelations: singleResults,
  absorptionRaw,
  absorptionClean,
  partialR_everything: pr_everything,
}, null, 2));

console.log(`\nSaved: ${outPath}`);
