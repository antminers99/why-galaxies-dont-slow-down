const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const stageA = require('../public/stage-A-master-table.json');
const s11 = require('../public/phase11-sensitivity-lab.json');
const d300 = require('../public/phase300-sample-salvage.json');

function pearsonR(x, y) {
  const n = x.length;
  if (n < 3) return NaN;
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) {
    const dx = x[i] - mx, dy = y[i] - my;
    num += dx * dy; dx2 += dx * dx; dy2 += dy * dy;
  }
  return dx2 > 0 && dy2 > 0 ? num / Math.sqrt(dx2 * dy2) : 0;
}

function spearmanR(x, y) {
  const n = x.length;
  if (n < 3) return NaN;
  function rank(arr) {
    const sorted = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
    const r = Array(n);
    for (let i = 0; i < n; i++) r[sorted[i].i] = i + 1;
    return r;
  }
  return pearsonR(rank(x), rank(y));
}

function linReg(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let num = 0, den = 0;
  for (let i = 0; i < n; i++) { num += (x[i] - mx) * (y[i] - my); den += (x[i] - mx) ** 2; }
  const slope = den > 0 ? num / den : 0;
  const intercept = my - slope * mx;
  let sse = 0;
  for (let i = 0; i < n; i++) sse += (y[i] - slope * x[i] - intercept) ** 2;
  const se = n > 2 ? Math.sqrt(sse / (n - 2) / (den || 1)) : 0;
  const R2 = 1 - sse / (y.reduce((a, v) => a + (v - my) ** 2, 0) || 1);
  return { slope, intercept, R2, se, n };
}

function solveLinear(A, b) {
  const n = b.length;
  const aug = A.map((row, i) => [...row, b[i]]);
  for (let col = 0; col < n; col++) {
    let maxRow = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row;
    }
    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
    if (Math.abs(aug[col][col]) < 1e-12) continue;
    for (let row = col + 1; row < n; row++) {
      const f = aug[row][col] / aug[col][col];
      for (let j = col; j <= n; j++) aug[row][j] -= f * aug[col][j];
    }
  }
  const x = Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = aug[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= aug[i][j] * x[j];
    x[i] /= aug[i][i] || 1;
  }
  return x;
}

function multiR2(X, y) {
  const n = y.length;
  const nv = X[0].length;
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
  const beta = solveLinear(XTX, XTy);
  let sse = 0, sst = 0;
  for (let i = 0; i < n; i++) {
    let pred = my;
    for (let j = 0; j < nv; j++) pred += beta[j] * (X[i][j] - mx[j]);
    sse += (y[i] - pred) ** 2;
    sst += (y[i] - my) ** 2;
  }
  return { R2: sst > 0 ? 1 - sse / sst : 0, beta, residuals: y.map((yi, i) => {
    let pred = my;
    for (let j = 0; j < nv; j++) pred += beta[j] * (X[i][j] - mx[j]);
    return yi - pred;
  })};
}


const sparcMap = {};
sparc.forEach(g => sparcMap[g.name] = g);

const stageAMap = {};
stageA.galaxies.forEach(g => stageAMap[g.name] = g);

const s11Map = {};
s11.galaxies.forEach(g => s11Map[g.name] = g);

const internal = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name];
  const sa = stageAMap[g.name];
  const ss = s11Map[g.name];
  if (!sp || sp.Vflat <= 0) continue;

  internal.push({
    name: g.name,
    logA0: g.logA0,
    Vflat: sp.Vflat,
    logVflat: Math.log10(sp.Vflat),
    L36: sp.L36,
    logL36: Math.log10(Math.max(sp.L36, 0.001)),
    MHI: sp.MHI,
    logMHI: g.logMHI,
    Rdisk: sp.Rdisk,
    logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)),
    T: sp.T,
    morphT: sp.T,
    SBdisk: sp.SBdisk,
    logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    SBeff: sp.SBeff,
    Reff: sp.Reff,
    inc: sp.inc,
    Q: sp.Q,
    eVflat: sp.eVflat,
    rcWig: g.rcWiggliness,
    envCode: g.envCode,
    logSig0: g.logSigma0 || (sa ? sa.logSigma0 : 0),
    logMR: g.logMeanRun || (sa ? sa.logMeanRun : 0),
    logMbar: Math.log10(Math.max(sp.L36 * 0.5 + sp.MHI * 1.33, 0.001) * 1e9),
    ncmAmp: sa ? sa.things_ncm_amp : null,
    ncmFrac: sa ? sa.things_ncm_frac : null,
    lopsidedness: sa ? sa.things_lopsidedness : null,
    bisymFlow: sa ? sa.things_bisymFlow : null,
    hiDef: sa ? sa.hi_deficiency : null,
    nPts: ss ? ss.n : null,
    Vmax: ss ? ss.Vmax : sp.Vflat,
    sample: 'internal'
  });
}

const external = [];
if (d300.augmentedGalaxies) {
  for (const g of d300.augmentedGalaxies) {
    if (!g.Vflat || g.Vflat <= 0) continue;
    const sp = sparcMap[g.name];
    external.push({
      name: g.name,
      logA0: g.logA0,
      Vflat: g.Vflat,
      logVflat: g.logVflat,
      L36: sp ? sp.L36 : Math.pow(10, g.logL36),
      logL36: g.logL36,
      logMHI: g.logMHI,
      Rdisk: sp ? sp.Rdisk : Math.pow(10, g.logRdisk),
      logRdisk: g.logRdisk,
      morphT: g.morphT,
      T: g.morphT,
      logMbar: g.logMbar,
      VfResid: g.VfResid,
      lhOuter: g.lhOuter,
      haloK: g.haloK,
      Q: g.Q,
      inc: g.inc,
      nRARpts: g.nRARpts,
      fitRMS: g.fitRMS,
      logMR: g.logMeanRun,
      sample: 'external'
    });
  }
}

console.log('='.repeat(70));
console.log('PHASE 403: ANATOMY OF VfResid');
console.log('='.repeat(70));
console.log(`\nInternal sample: N=${internal.length}`);
console.log(`External sample: N=${external.length}`);


console.log('\n' + '▓'.repeat(70));
console.log('PART A: CROSS-SAMPLE REPLICATION OF REGIME LAW');
console.log('    Is the strengthening real or a SPARC artifact?');
console.log('▓'.repeat(70));

function computeVfResid(gals) {
  const valid = gals.filter(g => g.logVflat && g.logMbar && g.logL36 && g.logRdisk);
  if (valid.length < 10) return;
  const X = valid.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT || 5]);
  const y = valid.map(g => g.logVflat);
  const result = multiR2(X, y);
  for (let i = 0; i < valid.length; i++) {
    valid[i].VfResid = result.residuals[i];
  }
}

computeVfResid(internal);

function regimeAnalysis(label, gals) {
  const valid = gals.filter(g => g.VfResid !== undefined && g.logA0 && isFinite(g.logA0));
  if (valid.length < 8) {
    console.log(`  ${label}: N=${valid.length} (insufficient)`);
    return null;
  }

  const rAll = pearsonR(valid.map(g => g.VfResid), valid.map(g => g.logA0));
  const sAll = spearmanR(valid.map(g => g.VfResid), valid.map(g => g.logA0));

  const lowV = valid.filter(g => g.Vflat < 120);
  const highV = valid.filter(g => g.Vflat >= 120);
  const vhighV = valid.filter(g => g.Vflat >= 180);

  const rLow = lowV.length >= 5 ? pearsonR(lowV.map(g => g.VfResid), lowV.map(g => g.logA0)) : null;
  const rHigh = highV.length >= 5 ? pearsonR(highV.map(g => g.VfResid), highV.map(g => g.logA0)) : null;
  const rVHigh = vhighV.length >= 5 ? pearsonR(vhighV.map(g => g.VfResid), vhighV.map(g => g.logA0)) : null;

  const delta = rLow !== null && rHigh !== null ? rHigh - rLow : null;
  const strengthens = delta !== null && delta > 0.05;

  console.log(`  ${label}: N=${valid.length}`);
  console.log(`    r(all)=${rAll.toFixed(3)}, rho(all)=${sAll.toFixed(3)}`);
  console.log(`    low-V (N=${lowV.length}): ${rLow !== null ? rLow.toFixed(3) : 'N/A'}`);
  console.log(`    high-V (N=${highV.length}): ${rHigh !== null ? rHigh.toFixed(3) : 'N/A'}`);
  console.log(`    v-high-V (N=${vhighV.length}): ${rVHigh !== null ? rVHigh.toFixed(3) : 'N/A'}`);
  console.log(`    delta = ${delta !== null ? (delta > 0 ? '+' : '') + delta.toFixed(3) : 'N/A'} → ${strengthens ? 'STRENGTHENS' : 'weakens/flat'}`);

  return { label, N: valid.length, rAll, sAll, rLow, rHigh, rVHigh, delta, strengthens,
    nLow: lowV.length, nHigh: highV.length, nVHigh: vhighV.length };
}

console.log('\n--- Internal (training) sample ---');
const intRegime = regimeAnalysis('Internal', internal);

console.log('\n--- External (validation) sample ---');
const extRegime = regimeAnalysis('External', external);

const pooled = [...internal.filter(g => g.VfResid !== undefined), ...external.filter(g => g.VfResid !== undefined)];
console.log('\n--- Pooled sample ---');
const poolRegime = regimeAnalysis('Pooled', pooled);

const highQ = internal.filter(g => g.Q === 1);
computeVfResid(highQ);
console.log('\n--- Internal Q=1 only (highest quality) ---');
const q1Regime = regimeAnalysis('Q=1', highQ);

const lowInc = internal.filter(g => g.inc > 40 && g.inc < 80);
computeVfResid(lowInc);
console.log('\n--- Moderate inclination (40°–80°) ---');
const incRegime = regimeAnalysis('ModInc', lowInc);


console.log('\n' + '▓'.repeat(70));
console.log('PART B: WHAT GENERATES VfResid?');
console.log('    Decompose VfResid into its physical drivers');
console.log('▓'.repeat(70));

const validInt = internal.filter(g => g.VfResid !== undefined);

const candidates = [
  { name: 'logMHI', extract: g => g.logMHI },
  { name: 'rcWiggliness', extract: g => g.rcWig },
  { name: 'envCode', extract: g => g.envCode },
  { name: 'logSigma0', extract: g => g.logSig0 },
  { name: 'logMeanRun', extract: g => g.logMR },
  { name: 'logSBdisk', extract: g => g.logSBdisk },
  { name: 'morphT', extract: g => g.morphT },
  { name: 'inc', extract: g => g.inc },
  { name: 'Q', extract: g => g.Q },
  { name: 'logReff', extract: g => Math.log10(Math.max(g.Reff || 0.01, 0.01)) },
  { name: 'logL36', extract: g => g.logL36 },
  { name: 'logMbar', extract: g => g.logMbar },
  { name: 'logA0', extract: g => g.logA0 },
];

console.log('\n--- Single predictors of VfResid ---');
console.log('    Which observable EXPLAINS where VfResid comes from?\n');

const singleResults = [];
for (const c of candidates) {
  const pairs = validInt.filter(g => {
    const v = c.extract(g);
    return v !== null && v !== undefined && isFinite(v);
  });
  if (pairs.length < 10) continue;
  const x = pairs.map(g => c.extract(g));
  const y = pairs.map(g => g.VfResid);
  const r = pearsonR(x, y);
  const rho = spearmanR(x, y);
  const reg = linReg(x, y);
  singleResults.push({ name: c.name, r, rho, R2: reg.R2, slope: reg.slope, N: pairs.length });
}

singleResults.sort((a, b) => Math.abs(b.r) - Math.abs(a.r));
console.log('  Rank  Variable       r       ρ       R²     slope    N');
console.log('  ' + '─'.repeat(60));
singleResults.forEach((s, i) => {
  console.log(`  ${String(i+1).padStart(2)}.   ${s.name.padEnd(14)} ${s.r >= 0 ? '+' : ''}${s.r.toFixed(3)}   ${s.rho >= 0 ? '+' : ''}${s.rho.toFixed(3)}   ${s.R2.toFixed(3)}   ${s.slope >= 0 ? '+' : ''}${s.slope.toFixed(3)}   ${s.N}`);
});

console.log('\n--- Multi-predictor models of VfResid ---');
console.log('    How much of VfResid can structural variables explain?\n');

const structVars = ['logMHI', 'rcWig', 'envCode', 'logSig0', 'logMR'];
const vfResidValid = validInt.filter(g => {
  return g.logMHI !== undefined && g.rcWig !== undefined && g.envCode !== undefined
    && g.logSig0 !== undefined && g.logMR !== undefined && isFinite(g.logSig0) && isFinite(g.logMR);
});

if (vfResidValid.length >= 15) {
  const models = [
    { name: 'logMHI only', vars: [g => g.logMHI] },
    { name: 'logMHI + rcWig', vars: [g => g.logMHI, g => g.rcWig] },
    { name: 'logMHI + rcWig + envCode', vars: [g => g.logMHI, g => g.rcWig, g => g.envCode] },
    { name: '5-var (MHI,rcWig,env,Sig0,MR)', vars: [g => g.logMHI, g => g.rcWig, g => g.envCode, g => g.logSig0, g => g.logMR] },
    { name: '5-var + logA0', vars: [g => g.logMHI, g => g.rcWig, g => g.envCode, g => g.logSig0, g => g.logMR, g => g.logA0] },
    { name: 'logA0 only', vars: [g => g.logA0] },
  ];

  for (const m of models) {
    const X = vfResidValid.map(g => m.vars.map(f => f(g)));
    const y = vfResidValid.map(g => g.VfResid);
    const result = multiR2(X, y);
    console.log(`  ${m.name.padEnd(40)} R²=${result.R2.toFixed(3)} (N=${vfResidValid.length})`);
  }
}

console.log('\n--- VfResid vs a₀ after controlling for structure ---');
const structModel = multiR2(
  vfResidValid.map(g => [g.logMHI, g.rcWig, g.envCode, g.logSig0, g.logMR]),
  vfResidValid.map(g => g.VfResid)
);
const structResid = structModel.residuals;
const a0vals = vfResidValid.map(g => g.logA0);
const rStructResid = pearsonR(structResid, a0vals);
console.log(`  r(VfResid_struct_residual, logA0) = ${rStructResid.toFixed(3)}`);
console.log(`  This is the VfResid-a₀ coupling AFTER removing all known structural drivers.`);


console.log('\n' + '▓'.repeat(70));
console.log('PART C: WHY HIGH-V ONLY? — Regime Anatomy');
console.log('    Map the activation pattern in detail');
console.log('▓'.repeat(70));

const validForRegime = internal.filter(g => g.VfResid !== undefined && isFinite(g.logA0));

const bins = [
  { label: 'V<80', min: 0, max: 80 },
  { label: '80-120', min: 80, max: 120 },
  { label: '120-160', min: 120, max: 160 },
  { label: '160-200', min: 160, max: 200 },
  { label: '200-250', min: 200, max: 250 },
  { label: 'V>250', min: 250, max: 999 },
];

console.log('\n--- r(VfResid, logA0) by Vflat bin ---\n');
console.log('  Bin          N     r(VfResid,a₀)   rho');
console.log('  ' + '─'.repeat(50));

const binResults = [];
for (const bin of bins) {
  const inBin = validForRegime.filter(g => g.Vflat >= bin.min && g.Vflat < bin.max);
  if (inBin.length < 4) {
    console.log(`  ${bin.label.padEnd(12)} ${String(inBin.length).padStart(3)}     —`);
    binResults.push({ bin: bin.label, N: inBin.length, r: null, rho: null });
    continue;
  }
  const r = pearsonR(inBin.map(g => g.VfResid), inBin.map(g => g.logA0));
  const rho = spearmanR(inBin.map(g => g.VfResid), inBin.map(g => g.logA0));
  console.log(`  ${bin.label.padEnd(12)} ${String(inBin.length).padStart(3)}     ${r >= 0 ? '+' : ''}${r.toFixed(3)}           ${rho >= 0 ? '+' : ''}${rho.toFixed(3)}`);
  binResults.push({ bin: bin.label, N: inBin.length, r, rho });
}

console.log('\n--- Running r(VfResid, logA0) with increasing Vflat threshold ---\n');
console.log('  Vflat≥    N     r(VfResid,a₀)   rho     R²(a₀→VfResid)');
console.log('  ' + '─'.repeat(60));

const thresholds = [0, 50, 80, 100, 120, 140, 160, 180, 200, 220, 250];
const threshResults = [];
for (const t of thresholds) {
  const above = validForRegime.filter(g => g.Vflat >= t);
  if (above.length < 5) continue;
  const r = pearsonR(above.map(g => g.VfResid), above.map(g => g.logA0));
  const rho = spearmanR(above.map(g => g.VfResid), above.map(g => g.logA0));
  const reg = linReg(above.map(g => g.logA0), above.map(g => g.VfResid));
  console.log(`  ≥${String(t).padStart(3)}     ${String(above.length).padStart(3)}     ${r >= 0 ? '+' : ''}${r.toFixed(3)}           ${rho >= 0 ? '+' : ''}${rho.toFixed(3)}     ${reg.R2.toFixed(3)}`);
  threshResults.push({ threshold: t, N: above.length, r, rho, R2: reg.R2 });
}

console.log('\n--- What changes at high-V? Property distributions ---\n');

const lowV = validForRegime.filter(g => g.Vflat < 120);
const highV = validForRegime.filter(g => g.Vflat >= 120);

function meanStd(arr) {
  const mu = arr.reduce((a, b) => a + b, 0) / arr.length;
  const std = Math.sqrt(arr.reduce((a, v) => a + (v - mu) ** 2, 0) / arr.length);
  return { mu, std };
}

const props = [
  { name: 'logMHI', extract: g => g.logMHI },
  { name: 'rcWiggliness', extract: g => g.rcWig },
  { name: 'envCode', extract: g => g.envCode },
  { name: 'logSigma0', extract: g => g.logSig0 },
  { name: 'logMeanRun', extract: g => g.logMR },
  { name: 'logA0', extract: g => g.logA0 },
  { name: 'VfResid', extract: g => g.VfResid },
  { name: 'inc', extract: g => g.inc },
  { name: 'Q', extract: g => g.Q },
  { name: 'nPts', extract: g => g.nPts },
  { name: 'logSBdisk', extract: g => g.logSBdisk },
];

console.log('  Property        low-V(N=' + lowV.length + ')           high-V(N=' + highV.length + ')');
console.log('  ' + '─'.repeat(60));
for (const p of props) {
  const lv = lowV.map(g => p.extract(g)).filter(v => v !== null && v !== undefined && isFinite(v));
  const hv = highV.map(g => p.extract(g)).filter(v => v !== null && v !== undefined && isFinite(v));
  if (lv.length < 3 || hv.length < 3) continue;
  const ls = meanStd(lv);
  const hs = meanStd(hv);
  console.log(`  ${p.name.padEnd(16)} ${ls.mu.toFixed(3)}±${ls.std.toFixed(3)}    ${hs.mu.toFixed(3)}±${hs.std.toFixed(3)}`);
}


console.log('\n' + '▓'.repeat(70));
console.log('PART D: VfResid-a₀ COUPLING CONDITIONAL ON OTHER VARIABLES');
console.log('    Does coupling survive after controlling for each variable?');
console.log('▓'.repeat(70));

console.log('\n--- Partial r(VfResid, a₀ | X) for each control variable ---\n');

function partialR(x, y, z) {
  if (x.length < 5) return NaN;
  const regXZ = linReg(z, x);
  const regYZ = linReg(z, y);
  const resX = x.map((v, i) => v - regXZ.slope * z[i] - regXZ.intercept);
  const resY = y.map((v, i) => v - regYZ.slope * z[i] - regYZ.intercept);
  return pearsonR(resX, resY);
}

const vfResidArr = validForRegime.map(g => g.VfResid);
const logA0Arr = validForRegime.map(g => g.logA0);

const controls = [
  { name: 'logMHI', vals: validForRegime.map(g => g.logMHI) },
  { name: 'rcWiggliness', vals: validForRegime.map(g => g.rcWig) },
  { name: 'envCode', vals: validForRegime.map(g => g.envCode) },
  { name: 'logSigma0', vals: validForRegime.map(g => g.logSig0) },
  { name: 'logMeanRun', vals: validForRegime.map(g => g.logMR) },
  { name: 'logVflat', vals: validForRegime.map(g => g.logVflat) },
  { name: 'morphT', vals: validForRegime.map(g => g.morphT) },
  { name: 'inc', vals: validForRegime.map(g => g.inc) },
  { name: 'logSBdisk', vals: validForRegime.map(g => g.logSBdisk) },
];

const rawR = pearsonR(vfResidArr, logA0Arr);
console.log(`  Raw r(VfResid, logA0) = ${rawR.toFixed(3)}\n`);
console.log('  Control variable    partial r    change');
console.log('  ' + '─'.repeat(45));

const partialResults = [];
for (const ctrl of controls) {
  const valid = ctrl.vals.every(v => v !== null && v !== undefined && isFinite(v));
  if (!valid) {
    console.log(`  ${ctrl.name.padEnd(20)} N/A (missing data)`);
    continue;
  }
  const pr = partialR(vfResidArr, logA0Arr, ctrl.vals);
  const change = pr - rawR;
  console.log(`  ${ctrl.name.padEnd(20)} ${pr >= 0 ? '+' : ''}${pr.toFixed(3)}       ${change >= 0 ? '+' : ''}${change.toFixed(3)}`);
  partialResults.push({ control: ctrl.name, partialR: pr, change });
}


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 403 SYNTHESIS');
console.log('═'.repeat(70));

console.log('\n1. CROSS-SAMPLE REPLICATION:');
if (intRegime && extRegime) {
  console.log(`   Internal: delta=${intRegime.delta !== null ? (intRegime.delta > 0 ? '+' : '') + intRegime.delta.toFixed(3) : 'N/A'} (${intRegime.strengthens ? 'STRENGTHENS' : 'flat/weakens'})`);
  console.log(`   External: delta=${extRegime.delta !== null ? (extRegime.delta > 0 ? '+' : '') + extRegime.delta.toFixed(3) : 'N/A'} (${extRegime.strengthens ? 'STRENGTHENS' : 'flat/weakens'})`);
  if (intRegime.strengthens && extRegime.strengthens) {
    console.log('   → REPLICATED: Regime strengthening seen in BOTH samples');
  } else if (intRegime.strengthens || extRegime.strengthens) {
    console.log('   → PARTIAL: Strengthening in one sample only');
  } else {
    console.log('   → NOTE: Pattern may be weaker than expected in split samples');
  }
}

console.log('\n2. VfResid DRIVERS:');
console.log('   Top single predictors: ' + singleResults.slice(0, 3).map(s => `${s.name}(r=${s.r.toFixed(2)})`).join(', '));

console.log('\n3. HIGH-V ACTIVATION:');
const thr160 = threshResults.find(t => t.threshold === 160);
const thrAll = threshResults.find(t => t.threshold === 0);
if (thr160 && thrAll) {
  console.log(`   Full sample: r=${thrAll.r.toFixed(3)}`);
  console.log(`   V≥160 only: r=${thr160.r.toFixed(3)}`);
  console.log(`   Activation ratio: ${(thr160.r / (thrAll.r || 0.001)).toFixed(2)}x`);
}

console.log('\n4. PARTIAL CORRELATIONS:');
console.log('   VfResid-a₀ coupling after structural control: r=' + rStructResid.toFixed(3));

const outPath = path.join(__dirname, '..', 'public', 'phase403-vfresid-anatomy.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '403',
  title: 'Anatomy of VfResid',
  timestamp: new Date().toISOString(),
  partA: {
    internal: intRegime,
    external: extRegime,
    pooled: poolRegime,
    q1: q1Regime,
    modInc: incRegime
  },
  partB: {
    singlePredictors: singleResults,
    structResidualVsA0: rStructResid
  },
  partC: {
    binned: binResults,
    thresholds: threshResults,
  },
  partD: {
    rawR: rawR,
    partials: partialResults
  }
}, null, 2));

console.log(`\nSaved: ${outPath}`);
