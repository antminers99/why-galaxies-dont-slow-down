const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const stageA = require('../public/stage-A-master-table.json');
const d300 = require('../public/phase300-sample-salvage.json');

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

function spearmanR(x, y) {
  const n = x.length;
  if (n < 4) return NaN;
  function rank(arr) {
    const s = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
    const r = Array(n);
    for (let i = 0; i < n; i++) r[s[i].i] = i + 1;
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
  const sst = y.reduce((a, v) => a + (v - my) ** 2, 0);
  return { slope, intercept, R2: sst > 0 ? 1 - sse / sst : 0, se: n > 2 ? Math.sqrt(sse / (n - 2) / (den || 1)) : 0, n };
}

function partialR(x, y, z) {
  if (x.length < 5) return NaN;
  const regXZ = linReg(z, x);
  const regYZ = linReg(z, y);
  const resX = x.map((v, i) => v - regXZ.slope * z[i] - regXZ.intercept);
  const resY = y.map((v, i) => v - regYZ.slope * z[i] - regYZ.intercept);
  return pearsonR(resX, resY);
}

function partialR2(x, y, zArr) {
  if (x.length < 5) return NaN;
  const n = x.length;
  const nz = zArr.length;
  function residualize(target, predictors) {
    const pn = predictors.length;
    const mt = target.reduce((a, b) => a + b, 0) / n;
    const mp = predictors.map(p => p.reduce((a, b) => a + b, 0) / n);
    const XTX = Array.from({ length: pn }, () => Array(pn).fill(0));
    const XTy = Array(pn).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < pn; j++) {
        XTy[j] += (predictors[j][i] - mp[j]) * (target[i] - mt);
        for (let k = 0; k < pn; k++) XTX[j][k] += (predictors[j][i] - mp[j]) * (predictors[k][i] - mp[k]);
      }
    }
    const aug = XTX.map((row, i) => [...row, XTy[i]]);
    for (let col = 0; col < pn; col++) {
      let maxRow = col;
      for (let row = col + 1; row < pn; row++) if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row;
      [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
      if (Math.abs(aug[col][col]) < 1e-12) continue;
      for (let row = col + 1; row < pn; row++) {
        const f = aug[row][col] / aug[col][col];
        for (let j = col; j <= pn; j++) aug[row][j] -= f * aug[col][j];
      }
    }
    const beta = Array(pn).fill(0);
    for (let i = pn - 1; i >= 0; i--) {
      beta[i] = aug[i][pn];
      for (let j = i + 1; j < pn; j++) beta[i] -= aug[i][j] * beta[j];
      beta[i] /= aug[i][i] || 1;
    }
    return target.map((v, i) => {
      let pred = mt;
      for (let j = 0; j < pn; j++) pred += beta[j] * (predictors[j][i] - mp[j]);
      return v - pred;
    });
  }
  const resX = residualize(x, zArr);
  const resY = residualize(y, zArr);
  return pearsonR(resX, resY);
}

function meanStd(arr) {
  if (arr.length === 0) return { mu: NaN, std: NaN };
  const mu = arr.reduce((a, b) => a + b, 0) / arr.length;
  const std = Math.sqrt(arr.reduce((a, v) => a + (v - mu) ** 2, 0) / arr.length);
  return { mu, std };
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
    const r = y[i] - pred;
    residuals.push(r);
    sse += r * r;
    sst += (y[i] - my) ** 2;
  }
  return { R2: sst > 0 ? 1 - sse / sst : 0, beta, residuals };
}


const sparcMap = {};
sparc.forEach(g => sparcMap[g.name] = g);
const stageAMap = {};
stageA.galaxies.forEach(g => stageAMap[g.name] = g);

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name];
  const sa = stageAMap[g.name];
  if (!sp || sp.Vflat <= 0) continue;

  const logMbar = Math.log10(Math.max(sp.L36 * 0.5 + sp.MHI * 1.33, 0.001) * 1e9);
  const logVflat = Math.log10(sp.Vflat);
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(sp.Rdisk, 0.01));

  gals.push({
    name: g.name,
    logA0: g.logA0,
    Vflat: sp.Vflat,
    logVflat,
    L36: sp.L36,
    logL36,
    MHI: sp.MHI,
    logMHI: g.logMHI,
    Rdisk: sp.Rdisk,
    logRdisk,
    morphT: sp.T,
    logMbar,
    SBdisk: sp.SBdisk,
    logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    SBeff: sp.SBeff,
    Reff: sp.Reff,
    inc: sp.inc,
    Q: sp.Q,
    eVflat: sp.eVflat,
    RHI: sp.RHI,

    rcWig: g.rcWiggliness,
    envCode: g.envCode,
    logSig0: g.logSigma0 || (sa ? sa.logSigma0 : 0),
    logMR: g.logMeanRun || (sa ? sa.logMeanRun : 0),

    ncmAmp: sa ? sa.things_ncm_amp : null,
    ncmFrac: sa ? sa.things_ncm_frac : null,
    lopsidedness: sa ? sa.things_lopsidedness : null,
    bisymFlow: sa ? sa.things_bisymFlow : null,
    c1: sa ? sa.things_c1 : null,
    s1: sa ? sa.things_s1 : null,
    c3: sa ? sa.things_c3 : null,
    s3: sa ? sa.things_s3 : null,
    hiDef: sa ? sa.hi_deficiency : null,
    inTHINGS: sa ? sa.in_THINGS : false,

    logGbarEff: Math.log10(Math.max(sp.Vflat ** 2 / (sp.Rdisk * 3.086e19) * 1e10, 1e-12)),
    logSurfDens: logMbar - 2 * logRdisk,
    dynTime: sp.Rdisk / sp.Vflat * 3.086e16 / 3.156e16,
    logDynTime: Math.log10(Math.max(sp.Rdisk / sp.Vflat * 3.086e16 / 3.156e16, 1e-3)),
    logAccScale: 2 * logVflat - logRdisk,
  });
}

const vfResidModel = multiR2(
  gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]),
  gals.map(g => g.logVflat)
);
for (let i = 0; i < gals.length; i++) {
  gals[i].VfResid = vfResidModel.residuals[i];
}

const valid = gals.filter(g => isFinite(g.VfResid) && isFinite(g.logA0) && isFinite(g.logSig0) && isFinite(g.logMR));

console.log('═'.repeat(70));
console.log('PHASE 404: SCENARIO DISCRIMINATION');
console.log('Which physical origin explains the irreducible VfResid–a₀ channel?');
console.log('═'.repeat(70));
console.log(`\nWorking sample: N=${valid.length}`);

const rawR = pearsonR(valid.map(g => g.VfResid), valid.map(g => g.logA0));
console.log(`Raw r(VfResid, logA0) = ${rawR.toFixed(3)}`);

const pr_Vflat = partialR(valid.map(g => g.VfResid), valid.map(g => g.logA0), valid.map(g => g.logVflat));
console.log(`Partial r(VfResid, a₀ | Vflat) = ${pr_Vflat.toFixed(3)}`);

const structResid = multiR2(
  valid.map(g => [g.logMHI, g.rcWig, g.envCode, g.logSig0, g.logMR]),
  valid.map(g => g.VfResid)
);
const rStructResid = pearsonR(structResid.residuals, valid.map(g => g.logA0));
console.log(`r(VfResid_struct_resid, a₀) = ${rStructResid.toFixed(3)}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('404B: DYNAMICAL COHERENCE HYPOTHESIS');
console.log('Does RC regularity / dynamical maturity explain the channel?');
console.log('▓'.repeat(70));

console.log('\n--- B1: Coherence metrics vs VfResid and a₀ ---\n');

const cohMetrics = [
  { name: 'logMeanRun', extract: g => g.logMR, desc: 'RC regularity' },
  { name: 'logSigma0', extract: g => g.logSig0, desc: 'RC scatter (inverse)' },
  { name: 'rcWiggliness', extract: g => g.rcWig, desc: 'RC wiggliness' },
  { name: 'logDynTime', extract: g => g.logDynTime, desc: 'Dynamical time (Rdisk/Vflat)' },
  { name: 'logAccScale', extract: g => g.logAccScale, desc: '2logV - logR (depth proxy)' },
  { name: 'logSurfDens', extract: g => g.logSurfDens, desc: 'Surface mass density' },
];

console.log('  Metric            r(.,VfRes)  r(.,a₀)   r(.,Vflat)  desc');
console.log('  ' + '─'.repeat(70));

const cohResults = [];
for (const m of cohMetrics) {
  const vals = valid.map(g => m.extract(g)).filter(v => isFinite(v));
  if (vals.length < valid.length * 0.8) continue;
  const vf = valid.filter(g => isFinite(m.extract(g)));
  const x = vf.map(g => m.extract(g));
  const rVfR = pearsonR(x, vf.map(g => g.VfResid));
  const rA0 = pearsonR(x, vf.map(g => g.logA0));
  const rVfl = pearsonR(x, vf.map(g => g.logVflat));
  console.log(`  ${m.name.padEnd(18)} ${(rVfR>=0?'+':'')+rVfR.toFixed(3)}     ${(rA0>=0?'+':'')+rA0.toFixed(3)}     ${(rVfl>=0?'+':'')+rVfl.toFixed(3)}     ${m.desc}`);
  cohResults.push({ name: m.name, rVfResid: rVfR, rA0, rVflat: rVfl });
}

console.log('\n--- B2: Does coherence explain the VfResid–a₀ coupling? ---');
console.log('    partial r(VfResid, a₀ | coherence metrics)\n');

const pr_MR = partialR(valid.map(g => g.VfResid), valid.map(g => g.logA0), valid.map(g => g.logMR));
const pr_Sig = partialR(valid.map(g => g.VfResid), valid.map(g => g.logA0), valid.map(g => g.logSig0));
const pr_Wig = partialR(valid.map(g => g.VfResid), valid.map(g => g.logA0), valid.map(g => g.rcWig));
const pr_DynT = partialR(valid.map(g => g.VfResid), valid.map(g => g.logA0), valid.map(g => g.logDynTime));

const pr_allCoh = partialR2(
  valid.map(g => g.VfResid),
  valid.map(g => g.logA0),
  [valid.map(g => g.logMR), valid.map(g => g.logSig0), valid.map(g => g.rcWig), valid.map(g => g.logDynTime)]
);

console.log(`  Raw:                                  r = ${rawR.toFixed(3)}`);
console.log(`  partial r | logMeanRun:               r = ${pr_MR.toFixed(3)} (change ${(pr_MR-rawR>=0?'+':'')+(pr_MR-rawR).toFixed(3)})`);
console.log(`  partial r | logSigma0:                r = ${pr_Sig.toFixed(3)} (change ${(pr_Sig-rawR>=0?'+':'')+(pr_Sig-rawR).toFixed(3)})`);
console.log(`  partial r | rcWiggliness:             r = ${pr_Wig.toFixed(3)} (change ${(pr_Wig-rawR>=0?'+':'')+(pr_Wig-rawR).toFixed(3)})`);
console.log(`  partial r | logDynTime:               r = ${pr_DynT.toFixed(3)} (change ${(pr_DynT-rawR>=0?'+':'')+(pr_DynT-rawR).toFixed(3)})`);
console.log(`  partial r | ALL 4 coherence:          r = ${pr_allCoh.toFixed(3)} (change ${(pr_allCoh-rawR>=0?'+':'')+(pr_allCoh-rawR).toFixed(3)})`);

console.log('\n--- B3: Does coherence explain the ACTIVATION RAMP? ---');
console.log('    If H4 is right, high-V should have higher coherence, explaining coupling\n');

const lowV = valid.filter(g => g.Vflat < 120);
const highV = valid.filter(g => g.Vflat >= 120);

console.log('  Metric            low-V              high-V             high/low diff');
console.log('  ' + '─'.repeat(70));
for (const m of cohMetrics) {
  const lv = meanStd(lowV.map(g => m.extract(g)).filter(v => isFinite(v)));
  const hv = meanStd(highV.map(g => m.extract(g)).filter(v => isFinite(v)));
  const diff = hv.mu - lv.mu;
  console.log(`  ${m.name.padEnd(18)} ${lv.mu.toFixed(3)}±${lv.std.toFixed(3)}       ${hv.mu.toFixed(3)}±${hv.std.toFixed(3)}       ${(diff>=0?'+':'')+diff.toFixed(3)}`);
}

console.log('\n--- B4: VfResid–a₀ coupling within coherence bins ---');
console.log('    Does coupling strength correlate with RC regularity?\n');

const mrMedian = valid.map(g => g.logMR).sort((a,b) => a-b)[Math.floor(valid.length/2)];
const lowCoh = valid.filter(g => g.logMR < mrMedian);
const highCoh = valid.filter(g => g.logMR >= mrMedian);

const rLowCoh = pearsonR(lowCoh.map(g => g.VfResid), lowCoh.map(g => g.logA0));
const rHighCoh = pearsonR(highCoh.map(g => g.VfResid), highCoh.map(g => g.logA0));

console.log(`  Low coherence (logMR < ${mrMedian.toFixed(2)}, N=${lowCoh.length}):  r = ${rLowCoh.toFixed(3)}`);
console.log(`  High coherence (logMR ≥ ${mrMedian.toFixed(2)}, N=${highCoh.length}): r = ${rHighCoh.toFixed(3)}`);
console.log(`  Delta: ${(rHighCoh-rLowCoh>=0?'+':'')+(rHighCoh-rLowCoh).toFixed(3)}`);

const sig0Median = valid.map(g => g.logSig0).sort((a,b) => a-b)[Math.floor(valid.length/2)];
const lowScatter = valid.filter(g => g.logSig0 < sig0Median);
const highScatter = valid.filter(g => g.logSig0 >= sig0Median);

const rLowScat = pearsonR(lowScatter.map(g => g.VfResid), lowScatter.map(g => g.logA0));
const rHighScat = pearsonR(highScatter.map(g => g.VfResid), highScatter.map(g => g.logA0));

console.log(`  Low scatter (logΣ₀ < ${sig0Median.toFixed(2)}, N=${lowScatter.length}):   r = ${rLowScat.toFixed(3)}`);
console.log(`  High scatter (logΣ₀ ≥ ${sig0Median.toFixed(2)}, N=${highScatter.length}):  r = ${rHighScat.toFixed(3)}`);
console.log(`  Delta: ${(rHighScat-rLowScat>=0?'+':'')+(rHighScat-rLowScat).toFixed(3)}`);

console.log('\n--- B5: Triple partial — coherence + Vflat controlled ---');
const pr_VfMR = partialR2(
  valid.map(g => g.VfResid), valid.map(g => g.logA0),
  [valid.map(g => g.logVflat), valid.map(g => g.logMR), valid.map(g => g.logSig0)]
);
console.log(`  partial r(VfResid, a₀ | Vflat, logMR, logΣ₀) = ${pr_VfMR.toFixed(3)}`);
console.log(`  Compare: partial r | Vflat alone = ${pr_Vflat.toFixed(3)}`);
console.log(`  Adding coherence metrics ${pr_VfMR < pr_Vflat ? 'REDUCES' : 'does NOT reduce'} the signal further.`);


const b404_verdict = {
  partialR_allCoherence: pr_allCoh,
  partialR_VflatMRSig: pr_VfMR,
  rLowCoh, rHighCoh,
  coherenceExplainsActivation: (rHighCoh - rLowCoh) > 0.1,
  coherenceAbsorbs: pr_allCoh < rawR - 0.15,
};

console.log('\n  B-VERDICT: Coherence ' + (b404_verdict.coherenceAbsorbs ? 'DOES' : 'does NOT') + ' absorb the channel.');
console.log('  Coherence ' + (b404_verdict.coherenceExplainsActivation ? 'DOES' : 'does NOT') + ' explain the activation ramp.');


console.log('\n\n' + '▓'.repeat(70));
console.log('404C: NON-CIRCULAR MOTION / ASYMMETRY HYPOTHESIS');
console.log('Do bars, warps, streaming motions explain VfResid?');
console.log('▓'.repeat(70));

const thingsGals = valid.filter(g =>
  g.ncmAmp !== null && g.ncmAmp !== undefined && isFinite(g.ncmAmp) &&
  g.lopsidedness !== null && g.lopsidedness !== undefined && isFinite(g.lopsidedness)
);

console.log(`\n  THINGS subset: N=${thingsGals.length}`);

if (thingsGals.length >= 8) {
  console.log('\n--- C1: Asymmetry metrics vs VfResid and a₀ ---\n');

  const asymMetrics = [
    { name: 'ncmAmp', extract: g => g.ncmAmp, desc: 'Non-circular motion amplitude' },
    { name: 'ncmFrac', extract: g => g.ncmFrac, desc: 'NCM fraction' },
    { name: 'lopsidedness', extract: g => g.lopsidedness, desc: 'Kinematic lopsidedness' },
    { name: 'bisymFlow', extract: g => g.bisymFlow, desc: 'Bisymmetric flow (bar/spiral)' },
    { name: 'c1', extract: g => g.c1, desc: 'Cos(1θ) harmonic' },
    { name: 's1', extract: g => g.s1, desc: 'Sin(1θ) harmonic' },
  ];

  console.log('  Metric            r(.,VfRes)  r(.,a₀)   r(.,Vflat)  N');
  console.log('  ' + '─'.repeat(65));

  const asymResults = [];
  for (const m of asymMetrics) {
    const vf = thingsGals.filter(g => {
      const v = m.extract(g);
      return v !== null && v !== undefined && isFinite(v);
    });
    if (vf.length < 5) continue;
    const x = vf.map(g => m.extract(g));
    const rVfR = pearsonR(x, vf.map(g => g.VfResid));
    const rA0 = pearsonR(x, vf.map(g => g.logA0));
    const rVfl = pearsonR(x, vf.map(g => g.logVflat));
    console.log(`  ${m.name.padEnd(18)} ${(rVfR>=0?'+':'')+rVfR.toFixed(3)}     ${(rA0>=0?'+':'')+rA0.toFixed(3)}     ${(rVfl>=0?'+':'')+rVfl.toFixed(3)}     ${vf.length}`);
    asymResults.push({ name: m.name, rVfResid: rVfR, rA0, rVflat: rVfl, N: vf.length });
  }

  console.log('\n--- C2: Does asymmetry absorb the VfResid–a₀ coupling? ---\n');

  const thVfR = thingsGals.map(g => g.VfResid);
  const thA0 = thingsGals.map(g => g.logA0);
  const rawTh = pearsonR(thVfR, thA0);

  console.log(`  Raw r(VfResid, a₀) in THINGS subset: ${rawTh.toFixed(3)} (N=${thingsGals.length})`);

  const pr_ncm = partialR(thVfR, thA0, thingsGals.map(g => g.ncmAmp));
  const pr_lop = partialR(thVfR, thA0, thingsGals.map(g => g.lopsidedness));
  const pr_bsf = partialR(thVfR, thA0, thingsGals.map(g => g.bisymFlow));

  console.log(`  partial r | ncmAmp:       ${pr_ncm.toFixed(3)} (change ${(pr_ncm-rawTh>=0?'+':'')+(pr_ncm-rawTh).toFixed(3)})`);
  console.log(`  partial r | lopsidedness: ${pr_lop.toFixed(3)} (change ${(pr_lop-rawTh>=0?'+':'')+(pr_lop-rawTh).toFixed(3)})`);
  console.log(`  partial r | bisymFlow:    ${pr_bsf.toFixed(3)} (change ${(pr_bsf-rawTh>=0?'+':'')+(pr_bsf-rawTh).toFixed(3)})`);

  const pr_allAsym = partialR2(thVfR, thA0, [
    thingsGals.map(g => g.ncmAmp),
    thingsGals.map(g => g.lopsidedness),
    thingsGals.map(g => g.bisymFlow),
  ]);
  console.log(`  partial r | ALL asym:     ${pr_allAsym.toFixed(3)} (change ${(pr_allAsym-rawTh>=0?'+':'')+(pr_allAsym-rawTh).toFixed(3)})`);

  console.log('\n--- C3: Asymmetry at low-V vs high-V ---');
  console.log('    If NCM explains suppression at low-V, asymmetry should be HIGHER there\n');

  const thLow = thingsGals.filter(g => g.Vflat < 120);
  const thHigh = thingsGals.filter(g => g.Vflat >= 120);
  console.log(`  low-V: N=${thLow.length}, high-V: N=${thHigh.length}`);

  for (const m of asymMetrics.slice(0, 4)) {
    const lv = thLow.map(g => m.extract(g)).filter(v => v !== null && isFinite(v));
    const hv = thHigh.map(g => m.extract(g)).filter(v => v !== null && isFinite(v));
    if (lv.length < 2 || hv.length < 2) continue;
    const ls = meanStd(lv);
    const hs = meanStd(hv);
    console.log(`  ${m.name.padEnd(16)} low=${ls.mu.toFixed(3)}±${ls.std.toFixed(3)}   high=${hs.mu.toFixed(3)}±${hs.std.toFixed(3)}`);
  }

  const c404_verdict = {
    partialR_allAsym: pr_allAsym,
    rawTHINGS: rawTh,
    asymAbsorbs: pr_allAsym < rawTh - 0.15,
  };
  console.log('\n  C-VERDICT: Non-circular motions ' + (c404_verdict.asymAbsorbs ? 'DO' : 'do NOT') + ' absorb the channel.');

} else {
  console.log('  Insufficient THINGS data for asymmetry tests.');
}


console.log('\n\n' + '▓'.repeat(70));
console.log('404A: EMERGENT ACCELERATION SCALE HYPOTHESIS');
console.log('Does the coupling reflect an acceleration-dependent phenomenon?');
console.log('▓'.repeat(70));

console.log('\n--- A1: Acceleration-scale proxies vs VfResid ---\n');

const accMetrics = [
  { name: 'logAccScale', extract: g => g.logAccScale, desc: 'V²/R proxy (2logV−logR)' },
  { name: 'logGbarEff', extract: g => g.logGbarEff, desc: 'Effective baryonic acceleration' },
  { name: 'logSurfDens', extract: g => g.logSurfDens, desc: 'Surface mass density (M/R²)' },
  { name: 'logSBdisk', extract: g => g.logSBdisk, desc: 'Disk surface brightness' },
  { name: 'logSBeff', extract: g => Math.log10(Math.max(g.SBeff || 0.01, 0.01)), desc: 'Effective surface brightness' },
];

console.log('  Metric            r(.,VfRes)  r(.,a₀)   partial r(VfR,a₀|metric)  desc');
console.log('  ' + '─'.repeat(80));

const accResults = [];
for (const m of accMetrics) {
  const vf = valid.filter(g => isFinite(m.extract(g)));
  if (vf.length < 10) continue;
  const x = vf.map(g => m.extract(g));
  const rVfR = pearsonR(x, vf.map(g => g.VfResid));
  const rA0 = pearsonR(x, vf.map(g => g.logA0));
  const pr = partialR(vf.map(g => g.VfResid), vf.map(g => g.logA0), x);
  console.log(`  ${m.name.padEnd(18)} ${(rVfR>=0?'+':'')+rVfR.toFixed(3)}     ${(rA0>=0?'+':'')+rA0.toFixed(3)}     ${(pr>=0?'+':'')+pr.toFixed(3)}                     ${m.desc}`);
  accResults.push({ name: m.name, rVfResid: rVfR, rA0, partialR: pr });
}

console.log('\n--- A2: Is VfResid an acceleration residual? ---');
console.log('    Test: after removing acceleration-scale dependence, does VfResid–a₀ remain?\n');

const pr_acc = partialR(valid.map(g => g.VfResid), valid.map(g => g.logA0), valid.map(g => g.logAccScale));
const pr_gbar = partialR(valid.map(g => g.VfResid), valid.map(g => g.logA0), valid.map(g => g.logGbarEff));
const pr_surfD = partialR(valid.map(g => g.VfResid), valid.map(g => g.logA0), valid.map(g => g.logSurfDens));

const pr_allAcc = partialR2(
  valid.map(g => g.VfResid), valid.map(g => g.logA0),
  [valid.map(g => g.logAccScale), valid.map(g => g.logSurfDens), valid.map(g => g.logSBdisk)]
);

console.log(`  Raw:                                 r = ${rawR.toFixed(3)}`);
console.log(`  partial r | logAccScale:             r = ${pr_acc.toFixed(3)}`);
console.log(`  partial r | logGbarEff:              r = ${pr_gbar.toFixed(3)}`);
console.log(`  partial r | logSurfDens:             r = ${pr_surfD.toFixed(3)}`);
console.log(`  partial r | ALL 3 acc proxies:       r = ${pr_allAcc.toFixed(3)}`);

console.log('\n--- A3: Acceleration regime bins ---');
console.log('    Does the coupling correlate with WHERE on the acc scale the galaxy sits?\n');

const accSorted = [...valid].sort((a, b) => a.logAccScale - b.logAccScale);
const nTertile = Math.floor(accSorted.length / 3);
const lowAcc = accSorted.slice(0, nTertile);
const midAcc = accSorted.slice(nTertile, 2 * nTertile);
const highAcc = accSorted.slice(2 * nTertile);

const rLowAcc = pearsonR(lowAcc.map(g => g.VfResid), lowAcc.map(g => g.logA0));
const rMidAcc = pearsonR(midAcc.map(g => g.VfResid), midAcc.map(g => g.logA0));
const rHighAcc = pearsonR(highAcc.map(g => g.VfResid), highAcc.map(g => g.logA0));

console.log(`  Low acc tertile  (N=${lowAcc.length}):  r = ${rLowAcc.toFixed(3)}  mean logAccScale = ${meanStd(lowAcc.map(g=>g.logAccScale)).mu.toFixed(2)}`);
console.log(`  Mid acc tertile  (N=${midAcc.length}):  r = ${rMidAcc.toFixed(3)}  mean logAccScale = ${meanStd(midAcc.map(g=>g.logAccScale)).mu.toFixed(2)}`);
console.log(`  High acc tertile (N=${highAcc.length}): r = ${rHighAcc.toFixed(3)}  mean logAccScale = ${meanStd(highAcc.map(g=>g.logAccScale)).mu.toFixed(2)}`);

const a404_verdict = {
  partialR_allAcc: pr_allAcc,
  accExplainsActivation: (rHighAcc - rLowAcc) > 0.1,
  accAbsorbs: pr_allAcc < rawR - 0.15,
};

console.log('\n  A-VERDICT: Acceleration-scale proxies ' + (a404_verdict.accAbsorbs ? 'DO' : 'do NOT') + ' absorb the channel.');
console.log('  Acceleration tertiles ' + (a404_verdict.accExplainsActivation ? 'DO' : 'do NOT') + ' show activation ramp.');


console.log('\n\n' + '▓'.repeat(70));
console.log('404D: STRUCTURAL HALO TRANSITION HYPOTHESIS');
console.log('Is there a qualitative halo profile change with mass?');
console.log('▓'.repeat(70));

const ext = d300.augmentedGalaxies.filter(g => g.Vflat > 0 && g.haloK !== undefined && g.lhOuter !== undefined);

console.log(`\n  External sample with haloK: N=${ext.length}`);

if (ext.length >= 10) {
  console.log('\n--- D1: Halo shape parameters vs mass ---\n');

  const rHK_V = pearsonR(ext.map(g => g.haloK), ext.map(g => g.logVflat));
  const rLH_V = pearsonR(ext.map(g => g.lhOuter), ext.map(g => g.logVflat));
  const rHK_A = pearsonR(ext.map(g => g.haloK), ext.map(g => g.logA0));
  const rLH_A = pearsonR(ext.map(g => g.lhOuter), ext.map(g => g.logA0));
  const rHK_VfR = pearsonR(ext.map(g => g.haloK), ext.map(g => g.VfResid));
  const rLH_VfR = pearsonR(ext.map(g => g.lhOuter), ext.map(g => g.VfResid));

  console.log(`  r(haloK, logVflat) = ${rHK_V.toFixed(3)}`);
  console.log(`  r(lhOuter, logVflat) = ${rLH_V.toFixed(3)}`);
  console.log(`  r(haloK, logA0) = ${rHK_A.toFixed(3)}`);
  console.log(`  r(lhOuter, logA0) = ${rLH_A.toFixed(3)}`);
  console.log(`  r(haloK, VfResid) = ${rHK_VfR.toFixed(3)}`);
  console.log(`  r(lhOuter, VfResid) = ${rLH_VfR.toFixed(3)}`);

  console.log('\n--- D2: Does halo shape absorb VfResid–a₀? ---\n');

  const extVfR = ext.map(g => g.VfResid);
  const extA0 = ext.map(g => g.logA0);
  const rawExt = pearsonR(extVfR, extA0);

  const pr_hk = partialR(extVfR, extA0, ext.map(g => g.haloK));
  const pr_lh = partialR(extVfR, extA0, ext.map(g => g.lhOuter));
  const pr_allHalo = partialR2(extVfR, extA0, [ext.map(g => g.haloK), ext.map(g => g.lhOuter)]);

  console.log(`  Raw r(VfResid, a₀) in ext sample: ${rawExt.toFixed(3)}`);
  console.log(`  partial r | haloK:     ${pr_hk.toFixed(3)}`);
  console.log(`  partial r | lhOuter:   ${pr_lh.toFixed(3)}`);
  console.log(`  partial r | both:      ${pr_allHalo.toFixed(3)}`);

  console.log('\n--- D3: Halo shape at low-V vs high-V ---\n');

  const extLow = ext.filter(g => g.Vflat < 120);
  const extHigh = ext.filter(g => g.Vflat >= 120);

  console.log(`  low-V (N=${extLow.length}): haloK=${meanStd(extLow.map(g=>g.haloK)).mu.toFixed(3)}±${meanStd(extLow.map(g=>g.haloK)).std.toFixed(3)}, lhOuter=${meanStd(extLow.map(g=>g.lhOuter)).mu.toFixed(3)}±${meanStd(extLow.map(g=>g.lhOuter)).std.toFixed(3)}`);
  console.log(`  high-V (N=${extHigh.length}): haloK=${meanStd(extHigh.map(g=>g.haloK)).mu.toFixed(3)}±${meanStd(extHigh.map(g=>g.haloK)).std.toFixed(3)}, lhOuter=${meanStd(extHigh.map(g=>g.lhOuter)).mu.toFixed(3)}±${meanStd(extHigh.map(g=>g.lhOuter)).std.toFixed(3)}`);

  const d404_verdict = {
    partialR_allHalo: pr_allHalo,
    haloShapeChanges: Math.abs(meanStd(extHigh.map(g=>g.haloK)).mu - meanStd(extLow.map(g=>g.haloK)).mu) > 0.1,
    haloAbsorbs: pr_allHalo < rawExt - 0.15,
  };

  console.log('\n  D-VERDICT: Halo structure ' + (d404_verdict.haloAbsorbs ? 'DOES' : 'does NOT') + ' absorb the channel.');
  console.log('  Halo shape ' + (d404_verdict.haloShapeChanges ? 'DOES' : 'does NOT') + ' change qualitatively with mass.');

} else {
  console.log('  Insufficient halo shape data.');
}


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 404: GRAND SCORECARD');
console.log('═'.repeat(70));

console.log('\n  Three empirical signatures each scenario must explain:');
console.log('  S1: Gradual activation ramp (r increases ~0.4→0.84 with Vflat)');
console.log('  S2: Partial r(VfResid,a₀|Vflat) = 0.788 (Vflat-independence)');
console.log('  S3: Irreducible residual r=0.538 after 5-var structural control');
console.log('');

function score(label, partialR_val, absorbsThreshold, activationTest, activationLabel) {
  const absorbs = partialR_val < rawR - 0.15;
  const explains = activationTest;
  let s1 = '?', s2 = '?', s3 = '?';

  if (absorbs) {
    s3 = 'YES — absorbs signal';
    s1 = explains ? 'YES' : 'PARTIAL';
  } else {
    s3 = 'NO — signal persists';
    s1 = explains ? 'PARTIAL' : 'NO';
  }
  s2 = absorbs ? 'YES' : 'NO';

  const total = [s1, s2, s3].filter(s => s.startsWith('YES')).length;
  return { label, s1, s2, s3, total, partialR: partialR_val };
}

const scores = [
  {
    label: '404B: Dynamical Coherence',
    partialR: pr_allCoh,
    absorbs: pr_allCoh < rawR - 0.15,
    activationTest: (rHighCoh - rLowCoh) > 0.1,
  },
  {
    label: '404A: Emergent Acceleration',
    partialR: pr_allAcc,
    absorbs: pr_allAcc < rawR - 0.15,
    activationTest: (rHighAcc - rLowAcc) > 0.1,
  },
];

console.log('  Scenario                      partial r(all)  Absorbs?  Activation?  Score');
console.log('  ' + '─'.repeat(75));

for (const s of scores) {
  const scoreVal = (s.absorbs ? 1 : 0) + (s.activationTest ? 1 : 0);
  console.log(`  ${s.label.padEnd(32)} ${s.partialR.toFixed(3)}          ${s.absorbs ? 'YES' : 'NO '}       ${s.activationTest ? 'YES' : 'NO '}          ${scoreVal}/2`);
}

if (thingsGals.length >= 8) {
  const pr_allAsym2 = partialR2(
    thingsGals.map(g => g.VfResid), thingsGals.map(g => g.logA0),
    [thingsGals.map(g => g.ncmAmp), thingsGals.map(g => g.lopsidedness), thingsGals.map(g => g.bisymFlow)]
  );
  const rawTh2 = pearsonR(thingsGals.map(g => g.VfResid), thingsGals.map(g => g.logA0));
  console.log(`  404C: Non-circular motions    ${pr_allAsym2.toFixed(3)}          ${pr_allAsym2 < rawTh2 - 0.15 ? 'YES' : 'NO '}       ?            ?`);
}

if (ext.length >= 10) {
  const pr_allHalo2 = partialR2(
    ext.map(g => g.VfResid), ext.map(g => g.logA0),
    [ext.map(g => g.haloK), ext.map(g => g.lhOuter)]
  );
  const rawExt2 = pearsonR(ext.map(g => g.VfResid), ext.map(g => g.logA0));
  const extLow2 = ext.filter(g => g.Vflat < 120);
  const extHigh2 = ext.filter(g => g.Vflat >= 120);
  const hkChange = Math.abs(meanStd(extHigh2.map(g=>g.haloK)).mu - meanStd(extLow2.map(g=>g.haloK)).mu) > 0.1;
  console.log(`  404D: Halo transition         ${pr_allHalo2.toFixed(3)}          ${pr_allHalo2 < rawExt2 - 0.15 ? 'YES' : 'NO '}       ${hkChange ? 'YES' : 'NO '}          ?`);
}


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 404 SYNTHESIS');
console.log('═'.repeat(70));

console.log('\nKey findings:');
console.log(`  1. Coherence metrics: partial r after controlling ALL = ${pr_allCoh.toFixed(3)} (raw ${rawR.toFixed(3)})`);
console.log(`  2. Acceleration proxies: partial r after ALL = ${pr_allAcc.toFixed(3)}`);
console.log(`  3. Coherence + Vflat + structure combined: partial r = ${pr_VfMR.toFixed(3)}`);
console.log(`  4. Vflat control alone: partial r = ${pr_Vflat.toFixed(3)}`);

const survivors = [];
if (pr_allCoh >= rawR - 0.15) survivors.push('Coherence does NOT absorb');
if (pr_allAcc >= rawR - 0.15) survivors.push('Acceleration does NOT absorb');
console.log('\n  ' + survivors.join('\n  '));

const outPath = path.join(__dirname, '..', 'public', 'phase404-scenario-discrimination.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '404',
  title: 'Scenario Discrimination',
  timestamp: new Date().toISOString(),
  sample: { N: valid.length, rawR },
  partialR_Vflat: pr_Vflat,
  structResidR: rStructResid,
  B_coherence: {
    metrics: cohResults,
    partialR_allCoherence: pr_allCoh,
    partialR_VflatMRSig: pr_VfMR,
    rLowCoh, rHighCoh,
    verdict: pr_allCoh < rawR - 0.15 ? 'ABSORBS' : 'DOES_NOT_ABSORB'
  },
  A_acceleration: {
    metrics: accResults,
    partialR_allAcc: pr_allAcc,
    rByTertile: { low: rLowAcc, mid: rMidAcc, high: rHighAcc },
    verdict: pr_allAcc < rawR - 0.15 ? 'ABSORBS' : 'DOES_NOT_ABSORB'
  },
}, null, 2));

console.log(`\nSaved: ${outPath}`);
