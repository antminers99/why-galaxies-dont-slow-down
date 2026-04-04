/**
 * Phase 58: Environmental Processing History
 * 
 * Two-pronged approach:
 * 58a: Karachentsev tidal index Θ₁ (N=15, pilot)
 * 58b: Morphology-gas deficit (N=56, full) — conditions on T, not L36
 * 
 * Key question: does a processing proxy explain a₀ variation
 * BEYOND envCode and logMHI?
 */
const fs = require('fs');

const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const p56 = JSON.parse(fs.readFileSync('public/phase56-frozen-baselines.json','utf8'));
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const gals = stageA.galaxies;
const N = gals.length;

function mean(arr) { return arr.reduce((s,v) => s+v, 0) / arr.length; }
function sd(arr) { const m = mean(arr); return Math.sqrt(arr.reduce((s,v) => s + (v-m)**2, 0) / (arr.length-1)); }

function pearson(x, y) {
  const n = x.length;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i]-mx)*(y[i]-my);
    sxx += (x[i]-mx)**2;
    syy += (y[i]-my)**2;
  }
  const r = sxy / Math.sqrt(sxx * syy);
  const t = r * Math.sqrt((n-2) / (1-r*r));
  return { r, t, n };
}

function spearman(x, y) {
  const n = x.length;
  function rank(arr) {
    const sorted = arr.map((v, i) => ({v, i})).sort((a, b) => a.v - b.v);
    const ranks = new Array(n);
    for (let i = 0; i < n; ) {
      let j = i;
      while (j < n && sorted[j].v === sorted[i].v) j++;
      const avgRank = (i + j - 1) / 2 + 1;
      for (let k = i; k < j; k++) ranks[sorted[k].i] = avgRank;
      i = j;
    }
    return ranks;
  }
  const rx = rank(x), ry = rank(y);
  const mx = mean(rx), my = mean(ry);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (rx[i]-mx)*(ry[i]-my);
    sxx += (rx[i]-mx)**2;
    syy += (ry[i]-my)**2;
  }
  return sxy / Math.sqrt(sxx * syy);
}

function ols(Y, X) {
  const n = Y.length, p = X[0].length + 1;
  const Xa = X.map(row => [1, ...row]);
  const XtX = Array.from({length: p}, () => new Array(p).fill(0));
  const XtY = new Array(p).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < p; j++) {
      XtY[j] += Xa[i][j] * Y[i];
      for (let l = 0; l < p; l++) XtX[j][l] += Xa[i][j] * Xa[i][l];
    }
  }
  const beta = solveLinear(XtX, XtY);
  const resid = Y.map((y, i) => y - Xa[i].reduce((s, x, j) => s + x * beta[j], 0));
  const rss = resid.reduce((s, r) => s + r*r, 0);
  const tss = Y.reduce((s, y) => s + (y - mean(Y))**2, 0);
  const r2adj = 1 - (1 - (1-rss/tss)) * (n-1) / (n-p);
  const se = Math.sqrt(rss / (n-p));
  const XtXinv = invertMatrix(XtX);
  const seBeta = beta.map((_, j) => Math.sqrt(Math.abs(XtXinv[j][j]) * rss / (n-p)));
  const tStats = beta.map((b, j) => b / (seBeta[j] || 1e-15));
  return { beta, seBeta, tStats, r2adj, se, resid, rss, tss, n, k: p };
}

function solveLinear(A, b) {
  const n = A.length;
  const M = A.map((row, i) => [...row, b[i]]);
  for (let i = 0; i < n; i++) {
    let maxRow = i;
    for (let j = i+1; j < n; j++) if (Math.abs(M[j][i]) > Math.abs(M[maxRow][i])) maxRow = j;
    [M[i], M[maxRow]] = [M[maxRow], M[i]];
    for (let j = i+1; j < n; j++) {
      const f = M[j][i] / M[i][i];
      for (let k = i; k <= n; k++) M[j][k] -= f * M[i][k];
    }
  }
  const x = new Array(n);
  for (let i = n-1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i+1; j < n; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}

function invertMatrix(A) {
  const n = A.length;
  const M = A.map((row, i) => [...row, ...Array.from({length: n}, (_, j) => i === j ? 1 : 0)]);
  for (let i = 0; i < n; i++) {
    let maxRow = i;
    for (let j = i+1; j < n; j++) if (Math.abs(M[j][i]) > Math.abs(M[maxRow][i])) maxRow = j;
    [M[i], M[maxRow]] = [M[maxRow], M[i]];
    const pivot = M[i][i];
    if (Math.abs(pivot) < 1e-15) { for (let j = 0; j < 2*n; j++) M[i][j] = 0; continue; }
    for (let j = 0; j < 2*n; j++) M[i][j] /= pivot;
    for (let j = 0; j < n; j++) {
      if (j === i) continue;
      const f = M[j][i];
      for (let k = 0; k < 2*n; k++) M[j][k] -= f * M[i][k];
    }
  }
  return M.map(row => row.slice(n));
}

function looCV(Y, X) {
  const n = Y.length;
  let ssPred = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0,i), ...Y.slice(i+1)];
    const Xt = [...X.slice(0,i), ...X.slice(i+1)];
    const fit = ols(Yt, Xt);
    const xi = [1, ...X[i]];
    const pred = xi.reduce((s, x, j) => s + x * fit.beta[j], 0);
    ssPred += (Y[i] - pred) ** 2;
  }
  return Math.sqrt(ssPred / n);
}

function permTest(Y, X, addIdx, nPerm) {
  const Xbase = X.map(row => row.filter((_, j) => j !== addIdx));
  const fitBase = ols(Y, Xbase);
  const fitFull = ols(Y, X);
  const obsImprove = fitBase.rss - fitFull.rss;
  let count = 0;
  for (let p = 0; p < nPerm; p++) {
    const newVar = X.map(row => row[addIdx]);
    for (let i = newVar.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [newVar[i], newVar[j]] = [newVar[j], newVar[i]];
    }
    const Xperm = X.map((row, k) => row.map((v, j) => j === addIdx ? newVar[k] : v));
    const fitPerm = ols(Y, Xperm);
    if (fitBase.rss - fitPerm.rss >= obsImprove) count++;
  }
  return count / nPerm;
}

function bootstrapCI(Y, X, varIdx, nBoot) {
  const coeffs = [];
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({length: Y.length}, () => Math.floor(Math.random() * Y.length));
    try {
      const fit = ols(idx.map(i => Y[i]), idx.map(i => X[i]));
      coeffs.push(fit.beta[varIdx + 1]);
    } catch(e) {}
  }
  coeffs.sort((a,b) => a-b);
  return { lo: coeffs[Math.floor(0.025*coeffs.length)], hi: coeffs[Math.floor(0.975*coeffs.length)], n: coeffs.length };
}

function jackknifeSigns(Y, X, varIdx) {
  const fullSign = Math.sign(ols(Y, X).beta[varIdx + 1]);
  let flips = 0;
  for (let i = 0; i < Y.length; i++) {
    const fit = ols([...Y.slice(0,i),...Y.slice(i+1)], [...X.slice(0,i),...X.slice(i+1)]);
    if (Math.sign(fit.beta[varIdx + 1]) !== fullSign) flips++;
  }
  return flips;
}

const sdY = sd(gals.map(g => g.logA0));
function gapClosedV(looRMS) { return 100 * (1 - looRMS**2 / sdY**2); }

// ======================================================================
// KARACHENTSEV TIDAL INDEX Θ₁ (N=15 within 11 Mpc)
// ======================================================================
const tidalData = {
  'NGC0024':  { theta1: -0.8, md: 'NGC0055' },
  'NGC0891':  { theta1: 2.9, md: 'NGC1023_group' },
  'NGC1705':  { theta1: -1.1, md: 'NGC1744' },
  'NGC2403':  { theta1: 1.5, md: 'M81_group' },
  'NGC2683':  { theta1: 0.0, md: 'NGC2903' },
  'NGC2903':  { theta1: -0.3, md: 'NGC3379_Leo' },
  'NGC2915':  { theta1: -0.7, md: 'NGC2784' },
  'NGC3521':  { theta1: 0.4, md: 'NGC3521_group' },
  'NGC3741':  { theta1: 1.2, md: 'NGC3031_M81' },
  'NGC4559':  { theta1: 1.5, md: 'NGC4565' },
  'NGC5055':  { theta1: 1.8, md: 'NGC5194_M51' },
  'NGC6503':  { theta1: -0.5, md: 'NGC6946' },
  'UGC01281': { theta1: 0.4, md: 'NGC672' },
  'UGC05721': { theta1: -0.9, md: 'NGC3198' },
  'UGC08490': { theta1: 1.7, md: 'NGC5457_M101' }
};

// ======================================================================
// MORPHOLOGY-GAS DEFICIT (N=56)
// Different from HI deficiency: conditions on morphological type T
// ======================================================================
const sparc56 = gals.map(g => {
  const sr = sparc.find(s => s.name === g.name);
  return { ...g, T: sr?.T, L36: sr?.L36, MHI: sr?.MHI, Vflat: sr?.Vflat };
});

// Fit expected logMHI from T (morphological type)
const validT = sparc56.filter(g => g.T !== undefined && g.T !== null);
const Tarr = validT.map(g => g.T);
const logMHI_arr = validT.map(g => g.logMHI);

// Linear fit: logMHI = a + b*T
const fitTMHI = ols(logMHI_arr, Tarr.map(t => [t]));
console.log('Morphology-gas scaling: logMHI = ' + fitTMHI.beta[0].toFixed(3) + ' + ' + fitTMHI.beta[1].toFixed(3) + '*T');
console.log('  r(T, logMHI) = ' + pearson(Tarr, logMHI_arr).r.toFixed(3));

// Compute morphology-gas deficit: DEF_T = expected_logMHI(T) - observed_logMHI
// DEF_T > 0 means less gas than expected for its type → processed/stripped
const morphGasDef = sparc56.map(g => {
  if (g.T === undefined || g.T === null) return null;
  const expected = fitTMHI.beta[0] + fitTMHI.beta[1] * g.T;
  return expected - g.logMHI;
});

// Also compute gas-fraction offset: f_gas relative to type
// f_gas = MHI / (MHI + Mstar) where Mstar ~ 0.5 * L36
const gasFrac = sparc56.map(g => {
  if (!g.L36 || !g.MHI) return null;
  const Mstar = 0.5 * g.L36 * 1e9;
  const Mgas = g.MHI * 1e9;
  return Math.log10(Mgas / (Mgas + Mstar));
});

// Fit gasFrac vs T
const validGF = sparc56.map((g, i) => ({ T: g.T, gf: gasFrac[i] })).filter(g => g.T !== null && g.gf !== null);
const fitGFT = ols(validGF.map(g => g.gf), validGF.map(g => [g.T]));
const gasFracOffset = sparc56.map((g, i) => {
  if (g.T === null || gasFrac[i] === null) return null;
  return gasFrac[i] - (fitGFT.beta[0] + fitGFT.beta[1] * g.T);
});
// gasFracOffset < 0 means LESS gas fraction than expected for type → processed

console.log('Gas fraction scaling: logFgas = ' + fitGFT.beta[0].toFixed(3) + ' + ' + fitGFT.beta[1].toFixed(3) + '*T');

// ======================================================================
// PRINT RESULTS
// ======================================================================
console.log('\n╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  PHASE 58: ENVIRONMENTAL PROCESSING HISTORY                                   ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

// ======================================================================
// 58a: TIDAL INDEX (N=15, pilot)
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  58a: KARACHENTSEV TIDAL INDEX Θ₁ (UNGC 2013)');
console.log('  N = 15 galaxies within 11 Mpc — PILOT ONLY');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const tidal15 = sparc56.filter(g => tidalData[g.name]).map(g => ({
  ...g, theta1: tidalData[g.name].theta1, md: tidalData[g.name].md
}));

console.log('  ┌──────────────┬────────┬────────┬─────────┬─────────┬─────────┬──────────────────┐');
console.log('  │ Galaxy       │ Θ₁     │envCode │ delta_a0│ logMHI  │ rcWig   │ Main Disturber   │');
console.log('  ├──────────────┼────────┼────────┼─────────┼─────────┼─────────┼──────────────────┤');
tidal15.sort((a,b) => b.theta1 - a.theta1).forEach(g => {
  console.log('  │ ' + g.name.padEnd(12) + ' │' +
    g.theta1.toFixed(1).padStart(6) + '  │   ' + g.envCode + '    │' +
    g.delta_a0.toFixed(3).padStart(7) + '  │' +
    g.logMHI.toFixed(2).padStart(7) + '  │' +
    g.rcWiggliness.toFixed(3).padStart(7) + '  │ ' + g.md.padEnd(16) + ' │');
});
console.log('  └──────────────┴────────┴────────┴─────────┴─────────┴─────────┴──────────────────┘\n');

// Correlations
const t15_theta = tidal15.map(g => g.theta1);
const t15_da0 = tidal15.map(g => g.delta_a0);
const t15_logA0 = tidal15.map(g => g.logA0);
const t15_env = tidal15.map(g => g.envCode);

const corrTheta_da0 = pearson(t15_theta, t15_da0);
const corrTheta_env = pearson(t15_theta, t15_env);
const spTheta_da0 = spearman(t15_theta, t15_da0);

console.log('  Pearson r(Θ₁, delta_a0) = ' + corrTheta_da0.r.toFixed(4) + ' (t=' + corrTheta_da0.t.toFixed(2) + ')');
console.log('  Spearman ρ(Θ₁, delta_a0) = ' + spTheta_da0.toFixed(4));
console.log('  Pearson r(Θ₁, envCode) = ' + corrTheta_env.r.toFixed(4));
console.log('  Direction: ' + (corrTheta_da0.r < 0 ? 'NEGATIVE (stronger tidal → lower a0)' : 'POSITIVE (stronger tidal → higher a0)'));
console.log();

// After envCode
const fitEnv15 = ols(t15_logA0, tidal15.map(g => [g.envCode]));
const corrTheta_residEnv = pearson(t15_theta, fitEnv15.resid);
console.log('  After envCode: r(Θ₁, resid) = ' + corrTheta_residEnv.r.toFixed(4) + ' (t=' + corrTheta_residEnv.t.toFixed(2) + ')');
console.log('  → Does Θ₁ add BEYOND envCode? ' + (Math.abs(corrTheta_residEnv.r) > 0.3 ? 'POSSIBLY' : 'NO'));
console.log();

// ======================================================================
// 58b: MORPHOLOGY-GAS DEFICIT (N=56, full)
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  58b: MORPHOLOGY-GAS DEFICIT (full, N=56)');
console.log('  DEF_T = expected_logMHI(T) - observed_logMHI');
console.log('  Conditions on morphological type T, NOT luminosity L3.6');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const Y = gals.map(g => g.logA0);
const mgd = morphGasDef;
const gfo = gasFracOffset;

// Check how DEF_T differs from HI deficiency
const hiDef = gals.map(g => g.hi_deficiency);
const corrDEFT_hiDef = pearson(mgd, hiDef);
const corrDEFT_logMHI = pearson(mgd, gals.map(g => g.logMHI));
console.log('  Cross-correlations:');
console.log('    r(DEF_T, HI_deficiency) = ' + corrDEFT_hiDef.r.toFixed(4));
console.log('    r(DEF_T, logMHI)        = ' + corrDEFT_logMHI.r.toFixed(4));
console.log('    r(HI_def, logMHI)       = ' + pearson(hiDef, gals.map(g => g.logMHI)).r.toFixed(4) + '  (for comparison)');
console.log('  → DEF_T is ' + (Math.abs(corrDEFT_logMHI.r) < 0.3 ? 'WEAKLY' : 'MODERATELY') + ' correlated with logMHI');
console.log('  → This is ' + (Math.abs(corrDEFT_logMHI.r) < Math.abs(pearson(hiDef, gals.map(g => g.logMHI)).r) ? 'LESS' : 'MORE') + ' correlated than HI_def with logMHI');
console.log();

// Raw correlation
const rawDEFT = pearson(mgd, Y);
console.log('  Raw r(DEF_T, logA0) = ' + rawDEFT.r.toFixed(4) + ' (t=' + rawDEFT.t.toFixed(3) + ')');

// After logMHI alone — THE KEY TEST
const fitMHI = ols(Y, gals.map(g => [g.logMHI]));
const corrDEFT_afterMHI = pearson(mgd, fitMHI.resid);
console.log('  |logMHI: r(DEF_T, resid) = ' + corrDEFT_afterMHI.r.toFixed(4) + ' (t=' + corrDEFT_afterMHI.t.toFixed(3) + ')');

// After envCode alone
const fitEnv = ols(Y, gals.map(g => [g.envCode]));
const corrDEFT_afterEnv = pearson(mgd, fitEnv.resid);
console.log('  |envCode: r(DEF_T, resid) = ' + corrDEFT_afterEnv.r.toFixed(4) + ' (t=' + corrDEFT_afterEnv.t.toFixed(3) + ')');

// After logMHI + envCode
const fitME = ols(Y, gals.map(g => [g.logMHI, g.envCode]));
const corrDEFT_afterME = pearson(mgd, fitME.resid);
console.log('  |logMHI+envCode: r(DEF_T, resid) = ' + corrDEFT_afterME.r.toFixed(4) + ' (t=' + corrDEFT_afterME.t.toFixed(3) + ')');
console.log();

// After Baseline A
const X_A = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0]);
const fitA = ols(Y, X_A);
const corrDEFT_afterA = pearson(mgd, fitA.resid);
console.log('  After Baseline A:');
console.log('    r(DEF_T, residA) = ' + corrDEFT_afterA.r.toFixed(4) + ' (t=' + corrDEFT_afterA.t.toFixed(3) + ')');

// Full model A + DEF_T
const X_A_deft = gals.map((g, i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, mgd[i]]);
const fitA_deft = ols(Y, X_A_deft);
console.log('    A + DEF_T: R2adj=' + fitA_deft.r2adj.toFixed(4) + ' (vs A: ' + fitA.r2adj.toFixed(4) + ')');
console.log('    DEF_T coeff=' + fitA_deft.beta[5].toFixed(4) + ' (t=' + fitA_deft.tStats[5].toFixed(3) + ')');

// After Baseline B
const X_B = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun]);
const fitB = ols(Y, X_B);
const corrDEFT_afterB = pearson(mgd, fitB.resid);
console.log('\n  After Baseline B:');
console.log('    r(DEF_T, residB) = ' + corrDEFT_afterB.r.toFixed(4) + ' (t=' + corrDEFT_afterB.t.toFixed(3) + ')');

const X_B_deft = gals.map((g, i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun, mgd[i]]);
const fitB_deft = ols(Y, X_B_deft);
console.log('    B + DEF_T: R2adj=' + fitB_deft.r2adj.toFixed(4) + ' (vs B: ' + fitB.r2adj.toFixed(4) + ')');
console.log('    DEF_T coeff=' + fitB_deft.beta[6].toFixed(4) + ' (t=' + fitB_deft.tStats[6].toFixed(3) + ')');

// envCode after adding DEF_T — does DEF_T absorb envCode?
console.log('\n  envCode after adding DEF_T:');
const X_noEnv_deft = gals.map((g, i) => [g.logMHI, g.rcWiggliness, g.logSigma0, mgd[i]]);
const fitNoEnv_deft = ols(Y, X_noEnv_deft);
const X_withEnv_deft = gals.map((g, i) => [g.logMHI, g.rcWiggliness, g.logSigma0, mgd[i], g.envCode]);
const fitWithEnv_deft = ols(Y, X_withEnv_deft);
console.log('    Without envCode (+ DEF_T): R2adj=' + fitNoEnv_deft.r2adj.toFixed(4));
console.log('    With envCode (+ DEF_T):    R2adj=' + fitWithEnv_deft.r2adj.toFixed(4));
console.log('    envCode t-stat in joint model: ' + fitWithEnv_deft.tStats[5].toFixed(3));
console.log('    → DEF_T absorbs envCode? ' + (Math.abs(fitWithEnv_deft.tStats[5]) < 2 ? 'PARTIALLY/YES' : 'NO'));

// LOO
console.log('\n  LOO Cross-Validation (variance-based gap%):');
const looA = looCV(Y, X_A);
const looA_deft = looCV(Y, X_A_deft);
const looB = looCV(Y, X_B);
const looB_deft = looCV(Y, X_B_deft);
console.log('    A:        ' + gapClosedV(looA).toFixed(1) + '%');
console.log('    A + DEF_T: ' + gapClosedV(looA_deft).toFixed(1) + '% (delta=' + (gapClosedV(looA_deft)-gapClosedV(looA)).toFixed(1) + 'pp)');
console.log('    B:        ' + gapClosedV(looB).toFixed(1) + '%');
console.log('    B + DEF_T: ' + gapClosedV(looB_deft).toFixed(1) + '% (delta=' + (gapClosedV(looB_deft)-gapClosedV(looB)).toFixed(1) + 'pp)');

// Permutation test
console.log('\n  Permutation (5000):');
const permA = permTest(Y, X_A_deft, 4, 5000);
const permB = permTest(Y, X_B_deft, 5, 5000);
console.log('    p(A + DEF_T) = ' + permA.toFixed(4));
console.log('    p(B + DEF_T) = ' + permB.toFixed(4));

// Bootstrap
console.log('\n  Bootstrap 95% CI (2000):');
const bootA = bootstrapCI(Y, X_A_deft, 4, 2000);
const bootB = bootstrapCI(Y, X_B_deft, 5, 2000);
console.log('    A + DEF_T: [' + bootA.lo.toFixed(4) + ', ' + bootA.hi.toFixed(4) + '] ' + (bootA.lo*bootA.hi > 0 ? 'excludes zero' : 'includes zero'));
console.log('    B + DEF_T: [' + bootB.lo.toFixed(4) + ', ' + bootB.hi.toFixed(4) + '] ' + (bootB.lo*bootB.hi > 0 ? 'excludes zero' : 'includes zero'));

// Jackknife
const jkA = jackknifeSigns(Y, X_A_deft, 4);
const jkB = jackknifeSigns(Y, X_B_deft, 5);
console.log('\n  Jackknife sign flips:');
console.log('    A + DEF_T: ' + jkA + '/' + N);
console.log('    B + DEF_T: ' + jkB + '/' + N);

// ======================================================================
// 58c: GAS-FRACTION OFFSET (alternative proxy)
// ======================================================================
console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  58c: GAS-FRACTION OFFSET (alternative processing proxy)');
console.log('  Residual of log(f_gas) vs T — negative = processed');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const corrGFO_logMHI = pearson(gfo, gals.map(g => g.logMHI));
const corrGFO_hiDef = pearson(gfo, hiDef);
const rawGFO = pearson(gfo, Y);
const corrGFO_afterMHI = pearson(gfo, fitMHI.resid);
const corrGFO_afterA = pearson(gfo, fitA.resid);
const corrGFO_afterB = pearson(gfo, fitB.resid);

console.log('  r(GFO, logMHI) = ' + corrGFO_logMHI.r.toFixed(4));
console.log('  r(GFO, HI_def) = ' + corrGFO_hiDef.r.toFixed(4));
console.log('  Raw r(GFO, logA0) = ' + rawGFO.r.toFixed(4) + ' (t=' + rawGFO.t.toFixed(3) + ')');
console.log('  |logMHI: r = ' + corrGFO_afterMHI.r.toFixed(4) + ' (t=' + corrGFO_afterMHI.t.toFixed(3) + ')');
console.log('  |Baseline A: r = ' + corrGFO_afterA.r.toFixed(4) + ' (t=' + corrGFO_afterA.t.toFixed(3) + ')');
console.log('  |Baseline B: r = ' + corrGFO_afterB.r.toFixed(4) + ' (t=' + corrGFO_afterB.t.toFixed(3) + ')');

// ======================================================================
// VERDICT
// ======================================================================
console.log('\n╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  FINAL VERDICT                                                                ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

// Check what survived
const survivesMHI_deft = Math.abs(corrDEFT_afterMHI.t) > 1.65;
const survivesA_deft = Math.abs(fitA_deft.tStats[5]) > 1.65;
const survivesB_deft = Math.abs(fitB_deft.tStats[6]) > 1.65;
const looImprovesA = gapClosedV(looA_deft) > gapClosedV(looA);
const looImprovesB = gapClosedV(looB_deft) > gapClosedV(looB);

let verdict;
if (survivesB_deft && looImprovesB && permB < 0.05) {
  verdict = 'CONFIRMED';
} else if (survivesA_deft && looImprovesA) {
  verdict = 'PARTIAL';
} else if (survivesMHI_deft) {
  verdict = 'PARTIAL (survives logMHI but not baselines)';
} else {
  verdict = 'FAIL';
}

console.log('  58a (Θ₁ pilot, N=15): r=' + corrTheta_da0.r.toFixed(3) + ', after envCode: r=' + corrTheta_residEnv.r.toFixed(3));
console.log('  58b (DEF_T, N=56):');
console.log('    Raw: r=' + rawDEFT.r.toFixed(4) + ', t=' + rawDEFT.t.toFixed(3));
console.log('    |logMHI: r=' + corrDEFT_afterMHI.r.toFixed(4) + ', t=' + corrDEFT_afterMHI.t.toFixed(3));
console.log('    |envCode: r=' + corrDEFT_afterEnv.r.toFixed(4) + ', t=' + corrDEFT_afterEnv.t.toFixed(3));
console.log('    |logMHI+envCode: r=' + corrDEFT_afterME.r.toFixed(4) + ', t=' + corrDEFT_afterME.t.toFixed(3));
console.log('    |Baseline A: t=' + fitA_deft.tStats[5].toFixed(3));
console.log('    |Baseline B: t=' + fitB_deft.tStats[6].toFixed(3));
console.log('    LOO: A delta=' + (gapClosedV(looA_deft)-gapClosedV(looA)).toFixed(1) + 'pp, B delta=' + (gapClosedV(looB_deft)-gapClosedV(looB)).toFixed(1) + 'pp');
console.log('    Perm: A=' + permA.toFixed(4) + ', B=' + permB.toFixed(4));
console.log('    Boot CI A: [' + bootA.lo.toFixed(4) + ', ' + bootA.hi.toFixed(4) + ']');
console.log('    Boot CI B: [' + bootB.lo.toFixed(4) + ', ' + bootB.hi.toFixed(4) + ']');
console.log('    JK flips: A=' + jkA + '/' + N + ', B=' + jkB + '/' + N);
console.log('  58c (GFO, N=56):');
console.log('    |logMHI: r=' + corrGFO_afterMHI.r.toFixed(4));
console.log('    |Baseline B: r=' + corrGFO_afterB.r.toFixed(4));
console.log();
console.log('  ═══════════════════════════════════════════════');
console.log('  VERDICT: ' + verdict);
console.log('  ═══════════════════════════════════════════════');
console.log();

if (verdict.startsWith('FAIL')) {
  console.log('  INTERPRETATION:');
  console.log('  Morphology-gas deficit (DEF_T) does NOT add information beyond');
  console.log('  the existing baselines. Despite conditioning on morphological type');
  console.log('  rather than luminosity, the processing signal is still absorbed');
  console.log('  by logMHI + envCode. The environmental effect captured by envCode');
  console.log('  remains a binary marker that cannot be refined with available data.');
}

// Save
const output = {
  phase: 58,
  program: 'New Program',
  subphases: {
    '58a_tidal': {
      type: 'pilot', n: 15,
      source: 'Karachentsev+2013 UNGC',
      variable: 'Theta_1',
      corrWithDeltaA0: corrTheta_da0,
      afterEnvCode: corrTheta_residEnv,
      tidalData: tidal15.map(g => ({ name: g.name, theta1: g.theta1, envCode: g.envCode, delta_a0: g.delta_a0, md: g.md }))
    },
    '58b_morphGasDef': {
      type: 'full', n: 56,
      definition: 'DEF_T = expected_logMHI(T) - observed_logMHI',
      scaling: { intercept: fitTMHI.beta[0], slope: fitTMHI.beta[1] },
      raw: rawDEFT,
      afterLogMHI: corrDEFT_afterMHI,
      afterEnvCode: corrDEFT_afterEnv,
      afterMHI_env: corrDEFT_afterME,
      afterBaselineA: { r: corrDEFT_afterA.r, coeff: fitA_deft.beta[5], t: fitA_deft.tStats[5] },
      afterBaselineB: { r: corrDEFT_afterB.r, coeff: fitB_deft.beta[6], t: fitB_deft.tStats[6] },
      loo: {
        A: gapClosedV(looA), A_deft: gapClosedV(looA_deft),
        B: gapClosedV(looB), B_deft: gapClosedV(looB_deft)
      },
      permutation: { A: permA, B: permB },
      bootstrap: { A: bootA, B: bootB },
      jackknife: { A: jkA, B: jkB }
    },
    '58c_gasFracOffset': {
      raw: rawGFO,
      afterLogMHI: corrGFO_afterMHI,
      afterBaselineB: corrGFO_afterB
    }
  },
  verdict
};

fs.writeFileSync('public/phase58-processing-history.json', JSON.stringify(output, null, 2));
console.log('  Results saved to public/phase58-processing-history.json');
