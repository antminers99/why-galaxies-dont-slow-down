/**
 * Phase 57 (New Program): HI Deficiency
 * Full formal test against frozen Baselines A and B
 * 
 * HI deficiency = log(MHI_expected) - log(MHI_observed)
 * where MHI_expected comes from L3.6-MHI scaling relation
 * 
 * Key question: Does gas DEFICIENCY (relative to expectations)
 * explain a0 variation BEYOND raw gas mass (logMHI)?
 */
const fs = require('fs');

const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const p56 = JSON.parse(fs.readFileSync('public/phase56-frozen-baselines.json','utf8'));

const gals = stageA.galaxies;
const N = gals.length;

// ======================================================================
// UTILITY FUNCTIONS
// ======================================================================
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
  const df = n - 2;
  const p = 2 * tDistCDF(-Math.abs(t), df);
  return { r, t, p, n };
}

function tDistCDF(t, df) {
  const x = df / (df + t*t);
  return 0.5 * betaInc(df/2, 0.5, x);
}

function betaInc(a, b, x) {
  if (x === 0 || x === 1) return x === 0 ? 1 : 0;
  const lbeta = lgamma(a) + lgamma(b) - lgamma(a+b);
  const front = Math.exp(Math.log(x)*a + Math.log(1-x)*b - lbeta) / a;
  let sum = 0, term = 1;
  for (let n = 0; n < 200; n++) {
    sum += term;
    term *= (a+n)*(a+b+n)/((a+1+n)*(n+1)) * x;
    if (Math.abs(term) < 1e-12) break;
  }
  return 1 - front * sum;
}

function lgamma(x) {
  const c = [76.18009172947146,-86.50532032941677,24.01409824083091,
    -1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5];
  let y = x, tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) ser += c[j] / ++y;
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

function ols(Y, X) {
  const n = Y.length, k = X[0].length;
  // Add intercept
  const Xa = X.map(row => [1, ...row]);
  const p = Xa[0].length;
  // X'X
  const XtX = Array.from({length: p}, () => new Array(p).fill(0));
  const XtY = new Array(p).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < p; j++) {
      XtY[j] += Xa[i][j] * Y[i];
      for (let l = 0; l < p; l++) {
        XtX[j][l] += Xa[i][j] * Xa[i][l];
      }
    }
  }
  const beta = solveLinear(XtX, XtY);
  const resid = Y.map((y, i) => y - Xa[i].reduce((s, x, j) => s + x * beta[j], 0));
  const rss = resid.reduce((s, r) => s + r*r, 0);
  const tss = Y.reduce((s, y) => s + (y - mean(Y))**2, 0);
  const r2 = 1 - rss / tss;
  const r2adj = 1 - (1-r2) * (n-1) / (n-p);
  const se = Math.sqrt(rss / (n-p));
  
  // Standard errors of coefficients
  const XtXinv = invertMatrix(XtX);
  const seBeta = beta.map((_, j) => Math.sqrt(XtXinv[j][j] * rss / (n-p)));
  const tStats = beta.map((b, j) => b / seBeta[j]);
  const pVals = tStats.map(t => 2 * tDistCDF(-Math.abs(t), n-p));
  
  return { beta, seBeta, tStats, pVals, r2, r2adj, se, resid, rss, tss, n, k: p };
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
    const Ytrain = [...Y.slice(0,i), ...Y.slice(i+1)];
    const Xtrain = [...X.slice(0,i), ...X.slice(i+1)];
    const fit = ols(Ytrain, Xtrain);
    const xi = [1, ...X[i]];
    const pred = xi.reduce((s, x, j) => s + x * fit.beta[j], 0);
    ssPred += (Y[i] - pred) ** 2;
  }
  const looRMS = Math.sqrt(ssPred / n);
  return looRMS;
}

function permTest(Y, X, addIdx, nPerm) {
  // Test significance of adding variable at addIdx
  const Xbase = X.map(row => row.filter((_, j) => j !== addIdx));
  const Xfull = X;
  
  const fitBase = ols(Y, Xbase);
  const fitFull = ols(Y, Xfull);
  const obsImprove = fitBase.rss - fitFull.rss;
  
  let count = 0;
  for (let p = 0; p < nPerm; p++) {
    // Permute the test variable
    const newVar = X.map(row => row[addIdx]);
    // Fisher-Yates shuffle
    for (let i = newVar.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [newVar[i], newVar[j]] = [newVar[j], newVar[i]];
    }
    const Xperm = X.map((row, k) => row.map((v, j) => j === addIdx ? newVar[k] : v));
    const fitPerm = ols(Y, Xperm);
    const permImprove = fitBase.rss - fitPerm.rss;
    if (permImprove >= obsImprove) count++;
  }
  return count / nPerm;
}

function bootstrapCI(Y, X, varIdx, nBoot) {
  const coeffs = [];
  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({length: Y.length}, () => Math.floor(Math.random() * Y.length));
    const Yb = idx.map(i => Y[i]);
    const Xb = idx.map(i => X[i]);
    try {
      const fit = ols(Yb, Xb);
      coeffs.push(fit.beta[varIdx + 1]); // +1 for intercept
    } catch(e) {}
  }
  coeffs.sort((a,b) => a-b);
  const lo = coeffs[Math.floor(0.025 * coeffs.length)];
  const hi = coeffs[Math.floor(0.975 * coeffs.length)];
  return { lo, hi, median: coeffs[Math.floor(0.5 * coeffs.length)], n: coeffs.length };
}

function jackknifeSigns(Y, X, varIdx) {
  let flips = 0;
  const fullFit = ols(Y, X);
  const fullSign = Math.sign(fullFit.beta[varIdx + 1]);
  for (let i = 0; i < Y.length; i++) {
    const Yj = [...Y.slice(0,i), ...Y.slice(i+1)];
    const Xj = [...X.slice(0,i), ...X.slice(i+1)];
    const fit = ols(Yj, Xj);
    if (Math.sign(fit.beta[varIdx + 1]) !== fullSign) flips++;
  }
  return flips;
}

function kfoldCV(Y, X, k) {
  const n = Y.length;
  const indices = Array.from({length: n}, (_, i) => i);
  // Shuffle
  for (let i = n-1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i+1));
    [indices[i], indices[j]] = [indices[j], indices[i]];
  }
  const foldSize = Math.floor(n / k);
  let ssPred = 0, count = 0;
  for (let f = 0; f < k; f++) {
    const testIdx = indices.slice(f * foldSize, (f+1) * foldSize);
    const trainIdx = indices.filter(i => !testIdx.includes(i));
    const Ytrain = trainIdx.map(i => Y[i]);
    const Xtrain = trainIdx.map(i => X[i]);
    const fit = ols(Ytrain, Xtrain);
    for (const i of testIdx) {
      const xi = [1, ...X[i]];
      const pred = xi.reduce((s, x, j) => s + x * fit.beta[j], 0);
      ssPred += (Y[i] - pred) ** 2;
      count++;
    }
  }
  return Math.sqrt(ssPred / count);
}

// Compute gap-closed %
function gapClosed(rms, rmsM0, rmsM2) {
  return 100 * (rmsM0 - rms) / (rmsM0 - rmsM2);
}

// ======================================================================
// PREPARE DATA
// ======================================================================
const Y = gals.map(g => g.logA0);
const hiDef = gals.map(g => g.hi_deficiency);

// Check for nulls
const validMask = hiDef.map(h => h !== null && !isNaN(h));
const nValid = validMask.filter(v => v).length;

// Baseline A variables: logMHI, rcWiggliness, envCode, logSigma0
const X_A = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0]);
// Baseline B: + logMeanRun
const X_B = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun]);
// A + hiDef
const X_A_hiDef = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.hi_deficiency]);
// B + hiDef
const X_B_hiDef = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun, g.hi_deficiency]);

// M0 and M2 reference RMS
const sdY = sd(Y);
const rmsM0 = sdY; // Universal a0
const rmsM2 = 0; // Per-galaxy (perfect)

console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  PHASE 57: HI DEFICIENCY                                                      ║');
console.log('║  Full formal test against frozen Baselines A and B                             ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

console.log('  Definition: DEF_HI = log(MHI_expected) - log(MHI_observed)');
console.log('  Expected MHI from L3.6-MHI scaling: log(MHI) = ' + 
  stageA.hiDefScaling.a.toFixed(3) + ' + ' + stageA.hiDefScaling.b.toFixed(3) + '*log(L36)');
console.log('  RMS of scaling relation: ' + stageA.hiDefScaling.rms.toFixed(3) + ' dex');
console.log('  DEF > 0 = less gas than expected (deficient)');
console.log('  DEF < 0 = more gas than expected (excess)');
console.log('  N galaxies: ' + nValid);
console.log('  HI deficiency range: [' + Math.min(...hiDef).toFixed(3) + ', ' + Math.max(...hiDef).toFixed(3) + ']');
console.log('  HI deficiency mean: ' + mean(hiDef).toFixed(4) + ', SD: ' + sd(hiDef).toFixed(4));
console.log();

// ======================================================================
// RAW CORRELATION
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  RAW CORRELATION: HI deficiency vs logA0');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const rawCorr = pearson(hiDef, Y);
console.log('  r = ' + rawCorr.r.toFixed(4));
console.log('  t = ' + rawCorr.t.toFixed(3));
console.log('  p = ' + rawCorr.p.toFixed(6));
console.log('  Direction: ' + (rawCorr.r > 0 ? 'POSITIVE (more deficient -> higher a0)' : 'NEGATIVE (more deficient -> lower a0)'));
console.log();

// Also check raw correlation with logMHI for comparison
const logMHI = gals.map(g => g.logMHI);
const corrMHI_hiDef = pearson(hiDef, logMHI);
console.log('  Cross-check: r(HIdef, logMHI) = ' + corrMHI_hiDef.r.toFixed(4));
console.log('  This tells us how much HIdef is just recoding logMHI.');
console.log();

// ======================================================================
// AFTER BASELINE A
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  AFTER BASELINE A (logMHI + wig + env + Sig0)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const fitA = ols(Y, X_A);
const residA = fitA.resid;
const corrAfterA = pearson(hiDef, residA);
console.log('  Baseline A: R2adj=' + fitA.r2adj.toFixed(4) + ', residSD=' + fitA.se.toFixed(4));
console.log('  r(HIdef, residA) = ' + corrAfterA.r.toFixed(4));
console.log('  t = ' + corrAfterA.t.toFixed(3));
console.log('  p = ' + corrAfterA.p.toFixed(6));
console.log();

// Full model A + hiDef
const fitA_hiDef = ols(Y, X_A_hiDef);
console.log('  Model A + HIdef:');
console.log('    R2adj = ' + fitA_hiDef.r2adj.toFixed(4) + ' (vs A: ' + fitA.r2adj.toFixed(4) + ')');
console.log('    residSD = ' + fitA_hiDef.se.toFixed(4) + ' (vs A: ' + fitA.se.toFixed(4) + ')');
console.log('    HIdef coeff = ' + fitA_hiDef.beta[5].toFixed(4) + ' (t=' + fitA_hiDef.tStats[5].toFixed(3) + ', p=' + fitA_hiDef.pVals[5].toFixed(4) + ')');
console.log();

// ======================================================================
// AFTER BASELINE B
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  AFTER BASELINE B (logMHI + wig + env + Sig0 + meanRun)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const fitB = ols(Y, X_B);
const residB = fitB.resid;
const corrAfterB = pearson(hiDef, residB);
console.log('  Baseline B: R2adj=' + fitB.r2adj.toFixed(4) + ', residSD=' + fitB.se.toFixed(4));
console.log('  r(HIdef, residB) = ' + corrAfterB.r.toFixed(4));
console.log('  t = ' + corrAfterB.t.toFixed(3));
console.log('  p = ' + corrAfterB.p.toFixed(6));
console.log();

// Full model B + hiDef
const fitB_hiDef = ols(Y, X_B_hiDef);
console.log('  Model B + HIdef:');
console.log('    R2adj = ' + fitB_hiDef.r2adj.toFixed(4) + ' (vs B: ' + fitB.r2adj.toFixed(4) + ')');
console.log('    residSD = ' + fitB_hiDef.se.toFixed(4) + ' (vs B: ' + fitB.se.toFixed(4) + ')');
console.log('    HIdef coeff = ' + fitB_hiDef.beta[6].toFixed(4) + ' (t=' + fitB_hiDef.tStats[6].toFixed(3) + ', p=' + fitB_hiDef.pVals[6].toFixed(4) + ')');
console.log();

// ======================================================================
// KEY TEST: AFTER logMHI ALONE
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  KEY TEST: HI deficiency after controlling for logMHI alone');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const fitMHI = ols(Y, gals.map(g => [g.logMHI]));
const residMHI = fitMHI.resid;
const corrAfterMHI = pearson(hiDef, residMHI);
console.log('  logMHI-only: R2adj=' + fitMHI.r2adj.toFixed(4));
console.log('  r(HIdef, resid|logMHI) = ' + corrAfterMHI.r.toFixed(4));
console.log('  t = ' + corrAfterMHI.t.toFixed(3));
console.log('  p = ' + corrAfterMHI.p.toFixed(6));
console.log('  INTERPRETATION: ' + (corrAfterMHI.p < 0.05 ? 
  'HIdef adds info BEYOND raw logMHI — not just a recoding!' : 
  'HIdef does NOT add info beyond logMHI — likely just a recoding.'));
console.log();

// Also test: logMHI after controlling for HIdef
const fitHIDef_only = ols(Y, gals.map(g => [g.hi_deficiency]));
const residHIDef = fitHIDef_only.resid;
const corrMHI_afterHIDef = pearson(logMHI, residHIDef);
console.log('  Reverse: r(logMHI, resid|HIdef) = ' + corrMHI_afterHIDef.r.toFixed(4));
console.log('  t = ' + corrMHI_afterHIDef.t.toFixed(3) + ', p = ' + corrMHI_afterHIDef.p.toFixed(6));
console.log('  This tells us if logMHI adds beyond HIdef.');
console.log();

// ======================================================================
// LOO CROSS-VALIDATION
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  LOO CROSS-VALIDATION');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const looA = looCV(Y, X_A);
const looA_hiDef = looCV(Y, X_A_hiDef);
const looB = looCV(Y, X_B);
const looB_hiDef = looCV(Y, X_B_hiDef);

const gcA = gapClosed(looA, sdY, 0);
const gcA_hiDef = gapClosed(looA_hiDef, sdY, 0);
const gcB = gapClosed(looB, sdY, 0);
const gcB_hiDef = gapClosed(looB_hiDef, sdY, 0);

console.log('  Baseline A LOO RMS: ' + looA.toFixed(4) + ' (' + gcA.toFixed(1) + '% gap closed)');
console.log('  A + HIdef LOO RMS:  ' + looA_hiDef.toFixed(4) + ' (' + gcA_hiDef.toFixed(1) + '% gap closed)');
console.log('  Delta:              ' + (gcA_hiDef - gcA).toFixed(1) + ' percentage points');
console.log('  Adds above A?      ' + (gcA_hiDef > gcA ? 'YES' : 'NO'));
console.log();
console.log('  Baseline B LOO RMS: ' + looB.toFixed(4) + ' (' + gcB.toFixed(1) + '% gap closed)');
console.log('  B + HIdef LOO RMS:  ' + looB_hiDef.toFixed(4) + ' (' + gcB_hiDef.toFixed(1) + '% gap closed)');
console.log('  Delta:              ' + (gcB_hiDef - gcB).toFixed(1) + ' percentage points');
console.log('  Adds above B?      ' + (gcB_hiDef > gcB ? 'YES' : 'NO'));
console.log();

// ======================================================================
// PERMUTATION TEST
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PERMUTATION TEST (5000 permutations)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const permP_A = permTest(Y, X_A_hiDef, 4, 5000);  // HIdef is index 4 in A+hiDef
const permP_B = permTest(Y, X_B_hiDef, 5, 5000);  // HIdef is index 5 in B+hiDef

console.log('  Permutation p (A + HIdef): ' + permP_A.toFixed(4));
console.log('  Permutation p (B + HIdef): ' + permP_B.toFixed(4));
console.log();

// ======================================================================
// BOOTSTRAP CI
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  BOOTSTRAP 95% CI (2000 resamples)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const bootA = bootstrapCI(Y, X_A_hiDef, 4, 2000);
const bootB = bootstrapCI(Y, X_B_hiDef, 5, 2000);

console.log('  Bootstrap CI (in A+HIdef): [' + bootA.lo.toFixed(4) + ', ' + bootA.hi.toFixed(4) + '] median=' + bootA.median.toFixed(4));
console.log('  Includes zero? ' + (bootA.lo <= 0 && bootA.hi >= 0 ? 'YES' : 'NO'));
console.log('  Bootstrap CI (in B+HIdef): [' + bootB.lo.toFixed(4) + ', ' + bootB.hi.toFixed(4) + '] median=' + bootB.median.toFixed(4));
console.log('  Includes zero? ' + (bootB.lo <= 0 && bootB.hi >= 0 ? 'YES' : 'NO'));
console.log();

// ======================================================================
// JACKKNIFE SIGN FLIPS
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  JACKKNIFE SIGN STABILITY');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const jkA = jackknifeSigns(Y, X_A_hiDef, 4);
const jkB = jackknifeSigns(Y, X_B_hiDef, 5);

console.log('  Sign flips (A + HIdef): ' + jkA + '/' + N + ' (' + (100*jkA/N).toFixed(1) + '%)');
console.log('  Sign flips (B + HIdef): ' + jkB + '/' + N + ' (' + (100*jkB/N).toFixed(1) + '%)');
console.log();

// ======================================================================
// EXTRA CONTROLS
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  EXTRA CONTROLS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// After logMHI only
console.log('  |logMHI: r(HIdef, resid|logMHI) = ' + corrAfterMHI.r.toFixed(4) + ', t=' + corrAfterMHI.t.toFixed(3) + ', p=' + corrAfterMHI.p.toFixed(4));

// After envCode only
const fitEnv = ols(Y, gals.map(g => [g.envCode]));
const corrAfterEnv = pearson(hiDef, fitEnv.resid);
console.log('  |envCode: r(HIdef, resid|envCode) = ' + corrAfterEnv.r.toFixed(4) + ', t=' + corrAfterEnv.t.toFixed(3) + ', p=' + corrAfterEnv.p.toFixed(4));

// After Sigma0 only
const fitSig = ols(Y, gals.map(g => [g.logSigma0]));
const corrAfterSig = pearson(hiDef, fitSig.resid);
console.log('  |Sigma0: r(HIdef, resid|Sigma0) = ' + corrAfterSig.r.toFixed(4) + ', t=' + corrAfterSig.t.toFixed(3) + ', p=' + corrAfterSig.p.toFixed(4));

// After logMHI + envCode jointly
const fitME = ols(Y, gals.map(g => [g.logMHI, g.envCode]));
const corrAfterME = pearson(hiDef, fitME.resid);
console.log('  |logMHI+envCode: r(HIdef, resid) = ' + corrAfterME.r.toFixed(4) + ', t=' + corrAfterME.t.toFixed(3) + ', p=' + corrAfterME.p.toFixed(4));
console.log();

// ======================================================================
// FINAL VERDICT
// ======================================================================
console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  FINAL VERDICT                                                                ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

// Determine verdict
const passesRaw = rawCorr.p < 0.05;
const passesAfterA = fitA_hiDef.pVals[5] < 0.05;
const passesAfterB = fitB_hiDef.pVals[6] < 0.05;
const passesAfterMHI = corrAfterMHI.p < 0.05;
const addsAboveA = gcA_hiDef > gcA;
const addsAboveB = gcB_hiDef > gcB;
const permPassA = permP_A < 0.05;
const permPassB = permP_B < 0.05;
const bootExclZeroA = !(bootA.lo <= 0 && bootA.hi >= 0);
const bootExclZeroB = !(bootB.lo <= 0 && bootB.hi >= 0);

let verdict;
if (passesAfterB && addsAboveB && permPassB && bootExclZeroB) {
  verdict = 'CONFIRMED';
} else if (passesAfterA && addsAboveA) {
  verdict = 'PARTIAL';
} else {
  verdict = 'FAIL';
}

console.log('  Raw: r=' + rawCorr.r.toFixed(4) + ', t=' + rawCorr.t.toFixed(3) + ', p=' + rawCorr.p.toFixed(4) + ' — ' + (passesRaw ? 'PASS' : 'FAIL'));
console.log('  |logMHI alone: r=' + corrAfterMHI.r.toFixed(4) + ', p=' + corrAfterMHI.p.toFixed(4) + ' — ' + (passesAfterMHI ? 'SURVIVES logMHI' : 'ABSORBED by logMHI'));
console.log('  |Baseline A: t=' + fitA_hiDef.tStats[5].toFixed(3) + ', p=' + fitA_hiDef.pVals[5].toFixed(4) + ' — ' + (passesAfterA ? 'PASS' : 'FAIL'));
console.log('  |Baseline B: t=' + fitB_hiDef.tStats[6].toFixed(3) + ', p=' + fitB_hiDef.pVals[6].toFixed(4) + ' — ' + (passesAfterB ? 'PASS' : 'FAIL'));
console.log('  LOO above A: ' + (addsAboveA ? 'YES (+' + (gcA_hiDef-gcA).toFixed(1) + 'pp)' : 'NO'));
console.log('  LOO above B: ' + (addsAboveB ? 'YES (+' + (gcB_hiDef-gcB).toFixed(1) + 'pp)' : 'NO'));
console.log('  Perm p (A): ' + permP_A.toFixed(4) + ' — ' + (permPassA ? 'PASS' : 'FAIL'));
console.log('  Perm p (B): ' + permP_B.toFixed(4) + ' — ' + (permPassB ? 'PASS' : 'FAIL'));
console.log('  Boot CI (A): [' + bootA.lo.toFixed(4) + ', ' + bootA.hi.toFixed(4) + '] — ' + (bootExclZeroA ? 'excludes zero' : 'includes zero'));
console.log('  Boot CI (B): [' + bootB.lo.toFixed(4) + ', ' + bootB.hi.toFixed(4) + '] — ' + (bootExclZeroB ? 'excludes zero' : 'includes zero'));
console.log('  JK sign flips A: ' + jkA + '/' + N);
console.log('  JK sign flips B: ' + jkB + '/' + N);
console.log();
console.log('  ═══════════════════════════════════════════════');
console.log('  VERDICT: ' + verdict);
console.log('  ═══════════════════════════════════════════════');
console.log();

if (verdict === 'FAIL') {
  console.log('  INTERPRETATION:');
  if (!passesAfterMHI) {
    console.log('  HI deficiency does NOT add information beyond raw logMHI.');
    console.log('  It is a RECODING of gas mass, not a new physical dimension.');
    console.log('  The question "how much gas relative to expectations" is');
    console.log('  already answered by "how much gas in absolute terms".');
  } else if (!passesAfterA) {
    console.log('  HI deficiency adds beyond logMHI alone, but is absorbed');
    console.log('  by the combination of baseline A variables.');
  } else if (!passesAfterB) {
    console.log('  HI deficiency adds beyond baseline A, but is absorbed');
    console.log('  by baseline B (including meanRun).');
  }
} else if (verdict === 'PARTIAL') {
  console.log('  INTERPRETATION:');
  console.log('  HI deficiency adds above baseline A but NOT above B.');
  console.log('  Partial information, likely absorbed by meanRun or other B vars.');
} else {
  console.log('  INTERPRETATION:');
  console.log('  HI deficiency adds NEW information above both baselines!');
  console.log('  Gas DEFICIENCY (history) adds beyond raw gas MASS.');
  console.log('  This opens the door to environmental processing tests.');
}

// Save results
const output = {
  phase: 57,
  program: 'New Program',
  variable: 'HI deficiency',
  definition: 'DEF = log(MHI_expected) - log(MHI_observed)',
  scaling: stageA.hiDefScaling,
  nGalaxies: nValid,
  raw: { r: rawCorr.r, t: rawCorr.t, p: rawCorr.p },
  afterLogMHI: { r: corrAfterMHI.r, t: corrAfterMHI.t, p: corrAfterMHI.p },
  afterBaselineA: { r: corrAfterA.r, t: corrAfterA.t, p: corrAfterA.p, 
    coeff: fitA_hiDef.beta[5], tStat: fitA_hiDef.tStats[5], pVal: fitA_hiDef.pVals[5] },
  afterBaselineB: { r: corrAfterB.r, t: corrAfterB.t, p: corrAfterB.p,
    coeff: fitB_hiDef.beta[6], tStat: fitB_hiDef.tStats[6], pVal: fitB_hiDef.pVals[6] },
  loo: {
    A: { rms: looA, gapClosed: gcA },
    A_hiDef: { rms: looA_hiDef, gapClosed: gcA_hiDef },
    B: { rms: looB, gapClosed: gcB },
    B_hiDef: { rms: looB_hiDef, gapClosed: gcB_hiDef }
  },
  permutation: { pA: permP_A, pB: permP_B },
  bootstrap: { A: bootA, B: bootB },
  jackknife: { flipsA: jkA, flipsB: jkB },
  extraControls: {
    afterLogMHI: { r: corrAfterMHI.r, p: corrAfterMHI.p },
    afterEnvCode: { r: corrAfterEnv.r, p: corrAfterEnv.p },
    afterSigma0: { r: corrAfterSig.r, p: corrAfterSig.p },
    afterMHI_env: { r: corrAfterME.r, p: corrAfterME.p }
  },
  verdict: verdict
};

fs.writeFileSync('public/phase57-hi-deficiency.json', JSON.stringify(output, null, 2));
console.log('  Results saved to public/phase57-hi-deficiency.json');
