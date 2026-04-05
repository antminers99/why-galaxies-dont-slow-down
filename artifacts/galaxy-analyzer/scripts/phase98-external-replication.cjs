#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 98: EXTERNAL REPLICATION OF THE DIRECT LAW');
console.log('  Does log(Vflat/InnerVmax) → outerSlope replicate on N=45?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function olsRegress(X, y) {
  const n = y.length, p = X[0].length;
  if (n <= p + 1) return null;
  const Xt = Array.from({ length: p }, (_, i) => X.map(row => row[i]));
  const XtX = Array.from({ length: p }, (_, i) =>
    Array.from({ length: p }, (_, j) =>
      Xt[i].reduce((s, _, k) => s + Xt[i][k] * Xt[j][k], 0)));
  const Xty = Xt.map(col => col.reduce((s, v, k) => s + v * y[k], 0));
  const aug = XtX.map((row, i) => [...row, Xty[i]]);
  for (let i = 0; i < p; i++) {
    let maxRow = i;
    for (let k = i + 1; k < p; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) maxRow = k;
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];
    if (Math.abs(aug[i][i]) < 1e-12) return null;
    for (let k = 0; k < p; k++) {
      if (k === i) continue;
      const f = aug[k][i] / aug[i][i];
      for (let j = i; j <= p; j++) aug[k][j] -= f * aug[i][j];
    }
  }
  const beta = aug.map((row, i) => row[p] / row[i]);
  const yhat = X.map(row => row.reduce((s, v, j) => s + v * beta[j], 0));
  const residuals = y.map((v, i) => v - yhat[i]);
  const ss_res = residuals.reduce((s, r) => s + r * r, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, residuals, yhat };
}

function looR2(X_fn, y_fn, data) {
  const n = data.length;
  let ss_res = 0, ss_tot = 0;
  const yAll = data.map(y_fn);
  const yMean = mean(yAll);
  for (let i = 0; i < n; i++) {
    const train = data.filter((_, j) => j !== i);
    const fit = olsRegress(train.map(X_fn), train.map(y_fn));
    if (!fit) return NaN;
    const yPred = X_fn(data[i]).reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss_res += (yAll[i] - yPred) ** 2;
    ss_tot += (yAll[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
}

const n45csv = fs.readFileSync(path.join(__dirname, '..', 'public', 'replication', 'N45_final_dataset.csv'), 'utf-8').trim().split('\n');
const n45header = n45csv[0].split(',');
const n45names = new Set();
for (let i = 1; i < n45csv.length; i++) {
  const cols = n45csv[i].split(',');
  n45names.add(cols[0]);
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { T, L36, Rdisk, MHI, Vflat, Q };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul });
}

function buildGalaxy(name, rcPoints, t1) {
  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    pts.push({ r: pt.rad, vobs: pt.vobs });
  }
  if (pts.length < 8) return null;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const innerPts = sorted.slice(0, half);
  const outerPts = sorted.slice(half);

  const innerVmax = Math.max(...innerPts.map(p => p.vobs));

  const logR = outerPts.map(p => Math.log10(p.r));
  const logV = outerPts.map(p => Math.log10(p.vobs));
  const mr = mean(logR), mv = mean(logV);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mr) * (logV[i] - mv); sxx += (logR[i] - mr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logRatio = Math.log10(t1.Vflat > 0 && innerVmax > 0 ? t1.Vflat / innerVmax : 1);

  return {
    name, nPts: pts.length,
    innerVmax, logInnerVmax: Math.log10(innerVmax > 0 ? innerVmax : 1),
    outerSlope, logRatio,
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    Vflat: t1.Vflat, fgas,
    logMHI: Math.log10(t1.MHI),
    T: t1.T, Q: t1.Q,
    isN45: n45names.has(name),
  };
}

const allGalaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0) continue;
  const g = buildGalaxy(name, rcPoints, t1);
  if (g) allGalaxies.push(g);
}

const n45hi = allGalaxies.filter(g => g.isN45 && g.Vflat >= VFLAT_BREAK);
const nonN45hi = allGalaxies.filter(g => !g.isN45 && g.Vflat >= VFLAT_BREAK);
const allHi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);

console.log('  N45 high-regime: N=' + n45hi.length);
console.log('  Non-N45 high-regime: N=' + nonN45hi.length);
console.log('  All high-regime: N=' + allHi.length + '\n');

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: TRAIN ON NON-N45, TEST ON N45');
console.log('  (True external validation — N45 never seen during training)');
console.log('══════════════════════════════════════════════════════════════\n');

const X_ratio = g => [1, g.logRatio];
const X_vmax = g => [1, g.logInnerVmax];
const X_full = g => [1, g.logInnerVmax, g.logSigma0, g.logVflat, g.fgas, g.logMHI];
const y_fn = g => g.outerSlope;

function trainTest(X_fn, trainSet, testSet, label) {
  const fitTrain = olsRegress(trainSet.map(X_fn), trainSet.map(y_fn));
  if (!fitTrain) { console.log('  ' + label + ': FAILED to fit'); return null; }

  const yTest = testSet.map(y_fn);
  const yPred = testSet.map(g => X_fn(g).reduce((s, v, j) => s + v * fitTrain.beta[j], 0));
  const yTestMean = mean(yTest);
  const ss_res = yTest.reduce((s, v, i) => s + (v - yPred[i]) ** 2, 0);
  const ss_tot = yTest.reduce((s, v) => s + (v - yTestMean) ** 2, 0);
  const testR2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0;

  const rPred = pearsonR(yPred, yTest);
  const rmse = Math.sqrt(ss_res / yTest.length);

  console.log('  ' + label + ':');
  console.log('    Train R² = ' + fitTrain.r2.toFixed(3) + ' (N=' + trainSet.length + ')');
  console.log('    Test R² = ' + testR2.toFixed(3) + ' (N=' + testSet.length + ')');
  console.log('    r(pred, actual) = ' + rPred.toFixed(3));
  console.log('    RMSE = ' + rmse.toFixed(4));
  console.log();

  return { label, trainR2: +fitTrain.r2.toFixed(3), testR2: +testR2.toFixed(3), rPred: +rPred.toFixed(3), rmse: +rmse.toFixed(4), beta: fitTrain.beta.map(b => +b.toFixed(4)) };
}

const r1 = trainTest(X_ratio, nonN45hi, n45hi, 'log(Vflat/InnerVmax) only');
const r2 = trainTest(X_vmax, nonN45hi, n45hi, 'logInnerVmax only');
const r3 = trainTest(X_full, nonN45hi, n45hi, 'Structural + logInnerVmax');

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: TRAIN ON N45, TEST ON NON-N45 (reverse)');
console.log('══════════════════════════════════════════════════════════════\n');

const r4 = trainTest(X_ratio, n45hi, nonN45hi, 'log(Vflat/InnerVmax) only');
const r5 = trainTest(X_vmax, n45hi, nonN45hi, 'logInnerVmax only');
const r6 = trainTest(X_full, n45hi, nonN45hi, 'Structural + logInnerVmax');

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: WITHIN-N45 LOO (gold-standard subset)');
console.log('══════════════════════════════════════════════════════════════\n');

const looN45ratio = looR2(X_ratio, y_fn, n45hi);
const looN45vmax = looR2(X_vmax, y_fn, n45hi);
const looN45full = looR2(X_full, y_fn, n45hi);

console.log('  N45 LOO R² (log ratio only):        ' + looN45ratio.toFixed(3));
console.log('  N45 LOO R² (logInnerVmax only):      ' + looN45vmax.toFixed(3));
console.log('  N45 LOO R² (structural + innerVmax): ' + looN45full.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: RAW CORRELATIONS WITHIN N45');
console.log('══════════════════════════════════════════════════════════════\n');

const rN45ratio = pearsonR(n45hi.map(g => g.logRatio), n45hi.map(g => g.outerSlope));
const rN45vmax = pearsonR(n45hi.map(g => g.logInnerVmax), n45hi.map(g => g.outerSlope));
console.log('  r(logRatio, outerSlope) within N45 = ' + rN45ratio.toFixed(3));
console.log('  r(logInnerVmax, outerSlope) within N45 = ' + rN45vmax.toFixed(3));

const rNonRatio = pearsonR(nonN45hi.map(g => g.logRatio), nonN45hi.map(g => g.outerSlope));
const rNonVmax = pearsonR(nonN45hi.map(g => g.logInnerVmax), nonN45hi.map(g => g.outerSlope));
console.log('  r(logRatio, outerSlope) within non-N45 = ' + rNonRatio.toFixed(3));
console.log('  r(logInnerVmax, outerSlope) within non-N45 = ' + rNonVmax.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: Q=1 ONLY within N45 (highest quality)');
console.log('══════════════════════════════════════════════════════════════\n');

const n45q1hi = n45hi.filter(g => g.Q === 1);
console.log('  N45 Q=1 high-regime: N=' + n45q1hi.length);

if (n45q1hi.length >= 8) {
  const rQ1ratio = pearsonR(n45q1hi.map(g => g.logRatio), n45q1hi.map(g => g.outerSlope));
  const rQ1vmax = pearsonR(n45q1hi.map(g => g.logInnerVmax), n45q1hi.map(g => g.outerSlope));
  console.log('  r(logRatio, outerSlope) Q=1 N45 = ' + rQ1ratio.toFixed(3));
  console.log('  r(logInnerVmax, outerSlope) Q=1 N45 = ' + rQ1vmax.toFixed(3));

  const looQ1ratio = looR2(X_ratio, y_fn, n45q1hi);
  const looQ1vmax = looR2(X_vmax, y_fn, n45q1hi);
  console.log('  LOO R² (logRatio) Q=1 N45: ' + looQ1ratio.toFixed(3));
  console.log('  LOO R² (logInnerVmax) Q=1 N45: ' + looQ1vmax.toFixed(3));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 6: COEFFICIENT STABILITY');
console.log('  Are the slopes consistent across N45 vs non-N45?');
console.log('══════════════════════════════════════════════════════════════\n');

const fitN45 = olsRegress(n45hi.map(X_ratio), n45hi.map(y_fn));
const fitNon = olsRegress(nonN45hi.map(X_ratio), nonN45hi.map(y_fn));
const fitAll = olsRegress(allHi.map(X_ratio), allHi.map(y_fn));

console.log('  Equation: outerSlope = a + b * log(Vflat/InnerVmax)');
console.log('    N45:     a=' + fitN45.beta[0].toFixed(3) + ', b=' + fitN45.beta[1].toFixed(3));
console.log('    Non-N45: a=' + fitNon.beta[0].toFixed(3) + ', b=' + fitNon.beta[1].toFixed(3));
console.log('    All:     a=' + fitAll.beta[0].toFixed(3) + ', b=' + fitAll.beta[1].toFixed(3));

const slopeDiff = Math.abs(fitN45.beta[1] - fitNon.beta[1]);
const slopeAvg = (Math.abs(fitN45.beta[1]) + Math.abs(fitNon.beta[1])) / 2;
const slopePctDiff = (slopeDiff / slopeAvg * 100).toFixed(1);
console.log('  Slope difference: ' + slopePctDiff + '%');

if (parseFloat(slopePctDiff) < 30) {
  console.log('  → Coefficients are STABLE across samples');
} else {
  console.log('  → Coefficients show INSTABILITY across samples');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 7: PERMUTATION on N45 alone');
console.log('══════════════════════════════════════════════════════════════\n');

const realRn45 = pearsonR(n45hi.map(g => g.logRatio), n45hi.map(g => g.outerSlope));
const nPerm = 5000;
let nBetter = 0;
const yN45 = n45hi.map(g => g.outerSlope);
for (let i = 0; i < nPerm; i++) {
  const shuffled = [...yN45].sort(() => Math.random() - 0.5);
  const pr = pearsonR(n45hi.map(g => g.logRatio), shuffled);
  if (Math.abs(pr) >= Math.abs(realRn45)) nBetter++;
}
console.log('  r(logRatio, outerSlope) in N45 = ' + realRn45.toFixed(3) + ', perm p = ' + (nBetter / nPerm).toFixed(4));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  FINAL SUMMARY');
console.log('══════════════════════════════════════════════════════════════\n');

const extTestR2 = r1 ? r1.testR2 : NaN;
const revTestR2 = r4 ? r4.testR2 : NaN;

console.log('  External replication (train non-N45, test N45):');
console.log('    logRatio:            test R² = ' + (r1 ? r1.testR2 : 'N/A'));
console.log('    logInnerVmax:        test R² = ' + (r2 ? r2.testR2 : 'N/A'));
console.log('    Structural+innerVmax: test R² = ' + (r3 ? r3.testR2 : 'N/A'));
console.log('  Reverse (train N45, test non-N45):');
console.log('    logRatio:            test R² = ' + (r4 ? r4.testR2 : 'N/A'));
console.log('    logInnerVmax:        test R² = ' + (r5 ? r5.testR2 : 'N/A'));
console.log('    Structural+innerVmax: test R² = ' + (r6 ? r6.testR2 : 'N/A'));
console.log('  Within-N45 LOO R²:');
console.log('    logRatio:            ' + looN45ratio.toFixed(3));
console.log('    logInnerVmax:        ' + looN45vmax.toFixed(3));

let verdict;
if (extTestR2 > 0.4 && revTestR2 > 0.4 && parseFloat(slopePctDiff) < 30) {
  verdict = 'REPLICATED';
  console.log('\n  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  EXTERNALLY REPLICATED: The direct law holds on N45        ║');
  console.log('  ║  log(Vflat/InnerVmax) → outerSlope is a robust physical   ║');
  console.log('  ║  law that generalizes across galaxy samples                ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else if (extTestR2 > 0.2) {
  verdict = 'PARTIALLY-REPLICATED';
  console.log('\n  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  PARTIALLY REPLICATED: signal present but weaker on N45    ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else {
  verdict = 'FAILED';
  console.log('\n  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  REPLICATION FAILED: law does not transfer to N45          ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
}

const output = {
  phase: '98',
  title: 'External Replication of the Direct Law',
  n45hi: n45hi.length,
  nonN45hi: nonN45hi.length,
  trainOnNon_testOnN45: { ratio: r1, vmax: r2, full: r3 },
  trainOnN45_testOnNon: { ratio: r4, vmax: r5, full: r6 },
  withinN45_LOO: { ratio: +looN45ratio.toFixed(3), vmax: +looN45vmax.toFixed(3), full: +looN45full.toFixed(3) },
  rawCorr: {
    n45_ratio: +rN45ratio.toFixed(3), n45_vmax: +rN45vmax.toFixed(3),
    non_ratio: +rNonRatio.toFixed(3), non_vmax: +rNonVmax.toFixed(3),
  },
  coeffStability: {
    n45_slope: +fitN45.beta[1].toFixed(3),
    non_slope: +fitNon.beta[1].toFixed(3),
    all_slope: +fitAll.beta[1].toFixed(3),
    pctDiff: +slopePctDiff,
  },
  permN45: { r: +realRn45.toFixed(3), p: +(nBetter / nPerm).toFixed(4) },
  verdict,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase98-external-replication.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase98-external-replication.json');
