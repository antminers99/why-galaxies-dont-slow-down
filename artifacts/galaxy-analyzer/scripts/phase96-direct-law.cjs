#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 96: DIRECT DYNAMICAL LAW WITHOUT a₀');
console.log('  Can we predict outerSlope from inner dynamics + structure?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0fromPts(pts) {
  if (pts.length < 3) return { logA0: NaN, a0: NaN };
  let lo = 0.5, hi = 5.0;
  for (let s = 0; s < 200; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gb = Math.pow(10, p.logGbar);
      c1 += (p.logGobs - Math.log10(mcgaughRAR(gb, Math.pow(10, m1)))) ** 2;
      c2 += (p.logGobs - Math.log10(mcgaughRAR(gb, Math.pow(10, m2)))) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  return { logA0: (lo + hi) / 2, a0: Math.pow(10, (lo + hi) / 2) };
}

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
    const Xtrain = train.map(X_fn);
    const ytrain = train.map(y_fn);
    const fit = olsRegress(Xtrain, ytrain);
    if (!fit) return NaN;
    const xTest = X_fn(data[i]);
    const yPred = xTest.reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss_res += (yAll[i] - yPred) ** 2;
    ss_tot += (yAll[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
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

const allGalaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    const gObs = pt.vobs * pt.vobs / pt.rad;
    const gBar = Math.abs(vBarSq) / pt.rad;
    if (gBar <= 0 || gObs <= 0 || !isFinite(gBar) || !isFinite(gObs)) continue;
    pts.push({ r: pt.rad, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), vobs: pt.vobs, vgas: pt.vgas, vdisk: pt.vdisk, vbul: pt.vbul || 0 });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const fitAll = fitA0fromPts(sorted);
  if (!isFinite(fitAll.logA0) || fitAll.logA0 < 0.5 || fitAll.logA0 > 5.5) continue;

  const half = Math.floor(sorted.length / 2);
  const innerPts = sorted.slice(0, half);
  const outerPts = sorted.slice(half);

  const fitInner = fitA0fromPts(innerPts);

  const logR_outer = outerPts.map(p => Math.log10(p.r));
  const logV_outer = outerPts.map(p => Math.log10(p.vobs));
  const mr = mean(logR_outer), mv = mean(logV_outer);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR_outer.length; i++) { sxy += (logR_outer[i] - mr) * (logV_outer[i] - mv); sxx += (logR_outer[i] - mr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const innerVobs = innerPts.map(p => p.vobs);
  const innerVmax = Math.max(...innerVobs);
  const logInnerVmax = Math.log10(innerVmax > 0 ? innerVmax : 1);

  const innerDiscrepancy = mean(innerPts.map(p => p.logGobs - p.logGbar));

  const innerCurvature = (() => {
    if (innerPts.length < 5) return NaN;
    const x = innerPts.map(p => Math.log10(p.r));
    const y = innerPts.map(p => Math.log10(p.vobs));
    const mx2 = mean(x);
    const X = x.map(xi => [1, xi - mx2, (xi - mx2) ** 2]);
    const fit2 = olsRegress(X, y);
    return fit2 ? fit2.beta[2] : NaN;
  })();

  const innerBaryonFrac = mean(innerPts.map(p => {
    const vBarSq2 = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) + UPSILON_BULGE * p.vbul * Math.abs(p.vbul) + p.vgas * Math.abs(p.vgas);
    return p.vobs > 0 ? Math.abs(vBarSq2) / (p.vobs * p.vobs) : NaN;
  }).filter(v => isFinite(v)));

  const innerSlope = (() => {
    const x = innerPts.map(p => Math.log10(p.r));
    const y = innerPts.map(p => Math.log10(p.vobs));
    const mx3 = mean(x), my3 = mean(y);
    let s3 = 0, sx3 = 0;
    for (let i = 0; i < x.length; i++) { s3 += (x[i] - mx3) * (y[i] - my3); sx3 += (x[i] - mx3) ** 2; }
    return sx3 > 0 ? s3 / sx3 : 0;
  })();

  const innerConc = (() => {
    if (innerPts.length < 4) return NaN;
    const rHalf = innerPts[Math.floor(innerPts.length / 2)].r;
    const rLast = innerPts[innerPts.length - 1].r;
    const vHalf = innerPts[Math.floor(innerPts.length / 2)].vobs;
    const vLast = innerPts[innerPts.length - 1].vobs;
    return rLast > 0 && rHalf > 0 ? (vLast * vLast * rLast) / (vHalf * vHalf * rHalf) : NaN;
  })();

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);

  allGalaxies.push({
    name, nPts: pts.length,
    logA0all: fitAll.logA0,
    logA0inner: fitInner.logA0,
    outerSlope,
    logInnerVmax,
    innerDiscrepancy: isFinite(innerDiscrepancy) ? innerDiscrepancy : NaN,
    innerCurvature: isFinite(innerCurvature) ? innerCurvature : NaN,
    innerBaryonFrac: isFinite(innerBaryonFrac) ? innerBaryonFrac : NaN,
    innerSlope,
    innerConc: isFinite(innerConc) ? innerConc : NaN,
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    Vflat: t1.Vflat, fgas,
    logMHI: Math.log10(t1.MHI),
    logL36: Math.log10(t1.L36),
    T: t1.T, Q: t1.Q,
  });
}

const hi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK &&
  isFinite(g.logA0inner) && isFinite(g.innerDiscrepancy) && isFinite(g.innerCurvature));
console.log('  High-regime valid: N=' + hi.length + '\n');

console.log('══════════════════════════════════════════════════════════════');
console.log('  MODEL 1: STRUCTURAL ONLY (no inner dynamics, no a₀)');
console.log('  outerSlope ~ logΣ₀ + logVflat + fgas + logMHI');
console.log('══════════════════════════════════════════════════════════════\n');

const M1_X = g => [1, g.logSigma0, g.logVflat, g.fgas, g.logMHI];
const y_fn = g => g.outerSlope;
const fitM1 = olsRegress(hi.map(M1_X), hi.map(y_fn));
const looM1 = looR2(M1_X, y_fn, hi);
console.log('  R² = ' + fitM1.r2.toFixed(3) + ', LOO R² = ' + looM1.toFixed(3));
console.log('  Gap = ' + ((1 - looM1 / fitM1.r2) * 100).toFixed(1) + '%');
console.log('  Coefficients: intercept=' + fitM1.beta[0].toFixed(3) + ', logΣ₀=' + fitM1.beta[1].toFixed(3) +
  ', logVflat=' + fitM1.beta[2].toFixed(3) + ', fgas=' + fitM1.beta[3].toFixed(3) + ', logMHI=' + fitM1.beta[4].toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  MODEL 2: STRUCTURAL + a₀ (the old approach)');
console.log('  outerSlope ~ logΣ₀ + logVflat + fgas + logMHI + logA0');
console.log('══════════════════════════════════════════════════════════════\n');

const M2_X = g => [1, g.logSigma0, g.logVflat, g.fgas, g.logMHI, g.logA0all];
const fitM2 = olsRegress(hi.map(M2_X), hi.map(y_fn));
const looM2 = looR2(M2_X, y_fn, hi);
console.log('  R² = ' + fitM2.r2.toFixed(3) + ', LOO R² = ' + looM2.toFixed(3));
console.log('  Gap = ' + ((1 - looM2 / fitM2.r2) * 100).toFixed(1) + '%');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  MODEL 3: STRUCTURAL + logInnerVmax (the new approach)');
console.log('  outerSlope ~ logΣ₀ + logVflat + fgas + logMHI + logInnerVmax');
console.log('══════════════════════════════════════════════════════════════\n');

const M3_X = g => [1, g.logSigma0, g.logVflat, g.fgas, g.logMHI, g.logInnerVmax];
const fitM3 = olsRegress(hi.map(M3_X), hi.map(y_fn));
const looM3 = looR2(M3_X, y_fn, hi);
console.log('  R² = ' + fitM3.r2.toFixed(3) + ', LOO R² = ' + looM3.toFixed(3));
console.log('  Gap = ' + ((1 - looM3 / fitM3.r2) * 100).toFixed(1) + '%');
console.log('  Coefficients: intercept=' + fitM3.beta[0].toFixed(3) + ', logΣ₀=' + fitM3.beta[1].toFixed(3) +
  ', logVflat=' + fitM3.beta[2].toFixed(3) + ', fgas=' + fitM3.beta[3].toFixed(3) + ', logMHI=' + fitM3.beta[4].toFixed(3) +
  ', logInnerVmax=' + fitM3.beta[5].toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  MODEL 4: FULL INNER DYNAMICS MODEL');
console.log('  outerSlope ~ logΣ₀ + logVflat + fgas + logMHI + logInnerVmax');
console.log('              + innerDiscrepancy + innerCurvature');
console.log('══════════════════════════════════════════════════════════════\n');

const M4_X = g => [1, g.logSigma0, g.logVflat, g.fgas, g.logMHI, g.logInnerVmax, g.innerDiscrepancy, g.innerCurvature];
const fitM4 = olsRegress(hi.map(M4_X), hi.map(y_fn));
const looM4 = looR2(M4_X, y_fn, hi);
console.log('  R² = ' + fitM4.r2.toFixed(3) + ', LOO R² = ' + looM4.toFixed(3));
console.log('  Gap = ' + ((1 - looM4 / fitM4.r2) * 100).toFixed(1) + '%');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  MODEL 5: INNER DYNAMICS ONLY (no structural catalog vars)');
console.log('  outerSlope ~ logInnerVmax + innerDiscrepancy + innerCurvature');
console.log('══════════════════════════════════════════════════════════════\n');

const M5_X = g => [1, g.logInnerVmax, g.innerDiscrepancy, g.innerCurvature];
const fitM5 = olsRegress(hi.map(M5_X), hi.map(y_fn));
const looM5 = looR2(M5_X, y_fn, hi);
console.log('  R² = ' + fitM5.r2.toFixed(3) + ', LOO R² = ' + looM5.toFixed(3));
console.log('  Gap = ' + ((1 - looM5 / fitM5.r2) * 100).toFixed(1) + '%');
console.log('  Coefficients: intercept=' + fitM5.beta[0].toFixed(3) + ', logInnerVmax=' + fitM5.beta[1].toFixed(3) +
  ', innerDisc=' + fitM5.beta[2].toFixed(3) + ', innerCurv=' + fitM5.beta[3].toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  MODEL 6: MINIMAL — logInnerVmax ONLY');
console.log('  outerSlope ~ logInnerVmax');
console.log('══════════════════════════════════════════════════════════════\n');

const M6_X = g => [1, g.logInnerVmax];
const fitM6 = olsRegress(hi.map(M6_X), hi.map(y_fn));
const looM6 = looR2(M6_X, y_fn, hi);
console.log('  R² = ' + fitM6.r2.toFixed(3) + ', LOO R² = ' + looM6.toFixed(3));
console.log('  Gap = ' + ((1 - looM6 / fitM6.r2) * 100).toFixed(1) + '%');
console.log('  outerSlope = ' + fitM6.beta[0].toFixed(3) + ' + ' + fitM6.beta[1].toFixed(3) + ' * logInnerVmax');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  MODEL 7: a₀(inner) ONLY for comparison');
console.log('  outerSlope ~ logA0_inner');
console.log('══════════════════════════════════════════════════════════════\n');

const M7_X = g => [1, g.logA0inner];
const fitM7 = olsRegress(hi.map(M7_X), hi.map(y_fn));
const looM7 = looR2(M7_X, y_fn, hi);
console.log('  R² = ' + fitM7.r2.toFixed(3) + ', LOO R² = ' + looM7.toFixed(3));
console.log('  Gap = ' + ((1 - looM7 / fitM7.r2) * 100).toFixed(1) + '%');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  MODEL COMPARISON SUMMARY');
console.log('══════════════════════════════════════════════════════════════\n');

const models = [
  { name: 'M1: Structural only', r2: fitM1.r2, loo: looM1, p: 4 },
  { name: 'M2: Structural + a₀', r2: fitM2.r2, loo: looM2, p: 5 },
  { name: 'M3: Structural + innerVmax', r2: fitM3.r2, loo: looM3, p: 5 },
  { name: 'M4: Structural + inner dynamics', r2: fitM4.r2, loo: looM4, p: 7 },
  { name: 'M5: Inner dynamics only', r2: fitM5.r2, loo: looM5, p: 3 },
  { name: 'M6: logInnerVmax only', r2: fitM6.r2, loo: looM6, p: 1 },
  { name: 'M7: logA0_inner only', r2: fitM7.r2, loo: looM7, p: 1 },
];

console.log('  ' + 'Model'.padEnd(38) + 'R²'.padStart(7) + 'LOO R²'.padStart(9) + '  Gap%'.padStart(7) + '  #p');
console.log('  ' + '-'.repeat(62));
for (const m of models) {
  const gap = ((1 - m.loo / m.r2) * 100).toFixed(1);
  console.log('  ' + m.name.padEnd(38) + m.r2.toFixed(3).padStart(7) + m.loo.toFixed(3).padStart(9) + gap.padStart(7) + m.p.toString().padStart(4));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  KEY QUESTION: Does adding a₀ improve the inner-dynamics model?');
console.log('══════════════════════════════════════════════════════════════\n');

const M8_X = g => [1, g.logInnerVmax, g.innerDiscrepancy, g.innerCurvature, g.logA0inner];
const fitM8 = olsRegress(hi.map(M8_X), hi.map(y_fn));
const looM8 = looR2(M8_X, y_fn, hi);
console.log('  M5 (inner only):         LOO R² = ' + looM5.toFixed(3));
console.log('  M5 + a₀(inner):          LOO R² = ' + looM8.toFixed(3));
const improvement = looM8 - looM5;
console.log('  Improvement from adding a₀: ' + (improvement > 0 ? '+' : '') + improvement.toFixed(3));
if (improvement < 0.01) {
  console.log('  → a₀ adds NOTHING to the inner dynamics model');
  console.log('  → a₀ is REDUNDANT when inner dynamics are included');
} else {
  console.log('  → a₀ adds marginal information beyond inner dynamics');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  REVERSE: Does innerVmax add to a₀?');
console.log('══════════════════════════════════════════════════════════════\n');

const M9_X = g => [1, g.logA0inner, g.logInnerVmax];
const fitM9 = olsRegress(hi.map(M9_X), hi.map(y_fn));
const looM9 = looR2(M9_X, y_fn, hi);
console.log('  a₀(inner) only:          LOO R² = ' + looM7.toFixed(3));
console.log('  a₀(inner) + innerVmax:   LOO R² = ' + looM9.toFixed(3));
const improvement2 = looM9 - looM7;
console.log('  Improvement from adding innerVmax: ' + (improvement2 > 0 ? '+' : '') + improvement2.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PERMUTATION TEST on best model');
console.log('══════════════════════════════════════════════════════════════\n');

const nPerm = 2000;
let nBetter = 0;
const yVals = hi.map(y_fn);
for (let i = 0; i < nPerm; i++) {
  const shuffled = [...yVals].sort(() => Math.random() - 0.5);
  const fitP = olsRegress(hi.map(M3_X), shuffled);
  if (fitP && fitP.r2 >= fitM3.r2) nBetter++;
}
console.log('  M3 (Structural + innerVmax): perm p = ' + (nBetter / nPerm).toFixed(4));

let nBetter2 = 0;
for (let i = 0; i < nPerm; i++) {
  const shuffled = [...yVals].sort(() => Math.random() - 0.5);
  const fitP = olsRegress(hi.map(M6_X), shuffled);
  if (fitP && fitP.r2 >= fitM6.r2) nBetter2++;
}
console.log('  M6 (innerVmax only): perm p = ' + (nBetter2 / nPerm).toFixed(4));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const bestDirect = models.reduce((a, b) => a.loo > b.loo ? a : b);
console.log('  Best model: ' + bestDirect.name + ' (LOO R² = ' + bestDirect.loo.toFixed(3) + ')');

const a0needed = improvement > 0.01;
if (!a0needed) {
  console.log('\n  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  a₀ is NOT NEEDED: inner dynamics fully subsume it          ║');
  console.log('  ║  The direct inner→outer law supersedes the a₀ framework    ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else {
  console.log('\n  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  a₀ retains MARGINAL value beyond inner dynamics            ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
}

const output = {
  phase: '96',
  title: 'Direct Dynamical Law Without a₀',
  hiN: hi.length,
  models: models.map(m => ({ name: m.name, r2: +m.r2.toFixed(3), looR2: +m.loo.toFixed(3), nPredictors: m.p })),
  a0_adds_to_inner: improvement > 0.01,
  a0_improvement: +improvement.toFixed(3),
  innerVmax_adds_to_a0: +improvement2.toFixed(3),
  bestModel: bestDirect.name,
  bestLooR2: +bestDirect.loo.toFixed(3),
  m3_coefficients: {
    intercept: +fitM3.beta[0].toFixed(3),
    logSigma0: +fitM3.beta[1].toFixed(3),
    logVflat: +fitM3.beta[2].toFixed(3),
    fgas: +fitM3.beta[3].toFixed(3),
    logMHI: +fitM3.beta[4].toFixed(3),
    logInnerVmax: +fitM3.beta[5].toFixed(3),
  },
  m6_equation: 'outerSlope = ' + fitM6.beta[0].toFixed(3) + ' + ' + fitM6.beta[1].toFixed(3) + ' * logInnerVmax',
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase96-direct-law.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase96-direct-law.json');
