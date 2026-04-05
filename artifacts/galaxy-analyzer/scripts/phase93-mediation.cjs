#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 93: MEDIATION TEST');
console.log('  Is Σ₀ the parent variable driving both a₀ and flatness?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0(pts) {
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
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, rms: Math.sqrt(ss_res / n), residuals, yhat };
}

function looCV(X, y) {
  const n = y.length; let ss = 0;
  for (let i = 0; i < n; i++) {
    const Xt = X.filter((_, j) => j !== i);
    const yt = y.filter((_, j) => j !== i);
    const fit = olsRegress(Xt, yt);
    if (!fit) return NaN;
    ss += (y[i] - X[i].reduce((s, v, j) => s + v * fit.beta[j], 0)) ** 2;
  }
  return Math.sqrt(ss / n);
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

const p25 = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase25-group-membership.json'), 'utf-8'));
const envLookup = {};
for (const g of p25.galaxyAssignments) envLookup[g.name] = g;

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
    pts.push({ r: pt.rad, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), vobs: pt.vobs, vgas: pt.vgas, vdisk: pt.vdisk });
  }
  if (pts.length < 5) continue;

  const fit = fitA0(pts);
  if (!isFinite(fit.logA0) || fit.logA0 < 0.5 || fit.logA0 > 5.5) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const outer = sorted.slice(half);
  const logR = outer.map(p => Math.log10(p.r));
  const logV = outer.map(p => Math.log10(p.vobs));
  const mrr = mean(logR), mvv = mean(logV);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mrr) * (logV[i] - mvv); sxx += (logR[i] - mrr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const rarResid = sorted.map(p => {
    const gb = Math.pow(10, p.logGbar);
    const pred = mcgaughRAR(gb, fit.a0);
    return p.logGobs - Math.log10(pred > 0 ? pred : 1e-20);
  });
  let currentSign = Math.sign(rarResid[0]);
  let runLen = 1;
  const residRun = [];
  for (let i = 1; i < rarResid.length; i++) {
    const s = Math.sign(rarResid[i]);
    if (s === currentSign) runLen++;
    else { residRun.push(runLen); runLen = 1; currentSign = s; }
  }
  residRun.push(runLen);
  const meanRunLen = residRun.reduce((s, v) => s + v, 0) / residRun.length;

  const ev = envLookup[name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;
  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const Mstar = t1.L36 * UPSILON_DISK;
  const diskDom = pts.length > 0 ? mean(pts.map(p => {
    const vd2 = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk);
    const vg2 = p.vgas * Math.abs(p.vgas);
    const total = Math.abs(vd2) + Math.abs(vg2);
    return total > 0 ? Math.abs(vd2) / total : 0.5;
  })) : 0.5;

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI: Math.log10(t1.MHI),
    logMeanRun: Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1),
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    Vflat: t1.Vflat, Q: t1.Q, T: t1.T,
    envCode, outerSlope, fgas, diskDom,
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    logFgas: Math.log10(fgas > 0 ? fgas : 1e-3),
    logMstar: Math.log10(Mstar),
  });
}

const hi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);
console.log('  High-regime: N=' + hi.length + '\n');

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: PARTIAL CORRELATION MATRIX');
console.log('  (controlling for each candidate parent)');
console.log('══════════════════════════════════════════════════════════════\n');

function partialCorr(x, y, z, data) {
  const xv = data.map(x), yv = data.map(y), zv = data.map(z);
  const fitXZ = olsRegress(zv.map(v => [1, v]), xv);
  const fitYZ = olsRegress(zv.map(v => [1, v]), yv);
  if (!fitXZ || !fitYZ) return NaN;
  return pearsonR(fitXZ.residuals, fitYZ.residuals);
}

function partialCorrMulti(xFn, yFn, zFns, data) {
  const xv = data.map(xFn), yv = data.map(yFn);
  const Z = data.map(g => [1, ...zFns.map(fn => fn(g))]);
  const fitXZ = olsRegress(Z, xv);
  const fitYZ = olsRegress(Z, yv);
  if (!fitXZ || !fitYZ) return NaN;
  return pearsonR(fitXZ.residuals, fitYZ.residuals);
}

const candidates = [
  { name: 'logSigma0', fn: g => g.logSigma0 },
  { name: 'logVflat', fn: g => g.logVflat },
  { name: 'fgas', fn: g => g.fgas },
  { name: 'logMHI', fn: g => g.logMHI },
  { name: 'logMstar', fn: g => g.logMstar },
  { name: 'diskDom', fn: g => g.diskDom },
  { name: 'T', fn: g => g.T },
];

console.log('  Raw correlations in high regime:');
console.log('  r(logA0, outerSlope) = ' + pearsonR(hi.map(g => g.logA0), hi.map(g => g.outerSlope)).toFixed(3));
for (const c of candidates) {
  const rA0 = pearsonR(hi.map(c.fn), hi.map(g => g.logA0));
  const rSlope = pearsonR(hi.map(c.fn), hi.map(g => g.outerSlope));
  console.log('  r(' + c.name + ', logA0) = ' + rA0.toFixed(3) + ',  r(' + c.name + ', slope) = ' + rSlope.toFixed(3));
}

console.log('\n  MEDIATION TEST: Does controlling for X eliminate r(A, B)?');
console.log('  A = logA0, B = outerSlope\n');

console.log('  ' + 'Control'.padEnd(15) + 'r(A,B|X)'.padStart(10) + 'Δ from raw'.padStart(12) + 'Mediation?'.padStart(12));
const rawAB = pearsonR(hi.map(g => g.logA0), hi.map(g => g.outerSlope));

const mediationResults = {};
for (const c of candidates) {
  const pr = partialCorr(g => g.logA0, g => g.outerSlope, c.fn, hi);
  const delta = pr - rawAB;
  const mediation = Math.abs(pr) < Math.abs(rawAB) * 0.5 ? 'STRONG' : Math.abs(pr) < Math.abs(rawAB) * 0.8 ? 'PARTIAL' : 'NONE';
  mediationResults[c.name] = { partial: +pr.toFixed(3), delta: +delta.toFixed(3), mediation };
  console.log('  ' + c.name.padEnd(15) + pr.toFixed(3).padStart(10) + delta.toFixed(3).padStart(12) + mediation.padStart(12));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: REVERSE MEDIATION');
console.log('  Does Σ₀ independently predict BOTH a₀ and slope?');
console.log('══════════════════════════════════════════════════════════════\n');

const rSigA0 = pearsonR(hi.map(g => g.logSigma0), hi.map(g => g.logA0));
const rSigSlope = pearsonR(hi.map(g => g.logSigma0), hi.map(g => g.outerSlope));
console.log('  Raw: r(Σ₀, logA0) = ' + rSigA0.toFixed(3) + ', r(Σ₀, slope) = ' + rSigSlope.toFixed(3));

const prSigA0givenSlope = partialCorr(g => g.logSigma0, g => g.logA0, g => g.outerSlope, hi);
const prSigSlopegivenA0 = partialCorr(g => g.logSigma0, g => g.outerSlope, g => g.logA0, hi);
console.log('  r(Σ₀, logA0 | slope) = ' + prSigA0givenSlope.toFixed(3) + '  (Σ₀→a₀ path survives controlling slope?)');
console.log('  r(Σ₀, slope | logA0) = ' + prSigSlopegivenA0.toFixed(3) + '  (Σ₀→slope path survives controlling a₀?)');

if (Math.abs(prSigA0givenSlope) > 0.1 && Math.abs(prSigSlopegivenA0) > 0.4) {
  console.log('  → Σ₀ has INDEPENDENT paths to BOTH a₀ and slope');
  console.log('    (consistent with Σ₀ as common cause)');
} else if (Math.abs(prSigSlopegivenA0) > 0.3) {
  console.log('  → Σ₀→slope path is robust; Σ₀→a₀ path weaker');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: FULL CANDIDATE PARENT SEARCH');
console.log('  Which variable, when controlled, most reduces r(X, slope)?');
console.log('  for each predictor X in the a₀ law');
console.log('══════════════════════════════════════════════════════════════\n');

const a0preds = [
  { name: 'logMHI', fn: g => g.logMHI },
  { name: 'logMeanRun', fn: g => g.logMeanRun },
  { name: 'logSigma0', fn: g => g.logSigma0 },
  { name: 'logVflat', fn: g => g.logVflat },
];

for (const pred of a0preds) {
  const rawR = pearsonR(hi.map(pred.fn), hi.map(g => g.outerSlope));
  console.log('  ' + pred.name + ' → slope: raw r = ' + rawR.toFixed(3));
  for (const ctrl of candidates) {
    if (ctrl.name === pred.name) continue;
    const pr = partialCorr(pred.fn, g => g.outerSlope, ctrl.fn, hi);
    const pct = rawR !== 0 ? ((rawR - pr) / rawR * 100).toFixed(0) : '---';
    console.log('    controlling ' + ctrl.name.padEnd(12) + ': r = ' + pr.toFixed(3) + '  (reduction: ' + pct + '%)');
  }
  console.log();
}

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 4: MULTIVARIATE PARTIAL CORRELATIONS');
console.log('  Control for ALL structural variables simultaneously');
console.log('══════════════════════════════════════════════════════════════\n');

const allControls = [g => g.logSigma0, g => g.logVflat, g => g.fgas, g => g.logMHI];

for (const pred of a0preds) {
  const controls = allControls.filter(fn => fn !== pred.fn);
  const pr = partialCorrMulti(pred.fn, g => g.outerSlope, controls, hi);
  console.log('  r(' + pred.name + ', slope | all others) = ' + (isNaN(pr) ? '---' : pr.toFixed(3)));
}

const prA0slope = partialCorrMulti(g => g.logA0, g => g.outerSlope, allControls, hi);
console.log('  r(logA0, slope | all structural) = ' + (isNaN(prA0slope) ? '---' : prA0slope.toFixed(3)));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: BARYON DOMINANCE AS THE ROOT');
console.log('  Composite: diskDom + logSigma0 + fgas');
console.log('══════════════════════════════════════════════════════════════\n');

const bdX = hi.map(g => [1, g.logSigma0, g.fgas, g.diskDom]);
const bdA0 = olsRegress(bdX, hi.map(g => g.logA0));
const bdSlope = olsRegress(bdX, hi.map(g => g.outerSlope));

console.log('  Baryon dominance composite → logA0:');
if (bdA0) console.log('    R² = ' + bdA0.r2.toFixed(4) + ', coefficients = ' + bdA0.beta.map(b => b.toFixed(4)).join(', '));

console.log('  Baryon dominance composite → outerSlope:');
if (bdSlope) console.log('    R² = ' + bdSlope.r2.toFixed(4) + ', coefficients = ' + bdSlope.beta.map(b => b.toFixed(4)).join(', '));

if (bdA0 && bdSlope) {
  const prResid = pearsonR(bdA0.residuals, bdSlope.residuals);
  console.log('\n  After removing baryon dominance from both:');
  console.log('    r(logA0 resid, slope resid | BD) = ' + prResid.toFixed(3));
  if (Math.abs(prResid) < 0.1) {
    console.log('    → Baryon dominance FULLY mediates the a₀-slope connection');
  } else {
    console.log('    → Residual connection remains: baryon dominance is partial mediator');
  }
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 6: DOES Σ₀ ALONE EXPLAIN FLATNESS?');
console.log('  Regression: slope = f(logSigma0) with LOO');
console.log('══════════════════════════════════════════════════════════════\n');

const sigX = hi.map(g => [1, g.logSigma0]);
const slopeY = hi.map(g => g.outerSlope);
const sigFit = olsRegress(sigX, slopeY);
const sigLOO = looCV(sigX, slopeY);
const m0sd = sd(slopeY);

console.log('  Σ₀ alone → outerSlope:');
console.log('    R² = ' + sigFit.r2.toFixed(4));
console.log('    LOO RMS = ' + sigLOO.toFixed(4) + ', M0 sd = ' + m0sd.toFixed(4));
console.log('    LOO gap = ' + ((m0sd - sigLOO) / m0sd * 100).toFixed(1) + '%');
console.log('    Coefficients: intercept=' + sigFit.beta[0].toFixed(4) + ', slope=' + sigFit.beta[1].toFixed(4));

const fgasX = hi.map(g => [1, g.fgas]);
const fgasFit = olsRegress(fgasX, slopeY);
const fgasLOO = looCV(fgasX, slopeY);
console.log('\n  fgas alone → outerSlope:');
console.log('    R² = ' + fgasFit.r2.toFixed(4));
console.log('    LOO gap = ' + ((m0sd - fgasLOO) / m0sd * 100).toFixed(1) + '%');

const sigFgasX = hi.map(g => [1, g.logSigma0, g.fgas]);
const sigFgasFit = olsRegress(sigFgasX, slopeY);
const sigFgasLOO = looCV(sigFgasX, slopeY);
console.log('\n  Σ₀ + fgas → outerSlope:');
console.log('    R² = ' + sigFgasFit.r2.toFixed(4));
console.log('    LOO gap = ' + ((m0sd - sigFgasLOO) / m0sd * 100).toFixed(1) + '%');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 7: CAUSAL DIRECTION — Σ₀ vs Vflat');
console.log('  Which is the more fundamental parent?');
console.log('══════════════════════════════════════════════════════════════\n');

const prSigSlopeCtrlV = partialCorr(g => g.logSigma0, g => g.outerSlope, g => g.logVflat, hi);
const prVSlopeCtrlSig = partialCorr(g => g.logVflat, g => g.outerSlope, g => g.logSigma0, hi);
const prSigA0ctrlV = partialCorr(g => g.logSigma0, g => g.logA0, g => g.logVflat, hi);
const prVA0ctrlSig = partialCorr(g => g.logVflat, g => g.logA0, g => g.logSigma0, hi);

console.log('  Predicting outerSlope:');
console.log('    r(Σ₀, slope | Vflat) = ' + prSigSlopeCtrlV.toFixed(3) + '  ← Σ₀ survives controlling Vflat?');
console.log('    r(Vflat, slope | Σ₀) = ' + prVSlopeCtrlSig.toFixed(3) + '  ← Vflat survives controlling Σ₀?');
console.log();
console.log('  Predicting logA0:');
console.log('    r(Σ₀, A0 | Vflat) = ' + prSigA0ctrlV.toFixed(3) + '  ← Σ₀ survives?');
console.log('    r(Vflat, A0 | Σ₀) = ' + prVA0ctrlSig.toFixed(3) + '  ← Vflat survives?');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 8: PATH ANALYSIS SUMMARY');
console.log('══════════════════════════════════════════════════════════════\n');

const paths = {};

paths['Sigma0_to_slope'] = { raw: +rSigSlope.toFixed(3), ctrlA0: +prSigSlopegivenA0.toFixed(3), ctrlVflat: +prSigSlopeCtrlV.toFixed(3) };
paths['Sigma0_to_A0'] = { raw: +rSigA0.toFixed(3), ctrlSlope: +prSigA0givenSlope.toFixed(3), ctrlVflat: +prSigA0ctrlV.toFixed(3) };
paths['Vflat_to_slope'] = { raw: +pearsonR(hi.map(g => g.logVflat), hi.map(g => g.outerSlope)).toFixed(3), ctrlSigma0: +prVSlopeCtrlSig.toFixed(3) };
paths['Vflat_to_A0'] = { raw: +pearsonR(hi.map(g => g.logVflat), hi.map(g => g.logA0)).toFixed(3), ctrlSigma0: +prVA0ctrlSig.toFixed(3) };
paths['A0_to_slope'] = { raw: +rawAB.toFixed(3), ctrlSigma0: +(partialCorr(g => g.logA0, g => g.outerSlope, g => g.logSigma0, hi)).toFixed(3) };

console.log('  Path strengths (raw → partial):');
for (const [key, val] of Object.entries(paths)) {
  console.log('    ' + key + ': ' + JSON.stringify(val));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  SYNTHESIS');
console.log('══════════════════════════════════════════════════════════════\n');

const sigSurvivesForSlope = Math.abs(prSigSlopeCtrlV) > 0.2;
const sigSurvivesForA0 = Math.abs(prSigA0ctrlV) > 0.1;
const vflatSurvivesForSlope = Math.abs(prVSlopeCtrlSig) > 0.2;

let verdict;
if (sigSurvivesForSlope && sigSurvivesForA0 && !vflatSurvivesForSlope) {
  verdict = 'SIGMA0-IS-PARENT';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  Σ₀ IS THE PRIMARY PARENT VARIABLE                         ║');
  console.log('  ║  It independently drives both a₀ and flatness              ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else if (sigSurvivesForSlope && vflatSurvivesForSlope) {
  verdict = 'DUAL-PARENTS';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  DUAL PARENT STRUCTURE: Σ₀ + Vflat both contribute         ║');
  console.log('  ║  independently to slope (and partially to a₀)              ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else {
  verdict = 'COMPLEX-MEDIATION';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  COMPLEX MEDIATION: no single parent dominates             ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
}

const output = {
  phase: '93',
  title: 'Mediation Test: Is Σ₀ the Parent Variable?',
  hiN: hi.length,
  mediationResults,
  paths,
  baryonDominance: {
    a0R2: bdA0 ? +bdA0.r2.toFixed(4) : null,
    slopeR2: bdSlope ? +bdSlope.r2.toFixed(4) : null,
    residualR: bdA0 && bdSlope ? +pearsonR(bdA0.residuals, bdSlope.residuals).toFixed(3) : null,
  },
  sigma0alone: { r2: +sigFit.r2.toFixed(4), looGap: +((m0sd - sigLOO) / m0sd * 100).toFixed(1) },
  verdict,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase93-mediation.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase93-mediation.json');
