#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 87: REGIME TRANSITION MAPPING');
console.log('  Where exactly does the MHI-a0 correlation flip?');
console.log('  Is the transition sharp or gradual?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

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
  const logA0 = (lo + hi) / 2, a0 = Math.pow(10, logA0);
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.logGbar);
    const pred = mcgaughRAR(gb, a0);
    ss += (p.logGobs - Math.log10(pred > 0 ? pred : 1e-20)) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length) };
}

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }

function pearsonR(x, y) {
  const n = x.length;
  if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function olsRegress(X, y) {
  const n = y.length, p = X[0].length;
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
  return { beta, r2: 1 - ss_res / ss_tot, rms: Math.sqrt(ss_res / n), residuals, yhat };
}

function looCV(X, y) {
  const n = y.length;
  let ss = 0;
  for (let i = 0; i < n; i++) {
    const Xtrain = X.filter((_, j) => j !== i);
    const ytrain = y.filter((_, j) => j !== i);
    const fit = olsRegress(Xtrain, ytrain);
    if (!fit) return NaN;
    const pred = X[i].reduce((s, v, j) => s + v * fit.beta[j], 0);
    ss += (y[i] - pred) ** 2;
  }
  return Math.sqrt(ss / n);
}

const csvPath = path.join(__dirname, '..', 'public', 'replication', 'N45_final_dataset.csv');
const csvLines = fs.readFileSync(csvPath, 'utf-8').trim().split('\n');
const headers = csvLines[0].split(',');
const trainCSV = csvLines.slice(1).map(l => {
  const vals = l.split(',');
  const obj = {};
  headers.forEach((h, i) => obj[h] = vals[i]);
  return obj;
});
const trainNames = new Set(trainCSV.map(d => d.galaxy_name));

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const D = parseFloat(line.substring(15, 21).trim());
  const fD = parseInt(line.substring(28, 29).trim());
  const inc = parseFloat(line.substring(30, 34).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { T, D, fD, inc, L36, Rdisk, MHI, Vflat, Q };
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
    pts.push({ r: pt.rad, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar) });
  }
  if (pts.length < 3) continue;
  const fit = fitA0(pts);
  if (!isFinite(fit.logA0) || fit.logA0 < 0.5 || fit.logA0 > 5.5) continue;
  const sorted = [...pts].sort((a, b) => a.r - b.r);
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
  const logMeanRun = Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1);
  const ev = envLookup[name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI: Math.log10(t1.MHI),
    logMeanRun,
    Vflat: t1.Vflat, Q: t1.Q, D: t1.D, fD: t1.fD, T: t1.T,
    L36: t1.L36, inc: t1.inc, envCode,
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    logL36: Math.log10(t1.L36),
    inTraining: trainNames.has(name),
  });
}

console.log('  Total SPARC galaxies with valid a0: ' + allGalaxies.length + '\n');

console.log('======================================================================');
console.log('  TEST 1: SLIDING WINDOW — r(MHI,a0) as function of Vflat');
console.log('  Window width = 30 galaxies, step = 5 km/s');
console.log('======================================================================\n');

const sortedByVflat = [...allGalaxies].sort((a, b) => a.Vflat - b.Vflat);
const windowSize = 30;
const slidingResults = [];

console.log('  ' + 'Vflat_center'.padEnd(14) + 'Vflat_range'.padEnd(18) + 'N'.padStart(4) +
  'r(MHI)'.padStart(9) + 'r(Run)'.padStart(9) + 'mean_logA0'.padStart(12));

for (let startIdx = 0; startIdx <= sortedByVflat.length - windowSize; startIdx += 3) {
  const window = sortedByVflat.slice(startIdx, startIdx + windowSize);
  const vMin = window[0].Vflat;
  const vMax = window[window.length - 1].Vflat;
  const vCenter = (vMin + vMax) / 2;
  const rMHI = pearsonR(window.map(g => g.logMHI), window.map(g => g.logA0));
  const rRun = pearsonR(window.map(g => g.logMeanRun), window.map(g => g.logA0));
  const mLogA0 = mean(window.map(g => g.logA0));

  slidingResults.push({
    vCenter: +vCenter.toFixed(1), vMin: +vMin.toFixed(1), vMax: +vMax.toFixed(1),
    N: window.length, rMHI: +rMHI.toFixed(3), rRun: +rRun.toFixed(3),
    meanLogA0: +mLogA0.toFixed(3),
  });

  console.log('  ' + vCenter.toFixed(1).padEnd(14) + (vMin.toFixed(0) + '-' + vMax.toFixed(0)).padEnd(18) +
    String(window.length).padStart(4) + rMHI.toFixed(3).padStart(9) + rRun.toFixed(3).padStart(9) +
    mLogA0.toFixed(3).padStart(12));
}

const flipIdx = slidingResults.findIndex(s => s.rMHI < 0);
const lastPos = flipIdx > 0 ? slidingResults[flipIdx - 1] : null;
const firstNeg = flipIdx >= 0 ? slidingResults[flipIdx] : null;

console.log('\n  Sign flip detected:');
if (lastPos && firstNeg) {
  console.log('    Last positive window: Vflat_center=' + lastPos.vCenter + ', r=' + lastPos.rMHI);
  console.log('    First negative window: Vflat_center=' + firstNeg.vCenter + ', r=' + firstNeg.rMHI);
  console.log('    Transition zone: ~' + lastPos.vCenter + ' to ~' + firstNeg.vCenter + ' km/s');
}

console.log('\n======================================================================');
console.log('  TEST 2: CUMULATIVE FROM BELOW — r(MHI,a0) for Vflat <= threshold');
console.log('======================================================================\n');

const cumBelow = [];
console.log('  ' + 'Vflat_max'.padEnd(12) + 'N'.padStart(4) + 'r(MHI)'.padStart(9) + 'r(Run)'.padStart(9));

for (let vmax = 40; vmax <= 250; vmax += 5) {
  const sub = allGalaxies.filter(g => g.Vflat <= vmax);
  if (sub.length < 10) continue;
  const r = pearsonR(sub.map(g => g.logMHI), sub.map(g => g.logA0));
  const rR = pearsonR(sub.map(g => g.logMeanRun), sub.map(g => g.logA0));
  cumBelow.push({ vMax: vmax, N: sub.length, rMHI: +r.toFixed(3), rRun: +rR.toFixed(3) });
  console.log('  ' + String(vmax).padEnd(12) + String(sub.length).padStart(4) +
    r.toFixed(3).padStart(9) + rR.toFixed(3).padStart(9));
}

console.log('\n======================================================================');
console.log('  TEST 3: CUMULATIVE FROM ABOVE — r(MHI,a0) for Vflat >= threshold');
console.log('======================================================================\n');

const cumAbove = [];
console.log('  ' + 'Vflat_min'.padEnd(12) + 'N'.padStart(4) + 'r(MHI)'.padStart(9) + 'r(Run)'.padStart(9));

for (let vmin = 20; vmin <= 180; vmin += 5) {
  const sub = allGalaxies.filter(g => g.Vflat >= vmin);
  if (sub.length < 10) continue;
  const r = pearsonR(sub.map(g => g.logMHI), sub.map(g => g.logA0));
  const rR = pearsonR(sub.map(g => g.logMeanRun), sub.map(g => g.logA0));
  cumAbove.push({ vMin: vmin, N: sub.length, rMHI: +r.toFixed(3), rRun: +rR.toFixed(3) });
  console.log('  ' + String(vmin).padEnd(12) + String(sub.length).padStart(4) +
    r.toFixed(3).padStart(9) + rR.toFixed(3).padStart(9));
}

console.log('\n======================================================================');
console.log('  TEST 4: BROKEN REGRESSION — find optimal Vflat breakpoint');
console.log('  Fit separate MHI slopes above/below each Vflat threshold');
console.log('======================================================================\n');

const brokenResults = [];
let bestBIC = Infinity, bestBreak = NaN;

console.log('  ' + 'Vflat_break'.padEnd(14) + 'N_lo'.padStart(5) + 'N_hi'.padStart(5) +
  'slope_lo'.padStart(10) + 'slope_hi'.padStart(10) + 'R2_total'.padStart(10) +
  'BIC'.padStart(10) + 'LOO_RMS'.padStart(10));

for (let vbreak = 40; vbreak <= 160; vbreak += 5) {
  const lo = allGalaxies.filter(g => g.Vflat < vbreak);
  const hi = allGalaxies.filter(g => g.Vflat >= vbreak);
  if (lo.length < 10 || hi.length < 10) continue;

  const indicator = allGalaxies.map(g => g.Vflat >= vbreak ? 1 : 0);
  const X = allGalaxies.map((g, i) => [1, g.logMHI, indicator[i], g.logMHI * indicator[i]]);
  const y = allGalaxies.map(g => g.logA0);

  const fit = olsRegress(X, y);
  if (!fit) continue;

  const slopeLo = fit.beta[1];
  const slopeHi = fit.beta[1] + fit.beta[3];
  const n = allGalaxies.length;
  const k = 4;
  const sse = fit.residuals.reduce((s, r) => s + r * r, 0);
  const bic = n * Math.log(sse / n) + k * Math.log(n);

  const looRMS = looCV(X, y);

  brokenResults.push({
    vBreak: vbreak, nLo: lo.length, nHi: hi.length,
    slopeLo: +slopeLo.toFixed(4), slopeHi: +slopeHi.toFixed(4),
    r2: +fit.r2.toFixed(4), bic: +bic.toFixed(2), looRMS: +looRMS.toFixed(4),
  });

  if (bic < bestBIC) { bestBIC = bic; bestBreak = vbreak; }

  console.log('  ' + String(vbreak).padEnd(14) + String(lo.length).padStart(5) +
    String(hi.length).padStart(5) + slopeLo.toFixed(4).padStart(10) +
    slopeHi.toFixed(4).padStart(10) + fit.r2.toFixed(4).padStart(10) +
    bic.toFixed(2).padStart(10) + looRMS.toFixed(4).padStart(10));
}

console.log('\n  Best BIC breakpoint: Vflat = ' + bestBreak + ' km/s (BIC = ' + bestBIC.toFixed(2) + ')');

const noBreakX = allGalaxies.map(g => [1, g.logMHI]);
const noBreakY = allGalaxies.map(g => g.logA0);
const noBreakFit = olsRegress(noBreakX, noBreakY);
const noBreakSSE = noBreakFit.residuals.reduce((s, r) => s + r * r, 0);
const noBreakBIC = allGalaxies.length * Math.log(noBreakSSE / allGalaxies.length) + 2 * Math.log(allGalaxies.length);
const noBreakLOO = looCV(noBreakX, noBreakY);
console.log('  No-break (single slope): BIC = ' + noBreakBIC.toFixed(2) + ', R2 = ' + noBreakFit.r2.toFixed(4) +
  ', LOO_RMS = ' + noBreakLOO.toFixed(4));
console.log('  ΔBIC = ' + (bestBIC - noBreakBIC).toFixed(2) + ' (negative = broken is better)\n');

console.log('======================================================================');
console.log('  TEST 5: INTERACTION MODEL — logMHI × logVflat');
console.log('  Does Vflat continuously modulate the MHI-a0 relationship?');
console.log('======================================================================\n');

const X_int = allGalaxies.map(g => [1, g.logMHI, g.logVflat, g.logMHI * g.logVflat]);
const y_int = allGalaxies.map(g => g.logA0);
const fit_int = olsRegress(X_int, y_int);
const sse_int = fit_int.residuals.reduce((s, r) => s + r * r, 0);
const bic_int = allGalaxies.length * Math.log(sse_int / allGalaxies.length) + 4 * Math.log(allGalaxies.length);
const loo_int = looCV(X_int, y_int);

console.log('  Interaction model: logA0 = b0 + b1*logMHI + b2*logVflat + b3*logMHI*logVflat');
console.log('  Coefficients: ' + fit_int.beta.map(b => b.toFixed(4)).join(', '));
console.log('  R2 = ' + fit_int.r2.toFixed(4) + ', BIC = ' + bic_int.toFixed(2) + ', LOO_RMS = ' + loo_int.toFixed(4));
console.log('  Interaction term (b3): ' + fit_int.beta[3].toFixed(4));
if (Math.abs(fit_int.beta[3]) > 0.05) {
  console.log('  → Interaction is substantial: Vflat continuously modulates MHI-a0 slope');
} else {
  console.log('  → Interaction is weak');
}

const zeroInteraction = -fit_int.beta[1] / fit_int.beta[3];
console.log('  MHI slope = 0 when logVflat = ' + zeroInteraction.toFixed(2) +
  ' → Vflat = ' + Math.pow(10, zeroInteraction).toFixed(0) + ' km/s');
console.log();

const X_add = allGalaxies.map(g => [1, g.logMHI, g.logVflat]);
const fit_add = olsRegress(X_add, y_int);
const sse_add = fit_add.residuals.reduce((s, r) => s + r * r, 0);
const bic_add = allGalaxies.length * Math.log(sse_add / allGalaxies.length) + 3 * Math.log(allGalaxies.length);
const loo_add = looCV(X_add, y_int);
console.log('  Additive model (no interaction): R2 = ' + fit_add.r2.toFixed(4) +
  ', BIC = ' + bic_add.toFixed(2) + ', LOO_RMS = ' + loo_add.toFixed(4));
console.log('  ΔBIC(int vs add) = ' + (bic_int - bic_add).toFixed(2) + ' (negative = interaction better)');

console.log('\n======================================================================');
console.log('  TEST 6: LOGISTIC TRANSITION MODEL');
console.log('  Fit r(MHI,a0) = r_hi + (r_lo - r_hi) / (1 + exp((Vflat - V0)/w))');
console.log('  to the sliding window results');
console.log('======================================================================\n');

if (slidingResults.length > 5) {
  let bestFitSSE = Infinity, bestParams = null;
  for (let V0 = 40; V0 <= 140; V0 += 2) {
    for (let w = 5; w <= 50; w += 2) {
      for (let rLo = 0.1; rLo <= 0.6; rLo += 0.05) {
        for (let rHi = -0.5; rHi <= -0.05; rHi += 0.05) {
          let sse = 0;
          for (const s of slidingResults) {
            const pred = rHi + (rLo - rHi) / (1 + Math.exp((s.vCenter - V0) / w));
            sse += (s.rMHI - pred) ** 2;
          }
          if (sse < bestFitSSE) {
            bestFitSSE = sse;
            bestParams = { V0, w, rLo, rHi };
          }
        }
      }
    }
  }

  if (bestParams) {
    console.log('  Best logistic fit:');
    console.log('    V0 (midpoint) = ' + bestParams.V0 + ' km/s');
    console.log('    w (width) = ' + bestParams.w + ' km/s');
    console.log('    r_lo (dwarf regime) = ' + bestParams.rLo.toFixed(2));
    console.log('    r_hi (massive regime) = ' + bestParams.rHi.toFixed(2));
    console.log('    SSE = ' + bestFitSSE.toFixed(4));
    console.log('    Transition range (10%-90%): ' +
      (bestParams.V0 - 2.2 * bestParams.w).toFixed(0) + ' to ' +
      (bestParams.V0 + 2.2 * bestParams.w).toFixed(0) + ' km/s');

    const sharpness = bestParams.w < 15 ? 'SHARP' : bestParams.w < 30 ? 'MODERATE' : 'GRADUAL';
    console.log('    Transition character: ' + sharpness + ' (w=' + bestParams.w + ' km/s)');
  }
}

console.log('\n======================================================================');
console.log('  TEST 7: ALTERNATIVE TRANSITION VARIABLES');
console.log('  Is Vflat the best separator, or is logMHI / logL36 / T better?');
console.log('======================================================================\n');

const altVars = [
  { name: 'logVflat', fn: g => g.logVflat },
  { name: 'logMHI', fn: g => g.logMHI },
  { name: 'logL36', fn: g => g.logL36 },
  { name: 'T (morphology)', fn: g => g.T },
];

for (const av of altVars) {
  const sorted = [...allGalaxies].sort((a, b) => av.fn(a) - av.fn(b));
  let bestBICalt = Infinity, bestCutAlt = NaN;

  for (let pct = 15; pct <= 85; pct += 5) {
    const cutIdx = Math.floor(sorted.length * pct / 100);
    const cutVal = av.fn(sorted[cutIdx]);
    const indicator = allGalaxies.map(g => av.fn(g) >= cutVal ? 1 : 0);
    const nHi = indicator.filter(v => v).length;
    const nLo = indicator.filter(v => !v).length;
    if (nLo < 10 || nHi < 10) continue;

    const X = allGalaxies.map((g, i) => [1, g.logMHI, indicator[i], g.logMHI * indicator[i]]);
    const y = allGalaxies.map(g => g.logA0);
    const fit = olsRegress(X, y);
    if (!fit) continue;

    const sse = fit.residuals.reduce((s, r) => s + r * r, 0);
    const bic = allGalaxies.length * Math.log(sse / allGalaxies.length) + 4 * Math.log(allGalaxies.length);
    if (bic < bestBICalt) { bestBICalt = bic; bestCutAlt = cutVal; }
  }

  console.log('  ' + av.name + ': best BIC = ' + bestBICalt.toFixed(2) +
    ' at cut = ' + (typeof bestCutAlt === 'number' ? bestCutAlt.toFixed(2) : bestCutAlt) +
    ', ΔBIC vs no-break = ' + (bestBICalt - noBreakBIC).toFixed(2));
}

console.log('\n======================================================================');
console.log('  FINAL SUMMARY');
console.log('======================================================================\n');

const bestBroken = brokenResults.reduce((best, r) => r.bic < best.bic ? r : best, { bic: Infinity });
console.log('  OPTIMAL BREAKPOINT: Vflat = ' + bestBroken.vBreak + ' km/s');
console.log('    Below: slope(MHI→a0) = ' + bestBroken.slopeLo + ' (positive = more gas → higher a0)');
console.log('    Above: slope(MHI→a0) = ' + bestBroken.slopeHi + ' (negative = more gas → lower a0)');
console.log('    Broken model R2 = ' + bestBroken.r2 + ', LOO_RMS = ' + bestBroken.looRMS);
console.log('    ΔBIC vs single slope: ' + (bestBroken.bic - noBreakBIC).toFixed(2));
console.log();

if (bestBroken.slopeLo > 0 && bestBroken.slopeHi < 0) {
  console.log('  ╔══════════════════════════════════════════════════════════════════════╗');
  console.log('  ║  CONFIRMED: MHI-a0 relationship REVERSES across Vflat regimes      ║');
  console.log('  ║  This is a genuine regime transition, not a selection artifact      ║');
  console.log('  ╚══════════════════════════════════════════════════════════════════════╝');
} else {
  console.log('  NOTE: Slopes do not cleanly reverse — transition may be more complex');
}

const output = {
  phase: '87',
  title: 'Regime Transition Mapping — Where does the MHI-a0 correlation flip?',
  nGalaxies: allGalaxies.length,
  slidingWindow: slidingResults,
  cumulativeBelow: cumBelow,
  cumulativeAbove: cumAbove,
  brokenRegression: brokenResults,
  bestBreakpoint: bestBroken,
  noBreakBIC: +noBreakBIC.toFixed(2),
  deltaBIC: +(bestBroken.bic - noBreakBIC).toFixed(2),
  interactionModel: {
    beta: fit_int.beta.map(b => +b.toFixed(4)),
    r2: +fit_int.r2.toFixed(4),
    bic: +bic_int.toFixed(2),
    looRMS: +loo_int.toFixed(4),
    zeroSlopeVflat: +Math.pow(10, zeroInteraction).toFixed(0),
  },
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase87-regime-transition.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase87-regime-transition.json');
