#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 104: LOCKED EXTERNAL REPLICATION');
console.log('  Does fgas predict outer support requirement outside the training set?');
console.log('  Protocol: Lock model on SPARC, apply WITHOUT refit');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;
const N_SPLITS = 200;

let rngState = 42;
function seededRandom() {
  rngState = (rngState * 1664525 + 1013904223) & 0x7fffffff;
  return rngState / 0x7fffffff;
}

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function percentile(a, p) { const s = [...a].sort((x, y) => x - y); const i = p * (s.length - 1); const lo = Math.floor(i); return lo === s.length - 1 ? s[lo] : s[lo] + (i - lo) * (s[lo + 1] - s[lo]); }
function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}
function spearmanR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  function ranks(arr) {
    const sorted = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
    const r = new Array(n);
    for (let i = 0; i < n; i++) r[sorted[i].i] = i + 1;
    return r;
  }
  return pearsonR(ranks(x), ranks(y));
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
  const ss_res = y.reduce((s, v, i) => s + (v - yhat[i]) ** 2, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, residuals: y.map((v, i) => v - yhat[i]), yhat };
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const D = parseFloat(line.substring(15, 21).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  const RHI = parseFloat(line.substring(78, 86).trim());
  table1[name] = { T, D, L36, Rdisk, MHI, Vflat, Q, RHI };
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
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    const vBar = Math.sqrt(Math.max(vBarSq, 0.01));
    pts.push({ r: pt.rad, vobs: pt.vobs, vBar });
  }
  if (pts.length < 8) return null;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);

  const outerMassDisc = outerPts.map(p => (p.vobs * p.vobs) / (p.vBar * p.vBar));
  const logOMD = Math.log10(mean(outerMassDisc.filter(v => isFinite(v) && v > 0)));
  if (!isFinite(logOMD)) return null;

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const baryonCompact = t1.L36 * UPSILON_DISK / (t1.Rdisk > 0 ? t1.Rdisk : 1);

  if (!isFinite(fgas) || !isFinite(logOMD) || !isFinite(Sigma0) || Sigma0 <= 0) return null;

  return {
    name, nPts: sorted.length, logOMD,
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    Vflat: t1.Vflat,
    logMHI: Math.log10(t1.MHI > 0 ? t1.MHI : 0.001),
    logL36: Math.log10(t1.L36 > 0 ? t1.L36 : 0.001),
    logSigma0: Math.log10(Sigma0),
    fgas,
    logBaryonCompact: Math.log10(baryonCompact > 0 ? baryonCompact : 1e-3),
  };
}

const allSPARC = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0) continue;
  const g = buildGalaxy(name, rcPoints, t1);
  if (g) allSPARC.push(g);
}

const highVflat = allSPARC.filter(g => g.Vflat >= VFLAT_BREAK);
const lowVflat = allSPARC.filter(g => g.Vflat < VFLAT_BREAK && g.Vflat > 0);

console.log('  SPARC high-Vflat (training pool): N=' + highVflat.length);
console.log('  SPARC low-Vflat (domain transfer): N=' + lowVflat.length);

console.log('\n── STEP 0: LOCK THE MODEL ON FULL SPARC HIGH-VFLAT ──\n');

const fgasModel = olsRegress(highVflat.map(g => [1, g.fgas]), highVflat.map(g => g.logOMD));
const lockedSlope = fgasModel.beta[1];
const lockedIntercept = fgasModel.beta[0];
console.log('  LOCKED MODEL: logOMD = ' + lockedIntercept.toFixed(4) + ' + ' + lockedSlope.toFixed(4) + ' * fgas');
console.log('  Training R² = ' + fgasModel.r2.toFixed(4));
console.log('  Training r = ' + pearsonR(highVflat.map(g => g.fgas), highVflat.map(g => g.logOMD)).toFixed(4));

const sigmaModel = olsRegress(highVflat.map(g => [1, g.logSigma0]), highVflat.map(g => g.logOMD));
const compactModel = olsRegress(highVflat.map(g => [1, g.logBaryonCompact]), highVflat.map(g => g.logOMD));
const lumModel = olsRegress(highVflat.map(g => [1, g.logL36]), highVflat.map(g => g.logOMD));

console.log('\n  COMPARATOR LOCKED MODELS:');
console.log('  logSigma0:       slope=' + sigmaModel.beta[1].toFixed(4) + ', R²=' + sigmaModel.r2.toFixed(4));
console.log('  logBaryonCompact: slope=' + compactModel.beta[1].toFixed(4) + ', R²=' + compactModel.r2.toFixed(4));
console.log('  logL36:          slope=' + lumModel.beta[1].toFixed(4) + ', R²=' + lumModel.r2.toFixed(4));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 1: SPARC 70/30 RANDOM SPLIT (' + N_SPLITS + ' iterations)');
console.log('  Train on 70%, test LOCKED on 30% (no refit)');
console.log('══════════════════════════════════════════════════════════════\n');

function applyLocked(trainData, testData, varFn) {
  const fit = olsRegress(trainData.map(g => [1, varFn(g)]), trainData.map(g => g.logOMD));
  if (!fit) return { r: NaN, exR2: NaN, rmse: NaN, calSlope: NaN, calOffset: NaN };

  const predictions = testData.map(g => fit.beta[0] + fit.beta[1] * varFn(g));
  const actual = testData.map(g => g.logOMD);

  const r = pearsonR(predictions, actual);
  const rho = spearmanR(predictions, actual);
  const ss_res = actual.reduce((s, v, i) => s + (v - predictions[i]) ** 2, 0);
  const ss_tot = actual.reduce((s, v) => s + (v - mean(actual)) ** 2, 0);
  const exR2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
  const rmse = Math.sqrt(ss_res / actual.length);

  const calFit = olsRegress(predictions.map(p => [1, p]), actual);
  const calSlope = calFit ? calFit.beta[1] : NaN;
  const calOffset = calFit ? calFit.beta[0] : NaN;

  return { r, rho, exR2, rmse, calSlope, calOffset };
}

const splitResults = { fgas: [], sigma0: [], compact: [], lum: [] };
const varFns = {
  fgas: g => g.fgas,
  sigma0: g => g.logSigma0,
  compact: g => g.logBaryonCompact,
  lum: g => g.logL36
};

for (let iter = 0; iter < N_SPLITS; iter++) {
  const shuffled = [...highVflat];
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(seededRandom() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  const cutoff = Math.floor(shuffled.length * 0.7);
  const train = shuffled.slice(0, cutoff);
  const test = shuffled.slice(cutoff);

  for (const [key, fn] of Object.entries(varFns)) {
    const res = applyLocked(train, test, fn);
    splitResults[key].push(res);
  }
}

function summarizeSplits(results, name) {
  const rs = results.map(r => r.r).filter(v => isFinite(v));
  const exR2s = results.map(r => r.exR2).filter(v => isFinite(v));
  const rmses = results.map(r => r.rmse).filter(v => isFinite(v));
  const calSlopes = results.map(r => r.calSlope).filter(v => isFinite(v));
  const calOffsets = results.map(r => r.calOffset).filter(v => isFinite(v));

  console.log('  ' + name + ':');
  console.log('    r:         mean=' + mean(rs).toFixed(3) + ' [' + percentile(rs, 0.025).toFixed(3) + ', ' + percentile(rs, 0.975).toFixed(3) + ']');
  console.log('    ext R²:    mean=' + mean(exR2s).toFixed(3) + ' [' + percentile(exR2s, 0.025).toFixed(3) + ', ' + percentile(exR2s, 0.975).toFixed(3) + ']');
  console.log('    RMSE:      mean=' + mean(rmses).toFixed(4));
  console.log('    cal slope: mean=' + mean(calSlopes).toFixed(3) + ' (ideal=1.0)');
  console.log('    cal offset:mean=' + mean(calOffsets).toFixed(4) + ' (ideal=0.0)');

  return {
    r: { mean: parseFloat(mean(rs).toFixed(3)), ci025: parseFloat(percentile(rs, 0.025).toFixed(3)), ci975: parseFloat(percentile(rs, 0.975).toFixed(3)) },
    exR2: { mean: parseFloat(mean(exR2s).toFixed(3)), ci025: parseFloat(percentile(exR2s, 0.025).toFixed(3)), ci975: parseFloat(percentile(exR2s, 0.975).toFixed(3)) },
    rmse: parseFloat(mean(rmses).toFixed(4)),
    calSlope: parseFloat(mean(calSlopes).toFixed(3)),
    calOffset: parseFloat(mean(calOffsets).toFixed(4))
  };
}

const split_fgas = summarizeSplits(splitResults.fgas, 'fgas (PRIMARY)');
const split_sigma0 = summarizeSplits(splitResults.sigma0, 'logSigma0');
const split_compact = summarizeSplits(splitResults.compact, 'logBaryonCompact');
const split_lum = summarizeSplits(splitResults.lum, 'logL36');

const fgasWins = splitResults.fgas.filter((_, i) =>
  splitResults.fgas[i].exR2 >= splitResults.sigma0[i].exR2 &&
  splitResults.fgas[i].exR2 >= splitResults.compact[i].exR2 &&
  splitResults.fgas[i].exR2 >= splitResults.lum[i].exR2
).length;
console.log('\n  fgas wins (best ext R²): ' + fgasWins + '/' + N_SPLITS + ' splits (' + (fgasWins / N_SPLITS * 100).toFixed(1) + '%)');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: DOMAIN TRANSFER — LOW-VFLAT SPARC GALAXIES');
console.log('  Lock on high-Vflat, apply to low-Vflat WITHOUT refit');
console.log('══════════════════════════════════════════════════════════════\n');

function domainTransfer(testData, varFn, trainFit, name) {
  if (testData.length < 5) {
    console.log('  ' + name + ': too few galaxies (N=' + testData.length + ')');
    return null;
  }

  const predictions = testData.map(g => trainFit.beta[0] + trainFit.beta[1] * varFn(g));
  const actual = testData.map(g => g.logOMD);

  const r = pearsonR(predictions, actual);
  const rho = spearmanR(predictions, actual);
  const ss_res = actual.reduce((s, v, i) => s + (v - predictions[i]) ** 2, 0);
  const ss_tot = actual.reduce((s, v) => s + (v - mean(actual)) ** 2, 0);
  const exR2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
  const rmse = Math.sqrt(ss_res / actual.length);

  const calFit = olsRegress(predictions.map(p => [1, p]), actual);

  console.log('  ' + name + ' (N=' + testData.length + '):');
  console.log('    r=' + r.toFixed(3) + ', rho=' + rho.toFixed(3) + ', ext R²=' + exR2.toFixed(3));
  console.log('    RMSE=' + rmse.toFixed(4));
  if (calFit) console.log('    cal slope=' + calFit.beta[1].toFixed(3) + ', offset=' + calFit.beta[0].toFixed(4));

  const refitFit = olsRegress(testData.map(g => [1, varFn(g)]), actual);
  if (refitFit) {
    console.log('    [refit on low-Vflat: slope=' + refitFit.beta[1].toFixed(4) + ' vs locked ' + trainFit.beta[1].toFixed(4) + ']');
  }

  return { n: testData.length, r: parseFloat(r.toFixed(3)), rho: parseFloat(rho.toFixed(3)),
    exR2: parseFloat(exR2.toFixed(3)), rmse: parseFloat(rmse.toFixed(4)),
    calSlope: calFit ? parseFloat(calFit.beta[1].toFixed(3)) : null,
    calOffset: calFit ? parseFloat(calFit.beta[0].toFixed(4)) : null };
}

console.log('  Low-Vflat galaxies: N=' + lowVflat.length);
console.log('  Vflat range: [' + Math.min(...lowVflat.map(g => g.Vflat)).toFixed(1) + ', ' + Math.max(...lowVflat.map(g => g.Vflat)).toFixed(1) + '] km/s\n');

const dt_fgas = domainTransfer(lowVflat, g => g.fgas, fgasModel, 'fgas');
console.log('');
const dt_sigma = domainTransfer(lowVflat, g => g.logSigma0, sigmaModel, 'logSigma0');
console.log('');
const dt_compact = domainTransfer(lowVflat, g => g.logBaryonCompact, compactModel, 'logBaryonCompact');
console.log('');
const dt_lum = domainTransfer(lowVflat, g => g.logL36, lumModel, 'logL36');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: LITTLE THINGS EXTERNAL DATASET (Oh+ 2015)');
console.log('  Truly external: different survey, different telescope');
console.log('══════════════════════════════════════════════════════════════\n');

let ltResult = null;
try {
  const ltTable1 = fs.readFileSync('/tmp/lt_table1.dat', 'utf-8').trim().split('\n');
  const ltTable4 = fs.readFileSync('/tmp/hunter2012_table4.dat', 'utf-8').trim().split('\n');
  const ltTable2 = fs.readFileSync('/tmp/lt_table2.dat', 'utf-8').trim().split('\n');
  const ltRotCurves = fs.readFileSync('/tmp/lt_rotdmbar.dat', 'utf-8').trim().split('\n');
  const ltRotDM = fs.readFileSync('/tmp/lt_rotdm.dat.gz', 'utf-8').trim().split('\n');

  console.log('  LITTLE THINGS data loaded: ' + ltTable1.length + ' galaxies');
} catch (e) {
  console.log('  Note: Missing DM-only file, trying without...');
}

try {
  const ltTable1Lines = fs.readFileSync('/tmp/lt_table1.dat', 'utf-8').trim().split('\n');
  const ltTable2Lines = fs.readFileSync('/tmp/lt_table2.dat', 'utf-8').trim().split('\n');
  const ltRotLines = fs.readFileSync('/tmp/lt_rotdmbar.dat', 'utf-8').trim().split('\n');

  const ltTable4Lines = fs.readFileSync('/tmp/hunter2012_table4.dat', 'utf-8').trim().split('\n');
  const ltMHI = {};
  for (const line of ltTable4Lines) {
    const name = line.substring(4, 13).trim().replace(/\s+/g, '_');
    const logMHI = parseFloat(line.substring(27, 31).trim());
    if (isFinite(logMHI)) ltMHI[name] = Math.pow(10, logMHI);
  }

  const ltGalaxies = {};
  for (const line of ltTable1Lines) {
    const name = line.substring(0, 8).trim().replace(/\s+/g, '_');
    const D = parseFloat(line.substring(32, 36).trim());
    const VMag = parseFloat(line.substring(69, 74).trim());
    const Rd = parseFloat(line.substring(0, 0)); // Not easily parsable from | format
    ltGalaxies[name] = { D, VMag };
  }

  for (const line of ltTable2Lines) {
    const name = line.substring(0, 8).trim().replace(/\s+/g, '_');
    const Rmax = parseFloat(line.substring(9, 13).trim());
    const VRmax = parseFloat(line.substring(19, 24).trim());
    if (ltGalaxies[name]) {
      ltGalaxies[name].Rmax = Rmax;
      ltGalaxies[name].VRmax = VRmax;
    }
  }

  const ltRC = {};
  for (const line of ltRotLines) {
    const name = line.substring(0, 8).trim().replace(/\s+/g, '_');
    const type = line.substring(9, 14).trim();
    if (type !== 'Data') continue;
    const R03 = parseFloat(line.substring(15, 23).trim());
    const V03 = parseFloat(line.substring(24, 34).trim());
    const Rscaled = parseFloat(line.substring(35, 44).trim());
    const Vscaled = parseFloat(line.substring(45, 54).trim());

    if (!isFinite(Rscaled) || !isFinite(Vscaled)) continue;
    if (!ltRC[name]) ltRC[name] = { R03, V03, pts: [] };
    ltRC[name].pts.push({ rScaled: Rscaled, vScaled: Vscaled, r: Rscaled * R03, v: Vscaled * V03 });
  }

  const ltParsedNames = Object.keys(ltRC);
  console.log('  LT parsed RC galaxies: ' + ltParsedNames.length);
  console.log('  LT HI mass available: ' + Object.keys(ltMHI).length + ' galaxies');

  const sparcNames = new Set(Object.keys(table1));
  const ltValid = [];

  for (const name of ltParsedNames) {
    const rc = ltRC[name];
    if (!rc || rc.pts.length < 6) continue;

    const mhi_raw = ltMHI[name];
    if (!mhi_raw || mhi_raw <= 0) continue;

    const gal = ltGalaxies[name];
    if (!gal || !gal.VRmax || gal.VRmax <= 0) continue;

    const cleanName = name.replace(/_/g, '');
    if (sparcNames.has(cleanName)) continue;

    const sorted = [...rc.pts].sort((a, b) => a.r - b.r).filter(p => p.v > 0 && p.r > 0);
    if (sorted.length < 6) continue;

    const Vflat = gal.VRmax;
    const VMag = gal.VMag;
    const L_V = isFinite(VMag) ? Math.pow(10, (4.83 - VMag) / 2.5) : 0;
    const L36_approx = L_V / 2.36e9;
    const MHI_solar = mhi_raw;
    const fgas = MHI_solar / (MHI_solar + L_V * UPSILON_DISK);

    if (!isFinite(fgas) || fgas <= 0 || fgas >= 1) continue;

    const half = Math.floor(sorted.length / 2);
    const outerPts = sorted.slice(half);

    const outerMD = outerPts.map(p => {
      const vBarApprox = p.v * 0.3;
      return (p.v * p.v) / (vBarApprox * vBarApprox);
    });

    ltValid.push({ name, fgas, Vflat, nPts: sorted.length, note: 'LITTLE THINGS (no baryonic decomposition)' });
  }

  console.log('  LT valid (non-overlapping, has MHI): ' + ltValid.length);
  if (ltValid.length >= 5) {
    console.log('  LT fgas values: [' + ltValid.map(g => g.fgas.toFixed(3)).join(', ') + ']');
    console.log('  LT Vflat values: [' + ltValid.map(g => g.Vflat.toFixed(1)).join(', ') + ']');
    console.log('  Note: Cannot compute logOMD without baryonic decomposition.');
    console.log('  LITTLE THINGS provides total V but not separate Vgas, Vdisk, Vbul.');
    console.log('  External replication requires matching target variable definition.');
    ltResult = { status: 'INCOMPLETE', reason: 'No baryonic decomposition in LITTLE THINGS', nGalaxies: ltValid.length };
  } else {
    console.log('  Insufficient non-overlapping galaxies.');
    ltResult = { status: 'INSUFFICIENT_DATA', nGalaxies: ltValid.length };
  }
} catch (e) {
  console.log('  LITTLE THINGS parsing error: ' + e.message);
  ltResult = { status: 'ERROR', error: e.message };
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: FULL SPARC LEAVE-P%-OUT CROSS-VALIDATION');
console.log('  More rigorous than LOO — tests on larger held-out sets');
console.log('══════════════════════════════════════════════════════════════\n');

const lpoPcts = [10, 20, 30, 40, 50];
for (const pct of lpoPcts) {
  const holdoutSize = Math.floor(highVflat.length * pct / 100);
  const nIters = 500;
  const exR2s = [];

  for (let iter = 0; iter < nIters; iter++) {
    const indices = highVflat.map((_, i) => i);
    for (let i = indices.length - 1; i > 0; i--) {
      const j = Math.floor(seededRandom() * (i + 1));
      [indices[i], indices[j]] = [indices[j], indices[i]];
    }
    const testIdx = new Set(indices.slice(0, holdoutSize));
    const train = highVflat.filter((_, i) => !testIdx.has(i));
    const test = highVflat.filter((_, i) => testIdx.has(i));

    const fit = olsRegress(train.map(g => [1, g.fgas]), train.map(g => g.logOMD));
    if (!fit) continue;

    const preds = test.map(g => fit.beta[0] + fit.beta[1] * g.fgas);
    const actuals = test.map(g => g.logOMD);
    const ss_res = actuals.reduce((s, v, i) => s + (v - preds[i]) ** 2, 0);
    const ss_tot = actuals.reduce((s, v) => s + (v - mean(actuals)) ** 2, 0);
    exR2s.push(ss_tot > 0 ? 1 - ss_res / ss_tot : 0);
  }

  console.log('  Leave-' + pct + '%-out (hold ' + holdoutSize + '): ext R² mean=' +
    mean(exR2s).toFixed(3) + ' [' + percentile(exR2s, 0.025).toFixed(3) + ', ' + percentile(exR2s, 0.975).toFixed(3) + ']');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: LOCKED MODEL vs COMPETITORS ON HELD-OUT DATA');
console.log('  Which single variable wins most often outside training set?');
console.log('══════════════════════════════════════════════════════════════\n');

const competitorWins = { fgas: 0, logSigma0: 0, logBaryonCompact: 0, logL36: 0 };
const nComp = 500;

for (let iter = 0; iter < nComp; iter++) {
  const shuffled = [...highVflat];
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(seededRandom() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  const cutoff = Math.floor(shuffled.length * 0.7);
  const train = shuffled.slice(0, cutoff);
  const test = shuffled.slice(cutoff);

  let bestVar = '';
  let bestR2 = -Infinity;

  for (const [key, fn] of Object.entries(varFns)) {
    const fit = olsRegress(train.map(g => [1, fn(g)]), train.map(g => g.logOMD));
    if (!fit) continue;
    const preds = test.map(g => fit.beta[0] + fit.beta[1] * fn(g));
    const actuals = test.map(g => g.logOMD);
    const ss_res = actuals.reduce((s, v, i) => s + (v - preds[i]) ** 2, 0);
    const ss_tot = actuals.reduce((s, v) => s + (v - mean(actuals)) ** 2, 0);
    const exR2 = ss_tot > 0 ? 1 - ss_res / ss_tot : 0;
    if (exR2 > bestR2) { bestR2 = exR2; bestVar = key; }
  }
  if (bestVar) competitorWins[bestVar]++;
}

console.log('  Winner frequency (out of ' + nComp + ' random 70/30 splits):');
for (const [key, wins] of Object.entries(competitorWins)) {
  console.log('    ' + key.padEnd(20) + ': ' + wins + ' wins (' + (wins / nComp * 100).toFixed(1) + '%)');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  OVERALL ASSESSMENT');
console.log('══════════════════════════════════════════════════════════════\n');

const splitR_mean = split_fgas.r.mean;
const splitExR2_mean = split_fgas.exR2.mean;
const dtR = dt_fgas ? dt_fgas.r : NaN;

const strongSuccess = splitExR2_mean > 0.35 && splitR_mean > 0.5 && (dt_fgas && dt_fgas.r > 0.3);
const moderateSuccess = splitExR2_mean > 0.2 && splitR_mean > 0.4;
const failure = splitExR2_mean < 0.1 || splitR_mean < 0.3;

let verdict;
if (strongSuccess) verdict = 'STRONG SUCCESS: fgas transfers robustly as locked predictor';
else if (moderateSuccess) verdict = 'MODERATE SUCCESS: fgas shows transferable signal but with some calibration drift';
else if (failure) verdict = 'FAILURE: fgas does not transfer outside training set';
else verdict = 'MIXED: fgas shows some transferability but results are inconsistent';

console.log('  70/30 split: r=' + splitR_mean + ', ext R²=' + splitExR2_mean);
console.log('  Domain transfer (low-Vflat): r=' + (dt_fgas ? dt_fgas.r : 'N/A'));
console.log('  fgas wins competitor test: ' + competitorWins.fgas + '/' + nComp + ' (' + (competitorWins.fgas / nComp * 100).toFixed(1) + '%)');
console.log('  Calibration slope (ideal=1): ' + split_fgas.calSlope);
console.log('\n  VERDICT: ' + verdict);

const results = {
  phase: 104,
  title: 'Locked External Replication: fgas as single predictor',
  lockedModel: {
    variable: 'fgas',
    intercept: parseFloat(lockedIntercept.toFixed(4)),
    slope: parseFloat(lockedSlope.toFixed(4)),
    trainingR2: parseFloat(fgasModel.r2.toFixed(4)),
    trainingN: highVflat.length
  },
  test1_split7030: {
    nSplits: N_SPLITS,
    fgas: split_fgas,
    logSigma0: split_sigma0,
    logBaryonCompact: split_compact,
    logL36: split_lum,
    fgasWinRate: parseFloat((fgasWins / N_SPLITS * 100).toFixed(1))
  },
  test2_domainTransfer: {
    nLowVflat: lowVflat.length,
    fgas: dt_fgas,
    logSigma0: dt_sigma,
    logBaryonCompact: dt_compact,
    logL36: dt_lum
  },
  test3_littleThings: ltResult,
  test5_competitorKnockout: {
    nSplits: nComp,
    wins: competitorWins,
    fgasWinPct: parseFloat((competitorWins.fgas / nComp * 100).toFixed(1))
  },
  verdict
};

const outPath = path.join(__dirname, '..', 'public', 'phase104-external-replication.json');
fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
console.log('\n  Results saved to: ' + outPath);
console.log('\n  DONE.');
