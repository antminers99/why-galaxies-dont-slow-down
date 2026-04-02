const fs = require('fs');
const path = require('path');

const G_KPC = 4.3009e-6;
const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const A0_KPC = 3.7032;

function mcGaughRAR(gbar) {
  const y = gbar / A0_KPC;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function linReg(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) {
    sxx += (x[i] - mx) ** 2;
    sxy += (x[i] - mx) * (y[i] - my);
  }
  const slope = sxy / sxx;
  const intercept = my - slope * mx;
  const yPred = x.map(xi => intercept + slope * xi);
  let ssRes = 0, ssTot = 0;
  for (let i = 0; i < n; i++) {
    ssRes += (y[i] - yPred[i]) ** 2;
    ssTot += (y[i] - my) ** 2;
  }
  const r2 = ssTot > 0 ? 1 - ssRes / ssTot : 0;
  const r = Math.sign(slope) * Math.sqrt(Math.abs(r2));
  const se = n > 2 ? Math.sqrt(ssRes / (n - 2) / sxx) : Infinity;
  return { slope, intercept, r2, r, se, n, ssRes, ssTot, rmse: Math.sqrt(ssRes / n) };
}

function multiLinReg(X, y) {
  const n = y.length;
  const p = X[0].length;
  const Xt = [];
  for (let j = 0; j < p; j++) {
    Xt.push(X.map(row => row[j]));
  }
  const XtX = [];
  for (let i = 0; i < p; i++) {
    XtX.push([]);
    for (let j = 0; j < p; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += Xt[i][k] * Xt[j][k];
      XtX[i].push(s);
    }
  }
  const Xty = [];
  for (let i = 0; i < p; i++) {
    let s = 0;
    for (let k = 0; k < n; k++) s += Xt[i][k] * y[k];
    Xty.push(s);
  }
  const beta = solveLinSys(XtX, Xty);
  let ssRes = 0;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let ssTot = 0;
  for (let i = 0; i < n; i++) {
    let pred = 0;
    for (let j = 0; j < p; j++) pred += X[i][j] * beta[j];
    ssRes += (y[i] - pred) ** 2;
    ssTot += (y[i] - my) ** 2;
  }
  const r2 = 1 - ssRes / ssTot;
  const adjR2 = 1 - (1 - r2) * (n - 1) / (n - p - 1);
  const rmse = Math.sqrt(ssRes / n);
  const aic = n * Math.log(ssRes / n) + 2 * p;
  const bic = n * Math.log(ssRes / n) + p * Math.log(n);
  return { beta, r2, adjR2, rmse, aic, bic, ssRes, n, p };
}

function solveLinSys(A, b) {
  const n = A.length;
  const M = A.map((row, i) => [...row, b[i]]);
  for (let col = 0; col < n; col++) {
    let maxRow = col;
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(M[row][col]) > Math.abs(M[maxRow][col])) maxRow = row;
    }
    [M[col], M[maxRow]] = [M[maxRow], M[col]];
    if (Math.abs(M[col][col]) < 1e-15) continue;
    for (let row = col + 1; row < n; row++) {
      const f = M[row][col] / M[col][col];
      for (let j = col; j <= n; j++) M[row][j] -= f * M[col][j];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}

function pearsonR(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b) / n;
  const my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return sxy / Math.sqrt(sxx * syy);
}

function rmsScatter(vals) {
  const m = vals.reduce((a, b) => a + b, 0) / vals.length;
  return Math.sqrt(vals.reduce((a, v) => a + (v - m) ** 2, 0) / vals.length);
}

function medAbsDev(vals) {
  const med = percentile(vals, 0.5);
  const devs = vals.map(v => Math.abs(v - med));
  return percentile(devs, 0.5);
}

function percentile(arr, p) {
  const s = [...arr].sort((a, b) => a - b);
  const idx = p * (s.length - 1);
  const lo = Math.floor(idx);
  const hi = Math.ceil(idx);
  return lo === hi ? s[lo] : s[lo] * (hi - idx) + s[hi] * (idx - lo);
}

console.log("=== BIVARIATE COLLAPSE RELATION ANALYSIS ===\n");

const sparcPoints = [];
const sparcDir = '/tmp/rotmod';
const files = fs.readdirSync(sparcDir).filter(f => f.endsWith('.dat'));

const rarJson = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const galaxyLookup = {};
for (const g of rarJson.perGalaxy) {
  galaxyLookup[g.name] = g;
}

let sparcGalCount = 0;
for (const file of files) {
  const name = file.replace('_rotmod.dat', '');
  const lines = fs.readFileSync(path.join(sparcDir, file), 'utf8').trim().split('\n');
  
  const gInfo = galaxyLookup[name];
  if (!gInfo) continue;
  
  const points = [];
  for (const line of lines) {
    if (line.trim().startsWith('#') || line.trim() === '') continue;
    const parts = line.trim().split(/\s+/).map(Number);
    if (parts.length < 6) continue;
    const [r, vobs, evobs, vgas, vdisk, vbul] = parts;
    
    if (!isFinite(r) || !isFinite(vobs) || r <= 0 || vobs <= 0) continue;
    
    const vBarSq = UPSILON_DISK * vdisk * Math.abs(vdisk) + UPSILON_BULGE * vbul * Math.abs(vbul) + vgas * Math.abs(vgas);
    
    const gObs = vobs * vobs / r;
    const gBar = Math.abs(vBarSq) / r;
    
    if (gBar <= 0 || gObs <= 0 || !isFinite(gBar) || !isFinite(gObs)) continue;
    
    const logGobs = Math.log10(gObs);
    const logGbar = Math.log10(gBar);
    
    const gRAR = mcGaughRAR(gBar);
    const deltaRAR = Math.log10(gObs) - Math.log10(gRAR);
    
    if (!isFinite(deltaRAR)) continue;
    
    const sigBar = gInfo.sigma_bar;
    const logSigBar = Math.log10(sigBar);
    
    if (!isFinite(logSigBar)) continue;
    
    points.push({
      galaxy: name,
      dataset: 'SPARC',
      r, vobs, 
      logGobs, logGbar, gObs, gBar,
      deltaRAR, sigBar, logSigBar,
      vmax: gInfo.Vmax,
      logVmax: Math.log10(gInfo.Vmax)
    });
  }
  
  if (points.length >= 3) {
    sparcGalCount++;
    sparcPoints.push(...points);
  }
}
console.log("SPARC: " + sparcGalCount + " galaxies, " + sparcPoints.length + " radial points");

const ltPoints = [];
function parseRotCurve(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const galaxies = {};
  for (const line of lines) {
    const name = line.substring(0, 8).trim();
    const type = line.substring(9, 14).trim();
    if (type !== 'Data') continue;
    const r03 = parseFloat(line.substring(15, 23).trim());
    const v03 = parseFloat(line.substring(24, 34).trim());
    const rScaled = parseFloat(line.substring(35, 44).trim());
    const vScaled = parseFloat(line.substring(45, 54).trim());
    const rKpc = rScaled * r03;
    const vKms = vScaled * v03;
    if (!galaxies[name]) galaxies[name] = [];
    galaxies[name].push({ r: rKpc, v: vKms });
  }
  return galaxies;
}

function parseTable2(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const results = {};
  for (const line of lines) {
    const name = line.substring(0, 8).trim();
    const rmax = parseFloat(line.substring(9, 13).trim());
    const vRmax = parseFloat(line.substring(19, 24).trim());
    const mgasStr = line.substring(139, 145).trim();
    const mstarSEDStr = line.substring(152, 157).trim();
    const mstarKStr = line.substring(146, 151).trim();
    const mgas = mgasStr ? parseFloat(mgasStr) * 1e7 : 0;
    const mstarSED = mstarSEDStr ? parseFloat(mstarSEDStr) * 1e7 : 0;
    const mstarK = mstarKStr ? parseFloat(mstarKStr) * 1e7 : 0;
    const mstar = mstarSED > 0 ? mstarSED : mstarK;
    results[name] = { rmax, vRmax, mgas, mstar, mbar: mstar + 1.33 * mgas };
  }
  return results;
}

function interpolateV(curve, r) {
  if (curve.length === 0) return NaN;
  if (r <= curve[0].r) return curve[0].v * (r / curve[0].r);
  if (r >= curve[curve.length - 1].r) return curve[curve.length - 1].v;
  for (let i = 0; i < curve.length - 1; i++) {
    if (r >= curve[i].r && r <= curve[i + 1].r) {
      const frac = (r - curve[i].r) / (curve[i + 1].r - curve[i].r);
      return curve[i].v + frac * (curve[i + 1].v - curve[i].v);
    }
  }
  return curve[curve.length - 1].v;
}

const rotTotal = parseRotCurve('/tmp/little_things/rotdmbar.dat');
const rotDM = parseRotCurve('/tmp/little_things/rotdm.dat');
const table2 = parseTable2('/tmp/little_things/table2.dat');

let ltGalCount = 0;
for (const name of Object.keys(rotTotal)) {
  if (!rotDM[name] || !table2[name]) continue;
  const t2 = table2[name];
  const curveTotal = rotTotal[name].sort((a, b) => a.r - b.r);
  const curveDM = rotDM[name].sort((a, b) => a.r - b.r);
  if (curveTotal.length < 5 || t2.mbar <= 0) continue;

  const pts = [];
  for (const pt of curveTotal) {
    const r = pt.r;
    const vObs = pt.v;
    const vDM = interpolateV(curveDM, r);
    if (isNaN(vDM) || vObs <= 0 || r <= 0) continue;

    let vBarSq = vObs * vObs - vDM * vDM;
    if (vBarSq < 0) vBarSq = 0;

    const gObs = vObs * vObs / r;
    const gBar = vBarSq / r;
    if (gBar <= 0 || gObs <= 0) continue;

    const gRAR = mcGaughRAR(gBar);
    const deltaRAR = Math.log10(gObs) - Math.log10(gRAR);

    const sigBar = t2.mbar / (Math.PI * r * r);
    const logSigBar = Math.log10(sigBar);

    pts.push({
      galaxy: name,
      dataset: 'LITTLE_THINGS',
      r, vobs: vObs,
      logGobs: Math.log10(gObs), logGbar: Math.log10(gBar),
      gObs, gBar,
      deltaRAR, sigBar, logSigBar,
      vmax: t2.vRmax,
      logVmax: Math.log10(t2.vRmax)
    });
  }
  if (pts.length >= 3) {
    ltGalCount++;
    ltPoints.push(...pts);
  }
}
console.log("LITTLE THINGS: " + ltGalCount + " galaxies, " + ltPoints.length + " radial points\n");

const allPoints = [...sparcPoints, ...ltPoints];
console.log("COMBINED: " + allPoints.length + " total radial points\n");

console.log("========================================");
console.log("TEST 1: STANDARD RAR (1-parameter McGaugh)");
console.log("========================================\n");

function rarResidualStats(points, label) {
  const deltas = points.map(p => p.deltaRAR);
  const rms = rmsScatter(deltas);
  const mad = medAbsDev(deltas);
  const mean = deltas.reduce((a, b) => a + b, 0) / deltas.length;
  console.log(label + ": n=" + deltas.length + " mean(ΔRAR)=" + mean.toFixed(4) + " rms=" + rms.toFixed(4) + " MAD=" + mad.toFixed(4));
  return { rms, mad, mean, n: deltas.length };
}

const sparcRAR = rarResidualStats(sparcPoints, "  SPARC standard RAR");
const ltRAR = rarResidualStats(ltPoints, "  LT    standard RAR");
const allRAR = rarResidualStats(allPoints, "  ALL   standard RAR");

console.log("\n========================================");
console.log("TEST 2: DOES Σ_bar EXPLAIN RAR RESIDUALS?");
console.log("  ΔRAR = a + b·log(Σ_bar)");
console.log("========================================\n");

function testDeltaVsSigma(points, label) {
  const x = points.map(p => p.logSigBar);
  const y = points.map(p => p.deltaRAR);
  const reg = linReg(x, y);
  const significance = Math.abs(reg.slope / reg.se);
  console.log("  " + label + ":");
  console.log("    slope = " + reg.slope.toFixed(5) + " ± " + reg.se.toFixed(5) + "  (" + significance.toFixed(1) + "σ)");
  console.log("    r = " + reg.r.toFixed(4) + "  R² = " + reg.r2.toFixed(4) + "  rmse = " + reg.rmse.toFixed(4));
  return reg;
}

const sparcDeltaSig = testDeltaVsSigma(sparcPoints, "SPARC");
const ltDeltaSig = testDeltaVsSigma(ltPoints, "LITTLE THINGS");
const allDeltaSig = testDeltaVsSigma(allPoints, "COMBINED");

console.log("\n  SIGN CONSISTENCY: SPARC slope " + (sparcDeltaSig.slope < 0 ? "negative" : "positive") + 
  ", LT slope " + (ltDeltaSig.slope < 0 ? "negative" : "positive") +
  " => " + (Math.sign(sparcDeltaSig.slope) === Math.sign(ltDeltaSig.slope) ? "CONSISTENT" : "INCONSISTENT"));

console.log("\n========================================");
console.log("TEST 3: BIVARIATE MODEL COMPARISON");
console.log("  Model A: log(g_obs) = a + b·log(g_bar)                    [standard RAR slope]");
console.log("  Model B: log(g_obs) = a + b·log(g_bar) + c·log(Σ_bar)     [bivariate]");
console.log("  Model C: log(g_obs) = a + b·log(g_bar) + c·log(Σ_bar) + d·log(g_bar)·log(Σ_bar)");
console.log("========================================\n");

function testBivariateModels(points, label) {
  const logGbar = points.map(p => p.logGbar);
  const logGobs = points.map(p => p.logGobs);
  const logSig = points.map(p => p.logSigBar);
  const n = points.length;

  const Xa = logGbar.map(g => [1, g]);
  const modelA = multiLinReg(Xa, logGobs);

  const Xb = logGbar.map((g, i) => [1, g, logSig[i]]);
  const modelB = multiLinReg(Xb, logGobs);

  const Xc = logGbar.map((g, i) => [1, g, logSig[i], g * logSig[i]]);
  const modelC = multiLinReg(Xc, logGobs);

  const fStatAB = ((modelA.ssRes - modelB.ssRes) / 1) / (modelB.ssRes / (n - 3));
  const fStatAC = ((modelA.ssRes - modelC.ssRes) / 2) / (modelC.ssRes / (n - 4));
  const fStatBC = ((modelB.ssRes - modelC.ssRes) / 1) / (modelC.ssRes / (n - 4));

  console.log("  " + label + " (n=" + n + "):");
  console.log("    Model A (g_bar only):           R²=" + modelA.r2.toFixed(5) + "  adjR²=" + modelA.adjR2.toFixed(5) + "  RMSE=" + modelA.rmse.toFixed(5) + "  AIC=" + modelA.aic.toFixed(1) + "  BIC=" + modelA.bic.toFixed(1));
  console.log("      coefficients: a=" + modelA.beta[0].toFixed(4) + " b=" + modelA.beta[1].toFixed(4));
  console.log("    Model B (+Σ_bar):               R²=" + modelB.r2.toFixed(5) + "  adjR²=" + modelB.adjR2.toFixed(5) + "  RMSE=" + modelB.rmse.toFixed(5) + "  AIC=" + modelB.aic.toFixed(1) + "  BIC=" + modelB.bic.toFixed(1));
  console.log("      coefficients: a=" + modelB.beta[0].toFixed(4) + " b=" + modelB.beta[1].toFixed(4) + " c(Σ)=" + modelB.beta[2].toFixed(5));
  console.log("    Model C (+interaction):         R²=" + modelC.r2.toFixed(5) + "  adjR²=" + modelC.adjR2.toFixed(5) + "  RMSE=" + modelC.rmse.toFixed(5) + "  AIC=" + modelC.aic.toFixed(1) + "  BIC=" + modelC.bic.toFixed(1));
  console.log("      coefficients: a=" + modelC.beta[0].toFixed(4) + " b=" + modelC.beta[1].toFixed(4) + " c(Σ)=" + modelC.beta[2].toFixed(5) + " d(g·Σ)=" + modelC.beta[3].toFixed(6));
  console.log("    F-test A vs B: F=" + fStatAB.toFixed(2) + " (Σ_bar adds predictive power: " + (fStatAB > 6.63 ? "YES p<0.01" : fStatAB > 3.84 ? "marginal p<0.05" : "NO") + ")");
  console.log("    F-test A vs C: F=" + fStatAC.toFixed(2));
  console.log("    F-test B vs C: F=" + fStatBC.toFixed(2) + " (interaction term needed: " + (fStatBC > 6.63 ? "YES" : "NO") + ")");
  
  const deltaAIC_AB = modelA.aic - modelB.aic;
  const deltaBIC_AB = modelA.bic - modelB.bic;
  console.log("    ΔAIC(A-B)=" + deltaAIC_AB.toFixed(1) + " ΔBIC(A-B)=" + deltaBIC_AB.toFixed(1) + " (positive = B wins)");

  return { modelA, modelB, modelC, fStatAB, deltaAIC_AB, deltaBIC_AB };
}

const sparcBV = testBivariateModels(sparcPoints, "SPARC");
console.log("");
const ltBV = testBivariateModels(ltPoints, "LITTLE THINGS");
console.log("");
const allBV = testBivariateModels(allPoints, "COMBINED");

console.log("\n========================================");
console.log("TEST 4: CROSS-VALIDATION (THE CRITICAL TEST)");
console.log("  Train on SPARC → predict LITTLE THINGS, and vice versa");
console.log("========================================\n");

function trainAndPredict(trainPts, testPts, trainLabel, testLabel) {
  const trainLogGbar = trainPts.map(p => p.logGbar);
  const trainLogGobs = trainPts.map(p => p.logGobs);
  const trainLogSig = trainPts.map(p => p.logSigBar);

  const XaT = trainLogGbar.map(g => [1, g]);
  const mA = multiLinReg(XaT, trainLogGobs);

  const XbT = trainLogGbar.map((g, i) => [1, g, trainLogSig[i]]);
  const mB = multiLinReg(XbT, trainLogGobs);

  let ssResA = 0, ssResB = 0;
  const testN = testPts.length;
  const residA = [], residB = [];

  for (const tp of testPts) {
    const predA = mA.beta[0] + mA.beta[1] * tp.logGbar;
    const predB = mB.beta[0] + mB.beta[1] * tp.logGbar + mB.beta[2] * tp.logSigBar;
    const rA = tp.logGobs - predA;
    const rB = tp.logGobs - predB;
    ssResA += rA * rA;
    ssResB += rB * rB;
    residA.push(rA);
    residB.push(rB);
  }

  const rmseA = Math.sqrt(ssResA / testN);
  const rmseB = Math.sqrt(ssResB / testN);
  const improvement = ((rmseA - rmseB) / rmseA * 100);

  console.log("  Train: " + trainLabel + " (" + trainPts.length + " pts) → Test: " + testLabel + " (" + testN + " pts)");
  console.log("    Model A (g_bar only) test RMSE: " + rmseA.toFixed(5));
  console.log("    Model B (+Σ_bar)     test RMSE: " + rmseB.toFixed(5));
  console.log("    Improvement: " + improvement.toFixed(2) + "%");
  console.log("    Trained coefficients: a=" + mB.beta[0].toFixed(4) + " b=" + mB.beta[1].toFixed(4) + " c(Σ)=" + mB.beta[2].toFixed(5));
  console.log("    Σ_bar coefficient sign: " + (mB.beta[2] < 0 ? "NEGATIVE" : "POSITIVE"));
  
  return { rmseA, rmseB, improvement, betaB: mB.beta, residA, residB };
}

const cv1 = trainAndPredict(sparcPoints, ltPoints, "SPARC", "LITTLE THINGS");
console.log("");
const cv2 = trainAndPredict(ltPoints, sparcPoints, "LITTLE THINGS", "SPARC");

console.log("\n  Cross-validation verdict:");
const signConsistent = Math.sign(cv1.betaB[2]) === Math.sign(cv2.betaB[2]);
console.log("    Σ_bar coefficient sign consistent: " + (signConsistent ? "YES" : "NO"));
console.log("    Both directions improve: " + (cv1.improvement > 0 && cv2.improvement > 0 ? "YES" : "NO"));

console.log("\n========================================");
console.log("TEST 5: COLLAPSE RELATION — RESCALED VARIABLES");
console.log("  Searching for combination that minimizes scatter");
console.log("========================================\n");

function testCollapse(points, formula, label) {
  const vals = points.map(formula).filter(v => isFinite(v) && !isNaN(v));
  if (vals.length < 10) return null;
  const rms = rmsScatter(vals);
  const mad = medAbsDev(vals);
  return { label, rms, mad, n: vals.length };
}

const collapseTests = [
  { label: "ΔRAR (standard)", fn: p => p.deltaRAR },
  { label: "ΔRAR - 0.05·log(Σ)", fn: p => p.deltaRAR - 0.05 * p.logSigBar },
  { label: "ΔRAR - 0.10·log(Σ)", fn: p => p.deltaRAR - 0.10 * p.logSigBar },
  { label: "ΔRAR - 0.15·log(Σ)", fn: p => p.deltaRAR - 0.15 * p.logSigBar },
  { label: "ΔRAR + 0.05·log(Σ)", fn: p => p.deltaRAR + 0.05 * p.logSigBar },
  { label: "ΔRAR + 0.10·log(Σ)", fn: p => p.deltaRAR + 0.10 * p.logSigBar },
];

console.log("  Searching for Σ correction that minimizes combined scatter:\n");
let bestCollapse = null;

for (let coeff = -0.30; coeff <= 0.30; coeff += 0.005) {
  const c = coeff;
  const sparcVals = sparcPoints.map(p => p.deltaRAR - c * p.logSigBar).filter(v => isFinite(v));
  const ltVals = ltPoints.map(p => p.deltaRAR - c * p.logSigBar).filter(v => isFinite(v));
  const allVals = [...sparcVals, ...ltVals];
  
  const sparcRMS = rmsScatter(sparcVals);
  const ltRMS = rmsScatter(ltVals);
  const allRMS = rmsScatter(allVals);
  
  const sparcMean = sparcVals.reduce((a, b) => a + b, 0) / sparcVals.length;
  const ltMean = ltVals.reduce((a, b) => a + b, 0) / ltVals.length;
  const meanDiff = Math.abs(sparcMean - ltMean);
  
  if (!bestCollapse || allRMS < bestCollapse.allRMS) {
    bestCollapse = { coeff: c, allRMS, sparcRMS, ltRMS, meanDiff, sparcMean, ltMean };
  }
}

const baselineSparcRMS = rmsScatter(sparcPoints.map(p => p.deltaRAR));
const baselineLtRMS = rmsScatter(ltPoints.map(p => p.deltaRAR));
const baselineAllRMS = rmsScatter([...sparcPoints, ...ltPoints].map(p => p.deltaRAR));
const baselineMeanDiff = Math.abs(
  sparcPoints.map(p => p.deltaRAR).reduce((a, b) => a + b, 0) / sparcPoints.length -
  ltPoints.map(p => p.deltaRAR).reduce((a, b) => a + b, 0) / ltPoints.length
);

console.log("  Baseline (no Σ correction):");
console.log("    SPARC rms = " + baselineSparcRMS.toFixed(5) + "  LT rms = " + baselineLtRMS.toFixed(5) + "  Combined rms = " + baselineAllRMS.toFixed(5));
console.log("    Mean gap (SPARC vs LT): " + baselineMeanDiff.toFixed(5));
console.log("");
console.log("  Best Σ correction: ΔRAR_corr = ΔRAR - " + bestCollapse.coeff.toFixed(3) + "·log(Σ_bar)");
console.log("    SPARC rms = " + bestCollapse.sparcRMS.toFixed(5) + "  LT rms = " + bestCollapse.ltRMS.toFixed(5) + "  Combined rms = " + bestCollapse.allRMS.toFixed(5));
console.log("    Mean gap (SPARC vs LT): " + bestCollapse.meanDiff.toFixed(5));
console.log("    Scatter reduction: " + ((baselineAllRMS - bestCollapse.allRMS) / baselineAllRMS * 100).toFixed(2) + "%");
console.log("    Mean gap reduction: " + ((baselineMeanDiff - bestCollapse.meanDiff) / baselineMeanDiff * 100).toFixed(2) + "%");

console.log("\n========================================");
console.log("TEST 6: STABILITY UNDER MULTIPLE DEFINITIONS");
console.log("  Testing with different Υ_disk, radial windows, Σ definitions");
console.log("========================================\n");

function computeSPARCwithParams(upsilon, rMinFrac, rMaxFrac) {
  const pts = [];
  for (const file of files) {
    const name = file.replace('_rotmod.dat', '');
    const gInfo = galaxyLookup[name];
    if (!gInfo) continue;
    const lines = fs.readFileSync(path.join(sparcDir, file), 'utf8').trim().split('\n');
    const allR = lines.map(l => parseFloat(l.trim().split(/\s+/)[0])).filter(r => r > 0);
    const rMax = Math.max(...allR);
    
    for (const line of lines) {
      if (line.trim().startsWith('#') || line.trim() === '') continue;
      const parts = line.trim().split(/\s+/).map(Number);
      if (parts.length < 6) continue;
      const [r, vobs, evobs, vgas, vdisk, vbul] = parts;
      if (!isFinite(r) || !isFinite(vobs) || r <= 0 || vobs <= 0) continue;
      if (r < rMinFrac * rMax || r > rMaxFrac * rMax) continue;
      
      const vBarSq = upsilon * vdisk * Math.abs(vdisk) + UPSILON_BULGE * vbul * Math.abs(vbul) + vgas * Math.abs(vgas);
      const gObs = vobs * vobs / r;
      const gBar = Math.abs(vBarSq) / r;
      if (gBar <= 0 || gObs <= 0 || !isFinite(gBar) || !isFinite(gObs)) continue;
      
      const gRAR = mcGaughRAR(gBar);
      const deltaRAR = Math.log10(gObs) - Math.log10(gRAR);
      if (!isFinite(deltaRAR)) continue;
      const sigBar = gInfo.sigma_bar;
      const logSigBar = Math.log10(sigBar);
      if (!isFinite(logSigBar)) continue;
      
      pts.push({ logGbar: Math.log10(gBar), logGobs: Math.log10(gObs), deltaRAR, logSigBar });
    }
  }
  return pts;
}

const paramGrid = [
  { upsilon: 0.3, rMin: 0, rMax: 1.0, label: "Υ=0.3, full range" },
  { upsilon: 0.5, rMin: 0, rMax: 1.0, label: "Υ=0.5, full range (baseline)" },
  { upsilon: 0.7, rMin: 0, rMax: 1.0, label: "Υ=0.7, full range" },
  { upsilon: 0.5, rMin: 0.2, rMax: 0.8, label: "Υ=0.5, inner 20-80%" },
  { upsilon: 0.5, rMin: 0.4, rMax: 1.0, label: "Υ=0.5, outer 40-100%" },
  { upsilon: 0.5, rMin: 0, rMax: 0.5, label: "Υ=0.5, inner 0-50%" },
  { upsilon: 0.3, rMin: 0.2, rMax: 0.8, label: "Υ=0.3, inner 20-80%" },
  { upsilon: 0.7, rMin: 0.2, rMax: 0.8, label: "Υ=0.7, inner 20-80%" },
];

console.log("  Testing SPARC ΔRAR vs log(Σ_bar) under varied definitions:\n");
console.log("  " + "Definition".padEnd(35) + "N".padStart(6) + "  slope".padStart(10) + "  ±se".padStart(8) + "  σ-sign".padStart(8) + "  sign");

const stabilityResults = [];
for (const p of paramGrid) {
  const pts = computeSPARCwithParams(p.upsilon, p.rMin, p.rMax);
  if (pts.length < 50) continue;
  const x = pts.map(p => p.logSigBar);
  const y = pts.map(p => p.deltaRAR);
  const reg = linReg(x, y);
  const sig = Math.abs(reg.slope / reg.se);
  const sign = reg.slope < 0 ? "NEG" : "POS";
  console.log("  " + p.label.padEnd(35) + String(pts.length).padStart(6) + ("  " + reg.slope.toFixed(5)).padStart(10) + ("  " + reg.se.toFixed(5)).padStart(8) + ("  " + sig.toFixed(1)).padStart(8) + "  " + sign);
  stabilityResults.push({ ...p, slope: reg.slope, se: reg.se, sig, n: pts.length, sign });
}

const allNeg = stabilityResults.every(s => s.slope < 0);
const allPos = stabilityResults.every(s => s.slope > 0);
console.log("\n  Sign stability: " + (allNeg ? "ALL NEGATIVE — ROBUST" : allPos ? "ALL POSITIVE — ROBUST" : "MIXED SIGNS — FRAGILE"));
const meanSig = stabilityResults.reduce((a, s) => a + s.sig, 0) / stabilityResults.length;
console.log("  Mean significance: " + meanSig.toFixed(1) + "σ across " + stabilityResults.length + " definitions");

console.log("\n========================================");
console.log("TEST 7: BIVARIATE FORMULA SEARCH");
console.log("  g_obs = g_bar / (1 - exp(-sqrt(g_bar/a0))) · (Σ_bar/Σ_0)^γ");
console.log("  Finding optimal γ and Σ_0");
console.log("========================================\n");

let bestGamma = null;
const logSig0 = 7.5;

for (let gamma = -0.30; gamma <= 0.30; gamma += 0.005) {
  const residuals = allPoints.map(p => {
    const gRAR = mcGaughRAR(p.gBar);
    const correction = Math.pow(10, gamma * (p.logSigBar - logSig0));
    const gPred = gRAR * correction;
    return Math.log10(p.gObs) - Math.log10(gPred);
  }).filter(v => isFinite(v));
  
  const rms = rmsScatter(residuals);
  if (!bestGamma || rms < bestGamma.rms) {
    bestGamma = { gamma, rms, n: residuals.length };
  }
}

const baseResiduals = allPoints.map(p => {
  const gRAR = mcGaughRAR(p.gBar);
  return Math.log10(p.gObs) - Math.log10(gRAR);
}).filter(v => isFinite(v));
const baseRMS = rmsScatter(baseResiduals);

console.log("  Standard RAR (γ=0):  rms = " + baseRMS.toFixed(5));
console.log("  Best bivariate:      rms = " + bestGamma.rms.toFixed(5) + "  γ = " + bestGamma.gamma.toFixed(3));
console.log("  Scatter reduction:   " + ((baseRMS - bestGamma.rms) / baseRMS * 100).toFixed(2) + "%");
console.log("  Formula: g_obs = RAR(g_bar) · (Σ_bar / 10^" + logSig0 + ")^" + bestGamma.gamma.toFixed(3));

const sparcResidBase = sparcPoints.map(p => {
  const gRAR = mcGaughRAR(p.gBar);
  return Math.log10(p.gObs) - Math.log10(gRAR);
}).filter(v => isFinite(v));
const sparcResidCorr = sparcPoints.map(p => {
  const gRAR = mcGaughRAR(p.gBar);
  const correction = Math.pow(10, bestGamma.gamma * (p.logSigBar - logSig0));
  return Math.log10(p.gObs) - Math.log10(gRAR * correction);
}).filter(v => isFinite(v));
const ltResidBase = ltPoints.map(p => {
  const gRAR = mcGaughRAR(p.gBar);
  return Math.log10(p.gObs) - Math.log10(gRAR);
}).filter(v => isFinite(v));
const ltResidCorr = ltPoints.map(p => {
  const gRAR = mcGaughRAR(p.gBar);
  const correction = Math.pow(10, bestGamma.gamma * (p.logSigBar - logSig0));
  return Math.log10(p.gObs) - Math.log10(gRAR * correction);
}).filter(v => isFinite(v));

console.log("\n  Per-dataset check:");
console.log("    SPARC:  base rms=" + rmsScatter(sparcResidBase).toFixed(5) + "  corrected rms=" + rmsScatter(sparcResidCorr).toFixed(5) + "  Δ=" + ((rmsScatter(sparcResidBase) - rmsScatter(sparcResidCorr)) / rmsScatter(sparcResidBase) * 100).toFixed(2) + "%");
console.log("    LT:     base rms=" + rmsScatter(ltResidBase).toFixed(5) + "  corrected rms=" + rmsScatter(ltResidCorr).toFixed(5) + "  Δ=" + ((rmsScatter(ltResidBase) - rmsScatter(ltResidCorr)) / rmsScatter(ltResidBase) * 100).toFixed(2) + "%");

console.log("\n========================================");
console.log("TEST 8: GALAXY-LEVEL AGGREGATION");
console.log("  Per-galaxy mean ΔRAR vs mean log(Σ_bar)");
console.log("========================================\n");

function galaxyLevelAnalysis(points, label) {
  const byGal = {};
  for (const p of points) {
    if (!byGal[p.galaxy]) byGal[p.galaxy] = { deltas: [], sigs: [] };
    byGal[p.galaxy].deltas.push(p.deltaRAR);
    byGal[p.galaxy].sigs.push(p.logSigBar);
  }
  const galNames = Object.keys(byGal);
  const galDeltas = galNames.map(n => byGal[n].deltas.reduce((a, b) => a + b) / byGal[n].deltas.length);
  const galSigs = galNames.map(n => byGal[n].sigs.reduce((a, b) => a + b) / byGal[n].sigs.length);
  
  const reg = linReg(galSigs, galDeltas);
  const sig = Math.abs(reg.slope / reg.se);
  console.log("  " + label + " (" + galNames.length + " galaxies):");
  console.log("    slope = " + reg.slope.toFixed(5) + " ± " + reg.se.toFixed(5) + "  (" + sig.toFixed(1) + "σ)");
  console.log("    r = " + reg.r.toFixed(4) + "  R² = " + reg.r2.toFixed(4));
  return { reg, sig, n: galNames.length, galDeltas, galSigs };
}

const sparcGalLevel = galaxyLevelAnalysis(sparcPoints, "SPARC galaxy-level");
const ltGalLevel = galaxyLevelAnalysis(ltPoints, "LT galaxy-level");

const combinedGalDeltas = [...sparcGalLevel.galDeltas, ...ltGalLevel.galDeltas];
const combinedGalSigs = [...sparcGalLevel.galSigs, ...ltGalLevel.galSigs];
const combinedGalReg = linReg(combinedGalSigs, combinedGalDeltas);
const combinedGalSig = Math.abs(combinedGalReg.slope / combinedGalReg.se);
console.log("  COMBINED galaxy-level (" + combinedGalDeltas.length + " galaxies):");
console.log("    slope = " + combinedGalReg.slope.toFixed(5) + " ± " + combinedGalReg.se.toFixed(5) + "  (" + combinedGalSig.toFixed(1) + "σ)");
console.log("    r = " + combinedGalReg.r.toFixed(4) + "  R² = " + combinedGalReg.r2.toFixed(4));

console.log("\n\n========================================");
console.log("FINAL VERDICT");
console.log("========================================\n");

const bivariateImprovesSPARC = sparcBV.fStatAB > 6.63;
const bivariateImprovesLT = ltBV.fStatAB > 6.63;
const crossValBothImprove = cv1.improvement > 0 && cv2.improvement > 0;
const signStable = allNeg || allPos;
const crossValSignConsistent = signConsistent;

console.log("  1. Σ_bar improves RAR fit (F-test p<0.01)?");
console.log("     SPARC:         " + (bivariateImprovesSPARC ? "YES (F=" + sparcBV.fStatAB.toFixed(1) + ")" : "NO (F=" + sparcBV.fStatAB.toFixed(1) + ")"));
console.log("     LITTLE THINGS: " + (bivariateImprovesLT ? "YES (F=" + ltBV.fStatAB.toFixed(1) + ")" : "NO (F=" + ltBV.fStatAB.toFixed(1) + ")"));
console.log("  2. Cross-validated (train→test improvement)?");
console.log("     SPARC→LT:      " + (cv1.improvement > 0 ? "YES (" + cv1.improvement.toFixed(1) + "%)" : "NO (" + cv1.improvement.toFixed(1) + "%)"));
console.log("     LT→SPARC:      " + (cv2.improvement > 0 ? "YES (" + cv2.improvement.toFixed(1) + "%)" : "NO (" + cv2.improvement.toFixed(1) + "%)"));
console.log("  3. Σ_bar coefficient sign consistent both directions? " + (crossValSignConsistent ? "YES" : "NO"));
console.log("  4. Sign stable under all definitions? " + (signStable ? "YES" : "NO") + " (mean " + meanSig.toFixed(1) + "σ)");
console.log("  5. Collapse reduces combined scatter? " + ((baselineAllRMS - bestCollapse.allRMS) / baselineAllRMS * 100).toFixed(1) + "%");
console.log("  6. Best bivariate formula: g_obs = RAR(g_bar) · (Σ_bar/10^" + logSig0 + ")^" + bestGamma.gamma.toFixed(3));

const treasure = bivariateImprovesSPARC && crossValBothImprove && crossValSignConsistent && signStable;
console.log("\n  === TREASURE? " + (treasure ? "YES — Bivariate relation is robust, cross-validated, and stable" : "FURTHER INVESTIGATION NEEDED — see details above") + " ===");

const output = {
  timestamp: new Date().toISOString(),
  datasets: {
    sparc: { nGalaxies: sparcGalCount, nPoints: sparcPoints.length },
    littleThings: { nGalaxies: ltGalCount, nPoints: ltPoints.length },
    combined: { nGalaxies: sparcGalCount + ltGalCount, nPoints: allPoints.length }
  },
  test1_standardRAR: {
    sparc: sparcRAR,
    littleThings: ltRAR,
    combined: allRAR
  },
  test2_deltaVsSigma: {
    sparc: { slope: sparcDeltaSig.slope, se: sparcDeltaSig.se, r: sparcDeltaSig.r, r2: sparcDeltaSig.r2, n: sparcDeltaSig.n },
    littleThings: { slope: ltDeltaSig.slope, se: ltDeltaSig.se, r: ltDeltaSig.r, r2: ltDeltaSig.r2, n: ltDeltaSig.n },
    combined: { slope: allDeltaSig.slope, se: allDeltaSig.se, r: allDeltaSig.r, r2: allDeltaSig.r2, n: allDeltaSig.n },
    signConsistent: Math.sign(sparcDeltaSig.slope) === Math.sign(ltDeltaSig.slope)
  },
  test3_bivariateModels: {
    sparc: {
      modelA: { r2: sparcBV.modelA.r2, adjR2: sparcBV.modelA.adjR2, rmse: sparcBV.modelA.rmse, aic: sparcBV.modelA.aic, bic: sparcBV.modelA.bic },
      modelB: { r2: sparcBV.modelB.r2, adjR2: sparcBV.modelB.adjR2, rmse: sparcBV.modelB.rmse, aic: sparcBV.modelB.aic, bic: sparcBV.modelB.bic, sigmaCoeff: sparcBV.modelB.beta[2] },
      fStatAB: sparcBV.fStatAB,
      deltaAIC: sparcBV.deltaAIC_AB,
      deltaBIC: sparcBV.deltaBIC_AB
    },
    littleThings: {
      modelA: { r2: ltBV.modelA.r2, adjR2: ltBV.modelA.adjR2, rmse: ltBV.modelA.rmse, aic: ltBV.modelA.aic, bic: ltBV.modelA.bic },
      modelB: { r2: ltBV.modelB.r2, adjR2: ltBV.modelB.adjR2, rmse: ltBV.modelB.rmse, aic: ltBV.modelB.aic, bic: ltBV.modelB.bic, sigmaCoeff: ltBV.modelB.beta[2] },
      fStatAB: ltBV.fStatAB,
      deltaAIC: ltBV.deltaAIC_AB,
      deltaBIC: ltBV.deltaBIC_AB
    },
    combined: {
      modelA: { r2: allBV.modelA.r2, adjR2: allBV.modelA.adjR2, rmse: allBV.modelA.rmse },
      modelB: { r2: allBV.modelB.r2, adjR2: allBV.modelB.adjR2, rmse: allBV.modelB.rmse, sigmaCoeff: allBV.modelB.beta[2] },
      fStatAB: allBV.fStatAB,
      deltaAIC: allBV.deltaAIC_AB
    }
  },
  test4_crossValidation: {
    sparcToLT: { rmseA: cv1.rmseA, rmseB: cv1.rmseB, improvement: cv1.improvement, sigmaCoeff: cv1.betaB[2] },
    ltToSPARC: { rmseA: cv2.rmseA, rmseB: cv2.rmseB, improvement: cv2.improvement, sigmaCoeff: cv2.betaB[2] },
    signConsistent: signConsistent,
    bothImprove: cv1.improvement > 0 && cv2.improvement > 0
  },
  test5_collapse: {
    baseline: { rms: baselineAllRMS, meanGap: baselineMeanDiff },
    best: { coeff: bestCollapse.coeff, rms: bestCollapse.allRMS, meanGap: bestCollapse.meanDiff },
    scatterReduction: (baselineAllRMS - bestCollapse.allRMS) / baselineAllRMS * 100,
    meanGapReduction: (baselineMeanDiff - bestCollapse.meanDiff) / baselineMeanDiff * 100
  },
  test6_stability: {
    results: stabilityResults.map(s => ({ label: s.label, slope: s.slope, se: s.se, sig: s.sig, n: s.n, sign: s.sign })),
    allSameSign: allNeg || allPos,
    signDirection: allNeg ? "negative" : allPos ? "positive" : "mixed",
    meanSignificance: meanSig
  },
  test7_bivariateFormula: {
    bestGamma: bestGamma.gamma,
    sigmaRef: Math.pow(10, logSig0),
    baseRMS: baseRMS,
    correctedRMS: bestGamma.rms,
    scatterReduction: (baseRMS - bestGamma.rms) / baseRMS * 100,
    formula: "g_obs = RAR(g_bar) * (Sigma_bar / " + Math.pow(10, logSig0).toExponential(1) + ")^" + bestGamma.gamma.toFixed(3),
    perDataset: {
      sparc: { baseRMS: rmsScatter(sparcResidBase), correctedRMS: rmsScatter(sparcResidCorr) },
      littleThings: { baseRMS: rmsScatter(ltResidBase), correctedRMS: rmsScatter(ltResidCorr) }
    }
  },
  test8_galaxyLevel: {
    sparc: { slope: sparcGalLevel.reg.slope, se: sparcGalLevel.reg.se, sig: sparcGalLevel.sig, r: sparcGalLevel.reg.r, n: sparcGalLevel.n },
    littleThings: { slope: ltGalLevel.reg.slope, se: ltGalLevel.reg.se, sig: ltGalLevel.sig, r: ltGalLevel.reg.r, n: ltGalLevel.n },
    combined: { slope: combinedGalReg.slope, se: combinedGalReg.se, sig: combinedGalSig, r: combinedGalReg.r, n: combinedGalDeltas.length }
  },
  verdict: {
    sigmaImprovesRAR_SPARC: bivariateImprovesSPARC,
    sigmaImprovesRAR_LT: bivariateImprovesLT,
    crossValidated: crossValBothImprove,
    signConsistent: crossValSignConsistent,
    definitionStable: signStable,
    isTreasure: treasure
  },
  pointData: {
    sparc: sparcPoints.map(p => ({
      galaxy: p.galaxy, r: p.r, logGobs: p.logGobs, logGbar: p.logGbar,
      deltaRAR: p.deltaRAR, logSigBar: p.logSigBar
    })),
    littleThings: ltPoints.map(p => ({
      galaxy: p.galaxy, r: p.r, logGobs: p.logGobs, logGbar: p.logGbar,
      deltaRAR: p.deltaRAR, logSigBar: p.logSigBar
    }))
  }
};

const outPath = path.join(__dirname, '..', 'public', 'bivariate-collapse.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log("\nResults saved to: " + outPath);
console.log("File size: " + (fs.statSync(outPath).size / 1024 / 1024).toFixed(1) + " MB");
