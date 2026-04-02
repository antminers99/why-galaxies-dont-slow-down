const fs = require('fs');
const path = require('path');

const A0 = 3.7032;
const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function rmsScatter(vals) {
  const m = vals.reduce((a, b) => a + b, 0) / vals.length;
  return Math.sqrt(vals.reduce((a, v) => a + (v - m) ** 2, 0) / vals.length);
}

function medAbsDev(vals) {
  const s = [...vals].sort((a, b) => a - b);
  const med = s[Math.floor(s.length / 2)];
  const devs = vals.map(v => Math.abs(v - med)).sort((a, b) => a - b);
  return devs[Math.floor(devs.length / 2)];
}

function stdRAR(gbar) {
  const y = gbar / A0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function solveLS(A, b) {
  const n = A.length, p = A[0].length;
  const Xt = [];
  for (let j = 0; j < p; j++) Xt.push(A.map(row => row[j]));
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
    for (let k = 0; k < n; k++) s += Xt[i][k] * b[k];
    Xty.push(s);
  }
  const M = XtX.map((row, i) => [...row, Xty[i]]);
  for (let col = 0; col < p; col++) {
    let maxR = col;
    for (let row = col + 1; row < p; row++) if (Math.abs(M[row][col]) > Math.abs(M[maxR][col])) maxR = row;
    [M[col], M[maxR]] = [M[maxR], M[col]];
    if (Math.abs(M[col][col]) < 1e-15) continue;
    for (let row = col + 1; row < p; row++) {
      const f = M[row][col] / M[col][col];
      for (let j = col; j <= p; j++) M[row][j] -= f * M[col][j];
    }
  }
  const x = new Array(p).fill(0);
  for (let i = p - 1; i >= 0; i--) {
    x[i] = M[i][p];
    for (let j = i + 1; j < p; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}

console.log("=== THE HONEST TEST: Is Sigma_bar actually helping? ===\n");
console.log("Critical question: How much improvement comes from BETTER g_bar fitting");
console.log("vs how much actually comes from ADDING Sigma_bar?\n");

const rarJson = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const galaxyLookup = {};
for (const g of rarJson.perGalaxy) galaxyLookup[g.name] = g;

const sparcPoints = [];
const sparcDir = '/tmp/rotmod';
const files = fs.readdirSync(sparcDir).filter(f => f.endsWith('.dat'));
for (const file of files) {
  const name = file.replace('_rotmod.dat', '');
  const gInfo = galaxyLookup[name];
  if (!gInfo) continue;
  const lines = fs.readFileSync(path.join(sparcDir, file), 'utf8').trim().split('\n');
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
    const sigBar = gInfo.sigma_bar;
    if (!isFinite(sigBar) || sigBar <= 0) continue;
    sparcPoints.push({ galaxy: name, gObs, gBar, sigBar, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), logSig: Math.log10(sigBar) });
  }
}

function parseRotCurve(fp) {
  const lines = fs.readFileSync(fp, 'utf8').trim().split('\n');
  const g = {};
  for (const line of lines) {
    const nm = line.substring(0, 8).trim();
    const tp = line.substring(9, 14).trim();
    if (tp !== 'Data') continue;
    if (!g[nm]) g[nm] = [];
    g[nm].push({ r: parseFloat(line.substring(35, 44).trim()) * parseFloat(line.substring(15, 23).trim()), v: parseFloat(line.substring(45, 54).trim()) * parseFloat(line.substring(24, 34).trim()) });
  }
  return g;
}
function parseT2(fp) {
  const lines = fs.readFileSync(fp, 'utf8').trim().split('\n');
  const r = {};
  for (const line of lines) {
    const nm = line.substring(0, 8).trim();
    const mgS = line.substring(139, 145).trim();
    const msS = line.substring(152, 157).trim();
    const mkS = line.substring(146, 151).trim();
    const mg = mgS ? parseFloat(mgS) * 1e7 : 0;
    const ms = msS ? parseFloat(msS) * 1e7 : (mkS ? parseFloat(mkS) * 1e7 : 0);
    r[nm] = { mbar: ms + 1.33 * mg };
  }
  return r;
}
function interpV(c, r) {
  if (c.length === 0) return NaN;
  if (r <= c[0].r) return c[0].v * (r / c[0].r);
  if (r >= c[c.length - 1].r) return c[c.length - 1].v;
  for (let i = 0; i < c.length - 1; i++) {
    if (r >= c[i].r && r <= c[i + 1].r) return c[i].v + (r - c[i].r) / (c[i + 1].r - c[i].r) * (c[i + 1].v - c[i].v);
  }
  return c[c.length - 1].v;
}

const rotT = parseRotCurve('/tmp/little_things/rotdmbar.dat');
const rotD = parseRotCurve('/tmp/little_things/rotdm.dat');
const t2 = parseT2('/tmp/little_things/table2.dat');
const ltPoints = [];
for (const nm of Object.keys(rotT)) {
  if (!rotD[nm] || !t2[nm]) continue;
  const cT = rotT[nm].sort((a, b) => a.r - b.r);
  const cD = rotD[nm].sort((a, b) => a.r - b.r);
  if (cT.length < 5 || t2[nm].mbar <= 0) continue;
  for (const pt of cT) {
    const r = pt.r, vO = pt.v;
    const vDM = interpV(cD, r);
    if (isNaN(vDM) || vO <= 0 || r <= 0) continue;
    let vBsq = vO * vO - vDM * vDM;
    if (vBsq < 0) vBsq = 0;
    const gObs = vO * vO / r, gBar = vBsq / r;
    if (gBar <= 0 || gObs <= 0) continue;
    const sigBar = t2[nm].mbar / (Math.PI * r * r);
    if (!isFinite(sigBar) || sigBar <= 0) continue;
    ltPoints.push({ galaxy: nm, gObs, gBar, sigBar, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), logSig: Math.log10(sigBar) });
  }
}

const allPoints = [...sparcPoints, ...ltPoints];
console.log("SPARC: " + sparcPoints.length + " points, LITTLE THINGS: " + ltPoints.length + " points, TOTAL: " + allPoints.length + "\n");

function fitAndEval(points, buildX) {
  const y = points.map(p => p.logGobs);
  const X = points.map(p => buildX(p));
  const beta = solveLS(X, y);
  const resids = points.map((p, i) => {
    let pred = 0;
    const row = buildX(p);
    for (let j = 0; j < beta.length; j++) pred += beta[j] * row[j];
    return p.logGobs - pred;
  });
  return { beta, rms: rmsScatter(resids), mad: medAbsDev(resids), n: points.length };
}

function evalWithBeta(points, buildX, beta) {
  const resids = points.map(p => {
    let pred = 0;
    const row = buildX(p);
    for (let j = 0; j < beta.length; j++) pred += beta[j] * row[j];
    return p.logGobs - pred;
  });
  return { rms: rmsScatter(resids), mad: medAbsDev(resids), n: points.length, mean: resids.reduce((a, b) => a + b, 0) / resids.length };
}

console.log("============================================");
console.log("CONTROLLED COMPARISON: g_bar ALONE vs g_bar + Sigma_bar");
console.log("============================================\n");

const modelDefs = [
  { name: "M0: RAR residuals only", buildX: p => [1] },
  { name: "M1: log(g_bar) linear", buildX: p => [1, p.logGbar], params: 2 },
  { name: "M2: log(g_bar) quadratic", buildX: p => [1, p.logGbar, p.logGbar * p.logGbar], params: 3 },
  { name: "M3: log(g_bar) cubic", buildX: p => [1, p.logGbar, p.logGbar * p.logGbar, p.logGbar ** 3], params: 4 },
  { name: "M4: log(g_bar) quad + log(Sig)", buildX: p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig], params: 4 },
  { name: "M5: FULL 5-param (quad + Sig + interaction)", buildX: p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig, p.logGbar * p.logSig], params: 5 },
  { name: "M6: log(Sig) only (control)", buildX: p => [1, p.logSig], params: 2 },
  { name: "M7: log(g_bar) + log(Sig) linear", buildX: p => [1, p.logGbar, p.logSig], params: 3 },
];

console.log("--- Fit on ALL data ---\n");
console.log("  " + "Model".padEnd(50) + "RMS".padStart(10) + "  MAD".padStart(10) + "  #params".padStart(8));

const allFits = [];
for (const m of modelDefs) {
  const res = fitAndEval(allPoints, m.buildX);
  console.log("  " + m.name.padEnd(50) + res.rms.toFixed(6).padStart(10) + ("  " + res.mad.toFixed(6)).padStart(10) + ("  " + res.beta.length).padStart(8));
  allFits.push({ ...m, ...res });
}

const baseRMS = evalWithBeta(allPoints, p => [1], [allPoints.reduce((a, p) => a + p.logGobs, 0) / allPoints.length]).rms;

console.log("\n--- KEY COMPARISON ---");
const m2 = allFits.find(f => f.name.includes("M2"));
const m4 = allFits.find(f => f.name.includes("M4"));
const m5 = allFits.find(f => f.name.includes("M5"));
const m3 = allFits.find(f => f.name.includes("M3"));

console.log("\n  Adding Sigma_bar to quadratic g_bar model:");
console.log("    M2 (g_bar quad only): rms=" + m2.rms.toFixed(6));
console.log("    M4 (+Sigma_bar):      rms=" + m4.rms.toFixed(6) + " (delta=" + ((m2.rms - m4.rms) / m2.rms * 100).toFixed(3) + "%)");
console.log("    M5 (+Sig+interaction): rms=" + m5.rms.toFixed(6) + " (delta=" + ((m2.rms - m5.rms) / m2.rms * 100).toFixed(3) + "%)");
console.log("    M3 (g_bar CUBIC, no Sig): rms=" + m3.rms.toFixed(6) + " (delta=" + ((m2.rms - m3.rms) / m2.rms * 100).toFixed(3) + "%)");

const sigmaContrib = ((m2.rms - m4.rms) / m2.rms * 100);
const totalContrib = ((m2.rms - m5.rms) / m2.rms * 100);
const cubicContrib = ((m2.rms - m3.rms) / m2.rms * 100);

console.log("\n  Verdict: Sigma_bar contributes " + sigmaContrib.toFixed(2) + "% improvement");
console.log("  Comparison: adding a CUBIC g_bar term contributes " + cubicContrib.toFixed(2) + "% improvement");
console.log("  -> " + (sigmaContrib > cubicContrib ? "Sigma_bar contributes MORE than cubic g_bar term" : "Adding another g_bar power contributes SAME or MORE than Sigma_bar"));

console.log("\n============================================");
console.log("F-TEST: Is Sigma_bar statistically significant?");
console.log("============================================\n");

const n = allPoints.length;
const ssRes2 = m2.rms ** 2 * n;
const ssRes4 = m4.rms ** 2 * n;
const ssRes5 = m5.rms ** 2 * n;
const ssRes3 = m3.rms ** 2 * n;

const fStat24 = ((ssRes2 - ssRes4) / 1) / (ssRes4 / (n - 4));
const fStat25 = ((ssRes2 - ssRes5) / 2) / (ssRes5 / (n - 5));
const fStat45 = ((ssRes4 - ssRes5) / 1) / (ssRes5 / (n - 5));

console.log("  F-test M2 vs M4 (adding log(Sig)):        F=" + fStat24.toFixed(1) + (fStat24 > 6.63 ? " (p<0.01 SIGNIFICANT)" : " (NOT significant)"));
console.log("  F-test M2 vs M5 (adding Sig+interaction): F=" + fStat25.toFixed(1) + (fStat25 > 6.63 ? " (p<0.01 SIGNIFICANT)" : " (NOT significant)"));
console.log("  F-test M4 vs M5 (adding interaction):     F=" + fStat45.toFixed(1) + (fStat45 > 6.63 ? " (p<0.01 SIGNIFICANT)" : " (NOT significant)"));

console.log("\n============================================");
console.log("TRUE CROSS-VALIDATION: Refit on each dataset");
console.log("============================================\n");

function trueCV(trainPts, testPts, buildX, modelName) {
  const trainFit = fitAndEval(trainPts, buildX);
  const testEval = evalWithBeta(testPts, buildX, trainFit.beta);
  return { trainRms: trainFit.rms, testRms: testEval.rms, testMean: testEval.mean, beta: trainFit.beta };
}

const models = [
  { name: "M1: log(g_bar) linear", buildX: p => [1, p.logGbar] },
  { name: "M2: log(g_bar) quadratic", buildX: p => [1, p.logGbar, p.logGbar * p.logGbar] },
  { name: "M4: quad + log(Sig)", buildX: p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig] },
  { name: "M5: FULL bivariate", buildX: p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig, p.logGbar * p.logSig] },
  { name: "M3: cubic g_bar only", buildX: p => [1, p.logGbar, p.logGbar * p.logGbar, p.logGbar ** 3] },
];

console.log("  Direction: Train on SPARC, test on LITTLE THINGS\n");
console.log("  " + "Model".padEnd(40) + "Train RMS".padStart(10) + "  Test RMS".padStart(10) + "  Test Mean".padStart(10));
for (const m of models) {
  const cv = trueCV(sparcPoints, ltPoints, m.buildX, m.name);
  console.log("  " + m.name.padEnd(40) + cv.trainRms.toFixed(6).padStart(10) + ("  " + cv.testRms.toFixed(6)).padStart(10) + ("  " + cv.testMean.toFixed(5)).padStart(10));
}

console.log("\n  Direction: Train on LITTLE THINGS, test on SPARC\n");
console.log("  " + "Model".padEnd(40) + "Train RMS".padStart(10) + "  Test RMS".padStart(10) + "  Test Mean".padStart(10));
for (const m of models) {
  const cv = trueCV(ltPoints, sparcPoints, m.buildX, m.name);
  console.log("  " + m.name.padEnd(40) + cv.trainRms.toFixed(6).padStart(10) + ("  " + cv.testRms.toFixed(6)).padStart(10) + ("  " + cv.testMean.toFixed(5)).padStart(10));
}

console.log("\n============================================");
console.log("THE DECISIVE COMPARISON");
console.log("============================================\n");

const cvM2_SL = trueCV(sparcPoints, ltPoints, p => [1, p.logGbar, p.logGbar * p.logGbar], "M2");
const cvM4_SL = trueCV(sparcPoints, ltPoints, p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig], "M4");
const cvM5_SL = trueCV(sparcPoints, ltPoints, p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig, p.logGbar * p.logSig], "M5");
const cvM3_SL = trueCV(sparcPoints, ltPoints, p => [1, p.logGbar, p.logGbar * p.logGbar, p.logGbar ** 3], "M3");

const cvM2_LS = trueCV(ltPoints, sparcPoints, p => [1, p.logGbar, p.logGbar * p.logGbar], "M2");
const cvM4_LS = trueCV(ltPoints, sparcPoints, p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig], "M4");
const cvM5_LS = trueCV(ltPoints, sparcPoints, p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig, p.logGbar * p.logSig], "M5");
const cvM3_LS = trueCV(ltPoints, sparcPoints, p => [1, p.logGbar, p.logGbar * p.logGbar, p.logGbar ** 3], "M3");

console.log("  SPARC -> LT cross-validation (test RMS):");
console.log("    M2 (g_bar quad only):      " + cvM2_SL.testRms.toFixed(6));
console.log("    M4 (quad + Sigma):         " + cvM4_SL.testRms.toFixed(6) + " delta=" + ((cvM2_SL.testRms - cvM4_SL.testRms) / cvM2_SL.testRms * 100).toFixed(2) + "%");
console.log("    M5 (full bivariate):       " + cvM5_SL.testRms.toFixed(6) + " delta=" + ((cvM2_SL.testRms - cvM5_SL.testRms) / cvM2_SL.testRms * 100).toFixed(2) + "%");
console.log("    M3 (cubic, no Sigma):      " + cvM3_SL.testRms.toFixed(6) + " delta=" + ((cvM2_SL.testRms - cvM3_SL.testRms) / cvM2_SL.testRms * 100).toFixed(2) + "%");

console.log("\n  LT -> SPARC cross-validation (test RMS):");
console.log("    M2 (g_bar quad only):      " + cvM2_LS.testRms.toFixed(6));
console.log("    M4 (quad + Sigma):         " + cvM4_LS.testRms.toFixed(6) + " delta=" + ((cvM2_LS.testRms - cvM4_LS.testRms) / cvM2_LS.testRms * 100).toFixed(2) + "%");
console.log("    M5 (full bivariate):       " + cvM5_LS.testRms.toFixed(6) + " delta=" + ((cvM2_LS.testRms - cvM5_LS.testRms) / cvM2_LS.testRms * 100).toFixed(2) + "%");
console.log("    M3 (cubic, no Sigma):      " + cvM3_LS.testRms.toFixed(6) + " delta=" + ((cvM2_LS.testRms - cvM3_LS.testRms) / cvM2_LS.testRms * 100).toFixed(2) + "%");

const sigmaHelpsSL = cvM4_SL.testRms < cvM2_SL.testRms;
const sigmaHelpsLS = cvM4_LS.testRms < cvM2_LS.testRms;
const sigmaBetterThanCubicSL = cvM4_SL.testRms < cvM3_SL.testRms;
const sigmaBetterThanCubicLS = cvM4_LS.testRms < cvM3_LS.testRms;

console.log("\n  Does Sigma_bar improve cross-dataset predictions?");
console.log("    SPARC->LT: " + (sigmaHelpsSL ? "YES" : "NO") + " (M4 vs M2: " + ((cvM2_SL.testRms - cvM4_SL.testRms) / cvM2_SL.testRms * 100).toFixed(2) + "%)");
console.log("    LT->SPARC: " + (sigmaHelpsLS ? "YES" : "NO") + " (M4 vs M2: " + ((cvM2_LS.testRms - cvM4_LS.testRms) / cvM2_LS.testRms * 100).toFixed(2) + "%)");
console.log("    Better than cubic (SL): " + (sigmaBetterThanCubicSL ? "YES" : "NO"));
console.log("    Better than cubic (LS): " + (sigmaBetterThanCubicLS ? "YES" : "NO"));

console.log("\n============================================");
console.log("COEFFICIENT STABILITY CHECK");
console.log("============================================\n");

const allBeta = fitAndEval(allPoints, p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig]).beta;
const sparcBeta = fitAndEval(sparcPoints, p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig]).beta;
const ltBeta = fitAndEval(ltPoints, p => [1, p.logGbar, p.logGbar * p.logGbar, p.logSig]).beta;

console.log("  M4 coefficients: log(g_obs) = a + b*log(g_bar) + c*log(g_bar)^2 + d*log(Sig)\n");
console.log("  " + "".padEnd(12) + "a".padStart(10) + "b".padStart(10) + "c".padStart(10) + "d(Sig)".padStart(10));
console.log("  " + "ALL".padEnd(12) + allBeta[0].toFixed(4).padStart(10) + allBeta[1].toFixed(4).padStart(10) + allBeta[2].toFixed(6).padStart(10) + allBeta[3].toFixed(6).padStart(10));
console.log("  " + "SPARC only".padEnd(12) + sparcBeta[0].toFixed(4).padStart(10) + sparcBeta[1].toFixed(4).padStart(10) + sparcBeta[2].toFixed(6).padStart(10) + sparcBeta[3].toFixed(6).padStart(10));
console.log("  " + "LT only".padEnd(12) + ltBeta[0].toFixed(4).padStart(10) + ltBeta[1].toFixed(4).padStart(10) + ltBeta[2].toFixed(6).padStart(10) + ltBeta[3].toFixed(6).padStart(10));

const sigCoeffStable = Math.sign(sparcBeta[3]) === Math.sign(ltBeta[3]);
console.log("\n  Sigma coefficient sign: SPARC=" + (sparcBeta[3] > 0 ? "+" : "-") + " LT=" + (ltBeta[3] > 0 ? "+" : "-") + " " + (sigCoeffStable ? "CONSISTENT" : "INCONSISTENT"));
console.log("  Sigma coefficient magnitude ratio: " + Math.abs(sparcBeta[3] / ltBeta[3]).toFixed(2));

console.log("\n============================================");
console.log("INFORMATION CRITERIA (accounting for model complexity)");
console.log("============================================\n");

function aic(n, k, ssRes) { return n * Math.log(ssRes / n) + 2 * k; }
function bic(n, k, ssRes) { return n * Math.log(ssRes / n) + k * Math.log(n); }

const ssM1 = allFits.find(f => f.name.includes("M1")).rms ** 2 * n;
const ssM2_ = m2.rms ** 2 * n;
const ssM3_ = m3.rms ** 2 * n;
const ssM4_ = m4.rms ** 2 * n;
const ssM5_ = m5.rms ** 2 * n;

console.log("  " + "Model".padEnd(35) + "AIC".padStart(10) + "  BIC".padStart(10) + "  #params".padStart(8));
console.log("  " + "M1 (linear g_bar)".padEnd(35) + aic(n, 2, ssM1).toFixed(1).padStart(10) + ("  " + bic(n, 2, ssM1).toFixed(1)).padStart(10) + "  2".padStart(8));
console.log("  " + "M2 (quad g_bar)".padEnd(35) + aic(n, 3, ssM2_).toFixed(1).padStart(10) + ("  " + bic(n, 3, ssM2_).toFixed(1)).padStart(10) + "  3".padStart(8));
console.log("  " + "M3 (cubic g_bar, no Sig)".padEnd(35) + aic(n, 4, ssM3_).toFixed(1).padStart(10) + ("  " + bic(n, 4, ssM3_).toFixed(1)).padStart(10) + "  4".padStart(8));
console.log("  " + "M4 (quad + Sig)".padEnd(35) + aic(n, 4, ssM4_).toFixed(1).padStart(10) + ("  " + bic(n, 4, ssM4_).toFixed(1)).padStart(10) + "  4".padStart(8));
console.log("  " + "M5 (full bivariate)".padEnd(35) + aic(n, 5, ssM5_).toFixed(1).padStart(10) + ("  " + bic(n, 5, ssM5_).toFixed(1)).padStart(10) + "  5".padStart(8));

console.log("\n  delta_AIC(M2->M4) = " + (aic(n, 3, ssM2_) - aic(n, 4, ssM4_)).toFixed(1) + " (positive = M4 wins)");
console.log("  delta_BIC(M2->M4) = " + (bic(n, 3, ssM2_) - bic(n, 4, ssM4_)).toFixed(1) + " (positive = M4 wins)");
console.log("  delta_AIC(M3->M4) = " + (aic(n, 4, ssM3_) - aic(n, 4, ssM4_)).toFixed(1) + " (positive = M4 wins, same #params)");

console.log("\n============================================");
console.log("FINAL HONEST VERDICT");
console.log("============================================\n");

const sigmaNetCV = ((cvM2_SL.testRms - cvM4_SL.testRms) + (cvM2_LS.testRms - cvM4_LS.testRms)) / 2;
const sigmaPercentCV = ((cvM2_SL.testRms - cvM4_SL.testRms) / cvM2_SL.testRms * 100 + (cvM2_LS.testRms - cvM4_LS.testRms) / cvM2_LS.testRms * 100) / 2;

console.log("  1. The 30% scatter reduction was MOSTLY from better g_bar fitting (polynomial vs McGaugh)");
console.log("  2. Sigma_bar's MARGINAL contribution (controlling for g_bar polynomial):");
console.log("     - In-sample (F-test): F=" + fStat24.toFixed(1) + (fStat24 > 6.63 ? " SIGNIFICANT" : " not significant"));
console.log("     - Cross-validated average: " + sigmaPercentCV.toFixed(2) + "% improvement");
console.log("     - Sign consistent: " + (sigCoeffStable ? "YES" : "NO"));
console.log("  3. Is the marginal Sigma contribution better than just adding g_bar^3?");
console.log("     - SPARC->LT: Sigma " + (sigmaBetterThanCubicSL ? ">" : "<=") + " cubic");
console.log("     - LT->SPARC: Sigma " + (sigmaBetterThanCubicLS ? ">" : "<=") + " cubic");

const genuineSigma = sigmaHelpsSL && sigmaHelpsLS && fStat24 > 6.63 && sigCoeffStable;
if (genuineSigma && sigmaPercentCV > 5) {
  console.log("\n  >>> TREASURE: Sigma_bar adds genuine, cross-validated predictive power beyond g_bar optimization <<<");
} else if (genuineSigma && sigmaPercentCV > 1) {
  console.log("\n  >>> REAL BUT SMALL: Sigma_bar adds genuine predictive power, but the effect is modest <<<");
} else if (genuineSigma) {
  console.log("\n  >>> DETECTABLE BUT MARGINAL: Sigma_bar contribution is statistically real but practically tiny <<<");
} else {
  console.log("\n  >>> NO INDEPENDENT SIGMA EFFECT: The improvement comes from better g_bar fitting, not from Sigma_bar <<<");
}

const output = {
  timestamp: new Date().toISOString(),
  nPoints: { sparc: sparcPoints.length, lt: ltPoints.length, all: allPoints.length },
  inSampleComparison: {
    M2_gbar_quad: { rms: m2.rms, params: 3 },
    M3_gbar_cubic: { rms: m3.rms, params: 4 },
    M4_quad_plus_sigma: { rms: m4.rms, params: 4 },
    M5_full_bivariate: { rms: m5.rms, params: 5 },
    fTest_M2_vs_M4: fStat24,
    sigmaContribution: sigmaContrib
  },
  crossValidation: {
    sparc_to_lt: {
      M2: cvM2_SL.testRms, M3: cvM3_SL.testRms, M4: cvM4_SL.testRms, M5: cvM5_SL.testRms,
      sigmaImprovement: ((cvM2_SL.testRms - cvM4_SL.testRms) / cvM2_SL.testRms * 100)
    },
    lt_to_sparc: {
      M2: cvM2_LS.testRms, M3: cvM3_LS.testRms, M4: cvM4_LS.testRms, M5: cvM5_LS.testRms,
      sigmaImprovement: ((cvM2_LS.testRms - cvM4_LS.testRms) / cvM2_LS.testRms * 100)
    }
  },
  coefficientStability: {
    sparc: { sigmaCoeff: sparcBeta[3] },
    lt: { sigmaCoeff: ltBeta[3] },
    signConsistent: sigCoeffStable
  },
  verdict: {
    genuineSigmaEffect: genuineSigma,
    crossValidatedImprovement: sigmaPercentCV,
    conclusion: genuineSigma && sigmaPercentCV > 5 ? "TREASURE" : genuineSigma && sigmaPercentCV > 1 ? "REAL_BUT_SMALL" : genuineSigma ? "DETECTABLE_BUT_MARGINAL" : "NO_INDEPENDENT_SIGMA"
  }
};

const outPath = path.join(__dirname, '..', 'public', 'formula-search.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log("\nResults saved to: " + outPath);
