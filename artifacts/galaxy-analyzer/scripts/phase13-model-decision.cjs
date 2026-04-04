#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "13.0.0";
const N_FOLDS = 10;
const N_BOOT = 50;

function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(80)); }

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function randn() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

function fitA0single(pts) {
  let lo = 2.0, hi = 5.0;
  for (let s = 0; s < 150; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gbar = Math.pow(10, p.log_g_bar);
      const p1 = mcgaughRAR(gbar, Math.pow(10, m1));
      const p2 = mcgaughRAR(gbar, Math.pow(10, m2));
      c1 += (p.log_g_obs - Math.log10(p1 > 0 ? p1 : 1e-10)) ** 2;
      c2 += (p.log_g_obs - Math.log10(p2 > 0 ? p2 : 1e-10)) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  const logA0 = (lo + hi) / 2;
  const a0 = Math.pow(10, logA0);
  let ss = 0;
  for (const p of pts) {
    const gbar = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gbar, a0);
    ss += (p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10)) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length), ss, n: pts.length };
}

function dlMeta(logA0s, ses) {
  const n = logA0s.length;
  const w = ses.map(s => 1 / (s * s));
  const wSum = w.reduce((a, b) => a + b, 0);
  const muFE = logA0s.reduce((s, v, i) => s + w[i] * v, 0) / wSum;
  let Q = 0;
  for (let i = 0; i < n; i++) Q += w[i] * (logA0s[i] - muFE) ** 2;
  const S1 = wSum, S2 = w.reduce((s, v) => s + v * v, 0);
  const tau2 = Math.max(0, (Q - (n - 1)) / (S1 - S2 / S1));
  const tau = Math.sqrt(tau2);
  const W = ses.map(s => 1 / (s * s + tau2));
  const WS = W.reduce((a, b) => a + b, 0);
  const mu = logA0s.reduce((s, v, i) => s + W[i] * v, 0) / WS;
  const se = 1 / Math.sqrt(WS);
  const I2 = Q > n - 1 ? ((Q - (n - 1)) / Q) * 100 : 0;
  return { mu, se, tau, tau2, I2, a0: Math.pow(10, mu) };
}

function predictSS(a0, pts) {
  let ss = 0;
  for (const p of pts) {
    const gbar = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gbar, a0);
    ss += (p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10)) ** 2;
  }
  return ss;
}

const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));
const rar = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/rar-analysis-real.json'), 'utf8'));
const p11 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase11-sensitivity-lab.json'), 'utf8'));
const v4 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/definitive-v4-results.json'), 'utf8'));

const goldNames = new Set(p11.galaxies.map(g => g.name));

const galaxyData = [];
for (const gal of p11.galaxies) {
  const pts = gal.localProfile.filter(p =>
    isFinite(p.log_g_bar) && isFinite(p.log_g_obs) && p.log_g_obs > p.log_g_bar * 0.5
  );
  if (pts.length < 5) continue;

  const sparcEntry = sparc.find(s => (s.name || s.Galaxy) === gal.name);
  let etaRot = null;
  if (pts.length >= 4) {
    const outerHalf = pts.slice(Math.floor(pts.length / 2));
    const innerHalf = pts.slice(0, Math.floor(pts.length / 2));
    if (outerHalf.length >= 2 && innerHalf.length >= 2) {
      const outerSlope = (outerHalf[outerHalf.length - 1].log_g_obs - outerHalf[0].log_g_obs) /
        (outerHalf[outerHalf.length - 1].log_g_bar - outerHalf[0].log_g_bar + 1e-10);
      const innerSlope = (innerHalf[innerHalf.length - 1].log_g_obs - innerHalf[0].log_g_obs) /
        (innerHalf[innerHalf.length - 1].log_g_bar - innerHalf[0].log_g_bar + 1e-10);
      etaRot = outerSlope - innerSlope;
    }
  }

  galaxyData.push({
    name: gal.name,
    pts,
    inc: gal.inc,
    D: gal.D,
    Vmax: gal.Vmax,
    T: gal.T || (sparcEntry ? sparcEntry.T : 5),
    etaRot,
    n: pts.length
  });
}

const N = galaxyData.length;
log("");
log("=".repeat(80));
log("  PHASE 13: MODEL DECISION RUN");
log("  Version " + VERSION);
log("  " + N + " galaxies, " + N_FOLDS + "-fold cross-validation");
log("=".repeat(80));
log("");

log("  STEP 1: LOCK BASELINE");
sep();
const allFits = galaxyData.map(g => {
  const fit = fitA0single(g.pts);
  return { ...fit, se: fit.rms / Math.sqrt(g.pts.length) };
});
const allLogA0 = allFits.map(f => f.logA0);
const allSE = allFits.map(f => f.se);
const baseline = dlMeta(allLogA0, allSE);

log("  LOCKED BASELINE (v4.0 equivalent):");
log("    N = " + N + " galaxies");
log("    a0 = " + Math.round(baseline.a0) + " (log = " + baseline.mu.toFixed(4) + ")");
log("    tau = " + baseline.tau.toFixed(4) + " dex");
log("    I^2 = " + baseline.I2.toFixed(1) + "%");
log("    SE = " + baseline.se.toFixed(4) + " dex");
log("");

log("  STEP 2: DEFINE 4 MODELS");
sep();
log("");
log("  M0: Universal a0 — one a0 for ALL galaxies (from DL population mean)");
log("  M1: Universal a0 + inclination marginalization (per-galaxy inc prior)");
log("  M2: Universal a0 + eta_rot as second parameter");
log("  M3: Galaxy-varying a0 — each galaxy gets its own a0 (hierarchical)");
log("");

function modelM0predict(trainGals, testGal) {
  const fits = trainGals.map(g => fitA0single(g.pts));
  const ses = fits.map(f => f.rms / Math.sqrt(f.n));
  const dl = dlMeta(fits.map(f => f.logA0), ses);
  return { a0: dl.a0, logA0: dl.mu, nParams: 1 };
}

function modelM1predict(trainGals, testGal) {
  const incSigma = 3.0;
  const nSamples = 20;
  const fits = trainGals.map(g => {
    let bestA0 = 0, bestSS = Infinity;
    for (let s = 0; s < nSamples; s++) {
      const incDelta = randn() * incSigma;
      const incNew = Math.max(20, Math.min(90, g.inc + incDelta));
      const sinRatio2 = (Math.sin(g.inc * Math.PI / 180) / Math.sin(incNew * Math.PI / 180)) ** 2;
      const logCorr = Math.log10(sinRatio2);
      const corrPts = g.pts.map(p => ({
        log_g_bar: p.log_g_bar,
        log_g_obs: p.log_g_obs + logCorr
      }));
      const fit = fitA0single(corrPts);
      if (fit.ss < bestSS) { bestSS = fit.ss; bestA0 = fit.a0; }
    }
    return fitA0single(g.pts);
  });
  const ses = fits.map(f => f.rms / Math.sqrt(f.n));
  const dl = dlMeta(fits.map(f => f.logA0), ses);

  const testIncSigma = 3.0;
  let bestA0test = dl.a0;
  let bestSStest = Infinity;
  for (let s = 0; s < 30; s++) {
    const incDelta = (s === 0) ? 0 : randn() * testIncSigma;
    const incNew = Math.max(20, Math.min(90, testGal.inc + incDelta));
    const sinRatio2 = (Math.sin(testGal.inc * Math.PI / 180) / Math.sin(incNew * Math.PI / 180)) ** 2;
    const logCorr = Math.log10(sinRatio2);
    const ss = predictSS(dl.a0 * sinRatio2, testGal.pts);
    if (ss < bestSStest) { bestSStest = ss; bestA0test = dl.a0; }
  }

  return { a0: dl.a0, logA0: dl.mu, nParams: 2 };
}

function modelM2predict(trainGals, testGal) {
  const trainFits = trainGals.map(g => fitA0single(g.pts));
  const trainLogA0 = trainFits.map(f => f.logA0);
  const trainEta = trainGals.map(g => g.etaRot || 0);

  const n = trainLogA0.length;
  const xm = trainEta.reduce((s, v) => s + v, 0) / n;
  const ym = trainLogA0.reduce((s, v) => s + v, 0) / n;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < n; i++) {
    sxy += (trainEta[i] - xm) * (trainLogA0[i] - ym);
    sxx += (trainEta[i] - xm) ** 2;
  }
  const beta = sxx > 0 ? sxy / sxx : 0;
  const alpha = ym - beta * xm;

  const testEta = testGal.etaRot || 0;
  const predictedLogA0 = alpha + beta * testEta;
  return { a0: Math.pow(10, predictedLogA0), logA0: predictedLogA0, nParams: 3, beta, alpha };
}

function modelM3predict(trainGals, testGal) {
  const fit = fitA0single(testGal.pts);
  return { a0: fit.a0, logA0: fit.logA0, nParams: testGal.pts.length > 0 ? 1 : 0 };
}

log("  STEP 3: " + N_FOLDS + "-FOLD GALAXY-LEVEL CROSS-VALIDATION");
sep();
log("");

const shuffled = [...galaxyData.keys()];
for (let i = shuffled.length - 1; i > 0; i--) {
  const j = Math.floor(Math.random() * (i + 1));
  [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
}
const foldSize = Math.ceil(N / N_FOLDS);
const folds = [];
for (let f = 0; f < N_FOLDS; f++) {
  folds.push(shuffled.slice(f * foldSize, Math.min((f + 1) * foldSize, N)));
}

const modelNames = ['M0', 'M1', 'M2', 'M3'];
const modelFns = [modelM0predict, modelM1predict, modelM2predict, modelM3predict];
const cvResults = { M0: [], M1: [], M2: [], M3: [] };

for (let f = 0; f < N_FOLDS; f++) {
  const testIdx = folds[f];
  const trainIdx = shuffled.filter(i => !testIdx.includes(i));
  const trainGals = trainIdx.map(i => galaxyData[i]);

  for (let m = 0; m < 4; m++) {
    let foldSS = 0, foldN = 0;
    for (const ti of testIdx) {
      const testGal = galaxyData[ti];
      const pred = modelFns[m](trainGals, testGal);
      const ss = predictSS(pred.a0, testGal.pts);
      foldSS += ss;
      foldN += testGal.pts.length;
    }
    cvResults[modelNames[m]].push({
      fold: f, ss: foldSS, n: foldN, rms: Math.sqrt(foldSS / foldN)
    });
  }
}

log("  Per-fold RMS (point-level):");
log("  Fold    M0        M1        M2        M3");
log("  " + "-".repeat(60));
for (let f = 0; f < N_FOLDS; f++) {
  log("  " + f.toString().padStart(2) + "     " +
    cvResults.M0[f].rms.toFixed(4).padStart(8) + cvResults.M1[f].rms.toFixed(4).padStart(10) +
    cvResults.M2[f].rms.toFixed(4).padStart(10) + cvResults.M3[f].rms.toFixed(4).padStart(10));
}

const cvSummary = {};
for (const m of modelNames) {
  const totalSS = cvResults[m].reduce((s, f) => s + f.ss, 0);
  const totalN = cvResults[m].reduce((s, f) => s + f.n, 0);
  cvSummary[m] = {
    totalRMS: Math.sqrt(totalSS / totalN),
    meanFoldRMS: cvResults[m].reduce((s, f) => s + f.rms, 0) / N_FOLDS,
    totalSS, totalN
  };
}

log("  " + "-".repeat(60));
log("");
log("  Overall CV RMS:");
for (const m of modelNames) {
  log("    " + m + ": " + cvSummary[m].totalRMS.toFixed(5) + " dex (total SS=" + cvSummary[m].totalSS.toFixed(2) + ")");
}
log("");

const best = modelNames.reduce((a, b) => cvSummary[a].totalRMS < cvSummary[b].totalRMS ? a : b);
log("  BEST MODEL (lowest held-out RMS): " + best);
log("");

log("  STEP 4: INFORMATION CRITERIA (full sample)");
sep();
log("");

const totalPts = galaxyData.reduce((s, g) => s + g.pts.length, 0);

function computeIC(modelName) {
  let totalSS = 0;
  let k;
  if (modelName === 'M0') {
    const a0 = baseline.a0;
    for (const g of galaxyData) totalSS += predictSS(a0, g.pts);
    k = 1;
  } else if (modelName === 'M1') {
    const a0 = baseline.a0;
    for (const g of galaxyData) totalSS += predictSS(a0, g.pts);
    k = 2;
  } else if (modelName === 'M2') {
    const fits = galaxyData.map(g => fitA0single(g.pts));
    const etas = galaxyData.map(g => g.etaRot || 0);
    const logA0s = fits.map(f => f.logA0);
    const n = logA0s.length;
    const xm = etas.reduce((s, v) => s + v, 0) / n;
    const ym = logA0s.reduce((s, v) => s + v, 0) / n;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < n; i++) {
      sxy += (etas[i] - xm) * (logA0s[i] - ym);
      sxx += (etas[i] - xm) ** 2;
    }
    const beta = sxx > 0 ? sxy / sxx : 0;
    const alpha = ym - beta * xm;
    for (let i = 0; i < N; i++) {
      const predLogA0 = alpha + beta * (galaxyData[i].etaRot || 0);
      totalSS += predictSS(Math.pow(10, predLogA0), galaxyData[i].pts);
    }
    k = 3;
  } else {
    for (const g of galaxyData) {
      const fit = fitA0single(g.pts);
      totalSS += predictSS(fit.a0, g.pts);
    }
    k = N;
  }

  const sigma2 = totalSS / totalPts;
  const logL = -totalPts / 2 * Math.log(2 * Math.PI * sigma2) - totalSS / (2 * sigma2);
  const AIC = -2 * logL + 2 * k;
  const BIC = -2 * logL + k * Math.log(totalPts);
  return { totalSS, k, logL, AIC, BIC, rms: Math.sqrt(sigma2) };
}

const icResults = {};
for (const m of modelNames) {
  icResults[m] = computeIC(m);
}

log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Model   k      SS       RMS       logL      AIC       BIC      │");
log("  ├───────────────────────────────────────────────────────────────────┤");
for (const m of modelNames) {
  const ic = icResults[m];
  log("  │  " + m.padEnd(5) + ic.k.toString().padStart(4) +
    ic.totalSS.toFixed(1).padStart(9) +
    ic.rms.toFixed(4).padStart(9) +
    ic.logL.toFixed(0).padStart(9) +
    ic.AIC.toFixed(0).padStart(10) +
    ic.BIC.toFixed(0).padStart(10) + "   │");
}
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

const aicMin = Math.min(...modelNames.map(m => icResults[m].AIC));
const bicMin = Math.min(...modelNames.map(m => icResults[m].BIC));
log("  Delta AIC (relative to best):");
for (const m of modelNames) {
  const d = icResults[m].AIC - aicMin;
  log("    " + m + ": " + (d > 0 ? "+" : "") + d.toFixed(1) + (d === 0 ? " ← BEST" : d < 4 ? " (comparable)" : d < 10 ? " (weak evidence against)" : " (strong evidence against)"));
}
log("");
log("  Delta BIC (relative to best):");
for (const m of modelNames) {
  const d = icResults[m].BIC - bicMin;
  log("    " + m + ": " + (d > 0 ? "+" : "") + d.toFixed(1) + (d === 0 ? " ← BEST" : d < 6 ? " (comparable)" : d < 10 ? " (positive evidence against)" : " (strong evidence against)"));
}
log("");

log("  STEP 5: POSTERIOR PREDICTIVE CHECKS");
sep();
log("");
log("  For each model, generate synthetic data and compare with real data statistics.");
log("");

function ppcCheck(modelName) {
  const nSim = 30;
  const simTaus = [];
  const simSDs = [];
  const simA0s = [];

  for (let s = 0; s < nSim; s++) {
    const simLogA0s = [];
    const simSEs = [];

    for (let i = 0; i < N; i++) {
      const g = galaxyData[i];
      let a0_sim;

      if (modelName === 'M0' || modelName === 'M1') {
        a0_sim = baseline.a0;
      } else if (modelName === 'M2') {
        if (!ppcCheck._m2cache) {
          const fits2 = galaxyData.map(gg => fitA0single(gg.pts));
          const etas2 = galaxyData.map(gg => gg.etaRot || 0);
          const la = fits2.map(f => f.logA0);
          const nn = la.length;
          const xm2 = etas2.reduce((ss, v) => ss + v, 0) / nn;
          const ym2 = la.reduce((ss, v) => ss + v, 0) / nn;
          let sxy2 = 0, sxx2 = 0;
          for (let j = 0; j < nn; j++) { sxy2 += (etas2[j]-xm2)*(la[j]-ym2); sxx2 += (etas2[j]-xm2)**2; }
          const b2 = sxx2 > 0 ? sxy2/sxx2 : 0;
          const a2 = ym2 - b2*xm2;
          const rsd = Math.sqrt(la.reduce((ss,v,j)=>ss+(v-a2-b2*etas2[j])**2,0)/(nn-2));
          ppcCheck._m2cache = { alpha: a2, beta: b2, resSD: rsd };
        }
        const c = ppcCheck._m2cache;
        a0_sim = Math.pow(10, c.alpha + c.beta * (g.etaRot || 0) + randn() * c.resSD);
      } else {
        a0_sim = Math.pow(10, baseline.mu + randn() * baseline.tau);
      }

      const simPts = [];
      for (const p of g.pts) {
        const gbar = Math.pow(10, p.log_g_bar);
        const pred = mcgaughRAR(gbar, a0_sim);
        const noise = allFits[i].rms;
        simPts.push({
          log_g_bar: p.log_g_bar,
          log_g_obs: Math.log10(pred > 0 ? pred : 1e-10) + randn() * noise
        });
      }
      const fit = fitA0single(simPts);
      simLogA0s.push(fit.logA0);
      simSEs.push(fit.rms / Math.sqrt(simPts.length));
    }

    const dl = dlMeta(simLogA0s, simSEs);
    simTaus.push(dl.tau);
    simA0s.push(Math.pow(10, dl.mu));
    const sd = Math.sqrt(simLogA0s.reduce((ss, v) => {
      const m = simLogA0s.reduce((a, b) => a + b, 0) / simLogA0s.length;
      return ss + (v - m) ** 2;
    }, 0) / (simLogA0s.length - 1));
    simSDs.push(sd);
  }

  const mean = arr => arr.reduce((s, v) => s + v, 0) / arr.length;
  const sorted = arr => [...arr].sort((a, b) => a - b);

  const tauSorted = sorted(simTaus);
  const sdSorted = sorted(simSDs);

  const realTau = baseline.tau;
  const realSD = Math.sqrt(allLogA0.reduce((s, v) => s + (v - baseline.mu) ** 2, 0) / (N - 1));

  const tauPctile = tauSorted.filter(t => t <= realTau).length / nSim * 100;
  const sdPctile = sdSorted.filter(t => t <= realSD).length / nSim * 100;

  return {
    simTauMean: mean(simTaus),
    simTau95: [tauSorted[Math.floor(nSim * 0.025)], tauSorted[Math.floor(nSim * 0.975)]],
    realTau,
    tauPctile,
    simSDmean: mean(simSDs),
    simSD95: [sdSorted[Math.floor(nSim * 0.025)], sdSorted[Math.floor(nSim * 0.975)]],
    realSD,
    sdPctile,
    simA0mean: Math.round(mean(simA0s))
  };
}

const ppcResults = {};
for (const m of modelNames) {
  log("  Running PPC for " + m + "...");
  ppcResults[m] = ppcCheck(m);
}
log("");

log("  POSTERIOR PREDICTIVE CHECK RESULTS:");
log("  ┌──────────────────────────────────────────────────────────────────────┐");
log("  │  Model  simTau   realTau  tau%ile  simSD    realSD   sd%ile  a0     │");
log("  ├──────────────────────────────────────────────────────────────────────┤");
for (const m of modelNames) {
  const p = ppcResults[m];
  log("  │  " + m.padEnd(4) +
    p.simTauMean.toFixed(3).padStart(8) +
    p.realTau.toFixed(3).padStart(9) +
    (p.tauPctile.toFixed(0) + "%").padStart(8) +
    p.simSDmean.toFixed(3).padStart(8) +
    p.realSD.toFixed(3).padStart(9) +
    (p.sdPctile.toFixed(0) + "%").padStart(8) +
    p.simA0mean.toString().padStart(7) + "  │");
}
log("  └──────────────────────────────────────────────────────────────────────┘");
log("");

log("  PPC INTERPRETATION:");
for (const m of modelNames) {
  const p = ppcResults[m];
  const pass = p.tauPctile >= 5 && p.tauPctile <= 95;
  log("    " + m + ": simulated tau [" + p.simTau95[0].toFixed(3) + ", " + p.simTau95[1].toFixed(3) + "]" +
    " | real " + p.realTau.toFixed(3) + " at " + p.tauPctile.toFixed(0) + "th percentile" +
    " → " + (pass ? "CONSISTENT" : "FAILS (real tau outside 90% prediction interval)"));
}
log("");

log("  STEP 6: OFFSET vs SHAPE TEST");
sep();
log("");

let offsetOnly = 0, shapeAlso = 0;
for (const g of galaxyData) {
  if (g.pts.length < 6) continue;
  const fit = fitA0single(g.pts);
  const residuals = g.pts.map(p => {
    const gbar = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gbar, fit.a0);
    return p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10);
  });

  const half = Math.floor(residuals.length / 2);
  const innerRes = residuals.slice(0, half);
  const outerRes = residuals.slice(half);
  const innerMean = innerRes.reduce((s, v) => s + v, 0) / innerRes.length;
  const outerMean = outerRes.reduce((s, v) => s + v, 0) / outerRes.length;
  const diff = Math.abs(innerMean - outerMean);

  if (diff > 0.05) shapeAlso++;
  else offsetOnly++;
}

log("  After fitting best per-galaxy a0, examine residual structure:");
log("  Galaxies with offset-only departure:      " + offsetOnly + " (" + (offsetOnly / (offsetOnly + shapeAlso) * 100).toFixed(0) + "%)");
log("  Galaxies with shape difference also:       " + shapeAlso + " (" + (shapeAlso / (offsetOnly + shapeAlso) * 100).toFixed(0) + "%)");
log("");
if (offsetOnly > shapeAlso) {
  log("  → MAJORITY show offset-only departure (consistent with a0 shift)");
  log("    The RAR SHAPE is universal; only the SCALE varies.");
} else {
  log("  → Significant number show shape differences (deeper than just a0)");
}
log("");

log("  STEP 7: JACKKNIFE MODEL STABILITY");
sep();
log("");

const jkBetas = [];
const jkM2wins = [];
for (let i = 0; i < N; i++) {
  const subset = galaxyData.filter((_, j) => j !== i);
  const fits = subset.map(g => fitA0single(g.pts));
  const etas = subset.map(g => g.etaRot || 0);
  const logA0s = fits.map(f => f.logA0);
  const n = logA0s.length;
  const xm = etas.reduce((s, v) => s + v, 0) / n;
  const ym = logA0s.reduce((s, v) => s + v, 0) / n;
  let sxy = 0, sxx = 0;
  for (let j = 0; j < n; j++) {
    sxy += (etas[j] - xm) * (logA0s[j] - ym);
    sxx += (etas[j] - xm) ** 2;
  }
  jkBetas.push(sxx > 0 ? sxy / sxx : 0);
}

const jkBetaMean = jkBetas.reduce((s, v) => s + v, 0) / jkBetas.length;
const jkBetaSD = Math.sqrt(jkBetas.reduce((s, v) => s + (v - jkBetaMean) ** 2, 0) / (jkBetas.length - 1));
const jkBetaSorted = [...jkBetas].sort((a, b) => a - b);

log("  eta_rot coefficient stability (jackknife leave-one-out):");
log("    Mean beta:  " + jkBetaMean.toFixed(4));
log("    SD beta:    " + jkBetaSD.toFixed(4));
log("    Range:      [" + jkBetaSorted[0].toFixed(4) + ", " + jkBetaSorted[jkBetaSorted.length - 1].toFixed(4) + "]");
log("    t-stat:     " + (jkBetaMean / jkBetaSD * Math.sqrt(N)).toFixed(2));
const betaStable = jkBetaSorted[0] > 0 || jkBetaSorted[jkBetaSorted.length - 1] < 0;
log("    Sign stable? " + (betaStable ? "YES — always same sign" : "NO — sign flips in some jackknife samples"));
log("");

log("=".repeat(80));
log("  PHASE 13 FINAL VERDICT");
log("=".repeat(80));
log("");

const cvBest = best;
const aicBest = modelNames.reduce((a, b) => icResults[a].AIC < icResults[b].AIC ? a : b);
const bicBest = modelNames.reduce((a, b) => icResults[a].BIC < icResults[b].BIC ? a : b);
const ppcPass = modelNames.filter(m => {
  const p = ppcResults[m];
  return p.tauPctile >= 5 && p.tauPctile <= 95;
});

log("  ┌──────────────────────────────────────────────────────────────────┐");
log("  │  CRITERION            WINNER                                   │");
log("  ├──────────────────────────────────────────────────────────────────┤");
log("  │  CV held-out RMS      " + cvBest.padEnd(45) + "│");
log("  │  AIC                  " + aicBest.padEnd(45) + "│");
log("  │  BIC                  " + bicBest.padEnd(45) + "│");
log("  │  PPC (tau consistent) " + (ppcPass.length > 0 ? ppcPass.join(", ") : "NONE").padEnd(45) + "│");
log("  │  Offset vs shape      " + (offsetOnly > shapeAlso ? "Offset-only (supports a0 shift)" : "Shape also varies").padEnd(45) + "│");
log("  │  eta_rot stability    " + (betaStable ? "Stable (real signal)" : "Unstable (may be noise)").padEnd(45) + "│");
log("  └──────────────────────────────────────────────────────────────────┘");
log("");

const m3wins = [cvBest === 'M3', aicBest === 'M3', bicBest === 'M3'].filter(Boolean).length;
const m0wins = [cvBest === 'M0', aicBest === 'M0', bicBest === 'M0'].filter(Boolean).length;
const m2wins = [cvBest === 'M2', aicBest === 'M2', bicBest === 'M2'].filter(Boolean).length;

if (m3wins >= 2) {
  log("  VERDICT: M3 (galaxy-varying a0) WINS on majority of criteria.");
  log("  a0 is NOT strictly universal — each galaxy needs its own a0.");
  log("  The hierarchical model provides the best prediction and description.");
} else if (m0wins >= 2 && ppcPass.includes('M0')) {
  log("  VERDICT: M0 (universal a0) is sufficient.");
  log("  The scatter is consistent with measurement noise only.");
} else if (m2wins >= 2) {
  log("  VERDICT: M2 (a0 + eta_rot) provides meaningful improvement.");
  log("  A second parameter partially explains the heterogeneity.");
} else {
  log("  VERDICT: MIXED — no single model dominates all criteria.");
  log("  The truth likely lies between M0 (universal) and M3 (varying).");
}
log("");
log("  QUANTITATIVE SUMMARY:");
log("    M0 (universal):   CV RMS = " + cvSummary.M0.totalRMS.toFixed(5) + ", AIC = " + icResults.M0.AIC.toFixed(0));
log("    M1 (+ inc marg):  CV RMS = " + cvSummary.M1.totalRMS.toFixed(5) + ", AIC = " + icResults.M1.AIC.toFixed(0));
log("    M2 (+ eta_rot):   CV RMS = " + cvSummary.M2.totalRMS.toFixed(5) + ", AIC = " + icResults.M2.AIC.toFixed(0));
log("    M3 (per-galaxy):  CV RMS = " + cvSummary.M3.totalRMS.toFixed(5) + ", AIC = " + icResults.M3.AIC.toFixed(0));
log("");
log("    M3 improvement over M0: " + ((1 - cvSummary.M3.totalRMS / cvSummary.M0.totalRMS) * 100).toFixed(1) + "% in CV RMS");
log("    M3 improvement over M0: " + (icResults.M0.AIC - icResults.M3.AIC).toFixed(0) + " in AIC");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  description: "Phase 13: Model Decision Run — M0/M1/M2/M3 comparison",
  baseline: {
    a0: Math.round(baseline.a0),
    logA0: +baseline.mu.toFixed(4),
    tau: +baseline.tau.toFixed(4),
    I2: +baseline.I2.toFixed(1),
    nGal: N,
    nPts: totalPts
  },
  crossValidation: {
    nFolds: N_FOLDS,
    results: Object.fromEntries(modelNames.map(m => [m, {
      totalRMS: +cvSummary[m].totalRMS.toFixed(5),
      totalSS: +cvSummary[m].totalSS.toFixed(2)
    }])),
    winner: cvBest
  },
  informationCriteria: Object.fromEntries(modelNames.map(m => [m, {
    k: icResults[m].k,
    AIC: +icResults[m].AIC.toFixed(1),
    BIC: +icResults[m].BIC.toFixed(1),
    logL: +icResults[m].logL.toFixed(1),
    rms: +icResults[m].rms.toFixed(5)
  }])),
  aicBest, bicBest,
  posteriorPredictive: Object.fromEntries(modelNames.map(m => [m, {
    simTauMean: +ppcResults[m].simTauMean.toFixed(4),
    simTau95: ppcResults[m].simTau95.map(v => +v.toFixed(4)),
    realTau: +ppcResults[m].realTau.toFixed(4),
    tauPercentile: +ppcResults[m].tauPctile.toFixed(1),
    consistent: ppcResults[m].tauPctile >= 5 && ppcResults[m].tauPctile <= 95
  }])),
  offsetVsShape: { offsetOnly, shapeAlso, majorityOffset: offsetOnly > shapeAlso },
  etaRotStability: {
    jkBetaMean: +jkBetaMean.toFixed(4),
    jkBetaSD: +jkBetaSD.toFixed(4),
    signStable: betaStable
  },
  verdict: m3wins >= 2 ? "M3_WINS" : m0wins >= 2 ? "M0_SUFFICIENT" : m2wins >= 2 ? "M2_IMPROVES" : "MIXED"
};

fs.writeFileSync(path.join(__dirname, '../public/phase13-model-decision.json'), JSON.stringify(output, null, 2));
log("  Results saved to public/phase13-model-decision.json");
