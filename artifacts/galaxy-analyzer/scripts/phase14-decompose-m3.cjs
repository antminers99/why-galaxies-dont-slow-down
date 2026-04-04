#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "14.0.0";

function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(80)); }

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function generalizedRAR(gbar, a0, gamma) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.pow(y, gamma / 2)));
}

function randn() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

function fitA0(pts) {
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

function fitA0gamma(pts) {
  let bestA0 = 3633, bestGamma = 0.5, bestSS = Infinity;
  for (let la = 2.5; la <= 4.5; la += 0.02) {
    for (let g = 0.2; g <= 1.5; g += 0.05) {
      let ss = 0;
      const a0 = Math.pow(10, la);
      for (const p of pts) {
        const gbar = Math.pow(10, p.log_g_bar);
        const pred = generalizedRAR(gbar, a0, g);
        if (!isFinite(pred) || pred <= 0) { ss += 100; continue; }
        ss += (p.log_g_obs - Math.log10(pred)) ** 2;
      }
      if (ss < bestSS) { bestSS = ss; bestA0 = a0; bestGamma = g; }
    }
  }
  let lo_a = Math.log10(bestA0) - 0.05, hi_a = Math.log10(bestA0) + 0.05;
  let lo_g = bestGamma - 0.1, hi_g = bestGamma + 0.1;
  lo_g = Math.max(0.1, lo_g); hi_g = Math.min(2.0, hi_g);
  for (let iter = 0; iter < 50; iter++) {
    const mids_a = [lo_a + (hi_a-lo_a)*0.33, lo_a + (hi_a-lo_a)*0.67];
    const mids_g = [lo_g + (hi_g-lo_g)*0.33, lo_g + (hi_g-lo_g)*0.67];
    let bss = Infinity, ba = mids_a[0], bg = mids_g[0];
    for (const la of mids_a) {
      for (const gg of mids_g) {
        let ss = 0;
        const a0 = Math.pow(10, la);
        for (const p of pts) {
          const gbar = Math.pow(10, p.log_g_bar);
          const pred = generalizedRAR(gbar, a0, gg);
          if (!isFinite(pred) || pred <= 0) { ss += 100; continue; }
          ss += (p.log_g_obs - Math.log10(pred)) ** 2;
        }
        if (ss < bss) { bss = ss; ba = la; bg = gg; }
      }
    }
    const ra = (hi_a - lo_a) * 0.2;
    const rg = (hi_g - lo_g) * 0.2;
    lo_a = ba - ra; hi_a = ba + ra;
    lo_g = Math.max(0.05, bg - rg); hi_g = Math.min(3.0, bg + rg);
    bestA0 = Math.pow(10, ba); bestGamma = bg; bestSS = bss;
  }
  return { a0: bestA0, logA0: Math.log10(bestA0), gamma: bestGamma, rms: Math.sqrt(bestSS / pts.length), ss: bestSS, n: pts.length };
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

function dlMeta(logA0s, ses) {
  const n = logA0s.length;
  const w = ses.map(s => 1 / (s * s));
  const wSum = w.reduce((a, b) => a + b, 0);
  const muFE = logA0s.reduce((s, v, i) => s + w[i] * v, 0) / wSum;
  let Q = 0;
  for (let i = 0; i < n; i++) Q += w[i] * (logA0s[i] - muFE) ** 2;
  const S1 = wSum, S2 = w.reduce((s, v) => s + v * v, 0);
  const tau2 = Math.max(0, (Q - (n - 1)) / (S1 - S2 / S1));
  return { mu: muFE, tau: Math.sqrt(tau2), tau2 };
}

const p11 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase11-sensitivity-lab.json'), 'utf8'));
const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));

const galaxyData = [];
for (const gal of p11.galaxies) {
  const pts = gal.localProfile.filter(p =>
    isFinite(p.log_g_bar) && isFinite(p.log_g_obs) && p.log_g_obs > p.log_g_bar * 0.5
  );
  if (pts.length < 5) continue;

  const sparcEntry = sparc.find(s => (s.name || s.Galaxy) === gal.name);

  let etaRot = 0;
  if (pts.length >= 4) {
    const oh = pts.slice(Math.floor(pts.length / 2));
    const ih = pts.slice(0, Math.floor(pts.length / 2));
    if (oh.length >= 2 && ih.length >= 2) {
      const os = (oh[oh.length-1].log_g_obs - oh[0].log_g_obs) / (oh[oh.length-1].log_g_bar - oh[0].log_g_bar + 1e-10);
      const is_ = (ih[ih.length-1].log_g_obs - ih[0].log_g_obs) / (ih[ih.length-1].log_g_bar - ih[0].log_g_bar + 1e-10);
      etaRot = os - is_;
    }
  }

  galaxyData.push({
    name: gal.name, pts, inc: gal.inc, D: gal.D, Vmax: gal.Vmax,
    T: gal.T || (sparcEntry ? sparcEntry.T : 5),
    etaRot, n: pts.length,
    logMHI: sparcEntry ? Math.log10((sparcEntry.MHI || 1e9)) : 9,
    Rdisk: sparcEntry ? (sparcEntry.Rdisk || 3) : 3,
    SBdisk: sparcEntry ? (sparcEntry.SBdisk || 100) : 100,
    SBeff: sparcEntry ? (sparcEntry.SBeff || 100) : 100,
    logL36: sparcEntry ? Math.log10((sparcEntry.L || 1e9)) : 9,
  });
}

const N = galaxyData.length;

log("");
log("=".repeat(80));
log("  PHASE 14: DECOMPOSE M3 — WHY DOES EACH GALAXY DIFFER?");
log("  Version " + VERSION);
log("  " + N + " galaxies");
log("=".repeat(80));
log("");

log("  STEP 1: FIT EACH GALAXY — OFFSET (a0) vs SHAPE (gamma)");
sep();
log("");

const results = [];
for (const g of galaxyData) {
  const fitSimple = fitA0(g.pts);
  const fitFull = fitA0gamma(g.pts);
  const ssM0 = predictSS(3656, g.pts);
  const rmsM0 = Math.sqrt(ssM0 / g.pts.length);

  results.push({
    name: g.name,
    logA0: fitSimple.logA0,
    a0: fitSimple.a0,
    rmsA0only: fitSimple.rms,
    gamma: fitFull.gamma,
    logA0_full: fitFull.logA0,
    rmsFull: fitFull.rms,
    rmsM0: rmsM0,
    n: g.n,
    inc: g.inc, D: g.D, Vmax: g.Vmax, T: g.T,
    etaRot: g.etaRot,
    logMHI: g.logMHI, Rdisk: g.Rdisk, SBdisk: g.SBdisk,
    logL36: g.logL36,
    deltaA0: fitSimple.logA0 - 3.563,
    deltaGamma: fitFull.gamma - 0.5,
    shapeImprovement: ((fitSimple.rms - fitFull.rms) / fitSimple.rms * 100)
  });
}

const meanGamma = results.reduce((s, r) => s + r.gamma, 0) / N;
const sdGamma = Math.sqrt(results.reduce((s, r) => s + (r.gamma - meanGamma) ** 2, 0) / (N - 1));

log("  OFFSET (a0) distribution:");
const deltaA0s = results.map(r => r.deltaA0);
const meanDA = deltaA0s.reduce((s, v) => s + v, 0) / N;
const sdDA = Math.sqrt(deltaA0s.reduce((s, v) => s + (v - meanDA) ** 2, 0) / (N - 1));
log("    Mean delta log(a0): " + meanDA.toFixed(4) + " dex");
log("    SD delta log(a0):   " + sdDA.toFixed(4) + " dex");
log("    Range: [" + Math.min(...deltaA0s).toFixed(3) + ", " + Math.max(...deltaA0s).toFixed(3) + "]");
log("");

log("  SHAPE (gamma) distribution:");
log("    Standard RAR has gamma = 0.5");
log("    Mean gamma: " + meanGamma.toFixed(3));
log("    SD gamma:   " + sdGamma.toFixed(3));
log("    Range: [" + Math.min(...results.map(r => r.gamma)).toFixed(3) + ", " + Math.max(...results.map(r => r.gamma)).toFixed(3) + "]");
log("");

const shapeImps = results.map(r => r.shapeImprovement);
const meanShapeImp = shapeImps.reduce((s, v) => s + v, 0) / N;
log("  Adding gamma (shape freedom) improves fit by " + meanShapeImp.toFixed(1) + "% on average");
log("");

const ssA0only = results.reduce((s, r) => s + r.rmsA0only ** 2 * r.n, 0);
const ssFull = results.reduce((s, r) => s + r.rmsFull ** 2 * r.n, 0);
const ssM0total = results.reduce((s, r) => s + r.rmsM0 ** 2 * r.n, 0);
const totalN = results.reduce((s, r) => s + r.n, 0);

const varM0 = ssM0total / totalN;
const varA0only = ssA0only / totalN;
const varFull = ssFull / totalN;

const offsetFraction = (varM0 - varA0only) / (varM0 - varFull) * 100;
const shapeFraction = (varA0only - varFull) / (varM0 - varFull) * 100;

log("  VARIANCE DECOMPOSITION:");
log("    M0 (universal a0):           variance = " + varM0.toFixed(6));
log("    M3-offset (per-galaxy a0):   variance = " + varA0only.toFixed(6));
log("    M3-full (a0 + gamma):        variance = " + varFull.toFixed(6));
log("");
log("  ┌────────────────────────────────────────────────────────┐");
log("  │  M3 improvement breakdown:                            │");
log("  │    Offset (a0 shift):  " + offsetFraction.toFixed(1).padStart(5) + "% of M3's advantage      │");
log("  │    Shape (gamma):      " + shapeFraction.toFixed(1).padStart(5) + "% of M3's advantage      │");
log("  └────────────────────────────────────────────────────────┘");
log("");

log("  STEP 2: WHAT GALAXY PROPERTIES PREDICT delta_a0?");
sep();
log("  Model: log(a0,i) = log(a0) + f(Xi)");
log("");

const predictors = [
  { name: 'log(Vmax)', values: results.map(r => Math.log10(r.Vmax)) },
  { name: 'inc', values: results.map(r => r.inc) },
  { name: 'log(D)', values: results.map(r => Math.log10(r.D)) },
  { name: 'T (Hubble)', values: results.map(r => r.T) },
  { name: 'eta_rot', values: results.map(r => r.etaRot) },
  { name: 'log(MHI)', values: results.map(r => r.logMHI) },
  { name: 'Rdisk', values: results.map(r => r.Rdisk) },
  { name: 'SBdisk', values: results.map(r => r.SBdisk) },
  { name: 'log(L36)', values: results.map(r => r.logL36) },
  { name: 'n_points', values: results.map(r => r.n) },
  { name: 'gamma', values: results.map(r => r.gamma) },
];

const target = results.map(r => r.deltaA0);
const corrResults = [];

for (const pred of predictors) {
  const x = pred.values;
  const y = target;
  const n = x.length;
  const xm = x.reduce((s, v) => s + v, 0) / n;
  const ym = y.reduce((s, v) => s + v, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - xm) * (y[i] - ym);
    sxx += (x[i] - xm) ** 2;
    syy += (y[i] - ym) ** 2;
  }
  const r = sxy / Math.sqrt(sxx * syy + 1e-20);
  const beta = sxx > 0 ? sxy / sxx : 0;
  const seBeta = Math.sqrt((syy - sxy * sxy / sxx) / ((n - 2) * sxx + 1e-20));
  const tStat = seBeta > 0 ? beta / seBeta : 0;
  const r2 = r * r;

  const residSS = y.reduce((s, v, i) => s + (v - (ym + beta * (x[i] - xm))) ** 2, 0);
  const baseSS = syy;
  const tauReduction = baseSS > 0 ? (1 - residSS / baseSS) * 100 : 0;

  corrResults.push({ name: pred.name, r: +r.toFixed(3), r2: +r2.toFixed(3), beta: +beta.toFixed(4), tStat: +tStat.toFixed(2), tauReduction: +tauReduction.toFixed(1) });
}

corrResults.sort((a, b) => Math.abs(b.r) - Math.abs(a.r));

log("  ┌──────────────────────────────────────────────────────────────────┐");
log("  │  Predictor X     r        r²     beta     |t|    tau-red(%)    │");
log("  ├──────────────────────────────────────────────────────────────────┤");
for (const c of corrResults) {
  const sig = Math.abs(c.tStat) > 2.0 ? " **" : Math.abs(c.tStat) > 1.65 ? " *" : "";
  log("  │  " + c.name.padEnd(14) + c.r.toFixed(3).padStart(7) + (c.r2 * 100).toFixed(1).padStart(8) + "%" +
    c.beta.toFixed(4).padStart(9) + Math.abs(c.tStat).toFixed(1).padStart(7) +
    c.tauReduction.toFixed(1).padStart(10) + "%" + sig.padEnd(4) + "│");
}
log("  └──────────────────────────────────────────────────────────────────┘");
log("  ** p < 0.05, * p < 0.10");
log("");

log("  STEP 3: MULTIVARIATE MODEL — a0,i = a0 + f(X1, X2, ...)");
sep();
log("");

const sigPreds = corrResults.filter(c => Math.abs(c.tStat) >= 1.65);
log("  Significant predictors (|t| >= 1.65):");
for (const p of sigPreds) {
  log("    " + p.name + " (r=" + p.r.toFixed(3) + ", tau-red=" + p.tauReduction.toFixed(1) + "%)");
}
log("");

if (sigPreds.length >= 2) {
  const top2names = sigPreds.slice(0, 2).map(p => p.name);
  log("  Top-2 multivariate model: delta_a0 = b1*" + top2names[0] + " + b2*" + top2names[1]);

  const x1 = predictors.find(p => p.name === top2names[0]).values;
  const x2 = predictors.find(p => p.name === top2names[1]).values;
  const y = target;
  const n = y.length;

  const mx1 = x1.reduce((s,v)=>s+v,0)/n, mx2 = x2.reduce((s,v)=>s+v,0)/n, my = y.reduce((s,v)=>s+v,0)/n;
  let s11=0,s12=0,s1y=0,s22=0,s2y=0;
  for (let i=0;i<n;i++) {
    const d1=x1[i]-mx1, d2=x2[i]-mx2, dy=y[i]-my;
    s11+=d1*d1; s12+=d1*d2; s1y+=d1*dy; s22+=d2*d2; s2y+=d2*dy;
  }
  const det = s11*s22-s12*s12;
  const b1 = det>0 ? (s22*s1y-s12*s2y)/det : 0;
  const b2 = det>0 ? (s11*s2y-s12*s1y)/det : 0;
  const a = my - b1*mx1 - b2*mx2;

  let residSS = 0, totalSS = 0;
  for (let i=0;i<n;i++) {
    const pred = a + b1*x1[i] + b2*x2[i];
    residSS += (y[i]-pred)**2;
    totalSS += (y[i]-my)**2;
  }
  const adjR2 = 1 - (residSS/(n-3)) / (totalSS/(n-1));
  const multiTauRed = (1 - residSS/totalSS) * 100;

  log("    b1 = " + b1.toFixed(4) + ", b2 = " + b2.toFixed(4) + ", intercept = " + a.toFixed(4));
  log("    R² = " + (1-residSS/totalSS).toFixed(3) + ", adj-R² = " + adjR2.toFixed(3));
  log("    Tau reduction: " + multiTauRed.toFixed(1) + "%");
  log("");
}

log("  STEP 4: CROSS-VALIDATED PREDICTIVE MODEL");
sep();
log("");
log("  Test: can we PREDICT a galaxy's a0 from its properties?");
log("  Leave-one-out cross-validation for top predictor.");
log("");

const bestPred = corrResults[0];
const bestX = predictors.find(p => p.name === bestPred.name).values;

let looPredSS = 0, looBaseSS = 0;
for (let i = 0; i < N; i++) {
  const trainX = bestX.filter((_, j) => j !== i);
  const trainY = target.filter((_, j) => j !== i);
  const n = trainX.length;
  const mx = trainX.reduce((s,v)=>s+v,0)/n;
  const my = trainY.reduce((s,v)=>s+v,0)/n;
  let sxy = 0, sxx = 0;
  for (let j = 0; j < n; j++) { sxy += (trainX[j]-mx)*(trainY[j]-my); sxx += (trainX[j]-mx)**2; }
  const beta = sxx > 0 ? sxy / sxx : 0;
  const alpha = my - beta * mx;
  const pred = alpha + beta * bestX[i];
  looPredSS += (target[i] - pred) ** 2;
  const baseMean = trainY.reduce((s,v)=>s+v,0)/n;
  looBaseSS += (target[i] - baseMean) ** 2;
}

const looR2 = 1 - looPredSS / looBaseSS;
log("  Best predictor: " + bestPred.name);
log("  LOO-CV R² = " + looR2.toFixed(3));
log("  LOO-CV RMS of delta_a0: " + Math.sqrt(looPredSS / N).toFixed(4) + " dex");
log("  Baseline RMS of delta_a0: " + Math.sqrt(looBaseSS / N).toFixed(4) + " dex");
log("");

if (looR2 > 0.05) {
  log("  → " + bestPred.name + " DOES predict a0 variation out-of-sample (R²=" + looR2.toFixed(3) + ")");
} else {
  log("  → " + bestPred.name + " does NOT predict a0 well out-of-sample");
}
log("");

log("  STEP 5: IS GAMMA CORRELATED WITH ANYTHING?");
sep();
log("");

const gammas = results.map(r => r.gamma);
const gammaCorrs = [];

for (const pred of predictors) {
  if (pred.name === 'gamma') continue;
  const x = pred.values;
  const n = x.length;
  const xm = x.reduce((s,v)=>s+v,0)/n;
  const gm = gammas.reduce((s,v)=>s+v,0)/n;
  let sxy=0,sxx=0,syy=0;
  for (let i=0;i<n;i++) { sxy+=(x[i]-xm)*(gammas[i]-gm); sxx+=(x[i]-xm)**2; syy+=(gammas[i]-gm)**2; }
  const r = sxy / Math.sqrt(sxx*syy+1e-20);
  gammaCorrs.push({ name: pred.name, r: +r.toFixed(3) });
}
gammaCorrs.sort((a,b) => Math.abs(b.r) - Math.abs(a.r));

log("  Correlations of gamma (shape parameter) with galaxy properties:");
for (const c of gammaCorrs) {
  const bar = "#".repeat(Math.round(Math.abs(c.r) * 30));
  log("    " + c.name.padEnd(14) + " r=" + c.r.toFixed(3).padStart(7) + "  " + bar);
}
log("");

const gammaA0corr = (() => {
  const x = results.map(r => r.deltaA0);
  const y = gammas;
  const n = x.length;
  const xm = x.reduce((s,v)=>s+v,0)/n;
  const ym = y.reduce((s,v)=>s+v,0)/n;
  let sxy=0,sxx=0,syy=0;
  for (let i=0;i<n;i++) { sxy+=(x[i]-xm)*(y[i]-ym); sxx+=(x[i]-xm)**2; syy+=(y[i]-ym)**2; }
  return sxy / Math.sqrt(sxx*syy+1e-20);
})();

log("  Correlation between delta_a0 and gamma: r = " + gammaA0corr.toFixed(3));
if (Math.abs(gammaA0corr) > 0.3) {
  log("  → a0 shift and shape change are CORRELATED — they may share a common cause");
} else {
  log("  → a0 shift and shape change are INDEPENDENT — two separate effects");
}
log("");

log("=".repeat(80));
log("  PHASE 14 CONCLUSIONS");
log("=".repeat(80));
log("");
log("  1. DECOMPOSITION OF M3's ADVANTAGE:");
log("     Offset (a0 shift):      " + offsetFraction.toFixed(1) + "% of improvement");
log("     Shape (gamma variation): " + shapeFraction.toFixed(1) + "% of improvement");
log("");

if (offsetFraction > 70) {
  log("     → M3 mostly wins through OFFSET — galaxies differ in scale, not shape.");
  log("     → The RAR transition is ~universal; only the acceleration scale shifts.");
} else if (shapeFraction > 70) {
  log("     → M3 mostly wins through SHAPE — galaxies differ in transition curvature.");
  log("     → This is DEEPER than just a0 variation.");
} else {
  log("     → BOTH offset and shape contribute significantly.");
  log("     → Galaxies differ in both scale AND transition curvature.");
}
log("");

log("  2. BEST PREDICTOR OF a0 VARIATION:");
log("     " + bestPred.name + " (r=" + bestPred.r.toFixed(3) + ", LOO-CV R²=" + looR2.toFixed(3) + ")");
if (sigPreds.length > 0) {
  log("     Equation: a0,i = a0 + " + bestPred.beta.toFixed(4) + " × (" + bestPred.name + " - mean)");
}
log("");

log("  3. THE EQUATION a0,i = a0 + f(Xi):");
if (looR2 > 0.1) {
  log("     f(Xi) is DETECTABLE — " + bestPred.name + " predicts a0 variation out-of-sample.");
  log("     This is the missing variable.");
} else if (looR2 > 0.02) {
  log("     f(Xi) is WEAK but present — " + bestPred.name + " has some predictive power.");
  log("     The missing variable is partially captured but not fully identified.");
} else {
  log("     f(Xi) is NOT FOUND — no measured property predicts a0 variation well.");
  log("     The variation is either random, physical, or driven by an unmeasured variable.");
}
log("");

log("  4. GAMMA (SHAPE) VARIATION:");
log("     Mean gamma = " + meanGamma.toFixed(3) + " ± " + sdGamma.toFixed(3) + " (standard = 0.5)");
log("     Gamma-a0 correlation: r = " + gammaA0corr.toFixed(3));
log("");

log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  description: "Phase 14: Decompose M3 — Why does each galaxy differ?",
  nGalaxies: N,
  varianceDecomposition: {
    M0_variance: +varM0.toFixed(6),
    M3offset_variance: +varA0only.toFixed(6),
    M3full_variance: +varFull.toFixed(6),
    offsetFraction: +offsetFraction.toFixed(1),
    shapeFraction: +shapeFraction.toFixed(1)
  },
  gammaDistribution: { mean: +meanGamma.toFixed(3), sd: +sdGamma.toFixed(3) },
  a0Predictors: corrResults,
  bestPredictor: { name: bestPred.name, r: bestPred.r, looR2: +looR2.toFixed(3) },
  gammaCorrelations: gammaCorrs,
  gammaA0correlation: +gammaA0corr.toFixed(3),
  perGalaxy: results.map(r => ({
    name: r.name, logA0: +r.logA0.toFixed(4), deltaA0: +r.deltaA0.toFixed(4),
    gamma: +r.gamma.toFixed(3), rmsA0only: +r.rmsA0only.toFixed(4), rmsFull: +r.rmsFull.toFixed(4),
    shapeImprovement: +r.shapeImprovement.toFixed(1),
    Vmax: r.Vmax, inc: r.inc, T: r.T, etaRot: +r.etaRot.toFixed(3)
  }))
};

fs.writeFileSync(path.join(__dirname, '../public/phase14-decompose-m3.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase14-decompose-m3.json");
