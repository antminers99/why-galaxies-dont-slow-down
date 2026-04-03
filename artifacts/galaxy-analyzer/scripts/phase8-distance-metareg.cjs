const fs = require('fs');
const path = require('path');

const VERSION = "8.0.0";
const TIMESTAMP = new Date().toISOString();

const UPSILON_FID = 0.5;
const MS2_CONV = 3.241e-14;
const C_MS = 299792458;
const H0_SI = 67.4e3 / 3.086e22;
const CH0_2PI = C_MS * H0_SI / (2 * Math.PI);
const CH0_2PI_GAL = CH0_2PI / MS2_CONV;

const rarData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const tsData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'transition-scale.json'), 'utf8'));
const sparcTable = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf8'));

function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(76)); }
function pad(s, n) { return String(s).padEnd(n); }
function padr(s, n) { return String(s).padStart(n); }

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function evalChi2(logGbar, logGobs, a0) {
  let chi2 = 0;
  for (let i = 0; i < logGbar.length; i++) {
    const gbar = Math.pow(10, logGbar[i]);
    const pred = mcgaughRAR(gbar, a0);
    if (!isFinite(pred) || pred <= 0) { chi2 += 100; continue; }
    const r = logGobs[i] - Math.log10(pred);
    chi2 += r * r;
  }
  return chi2;
}

const coarseGrid = [];
for (let logA = 2.0; logA <= 5.0; logA += 0.01) coarseGrid.push(logA);

function fitA0(logGbar, logGobs) {
  if (logGbar.length < 5) return null;
  let bestA0 = 3000, bestChi2 = Infinity;
  for (const logA of coarseGrid) {
    const a0 = Math.pow(10, logA);
    const chi2 = evalChi2(logGbar, logGobs, a0);
    if (chi2 < bestChi2) { bestChi2 = chi2; bestA0 = a0; }
  }
  let lo = Math.log10(bestA0) - 0.02, hi = Math.log10(bestA0) + 0.02;
  lo = Math.max(lo, 2.0); hi = Math.min(hi, 5.0);
  for (let step = 0; step < 100; step++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    const c1 = evalChi2(logGbar, logGobs, Math.pow(10, m1));
    const c2 = evalChi2(logGbar, logGobs, Math.pow(10, m2));
    if (c1 < c2) hi = m2; else lo = m1;
  }
  const finalA0 = Math.pow(10, (lo + hi) / 2);
  const finalChi2 = evalChi2(logGbar, logGobs, finalA0);
  const rms = Math.sqrt(finalChi2 / logGbar.length);
  return { a0: finalA0, logA0: Math.log10(finalA0), chi2: finalChi2, rms, n: logGbar.length };
}

function gaussPrior(x, mu, sigma) {
  return Math.exp(-0.5 * ((x - mu) / sigma) ** 2);
}

const upsilonGrid = [];
for (let u = 0.20; u <= 0.801; u += 0.05) upsilonGrid.push(+u.toFixed(2));
const dlogDGrid = [-0.10, -0.05, 0.00, 0.05, 0.10];

function marginalizeSingleGalaxy(gal) {
  if (gal.points.length < 5) return null;
  const baseLogGbar = gal.points.map(p => p.logGbar);
  const baseLogGobs = gal.points.map(p => p.logGobs);
  const gridResults = [];
  for (const upsilon of upsilonGrid) {
    const scale = upsilon / UPSILON_FID;
    const gbarAdj = (gal.isLT || gal.fgas >= 0.99)
      ? 0.0
      : Math.log10(scale * (1 - gal.fgas) + gal.fgas);
    const adjLogGbar = baseLogGbar.map(v => v + gbarAdj);
    for (const dlogD of dlogDGrid) {
      const adjLogGobs = baseLogGobs.map(v => v - dlogD);
      const fit = fitA0(adjLogGbar, adjLogGobs);
      if (!fit || fit.a0 > 1e5 || fit.a0 < 10) continue;
      const priorW = gaussPrior(upsilon, 0.50, 0.12) * gaussPrior(dlogD, 0.0, 0.07);
      gridResults.push({ ...fit, upsilon, dlogD, priorW });
    }
  }
  if (gridResults.length === 0) return null;
  const minChi2 = Math.min(...gridResults.map(r => r.chi2));
  const n = gal.points.length;
  const sigma2Hat = Math.max(minChi2 / (n - 1), 1e-6);
  const bestRMS = Math.sqrt(minChi2 / n);
  const withinSE2 = bestRMS * bestRMS / n;
  let sumW = 0, sumWlogA = 0, sumWlogA2 = 0;
  for (const r of gridResults) {
    const relChi2 = r.chi2 - minChi2;
    const w = Math.exp(-0.5 * relChi2 / sigma2Hat) * r.priorW;
    if (!isFinite(w)) continue;
    sumW += w;
    sumWlogA += w * r.logA0;
    sumWlogA2 += w * r.logA0 * r.logA0;
  }
  if (sumW <= 0) return null;
  const meanLogA = sumWlogA / sumW;
  const varLogA = sumWlogA2 / sumW - meanLogA * meanLogA;
  const totalSE2 = varLogA + withinSE2;
  const seLogA = Math.sqrt(Math.max(totalSE2, 1e-8));
  const gRange = Math.max(...baseLogGbar) - Math.min(...baseLogGbar);
  return {
    name: gal.name, logA0: meanLogA, a0: Math.pow(10, meanLogA),
    se: seLogA, n: gal.points.length, Vmax: gal.Vmax,
    distance: gal.distance, inc: gal.inc, fgas: gal.fgas,
    gRange, Q_sparc: gal.Q_sparc, fD: gal.fD, isLT: !!gal.isLT
  };
}

const sparcLookup = {};
for (const g of sparcTable) sparcLookup[g.name] = g;

const galaxyMap = {};
for (const g of rarData.perGalaxy) {
  const Mgas = 1.33 * g.MHI;
  const Mstar = UPSILON_FID * g.L36;
  const fgas = (Mgas + Mstar > 0) ? Mgas / (Mgas + Mstar) : 0.5;
  const sparc = sparcLookup[g.name] || {};
  galaxyMap[g.name] = {
    name: g.name, inc: g.inc, distance: g.distance,
    L36: g.L36, MHI: g.MHI, Vmax: g.Vmax, fgas: fgas,
    Q_sparc: sparc.Q || null, fD: sparc.fD || null,
    T: sparc.T !== undefined ? sparc.T : null,
    points: []
  };
}

for (const p of rarData.rarScatter) {
  if (galaxyMap[p.name]) {
    galaxyMap[p.name].points.push({ logGbar: p.log_g_bar, logGobs: p.log_g_obs });
  }
}

for (const p of tsData.plotPoints) {
  if (p.g.startsWith('LT_')) {
    if (!galaxyMap[p.g]) {
      galaxyMap[p.g] = {
        name: p.g, inc: 60, distance: 5.0, L36: 0, MHI: 0,
        Vmax: 50, fgas: 0.5, Q_sparc: null, fD: null, T: null,
        points: [], isLT: true
      };
    }
    galaxyMap[p.g].points.push({ logGbar: p.x, logGobs: p.x + p.y });
  }
}

const allGalaxies = Object.values(galaxyMap).filter(g => g.points.length >= 5);

log("================================================================================");
log("  PHASE 8: DISTANCE-METHOD MATCHED META-REGRESSION");
log("================================================================================");
log("  Version: " + VERSION);
log("  Date: " + TIMESTAMP);
log("  Purpose: Test whether precise vs Hubble-flow distance method effect on a0");
log("  persists after controlling for Vmax, inclination, n, dynamic range, Q.");
log("");

log("  Marginalizing all galaxies...");
const margAll = [];
for (const gal of allGalaxies) {
  const est = marginalizeSingleGalaxy(gal);
  if (est) margAll.push(est);
}
log("  Total marginalized: " + margAll.length);

const goldI45 = margAll.filter(g =>
  g.Vmax >= 50 && g.inc >= 45 && g.n >= 5 && g.gRange >= 1.0
);
log("  GOLD+i45 sample: " + goldI45.length);

const isPrecise = g => g.fD === 2 || g.fD === 3 || g.fD === 5;
const preciseGals = goldI45.filter(isPrecise);
const hubbleGals = goldI45.filter(g => g.fD === 1);
const otherGals = goldI45.filter(g => !isPrecise(g) && g.fD !== 1);

log("");
sep();
log("  SAMPLE BREAKDOWN BY DISTANCE METHOD");
sep();
log("");
log("  Precise (TRGB/Cepheid/SN, fD=2,3,5): " + preciseGals.length);
log("  Hubble-flow (fD=1):                   " + hubbleGals.length);
log("  Other (UMa cluster, fD=4):            " + otherGals.length);

function summaryStats(arr, label) {
  if (!arr.length) return;
  const mean = v => v.reduce((a, b) => a + b, 0) / v.length;
  const med = v => { const s = [...v].sort((a,b) => a-b); const m = Math.floor(s.length/2); return s.length%2 ? s[m] : (s[m-1]+s[m])/2; };
  const logA0s = arr.map(g => g.logA0);
  const vmaxs = arr.map(g => g.Vmax);
  const incs = arr.map(g => g.inc);
  const ns = arr.map(g => g.n);
  const grs = arr.map(g => g.gRange);
  const qs = arr.filter(g => g.Q_sparc === 1).length;
  log("  " + pad(label, 30) +
    " medLogA0=" + padr(med(logA0s).toFixed(3), 6) +
    " medVmax=" + padr(Math.round(med(vmaxs)), 4) +
    " medInc=" + padr(Math.round(med(incs)), 3) +
    " medN=" + padr(Math.round(med(ns)), 3) +
    " medGR=" + padr(med(grs).toFixed(2), 5) +
    " Q1%=" + padr((qs/arr.length*100).toFixed(0) + "%", 4));
}

log("");
summaryStats(preciseGals, "Precise");
summaryStats(hubbleGals, "Hubble-flow");
summaryStats(otherGals, "UMa cluster");

log("");
sep();
log("  STEP 1: UNADJUSTED COMPARISON");
sep();
log("");

function hierarchicalDL(galEsts, label) {
  if (galEsts.length < 3) return null;
  const logA0s = galEsts.map(g => g.logA0);
  const ses = galEsts.map(g => g.se);
  const ws = ses.map(s => 1 / (s * s));
  const sumW = ws.reduce((a, b) => a + b, 0);
  const muFE = ws.reduce((acc, w, i) => acc + w * logA0s[i], 0) / sumW;
  let Q = 0;
  for (let i = 0; i < logA0s.length; i++) Q += ws[i] * (logA0s[i] - muFE) ** 2;
  const k = galEsts.length;
  const sumW2 = ws.reduce((a, w) => a + w * w, 0);
  const C = sumW - sumW2 / sumW;
  let tau2 = Math.max(0, (Q - (k - 1)) / C);
  const wsRE = ses.map(s => 1 / (s * s + tau2));
  const sumWRE = wsRE.reduce((a, b) => a + b, 0);
  const muRE = wsRE.reduce((acc, w, i) => acc + w * logA0s[i], 0) / sumWRE;
  const seRE = Math.sqrt(1 / sumWRE);
  const tau = Math.sqrt(tau2);
  const I2 = Q > 0 ? Math.max(0, (Q - (k - 1)) / Q * 100) : 0;
  return {
    label, nGal: k, a0: Math.round(Math.pow(10, muRE)),
    logA0: +muRE.toFixed(4), se: +seRE.toFixed(4),
    tau: +tau.toFixed(4), I2: +I2.toFixed(1),
    lo1sig: Math.round(Math.pow(10, muRE - Math.sqrt(seRE*seRE + tau2))),
    hi1sig: Math.round(Math.pow(10, muRE + Math.sqrt(seRE*seRE + tau2)))
  };
}

const rPrec = hierarchicalDL(preciseGals, "Precise (fD=2,3,5)");
const rHub = hierarchicalDL(hubbleGals, "Hubble-flow (fD=1)");
const rOther = hierarchicalDL(otherGals, "UMa cluster (fD=4)");
const rAll = hierarchicalDL(goldI45, "GOLD+i45 (all)");

function printR(r) {
  if (!r) { log("    (insufficient data)"); return; }
  log("    " + pad(r.label, 30) + " n=" + padr(r.nGal, 3) +
    " a0=" + padr(r.a0, 5) + " logA0=" + padr(r.logA0, 7) +
    " tau=" + padr(r.tau, 6) + " 68%=[" + r.lo1sig + "," + r.hi1sig + "]");
}

printR(rAll);
printR(rPrec);
printR(rHub);
printR(rOther);

if (rPrec && rHub) {
  const rawDiff = rPrec.logA0 - rHub.logA0;
  const rawSE = Math.sqrt(rPrec.se * rPrec.se + rHub.se * rHub.se);
  const rawT = rawDiff / rawSE;
  log("");
  log("    Raw difference (precise - Hubble): " + (rawDiff >= 0 ? "+" : "") + rawDiff.toFixed(4) + " dex");
  log("    SE of difference: " + rawSE.toFixed(4));
  log("    t-statistic: " + rawT.toFixed(2));
  log("    Significance: " + (Math.abs(rawT) >= 2.0 ? "YES (|t| >= 2)" : "NO (|t| < 2)"));
}

log("");
sep();
log("  STEP 2: WEIGHTED LEAST SQUARES META-REGRESSION");
log("  Model: logA0_i = b0 + b1*isPrecise + b2*logVmax + b3*inc + b4*logN");
log("         + b5*gRange + b6*isQ1 + error_i");
log("  Weight: 1 / (se_i^2 + tau^2)");
sep();
log("");

function wlsMetaReg(gals, covariateNames, covariateFns, tau2) {
  const n = gals.length;
  const p = covariateNames.length + 1;
  if (n <= p + 1) return null;

  const y = gals.map(g => g.logA0);
  const W = gals.map(g => 1 / (g.se * g.se + tau2));

  const X = [];
  for (let i = 0; i < n; i++) {
    const row = [1];
    for (const fn of covariateFns) row.push(fn(gals[i]));
    X.push(row);
  }

  const XtWX = Array.from({ length: p }, () => new Float64Array(p));
  const XtWy = new Float64Array(p);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < p; j++) {
      XtWy[j] += X[i][j] * W[i] * y[i];
      for (let k = 0; k < p; k++) {
        XtWX[j][k] += X[i][j] * W[i] * X[i][k];
      }
    }
  }

  const beta = solveLinear(XtWX, XtWy, p);
  if (!beta) return null;

  let ssRes = 0;
  const resids = [];
  for (let i = 0; i < n; i++) {
    let pred = 0;
    for (let j = 0; j < p; j++) pred += X[i][j] * beta[j];
    const r = y[i] - pred;
    resids.push(r);
    ssRes += W[i] * r * r;
  }

  const invXtWX = invertMatrix(XtWX, p);
  if (!invXtWX) return null;

  const sigmaHat2 = ssRes / (n - p);
  const seBeta = [];
  for (let j = 0; j < p; j++) {
    seBeta.push(Math.sqrt(Math.max(invXtWX[j][j] * sigmaHat2, 1e-12)));
  }

  const results = [{ name: "intercept", coeff: beta[0], se: seBeta[0], t: beta[0] / seBeta[0] }];
  for (let j = 0; j < covariateNames.length; j++) {
    results.push({
      name: covariateNames[j],
      coeff: beta[j + 1],
      se: seBeta[j + 1],
      t: beta[j + 1] / seBeta[j + 1]
    });
  }

  let ssTot = 0;
  const yMean = y.reduce((a, b) => a + b, 0) / n;
  for (let i = 0; i < n; i++) ssTot += W[i] * (y[i] - yMean) ** 2;
  const R2 = 1 - ssRes / ssTot;

  return { coeffs: results, R2, n, p, sigmaHat2, resids };
}

function solveLinear(A, b, n) {
  const aug = [];
  for (let i = 0; i < n; i++) {
    aug.push([...A[i], b[i]]);
  }
  for (let col = 0; col < n; col++) {
    let maxRow = col, maxVal = Math.abs(aug[col][col]);
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(aug[row][col]) > maxVal) {
        maxVal = Math.abs(aug[row][col]);
        maxRow = row;
      }
    }
    if (maxVal < 1e-15) return null;
    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
    const pivot = aug[col][col];
    for (let j = col; j <= n; j++) aug[col][j] /= pivot;
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const factor = aug[row][col];
      for (let j = col; j <= n; j++) aug[row][j] -= factor * aug[col][j];
    }
  }
  return aug.map(row => row[n]);
}

function invertMatrix(A, n) {
  const aug = [];
  for (let i = 0; i < n; i++) {
    const row = [...A[i]];
    for (let j = 0; j < n; j++) row.push(i === j ? 1 : 0);
    aug.push(row);
  }
  for (let col = 0; col < n; col++) {
    let maxRow = col, maxVal = Math.abs(aug[col][col]);
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(aug[row][col]) > maxVal) {
        maxVal = Math.abs(aug[row][col]);
        maxRow = row;
      }
    }
    if (maxVal < 1e-15) return null;
    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
    const pivot = aug[col][col];
    for (let j = 0; j < 2 * n; j++) aug[col][j] /= pivot;
    for (let row = 0; row < n; row++) {
      if (row === col) continue;
      const factor = aug[row][col];
      for (let j = 0; j < 2 * n; j++) aug[row][j] -= factor * aug[col][j];
    }
  }
  return aug.map(row => {
    const inv = new Float64Array(n);
    for (let j = 0; j < n; j++) inv[j] = row[n + j];
    return inv;
  });
}

const tau2_baseline = rAll ? (rAll.tau * rAll.tau) : 0.08;

const covNames = ["isPrecise", "logVmax", "inc/90", "logN", "gRange", "isQ1"];
const covFns = [
  g => isPrecise(g) ? 1 : 0,
  g => Math.log10(g.Vmax),
  g => g.inc / 90,
  g => Math.log10(g.n),
  g => g.gRange,
  g => g.Q_sparc === 1 ? 1 : 0
];

const fullReg = wlsMetaReg(goldI45, covNames, covFns, tau2_baseline);

if (fullReg) {
  log("  Full model (n=" + fullReg.n + ", p=" + fullReg.p + ", R2=" + fullReg.R2.toFixed(4) + "):");
  log("");
  log("  " + pad("Variable", 16) + padr("Coeff", 10) + padr("SE", 10) + padr("t", 8) + "  Sig?");
  log("  " + "\u2500".repeat(55));
  for (const c of fullReg.coeffs) {
    const sig = Math.abs(c.t) >= 2.0 ? " **" : Math.abs(c.t) >= 1.5 ? " *" : "";
    log("  " + pad(c.name, 16) + padr(c.coeff.toFixed(4), 10) +
      padr(c.se.toFixed(4), 10) + padr(c.t.toFixed(2), 8) + sig);
  }
  log("");
  log("  ** = |t| >= 2.0 (significant at ~5%)");
  log("  *  = |t| >= 1.5 (marginal)");
}

log("");
sep();
log("  STEP 3: REDUCED MODELS (sequential variable elimination)");
sep();
log("");

const models = [
  {
    label: "isPrecise only",
    names: ["isPrecise"],
    fns: [g => isPrecise(g) ? 1 : 0]
  },
  {
    label: "isPrecise + logVmax",
    names: ["isPrecise", "logVmax"],
    fns: [g => isPrecise(g) ? 1 : 0, g => Math.log10(g.Vmax)]
  },
  {
    label: "isPrecise + logVmax + inc",
    names: ["isPrecise", "logVmax", "inc/90"],
    fns: [g => isPrecise(g) ? 1 : 0, g => Math.log10(g.Vmax), g => g.inc / 90]
  },
  {
    label: "isPrecise + logVmax + inc + logN",
    names: ["isPrecise", "logVmax", "inc/90", "logN"],
    fns: [g => isPrecise(g) ? 1 : 0, g => Math.log10(g.Vmax), g => g.inc / 90, g => Math.log10(g.n)]
  },
  {
    label: "isPrecise + logVmax + inc + logN + gRange",
    names: ["isPrecise", "logVmax", "inc/90", "logN", "gRange"],
    fns: [g => isPrecise(g) ? 1 : 0, g => Math.log10(g.Vmax), g => g.inc / 90, g => Math.log10(g.n), g => g.gRange]
  }
];

log("  " + pad("Model", 50) + padr("b(prec)", 9) + padr("SE", 9) + padr("t", 8) + padr("R2", 8));
log("  " + "\u2500".repeat(84));

const modelResults = [];
for (const m of models) {
  const reg = wlsMetaReg(goldI45, m.names, m.fns, tau2_baseline);
  if (!reg) continue;
  const precCoeff = reg.coeffs.find(c => c.name === "isPrecise");
  modelResults.push({ label: m.label, reg, precCoeff });
  log("  " + pad(m.label, 50) +
    padr(precCoeff.coeff.toFixed(4), 9) +
    padr(precCoeff.se.toFixed(4), 9) +
    padr(precCoeff.t.toFixed(2), 8) +
    padr(reg.R2.toFixed(4), 8));
}

if (fullReg) {
  const precCoeff = fullReg.coeffs.find(c => c.name === "isPrecise");
  log("  " + pad("FULL (+ gRange + isQ1)", 50) +
    padr(precCoeff.coeff.toFixed(4), 9) +
    padr(precCoeff.se.toFixed(4), 9) +
    padr(precCoeff.t.toFixed(2), 8) +
    padr(fullReg.R2.toFixed(4), 8));
}

log("");
sep();
log("  STEP 4: THREE-WAY DISTANCE SPLIT (precise / UMa / Hubble)");
sep();
log("");

const covNames3 = ["isPrecise", "isUMa", "logVmax", "inc/90", "logN"];
const covFns3 = [
  g => isPrecise(g) ? 1 : 0,
  g => g.fD === 4 ? 1 : 0,
  g => Math.log10(g.Vmax),
  g => g.inc / 90,
  g => Math.log10(g.n)
];

const reg3way = wlsMetaReg(goldI45, covNames3, covFns3, tau2_baseline);
if (reg3way) {
  log("  Three-way model (Hubble = reference):");
  log("");
  for (const c of reg3way.coeffs) {
    const sig = Math.abs(c.t) >= 2.0 ? " **" : Math.abs(c.t) >= 1.5 ? " *" : "";
    log("  " + pad(c.name, 16) + padr(c.coeff.toFixed(4), 10) +
      padr(c.se.toFixed(4), 10) + padr(c.t.toFixed(2), 8) + sig);
  }
  log("");
  const precC = reg3way.coeffs.find(c => c.name === "isPrecise");
  const umaC = reg3way.coeffs.find(c => c.name === "isUMa");
  if (precC && umaC) {
    log("  Precise vs Hubble: +" + precC.coeff.toFixed(4) + " dex (t=" + precC.t.toFixed(2) + ")");
    log("  UMa vs Hubble:     +" + umaC.coeff.toFixed(4) + " dex (t=" + umaC.t.toFixed(2) + ")");
    log("  Precise vs UMa:    +" + (precC.coeff - umaC.coeff).toFixed(4) + " dex");
  }
}

log("");
sep();
log("  STEP 5: PROPENSITY SCORE MATCHING");
log("  Match each precise-distance galaxy to closest Hubble-flow galaxy");
log("  in (logVmax, inc, logN, gRange) space, then compare logA0");
sep();
log("");

function propensityMatch(treated, control) {
  const normalize = (arr, fn) => {
    const vals = arr.map(fn);
    const mu = vals.reduce((a, b) => a + b, 0) / vals.length;
    const sd = Math.sqrt(vals.reduce((a, v) => a + (v - mu) ** 2, 0) / vals.length) || 1;
    return { mu, sd };
  };

  const allGals = [...treated, ...control];
  const normVmax = normalize(allGals, g => Math.log10(g.Vmax));
  const normInc = normalize(allGals, g => g.inc);
  const normN = normalize(allGals, g => Math.log10(g.n));
  const normGR = normalize(allGals, g => g.gRange);

  const feat = g => [
    (Math.log10(g.Vmax) - normVmax.mu) / normVmax.sd,
    (g.inc - normInc.mu) / normInc.sd,
    (Math.log10(g.n) - normN.mu) / normN.sd,
    (g.gRange - normGR.mu) / normGR.sd
  ];

  const dist = (a, b) => {
    const fa = feat(a), fb = feat(b);
    return Math.sqrt(fa.reduce((s, v, i) => s + (v - fb[i]) ** 2, 0));
  };

  const used = new Set();
  const pairs = [];
  for (const t of treated) {
    let bestD = Infinity, bestC = null;
    for (const c of control) {
      if (used.has(c.name)) continue;
      const d = dist(t, c);
      if (d < bestD) { bestD = d; bestC = c; }
    }
    if (bestC && bestD < 3.0) {
      used.add(bestC.name);
      pairs.push({ treated: t, control: bestC, distance: bestD });
    }
  }
  return pairs;
}

const pairs = propensityMatch(preciseGals, hubbleGals);
log("  Matched pairs: " + pairs.length + " / " + preciseGals.length + " precise galaxies");
log("");
log("  " + pad("Precise galaxy", 18) + pad("Matched HF galaxy", 18) +
  padr("logA0_p", 8) + padr("logA0_h", 8) + padr("diff", 8) + padr("matchD", 8));
log("  " + "\u2500".repeat(68));

let sumDiff = 0, sumDiff2 = 0;
for (const p of pairs) {
  const diff = p.treated.logA0 - p.control.logA0;
  sumDiff += diff;
  sumDiff2 += diff * diff;
  log("  " + pad(p.treated.name, 18) + pad(p.control.name, 18) +
    padr(p.treated.logA0.toFixed(3), 8) + padr(p.control.logA0.toFixed(3), 8) +
    padr((diff >= 0 ? "+" : "") + diff.toFixed(3), 8) +
    padr(p.distance.toFixed(2), 8));
}

if (pairs.length > 0) {
  const meanDiff = sumDiff / pairs.length;
  const varDiff = sumDiff2 / pairs.length - meanDiff * meanDiff;
  const seDiff = Math.sqrt(varDiff / pairs.length);
  const tDiff = meanDiff / seDiff;

  log("");
  log("  Mean paired difference: " + (meanDiff >= 0 ? "+" : "") + meanDiff.toFixed(4) + " dex");
  log("  SE of mean difference:  " + seDiff.toFixed(4));
  log("  Paired t-statistic:     " + tDiff.toFixed(2));
  log("  Significance:           " + (Math.abs(tDiff) >= 2.0 ? "YES" : "NO") +
    " (|t| " + (Math.abs(tDiff) >= 2.0 ? ">=" : "<") + " 2.0)");
  log("");
  log("  Effect size: " + (Math.abs(meanDiff) < 0.05 ? "NEGLIGIBLE" :
    Math.abs(meanDiff) < 0.10 ? "SMALL" :
    Math.abs(meanDiff) < 0.20 ? "MODERATE" : "LARGE") +
    " (" + Math.abs(meanDiff).toFixed(3) + " dex)");
}

log("");
sep();
log("  STEP 6: PERMUTATION TEST FOR isPrecise COEFFICIENT");
log("  Randomly shuffle distance-method labels 10000 times,");
log("  refit regression, measure b(isPrecise) each time.");
sep();
log("");

const NPERMS = 10000;
const obsCoeff = fullReg ? fullReg.coeffs.find(c => c.name === "isPrecise").coeff : 0;
let nMoreExtreme = 0;

const fDarray = goldI45.map(g => g.fD);

for (let perm = 0; perm < NPERMS; perm++) {
  const shuffled = [...fDarray];
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }

  const permGals = goldI45.map((g, idx) => ({
    ...g,
    fD_perm: shuffled[idx]
  }));

  const permCovFns = [
    g => (g.fD_perm === 2 || g.fD_perm === 3 || g.fD_perm === 5) ? 1 : 0,
    g => Math.log10(g.Vmax),
    g => g.inc / 90,
    g => Math.log10(g.n),
    g => g.gRange,
    g => g.Q_sparc === 1 ? 1 : 0
  ];

  const permReg = wlsMetaReg(permGals, covNames, permCovFns, tau2_baseline);
  if (permReg) {
    const permCoeff = permReg.coeffs.find(c => c.name === "isPrecise").coeff;
    if (Math.abs(permCoeff) >= Math.abs(obsCoeff)) nMoreExtreme++;
  }
}

const permP = nMoreExtreme / NPERMS;
log("  Observed b(isPrecise) = " + obsCoeff.toFixed(4));
log("  Permutations: " + NPERMS);
log("  More extreme: " + nMoreExtreme);
log("  Permutation p-value: " + permP.toFixed(4));
log("  Significance: " + (permP < 0.05 ? "YES (p < 0.05)" : "NO (p >= 0.05)"));

log("");
sep();
log("  STEP 7: LEAVE-ONE-OUT SENSITIVITY");
log("  Drop each galaxy, refit full model, track b(isPrecise)");
sep();
log("");

const looCoeffs = [];
let minCoeff = Infinity, maxCoeff = -Infinity, minName = "", maxName = "";

for (let i = 0; i < goldI45.length; i++) {
  const subset = goldI45.filter((_, idx) => idx !== i);
  const reg = wlsMetaReg(subset, covNames, covFns, tau2_baseline);
  if (!reg) continue;
  const c = reg.coeffs.find(c => c.name === "isPrecise");
  looCoeffs.push({ dropped: goldI45[i].name, coeff: c.coeff, t: c.t });
  if (c.coeff < minCoeff) { minCoeff = c.coeff; minName = goldI45[i].name; }
  if (c.coeff > maxCoeff) { maxCoeff = c.coeff; maxName = goldI45[i].name; }
}

const looMean = looCoeffs.reduce((a, c) => a + c.coeff, 0) / looCoeffs.length;
const looSD = Math.sqrt(looCoeffs.reduce((a, c) => a + (c.coeff - looMean) ** 2, 0) / looCoeffs.length);

log("  Full-model b(isPrecise) = " + obsCoeff.toFixed(4));
log("  LOO mean:  " + looMean.toFixed(4));
log("  LOO SD:    " + looSD.toFixed(4));
log("  LOO range: [" + minCoeff.toFixed(4) + ", " + maxCoeff.toFixed(4) + "]");
log("  Most influential drops: " + minName + " (min), " + maxName + " (max)");

const top5 = [...looCoeffs].sort((a, b) => Math.abs(b.coeff - obsCoeff) - Math.abs(a.coeff - obsCoeff)).slice(0, 5);
log("");
log("  Top 5 most influential galaxies:");
for (const c of top5) {
  log("    Drop " + pad(c.dropped, 18) + " => b=" + c.coeff.toFixed(4) + " (t=" + c.t.toFixed(2) + ")");
}

log("");
sep();
log("  PHASE 8 SUMMARY");
sep();
log("");

const results = {
  version: VERSION,
  timestamp: TIMESTAMP,
  description: "Phase 8: Distance-method matched meta-regression",
  unadjusted: {
    precise: rPrec,
    hubbleFlow: rHub,
    UMaCluster: rOther,
    all: rAll
  },
  fullRegression: fullReg ? {
    coefficients: fullReg.coeffs,
    R2: fullReg.R2,
    n: fullReg.n,
    p: fullReg.p
  } : null,
  propensityMatching: pairs.length > 0 ? {
    nPairs: pairs.length,
    meanDiff: +(sumDiff / pairs.length).toFixed(4),
    seDiff: +(Math.sqrt((sumDiff2 / pairs.length - (sumDiff / pairs.length) ** 2) / pairs.length)).toFixed(4),
    tStat: +((sumDiff / pairs.length) / Math.sqrt((sumDiff2 / pairs.length - (sumDiff / pairs.length) ** 2) / pairs.length)).toFixed(2),
    pairs: pairs.map(p => ({ precise: p.treated.name, hubble: p.control.name, diff: +(p.treated.logA0 - p.control.logA0).toFixed(4) }))
  } : null,
  permutationTest: {
    observedCoeff: +obsCoeff.toFixed(4),
    nPermutations: NPERMS,
    nMoreExtreme: nMoreExtreme,
    pValue: +permP.toFixed(4)
  },
  leaveOneOut: {
    mean: +looMean.toFixed(4),
    sd: +looSD.toFixed(4),
    range: [+minCoeff.toFixed(4), +maxCoeff.toFixed(4)],
    mostInfluential: top5.map(c => ({ dropped: c.dropped, coeff: +c.coeff.toFixed(4), t: +c.t.toFixed(2) }))
  }
};

if (fullReg) {
  const precC = fullReg.coeffs.find(c => c.name === "isPrecise");
  log("  ╔══════════════════════════════════════════════════════════════════╗");
  log("  ║  DISTANCE-METHOD META-REGRESSION RESULT                         ║");
  log("  ║                                                                  ║");
  log("  ║  Full model b(isPrecise) = " + pad(precC.coeff.toFixed(4) + " dex (t=" + precC.t.toFixed(2) + ")", 33) + "║");
  log("  ║  Permutation p-value = " + pad(permP.toFixed(4), 37) + "║");
  if (pairs.length > 0) {
    const mDiff = sumDiff / pairs.length;
    log("  ║  Propensity-matched diff = " + pad((mDiff >= 0 ? "+" : "") + mDiff.toFixed(4) + " dex", 33) + "║");
  }
  log("  ║  LOO range: [" + pad(minCoeff.toFixed(4) + ", " + maxCoeff.toFixed(4) + "]", 47) + "║");
  log("  ║                                                                  ║");

  if (permP < 0.05 && Math.abs(precC.t) >= 2.0) {
    log("  ║  VERDICT: Distance method effect is SIGNIFICANT.               ║");
    log("  ║  Precise distances give higher a0 even after controlling       ║");
    log("  ║  for Vmax, inclination, sample size, and dynamic range.        ║");
    results.verdict = "SIGNIFICANT";
  } else if (permP < 0.10 || Math.abs(precC.t) >= 1.5) {
    log("  ║  VERDICT: Distance method effect is SUGGESTIVE but not         ║");
    log("  ║  conclusive. The signal exists but does not reach              ║");
    log("  ║  conventional significance thresholds.                         ║");
    results.verdict = "SUGGESTIVE";
  } else {
    log("  ║  VERDICT: Distance method effect is NOT SIGNIFICANT            ║");
    log("  ║  after controlling for covariates. The raw split may be        ║");
    log("  ║  confounded by other galaxy properties.                        ║");
    results.verdict = "NOT SIGNIFICANT";
  }
  log("  ╚══════════════════════════════════════════════════════════════════╝");
}

const outPath = path.join(__dirname, '..', 'public', 'phase8-metareg-results.json');
fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
log("");
log("  Results saved to: public/phase8-metareg-results.json");
log("");
log("================================================================================");
log("  PHASE 8 COMPLETE");
log("================================================================================");
