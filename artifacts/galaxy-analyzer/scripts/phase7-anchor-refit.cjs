const fs = require('fs');
const path = require('path');

const VERSION = "7.0.0";
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

function median(arr) {
  if (!arr.length) return NaN;
  const s = [...arr].sort((a, b) => a - b);
  const mid = Math.floor(s.length / 2);
  return s.length % 2 !== 0 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
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
  const a0 = Math.pow(10, muRE);
  const ratio = a0 / CH0_2PI_GAL;
  const lo1sig = Math.pow(10, muRE - Math.sqrt(seRE * seRE + tau2));
  const hi1sig = Math.pow(10, muRE + Math.sqrt(seRE * seRE + tau2));
  return {
    label, nGal: k, a0: Math.round(a0),
    logA0: +muRE.toFixed(4), se: +seRE.toFixed(4),
    tau: +tau.toFixed(4), I2: +I2.toFixed(1),
    ratio: +ratio.toFixed(4),
    lo1sig: Math.round(lo1sig), hi1sig: Math.round(hi1sig)
  };
}

function pearsonR(xs, ys) {
  const n = xs.length;
  if (n < 5) return { r: NaN, p: NaN };
  const mx = xs.reduce((a, b) => a + b, 0) / n;
  const my = ys.reduce((a, b) => a + b, 0) / n;
  let sxx = 0, syy = 0, sxy = 0;
  for (let i = 0; i < n; i++) {
    sxx += (xs[i] - mx) ** 2;
    syy += (ys[i] - my) ** 2;
    sxy += (xs[i] - mx) * (ys[i] - my);
  }
  const r = sxy / Math.sqrt(sxx * syy);
  const t = r * Math.sqrt((n - 2) / (1 - r * r));
  return { r: +r.toFixed(4), t: +t.toFixed(2), n };
}

function printResult(r, indent) {
  indent = indent || "    ";
  if (!r) { log(indent + "(insufficient data)"); return; }
  log(indent + pad(r.label, 35) + " n=" + padr(r.nGal, 3) +
    " a0=" + padr(r.a0, 5) + " tau=" + padr(r.tau, 6) +
    " I2=" + padr(r.I2 + "%", 6) +
    " ratio=" + padr(r.ratio, 6) + " 68%=[" + r.lo1sig + "," + r.hi1sig + "]");
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
    Rdisk: g.Rdisk || sparc.Rdisk || 1.0,
    Q_kin: g.Q_kin, Q_sparc: sparc.Q || null,
    T: sparc.T !== undefined ? sparc.T : null,
    eD: sparc.eD || null, eInc: sparc.eInc || null,
    fD: sparc.fD || null, points: []
  };
}

for (const p of rarData.rarScatter) {
  if (galaxyMap[p.name]) {
    galaxyMap[p.name].points.push({
      logGbar: p.log_g_bar, logGobs: p.log_g_obs
    });
  }
}

for (const p of tsData.plotPoints) {
  if (p.g.startsWith('LT_')) {
    if (!galaxyMap[p.g]) {
      galaxyMap[p.g] = {
        name: p.g, inc: 60, distance: 5.0, L36: 0, MHI: 0,
        Vmax: 50, fgas: 0.5, Rdisk: 1.0, Q_kin: 0.7,
        Q_sparc: null, T: null, eD: null, eInc: null, fD: null,
        points: [], isLT: true
      };
    }
    galaxyMap[p.g].points.push({ logGbar: p.x, logGobs: p.x + p.y });
  }
}

const allGalaxies = Object.values(galaxyMap).filter(g => g.points.length >= 5);

log("================================================================================");
log("  PHASE 7: ANCHOR-SAMPLE REFIT");
log("================================================================================");
log("  Version: " + VERSION);
log("  Date: " + TIMESTAMP);
log("  Purpose: Independent test using highest-quality anchor galaxies only.");
log("  Key question: Does tau drop or remain high in a clean anchor sample?");
log("  If tau drops: remaining heterogeneity in GOLD+i45 was systematics.");
log("  If tau stays high: heterogeneity is real, not driven by data quality.");
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

log("");
sep();
log("  REFERENCE: GOLD+i45 BASELINE");
sep();
log("");
const baseline = hierarchicalDL(goldI45, "GOLD+i45 (full)");
printResult(baseline);
const baseCorrs = {
  distance: pearsonR(goldI45.map(g => Math.log10(g.distance)), goldI45.map(g => g.logA0)),
  inc: pearsonR(goldI45.map(g => g.inc), goldI45.map(g => g.logA0)),
  Vmax: pearsonR(goldI45.map(g => Math.log10(g.Vmax)), goldI45.map(g => g.logA0))
};
log("    r(logA0, logD) = " + baseCorrs.distance.r + " (t=" + baseCorrs.distance.t + ")");
log("    r(logA0, inc)  = " + baseCorrs.inc.r + " (t=" + baseCorrs.inc.t + ")");
log("    r(logA0, logV) = " + baseCorrs.Vmax.r + " (t=" + baseCorrs.Vmax.t + ")");

log("");
sep();
log("  TIER 1: ANCHOR SAMPLE (strictest)");
log("  Criteria: precise distance (TRGB/Cepheid/SN) + Q=1 + inc>=60");
sep();
log("");

const tier1 = goldI45.filter(g => {
  const isPrecise = g.fD === 2 || g.fD === 3 || g.fD === 5;
  return isPrecise && g.Q_sparc === 1 && g.inc >= 60;
});

log("  Galaxies: " + tier1.length);
for (const g of tier1.sort((a, b) => b.Vmax - a.Vmax)) {
  log("    " + pad(g.name, 20) +
    " Vmax=" + padr(Math.round(g.Vmax), 4) +
    " inc=" + padr(Math.round(g.inc), 3) +
    " D=" + padr(g.distance.toFixed(1), 5) +
    " fD=" + g.fD +
    " n=" + padr(g.n, 3) +
    " logA0=" + g.logA0.toFixed(3) +
    " a0=" + Math.round(g.a0));
}
log("");
const tier1R = hierarchicalDL(tier1, "ANCHOR (strict)");
printResult(tier1R);
if (tier1.length >= 5) {
  const t1corrs = {
    distance: pearsonR(tier1.map(g => Math.log10(g.distance)), tier1.map(g => g.logA0)),
    inc: pearsonR(tier1.map(g => g.inc), tier1.map(g => g.logA0)),
    Vmax: pearsonR(tier1.map(g => Math.log10(g.Vmax)), tier1.map(g => g.logA0))
  };
  log("    r(logA0, logD) = " + t1corrs.distance.r + " (t=" + t1corrs.distance.t + ")");
  log("    r(logA0, inc)  = " + t1corrs.inc.r + " (t=" + t1corrs.inc.t + ")");
  log("    r(logA0, logV) = " + t1corrs.Vmax.r + " (t=" + t1corrs.Vmax.t + ")");
}

log("");
sep();
log("  TIER 2: ANCHOR + UMa CLUSTER (relax to include fD=4)");
log("  Criteria: non-Hubble-flow (fD=2,3,4,5) + Q=1 + inc>=60");
sep();
log("");

const tier2 = goldI45.filter(g => {
  return g.fD !== 1 && g.Q_sparc === 1 && g.inc >= 60;
});

log("  Galaxies: " + tier2.length);
for (const g of tier2.sort((a, b) => b.Vmax - a.Vmax)) {
  log("    " + pad(g.name, 20) +
    " Vmax=" + padr(Math.round(g.Vmax), 4) +
    " inc=" + padr(Math.round(g.inc), 3) +
    " D=" + padr(g.distance.toFixed(1), 5) +
    " fD=" + g.fD +
    " n=" + padr(g.n, 3) +
    " logA0=" + g.logA0.toFixed(3));
}
log("");
const tier2R = hierarchicalDL(tier2, "ANCHOR+UMa (Q=1, fD!=1, i>=60)");
printResult(tier2R);
if (tier2.length >= 5) {
  const t2corrs = {
    distance: pearsonR(tier2.map(g => Math.log10(g.distance)), tier2.map(g => g.logA0)),
    inc: pearsonR(tier2.map(g => g.inc), tier2.map(g => g.logA0))
  };
  log("    r(logA0, logD) = " + t2corrs.distance.r + " (t=" + t2corrs.distance.t + ")");
  log("    r(logA0, inc)  = " + t2corrs.inc.r + " (t=" + t2corrs.inc.t + ")");
}

log("");
sep();
log("  TIER 3: RELAX INCLINATION (anchor distances + Q=1 + inc>=45)");
log("  Criteria: precise distance (fD=2,3,5) + Q=1 + inc>=45");
sep();
log("");

const tier3 = goldI45.filter(g => {
  const isPrecise = g.fD === 2 || g.fD === 3 || g.fD === 5;
  return isPrecise && g.Q_sparc === 1;
});

log("  Galaxies: " + tier3.length);
for (const g of tier3.sort((a, b) => b.Vmax - a.Vmax)) {
  log("    " + pad(g.name, 20) +
    " Vmax=" + padr(Math.round(g.Vmax), 4) +
    " inc=" + padr(Math.round(g.inc), 3) +
    " D=" + padr(g.distance.toFixed(1), 5) +
    " n=" + padr(g.n, 3) +
    " logA0=" + g.logA0.toFixed(3));
}
log("");
const tier3R = hierarchicalDL(tier3, "ANCHOR (precise+Q1, inc>=45)");
printResult(tier3R);
if (tier3.length >= 5) {
  const t3corrs = {
    distance: pearsonR(tier3.map(g => Math.log10(g.distance)), tier3.map(g => g.logA0)),
    inc: pearsonR(tier3.map(g => g.inc), tier3.map(g => g.logA0))
  };
  log("    r(logA0, logD) = " + t3corrs.distance.r + " (t=" + t3corrs.distance.t + ")");
  log("    r(logA0, inc)  = " + t3corrs.inc.r + " (t=" + t3corrs.inc.t + ")");
}

log("");
sep();
log("  TIER 4: Q=1 ONLY (any distance method, inc>=45)");
log("  Criteria: Q=1 + GOLD+i45 base cuts");
sep();
log("");

const tier4 = goldI45.filter(g => g.Q_sparc === 1);

log("  Galaxies: " + tier4.length);
const tier4R = hierarchicalDL(tier4, "Q=1 only (all dist methods)");
printResult(tier4R);
if (tier4.length >= 5) {
  const t4corrs = {
    distance: pearsonR(tier4.map(g => Math.log10(g.distance)), tier4.map(g => g.logA0)),
    inc: pearsonR(tier4.map(g => g.inc), tier4.map(g => g.logA0))
  };
  log("    r(logA0, logD) = " + t4corrs.distance.r + " (t=" + t4corrs.distance.t + ")");
  log("    r(logA0, inc)  = " + t4corrs.inc.r + " (t=" + t4corrs.inc.t + ")");
}

log("");
sep();
log("  TIER 5: PRECISE DISTANCE ONLY (any Q, inc>=45)");
log("  Criteria: fD=2,3,5 + GOLD+i45 base cuts");
sep();
log("");

const tier5 = goldI45.filter(g => g.fD === 2 || g.fD === 3 || g.fD === 5);

log("  Galaxies: " + tier5.length);
const tier5R = hierarchicalDL(tier5, "Precise dist only (any Q)");
printResult(tier5R);
if (tier5.length >= 5) {
  const t5corrs = {
    distance: pearsonR(tier5.map(g => Math.log10(g.distance)), tier5.map(g => g.logA0)),
    inc: pearsonR(tier5.map(g => g.inc), tier5.map(g => g.logA0))
  };
  log("    r(logA0, logD) = " + t5corrs.distance.r + " (t=" + t5corrs.distance.t + ")");
  log("    r(logA0, inc)  = " + t5corrs.inc.r + " (t=" + t5corrs.inc.t + ")");
}

log("");
sep();
log("  COMPARISON TABLE");
sep();
log("");
log("  " + pad("Sample", 40) + padr("n", 4) + padr("a0", 6) +
  padr("tau", 7) + padr("I2", 7) + padr("ratio", 7));
log("  " + "\u2500".repeat(71));

const allResults = [
  baseline, tier1R, tier2R, tier3R, tier4R, tier5R
].filter(r => r);

for (const r of allResults) {
  log("  " + pad(r.label, 40) + padr(r.nGal, 4) + padr(r.a0, 6) +
    padr(r.tau, 7) + padr(r.I2 + "%", 7) + padr(r.ratio, 7));
}

log("");
sep();
log("  KEY DIAGNOSTIC: TAU BEHAVIOR");
sep();
log("");

if (tier1R) {
  const tauChange = tier1R.tau - baseline.tau;
  const tauPct = (tauChange / baseline.tau * 100);
  log("  Baseline tau (GOLD+i45):  " + baseline.tau);
  log("  Anchor tau (strict):      " + tier1R.tau);
  log("  Change:                   " + (tauChange >= 0 ? "+" : "") + tauChange.toFixed(4) +
    " (" + (tauPct >= 0 ? "+" : "") + tauPct.toFixed(1) + "%)");
  log("");

  if (Math.abs(tauPct) < 15) {
    log("  INTERPRETATION: tau is STABLE across quality tiers.");
    log("  The between-galaxy heterogeneity (tau ~ 0.29 dex) is NOT driven");
    log("  by distance measurement errors or data quality differences.");
    log("  It appears to be an INTRINSIC property of the galaxy population,");
    log("  or reflects astrophysical diversity (disk structure, gas content,");
    log("  star formation history) rather than observational systematics.");
  } else if (tauPct < -15) {
    log("  INTERPRETATION: tau DROPS in the anchor sample.");
    log("  This suggests that some of the GOLD+i45 heterogeneity was");
    log("  driven by distance/quality systematics. The anchor sample");
    log("  provides a cleaner measurement of a0.");
  } else {
    log("  INTERPRETATION: tau INCREASES in the anchor sample.");
    log("  This could indicate that the strict cuts select a more");
    log("  heterogeneous subsample, or that the smaller sample size");
    log("  inflates the DL tau estimate.");
  }
}

log("");
sep();
log("  KEY DIAGNOSTIC: DISTANCE CORRELATION");
sep();
log("");

const baseDistR = pearsonR(goldI45.map(g => Math.log10(g.distance)), goldI45.map(g => g.logA0));
log("  GOLD+i45 r(logA0, logD) = " + baseDistR.r + " (n=" + baseDistR.n + ", t=" + baseDistR.t + ")");

if (tier1.length >= 5) {
  const anchorDistR = pearsonR(tier1.map(g => Math.log10(g.distance)), tier1.map(g => g.logA0));
  log("  Anchor   r(logA0, logD) = " + anchorDistR.r + " (n=" + anchorDistR.n + ", t=" + anchorDistR.t + ")");
  log("");

  if (Math.abs(anchorDistR.r) < 0.20) {
    log("  INTERPRETATION: Distance correlation DISAPPEARS in anchor sample.");
    log("  The r = " + baseDistR.r + " seen in GOLD+i45 was likely driven by");
    log("  Hubble-flow distance errors. With precise distances only,");
    log("  there is no significant distance-a0 correlation.");
  } else if (Math.abs(anchorDistR.r) > Math.abs(baseDistR.r)) {
    log("  INTERPRETATION: Distance correlation PERSISTS or STRENGTHENS.");
    log("  This would suggest a genuine distance-dependent effect,");
    log("  not merely Hubble-flow errors.");
  } else {
    log("  INTERPRETATION: Distance correlation WEAKENED but not eliminated.");
    log("  Mixed evidence — both distance errors and possible real effects.");
  }
}

log("");
sep();
log("  KEY DIAGNOSTIC: a0 VALUE SHIFT");
sep();
log("");

if (tier1R) {
  const a0shift = tier1R.logA0 - baseline.logA0;
  const a0pct = (tier1R.a0 / baseline.a0 - 1) * 100;
  log("  GOLD+i45 a0:  " + baseline.a0 + " (ratio=" + baseline.ratio + ")");
  log("  Anchor a0:    " + tier1R.a0 + " (ratio=" + tier1R.ratio + ")");
  log("  Shift:        " + (a0shift >= 0 ? "+" : "") + a0shift.toFixed(4) + " dex (" +
    (a0pct >= 0 ? "+" : "") + a0pct.toFixed(1) + "%)");
  log("");

  if (Math.abs(a0shift) < 0.05) {
    log("  a0 is STABLE: anchor and full GOLD+i45 agree within 0.05 dex.");
    log("  The headline value of " + baseline.a0 + " is not biased by");
    log("  low-quality distances or data.");
  } else if (a0shift > 0.05) {
    log("  a0 INCREASES in anchor sample. Precise-distance galaxies");
    log("  give higher a0, consistent with the Phase 6 distance split.");
    log("  The headline value may be biased LOW by Hubble-flow errors.");
  } else {
    log("  a0 DECREASES in anchor sample. This would suggest the");
    log("  headline value may be biased HIGH by selection effects.");
  }
}

log("");
sep();
log("  PHASE 7 VERDICT");
sep();
log("");

const results = {
  version: VERSION,
  timestamp: TIMESTAMP,
  description: "Phase 7: Anchor-sample refit with highest-quality galaxies",
  baseline: baseline,
  tiers: {
    tier1_strict: tier1R,
    tier1_galaxies: tier1.map(g => g.name),
    tier2_withUMa: tier2R,
    tier3_relaxInc: tier3R,
    tier4_Q1only: tier4R,
    tier5_precDist: tier5R
  },
  diagnostics: {}
};

if (tier1R) {
  const tauChange = ((tier1R.tau - baseline.tau) / baseline.tau * 100);
  const a0Shift = tier1R.logA0 - baseline.logA0;
  const anchorDistR = tier1.length >= 5 ?
    pearsonR(tier1.map(g => Math.log10(g.distance)), tier1.map(g => g.logA0)) : null;

  results.diagnostics = {
    tauChange_pct: +tauChange.toFixed(1),
    a0Shift_dex: +a0Shift.toFixed(4),
    baselineDistCorr: baseDistR.r,
    anchorDistCorr: anchorDistR ? anchorDistR.r : null,
    tauInterpretation: Math.abs(tauChange) < 15 ? "STABLE" :
      tauChange < -15 ? "DROPS" : "INCREASES",
    a0Interpretation: Math.abs(a0Shift) < 0.05 ? "STABLE" :
      a0Shift > 0.05 ? "INCREASES" : "DECREASES"
  };

  log("  ╔══════════════════════════════════════════════════════════════════╗");
  log("  ║  ANCHOR-SAMPLE REFIT RESULT                                     ║");
  log("  ║                                                                  ║");
  log("  ║  Anchor sample: " + pad(tier1R.nGal + " galaxies (Q=1, precise D, inc>=60)", 39) + "║");
  log("  ║  a0 = " + pad(tier1R.a0 + " (vs baseline " + baseline.a0 + ")", 50) + "║");
  log("  ║  tau = " + pad(tier1R.tau + " (vs baseline " + baseline.tau + ")", 49) + "║");
  log("  ║  ratio = " + pad(tier1R.ratio + " (vs baseline " + baseline.ratio + ")", 46) + "║");
  log("  ║                                                                  ║");
  log("  ║  tau change: " + pad((tauChange >= 0 ? "+" : "") + tauChange.toFixed(1) + "% => " + results.diagnostics.tauInterpretation, 43) + "║");
  log("  ║  a0 shift:   " + pad((a0Shift >= 0 ? "+" : "") + a0Shift.toFixed(3) + " dex => " + results.diagnostics.a0Interpretation, 43) + "║");
  if (anchorDistR) {
    log("  ║  r(a0,D):    " + pad(anchorDistR.r + " (vs baseline " + baseDistR.r + ")", 43) + "║");
  }
  log("  ╚══════════════════════════════════════════════════════════════════╝");
}

const outPath = path.join(__dirname, '..', 'public', 'phase7-anchor-results.json');
fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
log("");
log("  Results saved to: public/phase7-anchor-results.json");
log("");
log("================================================================================");
log("  PHASE 7 COMPLETE");
log("================================================================================");
