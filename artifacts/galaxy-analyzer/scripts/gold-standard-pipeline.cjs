const fs = require('fs');
const path = require('path');

const VERSION = "2.0.0";
const TIMESTAMP = new Date().toISOString();

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const INTERP_SET = ["McGaugh RAR", "Simple MOND", "Standard MOND"];
const UNIT_LABEL = "(km/s)^2/kpc";
const MS2_CONV = 3.241e-14;
const C_MS = 299792458;
const H0_SI = 67.4e3 / 3.086e22;
const CH0_2PI = C_MS * H0_SI / (2 * Math.PI);

const tsData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'transition-scale.json'), 'utf8'));
const rarData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));

const galaxyMeta = {};
for (const g of rarData.perGalaxy) galaxyMeta[g.name] = g;

const perGalaxyFits = {};
for (const g of tsData.perGalaxyA0.galaxies) perGalaxyFits[g.name] = g;

const ALL_POINTS = tsData.plotPoints;

function log(msg) { console.log(msg); }
function sep() { log("─".repeat(72)); }

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}
function simpleMOND(gbar, a0) {
  const y = gbar / a0;
  return gbar * (0.5 + Math.sqrt(0.25 + 1.0 / y));
}
function standardMOND(gbar, a0) {
  const y = gbar / a0;
  return gbar * (1 + Math.sqrt(1 + 4.0 / y)) / 2.0;
}
const INTERP_FUNCS = [mcgaughRAR, simpleMOND, standardMOND];

function rms(vals) {
  if (vals.length < 2) return Infinity;
  const m = vals.reduce((a, b) => a + b, 0) / vals.length;
  return Math.sqrt(vals.reduce((a, v) => a + (v - m) ** 2, 0) / vals.length);
}

function median(arr) {
  if (!arr.length) return NaN;
  const s = [...arr].sort((a, b) => a - b);
  const mid = Math.floor(s.length / 2);
  return s.length % 2 !== 0 ? s[mid] : (s[mid - 1] + s[mid]) / 2;
}

function mad(arr) {
  const med = median(arr);
  return median(arr.map(v => Math.abs(v - med)));
}

function pearsonR(x, y) {
  const n = x.length;
  if (n < 5) return NaN;
  const mx = x.reduce((a, b) => a + b) / n, my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return syy > 0 && sxx > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function fitA0(points, interpFunc) {
  let bestA0 = 3000, bestRMS = Infinity;
  for (let logA = 2.0; logA <= 5.0; logA += 0.01) {
    const a0 = Math.pow(10, logA);
    const r = evalRMS(points, interpFunc, a0);
    if (r < bestRMS) { bestRMS = r; bestA0 = a0; }
  }
  let lo = Math.log10(bestA0) - 0.02, hi = Math.log10(bestA0) + 0.02;
  for (let step = 0; step < 100; step++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    if (evalRMS(points, interpFunc, Math.pow(10, m1)) < evalRMS(points, interpFunc, Math.pow(10, m2))) hi = m2;
    else lo = m1;
  }
  const finalA0 = Math.pow(10, (lo + hi) / 2);
  const finalRMS = evalRMS(points, interpFunc, finalA0);
  const resids = getResids(points, interpFunc, finalA0);
  const meanBias = resids.reduce((a, b) => a + b, 0) / resids.length;
  return { a0: finalA0, logA0: Math.log10(finalA0), rms: finalRMS, bias: meanBias, n: resids.length };
}

function evalRMS(points, interpFunc, a0) {
  const resids = getResids(points, interpFunc, a0);
  return rms(resids);
}

function getResids(points, interpFunc, a0) {
  const resids = [];
  for (const p of points) {
    const gbar = Math.pow(10, p.x);
    const gobs = gbar * Math.pow(10, p.y);
    const pred = interpFunc(gbar, a0);
    if (!isFinite(pred) || pred <= 0) continue;
    resids.push(Math.log10(gobs) - Math.log10(pred));
  }
  return resids;
}

function getPointsForGalaxies(galaxyNames) {
  const shortSet = new Set(galaxyNames.map(n => n.substring(0, 12)));
  return ALL_POINTS.filter(p => shortSet.has(p.g));
}

function getGbarRange(galaxyName) {
  const pts = ALL_POINTS.filter(p => p.g === galaxyName.substring(0, 12));
  if (pts.length < 2) return 0;
  return Math.max(...pts.map(p => p.x)) - Math.min(...pts.map(p => p.x));
}

function blindAudit(sampleName, sampleGals, fitA0Val) {
  const residData = [];
  for (const gname of sampleGals) {
    const pts = ALL_POINTS.filter(p => p.g === gname.substring(0, 12));
    if (pts.length < 3) continue;
    const resids = [];
    for (const p of pts) {
      const gbar = Math.pow(10, p.x);
      const gobs = gbar * Math.pow(10, p.y);
      const pred = mcgaughRAR(gbar, fitA0Val);
      if (!isFinite(pred) || pred <= 0) continue;
      resids.push(Math.log10(gobs) - Math.log10(pred));
    }
    if (resids.length < 2) continue;
    const mean = resids.reduce((a, b) => a + b, 0) / resids.length;
    const meta = galaxyMeta[gname];
    if (!meta) continue;
    const gasFrac = meta.MHI > 0 && meta.L36 > 0 ? meta.MHI / meta.L36 : NaN;
    residData.push({
      name: gname, meanResid: mean,
      logVmax: Math.log10(meta.Vmax),
      logSB: meta.sigma_bar > 0 ? Math.log10(meta.sigma_bar) : NaN,
      logGasFrac: isFinite(gasFrac) && gasFrac > 0 ? Math.log10(gasFrac) : NaN,
      logL: meta.L36 > 0 ? Math.log10(meta.L36) : NaN,
      inc: meta.inc,
      logDist: meta.distance > 0 ? Math.log10(meta.distance) : NaN,
    });
  }
  const valid = residData.filter(g => isFinite(g.logVmax) && isFinite(g.logSB) && isFinite(g.logGasFrac) && isFinite(g.logL));
  const auditVars = [
    { name: "V_max", vals: valid.map(g => [g.meanResid, g.logVmax]) },
    { name: "Surface brightness", vals: valid.map(g => [g.meanResid, g.logSB]) },
    { name: "Gas fraction", vals: valid.map(g => [g.meanResid, g.logGasFrac]) },
    { name: "Luminosity", vals: valid.map(g => [g.meanResid, g.logL]) },
    { name: "Inclination", vals: valid.filter(g => g.inc > 0).map(g => [g.meanResid, g.inc]) },
    { name: "Distance", vals: valid.filter(g => isFinite(g.logDist)).map(g => [g.meanResid, g.logDist]) },
  ];
  const results = [];
  let clean = true;
  for (const av of auditVars) {
    if (av.vals.length < 10) continue;
    const r = pearsonR(av.vals.map(v => v[0]), av.vals.map(v => v[1]));
    const sig = Math.abs(r) > 0.3;
    if (sig) clean = false;
    results.push({ variable: av.name, r: +r.toFixed(3), n: av.vals.length, significant: sig });
  }
  return { sampleName, nGalaxies: valid.length, a0Used: fitA0Val, clean, correlations: results };
}

const allGalNames = Object.keys(galaxyMeta);

const cleanSample = allGalNames.filter(n => {
  const m = galaxyMeta[n]; const range = getGbarRange(n);
  return m.Vmax >= 80 && m.inc >= 30 && m.n >= 5 && range >= 0.5;
});
const goldSample = allGalNames.filter(n => {
  const m = galaxyMeta[n]; const range = getGbarRange(n);
  return m.Vmax >= 50 && m.inc >= 30 && m.n >= 5 && range >= 1.0;
});
const goldTight = allGalNames.filter(n => {
  const m = galaxyMeta[n]; const range = getGbarRange(n);
  return m.Vmax >= 50 && m.inc >= 45 && m.n >= 5 && range >= 1.0;
});
const highMass = allGalNames.filter(n => galaxyMeta[n].Vmax > 150);

log("╔════════════════════════════════════════════════════════════════════════╗");
log("║  GOLD-STANDARD a₀ MEASUREMENT PIPELINE                              ║");
log("║  Version: " + VERSION + "                                                      ║");
log("║  FINAL MEASUREMENT PROTOCOL                                         ║");
log("╚════════════════════════════════════════════════════════════════════════╝\n");

log("LOCKED PARAMETERS:");
log("  Υ_disk = " + UPSILON_DISK + " M☉/L☉");
log("  Υ_bulge = " + UPSILON_BULGE + " M☉/L☉");
log("  Interpolation functions: " + INTERP_SET.join(", "));
log("  Fitting: Golden section search, log(a₀) ∈ [2.0, 5.0], 100 refinement steps");
log("  Residual: log10(g_obs) - log10(g_pred)");
log("  Units: " + UNIT_LABEL);
log("  1 " + UNIT_LABEL + " = " + MS2_CONV + " m/s²");
log("");

log("DATASET:");
log("  Total points: " + ALL_POINTS.length);
log("  Total galaxies: " + allGalNames.length);
log("  Source: transition-scale.json (every-2nd-point subsample of 4123)");
log("");

log("SAMPLE DEFINITIONS:");
log("  ALL:       " + allGalNames.length + " galaxies, " + ALL_POINTS.length + " points (no cuts)");
log("  CLEAN:     " + cleanSample.length + " galaxies, " + getPointsForGalaxies(cleanSample).length + " points");
log("             (V≥80, inc≥30°, n≥5, range≥0.5 dex)");
log("  GOLD:      " + goldSample.length + " galaxies, " + getPointsForGalaxies(goldSample).length + " points");
log("             (V≥50, inc≥30°, n≥5, range≥1.0 dex)");
log("  GOLD+i45:  " + goldTight.length + " galaxies, " + getPointsForGalaxies(goldTight).length + " points");
log("             (V≥50, inc≥45°, n≥5, range≥1.0 dex)  ← NEW v2.0");
log("  HIGH-MASS: " + highMass.length + " galaxies, " + getPointsForGalaxies(highMass).length + " points");
log("             (V>150 km/s, no other cuts)");
log("");

sep();
log("LAYER 1: FIVE-LAYER a₀ MEASUREMENT");
sep();
log("");

const layers = [
  { name: "ALL", points: ALL_POINTS, gals: allGalNames },
  { name: "CLEAN", points: getPointsForGalaxies(cleanSample), gals: cleanSample },
  { name: "GOLD", points: getPointsForGalaxies(goldSample), gals: goldSample },
  { name: "GOLD+i45", points: getPointsForGalaxies(goldTight), gals: goldTight },
  { name: "HIGH-MASS", points: getPointsForGalaxies(highMass), gals: highMass },
];

const layerResults = {};
for (const layer of layers) {
  log("  " + layer.name + " (n=" + layer.points.length + ", " + layer.gals.length + " galaxies):");
  const results = [];
  for (let i = 0; i < INTERP_FUNCS.length; i++) {
    const fit = fitA0(layer.points, INTERP_FUNCS[i]);
    results.push({ func: INTERP_SET[i], ...fit });
    log("    " + INTERP_SET[i].padEnd(18) + "a₀ = " + fit.a0.toFixed(1).padEnd(10) +
      "log = " + fit.logA0.toFixed(4).padEnd(10) +
      "RMS = " + fit.rms.toFixed(4).padEnd(10) +
      "bias = " + fit.bias.toFixed(4));
  }
  const a0s = results.map(r => r.a0);
  const meanA0 = a0s.reduce((a, b) => a + b) / a0s.length;
  const spread = ((Math.max(...a0s) - Math.min(...a0s)) / meanA0 * 100);
  log("    Mean a₀: " + meanA0.toFixed(1) + " | spread: " + spread.toFixed(1) + "% | " +
    (meanA0 * MS2_CONV).toExponential(3) + " m/s² | ratio cH₀/2π: " + (meanA0 * MS2_CONV / CH0_2PI).toFixed(3));
  log("");
  layerResults[layer.name] = {
    results, nGal: layer.gals.length, nPts: layer.points.length,
    meanA0: +meanA0.toFixed(1), spread: +spread.toFixed(1),
    a0_ms2: +(meanA0 * MS2_CONV).toExponential(3),
    ratio_cH0: +(meanA0 * MS2_CONV / CH0_2PI).toFixed(3),
    bestRMS: +Math.min(...results.map(r => r.rms)).toFixed(4)
  };
}

sep();
log("LAYER 2: PER-GALAXY MEDIAN a₀ (McGaugh RAR)");
sep();
log("");

const perGalA0 = [];
for (const gname of allGalNames) {
  const pts = ALL_POINTS.filter(p => p.g === gname.substring(0, 12));
  if (pts.length < 5) continue;
  const range = getGbarRange(gname);
  if (range < 0.3) continue;
  const fit = fitA0(pts, mcgaughRAR);
  if (fit.a0 > 1e5 || fit.a0 < 10) continue;
  perGalA0.push({ name: gname, a0: fit.a0, logA0: fit.logA0, rms: fit.rms, range, Vmax: galaxyMeta[gname].Vmax, inc: galaxyMeta[gname].inc });
}

function getMedianForSample(sampleNames) {
  const vals = perGalA0.filter(g => sampleNames.includes(g.name)).map(g => g.logA0);
  if (vals.length < 3) return null;
  const med = median(vals);
  return { n: vals.length, medianLogA0: +med.toFixed(3), a0: +Math.pow(10, med).toFixed(0), mad: +mad(vals).toFixed(3) };
}

function reportMedian(label, sampleNames) {
  const r = getMedianForSample(sampleNames);
  if (!r) { log("  " + label.padEnd(16) + "— too few —"); return null; }
  log("  " + label.padEnd(16) + "n=" + String(r.n).padEnd(5) +
    "median log(a₀)=" + r.medianLogA0.toFixed(3).padEnd(8) +
    "a₀=" + String(r.a0).padEnd(8) +
    "MAD=" + r.mad.toFixed(3));
  return r;
}

const medians = {};
medians.all = reportMedian("ALL", allGalNames);
medians.clean = reportMedian("CLEAN", cleanSample);
medians.gold = reportMedian("GOLD", goldSample);
medians.goldTight = reportMedian("GOLD+i45", goldTight);
medians.highMass = reportMedian("HIGH-MASS", highMass);
log("");

sep();
log("LAYER 3: CONVERGENCE TABLE");
sep();
log("");

log("  ┌────────────┬──────┬───────┬────────────┬────────────┬──────────┐");
log("  │ Sample     │ n_gal│ n_pts │ Global fit │ Per-gal med│ G/M ratio│");
log("  ├────────────┼──────┼───────┼────────────┼────────────┼──────────┤");
for (const layer of layers) {
  const lr = layerResults[layer.name];
  const medKey = layer.name === "ALL" ? "all" : layer.name === "CLEAN" ? "clean" : layer.name === "GOLD" ? "gold" : layer.name === "GOLD+i45" ? "goldTight" : "highMass";
  const med = medians[medKey];
  const medA0 = med ? med.a0 : "—";
  const ratio = med ? (lr.meanA0 / med.a0).toFixed(2) : "—";
  log("  │ " + layer.name.padEnd(11) + "│ " + String(lr.nGal).padEnd(5) + "│ " + String(lr.nPts).padEnd(6) + "│ " +
    String(lr.meanA0).padEnd(11) + "│ " + String(medA0).padEnd(11) + "│ ×" + String(ratio).padEnd(9) + "│");
}
log("  └────────────┴──────┴───────┴────────────┴────────────┴──────────┘");
log("");

sep();
log("LAYER 4: BLIND RESIDUAL AUDIT — GOLD vs GOLD+i45");
sep();
log("");

const goldFit = fitA0(getPointsForGalaxies(goldSample), mcgaughRAR);
const goldTightFit = fitA0(getPointsForGalaxies(goldTight), mcgaughRAR);

log("  GOLD fit:      a₀ = " + goldFit.a0.toFixed(1) + ", RMS = " + goldFit.rms.toFixed(4) + ", n = " + getPointsForGalaxies(goldSample).length);
log("  GOLD+i45 fit:  a₀ = " + goldTightFit.a0.toFixed(1) + ", RMS = " + goldTightFit.rms.toFixed(4) + ", n = " + getPointsForGalaxies(goldTight).length);
log("");

const auditGold = blindAudit("GOLD", goldSample, goldFit.a0);
const auditGoldTight = blindAudit("GOLD+i45", goldTight, goldTightFit.a0);

function printAudit(audit) {
  log("  BLIND RESIDUAL AUDIT: " + audit.sampleName + " (n=" + audit.nGalaxies + ", a₀=" + audit.a0Used.toFixed(0) + ")");
  log("  " + "Variable".padEnd(22) + "r".padEnd(10) + "n".padEnd(6) + "|r|>0.3?  Verdict");
  log("  " + "─".repeat(60));
  for (const c of audit.correlations) {
    log("  " + c.variable.padEnd(22) + c.r.toFixed(3).padEnd(10) + String(c.n).padEnd(6) +
      (c.significant ? "YES ⚠    CONTAMINATED" : "no       CLEAN"));
  }
  log("  " + "─".repeat(60));
  log("  Verdict: " + (audit.clean ? "✓ ALL CLEAN" : "⚠ ISSUES (see above)"));
  log("");
}

printAudit(auditGold);
printAudit(auditGoldTight);

log("  ┌──────────────────────────────────────────────────────────────────┐");
if (auditGoldTight.clean) {
  log("  │ GOLD+i45 RESIDUAL AUDIT: ✓ FULLY PASSED                       │");
  log("  │ All 6 variables clean. Inclination tightening resolved the    │");
  log("  │ marginal signal from GOLD (inc≥30°).                          │");
  log("  │ GOLD+i45 is the recommended final measurement sample.         │");
} else {
  log("  │ GOLD+i45 RESIDUAL AUDIT: ⚠ ISSUES REMAIN                      │");
  for (const c of auditGoldTight.correlations.filter(c => c.significant)) {
    log("  │   " + c.variable.padEnd(18) + " r = " + c.r.toFixed(3).padEnd(20) + "│");
  }
  log("  │ Further investigation needed.                                  │");
}
log("  └──────────────────────────────────────────────────────────────────┘");

log("");
sep();
log("LAYER 5: INTERPOLATION STABILITY (GOLD+i45)");
sep();
log("");

const goldTightPts = getPointsForGalaxies(goldTight);
const tightInterpResults = [];
for (let i = 0; i < INTERP_FUNCS.length; i++) {
  const fit = fitA0(goldTightPts, INTERP_FUNCS[i]);
  tightInterpResults.push({ func: INTERP_SET[i], a0: fit.a0, logA0: fit.logA0, rms: fit.rms });
  log("  " + INTERP_SET[i].padEnd(18) + "a₀ = " + fit.a0.toFixed(1).padEnd(10) + "RMS = " + fit.rms.toFixed(4));
}
const tightA0vals = tightInterpResults.map(r => r.a0);
const tightMean = tightA0vals.reduce((a, b) => a + b) / tightA0vals.length;
const tightSpread = (Math.max(...tightA0vals) - Math.min(...tightA0vals)) / tightMean * 100;
log("  Mean: " + tightMean.toFixed(1) + " | Spread: " + tightSpread.toFixed(1) + "%");

log("");
sep();
log("LAYER 6: FORMAL UNCERTAINTY BUDGET");
sep();
log("");

const refA0 = layerResults["GOLD+i45"].meanA0;
const refLogA0 = Math.log10(refA0);

const interpA0s = tightInterpResults.map(r => r.logA0);
const interpRange = Math.max(...interpA0s) - Math.min(...interpA0s);
log("  Reference: GOLD+i45 mean a₀ = " + refA0.toFixed(1) + " (log = " + refLogA0.toFixed(4) + ")");
log("");

log("  1. INTERPOLATION FUNCTION UNCERTAINTY:");
log("     Range in log(a₀): " + interpRange.toFixed(4) + " dex");
log("     Half-range: ±" + (interpRange / 2).toFixed(4) + " dex");
log("     Fractional: ±" + (tightSpread / 2).toFixed(1) + "%");
log("");

const sampleA0s = [];
for (const name of ["GOLD", "GOLD+i45", "HIGH-MASS", "CLEAN"]) {
  if (layerResults[name]) sampleA0s.push(Math.log10(layerResults[name].meanA0));
}
const sampleRange = Math.max(...sampleA0s) - Math.min(...sampleA0s);
log("  2. SAMPLE SELECTION UNCERTAINTY:");
log("     Samples compared: GOLD, GOLD+i45, HIGH-MASS, CLEAN");
log("     Range in log(a₀): " + sampleRange.toFixed(4) + " dex");
log("     Half-range: ±" + (sampleRange / 2).toFixed(4) + " dex");
log("");

const goldTightMedian = medians.goldTight;
const globalVsMedian = goldTightMedian ? Math.abs(Math.log10(refA0) - goldTightMedian.medianLogA0) : NaN;
log("  3. GLOBAL vs PER-GALAXY UNCERTAINTY:");
if (goldTightMedian) {
  log("     GOLD+i45 global: log = " + refLogA0.toFixed(3));
  log("     GOLD+i45 median: log = " + goldTightMedian.medianLogA0.toFixed(3));
  log("     Difference: " + globalVsMedian.toFixed(3) + " dex");
  log("     Per-galaxy MAD: " + goldTightMedian.mad.toFixed(3) + " dex");
}
log("");

const upsilonTests = [];
for (const uFactor of [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15]) {
  const scaledPts = goldTightPts.map(p => ({ ...p, x: p.x + Math.log10(uFactor) }));
  const fit = fitA0(scaledPts, mcgaughRAR);
  upsilonTests.push({ factor: uFactor, a0: fit.a0, logA0: fit.logA0 });
}
const upsilonRange = Math.max(...upsilonTests.map(t => t.logA0)) - Math.min(...upsilonTests.map(t => t.logA0));
log("  4. MASS-TO-LIGHT RATIO (Υ★) UNCERTAINTY:");
log("     Tested: Υ★ × [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15]");
for (const t of upsilonTests) {
  log("       Υ★×" + t.factor.toFixed(2) + ":  a₀ = " + t.a0.toFixed(0).padEnd(8) + "log = " + t.logA0.toFixed(4));
}
log("     Range in log(a₀): " + upsilonRange.toFixed(4) + " dex (±" + (upsilonRange / 2).toFixed(4) + " dex)");
log("");

const incTests = [];
for (const incCut of [30, 35, 40, 45, 50, 55, 60]) {
  const sample = allGalNames.filter(n => {
    const m = galaxyMeta[n]; const range = getGbarRange(n);
    return m.Vmax >= 50 && m.inc >= incCut && m.n >= 5 && range >= 1.0;
  });
  const pts = getPointsForGalaxies(sample);
  if (pts.length < 50) continue;
  const fit = fitA0(pts, mcgaughRAR);
  incTests.push({ incCut, nGal: sample.length, nPts: pts.length, a0: fit.a0, logA0: fit.logA0, rms: fit.rms });
}
const incRange = incTests.length > 1 ? Math.max(...incTests.map(t => t.logA0)) - Math.min(...incTests.map(t => t.logA0)) : 0;
log("  5. INCLINATION CUT SENSITIVITY:");
for (const t of incTests) {
  log("     inc≥" + String(t.incCut).padEnd(3) + ": n=" + String(t.nGal).padEnd(4) + "pts=" + String(t.nPts).padEnd(5) +
    "a₀=" + t.a0.toFixed(0).padEnd(8) + "log=" + t.logA0.toFixed(4).padEnd(9) + "RMS=" + t.rms.toFixed(4));
}
log("     Range in log(a₀): " + incRange.toFixed(4) + " dex (±" + (incRange / 2).toFixed(4) + " dex)");
log("");

const rangeTests = [];
for (const minRange of [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]) {
  const sample = allGalNames.filter(n => {
    const m = galaxyMeta[n]; const range = getGbarRange(n);
    return m.Vmax >= 50 && m.inc >= 45 && m.n >= 5 && range >= minRange;
  });
  const pts = getPointsForGalaxies(sample);
  if (pts.length < 30) continue;
  const fit = fitA0(pts, mcgaughRAR);
  rangeTests.push({ minRange, nGal: sample.length, nPts: pts.length, a0: fit.a0, logA0: fit.logA0 });
}
const rangeTestRange = rangeTests.length > 1 ? Math.max(...rangeTests.map(t => t.logA0)) - Math.min(...rangeTests.map(t => t.logA0)) : 0;
log("  6. DYNAMIC RANGE CUT SENSITIVITY:");
for (const t of rangeTests) {
  log("     range≥" + t.minRange.toFixed(1) + " dex: n=" + String(t.nGal).padEnd(4) + "pts=" + String(t.nPts).padEnd(5) +
    "a₀=" + t.a0.toFixed(0).padEnd(8) + "log=" + t.logA0.toFixed(4));
}
log("     Range in log(a₀): " + rangeTestRange.toFixed(4) + " dex (±" + (rangeTestRange / 2).toFixed(4) + " dex)");
log("");

log("  ╔══════════════════════════════════════════════════════════════════╗");
log("  ║  FORMAL UNCERTAINTY BUDGET (all in dex of log a₀)              ║");
log("  ╠══════════════════════════════════════════════════════════════════╣");
log("  ║  Source                        │ ±σ (dex)                      ║");
log("  ╠════════════════════════════════╪═══════════════════════════════╣");
const budgetItems = [
  { source: "Interpolation function", val: interpRange / 2 },
  { source: "Sample selection", val: sampleRange / 2 },
  { source: "Global vs per-galaxy", val: goldTightMedian ? globalVsMedian : 0.1 },
  { source: "Υ★ (±15%)", val: upsilonRange / 2 },
  { source: "Inclination cut", val: incRange / 2 },
  { source: "Dynamic range cut", val: rangeTestRange / 2 },
];
let totalVar = 0;
for (const b of budgetItems) {
  log("  ║  " + b.source.padEnd(31) + "│ ±" + b.val.toFixed(4).padEnd(28) + "║");
  totalVar += b.val * b.val;
}
const totalSigma = Math.sqrt(totalVar);
log("  ╠════════════════════════════════╪═══════════════════════════════╣");
log("  ║  TOTAL (quadrature)            │ ±" + totalSigma.toFixed(4).padEnd(28) + "║");
log("  ╚════════════════════════════════╧═══════════════════════════════╝");
log("");

const a0Lower = Math.pow(10, refLogA0 - totalSigma);
const a0Upper = Math.pow(10, refLogA0 + totalSigma);
log("  FINAL RESULT:");
log("    log(a₀) = " + refLogA0.toFixed(3) + " ± " + totalSigma.toFixed(3) + " dex");
log("    a₀ = " + refA0.toFixed(0) + " [" + a0Lower.toFixed(0) + ", " + a0Upper.toFixed(0) + "] " + UNIT_LABEL);
log("    a₀ = " + (refA0 * MS2_CONV).toExponential(3) + " [" + (a0Lower * MS2_CONV).toExponential(2) + ", " + (a0Upper * MS2_CONV).toExponential(2) + "] m/s²");
log("    a₀/(cH₀/2π) = " + (refA0 * MS2_CONV / CH0_2PI).toFixed(3));

log("");
sep();
log("LAYER 7: LOW-MASS STRESS SAMPLE (excluded from measurement)");
sep();
log("");

const lowMass = allGalNames.filter(n => galaxyMeta[n].Vmax < 80 && galaxyMeta[n].Vmax > 0);
const lowMassPts = getPointsForGalaxies(lowMass);
const lowMassRanges = lowMass.map(n => getGbarRange(n));
log("  Low-mass galaxies (V < 80 km/s): " + lowMass.length);
log("  Points: " + lowMassPts.length);
log("  Median g_bar range: " + median(lowMassRanges).toFixed(2) + " dex");
log("  Fraction range < 0.5 dex: " + (lowMassRanges.filter(r => r < 0.5).length / lowMassRanges.length * 100).toFixed(0) + "%");
log("  Fraction range < 1.0 dex: " + (lowMassRanges.filter(r => r < 1.0).length / lowMassRanges.length * 100).toFixed(0) + "%");
if (lowMassPts.length >= 30) {
  const lmFit = fitA0(lowMassPts, mcgaughRAR);
  log("  Fit: a₀ = " + lmFit.a0.toFixed(0) + " (log = " + lmFit.logA0.toFixed(2) + "), RMS = " + lmFit.rms.toFixed(4));
  log("  ⚠ STRESS TEST ONLY. Cannot constrain a₀.");
}

log("");
sep();
log("FINAL MEASUREMENT SUMMARY");
sep();
log("");

const bestSample = auditGoldTight.clean ? "GOLD+i45" : "GOLD";
const bestLR = layerResults[bestSample];
const bestMed = bestSample === "GOLD+i45" ? medians.goldTight : medians.gold;

log("  ╔════════════════════════════════════════════════════════════════════╗");
log("  ║  FINAL a₀ MEASUREMENT TABLE                                      ║");
log("  ╠════════════════════════════════════════════════════════════════════╣");
log("  ║  Sample     │ Method      │ a₀ (km/s)²/kpc  │ a₀/(cH₀/2π)      ║");
log("  ╠════════════════════════════════════════════════════════════════════╣");
for (const layer of layers) {
  const lr = layerResults[layer.name];
  log("  ║  " + layer.name.padEnd(11) + "│ Global fit  │ " + String(lr.meanA0).padEnd(16) + "│ " + String(lr.ratio_cH0).padEnd(18) + "║");
}
log("  ╠════════════════════════════════════════════════════════════════════╣");
for (const [key, label] of [["all", "ALL"], ["clean", "CLEAN"], ["gold", "GOLD"], ["goldTight", "GOLD+i45"], ["highMass", "HIGH-MASS"]]) {
  if (medians[key]) {
    log("  ║  " + label.padEnd(11) + "│ Per-gal med │ " + String(medians[key].a0).padEnd(16) + "│ " + ((medians[key].a0 * MS2_CONV) / CH0_2PI).toFixed(3).padEnd(18) + "║");
  }
}
log("  ╠════════════════════════════════════════════════════════════════════╣");
log("  ║  RECOMMENDED: " + bestSample.padEnd(52) + "║");
log("  ║  a₀ = " + bestLR.meanA0 + " ± " + totalSigma.toFixed(3) + " dex = " + bestLR.a0_ms2 + " m/s²".padEnd(18) + "    ║");
log("  ║  a₀/(cH₀/2π) = " + String(bestLR.ratio_cH0).padEnd(50) + "  ║");
log("  ║  Range: [" + a0Lower.toFixed(0) + ", " + a0Upper.toFixed(0) + "] " + UNIT_LABEL.padEnd(35) + "    ║");
log("  ║  Residual audit: " + (auditGoldTight.clean ? "✓ ALL 6 CLEAN" : "⚠ ISSUES").padEnd(48) + "  ║");
log("  ║  Interp spread: " + tightSpread.toFixed(1) + "%".padEnd(49) + "  ║");
log("  ╚════════════════════════════════════════════════════════════════════╝");

log("");
sep();
log("REPRODUCIBILITY MANIFEST");
sep();
log("");
log("  Pipeline version: " + VERSION);
log("  Timestamp: " + TIMESTAMP);
log("  Input files: transition-scale.json, rar-analysis-real.json");
log("  Υ_disk = " + UPSILON_DISK + ", Υ_bulge = " + UPSILON_BULGE);
log("  Interpolations: " + INTERP_SET.join(", "));
log("  Excluded: Our ν(x)=x/(1+x) (ΔAIC=99)");
log("  Fit: log(a₀) ∈ [2.0, 5.0], step 0.01, 100 golden section refinements");
log("  Residual: log10(g_obs) - log10(interpFunc(g_bar, a₀))");
log("  GOLD cuts: V≥50, inc≥30°, n≥5, range≥1.0 dex");
log("  GOLD+i45 cuts: V≥50, inc≥45°, n≥5, range≥1.0 dex");
log("  Audit threshold: |r| > 0.3 = significant");
log("  LOW-MASS: stress sample only");
log("  Uncertainty: 6-component budget in quadrature");

const output = {
  version: VERSION,
  timestamp: TIMESTAMP,
  lockedParameters: {
    upsilonDisk: UPSILON_DISK, upsilonBulge: UPSILON_BULGE,
    interpolations: INTERP_SET,
    fitRange: "log(a0) in [2.0, 5.0]",
    residualDef: "log10(gobs) - log10(interpFunc(gbar, a0))",
    units: UNIT_LABEL, ms2Conv: MS2_CONV
  },
  samples: {
    all: { nGal: allGalNames.length, nPts: ALL_POINTS.length },
    clean: { nGal: cleanSample.length, nPts: getPointsForGalaxies(cleanSample).length, cuts: "V>=80, inc>=30, n>=5, range>=0.5" },
    gold: { nGal: goldSample.length, nPts: getPointsForGalaxies(goldSample).length, cuts: "V>=50, inc>=30, n>=5, range>=1.0" },
    goldTight: { nGal: goldTight.length, nPts: goldTightPts.length, cuts: "V>=50, inc>=45, n>=5, range>=1.0" },
    highMass: { nGal: highMass.length, nPts: getPointsForGalaxies(highMass).length, cuts: "V>150" }
  },
  globalFits: layerResults,
  perGalaxyMedians: medians,
  residualAudit: {
    gold: auditGold,
    goldTight: auditGoldTight
  },
  interpStability: {
    goldTight: tightInterpResults.map(r => ({ func: r.func, a0: +r.a0.toFixed(1), rms: +r.rms.toFixed(4) })),
    spreadPct: +tightSpread.toFixed(1)
  },
  uncertaintyBudget: {
    components: budgetItems.map(b => ({ source: b.source, halfRange_dex: +b.val.toFixed(4) })),
    total_dex: +totalSigma.toFixed(4),
    refLogA0: +refLogA0.toFixed(4),
    a0_lower: +a0Lower.toFixed(0),
    a0_upper: +a0Upper.toFixed(0)
  },
  sensitivityTests: {
    upsilon: upsilonTests.map(t => ({ factor: t.factor, a0: +t.a0.toFixed(0), logA0: +t.logA0.toFixed(4) })),
    inclination: incTests.map(t => ({ incCut: t.incCut, nGal: t.nGal, a0: +t.a0.toFixed(0), logA0: +t.logA0.toFixed(4), rms: +t.rms.toFixed(4) })),
    dynamicRange: rangeTests.map(t => ({ minRange: t.minRange, nGal: t.nGal, a0: +t.a0.toFixed(0), logA0: +t.logA0.toFixed(4) }))
  },
  recommendedMeasurement: {
    sample: bestSample,
    method: "global fit (mean of 3 interpolation functions)",
    a0_galUnits: +bestLR.meanA0.toFixed(1),
    a0_ms2: bestLR.a0_ms2,
    ratio_cH0_2pi: bestLR.ratio_cH0,
    logA0: +refLogA0.toFixed(4),
    uncertainty_dex: +totalSigma.toFixed(4),
    a0_lower: +a0Lower.toFixed(0),
    a0_upper: +a0Upper.toFixed(0),
    residualAuditPassed: auditGoldTight.clean,
    interpSpreadPct: +tightSpread.toFixed(1)
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'gold-standard-results.json'), JSON.stringify(output, null, 2));
log("\nResults saved to gold-standard-results.json");
log("Pipeline v" + VERSION + " complete.\n");
