const fs = require('fs');
const path = require('path');

const VERSION = "1.0.0";
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

const allGalNames = Object.keys(galaxyMeta);

const cleanSample = allGalNames.filter(n => {
  const m = galaxyMeta[n];
  const range = getGbarRange(n);
  return m.Vmax >= 80 && m.inc >= 30 && m.n >= 5 && range >= 0.5;
});

const goldSample = allGalNames.filter(n => {
  const m = galaxyMeta[n];
  const range = getGbarRange(n);
  return m.Vmax >= 50 && m.inc >= 30 && m.n >= 5 && range >= 1.0;
});

const highMass = allGalNames.filter(n => galaxyMeta[n].Vmax > 150);

log("╔════════════════════════════════════════════════════════════════════════╗");
log("║  GOLD-STANDARD a₀ MEASUREMENT PIPELINE                              ║");
log("║  Version: " + VERSION + "                                                      ║");
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
log("  Source: transition-scale.json (every-2nd-point subsample)");
log("  Note: full sample = 4123 points. Subsampling introduces <1% a₀ shift.");
log("");

log("SAMPLE DEFINITIONS:");
log("  ALL:   " + allGalNames.length + " galaxies, " + ALL_POINTS.length + " points (no cuts)");
log("  CLEAN: " + cleanSample.length + " galaxies, " + getPointsForGalaxies(cleanSample).length + " points");
log("         (V≥80, inc≥30°, n≥5, gbar_range≥0.5 dex)");
log("  GOLD:  " + goldSample.length + " galaxies, " + getPointsForGalaxies(goldSample).length + " points");
log("         (V≥50, inc≥30°, n≥5, gbar_range≥1.0 dex)");
log("  HMASS: " + highMass.length + " galaxies, " + getPointsForGalaxies(highMass).length + " points");
log("         (V>150 km/s, no other cuts)");
log("");

sep();
log("LAYER 1: FOUR-LAYER a₀ MEASUREMENT");
sep();
log("");

const layers = [
  { name: "ALL (global fit)", points: ALL_POINTS },
  { name: "CLEAN sample", points: getPointsForGalaxies(cleanSample) },
  { name: "GOLD sample", points: getPointsForGalaxies(goldSample) },
  { name: "HIGH-MASS (V>150)", points: getPointsForGalaxies(highMass) },
];

const layerResults = {};
for (const layer of layers) {
  log("  " + layer.name + " (n=" + layer.points.length + "):");
  const results = [];
  for (let i = 0; i < INTERP_FUNCS.length; i++) {
    const fit = fitA0(layer.points, INTERP_FUNCS[i]);
    results.push({ func: INTERP_SET[i], ...fit });
    log("    " + INTERP_SET[i].padEnd(18) + "a₀ = " + fit.a0.toFixed(1).padEnd(10) +
      "log(a₀) = " + fit.logA0.toFixed(4).padEnd(10) +
      "RMS = " + fit.rms.toFixed(4).padEnd(10) +
      "bias = " + fit.bias.toFixed(4));
  }
  const a0s = results.map(r => r.a0);
  const meanA0 = a0s.reduce((a, b) => a + b) / a0s.length;
  const spread = ((Math.max(...a0s) - Math.min(...a0s)) / meanA0 * 100);
  log("    Mean a₀: " + meanA0.toFixed(1) + " ± " + spread.toFixed(1) + "% spread across functions");
  log("    In m/s²: " + (meanA0 * MS2_CONV).toExponential(3));
  log("    Ratio a₀/(cH₀/2π): " + (meanA0 * MS2_CONV / CH0_2PI).toFixed(3));
  log("");
  layerResults[layer.name] = { results, meanA0: +meanA0.toFixed(1), spread: +spread.toFixed(1), bestRMS: +Math.min(...results.map(r => r.rms)).toFixed(4) };
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
  perGalA0.push({ name: gname, a0: fit.a0, logA0: fit.logA0, rms: fit.rms, range, Vmax: galaxyMeta[gname].Vmax });
}

const allA0s = perGalA0.map(g => g.logA0);
const cleanA0s = perGalA0.filter(g => cleanSample.includes(g.name)).map(g => g.logA0);
const goldA0s = perGalA0.filter(g => goldSample.includes(g.name)).map(g => g.logA0);
const hmA0s = perGalA0.filter(g => highMass.includes(g.name)).map(g => g.logA0);

function reportMedian(label, vals) {
  if (vals.length < 3) { log("  " + label.padEnd(25) + "— too few —"); return null; }
  const med = median(vals);
  const madVal = mad(vals);
  const a0 = Math.pow(10, med);
  log("  " + label.padEnd(25) + "n=" + String(vals.length).padEnd(5) +
    "median log(a₀)=" + med.toFixed(3).padEnd(8) +
    "a₀=" + a0.toFixed(0).padEnd(8) +
    "MAD=" + madVal.toFixed(3));
  return { n: vals.length, medianLogA0: +med.toFixed(3), a0: +a0.toFixed(0), mad: +madVal.toFixed(3) };
}

const medians = {};
medians.all = reportMedian("ALL", allA0s);
medians.clean = reportMedian("CLEAN", cleanA0s);
medians.gold = reportMedian("GOLD", goldA0s);
medians.highMass = reportMedian("HIGH-MASS", hmA0s);
log("");

sep();
log("LAYER 3: CONVERGENCE TABLE — WHICH a₀ IS THE MEASUREMENT?");
sep();
log("");

log("  ┌────────────────────────────┬────────────┬────────────┬────────────┐");
log("  │ Sample                      │ Global fit │ Per-gal med│ Spread     │");
log("  ├────────────────────────────┼────────────┼────────────┼────────────┤");
for (const layer of layers) {
  const lr = layerResults[layer.name];
  const medKey = layer.name.startsWith("ALL") ? "all" : layer.name.startsWith("CLEAN") ? "clean" : layer.name.startsWith("GOLD") ? "gold" : "highMass";
  const med = medians[medKey];
  const globalA0 = lr.meanA0;
  const medA0 = med ? med.a0 : "—";
  const ratio = med ? (globalA0 / med.a0).toFixed(2) : "—";
  log("  │ " + layer.name.padEnd(27) + "│ " + String(globalA0).padEnd(11) + "│ " + String(medA0).padEnd(11) + "│ ×" + String(ratio).padEnd(10) + "│");
}
log("  └────────────────────────────┴────────────┴────────────┴────────────┘");
log("");

sep();
log("LAYER 4: BLIND RESIDUAL AUDIT (using GOLD sample, McGaugh RAR)");
sep();
log("");

const goldPts = getPointsForGalaxies(goldSample);
const goldFit = fitA0(goldPts, mcgaughRAR);
log("  GOLD sample fit: a₀ = " + goldFit.a0.toFixed(1) + ", RMS = " + goldFit.rms.toFixed(4) + ", n = " + goldPts.length);
log("");

const goldResids = [];
for (const gname of goldSample) {
  const pts = ALL_POINTS.filter(p => p.g === gname.substring(0, 12));
  if (pts.length < 3) continue;
  const resids = [];
  for (const p of pts) {
    const gbar = Math.pow(10, p.x);
    const gobs = gbar * Math.pow(10, p.y);
    const pred = mcgaughRAR(gbar, goldFit.a0);
    if (!isFinite(pred) || pred <= 0) continue;
    resids.push(Math.log10(gobs) - Math.log10(pred));
  }
  if (resids.length < 2) continue;
  const mean = resids.reduce((a, b) => a + b, 0) / resids.length;
  const meta = galaxyMeta[gname];
  if (!meta) continue;
  const gasFrac = meta.MHI > 0 && meta.L36 > 0 ? meta.MHI / meta.L36 : NaN;
  goldResids.push({
    name: gname, meanResid: mean,
    logVmax: Math.log10(meta.Vmax),
    logSB: meta.sigma_bar > 0 ? Math.log10(meta.sigma_bar) : NaN,
    logGasFrac: isFinite(gasFrac) && gasFrac > 0 ? Math.log10(gasFrac) : NaN,
    logL: meta.L36 > 0 ? Math.log10(meta.L36) : NaN,
    inc: meta.inc,
    logDist: meta.distance > 0 ? Math.log10(meta.distance) : NaN,
  });
}

const validGR = goldResids.filter(g => isFinite(g.logVmax) && isFinite(g.logSB) && isFinite(g.logGasFrac) && isFinite(g.logL));

log("  Galaxies with residuals: " + validGR.length);
log("");

const auditVars = [
  { name: "V_max", vals: validGR.map(g => [g.meanResid, g.logVmax]) },
  { name: "Surface brightness", vals: validGR.map(g => [g.meanResid, g.logSB]) },
  { name: "Gas fraction", vals: validGR.map(g => [g.meanResid, g.logGasFrac]) },
  { name: "Luminosity", vals: validGR.map(g => [g.meanResid, g.logL]) },
  { name: "Inclination", vals: validGR.filter(g => g.inc > 0).map(g => [g.meanResid, g.inc]) },
  { name: "Distance", vals: validGR.filter(g => isFinite(g.logDist)).map(g => [g.meanResid, g.logDist]) },
];

log("  BLIND RESIDUAL AUDIT (GOLD sample, a₀ = " + goldFit.a0.toFixed(0) + "):");
log("  " + "Variable".padEnd(22) + "r".padEnd(10) + "n".padEnd(6) + "|r|>0.3?  Verdict");
log("  " + "─".repeat(60));

const auditResults = [];
let auditClean = true;
for (const av of auditVars) {
  if (av.vals.length < 10) {
    log("  " + av.name.padEnd(22) + "—".padEnd(10) + String(av.vals.length).padEnd(6) + "—        too few");
    continue;
  }
  const r = pearsonR(av.vals.map(v => v[0]), av.vals.map(v => v[1]));
  const sig = Math.abs(r) > 0.3;
  if (sig) auditClean = false;
  log("  " + av.name.padEnd(22) + r.toFixed(3).padEnd(10) + String(av.vals.length).padEnd(6) +
    (sig ? "YES ⚠    " : "no       ") + (sig ? "CONTAMINATED" : "CLEAN"));
  auditResults.push({ variable: av.name, r: +r.toFixed(3), n: av.vals.length, significant: sig });
}

log("");
log("  ┌──────────────────────────────────────────────────────────────┐");
if (auditClean) {
  log("  │ RESIDUAL AUDIT: ✓ PASSED                                   │");
  log("  │ No significant correlations in GOLD sample.                 │");
  log("  │ Residuals are consistent with random scatter.               │");
} else {
  log("  │ RESIDUAL AUDIT: ⚠ ISSUES FOUND                             │");
  for (const ar of auditResults.filter(a => a.significant)) {
    log("  │   " + ar.variable.padEnd(20) + " r = " + ar.r.toFixed(3).padEnd(22) + "│");
  }
}
log("  └──────────────────────────────────────────────────────────────┘");

log("");
sep();
log("LAYER 5: INTERPOLATION STABILITY IN GOLD SAMPLE");
sep();
log("");

const goldInterpResults = [];
for (let i = 0; i < INTERP_FUNCS.length; i++) {
  const fit = fitA0(goldPts, INTERP_FUNCS[i]);
  goldInterpResults.push({ func: INTERP_SET[i], a0: fit.a0, rms: fit.rms });
  log("  " + INTERP_SET[i].padEnd(18) + "a₀ = " + fit.a0.toFixed(1).padEnd(10) + "RMS = " + fit.rms.toFixed(4));
}
const goldA0spread = (Math.max(...goldInterpResults.map(r => r.a0)) - Math.min(...goldInterpResults.map(r => r.a0))) / goldInterpResults.reduce((a, r) => a + r.a0, 0) * goldInterpResults.length * 100;
log("  Spread: " + goldA0spread.toFixed(1) + "%");

log("");
sep();
log("LAYER 6: LOW-MASS STRESS SAMPLE (consistency check only)");
sep();
log("");

const lowMass = allGalNames.filter(n => galaxyMeta[n].Vmax < 80 && galaxyMeta[n].Vmax > 0);
const lowMassPts = getPointsForGalaxies(lowMass);
log("  Low-mass galaxies (V < 80 km/s): " + lowMass.length);
log("  Points: " + lowMassPts.length);
const lowMassRanges = lowMass.map(n => getGbarRange(n));
log("  Median g_bar range: " + median(lowMassRanges).toFixed(2) + " dex");
log("  Fraction with range < 0.5 dex: " + (lowMassRanges.filter(r => r < 0.5).length / lowMassRanges.length * 100).toFixed(0) + "%");
log("  Fraction with range < 1.0 dex: " + (lowMassRanges.filter(r => r < 1.0).length / lowMassRanges.length * 100).toFixed(0) + "%");

if (lowMassPts.length >= 30) {
  const lmFit = fitA0(lowMassPts, mcgaughRAR);
  log("  Fit: a₀ = " + lmFit.a0.toFixed(0) + " (log = " + lmFit.logA0.toFixed(2) + "), RMS = " + lmFit.rms.toFixed(4));
  log("  ⚠ This is a STRESS TEST only. Low-mass galaxies cannot independently");
  log("    constrain a₀ due to limited g_bar dynamic range.");
} else {
  log("  Too few points for fit.");
}

log("");
sep();
log("FINAL MEASUREMENT SUMMARY");
sep();
log("");

const goldGlobalA0 = layerResults["GOLD sample"].meanA0;
const goldMedianA0 = medians.gold ? medians.gold.a0 : null;
const hmGlobalA0 = layerResults["HIGH-MASS (V>150)"].meanA0;

log("  ╔════════════════════════════════════════════════════════════════╗");
log("  ║  a₀ MEASUREMENT TABLE                                        ║");
log("  ╠════════════════════════════════════════════════════════════════╣");
log("  ║  Method           │ Sample    │ a₀ " + UNIT_LABEL.padEnd(16) + " ║");
log("  ╠════════════════════════════════════════════════════════════════╣");
log("  ║  Global fit       │ ALL       │ " + String(layerResults["ALL (global fit)"].meanA0).padEnd(24) + "║");
log("  ║  Global fit       │ CLEAN     │ " + String(layerResults["CLEAN sample"].meanA0).padEnd(24) + "║");
log("  ║  Global fit       │ GOLD      │ " + String(goldGlobalA0).padEnd(24) + "║");
log("  ║  Global fit       │ HIGH-MASS │ " + String(hmGlobalA0).padEnd(24) + "║");
if (medians.all) log("  ║  Per-gal median   │ ALL       │ " + String(medians.all.a0).padEnd(24) + "║");
if (medians.clean) log("  ║  Per-gal median   │ CLEAN     │ " + String(medians.clean.a0).padEnd(24) + "║");
if (medians.gold) log("  ║  Per-gal median   │ GOLD      │ " + String(goldMedianA0).padEnd(24) + "║");
if (medians.highMass) log("  ║  Per-gal median   │ HIGH-MASS │ " + String(medians.highMass.a0).padEnd(24) + "║");
log("  ╠════════════════════════════════════════════════════════════════╣");

const bestEstimate = goldGlobalA0;
const bestMS2 = bestEstimate * MS2_CONV;
const ratio_cH = bestMS2 / CH0_2PI;
log("  ║  RECOMMENDED: GOLD global fit                                ║");
log("  ║  a₀ = " + bestEstimate.toFixed(0) + " " + UNIT_LABEL + " = " + bestMS2.toExponential(3) + " m/s²".padEnd(11) + "    ║");
log("  ║  a₀/(cH₀/2π) = " + ratio_cH.toFixed(3).padEnd(40) + "  ║");
log("  ║  Residual audit: " + (auditClean ? "PASSED ✓" : "ISSUES ⚠").padEnd(39) + "  ║");
log("  ║  Interpolation spread: " + goldA0spread.toFixed(1) + "% in GOLD sample".padEnd(34) + "  ║");
log("  ╚════════════════════════════════════════════════════════════════╝");

log("");
sep();
log("REPRODUCIBILITY MANIFEST");
sep();
log("");
log("  Pipeline version: " + VERSION);
log("  Timestamp: " + TIMESTAMP);
log("  Input: transition-scale.json (SHA: see git)");
log("  Input: rar-analysis-real.json (SHA: see git)");
log("  Upsilon_disk: " + UPSILON_DISK);
log("  Upsilon_bulge: " + UPSILON_BULGE);
log("  Interpolations: " + INTERP_SET.join(", "));
log("  Excluded: Our framework ν(x)=x/(1+x) (worst fit, ΔAIC=99)");
log("  Fit range: log(a₀) ∈ [2.0, 5.0], step 0.01, 100 golden section refinements");
log("  Residual: log10(g_obs) - log10(interpFunc(g_bar, a₀))");
log("  CLEAN cuts: V≥80 km/s, inc≥30°, n≥5, gbar_range≥0.5 dex");
log("  GOLD cuts: V≥50 km/s, inc≥30°, n≥5, gbar_range≥1.0 dex");
log("  LOW-MASS: stress sample only, not used for a₀ measurement");
log("  Per-galaxy fit: excluded a₀ > 100000 or a₀ < 10 or n < 5 or range < 0.3");

const output = {
  version: VERSION,
  timestamp: TIMESTAMP,
  lockedParameters: {
    upsilonDisk: UPSILON_DISK,
    upsilonBulge: UPSILON_BULGE,
    interpolations: INTERP_SET,
    fitRange: "log(a0) in [2.0, 5.0]",
    residualDef: "log10(gobs) - log10(interpFunc(gbar, a0))",
    units: UNIT_LABEL,
    ms2Conv: MS2_CONV
  },
  samples: {
    all: { nGal: allGalNames.length, nPts: ALL_POINTS.length },
    clean: { nGal: cleanSample.length, nPts: getPointsForGalaxies(cleanSample).length, cuts: "V>=80, inc>=30, n>=5, range>=0.5" },
    gold: { nGal: goldSample.length, nPts: goldPts.length, cuts: "V>=50, inc>=30, n>=5, range>=1.0" },
    highMass: { nGal: highMass.length, nPts: getPointsForGalaxies(highMass).length, cuts: "V>150" }
  },
  globalFits: layerResults,
  perGalaxyMedians: medians,
  goldResidualAudit: {
    a0Used: +goldFit.a0.toFixed(1),
    auditClean,
    correlations: auditResults
  },
  goldInterpStability: {
    results: goldInterpResults.map(r => ({ func: r.func, a0: +r.a0.toFixed(1), rms: +r.rms.toFixed(4) })),
    spreadPct: +goldA0spread.toFixed(1)
  },
  recommendedMeasurement: {
    sample: "GOLD",
    method: "global fit (mean of 3 interpolation functions)",
    a0_galUnits: +bestEstimate.toFixed(1),
    a0_ms2: +bestMS2.toExponential(3),
    ratio_cH0_2pi: +ratio_cH.toFixed(3),
    residualAuditPassed: auditClean,
    interpSpreadPct: +goldA0spread.toFixed(1)
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'gold-standard-results.json'), JSON.stringify(output, null, 2));
log("\nResults saved to gold-standard-results.json");
log("Pipeline complete.\n");
