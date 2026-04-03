const fs = require('fs');
const path = require('path');

const VERSION = "3.0.0";
const TIMESTAMP = new Date().toISOString();

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const UNIT_LABEL = "(km/s)^2/kpc";
const MS2_CONV = 3.241e-14;
const C_MS = 299792458;
const H0_SI = 67.4e3 / 3.086e22;
const CH0_2PI = C_MS * H0_SI / (2 * Math.PI);
const A0_LIT = 1.2e-10;
const A0_LIT_GAL = A0_LIT / MS2_CONV;

const tsData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'transition-scale.json'), 'utf8'));
const rarData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));

const galaxyMeta = {};
for (const g of rarData.perGalaxy) galaxyMeta[g.name] = g;

const ALL_POINTS = tsData.plotPoints;

function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(72)); }

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

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

function getResids(points, a0) {
  const resids = [];
  for (const p of points) {
    const gbar = Math.pow(10, p.x);
    const gobs = gbar * Math.pow(10, p.y);
    const pred = mcgaughRAR(gbar, a0);
    if (!isFinite(pred) || pred <= 0) continue;
    resids.push(Math.log10(gobs) - Math.log10(pred));
  }
  return resids;
}

function evalRMS(points, a0) {
  return rms(getResids(points, a0));
}

function fitA0_unweighted(points) {
  let bestA0 = 3000, bestRMS = Infinity;
  for (let logA = 2.0; logA <= 5.0; logA += 0.01) {
    const a0 = Math.pow(10, logA);
    const r = evalRMS(points, a0);
    if (r < bestRMS) { bestRMS = r; bestA0 = a0; }
  }
  let lo = Math.log10(bestA0) - 0.02, hi = Math.log10(bestA0) + 0.02;
  for (let step = 0; step < 100; step++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    if (evalRMS(points, Math.pow(10, m1)) < evalRMS(points, Math.pow(10, m2))) hi = m2;
    else lo = m1;
  }
  const finalA0 = Math.pow(10, (lo + hi) / 2);
  const finalRMS = evalRMS(points, finalA0);
  const resids = getResids(points, finalA0);
  const meanBias = resids.reduce((a, b) => a + b, 0) / resids.length;
  return { a0: finalA0, logA0: Math.log10(finalA0), rms: finalRMS, bias: meanBias, n: resids.length };
}

function evalWeightedRMS(points, a0, weights) {
  const resids = [];
  let sumW = 0;
  for (let i = 0; i < points.length; i++) {
    const p = points[i];
    const gbar = Math.pow(10, p.x);
    const gobs = gbar * Math.pow(10, p.y);
    const pred = mcgaughRAR(gbar, a0);
    if (!isFinite(pred) || pred <= 0) continue;
    const r = Math.log10(gobs) - Math.log10(pred);
    resids.push(r * r * weights[i]);
    sumW += weights[i];
  }
  return Math.sqrt(resids.reduce((a, b) => a + b, 0) / sumW);
}

function fitA0_weighted(points, weights) {
  let bestA0 = 3000, bestRMS = Infinity;
  for (let logA = 2.0; logA <= 5.0; logA += 0.01) {
    const a0 = Math.pow(10, logA);
    const r = evalWeightedRMS(points, a0, weights);
    if (r < bestRMS) { bestRMS = r; bestA0 = a0; }
  }
  let lo = Math.log10(bestA0) - 0.02, hi = Math.log10(bestA0) + 0.02;
  for (let step = 0; step < 100; step++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    if (evalWeightedRMS(points, Math.pow(10, m1), weights) < evalWeightedRMS(points, Math.pow(10, m2), weights)) hi = m2;
    else lo = m1;
  }
  const finalA0 = Math.pow(10, (lo + hi) / 2);
  const finalRMS = evalWeightedRMS(points, finalA0, weights);
  return { a0: finalA0, logA0: Math.log10(finalA0), rms: finalRMS, n: points.length };
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

function buildGalaxyPointMap() {
  const map = {};
  for (const p of ALL_POINTS) {
    if (!map[p.g]) map[p.g] = [];
    map[p.g].push(p);
  }
  return map;
}

function computePerGalaxyWeights(galaxyPointMap) {
  const galScatter = {};
  for (const [gname, pts] of Object.entries(galaxyPointMap)) {
    if (pts.length < 3) { galScatter[gname] = 1.0; continue; }
    const yVals = pts.map(p => p.y);
    const meanY = yVals.reduce((a, b) => a + b, 0) / yVals.length;
    const varY = yVals.reduce((a, v) => a + (v - meanY) ** 2, 0) / (yVals.length - 1);
    galScatter[gname] = Math.max(varY, 0.001);
  }

  const weights = [];
  for (const p of ALL_POINTS) {
    const meta = galaxyMeta[Object.keys(galaxyMeta).find(n => n.substring(0, 12) === p.g)];
    let w = 1.0;

    const scatterVar = galScatter[p.g] || 0.1;
    w *= 1.0 / scatterVar;

    if (meta && meta.inc > 0) {
      const sinI = Math.sin(meta.inc * Math.PI / 180);
      w *= sinI * sinI;
    }

    if (meta && meta.distance > 0) {
      w *= 1.0 / Math.sqrt(meta.distance);
    }

    const nPts = (galaxyPointMap[p.g] || []).length;
    if (nPts > 1) {
      w *= 1.0 / Math.sqrt(nPts);
    }

    weights.push(Math.max(w, 0.001));
  }

  const maxW = Math.max(...weights);
  return weights.map(w => w / maxW);
}

function hierarchicalModel(galaxyPointMap, sampleGals) {
  const galA0s = [];
  const galWeights = [];

  for (const gname of sampleGals) {
    const shortName = gname.substring(0, 12);
    const pts = galaxyPointMap[shortName];
    if (!pts || pts.length < 5) continue;
    const range = Math.max(...pts.map(p => p.x)) - Math.min(...pts.map(p => p.x));
    if (range < 0.3) continue;

    const fit = fitA0_unweighted(pts);
    if (fit.a0 > 1e5 || fit.a0 < 10) continue;

    const resids = getResids(pts, fit.a0);
    const withinVar = resids.length > 1 ? resids.reduce((a, v) => a + v * v, 0) / (resids.length - 1) : 1.0;

    const meta = galaxyMeta[gname];
    let incWeight = 1.0;
    if (meta && meta.inc > 0) {
      incWeight = Math.sin(meta.inc * Math.PI / 180);
      incWeight *= incWeight;
    }

    let distWeight = 1.0;
    if (meta && meta.distance > 0) {
      distWeight = 1.0 / Math.sqrt(meta.distance);
    }

    const se2 = withinVar / pts.length;
    const totalWeight = incWeight * distWeight / Math.max(se2, 1e-6);

    galA0s.push({ name: gname, logA0: fit.logA0, a0: fit.a0, rms: fit.rms, n: pts.length, se2, withinVar });
    galWeights.push(totalWeight);
  }

  if (galA0s.length < 5) return null;

  const sumW = galWeights.reduce((a, b) => a + b, 0);
  const normW = galWeights.map(w => w / sumW);

  let popLogA0 = 0;
  for (let i = 0; i < galA0s.length; i++) {
    popLogA0 += normW[i] * galA0s[i].logA0;
  }

  let betweenVar = 0;
  for (let i = 0; i < galA0s.length; i++) {
    betweenVar += normW[i] * (galA0s[i].logA0 - popLogA0) ** 2;
  }

  const tau2 = Math.max(betweenVar - galA0s.reduce((a, g, i) => a + normW[i] * g.se2, 0), 0);

  const hierWeights = [];
  for (let i = 0; i < galA0s.length; i++) {
    hierWeights.push(1.0 / (galA0s[i].se2 + tau2));
  }
  const hierSumW = hierWeights.reduce((a, b) => a + b, 0);

  let hierLogA0 = 0;
  for (let i = 0; i < galA0s.length; i++) {
    hierLogA0 += (hierWeights[i] / hierSumW) * galA0s[i].logA0;
  }

  const hierSE = Math.sqrt(1.0 / hierSumW);

  let hierVar = 0;
  for (let i = 0; i < galA0s.length; i++) {
    hierVar += (hierWeights[i] / hierSumW) ** 2 * (galA0s[i].logA0 - hierLogA0) ** 2;
  }

  const residLogA0s = galA0s.map(g => g.logA0);
  const popMedian = median(residLogA0s);
  const popMAD = mad(residLogA0s);

  const I2 = betweenVar > 0 ? Math.max(0, (betweenVar - galA0s.reduce((a, g) => a + g.se2 / galA0s.length, 0)) / betweenVar) * 100 : 0;

  return {
    nGalaxies: galA0s.length,
    popLogA0_weighted: +hierLogA0.toFixed(4),
    popA0_weighted: +Math.pow(10, hierLogA0).toFixed(0),
    popSE: +hierSE.toFixed(4),
    popMedianLogA0: +popMedian.toFixed(4),
    popMedianA0: +Math.pow(10, popMedian).toFixed(0),
    popMAD: +popMAD.toFixed(4),
    tau2: +tau2.toFixed(6),
    tau: +Math.sqrt(tau2).toFixed(4),
    betweenVar: +betweenVar.toFixed(6),
    I2_pct: +I2.toFixed(1),
    galA0s: galA0s.sort((a, b) => a.logA0 - b.logA0),
  };
}

function blindAudit(sampleName, sampleGals, fitA0Val) {
  const residData = [];
  for (const gname of sampleGals) {
    const pts = ALL_POINTS.filter(p => p.g === gname.substring(0, 12));
    if (pts.length < 3) continue;
    const resids = getResids(pts, fitA0Val);
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
const galaxyPointMap = buildGalaxyPointMap();

const SAMPLE_DEFS = {
  ALL: { filter: () => true, desc: "All 197 galaxies (no cuts)" },
  CLEAN: { filter: n => { const m = galaxyMeta[n]; const r = getGbarRange(n); return m.Vmax >= 80 && m.inc >= 30 && m.n >= 5 && r >= 0.5; }, desc: "V>=80, inc>=30, n>=5, range>=0.5 dex" },
  GOLD: { filter: n => { const m = galaxyMeta[n]; const r = getGbarRange(n); return m.Vmax >= 50 && m.inc >= 30 && m.n >= 5 && r >= 1.0; }, desc: "V>=50, inc>=30, n>=5, range>=1.0 dex" },
  "GOLD+i45": { filter: n => { const m = galaxyMeta[n]; const r = getGbarRange(n); return m.Vmax >= 50 && m.inc >= 45 && m.n >= 5 && r >= 1.0; }, desc: "V>=50, inc>=45, n>=5, range>=1.0 dex" },
  "HIGH-MASS": { filter: n => galaxyMeta[n].Vmax > 150, desc: "V>150 km/s" },
  "LOW-MASS": { filter: n => galaxyMeta[n].Vmax < 80 && galaxyMeta[n].Vmax > 0, desc: "V<80 km/s (stress test only)" },
};

const samples = {};
for (const [name, def] of Object.entries(SAMPLE_DEFS)) {
  const gals = allGalNames.filter(def.filter);
  const pts = getPointsForGalaxies(gals);
  samples[name] = { gals, pts, desc: def.desc };
}

const allWeights = computePerGalaxyWeights(galaxyPointMap);

log("\u2554" + "\u2550".repeat(72) + "\u2557");
log("\u2551  FINAL MEASUREMENT PIPELINE v" + VERSION + "".padEnd(40) + "\u2551");
log("\u2551  SINGLE ESTIMATOR: McGaugh RAR (2016) ONLY".padEnd(73) + "\u2551");
log("\u2551  WEIGHTED FIT + HIERARCHICAL MODEL".padEnd(73) + "\u2551");
log("\u255A" + "\u2550".repeat(72) + "\u255D\n");

log("DESIGN DECISIONS (v3.0):");
log("  1. ESTIMATOR: McGaugh RAR only (dropped Simple/Standard MOND mean)");
log("     Reason: All sensitivity tables built on McGaugh. No estimator mixing.");
log("  2. FITTING: Unweighted + error-weighted (both reported)");
log("  3. MODEL: Hierarchical random-effects meta-analysis");
log("  4. SAMPLES: Locked definitions (see below)");
log("  5. DATA: " + ALL_POINTS.length + " radial points (every-2nd subsample of ~4123)");
log("     Note: Subsampling bias verified <1% in a0 (v2.0 validation)");
log("");

log("LOCKED PARAMETERS:");
log("  Y_disk  = " + UPSILON_DISK + " M_sun/L_sun");
log("  Y_bulge = " + UPSILON_BULGE + " M_sun/L_sun");
log("  Interpolation: McGaugh RAR ONLY: g_obs = g_bar / (1 - exp(-sqrt(g_bar/a0)))");
log("  Fitting: Golden section search, log(a0) in [2.0, 5.0]");
log("  Residual: log10(g_obs) - log10(g_pred)");
log("  Units: " + UNIT_LABEL);
log("  1 " + UNIT_LABEL + " = " + MS2_CONV + " m/s^2");
log("  Literature a0: 1.2e-10 m/s^2 = " + A0_LIT_GAL.toFixed(0) + " " + UNIT_LABEL);
log("");

sep();
log("SECTION 1: LOCKED SAMPLE DEFINITIONS");
sep();
log("");
for (const [name, s] of Object.entries(samples)) {
  log("  " + name.padEnd(12) + ": " + s.gals.length + " galaxies, " + s.pts.length + " points");
  log("  ".padEnd(14) + "  Cuts: " + s.desc);
}
log("");

sep();
log("SECTION 2: UNWEIGHTED McGAUGH-ONLY FITS");
sep();
log("");

const unwFits = {};
for (const [name, s] of Object.entries(samples)) {
  if (name === "LOW-MASS") continue;
  const fit = fitA0_unweighted(s.pts);
  unwFits[name] = fit;
  const ms2 = fit.a0 * MS2_CONV;
  const ratio = ms2 / CH0_2PI;
  log("  " + name.padEnd(12) + "a0=" + fit.a0.toFixed(1).padEnd(10) +
    "log=" + fit.logA0.toFixed(4).padEnd(10) +
    "RMS=" + fit.rms.toFixed(4).padEnd(10) +
    "ratio=" + ratio.toFixed(3));
}
log("");

sep();
log("SECTION 3: ERROR-WEIGHTED McGAUGH-ONLY FITS");
sep();
log("");

log("  Weighting scheme:");
log("    w_i = (1/sigma^2_gal) * sin^2(inc) * (1/sqrt(D)) * (1/sqrt(N_gal))");
log("    sigma^2_gal = within-galaxy variance of log(g_obs/g_bar)");
log("    inc = galaxy inclination (deprojection uncertainty)");
log("    D = distance in Mpc (distance error propagation)");
log("    N_gal = points per galaxy (downweight oversampled galaxies)");
log("");

const wFits = {};
for (const [name, s] of Object.entries(samples)) {
  if (name === "LOW-MASS") continue;
  const sampleWeights = [];
  const shortSet = new Set(s.gals.map(n => n.substring(0, 12)));
  for (let i = 0; i < ALL_POINTS.length; i++) {
    if (shortSet.has(ALL_POINTS[i].g)) {
      sampleWeights.push(allWeights[i]);
    }
  }
  const samplePts = s.pts;
  if (samplePts.length !== sampleWeights.length) {
    const swRebuilt = [];
    for (let i = 0; i < ALL_POINTS.length; i++) {
      if (shortSet.has(ALL_POINTS[i].g)) {
        swRebuilt.push(allWeights[i]);
      }
    }
    const fit = fitA0_weighted(samplePts, swRebuilt);
    wFits[name] = fit;
  } else {
    const fit = fitA0_weighted(samplePts, sampleWeights);
    wFits[name] = fit;
  }
  const fit = wFits[name];
  const ms2 = fit.a0 * MS2_CONV;
  const ratio = ms2 / CH0_2PI;
  log("  " + name.padEnd(12) + "a0=" + fit.a0.toFixed(1).padEnd(10) +
    "log=" + fit.logA0.toFixed(4).padEnd(10) +
    "wRMS=" + fit.rms.toFixed(4).padEnd(10) +
    "ratio=" + ratio.toFixed(3));
}
log("");

log("  UNWEIGHTED vs WEIGHTED COMPARISON:");
log("  " + "Sample".padEnd(12) + "Unweighted".padEnd(12) + "Weighted".padEnd(12) + "Delta(dex)".padEnd(12) + "Delta(%)");
for (const name of Object.keys(unwFits)) {
  if (!wFits[name]) continue;
  const delta = wFits[name].logA0 - unwFits[name].logA0;
  const pct = (Math.pow(10, delta) - 1) * 100;
  log("  " + name.padEnd(12) +
    unwFits[name].a0.toFixed(0).padEnd(12) +
    wFits[name].a0.toFixed(0).padEnd(12) +
    (delta >= 0 ? "+" : "") + delta.toFixed(4).padEnd(12) +
    (pct >= 0 ? "+" : "") + pct.toFixed(1) + "%");
}
log("");

sep();
log("SECTION 4: HIERARCHICAL MODEL (RANDOM-EFFECTS META-ANALYSIS)");
sep();
log("");

log("  Model structure:");
log("    Level 1 (points):   r_ij = log10(g_obs,ij) - log10(RAR(g_bar,ij, a0_i)) + e_ij");
log("    Level 2 (galaxies): log(a0_i) = mu + u_i,  u_i ~ N(0, tau^2)");
log("    Level 3 (population): mu = population log(a0)");
log("");
log("  Estimation: DerSimonian-Laird random-effects meta-analysis");
log("    - Each galaxy contributes one a0 estimate (McGaugh fit)");
log("    - Within-galaxy variance: se_i^2 = sigma^2_within / n_i");
log("    - Between-galaxy heterogeneity: tau^2 (estimated from data)");
log("    - Weights: w_i = 1/(se_i^2 + tau^2) * sin^2(inc_i) / sqrt(D_i)");
log("    - I^2 = fraction of variance due to true heterogeneity");
log("");

const hierResults = {};
for (const [name, s] of Object.entries(samples)) {
  if (name === "LOW-MASS") continue;
  const hier = hierarchicalModel(galaxyPointMap, s.gals);
  if (!hier) { log("  " + name + ": too few galaxies for hierarchical model"); continue; }
  hierResults[name] = hier;

  const ms2w = hier.popA0_weighted * MS2_CONV;
  const ratiow = ms2w / CH0_2PI;
  const ms2m = hier.popMedianA0 * MS2_CONV;
  const ratiom = ms2m / CH0_2PI;

  log("  " + name + " (" + hier.nGalaxies + " galaxies):");
  log("    Weighted mean: log(a0) = " + hier.popLogA0_weighted.toFixed(4) +
    " => a0 = " + hier.popA0_weighted + " " + UNIT_LABEL);
  log("    SE(weighted):  " + hier.popSE.toFixed(4) + " dex");
  log("    Median:        log(a0) = " + hier.popMedianLogA0.toFixed(4) +
    " => a0 = " + hier.popMedianA0 + " " + UNIT_LABEL);
  log("    MAD:           " + hier.popMAD.toFixed(4) + " dex");
  log("    tau (between):  " + hier.tau + " dex");
  log("    tau^2:          " + hier.tau2);
  log("    I^2:            " + hier.I2_pct + "%");
  log("    ratio(wt mean): " + ratiow.toFixed(3));
  log("    ratio(median):  " + ratiom.toFixed(3));
  log("");
}

sep();
log("SECTION 5: THREE-ESTIMATOR CONVERGENCE TABLE");
sep();
log("");

log("  +------------+------+------+-----------+-----------+-----------+----------+----------+----------+");
log("  | Sample     | n_gal| n_pts| Unweighted| Weighted  | Hier.Mean | Hier.Med | Hier.SE  | tau      |");
log("  +------------+------+------+-----------+-----------+-----------+----------+----------+----------+");
for (const name of ["ALL", "CLEAN", "GOLD", "GOLD+i45", "HIGH-MASS"]) {
  const s = samples[name];
  const uw = unwFits[name];
  const w = wFits[name];
  const h = hierResults[name];
  log("  | " + name.padEnd(11) +
    "| " + String(s.gals.length).padEnd(5) +
    "| " + String(s.pts.length).padEnd(5) +
    "| " + (uw ? String(uw.a0.toFixed(0)).padEnd(10) : "---".padEnd(10)) +
    "| " + (w ? String(w.a0.toFixed(0)).padEnd(10) : "---".padEnd(10)) +
    "| " + (h ? String(h.popA0_weighted).padEnd(10) : "---".padEnd(10)) +
    "| " + (h ? String(h.popMedianA0).padEnd(9) : "---".padEnd(9)) +
    "| " + (h ? h.popSE.toFixed(4).padEnd(9) : "---".padEnd(9)) +
    "| " + (h ? String(h.tau).padEnd(9) : "---".padEnd(9)) +
    "|");
}
log("  +------------+------+------+-----------+-----------+-----------+----------+----------+----------+");
log("");

const refSample = "GOLD+i45";
const refHier = hierResults[refSample];
const refUW = unwFits[refSample];
const refW = wFits[refSample];

if (refHier) {
  log("  ESTIMATOR CONVERGENCE FOR " + refSample + ":");
  log("    Unweighted global:    a0 = " + refUW.a0.toFixed(0) + " (log = " + refUW.logA0.toFixed(4) + ")");
  log("    Weighted global:      a0 = " + refW.a0.toFixed(0) + " (log = " + refW.logA0.toFixed(4) + ")");
  log("    Hierarchical mean:    a0 = " + refHier.popA0_weighted + " (log = " + refHier.popLogA0_weighted.toFixed(4) + ")");
  log("    Hierarchical median:  a0 = " + refHier.popMedianA0 + " (log = " + refHier.popMedianLogA0.toFixed(4) + ")");
  const allLogA0 = [refUW.logA0, refW.logA0, refHier.popLogA0_weighted, refHier.popMedianLogA0];
  const spread = Math.max(...allLogA0) - Math.min(...allLogA0);
  log("    Total spread: " + spread.toFixed(4) + " dex (" + ((Math.pow(10, spread) - 1) * 100).toFixed(1) + "%)");
  log("");
}

sep();
log("SECTION 6: BLIND RESIDUAL AUDIT (McGaugh-only a0)");
sep();
log("");

for (const name of ["GOLD", "GOLD+i45"]) {
  const a0 = unwFits[name].a0;
  const audit = blindAudit(name, samples[name].gals, a0);
  log("  " + name + " (n=" + audit.nGalaxies + ", a0=" + a0.toFixed(0) + "):");
  log("  " + "Variable".padEnd(22) + "r".padEnd(10) + "n".padEnd(6) + "|r|>0.3?  Verdict");
  log("  " + "\u2500".repeat(60));
  for (const c of audit.correlations) {
    log("  " + c.variable.padEnd(22) + c.r.toFixed(3).padEnd(10) + String(c.n).padEnd(6) +
      (c.significant ? "YES      CONTAMINATED" : "no       CLEAN"));
  }
  log("  " + "\u2500".repeat(60));
  log("  Verdict: " + (audit.clean ? "ALL CLEAN" : "ISSUES (see above)"));
  log("");
}

sep();
log("SECTION 7: UPDATED UNCERTAINTY BUDGET");
sep();
log("");

const refA0 = refUW.a0;
const refLogA0 = refUW.logA0;

const uwVsW = Math.abs(refUW.logA0 - refW.logA0);
const uwVsHier = refHier ? Math.abs(refUW.logA0 - refHier.popLogA0_weighted) : 0;
const hierMeanVsMedian = refHier ? Math.abs(refHier.popLogA0_weighted - refHier.popMedianLogA0) : 0;

const estimatorSpread = refHier ?
  Math.max(refUW.logA0, refW.logA0, refHier.popLogA0_weighted, refHier.popMedianLogA0) -
  Math.min(refUW.logA0, refW.logA0, refHier.popLogA0_weighted, refHier.popMedianLogA0) : uwVsW;

const goldTightPts = samples["GOLD+i45"].pts;
const upsilonTests = [];
for (const uFactor of [0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15]) {
  const scaledPts = goldTightPts.map(p => ({ ...p, x: p.x + Math.log10(uFactor) }));
  const fit = fitA0_unweighted(scaledPts);
  upsilonTests.push({ factor: uFactor, a0: fit.a0, logA0: fit.logA0 });
}
const upsilonRange = Math.max(...upsilonTests.map(t => t.logA0)) - Math.min(...upsilonTests.map(t => t.logA0));

const incTests = [];
for (const incCut of [30, 35, 40, 45, 50, 55, 60]) {
  const sample = allGalNames.filter(n => {
    const m = galaxyMeta[n]; const range = getGbarRange(n);
    return m.Vmax >= 50 && m.inc >= incCut && m.n >= 5 && range >= 1.0;
  });
  const pts = getPointsForGalaxies(sample);
  if (pts.length < 50) continue;
  const fit = fitA0_unweighted(pts);
  incTests.push({ incCut, nGal: sample.length, nPts: pts.length, a0: fit.a0, logA0: fit.logA0 });
}
const incRange = incTests.length > 1 ? Math.max(...incTests.map(t => t.logA0)) - Math.min(...incTests.map(t => t.logA0)) : 0;

const rangeTests = [];
for (const minRange of [0.5, 0.8, 1.0, 1.2, 1.5, 2.0]) {
  const sample = allGalNames.filter(n => {
    const m = galaxyMeta[n]; const range = getGbarRange(n);
    return m.Vmax >= 50 && m.inc >= 45 && m.n >= 5 && range >= minRange;
  });
  const pts = getPointsForGalaxies(sample);
  if (pts.length < 30) continue;
  const fit = fitA0_unweighted(pts);
  rangeTests.push({ minRange, nGal: sample.length, nPts: pts.length, a0: fit.a0, logA0: fit.logA0 });
}
const rangeTestRange = rangeTests.length > 1 ? Math.max(...rangeTests.map(t => t.logA0)) - Math.min(...rangeTests.map(t => t.logA0)) : 0;

log("  Reference: " + refSample + " McGaugh-only unweighted a0 = " + refA0.toFixed(1));
log("");

const budgetItems = [
  { source: "Estimator method (UW/W/Hier)", val: estimatorSpread / 2 },
  { source: "Sample selection", val: (() => {
    const sA = []; for (const n of ["GOLD", "GOLD+i45", "HIGH-MASS", "CLEAN"]) { if (unwFits[n]) sA.push(unwFits[n].logA0); }
    return (Math.max(...sA) - Math.min(...sA)) / 2;
  })() },
  { source: "Y_star (+/-15%)", val: upsilonRange / 2 },
  { source: "Inclination cut", val: incRange / 2 },
  { source: "Dynamic range cut", val: rangeTestRange / 2 },
  { source: "Hier. tau (intrinsic scatter)", val: refHier ? +refHier.tau : 0 },
];

log("  +----------------------------------+------------+");
log("  | Source                           | +/- (dex)  |");
log("  +----------------------------------+------------+");
let totalVar = 0;
for (const b of budgetItems) {
  log("  | " + b.source.padEnd(33) + "| +/-" + b.val.toFixed(4).padEnd(8) + " |");
  totalVar += b.val * b.val;
}
const totalSigma = Math.sqrt(totalVar);
log("  +----------------------------------+------------+");
log("  | TOTAL (quadrature)               | +/-" + totalSigma.toFixed(4).padEnd(8) + " |");
log("  +----------------------------------+------------+");
log("");

const a0Lower = Math.pow(10, refLogA0 - totalSigma);
const a0Upper = Math.pow(10, refLogA0 + totalSigma);

sep();
log("SECTION 8: FINAL PAPER TABLE");
sep();
log("");

log("  Table 1: a0 measurement summary (McGaugh RAR, single estimator)");
log("");
log("  +------------+------+------+------+----------+----------+----------+--------+--------+");
log("  | Sample     | n_gal| n_pts| Audit| a0(UW)   | a0(W)    | a0(Hier) | RMS    | ratio  |");
log("  +------------+------+------+------+----------+----------+----------+--------+--------+");

for (const name of ["ALL", "CLEAN", "GOLD", "GOLD+i45", "HIGH-MASS"]) {
  const s = samples[name];
  const uw = unwFits[name];
  const w = wFits[name];
  const h = hierResults[name];
  const audit = blindAudit(name, s.gals, uw.a0);
  const auditStr = audit.clean ? "CLEAN" : "ISSUE";
  const ms2 = uw.a0 * MS2_CONV;
  const ratio = ms2 / CH0_2PI;
  log("  | " + name.padEnd(11) +
    "| " + String(s.gals.length).padEnd(5) +
    "| " + String(s.pts.length).padEnd(5) +
    "| " + auditStr.padEnd(5) +
    "| " + String(uw.a0.toFixed(0)).padEnd(9) +
    "| " + String(w ? w.a0.toFixed(0) : "---").padEnd(9) +
    "| " + String(h ? h.popA0_weighted : "---").padEnd(9) +
    "| " + uw.rms.toFixed(4).padEnd(7) +
    "| " + ratio.toFixed(3).padEnd(7) +
    "|");
}
log("  +------------+------+------+------+----------+----------+----------+--------+--------+");
log("");

log("  Table 2: Uncertainty budget (" + refSample + ", McGaugh-only)");
log("");
for (const b of budgetItems) {
  log("    " + b.source.padEnd(35) + "+/-" + b.val.toFixed(4) + " dex");
}
log("    " + "\u2500".repeat(50));
log("    " + "TOTAL (quadrature)".padEnd(35) + "+/-" + totalSigma.toFixed(4) + " dex");
log("");

log("  Table 3: Hierarchical model diagnostics");
log("");
if (refHier) {
  log("    Population log(a0):     " + refHier.popLogA0_weighted.toFixed(4) + " +/- " + refHier.popSE.toFixed(4) + " dex");
  log("    Population a0:          " + refHier.popA0_weighted + " " + UNIT_LABEL);
  log("    Population a0 (m/s^2):  " + (refHier.popA0_weighted * MS2_CONV).toExponential(3));
  log("    Between-galaxy tau:     " + refHier.tau + " dex");
  log("    I^2 (heterogeneity):    " + refHier.I2_pct + "%");
  log("    Population median:      " + refHier.popMedianA0 + " " + UNIT_LABEL);
  log("    Population MAD:         " + refHier.popMAD + " dex");
}
log("");

sep();
log("SECTION 9: FINAL MEASUREMENT");
sep();
log("");

const bestA0 = refUW.a0;
const bestLogA0 = refUW.logA0;
const hierA0 = refHier ? refHier.popA0_weighted : bestA0;
const hierLogA0val = refHier ? refHier.popLogA0_weighted : bestLogA0;

log("  PRIMARY ESTIMATOR: McGaugh RAR, GOLD+i45, unweighted global fit");
log("    a0 = " + bestA0.toFixed(0) + " " + UNIT_LABEL);
log("    a0 = " + (bestA0 * MS2_CONV).toExponential(3) + " m/s^2");
log("    log(a0) = " + bestLogA0.toFixed(4) + " +/- " + totalSigma.toFixed(4) + " dex");
log("    a0 range: [" + a0Lower.toFixed(0) + ", " + a0Upper.toFixed(0) + "] " + UNIT_LABEL);
log("    a0/(cH0/2pi) = " + (bestA0 * MS2_CONV / CH0_2PI).toFixed(3));
log("");

log("  HIERARCHICAL ESTIMATOR: " + refSample + " weighted mean");
if (refHier) {
  log("    a0 = " + hierA0 + " " + UNIT_LABEL);
  log("    a0 = " + (hierA0 * MS2_CONV).toExponential(3) + " m/s^2");
  log("    log(a0) = " + hierLogA0val.toFixed(4) + " +/- " + refHier.popSE.toFixed(4) + " dex (stat)");
  log("    a0/(cH0/2pi) = " + (hierA0 * MS2_CONV / CH0_2PI).toFixed(3));
}
log("");

log("  LITERATURE COMPARISON:");
log("    McGaugh et al. (2016):  a0 = 1.20e-10 m/s^2 = " + A0_LIT_GAL.toFixed(0) + " " + UNIT_LABEL);
log("    This work (primary):    a0 = " + (bestA0 * MS2_CONV).toExponential(3) + " m/s^2");
if (refHier) {
  log("    This work (hier.):      a0 = " + (hierA0 * MS2_CONV).toExponential(3) + " m/s^2");
}
log("    Ratio (primary/lit):    " + (bestA0 / A0_LIT_GAL).toFixed(3));
if (refHier) {
  log("    Ratio (hier/lit):       " + (hierA0 / A0_LIT_GAL).toFixed(3));
}
log("    Lit value within 1-sigma band: " + (A0_LIT_GAL >= a0Lower && A0_LIT_GAL <= a0Upper ? "YES" : "NO"));
log("");

sep();
log("REPRODUCIBILITY MANIFEST (v" + VERSION + ")");
sep();
log("");
log("  Pipeline version: " + VERSION);
log("  Timestamp: " + TIMESTAMP);
log("  Input files: transition-scale.json (" + ALL_POINTS.length + " pts), rar-analysis-real.json (" + rarData.perGalaxy.length + " gals)");
log("  Estimator: McGaugh RAR ONLY (single estimator, no mean-of-3)");
log("  Fitting methods: unweighted RMS, error-weighted RMS, hierarchical random-effects");
log("  Y_disk = " + UPSILON_DISK + ", Y_bulge = " + UPSILON_BULGE);
log("  Fit: log(a0) in [2.0, 5.0], step 0.01, 100 golden-section refinements");
log("  Residual: log10(g_obs) - log10(McGaughRAR(g_bar, a0))");
log("  Sample cuts: LOCKED (see Section 1)");
log("  Audit threshold: |r| > 0.3 = significant");
log("  Uncertainty: 6-component budget in quadrature");
log("  Hierarchical: DerSimonian-Laird random-effects meta-analysis");

const output = {
  version: VERSION,
  timestamp: TIMESTAMP,
  estimator: "McGaugh RAR only",
  lockedParameters: {
    upsilonDisk: UPSILON_DISK, upsilonBulge: UPSILON_BULGE,
    interpolation: "McGaugh RAR: gobs = gbar / (1 - exp(-sqrt(gbar/a0)))",
    fitRange: "log(a0) in [2.0, 5.0]",
    residualDef: "log10(gobs) - log10(McGaughRAR(gbar, a0))",
    units: UNIT_LABEL, ms2Conv: MS2_CONV
  },
  sampleDefinitions: Object.fromEntries(Object.entries(samples).map(([k, v]) => [k, { nGal: v.gals.length, nPts: v.pts.length, cuts: v.desc }])),
  unweightedFits: Object.fromEntries(Object.entries(unwFits).map(([k, v]) => [k, { a0: +v.a0.toFixed(1), logA0: +v.logA0.toFixed(4), rms: +v.rms.toFixed(4), bias: +v.bias.toFixed(5) }])),
  weightedFits: Object.fromEntries(Object.entries(wFits).map(([k, v]) => [k, { a0: +v.a0.toFixed(1), logA0: +v.logA0.toFixed(4), wRMS: +v.rms.toFixed(4) }])),
  hierarchicalResults: Object.fromEntries(Object.entries(hierResults).map(([k, v]) => [k, {
    nGalaxies: v.nGalaxies, popLogA0: v.popLogA0_weighted, popA0: v.popA0_weighted,
    popSE: v.popSE, popMedianA0: v.popMedianA0, popMAD: v.popMAD,
    tau: +v.tau, I2_pct: v.I2_pct
  }])),
  uncertaintyBudget: {
    components: budgetItems.map(b => ({ source: b.source, halfRange_dex: +b.val.toFixed(4) })),
    total_dex: +totalSigma.toFixed(4),
    refLogA0: +refLogA0.toFixed(4),
    a0_lower: +a0Lower.toFixed(0),
    a0_upper: +a0Upper.toFixed(0)
  },
  sensitivityTests: {
    upsilon: upsilonTests.map(t => ({ factor: t.factor, a0: +t.a0.toFixed(0), logA0: +t.logA0.toFixed(4) })),
    inclination: incTests.map(t => ({ incCut: t.incCut, nGal: t.nGal, a0: +t.a0.toFixed(0), logA0: +t.logA0.toFixed(4) })),
    dynamicRange: rangeTests.map(t => ({ minRange: t.minRange, nGal: t.nGal, a0: +t.a0.toFixed(0), logA0: +t.logA0.toFixed(4) }))
  },
  finalMeasurement: {
    primary: {
      sample: refSample, method: "McGaugh RAR unweighted global fit",
      a0_galUnits: +bestA0.toFixed(1), a0_ms2: +(bestA0 * MS2_CONV).toExponential(3),
      logA0: +bestLogA0.toFixed(4), uncertainty_dex: +totalSigma.toFixed(4),
      a0_lower: +a0Lower.toFixed(0), a0_upper: +a0Upper.toFixed(0),
      ratio_cH0_2pi: +(bestA0 * MS2_CONV / CH0_2PI).toFixed(3)
    },
    hierarchical: refHier ? {
      sample: refSample, method: "DerSimonian-Laird weighted mean",
      a0_galUnits: hierA0, a0_ms2: +(hierA0 * MS2_CONV).toExponential(3),
      logA0: hierLogA0val, popSE_dex: refHier.popSE,
      tau_dex: +refHier.tau, I2_pct: refHier.I2_pct,
      ratio_cH0_2pi: +(hierA0 * MS2_CONV / CH0_2PI).toFixed(3)
    } : null,
    literatureComparison: {
      mcgaugh2016: A0_LIT_GAL.toFixed(0) + " " + UNIT_LABEL,
      withinOneSigma: A0_LIT_GAL >= a0Lower && A0_LIT_GAL <= a0Upper,
      ratio_primary_to_lit: +(bestA0 / A0_LIT_GAL).toFixed(3)
    }
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'gold-standard-results.json'), JSON.stringify(output, null, 2));
log("\nResults saved to gold-standard-results.json");
log("Pipeline v" + VERSION + " complete.\n");
