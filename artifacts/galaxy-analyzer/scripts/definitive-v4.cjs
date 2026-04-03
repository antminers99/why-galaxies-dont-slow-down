const fs = require('fs');
const path = require('path');

const VERSION = "4.0.0";
const TIMESTAMP = new Date().toISOString();

const UPSILON_FID = 0.5;
const MS2_CONV = 3.241e-14;
const C_MS = 299792458;
const H0_SI = 67.4e3 / 3.086e22;
const CH0_2PI = C_MS * H0_SI / (2 * Math.PI);
const CH0_2PI_GAL = CH0_2PI / MS2_CONV;
const A0_LIT = 1.2e-10;
const A0_LIT_GAL = A0_LIT / MS2_CONV;

const rarData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const tsData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'transition-scale.json'), 'utf8'));

function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(72)); }
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

function mad(arr) {
  const med = median(arr);
  return median(arr.map(v => Math.abs(v - med)));
}

function fitA0(logGbar, logGobs, a0Grid) {
  let bestA0 = 3000, bestChi2 = Infinity;
  for (const logA of a0Grid) {
    const a0 = Math.pow(10, logA);
    let chi2 = 0;
    for (let i = 0; i < logGbar.length; i++) {
      const gbar = Math.pow(10, logGbar[i]);
      const pred = mcgaughRAR(gbar, a0);
      if (!isFinite(pred) || pred <= 0) { chi2 += 100; continue; }
      const r = logGobs[i] - Math.log10(pred);
      chi2 += r * r;
    }
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

log("================================================================================");
log("  DEFINITIVE PIPELINE v4.0 — FULL SAMPLE, NUISANCE MARGINALIZATION");
log("================================================================================");
log("  Version: " + VERSION);
log("  Date: " + TIMESTAMP);
log("");

const galaxyMap = {};
for (const g of rarData.perGalaxy) {
  const Mgas = 1.33 * g.MHI;
  const Mstar = UPSILON_FID * g.L36;
  const fgas = (Mgas + Mstar > 0) ? Mgas / (Mgas + Mstar) : 0.5;
  galaxyMap[g.name] = {
    name: g.name,
    inc: g.inc,
    distance: g.distance,
    L36: g.L36,
    MHI: g.MHI,
    Vmax: g.Vmax,
    fgas: fgas,
    points: []
  };
}

for (const p of rarData.rarScatter) {
  if (galaxyMap[p.name]) {
    galaxyMap[p.name].points.push({
      logGbar: p.log_g_bar,
      logGobs: p.log_g_obs
    });
  }
}

const ltGalNames = new Set();
for (const p of tsData.plotPoints) {
  if (p.g.startsWith('LT_')) {
    ltGalNames.add(p.g);
    if (!galaxyMap[p.g]) {
      galaxyMap[p.g] = {
        name: p.g,
        inc: 60,
        distance: 5.0,
        L36: 0,
        MHI: 0,
        Vmax: 50,
        fgas: 0.5,
        points: [],
        isLT: true
      };
    }
    galaxyMap[p.g].points.push({
      logGbar: p.x,
      logGobs: p.x + p.y
    });
  }
}

const allGalaxies = Object.values(galaxyMap).filter(g => g.points.length >= 3);

let totalPts = 0;
for (const g of allGalaxies) totalPts += g.points.length;
log("  SPARC galaxies: " + (allGalaxies.length - ltGalNames.size));
log("  LITTLE THINGS galaxies: " + ltGalNames.size);
log("  Total galaxies: " + allGalaxies.length);
log("  Total radial points: " + totalPts);
log("  (SPARC: " + rarData.rarScatter.length + " full, not subsampled)");
log("");

sep();
log("  STEP 1: PER-GALAXY NUISANCE MARGINALIZATION");
sep();
log("");
log("  Nuisance parameters:");
log("    Upsilon_star: grid [0.20, 0.25, ..., 0.80] (13 values)");
log("    Prior: N(0.50, 0.12^2), truncated to [0.20, 0.80]");
log("    Delta_logD: grid [-0.10, -0.05, ..., +0.10] (5 values)");
log("    Prior: N(0.00, 0.07^2) [~17% distance uncertainty]");
log("    Grid: 13 x 5 = 65 evaluations per galaxy");
log("");

const upsilonGrid = [];
for (let u = 0.20; u <= 0.801; u += 0.05) upsilonGrid.push(+u.toFixed(2));
const dlogDGrid = [-0.10, -0.05, 0.00, 0.05, 0.10];

function gaussPrior(x, mu, sigma) {
  return Math.exp(-0.5 * ((x - mu) / sigma) ** 2);
}

const margResults = [];
let nMarg = 0;

for (const gal of allGalaxies) {
  if (gal.points.length < 5) {
    margResults.push({ name: gal.name, skip: true, reason: 'too few points' });
    continue;
  }

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

      const fit = fitA0(adjLogGbar, adjLogGobs, coarseGrid);

      if (fit.a0 > 1e5 || fit.a0 < 10) continue;

      const priorW = gaussPrior(upsilon, 0.50, 0.12) * gaussPrior(dlogD, 0.0, 0.07);

      gridResults.push({
        upsilon, dlogD,
        logA0: fit.logA0,
        a0: fit.a0,
        rms: fit.rms,
        chi2: fit.chi2,
        priorW
      });
    }
  }

  if (gridResults.length < 3) {
    margResults.push({ name: gal.name, skip: true, reason: 'no valid grid points' });
    continue;
  }

  const minChi2 = Math.min(...gridResults.map(r => r.chi2));
  const sigma2Hat = Math.max(minChi2 / (gal.points.length - 1), 1e-6);
  for (const r of gridResults) {
    r.weight = Math.exp(-0.5 * (r.chi2 - minChi2) / sigma2Hat) * r.priorW;
  }

  const sumW = gridResults.reduce((a, r) => a + r.weight, 0);
  const normW = gridResults.map(r => r.weight / sumW);

  let margLogA0 = 0;
  for (let i = 0; i < gridResults.length; i++) {
    margLogA0 += normW[i] * gridResults[i].logA0;
  }

  let margVar = 0;
  for (let i = 0; i < gridResults.length; i++) {
    margVar += normW[i] * (gridResults[i].logA0 - margLogA0) ** 2;
  }

  const bestRMS = gridResults.reduce((best, r) => r.chi2 < best.chi2 ? r : best, gridResults[0]).rms;
  const withinSE2 = (bestRMS * bestRMS) / gal.points.length;

  const totalSE2 = margVar + withinSE2;

  const gRange = Math.max(...baseLogGbar) - Math.min(...baseLogGbar);

  const fixedFit = fitA0(baseLogGbar, baseLogGobs, coarseGrid);

  margResults.push({
    name: gal.name,
    skip: false,
    n: gal.points.length,
    inc: gal.inc,
    distance: gal.distance,
    Vmax: gal.Vmax,
    fgas: gal.fgas,
    gRange,
    fixedLogA0: fixedFit.logA0,
    fixedA0: fixedFit.a0,
    fixedRMS: fixedFit.rms,
    margLogA0,
    margA0: Math.pow(10, margLogA0),
    margSE: Math.sqrt(totalSE2),
    margVar,
    withinSE2,
    nGrid: gridResults.length,
    isLT: !!gal.isLT
  });

  nMarg++;
}

log("  Marginalized: " + nMarg + " galaxies");
log("  Skipped: " + margResults.filter(r => r.skip).length + " galaxies");
log("");

const valid = margResults.filter(r => !r.skip);
const deltaLogA0 = valid.map(r => r.margLogA0 - r.fixedLogA0);
log("  Effect of marginalization on log(a0):");
log("    Mean shift: " + (deltaLogA0.reduce((a, b) => a + b, 0) / deltaLogA0.length).toFixed(4) + " dex");
log("    Median shift: " + median(deltaLogA0).toFixed(4) + " dex");
log("    RMS shift: " + Math.sqrt(deltaLogA0.reduce((a, v) => a + v * v, 0) / deltaLogA0.length).toFixed(4) + " dex");
log("    Max |shift|: " + Math.max(...deltaLogA0.map(Math.abs)).toFixed(4) + " dex");
log("");

sep();
log("  STEP 2: SAMPLE DEFINITIONS");
sep();
log("");

function selectSample(name, filter) {
  const gals = valid.filter(filter);
  return { name, gals, n: gals.length, nPts: gals.reduce((a, g) => a + g.n, 0) };
}

const samples = {
  ALL: selectSample("ALL", g => true),
  CLEAN: selectSample("CLEAN", g => g.Vmax >= 80 && g.inc >= 30 && g.n >= 5 && g.gRange >= 0.5),
  GOLD: selectSample("GOLD", g => g.Vmax >= 50 && g.inc >= 30 && g.n >= 5 && g.gRange >= 1.0),
  GOLDi45: selectSample("GOLD+i45", g => g.Vmax >= 50 && g.inc >= 45 && g.n >= 5 && g.gRange >= 1.0),
  HIMASS: selectSample("HIGH-MASS", g => g.Vmax > 150),
  LOMASS: selectSample("LOW-MASS", g => g.Vmax < 80)
};

log("  " + pad("Sample", 12) + padr("n_gal", 8) + padr("n_pts", 8));
log("  " + "\u2500".repeat(28));
for (const s of Object.values(samples)) {
  log("  " + pad(s.name, 12) + padr(s.n, 8) + padr(s.nPts, 8));
}
log("");

sep();
log("  STEP 3: HIERARCHICAL MODEL (DerSimonian-Laird)");
sep();
log("");

function hierarchicalDL(sampleGals, label) {
  const galEsts = sampleGals.filter(g =>
    g.margA0 > 10 && g.margA0 < 1e5 && g.gRange >= 0.3 && g.n >= 5
  );

  if (galEsts.length < 5) return null;

  const logA0s = galEsts.map(g => g.margLogA0);
  const se2s = galEsts.map(g => Math.max(g.margSE * g.margSE, 1e-6));

  const w_fe = se2s.map(s => 1.0 / s);
  const sumW_fe = w_fe.reduce((a, b) => a + b, 0);
  const mu_fe = logA0s.reduce((a, v, i) => a + w_fe[i] * v, 0) / sumW_fe;

  const Q = logA0s.reduce((a, v, i) => a + w_fe[i] * (v - mu_fe) ** 2, 0);
  const k = galEsts.length;
  const C = sumW_fe - w_fe.reduce((a, w) => a + w * w, 0) / sumW_fe;
  const tau2 = Math.max((Q - (k - 1)) / C, 0);
  const tau = Math.sqrt(tau2);

  const w_re = se2s.map(s => 1.0 / (s + tau2));
  const sumW_re = w_re.reduce((a, b) => a + b, 0);
  const mu_re = logA0s.reduce((a, v, i) => a + w_re[i] * v, 0) / sumW_re;
  const se_re = Math.sqrt(1.0 / sumW_re);

  const I2 = Q > (k - 1) ? ((Q - (k - 1)) / Q) * 100 : 0;

  const medLogA0 = median(logA0s);
  const madLogA0 = mad(logA0s);

  const a0_re = Math.pow(10, mu_re);
  const a0_re_ms2 = a0_re * MS2_CONV;
  const ratio = a0_re_ms2 / CH0_2PI;

  const lo1sig = Math.pow(10, mu_re - Math.sqrt(se_re * se_re + tau2));
  const hi1sig = Math.pow(10, mu_re + Math.sqrt(se_re * se_re + tau2));

  return {
    label,
    nGal: galEsts.length,
    mu: +mu_re.toFixed(4),
    se: +se_re.toFixed(4),
    tau: +tau.toFixed(4),
    tau2: +tau2.toFixed(6),
    I2: +I2.toFixed(1),
    a0: +a0_re.toFixed(0),
    a0_ms2: a0_re_ms2,
    ratio: +ratio.toFixed(4),
    median: +Math.pow(10, medLogA0).toFixed(0),
    medianLog: +medLogA0.toFixed(4),
    mad: +madLogA0.toFixed(4),
    lo1sig: +lo1sig.toFixed(0),
    hi1sig: +hi1sig.toFixed(0),
    galEsts
  };
}

const hierResults = {};
for (const [key, sample] of Object.entries(samples)) {
  const h = hierarchicalDL(sample.gals, sample.name);
  hierResults[key] = h;
  if (h) {
    log("  " + pad(h.label, 12) + " n=" + pad(h.nGal, 4)
      + " mu=" + pad(h.mu.toFixed(4), 8)
      + " a0=" + padr(h.a0, 6)
      + " tau=" + pad(h.tau.toFixed(3), 6)
      + " I2=" + padr(h.I2.toFixed(1) + "%", 7)
      + " ratio=" + h.ratio.toFixed(3));
  } else {
    log("  " + pad(sample.name, 12) + " — insufficient data for hierarchical model");
  }
}
log("");

sep();
log("  STEP 4: FIXED vs MARGINALIZED COMPARISON");
sep();
log("");

function fixedHierarchical(sampleGals, label) {
  const galEsts = sampleGals.filter(g =>
    g.fixedA0 > 10 && g.fixedA0 < 1e5 && g.gRange >= 0.3 && g.n >= 5
  );
  if (galEsts.length < 5) return null;

  const logA0s = galEsts.map(g => g.fixedLogA0);
  const se2s = galEsts.map(g => {
    const withinRMS = g.fixedRMS || 0.15;
    return Math.max((withinRMS * withinRMS) / g.n, 1e-6);
  });

  const w_fe = se2s.map(s => 1.0 / s);
  const sumW_fe = w_fe.reduce((a, b) => a + b, 0);
  const mu_fe = logA0s.reduce((a, v, i) => a + w_fe[i] * v, 0) / sumW_fe;

  const Q = logA0s.reduce((a, v, i) => a + w_fe[i] * (v - mu_fe) ** 2, 0);
  const k = galEsts.length;
  const C = sumW_fe - w_fe.reduce((a, w) => a + w * w, 0) / sumW_fe;
  const tau2 = Math.max((Q - (k - 1)) / C, 0);
  const tau = Math.sqrt(tau2);

  const w_re = se2s.map(s => 1.0 / (s + tau2));
  const sumW_re = w_re.reduce((a, b) => a + b, 0);
  const mu_re = logA0s.reduce((a, v, i) => a + w_re[i] * v, 0) / sumW_re;
  const se_re = Math.sqrt(1.0 / sumW_re);
  const I2 = Q > (k - 1) ? ((Q - (k - 1)) / Q) * 100 : 0;

  return {
    label,
    nGal: galEsts.length,
    mu: +mu_re.toFixed(4),
    a0: +Math.pow(10, mu_re).toFixed(0),
    se: +se_re.toFixed(4),
    tau: +tau.toFixed(4),
    I2: +I2.toFixed(1),
    ratio: +(Math.pow(10, mu_re) * MS2_CONV / CH0_2PI).toFixed(4)
  };
}

log("  " + pad("Sample", 12) + pad("Method", 14) + padr("a0", 7) + padr("tau", 7) + padr("I2", 7) + padr("ratio", 8));
log("  " + "\u2500".repeat(55));

for (const [key, sample] of Object.entries(samples)) {
  if (key === 'LOMASS') continue;
  const hMarg = hierResults[key];
  const hFixed = fixedHierarchical(sample.gals, sample.name);
  if (hMarg) {
    log("  " + pad(sample.name, 12) + pad("marginalized", 14)
      + padr(hMarg.a0, 7) + padr(hMarg.tau.toFixed(3), 7)
      + padr(hMarg.I2.toFixed(1) + "%", 7) + padr(hMarg.ratio.toFixed(3), 8));
  }
  if (hFixed) {
    log("  " + pad("", 12) + pad("fixed Y*", 14)
      + padr(hFixed.a0, 7) + padr(hFixed.tau.toFixed(3), 7)
      + padr(hFixed.I2.toFixed(1) + "%", 7) + padr(hFixed.ratio.toFixed(3), 8));
  }
  if (hMarg && hFixed) {
    const delta = ((hMarg.a0 - hFixed.a0) / hFixed.a0 * 100).toFixed(1);
    const tauDelta = (hMarg.tau - hFixed.tau).toFixed(3);
    log("  " + pad("", 12) + pad("delta", 14)
      + padr(delta + "%", 7) + padr(tauDelta, 7) + padr("", 7) + padr("", 8));
  }
  log("");
}

sep();
log("  STEP 5: RATIO STABILITY TEST — a0/(cH0/2pi)");
sep();
log("");
log("  Question: Is the ratio a0/(cH0/2pi) stable across subsamples?");
log("  If yes: suggestive of cosmological connection.");
log("  If no: ratio is method-dependent, no claim possible.");
log("");

const gi45 = hierResults.GOLDi45;
if (gi45) {
  const gi45Gals = gi45.galEsts;

  const splits = [
    { name: "ALL valid", filter: g => true },
    { name: "Vmax>150", filter: g => g.Vmax > 150 },
    { name: "Vmax 80-150", filter: g => g.Vmax >= 80 && g.Vmax <= 150 },
    { name: "D < 15 Mpc", filter: g => g.distance < 15 },
    { name: "D >= 15 Mpc", filter: g => g.distance >= 15 },
    { name: "f_gas < 0.3", filter: g => g.fgas < 0.3 },
    { name: "f_gas >= 0.3", filter: g => g.fgas >= 0.3 },
    { name: "inc >= 60", filter: g => g.inc >= 60 },
    { name: "inc 45-60", filter: g => g.inc >= 45 && g.inc < 60 },
    { name: "n >= 20 pts", filter: g => g.n >= 20 },
    { name: "n < 20 pts", filter: g => g.n >= 5 && g.n < 20 },
  ];

  log("  GOLD+i45 subsample splits (hierarchical marginalized):");
  log("");
  log("  " + pad("Split", 16) + padr("n", 5) + padr("a0", 7) + padr("tau", 7) + padr("ratio", 8) + padr("SE", 7));
  log("  " + "\u2500".repeat(50));

  const ratios = [];

  for (const split of splits) {
    const subset = gi45Gals.filter(split.filter);
    if (subset.length < 5) {
      log("  " + pad(split.name, 16) + padr(subset.length, 5) + "  (too few)");
      continue;
    }

    const h = hierarchicalDL(subset, split.name);
    if (h) {
      log("  " + pad(split.name, 16) + padr(h.nGal, 5)
        + padr(h.a0, 7) + padr(h.tau.toFixed(3), 7)
        + padr(h.ratio.toFixed(3), 8) + padr(h.se.toFixed(3), 7));
      ratios.push({ name: split.name, ratio: h.ratio, a0: h.a0, n: h.nGal });
    }
  }

  log("");
  if (ratios.length >= 2) {
    const rs = ratios.map(r => r.ratio);
    const minR = Math.min(...rs), maxR = Math.max(...rs);
    const rangeR = maxR - minR;
    const meanR = rs.reduce((a, b) => a + b, 0) / rs.length;
    const stdR = Math.sqrt(rs.reduce((a, v) => a + (v - meanR) ** 2, 0) / (rs.length - 1));

    log("  RATIO STABILITY SUMMARY:");
    log("    Splits tested: " + ratios.length);
    log("    Ratio range: " + minR.toFixed(3) + " - " + maxR.toFixed(3));
    log("    Ratio spread: " + rangeR.toFixed(3));
    log("    Mean ratio: " + meanR.toFixed(3));
    log("    Std dev: " + stdR.toFixed(3));
    log("    CV (std/mean): " + (stdR / meanR * 100).toFixed(1) + "%");
    log("");

    if (rangeR < 0.3 && stdR < 0.15) {
      log("    VERDICT: Ratio is REASONABLY STABLE across splits.");
      log("    Suggestive of real physical connection, but substantial");
      log("    between-galaxy heterogeneity limits precision.");
    } else if (rangeR < 0.6) {
      log("    VERDICT: Ratio shows MODERATE variation across splits.");
      log("    Connection suggestive but not compelling.");
    } else {
      log("    VERDICT: Ratio is UNSTABLE across splits.");
      log("    No claim about cosmological connection defensible.");
    }
  }
}

log("");
sep();
log("  STEP 6: CROSS-SAMPLE HIERARCHICAL COMPARISON");
sep();
log("");

log("  " + pad("Sample", 12) + padr("n", 5) + padr("a0", 7) + padr("SE", 6)
  + padr("tau", 7) + padr("I2", 7) + padr("ratio", 8)
  + padr("[lo, hi]", 16));
log("  " + "\u2500".repeat(68));

for (const key of ['ALL', 'CLEAN', 'GOLD', 'GOLDi45', 'HIMASS']) {
  const h = hierResults[key];
  if (!h) continue;
  log("  " + pad(h.label, 12)
    + padr(h.nGal, 5) + padr(h.a0, 7) + padr(h.se.toFixed(3), 6)
    + padr(h.tau.toFixed(3), 7) + padr(h.I2.toFixed(1) + "%", 7)
    + padr(h.ratio.toFixed(3), 8)
    + padr("[" + h.lo1sig + ", " + h.hi1sig + "]", 16));
}
log("");
log("  Literature a0 = " + A0_LIT_GAL.toFixed(0) + " (km/s)^2/kpc = 1.20e-10 m/s^2");
log("  cH0/2pi = " + CH0_2PI_GAL.toFixed(0) + " (km/s)^2/kpc = " + CH0_2PI.toExponential(3) + " m/s^2");
log("");

sep();
log("  STEP 7: FINAL POSTERIOR");
sep();
log("");

const headline = hierResults.GOLDi45;
if (headline) {
  const totalVar = headline.se * headline.se + headline.tau * headline.tau;
  const totalSD = Math.sqrt(totalVar);
  const lo68 = Math.pow(10, headline.mu - totalSD);
  const hi68 = Math.pow(10, headline.mu + totalSD);
  const lo95 = Math.pow(10, headline.mu - 2 * totalSD);
  const hi95 = Math.pow(10, headline.mu + 2 * totalSD);

  log("  ╔══════════════════════════════════════════════════════════════════╗");
  log("  ║  DEFINITIVE a0 MEASUREMENT (v4.0, marginalized hierarchical)   ║");
  log("  ╠══════════════════════════════════════════════════════════════════╣");
  log("  ║  Sample: GOLD+i45 (" + headline.nGal + " galaxies)");
  log("  ║  Data: " + totalPts + " radial points (full, not subsampled)");
  log("  ║  Estimator: McGaugh RAR, per-galaxy Y* marginalization");
  log("  ║  Model: DerSimonian-Laird random-effects meta-analysis");
  log("  ║");
  log("  ║  HEADLINE:");
  log("  ║    a0 = " + headline.a0 + " (km/s)^2/kpc");
  log("  ║    a0 = " + headline.a0_ms2.toExponential(3) + " m/s^2");
  log("  ║    log(a0) = " + headline.mu.toFixed(4) + " +/- " + headline.se.toFixed(4) + " dex (stat)");
  log("  ║");
  log("  ║  POPULATION SCATTER:");
  log("  ║    tau = " + headline.tau.toFixed(4) + " dex (between-galaxy)");
  log("  ║    I^2 = " + headline.I2.toFixed(1) + "%");
  log("  ║");
  log("  ║  CREDIBLE INTERVALS:");
  log("  ║    68%: [" + lo68.toFixed(0) + ", " + hi68.toFixed(0) + "] (km/s)^2/kpc");
  log("  ║    95%: [" + lo95.toFixed(0) + ", " + hi95.toFixed(0) + "] (km/s)^2/kpc");
  log("  ║");
  log("  ║  COSMOLOGICAL RATIO:");
  log("  ║    a0/(cH0/2pi) = " + headline.ratio.toFixed(4));
  log("  ║");
  log("  ║  LITERATURE:");
  log("  ║    McGaugh+2016: a0 = 3703 (km/s)^2/kpc = 1.20e-10 m/s^2");
  log("  ║    Within 1-sigma: " + (A0_LIT_GAL >= lo68 && A0_LIT_GAL <= hi68 ? "YES" : "NO"));
  log("  ║");
  log("  ║  MEDIAN (robust):");
  log("  ║    a0 = " + headline.median + " (km/s)^2/kpc (MAD = " + headline.mad.toFixed(3) + " dex)");
  log("  ║");
  log("  ║  NOT CLAIMED: universal exact constant, dark matter solved,");
  log("  ║  MOND proved, cosmological origin established.");
  log("  ╚══════════════════════════════════════════════════════════════════╝");
  log("");

  const v3hier = { a0: 3374, tau: 0.245, I2: 89.4, ratio: 1.049 };
  log("  COMPARISON v3.0 (fixed Y*) vs v4.0 (marginalized):");
  log("    v3.0: a0=" + v3hier.a0 + ", tau=" + v3hier.tau + ", I2=" + v3hier.I2 + "%, ratio=" + v3hier.ratio);
  log("    v4.0: a0=" + headline.a0 + ", tau=" + headline.tau.toFixed(3) + ", I2=" + headline.I2.toFixed(1) + "%, ratio=" + headline.ratio.toFixed(3));
  const dA0 = ((headline.a0 - v3hier.a0) / v3hier.a0 * 100).toFixed(1);
  const dTau = (headline.tau - v3hier.tau).toFixed(3);
  log("    Delta a0: " + dA0 + "%");
  log("    Delta tau: " + dTau + " dex");
  log("    Interpretation: " + (Math.abs(+dTau) < 0.05
    ? "Y* marginalization has MINIMAL effect on tau — population heterogeneity is real, not Y*-driven."
    : Math.abs(+dTau) > 0.05
      ? "Y* marginalization REDUCES tau — some between-galaxy scatter was Y*-driven."
      : ""));

  const output = {
    version: VERSION,
    timestamp: TIMESTAMP,
    description: "Definitive pipeline: full sample, per-galaxy Y* marginalization, hierarchical DL model",
    data: {
      sparcGalaxies: allGalaxies.length - ltGalNames.size,
      ltGalaxies: ltGalNames.size,
      totalPoints: totalPts,
      sparcPoints: rarData.rarScatter.length
    },
    nuisanceGrid: {
      upsilon: upsilonGrid,
      deltaLogD: dlogDGrid,
      upsilonPrior: "N(0.50, 0.12^2)",
      distancePrior: "N(0.00, 0.07^2)"
    },
    samples: {},
    headline: {
      sample: "GOLD+i45",
      nGal: headline.nGal,
      a0_gal: headline.a0,
      a0_ms2: headline.a0_ms2,
      logA0: headline.mu,
      logA0_SE: headline.se,
      tau: headline.tau,
      I2: headline.I2,
      ratio_cH02pi: headline.ratio,
      CI68: [+lo68.toFixed(0), +hi68.toFixed(0)],
      CI95: [+lo95.toFixed(0), +hi95.toFixed(0)],
      median: headline.median,
      medianLog: headline.medianLog,
      mad: headline.mad
    },
    ratioStability: {},
    comparison_v3: {
      v3_a0: v3hier.a0,
      v4_a0: headline.a0,
      delta_pct: +dA0,
      v3_tau: v3hier.tau,
      v4_tau: +headline.tau.toFixed(3),
      delta_tau: +dTau
    }
  };

  for (const [key, h] of Object.entries(hierResults)) {
    if (!h) continue;
    output.samples[key] = {
      label: h.label,
      nGal: h.nGal,
      a0: h.a0,
      logA0: h.mu,
      se: h.se,
      tau: h.tau,
      I2: h.I2,
      ratio: h.ratio,
      lo1sig: h.lo1sig,
      hi1sig: h.hi1sig
    };
  }

  if (gi45) {
    const gi45Gals = gi45.galEsts;
    const splitTests = [
      { name: "ALL_valid", filter: g => true },
      { name: "Vmax_gt150", filter: g => g.Vmax > 150 },
      { name: "Vmax_80_150", filter: g => g.Vmax >= 80 && g.Vmax <= 150 },
      { name: "D_lt15", filter: g => g.distance < 15 },
      { name: "D_ge15", filter: g => g.distance >= 15 },
      { name: "fgas_lt0.3", filter: g => g.fgas < 0.3 },
      { name: "fgas_ge0.3", filter: g => g.fgas >= 0.3 },
      { name: "inc_ge60", filter: g => g.inc >= 60 },
      { name: "inc_45_60", filter: g => g.inc >= 45 && g.inc < 60 },
      { name: "n_ge20", filter: g => g.n >= 20 },
      { name: "n_lt20", filter: g => g.n >= 5 && g.n < 20 },
    ];

    for (const split of splitTests) {
      const subset = gi45Gals.filter(split.filter);
      if (subset.length < 5) continue;
      const h = hierarchicalDL(subset, split.name);
      if (h) {
        output.ratioStability[split.name] = {
          n: h.nGal, a0: h.a0, tau: h.tau, ratio: h.ratio, se: h.se
        };
      }
    }
  }

  const outPath = path.join(__dirname, '..', 'public', 'definitive-v4-results.json');
  fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
  log("");
  log("  Results saved to: public/definitive-v4-results.json");
}

log("");
log("================================================================================");
log("  PIPELINE COMPLETE");
log("================================================================================");
