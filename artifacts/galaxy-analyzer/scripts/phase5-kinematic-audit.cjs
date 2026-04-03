const fs = require('fs');
const path = require('path');

const VERSION = "5.0.0";
const TIMESTAMP = new Date().toISOString();

const UPSILON_FID = 0.5;
const MS2_CONV = 3.241e-14;
const C_MS = 299792458;
const H0_SI = 67.4e3 / 3.086e22;
const CH0_2PI = C_MS * H0_SI / (2 * Math.PI);
const CH0_2PI_GAL = CH0_2PI / MS2_CONV;
const SIGMA_HI = 10.0;

const rarData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const tsData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'transition-scale.json'), 'utf8'));

function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(72)); }

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

  const gridResults = [];
  for (const upsilon of upsilonGrid) {
    const scale = upsilon / UPSILON_FID;
    const gbarAdj = (gal.isLT || gal.fgas >= 0.99)
      ? 0.0
      : Math.log10(scale * (1 - gal.fgas) + gal.fgas);
    const adjLogGbar = gal.points.map(p => p.logGbar + gbarAdj);

    for (const dlogD of dlogDGrid) {
      const adjLogGobs = gal.points.map(p => p.logGobs - dlogD);
      const fit = fitA0(adjLogGbar, adjLogGobs);
      if (!fit || fit.a0 > 1e5 || fit.a0 < 10) continue;

      const priorW = gaussPrior(upsilon, 0.50, 0.12) * gaussPrior(dlogD, 0.0, 0.07);
      gridResults.push({ ...fit, upsilon, dlogD, priorW });
    }
  }

  if (gridResults.length === 0) return null;

  const minChi2 = Math.min(...gridResults.map(r => r.chi2));
  const n = gal.points.length;
  const bestRMS = Math.sqrt(minChi2 / n);
  const withinSE2 = bestRMS * bestRMS / n;

  let sumW = 0, sumWlogA = 0, sumWlogA2 = 0;
  for (const r of gridResults) {
    const relChi2 = r.chi2 - minChi2;
    const w = Math.exp(-0.5 * relChi2 / withinSE2) * r.priorW;
    if (!isFinite(w)) continue;
    sumW += w;
    sumWlogA += w * r.logA0;
    sumWlogA2 += w * r.logA0 * r.logA0;
  }

  if (sumW <= 0) return null;
  const meanLogA = sumWlogA / sumW;
  const varLogA = sumWlogA2 / sumW - meanLogA * meanLogA;
  const seLogA = Math.sqrt(Math.max(varLogA, 1e-8));

  return {
    name: gal.name,
    logA0: meanLogA,
    a0: Math.pow(10, meanLogA),
    se: seLogA,
    n: gal.points.length,
    Vmax: gal.Vmax,
    distance: gal.distance,
    inc: gal.inc,
    fgas: gal.fgas,
    Rdisk: gal.Rdisk || 1.0,
    Q_kin: gal.Q_kin || 0.7,
    eta_rot: gal.eta_rot || 0,
    S_out: gal.S_out || 1.0,
    isLT: !!gal.isLT
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
  return { label, nGal: k, a0: Math.round(a0), logA0: +muRE.toFixed(4), se: +seRE.toFixed(4), tau: +tau.toFixed(4), I2: +I2.toFixed(1), ratio: +ratio.toFixed(4) };
}

log("================================================================================");
log("  PHASE 5: KINEMATIC CONTAMINATION AUDIT");
log("================================================================================");
log("  Version: " + VERSION);
log("  Date: " + TIMESTAMP);
log("  Purpose: Test whether a0 survives non-circular motion and pressure");
log("           support corrections — the two remaining unchecked systematics.");
log("");

const galaxyMap = {};
for (const g of rarData.perGalaxy) {
  const Mgas = 1.33 * g.MHI;
  const Mstar = UPSILON_FID * g.L36;
  const fgas = (Mgas + Mstar > 0) ? Mgas / (Mgas + Mstar) : 0.5;
  galaxyMap[g.name] = {
    name: g.name, inc: g.inc, distance: g.distance,
    L36: g.L36, MHI: g.MHI, Vmax: g.Vmax, fgas: fgas,
    Rdisk: g.Rdisk || 1.0,
    Q_kin: g.Q_kin, eta_rot: g.eta_rot, eta_bar: g.eta_bar,
    S_out: g.S_out, sigma_bar: g.sigma_bar,
    points: [], pointsFull: []
  };
}

for (const p of rarData.rarScatter) {
  if (galaxyMap[p.name]) {
    galaxyMap[p.name].points.push({
      logGbar: p.log_g_bar, logGobs: p.log_g_obs,
      isInner: p.isInner, r: p.r, Q_kin: p.Q_kin
    });
    galaxyMap[p.name].pointsFull.push(p);
  }
}

for (const p of tsData.plotPoints) {
  if (p.g.startsWith('LT_')) {
    if (!galaxyMap[p.g]) {
      galaxyMap[p.g] = {
        name: p.g, inc: 60, distance: 5.0, L36: 0, MHI: 0,
        Vmax: 50, fgas: 0.5, Rdisk: 1.0, Q_kin: 0.7,
        eta_rot: 0, eta_bar: 0, S_out: 1.0, sigma_bar: 0,
        points: [], pointsFull: [], isLT: true
      };
    }
    galaxyMap[p.g].points.push({
      logGbar: p.x, logGobs: p.x + p.y,
      isInner: false, r: 0, Q_kin: 0.7
    });
  }
}

const allGalaxies = Object.values(galaxyMap).filter(g => g.points.length >= 5);
const sparcGals = allGalaxies.filter(g => !g.isLT);
const goldI45 = sparcGals.filter(g => g.Q_kin >= 0.7 && g.inc >= 45);

log("  Total galaxies (n>=5): " + allGalaxies.length);
log("  SPARC galaxies: " + sparcGals.length);
log("  GOLD+i45 baseline: " + goldI45.length);
log("");

sep();
log("  TEST 1: NON-CIRCULAR MOTION PROXIES");
sep();
log("");
log("  Without explicit bar/morphology classification, we use available");
log("  kinematic proxies from the SPARC data:");
log("    Q_kin  — kinematic data quality indicator");
log("    isInner — inner vs outer radial points");
log("    eta_rot — rotation curve quality metric");
log("    S_out  — outer slope symmetry");
log("");

log("  TEST 1A: Inner vs Outer radii (key non-circular motion diagnostic)");
log("  Non-circular motions (bars, spirals) are strongest in inner regions.");
log("  If they bias a0, inner-only and outer-only fits should disagree.");
log("");

function fitSubsetPoints(galaxies, pointFilter, label) {
  const galEsts = [];
  for (const gal of galaxies) {
    const pts = gal.points.filter(pointFilter);
    if (pts.length < 5) continue;
    const subGal = { ...gal, points: pts };
    const est = marginalizeSingleGalaxy(subGal);
    if (est) galEsts.push(est);
  }
  if (galEsts.length < 5) return null;
  return hierarchicalDL(galEsts, label);
}

const innerResult = fitSubsetPoints(goldI45, p => p.isInner, "Inner only");
const outerResult = fitSubsetPoints(goldI45, p => !p.isInner, "Outer only");

const margGoldI45 = [];
for (const gal of goldI45) {
  const est = marginalizeSingleGalaxy(gal);
  if (est) margGoldI45.push(est);
}
const baseline = hierarchicalDL(margGoldI45, "GOLD+i45 baseline");

function printResult(r) {
  if (!r) { log("    (insufficient data)"); return; }
  log("    " + r.label + ": n=" + r.nGal + " gal, a0=" + r.a0 +
    " (km/s)^2/kpc, tau=" + r.tau + ", ratio=" + r.ratio);
}

log("  Results (per-galaxy marginalized, then hierarchical DL):");
printResult(baseline);
printResult(innerResult);
printResult(outerResult);

if (baseline && innerResult && outerResult) {
  const innerDelta = Math.abs(innerResult.logA0 - baseline.logA0);
  const outerDelta = Math.abs(outerResult.logA0 - baseline.logA0);
  const innerOuterDelta = Math.abs(innerResult.logA0 - outerResult.logA0);
  log("");
  log("  Inner-Outer split: " + (innerOuterDelta * 100 / baseline.logA0).toFixed(1) + "% difference");
  log("  Delta log(a0): inner-outer = " + innerOuterDelta.toFixed(3) + " dex");
  log("  Interpretation: " + (innerOuterDelta < 0.1 ? "SMALL — non-circular motions not dominant" :
    innerOuterDelta < 0.2 ? "MODERATE — some non-circular contamination" :
      "LARGE — non-circular motions are a major concern"));
}

log("");
log("  TEST 1B: Kinematic quality Q_kin split");
log("  Higher Q_kin = better kinematic data, less contamination.");
log("");

const highQ = margGoldI45.filter(g => g.Q_kin >= 0.9);
const lowQ = margGoldI45.filter(g => g.Q_kin < 0.9);
const hqResult = hierarchicalDL(highQ, "Q_kin >= 0.9");
const lqResult = hierarchicalDL(lowQ, "Q_kin < 0.9");

printResult(hqResult);
printResult(lqResult);

if (hqResult && lqResult) {
  const qDelta = Math.abs(hqResult.logA0 - lqResult.logA0);
  log("  Q_kin split delta: " + qDelta.toFixed(3) + " dex (" +
    (Math.abs(hqResult.a0 - lqResult.a0) / baseline.a0 * 100).toFixed(1) + "% of baseline)");
}

log("");
log("  TEST 1C: Rotation quality eta_rot split");
log("");

const medEtaRot = median(margGoldI45.map(g => g.eta_rot));
const highEta = margGoldI45.filter(g => g.eta_rot >= medEtaRot);
const lowEta = margGoldI45.filter(g => g.eta_rot < medEtaRot);
const heResult = hierarchicalDL(highEta, "eta_rot >= median");
const leResult = hierarchicalDL(lowEta, "eta_rot < median");

printResult(heResult);
printResult(leResult);

if (heResult && leResult) {
  const eDelta = Math.abs(heResult.logA0 - leResult.logA0);
  log("  eta_rot split delta: " + eDelta.toFixed(3) + " dex");
}

log("");
log("  TEST 1D: Outer slope symmetry S_out split");
log("");

const medSout = median(sparcGals.map(g => g.S_out));
const highSout = margGoldI45.filter(g => g.S_out >= medSout);
const lowSout = margGoldI45.filter(g => g.S_out < medSout);
const hsResult = hierarchicalDL(highSout, "S_out >= median (symmetric)");
const lsResult = hierarchicalDL(lowSout, "S_out < median (asymmetric)");

printResult(hsResult);
printResult(lsResult);

if (hsResult && lsResult) {
  const sDelta = Math.abs(hsResult.logA0 - lsResult.logA0);
  log("  S_out split delta: " + sDelta.toFixed(3) + " dex");
  log("  Interpretation: " + (sDelta < 0.1 ? "SMALL — asymmetry not biasing a0" :
    "NOTABLE — asymmetric RCs yield different a0"));
}

sep();
log("");
sep();
log("  TEST 2: PRESSURE SUPPORT / ASYMMETRIC DRIFT CORRECTION");
sep();
log("");
log("  Applying first-order asymmetric drift correction:");
log("  g_obs_corrected = g_obs + sigma_HI^2 / h_gas");
log("  where sigma_HI = " + SIGMA_HI + " km/s (standard HI dispersion)");
log("  and h_gas = 2 * Rdisk (gas scale length ~ 2x stellar scale length)");
log("");
log("  This correction is largest for low-mass galaxies where");
log("  sigma^2/h_gas is comparable to g_obs.");
log("");

function applyPressureCorrection(gal) {
  const h_gas = 2.0 * (gal.Rdisk || 1.0);
  const correctionGal = SIGMA_HI * SIGMA_HI / h_gas;

  const correctedPoints = gal.points.map(p => {
    const gobs = Math.pow(10, p.logGobs);
    const gobs_corr = gobs + correctionGal;
    return {
      ...p,
      logGobs: Math.log10(gobs_corr),
      logGobs_orig: p.logGobs,
      correction_dex: Math.log10(gobs_corr) - p.logGobs
    };
  });

  return { ...gal, points: correctedPoints, pressureCorrGal: correctionGal };
}

const correctedGalaxies = goldI45.map(g => applyPressureCorrection(g));

const lowVmax = correctedGalaxies.filter(g => g.Vmax < 80);
const highVmax = correctedGalaxies.filter(g => g.Vmax >= 80);
const dwarfs = correctedGalaxies.filter(g => g.fgas >= 0.5);

log("  Correction magnitudes (sample statistics):");
const allCorr = correctedGalaxies.flatMap(g => g.points.map(p => p.correction_dex));
log("    Mean correction: " + (allCorr.reduce((a, b) => a + b, 0) / allCorr.length).toFixed(4) + " dex");
log("    Median correction: " + median(allCorr).toFixed(4) + " dex");
log("    Max correction: " + Math.max(...allCorr).toFixed(4) + " dex");
log("");

const lowVmaxCorr = correctedGalaxies.filter(g => g.Vmax < 80).flatMap(g => g.points.map(p => p.correction_dex));
const hiVmaxCorr = correctedGalaxies.filter(g => g.Vmax >= 80).flatMap(g => g.points.map(p => p.correction_dex));
if (lowVmaxCorr.length > 0) log("    Mean correction (Vmax<80): " + (lowVmaxCorr.reduce((a, b) => a + b, 0) / lowVmaxCorr.length).toFixed(4) + " dex");
if (hiVmaxCorr.length > 0) log("    Mean correction (Vmax>=80): " + (hiVmaxCorr.reduce((a, b) => a + b, 0) / hiVmaxCorr.length).toFixed(4) + " dex");
log("");

log("  TEST 2A: Full GOLD+i45 before vs after pressure correction");
log("");

const margCorrected = [];
for (const gal of correctedGalaxies) {
  const est = marginalizeSingleGalaxy(gal);
  if (est) margCorrected.push(est);
}
const corrBaseline = hierarchicalDL(margCorrected, "GOLD+i45 pressure-corrected");

log("  Before correction:");
printResult(baseline);
log("  After pressure correction:");
printResult(corrBaseline);

if (baseline && corrBaseline) {
  const shift = corrBaseline.logA0 - baseline.logA0;
  const shiftPct = (Math.pow(10, corrBaseline.logA0) / Math.pow(10, baseline.logA0) - 1) * 100;
  const tauChange = corrBaseline.tau - baseline.tau;
  log("");
  log("  Shift: " + (shiftPct > 0 ? "+" : "") + shiftPct.toFixed(1) + "% (" +
    (shift > 0 ? "+" : "") + shift.toFixed(4) + " dex)");
  log("  Tau change: " + (tauChange > 0 ? "+" : "") + tauChange.toFixed(4) + " dex");
  log("  Interpretation: " + (Math.abs(shiftPct) < 5 ? "SMALL shift" :
    Math.abs(shiftPct) < 15 ? "MODERATE shift" : "LARGE shift") +
    "; tau " + (tauChange < -0.02 ? "DECREASED (correction helps)" :
      tauChange > 0.02 ? "INCREASED (correction hurts)" : "STABLE"));
}

log("");
log("  TEST 2B: Pressure correction impact by mass regime");
log("");

const corrLowVmax = margCorrected.filter(g => g.Vmax < 80);
const corrMidVmax = margCorrected.filter(g => g.Vmax >= 80 && g.Vmax <= 150);
const corrHighVmax = margCorrected.filter(g => g.Vmax > 150);

const clv = hierarchicalDL(corrLowVmax, "Vmax<80 corrected");
const cmv = hierarchicalDL(corrMidVmax, "Vmax 80-150 corrected");
const chv = hierarchicalDL(corrHighVmax, "Vmax>150 corrected");

printResult(clv);
printResult(cmv);
printResult(chv);

log("");
log("  TEST 2C: Pressure correction impact on gas-rich galaxies");
log("");

const corrHighFgas = margCorrected.filter(g => g.fgas >= 0.3);
const corrLowFgas = margCorrected.filter(g => g.fgas < 0.3);
const cfh = hierarchicalDL(corrHighFgas, "fgas>=0.3 corrected");
const cfl = hierarchicalDL(corrLowFgas, "fgas<0.3 corrected");

printResult(cfh);
printResult(cfl);

sep();
log("");
sep();
log("  TEST 3: RED FLAG RE-EXAMINATION AFTER CORRECTIONS");
sep();
log("");
log("  Re-testing the three v4.0 red-flag splits with pressure correction applied.");
log("");

log("  RED FLAG 1: Distance split (v4.0: 35% difference)");
const corrDlt15 = margCorrected.filter(g => g.distance < 15);
const corrDge15 = margCorrected.filter(g => g.distance >= 15);
const cdl = hierarchicalDL(corrDlt15, "D<15 corrected");
const cdh = hierarchicalDL(corrDge15, "D>=15 corrected");
printResult(cdl);
printResult(cdh);
if (cdl && cdh) {
  const dPct = Math.abs(cdl.a0 - cdh.a0) / baseline.a0 * 100;
  log("  Distance split: " + dPct.toFixed(1) + "% (was 35% in v4.0)");
  log("  " + (dPct < 20 ? "IMPROVED" : dPct < 35 ? "SLIGHTLY IMPROVED" : "NO IMPROVEMENT"));
}

log("");
log("  RED FLAG 2: Inclination split (v4.0: 42% difference)");
const corrIge60 = margCorrected.filter(g => g.inc >= 60);
const corrI45_60 = margCorrected.filter(g => g.inc >= 45 && g.inc < 60);
const cih = hierarchicalDL(corrIge60, "inc>=60 corrected");
const cil = hierarchicalDL(corrI45_60, "inc 45-60 corrected");
printResult(cih);
printResult(cil);
if (cih && cil) {
  const iPct = Math.abs(cih.a0 - cil.a0) / baseline.a0 * 100;
  log("  Inclination split: " + iPct.toFixed(1) + "% (was 42% in v4.0)");
  log("  " + (iPct < 20 ? "IMPROVED" : iPct < 35 ? "SLIGHTLY IMPROVED" : "NO IMPROVEMENT"));
}

log("");
log("  RED FLAG 3: Velocity split (v4.0: 36% difference)");
const corrVhi = margCorrected.filter(g => g.Vmax > 150);
const corrVmid = margCorrected.filter(g => g.Vmax >= 80 && g.Vmax <= 150);
printResult(chv);
printResult(cmv);
if (chv && cmv) {
  const vPct = Math.abs(chv.a0 - cmv.a0) / baseline.a0 * 100;
  log("  Velocity split: " + vPct.toFixed(1) + "% (was 36% in v4.0)");
  log("  " + (vPct < 20 ? "IMPROVED" : vPct < 35 ? "SLIGHTLY IMPROVED" : "NO IMPROVEMENT"));
}

sep();
log("");
sep();
log("  TEST 4: COMBINED ANALYSIS — OUTER-ONLY + PRESSURE CORRECTED");
sep();
log("");
log("  The strongest test: restrict to outer radii (minimal bar/non-circular");
log("  contamination) AND apply pressure support correction.");
log("");

const correctedForOuter = goldI45.map(g => {
  const corrGal = applyPressureCorrection(g);
  const outerPts = corrGal.points.filter(p => !p.isInner);
  return { ...corrGal, points: outerPts };
});

const margOuterCorr = [];
for (const gal of correctedForOuter) {
  const est = marginalizeSingleGalaxy(gal);
  if (est) margOuterCorr.push(est);
}
const outerCorrResult = hierarchicalDL(margOuterCorr, "Outer + pressure corrected");

log("  v4.0 baseline (full, uncorrected):");
printResult(baseline);
log("  Outer-only (no pressure corr):");
printResult(outerResult);
log("  Outer + pressure corrected:");
printResult(outerCorrResult);

if (baseline && outerCorrResult) {
  const totalShift = Math.abs(outerCorrResult.logA0 - baseline.logA0);
  const totalPct = Math.abs(outerCorrResult.a0 - baseline.a0) / baseline.a0 * 100;
  const tauDelta = outerCorrResult.tau - baseline.tau;
  log("");
  log("  Total shift from baseline: " + totalPct.toFixed(1) + "% (" + totalShift.toFixed(3) + " dex)");
  log("  Tau change: " + (tauDelta > 0 ? "+" : "") + tauDelta.toFixed(4));
  log("  I2 change: " + baseline.I2 + "% -> " + outerCorrResult.I2 + "%");
}

sep();
log("");
sep();
log("  TEST 5: RESIDUAL ANALYSIS");
sep();
log("");

function computeResiduals(galEsts, a0) {
  return galEsts.map(g => ({
    name: g.name,
    logA0: g.logA0,
    residual: g.logA0 - Math.log10(a0),
    Vmax: g.Vmax,
    distance: g.distance,
    inc: g.inc,
    fgas: g.fgas,
    Q_kin: g.Q_kin
  }));
}

function pearsonR(xs, ys) {
  const n = xs.length;
  if (n < 3) return NaN;
  const mx = xs.reduce((a, b) => a + b, 0) / n;
  const my = ys.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) {
    const dx = xs[i] - mx, dy = ys[i] - my;
    num += dx * dy;
    dx2 += dx * dx;
    dy2 += dy * dy;
  }
  return num / Math.sqrt(dx2 * dy2);
}

const residsBefore = computeResiduals(margGoldI45, baseline.a0);
const residsAfter = computeResiduals(margCorrected, corrBaseline ? corrBaseline.a0 : baseline.a0);

log("  Correlations of per-galaxy a0 residuals with galaxy properties:");
log("");
log("  BEFORE pressure correction:");
const props = ['Vmax', 'distance', 'inc', 'fgas', 'Q_kin'];
for (const prop of props) {
  const vals = residsBefore.map(r => r[prop]).filter(v => v !== undefined);
  const res = residsBefore.filter(r => r[prop] !== undefined).map(r => r.residual);
  const r = pearsonR(vals, res);
  log("    r(residual, " + prop + ") = " + r.toFixed(3) +
    (Math.abs(r) > 0.3 ? " *** SIGNIFICANT" : Math.abs(r) > 0.15 ? " * moderate" : " (weak)"));
}

log("");
log("  AFTER pressure correction:");
for (const prop of props) {
  const vals = residsAfter.map(r => r[prop]).filter(v => v !== undefined);
  const res = residsAfter.filter(r => r[prop] !== undefined).map(r => r.residual);
  const r = pearsonR(vals, res);
  log("    r(residual, " + prop + ") = " + r.toFixed(3) +
    (Math.abs(r) > 0.3 ? " *** SIGNIFICANT" : Math.abs(r) > 0.15 ? " * moderate" : " (weak)"));
}

sep();
log("");
sep();
log("  PHASE 5 VERDICT");
sep();
log("");

const results = {
  version: VERSION,
  timestamp: TIMESTAMP,
  description: "Phase 5: Kinematic Contamination Audit — non-circular motions + pressure support",
  baseline: baseline,
  test1_innerOuter: {
    inner: innerResult,
    outer: outerResult,
    delta_dex: innerResult && outerResult ? Math.abs(innerResult.logA0 - outerResult.logA0) : null
  },
  test1_Qkin: {
    highQ: hqResult,
    lowQ: lqResult
  },
  test1_etaRot: {
    highEta: heResult,
    lowEta: leResult
  },
  test1_Sout: {
    highSout: hsResult,
    lowSout: lsResult
  },
  test2_pressureCorrection: {
    corrected: corrBaseline,
    shift_dex: baseline && corrBaseline ? corrBaseline.logA0 - baseline.logA0 : null,
    tau_change: baseline && corrBaseline ? corrBaseline.tau - baseline.tau : null,
    byMass: { lowVmax: clv, midVmax: cmv, highVmax: chv },
    byFgas: { highFgas: cfh, lowFgas: cfl }
  },
  test3_redFlags: {
    distance: { near: cdl, far: cdh, split_pct: cdl && cdh ? Math.abs(cdl.a0 - cdh.a0) / baseline.a0 * 100 : null },
    inclination: { high: cih, low: cil, split_pct: cih && cil ? Math.abs(cih.a0 - cil.a0) / baseline.a0 * 100 : null },
    velocity: { high: chv, mid: cmv, split_pct: chv && cmv ? Math.abs(chv.a0 - cmv.a0) / baseline.a0 * 100 : null }
  },
  test4_outerCorrected: outerCorrResult,
  test5_residualCorrelations: {
    before: {},
    after: {}
  }
};

for (const prop of props) {
  const valsBefore = residsBefore.map(r => r[prop]).filter(v => v !== undefined);
  const resBefore = residsBefore.filter(r => r[prop] !== undefined).map(r => r.residual);
  results.test5_residualCorrelations.before[prop] = +pearsonR(valsBefore, resBefore).toFixed(3);

  const valsAfter = residsAfter.map(r => r[prop]).filter(v => v !== undefined);
  const resAfter = residsAfter.filter(r => r[prop] !== undefined).map(r => r.residual);
  results.test5_residualCorrelations.after[prop] = +pearsonR(valsAfter, resAfter).toFixed(3);
}

let pass = true;
let verdictLines = [];

if (baseline && corrBaseline) {
  const shift = Math.abs(corrBaseline.logA0 - baseline.logA0);
  if (shift > 0.15) {
    pass = false;
    verdictLines.push("FAIL: Pressure correction shifts a0 by " + shift.toFixed(3) + " dex (>0.15)");
  } else {
    verdictLines.push("PASS: Pressure correction shifts a0 by only " + shift.toFixed(3) + " dex (<0.15)");
  }
}

if (innerResult && outerResult) {
  const ioSplit = Math.abs(innerResult.logA0 - outerResult.logA0);
  if (ioSplit > 0.2) {
    pass = false;
    verdictLines.push("FAIL: Inner-outer split = " + ioSplit.toFixed(3) + " dex (>0.2) — non-circular contamination");
  } else {
    verdictLines.push("PASS: Inner-outer split = " + ioSplit.toFixed(3) + " dex (<0.2) — non-circular motions controlled");
  }
}

const redFlagImproved = [];
if (results.test3_redFlags.distance.split_pct !== null) {
  const before = 35;
  const after = results.test3_redFlags.distance.split_pct;
  redFlagImproved.push("Distance: " + before + "% -> " + after.toFixed(1) + "%");
  if (after > 30) verdictLines.push("FLAG: Distance split still large (" + after.toFixed(1) + "%)");
}
if (results.test3_redFlags.inclination.split_pct !== null) {
  const before = 42;
  const after = results.test3_redFlags.inclination.split_pct;
  redFlagImproved.push("Inclination: " + before + "% -> " + after.toFixed(1) + "%");
  if (after > 30) verdictLines.push("FLAG: Inclination split still large (" + after.toFixed(1) + "%)");
}
if (results.test3_redFlags.velocity.split_pct !== null) {
  const before = 36;
  const after = results.test3_redFlags.velocity.split_pct;
  redFlagImproved.push("Velocity: " + before + "% -> " + after.toFixed(1) + "%");
  if (after > 30) verdictLines.push("FLAG: Velocity split still large (" + after.toFixed(1) + "%)");
}

results.verdict = {
  pass: pass,
  summary: pass ? "a0 survives kinematic contamination tests" : "a0 shows sensitivity to kinematic corrections",
  details: verdictLines,
  redFlagChanges: redFlagImproved
};

for (const line of verdictLines) {
  log("  " + line);
}
log("");
if (redFlagImproved.length > 0) {
  log("  Red flag changes after pressure correction:");
  for (const rf of redFlagImproved) log("    " + rf);
}
log("");
log("  OVERALL: " + results.verdict.summary);
log("");

const outPath = path.join(__dirname, '..', 'public', 'phase5-kinematic-results.json');
fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
log("  Results saved to: public/phase5-kinematic-results.json");
log("");
log("================================================================================");
log("  PHASE 5 COMPLETE");
log("================================================================================");
