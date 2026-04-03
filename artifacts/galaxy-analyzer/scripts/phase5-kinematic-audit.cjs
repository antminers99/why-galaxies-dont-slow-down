const fs = require('fs');
const path = require('path');

const VERSION = "5.1.0";
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

let sparcTable = null;
const sparcTablePath = path.join(__dirname, '..', 'public', 'sparc-table.json');
if (fs.existsSync(sparcTablePath)) {
  sparcTable = JSON.parse(fs.readFileSync(sparcTablePath, 'utf8'));
}

let hasRotmod = false;
const rotmodDir = '/tmp/rotmod';
if (fs.existsSync(rotmodDir)) {
  const rotFiles = fs.readdirSync(rotmodDir).filter(f => f.endsWith('_rotmod.dat'));
  hasRotmod = rotFiles.length >= 170;
}

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
    name: gal.name,
    logA0: meanLogA,
    a0: Math.pow(10, meanLogA),
    se: seLogA,
    n: gal.points.length,
    Vmax: gal.Vmax,
    distance: gal.distance,
    inc: gal.inc,
    fgas: gal.fgas,
    gRange: gRange,
    Rdisk: gal.Rdisk || 1.0,
    Q_kin: gal.Q_kin || 0.7,
    Q_sparc: gal.Q_sparc,
    T: gal.T,
    eD: gal.eD,
    eInc: gal.eInc,
    fD: gal.fD,
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

log("================================================================================");
log("  PHASE 5 v5.1: KINEMATIC CONTAMINATION AUDIT (CORRECTED)");
log("================================================================================");
log("  Version: " + VERSION);
log("  Date: " + TIMESTAMP);
log("  SPARC master table: " + (sparcTable ? "LOADED (" + sparcTable.length + " galaxies)" : "NOT FOUND"));
log("  SPARC rotmod files: " + (hasRotmod ? "AVAILABLE" : "NOT AVAILABLE"));
log("  Purpose: Test whether a0 survives non-circular motion and pressure");
log("           support corrections using CORRECT v4 sample definitions.");
log("");

const sparcLookup = {};
if (sparcTable) {
  for (const g of sparcTable) sparcLookup[g.name] = g;
}

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
    Q_kin: g.Q_kin, eta_rot: g.eta_rot, eta_bar: g.eta_bar,
    S_out: g.S_out, sigma_bar: g.sigma_bar,
    Q_sparc: sparc.Q || null,
    T: sparc.T !== undefined ? sparc.T : null,
    eD: sparc.eD || null,
    eInc: sparc.eInc || null,
    fD: sparc.fD || null,
    Vflat: sparc.Vflat || null,
    RHI: sparc.RHI || null,
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
        Q_sparc: null, T: null, eD: null, eInc: null, fD: null,
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

log("  Total galaxies (n>=5): " + allGalaxies.length);
log("  SPARC galaxies: " + sparcGals.length);

if (sparcTable) {
  log("  SPARC Q=1 (high): " + sparcGals.filter(g => g.Q_sparc === 1).length);
  log("  SPARC Q=2 (medium): " + sparcGals.filter(g => g.Q_sparc === 2).length);
  log("  SPARC Q=3 (low): " + sparcGals.filter(g => g.Q_sparc === 3).length);
  log("  Hubble types: " + JSON.stringify(sparcGals.reduce((a, g) => {
    if (g.T === null) return a;
    const m = { 0: 'S0', 1: 'Sa', 2: 'Sab', 3: 'Sb', 4: 'Sbc', 5: 'Sc', 6: 'Scd', 7: 'Sd', 8: 'Sdm', 9: 'Sm', 10: 'Im', 11: 'BCD' };
    const k = m[g.T] || 'T' + g.T;
    a[k] = (a[k] || 0) + 1; return a;
  }, {})));
}

log("");

sep();
log("  STEP 0: MARGINALIZE ALL GALAXIES (using v4 method)");
sep();
log("");

const margAll = [];
let nMarg = 0;
for (const gal of allGalaxies) {
  const est = marginalizeSingleGalaxy(gal);
  if (est) {
    margAll.push(est);
    nMarg++;
  }
}
log("  Marginalized: " + nMarg + " galaxies");
log("");

function selectAndHier(label, filterFn) {
  const gals = margAll.filter(filterFn);
  const result = hierarchicalDL(gals, label);
  return result;
}

const goldI45Filter = g => g.Vmax >= 50 && g.inc >= 45 && g.n >= 5 && g.gRange >= 1.0;
const goldI45 = margAll.filter(goldI45Filter);
const baseline = hierarchicalDL(goldI45, "GOLD+i45 (v4 def)");

log("  SAMPLE SIZES (v4 definitions):");
log("  ALL:        " + margAll.filter(g => true).length);
log("  CLEAN:      " + margAll.filter(g => g.Vmax >= 80 && g.inc >= 30 && g.n >= 5 && g.gRange >= 0.5).length);
log("  GOLD:       " + margAll.filter(g => g.Vmax >= 50 && g.inc >= 30 && g.n >= 5 && g.gRange >= 1.0).length);
log("  GOLD+i45:   " + goldI45.length + " galaxies  <-- THIS IS THE CORRECT SAMPLE");
log("  HIGH-MASS:  " + margAll.filter(g => g.Vmax > 150).length);
log("  LOW-MASS:   " + margAll.filter(g => g.Vmax < 80).length);
log("");

function printResult(r) {
  if (!r) { log("    (insufficient data)"); return; }
  log("    " + pad(r.label, 30) + " n=" + padr(r.nGal, 4) +
    " a0=" + padr(r.a0, 5) + " tau=" + padr(r.tau, 7) +
    " ratio=" + padr(r.ratio, 7) + " 68%=[" + r.lo1sig + "," + r.hi1sig + "]");
}

log("  BASELINE (matching v4.0):");
printResult(baseline);
log("");

sep();
log("  TEST 1: NON-CIRCULAR MOTION PROXIES");
sep();
log("");

log("  TEST 1A: Inner vs Outer radii");
log("  (Using v4 GOLD+i45 galaxies, split by isInner flag)");
log("");

function buildSubsetGalaxies(sourceGals, pointFilter) {
  const results = [];
  for (const est of sourceGals) {
    const gal = galaxyMap[est.name];
    if (!gal) continue;
    const filtPts = gal.points.filter(pointFilter);
    if (filtPts.length < 5) continue;
    const subGal = { ...gal, points: filtPts };
    const subEst = marginalizeSingleGalaxy(subGal);
    if (subEst && subEst.gRange >= 0.3) results.push(subEst);
  }
  return results;
}

const goldI45Names = new Set(goldI45.map(g => g.name));
const goldI45Gals = allGalaxies.filter(g => goldI45Names.has(g.name));

const innerEsts = buildSubsetGalaxies(goldI45, p => p.isInner);
const outerEsts = buildSubsetGalaxies(goldI45, p => !p.isInner);
const innerResult = hierarchicalDL(innerEsts, "Inner only (GOLD+i45)");
const outerResult = hierarchicalDL(outerEsts, "Outer only (GOLD+i45)");

printResult(baseline);
printResult(innerResult);
printResult(outerResult);

if (innerResult && outerResult) {
  const ioD = Math.abs(innerResult.logA0 - outerResult.logA0);
  log("");
  log("  Inner-outer delta: " + ioD.toFixed(3) + " dex (" +
    (Math.abs(innerResult.a0 - outerResult.a0) / baseline.a0 * 100).toFixed(1) + "%)");
  log("  " + (ioD < 0.1 ? "PASS: non-circular motions NOT dominant" :
    ioD < 0.2 ? "WARNING: moderate non-circular contamination" :
      "FAIL: non-circular motions are a major concern"));
  log("  Outer tau: " + outerResult.tau + " vs full tau: " + baseline.tau +
    (outerResult.tau < baseline.tau ? " (outer is cleaner)" : ""));
}

log("");
log("  TEST 1B: SPARC Quality flag Q split (from master table)");
log("");

if (sparcTable) {
  const q1 = goldI45.filter(g => g.Q_sparc === 1);
  const q2 = goldI45.filter(g => g.Q_sparc === 2);
  const q12 = goldI45.filter(g => g.Q_sparc <= 2);
  const q3 = goldI45.filter(g => g.Q_sparc === 3);
  const q1r = hierarchicalDL(q1, "Q=1 (high quality)");
  const q2r = hierarchicalDL(q2, "Q=2 (medium quality)");
  const q12r = hierarchicalDL(q12, "Q<=2 (high+medium)");
  const q3r = hierarchicalDL(q3, "Q=3 (low quality)");
  printResult(q1r);
  printResult(q2r);
  printResult(q12r);
  printResult(q3r);
  if (q1r && q2r) {
    log("  Q=1 vs Q=2 delta: " + Math.abs(q1r.logA0 - q2r.logA0).toFixed(3) + " dex");
  }
} else {
  log("  (SPARC master table not available)");
}

log("");
log("  TEST 1C: Hubble type split (early vs late spirals)");
log("");

if (sparcTable) {
  const early = goldI45.filter(g => g.T !== null && g.T <= 5);
  const late = goldI45.filter(g => g.T !== null && g.T > 5);
  const earlyR = hierarchicalDL(early, "Early type (T<=5, S0-Sc)");
  const lateR = hierarchicalDL(late, "Late type (T>5, Scd-BCD)");
  printResult(earlyR);
  printResult(lateR);
  if (earlyR && lateR) {
    const tDelta = Math.abs(earlyR.logA0 - lateR.logA0);
    log("  Morphology delta: " + tDelta.toFixed(3) + " dex (" +
      (Math.abs(earlyR.a0 - lateR.a0) / baseline.a0 * 100).toFixed(1) + "%)");
  }
}

log("");
log("  TEST 1D: Kinematic quality Q_kin split");
log("");

const highQk = goldI45.filter(g => g.Q_kin >= 0.9);
const lowQk = goldI45.filter(g => g.Q_kin < 0.9);
const hqR = hierarchicalDL(highQk, "Q_kin >= 0.9");
const lqR = hierarchicalDL(lowQk, "Q_kin < 0.9");
printResult(hqR);
printResult(lqR);

log("");
log("  TEST 1E: Distance method quality");
log("");

if (sparcTable) {
  const precise = goldI45.filter(g => g.fD === 2 || g.fD === 3 || g.fD === 5);
  const hubbleFlow = goldI45.filter(g => g.fD === 1);
  const precR = hierarchicalDL(precise, "TRGB/Cepheid/SN dist.");
  const hfR = hierarchicalDL(hubbleFlow, "Hubble flow dist.");
  printResult(precR);
  printResult(hfR);
  if (precR && hfR) {
    log("  Distance method delta: " + Math.abs(precR.logA0 - hfR.logA0).toFixed(3) + " dex");
  }
}

sep();
log("");
sep();
log("  TEST 2: PRESSURE SUPPORT / ASYMMETRIC DRIFT CORRECTION");
sep();
log("");
log("  sigma_HI = " + SIGMA_HI + " km/s, h_gas = 2*Rdisk");
log("  Correction: g_obs_corr = g_obs + sigma^2/h_gas (in (km/s)^2/kpc)");
log("");

function applyPressureCorrection(galName, points, Rdisk) {
  const h_gas = 2.0 * (Rdisk || 1.0);
  const corrGal = SIGMA_HI * SIGMA_HI / h_gas;
  return points.map(p => {
    const gobs = Math.pow(10, p.logGobs);
    const gobs_corr = gobs + corrGal;
    return { ...p, logGobs: Math.log10(gobs_corr), correction_dex: Math.log10(gobs_corr) - p.logGobs };
  });
}

const goldI45Corrected = [];
for (const est of goldI45) {
  const gal = galaxyMap[est.name];
  if (!gal) continue;
  const corrPts = applyPressureCorrection(est.name, gal.points, gal.Rdisk);
  const corrGal = { ...gal, points: corrPts };
  const corrEst = marginalizeSingleGalaxy(corrGal);
  if (corrEst && corrEst.gRange >= 1.0) goldI45Corrected.push(corrEst);
}

const corrBaseline = hierarchicalDL(goldI45Corrected, "GOLD+i45 pressure-corrected");

log("  Correction magnitudes:");
const allCorr = [];
for (const est of goldI45) {
  const gal = galaxyMap[est.name];
  if (!gal) continue;
  const corrPts = applyPressureCorrection(est.name, gal.points, gal.Rdisk);
  corrPts.forEach(p => allCorr.push(p.correction_dex));
}
log("    Mean: " + (allCorr.reduce((a, b) => a + b, 0) / allCorr.length).toFixed(4) + " dex");
log("    Median: " + median(allCorr).toFixed(4) + " dex");
log("    Max: " + Math.max(...allCorr).toFixed(4) + " dex");

const lowVmaxCorr = [];
const hiVmaxCorr = [];
for (const est of goldI45) {
  const gal = galaxyMap[est.name];
  if (!gal) continue;
  const corrPts = applyPressureCorrection(est.name, gal.points, gal.Rdisk);
  if (gal.Vmax < 100) corrPts.forEach(p => lowVmaxCorr.push(p.correction_dex));
  else corrPts.forEach(p => hiVmaxCorr.push(p.correction_dex));
}
if (lowVmaxCorr.length) log("    Mean (Vmax<100): " + (lowVmaxCorr.reduce((a, b) => a + b, 0) / lowVmaxCorr.length).toFixed(4) + " dex");
if (hiVmaxCorr.length) log("    Mean (Vmax>=100): " + (hiVmaxCorr.reduce((a, b) => a + b, 0) / hiVmaxCorr.length).toFixed(4) + " dex");

log("");
log("  TEST 2A: Before vs After (GOLD+i45)");
log("");
printResult(baseline);
printResult(corrBaseline);

if (baseline && corrBaseline) {
  const shift = corrBaseline.logA0 - baseline.logA0;
  const shiftPct = (corrBaseline.a0 / baseline.a0 - 1) * 100;
  log("");
  log("  Shift: " + (shiftPct > 0 ? "+" : "") + shiftPct.toFixed(1) + "% (" +
    (shift > 0 ? "+" : "") + shift.toFixed(4) + " dex)");
  log("  Tau: " + baseline.tau + " -> " + corrBaseline.tau +
    " (delta=" + (corrBaseline.tau - baseline.tau).toFixed(4) + ")");
}

log("");
log("  TEST 2B: By mass regime (corrected)");
log("");

const corrLV = goldI45Corrected.filter(g => g.Vmax < 100);
const corrMV = goldI45Corrected.filter(g => g.Vmax >= 100 && g.Vmax <= 150);
const corrHV = goldI45Corrected.filter(g => g.Vmax > 150);
printResult(hierarchicalDL(corrLV, "Vmax<100 corrected"));
printResult(hierarchicalDL(corrMV, "Vmax 100-150 corrected"));
printResult(hierarchicalDL(corrHV, "Vmax>150 corrected"));

log("");
log("  TEST 2C: By gas fraction (corrected)");
log("");
printResult(hierarchicalDL(goldI45Corrected.filter(g => g.fgas >= 0.3), "fgas>=0.3 corrected"));
printResult(hierarchicalDL(goldI45Corrected.filter(g => g.fgas < 0.3), "fgas<0.3 corrected"));

sep();
log("");
sep();
log("  TEST 3: RED FLAG RE-EXAMINATION");
sep();
log("");

function printRedFlag(name, v4pct, beforeFilter, afterFilter, ests, corrEsts) {
  const b1 = hierarchicalDL(ests.filter(beforeFilter), name + " before");
  const b2 = hierarchicalDL(ests.filter(afterFilter), name + " before (comp)");
  const a1 = hierarchicalDL(corrEsts.filter(beforeFilter), name + " corrected");
  const a2 = hierarchicalDL(corrEsts.filter(afterFilter), name + " corrected (comp)");

  log("  " + name.toUpperCase() + " (v4.0: " + v4pct + "% difference):");
  if (b1 && b2) {
    const bp = Math.abs(b1.a0 - b2.a0) / baseline.a0 * 100;
    log("    Before corr: " + b1.a0 + " vs " + b2.a0 + " (" + bp.toFixed(1) + "%)");
  }
  if (a1 && a2) {
    const ap = Math.abs(a1.a0 - a2.a0) / baseline.a0 * 100;
    log("    After corr:  " + a1.a0 + " vs " + a2.a0 + " (" + ap.toFixed(1) + "%)");
    log("    " + (ap < v4pct * 0.5 ? "IMPROVED SIGNIFICANTLY" :
      ap < v4pct * 0.8 ? "IMPROVED" : "NO SIGNIFICANT IMPROVEMENT"));
  }
  log("");
  return { before: b1 && b2 ? Math.abs(b1.a0 - b2.a0) / baseline.a0 * 100 : null,
           after: a1 && a2 ? Math.abs(a1.a0 - a2.a0) / baseline.a0 * 100 : null };
}

const rfDist = printRedFlag("Distance", 35,
  g => g.distance < 15, g => g.distance >= 15, goldI45, goldI45Corrected);
const rfInc = printRedFlag("Inclination", 42,
  g => g.inc >= 60, g => g.inc >= 45 && g.inc < 60, goldI45, goldI45Corrected);
const rfVel = printRedFlag("Velocity", 36,
  g => g.Vmax > 150, g => g.Vmax >= 80 && g.Vmax <= 150, goldI45, goldI45Corrected);

sep();
log("");
sep();
log("  TEST 4: COMBINED — OUTER + PRESSURE CORRECTED");
sep();
log("");

const outerCorrEsts = [];
for (const est of goldI45) {
  const gal = galaxyMap[est.name];
  if (!gal) continue;
  const outerPts = gal.points.filter(p => !p.isInner);
  if (outerPts.length < 5) continue;
  const corrPts = applyPressureCorrection(est.name, outerPts, gal.Rdisk);
  const corrGal = { ...gal, points: corrPts };
  const corrEst = marginalizeSingleGalaxy(corrGal);
  if (corrEst && corrEst.gRange >= 0.5) outerCorrEsts.push(corrEst);
}
const outerCorrResult = hierarchicalDL(outerCorrEsts, "Outer + pressure corr (GOLD)");

log("  v4.0 baseline (GOLD+i45):");
printResult(baseline);
log("  Outer-only:");
printResult(outerResult);
log("  Pressure-corrected:");
printResult(corrBaseline);
log("  Outer + pressure corrected:");
printResult(outerCorrResult);

if (baseline && outerCorrResult) {
  const totalShift = Math.abs(outerCorrResult.logA0 - baseline.logA0);
  log("");
  log("  Total shift from baseline: " + (Math.abs(outerCorrResult.a0 - baseline.a0) / baseline.a0 * 100).toFixed(1) + "% (" + totalShift.toFixed(3) + " dex)");
  log("  Tau: " + baseline.tau + " -> " + outerCorrResult.tau);
}

sep();
log("");
sep();
log("  TEST 5: RESIDUAL CORRELATIONS");
sep();
log("");

function pearsonR(xs, ys) {
  const n = xs.length;
  if (n < 3) return NaN;
  const mx = xs.reduce((a, b) => a + b, 0) / n;
  const my = ys.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) {
    const dx = xs[i] - mx, dy = ys[i] - my;
    num += dx * dy; dx2 += dx * dx; dy2 += dy * dy;
  }
  return (dx2 > 0 && dy2 > 0) ? num / Math.sqrt(dx2 * dy2) : 0;
}

function computeResidCorrs(ests, refA0, label) {
  const props = ['Vmax', 'distance', 'inc', 'fgas', 'Q_kin'];
  log("  " + label + ":");
  for (const prop of props) {
    const valid = ests.filter(g => g[prop] !== undefined && g[prop] !== null);
    const r = pearsonR(valid.map(g => g[prop]), valid.map(g => g.logA0 - Math.log10(refA0)));
    log("    r(residual, " + pad(prop, 10) + ") = " + padr(r.toFixed(3), 7) +
      (Math.abs(r) > 0.3 ? " *** SIGNIFICANT" : Math.abs(r) > 0.15 ? " * moderate" : ""));
  }
  if (sparcTable) {
    const withT = ests.filter(g => g.T !== null && g.T !== undefined);
    if (withT.length > 5) {
      const r = pearsonR(withT.map(g => g.T), withT.map(g => g.logA0 - Math.log10(refA0)));
      log("    r(residual, " + pad("HubbleT", 10) + ") = " + padr(r.toFixed(3), 7) +
        (Math.abs(r) > 0.3 ? " *** SIGNIFICANT" : Math.abs(r) > 0.15 ? " * moderate" : ""));
    }
  }
}

computeResidCorrs(goldI45, baseline.a0, "BEFORE pressure correction");
log("");
computeResidCorrs(goldI45Corrected, corrBaseline ? corrBaseline.a0 : baseline.a0, "AFTER pressure correction");

sep();
log("");
sep();
log("  PHASE 5 VERDICT");
sep();
log("");

const verdictLines = [];
let passCount = 0, failCount = 0;

if (innerResult && outerResult) {
  const ioD = Math.abs(innerResult.logA0 - outerResult.logA0);
  if (ioD < 0.15) {
    verdictLines.push("PASS: Inner-outer split = " + ioD.toFixed(3) + " dex — non-circular motions controlled");
    passCount++;
  } else {
    verdictLines.push("CONCERN: Inner-outer split = " + ioD.toFixed(3) + " dex");
    failCount++;
  }
}

if (baseline && corrBaseline) {
  const shift = Math.abs(corrBaseline.logA0 - baseline.logA0);
  if (shift < 0.15) {
    verdictLines.push("PASS: Pressure correction shift = " + shift.toFixed(3) + " dex (<0.15)");
    passCount++;
  } else {
    verdictLines.push("CONCERN: Pressure correction shift = " + shift.toFixed(3) + " dex");
    failCount++;
  }
  const tauD = corrBaseline.tau - baseline.tau;
  verdictLines.push("INFO: Tau change with pressure corr: " + (tauD > 0 ? "+" : "") + tauD.toFixed(4));
}

if (outerCorrResult) {
  const shift = Math.abs(outerCorrResult.logA0 - baseline.logA0);
  verdictLines.push("INFO: Combined (outer+pressure) shift = " + shift.toFixed(3) + " dex, tau = " + outerCorrResult.tau);
}

const redFlags = [];
if (rfDist.before !== null && rfDist.after !== null)
  redFlags.push("Distance: " + rfDist.before.toFixed(1) + "% -> " + rfDist.after.toFixed(1) + "%");
if (rfInc.before !== null && rfInc.after !== null)
  redFlags.push("Inclination: " + rfInc.before.toFixed(1) + "% -> " + rfInc.after.toFixed(1) + "%");
if (rfVel.before !== null && rfVel.after !== null)
  redFlags.push("Velocity: " + rfVel.before.toFixed(1) + "% -> " + rfVel.after.toFixed(1) + "%");

for (const v of verdictLines) log("  " + v);
log("");
if (redFlags.length > 0) {
  log("  Red flag changes:");
  for (const rf of redFlags) log("    " + rf);
}

const overallPass = passCount >= 2 && failCount === 0;
log("");
log("  OVERALL: " + (overallPass ? "PASS — a0 survives kinematic contamination tests" :
  "MIXED — some concerns remain"));

const output = {
  version: VERSION,
  timestamp: TIMESTAMP,
  description: "Phase 5 v5.1: Kinematic Contamination Audit (corrected sample definitions)",
  dataSource: {
    sparcTable: !!sparcTable,
    rotmodFiles: hasRotmod,
    sparcGalaxies: sparcGals.length,
    goldI45Count: goldI45.length
  },
  baseline: baseline,
  test1_innerOuter: {
    inner: innerResult,
    outer: outerResult,
    delta_dex: innerResult && outerResult ? +Math.abs(innerResult.logA0 - outerResult.logA0).toFixed(4) : null
  },
  test1_sparcQ: sparcTable ? {
    q1: hierarchicalDL(goldI45.filter(g => g.Q_sparc === 1), "Q=1"),
    q2: hierarchicalDL(goldI45.filter(g => g.Q_sparc === 2), "Q=2"),
    q12: hierarchicalDL(goldI45.filter(g => g.Q_sparc <= 2), "Q<=2"),
    q3: hierarchicalDL(goldI45.filter(g => g.Q_sparc === 3), "Q=3")
  } : null,
  test1_hubbleType: sparcTable ? {
    early: hierarchicalDL(goldI45.filter(g => g.T !== null && g.T <= 5), "Early (T<=5)"),
    late: hierarchicalDL(goldI45.filter(g => g.T !== null && g.T > 5), "Late (T>5)")
  } : null,
  test1_Qkin: { highQ: hqR, lowQ: lqR },
  test2_pressureCorrection: {
    corrected: corrBaseline,
    shift_dex: baseline && corrBaseline ? +(corrBaseline.logA0 - baseline.logA0).toFixed(4) : null,
    tau_change: baseline && corrBaseline ? +(corrBaseline.tau - baseline.tau).toFixed(4) : null
  },
  test3_redFlags: {
    distance: rfDist,
    inclination: rfInc,
    velocity: rfVel
  },
  test4_outerCorrected: outerCorrResult,
  verdict: {
    pass: overallPass,
    passCount, failCount,
    summary: overallPass ? "a0 survives kinematic contamination tests" : "Mixed results — some concerns remain",
    details: verdictLines,
    redFlagChanges: redFlags
  }
};

const outPath = path.join(__dirname, '..', 'public', 'phase5-kinematic-results.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
log("");
log("  Results saved to: public/phase5-kinematic-results.json");
log("");
log("================================================================================");
log("  PHASE 5 COMPLETE");
log("================================================================================");
