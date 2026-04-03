const fs = require('fs');
const path = require('path');

const VERSION = "6.0.0";
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
    gRange, Rdisk: gal.Rdisk || 1.0, Q_kin: gal.Q_kin || 0.7,
    Q_sparc: gal.Q_sparc, T: gal.T, eD: gal.eD, eInc: gal.eInc,
    fD: gal.fD, isLT: !!gal.isLT
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

function printResult(r, indent) {
  indent = indent || "    ";
  if (!r) { log(indent + "(insufficient data — fewer than 3 galaxies)"); return; }
  log(indent + pad(r.label, 35) + " n=" + padr(r.nGal, 3) +
    " a0=" + padr(r.a0, 5) + " tau=" + padr(r.tau, 6) +
    " ratio=" + padr(r.ratio, 6) + " 68%=[" + r.lo1sig + "," + r.hi1sig + "]");
}

function sampleStats(ests, label) {
  if (!ests.length) return null;
  return {
    label,
    n: ests.length,
    medianVmax: +median(ests.map(g => g.Vmax)).toFixed(1),
    medianInc: +median(ests.map(g => g.inc)).toFixed(1),
    medianDist: +median(ests.map(g => g.distance)).toFixed(1),
    medianFgas: +median(ests.map(g => g.fgas)).toFixed(3),
    medianN: +median(ests.map(g => g.n)).toFixed(0),
    medianGrange: +median(ests.map(g => g.gRange)).toFixed(2),
    meanVmax: +(ests.reduce((a, g) => a + g.Vmax, 0) / ests.length).toFixed(1),
    meanInc: +(ests.reduce((a, g) => a + g.inc, 0) / ests.length).toFixed(1),
    meanDist: +(ests.reduce((a, g) => a + g.distance, 0) / ests.length).toFixed(1),
  };
}

function printSampleComparison(s1, s2) {
  if (!s1 || !s2) return;
  log("      " + pad("", 18) + pad(s1.label, 16) + pad(s2.label, 16) + "Match?");
  log("      " + pad("n_gal", 18) + pad(s1.n, 16) + pad(s2.n, 16));
  log("      " + pad("median Vmax", 18) + pad(s1.medianVmax, 16) + pad(s2.medianVmax, 16) +
    (Math.abs(s1.medianVmax - s2.medianVmax) / Math.max(s1.medianVmax, s2.medianVmax) < 0.25 ? " OK" : " MISMATCH"));
  log("      " + pad("median inc", 18) + pad(s1.medianInc, 16) + pad(s2.medianInc, 16) +
    (Math.abs(s1.medianInc - s2.medianInc) < 10 ? " OK" : " MISMATCH"));
  log("      " + pad("median dist", 18) + pad(s1.medianDist, 16) + pad(s2.medianDist, 16) +
    (Math.abs(s1.medianDist - s2.medianDist) / Math.max(s1.medianDist, s2.medianDist) < 0.5 ? " OK" : " MISMATCH"));
  log("      " + pad("median fgas", 18) + pad(s1.medianFgas, 16) + pad(s2.medianFgas, 16) +
    (Math.abs(s1.medianFgas - s2.medianFgas) < 0.15 ? " OK" : " MISMATCH"));
  log("      " + pad("median n_pts", 18) + pad(s1.medianN, 16) + pad(s2.medianN, 16));
  log("      " + pad("median gRange", 18) + pad(s1.medianGrange, 16) + pad(s2.medianGrange, 16));
}

function nearestNeighborMatch(poolA, poolB, matchProps, weights) {
  const matched = [];
  const usedB = new Set();

  for (const a of poolA) {
    let bestIdx = -1, bestDist = Infinity;
    for (let j = 0; j < poolB.length; j++) {
      if (usedB.has(j)) continue;
      let d2 = 0;
      for (let k = 0; k < matchProps.length; k++) {
        const prop = matchProps[k];
        const va = a[prop], vb = poolB[j][prop];
        if (va === null || va === undefined || vb === null || vb === undefined) continue;
        const diff = va - vb;
        d2 += (diff * diff) * (weights[k] || 1);
      }
      if (d2 < bestDist) { bestDist = d2; bestIdx = j; }
    }
    if (bestIdx >= 0) {
      matched.push({ a: a, b: poolB[bestIdx], dist: bestDist });
      usedB.add(bestIdx);
    }
  }
  return matched;
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
      logGbar: p.log_g_bar, logGobs: p.log_g_obs, isInner: p.isInner
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
    galaxyMap[p.g].points.push({
      logGbar: p.x, logGobs: p.x + p.y, isInner: false
    });
  }
}

const allGalaxies = Object.values(galaxyMap).filter(g => g.points.length >= 5);

log("================================================================================");
log("  PHASE 6: MATCHED-SAMPLE SPLIT RESOLUTION");
log("================================================================================");
log("  Version: " + VERSION);
log("  Date: " + TIMESTAMP);
log("  Purpose: Resolve three remaining red-flag splits from v4.0/Phase 5");
log("           using nearest-neighbor matched-sample analysis.");
log("  Method: For each split, create matched pairs controlling for");
log("          confounding variables, then compare a0 within matched sets.");
log("  Threshold: delta < 0.05 dex (~12%) = RESOLVED");
log("             delta 0.05-0.10 dex = INCONCLUSIVE");
log("             delta > 0.10 dex = UNRESOLVED (real effect)");
log("");

log("  Marginalizing all galaxies...");
const margAll = [];
for (const gal of allGalaxies) {
  const est = marginalizeSingleGalaxy(gal);
  if (est) margAll.push(est);
}
log("  Marginalized: " + margAll.length + " galaxies");

const goldI45 = margAll.filter(g =>
  g.Vmax >= 50 && g.inc >= 45 && g.n >= 5 && g.gRange >= 1.0
);
const baseline = hierarchicalDL(goldI45, "GOLD+i45 baseline");

log("  GOLD+i45: " + goldI45.length + " galaxies");
log("");
log("  BASELINE:");
printResult(baseline);
log("");

sep();
log("  SPLIT 1: DISTANCE METHOD");
sep();
log("");
log("  Question: Do galaxies with precise distances (TRGB/Cepheid/SN)");
log("  give different a0 than Hubble-flow galaxies?");
log("  Phase 5 v5.1 found: 0.105 dex split (27%).");
log("");

const precise = goldI45.filter(g => g.fD === 2 || g.fD === 3 || g.fD === 5);
const hubbleFlow = goldI45.filter(g => g.fD === 1);
const umaCluster = goldI45.filter(g => g.fD === 4);

log("  1.1 RAW SPLIT (unmatched):");
log("");
const precR = hierarchicalDL(precise, "Precise (TRGB/Ceph/SN)");
const hfR = hierarchicalDL(hubbleFlow, "Hubble flow");
const umaR = hierarchicalDL(umaCluster, "UMa cluster");
printResult(precR);
printResult(hfR);
printResult(umaR);

if (precR && hfR) {
  const rawDelta = Math.abs(precR.logA0 - hfR.logA0);
  log("");
  log("    Raw delta: " + rawDelta.toFixed(3) + " dex (" +
    (Math.abs(precR.a0 - hfR.a0) / baseline.a0 * 100).toFixed(1) + "%)");
}

log("");
log("  1.2 SAMPLE COMPOSITION:");
log("");
const precStats = sampleStats(precise, "Precise");
const hfStats = sampleStats(hubbleFlow, "HubbleFlow");
printSampleComparison(precStats, hfStats);

log("");
log("  1.3 MATCHED ANALYSIS:");
log("  Matching precise to Hubble-flow on: Vmax, inc, fgas, n_pts");
log("");

const wDist = [1.0 / (50 * 50), 1.0 / (10 * 10), 1.0 / (0.15 * 0.15), 1.0 / (10 * 10)];
const distMatches = nearestNeighborMatch(precise, hubbleFlow,
  ['Vmax', 'inc', 'fgas', 'n'], wDist);

log("    Matched pairs: " + distMatches.length);

if (distMatches.length >= 3) {
  const matchedPrec = distMatches.map(m => m.a);
  const matchedHF = distMatches.map(m => m.b);

  const mPrecStats = sampleStats(matchedPrec, "Matched-Prec");
  const mHFStats = sampleStats(matchedHF, "Matched-HF");
  printSampleComparison(mPrecStats, mHFStats);

  log("");
  const mPrecR = hierarchicalDL(matchedPrec, "Matched precise");
  const mHFR = hierarchicalDL(matchedHF, "Matched Hubble-flow");
  printResult(mPrecR);
  printResult(mHFR);

  if (mPrecR && mHFR) {
    const matchedDelta = Math.abs(mPrecR.logA0 - mHFR.logA0);
    const matchedPct = Math.abs(mPrecR.a0 - mHFR.a0) / baseline.a0 * 100;
    log("");
    log("    Matched delta: " + matchedDelta.toFixed(3) + " dex (" + matchedPct.toFixed(1) + "%)");

    const pairDeltas = distMatches.map(m => m.a.logA0 - m.b.logA0);
    const meanPairDelta = pairDeltas.reduce((a, b) => a + b, 0) / pairDeltas.length;
    const sePairDelta = Math.sqrt(pairDeltas.reduce((a, d) => a + (d - meanPairDelta) ** 2, 0) / (pairDeltas.length * (pairDeltas.length - 1)));
    const tStat = Math.abs(meanPairDelta / sePairDelta);
    log("    Paired t-test: mean diff = " + meanPairDelta.toFixed(4) +
      " +/- " + sePairDelta.toFixed(4) + ", t = " + tStat.toFixed(2) +
      " (df=" + (pairDeltas.length - 1) + ")");
    log("    " + (tStat < 2.0 ? "NOT significant (t < 2.0)" : "SIGNIFICANT (t >= 2.0)"));

    log("");
    if (matchedDelta < 0.05) {
      log("    VERDICT: RESOLVED. Distance split disappears after matching.");
      log("    The raw split was driven by sample composition differences,");
      log("    not by distance measurement method per se.");
    } else if (matchedDelta < 0.10) {
      log("    VERDICT: INCONCLUSIVE. Split reduced but not eliminated.");
    } else {
      log("    VERDICT: UNRESOLVED. Split persists after matching.");
      log("    Distance measurement method may genuinely affect a0.");
    }
  }
}

log("");
log("  1.4 DISTANCE BINNING (continuous):");
log("");
const distBins = [
  { label: "D < 10 Mpc", filter: g => g.distance < 10 },
  { label: "10 <= D < 20", filter: g => g.distance >= 10 && g.distance < 20 },
  { label: "20 <= D < 40", filter: g => g.distance >= 20 && g.distance < 40 },
  { label: "D >= 40 Mpc", filter: g => g.distance >= 40 }
];
for (const bin of distBins) {
  const sub = goldI45.filter(bin.filter);
  const r = hierarchicalDL(sub, bin.label);
  printResult(r);
}

sep();
log("");
sep();
log("  SPLIT 2: INCLINATION");
sep();
log("");
log("  Question: Do high-inclination galaxies (inc >= 60) give different");
log("  a0 than moderate-inclination (45-60)?");
log("  Phase 5 v5.1 found: 32.6% split.");
log("");

const highInc = goldI45.filter(g => g.inc >= 60);
const modInc = goldI45.filter(g => g.inc >= 45 && g.inc < 60);

log("  2.1 RAW SPLIT (unmatched):");
log("");
const hiR = hierarchicalDL(highInc, "Inc >= 60");
const miR = hierarchicalDL(modInc, "Inc 45-59");
printResult(hiR);
printResult(miR);

if (hiR && miR) {
  const rawDelta = Math.abs(hiR.logA0 - miR.logA0);
  log("");
  log("    Raw delta: " + rawDelta.toFixed(3) + " dex (" +
    (Math.abs(hiR.a0 - miR.a0) / baseline.a0 * 100).toFixed(1) + "%)");
}

log("");
log("  2.2 SAMPLE COMPOSITION:");
log("");
const hiStats = sampleStats(highInc, "Inc>=60");
const miStats = sampleStats(modInc, "Inc45-59");
printSampleComparison(hiStats, miStats);

log("");
log("  2.3 MATCHED ANALYSIS:");
log("  Matching high-inc to moderate-inc on: Vmax, distance, fgas, n_pts");
log("");

const wInc = [1.0 / (50 * 50), 1.0 / (15 * 15), 1.0 / (0.15 * 0.15), 1.0 / (10 * 10)];

const smaller = highInc.length <= modInc.length ? highInc : modInc;
const larger = highInc.length <= modInc.length ? modInc : highInc;
const incMatches = nearestNeighborMatch(smaller, larger,
  ['Vmax', 'distance', 'fgas', 'n'], wInc);

log("    Matched pairs: " + incMatches.length);

if (incMatches.length >= 3) {
  const matchedSmall = incMatches.map(m => m.a);
  const matchedLarge = incMatches.map(m => m.b);

  const isHighFirst = highInc.length <= modInc.length;
  const matchedHi = isHighFirst ? matchedSmall : matchedLarge;
  const matchedMod = isHighFirst ? matchedLarge : matchedSmall;

  const mHiStats = sampleStats(matchedHi, "MatchHi");
  const mMiStats = sampleStats(matchedMod, "MatchMod");
  printSampleComparison(mHiStats, mMiStats);

  log("");
  const mHiR = hierarchicalDL(matchedHi, "Matched inc>=60");
  const mMiR = hierarchicalDL(matchedMod, "Matched inc 45-59");
  printResult(mHiR);
  printResult(mMiR);

  if (mHiR && mMiR) {
    const matchedDelta = Math.abs(mHiR.logA0 - mMiR.logA0);
    const matchedPct = Math.abs(mHiR.a0 - mMiR.a0) / baseline.a0 * 100;
    log("");
    log("    Matched delta: " + matchedDelta.toFixed(3) + " dex (" + matchedPct.toFixed(1) + "%)");

    const pairDeltas = incMatches.map(m => {
      const hi = isHighFirst ? m.a : m.b;
      const mod = isHighFirst ? m.b : m.a;
      return hi.logA0 - mod.logA0;
    });
    const meanPairDelta = pairDeltas.reduce((a, b) => a + b, 0) / pairDeltas.length;
    const sePairDelta = Math.sqrt(pairDeltas.reduce((a, d) => a + (d - meanPairDelta) ** 2, 0) / (pairDeltas.length * (pairDeltas.length - 1)));
    const tStat = Math.abs(meanPairDelta / sePairDelta);
    log("    Paired t-test: mean diff = " + meanPairDelta.toFixed(4) +
      " +/- " + sePairDelta.toFixed(4) + ", t = " + tStat.toFixed(2) +
      " (df=" + (pairDeltas.length - 1) + ")");
    log("    " + (tStat < 2.0 ? "NOT significant (t < 2.0)" : "SIGNIFICANT (t >= 2.0)"));

    log("");
    if (matchedDelta < 0.05) {
      log("    VERDICT: RESOLVED. Inclination split disappears after matching.");
    } else if (matchedDelta < 0.10) {
      log("    VERDICT: INCONCLUSIVE. Split reduced but not eliminated.");
    } else {
      log("    VERDICT: UNRESOLVED. Inclination split persists.");
    }
  }
}

log("");
log("  2.4 INCLINATION BINS (continuous):");
log("");
const incBins = [
  { label: "45 <= inc < 55", filter: g => g.inc >= 45 && g.inc < 55 },
  { label: "55 <= inc < 65", filter: g => g.inc >= 55 && g.inc < 65 },
  { label: "65 <= inc < 75", filter: g => g.inc >= 65 && g.inc < 75 },
  { label: "inc >= 75", filter: g => g.inc >= 75 }
];
for (const bin of incBins) {
  const sub = goldI45.filter(bin.filter);
  const r = hierarchicalDL(sub, bin.label);
  printResult(r);
}

sep();
log("");
sep();
log("  SPLIT 3: VELOCITY (MASS)");
sep();
log("");
log("  Question: Do high-Vmax and mid-Vmax galaxies give different a0?");
log("  Phase 5 v5.1 found: 34.7% split (Vmax>150 vs 80-150).");
log("");

const hiV = goldI45.filter(g => g.Vmax > 150);
const midV = goldI45.filter(g => g.Vmax >= 80 && g.Vmax <= 150);
const loV = goldI45.filter(g => g.Vmax >= 50 && g.Vmax < 80);

log("  3.1 RAW SPLIT (unmatched):");
log("");
const hvR = hierarchicalDL(hiV, "Vmax > 150");
const mvR = hierarchicalDL(midV, "Vmax 80-150");
const lvR = hierarchicalDL(loV, "Vmax 50-80");
printResult(hvR);
printResult(mvR);
printResult(lvR);

if (hvR && mvR) {
  const rawDelta = Math.abs(hvR.logA0 - mvR.logA0);
  log("");
  log("    Raw delta (hi vs mid): " + rawDelta.toFixed(3) + " dex (" +
    (Math.abs(hvR.a0 - mvR.a0) / baseline.a0 * 100).toFixed(1) + "%)");
}

log("");
log("  3.2 SAMPLE COMPOSITION:");
log("");
const hvStats = sampleStats(hiV, "Vmax>150");
const mvStats = sampleStats(midV, "Vmax80-150");
printSampleComparison(hvStats, mvStats);

log("");
log("  3.3 MATCHED ANALYSIS:");
log("  Matching Vmax>150 to Vmax 80-150 on: distance, inc, fgas, n_pts");
log("");

const wVel = [1.0 / (15 * 15), 1.0 / (10 * 10), 1.0 / (0.15 * 0.15), 1.0 / (10 * 10)];

const velSmaller = hiV.length <= midV.length ? hiV : midV;
const velLarger = hiV.length <= midV.length ? midV : hiV;
const velMatches = nearestNeighborMatch(velSmaller, velLarger,
  ['distance', 'inc', 'fgas', 'n'], wVel);

log("    Matched pairs: " + velMatches.length);

if (velMatches.length >= 3) {
  const matchedVSmall = velMatches.map(m => m.a);
  const matchedVLarge = velMatches.map(m => m.b);

  const isHiFirst = hiV.length <= midV.length;
  const matchedHiV = isHiFirst ? matchedVSmall : matchedVLarge;
  const matchedMidV = isHiFirst ? matchedVLarge : matchedVSmall;

  const mHvStats = sampleStats(matchedHiV, "MatchHi");
  const mMvStats = sampleStats(matchedMidV, "MatchMid");
  printSampleComparison(mHvStats, mMvStats);

  log("");
  const mHvR = hierarchicalDL(matchedHiV, "Matched Vmax>150");
  const mMvR = hierarchicalDL(matchedMidV, "Matched Vmax 80-150");
  printResult(mHvR);
  printResult(mMvR);

  if (mHvR && mMvR) {
    const matchedDelta = Math.abs(mHvR.logA0 - mMvR.logA0);
    const matchedPct = Math.abs(mHvR.a0 - mMvR.a0) / baseline.a0 * 100;
    log("");
    log("    Matched delta: " + matchedDelta.toFixed(3) + " dex (" + matchedPct.toFixed(1) + "%)");

    const pairDeltas = velMatches.map(m => {
      const hi = isHiFirst ? m.a : m.b;
      const mid = isHiFirst ? m.b : m.a;
      return hi.logA0 - mid.logA0;
    });
    const meanPairDelta = pairDeltas.reduce((a, b) => a + b, 0) / pairDeltas.length;
    const sePairDelta = Math.sqrt(pairDeltas.reduce((a, d) => a + (d - meanPairDelta) ** 2, 0) / (pairDeltas.length * (pairDeltas.length - 1)));
    const tStat = Math.abs(meanPairDelta / sePairDelta);
    log("    Paired t-test: mean diff = " + meanPairDelta.toFixed(4) +
      " +/- " + sePairDelta.toFixed(4) + ", t = " + tStat.toFixed(2) +
      " (df=" + (pairDeltas.length - 1) + ")");
    log("    " + (tStat < 2.0 ? "NOT significant (t < 2.0)" : "SIGNIFICANT (t >= 2.0)"));

    log("");
    if (matchedDelta < 0.05) {
      log("    VERDICT: RESOLVED. Velocity split disappears after matching.");
    } else if (matchedDelta < 0.10) {
      log("    VERDICT: INCONCLUSIVE. Split reduced but not eliminated.");
    } else {
      log("    VERDICT: UNRESOLVED. Velocity split persists.");
    }
  }
}

log("");
log("  3.4 VELOCITY BINS (continuous):");
log("");
const velBins = [
  { label: "50 <= Vmax < 80", filter: g => g.Vmax >= 50 && g.Vmax < 80 },
  { label: "80 <= Vmax < 120", filter: g => g.Vmax >= 80 && g.Vmax < 120 },
  { label: "120 <= Vmax < 180", filter: g => g.Vmax >= 120 && g.Vmax < 180 },
  { label: "180 <= Vmax < 250", filter: g => g.Vmax >= 180 && g.Vmax < 250 },
  { label: "Vmax >= 250", filter: g => g.Vmax >= 250 }
];
for (const bin of velBins) {
  const sub = goldI45.filter(bin.filter);
  const r = hierarchicalDL(sub, bin.label);
  printResult(r);
}

sep();
log("");
sep();
log("  TRIPLE-MATCHED GRAND TEST");
sep();
log("");
log("  Simultaneously match all galaxies into pairs that differ in ONLY");
log("  one property at a time, while controlling all others.");
log("");

const allProps = goldI45.map(g => ({
  ...g,
  logVmax: Math.log10(g.Vmax),
  logDist: Math.log10(g.distance)
}));

function jackknifeDelta(ests, splitFn, label) {
  const groupA = ests.filter(g => splitFn(g));
  const groupB = ests.filter(g => !splitFn(g));
  if (groupA.length < 3 || groupB.length < 3) return null;

  const rA = hierarchicalDL(groupA, label + " A");
  const rB = hierarchicalDL(groupB, label + " B");
  if (!rA || !rB) return null;

  const jkDeltas = [];
  for (let i = 0; i < ests.length; i++) {
    const jkEsts = ests.filter((_, j) => j !== i);
    const jkA = jkEsts.filter(g => splitFn(g));
    const jkB = jkEsts.filter(g => !splitFn(g));
    if (jkA.length < 3 || jkB.length < 3) continue;
    const jkRA = hierarchicalDL(jkA, "jk");
    const jkRB = hierarchicalDL(jkB, "jk");
    if (jkRA && jkRB) jkDeltas.push(jkRA.logA0 - jkRB.logA0);
  }

  const n = jkDeltas.length;
  if (n < 5) return { rA, rB, delta: Math.abs(rA.logA0 - rB.logA0), jkSE: null };
  const meanJK = jkDeltas.reduce((a, b) => a + b, 0) / n;
  const jkVar = ((n - 1) / n) * jkDeltas.reduce((a, d) => a + (d - meanJK) ** 2, 0);
  const jkSE = Math.sqrt(jkVar);

  return {
    rA, rB,
    delta: Math.abs(rA.logA0 - rB.logA0),
    jkSE,
    tStat: jkSE > 0 ? Math.abs(rA.logA0 - rB.logA0) / jkSE : Infinity
  };
}

log("  Jackknife significance test for each split within GOLD+i45:");
log("");

const jkDist = jackknifeDelta(goldI45,
  g => g.fD === 2 || g.fD === 3 || g.fD === 5, "DistMethod");
const jkInc = jackknifeDelta(goldI45,
  g => g.inc >= 60, "Inclination");
const jkVel = jackknifeDelta(goldI45,
  g => g.Vmax > 150, "Velocity");

function printJK(name, jk) {
  if (!jk) { log("    " + name + ": insufficient data"); return; }
  log("    " + pad(name, 15) + " delta=" + jk.delta.toFixed(3) + " dex" +
    (jk.jkSE ? " +/- " + jk.jkSE.toFixed(3) : "") +
    (jk.tStat ? " t=" + jk.tStat.toFixed(2) : "") +
    (jk.tStat && jk.tStat < 2.0 ? " (not significant)" :
      jk.tStat && jk.tStat >= 2.0 ? " (SIGNIFICANT)" : ""));
}

printJK("Distance", jkDist);
printJK("Inclination", jkInc);
printJK("Velocity", jkVel);

sep();
log("");
sep();
log("  PHASE 6 SUMMARY");
sep();
log("");

const results = {
  version: VERSION,
  timestamp: TIMESTAMP,
  description: "Phase 6: Matched-sample split resolution",
  baseline: baseline,
  threshold: "< 0.05 dex = RESOLVED, 0.05-0.10 = INCONCLUSIVE, > 0.10 = UNRESOLVED",
  splits: {}
};

function summarizeSplit(name, rawR1, rawR2, matchedR1, matchedR2, jk) {
  const raw = rawR1 && rawR2 ? {
    delta_dex: +Math.abs(rawR1.logA0 - rawR2.logA0).toFixed(4),
    delta_pct: +(Math.abs(rawR1.a0 - rawR2.a0) / baseline.a0 * 100).toFixed(1),
    r1: rawR1, r2: rawR2
  } : null;

  const matched = matchedR1 && matchedR2 ? {
    delta_dex: +Math.abs(matchedR1.logA0 - matchedR2.logA0).toFixed(4),
    delta_pct: +(Math.abs(matchedR1.a0 - matchedR2.a0) / baseline.a0 * 100).toFixed(1),
    r1: matchedR1, r2: matchedR2
  } : null;

  const effectiveDelta = matched ? matched.delta_dex : (raw ? raw.delta_dex : null);
  let verdict = "UNKNOWN";
  if (effectiveDelta !== null) {
    if (effectiveDelta < 0.05) verdict = "RESOLVED";
    else if (effectiveDelta < 0.10) verdict = "INCONCLUSIVE";
    else verdict = "UNRESOLVED";
  }

  const jackknife = jk ? {
    delta: +jk.delta.toFixed(4),
    se: jk.jkSE ? +jk.jkSE.toFixed(4) : null,
    tStat: jk.tStat ? +jk.tStat.toFixed(2) : null,
    significant: jk.tStat ? jk.tStat >= 2.0 : null
  } : null;

  log("  " + pad(name, 15) +
    " raw=" + (raw ? raw.delta_dex.toFixed(3) : "N/A") + " dex" +
    " matched=" + (matched ? matched.delta_dex.toFixed(3) : "N/A") + " dex" +
    " jk_t=" + (jackknife && jackknife.tStat ? jackknife.tStat.toFixed(2) : "N/A") +
    " => " + verdict);

  results.splits[name] = { raw, matched, jackknife, verdict };
  return verdict;
}

const distMatchedR1 = distMatches.length >= 3 ?
  hierarchicalDL(distMatches.map(m => m.a), "Matched precise") : null;
const distMatchedR2 = distMatches.length >= 3 ?
  hierarchicalDL(distMatches.map(m => m.b), "Matched Hubble-flow") : null;

let incMatchedR1 = null, incMatchedR2 = null;
if (incMatches.length >= 3) {
  const isHighFirst = highInc.length <= modInc.length;
  const mHi = isHighFirst ? incMatches.map(m => m.a) : incMatches.map(m => m.b);
  const mMod = isHighFirst ? incMatches.map(m => m.b) : incMatches.map(m => m.a);
  incMatchedR1 = hierarchicalDL(mHi, "Matched inc>=60");
  incMatchedR2 = hierarchicalDL(mMod, "Matched inc 45-59");
}

let velMatchedR1 = null, velMatchedR2 = null;
if (velMatches.length >= 3) {
  const isHiFirst = hiV.length <= midV.length;
  const mHi = isHiFirst ? velMatches.map(m => m.a) : velMatches.map(m => m.b);
  const mMid = isHiFirst ? velMatches.map(m => m.b) : velMatches.map(m => m.a);
  velMatchedR1 = hierarchicalDL(mHi, "Matched Vmax>150");
  velMatchedR2 = hierarchicalDL(mMid, "Matched Vmax 80-150");
}

const v1 = summarizeSplit("Distance", precR, hfR, distMatchedR1, distMatchedR2, jkDist);
const v2 = summarizeSplit("Inclination", hiR, miR, incMatchedR1, incMatchedR2, jkInc);
const v3 = summarizeSplit("Velocity", hvR, mvR, velMatchedR1, velMatchedR2, jkVel);

log("");

const allResolved = v1 === "RESOLVED" && v2 === "RESOLVED" && v3 === "RESOLVED";
const anyUnresolved = v1 === "UNRESOLVED" || v2 === "UNRESOLVED" || v3 === "UNRESOLVED";

if (allResolved) {
  log("  OVERALL: ALL THREE SPLITS RESOLVED.");
  log("  a0 is practically robust — remaining splits are sample composition effects.");
  log("  The precise value is limited only by irreducible tau (~0.29 dex).");
  results.overall = "ALL_RESOLVED";
} else if (anyUnresolved) {
  log("  OVERALL: AT LEAST ONE SPLIT UNRESOLVED.");
  log("  The precise a0 value is still limited by identified systematics.");
  log("  The transition scale exists but its exact value is method-dependent.");
  results.overall = "UNRESOLVED";
} else {
  log("  OVERALL: MIXED / INCONCLUSIVE.");
  log("  Splits reduced but not fully eliminated.");
  results.overall = "INCONCLUSIVE";
}

log("");
log("  RECOMMENDED STATEMENT:");
log("");
if (allResolved) {
  log("  \"The acceleration transition scale a0 = 3633 (km/s)^2/kpc = 1.18e-10 m/s^2");
  log("   is observationally robust. All three remaining systematic splits (distance,");
  log("   inclination, velocity) are resolved by matched-sample analysis, confirming");
  log("   they were confounding effects rather than genuine a0 variation.\"");
} else {
  log("  \"There is a real, robust transition acceleration scale in galaxy rotation");
  log("   curves. Its defensible value is approximately a0 ~ 3200-4200 (km/s)^2/kpc");
  log("   = (1.0-1.4) x 10^-10 m/s^2. The precise value remains limited by");
  log("   residual systematics tied to distance, inclination, and/or galaxy mass.\"");
}

const outPath = path.join(__dirname, '..', 'public', 'phase6-matched-results.json');
fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
log("");
log("  Results saved to: public/phase6-matched-results.json");
log("");
log("================================================================================");
log("  PHASE 6 COMPLETE");
log("================================================================================");
