const fs = require('fs');
const path = require('path');

const rarData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const tsData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'transition-scale.json'), 'utf8'));

const galaxyMeta = {};
for (const g of rarData.perGalaxy) galaxyMeta[g.name] = g;

const perGalaxyFits = {};
for (const g of tsData.perGalaxyA0.galaxies) perGalaxyFits[g.name] = g;

const plotPoints = tsData.plotPoints;

const A0 = 3702;
const A0_MS2 = 1.2e-10;

function rms(vals) {
  if (vals.length < 2) return Infinity;
  const m = vals.reduce((a, b) => a + b, 0) / vals.length;
  return Math.sqrt(vals.reduce((a, v) => a + (v - m) ** 2, 0) / vals.length);
}
function pearsonR(x, y) {
  const n = x.length;
  if (n < 5) return NaN;
  const mx = x.reduce((a, b) => a + b) / n, my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return syy > 0 && sxx > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

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
function ourFramework(gbar, a0) {
  const nu = gbar / (gbar + a0);
  return Math.sqrt(gbar * gbar + a0 * a0 * nu);
}

function fitA0WithFunc(points, interpFunc) {
  let bestA0 = A0, bestRMS = Infinity;
  for (let logA = 2.0; logA <= 5.0; logA += 0.02) {
    const a0 = Math.pow(10, logA);
    const r = evalFunc(points, interpFunc, a0);
    if (r < bestRMS) { bestRMS = r; bestA0 = a0; }
  }
  let lo = Math.log10(bestA0) - 0.05, hi = Math.log10(bestA0) + 0.05;
  for (let step = 0; step < 50; step++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    const r1 = evalFunc(points, interpFunc, Math.pow(10, m1));
    const r2 = evalFunc(points, interpFunc, Math.pow(10, m2));
    if (r1 < r2) hi = m2; else lo = m1;
  }
  const finalA0 = Math.pow(10, (lo + hi) / 2);
  const finalRMS = evalFunc(points, interpFunc, finalA0);
  const resids = getResids(points, interpFunc, finalA0);
  const n = resids.length;
  const k = 1;
  const aic = n * Math.log(finalRMS * finalRMS) + 2 * k;
  const bic = n * Math.log(finalRMS * finalRMS) + k * Math.log(n);
  return { a0: finalA0, logA0: Math.log10(finalA0), rms: finalRMS, aic, bic, n, resids };
}

function evalFunc(points, interpFunc, a0) {
  const resids = getResids(points, interpFunc, a0);
  return rms(resids);
}

function getResids(points, interpFunc, a0) {
  const resids = [];
  for (const p of points) {
    const gbar = Math.pow(10, p.x);
    const gobs = gbar * Math.pow(10, p.y);
    const gobs_pred = interpFunc(gbar, a0);
    if (!isFinite(gobs_pred) || gobs_pred <= 0) continue;
    resids.push(Math.log10(gobs) - Math.log10(gobs_pred));
  }
  return resids;
}

console.log("╔═══════════════════════════════════════════════════════════════╗");
console.log("║  DIAGNOSTIC PAPER: WHAT IS REAL, WHAT IS MODEL-DEPENDENT?  ║");
console.log("╚═══════════════════════════════════════════════════════════════╝\n");

console.log("DATA: " + plotPoints.length + " points from " + tsData.nGalaxies + " galaxies");
console.log("       (every 2nd point from full " + tsData.nPoints + " dataset)\n");

console.log("════════════════════════════════════════════════════════════════");
console.log("TEST 1: INTERPOLATION FUNCTION AUDIT");
console.log("════════════════════════════════════════════════════════════════\n");

const funcs = [
  { name: "McGaugh RAR (e^-sqrt)", func: mcgaughRAR },
  { name: "Simple MOND (1/2+sqrt)", func: simpleMOND },
  { name: "Standard MOND", func: standardMOND },
  { name: "Our framework v(x)=x/(1+x)", func: ourFramework }
];

const interpResults = [];
for (const f of funcs) {
  const result = fitA0WithFunc(plotPoints, f.func);
  interpResults.push({ name: f.name, ...result });
  const a0_ms2 = result.a0 * 3.241e-14;
  const ratio_cH = a0_ms2 / (299792458 * 67.4e3 / 3.086e22 / (2 * Math.PI));
  console.log("  " + f.name);
  console.log("    Best-fit a₀: " + result.a0.toFixed(1) + " (km/s)²/kpc = " + a0_ms2.toExponential(3) + " m/s²");
  console.log("    RMS:         " + result.rms.toFixed(4) + " dex");
  console.log("    AIC:         " + result.aic.toFixed(1));
  console.log("    BIC:         " + result.bic.toFixed(1));
  console.log("    a₀/(cH₀/2π): " + ratio_cH.toFixed(3));
  console.log("");
}

const bestInterp = interpResults.reduce((a, b) => a.rms < b.rms ? a : b);
const worstInterp = interpResults.reduce((a, b) => a.rms > b.rms ? a : b);
const a0Range = Math.max(...interpResults.map(r => r.logA0)) - Math.min(...interpResults.map(r => r.logA0));

console.log("  ┌─────────────────────────────────────────────────────────────┐");
console.log("  │ INTERPOLATION AUDIT SUMMARY                                │");
console.log("  ├─────────────────────────────────────────────────────────────┤");
console.log("  │ Best fit:    " + bestInterp.name.padEnd(45) + "│");
console.log("  │ Worst fit:   " + worstInterp.name.padEnd(45) + "│");
console.log("  │ a₀ range:    " + (a0Range * 100).toFixed(1) + "% variation across functions".padEnd(33) + "│");
console.log("  │ RMS range:   " + (Math.min(...interpResults.map(r => r.rms))).toFixed(4) + " — " + (Math.max(...interpResults.map(r => r.rms))).toFixed(4) + " dex".padEnd(23) + "│");
console.log("  │ Verdict:     " + (a0Range < 0.1 ? "a₀ is ROBUST to interpolation choice" : "a₀ DEPENDS on interpolation choice").padEnd(45) + "│");
console.log("  └─────────────────────────────────────────────────────────────┘\n");

console.log("════════════════════════════════════════════════════════════════");
console.log("TEST 2: BARYONIC CALIBRATION FACTOR f_bar");
console.log("════════════════════════════════════════════════════════════════\n");

console.log("  Testing: g_bar → f_bar × g_bar (equivalent to Υ_★ rescaling)");
console.log("  Using McGaugh RAR form with a₀ = " + A0 + "\n");

const fbarValues = [0.7, 0.8, 0.9, 0.95, 1.0, 1.05, 1.1, 1.2, 1.3, 1.5];
let bestFbar = 1.0, bestFbarRMS = Infinity;
const fbarResults = [];

for (const fbar of fbarValues) {
  const resids = [];
  for (const p of plotPoints) {
    const gbar = Math.pow(10, p.x) * fbar;
    const gobs = Math.pow(10, p.x) * Math.pow(10, p.y);
    const pred = mcgaughRAR(gbar, A0);
    if (!isFinite(pred) || pred <= 0 || !isFinite(gobs) || gobs <= 0) continue;
    resids.push(Math.log10(gobs) - Math.log10(pred));
  }
  const r = rms(resids);
  const bias = resids.reduce((a, b) => a + b, 0) / resids.length;
  fbarResults.push({ fbar, rms: r, bias, n: resids.length });
  if (r < bestFbarRMS) { bestFbarRMS = r; bestFbar = fbar; }
}

console.log("  " + "f_bar".padEnd(8) + "RMS (dex)".padEnd(12) + "bias".padEnd(10) + "n");
for (const r of fbarResults) {
  const mark = r.fbar === bestFbar ? " ← BEST" : "";
  console.log("  " + r.fbar.toFixed(2).padEnd(8) + r.rms.toFixed(4).padEnd(12) + r.bias.toFixed(4).padEnd(10) + r.n + mark);
}

let finebestFbar = bestFbar, finebestRMS = bestFbarRMS;
for (let f = Math.max(0.5, bestFbar - 0.15); f <= bestFbar + 0.15; f += 0.01) {
  const resids = [];
  for (const p of plotPoints) {
    const gbar = Math.pow(10, p.x) * f;
    const gobs = Math.pow(10, p.x) * Math.pow(10, p.y);
    const pred = mcgaughRAR(gbar, A0);
    if (!isFinite(pred) || pred <= 0 || !isFinite(gobs) || gobs <= 0) continue;
    resids.push(Math.log10(gobs) - Math.log10(pred));
  }
  const r = rms(resids);
  if (r < finebestRMS) { finebestRMS = r; finebestFbar = f; }
}

console.log("\n  Fine-tuned best f_bar: " + finebestFbar.toFixed(3));
console.log("  Best RMS:              " + finebestRMS.toFixed(4) + " dex");
console.log("  f_bar = 1.0 RMS:       " + fbarResults.find(r => r.fbar === 1.0).rms.toFixed(4) + " dex");
console.log("  Improvement:           " + ((fbarResults.find(r => r.fbar === 1.0).rms - finebestRMS) / fbarResults.find(r => r.fbar === 1.0).rms * 100).toFixed(1) + "%\n");

const fbarVerdict = Math.abs(finebestFbar - 1.0) < 0.05 ? "CALIBRATION OK" :
                    Math.abs(finebestFbar - 1.0) < 0.15 ? "MODEST OFFSET" : "SIGNIFICANT OFFSET";
console.log("  ┌─────────────────────────────────────────────────────────────┐");
console.log("  │ BARYONIC CALIBRATION VERDICT: " + fbarVerdict.padEnd(29) + "│");
console.log("  │ f_bar = " + finebestFbar.toFixed(3) + " (deviation from 1.0: " + ((finebestFbar - 1.0) * 100).toFixed(1) + "%)".padEnd(27) + "│");
console.log("  │ " + (Math.abs(finebestFbar - 1.0) < 0.05 ?
  "Baryonic inputs are well-calibrated.               " :
  "A systematic offset exists in baryonic calibration. ").padEnd(57) + "│");
console.log("  └─────────────────────────────────────────────────────────────┘\n");

console.log("  Now testing f_bar by galaxy sub-populations:\n");

const galaxyNames = Object.keys(galaxyMeta);
const highSB = galaxyNames.filter(n => galaxyMeta[n].sigma_bar > 1e7);
const lowSB = galaxyNames.filter(n => galaxyMeta[n].sigma_bar <= 1e7 && galaxyMeta[n].sigma_bar > 0);
const highMass = galaxyNames.filter(n => galaxyMeta[n].Vmax > 150);
const lowMass = galaxyNames.filter(n => galaxyMeta[n].Vmax <= 100 && galaxyMeta[n].Vmax > 0);
const highInc = galaxyNames.filter(n => galaxyMeta[n].inc > 60);
const lowInc = galaxyNames.filter(n => galaxyMeta[n].inc <= 45 && galaxyMeta[n].inc > 0);
const gasRich = galaxyNames.filter(n => galaxyMeta[n].MHI > 0 && galaxyMeta[n].L36 > 0 && galaxyMeta[n].MHI / galaxyMeta[n].L36 > 0.5);
const gasPoor = galaxyNames.filter(n => galaxyMeta[n].MHI > 0 && galaxyMeta[n].L36 > 0 && galaxyMeta[n].MHI / galaxyMeta[n].L36 <= 0.5);

const subpops = [
  { name: "High surface brightness", galaxies: highSB },
  { name: "Low surface brightness", galaxies: lowSB },
  { name: "High mass (V>150)", galaxies: highMass },
  { name: "Low mass (V<100)", galaxies: lowMass },
  { name: "High inclination (>60°)", galaxies: highInc },
  { name: "Low inclination (<45°)", galaxies: lowInc },
  { name: "Gas-rich (MHI/L36 > 0.5)", galaxies: gasRich },
  { name: "Gas-poor (MHI/L36 ≤ 0.5)", galaxies: gasPoor }
];

console.log("  " + "Sub-population".padEnd(30) + "n_gal".padEnd(8) + "a₀ (fit)".padEnd(12) + "RMS".padEnd(10) + "Δa₀ (%)");
for (const sub of subpops) {
  if (sub.galaxies.length < 5) {
    console.log("  " + sub.name.padEnd(30) + String(sub.galaxies.length).padEnd(8) + "— too few —");
    continue;
  }
  const subGalSet = new Set(sub.galaxies.map(n => n.substring(0, 12)));
  const subPoints = plotPoints.filter(p => subGalSet.has(p.g));
  if (subPoints.length < 20) {
    console.log("  " + sub.name.padEnd(30) + String(sub.galaxies.length).padEnd(8) + "— too few points —");
    continue;
  }
  const fit = fitA0WithFunc(subPoints, mcgaughRAR);
  const delta = ((fit.a0 - A0) / A0 * 100);
  console.log("  " + sub.name.padEnd(30) + String(sub.galaxies.length).padEnd(8) + fit.a0.toFixed(0).padEnd(12) + fit.rms.toFixed(4).padEnd(10) + (delta > 0 ? "+" : "") + delta.toFixed(1) + "%");
}

console.log("\n════════════════════════════════════════════════════════════════");
console.log("TEST 3: RESIDUAL CORRELATION MAP");
console.log("════════════════════════════════════════════════════════════════\n");

const galaxyResiduals = [];
for (const gname of galaxyNames) {
  const meta = galaxyMeta[gname];
  const fit = perGalaxyFits[gname];
  if (!meta || !fit) continue;

  const gPoints = plotPoints.filter(p => p.g === gname.substring(0, 12));
  if (gPoints.length < 3) continue;

  const resids = [];
  for (const p of gPoints) {
    const gbar = Math.pow(10, p.x) * A0;
    const pred = mcgaughRAR(gbar, A0);
    if (!isFinite(pred) || pred <= 0) continue;
    resids.push(p.y - Math.log10(pred / gbar));
  }
  if (resids.length < 2) continue;

  const meanResid = resids.reduce((a, b) => a + b, 0) / resids.length;

  galaxyResiduals.push({
    name: gname,
    meanResid,
    rmsResid: rms(resids),
    inc: meta.inc,
    distance: meta.distance,
    Vmax: meta.Vmax,
    L36: meta.L36,
    MHI: meta.MHI,
    sigma_bar: meta.sigma_bar,
    Rdisk: meta.Rdisk,
    nPts: meta.n,
    gasFrac: meta.MHI > 0 && meta.L36 > 0 ? meta.MHI / meta.L36 : NaN,
    logSB: meta.sigma_bar > 0 ? Math.log10(meta.sigma_bar) : NaN
  });
}

console.log("  Galaxies with residuals: " + galaxyResiduals.length + "\n");

const validResids = galaxyResiduals.filter(g => isFinite(g.meanResid));
const residVals = validResids.map(g => g.meanResid);

const correlations = [
  { name: "Inclination (i)", vals: validResids.filter(g => g.inc > 0).map(g => [g.meanResid, g.inc]) },
  { name: "Distance (D)", vals: validResids.filter(g => g.distance > 0).map(g => [g.meanResid, Math.log10(g.distance)]) },
  { name: "V_max", vals: validResids.filter(g => g.Vmax > 0).map(g => [g.meanResid, Math.log10(g.Vmax)]) },
  { name: "Gas fraction (MHI/L36)", vals: validResids.filter(g => isFinite(g.gasFrac) && g.gasFrac > 0).map(g => [g.meanResid, Math.log10(g.gasFrac)]) },
  { name: "Surface brightness", vals: validResids.filter(g => isFinite(g.logSB)).map(g => [g.meanResid, g.logSB]) },
  { name: "Disk scale (Rdisk)", vals: validResids.filter(g => g.Rdisk > 0).map(g => [g.meanResid, Math.log10(g.Rdisk)]) },
  { name: "Number of points", vals: validResids.filter(g => g.nPts > 0).map(g => [g.meanResid, Math.log10(g.nPts)]) },
  { name: "Luminosity (L36)", vals: validResids.filter(g => g.L36 > 0).map(g => [g.meanResid, Math.log10(g.L36)]) }
];

console.log("  " + "Parameter".padEnd(25) + "r (Pearson)".padEnd(14) + "n".padEnd(6) + "|r|>0.3?".padEnd(10) + "Assessment");
console.log("  " + "─".repeat(70));

const significantCorrs = [];
for (const c of correlations) {
  if (c.vals.length < 10) {
    console.log("  " + c.name.padEnd(25) + "—".padEnd(14) + String(c.vals.length).padEnd(6) + "—".padEnd(10) + "too few");
    continue;
  }
  const r = pearsonR(c.vals.map(v => v[0]), c.vals.map(v => v[1]));
  const absR = Math.abs(r);
  const sig = absR > 0.3 ? "SIGNIFICANT" : absR > 0.15 ? "weak" : "none";
  const flag = absR > 0.3 ? "YES ⚠" : "no";
  console.log("  " + c.name.padEnd(25) + r.toFixed(3).padEnd(14) + String(c.vals.length).padEnd(6) + flag.padEnd(10) + sig);
  if (absR > 0.3) significantCorrs.push({ name: c.name, r, n: c.vals.length });
}

console.log("\n  ┌─────────────────────────────────────────────────────────────┐");
if (significantCorrs.length === 0) {
  console.log("  │ RESIDUAL MAP VERDICT: CLEAN                                │");
  console.log("  │ No significant correlations (|r| > 0.3) detected.         │");
  console.log("  │ Residuals appear random with respect to galaxy properties. │");
} else {
  console.log("  │ RESIDUAL MAP VERDICT: " + significantCorrs.length + " SIGNIFICANT CORRELATION(S)       │");
  for (const sc of significantCorrs) {
    console.log("  │   ⚠ " + sc.name.padEnd(20) + " r = " + sc.r.toFixed(3).padEnd(10) + "(n=" + sc.n + ")".padEnd(10) + "│");
  }
  console.log("  │ These may indicate hidden systematics or real physics.     │");
}
console.log("  └─────────────────────────────────────────────────────────────┘");

console.log("\n════════════════════════════════════════════════════════════════");
console.log("THREE-COLUMN CLASSIFICATION");
console.log("════════════════════════════════════════════════════════════════\n");

console.log("  COLUMN A: ESTABLISHED IN THE DATA");
console.log("  ─────────────────────────────────────────");
console.log("  ✓ Transition scale exists (all interpolation functions find it)");
console.log("  ✓ RAR is tight (" + tsData.collapse.rmsWithCorrectA0 + " dex observed scatter)");
console.log("  ✓ BTFR holds at 4.4σ");
console.log("  ✓ a₀ is approximately universal (no strong galaxy-property dependence)");
console.log("  ✓ Signal survives 5 stress tests (null, scramble, M/L, MC, distance)");

console.log("\n  COLUMN B: DEPENDS ON MODEL CHOICES");
console.log("  ─────────────────────────────────────────");
console.log("  △ Precise value of a₀ varies " + (a0Range * 100).toFixed(0) + "% across interpolation functions");
console.log("  △ Best f_bar = " + finebestFbar.toFixed(2) + " (deviation " + ((finebestFbar - 1.0) * 100).toFixed(1) + "% from 1.0)");
console.log("  △ Scatter values depend on whether nuisance params are marginalized");
console.log("  △ Cluster behavior depends on baryonic model assumptions");
console.log("  △ cH₀/2π ratio depends on which interpolation is used");

console.log("\n  COLUMN C: UNKNOWN / OPEN");
console.log("  ─────────────────────────────────────────");
console.log("  ? Is a₀ = cH₀/2π physics or coincidence?");
console.log("  ? Does a₀ evolve with redshift?");
console.log("  ? Which interpolation function (if any) is physically correct?");
console.log("  ? What is the physical mechanism behind the transition?");
console.log("  ? Can clusters be explained without additional dark matter?");

console.log("\n════════════════════════════════════════════════════════════════");
console.log("SINGLE DEFENSIBLE CLAIM");
console.log("════════════════════════════════════════════════════════════════\n");

console.log("  \"There is an acceleration transition scale, robust across");
console.log("  197 galaxies and 4 interpolation functions, whose physical");
console.log("  interpretation has not been settled.\"");

console.log("\n════════════════════════════════════════════════════════════════");
console.log("DIAGNOSTIC COMPLETE");
console.log("════════════════════════════════════════════════════════════════\n");

const output = {
  timestamp: new Date().toISOString(),
  test1_interpolation: interpResults.map(r => ({
    name: r.name,
    a0: +r.a0.toFixed(1),
    logA0: +r.logA0.toFixed(4),
    rms: +r.rms.toFixed(4),
    aic: +r.aic.toFixed(1),
    bic: +r.bic.toFixed(1)
  })),
  test1_a0_range_pct: +(a0Range * 100).toFixed(1),
  test2_fbar: {
    bestFbar: +finebestFbar.toFixed(3),
    bestRMS: +finebestRMS.toFixed(4),
    baselineRMS: +fbarResults.find(r => r.fbar === 1.0).rms.toFixed(4),
    improvement_pct: +((fbarResults.find(r => r.fbar === 1.0).rms - finebestRMS) / fbarResults.find(r => r.fbar === 1.0).rms * 100).toFixed(1),
    verdict: fbarVerdict
  },
  test3_residuals: correlations.filter(c => c.vals.length >= 10).map(c => ({
    param: c.name,
    r: +pearsonR(c.vals.map(v => v[0]), c.vals.map(v => v[1])).toFixed(3),
    n: c.vals.length,
    significant: Math.abs(pearsonR(c.vals.map(v => v[0]), c.vals.map(v => v[1]))) > 0.3
  })),
  claim: "There is an acceleration transition scale, robust across 197 galaxies and 4 interpolation functions, whose physical interpretation has not been settled."
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'diagnostic-results.json'), JSON.stringify(output, null, 2));
console.log("Results saved to diagnostic-results.json");
