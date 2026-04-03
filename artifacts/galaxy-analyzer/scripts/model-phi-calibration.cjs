const fs = require('fs');
const path = require('path');

const A0_OBS = 3702;
const G_KPC = 4.3009e-6;
const UPSILON_DISK = 0.5, UPSILON_BULGE = 0.7;
const C_MS = 299792458;
const H0_PER_S = 70 / 3.0857e19;
const CH0_MS2 = C_MS * H0_PER_S;
const CH0_KPC = CH0_MS2 * 3.0857e16;

function rms(vals) {
  if (vals.length < 2) return Infinity;
  const m = vals.reduce((a, b) => a + b, 0) / vals.length;
  return Math.sqrt(vals.reduce((a, v) => a + (v - m) ** 2, 0) / vals.length);
}
function median(vals) {
  const s = [...vals].sort((a, b) => a - b);
  return s[Math.floor(s.length / 2)];
}
function pearsonR(x, y) {
  const n = x.length;
  if (n < 5) return 0;
  const mx = x.reduce((a, b) => a + b) / n, my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return syy > 0 && sxx > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}
function linreg(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b) / n, my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; }
  const slope = sxx > 0 ? sxy / sxx : 0;
  const intercept = my - slope * mx;
  return { slope, intercept };
}

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  if (y < 1e-10) return gbar;
  const denom = 1 - Math.exp(-Math.sqrt(y));
  if (denom < 1e-10) return gbar;
  return gbar / denom;
}

function fitA0(points) {
  let bestA0 = A0_OBS, bestRMS = Infinity;
  for (let logA = 2.5; logA <= 4.5; logA += 0.02) {
    const a0 = Math.pow(10, logA);
    const r = rms(points.map(p => Math.log10(p.gObs) - Math.log10(mcgaughRAR(p.gBar, a0))));
    if (r < bestRMS) { bestRMS = r; bestA0 = a0; }
  }
  let lo = Math.log10(bestA0) - 0.05, hi = Math.log10(bestA0) + 0.05;
  for (let step = 0; step < 50; step++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    const r1 = rms(points.map(p => Math.log10(p.gObs) - Math.log10(mcgaughRAR(p.gBar, Math.pow(10, m1)))));
    const r2 = rms(points.map(p => Math.log10(p.gObs) - Math.log10(mcgaughRAR(p.gBar, Math.pow(10, m2)))));
    if (r1 < r2) hi = m2; else lo = m1;
  }
  const a0 = Math.pow(10, (lo + hi) / 2);
  const r = rms(points.map(p => Math.log10(p.gObs) - Math.log10(mcgaughRAR(p.gBar, a0))));
  return { a0, logA0: Math.log10(a0), rms: r };
}

console.log("=======================================================");
console.log(" MODEL PHI CALIBRATION");
console.log(" Making Phi concrete: Phi = 1 + alpha*log(Sigma/Sigma0)");
console.log("=======================================================\n");

const rarJson = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const galaxyLookup = {};
for (const g of rarJson.perGalaxy) galaxyLookup[g.name] = g;

const galaxies = [];
const sparcDir = '/tmp/rotmod';
for (const file of fs.readdirSync(sparcDir).filter(f => f.endsWith('.dat'))) {
  const name = file.replace('_rotmod.dat', '');
  const info = galaxyLookup[name];
  if (!info) continue;
  
  const lines = fs.readFileSync(path.join(sparcDir, file), 'utf8').trim().split('\n');
  const pts = [];
  let maxR = 0;
  for (const line of lines) {
    if (line.trim().startsWith('#') || line.trim() === '') continue;
    const parts = line.trim().split(/\s+/).map(Number);
    if (parts.length < 6) continue;
    const [r, vobs, evobs, vgas, vdisk, vbul] = parts;
    if (r <= 0 || vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * vdisk * Math.abs(vdisk) + UPSILON_BULGE * (vbul || 0) * Math.abs(vbul || 0) + vgas * Math.abs(vgas);
    const gObs = vobs * vobs / r;
    const gBar = Math.abs(vBarSq) / r;
    if (gBar > 0 && gObs > 0 && isFinite(gBar) && isFinite(gObs)) {
      pts.push({ r, gObs, gBar });
      if (r > maxR) maxR = r;
    }
  }
  
  if (pts.length < 5) continue;
  
  const logRange = Math.log10(Math.max(...pts.map(p => p.gBar))) - Math.log10(Math.min(...pts.map(p => p.gBar)));
  if (logRange < 0.3) continue;
  
  const sigma_bar = info.sigma_bar || null;
  const maxV = info.maxV || Math.max(...pts.map(p => Math.sqrt(p.gObs * pts[0].r)));
  
  const fit = fitA0(pts);
  if (fit.rms > 0.25 || fit.logA0 < 2.5 || fit.logA0 > 4.5) continue;
  
  galaxies.push({
    name,
    nPts: pts.length,
    logRange,
    sigma_bar,
    logSigma: sigma_bar && sigma_bar > 0 ? Math.log10(sigma_bar) : null,
    maxV,
    logMaxV: Math.log10(maxV),
    maxR,
    a0_fit: fit.a0,
    logA0: fit.logA0,
    rms: fit.rms,
    phi: fit.a0 / A0_OBS,
    logPhi: fit.logA0 - Math.log10(A0_OBS),
    points: pts
  });
}

console.log("Galaxies with good a0 fits: " + galaxies.length);
console.log("Global a0 (known): " + A0_OBS + " (km/s)^2/kpc = " + (A0_OBS * 1e-10 / (3702 * 1e-10 / 1.2e-10)).toExponential(2) + " m/s^2");
console.log("cH0: " + CH0_MS2.toExponential(3) + " m/s^2");
console.log("epsilon = a0/cH0: " + (1.2e-10 / CH0_MS2).toFixed(4));
console.log("1/(2pi): " + (1 / (2 * Math.PI)).toFixed(4));
console.log();

const withSigma = galaxies.filter(g => g.logSigma !== null && isFinite(g.logSigma));
console.log("Galaxies with sigma_bar: " + withSigma.length);

const logPhiVals = withSigma.map(g => g.logPhi);
const logSigmaVals = withSigma.map(g => g.logSigma);

console.log("\nPhi statistics (all galaxies):");
console.log("  median(Phi): " + median(galaxies.map(g => g.phi)).toFixed(4));
console.log("  median(logPhi): " + median(galaxies.map(g => g.logPhi)).toFixed(4));
console.log("  RMS(logPhi): " + rms(galaxies.map(g => g.logPhi)).toFixed(4));

console.log("\n--- FITTING Phi(Sigma_bar) ---");

const rSigma = pearsonR(logSigmaVals, logPhiVals);
console.log("Correlation r(logSigma, logPhi): " + rSigma.toFixed(4));

const fitSigma = linreg(logSigmaVals, logPhiVals);
console.log("Linear fit: logPhi = " + fitSigma.slope.toFixed(4) + " * logSigma + " + fitSigma.intercept.toFixed(4));

const SIGMA0 = Math.pow(10, -fitSigma.intercept / fitSigma.slope);
const ALPHA = fitSigma.slope;
console.log("\nPHI MODEL: Phi = 1 + alpha*log10(Sigma/Sigma0)");
console.log("  alpha = " + ALPHA.toFixed(4));
console.log("  Sigma0 = " + SIGMA0.toFixed(1) + " M_sun/pc^2");
console.log("  log10(Sigma0) = " + Math.log10(SIGMA0).toFixed(3));

console.log("\n--- TESTING: Does Phi improve the collapse? ---");

const allPtsUniversal = [];
const allPtsPhiCorrected = [];

for (const g of withSigma) {
  const a0_universal = A0_OBS;
  const logPhi_pred = ALPHA * (g.logSigma - Math.log10(SIGMA0));
  const phi_pred = Math.pow(10, logPhi_pred);
  const a0_corrected = A0_OBS * phi_pred;
  
  for (const p of g.points) {
    const resid_univ = Math.log10(p.gObs) - Math.log10(mcgaughRAR(p.gBar, a0_universal));
    const resid_corr = Math.log10(p.gObs) - Math.log10(mcgaughRAR(p.gBar, a0_corrected));
    allPtsUniversal.push(resid_univ);
    allPtsPhiCorrected.push(resid_corr);
  }
}

const rmsUniversal = rms(allPtsUniversal);
const rmsCorrected = rms(allPtsPhiCorrected);
const improvement = ((rmsUniversal - rmsCorrected) / rmsUniversal * 100);

console.log("  RMS with universal a0: " + rmsUniversal.toFixed(4) + " dex");
console.log("  RMS with Phi-corrected a0: " + rmsCorrected.toFixed(4) + " dex");
console.log("  Improvement: " + improvement.toFixed(2) + "%");

console.log("\n--- CROSS-VALIDATION ---");

let cvSum = 0, cvN = 0;
for (let fold = 0; fold < 5; fold++) {
  const test = withSigma.filter((_, i) => i % 5 === fold);
  const train = withSigma.filter((_, i) => i % 5 !== fold);
  
  const trainFit = linreg(train.map(g => g.logSigma), train.map(g => g.logPhi));
  const trainSigma0 = Math.pow(10, -trainFit.intercept / trainFit.slope);
  
  let testResidsUniv = [], testResidsCorr = [];
  for (const g of test) {
    const logPhi_pred = trainFit.slope * (g.logSigma - Math.log10(trainSigma0));
    const a0c = A0_OBS * Math.pow(10, logPhi_pred);
    for (const p of g.points) {
      testResidsUniv.push(Math.log10(p.gObs) - Math.log10(mcgaughRAR(p.gBar, A0_OBS)));
      testResidsCorr.push(Math.log10(p.gObs) - Math.log10(mcgaughRAR(p.gBar, a0c)));
    }
  }
  const rU = rms(testResidsUniv), rC = rms(testResidsCorr);
  const imp = (rU - rC) / rU * 100;
  cvSum += imp;
  cvN++;
  console.log("  Fold " + fold + ": universal=" + rU.toFixed(4) + " corrected=" + rC.toFixed(4) + " improvement=" + imp.toFixed(2) + "%");
}
console.log("  Mean CV improvement: " + (cvSum / cvN).toFixed(2) + "%");

console.log("\n--- OTHER CANDIDATE VARIABLES FOR PHI ---");

const withMaxV = galaxies.filter(g => g.logMaxV > 0);
const rMaxV = pearsonR(withMaxV.map(g => g.logMaxV), withMaxV.map(g => g.logPhi));
const fitMaxV = linreg(withMaxV.map(g => g.logMaxV), withMaxV.map(g => g.logPhi));
console.log("r(logVmax, logPhi): " + rMaxV.toFixed(4) + " slope=" + fitMaxV.slope.toFixed(4));

const withMaxR = galaxies.filter(g => g.maxR > 0);
const rMaxR = pearsonR(withMaxR.map(g => Math.log10(g.maxR)), withMaxR.map(g => g.logPhi));
console.log("r(logRmax, logPhi): " + rMaxR.toFixed(4));

const withNPts = galaxies.filter(g => g.nPts >= 5);
const rNPts = pearsonR(withNPts.map(g => Math.log10(g.nPts)), withNPts.map(g => g.logPhi));
console.log("r(logNpts, logPhi): " + rNPts.toFixed(4));

console.log("\n--- WHAT HAPPENS IF WE VARY 'FEEDBACK' (proxy: Vmax bins) ---");

const bins = [
  { label: "Dwarfs (Vmax<80)", gals: galaxies.filter(g => g.maxV < 80) },
  { label: "Intermediate (80-150)", gals: galaxies.filter(g => g.maxV >= 80 && g.maxV < 150) },
  { label: "Massive (Vmax>150)", gals: galaxies.filter(g => g.maxV >= 150) }
];

for (const bin of bins) {
  if (bin.gals.length < 3) continue;
  const a0s = bin.gals.map(g => g.logA0);
  const med = median(a0s);
  const scatter = rms(a0s.map(v => v - median(a0s)));
  console.log("  " + bin.label + ": n=" + bin.gals.length + " median(logA0)=" + med.toFixed(3) + " (a0=" + Math.pow(10, med).toFixed(0) + ") scatter=" + scatter.toFixed(3) + " dex");
}

console.log("\n--- DOES a0 MOVE WITH GALAXY MASS? (feedback proxy) ---");

const sortedByV = [...galaxies].sort((a, b) => a.maxV - b.maxV);
const quartiles = [
  { label: "Q1 (lightest)", gals: sortedByV.slice(0, Math.floor(sortedByV.length / 4)) },
  { label: "Q2", gals: sortedByV.slice(Math.floor(sortedByV.length / 4), Math.floor(sortedByV.length / 2)) },
  { label: "Q3", gals: sortedByV.slice(Math.floor(sortedByV.length / 2), Math.floor(3 * sortedByV.length / 4)) },
  { label: "Q4 (heaviest)", gals: sortedByV.slice(Math.floor(3 * sortedByV.length / 4)) }
];

console.log("  " + "Quartile".padEnd(20) + "n".padEnd(6) + "med(logA0)".padEnd(14) + "a0".padEnd(10) + "scatter".padEnd(10) + "med(Vmax)");
for (const q of quartiles) {
  const a0s = q.gals.map(g => g.logA0);
  const m = median(a0s);
  const s = rms(a0s.map(v => v - m));
  const mV = median(q.gals.map(g => g.maxV));
  console.log("  " + q.label.padEnd(20) + q.gals.length.toString().padEnd(6) + m.toFixed(3).padEnd(14) + Math.pow(10, m).toFixed(0).padEnd(10) + s.toFixed(3).padEnd(10) + mV.toFixed(0));
}

const rVmaxA0 = pearsonR(galaxies.map(g => g.logMaxV), galaxies.map(g => g.logA0));
console.log("\n  Correlation r(logVmax, logA0): " + rVmaxA0.toFixed(4));
console.log("  " + (Math.abs(rVmaxA0) < 0.15 ? "WEAK — a0 does NOT systematically move with galaxy mass" : Math.abs(rVmaxA0) < 0.3 ? "MODERATE — some mass dependence" : "STRONG — a0 moves with mass!"));

console.log("\n=======================================================");
console.log(" FINAL MODEL SPECIFICATION");
console.log("=======================================================\n");

const epsilon = 1.2e-10 / CH0_MS2;
console.log("  a0 = epsilon * cH0 * Phi(Sigma_bar)\n");
console.log("  epsilon = " + epsilon.toFixed(4) + " ≈ 1/(2pi) = " + (1 / (2 * Math.PI)).toFixed(4));
console.log("  cH0 = " + CH0_MS2.toExponential(3) + " m/s^2");
console.log("  Phi = 10^[" + ALPHA.toFixed(4) + " * (log(Sigma) - " + Math.log10(SIGMA0).toFixed(3) + ")]");
console.log("  Sigma0 = " + SIGMA0.toFixed(1) + " M_sun/pc^2");
console.log("  alpha = " + ALPHA.toFixed(4));
console.log("\n  Or equivalently:");
console.log("  Phi ≈ 1 + " + (ALPHA * Math.log(10)).toFixed(3) + " * ln(Sigma/Sigma0)");
console.log("\n  For typical galaxy (Sigma_bar ~ 100 M_sun/pc^2):");
const typPhi = Math.pow(10, ALPHA * (Math.log10(100) - Math.log10(SIGMA0)));
console.log("  Phi = " + typPhi.toFixed(4));
console.log("  a0 = " + (epsilon * CH0_MS2 * typPhi).toExponential(3) + " m/s^2");

console.log("\n  IMPROVEMENT from Phi correction:");
console.log("  Without: " + rmsUniversal.toFixed(4) + " dex");
console.log("  With:    " + rmsCorrected.toFixed(4) + " dex");
console.log("  Gain:    " + improvement.toFixed(2) + "%");
console.log("  CV Gain: " + (cvSum / cvN).toFixed(2) + "%");

const output = {
  timestamp: new Date().toISOString(),
  model: {
    equation: "a0 = epsilon * c * H0 * Phi(Sigma_bar)",
    epsilon: +epsilon.toFixed(4),
    oneOver2Pi: +(1 / (2 * Math.PI)).toFixed(4),
    cH0_ms2: +CH0_MS2.toExponential(3),
    cH0_kpc: +CH0_KPC.toFixed(2)
  },
  phi: {
    formula: "Phi = 10^[alpha * (log10(Sigma_bar) - log10(Sigma0))]",
    formulaApprox: "Phi ≈ 1 + alpha_ln * ln(Sigma_bar/Sigma0)",
    alpha: +ALPHA.toFixed(4),
    alphaLn: +(ALPHA * Math.log(10)).toFixed(4),
    sigma0: +SIGMA0.toFixed(1),
    logSigma0: +Math.log10(SIGMA0).toFixed(3),
    correlationR: +rSigma.toFixed(4)
  },
  calibration: {
    nGalaxies: withSigma.length,
    rmsUniversal: +rmsUniversal.toFixed(4),
    rmsCorrected: +rmsCorrected.toFixed(4),
    improvement: +improvement.toFixed(2),
    cvImprovement: +(cvSum / cvN).toFixed(2)
  },
  otherVariables: {
    rVmax: +rMaxV.toFixed(4),
    rRmax: +rMaxR.toFixed(4),
    rNpts: +rNPts.toFixed(4)
  },
  massQuartiles: quartiles.map(q => {
    const a0s = q.gals.map(g => g.logA0);
    const m = median(a0s);
    return {
      label: q.label,
      n: q.gals.length,
      medLogA0: +m.toFixed(3),
      medA0: +Math.pow(10, m).toFixed(0),
      scatter: +rms(a0s.map(v => v - m)).toFixed(3),
      medVmax: +median(q.gals.map(g => g.maxV)).toFixed(0)
    };
  }),
  rVmaxA0: +rVmaxA0.toFixed(4),
  perGalaxyResults: galaxies.slice(0, 80).map(g => ({
    name: g.name,
    a0: +g.a0_fit.toFixed(1),
    logA0: +g.logA0.toFixed(3),
    phi: +g.phi.toFixed(4),
    logPhi: +g.logPhi.toFixed(4),
    sigma_bar: g.sigma_bar,
    logSigma: g.logSigma ? +g.logSigma.toFixed(3) : null,
    maxV: +g.maxV.toFixed(1),
    nPts: g.nPts,
    rms: +g.rms.toFixed(4)
  })),
  verdict: {
    phiDefined: true,
    phiPredictive: improvement > 0.5,
    phiSmall: Math.abs(ALPHA) < 0.5,
    a0MovesWithMass: Math.abs(rVmaxA0) > 0.15,
    summary: ""
  }
};

output.verdict.summary = Math.abs(rVmaxA0) < 0.15
  ? "a0 does NOT systematically depend on galaxy mass — the scale is truly universal. Phi correction gives " + improvement.toFixed(1) + "% improvement: detectable but the model works mostly because epsilon*cH0 already sets the right scale."
  : "a0 shows " + (Math.abs(rVmaxA0) < 0.3 ? "weak" : "moderate") + " dependence on galaxy mass (r=" + rVmaxA0.toFixed(3) + "). Phi correction helps by " + improvement.toFixed(1) + "%. The feedback-mass connection is real but small.";

fs.writeFileSync(path.join(__dirname, '..', 'public', 'model-phi.json'), JSON.stringify(output, null, 2));
console.log("\nSaved to public/model-phi.json");
