const fs = require('fs');
const path = require('path');

const dataPath = path.join(__dirname, '..', 'public', 'fdm-analysis.json');
const data = JSON.parse(fs.readFileSync(dataPath, 'utf8'));

const sparcGals = data.sparc.galaxies;
const simGals = data.simulation[0].galaxies;

function seededRNG(seed) {
  let s = seed;
  return () => { s = (s * 1664525 + 1013904223) & 0x7fffffff; return s / 0x7fffffff; };
}

function linearRegression(x, y) {
  const n = x.length;
  if (n < 3) return { slope: NaN, r: NaN };
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let sxx = 0, sxy = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxx += (x[i] - mx) ** 2;
    sxy += (x[i] - mx) * (y[i] - my);
    syy += (y[i] - my) ** 2;
  }
  const slope = sxy / sxx;
  const r2 = sxx > 0 && syy > 0 ? (sxy * sxy) / (sxx * syy) : 0;
  const r = Math.sign(slope) * Math.sqrt(Math.abs(r2));
  return { slope, r };
}

function resampleSlope(gals, xKey, yKey, rng) {
  const n = gals.length;
  const xs = [], ys = [];
  for (let i = 0; i < n; i++) {
    const idx = Math.floor(rng() * n);
    xs.push(gals[idx][xKey]);
    ys.push(gals[idx][yKey]);
  }
  return linearRegression(xs, ys).slope;
}

const N_BOOT = 10000;
const N_PERM = 10000;

console.log(`\n${'═'.repeat(70)}`);
console.log(`  STATISTICAL SIGNIFICANCE: OBSERVED vs SIMULATED`);
console.log(`${'═'.repeat(70)}\n`);

console.log(`  SPARC: ${sparcGals.length} galaxies, slope=${data.sparc.perGalaxy.slope.toFixed(4)}`);
console.log(`  ΛCDM:  ${simGals.length} mock galaxies, slope=${data.simulation[0].perGalaxy.slope.toFixed(4)}\n`);

console.log(`  ─── Bootstrap (${N_BOOT} iterations) ───\n`);

const rng1 = seededRNG(12345);
const sparcBootSlopes = [];
for (let i = 0; i < N_BOOT; i++) {
  sparcBootSlopes.push(resampleSlope(sparcGals, 'meanLogSigBar', 'meanFDM', rng1));
}
sparcBootSlopes.sort((a, b) => a - b);
const sparcCI = [sparcBootSlopes[Math.floor(N_BOOT * 0.025)], sparcBootSlopes[Math.floor(N_BOOT * 0.975)]];
const sparcMean = sparcBootSlopes.reduce((a, b) => a + b) / N_BOOT;
const sparcSD = Math.sqrt(sparcBootSlopes.reduce((s, v) => s + (v - sparcMean) ** 2, 0) / (N_BOOT - 1));

console.log(`  SPARC bootstrap: mean=${sparcMean.toFixed(5)}, SD=${sparcSD.toFixed(5)}`);
console.log(`  SPARC 95% CI: [${sparcCI[0].toFixed(5)}, ${sparcCI[1].toFixed(5)}]`);

const rng2 = seededRNG(54321);
const simBootSlopes = [];
for (let i = 0; i < N_BOOT; i++) {
  simBootSlopes.push(resampleSlope(simGals, 'logSigGal', 'meanFDM', rng2));
}
simBootSlopes.sort((a, b) => a - b);
const simCI = [simBootSlopes[Math.floor(N_BOOT * 0.025)], simBootSlopes[Math.floor(N_BOOT * 0.975)]];
const simMean = simBootSlopes.reduce((a, b) => a + b) / N_BOOT;
const simSD = Math.sqrt(simBootSlopes.reduce((s, v) => s + (v - simMean) ** 2, 0) / (N_BOOT - 1));

console.log(`  ΛCDM bootstrap:  mean=${simMean.toFixed(5)}, SD=${simSD.toFixed(5)}`);
console.log(`  ΛCDM 95% CI:  [${simCI[0].toFixed(5)}, ${simCI[1].toFixed(5)}]`);

const ciOverlap = sparcCI[0] <= simCI[1] && simCI[0] <= sparcCI[1];
console.log(`\n  CIs overlap? ${ciOverlap ? 'YES' : 'NO — statistically significant'}`);

console.log(`\n  ─── Permutation test (${N_PERM} shuffles) ───\n`);

const realDelta = data.sparc.perGalaxy.slope - data.simulation[0].perGalaxy.slope;
console.log(`  Real Δslope = ${realDelta.toFixed(5)}`);

const pooled = [];
for (const g of sparcGals) pooled.push({ x: g.meanLogSigBar, y: g.meanFDM, source: 'sparc' });
for (const g of simGals) pooled.push({ x: g.logSigGal, y: g.meanFDM, source: 'sim' });

const nSparc = sparcGals.length;
const rng3 = seededRNG(99999);
let moreExtreme = 0;
const permDeltas = [];

for (let i = 0; i < N_PERM; i++) {
  const shuffled = [...pooled];
  for (let j = shuffled.length - 1; j > 0; j--) {
    const k = Math.floor(rng3() * (j + 1));
    [shuffled[j], shuffled[k]] = [shuffled[k], shuffled[j]];
  }
  const groupA = shuffled.slice(0, nSparc);
  const groupB = shuffled.slice(nSparc);
  const regA = linearRegression(groupA.map(g => g.x), groupA.map(g => g.y));
  const regB = linearRegression(groupB.map(g => g.x), groupB.map(g => g.y));
  const delta = regA.slope - regB.slope;
  permDeltas.push(delta);
  if (Math.abs(delta) >= Math.abs(realDelta)) moreExtreme++;
}

const pValue = moreExtreme / N_PERM;
console.log(`  Permutation p-value (two-sided): ${pValue.toFixed(6)} (${moreExtreme}/${N_PERM})`);

permDeltas.sort((a, b) => a - b);
const permMean = permDeltas.reduce((a, b) => a + b) / N_PERM;
const permSD = Math.sqrt(permDeltas.reduce((s, v) => s + (v - permMean) ** 2, 0) / (N_PERM - 1));

console.log(`  Permutation null dist: mean=${permMean.toFixed(5)}, SD=${permSD.toFixed(5)}`);

const permHistBins = [];
const binMin = Math.min(...permDeltas, realDelta) - 0.005;
const binMax = Math.max(...permDeltas, realDelta) + 0.005;
const nBins = 40;
const binWidth = (binMax - binMin) / nBins;
for (let i = 0; i < nBins; i++) {
  const lo = binMin + i * binWidth;
  const hi = lo + binWidth;
  const count = permDeltas.filter(d => d >= lo && d < hi).length;
  permHistBins.push({ binCenter: +((lo + hi) / 2).toFixed(4), count });
}

console.log(`\n  ─── Fisher z-test ───\n`);

const rObs = data.sparc.perGalaxy.r;
const rSim = data.simulation[0].perGalaxy.r;
const nObs = data.sparc.perGalaxy.n;
const nSim = data.simulation[0].perGalaxy.n;

const zObs = 0.5 * Math.log((1 + Math.abs(rObs)) / (1 - Math.abs(rObs)));
const zSim = 0.5 * Math.log((1 + Math.abs(rSim)) / (1 - Math.abs(rSim)));
// Note: Using |r| since both correlations are negative; we compare magnitude of linear association
const seDiff = Math.sqrt(1 / (nObs - 3) + 1 / (nSim - 3));
const zScore = (zObs - zSim) / seDiff;

function normalCDF(z) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
  const p = 0.3275911;
  const sign = z < 0 ? -1 : 1;
  z = Math.abs(z) / Math.sqrt(2);
  const t = 1.0 / (1.0 + p * z);
  const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-z * z);
  return 0.5 * (1.0 + sign * y);
}

const pFisher = 2 * (1 - normalCDF(Math.abs(zScore)));

console.log(`  |r_obs| = ${Math.abs(rObs).toFixed(4)}, z_obs = ${zObs.toFixed(4)}`);
console.log(`  |r_sim| = ${Math.abs(rSim).toFixed(4)}, z_sim = ${zSim.toFixed(4)}`);
console.log(`  z-score = ${zScore.toFixed(4)}, p = ${pFisher.toFixed(6)}`);
console.log(`  → Correlations are ${pFisher < 0.05 ? 'SIGNIFICANTLY DIFFERENT' : 'not significantly different'}`);

console.log(`\n  ─── Effect size (Cohen's d) ───\n`);

const pooledSD = Math.sqrt((sparcSD ** 2 + simSD ** 2) / 2);
const cohensD = Math.abs(sparcMean - simMean) / pooledSD;
console.log(`  Cohen's d = ${cohensD.toFixed(3)}`);
console.log(`  Effect size: ${cohensD < 0.2 ? 'negligible' : cohensD < 0.5 ? 'small' : cohensD < 0.8 ? 'medium' : 'LARGE'}`);

const slopeRatio = sparcMean / simMean;
const slopeRatioErr = Math.abs(slopeRatio) * Math.sqrt((sparcSD / sparcMean) ** 2 + (simSD / simMean) ** 2);
console.log(`  Slope ratio: ${slopeRatio.toFixed(3)} ± ${slopeRatioErr.toFixed(3)}`);

console.log(`\n  ─── Welch's t-test on bootstrap distributions ───\n`);

const tStat = (sparcMean - simMean) / Math.sqrt(sparcSD ** 2 / N_BOOT + simSD ** 2 / N_BOOT);
const df = ((sparcSD ** 2 / N_BOOT + simSD ** 2 / N_BOOT) ** 2) /
  ((sparcSD ** 2 / N_BOOT) ** 2 / (N_BOOT - 1) + (simSD ** 2 / N_BOOT) ** 2 / (N_BOOT - 1));
const pWelch = 2 * (1 - normalCDF(Math.abs(tStat)));
console.log(`  t = ${tStat.toFixed(2)}, df = ${df.toFixed(0)}, p < ${pWelch < 1e-10 ? '1e-10' : pWelch.toFixed(8)}`);

console.log(`\n${'═'.repeat(70)}`);
console.log(`  SUMMARY`);
console.log(`${'═'.repeat(70)}`);
console.log(`  Bootstrap CIs: SPARC [${sparcCI[0].toFixed(4)}, ${sparcCI[1].toFixed(4)}] vs ΛCDM [${simCI[0].toFixed(4)}, ${simCI[1].toFixed(4)}]`);
console.log(`  CIs overlap: ${ciOverlap ? 'YES — not clearly separated' : 'NO — clearly separated'}`);
console.log(`  Permutation p: ${pValue.toFixed(6)}`);
console.log(`  Fisher z p: ${pFisher.toFixed(6)}`);
console.log(`  Cohen's d: ${cohensD.toFixed(3)} (${cohensD < 0.2 ? 'negligible' : cohensD < 0.5 ? 'small' : cohensD < 0.8 ? 'medium' : 'LARGE'})`);
console.log(`  Slope ratio: ${slopeRatio.toFixed(3)} ± ${slopeRatioErr.toFixed(3)}`);
console.log(`  Welch's t: ${tStat.toFixed(2)}, p < ${pWelch < 1e-10 ? '1e-10' : pWelch.toFixed(8)}`);
console.log(`  VERDICT: ${pValue < 0.05 && !ciOverlap ? 'STATISTICALLY SIGNIFICANT — observed slope is steeper than ΛCDM' : pValue < 0.05 ? 'SIGNIFICANT by permutation, but CIs overlap — marginal' : 'NOT clearly significant'}`);

const sparcBootHist = [];
const simBootHist = [];
const allSlopes = [...sparcBootSlopes, ...simBootSlopes];
const slopeMin = Math.min(...allSlopes) - 0.005;
const slopeMax = Math.max(...allSlopes) + 0.005;
const slopeBinW = (slopeMax - slopeMin) / 50;
for (let i = 0; i < 50; i++) {
  const lo = slopeMin + i * slopeBinW;
  const hi = lo + slopeBinW;
  const center = +((lo + hi) / 2).toFixed(5);
  sparcBootHist.push({ binCenter: center, count: sparcBootSlopes.filter(s => s >= lo && s < hi).length });
  simBootHist.push({ binCenter: center, count: simBootSlopes.filter(s => s >= lo && s < hi).length });
}

const significanceTest = {
  nBootstrap: N_BOOT,
  nPermutations: N_PERM,
  sparc: {
    observedSlope: data.sparc.perGalaxy.slope,
    bootstrapMean: sparcMean,
    bootstrapSD: sparcSD,
    ci95: sparcCI,
    n: sparcGals.length,
  },
  lcdm: {
    observedSlope: data.simulation[0].perGalaxy.slope,
    bootstrapMean: simMean,
    bootstrapSD: simSD,
    ci95: simCI,
    n: simGals.length,
  },
  ciOverlap,
  permutation: {
    realDelta,
    pValue,
    moreExtreme,
    nullMean: permMean,
    nullSD: permSD,
    histogram: permHistBins,
  },
  fisherZ: {
    rObs,
    rSim,
    zObs,
    zSim,
    zScore,
    pValue: pFisher,
  },
  effectSize: {
    cohensD,
    label: cohensD < 0.2 ? 'negligible' : cohensD < 0.5 ? 'small' : cohensD < 0.8 ? 'medium' : 'large',
    slopeRatio,
    slopeRatioErr,
  },
  welchT: {
    tStat,
    df,
    pValue: pWelch,
  },
  bootstrapHistograms: {
    sparc: sparcBootHist,
    lcdm: simBootHist,
  },
  whyStronger: [
    {
      mechanism: 'Adiabatic contraction',
      description: 'Real baryonic infall contracts the DM halo, creating a tighter density–DM coupling than simple NFW profiles predict.',
      prediction: 'High-Σ galaxies should show steeper inner DM profiles than NFW.',
      testable: true,
    },
    {
      mechanism: 'Feedback-driven core creation',
      description: 'Stellar feedback in low-Σ dwarfs creates DM cores, boosting f_DM at low density more than NFW expects.',
      prediction: 'Dwarf galaxies should show DM cores (not cusps) proportional to star formation history.',
      testable: true,
    },
    {
      mechanism: 'Abundance matching scatter',
      description: 'Our mock uses a deterministic M★–M₂₀₀ relation. Real scatter in this relation could steepen the slope.',
      prediction: 'Adding realistic scatter (0.2 dex) to M★–M₂₀₀ should produce slopes closer to observed.',
      testable: true,
    },
    {
      mechanism: 'MOND-like density threshold',
      description: 'In MOND/modified gravity, the transition from Newtonian to modified regime depends on surface density (Σ ∝ a₀/G), naturally creating a steeper f_DM–Σ relation.',
      prediction: 'The slope should correlate with a₀ = 1.2×10⁻¹⁰ m/s² and show a characteristic break.',
      testable: true,
    },
  ],
};

data.significanceTest = significanceTest;
fs.writeFileSync(dataPath, JSON.stringify(data, null, 2));
console.log(`\nSaved significanceTest to fdm-analysis.json`);
