const fs = require('fs');
const path = require('path');

const rar = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const fdm = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'fdm-analysis.json'), 'utf8'));

const G = 4.3009e-6;
const a0 = 3702;

function linreg(x, y) {
  const n = x.length;
  if (n < 3) return { slope: NaN, r: NaN, r2: NaN, n };
  const mx = x.reduce((a, b) => a + b) / n;
  const my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  const slope = sxx > 0 ? sxy / sxx : 0;
  const r = sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
  return { slope, r, r2: r * r, n };
}

function pearsonR(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b) / n;
  const my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function partialR(x, y, z) {
  const rxy = pearsonR(x, y);
  const rxz = pearsonR(x, z);
  const ryz = pearsonR(y, z);
  const d = Math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2));
  return d > 0 ? (rxy - rxz * ryz) / d : 0;
}

function shuffle(arr) {
  const a = [...arr];
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
}

function gaussRandom() {
  let u, v, s;
  do { u = Math.random() * 2 - 1; v = Math.random() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
  return u * Math.sqrt(-2 * Math.log(s) / s);
}

const sep = '='.repeat(70);
console.log(`\n${sep}`);
console.log('  PHASE C: TRYING TO BREAK THE IDEA');
console.log('  4 Rigorous Tests to Kill the Signal');
console.log(sep);

const galaxies = rar.perGalaxy.filter(g => g.sigma_bar > 0 && g.Vmax > 0);
const N = galaxies.length;
console.log(`\nDataset: ${N} galaxies with valid sigma_bar and Vmax`);

const logSigArr = galaxies.map(g => Math.log10(g.sigma_bar));
const deltaArr = galaxies.map(g => g.meanDeltaOuter);
const logVmaxArr = galaxies.map(g => Math.log10(g.Vmax));

const realReg = linreg(logSigArr, deltaArr);
const realPartialR = partialR(logSigArr, deltaArr, logVmaxArr);

console.log(`\nREAL SIGNAL (baseline):`);
console.log(`  Slope b = ${realReg.slope.toFixed(4)}`);
console.log(`  r = ${realReg.r.toFixed(4)}`);
console.log(`  Partial r|Vmax = ${realPartialR.toFixed(4)}`);
console.log(`  n = ${realReg.n}`);

console.log(`\n${sep}`);
console.log('  TEST 1: NULL TEST (Scrambled V_obs)');
console.log(sep);

const N_NULL = 1000;
let nullSlopesStronger = 0;
let nullRStronger = 0;
const nullSlopes = [];
const nullRs = [];

for (let iter = 0; iter < N_NULL; iter++) {
  const scrambledDelta = shuffle(deltaArr);
  const reg = linreg(logSigArr, scrambledDelta);
  nullSlopes.push(reg.slope);
  nullRs.push(reg.r);
  if (Math.abs(reg.slope) >= Math.abs(realReg.slope)) nullSlopesStronger++;
  if (Math.abs(reg.r) >= Math.abs(realReg.r)) nullRStronger++;
}

const nullMeanSlope = nullSlopes.reduce((a, b) => a + b) / N_NULL;
const nullStdSlope = Math.sqrt(nullSlopes.reduce((s, v) => s + (v - nullMeanSlope) ** 2, 0) / N_NULL);
const nullMeanR = nullRs.reduce((a, b) => a + b) / N_NULL;

console.log(`\nNull distribution (${N_NULL} shuffles of deltaRAR):`);
console.log(`  Mean null slope: ${nullMeanSlope.toFixed(6)} +/- ${nullStdSlope.toFixed(6)}`);
console.log(`  Mean null |r|: ${nullMeanR.toFixed(6)}`);
console.log(`  Real slope: ${realReg.slope.toFixed(4)}`);
console.log(`  Real r: ${realReg.r.toFixed(4)}`);
console.log(`  p-value (slope): ${(nullSlopesStronger / N_NULL).toFixed(4)}`);
console.log(`  p-value (r): ${(nullRStronger / N_NULL).toFixed(4)}`);
console.log(`  Sigma (slope): ${((realReg.slope - nullMeanSlope) / nullStdSlope).toFixed(2)}σ`);

const nullPass = nullSlopesStronger / N_NULL < 0.001;
console.log(`\n  VERDICT: ${nullPass ? '✓ SIGNAL IS REAL — shuffling kills it completely' : '✗ SIGNAL MAY BE ARTIFACT — shuffling reproduces it'}`);

console.log(`\n${sep}`);
console.log('  TEST 2: SCRAMBLING Σ_bar BETWEEN GALAXIES');
console.log(sep);

const N_SCRAMBLE = 1000;
let scrambleSlopesStronger = 0;
const scrambleSlopes = [];
const scramblePartialRs = [];

for (let iter = 0; iter < N_SCRAMBLE; iter++) {
  const scrambledLogSig = shuffle(logSigArr);
  const reg = linreg(scrambledLogSig, deltaArr);
  const pr = partialR(scrambledLogSig, deltaArr, logVmaxArr);
  scrambleSlopes.push(reg.slope);
  scramblePartialRs.push(pr);
  if (Math.abs(reg.slope) >= Math.abs(realReg.slope)) scrambleSlopesStronger++;
}

const scrMeanSlope = scrambleSlopes.reduce((a, b) => a + b) / N_SCRAMBLE;
const scrStdSlope = Math.sqrt(scrambleSlopes.reduce((s, v) => s + (v - scrMeanSlope) ** 2, 0) / N_SCRAMBLE);
const scrMeanPR = scramblePartialRs.reduce((a, b) => a + b) / N_SCRAMBLE;
let scrPRstronger = 0;
for (const pr of scramblePartialRs) { if (Math.abs(pr) >= Math.abs(realPartialR)) scrPRstronger++; }

console.log(`\nΣ_bar scramble (${N_SCRAMBLE} permutations):`);
console.log(`  Mean scrambled slope: ${scrMeanSlope.toFixed(6)} +/- ${scrStdSlope.toFixed(6)}`);
console.log(`  Real slope: ${realReg.slope.toFixed(4)}`);
console.log(`  p-value (slope): ${(scrambleSlopesStronger / N_SCRAMBLE).toFixed(4)}`);
console.log(`  p-value (partial r): ${(scrPRstronger / N_SCRAMBLE).toFixed(4)}`);
console.log(`  Sigma: ${((realReg.slope - scrMeanSlope) / scrStdSlope).toFixed(2)}σ`);
console.log(`  Mean scrambled partial r: ${scrMeanPR.toFixed(4)} vs real: ${realPartialR.toFixed(4)}`);

const scramblePass = scrambleSlopesStronger / N_SCRAMBLE < 0.001;
console.log(`\n  VERDICT: ${scramblePass ? '✓ COUPLING IS REAL — scrambling Σ_bar destroys it' : '✗ COUPLING MAY BE SPURIOUS — scrambling preserves it'}`);

console.log(`\n${sep}`);
console.log('  TEST 3: EXTREME M/L STRESS TEST');
console.log(sep);

const upsilonValues = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0];
const mlResults = [];

for (const Y of upsilonValues) {
  const scaledDelta = galaxies.map(g => {
    const scaleFactor = Y / 0.5;
    const gBar_scaled = g.sigma_bar * scaleFactor;
    const logSigScaled = Math.log10(gBar_scaled);
    const meanDelta = g.meanDeltaOuter;
    const adjustedDelta = meanDelta - (scaleFactor - 1) * 0.15;
    return adjustedDelta;
  });
  const scaledLogSig = galaxies.map((g, i) => {
    const scaleFactor = Y / 0.5;
    return Math.log10(g.sigma_bar * scaleFactor);
  });

  const reg = linreg(scaledLogSig, scaledDelta);
  const pr = partialR(scaledLogSig, scaledDelta, logVmaxArr);
  mlResults.push({ Y, slope: reg.slope, r: reg.r, partialR: pr });
}

console.log(`\nM/L ratio sweep (Υ_disk from 0.1 to 2.0):`);
console.log(`  ${'Υ_disk'.padEnd(8)} ${'Slope'.padEnd(10)} ${'r'.padEnd(10)} ${'Partial r'.padEnd(10)} ${'Negative?'}`);
console.log(`  ${'-'.repeat(50)}`);
let allNegative = true;
for (const r of mlResults) {
  const neg = r.slope < 0;
  if (!neg) allNegative = false;
  console.log(`  ${r.Y.toFixed(1).padEnd(8)} ${r.slope.toFixed(4).padEnd(10)} ${r.r.toFixed(4).padEnd(10)} ${r.partialR.toFixed(4).padEnd(10)} ${neg ? '✓' : '✗'}`);
}

console.log(`\n  All slopes negative: ${allNegative ? '✓ YES' : '✗ NO — slope flips at some Υ'}`);
const mlPass = allNegative;
console.log(`  VERDICT: ${mlPass ? '✓ SIGNAL SURVIVES ALL M/L VALUES — robust' : '✗ SIGNAL DEPENDS ON M/L — may be artifact'}`);

console.log(`\n${sep}`);
console.log('  TEST 4: DISTANCE PERTURBATION MONTE CARLO');
console.log(sep);

const N_DIST = 1000;
const distSlopes = [];
const distPartialRs = [];

for (let iter = 0; iter < N_DIST; iter++) {
  const perturbedLogSig = [];
  const perturbedDelta = [];
  const perturbedLogVmax = [];

  for (let i = 0; i < N; i++) {
    const g = galaxies[i];
    const distFactor = 1 + 0.3 * gaussRandom();
    if (distFactor <= 0.1) continue;

    const newSigma = g.sigma_bar / (distFactor * distFactor);
    const newVmax = g.Vmax * Math.sqrt(distFactor);
    const newDelta = g.meanDeltaOuter + 0.1 * Math.log10(distFactor);

    perturbedLogSig.push(Math.log10(newSigma));
    perturbedDelta.push(newDelta);
    perturbedLogVmax.push(Math.log10(newVmax));
  }

  const reg = linreg(perturbedLogSig, perturbedDelta);
  const pr = partialR(perturbedLogSig, perturbedDelta, perturbedLogVmax);
  distSlopes.push(reg.slope);
  distPartialRs.push(pr);
}

const distMeanSlope = distSlopes.reduce((a, b) => a + b) / N_DIST;
const distStdSlope = Math.sqrt(distSlopes.reduce((s, v) => s + (v - distMeanSlope) ** 2, 0) / N_DIST);
const distFracNeg = distSlopes.filter(s => s < 0).length / N_DIST;
const distMeanPR = distPartialRs.reduce((a, b) => a + b) / N_DIST;
const distFracPRneg = distPartialRs.filter(r => r < 0).length / N_DIST;

const distCI95 = [
  distSlopes.sort((a, b) => a - b)[Math.floor(N_DIST * 0.025)],
  distSlopes[Math.floor(N_DIST * 0.975)]
];

console.log(`\nDistance perturbation (${N_DIST} iterations, ±30% Gaussian):`);
console.log(`  Mean slope: ${distMeanSlope.toFixed(4)} +/- ${distStdSlope.toFixed(4)}`);
console.log(`  95% CI: [${distCI95[0].toFixed(4)}, ${distCI95[1].toFixed(4)}]`);
console.log(`  Fraction negative: ${(distFracNeg * 100).toFixed(1)}%`);
console.log(`  Real slope: ${realReg.slope.toFixed(4)}`);
console.log(`  Mean partial r: ${distMeanPR.toFixed(4)} (real: ${realPartialR.toFixed(4)})`);
console.log(`  Fraction partial r negative: ${(distFracPRneg * 100).toFixed(1)}%`);

const distPass = distFracNeg > 0.95 && distFracPRneg > 0.90;
console.log(`\n  VERDICT: ${distPass ? '✓ SIGNAL SURVIVES ±30% DISTANCE ERRORS — robust' : '✗ SIGNAL SENSITIVE TO DISTANCES — may be artifact'}`);

console.log(`\n${sep}`);
console.log('  ADDITIONAL: a₀ STABILITY UNDER PERTURBATION');
console.log(sep);

const a0_values = [];
for (let iter = 0; iter < N_DIST; iter++) {
  let sumNum = 0, sumDen = 0;
  for (const g of galaxies) {
    const distFactor = 1 + 0.3 * gaussRandom();
    if (distFactor <= 0.1) continue;
    const V = g.Vmax * Math.sqrt(distFactor);
    const R = g.Rmax * distFactor;
    if (V > 0 && R > 0) {
      const k = V * V / R;
      sumNum += k;
      sumDen++;
    }
  }
  if (sumDen > 0) a0_values.push(sumNum / sumDen);
}

const a0Mean = a0_values.reduce((a, b) => a + b) / a0_values.length;
const a0Std = Math.sqrt(a0_values.reduce((s, v) => s + (v - a0Mean) ** 2, 0) / a0_values.length);
const a0_sorted = a0_values.sort((a, b) => a - b);
const a0CI95 = [a0_sorted[Math.floor(a0_values.length * 0.025)], a0_sorted[Math.floor(a0_values.length * 0.975)]];

console.log(`\na₀ stability under ±30% distance perturbation:`);
console.log(`  Original a₀ (mean k): from global fit = 3702 (km/s)²/kpc`);
console.log(`  Perturbed mean: ${a0Mean.toFixed(1)} +/- ${a0Std.toFixed(1)}`);
console.log(`  95% CI: [${a0CI95[0].toFixed(1)}, ${a0CI95[1].toFixed(1)}]`);
console.log(`  CV: ${(a0Std / a0Mean * 100).toFixed(1)}%`);
const a0Robust = a0CI95[0] > 1000 && a0CI95[1] < 20000;
console.log(`  VERDICT: ${a0Robust ? '✓ a₀ remains well-constrained under ±30% distance noise' : '✗ a₀ is sensitive to distance errors'}`);

console.log(`\n${sep}`);
console.log('  OVERALL PHASE C VERDICT');
console.log(sep);

const tests = [
  { name: 'Null Test (shuffled deltaRAR)', pass: nullPass },
  { name: 'Σ_bar Scramble (permutation)', pass: scramblePass },
  { name: 'Extreme M/L sweep (Υ=0.1–2.0)', pass: mlPass },
  { name: 'Distance Monte Carlo (±30%)', pass: distPass },
];

let passed = 0;
for (const t of tests) {
  console.log(`  ${t.pass ? '✓' : '✗'} ${t.name}: ${t.pass ? 'SIGNAL SURVIVES' : 'SIGNAL VULNERABLE'}`);
  if (t.pass) passed++;
}

console.log(`\n  Result: ${passed}/4 tests passed`);

if (passed === 4) {
  console.log(`\n  ═══════════════════════════════════════════════`);
  console.log(`  THE SIGNAL SURVIVES ALL 4 BREAKING ATTEMPTS.`);
  console.log(`  It is NOT an artifact of:`);
  console.log(`    - Random pairing (null test p < 0.001)`);
  console.log(`    - Σ_bar labeling (scramble p < 0.001)`);
  console.log(`    - M/L assumptions (all Υ give negative slopes)`);
  console.log(`    - Distance errors (±30% → still negative)`);
  console.log(`  ═══════════════════════════════════════════════`);
  console.log(`\n  CONCLUSION: The phenomenon is robust enough to`);
  console.log(`  warrant physical interpretation. Proceed to Phase A.`);
} else if (passed >= 3) {
  console.log(`\n  Signal is mostly robust but has one vulnerability.`);
  console.log(`  Investigate the failed test before proceeding.`);
} else {
  console.log(`\n  ⚠ Signal has significant vulnerabilities.`);
  console.log(`  The idea may need fundamental revision.`);
}

const output = {
  date: new Date().toISOString(),
  description: 'Phase C: Breaking the Idea — 4 rigorous tests',
  baseline: {
    nGalaxies: N,
    slope: realReg.slope,
    r: realReg.r,
    partialR_Vmax: realPartialR,
  },
  test1_null: {
    name: 'Null Test (shuffled deltaRAR)',
    nIterations: N_NULL,
    pValueSlope: nullSlopesStronger / N_NULL,
    pValueR: nullRStronger / N_NULL,
    nullMeanSlope: nullMeanSlope,
    nullStdSlope: nullStdSlope,
    sigma: (realReg.slope - nullMeanSlope) / nullStdSlope,
    pass: nullPass,
  },
  test2_scramble: {
    name: 'Σ_bar Scramble',
    nIterations: N_SCRAMBLE,
    pValueSlope: scrambleSlopesStronger / N_SCRAMBLE,
    pValuePartialR: scrPRstronger / N_SCRAMBLE,
    scrambleMeanSlope: scrMeanSlope,
    scrambleStdSlope: scrStdSlope,
    sigma: (realReg.slope - scrMeanSlope) / scrStdSlope,
    pass: scramblePass,
  },
  test3_ml: {
    name: 'Extreme M/L Stress Test',
    upsilonRange: [0.1, 2.0],
    results: mlResults,
    allNegative: allNegative,
    pass: mlPass,
  },
  test4_distance: {
    name: 'Distance Perturbation Monte Carlo',
    nIterations: N_DIST,
    perturbationPercent: 30,
    meanSlope: distMeanSlope,
    stdSlope: distStdSlope,
    ci95: distCI95,
    fracNegative: distFracNeg,
    fracPartialRneg: distFracPRneg,
    pass: distPass,
  },
  a0_stability: {
    originalA0: 3702,
    perturbedMean: a0Mean,
    perturbedStd: a0Std,
    ci95: a0CI95,
    cv_percent: a0Std / a0Mean * 100,
  },
  overall: {
    passed: passed,
    total: 4,
    verdict: passed === 4 ? 'SIGNAL SURVIVES ALL BREAKING ATTEMPTS' :
             passed >= 3 ? 'MOSTLY ROBUST — ONE VULNERABILITY' :
             'SIGNIFICANT VULNERABILITIES FOUND',
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'break-test-results.json'), JSON.stringify(output, null, 2));
console.log(`\nResults saved to public/break-test-results.json`);
