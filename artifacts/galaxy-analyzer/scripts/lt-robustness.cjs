const fs = require('fs');
const path = require('path');

const G_DAG_KPC = 3.7032;

function parseTable1(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const galaxies = {};
  for (const line of lines) {
    const name = line.substring(0, 8).trim();
    const dist = parseFloat(line.substring(32, 36).trim());
    const inc = parseFloat(line.substring(59, 63).trim());
    const eInc = parseFloat(line.substring(64, 68).trim());
    galaxies[name] = { name, dist, inc, eInc };
  }
  return galaxies;
}

function parseTable2(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const results = {};
  for (const line of lines) {
    const name = line.substring(0, 8).trim();
    const rmax = parseFloat(line.substring(9, 13).trim());
    const vRmax = parseFloat(line.substring(19, 24).trim());
    const mgasStr = line.substring(139, 145).trim();
    const mstarKStr = line.substring(146, 151).trim();
    const mstarSEDStr = line.substring(152, 157).trim();
    const mgas = mgasStr ? parseFloat(mgasStr) * 1e7 : 0;
    const mstarK = mstarKStr ? parseFloat(mstarKStr) * 1e7 : 0;
    const mstarSED = mstarSEDStr ? parseFloat(mstarSEDStr) * 1e7 : 0;
    const mstar = mstarSED > 0 ? mstarSED : mstarK;
    results[name] = { rmax, vRmax, mgas, mstar, mstarK, mstarSED };
  }
  return results;
}

function parseRotCurve(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const galaxies = {};
  for (const line of lines) {
    const name = line.substring(0, 8).trim();
    const type = line.substring(9, 14).trim();
    if (type !== 'Data') continue;
    const r03 = parseFloat(line.substring(15, 23).trim());
    const v03 = parseFloat(line.substring(24, 34).trim());
    const rScaled = parseFloat(line.substring(35, 44).trim());
    const vScaled = parseFloat(line.substring(45, 54).trim());
    const eVScaled = parseFloat(line.substring(55, 63).trim());
    if (!galaxies[name]) galaxies[name] = [];
    galaxies[name].push({ r: rScaled * r03, v: vScaled * v03, ev: eVScaled * v03 });
  }
  return galaxies;
}

function interpolateV(curve, r) {
  if (curve.length === 0) return NaN;
  if (r <= curve[0].r) return curve[0].v * (r / curve[0].r);
  if (r >= curve[curve.length - 1].r) return curve[curve.length - 1].v;
  for (let i = 0; i < curve.length - 1; i++) {
    if (r >= curve[i].r && r <= curve[i + 1].r) {
      const frac = (r - curve[i].r) / (curve[i + 1].r - curve[i].r);
      return curve[i].v + frac * (curve[i + 1].v - curve[i].v);
    }
  }
  return curve[curve.length - 1].v;
}

function mcGillRAR(gbar) {
  const y = gbar / G_DAG_KPC;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function linearRegression(x, y) {
  const n = x.length;
  if (n < 3) return { slope: NaN, intercept: NaN, r2: NaN, r: NaN, sSlope: NaN, n };
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) { sxx += (x[i] - mx) ** 2; sxy += (x[i] - mx) * (y[i] - my); }
  const slope = sxy / sxx;
  const intercept = my - slope * mx;
  const yPred = x.map(xi => intercept + slope * xi);
  let ssRes = 0, ssTot = 0;
  for (let i = 0; i < n; i++) { ssRes += (y[i] - yPred[i]) ** 2; ssTot += (y[i] - my) ** 2; }
  const r2 = ssTot > 0 ? 1 - ssRes / ssTot : 0;
  const r = Math.sign(slope) * Math.sqrt(Math.abs(r2));
  const sSlope = n > 2 ? Math.sqrt(ssRes / (n - 2) / sxx) : NaN;
  return { slope, intercept, r2, r, sSlope, n };
}

function pearsonR(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b) / n;
  const my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxy / Math.sqrt(sxx * syy);
}

function partialR(x, y, z) {
  const rxy = pearsonR(x, y);
  const rxz = pearsonR(x, z);
  const ryz = pearsonR(y, z);
  return (rxy - rxz * ryz) / Math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2));
}

function gaussRandom() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

const table1 = parseTable1('/tmp/little_things/table1.dat');
const table2 = parseTable2('/tmp/little_things/table2.dat');
const rotTotal = parseRotCurve('/tmp/little_things/rotdmbar.dat');
const rotDM = parseRotCurve('/tmp/little_things/rotdm.dat');

function runAnalysis(upsilonStar, incPerts, distPerts) {
  const deltas = [];
  const logSigs = [];
  const logVmaxes = [];

  for (const name of Object.keys(rotTotal)) {
    if (!rotDM[name] || !table2[name] || !table1[name]) continue;
    const t1 = table1[name];
    const t2 = table2[name];
    const curveTotal = rotTotal[name].sort((a, b) => a.r - b.r);
    const curveDM = rotDM[name].sort((a, b) => a.r - b.r);
    if (curveTotal.length < 5 || t2.mstar <= 0 && t2.mgas <= 0) continue;

    const incDeg = t1.inc + (incPerts?.[name] || 0);
    const distFactor = 1 + (distPerts?.[name] || 0);
    const incRad = incDeg * Math.PI / 180;
    const incOrig = t1.inc * Math.PI / 180;
    const sinRatio = Math.sin(incOrig) / Math.sin(incRad);

    const mbar = upsilonStar * t2.mstar + 1.33 * t2.mgas;
    if (mbar <= 0) continue;
    const rmax = t2.rmax * distFactor;

    const outerPoints = curveTotal.filter(p => {
      const rAdj = p.r * distFactor;
      return rAdj >= rmax * 0.4 && rAdj <= rmax;
    });
    if (outerPoints.length < 3) continue;

    const ptDeltas = [];
    const ptLogSigs = [];

    for (const pt of outerPoints) {
      const r = pt.r * distFactor;
      const vObs = pt.v * sinRatio;
      const vDMraw = interpolateV(curveDM, pt.r);
      if (isNaN(vDMraw) || vObs <= 0 || r <= 0) continue;
      const vDM = vDMraw * sinRatio;

      let vBarSq = vObs * vObs - vDM * vDM;
      if (vBarSq < 0) vBarSq = 0;

      const gObs = vObs * vObs / r;
      const gBar = vBarSq / r;
      if (gBar <= 0 || gObs <= 0) continue;

      const gRAR = mcGillRAR(gBar);
      const deltaRAR = Math.log10(gObs) - Math.log10(gRAR);
      const sigBar = mbar / (Math.PI * r * r);
      ptDeltas.push(deltaRAR);
      ptLogSigs.push(Math.log10(sigBar));
    }

    if (ptDeltas.length >= 3) {
      const meanD = ptDeltas.reduce((a, b) => a + b) / ptDeltas.length;
      const meanS = ptLogSigs.reduce((a, b) => a + b) / ptLogSigs.length;
      deltas.push(meanD);
      logSigs.push(meanS);
      logVmaxes.push(Math.log10(t2.vRmax));
    }
  }

  if (deltas.length < 5) return null;
  const reg = linearRegression(logSigs, deltas);
  const pr = partialR(logSigs, deltas, logVmaxes);
  return { slope: reg.slope, slopeErr: reg.sSlope, r: reg.r, r2: reg.r2, n: reg.n, partialR: pr };
}

console.log(`\n${'='.repeat(70)}`);
console.log(`  LITTLE THINGS — Υ SENSITIVITY + MONTE CARLO ROBUSTNESS`);
console.log(`${'='.repeat(70)}`);

console.log(`\n╔══════════════════════════════════════════╗`);
console.log(`║  TEST A: Υ_star SENSITIVITY (7 values)  ║`);
console.log(`╚══════════════════════════════════════════╝`);

const upsilonValues = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8];
const upsilonResults = [];

console.log(`\n  Υ_star   slope b     ±σ_b      r       partial_r  n   sign`);
console.log(`  ${'─'.repeat(65)}`);

for (const ups of upsilonValues) {
  const res = runAnalysis(ups, null, null);
  const sign = res.slope < 0 ? '✓ NEG' : '✗ POS';
  console.log(`  ${ups.toFixed(1)}      ${res.slope.toFixed(4).padStart(8)}  ${res.slopeErr.toFixed(4).padStart(8)}  ${res.r.toFixed(4).padStart(7)}  ${res.partialR.toFixed(4).padStart(9)}  ${String(res.n).padStart(2)}  ${sign}`);
  upsilonResults.push({ upsilon: ups, slope: res.slope, slopeErr: res.slopeErr, r: res.r, r2: res.r2, partialR: res.partialR, n: res.n, negative: res.slope < 0 });
}

const allNeg = upsilonResults.every(r => r.negative);
const negCount = upsilonResults.filter(r => r.negative).length;
console.log(`\n  Slope negative in ${negCount}/${upsilonResults.length} configurations`);
console.log(`  Υ sensitivity: ${allNeg ? 'PASS ✓ — slope negative for ALL Υ values' : `PARTIAL — negative in ${negCount}/${upsilonResults.length}`}`);

console.log(`\n╔══════════════════════════════════════════════════╗`);
console.log(`║  TEST B: MONTE CARLO (1000 iterations)          ║`);
console.log(`║  Perturbing inclination (±eInc) and distance    ║`);
console.log(`║  Distance uncertainty: ±10% (typical for dIrr)  ║`);
console.log(`╚══════════════════════════════════════════════════╝`);

const N_MC = 1000;
const mcSlopes = [];
const mcPartialRs = [];
const mcRs = [];

for (let iter = 0; iter < N_MC; iter++) {
  const incPerts = {};
  const distPerts = {};
  for (const name of Object.keys(table1)) {
    const eInc = table1[name].eInc || 3;
    incPerts[name] = gaussRandom() * eInc;
    distPerts[name] = gaussRandom() * 0.10;
  }
  const res = runAnalysis(0.5, incPerts, distPerts);
  if (res) {
    mcSlopes.push(res.slope);
    mcPartialRs.push(res.partialR);
    mcRs.push(res.r);
  }
}

mcSlopes.sort((a, b) => a - b);
mcPartialRs.sort((a, b) => a - b);

const meanSlope = mcSlopes.reduce((a, b) => a + b) / mcSlopes.length;
const stdSlope = Math.sqrt(mcSlopes.reduce((a, b) => a + (b - meanSlope) ** 2, 0) / mcSlopes.length);
const ci025 = mcSlopes[Math.floor(mcSlopes.length * 0.025)];
const ci975 = mcSlopes[Math.floor(mcSlopes.length * 0.975)];
const fracNeg = mcSlopes.filter(s => s < 0).length / mcSlopes.length;

const meanPR = mcPartialRs.reduce((a, b) => a + b) / mcPartialRs.length;
const stdPR = Math.sqrt(mcPartialRs.reduce((a, b) => a + (b - meanPR) ** 2, 0) / mcPartialRs.length);
const prCI025 = mcPartialRs[Math.floor(mcPartialRs.length * 0.025)];
const prCI975 = mcPartialRs[Math.floor(mcPartialRs.length * 0.975)];
const fracPRNeg = mcPartialRs.filter(r => r < 0).length / mcPartialRs.length;

console.log(`\n  Completed ${mcSlopes.length}/${N_MC} valid iterations`);
console.log(`\n  --- Slope b Distribution ---`);
console.log(`  Mean:   ${meanSlope.toFixed(4)}`);
console.log(`  Std:    ${stdSlope.toFixed(4)}`);
console.log(`  95% CI: [${ci025.toFixed(4)}, ${ci975.toFixed(4)}]`);
console.log(`  Fraction negative: ${(fracNeg * 100).toFixed(1)}% (${mcSlopes.filter(s => s < 0).length}/${mcSlopes.length})`);
console.log(`  CI entirely below zero: ${ci975 < 0 ? 'YES ✓' : 'NO'}`);

console.log(`\n  --- Partial r Distribution ---`);
console.log(`  Mean:   ${meanPR.toFixed(4)}`);
console.log(`  Std:    ${stdPR.toFixed(4)}`);
console.log(`  95% CI: [${prCI025.toFixed(4)}, ${prCI975.toFixed(4)}]`);
console.log(`  Fraction negative: ${(fracPRNeg * 100).toFixed(1)}%`);

const mcPass = fracNeg >= 0.90;
console.log(`\n  Monte Carlo: ${mcPass ? 'PASS ✓' : 'PARTIAL'} — ${(fracNeg * 100).toFixed(1)}% of slopes negative`);

console.log(`\n╔══════════════════════════════════════════════════╗`);
console.log(`║  TEST C: COMBINED Υ × MC (7 × 200 = 1400 runs) ║`);
console.log(`╚══════════════════════════════════════════════════╝`);

const combinedResults = [];
const N_MC_COMB = 200;

for (const ups of upsilonValues) {
  let negCount = 0;
  const slopes = [];
  for (let iter = 0; iter < N_MC_COMB; iter++) {
    const incPerts = {};
    const distPerts = {};
    for (const name of Object.keys(table1)) {
      const eInc = table1[name].eInc || 3;
      incPerts[name] = gaussRandom() * eInc;
      distPerts[name] = gaussRandom() * 0.10;
    }
    const res = runAnalysis(ups, incPerts, distPerts);
    if (res) {
      slopes.push(res.slope);
      if (res.slope < 0) negCount++;
    }
  }
  const frac = slopes.length > 0 ? negCount / slopes.length : 0;
  const mean = slopes.reduce((a, b) => a + b, 0) / slopes.length;
  combinedResults.push({ upsilon: ups, fracNeg: frac, meanSlope: mean, nValid: slopes.length });
  console.log(`  Υ=${ups.toFixed(1)}: ${(frac * 100).toFixed(1)}% negative (mean b=${mean.toFixed(4)}, n=${slopes.length})`);
}

const allCombNeg = combinedResults.every(r => r.fracNeg >= 0.50);
console.log(`\n  All Υ have >50% negative: ${allCombNeg ? 'YES ✓' : 'NO ✗'}`);

console.log(`\n${'='.repeat(70)}`);
console.log(`  OVERALL ROBUSTNESS VERDICT`);
console.log(`${'='.repeat(70)}`);
console.log(`  A. Υ sensitivity:    ${allNeg ? 'PASS' : 'PARTIAL'} (${negCount}/${upsilonResults.length} negative)`);
console.log(`  B. Monte Carlo:      ${mcPass ? 'PASS' : 'PARTIAL'} (${(fracNeg * 100).toFixed(1)}% negative)`);
console.log(`  C. Combined Υ × MC: ${allCombNeg ? 'PASS' : 'PARTIAL'}`);
const overallPass = allNeg && mcPass;
console.log(`\n  >>> ${overallPass ? 'ROBUST: Result survives Υ variation AND measurement uncertainty' : 'FRAGILE: Result sensitive to systematics — interpret with caution'} <<<`);

const output = {
  upsilonSensitivity: {
    values: upsilonResults,
    allNegative: allNeg,
    negCount,
    total: upsilonResults.length
  },
  monteCarlo: {
    nIterations: mcSlopes.length,
    slope: { mean: meanSlope, std: stdSlope, ci95: [ci025, ci975], fracNegative: fracNeg },
    partialR: { mean: meanPR, std: stdPR, ci95: [prCI025, prCI975], fracNegative: fracPRNeg },
    allSlopes: mcSlopes,
    allPartialRs: mcPartialRs
  },
  combined: combinedResults,
  verdict: {
    upsilonPass: allNeg,
    mcPass,
    combinedPass: allCombNeg,
    overallRobust: overallPass
  }
};

const outPath = path.join(__dirname, '..', 'public', 'little-things-replication.json');
const existing = JSON.parse(fs.readFileSync(outPath, 'utf8'));
existing.robustness = output;
fs.writeFileSync(outPath, JSON.stringify(existing, null, 2));
console.log(`\nResults merged into: ${outPath}`);
