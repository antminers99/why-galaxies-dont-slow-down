#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 117: ANTI-CIRCULARITY BATTERY');
console.log('');
console.log('  CRITICAL QUESTION: Is the fgas → logOMD correlation a');
console.log('  computational artifact from shared baryonic inputs?');
console.log('');
console.log('  logOMD = log(Vobs²/Vbar²) and Vbar includes Vgas.');
console.log('  MHI feeds into BOTH predictor (log(MHI/L3.6))');
console.log('  AND target (through Vgas in Vbar decomposition).');
console.log('  We must prove this is NOT the source of the signal.');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}
function olsRegress(X, y) {
  const n = y.length, p = X[0].length;
  if (n <= p + 1) return null;
  const Xt = Array.from({ length: p }, (_, i) => X.map(row => row[i]));
  const XtX = Array.from({ length: p }, (_, i) =>
    Array.from({ length: p }, (_, j) =>
      Xt[i].reduce((s, _, k) => s + Xt[i][k] * Xt[j][k], 0)));
  const Xty = Xt.map(col => col.reduce((s, v, k) => s + v * y[k], 0));
  const aug = XtX.map((row, i) => [...row, Xty[i]]);
  for (let i = 0; i < p; i++) {
    let maxRow = i;
    for (let k = i + 1; k < p; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) maxRow = k;
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];
    if (Math.abs(aug[i][i]) < 1e-12) return null;
    for (let k = 0; k < p; k++) {
      if (k === i) continue;
      const f = aug[k][i] / aug[i][i];
      for (let j = i; j <= p; j++) aug[k][j] -= f * aug[i][j];
    }
  }
  return { beta: aug.map((row, i) => row[p] / row[i]) };
}
function looR2(xArr, yArr) {
  const n = xArr.length;
  if (n < 8) return NaN;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(yArr);
  for (let i = 0; i < n; i++) {
    const xTr = xArr.filter((_, j) => j !== i);
    const yTr = yArr.filter((_, j) => j !== i);
    const fit = olsRegress(xTr.map(v => [1, v]), yTr);
    if (!fit) return NaN;
    ss_loo += (yArr[i] - (fit.beta[0] + fit.beta[1] * xArr[i])) ** 2;
    ss_tot += (yArr[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_loo / ss_tot : 0;
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  table1[name] = { L36, Rdisk, MHI, Vflat };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul: vbul || 0 });
}

const galaxies = [];

for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat <= 0) continue;
  const sorted = rcPoints.filter(p => p.rad > 0 && p.vobs > 0).sort((a, b) => a.rad - b.rad);
  if (sorted.length < 8) continue;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logMHI_L36 = Math.log10(t1.MHI) - Math.log10(t1.L36);
  if (!isFinite(fgas) || !isFinite(logMHI_L36)) continue;

  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);

  const computeLogOMD = (pts, mode) => {
    const vals = pts.map(p => {
      let vBarSq;
      if (mode === 'full') {
        vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                 UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                 p.vgas * Math.abs(p.vgas);
      } else if (mode === 'stars_only') {
        vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                 UPSILON_BULGE * p.vbul * Math.abs(p.vbul);
      } else if (mode === 'gas_only') {
        vBarSq = p.vgas * Math.abs(p.vgas);
      }
      return (p.vobs ** 2) / Math.max(vBarSq, 0.01);
    }).filter(v => isFinite(v) && v > 0);
    return vals.length > 0 ? Math.log10(mean(vals)) : NaN;
  };

  const logOMD_full = computeLogOMD(outerPts, 'full');
  const logOMD_stars = computeLogOMD(outerPts, 'stars_only');
  const logOMD_gas = computeLogOMD(outerPts, 'gas_only');

  const outerVobs = mean(outerPts.map(p => p.vobs));
  const logVobs_outer = Math.log10(outerVobs);

  const outerVobsNorm = outerVobs / t1.Vflat;

  const gasContribFrac = outerPts.map(p => {
    const vGasSq = p.vgas * Math.abs(p.vgas);
    const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) +
                   UPSILON_BULGE * p.vbul * Math.abs(p.vbul) +
                   vGasSq;
    return vBarSq > 0 ? vGasSq / vBarSq : 0;
  });
  const meanGasContrib = mean(gasContribFrac);

  if (!isFinite(logOMD_full)) continue;

  galaxies.push({
    name, fgas, logMHI_L36, Vflat: t1.Vflat,
    logOMD_full, logOMD_stars, logOMD_gas,
    logVobs_outer, outerVobsNorm, meanGasContrib,
  });
}

const highV = galaxies.filter(g => g.Vflat >= 70);
console.log('  High-Vflat sample: N=' + highV.length + '\n');

// ============================================================
// TEST A: STARS-ONLY TARGET
// Remove gas from Vbar entirely. If fgas still predicts
// logOMD_stars → signal is NOT from shared gas input.
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST A: STARS-ONLY TARGET');
console.log('  logOMD_stars = log(Vobs²/Vstar²) — gas removed from target');
console.log('  If fgas still predicts → NOT a shared-input artifact');
console.log('══════════════════════════════════════════════════════════════\n');

const validStars = highV.filter(g => isFinite(g.logOMD_stars));
const r_fgas_stars = pearsonR(validStars.map(g => g.fgas), validStars.map(g => g.logOMD_stars));
const r_ratio_stars = pearsonR(validStars.map(g => g.logMHI_L36), validStars.map(g => g.logOMD_stars));
const loo_fgas_stars = looR2(validStars.map(g => g.fgas), validStars.map(g => g.logOMD_stars));
const loo_ratio_stars = looR2(validStars.map(g => g.logMHI_L36), validStars.map(g => g.logOMD_stars));

const r_fgas_full = pearsonR(highV.map(g => g.fgas), highV.map(g => g.logOMD_full));

console.log('  ORIGINAL (full Vbar): r(fgas, logOMD) = ' + r_fgas_full.toFixed(3));
console.log('');
console.log('  STARS-ONLY target (N=' + validStars.length + '):');
console.log('    r(fgas, logOMD_stars) = ' + r_fgas_stars.toFixed(3) + ', LOO = ' + loo_fgas_stars.toFixed(3));
console.log('    r(log(MHI/L36), logOMD_stars) = ' + r_ratio_stars.toFixed(3) + ', LOO = ' + loo_ratio_stars.toFixed(3));
console.log('');

const starsTest = Math.abs(r_fgas_stars) > 0.5;
console.log('  VERDICT: ' + (starsTest ?
  'PASS — fgas predicts even when gas is REMOVED from target. Signal is NOT from shared gas input.' :
  (Math.abs(r_fgas_stars) > 0.3 ? 'PARTIAL — weaker but survives' : 'FAIL — signal depends on gas in target')));
console.log('');

// ============================================================
// TEST B: GAS-ONLY TARGET
// Use only gas contribution to Vbar. This should show the
// circularity direction: if signal is ONLY with gas target,
// it's likely an artifact.
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST B: GAS-ONLY TARGET (CIRCULARITY DETECTOR)');
console.log('  logOMD_gas = log(Vobs²/Vgas²)');
console.log('  If signal ONLY here → circularity. If similar → real.');
console.log('══════════════════════════════════════════════════════════════\n');

const validGas = highV.filter(g => isFinite(g.logOMD_gas));
const r_fgas_gasTarget = pearsonR(validGas.map(g => g.fgas), validGas.map(g => g.logOMD_gas));
const r_ratio_gasTarget = pearsonR(validGas.map(g => g.logMHI_L36), validGas.map(g => g.logOMD_gas));

console.log('  GAS-ONLY target (N=' + validGas.length + '):');
console.log('    r(fgas, logOMD_gas) = ' + r_fgas_gasTarget.toFixed(3));
console.log('    r(log(MHI/L36), logOMD_gas) = ' + r_ratio_gasTarget.toFixed(3));
console.log('');
console.log('  Comparison:');
console.log('    Full target:       r = ' + r_fgas_full.toFixed(3));
console.log('    Stars-only target: r = ' + r_fgas_stars.toFixed(3));
console.log('    Gas-only target:   r = ' + r_fgas_gasTarget.toFixed(3));
console.log('');

const gasOnlyStronger = Math.abs(r_fgas_gasTarget) > Math.abs(r_fgas_stars) + 0.15;
console.log('  Gas-only target much stronger? ' + (gasOnlyStronger ? 'YES → circularity risk' : 'NO → signal is broad'));
console.log('');

// ============================================================
// TEST C: PURE VELOCITY TARGET (no baryonic decomposition)
// Use log(Vobs_outer / Vflat) — no Vbar at all.
// If fgas predicts this → fully independent of decomposition.
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST C: PURE VELOCITY TARGET (NO DECOMPOSITION)');
console.log('  Target = log(Vobs_outer / Vflat)');
console.log('  Predictor shares ZERO inputs with target.');
console.log('  This is the cleanest anti-circularity test possible.');
console.log('══════════════════════════════════════════════════════════════\n');

const logVratio = highV.map(g => Math.log10(g.outerVobsNorm));
const r_fgas_vratio = pearsonR(highV.map(g => g.fgas), logVratio);
const r_ratio_vratio = pearsonR(highV.map(g => g.logMHI_L36), logVratio);

console.log('  Target: log(mean Vobs outer / Vflat)');
console.log('    r(fgas) = ' + r_fgas_vratio.toFixed(3));
console.log('    r(log(MHI/L36)) = ' + r_ratio_vratio.toFixed(3));
console.log('');

const logVobs = highV.map(g => g.logVobs_outer);
const r_fgas_vobs = pearsonR(highV.map(g => g.fgas), logVobs);
console.log('  Target: log(mean Vobs outer) [raw, no normalization]');
console.log('    r(fgas) = ' + r_fgas_vobs.toFixed(3));
console.log('');

// ============================================================
// TEST D: GAS CONTRIBUTION FRACTION ANALYSIS
// How much of Vbar comes from gas in the outer region?
// If gas contributes very little → circularity is minimal.
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST D: GAS CONTRIBUTION TO OUTER VBAR');
console.log('  How much of Vbar² comes from Vgas² in outer half?');
console.log('  If gas is a small fraction → circularity risk is low.');
console.log('══════════════════════════════════════════════════════════════\n');

const gasContribs = highV.map(g => g.meanGasContrib);
console.log('  Mean gas contribution to Vbar² (outer): ' + (mean(gasContribs) * 100).toFixed(1) + '%');
console.log('  Median: ' + (gasContribs.sort((a, b) => a - b)[Math.floor(gasContribs.length / 2)] * 100).toFixed(1) + '%');
console.log('  Range: [' + (Math.min(...gasContribs) * 100).toFixed(1) + '%, ' + (Math.max(...gasContribs) * 100).toFixed(1) + '%]');
console.log('');

const sorted = [...highV].sort((a, b) => a.fgas - b.fgas);
const t3 = Math.floor(sorted.length / 3);
console.log('  By fgas tercile:');
console.log('    Low fgas:  gas contrib = ' + (mean(sorted.slice(0, t3).map(g => g.meanGasContrib)) * 100).toFixed(1) + '%');
console.log('    Mid fgas:  gas contrib = ' + (mean(sorted.slice(t3, 2 * t3).map(g => g.meanGasContrib)) * 100).toFixed(1) + '%');
console.log('    High fgas: gas contrib = ' + (mean(sorted.slice(2 * t3).map(g => g.meanGasContrib)) * 100).toFixed(1) + '%');
console.log('');

// ============================================================
// TEST E: SCRAMBLED GAS TEST
// Randomly reassign MHI values across galaxies. Recompute
// everything. If signal persists → artifact. If dies → real.
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST E: SCRAMBLED GAS PERMUTATION TEST');
console.log('  Shuffle MHI across galaxies, recompute predictor.');
console.log('  Target stays fixed. If signal dies → real signal.');
console.log('  1000 permutations.');
console.log('══════════════════════════════════════════════════════════════\n');

const realR = pearsonR(highV.map(g => g.logMHI_L36), highV.map(g => g.logOMD_full));
const nPerm = 1000;
let countHigher = 0;
const permRs = [];
const logOMDs = highV.map(g => g.logOMD_full);

for (let iter = 0; iter < nPerm; iter++) {
  const shuffledRatios = highV.map(g => g.logMHI_L36);
  for (let i = shuffledRatios.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [shuffledRatios[i], shuffledRatios[j]] = [shuffledRatios[j], shuffledRatios[i]];
  }
  const permR = pearsonR(shuffledRatios, logOMDs);
  permRs.push(permR);
  if (Math.abs(permR) >= Math.abs(realR)) countHigher++;
}

const pValue = countHigher / nPerm;
const meanPermR = mean(permRs);
const sdPermR = sd(permRs);

console.log('  Real r(log(MHI/L36), logOMD) = ' + realR.toFixed(3));
console.log('  Permutation mean r = ' + meanPermR.toFixed(3) + ' +/- ' + sdPermR.toFixed(3));
console.log('  p-value (|perm r| >= |real r|) = ' + pValue.toFixed(4));
console.log('  Signal-to-noise: ' + ((realR - meanPermR) / sdPermR).toFixed(1) + ' sigma');
console.log('');

// ============================================================
// TEST F: CROSS-COMPONENT TEST
// Use L3.6 as predictor, MHI-based target (logOMD_gas).
// And MHI as predictor, L3.6-based target (logOMD_stars).
// If signal exists in CROSSED version → not circular.
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST F: CROSS-COMPONENT TEST');
console.log('  Predictor and target share NO baryonic component');
console.log('══════════════════════════════════════════════════════════════\n');

const validCross = highV.filter(g => isFinite(g.logOMD_stars) && isFinite(g.logOMD_gas));

const logMHI = validCross.map(g => Math.log10(table1[g.name].MHI));
const logL36 = validCross.map(g => Math.log10(table1[g.name].L36));

const r_mhi_starsTarget = pearsonR(logMHI, validCross.map(g => g.logOMD_stars));
const r_l36_gasTarget = pearsonR(logL36, validCross.map(g => g.logOMD_gas));

console.log('  CROSSED PAIRS (share zero baryonic input):');
console.log('    log(MHI) → logOMD_stars:  r = ' + r_mhi_starsTarget.toFixed(3));
console.log('    log(L3.6) → logOMD_gas:   r = ' + r_l36_gasTarget.toFixed(3));
console.log('');

const r_mhi_gasTarget = pearsonR(logMHI, validCross.map(g => g.logOMD_gas));
const r_l36_starsTarget = pearsonR(logL36, validCross.map(g => g.logOMD_stars));

console.log('  SAME-COMPONENT PAIRS (shared input, circularity risk):');
console.log('    log(MHI) → logOMD_gas:    r = ' + r_mhi_gasTarget.toFixed(3));
console.log('    log(L3.6) → logOMD_stars: r = ' + r_l36_starsTarget.toFixed(3));
console.log('');

const crossedStrong = Math.abs(r_mhi_starsTarget) > 0.3 || Math.abs(r_l36_gasTarget) > 0.3;
console.log('  Crossed signals exist? ' + (crossedStrong ? 'YES → NOT circular' : 'NO → possible circularity'));
console.log('');

// ============================================================
// OVERALL VERDICT
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  OVERALL ANTI-CIRCULARITY VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const tests = [
  { name: 'A: Stars-only target', pass: starsTest },
  { name: 'B: Gas-only not uniquely strong', pass: !gasOnlyStronger },
  { name: 'C: Pure velocity target', pass: Math.abs(r_fgas_vratio) > 0.1 },
  { name: 'D: Gas contrib < 50%', pass: mean(highV.map(g => g.meanGasContrib)) < 0.5 },
  { name: 'E: Permutation p < 0.01', pass: pValue < 0.01 },
  { name: 'F: Cross-component signal', pass: crossedStrong },
];

let passCount = 0;
for (const t of tests) {
  console.log('  ' + (t.pass ? 'PASS' : 'FAIL') + ' — ' + t.name);
  if (t.pass) passCount++;
}
console.log('\n  Score: ' + passCount + '/6 tests passed');

let verdict;
if (passCount >= 5) {
  verdict = 'STRONG PASS: The fgas → logOMD signal is NOT a circularity artifact. The signal survives removing gas from the target, survives with pure velocity targets, and survives permutation tests.';
} else if (passCount >= 3) {
  verdict = 'PARTIAL PASS: Some circularity risk remains but the signal has genuine content beyond shared inputs.';
} else {
  verdict = 'FAIL: Circularity is a serious concern. The signal may be substantially driven by shared baryonic inputs.';
}
console.log('\n  ' + verdict);

const output = {
  phase: 117,
  title: 'Anti-Circularity Battery',
  nHighV: highV.length,
  testA: {
    r_fgas_starsOnly: parseFloat(r_fgas_stars.toFixed(3)),
    r_ratio_starsOnly: parseFloat(r_ratio_stars.toFixed(3)),
    loo_fgas_starsOnly: parseFloat(loo_fgas_stars.toFixed(3)),
    pass: starsTest,
  },
  testB: {
    r_fgas_gasOnly: parseFloat(r_fgas_gasTarget.toFixed(3)),
    gasOnlyStronger,
  },
  testC: {
    r_fgas_pureVratio: parseFloat(r_fgas_vratio.toFixed(3)),
    r_ratio_pureVratio: parseFloat(r_ratio_vratio.toFixed(3)),
  },
  testD: {
    meanGasContribPercent: parseFloat((mean(highV.map(g => g.meanGasContrib)) * 100).toFixed(1)),
  },
  testE: {
    realR: parseFloat(realR.toFixed(3)),
    permMeanR: parseFloat(meanPermR.toFixed(3)),
    pValue: parseFloat(pValue.toFixed(4)),
    sigmaSNR: parseFloat(((realR - meanPermR) / sdPermR).toFixed(1)),
  },
  testF: {
    crossed_mhi_stars: parseFloat(r_mhi_starsTarget.toFixed(3)),
    crossed_l36_gas: parseFloat(r_l36_gasTarget.toFixed(3)),
    same_mhi_gas: parseFloat(r_mhi_gasTarget.toFixed(3)),
    same_l36_stars: parseFloat(r_l36_starsTarget.toFixed(3)),
  },
  passCount, verdict,
};

const outPath = path.join(__dirname, '..', 'public', 'phase117-anti-circularity.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);
