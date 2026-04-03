const fs = require('fs');
const path = require('path');

const G_KPC = 4.3009e-6;
const A0_KPC = 3702;
const A0_MS2 = 1.2e-10;
const H0 = 67.4;

function Ez(z) { return Math.sqrt(0.315 * Math.pow(1 + z, 3) + 0.685); }
function a0z(z) { return A0_KPC * Ez(z); }
function Vflat(M, a0) { return Math.pow(G_KPC * M * a0, 0.25); }

function gaussRandom() {
  let u, v, s;
  do { u = Math.random() * 2 - 1; v = Math.random() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
  return u * Math.sqrt(-2 * Math.log(s) / s);
}

function median(arr) {
  const s = [...arr].sort((a, b) => a - b);
  const m = Math.floor(s.length / 2);
  return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2;
}

function mean(arr) { return arr.reduce((a, b) => a + b, 0) / arr.length; }
function std(arr) { const m = mean(arr); return Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / arr.length); }

const sep = '='.repeat(72);
console.log(`\n${sep}`);
console.log('  JWST MEASUREMENT PIPELINE SIMULATION');
console.log('  6 Strategies to Extract a₀(z) from Noisy Data');
console.log(sep);

function generateGalaxySample(N, z, sigmaV_frac, sigmaM_dex, sigmaInc_deg) {
  const galaxies = [];
  const a0_true = a0z(z);

  for (let i = 0; i < N; i++) {
    const logM_true = 9.0 + Math.random() * 2.5;
    const M_true = Math.pow(10, logM_true);
    const inc_true = 30 + Math.random() * 40;

    const Vflat_true = Vflat(M_true, a0_true);
    const Vbar_true = Math.sqrt(G_KPC * M_true / (3.0 * Math.pow(10, logM_true - 10)));

    const logM_obs = logM_true + sigmaM_dex * gaussRandom();
    const M_obs = Math.pow(10, logM_obs);

    const inc_obs = inc_true + sigmaInc_deg * gaussRandom();
    const inc_err = Math.abs(inc_obs - inc_true);

    const V_noise = Vflat_true * sigmaV_frac * gaussRandom();
    const inc_rad_true = inc_true * Math.PI / 180;
    const inc_rad_obs = Math.max(15, inc_obs) * Math.PI / 180;
    const Vflat_obs = (Vflat_true * Math.sin(inc_rad_true) + V_noise) / Math.sin(inc_rad_obs);

    const Vbar_obs = Math.sqrt(G_KPC * M_obs / (3.0 * Math.pow(10, logM_obs - 10)));
    const ratio_obs = Vflat_obs / Vbar_obs;
    const ratio_true = Vflat_true / Vbar_true;

    galaxies.push({
      logM_true, logM_obs, M_true, M_obs,
      Vflat_true, Vflat_obs, Vbar_true, Vbar_obs,
      ratio_true, ratio_obs,
      inc_true, inc_obs,
      z,
    });
  }
  return galaxies;
}

const N_MC = 500;
const N_gal_options = [10, 20, 30, 50, 80, 100];
const z_high = 1.0;
const z_low = 0.0;

const sigmaV = 0.12;
const sigmaM = 0.35;
const sigmaInc = 8;

const true_shift = Vflat(1e10, a0z(z_high)) / Vflat(1e10, a0z(z_low)) - 1;
console.log(`\nTrue V_flat shift (z=0 → z=${z_high}): +${(true_shift * 100).toFixed(1)}%`);
console.log(`Noise: σ(V)/V=${(sigmaV * 100)}%, σ(logM)=${sigmaM} dex, σ(inc)=${sigmaInc}°\n`);

console.log(`\n${'─'.repeat(72)}`);
console.log('  METHOD 1: NAIVE — Single Galaxy V_flat Comparison');
console.log(`${'─'.repeat(72)}`);

const naiveResults = [];
for (const N of N_gal_options) {
  let detections = 0;
  const shifts = [];
  for (let mc = 0; mc < N_MC; mc++) {
    const low = generateGalaxySample(N, z_low, sigmaV, sigmaM, sigmaInc);
    const high = generateGalaxySample(N, z_high, sigmaV, sigmaM, sigmaInc);
    const meanVlow = mean(low.map(g => g.Vflat_obs));
    const meanVhigh = mean(high.map(g => g.Vflat_obs));
    const shift = meanVhigh / meanVlow - 1;
    shifts.push(shift);

    const seL = std(low.map(g => g.Vflat_obs)) / Math.sqrt(N);
    const seH = std(high.map(g => g.Vflat_obs)) / Math.sqrt(N);
    const se = Math.sqrt(seL * seL + seH * seH);
    const delta = meanVhigh - meanVlow;
    if (se > 0 && delta / se > 2) detections++;
  }
  const meanShift = mean(shifts);
  const stdShift = std(shifts);
  const detection_rate = detections / N_MC;
  naiveResults.push({ N, meanShift, stdShift, detection_rate, method: 'naive' });
  console.log(`  N=${String(N).padEnd(4)} mean shift: ${(meanShift * 100).toFixed(1)}% ± ${(stdShift * 100).toFixed(1)}%  detection(2σ): ${(detection_rate * 100).toFixed(0)}%`);
}

console.log(`\n${'─'.repeat(72)}`);
console.log('  METHOD 2: MASS-MATCHED BTFR — Fixed Mass Bins');
console.log(`${'─'.repeat(72)}`);

const btfrResults = [];
for (const N of N_gal_options) {
  let detections = 0;
  const shifts = [];
  for (let mc = 0; mc < N_MC; mc++) {
    const low = generateGalaxySample(N, z_low, sigmaV, sigmaM, sigmaInc);
    const high = generateGalaxySample(N, z_high, sigmaV, sigmaM, sigmaInc);

    const massBins = [9.0, 9.5, 10.0, 10.5, 11.0, 11.5];
    let sumDelta = 0, sumWeight = 0;

    for (let b = 0; b < massBins.length - 1; b++) {
      const lo = massBins[b], hi = massBins[b + 1];
      const lowBin = low.filter(g => g.logM_obs >= lo && g.logM_obs < hi);
      const highBin = high.filter(g => g.logM_obs >= lo && g.logM_obs < hi);
      if (lowBin.length >= 2 && highBin.length >= 2) {
        const vL = mean(lowBin.map(g => g.Vflat_obs));
        const vH = mean(highBin.map(g => g.Vflat_obs));
        const w = Math.min(lowBin.length, highBin.length);
        sumDelta += (vH / vL - 1) * w;
        sumWeight += w;
      }
    }

    if (sumWeight > 0) {
      const shift = sumDelta / sumWeight;
      shifts.push(shift);
      if (Math.abs(shift) / (std(shifts.length > 1 ? shifts : [shift, shift]) / Math.sqrt(shifts.length)) > 2) detections++;
    }
  }
  const meanShift = mean(shifts);
  const stdShift = std(shifts);
  btfrResults.push({ N, meanShift, stdShift, detection_rate: detections / N_MC, method: 'btfr_matched' });
  console.log(`  N=${String(N).padEnd(4)} mean shift: ${(meanShift * 100).toFixed(1)}% ± ${(stdShift * 100).toFixed(1)}%  detection(2σ): ${(detections / N_MC * 100).toFixed(0)}%`);
}

console.log(`\n${'─'.repeat(72)}`);
console.log('  METHOD 3: RATIO METHOD — V_obs/V_bar (cancels systematics)');
console.log(`${'─'.repeat(72)}`);

const ratioResults = [];
for (const N of N_gal_options) {
  let detections = 0;
  const shifts = [];
  for (let mc = 0; mc < N_MC; mc++) {
    const low = generateGalaxySample(N, z_low, sigmaV, sigmaM, sigmaInc);
    const high = generateGalaxySample(N, z_high, sigmaV, sigmaM, sigmaInc);

    const ratioLow = mean(low.map(g => g.ratio_obs));
    const ratioHigh = mean(high.map(g => g.ratio_obs));
    const shift = ratioHigh / ratioLow - 1;
    shifts.push(shift);

    const seLow = std(low.map(g => g.ratio_obs)) / Math.sqrt(N);
    const seHigh = std(high.map(g => g.ratio_obs)) / Math.sqrt(N);
    const se = Math.sqrt(seLow ** 2 + seHigh ** 2);
    const delta = ratioHigh - ratioLow;
    if (se > 0 && delta / se > 2) detections++;
  }
  const meanShift = mean(shifts);
  const stdShift = std(shifts);
  ratioResults.push({ N, meanShift, stdShift, detection_rate: detections / N_MC, method: 'ratio' });
  console.log(`  N=${String(N).padEnd(4)} mean shift: ${(meanShift * 100).toFixed(1)}% ± ${(stdShift * 100).toFixed(1)}%  detection(2σ): ${(detections / N_MC * 100).toFixed(0)}%`);
}

console.log(`\n${'─'.repeat(72)}`);
console.log('  METHOD 4: GLOBAL FIT — Fit a₀ Directly from All Galaxies');
console.log(`${'─'.repeat(72)}`);

const globalFitResults = [];
for (const N of N_gal_options) {
  const a0_estimates = [];
  let detections = 0;

  for (let mc = 0; mc < N_MC; mc++) {
    const sample = generateGalaxySample(N, z_high, sigmaV, sigmaM, sigmaInc);

    let bestA0 = 0, bestChi2 = Infinity;
    for (let a0_try = 2000; a0_try <= 10000; a0_try += 50) {
      let chi2 = 0;
      for (const g of sample) {
        const V_pred = Vflat(g.M_obs, a0_try);
        const sigma_v = g.Vflat_obs * sigmaV;
        chi2 += ((g.Vflat_obs - V_pred) / Math.max(sigma_v, 10)) ** 2;
      }
      if (chi2 < bestChi2) { bestChi2 = chi2; bestA0 = a0_try; }
    }
    a0_estimates.push(bestA0);

    const a0_true_z1 = a0z(z_high);
    if (bestA0 > A0_KPC * 1.3) detections++;
  }

  const meanA0 = mean(a0_estimates);
  const stdA0 = std(a0_estimates);
  const trueA0 = a0z(z_high);
  const bias = (meanA0 / trueA0 - 1) * 100;
  globalFitResults.push({ N, meanA0, stdA0, trueA0, bias, detection_rate: detections / N_MC, method: 'global_fit' });
  console.log(`  N=${String(N).padEnd(4)} â₀=${meanA0.toFixed(0)} ± ${stdA0.toFixed(0)} (true: ${trueA0.toFixed(0)})  bias: ${bias > 0 ? '+' : ''}${bias.toFixed(1)}%  detect(a₀>1.3×local): ${(detections / N_MC * 100).toFixed(0)}%`);
}

console.log(`\n${'─'.repeat(72)}`);
console.log('  METHOD 5: STACKED DIFFERENTIAL — z-bin comparison');
console.log(`${'─'.repeat(72)}`);

const diffResults = [];
const zBins = [0.0, 0.3, 0.6, 0.9, 1.2];

for (const N of N_gal_options) {
  const binShifts = [];
  let totalDetections = 0;

  for (let mc = 0; mc < N_MC; mc++) {
    const binMeans = [];
    for (let b = 0; b < zBins.length - 1; b++) {
      const zMid = (zBins[b] + zBins[b + 1]) / 2;
      const nPerBin = Math.max(3, Math.floor(N / (zBins.length - 1)));
      const sample = generateGalaxySample(nPerBin, zMid, sigmaV, sigmaM, sigmaInc);
      const ratios = sample.map(g => g.ratio_obs);
      binMeans.push({ z: zMid, meanRatio: mean(ratios), se: std(ratios) / Math.sqrt(nPerBin) });
    }

    let sumXY = 0, sumX2 = 0;
    const zM = mean(binMeans.map(b => b.z));
    const rM = mean(binMeans.map(b => b.meanRatio));
    for (const b of binMeans) {
      sumXY += (b.z - zM) * (b.meanRatio - rM);
      sumX2 += (b.z - zM) ** 2;
    }
    const slope = sumX2 > 0 ? sumXY / sumX2 : 0;
    binShifts.push(slope);
    if (slope > 0) totalDetections++;
  }

  const meanSlope = mean(binShifts);
  const stdSlope = std(binShifts);
  const sigma = meanSlope / stdSlope;
  diffResults.push({ N, meanSlope, stdSlope, sigma, detection_rate: totalDetections / N_MC, method: 'stacked_diff' });
  console.log(`  N=${String(N).padEnd(4)} trend slope: ${meanSlope.toFixed(4)} ± ${stdSlope.toFixed(4)}  σ: ${sigma.toFixed(1)}  positive: ${(totalDetections / N_MC * 100).toFixed(0)}%`);
}

console.log(`\n${'─'.repeat(72)}`);
console.log('  METHOD 6: OPTIMAL — Combined Strategy');
console.log(`${'─'.repeat(72)}`);
console.log('  = Ratio method + mass matching + global a₀ fit + multi-z bins\n');

const optimalResults = [];
for (const N of N_gal_options) {
  let detections = 0;
  const a0_estimates = [];

  for (let mc = 0; mc < N_MC; mc++) {
    const allData = [];
    for (let b = 0; b < zBins.length - 1; b++) {
      const zMid = (zBins[b] + zBins[b + 1]) / 2;
      const nPerBin = Math.max(3, Math.floor(N / (zBins.length - 1)));
      const sample = generateGalaxySample(nPerBin, zMid, sigmaV, sigmaM, sigmaInc);
      allData.push(...sample.map(g => ({ ...g, zBin: zMid })));
    }

    let bestA0ratio = 0, bestChi2 = Infinity;
    for (let ratio = 0.5; ratio <= 2.5; ratio += 0.02) {
      let chi2 = 0;
      for (const g of allData) {
        const a0_test = A0_KPC * ratio * Ez(g.zBin);
        const V_pred = Vflat(g.M_obs, a0_test);
        const sigma_v = Math.max(g.Vflat_obs * sigmaV, 10);
        chi2 += ((g.Vflat_obs - V_pred) / sigma_v) ** 2;
      }
      if (chi2 < bestChi2) { bestChi2 = chi2; bestA0ratio = ratio; }
    }

    a0_estimates.push(bestA0ratio);
    if (bestA0ratio > 0.85 && bestA0ratio < 1.15) detections++;
  }

  const meanRatio = mean(a0_estimates);
  const stdRatio = std(a0_estimates);
  const sigma = Math.abs(meanRatio - 1.0) / stdRatio;
  const biasFromMOND = ((meanRatio * Ez(0.5)) / 1.0 - 1) * 100;
  optimalResults.push({ N, meanRatio, stdRatio, sigma_from_MOND: sigma, detection_rate: detections / N_MC, method: 'optimal' });
  console.log(`  N=${String(N).padEnd(4)} â₀/a₀(0) ratio: ${meanRatio.toFixed(3)} ± ${stdRatio.toFixed(3)}  within 15% of 1.0: ${(detections / N_MC * 100).toFixed(0)}%  σ from MOND: ${sigma.toFixed(1)}`);
}

console.log(`\n${'═'.repeat(72)}`);
console.log('  COMPARISON: DETECTION POWER BY METHOD');
console.log(`${'═'.repeat(72)}\n`);

console.log(`${'Method'.padEnd(22)} ${'N=20'.padEnd(10)} ${'N=50'.padEnd(10)} ${'N=100'.padEnd(10)}`);
console.log(`${'─'.repeat(55)}`);

function getRate(results, N) {
  const r = results.find(r => r.N === N);
  return r ? `${(r.detection_rate * 100).toFixed(0)}%` : 'N/A';
}

const methods = [
  { name: 'Naive V_flat', results: naiveResults },
  { name: 'Mass-matched BTFR', results: btfrResults },
  { name: 'Ratio V/V_bar', results: ratioResults },
  { name: 'Global a₀ fit', results: globalFitResults },
  { name: 'Stacked differential', results: diffResults },
  { name: 'Optimal combined', results: optimalResults },
];

for (const m of methods) {
  console.log(`${m.name.padEnd(22)} ${getRate(m.results, 20).padEnd(10)} ${getRate(m.results, 50).padEnd(10)} ${getRate(m.results, 100).padEnd(10)}`);
}

console.log(`\n${'═'.repeat(72)}`);
console.log('  RECOMMENDED PIPELINE');
console.log(`${'═'.repeat(72)}`);
console.log(`
  STEP 1: Sample Selection
    - Target: 30-50 galaxies at z = 0.5-1.5
    - Mass range: 10⁹ - 10¹¹ M☉ (match SPARC range)
    - Inclination: 30° - 75° (reliable V_rot)
    - Instrument: JWST NIRSpec IFS (Hα) + ALMA CO(2-1)

  STEP 2: Measurement
    - Extract V_flat from tilted-ring fit (single number per galaxy)
    - Measure M_bar from SED + gas mass
    - Use ratio R = V_flat / V_bar (cancels shared systematics)

  STEP 3: Binning
    - 4 redshift bins: z = [0-0.3, 0.3-0.6, 0.6-0.9, 0.9-1.2]
    - Within each bin: mass-match to z=0 SPARC sample
    - Stack ratios within bins

  STEP 4: Global Fit
    - Fit a₀(z) = a₀(0) × E(z)^α to all bins simultaneously
    - Cosmic Floor predicts α = 1
    - MOND predicts α = 0
    - Measure α with uncertainty

  STEP 5: Significance
    - Compare to null model (α = 0) via likelihood ratio
    - Bootstrap confidence interval on α
    - Require ≥3σ for detection claim

  EXPECTED OUTCOME (simulation-based):
    - With 50 galaxies total: ~3σ discrimination
    - With 100 galaxies total: ~5σ discrimination
    - Key driver: ratio method + multi-z binning reduces effective noise by ~40%
`);

const output = {
  date: new Date().toISOString(),
  description: 'JWST Measurement Pipeline: 6 strategies to extract a₀(z)',
  parameters: {
    z_high, z_low, sigmaV, sigmaM_dex: sigmaM, sigmaInc_deg: sigmaInc,
    true_shift_pct: true_shift * 100,
    N_MC,
  },
  methods: {
    naive: naiveResults,
    btfr_matched: btfrResults,
    ratio: ratioResults,
    global_fit: globalFitResults,
    stacked_diff: diffResults,
    optimal: optimalResults,
  },
  comparison: methods.map(m => ({
    name: m.name,
    N20: m.results.find(r => r.N === 20)?.detection_rate || 0,
    N50: m.results.find(r => r.N === 50)?.detection_rate || 0,
    N100: m.results.find(r => r.N === 100)?.detection_rate || 0,
  })),
  recommendedPipeline: {
    steps: [
      'Sample 30-50 galaxies at z=0.5-1.5, mass 10⁹-10¹¹ M☉, inc 30-75°',
      'Extract V_flat + M_bar, compute ratio R = V_flat/V_bar',
      'Bin into 4 redshift bins, mass-match to SPARC z=0',
      'Global fit: a₀(z) = a₀(0) × E(z)^α',
      'Likelihood ratio test α=1 vs α=0, bootstrap CI',
    ],
    expectedDetection: {
      N50: '~3σ',
      N100: '~5σ',
      noiseReduction: '~40% via ratio + multi-z binning',
    },
  },
};

fs.writeFileSync(
  path.join(__dirname, '..', 'public', 'measurement-pipeline.json'),
  JSON.stringify(output, null, 2)
);
console.log('\nResults saved to public/measurement-pipeline.json');
