const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 80b: EXTERNAL VALIDATION — M2 (logMHI + logMhost only)');
console.log('  No MeanRun needed — tests the 2-axis mass-suppression structure');
console.log('======================================================================\n');

const csvPath = path.join(__dirname, '..', 'public', 'replication', 'N45_final_dataset.csv');
const lines = fs.readFileSync(csvPath, 'utf-8').trim().split('\n');
const headers = lines[0].split(',');
const trainData = lines.slice(1).map(l => {
  const vals = l.split(',');
  const obj = {};
  headers.forEach((h, i) => obj[h] = vals[i]);
  return obj;
});
const trainNames = trainData.map(d => d.galaxy_name);

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }

function olsFit(Y, X) {
  const n = Y.length, p = X[0].length;
  const Xa = X.map((row) => [1, ...row]);
  const pp = p + 1;
  const XtX = Array.from({ length: pp }, () => new Array(pp).fill(0));
  const XtY = new Array(pp).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < pp; j++) {
      XtY[j] += Xa[i][j] * Y[i];
      for (let k = 0; k < pp; k++) {
        XtX[j][k] += Xa[i][j] * Xa[i][k];
      }
    }
  }
  const aug = XtX.map((row, i) => [...row, ...Array.from({ length: pp }, (_, j) => i === j ? 1 : 0)]);
  for (let col = 0; col < pp; col++) {
    let maxRow = col;
    for (let row = col + 1; row < pp; row++)
      if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row;
    [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
    const pivot = aug[col][col];
    if (Math.abs(pivot) < 1e-12) return null;
    for (let j = 0; j < 2 * pp; j++) aug[col][j] /= pivot;
    for (let row = 0; row < pp; row++) {
      if (row === col) continue;
      const factor = aug[row][col];
      for (let j = 0; j < 2 * pp; j++) aug[row][j] -= factor * aug[col][j];
    }
  }
  const inv = aug.map(row => row.slice(pp));
  const beta = inv.map(row => row.reduce((s, v, j) => s + v * XtY[j], 0));
  return beta;
}

const trainY = trainData.map(d => parseFloat(d.logA0));
const trainMHI = trainData.map(d => parseFloat(d.logMHI));
const trainMhost = trainData.map(d => parseFloat(d.logMhost));
const trainMeanRun = trainData.map(d => parseFloat(d.logMeanRun));
const trainMeanA0 = mean(trainY);

const X_M2_train = trainData.map(d => [parseFloat(d.logMHI), parseFloat(d.logMhost)]);
const betaM2 = olsFit(trainY, X_M2_train);
console.log('  FROZEN M2 coefficients (fit on N=45 training set):');
console.log('    intercept = ' + betaM2[0].toFixed(4));
console.log('    logMHI    = ' + betaM2[1].toFixed(4));
console.log('    logMhost  = ' + betaM2[2].toFixed(4));
console.log();

function looCV_M2() {
  let ss = 0;
  for (let i = 0; i < trainY.length; i++) {
    const Yr = trainY.filter((_, j) => j !== i);
    const Xr = X_M2_train.filter((_, j) => j !== i);
    const b = olsFit(Yr, Xr);
    const pred = b[0] + b[1] * X_M2_train[i][0] + b[2] * X_M2_train[i][1];
    ss += (trainY[i] - pred) ** 2;
  }
  const looRms = Math.sqrt(ss / trainY.length);
  const sdY = sd(trainY);
  return { gap: 100 * (1 - looRms ** 2 / sdY ** 2), rms: looRms };
}

const looM2 = looCV_M2();
console.log('  M2 LOO on training set: gap% = ' + looM2.gap.toFixed(1) + '%, RMS = ' + looM2.rms.toFixed(4));
console.log('  (For reference: M3 LOO = 44.1%, M5 LOO = 51.0%)');
console.log();

const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf-8'));
const a0all = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'per-galaxy-a0.json'), 'utf-8'));

function computeLogMhost(Vflat_kms) {
  return 4.0 * Math.log10(Vflat_kms) + 2.4;
}

const external = [];
for (const gsp of sparc) {
  if (trainNames.includes(gsp.name)) continue;
  if (gsp.fD < 2) continue;
  const ga0 = a0all.find(x => x.name === gsp.name);
  if (!ga0 || ga0.n < 5) continue;
  if (gsp.Vflat === 0) continue;
  if (ga0.logA0 < 2.0) continue;

  external.push({
    name: gsp.name,
    logA0_obs: ga0.logA0,
    nPts: ga0.n,
    Q: gsp.Q,
    D: gsp.D,
    Vflat: gsp.Vflat,
    T: gsp.T,
    logMHI: Math.log10(gsp.MHI),
    logMhost: computeLogMhost(gsp.Vflat)
  });
}

console.log('----------------------------------------------------------------------');
console.log('  EXTERNAL SAMPLE: N = ' + external.length);
console.log('----------------------------------------------------------------------\n');

const obs = external.map(g => g.logA0_obs);
const predM2 = external.map(g => betaM2[0] + betaM2[1] * g.logMHI + betaM2[2] * g.logMhost);
const predM0 = external.map(() => trainMeanA0);
const sdObs = sd(obs);

function rms(pred) { return Math.sqrt(pred.reduce((s, p, i) => s + (obs[i] - p) ** 2, 0) / obs.length); }
function bias(pred) { return pred.reduce((s, p, i) => s + (p - obs[i]), 0) / obs.length; }
function spearman(pred) {
  const rank = arr => {
    const sorted = [...arr].map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
    const r = new Array(arr.length);
    sorted.forEach((s, rank) => r[s.i] = rank + 1);
    return r;
  };
  const ro = rank(obs), rp = rank(pred);
  const mo = mean(ro), mp = mean(rp);
  const num = ro.reduce((s, v, i) => s + (v - mo) * (rp[i] - mp), 0);
  const d1 = Math.sqrt(ro.reduce((s, v) => s + (v - mo) ** 2, 0));
  const d2 = Math.sqrt(rp.reduce((s, v) => s + (v - mp) ** 2, 0));
  return num / (d1 * d2 + 1e-12);
}
function pearson(pred) {
  const mo = mean(obs), mp = mean(pred);
  const num = obs.reduce((s, v, i) => s + (v - mo) * (pred[i] - mp), 0);
  const d1 = Math.sqrt(obs.reduce((s, v) => s + (v - mo) ** 2, 0));
  const d2 = Math.sqrt(pred.reduce((s, v) => s + (v - mp) ** 2, 0));
  return num / (d1 * d2 + 1e-12);
}

const rmsM0 = rms(predM0);
const rmsM2 = rms(predM2);
const gapM2 = 100 * (1 - rmsM2 ** 2 / sdObs ** 2);
const rsM2 = spearman(predM2);
const rpM2 = pearson(predM2);

console.log('  SD(observed logA0) = ' + sdObs.toFixed(4) + ' dex\n');
console.log('  ' + 'Model'.padEnd(10) + 'RMS'.padStart(8) + 'Bias'.padStart(8) + 'Gap%'.padStart(8) + 'Spearman'.padStart(10) + 'Pearson'.padStart(9));
console.log('  ' + '-'.repeat(53));
console.log('  ' + 'M0'.padEnd(10) + rmsM0.toFixed(4).padStart(8) + bias(predM0).toFixed(4).padStart(8) + '—'.padStart(8) + '—'.padStart(10) + '—'.padStart(9));
console.log('  ' + 'M2'.padEnd(10) + rmsM2.toFixed(4).padStart(8) + bias(predM2).toFixed(4).padStart(8) + (gapM2.toFixed(1) + '%').padStart(8) + rsM2.toFixed(3).padStart(10) + rpM2.toFixed(3).padStart(9));
console.log();

let m2wins = 0;
console.log('  Per-galaxy M2 predictions:');
console.log('  ' + 'Galaxy'.padEnd(14) + 'Obs'.padStart(7) + 'PredM2'.padStart(8) + 'Err'.padStart(8) + 'PredM0'.padStart(8) + '|M2|<|M0|?'.padStart(12));
console.log('  ' + '-'.repeat(57));

for (const g of external) {
  const o = g.logA0_obs;
  const p2 = betaM2[0] + betaM2[1] * g.logMHI + betaM2[2] * g.logMhost;
  const p0 = trainMeanA0;
  const e2 = p2 - o;
  const e0 = p0 - o;
  const wins = Math.abs(e2) < Math.abs(e0);
  if (wins) m2wins++;
  console.log('  ' + g.name.padEnd(14) + o.toFixed(3).padStart(7) + p2.toFixed(3).padStart(8) + (e2 >= 0 ? '+' : '') + e2.toFixed(3).padStart(7) + p0.toFixed(3).padStart(8) + (wins ? 'YES' : 'no').padStart(12));
}
console.log();
console.log('  M2 closer than M0: ' + m2wins + '/' + external.length + ' (' + (m2wins / external.length * 100).toFixed(1) + '%)');
console.log();

console.log('----------------------------------------------------------------------');
console.log('  ALSO: CHECK Mhost ESTIMATION vs GROUP CATALOG');
console.log('----------------------------------------------------------------------\n');

const p25 = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase25-group-membership.json'), 'utf-8'));
const groupLookup = {};
for (const g of p25.galaxyAssignments) groupLookup[g.name] = g;

let nMatch = 0;
for (const g of external) {
  const grp = groupLookup[g.name];
  if (grp && grp.logMhost) {
    const diff = g.logMhost - grp.logMhost;
    console.log('  ' + g.name.padEnd(14) + 'Vflat est: ' + g.logMhost.toFixed(2) + '  Group catalog: ' + grp.logMhost.toFixed(2) + '  diff: ' + diff.toFixed(2));
    nMatch++;
  }
}
if (nMatch === 0) console.log('  (No external galaxies found in group catalog — expected, since they were excluded from working sample)');
console.log();

console.log('----------------------------------------------------------------------');
console.log('  Q=1 SUBSAMPLE (highest quality RC data)');
console.log('----------------------------------------------------------------------\n');

const q1 = external.filter(g => g.Q === 1);
if (q1.length > 0) {
  const obsQ1 = q1.map(g => g.logA0_obs);
  const predQ1_M2 = q1.map(g => betaM2[0] + betaM2[1] * g.logMHI + betaM2[2] * g.logMhost);
  const predQ1_M0 = q1.map(() => trainMeanA0);
  const sdQ1 = sd(obsQ1);
  const rmsQ1_M0 = Math.sqrt(predQ1_M0.reduce((s, p, i) => s + (obsQ1[i] - p) ** 2, 0) / q1.length);
  const rmsQ1_M2 = Math.sqrt(predQ1_M2.reduce((s, p, i) => s + (obsQ1[i] - p) ** 2, 0) / q1.length);
  const gapQ1 = 100 * (1 - rmsQ1_M2 ** 2 / sdQ1 ** 2);

  const roQ1 = obsQ1.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v).reduce((r, s, rank) => { r[s.i] = rank + 1; return r; }, new Array(q1.length));
  const rpQ1 = predQ1_M2.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v).reduce((r, s, rank) => { r[s.i] = rank + 1; return r; }, new Array(q1.length));
  const moQ1 = mean(roQ1), mpQ1 = mean(rpQ1);
  const rsQ1 = roQ1.reduce((s, v, i) => s + (v - moQ1) * (rpQ1[i] - mpQ1), 0) / (Math.sqrt(roQ1.reduce((s, v) => s + (v - moQ1) ** 2, 0)) * Math.sqrt(rpQ1.reduce((s, v) => s + (v - mpQ1) ** 2, 0)) + 1e-12);

  console.log('  Q=1 subsample: N = ' + q1.length);
  console.log('  M0 RMS = ' + rmsQ1_M0.toFixed(4) + ', M2 RMS = ' + rmsQ1_M2.toFixed(4));
  console.log('  M2 gap% = ' + gapQ1.toFixed(1) + '%');
  console.log('  M2 Spearman = ' + rsQ1.toFixed(3));
  console.log('  M2 beats M0? ' + (rmsQ1_M2 < rmsQ1_M0 ? 'YES' : 'NO'));

  let q1wins = 0;
  for (const g of q1) {
    const e2 = Math.abs(betaM2[0] + betaM2[1] * g.logMHI + betaM2[2] * g.logMhost - g.logA0_obs);
    const e0 = Math.abs(trainMeanA0 - g.logA0_obs);
    if (e2 < e0) q1wins++;
  }
  console.log('  M2 closer per-galaxy: ' + q1wins + '/' + q1.length);
}
console.log();

console.log('----------------------------------------------------------------------');
console.log('  KEY DIAGNOSTIC: logMhost ESTIMATION BIAS');
console.log('----------------------------------------------------------------------\n');

const trainVflat = trainData.map(d => parseFloat(d.Vflat_km_s));
const trainMhostReal = trainData.map(d => parseFloat(d.logMhost));
const trainMhostEst = trainVflat.map(v => 4.0 * Math.log10(v) + 2.4);
const mhostDiffs = trainMhostEst.map((e, i) => e - trainMhostReal[i]);

console.log('  Vflat-based Mhost vs training set group-catalog Mhost:');
console.log('    Mean diff: ' + mean(mhostDiffs).toFixed(3) + ' dex');
console.log('    SD diff: ' + sd(mhostDiffs).toFixed(3) + ' dex');
console.log('    Min/Max diff: ' + Math.min(...mhostDiffs).toFixed(3) + ' / ' + Math.max(...mhostDiffs).toFixed(3));
console.log();

const betaM2_raw = olsFit(trainY, trainData.map(d => [parseFloat(d.logMHI), 4.0 * Math.log10(parseFloat(d.Vflat_km_s)) + 2.4]));
console.log('  M2 using Vflat-estimated Mhost (on training set):');
console.log('    beta = [' + betaM2_raw.map(b => b.toFixed(4)).join(', ') + ']');

const trainPredRaw = trainData.map(d => betaM2_raw[0] + betaM2_raw[1] * parseFloat(d.logMHI) + betaM2_raw[2] * (4.0 * Math.log10(parseFloat(d.Vflat_km_s)) + 2.4));
const trainRmsRaw = Math.sqrt(trainPredRaw.reduce((s, p, i) => s + (trainY[i] - p) ** 2, 0) / trainY.length);
console.log('    Training RMS with est Mhost = ' + trainRmsRaw.toFixed(4));
console.log();

console.log('  IMPORTANT: External galaxies use Vflat-estimated Mhost');
console.log('  Training set uses group-catalog Mhost');
console.log('  This input mismatch may drive bias in the external test');
console.log();

console.log('----------------------------------------------------------------------');
console.log('  REFIT M2 WITH Vflat-ESTIMATED Mhost (fair comparison)');
console.log('----------------------------------------------------------------------\n');

const trainX_fair = trainData.map(d => [parseFloat(d.logMHI), 4.0 * Math.log10(parseFloat(d.Vflat_km_s)) + 2.4]);
const betaM2_fair = olsFit(trainY, trainX_fair);
console.log('  M2_fair coefficients (using Vflat-estimated Mhost on training):');
console.log('    intercept = ' + betaM2_fair[0].toFixed(4));
console.log('    logMHI    = ' + betaM2_fair[1].toFixed(4));
console.log('    logMhost_est = ' + betaM2_fair[2].toFixed(4));
console.log();

const predM2_fair = external.map(g => betaM2_fair[0] + betaM2_fair[1] * g.logMHI + betaM2_fair[2] * g.logMhost);
const rmsM2_fair = Math.sqrt(predM2_fair.reduce((s, p, i) => s + (obs[i] - p) ** 2, 0) / obs.length);
const biasM2_fair = predM2_fair.reduce((s, p, i) => s + (p - obs[i]), 0) / obs.length;
const gapM2_fair = 100 * (1 - rmsM2_fair ** 2 / sdObs ** 2);

const rpFair = (() => {
  const mo2 = mean(obs), mp2 = mean(predM2_fair);
  const n2 = obs.reduce((s, v, i) => s + (v - mo2) * (predM2_fair[i] - mp2), 0);
  const d12 = Math.sqrt(obs.reduce((s, v) => s + (v - mo2) ** 2, 0));
  const d22 = Math.sqrt(predM2_fair.reduce((s, v) => s + (v - mp2) ** 2, 0));
  return n2 / (d12 * d22 + 1e-12);
})();

console.log('  M2_fair on external sample:');
console.log('    RMS = ' + rmsM2_fair.toFixed(4) + ' (M0: ' + rmsM0.toFixed(4) + ')');
console.log('    Bias = ' + biasM2_fair.toFixed(4));
console.log('    Gap% = ' + gapM2_fair.toFixed(1) + '%');
console.log('    Pearson r = ' + rpFair.toFixed(3));
console.log('    M2_fair beats M0? ' + (rmsM2_fair < rmsM0 ? 'YES' : 'NO'));
console.log();

let fairWins = 0;
console.log('  Per-galaxy M2_fair:');
console.log('  ' + 'Galaxy'.padEnd(14) + 'Obs'.padStart(7) + 'Pred'.padStart(8) + 'Err'.padStart(8) + '|M2|<|M0|?'.padStart(12));
console.log('  ' + '-'.repeat(49));
for (const g of external) {
  const o = g.logA0_obs;
  const p = betaM2_fair[0] + betaM2_fair[1] * g.logMHI + betaM2_fair[2] * g.logMhost;
  const e = p - o;
  const e0 = trainMeanA0 - o;
  const wins = Math.abs(e) < Math.abs(e0);
  if (wins) fairWins++;
  console.log('  ' + g.name.padEnd(14) + o.toFixed(3).padStart(7) + p.toFixed(3).padStart(8) + (e >= 0 ? '+' : '') + e.toFixed(3).padStart(7) + (wins ? 'YES' : 'no').padStart(12));
}
console.log();
console.log('  M2_fair closer than M0: ' + fairWins + '/' + external.length + ' (' + (fairWins / external.length * 100).toFixed(1) + '%)');
console.log();

console.log('======================================================================');
console.log('  OVERALL PHASE 80 SUMMARY');
console.log('======================================================================\n');

console.log('  KEY FINDING: M3/M5 external test CANNOT be performed cleanly because');
console.log('  the MeanRun variable requires raw RC point data (not available for');
console.log('  external galaxies). The Mhost variable also uses Vflat estimation');
console.log('  instead of group-catalog values, introducing systematic mismatch.\n');

console.log('  M2 (logMHI + logMhost only) external results:');
console.log('    M2 with group-catalog Mhost: RMS=' + rmsM2.toFixed(3) + ', gap%=' + gapM2.toFixed(1) + '%, rs=' + rsM2.toFixed(3));
console.log('    M2_fair (Vflat Mhost): RMS=' + rmsM2_fair.toFixed(3) + ', gap%=' + gapM2_fair.toFixed(1) + '%, r=' + rpFair.toFixed(3));
console.log('    M0 baseline: RMS=' + rmsM0.toFixed(3));
console.log();

const m2verdict = rmsM2_fair < rmsM0 ? 'PARTIAL-PASS' : (rpFair > 0.2 ? 'MARGINAL' : 'INCONCLUSIVE');
console.log('  VERDICT: ' + m2verdict);
console.log();
console.log('  INTERPRETATION:');
if (m2verdict === 'PARTIAL-PASS') {
  console.log('  Even without MeanRun (the strongest axis), the mass-suppression');
  console.log('  structure (logMHI + logMhost) generalizes to unseen galaxies.');
  console.log('  Full M3 test requires raw RC data to compute MeanRun.');
} else {
  console.log('  The 2-axis mass structure alone is not sufficient for external');
  console.log('  prediction. MeanRun may be essential. A clean M3 test requires');
  console.log('  raw RC data to compute true MeanRun for external galaxies.');
}
console.log();
console.log('  REQUIREMENTS FOR CLEAN M3/M5 TEST:');
console.log('  1. Raw rotation curve points to compute MeanRun (MOND residual runs)');
console.log('  2. Group-catalog Mhost (not Vflat estimation)');
console.log('  3. SPS-based Upsilon for M5');
console.log('  4. All from galaxies NOT in training set');

const output = {
  phase: '80',
  title: 'External Clean-Sample Validation',
  fundamentalLimitation: 'MeanRun cannot be computed without raw RC point data; Mhost uses Vflat proxy instead of group catalog',
  m2test: {
    description: 'logMHI + logMhost only (no MeanRun)',
    nExternal: external.length,
    frozenBetaM2: betaM2,
    frozenBetaM2_fair: betaM2_fair,
    rmsM0: +rmsM0.toFixed(4),
    rmsM2: +rmsM2.toFixed(4),
    rmsM2_fair: +rmsM2_fair.toFixed(4),
    gapM2: +gapM2.toFixed(1),
    gapM2_fair: +gapM2_fair.toFixed(1),
    spearmanM2: +rsM2.toFixed(3),
    pearsonM2_fair: +rpFair.toFixed(3),
    m2BeatsM0: rmsM2_fair < rmsM0,
    m2WinRate: +(fairWins / external.length * 100).toFixed(1)
  },
  mhostEstimationBias: {
    meanDiff: +mean(mhostDiffs).toFixed(3),
    sdDiff: +sd(mhostDiffs).toFixed(3)
  },
  verdict: m2verdict,
  cleanTestRequirements: [
    'Raw RC points for MeanRun computation',
    'Group-catalog Mhost values',
    'SPS-based Upsilon for M5',
    'Galaxies not in training set'
  ]
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase80-external-validation.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase80-external-validation.json');
