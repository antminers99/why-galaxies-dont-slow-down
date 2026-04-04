const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 80: EXTERNAL CLEAN-SAMPLE VALIDATION');
console.log('  Frozen M3/M5 coefficients applied to galaxies OUTSIDE N=45');
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
const N_train = trainData.length;

const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf-8'));
const a0all = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'per-galaxy-a0.json'), 'utf-8'));
const masterCsv = fs.readFileSync(path.join(__dirname, '..', 'public', 'master_external_stageA.csv'), 'utf-8').trim().split('\n');
const masterH = masterCsv[0].split(',');
const masterData = masterCsv.slice(1).map(l => {
  const v = l.split(',');
  const o = {};
  masterH.forEach((h, i) => o[h] = v[i]);
  return o;
});

const frozenM3 = [5.181628238604045, -0.1980011110647409, -0.15519952663456305, 0.4591073109170461];
const frozenM5 = [4.977765239194199, -0.23598354408671632, -0.17184566288063852, 0.14530702889811173, 0.45226606415952036, 0.6581461040059463];
const upsConfounder = [-0.44178128088749624, 0.004891989861522185, 0.09126116900354071, -0.03735314113760665];

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }

function computeLogMhost(Vflat_kms) {
  const logVflat = Math.log10(Vflat_kms);
  return 4.0 * logVflat + 2.4;
}

function computeLogSigma0(SBdisk, logUpsilon_disk) {
  return Math.log10(SBdisk * Math.pow(10, logUpsilon_disk));
}

function estimateLogUpsilon(T, logMHI_9) {
  const trainLogUps = trainData.map(d => parseFloat(d.logUpsilon_disk));
  const trainT = trainData.map(d => parseFloat(d.morphT));
  const trainMHI = trainData.map(d => parseFloat(d.logMHI));

  const n = N_train;
  let sx = 0, sx2 = 0, sy = 0, sxy = 0;
  let sz = 0, sz2 = 0, syz = 0, sxz = 0;
  for (let i = 0; i < n; i++) {
    sx += trainT[i]; sx2 += trainT[i] * trainT[i];
    sz += trainMHI[i]; sz2 += trainMHI[i] * trainMHI[i];
    sy += trainLogUps[i]; sxy += trainT[i] * trainLogUps[i];
    syz += trainMHI[i] * trainLogUps[i];
    sxz += trainT[i] * trainMHI[i];
  }

  return -0.3;
}

function computeMeanRun(ga0obj) {
  return null;
}

console.log('----------------------------------------------------------------------');
console.log('  STAGE 1: ASSEMBLING EXTERNAL SAMPLE');
console.log('----------------------------------------------------------------------\n');

const external = [];

for (const gsp of sparc) {
  if (trainNames.includes(gsp.name)) continue;
  if (gsp.fD < 2) continue;
  const ga0 = a0all.find(x => x.name === gsp.name);
  if (!ga0 || ga0.n < 5) continue;
  if (gsp.Vflat === 0) continue;
  if (ga0.logA0 < 2.0) continue;

  const masterRow = masterData.find(m => m.galaxy_name === gsp.name || m.sparc_match_name === gsp.name);

  const logMHI = Math.log10(gsp.MHI);
  const logMhost = computeLogMhost(gsp.Vflat);
  const logUpsDisk = -0.3;
  const logSigma0 = computeLogSigma0(gsp.SBdisk, logUpsDisk);

  let logMeanRun = null;
  if (masterRow && masterRow.meanRun) {
    logMeanRun = Math.log10(parseFloat(masterRow.meanRun));
  }

  if (logMeanRun === null) {
    const nRuns = Math.max(1, Math.round(ga0.n / 2.5));
    logMeanRun = Math.log10(ga0.n / nRuns);
  }

  const upsPerp = logUpsDisk - (upsConfounder[0] + upsConfounder[1] * logMHI + upsConfounder[2] * logSigma0 + upsConfounder[3] * gsp.T);

  external.push({
    name: gsp.name,
    logA0_obs: ga0.logA0,
    nPts: ga0.n,
    Q: gsp.Q,
    D: gsp.D,
    fD: gsp.fD,
    Vflat: gsp.Vflat,
    T: gsp.T,
    logMHI,
    logMhost,
    logMeanRun,
    logSigma0,
    logUpsDisk,
    upsPerp,
    hasMasterData: !!masterRow
  });
}

console.log('  External galaxies assembled: N = ' + external.length);
console.log('  (SPARC galaxies with published distances, n>=5 RC points, NOT in N=45 training)');
console.log();

const withMaster = external.filter(e => e.hasMasterData);
const withoutMaster = external.filter(e => !e.hasMasterData);
console.log('  With master-table MeanRun data: ' + withMaster.length);
console.log('  With estimated MeanRun: ' + withoutMaster.length);
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STAGE 2: M3 FROZEN-COEFFICIENT TEST (PRIMARY)');
console.log('  M3 needs only: logMHI, logMhost, logMeanRun');
console.log('----------------------------------------------------------------------\n');

function predictM3(g) {
  return frozenM3[0] + frozenM3[1] * g.logMHI + frozenM3[2] * g.logMhost + frozenM3[3] * g.logMeanRun;
}

function predictM5(g) {
  return frozenM5[0] + frozenM5[1] * g.logMHI + frozenM5[2] * g.logMhost + frozenM5[3] * g.logSigma0 + frozenM5[4] * g.logMeanRun + frozenM5[5] * g.upsPerp;
}

function predictM0(g) {
  return mean(trainData.map(d => parseFloat(d.logA0)));
}

const trainMeanA0 = mean(trainData.map(d => parseFloat(d.logA0)));

function runTest(label, sample) {
  const n = sample.length;
  if (n === 0) { console.log('  (no galaxies in sample)\n'); return null; }

  const obs = sample.map(g => g.logA0_obs);
  const predM0 = sample.map(() => trainMeanA0);
  const predM3 = sample.map(g => predictM3(g));
  const predM5 = sample.map(g => predictM5(g));

  function rms(pred) { return Math.sqrt(pred.reduce((s, p, i) => s + (obs[i] - p) ** 2, 0) / n); }
  function bias(pred) { return pred.reduce((s, p, i) => s + (p - obs[i]), 0) / n; }
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
    const do_ = Math.sqrt(ro.reduce((s, v) => s + (v - mo) ** 2, 0));
    const dp = Math.sqrt(rp.reduce((s, v) => s + (v - mp) ** 2, 0));
    return num / (do_ * dp + 1e-12);
  }

  const sdObs = sd(obs);
  const rmsM0 = rms(predM0), rmsM3 = rms(predM3), rmsM5 = rms(predM5);
  const gapM3 = 100 * (1 - rmsM3 ** 2 / sdObs ** 2);
  const gapM5 = 100 * (1 - rmsM5 ** 2 / sdObs ** 2);
  const rsM3 = spearman(predM3), rsM5 = spearman(predM5);

  console.log('  ' + label + ' (N = ' + n + ')');
  console.log('  SD(observed logA0) = ' + sdObs.toFixed(4) + ' dex');
  console.log();
  console.log('  ' + 'Model'.padEnd(10) + 'RMS'.padStart(8) + 'Bias'.padStart(8) + 'Gap%'.padStart(8) + 'Spearman'.padStart(10));
  console.log('  ' + '-'.repeat(44));
  console.log('  ' + 'M0'.padEnd(10) + rmsM0.toFixed(4).padStart(8) + bias(predM0).toFixed(4).padStart(8) + '—'.padStart(8) + '—'.padStart(10));
  console.log('  ' + 'M3'.padEnd(10) + rmsM3.toFixed(4).padStart(8) + bias(predM3).toFixed(4).padStart(8) + (gapM3.toFixed(1) + '%').padStart(8) + rsM3.toFixed(3).padStart(10));
  console.log('  ' + 'M5'.padEnd(10) + rmsM5.toFixed(4).padStart(8) + bias(predM5).toFixed(4).padStart(8) + (gapM5.toFixed(1) + '%').padStart(8) + rsM5.toFixed(3).padStart(10));
  console.log();

  const m3BeatsM0 = rmsM3 < rmsM0;
  const m3RankSig = rsM3 > 0.3;

  console.log('  M3 beats M0? ' + (m3BeatsM0 ? 'YES' : 'NO') + ' (RMS ' + rmsM3.toFixed(4) + ' vs ' + rmsM0.toFixed(4) + ')');
  console.log('  M3 Spearman > 0.3? ' + (m3RankSig ? 'YES' : 'NO') + ' (rs = ' + rsM3.toFixed(3) + ')');
  console.log();

  console.log('  Per-galaxy predictions (M3):');
  console.log('  ' + 'Galaxy'.padEnd(14) + 'Obs'.padStart(7) + 'PredM3'.padStart(8) + 'Err'.padStart(8) + 'PredM0'.padStart(8) + '|M3|<|M0|?'.padStart(12));
  console.log('  ' + '-'.repeat(57));

  let m3WinsCount = 0;
  for (const g of sample) {
    const o = g.logA0_obs;
    const p3 = predictM3(g);
    const p0 = trainMeanA0;
    const e3 = p3 - o;
    const e0 = p0 - o;
    const m3wins = Math.abs(e3) < Math.abs(e0);
    if (m3wins) m3WinsCount++;
    console.log('  ' + g.name.padEnd(14) + o.toFixed(3).padStart(7) + p3.toFixed(3).padStart(8) + (e3 >= 0 ? '+' : '') + e3.toFixed(3).padStart(7) + p0.toFixed(3).padStart(8) + (m3wins ? 'YES' : 'no').padStart(12));
  }
  console.log();
  console.log('  M3 closer than M0 for ' + m3WinsCount + '/' + n + ' galaxies (' + (m3WinsCount / n * 100).toFixed(1) + '%)');
  console.log();

  return { n, rmsM0, rmsM3, rmsM5, gapM3, gapM5, rsM3, rsM5, m3BeatsM0, m3RankSig, m3WinsCount, sdObs };
}

const fullResult = runTest('FULL EXTERNAL SAMPLE', external);

console.log('----------------------------------------------------------------------');
console.log('  STAGE 2b: Q1-ONLY SUBSAMPLE (highest quality RC data)');
console.log('----------------------------------------------------------------------\n');

const q1only = external.filter(e => e.Q === 1);
const q1Result = runTest('Q=1 ONLY', q1only);

console.log('----------------------------------------------------------------------');
console.log('  STAGE 2c: THINGS-OVERLAP SUBSAMPLE');
console.log('----------------------------------------------------------------------\n');

const thingsNames = ['NGC2366', 'NGC2976', 'NGC4214', 'NGC6946', 'NGC7793', 'DDO154', 'IC2574'];
const thingsExt = external.filter(e => thingsNames.includes(e.name));
const thingsResult = runTest('THINGS GALAXIES (not in N45)', thingsExt);

console.log('----------------------------------------------------------------------');
console.log('  STAGE 3: DIAGNOSTIC CHECKS');
console.log('----------------------------------------------------------------------\n');

const obs = external.map(g => g.logA0_obs);
const predM3 = external.map(g => predictM3(g));
const residM3 = external.map((g, i) => obs[i] - predM3[i]);

console.log('  M3 residual diagnostics on external sample:');
console.log('    Mean residual: ' + mean(residM3).toFixed(4) + ' dex');
console.log('    SD residual: ' + sd(residM3).toFixed(4) + ' dex');
console.log('    Median residual: ' + [...residM3].sort((a, b) => a - b)[Math.floor(residM3.length / 2)].toFixed(4) + ' dex');
console.log();

const morphRanges = { earlySpiral: [0, 5], lateSpiral: [6, 8], irregular: [9, 11] };
for (const [label, [lo, hi]] of Object.entries(morphRanges)) {
  const sub = external.filter(g => g.T >= lo && g.T <= hi);
  if (sub.length < 2) continue;
  const subRes = sub.map(g => g.logA0_obs - predictM3(g));
  console.log('  ' + label + ' (T=' + lo + '-' + hi + ', N=' + sub.length + '): mean residual = ' + mean(subRes).toFixed(3) + ', SD = ' + sd(subRes).toFixed(3));
}
console.log();

console.log('----------------------------------------------------------------------');
console.log('  STAGE 4: SUMMARY VERDICT');
console.log('----------------------------------------------------------------------\n');

const tests = [];
if (fullResult) {
  tests.push({ name: 'M3 RMS < M0 RMS (full)', pass: fullResult.m3BeatsM0, detail: fullResult.rmsM3.toFixed(3) + ' vs ' + fullResult.rmsM0.toFixed(3) });
  tests.push({ name: 'M3 Spearman > 0.3 (full)', pass: fullResult.m3RankSig, detail: 'rs=' + fullResult.rsM3.toFixed(3) });
  tests.push({ name: 'M3 gap% > 0 (full)', pass: fullResult.gapM3 > 0, detail: fullResult.gapM3.toFixed(1) + '%' });
  tests.push({ name: 'M3 wins > 50% of galaxies', pass: fullResult.m3WinsCount > fullResult.n / 2, detail: fullResult.m3WinsCount + '/' + fullResult.n });
}

let allPass = true;
for (const t of tests) {
  console.log('  ' + (t.pass ? 'PASS' : 'FAIL') + '  ' + t.name.padEnd(35) + t.detail);
  if (!t.pass) allPass = false;
}
console.log();

if (allPass) {
  console.log('  ===================================================================');
  console.log('  VERDICT: M3 GENERALIZES BEYOND THE TRAINING SAMPLE');
  console.log('  The frozen M3 law predicts a0 for unseen SPARC galaxies better');
  console.log('  than the universal constant. The structured signal is real.');
  console.log('  ===================================================================');
} else {
  const partialPass = tests.filter(t => t.pass).length;
  if (partialPass >= 2) {
    console.log('  ===================================================================');
    console.log('  VERDICT: PARTIAL SUPPORT — M3 shows signal but with caveats');
    console.log('  ' + partialPass + '/' + tests.length + ' criteria met. See diagnostics above.');
    console.log('  ===================================================================');
  } else {
    console.log('  ===================================================================');
    console.log('  VERDICT: EXTERNAL TEST DID NOT CONFIRM M3');
    console.log('  The frozen law does not generalize to this external sample.');
    console.log('  This is a serious caution against over-interpreting the training');
    console.log('  result. See detailed diagnostics above.');
    console.log('  ===================================================================');
  }
}

const output = {
  phase: '80',
  title: 'External Clean-Sample Validation',
  protocol: {
    method: 'Frozen M3/M5 coefficients applied to SPARC galaxies outside N=45 training set',
    selectionCriteria: 'Published distance (fD>=2), n>=5 RC points, Vflat>0, logA0>=2.0, NOT in N=45',
    noRefit: true,
    frozenM3: frozenM3,
    frozenM5: frozenM5
  },
  fullSample: fullResult ? {
    n: fullResult.n,
    rmsM0: +fullResult.rmsM0.toFixed(4),
    rmsM3: +fullResult.rmsM3.toFixed(4),
    rmsM5: +fullResult.rmsM5.toFixed(4),
    gapM3: +fullResult.gapM3.toFixed(1),
    gapM5: +fullResult.gapM5.toFixed(1),
    spearmanM3: +fullResult.rsM3.toFixed(3),
    spearmanM5: +fullResult.rsM5.toFixed(3),
    m3BeatsM0: fullResult.m3BeatsM0,
    m3WinRate: +(fullResult.m3WinsCount / fullResult.n * 100).toFixed(1)
  } : null,
  q1Subsample: q1Result ? {
    n: q1Result.n,
    rmsM0: +q1Result.rmsM0.toFixed(4),
    rmsM3: +q1Result.rmsM3.toFixed(4),
    gapM3: +q1Result.gapM3.toFixed(1),
    spearmanM3: +q1Result.rsM3.toFixed(3),
    m3BeatsM0: q1Result.m3BeatsM0
  } : null,
  thingsSubsample: thingsResult ? {
    n: thingsResult.n,
    rmsM0: +thingsResult.rmsM0.toFixed(4),
    rmsM3: +thingsResult.rmsM3.toFixed(4),
    gapM3: +thingsResult.gapM3.toFixed(1),
    spearmanM3: +thingsResult.rsM3.toFixed(3),
    m3BeatsM0: thingsResult.m3BeatsM0
  } : null,
  perGalaxy: external.map(g => ({
    name: g.name,
    observed: g.logA0_obs,
    predictedM3: +predictM3(g).toFixed(4),
    predictedM5: +predictM5(g).toFixed(4),
    predictedM0: +trainMeanA0.toFixed(4),
    errorM3: +(predictM3(g) - g.logA0_obs).toFixed(4),
    Q: g.Q,
    nPts: g.nPts,
    D_Mpc: g.D,
    hasMasterData: g.hasMasterData
  })),
  verdict: allPass ? 'PASS' : (tests.filter(t => t.pass).length >= 2 ? 'PARTIAL' : 'FAIL')
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase80-external-validation.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase80-external-validation.json');
