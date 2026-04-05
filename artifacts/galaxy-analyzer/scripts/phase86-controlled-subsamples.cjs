#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 86: CONTROLLED SUBSAMPLES');
console.log('  Is N=45 just an extreme draw from Vflat-high regime?');
console.log('  Or does it need additional criteria beyond Vflat?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0(pts) {
  let lo = 0.5, hi = 5.0;
  for (let s = 0; s < 200; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gb = Math.pow(10, p.logGbar);
      c1 += (p.logGobs - Math.log10(mcgaughRAR(gb, Math.pow(10, m1)))) ** 2;
      c2 += (p.logGobs - Math.log10(mcgaughRAR(gb, Math.pow(10, m2)))) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  const logA0 = (lo + hi) / 2, a0 = Math.pow(10, logA0);
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.logGbar);
    const pred = mcgaughRAR(gb, a0);
    ss += (p.logGobs - Math.log10(pred > 0 ? pred : 1e-20)) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length) };
}

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function percentile(a, p) { const s = [...a].sort((x, y) => x - y); const i = Math.floor(s.length * p); return s[Math.min(i, s.length - 1)]; }

function pearsonR(x, y) {
  const n = x.length;
  if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

const csvPath = path.join(__dirname, '..', 'public', 'replication', 'N45_final_dataset.csv');
const csvLines = fs.readFileSync(csvPath, 'utf-8').trim().split('\n');
const headers = csvLines[0].split(',');
const trainCSV = csvLines.slice(1).map(l => {
  const vals = l.split(',');
  const obj = {};
  headers.forEach((h, i) => obj[h] = vals[i]);
  return obj;
});
const trainNames = new Set(trainCSV.map(d => d.galaxy_name));

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const D = parseFloat(line.substring(15, 21).trim());
  const fD = parseInt(line.substring(28, 29).trim());
  const inc = parseFloat(line.substring(30, 34).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { T, D, fD, inc, L36, Rdisk, MHI, Vflat, Q };
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
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul });
}

const p25 = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase25-group-membership.json'), 'utf-8'));
const envLookup = {};
for (const g of p25.galaxyAssignments) envLookup[g.name] = g;

const allGalaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0) continue;
  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    const gObs = pt.vobs * pt.vobs / pt.rad;
    const gBar = Math.abs(vBarSq) / pt.rad;
    if (gBar <= 0 || gObs <= 0 || !isFinite(gBar) || !isFinite(gObs)) continue;
    pts.push({ r: pt.rad, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar) });
  }
  if (pts.length < 3) continue;
  const fit = fitA0(pts);
  if (!isFinite(fit.logA0) || fit.logA0 < 0.5 || fit.logA0 > 5.5) continue;
  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const rarResid = sorted.map(p => {
    const gb = Math.pow(10, p.logGbar);
    const pred = mcgaughRAR(gb, fit.a0);
    return p.logGobs - Math.log10(pred > 0 ? pred : 1e-20);
  });
  let currentSign = Math.sign(rarResid[0]);
  let runLen = 1;
  const residRun = [];
  for (let i = 1; i < rarResid.length; i++) {
    const s = Math.sign(rarResid[i]);
    if (s === currentSign) runLen++;
    else { residRun.push(runLen); runLen = 1; currentSign = s; }
  }
  residRun.push(runLen);
  const meanRunLen = residRun.reduce((s, v) => s + v, 0) / residRun.length;
  const logMeanRun = Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1);
  const ev = envLookup[name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI: Math.log10(t1.MHI),
    logMeanRun,
    Vflat: t1.Vflat, Q: t1.Q, D: t1.D, fD: t1.fD, T: t1.T,
    L36: t1.L36, inc: t1.inc, envCode,
    inTraining: trainNames.has(name),
  });
}

const train45 = allGalaxies.filter(g => g.inTraining);
const TRAIN_R_MHI = pearsonR(train45.map(g => g.logMHI), train45.map(g => g.logA0));

console.log('  Total SPARC: ' + allGalaxies.length);
console.log('  Training N=45: r(MHI,a0) = ' + TRAIN_R_MHI.toFixed(3) + '\n');

const N_PERM = 5000;

function drawAndMeasure(pool, n, label) {
  const rs = [];
  for (let p = 0; p < N_PERM; p++) {
    const shuffled = [...pool].sort(() => Math.random() - 0.5);
    const sample = shuffled.slice(0, Math.min(n, shuffled.length));
    if (sample.length < 10) continue;
    const r = pearsonR(sample.map(g => g.logMHI), sample.map(g => g.logA0));
    if (!isNaN(r)) rs.push(r);
  }
  if (rs.length < 100) return null;
  const nBelow = rs.filter(r => r <= TRAIN_R_MHI).length;
  const pVal = nBelow / rs.length;
  return {
    label, poolSize: pool.length, drawSize: Math.min(n, pool.length),
    nDraws: rs.length,
    mean_r: +mean(rs).toFixed(3), sd_r: +sd(rs).toFixed(3),
    p5: +percentile(rs, 0.05).toFixed(3), p25: +percentile(rs, 0.25).toFixed(3),
    median_r: +percentile(rs, 0.50).toFixed(3),
    p75: +percentile(rs, 0.75).toFixed(3), p95: +percentile(rs, 0.95).toFixed(3),
    min_r: +Math.min(...rs).toFixed(3), max_r: +Math.max(...rs).toFixed(3),
    nBelowTraining: nBelow, pValue: +pVal.toFixed(4),
    training_r: +TRAIN_R_MHI.toFixed(3),
  };
}

function printResult(r) {
  if (!r) { console.log('    (insufficient pool size)'); return; }
  console.log('    Pool N=' + r.poolSize + ', draw N=' + r.drawSize + ', ' + r.nDraws + ' draws');
  console.log('    r(MHI,a0) distribution: mean=' + r.mean_r + ' ± ' + r.sd_r +
    '  [P5=' + r.p5 + ', median=' + r.median_r + ', P95=' + r.p95 + ']');
  console.log('    Range: [' + r.min_r + ', ' + r.max_r + ']');
  console.log('    P(r <= ' + r.training_r + '): ' + r.pValue + ' (' + r.nBelowTraining + '/' + r.nDraws + ')');
  const verdict = r.pValue > 0.10 ? 'TRAINING IS TYPICAL' :
    r.pValue > 0.01 ? 'TRAINING IS UNUSUAL' : 'TRAINING IS EXTREME';
  console.log('    → ' + verdict);
  console.log();
  return verdict;
}

console.log('======================================================================');
console.log('  TEST 1: Random N=45 from Vflat-only pools');
console.log('  Does the Vflat regime alone explain the training r?');
console.log('======================================================================\n');

const test1 = [];
const vflatPools = [
  { name: 'Vflat >= 60', pool: allGalaxies.filter(g => g.Vflat >= 60) },
  { name: 'Vflat >= 80', pool: allGalaxies.filter(g => g.Vflat >= 80) },
  { name: 'Vflat >= 100', pool: allGalaxies.filter(g => g.Vflat >= 100) },
  { name: 'Vflat >= 120', pool: allGalaxies.filter(g => g.Vflat >= 120) },
];

for (const vp of vflatPools) {
  console.log('  ── ' + vp.name + ' (pool N=' + vp.pool.length + ') ──');
  const r = drawAndMeasure(vp.pool, 45, vp.name);
  test1.push(r);
  printResult(r);
}

console.log('======================================================================');
console.log('  TEST 2: Vflat + ONE additional criterion');
console.log('  Which criterion on top of Vflat pushes r toward training?');
console.log('======================================================================\n');

const test2 = [];
const base = allGalaxies.filter(g => g.Vflat >= 80);
const addCriteria = [
  { name: 'Vflat>=80 only (baseline)', pool: base },
  { name: 'Vflat>=80 + Q<=2', pool: base.filter(g => g.Q <= 2) },
  { name: 'Vflat>=80 + Q==1', pool: base.filter(g => g.Q === 1) },
  { name: 'Vflat>=80 + nPts>=10', pool: base.filter(g => g.nPts >= 10) },
  { name: 'Vflat>=80 + nPts>=15', pool: base.filter(g => g.nPts >= 15) },
  { name: 'Vflat>=80 + nPts>=20', pool: base.filter(g => g.nPts >= 20) },
  { name: 'Vflat>=80 + fD>=2', pool: base.filter(g => g.fD >= 2) },
  { name: 'Vflat>=80 + T<=7', pool: base.filter(g => g.T <= 7) },
  { name: 'Vflat>=80 + T<=5', pool: base.filter(g => g.T <= 5) },
  { name: 'Vflat>=80 + env>=1', pool: base.filter(g => g.envCode >= 1) },
  { name: 'Vflat>=80 + D>=10', pool: base.filter(g => g.D >= 10) },
  { name: 'Vflat>=80 + L36>=5', pool: base.filter(g => g.L36 >= 5) },
];

for (const ac of addCriteria) {
  console.log('  ── ' + ac.name + ' (pool N=' + ac.pool.length + ') ──');
  const r = drawAndMeasure(ac.pool, 45, ac.name);
  test2.push(r);
  printResult(r);
}

console.log('======================================================================');
console.log('  TEST 3: Progressive matching toward N=45 profile');
console.log('  Match on multiple criteria simultaneously');
console.log('======================================================================\n');

const test3 = [];
const progressivePools = [
  { name: 'Vflat>=80', pool: allGalaxies.filter(g => g.Vflat >= 80) },
  { name: 'Vflat>=80 + nPts>=10', pool: allGalaxies.filter(g => g.Vflat >= 80 && g.nPts >= 10) },
  { name: 'Vflat>=80 + nPts>=10 + Q<=2', pool: allGalaxies.filter(g => g.Vflat >= 80 && g.nPts >= 10 && g.Q <= 2) },
  { name: 'Vflat>=80 + nPts>=10 + Q<=2 + T<=8', pool: allGalaxies.filter(g => g.Vflat >= 80 && g.nPts >= 10 && g.Q <= 2 && g.T <= 8) },
  { name: 'Vflat>=80 + nPts>=10 + Q<=2 + T<=8 + fD>=2', pool: allGalaxies.filter(g => g.Vflat >= 80 && g.nPts >= 10 && g.Q <= 2 && g.T <= 8 && g.fD >= 2) },
];

for (const pp of progressivePools) {
  console.log('  ── ' + pp.name + ' (pool N=' + pp.pool.length + ') ──');
  const r = drawAndMeasure(pp.pool, Math.min(45, pp.pool.length), pp.name);
  test3.push(r);
  printResult(r);
}

console.log('======================================================================');
console.log('  TEST 4: Draw from Vflat>=80, EXCLUDING training galaxies');
console.log('  Does the regime effect persist without the N=45 themselves?');
console.log('======================================================================\n');

const test4 = [];
const extOnly = allGalaxies.filter(g => g.Vflat >= 80 && !g.inTraining);
console.log('  ── Vflat>=80 EXCLUDING N=45 (pool N=' + extOnly.length + ') ──');
const r4 = drawAndMeasure(extOnly, 45, 'Vflat>=80 excl. training');
test4.push(r4);
printResult(r4);

const extOnlyStrict = allGalaxies.filter(g => g.Vflat >= 80 && g.nPts >= 10 && !g.inTraining);
console.log('  ── Vflat>=80 + nPts>=10 EXCLUDING N=45 (pool N=' + extOnlyStrict.length + ') ──');
const r4b = drawAndMeasure(extOnlyStrict, Math.min(45, extOnlyStrict.length), 'Vflat>=80+nPts>=10 excl.');
test4.push(r4b);
printResult(r4b);

console.log('======================================================================');
console.log('  TEST 5: Actual r(MHI,a0) of the Vflat>=80 pool itself');
console.log('  (not random draws, but the deterministic population value)');
console.log('======================================================================\n');

const pools = [
  { name: 'All SPARC', data: allGalaxies },
  { name: 'Vflat>=60', data: allGalaxies.filter(g => g.Vflat >= 60) },
  { name: 'Vflat>=80', data: allGalaxies.filter(g => g.Vflat >= 80) },
  { name: 'Vflat>=80, excl N=45', data: allGalaxies.filter(g => g.Vflat >= 80 && !g.inTraining) },
  { name: 'Vflat>=80+nPts>=10', data: allGalaxies.filter(g => g.Vflat >= 80 && g.nPts >= 10) },
  { name: 'Vflat>=80+nPts>=10, excl N=45', data: allGalaxies.filter(g => g.Vflat >= 80 && g.nPts >= 10 && !g.inTraining) },
  { name: 'Training N=45', data: train45 },
];

console.log('  ' + 'Pool'.padEnd(40) + 'N'.padStart(4) + 'r(MHI)'.padStart(9) + 'r(Run)'.padStart(9));
for (const p of pools) {
  if (p.data.length < 5) continue;
  const r_mhi = pearsonR(p.data.map(g => g.logMHI), p.data.map(g => g.logA0));
  const r_run = pearsonR(p.data.map(g => g.logMeanRun), p.data.map(g => g.logA0));
  console.log('  ' + p.name.padEnd(40) + String(p.data.length).padStart(4) +
    r_mhi.toFixed(3).padStart(9) + r_run.toFixed(3).padStart(9));
}
console.log();

console.log('======================================================================');
console.log('  FINAL SUMMARY');
console.log('======================================================================\n');

console.log('  Training r(MHI,a0) = ' + TRAIN_R_MHI.toFixed(3));
console.log();

const allResults = [...test1, ...test2, ...test3, ...test4].filter(r => r);
const typical = allResults.filter(r => r.pValue > 0.10);
const unusual = allResults.filter(r => r.pValue <= 0.10 && r.pValue > 0.01);
const extreme = allResults.filter(r => r.pValue <= 0.01);

console.log('  Pools where training r is TYPICAL (p>0.10): ' + typical.length);
for (const r of typical) console.log('    ' + r.label + ': p=' + r.pValue + ', mean_r=' + r.mean_r);
console.log();
console.log('  Pools where training r is UNUSUAL (0.01<p<=0.10): ' + unusual.length);
for (const r of unusual) console.log('    ' + r.label + ': p=' + r.pValue + ', mean_r=' + r.mean_r);
console.log();
console.log('  Pools where training r is EXTREME (p<=0.01): ' + extreme.length);
for (const r of extreme) console.log('    ' + r.label + ': p=' + r.pValue + ', mean_r=' + r.mean_r);
console.log();

if (typical.length > 0) {
  console.log('  ╔════════════════════════════════════════════════════════════════╗');
  console.log('  ║  SOME POOLS MAKE TRAINING r TYPICAL → regime effect real     ║');
  console.log('  ╚════════════════════════════════════════════════════════════════╝');
  console.log('  The pools that normalize the training r:');
  for (const r of typical) console.log('    ' + r.label);
} else {
  console.log('  ╔════════════════════════════════════════════════════════════════╗');
  console.log('  ║  TRAINING r REMAINS EXTREME IN ALL POOLS → pure artifact    ║');
  console.log('  ╚════════════════════════════════════════════════════════════════╝');
}
console.log();

const output = {
  phase: '86',
  title: 'Controlled Subsamples — Is N=45 extreme or typical within Vflat-high?',
  training_r_MHI: +TRAIN_R_MHI.toFixed(3),
  nPermutations: N_PERM,
  test1_vflatOnly: test1.filter(r => r),
  test2_vflatPlus: test2.filter(r => r),
  test3_progressive: test3.filter(r => r),
  test4_excludeTraining: test4.filter(r => r),
  nTypical: typical.length,
  nUnusual: unusual.length,
  nExtreme: extreme.length,
  typicalPools: typical.map(r => r.label),
  extremePools: extreme.map(r => r.label),
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase86-controlled-subsamples.json'), JSON.stringify(output, null, 2));
console.log('  Saved: public/phase86-controlled-subsamples.json');
