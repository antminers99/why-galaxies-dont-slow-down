#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 85: CRITERION-BY-CRITERION DISSECTION');
console.log('  Which selection criterion flips r(MHI, a0) from + to −?');
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
  const hasMhost = trainNames.has(name);

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI: Math.log10(t1.MHI),
    logMeanRun,
    Vflat: t1.Vflat, Q: t1.Q, D: t1.D, fD: t1.fD, T: t1.T,
    L36: t1.L36, Rdisk: t1.Rdisk, inc: t1.inc,
    envCode, hasMhost,
    inTraining: trainNames.has(name),
  });
}

console.log('  Total SPARC galaxies: ' + allGalaxies.length + '\n');

function reportCorrs(label, subset) {
  if (subset.length < 5) {
    console.log('  ' + label.padEnd(50) + 'N=' + String(subset.length).padStart(3) + '  (too few)');
    return { label, n: subset.length, r_MHI: NaN, r_Run: NaN, r_Vflat: NaN };
  }
  const r_mhi = pearsonR(subset.map(g => g.logMHI), subset.map(g => g.logA0));
  const r_run = pearsonR(subset.map(g => g.logMeanRun), subset.map(g => g.logA0));
  const r_vf = pearsonR(subset.map(g => Math.log10(g.Vflat > 0 ? g.Vflat : 50)), subset.map(g => g.logA0));
  const sign_mhi = r_mhi < -0.05 ? '(−)' : r_mhi > 0.05 ? '(+)' : '(~0)';
  console.log('  ' + label.padEnd(50) + 'N=' + String(subset.length).padStart(3) +
    '  r(MHI)=' + r_mhi.toFixed(3).padStart(7) + ' ' + sign_mhi +
    '  r(Run)=' + r_run.toFixed(3).padStart(7) +
    '  r(Vf)=' + r_vf.toFixed(3).padStart(7));
  return { label, n: subset.length, r_MHI: +r_mhi.toFixed(3), r_Run: +r_run.toFixed(3), r_Vflat: +r_vf.toFixed(3) };
}

console.log('======================================================================');
console.log('  TEST 1: Single criterion — does ONE criterion flip the sign?');
console.log('======================================================================\n');

const test1Results = [];
test1Results.push(reportCorrs('ALL SPARC (no cuts)', allGalaxies));
test1Results.push(reportCorrs('Training N=45', allGalaxies.filter(g => g.inTraining)));
console.log();

const criteria = [
  { name: 'fD >= 2 (published distance)', f: g => g.fD >= 2 },
  { name: 'fD == 1 (Hubble flow ONLY)', f: g => g.fD === 1 },
  { name: 'fD == 2 (TRGB only)', f: g => g.fD === 2 },
  { name: 'fD == 4 (Ursa Major only)', f: g => g.fD === 4 },
  { name: 'Q == 1 (best quality)', f: g => g.Q === 1 },
  { name: 'Q <= 2', f: g => g.Q <= 2 },
  { name: 'nPts >= 10', f: g => g.nPts >= 10 },
  { name: 'nPts >= 15', f: g => g.nPts >= 15 },
  { name: 'nPts >= 20', f: g => g.nPts >= 20 },
  { name: 'Vflat >= 60 km/s', f: g => g.Vflat >= 60 },
  { name: 'Vflat >= 80 km/s', f: g => g.Vflat >= 80 },
  { name: 'Vflat >= 100 km/s', f: g => g.Vflat >= 100 },
  { name: 'Vflat >= 120 km/s', f: g => g.Vflat >= 120 },
  { name: 'T <= 5 (Sa-Sc)', f: g => g.T <= 5 },
  { name: 'T <= 7 (Sa-Sd)', f: g => g.T <= 7 },
  { name: 'T <= 8 (not Sm/Im)', f: g => g.T <= 8 },
  { name: 'T >= 5 (Sc-Im late types)', f: g => g.T >= 5 },
  { name: 'envCode >= 1 (group/cluster)', f: g => g.envCode >= 1 },
  { name: 'envCode == 0 (field only)', f: g => g.envCode === 0 },
  { name: 'inc >= 50 deg', f: g => g.inc >= 50 },
  { name: 'inc >= 60 deg', f: g => g.inc >= 60 },
  { name: 'D >= 5 Mpc', f: g => g.D >= 5 },
  { name: 'D >= 10 Mpc', f: g => g.D >= 10 },
  { name: 'L36 >= 1 (1e9 Lsun)', f: g => g.L36 >= 1 },
  { name: 'L36 >= 5 (1e9 Lsun)', f: g => g.L36 >= 5 },
];

for (const c of criteria) {
  test1Results.push(reportCorrs(c.name, allGalaxies.filter(c.f)));
}
console.log();

const negMHI = test1Results.filter(r => r.r_MHI < -0.05 && r.n >= 5);
console.log('  Criteria that produce NEGATIVE r(MHI, a0):');
for (const r of negMHI) console.log('    ' + r.label + ': r=' + r.r_MHI);
console.log();

console.log('======================================================================');
console.log('  TEST 2: Two-criterion combinations — find the FLIP POINT');
console.log('======================================================================\n');

const test2Results = [];

const comboCriteria = [
  { name: 'fD>=2', f: g => g.fD >= 2 },
  { name: 'Q<=2', f: g => g.Q <= 2 },
  { name: 'nPts>=10', f: g => g.nPts >= 10 },
  { name: 'Vflat>=80', f: g => g.Vflat >= 80 },
  { name: 'T<=7', f: g => g.T <= 7 },
  { name: 'T<=8', f: g => g.T <= 8 },
  { name: 'env>=1', f: g => g.envCode >= 1 },
  { name: 'Vflat>=100', f: g => g.Vflat >= 100 },
  { name: 'nPts>=15', f: g => g.nPts >= 15 },
  { name: 'L36>=1', f: g => g.L36 >= 1 },
  { name: 'D>=10', f: g => g.D >= 10 },
];

for (let i = 0; i < comboCriteria.length; i++) {
  for (let j = i + 1; j < comboCriteria.length; j++) {
    const c1 = comboCriteria[i], c2 = comboCriteria[j];
    const subset = allGalaxies.filter(g => c1.f(g) && c2.f(g));
    const label = c1.name + ' + ' + c2.name;
    test2Results.push(reportCorrs(label, subset));
  }
}
console.log();

const negCombo = test2Results.filter(r => r.r_MHI < -0.10 && r.n >= 10);
console.log('  Two-criterion combos with r(MHI) < −0.10 and N >= 10:');
for (const r of negCombo.sort((a, b) => a.r_MHI - b.r_MHI)) {
  console.log('    ' + r.label + ': N=' + r.n + ', r=' + r.r_MHI);
}
console.log();

console.log('======================================================================');
console.log('  TEST 3: Cumulative — add criteria ONE AT A TIME');
console.log('  Track exactly where the sign flips');
console.log('======================================================================\n');

const cumulPaths = [
  {
    name: 'Path A: Vflat first',
    steps: [
      { label: 'Start: all', f: [] },
      { label: '+ Vflat>=80', f: [g => g.Vflat >= 80] },
      { label: '+ nPts>=10', f: [g => g.Vflat >= 80, g => g.nPts >= 10] },
      { label: '+ fD>=2', f: [g => g.Vflat >= 80, g => g.nPts >= 10, g => g.fD >= 2] },
      { label: '+ Q<=2', f: [g => g.Vflat >= 80, g => g.nPts >= 10, g => g.fD >= 2, g => g.Q <= 2] },
      { label: '+ T<=8', f: [g => g.Vflat >= 80, g => g.nPts >= 10, g => g.fD >= 2, g => g.Q <= 2, g => g.T <= 8] },
    ]
  },
  {
    name: 'Path B: fD first',
    steps: [
      { label: 'Start: all', f: [] },
      { label: '+ fD>=2', f: [g => g.fD >= 2] },
      { label: '+ Vflat>=80', f: [g => g.fD >= 2, g => g.Vflat >= 80] },
      { label: '+ nPts>=10', f: [g => g.fD >= 2, g => g.Vflat >= 80, g => g.nPts >= 10] },
      { label: '+ Q<=2', f: [g => g.fD >= 2, g => g.Vflat >= 80, g => g.nPts >= 10, g => g.Q <= 2] },
      { label: '+ T<=8', f: [g => g.fD >= 2, g => g.Vflat >= 80, g => g.nPts >= 10, g => g.Q <= 2, g => g.T <= 8] },
    ]
  },
  {
    name: 'Path C: nPts first',
    steps: [
      { label: 'Start: all', f: [] },
      { label: '+ nPts>=10', f: [g => g.nPts >= 10] },
      { label: '+ Vflat>=80', f: [g => g.nPts >= 10, g => g.Vflat >= 80] },
      { label: '+ fD>=2', f: [g => g.nPts >= 10, g => g.Vflat >= 80, g => g.fD >= 2] },
      { label: '+ Q<=2', f: [g => g.nPts >= 10, g => g.Vflat >= 80, g => g.fD >= 2, g => g.Q <= 2] },
      { label: '+ T<=8', f: [g => g.nPts >= 10, g => g.Vflat >= 80, g => g.fD >= 2, g => g.Q <= 2, g => g.T <= 8] },
    ]
  },
  {
    name: 'Path D: Luminosity path',
    steps: [
      { label: 'Start: all', f: [] },
      { label: '+ L36>=1', f: [g => g.L36 >= 1] },
      { label: '+ Vflat>=80', f: [g => g.L36 >= 1, g => g.Vflat >= 80] },
      { label: '+ nPts>=10', f: [g => g.L36 >= 1, g => g.Vflat >= 80, g => g.nPts >= 10] },
      { label: '+ fD>=2', f: [g => g.L36 >= 1, g => g.Vflat >= 80, g => g.nPts >= 10, g => g.fD >= 2] },
    ]
  },
];

const cumulResults = [];
for (const pathDef of cumulPaths) {
  console.log('  ── ' + pathDef.name + ' ──');
  for (const step of pathDef.steps) {
    let subset = [...allGalaxies];
    for (const fi of step.f) subset = subset.filter(fi);
    const r = reportCorrs(step.label, subset);
    cumulResults.push({ path: pathDef.name, ...r });
  }
  console.log();
}

console.log('======================================================================');
console.log('  TEST 4: The Vflat threshold scan');
console.log('  Scan Vflat threshold from 40 to 200 km/s');
console.log('======================================================================\n');

const vflatScan = [];
console.log('  Vflat_min    N   r(MHI,a0)  r(Run,a0)');
for (let vMin = 40; vMin <= 200; vMin += 10) {
  const subset = allGalaxies.filter(g => g.Vflat >= vMin);
  if (subset.length < 10) break;
  const r_mhi = pearsonR(subset.map(g => g.logMHI), subset.map(g => g.logA0));
  const r_run = pearsonR(subset.map(g => g.logMeanRun), subset.map(g => g.logA0));
  const marker = r_mhi < -0.05 ? ' ←FLIP' : '';
  console.log('  ' + String(vMin).padStart(4) + '       ' + String(subset.length).padStart(3) + '   ' +
    r_mhi.toFixed(3).padStart(7) + '    ' + r_run.toFixed(3).padStart(7) + marker);
  vflatScan.push({ vMin, n: subset.length, r_MHI: +r_mhi.toFixed(3), r_Run: +r_run.toFixed(3) });
}
console.log();

console.log('======================================================================');
console.log('  TEST 5: nPts threshold scan');
console.log('======================================================================\n');

const nptsScan = [];
console.log('  nPts_min    N   r(MHI,a0)  r(Run,a0)');
for (let nMin = 3; nMin <= 50; nMin += 3) {
  const subset = allGalaxies.filter(g => g.nPts >= nMin);
  if (subset.length < 10) break;
  const r_mhi = pearsonR(subset.map(g => g.logMHI), subset.map(g => g.logA0));
  const r_run = pearsonR(subset.map(g => g.logMeanRun), subset.map(g => g.logA0));
  const marker = r_mhi < -0.05 ? ' ←FLIP' : '';
  console.log('  ' + String(nMin).padStart(4) + '       ' + String(subset.length).padStart(3) + '   ' +
    r_mhi.toFixed(3).padStart(7) + '    ' + r_run.toFixed(3).padStart(7) + marker);
  nptsScan.push({ nMin, n: subset.length, r_MHI: +r_mhi.toFixed(3), r_Run: +r_run.toFixed(3) });
}
console.log();

console.log('======================================================================');
console.log('  TEST 6: Combined Vflat + nPts scan (the two key filters)');
console.log('======================================================================\n');

const comboScan = [];
console.log('  Vflat_min  nPts_min     N   r(MHI,a0)  r(Run,a0)');
for (const vMin of [50, 60, 70, 80, 90, 100, 120]) {
  for (const nMin of [5, 8, 10, 12, 15, 20]) {
    const subset = allGalaxies.filter(g => g.Vflat >= vMin && g.nPts >= nMin);
    if (subset.length < 10) continue;
    const r_mhi = pearsonR(subset.map(g => g.logMHI), subset.map(g => g.logA0));
    const r_run = pearsonR(subset.map(g => g.logMeanRun), subset.map(g => g.logA0));
    const marker = r_mhi < -0.10 ? ' ←FLIP' : '';
    console.log('  ' + String(vMin).padStart(4) + '         ' + String(nMin).padStart(3) + '       ' +
      String(subset.length).padStart(3) + '   ' + r_mhi.toFixed(3).padStart(7) + '    ' + r_run.toFixed(3).padStart(7) + marker);
    comboScan.push({ vMin, nMin, n: subset.length, r_MHI: +r_mhi.toFixed(3), r_Run: +r_run.toFixed(3) });
  }
}
console.log();

console.log('======================================================================');
console.log('  SUMMARY: Where does the sign flip?');
console.log('======================================================================\n');

const singleFlips = test1Results.filter(r => r.r_MHI < -0.05 && r.n >= 5);
console.log('  Single criteria that produce negative r(MHI):');
for (const r of singleFlips) console.log('    ' + r.label + ': N=' + r.n + ', r=' + r.r_MHI);
console.log();

const flipThreshold_vflat = vflatScan.find(v => v.r_MHI < -0.05);
const flipThreshold_npts = nptsScan.find(v => v.r_MHI < -0.05);
console.log('  Vflat threshold for flip: ' + (flipThreshold_vflat ? flipThreshold_vflat.vMin + ' km/s (N=' + flipThreshold_vflat.n + ')' : 'never flips'));
console.log('  nPts threshold for flip: ' + (flipThreshold_npts ? flipThreshold_npts.nMin + ' (N=' + flipThreshold_npts.n + ')' : 'never flips'));
console.log();

const comboFlips = comboScan.filter(v => v.r_MHI < -0.10);
if (comboFlips.length > 0) {
  console.log('  Combined Vflat+nPts that flip r < −0.10:');
  for (const c of comboFlips) console.log('    Vflat>=' + c.vMin + ' + nPts>=' + c.nMin + ': N=' + c.n + ', r=' + c.r_MHI);
} else {
  console.log('  No Vflat+nPts combo alone produces r < −0.10');
  console.log('  This means the flip requires MORE than just mass+data cuts');
}
console.log();

const output = {
  phase: '85',
  title: 'Criterion-by-Criterion Dissection',
  nTotal: allGalaxies.length,
  singleCriteria: test1Results,
  twoCriteriaCombos: test2Results.filter(r => !isNaN(r.r_MHI)),
  cumulativePaths: cumulResults,
  vflatScan,
  nptsScan,
  comboScan,
  singleFlips: singleFlips.map(r => r.label),
  flipVflat: flipThreshold_vflat ? flipThreshold_vflat.vMin : null,
  flipNpts: flipThreshold_npts ? flipThreshold_npts.nMin : null,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase85-criterion-dissection.json'), JSON.stringify(output, null, 2));
console.log('  Saved: public/phase85-criterion-dissection.json');
