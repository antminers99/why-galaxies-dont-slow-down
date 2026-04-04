#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 82: SAMPLE COMPARISON — Training N=45 vs External');
console.log('  Are the two populations similar or fundamentally different?');
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
    ss += (p.logGobs - Math.log10(mcgaughRAR(gb, a0))) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length) };
}

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function median(a) { const s = [...a].sort((x, y) => x - y); return s.length % 2 ? s[(s.length - 1) / 2] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function iqr(a) { const s = [...a].sort((x, y) => x - y); const q1 = s[Math.floor(s.length * 0.25)]; const q3 = s[Math.floor(s.length * 0.75)]; return [q1, q3]; }

function welchT(a, b) {
  const na = a.length, nb = b.length;
  const ma = mean(a), mb = mean(b);
  const va = a.reduce((s, v) => s + (v - ma) ** 2, 0) / (na - 1);
  const vb = b.reduce((s, v) => s + (v - mb) ** 2, 0) / (nb - 1);
  const se = Math.sqrt(va / na + vb / nb);
  if (se < 1e-12) return { t: 0, df: 2, p: 1.0, d: 0 };
  const t = (ma - mb) / se;
  const num = (va / na + vb / nb) ** 2;
  const den = (va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1);
  const df = num / den;
  const p = 2 * tCDF(-Math.abs(t), df);
  const pooledSD = Math.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2));
  const d = pooledSD > 0 ? (ma - mb) / pooledSD : 0;
  return { t, df, p, d };
}

function tCDF(t, df) {
  const x = df / (df + t * t);
  return 0.5 * betaInc(x, df / 2, 0.5);
}

function betaInc(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;
  const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log(1 - x));
  if (x < (a + 1) / (a + b + 2)) return bt * betaCF(x, a, b) / a;
  return 1 - bt * betaCF(1 - x, b, a) / b;
}

function betaCF(x, a, b) {
  const maxIter = 200;
  let qab = a + b, qap = a + 1, qam = a - 1;
  let c = 1, d = 1 - qab * x / qap;
  if (Math.abs(d) < 1e-30) d = 1e-30;
  d = 1 / d;
  let h = d;
  for (let m = 1; m <= maxIter; m++) {
    let m2 = 2 * m;
    let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1 / d;
    c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
    h *= d * c;
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1 + aa * d; if (Math.abs(d) < 1e-30) d = 1e-30; d = 1 / d;
    c = 1 + aa / c; if (Math.abs(c) < 1e-30) c = 1e-30;
    let del = d * c;
    h *= del;
    if (Math.abs(del - 1) < 3e-7) break;
  }
  return h;
}

function lgamma(x) {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let y = x, tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) ser += c[j] / ++y;
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

function mannWhitneyU(a, b) {
  const na = a.length, nb = b.length;
  const all = [...a.map(v => ({ v, g: 'a' })), ...b.map(v => ({ v, g: 'b' }))];
  all.sort((x, y) => x.v - y.v);
  let rank = 1;
  for (let i = 0; i < all.length;) {
    let j = i;
    while (j < all.length && all[j].v === all[i].v) j++;
    const avgRank = (rank + rank + j - i - 1) / 2;
    for (let k = i; k < j; k++) all[k].rank = avgRank;
    rank += j - i;
    i = j;
  }
  const Ra = all.filter(x => x.g === 'a').reduce((s, x) => s + x.rank, 0);
  const U = Ra - na * (na + 1) / 2;
  const mu = na * nb / 2;
  const sigma = Math.sqrt(na * nb * (na + nb + 1) / 12);
  const z = sigma > 0 ? (U - mu) / sigma : 0;
  const p = 2 * (1 - normalCDF(Math.abs(z)));
  return { U, z, p };
}

function normalCDF(x) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429;
  const p = 0.3275911;
  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.sqrt(2);
  const t = 1 / (1 + p * x);
  const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);
  return 0.5 * (1 + sign * y);
}

function ksTest(a, b) {
  const na = a.length, nb = b.length;
  const sa = [...a].sort((x, y) => x - y);
  const sb = [...b].sort((x, y) => x - y);
  let ia = 0, ib = 0, d = 0;
  while (ia < na && ib < nb) {
    const fa = (ia + 1) / na, fb = (ib + 1) / nb;
    if (sa[ia] <= sb[ib]) { d = Math.max(d, Math.abs(fa - ib / nb)); ia++; }
    else { d = Math.max(d, Math.abs(ia / na - fb)); ib++; }
  }
  while (ia < na) { d = Math.max(d, Math.abs((ia + 1) / na - 1)); ia++; }
  while (ib < nb) { d = Math.max(d, Math.abs(1 - (ib + 1) / nb)); ib++; }
  const ne = Math.sqrt(na * nb / (na + nb));
  const lambda = (ne + 0.12 + 0.11 / ne) * d;
  let p = 0;
  for (let k = 1; k <= 100; k++) p += 2 * (k % 2 ? 1 : -1) * Math.exp(-2 * k * k * lambda * lambda);
  p = Math.max(0, Math.min(1, 1 - p));
  return { D: d, p };
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

const sparcAll = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf-8'));
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

  const rMax = Math.max(...pts.map(p => p.r));

  const ev = envLookup[name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0,
    logMHI: Math.log10(t1.MHI),
    logMeanRun,
    Q: t1.Q, D: t1.D, fD: t1.fD, inc: t1.inc, T: t1.T,
    Vflat: t1.Vflat, L36: t1.L36, Rdisk: t1.Rdisk,
    rMax, envCode,
    inTraining: trainNames.has(name),
  });
}

const train = allGalaxies.filter(g => g.inTraining);
const ext = allGalaxies.filter(g => !g.inTraining && g.fD >= 2 && g.logA0 >= 2.0 && g.nPts >= 5);

console.log('  Training: N = ' + train.length);
console.log('  External (fD>=2, logA0>=2, nPts>=5): N = ' + ext.length);
console.log();

const variables = [
  { name: 'logMHI (gas mass)', key: 'logMHI', unit: 'dex' },
  { name: 'logMeanRun (kinematic coherence)', key: 'logMeanRun', unit: 'dex' },
  { name: 'logA0 (acceleration scale)', key: 'logA0', unit: 'dex' },
  { name: 'nPts (RC point count)', key: 'nPts', unit: '' },
  { name: 'Quality flag (Q)', key: 'Q', unit: '' },
  { name: 'Distance (D)', key: 'D', unit: 'Mpc' },
  { name: 'Inclination (inc)', key: 'inc', unit: 'deg' },
  { name: 'RC radial extent (rMax)', key: 'rMax', unit: 'kpc' },
  { name: 'Morphological type (T)', key: 'T', unit: '' },
  { name: 'Vflat', key: 'Vflat', unit: 'km/s' },
  { name: 'Luminosity L3.6', key: 'L36', unit: '1e9 Lsun' },
  { name: 'Disk scale length (Rdisk)', key: 'Rdisk', unit: 'kpc' },
  { name: 'Environment code', key: 'envCode', unit: '' },
];

console.log('======================================================================');
console.log('  VARIABLE-BY-VARIABLE COMPARISON');
console.log('======================================================================\n');

const results = [];

for (const v of variables) {
  const tVals = train.map(g => g[v.key]).filter(x => isFinite(x) && x !== null && x !== undefined);
  const eVals = ext.map(g => g[v.key]).filter(x => isFinite(x) && x !== null && x !== undefined);

  if (tVals.length < 3 || eVals.length < 3) {
    console.log('  ' + v.name + ': insufficient data\n');
    continue;
  }

  const tMean = mean(tVals), eMean = mean(eVals);
  const tMed = median(tVals), eMed = median(eVals);
  const tSD = sd(tVals), eSD = sd(eVals);
  const tIQR = iqr(tVals), eIQR = iqr(eVals);

  const welch = welchT(tVals, eVals);
  const mw = mannWhitneyU(tVals, eVals);
  const ks = ksTest(tVals, eVals);

  const sig = welch.p < 0.01 ? '***' : welch.p < 0.05 ? '**' : welch.p < 0.10 ? '*' : '';
  const effectSize = Math.abs(welch.d) > 0.8 ? 'LARGE' : Math.abs(welch.d) > 0.5 ? 'MEDIUM' : Math.abs(welch.d) > 0.2 ? 'SMALL' : 'NEGLIGIBLE';

  console.log('  ── ' + v.name + ' ' + sig + ' ──');
  console.log('               ' + 'Training'.padStart(12) + 'External'.padStart(12) + 'Diff'.padStart(10));
  console.log('    Mean       ' + tMean.toFixed(3).padStart(12) + eMean.toFixed(3).padStart(12) + (eMean - tMean >= 0 ? '+' : '') + (eMean - tMean).toFixed(3).padStart(9));
  console.log('    Median     ' + tMed.toFixed(3).padStart(12) + eMed.toFixed(3).padStart(12));
  console.log('    SD         ' + tSD.toFixed(3).padStart(12) + eSD.toFixed(3).padStart(12));
  console.log('    IQR        ' + ('[' + tIQR[0].toFixed(2) + ', ' + tIQR[1].toFixed(2) + ']').padStart(12) + ('[' + eIQR[0].toFixed(2) + ', ' + eIQR[1].toFixed(2) + ']').padStart(12));
  console.log('    Range      ' + ('[' + Math.min(...tVals).toFixed(2) + ', ' + Math.max(...tVals).toFixed(2) + ']').padStart(12) + ('[' + Math.min(...eVals).toFixed(2) + ', ' + Math.max(...eVals).toFixed(2) + ']').padStart(12));
  console.log('    Welch t    ' + welch.t.toFixed(2) + ', p=' + (welch.p < 0.001 ? '<0.001' : welch.p.toFixed(3)) + ', Cohen d=' + welch.d.toFixed(2) + ' (' + effectSize + ')');
  console.log('    Mann-Whit  ' + 'z=' + mw.z.toFixed(2) + ', p=' + (mw.p < 0.001 ? '<0.001' : mw.p.toFixed(3)));
  console.log('    KS test    ' + 'D=' + ks.D.toFixed(3) + ', p=' + (ks.p < 0.001 ? '<0.001' : ks.p.toFixed(3)));
  console.log();

  results.push({
    variable: v.name,
    key: v.key,
    training: { n: tVals.length, mean: +tMean.toFixed(3), median: +tMed.toFixed(3), sd: +tSD.toFixed(3), min: +Math.min(...tVals).toFixed(3), max: +Math.max(...tVals).toFixed(3) },
    external: { n: eVals.length, mean: +eMean.toFixed(3), median: +eMed.toFixed(3), sd: +eSD.toFixed(3), min: +Math.min(...eVals).toFixed(3), max: +Math.max(...eVals).toFixed(3) },
    welch: { t: +welch.t.toFixed(2), p: +welch.p.toFixed(4), d: +welch.d.toFixed(2) },
    mannWhitney: { z: +mw.z.toFixed(2), p: +mw.p.toFixed(4) },
    ks: { D: +ks.D.toFixed(3), p: +ks.p.toFixed(4) },
    effectSize,
    significant: welch.p < 0.05,
  });
}

console.log('======================================================================');
console.log('  QUALITY DISTRIBUTION COMPARISON');
console.log('======================================================================\n');

const qCounts = (arr) => {
  const c = { 1: 0, 2: 0, 3: 0 };
  arr.forEach(g => { if (c[g.Q] !== undefined) c[g.Q]++; });
  return c;
};
const tQ = qCounts(train), eQ = qCounts(ext);
console.log('  Quality  Training  External');
console.log('  Q=1      ' + String(tQ[1]).padStart(4) + ' (' + (tQ[1] / train.length * 100).toFixed(0) + '%)' + String(eQ[1]).padStart(6) + ' (' + (eQ[1] / ext.length * 100).toFixed(0) + '%)');
console.log('  Q=2      ' + String(tQ[2]).padStart(4) + ' (' + (tQ[2] / train.length * 100).toFixed(0) + '%)' + String(eQ[2]).padStart(6) + ' (' + (eQ[2] / ext.length * 100).toFixed(0) + '%)');
console.log('  Q=3      ' + String(tQ[3]).padStart(4) + ' (' + (tQ[3] / train.length * 100).toFixed(0) + '%)' + String(eQ[3]).padStart(6) + ' (' + (eQ[3] / ext.length * 100).toFixed(0) + '%)');
console.log();

console.log('======================================================================');
console.log('  ENVIRONMENT DISTRIBUTION COMPARISON');
console.log('======================================================================\n');

const envCounts = (arr) => {
  const c = { 0: 0, 1: 0, 2: 0 };
  arr.forEach(g => { if (c[g.envCode] !== undefined) c[g.envCode]++; });
  return c;
};
const tE = envCounts(train), eE = envCounts(ext);
console.log('  Env      Training  External');
console.log('  Field    ' + String(tE[0]).padStart(4) + ' (' + (tE[0] / train.length * 100).toFixed(0) + '%)' + String(eE[0]).padStart(6) + ' (' + (eE[0] / ext.length * 100).toFixed(0) + '%)');
console.log('  Group    ' + String(tE[1]).padStart(4) + ' (' + (tE[1] / train.length * 100).toFixed(0) + '%)' + String(eE[1]).padStart(6) + ' (' + (eE[1] / ext.length * 100).toFixed(0) + '%)');
console.log('  Cluster  ' + String(tE[2]).padStart(4) + ' (' + (tE[2] / train.length * 100).toFixed(0) + '%)' + String(eE[2]).padStart(6) + ' (' + (eE[2] / ext.length * 100).toFixed(0) + '%)');
console.log();

console.log('======================================================================');
console.log('  DISTANCE METHOD DISTRIBUTION');
console.log('======================================================================\n');

const fdCounts = (arr) => {
  const c = {};
  arr.forEach(g => { c[g.fD] = (c[g.fD] || 0) + 1; });
  return c;
};
const fdLabels = { 1: 'Hubble flow', 2: 'TRGB', 3: 'Cepheids', 4: 'Ursa Major', 5: 'Supernovae' };
const tFD = fdCounts(train), eFD = fdCounts(ext);
console.log('  Method        Training  External');
for (const k of [1, 2, 3, 4, 5]) {
  const tN = tFD[k] || 0, eN = eFD[k] || 0;
  console.log('  ' + (fdLabels[k] || 'fD=' + k).padEnd(14) +
    String(tN).padStart(4) + ' (' + (tN / train.length * 100).toFixed(0) + '%)' +
    String(eN).padStart(6) + ' (' + (eN / ext.length * 100).toFixed(0) + '%)');
}
console.log();

console.log('======================================================================');
console.log('  MORPHOLOGICAL TYPE DISTRIBUTION');
console.log('======================================================================\n');

const morphLabels = { 0: 'S0', 1: 'Sa', 2: 'Sab', 3: 'Sb', 4: 'Sbc', 5: 'Sc', 6: 'Scd', 7: 'Sd', 8: 'Sdm', 9: 'Sm', 10: 'Im', 11: 'BCD' };
const tT = {}, eT = {};
train.forEach(g => { tT[g.T] = (tT[g.T] || 0) + 1; });
ext.forEach(g => { eT[g.T] = (eT[g.T] || 0) + 1; });
console.log('  Type   Training  External');
for (let k = 0; k <= 11; k++) {
  const tN = tT[k] || 0, eN = eT[k] || 0;
  if (tN + eN === 0) continue;
  console.log('  ' + (k + ' ' + (morphLabels[k] || '')).padEnd(8) +
    String(tN).padStart(4) + ' (' + (tN / train.length * 100).toFixed(0) + '%)' +
    String(eN).padStart(6) + ' (' + (eN / ext.length * 100).toFixed(0) + '%)');
}
console.log();

console.log('======================================================================');
console.log('  SUMMARY DASHBOARD');
console.log('======================================================================\n');

const sigVars = results.filter(r => r.significant);
const largeEffect = results.filter(r => r.effectSize === 'LARGE');
const medEffect = results.filter(r => r.effectSize === 'MEDIUM');

console.log('  Total variables compared: ' + results.length);
console.log('  Significantly different (p<0.05): ' + sigVars.length + '/' + results.length);
console.log('    ' + sigVars.map(r => r.key + ' (d=' + r.welch.d + ')').join(', '));
console.log('  Large effect size (|d|>0.8): ' + largeEffect.length);
console.log('    ' + largeEffect.map(r => r.key + ' (d=' + r.welch.d + ')').join(', '));
console.log('  Medium effect size (|d|>0.5): ' + medEffect.length);
console.log('    ' + medEffect.map(r => r.key + ' (d=' + r.welch.d + ')').join(', '));
console.log();

const similar = results.filter(r => !r.significant && r.effectSize === 'NEGLIGIBLE' || r.effectSize === 'SMALL');
console.log('  Similar variables (p>0.05 or small effect): ' + similar.length);
console.log('    ' + similar.map(r => r.key).join(', '));
console.log();

const overallVerdict = sigVars.length >= 4 ? 'VERY DIFFERENT POPULATIONS' :
  sigVars.length >= 2 ? 'MODERATELY DIFFERENT' : 'SIMILAR POPULATIONS';

console.log('  ╔════════════════════════════════════════════════════════════╗');
console.log('  ║  OVERALL VERDICT: ' + overallVerdict.padEnd(39) + '║');
console.log('  ╚════════════════════════════════════════════════════════════╝');
console.log();

if (overallVerdict.includes('DIFFERENT')) {
  console.log('  INTERPRETATION: The training and external samples differ');
  console.log('  significantly on key variables. The law worked on a');
  console.log('  SPECIFIC TYPE of galaxy, not on galaxies in general.');
  console.log('  The biggest differences identify the selection effects:');
  for (const r of [...largeEffect, ...medEffect]) {
    const dir = r.welch.d > 0 ? 'higher' : 'lower';
    console.log('    - ' + r.variable + ': training has ' + dir + ' values (d=' + r.welch.d + ')');
  }
} else {
  console.log('  INTERPRETATION: The two samples are statistically similar.');
  console.log('  This means the law itself is not robust — it fails even');
  console.log('  on galaxies that look like the training set.');
}
console.log();

const output = {
  phase: '82',
  title: 'Sample Comparison — Training vs External',
  nTraining: train.length,
  nExternal: ext.length,
  variableComparisons: results,
  qualityDistribution: { training: tQ, external: eQ },
  environmentDistribution: { training: tE, external: eE },
  distanceMethod: { training: tFD, external: eFD },
  morphology: { training: tT, external: eT },
  nSignificant: sigVars.length,
  significantVariables: sigVars.map(r => r.key),
  overallVerdict,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase82-sample-comparison.json'), JSON.stringify(output, null, 2));
console.log('  Saved: public/phase82-sample-comparison.json');
