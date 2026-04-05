#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 89: LOW-REGIME AUDIT — Chaos or Noise?');
console.log('  Is the low-Vflat scatter real physics or measurement limitation?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

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

function fitA0Bootstrap(pts, nBoot) {
  const results = [];
  for (let b = 0; b < nBoot; b++) {
    const sample = [];
    for (let i = 0; i < pts.length; i++) {
      sample.push(pts[Math.floor(Math.random() * pts.length)]);
    }
    const fit = fitA0(sample);
    if (isFinite(fit.logA0)) results.push(fit.logA0);
  }
  return results;
}

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function median(a) { const s = [...a].sort((x, y) => x - y); return s[Math.floor(s.length / 2)]; }
function percentile(a, p) { const s = [...a].sort((x, y) => x - y); return s[Math.min(Math.floor(s.length * p), s.length - 1)]; }

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
  const eD = parseFloat(line.substring(22, 27).trim());
  const inc = parseFloat(line.substring(30, 34).trim());
  const eInc = parseFloat(line.substring(35, 39).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const eVflat = parseFloat(line.substring(106, 111).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { T, D, fD, eD, inc, eInc, L36, Rdisk, MHI, Vflat, eVflat, Q };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const evobs = parseFloat(line.substring(33, 38).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, evobs, vgas, vdisk, vbul });
}

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
    pts.push({ r: pt.rad, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), vobs: pt.vobs, evobs: pt.evobs || 0 });
  }
  if (pts.length < 3) continue;
  const fit = fitA0(pts);
  if (!isFinite(fit.logA0) || fit.logA0 < 0.5 || fit.logA0 > 5.5) continue;

  const logGbarRange = Math.max(...pts.map(p => p.logGbar)) - Math.min(...pts.map(p => p.logGbar));
  const logGobsRange = Math.max(...pts.map(p => p.logGobs)) - Math.min(...pts.map(p => p.logGobs));
  const medFracErr = median(pts.map(p => p.evobs > 0 && p.vobs > 0 ? p.evobs / p.vobs : 0));
  const maxRad = Math.max(...pts.map(p => p.r));

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

  allGalaxies.push({
    name, nPts: pts.length,
    logA0: fit.logA0, fitRMS: fit.rms,
    logMHI: Math.log10(t1.MHI),
    logMeanRun,
    Vflat: t1.Vflat, Q: t1.Q, D: t1.D, fD: t1.fD, T: t1.T,
    L36: t1.L36, inc: t1.inc,
    eD: t1.eD || 0, eInc: t1.eInc || 0, eVflat: t1.eVflat || 0,
    logGbarRange, logGobsRange, medFracErr, maxRad,
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    rawPts: pts,
  });
}

const lo = allGalaxies.filter(g => g.Vflat < VFLAT_BREAK);
const hi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);

console.log('  Total: ' + allGalaxies.length + ', Lo: ' + lo.length + ', Hi: ' + hi.length + '\n');

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: DATA QUALITY COMPARISON BY REGIME');
console.log('══════════════════════════════════════════════════════════════\n');

function regimeStats(data, label) {
  const stats = {};
  const vars = {
    nPts: g => g.nPts,
    logGbarRange: g => g.logGbarRange,
    logGobsRange: g => g.logGobsRange,
    medFracErr: g => g.medFracErr,
    maxRad: g => g.maxRad,
    fitRMS: g => g.fitRMS,
    D: g => g.D,
    inc: g => g.inc,
    Q: g => g.Q,
    logA0_sd: null,
  };

  for (const [v, fn] of Object.entries(vars)) {
    if (!fn) continue;
    const vals = data.map(fn);
    stats[v] = { mean: +mean(vals).toFixed(3), sd: +sd(vals).toFixed(3), median: +median(vals).toFixed(3) };
  }
  stats.logA0_sd = +sd(data.map(g => g.logA0)).toFixed(3);
  return stats;
}

const loStats = regimeStats(lo, 'Low');
const hiStats = regimeStats(hi, 'High');

console.log('  ' + 'Variable'.padEnd(18) + 'Lo-mean'.padStart(9) + 'Lo-sd'.padStart(9) +
  'Hi-mean'.padStart(9) + 'Hi-sd'.padStart(9) + 'Ratio'.padStart(9));

for (const v of Object.keys(loStats)) {
  if (v === 'logA0_sd') continue;
  const l = loStats[v], h = hiStats[v];
  const ratio = h.mean > 0 ? (l.mean / h.mean).toFixed(2) : '---';
  console.log('  ' + v.padEnd(18) + l.mean.toFixed(3).padStart(9) + l.sd.toFixed(3).padStart(9) +
    h.mean.toFixed(3).padStart(9) + h.sd.toFixed(3).padStart(9) + String(ratio).padStart(9));
}
console.log('  logA0_sd:        Lo=' + loStats.logA0_sd + ', Hi=' + hiStats.logA0_sd +
  ', Ratio=' + (loStats.logA0_sd / hiStats.logA0_sd).toFixed(2));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 2: a0 FIT STABILITY — Bootstrap uncertainty per galaxy');
console.log('══════════════════════════════════════════════════════════════\n');

const N_BOOT = 200;
const bootResults = { lo: [], hi: [] };

for (const g of allGalaxies) {
  if (g.rawPts.length < 5) continue;
  const boots = fitA0Bootstrap(g.rawPts, N_BOOT);
  if (boots.length < 50) continue;
  const bootSD = sd(boots);
  const regime = g.Vflat < VFLAT_BREAK ? 'lo' : 'hi';
  bootResults[regime].push({ name: g.name, logA0: g.logA0, bootSD, nPts: g.nPts, Vflat: g.Vflat });
}

const loBootSD = bootResults.lo.map(b => b.bootSD);
const hiBootSD = bootResults.hi.map(b => b.bootSD);

console.log('  Bootstrap SD of logA0 (N_boot=' + N_BOOT + '):');
console.log('    Low regime (N=' + loBootSD.length + '): mean=' + mean(loBootSD).toFixed(3) +
  ', median=' + median(loBootSD).toFixed(3) + ', P95=' + percentile(loBootSD, 0.95).toFixed(3));
console.log('    High regime (N=' + hiBootSD.length + '): mean=' + mean(hiBootSD).toFixed(3) +
  ', median=' + median(hiBootSD).toFixed(3) + ', P95=' + percentile(hiBootSD, 0.95).toFixed(3));

const loObsSD = sd(lo.map(g => g.logA0));
const hiObsSD = sd(hi.map(g => g.logA0));
const loMeasSD = mean(loBootSD);
const hiMeasSD = mean(hiBootSD);
const loIntrinsic = Math.sqrt(Math.max(0, loObsSD ** 2 - loMeasSD ** 2));
const hiIntrinsic = Math.sqrt(Math.max(0, hiObsSD ** 2 - hiMeasSD ** 2));

console.log('\n  Variance decomposition:');
console.log('    Low: observed SD=' + loObsSD.toFixed(3) + ', measurement SD=' + loMeasSD.toFixed(3) +
  ', intrinsic SD=' + loIntrinsic.toFixed(3) + ' (' + (loMeasSD ** 2 / loObsSD ** 2 * 100).toFixed(1) + '% measurement)');
console.log('    High: observed SD=' + hiObsSD.toFixed(3) + ', measurement SD=' + hiMeasSD.toFixed(3) +
  ', intrinsic SD=' + hiIntrinsic.toFixed(3) + ' (' + (hiMeasSD ** 2 / hiObsSD ** 2 * 100).toFixed(1) + '% measurement)');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: SENSITIVITY TO DISTANCE AND INCLINATION');
console.log('  Perturb D by ±20% and inc by ±5° and re-fit a0');
console.log('══════════════════════════════════════════════════════════════\n');

function perturbAndRefit(galaxies, perturbFn) {
  const deltas = [];
  for (const g of galaxies) {
    if (g.rawPts.length < 5) continue;
    const original = fitA0(g.rawPts);
    const perturbed = perturbFn(g);
    if (!perturbed || perturbed.length < 3) continue;
    const newFit = fitA0(perturbed);
    if (isFinite(newFit.logA0)) {
      deltas.push(Math.abs(newFit.logA0 - original.logA0));
    }
  }
  return deltas;
}

const distPertLo = perturbAndRefit(lo, g => {
  const factor = 1.2;
  return g.rawPts.map(p => ({
    ...p,
    logGobs: Math.log10(Math.pow(10, p.logGobs) * factor * factor),
    logGbar: Math.log10(Math.pow(10, p.logGbar) * factor * factor),
  }));
});

const distPertHi = perturbAndRefit(hi, g => {
  const factor = 1.2;
  return g.rawPts.map(p => ({
    ...p,
    logGobs: Math.log10(Math.pow(10, p.logGobs) * factor * factor),
    logGbar: Math.log10(Math.pow(10, p.logGbar) * factor * factor),
  }));
});

console.log('  Distance perturbation (+20%):');
console.log('    Low: mean|Δ logA0| = ' + mean(distPertLo).toFixed(3) + ', median = ' + median(distPertLo).toFixed(3));
console.log('    High: mean|Δ logA0| = ' + mean(distPertHi).toFixed(3) + ', median = ' + median(distPertHi).toFixed(3));
console.log('    Ratio (Lo/Hi): ' + (mean(distPertLo) / mean(distPertHi)).toFixed(2));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: DYNAMIC RANGE vs a0 SCATTER');
console.log('  Does limited gbar range in dwarfs cause a0 instability?');
console.log('══════════════════════════════════════════════════════════════\n');

const rGbarRange_lo = pearsonR(lo.map(g => g.logGbarRange), lo.map(g => Math.abs(g.logA0 - mean(lo.map(g2 => g2.logA0)))));
const rGbarRange_hi = pearsonR(hi.map(g => g.logGbarRange), hi.map(g => Math.abs(g.logA0 - mean(hi.map(g2 => g2.logA0)))));
const rNpts_lo = pearsonR(lo.map(g => g.nPts), lo.map(g => Math.abs(g.logA0 - mean(lo.map(g2 => g2.logA0)))));
const rNpts_hi = pearsonR(hi.map(g => g.nPts), hi.map(g => Math.abs(g.logA0 - mean(hi.map(g2 => g2.logA0)))));

console.log('  r(gbarRange, |residual_a0|):');
console.log('    Low: ' + rGbarRange_lo.toFixed(3) + '  High: ' + rGbarRange_hi.toFixed(3));
console.log('  r(nPts, |residual_a0|):');
console.log('    Low: ' + rNpts_lo.toFixed(3) + '  High: ' + rNpts_hi.toFixed(3));

const loWellConstrained = lo.filter(g => g.nPts >= 10 && g.logGbarRange >= 1.0);
const loOthers = lo.filter(g => g.nPts < 10 || g.logGbarRange < 1.0);

console.log('\n  Well-constrained low-Vflat (nPts>=10, gRange>=1.0 dex):');
console.log('    N=' + loWellConstrained.length + ', sd(logA0)=' + sd(loWellConstrained.map(g => g.logA0)).toFixed(3));
console.log('  Poorly-constrained low-Vflat:');
console.log('    N=' + loOthers.length + ', sd(logA0)=' + sd(loOthers.map(g => g.logA0)).toFixed(3));

if (loWellConstrained.length >= 10) {
  const rMHI_wc = pearsonR(loWellConstrained.map(g => g.logMHI), loWellConstrained.map(g => g.logA0));
  console.log('    Well-constrained r(MHI,a0) = ' + rMHI_wc.toFixed(3));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: QUALITY-STRATIFIED ANALYSIS IN LOW REGIME');
console.log('══════════════════════════════════════════════════════════════\n');

const qGroups = [
  { name: 'Q=1 (best)', data: lo.filter(g => g.Q === 1) },
  { name: 'Q=2', data: lo.filter(g => g.Q === 2) },
  { name: 'Q=3 (worst)', data: lo.filter(g => g.Q === 3) },
];

for (const qg of qGroups) {
  if (qg.data.length < 5) { console.log('  ' + qg.name + ': N=' + qg.data.length + ' (too few)'); continue; }
  console.log('  ' + qg.name + ': N=' + qg.data.length +
    ', sd(logA0)=' + sd(qg.data.map(g => g.logA0)).toFixed(3) +
    ', mean nPts=' + mean(qg.data.map(g => g.nPts)).toFixed(1) +
    ', mean gRange=' + mean(qg.data.map(g => g.logGbarRange)).toFixed(2));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 6: BIMODALITY / CLUSTERING IN LOW-REGIME logA0');
console.log('══════════════════════════════════════════════════════════════\n');

const loLogA0 = lo.map(g => g.logA0).sort((a, b) => a - b);
const loMean = mean(loLogA0);
const loSD = sd(loLogA0);

const bins = 10;
const binMin = loMean - 3 * loSD;
const binMax = loMean + 3 * loSD;
const binWidth = (binMax - binMin) / bins;
const histogram = Array(bins).fill(0);
for (const v of loLogA0) {
  const idx = Math.min(bins - 1, Math.max(0, Math.floor((v - binMin) / binWidth)));
  histogram[idx]++;
}

console.log('  Histogram of logA0 (low regime):');
for (let i = 0; i < bins; i++) {
  const center = binMin + (i + 0.5) * binWidth;
  const bar = '#'.repeat(histogram[i]);
  console.log('    ' + center.toFixed(2) + ' | ' + bar + ' (' + histogram[i] + ')');
}

const skewness = lo.reduce((s, g) => s + ((g.logA0 - loMean) / loSD) ** 3, 0) / lo.length;
const kurtosis = lo.reduce((s, g) => s + ((g.logA0 - loMean) / loSD) ** 4, 0) / lo.length - 3;
console.log('\n  Skewness: ' + skewness.toFixed(3) + ' (0 = symmetric)');
console.log('  Excess kurtosis: ' + kurtosis.toFixed(3) + ' (0 = Gaussian, <0 = flat, >0 = peaked)');

const dipGap = [];
for (let i = 1; i < loLogA0.length; i++) {
  dipGap.push({ gap: loLogA0[i] - loLogA0[i - 1], at: (loLogA0[i] + loLogA0[i - 1]) / 2 });
}
dipGap.sort((a, b) => b.gap - a.gap);
console.log('  Largest gaps in sorted logA0:');
for (let i = 0; i < Math.min(3, dipGap.length); i++) {
  console.log('    ' + dipGap[i].gap.toFixed(3) + ' dex at logA0=' + dipGap[i].at.toFixed(2));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 7: MORPHOLOGY SUB-REGIMES IN LOW-Vflat');
console.log('══════════════════════════════════════════════════════════════\n');

const morphGroups = [
  { name: 'T <= 5 (early-type)', data: lo.filter(g => g.T <= 5) },
  { name: 'T 6-8 (late spirals)', data: lo.filter(g => g.T >= 6 && g.T <= 8) },
  { name: 'T >= 9 (irregulars)', data: lo.filter(g => g.T >= 9) },
];

for (const mg of morphGroups) {
  if (mg.data.length < 5) { console.log('  ' + mg.name + ': N=' + mg.data.length + ' (too few)'); continue; }
  console.log('  ' + mg.name + ': N=' + mg.data.length +
    ', mean logA0=' + mean(mg.data.map(g => g.logA0)).toFixed(3) +
    ', sd=' + sd(mg.data.map(g => g.logA0)).toFixed(3) +
    ', mean Vflat=' + mean(mg.data.map(g => g.Vflat)).toFixed(1));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  FINAL SUMMARY');
console.log('══════════════════════════════════════════════════════════════\n');

const measFracLo = loMeasSD ** 2 / loObsSD ** 2;
const measFracHi = hiMeasSD ** 2 / hiObsSD ** 2;

console.log('  Observed scatter:    Lo=' + loObsSD.toFixed(3) + ' dex, Hi=' + hiObsSD.toFixed(3) + ' dex');
console.log('  Measurement scatter: Lo=' + loMeasSD.toFixed(3) + ' dex, Hi=' + hiMeasSD.toFixed(3) + ' dex');
console.log('  Intrinsic scatter:   Lo=' + loIntrinsic.toFixed(3) + ' dex, Hi=' + hiIntrinsic.toFixed(3) + ' dex');
console.log('  Measurement fraction: Lo=' + (measFracLo * 100).toFixed(1) + '%, Hi=' + (measFracHi * 100).toFixed(1) + '%');
console.log();

let verdict;
if (measFracLo > 0.5) {
  verdict = 'MEASUREMENT-DOMINATED';
  console.log('  ╔══════════════════════════════════════════════════════════════════╗');
  console.log('  ║  LOW REGIME IS MEASUREMENT-DOMINATED                           ║');
  console.log('  ║  The "chaos" is largely observational noise, not physics        ║');
  console.log('  ╚══════════════════════════════════════════════════════════════════╝');
} else if (measFracLo > 0.2) {
  verdict = 'MIXED';
  console.log('  ╔══════════════════════════════════════════════════════════════════╗');
  console.log('  ║  LOW REGIME: MIXED (measurement + intrinsic)                   ║');
  console.log('  ║  Substantial intrinsic scatter exists, but measurement is also  ║');
  console.log('  ║  a significant contributor                                     ║');
  console.log('  ╚══════════════════════════════════════════════════════════════════╝');
} else {
  verdict = 'INTRINSIC';
  console.log('  ╔══════════════════════════════════════════════════════════════════╗');
  console.log('  ║  LOW REGIME: INTRINSIC SCATTER DOMINATES                       ║');
  console.log('  ║  The large a0 variation in dwarfs is real physics              ║');
  console.log('  ╚══════════════════════════════════════════════════════════════════╝');
}

const output = {
  phase: '89',
  title: 'Low-Regime Audit — Chaos or Noise?',
  nLo: lo.length, nHi: hi.length,
  dataQuality: { lo: loStats, hi: hiStats },
  bootstrap: {
    lo: { N: loBootSD.length, meanSD: +mean(loBootSD).toFixed(3), medianSD: +median(loBootSD).toFixed(3) },
    hi: { N: hiBootSD.length, meanSD: +mean(hiBootSD).toFixed(3), medianSD: +median(hiBootSD).toFixed(3) },
  },
  varianceDecomposition: {
    lo: { observed: +loObsSD.toFixed(3), measurement: +loMeasSD.toFixed(3), intrinsic: +loIntrinsic.toFixed(3), measFrac: +(measFracLo * 100).toFixed(1) },
    hi: { observed: +hiObsSD.toFixed(3), measurement: +hiMeasSD.toFixed(3), intrinsic: +hiIntrinsic.toFixed(3), measFrac: +(measFracHi * 100).toFixed(1) },
  },
  dynamicRange: {
    rGbarRange_lo: +rGbarRange_lo.toFixed(3), rGbarRange_hi: +rGbarRange_hi.toFixed(3),
    rNpts_lo: +rNpts_lo.toFixed(3), rNpts_hi: +rNpts_hi.toFixed(3),
    wellConstrained: { N: loWellConstrained.length, sd: +sd(loWellConstrained.map(g => g.logA0)).toFixed(3) },
  },
  distribution: { skewness: +skewness.toFixed(3), excessKurtosis: +kurtosis.toFixed(3) },
  verdict,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase89-low-regime-audit.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase89-low-regime-audit.json');
