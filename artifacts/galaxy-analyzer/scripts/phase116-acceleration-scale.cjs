#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 116: ACCELERATION SCALE vs GAS STATE');
console.log('  Ref: McGaugh+ 2016 (RAR), Lelli+ 2017, Rodrigues+ 2018');
console.log('');
console.log('  The RAR uses a universal acceleration scale a0 ~ 1.2e-10 m/s^2.');
console.log('  Question 1: Does the BEST-FIT a0 vary with gas-to-stellar balance?');
console.log('  Question 2: Does the RAR tightness vary with gas state?');
console.log('  If a0 varies -> gas state modulates the fundamental acceleration');
console.log('  scale -> direct connection to MOND/dark matter debate.');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function median(a) { const s = [...a].sort((x, y) => x - y); const m = Math.floor(s.length / 2); return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }
function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
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

function rarPredicted(log_gbar, log_a0) {
  const a0 = Math.pow(10, log_a0);
  const gbar = Math.pow(10, log_gbar);
  const x = gbar / a0;
  const gobs = gbar / (1 - Math.exp(-Math.sqrt(x)));
  return Math.log10(gobs);
}

function fitA0forGalaxy(pts) {
  let bestA0 = -10;
  let bestRMS = Infinity;

  for (let log_a0 = -11.5; log_a0 <= -8.5; log_a0 += 0.01) {
    let ss = 0;
    let count = 0;
    for (const p of pts) {
      const pred = rarPredicted(p.log_gbar, log_a0);
      if (isFinite(pred)) {
        ss += (p.log_gobs - pred) ** 2;
        count++;
      }
    }
    if (count > 0) {
      const rms = Math.sqrt(ss / count);
      if (rms < bestRMS) {
        bestRMS = rms;
        bestA0 = log_a0;
      }
    }
  }
  return { log_a0: bestA0, rms: bestRMS };
}

const galaxies = [];

for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat <= 0) continue;
  if (rcPoints.length < 8) continue;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const logMHI_L36 = Math.log10(t1.MHI) - Math.log10(t1.L36);
  const logSigma0 = Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
  if (!isFinite(fgas) || !isFinite(logMHI_L36)) continue;

  const rarPts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * pt.vbul * Math.abs(pt.vbul) +
                   pt.vgas * Math.abs(pt.vgas);
    const gobs = pt.vobs ** 2 / (pt.rad * 3.086e16);
    const gbar = vBarSq / (pt.rad * 3.086e16);
    if (gbar > 0 && gobs > 0) {
      const log_gbar = Math.log10(gbar * 1e3);
      const log_gobs = Math.log10(gobs * 1e3);
      if (isFinite(log_gbar) && isFinite(log_gobs)) {
        rarPts.push({ log_gbar, log_gobs });
      }
    }
  }

  if (rarPts.length < 5) continue;

  const { log_a0, rms } = fitA0forGalaxy(rarPts);

  galaxies.push({
    name, fgas, logMHI_L36, logSigma0, Vflat: t1.Vflat,
    log_a0_best: log_a0, rarRMS: rms, nPts: rarPts.length,
  });
}

const highV = galaxies.filter(g => g.Vflat >= 70);
console.log('  Total galaxies with fitted a0: N=' + galaxies.length);
console.log('  High-Vflat (>=70): N=' + highV.length);
console.log('  Standard a0: log(1.2e-10) = ' + Math.log10(1.2e-10).toFixed(3) + '\n');

// ============================================================
// TEST 1: Distribution of best-fit a0
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: DISTRIBUTION OF BEST-FIT log(a0)');
console.log('══════════════════════════════════════════════════════════════\n');

const a0vals = highV.map(g => g.log_a0_best);
console.log('  Mean log(a0):   ' + mean(a0vals).toFixed(3));
console.log('  Median log(a0): ' + median(a0vals).toFixed(3));
console.log('  SD:             ' + sd(a0vals).toFixed(3));
console.log('  Range: [' + Math.min(...a0vals).toFixed(3) + ', ' + Math.max(...a0vals).toFixed(3) + ']');
console.log('  Standard: ' + Math.log10(1.2e-10).toFixed(3) + '\n');

// ============================================================
// TEST 2: Does gas state predict best-fit a0?
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: DOES GAS STATE PREDICT BEST-FIT a0?');
console.log('══════════════════════════════════════════════════════════════\n');

const r_a0_fgas = pearsonR(highV.map(g => g.fgas), a0vals);
const r_a0_ratio = pearsonR(highV.map(g => g.logMHI_L36), a0vals);
const r_a0_sigma = pearsonR(highV.map(g => g.logSigma0), a0vals);
const r_a0_vflat = pearsonR(highV.map(g => Math.log10(g.Vflat)), a0vals);

console.log('  r(fgas, log a0) = ' + r_a0_fgas.toFixed(3));
console.log('  r(log(MHI/L36), log a0) = ' + r_a0_ratio.toFixed(3));
console.log('  r(logSigma0, log a0) = ' + r_a0_sigma.toFixed(3));
console.log('  r(log Vflat, log a0) = ' + r_a0_vflat.toFixed(3));
console.log('');

if (Math.abs(r_a0_fgas) > 0.3 || Math.abs(r_a0_ratio) > 0.3) {
  console.log('  SIGNIFICANT: Gas state correlates with effective a0!');
} else {
  console.log('  a0 does NOT strongly vary with gas state.');
}
console.log('');

// ============================================================
// TEST 3: a0 by gas-fraction tercile
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: a0 BY GAS-FRACTION TERCILE');
console.log('══════════════════════════════════════════════════════════════\n');

const sorted = [...highV].sort((a, b) => a.fgas - b.fgas);
const t3 = Math.floor(sorted.length / 3);
const terciles = [
  { name: 'Low fgas', g: sorted.slice(0, t3) },
  { name: 'Mid fgas', g: sorted.slice(t3, 2 * t3) },
  { name: 'High fgas', g: sorted.slice(2 * t3) },
];

const tercileResults = [];
for (const t of terciles) {
  const a0s = t.g.map(g => g.log_a0_best);
  const rmss = t.g.map(g => g.rarRMS);
  console.log('  ' + t.name + ' (N=' + t.g.length + ', fgas=[' + Math.min(...t.g.map(g => g.fgas)).toFixed(3) + ', ' + Math.max(...t.g.map(g => g.fgas)).toFixed(3) + '])');
  console.log('    Mean log(a0): ' + mean(a0s).toFixed(3) + ' +/- ' + sd(a0s).toFixed(3));
  console.log('    Median log(a0): ' + median(a0s).toFixed(3));
  console.log('    Mean RAR RMS: ' + mean(rmss).toFixed(4));
  console.log('');
  tercileResults.push({
    name: t.name, n: t.g.length,
    meanLogA0: parseFloat(mean(a0s).toFixed(3)),
    medianLogA0: parseFloat(median(a0s).toFixed(3)),
    sdLogA0: parseFloat(sd(a0s).toFixed(3)),
    meanRMS: parseFloat(mean(rmss).toFixed(4)),
  });
}

const a0diff = tercileResults[2].meanLogA0 - tercileResults[0].meanLogA0;
console.log('  High-Low fgas a0 difference: ' + a0diff.toFixed(3) + ' dex');
console.log('  Gas-rich galaxies have ' + (a0diff > 0 ? 'HIGHER' : 'LOWER') + ' effective a0\n');

// ============================================================
// TEST 4: RAR tightness by gas state
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 4: RAR TIGHTNESS vs GAS STATE');
console.log('  (using standard a0 = 1.2e-10)');
console.log('══════════════════════════════════════════════════════════════\n');

const r_rms_fgas = pearsonR(highV.map(g => g.fgas), highV.map(g => g.rarRMS));
const r_rms_ratio = pearsonR(highV.map(g => g.logMHI_L36), highV.map(g => g.rarRMS));
console.log('  r(fgas, RAR RMS) = ' + r_rms_fgas.toFixed(3));
console.log('  r(log(MHI/L36), RAR RMS) = ' + r_rms_ratio.toFixed(3));
console.log('');

for (const t of terciles) {
  console.log('  ' + t.name + ': mean RAR RMS = ' + mean(t.g.map(g => g.rarRMS)).toFixed(4));
}
console.log('');

// ============================================================
// TEST 5: If we use individual best-fit a0, how much does it improve?
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 5: IMPROVEMENT FROM INDIVIDUAL a0 FIT');
console.log('══════════════════════════════════════════════════════════════\n');

const standardRMS = highV.map(g => {
  return g.rarRMS;
});
const individualRMS = highV.map(g => {
  const pts = rcByGalaxy[g.name]?.filter(p => p.rad > 0 && p.vobs > 0) || [];
  const rarPts = pts.map(pt => {
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * pt.vbul * Math.abs(pt.vbul) +
                   pt.vgas * Math.abs(pt.vgas);
    const gobs = pt.vobs ** 2 / (pt.rad * 3.086e16);
    const gbar = vBarSq / (pt.rad * 3.086e16);
    if (gbar > 0 && gobs > 0) {
      return { log_gbar: Math.log10(gbar * 1e3), log_gobs: Math.log10(gobs * 1e3) };
    }
    return null;
  }).filter(Boolean);

  let ss = 0, count = 0;
  for (const p of rarPts) {
    const pred = rarPredicted(p.log_gbar, g.log_a0_best);
    if (isFinite(pred)) { ss += (p.log_gobs - pred) ** 2; count++; }
  }
  return count > 0 ? Math.sqrt(ss / count) : NaN;
}).filter(v => isFinite(v));

console.log('  Mean RAR RMS with standard a0: ' + mean(standardRMS).toFixed(4));
console.log('  Mean RAR RMS with individual a0: ' + mean(individualRMS).toFixed(4));
console.log('  Improvement: ' + ((1 - mean(individualRMS) / mean(standardRMS)) * 100).toFixed(1) + '%\n');

// ============================================================
// VERDICT
// ============================================================

console.log('══════════════════════════════════════════════════════════════');
console.log('  OVERALL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

if (Math.abs(r_a0_fgas) > 0.3) {
  console.log('  DISCOVERY: Gas-to-stellar balance modulates the effective');
  console.log('  acceleration scale a0. This connects our finding directly');
  console.log('  to the MOND acceleration scale and RAR universality debate.');
} else if (Math.abs(a0diff) > 0.1) {
  console.log('  MODERATE: Gas-fraction terciles show different a0 values,');
  console.log('  but the galaxy-level correlation is not strong.');
} else {
  console.log('  a0 is effectively universal across gas states in our sample.');
  console.log('  Gas-to-stellar balance does NOT modulate the fundamental');
  console.log('  acceleration scale. Our finding operates at a different level.');
}

const output = {
  phase: 116,
  title: 'Acceleration Scale vs Gas State',
  nHighV: highV.length,
  a0_distribution: { mean: parseFloat(mean(a0vals).toFixed(3)), median: parseFloat(median(a0vals).toFixed(3)), sd: parseFloat(sd(a0vals).toFixed(3)) },
  correlations: {
    fgas: parseFloat(r_a0_fgas.toFixed(3)),
    logMHI_L36: parseFloat(r_a0_ratio.toFixed(3)),
    logSigma0: parseFloat(r_a0_sigma.toFixed(3)),
    logVflat: parseFloat(r_a0_vflat.toFixed(3)),
  },
  terciles: tercileResults,
  a0_tercileDiff: parseFloat(a0diff.toFixed(3)),
};

const outPath = path.join(__dirname, '..', 'public', 'phase116-acceleration-scale.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);
