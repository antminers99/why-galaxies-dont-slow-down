const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

function seededRNG(seed) { let s = seed | 0; return function() { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; }; }
const rng = seededRNG(20260407);

function normalize(n) { return n.replace(/\s+/g, '').replace(/^NGC0*/, 'NGC').replace(/^DDO0*/, 'DDO').replace(/^IC0*/, 'IC').replace(/^UGC0*/, 'UGC').replace(/^PGC0*/, 'PGC').toUpperCase(); }
function pearsonR(x, y) { const n = x.length; if (n < 3) return NaN; const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n; let sxy = 0, sxx = 0, syy = 0; for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; } return (sxx > 0 && syy > 0) ? sxy / Math.sqrt(sxx * syy) : 0; }
function ols(X, y) { const n = y.length, p = X[0].length; const my = y.reduce((a, b) => a + b, 0) / n; const mX = Array(p).fill(0); for (let j = 0; j < p; j++) { for (let i = 0; i < n; i++) mX[j] += X[i][j]; mX[j] /= n; } const Xc = X.map(r => r.map((v, j) => v - mX[j])); const yc = y.map(v => v - my); const XtX = Array.from({ length: p }, () => Array(p).fill(0)), Xty = Array(p).fill(0); for (let i = 0; i < n; i++) for (let j = 0; j < p; j++) { Xty[j] += Xc[i][j] * yc[i]; for (let k = 0; k < p; k++) XtX[j][k] += Xc[i][j] * Xc[i][k]; } const aug = XtX.map((r, i) => [...r, Xty[i]]); for (let c = 0; c < p; c++) { let mr = c; for (let r = c + 1; r < p; r++) if (Math.abs(aug[r][c]) > Math.abs(aug[mr][c])) mr = r; [aug[c], aug[mr]] = [aug[mr], aug[c]]; if (Math.abs(aug[c][c]) < 1e-14) continue; for (let r = c + 1; r < p; r++) { const f = aug[r][c] / aug[c][c]; for (let j = c; j <= p; j++) aug[r][j] -= f * aug[c][j]; } } const beta = Array(p).fill(0); for (let i = p - 1; i >= 0; i--) { beta[i] = aug[i][p]; for (let j = i + 1; j < p; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= Math.abs(aug[i][i]) > 1e-14 ? aug[i][i] : 1; } const res = []; for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < p; j++) pred += beta[j] * Xc[i][j]; res.push(y[i] - pred); } return { beta, residuals: res }; }
function zscore(arr) { const m = arr.reduce((a, b) => a + b, 0) / arr.length; const s = Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / arr.length); return s > 1e-10 ? arr.map(v => (v - m) / s) : arr.map(() => 0); }

const sparcMap = {}; sparc.forEach(g => { sparcMap[normalize(g.name)] = g; });
const resultsMap = {}; sparcResults.perGalaxy.forEach(g => { resultsMap[normalize(g.name)] = g; });

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[normalize(g.name)]; if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[normalize(g.name)]; if (!sr) continue;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  const hR = sr.models.newtonian.mse > 0 ? Math.log10(Math.max(sr.models.newtonian.mse / Math.max(sr.models.dark_halo_linear.mse, 0.001), 0.01)) : 0;
  gals.push({ name: g.name, logVflat: Math.log10(sp.Vflat), Vflat: sp.Vflat, logMbar: Math.log10(Math.max(Mbar, 1)), logL36: Math.log10(Math.max(sp.L36, 0.001)), logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)), morphT: sp.T, logMHI: g.logMHI, logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)), logA0: g.logA0, hR, dist: sp.D, inc: sp.inc || 0 });
}
const N = gals.length;
const vfR = ols(gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]), gals.map(g => g.logVflat)).residuals;
const a0R = ols(gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]), gals.map(g => g.logA0)).residuals;
const vfRz = zscore(vfR); const a0Rz = zscore(a0R);
const bilZ = zscore(vfRz.map((v, i) => v + a0Rz[i]));
for (let i = 0; i < N; i++) gals[i].DQ = bilZ[i];


function parseFITS(filePath) {
  const buf = fs.readFileSync(filePath);
  const headerBlocks = [];
  let offset = 0, headerDone = false;
  while (!headerDone && offset < buf.length) {
    const block = buf.slice(offset, offset + 2880).toString('ascii'); offset += 2880;
    for (let i = 0; i < 36; i++) { const card = block.substring(i * 80, (i + 1) * 80); headerBlocks.push(card); if (card.startsWith('END')) { headerDone = true; break; } }
  }
  const header = {};
  for (const card of headerBlocks) { if (card.startsWith('END')) break; if (!card.includes('=')) continue; const key = card.substring(0, 8).trim(); let valStr = card.substring(10, 80); const ci = valStr.indexOf('/'); if (ci >= 0 && !valStr.trimStart().startsWith("'")) valStr = valStr.substring(0, ci); valStr = valStr.trim().replace(/'/g, '').trim(); const nv = parseFloat(valStr); header[key] = isNaN(nv) ? valStr : nv; }
  const naxis = header['NAXIS'] || 0;
  const dims = []; for (let i = 1; i <= naxis; i++) dims.push(header['NAXIS' + i] || 0);
  const bitpix = header['BITPIX'] || -32;
  const totalPixels = dims.reduce((a, b) => a * b, 1);
  const bytesPerPixel = Math.abs(bitpix) / 8;
  const data = new Float64Array(totalPixels);
  for (let i = 0; i < totalPixels; i++) {
    const pos = offset + i * bytesPerPixel;
    if (pos + bytesPerPixel > buf.length) { data[i] = NaN; continue; }
    if (bitpix === -32) { const b = Buffer.alloc(4); b[0] = buf[pos]; b[1] = buf[pos+1]; b[2] = buf[pos+2]; b[3] = buf[pos+3]; data[i] = b.readFloatBE(0); }
    else if (bitpix === -64) { const b = Buffer.alloc(8); for (let j = 0; j < 8; j++) b[j] = buf[pos+j]; data[i] = b.readDoubleBE(0); }
    else data[i] = NaN;
    if (isNaN(data[i]) || !isFinite(data[i]) || Math.abs(data[i]) > 1e30) data[i] = NaN;
  }
  if (header['BSCALE'] && header['BSCALE'] !== 1) { const bs = header['BSCALE'], bz = header['BZERO'] || 0; for (let i = 0; i < data.length; i++) if (!isNaN(data[i])) data[i] = data[i] * bs + bz; }
  else if (header['BZERO'] && header['BZERO'] !== 0) { const bz = header['BZERO']; for (let i = 0; i < data.length; i++) if (!isNaN(data[i])) data[i] += bz; }
  return { header, dims, data };
}


function computeMSpectrum(fits, maxMode) {
  const nx = fits.dims[0], ny = fits.dims[1];
  const crpix1 = fits.header['CRPIX1'] || nx / 2, crpix2 = fits.header['CRPIX2'] || ny / 2;
  const cdelt1 = fits.header['CDELT1'] || 1, cdelt2 = fits.header['CDELT2'] || 1;

  const vel = [], coords = [];
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const v = fits.data[j * nx + i]; if (isNaN(v)) continue;
      vel.push(v);
      coords.push({ x: (i + 1 - crpix1) * cdelt1 * 3600, y: (j + 1 - crpix2) * cdelt2 * 3600 });
    }
  }
  if (vel.length < 100) return null;

  const sorted = vel.slice().sort((a, b) => a - b);
  const vsys = sorted[Math.floor(sorted.length / 2)];
  const vrel = vel.map(v => v - vsys);

  const radii = coords.map(c => Math.sqrt(c.x ** 2 + c.y ** 2));
  let maxR = 0; for (let k = 0; k < radii.length; k++) if (radii[k] > maxR) maxR = radii[k];

  const radialBins = 15;
  const azBins = Math.max(2 * maxMode + 2, 16);
  const binW = maxR / radialBins;

  const mPowerTotal = Array(maxMode + 1).fill(0);
  let validBins = 0;

  for (let rb = 0; rb < radialBins; rb++) {
    const rMin = rb * binW, rMax = (rb + 1) * binW;
    const azVals = Array.from({ length: azBins }, () => []);
    for (let k = 0; k < coords.length; k++) {
      const r = radii[k];
      if (r < rMin || r >= rMax) continue;
      const theta = Math.atan2(coords[k].y, coords[k].x);
      const azIdx = Math.floor(((theta + Math.PI) / (2 * Math.PI)) * azBins) % azBins;
      azVals[azIdx].push(vrel[k]);
    }
    const azMeans = azVals.map(a => a.length > 0 ? a.reduce((s, v) => s + v, 0) / a.length : NaN);
    const valid = azMeans.filter(v => !isNaN(v));
    if (valid.length < azBins * 0.5) continue;
    const mean = valid.reduce((a, b) => a + b, 0) / valid.length;

    let cnt = 0;
    for (let m = 0; m <= maxMode; m++) {
      let cm = 0, sm = 0;
      cnt = 0;
      for (let a = 0; a < azBins; a++) {
        if (isNaN(azMeans[a])) continue;
        const ang = (a + 0.5) * 2 * Math.PI / azBins;
        const dv = azMeans[a] - mean;
        if (m === 0) { cm += Math.abs(dv); }
        else { cm += dv * Math.cos(m * ang); sm += dv * Math.sin(m * ang); }
        cnt++;
      }
      if (m === 0) mPowerTotal[0] += cm / cnt;
      else mPowerTotal[m] += Math.sqrt(cm ** 2 + sm ** 2) / cnt;
    }
    validBins++;
  }

  if (validBins < 3) return null;

  const mPower = mPowerTotal.map(p => p / validBins);
  const totalPower = mPower.slice(1).reduce((a, b) => a + b, 0);
  const mFraction = mPower.map((p, i) => i === 0 ? 0 : (totalPower > 0 ? p / totalPower : 0));

  let dominantMode = 1;
  let maxPower = 0;
  for (let m = 1; m <= maxMode; m++) {
    if (mPower[m] > maxPower) { maxPower = mPower[m]; dominantMode = m; }
  }

  const m2Dominance = totalPower > 0 ? mPower[2] / totalPower : 0;
  const higherMFrac = totalPower > 0 ? (mPower.slice(3).reduce((a, b) => a + b, 0)) / totalPower : 0;
  const m2m1Ratio = mPower[1] > 0 ? mPower[2] / mPower[1] : Infinity;

  return {
    mPower,
    mFraction,
    dominantMode,
    m2Dominance,
    higherMFrac,
    m2m1Ratio,
    validBins,
    nValid: vel.length,
  };
}


console.log('='.repeat(72));
console.log('PROGRAM 11 — TEST B: FULL m-SPECTRUM (m=1..6)');
console.log('CDM vs Fuzzy DM discriminator');
console.log('CDM predicts: m=2 dominant from tidal/merger triaxiality');
console.log('FDM predicts: power spread to m=3,4 from quantum interference');
console.log('='.repeat(72));

const dataDir = path.join(__dirname, '..', 'data', 'things-2d');
const thingsGals = [
  { name: 'NGC2841' }, { name: 'NGC5055' }, { name: 'NGC3521' },
  { name: 'NGC6946' }, { name: 'NGC7331' }, { name: 'NGC2403' },
  { name: 'NGC2903' }, { name: 'NGC3198' }, { name: 'NGC4826' },
  { name: 'NGC4736' }, { name: 'DDO154' }, { name: 'IC2574' },
];

const maxMode = 6;
const results = [];

for (const tg of thingsGals) {
  const fitsPath = path.join(dataDir, tg.name + '_MOM1.FITS');
  if (!fs.existsSync(fitsPath)) continue;
  const fStat = fs.statSync(fitsPath);
  if (fStat.size < 1000) continue;

  const galData = gals.find(g => normalize(g.name) === normalize(tg.name));
  if (!galData) { console.log('SKIP ' + tg.name + ' (not in SPARC)'); continue; }

  const fits = parseFITS(fitsPath);
  const spec = computeMSpectrum(fits, maxMode);
  if (!spec) { console.log('SKIP ' + tg.name + ' (too few valid bins)'); continue; }

  results.push({
    name: tg.name,
    DQ: galData.DQ,
    Vflat: galData.Vflat,
    inc: galData.inc,
    logMbar: galData.logMbar,
    ...spec,
  });
}

console.log('\n  Galaxies analyzed: ' + results.length);


console.log('\n\n' + '#'.repeat(72));
console.log('B.1 — PER-GALAXY m-SPECTRUM');
console.log('#'.repeat(72));

console.log('\n  ' + 'Galaxy'.padEnd(12) + 'DQ'.padEnd(8) +
  'm=1'.padEnd(8) + 'm=2'.padEnd(8) + 'm=3'.padEnd(8) + 'm=4'.padEnd(8) + 'm=5'.padEnd(8) + 'm=6'.padEnd(8) + 'Dominant  m2_frac');
console.log('  ' + '-'.repeat(88));

for (const r of results) {
  let line = '  ' + r.name.padEnd(12) + r.DQ.toFixed(2).padEnd(8);
  for (let m = 1; m <= maxMode; m++) {
    line += r.mPower[m].toFixed(1).padEnd(8);
  }
  line += ('m=' + r.dominantMode).padEnd(10) + (r.m2Dominance * 100).toFixed(1) + '%';
  console.log(line);
}

console.log('\n  NORMALIZED FRACTIONS (% of total m>=1 power):');
console.log('  ' + 'Galaxy'.padEnd(12) + 'DQ'.padEnd(8) +
  'f1'.padEnd(8) + 'f2'.padEnd(8) + 'f3'.padEnd(8) + 'f4'.padEnd(8) + 'f5'.padEnd(8) + 'f6'.padEnd(8));
console.log('  ' + '-'.repeat(60));

for (const r of results) {
  let line = '  ' + r.name.padEnd(12) + r.DQ.toFixed(2).padEnd(8);
  for (let m = 1; m <= maxMode; m++) {
    line += ((r.mFraction[m] * 100).toFixed(1) + '%').padEnd(8);
  }
  console.log(line);
}


console.log('\n\n' + '#'.repeat(72));
console.log('B.2 — m=2 DOMINANCE STATISTICS');
console.log('#'.repeat(72));

const m2DomArr = results.map(r => r.m2Dominance);
const meanM2Dom = m2DomArr.reduce((a, b) => a + b, 0) / m2DomArr.length;
const nM2Dominant = results.filter(r => r.dominantMode === 2).length;
const nM1Dominant = results.filter(r => r.dominantMode === 1).length;
const nHigherDom = results.filter(r => r.dominantMode >= 3).length;
const meanHigherFrac = results.reduce((s, r) => s + r.higherMFrac, 0) / results.length;

console.log('\n  Mean m=2 fraction:          ' + (meanM2Dom * 100).toFixed(1) + '%');
console.log('  Mean higher-m (m>=3) frac:  ' + (meanHigherFrac * 100).toFixed(1) + '%');
console.log('  Galaxies with m=1 dominant: ' + nM1Dominant + '/' + results.length);
console.log('  Galaxies with m=2 dominant: ' + nM2Dominant + '/' + results.length);
console.log('  Galaxies with m>=3 dominant: ' + nHigherDom + '/' + results.length);


console.log('\n\n' + '#'.repeat(72));
console.log('B.3 — CORRELATIONS WITH DQ');
console.log('#'.repeat(72));

const validR = results.filter(r => !isNaN(r.DQ));
const dqArr = validR.map(r => r.DQ);

const corrResults = [];
for (let m = 1; m <= maxMode; m++) {
  const mArr = validR.map(r => r.mPower[m]);
  const fArr = validR.map(r => r.mFraction[m]);
  const r_power = pearsonR(dqArr, mArr);
  const r_frac = pearsonR(dqArr, fArr);
  corrResults.push({ m, r_power, r_frac });
  console.log('  r(DQ, m=' + m + ' power)    = ' + r_power.toFixed(4) + '    r(DQ, m=' + m + ' fraction) = ' + r_frac.toFixed(4));
}

const r_DQ_m2dom = pearsonR(dqArr, validR.map(r => r.m2Dominance));
const r_DQ_higher = pearsonR(dqArr, validR.map(r => r.higherMFrac));
console.log('\n  r(DQ, m2 dominance) = ' + r_DQ_m2dom.toFixed(4));
console.log('  r(DQ, higher-m frac) = ' + r_DQ_higher.toFixed(4));


console.log('\n\n' + '#'.repeat(72));
console.log('B.4 — PERMUTATION TESTS');
console.log('#'.repeat(72));

const nPerm = 10000;
const m2FracArr = validR.map(r => r.mFraction[2]);
const obs_r_m2frac = pearsonR(dqArr, m2FracArr);
const obs_r_higher = pearsonR(dqArr, validR.map(r => r.higherMFrac));

let countGE_m2frac = 0, countGE_higher = 0;
for (let p = 0; p < nPerm; p++) {
  const shuffled = dqArr.slice();
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  if (Math.abs(pearsonR(shuffled, m2FracArr)) >= Math.abs(obs_r_m2frac)) countGE_m2frac++;
  if (Math.abs(pearsonR(shuffled, validR.map(r => r.higherMFrac))) >= Math.abs(obs_r_higher)) countGE_higher++;
}

const pVal_m2frac = countGE_m2frac / nPerm;
const pVal_higher = countGE_higher / nPerm;
console.log('\n  r(DQ, m=2 fraction) = ' + obs_r_m2frac.toFixed(4) + ', p = ' + pVal_m2frac.toFixed(4));
console.log('  r(DQ, higher-m frac) = ' + obs_r_higher.toFixed(4) + ', p = ' + pVal_higher.toFixed(4));


console.log('\n\n' + '#'.repeat(72));
console.log('B.5 — HIGH-H vs LOW-H SPECTRUM COMPARISON');
console.log('#'.repeat(72));

const sortedByDQ = results.slice().sort((a, b) => b.DQ - a.DQ);
const nHalf = Math.floor(sortedByDQ.length / 2);
const highH = sortedByDQ.slice(0, nHalf);
const lowH = sortedByDQ.slice(sortedByDQ.length - nHalf);

console.log('\n  High-H (' + highH.length + '): ' + highH.map(g => g.name).join(', '));
console.log('  Low-H  (' + lowH.length + '): ' + lowH.map(g => g.name).join(', '));

console.log('\n  Mode  High-H frac  Low-H frac  Diff     Interpretation');
console.log('  ' + '-'.repeat(60));

for (let m = 1; m <= maxMode; m++) {
  const hMean = highH.reduce((s, g) => s + g.mFraction[m], 0) / highH.length;
  const lMean = lowH.reduce((s, g) => s + g.mFraction[m], 0) / lowH.length;
  const diff = hMean - lMean;
  let interp = '';
  if (m === 2 && diff > 0.05) interp = 'CDM-consistent (m=2 stronger in high-H)';
  else if (m === 2 && diff < -0.05) interp = 'Unexpected (m=2 weaker in high-H)';
  else if (m >= 3 && diff > 0.03) interp = 'FDM-hint? (higher modes stronger in high-H)';
  console.log('  m=' + m + '    ' +
    (hMean * 100).toFixed(1).padEnd(12) + '%' +
    (lMean * 100).toFixed(1).padEnd(12) + '%' +
    ((diff >= 0 ? '+' : '') + (diff * 100).toFixed(1) + 'pp').padEnd(9) +
    interp);
}

const hM2Dom = highH.reduce((s, g) => s + g.m2Dominance, 0) / highH.length;
const lM2Dom = lowH.reduce((s, g) => s + g.m2Dominance, 0) / lowH.length;
console.log('\n  m=2 dominance: High-H=' + (hM2Dom * 100).toFixed(1) + '%  Low-H=' + (lM2Dom * 100).toFixed(1) + '%');


console.log('\n\n' + '#'.repeat(72));
console.log('B.6 — TEST B VERDICT');
console.log('#'.repeat(72));

let testB_verdict, testB_favors;

if (nM2Dominant >= results.length * 0.5 && meanM2Dom > 0.3) {
  testB_verdict = 'M2-DOMINANT';
  testB_favors = 'CDM';
  console.log('\n  VERDICT: m=2 is the DOMINANT azimuthal mode.');
  console.log('  This is CONSISTENT with CDM triaxial halo origin.');
  console.log('  The m-spectrum is concentrated in m=2, not spread to higher modes.');
  console.log('  DISFAVORS Fuzzy DM quantum interference pattern.');
} else if (nHigherDom >= results.length * 0.4 || meanHigherFrac > 0.5) {
  testB_verdict = 'SPREAD-SPECTRUM';
  testB_favors = 'Fuzzy DM';
  console.log('\n  VERDICT: Power is SPREAD across multiple azimuthal modes.');
  console.log('  This is more consistent with Fuzzy DM interference patterns.');
  console.log('  DISFAVORS pure CDM triaxial halo interpretation.');
} else if (nM1Dominant >= results.length * 0.5) {
  testB_verdict = 'M1-DOMINANT';
  testB_favors = 'Lopsidedness/interactions';
  console.log('\n  VERDICT: m=1 (lopsidedness) dominates the spectrum.');
  console.log('  This suggests interaction/accretion history rather than halo shape.');
} else {
  testB_verdict = 'MIXED';
  testB_favors = 'Inconclusive';
  console.log('\n  VERDICT: Mixed m-spectrum — no single mode dominates across all galaxies.');
  console.log('  Need larger sample to distinguish CDM from FDM.');
}

console.log('\n  INTERPRETATION FOR DM IDENTITY:');
if (Math.abs(r_DQ_m2dom) > 0.4) {
  console.log('  m=2 dominance CORRELATES with DQ (r=' + r_DQ_m2dom.toFixed(3) + ').');
  console.log('  Galaxies with stronger H signal show more m=2 dominated spectra.');
  console.log('  This strengthens the halo triaxiality interpretation.');
} else {
  console.log('  m=2 dominance does NOT strongly correlate with DQ.');
  console.log('  The spectral shape is similar across high-H and low-H galaxies.');
}


console.log('\n\n' + '='.repeat(72));
console.log('COMBINED TEST A + B DISCRIMINATION TABLE');
console.log('='.repeat(72));

console.log('\n  ' + 'Test'.padEnd(12) + 'CDM/triaxial'.padEnd(20) + 'SIDM'.padEnd(20) + 'Fuzzy DM'.padEnd(20) + 'Result');
console.log('  ' + '-'.repeat(80));

const testAPath = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results', 'testA-radial-gradient.json');
let testA_favors_str = '(run Test A first)';
if (fs.existsSync(testAPath)) {
  const testAData = JSON.parse(fs.readFileSync(testAPath, 'utf8'));
  testA_favors_str = testAData.favors || testAData.verdict;
  console.log('  ' + 'Test A'.padEnd(12) + 'Flat gradient'.padEnd(20) + 'Core suppressed'.padEnd(20) + 'N/A'.padEnd(20) + testA_favors_str);
}

console.log('  ' + 'Test B'.padEnd(12) + 'm=2 dominant'.padEnd(20) + 'N/A'.padEnd(20) + 'Spread to m=3,4'.padEnd(20) + testB_favors);


const outDir = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results');
if (!fs.existsSync(outDir)) fs.mkdirSync(outDir, { recursive: true });

const outPath = path.join(outDir, 'testB-mspectrum.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: 11,
  test: 'B',
  title: 'Full m-Spectrum (CDM vs Fuzzy DM)',
  timestamp: new Date().toISOString(),
  N: results.length,
  maxMode,
  sampleStats: {
    meanM2Dominance: meanM2Dom,
    meanHigherMFrac: meanHigherFrac,
    nM1Dominant,
    nM2Dominant,
    nHigherDom,
  },
  correlations: corrResults,
  DQ_m2dominance: r_DQ_m2dom,
  DQ_higherFrac: r_DQ_higher,
  permutation: {
    nPerm,
    pVal_m2frac,
    pVal_higher,
  },
  groupComparison: {
    highH_m2dom: hM2Dom,
    lowH_m2dom: lM2Dom,
    highH_galaxies: highH.map(g => g.name),
    lowH_galaxies: lowH.map(g => g.name),
  },
  verdict: testB_verdict,
  favors: testB_favors,
  perGalaxy: results.map(r => ({
    name: r.name, DQ: r.DQ, Vflat: r.Vflat,
    mPower: r.mPower, mFraction: r.mFraction,
    dominantMode: r.dominantMode, m2Dominance: r.m2Dominance,
    higherMFrac: r.higherMFrac,
  })),
}, null, 2));
console.log('\nSaved: ' + outPath);
