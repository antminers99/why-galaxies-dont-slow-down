#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

function seededRNG(seed) { let s = seed | 0; return function() { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; }; }
const rng = seededRNG(20260410);

function normalize(n) { return n.replace(/\s+/g, '').replace(/^NGC0*/, 'NGC').replace(/^DDO0*/, 'DDO').replace(/^IC0*/, 'IC').replace(/^UGC0*/, 'UGC').replace(/^PGC0*/, 'PGC').toUpperCase(); }
function pearsonR(x, y) { const n = x.length; if (n < 3) return NaN; const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n; let sxy = 0, sxx = 0, syy = 0; for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; } return (sxx > 0 && syy > 0) ? sxy / Math.sqrt(sxx * syy) : 0; }
function ols(X, y) { const n = y.length, p = X[0].length; const my = y.reduce((a, b) => a + b, 0) / n; const mX = Array(p).fill(0); for (let j = 0; j < p; j++) { for (let i = 0; i < n; i++) mX[j] += X[i][j]; mX[j] /= n; } const Xc = X.map(r => r.map((v, j) => v - mX[j])); const yc = y.map(v => v - my); const XtX = Array.from({ length: p }, () => Array(p).fill(0)), Xty = Array(p).fill(0); for (let i = 0; i < n; i++) for (let j = 0; j < p; j++) { Xty[j] += Xc[i][j] * yc[i]; for (let k = 0; k < p; k++) XtX[j][k] += Xc[i][j] * Xc[i][k]; } const aug = XtX.map((r, i) => [...r, Xty[i]]); for (let c = 0; c < p; c++) { let mr = c; for (let r = c + 1; r < p; r++) if (Math.abs(aug[r][c]) > Math.abs(aug[mr][c])) mr = r; [aug[c], aug[mr]] = [aug[mr], aug[c]]; if (Math.abs(aug[c][c]) < 1e-14) continue; for (let r = c + 1; r < p; r++) { const f = aug[r][c] / aug[c][c]; for (let j = c; j <= p; j++) aug[r][j] -= f * aug[c][j]; } } const beta = Array(p).fill(0); for (let i = p - 1; i >= 0; i--) { beta[i] = aug[i][p]; for (let j = i + 1; j < p; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= Math.abs(aug[i][i]) > 1e-14 ? aug[i][i] : 1; } const res = []; for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < p; j++) pred += beta[j] * Xc[i][j]; res.push(y[i] - pred); } return { beta, residuals: res }; }
function zscore(arr) { const m = arr.reduce((a, b) => a + b, 0) / arr.length; const s = Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / arr.length); return s > 1e-10 ? arr.map(v => (v - m) / s) : arr.map(() => 0); }
function permutationP(x, y, nPerm, observed) { let count = 0; const yC = y.slice(); for (let p = 0; p < nPerm; p++) { for (let i = yC.length - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [yC[i], yC[j]] = [yC[j], yC[i]]; } if (Math.abs(pearsonR(x, yC)) >= Math.abs(observed)) count++; } return count / Math.max(nPerm, 1); }

const sparcMap = {}; sparc.forEach(g => { sparcMap[normalize(g.name)] = g; });
const resultsMap = {}; sparcResults.perGalaxy.forEach(g => { resultsMap[normalize(g.name)] = g; });
const gals = [];
for (const g of d56.perGalaxy) { const sp = sparcMap[normalize(g.name)]; if (!sp || sp.Vflat <= 0) continue; const sr = resultsMap[normalize(g.name)]; if (!sr) continue; const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9; const hR = sr.models.newtonian.mse > 0 ? Math.log10(Math.max(sr.models.newtonian.mse / Math.max(sr.models.dark_halo_linear.mse, 0.001), 0.01)) : 0; gals.push({ name: g.name, logVflat: Math.log10(sp.Vflat), Vflat: sp.Vflat, logMbar: Math.log10(Math.max(Mbar, 1)), logL36: Math.log10(Math.max(sp.L36, 0.001)), logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)), morphT: sp.T, logMHI: g.logMHI, logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)), logA0: g.logA0, hR, dist: sp.D, inc: sp.inc || 0 }); }
const N = gals.length;
const vfR = ols(gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]), gals.map(g => g.logVflat)).residuals;
const a0R = ols(gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]), gals.map(g => g.logA0)).residuals;
const vfRz = zscore(vfR); const a0Rz = zscore(a0R);
const bilZ = zscore(vfRz.map((v, i) => v + a0Rz[i]));
for (let i = 0; i < N; i++) gals[i].DQ = bilZ[i];

function parseFITS(filePath) {
  const buf = fs.readFileSync(filePath); const headerBlocks = []; let offset = 0, headerDone = false;
  while (!headerDone && offset < buf.length) { const block = buf.slice(offset, offset + 2880).toString('ascii'); offset += 2880; for (let i = 0; i < 36; i++) { const card = block.substring(i * 80, (i + 1) * 80); headerBlocks.push(card); if (card.startsWith('END')) { headerDone = true; break; } } }
  const header = {}; for (const card of headerBlocks) { if (card.startsWith('END')) break; if (!card.includes('=')) continue; const key = card.substring(0, 8).trim(); let valStr = card.substring(10, 80); const ci = valStr.indexOf('/'); if (ci >= 0 && !valStr.trimStart().startsWith("'")) valStr = valStr.substring(0, ci); valStr = valStr.trim().replace(/'/g, '').trim(); const nv = parseFloat(valStr); header[key] = isNaN(nv) ? valStr : nv; }
  const naxis = header['NAXIS'] || 0; const dims = []; for (let i = 1; i <= naxis; i++) dims.push(header['NAXIS' + i] || 0);
  const bitpix = header['BITPIX'] || -32; const totalPixels = dims.reduce((a, b) => a * b, 1); const bytesPerPixel = Math.abs(bitpix) / 8;
  const data = new Float64Array(totalPixels);
  for (let i = 0; i < totalPixels; i++) { const pos = offset + i * bytesPerPixel; if (pos + bytesPerPixel > buf.length) { data[i] = NaN; continue; } if (bitpix === -32) { const b = Buffer.alloc(4); b[0] = buf[pos]; b[1] = buf[pos+1]; b[2] = buf[pos+2]; b[3] = buf[pos+3]; data[i] = b.readFloatBE(0); } else if (bitpix === -64) { const b = Buffer.alloc(8); for (let j = 0; j < 8; j++) b[j] = buf[pos+j]; data[i] = b.readDoubleBE(0); } else data[i] = NaN; if (isNaN(data[i]) || !isFinite(data[i]) || Math.abs(data[i]) > 1e30) data[i] = NaN; }
  if (header['BSCALE'] && header['BSCALE'] !== 1) { const bs = header['BSCALE'], bz = header['BZERO'] || 0; for (let i = 0; i < data.length; i++) if (!isNaN(data[i])) data[i] = data[i] * bs + bz; }
  else if (header['BZERO'] && header['BZERO'] !== 0) { const bz = header['BZERO']; for (let i = 0; i < data.length; i++) if (!isNaN(data[i])) data[i] += bz; }
  return { header, dims, data };
}

function computeShapeAmplitude(fits, maxMode, nRadBins, azBins, rPctCut, excludeInnerFrac) {
  const nx = fits.dims[0], ny = fits.dims[1];
  const crpix1 = fits.header['CRPIX1'] || nx / 2, crpix2 = fits.header['CRPIX2'] || ny / 2;
  const cdelt1 = fits.header['CDELT1'] || 1, cdelt2 = fits.header['CDELT2'] || 1;
  const vel = [], coords = [];
  for (let j = 0; j < ny; j++) { for (let i = 0; i < nx; i++) { const v = fits.data[j * nx + i]; if (isNaN(v)) continue; vel.push(v); coords.push({ x: (i + 1 - crpix1) * cdelt1 * 3600, y: (j + 1 - crpix2) * cdelt2 * 3600 }); } }
  if (vel.length < 100) return null;
  const sorted = vel.slice().sort((a, b) => a - b);
  const vsys = sorted[Math.floor(sorted.length / 2)];
  const vrel = vel.map(v => v - vsys);
  const radii = coords.map(c => Math.sqrt(c.x ** 2 + c.y ** 2));
  const rSorted = radii.slice().sort((a, b) => a - b);
  const rMax = rSorted[Math.floor(rSorted.length * rPctCut)];
  const rMinCut = excludeInnerFrac * rMax;

  const binW = (rMax - rMinCut) / nRadBins;
  if (binW <= 0) return null;
  const bins = [];

  for (let rb = 0; rb < nRadBins; rb++) {
    const rMin = rMinCut + rb * binW, rMaxBin = rMinCut + (rb + 1) * binW;
    const azVals = Array.from({ length: azBins }, () => []);
    for (let k = 0; k < coords.length; k++) {
      if (radii[k] < rMin || radii[k] >= rMaxBin) continue;
      const theta = Math.atan2(coords[k].y, coords[k].x);
      const azIdx = Math.floor(((theta + Math.PI) / (2 * Math.PI)) * azBins) % azBins;
      azVals[azIdx].push(vrel[k]);
    }
    const azMeans = azVals.map(a => a.length > 0 ? a.reduce((s, v) => s + v, 0) / a.length : NaN);
    const valid = azMeans.filter(v => !isNaN(v));
    if (valid.length < azBins * 0.5) continue;
    const mean = valid.reduce((a, b) => a + b, 0) / valid.length;

    let nonAxiSum = 0;
    for (let m = 2; m <= maxMode; m++) {
      let cm = 0, sm = 0, cnt = 0;
      for (let a = 0; a < azBins; a++) {
        if (isNaN(azMeans[a])) continue;
        const ang = (a + 0.5) * 2 * Math.PI / azBins;
        const dv = azMeans[a] - mean;
        cm += dv * Math.cos(m * ang); sm += dv * Math.sin(m * ang);
        cnt++;
      }
      nonAxiSum += Math.sqrt(cm ** 2 + sm ** 2) / Math.max(cnt, 1);
    }
    bins.push(nonAxiSum);
  }
  if (bins.length < 3) return null;
  return bins.reduce((a, b) => a + b, 0) / bins.length;
}

function computeNormalizedShapeAmp(fits, maxMode, nRadBins, azBins, rPctCut, excludeInnerFrac, normMode) {
  const nx = fits.dims[0], ny = fits.dims[1];
  const crpix1 = fits.header['CRPIX1'] || nx / 2, crpix2 = fits.header['CRPIX2'] || ny / 2;
  const cdelt1 = fits.header['CDELT1'] || 1, cdelt2 = fits.header['CDELT2'] || 1;
  const vel = [], coords = [];
  for (let j = 0; j < ny; j++) { for (let i = 0; i < nx; i++) { const v = fits.data[j * nx + i]; if (isNaN(v)) continue; vel.push(v); coords.push({ x: (i + 1 - crpix1) * cdelt1 * 3600, y: (j + 1 - crpix2) * cdelt2 * 3600 }); } }
  if (vel.length < 100) return null;
  const sorted = vel.slice().sort((a, b) => a - b);
  const vsys = sorted[Math.floor(sorted.length / 2)];
  const vrel = vel.map(v => v - vsys);
  const radii = coords.map(c => Math.sqrt(c.x ** 2 + c.y ** 2));
  const rSorted = radii.slice().sort((a, b) => a - b);
  const rMax = rSorted[Math.floor(rSorted.length * rPctCut)];
  const rMinCut = excludeInnerFrac * rMax;
  const binW = (rMax - rMinCut) / nRadBins;
  if (binW <= 0) return null;
  const bins = [];

  for (let rb = 0; rb < nRadBins; rb++) {
    const rMin = rMinCut + rb * binW, rMaxBin = rMinCut + (rb + 1) * binW;
    const azVals = Array.from({ length: azBins }, () => []);
    for (let k = 0; k < coords.length; k++) {
      if (radii[k] < rMin || radii[k] >= rMaxBin) continue;
      const theta = Math.atan2(coords[k].y, coords[k].x);
      const azIdx = Math.floor(((theta + Math.PI) / (2 * Math.PI)) * azBins) % azBins;
      azVals[azIdx].push(vrel[k]);
    }
    const azMeans = azVals.map(a => a.length > 0 ? a.reduce((s, v) => s + v, 0) / a.length : NaN);
    const valid = azMeans.filter(v => !isNaN(v));
    if (valid.length < azBins * 0.5) continue;
    const mean = valid.reduce((a, b) => a + b, 0) / valid.length;
    const vrange = Math.max(...valid) - Math.min(...valid);

    let nonAxiSum = 0, totalSum = 0;
    for (let m = 1; m <= maxMode; m++) {
      let cm = 0, sm = 0, cnt = 0;
      for (let a = 0; a < azBins; a++) {
        if (isNaN(azMeans[a])) continue;
        const ang = (a + 0.5) * 2 * Math.PI / azBins;
        const dv = azMeans[a] - mean;
        cm += dv * Math.cos(m * ang); sm += dv * Math.sin(m * ang);
        cnt++;
      }
      const amp = Math.sqrt(cm ** 2 + sm ** 2) / Math.max(cnt, 1);
      totalSum += amp;
      if (m >= 2) nonAxiSum += amp;
    }

    if (normMode === 'fractional') {
      bins.push(totalSum > 0 ? nonAxiSum / totalSum : 0);
    } else if (normMode === 'velocity') {
      bins.push(vrange > 0 ? nonAxiSum / vrange : 0);
    } else {
      bins.push(nonAxiSum);
    }
  }
  if (bins.length < 3) return null;
  return bins.reduce((a, b) => a + b, 0) / bins.length;
}


console.log('='.repeat(72));
console.log('PROGRAM 12 / DM-4V: VERIFICATION OF shapeAmplitude');
console.log('4 stress tests: binning, normalization, LOO+pairs, bar exclusion');
console.log('='.repeat(72));

const dataDir = path.join(__dirname, '..', 'data', 'things-2d');
const thingsGals = [
  { name: 'NGC2841' }, { name: 'NGC5055' }, { name: 'NGC3521' },
  { name: 'NGC6946' }, { name: 'NGC7331' }, { name: 'NGC2403' },
  { name: 'NGC2903' }, { name: 'NGC3198' },
];
const maxMode = 6;

const fitsData = {};
const galDQ = {};
for (const tg of thingsGals) {
  const fitsPath = path.join(dataDir, tg.name + '_MOM1.FITS');
  if (!fs.existsSync(fitsPath) || fs.statSync(fitsPath).size < 1000) continue;
  const galData = gals.find(g => normalize(g.name) === normalize(tg.name));
  if (!galData) { console.log('  SKIP ' + tg.name); continue; }
  fitsData[tg.name] = parseFITS(fitsPath);
  galDQ[tg.name] = galData.DQ;
}

const galNames = Object.keys(fitsData);
console.log('\n  Galaxies: ' + galNames.length + ' — ' + galNames.join(', '));


console.log('\n' + '='.repeat(72));
console.log('V1: BINNING SENSITIVITY');
console.log('Does r(DQ, shapeAmp) survive different radial/azimuthal bin counts?');
console.log('='.repeat(72));

const binConfigs = [
  { nRad: 8, nAz: 12, label: '8rad x 12az' },
  { nRad: 10, nAz: 16, label: '10rad x 16az' },
  { nRad: 12, nAz: 16, label: '12rad x 16az (DM-4 default)' },
  { nRad: 12, nAz: 20, label: '12rad x 20az' },
  { nRad: 14, nAz: 16, label: '14rad x 16az' },
  { nRad: 16, nAz: 16, label: '16rad x 16az' },
  { nRad: 10, nAz: 12, label: '10rad x 12az' },
  { nRad: 14, nAz: 20, label: '14rad x 20az' },
];

const binResults = [];
for (const bc of binConfigs) {
  const amps = [], dqs = [];
  for (const name of galNames) {
    const sa = computeShapeAmplitude(fitsData[name], maxMode, bc.nRad, bc.nAz, 0.90, 0.0);
    if (sa !== null) { amps.push(sa); dqs.push(galDQ[name]); }
  }
  const r = pearsonR(dqs, amps);
  binResults.push({ label: bc.label, r, n: dqs.length });
  console.log('  ' + bc.label + ': r=' + r.toFixed(3) + ' (N=' + dqs.length + ')');
}

const binRValues = binResults.map(b => b.r);
const binRMean = binRValues.reduce((a, b) => a + b, 0) / binRValues.length;
const binRMin = Math.min(...binRValues);
const binRMax = Math.max(...binRValues);
const binRRange = binRMax - binRMin;
const binAllPositive = binRValues.every(r => r > 0.5);

console.log('\n  Mean r across configs: ' + binRMean.toFixed(3));
console.log('  Range: [' + binRMin.toFixed(3) + ', ' + binRMax.toFixed(3) + '] spread=' + binRRange.toFixed(3));
console.log('  All r > 0.5: ' + binAllPositive);

let V1_verdict;
if (binAllPositive && binRRange < 0.15) V1_verdict = 'ROBUST — shapeAmplitude insensitive to binning choice';
else if (binAllPositive && binRRange < 0.25) V1_verdict = 'STABLE — moderate sensitivity but sign preserved';
else V1_verdict = 'SENSITIVE — result depends on binning';

console.log('  V1 VERDICT: ' + V1_verdict);


console.log('\n' + '='.repeat(72));
console.log('V2: NORMALIZATION SENSITIVITY');
console.log('Does the result hold with raw, fractional, and velocity-normalized amps?');
console.log('='.repeat(72));

const normModes = ['raw', 'fractional', 'velocity'];
const normResults = [];

for (const nm of normModes) {
  const amps = [], dqs = [];
  for (const name of galNames) {
    let sa;
    if (nm === 'raw') {
      sa = computeShapeAmplitude(fitsData[name], maxMode, 12, 16, 0.90, 0.0);
    } else {
      sa = computeNormalizedShapeAmp(fitsData[name], maxMode, 12, 16, 0.90, 0.0, nm);
    }
    if (sa !== null) { amps.push(sa); dqs.push(galDQ[name]); }
  }
  const r = pearsonR(dqs, amps);
  const p = permutationP(dqs, amps, 2000, r);
  normResults.push({ mode: nm, r, p, n: dqs.length });
  console.log('  ' + nm + ': r=' + r.toFixed(3) + ' p=' + p.toFixed(4) + ' (N=' + dqs.length + ')');
}

const allNormPositive = normResults.every(n => n.r > 0.3);
const bestNorm = normResults.reduce((a, b) => Math.abs(a.r) > Math.abs(b.r) ? a : b);

console.log('\n  All normalizations r > 0.3: ' + allNormPositive);
console.log('  Best normalization: ' + bestNorm.mode + ' (r=' + bestNorm.r.toFixed(3) + ')');

let V2_verdict;
if (allNormPositive) V2_verdict = 'ROBUST — result survives all normalizations';
else V2_verdict = 'SENSITIVE — some normalizations break the correlation';

console.log('  V2 VERDICT: ' + V2_verdict);


console.log('\n' + '='.repeat(72));
console.log('V3: LOO + PAIR ROBUSTNESS (DEEP)');
console.log('='.repeat(72));

const baseAmps = [];
for (const name of galNames) {
  baseAmps.push({ name, DQ: galDQ[name], sa: computeShapeAmplitude(fitsData[name], maxMode, 12, 16, 0.90, 0.0) });
}
const validBase = baseAmps.filter(g => g.sa !== null);
const baseDQs = validBase.map(g => g.DQ);
const baseSAs = validBase.map(g => g.sa);
const baseR = pearsonR(baseDQs, baseSAs);

console.log('\n  Full sample r(DQ, shapeAmp): ' + baseR.toFixed(3));

console.log('\n  LOO:');
const looResults = [];
for (let drop = 0; drop < validBase.length; drop++) {
  const sub = validBase.filter((_, i) => i !== drop);
  const r = pearsonR(sub.map(g => g.DQ), sub.map(g => g.sa));
  looResults.push({ dropped: validBase[drop].name, r });
  console.log('    Drop ' + validBase[drop].name + ': r=' + r.toFixed(3));
}
const looAllPositive = looResults.every(l => l.r > 0);
const looMin = Math.min(...looResults.map(l => l.r));
const looMax = Math.max(...looResults.map(l => l.r));
const nLOOAbove05 = looResults.filter(l => l.r > 0.5).length;
console.log('    All positive: ' + looAllPositive);
console.log('    Range: [' + looMin.toFixed(3) + ', ' + looMax.toFixed(3) + ']');
console.log('    N with r > 0.5: ' + nLOOAbove05 + '/' + looResults.length);

console.log('\n  All-pairs rank consistency:');
let concordant = 0, discordant = 0;
for (let i = 0; i < validBase.length; i++) {
  for (let j = i + 1; j < validBase.length; j++) {
    const dqDiff = validBase[i].DQ - validBase[j].DQ;
    const saDiff = validBase[i].sa - validBase[j].sa;
    if (dqDiff * saDiff > 0) concordant++;
    else if (dqDiff * saDiff < 0) discordant++;
  }
}
const totalPairs = concordant + discordant;
const kendallTau = totalPairs > 0 ? (concordant - discordant) / totalPairs : 0;
console.log('    Concordant: ' + concordant + ', Discordant: ' + discordant);
console.log('    Kendall tau: ' + kendallTau.toFixed(3));

console.log('\n  Jackknife variance estimation:');
const jackR = looResults.map(l => l.r);
const jackMean = jackR.reduce((a, b) => a + b, 0) / jackR.length;
const jackVar = ((jackR.length - 1) / jackR.length) * jackR.reduce((s, r) => s + (r - jackMean) ** 2, 0);
const jackSE = Math.sqrt(jackVar);
console.log('    Jackknife mean r: ' + jackMean.toFixed(3));
console.log('    Jackknife SE: ' + jackSE.toFixed(3));
console.log('    Jackknife 95% CI: [' + (jackMean - 1.96 * jackSE).toFixed(3) + ', ' + (jackMean + 1.96 * jackSE).toFixed(3) + ']');

console.log('\n  Bootstrap (N=5000):');
const nBoot = 5000;
const bootR = [];
for (let b = 0; b < nBoot; b++) {
  const sample = [];
  for (let i = 0; i < validBase.length; i++) {
    sample.push(validBase[Math.floor(rng() * validBase.length)]);
  }
  bootR.push(pearsonR(sample.map(s => s.DQ), sample.map(s => s.sa)));
}
bootR.sort((a, b) => a - b);
const ci95 = [bootR[Math.floor(nBoot * 0.025)], bootR[Math.floor(nBoot * 0.975)]];
const bootMean = bootR.reduce((a, b) => a + b, 0) / nBoot;
const bootPctPositive = (bootR.filter(r => r > 0).length / nBoot * 100).toFixed(1);
const bootPctAbove05 = (bootR.filter(r => r > 0.5).length / nBoot * 100).toFixed(1);
console.log('    Mean r: ' + bootMean.toFixed(3));
console.log('    95% CI: [' + ci95[0].toFixed(3) + ', ' + ci95[1].toFixed(3) + ']');
console.log('    % positive: ' + bootPctPositive + '%');
console.log('    % above 0.5: ' + bootPctAbove05 + '%');

let V3_verdict;
if (looAllPositive && parseFloat(bootPctPositive) > 90 && kendallTau > 0.3) V3_verdict = 'ROBUST — correlation survives all LOO drops, bootstrap, and pair analysis';
else if (looAllPositive && parseFloat(bootPctPositive) > 80) V3_verdict = 'STABLE — positive under LOO but wide uncertainty';
else V3_verdict = 'FRAGILE — not robust under resampling';

console.log('  V3 VERDICT: ' + V3_verdict);


console.log('\n' + '='.repeat(72));
console.log('V4: BAR / INNER-REGION EXCLUSION');
console.log('Does shapeAmplitude survive when inner region (potential bar) is excluded?');
console.log('='.repeat(72));

const innerExclusionFracs = [0.0, 0.1, 0.15, 0.2, 0.25, 0.3];
const barResults = [];

for (const frac of innerExclusionFracs) {
  const amps = [], dqs = [];
  for (const name of galNames) {
    const sa = computeShapeAmplitude(fitsData[name], maxMode, 12, 16, 0.90, frac);
    if (sa !== null) { amps.push(sa); dqs.push(galDQ[name]); }
  }
  const r = pearsonR(dqs, amps);
  barResults.push({ innerExcluded: (frac * 100).toFixed(0) + '%', r, n: dqs.length });
  console.log('  Exclude inner ' + (frac * 100).toFixed(0) + '%: r=' + r.toFixed(3) + ' (N=' + dqs.length + ')');
}

const barRValues = barResults.map(b => b.r);
const barAllPositive = barRValues.every(r => r > 0.3);
const barMeanR = barRValues.reduce((a, b) => a + b, 0) / barRValues.length;
const noBarR = barRValues[0];
const maxExclR = barRValues[barRValues.length - 1];

console.log('\n  All exclusion levels r > 0.3: ' + barAllPositive);
console.log('  No exclusion r: ' + noBarR.toFixed(3));
console.log('  Max exclusion (30%) r: ' + maxExclR.toFixed(3));
console.log('  Mean r: ' + barMeanR.toFixed(3));

const barDriven = noBarR - maxExclR > 0.3;
console.log('  Bar-driven (r drops > 0.3 with exclusion): ' + barDriven);

let V4_verdict;
if (barAllPositive && !barDriven) V4_verdict = 'ROBUST — shapeAmplitude NOT driven by inner bar; survives exclusion';
else if (!barDriven) V4_verdict = 'STABLE — some weakening but not bar-dominated';
else V4_verdict = 'BAR-CONTAMINATED — inner region dominates shapeAmplitude';

console.log('  V4 VERDICT: ' + V4_verdict);


console.log('\n' + '='.repeat(72));
console.log('COMBINED DM-4V VERDICT');
console.log('='.repeat(72));

const verdicts = [
  { test: 'V1', name: 'Binning sensitivity', verdict: V1_verdict },
  { test: 'V2', name: 'Normalization sensitivity', verdict: V2_verdict },
  { test: 'V3', name: 'LOO + pair robustness', verdict: V3_verdict },
  { test: 'V4', name: 'Bar / inner exclusion', verdict: V4_verdict },
];

for (const v of verdicts) { console.log('  ' + v.test + ' (' + v.name + '): ' + v.verdict); }

const nRobust = verdicts.filter(v => v.verdict.startsWith('ROBUST')).length;
const nStable = verdicts.filter(v => v.verdict.startsWith('STABLE')).length;
const nFragile = verdicts.filter(v => !v.verdict.startsWith('ROBUST') && !v.verdict.startsWith('STABLE')).length;

console.log('\n  ROBUST: ' + nRobust + '/4');
console.log('  STABLE: ' + nStable + '/4');
console.log('  FRAGILE/SENSITIVE: ' + nFragile + '/4');

let overallVerdict;
if (nRobust >= 3) overallVerdict = 'shapeAmplitude CONFIRMED — robust quantitative carrier of H. Ready for DM-5.';
else if (nRobust + nStable >= 3) overallVerdict = 'shapeAmplitude CONFIRMED with caveats — stable but wide uncertainty (N=7). Proceed to DM-5 with caution.';
else if (nFragile >= 2) overallVerdict = 'shapeAmplitude FRAGILE — result may be artifact of analysis choices. DM-4 needs rethinking.';
else overallVerdict = 'shapeAmplitude MIXED — partially robust. Consider composite parameter.';

console.log('\n  OVERALL: ' + overallVerdict);

console.log('\n  Key numbers for claim:');
console.log('    DM-4 base r(DQ, shapeAmp): ' + baseR.toFixed(3));
console.log('    Binning range: [' + binRMin.toFixed(3) + ', ' + binRMax.toFixed(3) + ']');
console.log('    LOO range: [' + looMin.toFixed(3) + ', ' + looMax.toFixed(3) + ']');
console.log('    Kendall tau: ' + kendallTau.toFixed(3));
console.log('    Bootstrap 95% CI: [' + ci95[0].toFixed(3) + ', ' + ci95[1].toFixed(3) + ']');
console.log('    Bootstrap % positive: ' + bootPctPositive + '%');
console.log('    Bar exclusion stable: ' + !barDriven);


const output = {
  program: 12,
  phase: 'DM-4V',
  title: 'DM-4 Verification: shapeAmplitude stress tests',
  timestamp: new Date().toISOString(),
  N: validBase.length,
  galaxies: validBase.map(g => g.name),
  baseCorrelation: baseR,

  V1_binning: {
    configs: binResults,
    meanR: binRMean, minR: binRMin, maxR: binRMax, range: binRRange,
    allAbove05: binAllPositive,
    verdict: V1_verdict,
  },

  V2_normalization: {
    modes: normResults,
    allPositive: allNormPositive,
    bestMode: bestNorm.mode,
    verdict: V2_verdict,
  },

  V3_robustness: {
    loo: { results: looResults, allPositive: looAllPositive, min: looMin, max: looMax, nAbove05: nLOOAbove05 },
    kendallTau, concordant, discordant,
    jackknife: { mean: jackMean, se: jackSE, ci95: [jackMean - 1.96 * jackSE, jackMean + 1.96 * jackSE] },
    bootstrap: { mean: bootMean, ci95, pctPositive: parseFloat(bootPctPositive), pctAbove05: parseFloat(bootPctAbove05), nBoot },
    verdict: V3_verdict,
  },

  V4_barExclusion: {
    results: barResults,
    allPositive: barAllPositive,
    barDriven,
    verdict: V4_verdict,
  },

  combined: {
    nRobust, nStable, nFragile,
    overall: overallVerdict,
  },

  keyNumbers: {
    baseR, binRange: [binRMin, binRMax],
    looRange: [looMin, looMax],
    kendallTau,
    bootstrapCI95: ci95,
    bootstrapPctPositive: parseFloat(bootPctPositive),
    barExclStable: !barDriven,
  },
};

const outDir = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results');
fs.mkdirSync(outDir, { recursive: true });
const outPath = path.join(outDir, 'program12-DM4V-verification.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nResult saved: ' + outPath);
console.log('\n=== DM-4V COMPLETE ===');
