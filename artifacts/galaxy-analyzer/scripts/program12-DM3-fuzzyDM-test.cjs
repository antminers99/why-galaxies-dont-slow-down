#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

function seededRNG(seed) { let s = seed | 0; return function() { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; }; }
const rng = seededRNG(20260408);

function normalize(n) { return n.replace(/\s+/g, '').replace(/^NGC0*/, 'NGC').replace(/^DDO0*/, 'DDO').replace(/^IC0*/, 'IC').replace(/^UGC0*/, 'UGC').replace(/^PGC0*/, 'PGC').toUpperCase(); }
function pearsonR(x, y) { const n = x.length; if (n < 3) return NaN; const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n; let sxy = 0, sxx = 0, syy = 0; for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; } return (sxx > 0 && syy > 0) ? sxy / Math.sqrt(sxx * syy) : 0; }
function ols(X, y) { const n = y.length, p = X[0].length; const my = y.reduce((a, b) => a + b, 0) / n; const mX = Array(p).fill(0); for (let j = 0; j < p; j++) { for (let i = 0; i < n; i++) mX[j] += X[i][j]; mX[j] /= n; } const Xc = X.map(r => r.map((v, j) => v - mX[j])); const yc = y.map(v => v - my); const XtX = Array.from({ length: p }, () => Array(p).fill(0)), Xty = Array(p).fill(0); for (let i = 0; i < n; i++) for (let j = 0; j < p; j++) { Xty[j] += Xc[i][j] * yc[i]; for (let k = 0; k < p; k++) XtX[j][k] += Xc[i][j] * Xc[i][k]; } const aug = XtX.map((r, i) => [...r, Xty[i]]); for (let c = 0; c < p; c++) { let mr = c; for (let r = c + 1; r < p; r++) if (Math.abs(aug[r][c]) > Math.abs(aug[mr][c])) mr = r; [aug[c], aug[mr]] = [aug[mr], aug[c]]; if (Math.abs(aug[c][c]) < 1e-14) continue; for (let r = c + 1; r < p; r++) { const f = aug[r][c] / aug[c][c]; for (let j = c; j <= p; j++) aug[r][j] -= f * aug[c][j]; } } const beta = Array(p).fill(0); for (let i = p - 1; i >= 0; i--) { beta[i] = aug[i][p]; for (let j = i + 1; j < p; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= Math.abs(aug[i][i]) > 1e-14 ? aug[i][i] : 1; } const res = []; for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < p; j++) pred += beta[j] * Xc[i][j]; res.push(y[i] - pred); } return { beta, residuals: res }; }
function zscore(arr) { const m = arr.reduce((a, b) => a + b, 0) / arr.length; const s = Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / arr.length); return s > 1e-10 ? arr.map(v => (v - m) / s) : arr.map(() => 0); }
function permutationP(x, y, nPerm, observed) { let count = 0; const yC = y.slice(); for (let p = 0; p < nPerm; p++) { for (let i = yC.length - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [yC[i], yC[j]] = [yC[j], yC[i]]; } if (Math.abs(pearsonR(x, yC)) >= Math.abs(observed)) count++; } return count / nPerm; }

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

function analyzeGalaxy(fits, maxMode) {
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
  const rMax = rSorted[Math.floor(rSorted.length * 0.90)];

  const nRadBins = 12;
  const azBins = 16;
  const binW = rMax / nRadBins;
  const bins = [];

  for (let rb = 0; rb < nRadBins; rb++) {
    const rMin = rb * binW, rMaxBin = (rb + 1) * binW;
    const rNorm = ((rMin + rMaxBin) / 2) / rMax;
    const azVals = Array.from({ length: azBins }, () => []);
    for (let k = 0; k < coords.length; k++) {
      if (radii[k] < rMin || radii[k] >= rMaxBin) continue;
      const theta = Math.atan2(coords[k].y, coords[k].x);
      const azIdx = Math.floor(((theta + Math.PI) / (2 * Math.PI)) * azBins) % azBins;
      azVals[azIdx].push(vrel[k]);
    }
    const azMeans = azVals.map(a => a.length > 0 ? a.reduce((s, v) => s + v, 0) / a.length : NaN);
    const valid = azMeans.filter(v => !isNaN(v));
    if (valid.length < azBins * 0.6) continue;
    const mean = valid.reduce((a, b) => a + b, 0) / valid.length;

    const mAmp = Array(maxMode + 1).fill(0);
    const mPhase = Array(maxMode + 1).fill(0);
    for (let m = 0; m <= maxMode; m++) {
      let cm = 0, sm = 0, cnt = 0;
      for (let a = 0; a < azBins; a++) {
        if (isNaN(azMeans[a])) continue;
        const ang = (a + 0.5) * 2 * Math.PI / azBins;
        const dv = azMeans[a] - mean;
        if (m === 0) cm += Math.abs(dv);
        else { cm += dv * Math.cos(m * ang); sm += dv * Math.sin(m * ang); }
        cnt++;
      }
      mAmp[m] = m === 0 ? (cnt > 0 ? cm / cnt : 0) : Math.sqrt(cm ** 2 + sm ** 2) / Math.max(cnt, 1);
      mPhase[m] = m > 0 ? Math.atan2(sm, cm) : 0;
    }

    const totalP = mAmp.reduce((a, b) => a + b, 0);
    const nonRotP = mAmp.slice(2).reduce((a, b) => a + b, 0);
    const oddP = mAmp[3] + mAmp[5];
    const evenP = mAmp[2] + mAmp[4] + mAmp[6];
    const angC = totalP > 0 ? nonRotP / totalP : 0;

    const nonRotFracs = [];
    for (let m = 2; m <= maxMode; m++) nonRotFracs.push(nonRotP > 0 ? mAmp[m] / nonRotP : 0);
    let entropy = 0;
    for (const f of nonRotFracs) { if (f > 1e-10) entropy -= f * Math.log2(f); }
    const maxEntropy = Math.log2(nonRotFracs.length);
    const normalizedEntropy = maxEntropy > 0 ? entropy / maxEntropy : 0;

    let nEffective = 0;
    if (nonRotP > 0) {
      let sumSq = 0;
      for (let m = 2; m <= maxMode; m++) { const f = mAmp[m] / nonRotP; sumSq += f * f; }
      nEffective = sumSq > 0 ? 1 / sumSq : 0;
    }

    const domFrac = nonRotP > 0 ? Math.max(...mAmp.slice(2)) / nonRotP : 0;
    const outsideDomFrac = 1 - domFrac;

    bins.push({
      rNorm, mAmp, mPhase, totalP, nonRotP, angC, oddP, evenP,
      oddEvenRatio: evenP > 0 ? oddP / evenP : Infinity,
      spectralEntropy: normalizedEntropy,
      nEffectiveModes: nEffective,
      outsideDomFrac,
      nPix: azVals.reduce((s, a) => s + a.length, 0),
    });
  }

  if (bins.length < 4) return null;

  const globalMeanEntropy = bins.reduce((s, b) => s + b.spectralEntropy, 0) / bins.length;
  const globalMeanNEff = bins.reduce((s, b) => s + b.nEffectiveModes, 0) / bins.length;
  const globalOutsideDom = bins.reduce((s, b) => s + b.outsideDomFrac, 0) / bins.length;
  const globalOddEven = bins.reduce((s, b) => s + (isFinite(b.oddEvenRatio) ? b.oddEvenRatio : 5), 0) / bins.length;

  const phaseFlips = [];
  for (let m = 2; m <= maxMode; m++) {
    let flips = 0, totalPairs = 0;
    for (let i = 1; i < bins.length; i++) {
      const dp = bins[i].mPhase[m] - bins[i - 1].mPhase[m];
      let wrapped = ((dp + Math.PI) % (2 * Math.PI)) - Math.PI;
      if (wrapped < -Math.PI) wrapped += 2 * Math.PI;
      if (Math.abs(wrapped) > Math.PI / 2) flips++;
      totalPairs++;
    }
    phaseFlips.push({ m, flips, totalPairs, flipRate: totalPairs > 0 ? flips / totalPairs : 0 });
  }

  const meanFlipRate = phaseFlips.reduce((s, f) => s + f.flipRate, 0) / phaseFlips.length;

  const phaseVariance = [];
  for (let m = 2; m <= maxMode; m++) {
    const phases = bins.map(b => b.mPhase[m]);
    if (phases.length < 3) { phaseVariance.push({ m, variance: 0, isCoherent: true }); continue; }
    const diffs = [];
    for (let i = 1; i < phases.length; i++) {
      let dp = phases[i] - phases[i - 1];
      dp = ((dp + Math.PI) % (2 * Math.PI)) - Math.PI;
      if (dp < -Math.PI) dp += 2 * Math.PI;
      diffs.push(dp);
    }
    const meanDiff = diffs.reduce((a, b) => a + b, 0) / diffs.length;
    const variance = diffs.reduce((s, d) => s + (d - meanDiff) ** 2, 0) / diffs.length;
    phaseVariance.push({ m, variance, isCoherent: variance < 1.0 });
  }

  const meanPhaseVariance = phaseVariance.reduce((s, p) => s + p.variance, 0) / phaseVariance.length;
  const nCoherentModes = phaseVariance.filter(p => p.isCoherent).length;

  const ampVariability = [];
  for (let m = 2; m <= maxMode; m++) {
    const amps = bins.map(b => b.mAmp[m]);
    const mean = amps.reduce((a, b) => a + b, 0) / amps.length;
    const std = Math.sqrt(amps.reduce((s, a) => s + (a - mean) ** 2, 0) / amps.length);
    const cv = mean > 0 ? std / mean : 0;
    ampVariability.push({ m, mean, std, cv });
  }
  const meanAmpCV = ampVariability.reduce((s, a) => s + a.cv, 0) / ampVariability.length;

  const tercile = Math.floor(bins.length / 3);
  const innerBins = bins.slice(0, tercile);
  const outerBins = bins.slice(tercile * 2);

  const innerOddEven = innerBins.reduce((s, b) => s + (isFinite(b.oddEvenRatio) ? b.oddEvenRatio : 5), 0) / Math.max(innerBins.length, 1);
  const outerOddEven = outerBins.reduce((s, b) => s + (isFinite(b.oddEvenRatio) ? b.oddEvenRatio : 5), 0) / Math.max(outerBins.length, 1);
  const innerEntropy = innerBins.reduce((s, b) => s + b.spectralEntropy, 0) / Math.max(innerBins.length, 1);
  const outerEntropy = outerBins.reduce((s, b) => s + b.spectralEntropy, 0) / Math.max(outerBins.length, 1);

  return {
    globalMeanEntropy, globalMeanNEff, globalOutsideDom, globalOddEven,
    phaseFlips, meanFlipRate,
    phaseVariance, meanPhaseVariance, nCoherentModes,
    ampVariability, meanAmpCV,
    innerOddEven, outerOddEven, innerEntropy, outerEntropy,
    nBins: bins.length,
    radialProfile: bins,
  };
}


console.log('='.repeat(72));
console.log('PROGRAM 12 / DM-3: FUZZY DM vs CDM+SHAPE');
console.log('Question: Does the angular spectrum carry a wave-like fingerprint,');
console.log('          or is it ordinary halo shape complexity?');
console.log('='.repeat(72));

const dataDir = path.join(__dirname, '..', 'data', 'things-2d');
const thingsGals = [
  { name: 'NGC2841' }, { name: 'NGC5055' }, { name: 'NGC3521' },
  { name: 'NGC6946' }, { name: 'NGC7331' }, { name: 'NGC2403' },
  { name: 'NGC2903' }, { name: 'NGC3198' },
];

const maxMode = 6;
const analyzed = [];

for (const tg of thingsGals) {
  const fitsPath = path.join(dataDir, tg.name + '_MOM1.FITS');
  if (!fs.existsSync(fitsPath) || fs.statSync(fitsPath).size < 1000) continue;
  const galData = gals.find(g => normalize(g.name) === normalize(tg.name));
  if (!galData) { console.log('  SKIP ' + tg.name); continue; }
  const fits = parseFITS(fitsPath);
  const res = analyzeGalaxy(fits, maxMode);
  if (!res) { console.log('  SKIP ' + tg.name + ' (too few bins)'); continue; }
  analyzed.push({ name: tg.name, DQ: galData.DQ, Vflat: galData.Vflat, ...res });
}

console.log('\n  Galaxies analyzed: ' + analyzed.length);


console.log('\n' + '='.repeat(72));
console.log('TEST T1: SPECTRAL WIDTH');
console.log('Fuzzy DM predicts: broad spectrum (high entropy, many effective modes)');
console.log('CDM+shape predicts: moderate complexity, dominated by few modes');
console.log('='.repeat(72));

for (const g of analyzed) {
  console.log('\n  ' + g.name + ' (DQ=' + g.DQ.toFixed(2) + '):');
  console.log('    Spectral entropy (normalized): ' + g.globalMeanEntropy.toFixed(3));
  console.log('    Effective modes: ' + g.globalMeanNEff.toFixed(2));
  console.log('    Fraction outside dominant: ' + g.globalOutsideDom.toFixed(3));
}

const rDQ_entropy = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.globalMeanEntropy));
const rDQ_nEff = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.globalMeanNEff));
const rDQ_outsideDom = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.globalOutsideDom));
const pDQ_entropy = permutationP(analyzed.map(g => g.DQ), analyzed.map(g => g.globalMeanEntropy), 10000, rDQ_entropy);

const meanEntropy = analyzed.reduce((s, g) => s + g.globalMeanEntropy, 0) / analyzed.length;
const meanNEff = analyzed.reduce((s, g) => s + g.globalMeanNEff, 0) / analyzed.length;

console.log('\n  Mean spectral entropy: ' + meanEntropy.toFixed(3) + ' (max=1.0, Fuzzy expects >0.85)');
console.log('  Mean effective modes: ' + meanNEff.toFixed(2) + ' (max=5, Fuzzy expects >3.5)');
console.log('  r(DQ, entropy): ' + rDQ_entropy.toFixed(3) + ' (p=' + pDQ_entropy.toFixed(4) + ')');
console.log('  r(DQ, N_eff): ' + rDQ_nEff.toFixed(3));
console.log('  r(DQ, outside-dominant frac): ' + rDQ_outsideDom.toFixed(3));

let T1_verdict, T1_favors;
if (meanEntropy > 0.85 && meanNEff > 3.5) {
  T1_verdict = 'BROAD SPECTRUM — wave-like distribution of power across modes';
  T1_favors = 'Fuzzy DM';
} else if (meanEntropy < 0.7 || meanNEff < 2.5) {
  T1_verdict = 'NARROW SPECTRUM — dominated by few modes, consistent with shape complexity';
  T1_favors = 'CDM+shape';
} else {
  T1_verdict = 'MODERATE SPECTRUM — broader than pure m=2 but not wave-like';
  T1_favors = 'CDM+shape (not distinctively wave-like)';
}
console.log('\n  T1 VERDICT: ' + T1_verdict);
console.log('  T1 FAVORS: ' + T1_favors);


console.log('\n' + '='.repeat(72));
console.log('TEST T2: RADIAL PHASE COHERENCE');
console.log('Fuzzy DM predicts: oscillatory phase flips (wave interference)');
console.log('CDM+shape predicts: smooth/stable phase with radius');
console.log('='.repeat(72));

for (const g of analyzed) {
  console.log('\n  ' + g.name + ' (DQ=' + g.DQ.toFixed(2) + '):');
  console.log('    Mean phase flip rate: ' + g.meanFlipRate.toFixed(3));
  console.log('    Mean phase variance: ' + g.meanPhaseVariance.toFixed(3));
  console.log('    Coherent modes: ' + g.nCoherentModes + '/5');
  console.log('    Mean amplitude CV: ' + g.meanAmpCV.toFixed(3));
  for (const pf of g.phaseFlips) {
    console.log('      m=' + pf.m + ': flips=' + pf.flips + '/' + pf.totalPairs + ' rate=' + pf.flipRate.toFixed(2));
  }
}

const rDQ_flipRate = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.meanFlipRate));
const rDQ_phaseVar = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.meanPhaseVariance));
const rDQ_ampCV = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.meanAmpCV));
const rDQ_coherent = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.nCoherentModes));

const meanFlipRate = analyzed.reduce((s, g) => s + g.meanFlipRate, 0) / analyzed.length;
const meanPhaseVar = analyzed.reduce((s, g) => s + g.meanPhaseVariance, 0) / analyzed.length;
const meanCoherent = analyzed.reduce((s, g) => s + g.nCoherentModes, 0) / analyzed.length;

console.log('\n  Mean flip rate: ' + meanFlipRate.toFixed(3) + ' (Fuzzy expects >0.4)');
console.log('  Mean phase variance: ' + meanPhaseVar.toFixed(3) + ' (Fuzzy expects >1.5)');
console.log('  Mean coherent modes: ' + meanCoherent.toFixed(1) + '/5 (Fuzzy expects <2)');
console.log('  r(DQ, flip rate): ' + rDQ_flipRate.toFixed(3));
console.log('  r(DQ, phase variance): ' + rDQ_phaseVar.toFixed(3));
console.log('  r(DQ, amplitude CV): ' + rDQ_ampCV.toFixed(3));
console.log('  r(DQ, coherent modes): ' + rDQ_coherent.toFixed(3));

let T2_verdict, T2_favors;
if (meanFlipRate > 0.4 && meanPhaseVar > 1.5 && meanCoherent < 2) {
  T2_verdict = 'OSCILLATORY — wave-like phase structure detected';
  T2_favors = 'Fuzzy DM';
} else if (meanFlipRate < 0.25 && meanCoherent >= 3) {
  T2_verdict = 'SMOOTH/COHERENT — stable halo shape, no wave signature';
  T2_favors = 'CDM+shape';
} else {
  T2_verdict = 'INTERMEDIATE — some phase variation but not distinctively wave-like';
  T2_favors = 'CDM+shape (no clear wave fingerprint)';
}
console.log('\n  T2 VERDICT: ' + T2_verdict);
console.log('  T2 FAVORS: ' + T2_favors);


console.log('\n' + '='.repeat(72));
console.log('TEST T3: ODD/EVEN MODE ASYMMETRY');
console.log('Fuzzy DM predicts: wave patterns distribute to all modes, no parity bias');
console.log('CDM+shape predicts: odd dominance from spiral/warp structure');
console.log('='.repeat(72));

for (const g of analyzed) {
  console.log('\n  ' + g.name + ' (DQ=' + g.DQ.toFixed(2) + '):');
  console.log('    Global odd/even ratio: ' + g.globalOddEven.toFixed(2));
  console.log('    Inner odd/even: ' + g.innerOddEven.toFixed(2));
  console.log('    Outer odd/even: ' + g.outerOddEven.toFixed(2));
}

const rDQ_oddEven = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.globalOddEven));
const meanOddEven = analyzed.reduce((s, g) => s + g.globalOddEven, 0) / analyzed.length;
const nOddDominant = analyzed.filter(g => g.globalOddEven > 1.5).length;

console.log('\n  Mean odd/even ratio: ' + meanOddEven.toFixed(2) + ' (=1 for Fuzzy, >2 for CDM+shape)');
console.log('  N with strong odd dominance (ratio>1.5): ' + nOddDominant + '/' + analyzed.length);
console.log('  r(DQ, odd/even): ' + rDQ_oddEven.toFixed(3));

const innerMeanOE = analyzed.reduce((s, g) => s + g.innerOddEven, 0) / analyzed.length;
const outerMeanOE = analyzed.reduce((s, g) => s + g.outerOddEven, 0) / analyzed.length;
console.log('  Inner mean odd/even: ' + innerMeanOE.toFixed(2));
console.log('  Outer mean odd/even: ' + outerMeanOE.toFixed(2));

let T3_verdict, T3_favors;
if (meanOddEven > 2.0 && nOddDominant >= analyzed.length * 0.7) {
  T3_verdict = 'STRONG ODD DOMINANCE — parity asymmetry inconsistent with wave interference';
  T3_favors = 'CDM+shape';
} else if (meanOddEven < 1.3) {
  T3_verdict = 'BALANCED PARITY — consistent with wave-like mode distribution';
  T3_favors = 'Fuzzy DM';
} else {
  T3_verdict = 'MODERATE ODD BIAS — some parity asymmetry';
  T3_favors = 'CDM+shape (odd dominance argues against pure wave model)';
}
console.log('\n  T3 VERDICT: ' + T3_verdict);
console.log('  T3 FAVORS: ' + T3_favors);


console.log('\n' + '='.repeat(72));
console.log('TEST T4: HIDDENNESS PARADOX COMPATIBILITY');
console.log('CDM+shape: shape is 2D, hidden from 1D — naturally compatible');
console.log('Fuzzy DM: wave interference is 3D — but does it produce BILATERAL coupling?');
console.log('='.repeat(72));

const rDQ_angC_global = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => {
  const allC = g.radialProfile.reduce((s, b) => s + b.angC, 0) / g.radialProfile.length;
  return allC;
}));
const pDQ_angC = permutationP(
  analyzed.map(g => g.DQ),
  analyzed.map(g => g.radialProfile.reduce((s, b) => s + b.angC, 0) / g.radialProfile.length),
  10000, rDQ_angC_global);

console.log('\n  r(DQ, angular complexity): ' + rDQ_angC_global.toFixed(3) + ' (p=' + pDQ_angC.toFixed(4) + ')');

const waveUniqueness = meanFlipRate > 0.3 && meanEntropy > 0.8 && meanOddEven < 1.5;
console.log('\n  Wave-unique fingerprint present: ' + waveUniqueness);
console.log('  (requires: flip rate > 0.3 AND entropy > 0.8 AND odd/even < 1.5)');

let T4_verdict, T4_favors;
if (waveUniqueness) {
  T4_verdict = 'WAVE FINGERPRINT DETECTED — Fuzzy DM has unique signature beyond CDM+shape';
  T4_favors = 'Fuzzy DM';
} else {
  T4_verdict = 'NO WAVE FINGERPRINT — all observed features explainable by halo shape complexity';
  T4_favors = 'CDM+shape';
}
console.log('\n  T4 VERDICT: ' + T4_verdict);
console.log('  T4 FAVORS: ' + T4_favors);


console.log('\n' + '='.repeat(72));
console.log('COMBINED VERDICT TABLE');
console.log('='.repeat(72));

const tests = [
  { test: 'T1', name: 'Spectral width', verdict: T1_verdict, favors: T1_favors },
  { test: 'T2', name: 'Radial phase coherence', verdict: T2_verdict, favors: T2_favors },
  { test: 'T3', name: 'Odd/even asymmetry', verdict: T3_verdict, favors: T3_favors },
  { test: 'T4', name: 'Hiddenness + wave uniqueness', verdict: T4_verdict, favors: T4_favors },
];

for (const t of tests) { console.log('  ' + t.test + ' (' + t.name + '): ' + t.favors); }

const cdmWins = tests.filter(t => t.favors.startsWith('CDM')).length;
const fuzzyWins = tests.filter(t => t.favors.startsWith('Fuzzy')).length;
console.log('\n  CDM+shape: ' + cdmWins + '/4');
console.log('  Fuzzy DM: ' + fuzzyWins + '/4');

let overallVerdict;
if (cdmWins >= 3) overallVerdict = 'Fuzzy DM shows NO unique wave signature. CDM+shape explains all observed features. Fuzzy DM is effectively DEAD for this dataset.';
else if (fuzzyWins >= 3) overallVerdict = 'Fuzzy DM shows unique wave fingerprint not explainable by CDM+shape. Fuzzy DM SURVIVES as viable alternative.';
else if (cdmWins === 2 && fuzzyWins <= 1) overallVerdict = 'CDM+shape moderately favored. Fuzzy DM weakened but not killed.';
else overallVerdict = 'No decisive discrimination. Both models remain viable.';

console.log('\n  OVERALL: ' + overallVerdict);


const output = {
  program: 12,
  phase: 'DM-3',
  title: 'Fuzzy DM vs CDM+Shape — Wave Fingerprint Test',
  timestamp: new Date().toISOString(),
  N: analyzed.length,
  galaxies: analyzed.map(g => g.name),

  T1: {
    name: 'Spectral width',
    meanEntropy, meanNEff, rDQ_entropy, pDQ_entropy, rDQ_nEff, rDQ_outsideDom,
    verdict: T1_verdict, favors: T1_favors,
  },

  T2: {
    name: 'Radial phase coherence',
    meanFlipRate, meanPhaseVar, meanCoherent: Math.round(meanCoherent * 10) / 10,
    rDQ_flipRate, rDQ_phaseVar, rDQ_ampCV, rDQ_coherent,
    verdict: T2_verdict, favors: T2_favors,
  },

  T3: {
    name: 'Odd/even mode asymmetry',
    meanOddEven, nOddDominant, rDQ_oddEven,
    innerMeanOE, outerMeanOE,
    verdict: T3_verdict, favors: T3_favors,
  },

  T4: {
    name: 'Hiddenness + wave uniqueness',
    rDQ_angComplexity: rDQ_angC_global,
    pVal: pDQ_angC,
    waveUniqueness,
    verdict: T4_verdict, favors: T4_favors,
  },

  combinedVerdict: { cdmWins, fuzzyWins, overall: overallVerdict },

  perGalaxy: analyzed.map(g => ({
    name: g.name, DQ: g.DQ, Vflat: g.Vflat,
    spectralEntropy: g.globalMeanEntropy,
    nEffectiveModes: g.globalMeanNEff,
    outsideDomFrac: g.globalOutsideDom,
    oddEvenRatio: g.globalOddEven,
    meanFlipRate: g.meanFlipRate,
    meanPhaseVariance: g.meanPhaseVariance,
    nCoherentModes: g.nCoherentModes,
    meanAmpCV: g.meanAmpCV,
  })),
};

const outDir = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results');
fs.mkdirSync(outDir, { recursive: true });
const outPath = path.join(outDir, 'program12-DM3-fuzzyDM-test.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nResult saved: ' + outPath);
console.log('\n=== DM-3 COMPLETE ===');
