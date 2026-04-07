#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

function seededRNG(seed) { let s = seed | 0; return function() { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; }; }
const rng = seededRNG(20260412);

function normalize(n) { return n.replace(/\s+/g, '').replace(/^NGC0*/, 'NGC').replace(/^DDO0*/, 'DDO').replace(/^IC0*/, 'IC').replace(/^UGC0*/, 'UGC').replace(/^PGC0*/, 'PGC').toUpperCase(); }
function pearsonR(x, y) { const n = x.length; if (n < 3) return NaN; const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n; let sxy = 0, sxx = 0, syy = 0; for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; } return (sxx > 0 && syy > 0) ? sxy / Math.sqrt(sxx * syy) : 0; }
function ols(X, y) { const n = y.length, p = X[0].length; const my = y.reduce((a, b) => a + b, 0) / n; const mX = Array(p).fill(0); for (let j = 0; j < p; j++) { for (let i = 0; i < n; i++) mX[j] += X[i][j]; mX[j] /= n; } const Xc = X.map(r => r.map((v, j) => v - mX[j])); const yc = y.map(v => v - my); const XtX = Array.from({ length: p }, () => Array(p).fill(0)), Xty = Array(p).fill(0); for (let i = 0; i < n; i++) for (let j = 0; j < p; j++) { Xty[j] += Xc[i][j] * yc[i]; for (let k = 0; k < p; k++) XtX[j][k] += Xc[i][j] * Xc[i][k]; } const aug = XtX.map((r, i) => [...r, Xty[i]]); for (let c = 0; c < p; c++) { let mr = c; for (let r = c + 1; r < p; r++) if (Math.abs(aug[r][c]) > Math.abs(aug[mr][c])) mr = r; [aug[c], aug[mr]] = [aug[mr], aug[c]]; if (Math.abs(aug[c][c]) < 1e-14) continue; for (let r = c + 1; r < p; r++) { const f = aug[r][c] / aug[c][c]; for (let j = c; j <= p; j++) aug[r][j] -= f * aug[c][j]; } } const beta = Array(p).fill(0); for (let i = p - 1; i >= 0; i--) { beta[i] = aug[i][p]; for (let j = i + 1; j < p; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= Math.abs(aug[i][i]) > 1e-14 ? aug[i][i] : 1; } const res = []; for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < p; j++) pred += beta[j] * Xc[i][j]; res.push(y[i] - pred); } return { beta, residuals: res }; }
function zscore(arr) { const m = arr.reduce((a, b) => a + b, 0) / arr.length; const s = Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / arr.length); return s > 1e-10 ? arr.map(v => (v - m) / s) : arr.map(() => 0); }
function R2(residuals, y) { const my = y.reduce((a, b) => a + b, 0) / y.length; const ssTot = y.reduce((s, v) => s + (v - my) ** 2, 0); const ssRes = residuals.reduce((s, v) => s + v ** 2, 0); return ssTot > 0 ? 1 - ssRes / ssTot : 0; }

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

console.log('=== DM-5C: Extra Physics or Overfitting? ===');
console.log('Question: Is the +26.6pp R² improvement from extra parameters');
console.log('  (outerSupport, downstreamQuietness) a real structural signal,');
console.log('  or overfitting on N=7?\n');

const thingsDir = path.join(__dirname, '..', 'data', 'things-2d');
const thingsGalaxies = ['NGC2841', 'NGC5055', 'NGC3521', 'NGC7331', 'NGC2403', 'NGC2903', 'NGC3198'];
const maxMode = 6, nRadBins = 12, azBins = 16;

const galData = [];
for (const gName of thingsGalaxies) {
  const fitsPath = path.join(thingsDir, gName + '_MOM1.FITS');
  if (!fs.existsSync(fitsPath)) continue;
  const fits = parseFITS(fitsPath);
  const nx = fits.dims[0], ny = fits.dims[1];
  const crpix1 = fits.header['CRPIX1'] || nx / 2, crpix2 = fits.header['CRPIX2'] || ny / 2;
  const cdelt1 = fits.header['CDELT1'] || 1, cdelt2 = fits.header['CDELT2'] || 1;
  const vel = [], coords = [];
  for (let j = 0; j < ny; j++) { for (let i = 0; i < nx; i++) { const v = fits.data[j * nx + i]; if (isNaN(v)) continue; vel.push(v); coords.push({ x: (i + 1 - crpix1) * cdelt1 * 3600, y: (j + 1 - crpix2) * cdelt2 * 3600 }); } }
  if (vel.length < 100) continue;
  const sorted = vel.slice().sort((a, b) => a - b);
  const vsys = sorted[Math.floor(sorted.length / 2)];
  const vrel = vel.map(v => v - vsys);
  const radii = coords.map(c => Math.sqrt(c.x ** 2 + c.y ** 2));
  const rSorted = radii.slice().sort((a, b) => a - b);
  const rMax = rSorted[Math.floor(rSorted.length * 0.90)];
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
    }
    const nonAxiPower = mAmp.slice(2).reduce((a, b) => a + b, 0);
    bins.push({ rNorm, mAmp, nonAxiPower });
  }
  if (bins.length < 4) continue;

  const shapeAmplitude = bins.reduce((s, b) => s + b.nonAxiPower, 0) / bins.length;
  const outerBins = bins.filter(b => b.rNorm > 0.6);
  const outerSupport = outerBins.length > 0 ? outerBins.reduce((s, b) => s + b.nonAxiPower, 0) / outerBins.length : 0;
  const binPowers = bins.map(b => b.nonAxiPower);
  const bpMean = binPowers.reduce((a, b) => a + b, 0) / binPowers.length;
  const bpStd = Math.sqrt(binPowers.reduce((s, v) => s + (v - bpMean) ** 2, 0) / binPowers.length);
  const radialUniformity = bpMean > 0 ? 1 - bpStd / bpMean : 0;

  const sp = sparcMap[normalize(gName)];
  const sr = resultsMap[normalize(gName)];
  const rcResidualStd = sr && sr.models && sr.models.dark_halo_linear ? Math.sqrt(sr.models.dark_halo_linear.mse) : 0;
  const downstreamQuietness = rcResidualStd > 0 ? shapeAmplitude / rcResidualStd : 0;

  const galInfo = gals.find(g => normalize(g.name) === normalize(gName));

  galData.push({
    name: gName,
    DQ: galInfo ? galInfo.DQ : 0,
    shapeAmplitude,
    outerSupport,
    radialUniformity,
    downstreamQuietness,
    rcResidualStd,
  });
}

const DQs = galData.map(g => g.DQ);
const sAmps = galData.map(g => g.shapeAmplitude);
const oSups = galData.map(g => g.outerSupport);
const rUnis = galData.map(g => g.radialUniformity);
const dQuis = galData.map(g => g.downstreamQuietness);
const nG = galData.length;

console.log('Per-galaxy data:');
for (const g of galData) {
  console.log('  ' + g.name + ': DQ=' + g.DQ.toFixed(3) + ', sA=' + g.shapeAmplitude.toFixed(0) +
    ', oS=' + g.outerSupport.toFixed(0) + ', rU=' + g.radialUniformity.toFixed(3) +
    ', dQ=' + g.downstreamQuietness.toFixed(1));
}

console.log('\n\n========================================');
console.log('  TEST C1: INCREMENTAL LOO');
console.log('========================================');
console.log('Compare out-of-sample prediction error across 3 models.\n');

const models = [
  { label: 'M1: shapeAmplitude only', buildX: (g) => [g.shapeAmplitude] },
  { label: 'M2: shapeAmp + outerSupport', buildX: (g) => [g.shapeAmplitude, g.outerSupport] },
  { label: 'M3: shapeAmp + outerSup + quietness', buildX: (g) => [g.shapeAmplitude, g.outerSupport, g.downstreamQuietness] },
];

const looResults = [];
for (const model of models) {
  const looPredErrors = [];
  for (let i = 0; i < nG; i++) {
    const trainIdx = Array.from({ length: nG }, (_, j) => j).filter(j => j !== i);
    const trainX = trainIdx.map(j => model.buildX(galData[j]));
    const trainY = trainIdx.map(j => DQs[j]);
    const fit = ols(trainX, trainY);
    const testX = model.buildX(galData[i]);
    const mTrainY = trainY.reduce((a, b) => a + b, 0) / trainY.length;
    const mTrainX = trainX[0].map((_, c) => trainX.reduce((s, r) => s + r[c], 0) / trainX.length);
    let pred = mTrainY;
    for (let c = 0; c < testX.length; c++) pred += fit.beta[c] * (testX[c] - mTrainX[c]);
    const err = DQs[i] - pred;
    looPredErrors.push({ name: galData[i].name, actual: DQs[i], predicted: +pred.toFixed(3), error: +err.toFixed(3) });
  }
  const mse = looPredErrors.reduce((s, e) => s + e.error ** 2, 0) / nG;
  const rmse = Math.sqrt(mse);

  const fullX = galData.map(g => model.buildX(g));
  const fullFit = ols(fullX, DQs);
  const fullR2 = R2(fullFit.residuals, DQs);

  const looR2 = 1 - looPredErrors.reduce((s, e) => s + e.error ** 2, 0) /
    DQs.reduce((s, v) => s + (v - DQs.reduce((a, b) => a + b, 0) / nG) ** 2, 0);

  console.log('  ' + model.label);
  console.log('    In-sample R²: ' + fullR2.toFixed(3));
  console.log('    LOO R² (out-of-sample): ' + looR2.toFixed(3));
  console.log('    LOO RMSE: ' + rmse.toFixed(3));
  console.log('    R² drop (in→out): ' + ((fullR2 - looR2) * 100).toFixed(1) + ' pp');
  for (const e of looPredErrors) {
    console.log('      Drop ' + e.name + ': actual=' + e.actual.toFixed(3) + ', pred=' + e.predicted + ', err=' + e.error);
  }
  console.log('');
  looResults.push({ model: model.label, fullR2: +fullR2.toFixed(3), looR2: +looR2.toFixed(3), rmse: +rmse.toFixed(3), drop: +((fullR2 - looR2) * 100).toFixed(1), looPredictions: looPredErrors });
}

const m1LooR2 = looResults[0].looR2;
const m2LooR2 = looResults[1].looR2;
const m3LooR2 = looResults[2].looR2;
const bestLoo = Math.max(m1LooR2, m2LooR2, m3LooR2);
const bestModel = bestLoo === m1LooR2 ? 'M1' : bestLoo === m2LooR2 ? 'M2' : 'M3';

console.log('  LOO Summary:');
console.log('    M1 LOO R²: ' + m1LooR2.toFixed(3));
console.log('    M2 LOO R²: ' + m2LooR2.toFixed(3));
console.log('    M3 LOO R²: ' + m3LooR2.toFixed(3));
console.log('    Best out-of-sample: ' + bestModel);

const c1Overfitting = (looResults[2].fullR2 - looResults[2].looR2) > 0.3;
const c1ExtraHelps = m2LooR2 > m1LooR2 + 0.05 || m3LooR2 > m1LooR2 + 0.05;
const c1Verdict = c1Overfitting && !c1ExtraHelps ? 'OVERFITTING'
  : c1ExtraHelps && !c1Overfitting ? 'EXTRAS_GENUINE'
  : c1ExtraHelps ? 'EXTRAS_PARTIALLY_GENUINE'
  : 'SHAPE_ALONE_OPTIMAL';
console.log('  C1 VERDICT: ' + c1Verdict);

console.log('\n\n========================================');
console.log('  TEST C2: PERMUTATION NECESSITY');
console.log('========================================');
console.log('Fix shapeAmplitude, shuffle each extra parameter 5000 times.');
console.log('If R² collapses → that parameter is necessary.\n');

const nPerm = 5000;
const fullR2_M3 = looResults[2].fullR2;

function permTest(label, buildXreal, buildXshuffled, shuffleIdx) {
  const realX = galData.map(g => buildXreal(g));
  const realFit = ols(realX, DQs);
  const realR2 = R2(realFit.residuals, DQs);

  let countBetter = 0;
  const shuffledR2s = [];
  for (let p = 0; p < nPerm; p++) {
    const permVals = shuffleIdx.slice();
    for (let i = permVals.length - 1; i > 0; i--) {
      const j = Math.floor(rng() * (i + 1));
      [permVals[i], permVals[j]] = [permVals[j], permVals[i]];
    }
    const shuffX = galData.map((g, gi) => buildXshuffled(g, permVals[gi]));
    const shuffFit = ols(shuffX, DQs);
    const sR2 = R2(shuffFit.residuals, DQs);
    shuffledR2s.push(sR2);
    if (sR2 >= realR2) countBetter++;
  }
  const pVal = countBetter / nPerm;
  const meanShuffR2 = shuffledR2s.reduce((a, b) => a + b, 0) / nPerm;
  const r2Drop = realR2 - meanShuffR2;

  console.log('  ' + label);
  console.log('    Real R²: ' + realR2.toFixed(3));
  console.log('    Mean shuffled R²: ' + meanShuffR2.toFixed(3));
  console.log('    R² drop when shuffled: ' + r2Drop.toFixed(3));
  console.log('    p-value: ' + pVal.toFixed(4));
  console.log('    ' + (pVal < 0.1 ? 'NECESSARY (shuffling hurts)' : 'COSMETIC (shuffling doesnt hurt)'));

  return { label, realR2: +realR2.toFixed(3), meanShuffR2: +meanShuffR2.toFixed(3), r2Drop: +r2Drop.toFixed(3), pVal: +pVal.toFixed(4), necessary: pVal < 0.1 };
}

const oSupIdx = Array.from({ length: nG }, (_, i) => i);
const dQuiIdx = Array.from({ length: nG }, (_, i) => i);

const permOuterSupport = permTest(
  'Shuffle outerSupport (keep shapeAmp + quietness fixed)',
  g => [g.shapeAmplitude, g.outerSupport, g.downstreamQuietness],
  (g, pi) => [g.shapeAmplitude, galData[pi].outerSupport, g.downstreamQuietness],
  oSupIdx
);

const permQuietness = permTest(
  'Shuffle downstreamQuietness (keep shapeAmp + outerSupport fixed)',
  g => [g.shapeAmplitude, g.outerSupport, g.downstreamQuietness],
  (g, pi) => [g.shapeAmplitude, g.outerSupport, galData[pi].downstreamQuietness],
  dQuiIdx
);

const permShapeAmp = permTest(
  'Shuffle shapeAmplitude (keep outerSupport + quietness fixed)',
  g => [g.shapeAmplitude, g.outerSupport, g.downstreamQuietness],
  (g, pi) => [galData[pi].shapeAmplitude, g.outerSupport, g.downstreamQuietness],
  oSupIdx
);

const c2Verdict = permShapeAmp.necessary && !permOuterSupport.necessary && !permQuietness.necessary
  ? 'SHAPE_ONLY_NECESSARY'
  : permShapeAmp.necessary && (permOuterSupport.necessary || permQuietness.necessary)
  ? 'MULTIPLE_NECESSARY'
  : 'UNCLEAR';
console.log('\n  C2 VERDICT: ' + c2Verdict);

console.log('\n\n========================================');
console.log('  TEST C3: CAUSAL ORDERING');
console.log('========================================');
console.log('Test mediation: does outerSupport/quietness mediate shape→DQ,');
console.log('or is it the reverse (shape→DQ, extras are downstream)?\n');

const rShapeDQ = pearsonR(sAmps, DQs);
const rOuterDQ = pearsonR(oSups, DQs);
const rQuietDQ = pearsonR(dQuis, DQs);
const rShapeOuter = pearsonR(sAmps, oSups);
const rShapeQuiet = pearsonR(sAmps, dQuis);
const rOuterQuiet = pearsonR(oSups, dQuis);

console.log('  Pairwise correlations:');
console.log('    r(shape, DQ) = ' + rShapeDQ.toFixed(3));
console.log('    r(outer, DQ) = ' + rOuterDQ.toFixed(3));
console.log('    r(quiet, DQ) = ' + rQuietDQ.toFixed(3));
console.log('    r(shape, outer) = ' + rShapeOuter.toFixed(3));
console.log('    r(shape, quiet) = ' + rShapeQuiet.toFixed(3));
console.log('    r(outer, quiet) = ' + rOuterQuiet.toFixed(3));

const partialOuterDQ = (rOuterDQ - rShapeDQ * rShapeOuter) /
  (Math.sqrt(1 - rShapeDQ ** 2) * Math.sqrt(1 - rShapeOuter ** 2));
const partialQuietDQ = (rQuietDQ - rShapeDQ * rShapeQuiet) /
  (Math.sqrt(1 - rShapeDQ ** 2) * Math.sqrt(1 - rShapeQuiet ** 2));
const partialShapeDQ_outer = (rShapeDQ - rOuterDQ * rShapeOuter) /
  (Math.sqrt(1 - rOuterDQ ** 2) * Math.sqrt(1 - rShapeOuter ** 2));
const partialShapeDQ_quiet = (rShapeDQ - rQuietDQ * rShapeQuiet) /
  (Math.sqrt(1 - rQuietDQ ** 2) * Math.sqrt(1 - rShapeQuiet ** 2));

console.log('\n  Partial correlations (controlling for shape):');
console.log('    r(outer, DQ | shape) = ' + partialOuterDQ.toFixed(3));
console.log('    r(quiet, DQ | shape) = ' + partialQuietDQ.toFixed(3));
console.log('  Partial correlations (controlling for extras):');
console.log('    r(shape, DQ | outer) = ' + partialShapeDQ_outer.toFixed(3));
console.log('    r(shape, DQ | quiet) = ' + partialShapeDQ_quiet.toFixed(3));

const outerIsMediator = Math.abs(partialOuterDQ) < 0.3 && Math.abs(rOuterDQ) > 0.3;
const quietIsMediator = Math.abs(partialQuietDQ) < 0.3 && Math.abs(rQuietDQ) > 0.3;
const outerIsIndependent = Math.abs(partialOuterDQ) > 0.3;
const quietIsIndependent = Math.abs(partialQuietDQ) > 0.3;

console.log('\n  Mediation analysis:');
console.log('    outerSupport: ' + (outerIsMediator ? 'DOWNSTREAM (mediated by shape)' : outerIsIndependent ? 'INDEPENDENT contribution' : 'UNCLEAR'));
console.log('    quietness: ' + (quietIsMediator ? 'DOWNSTREAM (mediated by shape)' : quietIsIndependent ? 'INDEPENDENT contribution' : 'UNCLEAR'));

console.log('\n  Causal models:');
console.log('    Hypothesis A: shape → outer → DQ');
console.log('      r(shape,outer) = ' + rShapeOuter.toFixed(3) + ', r(outer,DQ|shape) = ' + partialOuterDQ.toFixed(3));
console.log('      ' + (rShapeOuter > 0.5 && Math.abs(partialOuterDQ) < 0.3 ? 'SUPPORTED — outer mediates' : 'NOT SUPPORTED'));
console.log('    Hypothesis B: shape → quiet → DQ');
console.log('      r(shape,quiet) = ' + rShapeQuiet.toFixed(3) + ', r(quiet,DQ|shape) = ' + partialQuietDQ.toFixed(3));
console.log('      ' + (rShapeQuiet > 0.5 && Math.abs(partialQuietDQ) < 0.3 ? 'SUPPORTED — quiet mediates' : 'NOT SUPPORTED'));
console.log('    Hypothesis C: shape → DQ directly, extras are byproducts');
console.log('      r(shape,DQ|outer) = ' + partialShapeDQ_outer.toFixed(3) + ', r(shape,DQ|quiet) = ' + partialShapeDQ_quiet.toFixed(3));
console.log('      ' + (Math.abs(partialShapeDQ_outer) > 0.3 || Math.abs(partialShapeDQ_quiet) > 0.3 ? 'SUPPORTED — shape has direct effect beyond extras' : 'NOT SUPPORTED'));

const hypA = rShapeOuter > 0.5 && Math.abs(partialOuterDQ) < 0.3;
const hypB = rShapeQuiet > 0.5 && Math.abs(partialQuietDQ) < 0.3;
const hypC = Math.abs(partialShapeDQ_outer) > 0.3 || Math.abs(partialShapeDQ_quiet) > 0.3;

const c3Verdict = hypC && !hypA && !hypB ? 'SHAPE_DIRECT_EXTRAS_DOWNSTREAM'
  : hypC && (hypA || hypB) ? 'SHAPE_DIRECT_SOME_MEDIATION'
  : hypA || hypB ? 'MEDIATION_DOMINANT'
  : 'CAUSAL_STRUCTURE_UNCLEAR';
console.log('\n  C3 VERDICT: ' + c3Verdict);

console.log('\n\n========================================');
console.log('  TEST C4: MINIMAL SUFFICIENT MODEL');
console.log('========================================');
console.log('Identify the simplest model that works out-of-sample.\n');

const adjustedR2 = (r2, n, p) => 1 - (1 - r2) * (n - 1) / (n - p - 1);

const modelComparison = [
  { label: 'M1: shapeAmplitude', params: 1, fullR2: looResults[0].fullR2, looR2: looResults[0].looR2, rmse: looResults[0].rmse },
  { label: 'M2: shape + outer', params: 2, fullR2: looResults[1].fullR2, looR2: looResults[1].looR2, rmse: looResults[1].rmse },
  { label: 'M3: shape + outer + quiet', params: 3, fullR2: looResults[2].fullR2, looR2: looResults[2].looR2, rmse: looResults[2].rmse },
];

console.log('  Model Comparison:');
console.log('  ┌────────────────────────────┬──────┬───────────┬──────────┬──────────┬──────────┐');
console.log('  │ Model                      │ Pars │ In-samp R²│ Adj R²   │ LOO R²   │ LOO RMSE │');
console.log('  ├────────────────────────────┼──────┼───────────┼──────────┼──────────┼──────────┤');
for (const m of modelComparison) {
  const adjR2 = adjustedR2(m.fullR2, nG, m.params);
  console.log('  │ ' + m.label.padEnd(27) + '│ ' + String(m.params).padEnd(5) + '│ ' +
    m.fullR2.toFixed(3).padEnd(10) + '│ ' + adjR2.toFixed(3).padEnd(9) + '│ ' +
    m.looR2.toFixed(3).padEnd(9) + '│ ' + m.rmse.toFixed(3).padEnd(9) + '│');
}
console.log('  └────────────────────────────┴──────┴───────────┴──────────┴──────────┴──────────┘');

const bestLooModel = modelComparison.reduce((best, m) => m.looR2 > best.looR2 ? m : best);
const parsimoniousModel = modelComparison.reduce((best, m) => {
  if (m.looR2 > best.looR2 - 0.05 && m.params < best.params) return m;
  if (m.looR2 > best.looR2 + 0.05) return m;
  return best;
});

console.log('\n  Best LOO R²: ' + bestLooModel.label + ' (R²=' + bestLooModel.looR2.toFixed(3) + ')');
console.log('  Most parsimonious (within 5pp of best): ' + parsimoniousModel.label);

const c4Verdict = parsimoniousModel.params === 1 ? 'SHAPE_ALONE_SUFFICIENT'
  : parsimoniousModel.params === 2 ? 'SHAPE_PLUS_OUTER'
  : 'FULL_MODEL_NEEDED';
console.log('  C4 VERDICT: ' + c4Verdict);

console.log('\n\n========================================');
console.log('  DM-5C OVERALL RESULTS');
console.log('========================================\n');

const verdicts = { C1: c1Verdict, C2: c2Verdict, C3: c3Verdict, C4: c4Verdict };
for (const [k, v] of Object.entries(verdicts)) {
  console.log('  ' + k + ': ' + v);
}

let overallVerdict;
const shapeOnly = c4Verdict === 'SHAPE_ALONE_SUFFICIENT' || c1Verdict === 'SHAPE_ALONE_OPTIMAL';
const extrasOverfit = c1Verdict === 'OVERFITTING' && c2Verdict === 'SHAPE_ONLY_NECESSARY';
const extrasReal = c1Verdict.includes('GENUINE') && c2Verdict === 'MULTIPLE_NECESSARY';

if (shapeOnly || extrasOverfit) {
  overallVerdict = 'SHAPE_AMPLITUDE_IS_THE_WHOLE_STORY';
  console.log('\n  OVERALL: shapeAmplitude IS the whole story');
  console.log('  The extra parameters do not add genuine out-of-sample power.');
  console.log('  The +26.6pp was overfitting on N=7.');
  console.log('  CDM + complex halo shape (captured by a single parameter) is SUFFICIENT.');
} else if (extrasReal) {
  overallVerdict = 'SHAPE_PLUS_STRUCTURE_IS_THE_STORY';
  console.log('\n  OVERALL: shapeAmplitude + spatial structure is the story');
  console.log('  Extra parameters add genuine information beyond total amplitude.');
  console.log('  The halo shape has internal structure (inner/outer, uniformity)');
  console.log('  that independently contributes to H.');
} else {
  overallVerdict = 'SHAPE_DOMINANT_EXTRAS_AMBIGUOUS';
  console.log('\n  OVERALL: shapeAmplitude dominates, extras are ambiguous');
  console.log('  shapeAmplitude is clearly the primary driver (R²=0.63).');
  console.log('  Extra parameters may add marginal information, but with N=7');
  console.log('  we cannot distinguish genuine structure from overfitting.');
  console.log('  Conclusion: shapeAmplitude is the safe, parsimonious answer.');
  console.log('  Whether spatial distribution adds real physics remains open');
  console.log('  until N >> 7 (e.g., WALLABY, MeerKAT).');
}

const output = {
  program: 'DM-5C',
  title: 'Extra Physics or Overfitting?',
  date: new Date().toISOString().slice(0, 10),
  question: 'Is the +26.6pp R² from extra params genuine or overfitting?',
  tests: {
    C1_incrementalLOO: { models: looResults, bestOutOfSample: bestModel, verdict: c1Verdict },
    C2_permutationNecessity: {
      shapeAmplitude: permShapeAmp,
      outerSupport: permOuterSupport,
      quietness: permQuietness,
      verdict: c2Verdict,
    },
    C3_causalOrdering: {
      pairwise: { rShapeDQ: +rShapeDQ.toFixed(3), rOuterDQ: +rOuterDQ.toFixed(3), rQuietDQ: +rQuietDQ.toFixed(3), rShapeOuter: +rShapeOuter.toFixed(3), rShapeQuiet: +rShapeQuiet.toFixed(3) },
      partial: { rOuterDQ_shape: +partialOuterDQ.toFixed(3), rQuietDQ_shape: +partialQuietDQ.toFixed(3), rShapeDQ_outer: +partialShapeDQ_outer.toFixed(3), rShapeDQ_quiet: +partialShapeDQ_quiet.toFixed(3) },
      hypotheses: { A_outerMediates: hypA, B_quietMediates: hypB, C_shapeDirectEffect: hypC },
      verdict: c3Verdict,
    },
    C4_minimalSufficient: { comparison: modelComparison, bestLOO: bestLooModel.label, parsimonious: parsimoniousModel.label, verdict: c4Verdict },
  },
  verdicts,
  overallVerdict,
  perGalaxy: galData.map(g => ({
    name: g.name, DQ: +g.DQ.toFixed(3),
    shapeAmplitude: +g.shapeAmplitude.toFixed(0),
    outerSupport: +g.outerSupport.toFixed(0),
    radialUniformity: +g.radialUniformity.toFixed(3),
    downstreamQuietness: +g.downstreamQuietness.toFixed(1),
  })),
};

const outPath = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results', 'program12-DM5C-extra-physics-test.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nResult saved to: ' + outPath);
