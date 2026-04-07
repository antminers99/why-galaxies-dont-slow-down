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


function computeRadialM2(fits) {
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

  const radialBins = 20;
  const azBins = 12;
  const binW = maxR / radialBins;

  const binResults = [];

  for (let rb = 0; rb < radialBins; rb++) {
    const rMin = rb * binW, rMax = (rb + 1) * binW;
    const rMid = (rMin + rMax) / 2;
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
    if (valid.length < 4) { binResults.push({ rMid, m2: NaN, nPix: 0 }); continue; }
    const mean = valid.reduce((a, b) => a + b, 0) / valid.length;
    let c2 = 0, s2 = 0, cnt = 0;
    for (let a = 0; a < azBins; a++) {
      if (isNaN(azMeans[a])) continue;
      const ang = (a + 0.5) * 2 * Math.PI / azBins;
      const dv = azMeans[a] - mean;
      c2 += dv * Math.cos(2 * ang); s2 += dv * Math.sin(2 * ang);
      cnt++;
    }
    const m2 = Math.sqrt(c2 ** 2 + s2 ** 2) / cnt;
    const nPix = azVals.reduce((s, a) => s + a.length, 0);
    binResults.push({ rMid, m2, nPix });
  }

  const validBins = binResults.filter(b => !isNaN(b.m2) && b.nPix > 0);
  if (validBins.length < 4) return null;

  const halfIdx = Math.floor(validBins.length / 2);
  const innerBins = validBins.slice(0, halfIdx);
  const outerBins = validBins.slice(halfIdx);

  const m2_inner = innerBins.reduce((s, b) => s + b.m2, 0) / innerBins.length;
  const m2_outer = outerBins.reduce((s, b) => s + b.m2, 0) / outerBins.length;
  const ratio_outer_inner = m2_inner > 0 ? m2_outer / m2_inner : NaN;

  const m2_total = validBins.reduce((s, b) => s + b.m2, 0) / validBins.length;

  return {
    m2_inner,
    m2_outer,
    ratio_outer_inner,
    m2_total,
    nInnerBins: innerBins.length,
    nOuterBins: outerBins.length,
    radialProfile: validBins.map(b => ({ r: b.rMid, m2: b.m2, nPix: b.nPix })),
    maxR,
    nValid: vel.length,
  };
}


console.log('='.repeat(72));
console.log('PROGRAM 11 — TEST A: RADIAL m=2 GRADIENT');
console.log('CDM vs SIDM discriminator');
console.log('CDM predicts: m=2 roughly constant or increasing outward');
console.log('SIDM predicts: m=2 suppressed in core (inner << outer)');
console.log('='.repeat(72));

const dataDir = path.join(__dirname, '..', 'data', 'things-2d');
const thingsGals = [
  { name: 'NGC2841' }, { name: 'NGC5055' }, { name: 'NGC3521' },
  { name: 'NGC6946' }, { name: 'NGC7331' }, { name: 'NGC2403' },
  { name: 'NGC2903' }, { name: 'NGC3198' }, { name: 'NGC4826' },
  { name: 'NGC4736' }, { name: 'DDO154' }, { name: 'IC2574' },
];

const results = [];
for (const tg of thingsGals) {
  const fitsPath = path.join(dataDir, tg.name + '_MOM1.FITS');
  if (!fs.existsSync(fitsPath)) continue;
  const fStat = fs.statSync(fitsPath);
  if (fStat.size < 1000) continue;

  const galData = gals.find(g => normalize(g.name) === normalize(tg.name));
  if (!galData) { console.log('SKIP ' + tg.name + ' (not in SPARC)'); continue; }

  const fits = parseFITS(fitsPath);
  const res = computeRadialM2(fits);
  if (!res) { console.log('SKIP ' + tg.name + ' (too few valid bins)'); continue; }

  results.push({
    name: tg.name,
    DQ: galData.DQ,
    Vflat: galData.Vflat,
    inc: galData.inc,
    logMbar: galData.logMbar,
    ...res,
  });
}

console.log('\n  Galaxies analyzed: ' + results.length);


console.log('\n\n' + '#'.repeat(72));
console.log('A.1 — PER-GALAXY RADIAL m=2 GRADIENT');
console.log('#'.repeat(72));

console.log('\n  ' + 'Galaxy'.padEnd(12) + 'DQ'.padEnd(8) + 'Vflat'.padEnd(7) + 'm2_inner'.padEnd(10) + 'm2_outer'.padEnd(10) + 'outer/inner'.padEnd(13) + 'Pattern');
console.log('  ' + '-'.repeat(72));

for (const r of results) {
  let pattern = '';
  if (r.ratio_outer_inner > 1.5) pattern = 'OUTER-DOMINANT (SIDM-like?)';
  else if (r.ratio_outer_inner > 0.8) pattern = 'FLAT (CDM-like)';
  else if (r.ratio_outer_inner > 0.3) pattern = 'INNER-DOMINANT';
  else pattern = 'CORE-CONCENTRATED';

  console.log('  ' + r.name.padEnd(12) +
    r.DQ.toFixed(2).padEnd(8) +
    ('' + Math.round(r.Vflat)).padEnd(7) +
    r.m2_inner.toFixed(1).padEnd(10) +
    r.m2_outer.toFixed(1).padEnd(10) +
    r.ratio_outer_inner.toFixed(2).padEnd(13) +
    pattern);
}


console.log('\n\n' + '#'.repeat(72));
console.log('A.2 — CORRELATIONS');
console.log('#'.repeat(72));

const validR = results.filter(r => !isNaN(r.ratio_outer_inner) && !isNaN(r.DQ));
const dqArr = validR.map(r => r.DQ);
const ratioArr = validR.map(r => r.ratio_outer_inner);
const m2InnerArr = validR.map(r => r.m2_inner);
const m2OuterArr = validR.map(r => r.m2_outer);
const m2TotalArr = validR.map(r => r.m2_total);

const r_DQ_ratio = pearsonR(dqArr, ratioArr);
const r_DQ_inner = pearsonR(dqArr, m2InnerArr);
const r_DQ_outer = pearsonR(dqArr, m2OuterArr);
const r_DQ_total = pearsonR(dqArr, m2TotalArr);

console.log('\n  N = ' + validR.length);
console.log('  r(DQ, outer/inner ratio) = ' + r_DQ_ratio.toFixed(4));
console.log('  r(DQ, m2_inner)          = ' + r_DQ_inner.toFixed(4));
console.log('  r(DQ, m2_outer)          = ' + r_DQ_outer.toFixed(4));
console.log('  r(DQ, m2_total)          = ' + r_DQ_total.toFixed(4));


console.log('\n\n' + '#'.repeat(72));
console.log('A.3 — PERMUTATION TEST for r(DQ, outer/inner)');
console.log('#'.repeat(72));

const nPerm = 10000;
const obsR = r_DQ_ratio;
let countGE = 0;
for (let p = 0; p < nPerm; p++) {
  const shuffled = dqArr.slice();
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(rng() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  if (Math.abs(pearsonR(shuffled, ratioArr)) >= Math.abs(obsR)) countGE++;
}
const pVal = countGE / nPerm;
console.log('\n  Observed r = ' + obsR.toFixed(4));
console.log('  p-value (two-tailed, 10000 permutations) = ' + pVal.toFixed(4));
console.log('  ' + (pVal < 0.05 ? 'SIGNIFICANT at p < 0.05' : 'NOT significant at p < 0.05'));


console.log('\n\n' + '#'.repeat(72));
console.log('A.4 — RADIAL PROFILE COMPARISON: High-H vs Low-H');
console.log('#'.repeat(72));

const sortedByDQ = results.slice().sort((a, b) => b.DQ - a.DQ);
const nHalf = Math.floor(sortedByDQ.length / 2);
const highH = sortedByDQ.slice(0, nHalf);
const lowH = sortedByDQ.slice(sortedByDQ.length - nHalf);

console.log('\n  High-H group (' + highH.length + '): ' + highH.map(g => g.name).join(', '));
console.log('  Low-H group  (' + lowH.length + '): ' + lowH.map(g => g.name).join(', '));

const maxBins = Math.max(...results.map(r => r.radialProfile.length));
const nProfileBins = Math.min(maxBins, 15);

console.log('\n  Radial bin  High-H m2    Low-H m2    Diff    Interpretation');
console.log('  ' + '-'.repeat(65));

for (let b = 0; b < nProfileBins; b++) {
  const hVals = highH.map(g => g.radialProfile[b]).filter(Boolean).filter(p => !isNaN(p.m2));
  const lVals = lowH.map(g => g.radialProfile[b]).filter(Boolean).filter(p => !isNaN(p.m2));
  if (hVals.length === 0 || lVals.length === 0) continue;
  const hMean = hVals.reduce((s, p) => s + p.m2, 0) / hVals.length;
  const lMean = lVals.reduce((s, p) => s + p.m2, 0) / lVals.length;
  const diff = hMean - lMean;
  const zone = b < nProfileBins / 3 ? 'INNER' : b < 2 * nProfileBins / 3 ? 'MID' : 'OUTER';
  console.log('  ' + (b + 1 + '').padEnd(14) +
    hMean.toFixed(1).padEnd(13) +
    lMean.toFixed(1).padEnd(12) +
    ((diff >= 0 ? '+' : '') + diff.toFixed(1)).padEnd(8) +
    zone);
}

const hMeanRatio = highH.reduce((s, g) => s + g.ratio_outer_inner, 0) / highH.length;
const lMeanRatio = lowH.reduce((s, g) => s + g.ratio_outer_inner, 0) / lowH.length;
console.log('\n  Mean outer/inner ratio:');
console.log('    High-H: ' + hMeanRatio.toFixed(3));
console.log('    Low-H:  ' + lMeanRatio.toFixed(3));
console.log('    Difference: ' + (hMeanRatio - lMeanRatio).toFixed(3));


console.log('\n\n' + '#'.repeat(72));
console.log('A.5 — TEST A VERDICT');
console.log('#'.repeat(72));

const meanRatio = validR.reduce((s, r) => s + r.ratio_outer_inner, 0) / validR.length;
const nOuterDominant = validR.filter(r => r.ratio_outer_inner > 1.2).length;
const nFlat = validR.filter(r => r.ratio_outer_inner >= 0.8 && r.ratio_outer_inner <= 1.2).length;
const nInnerDominant = validR.filter(r => r.ratio_outer_inner < 0.8).length;

console.log('\n  SAMPLE STATISTICS:');
console.log('  Mean outer/inner ratio: ' + meanRatio.toFixed(3));
console.log('  Outer-dominant (>1.2):  ' + nOuterDominant + '/' + validR.length);
console.log('  Flat (0.8-1.2):         ' + nFlat + '/' + validR.length);
console.log('  Inner-dominant (<0.8):  ' + nInnerDominant + '/' + validR.length);

console.log('\n  DQ CORRELATION:');
console.log('  r(DQ, outer/inner) = ' + r_DQ_ratio.toFixed(4) + ', p = ' + pVal.toFixed(4));

let testA_verdict, testA_favors;
if (meanRatio > 1.3 && nOuterDominant > validR.length * 0.6) {
  testA_verdict = 'SIDM-FAVORED';
  testA_favors = 'SIDM';
  console.log('\n  VERDICT: m=2 is systematically suppressed in cores.');
  console.log('  Pattern is CONSISTENT with SIDM core isotropization.');
} else if (meanRatio < 0.8 && nInnerDominant > validR.length * 0.5) {
  testA_verdict = 'INNER-CONCENTRATED';
  testA_favors = 'Neither (baryonic?)';
  console.log('\n  VERDICT: m=2 is concentrated in inner regions.');
  console.log('  This suggests baryonic (bar/spiral) origin rather than halo shape.');
} else {
  testA_verdict = 'CDM-CONSISTENT';
  testA_favors = 'CDM';
  console.log('\n  VERDICT: m=2 gradient is roughly flat or mixed.');
  console.log('  Pattern is CONSISTENT with CDM triaxial halos (no core suppression).');
  console.log('  Does NOT show the systematic core suppression expected from SIDM.');
}

console.log('\n  INTERPRETATION FOR DM IDENTITY:');
if (Math.abs(r_DQ_ratio) > 0.5 && pVal < 0.1) {
  console.log('  The radial gradient CORRELATES with DQ — high-H galaxies');
  console.log('  show a different radial m=2 profile than low-H galaxies.');
  console.log('  This radial structure is a potential DM model discriminator.');
} else {
  console.log('  The radial gradient does NOT strongly correlate with DQ.');
  console.log('  The m=2 radial profile is similar across high-H and low-H galaxies.');
  console.log('  This means inner/outer gradient alone does not discriminate DM models');
  console.log('  through the H channel.');
}


const outDir = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results');
if (!fs.existsSync(outDir)) fs.mkdirSync(outDir, { recursive: true });

const outPath = path.join(outDir, 'testA-radial-gradient.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: 11,
  test: 'A',
  title: 'Radial m=2 Gradient (CDM vs SIDM)',
  timestamp: new Date().toISOString(),
  N: validR.length,
  correlations: {
    DQ_ratio: r_DQ_ratio,
    DQ_m2_inner: r_DQ_inner,
    DQ_m2_outer: r_DQ_outer,
    DQ_m2_total: r_DQ_total,
  },
  permutation: { nPerm, pVal_ratio: pVal },
  sampleStats: {
    meanRatio,
    nOuterDominant,
    nFlat,
    nInnerDominant,
  },
  groupComparison: {
    highH_meanRatio: hMeanRatio,
    lowH_meanRatio: lMeanRatio,
    highH_galaxies: highH.map(g => g.name),
    lowH_galaxies: lowH.map(g => g.name),
  },
  verdict: testA_verdict,
  favors: testA_favors,
  perGalaxy: results.map(r => ({
    name: r.name, DQ: r.DQ, Vflat: r.Vflat, inc: r.inc,
    m2_inner: r.m2_inner, m2_outer: r.m2_outer,
    ratio_outer_inner: r.ratio_outer_inner, m2_total: r.m2_total,
    radialProfile: r.radialProfile,
  })),
}, null, 2));
console.log('\nSaved: ' + outPath);
