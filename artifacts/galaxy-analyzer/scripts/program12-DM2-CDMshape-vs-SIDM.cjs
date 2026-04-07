#!/usr/bin/env node
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

function permutationP(x, y, nPerm, observed) {
  let count = 0;
  const yCopy = y.slice();
  for (let p = 0; p < nPerm; p++) {
    for (let i = yCopy.length - 1; i > 0; i--) { const j = Math.floor(rng() * (i + 1)); [yCopy[i], yCopy[j]] = [yCopy[j], yCopy[i]]; }
    const rPerm = Math.abs(pearsonR(x, yCopy));
    if (rPerm >= Math.abs(observed)) count++;
  }
  return count / nPerm;
}

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


function computeRadialMSpectrum(fits, maxMode) {
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
  let maxR = 0; for (const r of radii) if (r > maxR) maxR = r;

  const nRadBins = 10;
  const azBins = 16;
  const binW = maxR / nRadBins;

  const innerZone = [];
  const outerZone = [];
  const allBins = [];

  for (let rb = 0; rb < nRadBins; rb++) {
    const rMin = rb * binW, rMax = (rb + 1) * binW;
    const rMid = (rMin + rMax) / 2;
    const rNorm = rMid / maxR;

    const azVals = Array.from({ length: azBins }, () => []);
    for (let k = 0; k < coords.length; k++) {
      if (radii[k] < rMin || radii[k] >= rMax) continue;
      const theta = Math.atan2(coords[k].y, coords[k].x);
      const azIdx = Math.floor(((theta + Math.PI) / (2 * Math.PI)) * azBins) % azBins;
      azVals[azIdx].push(vrel[k]);
    }
    const azMeans = azVals.map(a => a.length > 0 ? a.reduce((s, v) => s + v, 0) / a.length : NaN);
    const valid = azMeans.filter(v => !isNaN(v));
    if (valid.length < azBins * 0.5) continue;
    const mean = valid.reduce((a, b) => a + b, 0) / valid.length;

    const mPower = Array(maxMode + 1).fill(0);
    let cnt = 0;
    for (let m = 0; m <= maxMode; m++) {
      let cm = 0, sm = 0;
      cnt = 0;
      for (let a = 0; a < azBins; a++) {
        if (isNaN(azMeans[a])) continue;
        const ang = (a + 0.5) * 2 * Math.PI / azBins;
        const dv = azMeans[a] - mean;
        if (m === 0) cm += Math.abs(dv);
        else { cm += dv * Math.cos(m * ang); sm += dv * Math.sin(m * ang); }
        cnt++;
      }
      mPower[m] = m === 0 ? (cnt > 0 ? cm / cnt : 0) : Math.sqrt(cm ** 2 + sm ** 2) / Math.max(cnt, 1);
    }

    const totalPower = mPower.reduce((a, b) => a + b, 0);
    const nonRotPower = mPower.slice(2).reduce((a, b) => a + b, 0);
    const angularComplexity = totalPower > 0 ? nonRotPower / totalPower : 0;

    const binData = { rNorm, mPower, totalPower, nonRotPower, angularComplexity, nPixels: azVals.reduce((s, a) => s + a.length, 0) };
    allBins.push(binData);

    if (rNorm <= 0.5) innerZone.push(binData);
    else outerZone.push(binData);
  }

  if (allBins.length < 3) return null;

  const innerC = innerZone.length > 0 ? innerZone.reduce((s, b) => s + b.angularComplexity, 0) / innerZone.length : 0;
  const outerC = outerZone.length > 0 ? outerZone.reduce((s, b) => s + b.angularComplexity, 0) / outerZone.length : 0;
  const suppressionRatio = outerC > 0 ? innerC / outerC : NaN;

  const globalMPower = Array(maxMode + 1).fill(0);
  for (const b of allBins) for (let m = 0; m <= maxMode; m++) globalMPower[m] += b.mPower[m];
  const totalGlobal = globalMPower.reduce((a, b) => a + b, 0);
  const globalNonRot = globalMPower.slice(2).reduce((a, b) => a + b, 0);
  const globalAngComplexity = totalGlobal > 0 ? globalNonRot / totalGlobal : 0;

  const innerMPower = Array(maxMode + 1).fill(0);
  for (const b of innerZone) for (let m = 0; m <= maxMode; m++) innerMPower[m] += b.mPower[m];
  const outerMPower = Array(maxMode + 1).fill(0);
  for (const b of outerZone) for (let m = 0; m <= maxMode; m++) outerMPower[m] += b.mPower[m];

  return {
    innerC, outerC, suppressionRatio,
    globalAngComplexity,
    globalMPower,
    innerMPower,
    outerMPower,
    nBins: allBins.length,
    nInner: innerZone.length,
    nOuter: outerZone.length,
    radialProfile: allBins,
  };
}


console.log('='.repeat(72));
console.log('PROGRAM 12 / DM-2: CDM+SHAPE vs SIDM — DECISIVE KILL TEST');
console.log('Question: Is H from non-axisymmetric halo shape (CDM) or core isothermalisation (SIDM)?');
console.log('='.repeat(72));

const dataDir = path.join(__dirname, '..', 'data', 'things-2d');
const thingsGals = [
  { name: 'NGC2841' }, { name: 'NGC5055' }, { name: 'NGC3521' },
  { name: 'NGC6946' }, { name: 'NGC7331' }, { name: 'NGC2403' },
  { name: 'NGC2903' }, { name: 'NGC3198' }, { name: 'NGC4826' },
  { name: 'NGC4736' }, { name: 'DDO154' }, { name: 'IC2574' },
];

const maxMode = 6;
const analyzed = [];

for (const tg of thingsGals) {
  const fitsPath = path.join(dataDir, tg.name + '_MOM1.FITS');
  if (!fs.existsSync(fitsPath)) continue;
  if (fs.statSync(fitsPath).size < 1000) continue;

  const galData = gals.find(g => normalize(g.name) === normalize(tg.name));
  if (!galData) { console.log('  SKIP ' + tg.name + ' (not in SPARC)'); continue; }

  const fits = parseFITS(fitsPath);
  const spec = computeRadialMSpectrum(fits, maxMode);
  if (!spec) { console.log('  SKIP ' + tg.name + ' (too few bins)'); continue; }

  const sr = resultsMap[normalize(tg.name)];
  const haloConc = sr && sr.models && sr.models.dark_halo_linear ? sr.models.dark_halo_linear.mse : NaN;

  analyzed.push({
    name: tg.name,
    DQ: galData.DQ,
    Vflat: galData.Vflat,
    hR: galData.hR,
    ...spec,
    haloConc,
  });
}

console.log('\n  Galaxies analyzed: ' + analyzed.length);
if (analyzed.length < 5) { console.error('  Too few galaxies for meaningful analysis'); process.exit(1); }

for (const g of analyzed) {
  console.log('\n  ' + g.name + ' (DQ=' + g.DQ.toFixed(2) + ', Vflat=' + g.Vflat + ')');
  console.log('    Inner C = ' + g.innerC.toFixed(4) + '  Outer C = ' + g.outerC.toFixed(4) + '  Suppression = ' + g.suppressionRatio.toFixed(3));
  console.log('    Global angular complexity = ' + g.globalAngComplexity.toFixed(4));
}


console.log('\n' + '='.repeat(72));
console.log('TEST T1: NORMALIZED INNER SUPPRESSION');
console.log('SIDM predicts: systematic inner suppression (suppressionRatio << 1)');
console.log('CDM+shape predicts: no systematic inner suppression');
console.log('='.repeat(72));

const suppRatios = analyzed.map(g => g.suppressionRatio).filter(r => !isNaN(r));
const meanSuppression = suppRatios.reduce((a, b) => a + b, 0) / suppRatios.length;
const nSuppressed = suppRatios.filter(r => r < 0.8).length;
const nEnhanced = suppRatios.filter(r => r > 1.2).length;
const rDQ_supp = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.suppressionRatio));
const pDQ_supp = permutationP(analyzed.map(g => g.DQ), analyzed.map(g => g.suppressionRatio), 10000, rDQ_supp);

console.log('\n  Mean suppression ratio (inner/outer): ' + meanSuppression.toFixed(3));
console.log('  N with inner < outer (ratio < 0.8): ' + nSuppressed + '/' + analyzed.length);
console.log('  N with inner > outer (ratio > 1.2): ' + nEnhanced + '/' + analyzed.length);
console.log('  r(DQ, suppression ratio): ' + rDQ_supp.toFixed(3) + ' (p = ' + pDQ_supp.toFixed(4) + ')');

let T1_verdict, T1_favors;
if (meanSuppression < 0.6 && nSuppressed >= analyzed.length * 0.7) {
  T1_verdict = 'STRONG INNER SUPPRESSION — SIDM signal detected';
  T1_favors = 'SIDM';
} else if (meanSuppression < 0.8 && nSuppressed >= analyzed.length * 0.5) {
  T1_verdict = 'MODERATE INNER SUPPRESSION — weak SIDM signal';
  T1_favors = 'SIDM-LEAN';
} else if (meanSuppression > 0.8 && meanSuppression < 1.2) {
  T1_verdict = 'NO SYSTEMATIC SUPPRESSION — consistent with CDM+shape';
  T1_favors = 'CDM+shape';
} else {
  T1_verdict = 'MIXED — no clear pattern';
  T1_favors = 'INCONCLUSIVE';
}
console.log('\n  T1 VERDICT: ' + T1_verdict);
console.log('  T1 FAVORS: ' + T1_favors);


console.log('\n' + '='.repeat(72));
console.log('TEST T2: CORE-LINKED ANGULAR COMPLEXITY');
console.log('CDM+shape predicts: angular complexity correlates with halo shape state');
console.log('SIDM predicts: angular complexity correlates with core suppression');
console.log('='.repeat(72));

const rDQ_angC = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.globalAngComplexity));
const pDQ_angC = permutationP(analyzed.map(g => g.DQ), analyzed.map(g => g.globalAngComplexity), 10000, rDQ_angC);

const rDQ_outerC = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.outerC));
const rDQ_innerC = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.innerC));

const rhR_angC = pearsonR(analyzed.map(g => g.hR), analyzed.map(g => g.globalAngComplexity));

console.log('\n  r(DQ, global angular complexity): ' + rDQ_angC.toFixed(3) + ' (p = ' + pDQ_angC.toFixed(4) + ')');
console.log('  r(DQ, outer complexity): ' + rDQ_outerC.toFixed(3));
console.log('  r(DQ, inner complexity): ' + rDQ_innerC.toFixed(3));
console.log('  r(haloResponse, angular complexity): ' + rhR_angC.toFixed(3));
console.log('  Outer > Inner correlation with DQ: ' + (Math.abs(rDQ_outerC) > Math.abs(rDQ_innerC) ? 'YES — shape-led' : 'NO — core-led'));

let T2_verdict, T2_favors;
if (Math.abs(rDQ_outerC) > Math.abs(rDQ_innerC) + 0.1) {
  T2_verdict = 'OUTER-LED — angular complexity driven by outer halo shape';
  T2_favors = 'CDM+shape';
} else if (Math.abs(rDQ_innerC) > Math.abs(rDQ_outerC) + 0.1) {
  T2_verdict = 'INNER-LED — angular complexity driven by core region';
  T2_favors = 'SIDM';
} else {
  T2_verdict = 'NO CLEAR RADIAL PREFERENCE';
  T2_favors = 'INCONCLUSIVE';
}
console.log('\n  T2 VERDICT: ' + T2_verdict);
console.log('  T2 FAVORS: ' + T2_favors);


console.log('\n' + '='.repeat(72));
console.log('TEST T3: INACCESSIBILITY-STRENGTH PARADOX COMPATIBILITY');
console.log('CDM+shape: shape is 2D, destroyed by azimuthal average — strong + hidden');
console.log('SIDM: core changes are radial, partially visible in 1D — harder to hide');
console.log('='.repeat(72));

const strengthProxy = analyzed.map(g => Math.abs(g.DQ));
const hiddenProxy = analyzed.map(g => 1 - g.globalAngComplexity);

const rStr_DQ = pearsonR(strengthProxy, analyzed.map(g => g.globalAngComplexity));
const outerFracs = analyzed.map(g => {
  const outerTotal = g.outerMPower.reduce((a, b) => a + b, 0);
  const innerTotal = g.innerMPower.reduce((a, b) => a + b, 0);
  return outerTotal / (outerTotal + innerTotal + 1e-10);
});
const rDQ_outerFrac = pearsonR(analyzed.map(g => g.DQ), outerFracs);

const nonRotInner = analyzed.map(g => g.innerMPower.slice(2).reduce((a, b) => a + b, 0));
const nonRotOuter = analyzed.map(g => g.outerMPower.slice(2).reduce((a, b) => a + b, 0));
const outerDominanceFrac = analyzed.filter((g, i) => nonRotOuter[i] > nonRotInner[i]).length / analyzed.length;

console.log('\n  r(|DQ|, angular complexity): ' + rStr_DQ.toFixed(3));
console.log('  r(DQ, outer power fraction): ' + rDQ_outerFrac.toFixed(3));
console.log('  Fraction with outer > inner non-rot power: ' + (outerDominanceFrac * 100).toFixed(0) + '%');
console.log('  (SIDM would predict inner dominant for high-DQ if core rounding drives H)');

let T3_verdict, T3_favors;
if (outerDominanceFrac >= 0.7) {
  T3_verdict = 'OUTER-DOMINANT — non-rotational power concentrated outside, hidden from 1D';
  T3_favors = 'CDM+shape';
} else if (outerDominanceFrac <= 0.3) {
  T3_verdict = 'INNER-DOMINANT — non-rotational power concentrated in core';
  T3_favors = 'SIDM';
} else {
  T3_verdict = 'MIXED DISTRIBUTION';
  T3_favors = 'INCONCLUSIVE';
}
console.log('\n  T3 VERDICT: ' + T3_verdict);
console.log('  T3 FAVORS: ' + T3_favors);


console.log('\n' + '='.repeat(72));
console.log('TEST T4: MATCHED-PAIR DECISIVE CHECK');
console.log('Compare high-H vs low-H galaxies: shape-complexity or core-suppression?');
console.log('='.repeat(72));

const sortedByDQ = analyzed.slice().sort((a, b) => b.DQ - a.DQ);
const highH = sortedByDQ.slice(0, 3);
const lowH = sortedByDQ.slice(-3);

console.log('\n  High-H galaxies: ' + highH.map(g => g.name + '(' + g.DQ.toFixed(2) + ')').join(', '));
console.log('  Low-H galaxies: ' + lowH.map(g => g.name + '(' + g.DQ.toFixed(2) + ')').join(', '));

const highH_outerC = highH.reduce((s, g) => s + g.outerC, 0) / highH.length;
const lowH_outerC = lowH.reduce((s, g) => s + g.outerC, 0) / lowH.length;
const highH_innerC = highH.reduce((s, g) => s + g.innerC, 0) / highH.length;
const lowH_innerC = lowH.reduce((s, g) => s + g.innerC, 0) / lowH.length;
const highH_supp = highH.reduce((s, g) => s + g.suppressionRatio, 0) / highH.length;
const lowH_supp = lowH.reduce((s, g) => s + g.suppressionRatio, 0) / lowH.length;
const highH_globalC = highH.reduce((s, g) => s + g.globalAngComplexity, 0) / highH.length;
const lowH_globalC = lowH.reduce((s, g) => s + g.globalAngComplexity, 0) / lowH.length;

console.log('\n  High-H mean outer C: ' + highH_outerC.toFixed(4));
console.log('  Low-H  mean outer C: ' + lowH_outerC.toFixed(4));
console.log('  Ratio (high/low outer C): ' + (highH_outerC / lowH_outerC).toFixed(2) + 'x');

console.log('\n  High-H mean inner C: ' + highH_innerC.toFixed(4));
console.log('  Low-H  mean inner C: ' + lowH_innerC.toFixed(4));
console.log('  Ratio (high/low inner C): ' + (highH_innerC / lowH_innerC).toFixed(2) + 'x');

console.log('\n  High-H mean suppression: ' + highH_supp.toFixed(3));
console.log('  Low-H  mean suppression: ' + lowH_supp.toFixed(3));

console.log('\n  High-H mean global C: ' + highH_globalC.toFixed(4));
console.log('  Low-H  mean global C: ' + lowH_globalC.toFixed(4));

const outerCDiff = highH_outerC - lowH_outerC;
const innerCDiff = highH_innerC - lowH_innerC;
const shapeLed = Math.abs(outerCDiff) > Math.abs(innerCDiff);
const corePattern = highH_supp < lowH_supp * 0.8;

let T4_verdict, T4_favors;
if (shapeLed && !corePattern) {
  T4_verdict = 'SHAPE-LED — high-H have more outer complexity, no core suppression pattern';
  T4_favors = 'CDM+shape';
} else if (!shapeLed && corePattern) {
  T4_verdict = 'CORE-LED — high-H show systematic core suppression';
  T4_favors = 'SIDM';
} else if (shapeLed && corePattern) {
  T4_verdict = 'MIXED — both outer shape and core effects present';
  T4_favors = 'INCONCLUSIVE';
} else {
  T4_verdict = 'NO CLEAR PATTERN';
  T4_favors = 'INCONCLUSIVE';
}
console.log('\n  Shape-led (outer diff > inner diff): ' + shapeLed);
console.log('  Core suppression pattern: ' + corePattern);
console.log('\n  T4 VERDICT: ' + T4_verdict);
console.log('  T4 FAVORS: ' + T4_favors);


console.log('\n' + '='.repeat(72));
console.log('COMBINED VERDICT TABLE');
console.log('='.repeat(72));

const tests = [
  { test: 'T1', name: 'Inner suppression', verdict: T1_verdict, favors: T1_favors },
  { test: 'T2', name: 'Core-linked complexity', verdict: T2_verdict, favors: T2_favors },
  { test: 'T3', name: 'Hiddenness paradox', verdict: T3_verdict, favors: T3_favors },
  { test: 'T4', name: 'Matched pairs', verdict: T4_verdict, favors: T4_favors },
];

for (const t of tests) {
  console.log('  ' + t.test + ' (' + t.name + '): ' + t.favors);
}

const cdmWins = tests.filter(t => t.favors === 'CDM+shape').length;
const sidmWins = tests.filter(t => t.favors === 'SIDM' || t.favors === 'SIDM-LEAN').length;
const inconclusive = tests.filter(t => t.favors === 'INCONCLUSIVE').length;

console.log('\n  CDM+shape wins: ' + cdmWins + '/4');
console.log('  SIDM wins: ' + sidmWins + '/4');
console.log('  Inconclusive: ' + inconclusive + '/4');

let overallVerdict;
if (cdmWins >= 3) overallVerdict = 'CDM+shape is strongly favored over SIDM';
else if (cdmWins === 2 && sidmWins <= 1) overallVerdict = 'CDM+shape is moderately favored over SIDM';
else if (sidmWins >= 3) overallVerdict = 'SIDM is strongly favored over CDM+shape';
else if (sidmWins === 2 && cdmWins <= 1) overallVerdict = 'SIDM is moderately favored over CDM+shape';
else overallVerdict = 'No decisive discrimination between CDM+shape and SIDM';

console.log('\n  OVERALL: ' + overallVerdict);


const output = {
  program: 12,
  phase: 'DM-2',
  title: 'CDM+Shape vs SIDM — Decisive Kill Test',
  timestamp: new Date().toISOString(),
  N: analyzed.length,
  galaxies: analyzed.map(g => g.name),

  T1: {
    name: 'Normalized inner suppression',
    meanSuppressionRatio: meanSuppression,
    nSuppressed,
    nEnhanced,
    rDQ_suppression: rDQ_supp,
    pVal: pDQ_supp,
    verdict: T1_verdict,
    favors: T1_favors,
    perGalaxy: analyzed.map(g => ({ name: g.name, DQ: g.DQ, innerC: g.innerC, outerC: g.outerC, suppressionRatio: g.suppressionRatio })),
  },

  T2: {
    name: 'Core-linked angular complexity',
    rDQ_globalAngComplexity: rDQ_angC,
    pVal: pDQ_angC,
    rDQ_outerC: rDQ_outerC,
    rDQ_innerC: rDQ_innerC,
    rhR_angComplexity: rhR_angC,
    outerLed: Math.abs(rDQ_outerC) > Math.abs(rDQ_innerC),
    verdict: T2_verdict,
    favors: T2_favors,
  },

  T3: {
    name: 'Inaccessibility-strength paradox compatibility',
    rStrength_angComplexity: rStr_DQ,
    rDQ_outerFraction: rDQ_outerFrac,
    outerDominanceFraction: outerDominanceFrac,
    verdict: T3_verdict,
    favors: T3_favors,
  },

  T4: {
    name: 'Matched-pair decisive check',
    highH: highH.map(g => g.name),
    lowH: lowH.map(g => g.name),
    highH_outerC,
    lowH_outerC,
    highH_innerC,
    lowH_innerC,
    highH_suppression: highH_supp,
    lowH_suppression: lowH_supp,
    shapeLed,
    corePattern,
    verdict: T4_verdict,
    favors: T4_favors,
  },

  combinedVerdict: {
    cdmWins,
    sidmWins,
    inconclusive,
    overall: overallVerdict,
  },

  methodologyNotes: [
    'All angular complexity is NORMALIZED: C = sum(P_m2..m6) / P_total — not raw power',
    'Inner zone: r/R_max <= 0.5; Outer zone: r/R_max > 0.5',
    '16-bin azimuthal decomposition (Nyquist safe to m=7)',
    '10 radial bins per galaxy',
    'Suppression ratio = inner_C / outer_C (SIDM predicts < 0.6)',
    'No reliance on m=2 alone; full m=2..6 spectrum used',
  ],
};

const outDir = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results');
fs.mkdirSync(outDir, { recursive: true });
const outPath = path.join(outDir, 'program12-DM2-CDMshape-vs-SIDM.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nResult saved: ' + outPath);
console.log('\n=== DM-2 COMPLETE ===');
