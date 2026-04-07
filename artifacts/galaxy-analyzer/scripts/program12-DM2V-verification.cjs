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
  const headerBlocks = []; let offset = 0, headerDone = false;
  while (!headerDone && offset < buf.length) { const block = buf.slice(offset, offset + 2880).toString('ascii'); offset += 2880; for (let i = 0; i < 36; i++) { const card = block.substring(i * 80, (i + 1) * 80); headerBlocks.push(card); if (card.startsWith('END')) { headerDone = true; break; } } }
  const header = {};
  for (const card of headerBlocks) { if (card.startsWith('END')) break; if (!card.includes('=')) continue; const key = card.substring(0, 8).trim(); let valStr = card.substring(10, 80); const ci = valStr.indexOf('/'); if (ci >= 0 && !valStr.trimStart().startsWith("'")) valStr = valStr.substring(0, ci); valStr = valStr.trim().replace(/'/g, '').trim(); const nv = parseFloat(valStr); header[key] = isNaN(nv) ? valStr : nv; }
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

  const rSorted = radii.slice().sort((a, b) => a - b);
  const r90 = rSorted[Math.floor(rSorted.length * 0.90)];
  const rMax = r90;

  const nRadBins = 10;
  const azBins = 16;
  const binW = rMax / nRadBins;

  const bins = [];
  for (let rb = 0; rb < nRadBins; rb++) {
    const rMin = rb * binW, rMaxBin = (rb + 1) * binW;
    const rMid = (rMin + rMaxBin) / 2;
    const rNorm = rMid / rMax;

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

    const mPower = Array(maxMode + 1).fill(0);
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
      mPower[m] = m === 0 ? (cnt > 0 ? cm / cnt : 0) : Math.sqrt(cm ** 2 + sm ** 2) / Math.max(cnt, 1);
    }

    const totalP = mPower.reduce((a, b) => a + b, 0);
    const nonRotP = mPower.slice(2).reduce((a, b) => a + b, 0);
    const oddP = mPower[3] + mPower[5];
    const evenP = mPower[2] + mPower[4] + mPower[6];
    const angC = totalP > 0 ? nonRotP / totalP : 0;
    const oddFrac = nonRotP > 0 ? oddP / nonRotP : 0;
    const evenFrac = nonRotP > 0 ? evenP / nonRotP : 0;
    const nPix = azVals.reduce((s, a) => s + a.length, 0);

    bins.push({ rNorm, mPower, totalP, nonRotP, angC, oddP, evenP, oddFrac, evenFrac, nPix });
  }

  if (bins.length < 4) return null;

  const tercile = Math.floor(bins.length / 3);
  const innerBins = bins.slice(0, tercile);
  const midBins = bins.slice(tercile, tercile * 2);
  const outerBins = bins.slice(tercile * 2);

  const zoneStats = (zoneBins) => {
    if (zoneBins.length === 0) return { angC: 0, oddFrac: 0, evenFrac: 0, nonRotP: 0, nBins: 0 };
    const avgAngC = zoneBins.reduce((s, b) => s + b.angC, 0) / zoneBins.length;
    const avgOdd = zoneBins.reduce((s, b) => s + b.oddFrac, 0) / zoneBins.length;
    const avgEven = zoneBins.reduce((s, b) => s + b.evenFrac, 0) / zoneBins.length;
    const totalNonRot = zoneBins.reduce((s, b) => s + b.nonRotP, 0);
    return { angC: avgAngC, oddFrac: avgOdd, evenFrac: avgEven, nonRotP: totalNonRot, nBins: zoneBins.length };
  };

  const inner = zoneStats(innerBins);
  const mid = zoneStats(midBins);
  const outer = zoneStats(outerBins);
  const all = zoneStats(bins);

  const suppressionInOut = outer.angC > 0 ? inner.angC / outer.angC : null;

  const innerNonRotFrac = (inner.nonRotP + mid.nonRotP + outer.nonRotP) > 0 ?
    inner.nonRotP / (inner.nonRotP + mid.nonRotP + outer.nonRotP) : 0;
  const outerNonRotFrac = (inner.nonRotP + mid.nonRotP + outer.nonRotP) > 0 ?
    outer.nonRotP / (inner.nonRotP + mid.nonRotP + outer.nonRotP) : 0;

  return {
    inner, mid, outer, all,
    suppressionInOut,
    innerNonRotFrac,
    outerNonRotFrac,
    nBins: bins.length,
    hasValidOuter: outer.nBins >= 2,
    radialProfile: bins,
  };
}


console.log('='.repeat(72));
console.log('PROGRAM 12 / DM-2V: VERIFICATION — CDM+SHAPE vs SIDM');
console.log('Fixes: fractional radial coverage, outer-coverage filter,');
console.log('       odd/even mode separation, LOO + bootstrap robustness');
console.log('='.repeat(72));

const dataDir = path.join(__dirname, '..', 'data', 'things-2d');
const thingsGals = [
  { name: 'NGC2841' }, { name: 'NGC5055' }, { name: 'NGC3521' },
  { name: 'NGC6946' }, { name: 'NGC7331' }, { name: 'NGC2403' },
  { name: 'NGC2903' }, { name: 'NGC3198' }, { name: 'NGC4826' },
  { name: 'NGC4736' }, { name: 'DDO154' }, { name: 'IC2574' },
];

const maxMode = 6;
const allAnalyzed = [];

for (const tg of thingsGals) {
  const fitsPath = path.join(dataDir, tg.name + '_MOM1.FITS');
  if (!fs.existsSync(fitsPath)) continue;
  if (fs.statSync(fitsPath).size < 1000) continue;
  const galData = gals.find(g => normalize(g.name) === normalize(tg.name));
  if (!galData) { console.log('  SKIP ' + tg.name + ' (not in SPARC)'); continue; }

  const fits = parseFITS(fitsPath);
  const res = analyzeGalaxy(fits, maxMode);
  if (!res) { console.log('  SKIP ' + tg.name + ' (too few bins)'); continue; }

  allAnalyzed.push({ name: tg.name, DQ: galData.DQ, Vflat: galData.Vflat, hR: galData.hR, ...res });
}

console.log('\n  Total galaxies analyzed: ' + allAnalyzed.length);
for (const g of allAnalyzed) {
  console.log('  ' + g.name + ' DQ=' + g.DQ.toFixed(2) + ' bins=' + g.nBins +
    ' inner=' + g.inner.angC.toFixed(4) + ' mid=' + g.mid.angC.toFixed(4) +
    ' outer=' + g.outer.angC.toFixed(4) + ' hasOuter=' + g.hasValidOuter +
    ' supp=' + (g.suppressionInOut !== null ? g.suppressionInOut.toFixed(3) : 'N/A'));
}

const clean = allAnalyzed.filter(g => g.hasValidOuter);
console.log('\n  Clean sample (valid outer coverage): ' + clean.length + ' / ' + allAnalyzed.length);
console.log('  Removed: ' + allAnalyzed.filter(g => !g.hasValidOuter).map(g => g.name).join(', '));
console.log('  Kept: ' + clean.map(g => g.name).join(', '));

if (clean.length < 4) {
  console.error('\n  ERROR: Too few galaxies with valid outer coverage. Cannot run clean analysis.');
  process.exit(1);
}

console.log('\n' + '='.repeat(72));
console.log('SECTION 1: CLEAN INNER vs OUTER SUPPRESSION (T1 redo)');
console.log('='.repeat(72));

const suppRatios = clean.map(g => g.suppressionInOut).filter(s => s !== null);
const meanSupp = suppRatios.reduce((a, b) => a + b, 0) / suppRatios.length;
const nSuppressed = suppRatios.filter(r => r < 0.8).length;
const rDQ_supp = pearsonR(clean.filter(g => g.suppressionInOut !== null).map(g => g.DQ),
  clean.filter(g => g.suppressionInOut !== null).map(g => g.suppressionInOut));

console.log('\n  Mean suppression ratio (inner/outer): ' + meanSupp.toFixed(3));
console.log('  N inner < outer: ' + nSuppressed + '/' + suppRatios.length);
console.log('  r(DQ, suppression): ' + rDQ_supp.toFixed(3));

for (const g of clean) {
  if (g.suppressionInOut !== null) {
    console.log('    ' + g.name + ': inner=' + g.inner.angC.toFixed(4) + ' outer=' + g.outer.angC.toFixed(4) + ' ratio=' + g.suppressionInOut.toFixed(3) + ' DQ=' + g.DQ.toFixed(2));
  }
}


console.log('\n' + '='.repeat(72));
console.log('SECTION 2: ODD vs EVEN MODE STRUCTURE');
console.log('Odd modes (m=3,5) = spiral/warp/asymmetry');
console.log('Even modes (m=2,4,6) = bar/triaxiality/ellipticity');
console.log('='.repeat(72));

for (const g of clean) {
  console.log('\n  ' + g.name + ' (DQ=' + g.DQ.toFixed(2) + '):');
  console.log('    Inner: odd=' + g.inner.oddFrac.toFixed(3) + ' even=' + g.inner.evenFrac.toFixed(3));
  console.log('    Mid:   odd=' + g.mid.oddFrac.toFixed(3) + ' even=' + g.mid.evenFrac.toFixed(3));
  console.log('    Outer: odd=' + g.outer.oddFrac.toFixed(3) + ' even=' + g.outer.evenFrac.toFixed(3));
}

const rDQ_innerOdd = pearsonR(clean.map(g => g.DQ), clean.map(g => g.inner.oddFrac));
const rDQ_outerOdd = pearsonR(clean.map(g => g.DQ), clean.map(g => g.outer.oddFrac));
const rDQ_innerEven = pearsonR(clean.map(g => g.DQ), clean.map(g => g.inner.evenFrac));
const rDQ_outerEven = pearsonR(clean.map(g => g.DQ), clean.map(g => g.outer.evenFrac));

console.log('\n  r(DQ, inner odd frac): ' + rDQ_innerOdd.toFixed(3));
console.log('  r(DQ, outer odd frac): ' + rDQ_outerOdd.toFixed(3));
console.log('  r(DQ, inner even frac): ' + rDQ_innerEven.toFixed(3));
console.log('  r(DQ, outer even frac): ' + rDQ_outerEven.toFixed(3));


console.log('\n' + '='.repeat(72));
console.log('SECTION 3: INNER-ONLY vs OUTER-ONLY COMPLEXITY (T2/T3 redo)');
console.log('='.repeat(72));

const rDQ_innerC = pearsonR(clean.map(g => g.DQ), clean.map(g => g.inner.angC));
const rDQ_midC = pearsonR(clean.map(g => g.DQ), clean.map(g => g.mid.angC));
const rDQ_outerC = pearsonR(clean.map(g => g.DQ), clean.map(g => g.outer.angC));
const rDQ_allC = pearsonR(clean.map(g => g.DQ), clean.map(g => g.all.angC));

console.log('\n  r(DQ, inner angular C): ' + rDQ_innerC.toFixed(3));
console.log('  r(DQ, mid angular C):   ' + rDQ_midC.toFixed(3));
console.log('  r(DQ, outer angular C): ' + rDQ_outerC.toFixed(3));
console.log('  r(DQ, global angular C): ' + rDQ_allC.toFixed(3));

const rDQ_innerNRF = pearsonR(clean.map(g => g.DQ), clean.map(g => g.innerNonRotFrac));
const rDQ_outerNRF = pearsonR(clean.map(g => g.DQ), clean.map(g => g.outerNonRotFrac));
console.log('\n  r(DQ, inner non-rot fraction of total): ' + rDQ_innerNRF.toFixed(3));
console.log('  r(DQ, outer non-rot fraction of total): ' + rDQ_outerNRF.toFixed(3));

const innerLed = Math.abs(rDQ_innerC) > Math.abs(rDQ_outerC);
console.log('\n  Inner-led (|r_inner| > |r_outer|): ' + innerLed);
console.log('  Difference: |r_inner| - |r_outer| = ' + (Math.abs(rDQ_innerC) - Math.abs(rDQ_outerC)).toFixed(3));


console.log('\n' + '='.repeat(72));
console.log('SECTION 4: LOO ROBUSTNESS');
console.log('='.repeat(72));

const looResults = [];
for (let drop = 0; drop < clean.length; drop++) {
  const sub = clean.filter((_, i) => i !== drop);
  const rInner = pearsonR(sub.map(g => g.DQ), sub.map(g => g.inner.angC));
  const rOuter = pearsonR(sub.map(g => g.DQ), sub.map(g => g.outer.angC));
  const rGlobal = pearsonR(sub.map(g => g.DQ), sub.map(g => g.all.angC));
  const rSupp = pearsonR(sub.filter(g => g.suppressionInOut !== null).map(g => g.DQ),
    sub.filter(g => g.suppressionInOut !== null).map(g => g.suppressionInOut));
  const innerLedLoo = Math.abs(rInner) > Math.abs(rOuter);
  looResults.push({ dropped: clean[drop].name, rInner, rOuter, rGlobal, rSupp, innerLed: innerLedLoo });
  console.log('  Drop ' + clean[drop].name + ': r_inner=' + rInner.toFixed(3) + ' r_outer=' + rOuter.toFixed(3) + ' r_global=' + rGlobal.toFixed(3) + ' inner-led=' + innerLedLoo);
}

const nInnerLedLoo = looResults.filter(l => l.innerLed).length;
const looInnerStable = nInnerLedLoo >= looResults.length * 0.6;
console.log('\n  Inner-led in LOO: ' + nInnerLedLoo + '/' + looResults.length + ' (' + (looInnerStable ? 'STABLE' : 'UNSTABLE') + ')');


console.log('\n' + '='.repeat(72));
console.log('SECTION 5: BOOTSTRAP CONFIDENCE');
console.log('='.repeat(72));

const nBoot = 5000;
let bootInnerLed = 0, bootSidmFav = 0;
const bootRInner = [], bootROuter = [];

for (let b = 0; b < nBoot; b++) {
  const sample = [];
  for (let i = 0; i < clean.length; i++) {
    sample.push(clean[Math.floor(rng() * clean.length)]);
  }
  const rI = pearsonR(sample.map(g => g.DQ), sample.map(g => g.inner.angC));
  const rO = pearsonR(sample.map(g => g.DQ), sample.map(g => g.outer.angC));
  bootRInner.push(rI);
  bootROuter.push(rO);
  if (Math.abs(rI) > Math.abs(rO)) bootInnerLed++;
  const suppVals = sample.filter(g => g.suppressionInOut !== null);
  if (suppVals.length >= 3) {
    const meanS = suppVals.reduce((s, g) => s + g.suppressionInOut, 0) / suppVals.length;
    if (meanS < 0.85 && Math.abs(rI) > Math.abs(rO)) bootSidmFav++;
  }
}

const bootInnerPct = (bootInnerLed / nBoot * 100).toFixed(1);
const bootSidmPct = (bootSidmFav / nBoot * 100).toFixed(1);
bootRInner.sort((a, b) => a - b);
bootROuter.sort((a, b) => a - b);
const ciInner = [bootRInner[Math.floor(nBoot * 0.025)], bootRInner[Math.floor(nBoot * 0.975)]];
const ciOuter = [bootROuter[Math.floor(nBoot * 0.025)], bootROuter[Math.floor(nBoot * 0.975)]];

console.log('\n  Bootstrap (N=' + nBoot + '):');
console.log('  Inner-led fraction: ' + bootInnerPct + '%');
console.log('  SIDM-favored fraction (inner-led + suppression < 0.85): ' + bootSidmPct + '%');
console.log('  r(DQ, inner C) 95% CI: [' + ciInner[0].toFixed(3) + ', ' + ciInner[1].toFixed(3) + ']');
console.log('  r(DQ, outer C) 95% CI: [' + ciOuter[0].toFixed(3) + ', ' + ciOuter[1].toFixed(3) + ']');


console.log('\n' + '='.repeat(72));
console.log('SECTION 6: FINAL CLEAN VERDICT');
console.log('='.repeat(72));

const T1clean = meanSupp < 0.7 ? 'SIDM' : (meanSupp > 0.9 ? 'CDM+shape' : 'INCONCLUSIVE');
const T2clean = Math.abs(rDQ_innerC) > Math.abs(rDQ_outerC) + 0.05 ? 'SIDM' :
  (Math.abs(rDQ_outerC) > Math.abs(rDQ_innerC) + 0.05 ? 'CDM+shape' : 'INCONCLUSIVE');
const T3clean = clean.filter(g => g.innerNonRotFrac > g.outerNonRotFrac).length > clean.length * 0.6 ? 'SIDM' :
  (clean.filter(g => g.outerNonRotFrac > g.innerNonRotFrac).length > clean.length * 0.6 ? 'CDM+shape' : 'INCONCLUSIVE');
const T4boot = parseFloat(bootInnerPct) > 65 ? 'SIDM' : (parseFloat(bootInnerPct) < 35 ? 'CDM+shape' : 'INCONCLUSIVE');

const cleanTests = [
  { test: 'T1-clean', name: 'Suppression (filtered)', favors: T1clean },
  { test: 'T2-clean', name: 'Inner vs outer correlation', favors: T2clean },
  { test: 'T3-clean', name: 'Non-rot power distribution', favors: T3clean },
  { test: 'T4-boot', name: 'Bootstrap robustness', favors: T4boot },
];

for (const t of cleanTests) {
  console.log('  ' + t.test + ' (' + t.name + '): ' + t.favors);
}

const cdmClean = cleanTests.filter(t => t.favors === 'CDM+shape').length;
const sidmClean = cleanTests.filter(t => t.favors === 'SIDM').length;
const incClean = cleanTests.filter(t => t.favors === 'INCONCLUSIVE').length;

console.log('\n  CDM+shape: ' + cdmClean + '/4');
console.log('  SIDM: ' + sidmClean + '/4');
console.log('  Inconclusive: ' + incClean + '/4');

let finalVerdict;
if (sidmClean >= 3) finalVerdict = 'SIDM SIGNAL SURVIVES CLEANING — genuine inner-dominated angular complexity';
else if (sidmClean === 2 && cdmClean <= 1 && looInnerStable) finalVerdict = 'SIDM LEANS — inner-led persists but marginal with N=' + clean.length;
else if (cdmClean >= 3) finalVerdict = 'CDM+shape CONFIRMED — DM-2 SIDM signal was coverage artifact';
else if (cdmClean === 2 && sidmClean <= 1) finalVerdict = 'CDM+shape LEANS — outer-led after cleaning';
else finalVerdict = 'NO CLEAR WINNER — DM-2 result was fragile, neither model decisively favored';

console.log('\n  FINAL VERDICT: ' + finalVerdict);

const fragile = (sidmClean <= 1 && cdmClean <= 1) || !looInnerStable;
console.log('  DM-2 original SIDM advantage was: ' + (fragile ? 'FRAGILE' : 'ROBUST'));


const output = {
  program: 12,
  phase: 'DM-2V',
  title: 'DM-2 Verification: CDM+Shape vs SIDM after cleaning',
  timestamp: new Date().toISOString(),
  totalGalaxies: allAnalyzed.length,
  cleanGalaxies: clean.length,
  removedForCoverage: allAnalyzed.filter(g => !g.hasValidOuter).map(g => g.name),
  keptGalaxies: clean.map(g => g.name),

  cleanSuppressionTest: {
    meanSuppression: meanSupp,
    nSuppressed,
    rDQ_suppression: rDQ_supp,
    verdict: T1clean,
  },

  radialCorrelations: {
    rDQ_innerC, rDQ_midC, rDQ_outerC, rDQ_allC,
    innerLed,
    difference: Math.abs(rDQ_innerC) - Math.abs(rDQ_outerC),
    verdict: T2clean,
  },

  oddEvenStructure: {
    rDQ_innerOdd, rDQ_outerOdd,
    rDQ_innerEven, rDQ_outerEven,
  },

  nonRotDistribution: {
    perGalaxy: clean.map(g => ({ name: g.name, DQ: g.DQ, innerFrac: g.innerNonRotFrac, outerFrac: g.outerNonRotFrac })),
    verdict: T3clean,
  },

  loo: {
    nInnerLed: nInnerLedLoo,
    total: looResults.length,
    stable: looInnerStable,
    perDrop: looResults,
  },

  bootstrap: {
    nBoot,
    innerLedPct: parseFloat(bootInnerPct),
    sidmFavoredPct: parseFloat(bootSidmPct),
    ciRInner: ciInner,
    ciROuter: ciOuter,
    verdict: T4boot,
  },

  cleanVerdict: {
    T1: T1clean, T2: T2clean, T3: T3clean, T4: T4boot,
    cdmWins: cdmClean, sidmWins: sidmClean, inconclusive: incClean,
    overall: finalVerdict,
    dm2WasFragile: fragile,
  },

  perGalaxy: clean.map(g => ({
    name: g.name, DQ: g.DQ, Vflat: g.Vflat,
    innerAngC: g.inner.angC, midAngC: g.mid.angC, outerAngC: g.outer.angC,
    globalAngC: g.all.angC,
    suppressionRatio: g.suppressionInOut,
    innerOddFrac: g.inner.oddFrac, outerOddFrac: g.outer.oddFrac,
    innerEvenFrac: g.inner.evenFrac, outerEvenFrac: g.outer.evenFrac,
    innerNonRotFrac: g.innerNonRotFrac, outerNonRotFrac: g.outerNonRotFrac,
  })),
};

const outDir = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results');
fs.mkdirSync(outDir, { recursive: true });
const outPath = path.join(outDir, 'program12-DM2V-verification.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nResult saved: ' + outPath);
console.log('\n=== DM-2V COMPLETE ===');
