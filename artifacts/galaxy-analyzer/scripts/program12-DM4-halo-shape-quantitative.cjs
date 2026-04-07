#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

function seededRNG(seed) { let s = seed | 0; return function() { s = (s * 1103515245 + 12345) & 0x7fffffff; return s / 0x7fffffff; }; }
const rng = seededRNG(20260409);

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
    const angC = totalP > 0 ? nonRotP / totalP : 0;

    bins.push({ rNorm, mAmp, mPhase, totalP, nonRotP, angC });
  }
  if (bins.length < 4) return null;

  const m2Amps = bins.map(b => b.mAmp[2]);
  const m2Mean = m2Amps.reduce((a, b) => a + b, 0) / m2Amps.length;
  const m2Std = Math.sqrt(m2Amps.reduce((s, a) => s + (a - m2Mean) ** 2, 0) / m2Amps.length);

  const ellipticityProxy = m2Mean;
  const ellipticityCV = m2Mean > 0 ? m2Std / m2Mean : 0;

  const nonAxiAmps = bins.map(b => {
    let sum = 0;
    for (let m = 2; m <= maxMode; m++) sum += b.mAmp[m];
    return sum;
  });
  const shapeAmplitude = nonAxiAmps.reduce((a, b) => a + b, 0) / nonAxiAmps.length;
  const shapeAmpStd = Math.sqrt(nonAxiAmps.reduce((s, a) => s + (a - shapeAmplitude) ** 2, 0) / nonAxiAmps.length);
  const shapeAmpCV = shapeAmplitude > 0 ? shapeAmpStd / shapeAmplitude : 0;

  const phaseCoherence = [];
  for (let m = 2; m <= maxMode; m++) {
    const phases = bins.map(b => b.mPhase[m]);
    let sumCos = 0, sumSin = 0;
    for (const p of phases) { sumCos += Math.cos(p); sumSin += Math.sin(p); }
    const R = Math.sqrt(sumCos ** 2 + sumSin ** 2) / phases.length;
    phaseCoherence.push({ m, R, meanAngle: Math.atan2(sumSin, sumCos) });
  }
  const meanPhaseCoherence = phaseCoherence.reduce((s, p) => s + p.R, 0) / phaseCoherence.length;

  const paGradients = [];
  for (let m = 2; m <= Math.min(4, maxMode); m++) {
    const phases = bins.map(b => b.mPhase[m]);
    const rNorms = bins.map(b => b.rNorm);
    if (phases.length < 3) continue;
    const unwrapped = [phases[0]];
    for (let i = 1; i < phases.length; i++) {
      let dp = phases[i] - unwrapped[i - 1];
      while (dp > Math.PI) dp -= 2 * Math.PI;
      while (dp < -Math.PI) dp += 2 * Math.PI;
      unwrapped.push(unwrapped[i - 1] + dp);
    }
    const totalTwist = Math.abs(unwrapped[unwrapped.length - 1] - unwrapped[0]);
    const twistPerR = totalTwist / (rNorms[rNorms.length - 1] - rNorms[0] + 0.01);
    paGradients.push({ m, totalTwist: totalTwist * 180 / Math.PI, twistPerR: twistPerR * 180 / Math.PI });
  }
  const meanTwist = paGradients.length > 0 ? paGradients.reduce((s, p) => s + p.totalTwist, 0) / paGradients.length : 0;

  const isotropyIndex = (() => {
    const globalAmps = Array(maxMode + 1).fill(0);
    for (const b of bins) { for (let m = 0; m <= maxMode; m++) globalAmps[m] += b.mAmp[m]; }
    const nonRot = globalAmps.slice(2).reduce((a, b) => a + b, 0);
    if (nonRot <= 0) return 0;
    const fracs = [];
    for (let m = 2; m <= maxMode; m++) fracs.push(globalAmps[m] / nonRot);
    let H = 0;
    for (const f of fracs) if (f > 1e-10) H -= f * Math.log2(f);
    return H / Math.log2(fracs.length);
  })();

  const oddPowerFrac = (() => {
    let odd = 0, total = 0;
    for (const b of bins) { odd += b.mAmp[3] + b.mAmp[5]; for (let m = 2; m <= maxMode; m++) total += b.mAmp[m]; }
    return total > 0 ? odd / total : 0;
  })();

  const shapeCoherence = meanPhaseCoherence * (1 - shapeAmpCV);

  const globalAngC = bins.reduce((s, b) => s + b.angC, 0) / bins.length;

  return {
    ellipticityProxy, ellipticityCV,
    shapeAmplitude, shapeAmpCV,
    phaseCoherence, meanPhaseCoherence,
    paGradients, meanTwist,
    isotropyIndex, oddPowerFrac,
    shapeCoherence,
    globalAngC,
    nBins: bins.length,
    radialProfile: bins,
  };
}


console.log('='.repeat(72));
console.log('PROGRAM 12 / DM-4: QUANTITATIVE HALO SHAPE TEST');
console.log('Question: What shape parameter(s) best capture H?');
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
  analyzed.push({ name: tg.name, DQ: galData.DQ, Vflat: galData.Vflat, hR: galData.hR, ...res });
}

console.log('\n  Galaxies analyzed: ' + analyzed.length);


console.log('\n' + '='.repeat(72));
console.log('TEST T1: NON-AXISYMMETRY STRENGTH');
console.log('Does H increase with the magnitude of non-axisymmetric structure?');
console.log('='.repeat(72));

for (const g of analyzed) {
  console.log('\n  ' + g.name + ' (DQ=' + g.DQ.toFixed(2) + '):');
  console.log('    Shape amplitude (total non-axi): ' + g.shapeAmplitude.toFixed(4));
  console.log('    Ellipticity proxy (m=2 amp): ' + g.ellipticityProxy.toFixed(4));
  console.log('    Angular complexity: ' + g.globalAngC.toFixed(4));
}

const rDQ_shapeAmp = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.shapeAmplitude));
const rDQ_ellip = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.ellipticityProxy));
const rDQ_angC = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.globalAngC));
const rDQ_hR = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.hR));

const pDQ_shapeAmp = permutationP(analyzed.map(g => g.DQ), analyzed.map(g => g.shapeAmplitude), 10000, rDQ_shapeAmp);
const pDQ_angC = permutationP(analyzed.map(g => g.DQ), analyzed.map(g => g.globalAngC), 10000, rDQ_angC);

console.log('\n  Correlations with DQ (H):');
console.log('    r(DQ, shape amplitude): ' + rDQ_shapeAmp.toFixed(3) + ' (p=' + pDQ_shapeAmp.toFixed(4) + ')');
console.log('    r(DQ, ellipticity proxy): ' + rDQ_ellip.toFixed(3));
console.log('    r(DQ, angular complexity): ' + rDQ_angC.toFixed(3) + ' (p=' + pDQ_angC.toFixed(4) + ')');
console.log('    r(DQ, haloResponse): ' + rDQ_hR.toFixed(3) + ' [reference]');

const T1_bestParam = Math.abs(rDQ_shapeAmp) > Math.abs(rDQ_angC) ? 'shapeAmplitude' : 'angularComplexity';
const T1_bestR = Math.max(Math.abs(rDQ_shapeAmp), Math.abs(rDQ_angC));
const beatsHR = T1_bestR > Math.abs(rDQ_hR);

console.log('\n  Best shape proxy: ' + T1_bestParam + ' (r=' + T1_bestR.toFixed(3) + ')');
console.log('  Beats haloResponse (r=' + Math.abs(rDQ_hR).toFixed(3) + '): ' + beatsHR);


console.log('\n' + '='.repeat(72));
console.log('TEST T2: DIRECTIONAL STABILITY (SHAPE COHERENCE)');
console.log('Is the shape orientation stable with radius, or twisted/random?');
console.log('='.repeat(72));

for (const g of analyzed) {
  console.log('\n  ' + g.name + ' (DQ=' + g.DQ.toFixed(2) + '):');
  console.log('    Mean phase coherence (R): ' + g.meanPhaseCoherence.toFixed(3));
  console.log('    Shape amplitude CV: ' + g.shapeAmpCV.toFixed(3));
  console.log('    Combined shape coherence: ' + g.shapeCoherence.toFixed(3));
  console.log('    Mean PA twist (deg): ' + g.meanTwist.toFixed(1));
  for (const pa of g.paGradients) {
    console.log('      m=' + pa.m + ': total twist=' + pa.totalTwist.toFixed(1) + ' deg');
  }
}

const rDQ_phaseCoherence = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.meanPhaseCoherence));
const rDQ_shapeCoherence = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.shapeCoherence));
const rDQ_twist = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.meanTwist));
const rDQ_ampCV = pearsonR(analyzed.map(g => g.DQ), analyzed.map(g => g.shapeAmpCV));

const pDQ_shapeCoherence = permutationP(analyzed.map(g => g.DQ), analyzed.map(g => g.shapeCoherence), 10000, rDQ_shapeCoherence);

console.log('\n  Correlations with DQ:');
console.log('    r(DQ, phase coherence): ' + rDQ_phaseCoherence.toFixed(3));
console.log('    r(DQ, shape coherence): ' + rDQ_shapeCoherence.toFixed(3) + ' (p=' + pDQ_shapeCoherence.toFixed(4) + ')');
console.log('    r(DQ, PA twist): ' + rDQ_twist.toFixed(3));
console.log('    r(DQ, amplitude CV): ' + rDQ_ampCV.toFixed(3));

const coherenceMatters = Math.abs(rDQ_shapeCoherence) > Math.abs(rDQ_shapeAmp) * 0.5;
console.log('\n  Coherence matters (|r_coherence| > 0.5*|r_amplitude|): ' + coherenceMatters);
console.log('  Interpretation: H driven by ' + (coherenceMatters ? 'ORGANIZED shape (amplitude + coherence)' : 'shape AMPLITUDE alone'));


console.log('\n' + '='.repeat(72));
console.log('TEST T3: WHAT COUPLES TO H — AMPLITUDE, COHERENCE, OR BOTH?');
console.log('='.repeat(72));

const candidates = [
  { name: 'shapeAmplitude', vals: analyzed.map(g => g.shapeAmplitude) },
  { name: 'ellipticityProxy', vals: analyzed.map(g => g.ellipticityProxy) },
  { name: 'angularComplexity', vals: analyzed.map(g => g.globalAngC) },
  { name: 'shapeCoherence', vals: analyzed.map(g => g.shapeCoherence) },
  { name: 'phaseCoherence', vals: analyzed.map(g => g.meanPhaseCoherence) },
  { name: 'isotropyIndex', vals: analyzed.map(g => g.isotropyIndex) },
  { name: 'oddPowerFraction', vals: analyzed.map(g => g.oddPowerFrac) },
  { name: 'meanTwist', vals: analyzed.map(g => g.meanTwist) },
  { name: 'haloResponse', vals: analyzed.map(g => g.hR) },
];

const DQs = analyzed.map(g => g.DQ);
const ranked = candidates.map(c => {
  const r = pearsonR(DQs, c.vals);
  return { name: c.name, r, absR: Math.abs(r) };
}).sort((a, b) => b.absR - a.absR);

console.log('\n  Parameter ranking by |r(DQ, param)|:');
for (const c of ranked) {
  const marker = c.name === 'haloResponse' ? ' <-- reference' : '';
  console.log('    ' + c.name + ': r=' + c.r.toFixed(3) + ' |r|=' + c.absR.toFixed(3) + marker);
}

const top1 = ranked[0];
const top2 = ranked[1];
console.log('\n  Best single parameter: ' + top1.name + ' (r=' + top1.r.toFixed(3) + ')');
console.log('  Second best: ' + top2.name + ' (r=' + top2.r.toFixed(3) + ')');

const top1Vals = candidates.find(c => c.name === top1.name).vals;
const top2Vals = candidates.find(c => c.name === top2.name).vals;

const r12 = pearsonR(top1Vals, top2Vals);
console.log('  Correlation between top-1 and top-2: r=' + r12.toFixed(3));
const independent = Math.abs(r12) < 0.8;
console.log('  Independent (|r|<0.8): ' + independent);

if (independent) {
  const combined = analyzed.map((_, i) => (zscore(top1Vals)[i] + zscore(top2Vals)[i]) / 2);
  const rCombined = pearsonR(DQs, combined);
  const pCombined = permutationP(DQs, combined, 10000, rCombined);
  console.log('\n  Combined parameter (z-avg of top-1 + top-2):');
  console.log('    r(DQ, combined): ' + rCombined.toFixed(3) + ' (p=' + pCombined.toFixed(4) + ')');
  console.log('    Improvement over best single: ' + (Math.abs(rCombined) - top1.absR > 0.02 ? 'YES' : 'NO'));
}


console.log('\n' + '='.repeat(72));
console.log('TEST T4: LOO + BOOTSTRAP STABILITY');
console.log('='.repeat(72));

const bestParamVals = candidates.find(c => c.name === top1.name).vals;
const bestParamName = top1.name;

console.log('\n  LOO for best parameter (' + bestParamName + '):');
const looR = [];
for (let drop = 0; drop < analyzed.length; drop++) {
  const subDQ = DQs.filter((_, i) => i !== drop);
  const subVals = bestParamVals.filter((_, i) => i !== drop);
  const r = pearsonR(subDQ, subVals);
  looR.push(r);
  console.log('    Drop ' + analyzed[drop].name + ': r=' + r.toFixed(3));
}
const looMin = Math.min(...looR);
const looMax = Math.max(...looR);
const looMean = looR.reduce((a, b) => a + b, 0) / looR.length;
const looStable = looR.every(r => Math.sign(r) === Math.sign(top1.r));
console.log('    LOO range: [' + looMin.toFixed(3) + ', ' + looMax.toFixed(3) + '] mean=' + looMean.toFixed(3));
console.log('    Sign-stable: ' + looStable);

console.log('\n  Bootstrap (N=5000) for best parameter:');
const nBoot = 5000;
const bootR = [];
for (let b = 0; b < nBoot; b++) {
  const sample = [];
  for (let i = 0; i < analyzed.length; i++) {
    const idx = Math.floor(rng() * analyzed.length);
    sample.push({ dq: DQs[idx], val: bestParamVals[idx] });
  }
  const r = pearsonR(sample.map(s => s.dq), sample.map(s => s.val));
  bootR.push(r);
}
bootR.sort((a, b) => a - b);
const ci95 = [bootR[Math.floor(nBoot * 0.025)], bootR[Math.floor(nBoot * 0.975)]];
const bootMean = bootR.reduce((a, b) => a + b, 0) / nBoot;
const bootPctPositive = bootR.filter(r => r > 0).length / nBoot * 100;

console.log('    Bootstrap mean r: ' + bootMean.toFixed(3));
console.log('    95% CI: [' + ci95[0].toFixed(3) + ', ' + ci95[1].toFixed(3) + ']');
console.log('    % positive: ' + bootPctPositive.toFixed(1) + '%');

const T4_stable = looStable && !ci95.includes(0) && (ci95[0] > 0 || ci95[1] < 0);
console.log('    Overall stability: ' + (T4_stable ? 'STABLE' : 'MARGINAL (wide CI expected with N=7)'));


console.log('\n' + '='.repeat(72));
console.log('TEST T5: MATCHED PAIRS');
console.log('='.repeat(72));

const sortedByDQ = analyzed.slice().sort((a, b) => b.DQ - a.DQ);
const half = Math.floor(sortedByDQ.length / 2);
const highH = sortedByDQ.slice(0, half);
const lowH = sortedByDQ.slice(-half);

console.log('\n  High-H group: ' + highH.map(g => g.name + '(' + g.DQ.toFixed(2) + ')').join(', '));
console.log('  Low-H group: ' + lowH.map(g => g.name + '(' + g.DQ.toFixed(2) + ')').join(', '));

const metrics = ['shapeAmplitude', 'globalAngC', 'shapeCoherence', 'meanPhaseCoherence', 'ellipticityProxy', 'oddPowerFrac'];
for (const metric of metrics) {
  const highMean = highH.reduce((s, g) => s + g[metric], 0) / highH.length;
  const lowMean = lowH.reduce((s, g) => s + g[metric], 0) / lowH.length;
  const diff = highMean - lowMean;
  const pctDiff = lowMean !== 0 ? (diff / Math.abs(lowMean) * 100).toFixed(1) : 'Inf';
  console.log('  ' + metric + ': high=' + highMean.toFixed(4) + ' low=' + lowMean.toFixed(4) + ' diff=' + diff.toFixed(4) + ' (' + pctDiff + '%)');
}


console.log('\n' + '='.repeat(72));
console.log('FINAL SUMMARY: BEST SHAPE PARAMETERS FOR H');
console.log('='.repeat(72));

console.log('\n  Parameter ranking (top 3):');
for (let i = 0; i < Math.min(3, ranked.length); i++) {
  console.log('    ' + (i + 1) + '. ' + ranked[i].name + ': r(DQ)=' + ranked[i].r.toFixed(3));
}

const bestBeatsHR = ranked[0].absR > Math.abs(rDQ_hR);
const hrRank = ranked.findIndex(r => r.name === 'haloResponse') + 1;

console.log('\n  haloResponse rank: ' + hrRank + '/' + ranked.length + ' (r=' + rDQ_hR.toFixed(3) + ')');
console.log('  Best shape param beats haloResponse: ' + bestBeatsHR);

let overallVerdict;
if (ranked[0].absR > 0.6 && looStable) {
  overallVerdict = 'SUCCESS — shape parameter "' + ranked[0].name + '" captures H with r=' + ranked[0].r.toFixed(3) + ', LOO-stable. Halo non-axisymmetry is the quantitative carrier of H.';
} else if (ranked[0].absR > 0.5) {
  overallVerdict = 'PARTIAL SUCCESS — shape parameter correlates with H but not dominant. H may require multi-parameter description.';
} else {
  overallVerdict = 'H deeper than shape alone — shape parameters do not fully capture H. Need to investigate dark sector coupling or other mechanisms.';
}

console.log('\n  OVERALL VERDICT: ' + overallVerdict);

const dm2vResult = (() => { try { return JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results', 'program12-DM2V-verification.json'), 'utf8')); } catch { return null; } })();
const dm2v_rOuter = dm2vResult ? dm2vResult.radialCorrelations.rDQ_outerC : null;

console.log('\n  Cross-check with DM-2V:');
if (dm2v_rOuter !== null) {
  console.log('    DM-2V r(DQ, outer angC): ' + dm2v_rOuter.toFixed(3));
  console.log('    DM-4 r(DQ, angC): ' + rDQ_angC.toFixed(3));
  console.log('    Consistent: ' + (Math.abs(dm2v_rOuter - rDQ_angC) < 0.15 ? 'YES' : 'CHECK'));
}


const output = {
  program: 12,
  phase: 'DM-4',
  title: 'Quantitative Halo Shape Test — What shape parameter best captures H?',
  timestamp: new Date().toISOString(),
  N: analyzed.length,
  galaxies: analyzed.map(g => g.name),

  T1: {
    name: 'Non-axisymmetry strength',
    rDQ_shapeAmp, pDQ_shapeAmp,
    rDQ_ellipticity: rDQ_ellip,
    rDQ_angularComplexity: rDQ_angC, pDQ_angC,
    rDQ_haloResponse: rDQ_hR,
    bestProxy: T1_bestParam,
    bestR: T1_bestR,
    beatsHaloResponse: beatsHR,
  },

  T2: {
    name: 'Directional stability',
    rDQ_phaseCoherence, rDQ_shapeCoherence, pDQ_shapeCoherence,
    rDQ_twist, rDQ_ampCV,
    coherenceMatters,
  },

  T3: {
    name: 'Parameter ranking',
    ranking: ranked,
    bestSingle: top1,
    secondBest: top2,
    correlation_top1_top2: r12,
    independent,
  },

  T4: {
    name: 'LOO + Bootstrap stability',
    bestParameter: bestParamName,
    loo: { min: looMin, max: looMax, mean: looMean, signStable: looStable },
    bootstrap: { mean: bootMean, ci95, pctPositive: bootPctPositive, nBoot },
    stable: T4_stable,
  },

  T5: {
    name: 'Matched pairs',
    highH: highH.map(g => g.name),
    lowH: lowH.map(g => g.name),
  },

  overallVerdict,

  perGalaxy: analyzed.map(g => ({
    name: g.name, DQ: g.DQ, Vflat: g.Vflat, hR: g.hR,
    shapeAmplitude: g.shapeAmplitude,
    ellipticityProxy: g.ellipticityProxy,
    angularComplexity: g.globalAngC,
    shapeCoherence: g.shapeCoherence,
    phaseCoherence: g.meanPhaseCoherence,
    isotropyIndex: g.isotropyIndex,
    oddPowerFrac: g.oddPowerFrac,
    meanTwist: g.meanTwist,
  })),
};

const outDir = path.join(__dirname, '..', 'public', 'replication', 'dark-matter-program', 'results');
fs.mkdirSync(outDir, { recursive: true });
const outPath = path.join(outDir, 'program12-DM4-halo-shape-quantitative.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\nResult saved: ' + outPath);
console.log('\n=== DM-4 COMPLETE ===');
