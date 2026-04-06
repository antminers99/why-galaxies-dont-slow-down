const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

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
  gals.push({ name: g.name, logVflat: Math.log10(sp.Vflat), Vflat: sp.Vflat, logMbar: Math.log10(Math.max(Mbar, 1)), logL36: Math.log10(Math.max(sp.L36, 0.001)), logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)), morphT: sp.T, logMHI: g.logMHI, logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)), logA0: g.logA0, hR, dist: sp.D, inc: sp.Inc || 0 });
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

function computeM2Power(fits) {
  const nx = fits.dims[0], ny = fits.dims[1];
  const crpix1 = fits.header['CRPIX1'] || nx/2, crpix2 = fits.header['CRPIX2'] || ny/2;
  const cdelt1 = fits.header['CDELT1'] || 1;
  const pixScale = Math.abs(cdelt1) * 3600;

  const vel = [], coords = [];
  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const v = fits.data[j * nx + i]; if (isNaN(v)) continue;
      vel.push(v);
      coords.push({ x: (i + 1 - crpix1) * cdelt1 * 3600, y: (j + 1 - crpix2) * (fits.header['CDELT2'] || 1) * 3600 });
    }
  }
  if (vel.length < 100) return { m2Power: NaN, m1Power: NaN, m2m0: NaN, paTwist: NaN, coherenceOuter: NaN, nValid: vel.length };

  const sorted = vel.slice().sort((a, b) => a - b);
  const vsys = sorted[Math.floor(sorted.length / 2)];
  const vrel = vel.map(v => v - vsys);

  const maxR = Math.max(nx, ny) * Math.abs(cdelt1) * 3600 / 2;
  const radialBins = 15, azBins = 8;
  const binW = maxR / radialBins;
  let m2Total = 0, m1Total = 0, m0Total = 0, validBins = 0;

  for (let rb = 0; rb < radialBins; rb++) {
    const rMin = rb * binW, rMax = (rb + 1) * binW;
    const azVals = Array.from({ length: azBins }, () => []);
    for (let k = 0; k < coords.length; k++) {
      const r = Math.sqrt(coords[k].x ** 2 + coords[k].y ** 2);
      if (r < rMin || r >= rMax) continue;
      const theta = Math.atan2(coords[k].y, coords[k].x);
      const azIdx = Math.floor(((theta + Math.PI) / (2 * Math.PI)) * azBins) % azBins;
      azVals[azIdx].push(vrel[k]);
    }
    const azMeans = azVals.map(a => a.length > 0 ? a.reduce((s, v) => s + v, 0) / a.length : NaN);
    const valid = azMeans.filter(v => !isNaN(v));
    if (valid.length < 4) continue;
    const mean = valid.reduce((a, b) => a + b, 0) / valid.length;
    let c1 = 0, s1 = 0, c2 = 0, s2 = 0, cnt = 0;
    for (let a = 0; a < azBins; a++) {
      if (isNaN(azMeans[a])) continue;
      const ang = (a + 0.5) * 2 * Math.PI / azBins;
      const dv = azMeans[a] - mean;
      c1 += dv * Math.cos(ang); s1 += dv * Math.sin(ang);
      c2 += dv * Math.cos(2 * ang); s2 += dv * Math.sin(2 * ang);
      cnt++;
    }
    m2Total += Math.sqrt(c2 ** 2 + s2 ** 2) / cnt;
    m1Total += Math.sqrt(c1 ** 2 + s1 ** 2) / cnt;
    m0Total += Math.abs(mean);
    validBins++;
  }
  const m2Power = validBins > 0 ? m2Total / validBins : NaN;
  const m1Power = validBins > 0 ? m1Total / validBins : NaN;
  const m0Power = validBins > 0 ? m0Total / validBins : NaN;
  const m2m0 = m0Power > 0 ? m2Power / m0Power : NaN;

  const innerR = maxR * 0.3, outerR = maxR * 0.7;
  let inPA_x = 0, inPA_y = 0, outPA_x = 0, outPA_y = 0;
  for (let k = 0; k < coords.length; k++) {
    const r = Math.sqrt(coords[k].x ** 2 + coords[k].y ** 2);
    if (r < innerR && r > 5) { inPA_x += vrel[k] * coords[k].x; inPA_y += vrel[k] * coords[k].y; }
    else if (r > outerR && r < maxR * 0.9) { outPA_x += vrel[k] * coords[k].x; outPA_y += vrel[k] * coords[k].y; }
  }
  const inPA = Math.atan2(inPA_y, inPA_x);
  const outPA = Math.atan2(outPA_y, outPA_x);
  let twist = Math.abs(outPA - inPA) * 180 / Math.PI;
  if (twist > 180) twist = 360 - twist;

  let cohOuter = 0, nOut = 0;
  const kinPA = (Math.atan2(inPA_y + outPA_y, inPA_x + outPA_x));
  for (let k = 0; k < coords.length; k++) {
    const r = Math.sqrt(coords[k].x ** 2 + coords[k].y ** 2);
    if (r > outerR) {
      const es = Math.sign(coords[k].x * Math.cos(kinPA) + coords[k].y * Math.sin(kinPA));
      cohOuter += (Math.sign(vrel[k]) === es || vrel[k] === 0) ? 1 : 0;
      nOut++;
    }
  }
  cohOuter = nOut > 0 ? cohOuter / nOut : NaN;

  return { m2Power, m1Power, m0Power, m2m0, paTwist: twist, coherenceOuter: cohOuter, nValid: vel.length };
}


console.log('='.repeat(72));
console.log('PROGRAM 9 — PHASE 903: DECISIVE MATCHED IFU TEST');
console.log('Expanded sample + permutation validation');
console.log('='.repeat(72));

const dataDir = path.join(__dirname, '..', 'data', 'things-2d');
const thingsGals = [
  { name: 'NGC2841', matchName: 'NGC2841' },
  { name: 'NGC5055', matchName: 'NGC5055' },
  { name: 'NGC3521', matchName: 'NGC3521' },
  { name: 'NGC6946', matchName: 'NGC6946' },
  { name: 'NGC7331', matchName: 'NGC7331' },
  { name: 'NGC2403', matchName: 'NGC2403' },
  { name: 'NGC2903', matchName: 'NGC2903' },
  { name: 'NGC3198', matchName: 'NGC3198' },
  { name: 'NGC4826', matchName: 'NGC4826' },
  { name: 'NGC4736', matchName: 'NGC4736' },
];

const results = [];
for (const tg of thingsGals) {
  const fitsPath = path.join(dataDir, tg.name + '_MOM1.FITS');
  if (!fs.existsSync(fitsPath)) { console.log('SKIP ' + tg.name + ' (no file)'); continue; }
  const fStat = fs.statSync(fitsPath);
  if (fStat.size < 1000) { console.log('SKIP ' + tg.name + ' (file too small: ' + fStat.size + ' bytes)'); continue; }

  const galData = gals.find(g => normalize(g.name) === normalize(tg.matchName));
  if (!galData) { console.log('SKIP ' + tg.name + ' (not in N=55)'); continue; }

  const fits = parseFITS(fitsPath);
  const metrics = computeM2Power(fits);

  console.log(tg.name.padEnd(12) + 'DQ=' + galData.DQ.toFixed(2).padEnd(8) + 'Vflat=' + ('' + Math.round(galData.Vflat)).padEnd(6) + 'm2=' + metrics.m2Power.toFixed(0).padEnd(10) + 'm1=' + metrics.m1Power.toFixed(0).padEnd(10) + 'twist=' + metrics.paTwist.toFixed(1).padEnd(8) + 'cohOut=' + metrics.coherenceOuter.toFixed(3) + '  nPix=' + metrics.nValid);

  results.push({ name: tg.name, DQ: galData.DQ, Vflat: galData.Vflat, dist: galData.dist, inc: galData.inc, logMbar: galData.logMbar, m2Power: metrics.m2Power, m1Power: metrics.m1Power, m0Power: metrics.m0Power, m2m0: metrics.m2m0, paTwist: metrics.paTwist, coherenceOuter: metrics.coherenceOuter, nValid: metrics.nValid });
}

console.log('\n  Total galaxies analyzed: ' + results.length);


console.log('\n\n' + '#'.repeat(72));
console.log('903.1 — DQ-m2 CORRELATION (EXPANDED SAMPLE)');
console.log('#'.repeat(72));

const validR = results.filter(r => !isNaN(r.m2Power) && !isNaN(r.DQ));
const dqArr = validR.map(r => r.DQ);
const m2Arr = validR.map(r => r.m2Power);
const logM2 = validR.map(r => Math.log10(Math.max(r.m2Power, 1)));
const m1Arr = validR.map(r => r.m1Power);
const twistArr = validR.map(r => r.paTwist);
const cohArr = validR.map(r => r.coherenceOuter);

const rDQ_m2 = pearsonR(dqArr, m2Arr);
const rDQ_logM2 = pearsonR(dqArr, logM2);
const rDQ_m1 = pearsonR(dqArr, m1Arr);
const rDQ_twist = pearsonR(dqArr, twistArr);
const rDQ_coh = pearsonR(dqArr, cohArr);

console.log('\n  N = ' + validR.length);
console.log('  r(DQ, m2_power) = ' + rDQ_m2.toFixed(4));
console.log('  r(DQ, log_m2)   = ' + rDQ_logM2.toFixed(4));
console.log('  r(DQ, m1_power) = ' + rDQ_m1.toFixed(4));
console.log('  r(DQ, PA_twist) = ' + rDQ_twist.toFixed(4));
console.log('  r(DQ, coh_out)  = ' + rDQ_coh.toFixed(4));

const rVf_m2 = pearsonR(validR.map(r => Math.log10(r.Vflat)), m2Arr);
const rMbar_m2 = pearsonR(validR.map(r => r.logMbar), m2Arr);
console.log('\n  Confounders:');
console.log('  r(logVflat, m2) = ' + rVf_m2.toFixed(4));
console.log('  r(logMbar, m2)  = ' + rMbar_m2.toFixed(4));


console.log('\n\n' + '#'.repeat(72));
console.log('903.2 — PERMUTATION TEST');
console.log('#'.repeat(72));

const nPerm = 10000;
let countGE_m2 = 0, countGE_logM2 = 0;
const observedR_m2 = rDQ_m2;
const observedR_logM2 = rDQ_logM2;

for (let p = 0; p < nPerm; p++) {
  const shuffledDQ = dqArr.slice();
  for (let i = shuffledDQ.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    const tmp = shuffledDQ[i]; shuffledDQ[i] = shuffledDQ[j]; shuffledDQ[j] = tmp;
  }
  const rPerm_m2 = pearsonR(shuffledDQ, m2Arr);
  const rPerm_logM2 = pearsonR(shuffledDQ, logM2);
  if (rPerm_m2 >= observedR_m2) countGE_m2++;
  if (rPerm_logM2 >= observedR_logM2) countGE_logM2++;
}

const pVal_m2 = countGE_m2 / nPerm;
const pVal_logM2 = countGE_logM2 / nPerm;
console.log('\n  Permutation test (N=' + nPerm + '):');
console.log('  r(DQ, m2) = ' + observedR_m2.toFixed(4) + ', p = ' + pVal_m2.toFixed(4) + (pVal_m2 < 0.05 ? ' * SIGNIFICANT' : ''));
console.log('  r(DQ, log_m2) = ' + observedR_logM2.toFixed(4) + ', p = ' + pVal_logM2.toFixed(4) + (pVal_logM2 < 0.05 ? ' * SIGNIFICANT' : ''));


console.log('\n\n' + '#'.repeat(72));
console.log('903.3 — PARTIAL CORRELATION (controlling Vflat)');
console.log('#'.repeat(72));

const logVf = validR.map(r => Math.log10(r.Vflat));
const dqResid = ols(logVf.map(v => [v]), dqArr).residuals;
const m2Resid = ols(logVf.map(v => [v]), m2Arr).residuals;
const logM2Resid = ols(logVf.map(v => [v]), logM2).residuals;
const rPartial_m2 = pearsonR(dqResid, m2Resid);
const rPartial_logM2 = pearsonR(dqResid, logM2Resid);

console.log('\n  Partial correlation (controlling logVflat):');
console.log('  r_partial(DQ, m2 | logVflat) = ' + rPartial_m2.toFixed(4));
console.log('  r_partial(DQ, log_m2 | logVflat) = ' + rPartial_logM2.toFixed(4));

const mbarResid_dq = ols(validR.map(r => [Math.log10(r.Vflat), r.logMbar]), dqArr).residuals;
const mbarResid_m2 = ols(validR.map(r => [Math.log10(r.Vflat), r.logMbar]), m2Arr).residuals;
const rPartial_2 = pearsonR(mbarResid_dq, mbarResid_m2);
console.log('  r_partial(DQ, m2 | logVflat, logMbar) = ' + rPartial_2.toFixed(4));


console.log('\n\n' + '#'.repeat(72));
console.log('903.4 — MATCHED PAIR DIRECT COMPARISON');
console.log('#'.repeat(72));

const byName = {};
for (const r of results) byName[r.name] = r;

const pairs = [
  { high: 'NGC2841', low: 'NGC5055', label: 'GOLD PAIR' },
  { high: 'NGC3521', low: 'NGC6946', label: 'PAIR 2' },
  { high: 'NGC3521', low: 'NGC5055', label: 'PAIR 3' },
];

for (const pair of pairs) {
  const h = byName[pair.high], l = byName[pair.low];
  if (!h || !l) { console.log('\n  ' + pair.label + ': missing data'); continue; }
  console.log('\n  ' + pair.label + ': ' + pair.high + ' (DQ=' + h.DQ.toFixed(2) + ') vs ' + pair.low + ' (DQ=' + l.DQ.toFixed(2) + ')');
  console.log('    m2 power: ' + h.m2Power.toFixed(0) + ' vs ' + l.m2Power.toFixed(0) + ' (ratio ' + (h.m2Power / Math.max(l.m2Power, 1)).toFixed(1) + 'x)');
  console.log('    m1 power: ' + h.m1Power.toFixed(0) + ' vs ' + l.m1Power.toFixed(0) + ' (ratio ' + (h.m1Power / Math.max(l.m1Power, 1)).toFixed(1) + 'x)');
  console.log('    PA twist: ' + h.paTwist.toFixed(1) + '° vs ' + l.paTwist.toFixed(1) + '°');
  console.log('    coh outer: ' + h.coherenceOuter.toFixed(3) + ' vs ' + l.coherenceOuter.toFixed(3));
  console.log('    High-H m2 > Low-H m2: ' + (h.m2Power > l.m2Power ? 'YES ✓' : 'NO ✗'));
}


console.log('\n\n' + '#'.repeat(72));
console.log('903.5 — LEAVE-ONE-OUT STABILITY');
console.log('#'.repeat(72));

const looResults = [];
for (let i = 0; i < validR.length; i++) {
  const subDQ = dqArr.filter((_, j) => j !== i);
  const subM2 = m2Arr.filter((_, j) => j !== i);
  const subLogM2 = logM2.filter((_, j) => j !== i);
  looResults.push({ left: validR[i].name, r_m2: pearsonR(subDQ, subM2), r_logM2: pearsonR(subDQ, subLogM2) });
}

console.log('\n  LOO r(DQ, m2):');
for (const l of looResults) {
  console.log('    leave ' + l.left.padEnd(12) + 'r_m2=' + l.r_m2.toFixed(4) + '  r_logM2=' + l.r_logM2.toFixed(4));
}
const meanLOO_m2 = looResults.reduce((s, l) => s + l.r_m2, 0) / looResults.length;
const meanLOO_logM2 = looResults.reduce((s, l) => s + l.r_logM2, 0) / looResults.length;
const minLOO = looResults.reduce((m, l) => Math.min(m, l.r_m2), Infinity);
console.log('\n  Mean LOO r(DQ, m2) = ' + meanLOO_m2.toFixed(4));
console.log('  Mean LOO r(DQ, logM2) = ' + meanLOO_logM2.toFixed(4));
console.log('  Min LOO r(DQ, m2) = ' + minLOO.toFixed(4));
console.log('  All positive: ' + (looResults.every(l => l.r_m2 > 0) ? 'YES' : 'NO'));


console.log('\n\n' + '#'.repeat(72));
console.log('903.6 — PHASE 903 VERDICT');
console.log('#'.repeat(72));

const isSignificant = pVal_m2 < 0.05 || pVal_logM2 < 0.05;
const isStrong = Math.abs(rDQ_m2) > 0.5 || Math.abs(rDQ_logM2) > 0.5;
const isStable = minLOO > 0;
const pairsConfirm = pairs.filter(p => byName[p.high] && byName[p.low] && byName[p.high].m2Power > byName[p.low].m2Power).length;

console.log('\n  CRITERIA CHECK:');
console.log('  1. r(DQ, m2) significant by permutation (p<0.05): ' + (isSignificant ? 'PASS' : 'FAIL'));
console.log('  2. r(DQ, m2) strong (|r| > 0.5): ' + (isStrong ? 'PASS' : 'FAIL'));
console.log('  3. LOO stable (all positive): ' + (isStable ? 'PASS' : 'FAIL'));
console.log('  4. Matched pairs confirm (' + pairsConfirm + '/' + pairs.length + '): ' + (pairsConfirm >= 2 ? 'PASS' : 'FAIL'));

const passed = [isSignificant, isStrong, isStable, pairsConfirm >= 2].filter(Boolean).length;

if (passed >= 3) {
  console.log('\n  PHASE 903 VERDICT: ***PASS*** (' + passed + '/4 criteria met)');
  console.log('  The m=2 azimuthal power in THINGS velocity fields is a robust');
  console.log('  2D discriminant of the hidden state H.');
  console.log('  PROVISIONAL CARRIER: Triaxial/oval halo structure.');
} else if (passed >= 2) {
  console.log('\n  PHASE 903 VERDICT: PARTIAL PASS (' + passed + '/4 criteria met)');
  console.log('  Suggestive but not decisive. Need larger sample or better resolution.');
} else {
  console.log('\n  PHASE 903 VERDICT: FAIL (' + passed + '/4 criteria met)');
  console.log('  2D velocity fields do not robustly distinguish High-H from Low-H.');
}


const outPath = path.join(__dirname, '..', 'public', 'program9-phase903.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: 9, phase: 903,
  title: 'Decisive Matched IFU Test',
  timestamp: new Date().toISOString(),
  N: validR.length,
  correlations: { DQ_m2: rDQ_m2, DQ_logM2: rDQ_logM2, DQ_m1: rDQ_m1, DQ_twist: rDQ_twist, DQ_coh: rDQ_coh },
  partialCorrelations: { DQ_m2_logVf: rPartial_m2, DQ_logM2_logVf: rPartial_logM2, DQ_m2_logVf_logMbar: rPartial_2 },
  permutation: { nPerm, pVal_m2, pVal_logM2 },
  loo: { mean_m2: meanLOO_m2, mean_logM2: meanLOO_logM2, min_m2: minLOO, allPositive: looResults.every(l => l.r_m2 > 0) },
  pairsConfirm,
  passed,
  verdict: passed >= 3 ? 'PASS' : passed >= 2 ? 'PARTIAL' : 'FAIL',
  results,
}, null, 2));
console.log('\nSaved: ' + outPath);
