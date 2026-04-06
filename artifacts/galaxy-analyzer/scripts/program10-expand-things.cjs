const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');
const p903 = require('../public/program9-phase903.json');

function normalize(n) {
  return n.replace(/\s+/g, '').replace(/^NGC0*/, 'NGC').replace(/^DDO0*/, 'DDO').replace(/^IC0*/, 'IC').replace(/^UGC0*/, 'UGC').replace(/^PGC0*/, 'PGC').toUpperCase();
}
function pearsonR(x, y) {
  const n = x.length; if (n < 3) return NaN;
  const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
  return (sxx > 0 && syy > 0) ? sxy / Math.sqrt(sxx * syy) : 0;
}
function ols(X, y) {
  const n = y.length, p = X[0].length;
  const my = y.reduce((a, b) => a + b, 0) / n;
  const mX = Array(p).fill(0);
  for (let j = 0; j < p; j++) { for (let i = 0; i < n; i++) mX[j] += X[i][j]; mX[j] /= n; }
  const Xc = X.map(r => r.map((v, j) => v - mX[j])), yc = y.map(v => v - my);
  const XtX = Array.from({ length: p }, () => Array(p).fill(0)), Xty = Array(p).fill(0);
  for (let i = 0; i < n; i++) for (let j = 0; j < p; j++) { Xty[j] += Xc[i][j] * yc[i]; for (let k = 0; k < p; k++) XtX[j][k] += Xc[i][j] * Xc[i][k]; }
  const aug = XtX.map((r, i) => [...r, Xty[i]]);
  for (let c = 0; c < p; c++) { let mr = c; for (let r = c + 1; r < p; r++) if (Math.abs(aug[r][c]) > Math.abs(aug[mr][c])) mr = r; [aug[c], aug[mr]] = [aug[mr], aug[c]]; if (Math.abs(aug[c][c]) < 1e-14) continue; for (let r = c + 1; r < p; r++) { const f = aug[r][c] / aug[c][c]; for (let j = c; j <= p; j++) aug[r][j] -= f * aug[c][j]; } }
  const beta = Array(p).fill(0);
  for (let i = p - 1; i >= 0; i--) { beta[i] = aug[i][p]; for (let j = i + 1; j < p; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= Math.abs(aug[i][i]) > 1e-14 ? aug[i][i] : 1; }
  const res = []; for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < p; j++) pred += beta[j] * Xc[i][j]; res.push(y[i] - pred); }
  return res;
}
function zscore(arr) {
  const m = arr.reduce((a, b) => a + b, 0) / arr.length;
  const s = Math.sqrt(arr.reduce((ss, v) => ss + (v - m) ** 2, 0) / arr.length);
  return s > 1e-10 ? arr.map(v => (v - m) / s) : arr.map(() => 0);
}
function permTest(x, y, nPerm) {
  const rObs = Math.abs(pearsonR(x, y));
  let count = 0;
  for (let p = 0; p < nPerm; p++) {
    const yp = y.slice();
    for (let i = yp.length - 1; i > 0; i--) { const j = Math.floor(Math.random() * (i + 1)); [yp[i], yp[j]] = [yp[j], yp[i]]; }
    if (Math.abs(pearsonR(x, yp)) >= rObs) count++;
  }
  return count / nPerm;
}


console.log('='.repeat(72));
console.log('PROGRAM 10.2 — EXPAND THINGS SAMPLE');
console.log('Include all SPARC galaxies with THINGS FITS data');
console.log('='.repeat(72));


const sparcMap = {}; sparc.forEach(g => { sparcMap[normalize(g.name)] = g; });
const resultsMap = {}; sparcResults.perGalaxy.forEach(g => { resultsMap[normalize(g.name)] = g; });
const d56Map = {}; d56.perGalaxy.forEach(g => { d56Map[normalize(g.name)] = g; });


const thingsDir = path.join(__dirname, '..', 'data', 'things-2d');
const fitsFiles = fs.readdirSync(thingsDir).filter(f => f.endsWith('.FITS') && f.includes('MOM1'));

console.log('\n  THINGS MOM1 FITS on disk: ' + fitsFiles.length);


console.log('\n\n' + '#'.repeat(72));
console.log('10.2.1 — CHECK ALL FITS AGAINST FULL SPARC (N=175)');
console.log('#'.repeat(72));

const fitsGalaxies = [];

for (const f of fitsFiles) {
  const rawName = f.replace('_MOM1.FITS', '');
  const nn = normalize(rawName);

  const sp = sparcMap[nn];
  const sr = resultsMap[nn];
  const inD56 = d56Map[nn];

  console.log('\n  ' + rawName.padEnd(12));
  if (sp) {
    console.log('    SPARC: YES  Vflat=' + sp.Vflat + '  D=' + sp.D + '  T=' + sp.T + '  Q=' + sp.Q + '  Inc=' + sp.inc);
    if (sr) {
      console.log('    Results: YES  a0=' + (sr.a0 || 'N/A'));
    } else {
      console.log('    Results: NO (no RAR fit available)');
    }
    if (inD56) {
      console.log('    In d56 (published quality): YES  logA0=' + inD56.logA0.toFixed(4));
    } else {
      console.log('    In d56 (published quality): NO');
      if (sp.Vflat > 0 && sr) {
        console.log('    >>> CANDIDATE FOR EXPANSION');
      }
    }
    fitsGalaxies.push({ name: rawName, nn, sp, sr, inD56: !!inD56 });
  } else {
    console.log('    SPARC: NO — not in database');
  }
}


console.log('\n\n' + '#'.repeat(72));
console.log('10.2.2 — BUILD EXPANDED BASELINE');
console.log('#'.repeat(72));

const expandedGals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[normalize(g.name)]; if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[normalize(g.name)]; if (!sr) continue;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  expandedGals.push({
    name: g.name, logVflat: Math.log10(sp.Vflat), Vflat: sp.Vflat,
    logMbar: Math.log10(Math.max(Mbar, 1)), logL36: Math.log10(Math.max(sp.L36, 0.001)),
    logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)), morphT: sp.T,
    logMHI: g.logMHI, logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    logA0: g.logA0, dist: sp.D, inc: sp.inc || 0,
    inD56: true,
  });
}

for (const fg of fitsGalaxies) {
  if (fg.inD56) continue;
  if (!fg.sp || fg.sp.Vflat <= 0) continue;
  if (!fg.sr) continue;
  if (!fg.sr.a0 || fg.sr.a0 <= 0) continue;
  const sp = fg.sp;
  const sr = fg.sr;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  if (Mbar <= 0) continue;
  const logA0 = Math.log10(sr.a0);
  if (!isFinite(logA0)) continue;
  const logMHI = Math.log10(Math.max(sp.MHI * 1e9, 1));

  expandedGals.push({
    name: fg.name, logVflat: Math.log10(sp.Vflat), Vflat: sp.Vflat,
    logMbar: Math.log10(Math.max(Mbar, 1)), logL36: Math.log10(Math.max(sp.L36, 0.001)),
    logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)), morphT: sp.T,
    logMHI: logMHI, logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    logA0: logA0, dist: sp.D, inc: sp.inc || 0,
    inD56: false,
  });
  console.log('\n  ADDED: ' + fg.name + '  Vflat=' + sp.Vflat + '  logA0=' + logA0.toFixed(4) + '  D=' + sp.D);
}

console.log('\n  Expanded baseline: N = ' + expandedGals.length + ' (was ' + d56.perGalaxy.length + ')');

const vfR = ols(expandedGals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]), expandedGals.map(g => g.logVflat));
const a0R = ols(expandedGals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]), expandedGals.map(g => g.logA0));
const bilZ = zscore(zscore(vfR).map((v, i) => v + zscore(a0R)[i]));
for (let i = 0; i < expandedGals.length; i++) {
  expandedGals[i].DQ = bilZ[i]; expandedGals[i].VfResid = vfR[i]; expandedGals[i].a0Resid = a0R[i];
}
const galMap = {}; expandedGals.forEach(g => { galMap[normalize(g.name)] = g; });

const rBilateral = pearsonR(vfR, a0R);
console.log('  r(VfResid, a0Resid) on expanded baseline: ' + rBilateral.toFixed(4));


console.log('\n\n' + '#'.repeat(72));
console.log('10.2.3 — PROCESS ALL THINGS FITS');
console.log('#'.repeat(72));

function readFitsMom1(filePath) {
  const buf = fs.readFileSync(filePath);
  let headerEnd = -1;
  for (let i = 0; i < buf.length - 3; i++) {
    if (buf[i] === 0x45 && buf[i+1] === 0x4E && buf[i+2] === 0x44 && buf[i+3] === 0x20) { headerEnd = i; break; }
  }
  if (headerEnd < 0) throw new Error('No END in FITS header');
  const headerStr = buf.slice(0, headerEnd + 80).toString('ascii');
  const cards = {};
  for (let i = 0; i < headerStr.length; i += 80) {
    const card = headerStr.slice(i, i + 80);
    const key = card.slice(0, 8).trim();
    if (card[8] === '=' && key) {
      let valStr = card.slice(10, 30).trim();
      if (valStr.startsWith("'")) valStr = valStr.replace(/'/g, '').trim();
      const num = parseFloat(valStr); cards[key] = isNaN(num) ? valStr : num;
    }
  }
  const naxis1 = cards.NAXIS1 || 0, naxis2 = cards.NAXIS2 || 0, bitpix = cards.BITPIX || -32;
  const bscale = cards.BSCALE || 1, bzero = cards.BZERO || 0;
  const crpix1 = cards.CRPIX1 || naxis1/2, crpix2 = cards.CRPIX2 || naxis2/2;
  const cdelt1 = cards.CDELT1 || 1;
  const dataStart = Math.ceil((headerEnd + 80) / 2880) * 2880;
  const bytesPerPix = Math.abs(bitpix) / 8;
  const pixels = [];
  for (let j = 0; j < naxis2; j++) {
    for (let i = 0; i < naxis1; i++) {
      const pos = dataStart + (j * naxis1 + i) * bytesPerPix;
      if (pos + bytesPerPix > buf.length) continue;
      let val;
      if (bitpix === -32) val = buf.readFloatBE(pos);
      else if (bitpix === -64) val = buf.readDoubleBE(pos);
      else if (bitpix === 16) val = buf.readInt16BE(pos);
      else if (bitpix === 32) val = buf.readInt32BE(pos);
      else continue;
      val = val * bscale + bzero;
      if (!isFinite(val) || Math.abs(val) > 1e6) continue;
      pixels.push({ x: i, y: j, val });
    }
  }
  return { naxis1, naxis2, crpix1, crpix2, cdelt1, pixels, cards };
}

function computeM2Power(fits, vsys) {
  const { pixels, crpix1, crpix2, cdelt1 } = fits;
  const pixScale = Math.abs(cdelt1) * 3600;
  const cx = crpix1 - 1, cy = crpix2 - 1;
  const nBins = 15;
  const rads = pixels.map(p => Math.sqrt((p.x - cx) ** 2 + (p.y - cy) ** 2) * pixScale);
  const maxR = rads.reduce((a, b) => a > b ? a : b, 0);
  const binW = maxR / nBins;
  if (binW < 1) return null;
  const bins = Array.from({ length: nBins }, () => ({ a0: 0, a1c: 0, a1s: 0, a2c: 0, a2s: 0, n: 0 }));
  for (let k = 0; k < pixels.length; k++) {
    const p = pixels[k], dx = p.x - cx, dy = p.y - cy;
    const bi = Math.min(Math.floor(rads[k] / binW), nBins - 1);
    const theta = Math.atan2(dy, dx), vCorr = p.val - vsys;
    bins[bi].a0 += Math.abs(vCorr);
    bins[bi].a1c += vCorr * Math.cos(theta); bins[bi].a1s += vCorr * Math.sin(theta);
    bins[bi].a2c += vCorr * Math.cos(2 * theta); bins[bi].a2s += vCorr * Math.sin(2 * theta);
    bins[bi].n++;
  }
  let m0S = 0, m1S = 0, m2S = 0, nU = 0;
  const paInner = [], paOuter = [];
  for (let i = 0; i < nBins; i++) {
    if (bins[i].n < 10) continue;
    const n = bins[i].n;
    m0S += (bins[i].a0 / n) ** 2;
    m1S += ((bins[i].a1c / n) ** 2 + (bins[i].a1s / n) ** 2);
    m2S += ((bins[i].a2c / n) ** 2 + (bins[i].a2s / n) ** 2);
    nU++;
    const pa = 0.5 * Math.atan2(bins[i].a2s / n, bins[i].a2c / n) * 180 / Math.PI;
    if (i < nBins / 2) paInner.push(pa); else paOuter.push(pa);
  }
  if (nU === 0) return null;
  const meanI = paInner.length > 0 ? paInner.reduce((a, b) => a + b, 0) / paInner.length : 0;
  const meanO = paOuter.length > 0 ? paOuter.reduce((a, b) => a + b, 0) / paOuter.length : 0;
  const outerBins = bins.slice(Math.floor(nBins * 0.5));
  let coh = 0, totO = 0;
  const refSign = outerBins.find(b => b.n >= 10)?.a2c || 1;
  for (const b of outerBins) { if (b.n < 10) continue; totO++; if (Math.sign(b.a2c) === Math.sign(refSign)) coh++; }
  return { m0: m0S / nU, m1: m1S / nU, m2: m2S / nU, m2m0: m0S > 0 ? (m2S / nU) / (m0S / nU) : 0,
    paTwist: Math.abs(meanO - meanI), coherence: totO > 0 ? coh / totO : 0, nValid: pixels.length };
}


const allResults = [];
const p903Names = new Set(p903.results.map(r => normalize(r.name)));

for (const f of fitsFiles) {
  const rawName = f.replace('_MOM1.FITS', '');
  const nn = normalize(rawName);
  const g = galMap[nn];

  if (!g) {
    console.log('\n  SKIP: ' + rawName + ' (not in expanded SPARC baseline)');
    continue;
  }

  const fromP9 = p903Names.has(nn);
  if (fromP9) {
    const orig = p903.results.find(r => normalize(r.name) === nn);
    allResults.push({
      name: g.name, DQ: g.DQ, Vflat: g.Vflat, dist: g.dist, logMbar: g.logMbar, inc: g.inc,
      m2Power: orig.m2Power, m1Power: orig.m1Power, m0Power: orig.m0Power, m2m0: orig.m2m0,
      paTwist: orig.paTwist, coherenceOuter: orig.coherenceOuter, nValid: orig.nValid,
      survey: 'THINGS', source: 'Program9', inD56: g.inD56,
    });
    console.log('\n  ' + g.name.padEnd(12) + ' DQ=' + g.DQ.toFixed(2).padEnd(8) + ' [Program 9 — reusing]');
    continue;
  }

  console.log('\n  Processing NEW: ' + rawName + '  DQ=' + g.DQ.toFixed(2) + '  Vflat=' + g.Vflat);
  try {
    const fitsPath = path.join(thingsDir, f);
    const fits = readFitsMom1(fitsPath);
    const vals = fits.pixels.map(p => p.val);
    const vsys = vals.reduce((a, b) => a + b, 0) / vals.length;
    const m2data = computeM2Power(fits, vsys);
    if (!m2data) { console.log('    SKIP: insufficient bins'); continue; }
    allResults.push({
      name: g.name, DQ: g.DQ, Vflat: g.Vflat, dist: g.dist, logMbar: g.logMbar, inc: g.inc,
      m2Power: m2data.m2, m1Power: m2data.m1, m0Power: m2data.m0, m2m0: m2data.m2m0,
      paTwist: m2data.paTwist, coherenceOuter: m2data.coherence, nValid: m2data.nValid,
      survey: 'THINGS', source: 'Program10-expanded', inD56: g.inD56,
    });
    console.log('    m2=' + m2data.m2.toFixed(1) + '  m2/m0=' + m2data.m2m0.toFixed(3) + '  twist=' + m2data.paTwist.toFixed(1) + '  nPix=' + m2data.nValid);
  } catch(e) {
    console.log('    ERROR: ' + e.message);
  }
}


console.log('\n\n' + '#'.repeat(72));
console.log('10.2.4 — EXPANDED CORRELATION ANALYSIS');
console.log('#'.repeat(72));

const valid = allResults.filter(r => r.m2Power > 0 && isFinite(r.m2Power) && isFinite(r.DQ));
const fromP9 = valid.filter(r => r.source === 'Program9');
const newGals = valid.filter(r => r.source !== 'Program9');

console.log('\n  Total valid galaxies: ' + valid.length);
console.log('  From Program 9: ' + fromP9.length);
console.log('  NEW (Program 10): ' + newGals.length);
console.log('  From d56 baseline: ' + valid.filter(r => r.inD56).length);
console.log('  Expanded beyond d56: ' + valid.filter(r => !r.inD56).length);

if (valid.length >= 4) {
  const dqs = valid.map(r => r.DQ), logM2s = valid.map(r => Math.log10(r.m2Power));

  const rAll = pearsonR(dqs, logM2s);
  const pAll = permTest(dqs, logM2s, 10000);

  console.log('\n  FULL SAMPLE (N=' + valid.length + '):');
  console.log('    r(DQ, log m2) = ' + rAll.toFixed(4));
  console.log('    p-value = ' + pAll.toFixed(4));

  if (fromP9.length >= 4) {
    const rP9 = pearsonR(fromP9.map(r => r.DQ), fromP9.map(r => Math.log10(r.m2Power)));
    console.log('\n  PROGRAM 9 SUBSET (N=' + fromP9.length + '):');
    console.log('    r(DQ, log m2) = ' + rP9.toFixed(4));
  }

  if (newGals.length >= 3) {
    const rNew = pearsonR(newGals.map(r => r.DQ), newGals.map(r => Math.log10(r.m2Power)));
    console.log('\n  NEW GALAXIES ONLY (N=' + newGals.length + '):');
    console.log('    r(DQ, log m2) = ' + rNew.toFixed(4));
  }

  console.log('\n  Per-galaxy:');
  console.log('  ' + 'Galaxy'.padEnd(14) + 'DQ'.padEnd(8) + 'log(m2)'.padEnd(10) + 'Vflat'.padEnd(8) + 'Source'.padEnd(18) + 'In d56');
  console.log('  ' + '-'.repeat(70));
  for (const r of valid.sort((a, b) => b.DQ - a.DQ)) {
    console.log('  ' + r.name.padEnd(14) + r.DQ.toFixed(2).padEnd(8) + Math.log10(r.m2Power).toFixed(3).padEnd(10) + Math.round(r.Vflat).toString().padEnd(8) + r.source.padEnd(18) + (r.inD56 ? 'YES' : 'NO'));
  }

  console.log('\n  LOO:');
  const looRs = [];
  for (let i = 0; i < valid.length; i++) {
    const xL = dqs.filter((_, j) => j !== i), yL = logM2s.filter((_, j) => j !== i);
    looRs.push(pearsonR(xL, yL));
  }
  const allPos = looRs.every(r => r > 0);
  console.log('    All positive: ' + allPos + '  (' + looRs.filter(r => r > 0).length + '/' + looRs.length + ')');
  console.log('    Mean: ' + (looRs.reduce((a, b) => a + b, 0) / looRs.length).toFixed(4));
  console.log('    Range: [' + Math.min(...looRs).toFixed(4) + ', ' + Math.max(...looRs).toFixed(4) + ']');

  for (let i = 0; i < valid.length; i++) {
    const dropped = valid[i].name;
    console.log('      Drop ' + dropped.padEnd(12) + ' => r=' + looRs[i].toFixed(4) + '  (' + (looRs[i] > rAll ? '+' : '') + (looRs[i] - rAll).toFixed(4) + ')' + (valid[i].source !== 'Program9' ? '  [NEW]' : ''));
  }

  const d56Only = valid.filter(r => r.inD56);
  if (d56Only.length >= 4 && d56Only.length < valid.length) {
    const rD56 = pearsonR(d56Only.map(r => r.DQ), d56Only.map(r => Math.log10(r.m2Power)));
    console.log('\n  d56-only subset (N=' + d56Only.length + '): r = ' + rD56.toFixed(4));
  }

  const result = {
    program: 10, phase: '10.2', title: 'Expanded THINGS Sample',
    timestamp: new Date().toISOString(),
    N_total: valid.length, N_program9: fromP9.length, N_new: newGals.length,
    N_d56: valid.filter(r => r.inD56).length, N_expanded: valid.filter(r => !r.inD56).length,
    correlation: { r: rAll, pValue: pAll, r_program9: 0.847 },
    loo: { allPositive: allPos, mean: looRs.reduce((a, b) => a + b, 0) / looRs.length,
      min: Math.min(...looRs), max: Math.max(...looRs) },
    results: valid.sort((a, b) => b.DQ - a.DQ),
    expandedBaseline: { N: expandedGals.length, rBilateral },
  };
  fs.writeFileSync(path.join(__dirname, '..', 'public', 'program10-expanded.json'), JSON.stringify(result, null, 2));
  console.log('\n  Saved: program10-expanded.json');
} else {
  console.log('\n  Insufficient data (N=' + valid.length + ')');
}

console.log('\n' + '='.repeat(72));
console.log('PROGRAM 10.2 COMPLETE');
console.log('='.repeat(72));
