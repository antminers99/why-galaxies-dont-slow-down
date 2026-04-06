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


const sparcMap = {}; sparc.forEach(g => { sparcMap[normalize(g.name)] = g; });
const resultsMap = {}; sparcResults.perGalaxy.forEach(g => { resultsMap[normalize(g.name)] = g; });

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[normalize(g.name)]; if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[normalize(g.name)]; if (!sr) continue;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  gals.push({
    name: g.name, logVflat: Math.log10(sp.Vflat), Vflat: sp.Vflat,
    logMbar: Math.log10(Math.max(Mbar, 1)), logL36: Math.log10(Math.max(sp.L36, 0.001)),
    logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)), morphT: sp.T,
    logMHI: g.logMHI, logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    logA0: g.logA0, dist: sp.D, inc: sp.Inc || 0,
  });
}
const vfR = ols(gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]), gals.map(g => g.logVflat));
const a0R = ols(gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]), gals.map(g => g.logA0));
const bilZ = zscore(zscore(vfR).map((v, i) => v + zscore(a0R)[i]));
for (let i = 0; i < gals.length; i++) { gals[i].DQ = bilZ[i]; gals[i].VfResid = vfR[i]; gals[i].a0Resid = a0R[i]; }
const galMap = {}; gals.forEach(g => { galMap[normalize(g.name)] = g; });


const dataDir = path.join(__dirname, '..', 'data');
const thingsDir = path.join(dataDir, 'things-2d');
const multiDir = path.join(dataDir, 'multi-survey-2d');
fs.mkdirSync(multiDir, { recursive: true });


console.log('='.repeat(72));
console.log('PROGRAM 10 — MULTI-SURVEY 2D CARRIER VERIFICATION');
console.log('From observational carrier to physical identification');
console.log('='.repeat(72));


console.log('\n' + '#'.repeat(72));
console.log('10.1 — EXPANDED SURVEY CENSUS');
console.log('#'.repeat(72));

const surveyData = {
  THINGS: { dir: thingsDir, galaxies: new Map(), type: 'HI-21cm', resolution: '~6-12"', fov: 'full' },
  MaNGA: { dir: path.join(multiDir, 'manga'), galaxies: new Map(), type: 'optical-IFU', resolution: '~2.5"', fov: '1-2 Re' },
  CALIFA: { dir: path.join(multiDir, 'califa'), galaxies: new Map(), type: 'optical-IFU', resolution: '~2.5"', fov: '2.5 Re' },
  PHANGS: { dir: path.join(multiDir, 'phangs'), galaxies: new Map(), type: 'ALMA-CO/MUSE', resolution: '~1"', fov: 'inner disc' },
};

const thingsSurvey = ['NGC628','NGC925','NGC2366','NGC2403','NGC2841','NGC2903','NGC2976','NGC3031','NGC3184','NGC3198','NGC3351','NGC3521','NGC3621','NGC3627','NGC4214','NGC4449','NGC4736','NGC4826','NGC5055','NGC5194','NGC5236','NGC5457','NGC6946','NGC7331','NGC7793','DDO50','DDO53','DDO63','DDO154','DDO165','IC2574','NGC1569','HOII'].map(normalize);
const mangaGals = ['NGC2841','NGC3198','NGC3521','NGC4051','NGC5055','NGC7331','UGC3546','UGC6786','UGC9037','NGC2403','NGC3031','NGC4088','NGC4559','UGC2953','NGC4826','IC2574'].map(normalize);
const califaGals = ['NGC2841','NGC3521','NGC5055','NGC7331','NGC4826','NGC6946','UGC3580','NGC3198','NGC2403','NGC3031'].map(normalize);
const phangsGals = ['NGC628','NGC1087','NGC1300','NGC1365','NGC1385','NGC1433','NGC1512','NGC1566','NGC1672','NGC2835','NGC2903','NGC3351','NGC3521','NGC3627','NGC4254','NGC4303','NGC4321','NGC4535','NGC4548','NGC4571','NGC4654','NGC4689','NGC4826','NGC5068','NGC5248','NGC6744','NGC7496'].map(normalize);

const allTargets = new Map();

function addTarget(name, survey) {
  const nn = normalize(name);
  const g = galMap[nn];
  if (!g) return;
  if (!allTargets.has(nn)) {
    allTargets.set(nn, { name: g.name, DQ: g.DQ, Vflat: g.Vflat, dist: g.dist, logMbar: g.logMbar, inc: g.inc, surveys: new Set() });
  }
  allTargets.get(nn).surveys.add(survey);
}

thingsSurvey.forEach(n => addTarget(n, 'THINGS'));
mangaGals.forEach(n => addTarget(n, 'MaNGA'));
califaGals.forEach(n => addTarget(n, 'CALIFA'));
phangsGals.forEach(n => addTarget(n, 'PHANGS'));

const targetList = Array.from(allTargets.values()).sort((a, b) => b.DQ - a.DQ);

console.log('\nSPARC galaxies with ANY 2D survey coverage: ' + targetList.length);
console.log('\n  ' + 'Galaxy'.padEnd(15) + 'DQ'.padEnd(8) + 'Vflat'.padEnd(8) + 'Dist'.padEnd(8) + 'Surveys'.padEnd(30) + 'Status');
console.log('  ' + '-'.repeat(85));

const thingsProcessed = new Set(p903.results.map(r => normalize(r.name)));
let newTargets = 0;
const priorityTargets = [];

for (const t of targetList) {
  const nn = normalize(t.name);
  const survList = Array.from(t.surveys).join(', ');
  const processed = thingsProcessed.has(nn);
  const status = processed ? 'DONE (Program 9)' : 'NEW TARGET';
  if (!processed) {
    newTargets++;
    priorityTargets.push(t);
  }
  const marker = Math.abs(t.DQ) > 1.0 ? (t.DQ > 0 ? ' ***' : ' ---') : '';
  console.log('  ' + t.name.padEnd(15) + t.DQ.toFixed(2).padEnd(8) + Math.round(t.Vflat).toString().padEnd(8) + t.dist.toFixed(1).padEnd(8) + survList.padEnd(30) + status + marker);
}

console.log('\n  Already processed (THINGS, Program 9): ' + thingsProcessed.size);
console.log('  New targets with 2D coverage: ' + newTargets);
console.log('  Total available for Program 10: ' + (thingsProcessed.size + newTargets));


console.log('\n\n' + '#'.repeat(72));
console.log('10.1.1 — MULTI-SURVEY STRATEGY');
console.log('#'.repeat(72));

const thingsOnly = [];
const mangaOnly = [];
const califaOnly = [];
const multiSurvey = [];
const allNew = [];

for (const t of priorityTargets) {
  const nn = normalize(t.name);
  const s = t.surveys;
  if (s.size > 1) multiSurvey.push(t);
  if (s.has('THINGS') && !thingsProcessed.has(nn)) thingsOnly.push(t);
  if (s.has('MaNGA')) mangaOnly.push(t);
  if (s.has('CALIFA')) califaOnly.push(t);
  allNew.push(t);
}

console.log('\n  Survey breakdown of NEW targets:');
console.log('    With THINGS data (not yet processed): ' + thingsOnly.length);
console.log('    With MaNGA data: ' + mangaOnly.length);
console.log('    With CALIFA data: ' + califaOnly.length);
console.log('    With multi-survey overlap: ' + multiSurvey.length);

console.log('\n  KEY: Multi-survey overlap galaxies (cross-validation):');
for (const t of multiSurvey) {
  console.log('    ' + t.name.padEnd(15) + 'DQ=' + t.DQ.toFixed(2).padEnd(8) + Array.from(t.surveys).join(', '));
}


console.log('\n\n' + '#'.repeat(72));
console.log('10.1.2 — EXISTING THINGS FITS: EXPAND TO ALL AVAILABLE');
console.log('#'.repeat(72));

const fitsFiles = fs.existsSync(thingsDir) ? fs.readdirSync(thingsDir).filter(f => f.endsWith('.FITS') && f.includes('MOM1')) : [];
console.log('\n  THINGS MOM1 FITS files on disk: ' + fitsFiles.length);
for (const f of fitsFiles) {
  const nn = normalize(f.replace('_MOM1.FITS', ''));
  const g = galMap[nn];
  const processed = thingsProcessed.has(nn);
  console.log('    ' + f.padEnd(25) + (g ? 'SPARC: YES  DQ=' + g.DQ.toFixed(2) : 'SPARC: NO') + (processed ? '  [PROCESSED]' : '  [NOT YET PROCESSED]'));
}

const unprocessedThings = fitsFiles.map(f => normalize(f.replace('_MOM1.FITS', ''))).filter(nn => galMap[nn] && !thingsProcessed.has(nn));
console.log('\n  THINGS FITS available but NOT processed: ' + unprocessedThings.length);
for (const nn of unprocessedThings) {
  const g = galMap[nn];
  if (g) console.log('    ' + g.name.padEnd(15) + 'DQ=' + g.DQ.toFixed(2));
}


console.log('\n\n' + '#'.repeat(72));
console.log('10.1.3 — PROCESS UNPROCESSED THINGS GALAXIES');
console.log('#'.repeat(72));

let fitsReader;
try { fitsReader = require('../src/fits-reader.cjs'); } catch(e) { fitsReader = null; }
if (!fitsReader) {
  try { fitsReader = require('./fits-utils.cjs'); } catch(e) { fitsReader = null; }
}

function readFitsMom1(filePath) {
  const buf = fs.readFileSync(filePath);
  let headerEnd = -1;
  for (let i = 0; i < buf.length - 3; i++) {
    if (buf[i] === 0x45 && buf[i+1] === 0x4E && buf[i+2] === 0x44 && buf[i+3] === 0x20) {
      headerEnd = i;
      break;
    }
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
      const num = parseFloat(valStr);
      cards[key] = isNaN(num) ? valStr : num;
    }
  }

  const naxis1 = cards.NAXIS1 || 0, naxis2 = cards.NAXIS2 || 0;
  const bitpix = cards.BITPIX || -32;
  const bscale = cards.BSCALE || 1, bzero = cards.BZERO || 0;
  const crpix1 = cards.CRPIX1 || naxis1/2, crpix2 = cards.CRPIX2 || naxis2/2;
  const cdelt1 = cards.CDELT1 || 1, cdelt2 = cards.CDELT2 || 1;
  const crval1 = cards.CRVAL1 || 0, crval2 = cards.CRVAL2 || 0;

  const dataStart = Math.ceil((headerEnd + 80) / 2880) * 2880;
  const naxis3 = cards.NAXIS3 || 1;
  const planeSize = naxis1 * naxis2;
  const bytesPerPix = Math.abs(bitpix) / 8;
  const offset = dataStart;

  const pixels = [];
  for (let j = 0; j < naxis2; j++) {
    for (let i = 0; i < naxis1; i++) {
      const idx = j * naxis1 + i;
      const pos = offset + idx * bytesPerPix;
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

  return { naxis1, naxis2, crpix1, crpix2, cdelt1, cdelt2, crval1, crval2, pixels, cards };
}

function computeM2Power(fits, vsys) {
  const { pixels, crpix1, crpix2, cdelt1 } = fits;
  const pixScale = Math.abs(cdelt1) * 3600;

  const cx = crpix1 - 1, cy = crpix2 - 1;
  const radii = pixels.map(p => Math.sqrt((p.x - cx) ** 2 + (p.y - cy) ** 2) * pixScale);
  const maxR = radii.reduce((a, b) => Math.max(a, b), 0);

  const nBins = 15;
  const binW = maxR / nBins;
  if (binW < 1) return { m0: 0, m1: 0, m2: 0, m2m0: 0, paTwist: 0, coherence: 0, nValid: pixels.length };

  const bins = Array.from({ length: nBins }, () => ({ a0: 0, a1c: 0, a1s: 0, a2c: 0, a2s: 0, n: 0 }));

  for (let k = 0; k < pixels.length; k++) {
    const p = pixels[k];
    const dx = p.x - cx, dy = p.y - cy;
    const r = radii[k];
    const bi = Math.min(Math.floor(r / binW), nBins - 1);
    const theta = Math.atan2(dy, dx);
    const vCorr = p.val - vsys;
    bins[bi].a0 += Math.abs(vCorr);
    bins[bi].a1c += vCorr * Math.cos(theta);
    bins[bi].a1s += vCorr * Math.sin(theta);
    bins[bi].a2c += vCorr * Math.cos(2 * theta);
    bins[bi].a2s += vCorr * Math.sin(2 * theta);
    bins[bi].n++;
  }

  let m0Sum = 0, m1Sum = 0, m2Sum = 0, nBinsUsed = 0;
  const paInner = [], paOuter = [];

  for (let i = 0; i < nBins; i++) {
    if (bins[i].n < 10) continue;
    const n = bins[i].n;
    const a0 = bins[i].a0 / n;
    const a1 = Math.sqrt((bins[i].a1c / n) ** 2 + (bins[i].a1s / n) ** 2);
    const a2 = Math.sqrt((bins[i].a2c / n) ** 2 + (bins[i].a2s / n) ** 2);
    m0Sum += a0 ** 2;
    m1Sum += a1 ** 2;
    m2Sum += a2 ** 2;
    nBinsUsed++;

    const pa = 0.5 * Math.atan2(bins[i].a2s / n, bins[i].a2c / n) * 180 / Math.PI;
    if (i < nBins / 2) paInner.push(pa);
    else paOuter.push(pa);
  }

  if (nBinsUsed === 0) return { m0: 0, m1: 0, m2: 0, m2m0: 0, paTwist: 0, coherence: 0, nValid: pixels.length };

  const m0 = m0Sum / nBinsUsed;
  const m1 = m1Sum / nBinsUsed;
  const m2 = m2Sum / nBinsUsed;

  const meanInner = paInner.length > 0 ? paInner.reduce((a, b) => a + b, 0) / paInner.length : 0;
  const meanOuter = paOuter.length > 0 ? paOuter.reduce((a, b) => a + b, 0) / paOuter.length : 0;
  const paTwist = Math.abs(meanOuter - meanInner);

  const outerBins = bins.slice(Math.floor(nBins * 0.5));
  let coherent = 0, totalOuter = 0;
  for (const b of outerBins) {
    if (b.n < 10) continue;
    totalOuter++;
    const sign = Math.sign(b.a2c);
    if (sign === Math.sign(outerBins.find(ob => ob.n >= 10)?.a2c || 1)) coherent++;
  }
  const coherence = totalOuter > 0 ? coherent / totalOuter : 0;

  return { m0, m1, m2, m2m0: m0 > 0 ? m2 / m0 : 0, paTwist, coherence, nValid: pixels.length };
}


const allResults = [];

for (const r of p903.results) {
  allResults.push({
    name: r.name, DQ: r.DQ, Vflat: r.Vflat, dist: r.dist, logMbar: r.logMbar, inc: r.inc,
    m2Power: r.m2Power, m1Power: r.m1Power, m0Power: r.m0Power, m2m0: r.m2m0,
    paTwist: r.paTwist, coherenceOuter: r.coherenceOuter, nValid: r.nValid,
    survey: 'THINGS', source: 'Program9',
  });
}
console.log('\nExisting Program 9 results: ' + allResults.length + ' galaxies');

for (const nn of unprocessedThings) {
  const g = galMap[nn];
  if (!g) continue;
  const fitsFile = path.join(thingsDir, nn.replace('NGC', 'NGC') + '_MOM1.FITS');
  const altFile = path.join(thingsDir, g.name.replace(/\s+/g, '') + '_MOM1.FITS');
  const actualFile = fs.existsSync(fitsFile) ? fitsFile : fs.existsSync(altFile) ? altFile : null;
  if (!actualFile) {
    const candidates = fitsFiles.filter(f => normalize(f.replace('_MOM1.FITS', '')) === nn);
    if (candidates.length > 0) {
      const cf = path.join(thingsDir, candidates[0]);
      if (fs.existsSync(cf)) {
        console.log('\n  Processing: ' + g.name + ' (THINGS, DQ=' + g.DQ.toFixed(2) + ')');
        try {
          const fits = readFitsMom1(cf);
          const vals = fits.pixels.map(p => p.val);
          const vsys = vals.reduce((a, b) => a + b, 0) / vals.length;
          const m2data = computeM2Power(fits, vsys);
          allResults.push({
            name: g.name, DQ: g.DQ, Vflat: g.Vflat, dist: g.dist, logMbar: g.logMbar, inc: g.inc,
            m2Power: m2data.m2, m1Power: m2data.m1, m0Power: m2data.m0, m2m0: m2data.m2m0,
            paTwist: m2data.paTwist, coherenceOuter: m2data.coherence, nValid: m2data.nValid,
            survey: 'THINGS', source: 'Program10-new',
          });
          console.log('    m2=' + m2data.m2.toFixed(1) + '  m2/m0=' + m2data.m2m0.toFixed(3) + '  PA twist=' + m2data.paTwist.toFixed(1) + '°  nPix=' + m2data.nValid);
        } catch(e) {
          console.log('    ERROR: ' + e.message);
        }
      }
    }
    continue;
  }

  console.log('\n  Processing: ' + g.name + ' (THINGS, DQ=' + g.DQ.toFixed(2) + ')');
  try {
    const fits = readFitsMom1(actualFile);
    const vals = fits.pixels.map(p => p.val);
    const vsys = vals.reduce((a, b) => a + b, 0) / vals.length;
    const m2data = computeM2Power(fits, vsys);
    allResults.push({
      name: g.name, DQ: g.DQ, Vflat: g.Vflat, dist: g.dist, logMbar: g.logMbar, inc: g.inc,
      m2Power: m2data.m2, m1Power: m2data.m1, m0Power: m2data.m0, m2m0: m2data.m2m0,
      paTwist: m2data.paTwist, coherenceOuter: m2data.coherence, nValid: m2data.nValid,
      survey: 'THINGS', source: 'Program10-new',
    });
    console.log('    m2=' + m2data.m2.toFixed(1) + '  m2/m0=' + m2data.m2m0.toFixed(3) + '  PA twist=' + m2data.paTwist.toFixed(1) + '°  nPix=' + m2data.nValid);
  } catch(e) {
    console.log('    ERROR: ' + e.message);
  }
}


console.log('\n\n' + '#'.repeat(72));
console.log('10.1.4 — EXPANDED CORRELATION ANALYSIS');
console.log('#'.repeat(72));

const valid = allResults.filter(r => r.m2Power > 0 && isFinite(r.m2Power) && isFinite(r.DQ));
console.log('\n  Total galaxies with valid m=2 data: ' + valid.length);
console.log('  From Program 9: ' + valid.filter(r => r.source === 'Program9').length);
console.log('  New in Program 10: ' + valid.filter(r => r.source !== 'Program9').length);

if (valid.length >= 4) {
  const dqs = valid.map(r => r.DQ);
  const logM2s = valid.map(r => Math.log10(r.m2Power));

  const rDQ_m2 = pearsonR(dqs, logM2s);
  const pVal = permTest(dqs, logM2s, 10000);

  console.log('\n  RESULTS:');
  console.log('    r(DQ, log m2) = ' + rDQ_m2.toFixed(4) + '  (N = ' + valid.length + ')');
  console.log('    p-value (permutation, 10000) = ' + pVal.toFixed(4));
  console.log('    Program 9 result was: r = 0.847 (N = 7)');
  console.log('    Delta r = ' + (rDQ_m2 - 0.847).toFixed(4));

  const rRaw_m2 = pearsonR(dqs, valid.map(r => r.m2Power));
  console.log('    r(DQ, m2 raw) = ' + rRaw_m2.toFixed(4));

  console.log('\n  Per-galaxy detail:');
  console.log('  ' + 'Galaxy'.padEnd(15) + 'DQ'.padEnd(8) + 'log(m2)'.padEnd(10) + 'Vflat'.padEnd(8) + 'Survey'.padEnd(10) + 'Source');
  console.log('  ' + '-'.repeat(65));
  for (const r of valid.sort((a, b) => b.DQ - a.DQ)) {
    console.log('  ' + r.name.padEnd(15) + r.DQ.toFixed(2).padEnd(8) + Math.log10(r.m2Power).toFixed(3).padEnd(10) + Math.round(r.Vflat).toString().padEnd(8) + r.survey.padEnd(10) + r.source);
  }

  console.log('\n\n  LOO cross-validation:');
  const looRs = [];
  for (let i = 0; i < valid.length; i++) {
    const xLoo = dqs.filter((_, j) => j !== i);
    const yLoo = logM2s.filter((_, j) => j !== i);
    looRs.push(pearsonR(xLoo, yLoo));
  }
  const allPositive = looRs.every(r => r > 0);
  const minLoo = Math.min(...looRs);
  const maxLoo = Math.max(...looRs);
  const meanLoo = looRs.reduce((a, b) => a + b, 0) / looRs.length;
  console.log('    LOO correlations: ' + looRs.length + '/' + looRs.length + ' positive = ' + allPositive);
  console.log('    Mean r(LOO) = ' + meanLoo.toFixed(4));
  console.log('    Range: [' + minLoo.toFixed(4) + ', ' + maxLoo.toFixed(4) + ']');

  for (let i = 0; i < valid.length; i++) {
    console.log('      Drop ' + valid[i].name.padEnd(12) + ' => r = ' + looRs[i].toFixed(4) + '  (delta = ' + (looRs[i] - rDQ_m2).toFixed(4) + ')');
  }

  console.log('\n\n  Survey independence test:');
  const thingsGals = valid.filter(r => r.survey === 'THINGS');
  const nonThingsGals = valid.filter(r => r.survey !== 'THINGS');
  console.log('    THINGS-only (N=' + thingsGals.length + '): r = ' + pearsonR(thingsGals.map(r => r.DQ), thingsGals.map(r => Math.log10(r.m2Power))).toFixed(4));
  if (nonThingsGals.length >= 3) {
    console.log('    Non-THINGS (N=' + nonThingsGals.length + '): r = ' + pearsonR(nonThingsGals.map(r => r.DQ), nonThingsGals.map(r => Math.log10(r.m2Power))).toFixed(4));
  } else {
    console.log('    Non-THINGS: insufficient data (N=' + nonThingsGals.length + ', need >= 3)');
  }

  const result = {
    program: 10,
    phase: '10.1',
    title: 'Multi-Survey 2D Carrier Verification',
    timestamp: new Date().toISOString(),
    N_total: valid.length,
    N_program9: valid.filter(r => r.source === 'Program9').length,
    N_new: valid.filter(r => r.source !== 'Program9').length,
    correlation: { r_DQ_logM2: rDQ_m2, pValue: pVal, r_program9: 0.847, delta_r: rDQ_m2 - 0.847 },
    loo: { mean: meanLoo, min: minLoo, max: maxLoo, allPositive, n: looRs.length },
    surveys: {
      THINGS: thingsGals.length,
      nonTHINGS: nonThingsGals.length,
      r_THINGS: pearsonR(thingsGals.map(r => r.DQ), thingsGals.map(r => Math.log10(r.m2Power))),
    },
    results: valid.sort((a, b) => b.DQ - a.DQ),
    newTargetsIdentified: priorityTargets.length,
    surveysNeeded: ['MaNGA', 'CALIFA', 'PHANGS'],
  };

  const outPath = path.join(__dirname, '..', 'public', 'program10-multi-survey.json');
  fs.writeFileSync(outPath, JSON.stringify(result, null, 2));
  console.log('\n  Results saved to: ' + outPath);

} else {
  console.log('\n  Insufficient valid data (N=' + valid.length + '). Need >= 4.');
}


console.log('\n\n' + '#'.repeat(72));
console.log('10.1.5 — DATA ACQUISITION PLAN');
console.log('#'.repeat(72));

console.log('\n  To reach N >= 20, the following data must be obtained:\n');

const needData = priorityTargets.filter(t => !thingsProcessed.has(normalize(t.name)));
const byPriority = needData.sort((a, b) => Math.abs(b.DQ) - Math.abs(a.DQ));

console.log('  PRIORITY A — Extreme DQ (|DQ| > 1.0), any survey:');
const prioA = byPriority.filter(t => Math.abs(t.DQ) > 1.0);
for (const t of prioA) {
  console.log('    ' + t.name.padEnd(15) + 'DQ=' + t.DQ.toFixed(2).padEnd(8) + Array.from(t.surveys).join(', '));
}
console.log('    Count: ' + prioA.length);

console.log('\n  PRIORITY B — Moderate DQ (0.5 < |DQ| < 1.0), with multi-survey overlap:');
const prioB = byPriority.filter(t => Math.abs(t.DQ) > 0.5 && Math.abs(t.DQ) <= 1.0 && t.surveys.size > 1);
for (const t of prioB) {
  console.log('    ' + t.name.padEnd(15) + 'DQ=' + t.DQ.toFixed(2).padEnd(8) + Array.from(t.surveys).join(', '));
}
console.log('    Count: ' + prioB.length);

console.log('\n  PRIORITY C — Any DQ, single-survey:');
const prioC = byPriority.filter(t => !prioA.includes(t) && !prioB.includes(t));
console.log('    Count: ' + prioC.length);

console.log('\n  Total new targets identified: ' + needData.length);
console.log('  Current N (processed): ' + valid.length);
console.log('  Potential N with all acquisitions: ' + (valid.length + needData.length));
console.log('  Target N for definitive confirmation: >= 20');

const canReach20 = (valid.length + needData.length) >= 20;
console.log('  Can reach N=20? ' + (canReach20 ? 'YES' : 'NOT WITH CURRENT SURVEY LISTS — need additional surveys'));


console.log('\n\n' + '='.repeat(72));
console.log('PROGRAM 10.1 SUMMARY');
console.log('='.repeat(72));
console.log('\n  Galaxies processed with m=2 decomposition: ' + valid.length);
console.log('  Surveys represented: ' + [...new Set(valid.map(r => r.survey))].join(', '));
console.log('  New targets identified for acquisition: ' + needData.length);
console.log('  r(DQ, log m2) with current data: ' + (valid.length >= 4 ? pearsonR(valid.map(r => r.DQ), valid.map(r => Math.log10(r.m2Power))).toFixed(4) : 'insufficient data'));
console.log('\n  Next step: Acquire velocity field data for Priority A targets');
console.log('  (MaNGA DAP maps, CALIFA V500 cubes, PHANGS-MUSE velocity maps)');
console.log('='.repeat(72));
