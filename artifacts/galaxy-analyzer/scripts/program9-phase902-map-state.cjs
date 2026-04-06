const fs = require('fs');
const path = require('path');

function parseFITS(filePath) {
  const buf = fs.readFileSync(filePath);
  const headerBlocks = [];
  let offset = 0;
  let headerDone = false;
  while (!headerDone && offset < buf.length) {
    const block = buf.slice(offset, offset + 2880).toString('ascii');
    offset += 2880;
    for (let i = 0; i < 36; i++) {
      const card = block.substring(i * 80, (i + 1) * 80);
      headerBlocks.push(card);
      if (card.startsWith('END')) { headerDone = true; break; }
    }
  }

  const header = {};
  for (const card of headerBlocks) {
    if (card.startsWith('END')) break;
    if (!card.includes('=')) continue;
    const key = card.substring(0, 8).trim();
    let valStr = card.substring(10, 80);
    const commentIdx = valStr.indexOf('/');
    if (commentIdx >= 0 && !valStr.trimStart().startsWith("'")) valStr = valStr.substring(0, commentIdx);
    valStr = valStr.trim().replace(/'/g, '').trim();
    const numVal = parseFloat(valStr);
    header[key] = isNaN(numVal) ? valStr : numVal;
  }

  const naxis = header['NAXIS'] || 0;
  const dims = [];
  for (let i = 1; i <= naxis; i++) dims.push(header['NAXIS' + i] || 0);
  const bitpix = header['BITPIX'] || -32;
  const blank = header['BLANK'];

  const totalPixels = dims.reduce((a, b) => a * b, 1);
  const bytesPerPixel = Math.abs(bitpix) / 8;
  const dataSize = totalPixels * bytesPerPixel;

  const data = new Float64Array(totalPixels);
  for (let i = 0; i < totalPixels; i++) {
    const pos = offset + i * bytesPerPixel;
    if (pos + bytesPerPixel > buf.length) { data[i] = NaN; continue; }
    if (bitpix === -32) {
      const b = Buffer.alloc(4);
      b[0] = buf[pos]; b[1] = buf[pos + 1]; b[2] = buf[pos + 2]; b[3] = buf[pos + 3];
      data[i] = b.readFloatBE(0);
    } else if (bitpix === -64) {
      const b = Buffer.alloc(8);
      for (let j = 0; j < 8; j++) b[j] = buf[pos + j];
      data[i] = b.readDoubleBE(0);
    } else if (bitpix === 16) {
      let val = (buf[pos] << 8) | buf[pos + 1];
      if (val >= 32768) val -= 65536;
      if (blank !== undefined && val === blank) data[i] = NaN;
      else data[i] = val;
    } else if (bitpix === 32) {
      let val = (buf[pos] << 24) | (buf[pos + 1] << 16) | (buf[pos + 2] << 8) | buf[pos + 3];
      if (blank !== undefined && val === blank) data[i] = NaN;
      else data[i] = val;
    } else {
      data[i] = NaN;
    }
    if (isNaN(data[i]) || !isFinite(data[i]) || Math.abs(data[i]) > 1e30) data[i] = NaN;
  }

  if (header['BSCALE'] && header['BSCALE'] !== 1) {
    const bs = header['BSCALE'], bz = header['BZERO'] || 0;
    for (let i = 0; i < data.length; i++) if (!isNaN(data[i])) data[i] = data[i] * bs + bz;
  } else if (header['BZERO'] && header['BZERO'] !== 0) {
    const bz = header['BZERO'];
    for (let i = 0; i < data.length; i++) if (!isNaN(data[i])) data[i] += bz;
  }

  return { header, dims, data };
}

function analyzeVelocityField(fits, galaxyParams) {
  const nx = fits.dims[0], ny = fits.dims[1];
  const nz = fits.dims.length > 2 ? fits.dims[2] : 1;
  const crpix1 = fits.header['CRPIX1'] || nx / 2;
  const crpix2 = fits.header['CRPIX2'] || ny / 2;
  const cdelt1 = fits.header['CDELT1'] || 1;
  const cdelt2 = fits.header['CDELT2'] || 1;
  const crval1 = fits.header['CRVAL1'] || 0;
  const crval2 = fits.header['CRVAL2'] || 0;

  const pixScale = Math.abs(cdelt1) * 3600;

  const vel = [];
  const coords = [];
  let validCount = 0;

  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const idx = j * nx + i;
      const v = fits.data[idx];
      if (isNaN(v)) continue;
      const x = (i + 1 - crpix1) * cdelt1 * 3600;
      const y = (j + 1 - crpix2) * cdelt2 * 3600;
      vel.push(v);
      coords.push({ x, y, i, j });
      validCount++;
    }
  }

  const vsys = vel.slice().sort((a, b) => a - b)[Math.floor(vel.length / 2)];
  const vrel = vel.map(v => v - vsys);

  const vMean = vel.reduce((a, b) => a + b, 0) / vel.length;
  const vStd = Math.sqrt(vel.reduce((s, v) => s + (v - vMean) ** 2, 0) / vel.length);

  let asymmetry180 = 0, asymmetryCount = 0;
  const cx = crpix1 - 1, cy = crpix2 - 1;

  for (let j = 0; j < ny; j++) {
    for (let i = 0; i < nx; i++) {
      const v1 = fits.data[j * nx + i];
      if (isNaN(v1)) continue;
      const oi = Math.round(2 * cx - i), oj = Math.round(2 * cy - j);
      if (oi < 0 || oi >= nx || oj < 0 || oj >= ny) continue;
      const v2 = fits.data[oj * nx + oi];
      if (isNaN(v2)) continue;
      const expected = -(v2 - vsys);
      const actual = v1 - vsys;
      if (Math.abs(actual) > 5) {
        asymmetry180 += ((actual - expected) / Math.max(Math.abs(actual), 1)) ** 2;
        asymmetryCount++;
      }
    }
  }
  asymmetry180 = asymmetryCount > 0 ? Math.sqrt(asymmetry180 / asymmetryCount) : NaN;

  const radialBins = 20;
  const maxR = Math.max(nx, ny) * Math.abs(cdelt1) * 3600 / 2;
  const binWidth = maxR / radialBins;
  const azBins = 8;
  const azimuthalPower = [];
  const radialAsymmetry = [];

  for (let rb = 0; rb < radialBins; rb++) {
    const rMin = rb * binWidth, rMax = (rb + 1) * binWidth;
    const azVals = Array.from({ length: azBins }, () => []);

    for (let k = 0; k < coords.length; k++) {
      const c = coords[k];
      const r = Math.sqrt(c.x ** 2 + c.y ** 2);
      if (r < rMin || r >= rMax) continue;
      const theta = Math.atan2(c.y, c.x);
      const azIdx = Math.floor(((theta + Math.PI) / (2 * Math.PI)) * azBins) % azBins;
      azVals[azIdx].push(vrel[k]);
    }

    const azMeans = azVals.map(arr => arr.length > 0 ? arr.reduce((a, b) => a + b, 0) / arr.length : NaN);
    const validAz = azMeans.filter(v => !isNaN(v));
    if (validAz.length < 4) {
      azimuthalPower.push({ r: (rMin + rMax) / 2, m0: NaN, m1: NaN, m2: NaN, totalPower: NaN });
      radialAsymmetry.push(NaN);
      continue;
    }

    const azMean = validAz.reduce((a, b) => a + b, 0) / validAz.length;
    let m0 = azMean;
    let cos1 = 0, sin1 = 0, cos2 = 0, sin2 = 0, count = 0;
    for (let a = 0; a < azBins; a++) {
      if (isNaN(azMeans[a])) continue;
      const angle = (a + 0.5) * 2 * Math.PI / azBins;
      const dv = azMeans[a] - azMean;
      cos1 += dv * Math.cos(angle);
      sin1 += dv * Math.sin(angle);
      cos2 += dv * Math.cos(2 * angle);
      sin2 += dv * Math.sin(2 * angle);
      count++;
    }
    const m1power = Math.sqrt(cos1 ** 2 + sin1 ** 2) / count;
    const m2power = Math.sqrt(cos2 ** 2 + sin2 ** 2) / count;

    azimuthalPower.push({ r: (rMin + rMax) / 2, m0: Math.abs(azMean), m1: m1power, m2: m2power, totalPower: m1power + m2power });

    const azStd = Math.sqrt(validAz.reduce((s, v) => s + (v - azMean) ** 2, 0) / validAz.length);
    radialAsymmetry.push(Math.abs(azMean) > 1 ? azStd / Math.abs(azMean) : NaN);
  }

  const validPower = azimuthalPower.filter(p => !isNaN(p.totalPower) && p.m0 > 0);
  const m1Total = validPower.length > 0 ? validPower.reduce((s, p) => s + p.m1, 0) / validPower.length : NaN;
  const m2Total = validPower.length > 0 ? validPower.reduce((s, p) => s + p.m2, 0) / validPower.length : NaN;
  const m1m0ratio = validPower.length > 0 ? validPower.reduce((s, p) => s + (p.m0 > 5 ? p.m1 / p.m0 : 0), 0) / validPower.filter(p => p.m0 > 5).length : NaN;
  const m2m0ratio = validPower.length > 0 ? validPower.reduce((s, p) => s + (p.m0 > 5 ? p.m2 / p.m0 : 0), 0) / validPower.filter(p => p.m0 > 5).length : NaN;

  let quadrantVels = [[], [], [], []];
  for (let k = 0; k < coords.length; k++) {
    const c = coords[k];
    const qx = c.x >= 0 ? 1 : 0;
    const qy = c.y >= 0 ? 1 : 0;
    const q = qy * 2 + qx;
    quadrantVels[q].push(vrel[k]);
  }
  const qMeans = quadrantVels.map(arr => arr.length > 0 ? arr.reduce((a, b) => a + b, 0) / arr.length : 0);
  const qStds = quadrantVels.map((arr, idx) => {
    if (arr.length < 2) return 0;
    const m = qMeans[idx];
    return Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / arr.length);
  });
  const quadrantBalance = Math.abs(qMeans[0] + qMeans[3]) + Math.abs(qMeans[1] + qMeans[2]);
  const maxQdiff = Math.max(...qMeans) - Math.min(...qMeans);

  let halfFluxAsym = 0;
  let upperSum = 0, lowerSum = 0, upperN = 0, lowerN = 0;
  for (let k = 0; k < coords.length; k++) {
    if (coords[k].y >= 0) { upperSum += Math.abs(vrel[k]); upperN++; }
    else { lowerSum += Math.abs(vrel[k]); lowerN++; }
  }
  if (upperN > 0 && lowerN > 0) {
    halfFluxAsym = Math.abs(upperSum / upperN - lowerSum / lowerN) / Math.max(upperSum / upperN, lowerSum / lowerN, 1);
  }

  let kinPA_residual = NaN;
  let maxVel = -Infinity, minVel = Infinity, maxAngle = 0, minAngle = 0;
  for (let k = 0; k < coords.length; k++) {
    const r = Math.sqrt(coords[k].x ** 2 + coords[k].y ** 2);
    if (r < 20 || r > maxR * 0.8) continue;
    if (vrel[k] > maxVel) { maxVel = vrel[k]; maxAngle = Math.atan2(coords[k].y, coords[k].x); }
    if (vrel[k] < minVel) { minVel = vrel[k]; minAngle = Math.atan2(coords[k].y, coords[k].x); }
  }
  const kinPA = (maxAngle + minAngle + Math.PI) / 2;

  const innerR = maxR * 0.3, outerR = maxR * 0.7;
  let innerPA_x = 0, innerPA_y = 0, outerPA_x = 0, outerPA_y = 0;
  for (let k = 0; k < coords.length; k++) {
    const r = Math.sqrt(coords[k].x ** 2 + coords[k].y ** 2);
    const w = vrel[k];
    if (r < innerR && r > 10) {
      innerPA_x += w * coords[k].x;
      innerPA_y += w * coords[k].y;
    } else if (r > outerR && r < maxR * 0.9) {
      outerPA_x += w * coords[k].x;
      outerPA_y += w * coords[k].y;
    }
  }
  const innerPA = Math.atan2(innerPA_y, innerPA_x);
  const outerPA = Math.atan2(outerPA_y, outerPA_x);
  const paTwist = Math.abs(outerPA - innerPA) * 180 / Math.PI;
  const paTwistNorm = paTwist > 180 ? 360 - paTwist : paTwist;

  let coherenceInner = 0, coherenceOuter = 0, nInner = 0, nOuter = 0;
  for (let k = 0; k < coords.length; k++) {
    const r = Math.sqrt(coords[k].x ** 2 + coords[k].y ** 2);
    const sign = Math.sign(vrel[k]);
    const expectedSign = Math.sign(coords[k].x * Math.cos(kinPA) + coords[k].y * Math.sin(kinPA));
    if (r < innerR && r > 5) {
      coherenceInner += (sign === expectedSign || sign === 0) ? 1 : 0;
      nInner++;
    } else if (r > outerR) {
      coherenceOuter += (sign === expectedSign || sign === 0) ? 1 : 0;
      nOuter++;
    }
  }
  coherenceInner = nInner > 0 ? coherenceInner / nInner : NaN;
  coherenceOuter = nOuter > 0 ? coherenceOuter / nOuter : NaN;

  return {
    validPixels: validCount,
    imageSize: [nx, ny],
    pixelScale_arcsec: pixScale,
    vsys: vsys,
    vRange: [vel.reduce((a, b) => Math.min(a, b), Infinity), vel.reduce((a, b) => Math.max(a, b), -Infinity)],
    vStd,
    asymmetry180,
    m1_mean: m1Total,
    m2_mean: m2Total,
    m1m0_ratio: m1m0ratio,
    m2m0_ratio: m2m0ratio,
    quadrantBalance,
    quadrantMeans: qMeans,
    halfFluxAsymmetry: halfFluxAsym,
    paTwist_deg: paTwistNorm,
    coherenceInner,
    coherenceOuter,
    coherenceRatio: coherenceOuter > 0 ? coherenceInner / coherenceOuter : NaN,
    azimuthalProfile: azimuthalPower,
  };
}

console.log('='.repeat(72));
console.log('PROGRAM 9 — PHASE 902: MAP-LEVEL STATE TEST');
console.log('2D velocity field analysis of matched High-H / Low-H pairs');
console.log('='.repeat(72));

const dataDir = path.join(__dirname, '..', 'data', 'things-2d');

const galaxies = [
  { name: 'NGC2841', DQ: 2.25, role: 'HIGH-H', dist: 14.1, Vflat: 285 },
  { name: 'NGC5055', DQ: -1.52, role: 'LOW-H', dist: 9.9, Vflat: 179 },
  { name: 'NGC3521', DQ: 0.52, role: 'MID-HIGH', dist: 7.7, Vflat: 214 },
  { name: 'NGC6946', DQ: -1.07, role: 'LOW-H', dist: 5.9, Vflat: 158 },
  { name: 'NGC7331', DQ: -0.14, role: 'NEUTRAL', dist: 14.7, Vflat: 239 },
  { name: 'NGC2403', DQ: 0.14, role: 'NEUTRAL', dist: 3.2, Vflat: 131 },
];

const results = [];

for (const gal of galaxies) {
  const mom1File = path.join(dataDir, gal.name + '_MOM1.FITS');
  if (!fs.existsSync(mom1File)) {
    console.log('\n  SKIPPED ' + gal.name + ': no MOM1 file');
    continue;
  }

  console.log('\n' + '-'.repeat(72));
  console.log('Analyzing: ' + gal.name + ' (' + gal.role + ', DQ=' + gal.DQ.toFixed(2) + ')');
  console.log('-'.repeat(72));

  const fits = parseFITS(mom1File);
  console.log('  Image dimensions: ' + fits.dims.join(' x '));
  console.log('  BITPIX: ' + fits.header['BITPIX']);
  console.log('  CRPIX: [' + fits.header['CRPIX1'] + ', ' + fits.header['CRPIX2'] + ']');
  console.log('  CDELT: [' + fits.header['CDELT1'] + ', ' + fits.header['CDELT2'] + ']');

  const analysis = analyzeVelocityField(fits, gal);

  console.log('\n  Valid pixels: ' + analysis.validPixels);
  console.log('  Pixel scale: ' + analysis.pixelScale_arcsec.toFixed(2) + ' arcsec');
  console.log('  Vsys: ' + analysis.vsys.toFixed(1) + ' km/s');
  console.log('  V range: [' + analysis.vRange[0].toFixed(1) + ', ' + analysis.vRange[1].toFixed(1) + '] km/s');
  console.log('  V std: ' + analysis.vStd.toFixed(1) + ' km/s');
  console.log('\n  --- 2D ASYMMETRY METRICS ---');
  console.log('  180° rotational asymmetry: ' + analysis.asymmetry180.toFixed(4));
  console.log('  m=1 azimuthal power (mean): ' + analysis.m1_mean.toFixed(3) + ' km/s');
  console.log('  m=2 azimuthal power (mean): ' + analysis.m2_mean.toFixed(3) + ' km/s');
  console.log('  m1/m0 ratio: ' + analysis.m1m0_ratio.toFixed(4));
  console.log('  m2/m0 ratio: ' + analysis.m2m0_ratio.toFixed(4));
  console.log('  Quadrant balance: ' + analysis.quadrantBalance.toFixed(2) + ' km/s');
  console.log('  Half-flux asymmetry: ' + analysis.halfFluxAsymmetry.toFixed(4));
  console.log('  PA twist (inner-outer): ' + analysis.paTwist_deg.toFixed(1) + '°');
  console.log('  Coherence inner: ' + analysis.coherenceInner.toFixed(3));
  console.log('  Coherence outer: ' + analysis.coherenceOuter.toFixed(3));
  console.log('  Coherence ratio (inner/outer): ' + analysis.coherenceRatio.toFixed(3));

  results.push({
    name: gal.name,
    DQ: gal.DQ,
    role: gal.role,
    Vflat: gal.Vflat,
    dist: gal.dist,
    validPixels: analysis.validPixels,
    asymmetry180: analysis.asymmetry180,
    m1_mean: analysis.m1_mean,
    m2_mean: analysis.m2_mean,
    m1m0_ratio: analysis.m1m0_ratio,
    m2m0_ratio: analysis.m2m0_ratio,
    quadrantBalance: analysis.quadrantBalance,
    halfFluxAsymmetry: analysis.halfFluxAsymmetry,
    paTwist_deg: analysis.paTwist_deg,
    coherenceInner: analysis.coherenceInner,
    coherenceOuter: analysis.coherenceOuter,
    coherenceRatio: analysis.coherenceRatio,
  });
}

console.log('\n\n' + '#'.repeat(72));
console.log('902.1 — COMPARATIVE ANALYSIS');
console.log('#'.repeat(72));

const metrics = ['asymmetry180', 'm1_mean', 'm2_mean', 'm1m0_ratio', 'm2m0_ratio', 'quadrantBalance', 'halfFluxAsymmetry', 'paTwist_deg', 'coherenceInner', 'coherenceOuter', 'coherenceRatio'];

console.log('\n  Summary table:\n');
console.log('  ' + 'Galaxy'.padEnd(12) + 'Role'.padEnd(12) + 'DQ'.padEnd(8) + 'asym180'.padEnd(10) + 'm1/m0'.padEnd(10) + 'm2/m0'.padEnd(10) + 'quadBal'.padEnd(10) + 'halfAsym'.padEnd(10) + 'paTwist'.padEnd(10) + 'cohRatio'.padEnd(10));
console.log('  ' + '-'.repeat(100));

for (const r of results) {
  console.log('  ' + r.name.padEnd(12) + r.role.padEnd(12) + r.DQ.toFixed(2).padEnd(8) + r.asymmetry180.toFixed(4).padEnd(10) + r.m1m0_ratio.toFixed(4).padEnd(10) + r.m2m0_ratio.toFixed(4).padEnd(10) + r.quadrantBalance.toFixed(1).padEnd(10) + r.halfFluxAsymmetry.toFixed(4).padEnd(10) + r.paTwist_deg.toFixed(1).padEnd(10) + r.coherenceRatio.toFixed(3).padEnd(10));
}


console.log('\n\n' + '#'.repeat(72));
console.log('902.2 — HIGH-H vs LOW-H METRIC COMPARISON');
console.log('#'.repeat(72));

const highHResults = results.filter(r => r.DQ > 0.5);
const lowHResults = results.filter(r => r.DQ < -0.5);

console.log('\n  High-H galaxies: ' + highHResults.map(r => r.name).join(', '));
console.log('  Low-H galaxies: ' + lowHResults.map(r => r.name).join(', '));
console.log('');

for (const m of metrics) {
  const highVals = highHResults.map(r => r[m]).filter(v => !isNaN(v));
  const lowVals = lowHResults.map(r => r[m]).filter(v => !isNaN(v));
  if (highVals.length === 0 || lowVals.length === 0) continue;
  const highMean = highVals.reduce((a, b) => a + b, 0) / highVals.length;
  const lowMean = lowVals.reduce((a, b) => a + b, 0) / lowVals.length;
  const ratio = lowMean !== 0 ? highMean / lowMean : NaN;
  const diff = highMean - lowMean;
  const allVals = highVals.concat(lowVals);
  const globalMean = allVals.reduce((a, b) => a + b, 0) / allVals.length;
  const globalStd = Math.sqrt(allVals.reduce((s, v) => s + (v - globalMean) ** 2, 0) / allVals.length);
  const effectSize = globalStd > 0 ? diff / globalStd : 0;
  console.log('  ' + m.padEnd(20) + 'High=' + highMean.toFixed(4).padEnd(10) + 'Low=' + lowMean.toFixed(4).padEnd(10) + 'Δ=' + diff.toFixed(4).padEnd(10) + 'ratio=' + ratio.toFixed(2).padEnd(8) + 'effect=' + effectSize.toFixed(2));
}


console.log('\n\n' + '#'.repeat(72));
console.log('902.3 — GOLD PAIR DEEP COMPARISON: NGC 2841 vs NGC 5055');
console.log('#'.repeat(72));

const ngc2841 = results.find(r => r.name === 'NGC2841');
const ngc5055 = results.find(r => r.name === 'NGC5055');

if (ngc2841 && ngc5055) {
  console.log('\n  NGC 2841 (HIGH-H, DQ=+2.25) vs NGC 5055 (LOW-H, DQ=-1.52)');
  console.log('  DQ difference: 3.78 sigma\n');

  for (const m of metrics) {
    const v1 = ngc2841[m], v2 = ngc5055[m];
    if (isNaN(v1) || isNaN(v2)) continue;
    const ratio = v2 !== 0 ? v1 / v2 : NaN;
    const winner = v1 > v2 ? 'NGC2841 higher' : v1 < v2 ? 'NGC5055 higher' : 'equal';
    console.log('  ' + m.padEnd(22) + 'NGC2841=' + v1.toFixed(4).padEnd(12) + 'NGC5055=' + v2.toFixed(4).padEnd(12) + 'ratio=' + ratio.toFixed(3).padEnd(8) + winner);
  }
}


console.log('\n\n' + '#'.repeat(72));
console.log('902.4 — DQ CORRELATION WITH 2D METRICS');
console.log('#'.repeat(72));

function pearsonR(x, y) {
  const n = x.length; if (n < 3) return NaN;
  const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; sxy += dx * dy; sxx += dx * dx; syy += dy * dy; }
  return (sxx > 0 && syy > 0) ? sxy / Math.sqrt(sxx * syy) : 0;
}

const dqVals = results.map(r => r.DQ);
console.log('\n  Correlations between DQ and 2D metrics (N=' + results.length + '):');
for (const m of metrics) {
  const mVals = results.map(r => r[m]);
  const validPairs = [];
  for (let i = 0; i < results.length; i++) {
    if (!isNaN(mVals[i]) && !isNaN(dqVals[i])) validPairs.push({ dq: dqVals[i], mv: mVals[i] });
  }
  if (validPairs.length < 3) continue;
  const r = pearsonR(validPairs.map(p => p.dq), validPairs.map(p => p.mv));
  const sig = Math.abs(r) > 0.7 ? ' *** STRONG' : Math.abs(r) > 0.4 ? ' ** MODERATE' : Math.abs(r) > 0.2 ? ' * WEAK' : '';
  console.log('  r(DQ, ' + m.padEnd(20) + ') = ' + r.toFixed(3) + ' (N=' + validPairs.length + ')' + sig);
}


console.log('\n\n' + '#'.repeat(72));
console.log('902.5 — PHASE 902 VERDICT');
console.log('#'.repeat(72));

const strongMetrics = [];
const weakMetrics = [];
for (const m of metrics) {
  const mVals = results.map(r => r[m]);
  const validPairs = [];
  for (let i = 0; i < results.length; i++) {
    if (!isNaN(mVals[i]) && !isNaN(dqVals[i])) validPairs.push({ dq: dqVals[i], mv: mVals[i] });
  }
  if (validPairs.length < 3) continue;
  const r = pearsonR(validPairs.map(p => p.dq), validPairs.map(p => p.mv));
  if (Math.abs(r) > 0.5) strongMetrics.push({ metric: m, r });
  else if (Math.abs(r) > 0.3) weakMetrics.push({ metric: m, r });
}

console.log('\n  Strong DQ correlations (|r| > 0.5): ' + strongMetrics.length);
for (const s of strongMetrics) console.log('    ' + s.metric + ': r = ' + s.r.toFixed(3));

console.log('  Moderate DQ correlations (0.3 < |r| < 0.5): ' + weakMetrics.length);
for (const w of weakMetrics) console.log('    ' + w.metric + ': r = ' + w.r.toFixed(3));

if (strongMetrics.length >= 2) {
  console.log('\n  VERDICT: STRONG 2D SIGNAL DETECTED');
  console.log('  Multiple 2D metrics show significant correlation with DQ.');
  console.log('  H has a detectable 2D kinematic signature.');
} else if (strongMetrics.length >= 1 || weakMetrics.length >= 3) {
  console.log('\n  VERDICT: PARTIAL 2D SIGNAL');
  console.log('  Some 2D metrics correlate with DQ but not all.');
  console.log('  H may have a partial 2D component.');
} else {
  console.log('\n  VERDICT: NO CLEAR 2D SIGNAL (from THINGS velocity fields alone)');
  console.log('  DQ does not strongly predict 2D kinematic asymmetry metrics.');
  console.log('  This could mean: H is NOT angular, OR sample too small, OR metrics need refinement.');
}


const outPath = path.join(__dirname, '..', 'public', 'program9-phase902.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: 9,
  phase: 902,
  title: 'Map-Level State Test — THINGS Velocity Fields',
  timestamp: new Date().toISOString(),
  galaxiesAnalyzed: results.length,
  results,
  strongCorrelations: strongMetrics,
  moderateCorrelations: weakMetrics,
}, null, 2));
console.log('\nSaved: ' + outPath);
