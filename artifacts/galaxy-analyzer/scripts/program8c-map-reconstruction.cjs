const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

const G = 4.3009e-6;
const Nradial = 20;

function normalize(n) {
  return n.replace(/\s+/g, '').replace(/^NGC0*/, 'NGC').replace(/^DDO0*/, 'DDO').replace(/^IC0*/, 'IC').replace(/^UGC0*/, 'UGC').replace(/^PGC0*/, 'PGC').toUpperCase();
}

function pearsonR(x, y) {
  const n = x.length; if (n < 4) return NaN;
  const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; num += dx * dy; dx2 += dx * dx; dy2 += dy * dy; }
  return dx2 > 0 && dy2 > 0 ? num / Math.sqrt(dx2 * dy2) : 0;
}

function multiR2(X, y) {
  const n = y.length, nv = X[0].length;
  const my = y.reduce((a, b) => a + b, 0) / n;
  const mx = Array(nv).fill(0);
  for (let j = 0; j < nv; j++) { for (let i = 0; i < n; i++) mx[j] += X[i][j]; mx[j] /= n; }
  const XTX = Array.from({ length: nv }, () => Array(nv).fill(0)), XTy = Array(nv).fill(0);
  for (let i = 0; i < n; i++) { for (let j = 0; j < nv; j++) { XTy[j] += (X[i][j] - mx[j]) * (y[i] - my); for (let k = 0; k < nv; k++) XTX[j][k] += (X[i][j] - mx[j]) * (X[i][k] - mx[k]); } }
  const aug = XTX.map((row, i) => [...row, XTy[i]]);
  for (let col = 0; col < nv; col++) { let maxRow = col; for (let row = col + 1; row < nv; row++) if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row; [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]]; if (Math.abs(aug[col][col]) < 1e-12) continue; for (let row = col + 1; row < nv; row++) { const f = aug[row][col] / aug[col][col]; for (let j = col; j <= nv; j++) aug[row][j] -= f * aug[col][j]; } }
  const beta = Array(nv).fill(0);
  for (let i = nv - 1; i >= 0; i--) { beta[i] = aug[i][nv]; for (let j = i + 1; j < nv; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= aug[i][i] || 1; }
  const resid = [];
  for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < nv; j++) pred += beta[j] * (X[i][j] - mx[j]); resid.push(y[i] - pred); }
  return { residuals: resid, beta };
}

function partialR(x, y, controls) {
  const n = x.length;
  if (n < controls[0].length + 4) return NaN;
  const residualise = (v, C) => {
    const nv2 = C[0].length;
    const my2 = v.reduce((a, b) => a + b, 0) / n;
    const mc2 = Array(nv2).fill(0);
    for (let j = 0; j < nv2; j++) { for (let i = 0; i < n; i++) mc2[j] += C[i][j]; mc2[j] /= n; }
    const XTX2 = Array.from({ length: nv2 }, () => Array(nv2).fill(0));
    const XTy2 = Array(nv2).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < nv2; j++) {
        XTy2[j] += (C[i][j] - mc2[j]) * (v[i] - my2);
        for (let k = 0; k < nv2; k++) XTX2[j][k] += (C[i][j] - mc2[j]) * (C[i][k] - mc2[k]);
      }
    }
    const aug2 = XTX2.map((row, i) => [...row, XTy2[i]]);
    for (let col = 0; col < nv2; col++) {
      let maxRow = col;
      for (let row = col + 1; row < nv2; row++) if (Math.abs(aug2[row][col]) > Math.abs(aug2[maxRow][col])) maxRow = row;
      [aug2[col], aug2[maxRow]] = [aug2[maxRow], aug2[col]];
      if (Math.abs(aug2[col][col]) < 1e-12) continue;
      for (let row = col + 1; row < nv2; row++) { const f = aug2[row][col] / aug2[col][col]; for (let j = col; j <= nv2; j++) aug2[row][j] -= f * aug2[col][j]; }
    }
    const beta2 = Array(nv2).fill(0);
    for (let i = nv2 - 1; i >= 0; i--) { beta2[i] = aug2[i][nv2]; for (let j = i + 1; j < nv2; j++) beta2[i] -= aug2[i][j] * beta2[j]; beta2[i] /= aug2[i][i] || 1; }
    const r = [];
    for (let i = 0; i < n; i++) { let pred = my2; for (let j = 0; j < nv2; j++) pred += beta2[j] * (C[i][j] - mc2[j]); r.push(v[i] - pred); }
    return r;
  };
  return pearsonR(residualise(x, controls), residualise(y, controls));
}

const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[g.name] = g; resultsMap[normalize(g.name)] = g; });

const tsContent = fs.readFileSync(path.join(__dirname, '..', 'src', 'data', 'sparc-datasets.ts'), 'utf8');
const rcMapRaw = {};
const re = /"([^"]+)":\s*\{[^[]*data:\s*\[([\s\S]*?)\]/g;
let m2;
while ((m2 = re.exec(tsContent)) !== null) {
  const pts = [];
  const ptRe = /r:\s*([\d.]+)\s*,\s*v:\s*([\d.]+)/g;
  let pm;
  while ((pm = ptRe.exec(m2[2])) !== null) pts.push({ r: parseFloat(pm[1]), v: parseFloat(pm[2]) });
  if (pts.length >= 3) rcMapRaw[m2[1]] = pts;
}
const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) rcMap[normalize(name)] = pts;


function interpolateRC(rc, Rmax, Npts) {
  const profile = [];
  for (let i = 0; i < Npts; i++) {
    const rTarget = (i + 0.5) / Npts * Rmax;
    let lo = 0, hi = rc.length - 1;
    for (let j = 0; j < rc.length; j++) {
      if (rc[j].r <= rTarget) lo = j;
      if (rc[j].r >= rTarget && hi === rc.length - 1) hi = j;
    }
    if (lo === hi) { profile.push(rc[lo].v); continue; }
    const frac = (rTarget - rc[lo].r) / (rc[hi].r - rc[lo].r || 1);
    profile.push(rc[lo].v + frac * (rc[hi].v - rc[lo].v));
  }
  return profile;
}

function buildFullState(rc, Vflat, Rdisk, Mbar, Rmax) {
  const Npts = Nradial;
  const vProfile = interpolateRC(rc, Rmax, Npts);
  const vNorm = vProfile.map(v => v / Vflat);

  const haloFracProfile = [];
  for (let i = 0; i < Npts; i++) {
    const rMid = (i + 0.5) / Npts * Rmax;
    const encMass = Mbar * Math.min(rMid / Rmax, 1);
    const V_bar = Math.sqrt(G * encMass / Math.max(rMid, 0.01));
    const V_obs = vProfile[i];
    const haloFrac = Math.max(0, 1 - (V_bar / Math.max(V_obs, 1)) ** 2);
    haloFracProfile.push(haloFrac);
  }

  const residProfile = vNorm.map(v => v - 1);

  const gradProfile = [];
  for (let i = 1; i < Npts; i++) {
    gradProfile.push(vNorm[i] - vNorm[i - 1]);
  }

  const curvProfile = [];
  for (let i = 1; i < Npts - 1; i++) {
    curvProfile.push(vNorm[i + 1] - 2 * vNorm[i] + vNorm[i - 1]);
  }

  const haloGradProfile = [];
  for (let i = 1; i < Npts; i++) {
    haloGradProfile.push(haloFracProfile[i] - haloFracProfile[i - 1]);
  }

  const cumHaloProfile = [];
  let cumSum = 0;
  for (let i = 0; i < Npts; i++) {
    cumSum += haloFracProfile[i];
    cumHaloProfile.push(cumSum / (i + 1));
  }

  const fourierCoeffs = [];
  for (let k = 1; k <= 5; k++) {
    let cosSum = 0, sinSum = 0;
    for (let i = 0; i < Npts; i++) {
      const theta = 2 * Math.PI * k * i / Npts;
      cosSum += residProfile[i] * Math.cos(theta);
      sinSum += residProfile[i] * Math.sin(theta);
    }
    fourierCoeffs.push(Math.sqrt(cosSum ** 2 + sinSum ** 2) / Npts);
  }

  const turnoverBin = vNorm.indexOf(Math.max(...vNorm));
  const turnoverFrac = turnoverBin / Npts;

  const innerSlope = Npts >= 4 ? (vNorm[3] - vNorm[0]) / 3 : 0;
  const outerSlope = Npts >= 4 ? (vNorm[Npts - 1] - vNorm[Npts - 4]) / 3 : 0;

  let asymmetry = 0;
  const half = Math.floor(Npts / 2);
  for (let i = 0; i < half; i++) {
    const j = Npts - 1 - i;
    asymmetry += (haloFracProfile[i] - haloFracProfile[j]) ** 2;
  }
  asymmetry = Math.sqrt(asymmetry / half);

  let roughness = 0;
  for (let i = 1; i < Npts; i++) {
    roughness += (vNorm[i] - vNorm[i - 1]) ** 2;
  }
  roughness = Math.sqrt(roughness / (Npts - 1));

  const innerHaloMean = haloFracProfile.slice(0, Math.floor(Npts / 3)).reduce((a, b) => a + b, 0) / Math.floor(Npts / 3);
  const outerHaloMean = haloFracProfile.slice(Math.floor(2 * Npts / 3)).reduce((a, b) => a + b, 0) / (Npts - Math.floor(2 * Npts / 3));
  const haloContrast = outerHaloMean - innerHaloMean;

  return {
    vNorm, haloFracProfile, residProfile, gradProfile, curvProfile,
    haloGradProfile, cumHaloProfile, fourierCoeffs,
    turnoverFrac, innerSlope, outerSlope, asymmetry, roughness, haloContrast,
    fullVector: [
      ...vNorm, ...haloFracProfile, ...residProfile,
      ...fourierCoeffs,
      turnoverFrac, innerSlope, outerSlope, asymmetry, roughness, haloContrast
    ],
  };
}

function vectorDist(a, b) {
  let sum = 0;
  const len = Math.min(a.length, b.length);
  for (let i = 0; i < len; i++) sum += (a[i] - b[i]) ** 2;
  return Math.sqrt(sum / len);
}

function cosineSim(a, b) {
  let dot = 0, na = 0, nb = 0;
  const len = Math.min(a.length, b.length);
  for (let i = 0; i < len; i++) { dot += a[i] * b[i]; na += a[i] * a[i]; nb += b[i] * b[i]; }
  return na > 0 && nb > 0 ? dot / (Math.sqrt(na) * Math.sqrt(nb)) : 0;
}


const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sr) continue;
  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 10) continue;
  const Vflat = sp.Vflat, Rdisk = sp.Rdisk;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  const Rmax = rc[rc.length - 1].r;
  const logMbar = Math.log10(Math.max(Mbar, 1));
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(Rdisk, 0.01));
  const logSBdisk = Math.log10(Math.max(sp.SBdisk, 0.01));
  const logK = sr.models.dark_halo_linear.k > 0 ? Math.log10(sr.models.dark_halo_linear.k) : -5;
  const hR = sr.models.newtonian.mse > 0 ? Math.log10(Math.max(sr.models.newtonian.mse / Math.max(sr.models.dark_halo_linear.mse, 0.001), 0.01)) : 0;
  const V_Newt_Rmax = Math.sqrt(G * Mbar / Rmax);
  const dmFrac = Math.max(0, 1 - (V_Newt_Rmax / Vflat) ** 2);

  const state = buildFullState(rc, Vflat, Rdisk, Mbar, Rmax);
  if (!state) continue;

  gals.push({
    name: g.name, Vflat, Rdisk, Rmax, Mbar, rc,
    logVflat: Math.log10(Vflat), logL36, logRdisk, logMbar,
    logMHI: g.logMHI, morphT: sp.T, logA0: g.logA0,
    logSBdisk, envCode: g.envCode,
    logK, dmFrac, hR,
    state,
  });
}

const N = gals.length;

const struct4 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
const struct6 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
const vfM = multiR2(struct4, gals.map(g => g.logVflat));
const a0M = multiR2(struct6, gals.map(g => g.logA0));
for (let i = 0; i < N; i++) { gals[i].VfR = vfM.residuals[i]; gals[i].a0R = a0M.residuals[i]; }
const sdVf = Math.sqrt(gals.reduce((a, g) => a + g.VfR ** 2, 0) / N);
const sdA0 = Math.sqrt(gals.reduce((a, g) => a + g.a0R ** 2, 0) / N);
for (let i = 0; i < N; i++) {
  gals[i].VfR_z = gals[i].VfR / sdVf;
  gals[i].a0R_z = gals[i].a0R / sdA0;
  gals[i].L_sum = gals[i].VfR_z + gals[i].a0R_z;
}
const ctrlV = gals.map(g => [g.logK, g.dmFrac, g.envCode]);
const Lresid = multiR2(ctrlV, gals.map(g => g.L_sum));
for (let i = 0; i < N; i++) gals[i].dq = Lresid.residuals[i];
gals.sort((a, b) => b.dq - a.dq);

const r_VfR_a0R = pearsonR(gals.map(g => g.VfR), gals.map(g => g.a0R));
const vfArr = gals.map(g => g.VfR);
const a0Arr = gals.map(g => g.a0R);
const dqArr = gals.map(g => g.dq);


console.log('='.repeat(70));
console.log('PROGRAM 8C: MATCHED 2D MAP RECONSTRUCTION');
console.log('Can the FULL rotation-curve state — not scalar summaries —');
console.log('discriminate high-H from low-H galaxies?');
console.log('='.repeat(70));
console.log('\nN = ' + N + ' galaxies with >= 10 RC points');
console.log('Full state vector: ' + gals[0].state.fullVector.length + ' dimensions');
console.log('  - ' + Nradial + '-bin normalised velocity profile');
console.log('  - ' + Nradial + '-bin halo fraction profile');
console.log('  - ' + Nradial + '-bin residual profile');
console.log('  - 5 Fourier coefficients');
console.log('  - 6 shape scalars (turnover, slopes, asymmetry, roughness, contrast)');
console.log('r(VfR, a0R) = ' + r_VfR_a0R.toFixed(3));


console.log('\n\n' + '#'.repeat(70));
console.log('8C-1: MATCHED-PAIR FULL-STATE COMPARISON');
console.log('#'.repeat(70));

const targetNames = ['NGC2841', 'NGC3741', 'ESO563-G021'];
function matchScore(t, c) {
  const dVf = Math.abs(Math.log10(t.Vflat) - Math.log10(c.Vflat)) / 0.3;
  const dMb = Math.abs(t.logMbar - c.logMbar) / 0.5;
  const dRd = Math.abs(Math.log10(t.Rdisk) - Math.log10(c.Rdisk)) / 0.3;
  const dT = Math.abs(t.morphT - c.morphT) / 3;
  return Math.sqrt(dVf ** 2 + dMb ** 2 + dRd ** 2 + dT ** 2);
}

const targets = targetNames.map(n => gals.find(g => g.name === n)).filter(Boolean);
const matchedControls = [];
const used = new Set(targetNames);
for (const target of targets) {
  const lowDQ = gals.filter(g => g.dq < 0 && !used.has(g.name));
  const scored = lowDQ.map(g => ({ g, score: matchScore(target, g) })).sort((a, b) => a.score - b.score);
  matchedControls.push(scored[0].g);
  used.add(scored[0].g.name);
}

for (let p = 0; p < targets.length; p++) {
  const t = targets[p], c = matchedControls[p];
  console.log('\n  PAIR ' + (p + 1) + ': ' + t.name + ' (DQ=' + t.dq.toFixed(2) + ') vs ' + c.name + ' (DQ=' + c.dq.toFixed(2) + ')');
  console.log('  ' + '-'.repeat(80));

  const eucDist = vectorDist(t.state.fullVector, c.state.fullVector);
  const cosSim = cosineSim(t.state.fullVector, c.state.fullVector);
  const vProfileDist = vectorDist(t.state.vNorm, c.state.vNorm);
  const haloProfileDist = vectorDist(t.state.haloFracProfile, c.state.haloFracProfile);
  const fourierDist = vectorDist(t.state.fourierCoeffs, c.state.fourierCoeffs);

  console.log('  Full-state distance:  ' + eucDist.toFixed(4));
  console.log('  Cosine similarity:    ' + cosSim.toFixed(4));
  console.log('  V-profile distance:   ' + vProfileDist.toFixed(4));
  console.log('  Halo-frac distance:   ' + haloProfileDist.toFixed(4));
  console.log('  Fourier distance:     ' + fourierDist.toFixed(4));

  console.log('\n  NORMALISED VELOCITY PROFILE (20-bin):');
  let line1 = '  T: ';
  let line2 = '  C: ';
  let line3 = '  d: ';
  for (let b = 0; b < Nradial; b++) {
    line1 += t.state.vNorm[b].toFixed(2).padEnd(6);
    line2 += c.state.vNorm[b].toFixed(2).padEnd(6);
    line3 += ((t.state.vNorm[b] - c.state.vNorm[b]) >= 0 ? '+' : '') + (t.state.vNorm[b] - c.state.vNorm[b]).toFixed(2) + ' ';
    if (b === 9) { console.log(line1); console.log(line2); console.log(line3); line1 = '  T: '; line2 = '  C: '; line3 = '  d: '; }
  }
  console.log(line1); console.log(line2); console.log(line3);

  console.log('\n  HALO FRACTION PROFILE:');
  let h1 = '  T: ', h2 = '  C: ', h3 = '  d: ';
  for (let b = 0; b < Nradial; b++) {
    h1 += t.state.haloFracProfile[b].toFixed(2).padEnd(6);
    h2 += c.state.haloFracProfile[b].toFixed(2).padEnd(6);
    h3 += ((t.state.haloFracProfile[b] - c.state.haloFracProfile[b]) >= 0 ? '+' : '') + (t.state.haloFracProfile[b] - c.state.haloFracProfile[b]).toFixed(2) + ' ';
    if (b === 9) { console.log(h1); console.log(h2); console.log(h3); h1 = '  T: '; h2 = '  C: '; h3 = '  d: '; }
  }
  console.log(h1); console.log(h2); console.log(h3);

  console.log('\n  SHAPE SCALARS:');
  const scalars = ['turnoverFrac', 'innerSlope', 'outerSlope', 'asymmetry', 'roughness', 'haloContrast'];
  for (const s of scalars) {
    console.log('  ' + s.padEnd(20) + 'T=' + t.state[s].toFixed(4).padEnd(10) + 'C=' + c.state[s].toFixed(4).padEnd(10) + 'd=' + (t.state[s] - c.state[s]).toFixed(4));
  }

  console.log('\n  FOURIER SPECTRUM:');
  for (let k = 0; k < 5; k++) {
    console.log('  k=' + (k + 1) + '  T=' + t.state.fourierCoeffs[k].toFixed(4).padEnd(10) + 'C=' + c.state.fourierCoeffs[k].toFixed(4).padEnd(10) + 'd=' + (t.state.fourierCoeffs[k] - c.state.fourierCoeffs[k]).toFixed(4));
  }
}


console.log('\n\n' + '#'.repeat(70));
console.log('8C-2: QUINTILE PROFILE SEPARATION');
console.log('#'.repeat(70));

const quintileSize = Math.floor(N / 5);
const Q1 = gals.slice(0, quintileSize);
const Q5 = gals.slice(N - quintileSize);

const q1MeanV = Array(Nradial).fill(0);
const q5MeanV = Array(Nradial).fill(0);
const q1MeanH = Array(Nradial).fill(0);
const q5MeanH = Array(Nradial).fill(0);
for (const g of Q1) for (let b = 0; b < Nradial; b++) { q1MeanV[b] += g.state.vNorm[b] / Q1.length; q1MeanH[b] += g.state.haloFracProfile[b] / Q1.length; }
for (const g of Q5) for (let b = 0; b < Nradial; b++) { q5MeanV[b] += g.state.vNorm[b] / Q5.length; q5MeanH[b] += g.state.haloFracProfile[b] / Q5.length; }

console.log('\n  MEAN NORMALISED VELOCITY (Q1 high-H vs Q5 low-H):');
console.log('  ' + 'Bin'.padEnd(6) + 'Q1'.padEnd(10) + 'Q5'.padEnd(10) + 'Diff'.padEnd(10) + '| Halo Q1'.padEnd(12) + 'Halo Q5'.padEnd(12) + 'Halo diff');
console.log('  ' + '-'.repeat(75));
for (let b = 0; b < Nradial; b++) {
  const vd = q1MeanV[b] - q5MeanV[b];
  const hd = q1MeanH[b] - q5MeanH[b];
  console.log('  ' + (b + 1 + '').padEnd(6) + q1MeanV[b].toFixed(3).padEnd(10) + q5MeanV[b].toFixed(3).padEnd(10) + ((vd >= 0 ? '+' : '') + vd.toFixed(3)).padEnd(10) + '| ' + q1MeanH[b].toFixed(3).padEnd(12) + q5MeanH[b].toFixed(3).padEnd(12) + (hd >= 0 ? '+' : '') + hd.toFixed(3));
}

const vProfileSep = vectorDist(q1MeanV, q5MeanV);
const hProfileSep = vectorDist(q1MeanH, q5MeanH);
console.log('\n  V-profile quintile separation: ' + vProfileSep.toFixed(4));
console.log('  Halo-profile quintile separation: ' + hProfileSep.toFixed(4));


console.log('\n\n' + '#'.repeat(70));
console.log('8C-3: LEAVE-ONE-OUT NEAREST-NEIGHBOUR CLASSIFICATION');
console.log('Can full-state distance assign galaxies to correct H-class?');
console.log('#'.repeat(70));

const medianDQ = gals.slice().sort((a, b) => a.dq - b.dq)[Math.floor(N / 2)].dq;
for (let i = 0; i < N; i++) gals[i].highH = gals[i].dq > medianDQ;

let looCorrect = 0;
let looScalarCorrect = 0;
for (let i = 0; i < N; i++) {
  let nearestDist = Infinity, nearestClass = false;
  let nearestScalarDist = Infinity, nearestScalarClass = false;
  for (let j = 0; j < N; j++) {
    if (i === j) continue;
    const dist = vectorDist(gals[i].state.fullVector, gals[j].state.fullVector);
    if (dist < nearestDist) { nearestDist = dist; nearestClass = gals[j].highH; }
    const scalarDist = Math.abs(gals[i].hR - gals[j].hR);
    if (scalarDist < nearestScalarDist) { nearestScalarDist = scalarDist; nearestScalarClass = gals[j].highH; }
  }
  if (nearestClass === gals[i].highH) looCorrect++;
  if (nearestScalarClass === gals[i].highH) looScalarCorrect++;
}

const looAcc = looCorrect / N;
const looScalarAcc = looScalarCorrect / N;
console.log('\n  Full-state LOO-1NN accuracy: ' + (looAcc * 100).toFixed(1) + '% (' + looCorrect + '/' + N + ')');
console.log('  Scalar (hR) LOO-1NN accuracy: ' + (looScalarAcc * 100).toFixed(1) + '% (' + looScalarCorrect + '/' + N + ')');
console.log('  Chance level: 50%');
console.log('  Full-state vs scalar: ' + (looAcc > looScalarAcc ? 'FULL-STATE WINS' : 'SCALAR WINS'));

let looK3correct = 0;
let looK3scalarCorrect = 0;
for (let i = 0; i < N; i++) {
  const dists = [];
  const scalarDists = [];
  for (let j = 0; j < N; j++) {
    if (i === j) continue;
    dists.push({ dist: vectorDist(gals[i].state.fullVector, gals[j].state.fullVector), highH: gals[j].highH });
    scalarDists.push({ dist: Math.abs(gals[i].hR - gals[j].hR), highH: gals[j].highH });
  }
  dists.sort((a, b) => a.dist - b.dist);
  scalarDists.sort((a, b) => a.dist - b.dist);
  const k3vote = dists.slice(0, 3).filter(d => d.highH).length >= 2;
  const k3scalarVote = scalarDists.slice(0, 3).filter(d => d.highH).length >= 2;
  if (k3vote === gals[i].highH) looK3correct++;
  if (k3scalarVote === gals[i].highH) looK3scalarCorrect++;
}

const looK3acc = looK3correct / N;
const looK3scalarAcc = looK3scalarCorrect / N;
console.log('\n  Full-state LOO-3NN accuracy: ' + (looK3acc * 100).toFixed(1) + '%');
console.log('  Scalar LOO-3NN accuracy: ' + (looK3scalarAcc * 100).toFixed(1) + '%');


console.log('\n\n' + '#'.repeat(70));
console.log('8C-4: PROFILE-BIN CORRELATION SCAN');
console.log('Which specific radial bins carry DQ information?');
console.log('#'.repeat(70));

console.log('\n  V-PROFILE BIN CORRELATIONS WITH DQ:');
console.log('  ' + 'Bin'.padEnd(6) + 'r(DQ, vNorm[b])'.padEnd(20) + 't'.padEnd(8) + 'Sig?'.padEnd(8) + '| r(DQ, hFrac[b])'.padEnd(20) + 't'.padEnd(8) + 'Sig?');
console.log('  ' + '-'.repeat(80));

let sigVbins = 0, sigHbins = 0;
for (let b = 0; b < Nradial; b++) {
  const vBin = gals.map(g => g.state.vNorm[b]);
  const hBin = gals.map(g => g.state.haloFracProfile[b]);
  const rv = pearsonR(dqArr, vBin);
  const rh = pearsonR(dqArr, hBin);
  const tv = rv * Math.sqrt((N - 2) / (1 - rv * rv));
  const th = rh * Math.sqrt((N - 2) / (1 - rh * rh));
  const sigV = Math.abs(tv) > 2;
  const sigH = Math.abs(th) > 2;
  if (sigV) sigVbins++;
  if (sigH) sigHbins++;
  console.log('  ' + (b + 1 + '').padEnd(6) + ((rv >= 0 ? '+' : '') + rv.toFixed(3)).padEnd(20) + tv.toFixed(2).padEnd(8) + (sigV ? 'YES' : 'no').padEnd(8) + '| ' + ((rh >= 0 ? '+' : '') + rh.toFixed(3)).padEnd(20) + th.toFixed(2).padEnd(8) + (sigH ? 'YES' : 'no'));
}
console.log('\n  Significant V bins: ' + sigVbins + '/' + Nradial);
console.log('  Significant halo bins: ' + sigHbins + '/' + Nradial);


console.log('\n\n' + '#'.repeat(70));
console.log('8C-5: FULL-STATE CHANNEL ABSORPTION');
console.log('#'.repeat(70));

const top8bins = [];
for (let b = 0; b < Nradial; b++) {
  const rv = Math.abs(pearsonR(dqArr, gals.map(g => g.state.vNorm[b])));
  const rh = Math.abs(pearsonR(dqArr, gals.map(g => g.state.haloFracProfile[b])));
  top8bins.push({ b, r: Math.max(rv, rh), type: rv > rh ? 'v' : 'h' });
}
top8bins.sort((a, b) => b.r - a.r);
const selected = top8bins.slice(0, 8);

const selectedVec = gals.map(g => selected.map(s => s.type === 'v' ? g.state.vNorm[s.b] : g.state.haloFracProfile[s.b]));
const selectedPlusCtrl = gals.map((g, i) => [g.logK, g.dmFrac, g.envCode, ...selectedVec[i]]);
const selectedPlusCtrlHR = gals.map((g, i) => [g.logK, g.dmFrac, g.envCode, g.hR, ...selectedVec[i]]);

const absorptions = [
  { name: 'raw', r: r_VfR_a0R },
  { name: 'ctrl: logK+dmFrac+env', r: partialR(vfArr, a0Arr, gals.map(g => [g.logK, g.dmFrac, g.envCode])) },
  { name: 'ctrl: haloResponse only', r: partialR(vfArr, a0Arr, gals.map(g => [g.hR])) },
  { name: 'ctrl: top-8 profile bins', r: partialR(vfArr, a0Arr, selectedVec) },
  { name: 'ctrl: full+hR', r: partialR(vfArr, a0Arr, gals.map(g => [g.logK, g.dmFrac, g.envCode, g.hR])) },
  { name: 'ctrl: full+top8', r: partialR(vfArr, a0Arr, selectedPlusCtrl) },
  { name: 'ctrl: full+hR+top8', r: partialR(vfArr, a0Arr, selectedPlusCtrlHR) },
];

console.log('\n  ' + 'Controls'.padEnd(35) + 'partial r'.padEnd(12) + 'Absorbed');
console.log('  ' + '-'.repeat(60));
for (const a of absorptions) {
  const red = (1 - Math.abs(a.r) / Math.abs(r_VfR_a0R)) * 100;
  console.log('  ' + a.name.padEnd(35) + ((a.r >= 0 ? '+' : '') + a.r.toFixed(3)).padEnd(12) + red.toFixed(1) + '%');
}

const top8Abs = absorptions.find(t => t.name === 'ctrl: top-8 profile bins');
const top8Red = (1 - Math.abs(top8Abs.r) / Math.abs(r_VfR_a0R)) * 100;
const hRRed = (1 - Math.abs(absorptions[2].r) / Math.abs(r_VfR_a0R)) * 100;
const profileBeatsHR = top8Red > hRRed;


console.log('\n\n' + '#'.repeat(70));
console.log('8C-6: SPECTRAL SIGNATURE ANALYSIS');
console.log('Do high-H galaxies have a distinct Fourier spectrum?');
console.log('#'.repeat(70));

const q1Fourier = Array(5).fill(0);
const q5Fourier = Array(5).fill(0);
for (const g of Q1) for (let k = 0; k < 5; k++) q1Fourier[k] += g.state.fourierCoeffs[k] / Q1.length;
for (const g of Q5) for (let k = 0; k < 5; k++) q5Fourier[k] += g.state.fourierCoeffs[k] / Q5.length;

console.log('\n  FOURIER SPECTRUM (Q1 vs Q5):');
console.log('  ' + 'k'.padEnd(6) + 'Q1'.padEnd(12) + 'Q5'.padEnd(12) + 'Ratio'.padEnd(12) + 'r(DQ,Fk)');
console.log('  ' + '-'.repeat(55));
for (let k = 0; k < 5; k++) {
  const ratio = q5Fourier[k] > 0 ? q1Fourier[k] / q5Fourier[k] : 999;
  const rk = pearsonR(dqArr, gals.map(g => g.state.fourierCoeffs[k]));
  console.log('  ' + (k + 1 + '').padEnd(6) + q1Fourier[k].toFixed(4).padEnd(12) + q5Fourier[k].toFixed(4).padEnd(12) + ratio.toFixed(2).padEnd(12) + (rk >= 0 ? '+' : '') + rk.toFixed(3));
}


console.log('\n\n' + '='.repeat(70));
console.log('PROGRAM 8C GRAND VERDICT');
console.log('='.repeat(70));

const test1 = looAcc > 0.6 || looK3acc > 0.6;
const test2 = looAcc > looScalarAcc || looK3acc > looK3scalarAcc;
const test3 = profileBeatsHR;
const test4 = sigVbins >= 3 || sigHbins >= 3;

const allTests = [
  { name: 'T1: Full-state LOO accuracy > 60%', pass: test1 },
  { name: 'T2: Full-state beats scalar classification', pass: test2 },
  { name: 'T3: Profile bins absorb more channel than hR', pass: test3 },
  { name: 'T4: >= 3 profile bins individually significant', pass: test4 },
];

let totalPass = 0;
for (const t of allTests) {
  if (t.pass) totalPass++;
  console.log('  ' + (t.pass ? 'PASS' : 'FAIL').padEnd(6) + t.name);
}
console.log('  Total: ' + totalPass + '/4');

console.log('\n  CRITICAL COMPARISON:');
console.log('  Full-state LOO: ' + (looAcc * 100).toFixed(1) + '% vs scalar: ' + (looScalarAcc * 100).toFixed(1) + '%');
console.log('  Profile absorption: ' + top8Red.toFixed(1) + '% vs hR: ' + hRRed.toFixed(1) + '%');
console.log('  Significant bins: V=' + sigVbins + ', H=' + sigHbins);

if (totalPass >= 3) {
  console.log('\n  VERDICT: MAP-LEVEL BREAKTHROUGH');
  console.log('  The full rotation-curve state carries information about H');
  console.log('  that scalar summaries cannot capture. The hidden variable');
  console.log('  lives in the SHAPE of the curve, not in any single number.');
} else if (totalPass >= 2) {
  console.log('\n  VERDICT: PARTIAL MAP SIGNAL');
  console.log('  Some map-level information exists but does not decisively');
  console.log('  break the scalar ceiling.');
} else {
  console.log('\n  VERDICT: MAP CEILING = SCALAR CEILING');
  console.log('  The full rotation-curve state carries no more information');
  console.log('  about H than scalar summaries. The information barrier is');
  console.log('  not about summary statistics vs. profiles — it is about');
  console.log('  the fundamental dimensionality of 1D rotation curves.');
  console.log('  H requires TRULY 2D data: IFU velocity fields, lensing');
  console.log('  maps, or cosmological assembly histories.');
}


const outPath = path.join(__dirname, '..', 'public', 'program8c-map-reconstruction.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: '8C',
  title: 'Matched 2D Map Reconstruction',
  timestamp: new Date().toISOString(),
  N,
  stateDimension: gals[0].state.fullVector.length,
  r_VfR_a0R,
  classification: {
    loo1nn: { fullState: looAcc, scalar: looScalarAcc },
    loo3nn: { fullState: looK3acc, scalar: looK3scalarAcc },
  },
  absorption: absorptions.map(a => ({ name: a.name, r: a.r })),
  significantBins: { V: sigVbins, H: sigHbins },
  quintileSeparation: { vProfile: vProfileSep, haloProfile: hProfileSep },
  tests: { t1: test1, t2: test2, t3: test3, t4: test4, totalPass },
  verdict: totalPass >= 3 ? 'MAP BREAKTHROUGH' : totalPass >= 2 ? 'PARTIAL' : 'MAP CEILING = SCALAR CEILING',
}, null, 2));
console.log('\nSaved: ' + outPath);
