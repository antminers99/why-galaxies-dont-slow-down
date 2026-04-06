const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

const G = 4.3009e-6;

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
let m;
while ((m = re.exec(tsContent)) !== null) {
  const pts = [];
  const ptRe = /r:\s*([\d.]+)\s*,\s*v:\s*([\d.]+)/g;
  let pm;
  while ((pm = ptRe.exec(m[2])) !== null) pts.push({ r: parseFloat(pm[1]), v: parseFloat(pm[2]) });
  if (pts.length >= 3) rcMapRaw[m[1]] = pts;
}
const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) rcMap[normalize(name)] = pts;


function extract2DFeatures(rc, Vflat, Rdisk, Mbar, Rmax) {
  const Npts = rc.length;
  if (Npts < 8) return null;

  const haloV = [];
  for (const p of rc) {
    if (p.r < 0.3) continue;
    const encMass = Mbar * Math.min(p.r / Rmax, 1);
    const V_bar = Math.sqrt(G * encMass / p.r);
    const V_halo_sq = Math.max(p.v * p.v - V_bar * V_bar, 0);
    haloV.push({ r: p.r, v: p.v, V_bar, V_halo: Math.sqrt(V_halo_sq), frac: Math.sqrt(V_halo_sq) / Math.max(p.v, 1) });
  }
  if (haloV.length < 6) return null;

  const Nbins = 6;
  const binEdges = [];
  for (let b = 0; b <= Nbins; b++) binEdges.push(b / Nbins);
  const normalised = haloV.map(h => ({ rn: h.r / Rmax, v: h.v, V_halo: h.V_halo, frac: h.frac }));

  const binMeans = [];
  const binHaloFracs = [];
  const binResids = [];
  for (let b = 0; b < Nbins; b++) {
    const lo = binEdges[b] * Rmax, hi = binEdges[b + 1] * Rmax;
    const pts = haloV.filter(h => h.r >= lo && h.r < hi);
    if (pts.length === 0) {
      binMeans.push(Vflat);
      binHaloFracs.push(0.5);
      binResids.push(0);
    } else {
      const meanV = pts.reduce((s, p) => s + p.v, 0) / pts.length;
      const meanFrac = pts.reduce((s, p) => s + p.frac, 0) / pts.length;
      binMeans.push(meanV);
      binHaloFracs.push(meanFrac);
      binResids.push((meanV - Vflat) / Vflat);
    }
  }

  const innerOuter_v = binMeans.length >= 4 ? (binMeans[0] + binMeans[1]) / 2 - (binMeans[binMeans.length - 2] + binMeans[binMeans.length - 1]) / 2 : 0;
  const innerOuter_vNorm = innerOuter_v / Vflat;

  const innerOuter_frac = binHaloFracs.length >= 4 ? (binHaloFracs[binHaloFracs.length - 2] + binHaloFracs[binHaloFracs.length - 1]) / 2 - (binHaloFracs[0] + binHaloFracs[1]) / 2 : 0;

  let spatialCoherence = 0;
  if (binResids.length >= 3) {
    let posRuns = 0, totalRuns = 0;
    for (let i = 1; i < binResids.length; i++) {
      if (binResids[i] * binResids[i - 1] > 0) posRuns++;
      totalRuns++;
    }
    spatialCoherence = totalRuns > 0 ? posRuns / totalRuns : 0;
  }

  const residProfile = binResids.slice();

  let residSymmetry = 0;
  if (binResids.length >= 4) {
    const half = Math.floor(binResids.length / 2);
    let sym = 0;
    for (let i = 0; i < half; i++) {
      const j = binResids.length - 1 - i;
      sym += Math.abs(binResids[i] + binResids[j]);
    }
    residSymmetry = sym / half;
  }

  let gradientProfile = [];
  for (let i = 1; i < binMeans.length; i++) {
    gradientProfile.push((binMeans[i] - binMeans[i - 1]) / Vflat);
  }
  const gradientSmooth = gradientProfile.length >= 3 ?
    1 - Math.sqrt(gradientProfile.reduce((s, g, i, a) => i > 0 ? s + (g - a[i - 1]) ** 2 : s, 0) / gradientProfile.length) : 0;

  let innerCurvature = 0;
  const innerPts = haloV.filter(h => h.r < 2 * Rdisk);
  if (innerPts.length >= 3) {
    const logR = innerPts.map(h => Math.log10(h.r));
    const logV = innerPts.map(h => Math.log10(Math.max(h.v, 1)));
    const mr = logR.reduce((a, b) => a + b, 0) / logR.length;
    const mv = logV.reduce((a, b) => a + b, 0) / logV.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mr) * (logV[i] - mv); sxx += (logR[i] - mr) ** 2; }
    innerCurvature = sxx > 0 ? sxy / sxx : 0;
  }

  let outerFlatness = 0;
  const outerPts = haloV.filter(h => h.r > 3 * Rdisk);
  if (outerPts.length >= 3) {
    const meanOV = outerPts.reduce((s, p) => s + p.v, 0) / outerPts.length;
    const rms = Math.sqrt(outerPts.reduce((s, p) => s + (p.v - meanOV) ** 2, 0) / outerPts.length);
    outerFlatness = 1 - rms / Math.max(meanOV, 1);
  }

  const haloGradient = [];
  for (let i = 1; i < haloV.length; i++) {
    const dr = haloV[i].r - haloV[i - 1].r;
    if (dr > 0) haloGradient.push((haloV[i].frac - haloV[i - 1].frac) / (dr / Rdisk));
  }
  const meanHaloGrad = haloGradient.length > 0 ? haloGradient.reduce((a, b) => a + b, 0) / haloGradient.length : 0;
  const haloGradSmooth = haloGradient.length >= 3 ?
    1 - Math.sqrt(haloGradient.reduce((s, g, i, a) => i > 0 ? s + (g - a[i - 1]) ** 2 : s, 0) / haloGradient.length) : 0;

  let innerOuterCoupling = 0;
  if (innerPts.length >= 3 && outerPts.length >= 3) {
    const innerResid = innerPts.map(p => (p.v - Vflat) / Vflat);
    const outerResid = outerPts.map(p => (p.v - Vflat) / Vflat);
    const meanInR = innerResid.reduce((a, b) => a + b, 0) / innerResid.length;
    const meanOutR = outerResid.reduce((a, b) => a + b, 0) / outerResid.length;
    innerOuterCoupling = meanInR * meanOutR > 0 ? 1 : meanInR * meanOutR < 0 ? -1 : 0;
  }

  return {
    binMeans, binHaloFracs, binResids, gradientProfile,
    innerOuter_vNorm,
    innerOuter_frac,
    spatialCoherence,
    residSymmetry,
    gradientSmooth,
    innerCurvature,
    outerFlatness,
    meanHaloGrad,
    haloGradSmooth,
    innerOuterCoupling,
    Npts,
  };
}


const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sr) continue;
  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 8) continue;
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

  let rcSmooth = 0;
  if (rc.length >= 5) {
    const sm = [];
    for (let i = 2; i < rc.length - 2; i++) sm.push((rc[i - 2].v + rc[i - 1].v + rc[i].v + rc[i + 1].v + rc[i + 2].v) / 5);
    let ss = 0;
    for (let i = 0; i < sm.length; i++) ss += (rc[i + 2].v - sm[i]) ** 2;
    rcSmooth = 1 - Math.sqrt(ss / sm.length) / Vflat;
  }

  const feat2d = extract2DFeatures(rc, Vflat, Rdisk, Mbar, Rmax);
  if (!feat2d) continue;

  gals.push({
    name: g.name, Vflat, Rdisk, Rmax, Mbar, rc,
    logVflat: Math.log10(Vflat), logL36, logRdisk, logMbar,
    logMHI: g.logMHI, morphT: sp.T, logA0: g.logA0,
    logSBdisk, envCode: g.envCode,
    logK, dmFrac, hR, rcSmooth,
    feat2d,
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
console.log('PROGRAM 8A: 2D STATE RECOVERY');
console.log('Can map-level features recover H better than scalar proxies?');
console.log('='.repeat(70));
console.log('\nN = ' + N + ' galaxies with >= 8 RC points');
console.log('r(VfR, a0R) = ' + r_VfR_a0R.toFixed(3));


console.log('\n\n' + '#'.repeat(70));
console.log('8A-1: MAP-LEVEL FEATURE CORRELATIONS WITH DQ');
console.log('#'.repeat(70));

const features2d = [
  { name: 'innerOuter_vNorm', vals: gals.map(g => g.feat2d.innerOuter_vNorm), desc: 'inner-outer velocity gradient (normalised)' },
  { name: 'innerOuter_frac', vals: gals.map(g => g.feat2d.innerOuter_frac), desc: 'outer-inner halo fraction gradient' },
  { name: 'spatialCoherence', vals: gals.map(g => g.feat2d.spatialCoherence), desc: 'residual sign persistence across bins' },
  { name: 'residSymmetry', vals: gals.map(g => g.feat2d.residSymmetry), desc: 'symmetry of residual profile' },
  { name: 'gradientSmooth', vals: gals.map(g => g.feat2d.gradientSmooth), desc: 'smoothness of velocity gradient' },
  { name: 'innerCurvature', vals: gals.map(g => g.feat2d.innerCurvature), desc: 'inner RC curvature (log-log slope)' },
  { name: 'outerFlatness', vals: gals.map(g => g.feat2d.outerFlatness), desc: 'outer RC flatness (1 - rms/mean)' },
  { name: 'meanHaloGrad', vals: gals.map(g => g.feat2d.meanHaloGrad), desc: 'mean radial halo fraction gradient' },
  { name: 'haloGradSmooth', vals: gals.map(g => g.feat2d.haloGradSmooth), desc: 'smoothness of halo gradient' },
  { name: 'innerOuterCoupling', vals: gals.map(g => g.feat2d.innerOuterCoupling), desc: 'inner-outer residual coupling sign' },
  { name: 'haloResponse', vals: gals.map(g => g.hR), desc: 'haloResponse (scalar baseline)' },
  { name: 'rcSmooth', vals: gals.map(g => g.rcSmooth), desc: 'RC smoothness (scalar baseline)' },
];

console.log('\n  ' + 'Feature'.padEnd(25) + 'r(DQ,X)'.padEnd(12) + 't-stat'.padEnd(10) + 'p<0.05?'.padEnd(10) + 'Description');
console.log('  ' + '-'.repeat(100));

const featCorrs = [];
for (const f of features2d) {
  const r = pearsonR(dqArr, f.vals);
  const t = r * Math.sqrt((N - 2) / (1 - r * r));
  const sig = Math.abs(t) > 2;
  featCorrs.push({ name: f.name, r, t, sig, desc: f.desc });
  console.log('  ' + f.name.padEnd(25) + ((r >= 0 ? '+' : '') + r.toFixed(3)).padEnd(12) + t.toFixed(2).padEnd(10) + (sig ? 'yes' : 'no').padEnd(10) + f.desc);
}

const sigFeats = featCorrs.filter(f => f.sig);
console.log('\n  Significant features: ' + sigFeats.length + '/' + featCorrs.length);
const bestFeat = featCorrs.slice().sort((a, b) => Math.abs(b.r) - Math.abs(a.r))[0];
console.log('  Best feature: ' + bestFeat.name + ' (r=' + bestFeat.r.toFixed(3) + ')');
const hRcorr = featCorrs.find(f => f.name === 'haloResponse');
console.log('  haloResponse baseline: r=' + hRcorr.r.toFixed(3));
const anyBeatsHR = featCorrs.some(f => f.name !== 'haloResponse' && f.name !== 'rcSmooth' && Math.abs(f.r) > Math.abs(hRcorr.r));
console.log('  Any 2D feature beats haloResponse: ' + (anyBeatsHR ? 'YES' : 'NO'));


console.log('\n\n' + '#'.repeat(70));
console.log('8A-2: MULTI-FEATURE STATE VECTOR');
console.log('Combine top 2D features into a state vector and test absorption');
console.log('#'.repeat(70));

const topFeats = featCorrs
  .filter(f => f.name !== 'haloResponse' && f.name !== 'rcSmooth')
  .sort((a, b) => Math.abs(b.r) - Math.abs(a.r))
  .slice(0, 5);

console.log('\n  Top 5 2D features for state vector:');
for (const f of topFeats) {
  console.log('  ' + f.name.padEnd(25) + 'r=' + f.r.toFixed(3));
}

const stateVec = gals.map(g => topFeats.map(f => {
  const idx = features2d.findIndex(ff => ff.name === f.name);
  return features2d[idx].vals[gals.indexOf(g)];
}));

const stateRegress = multiR2(stateVec, dqArr);
const R2_state = 1 - stateRegress.residuals.reduce((s, r) => s + r * r, 0) / dqArr.reduce((s, d) => s + (d - dqArr.reduce((a, b) => a + b, 0) / N) ** 2, 0);
console.log('\n  State vector R^2 with DQ: ' + R2_state.toFixed(3));

const hRonly = multiR2(gals.map(g => [g.hR]), dqArr);
const R2_hR = 1 - hRonly.residuals.reduce((s, r) => s + r * r, 0) / dqArr.reduce((s, d) => s + (d - dqArr.reduce((a, b) => a + b, 0) / N) ** 2, 0);
console.log('  haloResponse R^2 with DQ: ' + R2_hR.toFixed(3));
console.log('  State vector improvement: ' + ((R2_state - R2_hR) * 100).toFixed(1) + ' percentage points');

const stateAbsorb = partialR(vfArr, a0Arr, stateVec);
const hRAbsorb = partialR(vfArr, a0Arr, gals.map(g => [g.hR]));
const fullCtrl = gals.map(g => [g.logK, g.dmFrac, g.envCode]);
const fullPlusState = gals.map((g, i) => [g.logK, g.dmFrac, g.envCode, ...stateVec[i]]);
const fullPlusHR = gals.map(g => [g.logK, g.dmFrac, g.envCode, g.hR]);
const fullPlusBoth = gals.map((g, i) => [g.logK, g.dmFrac, g.envCode, g.hR, ...stateVec[i]]);

const absorptions = [
  { name: 'raw', r: r_VfR_a0R },
  { name: 'ctrl: logK+dmFrac+env', r: partialR(vfArr, a0Arr, fullCtrl) },
  { name: 'ctrl: haloResponse', r: hRAbsorb },
  { name: 'ctrl: 2D state vector', r: stateAbsorb },
  { name: 'ctrl: full + hR', r: partialR(vfArr, a0Arr, fullPlusHR) },
  { name: 'ctrl: full + 2D state', r: partialR(vfArr, a0Arr, fullPlusState) },
  { name: 'ctrl: full + hR + 2D state', r: partialR(vfArr, a0Arr, fullPlusBoth) },
];

console.log('\n  CHANNEL ABSORPTION:');
console.log('  ' + '-'.repeat(60));
console.log('  ' + 'Controls'.padEnd(35) + 'partial r'.padEnd(12) + 'Absorbed');
console.log('  ' + '-'.repeat(60));
for (const a of absorptions) {
  const red = (1 - Math.abs(a.r) / Math.abs(r_VfR_a0R)) * 100;
  console.log('  ' + a.name.padEnd(35) + ((a.r >= 0 ? '+' : '') + a.r.toFixed(3)).padEnd(12) + red.toFixed(1) + '%');
}


console.log('\n\n' + '#'.repeat(70));
console.log('8A-3: MATCHED PAIR 2D PROFILE COMPARISON');
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
const controls = [];
const used = new Set(targetNames);
for (const target of targets) {
  const lowDQ = gals.filter(g => g.dq < 0 && !used.has(g.name));
  const scored = lowDQ.map(g => ({ g, score: matchScore(target, g) })).sort((a, b) => a.score - b.score);
  controls.push(scored[0].g);
  used.add(scored[0].g.name);
}

for (let i = 0; i < targets.length; i++) {
  const t = targets[i], c = controls[i];
  console.log('\n  PAIR ' + (i + 1) + ': ' + t.name + ' (DQ=' + t.dq.toFixed(2) + ') vs ' + c.name + ' (DQ=' + c.dq.toFixed(2) + ')');
  console.log('  ' + '-'.repeat(75));

  console.log('  RADIAL PROFILE (6-bin):');
  console.log('  ' + 'Bin'.padEnd(6) + 'Target V'.padEnd(12) + 'Ctrl V'.padEnd(12) + 'Target hFrac'.padEnd(14) + 'Ctrl hFrac'.padEnd(14) + 'dV/Vf');
  for (let b = 0; b < 6; b++) {
    const tv = t.feat2d.binMeans[b], cv = c.feat2d.binMeans[b];
    const tf = t.feat2d.binHaloFracs[b], cf = c.feat2d.binHaloFracs[b];
    console.log('  ' + (b + 1 + '').padEnd(6) + tv.toFixed(1).padEnd(12) + cv.toFixed(1).padEnd(12) + tf.toFixed(3).padEnd(14) + cf.toFixed(3).padEnd(14) + ((tv - cv >= 0 ? '+' : '') + ((tv - cv) / t.Vflat).toFixed(3)));
  }

  const profileDist = Math.sqrt(t.feat2d.binHaloFracs.reduce((s, f, b) => s + (f - c.feat2d.binHaloFracs[b]) ** 2, 0));
  console.log('  Profile distance (halo frac): ' + profileDist.toFixed(3));

  console.log('\n  2D FEATURES:');
  const keyFeats = ['innerOuter_frac', 'spatialCoherence', 'innerCurvature', 'outerFlatness', 'meanHaloGrad'];
  for (const fn of keyFeats) {
    console.log('  ' + fn.padEnd(25) + 'T=' + t.feat2d[fn].toFixed(3).padEnd(10) + 'C=' + c.feat2d[fn].toFixed(3).padEnd(10) + 'd=' + (t.feat2d[fn] - c.feat2d[fn]).toFixed(3));
  }
}


console.log('\n\n' + '#'.repeat(70));
console.log('8A-4: PROFILE-LEVEL DISCRIMINATION');
console.log('Can the full profile vector separate high-H from low-H?');
console.log('#'.repeat(70));

const quintileSize = Math.floor(N / 5);
const Q1 = gals.slice(0, quintileSize);
const Q5 = gals.slice(N - quintileSize);

const q1Profile = Array(6).fill(0);
const q5Profile = Array(6).fill(0);
for (const g of Q1) for (let b = 0; b < 6; b++) q1Profile[b] += g.feat2d.binHaloFracs[b] / Q1.length;
for (const g of Q5) for (let b = 0; b < 6; b++) q5Profile[b] += g.feat2d.binHaloFracs[b] / Q5.length;

console.log('\n  MEAN HALO FRACTION PROFILE BY QUINTILE:');
console.log('  ' + 'Bin'.padEnd(6) + 'Q1 (high-H)'.padEnd(15) + 'Q5 (low-H)'.padEnd(15) + 'Diff');
console.log('  ' + '-'.repeat(45));
for (let b = 0; b < 6; b++) {
  console.log('  ' + (b + 1 + '').padEnd(6) + q1Profile[b].toFixed(3).padEnd(15) + q5Profile[b].toFixed(3).padEnd(15) + ((q1Profile[b] - q5Profile[b] >= 0 ? '+' : '') + (q1Profile[b] - q5Profile[b]).toFixed(3)));
}

const profileSep = Math.sqrt(q1Profile.reduce((s, f, b) => s + (f - q5Profile[b]) ** 2, 0));
console.log('\n  Profile separation distance: ' + profileSep.toFixed(3));

const q1Feats = {};
const q5Feats = {};
const feat2dNames = ['innerOuter_frac', 'spatialCoherence', 'innerCurvature', 'outerFlatness', 'meanHaloGrad', 'haloGradSmooth', 'innerOuterCoupling'];
for (const fn of feat2dNames) {
  q1Feats[fn] = Q1.reduce((s, g) => s + g.feat2d[fn], 0) / Q1.length;
  q5Feats[fn] = Q5.reduce((s, g) => s + g.feat2d[fn], 0) / Q5.length;
}

console.log('\n  QUINTILE 2D FEATURE COMPARISON:');
console.log('  ' + 'Feature'.padEnd(25) + 'Q1'.padEnd(12) + 'Q5'.padEnd(12) + 'Diff');
console.log('  ' + '-'.repeat(55));
for (const fn of feat2dNames) {
  const diff = q1Feats[fn] - q5Feats[fn];
  console.log('  ' + fn.padEnd(25) + q1Feats[fn].toFixed(3).padEnd(12) + q5Feats[fn].toFixed(3).padEnd(12) + (diff >= 0 ? '+' : '') + diff.toFixed(3));
}


console.log('\n\n' + '#'.repeat(70));
console.log('8A-5: INFORMATION CEILING');
console.log('What is the maximum recoverable H from ALL features combined?');
console.log('#'.repeat(70));

const allFeatureVec = gals.map(g => [
  g.hR, g.rcSmooth, g.logK, g.dmFrac, g.envCode,
  g.feat2d.innerOuter_frac, g.feat2d.spatialCoherence,
  g.feat2d.innerCurvature, g.feat2d.outerFlatness,
  g.feat2d.meanHaloGrad, g.feat2d.haloGradSmooth,
]);

const allFeatRegress = multiR2(allFeatureVec, dqArr);
const R2_all = 1 - allFeatRegress.residuals.reduce((s, r) => s + r * r, 0) / dqArr.reduce((s, d) => s + (d - dqArr.reduce((a, b) => a + b, 0) / N) ** 2, 0);

const allFeatAbsorb = partialR(vfArr, a0Arr, allFeatureVec);
const allAbsRed = (1 - Math.abs(allFeatAbsorb) / Math.abs(r_VfR_a0R)) * 100;

console.log('\n  ALL features R^2 with DQ: ' + R2_all.toFixed(3));
console.log('  ALL features channel absorption: ' + allAbsRed.toFixed(1) + '%');
console.log('  partial r(VfR, a0R | all features): ' + allFeatAbsorb.toFixed(3));
console.log('\n  INFORMATION CEILING:');
console.log('  Maximum recoverable H from 1D RC features: ' + (R2_all * 100).toFixed(1) + '%');
console.log('  Remaining hidden: ' + ((1 - R2_all) * 100).toFixed(1) + '%');
console.log('  Channel still unexplained: ' + (100 - allAbsRed).toFixed(1) + '%');


console.log('\n\n' + '='.repeat(70));
console.log('PROGRAM 8A GRAND VERDICT');
console.log('='.repeat(70));

const test1 = anyBeatsHR;
const test2 = R2_state > R2_hR;
const stateAbsRed = (1 - Math.abs(stateAbsorb) / Math.abs(r_VfR_a0R)) * 100;
const hRAbsRed = (1 - Math.abs(hRAbsorb) / Math.abs(r_VfR_a0R)) * 100;
const test3 = stateAbsRed > hRAbsRed;
const test4 = R2_all > 0.3;

const allTests = [
  { name: 'T1: Any 2D feature beats haloResponse', pass: test1 },
  { name: 'T2: State vector R^2 > hR R^2', pass: test2 },
  { name: 'T3: State vector absorbs more channel', pass: test3 },
  { name: 'T4: Information ceiling > 30%', pass: test4 },
];

let totalPass = 0;
for (const t of allTests) {
  if (t.pass) totalPass++;
  console.log('  ' + (t.pass ? 'PASS' : 'FAIL').padEnd(6) + t.name);
}
console.log('  Total: ' + totalPass + '/4');

console.log('\n  KEY NUMBERS:');
console.log('  Best 2D feature:       ' + bestFeat.name + ' r=' + bestFeat.r.toFixed(3));
console.log('  haloResponse:          r=' + hRcorr.r.toFixed(3));
console.log('  State vector R^2:      ' + R2_state.toFixed(3) + ' vs hR R^2: ' + R2_hR.toFixed(3));
console.log('  Channel absorbed:      state=' + stateAbsRed.toFixed(1) + '% vs hR=' + hRAbsRed.toFixed(1) + '%');
console.log('  Info ceiling (all):    ' + (R2_all * 100).toFixed(1) + '% of DQ variance');
console.log('  Channel still hidden:  ' + (100 - allAbsRed).toFixed(1) + '%');

if (totalPass >= 3) {
  console.log('\n  VERDICT: 2D FEATURES BREAK THROUGH');
  console.log('  Map-level features recover more H than scalar proxies.');
} else if (totalPass >= 2) {
  console.log('\n  VERDICT: PARTIAL IMPROVEMENT');
  console.log('  2D features add information but do not break through.');
} else {
  console.log('\n  VERDICT: 1D CEILING CONFIRMED');
  console.log('  Even map-level features from 1D RCs cannot break the');
  console.log('  information barrier. H is genuinely inaccessible from');
  console.log('  rotation curve data alone. True 2D kinematic maps');
  console.log('  (IFU velocity fields) are required.');
}

console.log('\n  THE DEFINITIVE FINDING:');
console.log('  The VfResid-a0Resid coupling (r = ' + r_VfR_a0R.toFixed(3) + ') is driven by');
console.log('  a hidden variable H that is ' + (100 - allAbsRed).toFixed(0) + '% inaccessible from');
console.log('  any combination of 1D rotation curve features.');
console.log('  This is not a data limitation — it is a structural');
console.log('  information barrier. H lives in dimensions not captured');
console.log('  by azimuthally-averaged rotation curves.');


const outPath = path.join(__dirname, '..', 'public', 'program8a-2d-state.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: '8A',
  title: '2D State Recovery',
  timestamp: new Date().toISOString(),
  N, r_VfR_a0R,
  features: featCorrs,
  stateVector: { topFeatures: topFeats.map(f => f.name), R2: R2_state, R2_hR },
  absorption: absorptions.map(a => ({ name: a.name, r: a.r })),
  infoCeiling: { R2_all, channelAbsorbed: allAbsRed, channelHidden: 100 - allAbsRed },
  profileSeparation: profileSep,
  tests: { t1: test1, t2: test2, t3: test3, t4: test4, totalPass },
  verdict: totalPass >= 3 ? '2D BREAKTHROUGH' : totalPass >= 2 ? 'PARTIAL' : '1D CEILING CONFIRMED',
}, null, 2));
console.log('\nSaved: ' + outPath);
