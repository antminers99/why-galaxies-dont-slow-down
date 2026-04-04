#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "31.0.0";
function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(80)); }

const H0 = 73.0;
const c_km = 299792.458;

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0(pts) {
  let lo = 2.0, hi = 5.0;
  for (let s = 0; s < 150; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gb = Math.pow(10, p.log_g_bar);
      c1 += (p.log_g_obs - Math.log10(mcgaughRAR(gb, Math.pow(10, m1)))) ** 2;
      c2 += (p.log_g_obs - Math.log10(mcgaughRAR(gb, Math.pow(10, m2)))) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  const logA0 = (lo + hi) / 2, a0 = Math.pow(10, logA0);
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.log_g_bar);
    ss += (p.log_g_obs - Math.log10(mcgaughRAR(gb, a0))) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length), ss, n: pts.length };
}

function predictSS(a0, pts) {
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gb, a0);
    ss += (p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10)) ** 2;
  }
  return ss;
}

function corrWith(x, y) {
  const n = x.length;
  if (n < 4) return 0;
  const mx = x.reduce((s, v) => s + v, 0) / n;
  const my = y.reduce((s, v) => s + v, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i]-mx)*(y[i]-my); sxx += (x[i]-mx)**2; syy += (y[i]-my)**2; }
  return sxy / Math.sqrt(sxx * syy + 1e-20);
}

function linReg(X, y) {
  const n = y.length, p = X[0] ? X[0].length : 0;
  if (p === 0) return { coefs: [], intercept: y.reduce((s,v)=>s+v,0)/n, residuals: y.map(v => v - y.reduce((s,v)=>s+v,0)/n) };
  const means = [];
  for (let j = 0; j < p; j++) { let s = 0; for (let i = 0; i < n; i++) s += X[i][j]; means.push(s / n); }
  const my = y.reduce((s,v)=>s+v,0) / n;
  const Xc = X.map(row => row.map((v,j) => v - means[j]));
  const yc = y.map(v => v - my);
  const XtX = Array.from({length:p}, () => Array(p).fill(0));
  const Xty = Array(p).fill(0);
  for (let i=0;i<n;i++) { for (let j=0;j<p;j++) { Xty[j]+=Xc[i][j]*yc[i]; for (let k=0;k<p;k++) XtX[j][k]+=Xc[i][j]*Xc[i][k]; } }
  for (let j=0;j<p;j++) XtX[j][j] += 1e-10;
  const A = XtX.map((r,i) => [...r, Xty[i]]);
  for (let j=0;j<p;j++) { let mx=j; for(let i=j+1;i<p;i++) if(Math.abs(A[i][j])>Math.abs(A[mx][j]))mx=i; [A[j],A[mx]]=[A[mx],A[j]]; for(let i=j+1;i<p;i++){const f=A[i][j]/A[j][j];for(let k=j;k<=p;k++)A[i][k]-=f*A[j][k];} }
  const b = Array(p).fill(0);
  for (let j=p-1;j>=0;j--) { b[j]=A[j][p]; for(let k=j+1;k<p;k++)b[j]-=A[j][k]*b[k]; b[j]/=A[j][j]; }
  const intercept = my - b.reduce((s,v,j) => s + v * means[j], 0);
  const residuals = y.map((yi, i) => yi - (intercept + b.reduce((s, c, j) => s + c * X[i][j], 0)));
  return { coefs: b, intercept, residuals };
}

function residualize(y, controls) {
  if (controls.length === 0) return [...y];
  const n = y.length;
  const X = [];
  for (let i = 0; i < n; i++) X.push(controls.map(c => c[i]));
  return linReg(X, y).residuals;
}

const p11 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase11-sensitivity-lab.json'), 'utf8'));
const sparcAll = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));
const p25 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase25-group-membership.json'), 'utf8'));

const envLookup = {};
for (const g of p25.galaxyAssignments) envLookup[g.name] = g;

const galaxyData = [];

for (const gal of p11.galaxies) {
  const pts = gal.localProfile.filter(p =>
    isFinite(p.log_g_bar) && isFinite(p.log_g_obs) && p.log_g_obs > p.log_g_bar * 0.5 && p.r > 0
  );
  if (pts.length < 5) continue;
  const sp = sparcAll.find(s => s.name === gal.name);
  if (!sp || !sp.MHI || sp.MHI <= 0) continue;
  const L = sp.L36 || sp.L || 0;
  if (L <= 0) continue;

  const fit = fitA0(pts);
  const logMHI = Math.log10(sp.MHI);
  const n = pts.length;
  const logVflat = Math.log10(sp.Vflat || 100);
  const D = sp.D;
  const eD = sp.eD || 0;

  const z = (H0 * D) / c_km;
  const logD = Math.log10(D);

  const rArr = pts.map(p => p.r);
  const rMax = Math.max(...rArr);
  const rMid = rMax / 2;
  const innerPts = pts.filter(p => p.r <= rMid);
  const outerPts = pts.filter(p => p.r > rMid);

  function residualVariability(subset) {
    if (subset.length < 3) return 0;
    const xv = subset.map(p => p.log_g_bar);
    const yv = subset.map(p => p.log_g_obs);
    const mx = xv.reduce((s,v)=>s+v,0)/xv.length;
    const my = yv.reduce((s,v)=>s+v,0)/yv.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < xv.length; i++) { sxy += (xv[i]-mx)*(yv[i]-my); sxx += (xv[i]-mx)**2; }
    const b = sxx > 0 ? sxy / sxx : 0;
    const a = my - b * mx;
    let ss = 0;
    for (let i = 0; i < xv.length; i++) ss += (yv[i] - (a + b * xv[i])) ** 2;
    return Math.sqrt(ss / xv.length);
  }

  const rcWiggliness = residualVariability(innerPts) + residualVariability(outerPts);

  const ev = envLookup[gal.name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms,
    n, logMHI, logVflat, rcWiggliness, envCode,
    D, eD, z, logD,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const zArr = galaxyData.map(g => g.z);
const logDArr = galaxyData.map(g => g.logD);
const DArr = galaxyData.map(g => g.D);
const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const vflatArr = galaxyData.map(g => g.logVflat);
const nPtsArr = galaxyData.map(g => g.n);

log("");
log("=".repeat(80));
log("  PHASE 31: REDSHIFT / COSMIC TIME / DISTANCE");
log("  Does distance (as proxy for cosmic time) predict a0 variation?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  NEW DOOR: COSMIC CONTEXT");
sep();
log("  Variables in this door:");
log("    1. redshift / cosmic time / distance  ← TESTING NOW");
log("    2. cosmic web environment             ← next");
log("");

log("  NOTE: SPARC is a LOCAL sample (z < 0.03).");
log("  Any redshift/distance signal here is NOT about cosmic evolution");
log("  of a0 over billions of years. Rather, it tests whether:");
log("  - Distance-dependent systematics contaminate the a0 measurement");
log("  - Selection effects at different distances bias the result");
log("  - The local Hubble flow introduces artifacts");
log("");

const Dmin = Math.min(...DArr), Dmax = Math.max(...DArr);
const zMin = Math.min(...zArr), zMax = Math.max(...zArr);
log("  DISTANCE RANGE:");
log("    D: " + Dmin.toFixed(1) + " — " + Dmax.toFixed(1) + " Mpc");
log("    z: " + zMin.toFixed(5) + " — " + zMax.toFixed(5));
log("    (all z << 0.1, deeply local universe)");
log("");

log("  PER-GALAXY (sorted by distance):");
const sorted = [...galaxyData].sort((a,b) => a.D - b.D);
for (const g of sorted) {
  log("    " + g.name.padEnd(16) + "D=" + g.D.toFixed(1).padStart(6) + " Mpc  z=" + g.z.toFixed(5) +
    "  logA0=" + g.logA0.toFixed(3) + "  n=" + g.n.toString().padStart(2));
}
log("");

log("=".repeat(80));
log("  TEST 1: RAW CORRELATIONS");
log("=".repeat(80));
log("");

for (const [label, arr] of [['D (Mpc)', DArr], ['log(D)', logDArr], ['z', zArr]]) {
  const r = corrWith(arr, deltaA0);
  const t = Math.abs(r) * Math.sqrt(N-2) / Math.sqrt(1-r**2+1e-10);
  log("  " + label.padEnd(14) + " vs delta_a0: r=" + r.toFixed(3) + " t=" + t.toFixed(1));
}
log("");

const rDnPts = corrWith(logDArr, nPtsArr);
const rDmhi = corrWith(logDArr, mhiArr);
const rDvflat = corrWith(logDArr, vflatArr);
const rDenv = corrWith(logDArr, envArr);
const rDwig = corrWith(logDArr, wigArr);

log("  CONFOUNDERS — log(D) correlates with:");
log("    n_points:     r=" + rDnPts.toFixed(3) + " (selection: distant → fewer points)");
log("    logMHI:       r=" + rDmhi.toFixed(3) + " (Malmquist: distant → more massive)");
log("    logVflat:     r=" + rDvflat.toFixed(3) + " (Malmquist: distant → faster)");
log("    envCode:      r=" + rDenv.toFixed(3));
log("    rcWiggliness: r=" + rDwig.toFixed(3));
log("");

log("=".repeat(80));
log("  TEST 2: CONFOUNDER STRIPPING");
log("=".repeat(80));
log("");

const controlSets = [
  { name: 'raw', vals: [] },
  { name: '|MHI', vals: [mhiArr] },
  { name: '|MHI+Vflat', vals: [mhiArr, vflatArr] },
  { name: '|MHI+wig', vals: [mhiArr, wigArr] },
  { name: '|MHI+wig+env', vals: [mhiArr, wigArr, envArr] },
  { name: '|MHI+Vflat+wig', vals: [mhiArr, vflatArr, wigArr] },
  { name: '|MHI+Vflat+wig+env', vals: [mhiArr, vflatArr, wigArr, envArr] },
  { name: '|nPts', vals: [nPtsArr] },
  { name: '|nPts+MHI', vals: [nPtsArr, mhiArr] },
];

for (const [label, arr] of [['log(D)', logDArr], ['z', zArr]]) {
  log("  " + label + " vs delta_a0:");
  log("  \u250c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510");
  log("  \u2502  Controls                 r_partial   |t|    status       \u2502");
  log("  \u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524");
  for (const cs of controlSets) {
    const residY = residualize(deltaA0, cs.vals);
    const residX = residualize(arr, cs.vals);
    const r = corrWith(residX, residY);
    const df = N - 2 - cs.vals.length;
    const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
    const status = t >= 2.0 ? "SURVIVES" : t >= 1.65 ? "MARGINAL" : "FAILS";
    log("  \u2502  " + cs.name.padEnd(24) + r.toFixed(3).padStart(8) + t.toFixed(1).padStart(7) + "    " + status.padEnd(10) + "\u2502");
  }
  log("  \u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518");
  log("");
}

log("=".repeat(80));
log("  TEST 3: DISTANCE AS SYSTEMATIC — n_points INTERACTION");
log("  (Distant galaxies have fewer data points → noisier a0 fits)");
log("=".repeat(80));
log("");

const nearIdx = galaxyData.map((g,i) => g.D < 20 ? i : -1).filter(i => i >= 0);
const midIdx = galaxyData.map((g,i) => g.D >= 20 && g.D < 60 ? i : -1).filter(i => i >= 0);
const farIdx = galaxyData.map((g,i) => g.D >= 60 ? i : -1).filter(i => i >= 0);

function groupStats(idx) {
  const la = idx.map(i => allLogA0[i]);
  const mu = la.reduce((s,v)=>s+v,0)/la.length;
  const sd = Math.sqrt(la.reduce((s,v)=>s+(v-mu)**2,0)/(la.length-1));
  const npts = idx.map(i => galaxyData[i].n);
  const avgN = npts.reduce((s,v)=>s+v,0)/npts.length;
  return { n: la.length, meanLogA0: mu, sdLogA0: sd, avgNpts: avgN };
}

const nearS = groupStats(nearIdx), midS = groupStats(midIdx), farS = groupStats(farIdx);

log("  Distance bins:");
log("  \u250c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510");
log("  \u2502  Bin           N    <logA0>   SD     <nPts>          \u2502");
log("  \u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524");
log("  \u2502  D<20 Mpc     " + nearS.n.toString().padStart(2) + "    " + nearS.meanLogA0.toFixed(3) + "   " + nearS.sdLogA0.toFixed(3) + "   " + nearS.avgNpts.toFixed(0).padStart(4) + "            \u2502");
log("  \u2502  20-60 Mpc    " + midS.n.toString().padStart(2) + "    " + midS.meanLogA0.toFixed(3) + "   " + midS.sdLogA0.toFixed(3) + "   " + midS.avgNpts.toFixed(0).padStart(4) + "            \u2502");
log("  \u2502  D>60 Mpc     " + farS.n.toString().padStart(2) + "    " + farS.meanLogA0.toFixed(3) + "   " + farS.sdLogA0.toFixed(3) + "   " + farS.avgNpts.toFixed(0).padStart(4) + "            \u2502");
log("  \u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518");
log("");

const nearMean = nearS.meanLogA0;
const farMean = farS.meanLogA0;
const tBins = Math.abs(nearMean - farMean) / Math.sqrt(nearS.sdLogA0**2/nearS.n + farS.sdLogA0**2/farS.n);
log("  Near vs Far mean difference: " + (nearMean - farMean).toFixed(3) + " dex, t=" + tBins.toFixed(1));
log("  " + (tBins >= 2.0 ? "SIGNIFICANT" : "NOT SIGNIFICANT"));
log("");

log("=".repeat(80));
log("  TEST 4: PERMUTATION TESTS");
log("=".repeat(80));
log("");

const NPERMS = 10000;
for (const [label, arr, ctrls] of [
  ['logD|raw', logDArr, []],
  ['logD|MHI', logDArr, [mhiArr]],
  ['logD|MHI+wig', logDArr, [mhiArr, wigArr]],
  ['logD|MHI+wig+env', logDArr, [mhiArr, wigArr, envArr]],
  ['logD|nPts', logDArr, [nPtsArr]],
]) {
  const residY = residualize(deltaA0, ctrls);
  const residX = residualize(arr, ctrls);
  const obsR = Math.abs(corrWith(residX, residY));
  let cnt = 0;
  const sh = [...residY];
  for (let p = 0; p < NPERMS; p++) {
    for (let i = sh.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [sh[i], sh[j]] = [sh[j], sh[i]];
    }
    if (Math.abs(corrWith(residX, sh)) >= obsR) cnt++;
  }
  log("  " + label.padEnd(22) + " |r|=" + obsR.toFixed(3) + " perm_p=" + (cnt/NPERMS).toFixed(4) +
    (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "**" : ""));
}
log("");

log("=".repeat(80));
log("  TEST 5: LOO CROSS-VALIDATION");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'envCode': envArr,
  'logD': logDArr, 'logVflat': vflatArr, 'nPts': nPtsArr,
};

function getVal(gIdx, fname) { return allFeatures[fname][gIdx]; }

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_D: log(D)', features: ['logD'] },
  { name: 'M_MHI', features: ['logMHI'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
  { name: 'M_MHI+wig+D', features: ['logMHI', 'rcWiggliness', 'logD'] },
  { name: 'M_MHI+wig+env', features: ['logMHI', 'rcWiggliness', 'envCode'] },
  { name: 'M_MHI+wig+env+D', features: ['logMHI', 'rcWiggliness', 'envCode', 'logD'] },
  { name: 'M6: Per-galaxy', features: ['__free__'] },
];

const cvResults = [];
for (const model of modelDefs) {
  let totalSS = 0, totalN = 0;
  if (model.features[0] === '__free__') {
    for (let i = 0; i < N; i++) { totalSS += predictSS(galaxyData[i].a0, galaxyData[i].pts); totalN += galaxyData[i].pts.length; }
    cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), k: N }); continue;
  }
  if (model.features.length === 0) {
    for (let i = 0; i < N; i++) { const trainY = allLogA0.filter((_, j) => j !== i); const mu = trainY.reduce((s,v) => s+v, 0) / trainY.length; totalSS += predictSS(Math.pow(10, mu), galaxyData[i].pts); totalN += galaxyData[i].pts.length; }
    cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), k: 1 }); continue;
  }
  for (let i = 0; i < N; i++) {
    const trainIdx = [...Array(N).keys()].filter(j => j !== i);
    const trainX = trainIdx.map(j => model.features.map(fn => getVal(j, fn)));
    const trainY = trainIdx.map(j => allLogA0[j]);
    const reg = linReg(trainX, trainY);
    const testX = model.features.map(fn => getVal(i, fn));
    let predLogA0 = reg.intercept + reg.coefs.reduce((s, c, j) => s + c * testX[j], 0);
    predLogA0 = Math.max(2.5, Math.min(4.5, predLogA0));
    totalSS += predictSS(Math.pow(10, predLogA0), galaxyData[i].pts);
    totalN += galaxyData[i].pts.length;
  }
  cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), k: model.features.length + 1 });
}

const m0rms = cvResults[0].cvRMS;
const m6rms = cvResults[cvResults.length - 1].cvRMS;
const gap = m0rms - m6rms;

log("  \u250c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510");
log("  \u2502  Model                    k    CV-RMS   vs M0   gap-closed             \u2502");
log("  \u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const gc = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  \u2502  " + r.name.padEnd(24) + r.k.toString().padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(8) + (gc + "%").padStart(13) + "            \u2502");
}
log("  \u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518");
log("");

const fn = (name) => cvResults.find(r => r.name === name);
const mhiWigGap = gap > 0 ? (m0rms - fn('M_MHI+wig').cvRMS) / gap * 100 : 0;
const mhiWigEnvGap = gap > 0 ? (m0rms - fn('M_MHI+wig+env').cvRMS) / gap * 100 : 0;
const mhiWigDGap = gap > 0 ? (m0rms - fn('M_MHI+wig+D').cvRMS) / gap * 100 : 0;
const mhiWigEnvDGap = gap > 0 ? (m0rms - fn('M_MHI+wig+env+D').cvRMS) / gap * 100 : 0;
const dAloneGap = gap > 0 ? (m0rms - fn('M_D: log(D)').cvRMS) / gap * 100 : 0;

log("  KEY COMPARISONS:");
log("    log(D) alone:               " + dAloneGap.toFixed(1) + "%");
log("    MHI + wig:                  " + mhiWigGap.toFixed(1) + "%");
log("    MHI + wig + D:              " + mhiWigDGap.toFixed(1) + "%");
log("    MHI + wig + env:            " + mhiWigEnvGap.toFixed(1) + "%");
log("    MHI + wig + env + D:        " + mhiWigEnvDGap.toFixed(1) + "%");
log("");

const dAdds = mhiWigDGap > mhiWigGap + 1;
const dAddsToEnv = mhiWigEnvDGap > mhiWigEnvGap + 1;
log("  D adds beyond MHI+wig? " + (dAdds ? "YES" : "NO"));
log("  D adds beyond MHI+wig+env? " + (dAddsToEnv ? "YES" : "NO"));
log("");

log("=".repeat(80));
log("  TEST 6: DISTANCE-DEPENDENT SYSTEMATICS CHECK");
log("  Does distance predict a0 SCATTER (not just mean)?");
log("=".repeat(80));
log("");

const rmsArr = galaxyData.map(g => g.rms);
const rDrms = corrWith(logDArr, rmsArr);
const tDrms = Math.abs(rDrms) * Math.sqrt(N-2) / Math.sqrt(1-rDrms**2+1e-10);
log("  log(D) vs RMS residual: r=" + rDrms.toFixed(3) + " t=" + tDrms.toFixed(1));
log("  " + (tDrms >= 2 ? "Distant galaxies have systematically different fit quality" : "No distance-dependent fit quality bias"));
log("");

const absDA0 = deltaA0.map(v => Math.abs(v));
const rDabsDA0 = corrWith(logDArr, absDA0);
const tDabsDA0 = Math.abs(rDabsDA0) * Math.sqrt(N-2) / Math.sqrt(1-rDabsDA0**2+1e-10);
log("  log(D) vs |delta_a0|: r=" + rDabsDA0.toFixed(3) + " t=" + tDabsDA0.toFixed(1));
log("  " + (tDabsDA0 >= 2 ? "Distance increases a0 SCATTER" : "No distance-dependent scatter inflation"));
log("");

log("=".repeat(80));
log("  PHASE 31 — STRICT CRITERIA");
log("=".repeat(80));
log("");

const rRaw = corrWith(logDArr, deltaA0);
const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);
const rAfterMHI = corrWith(residualize(logDArr, [mhiArr]), residualize(deltaA0, [mhiArr]));
const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);
const rAfterMHIWig = corrWith(residualize(logDArr, [mhiArr, wigArr]), residualize(deltaA0, [mhiArr, wigArr]));
const tAfterMHIWig = Math.abs(rAfterMHIWig) * Math.sqrt(N-4) / Math.sqrt(1-rAfterMHIWig**2+1e-10);

const criteria = [
  { name: 'Clear raw |t|>=2', pass: tRaw >= 2.0, val: 'r=' + rRaw.toFixed(3) + ' t=' + tRaw.toFixed(1) },
  { name: 'Independent of MHI', pass: tAfterMHI >= 1.65, val: 'r=' + rAfterMHI.toFixed(3) + ' t=' + tAfterMHI.toFixed(1) },
  { name: 'Independent of MHI+wig', pass: tAfterMHIWig >= 1.65, val: 'r=' + rAfterMHIWig.toFixed(3) + ' t=' + tAfterMHIWig.toFixed(1) },
  { name: 'LOO-CV > M0', pass: dAloneGap > 2, val: dAloneGap.toFixed(1) + '%' },
  { name: 'Adds to MHI+wig', pass: dAdds, val: (mhiWigDGap - mhiWigGap).toFixed(1) + '%' },
  { name: 'Adds beyond MHI+wig+env', pass: dAddsToEnv, val: (mhiWigEnvDGap - mhiWigEnvGap).toFixed(1) + '%' },
];

const nPassed = criteria.filter(c => c.pass).length;

log("  \u250c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2510");
log("  \u2502  Criterion                   Pass?    Value                     \u2502");
log("  \u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524");
for (const c of criteria) {
  log("  \u2502  " + c.name.padEnd(29) + (c.pass ? " YES " : " NO  ") + "  " + c.val.padEnd(27) + "\u2502");
}
log("  \u251c\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2524");
log("  \u2502  TOTAL: " + nPassed + "/6 criteria met" + " ".repeat(39) + "\u2502");
log("  \u2514\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2518");
log("");

log("=".repeat(80));
log("  PHASE 31 — VERDICT");
log("=".repeat(80));
log("");

if (nPassed >= 4) {
  log("  \u2554\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2557");
  log("  \u2551  REDSHIFT/DISTANCE shows significant independent signal.               \u2551");
  log("  \u255a\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u255d");
} else if (nPassed >= 2) {
  log("  \u2554\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2557");
  log("  \u2551  REDSHIFT/DISTANCE shows PARTIAL signal.                               \u2551");
  log("  \u255a\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u255d");
} else {
  log("  \u2554\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2557");
  log("  \u2551  REDSHIFT/DISTANCE FAILS.                                             \u2551");
  log("  \u2551  No cosmic-time / distance effect on a0 in this local sample.         \u2551");
  log("  \u255a\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u2550\u255d");
}
log("");

log("  IMPORTANT CONTEXT:");
log("  This SPARC sample spans z ~ 0.0001 to 0.025 (all local).");
log("  Over this tiny redshift range, cosmic evolution of a0 would be");
log("  negligible even in theories where a0 varies with time.");
log("  A null result here means no distance-dependent SYSTEMATICS,");
log("  NOT that a0 is constant over cosmic time.");
log("");
log("  Best model: MHI + rcWiggliness + envCode = 14.9% gap closure");
log("  (unchanged)");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  distanceRange: { Dmin: +Dmin.toFixed(1), Dmax: +Dmax.toFixed(1), zMin: +zMin.toFixed(5), zMax: +zMax.toFixed(5) },
  rawCorrelation: { r: +rRaw.toFixed(3), t: +tRaw.toFixed(1) },
  afterMHI: { r: +rAfterMHI.toFixed(3), t: +tAfterMHI.toFixed(1) },
  afterMHIWig: { r: +rAfterMHIWig.toFixed(3), t: +tAfterMHIWig.toFixed(1) },
  distanceBins: {
    near: { n: nearS.n, meanLogA0: +nearS.meanLogA0.toFixed(3), sd: +nearS.sdLogA0.toFixed(3) },
    mid: { n: midS.n, meanLogA0: +midS.meanLogA0.toFixed(3), sd: +midS.sdLogA0.toFixed(3) },
    far: { n: farS.n, meanLogA0: +farS.meanLogA0.toFixed(3), sd: +farS.sdLogA0.toFixed(3) },
  },
  distanceSystematics: {
    rDvsRMS: +rDrms.toFixed(3), rDvsAbsDeltaA0: +rDabsDA0.toFixed(3),
  },
  looGapClosed: {
    D: +dAloneGap.toFixed(1), mhiWig: +mhiWigGap.toFixed(1),
    mhiWigD: +mhiWigDGap.toFixed(1), mhiWigEnv: +mhiWigEnvGap.toFixed(1),
    mhiWigEnvD: +mhiWigEnvDGap.toFixed(1),
  },
  successCriteria: { passed: nPassed, total: 6 },
  verdict: nPassed >= 4 ? 'SIGNIFICANT' : nPassed >= 2 ? 'PARTIAL' : 'FAILED',
};

fs.writeFileSync(path.join(__dirname, '../public/phase31-redshift.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase31-redshift.json");
