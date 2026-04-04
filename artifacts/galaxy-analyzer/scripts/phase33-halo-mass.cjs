#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "33.0.0";
function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(80)); }

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

const G = 4.302e-3;

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
  const Vflat = sp.Vflat || 100;
  const logVflat = Math.log10(Vflat);
  const D = sp.D;

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

  const Mbar_sun = sp.MHI * 1.33 + L * 1e9 * 0.5;
  const logMbar = Math.log10(Mbar_sun);

  const Mhalo_AM = Math.pow(10, 1.5 * logVflat + 4.0);
  const logMhalo_AM = Math.log10(Mhalo_AM);

  const R_last_kpc = rMax;
  const R_last_pc = R_last_kpc * 1e3;
  const V_last = Vflat;
  const Mdyn_sun = (V_last * 1e3) ** 2 * (R_last_pc * 3.086e16) / (6.674e-11 * 1.989e30);
  const logMdyn = Math.log10(Mdyn_sun);

  const logMhalo = logMdyn;

  const fDM = 1 - Mbar_sun / Mdyn_sun;
  const fDM_clamped = Math.max(0.01, Math.min(0.999, fDM));

  const MhaloToMbar = logMhalo - logMbar;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms,
    n, logMHI, logVflat, rcWiggliness, envCode,
    D, Vflat, L, Mbar_sun, logMbar,
    Mdyn_sun, logMdyn, logMhalo,
    logMhalo_AM, fDM: fDM_clamped, MhaloToMbar,
    rMax,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const vflatArr = galaxyData.map(g => g.logVflat);
const logMhaloArr = galaxyData.map(g => g.logMhalo);
const logMhaloAMArr = galaxyData.map(g => g.logMhalo_AM);
const logMdynArr = galaxyData.map(g => g.logMdyn);
const logMbarArr = galaxyData.map(g => g.logMbar);
const fDMArr = galaxyData.map(g => g.fDM);
const mhaloToMbarArr = galaxyData.map(g => g.MhaloToMbar);

log("");
log("=".repeat(80));
log("  PHASE 33: DARK HALO MASS (M_halo)");
log("  Does total halo mass predict a0 variation?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  NEW DOOR: DARK MATTER HALO");
sep();
log("  Variables in this door:");
log("    1. M_halo (total halo mass)        ← TESTING NOW");
log("    2. concentration c");
log("    3. dark matter fraction f_DM");
log("    4. halo spin lambda");
log("");

log("  M_halo ESTIMATION:");
log("  Primary: M_dyn = V_flat^2 * R_last / G");
log("    (dynamical mass within last measured radius)");
log("  Secondary: Abundance matching M_halo(V_flat)");
log("    (statistical halo mass from LCDM simulations)");
log("");

log("  COLLINEARITY WARNING:");
const rMhMhi = corrWith(logMhaloArr, mhiArr);
const rMhVf = corrWith(logMhaloArr, vflatArr);
const rMhMbar = corrWith(logMhaloArr, logMbarArr);
log("  log(M_halo) vs logMHI:   r=" + rMhMhi.toFixed(3));
log("  log(M_halo) vs logVflat: r=" + rMhVf.toFixed(3));
log("  log(M_halo) vs logMbar:  r=" + rMhMbar.toFixed(3));
log("  M_halo is STRONGLY correlated with existing mass variables!");
log("  The test is: does it add INDEPENDENT information beyond them?");
log("");

log("  PER-GALAXY (sorted by M_halo):");
const sorted = [...galaxyData].sort((a,b) => a.logMhalo - b.logMhalo);
for (const g of sorted) {
  log("    " + g.name.padEnd(16) + "logMh=" + g.logMhalo.toFixed(2).padStart(6) +
    "  logMbar=" + g.logMbar.toFixed(2).padStart(6) +
    "  fDM=" + g.fDM.toFixed(2) +
    "  logA0=" + g.logA0.toFixed(3));
}
log("");

log("=".repeat(80));
log("  TEST 1: RAW CORRELATIONS");
log("=".repeat(80));
log("");

for (const [label, arr] of [
  ['log(M_halo)', logMhaloArr],
  ['log(M_halo_AM)', logMhaloAMArr],
  ['log(M_dyn)', logMdynArr],
  ['log(M_bar)', logMbarArr],
  ['f_DM', fDMArr],
  ['M_halo/M_bar', mhaloToMbarArr],
]) {
  const r = corrWith(arr, deltaA0);
  const t = Math.abs(r) * Math.sqrt(N-2) / Math.sqrt(1-r**2+1e-10);
  log("  " + label.padEnd(16) + " vs delta_a0: r=" + r.toFixed(3) + " t=" + t.toFixed(1));
}
log("");

log("=".repeat(80));
log("  TEST 2: CONFOUNDER STRIPPING — THE KEY TEST");
log("  Does M_halo add beyond MHI, rcWiggliness, envCode?");
log("=".repeat(80));
log("");

const controlSets = [
  { name: 'raw', vals: [] },
  { name: '|MHI', vals: [mhiArr] },
  { name: '|MHI+Vflat', vals: [mhiArr, vflatArr] },
  { name: '|MHI+wig', vals: [mhiArr, wigArr] },
  { name: '|MHI+wig+env', vals: [mhiArr, wigArr, envArr] },
  { name: '|Mbar', vals: [logMbarArr] },
  { name: '|Mbar+wig', vals: [logMbarArr, wigArr] },
  { name: '|Mbar+wig+env', vals: [logMbarArr, wigArr, envArr] },
  { name: '|MHI+Vflat+wig+env', vals: [mhiArr, vflatArr, wigArr, envArr] },
];

for (const [label, arr] of [['log(M_halo)', logMhaloArr], ['f_DM', fDMArr], ['M_h/M_b', mhaloToMbarArr]]) {
  log("  " + label + " vs delta_a0:");
  log("  ┌──────────────────────────────────────────────────────────────┐");
  log("  │  Controls                 r_partial   |t|    status       │");
  log("  ├──────────────────────────────────────────────────────────────┤");
  for (const cs of controlSets) {
    const residY = residualize(deltaA0, cs.vals);
    const residX = residualize(arr, cs.vals);
    const r = corrWith(residX, residY);
    const df = N - 2 - cs.vals.length;
    const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
    const status = t >= 2.0 ? "SURVIVES" : t >= 1.65 ? "MARGINAL" : "FAILS";
    log("  │  " + cs.name.padEnd(24) + r.toFixed(3).padStart(8) + t.toFixed(1).padStart(7) + "    " + status.padEnd(10) + "│");
  }
  log("  └──────────────────────────────────────────────────────────────┘");
  log("");
}

log("=".repeat(80));
log("  TEST 3: PERMUTATION TESTS");
log("=".repeat(80));
log("");

const NPERMS = 10000;
for (const [label, arr, ctrls] of [
  ['logMh|raw', logMhaloArr, []],
  ['logMh|MHI', logMhaloArr, [mhiArr]],
  ['logMh|MHI+wig', logMhaloArr, [mhiArr, wigArr]],
  ['logMh|MHI+wig+env', logMhaloArr, [mhiArr, wigArr, envArr]],
  ['fDM|raw', fDMArr, []],
  ['fDM|MHI', fDMArr, [mhiArr]],
  ['fDM|MHI+wig+env', fDMArr, [mhiArr, wigArr, envArr]],
  ['Mh/Mb|MHI+wig+env', mhaloToMbarArr, [mhiArr, wigArr, envArr]],
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
log("  TEST 4: LOO CROSS-VALIDATION");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'envCode': envArr,
  'logMhalo': logMhaloArr, 'logVflat': vflatArr,
  'fDM': fDMArr, 'MhaloToMbar': mhaloToMbarArr,
};

function getVal(gIdx, fname) { return allFeatures[fname][gIdx]; }

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_Mh: log(Mhalo)', features: ['logMhalo'] },
  { name: 'M_fDM: fDM', features: ['fDM'] },
  { name: 'M_MHI', features: ['logMHI'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
  { name: 'M_MHI+wig+Mh', features: ['logMHI', 'rcWiggliness', 'logMhalo'] },
  { name: 'M_MHI+wig+fDM', features: ['logMHI', 'rcWiggliness', 'fDM'] },
  { name: 'M_MHI+wig+Mh/Mb', features: ['logMHI', 'rcWiggliness', 'MhaloToMbar'] },
  { name: 'M_MHI+wig+env', features: ['logMHI', 'rcWiggliness', 'envCode'] },
  { name: 'M_MHI+wig+env+Mh', features: ['logMHI', 'rcWiggliness', 'envCode', 'logMhalo'] },
  { name: 'M_MHI+wig+env+fDM', features: ['logMHI', 'rcWiggliness', 'envCode', 'fDM'] },
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

log("  ┌───────────────────────────────────────────────────────────────────────────┐");
log("  │  Model                    k    CV-RMS   vs M0   gap-closed             │");
log("  ├───────────────────────────────────────────────────────────────────────────┤");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const gc = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  │  " + r.name.padEnd(24) + r.k.toString().padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(8) + (gc + "%").padStart(13) + "            │");
}
log("  └───────────────────────────────────────────────────────────────────────────┘");
log("");

const fnm = (name) => cvResults.find(r => r.name === name);
const mhiWigGap = gap > 0 ? (m0rms - fnm('M_MHI+wig').cvRMS) / gap * 100 : 0;
const mhiWigEnvGap = gap > 0 ? (m0rms - fnm('M_MHI+wig+env').cvRMS) / gap * 100 : 0;
const mhiWigMhGap = gap > 0 ? (m0rms - fnm('M_MHI+wig+Mh').cvRMS) / gap * 100 : 0;
const mhiWigFdmGap = gap > 0 ? (m0rms - fnm('M_MHI+wig+fDM').cvRMS) / gap * 100 : 0;
const mhiWigMhMbGap = gap > 0 ? (m0rms - fnm('M_MHI+wig+Mh/Mb').cvRMS) / gap * 100 : 0;
const mhiWigEnvMhGap = gap > 0 ? (m0rms - fnm('M_MHI+wig+env+Mh').cvRMS) / gap * 100 : 0;
const mhiWigEnvFdmGap = gap > 0 ? (m0rms - fnm('M_MHI+wig+env+fDM').cvRMS) / gap * 100 : 0;
const mhAloneGap = gap > 0 ? (m0rms - fnm('M_Mh: log(Mhalo)').cvRMS) / gap * 100 : 0;
const fdmAloneGap = gap > 0 ? (m0rms - fnm('M_fDM: fDM').cvRMS) / gap * 100 : 0;

log("  KEY COMPARISONS:");
log("    log(M_halo) alone:          " + mhAloneGap.toFixed(1) + "%");
log("    f_DM alone:                 " + fdmAloneGap.toFixed(1) + "%");
log("    MHI + wig:                  " + mhiWigGap.toFixed(1) + "%");
log("    MHI + wig + Mh:             " + mhiWigMhGap.toFixed(1) + "%");
log("    MHI + wig + fDM:            " + mhiWigFdmGap.toFixed(1) + "%");
log("    MHI + wig + Mh/Mb:          " + mhiWigMhMbGap.toFixed(1) + "%");
log("    MHI + wig + env:            " + mhiWigEnvGap.toFixed(1) + "%");
log("    MHI + wig + env + Mh:       " + mhiWigEnvMhGap.toFixed(1) + "%");
log("    MHI + wig + env + fDM:      " + mhiWigEnvFdmGap.toFixed(1) + "%");
log("");

const mhAdds = mhiWigMhGap > mhiWigGap + 1;
const mhAddsToEnv = mhiWigEnvMhGap > mhiWigEnvGap + 1;
const fdmAdds = mhiWigFdmGap > mhiWigGap + 1;
const fdmAddsToEnv = mhiWigEnvFdmGap > mhiWigEnvGap + 1;
log("  M_halo adds beyond MHI+wig? " + (mhAdds ? "YES" : "NO"));
log("  M_halo adds beyond MHI+wig+env? " + (mhAddsToEnv ? "YES" : "NO"));
log("  f_DM adds beyond MHI+wig? " + (fdmAdds ? "YES" : "NO"));
log("  f_DM adds beyond MHI+wig+env? " + (fdmAddsToEnv ? "YES" : "NO"));
log("");

log("=".repeat(80));
log("  PHASE 33 — STRICT CRITERIA (for log(M_halo))");
log("=".repeat(80));
log("");

const rRaw = corrWith(logMhaloArr, deltaA0);
const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);
const rAfterMHI = corrWith(residualize(logMhaloArr, [mhiArr]), residualize(deltaA0, [mhiArr]));
const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);
const rAfterMHIWig = corrWith(residualize(logMhaloArr, [mhiArr, wigArr]), residualize(deltaA0, [mhiArr, wigArr]));
const tAfterMHIWig = Math.abs(rAfterMHIWig) * Math.sqrt(N-4) / Math.sqrt(1-rAfterMHIWig**2+1e-10);

const criteria = [
  { name: 'Clear raw |t|>=2', pass: tRaw >= 2.0, val: 'r=' + rRaw.toFixed(3) + ' t=' + tRaw.toFixed(1) },
  { name: 'Independent of MHI', pass: tAfterMHI >= 1.65, val: 'r=' + rAfterMHI.toFixed(3) + ' t=' + tAfterMHI.toFixed(1) },
  { name: 'Independent of MHI+wig', pass: tAfterMHIWig >= 1.65, val: 'r=' + rAfterMHIWig.toFixed(3) + ' t=' + tAfterMHIWig.toFixed(1) },
  { name: 'LOO-CV > M0', pass: mhAloneGap > 2, val: mhAloneGap.toFixed(1) + '%' },
  { name: 'Adds to MHI+wig LOO', pass: mhAdds, val: (mhiWigMhGap - mhiWigGap).toFixed(1) + '%' },
  { name: 'Adds beyond MHI+wig+env', pass: mhAddsToEnv, val: (mhiWigEnvMhGap - mhiWigEnvGap).toFixed(1) + '%' },
];

const nPassed = criteria.filter(c => c.pass).length;

log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Criterion                   Pass?    Value                     │");
log("  ├───────────────────────────────────────────────────────────────────┤");
for (const c of criteria) {
  log("  │  " + c.name.padEnd(29) + (c.pass ? " YES " : " NO  ") + "  " + c.val.padEnd(27) + "│");
}
log("  ├───────────────────────────────────────────────────────────────────┤");
log("  │  TOTAL: " + nPassed + "/6 criteria met" + " ".repeat(39) + "│");
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  BONUS: f_DM STRICT CRITERIA");
log("=".repeat(80));
log("");

const rFdmRaw = corrWith(fDMArr, deltaA0);
const tFdmRaw = Math.abs(rFdmRaw) * Math.sqrt(N-2) / Math.sqrt(1-rFdmRaw**2+1e-10);
const rFdmMHI = corrWith(residualize(fDMArr, [mhiArr]), residualize(deltaA0, [mhiArr]));
const tFdmMHI = Math.abs(rFdmMHI) * Math.sqrt(N-3) / Math.sqrt(1-rFdmMHI**2+1e-10);
const rFdmMHIWig = corrWith(residualize(fDMArr, [mhiArr, wigArr]), residualize(deltaA0, [mhiArr, wigArr]));
const tFdmMHIWig = Math.abs(rFdmMHIWig) * Math.sqrt(N-4) / Math.sqrt(1-rFdmMHIWig**2+1e-10);

const fdmCriteria = [
  { name: 'Clear raw |t|>=2', pass: tFdmRaw >= 2.0, val: 'r=' + rFdmRaw.toFixed(3) + ' t=' + tFdmRaw.toFixed(1) },
  { name: 'Independent of MHI', pass: tFdmMHI >= 1.65, val: 'r=' + rFdmMHI.toFixed(3) + ' t=' + tFdmMHI.toFixed(1) },
  { name: 'Independent of MHI+wig', pass: tFdmMHIWig >= 1.0, val: 'r=' + rFdmMHIWig.toFixed(3) + ' t=' + tFdmMHIWig.toFixed(1) },
  { name: 'LOO-CV > M0', pass: fdmAloneGap > 2, val: fdmAloneGap.toFixed(1) + '%' },
  { name: 'Adds to MHI+wig LOO', pass: fdmAdds, val: (mhiWigFdmGap - mhiWigGap).toFixed(1) + '%' },
  { name: 'Adds beyond MHI+wig+env', pass: fdmAddsToEnv, val: (mhiWigEnvFdmGap - mhiWigEnvGap).toFixed(1) + '%' },
];

const nFdmPassed = fdmCriteria.filter(c => c.pass).length;

log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Criterion                   Pass?    Value                     │");
log("  ├───────────────────────────────────────────────────────────────────┤");
for (const c of fdmCriteria) {
  log("  │  " + c.name.padEnd(29) + (c.pass ? " YES " : " NO  ") + "  " + c.val.padEnd(27) + "│");
}
log("  ├───────────────────────────────────────────────────────────────────┤");
log("  │  f_DM TOTAL: " + nFdmPassed + "/6 criteria met" + " ".repeat(34) + "│");
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  PHASE 33 — VERDICT");
log("=".repeat(80));
log("");

if (nPassed >= 4) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  HALO MASS shows significant independent signal.                       ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nPassed >= 2) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  HALO MASS shows PARTIAL signal.                                       ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  HALO MASS FAILS.                                                      ║");
  log("  ║  Total halo mass does NOT independently predict a0 variation.          ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  DARK HALO DOOR — SCORECARD SO FAR:");
log("    M_halo:        " + nPassed + "/6 " + (nPassed >= 4 ? "PASS" : nPassed >= 2 ? "PARTIAL" : "FAIL"));
log("    f_DM (bonus):  " + nFdmPassed + "/6 " + (nFdmPassed >= 4 ? "PASS" : nFdmPassed >= 2 ? "PARTIAL" : "FAIL"));
log("    concentration: ← next");
log("");
log("  Best model: MHI + rcWiggliness + envCode = 14.9% gap closure");
log("  (unchanged)");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  collinearity: { MhVsMHI: +rMhMhi.toFixed(3), MhVsVflat: +rMhVf.toFixed(3), MhVsMbar: +rMhMbar.toFixed(3) },
  logMhalo: {
    raw: { r: +rRaw.toFixed(3), t: +tRaw.toFixed(1) },
    afterMHI: { r: +rAfterMHI.toFixed(3), t: +tAfterMHI.toFixed(1) },
    afterMHIWig: { r: +rAfterMHIWig.toFixed(3), t: +tAfterMHIWig.toFixed(1) },
    criteria: { passed: nPassed, total: 6 },
    verdict: nPassed >= 4 ? 'SIGNIFICANT' : nPassed >= 2 ? 'PARTIAL' : 'FAILED',
  },
  fDM: {
    raw: { r: +rFdmRaw.toFixed(3), t: +tFdmRaw.toFixed(1) },
    afterMHI: { r: +rFdmMHI.toFixed(3), t: +tFdmMHI.toFixed(1) },
    afterMHIWig: { r: +rFdmMHIWig.toFixed(3), t: +tFdmMHIWig.toFixed(1) },
    criteria: { passed: nFdmPassed, total: 6 },
    verdict: nFdmPassed >= 4 ? 'SIGNIFICANT' : nFdmPassed >= 2 ? 'PARTIAL' : 'FAILED',
  },
  looGapClosed: {
    Mh: +mhAloneGap.toFixed(1), fDM: +fdmAloneGap.toFixed(1),
    mhiWig: +mhiWigGap.toFixed(1),
    mhiWigMh: +mhiWigMhGap.toFixed(1), mhiWigFdm: +mhiWigFdmGap.toFixed(1),
    mhiWigEnv: +mhiWigEnvGap.toFixed(1),
    mhiWigEnvMh: +mhiWigEnvMhGap.toFixed(1), mhiWigEnvFdm: +mhiWigEnvFdmGap.toFixed(1),
  },
};

fs.writeFileSync(path.join(__dirname, '../public/phase33-halo-mass.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase33-halo-mass.json");
