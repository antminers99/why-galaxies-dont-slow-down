#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "33b.0.0";
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
  const Mdyn_sun = (Vflat * 1e3) ** 2 * (rMax * 1e3 * 3.086e16) / (6.674e-11 * 1.989e30);
  const logMhalo = Math.log10(Mdyn_sun);
  const logRlast = Math.log10(rMax);

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms,
    n, logMHI, logVflat, rcWiggliness, envCode,
    logMhalo, logMbar, logRlast, Vflat, rMax, D,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const mhaloArr = galaxyData.map(g => g.logMhalo);
const vflatArr = galaxyData.map(g => g.logVflat);
const rlastArr = galaxyData.map(g => g.logRlast);
const mbarArr = galaxyData.map(g => g.logMbar);

const BASELINE = [mhiArr, wigArr, envArr];

log("");
log("=".repeat(80));
log("  PHASE 33b: M_halo CIRCULARITY & STABILITY CHECK");
log("  Is the M_halo suppressor signal real or a measurement artifact?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  THE CIRCULARITY CONCERN:");
log("  M_halo = V_flat^2 * R_last / G");
log("  a0 is fitted from the SAME rotation curve");
log("  Both share the velocity measurements at large radii.");
log("  If V_flat has noise, it could inflate BOTH M_halo AND a0.");
log("");

log("=".repeat(80));
log("  TEST 1: DECOMPOSITION — V_flat vs R_last");
log("  logMhalo = 2*logVflat + logRlast + const");
log("  Which component carries the signal?");
log("=".repeat(80));
log("");

const rVfRl = corrWith(vflatArr, rlastArr);
log("  Collinearity: logVflat vs logRlast: r=" + rVfRl.toFixed(3));
log("");

for (const [label, arr] of [['logMhalo', mhaloArr], ['logVflat', vflatArr], ['logRlast', rlastArr]]) {
  log("  " + label + " vs delta_a0:");
  for (const [cname, ctrls] of [
    ['raw', []], ['|MHI', [mhiArr]], ['|MHI+wig', [mhiArr, wigArr]], ['|baseline', BASELINE],
  ]) {
    const rx = residualize(arr, ctrls);
    const ry = residualize(deltaA0, ctrls);
    const r = corrWith(rx, ry);
    const df = N - 2 - ctrls.length;
    const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
    log("    " + cname.padEnd(14) + " r=" + r.toFixed(3) + "  t=" + t.toFixed(1) +
      (t >= 2.0 ? " SURVIVES" : t >= 1.65 ? " MARGINAL" : " fails"));
  }
  log("");
}

log("  KEY QUESTION: If the signal is ONLY in V_flat, circularity is more");
log("  concerning (V_flat enters a0 directly). If R_last also carries signal,");
log("  it's less circular (R_last is geometric, not dynamical).");
log("");

log("=".repeat(80));
log("  TEST 2: CONTROLLING FOR V_flat DIRECTLY");
log("  If M_halo signal is just V_flat noise, it should vanish when");
log("  we control for V_flat.");
log("=".repeat(80));
log("");

for (const [label, ctrls] of [
  ['|Vflat', [vflatArr]],
  ['|Vflat+MHI', [vflatArr, mhiArr]],
  ['|Vflat+MHI+wig', [vflatArr, mhiArr, wigArr]],
  ['|Vflat+MHI+wig+env', [vflatArr, mhiArr, wigArr, envArr]],
]) {
  const rx = residualize(mhaloArr, ctrls);
  const ry = residualize(deltaA0, ctrls);
  const r = corrWith(rx, ry);
  const df = N - 2 - ctrls.length;
  const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
  log("  logMhalo " + label.padEnd(26) + " r=" + r.toFixed(3) + " t=" + t.toFixed(1) +
    (t >= 2.0 ? " SURVIVES" : t >= 1.65 ? " MARGINAL" : " fails"));
}
log("");
log("  If M_halo still adds after controlling V_flat,");
log("  the signal is NOT just shared V_flat noise.");
log("");

log("=".repeat(80));
log("  TEST 3: R_last ALONE (GEOMETRIC, NON-CIRCULAR)");
log("  R_last is the extent of the rotation curve, not a velocity.");
log("  If it carries independent signal, circularity is less worrying.");
log("=".repeat(80));
log("");

for (const [label, ctrls] of [
  ['|MHI', [mhiArr]],
  ['|MHI+Vflat', [mhiArr, vflatArr]],
  ['|MHI+wig', [mhiArr, wigArr]],
  ['|baseline', BASELINE],
  ['|baseline+Vflat', [...BASELINE, vflatArr]],
]) {
  const rx = residualize(rlastArr, ctrls);
  const ry = residualize(deltaA0, ctrls);
  const r = corrWith(rx, ry);
  const df = N - 2 - ctrls.length;
  const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
  log("  logRlast " + label.padEnd(26) + " r=" + r.toFixed(3) + " t=" + t.toFixed(1) +
    (t >= 2.0 ? " SURVIVES" : t >= 1.65 ? " MARGINAL" : " fails"));
}
log("");

log("=".repeat(80));
log("  TEST 4: JACKKNIFE STABILITY");
log("  Remove one galaxy at a time. Is the suppressor stable?");
log("=".repeat(80));
log("");

const jackR = [];
const jackNames = [];
for (let i = 0; i < N; i++) {
  const idx = [...Array(N).keys()].filter(j => j !== i);
  const subMh = idx.map(j => mhaloArr[j]);
  const subDa = idx.map(j => deltaA0[j]);
  const subBL = BASELINE.map(arr => idx.map(j => arr[j]));
  const rx = residualize(subMh, subBL);
  const ry = residualize(subDa, subBL);
  const r = corrWith(rx, ry);
  jackR.push(r);
  jackNames.push(galaxyData[i].name);
}

const jackMean = jackR.reduce((s,v)=>s+v,0)/N;
const jackSD = Math.sqrt(jackR.map(v => (v - jackMean)**2).reduce((s,v)=>s+v,0)/(N-1));

log("  Jackknife partial r (M_halo|baseline): N=" + N + " iterations");
log("  Mean r = " + jackMean.toFixed(3));
log("  SD = " + jackSD.toFixed(3));
log("  Range: [" + Math.min(...jackR).toFixed(3) + ", " + Math.max(...jackR).toFixed(3) + "]");
log("");

const fullR = corrWith(residualize(mhaloArr, BASELINE), residualize(deltaA0, BASELINE));
const influential = [];
for (let i = 0; i < N; i++) {
  const shift = Math.abs(jackR[i] - fullR);
  if (shift > 2 * jackSD) {
    influential.push({ name: jackNames[i], r: jackR[i], shift });
  }
}

if (influential.length === 0) {
  log("  No galaxy shifts r by >2 SD. Signal is STABLE.");
} else {
  log("  Influential galaxies (shift > 2 SD):");
  for (const g of influential) {
    log("    " + g.name.padEnd(16) + " r_without=" + g.r.toFixed(3) + " shift=" + g.shift.toFixed(3));
  }
}

const signFlips = jackR.filter(r => Math.sign(r) !== Math.sign(fullR)).length;
log("  Sign flips: " + signFlips + "/" + N + (signFlips === 0 ? " (perfectly stable)" : " (UNSTABLE!)"));
log("");

log("=".repeat(80));
log("  TEST 5: BOOTSTRAP CI ON PARTIAL CORRELATION");
log("  Does the suppressor effect's 95% CI exclude zero?");
log("=".repeat(80));
log("");

const NBOOT = 2000;
const bootR_base = [];
const bootR_vflat = [];
for (let b = 0; b < NBOOT; b++) {
  const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
  const bMh = idx.map(i => mhaloArr[i]);
  const bDa = idx.map(i => deltaA0[i]);

  const bBL = BASELINE.map(arr => idx.map(i => arr[i]));
  const rx1 = residualize(bMh, bBL);
  const ry1 = residualize(bDa, bBL);
  bootR_base.push(corrWith(rx1, ry1));

  const bBLV = [...BASELINE, vflatArr].map(arr => idx.map(i => arr[i]));
  const rx2 = residualize(bMh, bBLV);
  const ry2 = residualize(bDa, bBLV);
  bootR_vflat.push(corrWith(rx2, ry2));
}

bootR_base.sort((a,b) => a - b);
bootR_vflat.sort((a,b) => a - b);
const ci_base_lo = bootR_base[Math.floor(NBOOT * 0.025)];
const ci_base_hi = bootR_base[Math.floor(NBOOT * 0.975)];
const ci_vf_lo = bootR_vflat[Math.floor(NBOOT * 0.025)];
const ci_vf_hi = bootR_vflat[Math.floor(NBOOT * 0.975)];

log("  logMhalo|baseline:      95% CI = [" + ci_base_lo.toFixed(3) + ", " + ci_base_hi.toFixed(3) + "]" +
  (ci_base_lo > 0 ? " → EXCLUDES ZERO ✓" : " → crosses zero ✗"));
log("  logMhalo|baseline+Vflat: 95% CI = [" + ci_vf_lo.toFixed(3) + ", " + ci_vf_hi.toFixed(3) + "]" +
  (ci_vf_lo > 0 ? " → EXCLUDES ZERO ✓" : ci_vf_hi < 0 ? " → EXCLUDES ZERO ✓" : " → crosses zero ✗"));
log("");

log("=".repeat(80));
log("  TEST 6: PERMUTATION p — AFTER CONTROLLING V_flat");
log("  Most stringent test: does M_halo add beyond baseline+Vflat?");
log("=".repeat(80));
log("");

const NPERMS = 10000;
for (const [label, ctrls] of [
  ['|baseline', BASELINE],
  ['|baseline+Vflat', [...BASELINE, vflatArr]],
  ['|MHI+Vflat', [mhiArr, vflatArr]],
  ['|MHI+wig+Vflat', [mhiArr, wigArr, vflatArr]],
]) {
  const rx = residualize(mhaloArr, ctrls);
  const ry = residualize(deltaA0, ctrls);
  const obsR = Math.abs(corrWith(rx, ry));
  let cnt = 0;
  const sh = [...ry];
  for (let p = 0; p < NPERMS; p++) {
    for (let i = sh.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [sh[i], sh[j]] = [sh[j], sh[i]];
    }
    if (Math.abs(corrWith(rx, sh)) >= obsR) cnt++;
  }
  log("  logMhalo " + label.padEnd(22) + " |r|=" + obsR.toFixed(3) + " perm_p=" + (cnt/NPERMS).toFixed(4) +
    (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "**" : ""));
}
log("");

log("=".repeat(80));
log("  TEST 7: STRATIFICATION");
log("  Does the signal hold in subsamples?");
log("=".repeat(80));
log("");

const medianMHI = [...mhiArr].sort((a,b) => a - b)[Math.floor(N/2)];
const lowMHI = galaxyData.map((g,i) => i).filter(i => mhiArr[i] < medianMHI);
const highMHI = galaxyData.map((g,i) => i).filter(i => mhiArr[i] >= medianMHI);

for (const [label, idx] of [['Low MHI (N=' + lowMHI.length + ')', lowMHI], ['High MHI (N=' + highMHI.length + ')', highMHI]]) {
  const subMh = idx.map(i => mhaloArr[i]);
  const subDa = idx.map(i => deltaA0[i]);
  const subBL = BASELINE.map(arr => idx.map(i => arr[i]));
  const rx = residualize(subMh, subBL);
  const ry = residualize(subDa, subBL);
  const r = corrWith(rx, ry);
  const df = idx.length - 5;
  const t = df > 0 ? Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10) : 0;
  log("  " + label.padEnd(20) + " r=" + r.toFixed(3) + " t=" + t.toFixed(1) +
    (r > 0 ? " (same sign ✓)" : " (OPPOSITE sign ✗)"));
}

const fieldIdx = galaxyData.map((g,i) => i).filter(i => envArr[i] === 0);
const groupIdx = galaxyData.map((g,i) => i).filter(i => envArr[i] > 0);

for (const [label, idx] of [['Field (N=' + fieldIdx.length + ')', fieldIdx], ['Group/cluster (N=' + groupIdx.length + ')', groupIdx]]) {
  if (idx.length < 8) { log("  " + label + " — too few galaxies"); continue; }
  const subMh = idx.map(i => mhaloArr[i]);
  const subDa = idx.map(i => deltaA0[i]);
  const subMHI = idx.map(i => mhiArr[i]);
  const subWig = idx.map(i => wigArr[i]);
  const rx = residualize(subMh, [subMHI, subWig]);
  const ry = residualize(subDa, [subMHI, subWig]);
  const r = corrWith(rx, ry);
  const df = idx.length - 4;
  const t = df > 0 ? Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10) : 0;
  log("  " + label.padEnd(26) + " r=" + r.toFixed(3) + " t=" + t.toFixed(1) +
    (r > 0 ? " (same sign ✓)" : " (OPPOSITE sign ✗)"));
}
log("");

log("=".repeat(80));
log("  TEST 8: LOO — M_halo vs V_flat as predictor");
log("  If circularity drives everything, V_flat should do just as well.");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'envCode': envArr,
  'logMhalo': mhaloArr, 'logVflat': vflatArr, 'logRlast': rlastArr,
};
function getVal(gIdx, fname) { return allFeatures[fname][gIdx]; }

function looGapClosed(featureNames) {
  let totalSS_model = 0, totalSS_m0 = 0, totalSS_free = 0, totalN = 0;
  for (let i = 0; i < N; i++) {
    const trainIdx = [...Array(N).keys()].filter(j => j !== i);
    const trainY = trainIdx.map(j => allLogA0[j]);
    const mu = trainY.reduce((s,v)=>s+v,0)/trainY.length;
    totalSS_m0 += predictSS(Math.pow(10, mu), galaxyData[i].pts);
    totalSS_free += predictSS(galaxyData[i].a0, galaxyData[i].pts);
    totalN += galaxyData[i].pts.length;
    if (featureNames.length > 0) {
      const trainX = trainIdx.map(j => featureNames.map(fn => getVal(j, fn)));
      const reg = linReg(trainX, trainY);
      const testX = featureNames.map(fn => getVal(i, fn));
      let pred = reg.intercept + reg.coefs.reduce((s, c, j) => s + c * testX[j], 0);
      pred = Math.max(2.5, Math.min(4.5, pred));
      totalSS_model += predictSS(Math.pow(10, pred), galaxyData[i].pts);
    }
  }
  const m0rms = Math.sqrt(totalSS_m0 / totalN);
  const m6rms = Math.sqrt(totalSS_free / totalN);
  const gap = m0rms - m6rms;
  if (featureNames.length === 0) return { gc: 0, m0rms, m6rms, gap };
  const modelRms = Math.sqrt(totalSS_model / totalN);
  return { gc: gap > 0 ? (m0rms - modelRms) / gap * 100 : 0 };
}

const models = [
  ['baseline (MHI+wig+env)', ['logMHI', 'rcWiggliness', 'envCode']],
  ['baseline+Mhalo', ['logMHI', 'rcWiggliness', 'envCode', 'logMhalo']],
  ['baseline+Vflat', ['logMHI', 'rcWiggliness', 'envCode', 'logVflat']],
  ['baseline+Rlast', ['logMHI', 'rcWiggliness', 'envCode', 'logRlast']],
  ['baseline+Vflat+Rlast', ['logMHI', 'rcWiggliness', 'envCode', 'logVflat', 'logRlast']],
];

log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Model                         LOO gap-closed                  │");
log("  ├───────────────────────────────────────────────────────────────────┤");
for (const [name, feats] of models) {
  const gc = looGapClosed(feats);
  log("  │  " + name.padEnd(32) + (gc.gc.toFixed(1) + "%").padStart(8) + " ".repeat(24) + "│");
}
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

log("  INTERPRETATION:");
log("  If baseline+Vflat ≈ baseline+Mhalo → signal is mostly V_flat (MORE circular)");
log("  If baseline+Rlast adds substantially → signal has geometric component (LESS circular)");
log("  If baseline+Mhalo >> baseline+Vflat → combination V²R carries new info");
log("");

log("=".repeat(80));
log("  PHASE 33b — CIRCULARITY VERDICT");
log("=".repeat(80));
log("");

const vflatAfterBL = (() => {
  const rx = residualize(vflatArr, BASELINE);
  const ry = residualize(deltaA0, BASELINE);
  return corrWith(rx, ry);
})();
const rlastAfterBL = (() => {
  const rx = residualize(rlastArr, BASELINE);
  const ry = residualize(deltaA0, BASELINE);
  return corrWith(rx, ry);
})();
const mhAfterBLVf = (() => {
  const rx = residualize(mhaloArr, [...BASELINE, vflatArr]);
  const ry = residualize(deltaA0, [...BASELINE, vflatArr]);
  return corrWith(rx, ry);
})();

const vflatCarriesSignal = Math.abs(vflatAfterBL) > 0.15;
const rlastCarriesSignal = Math.abs(rlastAfterBL) > 0.15;
const mhSurvivesVflat = Math.abs(mhAfterBLVf) > 0.15;
const jackStable = signFlips === 0 && jackSD < 0.10;
const bootstrapExcludes = ci_base_lo > 0;
const bootstrapVfExcludes = ci_vf_lo > 0 || ci_vf_hi < 0;

const circTests = [
  { name: 'Jackknife stable (no sign flips)', pass: jackStable },
  { name: 'Bootstrap CI|baseline excludes 0', pass: bootstrapExcludes },
  { name: 'Signal NOT only Vflat (Mh|BL+Vf)', pass: mhSurvivesVflat },
  { name: 'R_last carries independent signal', pass: rlastCarriesSignal },
  { name: 'Bootstrap CI|BL+Vf excludes 0', pass: bootstrapVfExcludes },
];

const nCircPass = circTests.filter(c => c.pass).length;

log("  CIRCULARITY / STABILITY SCORECARD:");
log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Test                             Pass?                        │");
log("  ├───────────────────────────────────────────────────────────────────┤");
for (const c of circTests) {
  log("  │  " + c.name.padEnd(37) + (c.pass ? " YES" : " NO ") + " ".repeat(23) + "│");
}
log("  ├───────────────────────────────────────────────────────────────────┤");
log("  │  CIRCULARITY SCORE: " + nCircPass + "/5" + " ".repeat(40) + "│");
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

if (nCircPass >= 4) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  M_halo: CONFIRMED                                                    ║");
  log("  ║  Suppressor signal is REAL, survives circularity checks.               ║");
  log("  ║  Can be added to baseline with confidence.                             ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nCircPass >= 3) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  M_halo: CONFIRMED WITH CAVEAT                                        ║");
  log("  ║  Signal is likely real but some circularity concern remains.            ║");
  log("  ║  Add to baseline but flag the caveat.                                  ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nCircPass >= 2) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  M_halo: PARTIAL — circularity NOT resolved                            ║");
  log("  ║  Signal may be partially or fully driven by shared measurement error.   ║");
  log("  ║  Do NOT add to baseline.                                               ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  M_halo: LIKELY ARTIFACT                                               ║");
  log("  ║  Signal is probably driven by circularity.                              ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  decomposition: {
    VflatAfterBaseline: { r: +vflatAfterBL.toFixed(3) },
    RlastAfterBaseline: { r: +rlastAfterBL.toFixed(3) },
    MhAfterBaselineVflat: { r: +mhAfterBLVf.toFixed(3) },
  },
  jackknife: { mean: +jackMean.toFixed(3), sd: +jackSD.toFixed(3), signFlips, nInfluential: influential.length },
  bootstrap: {
    baseline: { lo: +ci_base_lo.toFixed(3), hi: +ci_base_hi.toFixed(3), excludesZero: bootstrapExcludes },
    baselineVflat: { lo: +ci_vf_lo.toFixed(3), hi: +ci_vf_hi.toFixed(3), excludesZero: bootstrapVfExcludes },
  },
  circularityScore: { passed: nCircPass, total: 5 },
  verdict: nCircPass >= 4 ? 'CONFIRMED' : nCircPass >= 3 ? 'CONFIRMED_WITH_CAVEAT' : nCircPass >= 2 ? 'PARTIAL' : 'LIKELY_ARTIFACT',
};

fs.writeFileSync(path.join(__dirname, '../public/phase33b-circularity.json'), JSON.stringify(output, null, 2));
log("  Results saved to public/phase33b-circularity.json");
