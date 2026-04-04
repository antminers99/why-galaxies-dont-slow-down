#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "35.0.0";
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
  const Vflat = sp.Vflat || 100;
  const logVflat = Math.log10(Vflat);

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
    logA0: fit.logA0, a0: fit.a0,
    logMHI, logVflat, rcWiggliness, envCode, Vflat, rMax,
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

const BASELINE = [mhiArr, wigArr, envArr];

log("");
log("=".repeat(80));
log("  PHASE 35: V_flat (Flat Rotation Velocity)");
log("  FAST PROTOCOL — Baseline: MHI + wig + env");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  CONTEXT:");
log("  Phase 33b showed: the M_halo suppressor signal = V_flat.");
log("  V_flat|baseline: r=0.472, t=3.8 (strongest partial in project).");
log("  Now testing V_flat DIRECTLY as a halo-door variable.");
log("");
log("  CRITICAL CIRCULARITY FLAG:");
log("  V_flat comes from the SAME rotation curve that determines a0.");
log("  Both are derived from outer velocities.");
log("  This is the MOST circular variable in the entire project.");
log("  Even if it passes all tests, the circularity caveat is SEVERE.");
log("");

log("  COLLINEARITY:");
const rVfMhi = corrWith(vflatArr, mhiArr);
const rVfWig = corrWith(vflatArr, wigArr);
const rVfEnv = corrWith(vflatArr, envArr);
log("  logVflat vs logMHI: r=" + rVfMhi.toFixed(3));
log("  logVflat vs wig:    r=" + rVfWig.toFixed(3));
log("  logVflat vs env:    r=" + rVfEnv.toFixed(3));
log("");

log("=".repeat(80));
log("  STEP 1: FAST SCREEN");
log("=".repeat(80));
log("");

const rRaw = corrWith(vflatArr, deltaA0);
const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);

const residX_mhi = residualize(vflatArr, [mhiArr]);
const residY_mhi = residualize(deltaA0, [mhiArr]);
const rMHI = corrWith(residX_mhi, residY_mhi);
const tMHI = Math.abs(rMHI) * Math.sqrt(N-3) / Math.sqrt(1-rMHI**2+1e-10);

const residX_mhiw = residualize(vflatArr, [mhiArr, wigArr]);
const residY_mhiw = residualize(deltaA0, [mhiArr, wigArr]);
const rMHIW = corrWith(residX_mhiw, residY_mhiw);
const tMHIW = Math.abs(rMHIW) * Math.sqrt(N-4) / Math.sqrt(1-rMHIW**2+1e-10);

const residX_bl = residualize(vflatArr, BASELINE);
const residY_bl = residualize(deltaA0, BASELINE);
const rBL = corrWith(residX_bl, residY_bl);
const tBL = Math.abs(rBL) * Math.sqrt(N-5) / Math.sqrt(1-rBL**2+1e-10);

log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Test                    r         t       Status              │");
log("  ├───────────────────────────────────────────────────────────────────┤");
log("  │  Raw                  " + rRaw.toFixed(3).padStart(7) + tRaw.toFixed(1).padStart(8) + "       " + (tRaw >= 2 ? "PASS" : "FAIL").padEnd(10) + "│");
log("  │  After MHI            " + rMHI.toFixed(3).padStart(7) + tMHI.toFixed(1).padStart(8) + "       " + (tMHI >= 1.65 ? "PASS" : "FAIL").padEnd(10) + "│");
log("  │  After MHI+wig        " + rMHIW.toFixed(3).padStart(7) + tMHIW.toFixed(1).padStart(8) + "       " + (tMHIW >= 1.65 ? "PASS" : "FAIL").padEnd(10) + "│");
log("  │  After baseline       " + rBL.toFixed(3).padStart(7) + tBL.toFixed(1).padStart(8) + "       " + (tBL >= 1.65 ? "PASS" : "FAIL").padEnd(10) + "│");
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'envCode': envArr, 'logVflat': vflatArr,
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

const baselineGC = looGapClosed(['logMHI', 'rcWiggliness', 'envCode']);
const vflatGC = looGapClosed(['logVflat']);
const blVfGC = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', 'logVflat']);

log("  LOO CROSS-VALIDATION:");
log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Model                         LOO gap-closed                  │");
log("  ├───────────────────────────────────────────────────────────────────┤");
log("  │  Vflat alone                   " + (vflatGC.gc.toFixed(1) + "%").padStart(8) + " ".repeat(24) + "│");
log("  │  baseline (MHI+wig+env)        " + (baselineGC.gc.toFixed(1) + "%").padStart(8) + " ".repeat(24) + "│");
log("  │  baseline + Vflat              " + (blVfGC.gc.toFixed(1) + "%").padStart(8) + " ".repeat(24) + "│");
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

const adds = blVfGC.gc > baselineGC.gc + 1;
log("  V_flat adds beyond baseline? " + (adds ? "YES (+" + (blVfGC.gc - baselineGC.gc).toFixed(1) + "%)" : "NO"));
log("");

const fastPass = tBL >= 1.65 && adds;
const failFast = tRaw < 1.0 && tBL < 1.0;

if (failFast) {
  log("  ══════════════════════════════════════════════════════════");
  log("  FAIL-FAST: V_flat does not pass basic screen.");
  log("  ══════════════════════════════════════════════════════════");
} else if (fastPass) {
  log("  PROMISING — Proceeding to deep tests...");
  log("");

  log("=".repeat(80));
  log("  STEP 2: DEEP TESTS");
  log("=".repeat(80));
  log("");

  log("  Permutation tests (10000 perms):");
  const NPERMS = 10000;
  for (const [label, ctrls] of [
    ['raw', []], ['|MHI', [mhiArr]], ['|MHI+wig', [mhiArr, wigArr]], ['|baseline', BASELINE],
  ]) {
    const rx = residualize(vflatArr, ctrls);
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
    log("    logVflat " + label.padEnd(14) + " |r|=" + obsR.toFixed(3) + " p=" + (cnt/NPERMS).toFixed(4) +
      (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "**" : ""));
  }
  log("");

  log("  Jackknife (N=" + N + "):");
  const jackR = [];
  for (let i = 0; i < N; i++) {
    const idx = [...Array(N).keys()].filter(j => j !== i);
    const subVf = idx.map(j => vflatArr[j]);
    const subDa = idx.map(j => deltaA0[j]);
    const subBL = BASELINE.map(arr => idx.map(j => arr[j]));
    const rx = residualize(subVf, subBL);
    const ry = residualize(subDa, subBL);
    jackR.push(corrWith(rx, ry));
  }
  const jackMean = jackR.reduce((s,v)=>s+v,0)/N;
  const jackSD = Math.sqrt(jackR.map(v => (v - jackMean)**2).reduce((s,v)=>s+v,0)/(N-1));
  const signFlips = jackR.filter(r => Math.sign(r) !== Math.sign(rBL)).length;
  log("    Mean r=" + jackMean.toFixed(3) + " SD=" + jackSD.toFixed(3));
  log("    Range: [" + Math.min(...jackR).toFixed(3) + ", " + Math.max(...jackR).toFixed(3) + "]");
  log("    Sign flips: " + signFlips + "/" + N);
  log("");

  log("  Bootstrap CI (2000 resamples):");
  const NBOOT = 2000;
  const bootR = [];
  for (let b = 0; b < NBOOT; b++) {
    const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
    const bVf = idx.map(i => vflatArr[i]);
    const bDa = idx.map(i => deltaA0[i]);
    const bBL = BASELINE.map(arr => idx.map(i => arr[i]));
    const rx = residualize(bVf, bBL);
    const ry = residualize(bDa, bBL);
    bootR.push(corrWith(rx, ry));
  }
  bootR.sort((a,b) => a - b);
  const ciLo = bootR[Math.floor(NBOOT * 0.025)];
  const ciHi = bootR[Math.floor(NBOOT * 0.975)];
  log("    95% CI: [" + ciLo.toFixed(3) + ", " + ciHi.toFixed(3) + "]" +
    (ciLo > 0 ? " → EXCLUDES ZERO" : " → crosses zero"));
  log("");

  log("  Stratification:");
  const medianMHI = [...mhiArr].sort((a,b) => a - b)[Math.floor(N/2)];
  for (const [label, idx] of [
    ['Low MHI', galaxyData.map((_,i) => i).filter(i => mhiArr[i] < medianMHI)],
    ['High MHI', galaxyData.map((_,i) => i).filter(i => mhiArr[i] >= medianMHI)],
    ['Field', galaxyData.map((_,i) => i).filter(i => envArr[i] === 0)],
    ['Group/cluster', galaxyData.map((_,i) => i).filter(i => envArr[i] > 0)],
  ]) {
    if (idx.length < 8) { log("    " + label + ": too few"); continue; }
    const subVf = idx.map(i => vflatArr[i]);
    const subDa = idx.map(i => deltaA0[i]);
    const subMHI = idx.map(i => mhiArr[i]);
    const subWig = idx.map(i => wigArr[i]);
    const rx = residualize(subVf, [subMHI, subWig]);
    const ry = residualize(subDa, [subMHI, subWig]);
    const r = corrWith(rx, ry);
    const df = idx.length - 4;
    const t = df > 0 ? Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10) : 0;
    log("    " + (label + " (N=" + idx.length + ")").padEnd(24) + " r=" + r.toFixed(3) + " t=" + t.toFixed(1) +
      (r > 0 ? " same sign" : " OPPOSITE"));
  }
  log("");
} else {
  log("  WEAK — V_flat has signal but doesn't clearly pass all fast criteria.");
  log("  Proceeding to deep tests anyway (V_flat is central to project)...");
  log("");

  log("=".repeat(80));
  log("  STEP 2: DEEP TESTS (for completeness)");
  log("=".repeat(80));
  log("");

  log("  Permutation tests (10000 perms):");
  const NPERMS = 10000;
  for (const [label, ctrls] of [
    ['raw', []], ['|MHI', [mhiArr]], ['|MHI+wig', [mhiArr, wigArr]], ['|baseline', BASELINE],
  ]) {
    const rx = residualize(vflatArr, ctrls);
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
    log("    logVflat " + label.padEnd(14) + " |r|=" + obsR.toFixed(3) + " p=" + (cnt/NPERMS).toFixed(4) +
      (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "**" : ""));
  }
  log("");

  log("  Jackknife (N=" + N + "):");
  const jackR = [];
  for (let i = 0; i < N; i++) {
    const idx = [...Array(N).keys()].filter(j => j !== i);
    const subVf = idx.map(j => vflatArr[j]);
    const subDa = idx.map(j => deltaA0[j]);
    const subBL = BASELINE.map(arr => idx.map(j => arr[j]));
    const rx = residualize(subVf, subBL);
    const ry = residualize(subDa, subBL);
    jackR.push(corrWith(rx, ry));
  }
  const jackMean = jackR.reduce((s,v)=>s+v,0)/N;
  const jackSD = Math.sqrt(jackR.map(v => (v - jackMean)**2).reduce((s,v)=>s+v,0)/(N-1));
  const signFlips = jackR.filter(r => Math.sign(r) !== Math.sign(rBL)).length;
  log("    Mean r=" + jackMean.toFixed(3) + " SD=" + jackSD.toFixed(3));
  log("    Range: [" + Math.min(...jackR).toFixed(3) + ", " + Math.max(...jackR).toFixed(3) + "]");
  log("    Sign flips: " + signFlips + "/" + N);
  log("");

  log("  Bootstrap CI (2000 resamples):");
  const NBOOT = 2000;
  const bootR = [];
  for (let b = 0; b < NBOOT; b++) {
    const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
    const bVf = idx.map(i => vflatArr[i]);
    const bDa = idx.map(i => deltaA0[i]);
    const bBL = BASELINE.map(arr => idx.map(i => arr[i]));
    const rx = residualize(bVf, bBL);
    const ry = residualize(bDa, bBL);
    bootR.push(corrWith(rx, ry));
  }
  bootR.sort((a,b) => a - b);
  const ciLo = bootR[Math.floor(NBOOT * 0.025)];
  const ciHi = bootR[Math.floor(NBOOT * 0.975)];
  log("    95% CI: [" + ciLo.toFixed(3) + ", " + ciHi.toFixed(3) + "]" +
    (ciLo > 0 ? " → EXCLUDES ZERO" : " → crosses zero"));
  log("");

  log("  Stratification:");
  const medianMHI = [...mhiArr].sort((a,b) => a - b)[Math.floor(N/2)];
  for (const [label, idx] of [
    ['Low MHI', galaxyData.map((_,i) => i).filter(i => mhiArr[i] < medianMHI)],
    ['High MHI', galaxyData.map((_,i) => i).filter(i => mhiArr[i] >= medianMHI)],
    ['Field', galaxyData.map((_,i) => i).filter(i => envArr[i] === 0)],
    ['Group/cluster', galaxyData.map((_,i) => i).filter(i => envArr[i] > 0)],
  ]) {
    if (idx.length < 8) { log("    " + label + ": too few"); continue; }
    const subVf = idx.map(i => vflatArr[i]);
    const subDa = idx.map(i => deltaA0[i]);
    const subMHI = idx.map(i => mhiArr[i]);
    const subWig = idx.map(i => wigArr[i]);
    const rx = residualize(subVf, [subMHI, subWig]);
    const ry = residualize(subDa, [subMHI, subWig]);
    const r = corrWith(rx, ry);
    const df = idx.length - 4;
    const t = df > 0 ? Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10) : 0;
    log("    " + (label + " (N=" + idx.length + ")").padEnd(24) + " r=" + r.toFixed(3) + " t=" + t.toFixed(1) +
      (r > 0 ? " same sign" : " OPPOSITE"));
  }
  log("");
}

log("=".repeat(80));
log("  PHASE 35 — FAST PROTOCOL SUMMARY");
log("=".repeat(80));
log("");

log("  ┌──────────────────────────────────────────────────────────────────────┐");
log("  │  Test                        Result              Pass/Fail        │");
log("  ├──────────────────────────────────────────────────────────────────────┤");
log("  │  Raw correlation             r=" + rRaw.toFixed(3) + " t=" + tRaw.toFixed(1) + "          " + (tRaw >= 2 ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  After baseline control      r=" + rBL.toFixed(3) + " t=" + tBL.toFixed(1) + "          " + (tBL >= 1.65 ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  LOO > baseline              " + blVfGC.gc.toFixed(1) + "% vs " + baselineGC.gc.toFixed(1) + "%" + "           " + (adds ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  Adds above baseline         " + (adds ? "YES (+"+((blVfGC.gc-baselineGC.gc).toFixed(1))+"%)" : "NO").padEnd(24) + (adds ? "PASS" : "FAIL").padEnd(12) + "│");
log("  ├──────────────────────────────────────────────────────────────────────┤");

const nFastPass = [tRaw >= 2, tBL >= 1.65, adds, blVfGC.gc > baselineGC.gc + 1].filter(Boolean).length;

let verdict;
if (nFastPass >= 3) {
  verdict = "CONFIRMED — but CIRCULARITY CAVEAT remains SEVERE";
} else if (nFastPass >= 2) {
  verdict = "PARTIAL — signal present but circularity unresolved";
} else {
  verdict = "FAIL";
}

log("  │  FINAL VERDICT: " + verdict.padEnd(50) + "│");
log("  └──────────────────────────────────────────────────────────────────────┘");
log("");

log("  ⚠  MANDATORY CIRCULARITY NOTE:");
log("  V_flat enters BOTH sides of this analysis:");
log("    - a0 is fitted from rotation curve (includes V_flat)");
log("    - V_flat IS the outer rotation velocity");
log("  Even if statistically strong, this correlation may be");
log("  MECHANICALLY inevitable rather than physically meaningful.");
log("  This caveat applies equally to M_halo (which IS V_flat).");
log("");

log("  DARK HALO DOOR — SCORECARD:");
log("    M_halo:        PARTIAL (= V_flat, circularity 2/5)");
log("    concentration: FAIL");
log("    V_flat:        " + (nFastPass >= 3 ? "PASS but CIRCULAR" : nFastPass >= 2 ? "PARTIAL" : "FAIL"));
log("    Next: R_s (scale radius)");
log("");
log("  Baseline UNCHANGED: MHI + rcWiggliness + envCode = 14.9% gap closure");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  collinearity: { VfMHI: +rVfMhi.toFixed(3), VfWig: +rVfWig.toFixed(3), VfEnv: +rVfEnv.toFixed(3) },
  fastScreen: {
    rawR: +rRaw.toFixed(3), rawT: +tRaw.toFixed(1),
    afterMHI: { r: +rMHI.toFixed(3), t: +tMHI.toFixed(1) },
    afterMHIWig: { r: +rMHIW.toFixed(3), t: +tMHIW.toFixed(1) },
    afterBaseline: { r: +rBL.toFixed(3), t: +tBL.toFixed(1) },
  },
  loo: { vflatAlone: +vflatGC.gc.toFixed(1), baseline: +baselineGC.gc.toFixed(1), baselineVflat: +blVfGC.gc.toFixed(1) },
  addsAboveBaseline: adds,
  nFastPass,
  verdict: nFastPass >= 3 ? 'CONFIRMED_WITH_CIRCULARITY' : nFastPass >= 2 ? 'PARTIAL' : 'FAIL',
  circularitySeverity: 'SEVERE',
};

fs.writeFileSync(path.join(__dirname, '../public/phase35-vflat.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase35-vflat.json");
