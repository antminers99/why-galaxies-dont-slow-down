#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "37.0.0";
function log(msg) { console.log(msg); }

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

  const Rdisk = sp.Rdisk || sp.Reff || 3.0;
  const Mbar_sun = sp.MHI * 1.33 + L * 1e9 * 0.5;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const innermost = sorted.slice(0, Math.max(3, Math.floor(sorted.length * 0.2)));
  const gObsInner = innermost.map(p => Math.pow(10, p.log_g_obs));
  const rInner = innermost.map(p => p.r);

  const rho0_obs = (() => {
    if (innermost.length < 2) return 1e-3;
    const gMax = Math.max(...gObsInner);
    const rMin = Math.min(...rInner);
    const rKpc = Math.max(rMin, 0.5);
    const rCm = rKpc * 3.086e21;
    const gCgs = gMax * 1e-10;
    const rho = 3 * gCgs / (4 * Math.PI * 6.674e-8 * rCm);
    return rho > 0 ? rho : 1e-30;
  })();
  const logRho0 = Math.log10(rho0_obs);

  const Sigma0_bar = (() => {
    const rKpc = Rdisk;
    const rPc = rKpc * 1e3;
    const areaPc2 = Math.PI * rPc * rPc;
    const sigmaM = Mbar_sun / areaPc2;
    return sigmaM > 0 ? sigmaM : 1e-3;
  })();
  const logSigma0 = Math.log10(Sigma0_bar);

  const gObsCentral = Math.max(...gObsInner);
  const logGobsCentral = Math.log10(gObsCentral);

  const rho0_pseudo = (() => {
    if (innermost.length < 3) return 1e-3;
    const innerSlope = (() => {
      const xv = innermost.map(p => Math.log10(p.r));
      const yv = innermost.map(p => p.log_g_obs);
      const mx = xv.reduce((s,v)=>s+v,0)/xv.length;
      const my = yv.reduce((s,v)=>s+v,0)/yv.length;
      let sxy = 0, sxx = 0;
      for (let i = 0; i < xv.length; i++) { sxy += (xv[i]-mx)*(yv[i]-my); sxx += (xv[i]-mx)**2; }
      return sxx > 0 ? sxy / sxx : 0;
    })();
    const gAtR1 = gObsInner[0];
    const r1kpc = Math.max(rInner[0], 0.5);
    const V1 = Math.sqrt(gAtR1 * r1kpc * 3.086e19) / 1e3;
    const rho = V1 * V1 / (4 * Math.PI * 6.674e-11 * (r1kpc * 3.086e19) ** 2) * 1.989e30;
    return rho > 0 ? rho : 1e-30;
  })();
  const logRho0Pseudo = Math.log10(rho0_pseudo > 0 ? rho0_pseudo : 1e-30);

  const MdynInner = (() => {
    if (innermost.length < 2) return 1e8;
    const rOut = Math.max(...rInner);
    const gOut = gObsInner[gObsInner.length - 1];
    const VkmsOut = Math.sqrt(gOut * rOut * 3.086e19) / 1e3;
    const Mdyn = (VkmsOut * 1e3) ** 2 * (rOut * 1e3 * 3.086e16) / (6.674e-11 * 1.989e30);
    return Mdyn > 0 ? Mdyn : 1e8;
  })();
  const logMdynInner = Math.log10(MdynInner);

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, logVflat, rcWiggliness, envCode,
    logRho0, logSigma0, logGobsCentral, logRho0Pseudo, logMdynInner,
    Rdisk, rMax,
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
const rho0Arr = galaxyData.map(g => g.logRho0);
const sigma0Arr = galaxyData.map(g => g.logSigma0);
const gCentralArr = galaxyData.map(g => g.logGobsCentral);
const rho0PseudoArr = galaxyData.map(g => g.logRho0Pseudo);
const mdynInnerArr = galaxyData.map(g => g.logMdynInner);

const BASELINE = [mhiArr, wigArr, envArr];

log("");
log("=".repeat(80));
log("  PHASE 37: rho_0 (Central Density)");
log("  FAST PROTOCOL — Baseline: MHI + wig + env");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  CENTRAL DENSITY PROXIES:");
log("  - rho0_obs: central mass density from innermost g_obs");
log("  - Sigma0_bar: baryonic surface density (Mbar/pi*Rdisk^2)");
log("  - gobs_central: peak observed acceleration in inner region");
log("  - rho0_pseudo: pseudo isothermal density from inner V(r)");
log("  - Mdyn_inner: dynamical mass within inner RC region");
log("");

log("  CIRCULARITY:");
log("  - Sigma0_bar: EXTERNAL (photometric) — LOW circularity");
log("  - rho0_obs: from inner g_obs — HIGH circularity (g_obs → a0)");
log("  - gobs_central: DIRECTLY from RC — HIGHEST circularity");
log("  - rho0_pseudo, Mdyn_inner: from inner velocities — HIGH");
log("");

const candidates = [
  { name: 'rho0_obs', arr: rho0Arr, feat: 'logRho0', circ: 'HIGH' },
  { name: 'Sigma0_bar', arr: sigma0Arr, feat: 'logSigma0', circ: 'LOW' },
  { name: 'gobs_central', arr: gCentralArr, feat: 'logGobsCentral', circ: 'HIGHEST' },
  { name: 'rho0_pseudo', arr: rho0PseudoArr, feat: 'logRho0Pseudo', circ: 'HIGH' },
  { name: 'Mdyn_inner', arr: mdynInnerArr, feat: 'logMdynInner', circ: 'HIGH' },
];

log("  COLLINEARITY:");
for (const c of candidates) {
  const rMhi = corrWith(c.arr, mhiArr);
  const rWig = corrWith(c.arr, wigArr);
  const rVf = corrWith(c.arr, vflatArr);
  log("  " + c.name.padEnd(14) + " vs MHI:" + rMhi.toFixed(2).padStart(6) + "  wig:" + rWig.toFixed(2).padStart(6) + "  Vf:" + rVf.toFixed(2).padStart(6) + "  circ:" + c.circ);
}
log("");

log("=".repeat(80));
log("  STEP 1: FAST SCREEN");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'envCode': envArr,
  'logVflat': vflatArr,
  'logRho0': rho0Arr, 'logSigma0': sigma0Arr,
  'logGobsCentral': gCentralArr, 'logRho0Pseudo': rho0PseudoArr,
  'logMdynInner': mdynInnerArr,
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
  if (featureNames.length === 0) return { gc: 0 };
  const modelRms = Math.sqrt(totalSS_model / totalN);
  return { gc: gap > 0 ? (m0rms - modelRms) / gap * 100 : 0 };
}

const baselineGC = looGapClosed(['logMHI', 'rcWiggliness', 'envCode']);
log("  BASELINE: " + baselineGC.gc.toFixed(1) + "% gap closure");
log("");

log("  ┌──────────────────────────────────────────────────────────────────────────────────────────┐");
log("  │  Proxy         Circ     Raw r    Raw t   |base r   |base t   LOO+base  Adds?  Verdict │");
log("  ├──────────────────────────────────────────────────────────────────────────────────────────┤");

const results = [];
for (const cand of candidates) {
  const rRaw = corrWith(cand.arr, deltaA0);
  const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);

  const rx = residualize(cand.arr, BASELINE);
  const ry = residualize(deltaA0, BASELINE);
  const rBase = corrWith(rx, ry);
  const tBase = Math.abs(rBase) * Math.sqrt(N-5) / Math.sqrt(1-rBase**2+1e-10);

  const gc = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', cand.feat]);
  const adds = gc.gc > baselineGC.gc + 1;

  const failFast = tRaw < 1.0 && tBase < 1.0;
  const promising = tBase >= 1.65 || adds;
  const verdict = failFast ? "FAIL-FAST" : promising ? "PROMISING" : "WEAK";

  results.push({ ...cand, rRaw, tRaw, rBase, tBase, gc: gc.gc, adds, verdict });

  log("  │  " + cand.name.padEnd(13) + cand.circ.padEnd(8) +
    rRaw.toFixed(3).padStart(7) + tRaw.toFixed(1).padStart(7) +
    rBase.toFixed(3).padStart(9) + tBase.toFixed(1).padStart(9) +
    (gc.gc.toFixed(1) + "%").padStart(10) +
    (adds ? "  YES " : "  NO  ") + verdict.padStart(10) + " │");
}
log("  └──────────────────────────────────────────────────────────────────────────────────────────┘");
log("");

const promising = results.filter(r => r.verdict === "PROMISING");

if (promising.length === 0) {
  log("  ══════════════════════════════════════════════════════════");
  log("  NO central density proxy passes fast screen.");
  log("  ALL FAIL or WEAK. Skipping deep tests.");
  log("  ══════════════════════════════════════════════════════════");
} else {
  log("  PROMISING: " + promising.map(p => p.name).join(", "));
  log("  Proceeding to deep tests...");
  log("");

  for (const cand of promising) {
    log("=".repeat(80));
    log("  DEEP TEST: " + cand.name + " (circularity: " + cand.circ + ")");
    log("=".repeat(80));
    log("");

    log("  Confounder stripping:");
    for (const [label, ctrls] of [
      ['raw', []], ['|MHI', [mhiArr]], ['|MHI+wig', [mhiArr, wigArr]],
      ['|baseline', BASELINE], ['|baseline+Vflat', [...BASELINE, vflatArr]],
    ]) {
      const rx2 = residualize(cand.arr, ctrls);
      const ry2 = residualize(deltaA0, ctrls);
      const r = corrWith(rx2, ry2);
      const df = N - 2 - ctrls.length;
      const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
      log("    " + label.padEnd(20) + " r=" + r.toFixed(3) + " t=" + t.toFixed(1) +
        (t >= 2.0 ? " SURVIVES" : t >= 1.65 ? " MARGINAL" : " fails"));
    }
    log("");

    log("  Permutation (10000):");
    const NPERMS = 10000;
    for (const [label, ctrls] of [['|baseline', BASELINE], ['|baseline+Vflat', [...BASELINE, vflatArr]]]) {
      const rx2 = residualize(cand.arr, ctrls);
      const ry2 = residualize(deltaA0, ctrls);
      const obsR = Math.abs(corrWith(rx2, ry2));
      let cnt = 0;
      const sh = [...ry2];
      for (let p = 0; p < NPERMS; p++) {
        for (let i = sh.length - 1; i > 0; i--) {
          const j = Math.floor(Math.random() * (i + 1));
          [sh[i], sh[j]] = [sh[j], sh[i]];
        }
        if (Math.abs(corrWith(rx2, sh)) >= obsR) cnt++;
      }
      log("    " + label.padEnd(20) + " |r|=" + obsR.toFixed(3) + " p=" + (cnt/NPERMS).toFixed(4) +
        (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "**" : ""));
    }
    log("");

    log("  Jackknife:");
    const jackR = [];
    for (let i = 0; i < N; i++) {
      const idx = [...Array(N).keys()].filter(j => j !== i);
      const subX = idx.map(j => cand.arr[j]);
      const subY = idx.map(j => deltaA0[j]);
      const subBL = BASELINE.map(arr => idx.map(j => arr[j]));
      const rx2 = residualize(subX, subBL);
      const ry2 = residualize(subY, subBL);
      jackR.push(corrWith(rx2, ry2));
    }
    const jackMean = jackR.reduce((s,v)=>s+v,0)/N;
    const jackSD = Math.sqrt(jackR.map(v => (v - jackMean)**2).reduce((s,v)=>s+v,0)/(N-1));
    const signFlips = jackR.filter(r => Math.sign(r) !== Math.sign(cand.rBase)).length;
    log("    Mean=" + jackMean.toFixed(3) + " SD=" + jackSD.toFixed(3) + " flips=" + signFlips + "/" + N);
    log("");

    log("  Bootstrap CI (2000):");
    const NBOOT = 2000;
    const bootR = [];
    for (let b = 0; b < NBOOT; b++) {
      const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
      const bX = idx.map(i => cand.arr[i]);
      const bY = idx.map(i => deltaA0[i]);
      const bBL = BASELINE.map(arr => idx.map(i => arr[i]));
      const rx2 = residualize(bX, bBL);
      const ry2 = residualize(bY, bBL);
      bootR.push(corrWith(rx2, ry2));
    }
    bootR.sort((a,b) => a - b);
    const lo = bootR[Math.floor(NBOOT * 0.025)];
    const hi = bootR[Math.floor(NBOOT * 0.975)];
    log("    95% CI: [" + lo.toFixed(3) + ", " + hi.toFixed(3) + "]" +
      (lo > 0 || hi < 0 ? " excludes zero" : " crosses zero"));
    log("");
  }
}

log("=".repeat(80));
log("  PHASE 37 — FAST PROTOCOL SUMMARY");
log("=".repeat(80));
log("");

const best = results.reduce((a,b) => Math.abs(b.tBase) > Math.abs(a.tBase) ? b : a, results[0]);
const nFP = [best.tRaw >= 2, best.tBase >= 1.65, best.adds, best.gc > baselineGC.gc + 1].filter(Boolean).length;

log("  ┌──────────────────────────────────────────────────────────────────────┐");
log("  │  Test                        Result              Pass/Fail        │");
log("  ├──────────────────────────────────────────────────────────────────────┤");
log("  │  Raw correlation             r=" + best.rRaw.toFixed(3) + " t=" + best.tRaw.toFixed(1) + "          " + (best.tRaw >= 2 ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  After baseline control      r=" + best.rBase.toFixed(3) + " t=" + best.tBase.toFixed(1) + "          " + (best.tBase >= 1.65 ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  LOO > baseline              " + best.gc.toFixed(1) + "% vs " + baselineGC.gc.toFixed(1) + "%" + "           " + (best.adds ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  Adds above baseline         " + (best.adds ? "YES" : "NO").padEnd(24) + (best.adds ? "PASS" : "FAIL").padEnd(12) + "│");
log("  ├──────────────────────────────────────────────────────────────────────┤");

const verdict = nFP >= 3 ? "CONFIRMED" : nFP >= 2 ? "PARTIAL" : "FAIL";
log("  │  FINAL VERDICT: " + verdict.padEnd(12) + " (best: " + best.name + ", circ: " + best.circ + ")" + " ".repeat(10) + "│");
log("  └──────────────────────────────────────────────────────────────────────┘");
log("");

log("  DARK HALO DOOR — SCORECARD:");
log("    M_halo/Mdyn:   PARTIAL (= V_flat, circularity 2/5)");
log("    concentration: FAIL");
log("    V_flat:        PASS but CIRCULAR");
log("    R_s:           PARTIAL (Rs_Vrise)");
log("    rho_0:         " + verdict + " (best=" + best.name + ")");
log("    Next: f_DM (dark matter fraction)");
log("");
log("  Baseline UNCHANGED: MHI + rcWiggliness + envCode = " + baselineGC.gc.toFixed(1) + "%");
log("");
log("=".repeat(80));

const output = {
  version: VERSION, timestamp: new Date().toISOString(), nGalaxies: N,
  baseline: { gc: +baselineGC.gc.toFixed(1) },
  proxies: results.map(r => ({
    name: r.name, circ: r.circ, rawR: +r.rRaw.toFixed(3), rawT: +r.tRaw.toFixed(1),
    baseR: +r.rBase.toFixed(3), baseT: +r.tBase.toFixed(1),
    looGC: +r.gc.toFixed(1), adds: r.adds, verdict: r.verdict,
  })),
  bestProxy: best.name, finalVerdict: verdict,
};

fs.writeFileSync(path.join(__dirname, '../public/phase37-central-density.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase37-central-density.json");
