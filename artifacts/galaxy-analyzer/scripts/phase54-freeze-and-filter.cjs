#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "54.0.0";
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

function looGapClosed(featArrays, allLogA0, galaxyData, N) {
  let totalSS_model = 0, totalSS_m0 = 0, totalSS_free = 0, totalN = 0;
  for (let i = 0; i < N; i++) {
    const trainIdx = [...Array(N).keys()].filter(j => j !== i);
    const trainY = trainIdx.map(j => allLogA0[j]);
    const mu = trainY.reduce((s,v)=>s+v,0)/trainY.length;
    totalSS_m0 += predictSS(Math.pow(10, mu), galaxyData[i].pts);
    totalSS_free += predictSS(galaxyData[i].a0, galaxyData[i].pts);
    totalN += galaxyData[i].pts.length;
    if (featArrays.length > 0) {
      const trainX = trainIdx.map(j => featArrays.map(arr => arr[j]));
      const reg = linReg(trainX, trainY);
      const testX = featArrays.map(arr => arr[i]);
      let pred = reg.intercept + reg.coefs.reduce((s, c, j) => s + c * testX[j], 0);
      pred = Math.max(2.5, Math.min(4.5, pred));
      totalSS_model += predictSS(Math.pow(10, pred), galaxyData[i].pts);
    } else {
      totalSS_model = totalSS_m0;
    }
  }
  const m0rms = Math.sqrt(totalSS_m0 / totalN);
  const m6rms = Math.sqrt(totalSS_free / totalN);
  const gap = m0rms - m6rms;
  const modelRms = Math.sqrt(totalSS_model / totalN);
  return { gc: gap > 0 ? (m0rms - modelRms) / gap * 100 : 0 };
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
  const Sigma0_bar = Mbar_sun / (Math.PI * (Rdisk * 1e3) ** 2);
  const logSigma0 = Math.log10(Sigma0_bar > 0 ? Sigma0_bar : 1e-3);

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const rarResid = sorted.map(p => {
    const gb = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gb, fit.a0);
    return p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10);
  });

  let currentSign = Math.sign(rarResid[0]);
  let runLen = 1;
  const residRun = [];
  for (let i = 1; i < rarResid.length; i++) {
    const s = Math.sign(rarResid[i]);
    if (s === currentSign) { runLen++; }
    else { residRun.push(runLen); runLen = 1; currentSign = s; }
  }
  residRun.push(runLen);
  const meanRunLen = residRun.reduce((s,v)=>s+v,0) / residRun.length;
  const logMeanRun = Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1);

  const nPts = pts.length;
  const outerFrac = outerPts.length / pts.length;

  const outerResid = sorted.slice(Math.floor(sorted.length * 0.5)).map(p => {
    const gb = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gb, fit.a0);
    return p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10);
  });
  let outerZC = 0;
  for (let i = 1; i < outerResid.length; i++) {
    if (Math.sign(outerResid[i]) !== Math.sign(outerResid[i-1])) outerZC++;
  }
  const outerOsc = outerResid.length > 1 ? outerZC / (outerResid.length - 1) : 0;

  const tBstRd = (() => {
    let maxBoost = -Infinity, maxR = rMax * 0.5;
    for (const p of sorted) {
      const gb = Math.pow(10, p.log_g_bar);
      const predUniv = mcgaughRAR(gb, 3633);
      const boost = p.log_g_obs - Math.log10(predUniv > 0 ? predUniv : 1e-10);
      if (boost > maxBoost) { maxBoost = boost; maxR = p.r; }
    }
    return maxR / Rdisk;
  })();
  const logTBstRd = Math.log10(tBstRd > 0.01 ? tBstRd : 0.01);

  const innerRise = (() => {
    if (sorted.length < 5) return 0;
    const n20 = Math.max(3, Math.floor(sorted.length * 0.2));
    const early = sorted.slice(0, n20).map(p => p.log_g_obs);
    const earlyMean = early.reduce((s,v)=>s+v,0)/early.length;
    const mid = sorted.slice(Math.floor(sorted.length * 0.3), Math.floor(sorted.length * 0.5));
    if (mid.length === 0) return 0;
    const midMean = mid.map(p => p.log_g_obs).reduce((s,v)=>s+v,0)/mid.length;
    return midMean - earlyMean;
  })();

  const RsVrise = (() => {
    let peakV = -Infinity, peakR = 1;
    for (const p of sorted) {
      if (p.log_g_obs > peakV) { peakV = p.log_g_obs; peakR = p.r; }
    }
    return peakR / Rdisk;
  })();
  const logRsVrise = Math.log10(RsVrise > 0.01 ? RsVrise : 0.01);

  const maxJump = (() => {
    let mj = 0;
    for (let i = 1; i < sorted.length; i++) {
      const dv = Math.abs(sorted[i].log_g_obs - sorted[i-1].log_g_obs);
      const dr = Math.abs(sorted[i].r - sorted[i-1].r);
      if (dr > 0) mj = Math.max(mj, dv / dr);
    }
    return mj;
  })();

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, logVflat, rcWiggliness, envCode, logSigma0,
    logMeanRun, nPts, outerFrac, outerOsc,
    logTBstRd, innerRise, logRsVrise, maxJump,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const sigma0Arr = galaxyData.map(g => g.logSigma0);
const meanRunArr = galaxyData.map(g => g.logMeanRun);
const vflatArr = galaxyData.map(g => g.logVflat);

const CONS_BL = [mhiArr, wigArr, envArr, sigma0Arr];
const WORK_BL = [mhiArr, wigArr, envArr, sigma0Arr, meanRunArr];

log("");
log("╔══════════════════════════════════════════════════════════════════════════════════╗");
log("║  PHASE 54: FREEZE BASELINES + FILTER PARTIAL VARIABLES                        ║");
log("║  Version " + VERSION + ", " + N + " galaxies                                           ║");
log("╚══════════════════════════════════════════════════════════════════════════════════╝");
log("");

log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  PART A: BASELINE FREEZE");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

const gcCons = looGapClosed(CONS_BL, allLogA0, galaxyData, N);
const gcWork = looGapClosed(WORK_BL, allLogA0, galaxyData, N);

log("  ┌─────────────────────────────────────────────────────────────────┐");
log("  │  BASELINE A — CONSERVATIVE (4 variables)                      │");
log("  │  Variables: MHI + rcWiggliness + envCode + Sigma0_bar         │");
log("  │  Circularity: 2 ZERO + 1 MED + 1 ZERO = mostly clean        │");
log("  │  LOO gap-closed: " + gcCons.gc.toFixed(1) + "%                                       │");
log("  │  Use for: cross-sample claims, publication-grade results      │");
log("  ├─────────────────────────────────────────────────────────────────┤");
log("  │  BASELINE B — EXTENDED (5 variables)                          │");
log("  │  Variables: conservative + meanRun                            │");
log("  │  Circularity: adds MED (RAR residual structure)              │");
log("  │  LOO gap-closed: " + gcWork.gc.toFixed(1) + "%                                       │");
log("  │  Use for: internal analysis, hypothesis generation            │");
log("  └─────────────────────────────────────────────────────────────────┘");
log("");

log("  Confirmed variables registry:");
log("  ┌────────────────┬──────────┬───────┬────────────┬───────────────────┐");
log("  │ Variable       │ Status   │ Circ  │ In BL      │ Door              │");
log("  ├────────────────┼──────────┼───────┼────────────┼───────────────────┤");
log("  │ log(MHI)       │ CONFIRM  │ ZERO  │ A + B      │ Galaxy History    │");
log("  │ rcWiggliness   │ CONFIRM  │ MED   │ A + B      │ Gas/Kinematics   │");
log("  │ envCode        │ CONFIRM  │ ZERO  │ A + B      │ Environment      │");
log("  │ Sigma0_bar     │ CONFIRM  │ ZERO  │ A + B      │ DM Halo          │");
log("  │ meanRun        │ CONFIRM  │ MED   │ B only     │ Internal Struct  │");
log("  │ logVflat       │ CONFIRM* │ ZERO  │ neither    │ History (collin)  │");
log("  │ tBst/Rd        │ CONFIRM* │ MED   │ neither    │ Struct (no LOO)   │");
log("  │ wigRatio       │ CONFIRM* │ MED   │ neither    │ Kin (=wiggliness) │");
log("  └────────────────┴──────────┴───────┴────────────┴───────────────────┘");
log("  * = confirmed signal, not in baseline (collinear or no LOO improvement)");
log("");

log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  PART B: FILTER PARTIAL VARIABLES");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

const partials = [
  { name: 'outerOsc', arr: galaxyData.map(g => g.outerOsc), circ: 'MED', note: 'outer RAR zero-crossings' },
  { name: 'tBst/Rd', arr: galaxyData.map(g => g.logTBstRd), circ: 'MED', note: 'boost transition radius' },
  { name: 'RsVrise', arr: galaxyData.map(g => g.logRsVrise), circ: 'LOW', note: 'peak-V radius / Rdisk' },
  { name: 'innerRise', arr: galaxyData.map(g => g.innerRise), circ: 'LOW', note: 'inner RC rise rate' },
  { name: 'nPts', arr: galaxyData.map(g => g.nPts), circ: 'LOW', note: 'number of data points' },
  { name: 'outerFrac', arr: galaxyData.map(g => g.outerFrac), circ: 'LOW', note: 'outer disk fraction' },
  { name: 'maxJump', arr: galaxyData.map(g => g.maxJump), circ: 'MED', note: 'max RC gradient' },
];

function deepTest(candArr, baseArrs, label) {
  const rx = residualize(candArr, baseArrs);
  const ry = residualize(deltaA0, baseArrs);
  const rObs = corrWith(rx, ry);
  const df = N - 2 - baseArrs.length;
  const tObs = Math.abs(rObs) * Math.sqrt(df) / Math.sqrt(1 - rObs**2 + 1e-10);

  const gc = looGapClosed([...baseArrs, candArr], allLogA0, galaxyData, N);

  const NPERMS = 10000;
  const sh = [...ry];
  let cnt = 0;
  for (let p = 0; p < NPERMS; p++) {
    for (let i = sh.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [sh[i], sh[j]] = [sh[j], sh[i]];
    }
    if (Math.abs(corrWith(rx, sh)) >= Math.abs(rObs)) cnt++;
  }
  const permP = cnt / NPERMS;

  const NBOOT = 2000;
  const bootR = [];
  for (let b = 0; b < NBOOT; b++) {
    const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
    const bX = idx.map(i => candArr[i]);
    const bY = idx.map(i => deltaA0[i]);
    const bBL = baseArrs.map(arr => idx.map(i => arr[i]));
    const brx = residualize(bX, bBL);
    const bry = residualize(bY, bBL);
    bootR.push(corrWith(brx, bry));
  }
  bootR.sort((a,b) => a - b);
  const lo = bootR[Math.floor(NBOOT * 0.025)];
  const hi = bootR[Math.floor(NBOOT * 0.975)];
  const excludesZero = lo > 0 || hi < 0;

  return { rObs, tObs, gc: gc.gc, permP, bootLo: lo, bootHi: hi, excludesZero };
}

log("  Testing each partial against BOTH baselines:");
log("");

const allResults = [];

for (const cand of partials) {
  log("  ── " + cand.name + " (circ=" + cand.circ + ") ──");

  const consTest = deepTest(cand.arr, CONS_BL, 'cons');
  const workTest = deepTest(cand.arr, WORK_BL, 'work');

  const colMR = corrWith(cand.arr, meanRunArr);
  const colMHI = corrWith(cand.arr, mhiArr);
  const colWig = corrWith(cand.arr, wigArr);

  log("    vs Conservative BL (" + gcCons.gc.toFixed(1) + "%):");
  log("      |BL r=" + consTest.rObs.toFixed(3) + " t=" + consTest.tObs.toFixed(1) +
    " LOO=" + consTest.gc.toFixed(1) + "% perm=" + consTest.permP.toFixed(4) +
    " boot=[" + consTest.bootLo.toFixed(3) + "," + consTest.bootHi.toFixed(3) + "]");
  log("    vs Extended BL (" + gcWork.gc.toFixed(1) + "%):");
  log("      |BL r=" + workTest.rObs.toFixed(3) + " t=" + workTest.tObs.toFixed(1) +
    " LOO=" + workTest.gc.toFixed(1) + "% perm=" + workTest.permP.toFixed(4) +
    " boot=[" + workTest.bootLo.toFixed(3) + "," + workTest.bootHi.toFixed(3) + "]");
  log("    Collin: mRun=" + colMR.toFixed(2) + " MHI=" + colMHI.toFixed(2) + " wig=" + colWig.toFixed(2));

  const passesConsPerm = consTest.permP < 0.05;
  const passesConsBoot = consTest.excludesZero;
  const addsOverCons = consTest.gc > gcCons.gc + 1;

  const passesWorkPerm = workTest.permP < 0.05;
  const passesWorkBoot = workTest.excludesZero;
  const addsOverWork = workTest.gc > gcWork.gc + 1;

  let verdict;
  if (passesWorkPerm && passesWorkBoot && addsOverWork) {
    verdict = 'SECONDARY SURVIVOR (adds over extended)';
  } else if (passesConsPerm && passesConsBoot && addsOverCons) {
    verdict = 'SECONDARY SURVIVOR (adds over conservative only)';
  } else if (passesConsPerm && passesConsBoot) {
    verdict = 'MARGINAL (passes tests but no LOO gain)';
  } else if (passesConsPerm || passesConsBoot) {
    verdict = 'WEAK (partial pass only)';
  } else {
    verdict = 'DISCARD';
  }

  log("    VERDICT: " + verdict);
  log("");

  allResults.push({
    name: cand.name, circ: cand.circ, note: cand.note,
    cons: {
      r: +consTest.rObs.toFixed(3), t: +consTest.tObs.toFixed(1),
      loo: +consTest.gc.toFixed(1), permP: +consTest.permP.toFixed(4),
      boot: [+consTest.bootLo.toFixed(3), +consTest.bootHi.toFixed(3)],
      excludesZero: consTest.excludesZero,
    },
    work: {
      r: +workTest.rObs.toFixed(3), t: +workTest.tObs.toFixed(1),
      loo: +workTest.gc.toFixed(1), permP: +workTest.permP.toFixed(4),
      boot: [+workTest.bootLo.toFixed(3), +workTest.bootHi.toFixed(3)],
      excludesZero: workTest.excludesZero,
    },
    collinearity: { meanRun: +colMR.toFixed(2), MHI: +colMHI.toFixed(2), wig: +colWig.toFixed(2) },
    verdict,
  });
}

log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  PART C: SUMMARY TABLE");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");
log("  ┌─────────────┬───────┬────────┬────────┬────────┬────────┬───────────────────────────────┐");
log("  │ Variable    │ Circ  │ Cons t │ Cons p │ Work t │ Work p │ Verdict                       │");
log("  ├─────────────┼───────┼────────┼────────┼────────┼────────┼───────────────────────────────┤");
for (const r of allResults) {
  log("  │ " + r.name.padEnd(11) + " │ " + r.circ.padEnd(5) + " │ " +
    r.cons.t.toFixed(1).padStart(6) + " │ " + r.cons.permP.toFixed(3).padStart(6) + " │ " +
    r.work.t.toFixed(1).padStart(6) + " │ " + r.work.permP.toFixed(3).padStart(6) + " │ " +
    r.verdict.substring(0, 29).padEnd(29) + " │");
}
log("  └─────────────┴───────┴────────┴────────┴────────┴────────┴───────────────────────────────┘");
log("");

const survivors = allResults.filter(r => r.verdict.includes('SURVIVOR'));
const marginals = allResults.filter(r => r.verdict.includes('MARGINAL'));
const discards = allResults.filter(r => r.verdict === 'DISCARD' || r.verdict.includes('WEAK'));

log("  Secondary Survivors: " + (survivors.length > 0 ? survivors.map(s => s.name).join(', ') : 'NONE'));
log("  Marginals:           " + (marginals.length > 0 ? marginals.map(s => s.name).join(', ') : 'NONE'));
log("  Discarded:           " + discards.map(s => s.name).join(', '));
log("");
log("  BASELINES FROZEN. No further baseline changes permitted.");
log("");

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  baselines: {
    conservative: { vars: ['MHI', 'rcWiggliness', 'envCode', 'Sigma0_bar'], gc: +gcCons.gc.toFixed(1) },
    extended: { vars: ['MHI', 'rcWiggliness', 'envCode', 'Sigma0_bar', 'meanRun'], gc: +gcWork.gc.toFixed(1) },
  },
  confirmedRegistry: [
    { name: 'log(MHI)', status: 'CONFIRMED', circ: 'ZERO', baseline: 'A+B', door: 'Galaxy History' },
    { name: 'rcWiggliness', status: 'CONFIRMED', circ: 'MED', baseline: 'A+B', door: 'Gas/Kinematics' },
    { name: 'envCode', status: 'CONFIRMED', circ: 'ZERO', baseline: 'A+B', door: 'Environment' },
    { name: 'Sigma0_bar', status: 'CONFIRMED', circ: 'ZERO', baseline: 'A+B', door: 'DM Halo' },
    { name: 'meanRun', status: 'CONFIRMED', circ: 'MED', baseline: 'B only', door: 'Internal Structure' },
    { name: 'logVflat', status: 'CONFIRMED*', circ: 'ZERO', baseline: 'neither', door: 'History (collinear)' },
    { name: 'tBst/Rd', status: 'CONFIRMED*', circ: 'MED', baseline: 'neither', door: 'Structure (no LOO)' },
    { name: 'wigRatio', status: 'CONFIRMED*', circ: 'MED', baseline: 'neither', door: 'Kin (=wiggliness)' },
  ],
  partialResults: allResults,
  survivors: survivors.map(s => s.name),
  marginals: marginals.map(s => s.name),
  discarded: discards.map(s => s.name),
};

fs.writeFileSync(path.join(__dirname, '../public/phase54-freeze-filter.json'), JSON.stringify(output, null, 2));
log("  Results saved to public/phase54-freeze-filter.json");
