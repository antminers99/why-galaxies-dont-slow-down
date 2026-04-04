#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "48.0.0";
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

  const halfLen = Math.floor(sorted.length / 2);
  const innerResid = sorted.slice(0, halfLen).map(p => {
    const gb = Math.pow(10, p.log_g_bar);
    return p.log_g_obs - Math.log10(mcgaughRAR(gb, fit.a0));
  });
  const outerResid = sorted.slice(halfLen).map(p => {
    const gb = Math.pow(10, p.log_g_bar);
    return p.log_g_obs - Math.log10(mcgaughRAR(gb, fit.a0));
  });

  const innerMean = innerResid.length > 0 ? innerResid.reduce((s,v)=>s+v,0)/innerResid.length : 0;
  const outerMean = outerResid.length > 0 ? outerResid.reduce((s,v)=>s+v,0)/outerResid.length : 0;
  const inOutDiff = innerMean - outerMean;

  const innerRMS = innerResid.length > 0 ? Math.sqrt(innerResid.map(v=>v**2).reduce((s,v)=>s+v,0)/innerResid.length) : 0;
  const outerRMS = outerResid.length > 0 ? Math.sqrt(outerResid.map(v=>v**2).reduce((s,v)=>s+v,0)/outerResid.length) : 0;
  const logInnerRMS = Math.log10(innerRMS > 0.001 ? innerRMS : 0.001);
  const logOuterRMS = Math.log10(outerRMS > 0.001 ? outerRMS : 0.001);
  const rmsRatio = outerRMS > 0.001 ? innerRMS / outerRMS : 1;
  const logRmsRatio = Math.log10(rmsRatio > 0.01 ? rmsRatio : 0.01);

  const innerSlope = (() => {
    if (innerPts.length < 3) return 0;
    const x = innerPts.map(p => p.log_g_bar);
    const y = innerPts.map(p => p.log_g_obs);
    const mx = x.reduce((s,v)=>s+v,0)/x.length;
    const my = y.reduce((s,v)=>s+v,0)/y.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < x.length; i++) { sxy += (x[i]-mx)*(y[i]-my); sxx += (x[i]-mx)**2; }
    return sxx > 0 ? sxy / sxx : 0;
  })();

  const outerSlope = (() => {
    if (outerPts.length < 3) return 0;
    const x = outerPts.map(p => p.log_g_bar);
    const y = outerPts.map(p => p.log_g_obs);
    const mx = x.reduce((s,v)=>s+v,0)/x.length;
    const my = y.reduce((s,v)=>s+v,0)/y.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < x.length; i++) { sxy += (x[i]-mx)*(y[i]-my); sxx += (x[i]-mx)**2; }
    return sxx > 0 ? sxy / sxx : 0;
  })();

  const slopeDiff = innerSlope - outerSlope;

  const SBdisk = sp.SBdisk || 21;
  const inc = sp.inc || 60;

  const innerGbar = innerPts.length > 0 ? innerPts.map(p=>p.log_g_bar).reduce((s,v)=>s+v,0)/innerPts.length : -10;
  const outerGbar = outerPts.length > 0 ? outerPts.map(p=>p.log_g_bar).reduce((s,v)=>s+v,0)/outerPts.length : -10;
  const gbarContrast = innerGbar - outerGbar;

  const innerGobs = innerPts.length > 0 ? innerPts.map(p=>p.log_g_obs).reduce((s,v)=>s+v,0)/innerPts.length : -10;
  const outerGobs = outerPts.length > 0 ? outerPts.map(p=>p.log_g_obs).reduce((s,v)=>s+v,0)/outerPts.length : -10;
  const gobsContrast = innerGobs - outerGobs;

  const boostInner = innerPts.length > 0
    ? innerPts.map(p => p.log_g_obs - p.log_g_bar).reduce((s,v)=>s+v,0)/innerPts.length : 0;
  const boostOuter = outerPts.length > 0
    ? outerPts.map(p => p.log_g_obs - p.log_g_bar).reduce((s,v)=>s+v,0)/outerPts.length : 0;
  const boostGradient = boostOuter - boostInner;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, logVflat, rcWiggliness, envCode, logSigma0,
    logMeanRun,
    inOutDiff, logInnerRMS, logOuterRMS, logRmsRatio,
    slopeDiff, gbarContrast, gobsContrast, boostGradient,
    innerSlope, outerSlope,
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
const vflatArr = galaxyData.map(g => g.logVflat);
const meanRunArr = galaxyData.map(g => g.logMeanRun);

const CONS_BASELINE = [mhiArr, wigArr, envArr, sigma0Arr];
const WORK_BASELINE = [mhiArr, wigArr, envArr, sigma0Arr, meanRunArr];

log("");
log("=".repeat(80));
log("  PHASE 48: Internal Structure — Inner-Outer Difference");
log("  FAST PROTOCOL");
log("  Conservative BL: MHI + wig + env + Sigma0_bar");
log("  Working BL:      + meanRun (confirmed Phase 47)");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  CONTEXT:");
log("  Do inner vs outer halves of the RC behave differently");
log("  in ways that predict a0 variation, beyond what meanRun");
log("  already captures? If yes → independent structural signal.");
log("  If no → redundant with meanRun.");
log("");

log("  PROXIES:");
log("  1) inOutDiff: inner - outer mean RAR residual");
log("  2) innerRMS: RMS of inner-half residuals");
log("  3) outerRMS: RMS of outer-half residuals");
log("  4) rmsRatio: inner/outer RMS ratio");
log("  5) slopeDiff: inner RAR slope - outer RAR slope");
log("  6) gbarContr: inner - outer mean log(g_bar)");
log("  7) gobsContr: inner - outer mean log(g_obs)");
log("  8) boostGrad: outer boost - inner boost (gobs-gbar gradient)");
log("");

const candidates = [
  { name: 'inOutDiff', arr: galaxyData.map(g => g.inOutDiff), circ: 'HIGH', note: 'inner-outer mean resid' },
  { name: 'innerRMS', arr: galaxyData.map(g => g.logInnerRMS), circ: 'HIGH', note: 'inner half residual RMS' },
  { name: 'outerRMS', arr: galaxyData.map(g => g.logOuterRMS), circ: 'HIGH', note: 'outer half residual RMS' },
  { name: 'rmsRatio', arr: galaxyData.map(g => g.logRmsRatio), circ: 'HIGH', note: 'inner/outer RMS ratio' },
  { name: 'slopeDif', arr: galaxyData.map(g => g.slopeDiff), circ: 'MED', note: 'inner-outer RAR slope diff' },
  { name: 'gbarCntr', arr: galaxyData.map(g => g.gbarContrast), circ: 'LOW', note: 'inner-outer g_bar contrast' },
  { name: 'gobsCntr', arr: galaxyData.map(g => g.gobsContrast), circ: 'MED', note: 'inner-outer g_obs contrast' },
  { name: 'boostGrd', arr: galaxyData.map(g => g.boostGradient), circ: 'HIGH', note: 'outer-inner boost gradient' },
];

function looGapClosed(baseArrays, extraArr) {
  const featArrays = [...baseArrays];
  if (extraArr) featArrays.push(extraArr);

  let totalSS_model = 0, totalSS_m0 = 0, totalSS_free = 0, totalN = 0;
  for (let i = 0; i < N; i++) {
    const trainIdx = [...Array(N).keys()].filter(j => j !== i);
    const trainY = trainIdx.map(j => allLogA0[j]);
    const mu = trainY.reduce((s,v)=>s+v,0)/trainY.length;
    totalSS_m0 += predictSS(Math.pow(10, mu), galaxyData[i].pts);
    totalSS_free += predictSS(galaxyData[i].a0, galaxyData[i].pts);
    totalN += galaxyData[i].pts.length;

    const trainX = trainIdx.map(j => featArrays.map(arr => arr[j]));
    const reg = linReg(trainX, trainY);
    const testX = featArrays.map(arr => arr[i]);
    let pred = reg.intercept + reg.coefs.reduce((s, c, j) => s + c * testX[j], 0);
    pred = Math.max(2.5, Math.min(4.5, pred));
    totalSS_model += predictSS(Math.pow(10, pred), galaxyData[i].pts);
  }
  const m0rms = Math.sqrt(totalSS_m0 / totalN);
  const m6rms = Math.sqrt(totalSS_free / totalN);
  const gap = m0rms - m6rms;
  const modelRms = Math.sqrt(totalSS_model / totalN);
  return { gc: gap > 0 ? (m0rms - modelRms) / gap * 100 : 0 };
}

const consGC = looGapClosed(CONS_BASELINE, null);
const workGC = looGapClosed(WORK_BASELINE, null);
log("  CONSERVATIVE BASELINE: " + consGC.gc.toFixed(1) + "% gap closure");
log("  WORKING BASELINE (+meanRun): " + workGC.gc.toFixed(1) + "% gap closure");
log("");

log("  COLLINEARITY with meanRun:");
for (const c of candidates) {
  const rMR = corrWith(c.arr, meanRunArr);
  const rMhi = corrWith(c.arr, mhiArr);
  const rWig = corrWith(c.arr, wigArr);
  const rVf = corrWith(c.arr, vflatArr);
  log("  " + c.name.padEnd(12) + " vs mRun:" + rMR.toFixed(2).padStart(6) +
    "  MHI:" + rMhi.toFixed(2).padStart(6) + "  wig:" + rWig.toFixed(2).padStart(6) +
    "  Vf:" + rVf.toFixed(2).padStart(6));
}
log("");

log("=".repeat(80));
log("  STEP 1: FAST SCREEN (against BOTH baselines)");
log("=".repeat(80));
log("");

log("  ┌─────────────────────────────────────────────────────────────────────────────────────────────────────────────────┐");
log("  │  Proxy        Circ   |consBL t  consLOO  |workBL t  workLOO  Adds/c  Adds/w  Verdict                        │");
log("  ├─────────────────────────────────────────────────────────────────────────────────────────────────────────────────┤");

const results = [];
for (const cand of candidates) {
  const rxC = residualize(cand.arr, CONS_BASELINE);
  const ryC = residualize(deltaA0, CONS_BASELINE);
  const rCons = corrWith(rxC, ryC);
  const tCons = Math.abs(rCons) * Math.sqrt(N-6) / Math.sqrt(1-rCons**2+1e-10);
  const gcCons = looGapClosed(CONS_BASELINE, cand.arr);

  const rxW = residualize(cand.arr, WORK_BASELINE);
  const ryW = residualize(deltaA0, WORK_BASELINE);
  const rWork = corrWith(rxW, ryW);
  const tWork = Math.abs(rWork) * Math.sqrt(N-7) / Math.sqrt(1-rWork**2+1e-10);
  const gcWork = looGapClosed(WORK_BASELINE, cand.arr);

  const addsCons = gcCons.gc > consGC.gc + 1;
  const addsWork = gcWork.gc > workGC.gc + 1;

  const rRaw = corrWith(cand.arr, deltaA0);
  const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);

  const failFast = tRaw < 1.0 && tCons < 1.0 && tWork < 1.0;
  const promising = tWork >= 1.65 || addsWork;
  let verdict = failFast ? "FAIL-FAST" : promising ? "PROMISING" : "WEAK";

  if (cand.circ === 'HIGH' && promising) verdict = "PROMISING(circ!)";

  results.push({
    ...cand, rRaw, tRaw,
    rCons, tCons, gcCons: gcCons.gc, addsCons,
    rWork, tWork, gcWork: gcWork.gc, addsWork,
    verdict,
  });

  log("  │  " + cand.name.padEnd(12) + cand.circ.padEnd(5) +
    tCons.toFixed(1).padStart(10) + (gcCons.gc.toFixed(1)+"%").padStart(8) +
    tWork.toFixed(1).padStart(10) + (gcWork.gc.toFixed(1)+"%").padStart(8) +
    (addsCons ? "  YES " : "  NO  ") +
    (addsWork ? "  YES " : "  NO  ") +
    verdict.padStart(22) + "   │");
}
log("  └─────────────────────────────────────────────────────────────────────────────────────────────────────────────────┘");
log("");

const promising2 = results.filter(r => r.verdict.startsWith("PROMISING"));

if (promising2.length === 0) {
  log("  ══════════════════════════════════════════════════════════");
  log("  NO inner-outer proxy passes fast screen against working BL.");
  log("  All are redundant with meanRun or have no signal.");
  log("  ══════════════════════════════════════════════════════════");
} else {
  const lowCirc = promising2.filter(r => r.circ !== 'HIGH');
  const deepCandidates = lowCirc.length > 0 ? lowCirc : promising2;

  log("  PROMISING: " + promising2.map(p => p.name + "(circ=" + p.circ + ")").join(", "));
  log("  Deep testing: " + deepCandidates.map(p => p.name).join(", "));
  log("");

  for (const cand of deepCandidates) {
    log("=".repeat(80));
    log("  DEEP TEST: " + cand.name + " (circ=" + cand.circ + ") against WORKING BL");
    log("=".repeat(80));
    log("");

    log("  Confounder stripping:");
    for (const [label, ctrls] of [
      ['raw', []],
      ['|MHI', [mhiArr]],
      ['|consBL', CONS_BASELINE],
      ['|workBL', WORK_BASELINE],
      ['|workBL+Vflat', [...WORK_BASELINE, vflatArr]],
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

    log("  Permutation (10000) against workBL:");
    const NPERMS = 10000;
    const rx2 = residualize(cand.arr, WORK_BASELINE);
    const ry2 = residualize(deltaA0, WORK_BASELINE);
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
    log("    |r|=" + obsR.toFixed(3) + " p=" + (cnt/NPERMS).toFixed(4) +
      (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "**" : ""));
    log("");

    log("  Jackknife:");
    const jackR = [];
    for (let i = 0; i < N; i++) {
      const idx = [...Array(N).keys()].filter(j => j !== i);
      const subX = idx.map(j => cand.arr[j]);
      const subY = idx.map(j => deltaA0[j]);
      const subBL = WORK_BASELINE.map(arr => idx.map(j => arr[j]));
      const rx3 = residualize(subX, subBL);
      const ry3 = residualize(subY, subBL);
      jackR.push(corrWith(rx3, ry3));
    }
    const jackMean = jackR.reduce((s,v)=>s+v,0)/N;
    const jackSD = Math.sqrt(jackR.map(v => (v - jackMean)**2).reduce((s,v)=>s+v,0)/(N-1));
    const signFlips = jackR.filter(r => Math.sign(r) !== Math.sign(cand.rWork)).length;
    log("    Mean=" + jackMean.toFixed(3) + " SD=" + jackSD.toFixed(3) + " flips=" + signFlips + "/" + N);
    log("");

    log("  Bootstrap CI (2000):");
    const NBOOT = 2000;
    const bootR = [];
    for (let b = 0; b < NBOOT; b++) {
      const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
      const bX = idx.map(i => cand.arr[i]);
      const bY = idx.map(i => deltaA0[i]);
      const bBL = WORK_BASELINE.map(arr => idx.map(i => arr[i]));
      const rx3 = residualize(bX, bBL);
      const ry3 = residualize(bY, bBL);
      bootR.push(corrWith(rx3, ry3));
    }
    bootR.sort((a,b) => a - b);
    const lo = bootR[Math.floor(NBOOT * 0.025)];
    const hi = bootR[Math.floor(NBOOT * 0.975)];
    log("    95% CI: [" + lo.toFixed(3) + ", " + hi.toFixed(3) + "]" +
      (lo > 0 || hi < 0 ? " excludes zero" : " crosses zero"));
    log("");
  }
}

log("");
log("=".repeat(80));
log("  PHASE 48 — SUMMARY");
log("=".repeat(80));
log("");

log("  Results summary:");
for (const r of results) {
  log("    " + r.name.padEnd(12) + " circ=" + r.circ.padEnd(5) +
    " |consBL t=" + r.tCons.toFixed(1).padStart(4) +
    " |workBL t=" + r.tWork.toFixed(1).padStart(4) +
    " workLOO=" + r.gcWork.toFixed(1).padStart(5) + "% " + r.verdict);
}
log("");
log("  Conservative BL: " + consGC.gc.toFixed(1) + "%");
log("  Working BL (+meanRun): " + workGC.gc.toFixed(1) + "%");
log("");
log("=".repeat(80));

const output = {
  version: VERSION, timestamp: new Date().toISOString(), nGalaxies: N,
  conservativeBL: { gc: +consGC.gc.toFixed(1) },
  workingBL: { gc: +workGC.gc.toFixed(1), includes: 'meanRun' },
  proxies: results.map(r => ({
    name: r.name, circ: r.circ, note: r.note,
    consBLr: +r.rCons.toFixed(3), consBLt: +r.tCons.toFixed(1), consLOO: +r.gcCons.toFixed(1),
    workBLr: +r.rWork.toFixed(3), workBLt: +r.tWork.toFixed(1), workLOO: +r.gcWork.toFixed(1),
    addsCons: r.addsCons, addsWork: r.addsWork, verdict: r.verdict,
  })),
};

fs.writeFileSync(path.join(__dirname, '../public/phase48-inner-outer.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase48-inner-outer.json");
