#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "57.0.0";
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
  const rNorm = sorted.map(p => p.r / rMax);

  const residMean = rarResid.reduce((s,v)=>s+v,0)/nPts;
  const residVar = rarResid.map(v=>(v-residMean)**2).reduce((s,v)=>s+v,0)/nPts;
  const residSkew = rarResid.map(v=>((v-residMean)/Math.sqrt(residVar+1e-20))**3).reduce((s,v)=>s+v,0)/nPts;
  const residKurt = rarResid.map(v=>((v-residMean)/Math.sqrt(residVar+1e-20))**4).reduce((s,v)=>s+v,0)/nPts - 3;

  const innerResid = sorted.filter((_,i) => rNorm[i] <= 0.5).map((_,i) => rarResid[i]);
  const outerResid = sorted.filter((_,i) => rNorm[i] > 0.5).map((_,idx) => {
    const actualIdx = sorted.findIndex((_, i) => rNorm[i] > 0.5) + idx;
    return rarResid[actualIdx] || 0;
  });

  const innerRMS = innerResid.length > 1 ? Math.sqrt(innerResid.map(v=>v**2).reduce((s,v)=>s+v,0)/innerResid.length) : 0;
  const outerRMS2 = [];
  for (let i = 0; i < nPts; i++) {
    if (rNorm[i] > 0.5) outerRMS2.push(rarResid[i]);
  }
  const outerRMS = outerRMS2.length > 1 ? Math.sqrt(outerRMS2.map(v=>v**2).reduce((s,v)=>s+v,0)/outerRMS2.length) : 0;
  const innerOuterRatio = outerRMS > 0.001 ? innerRMS / outerRMS : 1;

  let autocorr1 = 0, autocorr2 = 0;
  if (nPts > 3) {
    let s1 = 0, s2 = 0, ss0 = 0;
    for (let i = 0; i < nPts; i++) ss0 += (rarResid[i] - residMean) ** 2;
    for (let i = 0; i < nPts - 1; i++) s1 += (rarResid[i] - residMean) * (rarResid[i+1] - residMean);
    for (let i = 0; i < nPts - 2; i++) s2 += (rarResid[i] - residMean) * (rarResid[i+2] - residMean);
    autocorr1 = ss0 > 0 ? s1 / ss0 : 0;
    autocorr2 = ss0 > 0 ? s2 / ss0 : 0;
  }

  const jumps = [];
  for (let i = 1; i < nPts; i++) jumps.push(Math.abs(rarResid[i] - rarResid[i-1]));
  const rmsJump = jumps.length > 0 ? Math.sqrt(jumps.map(v=>v**2).reduce((s,v)=>s+v,0)/jumps.length) : 0;
  const maxJump = jumps.length > 0 ? Math.max(...jumps) : 0;
  const smoothness = residVar > 0 ? rmsJump / Math.sqrt(residVar) : 1;

  let harmonicPower2 = 0, harmonicPower3 = 0;
  if (nPts >= 6) {
    let cos2 = 0, sin2 = 0, cos3 = 0, sin3 = 0;
    for (let i = 0; i < nPts; i++) {
      const phi = 2 * Math.PI * rNorm[i];
      cos2 += rarResid[i] * Math.cos(2 * phi);
      sin2 += rarResid[i] * Math.sin(2 * phi);
      cos3 += rarResid[i] * Math.cos(3 * phi);
      sin3 += rarResid[i] * Math.sin(3 * phi);
    }
    harmonicPower2 = Math.sqrt(cos2**2 + sin2**2) / nPts;
    harmonicPower3 = Math.sqrt(cos3**2 + sin3**2) / nPts;
  }

  let radialGradient = 0;
  if (nPts > 3) {
    let sxy2 = 0, sxx2 = 0;
    const mr = rNorm.reduce((s,v)=>s+v,0)/nPts;
    for (let i = 0; i < nPts; i++) {
      sxy2 += (rNorm[i] - mr) * rarResid[i];
      sxx2 += (rNorm[i] - mr) ** 2;
    }
    radialGradient = sxx2 > 0 ? sxy2 / sxx2 : 0;
  }

  const absResid = rarResid.map(Math.abs);
  const medianAbsResid = [...absResid].sort((a,b)=>a-b)[Math.floor(nPts/2)];
  const robustRMS = medianAbsResid * 1.4826;

  let runningVar = 0;
  if (nPts >= 5) {
    const winSize = Math.max(3, Math.floor(nPts / 4));
    const localVars = [];
    for (let i = 0; i <= nPts - winSize; i++) {
      const win = rarResid.slice(i, i + winSize);
      const wm = win.reduce((s,v)=>s+v,0)/win.length;
      localVars.push(win.map(v=>(v-wm)**2).reduce((s,v)=>s+v,0)/win.length);
    }
    const meanLV = localVars.reduce((s,v)=>s+v,0)/localVars.length;
    runningVar = Math.sqrt(localVars.map(v=>(v-meanLV)**2).reduce((s,v)=>s+v,0)/localVars.length);
  }

  let curvature = 0;
  if (nPts >= 3) {
    let totalCurv = 0;
    for (let i = 1; i < nPts - 1; i++) {
      totalCurv += Math.abs(rarResid[i+1] - 2*rarResid[i] + rarResid[i-1]);
    }
    curvature = totalCurv / (nPts - 2);
  }

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, rcWiggliness, envCode, logSigma0, logMeanRun,
    residSkew, residKurt,
    innerOuterRatio: Math.log10(Math.max(0.01, innerOuterRatio)),
    autocorr1, autocorr2,
    rmsJump, maxJump, smoothness,
    harmonicPower2: Math.log10(Math.max(1e-4, harmonicPower2)),
    harmonicPower3: Math.log10(Math.max(1e-4, harmonicPower3)),
    radialGradient,
    robustRMS,
    runningVar: Math.log10(Math.max(1e-6, runningVar)),
    curvature,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);

const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const sigma0Arr = galaxyData.map(g => g.logSigma0);
const meanRunArr = galaxyData.map(g => g.logMeanRun);
const CONS_BL = [mhiArr, wigArr, envArr, sigma0Arr];
const WORK_BL = [mhiArr, wigArr, envArr, sigma0Arr, meanRunArr];

const gcCons = looGapClosed(CONS_BL, allLogA0, galaxyData, N);
const gcWork = looGapClosed(WORK_BL, allLogA0, galaxyData, N);

log("");
log("╔══════════════════════════════════════════════════════════════════════════════════╗");
log("║  PHASE 57: 2D KINEMATICS DECISIVE TEST                                        ║");
log("║  Version " + VERSION + ", " + N + " galaxies                                           ║");
log("╚══════════════════════════════════════════════════════════════════════════════════╝");
log("");
log("  NOTE: True 2D velocity fields not available in SPARC.");
log("  This phase extracts ALL possible kinematic proxies from the 1D rotation curves");
log("  that correlate with 2D kinematic properties (asymmetry, harmonic structure,");
log("  non-circular motions, coherence patterns).");
log("");
log("  Conservative BL: " + gcCons.gc.toFixed(1) + "% | Extended BL: " + gcWork.gc.toFixed(1) + "%");
log("");

const kinVars = [
  { name: 'residSkew', arr: galaxyData.map(g => g.residSkew), note: '1D proxy for kinematic lopsidedness', circ: 'MED' },
  { name: 'residKurt', arr: galaxyData.map(g => g.residKurt), note: '1D proxy for non-Gaussian velocity structure', circ: 'MED' },
  { name: 'inOutRatio', arr: galaxyData.map(g => g.innerOuterRatio), note: '1D proxy for radial non-circular power ratio', circ: 'MED' },
  { name: 'autocorr1', arr: galaxyData.map(g => g.autocorr1), note: '1D proxy for coherent kinematic structures', circ: 'MED' },
  { name: 'autocorr2', arr: galaxyData.map(g => g.autocorr2), note: '1D proxy for bistable kinematic patterns', circ: 'MED' },
  { name: 'rmsJump', arr: galaxyData.map(g => g.rmsJump), note: '1D proxy for local velocity irregularity', circ: 'MED' },
  { name: 'maxJump', arr: galaxyData.map(g => g.maxJump), note: '1D proxy for strongest kinematic discontinuity', circ: 'MED' },
  { name: 'smoothness', arr: galaxyData.map(g => g.smoothness), note: 'jump-to-scatter ratio (non-circular proxy)', circ: 'MED' },
  { name: 'harm2', arr: galaxyData.map(g => g.harmonicPower2), note: '1D proxy for m=2 (bar/bisymmetric flow)', circ: 'MED' },
  { name: 'harm3', arr: galaxyData.map(g => g.harmonicPower3), note: '1D proxy for m=3 harmonic distortion', circ: 'MED' },
  { name: 'radGrad', arr: galaxyData.map(g => g.radialGradient), note: '1D proxy for radial flow / systematic tilt', circ: 'MED' },
  { name: 'robustRMS', arr: galaxyData.map(g => g.robustRMS), note: 'MAD-based scatter (outlier-resistant)', circ: 'MED' },
  { name: 'runningVar', arr: galaxyData.map(g => g.runningVar), note: '1D proxy for spatially varying non-circular motion', circ: 'HIGH' },
  { name: 'curvature', arr: galaxyData.map(g => g.curvature), note: '1D proxy for second-derivative structure', circ: 'MED' },
];

log("  ┌──────────────┬──────┬────────┬────────┬────────┬────────┬────────┬────────┬───────────┐");
log("  │ Variable     │ Circ │ |Cs r  │ Cs LOO │ Cs p   │ |Wk r  │ Wk LOO │ Wk p   │ Verdict   │");
log("  ├──────────────┼──────┼────────┼────────┼────────┼────────┼────────┼────────┼───────────┤");

const deltaA0 = allLogA0.map(v => v - allLogA0.reduce((s,v)=>s+v,0)/N);
const NPERMS = 10000;

const kinResults = [];

for (const kv of kinVars) {
  const rxC = residualize(kv.arr, CONS_BL);
  const ryC = residualize(deltaA0, CONS_BL);
  const rC = corrWith(rxC, ryC);
  const gcC = looGapClosed([...CONS_BL, kv.arr], allLogA0, galaxyData, N);

  const rxW = residualize(kv.arr, WORK_BL);
  const ryW = residualize(deltaA0, WORK_BL);
  const rW = corrWith(rxW, ryW);
  const gcW = looGapClosed([...WORK_BL, kv.arr], allLogA0, galaxyData, N);

  let cntC = 0, cntW = 0;
  const shC = [...ryC], shW = [...ryW];
  for (let p = 0; p < NPERMS; p++) {
    for (let i = shC.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [shC[i], shC[j]] = [shC[j], shC[i]];
      [shW[i], shW[j]] = [shW[j], shW[i]];
    }
    if (Math.abs(corrWith(rxC, shC)) >= Math.abs(rC)) cntC++;
    if (Math.abs(corrWith(rxW, shW)) >= Math.abs(rW)) cntW++;
  }
  const permC = cntC / NPERMS;
  const permW = cntW / NPERMS;

  const collinWig = corrWith(kv.arr, wigArr);
  const collinMR = corrWith(kv.arr, meanRunArr);

  let verdict;
  if (permW < 0.05 && gcW.gc > gcWork.gc + 1) verdict = 'SURVIVOR';
  else if (permC < 0.05 && gcC.gc > gcCons.gc + 1) verdict = 'PARTIAL';
  else verdict = 'FAIL';

  kinResults.push({
    name: kv.name, circ: kv.circ, note: kv.note,
    cons: { r: +rC.toFixed(3), loo: +gcC.gc.toFixed(1), permP: +permC.toFixed(4) },
    work: { r: +rW.toFixed(3), loo: +gcW.gc.toFixed(1), permP: +permW.toFixed(4) },
    collinWig: +collinWig.toFixed(2), collinMR: +collinMR.toFixed(2),
    verdict,
  });

  log("  │ " + kv.name.padEnd(12) + " │ " + kv.circ.padEnd(4) + " │ " +
    rC.toFixed(3).padStart(6) + " │ " + (gcC.gc.toFixed(1)+"%").padStart(6) + " │ " +
    permC.toFixed(3).padStart(6) + " │ " +
    rW.toFixed(3).padStart(6) + " │ " + (gcW.gc.toFixed(1)+"%").padStart(6) + " │ " +
    permW.toFixed(3).padStart(6) + " │ " + verdict.padEnd(9) + " │");
}
log("  └──────────────┴──────┴────────┴────────┴────────┴────────┴────────┴────────┴───────────┘");

log("");
log("  Collinearity with existing kinematic variables:");
log("  ┌──────────────┬────────┬────────┐");
log("  │ Variable     │ r(wig) │ r(mRun)│");
log("  ├──────────────┼────────┼────────┤");
for (const kr of kinResults) {
  log("  │ " + kr.name.padEnd(12) + " │ " + kr.collinWig.toFixed(2).padStart(6) + " │ " + kr.collinMR.toFixed(2).padStart(6) + " │");
}
log("  └──────────────┴────────┴────────┘");

const survivors = kinResults.filter(r => r.verdict === 'SURVIVOR');
const partials = kinResults.filter(r => r.verdict === 'PARTIAL');

log("");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  MODEL COMPARISON: M57-0, M57-1, M57-2");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

const m570 = looGapClosed(CONS_BL, allLogA0, galaxyData, N);
const m571 = looGapClosed([...CONS_BL, wigArr, meanRunArr], allLogA0, galaxyData, N);

const kinSurvivorArrs = survivors.map(s => kinVars.find(kv => kv.name === s.name).arr);
const kinPartialArrs = partials.map(s => kinVars.find(kv => kv.name === s.name).arr);
const allKinArrs = [...kinSurvivorArrs, ...kinPartialArrs];

const m572 = allKinArrs.length > 0
  ? looGapClosed([...CONS_BL, ...allKinArrs], allLogA0, galaxyData, N)
  : { gc: m570.gc };

const m57combo = allKinArrs.length > 0
  ? looGapClosed([...CONS_BL, wigArr, meanRunArr, ...allKinArrs], allLogA0, galaxyData, N)
  : m571;

log("  M57-0 (conservative baseline only):                LOO = " + m570.gc.toFixed(1) + "%");
log("  M57-1 (conservative + rcWig + meanRun):            LOO = " + m571.gc.toFixed(1) + "%");
log("  M57-2 (conservative + 2D kinematic proxies):       LOO = " + m572.gc.toFixed(1) + "%");
if (allKinArrs.length > 0) {
  log("  M57-combo (conservative + wig + mRun + 2D kin):    LOO = " + m57combo.gc.toFixed(1) + "%");
}
log("");

log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  ABSORPTION TEST: Do 2D proxies absorb rcWiggliness / meanRun?");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

if (partials.length > 0 || survivors.length > 0) {
  const sigKin = [...survivors, ...partials];
  for (const sk of sigKin) {
    const skArr = kinVars.find(kv => kv.name === sk.name).arr;
    const withKin = [...CONS_BL, skArr];
    const rxWig = residualize(wigArr, withKin);
    const ryWig = residualize(deltaA0, withKin);
    const wigAfterKin = corrWith(rxWig, ryWig);
    const rxMR = residualize(meanRunArr, withKin);
    const ryMR = residualize(deltaA0, withKin);
    const mrAfterKin = corrWith(rxMR, ryMR);

    log("  After controlling for " + sk.name + ":");
    log("    rcWiggliness partial r: " + wigAfterKin.toFixed(3) + " (was 0.345 in baseline A)");
    log("    meanRun partial r:      " + mrAfterKin.toFixed(3) + " (was 0.368 in baseline B)");
    log("    Absorption: wig=" + ((1 - Math.abs(wigAfterKin)/0.345)*100).toFixed(0) + "%, mRun=" + ((1 - Math.abs(mrAfterKin)/0.368)*100).toFixed(0) + "%");
    log("");
  }
} else {
  log("  No significant kinematic proxies survived. Cannot test absorption.");
  log("  rcWiggliness and meanRun remain the best available kinematic metrics.");
  log("");
}

log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  PHASE 57 — VERDICT");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

const hasSurvivor = survivors.length > 0;
const hasPartial = partials.length > 0;

if (hasSurvivor) {
  log("  RESULT: 2D kinematic proxies ADD beyond baseline.");
  log("  Survivors: " + survivors.map(s => s.name).join(", "));
  log("  This suggests the signal may partly be non-circular motions.");
} else if (hasPartial) {
  log("  RESULT: Some 2D kinematic proxies show partial signal.");
  log("  Partials: " + partials.map(s => s.name).join(", "));
  log("  But they do not improve beyond the extended baseline.");
  log("  rcWiggliness and meanRun already capture most 1D-accessible kinematic info.");
} else {
  log("  RESULT: NO 2D kinematic proxy adds independent information.");
  log("  All 14 metrics either fail outright or are absorbed by rcWiggliness/meanRun.");
  log("  INTERPRETATION:");
  log("    rcWiggliness + meanRun already capture all kinematic information");
  log("    extractable from 1D rotation curves.");
  log("    Whether TRUE 2D velocity fields add more requires actual IFU/HI cube data.");
  log("    The kinematic signal is DEEPER than non-circular motions visible in 1D RCs.");
}
log("");
log("  DATA LIMITATION:");
log("  Without actual 2D velocity fields (THINGS/LITTLE THINGS HI cubes or IFU data),");
log("  we cannot test true kinematic asymmetry, harmonic A amplitudes, or");
log("  radial/bisymmetric flow components directly.");
log("  This remains an open question for external data.");
log("");

const output = {
  version: VERSION, timestamp: new Date().toISOString(), nGalaxies: N,
  conservativeBL: +gcCons.gc.toFixed(1),
  extendedBL: +gcWork.gc.toFixed(1),
  dataLimitation: "No true 2D velocity fields available; all metrics are 1D RC proxies",
  kinematicVars: kinResults,
  modelComparison: {
    M570_consBaseline: +m570.gc.toFixed(1),
    M571_consWithKin: +m571.gc.toFixed(1),
    M572_consWith2D: +m572.gc.toFixed(1),
  },
  survivors: survivors.map(s => s.name),
  partials: partials.map(s => s.name),
  verdict: hasSurvivor ? '2D_ADDS' : hasPartial ? 'PARTIAL_SIGNAL' : 'NO_INDEPENDENT_INFO',
};

fs.writeFileSync(path.join(__dirname, '../public/phase57-2d-kinematics.json'), JSON.stringify(output, null, 2));
log("  Results saved to public/phase57-2d-kinematics.json");
