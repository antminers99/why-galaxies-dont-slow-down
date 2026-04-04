#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "59.0.0";
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
  const logVflat = sp.Vflat > 0 ? Math.log10(sp.Vflat) : 2.0;

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
  let rl = 1;
  const residRun = [];
  for (let i = 1; i < rarResid.length; i++) {
    const s = Math.sign(rarResid[i]);
    if (s === currentSign) { rl++; }
    else { residRun.push(rl); rl = 1; currentSign = s; }
  }
  residRun.push(rl);
  const meanRunLen = residRun.reduce((s,v)=>s+v,0) / residRun.length;
  const logMeanRun = Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1);

  const Mstar_sun = L * 1e9 * 0.5;
  const logMstar = Math.log10(Mstar_sun > 0 ? Mstar_sun : 1e6);
  const logL36 = Math.log10(L > 0 ? L : 0.01);
  const SBdisk = sp.SBdisk || 21.0;
  const SBeff = sp.SBeff || 22.0;

  const fgas = sp.MHI * 1.33 / Mbar_sun;
  const T = sp.T || gal.T || 5;

  const mlSB = -0.4 * (SBdisk - 21.0);
  const mlMstar = logMstar / 10;
  const mlT = T < 3 ? 0.7 : (T < 6 ? 0.5 : 0.3);
  const mlVflat = logVflat / 2.5;
  const mlColor = T < 2 ? 1.0 : (T < 4 ? 0.7 : (T < 7 ? 0.4 : 0.2));

  const logMLsb = Math.log10(Math.max(0.01, mlSB + 1));
  const logMstar_L = logMstar - logL36;
  const sbGrad = SBeff - SBdisk;
  const logStellarDensity = logMstar - 2 * Math.log10(Rdisk > 0 ? Rdisk : 1);

  const dynamicML = sp.Vflat > 0 ? Math.log10(sp.Vflat**2 * rMax / (L * 1e9 * 0.5 > 0 ? L * 1e9 * 0.5 : 1e6)) : 0;

  const barGravDom = (() => {
    let count = 0;
    for (const p of pts) {
      if (p.log_g_bar > p.log_g_obs - 0.3) count++;
    }
    return count / pts.length;
  })();

  const innerBarFrac = (() => {
    let count = 0;
    for (const p of innerPts) {
      if (p.log_g_bar > p.log_g_obs - 0.3) count++;
    }
    return innerPts.length > 0 ? count / innerPts.length : 0;
  })();

  const gbarSpread = (() => {
    const gbs = pts.map(p => p.log_g_bar);
    return Math.max(...gbs) - Math.min(...gbs);
  })();

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, rcWiggliness, envCode, logSigma0, logMeanRun,
    logMstar, logL36, SBdisk, SBeff, T, logVflat,
    logMLsb, logMstar_L, sbGrad, logStellarDensity,
    dynamicML, barGravDom, innerBarFrac, gbarSpread,
    mlColor: Math.log10(mlColor > 0 ? mlColor : 0.1),
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const deltaA0 = allLogA0.map(v => v - allLogA0.reduce((s,v)=>s+v,0)/N);

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
log("║  PHASE 59: TRUE STELLAR MASS-TO-LIGHT RATIO TEST                              ║");
log("║  Version " + VERSION + ", " + N + " galaxies                                           ║");
log("╚══════════════════════════════════════════════════════════════════════════════════╝");
log("");
log("  HYPOTHESIS: The Sigma0 signal and part of the unexplained gap come from");
log("  stellar population variations (Y*) that a constant M/L = 0.5 misses.");
log("  If true, Y*-sensitive proxies should absorb Sigma0 and close more gap.");
log("");
log("  NOTE: True SPS-derived Y* not available. Using SPARC observables as proxies:");
log("  surface brightness, color (T-type), M*/L ratio, stellar density, etc.");
log("");
log("  Conservative BL: " + gcCons.gc.toFixed(1) + "% | Extended BL: " + gcWork.gc.toFixed(1) + "%");
log("");

const mlVars = [
  { name: 'logMLsb', arr: galaxyData.map(g => g.logMLsb), note: 'SB-based M/L proxy', circ: 'LOW' },
  { name: 'logM*/L', arr: galaxyData.map(g => g.logMstar_L), note: 'stellar mass to 3.6um luminosity ratio', circ: 'LOW' },
  { name: 'sbGrad', arr: galaxyData.map(g => g.sbGrad), note: 'SB gradient (SBeff - SBdisk)', circ: 'ZERO' },
  { name: 'logStelD', arr: galaxyData.map(g => g.logStellarDensity), note: 'stellar surface density M*/Rd^2', circ: 'ZERO' },
  { name: 'dynML', arr: galaxyData.map(g => g.dynamicML), note: 'dynamic M/L from Vflat^2*R/Mstar', circ: 'LOW' },
  { name: 'barGravD', arr: galaxyData.map(g => g.barGravDom), note: 'baryon-dominated fraction of RC', circ: 'MED' },
  { name: 'innBarF', arr: galaxyData.map(g => g.innerBarFrac), note: 'inner baryon-dominated fraction', circ: 'MED' },
  { name: 'gbarSprd', arr: galaxyData.map(g => g.gbarSpread), note: 'dynamic range of gbar', circ: 'LOW' },
  { name: 'mlColor', arr: galaxyData.map(g => g.mlColor), note: 'T-type based color/M-L proxy', circ: 'ZERO' },
  { name: 'SBdisk', arr: galaxyData.map(g => g.SBdisk), note: 'disk central surface brightness', circ: 'ZERO' },
  { name: 'logMstar', arr: galaxyData.map(g => g.logMstar), note: 'total stellar mass (from L*0.5)', circ: 'ZERO' },
];

log("  ┌──────────────┬──────┬────────┬────────┬────────┬────────┬────────┬────────┬──────────┬───────────┐");
log("  │ Variable     │ Circ │ |Cs r  │ Cs LOO │ Cs p   │ |Wk r  │ Wk LOO │ Wk p   │ r(Sig0)  │ Verdict   │");
log("  ├──────────────┼──────┼────────┼────────┼────────┼────────┼────────┼────────┼──────────┼───────────┤");

const NPERMS = 10000;
const mlResults = [];

for (const mv of mlVars) {
  const sdV = Math.sqrt(mv.arr.map(v => (v - mv.arr.reduce((s,v)=>s+v,0)/N)**2).reduce((s,v)=>s+v,0)/(N-1));
  if (sdV < 1e-8) {
    mlResults.push({ name: mv.name, circ: mv.circ, note: mv.note,
      cons: { r: 0, loo: gcCons.gc, permP: 1 }, work: { r: 0, loo: gcWork.gc, permP: 1 },
      rSigma0: 0, verdict: 'FAIL (no var)' });
    log("  │ " + mv.name.padEnd(12) + " │ " + mv.circ.padEnd(4) + " │    0.0 │  " +
      gcCons.gc.toFixed(1) + "% │  1.000 │    0.0 │  " + gcWork.gc.toFixed(1) + "% │  1.000 │    0.000 │ FAIL      │");
    continue;
  }

  const rxC = residualize(mv.arr, CONS_BL);
  const ryC = residualize(deltaA0, CONS_BL);
  const rC = corrWith(rxC, ryC);
  const gcC = looGapClosed([...CONS_BL, mv.arr], allLogA0, galaxyData, N);

  const rxW = residualize(mv.arr, WORK_BL);
  const ryW = residualize(deltaA0, WORK_BL);
  const rW = corrWith(rxW, ryW);
  const gcW = looGapClosed([...WORK_BL, mv.arr], allLogA0, galaxyData, N);

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

  const rSig = corrWith(mv.arr, sigma0Arr);

  let verdict;
  if (permW < 0.05 && gcW.gc > gcWork.gc + 1) verdict = 'SURVIVOR';
  else if (permC < 0.05 && gcC.gc > gcCons.gc + 1) verdict = 'PARTIAL';
  else verdict = 'FAIL';

  mlResults.push({
    name: mv.name, circ: mv.circ, note: mv.note,
    cons: { r: +rC.toFixed(3), loo: +gcC.gc.toFixed(1), permP: +permC.toFixed(4) },
    work: { r: +rW.toFixed(3), loo: +gcW.gc.toFixed(1), permP: +permW.toFixed(4) },
    rSigma0: +rSig.toFixed(3),
    verdict,
  });

  log("  │ " + mv.name.padEnd(12) + " │ " + mv.circ.padEnd(4) + " │ " +
    rC.toFixed(3).padStart(6) + " │ " + (gcC.gc.toFixed(1)+"%").padStart(6) + " │ " +
    permC.toFixed(3).padStart(6) + " │ " +
    rW.toFixed(3).padStart(6) + " │ " + (gcW.gc.toFixed(1)+"%").padStart(6) + " │ " +
    permW.toFixed(3).padStart(6) + " │ " + rSig.toFixed(3).padStart(8) + " │ " + verdict.padEnd(9) + " │");
}
log("  └──────────────┴──────┴────────┴────────┴────────┴────────┴────────┴────────┴──────────┴───────────┘");

const survivors = mlResults.filter(r => r.verdict === 'SURVIVOR');
const partials = mlResults.filter(r => r.verdict === 'PARTIAL');

log("");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  ABSORPTION TEST: Does Y* absorb Sigma0_bar?");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

const sigVars = [...survivors, ...partials];
if (sigVars.length > 0) {
  for (const sv of sigVars) {
    const svArr = mlVars.find(mv => mv.name === sv.name).arr;
    const noSig = [mhiArr, wigArr, envArr];
    const withML = [...noSig, svArr];

    const rxSig = residualize(sigma0Arr, withML);
    const rySig = residualize(deltaA0, withML);
    const sigAfter = corrWith(rxSig, rySig);

    const rxML = residualize(svArr, [...noSig, sigma0Arr]);
    const ryML = residualize(deltaA0, [...noSig, sigma0Arr]);
    const mlAfterSig = corrWith(rxML, ryML);

    log("  After controlling for " + sv.name + ":");
    log("    Sigma0_bar partial r: " + sigAfter.toFixed(3) + " (was 0.351 in baseline)");
    log("    " + sv.name + " partial r after Sigma0: " + mlAfterSig.toFixed(3));
    log("    Absorption of Sigma0: " + ((1 - Math.abs(sigAfter)/0.351)*100).toFixed(0) + "%");

    const gcReplace = looGapClosed([mhiArr, wigArr, envArr, svArr], allLogA0, galaxyData, N);
    log("    LOO with " + sv.name + " instead of Sigma0: " + gcReplace.gc.toFixed(1) + "% (was " + gcCons.gc.toFixed(1) + "%)");
    log("");
  }
} else {
  log("  No Y* proxy survived. Sigma0_bar cannot be absorbed.");
  log("  Sigma0_bar is NOT simply a stellar population effect.");
  log("  It captures something about baryon concentration that is");
  log("  orthogonal to M/L variations available from SPARC photometry.");
  log("");
}

log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  CRITICAL QUESTION: Could true SPS Y* explain the 61.9% gap?");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

const rSig0_delta = corrWith(sigma0Arr, deltaA0);
log("  Raw Sigma0-deltaA0 correlation: r = " + rSig0_delta.toFixed(3));
log("  Sigma0 baseline contribution: 0.351 (partial r in conservative BL)");
log("");
log("  The question is: does the 61.9% unexplained gap correlate with");
log("  any available stellar population proxy?");
log("");

const residBL = linReg(
  galaxyData.map((_, i) => WORK_BL.map(arr => arr[i])),
  allLogA0
).residuals;

log("  Correlation of EXTENDED BL RESIDUALS with Y*-related proxies:");
log("  ┌──────────────┬──────────┐");
log("  │ Proxy        │ r(resid) │");
log("  ├──────────────┼──────────┤");
for (const mv of mlVars) {
  const r = corrWith(mv.arr, residBL);
  log("  │ " + mv.name.padEnd(12) + " │ " + r.toFixed(3).padStart(8) + " │");
}
log("  └──────────────┴──────────┘");
log("");
log("  If all correlations are weak: the gap is NOT a Y* problem.");
log("  If some are moderate: true SPS data might help.");
log("");

log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  PHASE 59 — VERDICT");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

if (survivors.length > 0) {
  log("  RESULT: Y* proxies ADD beyond baseline.");
  log("  The gap is PARTLY a stellar population problem.");
} else if (partials.length > 0) {
  log("  RESULT: Some Y* proxies show partial signal, but not decisive.");
} else {
  log("  RESULT: NO Y* proxy adds independent information.");
  log("  INTERPRETATION:");
  log("    1. SPARC photometric proxies for Y* have been exhausted");
  log("    2. Sigma0_bar is NOT simply a Y* effect — it captures baryon geometry");
  log("    3. True SPS-derived Y* (from multi-band SED fitting) remains untested");
  log("    4. The 61.9% unexplained gap is NOT attributable to Y* variations");
  log("       accessible from single-band 3.6um photometry");
  log("    5. Resolution requires: IFU spectroscopy, multi-band imaging,");
  log("       or independent stellar mass estimates from dynamics");
}
log("");

const output = {
  version: VERSION, timestamp: new Date().toISOString(), nGalaxies: N,
  conservativeBL: +gcCons.gc.toFixed(1),
  extendedBL: +gcWork.gc.toFixed(1),
  stellarMLvars: mlResults,
  survivors: survivors.map(s => s.name),
  partials: partials.map(s => s.name),
  verdict: survivors.length > 0 ? 'ML_ADDS' : partials.length > 0 ? 'PARTIAL' : 'SIGMA0_IRREDUCIBLE',
};

fs.writeFileSync(path.join(__dirname, '../public/phase59-stellar-ml.json'), JSON.stringify(output, null, 2));
log("  Results saved to public/phase59-stellar-ml.json");
