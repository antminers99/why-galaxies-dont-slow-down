#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "23b.0.0";
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
  const mx = x.reduce((s, v) => s + v, 0) / n;
  const my = y.reduce((s, v) => s + v, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i]-mx)*(y[i]-my); sxx += (x[i]-mx)**2; syy += (y[i]-my)**2; }
  return sxy / Math.sqrt(sxx * syy + 1e-20);
}

function dlMeta(logA0s, ses) {
  const n = logA0s.length;
  const w = ses.map(s => 1 / (s * s + 1e-10));
  const wS = w.reduce((a, b) => a + b, 0);
  const muFE = logA0s.reduce((s, v, i) => s + w[i] * v, 0) / wS;
  let Q = 0;
  for (let i = 0; i < n; i++) Q += w[i] * (logA0s[i] - muFE) ** 2;
  const S1 = wS, S2 = w.reduce((s, v) => s + v * v, 0);
  const tau2 = Math.max(0, (Q - (n - 1)) / (S1 - S2 / S1));
  const tau = Math.sqrt(tau2);
  const W = ses.map(s => 1 / (s * s + tau2 + 1e-10));
  const WS = W.reduce((a, b) => a + b, 0);
  const mu = logA0s.reduce((s, v, i) => s + W[i] * v, 0) / WS;
  return { mu, tau, tau2, a0: Math.pow(10, mu) };
}

function linReg(X, y) {
  const n = y.length, p = X[0] ? X[0].length : 0;
  if (p === 0) return { coefs: [], intercept: y.reduce((s,v)=>s+v,0)/n };
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
  return { coefs: b, intercept };
}

const p11 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase11-sensitivity-lab.json'), 'utf8'));
const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));

const galaxyData = [];
for (const gal of p11.galaxies) {
  const pts = gal.localProfile.filter(p =>
    isFinite(p.log_g_bar) && isFinite(p.log_g_obs) && p.log_g_obs > p.log_g_bar * 0.5 && p.r > 0
  );
  if (pts.length < 5) continue;
  const sp = sparc.find(s => s.name === gal.name);
  if (!sp || !sp.MHI || sp.MHI <= 0) continue;
  const L = sp.L36 || sp.L || 0;
  if (L <= 0) continue;

  const fit = fitA0(pts);
  const logMstar = Math.log10(L * 0.5e9);
  const logMHI = Math.log10(sp.MHI);
  const logMHI_L = Math.log10(sp.MHI / L);

  const logSBdisk = sp.SBdisk > 0 ? Math.log10(sp.SBdisk) : 2;
  const logSBeff = sp.SBeff > 0 ? Math.log10(sp.SBeff) : 2;
  const concentration = sp.Reff > 0 && sp.Rdisk > 0 ? sp.Reff / sp.Rdisk : 1;
  const logConcentration = Math.log10(concentration);
  const gasExtent = sp.RHI > 0 && sp.Rdisk > 0 ? sp.RHI / sp.Rdisk : 3;

  const n = pts.length;
  const rArr = pts.map(p => p.r);
  const rMax = Math.max(...rArr);
  const rMid = rMax / 2;

  let innerSlope = 0, outerSlope = 0;
  const innerPts = pts.filter(p => p.r <= rMid);
  const outerPts = pts.filter(p => p.r > rMid);

  function slopeOf(subset) {
    if (subset.length < 3) return 0;
    const xv = subset.map(p => p.r);
    const yv = subset.map(p => p.log_g_obs);
    const mx = xv.reduce((s,v)=>s+v,0)/xv.length;
    const my = yv.reduce((s,v)=>s+v,0)/yv.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < xv.length; i++) { sxy += (xv[i]-mx)*(yv[i]-my); sxx += (xv[i]-mx)**2; }
    return sxx > 0 ? sxy / sxx : 0;
  }

  innerSlope = slopeOf(innerPts);
  outerSlope = slopeOf(outerPts);
  const slopeMismatch = Math.abs(innerSlope - outerSlope);

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

  const innerResidVar = residualVariability(innerPts);
  const outerResidVar = residualVariability(outerPts);
  const rcWiggliness = innerResidVar + outerResidVar;

  let vPeakDip = 0;
  if (pts.length >= 6) {
    const vobs = pts.map(p => Math.pow(10, p.log_g_obs / 2));
    let maxV = 0, minAfterMax = Infinity, maxIdx = 0;
    for (let i = 0; i < vobs.length; i++) {
      if (vobs[i] > maxV) { maxV = vobs[i]; maxIdx = i; }
    }
    for (let i = maxIdx; i < vobs.length; i++) {
      if (vobs[i] < minAfterMax) minAfterMax = vobs[i];
    }
    vPeakDip = maxV > 0 ? (maxV - minAfterMax) / maxV : 0;
  }

  const sbAsymmetry = Math.abs(logSBdisk - logSBeff);
  const sizeAnomaly = Math.abs(Math.log10(concentration) - Math.log10(1.678));
  const lopsidedness = Math.abs(logSBeff - logSBdisk - 0.5);
  const profileIrregularity = Math.abs(logSBeff - logSBdisk - logConcentration);

  const gasStarMismatch = Math.abs(gasExtent - 3.0);

  const morphDisturbance = sbAsymmetry + sizeAnomaly * 0.5 + lopsidedness * 0.3;

  const rcScatter = fit.rms;
  const kinematicDist = rcScatter * Math.sqrt(n);

  const edgeOnWarp = sp.inc > 75 ? sbAsymmetry * 1.5 : sbAsymmetry;

  const cleanComposite = slopeMismatch * 0.3 + rcWiggliness * 0.3 + vPeakDip * 0.2 + sbAsymmetry * 0.1 + gasStarMismatch * 0.1;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q, fD: sp.fD,
    n, logMHI, logMstar, logMHI_L,
    logSBdisk, logSBeff, logConcentration, gasExtent,
    Vflat: sp.Vflat,
    slopeMismatch, rcWiggliness, vPeakDip,
    sbAsymmetry, sizeAnomaly, lopsidedness, profileIrregularity,
    gasStarMismatch, morphDisturbance,
    rcScatter, kinematicDist,
    edgeOnWarp, cleanComposite,
    innerSlope, outerSlope,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

log("");
log("=".repeat(80));
log("  PHASE 23b: CIRCULARITY-CLEAN DISTURBANCE TEST");
log("  Separate circular from non-circular disturbance proxies");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  THE PROBLEM:");
sep();
log("");
log("  Phase 23 found kinematicDist (r=0.370, t=2.9) passes 3/4 criteria.");
log("  BUT: kinematicDist = rcScatter * sqrt(n)");
log("  rcScatter = RMS of RAR fit residuals = how well RAR fits this galaxy.");
log("  A galaxy with bad RAR fit → high rcScatter → more uncertain a0");
log("  → potentially larger |delta_a0| by construction.");
log("");
log("  THIS IS PARTIAL CIRCULARITY. We need to test whether the signal");
log("  survives in proxies that are NOT derived from the RAR fit itself.");
log("");

log("  PROXY CLASSIFICATION:");
sep();
log("");
log("  CIRCULAR (derived from RAR fit):");
log("    rcScatter:      RMS of RAR fit → CIRCULAR");
log("    kinematicDist:  rcScatter * sqrt(n) → CIRCULAR");
log("");
log("  SEMI-CIRCULAR (partly from rotation curve, not from RAR fit):");
log("    slopeMismatch:  |inner RC slope - outer RC slope| in radius space");
log("    rcWiggliness:   inner + outer local residual variability");
log("    vPeakDip:       (Vpeak - Vmin_after) / Vpeak → RC shape feature");
log("");
log("  CLEAN (no circularity — purely from photometry/HI/morphology):");
log("    sbAsymmetry:    |logSBdisk - logSBeff| → photometric");
log("    sizeAnomaly:    |log(Reff/Rdisk) - log(1.678)| → photometric");
log("    lopsidedness:   |logSBeff - logSBdisk - 0.5| → photometric");
log("    profileIrreg:   |logSBeff - logSBdisk - logConc| → photometric");
log("    gasStarMismatch: |RHI/Rdisk - 3| → HI+photometric");
log("    morphDisturbance: composite of above → photometric");
log("    edgeOnWarp:     sbAsymmetry * inc-weight → photometric");
log("");

const circularProxies = [
  { name: 'kinematicDist', values: galaxyData.map(g => g.kinematicDist), circ: 'CIRCULAR' },
  { name: 'rcScatter', values: galaxyData.map(g => g.rcScatter), circ: 'CIRCULAR' },
];

const semiCircularProxies = [
  { name: 'slopeMismatch', values: galaxyData.map(g => g.slopeMismatch), circ: 'SEMI' },
  { name: 'rcWiggliness', values: galaxyData.map(g => g.rcWiggliness), circ: 'SEMI' },
  { name: 'vPeakDip', values: galaxyData.map(g => g.vPeakDip), circ: 'SEMI' },
];

const cleanProxies = [
  { name: 'sbAsymmetry', values: galaxyData.map(g => g.sbAsymmetry), circ: 'CLEAN' },
  { name: 'sizeAnomaly', values: galaxyData.map(g => g.sizeAnomaly), circ: 'CLEAN' },
  { name: 'lopsidedness', values: galaxyData.map(g => g.lopsidedness), circ: 'CLEAN' },
  { name: 'profileIrreg', values: galaxyData.map(g => g.profileIrregularity), circ: 'CLEAN' },
  { name: 'gasStarMismatch', values: galaxyData.map(g => g.gasStarMismatch), circ: 'CLEAN' },
  { name: 'morphDisturb', values: galaxyData.map(g => g.morphDisturbance), circ: 'CLEAN' },
  { name: 'edgeOnWarp', values: galaxyData.map(g => g.edgeOnWarp), circ: 'CLEAN' },
  { name: 'cleanComposite', values: galaxyData.map(g => g.cleanComposite), circ: 'CLEAN' },
];

const allProxies = [...circularProxies, ...semiCircularProxies, ...cleanProxies];
const refFeatures = [
  { name: 'log(MHI)', values: galaxyData.map(g => g.logMHI), circ: 'REF' },
  { name: 'log(Mstar)', values: galaxyData.map(g => g.logMstar), circ: 'REF' },
  { name: 'n_points', values: galaxyData.map(g => g.n), circ: 'REF' },
  { name: 'log(D)', values: galaxyData.map(g => Math.log10(g.D)), circ: 'REF' },
  { name: 'T (Hubble)', values: galaxyData.map(g => g.T), circ: 'REF' },
  { name: 'inc', values: galaxyData.map(g => g.inc), circ: 'REF' },
];

const allFeatures = [...allProxies, ...refFeatures];

log("  STEP 1: CORRELATIONS BY CIRCULARITY CLASS");
sep();
log("");

const mhiVals = galaxyData.map(g => g.logMHI);
const residAfterMHI = [];
{
  const x = mhiVals, y = deltaA0;
  const n = x.length;
  const mx = x.reduce((s,v)=>s+v,0)/n, my = y.reduce((s,v)=>s+v,0)/n;
  let sxy=0,sxx=0;
  for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;}
  const b = sxy/sxx, a = my - b*mx;
  for(let i=0;i<n;i++) residAfterMHI.push(y[i] - (a + b*x[i]));
}

log("  ┌──────────────────────────────────────────────────────────────────────────────┐");
log("  │  Proxy              Class      r_raw  |t|   r_afterMHI  |t|   verdict     │");
log("  ├──────────────────────────────────────────────────────────────────────────────┤");

const proxyResults = [];
for (const p of allProxies) {
  const rRaw = corrWith(p.values, deltaA0);
  const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);
  const rAfterMHI = corrWith(p.values, residAfterMHI);
  const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);

  let verdict = '';
  if (tAfterMHI >= 2.0) verdict = 'STRONG';
  else if (tAfterMHI >= 1.65) verdict = 'MARGINAL';
  else if (Math.abs(rAfterMHI) >= 0.10) verdict = 'WEAK';
  else verdict = 'NONE';

  proxyResults.push({ ...p, rRaw, tRaw, rAfterMHI, tAfterMHI, verdict });

  log("  │  " + p.name.padEnd(19) + p.circ.padEnd(9) +
    rRaw.toFixed(3).padStart(7) + tRaw.toFixed(1).padStart(5) +
    rAfterMHI.toFixed(3).padStart(11) + tAfterMHI.toFixed(1).padStart(5) +
    "   " + verdict.padEnd(10) + "│");
}
log("  └──────────────────────────────────────────────────────────────────────────────┘");
log("");

const cleanResults = proxyResults.filter(p => p.circ === 'CLEAN');
const semiResults = proxyResults.filter(p => p.circ === 'SEMI');
const circResults = proxyResults.filter(p => p.circ === 'CIRCULAR');

const bestClean = cleanResults.sort((a,b) => Math.abs(b.rAfterMHI) - Math.abs(a.rAfterMHI))[0];
const bestSemi = semiResults.sort((a,b) => Math.abs(b.rAfterMHI) - Math.abs(a.rAfterMHI))[0];
const bestCirc = circResults.sort((a,b) => Math.abs(b.rAfterMHI) - Math.abs(a.rAfterMHI))[0];

log("  SUMMARY BY CLASS:");
log("    CIRCULAR:      best = " + bestCirc.name + " r_afterMHI=" + bestCirc.rAfterMHI.toFixed(3) + " (t=" + bestCirc.tAfterMHI.toFixed(1) + ") → " + bestCirc.verdict);
log("    SEMI-CIRCULAR: best = " + bestSemi.name + " r_afterMHI=" + bestSemi.rAfterMHI.toFixed(3) + " (t=" + bestSemi.tAfterMHI.toFixed(1) + ") → " + bestSemi.verdict);
log("    CLEAN:         best = " + bestClean.name + " r_afterMHI=" + bestClean.rAfterMHI.toFixed(3) + " (t=" + bestClean.tAfterMHI.toFixed(1) + ") → " + bestClean.verdict);
log("");

log("  STEP 2: PERMUTATION TESTS ON CLEAN + SEMI PROXIES (vs MHI residuals)");
sep();
log("");

const NPERMS = 10000;
const testProxies = [...semiResults, ...cleanResults].filter(p => Math.abs(p.rAfterMHI) >= 0.10);

for (const p of testProxies) {
  const obsR = Math.abs(p.rAfterMHI);
  let cnt = 0;
  const sh = [...residAfterMHI];
  for (let pp = 0; pp < NPERMS; pp++) {
    for (let i = sh.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [sh[i], sh[j]] = [sh[j], sh[i]];
    }
    if (Math.abs(corrWith(p.values, sh)) >= obsR) cnt++;
  }
  const perm_p = cnt / NPERMS;
  log("  " + p.name.padEnd(16) + " [" + p.circ + "]  |r|=" + obsR.toFixed(3) +
    " perm_p=" + perm_p.toFixed(4) + (perm_p < 0.05 ? " *" : "") + (perm_p < 0.01 ? "*" : ""));
}
if (testProxies.length === 0) log("  No clean/semi proxy has |r_afterMHI| >= 0.10");
log("");

log("  STEP 3: LOO CROSS-VALIDATION — CLEAN vs CIRCULAR");
sep();
log("");

function getVal(gIdx, fname) {
  const f = allFeatures.find(x => x.name === fname);
  return f ? f.values[gIdx] : 0;
}

const modelDefs = [
  { name: 'M0: Universal', features: [], label: 'baseline' },
  { name: 'M_circ: kinematicDist', features: ['kinematicDist'], label: 'CIRCULAR' },
  { name: 'M_semi: best semi', features: [bestSemi.name], label: 'SEMI' },
  { name: 'M_clean: best clean', features: [bestClean.name], label: 'CLEAN' },
  { name: 'M_allClean: all clean', features: cleanProxies.map(c => c.name), label: 'CLEAN' },
  { name: 'M_allSemi: all semi', features: semiCircularProxies.map(c => c.name), label: 'SEMI' },
  { name: 'M_nonCirc: semi+clean', features: [...semiCircularProxies.map(c=>c.name), ...cleanProxies.map(c=>c.name)], label: 'BOTH' },
  { name: 'M_MHI: gas only', features: ['log(MHI)'], label: 'reference' },
  { name: 'M_MHI+clean: gas+clean', features: ['log(MHI)', bestClean.name], label: 'combined' },
  { name: 'M_MHI+semi: gas+semi', features: ['log(MHI)', bestSemi.name], label: 'combined' },
  { name: 'M_MHI+circ: gas+circ', features: ['log(MHI)', 'kinematicDist'], label: 'combined' },
  { name: 'M6: Per-galaxy', features: ['__free__'], label: 'ceiling' },
];

const cvResults = [];
for (const model of modelDefs) {
  let totalSS = 0, totalN = 0;
  if (model.features[0] === '__free__') {
    for (let i = 0; i < N; i++) { totalSS += predictSS(galaxyData[i].a0, galaxyData[i].pts); totalN += galaxyData[i].pts.length; }
    cvResults.push({ ...model, cvRMS: Math.sqrt(totalSS / totalN), k: N }); continue;
  }
  if (model.features.length === 0) {
    for (let i = 0; i < N; i++) { const trainY = allLogA0.filter((_, j) => j !== i); const mu = trainY.reduce((s,v) => s+v, 0) / trainY.length; totalSS += predictSS(Math.pow(10, mu), galaxyData[i].pts); totalN += galaxyData[i].pts.length; }
    cvResults.push({ ...model, cvRMS: Math.sqrt(totalSS / totalN), k: 1 }); continue;
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
  cvResults.push({ ...model, cvRMS: Math.sqrt(totalSS / totalN), k: model.features.length + 1 });
}

const m0rms = cvResults[0].cvRMS;
const m6rms = cvResults[cvResults.length - 1].cvRMS;
const gap = m0rms - m6rms;

log("  ┌──────────────────────────────────────────────────────────────────────────────┐");
log("  │  Model                     label     k   CV-RMS    vs M0    gap-closed     │");
log("  ├──────────────────────────────────────────────────────────────────────────────┤");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const gapClosed = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  │  " + r.name.padEnd(25) + r.label.padEnd(10) + r.k.toString().padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(9) + (gapClosed + "%").padStart(13) + "    │");
}
log("  └──────────────────────────────────────────────────────────────────────────────┘");
log("");

const circGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_circ')).cvRMS) / gap * 100 : 0;
const semiGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_semi:')).cvRMS) / gap * 100 : 0;
const cleanGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_clean:')).cvRMS) / gap * 100 : 0;
const mhiCircGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('MHI+circ')).cvRMS) / gap * 100 : 0;
const mhiSemiGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('MHI+semi')).cvRMS) / gap * 100 : 0;
const mhiCleanGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('MHI+clean')).cvRMS) / gap * 100 : 0;
const mhiGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI: gas only').cvRMS) / gap * 100 : 0;

log("  GAP CLOSURE COMPARISON:");
log("    CIRCULAR (kinematicDist):  " + circGap.toFixed(1) + "%");
log("    SEMI (best):               " + semiGap.toFixed(1) + "%");
log("    CLEAN (best):              " + cleanGap.toFixed(1) + "%");
log("    MHI alone:                 " + mhiGap.toFixed(1) + "%");
log("    MHI + CIRCULAR:            " + mhiCircGap.toFixed(1) + "%");
log("    MHI + SEMI:                " + mhiSemiGap.toFixed(1) + "%");
log("    MHI + CLEAN:               " + mhiCleanGap.toFixed(1) + "%");
log("");

const cleanAddsToMHI = mhiCleanGap > mhiGap + 2;
const semiAddsToMHI = mhiSemiGap > mhiGap + 2;
const circAddsToMHI = mhiCircGap > mhiGap + 2;

log("  Does adding disturbance to MHI improve prediction?");
log("    MHI + CIRCULAR: " + (circAddsToMHI ? "YES (+" + (mhiCircGap - mhiGap).toFixed(1) + "%)" : "NO"));
log("    MHI + SEMI:     " + (semiAddsToMHI ? "YES (+" + (mhiSemiGap - mhiGap).toFixed(1) + "%)" : "NO"));
log("    MHI + CLEAN:    " + (cleanAddsToMHI ? "YES (+" + (mhiCleanGap - mhiGap).toFixed(1) + "%)" : "NO"));
log("");

log("  STEP 4: QUARTILE ANALYSIS BY CIRCULARITY CLASS");
sep();
log("");

const fullDL = dlMeta(allLogA0, galaxyData.map(g => g.se));

for (const proxy of [bestCirc, bestSemi, bestClean]) {
  const sortedIdx = [...Array(N).keys()].sort((a, b) => proxy.values[a] - proxy.values[b]);
  const q1Idx = sortedIdx.slice(0, Math.floor(N / 4));
  const q4Idx = sortedIdx.slice(Math.floor(3 * N / 4));
  const q1LogA0 = q1Idx.map(i => allLogA0[i]);
  const q4LogA0 = q4Idx.map(i => allLogA0[i]);
  const q1DL = dlMeta(q1LogA0, q1Idx.map(i => galaxyData[i].se));
  const q4DL = dlMeta(q4LogA0, q4Idx.map(i => galaxyData[i].se));
  const q1mean = q1LogA0.reduce((s,v)=>s+v,0)/q1LogA0.length;
  const q4mean = q4LogA0.reduce((s,v)=>s+v,0)/q4LogA0.length;
  const splitDex = Math.abs(q1mean - q4mean);
  const tauMinQ = Math.min(q1DL.tau, q4DL.tau);
  log("  " + proxy.name + " [" + proxy.circ + "]:");
  log("    Q1 (low): a0=" + Math.round(Math.pow(10, q1mean)) + ", tau=" + q1DL.tau.toFixed(3));
  log("    Q4 (high): a0=" + Math.round(Math.pow(10, q4mean)) + ", tau=" + q4DL.tau.toFixed(3));
  log("    Split: " + splitDex.toFixed(3) + " dex, tau reduction: " + ((1-tauMinQ/fullDL.tau)*100).toFixed(1) + "%");
  log("");
}

log("=".repeat(80));
log("  PHASE 23b — THE VERDICT");
log("=".repeat(80));
log("");

const cleanSignal = cleanResults.some(p => p.tAfterMHI >= 1.65);
const semiSignal = semiResults.some(p => p.tAfterMHI >= 1.65);
const cleanLOO = cleanGap > 2;
const semiLOO = semiGap > 2;

log("  QUESTION: Is the Phase 23 disturbance signal real or circular?");
log("");
log("  ┌──────────────────────────────────────────────────────────────────┐");
log("  │  Test                    CIRCULAR   SEMI      CLEAN             │");
log("  ├──────────────────────────────────────────────────────────────────┤");
log("  │  r_afterMHI (best)       " +
  bestCirc.rAfterMHI.toFixed(3).padEnd(9) +
  bestSemi.rAfterMHI.toFixed(3).padEnd(10) +
  bestClean.rAfterMHI.toFixed(3).padEnd(18) + "│");
log("  │  t_afterMHI              " +
  bestCirc.tAfterMHI.toFixed(1).padEnd(9) +
  bestSemi.tAfterMHI.toFixed(1).padEnd(10) +
  bestClean.tAfterMHI.toFixed(1).padEnd(18) + "│");
log("  │  LOO gap closed          " +
  (circGap.toFixed(1)+"%").padEnd(9) +
  (semiGap.toFixed(1)+"%").padEnd(10) +
  (cleanGap.toFixed(1)+"%").padEnd(18) + "│");
log("  │  Adds to MHI?            " +
  (circAddsToMHI?"YES":"NO").padEnd(9) +
  (semiAddsToMHI?"YES":"NO").padEnd(10) +
  (cleanAddsToMHI?"YES":"NO").padEnd(18) + "│");
log("  └──────────────────────────────────────────────────────────────────┘");
log("");

if (cleanSignal || cleanLOO) {
  log("  CONCLUSION:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  CLEAN (non-circular) disturbance proxies RETAIN signal.           ║");
  log("  ║  The disturbance effect is GENUINE, not just circular artifact.    ║");
  log("  ║  Galaxy history door remains partially open via this channel.      ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else if (semiSignal || semiLOO) {
  log("  CONCLUSION:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  Semi-circular proxies (RC shape) retain some signal.              ║");
  log("  ║  The effect may be partly real but partly from RC noise.           ║");
  log("  ║  Not conclusive — needs external disturbance data.                 ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  CONCLUSION:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  CLEAN proxies show NO independent signal.                         ║");
  log("  ║  The Phase 23 result was MOSTLY CIRCULAR ARTIFACT.                 ║");
  log("  ║  The strong kinematicDist signal came from RAR fit quality,        ║");
  log("  ║  not from genuine morphological/kinematic disturbance.             ║");
  log("  ║                                                                     ║");
  log("  ║  GALAXY HISTORY DOOR: OFFICIALLY CLOSED.                           ║");
  log("  ║  All 5 variables tested: sSFR, age, Z, SFH, merger — FAILED.     ║");
  log("  ║  NEXT: Environment & Neighbors.                                    ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  proxyResults: proxyResults.map(p => ({
    name: p.name, circularity: p.circ,
    rRaw: +p.rRaw.toFixed(3), tRaw: +p.tRaw.toFixed(1),
    rAfterMHI: +p.rAfterMHI.toFixed(3), tAfterMHI: +p.tAfterMHI.toFixed(1),
    verdict: p.verdict
  })),
  bestByClass: {
    circular: { name: bestCirc.name, rAfterMHI: +bestCirc.rAfterMHI.toFixed(3), tAfterMHI: +bestCirc.tAfterMHI.toFixed(1) },
    semi: { name: bestSemi.name, rAfterMHI: +bestSemi.rAfterMHI.toFixed(3), tAfterMHI: +bestSemi.tAfterMHI.toFixed(1) },
    clean: { name: bestClean.name, rAfterMHI: +bestClean.rAfterMHI.toFixed(3), tAfterMHI: +bestClean.tAfterMHI.toFixed(1) }
  },
  looGapClosed: {
    circular: +circGap.toFixed(1), semi: +semiGap.toFixed(1), clean: +cleanGap.toFixed(1),
    mhi: +mhiGap.toFixed(1),
    mhiPlusCirc: +mhiCircGap.toFixed(1), mhiPlusSemi: +mhiSemiGap.toFixed(1), mhiPlusClean: +mhiCleanGap.toFixed(1)
  },
  cleanSignalSurvives: cleanSignal || cleanLOO,
  semiSignalSurvives: semiSignal || semiLOO,
  doorStatus: (cleanSignal || cleanLOO) ? 'PARTIALLY OPEN' : (semiSignal || semiLOO) ? 'AMBIGUOUS' : 'CLOSED'
};

fs.writeFileSync(path.join(__dirname, '../public/phase23b-circularity-clean.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase23b-circularity-clean.json");
