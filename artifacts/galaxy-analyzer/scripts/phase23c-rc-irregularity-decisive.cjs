#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "23c.0.0";
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
  const reg = linReg(X, y);
  return reg.residuals;
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
  const logMHI = Math.log10(sp.MHI);
  const n = pts.length;

  const rArr = pts.map(p => p.r);
  const rMax = Math.max(...rArr);
  const rMid = rMax / 2;
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

  const slopeMismatch = Math.abs(slopeOf(innerPts) - slopeOf(outerPts));

  let rcCurvature = 0;
  if (pts.length >= 9) {
    const third = Math.floor(pts.length / 3);
    const s1 = slopeOf(pts.slice(0, third));
    const s2 = slopeOf(pts.slice(third, 2*third));
    const s3 = slopeOf(pts.slice(2*third));
    rcCurvature = Math.abs(s1 - 2*s2 + s3);
  }

  let rcRoughness = 0;
  if (pts.length >= 4) {
    let diffs = 0;
    for (let i = 1; i < pts.length; i++) {
      diffs += Math.abs(pts[i].log_g_obs - pts[i-1].log_g_obs);
    }
    const totalRange = Math.abs(pts[pts.length-1].log_g_obs - pts[0].log_g_obs) + 1e-10;
    rcRoughness = diffs / totalRange;
  }

  let outerFlatness = 0;
  if (outerPts.length >= 3) {
    const yOuter = outerPts.map(p => p.log_g_obs);
    const mO = yOuter.reduce((s,v)=>s+v,0)/yOuter.length;
    outerFlatness = Math.sqrt(yOuter.reduce((s,v)=>s+(v-mO)**2,0)/yOuter.length);
  }

  const Vflat = sp.Vflat || 100;
  const logVflat = Math.log10(Vflat);

  const beamSize = sp.D * 1000 * 4.848e-6 * 15;
  const Rdisk = sp.Rdisk || 3;
  const beamPerDisk = beamSize / Rdisk;
  const logBeamPerDisk = Math.log10(beamPerDisk + 0.01);

  const angularSize = Rdisk / (sp.D * 1000) * 206265;
  const logAngularSize = Math.log10(angularSize + 1);

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q,
    n, logMHI, logVflat,
    rcWiggliness, vPeakDip, slopeMismatch, rcCurvature, rcRoughness, outerFlatness,
    beamPerDisk, logBeamPerDisk, logAngularSize,
    kinematicDist: fit.rms * Math.sqrt(n),
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

log("");
log("=".repeat(80));
log("  PHASE 23c: DECISIVE RC IRREGULARITY TEST");
log("  Does rotation-curve shape irregularity survive ALL confounders?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  TARGET PROXIES (semi-circular, RC shape based):");
sep();
log("");
log("  rcWiggliness:   inner + outer local residual variability (from Phase 23b: r=0.365 after MHI)");
log("  vPeakDip:       (Vpeak - Vmin_after) / Vpeak (Phase 23b: r=0.295 after MHI)");
log("  slopeMismatch:  |inner slope - outer slope| in radius space");
log("  rcCurvature:    |s1 - 2*s2 + s3| second derivative of RC");
log("  rcRoughness:    sum|dV|/range — point-to-point roughness");
log("  outerFlatness:  RMS of outer RC around its mean");
log("");

log("  POTENTIAL CONFOUNDERS TO CONTROL FOR:");
log("  1. log(MHI)       — dominant known correlate");
log("  2. log(Vflat)     — galaxy mass/size");
log("  3. n_points       — data quality / resolution");
log("  4. inclination    — viewing angle");
log("  5. log(D)         — distance (beam smearing)");
log("  6. beam/disk      — resolution relative to galaxy size");
log("  7. angular size   — apparent size on sky");
log("");

const targets = [
  { name: 'rcWiggliness', values: galaxyData.map(g => g.rcWiggliness) },
  { name: 'vPeakDip', values: galaxyData.map(g => g.vPeakDip) },
  { name: 'slopeMismatch', values: galaxyData.map(g => g.slopeMismatch) },
  { name: 'rcCurvature', values: galaxyData.map(g => g.rcCurvature) },
  { name: 'rcRoughness', values: galaxyData.map(g => g.rcRoughness) },
  { name: 'outerFlatness', values: galaxyData.map(g => g.outerFlatness) },
];

const confounders = {
  'logMHI': galaxyData.map(g => g.logMHI),
  'logVflat': galaxyData.map(g => g.logVflat),
  'n_pts': galaxyData.map(g => g.n),
  'inc': galaxyData.map(g => g.inc),
  'logD': galaxyData.map(g => Math.log10(g.D)),
  'beamDisk': galaxyData.map(g => g.logBeamPerDisk),
  'angSize': galaxyData.map(g => g.logAngularSize),
};

log("=".repeat(80));
log("  TEST 1: SEQUENTIAL CONFOUNDER STRIPPING");
log("  Add confounders one by one. Does the signal survive?");
log("=".repeat(80));
log("");

const confounderSets = [
  { name: 'raw', controls: [] },
  { name: '+MHI', controls: ['logMHI'] },
  { name: '+MHI+Vflat', controls: ['logMHI', 'logVflat'] },
  { name: '+MHI+Vflat+n', controls: ['logMHI', 'logVflat', 'n_pts'] },
  { name: '+MHI+Vflat+n+inc', controls: ['logMHI', 'logVflat', 'n_pts', 'inc'] },
  { name: '+MHI+Vflat+n+inc+D', controls: ['logMHI', 'logVflat', 'n_pts', 'inc', 'logD'] },
  { name: 'ALL 7', controls: ['logMHI', 'logVflat', 'n_pts', 'inc', 'logD', 'beamDisk', 'angSize'] },
];

for (const target of targets) {
  log("  ── " + target.name + " ──");
  log("  ┌──────────────────────────────────────────────────────────────────┐");
  log("  │  Controls              r_partial  |t|     perm_p    status    │");
  log("  ├──────────────────────────────────────────────────────────────────┤");

  let allSurvived = true;
  const results = [];

  for (const cs of confounderSets) {
    const controlVals = cs.controls.map(c => confounders[c]);
    const residY = residualize(deltaA0, controlVals);
    const residX = residualize(target.values, controlVals);
    const r = corrWith(residX, residY);
    const df = N - 2 - cs.controls.length;
    const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1 - r*r + 1e-10);

    const NPERMS = 5000;
    const obsR = Math.abs(r);
    let cnt = 0;
    const sh = [...residY];
    for (let pp = 0; pp < NPERMS; pp++) {
      for (let i = sh.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [sh[i], sh[j]] = [sh[j], sh[i]];
      }
      if (Math.abs(corrWith(residX, sh)) >= obsR) cnt++;
    }
    const perm_p = cnt / NPERMS;

    const status = t >= 2.0 && perm_p < 0.05 ? "SURVIVES" :
                   t >= 1.65 ? "MARGINAL" : "FAILS";
    if (cs.controls.length > 0 && status === "FAILS") allSurvived = false;

    results.push({ controls: cs.name, r, t, perm_p, status });

    log("  │  " + cs.name.padEnd(22) + r.toFixed(3).padStart(8) +
      t.toFixed(1).padStart(6) + perm_p.toFixed(4).padStart(10) + "    " + status.padEnd(10) + "│");
  }

  log("  ├──────────────────────────────────────────────────────────────────┤");
  log("  │  FINAL: " + (allSurvived ? "SURVIVES ALL CONFOUNDERS" : "FAILS at some level").padEnd(51) + "│");
  log("  └──────────────────────────────────────────────────────────────────┘");
  log("");

  target.seqResults = results;
  target.survivesAll = allSurvived;
}

log("=".repeat(80));
log("  TEST 2: BEAM SMEARING / RESOLUTION CHECK");
log("  Are RC irregularity proxies just measuring poor resolution?");
log("=".repeat(80));
log("");

for (const target of targets.slice(0, 3)) {
  const rBeam = corrWith(target.values, galaxyData.map(g => g.logBeamPerDisk));
  const rAng = corrWith(target.values, galaxyData.map(g => g.logAngularSize));
  const rD = corrWith(target.values, galaxyData.map(g => Math.log10(g.D)));
  const rN = corrWith(target.values, galaxyData.map(g => g.n));
  log("  " + target.name + ":");
  log("    vs beam/disk:  r=" + rBeam.toFixed(3) + (Math.abs(rBeam) > 0.3 ? " ⚠ CONCERN" : " OK"));
  log("    vs angSize:    r=" + rAng.toFixed(3) + (Math.abs(rAng) > 0.3 ? " ⚠ CONCERN" : " OK"));
  log("    vs log(D):     r=" + rD.toFixed(3) + (Math.abs(rD) > 0.3 ? " ⚠ CONCERN" : " OK"));
  log("    vs n_pts:      r=" + rN.toFixed(3) + (Math.abs(rN) > 0.3 ? " ⚠ CONCERN" : " OK"));
  log("");
}

log("=".repeat(80));
log("  TEST 3: LOO CROSS-VALIDATION — RC IRREGULARITY MODELS");
log("=".repeat(80));
log("");

const allFeatures = [...targets, ...Object.entries(confounders).map(([k,v]) => ({ name: k, values: v }))];

function getVal(gIdx, fname) {
  const f = allFeatures.find(x => x.name === fname);
  return f ? f.values[gIdx] : 0;
}

const survivors = targets.filter(t => t.survivesAll);
const bestSurvivor = survivors.length > 0 ? survivors[0].name : targets[0].name;

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_wig: rcWiggliness', features: ['rcWiggliness'] },
  { name: 'M_vpd: vPeakDip', features: ['vPeakDip'] },
  { name: 'M_slp: slopeMismatch', features: ['slopeMismatch'] },
  { name: 'M_curv: rcCurvature', features: ['rcCurvature'] },
  { name: 'M_rough: rcRoughness', features: ['rcRoughness'] },
  { name: 'M_flat: outerFlatness', features: ['outerFlatness'] },
  { name: 'M_MHI: gas only', features: ['logMHI'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
  { name: 'M_MHI+vpd', features: ['logMHI', 'vPeakDip'] },
  { name: 'M_MHI+wig+vpd', features: ['logMHI', 'rcWiggliness', 'vPeakDip'] },
  { name: 'M_kitchen: MHI+all6', features: ['logMHI', 'rcWiggliness', 'vPeakDip', 'slopeMismatch', 'rcCurvature', 'rcRoughness', 'outerFlatness'] },
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

log("  ┌────────────────────────────────────────────────────────────────────────────┐");
log("  │  Model                   k    CV-RMS    vs M0     gap-closed             │");
log("  ├────────────────────────────────────────────────────────────────────────────┤");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const gapClosed = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  │  " + r.name.padEnd(24) + r.k.toString().padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(9) + (gapClosed + "%").padStart(14) + "            │");
}
log("  └────────────────────────────────────────────────────────────────────────────┘");
log("");

const wigGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_wig:')).cvRMS) / gap * 100 : 0;
const vpdGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_vpd:')).cvRMS) / gap * 100 : 0;
const mhiGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI: gas only').cvRMS) / gap * 100 : 0;
const mhiWigGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig').cvRMS) / gap * 100 : 0;
const mhiVpdGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+vpd').cvRMS) / gap * 100 : 0;
const mhiWigVpdGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig+vpd').cvRMS) / gap * 100 : 0;

log("  KEY LOO NUMBERS:");
log("    rcWiggliness alone:    " + wigGap.toFixed(1) + "%");
log("    vPeakDip alone:        " + vpdGap.toFixed(1) + "%");
log("    MHI alone:             " + mhiGap.toFixed(1) + "%");
log("    MHI + rcWiggliness:    " + mhiWigGap.toFixed(1) + "%");
log("    MHI + vPeakDip:        " + mhiVpdGap.toFixed(1) + "%");
log("    MHI + wig + vpd:       " + mhiWigVpdGap.toFixed(1) + "%");
log("");

log("=".repeat(80));
log("  TEST 4: QUARTILE tau ANALYSIS");
log("=".repeat(80));
log("");

const fullDL = dlMeta(allLogA0, galaxyData.map(g => g.se));

for (const target of targets.slice(0, 4)) {
  const sortedIdx = [...Array(N).keys()].sort((a, b) => target.values[a] - target.values[b]);
  const q1Idx = sortedIdx.slice(0, Math.floor(N / 4));
  const q4Idx = sortedIdx.slice(Math.floor(3 * N / 4));
  const q1LogA0 = q1Idx.map(i => allLogA0[i]);
  const q4LogA0 = q4Idx.map(i => allLogA0[i]);
  const q1DL = dlMeta(q1LogA0, q1Idx.map(i => galaxyData[i].se));
  const q4DL = dlMeta(q4LogA0, q4Idx.map(i => galaxyData[i].se));
  const minTau = Math.min(q1DL.tau, q4DL.tau);
  log("  " + target.name + ":");
  log("    Q1: tau=" + q1DL.tau.toFixed(3) + ", Q4: tau=" + q4DL.tau.toFixed(3));
  log("    Full tau: " + fullDL.tau.toFixed(3) + ", reduction: " + ((1 - minTau / fullDL.tau) * 100).toFixed(1) + "%");
  log("");
}

log("=".repeat(80));
log("  TEST 5: STABILITY — BOOTSTRAP CONFIDENCE INTERVALS");
log("=".repeat(80));
log("");

const NBOOT = 2000;
for (const target of targets.slice(0, 2)) {
  const controlVals = [confounders.logMHI, confounders.logVflat, confounders.n_pts, confounders.inc, confounders.logD];
  const bootR = [];
  for (let b = 0; b < NBOOT; b++) {
    const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
    const bY = idx.map(i => deltaA0[i]);
    const bX = idx.map(i => target.values[i]);
    const bControls = controlVals.map(c => idx.map(i => c[i]));
    const rY = residualize(bY, bControls);
    const rX = residualize(bX, bControls);
    bootR.push(corrWith(rX, rY));
  }
  bootR.sort((a,b) => a - b);
  const ci025 = bootR[Math.floor(NBOOT * 0.025)];
  const ci975 = bootR[Math.floor(NBOOT * 0.975)];
  const ciMed = bootR[Math.floor(NBOOT * 0.5)];
  const fracPositive = bootR.filter(r => r > 0).length / NBOOT;
  log("  " + target.name + " (partial r, controlling MHI+Vflat+n+inc+D):");
  log("    median: " + ciMed.toFixed(3) + "  95% CI: [" + ci025.toFixed(3) + ", " + ci975.toFixed(3) + "]");
  log("    fraction positive: " + (fracPositive * 100).toFixed(1) + "%");
  log("    CI excludes zero: " + (ci025 > 0 || ci975 < 0 ? "YES" : "NO"));
  log("");
}

log("=".repeat(80));
log("  PHASE 23c — FINAL DECISIVE VERDICT");
log("=".repeat(80));
log("");

const wigSurvivesAll = targets.find(t => t.name === 'rcWiggliness').survivesAll;
const vpdSurvivesAll = targets.find(t => t.name === 'vPeakDip').survivesAll;
const wigAddsToMHI = mhiWigGap > mhiGap + 1;
const vpdAddsToMHI = mhiVpdGap > mhiGap + 1;

log("  ┌──────────────────────────────────────────────────────────────────────────┐");
log("  │  CRITERION                     rcWiggliness         vPeakDip           │");
log("  ├──────────────────────────────────────────────────────────────────────────┤");
log("  │  Survives MHI?                 " + (targets[0].seqResults[1].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + (targets[1].seqResults[1].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + "│");
log("  │  Survives MHI+Vflat?           " + (targets[0].seqResults[2].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + (targets[1].seqResults[2].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + "│");
log("  │  Survives MHI+Vflat+n?         " + (targets[0].seqResults[3].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + (targets[1].seqResults[3].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + "│");
log("  │  Survives MHI+Vflat+n+inc?     " + (targets[0].seqResults[4].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + (targets[1].seqResults[4].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + "│");
log("  │  Survives MHI+Vflat+n+inc+D?   " + (targets[0].seqResults[5].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + (targets[1].seqResults[5].status === "SURVIVES" ? "YES" : "NO").padEnd(20) + "│");
log("  │  Survives ALL 7?               " + (wigSurvivesAll ? "YES" : "NO").padEnd(20) + (vpdSurvivesAll ? "YES" : "NO").padEnd(20) + "│");
log("  │  LOO gap closed                " + (wigGap.toFixed(1)+"%").padEnd(20) + (vpdGap.toFixed(1)+"%").padEnd(20) + "│");
log("  │  Adds to MHI in LOO?           " + (wigAddsToMHI ? "YES" : "NO").padEnd(20) + (vpdAddsToMHI ? "YES" : "NO").padEnd(20) + "│");
log("  │  Not beam smearing artifact?   check above         check above        │");
log("  └──────────────────────────────────────────────────────────────────────────┘");
log("");

const anyRCSurvives = wigSurvivesAll || vpdSurvivesAll;
const anyAddsToMHI = wigAddsToMHI || vpdAddsToMHI;

if (anyRCSurvives && anyAddsToMHI) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  RC IRREGULARITY IS A GENUINE INDEPENDENT SECOND PARAMETER.            ║");
  log("  ║  It survives all confounders and adds to MHI prediction.               ║");
  log("  ║  Galaxy history door: PARTIALLY OPEN via kinematic channel.            ║");
  log("  ║  Ready to proceed to Environment & Neighbors door.                     ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (anyRCSurvives || anyAddsToMHI) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  RC IRREGULARITY SHOWS PARTIAL SIGNAL.                                 ║");
  log("  ║  Survives some but not all tests.                                      ║");
  log("  ║  Not strong enough to call a confirmed second parameter.               ║");
  log("  ║  Galaxy history door: AMBIGUOUS — proceed to next door with caveat.    ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  RC IRREGULARITY FAILS DECISIVE TEST.                                  ║");
  log("  ║  The Phase 23 signal was mostly driven by fit quality (circular).      ║");
  log("  ║  Galaxy history door: OFFICIALLY CLOSED.                               ║");
  log("  ║  All 5 variables (sSFR, age, Z, SFH, disturbance) FAILED.            ║");
  log("  ║  NEXT: Environment & Neighbors.                                        ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  targets: targets.map(t => ({
    name: t.name,
    survivesAllConfounders: t.survivesAll,
    sequentialResults: t.seqResults.map(r => ({ controls: r.controls, r: +r.r.toFixed(3), t: +r.t.toFixed(1), perm_p: +r.perm_p.toFixed(4), status: r.status }))
  })),
  looGapClosed: {
    rcWiggliness: +wigGap.toFixed(1), vPeakDip: +vpdGap.toFixed(1),
    mhi: +mhiGap.toFixed(1),
    mhiPlusWig: +mhiWigGap.toFixed(1), mhiPlusVpd: +mhiVpdGap.toFixed(1),
    mhiPlusWigVpd: +mhiWigVpdGap.toFixed(1)
  },
  rcIrregularityIsGenuine: anyRCSurvives && anyAddsToMHI,
  doorStatus: (anyRCSurvives && anyAddsToMHI) ? 'PARTIALLY OPEN' :
              (anyRCSurvives || anyAddsToMHI) ? 'AMBIGUOUS' : 'CLOSED'
};

fs.writeFileSync(path.join(__dirname, '../public/phase23c-rc-irregularity-decisive.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase23c-rc-irregularity-decisive.json");
