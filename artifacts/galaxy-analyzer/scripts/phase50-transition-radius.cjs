#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "50.0.0";
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
  let runLen2 = 1;
  const residRun = [];
  for (let i = 1; i < rarResid.length; i++) {
    const s = Math.sign(rarResid[i]);
    if (s === currentSign) { runLen2++; }
    else { residRun.push(runLen2); runLen2 = 1; currentSign = s; }
  }
  residRun.push(runLen2);
  const meanRunLen = residRun.reduce((s,v)=>s+v,0) / residRun.length;
  const logMeanRun = Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1);

  const gbarArr = sorted.map(p => p.log_g_bar);
  const gobsArr = sorted.map(p => p.log_g_obs);
  const boostArr = sorted.map(p => p.log_g_obs - p.log_g_bar);

  let transR_gbar = rMax * 0.5;
  const a0_threshold = Math.log10(fit.a0);
  for (let i = 0; i < sorted.length; i++) {
    if (gbarArr[i] <= a0_threshold) {
      transR_gbar = sorted[i].r;
      break;
    }
  }
  const logTransR = Math.log10(transR_gbar > 0.1 ? transR_gbar : 0.1);
  const transR_norm = transR_gbar / Rdisk;
  const logTransRnorm = Math.log10(transR_norm > 0.01 ? transR_norm : 0.01);

  let transR_boost = rMax * 0.5;
  for (let i = 1; i < sorted.length; i++) {
    if (boostArr[i] > 0.3) {
      transR_boost = sorted[i].r;
      break;
    }
  }
  const logTransBoost = Math.log10(transR_boost > 0.1 ? transR_boost : 0.1);
  const transBoost_norm = transR_boost / Rdisk;
  const logTransBoostNorm = Math.log10(transBoost_norm > 0.01 ? transBoost_norm : 0.01);

  const fracBelow_a0 = gbarArr.filter(g => g <= a0_threshold).length / gbarArr.length;

  let peakR = sorted[0].r;
  let peakGobs = gobsArr[0];
  for (let i = 1; i < sorted.length; i++) {
    if (gobsArr[i] > peakGobs) { peakGobs = gobsArr[i]; peakR = sorted[i].r; }
  }
  const logPeakR = Math.log10(peakR > 0.1 ? peakR : 0.1);
  const peakR_norm = peakR / Rdisk;
  const logPeakRnorm = Math.log10(peakR_norm > 0.01 ? peakR_norm : 0.01);

  const Rflat_idx = (() => {
    if (sorted.length < 5) return sorted.length - 1;
    for (let i = Math.floor(sorted.length * 0.3); i < sorted.length - 2; i++) {
      const window = gobsArr.slice(i, Math.min(i + 3, sorted.length));
      const wMean = window.reduce((s,v)=>s+v,0)/window.length;
      const wVar = window.map(v=>(v-wMean)**2).reduce((s,v)=>s+v,0)/window.length;
      if (Math.sqrt(wVar) < 0.05) return i;
    }
    return sorted.length - 1;
  })();
  const Rflat = sorted[Rflat_idx].r;
  const logRflat = Math.log10(Rflat > 0.1 ? Rflat : 0.1);
  const Rflat_norm = Rflat / Rdisk;
  const logRflatNorm = Math.log10(Rflat_norm > 0.01 ? Rflat_norm : 0.01);

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, logVflat, rcWiggliness, envCode, logSigma0,
    logMeanRun,
    logTransR, logTransRnorm,
    logTransBoost, logTransBoostNorm,
    fracBelow_a0,
    logPeakR, logPeakRnorm,
    logRflat, logRflatNorm,
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
log("  PHASE 50: Internal Structure — Transition Radius");
log("  FAST PROTOCOL (dual baseline)");
log("  Working BL: MHI + wig + env + Sigma0_bar + meanRun");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  PROXIES:");
log("  1) transR: radius where g_bar drops below a0");
log("  2) transR/Rd: transition radius normalized by Rdisk");
log("  3) transBoost: radius where boost > 0.3 dex");
log("  4) transBst/Rd: boost transition normalized");
log("  5) fracBelow: fraction of RC points below a0");
log("  6) peakR: radius of peak g_obs");
log("  7) peakR/Rd: peak radius normalized");
log("  8) Rflat: radius where RC becomes flat");
log("  9) Rflat/Rd: flat radius normalized");
log("");

const candidates = [
  { name: 'transR', arr: galaxyData.map(g => g.logTransR), circ: 'HIGH', note: 'g_bar=a0 radius (uses a0!)' },
  { name: 'transR/Rd', arr: galaxyData.map(g => g.logTransRnorm), circ: 'HIGH', note: 'normalized (uses a0!)' },
  { name: 'transBst', arr: galaxyData.map(g => g.logTransBoost), circ: 'MED', note: 'boost>0.3 radius' },
  { name: 'tBst/Rd', arr: galaxyData.map(g => g.logTransBoostNorm), circ: 'MED', note: 'boost transition norm' },
  { name: 'fracBlow', arr: galaxyData.map(g => g.fracBelow_a0), circ: 'HIGH', note: 'frac below a0 (uses a0!)' },
  { name: 'peakR', arr: galaxyData.map(g => g.logPeakR), circ: 'LOW', note: 'peak g_obs radius' },
  { name: 'peakR/Rd', arr: galaxyData.map(g => g.logPeakRnorm), circ: 'LOW', note: 'peak radius normalized' },
  { name: 'Rflat', arr: galaxyData.map(g => g.logRflat), circ: 'LOW', note: 'RC flattening radius' },
  { name: 'Rflat/Rd', arr: galaxyData.map(g => g.logRflatNorm), circ: 'LOW', note: 'flat radius normalized' },
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
log("  CONSERVATIVE BL: " + consGC.gc.toFixed(1) + "%");
log("  WORKING BL: " + workGC.gc.toFixed(1) + "%");
log("");

log("  COLLINEARITY:");
for (const c of candidates) {
  const rMR = corrWith(c.arr, meanRunArr);
  const rMhi = corrWith(c.arr, mhiArr);
  const rVf = corrWith(c.arr, vflatArr);
  log("  " + c.name.padEnd(12) + " vs mRun:" + rMR.toFixed(2).padStart(6) +
    "  MHI:" + rMhi.toFixed(2).padStart(6) + "  Vf:" + rVf.toFixed(2).padStart(6));
}
log("");

log("=".repeat(80));
log("  FAST SCREEN");
log("=".repeat(80));
log("");

log("  ┌─────────────────────────────────────────────────────────────────────────────────────────────────────────┐");
log("  │  Proxy        Circ   |consBL t  consLOO  |workBL t  workLOO  Adds/w  Verdict                        │");
log("  ├─────────────────────────────────────────────────────────────────────────────────────────────────────────┤");

const results = [];
for (const cand of candidates) {
  const rRaw = corrWith(cand.arr, deltaA0);
  const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);

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

  const addsWork = gcWork.gc > workGC.gc + 1;

  const failFast = tRaw < 1.0 && tCons < 1.0 && tWork < 1.0;
  const promising = tWork >= 1.65 || addsWork;
  let verdict = failFast ? "FAIL-FAST" : promising ? "PROMISING" : "WEAK";
  if (cand.circ === 'HIGH' && promising) verdict = "PROMISING(circ!)";

  results.push({
    ...cand, rRaw, tRaw,
    rCons, tCons, gcCons: gcCons.gc,
    rWork, tWork, gcWork: gcWork.gc, addsWork,
    verdict,
  });

  log("  │  " + cand.name.padEnd(12) + cand.circ.padEnd(5) +
    tCons.toFixed(1).padStart(10) + (gcCons.gc.toFixed(1)+"%").padStart(8) +
    tWork.toFixed(1).padStart(10) + (gcWork.gc.toFixed(1)+"%").padStart(8) +
    (addsWork ? "  YES " : "  NO  ") +
    verdict.padStart(22) + "   │");
}
log("  └─────────────────────────────────────────────────────────────────────────────────────────────────────────┘");
log("");

const promising2 = results.filter(r => r.verdict.startsWith("PROMISING"));

if (promising2.length === 0) {
  log("  ══════════════════════════════════════════════════════════");
  log("  NO transition radius proxy passes fast screen.");
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
log("  PHASE 50 — SUMMARY");
log("=".repeat(80));
log("");
log("  Results summary:");
for (const r of results) {
  log("    " + r.name.padEnd(12) + " circ=" + r.circ.padEnd(5) +
    " |workBL t=" + r.tWork.toFixed(1).padStart(4) +
    " workLOO=" + r.gcWork.toFixed(1).padStart(5) + "% " + r.verdict);
}
log("");
log("  Working BL: " + workGC.gc.toFixed(1) + "%");
log("");
log("=".repeat(80));

const output = {
  version: VERSION, timestamp: new Date().toISOString(), nGalaxies: N,
  conservativeBL: { gc: +consGC.gc.toFixed(1) },
  workingBL: { gc: +workGC.gc.toFixed(1) },
  proxies: results.map(r => ({
    name: r.name, circ: r.circ, note: r.note,
    consBLt: +r.tCons.toFixed(1), consLOO: +r.gcCons.toFixed(1),
    workBLt: +r.tWork.toFixed(1), workLOO: +r.gcWork.toFixed(1),
    addsWork: r.addsWork, verdict: r.verdict,
  })),
};

fs.writeFileSync(path.join(__dirname, '../public/phase50-transition-radius.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase50-transition-radius.json");
