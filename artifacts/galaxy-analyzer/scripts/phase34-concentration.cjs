#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "34.0.0";
function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(80)); }

const G_pc = 4.302e-3;

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
  const R_last_kpc = rMax;
  const Mdyn_sun = (Vflat * 1e3) ** 2 * (R_last_kpc * 1e3 * 3.086e16) / (6.674e-11 * 1.989e30);
  const logMhalo = Math.log10(Mdyn_sun);

  const Rdisk = sp.Rdisk || sp.Reff || 3.0;

  const outerV = outerPts.map(p => Math.pow(10, p.log_g_obs) * p.r * 3.086e19);
  const innerV = innerPts.length > 2 ? innerPts.map(p => Math.pow(10, p.log_g_obs) * p.r * 3.086e19) : [];

  const Vmax = Vflat;
  const rHalf = rMax * 0.5;
  let Vhalf = Vflat;
  const ptsNearHalf = pts.filter(p => Math.abs(p.r - rHalf) < rHalf * 0.3);
  if (ptsNearHalf.length > 0) {
    const gObs = ptsNearHalf.map(p => Math.pow(10, p.log_g_obs));
    const rr = ptsNearHalf.map(p => p.r);
    const vv = gObs.map((g, i) => Math.sqrt(g * rr[i] * 3.086e19));
    Vhalf = vv.reduce((s,v) => s + v, 0) / vv.length / 1e3;
  }

  const cProxy = Vmax / (Vhalf > 0 ? Vhalf : Vmax);

  const rcShape = Vflat > 0 && innerPts.length > 2 ? (() => {
    const r2d = rMax * 0.2;
    const ptsAt2d = pts.filter(p => Math.abs(p.r - r2d) < r2d * 0.5);
    if (ptsAt2d.length === 0) return 1.0;
    const gObs2d = ptsAt2d.map(p => Math.pow(10, p.log_g_obs));
    const rr2d = ptsAt2d.map(p => p.r);
    const v2d = gObs2d.map((g, i) => Math.sqrt(g * rr2d[i] * 3.086e19)) ;
    const V2d = v2d.reduce((s,v) => s + v, 0) / v2d.length / 1e3;
    return Vflat / (V2d > 0 ? V2d : Vflat);
  })() : 1.0;

  const innerSlope = (() => {
    if (innerPts.length < 3) return 0;
    const xv = innerPts.map(p => Math.log10(p.r));
    const yv = innerPts.map(p => p.log_g_obs);
    const mx = xv.reduce((s,v)=>s+v,0)/xv.length;
    const my = yv.reduce((s,v)=>s+v,0)/yv.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < xv.length; i++) { sxy += (xv[i]-mx)*(yv[i]-my); sxx += (xv[i]-mx)**2; }
    return sxx > 0 ? sxy / sxx : 0;
  })();

  const RdiskToRlast = Rdisk / rMax;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms,
    n, logMHI, logVflat, rcWiggliness, envCode,
    logMhalo, logMbar, Rdisk, rMax,
    cProxy, rcShape, innerSlope, RdiskToRlast,
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
const cProxyArr = galaxyData.map(g => g.cProxy);
const rcShapeArr = galaxyData.map(g => g.rcShape);
const innerSlopeArr = galaxyData.map(g => g.innerSlope);
const rdiskRatioArr = galaxyData.map(g => g.RdiskToRlast);

const BASELINE = [mhiArr, wigArr, envArr];
const BASELINE_PLUS_MH = [mhiArr, wigArr, envArr, mhaloArr];

log("");
log("=".repeat(80));
log("  PHASE 34: HALO CONCENTRATION (c)");
log("  FAST PROTOCOL — Baseline: MHI + wig + env");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  CONCENTRATION PROXIES:");
log("  - cProxy = Vmax / V(r_half): rotation curve shape concentration");
log("  - rcShape = Vflat / V(0.2*Rmax): inner RC steepness");
log("  - innerSlope: d(log g_obs)/d(log r) in inner half");
log("  - Rdisk/Rlast: compactness ratio");
log("  All derived from rotation curve shape (no external data needed).");
log("");

log("  COLLINEARITY with existing variables:");
for (const [label, arr] of [['cProxy', cProxyArr], ['rcShape', rcShapeArr], ['innerSlope', innerSlopeArr], ['Rd/Rl', rdiskRatioArr]]) {
  const rMhi = corrWith(arr, mhiArr);
  const rWig = corrWith(arr, wigArr);
  const rMh = corrWith(arr, mhaloArr);
  log("  " + label.padEnd(14) + " vs MHI:" + rMhi.toFixed(2).padStart(6) + "  wig:" + rWig.toFixed(2).padStart(6) + "  Mh:" + rMh.toFixed(2).padStart(6));
}
log("");

log("=".repeat(80));
log("  STEP 1: FAST SCREEN — Raw + Baseline Control + LOO");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'envCode': envArr,
  'logMhalo': mhaloArr, 'cProxy': cProxyArr, 'rcShape': rcShapeArr,
  'innerSlope': innerSlopeArr, 'RdRl': rdiskRatioArr,
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
  return { gc: gap > 0 ? (m0rms - modelRms) / gap * 100 : 0, m0rms, m6rms, gap, modelRms };
}

const baselineGC = looGapClosed(['logMHI', 'rcWiggliness', 'envCode']);
const baselineMhGC = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', 'logMhalo']);

log("  BASELINE REFERENCE:");
log("  MHI+wig+env:       " + baselineGC.gc.toFixed(1) + "% gap closure");
log("  MHI+wig+env+Mh:    " + baselineMhGC.gc.toFixed(1) + "% gap closure");
log("");

const candidates = [
  { name: 'cProxy', arr: cProxyArr, features: 'cProxy' },
  { name: 'rcShape', arr: rcShapeArr, features: 'rcShape' },
  { name: 'innerSlope', arr: innerSlopeArr, features: 'innerSlope' },
  { name: 'Rdisk/Rlast', arr: rdiskRatioArr, features: 'RdRl' },
];

const results = [];

log("  ┌──────────────────────────────────────────────────────────────────────────────────┐");
log("  │  Proxy         Raw r    Raw t   |base r   |base t   LOO+base  Adds?  Verdict  │");
log("  ├──────────────────────────────────────────────────────────────────────────────────┤");

for (const cand of candidates) {
  const rRaw = corrWith(cand.arr, deltaA0);
  const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);

  const residX = residualize(cand.arr, BASELINE);
  const residY = residualize(deltaA0, BASELINE);
  const rBase = corrWith(residX, residY);
  const tBase = Math.abs(rBase) * Math.sqrt(N-5) / Math.sqrt(1-rBase**2+1e-10);

  const gc = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', cand.features]);
  const adds = gc.gc > baselineGC.gc + 1;

  const failFast = tRaw < 1.0 && tBase < 1.0;
  const promising = tBase >= 1.65 || adds;
  const verdict = failFast ? "FAIL-FAST" : promising ? "PROMISING" : "WEAK";

  results.push({ ...cand, rRaw, tRaw, rBase, tBase, gc: gc.gc, adds, verdict });

  log("  │  " + cand.name.padEnd(13) +
    rRaw.toFixed(3).padStart(7) + tRaw.toFixed(1).padStart(7) +
    rBase.toFixed(3).padStart(9) + tBase.toFixed(1).padStart(9) +
    (gc.gc.toFixed(1) + "%").padStart(10) +
    (adds ? "  YES " : "  NO  ") + verdict.padStart(10) + "  │");
}
log("  └──────────────────────────────────────────────────────────────────────────────────┘");
log("");

const promising = results.filter(r => r.verdict === "PROMISING");
const anyPromising = promising.length > 0;

if (!anyPromising) {
  log("  ══════════════════════════════════════════════════════════");
  log("  NO concentration proxy passes fast screen.");
  log("  All FAIL or WEAK. Skipping deep tests.");
  log("  ══════════════════════════════════════════════════════════");
  log("");
} else {
  log("  PROMISING candidates found: " + promising.map(p => p.name).join(", "));
  log("  Proceeding to DEEP TESTS...");
  log("");

  for (const cand of promising) {
    log("=".repeat(80));
    log("  DEEP TEST: " + cand.name);
    log("=".repeat(80));
    log("");

    log("  Confounder stripping:");
    for (const [label, ctrls] of [
      ['raw', []], ['|MHI', [mhiArr]], ['|MHI+wig', [mhiArr, wigArr]],
      ['|MHI+wig+env', BASELINE], ['|MHI+wig+env+Mh', BASELINE_PLUS_MH],
    ]) {
      const rx = residualize(cand.arr, ctrls);
      const ry = residualize(deltaA0, ctrls);
      const r = corrWith(rx, ry);
      const df = N - 2 - ctrls.length;
      const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
      log("    " + label.padEnd(22) + "r=" + r.toFixed(3) + " t=" + t.toFixed(1) +
        (t >= 2.0 ? " SURVIVES" : t >= 1.65 ? " MARGINAL" : " fails"));
    }
    log("");

    log("  Permutation test (10000 perms):");
    const NPERMS = 10000;
    for (const [label, ctrls] of [['|baseline', BASELINE], ['|baseline+Mh', BASELINE_PLUS_MH]]) {
      const rx = residualize(cand.arr, ctrls);
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
      log("    " + label.padEnd(18) + " |r|=" + obsR.toFixed(3) + " p=" + (cnt/NPERMS).toFixed(4) +
        (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "**" : ""));
    }
    log("");

    log("  LOO with expanded baseline:");
    const gcWithMh = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', 'logMhalo', cand.features]);
    log("    MHI+wig+env+Mh:        " + baselineMhGC.gc.toFixed(1) + "%");
    log("    MHI+wig+env+Mh+" + cand.name + ": " + gcWithMh.gc.toFixed(1) + "%");
    log("    Adds beyond Mh? " + (gcWithMh.gc > baselineMhGC.gc + 1 ? "YES (+" + (gcWithMh.gc - baselineMhGC.gc).toFixed(1) + "%)" : "NO"));
    log("");

    log("  Bootstrap CI (1000 resamples):");
    const NBOOT = 1000;
    const bootR = [];
    for (let b = 0; b < NBOOT; b++) {
      const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
      const bx = residualize(idx.map(i => cand.arr[i]), idx.map(i => BASELINE.map(c => c[i])).reduce((acc, row) => {
        row.forEach((v, j) => { if (!acc[j]) acc[j] = []; acc[j].push(v); });
        return acc;
      }, []));
      const by = residualize(idx.map(i => deltaA0[i]), idx.map(i => BASELINE.map(c => c[i])).reduce((acc, row) => {
        row.forEach((v, j) => { if (!acc[j]) acc[j] = []; acc[j].push(v); });
        return acc;
      }, []));
      bootR.push(corrWith(bx, by));
    }
    bootR.sort((a,b) => a - b);
    const lo = bootR[Math.floor(NBOOT * 0.025)];
    const hi = bootR[Math.floor(NBOOT * 0.975)];
    const crossesZero = lo <= 0 && hi >= 0;
    log("    95% CI: [" + lo.toFixed(3) + ", " + hi.toFixed(3) + "]" +
      (crossesZero ? " → crosses zero (WEAK)" : " → excludes zero (STRONG)"));
    log("");
  }
}

log("=".repeat(80));
log("  PHASE 34 — FAST PROTOCOL SUMMARY");
log("=".repeat(80));
log("");
log("  ┌──────────────────────────────────────────────────────────────────────┐");
log("  │  Test                        Result              Pass/Fail        │");
log("  ├──────────────────────────────────────────────────────────────────────┤");

const best = results.reduce((a,b) => Math.abs(b.tBase) > Math.abs(a.tBase) ? b : a, results[0]);

log("  │  Raw correlation             r=" + best.rRaw.toFixed(3) + " t=" + best.tRaw.toFixed(1) + "          " + (best.tRaw >= 2 ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  After baseline control      r=" + best.rBase.toFixed(3) + " t=" + best.tBase.toFixed(1) + "          " + (best.tBase >= 1.65 ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  LOO > baseline              " + best.gc.toFixed(1) + "% vs " + baselineGC.gc.toFixed(1) + "%" + "             " + (best.adds ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  Adds above baseline         " + (best.adds ? "YES" : "NO") + "                       " + (best.adds ? "PASS" : "FAIL").padEnd(12) + "│");

if (anyPromising) {
  const permRx = residualize(best.arr, BASELINE);
  const permRy = residualize(deltaA0, BASELINE);
  const obsR = Math.abs(corrWith(permRx, permRy));
  log("  │  Permutation p               see deep tests           " + "     ".padEnd(12) + "│");
  log("  │  Bootstrap CI                see deep tests           " + "     ".padEnd(12) + "│");
}
log("  ├──────────────────────────────────────────────────────────────────────┤");

const bestVerdict = best.verdict === "PROMISING" ? (best.adds && best.tBase >= 2.0 ? "CONFIRMED" : "PARTIAL") : "FAIL";
const nPass = [best.tRaw >= 2, best.tBase >= 1.65, best.adds, best.gc > baselineGC.gc + 1].filter(Boolean).length;

log("  │  FINAL VERDICT: " + bestVerdict.padEnd(12) + " (best proxy: " + best.name + ")" + " ".repeat(14) + "│");
log("  └──────────────────────────────────────────────────────────────────────┘");
log("");

log("  DARK HALO DOOR — SCORECARD:");
log("    M_halo:        4/6 PASS (suppressor, +18% LOO)");
log("    concentration: " + bestVerdict + " (best=" + best.name + ", " + nPass + "/4 fast)");
log("    f_DM:          next (or skip if absorbed)");
log("    halo spin:     next");
log("");
log("  Best model: MHI + rcWiggliness + envCode = 14.9% gap closure");
log("  Candidate: MHI + wig + env + logMhalo = " + baselineMhGC.gc.toFixed(1) + "% gap closure");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  baseline: { gc: +baselineGC.gc.toFixed(1) },
  baselinePlusMh: { gc: +baselineMhGC.gc.toFixed(1) },
  proxies: results.map(r => ({
    name: r.name, rawR: +r.rRaw.toFixed(3), rawT: +r.tRaw.toFixed(1),
    baseR: +r.rBase.toFixed(3), baseT: +r.tBase.toFixed(1),
    looGC: +r.gc.toFixed(1), adds: r.adds, verdict: r.verdict,
  })),
  bestProxy: best.name,
  finalVerdict: bestVerdict,
};

fs.writeFileSync(path.join(__dirname, '../public/phase34-concentration.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase34-concentration.json");
