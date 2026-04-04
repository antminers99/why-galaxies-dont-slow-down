#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "38.0.0";
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
  const logMbar = Math.log10(Mbar_sun);
  const Sigma0_bar = Mbar_sun / (Math.PI * (Rdisk * 1e3) ** 2);
  const logSigma0 = Math.log10(Sigma0_bar > 0 ? Sigma0_bar : 1e-3);

  const Mdyn_sun = (Vflat * 1e3) ** 2 * (rMax * 1e3 * 3.086e16) / (6.674e-11 * 1.989e30);
  const logMhalo = Math.log10(Mdyn_sun);

  const fDM = Math.max(0.01, Math.min(0.999, 1 - Mbar_sun / Mdyn_sun));
  const logFDM = Math.log10(fDM);
  const MhaloToMbar = logMhalo - logMbar;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, logVflat, rcWiggliness, envCode, logSigma0,
    fDM, logFDM, logMhalo, logMbar, MhaloToMbar,
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
const fDMArr = galaxyData.map(g => g.fDM);
const logFDMArr = galaxyData.map(g => g.logFDM);
const mhaloArr = galaxyData.map(g => g.logMhalo);
const mhaloToMbarArr = galaxyData.map(g => g.MhaloToMbar);

const OLD_BASELINE = [mhiArr, wigArr, envArr];
const NEW_BASELINE = [mhiArr, wigArr, envArr, sigma0Arr];

log("");
log("=".repeat(80));
log("  PHASE 38: f_DM (Dark Matter Fraction)");
log("  FAST PROTOCOL — NEW Baseline: MHI + wig + env + Sigma0_bar");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  BASELINE UPDATE:");
log("  OLD: MHI + wig + env = 14.9%");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'envCode': envArr,
  'logSigma0': sigma0Arr, 'logVflat': vflatArr,
  'fDM': fDMArr, 'logFDM': logFDMArr, 'logMhalo': mhaloArr,
  'MhaloToMbar': mhaloToMbarArr,
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

const oldBLgc = looGapClosed(['logMHI', 'rcWiggliness', 'envCode']);
const newBLgc = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', 'logSigma0']);
log("  NEW: MHI + wig + env + Sigma0_bar = " + newBLgc.gc.toFixed(1) + "%");
log("");

log("  f_DM PROXIES:");
log("  - fDM: 1 - Mbar/Mdyn (dark matter fraction, linear)");
log("  - log(fDM): log10 of above");
log("  - Mhalo/Mbar: ratio of total to baryonic mass");
log("");

log("  CIRCULARITY: HIGH");
log("  fDM uses Mdyn = Vflat^2 * Rlast / G (from rotation curve).");
log("  Same circularity concern as M_halo.");
log("");

const candidates = [
  { name: 'fDM', arr: fDMArr, feat: 'fDM', circ: 'HIGH' },
  { name: 'log(fDM)', arr: logFDMArr, feat: 'logFDM', circ: 'HIGH' },
  { name: 'Mh/Mb', arr: mhaloToMbarArr, feat: 'MhaloToMbar', circ: 'HIGH' },
];

log("  COLLINEARITY with NEW baseline vars:");
for (const c of candidates) {
  const rMhi = corrWith(c.arr, mhiArr);
  const rWig = corrWith(c.arr, wigArr);
  const rSig = corrWith(c.arr, sigma0Arr);
  const rVf = corrWith(c.arr, vflatArr);
  log("  " + c.name.padEnd(12) + " vs MHI:" + rMhi.toFixed(2).padStart(6) + "  wig:" + rWig.toFixed(2).padStart(6) +
    "  Sig0:" + rSig.toFixed(2).padStart(6) + "  Vf:" + rVf.toFixed(2).padStart(6));
}
log("");

log("=".repeat(80));
log("  STEP 1: FAST SCREEN — against NEW baseline");
log("=".repeat(80));
log("");

log("  ┌──────────────────────────────────────────────────────────────────────────────────────────┐");
log("  │  Proxy       Circ     Raw r    Raw t   |newBL r  |newBL t  LOO+newBL Adds?  Verdict   │");
log("  ├──────────────────────────────────────────────────────────────────────────────────────────┤");

const results = [];
for (const cand of candidates) {
  const rRaw = corrWith(cand.arr, deltaA0);
  const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);

  const rx = residualize(cand.arr, NEW_BASELINE);
  const ry = residualize(deltaA0, NEW_BASELINE);
  const rNewBL = corrWith(rx, ry);
  const tNewBL = Math.abs(rNewBL) * Math.sqrt(N-6) / Math.sqrt(1-rNewBL**2+1e-10);

  const gc = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', 'logSigma0', cand.feat]);
  const adds = gc.gc > newBLgc.gc + 1;

  const failFast = tRaw < 1.0 && tNewBL < 1.0;
  const promising = tNewBL >= 1.65 || adds;
  const verdict = failFast ? "FAIL-FAST" : promising ? "PROMISING" : "WEAK";

  results.push({ ...cand, rRaw, tRaw, rNewBL, tNewBL, gc: gc.gc, adds, verdict });

  log("  │  " + cand.name.padEnd(11) + cand.circ.padEnd(8) +
    rRaw.toFixed(3).padStart(7) + tRaw.toFixed(1).padStart(7) +
    rNewBL.toFixed(3).padStart(9) + tNewBL.toFixed(1).padStart(9) +
    (gc.gc.toFixed(1) + "%").padStart(10) +
    (adds ? "  YES " : "  NO  ") + verdict.padStart(10) + "   │");
}
log("  └──────────────────────────────────────────────────────────────────────────────────────────┘");
log("");

log("  ALSO: Test against OLD baseline for comparison:");
for (const cand of candidates) {
  const rx = residualize(cand.arr, OLD_BASELINE);
  const ry = residualize(deltaA0, OLD_BASELINE);
  const r = corrWith(rx, ry);
  const t = Math.abs(r) * Math.sqrt(N-5) / Math.sqrt(1-r**2+1e-10);
  const gc = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', cand.feat]);
  log("  " + cand.name.padEnd(12) + " |oldBL: r=" + r.toFixed(3) + " t=" + t.toFixed(1) + " LOO=" + gc.gc.toFixed(1) + "%");
}
log("");

const promising = results.filter(r => r.verdict === "PROMISING");

if (promising.length === 0) {
  log("  ══════════════════════════════════════════════════════════");
  log("  NO f_DM proxy passes fast screen against NEW baseline.");
  log("  Skipping deep tests.");
  log("  ══════════════════════════════════════════════════════════");
} else {
  log("  PROMISING: " + promising.map(p => p.name).join(", "));
  log("  Proceeding to deep tests...");
  log("");

  for (const cand of promising) {
    log("=".repeat(80));
    log("  DEEP TEST: " + cand.name);
    log("=".repeat(80));
    log("");

    log("  Confounder stripping:");
    for (const [label, ctrls] of [
      ['raw', []], ['|MHI', [mhiArr]],
      ['|oldBL', OLD_BASELINE],
      ['|newBL', NEW_BASELINE],
      ['|newBL+Vflat', [...NEW_BASELINE, vflatArr]],
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
    for (const [label, ctrls] of [['|newBL', NEW_BASELINE], ['|newBL+Vflat', [...NEW_BASELINE, vflatArr]]]) {
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
      const subBL = NEW_BASELINE.map(arr => idx.map(j => arr[j]));
      const rx2 = residualize(subX, subBL);
      const ry2 = residualize(subY, subBL);
      jackR.push(corrWith(rx2, ry2));
    }
    const jackMean = jackR.reduce((s,v)=>s+v,0)/N;
    const jackSD = Math.sqrt(jackR.map(v => (v - jackMean)**2).reduce((s,v)=>s+v,0)/(N-1));
    const signFlips = jackR.filter(r => Math.sign(r) !== Math.sign(cand.rNewBL)).length;
    log("    Mean=" + jackMean.toFixed(3) + " SD=" + jackSD.toFixed(3) + " flips=" + signFlips + "/" + N);
    log("");

    log("  Bootstrap CI (2000):");
    const NBOOT = 2000;
    const bootR = [];
    for (let b = 0; b < NBOOT; b++) {
      const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
      const bX = idx.map(i => cand.arr[i]);
      const bY = idx.map(i => deltaA0[i]);
      const bBL = NEW_BASELINE.map(arr => idx.map(i => arr[i]));
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
log("  PHASE 38 — FAST PROTOCOL SUMMARY");
log("=".repeat(80));
log("");

const best = results.reduce((a,b) => Math.abs(b.tNewBL) > Math.abs(a.tNewBL) ? b : a, results[0]);
const nFP = [best.tRaw >= 2, best.tNewBL >= 1.65, best.adds, best.gc > newBLgc.gc + 1].filter(Boolean).length;

log("  ┌──────────────────────────────────────────────────────────────────────┐");
log("  │  Test                        Result              Pass/Fail        │");
log("  ├──────────────────────────────────────────────────────────────────────┤");
log("  │  Raw correlation             r=" + best.rRaw.toFixed(3) + " t=" + best.tRaw.toFixed(1) + "          " + (best.tRaw >= 2 ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  After NEW baseline control  r=" + best.rNewBL.toFixed(3) + " t=" + best.tNewBL.toFixed(1) + "          " + (best.tNewBL >= 1.65 ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  LOO > NEW baseline          " + best.gc.toFixed(1) + "% vs " + newBLgc.gc.toFixed(1) + "%" + "          " + (best.adds ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  Adds above NEW baseline     " + (best.adds ? "YES" : "NO").padEnd(24) + (best.adds ? "PASS" : "FAIL").padEnd(12) + "│");
log("  ├──────────────────────────────────────────────────────────────────────┤");

const verdict = nFP >= 3 ? "CONFIRMED" : nFP >= 2 ? "PARTIAL" : "FAIL";
log("  │  FINAL VERDICT: " + verdict.padEnd(50) + "│");
log("  └──────────────────────────────────────────────────────────────────────┘");
log("");

log("  DARK HALO DOOR — SCORECARD:");
log("    M_halo/Mdyn:   PARTIAL (= V_flat)");
log("    concentration: FAIL");
log("    V_flat:        PASS but CIRCULAR");
log("    R_s:           PARTIAL (Rs_Vrise)");
log("    rho_0:         Sigma0_bar CONFIRMED (reclassified: baryon structure)");
log("    f_DM:          " + verdict);
log("    Next: f_DM at reference radius");
log("");
log("  Baseline: MHI + wig + env + Sigma0_bar = " + newBLgc.gc.toFixed(1) + "%");
log("");
log("=".repeat(80));

const output = {
  version: VERSION, timestamp: new Date().toISOString(), nGalaxies: N,
  oldBaseline: { gc: +oldBLgc.gc.toFixed(1) },
  newBaseline: { gc: +newBLgc.gc.toFixed(1) },
  proxies: results.map(r => ({
    name: r.name, circ: r.circ, rawR: +r.rRaw.toFixed(3), rawT: +r.tRaw.toFixed(1),
    newBLr: +r.rNewBL.toFixed(3), newBLt: +r.tNewBL.toFixed(1),
    looGC: +r.gc.toFixed(1), adds: r.adds, verdict: r.verdict,
  })),
  bestProxy: best.name, finalVerdict: verdict,
};

fs.writeFileSync(path.join(__dirname, '../public/phase38-fdm.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase38-fdm.json");
