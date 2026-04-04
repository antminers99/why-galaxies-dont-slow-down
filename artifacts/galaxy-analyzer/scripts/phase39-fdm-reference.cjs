#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "39.0.0";
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

  function fDMAtRadius(rRef) {
    const nearby = sorted.filter(p => Math.abs(p.r - rRef) < rRef * 0.4 || Math.abs(p.r - rRef) < 1.0);
    if (nearby.length === 0) return null;
    const closest = nearby.reduce((a, b) => Math.abs(a.r - rRef) < Math.abs(b.r - rRef) ? a : b);
    const gObs = Math.pow(10, closest.log_g_obs);
    const gBar = Math.pow(10, closest.log_g_bar);
    const ratio = gBar / gObs;
    const fDMlocal = 1 - ratio;
    return Math.max(-0.5, Math.min(0.999, fDMlocal));
  }

  const fDM_2Rd = fDMAtRadius(2 * Rdisk);
  const fDM_3Rd = fDMAtRadius(3 * Rdisk);
  const fDM_Rmax50 = fDMAtRadius(rMax * 0.5);
  const fDM_5kpc = fDMAtRadius(5.0);
  const fDM_10kpc = fDMAtRadius(10.0);

  const dmDominanceR = (() => {
    for (const p of sorted) {
      const gObs = Math.pow(10, p.log_g_obs);
      const gBar = Math.pow(10, p.log_g_bar);
      if (gBar / gObs < 0.5) return p.r;
    }
    return rMax;
  })();
  const logDmDomR = Math.log10(dmDominanceR > 0.1 ? dmDominanceR : 0.1);
  const dmDomRnorm = dmDominanceR / Rdisk;
  const logDmDomRnorm = Math.log10(dmDomRnorm > 0.01 ? dmDomRnorm : 0.01);

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, logVflat, rcWiggliness, envCode, logSigma0,
    fDM_2Rd, fDM_3Rd, fDM_Rmax50, fDM_5kpc, fDM_10kpc,
    logDmDomR, logDmDomRnorm,
    Rdisk, rMax,
  });
}

const valid2Rd = galaxyData.filter(g => g.fDM_2Rd !== null);
const valid3Rd = galaxyData.filter(g => g.fDM_3Rd !== null);
const valid50 = galaxyData.filter(g => g.fDM_Rmax50 !== null);
const valid5k = galaxyData.filter(g => g.fDM_5kpc !== null);
const valid10k = galaxyData.filter(g => g.fDM_10kpc !== null);

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const sigma0Arr = galaxyData.map(g => g.logSigma0);
const vflatArr = galaxyData.map(g => g.logVflat);
const dmDomRArr = galaxyData.map(g => g.logDmDomR);
const dmDomRnormArr = galaxyData.map(g => g.logDmDomRnorm);

const NEW_BASELINE = [mhiArr, wigArr, envArr, sigma0Arr];

log("");
log("=".repeat(80));
log("  PHASE 39: f_DM at Reference Radii");
log("  FAST PROTOCOL — Baseline: MHI + wig + env + Sigma0_bar");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  PROXIES:");
log("  - f_DM(2Rd): DM fraction at 2 disk scale lengths (N=" + valid2Rd.length + ")");
log("  - f_DM(3Rd): DM fraction at 3 disk scale lengths (N=" + valid3Rd.length + ")");
log("  - f_DM(Rmax/2): DM fraction at half last radius (N=" + valid50.length + ")");
log("  - f_DM(5kpc): DM fraction at 5 kpc (N=" + valid5k.length + ")");
log("  - f_DM(10kpc): DM fraction at 10 kpc (N=" + valid10k.length + ")");
log("  - R_DM_dom: radius where DM first dominates (g_bar/g_obs < 0.5)");
log("  - R_DM_dom/Rdisk: normalized DM dominance radius");
log("");

log("  CIRCULARITY: HIGH for all (from g_obs/g_bar ratio in RC)");
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'envCode': envArr,
  'logSigma0': sigma0Arr, 'logVflat': vflatArr,
  'logDmDomR': dmDomRArr, 'logDmDomRnorm': dmDomRnormArr,
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

const baselineGC = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', 'logSigma0']);
log("  BASELINE: " + baselineGC.gc.toFixed(1) + "% gap closure");
log("");

log("=".repeat(80));
log("  STEP 1: FAST SCREEN — f_DM at reference radii");
log("=".repeat(80));
log("");

const refRadiusTests = [
  { name: 'f_DM(2Rd)', data: valid2Rd, key: 'fDM_2Rd' },
  { name: 'f_DM(3Rd)', data: valid3Rd, key: 'fDM_3Rd' },
  { name: 'f_DM(Rmax/2)', data: valid50, key: 'fDM_Rmax50' },
  { name: 'f_DM(5kpc)', data: valid5k, key: 'fDM_5kpc' },
  { name: 'f_DM(10kpc)', data: valid10k, key: 'fDM_10kpc' },
];

log("  ┌──────────────────────────────────────────────────────────────────────────────┐");
log("  │  Proxy           N    Raw r    Raw t   |newBL r  |newBL t   Verdict        │");
log("  ├──────────────────────────────────────────────────────────────────────────────┤");

const refResults = [];
for (const test of refRadiusTests) {
  const subset = test.data;
  const n = subset.length;
  if (n < 15) {
    log("  │  " + test.name.padEnd(15) + n.toString().padStart(4) + "   — too few galaxies —".padEnd(45) + "SKIP    │");
    refResults.push({ name: test.name, n, verdict: 'SKIP' });
    continue;
  }

  const subIdx = subset.map(g => galaxyData.indexOf(g));
  const subFdm = subIdx.map(i => galaxyData[i][test.key]);
  const subDa = subIdx.map(i => deltaA0[i]);
  const subBL = NEW_BASELINE.map(arr => subIdx.map(i => arr[i]));

  const rRaw = corrWith(subFdm, subDa);
  const tRaw = Math.abs(rRaw) * Math.sqrt(n-2) / Math.sqrt(1-rRaw**2+1e-10);

  const rx = residualize(subFdm, subBL);
  const ry = residualize(subDa, subBL);
  const rBL = corrWith(rx, ry);
  const df = n - 6;
  const tBL = df > 0 ? Math.abs(rBL) * Math.sqrt(df) / Math.sqrt(1-rBL**2+1e-10) : 0;

  const failFast = tRaw < 1.0 && tBL < 1.0;
  const promising = tBL >= 1.65;
  const verdict = failFast ? "FAIL-FAST" : promising ? "PROMISING" : "WEAK";

  refResults.push({ name: test.name, n, rRaw, tRaw, rBL, tBL, verdict });

  log("  │  " + test.name.padEnd(15) + n.toString().padStart(4) +
    rRaw.toFixed(3).padStart(8) + tRaw.toFixed(1).padStart(7) +
    rBL.toFixed(3).padStart(9) + tBL.toFixed(1).padStart(9) + "   " + verdict.padEnd(12) + "│");
}
log("  └──────────────────────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  STEP 2: DM DOMINANCE RADIUS (full sample, N=" + N + ")");
log("=".repeat(80));
log("");

const domCandidates = [
  { name: 'R_DM_dom', arr: dmDomRArr, feat: 'logDmDomR', circ: 'MOD' },
  { name: 'R_DM_dom/Rd', arr: dmDomRnormArr, feat: 'logDmDomRnorm', circ: 'MOD' },
];

log("  ┌──────────────────────────────────────────────────────────────────────────────────────┐");
log("  │  Proxy         Circ   Raw r    Raw t   |newBL r  |newBL t  LOO+newBL Adds? Verdict │");
log("  ├──────────────────────────────────────────────────────────────────────────────────────┤");

const domResults = [];
for (const cand of domCandidates) {
  const rRaw = corrWith(cand.arr, deltaA0);
  const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);

  const rx = residualize(cand.arr, NEW_BASELINE);
  const ry = residualize(deltaA0, NEW_BASELINE);
  const rBL = corrWith(rx, ry);
  const tBL = Math.abs(rBL) * Math.sqrt(N-6) / Math.sqrt(1-rBL**2+1e-10);

  const gc = looGapClosed(['logMHI', 'rcWiggliness', 'envCode', 'logSigma0', cand.feat]);
  const adds = gc.gc > baselineGC.gc + 1;

  const failFast = tRaw < 1.0 && tBL < 1.0;
  const promising = tBL >= 1.65 || adds;
  const verdict = failFast ? "FAIL-FAST" : promising ? "PROMISING" : "WEAK";

  domResults.push({ ...cand, rRaw, tRaw, rBL, tBL, gc: gc.gc, adds, verdict });

  log("  │  " + cand.name.padEnd(13) + cand.circ.padEnd(5) +
    rRaw.toFixed(3).padStart(7) + tRaw.toFixed(1).padStart(7) +
    rBL.toFixed(3).padStart(9) + tBL.toFixed(1).padStart(9) +
    (gc.gc.toFixed(1) + "%").padStart(10) +
    (adds ? " YES " : " NO  ") + verdict.padStart(10) + " │");
}
log("  └──────────────────────────────────────────────────────────────────────────────────────┘");
log("");

const allPromising = [...refResults.filter(r => r.verdict === "PROMISING"), ...domResults.filter(r => r.verdict === "PROMISING")];

if (allPromising.length === 0) {
  log("  ══════════════════════════════════════════════════════════");
  log("  NO f_DM-at-reference-radius proxy passes fast screen.");
  log("  ALL FAIL or WEAK. Skipping deep tests.");
  log("  ══════════════════════════════════════════════════════════");
} else {
  log("  PROMISING: " + allPromising.map(p => p.name).join(", "));
  log("  Proceeding to deep tests...");
  log("");

  for (const cand of domResults.filter(r => r.verdict === "PROMISING")) {
    log("=".repeat(80));
    log("  DEEP TEST: " + cand.name);
    log("=".repeat(80));
    log("");

    log("  Confounder stripping:");
    for (const [label, ctrls] of [
      ['raw', []], ['|MHI', [mhiArr]],
      ['|newBL', NEW_BASELINE], ['|newBL+Vflat', [...NEW_BASELINE, vflatArr]],
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
    for (const [label, ctrls] of [['|newBL', NEW_BASELINE]]) {
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
  }
}

log("=".repeat(80));
log("  PHASE 39 — FAST PROTOCOL SUMMARY");
log("=".repeat(80));
log("");

const allResults = [...refResults, ...domResults];
const bestResult = domResults.length > 0
  ? domResults.reduce((a,b) => Math.abs(b.tBL || 0) > Math.abs(a.tBL || 0) ? b : a)
  : { name: 'none', tRaw: 0, tBL: 0, adds: false, gc: baselineGC.gc, verdict: 'FAIL' };

const anyPassed = allResults.some(r => r.verdict === "PROMISING");

log("  ┌──────────────────────────────────────────────────────────────────────┐");
log("  │  Test                        Result              Pass/Fail        │");
log("  ├──────────────────────────────────────────────────────────────────────┤");
log("  │  Any ref-radius f_DM passes  " + (anyPassed ? "YES" : "NO").padEnd(24) + (anyPassed ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  Best DM-dominance proxy     " + bestResult.name.padEnd(24) + " ".repeat(12) + "│");
log("  │  After NEW baseline          t=" + (bestResult.tBL||0).toFixed(1) + "                   " + ((bestResult.tBL||0) >= 1.65 ? "PASS" : "FAIL").padEnd(12) + "│");
log("  │  LOO adds                    " + (bestResult.adds ? "YES" : "NO").padEnd(24) + (bestResult.adds ? "PASS" : "FAIL").padEnd(12) + "│");
log("  ├──────────────────────────────────────────────────────────────────────┤");

const finalVerdict = anyPassed ? "PARTIAL" : "FAIL";
log("  │  FINAL VERDICT: " + finalVerdict.padEnd(50) + "│");
log("  └──────────────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  DARK HALO DOOR — FINAL SCORECARD");
log("=".repeat(80));
log("");
log("  ┌───────────────────────────────────────────────────────────────────────┐");
log("  │  Variable              Verdict          Notes                       │");
log("  ├───────────────────────────────────────────────────────────────────────┤");
log("  │  M_halo / M_dyn        PARTIAL          = V_flat, circularity 2/5  │");
log("  │  concentration c       FAIL             0/4 fast screen            │");
log("  │  V_flat                PASS+CIRCULAR    t=3.8 but mechanically     │");
log("  │  R_s                   PARTIAL          Rs_Vrise, perm p=0.08      │");
log("  │  rho_0 (halo proper)   FAIL             all HIGH circularity       │");
log("  │  Sigma0_bar            CONFIRMED        baryon structure, p=0.007  │");
log("  │  f_DM (global)         FAIL             absorbed by Sigma0_bar     │");
log("  │  f_DM (ref radius)     " + finalVerdict.padEnd(17) + "                             │");
log("  └───────────────────────────────────────────────────────────────────────┘");
log("");

log("  BASELINE: MHI + rcWiggliness + envCode + Sigma0_bar = " + baselineGC.gc.toFixed(1) + "%");
log("");
log("=".repeat(80));

const output = {
  version: VERSION, timestamp: new Date().toISOString(), nGalaxies: N,
  baseline: { gc: +baselineGC.gc.toFixed(1) },
  refRadiusResults: refResults.map(r => ({ name: r.name, n: r.n, rawR: r.rRaw ? +r.rRaw.toFixed(3) : null, rawT: r.tRaw ? +r.tRaw.toFixed(1) : null, blR: r.rBL ? +r.rBL.toFixed(3) : null, blT: r.tBL ? +r.tBL.toFixed(1) : null, verdict: r.verdict })),
  dmDominanceResults: domResults.map(r => ({ name: r.name, rawR: +r.rRaw.toFixed(3), rawT: +r.tRaw.toFixed(1), blR: +r.rBL.toFixed(3), blT: +r.tBL.toFixed(1), looGC: +r.gc.toFixed(1), adds: r.adds, verdict: r.verdict })),
  finalVerdict,
};

fs.writeFileSync(path.join(__dirname, '../public/phase39-fdm-reference.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase39-fdm-reference.json");
