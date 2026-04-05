#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "58.0.0";
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
  const fgas = sp.MHI * 1.33 / Mbar_sun;
  const logFgas = Math.log10(fgas > 0 ? fgas : 1e-3);

  const RHI = sp.RHI || 0;
  const RHI_Rdisk = RHI > 0 && Rdisk > 0 ? RHI / Rdisk : 2.0;
  const logRHI_Rdisk = Math.log10(RHI_Rdisk > 0 ? RHI_Rdisk : 0.5);

  const MHI_Mstar = sp.MHI * 1.33 / (Mstar_sun > 0 ? Mstar_sun : 1e6);
  const logHIrichness = Math.log10(MHI_Mstar > 0 ? MHI_Mstar : 1e-3);

  const hiDef = logHIrichness < -0.3 ? 1 : 0;

  const T = sp.T || gal.T ?? 5;
  const colorProxy = T < 3 ? 1 : (T < 6 ? 0.5 : 0);

  const rpSusceptibility = fgas > 0.3 ? Math.log10(fgas * (sp.Vflat > 0 ? 1/sp.Vflat : 0.01)) : -2;

  const outerSlope = (() => {
    const outers = sorted.filter((_, i) => rArr[i] / rMax > 0.7);
    if (outers.length < 3) return 0;
    const ox = outers.map(p => p.r);
    const oy = outers.map(p => p.log_g_obs);
    const mx2 = ox.reduce((s,v)=>s+v,0)/ox.length;
    const my2 = oy.reduce((s,v)=>s+v,0)/oy.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < ox.length; i++) { sxy += (ox[i]-mx2)*(oy[i]-my2); sxx += (ox[i]-mx2)**2; }
    return sxx > 0 ? sxy / sxx : 0;
  })();

  const envGasCombined = envCode * logFgas;
  const envMassCombined = envCode * logMstar;

  const preprocessFlag = envCode >= 1 && fgas < 0.3 ? 1 : 0;
  const strippedFlag = envCode >= 1 && RHI_Rdisk < 1.5 ? 1 : 0;
  const harassedFlag = envCode >= 1 && rcWiggliness > 0.06 ? 1 : 0;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, rcWiggliness, envCode, logSigma0, logMeanRun,
    logFgas, logRHI_Rdisk, logHIrichness, hiDef,
    colorProxy, rpSusceptibility, outerSlope,
    envGasCombined, envMassCombined,
    preprocessFlag, strippedFlag, harassedFlag,
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
log("║  PHASE 58: ENVIRONMENTAL PROCESSING HISTORY                                   ║");
log("║  Version " + VERSION + ", " + N + " galaxies                                           ║");
log("╚══════════════════════════════════════════════════════════════════════════════════╝");
log("");
log("  HYPOTHESIS: envCode works because what matters is not current density,");
log("  but processing history — time since infall, stripping, harassment.");
log("");
log("  Conservative BL: " + gcCons.gc.toFixed(1) + "% | Extended BL: " + gcWork.gc.toFixed(1) + "%");
log("");

const procVars = [
  { name: 'logFgas', arr: galaxyData.map(g => g.logFgas), note: 'gas fraction (stripping proxy)', circ: 'ZERO' },
  { name: 'logRHI/Rd', arr: galaxyData.map(g => g.logRHI_Rdisk), note: 'HI extent relative to disk (truncation proxy)', circ: 'ZERO' },
  { name: 'logHIrich', arr: galaxyData.map(g => g.logHIrichness), note: 'HI richness = MHI/Mstar (depletion proxy)', circ: 'ZERO' },
  { name: 'hiDef', arr: galaxyData.map(g => g.hiDef), note: 'HI deficiency flag', circ: 'ZERO' },
  { name: 'colorPrx', arr: galaxyData.map(g => g.colorProxy), note: 'T-type color proxy (quenching proxy)', circ: 'ZERO' },
  { name: 'rpSusc', arr: galaxyData.map(g => g.rpSusceptibility), note: 'ram-pressure susceptibility (fgas/Vflat)', circ: 'LOW' },
  { name: 'outerSlp', arr: galaxyData.map(g => g.outerSlope), note: 'outer RC decline (gas removal signature)', circ: 'MED' },
  { name: 'env*gas', arr: galaxyData.map(g => g.envGasCombined), note: 'environment x gas fraction', circ: 'ZERO' },
  { name: 'env*mass', arr: galaxyData.map(g => g.envMassCombined), note: 'environment x stellar mass', circ: 'ZERO' },
  { name: 'preProcF', arr: galaxyData.map(g => g.preprocessFlag), note: 'preprocessing flag (env>=1 & fgas<0.3)', circ: 'ZERO' },
  { name: 'strippedF', arr: galaxyData.map(g => g.strippedFlag), note: 'stripped flag (env>=1 & RHI/Rd<1.5)', circ: 'ZERO' },
  { name: 'harassedF', arr: galaxyData.map(g => g.harassedFlag), note: 'harassed flag (env>=1 & wig>0.06)', circ: 'MED' },
];

log("  ┌──────────────┬──────┬────────┬────────┬────────┬────────┬────────┬────────┬──────────┬───────────┐");
log("  │ Variable     │ Circ │ |Cs r  │ Cs LOO │ Cs p   │ |Wk r  │ Wk LOO │ Wk p   │ r(env)   │ Verdict   │");
log("  ├──────────────┼──────┼────────┼────────┼────────┼────────┼────────┼────────┼──────────┼───────────┤");

const NPERMS = 10000;
const procResults = [];

for (const pv of procVars) {
  const sdV = Math.sqrt(pv.arr.map(v => (v - pv.arr.reduce((s,v)=>s+v,0)/N)**2).reduce((s,v)=>s+v,0)/(N-1));
  if (sdV < 1e-8) {
    procResults.push({ name: pv.name, circ: pv.circ, note: pv.note,
      cons: { r: 0, loo: gcCons.gc, permP: 1 }, work: { r: 0, loo: gcWork.gc, permP: 1 },
      rEnv: 0, verdict: 'FAIL (no variance)' });
    log("  │ " + pv.name.padEnd(12) + " │ " + pv.circ.padEnd(4) + " │    0.0 │  " +
      gcCons.gc.toFixed(1) + "% │  1.000 │    0.0 │  " + gcWork.gc.toFixed(1) + "% │  1.000 │    0.000 │ FAIL      │");
    continue;
  }

  const rxC = residualize(pv.arr, CONS_BL);
  const ryC = residualize(deltaA0, CONS_BL);
  const rC = corrWith(rxC, ryC);
  const gcC = looGapClosed([...CONS_BL, pv.arr], allLogA0, galaxyData, N);

  const rxW = residualize(pv.arr, WORK_BL);
  const ryW = residualize(deltaA0, WORK_BL);
  const rW = corrWith(rxW, ryW);
  const gcW = looGapClosed([...WORK_BL, pv.arr], allLogA0, galaxyData, N);

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

  const rEnv = corrWith(pv.arr, envArr);

  let verdict;
  if (permW < 0.05 && gcW.gc > gcWork.gc + 1) verdict = 'SURVIVOR';
  else if (permC < 0.05 && gcC.gc > gcCons.gc + 1) verdict = 'PARTIAL';
  else verdict = 'FAIL';

  procResults.push({
    name: pv.name, circ: pv.circ, note: pv.note,
    cons: { r: +rC.toFixed(3), loo: +gcC.gc.toFixed(1), permP: +permC.toFixed(4) },
    work: { r: +rW.toFixed(3), loo: +gcW.gc.toFixed(1), permP: +permW.toFixed(4) },
    rEnv: +rEnv.toFixed(3),
    verdict,
  });

  log("  │ " + pv.name.padEnd(12) + " │ " + pv.circ.padEnd(4) + " │ " +
    rC.toFixed(3).padStart(6) + " │ " + (gcC.gc.toFixed(1)+"%").padStart(6) + " │ " +
    permC.toFixed(3).padStart(6) + " │ " +
    rW.toFixed(3).padStart(6) + " │ " + (gcW.gc.toFixed(1)+"%").padStart(6) + " │ " +
    permW.toFixed(3).padStart(6) + " │ " + rEnv.toFixed(3).padStart(8) + " │ " + verdict.padEnd(9) + " │");
}
log("  └──────────────┴──────┴────────┴────────┴────────┴────────┴────────┴────────┴──────────┴───────────┘");

const survivors = procResults.filter(r => r.verdict === 'SURVIVOR');
const partials = procResults.filter(r => r.verdict === 'PARTIAL');

log("");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  ABSORPTION TEST: Does processing history absorb envCode?");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

const sigVars = [...survivors, ...partials];
if (sigVars.length > 0) {
  for (const sv of sigVars) {
    const svArr = procVars.find(pv => pv.name === sv.name).arr;
    const noEnvBL = [mhiArr, wigArr, sigma0Arr];
    const withProc = [...noEnvBL, svArr];

    const rxEnv = residualize(envArr, withProc);
    const ryEnv = residualize(deltaA0, withProc);
    const envAfter = corrWith(rxEnv, ryEnv);

    const rxProc = residualize(svArr, [...noEnvBL, envArr]);
    const ryProc = residualize(deltaA0, [...noEnvBL, envArr]);
    const procAfterEnv = corrWith(rxProc, ryProc);

    log("  After controlling for " + sv.name + " (replacing envCode):");
    log("    envCode partial r: " + envAfter.toFixed(3) + " (was -0.385 in baseline)");
    log("    " + sv.name + " partial r after envCode: " + procAfterEnv.toFixed(3));
    log("    Absorption of envCode: " + ((1 - Math.abs(envAfter)/0.385)*100).toFixed(0) + "%");

    const gcReplace = looGapClosed([mhiArr, wigArr, svArr, sigma0Arr], allLogA0, galaxyData, N);
    log("    LOO with " + sv.name + " instead of envCode: " + gcReplace.gc.toFixed(1) + "% (was " + gcCons.gc.toFixed(1) + "%)");
    log("");
  }
} else {
  log("  No processing variables survived. envCode cannot be absorbed.");
  log("  envCode remains a simple threshold marker (field/group/cluster).");
  log("  True processing history proxies (time since infall, tidal history)");
  log("  are not available in SPARC data.");
  log("");
}

log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  PHASE 58 — VERDICT");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");

if (survivors.length > 0) {
  log("  RESULT: Processing history variables ADD beyond baseline.");
  log("  Survivors: " + survivors.map(s => s.name).join(", "));
  log("  envCode may be partially explained by processing history.");
} else if (partials.length > 0) {
  log("  RESULT: Some processing proxies show partial signal.");
  log("  Partials: " + partials.map(s => s.name).join(", "));
  log("  But they don't fully replace envCode.");
} else {
  log("  RESULT: NO processing proxy adds independent information beyond envCode.");
  log("  INTERPRETATION:");
  log("    1. envCode is already the optimal available proxy for environment");
  log("    2. The environmental signal is simple: field vs non-field binary");
  log("    3. True processing history (infall time, stripping history) is not");
  log("       available in SPARC and would require complementary surveys");
  log("    4. The environmental effect on a0 is a THRESHOLD, not a gradient");
}
log("");

const output = {
  version: VERSION, timestamp: new Date().toISOString(), nGalaxies: N,
  conservativeBL: +gcCons.gc.toFixed(1),
  extendedBL: +gcWork.gc.toFixed(1),
  processingVars: procResults,
  survivors: survivors.map(s => s.name),
  partials: partials.map(s => s.name),
  verdict: survivors.length > 0 ? 'PROCESSING_ADDS' : partials.length > 0 ? 'PARTIAL' : 'ENVCODE_IRREDUCIBLE',
};

fs.writeFileSync(path.join(__dirname, '../public/phase58-env-processing.json'), JSON.stringify(output, null, 2));
log("  Results saved to public/phase58-env-processing.json");
