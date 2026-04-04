#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "55.0.0";
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

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0,
    logMHI, rcWiggliness, envCode, logSigma0, logMeanRun,
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
log("║  PHASE 55: INTERACTION TERMS                                                   ║");
log("║  Version " + VERSION + ", " + N + " galaxies                                           ║");
log("╚══════════════════════════════════════════════════════════════════════════════════╝");
log("");
log("  Conservative BL: " + gcCons.gc.toFixed(1) + "% | Extended BL: " + gcWork.gc.toFixed(1) + "%");
log("");

function stdize(arr) {
  const m = arr.reduce((s,v)=>s+v,0)/arr.length;
  const sd = Math.sqrt(arr.map(v=>(v-m)**2).reduce((s,v)=>s+v,0)/(arr.length-1));
  return arr.map(v => sd > 0 ? (v - m) / sd : 0);
}

const sMHI = stdize(mhiArr);
const sWig = stdize(wigArr);
const sEnv = stdize(envArr);
const sSig = stdize(sigma0Arr);
const sMR = stdize(meanRunArr);

const interactions = [
  { name: 'MHI x env', arr: sMHI.map((v,i) => v * sEnv[i]), note: 'gas mass modulated by environment' },
  { name: 'MHI x Sig0', arr: sMHI.map((v,i) => v * sSig[i]), note: 'gas mass modulated by density' },
  { name: 'wig x env', arr: sWig.map((v,i) => v * sEnv[i]), note: 'wiggliness modulated by environment' },
  { name: 'wig x mRun', arr: sWig.map((v,i) => v * sMR[i]), note: 'wiggliness modulated by run coherence' },
  { name: 'Sig0 x mRun', arr: sSig.map((v,i) => v * sMR[i]), note: 'density modulated by run coherence' },
  { name: 'MHI x mRun', arr: sMHI.map((v,i) => v * sMR[i]), note: 'gas mass modulated by run coherence' },
  { name: 'MHI x wig', arr: sMHI.map((v,i) => v * sWig[i]), note: 'gas mass modulated by wiggliness' },
  { name: 'env x Sig0', arr: sEnv.map((v,i) => v * sSig[i]), note: 'environment modulated by density' },
  { name: 'env x mRun', arr: sEnv.map((v,i) => v * sMR[i]), note: 'environment modulated by run coherence' },
  { name: 'wig x Sig0', arr: sWig.map((v,i) => v * sSig[i]), note: 'wiggliness modulated by density' },
];

log("  Testing " + interactions.length + " interaction terms:");
log("");

log("  ┌───────────────┬────────┬────────┬────────┬────────┬────────┬────────┬───────────────────────┐");
log("  │ Interaction   │ |Cs t  │ Cs LOO │ Cs p   │ |Wk t  │ Wk LOO │ Wk p   │ Verdict               │");
log("  ├───────────────┼────────┼────────┼────────┼────────┼────────┼────────┼───────────────────────┤");

const intResults = [];

for (const ix of interactions) {
  const rxC = residualize(ix.arr, CONS_BL);
  const ryC = residualize(deltaA0, CONS_BL);
  const rC = corrWith(rxC, ryC);
  const tC = Math.abs(rC) * Math.sqrt(N-6) / Math.sqrt(1-rC**2+1e-10);
  const gcC = looGapClosed([...CONS_BL, ix.arr], allLogA0, galaxyData, N);

  const rxW = residualize(ix.arr, WORK_BL);
  const ryW = residualize(deltaA0, WORK_BL);
  const rW = corrWith(rxW, ryW);
  const tW = Math.abs(rW) * Math.sqrt(N-7) / Math.sqrt(1-rW**2+1e-10);
  const gcW = looGapClosed([...WORK_BL, ix.arr], allLogA0, galaxyData, N);

  const NPERMS = 10000;
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

  const addsC = gcC.gc > gcCons.gc + 1;
  const addsW = gcW.gc > gcWork.gc + 1;

  let verdict;
  if (permW < 0.05 && addsW) verdict = 'SURVIVOR (extended)';
  else if (permC < 0.05 && addsC) verdict = 'SURVIVOR (conservative)';
  else if (permC < 0.05) verdict = 'MARGINAL';
  else verdict = 'FAIL';

  intResults.push({
    name: ix.name, note: ix.note,
    cons: { r: +rC.toFixed(3), t: +tC.toFixed(1), loo: +gcC.gc.toFixed(1), permP: +permC.toFixed(4) },
    work: { r: +rW.toFixed(3), t: +tW.toFixed(1), loo: +gcW.gc.toFixed(1), permP: +permW.toFixed(4) },
    addsC, addsW, verdict,
  });

  log("  │ " + ix.name.padEnd(13) + " │ " +
    tC.toFixed(1).padStart(6) + " │ " + (gcC.gc.toFixed(1)+"%").padStart(6) + " │ " +
    permC.toFixed(3).padStart(6) + " │ " +
    tW.toFixed(1).padStart(6) + " │ " + (gcW.gc.toFixed(1)+"%").padStart(6) + " │ " +
    permW.toFixed(3).padStart(6) + " │ " +
    verdict.padEnd(21) + " │");
}

log("  └───────────────┴────────┴────────┴────────┴────────┴────────┴────────┴───────────────────────┘");
log("");

const ixSurvivors = intResults.filter(r => r.verdict.includes('SURVIVOR'));
const ixMarginals = intResults.filter(r => r.verdict === 'MARGINAL');

if (ixSurvivors.length > 0) {
  log("  SURVIVORS:");
  for (const s of ixSurvivors) {
    log("    " + s.name + ": " + s.note);
    log("      cons LOO=" + s.cons.loo + "% (perm=" + s.cons.permP + ")");
    log("      work LOO=" + s.work.loo + "% (perm=" + s.work.permP + ")");

    log("    Deep validation:");
    const ixArr = interactions.find(i => i.name === s.name).arr;
    const bl = s.verdict.includes('extended') ? WORK_BL : CONS_BL;
    const rx = residualize(ixArr, bl);
    const ry = residualize(deltaA0, bl);

    const NBOOT = 2000;
    const bootR = [];
    for (let b = 0; b < NBOOT; b++) {
      const idx = Array.from({length: N}, () => Math.floor(Math.random() * N));
      const bX = idx.map(i => ixArr[i]);
      const bY = idx.map(i => deltaA0[i]);
      const bBL = bl.map(arr => idx.map(i => arr[i]));
      const brx = residualize(bX, bBL);
      const bry = residualize(bY, bBL);
      bootR.push(corrWith(brx, bry));
    }
    bootR.sort((a,b) => a - b);
    const lo = bootR[Math.floor(NBOOT * 0.025)];
    const hi = bootR[Math.floor(NBOOT * 0.975)];
    log("      Bootstrap 95% CI: [" + lo.toFixed(3) + ", " + hi.toFixed(3) + "]" +
      (lo > 0 || hi < 0 ? " excludes zero" : " crosses zero"));

    const jackR = [];
    for (let i = 0; i < N; i++) {
      const idx = [...Array(N).keys()].filter(j => j !== i);
      const subX = idx.map(j => ixArr[j]);
      const subY = idx.map(j => deltaA0[j]);
      const subBL = bl.map(arr => idx.map(j => arr[j]));
      const brx = residualize(subX, subBL);
      const bry = residualize(subY, subBL);
      jackR.push(corrWith(brx, bry));
    }
    const jackMean = jackR.reduce((s,v)=>s+v,0)/N;
    const signFlips = jackR.filter(r => Math.sign(r) !== Math.sign(corrWith(rx, ry))).length;
    log("      Jackknife: mean=" + jackMean.toFixed(3) + " flips=" + signFlips + "/" + N);
    log("");
  }
}

if (ixMarginals.length > 0) {
  log("  MARGINALS: " + ixMarginals.map(m => m.name).join(', '));
}

const ixFails = intResults.filter(r => r.verdict === 'FAIL');
if (ixFails.length > 0) {
  log("  FAILED: " + ixFails.map(f => f.name).join(', '));
}

log("");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("  PHASE 55 — CONCLUSION");
log("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━");
log("");
log("  Conservative BL: " + gcCons.gc.toFixed(1) + "% (frozen)");
log("  Extended BL:     " + gcWork.gc.toFixed(1) + "% (frozen)");
log("  Interaction survivors: " + (ixSurvivors.length > 0 ? ixSurvivors.map(s => s.name).join(', ') : 'NONE'));
log("  Interaction marginals: " + (ixMarginals.length > 0 ? ixMarginals.map(s => s.name).join(', ') : 'NONE'));
log("");

const output = {
  version: VERSION, timestamp: new Date().toISOString(), nGalaxies: N,
  conservativeBL: +gcCons.gc.toFixed(1),
  extendedBL: +gcWork.gc.toFixed(1),
  interactions: intResults,
  survivors: ixSurvivors.map(s => s.name),
  marginals: ixMarginals.map(s => s.name),
};

fs.writeFileSync(path.join(__dirname, '../public/phase55-interactions.json'), JSON.stringify(output, null, 2));
log("  Results saved to public/phase55-interactions.json");
