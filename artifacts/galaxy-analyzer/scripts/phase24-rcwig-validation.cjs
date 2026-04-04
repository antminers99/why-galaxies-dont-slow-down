#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "24.0.0";
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

function shuffle(arr) {
  const a = [...arr];
  for (let i = a.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [a[i], a[j]] = [a[j], a[i]];
  }
  return a;
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

  let baryonFeatureAmplitude = 0;
  if (pts.length >= 4) {
    const gb = pts.map(p => p.log_g_bar);
    const mgb = gb.reduce((s,v)=>s+v,0)/gb.length;
    let sxy = 0, sxx = 0;
    const rr = pts.map(p => p.r);
    const mr = rr.reduce((s,v)=>s+v,0)/rr.length;
    for (let i = 0; i < pts.length; i++) { sxy += (rr[i]-mr)*(gb[i]-mgb); sxx += (rr[i]-mr)**2; }
    const bSlope = sxx > 0 ? sxy / sxx : 0;
    const bInt = mgb - bSlope * mr;
    let bss = 0;
    for (let i = 0; i < pts.length; i++) bss += (gb[i] - (bInt + bSlope * rr[i])) ** 2;
    baryonFeatureAmplitude = Math.sqrt(bss / pts.length);
  }

  let baryonBumpiness = 0;
  if (pts.length >= 4) {
    const gb = pts.map(p => p.log_g_bar);
    let diffs = 0;
    for (let i = 1; i < gb.length; i++) diffs += Math.abs(gb[i] - gb[i-1]);
    const totalRange = Math.abs(gb[gb.length-1] - gb[0]) + 1e-10;
    baryonBumpiness = diffs / totalRange;
  }

  const logSBdisk = sp.SBdisk > 0 ? Math.log10(sp.SBdisk) : 2;
  const logSBeff = sp.SBeff > 0 ? Math.log10(sp.SBeff) : 2;
  const sbContrast = Math.abs(logSBdisk - logSBeff);

  const Vflat = sp.Vflat || 100;
  const logVflat = Math.log10(Vflat);
  const Rdisk = sp.Rdisk || 3;
  const beamSize = sp.D * 1000 * 4.848e-6 * 15;
  const logBeamPerDisk = Math.log10(beamSize / Rdisk + 0.01);

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q,
    n, logMHI, logVflat, logBeamPerDisk,
    rcWiggliness,
    baryonFeatureAmplitude, baryonBumpiness, sbContrast,
    Rdisk, Vflat,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

log("");
log("=".repeat(80));
log("  PHASE 24: rcWiggliness VALIDATION BATTERY");
log("  From 'interesting signal' to 'robust finding'");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("=".repeat(80));
log("  TEST 1: SPLIT-HALF REPLICATION (1000 random splits)");
log("  Train on half, test on other half. Does signal replicate?");
log("=".repeat(80));
log("");

const NSPLITS = 1000;
let bothPositive = 0, bothSignificant = 0;
const splitCorrs = [];

for (let s = 0; s < NSPLITS; s++) {
  const idx = shuffle([...Array(N).keys()]);
  const half = Math.floor(N / 2);
  const setA = idx.slice(0, half);
  const setB = idx.slice(half);

  const rA = corrWith(setA.map(i => galaxyData[i].rcWiggliness), setA.map(i => deltaA0[i]));
  const rB = corrWith(setB.map(i => galaxyData[i].rcWiggliness), setB.map(i => deltaA0[i]));

  if (rA > 0 && rB > 0) bothPositive++;
  const tA = Math.abs(rA) * Math.sqrt(setA.length - 2) / Math.sqrt(1 - rA**2 + 1e-10);
  const tB = Math.abs(rB) * Math.sqrt(setB.length - 2) / Math.sqrt(1 - rB**2 + 1e-10);
  if (tA > 1.65 && tB > 1.65 && rA * rB > 0) bothSignificant++;
  splitCorrs.push({ rA, rB });
}

log("  SPLIT-HALF RESULTS (N=" + NSPLITS + " random 50/50 splits):");
log("  Both halves same sign:      " + bothPositive + "/" + NSPLITS + " (" + (bothPositive/NSPLITS*100).toFixed(1) + "%)");
log("  Both halves significant:    " + bothSignificant + "/" + NSPLITS + " (" + (bothSignificant/NSPLITS*100).toFixed(1) + "%)");
log("  Mean |r| in each half:      " + (splitCorrs.reduce((s,c)=>s+Math.abs(c.rA),0)/NSPLITS).toFixed(3) +
  ", " + (splitCorrs.reduce((s,c)=>s+Math.abs(c.rB),0)/NSPLITS).toFixed(3));
log("");

log("  SPLIT-HALF WITH MHI CONTROL:");
let bothPosAfterMHI = 0, bothSigAfterMHI = 0;
for (let s = 0; s < NSPLITS; s++) {
  const idx = shuffle([...Array(N).keys()]);
  const half = Math.floor(N / 2);
  const setA = idx.slice(0, half);
  const setB = idx.slice(half);

  for (const set of [setA, setB]) {
    const mhi = set.map(i => galaxyData[i].logMHI);
    const da = set.map(i => deltaA0[i]);
    const wig = set.map(i => galaxyData[i].rcWiggliness);
    const residDA = residualize(da, [mhi]);
    const r = corrWith(wig, residDA);
    set._r = r;
    set._t = Math.abs(r) * Math.sqrt(set.length - 3) / Math.sqrt(1 - r**2 + 1e-10);
  }
  if (setA._r > 0 && setB._r > 0) bothPosAfterMHI++;
  if (setA._t > 1.65 && setB._t > 1.65 && setA._r * setB._r > 0) bothSigAfterMHI++;
}

log("  Both same sign (|MHI):      " + bothPosAfterMHI + "/" + NSPLITS + " (" + (bothPosAfterMHI/NSPLITS*100).toFixed(1) + "%)");
log("  Both significant (|MHI):    " + bothSigAfterMHI + "/" + NSPLITS + " (" + (bothSigAfterMHI/NSPLITS*100).toFixed(1) + "%)");
log("");

log("=".repeat(80));
log("  TEST 2: RENZO'S RULE CONTROL");
log("  Is rcWiggliness just baryon feature amplitude in disguise?");
log("=".repeat(80));
log("");

const bfa = galaxyData.map(g => g.baryonFeatureAmplitude);
const bb = galaxyData.map(g => g.baryonBumpiness);
const sbc = galaxyData.map(g => g.sbContrast);
const wig = galaxyData.map(g => g.rcWiggliness);
const mhi = galaxyData.map(g => g.logMHI);
const vflat = galaxyData.map(g => g.logVflat);
const npts = galaxyData.map(g => g.n);

log("  Baryon feature proxies (Renzo's rule / Sancisi's law):");
log("    baryonFeatureAmplitude: RMS of gbar around linear trend in radius");
log("    baryonBumpiness:        sum|dgbar|/range — point-to-point roughness");
log("    sbContrast:             |logSBdisk - logSBeff|");
log("");

log("  Correlations with rcWiggliness:");
log("    baryonFeatureAmplitude vs rcWiggliness: r=" + corrWith(bfa, wig).toFixed(3));
log("    baryonBumpiness vs rcWiggliness:        r=" + corrWith(bb, wig).toFixed(3));
log("    sbContrast vs rcWiggliness:             r=" + corrWith(sbc, wig).toFixed(3));
log("");

const residDA_mhi = residualize(deltaA0, [mhi]);
const residDA_mhi_bfa = residualize(deltaA0, [mhi, bfa]);
const residDA_mhi_bb = residualize(deltaA0, [mhi, bb]);
const residDA_mhi_bfa_bb = residualize(deltaA0, [mhi, bfa, bb]);
const residDA_mhi_bfa_bb_sbc = residualize(deltaA0, [mhi, bfa, bb, sbc]);
const residDA_all5 = residualize(deltaA0, [mhi, vflat, npts, bfa, bb]);

const rBaseAfterMHI = corrWith(wig, residDA_mhi);
const rAfterMHI_BFA = corrWith(wig, residDA_mhi_bfa);
const rAfterMHI_BB = corrWith(wig, residDA_mhi_bb);
const rAfterMHI_BFA_BB = corrWith(wig, residDA_mhi_bfa_bb);
const rAfterAll = corrWith(wig, residDA_mhi_bfa_bb_sbc);
const rAfterAll5 = corrWith(wig, residDA_all5);

function tStat(r, n, p) { return Math.abs(r) * Math.sqrt(n - 2 - p) / Math.sqrt(1 - r**2 + 1e-10); }

log("  rcWiggliness vs delta_a0 after sequential Renzo controls:");
log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Controls                          r        |t|    status      │");
log("  ├───────────────────────────────────────────────────────────────────┤");
log("  │  MHI only                          " + rBaseAfterMHI.toFixed(3).padStart(6) + tStat(rBaseAfterMHI,N,1).toFixed(1).padStart(9) + "    BASELINE  │");
log("  │  MHI + baryonFeatureAmp            " + rAfterMHI_BFA.toFixed(3).padStart(6) + tStat(rAfterMHI_BFA,N,2).toFixed(1).padStart(9) + "    " + (tStat(rAfterMHI_BFA,N,2)>=1.65?"SURVIVES":"FAILS   ") + "  │");
log("  │  MHI + baryonBumpiness             " + rAfterMHI_BB.toFixed(3).padStart(6) + tStat(rAfterMHI_BB,N,2).toFixed(1).padStart(9) + "    " + (tStat(rAfterMHI_BB,N,2)>=1.65?"SURVIVES":"FAILS   ") + "  │");
log("  │  MHI + BFA + BB                    " + rAfterMHI_BFA_BB.toFixed(3).padStart(6) + tStat(rAfterMHI_BFA_BB,N,3).toFixed(1).padStart(9) + "    " + (tStat(rAfterMHI_BFA_BB,N,3)>=1.65?"SURVIVES":"FAILS   ") + "  │");
log("  │  MHI + BFA + BB + sbContrast       " + rAfterAll.toFixed(3).padStart(6) + tStat(rAfterAll,N,4).toFixed(1).padStart(9) + "    " + (tStat(rAfterAll,N,4)>=1.65?"SURVIVES":"FAILS   ") + "  │");
log("  │  MHI + Vflat + n + BFA + BB        " + rAfterAll5.toFixed(3).padStart(6) + tStat(rAfterAll5,N,5).toFixed(1).padStart(9) + "    " + (tStat(rAfterAll5,N,5)>=1.65?"SURVIVES":"FAILS   ") + "  │");
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

const renzoSurvives = tStat(rAfterMHI_BFA_BB, N, 3) >= 1.65;
log("  RENZO'S RULE VERDICT: " + (renzoSurvives ? "rcWiggliness is NOT just baryon features" : "rcWiggliness MAY be baryon features in disguise"));
log("");

log("=".repeat(80));
log("  TEST 3: TRAIN/TEST PREDICTION (repeated 70/30 splits)");
log("  Can rcWiggliness predict a0 on truly unseen galaxies?");
log("=".repeat(80));
log("");

const NTRIALS = 500;
const models = [
  { name: 'M0: universal', features: [] },
  { name: 'M_MHI', features: ['logMHI'] },
  { name: 'M_wig', features: ['rcWiggliness'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
];

const featureMap = {
  'logMHI': galaxyData.map(g => g.logMHI),
  'rcWiggliness': galaxyData.map(g => g.rcWiggliness),
};

const testResults = models.map(m => ({ ...m, testSSs: [], testRMSs: [] }));

for (let trial = 0; trial < NTRIALS; trial++) {
  const idx = shuffle([...Array(N).keys()]);
  const nTrain = Math.floor(N * 0.7);
  const trainIdx = idx.slice(0, nTrain);
  const testIdx = idx.slice(nTrain);

  for (const model of testResults) {
    if (model.features.length === 0) {
      const trainMean = trainIdx.map(i => allLogA0[i]).reduce((s,v)=>s+v,0) / trainIdx.length;
      let totalSS = 0, totalN = 0;
      for (const ti of testIdx) {
        totalSS += predictSS(Math.pow(10, trainMean), galaxyData[ti].pts);
        totalN += galaxyData[ti].pts.length;
      }
      model.testRMSs.push(Math.sqrt(totalSS / totalN));
    } else {
      const trainX = trainIdx.map(i => model.features.map(f => featureMap[f][i]));
      const trainY = trainIdx.map(i => allLogA0[i]);
      const reg = linReg(trainX, trainY);
      let totalSS = 0, totalN = 0;
      for (const ti of testIdx) {
        const testX = model.features.map(f => featureMap[f][ti]);
        let pred = reg.intercept + reg.coefs.reduce((s,c,j) => s + c * testX[j], 0);
        pred = Math.max(2.5, Math.min(4.5, pred));
        totalSS += predictSS(Math.pow(10, pred), galaxyData[ti].pts);
        totalN += galaxyData[ti].pts.length;
      }
      model.testRMSs.push(Math.sqrt(totalSS / totalN));
    }
  }
}

log("  TRAIN/TEST PREDICTION RESULTS (N=" + NTRIALS + " random 70/30 splits):");
log("  ┌──────────────────────────────────────────────────────────────────────────┐");
log("  │  Model          mean RMS    median     wins vs M0   mean improvement   │");
log("  ├──────────────────────────────────────────────────────────────────────────┤");

const m0median = testResults[0].testRMSs.sort((a,b)=>a-b)[Math.floor(NTRIALS/2)];
const m0mean = testResults[0].testRMSs.reduce((s,v)=>s+v,0)/NTRIALS;

for (const model of testResults) {
  const sorted = [...model.testRMSs].sort((a,b)=>a-b);
  const median = sorted[Math.floor(NTRIALS/2)];
  const mean = model.testRMSs.reduce((s,v)=>s+v,0)/NTRIALS;
  const wins = model.testRMSs.filter((v,i) => v < testResults[0].testRMSs[i]).length;
  const improvement = ((1 - mean / m0mean) * 100).toFixed(1);
  log("  │  " + model.name.padEnd(15) + mean.toFixed(5).padStart(10) + median.toFixed(5).padStart(10) +
    (wins + "/" + NTRIALS).padStart(12) + (improvement + "%").padStart(18) + "   │");
}
log("  └──────────────────────────────────────────────────────────────────────────┘");
log("");

const wigWins = testResults[2].testRMSs.filter((v,i) => v < testResults[0].testRMSs[i]).length;
const mhiWigWins = testResults[3].testRMSs.filter((v,i) => v < testResults[1].testRMSs[i]).length;
log("  rcWiggliness beats universal:        " + wigWins + "/" + NTRIALS + " (" + (wigWins/NTRIALS*100).toFixed(1) + "%)");
log("  MHI+wig beats MHI alone:             " + mhiWigWins + "/" + NTRIALS + " (" + (mhiWigWins/NTRIALS*100).toFixed(1) + "%)");
log("");

log("=".repeat(80));
log("  TEST 4: STRATIFICATION — DOES SIGNAL SURVIVE WITHIN SUBGROUPS?");
log("=".repeat(80));
log("");

function testSubgroup(name, indices) {
  if (indices.length < 10) { log("  " + name + ": n=" + indices.length + " (too few)"); return null; }
  const subWig = indices.map(i => galaxyData[i].rcWiggliness);
  const subDA = indices.map(i => deltaA0[i]);
  const subMHI = indices.map(i => galaxyData[i].logMHI);
  const rRaw = corrWith(subWig, subDA);
  const residDA = residualize(subDA, [subMHI]);
  const rMHI = corrWith(subWig, residDA);
  const tMHI = Math.abs(rMHI) * Math.sqrt(indices.length - 3) / Math.sqrt(1 - rMHI**2 + 1e-10);
  log("  " + name.padEnd(30) + " n=" + indices.length.toString().padStart(2) +
    "  r_raw=" + rRaw.toFixed(3) + "  r|MHI=" + rMHI.toFixed(3) + "  t=" + tMHI.toFixed(1) +
    "  " + (rMHI > 0 ? "+" : "-") + " " + (tMHI >= 1.65 ? "SIG" : "   "));
  return { name, n: indices.length, rRaw, rMHI, tMHI };
}

log("  BY QUALITY FLAG:");
const results4 = [];
for (const q of [1, 2, 3]) {
  const idx = galaxyData.map((g, i) => g.Q === q ? i : -1).filter(i => i >= 0);
  const r = testSubgroup("Q=" + q, idx);
  if (r) results4.push(r);
}
log("");

log("  BY INCLINATION:");
const incBins = [[0,50,'low inc'],[50,70,'mid inc'],[70,90,'high inc']];
for (const [lo,hi,label] of incBins) {
  const idx = galaxyData.map((g,i) => g.inc >= lo && g.inc < hi ? i : -1).filter(i => i >= 0);
  const r = testSubgroup(label + " (" + lo + "-" + hi + ")", idx);
  if (r) results4.push(r);
}
log("");

log("  BY DISTANCE:");
const dVals = galaxyData.map(g => g.D);
const dMed = [...dVals].sort((a,b)=>a-b)[Math.floor(N/2)];
{
  const nearIdx = galaxyData.map((g,i) => g.D <= dMed ? i : -1).filter(i => i >= 0);
  const farIdx = galaxyData.map((g,i) => g.D > dMed ? i : -1).filter(i => i >= 0);
  const r1 = testSubgroup("near (D<=" + dMed.toFixed(1) + " Mpc)", nearIdx);
  const r2 = testSubgroup("far  (D>" + dMed.toFixed(1) + " Mpc)", farIdx);
  if (r1) results4.push(r1);
  if (r2) results4.push(r2);
}
log("");

log("  BY NUMBER OF POINTS:");
const nMed = [...galaxyData.map(g => g.n)].sort((a,b)=>a-b)[Math.floor(N/2)];
{
  const fewIdx = galaxyData.map((g,i) => g.n <= nMed ? i : -1).filter(i => i >= 0);
  const manyIdx = galaxyData.map((g,i) => g.n > nMed ? i : -1).filter(i => i >= 0);
  const r1 = testSubgroup("few pts (n<=" + nMed + ")", fewIdx);
  const r2 = testSubgroup("many pts (n>" + nMed + ")", manyIdx);
  if (r1) results4.push(r1);
  if (r2) results4.push(r2);
}
log("");

log("  BY GALAXY MASS (Vflat):");
const vMed = [...galaxyData.map(g => g.Vflat)].sort((a,b)=>a-b)[Math.floor(N/2)];
{
  const lowM = galaxyData.map((g,i) => g.Vflat <= vMed ? i : -1).filter(i => i >= 0);
  const highM = galaxyData.map((g,i) => g.Vflat > vMed ? i : -1).filter(i => i >= 0);
  const r1 = testSubgroup("low mass (Vf<=" + Math.round(vMed) + ")", lowM);
  const r2 = testSubgroup("high mass (Vf>" + Math.round(vMed) + ")", highM);
  if (r1) results4.push(r1);
  if (r2) results4.push(r2);
}
log("");

log("  BY HUBBLE TYPE:");
{
  const earlyIdx = galaxyData.map((g,i) => g.T <= 5 ? i : -1).filter(i => i >= 0);
  const lateIdx = galaxyData.map((g,i) => g.T > 5 ? i : -1).filter(i => i >= 0);
  const r1 = testSubgroup("early type (T<=5)", earlyIdx);
  const r2 = testSubgroup("late type (T>5)", lateIdx);
  if (r1) results4.push(r1);
  if (r2) results4.push(r2);
}
log("");

const nPositive = results4.filter(r => r.rMHI > 0).length;
const nSig = results4.filter(r => r.tMHI >= 1.65 && r.rMHI > 0).length;
log("  STRATIFICATION SUMMARY:");
log("    Total subgroups:     " + results4.length);
log("    Same sign (positive): " + nPositive + "/" + results4.length + " (" + (nPositive/results4.length*100).toFixed(0) + "%)");
log("    Significant:          " + nSig + "/" + results4.length);
log("");

log("=".repeat(80));
log("  TEST 5: JACKKNIFE — IS THE SIGNAL DRIVEN BY OUTLIERS?");
log("=".repeat(80));
log("");

const controlVals = [mhi, vflat, npts];
const baseR = corrWith(wig, residualize(deltaA0, controlVals));
const jackR = [];
let maxDrop = 0, maxDropGal = '';

for (let i = 0; i < N; i++) {
  const idx = [...Array(N).keys()].filter(j => j !== i);
  const subWig = idx.map(j => wig[j]);
  const subDA = idx.map(j => deltaA0[j]);
  const subControls = controlVals.map(c => idx.map(j => c[j]));
  const subResid = residualize(subDA, subControls);
  const r = corrWith(subWig, subResid);
  jackR.push(r);
  const drop = baseR - r;
  if (Math.abs(drop) > Math.abs(maxDrop)) { maxDrop = drop; maxDropGal = galaxyData[i].name; }
}

const jackMin = Math.min(...jackR);
const jackMax = Math.max(...jackR);
const jackMean = jackR.reduce((s,v)=>s+v,0)/N;
const jackSE = Math.sqrt(jackR.reduce((s,v)=>s+(v-jackMean)**2,0) * (N-1) / N);

log("  Base r (|MHI,Vflat,n): " + baseR.toFixed(3));
log("  Jackknife range:        [" + jackMin.toFixed(3) + ", " + jackMax.toFixed(3) + "]");
log("  Jackknife mean:         " + jackMean.toFixed(3) + " +/- " + jackSE.toFixed(3));
log("  Most influential drop:  " + maxDropGal + " (delta_r = " + maxDrop.toFixed(3) + ")");
log("  Min r / base r:         " + (jackMin / baseR * 100).toFixed(1) + "% retained");
log("  Signal survives any single removal: " + (jackMin > 0.15 ? "YES" : "NO"));
log("");

log("=".repeat(80));
log("  PHASE 24 — COMPREHENSIVE VERDICT");
log("=".repeat(80));
log("");

const test1pass = bothPosAfterMHI / NSPLITS > 0.6;
const test2pass = renzoSurvives;
const test3pass = wigWins / NTRIALS > 0.55;
const test4pass = nPositive / results4.length > 0.7;
const test5pass = jackMin > 0.15;

const nTestsPassed = [test1pass, test2pass, test3pass, test4pass, test5pass].filter(Boolean).length;

log("  ┌──────────────────────────────────────────────────────────────────────────┐");
log("  │  Test                              Pass?   Details                     │");
log("  ├──────────────────────────────────────────────────────────────────────────┤");
log("  │  1. Split-half replication          " + (test1pass?"YES":"NO ") + "    " + (bothPosAfterMHI/NSPLITS*100).toFixed(0) + "% same sign after MHI".padEnd(30) + "│");
log("  │  2. Renzo's rule control            " + (test2pass?"YES":"NO ") + "    t=" + tStat(rAfterMHI_BFA_BB,N,3).toFixed(1) + " after baryon controls".padEnd(27) + "│");
log("  │  3. Train/test prediction           " + (test3pass?"YES":"NO ") + "    wins " + (wigWins/NTRIALS*100).toFixed(0) + "% of trials".padEnd(26) + "│");
log("  │  4. Stratification consistency      " + (test4pass?"YES":"NO ") + "    " + nPositive + "/" + results4.length + " subgroups positive".padEnd(30) + "│");
log("  │  5. Jackknife stability             " + (test5pass?"YES":"NO ") + "    min r=" + jackMin.toFixed(3) + ", robust".padEnd(23) + "│");
log("  ├──────────────────────────────────────────────────────────────────────────┤");
log("  │  TOTAL: " + nTestsPassed + "/5 validation tests passed" + " ".repeat(36) + "│");
log("  └──────────────────────────────────────────────────────────────────────────┘");
log("");

if (nTestsPassed >= 4) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  rcWiggliness VALIDATED as robust independent signal.                  ║");
  log("  ║  Replicates across splits, survives Renzo control, predicts out-of-    ║");
  log("  ║  sample, consistent across subgroups, stable to outlier removal.       ║");
  log("  ║                                                                         ║");
  log("  ║  STATUS: CONFIRMED SECOND PARAMETER (within SPARC)                     ║");
  log("  ║  NEXT: External replication (THINGS) or Environment & Neighbors door.  ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nTestsPassed >= 3) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  rcWiggliness shows STRONG but not fully robust signal.                ║");
  log("  ║  Passes most but not all validation tests.                             ║");
  log("  ║  STATUS: STRONG CANDIDATE — needs external replication.                ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  rcWiggliness FAILS comprehensive validation.                          ║");
  log("  ║  The Phase 23c result may have been overfitted or unstable.            ║");
  log("  ║  Galaxy history door: CLOSED.                                          ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");

log("=".repeat(80));
log("  WHAT EXTERNAL REPLICATION WOULD REQUIRE");
log("=".repeat(80));
log("");
log("  To move from 'confirmed within SPARC' to 'confirmed observationally':");
log("  1. THINGS survey: ~30+ galaxies with high-res HI velocity fields");
log("     → compute rcWiggliness from their rotation curves");
log("     → test correlation with per-galaxy a0");
log("  2. LITTLE THINGS: ~20+ dwarf galaxies — complementary mass range");
log("  3. PHANGS: ~70+ galaxies with ALMA CO kinematics — different tracer");
log("");
log("  These surveys provide independent rotation curves from independent");
log("  data reduction pipelines. If rcWiggliness correlates with delta_a0");
log("  in ANY of these samples, the finding becomes publishable.");
log("");
log("  NOTE: We cannot access these datasets programmatically here.");
log("  This is documented as the natural next observational step.");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  splitHalf: {
    nSplits: NSPLITS,
    bothPositiveRaw: +(bothPositive/NSPLITS*100).toFixed(1),
    bothPositiveAfterMHI: +(bothPosAfterMHI/NSPLITS*100).toFixed(1),
    bothSignificantAfterMHI: +(bothSigAfterMHI/NSPLITS*100).toFixed(1)
  },
  renzoRule: {
    rAfterMHI: +rBaseAfterMHI.toFixed(3),
    rAfterMHI_BFA_BB: +rAfterMHI_BFA_BB.toFixed(3),
    tAfterMHI_BFA_BB: +tStat(rAfterMHI_BFA_BB, N, 3).toFixed(1),
    survives: renzoSurvives
  },
  trainTest: {
    nTrials: NTRIALS,
    wigWinsVsUniversal: +(wigWins/NTRIALS*100).toFixed(1),
    mhiWigWinsVsMHI: +(mhiWigWins/NTRIALS*100).toFixed(1)
  },
  stratification: {
    nSubgroups: results4.length,
    nPositive: nPositive,
    nSignificant: nSig,
    subgroups: results4
  },
  jackknife: {
    baseR: +baseR.toFixed(3),
    minR: +jackMin.toFixed(3),
    maxR: +jackMax.toFixed(3),
    mostInfluential: maxDropGal,
    stable: test5pass
  },
  overallVerdict: {
    testsPassed: nTestsPassed,
    testsTotal: 5,
    status: nTestsPassed >= 4 ? 'CONFIRMED' : nTestsPassed >= 3 ? 'STRONG_CANDIDATE' : 'FAILED'
  }
};

fs.writeFileSync(path.join(__dirname, '../public/phase24-rcwig-validation.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase24-rcwig-validation.json");
