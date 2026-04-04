#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "23.0.0";
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

  let etaRot = 0;
  if (pts.length >= 4) {
    const oh = pts.slice(Math.floor(pts.length / 2)), ih = pts.slice(0, Math.floor(pts.length / 2));
    if (oh.length >= 2 && ih.length >= 2) {
      etaRot = (oh[oh.length-1].log_g_obs-oh[0].log_g_obs)/(oh[oh.length-1].log_g_bar-oh[0].log_g_bar+1e-10) -
               (ih[ih.length-1].log_g_obs-ih[0].log_g_obs)/(ih[ih.length-1].log_g_bar-ih[0].log_g_bar+1e-10);
    }
  }

  const logSBdisk = sp.SBdisk > 0 ? Math.log10(sp.SBdisk) : 2;
  const logSBeff = sp.SBeff > 0 ? Math.log10(sp.SBeff) : 2;
  const concentration = sp.Reff > 0 && sp.Rdisk > 0 ? sp.Reff / sp.Rdisk : 1;
  const logConcentration = Math.log10(concentration);
  const gasExtent = sp.RHI > 0 && sp.Rdisk > 0 ? sp.RHI / sp.Rdisk : 3;

  const sbAsymmetry = Math.abs(logSBdisk - logSBeff);
  const profileIrregularity = sp.SBeff > 0 && sp.SBdisk > 0 ?
    Math.abs(Math.log10(sp.SBeff / sp.SBdisk) - Math.log10(concentration)) : 0;

  const gasStarMismatch = Math.abs(gasExtent - 3.0);

  const rcScatter = fit.rms;

  let rcAsymmetry = 0;
  if (pts.length >= 6) {
    const half = Math.floor(pts.length / 2);
    const inner = pts.slice(0, half);
    const outer = pts.slice(half);
    const innerSlope = inner.length >= 2 ?
      (inner[inner.length-1].log_g_obs - inner[0].log_g_obs) / (inner[inner.length-1].log_g_bar - inner[0].log_g_bar + 1e-10) : 1;
    const outerSlope = outer.length >= 2 ?
      (outer[outer.length-1].log_g_obs - outer[0].log_g_obs) / (outer[outer.length-1].log_g_bar - outer[0].log_g_bar + 1e-10) : 1;
    rcAsymmetry = Math.abs(innerSlope - outerSlope);
  }

  const sizeRatio = sp.Reff > 0 && sp.Rdisk > 0 ? sp.Reff / sp.Rdisk : 1;
  const sizeAnomaly = Math.abs(Math.log10(sizeRatio) - Math.log10(1.678));

  const lopsidedness = sp.SBeff > 0 && sp.SBdisk > 0 ?
    Math.abs(logSBeff - logSBdisk - 0.5) : 0;

  const kinematicDisturbance = rcScatter * Math.sqrt(pts.length);

  const morphDisturbance = sbAsymmetry + sizeAnomaly * 0.5 + lopsidedness * 0.3;

  const totalDisturbance = morphDisturbance + kinematicDisturbance * 0.3 + gasStarMismatch * 0.2;

  const interactionProxy = gasExtent > 5 ? 1 : 0;
  const warpProxy = sp.inc > 80 ? sbAsymmetry * 1.5 : sbAsymmetry;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q, fD: sp.fD,
    n: pts.length, etaRot,
    logMHI, logL36: Math.log10(L), logMstar, logMHI_L,
    logSBdisk, logSBeff, logConcentration, gasExtent,
    Vflat: sp.Vflat,
    sbAsymmetry, profileIrregularity, gasStarMismatch,
    rcScatter, rcAsymmetry, sizeAnomaly, lopsidedness,
    kinematicDisturbance, morphDisturbance, totalDisturbance,
    interactionProxy, warpProxy,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

log("");
log("=".repeat(80));
log("  PHASE 23: MERGER / DISTURBANCE PROXIES");
log("  Does morphological or kinematic disturbance explain a0?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  THIS IS THE FINAL TEST IN THE 'GALAXY HISTORY' DOOR.");
log("  If it fails, the entire door closes.");
log("");

log("  DISTURBANCE PROXIES CONSTRUCTED:");
sep();
log("");
log("  Morphological disturbance:");
log("    1. sbAsymmetry:       |logSBdisk - logSBeff|");
log("    2. profileIrregularity: deviation from expected SB-concentration relation");
log("    3. sizeAnomaly:       |log(Reff/Rdisk) - log(1.678)| (exponential disk)");
log("    4. lopsidedness:      |logSBeff - logSBdisk - 0.5| (asymmetric light)");
log("    5. morphDisturbance:  composite morphological index");
log("");
log("  Kinematic disturbance:");
log("    6. rcScatter:         RMS residual of RAR fit (noisy RC)");
log("    7. rcAsymmetry:       |inner slope - outer slope| in RAR");
log("    8. kinematicDist:     rcScatter * sqrt(n_points)");
log("");
log("  Gas-structure disturbance:");
log("    9. gasStarMismatch:   |RHI/Rdisk - 3| (abnormal gas extent)");
log("    10. interactionProxy: flag for very extended gas (RHI/Rdisk > 5)");
log("    11. warpProxy:        sbAsymmetry boosted for edge-on galaxies");
log("");
log("  Combined:");
log("    12. totalDisturbance: weighted sum of morph + kin + gas");
log("");

const allFeatures = [
  { name: 'sbAsymmetry', values: galaxyData.map(g => g.sbAsymmetry), type: 'morph' },
  { name: 'profileIrreg', values: galaxyData.map(g => g.profileIrregularity), type: 'morph' },
  { name: 'sizeAnomaly', values: galaxyData.map(g => g.sizeAnomaly), type: 'morph' },
  { name: 'lopsidedness', values: galaxyData.map(g => g.lopsidedness), type: 'morph' },
  { name: 'morphDisturb', values: galaxyData.map(g => g.morphDisturbance), type: 'morph' },
  { name: 'rcScatter', values: galaxyData.map(g => g.rcScatter), type: 'kinematic' },
  { name: 'rcAsymmetry', values: galaxyData.map(g => g.rcAsymmetry), type: 'kinematic' },
  { name: 'kinematicDist', values: galaxyData.map(g => g.kinematicDisturbance), type: 'kinematic' },
  { name: 'gasStarMismatch', values: galaxyData.map(g => g.gasStarMismatch), type: 'gas' },
  { name: 'warpProxy', values: galaxyData.map(g => g.warpProxy), type: 'morph' },
  { name: 'totalDisturb', values: galaxyData.map(g => g.totalDisturbance), type: 'combined' },
  { name: 'log(MHI)', values: galaxyData.map(g => g.logMHI), type: 'reference' },
  { name: 'log(Mstar)', values: galaxyData.map(g => g.logMstar), type: 'reference' },
  { name: 'eta_rot', values: galaxyData.map(g => g.etaRot), type: 'reference' },
  { name: 'T (Hubble)', values: galaxyData.map(g => g.T), type: 'reference' },
  { name: 'n_points', values: galaxyData.map(g => g.n), type: 'quality' },
  { name: 'log(D)', values: galaxyData.map(g => Math.log10(g.D)), type: 'quality' },
  { name: 'inc', values: galaxyData.map(g => g.inc), type: 'quality' },
];

log("  STEP 1: CORRELATIONS — DISTURBANCE PROXIES vs delta_a0");
sep();
log("");

const corrs = allFeatures.map(f => {
  const r = corrWith(f.values, deltaA0);
  const t = Math.abs(r) * Math.sqrt(N - 2) / Math.sqrt(1 - r * r + 1e-10);
  return { ...f, r, t };
}).sort((a, b) => Math.abs(b.r) - Math.abs(a.r));

const distTypes = ['morph','kinematic','gas','combined'];

log("  ┌───────────────────────────────────────────────────────────────────────┐");
log("  │  Feature            Type              r        |t|    sig           │");
log("  ├───────────────────────────────────────────────────────────────────────┤");
for (const c of corrs) {
  const sig = c.t > 2.58 ? "***" : c.t > 2.0 ? "**" : c.t > 1.65 ? "*" : "";
  const marker = distTypes.includes(c.type) ? " [D]" : "";
  log("  │  " + c.name.padEnd(19) + c.type.padEnd(12) +
    c.r.toFixed(3).padStart(7) + c.t.toFixed(1).padStart(8) + "    " + (sig + marker).padEnd(10) + "│");
}
log("  └───────────────────────────────────────────────────────────────────────┘");
log("");

const distCorrs = corrs.filter(c => distTypes.includes(c.type));
const bestDist = distCorrs[0];
log("  BEST DISTURBANCE PROXY: " + bestDist.name + " (r=" + bestDist.r.toFixed(3) + ", t=" + bestDist.t.toFixed(1) + ")");
log("");

log("  STEP 2: INDEPENDENCE TEST — DOES DISTURBANCE ADD TO MHI?");
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

const residAfterBoth = [];
{
  const X = galaxyData.map(g => [g.logMHI, g.logMstar]);
  const reg = linReg(X, deltaA0);
  for(let i=0;i<N;i++) residAfterBoth.push(deltaA0[i] - (reg.intercept + reg.coefs[0]*X[i][0] + reg.coefs[1]*X[i][1]));
}

const independentDist = [];
for (const mc of distCorrs) {
  const rRaw = corrWith(mc.values, deltaA0);
  const rAfterMHI = corrWith(mc.values, residAfterMHI);
  const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);
  const rAfterBoth = corrWith(mc.values, residAfterBoth);
  const tAfterBoth = Math.abs(rAfterBoth) * Math.sqrt(N-4) / Math.sqrt(1-rAfterBoth**2+1e-10);
  const status = tAfterMHI >= 1.65 ? "INDEPENDENT" : Math.abs(rAfterMHI) >= 0.10 ? "WEAK" : "ABSORBED";
  if (Math.abs(rAfterMHI) >= 0.10) independentDist.push({ ...mc, rAfterMHI, tAfterMHI });
  log("  " + mc.name.padEnd(16) +
    " raw=" + rRaw.toFixed(3) +
    " |MHI: r=" + rAfterMHI.toFixed(3) + " t=" + tAfterMHI.toFixed(1) +
    " |both: r=" + rAfterBoth.toFixed(3) + " t=" + tAfterBoth.toFixed(1) +
    "  " + status);
}
log("");

if (independentDist.length > 0) {
  log("  DISTURBANCE PROXIES WITH RESIDUAL SIGNAL:");
  for (const mc of independentDist) {
    log("    " + mc.name + ": r_afterMHI=" + mc.rAfterMHI.toFixed(3) + " (t=" + mc.tAfterMHI.toFixed(1) + ")");
  }
} else {
  log("  NO disturbance proxy retains signal after MHI.");
}
log("");

log("  STEP 3: PERMUTATION TEST");
sep();
log("");

const NPERMS = 5000;
for (const mc of distCorrs.slice(0, 3)) {
  const obsR = Math.abs(corrWith(mc.values, deltaA0));
  let cnt = 0;
  const sh = [...deltaA0];
  for (let p = 0; p < NPERMS; p++) {
    for (let i = sh.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [sh[i], sh[j]] = [sh[j], sh[i]];
    }
    if (Math.abs(corrWith(mc.values, sh)) >= obsR) cnt++;
  }
  log("  " + mc.name.padEnd(16) + " |r|=" + obsR.toFixed(3) + " perm_p=" + (cnt/NPERMS).toFixed(4) +
    (cnt/NPERMS < 0.05 ? " *" : ""));
}

if (independentDist.length > 0) {
  log("");
  log("  Permutation on MHI residuals:");
  for (const mc of independentDist.slice(0, 3)) {
    const obsR = Math.abs(mc.rAfterMHI);
    let cnt = 0;
    const sh = [...residAfterMHI];
    for (let p = 0; p < NPERMS; p++) {
      for (let i = sh.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [sh[i], sh[j]] = [sh[j], sh[i]];
      }
      if (Math.abs(corrWith(mc.values, sh)) >= obsR) cnt++;
    }
    log("  " + mc.name.padEnd(16) + " |r|=" + obsR.toFixed(3) + " perm_p=" + (cnt/NPERMS).toFixed(4));
  }
}
log("");

log("  STEP 4: LOO CROSS-VALIDATION");
sep();
log("");

function getVal(gIdx, fname) {
  const f = allFeatures.find(x => x.name === fname);
  return f ? f.values[gIdx] : 0;
}

const topDist = distCorrs.slice(0, 3).map(c => c.name);
const modelDefs = [
  { name: 'M0: Universal a0', features: [] },
  { name: 'M_D1: best disturb', features: [bestDist.name] },
  { name: 'M_D2: top-2 disturb', features: topDist.slice(0, 2) },
  { name: 'M_D3: top-3 disturb', features: topDist },
  { name: 'M_Dall: all disturb', features: distCorrs.map(c => c.name) },
  { name: 'M_MHI: gas only', features: ['log(MHI)'] },
  { name: 'M_D+MHI: best+gas', features: [bestDist.name, 'log(MHI)'] },
  { name: 'M_kitchen: everything', features: allFeatures.map(f => f.name) },
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

log("  ┌──────────────────────────────────────────────────────────────────────────┐");
log("  │  Model                       k    CV-RMS    vs M0    gap-closed        │");
log("  ├──────────────────────────────────────────────────────────────────────────┤");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const gapClosed = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  │  " + r.name.padEnd(27) + r.k.toString().padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(9) + (gapClosed + "%").padStart(13) + "       │");
}
log("  └──────────────────────────────────────────────────────────────────────────┘");
log("");

const bestDistGap = Math.max(...cvResults.filter(r => r.name.includes('M_D')).map(r => gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0));
const mhiGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI: gas only').cvRMS) / gap * 100 : 0;

log("  STEP 5: QUARTILE ANALYSIS — DISTURBED vs UNDISTURBED");
sep();
log("");

const bestVals = bestDist.values;
const sortedIdx = [...Array(N).keys()].sort((a, b) => bestVals[a] - bestVals[b]);
const q1Idx = sortedIdx.slice(0, Math.floor(N / 4));
const q4Idx = sortedIdx.slice(Math.floor(3 * N / 4));

const q1LogA0 = q1Idx.map(i => allLogA0[i]);
const q4LogA0 = q4Idx.map(i => allLogA0[i]);
const q1mean = q1LogA0.reduce((s,v)=>s+v,0)/q1LogA0.length;
const q4mean = q4LogA0.reduce((s,v)=>s+v,0)/q4LogA0.length;
const q1sd = Math.sqrt(q1LogA0.reduce((s,v)=>s+(v-q1mean)**2,0)/(q1LogA0.length-1));
const q4sd = Math.sqrt(q4LogA0.reduce((s,v)=>s+(v-q4mean)**2,0)/(q4LogA0.length-1));
const q1DL = dlMeta(q1LogA0, q1Idx.map(i => galaxyData[i].se));
const q4DL = dlMeta(q4LogA0, q4Idx.map(i => galaxyData[i].se));
const fullDL = dlMeta(allLogA0, galaxyData.map(g => g.se));

const splitDex = Math.abs(q1mean - q4mean);
const splitT = splitDex / Math.sqrt(q1sd**2/q1Idx.length + q4sd**2/q4Idx.length);
const tauReduction = Math.min(q1DL.tau, q4DL.tau) / fullDL.tau;

log("  Split by " + bestDist.name + ":");
log("  Q1 (least disturbed): n=" + q1Idx.length + ", a0=" + Math.round(Math.pow(10, q1mean)) + ", tau=" + q1DL.tau.toFixed(3));
log("  Q4 (most disturbed):  n=" + q4Idx.length + ", a0=" + Math.round(Math.pow(10, q4mean)) + ", tau=" + q4DL.tau.toFixed(3));
log("  Split: " + splitDex.toFixed(3) + " dex (t=" + splitT.toFixed(2) + ")");
log("  Full tau: " + fullDL.tau.toFixed(3));
log("  Within-quartile tau reduction: " + ((1-tauReduction)*100).toFixed(1) + "%");
log("");

log("  STEP 6: STRICT SUCCESS CRITERIA");
sep();
log("");

const rAfterMHI_best = corrWith(bestDist.values, residAfterMHI);
const tAfterMHI_best = Math.abs(rAfterMHI_best) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI_best**2+1e-10);

const criteria = [
  { name: 'Clear correlation |t|>=2', pass: bestDist.t >= 2.0, val: 't=' + bestDist.t.toFixed(1) },
  { name: 'Independent of MHI', pass: tAfterMHI_best >= 1.65, val: 'r=' + rAfterMHI_best.toFixed(3) + ' t=' + tAfterMHI_best.toFixed(1) },
  { name: 'LOO-CV > M0', pass: bestDistGap > 0, val: bestDistGap.toFixed(1) + '%' },
  { name: 'tau reduction > 10%', pass: (1-tauReduction)*100 > 10, val: ((1-tauReduction)*100).toFixed(1) + '%' },
];

const nPassed = criteria.filter(c => c.pass).length;

log("  ┌─────────────────────────────────────────────────────────────────┐");
log("  │  Criterion              Pass?    Value                        │");
log("  ├─────────────────────────────────────────────────────────────────┤");
for (const c of criteria) {
  log("  │  " + c.name.padEnd(24) + (c.pass ? " YES " : " NO  ") + "  " + c.val.padEnd(28) + "│");
}
log("  ├─────────────────────────────────────────────────────────────────┤");
log("  │  TOTAL: " + nPassed + "/4 criteria met" + " ".repeat(37) + "│");
log("  └─────────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  PHASE 23 — CONCLUSIONS");
log("=".repeat(80));
log("");

log("  1. BEST DISTURBANCE PROXY: " + bestDist.name);
log("     r = " + bestDist.r.toFixed(3) + " (t=" + bestDist.t.toFixed(1) + ")");
log("");

log("  2. INDEPENDENCE FROM MHI: " + (tAfterMHI_best >= 1.65 ? "YES" : "NO"));
log("     After MHI: r=" + rAfterMHI_best.toFixed(3) + " (t=" + tAfterMHI_best.toFixed(1) + ")");
log("");

log("  3. LOO-CV GAP CLOSURE:");
log("     Best disturb model:  " + bestDistGap.toFixed(1) + "%");
log("     MHI alone:           " + mhiGap.toFixed(1) + "%");
log("     Phase 22 SFH:        4.0%");
log("     Phase 21 metal:      4.2%");
log("     Phase 20 age:        0.9%");
log("     Phase 19 sSFR:       -1.8%");
log("");

log("  4. SUCCESS CRITERIA: " + nPassed + "/4 → " + (nPassed >= 3 ? "PASS" : "FAIL"));
log("");

if (nPassed >= 3) {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  Disturbance/merger proxies show promise.                           ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  Disturbance/merger proxies FAIL.                                   ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}
log("");

log("=".repeat(80));
log("  GALAXY HISTORY DOOR — FINAL FINAL SCORECARD");
log("=".repeat(80));
log("");
log("  ┌────────────────────────────────────────────────────────────────────────┐");
log("  │  Variable       Status          Gap Closed  Independent of MHI?      │");
log("  ├────────────────────────────────────────────────────────────────────────┤");
log("  │  sSFR           FAILED          -1.8%       partially                │");
log("  │  Stellar age    FAILED           0.9%       NO                       │");
log("  │  Metallicity    FAILED           4.2%       NO (absorbed)            │");
log("  │  SFH            FAILED           4.0%       NO (absorbed)            │");
log("  │  Merger/disturb " + (nPassed >= 3 ? "PARTIAL " : "FAILED  ") + "        " + bestDistGap.toFixed(1).padStart(5) + "%       " + (tAfterMHI_best >= 1.65 ? "YES" : "NO").padEnd(25) + "│");
log("  ├────────────────────────────────────────────────────────────────────────┤");

const allFailed = nPassed < 3;
if (allFailed) {
  log("  │                                                                      │");
  log("  │  ══════════════════════════════════════════════════════════════════  │");
  log("  │  GALAXY HISTORY DOOR: OFFICIALLY CLOSED                             │");
  log("  │  All 5 variables tested. None pass strict criteria.                 │");
  log("  │  MHI remains the dominant, unexplained correlate.                   │");
  log("  │  NEXT DOOR: Environment & Neighbors                                │");
  log("  │  ══════════════════════════════════════════════════════════════════  │");
}
log("  └────────────────────────────────────────────────────────────────────────┘");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  bestDistProxy: { name: bestDist.name, r: +bestDist.r.toFixed(3), t: +bestDist.t.toFixed(1) },
  correlations: corrs.map(c => ({ name: c.name, type: c.type, r: +c.r.toFixed(3), t: +c.t.toFixed(1) })),
  independenceAfterMHI: distCorrs.map(c => {
    const rA = corrWith(c.values, residAfterMHI);
    return { name: c.name, r: +rA.toFixed(3) };
  }),
  cvResults: cvResults.map(r => ({
    name: r.name, k: r.k, cvRMS: +r.cvRMS.toFixed(5),
    gapClosed: +(gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0).toFixed(1)
  })),
  bestDistGapClosed: +bestDistGap.toFixed(1),
  mhiGapClosed: +mhiGap.toFixed(1),
  successCriteria: { passed: nPassed, total: 4 },
  quartileSplit: {
    q1: { n: q1Idx.length, meanA0: Math.round(Math.pow(10, q1mean)), tau: +q1DL.tau.toFixed(3) },
    q4: { n: q4Idx.length, meanA0: Math.round(Math.pow(10, q4mean)), tau: +q4DL.tau.toFixed(3) },
    splitDex: +splitDex.toFixed(3), splitT: +splitT.toFixed(2)
  },
  doorStatus: allFailed ? 'CLOSED' : 'PARTIALLY OPEN',
  galaxyHistoryScorecard: {
    sSFR: 'FAILED', stellarAge: 'FAILED', metallicity: 'FAILED',
    SFH: 'FAILED', mergerDisturb: nPassed >= 3 ? 'PARTIAL' : 'FAILED'
  }
};

fs.writeFileSync(path.join(__dirname, '../public/phase23-merger-disturbance.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase23-merger-disturbance.json");
