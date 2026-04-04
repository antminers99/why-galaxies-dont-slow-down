#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "22.0.0";
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

  const gasFraction = sp.MHI / (sp.MHI + L * 0.5);
  const logGasFrac = Math.log10(gasFraction);
  const logSBdisk = sp.SBdisk > 0 ? Math.log10(sp.SBdisk) : 2;
  const logSBeff = sp.SBeff > 0 ? Math.log10(sp.SBeff) : 2;
  const concentration = sp.Reff > 0 && sp.Rdisk > 0 ? sp.Reff / sp.Rdisk : 1;
  const logConcentration = Math.log10(concentration);
  const gasExtent = sp.RHI > 0 && sp.Rdisk > 0 ? sp.RHI / sp.Rdisk : 3;

  const recentToPast = logMHI_L + 0.5 * sp.T / 10;

  const burstiness = Math.abs(logSBdisk - logSBeff);

  const youngFracProxy = gasFraction * (sp.T > 5 ? 1.5 : sp.T > 3 ? 1.0 : 0.5);

  const depletionTime = sp.MHI * 1e9 / (Math.pow(10, logSBdisk) * sp.Rdisk * sp.Rdisk * Math.PI * 1e-3 + 1e-5);
  const logDepletionTime = Math.log10(Math.max(depletionTime, 0.01));

  const massBuildupRate = logMstar - logMHI - Math.log10(Math.max(sp.D, 1)) * 0.3;

  const diskSettling = logSBdisk - logConcentration;

  const SFE = logMstar - logMHI;

  const gasConsumption = logMstar - logGasFrac;

  const birthParameter = logMHI_L + 0.4 * sp.T / 10 - 0.2 * logMstar / 10;

  const baryonConversion = logMstar - Math.log10(sp.MHI * 1e9 + L * 0.5e9);

  const diskAge = -sp.T * 0.1 + logSBeff * 0.3 - logMHI_L * 0.4 + logConcentration * 0.2;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q, fD: sp.fD,
    n: pts.length, etaRot,
    logMHI, logL36: Math.log10(L), logMstar, logMHI_L, logGasFrac,
    logSBdisk, logSBeff, logConcentration, gasExtent,
    Vflat: sp.Vflat,
    recentToPast, burstiness, youngFracProxy, logDepletionTime,
    massBuildupRate, diskSettling, SFE, gasConsumption,
    birthParameter, baryonConversion, diskAge,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

log("");
log("=".repeat(80));
log("  PHASE 22: STAR FORMATION HISTORY (SFH)");
log("  Does the PATTERN of star formation over time explain a0?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  KEY DISTINCTION FROM PHASE 19 (sSFR):");
sep();
log("");
log("  Phase 19 asked: is the galaxy forming stars NOW?");
log("  Phase 22 asks:  HOW did it form stars OVER TIME?");
log("");
log("  SFH captures:");
log("    - recent-to-past SFR ratio (declining vs rising history)");
log("    - burstiness (smooth vs bursty formation)");
log("    - star formation efficiency (how much gas → stars)");
log("    - depletion timescale (how fast gas is consumed)");
log("    - birth parameter (overall evolutionary state)");
log("    - disk settling age (how dynamically settled is the disk)");
log("");

log("  SFH PROXIES CONSTRUCTED:");
sep();
log("");
log("  1. recentToPast:     log(MHI/L) + 0.5*T/10     recent/past SFR ratio");
log("  2. burstiness:       |logSBdisk - logSBeff|     SB profile irregularity");
log("  3. youngFracProxy:   gasFrac * T-weight         young stellar fraction");
log("  4. logDepletionTime: MHI / (SBdisk * Rdisk^2)   gas consumption timescale");
log("  5. massBuildupRate:  logM* - logMHI - 0.3*logD  mass assembly rate");
log("  6. diskSettling:     logSBdisk - logConc         disk dynamical maturity");
log("  7. SFE:              logM* - logMHI              star formation efficiency");
log("  8. gasConsumption:   logM* - logGasFrac          gas consumption fraction");
log("  9. birthParameter:   logMHI/L + 0.4*T/10 - ...  evolutionary stage (b)");
log("  10. baryonConversion: logM* - log(Mbar)          baryon conversion fraction");
log("  11. diskAge:          composite disk age          disk formation time proxy");
log("");

const allFeatures = [
  { name: 'recentToPast', values: galaxyData.map(g => g.recentToPast), type: 'SFH' },
  { name: 'burstiness', values: galaxyData.map(g => g.burstiness), type: 'SFH' },
  { name: 'youngFracProxy', values: galaxyData.map(g => g.youngFracProxy), type: 'SFH' },
  { name: 'logDepletion', values: galaxyData.map(g => g.logDepletionTime), type: 'SFH' },
  { name: 'massBuildup', values: galaxyData.map(g => g.massBuildupRate), type: 'SFH' },
  { name: 'diskSettling', values: galaxyData.map(g => g.diskSettling), type: 'SFH' },
  { name: 'SFE', values: galaxyData.map(g => g.SFE), type: 'SFH' },
  { name: 'gasConsumption', values: galaxyData.map(g => g.gasConsumption), type: 'SFH' },
  { name: 'birthParameter', values: galaxyData.map(g => g.birthParameter), type: 'SFH' },
  { name: 'baryonConversion', values: galaxyData.map(g => g.baryonConversion), type: 'SFH' },
  { name: 'diskAge', values: galaxyData.map(g => g.diskAge), type: 'SFH' },
  { name: 'log(MHI)', values: galaxyData.map(g => g.logMHI), type: 'reference' },
  { name: 'log(Mstar)', values: galaxyData.map(g => g.logMstar), type: 'reference' },
  { name: 'eta_rot', values: galaxyData.map(g => g.etaRot), type: 'reference' },
  { name: 'T (Hubble)', values: galaxyData.map(g => g.T), type: 'reference' },
  { name: 'n_points', values: galaxyData.map(g => g.n), type: 'quality' },
  { name: 'log(D)', values: galaxyData.map(g => Math.log10(g.D)), type: 'quality' },
  { name: 'inc', values: galaxyData.map(g => g.inc), type: 'quality' },
];

log("  STEP 1: CORRELATIONS — SFH PROXIES vs delta_a0");
sep();
log("");

const corrs = allFeatures.map(f => {
  const r = corrWith(f.values, deltaA0);
  const t = Math.abs(r) * Math.sqrt(N - 2) / Math.sqrt(1 - r * r + 1e-10);
  return { ...f, r, t };
}).sort((a, b) => Math.abs(b.r) - Math.abs(a.r));

log("  ┌───────────────────────────────────────────────────────────────────────┐");
log("  │  Feature            Type              r        |t|    sig           │");
log("  ├───────────────────────────────────────────────────────────────────────┤");
for (const c of corrs) {
  const sig = c.t > 2.58 ? "***" : c.t > 2.0 ? "**" : c.t > 1.65 ? "*" : "";
  const marker = c.type === 'SFH' ? " [SFH]" : "";
  log("  │  " + c.name.padEnd(19) + c.type.padEnd(12) +
    c.r.toFixed(3).padStart(7) + c.t.toFixed(1).padStart(8) + "    " + (sig + marker).padEnd(10) + "│");
}
log("  └───────────────────────────────────────────────────────────────────────┘");
log("");

const sfhCorrs = corrs.filter(c => c.type === 'SFH');
const bestSFH = sfhCorrs[0];
log("  BEST SFH PROXY: " + bestSFH.name + " (r=" + bestSFH.r.toFixed(3) + ", t=" + bestSFH.t.toFixed(1) + ")");
log("");

log("  STEP 2: INDEPENDENCE TEST — DOES SFH ADD TO MHI?");
sep();
log("");

const mhiVals = galaxyData.map(g => g.logMHI);
const mstarVals = galaxyData.map(g => g.logMstar);

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

for (const mc of sfhCorrs) {
  const rRaw = corrWith(mc.values, deltaA0);
  const rAfterMHI = corrWith(mc.values, residAfterMHI);
  const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);
  const rAfterBoth = corrWith(mc.values, residAfterBoth);
  const tAfterBoth = Math.abs(rAfterBoth) * Math.sqrt(N-4) / Math.sqrt(1-rAfterBoth**2+1e-10);
  const status = Math.abs(rAfterMHI) < 0.10 ? "ABSORBED" : tAfterMHI >= 1.65 ? "INDEPENDENT" : "WEAK";
  log("  " + mc.name.padEnd(16) +
    " raw=" + rRaw.toFixed(3) +
    " |MHI: r=" + rAfterMHI.toFixed(3) + " t=" + tAfterMHI.toFixed(1) +
    " |both: r=" + rAfterBoth.toFixed(3) + " t=" + tAfterBoth.toFixed(1) +
    "  " + status);
}
log("");

const independentSFH = sfhCorrs.filter(mc => {
  const rAfterMHI = corrWith(mc.values, residAfterMHI);
  return Math.abs(rAfterMHI) >= 0.10;
});

if (independentSFH.length > 0) {
  log("  INDEPENDENT SFH signals (after MHI):");
  for (const mc of independentSFH) {
    const rAfterMHI = corrWith(mc.values, residAfterMHI);
    const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);
    log("    " + mc.name + ": r=" + rAfterMHI.toFixed(3) + " (t=" + tAfterMHI.toFixed(1) + ")");
  }
} else {
  log("  NO SFH proxy has independent signal after MHI.");
}
log("");

log("  STEP 3: PERMUTATION TEST — BEST SFH PROXY");
sep();
log("");

const NPERMS = 5000;
const obsR = Math.abs(corrWith(bestSFH.values, deltaA0));
let permCount = 0;
const shuffled = [...deltaA0];
for (let p = 0; p < NPERMS; p++) {
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  if (Math.abs(corrWith(bestSFH.values, shuffled)) >= obsR) permCount++;
}
const permP = permCount / NPERMS;
log("  " + bestSFH.name + ": |r|=" + obsR.toFixed(3) + " perm_p=" + permP.toFixed(4) +
  (permP < 0.05 ? " *" : "") + (permP < 0.01 ? "*" : ""));
log("");

if (independentSFH.length > 0) {
  log("  Permutation for independent proxies:");
  for (const mc of independentSFH.slice(0, 3)) {
    const obsR2 = Math.abs(corrWith(mc.values, residAfterMHI));
    let cnt = 0;
    const sh = [...residAfterMHI];
    for (let p = 0; p < NPERMS; p++) {
      for (let i = sh.length - 1; i > 0; i--) {
        const j = Math.floor(Math.random() * (i + 1));
        [sh[i], sh[j]] = [sh[j], sh[i]];
      }
      if (Math.abs(corrWith(mc.values, sh)) >= obsR2) cnt++;
    }
    log("  " + mc.name + " (vs MHI residuals): |r|=" + obsR2.toFixed(3) + " perm_p=" + (cnt/NPERMS).toFixed(4));
  }
  log("");
}

log("  STEP 4: LOO CROSS-VALIDATION — SFH MODELS");
sep();
log("");

function getVal(gIdx, fname) {
  const f = allFeatures.find(x => x.name === fname);
  return f ? f.values[gIdx] : 0;
}

const topSFH = sfhCorrs.slice(0, 3).map(c => c.name);
const modelDefs = [
  { name: 'M0: Universal a0', features: [] },
  { name: 'M_SFH1: best SFH', features: [bestSFH.name] },
  { name: 'M_SFH2: top-2 SFH', features: topSFH.slice(0, 2) },
  { name: 'M_SFH3: top-3 SFH', features: topSFH },
  { name: 'M_SFH_all: all SFH', features: sfhCorrs.map(c => c.name) },
  { name: 'M_MHI: gas only', features: ['log(MHI)'] },
  { name: 'M_SFH+MHI: best+gas', features: [bestSFH.name, 'log(MHI)'] },
  { name: 'M_kitchen: everything', features: allFeatures.map(f => f.name) },
  { name: 'M6: Per-galaxy', features: ['__free__'] },
];

const cvResults = [];
for (const model of modelDefs) {
  let totalSS = 0, totalN = 0;

  if (model.features[0] === '__free__') {
    for (let i = 0; i < N; i++) {
      totalSS += predictSS(galaxyData[i].a0, galaxyData[i].pts);
      totalN += galaxyData[i].pts.length;
    }
    cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), k: N });
    continue;
  }

  if (model.features.length === 0) {
    for (let i = 0; i < N; i++) {
      const trainY = allLogA0.filter((_, j) => j !== i);
      const mu = trainY.reduce((s,v) => s+v, 0) / trainY.length;
      totalSS += predictSS(Math.pow(10, mu), galaxyData[i].pts);
      totalN += galaxyData[i].pts.length;
    }
    cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), k: 1 });
    continue;
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

const bestSFHGap = Math.max(...cvResults.filter(r => r.name.includes('M_SFH')).map(r => gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0));
const mhiGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI: gas only').cvRMS) / gap * 100 : 0;

log("  STEP 5: QUARTILE ANALYSIS");
sep();
log("");

const bestVals = bestSFH.values;
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

log("  Split by " + bestSFH.name + ":");
log("  Q1: n=" + q1Idx.length + ", a0=" + Math.round(Math.pow(10, q1mean)) + ", tau=" + q1DL.tau.toFixed(3));
log("  Q4: n=" + q4Idx.length + ", a0=" + Math.round(Math.pow(10, q4mean)) + ", tau=" + q4DL.tau.toFixed(3));
log("  Split: " + splitDex.toFixed(3) + " dex (t=" + splitT.toFixed(2) + ")");
log("  Full tau: " + fullDL.tau.toFixed(3));
log("");

log("  STEP 6: STRICT SUCCESS CRITERIA");
sep();
log("");

const rAfterMHI_best = corrWith(bestSFH.values, residAfterMHI);
const tAfterMHI_best = Math.abs(rAfterMHI_best) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI_best**2+1e-10);

const tauReduction = Math.min(q1DL.tau, q4DL.tau) / fullDL.tau;

const criteria = [
  { name: 'Clear correlation |t|>=2', pass: bestSFH.t >= 2.0, val: 't=' + bestSFH.t.toFixed(1) },
  { name: 'Independent of MHI', pass: tAfterMHI_best >= 1.65, val: 'r=' + rAfterMHI_best.toFixed(3) + ' t=' + tAfterMHI_best.toFixed(1) },
  { name: 'LOO-CV > M0', pass: bestSFHGap > 0, val: bestSFHGap.toFixed(1) + '%' },
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
log("  PHASE 22 — CONCLUSIONS");
log("=".repeat(80));
log("");

log("  1. BEST SFH PROXY: " + bestSFH.name);
log("     r = " + bestSFH.r.toFixed(3) + " (t=" + bestSFH.t.toFixed(1) + ")");
log("");

log("  2. INDEPENDENCE FROM MHI: " + (tAfterMHI_best >= 1.65 ? "YES" : "NO"));
log("     After MHI: r=" + rAfterMHI_best.toFixed(3) + " (t=" + tAfterMHI_best.toFixed(1) + ")");
log("");

log("  3. LOO-CV GAP CLOSURE:");
log("     Best SFH model:  " + bestSFHGap.toFixed(1) + "%");
log("     MHI alone:       " + mhiGap.toFixed(1) + "%");
log("     Phase 21 metal:  4.2%");
log("     Phase 20 age:    0.9%");
log("     Phase 19 sSFR:   -1.8%");
log("     Phase 15 baryon: 14-16%");
log("");

log("  4. SUCCESS CRITERIA: " + nPassed + "/4 → " + (nPassed >= 3 ? "PASS" : "FAIL"));
log("");

if (nPassed >= 3) {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  SFH proxies show promising signal. Star formation HISTORY         ║");
  log("  ║  (not just current rate) contributes to a0 variation.              ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  SFH proxies FAIL the strict success criteria.                     ║");
  log("  ║  Like sSFR, age, and metallicity — all galaxy history variables    ║");
  log("  ║  are either too weak or absorbed by MHI.                           ║");
  log("  ║                                                                     ║");
  log("  ║  THE ENTIRE 'GALAXY HISTORY' DOOR IS NOW CLOSED.                   ║");
  log("  ║  MHI remains the dominant and unexplained correlate.               ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  GALAXY HISTORY — FINAL SCORECARD:");
log("  ┌────────────────────────────────────────────────────────────────────────┐");
log("  │  Variable       Status          Gap Closed  Independent of MHI?      │");
log("  ├────────────────────────────────────────────────────────────────────────┤");
log("  │  sSFR           FAILED          -1.8%       partially                │");
log("  │  Stellar age    FAILED           0.9%       NO                       │");
log("  │  Metallicity    FAILED           4.2%       NO (completely absorbed) │");
log("  │  SFH            " + (nPassed >= 3 ? "PARTIAL " : "FAILED  ") + "        " + bestSFHGap.toFixed(1).padStart(5) + "%       " + (tAfterMHI_best >= 1.65 ? "YES" : "NO").padEnd(25) + "│");
log("  │  Merger proxy   SKIPPED         —           —                        │");
log("  ├────────────────────────────────────────────────────────────────────────┤");
log("  │  DOOR STATUS:   " + (nPassed >= 3 ? "PARTIALLY OPEN" : "CLOSED") + " ".repeat(nPassed >= 3 ? 38 : 46) + "│");
log("  └────────────────────────────────────────────────────────────────────────┘");
log("");

log("  THE MHI QUESTION:");
log("  Every galaxy history variable is either too weak or absorbed by MHI.");
log("  MHI (log gas mass) remains the dominant correlate of delta_a0.");
log("  The next investigation should ask: WHY does MHI predict a0?");
log("  Candidate answers:");
log("    a) MHI affects mass modeling (gas disk contribution to gbar)");
log("    b) MHI traces halo properties (gas retention = halo potential)");
log("    c) MHI is a proxy for an unmeasured structural parameter");
log("    d) MHI relates to distance/selection biases");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  bestSFHProxy: { name: bestSFH.name, r: +bestSFH.r.toFixed(3), t: +bestSFH.t.toFixed(1) },
  correlations: corrs.map(c => ({ name: c.name, type: c.type, r: +c.r.toFixed(3), t: +c.t.toFixed(1) })),
  independenceAfterMHI: sfhCorrs.map(c => {
    const rA = corrWith(c.values, residAfterMHI);
    return { name: c.name, r: +rA.toFixed(3) };
  }),
  cvResults: cvResults.map(r => ({
    name: r.name, k: r.k, cvRMS: +r.cvRMS.toFixed(5),
    gapClosed: +(gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0).toFixed(1)
  })),
  bestSFHGapClosed: +bestSFHGap.toFixed(1),
  mhiGapClosed: +mhiGap.toFixed(1),
  permutationP: +permP.toFixed(4),
  successCriteria: { passed: nPassed, total: 4 },
  quartileSplit: {
    q1: { n: q1Idx.length, meanA0: Math.round(Math.pow(10, q1mean)), tau: +q1DL.tau.toFixed(3) },
    q4: { n: q4Idx.length, meanA0: Math.round(Math.pow(10, q4mean)), tau: +q4DL.tau.toFixed(3) },
    splitDex: +splitDex.toFixed(3), splitT: +splitT.toFixed(2)
  },
  doorStatus: nPassed >= 3 ? 'partially open' : 'CLOSED',
  galaxyHistoryScorecard: {
    sSFR: 'FAILED', stellarAge: 'FAILED', metallicity: 'FAILED',
    SFH: nPassed >= 3 ? 'PARTIAL' : 'FAILED',
    mergerProxy: 'SKIPPED'
  }
};

fs.writeFileSync(path.join(__dirname, '../public/phase22-sfh.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase22-sfh.json");
