#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "19.0.0";
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

  let etaRot = 0;
  if (pts.length >= 4) {
    const oh = pts.slice(Math.floor(pts.length / 2)), ih = pts.slice(0, Math.floor(pts.length / 2));
    if (oh.length >= 2 && ih.length >= 2) {
      etaRot = (oh[oh.length-1].log_g_obs-oh[0].log_g_obs)/(oh[oh.length-1].log_g_bar-oh[0].log_g_bar+1e-10) -
               (ih[ih.length-1].log_g_obs-ih[0].log_g_obs)/(ih[ih.length-1].log_g_bar-ih[0].log_g_bar+1e-10);
    }
  }

  const logMstar = Math.log10(L * 0.5e9);
  const logSFR_proxy = Math.log10(sp.MHI * 1e9) - 0.32 * logMstar + 0.5 * Math.log10(sp.SBdisk > 0 ? sp.SBdisk : 100);
  const logsSFR_proxy = logSFR_proxy - logMstar;

  const gasFraction = sp.MHI / (sp.MHI + L * 0.5);
  const logGasFrac = Math.log10(gasFraction);
  const logMHI_L = Math.log10(sp.MHI / L);

  const compactness = sp.SBdisk > 0 ? Math.log10(sp.SBdisk) : 2;
  const gasExtent = sp.RHI > 0 && sp.Rdisk > 0 ? sp.RHI / sp.Rdisk : 3;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q, fD: sp.fD,
    Vmax: gal.Vmax, n: pts.length, etaRot,
    logMHI: Math.log10(sp.MHI), logL36: Math.log10(L),
    logMstar, logMHI_L, logGasFrac, gasFraction,
    logsSFR_proxy, logSFR_proxy,
    SBdisk: sp.SBdisk, Rdisk: sp.Rdisk, Reff: sp.Reff,
    RHI: sp.RHI, gasExtent, compactness,
    Vflat: sp.Vflat,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

log("");
log("=".repeat(80));
log("  PHASE 19: GALAXY HISTORY — sSFR AS THE MISSING VARIABLE");
log("  Does star formation history explain the a0 variation?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  STAR FORMATION PROXIES FROM SPARC DATA:");
sep();
log("");
log("  Available proxies (no external data needed):");
log("    log(MHI/L)     gas-to-stellar ratio = gas fraction proxy");
log("    log(gasFrac)   MHI/(MHI + 0.5*L)   = gas mass fraction");
log("    log(sSFR)      composite proxy from MHI + SBdisk + Mstar");
log("    T (Hubble)     morphological type (late = more star-forming)");
log("    SBdisk          disk surface brightness (star formation surface density)");
log("    RHI/Rdisk      gas extent / stellar extent (gas reservoir indicator)");
log("");

log("  WHY sSFR COULD BE THE MISSING LINK:");
log("    log(MHI) was our best predictor (r=-0.388, Phase 14)");
log("    eta_rot was significant (r=0.334, Phase 10)");
log("    T was significant (r=0.278, Phase 10)");
log("    All three track star formation activity.");
log("    sSFR may be the common underlying variable.");
log("");

log("  STEP 1: CORRELATIONS WITH delta_a0");
sep();
log("");

const historyFeatures = [
  { name: 'log(MHI/L)', values: galaxyData.map(g => g.logMHI_L), type: 'sSFR proxy' },
  { name: 'log(gasFrac)', values: galaxyData.map(g => g.logGasFrac), type: 'sSFR proxy' },
  { name: 'log(sSFR)', values: galaxyData.map(g => g.logsSFR_proxy), type: 'sSFR composite' },
  { name: 'T (Hubble)', values: galaxyData.map(g => g.T), type: 'morphology' },
  { name: 'log(MHI)', values: galaxyData.map(g => g.logMHI), type: 'gas mass' },
  { name: 'eta_rot', values: galaxyData.map(g => g.etaRot), type: 'RC shape' },
  { name: 'compactness', values: galaxyData.map(g => g.compactness), type: 'structure' },
  { name: 'gasExtent', values: galaxyData.map(g => g.gasExtent), type: 'gas structure' },
  { name: 'log(Mstar)', values: galaxyData.map(g => g.logMstar), type: 'stellar mass' },
  { name: 'log(D)', values: galaxyData.map(g => Math.log10(g.D)), type: 'distance' },
  { name: 'inc', values: galaxyData.map(g => g.inc), type: 'geometry' },
  { name: 'n_points', values: galaxyData.map(g => g.n), type: 'data quality' },
];

const corrs = historyFeatures.map(f => {
  const r = corrWith(f.values, deltaA0);
  const t = Math.abs(r) * Math.sqrt(N - 2) / Math.sqrt(1 - r * r + 1e-10);
  return { ...f, r, t };
}).sort((a, b) => Math.abs(b.r) - Math.abs(a.r));

log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Feature          Type            r        |t|    sig           │");
log("  ├───────────────────────────────────────────────────────────────────┤");
for (const c of corrs) {
  const sig = c.t > 2.58 ? "***" : c.t > 2.0 ? "**" : c.t > 1.65 ? "*" : "";
  const bar = (c.r > 0 ? "+" : "-").repeat(Math.round(Math.abs(c.r) * 20));
  log("  │  " + c.name.padEnd(17) + c.type.padEnd(16) +
    c.r.toFixed(3).padStart(7) + c.t.toFixed(1).padStart(8) + "    " + sig.padEnd(4) + "│");
}
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

const sSFRcorrs = corrs.filter(c => c.type.includes('sSFR'));
const bestsSFR = sSFRcorrs[0];

log("  BEST sSFR PROXY: " + bestsSFR.name + " (r=" + bestsSFR.r.toFixed(3) + ", t=" + bestsSFR.t.toFixed(1) + ")");
log("");

log("  STEP 2: IS sSFR THE COMMON THREAD?");
sep();
log("");

log("  Correlations among the clue variables:");
const clueVars = ['log(MHI/L)', 'log(MHI)', 'T (Hubble)', 'eta_rot'];
for (const v1 of clueVars) {
  const f1 = historyFeatures.find(f => f.name === v1);
  const row = clueVars.map(v2 => {
    if (v1 === v2) return " 1.00";
    const f2 = historyFeatures.find(f => f.name === v2);
    return corrWith(f1.values, f2.values).toFixed(2).padStart(5);
  });
  log("    " + v1.padEnd(14) + row.join("  "));
}
log("");

const residAftersSFR = [];
{
  const x = bestsSFR.values;
  const y = deltaA0;
  const n = x.length;
  const mx = x.reduce((s,v)=>s+v,0)/n, my = y.reduce((s,v)=>s+v,0)/n;
  let sxy=0,sxx=0;
  for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;}
  const b = sxy/sxx, a = my - b*mx;
  for(let i=0;i<n;i++) residAftersSFR.push(y[i] - (a + b*x[i]));
}

log("  After removing sSFR effect, do other clues still correlate?");
for (const v of clueVars) {
  if (v === bestsSFR.name) continue;
  const f = historyFeatures.find(ff => ff.name === v);
  const rBefore = corrWith(f.values, deltaA0);
  const rAfter = corrWith(f.values, residAftersSFR);
  log("    " + v.padEnd(14) + " r_before=" + rBefore.toFixed(3) + "  r_after_sSFR=" + rAfter.toFixed(3) +
    "  (" + (Math.abs(rAfter) < Math.abs(rBefore) * 0.5 ? "ABSORBED" : "INDEPENDENT") + ")");
}
log("");

log("  STEP 3: LOO CROSS-VALIDATION — sSFR-BASED MODELS");
sep();
log("");

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_sSFR: gas frac only', features: [bestsSFR.name] },
  { name: 'M_hist2: gasFrac+T', features: [bestsSFR.name, 'T (Hubble)'] },
  { name: 'M_hist3: gasFrac+T+eta', features: [bestsSFR.name, 'T (Hubble)', 'eta_rot'] },
  { name: 'M_full: all history', features: sSFRcorrs.map(c => c.name).concat(['T (Hubble)', 'eta_rot', 'gasExtent', 'compactness']) },
  { name: 'M_kitchen: everything', features: historyFeatures.map(f => f.name) },
  { name: 'M6: Per-galaxy', features: ['__free__'] },
];

function getVal(gIdx, fname) {
  const f = historyFeatures.find(x => x.name === fname);
  return f ? f.values[gIdx] : 0;
}

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

const bestHistGap = Math.max(...cvResults.filter(r => !r.name.includes('M0') && !r.name.includes('M6')).map(r => gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0));

log("  STEP 4: GALAXY QUARTILE ANALYSIS — HIGH vs LOW sSFR");
sep();
log("");

const sortedBysSFR = [...galaxyData].map((g, i) => ({ ...g, idx: i })).sort((a, b) => a[bestsSFR.name.replace('log(','').replace(')','').replace('/','_')] || a.logMHI_L - (b[bestsSFR.name.replace('log(','').replace(')','').replace('/','_')] || b.logMHI_L));
const sSFRvals = galaxyData.map(g => g.logMHI_L);
const sortedIdx = [...Array(N).keys()].sort((a, b) => sSFRvals[a] - sSFRvals[b]);

const q1Idx = sortedIdx.slice(0, Math.floor(N / 4));
const q4Idx = sortedIdx.slice(Math.floor(3 * N / 4));

const q1LogA0 = q1Idx.map(i => allLogA0[i]);
const q4LogA0 = q4Idx.map(i => allLogA0[i]);
const q1se = q1Idx.map(i => galaxyData[i].se);
const q4se = q4Idx.map(i => galaxyData[i].se);

const q1DL = dlMeta(q1LogA0, q1se);
const q4DL = dlMeta(q4LogA0, q4se);

const q1mean = q1LogA0.reduce((s,v)=>s+v,0)/q1LogA0.length;
const q4mean = q4LogA0.reduce((s,v)=>s+v,0)/q4LogA0.length;
const q1sd = Math.sqrt(q1LogA0.reduce((s,v)=>s+(v-q1mean)**2,0)/(q1LogA0.length-1));
const q4sd = Math.sqrt(q4LogA0.reduce((s,v)=>s+(v-q4mean)**2,0)/(q4LogA0.length-1));

log("  Split by " + bestsSFR.name + " quartiles:");
log("");
log("  Q1 (lowest sSFR — quiescent/old):");
log("    n=" + q1Idx.length + ", mean a0=" + Math.round(Math.pow(10, q1mean)));
log("    SD=" + q1sd.toFixed(3) + " dex, tau=" + q1DL.tau.toFixed(3) + " dex");
log("    Galaxies: " + q1Idx.map(i => galaxyData[i].name).join(", "));
log("");
log("  Q4 (highest sSFR — active/young):");
log("    n=" + q4Idx.length + ", mean a0=" + Math.round(Math.pow(10, q4mean)));
log("    SD=" + q4sd.toFixed(3) + " dex, tau=" + q4DL.tau.toFixed(3) + " dex");
log("    Galaxies: " + q4Idx.map(i => galaxyData[i].name).join(", "));
log("");

const splitDex = Math.abs(q1mean - q4mean);
const splitT = splitDex / Math.sqrt(q1sd**2/q1Idx.length + q4sd**2/q4Idx.length);

log("  Q1 vs Q4 split: " + splitDex.toFixed(3) + " dex (t=" + splitT.toFixed(2) + ")");
log("  Direction: " + (q4mean < q1mean ? "HIGH sSFR → LOWER a0" : "HIGH sSFR → HIGHER a0"));
log("");

log("  Within-quartile tau:");
log("    Q1 tau=" + q1DL.tau.toFixed(3) + " vs full tau=" + dlMeta(allLogA0, galaxyData.map(g=>g.se)).tau.toFixed(3));
log("    Q4 tau=" + q4DL.tau.toFixed(3));
const tauReductionQ = Math.min(q1DL.tau, q4DL.tau) / dlMeta(allLogA0, galaxyData.map(g=>g.se)).tau;
log("    Within-quartile tau reduction: " + ((1-tauReductionQ)*100).toFixed(1) + "%");
log("");

log("  STEP 5: PRINCIPAL COMPONENT OF GALAXY HISTORY");
sep();
log("");

const histVars = [
  { name: 'logMHI_L', values: galaxyData.map(g => g.logMHI_L) },
  { name: 'logGasFrac', values: galaxyData.map(g => g.logGasFrac) },
  { name: 'T', values: galaxyData.map(g => g.T) },
  { name: 'etaRot', values: galaxyData.map(g => g.etaRot) },
  { name: 'gasExtent', values: galaxyData.map(g => g.gasExtent) },
];

const means = histVars.map(v => v.values.reduce((s,x)=>s+x,0)/N);
const stds = histVars.map((v,j) => Math.sqrt(v.values.reduce((s,x)=>s+(x-means[j])**2,0)/N) || 1);
const Z = galaxyData.map((g, i) => histVars.map((v, j) => (v.values[i] - means[j]) / stds[j]));

const pc1 = galaxyData.map((g, i) => {
  return 0.5 * Z[i][0] + 0.5 * Z[i][1] + 0.3 * Z[i][2] + 0.2 * Z[i][3] + 0.2 * Z[i][4];
});

const pc1Corr = corrWith(pc1, deltaA0);
log("  Composite 'galaxy history index' (weighted MHI/L + gasFrac + T + eta + gasExt):");
log("  Correlation with delta_a0: r = " + pc1Corr.toFixed(3));
const pc1t = Math.abs(pc1Corr) * Math.sqrt(N-2) / Math.sqrt(1-pc1Corr**2+1e-10);
log("  |t| = " + pc1t.toFixed(1));
log("");

let looSS_pc1 = 0, looN_pc1 = 0;
for (let i = 0; i < N; i++) {
  const trainIdx = [...Array(N).keys()].filter(j => j !== i);
  const trainX = trainIdx.map(j => [[pc1[j]]]).map(x => x[0]);
  const trainY = trainIdx.map(j => allLogA0[j]);
  const reg = linReg(trainX, trainY);
  let predLogA0 = reg.intercept + reg.coefs[0] * pc1[i];
  predLogA0 = Math.max(2.5, Math.min(4.5, predLogA0));
  looSS_pc1 += predictSS(Math.pow(10, predLogA0), galaxyData[i].pts);
  looN_pc1 += galaxyData[i].pts.length;
}
const pc1GapClosed = gap > 0 ? (m0rms - Math.sqrt(looSS_pc1 / looN_pc1)) / gap * 100 : 0;
log("  LOO-CV gap closed by composite index: " + pc1GapClosed.toFixed(1) + "%");
log("");

log("=".repeat(80));
log("  PHASE 19 — CONCLUSIONS");
log("=".repeat(80));
log("");

log("  1. BEST sSFR PROXY: " + bestsSFR.name);
log("     r = " + bestsSFR.r.toFixed(3) + " with delta_a0 (t=" + bestsSFR.t.toFixed(1) + ")");
log("");
log("  2. DIRECTION: " + (bestsSFR.r < 0 ?
  "Gas-rich (high sSFR) galaxies → LOWER apparent a0" :
  "Gas-rich (high sSFR) galaxies → HIGHER apparent a0"));
log("");
log("  3. DOES sSFR ABSORB THE OTHER CLUES?");

const mhiF = historyFeatures.find(f => f.name === 'log(MHI)');
const rMHI_after = corrWith(mhiF.values, residAftersSFR);
const etaF = historyFeatures.find(f => f.name === 'eta_rot');
const rEta_after = corrWith(etaF.values, residAftersSFR);

if (Math.abs(rMHI_after) < 0.15 && Math.abs(rEta_after) < 0.15) {
  log("     YES — after removing sSFR effect, MHI and eta_rot are absorbed.");
  log("     sSFR IS the common underlying variable.");
} else {
  log("     PARTIALLY — MHI residual r=" + rMHI_after.toFixed(3) + ", eta_rot residual r=" + rEta_after.toFixed(3));
  log("     sSFR captures some but not all of the signal from MHI/eta.");
}
log("");

log("  4. GAP CLOSURE (LOO-CV):");
log("     sSFR alone:        " + (gap > 0 ? ((m0rms - cvResults[1].cvRMS) / gap * 100).toFixed(1) : '0') + "%");
log("     Best history model: " + bestHistGap.toFixed(1) + "%");
log("     Composite index:    " + pc1GapClosed.toFixed(1) + "%");
log("     Phase 15 baryon:   14-16%");
log("     M6 per-galaxy:     100%");
log("");

log("  5. QUARTILE SPLIT:");
log("     Q1 (quiescent): a0=" + Math.round(Math.pow(10, q1mean)) + ", tau=" + q1DL.tau.toFixed(3));
log("     Q4 (active):    a0=" + Math.round(Math.pow(10, q4mean)) + ", tau=" + q4DL.tau.toFixed(3));
log("     Split: " + splitDex.toFixed(3) + " dex (t=" + splitT.toFixed(2) + ")");
log("     Within-quartile tau reduction: " + ((1-tauReductionQ)*100).toFixed(1) + "%");
log("");

if (bestHistGap > 25) {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  sSFR / GALAXY HISTORY EXPLAINS SUBSTANTIAL a0 VARIATION.           ║");
  log("  ║  Star formation history IS the missing variable.                     ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else if (bestHistGap > 12) {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  sSFR / GALAXY HISTORY CONTRIBUTES to explaining a0 variation.      ║");
  log("  ║  It captures the MHI/T/eta clues but doesn't resolve everything.    ║");
  log("  ║  Better sSFR estimates (UV/IR) might close more of the gap.         ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  sSFR proxy does NOT substantially improve over individual clues.   ║");
  log("  ║  Either the proxy is too crude, or the variation has deeper roots.  ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  bestsSFRproxy: { name: bestsSFR.name, r: +bestsSFR.r.toFixed(3), t: +bestsSFR.t.toFixed(1) },
  correlations: corrs.map(c => ({ name: c.name, type: c.type, r: +c.r.toFixed(3), t: +c.t.toFixed(1) })),
  cvResults: cvResults.map(r => ({
    name: r.name, k: r.k, cvRMS: +r.cvRMS.toFixed(5),
    gapClosed: +(gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0).toFixed(1)
  })),
  bestHistGapClosed: +bestHistGap.toFixed(1),
  quartileSplit: {
    q1: { n: q1Idx.length, meanA0: Math.round(Math.pow(10, q1mean)), tau: +q1DL.tau.toFixed(3), sd: +q1sd.toFixed(3) },
    q4: { n: q4Idx.length, meanA0: Math.round(Math.pow(10, q4mean)), tau: +q4DL.tau.toFixed(3), sd: +q4sd.toFixed(3) },
    splitDex: +splitDex.toFixed(3), splitT: +splitT.toFixed(2)
  },
  compositeIndex: { r: +pc1Corr.toFixed(3), gapClosed: +pc1GapClosed.toFixed(1) }
};

fs.writeFileSync(path.join(__dirname, '../public/phase19-galaxy-history.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase19-galaxy-history.json");
