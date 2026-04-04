#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "21.0.0";
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

  let etaRot = 0;
  if (pts.length >= 4) {
    const oh = pts.slice(Math.floor(pts.length / 2)), ih = pts.slice(0, Math.floor(pts.length / 2));
    if (oh.length >= 2 && ih.length >= 2) {
      etaRot = (oh[oh.length-1].log_g_obs-oh[0].log_g_obs)/(oh[oh.length-1].log_g_bar-oh[0].log_g_bar+1e-10) -
               (ih[ih.length-1].log_g_obs-ih[0].log_g_obs)/(ih[ih.length-1].log_g_bar-ih[0].log_g_bar+1e-10);
    }
  }

  const logMHI_L = Math.log10(sp.MHI / L);
  const logSBdisk = sp.SBdisk > 0 ? Math.log10(sp.SBdisk) : 2;
  const logSBeff = sp.SBeff > 0 ? Math.log10(sp.SBeff) : 2;
  const concentration = sp.Reff > 0 && sp.Rdisk > 0 ? sp.Reff / sp.Rdisk : 1;
  const logConcentration = Math.log10(concentration);
  const gasExtent = sp.RHI > 0 && sp.Rdisk > 0 ? sp.RHI / sp.Rdisk : 3;

  const MZR_proxy = 8.9 + 0.3 * (logMstar - 10) - 0.1 * (logMstar - 10) ** 2;
  const FMR_proxy = MZR_proxy - 0.32 * logMHI_L;
  const yieldProxy = logMstar - Math.log10(sp.MHI * 1e9) * 0.6;
  const enrichmentProxy = logMstar - logMHI_L * 0.3;

  const metalGrad_proxy = -0.1 * (logMstar - 9.5) + 0.05 * sp.T;
  const effectiveYield = logMstar - Math.log10(sp.MHI * 1e9 + L * 0.5e9) + 9;
  const closedBoxDeviation = logMHI_L + 0.5 * logMstar - 4.5;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q, fD: sp.fD,
    n: pts.length, etaRot,
    logMHI: Math.log10(sp.MHI), logL36: Math.log10(L),
    logMstar, logMHI_L, logSBdisk, logSBeff,
    logConcentration, gasExtent,
    Vflat: sp.Vflat,
    MZR_proxy, FMR_proxy, yieldProxy, enrichmentProxy,
    metalGrad_proxy, effectiveYield, closedBoxDeviation,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

log("");
log("=".repeat(80));
log("  PHASE 21: METALLICITY — CHEMICAL ENRICHMENT HISTORY");
log("  Does metallicity explain a0 variation?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  METALLICITY STATUS IN SPARC:");
sep();
log("");
log("  Direct metallicity measurements: NOT in SPARC catalog");
log("  Gas-phase O/H: NOT available");
log("  Stellar [Fe/H]: NOT available");
log("");
log("  PROXIES CONSTRUCTED:");
log("    1. MZR proxy:   12+log(O/H) ~ 8.9 + 0.3*(logM*-10) - 0.1*(logM*-10)^2");
log("       The mass-metallicity relation (Tremonti+04, Zahid+14)");
log("    2. FMR proxy:   MZR - 0.32*log(MHI/L)");
log("       Fundamental metallicity relation (Mannucci+10): Z(M*, SFR)");
log("    3. Yield proxy: logM* - 0.6*log(MHI)");
log("       Effective nucleosynthetic yield indicator");
log("    4. Enrichment:  logM* - 0.3*log(MHI/L)");
log("       Overall enrichment level");
log("    5. Metal gradient: -0.1*(logM*-9.5) + 0.05*T");
log("       Radial metallicity gradient proxy");
log("    6. Effective yield: logM* - log(MHI+M*) + 9");
log("       Closed-box-like effective yield");
log("    7. Closed-box deviation: logMHI/L + 0.5*logM* - 4.5");
log("       Deviation from simple closed-box enrichment");
log("");

log("  IMPORTANT CAVEAT:");
log("  All proxies are constructed from Mstar and MHI — the same");
log("  variables already tested. If metallicity proxies correlate,");
log("  it may simply reflect the known Mstar/MHI correlations.");
log("  We need to check for INDEPENDENT signal beyond Mstar and MHI.");
log("");

log("  STEP 1: CORRELATIONS — METALLICITY PROXIES vs delta_a0");
sep();
log("");

const allFeatures = [
  { name: 'log(MHI)', values: galaxyData.map(g => g.logMHI), type: 'gas mass' },
  { name: 'log(Mstar)', values: galaxyData.map(g => g.logMstar), type: 'stellar mass' },
  { name: 'MZR proxy', values: galaxyData.map(g => g.MZR_proxy), type: 'metallicity' },
  { name: 'FMR proxy', values: galaxyData.map(g => g.FMR_proxy), type: 'metallicity' },
  { name: 'yield proxy', values: galaxyData.map(g => g.yieldProxy), type: 'metallicity' },
  { name: 'enrichment', values: galaxyData.map(g => g.enrichmentProxy), type: 'metallicity' },
  { name: 'metal gradient', values: galaxyData.map(g => g.metalGrad_proxy), type: 'metallicity' },
  { name: 'eff. yield', values: galaxyData.map(g => g.effectiveYield), type: 'metallicity' },
  { name: 'closed-box dev', values: galaxyData.map(g => g.closedBoxDeviation), type: 'metallicity' },
  { name: 'eta_rot', values: galaxyData.map(g => g.etaRot), type: 'RC shape' },
  { name: 'T (Hubble)', values: galaxyData.map(g => g.T), type: 'morphology' },
  { name: 'log(MHI/L)', values: galaxyData.map(g => g.logMHI_L), type: 'gas fraction' },
  { name: 'n_points', values: galaxyData.map(g => g.n), type: 'data quality' },
  { name: 'log(D)', values: galaxyData.map(g => Math.log10(g.D)), type: 'distance' },
  { name: 'inc', values: galaxyData.map(g => g.inc), type: 'geometry' },
];

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
  const marker = c.type === 'metallicity' ? " [Z]" : "";
  log("  │  " + c.name.padEnd(19) + c.type.padEnd(18) +
    c.r.toFixed(3).padStart(7) + c.t.toFixed(1).padStart(8) + "    " + (sig + marker).padEnd(8) + "│");
}
log("  └───────────────────────────────────────────────────────────────────────┘");
log("");

const metalCorrs = corrs.filter(c => c.type === 'metallicity');
const bestMetal = metalCorrs[0];
log("  BEST METALLICITY PROXY: " + bestMetal.name + " (r=" + bestMetal.r.toFixed(3) + ", t=" + bestMetal.t.toFixed(1) + ")");
log("");

log("  STEP 2: INDEPENDENCE TEST — DOES METALLICITY ADD TO MHI + MSTAR?");
sep();
log("");

const mhiVals = galaxyData.map(g => g.logMHI);
const mstarVals = galaxyData.map(g => g.logMstar);

for (const mc of metalCorrs) {
  const rZ = corrWith(mc.values, deltaA0);
  const rZ_MHI = corrWith(mc.values, mhiVals);
  const rZ_Mstar = corrWith(mc.values, mstarVals);

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
  const rZ_afterMHI = corrWith(mc.values, residAfterMHI);
  const tZ_afterMHI = Math.abs(rZ_afterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rZ_afterMHI**2+1e-10);

  const residAfterBoth = [];
  {
    const X = galaxyData.map(g => [g.logMHI, g.logMstar]);
    const y = deltaA0;
    const reg = linReg(X, y);
    for(let i=0;i<N;i++) residAfterBoth.push(y[i] - (reg.intercept + reg.coefs[0]*X[i][0] + reg.coefs[1]*X[i][1]));
  }
  const rZ_afterBoth = corrWith(mc.values, residAfterBoth);
  const tZ_afterBoth = Math.abs(rZ_afterBoth) * Math.sqrt(N-4) / Math.sqrt(1-rZ_afterBoth**2+1e-10);

  const absorbed = Math.abs(rZ_afterBoth) < 0.10;
  log("  " + mc.name.padEnd(16) +
    " raw r=" + rZ.toFixed(3) +
    " | after MHI: r=" + rZ_afterMHI.toFixed(3) + " (t=" + tZ_afterMHI.toFixed(1) + ")" +
    " | after MHI+M*: r=" + rZ_afterBoth.toFixed(3) + " (t=" + tZ_afterBoth.toFixed(1) + ")" +
    (absorbed ? " ABSORBED" : " INDEPENDENT"));
}
log("");

log("  STEP 3: PERMUTATION TEST — IS ANY METALLICITY SIGNAL REAL?");
sep();
log("");

const NPERMS = 5000;
for (const mc of metalCorrs.slice(0, 3)) {
  const obsR = Math.abs(corrWith(mc.values, deltaA0));
  let count = 0;
  const delta = [...deltaA0];
  for (let p = 0; p < NPERMS; p++) {
    for (let i = delta.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [delta[i], delta[j]] = [delta[j], delta[i]];
    }
    if (Math.abs(corrWith(mc.values, delta)) >= obsR) count++;
  }
  const perm_p = count / NPERMS;
  log("  " + mc.name.padEnd(16) + " |r|=" + obsR.toFixed(3) + " perm_p=" + perm_p.toFixed(4) +
    (perm_p < 0.05 ? " *" : "") + (perm_p < 0.01 ? "*" : ""));
}
log("");

log("  STEP 4: LOO CROSS-VALIDATION — METALLICITY MODELS");
sep();
log("");

function getVal(gIdx, fname) {
  const f = allFeatures.find(x => x.name === fname);
  return f ? f.values[gIdx] : 0;
}

const modelDefs = [
  { name: 'M0: Universal a0', features: [] },
  { name: 'M_Z1: best Z proxy', features: [bestMetal.name] },
  { name: 'M_Z2: MZR+FMR', features: ['MZR proxy', 'FMR proxy'] },
  { name: 'M_Z3: all Z proxies', features: metalCorrs.map(c => c.name) },
  { name: 'M_Z+MHI: Z+gas', features: [bestMetal.name, 'log(MHI)'] },
  { name: 'M_MHI: gas only', features: ['log(MHI)'] },
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

const bestZGap = Math.max(...cvResults.filter(r => r.name.includes('M_Z')).map(r => gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0));
const mhiGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_MHI')).cvRMS) / gap * 100 : 0;

log("  STEP 5: QUARTILE ANALYSIS — HIGH vs LOW METALLICITY");
sep();
log("");

const bestMetalVals = bestMetal.values;
const sortedIdx = [...Array(N).keys()].sort((a, b) => bestMetalVals[a] - bestMetalVals[b]);
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

log("  Split by " + bestMetal.name + " (best metallicity proxy):");
log("");
log("  Q1 (lowest Z — metal-poor):");
log("    n=" + q1Idx.length + ", mean a0=" + Math.round(Math.pow(10, q1mean)));
log("    SD=" + q1sd.toFixed(3) + " dex, tau=" + q1DL.tau.toFixed(3) + " dex");
log("    Galaxies: " + q1Idx.map(i => galaxyData[i].name).join(", "));
log("");
log("  Q4 (highest Z — metal-rich):");
log("    n=" + q4Idx.length + ", mean a0=" + Math.round(Math.pow(10, q4mean)));
log("    SD=" + q4sd.toFixed(3) + " dex, tau=" + q4DL.tau.toFixed(3) + " dex");
log("    Galaxies: " + q4Idx.map(i => galaxyData[i].name).join(", "));
log("");
log("  Split: " + splitDex.toFixed(3) + " dex (t=" + splitT.toFixed(2) + ")");
log("  Direction: " + (q4mean < q1mean ? "HIGH Z → LOWER a0" : "HIGH Z → HIGHER a0"));
log("");
log("  Within-quartile tau:");
log("    Q1 (metal-poor): tau=" + q1DL.tau.toFixed(3) + " (" + ((1-q1DL.tau/fullDL.tau)*100).toFixed(0) + "% reduction)");
log("    Q4 (metal-rich): tau=" + q4DL.tau.toFixed(3) + " (" + ((1-q4DL.tau/fullDL.tau)*100).toFixed(0) + "% reduction)");
log("    Full:            tau=" + fullDL.tau.toFixed(3));
log("");

log("  STEP 6: PHYSICS — METALLICITY CHANNELS TO a0");
sep();
log("");
log("  How could metallicity physically affect measured a0?");
log("");
log("  Channel 1: M*/L at 3.6um");
log("    Metal-rich stars → slightly higher M*/L at 3.6um");
log("    But effect is SMALL at 3.6um (dominated by AGB/RGB stars)");
log("    Expected delta(M*/L) ~ 0.05-0.1 per dex in Z");
log("    → Would produce delta(log a0) ~ 0.03-0.05 per dex in Z");
log("    → TOO SMALL to explain tau=0.22");
log("");
log("  Channel 2: Gas cooling and disk structure");
log("    Higher Z → more efficient cooling → denser gas → different dynamics");
log("    But this affects gas distribution, not directly a0");
log("");
log("  Channel 3: Correlated with mass (MZR)");
log("    Massive galaxies are metal-rich AND have different dynamics");
log("    This is the DOMINANT channel but it's already captured by Mstar");
log("");

log("  STEP 7: STRICT SUCCESS CRITERIA EVALUATION");
sep();
log("");

const criteria = [
  { name: '|t| >= 2.0', pass: bestMetal.t >= 2.0, val: bestMetal.t.toFixed(1) },
  { name: 'perm p < 0.05', pass: false, val: 'see above' },
  { name: 'LOO > M0', pass: bestZGap > 0, val: bestZGap.toFixed(1) + '%' },
  { name: 'gap > age (0.9%)', pass: bestZGap > 0.9, val: bestZGap.toFixed(1) + '%' },
  { name: 'Independent of MHI', pass: false, val: 'see step 2' },
];

const obsR_best = Math.abs(corrWith(bestMetal.values, deltaA0));
let permCount = 0;
const shuffled = [...deltaA0];
for (let p = 0; p < NPERMS; p++) {
  for (let i = shuffled.length - 1; i > 0; i--) {
    const j = Math.floor(Math.random() * (i + 1));
    [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]];
  }
  if (Math.abs(corrWith(bestMetal.values, shuffled)) >= obsR_best) permCount++;
}
const bestPerm = permCount / NPERMS;
criteria[1].pass = bestPerm < 0.05;
criteria[1].val = 'p=' + bestPerm.toFixed(4);

const residAfterMHI_best = [];
{
  const x = mhiVals, y = deltaA0;
  const n = x.length;
  const mx = x.reduce((s,v)=>s+v,0)/n, my = y.reduce((s,v)=>s+v,0)/n;
  let sxy=0,sxx=0;
  for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;}
  const b = sxy/sxx, a = my - b*mx;
  for(let i=0;i<n;i++) residAfterMHI_best.push(y[i] - (a + b*x[i]));
}
const rZ_afterMHI_best = corrWith(bestMetal.values, residAfterMHI_best);
const tZ_afterMHI_best = Math.abs(rZ_afterMHI_best) * Math.sqrt(N-3) / Math.sqrt(1-rZ_afterMHI_best**2+1e-10);
criteria[4].pass = tZ_afterMHI_best >= 1.65;
criteria[4].val = 'r=' + rZ_afterMHI_best.toFixed(3) + ' t=' + tZ_afterMHI_best.toFixed(1);

const nPassed = criteria.filter(c => c.pass).length;

log("  ┌─────────────────────────────────────────────────────────────────┐");
log("  │  Criterion              Pass?    Value                        │");
log("  ├─────────────────────────────────────────────────────────────────┤");
for (const c of criteria) {
  log("  │  " + c.name.padEnd(22) + (c.pass ? "  YES " : "  NO  ") + "  " + c.val.padEnd(30) + "│");
}
log("  ├─────────────────────────────────────────────────────────────────┤");
log("  │  TOTAL: " + nPassed + "/5 criteria met" + " ".repeat(37) + "│");
log("  └─────────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  PHASE 21 — CONCLUSIONS");
log("=".repeat(80));
log("");

log("  1. BEST METALLICITY PROXY: " + bestMetal.name);
log("     r = " + bestMetal.r.toFixed(3) + " with delta_a0 (t=" + bestMetal.t.toFixed(1) + ")");
log("");

log("  2. INDEPENDENCE: " + (tZ_afterMHI_best >= 1.65 ? "YES" : "NO"));
log("     After MHI: r=" + rZ_afterMHI_best.toFixed(3) + " (t=" + tZ_afterMHI_best.toFixed(1) + ")");
log("");

log("  3. LOO-CV GAP CLOSURE:");
log("     Best Z model: " + bestZGap.toFixed(1) + "%");
log("     MHI alone:    " + mhiGap.toFixed(1) + "%");
log("     Phase 20 age: 0.9%");
log("     Phase 19 sSFR: -1.8%");
log("     Phase 15 baryon: 14-16%");
log("");

log("  4. QUARTILE SPLIT:");
log("     Metal-poor: a0=" + Math.round(Math.pow(10, q1mean)) + ", tau=" + q1DL.tau.toFixed(3));
log("     Metal-rich: a0=" + Math.round(Math.pow(10, q4mean)) + ", tau=" + q4DL.tau.toFixed(3));
log("     Split: " + splitDex.toFixed(3) + " dex (t=" + splitT.toFixed(2) + ")");
log("");

log("  5. SUCCESS CRITERIA: " + nPassed + "/5 → " + (nPassed >= 4 ? "CONFIRMED" : nPassed >= 2 ? "PARTIAL" : "FAILED"));
log("");

if (nPassed >= 4) {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  METALLICITY IS A CONFIRMED PREDICTOR of a0 variation.             ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else if (nPassed >= 2) {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  Metallicity shows PARTIAL evidence but does not meet the strict   ║");
  log("  ║  success criteria. Like age, it tracks mass-related properties     ║");
  log("  ║  but adds no independent explanatory power beyond MHI/Mstar.      ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  Metallicity proxies FAIL. No independent signal beyond MHI/Mstar. ║");
  log("  ║  All apparent correlations are driven by the mass-metallicity      ║");
  log("  ║  relation — metallicity is a REDUNDANT label for stellar mass.     ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  GALAXY HISTORY SCORECARD (cumulative):");
log("  ┌────────────────────────────────────────────────────────────────────┐");
log("  │  Variable       Status        Gap Closed  Independent?           │");
log("  ├────────────────────────────────────────────────────────────────────┤");
log("  │  sSFR           CHECKED       -1.8%       partially              │");
log("  │  Stellar age    CHECKED        0.9%       NO (absorbed by MHI)   │");
log("  │  Metallicity    CHECKED       " + bestZGap.toFixed(1).padStart(5) + "%       " + (tZ_afterMHI_best >= 1.65 ? "YES" : "NO (absorbed by MHI)").padEnd(23) + "│");
log("  │  SFH            PENDING       —           —                      │");
log("  │  Merger proxy   PENDING       —           —                      │");
log("  └────────────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  bestMetalProxy: { name: bestMetal.name, r: +bestMetal.r.toFixed(3), t: +bestMetal.t.toFixed(1) },
  correlations: corrs.map(c => ({ name: c.name, type: c.type, r: +c.r.toFixed(3), t: +c.t.toFixed(1) })),
  cvResults: cvResults.map(r => ({
    name: r.name, k: r.k, cvRMS: +r.cvRMS.toFixed(5),
    gapClosed: +(gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0).toFixed(1)
  })),
  bestZGapClosed: +bestZGap.toFixed(1),
  independence: {
    afterMHI: { r: +rZ_afterMHI_best.toFixed(3), t: +tZ_afterMHI_best.toFixed(1) },
  },
  permutationP: +bestPerm.toFixed(4),
  successCriteria: { passed: nPassed, total: 5 },
  quartileSplit: {
    q1: { n: q1Idx.length, meanA0: Math.round(Math.pow(10, q1mean)), tau: +q1DL.tau.toFixed(3) },
    q4: { n: q4Idx.length, meanA0: Math.round(Math.pow(10, q4mean)), tau: +q4DL.tau.toFixed(3) },
    splitDex: +splitDex.toFixed(3), splitT: +splitT.toFixed(2)
  },
  comparison: {
    phase19_sSFR: -1.8,
    phase20_age: 0.9,
    phase21_metal: +bestZGap.toFixed(1),
    phase15_baryon: 15.7
  }
};

fs.writeFileSync(path.join(__dirname, '../public/phase21-metallicity.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase21-metallicity.json");
