#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "20.0.0";
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

  const concentration = sp.Reff > 0 && sp.Rdisk > 0 ? sp.Reff / sp.Rdisk : 1;
  const logConcentration = Math.log10(concentration);
  const bulgeProxy = sp.SBeff > 0 && sp.SBdisk > 0 ? sp.SBeff / sp.SBdisk : 0.1;
  const logBulgeProxy = Math.log10(bulgeProxy > 0 ? bulgeProxy : 0.01);

  const earlyType = sp.T <= 3 ? 1 : 0;
  const logSBdisk = sp.SBdisk > 0 ? Math.log10(sp.SBdisk) : 2;
  const logSBeff = sp.SBeff > 0 ? Math.log10(sp.SBeff) : 2;

  const logMHI_L = Math.log10(sp.MHI / L);
  const gasRichness = Math.log10(sp.MHI / L);
  const gasPoor = logMHI_L < -1.0 ? 1 : 0;

  const ageProxy1 = -sp.T * 0.4 + logMstar * 0.3 + logConcentration * 0.3;
  const ageProxy2 = -gasRichness * 0.5 + logSBeff * 0.3 + logConcentration * 0.2;
  const downsizingAge = logMstar;
  const morphAge = -sp.T;
  const concAge = logConcentration;
  const sbAge = logSBeff;

  const massWeightedAgeProxy = 0.4 * (logMstar - 9) / 2 + 0.3 * (-sp.T + 5) / 10 + 0.2 * logConcentration + 0.1 * (-logMHI_L + 1) / 3;
  const compositeAgeIndex = massWeightedAgeProxy;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q, fD: sp.fD,
    n: pts.length, etaRot,
    logMHI: Math.log10(sp.MHI), logL36: Math.log10(L),
    logMstar, logMHI_L,
    concentration, logConcentration,
    bulgeProxy, logBulgeProxy,
    earlyType, logSBdisk, logSBeff,
    gasPoor, gasRichness,
    ageProxy1, ageProxy2,
    downsizingAge, morphAge, concAge, sbAge,
    compositeAgeIndex,
    Vflat: sp.Vflat, SBdisk: sp.SBdisk, SBeff: sp.SBeff,
    Reff: sp.Reff, Rdisk: sp.Rdisk, RHI: sp.RHI, MHI: sp.MHI, L36: L,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

log("");
log("=".repeat(80));
log("  PHASE 20: STELLAR AGE — THE ACCUMULATED HISTORY");
log("  Does the age of the stellar population explain a0 variation?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  WHY STELLAR AGE AFTER sSFR:");
sep();
log("");
log("  Phase 19 showed: sSFR (current star formation) is NOT the answer.");
log("  But sSFR measures only the PRESENT. Stellar age measures the");
log("  ACCUMULATED HISTORY — how old the bulk of stars actually are.");
log("");
log("  Key physics pathway:");
log("    old population → higher M*/L at 3.6um → different mass model →");
log("    different gbar → different fitted a0");
log("");
log("  At 3.6um, M*/L varies from ~0.3 (young) to ~0.8 (old).");
log("  We assumed 0.5 universally. If age-dependent M*/L is the");
log("  unmodeled systematic, it would manifest as a0 variation.");
log("");

log("  STELLAR AGE PROXIES FROM SPARC:");
sep();
log("");
log("  Direct age measurements: NOT available in SPARC");
log("  Proxies constructed from available data:");
log("    1. log(Mstar)        downsizing: massive = older");
log("    2. -T (morphology)   early-type = older");
log("    3. log(Reff/Rdisk)   concentration: bulge-dominated = older");
log("    4. log(SBeff)        effective SB: high = older/denser");
log("    5. -log(MHI/L)       gas-poor = older (inverted gas fraction)");
log("    6. Composite age     weighted combination of 1-5");
log("");

const ageFeatures = [
  { name: 'log(Mstar)', values: galaxyData.map(g => g.logMstar), type: 'downsizing' },
  { name: '-T (morph)', values: galaxyData.map(g => -g.T), type: 'morphology' },
  { name: 'log(Reff/Rdisk)', values: galaxyData.map(g => g.logConcentration), type: 'concentration' },
  { name: 'log(SBeff)', values: galaxyData.map(g => g.logSBeff), type: 'surface brightness' },
  { name: '-log(MHI/L)', values: galaxyData.map(g => -g.logMHI_L), type: 'gas depletion' },
  { name: 'compositeAge', values: galaxyData.map(g => g.compositeAgeIndex), type: 'composite' },
  { name: 'log(SBdisk)', values: galaxyData.map(g => g.logSBdisk), type: 'disk SB' },
  { name: 'log(bulge)', values: galaxyData.map(g => g.logBulgeProxy), type: 'bulge fraction' },
  { name: 'eta_rot', values: galaxyData.map(g => g.etaRot), type: 'RC shape' },
  { name: 'log(MHI)', values: galaxyData.map(g => g.logMHI), type: 'gas mass' },
  { name: 'log(D)', values: galaxyData.map(g => Math.log10(g.D)), type: 'distance' },
  { name: 'inc', values: galaxyData.map(g => g.inc), type: 'geometry' },
  { name: 'n_points', values: galaxyData.map(g => g.n), type: 'data quality' },
];

log("  STEP 1: CORRELATIONS — AGE PROXIES vs delta_a0");
sep();
log("");

const corrs = ageFeatures.map(f => {
  const r = corrWith(f.values, deltaA0);
  const t = Math.abs(r) * Math.sqrt(N - 2) / Math.sqrt(1 - r * r + 1e-10);
  return { ...f, r, t };
}).sort((a, b) => Math.abs(b.r) - Math.abs(a.r));

log("  ┌───────────────────────────────────────────────────────────────────────┐");
log("  │  Feature            Type              r        |t|    sig           │");
log("  ├───────────────────────────────────────────────────────────────────────┤");
for (const c of corrs) {
  const sig = c.t > 2.58 ? "***" : c.t > 2.0 ? "**" : c.t > 1.65 ? "*" : "";
  log("  │  " + c.name.padEnd(19) + c.type.padEnd(18) +
    c.r.toFixed(3).padStart(7) + c.t.toFixed(1).padStart(8) + "    " + sig.padEnd(4) + "│");
}
log("  └───────────────────────────────────────────────────────────────────────┘");
log("");

const ageProxies = corrs.filter(c => ['downsizing','morphology','concentration','surface brightness','gas depletion','composite'].includes(c.type));
const bestAge = ageProxies[0];
log("  BEST AGE PROXY: " + bestAge.name + " (r=" + bestAge.r.toFixed(3) + ", t=" + bestAge.t.toFixed(1) + ")");
log("");

log("  STEP 2: MUTUAL CORRELATIONS AMONG AGE PROXIES");
sep();
log("");

const ageNames = ['log(Mstar)', '-T (morph)', 'log(Reff/Rdisk)', 'log(SBeff)', '-log(MHI/L)'];
log("  " + "".padEnd(16) + ageNames.map(n => n.substring(0,8).padStart(9)).join(""));
for (const v1 of ageNames) {
  const f1 = ageFeatures.find(f => f.name === v1);
  const row = ageNames.map(v2 => {
    if (v1 === v2) return "    1.00";
    const f2 = ageFeatures.find(f => f.name === v2);
    return corrWith(f1.values, f2.values).toFixed(2).padStart(8);
  });
  log("  " + v1.padEnd(16) + row.join(""));
}
log("");

log("  STEP 3: THE M*/L PATHWAY — AGE → M*/L → a0?");
sep();
log("");

log("  Theory: If stellar age varies, then M*/L at 3.6um varies.");
log("  We assumed M*/L = 0.5 everywhere. True range: 0.3-0.8.");
log("  If actual M*/L differs from 0.5, gbar shifts, changing fitted a0.");
log("");
log("  Prediction: if M*/L is the mechanism, then");
log("    old (high M*/L) → gbar overestimated → a0 shifts DOWN");
log("    young (low M*/L) → gbar underestimated → a0 shifts UP");
log("");

const mstarCorr = corrWith(galaxyData.map(g => g.logMstar), deltaA0);
const morphCorr = corrWith(galaxyData.map(g => -g.T), deltaA0);
log("  Observed:");
log("    log(Mstar) vs delta_a0: r = " + mstarCorr.toFixed(3));
log("    Direction: massive (old) → " + (mstarCorr < 0 ? "LOWER" : "HIGHER") + " a0");
log("    -T (early=old) vs delta_a0: r = " + morphCorr.toFixed(3));
log("    Direction: early-type (old) → " + (morphCorr < 0 ? "LOWER" : "HIGHER") + " a0");
log("");

const mlPrediction = mstarCorr < 0;
if (mlPrediction) {
  log("  ✓ CONSISTENT with M*/L pathway: old/massive → lower a0");
  log("    This is what we'd expect if M*/L > 0.5 for old galaxies");
  log("    → gbar too high → a0 pushed down to compensate");
} else {
  log("  ✗ INCONSISTENT with M*/L pathway: old/massive → higher a0");
  log("    This goes AGAINST the M*/L prediction");
}
log("");

log("  Quantifying the M*/L effect:");
log("  If delta(log M*/L) ~ " + (mstarCorr * 0.2).toFixed(2) + " per dex in Mstar,");
log("  this would produce delta(log a0) ~ " + (mstarCorr * 0.2 * 0.5).toFixed(3) + " per dex in Mstar");
log("  Observed slope: " + (corrs.find(c => c.name === 'log(Mstar)').r * Math.sqrt(deltaA0.reduce((s,v)=>s+v*v,0)/N) / Math.sqrt(galaxyData.map(g=>g.logMstar).reduce((s,v)=>s+(v-galaxyData.map(g=>g.logMstar).reduce((a,b)=>a+b,0)/N)**2,0)/N)).toFixed(3) + " dex per dex Mstar");
log("");

log("  STEP 4: LOO CROSS-VALIDATION — AGE-BASED MODELS");
sep();
log("");

function getVal(gIdx, fname) {
  const f = ageFeatures.find(x => x.name === fname);
  return f ? f.values[gIdx] : 0;
}

const modelDefs = [
  { name: 'M0: Universal a0', features: [] },
  { name: 'M_age1: best age proxy', features: [bestAge.name] },
  { name: 'M_age2: Mstar+morph', features: ['log(Mstar)', '-T (morph)'] },
  { name: 'M_age3: Mstar+morph+conc', features: ['log(Mstar)', '-T (morph)', 'log(Reff/Rdisk)'] },
  { name: 'M_age4: all age proxies', features: ageNames },
  { name: 'M_age+eta: age+RC shape', features: [...ageNames, 'eta_rot'] },
  { name: 'M_kitchen: everything', features: ageFeatures.map(f => f.name) },
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

const bestAgeGap = Math.max(...cvResults.filter(r => !r.name.includes('M0') && !r.name.includes('M6') && !r.name.includes('kitchen')).map(r => gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0));
const kitchenGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('kitchen')).cvRMS) / gap * 100 : 0;

log("  STEP 5: QUARTILE ANALYSIS — OLD vs YOUNG GALAXIES");
sep();
log("");

const bestAgeVals = bestAge.values;
const sortedIdx = [...Array(N).keys()].sort((a, b) => bestAgeVals[a] - bestAgeVals[b]);
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

log("  Split by " + bestAge.name + " (best age proxy):");
log("");
log("  Q1 (youngest population):");
log("    n=" + q1Idx.length + ", mean a0=" + Math.round(Math.pow(10, q1mean)) + ", SD=" + q1sd.toFixed(3) + ", tau=" + q1DL.tau.toFixed(3));
log("    Galaxies: " + q1Idx.map(i => galaxyData[i].name).join(", "));
log("");
log("  Q4 (oldest population):");
log("    n=" + q4Idx.length + ", mean a0=" + Math.round(Math.pow(10, q4mean)) + ", SD=" + q4sd.toFixed(3) + ", tau=" + q4DL.tau.toFixed(3));
log("    Galaxies: " + q4Idx.map(i => galaxyData[i].name).join(", "));
log("");
log("  Split: " + splitDex.toFixed(3) + " dex (t=" + splitT.toFixed(2) + ")");
log("  Direction: " + (q4mean > q1mean ? "OLD → HIGHER a0" : "OLD → LOWER a0"));
log("");
log("  Within-quartile tau:");
log("    Q1 (young): tau=" + q1DL.tau.toFixed(3) + " (" + ((1-q1DL.tau/fullDL.tau)*100).toFixed(0) + "% reduction)");
log("    Q4 (old):   tau=" + q4DL.tau.toFixed(3) + " (" + ((1-q4DL.tau/fullDL.tau)*100).toFixed(0) + "% reduction)");
log("    Full:       tau=" + fullDL.tau.toFixed(3));
log("");

log("  STEP 6: PARTIAL CORRELATION — AGE vs MHI");
sep();
log("");

const mhiVals = galaxyData.map(g => g.logMHI);
const rMHI = corrWith(mhiVals, deltaA0);
const rAge = corrWith(bestAgeVals, deltaA0);
const rMHI_Age = corrWith(mhiVals, bestAgeVals);

const rMHI_partial = (rMHI - rAge * rMHI_Age) / Math.sqrt((1 - rAge**2) * (1 - rMHI_Age**2) + 1e-10);
const rAge_partial = (rAge - rMHI * rMHI_Age) / Math.sqrt((1 - rMHI**2) * (1 - rMHI_Age**2) + 1e-10);
const tMHI_partial = Math.abs(rMHI_partial) * Math.sqrt(N - 3) / Math.sqrt(1 - rMHI_partial**2 + 1e-10);
const tAge_partial = Math.abs(rAge_partial) * Math.sqrt(N - 3) / Math.sqrt(1 - rAge_partial**2 + 1e-10);

log("  Correlation matrix:");
log("    delta_a0 vs MHI: r = " + rMHI.toFixed(3));
log("    delta_a0 vs age: r = " + rAge.toFixed(3));
log("    MHI vs age:      r = " + rMHI_Age.toFixed(3));
log("");
log("  Partial correlations (controlling for the other):");
log("    MHI | age:  r_partial = " + rMHI_partial.toFixed(3) + " (t=" + tMHI_partial.toFixed(1) + ")");
log("    age | MHI:  r_partial = " + rAge_partial.toFixed(3) + " (t=" + tAge_partial.toFixed(1) + ")");
log("");

if (Math.abs(rMHI_partial) > Math.abs(rAge_partial)) {
  log("  → MHI remains DOMINANT after controlling for age");
  if (tAge_partial < 1.65) {
    log("  → Age has NO independent signal beyond MHI");
  } else {
    log("  → Age has some independent signal alongside MHI");
  }
} else {
  log("  → Age is DOMINANT after controlling for MHI");
}
log("");

log("  STEP 7: SENSITIVITY — WHAT IF M*/L ACTUALLY VARIES?");
sep();
log("");

log("  Thought experiment: What M*/L variation is needed to explain tau?");
log("  tau = " + fullDL.tau.toFixed(3) + " dex means delta(log a0) ~ +-" + fullDL.tau.toFixed(3));
log("  In RAR: a0 ~ gbar / f(gbar/a0)");
log("  If gbar = (M*/L) * L * G/r^2, then delta(log gbar) = delta(log M*/L)");
log("  Roughly: delta(log a0) ~ delta(log M*/L) (in deep MOND regime)");
log("  So we need delta(log M*/L) ~ " + fullDL.tau.toFixed(3) + " dex = factor " + Math.pow(10, fullDL.tau).toFixed(2));
log("  M*/L range needed: 0.5 / " + Math.pow(10, fullDL.tau).toFixed(2) + " to 0.5 * " + Math.pow(10, fullDL.tau).toFixed(2));
log("  = " + (0.5 / Math.pow(10, fullDL.tau)).toFixed(2) + " to " + (0.5 * Math.pow(10, fullDL.tau)).toFixed(2));
log("  Known 3.6um M*/L range: ~0.2 to ~0.8");
log("  Required range: " + (0.5/Math.pow(10, fullDL.tau)).toFixed(2) + " to " + (0.5*Math.pow(10, fullDL.tau)).toFixed(2));
const mlRangeOK = (0.5/Math.pow(10, fullDL.tau)) >= 0.15 && (0.5*Math.pow(10, fullDL.tau)) <= 1.0;
log("  " + (mlRangeOK ? "PLAUSIBLE" : "IMPLAUSIBLE") + " — the required M*/L range " +
  (mlRangeOK ? "falls within known bounds" : "exceeds known bounds"));
log("");

log("  STEP 8: COMPARISON WITH PHASE 19 (sSFR)");
sep();
log("");
log("  ┌────────────────────────────────────────────────────────────────────┐");
log("  │                    Phase 19 (sSFR)    Phase 20 (age)             │");
log("  ├────────────────────────────────────────────────────────────────────┤");
log("  │  Best proxy r:      0.215             " + bestAge.r.toFixed(3).padEnd(22) + "│");
log("  │  Best proxy t:      1.6               " + bestAge.t.toFixed(1).padEnd(22) + "│");
log("  │  LOO gap closed:    -1.8%             " + bestAgeGap.toFixed(1).padEnd(21) + "%│");
log("  │  Kitchen gap:       10.4%             " + kitchenGap.toFixed(1).padEnd(21) + "%│");
log("  │  Q split (dex):     0.030             " + splitDex.toFixed(3).padEnd(22) + "│");
log("  │  Q1 tau:            0.161             " + q1DL.tau.toFixed(3).padEnd(22) + "│");
log("  │  Q4 tau:            0.257             " + q4DL.tau.toFixed(3).padEnd(22) + "│");
log("  └────────────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  PHASE 20 — CONCLUSIONS");
log("=".repeat(80));
log("");

log("  1. BEST AGE PROXY: " + bestAge.name);
log("     r = " + bestAge.r.toFixed(3) + " with delta_a0 (t=" + bestAge.t.toFixed(1) + ")");
log("");

const ageDirection = bestAge.r < 0 ? "older" : "younger";
log("  2. DIRECTION: " + ageDirection + " stellar populations → " + (bestAge.r < 0 ? "LOWER" : "HIGHER") + " a0");
log("");

log("  3. PARTIAL CORRELATIONS:");
log("     MHI after controlling for age: r=" + rMHI_partial.toFixed(3) + " (t=" + tMHI_partial.toFixed(1) + ")");
log("     Age after controlling for MHI:  r=" + rAge_partial.toFixed(3) + " (t=" + tAge_partial.toFixed(1) + ")");
log("     → " + (Math.abs(rMHI_partial) > Math.abs(rAge_partial) ? "MHI is DOMINANT" : "Age is DOMINANT"));
log("");

log("  4. M*/L PATHWAY: " + (mlRangeOK ? "PLAUSIBLE" : "IMPLAUSIBLE") + " as full explanation");
log("     Required M*/L range: " + (0.5/Math.pow(10, fullDL.tau)).toFixed(2) + "-" + (0.5*Math.pow(10, fullDL.tau)).toFixed(2));
log("     Known range: 0.2-0.8");
log("");

log("  5. GAP CLOSURE:");
log("     Best age model (LOO): " + bestAgeGap.toFixed(1) + "%");
log("     Phase 19 sSFR: -1.8%");
log("     Phase 15 baryon: 14-16%");
log("");

if (bestAgeGap > 20) {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  STELLAR AGE PROXIES EXPLAIN SUBSTANTIAL a0 VARIATION.             ║");
  log("  ║  The accumulated history matters more than current activity.        ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else if (bestAgeGap > 5) {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  Age proxies show MODEST improvement over sSFR.                    ║");
  log("  ║  Some signal exists but does not resolve the heterogeneity.        ║");
  log("  ║  M*/L variation remains a plausible unmodeled systematic.          ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  VERDICT:");
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  Age proxies (from SPARC) do NOT explain the heterogeneity.        ║");
  log("  ║  Like sSFR, they fail out-of-sample prediction.                    ║");
  log("  ║  The SPARC proxies may be too crude to capture real age effects.   ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}

log("");

log("  GALAXY HISTORY SCORECARD (cumulative):");
log("  ┌────────────────────────────────────────────────┐");
log("  │  Variable       Status        Gap Closed       │");
log("  ├────────────────────────────────────────────────┤");
log("  │  sSFR           CHECKED       -1.8% (fail)     │");
log("  │  Stellar age    CHECKED       " + bestAgeGap.toFixed(1).padStart(5) + "%" + "          │");
log("  │  Metallicity    PENDING       —                │");
log("  │  SFH            PENDING       —                │");
log("  │  Merger proxy   PENDING       —                │");
log("  └────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  bestAgeProxy: { name: bestAge.name, r: +bestAge.r.toFixed(3), t: +bestAge.t.toFixed(1) },
  correlations: corrs.map(c => ({ name: c.name, type: c.type, r: +c.r.toFixed(3), t: +c.t.toFixed(1) })),
  cvResults: cvResults.map(r => ({
    name: r.name, k: r.k, cvRMS: +r.cvRMS.toFixed(5),
    gapClosed: +(gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0).toFixed(1)
  })),
  bestAgeGapClosed: +bestAgeGap.toFixed(1),
  partialCorrelations: {
    mhi_given_age: { r: +rMHI_partial.toFixed(3), t: +tMHI_partial.toFixed(1) },
    age_given_mhi: { r: +rAge_partial.toFixed(3), t: +tAge_partial.toFixed(1) }
  },
  mlPathway: {
    plausible: mlRangeOK,
    requiredRange: [(0.5/Math.pow(10, fullDL.tau)).toFixed(2), (0.5*Math.pow(10, fullDL.tau)).toFixed(2)].map(Number)
  },
  quartileSplit: {
    q1: { n: q1Idx.length, meanA0: Math.round(Math.pow(10, q1mean)), tau: +q1DL.tau.toFixed(3), sd: +q1sd.toFixed(3) },
    q4: { n: q4Idx.length, meanA0: Math.round(Math.pow(10, q4mean)), tau: +q4DL.tau.toFixed(3), sd: +q4sd.toFixed(3) },
    splitDex: +splitDex.toFixed(3), splitT: +splitT.toFixed(2)
  },
  comparison: {
    phase19_sSFR_gapClosed: -1.8,
    phase20_age_gapClosed: +bestAgeGap.toFixed(1),
    phase15_baryon_gapClosed: 15.7
  }
};

fs.writeFileSync(path.join(__dirname, '../public/phase20-stellar-age.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase20-stellar-age.json");
