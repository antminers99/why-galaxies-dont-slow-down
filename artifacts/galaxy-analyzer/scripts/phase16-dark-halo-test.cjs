#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "16.0.0";
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

function hNFW(x) {
  if (x <= 0) return 0;
  return Math.log(1 + x) - x / (1 + x);
}

function fitNFW(pts) {
  let bestA = 1e4, bestRs = 5, bestSS = Infinity;

  for (let logA = 2; logA <= 7; logA += 0.25) {
    for (let logRs = -0.5; logRs <= 2.5; logRs += 0.15) {
      const A = Math.pow(10, logA);
      const rs = Math.pow(10, logRs);
      let ss = 0;
      for (const p of pts) {
        const r = p.r;
        const gobs = Math.pow(10, p.log_g_obs);
        const gbar = Math.pow(10, p.log_g_bar);
        const gDM_obs = Math.max(0, gobs - gbar);
        const x = r / rs;
        const gDM_nfw = A / (r * r) * hNFW(x);
        if (gDM_obs > 0 && gDM_nfw > 0) {
          ss += (Math.log10(gDM_obs) - Math.log10(gDM_nfw)) ** 2;
        } else {
          ss += (gDM_obs - gDM_nfw) ** 2 / (gobs * gobs + 1);
        }
      }
      if (ss < bestSS) { bestSS = ss; bestA = A; bestRs = rs; }
    }
  }

  let logA = Math.log10(bestA), logRs = Math.log10(bestRs);
  for (let iter = 0; iter < 80; iter++) {
    const step = 0.05 / (1 + iter * 0.1);
    let improved = false;
    for (const dA of [-step, 0, step]) {
      for (const dR of [-step, 0, step]) {
        if (dA === 0 && dR === 0) continue;
        const A = Math.pow(10, logA + dA);
        const rs = Math.pow(10, logRs + dR);
        let ss = 0;
        for (const p of pts) {
          const r = p.r;
          const gobs = Math.pow(10, p.log_g_obs);
          const gbar = Math.pow(10, p.log_g_bar);
          const gDM_obs = Math.max(0, gobs - gbar);
          const x = r / rs;
          const gDM_nfw = A / (r * r) * hNFW(x);
          if (gDM_obs > 0 && gDM_nfw > 0) {
            ss += (Math.log10(gDM_obs) - Math.log10(gDM_nfw)) ** 2;
          } else {
            ss += (gDM_obs - gDM_nfw) ** 2 / (gobs * gobs + 1);
          }
        }
        if (ss < bestSS) { bestSS = ss; logA += dA; logRs += dR; improved = true; }
      }
    }
    if (!improved) break;
  }

  const A = Math.pow(10, logA);
  const rs = Math.pow(10, logRs);

  const G_kpc = 4.302e-6;
  const rhoS = A / (4 * Math.PI * G_kpc * rs * rs * rs);

  let maxVDM2 = 0, rmax = 0;
  for (const p of pts) {
    const vdm2 = A / p.r * hNFW(p.r / rs);
    if (vdm2 > maxVDM2) { maxVDM2 = vdm2; rmax = p.r; }
  }
  const VmaxDM = Math.sqrt(maxVDM2);

  const rOut = pts[pts.length - 1].r;
  const rHalf = pts[Math.floor(pts.length / 2)].r;

  let fDM_half = 0, fDM_out = 0;
  {
    const gobs_h = Math.pow(10, pts[Math.floor(pts.length / 2)].log_g_obs);
    const gbar_h = Math.pow(10, pts[Math.floor(pts.length / 2)].log_g_bar);
    fDM_half = Math.max(0, 1 - gbar_h / gobs_h);
  }
  {
    const gobs_o = Math.pow(10, pts[pts.length - 1].log_g_obs);
    const gbar_o = Math.pow(10, pts[pts.length - 1].log_g_bar);
    fDM_out = Math.max(0, 1 - gbar_o / gobs_o);
  }

  const r200_approx = rs * Math.pow(200 * rhoS / (3 * 1.879e-29 * 3.086e19 * 3.086e19 * 3.086e19 / 1.989e30), 1 / 3);
  const c_approx = Math.max(1, Math.min(50, rOut / rs * 2));

  return {
    logA: +logA.toFixed(3),
    logRs: +logRs.toFixed(3),
    rs: +rs.toFixed(2),
    logRhoS: +Math.log10(rhoS > 0 ? rhoS : 1).toFixed(3),
    VmaxDM: +VmaxDM.toFixed(1),
    fDM_half: +fDM_half.toFixed(3),
    fDM_out: +fDM_out.toFixed(3),
    c_approx: +c_approx.toFixed(1),
    rms: +Math.sqrt(bestSS / pts.length).toFixed(4),
    n: pts.length
  };
}

const p11 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase11-sensitivity-lab.json'), 'utf8'));
const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));

const galaxyData = [];
for (const gal of p11.galaxies) {
  const pts = gal.localProfile.filter(p =>
    isFinite(p.log_g_bar) && isFinite(p.log_g_obs) && p.log_g_obs > p.log_g_bar * 0.5 && p.r > 0
  );
  if (pts.length < 5) continue;
  const sp = sparc.find(s => (s.name || s.Galaxy) === gal.name);

  let etaRot = 0;
  if (pts.length >= 4) {
    const oh = pts.slice(Math.floor(pts.length / 2)), ih = pts.slice(0, Math.floor(pts.length / 2));
    if (oh.length >= 2 && ih.length >= 2) {
      etaRot = (oh[oh.length-1].log_g_obs-oh[0].log_g_obs)/(oh[oh.length-1].log_g_bar-oh[0].log_g_bar+1e-10) -
               (ih[ih.length-1].log_g_obs-ih[0].log_g_obs)/(ih[ih.length-1].log_g_bar-ih[0].log_g_bar+1e-10);
    }
  }

  galaxyData.push({
    name: gal.name, pts, inc: gal.inc, D: gal.D, Vmax: gal.Vmax,
    T: gal.T || (sp ? sp.T : 5), etaRot, n: pts.length,
    logMHI: sp ? Math.log10(sp.MHI || 1e9) : 9,
    Rdisk: sp ? (sp.Rdisk || 3) : 3,
    SBdisk: sp ? (sp.SBdisk || 100) : 100,
    logL36: sp ? Math.log10(sp.L || 1e9) : 9,
  });
}

const N = galaxyData.length;

log("");
log("=".repeat(80));
log("  PHASE 16: DARK HALO EXPLANATORY TEST");
log("  Does halo structure explain the a0 variation?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  STEP 1: FIT NFW HALOS");
sep();
log("");

const haloFits = [];
const perGalA0 = [];
for (let i = 0; i < N; i++) {
  const g = galaxyData[i];
  const halo = fitNFW(g.pts);
  const a0fit = fitA0(g.pts);
  haloFits.push(halo);
  perGalA0.push(a0fit);
}

const meanLogRs = haloFits.reduce((s, h) => s + h.logRs, 0) / N;
const meanLogRhoS = haloFits.reduce((s, h) => s + h.logRhoS, 0) / N;
const meanFDM = haloFits.reduce((s, h) => s + h.fDM_out, 0) / N;
const meanVmaxDM = haloFits.reduce((s, h) => s + h.VmaxDM, 0) / N;

log("  NFW fit summary (" + N + " galaxies):");
log("    log(rs):      mean=" + meanLogRs.toFixed(2) + " kpc");
log("    log(rho_s):   mean=" + meanLogRhoS.toFixed(2) + " M_sun/kpc^3");
log("    f_DM(outer):  mean=" + meanFDM.toFixed(3));
log("    V_max,DM:     mean=" + meanVmaxDM.toFixed(1) + " km/s");
log("");

log("  STEP 2: CORRELATE HALO PROPERTIES WITH delta_a0");
sep();
log("");

const allLogA0 = perGalA0.map(f => f.logA0);
const meanA0 = allLogA0.reduce((s, v) => s + v, 0) / N;
const deltaA0 = allLogA0.map(v => v - meanA0);

const dmFeatures = [
  { name: 'log(rs)', values: haloFits.map(h => h.logRs) },
  { name: 'log(rho_s)', values: haloFits.map(h => h.logRhoS) },
  { name: 'f_DM(half)', values: haloFits.map(h => h.fDM_half) },
  { name: 'f_DM(out)', values: haloFits.map(h => h.fDM_out) },
  { name: 'V_max,DM', values: haloFits.map(h => h.VmaxDM) },
  { name: 'log(V_maxDM)', values: haloFits.map(h => Math.log10(Math.max(1, h.VmaxDM))) },
  { name: 'c_approx', values: haloFits.map(h => h.c_approx) },
  { name: 'log(A)', values: haloFits.map(h => h.logA) },
];

const baryonFeatures = [
  { name: 'log(MHI)', values: galaxyData.map(g => g.logMHI) },
  { name: 'log(D)', values: galaxyData.map(g => Math.log10(g.D)) },
  { name: 'inc', values: galaxyData.map(g => g.inc) },
  { name: 'T', values: galaxyData.map(g => g.T) },
  { name: 'Rdisk', values: galaxyData.map(g => g.Rdisk) },
  { name: 'log(Vmax)', values: galaxyData.map(g => Math.log10(g.Vmax)) },
  { name: 'eta_rot', values: galaxyData.map(g => g.etaRot) },
  { name: 'n_points', values: galaxyData.map(g => g.n) },
];

function corrWith(x, y) {
  const n = x.length;
  const mx = x.reduce((s, v) => s + v, 0) / n;
  const my = y.reduce((s, v) => s + v, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i]-mx)*(y[i]-my); sxx += (x[i]-mx)**2; syy += (y[i]-my)**2; }
  return sxy / Math.sqrt(sxx * syy + 1e-20);
}

log("  DM halo properties vs delta_a0:");
log("  ┌─────────────────────────────────────────────────┐");
for (const f of dmFeatures) {
  const r = corrWith(f.values, deltaA0);
  const bar = "#".repeat(Math.round(Math.abs(r) * 30));
  log("  │  " + f.name.padEnd(14) + " r=" + r.toFixed(3).padStart(7) + "  " + bar.padEnd(20) + "│");
}
log("  └─────────────────────────────────────────────────┘");
log("");
log("  Baryon properties vs delta_a0 (for comparison):");
log("  ┌─────────────────────────────────────────────────┐");
for (const f of baryonFeatures) {
  const r = corrWith(f.values, deltaA0);
  const bar = "#".repeat(Math.round(Math.abs(r) * 30));
  log("  │  " + f.name.padEnd(14) + " r=" + r.toFixed(3).padStart(7) + "  " + bar.padEnd(20) + "│");
}
log("  └─────────────────────────────────────────────────┘");
log("");

log("  STEP 3: MODEL LADDER WITH GALAXY-LEVEL LOO-CV");
sep();
log("");
log("  Train on N-1 galaxies, predict a0 for held-out galaxy.");
log("  Compare predicted a0 vs actual a0 on that galaxy's data points.");
log("");

function linReg(X, y) {
  const n = y.length, p = X[0] ? X[0].length : 0;
  if (p === 0) { const m = y.reduce((s,v)=>s+v,0)/n; return { coefs: [], intercept: m }; }
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

function getFeatureVec(gIdx, featureList) {
  return featureList.map(f => f.values[gIdx]);
}

const allFeatures = [...baryonFeatures, ...dmFeatures];

const modelDefs = [
  { name: 'M0: Universal a0', features: [] },
  { name: 'M1: Baryon only', features: baryonFeatures.map(f => f.name) },
  { name: 'M2: DM halo only', features: dmFeatures.map(f => f.name) },
  { name: 'M3: Baryon + DM', features: allFeatures.map(f => f.name) },
  { name: 'M6: Per-galaxy', features: ['__free__'] },
];

const cvResults = [];

for (const model of modelDefs) {
  let totalSS = 0, totalN = 0;

  if (model.features[0] === '__free__') {
    for (let i = 0; i < N; i++) {
      totalSS += predictSS(perGalA0[i].a0, galaxyData[i].pts);
      totalN += galaxyData[i].pts.length;
    }
    cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), totalSS, totalN, k: N });
    continue;
  }

  const featIdx = model.features.map(fn => allFeatures.findIndex(f => f.name === fn));

  for (let i = 0; i < N; i++) {
    const trainIdx = [...Array(N).keys()].filter(j => j !== i);
    const trainY = trainIdx.map(j => allLogA0[j]);

    if (model.features.length === 0) {
      const mu = trainY.reduce((s,v)=>s+v,0) / trainY.length;
      totalSS += predictSS(Math.pow(10, mu), galaxyData[i].pts);
    } else {
      const trainX = trainIdx.map(j => featIdx.map(fi => allFeatures[fi].values[j]));
      const reg = linReg(trainX, trainY);
      const testX = featIdx.map(fi => allFeatures[fi].values[i]);
      const predLogA0 = reg.intercept + reg.coefs.reduce((s,c,j) => s + c * testX[j], 0);
      const clampedLogA0 = Math.max(2.5, Math.min(4.5, predLogA0));
      totalSS += predictSS(Math.pow(10, clampedLogA0), galaxyData[i].pts);
    }
    totalN += galaxyData[i].pts.length;
  }
  cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), totalSS, totalN, k: model.features.length + 1 });
}

const m0rms = cvResults[0].cvRMS;
const m6rms = cvResults[cvResults.length - 1].cvRMS;
const gap = m0rms - m6rms;

log("  ┌──────────────────────────────────────────────────────────────────────────┐");
log("  │  Model                    k    CV-RMS    vs M0     vs M6    gap-closed  │");
log("  ├──────────────────────────────────────────────────────────────────────────┤");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const vsM6 = ((r.cvRMS / m6rms - 1) * 100).toFixed(1);
  const gapClosed = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  │  " + r.name.padEnd(25) + (r.k.toString()).padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(9) + ("+" + vsM6 + "%").padStart(9) +
    (gapClosed + "%").padStart(12) + "  │");
}
log("  └──────────────────────────────────────────────────────────────────────────┘");
log("");

const dmGapClosed = gap > 0 ? (m0rms - cvResults[2].cvRMS) / gap * 100 : 0;
const baryonGapClosed = gap > 0 ? (m0rms - cvResults[1].cvRMS) / gap * 100 : 0;
const combinedGapClosed = gap > 0 ? (m0rms - cvResults[3].cvRMS) / gap * 100 : 0;

log("  STEP 4: BEST DM-ONLY PREDICTORS (stepwise)");
sep();
log("");

const dmCorrs = dmFeatures.map(f => ({
  name: f.name,
  r: corrWith(f.values, deltaA0),
  values: f.values
})).sort((a,b) => Math.abs(b.r) - Math.abs(a.r));

log("  Top DM predictors of delta_a0:");
for (const c of dmCorrs) {
  const sig = Math.abs(c.r) * Math.sqrt(N - 2) / Math.sqrt(1 - c.r * c.r + 1e-10);
  log("    " + c.name.padEnd(14) + " r=" + c.r.toFixed(3).padStart(7) +
    " |t|=" + sig.toFixed(1).padStart(5) + (sig > 2 ? " **" : sig > 1.65 ? " *" : ""));
}
log("");

log("  STEP 5: RIDGE-REGULARIZED DM MODEL (prevent overfitting)");
sep();
log("");

function ridgeReg(X, y, lambda) {
  const n = y.length, p = X[0].length;
  const means = [];
  for (let j = 0; j < p; j++) { let s = 0; for (let i = 0; i < n; i++) s += X[i][j]; means.push(s / n); }
  const my = y.reduce((s,v) => s + v, 0) / n;
  const stds = [];
  for (let j = 0; j < p; j++) {
    let s = 0; for (let i = 0; i < n; i++) s += (X[i][j] - means[j]) ** 2;
    stds.push(Math.sqrt(s / n) || 1);
  }
  const Xs = X.map(row => row.map((v, j) => (v - means[j]) / stds[j]));
  const ys = y.map(v => v - my);
  const XtX = Array.from({length:p}, () => Array(p).fill(0));
  const Xty = Array(p).fill(0);
  for (let i = 0; i < n; i++) { for (let j = 0; j < p; j++) { Xty[j] += Xs[i][j] * ys[i]; for (let k = 0; k < p; k++) XtX[j][k] += Xs[i][j] * Xs[i][k]; } }
  for (let j = 0; j < p; j++) XtX[j][j] += lambda;
  const A = XtX.map((r,i) => [...r, Xty[i]]);
  for (let j = 0; j < p; j++) { let mx = j; for (let i=j+1;i<p;i++) if(Math.abs(A[i][j])>Math.abs(A[mx][j]))mx=i; [A[j],A[mx]]=[A[mx],A[j]]; for(let i=j+1;i<p;i++){const f=A[i][j]/A[j][j];for(let k=j;k<=p;k++)A[i][k]-=f*A[j][k];} }
  const bStd = Array(p).fill(0);
  for (let j = p-1; j >= 0; j--) { bStd[j] = A[j][p]; for (let k=j+1;k<p;k++) bStd[j] -= A[j][k]*bStd[k]; bStd[j] /= A[j][j]; }
  const coefs = bStd.map((b, j) => b / stds[j]);
  const intercept = my - coefs.reduce((s, c, j) => s + c * means[j], 0);
  return { coefs, intercept };
}

const bestLambdas = [0.1, 0.5, 1, 2, 5, 10, 20, 50];
let bestLambda = 1, bestRidgeRMS = Infinity;

for (const lambda of bestLambdas) {
  let totalSS = 0, totalN = 0;
  for (let i = 0; i < N; i++) {
    const trainIdx = [...Array(N).keys()].filter(j => j !== i);
    const trainX = trainIdx.map(j => dmFeatures.map(f => f.values[j]));
    const trainY = trainIdx.map(j => allLogA0[j]);
    const reg = ridgeReg(trainX, trainY, lambda);
    const testX = dmFeatures.map(f => f.values[i]);
    const predLogA0 = reg.intercept + reg.coefs.reduce((s, c, j) => s + c * testX[j], 0);
    const clamped = Math.max(2.5, Math.min(4.5, predLogA0));
    totalSS += predictSS(Math.pow(10, clamped), galaxyData[i].pts);
    totalN += galaxyData[i].pts.length;
  }
  const rms = Math.sqrt(totalSS / totalN);
  if (rms < bestRidgeRMS) { bestRidgeRMS = rms; bestLambda = lambda; }
}

log("  Best ridge lambda: " + bestLambda);
log("  Ridge DM-only CV-RMS: " + bestRidgeRMS.toFixed(5));
const ridgeGapClosed = gap > 0 ? (m0rms - bestRidgeRMS) / gap * 100 : 0;
log("  Gap closed: " + ridgeGapClosed.toFixed(1) + "%");
log("");

let combinedRidgeRMS = Infinity;
for (const lambda of bestLambdas) {
  let totalSS = 0, totalN = 0;
  const allF = [...baryonFeatures, ...dmFeatures];
  for (let i = 0; i < N; i++) {
    const trainIdx = [...Array(N).keys()].filter(j => j !== i);
    const trainX = trainIdx.map(j => allF.map(f => f.values[j]));
    const trainY = trainIdx.map(j => allLogA0[j]);
    const reg = ridgeReg(trainX, trainY, lambda);
    const testX = allF.map(f => f.values[i]);
    const predLogA0 = reg.intercept + reg.coefs.reduce((s, c, j) => s + c * testX[j], 0);
    const clamped = Math.max(2.5, Math.min(4.5, predLogA0));
    totalSS += predictSS(Math.pow(10, clamped), galaxyData[i].pts);
    totalN += galaxyData[i].pts.length;
  }
  const rms = Math.sqrt(totalSS / totalN);
  if (rms < combinedRidgeRMS) combinedRidgeRMS = rms;
}

const combinedRidgeGap = gap > 0 ? (m0rms - combinedRidgeRMS) / gap * 100 : 0;
log("  Ridge Baryon+DM CV-RMS: " + combinedRidgeRMS.toFixed(5));
log("  Gap closed (baryon+DM, ridge): " + combinedRidgeGap.toFixed(1) + "%");
log("");

log("=".repeat(80));
log("  PHASE 16 — DECISION");
log("=".repeat(80));
log("");
log("  HIERARCHY OF MODELS (gap closed by each):");
log("");

const allModelsForViz = [
  { name: 'M0: Universal a0', gap: 0 },
  { name: 'M1: Baryon only', gap: baryonGapClosed },
  { name: 'M2: DM halo only', gap: dmGapClosed },
  { name: 'M2r: DM ridge', gap: ridgeGapClosed },
  { name: 'M3: Baryon+DM', gap: combinedGapClosed },
  { name: 'M3r: Baryon+DM ridge', gap: combinedRidgeGap },
  { name: 'M6: Per-galaxy', gap: 100 },
];

for (const m of allModelsForViz) {
  const barLen = Math.max(0, Math.round(m.gap / 2));
  const bar = "\u2588".repeat(barLen);
  log("  " + m.name.padEnd(24) + (m.gap.toFixed(1) + "%").padStart(7) + "  " + bar);
}
log("");
log("  " + "\u2500".repeat(40));
log("  0%         25%        50%        75%        100%");
log("  M0                                           M6");
log("");

const bestExplanatory = Math.max(dmGapClosed, ridgeGapClosed, combinedRidgeGap);

if (bestExplanatory > 60) {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  DM HALO STRUCTURE EXPLAINS MOST OF THE a0 VARIATION.               ║");
  log("  ║  The 'variation' in a0 is largely a reflection of halo diversity.    ║");
  log("  ║  a0 may be effectively universal after halo corrections.             ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else if (bestExplanatory > 30) {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  DM HALO STRUCTURE PARTIALLY EXPLAINS a0 VARIATION.                 ║");
  log("  ║  Halo properties contribute, but a large fraction remains.           ║");
  log("  ║  Neither baryons alone nor halos alone resolve the puzzle.           ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  DM HALO PROPERTIES DO NOT EXPLAIN a0 VARIATION.                    ║");
  log("  ║  Even with NFW halo fits, most variation remains unexplained.        ║");
  log("  ║  The variation is either deeper physics or unmeasured variables.     ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  DETAILED INTERPRETATION:");
log("");
log("  DM-only (M2) gap closed:        " + dmGapClosed.toFixed(1) + "%");
log("  DM ridge gap closed:            " + ridgeGapClosed.toFixed(1) + "%");
log("  Baryon-only (M1) gap closed:    " + baryonGapClosed.toFixed(1) + "%");
log("  Baryon+DM combined gap closed:  " + combinedRidgeGap.toFixed(1) + "%");
log("  Remaining unexplained:          " + (100 - combinedRidgeGap).toFixed(1) + "%");
log("");

if (dmGapClosed > baryonGapClosed * 1.5) {
  log("  → DM halo is a BETTER predictor of a0 variation than baryonic properties.");
  log("  → This supports: apparent a0 variation reflects halo structure diversity.");
} else if (baryonGapClosed > dmGapClosed * 1.5) {
  log("  → Baryonic properties are BETTER predictors than halo structure.");
  log("  → DM halo fitting does not add substantial explanatory power.");
} else {
  log("  → DM and baryon properties have SIMILAR explanatory power.");
  log("  → Both contribute partially to explaining a0 variation.");
}
log("");

log("  CORRELATION BETWEEN DM AND BARYON FEATURES:");
const dmBaryonCorr = corrWith(
  haloFits.map(h => h.fDM_out),
  galaxyData.map(g => g.logMHI)
);
log("  f_DM(out) vs log(MHI): r = " + dmBaryonCorr.toFixed(3));
const dmVmaxCorr = corrWith(
  haloFits.map(h => Math.log10(Math.max(1, h.VmaxDM))),
  galaxyData.map(g => Math.log10(g.Vmax))
);
log("  log(V_max,DM) vs log(Vmax): r = " + dmVmaxCorr.toFixed(3));
log("");

log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  description: "Phase 16: Dark Halo Explanatory Test",
  nGalaxies: N,
  haloFitSummary: {
    meanLogRs: +meanLogRs.toFixed(3),
    meanLogRhoS: +meanLogRhoS.toFixed(3),
    meanFDM_out: +meanFDM.toFixed(3),
    meanVmaxDM: +meanVmaxDM.toFixed(1)
  },
  dmCorrelations: dmFeatures.map(f => ({
    name: f.name,
    r: +corrWith(f.values, deltaA0).toFixed(3)
  })),
  cvResults: cvResults.map(r => ({
    name: r.name, k: r.k,
    cvRMS: +r.cvRMS.toFixed(5),
    gapClosed: +(gap > 0 ? (m0rms - r.cvRMS) / gap * 100 : 0).toFixed(1)
  })),
  ridgeDM: { cvRMS: +bestRidgeRMS.toFixed(5), gapClosed: +ridgeGapClosed.toFixed(1), lambda: bestLambda },
  ridgeCombined: { cvRMS: +combinedRidgeRMS.toFixed(5), gapClosed: +combinedRidgeGap.toFixed(1) },
  verdict: bestExplanatory > 60 ? "DM_EXPLAINS" : bestExplanatory > 30 ? "PARTIAL" : "DM_INSUFFICIENT",
  perGalaxy: galaxyData.map((g, i) => ({
    name: g.name,
    logA0: +allLogA0[i].toFixed(4),
    deltaA0: +deltaA0[i].toFixed(4),
    halo: haloFits[i]
  }))
};

fs.writeFileSync(path.join(__dirname, '../public/phase16-dark-halo.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase16-dark-halo.json");
