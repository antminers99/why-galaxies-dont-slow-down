#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 133C: WHAT DRIVES VfResid? — Coupling Law Derivation');
console.log('');
console.log('  Instead of predicting a₀, now predict VfResid itself.');
console.log('  What halo/environment/dynamics properties generate the');
console.log('  coupling residual? Especially in the massive galaxy regime.');
console.log('======================================================================\n');

const stageA = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'stage-A-master-table.json'), 'utf8'));
const tidalData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase58a2-tidal-expansion.json'), 'utf8')).tidalData;
const sparcTable = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-table.json'), 'utf8'));
const sparcResults = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'sparc-results.json'), 'utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparcTable.forEach(s => { sparcMap[s.name] = s; });
const resMap = {}; sparcResults.perGalaxy.forEach(g => { if (g.rawName) resMap[g.rawName] = g; });
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
const N = gals45.length;

function mean(a) { return a.reduce((s, v) => s + v, 0) / a.length; }
function sd(a) { const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function pearsonR(x, y) {
  const n = x.length, mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}
function solveLinear(A, b) {
  const n = A.length;
  const M = A.map((r, i) => [...r, b[i]]);
  for (let i = 0; i < n; i++) {
    let mx = i;
    for (let j = i + 1; j < n; j++) if (Math.abs(M[j][i]) > Math.abs(M[mx][i])) mx = j;
    [M[i], M[mx]] = [M[mx], M[i]];
    if (Math.abs(M[i][i]) < 1e-15) continue;
    for (let j = i + 1; j < n; j++) {
      const f = M[j][i] / M[i][i];
      for (let k = i; k <= n; k++) M[j][k] -= f * M[i][k];
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = M[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= M[i][j] * x[j];
    x[i] /= M[i][i];
  }
  return x;
}
function ols(Y, X) {
  const n = Y.length, p = X[0].length + 1;
  const Xa = X.map(r => [1, ...r]);
  const XtX = Array.from({ length: p }, () => new Array(p).fill(0));
  const XtY = new Array(p).fill(0);
  for (let i = 0; i < n; i++)
    for (let j = 0; j < p; j++) {
      XtY[j] += Xa[i][j] * Y[i];
      for (let l = 0; l < p; l++) XtX[j][l] += Xa[i][j] * Xa[i][l];
    }
  const beta = solveLinear(XtX, XtY);
  const resid = Y.map((y, i) => y - Xa[i].reduce((s, x, j) => s + x * beta[j], 0));
  const rss = resid.reduce((s, r) => s + r * r, 0);
  const tss = Y.reduce((s, y) => s + (y - mean(Y)) ** 2, 0);
  return { beta, resid, rss, tss, r2: tss > 0 ? 1 - rss / tss : 0 };
}
function looR2(Y, X) {
  const n = Y.length; let ss = 0;
  const tss = Y.reduce((s, y) => s + (y - mean(Y)) ** 2, 0);
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xt = [...X.slice(0, i), ...X.slice(i + 1)];
    try {
      const f = ols(Yt, Xt);
      const xi = [1, ...X[i]];
      ss += (Y[i] - xi.reduce((s, x, j) => s + x * f.beta[j], 0)) ** 2;
    } catch (e) { ss += (Y[i] - mean(Y)) ** 2; }
  }
  return 1 - ss / tss;
}

const logVflat = gals45.map(g => Math.log10(sparcMap[g.name].Vflat));
const Vflat_raw = gals45.map(g => sparcMap[g.name].Vflat);
const logMbar = gals45.map(g => {
  const s = sparcMap[g.name]; return Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
});
const logL36 = gals45.map(g => Math.log10(Math.max(sparcMap[g.name].L36, 0.01)));
const logRdisk = gals45.map(g => Math.log10(sparcMap[g.name].Rdisk));
const morphT = gals45.map(g => sparcMap[g.name].T);
const structX = gals45.map((_, i) => [logMbar[i], logL36[i], logRdisk[i], morphT[i]]);
const VfResid = ols(logVflat, structX).resid;

function safeLog(v) { return v > 0 ? Math.log10(v) : Math.log10(Math.max(Math.abs(v), 1)); }

const predictors = {};
gals45.forEach((g, i) => {
  const r = resMap[g.name]; const s = sparcMap[g.name];
  const dhl = r.models.dark_halo_linear; const lh = r.models.log_halo;
  const mo = r.models.mond; const newt = r.models.newtonian;
  const mgh = r.models.modified_gravity_halo;
  const dhf = r.models.dark_halo_flat;

  const Vflat2 = s.Vflat * s.Vflat;
  const Mstar = s.L36 * 0.5 * 1e9;
  const Mgas = Math.pow(10, g.logMHI) * 1.33 * 1e9;
  const Vbar2 = 6.674e-11 * (Mstar + Mgas) * 1.989e30 / (s.Rdisk * 3.086e19) / 1e6;
  const fDM = Math.max(0, 1 - Math.min(1, Vbar2 / Vflat2));

  if (i === 0) {
    Object.keys({
      haloK_linear: 0, dhl_mse: 0, dhl_improve: 0,
      lh_improve: 0, lh_outerImprove: 0, lh_innerImprove: 0,
      haloK_log: 0, mond_improve: 0, mond_mse_ratio: 0,
      mgh_k: 0, dhf_improve: 0, newt_mse: 0,
      fDM_flat: 0, logExtent: 0, logSBeff: 0,
      concentration: 0, envCode: 0, logMhost: 0,
      logMeanRun: 0, logMHI: 0, logSig0: 0
    }).forEach(k => { predictors[k] = []; });
  }

  predictors.haloK_linear.push(safeLog(Math.max(dhl.k, 1)));
  predictors.dhl_mse.push(Math.log10(Math.max(dhl.mse, 0.1)));
  predictors.dhl_improve.push(dhl.improvementVsNewton);
  predictors.lh_improve.push(lh.improvementVsNewton);
  predictors.lh_outerImprove.push(lh.outerImprovement);
  predictors.lh_innerImprove.push(lh.innerImprovement);
  predictors.haloK_log.push(safeLog(Math.max(lh.k, 1)));
  predictors.mond_improve.push(mo.improvementVsNewton);
  predictors.mond_mse_ratio.push(Math.log10(Math.max(mo.mse / Math.max(newt.mse, 0.1), 0.001)));
  predictors.mgh_k.push(safeLog(Math.max(mgh.k, 1)));
  predictors.dhf_improve.push(dhf.improvementVsNewton);
  predictors.newt_mse.push(Math.log10(Math.max(newt.mse, 0.1)));
  predictors.fDM_flat.push(fDM);
  predictors.logExtent.push(Math.log10(r.maxR / s.Rdisk));
  predictors.logSBeff.push(Math.log10(Math.max(s.SBeff, 0.01)));
  predictors.concentration.push(Math.log10(s.Reff / s.Rdisk));
  predictors.envCode.push(g.envCode);
  predictors.logMhost.push(tdMap[g.name].logMhost);
  predictors.logMeanRun.push(g.logMeanRun);
  predictors.logMHI.push(g.logMHI);
  predictors.logSig0.push(g.logSigma0);
});

console.log('  VfResid: mean=' + mean(VfResid).toFixed(4) + ' sd=' + sd(VfResid).toFixed(4) + '\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 1: WHAT CORRELATES WITH VfResid? — Full sample (N=45)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const predNames = Object.keys(predictors);
const corrs = predNames.map(name => {
  const r = pearsonR(VfResid, predictors[name]);
  return { name, r, rAbs: Math.abs(r) };
}).sort((a, b) => b.rAbs - a.rAbs);

for (const c of corrs) {
  const flag = c.rAbs > 0.5 ? '★★★' : c.rAbs > 0.3 ? '★★ ' : c.rAbs > 0.15 ? '★  ' : '   ';
  console.log('  ' + flag + ' r=' + c.r.toFixed(3).padStart(7) + '  ' + c.name);
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 2: PREDICTING VfResid — Single-variable models');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const topPreds = corrs.filter(c => c.rAbs > 0.15);
const singleResults = [];
console.log('  ' + 'Predictor'.padEnd(22) + '  R²(fit)  R²(LOO)  r(VfR)');
console.log('  ' + '─'.repeat(55));

for (const tp of topPreds) {
  const X = predictors[tp.name].map(v => [v]);
  const fit = ols(VfResid, X);
  const looR2val = looR2(VfResid, X);
  singleResults.push({ name: tp.name, r2Fit: +fit.r2.toFixed(3), r2LOO: +looR2val.toFixed(3), r: +tp.r.toFixed(3) });
  console.log('  ' + tp.name.padEnd(22) + fit.r2.toFixed(3).padStart(7) + '  ' + looR2val.toFixed(3).padStart(7) + '  ' + tp.r.toFixed(3));
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 3: MULTI-PREDICTOR MODELS FOR VfResid');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const modelSpecs = [
  { name: 'haloK only', vars: ['haloK_linear'] },
  { name: 'haloK + mondImprove', vars: ['haloK_linear', 'mond_improve'] },
  { name: 'haloK + mondMSEratio', vars: ['haloK_linear', 'mond_mse_ratio'] },
  { name: 'haloK + env + MR', vars: ['haloK_linear', 'envCode', 'logMeanRun'] },
  { name: 'haloK + mghK', vars: ['haloK_linear', 'mgh_k'] },
  { name: 'haloK + lhOuter', vars: ['haloK_linear', 'lh_outerImprove'] },
  { name: 'haloK + mondImprove + env', vars: ['haloK_linear', 'mond_improve', 'envCode'] },
  { name: 'haloK + mondImprove + lhOuter', vars: ['haloK_linear', 'mond_improve', 'lh_outerImprove'] },
  { name: 'ALL top 5', vars: ['haloK_linear', 'mond_improve', 'mgh_k', 'envCode', 'logMeanRun'] }
];

console.log('  ' + 'Model'.padEnd(35) + '  R²(fit)  R²(LOO)  p');
console.log('  ' + '─'.repeat(60));

const multiResults = [];
for (const spec of modelSpecs) {
  const X = gals45.map((_, i) => spec.vars.map(v => predictors[v][i]));
  const fit = ols(VfResid, X);
  const looR2val = looR2(VfResid, X);
  multiResults.push({ name: spec.name, vars: spec.vars, r2Fit: +fit.r2.toFixed(3), r2LOO: +looR2val.toFixed(3), p: spec.vars.length });
  console.log('  ' + spec.name.padEnd(35) + fit.r2.toFixed(3).padStart(7) + '  ' + looR2val.toFixed(3).padStart(7) + '  ' + spec.vars.length);
}
console.log('');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 4: REGIME-SPECIFIC — High-Vflat galaxies only');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const highIdx = gals45.map((_, i) => i).filter(i => Vflat_raw[i] >= 120);
const VfResid_high = highIdx.map(i => VfResid[i]);
const nHigh = highIdx.length;

console.log('  High-Vflat regime (Vflat >= 120): N=' + nHigh + '\n');

if (nHigh >= 12) {
  const highCorrs = predNames.map(name => {
    const r = pearsonR(VfResid_high, highIdx.map(i => predictors[name][i]));
    return { name, r, rAbs: Math.abs(r) };
  }).sort((a, b) => b.rAbs - a.rAbs);

  console.log('  Top correlates in high-Vflat regime:');
  for (const c of highCorrs.slice(0, 10)) {
    const flag = c.rAbs > 0.5 ? '★★★' : c.rAbs > 0.3 ? '★★ ' : '★  ';
    console.log('  ' + flag + ' r=' + c.r.toFixed(3).padStart(7) + '  ' + c.name);
  }
  console.log('');

  console.log('  Multi-predictor models (high-Vflat):');
  console.log('  ' + 'Model'.padEnd(35) + '  R²(fit)  R²(LOO)');
  console.log('  ' + '─'.repeat(55));

  for (const spec of modelSpecs.slice(0, 6)) {
    const X = highIdx.map(i => spec.vars.map(v => predictors[v][i]));
    if (X.length < spec.vars.length + 3) continue;
    try {
      const fit = ols(VfResid_high, X);
      const looR2val = looR2(VfResid_high, X);
      console.log('  ' + spec.name.padEnd(35) + fit.r2.toFixed(3).padStart(7) + '  ' + looR2val.toFixed(3).padStart(7));
    } catch (e) { }
  }
  console.log('');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROUND 5: WHAT REMAINS UNEXPLAINED?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const bestModel = multiResults.sort((a, b) => b.r2LOO - a.r2LOO)[0];
console.log('  Best model for VfResid: ' + bestModel.name + ' (LOO R²=' + bestModel.r2LOO + ')\n');

const bestX = gals45.map((_, i) => bestModel.vars.map(v => predictors[v][i]));
const bestFit = ols(VfResid, bestX);
const unexplained = bestFit.resid;
const unexpSD = sd(unexplained);
const origSD = sd(VfResid);

console.log('  Original VfResid SD: ' + origSD.toFixed(4));
console.log('  Unexplained residual SD: ' + unexpSD.toFixed(4));
console.log('  Fraction explained: ' + (bestFit.r2 * 100).toFixed(1) + '%');
console.log('  Fraction unexplained: ' + ((1 - bestFit.r2) * 100).toFixed(1) + '%\n');

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const core3X = gals45.map((_, i) => [gals45[i].logMHI, tdMap[gals45[i].name].logMhost, gals45[i].logMeanRun]);
const coreUnexplainedX = gals45.map((_, j) => [...core3X[j], unexplained[j]]);

function gapPct(rms, sdy) { return sdy > 0 ? 100 * (1 - rms ** 2 / sdy ** 2) : 0; }
function looRMS(Y, X) {
  const n = Y.length; let ss = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xt = [...X.slice(0, i), ...X.slice(i + 1)];
    try {
      const f = ols(Yt, Xt);
      const xi = [1, ...X[i]];
      ss += (Y[i] - xi.reduce((s, x, j) => s + x * f.beta[j], 0)) ** 2;
    } catch (e) { ss += (Y[i] - mean(Y)) ** 2; }
  }
  return Math.sqrt(ss / n);
}

const coreGap = gapPct(looRMS(Y, core3X), sdY);
const coreVfResidGap = gapPct(looRMS(Y, gals45.map((_, i) => [...core3X[i], VfResid[i]])), sdY);
const coreUnexGap = gapPct(looRMS(Y, coreUnexplainedX), sdY);

console.log('  Can the unexplained part of VfResid still predict a₀?');
console.log('    Core gap: ' + coreGap.toFixed(1) + '%');
console.log('    Core + full VfResid gap: ' + coreVfResidGap.toFixed(1) + '%');
console.log('    Core + unexplained-VfResid gap: ' + coreUnexGap.toFixed(1) + '%');
console.log('    Unexplained still adds: ' + (coreUnexGap > coreGap ? '+' : '') + (coreUnexGap - coreGap).toFixed(1) + 'pp');
console.log('');

if (coreUnexGap > coreGap + 3) {
  console.log('  → DEEP RESIDUAL: Even after removing all known halo diagnostics,');
  console.log('    VfResid still predicts a₀. There is physics in VfResid that');
  console.log('    goes beyond what any available halo proxy captures.');
} else {
  console.log('  → FULLY EXPLAINED: The halo diagnostics capture essentially all');
  console.log('    of VfResid\'s predictive content for a₀.');
}
console.log('');

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 133C: VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const bestSingleR2 = singleResults.sort((a, b) => b.r2LOO - a.r2LOO)[0];
const bestMultiR2 = bestModel.r2LOO;
const unexpStillPredicts = coreUnexGap > coreGap + 3;

let verdict;
if (bestMultiR2 > 0.6 && !unexpStillPredicts) {
  verdict = 'COUPLING_LAW_FOUND';
} else if (bestMultiR2 > 0.4) {
  verdict = 'PARTIAL_COUPLING_LAW';
} else if (bestSingleR2.r2LOO > 0.2) {
  verdict = 'DOMINANT_DRIVER_IDENTIFIED';
} else {
  verdict = 'COUPLING_REMAINS_OPAQUE';
}

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + verdict.padEnd(55) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════════╝\n');

console.log('  COUPLING LAW SCORECARD:');
console.log('    Best single predictor: ' + bestSingleR2.name + ' (LOO R²=' + bestSingleR2.r2LOO + ')');
console.log('    Best multi-predictor: ' + bestModel.name + ' (LOO R²=' + bestMultiR2 + ')');
console.log('    Unexplained part still predicts a₀: ' + (unexpStillPredicts ? 'YES (+' + (coreUnexGap - coreGap).toFixed(1) + 'pp)' : 'NO'));
console.log('');

const output = {
  phase: '133C',
  title: 'What Drives VfResid? — Coupling Law Derivation',
  verdict,
  singlePredictors: singleResults,
  multiModels: multiResults,
  bestModel: bestModel,
  unexplainedAnalysis: {
    fractionExplained: +bestFit.r2.toFixed(3),
    coreGap: +coreGap.toFixed(1),
    coreVfResidGap: +coreVfResidGap.toFixed(1),
    coreUnexplainedGap: +coreUnexGap.toFixed(1),
    unexplainedStillPredicts: unexpStillPredicts
  }
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase133c-coupling-drivers.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase133c-coupling-drivers.json');
