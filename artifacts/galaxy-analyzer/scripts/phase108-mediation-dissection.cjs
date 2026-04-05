#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 108: MEDIATION / PROXY DISSECTION');
console.log('  Is fgas the real driver, or a proxy for something deeper?');
console.log('  Residualize each variable against the others, then test.');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }

function olsRegress(X, y) {
  const n = y.length, p = X[0].length;
  if (n <= p + 1) return null;
  const Xt = Array.from({ length: p }, (_, i) => X.map(row => row[i]));
  const XtX = Array.from({ length: p }, (_, i) =>
    Array.from({ length: p }, (_, j) =>
      Xt[i].reduce((s, _, k) => s + Xt[i][k] * Xt[j][k], 0)));
  const Xty = Xt.map(col => col.reduce((s, v, k) => s + v * y[k], 0));
  const aug = XtX.map((row, i) => [...row, Xty[i]]);
  for (let i = 0; i < p; i++) {
    let maxRow = i;
    for (let k = i + 1; k < p; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) maxRow = k;
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];
    if (Math.abs(aug[i][i]) < 1e-12) return null;
    for (let k = 0; k < p; k++) {
      if (k === i) continue;
      const f = aug[k][i] / aug[i][i];
      for (let j = i; j <= p; j++) aug[k][j] -= f * aug[i][j];
    }
  }
  const beta = aug.map((row, i) => row[p] / row[i]);
  const yhat = X.map(row => row.reduce((s, v, j) => s + v * beta[j], 0));
  return { beta, yhat, residuals: y.map((v, i) => v - yhat[i]) };
}

function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function looR2_arrays(xArr, yArr) {
  const n = xArr.length;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(yArr);
  for (let i = 0; i < n; i++) {
    const xTr = xArr.filter((_, j) => j !== i);
    const yTr = yArr.filter((_, j) => j !== i);
    const fit = olsRegress(xTr.map(v => [1, v]), yTr);
    if (!fit) return NaN;
    const pred = fit.beta[0] + fit.beta[1] * xArr[i];
    ss_loo += (yArr[i] - pred) ** 2;
    ss_tot += (yArr[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_loo / ss_tot : 0;
}

function residualize(target, predictors) {
  const X = predictors[0].map((_, i) => [1, ...predictors.map(p => p[i])]);
  const fit = olsRegress(X, target);
  return fit ? fit.residuals : target;
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  table1[name] = { L36, Rdisk, MHI, Vflat };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul });
}

const galaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat < 70) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    pts.push({ r: pt.rad, vobs: pt.vobs, vBar: Math.sqrt(Math.max(vBarSq, 0.01)) });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);
  const outerMD = outerPts.map(p => (p.vobs ** 2) / (p.vBar ** 2)).filter(v => isFinite(v) && v > 0);
  const logOMD = Math.log10(mean(outerMD));
  if (!isFinite(logOMD)) continue;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1);
  const baryonCompact = t1.L36 * UPSILON_DISK / (t1.Rdisk > 0 ? t1.Rdisk : 1);

  if (!isFinite(fgas) || !isFinite(Sigma0) || Sigma0 <= 0) continue;

  galaxies.push({
    name, logOMD, fgas,
    logSigma0: Math.log10(Sigma0),
    logL36: Math.log10(t1.L36),
    logRdisk: Math.log10(t1.Rdisk > 0 ? t1.Rdisk : 0.1),
    logMHI: Math.log10(t1.MHI),
    logBaryonCompact: Math.log10(baryonCompact > 0 ? baryonCompact : 1e-3),
  });
}

console.log('  N = ' + galaxies.length + '\n');

const vars = {
  fgas: galaxies.map(g => g.fgas),
  logSigma0: galaxies.map(g => g.logSigma0),
  logL36: galaxies.map(g => g.logL36),
  logRdisk: galaxies.map(g => g.logRdisk),
  logMHI: galaxies.map(g => g.logMHI),
  logBaryonCompact: galaxies.map(g => g.logBaryonCompact),
};
const yArr = galaxies.map(g => g.logOMD);

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: RESIDUALIZE fgas AGAINST EACH VARIABLE');
console.log('  fgas_perp = fgas after removing shared info with X');
console.log('  If fgas_perp still predicts logOMD -> fgas has unique content');
console.log('══════════════════════════════════════════════════════════════\n');

const rawR_fgas = pearsonR(vars.fgas, yArr);
const rawLOO_fgas = looR2_arrays(vars.fgas, yArr);
console.log('  RAW fgas: r=' + rawR_fgas.toFixed(3) + ', LOO=' + rawLOO_fgas.toFixed(3) + '\n');

const fgasPerpResults = {};
const controls = ['logSigma0', 'logL36', 'logRdisk', 'logMHI', 'logBaryonCompact'];

for (const cName of controls) {
  const fgas_perp = residualize(vars.fgas, [vars[cName]]);
  const r = pearsonR(fgas_perp, yArr);
  const loo = looR2_arrays(fgas_perp, yArr);
  const retained = (r / rawR_fgas * 100);

  console.log('  fgas | ' + cName + ':');
  console.log('    r(fgas_perp, logOMD) = ' + r.toFixed(3) + '  LOO = ' + loo.toFixed(3));
  console.log('    Signal retained: ' + retained.toFixed(1) + '% of raw r');
  console.log('');

  fgasPerpResults[cName] = { r: parseFloat(r.toFixed(3)), loo: parseFloat(loo.toFixed(3)), retained: parseFloat(retained.toFixed(1)) };
}

const fgas_perp_all = residualize(vars.fgas, controls.map(c => vars[c]));
const r_perp_all = pearsonR(fgas_perp_all, yArr);
const loo_perp_all = looR2_arrays(fgas_perp_all, yArr);
console.log('  fgas | ALL controls:');
console.log('    r(fgas_perp, logOMD) = ' + r_perp_all.toFixed(3) + '  LOO = ' + loo_perp_all.toFixed(3));
console.log('    Signal retained: ' + (r_perp_all / rawR_fgas * 100).toFixed(1) + '% of raw r\n');
fgasPerpResults['all'] = { r: parseFloat(r_perp_all.toFixed(3)), loo: parseFloat(loo_perp_all.toFixed(3)), retained: parseFloat((r_perp_all / rawR_fgas * 100).toFixed(1)) };

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: RESIDUALIZE EACH COMPETITOR AGAINST fgas');
console.log('  X_perp = X after removing shared info with fgas');
console.log('  If X_perp still predicts logOMD -> X has unique content');
console.log('══════════════════════════════════════════════════════════════\n');

const compPerpResults = {};
for (const cName of controls) {
  const rawR = pearsonR(vars[cName], yArr);
  const x_perp = residualize(vars[cName], [vars.fgas]);
  const r = pearsonR(x_perp, yArr);
  const loo = looR2_arrays(x_perp, yArr);
  const retained = rawR !== 0 ? (r / rawR * 100) : 0;

  console.log('  ' + cName + ' | fgas:');
  console.log('    raw r = ' + rawR.toFixed(3));
  console.log('    r(' + cName + '_perp, logOMD) = ' + r.toFixed(3) + '  LOO = ' + loo.toFixed(3));
  console.log('    Signal retained: ' + retained.toFixed(1) + '% of raw r');
  console.log('');

  compPerpResults[cName] = { rawR: parseFloat(rawR.toFixed(3)), perpR: parseFloat(r.toFixed(3)), loo: parseFloat(loo.toFixed(3)), retained: parseFloat(retained.toFixed(1)) };
}

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 3: SYMMETRIC MEDIATION — fgas vs logSigma0');
console.log('  The critical pair (r = -0.81 between them)');
console.log('══════════════════════════════════════════════════════════════\n');

const fgas_perpSig = residualize(vars.fgas, [vars.logSigma0]);
const sig_perpFgas = residualize(vars.logSigma0, [vars.fgas]);

const r_fgas_perpSig = pearsonR(fgas_perpSig, yArr);
const r_sig_perpFgas = pearsonR(sig_perpFgas, yArr);
const loo_fgas_perpSig = looR2_arrays(fgas_perpSig, yArr);
const loo_sig_perpFgas = looR2_arrays(sig_perpFgas, yArr);

console.log('  fgas | logSigma0:  r=' + r_fgas_perpSig.toFixed(3) + ', LOO=' + loo_fgas_perpSig.toFixed(3));
console.log('  logSigma0 | fgas:  r=' + r_sig_perpFgas.toFixed(3) + ', LOO=' + loo_sig_perpFgas.toFixed(3));
console.log('');
const sigDominates = Math.abs(r_sig_perpFgas) > Math.abs(r_fgas_perpSig);
console.log('  Winner: ' + (sigDominates ? 'logSigma0 has MORE unique content than fgas' : 'fgas has MORE unique content than logSigma0'));
console.log('  Ratio of unique signals: |r_fgas_perp|/|r_sig_perp| = ' + (Math.abs(r_fgas_perpSig) / Math.abs(r_sig_perpFgas)).toFixed(2));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: TRIPLE MEDIATION — fgas vs logL36 + logRdisk');
console.log('  fgas = MHI/(MHI + 0.5*L36). Is it just L36 and size?');
console.log('══════════════════════════════════════════════════════════════\n');

const fgas_perpLR = residualize(vars.fgas, [vars.logL36, vars.logRdisk]);
const r_fgas_perpLR = pearsonR(fgas_perpLR, yArr);
const loo_fgas_perpLR = looR2_arrays(fgas_perpLR, yArr);

console.log('  fgas | (logL36 + logRdisk):');
console.log('    r = ' + r_fgas_perpLR.toFixed(3) + ', LOO = ' + loo_fgas_perpLR.toFixed(3));
console.log('    Signal retained: ' + (r_fgas_perpLR / rawR_fgas * 100).toFixed(1) + '%');

const fgas_perpLRM = residualize(vars.fgas, [vars.logL36, vars.logRdisk, vars.logMHI]);
const r_fgas_perpLRM = pearsonR(fgas_perpLRM, yArr);
const loo_fgas_perpLRM = looR2_arrays(fgas_perpLRM, yArr);

console.log('\n  fgas | (logL36 + logRdisk + logMHI):');
console.log('    r = ' + r_fgas_perpLRM.toFixed(3) + ', LOO = ' + loo_fgas_perpLRM.toFixed(3));
console.log('    Signal retained: ' + (r_fgas_perpLRM / rawR_fgas * 100).toFixed(1) + '%');

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: DECOMPOSE fgas — MHI component vs L36 component');
console.log('  Which part of fgas = MHI/(MHI+0.5*L36) carries the signal?');
console.log('══════════════════════════════════════════════════════════════\n');

const r_MHI = pearsonR(vars.logMHI, yArr);
const r_L36 = pearsonR(vars.logL36, yArr);
const loo_MHI = looR2_arrays(vars.logMHI, yArr);
const loo_L36 = looR2_arrays(vars.logL36, yArr);

const logMHI_L36_ratio = galaxies.map(g => g.logMHI - g.logL36);
const r_ratio = pearsonR(logMHI_L36_ratio, yArr);
const loo_ratio = looR2_arrays(logMHI_L36_ratio, yArr);

console.log('  logMHI alone:        r=' + r_MHI.toFixed(3) + ', LOO=' + loo_MHI.toFixed(3));
console.log('  logL36 alone:        r=' + r_L36.toFixed(3) + ', LOO=' + loo_L36.toFixed(3));
console.log('  log(MHI/L36) ratio:  r=' + r_ratio.toFixed(3) + ', LOO=' + loo_ratio.toFixed(3));
console.log('  fgas:                r=' + rawR_fgas.toFixed(3) + ', LOO=' + rawLOO_fgas.toFixed(3));
console.log('');

const MHI_perpL36 = residualize(vars.logMHI, [vars.logL36]);
const L36_perpMHI = residualize(vars.logL36, [vars.logMHI]);
const r_MHI_perpL36 = pearsonR(MHI_perpL36, yArr);
const r_L36_perpMHI = pearsonR(L36_perpMHI, yArr);

console.log('  logMHI | logL36:  r=' + r_MHI_perpL36.toFixed(3) + ' (unique gas content)');
console.log('  logL36 | logMHI:  r=' + r_L36_perpMHI.toFixed(3) + ' (unique stellar luminosity)');
console.log('  -> ' + (Math.abs(r_MHI_perpL36) > Math.abs(r_L36_perpMHI) ? 'Gas mass carries more unique signal' : 'Stellar luminosity carries more unique signal'));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  OVERALL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const fgasSurvivesAll = Math.abs(r_perp_all) > 0.2;
const fgasSurvivesSigma = Math.abs(r_fgas_perpSig) > 0.3;
const sigmaDies = Math.abs(r_sig_perpFgas) < 0.15;
const fgasSurvivesLR = Math.abs(r_fgas_perpLR) > 0.2;

let verdict;
if (fgasSurvivesSigma && sigmaDies) {
  verdict = 'CAUSAL CANDIDATE: fgas has unique predictive content that survives controlling for Sigma0. Sigma0 is fully mediated by fgas. fgas is not a proxy for surface density.';
} else if (fgasSurvivesSigma && !sigmaDies) {
  verdict = 'INDEPENDENT SIGNALS: Both fgas and Sigma0 carry unique information. The outer support requirement has two independent drivers.';
} else if (!fgasSurvivesSigma && sigmaDies) {
  verdict = 'MUTUAL PROXY: Both fgas and Sigma0 are proxies for the same underlying variable. Neither has unique content.';
} else {
  verdict = 'PROXY DETECTED: fgas loses its signal when Sigma0 is controlled. fgas is likely a proxy for diffuse baryonic structure (low surface density).';
}

console.log('  fgas_perp survives Sigma0? ' + (fgasSurvivesSigma ? 'YES (r=' + r_fgas_perpSig.toFixed(3) + ')' : 'NO (r=' + r_fgas_perpSig.toFixed(3) + ')'));
console.log('  Sigma0_perp survives fgas? ' + (!sigmaDies ? 'YES (r=' + r_sig_perpFgas.toFixed(3) + ')' : 'NO (r=' + r_sig_perpFgas.toFixed(3) + ')'));
console.log('  fgas survives ALL controls? ' + (fgasSurvivesAll ? 'YES' : 'NO') + ' (r=' + r_perp_all.toFixed(3) + ')');
console.log('  fgas survives L36+Rdisk? ' + (fgasSurvivesLR ? 'YES' : 'NO') + ' (r=' + r_fgas_perpLR.toFixed(3) + ')');
console.log('');
console.log('  VERDICT: ' + verdict);

const output = {
  phase: 108,
  title: 'Mediation / Proxy Dissection',
  n: galaxies.length,
  rawFgas: { r: parseFloat(rawR_fgas.toFixed(3)), loo: parseFloat(rawLOO_fgas.toFixed(3)) },
  fgasPerpResults,
  competitorPerpResults: compPerpResults,
  symmetricMediation: {
    fgas_perpSigma0: { r: parseFloat(r_fgas_perpSig.toFixed(3)), loo: parseFloat(loo_fgas_perpSig.toFixed(3)) },
    sigma0_perpFgas: { r: parseFloat(r_sig_perpFgas.toFixed(3)), loo: parseFloat(loo_sig_perpFgas.toFixed(3)) },
  },
  fgasDecomposition: {
    logMHI: { r: parseFloat(r_MHI.toFixed(3)), loo: parseFloat(loo_MHI.toFixed(3)) },
    logL36: { r: parseFloat(r_L36.toFixed(3)), loo: parseFloat(loo_L36.toFixed(3)) },
    logMHI_L36_ratio: { r: parseFloat(r_ratio.toFixed(3)), loo: parseFloat(loo_ratio.toFixed(3)) },
    logMHI_perpL36: parseFloat(r_MHI_perpL36.toFixed(3)),
    logL36_perpMHI: parseFloat(r_L36_perpMHI.toFixed(3)),
  },
  verdict
};

const outPath = path.join(__dirname, '..', 'public', 'phase108-mediation-dissection.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);
