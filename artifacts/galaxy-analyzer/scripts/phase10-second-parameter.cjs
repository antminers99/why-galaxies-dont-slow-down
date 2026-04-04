#!/usr/bin/env node
/**
 * Phase 10: Second-Parameter Search
 * v10.0.0
 *
 * Tests whether adding a second variable X to the RAR reduces galaxy-level tau.
 *
 * Model 0: logA0_eff = mu                    (single universal a0)
 * Model 1: logA0_eff = mu + beta * X_norm    (a0 depends on second parameter)
 *
 * Candidates:
 *   1. SBdisk   — disk central surface brightness (structure/EFE proxy)
 *   2. SBeff    — effective surface brightness
 *   3. inc      — inclination (kinematic correction proxy)
 *   4. logMHI   — log HI gas mass (gas disequilibrium proxy)
 *   5. T        — Hubble type (morphology/bar proxy)
 *   6. logL36   — log 3.6um luminosity (mass proxy)
 *   7. Rdisk    — disk scale length (size proxy)
 *   8. fD       — distance method flag (systematic proxy)
 *   9. Q        — SPARC quality flag
 *  10. logSigmaBar — log baryonic surface density
 *
 * For each: fit beta, compute tau_new, delta-tau, permutation p-value,
 *           AIC improvement, jackknife stability.
 */

const fs = require('fs');
const path = require('path');

const rar = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/rar-analysis-real.json'), 'utf8'));
const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));

const A0_REF = 3703;

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0(pts) {
  let bestA0 = 3000, bestChi2 = Infinity;
  for (let logA = 2.5; logA <= 4.5; logA += 0.005) {
    const a0 = Math.pow(10, logA);
    let chi2 = 0;
    for (const p of pts) {
      const gbar = Math.pow(10, p.log_g_bar);
      const pred = mcgaughRAR(gbar, a0);
      if (!isFinite(pred) || pred <= 0) { chi2 += 100; continue; }
      chi2 += (p.log_g_obs - Math.log10(pred)) ** 2;
    }
    if (chi2 < bestChi2) { bestChi2 = chi2; bestA0 = a0; }
  }
  let lo = Math.log10(bestA0) - 0.01, hi = Math.log10(bestA0) + 0.01;
  for (let s = 0; s < 80; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gbar = Math.pow(10, p.log_g_bar);
      c1 += (p.log_g_obs - Math.log10(mcgaughRAR(gbar, Math.pow(10, m1)))) ** 2;
      c2 += (p.log_g_obs - Math.log10(mcgaughRAR(gbar, Math.pow(10, m2)))) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  return Math.pow(10, (lo + hi) / 2);
}

function fitPointLevel(pts, a0) {
  let ss = 0;
  for (const p of pts) {
    const gbar = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gbar, a0);
    if (!isFinite(pred) || pred <= 0) continue;
    ss += (p.log_g_obs - Math.log10(pred)) ** 2;
  }
  return ss;
}

const galPts = {};
for (const p of rar.rarScatter) {
  if (!galPts[p.name]) galPts[p.name] = [];
  galPts[p.name].push(p);
}

const sparcLookup = {};
for (const g of sparc) sparcLookup[g.name] = g;

const goldI45 = [];
for (const g of rar.perGalaxy) {
  const pts = galPts[g.name] || [];
  if (pts.length < 5 || g.Vmax < 50 || g.inc < 45) continue;
  const gRange = Math.max(...pts.map(p => p.log_g_bar)) - Math.min(...pts.map(p => p.log_g_bar));
  if (gRange < 1.0) continue;
  const a0 = fitA0(pts);
  if (a0 < 10 || a0 > 1e5) continue;
  const sp = sparcLookup[g.name] || {};
  goldI45.push({
    name: g.name,
    a0,
    logA0: Math.log10(a0),
    pts,
    n: pts.length,
    Vmax: g.Vmax,
    inc: g.inc,
    SBdisk: sp.SBdisk || null,
    SBeff: sp.SBeff || null,
    T: sp.T != null ? sp.T : null,
    L36: sp.L36 || null,
    MHI: sp.MHI || null,
    Rdisk: sp.Rdisk || null,
    fD: sp.fD || null,
    Q: sp.Q || null,
    sigma_bar: g.sigma_bar || null,
    eta_rot: g.eta_rot || null,
    S_out: g.S_out || null,
  });
}

const logA0s = goldI45.map(g => g.logA0);
const meanLogA0 = logA0s.reduce((a, b) => a + b, 0) / logA0s.length;
const baselineTau = Math.sqrt(logA0s.reduce((a, v) => a + (v - meanLogA0) ** 2, 0) / logA0s.length);

console.log(`\n${'='.repeat(70)}`);
console.log(`  PHASE 10: SECOND-PARAMETER SEARCH`);
console.log(`${'='.repeat(70)}`);
console.log(`  Galaxies: ${goldI45.length} GOLD+i45 (benchmark)`);
console.log(`  Baseline tau (SD logA0): ${baselineTau.toFixed(4)} dex`);
console.log(`  Baseline mean logA0: ${meanLogA0.toFixed(4)}\n`);

function normalize(vals) {
  const valid = vals.filter(v => v !== null && isFinite(v));
  const mu = valid.reduce((a, b) => a + b, 0) / valid.length;
  const sd = Math.sqrt(valid.reduce((a, v) => a + (v - mu) ** 2, 0) / valid.length);
  return { mu, sd };
}

function buildCandidate(name, extractor) {
  const vals = goldI45.map(g => {
    const v = extractor(g);
    return (v !== null && v !== undefined && isFinite(v) && v > 0) ? v : null;
  });
  const nValid = vals.filter(v => v !== null).length;
  if (nValid < 20) return null;
  return { name, vals, nValid };
}

const candidates = [
  buildCandidate('SBdisk', g => g.SBdisk),
  buildCandidate('SBeff', g => g.SBeff),
  buildCandidate('inc', g => g.inc),
  buildCandidate('logMHI', g => g.MHI > 0 ? Math.log10(g.MHI) : null),
  buildCandidate('T', g => g.T),
  buildCandidate('logL36', g => g.L36 > 0 ? Math.log10(g.L36) : null),
  buildCandidate('Rdisk', g => g.Rdisk),
  buildCandidate('logSigmaBar', g => g.sigma_bar > 0 ? Math.log10(g.sigma_bar) : null),
  buildCandidate('eta_rot', g => g.eta_rot),
  buildCandidate('S_out', g => g.S_out),
].filter(Boolean);

function fitLinear(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let ssxy = 0, ssxx = 0;
  for (let i = 0; i < n; i++) {
    ssxy += (x[i] - mx) * (y[i] - my);
    ssxx += (x[i] - mx) ** 2;
  }
  const beta = ssxy / ssxx;
  const alpha = my - beta * mx;
  const residuals = y.map((yi, i) => yi - alpha - beta * x[i]);
  const ssRes = residuals.reduce((a, r) => a + r * r, 0);
  const ssTot = y.reduce((a, yi) => a + (yi - my) ** 2, 0);
  const r2 = 1 - ssRes / ssTot;
  const tauNew = Math.sqrt(ssRes / n);
  const seBeta = Math.sqrt(ssRes / (n - 2) / ssxx);
  const tStat = beta / seBeta;
  const r = Math.sqrt(Math.max(0, r2)) * Math.sign(beta);
  return { alpha, beta, seBeta, tStat, r2, r, tauNew, residuals };
}

function permutationTest(x, y, observedR2, nPerm) {
  let count = 0;
  const yArr = [...y];
  for (let p = 0; p < nPerm; p++) {
    for (let i = yArr.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [yArr[i], yArr[j]] = [yArr[j], yArr[i]];
    }
    const fit = fitLinear(x, yArr);
    if (fit.r2 >= observedR2) count++;
  }
  return count / nPerm;
}

function jackknifeBeta(x, y) {
  const n = x.length;
  const betas = [];
  for (let i = 0; i < n; i++) {
    const xj = x.filter((_, j) => j !== i);
    const yj = y.filter((_, j) => j !== i);
    betas.push(fitLinear(xj, yj).beta);
  }
  const meanB = betas.reduce((a, b) => a + b, 0) / n;
  const seJK = Math.sqrt((n - 1) / n * betas.reduce((a, b) => a + (b - meanB) ** 2, 0));
  return { betas, meanB, seJK, range: [Math.min(...betas), Math.max(...betas)] };
}

function computeAIC(n, ssRes, k) {
  return n * Math.log(ssRes / n) + 2 * k;
}

function fitPointLevelModel1(galaxies, x, beta, globalMu) {
  let ssTotalM0 = 0, ssTotalM1 = 0, nPts = 0;
  for (let i = 0; i < galaxies.length; i++) {
    if (x[i] === null) continue;
    const g = galaxies[i];
    const a0_m0 = Math.pow(10, globalMu);
    const a0_m1 = Math.pow(10, globalMu + beta * x[i]);
    for (const p of g.pts) {
      const gbar = Math.pow(10, p.log_g_bar);
      const pred0 = mcgaughRAR(gbar, a0_m0);
      const pred1 = mcgaughRAR(gbar, a0_m1);
      if (!isFinite(pred0) || pred0 <= 0 || !isFinite(pred1) || pred1 <= 0) continue;
      ssTotalM0 += (p.log_g_obs - Math.log10(pred0)) ** 2;
      ssTotalM1 += (p.log_g_obs - Math.log10(pred1)) ** 2;
      nPts++;
    }
  }
  return { rmsM0: Math.sqrt(ssTotalM0 / nPts), rmsM1: Math.sqrt(ssTotalM1 / nPts), nPts };
}

console.log(`${'─'.repeat(70)}`);
console.log('  CANDIDATE SECOND PARAMETERS');
console.log(`${'─'.repeat(70)}\n`);

const results = [];

for (const cand of candidates) {
  const validIdx = [];
  const xRaw = [];
  const yRaw = [];
  for (let i = 0; i < goldI45.length; i++) {
    if (cand.vals[i] !== null) {
      validIdx.push(i);
      xRaw.push(cand.vals[i]);
      yRaw.push(goldI45[i].logA0);
    }
  }

  const { mu: xMu, sd: xSD } = normalize(xRaw);
  const xNorm = xRaw.map(v => (v - xMu) / xSD);

  const fit = fitLinear(xNorm, yRaw);
  const pPerm = permutationTest(xNorm, yRaw, fit.r2, 2000);
  const jk = jackknifeBeta(xNorm, yRaw);

  const n = xRaw.length;
  const ssTotY = yRaw.reduce((a, y) => a + (y - yRaw.reduce((s, v) => s + v, 0) / n) ** 2, 0);
  const ssResM0 = ssTotY;
  const ssResM1 = fit.residuals.reduce((a, r) => a + r * r, 0);
  const aicM0 = computeAIC(n, ssResM0, 1);
  const aicM1 = computeAIC(n, ssResM1, 2);
  const deltaAIC = aicM1 - aicM0;

  const betaRaw = fit.beta / xSD;
  const xNormForPtLevel = cand.vals.map(v => v !== null ? (v - xMu) / xSD : null);
  const ptLevel = fitPointLevelModel1(goldI45, xNormForPtLevel, fit.beta, fit.alpha);

  const deltaTau = fit.tauNew - baselineTau;
  const deltaTauPct = (deltaTau / baselineTau) * 100;

  const result = {
    name: cand.name,
    nValid: n,
    beta: +fit.beta.toFixed(4),
    seBeta: +fit.seBeta.toFixed(4),
    tStat: +fit.tStat.toFixed(3),
    r: +fit.r.toFixed(4),
    r2: +fit.r2.toFixed(4),
    pPerm: +pPerm.toFixed(4),
    baselineTau: +baselineTau.toFixed(4),
    tauNew: +fit.tauNew.toFixed(4),
    deltaTauPct: +deltaTauPct.toFixed(1),
    deltaAIC: +deltaAIC.toFixed(2),
    jkMeanBeta: +jk.meanB.toFixed(4),
    jkSEBeta: +jk.seJK.toFixed(4),
    jkRange: [+jk.range[0].toFixed(4), +jk.range[1].toFixed(4)],
    pointRMS_M0: +ptLevel.rmsM0.toFixed(4),
    pointRMS_M1: +ptLevel.rmsM1.toFixed(4),
    pointImprPct: +(((ptLevel.rmsM0 - ptLevel.rmsM1) / ptLevel.rmsM0) * 100).toFixed(2),
    significant: pPerm < 0.05 && Math.abs(fit.tStat) >= 2.0 && deltaAIC < -2,
  };
  results.push(result);

  const sig = result.significant ? '*** SIGNIFICANT ***' : '';
  console.log(`  ${cand.name} (n=${n})`);
  console.log(`    beta = ${fit.beta.toFixed(4)} +/- ${fit.seBeta.toFixed(4)}, t = ${fit.tStat.toFixed(2)}`);
  console.log(`    r = ${fit.r.toFixed(3)}, R2 = ${(fit.r2 * 100).toFixed(1)}%`);
  console.log(`    p(perm) = ${pPerm.toFixed(4)} ${sig}`);
  console.log(`    tau: ${baselineTau.toFixed(4)} -> ${fit.tauNew.toFixed(4)} (${deltaTauPct.toFixed(1)}%)`);
  console.log(`    deltaAIC = ${deltaAIC.toFixed(2)}`);
  console.log(`    JK beta: ${jk.meanB.toFixed(4)} +/- ${jk.seJK.toFixed(4)} [${jk.range[0].toFixed(3)}, ${jk.range[1].toFixed(3)}]`);
  console.log(`    Point RMS: ${ptLevel.rmsM0.toFixed(4)} -> ${ptLevel.rmsM1.toFixed(4)} (${result.pointImprPct}%)`);
  console.log('');
}

results.sort((a, b) => a.deltaTauPct - b.deltaTauPct);

console.log(`${'─'.repeat(70)}`);
console.log('  RANKING BY TAU REDUCTION');
console.log(`${'─'.repeat(70)}\n`);

console.log('  Rank  Candidate      tau_new   delta%   r       p(perm)  dAIC   Sig?');
console.log('  ' + '─'.repeat(68));
results.forEach((r, i) => {
  console.log(`  ${String(i + 1).padStart(2)}.   ${r.name.padEnd(15)} ${r.tauNew.toFixed(4)}   ${(r.deltaTauPct >= 0 ? '+' : '') + r.deltaTauPct.toFixed(1).padStart(5)}%  ${r.r.toFixed(3).padStart(6)}  ${r.pPerm.toFixed(3).padStart(6)}   ${r.deltaAIC.toFixed(1).padStart(5)}  ${r.significant ? 'YES' : ' no'}`);
});

const bestCandidate = results[0];
const anySignificant = results.filter(r => r.significant);

console.log(`\n${'─'.repeat(70)}`);
console.log('  MULTI-PARAMETER MODEL TEST');
console.log(`${'─'.repeat(70)}\n`);

const top3 = results.slice(0, 3);
const multiX = [];
const multiY = [];
const multiNames = top3.map(r => r.name);

for (let i = 0; i < goldI45.length; i++) {
  const g = goldI45[i];
  const sp = sparcLookup[g.name] || {};
  const xVals = multiNames.map(name => {
    const cand = candidates.find(c => c.name === name);
    return cand ? cand.vals[i] : null;
  });
  if (xVals.some(v => v === null)) continue;
  multiX.push(xVals);
  multiY.push(g.logA0);
}

if (multiX.length >= 20) {
  const nM = multiX.length;
  const k = multiNames.length;
  const means = multiNames.map((_, j) => multiX.reduce((s, row) => s + row[j], 0) / nM);
  const sds = multiNames.map((_, j) => {
    const mu = means[j];
    return Math.sqrt(multiX.reduce((s, row) => s + (row[j] - mu) ** 2, 0) / nM);
  });
  const xNormed = multiX.map(row => row.map((v, j) => (v - means[j]) / sds[j]));

  const my = multiY.reduce((a, b) => a + b, 0) / nM;
  const XtX = Array.from({ length: k + 1 }, () => new Array(k + 1).fill(0));
  const XtY = new Array(k + 1).fill(0);
  for (let i = 0; i < nM; i++) {
    const xi = [1, ...xNormed[i]];
    for (let a = 0; a <= k; a++) {
      for (let b = 0; b <= k; b++) XtX[a][b] += xi[a] * xi[b];
      XtY[a] += xi[a] * multiY[i];
    }
  }

  function solveGauss(A, b) {
    const n = b.length;
    const aug = A.map((row, i) => [...row, b[i]]);
    for (let col = 0; col < n; col++) {
      let maxRow = col;
      for (let row = col + 1; row < n; row++) {
        if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row;
      }
      [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]];
      if (Math.abs(aug[col][col]) < 1e-12) continue;
      for (let row = col + 1; row < n; row++) {
        const f = aug[row][col] / aug[col][col];
        for (let j = col; j <= n; j++) aug[row][j] -= f * aug[col][j];
      }
    }
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      x[i] = aug[i][n];
      for (let j = i + 1; j < n; j++) x[i] -= aug[i][j] * x[j];
      x[i] /= aug[i][i];
    }
    return x;
  }

  const betas = solveGauss(XtX, XtY);
  const residuals = multiY.map((yi, i) => {
    const xi = [1, ...xNormed[i]];
    const pred = xi.reduce((s, v, j) => s + v * betas[j], 0);
    return yi - pred;
  });
  const ssResMulti = residuals.reduce((a, r) => a + r * r, 0);
  const tauMulti = Math.sqrt(ssResMulti / nM);
  const r2Multi = 1 - ssResMulti / multiY.reduce((a, y) => a + (y - my) ** 2, 0);
  const aicMulti = computeAIC(nM, ssResMulti, k + 1);
  const aicM0Multi = computeAIC(nM, multiY.reduce((a, y) => a + (y - my) ** 2, 0), 1);

  console.log(`  Top-3 combined model: ${multiNames.join(' + ')}`);
  console.log(`  n = ${nM}, k = ${k}`);
  console.log(`  Betas: ${betas.map((b, i) => (i === 0 ? 'intercept' : multiNames[i - 1]) + '=' + b.toFixed(4)).join(', ')}`);
  console.log(`  R2 = ${(r2Multi * 100).toFixed(1)}%`);
  console.log(`  tau: ${baselineTau.toFixed(4)} -> ${tauMulti.toFixed(4)} (${((tauMulti - baselineTau) / baselineTau * 100).toFixed(1)}%)`);
  console.log(`  deltaAIC = ${(aicMulti - aicM0Multi).toFixed(2)}`);
}

console.log(`\n${'═'.repeat(70)}`);
console.log('  PHASE 10 VERDICT');
console.log(`${'═'.repeat(70)}\n`);

if (anySignificant.length === 0) {
  console.log('  No single second parameter significantly reduces tau.');
  console.log('  All candidates: p(perm) > 0.05 or |t| < 2 or deltaAIC > -2.');
  console.log('  The galaxy-level heterogeneity is NOT explained by any');
  console.log('  measured structural or observational property.');
  console.log('');
  console.log('  IMPLICATION: Either:');
  console.log('  (a) The second parameter is unmeasured (e.g., EFE, halo properties)');
  console.log('  (b) The heterogeneity is genuinely intrinsic');
  console.log('  (c) Multiple small effects combine nonlinearly');
} else {
  console.log(`  ${anySignificant.length} candidate(s) significantly reduce tau:`);
  for (const r of anySignificant) {
    console.log(`    ${r.name}: tau ${r.baselineTau.toFixed(4)} -> ${r.tauNew.toFixed(4)} (${r.deltaTauPct.toFixed(1)}%), p=${r.pPerm.toFixed(4)}`);
  }
}

const output = {
  version: '10.0.0',
  timestamp: new Date().toISOString(),
  description: 'Phase 10: Second-parameter search for galaxy-level tau reduction',
  nGalaxies: goldI45.length,
  baselineTau: +baselineTau.toFixed(4),
  baselineMeanLogA0: +meanLogA0.toFixed(4),
  candidates: results,
  anySignificant: anySignificant.length > 0,
  significantCandidates: anySignificant.map(r => r.name),
  ranking: results.map(r => r.name),
  verdict: anySignificant.length === 0
    ? 'No single measured second parameter significantly reduces galaxy-level tau. Heterogeneity remains unexplained by available observational properties.'
    : `${anySignificant.length} candidate(s) significantly reduce tau: ${anySignificant.map(r => r.name).join(', ')}`,
};

fs.writeFileSync(
  path.join(__dirname, '../public/phase10-second-param-results.json'),
  JSON.stringify(output, null, 2)
);
console.log('\n  Results written to public/phase10-second-param-results.json');
console.log(`${'═'.repeat(70)}\n`);
