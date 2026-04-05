#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 132B: MEDIATION / CAUSAL ORDERING TEST');
console.log('');
console.log('  Does VfResid mediate the halo → a₀ effect?');
console.log('  Or are they partially independent channels?');
console.log('  Testing: Halo → VfResid → a₀  vs  Halo → a₀ + VfResid → a₀');
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
function looCV(Y, X) {
  const n = Y.length; let ss = 0;
  for (let i = 0; i < n; i++) {
    const Yt = [...Y.slice(0, i), ...Y.slice(i + 1)];
    const Xt = [...X.slice(0, i), ...X.slice(i + 1)];
    const f = ols(Yt, Xt);
    const xi = [1, ...X[i]];
    ss += (Y[i] - xi.reduce((s, x, j) => s + x * f.beta[j], 0)) ** 2;
  }
  return Math.sqrt(ss / n);
}
function gapPct(rms, sdy) { return 100 * (1 - rms ** 2 / sdy ** 2); }

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const logMHI = gals45.map(g => g.logMHI);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logMR = gals45.map(g => g.logMeanRun);
const logVflat = gals45.map(g => Math.log10(sparcMap[g.name].Vflat));
const logMbar = gals45.map(g => {
  const s = sparcMap[g.name];
  return Math.log10(s.L36 * 0.5 * 1e9 + Math.pow(10, g.logMHI) * 1.33 * 1e9);
});
const logL36 = gals45.map(g => Math.log10(Math.max(sparcMap[g.name].L36, 0.01)));
const logRdisk = gals45.map(g => Math.log10(sparcMap[g.name].Rdisk));
const morphT = gals45.map(g => sparcMap[g.name].T);

const structX = gals45.map((_, i) => [logMbar[i], logL36[i], logRdisk[i], morphT[i]]);
const VfResid = ols(logVflat, structX).resid;

function safeLog(v) { return v > 0 ? Math.log10(v) : Math.log10(Math.max(Math.abs(v), 1)); }

function getHaloProxy(g, name) {
  const r = resMap[g.name];
  const s = sparcMap[g.name];
  const dhl = r.models.dark_halo_linear;
  const lh = r.models.log_halo;
  const mo = r.models.mond;
  const newt = r.models.newtonian;
  const mgh = r.models.modified_gravity_halo;
  const dhf = r.models.dark_halo_flat;
  switch (name) {
    case 'haloK_linear': return safeLog(Math.max(dhl.k, 1));
    case 'dhl_mse': return Math.log10(Math.max(dhl.mse, 0.1));
    case 'mond_improve': return mo.improvementVsNewton;
    case 'mond_mse_ratio': return Math.log10(Math.max(mo.mse / Math.max(newt.mse, 0.1), 0.001));
    case 'mgh_k': return safeLog(Math.max(mgh.k, 1));
    case 'lh_outerImprove': return lh.outerImprovement;
    case 'lh_improve': return lh.improvementVsNewton;
    case 'mond_mse': return Math.log10(Math.max(mo.mse, 0.1));
    case 'newt_mse': return Math.log10(Math.max(newt.mse, 0.1));
    case 'dhf_improve': return dhf.improvementVsNewton;
    case 'logSBeff': return Math.log10(Math.max(s.SBeff, 0.01));
    case 'concentration': return Math.log10(s.Reff / s.Rdisk);
    default: return 0;
  }
}

const core3X = gals45.map((_, i) => [logMHI[i], logMhost[i], logMR[i]]);
const coreGap = gapPct(looCV(Y, core3X), sdY);
const coreVfResidX = gals45.map((_, i) => [...core3X[i], VfResid[i]]);
const vfResidGap = gapPct(looCV(Y, coreVfResidX), sdY);

const haloCandidates = [
  'haloK_linear', 'dhl_mse', 'mond_improve', 'mond_mse_ratio',
  'mgh_k', 'lh_outerImprove', 'lh_improve', 'mond_mse',
  'newt_mse', 'dhf_improve', 'logSBeff', 'concentration'
];

console.log('  BASELINE: Core gap=' + coreGap.toFixed(1) + '%  Core+VfResid gap=' + vfResidGap.toFixed(1) + '%\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: MEDIATION QUADRUPLET — For each halo proxy');
console.log('  A: Core + halo    B: Core + VfResid    C: Core + halo + VfResid');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const mediationResults = [];

console.log('  ' + 'Halo Proxy'.padEnd(22) + '  A(+halo)  B(+VfR)  C(both)  absorbA  absorbB  Verdict');
console.log('  ' + '─'.repeat(85));

for (const hName of haloCandidates) {
  const hArr = gals45.map(g => getHaloProxy(g, hName));

  const coreHaloX = gals45.map((_, i) => [...core3X[i], hArr[i]]);
  const coreBothX = gals45.map((_, i) => [...core3X[i], hArr[i], VfResid[i]]);

  const gapA = gapPct(looCV(Y, coreHaloX), sdY);
  const gapB = vfResidGap;
  const gapC = gapPct(looCV(Y, coreBothX), sdY);

  const addedByHalo = gapA - coreGap;
  const addedByVfResid = gapB - coreGap;
  const absorbedA = gapC > gapB ? (gapC - gapB) : 0;
  const absorbedB = gapC > gapA ? (gapC - gapA) : 0;
  const absorbPctA = addedByHalo > 0 ? 100 * (1 - absorbedA / addedByHalo) : 0;
  const absorbPctB = addedByVfResid > 0 ? 100 * (1 - absorbedB / addedByVfResid) : 0;

  let verdict;
  if (absorbPctA > 70 && absorbPctB < 30) {
    verdict = 'FULL_MEDIATION';
  } else if (absorbPctA > 50 && absorbPctB < 50) {
    verdict = 'PARTIAL_MED';
  } else if (absorbPctA > 30 || absorbPctB > 30) {
    verdict = 'MUTUAL_PARTIAL';
  } else {
    verdict = 'INDEPENDENT';
  }

  mediationResults.push({
    name: hName, gapA, gapB, gapC,
    addedByHalo, addedByVfResid,
    absorbPctA: +absorbPctA.toFixed(1),
    absorbPctB: +absorbPctB.toFixed(1),
    verdict
  });

  console.log('  ' + hName.padEnd(22) + gapA.toFixed(1).padStart(7) + '%  ' +
    gapB.toFixed(1).padStart(6) + '%  ' + gapC.toFixed(1).padStart(6) + '%  ' +
    absorbPctA.toFixed(0).padStart(5) + '%  ' + absorbPctB.toFixed(0).padStart(6) + '%  ' + verdict);
}
console.log('');
console.log('  absorbA = how much of halo\'s added value is absorbed by VfResid');
console.log('  absorbB = how much of VfResid\'s added value is absorbed by halo');
console.log('  FULL_MEDIATION = VfResid absorbs >70% of halo, halo absorbs <30% of VfResid\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: SOBEL-LIKE MEDIATION PATH — a→b→c effect sizes');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const topHalo = ['haloK_linear', 'dhl_mse', 'mond_improve', 'mgh_k', 'lh_outerImprove'];

console.log('  Path: Halo → VfResid → logA0 (controlling for Core)\n');

for (const hName of topHalo) {
  const hArr = gals45.map(g => getHaloProxy(g, hName));

  const hResid = ols(hArr, core3X).resid;
  const vfResidC = ols(VfResid, core3X).resid;
  const yResid = ols(Y, core3X).resid;

  const rAB = pearsonR(hResid, vfResidC);
  const rBC = pearsonR(vfResidC, yResid);
  const rAC = pearsonR(hResid, yResid);

  const rAC_afterB = pearsonR(
    ols(hResid, vfResidC.map(v => [v])).resid,
    ols(yResid, vfResidC.map(v => [v])).resid
  );

  const indirect = rAB * rBC;
  const direct = rAC_afterB;
  const total = rAC;
  const medPct = total !== 0 ? 100 * indirect / total : 0;

  console.log('  ' + hName + ':');
  console.log('    Path a (halo→VfResid|Core):   r = ' + rAB.toFixed(3));
  console.log('    Path b (VfResid→a₀|Core):     r = ' + rBC.toFixed(3));
  console.log('    Total (halo→a₀|Core):         r = ' + rAC.toFixed(3));
  console.log('    Direct (halo→a₀|Core+VfResid): r = ' + direct.toFixed(3));
  console.log('    Indirect (a×b):                r = ' + indirect.toFixed(3));
  console.log('    Mediation %:                   ' + medPct.toFixed(1) + '%');
  console.log('    → ' + (medPct > 70 ? 'FULL MEDIATION' : medPct > 40 ? 'PARTIAL MEDIATION' : 'WEAK/NO MEDIATION'));
  console.log('');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: ENTRY-ORDER SENSITIVITY — Does adding order matter?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

for (const hName of topHalo) {
  const hArr = gals45.map(g => getHaloProxy(g, hName));

  const gapHaloOnly = gapPct(looCV(Y, gals45.map((_, i) => [...core3X[i], hArr[i]])), sdY);
  const gapVfResidOnly = vfResidGap;
  const gapHaloFirst = gapPct(looCV(Y, gals45.map((_, i) => [...core3X[i], hArr[i], VfResid[i]])), sdY);
  const gapVfResidFirst = gapPct(looCV(Y, gals45.map((_, i) => [...core3X[i], VfResid[i], hArr[i]])), sdY);

  const haloAddsToVfResid = gapVfResidFirst - gapVfResidOnly;
  const vfResidAddsToHalo = gapHaloFirst - gapHaloOnly;

  console.log('  ' + hName + ':');
  console.log('    Core + halo alone:        ' + gapHaloOnly.toFixed(1) + '% gap');
  console.log('    Core + VfResid alone:     ' + gapVfResidOnly.toFixed(1) + '% gap');
  console.log('    Core + halo + VfResid:    ' + gapHaloFirst.toFixed(1) + '% gap');
  console.log('    VfResid adds to halo:     +' + vfResidAddsToHalo.toFixed(1) + 'pp');
  console.log('    Halo adds to VfResid:     ' + (haloAddsToVfResid > 0 ? '+' : '') + haloAddsToVfResid.toFixed(1) + 'pp');
  console.log('    → ' + (vfResidAddsToHalo > 5 && haloAddsToVfResid < 3 ? 'VfResid ABSORBS halo signal'
    : vfResidAddsToHalo > haloAddsToVfResid * 2 ? 'VfResid DOMINANT channel'
      : Math.abs(vfResidAddsToHalo - haloAddsToVfResid) < 3 ? 'ROUGHLY SYMMETRIC'
        : 'HALO contributes independently'));
  console.log('');
}

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: BOOTSTRAP MEDIATION STABILITY');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const nBoot = 1000;
for (const hName of ['haloK_linear', 'dhl_mse', 'mgh_k']) {
  const hArr = gals45.map(g => getHaloProxy(g, hName));
  const medPcts = [];

  for (let b = 0; b < nBoot; b++) {
    const idx = Array.from({ length: N }, () => Math.floor(Math.random() * N));
    const Yb = idx.map(i => Y[i]);
    const coreXb = idx.map(i => core3X[i]);
    const hArrB = idx.map(i => hArr[i]);
    const structXb = idx.map(i => structX[i]);
    const logVflatB = idx.map(i => logVflat[i]);

    try {
      const vfResB = ols(logVflatB, structXb).resid;
      const hResid = ols(hArrB, coreXb).resid;
      const vfRC = ols(vfResB, coreXb).resid;
      const yResid = ols(Yb, coreXb).resid;

      const rAB = pearsonR(hResid, vfRC);
      const rBC = pearsonR(vfRC, yResid);
      const rAC = pearsonR(hResid, yResid);
      if (Math.abs(rAC) > 0.05) {
        medPcts.push(100 * rAB * rBC / rAC);
      }
    } catch (e) { }
  }

  medPcts.sort((a, b) => a - b);
  const p25 = medPcts[Math.floor(medPcts.length * 0.25)];
  const p50 = medPcts[Math.floor(medPcts.length * 0.5)];
  const p75 = medPcts[Math.floor(medPcts.length * 0.75)];
  const above50 = medPcts.filter(m => m > 50).length / medPcts.length;

  console.log('  ' + hName + ' mediation % bootstrap (' + medPcts.length + ' valid):');
  console.log('    25th: ' + p25.toFixed(1) + '%  50th: ' + p50.toFixed(1) + '%  75th: ' + p75.toFixed(1) + '%');
  console.log('    P(mediation > 50%): ' + (above50 * 100).toFixed(1) + '%');
  console.log('');
}

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 132B: VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const fullMed = mediationResults.filter(m => m.verdict === 'FULL_MEDIATION').length;
const partMed = mediationResults.filter(m => m.verdict === 'PARTIAL_MED').length;
const mutual = mediationResults.filter(m => m.verdict === 'MUTUAL_PARTIAL').length;
const indep = mediationResults.filter(m => m.verdict === 'INDEPENDENT').length;

let overallVerdict;
if (fullMed >= haloCandidates.length * 0.5) {
  overallVerdict = 'VFRESID_IS_PRIMARY_MEDIATOR';
} else if (fullMed + partMed >= haloCandidates.length * 0.5) {
  overallVerdict = 'VFRESID_DOMINANT_CHANNEL';
} else if (mutual >= haloCandidates.length * 0.5) {
  overallVerdict = 'PARALLEL_CHANNELS';
} else {
  overallVerdict = 'MIXED_PICTURE';
}

console.log('  ╔══════════════════════════════════════════════════════════════════╗');
console.log('  ║  VERDICT: ' + overallVerdict.padEnd(55) + '║');
console.log('  ╚══════════════════════════════════════════════════════════════════╝\n');

console.log('  MEDIATION TALLY:');
console.log('    Full mediation:     ' + fullMed + '/' + haloCandidates.length);
console.log('    Partial mediation:  ' + partMed + '/' + haloCandidates.length);
console.log('    Mutual partial:     ' + mutual + '/' + haloCandidates.length);
console.log('    Independent:        ' + indep + '/' + haloCandidates.length);
console.log('');

const output = {
  phase: '132B',
  title: 'Mediation / Causal Ordering Test',
  verdict: overallVerdict,
  mediationTally: { fullMed, partMed, mutual, indep, total: haloCandidates.length },
  mediationResults: mediationResults.map(m => ({
    name: m.name,
    gapA: +m.gapA.toFixed(1), gapB: +m.gapB.toFixed(1), gapC: +m.gapC.toFixed(1),
    absorbPctA: m.absorbPctA, absorbPctB: m.absorbPctB, verdict: m.verdict
  }))
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase132b-mediation-causal.json'),
  JSON.stringify(output, null, 2));
console.log('  Saved to public/phase132b-mediation-causal.json');
