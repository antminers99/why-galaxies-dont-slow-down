const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

const G = 4.3009e-6;

function normalize(n) {
  return n.replace(/\s+/g, '').replace(/^NGC0*/, 'NGC').replace(/^DDO0*/, 'DDO').replace(/^IC0*/, 'IC').replace(/^UGC0*/, 'UGC').replace(/^PGC0*/, 'PGC').toUpperCase();
}

function pearsonR(x, y) {
  const n = x.length; if (n < 4) return NaN;
  const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; num += dx * dy; dx2 += dx * dx; dy2 += dy * dy; }
  return dx2 > 0 && dy2 > 0 ? num / Math.sqrt(dx2 * dy2) : 0;
}

function multiR2(X, y) {
  const n = y.length, nv = X[0].length;
  const my = y.reduce((a, b) => a + b, 0) / n;
  const mx = Array(nv).fill(0);
  for (let j = 0; j < nv; j++) { for (let i = 0; i < n; i++) mx[j] += X[i][j]; mx[j] /= n; }
  const XTX = Array.from({ length: nv }, () => Array(nv).fill(0)), XTy = Array(nv).fill(0);
  for (let i = 0; i < n; i++) { for (let j = 0; j < nv; j++) { XTy[j] += (X[i][j] - mx[j]) * (y[i] - my); for (let k = 0; k < nv; k++) XTX[j][k] += (X[i][j] - mx[j]) * (X[i][k] - mx[k]); } }
  const aug = XTX.map((row, i) => [...row, XTy[i]]);
  for (let col = 0; col < nv; col++) { let maxRow = col; for (let row = col + 1; row < nv; row++) if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row; [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]]; if (Math.abs(aug[col][col]) < 1e-12) continue; for (let row = col + 1; row < nv; row++) { const f = aug[row][col] / aug[col][col]; for (let j = col; j <= nv; j++) aug[row][j] -= f * aug[col][j]; } }
  const beta = Array(nv).fill(0);
  for (let i = nv - 1; i >= 0; i--) { beta[i] = aug[i][nv]; for (let j = i + 1; j < nv; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= aug[i][i] || 1; }
  const resid = [];
  for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < nv; j++) pred += beta[j] * (X[i][j] - mx[j]); resid.push(y[i] - pred); }
  return { residuals: resid, beta };
}

function partialR(x, y, controls) {
  const n = x.length;
  if (n < controls[0].length + 4) return NaN;
  const residualise = (v, C) => {
    const nv = C[0].length;
    const my = v.reduce((a, b) => a + b, 0) / n;
    const mc = Array(nv).fill(0);
    for (let j = 0; j < nv; j++) { for (let i = 0; i < n; i++) mc[j] += C[i][j]; mc[j] /= n; }
    const XTX = Array.from({ length: nv }, () => Array(nv).fill(0));
    const XTy = Array(nv).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < nv; j++) {
        XTy[j] += (C[i][j] - mc[j]) * (v[i] - my);
        for (let k = 0; k < nv; k++) XTX[j][k] += (C[i][j] - mc[j]) * (C[i][k] - mc[k]);
      }
    }
    const aug2 = XTX.map((row, i) => [...row, XTy[i]]);
    for (let col = 0; col < nv; col++) {
      let maxRow = col;
      for (let row = col + 1; row < nv; row++) if (Math.abs(aug2[row][col]) > Math.abs(aug2[maxRow][col])) maxRow = row;
      [aug2[col], aug2[maxRow]] = [aug2[maxRow], aug2[col]];
      if (Math.abs(aug2[col][col]) < 1e-12) continue;
      for (let row = col + 1; row < nv; row++) { const f = aug2[row][col] / aug2[col][col]; for (let j = col; j <= nv; j++) aug2[row][j] -= f * aug2[col][j]; }
    }
    const beta2 = Array(nv).fill(0);
    for (let i = nv - 1; i >= 0; i--) { beta2[i] = aug2[i][nv]; for (let j = i + 1; j < nv; j++) beta2[i] -= aug2[i][j] * beta2[j]; beta2[i] /= aug2[i][i] || 1; }
    const r = [];
    for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < nv; j++) pred += beta2[j] * (C[i][j] - mc[j]); r.push(v[i] - pred); }
    return r;
  };
  return pearsonR(residualise(x, controls), residualise(y, controls));
}

const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[g.name] = g; resultsMap[normalize(g.name)] = g; });

const tsContent = fs.readFileSync(path.join(__dirname, '..', 'src', 'data', 'sparc-datasets.ts'), 'utf8');
const rcMapRaw = {};
const re = /"([^"]+)":\s*\{[^[]*data:\s*\[([\s\S]*?)\]/g;
let m;
while ((m = re.exec(tsContent)) !== null) {
  const pts = [];
  const ptRe = /r:\s*([\d.]+)\s*,\s*v:\s*([\d.]+)/g;
  let pm;
  while ((pm = ptRe.exec(m[2])) !== null) pts.push({ r: parseFloat(pm[1]), v: parseFloat(pm[2]) });
  if (pts.length >= 3) rcMapRaw[m[1]] = pts;
}
const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) rcMap[normalize(name)] = pts;

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sr) continue;
  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 5) continue;
  const Vflat = sp.Vflat, Rdisk = sp.Rdisk;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  const Rmax = rc[rc.length - 1].r;
  const V_Newt = Math.sqrt(G * Mbar / Rmax);
  const dmFrac = Math.max(0, 1 - (V_Newt / Vflat) ** 2);
  const logK = sr.models.dark_halo_linear.k > 0 ? Math.log10(sr.models.dark_halo_linear.k) : -5;
  const hR = sr.models.newtonian.mse > 0 ? Math.log10(Math.max(sr.models.newtonian.mse / Math.max(sr.models.dark_halo_linear.mse, 0.001), 0.01)) : 0;
  const logMbar = Math.log10(Math.max(Mbar, 1));
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(Rdisk, 0.01));
  const logSBdisk = Math.log10(Math.max(sp.SBdisk, 0.01));
  const outerPts = rc.filter(p => p.r > 3 * Rdisk);
  let outerSlope = 0;
  if (outerPts.length >= 3) {
    const mx2 = outerPts.reduce((a, p) => a + p.r, 0) / outerPts.length;
    const my2 = outerPts.reduce((a, p) => a + p.v, 0) / outerPts.length;
    let n2 = 0, d2 = 0;
    for (const p of outerPts) { n2 += (p.r - mx2) * (p.v - my2); d2 += (p.r - mx2) ** 2; }
    outerSlope = d2 > 0 ? n2 / d2 : 0;
  }
  let rcSmooth = 0;
  if (rc.length >= 5) {
    const sm = [];
    for (let i = 2; i < rc.length - 2; i++) sm.push((rc[i - 2].v + rc[i - 1].v + rc[i].v + rc[i + 1].v + rc[i + 2].v) / 5);
    let ss = 0;
    for (let i = 0; i < sm.length; i++) ss += (rc[i + 2].v - sm[i]) ** 2;
    rcSmooth = 1 - Math.sqrt(ss / sm.length) / Vflat;
  }
  const MHI = Math.pow(10, g.logMHI || 8);
  const Mstar = Math.pow(10, logL36) * 0.5e9;
  const gasFrac = MHI / (MHI + Mstar);

  gals.push({
    name: g.name, Vflat, Rdisk, Rmax, Mbar, rc,
    logVflat: Math.log10(Vflat), logL36, logRdisk, logMbar,
    logMHI: g.logMHI, morphT: sp.T,
    logSBdisk, envCode: g.envCode,
    logK, dmFrac, hR, outerSlope, rcSmooth, gasFrac,
    Npts: rc.length, D: sp.D,
  });
}

const N = gals.length;
const struct4 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
const struct6 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
const vfM = multiR2(struct4, gals.map(g => g.logVflat));
const a0M = multiR2(struct6, gals.map(g => g.logA0));
for (let i = 0; i < N; i++) { gals[i].VfR = vfM.residuals[i]; gals[i].a0R = a0M.residuals[i]; }
const sdVf = Math.sqrt(gals.reduce((a, g) => a + g.VfR ** 2, 0) / N);
const sdA0 = Math.sqrt(gals.reduce((a, g) => a + g.a0R ** 2, 0) / N);
for (let i = 0; i < N; i++) {
  gals[i].VfR_z = gals[i].VfR / sdVf;
  gals[i].a0R_z = gals[i].a0R / sdA0;
  gals[i].L_sum = gals[i].VfR_z + gals[i].a0R_z;
}
const ctrlVars = gals.map(g => [g.logK, g.dmFrac, g.envCode]);
const Lresid = multiR2(ctrlVars, gals.map(g => g.L_sum));
for (let i = 0; i < N; i++) gals[i].dq = Lresid.residuals[i];

const r_VfR_a0R_SPARC = pearsonR(gals.map(g => g.VfR), gals.map(g => g.a0R));
const r_DQ_hR_SPARC = pearsonR(gals.map(g => g.dq), gals.map(g => g.hR));

console.log('='.repeat(70));
console.log('PROGRAM 6B: DECISIVE GENERATIVE MODEL');
console.log('Make → Test → Iterate');
console.log('='.repeat(70));
console.log('\nN = ' + N + ' SPARC galaxies');
console.log('SPARC targets: r(VfR,a0R) = ' + r_VfR_a0R_SPARC.toFixed(3) + ', r(DQ,hR) = ' + r_DQ_hR_SPARC.toFixed(3));


function seedRng(seed) {
  let s = seed;
  return () => { s = (s * 16807 + 0) % 2147483647; return s / 2147483647; };
}

function gaussRng(rng) {
  let u, v, s2;
  do { u = 2 * rng() - 1; v = 2 * rng() - 1; s2 = u * u + v * v; } while (s2 >= 1 || s2 === 0);
  return u * Math.sqrt(-2 * Math.log(s2) / s2);
}

function generateModel(params, Nsim, seed) {
  const rng = seedRng(seed);
  const gs = [];
  for (let i = 0; i < Nsim; i++) {
    const logMbar = 8 + rng() * 3.5;
    const morphT = Math.floor(rng() * 10);
    const envCode = rng() < 0.3 ? 1 : 0;
    const logRdisk = -0.5 + rng() * 1.5;
    const logL36 = logMbar - 9 + gaussRng(rng) * 0.1;
    const logMHI = logMbar - 1 + gaussRng(rng) * 0.3;
    const logSBdisk = 1 + gaussRng(rng) * 0.3;

    const H = gaussRng(rng);

    let VfR = params.alpha_Vf * H + gaussRng(rng) * params.sigma_obs;
    let a0R = params.alpha_a0 * H + gaussRng(rng) * params.sigma_obs;

    let hR_base = 0.5 + gaussRng(rng) * 0.3;
    let dmFrac_base = 0.3 + rng() * 0.5;

    if (params.gamma_hR) {
      hR_base += params.gamma_hR * H;
    }

    if (params.gamma_dmFrac) {
      dmFrac_base += params.gamma_dmFrac * H;
    }

    let envEffect = 0;
    if (params.gamma_env) {
      envEffect = params.gamma_env * H * envCode;
      VfR += envEffect * 0.1;
    }

    let dq_component = 0;
    if (params.delta_dq) {
      dq_component = params.delta_dq * H;
    }

    const logK = -1 + gaussRng(rng) * 0.5;
    const outerSlope = gaussRng(rng) * 2;
    const rcSmooth = 0.95 + gaussRng(rng) * 0.02;
    const gasFrac = 0.1 + rng() * 0.5;

    gs.push({
      logMbar, morphT, envCode, logRdisk, logL36, logMHI, logSBdisk,
      logVflat: 0.25 * logMbar + gaussRng(rng) * 0.05,
      logA0: 3.5 + gaussRng(rng) * 0.1,
      H, VfR, a0R,
      hR: hR_base, dmFrac: dmFrac_base, logK, outerSlope, rcSmooth, gasFrac,
      dq_raw: dq_component,
    });
  }

  const struct4s = gs.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
  const struct6s = gs.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
  const vfMs = multiR2(struct4s, gs.map(g => g.logVflat));
  const a0Ms = multiR2(struct6s, gs.map(g => g.logA0));

  for (let i = 0; i < Nsim; i++) {
    gs[i].VfR_clean = vfMs.residuals[i] + gs[i].VfR;
    gs[i].a0R_clean = a0Ms.residuals[i] + gs[i].a0R;
  }

  const sdV = Math.sqrt(gs.reduce((a, g) => a + g.VfR_clean ** 2, 0) / Nsim);
  const sdA = Math.sqrt(gs.reduce((a, g) => a + g.a0R_clean ** 2, 0) / Nsim);
  for (let i = 0; i < Nsim; i++) {
    gs[i].VfR_z = gs[i].VfR_clean / (sdV || 1);
    gs[i].a0R_z = gs[i].a0R_clean / (sdA || 1);
    gs[i].L_sum = gs[i].VfR_z + gs[i].a0R_z;
  }
  const ctrlS = gs.map(g => [g.logK, g.dmFrac, g.envCode]);
  const LresS = multiR2(ctrlS, gs.map(g => g.L_sum));
  for (let i = 0; i < Nsim; i++) gs[i].dq = LresS.residuals[i];

  return gs;
}


function testModel(gs, params) {
  const Nsim = gs.length;
  const results = {};

  const r_vf_a0 = pearsonR(gs.map(g => g.VfR_clean), gs.map(g => g.a0R_clean));
  results.T1_bilateral = {
    r: r_vf_a0,
    target: 0.804,
    pass: r_vf_a0 > 0.3 && r_vf_a0 > 0,
    strong: Math.abs(r_vf_a0 - 0.804) < 0.15,
  };

  const r_dq_hr = pearsonR(gs.map(g => g.dq), gs.map(g => g.hR));
  const dqArr = gs.map(g => g.dq);
  const hrArr = gs.map(g => g.hR);
  const sorted = gs.slice().sort((a, b) => b.dq - a.dq);
  const topQ = sorted.slice(0, Math.floor(Nsim / 4));
  const botQ = sorted.slice(Math.floor(3 * Nsim / 4));
  const topHR = topQ.reduce((s, g) => s + g.hR, 0) / topQ.length;
  const botHR = botQ.reduce((s, g) => s + g.hR, 0) / botQ.length;
  results.T2_haloRespSign = {
    r_dq_hr,
    topQ_hR: topHR,
    botQ_hR: botHR,
    signPositive: topHR > botHR,
    pass: topHR > botHR,
  };

  const dqSD = Math.sqrt(gs.reduce((a, g) => a + g.dq ** 2, 0) / Nsim);
  const dqSkew = gs.reduce((a, g) => a + (g.dq / (dqSD || 1)) ** 3, 0) / Nsim;
  results.T3_dqDistribution = {
    dqSD,
    dqSkew,
    bilateral: Math.abs(dqSkew) < 1.0,
    pass: dqSD > 0.3 && Math.abs(dqSkew) < 1.5,
  };

  const halfN = Math.floor(Nsim / 2);
  const gs1 = gs.slice(0, halfN), gs2 = gs.slice(halfN);
  const r1 = pearsonR(gs1.map(g => g.VfR_clean), gs1.map(g => g.a0R_clean));
  const r2 = pearsonR(gs2.map(g => g.VfR_clean), gs2.map(g => g.a0R_clean));
  const shrinkage = Math.abs(r1 - r2);
  results.T4_constructionIndep = {
    r_half1: r1, r_half2: r2,
    shrinkage,
    pass: shrinkage < 0.15,
  };

  const sortedDQ = gs.slice().sort((a, b) => b.dq - a.dq);
  const highH = sortedDQ.slice(0, Math.floor(Nsim / 5));
  const lowH = sortedDQ.slice(Math.floor(4 * Nsim / 5));
  const highVfR = highH.reduce((s, g) => s + g.VfR_clean, 0) / highH.length;
  const lowVfR = lowH.reduce((s, g) => s + g.VfR_clean, 0) / lowH.length;
  const higha0R = highH.reduce((s, g) => s + g.a0R_clean, 0) / highH.length;
  const lowa0R = lowH.reduce((s, g) => s + g.a0R_clean, 0) / lowH.length;
  const highSmooth = highH.reduce((s, g) => s + g.rcSmooth, 0) / highH.length;
  const lowSmooth = lowH.reduce((s, g) => s + g.rcSmooth, 0) / lowH.length;
  const bilateralPass = highVfR > lowVfR && higha0R > lowa0R;
  const quietnessNotRequired = Math.abs(highSmooth - lowSmooth) < 0.05;
  results.T5_matchedPair = {
    highVfR, lowVfR, higha0R, lowa0R,
    highSmooth, lowSmooth,
    bilateralPass,
    quietnessNotRequired,
    pass: bilateralPass,
  };

  const ctrlMed = gs.map(g => [g.hR, g.rcSmooth, g.gasFrac]);
  const pr_vf = partialR(gs.map(g => g.dq), gs.map(g => g.VfR_clean), ctrlMed);
  const pr_a0 = partialR(gs.map(g => g.dq), gs.map(g => g.a0R_clean), ctrlMed);
  const r_dq_vf_zero = pearsonR(gs.map(g => g.dq), gs.map(g => g.VfR_clean));
  const r_dq_a0_zero = pearsonR(gs.map(g => g.dq), gs.map(g => g.a0R_clean));
  const red_vf = 1 - Math.abs(pr_vf) / Math.abs(r_dq_vf_zero || 0.001);
  const red_a0 = 1 - Math.abs(pr_a0) / Math.abs(r_dq_a0_zero || 0.001);
  results.T6_noMediator = {
    r_dq_vf_zero, r_dq_a0_zero,
    partial_vf: pr_vf, partial_a0: pr_a0,
    reduction_vf: red_vf, reduction_a0: red_a0,
    directSurvives: Math.abs(pr_vf) > 0.05 || Math.abs(pr_a0) > 0.05,
    pass: red_vf < 0.5 && red_a0 < 0.5,
  };

  let totalPass = 0;
  let criticalPass = 0;
  const tests = ['T1_bilateral', 'T2_haloRespSign', 'T3_dqDistribution', 'T4_constructionIndep', 'T5_matchedPair', 'T6_noMediator'];
  const criticalTests = ['T2_haloRespSign', 'T5_matchedPair', 'T6_noMediator'];
  for (const t of tests) if (results[t].pass) totalPass++;
  for (const t of criticalTests) if (results[t].pass) criticalPass++;

  results.summary = {
    totalPass,
    criticalPass,
    verdict: totalPass >= 5 && criticalPass >= 3 ? 'LEAD MODEL' :
             totalPass >= 4 ? 'PROMISING' :
             criticalPass < 2 ? 'KILLED' : 'WEAK',
  };

  return results;
}


console.log('\n\n' + '#'.repeat(70));
console.log('STAGE 1: MAKE — Three minimal models');
console.log('#'.repeat(70));

const models = {
  M1: {
    name: 'M1: Pure Common-Cause',
    desc: 'H → VfResid, H → a0Resid, nothing else',
    params: { alpha_Vf: 0.125, alpha_a0: 0.10, sigma_obs: 0.02 },
  },
  M2: {
    name: 'M2: Common-Cause + Halo Coupling',
    desc: 'M1 + weak H → haloResponse, H → dmFrac (downstream, not mediating)',
    params: { alpha_Vf: 0.125, alpha_a0: 0.10, sigma_obs: 0.02, gamma_hR: 0.08, gamma_dmFrac: 0.05, gamma_env: 0.03 },
  },
  M3: {
    name: 'M3: Common-Cause + Dark-Quarter Term',
    desc: 'M2 + independent DQ component from H',
    params: { alpha_Vf: 0.125, alpha_a0: 0.10, sigma_obs: 0.02, gamma_hR: 0.08, gamma_dmFrac: 0.05, gamma_env: 0.03, delta_dq: 0.15 },
  },
};

for (const [key, model] of Object.entries(models)) {
  console.log('\n  ' + model.name);
  console.log('  ' + model.desc);
  const pStr = Object.entries(model.params).map(([k, v]) => k + '=' + v).join(', ');
  console.log('  Params: ' + pStr);
}


console.log('\n\n' + '#'.repeat(70));
console.log('STAGE 2: TEST — Iteration 1 (minimal versions)');
console.log('#'.repeat(70));

const Nsim = 5000;
const seed = 42;
const iter1Results = {};

for (const [key, model] of Object.entries(models)) {
  const gs = generateModel(model.params, Nsim, seed);
  const results = testModel(gs, model.params);
  iter1Results[key] = results;

  console.log('\n  ' + '-'.repeat(65));
  console.log('  ' + model.name);
  console.log('  ' + '-'.repeat(65));
  const tests = ['T1_bilateral', 'T2_haloRespSign', 'T3_dqDistribution', 'T4_constructionIndep', 'T5_matchedPair', 'T6_noMediator'];
  const labels = ['T1: Bilateral channel', 'T2: haloResp sign (+)', 'T3: DQ distribution', 'T4: Construction-indep', 'T5: Matched-pair pattern', 'T6: No mediator needed'];
  for (let i = 0; i < tests.length; i++) {
    const r = results[tests[i]];
    const critical = ['T2_haloRespSign', 'T5_matchedPair', 'T6_noMediator'].includes(tests[i]);
    let detail = '';
    if (tests[i] === 'T1_bilateral') detail = 'r=' + r.r.toFixed(3) + ' (target: 0.804)';
    if (tests[i] === 'T2_haloRespSign') detail = 'r(DQ,hR)=' + r.r_dq_hr.toFixed(3) + ', top-bot=' + (r.topQ_hR - r.botQ_hR).toFixed(3);
    if (tests[i] === 'T3_dqDistribution') detail = 'SD=' + r.dqSD.toFixed(3) + ', skew=' + r.dqSkew.toFixed(3);
    if (tests[i] === 'T4_constructionIndep') detail = 'r1=' + r.r_half1.toFixed(3) + ', r2=' + r.r_half2.toFixed(3) + ', shrink=' + r.shrinkage.toFixed(3);
    if (tests[i] === 'T5_matchedPair') detail = 'VfR: ' + r.highVfR.toFixed(3) + '>' + r.lowVfR.toFixed(3) + ', a0R: ' + r.higha0R.toFixed(3) + '>' + r.lowa0R.toFixed(3);
    if (tests[i] === 'T6_noMediator') detail = 'red_vf=' + (r.reduction_vf * 100).toFixed(1) + '%, red_a0=' + (r.reduction_a0 * 100).toFixed(1) + '%';
    console.log('  ' + (r.pass ? 'PASS' : 'FAIL').padEnd(6) + (critical ? '*' : ' ') + labels[i].padEnd(28) + detail);
  }
  console.log('  SCORE: ' + results.summary.totalPass + '/6 (critical: ' + results.summary.criticalPass + '/3)  →  ' + results.summary.verdict);
}


console.log('\n\n' + '#'.repeat(70));
console.log('STAGE 3: ITERATE — Tune the best model');
console.log('#'.repeat(70));

const bestKey = Object.entries(iter1Results).sort((a, b) => {
  const sa = a[1].summary; const sb = b[1].summary;
  if (sa.totalPass !== sb.totalPass) return sb.totalPass - sa.totalPass;
  return sb.criticalPass - sa.criticalPass;
})[0][0];

console.log('\n  Best from Iteration 1: ' + bestKey + ' (' + models[bestKey].name + ')');
console.log('  Score: ' + iter1Results[bestKey].summary.totalPass + '/6');
console.log('\n  Tuning only 3 knobs: alpha_Vf, alpha_a0, delta_dq');

const baseParams = Object.assign({}, models[bestKey].params);

const tuningGrid = [
  Object.assign({}, baseParams, { name: 'baseline' }),
  Object.assign({}, baseParams, { name: 'alpha x1.5', alpha_Vf: baseParams.alpha_Vf * 1.5, alpha_a0: baseParams.alpha_a0 * 1.5 }),
  Object.assign({}, baseParams, { name: 'alpha x2.0', alpha_Vf: baseParams.alpha_Vf * 2.0, alpha_a0: baseParams.alpha_a0 * 2.0 }),
  Object.assign({}, baseParams, { name: 'alpha x2.5', alpha_Vf: baseParams.alpha_Vf * 2.5, alpha_a0: baseParams.alpha_a0 * 2.5 }),
  Object.assign({}, baseParams, { name: 'alpha x3.0', alpha_Vf: baseParams.alpha_Vf * 3.0, alpha_a0: baseParams.alpha_a0 * 3.0 }),
  Object.assign({}, baseParams, { name: 'dq x0', delta_dq: 0 }),
  Object.assign({}, baseParams, { name: 'dq x2', delta_dq: (baseParams.delta_dq || 0.15) * 2 }),
  Object.assign({}, baseParams, { name: 'alpha x2.5 + dq x2', alpha_Vf: baseParams.alpha_Vf * 2.5, alpha_a0: baseParams.alpha_a0 * 2.5, delta_dq: (baseParams.delta_dq || 0.15) * 2 }),
];

console.log('\n  ITERATION 2 — AMPLITUDE TUNING:');
console.log('  ' + '-'.repeat(100));
console.log('  ' + 'Variant'.padEnd(25) + 'r(VfR,a0R)'.padEnd(12) + 'hR sign'.padEnd(10) + 'DQ SD'.padEnd(10) + 'CI shrink'.padEnd(12) + 'bilateral'.padEnd(12) + 'no-med'.padEnd(10) + 'Score'.padEnd(8) + 'Verdict');
console.log('  ' + '-'.repeat(100));

let bestTune = null;
let bestTuneScore = 0;
const iter2Results = [];

for (const tune of tuningGrid) {
  const gs = generateModel(tune, Nsim, seed);
  const results = testModel(gs, tune);
  const score = results.summary.totalPass;
  iter2Results.push({ name: tune.name, results, params: tune });

  if (score > bestTuneScore || (score === bestTuneScore && results.summary.criticalPass > (bestTune ? bestTune.results.summary.criticalPass : 0))) {
    bestTuneScore = score;
    bestTune = { name: tune.name, results, params: tune };
  }

  console.log('  ' + tune.name.padEnd(25) +
    results.T1_bilateral.r.toFixed(3).padEnd(12) +
    (results.T2_haloRespSign.pass ? '+' : '-').padEnd(10) +
    results.T3_dqDistribution.dqSD.toFixed(3).padEnd(10) +
    results.T4_constructionIndep.shrinkage.toFixed(3).padEnd(12) +
    (results.T5_matchedPair.pass ? 'YES' : 'no').padEnd(12) +
    (results.T6_noMediator.pass ? 'YES' : 'no').padEnd(10) +
    (score + '/6').padEnd(8) +
    results.summary.verdict);
}


console.log('\n\n' + '#'.repeat(70));
console.log('ITERATION 3: FINAL VALIDATION OF BEST MODEL');
console.log('#'.repeat(70));

console.log('\n  Best tuned variant: ' + bestTune.name);
const finalParams = bestTune.params;
console.log('  Params: ' + Object.entries(finalParams).filter(([k]) => k !== 'name').map(([k, v]) => k + '=' + v).join(', '));

const seeds = [42, 137, 271, 314, 577];
console.log('\n  MULTI-SEED STABILITY TEST:');
console.log('  ' + '-'.repeat(80));
console.log('  ' + 'Seed'.padEnd(8) + 'r(VfR,a0R)'.padEnd(12) + 'hR+'.padEnd(6) + 'DQ_SD'.padEnd(10) + 'CI'.padEnd(10) + 'bilateral'.padEnd(12) + 'no-med'.padEnd(10) + 'Score');
console.log('  ' + '-'.repeat(80));

let allScores = [];
for (const s of seeds) {
  const gs = generateModel(finalParams, Nsim, s);
  const results = testModel(gs, finalParams);
  allScores.push(results.summary.totalPass);
  console.log('  ' + (s + '').padEnd(8) +
    results.T1_bilateral.r.toFixed(3).padEnd(12) +
    (results.T2_haloRespSign.pass ? '+' : '-').padEnd(6) +
    results.T3_dqDistribution.dqSD.toFixed(3).padEnd(10) +
    results.T4_constructionIndep.shrinkage.toFixed(3).padEnd(10) +
    (results.T5_matchedPair.pass ? 'YES' : 'no').padEnd(12) +
    (results.T6_noMediator.pass ? 'YES' : 'no').padEnd(10) +
    results.summary.totalPass + '/6');
}

const meanScore = allScores.reduce((a, b) => a + b, 0) / allScores.length;
const minScore = Math.min(...allScores);
console.log('\n  Mean score: ' + meanScore.toFixed(1) + '/6');
console.log('  Min score:  ' + minScore + '/6');
console.log('  Stable (all >= 4): ' + (minScore >= 4 ? 'YES' : 'NO'));


console.log('\n\n' + '#'.repeat(70));
console.log('COMPARISON: WHY OTHER MODELS FAILED');
console.log('#'.repeat(70));

for (const [key, model] of Object.entries(models)) {
  if (key === bestKey) continue;
  const r = iter1Results[key];
  const failures = [];
  if (!r.T1_bilateral.pass) failures.push('T1: bilateral channel too weak (r=' + r.T1_bilateral.r.toFixed(3) + ')');
  if (!r.T2_haloRespSign.pass) failures.push('T2: haloResp sign wrong');
  if (!r.T3_dqDistribution.pass) failures.push('T3: DQ distribution bad');
  if (!r.T4_constructionIndep.pass) failures.push('T4: construction-dependent');
  if (!r.T5_matchedPair.pass) failures.push('T5: matched-pair pattern wrong');
  if (!r.T6_noMediator.pass) failures.push('T6: needs mediator');
  console.log('\n  ' + model.name + ': ' + r.summary.totalPass + '/6 → ' + r.summary.verdict);
  if (failures.length > 0) {
    for (const f of failures) console.log('    FAILED: ' + f);
  } else {
    console.log('    (all tests passed, but outperformed by ' + bestKey + ')');
  }
}


console.log('\n\n' + '='.repeat(70));
console.log('PROGRAM 6B GRAND VERDICT');
console.log('='.repeat(70));

console.log('\n  LEAD MODEL: ' + bestTune.name + ' (based on ' + bestKey + ')');
console.log('  Final score: ' + bestTuneScore + '/6');

const finalGS = generateModel(finalParams, Nsim, 42);
const finalResults = testModel(finalGS, finalParams);

console.log('\n  FULL TEST RESULTS:');
const fTests = ['T1_bilateral', 'T2_haloRespSign', 'T3_dqDistribution', 'T4_constructionIndep', 'T5_matchedPair', 'T6_noMediator'];
const fLabels = ['T1: Bilateral channel', 'T2: haloResp sign (+)', 'T3: DQ distribution', 'T4: Construction-indep', 'T5: Matched-pair pattern', 'T6: No mediator needed'];
for (let i = 0; i < fTests.length; i++) {
  const r = finalResults[fTests[i]];
  console.log('  ' + (r.pass ? 'PASS' : 'FAIL').padEnd(6) + fLabels[i]);
}

console.log('\n  KEY METRICS:');
console.log('    r(VfR,a0R) = ' + finalResults.T1_bilateral.r.toFixed(3) + ' (SPARC: ' + r_VfR_a0R_SPARC.toFixed(3) + ')');
console.log('    haloResp sign: ' + (finalResults.T2_haloRespSign.signPositive ? 'POSITIVE (correct)' : 'NEGATIVE (wrong)'));
console.log('    DQ SD = ' + finalResults.T3_dqDistribution.dqSD.toFixed(3));
console.log('    Construction shrinkage = ' + finalResults.T4_constructionIndep.shrinkage.toFixed(3));
console.log('    Bilateral matched: ' + (finalResults.T5_matchedPair.bilateralPass ? 'YES' : 'NO'));
console.log('    Direct path survives: ' + (finalResults.T6_noMediator.directSurvives ? 'YES' : 'NO'));

console.log('\n  DECISIVE NEW PREDICTION:');
console.log('  If H is a genuine common-cause hidden state with bilateral drive:');
console.log('  1. IFU surveys should find that galaxies with elevated VfResid AND a0Resid');
console.log('     simultaneously have systematically DIFFERENT halo density profiles');
console.log('     (higher inner density at fixed concentration) than matched controls.');
console.log('  2. Cosmological simulations (NIHAO/FIRE) should produce a latent variable');
console.log('     in halo assembly history that correlates with both Vflat excess and');
console.log('     acceleration-scale excess, but NOT with kinematic morphology metrics.');
console.log('  3. The bilateral residual excess should be ABSENT in isolated WDM halos');
console.log('     but PRESENT in CDM halos with standard assembly histories.');


const outPath = path.join(__dirname, '..', 'public', 'program6b-decisive-model.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: '6B',
  title: 'Decisive Generative Model',
  timestamp: new Date().toISOString(),
  N_sparc: N, N_sim: Nsim,
  sparc_targets: { r_VfR_a0R: r_VfR_a0R_SPARC, r_DQ_hR: r_DQ_hR_SPARC },
  models: Object.fromEntries(Object.entries(models).map(([k, v]) => [k, { name: v.name, params: v.params, iter1: iter1Results[k].summary }])),
  bestModel: bestKey,
  bestTune: { name: bestTune.name, params: Object.fromEntries(Object.entries(finalParams).filter(([k]) => k !== 'name')), score: bestTuneScore },
  finalResults: {
    T1: finalResults.T1_bilateral,
    T2: finalResults.T2_haloRespSign,
    T3: finalResults.T3_dqDistribution,
    T4: finalResults.T4_constructionIndep,
    T5: finalResults.T5_matchedPair,
    T6: finalResults.T6_noMediator,
    summary: finalResults.summary,
  },
  multiSeedStability: { scores: allScores, meanScore, minScore },
}, null, 2));
console.log('\nSaved: ' + outPath);
