const fs = require('fs');
const path = require('path');

const G = 4.3009e-6;
const rhoCrit = 136.18;

function nfwMenc(r, M200, c) {
  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const rs = r200 / c;
  const x = r / rs;
  const norm = Math.log(1 + c) - c / (1 + c);
  return M200 * (Math.log(1 + x) - x / (1 + x)) / norm;
}

function coredNFWMenc(r, M200, c, rCore, alpha) {
  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const rs = r200 / c;
  const x = r / rs;
  const norm = Math.log(1 + c) - c / (1 + c);
  const coreSuppress = Math.pow(1 + Math.pow(rCore / r, 2), -alpha / 2);
  return M200 * (Math.log(1 + x) - x / (1 + x)) / norm * coreSuppress;
}

function contractedNFWMenc(r, M200, c, Mbar, Rd) {
  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const fb = Mbar / M200;
  const rContract = r * (1 - fb * Math.exp(-r / (2 * Rd)));
  const rEff = Math.max(rContract, r * 0.3);
  const rs = r200 / c;
  const x = rEff / rs;
  const norm = Math.log(1 + c) - c / (1 + c);
  return M200 * (1 - fb) * (Math.log(1 + x) - x / (1 + x)) / norm;
}

function expDiskMenc(r, Mdisk, Rd) {
  const y = r / Rd;
  return Mdisk * (1 - (1 + y) * Math.exp(-y));
}

function gasDiskMenc(r, Mgas, Rg) {
  const y = r / Rg;
  return Mgas * (1 - (1 + y) * Math.exp(-y));
}

function rarPredict(gbar, a0) {
  return gbar / (1 - Math.exp(-Math.sqrt(gbar / a0)));
}

function fitA0(points) {
  if (points.length < 3) return { a0: 3703, rmse: 1 };
  let bestA0 = 3703, bestRMSE = Infinity;
  for (let logA0 = 2.5; logA0 <= 4.5; logA0 += 0.01) {
    const a0 = Math.pow(10, logA0);
    let sse = 0;
    for (const p of points) {
      const pred = rarPredict(p.gbar, a0);
      const diff = Math.log10(p.gobs) - Math.log10(pred);
      sse += diff * diff;
    }
    const rmse = Math.sqrt(sse / points.length);
    if (rmse < bestRMSE) { bestRMSE = rmse; bestA0 = a0; }
  }
  return { a0: bestA0, rmse: bestRMSE };
}

function pearsonR(x, y) {
  const n = x.length;
  if (n < 3) return 0;
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) {
    const dx = x[i] - mx, dy = y[i] - my;
    num += dx * dy; dx2 += dx * dx; dy2 += dy * dy;
  }
  return dx2 > 0 && dy2 > 0 ? num / Math.sqrt(dx2 * dy2) : 0;
}

function linReg(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let num = 0, den = 0;
  for (let i = 0; i < n; i++) { num += (x[i] - mx) * (y[i] - my); den += (x[i] - mx) ** 2; }
  const slope = den > 0 ? num / den : 0;
  const intercept = my - slope * mx;
  return { slope, intercept };
}

function rand(min, max) { return min + Math.random() * (max - min); }
function randNorm(mu, sigma) {
  const u1 = Math.random(), u2 = Math.random();
  return mu + sigma * Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}

console.log('='.repeat(70));
console.log('PHASE 400: Origin of VfResid');
console.log('What physically generates the dominant a₀-coupling channel?');
console.log('='.repeat(70));

const NGAL = 300;
const A0_TRUE = 3703;

function generateGalaxy(logM200, scenario) {
  const M200 = Math.pow(10, logM200);

  const cBase = 10 * Math.pow(M200 / 1e12, -0.1);
  const cScatter = randNorm(0, 0.12);
  const c = cBase * Math.pow(10, cScatter);

  const fbBase = 0.05 * Math.pow(M200 / 1e12, 0.3);
  const fb = Math.min(fbBase * (1 + randNorm(0, 0.2)), 0.16);
  const Mbar = M200 * fb;
  const fgas = 0.1 + 0.5 * Math.pow(M200 / 1e10, -0.3) + randNorm(0, 0.05);
  const fgasClamped = Math.max(0.05, Math.min(0.9, fgas));
  const Mstar = Mbar * (1 - fgasClamped);
  const Mgas = Mbar * fgasClamped;

  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const lambda = 0.035 * (1 + randNorm(0, 0.3));
  const Rd = lambda * r200 / Math.sqrt(2);
  const Rg = Rd * (1.5 + rand(0, 1.5));

  let dynTimeProxy = 0;
  let assemblyProxy = 0;
  let feedbackStrength = 0;
  let contractionStrength = 0;

  switch (scenario) {
    case 'baseline':
      break;

    case 'halo_response':
      contractionStrength = 0.3 + 0.5 * (logM200 - 10) / 2.5;
      contractionStrength = Math.max(0, Math.min(1, contractionStrength));
      break;

    case 'assembly_history':
      assemblyProxy = randNorm(0, 0.3);
      break;

    case 'feedback':
      const mstarMhalo = Mstar / M200;
      if (mstarMhalo < 0.001) feedbackStrength = 0.8;
      else if (mstarMhalo < 0.01) feedbackStrength = 0.5;
      else feedbackStrength = 0.1;
      break;

    case 'dynamical_integration':
      dynTimeProxy = (logM200 - 10) / 2.5;
      dynTimeProxy = Math.max(0, Math.min(1, dynTimeProxy));
      break;

    case 'combined':
      contractionStrength = 0.2 + 0.3 * (logM200 - 10) / 2.5;
      contractionStrength = Math.max(0, Math.min(0.5, contractionStrength));
      const msr = Mstar / M200;
      if (msr < 0.001) feedbackStrength = 0.6;
      else if (msr < 0.01) feedbackStrength = 0.3;
      else feedbackStrength = 0.05;
      dynTimeProxy = (logM200 - 10) / 2.5;
      dynTimeProxy = Math.max(0, Math.min(1, dynTimeProxy));
      assemblyProxy = randNorm(0, 0.15);
      break;
  }

  let a0Effective = A0_TRUE;
  if (scenario === 'assembly_history') {
    a0Effective = A0_TRUE * Math.pow(10, 0.15 * assemblyProxy);
  }
  if (scenario === 'dynamical_integration') {
    a0Effective = A0_TRUE * Math.pow(10, 0.2 * dynTimeProxy + randNorm(0, 0.05));
  }
  if (scenario === 'combined') {
    const shift = 0.12 * dynTimeProxy + 0.05 * contractionStrength - 0.08 * feedbackStrength + 0.05 * assemblyProxy;
    a0Effective = A0_TRUE * Math.pow(10, shift + randNorm(0, 0.04));
  }

  const nPts = 20;
  const rMin = 0.3 * Rd;
  const rMax = 5 * Rd;
  const points = [];
  let Vflat = 0;

  for (let i = 0; i < nPts; i++) {
    const r = rMin + (rMax - rMin) * (i / (nPts - 1));

    const MstarEnc = expDiskMenc(r, Mstar, Rd);
    const MgasEnc = gasDiskMenc(r, Mgas * 1.33, Rg);
    const MbarEnc = MstarEnc + MgasEnc;

    let MdmEnc;
    if (scenario === 'halo_response' || (scenario === 'combined' && contractionStrength > 0.1)) {
      MdmEnc = contractedNFWMenc(r, M200, c, Mbar, Rd) * (1 - contractionStrength * 0.3);
      MdmEnc += contractionStrength * nfwMenc(r, M200, c) * 0.3;
    } else if (scenario === 'feedback' || (scenario === 'combined' && feedbackStrength > 0.1)) {
      const rCore = feedbackStrength * Rd;
      MdmEnc = coredNFWMenc(r, M200, c, rCore, feedbackStrength * 2);
    } else {
      MdmEnc = nfwMenc(r, M200, c);
    }

    const Mtot = MbarEnc + MdmEnc;
    let Vobs = Math.sqrt(G * Mtot / r);
    const Vbar = Math.sqrt(G * MbarEnc / r);

    const gbar = Vbar * Vbar / r;

    if (a0Effective !== A0_TRUE && gbar > 0) {
      const gobsRarTrue = rarPredict(gbar, A0_TRUE);
      const gobsRarShifted = rarPredict(gbar, a0Effective);
      if (gobsRarTrue > 0) {
        const correction = gobsRarShifted / gobsRarTrue;
        Vobs = Vobs * Math.sqrt(correction);
      }
    }

    if (Vobs > Vflat) Vflat = Vobs;

    const gobs = Vobs * Vobs / r;

    if (gbar > 0 && gobs > 0) {
      const noise = 1 + randNorm(0, 0.03);
      points.push({ r, gobs: gobs * noise, gbar, Vobs, Vbar });
    }
  }

  const L36 = Mstar / 0.5;
  const logMbar = Math.log10(Mbar);
  const logL36 = Math.log10(L36);
  const logRd = Math.log10(Rd);
  const morphT = 3 + (logM200 - 11) * 2 + randNorm(0, 1.5);

  const fit = fitA0(points);

  const haloK = Math.log10(M200 * fb > 0 ? (nfwMenc(2 * Rd, M200, c) / MbarEnc_at2Rd(Mstar, Mgas, Rd, Rg)) : 1);

  function MbarEnc_at2Rd(ms, mg, rd, rg) {
    return expDiskMenc(2*rd, ms, rd) + gasDiskMenc(2*rd, mg*1.33, rg);
  }

  return {
    logM200,
    c,
    fb,
    fgas: fgasClamped,
    Mstar,
    Mgas,
    Mbar,
    Rd,
    Rg,
    Vflat,
    logVflat: Math.log10(Math.max(Vflat, 1)),
    logMbar,
    logL36,
    logRd,
    morphT,
    logMHI: Math.log10(Mgas / 1e9),
    a0: fit.a0,
    logA0: Math.log10(fit.a0),
    a0True: a0Effective,
    rmse: fit.rmse,
    haloK,
    nPts: points.length,
    dynTimeProxy,
    assemblyProxy,
    feedbackStrength,
    contractionStrength,
    scenario
  };
}

function computeVfResid(galaxies) {
  const x = galaxies.map(g => [g.logMbar, g.logL36, g.logRd, g.morphT]);
  const y = galaxies.map(g => g.logVflat);

  const n = y.length;
  const nv = 4;
  const mx = Array(nv).fill(0);
  const my = y.reduce((a, b) => a + b, 0) / n;
  for (let j = 0; j < nv; j++) {
    for (let i = 0; i < n; i++) mx[j] += x[i][j];
    mx[j] /= n;
  }

  const XTX = Array.from({length: nv}, () => Array(nv).fill(0));
  const XTy = Array(nv).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < nv; j++) {
      XTy[j] += (x[i][j] - mx[j]) * (y[i] - my);
      for (let k = 0; k < nv; k++) {
        XTX[j][k] += (x[i][j] - mx[j]) * (x[i][k] - mx[k]);
      }
    }
  }

  const beta = solveLinear(XTX, XTy);

  let ssRes = 0, ssTot = 0;
  for (let i = 0; i < n; i++) {
    let pred = my;
    for (let j = 0; j < nv; j++) pred += beta[j] * (x[i][j] - mx[j]);
    galaxies[i].vfPred = pred;
    galaxies[i].VfResid = y[i] - pred;
    ssRes += (y[i] - pred) ** 2;
    ssTot += (y[i] - my) ** 2;
  }

  return { R2: 1 - ssRes / ssTot, beta, mx, my };
}

function solveLinear(A, b) {
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
  const x = Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = aug[i][n];
    for (let j = i + 1; j < n; j++) x[i] -= aug[i][j] * x[j];
    x[i] /= aug[i][i] || 1;
  }
  return x;
}

function analyzeScenario(name, galaxies) {
  console.log(`\n${'─'.repeat(60)}`);
  console.log(`Scenario: ${name} (N=${galaxies.length})`);
  console.log('─'.repeat(60));

  const valid = galaxies.filter(g => g.rmse < 0.3 && isFinite(g.logA0) && isFinite(g.logVflat) && g.Vflat > 20);
  console.log(`  Valid galaxies: ${valid.length}/${galaxies.length}`);

  const structFit = computeVfResid(valid);
  console.log(`  Structural Vflat model R²: ${structFit.R2.toFixed(3)}`);

  const logA0 = valid.map(g => g.logA0);
  const VfR = valid.map(g => g.VfResid);
  const logVf = valid.map(g => g.logVflat);

  const rVfResidA0 = pearsonR(VfR, logA0);
  const rVflatA0 = pearsonR(logVf, logA0);
  console.log(`  r(VfResid, logA0) = ${rVfResidA0.toFixed(3)}`);
  console.log(`  r(logVflat, logA0) = ${rVflatA0.toFixed(3)}`);

  const rHaloK = pearsonR(valid.map(g => g.haloK), VfR);
  const rLogM200 = pearsonR(valid.map(g => g.logM200), VfR);
  const rConc = pearsonR(valid.map(g => Math.log10(g.c)), VfR);
  const rFgas = pearsonR(valid.map(g => g.fgas), VfR);

  console.log(`\n  VfResid correlations with known halo truth:`);
  console.log(`    r(haloK, VfResid) = ${rHaloK.toFixed(3)}`);
  console.log(`    r(logM200, VfResid) = ${rLogM200.toFixed(3)}`);
  console.log(`    r(log c, VfResid) = ${rConc.toFixed(3)}`);
  console.log(`    r(fgas, VfResid) = ${rFgas.toFixed(3)}`);

  const highV = valid.filter(g => g.Vflat >= 120);
  const lowV = valid.filter(g => g.Vflat < 120);
  const vhighV = valid.filter(g => g.Vflat >= 180);

  console.log(`\n  Regime analysis:`);
  if (lowV.length >= 5) {
    const rLow = pearsonR(lowV.map(g => g.VfResid), lowV.map(g => g.logA0));
    console.log(`    Low-V (N=${lowV.length}): r(VfResid, logA0) = ${rLow.toFixed(3)}`);
  }
  if (highV.length >= 5) {
    const rHigh = pearsonR(highV.map(g => g.VfResid), highV.map(g => g.logA0));
    console.log(`    High-V (N=${highV.length}): r(VfResid, logA0) = ${rHigh.toFixed(3)}`);
  }
  if (vhighV.length >= 5) {
    const rVH = pearsonR(vhighV.map(g => g.VfResid), vhighV.map(g => g.logA0));
    console.log(`    V-High-V (N=${vhighV.length}): r(VfResid, logA0) = ${rVH.toFixed(3)}`);
  }

  const trueParams = {};
  if (name.includes('response')) {
    const rContract = pearsonR(valid.map(g => g.contractionStrength), VfR);
    console.log(`    r(contractionStrength, VfResid) = ${rContract.toFixed(3)}`);
    trueParams.contractionStrength = rContract;
  }
  if (name.includes('feedback')) {
    const rFB = pearsonR(valid.map(g => g.feedbackStrength), VfR);
    console.log(`    r(feedbackStrength, VfResid) = ${rFB.toFixed(3)}`);
    trueParams.feedbackStrength = rFB;
  }
  if (name.includes('integration')) {
    const rDyn = pearsonR(valid.map(g => g.dynTimeProxy), VfR);
    console.log(`    r(dynTimeProxy, VfResid) = ${rDyn.toFixed(3)}`);
    trueParams.dynTimeProxy = rDyn;
  }
  if (name.includes('assembly')) {
    const rAsm = pearsonR(valid.map(g => g.assemblyProxy), VfR);
    console.log(`    r(assemblyProxy, VfResid) = ${rAsm.toFixed(3)}`);
    trueParams.assemblyProxy = rAsm;
  }

  const rDynTime = pearsonR(valid.map(g => g.dynTimeProxy), VfR);
  const rFeedback = pearsonR(valid.map(g => g.feedbackStrength), VfR);
  const rAssembly = pearsonR(valid.map(g => g.assemblyProxy), VfR);
  const rContraction = pearsonR(valid.map(g => g.contractionStrength), VfR);

  console.log(`\n  All mechanism correlations with VfResid:`);
  console.log(`    dynTimeProxy:       ${rDynTime.toFixed(3)}`);
  console.log(`    contractionStr:     ${rContraction.toFixed(3)}`);
  console.log(`    feedbackStrength:   ${rFeedback.toFixed(3)}`);
  console.log(`    assemblyProxy:      ${rAssembly.toFixed(3)}`);

  const rHaloKA0 = pearsonR(valid.map(g => g.haloK), logA0);
  const rLogM200A0 = pearsonR(valid.map(g => g.logM200), logA0);
  const rConcA0 = pearsonR(valid.map(g => Math.log10(g.c)), logA0);

  console.log(`\n  True halo parameters → a₀:`);
  console.log(`    r(haloK, logA0) = ${rHaloKA0.toFixed(3)}`);
  console.log(`    r(logM200, logA0) = ${rLogM200A0.toFixed(3)}`);
  console.log(`    r(log c, logA0) = ${rConcA0.toFixed(3)}`);

  let structA0R2 = -1;
  let structPlusVfR_A0R2 = -1;
  {
    const sx = valid.map(g => [g.logMbar, g.logL36, g.logRd, g.morphT]);
    const sy = logA0;
    const sn = sy.length;
    const smx = Array(4).fill(0);
    const smy = sy.reduce((a,b) => a+b, 0) / sn;
    for (let j = 0; j < 4; j++) {
      for (let i = 0; i < sn; i++) smx[j] += sx[i][j];
      smx[j] /= sn;
    }
    const sXTX = Array.from({length:4}, ()=>Array(4).fill(0));
    const sXTy = Array(4).fill(0);
    for (let i = 0; i < sn; i++) {
      for (let j = 0; j < 4; j++) {
        sXTy[j] += (sx[i][j] - smx[j]) * (sy[i] - smy);
        for (let k = 0; k < 4; k++) sXTX[j][k] += (sx[i][j] - smx[j]) * (sx[i][k] - smx[k]);
      }
    }
    const sb = solveLinear(sXTX, sXTy);
    let ssr = 0, sst = 0;
    for (let i = 0; i < sn; i++) {
      let p = smy;
      for (let j = 0; j < 4; j++) p += sb[j] * (sx[i][j] - smx[j]);
      ssr += (sy[i] - p) ** 2;
      sst += (sy[i] - smy) ** 2;
    }
    structA0R2 = sst > 0 ? 1 - ssr / sst : 0;

    const sx5 = valid.map((g, i) => [...sx[i], VfR[i]]);
    const smx5 = Array(5).fill(0);
    for (let j = 0; j < 5; j++) {
      for (let i = 0; i < sn; i++) smx5[j] += sx5[i][j];
      smx5[j] /= sn;
    }
    const sXTX5 = Array.from({length:5}, ()=>Array(5).fill(0));
    const sXTy5 = Array(5).fill(0);
    for (let i = 0; i < sn; i++) {
      for (let j = 0; j < 5; j++) {
        sXTy5[j] += (sx5[i][j] - smx5[j]) * (sy[i] - smy);
        for (let k = 0; k < 5; k++) sXTX5[j][k] += (sx5[i][j] - smx5[j]) * (sx5[i][k] - smx5[k]);
      }
    }
    const sb5 = solveLinear(sXTX5, sXTy5);
    let ssr5 = 0;
    for (let i = 0; i < sn; i++) {
      let p = smy;
      for (let j = 0; j < 5; j++) p += sb5[j] * (sx5[i][j] - smx5[j]);
      ssr5 += (sy[i] - p) ** 2;
    }
    structPlusVfR_A0R2 = sst > 0 ? 1 - ssr5 / sst : 0;
  }
  const vfResidLift = structPlusVfR_A0R2 - structA0R2;
  console.log(`\n  Structural → logA0: R² = ${structA0R2.toFixed(3)}`);
  console.log(`  Structural+VfResid → logA0: R² = ${structPlusVfR_A0R2.toFixed(3)}`);
  console.log(`  VfResid lift: +${(vfResidLift * 100).toFixed(1)}pp (non-structural information in VfResid)`);

  return {
    name,
    N: valid.length,
    structR2: structFit.R2,
    rVfResidA0,
    rVflatA0,
    rHaloK,
    rLogM200,
    rConc,
    rFgas,
    rHaloKA0,
    rLogM200A0,
    rConcA0,
    structA0R2,
    structPlusVfR_A0R2,
    vfResidLift,
    regimeLow: lowV.length >= 5 ? pearsonR(lowV.map(g=>g.VfResid), lowV.map(g=>g.logA0)) : null,
    regimeHigh: highV.length >= 5 ? pearsonR(highV.map(g=>g.VfResid), highV.map(g=>g.logA0)) : null,
    regimeVHigh: vhighV.length >= 5 ? pearsonR(vhighV.map(g=>g.VfResid), vhighV.map(g=>g.logA0)) : null,
    trueParams,
    mechanismCorrs: { rDynTime, rFeedback, rAssembly, rContraction }
  };
}

const scenarios = [
  { name: 'baseline', desc: 'Pure NFW + exponential disk, universal a₀' },
  { name: 'halo_response', desc: 'Adiabatic contraction (mass-dependent)' },
  { name: 'assembly_history', desc: 'Assembly-dependent a₀ scatter' },
  { name: 'feedback', desc: 'SN feedback core formation (Di Cintio)' },
  { name: 'dynamical_integration', desc: 'Accumulated baryon-halo processing' },
  { name: 'combined', desc: 'All mechanisms active simultaneously' }
];

const allResults = [];

for (const sc of scenarios) {
  console.log(`\n${'═'.repeat(70)}`);
  console.log(`GENERATING: ${sc.name} — ${sc.desc}`);

  Math.seedrandom = null;
  const galaxies = [];
  for (let i = 0; i < NGAL; i++) {
    const logM200 = rand(10.0, 12.5);
    galaxies.push(generateGalaxy(logM200, sc.name));
  }

  const result = analyzeScenario(sc.name, galaxies);
  result.desc = sc.desc;
  allResults.push(result);
}

console.log('\n\n');
console.log('═'.repeat(70));
console.log('PHASE 400 SUMMARY: Which scenario reproduces the SPARC VfResid pattern?');
console.log('═'.repeat(70));

console.log('\nSPARC observed pattern (for comparison):');
console.log('  r(VfResid, logA0) = 0.70 (full), 0.78 (high-V), 0.83 (V-high-V)');
console.log('  r(haloK, VfResid) = 0.60');
console.log('  Structural → VfResid: R² < 0 (FAILS)');
console.log('  Regime strengthening: low→high r increases');

console.log('\n' + '─'.repeat(70));
console.log(padR('Scenario', 25) + padR('r(VfR,A0)', 12) + padR('r(hK,VfR)', 12) +
            padR('VfR lift', 10) + padR('Low-V', 10) + padR('Hi-V', 10) + 'VHi-V');
console.log('─'.repeat(70));

function padR(s, n) { return String(s).padEnd(n); }

for (const r of allResults) {
  console.log(
    padR(r.name, 25) +
    padR(r.rVfResidA0.toFixed(3), 12) +
    padR(r.rHaloK.toFixed(3), 12) +
    padR((r.vfResidLift * 100).toFixed(1) + 'pp', 10) +
    padR(r.regimeLow !== null ? r.regimeLow.toFixed(3) : 'N/A', 10) +
    padR(r.regimeHigh !== null ? r.regimeHigh.toFixed(3) : 'N/A', 10) +
    (r.regimeVHigh !== null ? r.regimeVHigh.toFixed(3) : 'N/A')
  );
}

console.log('─'.repeat(70));

console.log('\nDiagnostic: Which scenario best matches SPARC observations?');
const sparcPattern = { rVfRA0: 0.70, rHaloK: 0.60, structR2: -0.14, regimeStrength: true };

for (const r of allResults) {
  let score = 0;
  let notes = [];

  if (Math.abs(r.rVfResidA0) > 0.3) { score++; notes.push('VfR-A0 signal present'); }
  if (Math.abs(r.rVfResidA0) > 0.5) { score++; notes.push('VfR-A0 signal strong'); }
  if (Math.abs(r.rHaloK) > 0.3) { score++; notes.push('haloK-VfR link'); }
  if (r.vfResidLift > 0.3) { score++; notes.push('VfResid carries non-structural info'); }
  if (r.regimeLow !== null && r.regimeHigh !== null) {
    if (Math.abs(r.regimeHigh) > Math.abs(r.regimeLow) + 0.1) {
      score++; notes.push('regime strengthening');
    }
  }
  if (Math.abs(r.rVfResidA0 - 0.70) < 0.2) { score++; notes.push('magnitude match'); }

  console.log(`  ${r.name}: score=${score}/6 [${notes.join(', ')}]`);
}

console.log('\n' + '═'.repeat(70));
console.log('KEY QUESTION ANSWERS:');
console.log('═'.repeat(70));

const baselineResult = allResults.find(r => r.name === 'baseline');
const dynResult = allResults.find(r => r.name === 'dynamical_integration');
const combinedResult = allResults.find(r => r.name === 'combined');
const feedbackResult = allResults.find(r => r.name === 'feedback');
const responseResult = allResults.find(r => r.name === 'halo_response');

console.log('\n1. Does VfResid appear in pure ΛCDM (no extra physics)?');
console.log(`   r(VfResid, logA0) = ${baselineResult.rVfResidA0.toFixed(3)}`);
console.log(`   → ${Math.abs(baselineResult.rVfResidA0) < 0.2 ? 'NO — baseline NFW does NOT produce VfResid channel' :
              Math.abs(baselineResult.rVfResidA0) < 0.4 ? 'WEAK — some signal from NFW diversity alone' :
              'YES — NFW diversity already produces VfResid channel'}`);

console.log('\n2. Which mechanism generates VfResid most like SPARC?');
const byMatch = allResults.slice().sort((a, b) => {
  const dA = Math.abs(Math.abs(a.rVfResidA0) - 0.70);
  const dB = Math.abs(Math.abs(b.rVfResidA0) - 0.70);
  return dA - dB;
});
console.log(`   Best match: ${byMatch[0].name} (r = ${byMatch[0].rVfResidA0.toFixed(3)} vs SPARC 0.70)`);
console.log(`   Runner-up:  ${byMatch[1].name} (r = ${byMatch[1].rVfResidA0.toFixed(3)})`);

console.log('\n3. Does the regime pattern (low→high strengthening) appear?');
for (const r of allResults) {
  if (r.regimeLow !== null && r.regimeHigh !== null) {
    const strengthens = Math.abs(r.regimeHigh) > Math.abs(r.regimeLow) + 0.05;
    console.log(`   ${r.name}: low=${r.regimeLow.toFixed(3)}, high=${r.regimeHigh.toFixed(3)} → ${strengthens ? 'YES' : 'NO'}`);
  }
}

console.log('\n4. Does VfResid carry non-structural information?');
for (const r of allResults) {
  console.log(`   ${r.name}: struct→A0 R²=${r.structA0R2.toFixed(3)}, +VfResid R²=${r.structPlusVfR_A0R2.toFixed(3)}, lift=+${(r.vfResidLift*100).toFixed(1)}pp`);
}

const outPath = path.join(__dirname, '..', 'public', 'phase400-vfresid-origin.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '400',
  title: 'Origin of VfResid — Mock Galaxy Simulation Program',
  timestamp: new Date().toISOString(),
  NGAL,
  scenarios: allResults.map(r => ({
    name: r.name,
    desc: r.desc,
    N: r.N,
    rVfResidA0: +r.rVfResidA0.toFixed(4),
    rVflatA0: +r.rVflatA0.toFixed(4),
    rHaloK_VfResid: +r.rHaloK.toFixed(4),
    rLogM200_VfResid: +r.rLogM200.toFixed(4),
    rConc_VfResid: +r.rConc.toFixed(4),
    structA0R2: +r.structA0R2.toFixed(4),
    structPlusVfR_A0R2: +r.structPlusVfR_A0R2.toFixed(4),
    vfResidLift: +r.vfResidLift.toFixed(4),
    rHaloK_A0: +r.rHaloKA0.toFixed(4),
    rLogM200_A0: +r.rLogM200A0.toFixed(4),
    regime: {
      lowV: r.regimeLow !== null ? +r.regimeLow.toFixed(4) : null,
      highV: r.regimeHigh !== null ? +r.regimeHigh.toFixed(4) : null,
      vHighV: r.regimeVHigh !== null ? +r.regimeVHigh.toFixed(4) : null
    },
    mechanismCorrs: {
      dynTimeProxy: +r.mechanismCorrs.rDynTime.toFixed(4),
      feedbackStrength: +r.mechanismCorrs.rFeedback.toFixed(4),
      assemblyProxy: +r.mechanismCorrs.rAssembly.toFixed(4),
      contractionStrength: +r.mechanismCorrs.rContraction.toFixed(4)
    }
  })),
  sparcComparison: {
    rVfResidA0: 0.70,
    rHaloK: 0.60,
    vfResidLift: 0.171,
    regimeLow: 0.30,
    regimeHigh: 0.75,
    regimeVHigh: 0.83
  }
}, null, 2));

console.log(`\nResults saved to: ${outPath}`);
console.log('\nPhase 400A complete.');
