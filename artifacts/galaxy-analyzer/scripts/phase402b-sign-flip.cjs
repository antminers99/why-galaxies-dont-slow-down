const fs = require('fs');
const path = require('path');

const G = 4.3009e-6;
const rhoCrit = 136.18;
const A0_TRUE = 3703;

function nfwMenc(r, M200, c) {
  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const rs = r200 / c;
  const x = r / rs;
  const norm = Math.log(1 + c) - c / (1 + c);
  return M200 * (Math.log(1 + x) - x / (1 + x)) / norm;
}

function burkertMenc(r, rho0, r0) {
  const x = r / r0;
  const enclosed = Math.log(1 + x) + 0.5 * Math.log(1 + x * x) - Math.atan(x);
  return 4 * Math.PI * rho0 * r0 * r0 * r0 * enclosed;
}

function burkertFromM200(r, M200, c, coreBoost) {
  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const r0 = r200 / c * (coreBoost || 1.0);
  const x200 = r200 / r0;
  const norm200 = Math.log(1 + x200) + 0.5 * Math.log(1 + x200 * x200) - Math.atan(x200);
  const rho0 = M200 / (4 * Math.PI * r0 * r0 * r0 * norm200);
  return burkertMenc(r, rho0, r0);
}

function feedbackCoredMenc(r, M200, c, Mstar, coreFrac) {
  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const rs = r200 / c;
  const logMstarM200 = Math.log10(Mstar / M200);
  const feedbackEff = Math.max(0, Math.min(1, -0.5 * (logMstarM200 + 2.5)));
  const effectiveCoreFrac = coreFrac * feedbackEff;

  const MdmNFW = nfwMenc(r, M200 * (1 - effectiveCoreFrac), c);
  const rCore = rs * (1 + 2 * feedbackEff);
  const x = r / rCore;
  const coreEnclosed = x - Math.atan(x);
  const x200 = r200 / rCore;
  const coreNorm = x200 - Math.atan(x200);
  const MdmCore = coreNorm > 0 ? M200 * effectiveCoreFrac * coreEnclosed / coreNorm : 0;

  return MdmNFW + MdmCore;
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
  if (gbar <= 0) return gbar;
  return gbar / (1 - Math.exp(-Math.sqrt(gbar / a0)));
}

function fitA0(points) {
  if (points.length < 3) return { a0: 3703, rmse: 1 };
  let bestA0 = 3703, bestRMSE = Infinity;
  for (let logA0 = 2.5; logA0 <= 4.5; logA0 += 0.005) {
    const a0 = Math.pow(10, logA0);
    let sse = 0, cnt = 0;
    for (const p of points) {
      if (p.gbar <= 0 || p.gobs <= 0) continue;
      const pred = rarPredict(p.gbar, a0);
      if (pred <= 0) continue;
      const diff = Math.log10(p.gobs) - Math.log10(pred);
      sse += diff * diff; cnt++;
    }
    if (cnt < 3) continue;
    const rmse = Math.sqrt(sse / cnt);
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

function rand(min, max) { return min + Math.random() * (max - min); }
function randNorm(mu, sigma) {
  const u1 = Math.random(), u2 = Math.random();
  return mu + sigma * Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
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

function computeVfResid(gals) {
  const x = gals.map(g => [g.logMbar, g.logL36, g.logRd, g.morphT]);
  const y = gals.map(g => g.logVflat);
  const n = y.length;
  const nv = 4;
  const mx = Array(nv).fill(0);
  const my = y.reduce((a, b) => a + b, 0) / n;
  for (let j = 0; j < nv; j++) { for (let i = 0; i < n; i++) mx[j] += x[i][j]; mx[j] /= n; }
  const XTX = Array.from({ length: nv }, () => Array(nv).fill(0));
  const XTy = Array(nv).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < nv; j++) {
      XTy[j] += (x[i][j] - mx[j]) * (y[i] - my);
      for (let k = 0; k < nv; k++) XTX[j][k] += (x[i][j] - mx[j]) * (x[i][k] - mx[k]);
    }
  }
  const beta = solveLinear(XTX, XTy);
  for (let i = 0; i < n; i++) {
    let pred = my;
    for (let j = 0; j < nv; j++) pred += beta[j] * (x[i][j] - mx[j]);
    gals[i].VfResid = y[i] - pred;
  }
}

console.log('='.repeat(70));
console.log('PHASE 402b: FOCUSED SIGN-FLIP PROGRAM');
console.log('='.repeat(70));
console.log('\nKey insight from 402: DC14 transition was the ONLY model to flip.');
console.log('Why? Feedback-driven cores at low mass break the c→a₀ link.');
console.log('');
console.log('This phase tests the critical mechanism:');
console.log('  LOW MASS: stochastic core formation → decorrelates c from a₀');
console.log('  HIGH MASS: cuspy NFW preserved → c→a₀ link remains');
console.log('  RESULT: low-V coupling drops, high-V stays → STRENGTHENING');
console.log('');

const NGAL = 800;

function generateGalaxy(logM200, opts) {
  const M200 = Math.pow(10, logM200);
  const massFrac = Math.max(0, Math.min(1, (logM200 - 10) / 2.5));

  const cBase = 10 * Math.pow(M200 / 1e12, -0.1);
  const concScatter = randNorm(0, 0.12);
  const c = cBase * Math.pow(10, concScatter);

  const fbBase = 0.05 * Math.pow(M200 / 1e12, 0.3);
  const fb = Math.min(fbBase * (1 + randNorm(0, 0.15)), 0.16);
  const Mbar = M200 * fb;
  const fgas = 0.1 + 0.5 * Math.pow(M200 / 1e10, -0.3) + randNorm(0, 0.05);
  const fgasClamped = Math.max(0.05, Math.min(0.9, fgas));
  const Mstar = Mbar * (1 - fgasClamped);
  const Mgas = Mbar * fgasClamped;

  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const lambda = 0.035 * (1 + randNorm(0, 0.25));
  const Rd = lambda * r200 / Math.sqrt(2);
  const Rg = Rd * (1.5 + rand(0, 1.5));

  const haloType = opts.haloType || 'nfw';
  const coreScatterStr = opts.coreScatterStrength || 0;
  const coreFrac = opts.coreFraction || 0.3;
  const coreSizeScatter = opts.coreSizeScatter || 0;

  let useBurkert = false;
  let coreBoost = 1.0;

  if (haloType === 'feedbackTransition') {
    const transLogM = opts.transitionLogM || 11.0;
    const transWidth = opts.transitionWidth || 0.5;
    const sigmoid = 1 / (1 + Math.exp(-(logM200 - transLogM) / transWidth));
    useBurkert = Math.random() > sigmoid;

    if (useBurkert && coreSizeScatter > 0) {
      coreBoost = Math.pow(10, coreSizeScatter * randNorm(0, 1) * (1 - massFrac));
      coreBoost = Math.max(0.3, Math.min(5.0, coreBoost));
    }
  } else if (haloType === 'burkertAll') {
    useBurkert = true;
    if (coreSizeScatter > 0) {
      coreBoost = Math.pow(10, coreSizeScatter * randNorm(0, 1));
      coreBoost = Math.max(0.3, Math.min(5.0, coreBoost));
    }
  } else if (haloType === 'feedbackCored') {
    // pass
  }

  let a0Phys = A0_TRUE;
  if (opts.genuineA0) {
    const str = opts.a0Strength || 0.15;
    if (opts.a0PotDepth) {
      const Vcirc200 = Math.sqrt(G * M200 / r200);
      const VcNorm = Vcirc200 / 150;
      a0Phys = A0_TRUE * Math.pow(10, str * VcNorm * VcNorm * (1 + randNorm(0, opts.a0Scatter || 0.15)));
    } else if (opts.a0MassDep) {
      a0Phys = A0_TRUE * Math.pow(10, str * massFrac * massFrac * (1 + randNorm(0, opts.a0Scatter || 0.15)));
    }
  }

  let lowMassNoise = 0;
  if (opts.lowMassNoise) {
    lowMassNoise = (opts.noiseStrength || 0.3) * (1 - massFrac) * randNorm(0, 1);
  }

  const nPts = Math.max(5, Math.min(30, Math.round(6 + 20 * massFrac + randNorm(0, 2))));
  const rMin = 0.3 * Rd;
  const rMax = (3 + 3 * massFrac) * Rd;
  const points = [];
  let Vflat = 0;

  const YstarErr = randNorm(0, 0.10);
  const velNoise = 0.08 * (1 + 1.5 * (1 - massFrac));

  for (let i = 0; i < nPts; i++) {
    const r = rMin + (rMax - rMin) * (i / (nPts - 1));
    const MstarEnc = expDiskMenc(r, Mstar, Rd);
    const MgasEnc = gasDiskMenc(r, Mgas * 1.33, Rg);
    const MbarEnc = MstarEnc + MgasEnc;

    let MdmEnc;
    if (haloType === 'feedbackCored') {
      MdmEnc = feedbackCoredMenc(r, M200, c, Mstar, coreFrac);
    } else if (useBurkert) {
      MdmEnc = burkertFromM200(r, M200, c, coreBoost);
    } else {
      MdmEnc = nfwMenc(r, M200, c);
    }

    if (coreScatterStr > 0 && !useBurkert) {
      const randCore = coreScatterStr * (1 - massFrac) * randNorm(0, 1);
      const coreSup = Math.pow(1 + Math.pow(Rd / r, 2), -Math.abs(randCore));
      MdmEnc *= coreSup;
    }

    const Mtot = MbarEnc + Math.max(MdmEnc, 0);
    let VobsTrue = Math.sqrt(G * Math.max(Mtot, 1) / r);

    if (a0Phys !== A0_TRUE) {
      const gbarBase = G * MbarEnc / (r * r);
      if (gbarBase > 0) {
        const rarBase = rarPredict(gbarBase, A0_TRUE);
        const rarPhys = rarPredict(gbarBase, a0Phys);
        if (rarBase > 0) {
          VobsTrue *= Math.sqrt(rarPhys / rarBase);
        }
      }
    }

    if (lowMassNoise !== 0) {
      VobsTrue *= (1 + 0.15 * lowMassNoise * Math.exp(-r / (2 * Rd)));
    }

    const Vobs = VobsTrue * (1 + randNorm(0, velNoise));
    const MbarObs = MbarEnc * Math.pow(10, YstarErr);
    const Vbar = Math.sqrt(G * Math.max(MbarObs, 1) / r);

    if (Math.abs(Vobs) > Vflat) Vflat = Math.abs(Vobs);
    const gobs = Vobs * Math.abs(Vobs) / r;
    const gbar = Vbar * Vbar / r;
    if (gbar > 0 && gobs > 0) {
      points.push({ r, gobs, gbar });
    }
  }

  const L36 = Mstar / 0.5;
  const morphT = 3 + (logM200 - 11) * 2 + randNorm(0, 1.5);
  const fit = fitA0(points);

  return {
    logM200, c, fb, Mstar, Mgas, Mbar, Rd,
    Vflat: Math.max(Vflat, 1), logVflat: Math.log10(Math.max(Vflat, 1)),
    logMbar: Math.log10(Mbar), logL36: Math.log10(L36), logRd: Math.log10(Rd),
    morphT, logMHI: Math.log10(Mgas / 1e9),
    a0: fit.a0, logA0: Math.log10(fit.a0), rmse: fit.rmse,
    logC: Math.log10(c), concScatter, massFrac,
    useBurkert, coreBoost, a0Phys,
    nPts: points.length
  };
}

function analyzeRegime(label, galaxies) {
  const valid = galaxies.filter(g => g.rmse < 0.5 && isFinite(g.logA0) && g.Vflat > 15 && g.nPts >= 4);
  if (valid.length < 20) {
    console.log(`  ${label}: INSUFFICIENT DATA (N=${valid.length})`);
    return { label, N: valid.length, rAll: 0, rLow: null, rHigh: null, rVHigh: null, delta: null, strengthens: false };
  }
  computeVfResid(valid);

  const rAll = pearsonR(valid.map(g => g.VfResid), valid.map(g => g.logA0));

  const lowV = valid.filter(g => g.Vflat < 120);
  const highV = valid.filter(g => g.Vflat >= 120);
  const vhighV = valid.filter(g => g.Vflat >= 180);

  const rLow = lowV.length >= 8 ? pearsonR(lowV.map(g => g.VfResid), lowV.map(g => g.logA0)) : null;
  const rHigh = highV.length >= 8 ? pearsonR(highV.map(g => g.VfResid), highV.map(g => g.logA0)) : null;
  const rVHigh = vhighV.length >= 8 ? pearsonR(vhighV.map(g => g.VfResid), vhighV.map(g => g.logA0)) : null;

  const strengthens = rLow !== null && rHigh !== null && rHigh > rLow + 0.05;
  const delta = rLow !== null && rHigh !== null ? rHigh - rLow : null;

  console.log(`  ${label}: N=${valid.length}, r(all)=${rAll.toFixed(3)}`);
  console.log(`    low-V(N=${lowV.length}): ${rLow !== null ? rLow.toFixed(3) : 'N/A'}`);
  console.log(`    high-V(N=${highV.length}): ${rHigh !== null ? rHigh.toFixed(3) : 'N/A'}`);
  console.log(`    v-high-V(N=${vhighV.length}): ${rVHigh !== null ? rVHigh.toFixed(3) : 'N/A'}`);
  console.log(`    delta = ${delta !== null ? (delta > 0 ? '+' : '') + delta.toFixed(3) : 'N/A'}`);
  console.log(`    → ${strengthens ? 'STRENGTHENS ✓✓✓' : 'weakens ✗'}`);

  return { label, N: valid.length, rAll, rLow, rHigh, rVHigh, delta, strengthens };
}

function run(label, desc, opts) {
  console.log(`\n${'─'.repeat(60)}`);
  console.log(`${label}: ${desc}`);
  console.log('─'.repeat(60));
  const gals = [];
  for (let i = 0; i < NGAL; i++) {
    gals.push(generateGalaxy(rand(10.0, 12.5), opts));
  }
  return analyzeRegime(label, gals);
}

const R = [];

console.log('\n' + '▓'.repeat(70));
console.log('TEST 1: NFW BASELINE (control)');
console.log('▓'.repeat(70));
R.push(run('T0_nfw_baseline', 'Pure NFW — Phase 400 replication', { haloType: 'nfw' }));

console.log('\n' + '▓'.repeat(70));
console.log('TEST 2: FEEDBACK-CORED HALOES');
console.log('    Mass-dependent core: stronger feedback at low Mstar/Mhalo');
console.log('    Core fraction and size depend on stellar-to-halo mass ratio');
console.log('▓'.repeat(70));

for (const cf of [0.2, 0.4, 0.6]) {
  R.push(run(`T1_fbCore_${cf}`, `Feedback cores, coreFrac=${cf}`, {
    haloType: 'feedbackCored', coreFraction: cf
  }));
}

console.log('\n' + '▓'.repeat(70));
console.log('TEST 3: BURKERT TRANSITION + STOCHASTIC CORE SIZE');
console.log('    Low-mass: Burkert with RANDOM core sizes (uncorrelated with c)');
console.log('    High-mass: NFW (cuspy)');
console.log('    The stochastic cores should break c→a₀ link at low mass');
console.log('▓'.repeat(70));

for (const css of [0.2, 0.4, 0.6, 0.8]) {
  R.push(run(`T2_burkTrans_css${css}`, `Burkert→NFW, core size scatter=${css}`, {
    haloType: 'feedbackTransition', transitionLogM: 11.0, transitionWidth: 0.5,
    coreSizeScatter: css
  }));
}

console.log('\n' + '▓'.repeat(70));
console.log('TEST 4: CORE SCATTER IN NFW (random core injection)');
console.log('    NFW baseline + random mass-dependent core perturbation');
console.log('    Stronger at low mass, absent at high mass');
console.log('▓'.repeat(70));

for (const cs of [0.3, 0.5, 0.8, 1.2]) {
  R.push(run(`T3_coreScatter_${cs}`, `NFW + random core scatter=${cs}`, {
    haloType: 'nfw', coreScatterStrength: cs
  }));
}

console.log('\n' + '▓'.repeat(70));
console.log('TEST 5: HYBRID — Cores + genuine a₀');
console.log('    Feedback cores at low mass + potential-depth a₀ at high mass');
console.log('▓'.repeat(70));

R.push(run('T4a_fbCore+a0', 'FeedbackCore + pot-depth a₀', {
  haloType: 'feedbackCored', coreFraction: 0.4,
  genuineA0: true, a0PotDepth: true, a0Strength: 0.15, a0Scatter: 0.15
}));

R.push(run('T4b_fbCore+a0_strong', 'FeedbackCore + strong pot-depth a₀', {
  haloType: 'feedbackCored', coreFraction: 0.5,
  genuineA0: true, a0PotDepth: true, a0Strength: 0.25, a0Scatter: 0.1
}));

R.push(run('T4c_burkTrans+a0', 'Burkert→NFW + pot-depth a₀ + core scatter', {
  haloType: 'feedbackTransition', transitionLogM: 11.0, coreSizeScatter: 0.5,
  genuineA0: true, a0PotDepth: true, a0Strength: 0.20, a0Scatter: 0.15
}));

R.push(run('T4d_coreScatter+a0', 'NFW + core scatter + pot-depth a₀', {
  haloType: 'nfw', coreScatterStrength: 0.8,
  genuineA0: true, a0PotDepth: true, a0Strength: 0.20, a0Scatter: 0.15
}));

console.log('\n' + '▓'.repeat(70));
console.log('TEST 6: TRIPLE COMBO — Cores + a₀ + low-mass noise');
console.log('▓'.repeat(70));

R.push(run('T5a_triple_mild', 'FbCore + a₀ + noise (mild)', {
  haloType: 'feedbackCored', coreFraction: 0.4,
  genuineA0: true, a0PotDepth: true, a0Strength: 0.15, a0Scatter: 0.15,
  lowMassNoise: true, noiseStrength: 0.3
}));

R.push(run('T5b_triple_mod', 'FbCore + a₀ + noise (moderate)', {
  haloType: 'feedbackCored', coreFraction: 0.5,
  genuineA0: true, a0PotDepth: true, a0Strength: 0.20, a0Scatter: 0.12,
  lowMassNoise: true, noiseStrength: 0.5
}));

R.push(run('T5c_triple_strong', 'FbCore + strong a₀ + strong noise', {
  haloType: 'feedbackCored', coreFraction: 0.6,
  genuineA0: true, a0PotDepth: true, a0Strength: 0.30, a0Scatter: 0.1,
  lowMassNoise: true, noiseStrength: 0.8
}));

R.push(run('T5d_burkTrans_full', 'Burkert→NFW + a₀ + noise + core scatter', {
  haloType: 'feedbackTransition', transitionLogM: 10.8, coreSizeScatter: 0.6,
  genuineA0: true, a0PotDepth: true, a0Strength: 0.25, a0Scatter: 0.12,
  lowMassNoise: true, noiseStrength: 0.5
}));

R.push(run('T5e_coreScatter_full', 'NFW + core scatter + a₀ + noise', {
  haloType: 'nfw', coreScatterStrength: 1.0,
  genuineA0: true, a0PotDepth: true, a0Strength: 0.25, a0Scatter: 0.12,
  lowMassNoise: true, noiseStrength: 0.5
}));

R.push(run('T5f_massDep_a0', 'FbCore + mass-dep a₀ + noise', {
  haloType: 'feedbackCored', coreFraction: 0.5,
  genuineA0: true, a0MassDep: true, a0Strength: 0.25, a0Scatter: 0.12,
  lowMassNoise: true, noiseStrength: 0.5
}));


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 402b RESULTS SUMMARY');
console.log('═'.repeat(70));

function pad(s, n) { return String(s).padEnd(n); }

console.log('\nSPARC:            r(all)=0.70  low=0.30  high=0.75  delta=+0.45');
console.log('');
console.log(pad('Test', 30) + pad('r(all)', 9) + pad('low', 9) + pad('high', 9) + pad('v-high', 9) + pad('delta', 9) + 'Pattern');
console.log('─'.repeat(84));

const strengthened = [];
for (const r of R) {
  const pattern = r.strengthens ? 'STRENGTHENS ✓✓✓' : 'weakens';
  console.log(
    pad(r.label, 30) +
    pad(r.rAll.toFixed(3), 9) +
    pad(r.rLow !== null ? r.rLow.toFixed(3) : 'N/A', 9) +
    pad(r.rHigh !== null ? r.rHigh.toFixed(3) : 'N/A', 9) +
    pad(r.rVHigh !== null ? r.rVHigh.toFixed(3) : 'N/A', 9) +
    pad(r.delta !== null ? (r.delta > 0 ? '+' : '') + r.delta.toFixed(3) : 'N/A', 9) +
    pattern
  );
  if (r.strengthens) strengthened.push(r);
}
console.log('─'.repeat(84));

console.log(`\n${'═'.repeat(70)}`);
console.log('VERDICT');
console.log('═'.repeat(70));

if (strengthened.length > 0) {
  console.log(`\n✓ REGIME SIGN FLIPPED in ${strengthened.length}/${R.length} configs!`);
  const best = strengthened.sort((a, b) => (b.delta || 0) - (a.delta || 0))[0];
  console.log(`\nBest: ${best.label}`);
  console.log(`  delta=${best.delta > 0 ? '+' : ''}${best.delta.toFixed(3)} (SPARC: +0.45)`);
  console.log(`  low=${best.rLow.toFixed(3)} (SPARC: 0.30)`);
  console.log(`  high=${best.rHigh.toFixed(3)} (SPARC: 0.75)`);

  console.log('\nKey ingredients:');
  const hasFb = strengthened.some(s => s.label.includes('fbCore'));
  const hasBurk = strengthened.some(s => s.label.includes('burk'));
  const hasCS = strengthened.some(s => s.label.includes('coreScatter'));
  const hasA0 = strengthened.some(s => s.label.includes('a0'));
  const hasNoise = strengthened.some(s => s.label.includes('triple') || s.label.includes('full'));
  if (hasFb) console.log('  ✓ Feedback-driven cores (mass-dependent)');
  if (hasBurk) console.log('  ✓ Burkert profile transition');
  if (hasCS) console.log('  ✓ Stochastic core size scatter');
  if (hasA0) console.log('  ✓ Genuine a₀ variation (potential-depth dependent)');
  if (hasNoise) console.log('  ✓ Low-mass non-circular motions');
} else {
  console.log('\n✗ NO configuration achieved regime strengthening.');
  const sorted = R.filter(r => r.delta !== null).sort((a, b) => b.delta - a.delta);
  console.log('\nClosest:');
  for (let i = 0; i < Math.min(5, sorted.length); i++) {
    console.log(`  ${sorted[i].label}: delta=${sorted[i].delta > 0 ? '+' : ''}${sorted[i].delta.toFixed(3)}`);
  }
}

const outPath = path.join(__dirname, '..', 'public', 'phase402-beyond-nfw.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '402b',
  title: 'Focused Sign-Flip: Cored Haloes + Genuine Physics',
  timestamp: new Date().toISOString(),
  NGAL,
  sparcReference: { rAll: 0.70, rLow: 0.30, rHigh: 0.75, delta: 0.45 },
  experiments: R.map(r => ({
    label: r.label, N: r.N,
    rAll: +r.rAll.toFixed(4),
    rLow: r.rLow !== null ? +r.rLow.toFixed(4) : null,
    rHigh: r.rHigh !== null ? +r.rHigh.toFixed(4) : null,
    rVHigh: r.rVHigh !== null ? +r.rVHigh.toFixed(4) : null,
    delta: r.delta !== null ? +r.delta.toFixed(4) : null,
    strengthens: r.strengthens
  })),
  nStrengthened: strengthened.length
}, null, 2));

console.log(`\nSaved: ${outPath}`);
