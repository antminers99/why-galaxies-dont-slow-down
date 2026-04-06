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

function burkertMenc(r, M200, c) {
  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const r0 = r200 / c;
  const x = r / r0;
  const enclosed = Math.log(1 + x) + 0.5 * Math.log(1 + x * x) - Math.atan(x);
  const cNorm = Math.log(1 + c) + 0.5 * Math.log(1 + c * c) - Math.atan(c);
  return M200 * enclosed / cNorm;
}

function pseudoIsoMenc(r, M200, c) {
  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const rc = r200 / c;
  const enclosed = r - rc * Math.atan(r / rc);
  const cNorm = r200 - rc * Math.atan(r200 / rc);
  return M200 * enclosed / Math.max(cNorm, 1e-10);
}

function dc14Menc(r, M200, c, Mstar) {
  const logMstarMhalo = Math.log10(Mstar / M200);
  const X = logMstarMhalo + 4.1;
  const alpha = 2.94 - Math.log10(Math.pow(10, 4.899 * X) + Math.pow(10, 2.641 * X));
  const beta = 4.23 + 1.309 * X + 0.0 * X * X;
  const gamma = Math.max(0, Math.min(1.3, -0.06 + Math.log10(Math.pow(10, -0.68 * X) + Math.pow(10, 0.0 * X))));

  const r200 = Math.pow(M200 / (4/3 * Math.PI * 200 * rhoCrit), 1/3);
  const rs = r200 / c;
  const x = r / rs;

  const nSteps = 50;
  const dr = r / nSteps;
  let mass = 0;
  for (let i = 0; i < nSteps; i++) {
    const ri = (i + 0.5) * dr;
    const xi = ri / rs;
    const rho = 1 / (Math.pow(xi, gamma) * Math.pow(1 + Math.pow(xi, alpha), (beta - gamma) / alpha));
    mass += 4 * Math.PI * ri * ri * rho * dr;
  }

  const drN = r200 / nSteps;
  let massNorm = 0;
  for (let i = 0; i < nSteps; i++) {
    const ri = (i + 0.5) * drN;
    const xi = ri / rs;
    const rho = 1 / (Math.pow(xi, gamma) * Math.pow(1 + Math.pow(xi, alpha), (beta - gamma) / alpha));
    massNorm += 4 * Math.PI * ri * ri * rho * drN;
  }

  return massNorm > 0 ? M200 * mass / massNorm : nfwMenc(r, M200, c);
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
console.log('PHASE 402: BEYOND NFW — Cored Haloes + Genuine Physics');
console.log('='.repeat(70));
console.log('\nPhase 401 proved: NFW+disk CANNOT produce regime strengthening.');
console.log('Now testing qualitatively different halo physics.');
console.log('\nSuccess = r(high-V) > r(low-V)  [regime STRENGTHENING]');
console.log('SPARC: low=0.30, high=0.75');
console.log('');

const NGAL = 600;

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

  let haloModel = opts.haloProfile || 'nfw';
  if (opts.massTransition) {
    const transitionMass = opts.transitionLogM || 11.0;
    const transWidth = opts.transitionWidth || 0.5;
    const sigmoid = 1 / (1 + Math.exp(-(logM200 - transitionMass) / transWidth));
    haloModel = Math.random() < sigmoid ? 'nfw' : opts.lowMassProfile || 'burkert';
  }

  let a0Phys = A0_TRUE;
  if (opts.genuineA0Variation) {
    const a0Strength = opts.a0Strength || 0.15;
    const a0MassDep = opts.a0MassDep || 0;
    const a0Scatter = opts.a0Scatter || 0.1;

    let a0Shift = 0;

    if (opts.a0ActivatesHighMass) {
      a0Shift = a0Strength * massFrac * massFrac * (1 + randNorm(0, a0Scatter));
    } else if (opts.a0ActivatesPotentialDepth) {
      const Vcirc200 = Math.sqrt(G * M200 / r200);
      const VcircNorm = Vcirc200 / 150;
      a0Shift = a0Strength * VcircNorm * VcircNorm * (1 + randNorm(0, a0Scatter));
    } else {
      a0Shift = a0Strength * randNorm(0, 1) + a0MassDep * massFrac;
    }

    a0Phys = A0_TRUE * Math.pow(10, a0Shift);
  }

  let lowMassNonCirc = 0;
  if (opts.stochasticLowMass) {
    const ncStr = opts.nonCircStrength || 0.3;
    lowMassNonCirc = ncStr * (1 - massFrac) * randNorm(0, 1);
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
    if (haloModel === 'burkert') {
      MdmEnc = burkertMenc(r, M200, c);
    } else if (haloModel === 'piso') {
      MdmEnc = pseudoIsoMenc(r, M200, c);
    } else if (haloModel === 'dc14') {
      MdmEnc = dc14Menc(r, M200, c, Mstar);
    } else {
      MdmEnc = nfwMenc(r, M200, c);
    }

    const Mtot = MbarEnc + MdmEnc;
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

    if (lowMassNonCirc !== 0) {
      VobsTrue *= (1 + 0.15 * lowMassNonCirc * Math.exp(-r / (2 * Rd)));
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
    logM200, c, fb, fgas: fgasClamped, Mstar, Mgas, Mbar, Rd, Rg,
    Vflat: Math.max(Vflat, 1), logVflat: Math.log10(Math.max(Vflat, 1)),
    logMbar: Math.log10(Mbar), logL36: Math.log10(L36), logRd: Math.log10(Rd),
    morphT, logMHI: Math.log10(Mgas / 1e9),
    a0: fit.a0, logA0: Math.log10(fit.a0), rmse: fit.rmse,
    logC: Math.log10(c), concScatter, massFrac, haloModel,
    a0Phys, logA0Phys: Math.log10(a0Phys),
    nPts: points.length
  };
}

function analyzeRegime(label, galaxies) {
  const valid = galaxies.filter(g => g.rmse < 0.5 && isFinite(g.logA0) && g.Vflat > 15 && g.nPts >= 4);
  computeVfResid(valid);

  const rAll = pearsonR(valid.map(g => g.VfResid), valid.map(g => g.logA0));

  const lowV = valid.filter(g => g.Vflat < 120);
  const highV = valid.filter(g => g.Vflat >= 120);
  const vhighV = valid.filter(g => g.Vflat >= 180);

  const rLow = lowV.length >= 8 ? pearsonR(lowV.map(g => g.VfResid), lowV.map(g => g.logA0)) : null;
  const rHigh = highV.length >= 8 ? pearsonR(highV.map(g => g.VfResid), highV.map(g => g.logA0)) : null;
  const rVHigh = vhighV.length >= 8 ? pearsonR(vhighV.map(g => g.VfResid), vhighV.map(g => g.logA0)) : null;

  const strengthens = rLow !== null && rHigh !== null && rHigh > rLow + 0.05;
  const sign = strengthens ? 'STRENGTHENS' : 'weakens';

  const delta = rLow !== null && rHigh !== null ? rHigh - rLow : null;

  console.log(`  ${label}: N=${valid.length}, r(all)=${rAll.toFixed(3)}`);
  console.log(`    low-V(N=${lowV.length}): ${rLow !== null ? rLow.toFixed(3) : 'N/A'}`);
  console.log(`    high-V(N=${highV.length}): ${rHigh !== null ? rHigh.toFixed(3) : 'N/A'}`);
  console.log(`    v-high-V(N=${vhighV.length}): ${rVHigh !== null ? rVHigh.toFixed(3) : 'N/A'}`);
  console.log(`    delta(high-low) = ${delta !== null ? (delta > 0 ? '+' : '') + delta.toFixed(3) : 'N/A'}`);
  console.log(`    → ${sign} ${strengthens ? '✓✓✓' : '✗'}`);

  return { label, N: valid.length, rAll, rLow, rHigh, rVHigh, delta, strengthens };
}

function runExperiment(label, desc, opts) {
  console.log(`\n${'═'.repeat(65)}`);
  console.log(`${label}`);
  console.log(`  ${desc}`);
  console.log('═'.repeat(65));

  const galaxies = [];
  for (let i = 0; i < NGAL; i++) {
    const logM200 = rand(10.0, 12.5);
    galaxies.push(generateGalaxy(logM200, opts));
  }

  return analyzeRegime(label, galaxies);
}

const allResults = [];


console.log('\n' + '▓'.repeat(70));
console.log('BRANCH A: CORED LOW-MASS HALOES');
console.log('    Break concentration dominance at low mass');
console.log('▓'.repeat(70));

allResults.push(runExperiment('A1_burkert_all', 'All galaxies use Burkert profile', { haloProfile: 'burkert' }));
allResults.push(runExperiment('A2_piso_all', 'All galaxies use pseudo-isothermal', { haloProfile: 'piso' }));
allResults.push(runExperiment('A3_dc14_all', 'All galaxies use DC14 profile', { haloProfile: 'dc14' }));
allResults.push(runExperiment('A4_burkert_transition', 'Burkert at low mass, NFW at high mass (logM=11)', {
  massTransition: true, lowMassProfile: 'burkert', transitionLogM: 11.0
}));
allResults.push(runExperiment('A5_piso_transition', 'PseudoIso at low mass, NFW at high mass', {
  massTransition: true, lowMassProfile: 'piso', transitionLogM: 11.0
}));
allResults.push(runExperiment('A6_dc14_transition', 'DC14 at low mass, NFW at high mass', {
  massTransition: true, lowMassProfile: 'dc14', transitionLogM: 11.0
}));


console.log('\n' + '▓'.repeat(70));
console.log('BRANCH B: GENUINE HIGH-MASS a₀ ACTIVATION');
console.log('    Physical a₀ variation that activates in deep potential wells');
console.log('▓'.repeat(70));

allResults.push(runExperiment('B1_a0_highMass_weak', 'a₀ activates at high mass (weak, σ=0.10)', {
  genuineA0Variation: true, a0ActivatesHighMass: true, a0Strength: 0.10, a0Scatter: 0.2
}));
allResults.push(runExperiment('B2_a0_highMass_mod', 'a₀ activates at high mass (moderate, σ=0.15)', {
  genuineA0Variation: true, a0ActivatesHighMass: true, a0Strength: 0.15, a0Scatter: 0.2
}));
allResults.push(runExperiment('B3_a0_highMass_strong', 'a₀ activates at high mass (strong, σ=0.20)', {
  genuineA0Variation: true, a0ActivatesHighMass: true, a0Strength: 0.20, a0Scatter: 0.2
}));
allResults.push(runExperiment('B4_a0_potDepth', 'a₀ activates with potential depth (Vcirc²)', {
  genuineA0Variation: true, a0ActivatesPotentialDepth: true, a0Strength: 0.15, a0Scatter: 0.2
}));
allResults.push(runExperiment('B5_a0_potDepth_strong', 'a₀ activates with potential depth (strong)', {
  genuineA0Variation: true, a0ActivatesPotentialDepth: true, a0Strength: 0.25, a0Scatter: 0.15
}));


console.log('\n' + '▓'.repeat(70));
console.log('BRANCH C: HYBRID — Cored low-mass + Genuine high-mass physics');
console.log('    The leading hypothesis for regime sign flip');
console.log('▓'.repeat(70));

allResults.push(runExperiment('C1_burkert+a0_weak', 'Burkert→NFW transition + weak a₀ activation', {
  massTransition: true, lowMassProfile: 'burkert', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesHighMass: true, a0Strength: 0.10, a0Scatter: 0.2
}));
allResults.push(runExperiment('C2_burkert+a0_mod', 'Burkert→NFW + moderate a₀ activation', {
  massTransition: true, lowMassProfile: 'burkert', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesHighMass: true, a0Strength: 0.15, a0Scatter: 0.2
}));
allResults.push(runExperiment('C3_burkert+a0_strong', 'Burkert→NFW + strong a₀ activation', {
  massTransition: true, lowMassProfile: 'burkert', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesHighMass: true, a0Strength: 0.20, a0Scatter: 0.2
}));
allResults.push(runExperiment('C4_burkert+potDepth', 'Burkert→NFW + potential-depth a₀', {
  massTransition: true, lowMassProfile: 'burkert', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesPotentialDepth: true, a0Strength: 0.20, a0Scatter: 0.15
}));
allResults.push(runExperiment('C5_piso+a0_strong', 'PseudoIso→NFW + strong a₀ activation', {
  massTransition: true, lowMassProfile: 'piso', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesHighMass: true, a0Strength: 0.20, a0Scatter: 0.2
}));
allResults.push(runExperiment('C6_dc14+a0_strong', 'DC14→NFW + strong a₀ activation', {
  massTransition: true, lowMassProfile: 'dc14', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesHighMass: true, a0Strength: 0.20, a0Scatter: 0.2
}));


console.log('\n' + '▓'.repeat(70));
console.log('BRANCH C+: HYBRID + LOW-MASS DISRUPTION');
console.log('    Maximum realism: cores + noise + genuine physics');
console.log('▓'.repeat(70));

allResults.push(runExperiment('C+1_burkert+a0+nc', 'Burkert→NFW + a₀ activation + non-circ', {
  massTransition: true, lowMassProfile: 'burkert', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesHighMass: true, a0Strength: 0.15, a0Scatter: 0.2,
  stochasticLowMass: true, nonCircStrength: 0.5
}));
allResults.push(runExperiment('C+2_burkert+potDepth+nc', 'Burkert→NFW + potential-depth + non-circ', {
  massTransition: true, lowMassProfile: 'burkert', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesPotentialDepth: true, a0Strength: 0.20, a0Scatter: 0.15,
  stochasticLowMass: true, nonCircStrength: 0.5
}));
allResults.push(runExperiment('C+3_piso+potDepth+nc', 'PseudoIso→NFW + potential-depth + non-circ', {
  massTransition: true, lowMassProfile: 'piso', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesPotentialDepth: true, a0Strength: 0.25, a0Scatter: 0.15,
  stochasticLowMass: true, nonCircStrength: 0.5
}));
allResults.push(runExperiment('C+4_full_best', 'Burkert→NFW + strong potDepth + strong non-circ', {
  massTransition: true, lowMassProfile: 'burkert', transitionLogM: 10.8,
  genuineA0Variation: true, a0ActivatesPotentialDepth: true, a0Strength: 0.30, a0Scatter: 0.1,
  stochasticLowMass: true, nonCircStrength: 0.8
}));
allResults.push(runExperiment('C+5_dc14+full', 'DC14→NFW + potDepth + non-circ (max realism)', {
  massTransition: true, lowMassProfile: 'dc14', transitionLogM: 11.0,
  genuineA0Variation: true, a0ActivatesPotentialDepth: true, a0Strength: 0.25, a0Scatter: 0.15,
  stochasticLowMass: true, nonCircStrength: 0.5
}));


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 402 RESULTS SUMMARY');
console.log('═'.repeat(70));

function pad(s, n) { return String(s).padEnd(n); }

console.log('\nSPARC reference:  r(all)=0.70  low=0.30  high=0.75  delta=+0.45  STRENGTHENS');
console.log('Phase401 best:    r(all)=0.60  low=0.64  high=0.62  delta=-0.03  weakens');
console.log('');
console.log(pad('Experiment', 35) + pad('r(all)', 9) + pad('low', 9) + pad('high', 9) + pad('v-high', 9) + pad('delta', 9) + 'Pattern');
console.log('─'.repeat(90));

const strengthened = [];

for (const r of allResults) {
  const pattern = r.strengthens ? 'STRENGTHENS ✓✓✓' : 'weakens';
  console.log(
    pad(r.label, 35) +
    pad(r.rAll.toFixed(3), 9) +
    pad(r.rLow !== null ? r.rLow.toFixed(3) : 'N/A', 9) +
    pad(r.rHigh !== null ? r.rHigh.toFixed(3) : 'N/A', 9) +
    pad(r.rVHigh !== null ? r.rVHigh.toFixed(3) : 'N/A', 9) +
    pad(r.delta !== null ? (r.delta > 0 ? '+' : '') + r.delta.toFixed(3) : 'N/A', 9) +
    pattern
  );
  if (r.strengthens) strengthened.push(r);
}

console.log('─'.repeat(90));

console.log(`\n${'═'.repeat(70)}`);
console.log('VERDICT');
console.log('═'.repeat(70));

if (strengthened.length > 0) {
  console.log(`\n✓ REGIME STRENGTHENING ACHIEVED in ${strengthened.length}/${allResults.length} configurations!`);
  console.log('\nConfigurations that FLIP the sign:');
  for (const s of strengthened) {
    console.log(`  ${s.label}: delta=${s.delta > 0 ? '+' : ''}${s.delta.toFixed(3)}`);
    console.log(`    low=${s.rLow !== null ? s.rLow.toFixed(3) : '?'}, high=${s.rHigh !== null ? s.rHigh.toFixed(3) : '?'}, v-high=${s.rVHigh !== null ? s.rVHigh.toFixed(3) : '?'}`);
  }

  const best = strengthened.sort((a, b) => (b.delta || 0) - (a.delta || 0))[0];
  console.log(`\nBest match to SPARC:`);
  console.log(`  ${best.label}`);
  console.log(`  r(all)=${best.rAll.toFixed(3)} (SPARC: 0.70)`);
  console.log(`  low=${best.rLow !== null ? best.rLow.toFixed(3) : '?'} (SPARC: 0.30)`);
  console.log(`  high=${best.rHigh !== null ? best.rHigh.toFixed(3) : '?'} (SPARC: 0.75)`);
  console.log(`  delta=${best.delta !== null ? (best.delta > 0 ? '+' : '') + best.delta.toFixed(3) : '?'} (SPARC: +0.45)`);

  console.log('\nRequired ingredients for regime flip:');
  const hasCored = strengthened.some(s => s.label.includes('burkert') || s.label.includes('piso') || s.label.includes('dc14'));
  const hasA0 = strengthened.some(s => s.label.includes('a0') || s.label.includes('potDepth'));
  const hasNC = strengthened.some(s => s.label.includes('nc'));
  const coreOnly = strengthened.some(s => s.label.startsWith('A'));
  const a0Only = strengthened.some(s => s.label.startsWith('B'));
  const hybrid = strengthened.some(s => s.label.startsWith('C'));

  if (coreOnly) console.log('  ✓ Cored haloes ALONE can flip the sign (Branch A)');
  if (a0Only) console.log('  ✓ Genuine a₀ variation ALONE can flip the sign (Branch B)');
  if (hybrid) console.log('  ✓ Hybrid model flips the sign (Branch C)');
  if (hasCored) console.log('  ✓ Cored low-mass haloes break concentration dominance');
  if (hasA0) console.log('  ✓ Genuine a₀ variation adds signal at high mass');
  if (hasNC) console.log('  ✓ Low-mass non-circular motions further weaken low-V signal');
} else {
  console.log('\n✗ REGIME STRENGTHENING NOT ACHIEVED.');
  console.log('  Even cored haloes + genuine physics fail to flip the sign.');
  console.log('  The SPARC pattern may require:');
  console.log('  - Fundamentally non-CDM physics');
  console.log('  - MOND-like acceleration scale');
  console.log('  - Selection effects in SPARC sample construction');

  const sorted = allResults.slice().sort((a, b) => (b.delta || -999) - (a.delta || -999));
  console.log('\nClosest approaches:');
  for (let i = 0; i < Math.min(5, sorted.length); i++) {
    const s = sorted[i];
    console.log(`  ${s.label}: delta=${s.delta !== null ? (s.delta > 0 ? '+' : '') + s.delta.toFixed(3) : 'N/A'}`);
  }
}

const outPath = path.join(__dirname, '..', 'public', 'phase402-beyond-nfw.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '402',
  title: 'Beyond NFW — Cored Haloes + Genuine Physics',
  timestamp: new Date().toISOString(),
  NGAL,
  sparcReference: { rAll: 0.70, rLow: 0.30, rHigh: 0.75, delta: 0.45 },
  experiments: allResults.map(r => ({
    label: r.label, N: r.N,
    rAll: r.rAll !== null ? +r.rAll.toFixed(4) : null,
    rLow: r.rLow !== null ? +r.rLow.toFixed(4) : null,
    rHigh: r.rHigh !== null ? +r.rHigh.toFixed(4) : null,
    rVHigh: r.rVHigh !== null ? +r.rVHigh.toFixed(4) : null,
    delta: r.delta !== null ? +r.delta.toFixed(4) : null,
    strengthens: r.strengthens
  })),
  nStrengthened: strengthened.length,
  branches: {
    A: allResults.filter(r => r.label.startsWith('A')).map(r => r.label),
    B: allResults.filter(r => r.label.startsWith('B')).map(r => r.label),
    C: allResults.filter(r => r.label.startsWith('C')).map(r => r.label)
  }
}, null, 2));

console.log(`\nResults saved to: ${outPath}`);
