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
console.log('PHASE 401: What Physics Flips the Regime Sign?');
console.log('='.repeat(70));
console.log('\nGoal: Find a mechanism that produces STRENGTHENING at high-V');
console.log('      (not just coupling — the CORRECT regime direction)');
console.log('\nSuccess criterion: r(VfResid,a0) at high-V > r at low-V');
console.log('SPARC reference: low=0.30, high=0.75 (strengthens)');
console.log('Phase 400 baseline: low=0.93, high=0.87 (weakens)');
console.log('');

const NGAL = 600;

function generateGalaxy(logM200, opts) {
  const M200 = Math.pow(10, logM200);
  const massFrac = Math.max(0, Math.min(1, (logM200 - 10) / 2.5));

  const cBase = 10 * Math.pow(M200 / 1e12, -0.1);
  let concScatter = randNorm(0, 0.12);

  if (opts.concentrationScatterReduction) {
    const reductionAtHighMass = opts.concScatterReductionStrength || 0.6;
    const effectiveScatter = 0.12 * (1 - massFrac * reductionAtHighMass);
    concScatter = randNorm(0, effectiveScatter);
  }

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

  let contractionFactor = 1.0;
  if (opts.massDepContraction) {
    const contrStrength = opts.contractionStrength || 0.4;
    const fbNorm = fb / 0.05;
    const cNorm = c / cBase;
    contractionFactor = 1 + contrStrength * massFrac * fbNorm * (1 + 0.3 * (cNorm - 1));
    contractionFactor = Math.max(1.0, Math.min(contractionFactor, 2.5));
  }

  let lowMassNonCirc = 0;
  if (opts.stochasticLowMass) {
    const ncStrength = opts.nonCircStrength || 0.3;
    lowMassNonCirc = ncStrength * (1 - massFrac) * randNorm(0, 1);
  }

  let responseAntiCorr = 0;
  if (opts.responseAntiCorrelation) {
    const respStr = opts.responseStrength || 0.3;
    responseAntiCorr = -respStr * concScatter * massFrac;
  }

  let dynIntegration = 0;
  if (opts.deepWellIntegration) {
    const intStr = opts.integrationStrength || 0.25;
    dynIntegration = intStr * massFrac * massFrac * (1 + 0.5 * randNorm(0, 0.3));
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

    const rContract = r / contractionFactor;
    const MdmEnc = nfwMenc(rContract, M200, c);

    const Mtot = MbarEnc + MdmEnc;
    let VobsTrue = Math.sqrt(G * Mtot / r);

    if (lowMassNonCirc !== 0) {
      VobsTrue *= (1 + 0.15 * lowMassNonCirc * Math.exp(-r / (2 * Rd)));
    }

    if (responseAntiCorr !== 0) {
      const gobsBase = VobsTrue * VobsTrue / r;
      const gbarBase = G * MbarEnc / (r * r);
      if (gbarBase > 0) {
        const rarBase = rarPredict(gbarBase, A0_TRUE);
        const a0Shifted = A0_TRUE * Math.pow(10, responseAntiCorr);
        const rarShifted = rarPredict(gbarBase, a0Shifted);
        if (rarBase > 0) {
          VobsTrue *= Math.sqrt(rarShifted / rarBase);
        }
      }
    }

    if (dynIntegration !== 0) {
      const gbarBase = G * MbarEnc / (r * r);
      if (gbarBase > 0) {
        const rarBase = rarPredict(gbarBase, A0_TRUE);
        const a0Shifted = A0_TRUE * Math.pow(10, dynIntegration);
        const rarShifted = rarPredict(gbarBase, a0Shifted);
        if (rarBase > 0) {
          VobsTrue *= Math.sqrt(rarShifted / rarBase);
        }
      }
    }

    const Vobs = VobsTrue * (1 + randNorm(0, velNoise));
    const MbarObs = MbarEnc * Math.pow(10, YstarErr);
    const Vbar = Math.sqrt(G * MbarObs / r);

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
    logC: Math.log10(c), concScatter, massFrac,
    contractionFactor, lowMassNonCirc, responseAntiCorr, dynIntegration,
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

  console.log(`  ${label}: N=${valid.length}, r(all)=${rAll.toFixed(3)}`);
  console.log(`    low-V(N=${lowV.length}): ${rLow !== null ? rLow.toFixed(3) : 'N/A'}`);
  console.log(`    high-V(N=${highV.length}): ${rHigh !== null ? rHigh.toFixed(3) : 'N/A'}`);
  console.log(`    v-high-V(N=${vhighV.length}): ${rVHigh !== null ? rVHigh.toFixed(3) : 'N/A'}`);
  console.log(`    → ${sign} ${strengthens ? '✓✓✓' : '✗'}`);

  return { label, N: valid.length, rAll, rLow, rHigh, rVHigh, strengthens };
}

function runExperiment(label, desc, opts) {
  console.log(`\n${'═'.repeat(65)}`);
  console.log(`EXPERIMENT: ${label}`);
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
console.log('E1: STOCHASTIC LOW-MASS DISRUPTION');
console.log('    Feedback/non-circular motions scramble concentration');
console.log('    signal at low mass → should weaken low-V r');
console.log('▓'.repeat(70));

for (const str of [0.1, 0.2, 0.3, 0.5, 0.8]) {
  allResults.push(runExperiment(
    `E1_nonCirc_${str}`,
    `Low-mass non-circular motions strength=${str}`,
    { stochasticLowMass: true, nonCircStrength: str }
  ));
}


console.log('\n' + '▓'.repeat(70));
console.log('E2: MASS-DEPENDENT ADIABATIC CONTRACTION');
console.log('    Baryons contract DM halo more in massive galaxies');
console.log('    → should add new signal at high mass');
console.log('▓'.repeat(70));

for (const str of [0.2, 0.4, 0.6, 0.8, 1.0]) {
  allResults.push(runExperiment(
    `E2_contract_${str}`,
    `Mass-dep contraction strength=${str}`,
    { massDepContraction: true, contractionStrength: str }
  ));
}


console.log('\n' + '▓'.repeat(70));
console.log('E3: COMBINED — Low-mass disruption + High-mass contraction');
console.log('    The hypothesis: both effects together flip the sign');
console.log('▓'.repeat(70));

const combos = [
  { nc: 0.2, ct: 0.3 },
  { nc: 0.3, ct: 0.4 },
  { nc: 0.5, ct: 0.5 },
  { nc: 0.5, ct: 0.8 },
  { nc: 0.8, ct: 0.6 },
  { nc: 0.8, ct: 1.0 },
  { nc: 1.0, ct: 0.8 },
];

for (const combo of combos) {
  allResults.push(runExperiment(
    `E3_nc${combo.nc}_ct${combo.ct}`,
    `NonCirc=${combo.nc} + Contraction=${combo.ct}`,
    {
      stochasticLowMass: true, nonCircStrength: combo.nc,
      massDepContraction: true, contractionStrength: combo.ct
    }
  ));
}


console.log('\n' + '▓'.repeat(70));
console.log('E4: CONCENTRATION-RESPONSE ANTI-CORRELATION');
console.log('    Low-c halos get MORE contracted by baryons at high mass');
console.log('    → reduces effective c scatter → destroys baseline signal');
console.log('    + adds new physics signal that grows with mass');
console.log('▓'.repeat(70));

for (const str of [0.2, 0.4, 0.6, 0.8, 1.0]) {
  allResults.push(runExperiment(
    `E4_antiCorr_${str}`,
    `Conc-response anti-correlation strength=${str}`,
    {
      responseAntiCorrelation: true, responseStrength: str,
      stochasticLowMass: true, nonCircStrength: 0.3,
      massDepContraction: true, contractionStrength: 0.3
    }
  ));
}


console.log('\n' + '▓'.repeat(70));
console.log('E5: DEEP-WELL DYNAMICAL INTEGRATION');
console.log('    Accumulated baryon-halo processing scales as massFrac^2');
console.log('    → strongly nonlinear growth with potential depth');
console.log('▓'.repeat(70));

for (const str of [0.15, 0.25, 0.35, 0.5]) {
  allResults.push(runExperiment(
    `E5_dynInt_${str}`,
    `Deep-well integration strength=${str}`,
    {
      deepWellIntegration: true, integrationStrength: str,
      stochasticLowMass: true, nonCircStrength: 0.5,
      massDepContraction: true, contractionStrength: 0.4
    }
  ));
}


console.log('\n' + '▓'.repeat(70));
console.log('E6: FULL KITCHEN SINK');
console.log('    All mechanisms + concentration scatter reduction at high mass');
console.log('▓'.repeat(70));

const kitchenSink = [
  { nc: 0.5, ct: 0.5, ac: 0.4, di: 0.25, csr: 0.3 },
  { nc: 0.5, ct: 0.6, ac: 0.5, di: 0.3, csr: 0.5 },
  { nc: 0.8, ct: 0.8, ac: 0.6, di: 0.35, csr: 0.6 },
  { nc: 1.0, ct: 1.0, ac: 0.8, di: 0.5, csr: 0.7 },
];

for (const ks of kitchenSink) {
  allResults.push(runExperiment(
    `E6_full_nc${ks.nc}_ct${ks.ct}_di${ks.di}`,
    `Full: nc=${ks.nc} ct=${ks.ct} ac=${ks.ac} di=${ks.di} csr=${ks.csr}`,
    {
      stochasticLowMass: true, nonCircStrength: ks.nc,
      massDepContraction: true, contractionStrength: ks.ct,
      responseAntiCorrelation: true, responseStrength: ks.ac,
      deepWellIntegration: true, integrationStrength: ks.di,
      concentrationScatterReduction: true, concScatterReductionStrength: ks.csr
    }
  ));
}


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 401 RESULTS SUMMARY');
console.log('═'.repeat(70));

function pad(s, n) { return String(s).padEnd(n); }

console.log('\nSPARC reference:  r(all)=0.70  low=0.30  high=0.75  STRENGTHENS');
console.log('Phase400 baseline: r(all)=0.88  low=0.93  high=0.87  weakens');
console.log('');
console.log(pad('Experiment', 40) + pad('r(all)', 10) + pad('low', 10) + pad('high', 10) + pad('v-high', 10) + 'Pattern');
console.log('─'.repeat(90));

const strengthened = [];
const bestMatch = { score: -Infinity, result: null };

for (const r of allResults) {
  const pattern = r.strengthens ? 'STRENGTHENS ✓✓✓' : 'weakens';
  console.log(
    pad(r.label, 40) +
    pad(r.rAll.toFixed(3), 10) +
    pad(r.rLow !== null ? r.rLow.toFixed(3) : 'N/A', 10) +
    pad(r.rHigh !== null ? r.rHigh.toFixed(3) : 'N/A', 10) +
    pad(r.rVHigh !== null ? r.rVHigh.toFixed(3) : 'N/A', 10) +
    pattern
  );

  if (r.strengthens) strengthened.push(r);

  let score = 0;
  if (r.strengthens) score += 3;
  if (r.rAll !== null && Math.abs(r.rAll - 0.70) < 0.15) score += 1;
  if (r.rLow !== null && Math.abs(r.rLow - 0.30) < 0.2) score += 1;
  if (r.rHigh !== null && Math.abs(r.rHigh - 0.75) < 0.2) score += 1;
  if (score > bestMatch.score) { bestMatch.score = score; bestMatch.result = r; }
}

console.log('─'.repeat(90));

console.log(`\n${'═'.repeat(70)}`);
console.log('VERDICT');
console.log('═'.repeat(70));

if (strengthened.length > 0) {
  console.log(`\n✓ REGIME STRENGTHENING ACHIEVED in ${strengthened.length} configurations!`);
  console.log('\nSuccessful configurations:');
  for (const s of strengthened) {
    console.log(`  ${s.label}: low=${s.rLow !== null ? s.rLow.toFixed(3) : '?'}, high=${s.rHigh !== null ? s.rHigh.toFixed(3) : '?'}, v-high=${s.rVHigh !== null ? s.rVHigh.toFixed(3) : '?'}`);
  }
  console.log(`\nBest SPARC match: ${bestMatch.result.label} (score=${bestMatch.score}/6)`);
  if (bestMatch.result.rLow !== null && bestMatch.result.rHigh !== null) {
    console.log(`  low=${bestMatch.result.rLow.toFixed(3)} (SPARC: 0.30)`);
    console.log(`  high=${bestMatch.result.rHigh.toFixed(3)} (SPARC: 0.75)`);
    console.log(`  v-high=${bestMatch.result.rVHigh !== null ? bestMatch.result.rVHigh.toFixed(3) : 'N/A'} (SPARC: 0.83)`);
  }

  console.log('\nRequired ingredients for regime flip:');
  const hasNC = strengthened.some(s => s.label.includes('nonCirc') || s.label.includes('nc'));
  const hasCT = strengthened.some(s => s.label.includes('contract') || s.label.includes('ct'));
  const hasAC = strengthened.some(s => s.label.includes('antiCorr') || s.label.includes('ac'));
  const hasDI = strengthened.some(s => s.label.includes('dynInt') || s.label.includes('di'));
  if (hasNC) console.log('  ✓ Stochastic low-mass disruption (non-circular motions)');
  if (hasCT) console.log('  ✓ Mass-dependent adiabatic contraction');
  if (hasAC) console.log('  ✓ Concentration-response anti-correlation');
  if (hasDI) console.log('  ✓ Deep-well dynamical integration');
} else {
  console.log('\n✗ REGIME STRENGTHENING NOT ACHIEVED in any configuration.');
  console.log('  The concentration baseline remains dominant.');
  console.log('  Either the physics tested is insufficient, or the');
  console.log('  SPARC pattern has a qualitatively different origin.');

  console.log('\nClosest approaches:');
  const sorted = allResults.slice().sort((a, b) => {
    const da = a.rHigh !== null && a.rLow !== null ? a.rHigh - a.rLow : -999;
    const db = b.rHigh !== null && b.rLow !== null ? b.rHigh - b.rLow : -999;
    return db - da;
  });
  for (let i = 0; i < Math.min(5, sorted.length); i++) {
    const s = sorted[i];
    if (s.rLow !== null && s.rHigh !== null) {
      const delta = s.rHigh - s.rLow;
      console.log(`  ${s.label}: delta=${delta > 0 ? '+' : ''}${delta.toFixed(3)} (low=${s.rLow.toFixed(3)}, high=${s.rHigh.toFixed(3)})`);
    }
  }
}

const outPath = path.join(__dirname, '..', 'public', 'phase401-regime-flip.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '401',
  title: 'What Physics Flips the Regime Sign?',
  timestamp: new Date().toISOString(),
  NGAL,
  successCriterion: 'r(VfResid,a0) at high-V > r at low-V (regime strengthening)',
  sparcReference: { rAll: 0.70, rLow: 0.30, rHigh: 0.75, rVHigh: 0.83 },
  experiments: allResults.map(r => ({
    label: r.label,
    N: r.N,
    rAll: r.rAll !== null ? +r.rAll.toFixed(4) : null,
    rLow: r.rLow !== null ? +r.rLow.toFixed(4) : null,
    rHigh: r.rHigh !== null ? +r.rHigh.toFixed(4) : null,
    rVHigh: r.rVHigh !== null ? +r.rVHigh.toFixed(4) : null,
    strengthens: r.strengthens
  })),
  nStrengthened: strengthened.length,
  bestMatch: bestMatch.result ? {
    label: bestMatch.result.label,
    score: bestMatch.score,
    rAll: bestMatch.result.rAll !== null ? +bestMatch.result.rAll.toFixed(4) : null,
    rLow: bestMatch.result.rLow !== null ? +bestMatch.result.rLow.toFixed(4) : null,
    rHigh: bestMatch.result.rHigh !== null ? +bestMatch.result.rHigh.toFixed(4) : null,
    rVHigh: bestMatch.result.rVHigh !== null ? +bestMatch.result.rVHigh.toFixed(4) : null
  } : null
}, null, 2));

console.log(`\nResults saved to: ${outPath}`);
