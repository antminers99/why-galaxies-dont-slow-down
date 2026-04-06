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
  const my = y.reduce((a,b)=>a+b,0)/n;
  for (let j=0;j<nv;j++) { for (let i=0;i<n;i++) mx[j]+=x[i][j]; mx[j]/=n; }
  const XTX = Array.from({length:nv},()=>Array(nv).fill(0));
  const XTy = Array(nv).fill(0);
  for (let i=0;i<n;i++) {
    for (let j=0;j<nv;j++) {
      XTy[j]+=(x[i][j]-mx[j])*(y[i]-my);
      for (let k=0;k<nv;k++) XTX[j][k]+=(x[i][j]-mx[j])*(x[i][k]-mx[k]);
    }
  }
  const beta = solveLinear(XTX, XTy);
  let ssRes=0, ssTot=0;
  for (let i=0;i<n;i++) {
    let pred=my;
    for (let j=0;j<nv;j++) pred+=beta[j]*(x[i][j]-mx[j]);
    gals[i].VfResid=y[i]-pred;
    ssRes+=(y[i]-pred)**2; ssTot+=(y[i]-my)**2;
  }
  return { R2: 1-ssRes/ssTot };
}

console.log('═'.repeat(70));
console.log('PHASE 400C: Realistic Mock — Noise + Sparse Data + Physical Signal');
console.log('═'.repeat(70));
console.log('\nThe question: Can we reproduce the SPARC regime pattern');
console.log('(strengthening at high-V) with realistic observational conditions?');
console.log('');

function generateRealisticGalaxy(logM200, opts) {
  const M200 = Math.pow(10, logM200);
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

  const massFrac = Math.max(0, Math.min(1, (logM200 - 10) / 2.5));

  let physicalShift = 0;
  if (opts.dynIntegration) {
    physicalShift += opts.dynStrength * massFrac * (1 + 0.3 * concScatter);
  }
  if (opts.lowMassNoise) {
    const noiseScale = opts.lowNoiseAmp * (1 - massFrac) + 0.02 * massFrac;
    physicalShift += randNorm(0, noiseScale);
  }

  const nPtsBase = opts.sparseLowMass ?
    Math.round(5 + 20 * massFrac + randNorm(0, 2)) :
    20;
  const nPts = Math.max(4, Math.min(30, nPtsBase));

  const Vobs_noise = opts.velocityNoise || 0.05;
  const lowMassNoiseMult = opts.sparseLowMass ? (1 + 2 * (1 - massFrac)) : 1;
  const effectiveNoise = Vobs_noise * lowMassNoiseMult;

  const rMin = 0.3 * Rd;
  const rMax = (opts.sparseLowMass ? (3 + 3 * massFrac) : 5) * Rd;
  const points = [];
  let Vflat = 0;

  const YstarError = randNorm(0, opts.mlNoise || 0.1);

  for (let i = 0; i < nPts; i++) {
    const r = rMin + (rMax - rMin) * (i / (nPts - 1));
    const MstarEnc = expDiskMenc(r, Mstar, Rd);
    const MgasEnc = gasDiskMenc(r, Mgas * 1.33, Rg);
    const MbarEnc = MstarEnc + MgasEnc;
    const MdmEnc = nfwMenc(r, M200, c);
    const Mtot = MbarEnc + MdmEnc;
    let VobsTrue = Math.sqrt(G * Mtot / r);
    const MbarObs = MbarEnc * Math.pow(10, YstarError);
    const Vbar = Math.sqrt(G * MbarObs / r);
    const gbar = Vbar * Vbar / r;

    if (physicalShift !== 0 && gbar > 0) {
      const a0Eff = 3703 * Math.pow(10, physicalShift);
      const gobsRarTrue = rarPredict(gbar, 3703);
      const gobsRarShifted = rarPredict(gbar, a0Eff);
      if (gobsRarTrue > 0) {
        VobsTrue = VobsTrue * Math.sqrt(gobsRarShifted / gobsRarTrue);
      }
    }

    const Vobs = VobsTrue * (1 + randNorm(0, effectiveNoise));
    if (Math.abs(Vobs) > Vflat) Vflat = Math.abs(Vobs);
    const gobs = Vobs * Math.abs(Vobs) / r;
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
    logC: Math.log10(c), concScatter, physicalShift,
    nPts: points.length, massFrac
  };
}

function runTest(label, opts, NGAL = 500) {
  console.log(`\n${'─'.repeat(60)}`);
  console.log(`${label}`);
  console.log('─'.repeat(60));

  const gals = [];
  for (let i = 0; i < NGAL; i++) {
    const logM200 = rand(10.0, 12.5);
    gals.push(generateRealisticGalaxy(logM200, opts));
  }

  const valid = gals.filter(g => g.rmse < 0.5 && isFinite(g.logA0) && g.Vflat > 15);
  computeVfResid(valid);

  const rAll = pearsonR(valid.map(g=>g.VfResid), valid.map(g=>g.logA0));
  console.log(`  N=${valid.length}, r(VfResid,logA0) = ${rAll.toFixed(3)}`);

  const thresholds = [50, 80, 100, 120, 150, 180];
  const regimeData = [];
  for (const th of thresholds) {
    const sub = valid.filter(g => g.Vflat >= th);
    if (sub.length < 8) continue;
    const r = pearsonR(sub.map(g=>g.VfResid), sub.map(g=>g.logA0));
    regimeData.push({ threshold: th, N: sub.length, r });
    console.log(`  Vflat >= ${th}: N=${sub.length}, r=${r.toFixed(3)}`);
  }

  const lowV = valid.filter(g => g.Vflat < 120);
  const highV = valid.filter(g => g.Vflat >= 120);
  let rLow = null, rHigh = null;
  if (lowV.length >= 5 && highV.length >= 5) {
    rLow = pearsonR(lowV.map(g=>g.VfResid), lowV.map(g=>g.logA0));
    rHigh = pearsonR(highV.map(g=>g.VfResid), highV.map(g=>g.logA0));
    const strengthens = rHigh > rLow + 0.05;
    console.log(`  Regime: low=${rLow.toFixed(3)}, high=${rHigh.toFixed(3)} → ${strengthens ? 'STRENGTHENS ✓' : 'weakens ✗'}`);
  }

  return { label, N: valid.length, rAll, rLow, rHigh, regimeData };
}

const results = [];

results.push(runTest('A: Clean mock (no noise, no physics)', {
  dynIntegration: false, lowMassNoise: false, sparseLowMass: false,
  velocityNoise: 0.02, mlNoise: 0.02
}));

results.push(runTest('B: Realistic noise only (SPARC-like errors)', {
  dynIntegration: false, lowMassNoise: false, sparseLowMass: true,
  velocityNoise: 0.08, mlNoise: 0.12
}));

results.push(runTest('C: Noise + low-mass a₀ scatter', {
  dynIntegration: false, lowMassNoise: true, lowNoiseAmp: 0.35,
  sparseLowMass: true, velocityNoise: 0.08, mlNoise: 0.12
}));

results.push(runTest('D: Noise + dynamical integration (weak)', {
  dynIntegration: true, dynStrength: 0.15, lowMassNoise: true, lowNoiseAmp: 0.2,
  sparseLowMass: true, velocityNoise: 0.08, mlNoise: 0.12
}));

results.push(runTest('E: Noise + dynamical integration (strong)', {
  dynIntegration: true, dynStrength: 0.3, lowMassNoise: true, lowNoiseAmp: 0.35,
  sparseLowMass: true, velocityNoise: 0.08, mlNoise: 0.12
}));

results.push(runTest('F: Noise + strong dynInt + heavy low-mass scatter', {
  dynIntegration: true, dynStrength: 0.4, lowMassNoise: true, lowNoiseAmp: 0.5,
  sparseLowMass: true, velocityNoise: 0.10, mlNoise: 0.15
}));

results.push(runTest('G: SPARC-calibrated (tuned to match observed pattern)', {
  dynIntegration: true, dynStrength: 0.35, lowMassNoise: true, lowNoiseAmp: 0.45,
  sparseLowMass: true, velocityNoise: 0.12, mlNoise: 0.15
}));

console.log('\n\n' + '═'.repeat(70));
console.log('SUMMARY: Regime pattern reproduction');
console.log('═'.repeat(70));

function pad(s, n) { return String(s).padEnd(n); }
console.log(pad('Test', 55) + pad('r(all)', 10) + pad('r(low)', 10) + pad('r(high)', 10) + 'Pattern');
console.log('─'.repeat(95));
console.log(pad('SPARC OBSERVED', 55) + pad('0.70', 10) + pad('0.30', 10) + pad('0.75', 10) + 'STRENGTHENS');
console.log('─'.repeat(95));
for (const r of results) {
  const pattern = r.rLow !== null && r.rHigh !== null ?
    (r.rHigh > r.rLow + 0.05 ? 'STRENGTHENS ✓' : 'weakens ✗') : 'N/A';
  console.log(
    pad(r.label, 55) +
    pad(r.rAll.toFixed(3), 10) +
    pad(r.rLow !== null ? r.rLow.toFixed(3) : 'N/A', 10) +
    pad(r.rHigh !== null ? r.rHigh.toFixed(3) : 'N/A', 10) +
    pattern
  );
}

console.log('\n' + '═'.repeat(70));
console.log('PHASE 400C VERDICT');
console.log('═'.repeat(70));

const testE = results.find(r => r.label.startsWith('E'));
const testF = results.find(r => r.label.startsWith('F'));
const testG = results.find(r => r.label.startsWith('G'));

const anyStrengthen = results.some(r => r.rHigh !== null && r.rLow !== null && r.rHigh > r.rLow + 0.05);

if (anyStrengthen) {
  const match = results.filter(r => r.rHigh !== null && r.rLow !== null && r.rHigh > r.rLow + 0.05);
  console.log('\nRegime strengthening REPRODUCED in:');
  for (const m of match) {
    console.log(`  ${m.label}: low=${m.rLow.toFixed(3)}, high=${m.rHigh.toFixed(3)}`);
  }
  console.log('\nRequired ingredients:');
  console.log('  1. Realistic observational noise (velocity errors, M/L uncertainty)');
  console.log('  2. Sparse low-mass data (fewer points, shorter radial extent)');
  console.log('  3. Physical a₀ variation that strengthens with mass (dynIntegration)');
  console.log('  4. Large intrinsic scatter at low mass (diversity)');
  console.log('\nInterpretation: The SPARC regime pattern is a COMBINATION of:');
  console.log('  - Observational suppression of the concentration signal at low-V');
  console.log('  - Physical amplification of a genuine coupling at high-V');
  console.log('  - Both effects working together create the observed strengthening');
} else {
  console.log('\nRegime strengthening NOT reproduced in any test.');
  console.log('The concentration baseline is too strong even with realistic noise.');
  console.log('This means the SPARC pattern requires QUALITATIVELY different physics');
  console.log('from what simple NFW+disk models can produce.');
  console.log('\nPossible explanations:');
  console.log('  1. Non-circular motions at low-V that scramble the signal');
  console.log('  2. Baryon-dependent halo profiles (not simple NFW)');
  console.log('  3. A genuinely mass-dependent acceleration scale');
  console.log('  4. Selection effects in SPARC sample construction');
}

const outPath = path.join(__dirname, '..', 'public', 'phase400c-noise-physics.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '400C',
  title: 'Realistic Mock with Noise and Physical Signal',
  timestamp: new Date().toISOString(),
  tests: results.map(r => ({
    label: r.label,
    N: r.N,
    rAll: +r.rAll.toFixed(4),
    rLow: r.rLow !== null ? +r.rLow.toFixed(4) : null,
    rHigh: r.rHigh !== null ? +r.rHigh.toFixed(4) : null,
    strengthens: r.rHigh !== null && r.rLow !== null && r.rHigh > r.rLow + 0.05,
    regimeData: r.regimeData
  })),
  sparcComparison: { rAll: 0.70, rLow: 0.30, rHigh: 0.75, pattern: 'strengthens' },
  anyStrengthens: anyStrengthen
}, null, 2));

console.log(`\nResults saved to: ${outPath}`);
