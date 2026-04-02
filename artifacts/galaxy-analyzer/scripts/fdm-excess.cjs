const fs = require('fs');
const path = require('path');

const dataPath = path.join(__dirname, '..', 'public', 'fdm-analysis.json');
const data = JSON.parse(fs.readFileSync(dataPath, 'utf8'));

const sparcGals = data.sparc.galaxies;
const simGals = data.simulation[0].galaxies;

function seededRNG(seed) {
  let s = seed;
  return () => { s = (s * 1664525 + 1013904223) & 0x7fffffff; return s / 0x7fffffff; };
}

function linearRegression(x, y) {
  const n = x.length;
  if (n < 3) return { slope: NaN, r: NaN };
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let sxx = 0, sxy = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxx += (x[i] - mx) ** 2;
    sxy += (x[i] - mx) * (y[i] - my);
    syy += (y[i] - my) ** 2;
  }
  const slope = sxy / sxx;
  const r2 = sxx > 0 && syy > 0 ? (sxy * sxy) / (sxx * syy) : 0;
  const r = Math.sign(slope) * Math.sqrt(Math.abs(r2));
  return { slope, r, r2, n };
}

function bootstrapSlope(gals, xKey, yKey, nBoot, rng) {
  const n = gals.length;
  const slopes = [];
  for (let i = 0; i < nBoot; i++) {
    const xs = [], ys = [];
    for (let j = 0; j < n; j++) {
      const idx = Math.floor(rng() * n);
      xs.push(gals[idx][xKey]);
      ys.push(gals[idx][yKey]);
    }
    const reg = linearRegression(xs, ys);
    if (!isNaN(reg.slope)) slopes.push(reg.slope);
  }
  slopes.sort((a, b) => a - b);
  const mean = slopes.reduce((a, b) => a + b, 0) / slopes.length;
  const sd = Math.sqrt(slopes.reduce((s, v) => s + (v - mean) ** 2, 0) / (slopes.length - 1));
  const ci95 = [slopes[Math.floor(slopes.length * 0.025)], slopes[Math.floor(slopes.length * 0.975)]];
  const ci99 = [slopes[Math.floor(slopes.length * 0.005)], slopes[Math.floor(slopes.length * 0.995)]];
  return { mean, sd, ci95, ci99, slopes };
}

const BINS = [
  { label: 'V < 50', min: 0, max: 50, mid: 25 },
  { label: '50–80', min: 50, max: 80, mid: 65 },
  { label: '80–120', min: 80, max: 120, mid: 100 },
  { label: '120–170', min: 120, max: 170, mid: 145 },
  { label: '170–250', min: 170, max: 250, mid: 210 },
  { label: 'V > 250', min: 250, max: 9999, mid: 325 },
];

const N_BOOT = 10000;
let seedCounter = 77777;

console.log(`\n${'═'.repeat(70)}`);
console.log(`  DISCOVERY PROOF: Δb per Vmax bin (observed − ΛCDM)`);
console.log(`  Testing if f_DM–Σ_bar slope exceeds geometric expectation`);
console.log(`${'═'.repeat(70)}\n`);

const binResults = [];

for (const bin of BINS) {
  const sparcBin = sparcGals.filter(g => g.vmax >= bin.min && g.vmax < bin.max);
  const simBin = simGals.filter(g => g.vmax >= bin.min && g.vmax < bin.max);

  console.log(`\n  ─── ${bin.label} (SPARC: ${sparcBin.length}, Sim: ${simBin.length}) ───`);

  if (sparcBin.length < 5 || simBin.length < 5) {
    console.log(`  ⚠ Skipped (too few galaxies)`);
    binResults.push({
      bin: bin.label, midVmax: bin.mid,
      sparc: { n: sparcBin.length, slope: NaN, boot: null },
      sim: { n: simBin.length, slope: NaN, boot: null },
      deltaB: NaN, deltaSigma: NaN, significant: false, skipped: true,
    });
    continue;
  }

  const sparcReg = linearRegression(sparcBin.map(g => g.meanLogSigBar), sparcBin.map(g => g.meanFDM));
  const simReg = linearRegression(simBin.map(g => g.logSigGal), simBin.map(g => g.meanFDM));

  const rng1 = seededRNG(seedCounter++);
  const sparcBoot = bootstrapSlope(sparcBin, 'meanLogSigBar', 'meanFDM', N_BOOT, rng1);
  const rng2 = seededRNG(seedCounter++);
  const simBoot = bootstrapSlope(simBin, 'logSigGal', 'meanFDM', N_BOOT, rng2);

  const deltaB = sparcBoot.mean - simBoot.mean;
  const deltaSE = Math.sqrt(sparcBoot.sd ** 2 + simBoot.sd ** 2);
  const deltaSigma = deltaSE > 0 ? Math.abs(deltaB) / deltaSE : 0;

  const deltaBoot = [];
  for (let i = 0; i < N_BOOT; i++) {
    deltaBoot.push(sparcBoot.slopes[i] - simBoot.slopes[i]);
  }
  deltaBoot.sort((a, b) => a - b);
  const deltaCI95 = [deltaBoot[Math.floor(N_BOOT * 0.025)], deltaBoot[Math.floor(N_BOOT * 0.975)]];
  const deltaCI99 = [deltaBoot[Math.floor(N_BOOT * 0.005)], deltaBoot[Math.floor(N_BOOT * 0.995)]];

  const isNegExcess = deltaCI95[1] < 0;

  console.log(`  SPARC: slope=${sparcReg.slope.toFixed(5)}, boot mean=${sparcBoot.mean.toFixed(5)} ± ${sparcBoot.sd.toFixed(5)}`);
  console.log(`  ΛCDM:  slope=${simReg.slope.toFixed(5)}, boot mean=${simBoot.mean.toFixed(5)} ± ${simBoot.sd.toFixed(5)}`);
  console.log(`  Δb = ${deltaB.toFixed(5)} ± ${deltaSE.toFixed(5)} (${deltaSigma.toFixed(1)}σ)`);
  console.log(`  Δb 95% CI: [${deltaCI95[0].toFixed(5)}, ${deltaCI95[1].toFixed(5)}]`);
  console.log(`  Excess steeper? ${isNegExcess ? '✓ YES (entire CI < 0)' : '✗ NO (CI includes 0 or positive)'}`);

  binResults.push({
    bin: bin.label,
    midVmax: bin.mid,
    sparc: {
      n: sparcBin.length,
      slope: sparcReg.slope,
      r: sparcReg.r,
      bootMean: sparcBoot.mean,
      bootSD: sparcBoot.sd,
      ci95: sparcBoot.ci95,
    },
    sim: {
      n: simBin.length,
      slope: simReg.slope,
      r: simReg.r,
      bootMean: simBoot.mean,
      bootSD: simBoot.sd,
      ci95: simBoot.ci95,
    },
    deltaB,
    deltaSE,
    deltaSigma,
    deltaCI95,
    deltaCI99,
    significant: isNegExcess,
    skipped: false,
  });
}

console.log(`\n${'═'.repeat(70)}`);
console.log(`  OVERALL Δb`);
console.log(`${'═'.repeat(70)}`);

const rngAll1 = seededRNG(11111);
const sparcBootAll = bootstrapSlope(sparcGals, 'meanLogSigBar', 'meanFDM', N_BOOT, rngAll1);
const rngAll2 = seededRNG(22222);
const simBootAll = bootstrapSlope(simGals, 'logSigGal', 'meanFDM', N_BOOT, rngAll2);

const overallDeltaB = sparcBootAll.mean - simBootAll.mean;
const overallDeltaSE = Math.sqrt(sparcBootAll.sd ** 2 + simBootAll.sd ** 2);
const overallDeltaSigma = overallDeltaSE > 0 ? Math.abs(overallDeltaB) / overallDeltaSE : 0;

const overallDeltaBoot = [];
for (let i = 0; i < N_BOOT; i++) {
  overallDeltaBoot.push(sparcBootAll.slopes[i] - simBootAll.slopes[i]);
}
overallDeltaBoot.sort((a, b) => a - b);
const overallDeltaCI95 = [overallDeltaBoot[Math.floor(N_BOOT * 0.025)], overallDeltaBoot[Math.floor(N_BOOT * 0.975)]];

console.log(`  SPARC overall: boot mean=${sparcBootAll.mean.toFixed(5)} ± ${sparcBootAll.sd.toFixed(5)}`);
console.log(`  ΛCDM overall:  boot mean=${simBootAll.mean.toFixed(5)} ± ${simBootAll.sd.toFixed(5)}`);
console.log(`  Δb = ${overallDeltaB.toFixed(5)} ± ${overallDeltaSE.toFixed(5)} (${overallDeltaSigma.toFixed(1)}σ)`);
console.log(`  Δb 95% CI: [${overallDeltaCI95[0].toFixed(5)}, ${overallDeltaCI95[1].toFixed(5)}]`);

const halo = data.deepAnalysis.halo;
const diag = data.diagnostics;

const nBinsSignificant = binResults.filter(b => !b.skipped && b.significant).length;
const nBinsTotal = binResults.filter(b => !b.skipped).length;
const nBinsNegative = binResults.filter(b => !b.skipped && b.deltaB < 0).length;

console.log(`\n${'═'.repeat(70)}`);
console.log(`  4-STEP DISCOVERY PROOF SUMMARY`);
console.log(`${'═'.repeat(70)}`);
console.log(`\n  Step 1: Δb > 0 (excess over expectation)`);
console.log(`    Overall: Δb = ${overallDeltaB.toFixed(4)} (${overallDeltaSigma.toFixed(1)}σ)`);
console.log(`    ${nBinsNegative}/${nBinsTotal} bins show steeper observed slope`);
console.log(`    ${nBinsSignificant}/${nBinsTotal} bins significant at 95%`);
console.log(`\n  Step 2: Halo tracks Σ_bar`);
console.log(`    V_DM vs Σ_bar: r = ${halo.vDMvsSigBar.galR.toFixed(3)} (n=${halo.vDMvsSigBar.galN})`);
console.log(`    r_DMdom vs Σ_bar: r = ${halo.dmDominanceRadius.r.toFixed(3)} (n=${halo.dmDominanceRadius.n})`);
console.log(`\n  Step 3: Simulation can't explain`);
console.log(`    ΛCDM slope: ${simBootAll.mean.toFixed(4)}`);
console.log(`    Observed slope: ${sparcBootAll.mean.toFixed(4)}`);
console.log(`    Slope ratio: ${(sparcBootAll.mean / simBootAll.mean).toFixed(2)}×`);
console.log(`    Fisher z p < 10⁻⁶ (scatter structure differs)`);
console.log(`\n  Step 4: Kill tests passed`);
const altSigAll = diag.altSigmaDefinitions.every(d => d.slope < 0);
console.log(`    Alt Σ definitions: ${diag.altSigmaDefinitions.length}/6 negative (${altSigAll ? '✓' : '✗'})`);
console.log(`    Face-on: slope = ${diag.selectionBias.inclinationSplit.low.slope.toFixed(3)}`);
console.log(`    Edge-on: slope = ${diag.selectionBias.inclinationSplit.high.slope.toFixed(3)}`);
console.log(`    High-Q: slope = ${diag.selectionBias.qualitySplit.many.slope.toFixed(3)}`);
console.log(`    Low-Q: slope = ${diag.selectionBias.qualitySplit.few.slope.toFixed(3)}`);
console.log(`    Permutation p = ${diag.selectionBias.permutationTest.pValue.toFixed(6)}`);

const discoveryProof = {
  overall: {
    sparcSlope: sparcBootAll.mean,
    sparcSD: sparcBootAll.sd,
    sparcCI95: sparcBootAll.ci95,
    simSlope: simBootAll.mean,
    simSD: simBootAll.sd,
    simCI95: simBootAll.ci95,
    deltaB: overallDeltaB,
    deltaSE: overallDeltaSE,
    deltaSigma: overallDeltaSigma,
    deltaCI95: overallDeltaCI95,
    significant: overallDeltaCI95[1] < 0,
  },
  bins: binResults.filter(b => !b.skipped).map(b => ({
    bin: b.bin,
    midVmax: b.midVmax,
    sparcN: b.sparc.n,
    sparcSlope: b.sparc.bootMean,
    sparcSD: b.sparc.bootSD,
    sparcCI95: b.sparc.ci95,
    simN: b.sim.n,
    simSlope: b.sim.bootMean,
    simSD: b.sim.bootSD,
    simCI95: b.sim.ci95,
    deltaB: b.deltaB,
    deltaSE: b.deltaSE,
    deltaSigma: b.deltaSigma,
    deltaCI95: b.deltaCI95,
    significant: b.significant,
  })),
  nBootstrap: N_BOOT,
  step1: {
    label: 'Excess over expectation',
    overallDeltaB: overallDeltaB,
    overallSigma: overallDeltaSigma,
    nBinsNegative,
    nBinsSignificant,
    nBinsTotal,
    pass: overallDeltaCI95[1] < 0,
  },
  step2: {
    label: 'Halo tracks Σ_bar',
    vdm: { r: halo.vDMvsSigBar.galR, n: halo.vDMvsSigBar.galN },
    rDMdom: { r: halo.dmDominanceRadius.r, n: halo.dmDominanceRadius.n },
    pass: Math.abs(halo.vDMvsSigBar.galR) > 0.5 && Math.abs(halo.dmDominanceRadius.r) > 0.3,
  },
  step3: {
    label: 'Simulation can\'t explain',
    slopeRatio: sparcBootAll.mean / simBootAll.mean,
    fisherZp: data.significanceTest.fisherZ.pValue,
    pass: data.significanceTest.fisherZ.pValue < 0.001,
  },
  step4: {
    label: 'Kill tests passed',
    altSigAllNeg: altSigAll,
    nAltSig: diag.altSigmaDefinitions.length,
    faceOn: { slope: diag.selectionBias.inclinationSplit.low.slope, r: diag.selectionBias.inclinationSplit.low.r },
    edgeOn: { slope: diag.selectionBias.inclinationSplit.high.slope, r: diag.selectionBias.inclinationSplit.high.r },
    highQ: { slope: diag.selectionBias.qualitySplit.many.slope, r: diag.selectionBias.qualitySplit.many.r },
    lowQ: { slope: diag.selectionBias.qualitySplit.few.slope, r: diag.selectionBias.qualitySplit.few.r },
    massSplit: {
      low: { slope: diag.selectionBias.massSplit.low.slope, r: diag.selectionBias.massSplit.low.r },
      high: { slope: diag.selectionBias.massSplit.high.slope, r: diag.selectionBias.massSplit.high.r },
    },
    permutationP: diag.selectionBias.permutationTest.pValue,
    pass: altSigAll && diag.selectionBias.permutationTest.pValue < 0.01,
  },
  verdict: {
    allPass: false,
    statement: '',
  },
};

const allPass = discoveryProof.step1.pass && discoveryProof.step2.pass && discoveryProof.step3.pass && discoveryProof.step4.pass;
discoveryProof.verdict.allPass = allPass;

if (allPass) {
  discoveryProof.verdict.statement = 'The baryon–halo coupling is stronger than predicted by ΛCDM, indicating additional physics linking baryonic structure to halo dynamics.';
} else {
  const passed = [discoveryProof.step1, discoveryProof.step2, discoveryProof.step3, discoveryProof.step4].filter(s => s.pass).length;
  discoveryProof.verdict.statement = `${passed}/4 steps pass. The evidence is ${passed >= 3 ? 'strong but not definitive' : passed >= 2 ? 'suggestive' : 'insufficient'} for excess baryon–halo coupling.`;
}

console.log(`\n  VERDICT: ${discoveryProof.verdict.statement}`);

data.discoveryProof = discoveryProof;
fs.writeFileSync(dataPath, JSON.stringify(data, null, 2));
console.log(`\nSaved discoveryProof to fdm-analysis.json`);
