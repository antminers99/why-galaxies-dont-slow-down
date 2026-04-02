const fs = require('fs');
const path = require('path');

const G = 4.3009e-6;
const RNG = (() => { let s = 42; return () => { s = (s * 1664525 + 1013904223) & 0x7fffffff; return s / 0x7fffffff; }; })();
function gaussRNG(mu, sig) { const u1 = RNG(), u2 = RNG(); return mu + sig * Math.sqrt(-2 * Math.log(u1 + 1e-10)) * Math.cos(2 * Math.PI * u2); }

function linearRegression(x, y) {
  const n = x.length;
  if (n < 3) return { slope: NaN, intercept: NaN, r2: NaN, r: NaN, sSlope: NaN, n };
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) { sxx += (x[i] - mx) ** 2; sxy += (x[i] - mx) * (y[i] - my); }
  const slope = sxy / sxx;
  const intercept = my - slope * mx;
  let ssRes = 0, ssTot = 0;
  for (let i = 0; i < n; i++) { ssRes += (y[i] - (intercept + slope * x[i])) ** 2; ssTot += (y[i] - my) ** 2; }
  const r2 = ssTot > 0 ? 1 - ssRes / ssTot : 0;
  const r = Math.sign(slope) * Math.sqrt(Math.abs(r2));
  const sSlope = n > 2 && sxx > 0 ? Math.sqrt(ssRes / (n - 2) / sxx) : NaN;
  return { slope, intercept, r2, r, sSlope, n };
}

function pearsonR(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b) / n;
  const my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxy / Math.sqrt(sxx * syy);
}

function partialR(x, y, z) {
  const rxy = pearsonR(x, y);
  const rxz = pearsonR(x, z);
  const ryz = pearsonR(y, z);
  return (rxy - rxz * ryz) / Math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2));
}

function nfwMenc(r, M200, c) {
  const rho_crit = 136.18;
  const R200 = Math.cbrt(M200 / (4 / 3 * Math.PI * 200 * rho_crit));
  const rs = R200 / c;
  const x = r / rs;
  const gc = Math.log(1 + c) - c / (1 + c);
  return M200 * (Math.log(1 + x) - x / (1 + x)) / gc;
}

function expDiskMenc(r, Mstar, Rd) {
  const x = r / Rd;
  return Mstar * (1 - (1 + x) * Math.exp(-x));
}

function gasDiskMenc(r, Mgas, Rg) {
  const x = r / Rg;
  return Mgas * (1 - (1 + x) * Math.exp(-x));
}

console.log(`\n${'═'.repeat(70)}`);
console.log(`  SIMULATION TEST: NFW + Exponential Disk Mock Galaxies`);
console.log(`  Testing f_DM vs Σ_bar on controlled mock data`);
console.log(`${'═'.repeat(70)}\n`);

const scenarios = [
  {
    name: 'ΛCDM (abundance matching, no scatter)',
    description: 'Standard NFW halos with abundance-matched baryon fraction',
    generate: (n) => {
      const gals = [];
      for (let i = 0; i < n; i++) {
        const logM200 = 10 + RNG() * 2.5;
        const M200 = Math.pow(10, logM200);
        const c = 10 * Math.pow(M200 / 1e12, -0.1) * (0.8 + 0.4 * RNG());
        const fb = 0.05 * Math.pow(M200 / 1e12, 0.3);
        const fstar = 0.6 + 0.2 * RNG();
        const Mstar = M200 * fb * fstar;
        const Mgas = M200 * fb * (1 - fstar) * 1.33;
        const Rd = 0.015 * Math.cbrt(M200 / 1e12) * (0.7 + 0.6 * RNG()) * 1000;
        const Rg = Rd * (1.5 + RNG());
        gals.push({ M200, c, Mstar, Mgas, Rd, Rg, type: logM200 < 11 ? 'dwarf' : logM200 < 12 ? 'spiral' : 'massive' });
      }
      return gals;
    },
  },
  {
    name: 'ΛCDM + baryonic feedback (EAGLE-like)',
    description: 'Halo concentration anti-correlates with baryon density (feedback expands halos)',
    generate: (n) => {
      const gals = [];
      for (let i = 0; i < n; i++) {
        const logM200 = 10 + RNG() * 2.5;
        const M200 = Math.pow(10, logM200);
        const fb = 0.05 * Math.pow(M200 / 1e12, 0.3) * (0.5 + RNG());
        const fstar = 0.5 + 0.3 * RNG();
        const Mstar = M200 * fb * fstar;
        const Mgas = M200 * fb * (1 - fstar) * 1.33;
        const Rd = 0.015 * Math.cbrt(M200 / 1e12) * (0.5 + RNG()) * 1000;
        const Rg = Rd * (1.5 + RNG());
        const SigmaBar = (Mstar + Mgas) / (Math.PI * Rd * Rd);
        const logSig = Math.log10(SigmaBar);
        const cFeedback = 10 * Math.pow(M200 / 1e12, -0.1) * (1 + 0.15 * (logSig - 7));
        const c = Math.max(2, cFeedback * (0.8 + 0.4 * RNG()));
        gals.push({ M200, c, Mstar, Mgas, Rd, Rg, type: logM200 < 11 ? 'dwarf' : logM200 < 12 ? 'spiral' : 'massive' });
      }
      return gals;
    },
  },
  {
    name: 'Independent halos (null hypothesis)',
    description: 'Halo params randomized independently of baryons',
    generate: (n) => {
      const gals = [];
      for (let i = 0; i < n; i++) {
        const logM200 = 10 + RNG() * 2.5;
        const M200 = Math.pow(10, logM200);
        const c = 4 + RNG() * 16;
        const Mstar = Math.pow(10, 7 + RNG() * 4);
        const Mgas = Mstar * (0.5 + RNG() * 2);
        const Rd = 0.5 + RNG() * 8;
        const Rg = Rd * (1.5 + RNG());
        gals.push({ M200, c, Mstar, Mgas, Rd, Rg, type: logM200 < 11 ? 'dwarf' : logM200 < 12 ? 'spiral' : 'massive' });
      }
      return gals;
    },
  },
];

const allResults = [];

for (const scenario of scenarios) {
  console.log(`\n╔══════════════════════════════════════════════════════════════╗`);
  console.log(`║  ${scenario.name.padEnd(58)} ║`);
  console.log(`╚══════════════════════════════════════════════════════════════╝`);
  console.log(`  ${scenario.description}\n`);

  const mockGals = scenario.generate(300);
  const galaxyData = [];
  const allPoints = [];

  for (const gal of mockGals) {
    const radii = [];
    for (let i = 1; i <= 25; i++) radii.push(gal.Rd * 0.2 * i);
    const rmax = radii[radii.length - 1];

    const pts = [];
    for (const r of radii) {
      const Mbar_enc = expDiskMenc(r, gal.Mstar, gal.Rd) + gasDiskMenc(r, gal.Mgas, gal.Rg);
      const Mdm_enc = nfwMenc(r, gal.M200, gal.c);
      const Mtot = Mbar_enc + Mdm_enc;

      const vObs = Math.sqrt(G * Mtot / r);
      const vBar = Math.sqrt(G * Mbar_enc / r);

      const gObs = vObs * vObs / r;
      const gBar = vBar * vBar / r;

      if (gObs <= 0 || gBar <= 0) continue;
      const fDM = (gObs - gBar) / gObs;
      if (fDM < 0 || fDM > 1) continue;

      const sigBarEnc = Mbar_enc / (Math.PI * r * r);
      const sigBarLocal = (gal.Mstar * Math.exp(-r / gal.Rd) / (2 * Math.PI * gal.Rd * gal.Rd) +
                            gal.Mgas * Math.exp(-r / gal.Rg) / (2 * Math.PI * gal.Rg * gal.Rg));

      pts.push({ r, fDM, logSigEnc: Math.log10(sigBarEnc), logSigLocal: Math.log10(sigBarLocal + 1e-10), rNorm: r / rmax, vObs });
    }

    if (pts.length < 5) continue;

    const vmax = pts.reduce((mx, p) => Math.max(mx, p.vObs), 0);
    const sigBarGal = (gal.Mstar + gal.Mgas) / (Math.PI * gal.Rd * gal.Rd);
    const logSigGal = Math.log10(sigBarGal);

    const meanFDM = pts.reduce((s, p) => s + p.fDM, 0) / pts.length;

    const sigEncHalf = pts.filter(p => p.rNorm <= 0.55 && p.rNorm >= 0.45);
    const sigAtHalf = sigEncHalf.length > 0 ? sigEncHalf[0].logSigEnc : logSigGal;

    galaxyData.push({
      name: `mock_${galaxyData.length}`, vmax, logVmax: Math.log10(vmax),
      meanFDM, logSigGal, logSigAtHalf: sigAtHalf, type: gal.type,
      M200: gal.M200, c: gal.c, Mstar: gal.Mstar, Rd: gal.Rd,
    });

    for (const pt of pts) {
      allPoints.push({ ...pt, galaxy: `mock_${galaxyData.length - 1}`, vmax, type: gal.type });
    }
  }

  console.log(`  Mock galaxies: ${galaxyData.length}, Radial points: ${allPoints.length}\n`);

  const regPtEnc = linearRegression(allPoints.map(p => p.logSigEnc), allPoints.map(p => p.fDM));
  const regPtLocal = linearRegression(allPoints.map(p => p.logSigLocal), allPoints.map(p => p.fDM));

  const regGal = linearRegression(galaxyData.map(g => g.logSigGal), galaxyData.map(g => g.meanFDM));
  const partRGal = partialR(galaxyData.map(g => g.logSigGal), galaxyData.map(g => g.meanFDM), galaxyData.map(g => g.logVmax));

  console.log(`  Point-level (Σ_enc):   b=${regPtEnc.slope.toFixed(5)}, r=${regPtEnc.r.toFixed(4)}, R²=${regPtEnc.r2.toFixed(4)}`);
  console.log(`  Point-level (Σ_local): b=${regPtLocal.slope.toFixed(5)}, r=${regPtLocal.r.toFixed(4)}, R²=${regPtLocal.r2.toFixed(4)}`);
  console.log(`  Per-galaxy:            b=${regGal.slope.toFixed(5)}, r=${regGal.r.toFixed(4)}, partial r|V=${partRGal.toFixed(4)}\n`);

  const typeResults = [];
  for (const type of ['dwarf', 'spiral', 'massive']) {
    const gals = galaxyData.filter(g => g.type === type);
    if (gals.length < 10) continue;
    const reg = linearRegression(gals.map(g => g.logSigGal), gals.map(g => g.meanFDM));
    console.log(`  ${type.padEnd(10)} (${gals.length}): b=${reg.slope.toFixed(5)}, r=${reg.r.toFixed(4)}`);
    typeResults.push({ type, n: gals.length, slope: reg.slope, r: reg.r });
  }

  const innerPts = allPoints.filter(p => p.rNorm < 0.5);
  const outerPts = allPoints.filter(p => p.rNorm >= 0.5);
  const regInner = linearRegression(innerPts.map(p => p.logSigEnc), innerPts.map(p => p.fDM));
  const regOuter = linearRegression(outerPts.map(p => p.logSigEnc), outerPts.map(p => p.fDM));
  console.log(`\n  Inner (r<0.5): b=${regInner.slope.toFixed(5)}, r=${regInner.r.toFixed(4)}`);
  console.log(`  Outer (r≥0.5): b=${regOuter.slope.toFixed(5)}, r=${regOuter.r.toFixed(4)}`);

  const samplePts = [];
  for (let i = 0; i < allPoints.length; i += Math.max(1, Math.floor(allPoints.length / 400))) {
    samplePts.push({ logSigEnc: +allPoints[i].logSigEnc.toFixed(3), fDM: +allPoints[i].fDM.toFixed(4), type: allPoints[i].type });
  }

  allResults.push({
    name: scenario.name,
    description: scenario.description,
    nGalaxies: galaxyData.length,
    nPoints: allPoints.length,
    pointLevel: { encSlope: regPtEnc.slope, encR: regPtEnc.r, encR2: regPtEnc.r2, localSlope: regPtLocal.slope, localR: regPtLocal.r },
    perGalaxy: { slope: regGal.slope, r: regGal.r, r2: regGal.r2, partialR: partRGal, n: regGal.n },
    byType: typeResults,
    innerOuter: {
      inner: { slope: regInner.slope, r: regInner.r, n: regInner.n },
      outer: { slope: regOuter.slope, r: regOuter.r, n: regOuter.n },
    },
    samplePoints: samplePts,
    galaxies: galaxyData.map(g => ({ meanFDM: +g.meanFDM.toFixed(4), logSigGal: +g.logSigGal.toFixed(3), vmax: +g.vmax.toFixed(1), type: g.type })),
  });
}

console.log(`\n${'═'.repeat(70)}`);
console.log(`  COMPARISON: SIMULATION vs OBSERVATION`);
console.log(`${'═'.repeat(70)}`);
console.log(`\n  ${'Scenario'.padEnd(45)} ${'slope'.padStart(8)} ${'r'.padStart(8)} ${'pr|V'.padStart(8)}`);
console.log(`  ${'-'.repeat(71)}`);

const sparcData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'fdm-analysis.json'), 'utf8'));
const sparcSlope = sparcData.sparc.perGalaxy.slope;
const sparcR = sparcData.sparc.perGalaxy.r;
const sparcPR = sparcData.sparc.perGalaxy.partialR;
console.log(`  ${'SPARC (OBSERVED)'.padEnd(45)} ${sparcSlope.toFixed(4).padStart(8)} ${sparcR.toFixed(4).padStart(8)} ${sparcPR.toFixed(4).padStart(8)}`);
for (const r of allResults) {
  console.log(`  ${r.name.padEnd(45)} ${r.perGalaxy.slope.toFixed(4).padStart(8)} ${r.perGalaxy.r.toFixed(4).padStart(8)} ${r.perGalaxy.partialR.toFixed(4).padStart(8)}`);
}

console.log(`\n  Key question: Does standard ΛCDM produce the observed slope?`);
const cdmSlope = allResults[0].perGalaxy.slope;
const obsSlope = sparcSlope;
console.log(`  ΛCDM slope:    ${cdmSlope.toFixed(4)}`);
console.log(`  Observed slope: ${obsSlope.toFixed(4)}`);
console.log(`  Ratio obs/sim:  ${(obsSlope / cdmSlope).toFixed(2)}`);
console.log(`  → ${Math.abs(obsSlope) > Math.abs(cdmSlope) * 1.5 ? 'Observed effect is STRONGER than ΛCDM predicts' : 'Comparable to ΛCDM'}`);

sparcData.simulation = allResults;
fs.writeFileSync(path.join(__dirname, '..', 'public', 'fdm-analysis.json'), JSON.stringify(sparcData, null, 2));
console.log(`\nSaved simulation results to fdm-analysis.json`);
