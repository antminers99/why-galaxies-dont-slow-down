const fs = require('fs');
const path = require('path');

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
  if (n < 3) return NaN;
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
  const denom = Math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2));
  return denom > 0 ? (rxy - rxz * ryz) / denom : NaN;
}

function seededRNG(seed) {
  let s = seed;
  return () => { s = (s * 1664525 + 1013904223) & 0x7fffffff; return s / 0x7fffffff; };
}

function bootstrapSlope(xs, ys, nBoot, rng) {
  const n = xs.length;
  const slopes = [];
  for (let i = 0; i < nBoot; i++) {
    const bx = [], by = [];
    for (let j = 0; j < n; j++) {
      const idx = Math.floor(rng() * n);
      bx.push(xs[idx]); by.push(ys[idx]);
    }
    const reg = linearRegression(bx, by);
    if (!isNaN(reg.slope)) slopes.push(reg.slope);
  }
  slopes.sort((a, b) => a - b);
  const mean = slopes.reduce((a, b) => a + b, 0) / slopes.length;
  const sd = Math.sqrt(slopes.reduce((s, v) => s + (v - mean) ** 2, 0) / (slopes.length - 1));
  return { mean, sd, ci95: [slopes[Math.floor(slopes.length * 0.025)], slopes[Math.floor(slopes.length * 0.975)]] };
}

const rotmodDir = '/tmp/rotmod';
const sparcTablePath = '/tmp/sparc_table.mrt';
const UPSILON_D = 0.5, UPSILON_B = 0.7;
const N_BOOT = 5000;

const sparcTable = {};
const tableLines = fs.readFileSync(sparcTablePath, 'utf8').split('\n');
for (const line of tableLines) {
  const parts = line.trim().split(/\s+/);
  if (parts.length < 7) continue;
  const name = parts[0];
  if (name === 'Galaxy' || name.startsWith('#') || name.startsWith('-')) continue;
  sparcTable[name] = { dist: parseFloat(parts[2]), inc: parseFloat(parts[5]) };
}

const sparcReal = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));

const allPoints = [];
const galaxies = [];
const rotmodFiles = fs.readdirSync(rotmodDir).filter(f => f.endsWith('_rotmod.dat'));

for (const file of rotmodFiles) {
  const name = file.replace('_rotmod.dat', '');
  const info = sparcTable[name];
  if (!info) continue;

  const perGalaxy = sparcReal.perGalaxy?.find((g) => g.name === name);
  if (!perGalaxy || perGalaxy.sigma_bar <= 0) continue;

  const lines = fs.readFileSync(path.join(rotmodDir, file), 'utf8').trim().split('\n');
  const dataLines = lines.filter(l => !l.startsWith('#') && l.trim());

  const pts = [];
  for (const line of dataLines) {
    const p = line.trim().split(/\s+/).map(Number);
    if (p.length < 7) continue;
    const [r, vObs, eV, vGas, vDisk, vBulge] = [p[0], p[1], p[2], p[3], p[4], p[5]];
    if (r <= 0 || vObs <= 0) continue;
    const vBarSq = vGas * Math.abs(vGas) + UPSILON_D * vDisk * Math.abs(vDisk) + UPSILON_B * vBulge * Math.abs(vBulge);
    if (vBarSq <= 0) continue;
    const gObs = vObs * vObs / r;
    const gBar = vBarSq / r;
    const fDM = (gObs - gBar) / gObs;
    const vDMsq = vObs * vObs - vBarSq;
    if (fDM < 0 || fDM > 1 || gObs <= 0 || gBar <= 0) continue;
    pts.push({ r, vObs, vBar: Math.sqrt(vBarSq), fDM, vDMsq, vDM: vDMsq > 0 ? Math.sqrt(vDMsq) : 0 });
  }

  if (pts.length < 5) continue;

  const rmax = pts[pts.length - 1].r;
  const vmax = pts.reduce((mx, p) => Math.max(mx, p.vObs), 0);
  const logVmax = Math.log10(vmax);
  const sigBarGal = perGalaxy.sigma_bar;
  const logSigBar = Math.log10(sigBarGal);

  let rDMdom = NaN;
  for (const pt of pts) {
    if (pt.vDMsq > pt.vBar * pt.vBar && isNaN(rDMdom)) rDMdom = pt.r;
  }

  const innerPts = pts.filter(p => p.r / rmax < 0.3 && p.vDMsq > 0);
  let alpha = NaN;
  if (innerPts.length >= 3) {
    const reg = linearRegression(innerPts.map(p => Math.log10(p.r)), innerPts.map(p => Math.log10(Math.sqrt(p.vDMsq))));
    alpha = reg.slope;
  }

  const meanVDM = pts.filter(p => p.vDMsq > 0).reduce((s, p) => s + Math.sqrt(p.vDMsq), 0) / Math.max(1, pts.filter(p => p.vDMsq > 0).length);
  const meanFDM = pts.reduce((s, p) => s + p.fDM, 0) / pts.length;

  const innerMeanVDM = pts.filter(p => p.r / rmax < 0.5 && p.vDMsq > 0);
  const outerMeanVDM = pts.filter(p => p.r / rmax >= 0.5 && p.vDMsq > 0);

  galaxies.push({
    name, vmax, logVmax, logSigBar, meanFDM,
    meanVDM, logMeanVDM: meanVDM > 0 ? Math.log10(meanVDM) : NaN,
    rDMdomNorm: isNaN(rDMdom) ? null : rDMdom / rmax,
    alpha: isNaN(alpha) ? null : alpha,
    innerMeanVDM: innerMeanVDM.length >= 3 ? innerMeanVDM.reduce((s, p) => s + Math.sqrt(p.vDMsq), 0) / innerMeanVDM.length : null,
    outerMeanVDM: outerMeanVDM.length >= 3 ? outerMeanVDM.reduce((s, p) => s + Math.sqrt(p.vDMsq), 0) / outerMeanVDM.length : null,
    innerFDM: pts.filter(p => p.r / rmax < 0.5).length >= 3 ? pts.filter(p => p.r / rmax < 0.5).reduce((s, p) => s + p.fDM, 0) / pts.filter(p => p.r / rmax < 0.5).length : null,
    outerFDM: pts.filter(p => p.r / rmax >= 0.5).length >= 3 ? pts.filter(p => p.r / rmax >= 0.5).reduce((s, p) => s + p.fDM, 0) / pts.filter(p => p.r / rmax >= 0.5).length : null,
  });
}

console.log(`\n${'═'.repeat(70)}`);
console.log(`  BARYON–HALO COUPLING ANALYSIS`);
console.log(`  ${galaxies.length} galaxies from SPARC`);
console.log(`${'═'.repeat(70)}\n`);

console.log(`╔══════════════════════════════════════════════════════════════╗`);
console.log(`║  PHASE 1: COUPLING FINGERPRINT (raw + partial r)           ║`);
console.log(`╚══════════════════════════════════════════════════════════════╝\n`);

const vdmGals = galaxies.filter(g => g.logMeanVDM && !isNaN(g.logMeanVDM));
const rdomGals = galaxies.filter(g => g.rDMdomNorm !== null);
const alphaGals = galaxies.filter(g => g.alpha !== null);

const metrics = [
  { name: '⟨V_DM⟩', gals: vdmGals, xKey: 'logSigBar', yKey: 'logMeanVDM', zKey: 'logVmax', expectedSign: '+' },
  { name: 'r_DMdom/r_max', gals: rdomGals, xKey: 'logSigBar', yKey: 'rDMdomNorm', zKey: 'logVmax', expectedSign: '+' },
  { name: 'inner α', gals: alphaGals, xKey: 'logSigBar', yKey: 'alpha', zKey: 'logVmax', expectedSign: '−' },
];

const fingerprint = [];

for (const m of metrics) {
  const xs = m.gals.map(g => g[m.xKey]);
  const ys = m.gals.map(g => g[m.yKey]);
  const zs = m.gals.map(g => g[m.zKey]);

  const reg = linearRegression(xs, ys);
  const pr = partialR(xs, ys, zs);

  const rng = seededRNG(fingerprint.length * 11111 + 42);
  const boot = bootstrapSlope(xs, ys, N_BOOT, rng);

  console.log(`  ${m.name} vs log(Σ_bar):`);
  console.log(`    slope = ${reg.slope.toFixed(5)}, r = ${reg.r.toFixed(4)}, n = ${m.gals.length}`);
  console.log(`    partial r|Vmax = ${pr.toFixed(4)}`);
  console.log(`    bootstrap: ${boot.mean.toFixed(5)} [${boot.ci95[0].toFixed(5)}, ${boot.ci95[1].toFixed(5)}]`);
  console.log(`    Expected sign: ${m.expectedSign}, Actual: ${reg.slope > 0 ? '+' : '−'} ${m.expectedSign === (reg.slope > 0 ? '+' : '−') ? '✓' : '✗'}\n`);

  fingerprint.push({
    metric: m.name,
    n: m.gals.length,
    slope: reg.slope,
    slopeErr: reg.sSlope,
    r: reg.r,
    r2: reg.r2,
    partialR: pr,
    bootMean: boot.mean,
    bootCI95: boot.ci95,
    expectedSign: m.expectedSign,
    matchesExpected: m.expectedSign === (reg.slope > 0 ? '+' : '−'),
  });
}

console.log(`╔══════════════════════════════════════════════════════════════╗`);
console.log(`║  PHASE 2: INNER vs OUTER COUPLING                         ║`);
console.log(`╚══════════════════════════════════════════════════════════════╝\n`);

const innerVDMgals = galaxies.filter(g => g.innerMeanVDM !== null && g.innerMeanVDM > 0);
const outerVDMgals = galaxies.filter(g => g.outerMeanVDM !== null && g.outerMeanVDM > 0);
const innerFDMgals = galaxies.filter(g => g.innerFDM !== null);
const outerFDMgals = galaxies.filter(g => g.outerFDM !== null);

const innerVDMreg = linearRegression(innerVDMgals.map(g => g.logSigBar), innerVDMgals.map(g => Math.log10(g.innerMeanVDM)));
const outerVDMreg = linearRegression(outerVDMgals.map(g => g.logSigBar), outerVDMgals.map(g => Math.log10(g.outerMeanVDM)));
const innerFDMreg = linearRegression(innerFDMgals.map(g => g.logSigBar), innerFDMgals.map(g => g.innerFDM));
const outerFDMreg = linearRegression(outerFDMgals.map(g => g.logSigBar), outerFDMgals.map(g => g.outerFDM));

const innerVDMpr = partialR(innerVDMgals.map(g => g.logSigBar), innerVDMgals.map(g => Math.log10(g.innerMeanVDM)), innerVDMgals.map(g => g.logVmax));
const outerVDMpr = partialR(outerVDMgals.map(g => g.logSigBar), outerVDMgals.map(g => Math.log10(g.outerMeanVDM)), outerVDMgals.map(g => g.logVmax));

console.log(`  Inner V_DM vs Σ: slope=${innerVDMreg.slope.toFixed(5)}, r=${innerVDMreg.r.toFixed(4)}, pr|V=${innerVDMpr.toFixed(4)}, n=${innerVDMreg.n}`);
console.log(`  Outer V_DM vs Σ: slope=${outerVDMreg.slope.toFixed(5)}, r=${outerVDMreg.r.toFixed(4)}, pr|V=${outerVDMpr.toFixed(4)}, n=${outerVDMreg.n}`);
console.log(`  Inner f_DM vs Σ: slope=${innerFDMreg.slope.toFixed(5)}, r=${innerFDMreg.r.toFixed(4)}, n=${innerFDMreg.n}`);
console.log(`  Outer f_DM vs Σ: slope=${outerFDMreg.slope.toFixed(5)}, r=${outerFDMreg.r.toFixed(4)}, n=${outerFDMreg.n}`);
console.log(`  → Inner coupling stronger? V_DM: ${Math.abs(innerVDMreg.r) > Math.abs(outerVDMreg.r) ? 'YES' : 'NO'}, f_DM: ${Math.abs(innerFDMreg.r) > Math.abs(outerFDMreg.r) ? 'YES' : 'NO'}\n`);

const innerOuter = {
  innerVDM: { slope: innerVDMreg.slope, r: innerVDMreg.r, partialR: innerVDMpr, n: innerVDMreg.n },
  outerVDM: { slope: outerVDMreg.slope, r: outerVDMreg.r, partialR: outerVDMpr, n: outerVDMreg.n },
  innerFDM: { slope: innerFDMreg.slope, r: innerFDMreg.r, n: innerFDMreg.n },
  outerFDM: { slope: outerFDMreg.slope, r: outerFDMreg.r, n: outerFDMreg.n },
  innerStronger: {
    vdm: Math.abs(innerVDMreg.r) > Math.abs(outerVDMreg.r),
    fdm: Math.abs(innerFDMreg.r) > Math.abs(outerFDMreg.r),
  },
};

console.log(`╔══════════════════════════════════════════════════════════════╗`);
console.log(`║  PHASE 3: RESIDUAL TESTS (remove Vmax dependence)          ║`);
console.log(`╚══════════════════════════════════════════════════════════════╝\n`);

function computeResiduals(gals, yKey, yTransform) {
  const xs = gals.map(g => g.logVmax);
  const ys = gals.map(g => yTransform ? yTransform(g[yKey]) : g[yKey]);
  const baseline = linearRegression(xs, ys);
  const residuals = gals.map((g, i) => {
    const predicted = baseline.intercept + baseline.slope * g.logVmax;
    return { ...g, residual: ys[i] - predicted };
  });
  const resReg = linearRegression(residuals.map(r => r.logSigBar), residuals.map(r => r.residual));
  return { baseline, residuals, residualReg: resReg };
}

const vdmResid = computeResiduals(vdmGals, 'logMeanVDM', null);
const rdomResid = computeResiduals(rdomGals, 'rDMdomNorm', null);
const alphaResid = computeResiduals(alphaGals, 'alpha', null);

console.log(`  ΔV_DM (residual after Vmax regression) vs log(Σ_bar):`);
console.log(`    slope = ${vdmResid.residualReg.slope.toFixed(5)}, r = ${vdmResid.residualReg.r.toFixed(4)}, n = ${vdmResid.residualReg.n}`);
console.log(`    Baseline R² (Vmax only): ${vdmResid.baseline.r2.toFixed(4)}`);
console.log(`    Σ_bar still explains residuals? ${Math.abs(vdmResid.residualReg.r) > 0.15 ? 'YES ✓' : 'Weak'}\n`);

console.log(`  Δr_DMdom (residual after Vmax regression) vs log(Σ_bar):`);
console.log(`    slope = ${rdomResid.residualReg.slope.toFixed(5)}, r = ${rdomResid.residualReg.r.toFixed(4)}, n = ${rdomResid.residualReg.n}`);
console.log(`    Baseline R² (Vmax only): ${rdomResid.baseline.r2.toFixed(4)}`);
console.log(`    Σ_bar still explains residuals? ${Math.abs(rdomResid.residualReg.r) > 0.15 ? 'YES ✓' : 'Weak'}\n`);

console.log(`  Δα (residual after Vmax regression) vs log(Σ_bar):`);
console.log(`    slope = ${alphaResid.residualReg.slope.toFixed(5)}, r = ${alphaResid.residualReg.r.toFixed(4)}, n = ${alphaResid.residualReg.n}`);
console.log(`    Baseline R² (Vmax only): ${alphaResid.baseline.r2.toFixed(4)}`);
console.log(`    Σ_bar still explains residuals? ${Math.abs(alphaResid.residualReg.r) > 0.15 ? 'YES ✓' : 'Weak'}\n`);

const residualTests = {
  vdm: {
    baselineR2: vdmResid.baseline.r2,
    residualSlope: vdmResid.residualReg.slope,
    residualIntercept: vdmResid.residualReg.intercept,
    residualR: vdmResid.residualReg.r,
    residualR2: vdmResid.residualReg.r2,
    n: vdmResid.residualReg.n,
    sigBarExplains: Math.abs(vdmResid.residualReg.r) > 0.15,
  },
  rDMdom: {
    baselineR2: rdomResid.baseline.r2,
    residualSlope: rdomResid.residualReg.slope,
    residualIntercept: rdomResid.residualReg.intercept,
    residualR: rdomResid.residualReg.r,
    residualR2: rdomResid.residualReg.r2,
    n: rdomResid.residualReg.n,
    sigBarExplains: Math.abs(rdomResid.residualReg.r) > 0.15,
  },
  alpha: {
    baselineR2: alphaResid.baseline.r2,
    residualSlope: alphaResid.residualReg.slope,
    residualIntercept: alphaResid.residualReg.intercept,
    residualR: alphaResid.residualReg.r,
    residualR2: alphaResid.residualReg.r2,
    n: alphaResid.residualReg.n,
    sigBarExplains: Math.abs(alphaResid.residualReg.r) > 0.15,
  },
};

console.log(`╔══════════════════════════════════════════════════════════════╗`);
console.log(`║  PHASE 4: COUPLING LAW FITS                               ║`);
console.log(`╚══════════════════════════════════════════════════════════════╝\n`);

const vdmLaw = linearRegression(vdmGals.map(g => g.logSigBar), vdmGals.map(g => g.logMeanVDM));
const rdomLaw = linearRegression(rdomGals.map(g => g.logSigBar), rdomGals.map(g => g.rDMdomNorm));
const alphaLaw = linearRegression(alphaGals.map(g => g.logSigBar), alphaGals.map(g => g.alpha));

console.log(`  ⟨V_DM⟩ = 10^(${vdmLaw.intercept.toFixed(3)} + ${vdmLaw.slope.toFixed(4)} · log Σ_bar)   R²=${vdmLaw.r2.toFixed(3)}`);
console.log(`  r_DMdom = ${rdomLaw.intercept.toFixed(3)} + ${rdomLaw.slope.toFixed(4)} · log Σ_bar   R²=${rdomLaw.r2.toFixed(3)}`);
console.log(`  α = ${alphaLaw.intercept.toFixed(3)} + ${alphaLaw.slope.toFixed(4)} · log Σ_bar   R²=${alphaLaw.r2.toFixed(3)}\n`);

const couplingLaws = {
  vdm: { intercept: vdmLaw.intercept, slope: vdmLaw.slope, r2: vdmLaw.r2, r: vdmLaw.r, formula: 'log⟨V_DM⟩ = A + B·logΣ_bar' },
  rDMdom: { intercept: rdomLaw.intercept, slope: rdomLaw.slope, r2: rdomLaw.r2, r: rdomLaw.r, formula: 'r_DMdom/r_max = C + D·logΣ_bar' },
  alpha: { intercept: alphaLaw.intercept, slope: alphaLaw.slope, r2: alphaLaw.r2, r: alphaLaw.r, formula: 'α = E + F·logΣ_bar' },
};

console.log(`╔══════════════════════════════════════════════════════════════╗`);
console.log(`║  PHASE 5: MASS-BINNED COUPLING                            ║`);
console.log(`╚══════════════════════════════════════════════════════════════╝\n`);

const massBins = [
  { label: 'Low mass (V<80)', min: 0, max: 80 },
  { label: 'High mass (V≥80)', min: 80, max: 9999 },
];

const massBinResults = [];
for (const bin of massBins) {
  const binGals = vdmGals.filter(g => g.vmax >= bin.min && g.vmax < bin.max);
  if (binGals.length < 5) { massBinResults.push({ label: bin.label, n: binGals.length, skipped: true }); continue; }

  const reg = linearRegression(binGals.map(g => g.logSigBar), binGals.map(g => g.logMeanVDM));
  console.log(`  ${bin.label}: V_DM slope=${reg.slope.toFixed(4)}, r=${reg.r.toFixed(4)}, n=${binGals.length}`);
  massBinResults.push({ label: bin.label, n: binGals.length, slope: reg.slope, r: reg.r, skipped: false });
}

console.log(`\n${'═'.repeat(70)}`);
console.log(`  COUPLING VERDICT`);
console.log(`${'═'.repeat(70)}\n`);

const nResidualsPassing = Object.values(residualTests).filter(r => r.sigBarExplains).length;

const tests = [
  { name: 'V_DM tracks Σ_bar', pass: fingerprint[0].matchesExpected && Math.abs(fingerprint[0].r) > 0.3 },
  { name: 'r_DMdom tracks Σ_bar', pass: fingerprint[1].matchesExpected && Math.abs(fingerprint[1].r) > 0.2 },
  { name: 'α tracks Σ_bar (sign + |r|>0.1)', pass: fingerprint.length > 2 && fingerprint[2].matchesExpected && Math.abs(fingerprint[2].r) > 0.1 },
  { name: 'Partial r|Vmax ≥2 metrics >0.2', pass: fingerprint.filter(f => Math.abs(f.partialR) > 0.2).length >= 2 },
  { name: 'Residuals ≥2/3 pass', pass: nResidualsPassing >= 2 },
  { name: 'Inner coupling stronger', pass: innerOuter.innerStronger.vdm || innerOuter.innerStronger.fdm },
];

let passCount = 0;
for (const t of tests) {
  console.log(`  ${t.pass ? '✓' : '✗'} ${t.name}`);
  if (t.pass) passCount++;
}

const couplingWins = passCount >= 4;
console.log(`\n  Score: ${passCount}/6`);
console.log(`  ${couplingWins ? '→ COUPLING WINS: Baryonic surface density tracks halo structure itself' : '→ Evidence is mixed'}`);

const perGalaxyData = vdmGals.map(g => ({
  name: g.name,
  logSigBar: +g.logSigBar.toFixed(3),
  logVDM: +g.logMeanVDM.toFixed(3),
  vmax: +g.vmax.toFixed(1),
}));

const rdomPlotData = rdomGals.map(g => ({
  name: g.name,
  logSigBar: +g.logSigBar.toFixed(3),
  rDMdomNorm: +g.rDMdomNorm.toFixed(4),
  vmax: +g.vmax.toFixed(1),
}));

const alphaPlotData = alphaGals.map(g => ({
  name: g.name,
  logSigBar: +g.logSigBar.toFixed(3),
  alpha: +g.alpha.toFixed(4),
  vmax: +g.vmax.toFixed(1),
}));

const residualPlotData = {
  vdm: vdmResid.residuals.map(g => ({ logSigBar: +g.logSigBar.toFixed(3), residual: +g.residual.toFixed(4) })),
  rDMdom: rdomResid.residuals.map(g => ({ logSigBar: +g.logSigBar.toFixed(3), residual: +g.residual.toFixed(4) })),
  alpha: alphaResid.residuals.map(g => ({ logSigBar: +g.logSigBar.toFixed(3), residual: +g.residual.toFixed(4) })),
};

const couplingAnalysis = {
  fingerprint,
  innerOuter,
  residualTests,
  couplingLaws,
  massBins: massBinResults,
  tests: tests.map(t => ({ name: t.name, pass: t.pass })),
  passCount,
  couplingWins,
  verdict: couplingWins
    ? 'Baryonic surface density does not merely correlate with the apparent dark matter fraction; it tracks halo structure itself, consistent with direct baryon–halo coupling.'
    : `Evidence is mixed (${passCount}/6 tests pass). Some halo properties track Σ_bar beyond mass dependence, but the signal is not uniformly strong.`,
  plotData: {
    vdm: perGalaxyData,
    rDMdom: rdomPlotData,
    alpha: alphaPlotData,
    residuals: residualPlotData,
  },
};

const dataPath = path.join(__dirname, '..', 'public', 'fdm-analysis.json');
const data = JSON.parse(fs.readFileSync(dataPath, 'utf8'));
data.couplingAnalysis = couplingAnalysis;
fs.writeFileSync(dataPath, JSON.stringify(data, null, 2));
console.log(`\nSaved couplingAnalysis to fdm-analysis.json`);
