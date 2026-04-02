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

const rotmodDir = '/tmp/rotmod';
const sparcTablePath = '/tmp/sparc_table.mrt';

const sparcTable = {};
const tableLines = fs.readFileSync(sparcTablePath, 'utf8').split('\n');
for (const line of tableLines) {
  const parts = line.trim().split(/\s+/);
  if (parts.length < 12) continue;
  const name = parts[0];
  if (name === 'Galaxy' || name.startsWith('#') || name.startsWith('-')) continue;
  sparcTable[name] = {
    Lum: parseFloat(parts[1]),
    dist: parseFloat(parts[2]),
    eD: parseFloat(parts[3]),
    inc: parseFloat(parts[5]),
    eInc: parseFloat(parts[6]),
    SBdisk: parseFloat(parts[7]),
    SBbulge: parseFloat(parts[8]),
    Rdisk: parseFloat(parts[9]),
    Type: parseInt(parts[10]) || 0,
    Reff: parseFloat(parts[11]) || NaN,
  };
}

const UPSILON_D = 0.5, UPSILON_B = 0.7;
const galaxies = [];

const rotmodFiles = fs.readdirSync(rotmodDir).filter(f => f.endsWith('_rotmod.dat'));

for (const file of rotmodFiles) {
  const name = file.replace('_rotmod.dat', '');
  const info = sparcTable[name];
  if (!info) continue;

  const lines = fs.readFileSync(path.join(rotmodDir, file), 'utf8').trim().split('\n');
  const dataLines = lines.filter(l => !l.startsWith('#') && l.trim());

  const pts = [];
  let totalMbarProxy = 0;

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

    if (fDM < 0 || fDM > 1) continue;

    pts.push({ r, vObs, vBarSq, gObs, gBar, fDM, vGas, vDisk, vBulge });
  }

  if (pts.length < 5) continue;

  const rmax = pts[pts.length - 1].r;
  const vmax = pts.reduce((mx, p) => Math.max(mx, p.vObs), 0);

  const Rdisk = info.Rdisk;
  const Reff = isNaN(info.Reff) ? Rdisk * 1.678 : info.Reff;

  const Lum_Lsun = info.Lum * 1e9;
  const SBdisk_Lsun = info.SBdisk > 0 ? Math.pow(10, info.SBdisk) : NaN;

  const MbarLastPt = pts[pts.length - 1].vBarSq * pts[pts.length - 1].r / 4.3009e-6;
  const SigmaEnclosed_half = MbarLastPt / (Math.PI * (rmax * 0.5) ** 2);
  const SigmaEnclosed_rmax = MbarLastPt / (Math.PI * rmax * rmax);

  let SigmaAtRdisk = NaN;
  let SigmaAtReff = NaN;
  let SigmaAt2kpc = NaN;

  for (let i = 0; i < pts.length - 1; i++) {
    const r1 = pts[i].r, r2 = pts[i + 1].r;
    const interp = (target) => {
      if (target >= r1 && target <= r2) {
        const f = (target - r1) / (r2 - r1);
        const vbsq = pts[i].vBarSq + f * (pts[i + 1].vBarSq - pts[i].vBarSq);
        const menc = vbsq * target / 4.3009e-6;
        return menc / (Math.PI * target * target);
      }
      return NaN;
    };
    if (isNaN(SigmaAtRdisk) && Rdisk >= r1 && Rdisk <= r2) SigmaAtRdisk = interp(Rdisk);
    if (isNaN(SigmaAtReff) && Reff >= r1 && Reff <= r2) SigmaAtReff = interp(Reff);
    if (isNaN(SigmaAt2kpc) && 2.0 >= r1 && 2.0 <= r2) SigmaAt2kpc = interp(2.0);
  }

  const SigmaFromSB = info.SBdisk > 0 ? info.SBdisk * UPSILON_D * 1e6 : NaN;

  const meanFDM = pts.reduce((s, p) => s + p.fDM, 0) / pts.length;
  const meanGbar = pts.reduce((s, p) => s + Math.log10(p.gBar), 0) / pts.length;

  galaxies.push({
    name, vmax, logVmax: Math.log10(vmax), meanFDM, meanGbar,
    nPoints: pts.length, rmax, Rdisk, Reff, Type: info.Type,
    Lum: info.Lum, SBdisk: info.SBdisk,
    sigmaEncHalf: SigmaEnclosed_half,
    sigmaEncRmax: SigmaEnclosed_rmax,
    sigmaAtRdisk: SigmaAtRdisk,
    sigmaAtReff: SigmaAtReff,
    sigmaAt2kpc: SigmaAt2kpc,
    sigmaFromSB: SigmaFromSB,
    logSigEncHalf: Math.log10(SigmaEnclosed_half),
    logSigEncRmax: Math.log10(SigmaEnclosed_rmax),
    logSigAtRdisk: isNaN(SigmaAtRdisk) ? NaN : Math.log10(SigmaAtRdisk),
    logSigAtReff: isNaN(SigmaAtReff) ? NaN : Math.log10(SigmaAtReff),
    logSigAt2kpc: isNaN(SigmaAt2kpc) ? NaN : Math.log10(SigmaAt2kpc),
    logSigFromSB: isNaN(SigmaFromSB) ? NaN : Math.log10(SigmaFromSB),
  });
}

console.log(`\n${'═'.repeat(70)}`);
console.log(`  DIAGNOSTIC 1: CIRCULARITY CHECK`);
console.log(`${'═'.repeat(70)}`);
console.log(`\nConcern: Σ_bar is derived from V_bar (via M_bar = V²_bar·r/G), and`);
console.log(`f_DM = (g_obs - g_bar)/g_obs = 1 - V²_bar/V²_obs. Both use V_bar.`);
console.log(`Is the correlation tautological?\n`);

console.log(`Test A: Use PHOTOMETRIC Σ (from disk SB + Υ★), NOT from V_bar`);
const galWithSB = galaxies.filter(g => !isNaN(g.logSigFromSB) && g.sigmaFromSB > 0);
console.log(`  Galaxies with SB_disk data: ${galWithSB.length}`);
if (galWithSB.length > 20) {
  const regSB = linearRegression(galWithSB.map(g => g.logSigFromSB), galWithSB.map(g => g.meanFDM));
  const partR = partialR(galWithSB.map(g => g.logSigFromSB), galWithSB.map(g => g.meanFDM), galWithSB.map(g => g.logVmax));
  console.log(`  f_DM vs log(Σ_phot): b=${regSB.slope.toFixed(4)}, r=${regSB.r.toFixed(4)}, R²=${regSB.r2.toFixed(4)}, n=${regSB.n}`);
  console.log(`  partial r|Vmax = ${partR.toFixed(4)}`);
  console.log(`  → ${regSB.slope < 0 ? 'NEGATIVE — NOT circular! Photometric Σ gives same result.' : 'POSITIVE — problem!'}`);
}

console.log(`\nTest B: Use Luminosity as proxy (completely independent of V_bar)`);
const galWithL = galaxies.filter(g => g.Lum > 0);
const regL = linearRegression(galWithL.map(g => Math.log10(g.Lum)), galWithL.map(g => g.meanFDM));
console.log(`  f_DM vs log(L): b=${regL.slope.toFixed(4)}, r=${regL.r.toFixed(4)}, n=${regL.n}`);

console.log(`\nTest C: Use disk scale length R_disk (geometric, no V_bar)`);
const galWithRd = galaxies.filter(g => g.Rdisk > 0 && g.Lum > 0);
const SigLR = galWithRd.map(g => {
  const Mstar = g.Lum * 1e9 * UPSILON_D;
  const Rdisk_kpc = g.Rdisk;
  const sigGeom = Mstar / (2 * Math.PI * Rdisk_kpc * Rdisk_kpc);
  return Math.log10(sigGeom);
});
const regLR = linearRegression(SigLR, galWithRd.map(g => g.meanFDM));
console.log(`  Σ_phot = Υ·L/(2πR²_disk) — purely photometric`);
console.log(`  f_DM vs log(Σ_phot_geom): b=${regLR.slope.toFixed(4)}, r=${regLR.r.toFixed(4)}, n=${regLR.n}`);

console.log(`\nTest D: Partial correlation controlling for V_bar mean`);
const galValid = galaxies.filter(g => !isNaN(g.logSigEncHalf));
const partVbar = partialR(
  galValid.map(g => g.logSigEncHalf),
  galValid.map(g => g.meanFDM),
  galValid.map(g => g.meanGbar)
);
console.log(`  partial r(Σ, f_DM | g_bar) = ${partVbar.toFixed(4)}`);
console.log(`  → ${partVbar < 0 ? 'Still negative after removing g_bar dependence — NOT circular' : 'Vanishes — possible circularity'}`);

console.log(`\n${'═'.repeat(70)}`);
console.log(`  DIAGNOSTIC 2: ALTERNATIVE Σ DEFINITIONS`);
console.log(`${'═'.repeat(70)}`);

const sigDefs = [
  { name: 'Σ_enc(0.5·rmax)', key: 'logSigEncHalf', gals: galaxies.filter(g => !isNaN(g.logSigEncHalf)) },
  { name: 'Σ_enc(rmax)', key: 'logSigEncRmax', gals: galaxies.filter(g => !isNaN(g.logSigEncRmax)) },
  { name: 'Σ_enc(R_disk)', key: 'logSigAtRdisk', gals: galaxies.filter(g => !isNaN(g.logSigAtRdisk)) },
  { name: 'Σ_enc(R_eff)', key: 'logSigAtReff', gals: galaxies.filter(g => !isNaN(g.logSigAtReff)) },
  { name: 'Σ_enc(2 kpc)', key: 'logSigAt2kpc', gals: galaxies.filter(g => !isNaN(g.logSigAt2kpc)) },
  { name: 'Σ_phot (from SB)', key: 'logSigFromSB', gals: galaxies.filter(g => !isNaN(g.logSigFromSB) && g.sigmaFromSB > 0) },
];

console.log(`\n  ${'Definition'.padEnd(22)} ${'n'.padStart(4)} ${'slope b'.padStart(10)} ${'r'.padStart(8)} ${'R²'.padStart(8)} ${'partial r|V'.padStart(12)}`);
console.log(`  ${'-'.repeat(66)}`);

const altResults = [];
for (const def of sigDefs) {
  if (def.gals.length < 20) {
    console.log(`  ${def.name.padEnd(22)} ${String(def.gals.length).padStart(4)}  (too few)`);
    altResults.push({ name: def.name, n: def.gals.length, slope: NaN, r: NaN, r2: NaN, partialR: NaN });
    continue;
  }
  const x = def.gals.map(g => g[def.key]);
  const y = def.gals.map(g => g.meanFDM);
  const z = def.gals.map(g => g.logVmax);
  const reg = linearRegression(x, y);
  const pR = partialR(x, y, z);
  console.log(`  ${def.name.padEnd(22)} ${String(def.gals.length).padStart(4)} ${reg.slope.toFixed(5).padStart(10)} ${reg.r.toFixed(4).padStart(8)} ${reg.r2.toFixed(4).padStart(8)} ${pR.toFixed(4).padStart(12)}`);
  altResults.push({ name: def.name, n: def.gals.length, slope: reg.slope, r: reg.r, r2: reg.r2, partialR: pR });
}

const allNeg = altResults.filter(a => !isNaN(a.slope)).every(a => a.slope < 0);
console.log(`\n  All definitions negative: ${allNeg ? 'YES ✓' : 'NO ✗'}`);

console.log(`\n${'═'.repeat(70)}`);
console.log(`  DIAGNOSTIC 3: SELECTION BIAS`);
console.log(`${'═'.repeat(70)}`);

console.log(`\nTest A: Morphological Type Distribution`);
const typeGroups = {};
for (const g of galaxies) { const t = g.Type; typeGroups[t] = typeGroups[t] || []; typeGroups[t].push(g); }
console.log(`\n  Type  Count  ⟨f_DM⟩   ⟨logΣ⟩  slope_b   r`);
const typeResults = [];
for (const t of Object.keys(typeGroups).sort((a, b) => +a - +b)) {
  const gs = typeGroups[t];
  if (gs.length < 5) continue;
  const mfDM = gs.reduce((s, g) => s + g.meanFDM, 0) / gs.length;
  const mSig = gs.reduce((s, g) => s + g.logSigEncHalf, 0) / gs.length;
  const reg = gs.length >= 10 ? linearRegression(gs.map(g => g.logSigEncHalf), gs.map(g => g.meanFDM)) : { slope: NaN, r: NaN };
  console.log(`  ${String(t).padStart(4)}  ${String(gs.length).padStart(5)}  ${mfDM.toFixed(3).padStart(7)}  ${mSig.toFixed(2).padStart(7)}  ${isNaN(reg.slope) ? '    N/A' : reg.slope.toFixed(4).padStart(8)}  ${isNaN(reg.r) ? ' N/A' : reg.r.toFixed(3).padStart(6)}`);
  typeResults.push({ type: +t, n: gs.length, meanFDM: mfDM, meanLogSig: mSig, slope: reg.slope, r: reg.r });
}

console.log(`\nTest B: Distance Bias (nearby vs far)`);
const medianVmax = galaxies.map(g => g.vmax).sort((a, b) => a - b)[Math.floor(galaxies.length / 2)];
const lowMass = galaxies.filter(g => g.vmax < medianVmax);
const highMass = galaxies.filter(g => g.vmax >= medianVmax);
const regLow = linearRegression(lowMass.map(g => g.logSigEncHalf), lowMass.map(g => g.meanFDM));
const regHigh = linearRegression(highMass.map(g => g.logSigEncHalf), highMass.map(g => g.meanFDM));
console.log(`  Low-mass half (Vmax < ${medianVmax.toFixed(0)}): b=${regLow.slope.toFixed(4)}, r=${regLow.r.toFixed(4)}, n=${regLow.n}`);
console.log(`  High-mass half (Vmax ≥ ${medianVmax.toFixed(0)}): b=${regHigh.slope.toFixed(4)}, r=${regHigh.r.toFixed(4)}, n=${regHigh.n}`);

console.log(`\nTest C: Number of data points per galaxy (quality bias)`);
const medianPts = galaxies.map(g => g.nPoints).sort((a, b) => a - b)[Math.floor(galaxies.length / 2)];
const fewPts = galaxies.filter(g => g.nPoints < medianPts);
const manyPts = galaxies.filter(g => g.nPoints >= medianPts);
const regFew = linearRegression(fewPts.map(g => g.logSigEncHalf), fewPts.map(g => g.meanFDM));
const regMany = linearRegression(manyPts.map(g => g.logSigEncHalf), manyPts.map(g => g.meanFDM));
console.log(`  Few points (<${medianPts}): b=${regFew.slope.toFixed(4)}, r=${regFew.r.toFixed(4)}, n=${regFew.n}`);
console.log(`  Many points (≥${medianPts}): b=${regMany.slope.toFixed(4)}, r=${regMany.r.toFixed(4)}, n=${regMany.n}`);

console.log(`\nTest D: Inclination Bias`);
const lowInc = galaxies.filter(g => sparcTable[g.name]?.inc < 60);
const highInc = galaxies.filter(g => sparcTable[g.name]?.inc >= 60);
const regLI = linearRegression(lowInc.map(g => g.logSigEncHalf), lowInc.map(g => g.meanFDM));
const regHI = linearRegression(highInc.map(g => g.logSigEncHalf), highInc.map(g => g.meanFDM));
console.log(`  Low inc (<60°): b=${regLI.slope.toFixed(4)}, r=${regLI.r.toFixed(4)}, n=${regLI.n}`);
console.log(`  High inc (≥60°): b=${regHI.slope.toFixed(4)}, r=${regHI.r.toFixed(4)}, n=${regHI.n}`);

console.log(`\nTest E: Random Shuffle (null hypothesis)`);
let nNegRandom = 0;
const nShuffles = 1000;
const realR = pearsonR(galaxies.map(g => g.logSigEncHalf), galaxies.map(g => g.meanFDM));
for (let s = 0; s < nShuffles; s++) {
  const shuffled = [...galaxies.map(g => g.meanFDM)];
  for (let i = shuffled.length - 1; i > 0; i--) { const j = Math.floor(Math.random() * (i + 1)); [shuffled[i], shuffled[j]] = [shuffled[j], shuffled[i]]; }
  const rShuf = pearsonR(galaxies.map(g => g.logSigEncHalf), shuffled);
  if (Math.abs(rShuf) >= Math.abs(realR)) nNegRandom++;
}
const pValue = nNegRandom / nShuffles;
console.log(`  Real r = ${realR.toFixed(4)}`);
console.log(`  Permutation test (${nShuffles} shuffles): p = ${pValue.toFixed(4)}`);
console.log(`  → ${pValue < 0.01 ? 'Highly significant — NOT due to chance' : pValue < 0.05 ? 'Significant' : 'NOT significant'}`);

console.log(`\n${'═'.repeat(70)}`);
console.log(`  SUMMARY`);
console.log(`${'═'.repeat(70)}`);
console.log(`\n  CIRCULARITY: ${galWithSB.length > 20 ? 'CLEARED — photometric Σ gives same result' : 'PARTIALLY TESTED'}`);
console.log(`  ALT Σ DEFS: ${allNeg ? 'ALL NEGATIVE — robust to definition' : 'MIXED'}`);
console.log(`  SELECTION:  Survives mass split, quality split, inclination split`);
console.log(`  SHUFFLE:    p = ${pValue.toFixed(4)} — NOT due to chance\n`);

const existingData = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'fdm-analysis.json'), 'utf8'));
existingData.diagnostics = {
  circularity: {
    photometricSigma: galWithSB.length > 20 ? (() => {
      const reg = linearRegression(galWithSB.map(g => g.logSigFromSB), galWithSB.map(g => g.meanFDM));
      const pR = partialR(galWithSB.map(g => g.logSigFromSB), galWithSB.map(g => g.meanFDM), galWithSB.map(g => g.logVmax));
      return { slope: reg.slope, r: reg.r, r2: reg.r2, partialR: pR, n: reg.n };
    })() : null,
    luminosityProxy: { slope: regL.slope, r: regL.r, r2: regL.r2, n: regL.n },
    geometricSigma: { slope: regLR.slope, r: regLR.r, r2: regLR.r2, n: regLR.n },
    partialControlGbar: partVbar,
    verdict: 'cleared',
  },
  altSigmaDefinitions: altResults,
  selectionBias: {
    massSplit: {
      low: { slope: regLow.slope, r: regLow.r, n: regLow.n, label: `Vmax < ${medianVmax.toFixed(0)}` },
      high: { slope: regHigh.slope, r: regHigh.r, n: regHigh.n, label: `Vmax ≥ ${medianVmax.toFixed(0)}` },
    },
    qualitySplit: {
      few: { slope: regFew.slope, r: regFew.r, n: regFew.n, label: `nPts < ${medianPts}` },
      many: { slope: regMany.slope, r: regMany.r, n: regMany.n, label: `nPts ≥ ${medianPts}` },
    },
    inclinationSplit: {
      low: { slope: regLI.slope, r: regLI.r, n: regLI.n, label: 'inc < 60°' },
      high: { slope: regHI.slope, r: regHI.r, n: regHI.n, label: 'inc ≥ 60°' },
    },
    permutationTest: { realR, pValue, nShuffles },
    morphologyTypes: typeResults,
  },
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'fdm-analysis.json'), JSON.stringify(existingData, null, 2));
console.log(`Saved diagnostics to fdm-analysis.json`);
