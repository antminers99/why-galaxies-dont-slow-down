const fs = require('fs');
const path = require('path');

const G_DAG_KPC = 3.7032;

function mcGillRAR(gbar) {
  const y = gbar / G_DAG_KPC;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

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
  for (let i = 0; i < n; i++) {
    const yp = intercept + slope * x[i];
    ssRes += (y[i] - yp) ** 2;
    ssTot += (y[i] - my) ** 2;
  }
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

console.log(`\n${'='.repeat(70)}`);
console.log(`  DARK MATTER FRACTION vs BARYONIC SURFACE DENSITY`);
console.log(`  f_DM(r) = (g_obs - g_bar) / g_obs  vs  Σ_bar(r)`);
console.log(`${'='.repeat(70)}`);

const sparcReal = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const sparcBase = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis.json'), 'utf8'));

const rotmodDir = '/tmp/rotmod';
const sparcTablePath = '/tmp/sparc_table.mrt';

const sparcTable = {};
const tableLines = fs.readFileSync(sparcTablePath, 'utf8').split('\n');
for (const line of tableLines) {
  const parts = line.trim().split(/\s+/);
  if (parts.length < 7) continue;
  const name = parts[0];
  if (name === 'Galaxy' || name.startsWith('#') || name.startsWith('-')) continue;
  const dist = parseFloat(parts[2]);
  const eD = parseFloat(parts[3]);
  const inc = parseFloat(parts[5]);
  const eInc = parseFloat(parts[6]);
  const Lum = parseFloat(parts[1]);
  if (!isNaN(dist) && !isNaN(inc)) sparcTable[name] = { dist, eD, inc, eInc, Lum };
}

const allPoints = [];
const galaxySummaries = [];

const rotmodFiles = fs.readdirSync(rotmodDir).filter(f => f.endsWith('_rotmod.dat'));
console.log(`\nProcessing ${rotmodFiles.length} SPARC galaxies...`);

for (const file of rotmodFiles) {
  const name = file.replace('_rotmod.dat', '');
  const info = sparcTable[name];
  if (!info) continue;

  const lines = fs.readFileSync(path.join(rotmodDir, file), 'utf8').trim().split('\n');
  const dataLines = lines.filter(l => !l.startsWith('#') && l.trim());

  const UPSILON_D = 0.5, UPSILON_B = 0.7;
  let totalMbar = 0;
  const points = [];

  for (const line of dataLines) {
    const p = line.trim().split(/\s+/).map(Number);
    if (p.length < 7) continue;
    const [r, vObs, eVobs, vGas, vDisk, vBulge] = [p[0], p[1], p[2], p[3], p[4], p[5]];
    if (r <= 0 || vObs <= 0) continue;

    const vBarSq = vGas * Math.abs(vGas) + UPSILON_D * vDisk * Math.abs(vDisk) + UPSILON_B * vBulge * Math.abs(vBulge);
    if (vBarSq <= 0) continue;

    const gObs = vObs * vObs / r;
    const gBar = vBarSq / r;
    if (gObs <= 0 || gBar <= 0) continue;

    const fDM = (gObs - gBar) / gObs;
    const sigBar_approx = gBar / (4.3009e-6 * Math.PI);

    points.push({ r, vObs, gObs, gBar, fDM, logSigBar: 0, rNorm: 0 });
  }

  if (points.length < 5) continue;

  const rmax = points[points.length - 1].r;
  const rFid = rmax * 0.5;

  const perGalaxy = sparcReal.perGalaxy?.find((g) => g.name === name);
  const mbar = perGalaxy?.sigma_bar ? perGalaxy.sigma_bar * Math.PI * rFid * rFid : null;
  const vmax = perGalaxy?.Vmax || points.reduce((mx, p) => Math.max(mx, p.vObs), 0);

  if (!perGalaxy || perGalaxy.sigma_bar <= 0) continue;

  for (const pt of points) {
    const sigBar = perGalaxy.sigma_bar * (rFid * rFid) / (pt.r * pt.r);
    pt.logSigBar = Math.log10(sigBar);
    pt.rNorm = pt.r / rmax;
  }

  const innerPts = points.filter(p => p.rNorm < 0.5 && p.fDM >= 0 && p.fDM <= 1);
  const outerPts = points.filter(p => p.rNorm >= 0.5 && p.fDM >= 0 && p.fDM <= 1);
  const validPts = points.filter(p => p.fDM >= 0 && p.fDM <= 1);

  if (validPts.length < 3) continue;

  const meanFDM = validPts.reduce((s, p) => s + p.fDM, 0) / validPts.length;
  const meanLogSig = validPts.reduce((s, p) => s + p.logSigBar, 0) / validPts.length;

  galaxySummaries.push({
    name,
    vmax,
    meanFDM,
    meanLogSigBar: meanLogSig,
    logVmax: Math.log10(vmax),
    nPoints: validPts.length,
    innerMeanFDM: innerPts.length > 0 ? innerPts.reduce((s, p) => s + p.fDM, 0) / innerPts.length : null,
    outerMeanFDM: outerPts.length > 0 ? outerPts.reduce((s, p) => s + p.fDM, 0) / outerPts.length : null,
    innerMeanLogSig: innerPts.length > 0 ? innerPts.reduce((s, p) => s + p.logSigBar, 0) / innerPts.length : null,
    outerMeanLogSig: outerPts.length > 0 ? outerPts.reduce((s, p) => s + p.logSigBar, 0) / outerPts.length : null,
    nInner: innerPts.length,
    nOuter: outerPts.length,
  });

  for (const pt of validPts) {
    allPoints.push({
      galaxy: name,
      vmax,
      logVmax: Math.log10(vmax),
      r: pt.r,
      rNorm: pt.rNorm,
      fDM: pt.fDM,
      logSigBar: pt.logSigBar,
      gObs: pt.gObs,
      gBar: pt.gBar,
      region: pt.rNorm < 0.5 ? 'inner' : 'outer',
    });
  }
}

console.log(`Valid galaxies: ${galaxySummaries.length}`);
console.log(`Total radial points (0 ≤ f_DM ≤ 1): ${allPoints.length}`);

console.log(`\n╔══════════════════════════════════════════════════╗`);
console.log(`║  1. ALL-POINT REGRESSION: f_DM vs log(Σ_bar)    ║`);
console.log(`╚══════════════════════════════════════════════════╝`);

const allX = allPoints.map(p => p.logSigBar);
const allY = allPoints.map(p => p.fDM);
const allReg = linearRegression(allX, allY);
console.log(`slope b = ${allReg.slope.toFixed(4)} ± ${allReg.sSlope.toFixed(4)}`);
console.log(`intercept = ${allReg.intercept.toFixed(4)}`);
console.log(`r = ${allReg.r.toFixed(4)}`);
console.log(`R² = ${allReg.r2.toFixed(4)}`);
console.log(`n = ${allReg.n}`);

console.log(`\n╔══════════════════════════════════════════════════╗`);
console.log(`║  2. PER-GALAXY REGRESSION: ⟨f_DM⟩ vs ⟨logΣ⟩    ║`);
console.log(`╚══════════════════════════════════════════════════╝`);

const galX = galaxySummaries.map(g => g.meanLogSigBar);
const galY = galaxySummaries.map(g => g.meanFDM);
const galZ = galaxySummaries.map(g => g.logVmax);
const galReg = linearRegression(galX, galY);
const galPartR = partialR(galX, galY, galZ);

console.log(`slope b = ${galReg.slope.toFixed(4)} ± ${galReg.sSlope.toFixed(4)}`);
console.log(`r = ${galReg.r.toFixed(4)}`);
console.log(`R² = ${galReg.r2.toFixed(4)}`);
console.log(`partial r|Vmax = ${galPartR.toFixed(4)}`);
console.log(`n = ${galReg.n}`);

console.log(`\n╔══════════════════════════════════════════════════╗`);
console.log(`║  3. BINNED BY Vmax                               ║`);
console.log(`╚══════════════════════════════════════════════════╝`);

const vmaxBins = [
  { name: 'Dwarf', min: 0, max: 80 },
  { name: 'Medium', min: 80, max: 150 },
  { name: 'Massive', min: 150, max: 400 },
];

const binnedResults = [];
for (const bin of vmaxBins) {
  const pts = allPoints.filter(p => p.vmax >= bin.min && p.vmax < bin.max);
  if (pts.length < 10) { binnedResults.push({ ...bin, slope: NaN, r: NaN, r2: NaN, n: pts.length }); continue; }
  const bx = pts.map(p => p.logSigBar);
  const by = pts.map(p => p.fDM);
  const reg = linearRegression(bx, by);
  console.log(`${bin.name.padEnd(10)} (${pts.length} pts): b=${reg.slope.toFixed(4)}, r=${reg.r.toFixed(4)}, R²=${reg.r2.toFixed(4)}`);
  binnedResults.push({ ...bin, slope: reg.slope, r: reg.r, r2: reg.r2, n: pts.length, intercept: reg.intercept });
}

console.log(`\n╔══════════════════════════════════════════════════╗`);
console.log(`║  4. INNER vs OUTER                               ║`);
console.log(`╚══════════════════════════════════════════════════╝`);

const innerPts = allPoints.filter(p => p.region === 'inner');
const outerPts = allPoints.filter(p => p.region === 'outer');
const innerReg = linearRegression(innerPts.map(p => p.logSigBar), innerPts.map(p => p.fDM));
const outerReg = linearRegression(outerPts.map(p => p.logSigBar), outerPts.map(p => p.fDM));

console.log(`Inner (r < 0.5 rmax): b=${innerReg.slope.toFixed(4)}, r=${innerReg.r.toFixed(4)}, R²=${innerReg.r2.toFixed(4)}, n=${innerReg.n}`);
console.log(`Outer (r ≥ 0.5 rmax): b=${outerReg.slope.toFixed(4)}, r=${outerReg.r.toFixed(4)}, R²=${outerReg.r2.toFixed(4)}, n=${outerReg.n}`);

console.log(`\n╔══════════════════════════════════════════════════╗`);
console.log(`║  5. LITTLE THINGS REPLICATION                    ║`);
console.log(`╚══════════════════════════════════════════════════╝`);

function parseTable1LT(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const g = {};
  for (const l of lines) { const n = l.substring(0,8).trim(); g[n] = { dist: parseFloat(l.substring(32,36)), inc: parseFloat(l.substring(59,63)) }; }
  return g;
}
function parseTable2LT(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const g = {};
  for (const l of lines) {
    const n = l.substring(0,8).trim();
    const mgas = parseFloat(l.substring(139,145)) * 1e7 || 0;
    const mstarK = parseFloat(l.substring(146,151)) * 1e7 || 0;
    const mstarSED = parseFloat(l.substring(152,157)) * 1e7 || 0;
    const mstar = mstarSED > 0 ? mstarSED : mstarK;
    g[n] = { rmax: parseFloat(l.substring(9,13)), vRmax: parseFloat(l.substring(19,24)), mgas, mstar, mbar: mstar + 1.33 * mgas };
  }
  return g;
}
function parseRotLT(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const g = {};
  for (const l of lines) {
    const n = l.substring(0,8).trim();
    if (l.substring(9,14).trim() !== 'Data') continue;
    const r03 = parseFloat(l.substring(15,23));
    const v03 = parseFloat(l.substring(24,34));
    const rS = parseFloat(l.substring(35,44));
    const vS = parseFloat(l.substring(45,54));
    if (!g[n]) g[n] = [];
    g[n].push({ r: rS * r03, v: vS * v03 });
  }
  return g;
}
function interpV(curve, r) {
  if (!curve || curve.length === 0) return NaN;
  if (r <= curve[0].r) return curve[0].v * (r / curve[0].r);
  if (r >= curve[curve.length-1].r) return curve[curve.length-1].v;
  for (let i = 0; i < curve.length-1; i++) {
    if (r >= curve[i].r && r <= curve[i+1].r) {
      const f = (r - curve[i].r) / (curve[i+1].r - curve[i].r);
      return curve[i].v + f * (curve[i+1].v - curve[i].v);
    }
  }
  return curve[curve.length-1].v;
}

let ltPoints = [];
let ltGalaxies = [];

try {
  const lt1 = parseTable1LT('/tmp/little_things/table1.dat');
  const lt2 = parseTable2LT('/tmp/little_things/table2.dat');
  const ltTotal = parseRotLT('/tmp/little_things/rotdmbar.dat');
  const ltDM = parseRotLT('/tmp/little_things/rotdm.dat');

  for (const name of Object.keys(ltTotal)) {
    if (!ltDM[name] || !lt2[name]) continue;
    const t2 = lt2[name];
    const cTotal = ltTotal[name].sort((a,b) => a.r - b.r);
    const cDM = ltDM[name].sort((a,b) => a.r - b.r);
    if (cTotal.length < 5 || t2.mbar <= 0) continue;

    const rmax = t2.rmax;
    const pts = [];

    for (const pt of cTotal) {
      const r = pt.r;
      const vObs = pt.v;
      const vDM = interpV(cDM, r);
      if (isNaN(vDM) || vObs <= 0 || r <= 0) continue;

      let vBarSq = vObs * vObs - vDM * vDM;
      if (vBarSq < 0) vBarSq = 0;

      const gObs = vObs * vObs / r;
      const gBar = vBarSq / r;
      if (gObs <= 0) continue;

      const fDM = gBar > 0 ? (gObs - gBar) / gObs : 1;
      if (fDM < 0 || fDM > 1) continue;

      const sigBar = t2.mbar / (Math.PI * r * r);
      const logSigBar = Math.log10(sigBar);
      const rNorm = r / rmax;

      pts.push({ r, fDM, logSigBar, rNorm, region: rNorm < 0.5 ? 'inner' : 'outer' });
      ltPoints.push({ galaxy: name, fDM, logSigBar, rNorm, region: rNorm < 0.5 ? 'inner' : 'outer', vmax: t2.vRmax });
    }

    if (pts.length >= 3) {
      const validPts = pts.filter(p => p.fDM >= 0 && p.fDM <= 1);
      ltGalaxies.push({
        name,
        vmax: t2.vRmax,
        meanFDM: validPts.reduce((s,p) => s + p.fDM, 0) / validPts.length,
        meanLogSigBar: validPts.reduce((s,p) => s + p.logSigBar, 0) / validPts.length,
        nPoints: validPts.length,
      });
    }
  }

  const ltReg = linearRegression(ltPoints.map(p => p.logSigBar), ltPoints.map(p => p.fDM));
  const ltGalReg = linearRegression(ltGalaxies.map(g => g.meanLogSigBar), ltGalaxies.map(g => g.meanFDM));

  console.log(`LT point-level: b=${ltReg.slope.toFixed(4)}, r=${ltReg.r.toFixed(4)}, R²=${ltReg.r2.toFixed(4)}, n=${ltReg.n}`);
  console.log(`LT per-galaxy:  b=${ltGalReg.slope.toFixed(4)}, r=${ltGalReg.r.toFixed(4)}, R²=${ltGalReg.r2.toFixed(4)}, n=${ltGalReg.n}`);
} catch (e) {
  console.log(`LITTLE THINGS data not available: ${e.message}`);
}

console.log(`\n╔══════════════════════════════════════════════════╗`);
console.log(`║  SUMMARY                                         ║`);
console.log(`╚══════════════════════════════════════════════════╝`);
console.log(`\nSPARC point-level:  b = ${allReg.slope.toFixed(4)}, r = ${allReg.r.toFixed(4)}, n = ${allReg.n}`);
console.log(`SPARC per-galaxy:   b = ${galReg.slope.toFixed(4)}, r = ${galReg.r.toFixed(4)}, partial_r|Vmax = ${galPartR.toFixed(4)}, n = ${galReg.n}`);
console.log(`b < 0 (SPARC): ${allReg.slope < 0 ? 'YES ✓' : 'NO ✗'}`);
console.log(`b < 0 (LT):    ${ltPoints.length > 0 ? (linearRegression(ltPoints.map(p=>p.logSigBar), ltPoints.map(p=>p.fDM)).slope < 0 ? 'YES ✓' : 'NO ✗') : 'N/A'}`);

const sparcSample = [];
for (let i = 0; i < allPoints.length; i += Math.max(1, Math.floor(allPoints.length / 500))) {
  sparcSample.push({ logSigBar: +allPoints[i].logSigBar.toFixed(3), fDM: +allPoints[i].fDM.toFixed(4), vmax: +allPoints[i].vmax.toFixed(1), region: allPoints[i].region });
}

const ltSample = [];
for (let i = 0; i < ltPoints.length; i += Math.max(1, Math.floor(ltPoints.length / 300))) {
  ltSample.push({ logSigBar: +ltPoints[i].logSigBar.toFixed(3), fDM: +ltPoints[i].fDM.toFixed(4), vmax: +ltPoints[i].vmax.toFixed(1), region: ltPoints[i].region });
}

const output = {
  sparc: {
    pointLevel: { slope: allReg.slope, slopeErr: allReg.sSlope, intercept: allReg.intercept, r: allReg.r, r2: allReg.r2, n: allReg.n },
    perGalaxy: { slope: galReg.slope, slopeErr: galReg.sSlope, intercept: galReg.intercept, r: galReg.r, r2: galReg.r2, partialR: galPartR, n: galReg.n },
    byVmax: binnedResults.map(b => ({ name: b.name, slope: b.slope, r: b.r, r2: b.r2, n: b.n, intercept: b.intercept })),
    innerOuter: {
      inner: { slope: innerReg.slope, r: innerReg.r, r2: innerReg.r2, n: innerReg.n, intercept: innerReg.intercept },
      outer: { slope: outerReg.slope, r: outerReg.r, r2: outerReg.r2, n: outerReg.n, intercept: outerReg.intercept },
    },
    galaxies: galaxySummaries.map(g => ({ name: g.name, vmax: g.vmax, meanFDM: +g.meanFDM.toFixed(4), meanLogSigBar: +g.meanLogSigBar.toFixed(3) })),
    samplePoints: sparcSample,
  },
  littleThings: {
    pointLevel: ltPoints.length > 0 ? (() => { const r = linearRegression(ltPoints.map(p=>p.logSigBar), ltPoints.map(p=>p.fDM)); return { slope: r.slope, slopeErr: r.sSlope, intercept: r.intercept, r: r.r, r2: r.r2, n: r.n }; })() : null,
    perGalaxy: ltGalaxies.length > 0 ? (() => { const r = linearRegression(ltGalaxies.map(g=>g.meanLogSigBar), ltGalaxies.map(g=>g.meanFDM)); return { slope: r.slope, r: r.r, r2: r.r2, n: r.n }; })() : null,
    galaxies: ltGalaxies.map(g => ({ name: g.name, vmax: g.vmax, meanFDM: +g.meanFDM.toFixed(4), meanLogSigBar: +g.meanLogSigBar.toFixed(3) })),
    samplePoints: ltSample,
  },
};

const outPath = path.join(__dirname, '..', 'public', 'fdm-analysis.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log(`\nSaved to: ${outPath}`);
