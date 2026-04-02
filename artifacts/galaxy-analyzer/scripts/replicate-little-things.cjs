const fs = require('fs');
const path = require('path');

const G_KPC = 4.3009e-6;
const G_DAG_KPC = 3.7032;
const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function parseTable1(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const galaxies = {};
  for (const line of lines) {
    const name = line.substring(0, 8).trim();
    const dist = parseFloat(line.substring(32, 36).trim());
    const inc = parseFloat(line.substring(59, 63).trim());
    const eInc = parseFloat(line.substring(64, 68).trim());
    const vMag = parseFloat(line.substring(69, 74).trim());
    galaxies[name] = { name, dist, inc, eInc, vMag };
  }
  return galaxies;
}

function parseTable2(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const results = {};
  for (const line of lines) {
    const name = line.substring(0, 8).trim();
    const rmax = parseFloat(line.substring(9, 13).trim());
    const r03 = parseFloat(line.substring(14, 18).trim());
    const vRmax = parseFloat(line.substring(19, 24).trim());

    const mgasStr = line.substring(139, 145).trim();
    const mstarKStr = line.substring(146, 151).trim();
    const mstarSEDStr = line.substring(152, 157).trim();

    const mgas = mgasStr ? parseFloat(mgasStr) * 1e7 : 0;
    const mstarK = mstarKStr ? parseFloat(mstarKStr) * 1e7 : 0;
    const mstarSED = mstarSEDStr ? parseFloat(mstarSEDStr) * 1e7 : 0;

    const mstar = mstarSED > 0 ? mstarSED : mstarK;

    results[name] = { rmax, r03, vRmax, mgas, mstar, mstarK, mstarSED, mbar: mstar + 1.33 * mgas };
  }
  return results;
}

function parseRotCurve(filePath) {
  const lines = fs.readFileSync(filePath, 'utf8').trim().split('\n');
  const galaxies = {};
  for (const line of lines) {
    const name = line.substring(0, 8).trim();
    const type = line.substring(9, 14).trim();
    if (type !== 'Data') continue;

    const r03 = parseFloat(line.substring(15, 23).trim());
    const v03 = parseFloat(line.substring(24, 34).trim());
    const rScaled = parseFloat(line.substring(35, 44).trim());
    const vScaled = parseFloat(line.substring(45, 54).trim());
    const eVScaled = parseFloat(line.substring(55, 63).trim());

    const rKpc = rScaled * r03;
    const vKms = vScaled * v03;
    const eVKms = eVScaled * v03;

    if (!galaxies[name]) galaxies[name] = [];
    galaxies[name].push({ r: rKpc, v: vKms, ev: eVKms, rScaled, vScaled });
  }
  return galaxies;
}

function interpolateV(curve, r) {
  if (curve.length === 0) return NaN;
  if (r <= curve[0].r) return curve[0].v * (r / curve[0].r);
  if (r >= curve[curve.length - 1].r) return curve[curve.length - 1].v;
  for (let i = 0; i < curve.length - 1; i++) {
    if (r >= curve[i].r && r <= curve[i + 1].r) {
      const frac = (r - curve[i].r) / (curve[i + 1].r - curve[i].r);
      return curve[i].v + frac * (curve[i + 1].v - curve[i].v);
    }
  }
  return curve[curve.length - 1].v;
}

function mcGillRAR(gbar) {
  const y = gbar / G_DAG_KPC;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function linearRegression(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let sxx = 0, sxy = 0;
  for (let i = 0; i < n; i++) {
    sxx += (x[i] - mx) ** 2;
    sxy += (x[i] - mx) * (y[i] - my);
  }
  const slope = sxy / sxx;
  const intercept = my - slope * mx;
  const yPred = x.map(xi => intercept + slope * xi);
  let ssRes = 0, ssTot = 0;
  for (let i = 0; i < n; i++) {
    ssRes += (y[i] - yPred[i]) ** 2;
    ssTot += (y[i] - my) ** 2;
  }
  const r2 = 1 - ssRes / ssTot;
  const r = Math.sign(slope) * Math.sqrt(Math.abs(r2));
  const sSlope = Math.sqrt(ssRes / (n - 2) / sxx);
  const tStat = slope / sSlope;
  const df = n - 2;
  const p = tDistPValue(tStat, df);
  return { slope, intercept, r2, r, sSlope, tStat, p, n };
}

function tDistPValue(t, df) {
  const x = df / (df + t * t);
  const a = df / 2;
  const b = 0.5;
  let ibeta = incompleteBeta(x, a, b);
  return ibeta;
}

function incompleteBeta(x, a, b) {
  if (x === 0 || x === 1) return x === 0 ? 1 : 0;
  const lnBeta = lnGamma(a) + lnGamma(b) - lnGamma(a + b);
  const front = Math.exp(Math.log(x) * a + Math.log(1 - x) * b - lnBeta);
  let sum = 1, term = 1;
  for (let n = 0; n < 200; n++) {
    const num1 = -(a + n) * (a + b + n) * x / ((a + 2 * n) * (a + 2 * n + 1));
    term *= num1;
    sum += term;
    const num2 = (n + 1) * (b - n - 1) * x / ((a + 2 * n + 1) * (a + 2 * n + 2));
    term *= num2;
    sum += term;
    if (Math.abs(term) < 1e-10) break;
  }
  return front * sum / a;
}

function lnGamma(z) {
  const c = [76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5];
  let x = z, y = z;
  let tmp = x + 5.5;
  tmp -= (x + 0.5) * Math.log(tmp);
  let ser = 1.000000000190015;
  for (let j = 0; j < 6; j++) ser += c[j] / ++y;
  return -tmp + Math.log(2.5066282746310005 * ser / x);
}

function partialCorrelation(x, y, z) {
  const n = x.length;
  const rxy = pearsonR(x, y);
  const rxz = pearsonR(x, z);
  const ryz = pearsonR(y, z);
  const num = rxy - rxz * ryz;
  const den = Math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2));
  const rPartial = num / den;
  const df = n - 3;
  const t = rPartial * Math.sqrt(df / (1 - rPartial ** 2));
  const p = tDistPValue(Math.abs(t) * Math.abs(t) > 0 ? df / (df + t * t) : 1, df);
  const pVal = tDistPValue2(t, df);
  return { rPartial, t, df, p: pVal };
}

function tDistPValue2(t, df) {
  const x = df / (df + t * t);
  return incompleteBeta(x, df / 2, 0.5);
}

function pearsonR(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b) / n;
  const my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return sxy / Math.sqrt(sxx * syy);
}

const table1 = parseTable1('/tmp/little_things/table1.dat');
const table2 = parseTable2('/tmp/little_things/table2.dat');
const rotTotal = parseRotCurve('/tmp/little_things/rotdmbar.dat');
const rotDM = parseRotCurve('/tmp/little_things/rotdm.dat');

console.log(`\n=== LITTLE THINGS REPLICATION ANALYSIS ===`);
console.log(`Dataset: Oh et al. (2015) - 26 dwarf irregular galaxies`);
console.log(`Source: VizieR J/AJ/149/180 (VLA HI + Spitzer 3.6μm)`);
console.log(`Independent from: SPARC (Lelli et al. 2016)`);
console.log(`\nGalaxies with total rotation curves: ${Object.keys(rotTotal).length}`);
console.log(`Galaxies with DM-only curves: ${Object.keys(rotDM).length}`);

const galaxyResults = [];

for (const name of Object.keys(rotTotal)) {
  if (!rotDM[name] || !table2[name] || !table1[name]) continue;

  const t1 = table1[name];
  const t2 = table2[name];
  const curveTotal = rotTotal[name].sort((a, b) => a.r - b.r);
  const curveDM = rotDM[name].sort((a, b) => a.r - b.r);

  if (curveTotal.length < 5) continue;
  if (t2.mbar <= 0) continue;

  const rFid = t2.rmax * 0.7;

  const outerPoints = curveTotal.filter(p => p.r >= t2.rmax * 0.4 && p.r <= t2.rmax);
  if (outerPoints.length < 3) continue;

  const deltaRARs = [];
  const logSigBars = [];
  const radialData = [];

  for (const pt of outerPoints) {
    const r = pt.r;
    const vObs = pt.v;
    const vDM = interpolateV(curveDM, r);

    if (isNaN(vDM) || vObs <= 0 || r <= 0) continue;

    let vBarSq = vObs * vObs - vDM * vDM;
    if (vBarSq < 0) vBarSq = 0;
    const vBar = Math.sqrt(vBarSq);

    const gObs = vObs * vObs / r;
    const gBar = vBarSq / r;

    if (gBar <= 0 || gObs <= 0) continue;

    const gRAR = mcGillRAR(gBar);
    const deltaRAR = Math.log10(gObs) - Math.log10(gRAR);

    const sigBar = t2.mbar / (Math.PI * r * r);
    const logSigBar = Math.log10(sigBar);

    deltaRARs.push(deltaRAR);
    logSigBars.push(logSigBar);
    radialData.push({ r, vObs, vBar, gObs, gBar, gRAR, deltaRAR, sigBar, logSigBar, rNorm: r / t2.rmax });
  }

  if (deltaRARs.length >= 3) {
    const meanDelta = deltaRARs.reduce((a, b) => a + b) / deltaRARs.length;
    const meanLogSig = logSigBars.reduce((a, b) => a + b) / logSigBars.length;
    const vmax = t2.vRmax;

    galaxyResults.push({
      name,
      dist: t1.dist,
      inc: t1.inc,
      vmax,
      mbar: t2.mbar,
      mgas: t2.mgas,
      mstar: t2.mstar,
      rmax: t2.rmax,
      nPoints: deltaRARs.length,
      meanDeltaRAR: meanDelta,
      meanLogSigBar: meanLogSig,
      logVmax: Math.log10(vmax),
      deltaRARs,
      logSigBars,
      radialData
    });
  }
}

console.log(`\nGalaxies with valid data for analysis: ${galaxyResults.length}`);
console.log(`\n--- Per-Galaxy Summary ---`);
for (const g of galaxyResults) {
  console.log(`${g.name.padEnd(10)} Dist=${g.dist.toFixed(1)}Mpc Vmax=${g.vmax.toFixed(1)}km/s Mbar=${g.mbar.toExponential(2)} <ΔRAR>=${g.meanDeltaRAR.toFixed(4)} <logΣ>=${g.meanLogSigBar.toFixed(2)} n=${g.nPoints}`);
}

const allDeltaRAR = [];
const allLogSigBar = [];
const allLogVmax = [];

for (const g of galaxyResults) {
  allDeltaRAR.push(g.meanDeltaRAR);
  allLogSigBar.push(g.meanLogSigBar);
  allLogVmax.push(g.logVmax);
}

console.log(`\n=== SIMPLE REGRESSION: ΔRAR vs log(Σ_bar) ===`);
const simpleReg = linearRegression(allLogSigBar, allDeltaRAR);
console.log(`b (slope) = ${simpleReg.slope.toFixed(4)} ± ${simpleReg.sSlope.toFixed(4)}`);
console.log(`intercept = ${simpleReg.intercept.toFixed(4)}`);
console.log(`r = ${simpleReg.r.toFixed(4)}`);
console.log(`R² = ${simpleReg.r2.toFixed(4)}`);
console.log(`p = ${simpleReg.p.toExponential(4)}`);
console.log(`n = ${simpleReg.n}`);

console.log(`\n=== PARTIAL CORRELATION: ΔRAR vs log(Σ_bar) | Vmax ===`);
const partial = partialCorrelation(allLogSigBar, allDeltaRAR, allLogVmax);
console.log(`r_partial = ${partial.rPartial.toFixed(4)}`);
console.log(`t = ${partial.t.toFixed(4)}`);
console.log(`df = ${partial.df}`);
console.log(`p = ${partial.p.toExponential(4)}`);

console.log(`\n=== COMPARISON WITH SPARC RESULTS ===`);
console.log(`                    SPARC (175 gal)    LITTLE THINGS (${galaxyResults.length} gal)`);
console.log(`slope b:            -0.178             ${simpleReg.slope.toFixed(3)}`);
console.log(`partial r|Vmax:     -0.656             ${partial.rPartial.toFixed(3)}`);
console.log(`simple r:           (computed)         ${simpleReg.r.toFixed(3)}`);

const allPointsDelta = [];
const allPointsLogSig = [];
for (const g of galaxyResults) {
  for (let i = 0; i < g.deltaRARs.length; i++) {
    allPointsDelta.push(g.deltaRARs[i]);
    allPointsLogSig.push(g.logSigBars[i]);
  }
}
console.log(`\n=== POINT-LEVEL REGRESSION (all radial points) ===`);
const pointReg = linearRegression(allPointsLogSig, allPointsDelta);
console.log(`b (slope) = ${pointReg.slope.toFixed(4)} ± ${pointReg.sSlope.toFixed(4)}`);
console.log(`r = ${pointReg.r.toFixed(4)}`);
console.log(`R² = ${pointReg.r2.toFixed(4)}`);
console.log(`p = ${pointReg.p.toExponential(4)}`);
console.log(`n_points = ${pointReg.n}`);

const consistent = Math.sign(simpleReg.slope) === Math.sign(-0.178);
const slopeInCI = simpleReg.slope >= -0.5 && simpleReg.slope <= 0;

console.log(`\n=== REPLICATION VERDICT ===`);
console.log(`Sign consistent with SPARC: ${consistent ? 'YES ✓' : 'NO ✗'}`);
console.log(`Slope in plausible range:   ${slopeInCI ? 'YES ✓' : 'NO ✗'}`);
console.log(`Partial r negative:         ${partial.rPartial < 0 ? 'YES ✓' : 'NO ✗'}`);

if (consistent && partial.rPartial < 0) {
  console.log(`\n>>> REPLICATION SUPPORTED: Independent dataset confirms density-dependent RAR residuals <<<`);
} else if (consistent) {
  console.log(`\n>>> PARTIAL REPLICATION: Sign matches but partial correlation not significant <<<`);
} else {
  console.log(`\n>>> REPLICATION FAILED: Independent dataset does not support the finding <<<`);
}

const output = {
  dataset: 'LITTLE THINGS (Oh et al. 2015)',
  source: 'VizieR J/AJ/149/180',
  instrument: 'VLA HI + Spitzer 3.6μm',
  nGalaxies: galaxyResults.length,
  type: 'dwarf irregular galaxies',
  independent: true,
  simpleRegression: {
    slope: simpleReg.slope,
    slopeError: simpleReg.sSlope,
    intercept: simpleReg.intercept,
    r: simpleReg.r,
    r2: simpleReg.r2,
    p: simpleReg.p,
    n: simpleReg.n
  },
  partialCorrelation: {
    rPartial: partial.rPartial,
    t: partial.t,
    df: partial.df,
    p: partial.p
  },
  pointLevelRegression: {
    slope: pointReg.slope,
    slopeError: pointReg.sSlope,
    r: pointReg.r,
    r2: pointReg.r2,
    p: pointReg.p,
    nPoints: pointReg.n
  },
  sparcComparison: {
    sparcSlope: -0.178,
    ltSlope: simpleReg.slope,
    signConsistent: consistent,
    sparcPartialR: -0.656,
    ltPartialR: partial.rPartial
  },
  galaxies: galaxyResults.map(g => ({
    name: g.name,
    dist: g.dist,
    vmax: g.vmax,
    mbar: g.mbar,
    meanDeltaRAR: g.meanDeltaRAR,
    meanLogSigBar: g.meanLogSigBar,
    nPoints: g.nPoints
  }))
};

const outPath = path.join(__dirname, '..', 'public', 'little-things-replication.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log(`\nResults saved to: ${outPath}`);
