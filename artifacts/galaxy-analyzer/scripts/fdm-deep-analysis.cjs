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
  const aic = n * Math.log(ssRes / n) + 4;
  const bic = n * Math.log(ssRes / n) + 2 * Math.log(n);
  return { slope, intercept, r2, r, sSlope, n, aic, bic, ssRes };
}

function pearsonR(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b) / n;
  const my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxy / Math.sqrt(sxx * syy);
}

const rotmodDir = '/tmp/rotmod';
const sparcTablePath = '/tmp/sparc_table.mrt';

const sparcTable = {};
const tableLines = fs.readFileSync(sparcTablePath, 'utf8').split('\n');
for (const line of tableLines) {
  const parts = line.trim().split(/\s+/);
  if (parts.length < 7) continue;
  const name = parts[0];
  if (name === 'Galaxy' || name.startsWith('#') || name.startsWith('-')) continue;
  sparcTable[name] = { dist: parseFloat(parts[2]), eD: parseFloat(parts[3]), inc: parseFloat(parts[5]), eInc: parseFloat(parts[6]) };
}

const UPSILON_D = 0.5, UPSILON_B = 0.7;
const allPoints = [];
const galaxies = [];

const rotmodFiles = fs.readdirSync(rotmodDir).filter(f => f.endsWith('_rotmod.dat'));

for (const file of rotmodFiles) {
  const name = file.replace('_rotmod.dat', '');
  const info = sparcTable[name];
  if (!info) continue;

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

    const vBar = Math.sqrt(vBarSq);
    const gObs = vObs * vObs / r;
    const gBar = vBarSq / r;
    const fDM = (gObs - gBar) / gObs;
    const vDMsq = vObs * vObs - vBarSq;

    if (fDM < 0 || fDM > 1) continue;
    if (gObs <= 0 || gBar <= 0) continue;

    pts.push({ r, vObs, vBar, gObs, gBar, fDM, vDMsq, logGratio: Math.log10(gObs / gBar), logMratio: Math.log10(vObs * vObs / vBarSq) });
  }

  if (pts.length < 5) continue;

  const rmax = pts[pts.length - 1].r;
  const vmax = pts.reduce((mx, p) => Math.max(mx, p.vObs), 0);
  const logVmax = Math.log10(vmax);

  const sparcReal = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
  const perGalaxy = sparcReal.perGalaxy?.find((g) => g.name === name);
  if (!perGalaxy || perGalaxy.sigma_bar <= 0) continue;

  let rDMdom = NaN;
  for (const pt of pts) {
    if (pt.vDMsq > pt.vBar * pt.vBar && isNaN(rDMdom)) {
      rDMdom = pt.r;
    }
  }

  const sigBarGal = perGalaxy.sigma_bar;
  const logSigBarGal = Math.log10(sigBarGal);

  for (const pt of pts) {
    const sigBar = sigBarGal * (rmax * 0.5) ** 2 / (pt.r * pt.r);
    const logSigBar = Math.log10(sigBar);
    allPoints.push({
      galaxy: name, vmax, logVmax, r: pt.r, rNorm: pt.r / rmax,
      fDM: pt.fDM, logSigBar, gObs: pt.gObs, gBar: pt.gBar,
      logGratio: pt.logGratio, logMratio: pt.logMratio,
      vDMsq: pt.vDMsq, vDM: pt.vDMsq > 0 ? Math.sqrt(pt.vDMsq) : 0,
      region: pt.r / rmax < 0.5 ? 'inner' : 'outer',
    });
  }

  const validPts = pts;
  const meanFDM = validPts.reduce((s, p) => s + p.fDM, 0) / validPts.length;
  const meanLogGratio = validPts.reduce((s, p) => s + p.logGratio, 0) / validPts.length;
  const meanLogMratio = validPts.reduce((s, p) => s + p.logMratio, 0) / validPts.length;
  const meanVDMsq = validPts.reduce((s, p) => s + p.vDMsq, 0) / validPts.length;

  galaxies.push({
    name, vmax, logVmax, logSigBar: logSigBarGal,
    meanFDM, meanLogGratio, meanLogMratio, meanVDMsq,
    rDMdom: isNaN(rDMdom) ? null : rDMdom, rDMdomNorm: isNaN(rDMdom) ? null : rDMdom / rmax,
    nPoints: validPts.length,
  });
}

console.log(`\n${'═'.repeat(70)}`);
console.log(`  DEEP ANALYSIS: FUNCTIONAL FORMS + MASS SCALING + HALO CONNECTION`);
console.log(`${'═'.repeat(70)}`);
console.log(`Galaxies: ${galaxies.length}, Points: ${allPoints.length}\n`);

console.log(`╔══════════════════════════════════════════════════════════════╗`);
console.log(`║  TEST 1: FUNCTIONAL FORM COMPARISON                        ║`);
console.log(`╚══════════════════════════════════════════════════════════════╝\n`);

const forms = [
  { name: 'f_DM = a + b·log(Σ)', yKey: 'fDM', yLabel: 'f_DM', galYKey: 'meanFDM' },
  { name: 'log(g_obs/g_bar) = a + b·log(Σ)', yKey: 'logGratio', yLabel: 'log(g_obs/g_bar)', galYKey: 'meanLogGratio' },
  { name: 'log(M_dyn/M_bar) = a + b·log(Σ)', yKey: 'logMratio', yLabel: 'log(V²_obs/V²_bar)', galYKey: 'meanLogMratio' },
];

const formResults = [];
for (const form of forms) {
  const ptX = allPoints.map(p => p.logSigBar);
  const ptY = allPoints.map(p => p[form.yKey]);
  const ptReg = linearRegression(ptX, ptY);

  const galX = galaxies.map(g => g.logSigBar);
  const galY = galaxies.map(g => g[form.galYKey]);
  const galReg = linearRegression(galX, galY);

  console.log(`${form.name}`);
  console.log(`  Point:  b=${ptReg.slope.toFixed(5)} ± ${ptReg.sSlope.toFixed(5)}, r=${ptReg.r.toFixed(4)}, R²=${ptReg.r2.toFixed(4)}, AIC=${ptReg.aic.toFixed(1)}`);
  console.log(`  Galaxy: b=${galReg.slope.toFixed(5)} ± ${galReg.sSlope.toFixed(5)}, r=${galReg.r.toFixed(4)}, R²=${galReg.r2.toFixed(4)}, AIC=${galReg.aic.toFixed(1)}\n`);

  formResults.push({
    name: form.name,
    label: form.yLabel,
    point: { slope: ptReg.slope, slopeErr: ptReg.sSlope, intercept: ptReg.intercept, r: ptReg.r, r2: ptReg.r2, aic: ptReg.aic, bic: ptReg.bic, n: ptReg.n },
    galaxy: { slope: galReg.slope, slopeErr: galReg.sSlope, intercept: galReg.intercept, r: galReg.r, r2: galReg.r2, aic: galReg.aic, bic: galReg.bic, n: galReg.n },
  });
}

const bestForm = formResults.reduce((best, f) => Math.abs(f.point.r) > Math.abs(best.point.r) ? f : best);
console.log(`BEST FORM (by |r|): ${bestForm.name} with r = ${bestForm.point.r.toFixed(4)}\n`);

console.log(`╔══════════════════════════════════════════════════════════════╗`);
console.log(`║  TEST 2: SLOPE b(Vmax) — MASS DEPENDENCE                   ║`);
console.log(`╚══════════════════════════════════════════════════════════════╝\n`);

const fineBins = [
  { name: 'V < 50', min: 0, max: 50 },
  { name: '50–80', min: 50, max: 80 },
  { name: '80–120', min: 80, max: 120 },
  { name: '120–170', min: 120, max: 170 },
  { name: '170–250', min: 170, max: 250 },
  { name: 'V > 250', min: 250, max: 9999 },
];

const slopeVsVmax = [];
for (const bin of fineBins) {
  const pts = allPoints.filter(p => p.vmax >= bin.min && p.vmax < bin.max);
  if (pts.length < 20) { 
    slopeVsVmax.push({ ...bin, slope: NaN, r: NaN, r2: NaN, n: pts.length, midVmax: (bin.min + Math.min(bin.max, 400)) / 2 });
    continue;
  }
  const reg = linearRegression(pts.map(p => p.logSigBar), pts.map(p => p.fDM));
  console.log(`  ${bin.name.padEnd(10)} (${pts.length} pts): b=${reg.slope.toFixed(5)}, r=${reg.r.toFixed(4)}, R²=${reg.r2.toFixed(4)}`);
  slopeVsVmax.push({ ...bin, slope: reg.slope, r: reg.r, r2: reg.r2, n: pts.length, midVmax: (bin.min + Math.min(bin.max, 400)) / 2 });
}

const validBins = slopeVsVmax.filter(b => !isNaN(b.slope));
const bVsV = linearRegression(validBins.map(b => Math.log10(b.midVmax)), validBins.map(b => b.slope));
console.log(`\n  b vs log(Vmax): slope=${bVsV.slope.toFixed(5)}, r=${bVsV.r.toFixed(4)}, R²=${bVsV.r2.toFixed(4)}`);
console.log(`  → Slope b ${bVsV.slope < 0 ? 'becomes more negative (stronger)' : 'becomes less negative (weaker)'} for more massive galaxies\n`);

console.log(`╔══════════════════════════════════════════════════════════════╗`);
console.log(`║  TEST 3: HALO CONNECTION                                    ║`);
console.log(`╚══════════════════════════════════════════════════════════════╝\n`);

const ptsWithVDM = allPoints.filter(p => p.vDMsq > 0);
const logVDM = ptsWithVDM.map(p => Math.log10(Math.sqrt(p.vDMsq)));
const logSigVDM = ptsWithVDM.map(p => p.logSigBar);
const vdmReg = linearRegression(logSigVDM, logVDM);
console.log(`  log(V_DM) vs log(Σ_bar):`);
console.log(`    b=${vdmReg.slope.toFixed(5)}, r=${vdmReg.r.toFixed(4)}, R²=${vdmReg.r2.toFixed(4)}, n=${vdmReg.n}\n`);

const galVDMsq = galaxies.filter(g => g.meanVDMsq > 0);
const galLogVDM = galVDMsq.map(g => Math.log10(Math.sqrt(g.meanVDMsq)));
const galLogSigVDM = galVDMsq.map(g => g.logSigBar);
const galVdmReg = linearRegression(galLogSigVDM, galLogVDM);
console.log(`  Per-galaxy log(⟨V_DM⟩) vs log(Σ_bar):`);
console.log(`    b=${galVdmReg.slope.toFixed(5)}, r=${galVdmReg.r.toFixed(4)}, R²=${galVdmReg.r2.toFixed(4)}, n=${galVdmReg.n}\n`);

const galWithRdom = galaxies.filter(g => g.rDMdomNorm !== null);
console.log(`  Galaxies with DM dominance radius: ${galWithRdom.length}`);
if (galWithRdom.length >= 10) {
  const rdomReg = linearRegression(galWithRdom.map(g => g.logSigBar), galWithRdom.map(g => g.rDMdomNorm));
  console.log(`  r_DMdom/r_max vs log(Σ_bar): b=${rdomReg.slope.toFixed(5)}, r=${rdomReg.r.toFixed(4)}, R²=${rdomReg.r2.toFixed(4)}, n=${rdomReg.n}`);
  
  const logRdom = galWithRdom.map(g => Math.log10(g.rDMdomNorm));
  const rdomLogReg = linearRegression(galWithRdom.map(g => g.logSigBar), logRdom);
  console.log(`  log(r_DMdom/r_max) vs log(Σ_bar): b=${rdomLogReg.slope.toFixed(5)}, r=${rdomLogReg.r.toFixed(4)}`);
}

const coreCuspProxy = allPoints.filter(p => p.rNorm < 0.3 && p.vDMsq > 0);
if (coreCuspProxy.length > 50) {
  const innerSlope = [];
  const galNames = [...new Set(coreCuspProxy.map(p => p.galaxy))];
  for (const gn of galNames) {
    const gpts = coreCuspProxy.filter(p => p.galaxy === gn).sort((a, b) => a.r - b.r);
    if (gpts.length < 3) continue;
    const reg = linearRegression(gpts.map(p => Math.log10(p.r)), gpts.map(p => Math.log10(Math.sqrt(p.vDMsq))));
    const gal = galaxies.find(g => g.name === gn);
    if (gal && !isNaN(reg.slope)) {
      innerSlope.push({ name: gn, alpha: reg.slope, logSigBar: gal.logSigBar, vmax: gal.vmax });
    }
  }
  if (innerSlope.length >= 10) {
    const alphaVsSig = linearRegression(innerSlope.map(s => s.logSigBar), innerSlope.map(s => s.alpha));
    console.log(`\n  Inner DM slope α (V_DM ∝ r^α) vs log(Σ_bar):`);
    console.log(`    b=${alphaVsSig.slope.toFixed(5)}, r=${alphaVsSig.r.toFixed(4)}, n=${alphaVsSig.n}`);
    console.log(`    Interpretation: ${alphaVsSig.slope > 0 ? 'Higher Σ → steeper inner rise (more cuspy)' : 'Higher Σ → shallower inner rise (more cored)'}`);
  }
}

console.log(`\n${'═'.repeat(70)}`);
console.log(`  SUMMARY`);
console.log(`${'═'.repeat(70)}`);
console.log(`\n  BEST FUNCTIONAL FORM: ${bestForm.name}`);
console.log(`    r = ${bestForm.point.r.toFixed(4)}, R² = ${bestForm.point.r2.toFixed(4)}`);
console.log(`  MASS SCALING: b becomes ${bVsV.slope < 0 ? 'MORE' : 'LESS'} negative for massive galaxies`);
console.log(`    b(Vmax) correlation: r = ${bVsV.r.toFixed(4)}`);
console.log(`  HALO: V_DM correlates with Σ_bar: r = ${vdmReg.r.toFixed(4)}`);

const fdmExisting = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'fdm-analysis.json'), 'utf8'));

fdmExisting.deepAnalysis = {
  functionalForms: formResults.map(f => ({
    name: f.name,
    label: f.label,
    pointSlope: f.point.slope,
    pointR: f.point.r,
    pointR2: f.point.r2,
    pointAIC: f.point.aic,
    pointN: f.point.n,
    galaxySlope: f.galaxy.slope,
    galaxyR: f.galaxy.r,
    galaxyR2: f.galaxy.r2,
    galaxyAIC: f.galaxy.aic,
    galaxyN: f.galaxy.n,
  })),
  bestForm: bestForm.name,
  slopeVsVmax: slopeVsVmax.filter(b => !isNaN(b.slope)).map(b => ({
    bin: b.name, midVmax: b.midVmax, slope: b.slope, r: b.r, r2: b.r2, n: b.n,
  })),
  slopeMassScaling: { slope: bVsV.slope, r: bVsV.r, r2: bVsV.r2, direction: bVsV.slope < 0 ? 'stronger for massive' : 'weaker for massive' },
  halo: {
    vDMvsSigBar: { pointSlope: vdmReg.slope, pointR: vdmReg.r, pointR2: vdmReg.r2, pointN: vdmReg.n, galSlope: galVdmReg.slope, galR: galVdmReg.r, galR2: galVdmReg.r2, galN: galVdmReg.n },
    dmDominanceRadius: galWithRdom.length >= 10 ? (() => {
      const reg = linearRegression(galWithRdom.map(g => g.logSigBar), galWithRdom.map(g => g.rDMdomNorm));
      return { slope: reg.slope, r: reg.r, r2: reg.r2, n: reg.n };
    })() : null,
    dmDomGalaxies: galWithRdom.map(g => ({ name: g.name, logSigBar: +g.logSigBar.toFixed(3), rDMdomNorm: +g.rDMdomNorm.toFixed(4), vmax: +g.vmax.toFixed(1) })),
    innerSlopeData: (() => {
      const coreCusp = allPoints.filter(p => p.rNorm < 0.3 && p.vDMsq > 0);
      const galNames = [...new Set(coreCusp.map(p => p.galaxy))];
      const slopes = [];
      for (const gn of galNames) {
        const gpts = coreCusp.filter(p => p.galaxy === gn).sort((a, b) => a.r - b.r);
        if (gpts.length < 3) continue;
        const reg = linearRegression(gpts.map(p => Math.log10(p.r)), gpts.map(p => Math.log10(Math.sqrt(p.vDMsq))));
        const gal = galaxies.find(g => g.name === gn);
        if (gal && !isNaN(reg.slope)) slopes.push({ name: gn, alpha: +reg.slope.toFixed(4), logSigBar: +gal.logSigBar.toFixed(3), vmax: +gal.vmax.toFixed(1) });
      }
      if (slopes.length < 10) return null;
      const alphaReg = linearRegression(slopes.map(s => s.logSigBar), slopes.map(s => s.alpha));
      return { slope: alphaReg.slope, r: alphaReg.r, r2: alphaReg.r2, n: alphaReg.n, galaxies: slopes };
    })(),
  },
};

const outPath = path.join(__dirname, '..', 'public', 'fdm-analysis.json');
fs.writeFileSync(outPath, JSON.stringify(fdmExisting, null, 2));
console.log(`\nSaved to: ${outPath}`);
