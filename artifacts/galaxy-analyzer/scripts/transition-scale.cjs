const fs = require('fs');
const path = require('path');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const A0_KNOWN = 3.7032;

function rmsScatter(vals) {
  if (vals.length < 2) return Infinity;
  const m = vals.reduce((a, b) => a + b, 0) / vals.length;
  return Math.sqrt(vals.reduce((a, v) => a + (v - m) ** 2, 0) / vals.length);
}

function medAbsDev(vals) {
  if (vals.length < 2) return Infinity;
  const s = [...vals].sort((a, b) => a - b);
  const med = s[Math.floor(s.length / 2)];
  const devs = vals.map(v => Math.abs(v - med)).sort((a, b) => a - b);
  return devs[Math.floor(devs.length / 2)];
}

console.log("=== TRANSITION SCALE ANALYSIS ===");
console.log("Is there a hidden scale g0 where ALL galaxies transition?\n");

const rarJson = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const galaxyLookup = {};
for (const g of rarJson.perGalaxy) galaxyLookup[g.name] = g;

const allPoints = [];
const perGalaxy = {};
const sparcDir = '/tmp/rotmod';
const files = fs.readdirSync(sparcDir).filter(f => f.endsWith('.dat'));

for (const file of files) {
  const name = file.replace('_rotmod.dat', '');
  const gInfo = galaxyLookup[name];
  if (!gInfo) continue;
  const lines = fs.readFileSync(path.join(sparcDir, file), 'utf8').trim().split('\n');
  for (const line of lines) {
    if (line.trim().startsWith('#') || line.trim() === '') continue;
    const parts = line.trim().split(/\s+/).map(Number);
    if (parts.length < 6) continue;
    const [r, vobs, evobs, vgas, vdisk, vbul] = parts;
    if (!isFinite(r) || !isFinite(vobs) || r <= 0 || vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * vdisk * Math.abs(vdisk) + UPSILON_BULGE * vbul * Math.abs(vbul) + vgas * Math.abs(vgas);
    const gObs = vobs * vobs / r;
    const gBar = Math.abs(vBarSq) / r;
    if (gBar <= 0 || gObs <= 0 || !isFinite(gBar) || !isFinite(gObs)) continue;
    const ratio = gObs / gBar;
    const pt = { galaxy: name, gObs, gBar, ratio, logGbar: Math.log10(gBar), logRatio: Math.log10(ratio) };
    allPoints.push(pt);
    if (!perGalaxy[name]) perGalaxy[name] = [];
    perGalaxy[name].push(pt);
  }
}

function parseRotCurve(fp) {
  const lines = fs.readFileSync(fp, 'utf8').trim().split('\n');
  const g = {};
  for (const line of lines) {
    const nm = line.substring(0, 8).trim();
    const tp = line.substring(9, 14).trim();
    if (tp !== 'Data') continue;
    if (!g[nm]) g[nm] = [];
    g[nm].push({ r: parseFloat(line.substring(35, 44).trim()) * parseFloat(line.substring(15, 23).trim()), v: parseFloat(line.substring(45, 54).trim()) * parseFloat(line.substring(24, 34).trim()) });
  }
  return g;
}
function parseT2(fp) {
  const lines = fs.readFileSync(fp, 'utf8').trim().split('\n');
  const r = {};
  for (const line of lines) {
    const nm = line.substring(0, 8).trim();
    const mgS = line.substring(139, 145).trim();
    const msS = line.substring(152, 157).trim();
    const mkS = line.substring(146, 151).trim();
    const mg = mgS ? parseFloat(mgS) * 1e7 : 0;
    const ms = msS ? parseFloat(msS) * 1e7 : (mkS ? parseFloat(mkS) * 1e7 : 0);
    r[nm] = { mbar: ms + 1.33 * mg };
  }
  return r;
}
function interpV(c, r) {
  if (c.length === 0) return NaN;
  if (r <= c[0].r) return c[0].v * (r / c[0].r);
  if (r >= c[c.length - 1].r) return c[c.length - 1].v;
  for (let i = 0; i < c.length - 1; i++) {
    if (r >= c[i].r && r <= c[i + 1].r) return c[i].v + (r - c[i].r) / (c[i + 1].r - c[i].r) * (c[i + 1].v - c[i].v);
  }
  return c[c.length - 1].v;
}

const rotT = parseRotCurve('/tmp/little_things/rotdmbar.dat');
const rotD = parseRotCurve('/tmp/little_things/rotdm.dat');
const t2 = parseT2('/tmp/little_things/table2.dat');
let ltCount = 0;
for (const nm of Object.keys(rotT)) {
  if (!rotD[nm] || !t2[nm]) continue;
  const cT = rotT[nm].sort((a, b) => a.r - b.r);
  const cD = rotD[nm].sort((a, b) => a.r - b.r);
  if (cT.length < 5 || t2[nm].mbar <= 0) continue;
  ltCount++;
  for (const pt of cT) {
    const r = pt.r, vO = pt.v;
    const vDM = interpV(cD, r);
    if (isNaN(vDM) || vO <= 0 || r <= 0) continue;
    let vBsq = vO * vO - vDM * vDM;
    if (vBsq < 0) vBsq = 0;
    const gObs = vO * vO / r, gBar = vBsq / r;
    if (gBar <= 0 || gObs <= 0) continue;
    const ratio = gObs / gBar;
    const p = { galaxy: "LT_" + nm, gObs, gBar, ratio, logGbar: Math.log10(gBar), logRatio: Math.log10(ratio) };
    allPoints.push(p);
    if (!perGalaxy["LT_" + nm]) perGalaxy["LT_" + nm] = [];
    perGalaxy["LT_" + nm].push(p);
  }
}

const galaxyNames = Object.keys(perGalaxy).filter(n => perGalaxy[n].length >= 5);
console.log("Total: " + allPoints.length + " points, " + galaxyNames.length + " galaxies (" + (galaxyNames.length - ltCount) + " SPARC + " + ltCount + " LT)\n");

console.log("============================================");
console.log("STEP 1: g_obs/g_bar vs g_bar (the transition plot)");
console.log("============================================\n");

allPoints.sort((a, b) => a.logGbar - b.logGbar);
const nBins = 30;
const minLog = allPoints[0].logGbar;
const maxLog = allPoints[allPoints.length - 1].logGbar;
const binWidth = (maxLog - minLog) / nBins;

const bins = [];
for (let i = 0; i < nBins; i++) {
  const lo = minLog + i * binWidth;
  const hi = lo + binWidth;
  const inBin = allPoints.filter(p => p.logGbar >= lo && p.logGbar < hi);
  if (inBin.length < 5) continue;
  const ratios = inBin.map(p => p.logRatio);
  const mean = ratios.reduce((a, b) => a + b, 0) / ratios.length;
  const median = [...ratios].sort((a, b) => a - b)[Math.floor(ratios.length / 2)];
  const scatter = rmsScatter(ratios);
  bins.push({ logGbarMid: (lo + hi) / 2, gBarMid: Math.pow(10, (lo + hi) / 2), meanLogRatio: mean, medianLogRatio: median, scatter, n: inBin.length });
}

console.log("  " + "log(g_bar)".padEnd(12) + "g_bar".padEnd(14) + "log(ratio)".padEnd(12) + "ratio".padEnd(10) + "scatter".padEnd(10) + "n".padEnd(6));
for (const b of bins) {
  console.log("  " + b.logGbarMid.toFixed(2).padEnd(12) + b.gBarMid.toExponential(2).padEnd(14) + b.medianLogRatio.toFixed(3).padEnd(12) + Math.pow(10, b.medianLogRatio).toFixed(3).padEnd(10) + b.scatter.toFixed(4).padEnd(10) + String(b.n).padEnd(6));
}

console.log("\n  Key observations:");
const highAccBins = bins.filter(b => b.logGbarMid > 1);
const lowAccBins = bins.filter(b => b.logGbarMid < -1);
if (highAccBins.length > 0 && lowAccBins.length > 0) {
  const highRatio = Math.pow(10, highAccBins.reduce((a, b) => a + b.medianLogRatio, 0) / highAccBins.length);
  const lowRatio = Math.pow(10, lowAccBins.reduce((a, b) => a + b.medianLogRatio, 0) / lowAccBins.length);
  console.log("  High acceleration (g_bar >> a0): median ratio = " + highRatio.toFixed(3) + " (should be ~1)");
  console.log("  Low acceleration  (g_bar << a0): median ratio = " + lowRatio.toFixed(1) + " (DM dominates)");
}

console.log("\n============================================");
console.log("STEP 2: Finding the transition scale g0");
console.log("============================================\n");

console.log("  Method 1: Where does log(g_obs/g_bar) cross 0.1 dex (ratio = 1.26)?");
let transitionBin1 = null;
for (let i = bins.length - 1; i >= 0; i--) {
  if (bins[i].medianLogRatio > 0.1) {
    transitionBin1 = bins[i];
    if (i < bins.length - 1) {
      const next = bins[i + 1];
      const frac = (0.1 - next.medianLogRatio) / (bins[i].medianLogRatio - next.medianLogRatio);
      transitionBin1 = { logGbarMid: next.logGbarMid + frac * (bins[i].logGbarMid - next.logGbarMid) };
    }
    break;
  }
}
if (transitionBin1) {
  console.log("  Transition at log(g_bar) = " + transitionBin1.logGbarMid.toFixed(3) + " -> g0 = " + Math.pow(10, transitionBin1.logGbarMid).toFixed(3) + " (km/s)^2/kpc");
}

console.log("\n  Method 2: Maximum curvature (second derivative of ratio curve)");
let maxCurv = 0, maxCurvIdx = -1;
for (let i = 1; i < bins.length - 1; i++) {
  const d2 = bins[i + 1].medianLogRatio - 2 * bins[i].medianLogRatio + bins[i - 1].medianLogRatio;
  if (Math.abs(d2) > maxCurv) {
    maxCurv = Math.abs(d2);
    maxCurvIdx = i;
  }
}
if (maxCurvIdx >= 0) {
  console.log("  Max curvature at log(g_bar) = " + bins[maxCurvIdx].logGbarMid.toFixed(3) + " -> g0 = " + Math.pow(10, bins[maxCurvIdx].logGbarMid).toFixed(3) + " (km/s)^2/kpc");
}

console.log("\n  Method 3: Per-galaxy transition point extraction");
const perGalaxyG0 = [];
for (const name of galaxyNames) {
  const pts = perGalaxy[name].sort((a, b) => a.gBar - b.gBar);
  if (pts.length < 5) continue;
  for (let i = 0; i < pts.length - 1; i++) {
    if (pts[i].logRatio > 0.1 && pts[i + 1].logRatio <= 0.1) {
      const frac = (0.1 - pts[i + 1].logRatio) / (pts[i].logRatio - pts[i + 1].logRatio);
      const logG0 = pts[i + 1].logGbar + frac * (pts[i].logGbar - pts[i + 1].logGbar);
      perGalaxyG0.push({ name, g0: Math.pow(10, logG0), logG0 });
      break;
    }
  }
}
if (perGalaxyG0.length > 0) {
  const g0vals = perGalaxyG0.map(g => g.logG0);
  const meanG0 = g0vals.reduce((a, b) => a + b, 0) / g0vals.length;
  const medG0 = [...g0vals].sort((a, b) => a - b)[Math.floor(g0vals.length / 2)];
  const scatterG0 = rmsScatter(g0vals);
  const madG0 = medAbsDev(g0vals);
  console.log("  Galaxies with measurable transition: " + perGalaxyG0.length + "/" + galaxyNames.length);
  console.log("  Mean log(g0)   = " + meanG0.toFixed(3) + " -> g0 = " + Math.pow(10, meanG0).toFixed(3));
  console.log("  Median log(g0) = " + medG0.toFixed(3) + " -> g0 = " + Math.pow(10, medG0).toFixed(3));
  console.log("  Scatter (RMS)  = " + scatterG0.toFixed(3) + " dex");
  console.log("  MAD            = " + madG0.toFixed(3) + " dex");
  console.log("  Known a0 = " + A0_KNOWN.toFixed(3) + " (km/s)^2/kpc [log = " + Math.log10(A0_KNOWN).toFixed(3) + "]");
  console.log("  Consistency with a0: " + (Math.abs(medG0 - Math.log10(A0_KNOWN)) < 0.3 ? "CONSISTENT" : "DIFFERENT") + " (offset = " + (medG0 - Math.log10(A0_KNOWN)).toFixed(3) + " dex)");
}

console.log("\n  Method 4: Fitting the McGaugh function to extract a0 per galaxy");
function fitA0(points) {
  let bestA0 = 0, bestRMS = Infinity;
  for (let logA = -0.5; logA <= 1.5; logA += 0.005) {
    const a0 = Math.pow(10, logA);
    const resids = points.map(p => {
      const y = p.gBar / a0;
      const pred = p.gBar / (1 - Math.exp(-Math.sqrt(y)));
      return Math.log10(p.gObs) - Math.log10(pred);
    });
    const rms = rmsScatter(resids);
    if (rms < bestRMS) { bestRMS = rms; bestA0 = a0; }
  }
  return { a0: bestA0, logA0: Math.log10(bestA0), rms: bestRMS };
}

const perGalaxyA0 = [];
for (const name of galaxyNames) {
  const pts = perGalaxy[name];
  if (pts.length < 5) continue;
  const fit = fitA0(pts);
  if (fit.rms < 0.5) perGalaxyA0.push({ name, ...fit });
}

if (perGalaxyA0.length > 0) {
  const a0vals = perGalaxyA0.map(g => g.logA0);
  const meanA0 = a0vals.reduce((a, b) => a + b, 0) / a0vals.length;
  const medA0 = [...a0vals].sort((a, b) => a - b)[Math.floor(a0vals.length / 2)];
  const scatterA0 = rmsScatter(a0vals);
  const madA0 = medAbsDev(a0vals);
  console.log("\n  Per-galaxy a0 fit (McGaugh interpolation):");
  console.log("  Galaxies with good fit: " + perGalaxyA0.length + "/" + galaxyNames.length);
  console.log("  Mean log(a0)   = " + meanA0.toFixed(3) + " -> a0 = " + Math.pow(10, meanA0).toFixed(3));
  console.log("  Median log(a0) = " + medA0.toFixed(3) + " -> a0 = " + Math.pow(10, medA0).toFixed(3));
  console.log("  Scatter (RMS)  = " + scatterA0.toFixed(3) + " dex");
  console.log("  MAD            = " + madA0.toFixed(3) + " dex");
  console.log("  Known a0       = " + A0_KNOWN.toFixed(3) + " [log = " + Math.log10(A0_KNOWN).toFixed(3) + "]");
  
  const fractionalScatter = (Math.pow(10, scatterA0) - 1) * 100;
  console.log("  Fractional scatter: " + fractionalScatter.toFixed(1) + "%");
  console.log("  -> a0 is " + (scatterA0 < 0.15 ? "REMARKABLY TIGHT" : scatterA0 < 0.3 ? "FAIRLY CONSISTENT" : "VARIABLE") + " across galaxies");

  const sparcA0 = perGalaxyA0.filter(g => !g.name.startsWith("LT_"));
  const ltA0 = perGalaxyA0.filter(g => g.name.startsWith("LT_"));
  if (sparcA0.length > 0 && ltA0.length > 0) {
    const sparcMed = [...sparcA0.map(g => g.logA0)].sort((a, b) => a - b)[Math.floor(sparcA0.length / 2)];
    const ltMed = [...ltA0.map(g => g.logA0)].sort((a, b) => a - b)[Math.floor(ltA0.length / 2)];
    console.log("\n  Cross-dataset consistency:");
    console.log("    SPARC median a0 = " + Math.pow(10, sparcMed).toFixed(3) + " (" + sparcA0.length + " galaxies)");
    console.log("    LT median a0    = " + Math.pow(10, ltMed).toFixed(3) + " (" + ltA0.length + " galaxies)");
    console.log("    Offset: " + (sparcMed - ltMed).toFixed(3) + " dex");
  }
  
  console.log("\n  Distribution of per-galaxy a0:");
  const sorted = [...a0vals].sort((a, b) => a - b);
  const p10 = sorted[Math.floor(sorted.length * 0.1)];
  const p25 = sorted[Math.floor(sorted.length * 0.25)];
  const p50 = sorted[Math.floor(sorted.length * 0.5)];
  const p75 = sorted[Math.floor(sorted.length * 0.75)];
  const p90 = sorted[Math.floor(sorted.length * 0.9)];
  console.log("    10th percentile: " + Math.pow(10, p10).toFixed(3));
  console.log("    25th percentile: " + Math.pow(10, p25).toFixed(3));
  console.log("    50th (median):   " + Math.pow(10, p50).toFixed(3));
  console.log("    75th percentile: " + Math.pow(10, p75).toFixed(3));
  console.log("    90th percentile: " + Math.pow(10, p90).toFixed(3));
  console.log("    IQR factor:      " + (Math.pow(10, p75) / Math.pow(10, p25)).toFixed(2) + "x");
}

console.log("\n============================================");
console.log("STEP 3: COLLAPSE with normalized x = g_bar / g0");
console.log("============================================\n");

const g0_global = Math.pow(10, perGalaxyA0.length > 0 ? [...perGalaxyA0.map(g => g.logA0)].sort((a, b) => a - b)[Math.floor(perGalaxyA0.length / 2)] : Math.log10(A0_KNOWN));

console.log("  Using g0 = " + g0_global.toFixed(4) + " (median per-galaxy a0)\n");

function collapseScatter(points, a0) {
  const resids = [];
  for (const p of points) {
    const x = p.gBar / a0;
    const predRatio = 1 / (1 - Math.exp(-Math.sqrt(x)));
    const pred = Math.log10(predRatio);
    const obs = p.logRatio;
    resids.push(obs - pred);
  }
  return { rms: rmsScatter(resids), mad: medAbsDev(resids), n: resids.length, resids };
}

const collapseGlobal = collapseScatter(allPoints, g0_global);
const collapseKnown = collapseScatter(allPoints, A0_KNOWN);

console.log("  Collapse with median per-galaxy g0: rms=" + collapseGlobal.rms.toFixed(6) + " MAD=" + collapseGlobal.mad.toFixed(6));
console.log("  Collapse with known a0 (3.7032):    rms=" + collapseKnown.rms.toFixed(6) + " MAD=" + collapseKnown.mad.toFixed(6));

console.log("\n  Collapse quality in bins (using median g0):");
const normBins = [];
for (let i = 0; i < 25; i++) {
  const lo = -2.5 + i * 0.3;
  const hi = lo + 0.3;
  const inBin = allPoints.filter(p => Math.log10(p.gBar / g0_global) >= lo && Math.log10(p.gBar / g0_global) < hi);
  if (inBin.length < 10) continue;
  const logRatios = inBin.map(p => p.logRatio);
  const predRatios = inBin.map(p => {
    const x = p.gBar / g0_global;
    return Math.log10(1 / (1 - Math.exp(-Math.sqrt(x))));
  });
  const resids = logRatios.map((lr, j) => lr - predRatios[j]);
  normBins.push({ logX: (lo + hi) / 2, scatter: rmsScatter(resids), n: inBin.length, meanResid: resids.reduce((a, b) => a + b, 0) / resids.length });
}

console.log("  " + "log(g_bar/g0)".padEnd(15) + "scatter".padEnd(10) + "bias".padEnd(10) + "n".padEnd(6));
for (const b of normBins) {
  console.log("  " + b.logX.toFixed(2).padEnd(15) + b.scatter.toFixed(4).padEnd(10) + b.meanResid.toFixed(4).padEnd(10) + String(b.n).padEnd(6));
}

console.log("\n============================================");
console.log("STEP 4: Is g0 truly UNIVERSAL? The critical test");
console.log("============================================\n");

if (perGalaxyA0.length > 0) {
  console.log("  Comparing per-galaxy a0 with galaxy properties:\n");
  
  const galProps = [];
  for (const fit of perGalaxyA0) {
    const gInfo = galaxyLookup[fit.name] || galaxyLookup[fit.name.replace("LT_", "")];
    const pts = perGalaxy[fit.name];
    if (!pts || pts.length < 5) continue;
    const maxGbar = Math.max(...pts.map(p => p.gBar));
    const minGbar = Math.min(...pts.map(p => p.gBar));
    const medGbar = [...pts.map(p => p.logGbar)].sort((a, b) => a - b)[Math.floor(pts.length / 2)];
    galProps.push({
      name: fit.name,
      logA0: fit.logA0,
      logMaxGbar: Math.log10(maxGbar),
      logMinGbar: Math.log10(minGbar),
      medLogGbar: medGbar,
      nPts: pts.length,
      sigBar: gInfo ? gInfo.sigma_bar : null,
      isLT: fit.name.startsWith("LT_")
    });
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
    return syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
  }
  
  const withSig = galProps.filter(g => g.sigBar && g.sigBar > 0);
  
  const rMaxGbar = pearsonR(galProps.map(g => g.logA0), galProps.map(g => g.logMaxGbar));
  const rMinGbar = pearsonR(galProps.map(g => g.logA0), galProps.map(g => g.logMinGbar));
  const rMedGbar = pearsonR(galProps.map(g => g.logA0), galProps.map(g => g.medLogGbar));
  const rSigBar = withSig.length > 10 ? pearsonR(withSig.map(g => g.logA0), withSig.map(g => Math.log10(g.sigBar))) : NaN;
  
  console.log("  Correlation of per-galaxy log(a0) with:");
  console.log("    max(g_bar):     r = " + rMaxGbar.toFixed(3) + (Math.abs(rMaxGbar) > 0.3 ? " CORRELATED" : " weak"));
  console.log("    min(g_bar):     r = " + rMinGbar.toFixed(3) + (Math.abs(rMinGbar) > 0.3 ? " CORRELATED" : " weak"));
  console.log("    median(g_bar):  r = " + rMedGbar.toFixed(3) + (Math.abs(rMedGbar) > 0.3 ? " CORRELATED" : " weak"));
  if (!isNaN(rSigBar)) console.log("    log(Sigma_bar): r = " + rSigBar.toFixed(3) + (Math.abs(rSigBar) > 0.3 ? " CORRELATED" : " weak"));
  
  console.log("\n  If a0 is truly universal, ALL correlations should be weak (|r| < 0.2).");
  const allWeak = Math.abs(rMaxGbar) < 0.2 && Math.abs(rMinGbar) < 0.2 && Math.abs(rMedGbar) < 0.2;
  if (allWeak) {
    console.log("  >>> a0 is REMARKABLY INDEPENDENT of galaxy properties <<<");
  } else {
    console.log("  >>> a0 shows SOME dependence on galaxy properties (not perfectly universal) <<<");
    console.log("  This could be: real variation, or artifacts of limited radial coverage");
  }
}

console.log("\n============================================");
console.log("STEP 5: COLLAPSE QUALITY — per-galaxy vs universal a0");
console.log("============================================\n");

let rmsPerGal = 0, countPerGal = 0;
let rmsUniversal = 0, countUniversal = 0;
for (const fit of perGalaxyA0) {
  const pts = perGalaxy[fit.name];
  if (!pts || pts.length < 5) continue;
  const cPer = collapseScatter(pts, fit.a0);
  const cUni = collapseScatter(pts, g0_global);
  rmsPerGal += cPer.rms ** 2 * cPer.n;
  countPerGal += cPer.n;
  rmsUniversal += cUni.rms ** 2 * cUni.n;
  countUniversal += cUni.n;
}
rmsPerGal = Math.sqrt(rmsPerGal / countPerGal);
rmsUniversal = Math.sqrt(rmsUniversal / countUniversal);

console.log("  Weighted RMS with per-galaxy a0: " + rmsPerGal.toFixed(6));
console.log("  Weighted RMS with universal a0:  " + rmsUniversal.toFixed(6));
console.log("  Improvement from per-galaxy a0:  " + ((rmsUniversal - rmsPerGal) / rmsUniversal * 100).toFixed(2) + "%");
console.log("  -> " + ((rmsUniversal - rmsPerGal) / rmsUniversal * 100 < 10 ? "Universal a0 works almost as well as per-galaxy (a0 IS universal)" : "Per-galaxy a0 significantly better (a0 varies between galaxies)"));

console.log("\n============================================");
console.log("STEP 6: THE SMOKING GUN — a0 vs cosmological coincidence");
console.log("============================================\n");

const H0 = 67.4;
const c_km = 299792.458;
const cH0 = c_km * H0 / 3.086e19;
const a0_si = A0_KNOWN * 1e3 * 1e3 / (3.086e19);

console.log("  Known acceleration scales:");
console.log("    a0 (MOND/RAR) = " + a0_si.toExponential(3) + " m/s^2 = " + A0_KNOWN.toFixed(4) + " (km/s)^2/kpc");
console.log("    cH0            = " + (c_km * H0 / 1e3 / 3.086e19).toExponential(3) + " m/s^2");
console.log("    a0 / cH0       ~ " + (a0_si / (c_km * 1e3 * H0 / 3.086e19)).toFixed(2));
console.log("    c^2 * sqrt(Lambda/3) ~ a0  [cosmological coincidence]");
console.log("\n  This coincidence (a0 ~ cH0) is one of the most puzzling features.");
console.log("  If a0 were random feedback, it should correlate with galaxy properties.");
console.log("  If it's fundamental, it should be universal and match a cosmological scale.");

console.log("\n============================================");
console.log("FINAL VERDICT: Is there a hidden scale?");
console.log("============================================\n");

const a0Scatter = perGalaxyA0.length > 0 ? rmsScatter(perGalaxyA0.map(g => g.logA0)) : Infinity;
const a0MAD = perGalaxyA0.length > 0 ? medAbsDev(perGalaxyA0.map(g => g.logA0)) : Infinity;
const medA0 = perGalaxyA0.length > 0 ? Math.pow(10, [...perGalaxyA0.map(g => g.logA0)].sort((a, b) => a - b)[Math.floor(perGalaxyA0.length / 2)]) : A0_KNOWN;
const collapseImprovement = ((rmsUniversal - rmsPerGal) / rmsUniversal * 100);

console.log("  1. TRANSITION EXISTS: " + (bins.some(b => b.medianLogRatio > 0.3) && bins.some(b => b.medianLogRatio < 0.05) ? "YES" : "UNCLEAR"));
console.log("  2. TRANSITION SCALE: g0 = " + medA0.toFixed(3) + " (km/s)^2/kpc [= " + (medA0 * 1e6 / 3.086e19).toExponential(2) + " m/s^2]");
console.log("  3. UNIVERSAL (same for all galaxies):");
console.log("     Scatter: " + a0Scatter.toFixed(3) + " dex (RMS), " + a0MAD.toFixed(3) + " dex (MAD)");
console.log("     Per-galaxy vs universal: " + collapseImprovement.toFixed(1) + "% worse with universal a0");
console.log("  4. CONSISTENT WITH KNOWN a0: offset = " + Math.abs(Math.log10(medA0) - Math.log10(A0_KNOWN)).toFixed(3) + " dex from Milgrom/McGaugh value");
console.log("  5. COLLAPSE QUALITY: RMS = " + collapseKnown.rms.toFixed(4) + " dex around predicted curve");

const isUniversal = a0MAD < 0.2 && collapseImprovement < 15;
const transitionClear = bins.some(b => b.medianLogRatio > 0.3) && bins.some(b => b.medianLogRatio < 0.05);

if (transitionClear && isUniversal) {
  console.log("\n  >>> YES: There IS a hidden scale. It is REAL, SHARP, and REMARKABLY UNIVERSAL <<<");
  console.log("  >>> g0 ~ " + medA0.toFixed(2) + " (km/s)^2/kpc ~ 1.2e-10 m/s^2 ~ cH0 <<<");
  console.log("  >>> All " + perGalaxyA0.length + " galaxies collapse onto ONE curve when normalized by g0 <<<");
} else if (transitionClear) {
  console.log("\n  >>> TRANSITION EXISTS but a0 varies between galaxies (not perfectly universal) <<<");
} else {
  console.log("\n  >>> NO clear transition found <<<");
}

const plotData = {
  bins: bins.map(b => ({ logGbar: b.logGbarMid, medianLogRatio: b.medianLogRatio, scatter: b.scatter, n: b.n })),
  normBins: normBins.map(b => ({ logX: b.logX, scatter: b.scatter, n: b.n, meanResid: b.meanResid })),
  rawPoints: allPoints.filter((_, i) => i % 3 === 0).map(p => ({ logGbar: p.logGbar, logRatio: p.logRatio, galaxy: p.galaxy })),
};

const output = {
  timestamp: new Date().toISOString(),
  nPoints: allPoints.length,
  nGalaxies: galaxyNames.length,
  transitionScale: {
    g0_median: medA0,
    g0_known: A0_KNOWN,
    log_g0_median: Math.log10(medA0),
    log_g0_known: Math.log10(A0_KNOWN),
    offset_dex: Math.abs(Math.log10(medA0) - Math.log10(A0_KNOWN)),
    a0_si: medA0 * 1e6 / 3.086e19
  },
  universality: {
    scatterRMS: a0Scatter,
    scatterMAD: a0MAD,
    nGalaxiesWithFit: perGalaxyA0.length,
    perGalaxyVsUniversal: collapseImprovement,
    isUniversal: isUniversal
  },
  collapse: {
    rmsGlobal: collapseKnown.rms,
    madGlobal: collapseKnown.mad,
    rmsPerGalaxy: rmsPerGal,
    rmsUniversal: rmsUniversal
  },
  correlations: {
    withMaxGbar: typeof rMaxGbar !== 'undefined' ? rMaxGbar : null,
    withMinGbar: typeof rMinGbar !== 'undefined' ? rMinGbar : null,
    withMedGbar: typeof rMedGbar !== 'undefined' ? rMedGbar : null,
    withSigmaBar: typeof rSigBar !== 'undefined' && !isNaN(rSigBar) ? rSigBar : null
  },
  perGalaxyA0: perGalaxyA0.map(g => ({
    name: g.name,
    a0: g.a0,
    logA0: g.logA0,
    rms: g.rms,
    isLT: g.name.startsWith("LT_")
  })),
  plotData,
  verdict: transitionClear && isUniversal ? "UNIVERSAL_SCALE_FOUND" : transitionClear ? "TRANSITION_EXISTS_BUT_VARIABLE" : "NO_CLEAR_TRANSITION"
};

const outPath = path.join(__dirname, '..', 'public', 'transition-scale.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log("\nResults saved to: " + outPath);
