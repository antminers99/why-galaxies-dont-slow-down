#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "32.0.0";
function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(80)); }

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0(pts) {
  let lo = 2.0, hi = 5.0;
  for (let s = 0; s < 150; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gb = Math.pow(10, p.log_g_bar);
      c1 += (p.log_g_obs - Math.log10(mcgaughRAR(gb, Math.pow(10, m1)))) ** 2;
      c2 += (p.log_g_obs - Math.log10(mcgaughRAR(gb, Math.pow(10, m2)))) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  const logA0 = (lo + hi) / 2, a0 = Math.pow(10, logA0);
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.log_g_bar);
    ss += (p.log_g_obs - Math.log10(mcgaughRAR(gb, a0))) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length), ss, n: pts.length };
}

function predictSS(a0, pts) {
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gb, a0);
    ss += (p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10)) ** 2;
  }
  return ss;
}

function corrWith(x, y) {
  const n = x.length;
  if (n < 4) return 0;
  const mx = x.reduce((s, v) => s + v, 0) / n;
  const my = y.reduce((s, v) => s + v, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i]-mx)*(y[i]-my); sxx += (x[i]-mx)**2; syy += (y[i]-my)**2; }
  return sxy / Math.sqrt(sxx * syy + 1e-20);
}

function linReg(X, y) {
  const n = y.length, p = X[0] ? X[0].length : 0;
  if (p === 0) return { coefs: [], intercept: y.reduce((s,v)=>s+v,0)/n, residuals: y.map(v => v - y.reduce((s,v)=>s+v,0)/n) };
  const means = [];
  for (let j = 0; j < p; j++) { let s = 0; for (let i = 0; i < n; i++) s += X[i][j]; means.push(s / n); }
  const my = y.reduce((s,v)=>s+v,0) / n;
  const Xc = X.map(row => row.map((v,j) => v - means[j]));
  const yc = y.map(v => v - my);
  const XtX = Array.from({length:p}, () => Array(p).fill(0));
  const Xty = Array(p).fill(0);
  for (let i=0;i<n;i++) { for (let j=0;j<p;j++) { Xty[j]+=Xc[i][j]*yc[i]; for (let k=0;k<p;k++) XtX[j][k]+=Xc[i][j]*Xc[i][k]; } }
  for (let j=0;j<p;j++) XtX[j][j] += 1e-10;
  const A = XtX.map((r,i) => [...r, Xty[i]]);
  for (let j=0;j<p;j++) { let mx=j; for(let i=j+1;i<p;i++) if(Math.abs(A[i][j])>Math.abs(A[mx][j]))mx=i; [A[j],A[mx]]=[A[mx],A[j]]; for(let i=j+1;i<p;i++){const f=A[i][j]/A[j][j];for(let k=j;k<=p;k++)A[i][k]-=f*A[j][k];} }
  const b = Array(p).fill(0);
  for (let j=p-1;j>=0;j--) { b[j]=A[j][p]; for(let k=j+1;k<p;k++)b[j]-=A[j][k]*b[k]; b[j]/=A[j][j]; }
  const intercept = my - b.reduce((s,v,j) => s + v * means[j], 0);
  const residuals = y.map((yi, i) => yi - (intercept + b.reduce((s, c, j) => s + c * X[i][j], 0)));
  return { coefs: b, intercept, residuals };
}

function residualize(y, controls) {
  if (controls.length === 0) return [...y];
  const n = y.length;
  const X = [];
  for (let i = 0; i < n; i++) X.push(controls.map(c => c[i]));
  return linReg(X, y).residuals;
}

const cosmicWebData = {
  'NGC3726':     { cwType: 'filament', cwCode: 2, conf: 'high',   note: 'UMa cluster on Ursa Major filament' },
  'NGC3769':     { cwType: 'filament', cwCode: 2, conf: 'high',   note: 'UMa cluster on Ursa Major filament' },
  'NGC3893':     { cwType: 'filament', cwCode: 2, conf: 'high',   note: 'UMa cluster on Ursa Major filament' },
  'NGC4013':     { cwType: 'filament', cwCode: 2, conf: 'high',   note: 'UMa cluster on Ursa Major filament' },
  'NGC4100':     { cwType: 'filament', cwCode: 2, conf: 'high',   note: 'UMa cluster on Ursa Major filament' },
  'NGC4138':     { cwType: 'filament', cwCode: 2, conf: 'medium', note: 'UMa/CVn boundary on filament' },
  'NGC4157':     { cwType: 'filament', cwCode: 2, conf: 'high',   note: 'UMa cluster on filament' },
  'NGC4217':     { cwType: 'filament', cwCode: 2, conf: 'high',   note: 'UMa cluster on filament' },
  'UGC06973':    { cwType: 'filament', cwCode: 2, conf: 'medium', note: 'UMa periphery, filament' },
  'NGC4559':     { cwType: 'filament', cwCode: 2, conf: 'medium', note: 'Coma I cloud, Virgo-Coma filament' },

  'NGC5005':     { cwType: 'filament', cwCode: 2, conf: 'medium', note: 'CVn I cloud, CfA Great Wall' },
  'NGC5033':     { cwType: 'filament', cwCode: 2, conf: 'medium', note: 'CVn I cloud, CfA Great Wall' },
  'NGC5055':     { cwType: 'filament', cwCode: 2, conf: 'medium', note: 'M51 group, CVn spur' },
  'NGC5371':     { cwType: 'filament', cwCode: 2, conf: 'low',    note: 'CVn I cloud, filament' },
  'NGC5907':     { cwType: 'filament', cwCode: 2, conf: 'medium', note: 'Draco/NGC5866 group, filament' },

  'NGC2403':     { cwType: 'filament', cwCode: 2, conf: 'high',   note: 'M81 group, on M81-M82 filament to Virgo' },
  'NGC0891':     { cwType: 'filament', cwCode: 2, conf: 'medium', note: 'NGC1023 group, Perseus-Pisces filament' },
  'NGC7331':     { cwType: 'filament', cwCode: 2, conf: 'medium', note: 'Small group, Pegasus filament' },

  'NGC0024':     { cwType: 'wall', cwCode: 3, conf: 'medium', note: 'Sculptor group, part of Sculptor Wall / Southern Wall' },
  'NGC0289':     { cwType: 'wall', cwCode: 3, conf: 'medium', note: 'Sculptor Wall region' },
  'NGC7814':     { cwType: 'wall', cwCode: 3, conf: 'low',    note: 'Pegasus / local sheet' },

  'NGC2841':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Local Sheet / Supergalactic plane, relatively isolated' },
  'NGC2903':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Leo spur, between filaments' },
  'NGC3198':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'UMa periphery, sheet/inter-filament' },
  'NGC3521':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Leo region, sheet' },
  'NGC2683':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Lynx, local sheet' },
  'NGC1003':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Perseus region, sheet/inter-filament' },
  'NGC1090':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Cetus, sheet' },
  'NGC6015':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Draco field, sheet' },
  'NGC6503':     { cwType: 'void',  cwCode: 0, conf: 'high',   note: 'Local Void edge — classic void galaxy' },

  'NGC2915':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Centaurus region, likely sheet' },
  'NGC3741':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'CVn field, sheet' },
  'NGC1705':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Dorado/Eridanus region, sheet' },
  'UGC01281':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Perseus-Pisces foreground, sheet' },
  'UGC05721':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Leo region, sheet' },
  'UGC08490':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Bootes, sheet' },

  'NGC2955':     { cwType: 'filament', cwCode: 2, conf: 'low', note: 'D=98 Mpc, likely filament (CfA2 Great Wall region)' },
  'NGC2998':     { cwType: 'filament', cwCode: 2, conf: 'low', note: 'D=68, Leo/UMa filament' },
  'NGC6674':     { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Hercules, inter-filament' },
  'NGC6195':     { cwType: 'filament', cwCode: 2, conf: 'low', note: 'D=128, Hercules supercluster filament' },
  'ESO563-G021': { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Southern, sheet/inter-filament' },
  'F571-8':      { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'LSB D=53, sheet' },
  'IC4202':      { cwType: 'filament', cwCode: 2, conf: 'low', note: 'D=100, filament' },
  'NGC0801':     { cwType: 'filament', cwCode: 2, conf: 'low', note: 'D=81, Perseus-Pisces filament' },
  'UGC00128':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Pisces, sheet' },
  'UGC02885':    { cwType: 'void',  cwCode: 0, conf: 'low',    note: 'Very isolated giant, possible void/underdense' },
  'UGC02916':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Cam/Lynx, sheet' },
  'UGC02953':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Cam, sheet' },
  'UGC03205':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Lynx, sheet' },
  'UGC03546':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Cam/Lynx, sheet' },
  'UGC03580':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Cam, sheet' },
  'UGC06786':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'UMa outskirts, sheet' },
  'UGC06787':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'UMa outskirts, sheet' },
  'UGC08699':    { cwType: 'filament', cwCode: 2, conf: 'low', note: 'Bootes, CfA2 region' },
  'UGC09037':    { cwType: 'sheet', cwCode: 1, conf: 'low',    note: 'Bootes, sheet' },
  'UGC09133':    { cwType: 'filament', cwCode: 2, conf: 'low', note: 'Corona Borealis region, filament' },
};

const p11 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase11-sensitivity-lab.json'), 'utf8'));
const sparcAll = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));
const p25 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase25-group-membership.json'), 'utf8'));

const envLookup = {};
for (const g of p25.galaxyAssignments) envLookup[g.name] = g;

const galaxyData = [];

for (const gal of p11.galaxies) {
  const pts = gal.localProfile.filter(p =>
    isFinite(p.log_g_bar) && isFinite(p.log_g_obs) && p.log_g_obs > p.log_g_bar * 0.5 && p.r > 0
  );
  if (pts.length < 5) continue;
  const sp = sparcAll.find(s => s.name === gal.name);
  if (!sp || !sp.MHI || sp.MHI <= 0) continue;
  const L = sp.L36 || sp.L || 0;
  if (L <= 0) continue;

  const fit = fitA0(pts);
  const logMHI = Math.log10(sp.MHI);
  const n = pts.length;
  const logVflat = Math.log10(sp.Vflat || 100);

  const rArr = pts.map(p => p.r);
  const rMax = Math.max(...rArr);
  const rMid = rMax / 2;
  const innerPts = pts.filter(p => p.r <= rMid);
  const outerPts = pts.filter(p => p.r > rMid);

  function residualVariability(subset) {
    if (subset.length < 3) return 0;
    const xv = subset.map(p => p.log_g_bar);
    const yv = subset.map(p => p.log_g_obs);
    const mx = xv.reduce((s,v)=>s+v,0)/xv.length;
    const my = yv.reduce((s,v)=>s+v,0)/yv.length;
    let sxy = 0, sxx = 0;
    for (let i = 0; i < xv.length; i++) { sxy += (xv[i]-mx)*(yv[i]-my); sxx += (xv[i]-mx)**2; }
    const b = sxx > 0 ? sxy / sxx : 0;
    const a = my - b * mx;
    let ss = 0;
    for (let i = 0; i < xv.length; i++) ss += (yv[i] - (a + b * xv[i])) ** 2;
    return Math.sqrt(ss / xv.length);
  }

  const rcWiggliness = residualVariability(innerPts) + residualVariability(outerPts);

  const ev = envLookup[gal.name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  const cw = cosmicWebData[gal.name];
  if (!cw) continue;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms,
    n, logMHI, logVflat, rcWiggliness, envCode,
    D: sp.D,
    cwType: cw.cwType, cwCode: cw.cwCode, cwConf: cw.conf, cwNote: cw.note,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const cwArr = galaxyData.map(g => g.cwCode);
const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const vflatArr = galaxyData.map(g => g.logVflat);

const nVoid = galaxyData.filter(g => g.cwCode === 0).length;
const nSheet = galaxyData.filter(g => g.cwCode === 1).length;
const nFilament = galaxyData.filter(g => g.cwCode === 2).length;
const nWall = galaxyData.filter(g => g.cwCode === 3).length;
const nHigh = galaxyData.filter(g => g.cwConf === 'high').length;
const nMed = galaxyData.filter(g => g.cwConf === 'medium').length;
const nLow = galaxyData.filter(g => g.cwConf === 'low').length;

log("");
log("=".repeat(80));
log("  PHASE 32: COSMIC WEB ENVIRONMENT");
log("  Does large-scale cosmic web position predict a0 variation?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  COSMIC CONTEXT DOOR — FINAL VARIABLE:");
sep();
log("  1. redshift / distance:        2/6 FAIL (Malmquist artifact)");
log("  2. cosmic web environment:     ← TESTING NOW (LAST)");
log("");

log("  DEFINITION:");
log("  Cosmic web classification based on large-scale structure:");
log("    0 = void (underdense regions between filaments)");
log("    1 = sheet/wall (planar structures between filaments)");
log("    2 = filament (elongated overdensities connecting clusters)");
log("    3 = wall/node (dense walls or near cluster nodes)");
log("  Scale: ~10-100 Mpc, much larger than group membership (~1 Mpc)");
log("");

log("  DATA QUALITY:");
log("  ┌─────────────────────────────────────────────────────────────┐");
log("  │  Confidence   N    fraction                               │");
log("  ├─────────────────────────────────────────────────────────────┤");
log("  │  high          " + nHigh.toString().padStart(2) + "   " + (nHigh/N*100).toFixed(0) + "%".padEnd(40) + "│");
log("  │  medium        " + nMed.toString().padStart(2) + "   " + (nMed/N*100).toFixed(0) + "%".padEnd(40) + "│");
log("  │  low           " + nLow.toString().padStart(2) + "   " + (nLow/N*100).toFixed(0) + "%".padEnd(40) + "│");
log("  └─────────────────────────────────────────────────────────────┘");
log("");

log("  COSMIC WEB DISTRIBUTION:");
log("  ┌─────────────────────────────────────────────────────────────┐");
log("  │  Type        Code  N    <logA0>  SD                       │");
log("  ├─────────────────────────────────────────────────────────────┤");

for (const [label, code] of [['void', 0], ['sheet', 1], ['filament', 2], ['wall/node', 3]]) {
  const idx = galaxyData.filter(g => g.cwCode === code);
  if (idx.length === 0) { log("  │  " + label.padEnd(11) + code.toString().padStart(3) + "    0    N/A     N/A" + " ".repeat(21) + "│"); continue; }
  const la = idx.map(g => g.logA0);
  const mu = la.reduce((s,v)=>s+v,0)/la.length;
  const sd = la.length > 1 ? Math.sqrt(la.reduce((s,v)=>s+(v-mu)**2,0)/(la.length-1)) : 0;
  log("  │  " + label.padEnd(11) + code.toString().padStart(3) + idx.length.toString().padStart(5) + "    " + mu.toFixed(3) + "  " + sd.toFixed(3) + " ".repeat(21) + "│");
}
log("  └─────────────────────────────────────────────────────────────┘");
log("");

log("  COLLINEARITY CHECK:");
const rCwEnv = corrWith(cwArr, envArr);
log("  cwCode vs envCode: r=" + rCwEnv.toFixed(3));
log("  " + (Math.abs(rCwEnv) > 0.8 ? "VERY HIGH" : Math.abs(rCwEnv) > 0.6 ? "HIGH" : Math.abs(rCwEnv) > 0.4 ? "MODERATE" : "LOW") + " collinearity with local environment");
log("");

log("  PER-GALAXY (sorted by cwCode):");
const sorted = [...galaxyData].sort((a,b) => a.cwCode - b.cwCode || a.name.localeCompare(b.name));
for (const g of sorted) {
  log("    " + g.name.padEnd(16) + "cw=" + g.cwCode + " " + g.cwType.padEnd(9) +
    " env=" + g.envCode + " " + g.cwConf.padEnd(7) + " logA0=" + g.logA0.toFixed(3));
}
log("");

log("=".repeat(80));
log("  TEST 1: RAW CORRELATIONS");
log("=".repeat(80));
log("");

const rRaw = corrWith(cwArr, deltaA0);
const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);
log("  cwCode vs delta_a0: r=" + rRaw.toFixed(3) + " t=" + tRaw.toFixed(1));
log("");

const isFilament = galaxyData.map(g => g.cwCode >= 2 ? 1 : 0);
const rFil = corrWith(isFilament, deltaA0);
const tFil = Math.abs(rFil) * Math.sqrt(N-2) / Math.sqrt(1-rFil**2+1e-10);
log("  isFilament/wall (binary) vs delta_a0: r=" + rFil.toFixed(3) + " t=" + tFil.toFixed(1));
log("");

log("=".repeat(80));
log("  TEST 2: CONFOUNDER STRIPPING — THE CRITICAL TEST");
log("  Does cwCode survive after controlling envCode?");
log("=".repeat(80));
log("");

const controlSets = [
  { name: 'raw', vals: [] },
  { name: '|MHI', vals: [mhiArr] },
  { name: '|MHI+Vflat', vals: [mhiArr, vflatArr] },
  { name: '|MHI+wig', vals: [mhiArr, wigArr] },
  { name: '|MHI+wig+env', vals: [mhiArr, wigArr, envArr] },
  { name: '|env', vals: [envArr] },
  { name: '|MHI+env', vals: [mhiArr, envArr] },
  { name: '|MHI+Vflat+wig+env', vals: [mhiArr, vflatArr, wigArr, envArr] },
];

for (const [label, arr] of [['cwCode', cwArr], ['isFilament', isFilament]]) {
  log("  " + label + " vs delta_a0:");
  log("  ┌──────────────────────────────────────────────────────────────┐");
  log("  │  Controls                 r_partial   |t|    status       │");
  log("  ├──────────────────────────────────────────────────────────────┤");
  for (const cs of controlSets) {
    const residY = residualize(deltaA0, cs.vals);
    const residX = residualize(arr, cs.vals);
    const r = corrWith(residX, residY);
    const df = N - 2 - cs.vals.length;
    const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
    const status = t >= 2.0 ? "SURVIVES" : t >= 1.65 ? "MARGINAL" : "FAILS";
    log("  │  " + cs.name.padEnd(24) + r.toFixed(3).padStart(8) + t.toFixed(1).padStart(7) + "    " + status.padEnd(10) + "│");
  }
  log("  └──────────────────────────────────────────────────────────────┘");
  log("");
}

log("=".repeat(80));
log("  TEST 3: cwCode AFTER envCode — ABSORPTION CHECK");
log("=".repeat(80));
log("");

const residAfterEnv = residualize(deltaA0, [mhiArr, wigArr, envArr]);
const rCwAfterEnv = corrWith(residualize(cwArr, [mhiArr, wigArr, envArr]), residAfterEnv);
const tCwAfterEnv = Math.abs(rCwAfterEnv) * Math.sqrt(N-5) / Math.sqrt(1-rCwAfterEnv**2+1e-10);

const residAfterCw = residualize(deltaA0, [mhiArr, wigArr, cwArr]);
const rEnvAfterCw = corrWith(residualize(envArr, [mhiArr, wigArr, cwArr]), residAfterCw);
const tEnvAfterCw = Math.abs(rEnvAfterCw) * Math.sqrt(N-5) / Math.sqrt(1-rEnvAfterCw**2+1e-10);

log("  cwCode AFTER envCode+MHI+wig: r=" + rCwAfterEnv.toFixed(3) + " t=" + tCwAfterEnv.toFixed(1) +
  " → " + (tCwAfterEnv >= 1.65 ? "ADDS beyond envCode" : "ABSORBED"));
log("  envCode AFTER cwCode+MHI+wig: r=" + rEnvAfterCw.toFixed(3) + " t=" + tEnvAfterCw.toFixed(1) +
  " → " + (tEnvAfterCw >= 1.65 ? "envCode STILL adds" : "envCode ABSORBED"));
log("");

if (tCwAfterEnv < 1.65 && tEnvAfterCw >= 1.65) {
  log("  → envCode ABSORBS cosmic web. Local environment is deeper.");
} else if (tCwAfterEnv >= 1.65 && tEnvAfterCw < 1.65) {
  log("  → Cosmic web ABSORBS envCode! Large-scale structure is deeper.");
} else if (tCwAfterEnv >= 1.65 && tEnvAfterCw >= 1.65) {
  log("  → BOTH carry independent info at different scales.");
} else {
  log("  → Neither adds when the other is present (same information).");
}
log("");

log("=".repeat(80));
log("  TEST 4: PERMUTATION TESTS");
log("=".repeat(80));
log("");

const NPERMS = 10000;
for (const [label, arr, ctrls] of [
  ['cwCode|raw', cwArr, []],
  ['cwCode|MHI', cwArr, [mhiArr]],
  ['cwCode|MHI+wig', cwArr, [mhiArr, wigArr]],
  ['cwCode|MHI+wig+env', cwArr, [mhiArr, wigArr, envArr]],
  ['cwCode|env', cwArr, [envArr]],
]) {
  const residY = residualize(deltaA0, ctrls);
  const residX = residualize(arr, ctrls);
  const obsR = Math.abs(corrWith(residX, residY));
  let cnt = 0;
  const sh = [...residY];
  for (let p = 0; p < NPERMS; p++) {
    for (let i = sh.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [sh[i], sh[j]] = [sh[j], sh[i]];
    }
    if (Math.abs(corrWith(residX, sh)) >= obsR) cnt++;
  }
  log("  " + label.padEnd(22) + " |r|=" + obsR.toFixed(3) + " perm_p=" + (cnt/NPERMS).toFixed(4) +
    (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "**" : ""));
}
log("");

log("=".repeat(80));
log("  TEST 5: LOO CROSS-VALIDATION");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'envCode': envArr,
  'cwCode': cwArr, 'logVflat': vflatArr,
};

function getVal(gIdx, fname) { return allFeatures[fname][gIdx]; }

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_cw: cwCode', features: ['cwCode'] },
  { name: 'M_env', features: ['envCode'] },
  { name: 'M_MHI', features: ['logMHI'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
  { name: 'M_MHI+wig+cw', features: ['logMHI', 'rcWiggliness', 'cwCode'] },
  { name: 'M_MHI+wig+env', features: ['logMHI', 'rcWiggliness', 'envCode'] },
  { name: 'M_MHI+wig+env+cw', features: ['logMHI', 'rcWiggliness', 'envCode', 'cwCode'] },
  { name: 'M6: Per-galaxy', features: ['__free__'] },
];

const cvResults = [];
for (const model of modelDefs) {
  let totalSS = 0, totalN = 0;
  if (model.features[0] === '__free__') {
    for (let i = 0; i < N; i++) { totalSS += predictSS(galaxyData[i].a0, galaxyData[i].pts); totalN += galaxyData[i].pts.length; }
    cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), k: N }); continue;
  }
  if (model.features.length === 0) {
    for (let i = 0; i < N; i++) { const trainY = allLogA0.filter((_, j) => j !== i); const mu = trainY.reduce((s,v) => s+v, 0) / trainY.length; totalSS += predictSS(Math.pow(10, mu), galaxyData[i].pts); totalN += galaxyData[i].pts.length; }
    cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), k: 1 }); continue;
  }
  for (let i = 0; i < N; i++) {
    const trainIdx = [...Array(N).keys()].filter(j => j !== i);
    const trainX = trainIdx.map(j => model.features.map(fn => getVal(j, fn)));
    const trainY = trainIdx.map(j => allLogA0[j]);
    const reg = linReg(trainX, trainY);
    const testX = model.features.map(fn => getVal(i, fn));
    let predLogA0 = reg.intercept + reg.coefs.reduce((s, c, j) => s + c * testX[j], 0);
    predLogA0 = Math.max(2.5, Math.min(4.5, predLogA0));
    totalSS += predictSS(Math.pow(10, predLogA0), galaxyData[i].pts);
    totalN += galaxyData[i].pts.length;
  }
  cvResults.push({ name: model.name, cvRMS: Math.sqrt(totalSS / totalN), k: model.features.length + 1 });
}

const m0rms = cvResults[0].cvRMS;
const m6rms = cvResults[cvResults.length - 1].cvRMS;
const gap = m0rms - m6rms;

log("  ┌───────────────────────────────────────────────────────────────────────────┐");
log("  │  Model                    k    CV-RMS   vs M0   gap-closed             │");
log("  ├───────────────────────────────────────────────────────────────────────────┤");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const gc = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  │  " + r.name.padEnd(24) + r.k.toString().padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(8) + (gc + "%").padStart(13) + "            │");
}
log("  └───────────────────────────────────────────────────────────────────────────┘");
log("");

const fnm = (name) => cvResults.find(r => r.name === name);
const mhiWigGap = gap > 0 ? (m0rms - fnm('M_MHI+wig').cvRMS) / gap * 100 : 0;
const mhiWigEnvGap = gap > 0 ? (m0rms - fnm('M_MHI+wig+env').cvRMS) / gap * 100 : 0;
const mhiWigCwGap = gap > 0 ? (m0rms - fnm('M_MHI+wig+cw').cvRMS) / gap * 100 : 0;
const mhiWigEnvCwGap = gap > 0 ? (m0rms - fnm('M_MHI+wig+env+cw').cvRMS) / gap * 100 : 0;
const cwAloneGap = gap > 0 ? (m0rms - fnm('M_cw: cwCode').cvRMS) / gap * 100 : 0;
const envAloneGap = gap > 0 ? (m0rms - fnm('M_env').cvRMS) / gap * 100 : 0;

log("  KEY COMPARISONS:");
log("    cwCode alone:               " + cwAloneGap.toFixed(1) + "%");
log("    envCode alone:              " + envAloneGap.toFixed(1) + "%");
log("    MHI + wig:                  " + mhiWigGap.toFixed(1) + "%");
log("    MHI + wig + cw:             " + mhiWigCwGap.toFixed(1) + "%");
log("    MHI + wig + env:            " + mhiWigEnvGap.toFixed(1) + "%");
log("    MHI + wig + env + cw:       " + mhiWigEnvCwGap.toFixed(1) + "%");
log("");

const cwAdds = mhiWigCwGap > mhiWigGap + 1;
const cwAddsToEnv = mhiWigEnvCwGap > mhiWigEnvGap + 1;
log("  cw adds beyond MHI+wig? " + (cwAdds ? "YES (+" + (mhiWigCwGap-mhiWigGap).toFixed(1) + "%)" : "NO"));
log("  cw adds beyond MHI+wig+env? " + (cwAddsToEnv ? "YES" : "NO"));
log("");

log("=".repeat(80));
log("  PHASE 32 — STRICT CRITERIA");
log("=".repeat(80));
log("");

const rAfterMHI = corrWith(residualize(cwArr, [mhiArr]), residualize(deltaA0, [mhiArr]));
const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);
const rAfterMHIWig = corrWith(residualize(cwArr, [mhiArr, wigArr]), residualize(deltaA0, [mhiArr, wigArr]));
const tAfterMHIWig = Math.abs(rAfterMHIWig) * Math.sqrt(N-4) / Math.sqrt(1-rAfterMHIWig**2+1e-10);

const criteria = [
  { name: 'Clear raw |t|>=2', pass: tRaw >= 2.0, val: 'r=' + rRaw.toFixed(3) + ' t=' + tRaw.toFixed(1) },
  { name: 'Independent of MHI', pass: tAfterMHI >= 1.65, val: 'r=' + rAfterMHI.toFixed(3) + ' t=' + tAfterMHI.toFixed(1) },
  { name: 'Independent of MHI+wig', pass: tAfterMHIWig >= 1.65, val: 'r=' + rAfterMHIWig.toFixed(3) + ' t=' + tAfterMHIWig.toFixed(1) },
  { name: 'LOO-CV > M0', pass: cwAloneGap > 2, val: cwAloneGap.toFixed(1) + '%' },
  { name: 'Adds to MHI+wig', pass: cwAdds, val: (mhiWigCwGap - mhiWigGap).toFixed(1) + '%' },
  { name: 'Adds beyond MHI+wig+env', pass: cwAddsToEnv, val: (mhiWigEnvCwGap - mhiWigEnvGap).toFixed(1) + '%' },
];

const nPassed = criteria.filter(c => c.pass).length;

log("  ┌───────────────────────────────────────────────────────────────────┐");
log("  │  Criterion                   Pass?    Value                     │");
log("  ├───────────────────────────────────────────────────────────────────┤");
for (const c of criteria) {
  log("  │  " + c.name.padEnd(29) + (c.pass ? " YES " : " NO  ") + "  " + c.val.padEnd(27) + "│");
}
log("  ├───────────────────────────────────────────────────────────────────┤");
log("  │  TOTAL: " + nPassed + "/6 criteria met" + " ".repeat(39) + "│");
log("  └───────────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  PHASE 32 — VERDICT");
log("=".repeat(80));
log("");

if (nPassed >= 4) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  COSMIC WEB shows significant independent signal.                      ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nPassed >= 2) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  COSMIC WEB shows PARTIAL signal.                                      ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  COSMIC WEB FAILS.                                                     ║");
  log("  ║  Large-scale structure does NOT independently predict a0.              ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");

log("=".repeat(80));
log("  COSMIC CONTEXT DOOR — FINAL SCORECARD (DOOR NOW CLOSED)");
log("=".repeat(80));
log("");
log("  ┌──────────────────────────────────────────────────────────────────┐");
log("  │  Variable                    Result   Adds to best model?     │");
log("  ├──────────────────────────────────────────────────────────────────┤");
log("  │  redshift / distance         2/6 FAIL  NO (Malmquist bias)    │");
log("  │  cosmic web environment      " + nPassed + "/6 " + (nPassed >= 4 ? "PASS" : "FAIL").padEnd(6) + (cwAddsToEnv ? "YES" : "NO").padEnd(25) + "│");
log("  └──────────────────────────────────────────────────────────────────┘");
log("");
log("  COSMIC CONTEXT DOOR: OFFICIALLY CLOSED.");
log("");

log("=".repeat(80));
log("  PROJECT STATUS — ALL CLOSED DOORS");
log("=".repeat(80));
log("");
log("  DOOR 1: GALAXY HISTORY ← CLOSED");
log("    Survivor: rcWiggliness (RC shape complexity)");
log("    All intrinsic history variables (sSFR, age, metallicity, SFH, merger) failed");
log("");
log("  DOOR 2: ENVIRONMENT & NEIGHBORS ← CLOSED");
log("    Survivor: envCode (group membership)");
log("    satellite/central, dNN, tidal, density, group-distance all failed/redundant");
log("    Signal is CATEGORICAL (threshold effect, not gradient)");
log("");
log("  DOOR 3: COSMIC CONTEXT ← NOW CLOSED");
log("    redshift/distance: FAIL (Malmquist artifact)");
log("    cosmic web: " + nPassed + "/6 " + (nPassed >= 4 ? "PASS" : "FAIL"));
log("");
log("  BEST MODEL (ENTIRE PROJECT):");
log("  MHI + rcWiggliness + envCode = 14.9% gap closure");
log("  tau = 0.22 dex (robust, unchanged)");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  cosmicWebDistribution: { void: nVoid, sheet: nSheet, filament: nFilament, wall: nWall },
  dataQuality: { high: nHigh, medium: nMed, low: nLow },
  cwVsEnvCollinearity: +rCwEnv.toFixed(3),
  rawCorrelation: { r: +rRaw.toFixed(3), t: +tRaw.toFixed(1) },
  afterMHI: { r: +rAfterMHI.toFixed(3), t: +tAfterMHI.toFixed(1) },
  afterMHIWig: { r: +rAfterMHIWig.toFixed(3), t: +tAfterMHIWig.toFixed(1) },
  cwAfterEnv: { r: +rCwAfterEnv.toFixed(3), t: +tCwAfterEnv.toFixed(1) },
  envAfterCw: { r: +rEnvAfterCw.toFixed(3), t: +tEnvAfterCw.toFixed(1) },
  looGapClosed: {
    cw: +cwAloneGap.toFixed(1), env: +envAloneGap.toFixed(1),
    mhiWig: +mhiWigGap.toFixed(1), mhiWigCw: +mhiWigCwGap.toFixed(1),
    mhiWigEnv: +mhiWigEnvGap.toFixed(1), mhiWigEnvCw: +mhiWigEnvCwGap.toFixed(1),
  },
  successCriteria: { passed: nPassed, total: 6 },
  verdict: nPassed >= 4 ? 'SIGNIFICANT' : nPassed >= 2 ? 'PARTIAL' : 'FAILED',
  doorClosed: true,
};

fs.writeFileSync(path.join(__dirname, '../public/phase32-cosmic-web.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase32-cosmic-web.json");
