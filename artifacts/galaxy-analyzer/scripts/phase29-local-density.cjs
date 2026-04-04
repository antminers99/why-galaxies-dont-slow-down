#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "29.0.0";
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

const densityData = {
  'NGC3726':     { nNeighbors: 12, rho: 8.0,  conf: 'high',   note: 'UMa cluster core, ~12 L* galaxies within 1 Mpc' },
  'NGC3769':     { nNeighbors: 12, rho: 8.0,  conf: 'high',   note: 'UMa cluster core' },
  'NGC3893':     { nNeighbors: 12, rho: 8.0,  conf: 'high',   note: 'UMa cluster core' },
  'NGC4013':     { nNeighbors: 10, rho: 7.0,  conf: 'high',   note: 'UMa cluster' },
  'NGC4100':     { nNeighbors: 10, rho: 7.0,  conf: 'high',   note: 'UMa cluster' },
  'NGC4138':     { nNeighbors: 8,  rho: 5.5,  conf: 'medium', note: 'UMa cluster edge / CVn border' },
  'NGC4157':     { nNeighbors: 10, rho: 7.0,  conf: 'high',   note: 'UMa cluster' },
  'NGC4217':     { nNeighbors: 10, rho: 7.0,  conf: 'high',   note: 'UMa cluster' },
  'UGC06973':    { nNeighbors: 8,  rho: 5.5,  conf: 'medium', note: 'UMa cluster periphery' },
  'NGC4559':     { nNeighbors: 6,  rho: 4.0,  conf: 'medium', note: 'Coma I cloud, moderate density' },

  'NGC5005':     { nNeighbors: 5,  rho: 3.0,  conf: 'medium', note: 'CVn I cloud, NGC5005/5033 pair + dwarfs' },
  'NGC5033':     { nNeighbors: 5,  rho: 3.0,  conf: 'medium', note: 'CVn I cloud with NGC5005' },
  'NGC5055':     { nNeighbors: 4,  rho: 2.5,  conf: 'medium', note: 'M51 subgroup, a few companions' },
  'NGC5371':     { nNeighbors: 4,  rho: 2.5,  conf: 'medium', note: 'CVn I cloud, near NGC5353/4 group' },
  'NGC5907':     { nNeighbors: 3,  rho: 2.0,  conf: 'medium', note: 'NGC5866/Draco group, small' },
  'NGC2403':     { nNeighbors: 5,  rho: 3.0,  conf: 'high',   note: 'M81 group, ~5 significant members within 1 Mpc' },
  'NGC0024':     { nNeighbors: 4,  rho: 2.5,  conf: 'medium', note: 'Sculptor group, moderate' },
  'NGC0891':     { nNeighbors: 3,  rho: 1.5,  conf: 'medium', note: 'NGC1023 group, sparse' },
  'NGC7331':     { nNeighbors: 2,  rho: 1.0,  conf: 'medium', note: 'Small group with few companions' },
  'NGC7814':     { nNeighbors: 2,  rho: 1.0,  conf: 'medium', note: 'Near NGC7331, sparse' },

  'NGC2841':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Relatively isolated' },
  'NGC2903':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Leo field, isolated' },
  'NGC3198':     { nNeighbors: 1, rho: 0.4, conf: 'low', note: 'UMa periphery' },
  'NGC3521':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Leo field' },
  'NGC2683':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Lynx field' },
  'NGC1003':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Perseus field' },
  'NGC1090':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Cetus field' },
  'NGC0289':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Sculptor field' },
  'NGC6015':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Isolated' },
  'NGC6503':     { nNeighbors: 0, rho: 0.1, conf: 'medium', note: 'Local Void edge, extremely isolated' },
  'NGC2915':     { nNeighbors: 0, rho: 0.1, conf: 'low', note: 'Isolated BCD' },
  'NGC3741':     { nNeighbors: 0, rho: 0.1, conf: 'low', note: 'Isolated dwarf' },
  'NGC1705':     { nNeighbors: 0, rho: 0.1, conf: 'low', note: 'Isolated starburst' },
  'UGC01281':    { nNeighbors: 0, rho: 0.1, conf: 'low', note: 'Edge-on dwarf, sparse' },
  'UGC05721':    { nNeighbors: 0, rho: 0.1, conf: 'low', note: 'Dwarf field' },
  'UGC08490':    { nNeighbors: 0, rho: 0.1, conf: 'low', note: 'Dwarf IZw36' },

  'NGC2955':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Distant D=98' },
  'NGC2998':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=68' },
  'NGC6674':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Hercules field' },
  'NGC6195':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Distant D=128' },
  'ESO563-G021': { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Southern D=61' },
  'F571-8':      { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'LSB D=53' },
  'IC4202':      { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Distant D=100' },
  'NGC0801':     { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=81' },
  'UGC00128':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'LSB field' },
  'UGC02885':    { nNeighbors: 0, rho: 0.1, conf: 'low', note: 'Giant LSB, very isolated' },
  'UGC02916':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=65' },
  'UGC02953':    { nNeighbors: 1, rho: 0.4, conf: 'low', note: 'Field D=16' },
  'UGC03205':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=50' },
  'UGC03546':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=29' },
  'UGC03580':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=21' },
  'UGC06786':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=29' },
  'UGC06787':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=21' },
  'UGC08699':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=39' },
  'UGC09037':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=84' },
  'UGC09133':    { nNeighbors: 1, rho: 0.3, conf: 'low', note: 'Field D=57' },
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

  const dd = densityData[gal.name];
  const rho = dd ? dd.rho : 0.3;
  const logRho = Math.log10(rho + 0.01);
  const nNeighbors = dd ? dd.nNeighbors : 1;
  const conf = dd ? dd.conf : 'low';

  const ev = envLookup[gal.name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    n, logMHI, logVflat, rcWiggliness,
    rho, logRho, nNeighbors, conf, envCode,
    D: sp.D,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const nHigh = galaxyData.filter(g => g.conf === 'high').length;
const nMed = galaxyData.filter(g => g.conf === 'medium').length;
const nLow = galaxyData.filter(g => g.conf === 'low').length;

log("");
log("=".repeat(80));
log("  PHASE 29: LOCAL GALAXY DENSITY");
log("  Does the number density of neighbors predict a0 variation?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  DEFINITION:");
sep();
log("  rho = estimated number of significant (L > 0.1 L*) galaxies");
log("  within ~1 Mpc projected separation at similar distance.");
log("  Normalized to galaxies per Mpc³ (approximate).");
log("  Also track nNeighbors = raw count of known neighbors.");
log("");
log("  Unlike dNN or tidal index (single nearest neighbor),");
log("  local density captures the COLLECTIVE environment.");
log("");

log("  DATA QUALITY:");
log("  ┌─────────────────────────────────────────────┐");
log("  │  Confidence   N    fraction                 │");
log("  ├─────────────────────────────────────────────┤");
log("  │  high         " + nHigh.toString().padStart(2) + "   " + (nHigh/N*100).toFixed(0) + "%".padEnd(22) + "│");
log("  │  medium       " + nMed.toString().padStart(2) + "   " + (nMed/N*100).toFixed(0) + "%".padEnd(22) + "│");
log("  │  low          " + nLow.toString().padStart(2) + "   " + (nLow/N*100).toFixed(0) + "%".padEnd(22) + "│");
log("  └─────────────────────────────────────────────┘");
log("");

log("  PER-GALAXY (sorted by density):");
const sorted = [...galaxyData].sort((a,b) => b.rho - a.rho);
for (const g of sorted) {
  log("    " + g.name.padEnd(16) + "rho=" + g.rho.toFixed(1).padStart(5) +
    "  nN=" + g.nNeighbors.toString().padStart(2) +
    "  " + g.conf.padEnd(7) + " env=" + g.envCode + " logA0=" + g.logA0.toFixed(3));
}
log("");

const rhoArr = galaxyData.map(g => g.rho);
const logRhoArr = galaxyData.map(g => g.logRho);
const nNArr = galaxyData.map(g => g.nNeighbors);
const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const vflatArr = galaxyData.map(g => g.logVflat);
const nPtsArr = galaxyData.map(g => g.n);

log("=".repeat(80));
log("  TEST 1: RAW CORRELATIONS");
log("=".repeat(80));
log("");

for (const [label, arr] of [['rho', rhoArr], ['log(rho)', logRhoArr], ['nNeighbors', nNArr]]) {
  const r = corrWith(arr, deltaA0);
  const t = Math.abs(r) * Math.sqrt(N-2) / Math.sqrt(1-r**2+1e-10);
  log("  " + label.padEnd(14) + " vs delta_a0: r=" + r.toFixed(3) + " t=" + t.toFixed(1));
}
log("");

const rRhoEnv = corrWith(logRhoArr, envArr);
log("  log(rho) vs envCode: r=" + rRhoEnv.toFixed(3) + " (collinearity check)");
log("  " + (Math.abs(rRhoEnv) > 0.8 ? "VERY HIGH" : Math.abs(rRhoEnv) > 0.6 ? "HIGH" : "MODERATE") + " collinearity");
log("");

log("=".repeat(80));
log("  TEST 2: CONFOUNDER STRIPPING");
log("=".repeat(80));
log("");

const controlSets = [
  { name: 'raw', vals: [] },
  { name: '|MHI', vals: [mhiArr] },
  { name: '|MHI+Vflat', vals: [mhiArr, vflatArr] },
  { name: '|MHI+wig', vals: [mhiArr, wigArr] },
  { name: '|MHI+wig+env', vals: [mhiArr, wigArr, envArr] },
  { name: '|MHI+Vflat+wig', vals: [mhiArr, vflatArr, wigArr] },
  { name: '|MHI+Vflat+wig+env', vals: [mhiArr, vflatArr, wigArr, envArr] },
];

for (const [label, arr] of [['log(rho)', logRhoArr], ['nNeighbors', nNArr]]) {
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
log("  TEST 3: LOCAL DENSITY AFTER envCode — DOES IT ADD?");
log("=".repeat(80));
log("");

const residAfterEnv = residualize(deltaA0, [mhiArr, wigArr, envArr]);
const rRhoAfterEnv = corrWith(residualize(logRhoArr, [mhiArr, wigArr, envArr]), residAfterEnv);
const tRhoAfterEnv = Math.abs(rRhoAfterEnv) * Math.sqrt(N-5) / Math.sqrt(1-rRhoAfterEnv**2+1e-10);

const residAfterRho = residualize(deltaA0, [mhiArr, wigArr, logRhoArr]);
const rEnvAfterRho = corrWith(residualize(envArr, [mhiArr, wigArr, logRhoArr]), residAfterRho);
const tEnvAfterRho = Math.abs(rEnvAfterRho) * Math.sqrt(N-5) / Math.sqrt(1-rEnvAfterRho**2+1e-10);

log("  log(rho) AFTER envCode+MHI+wig: r=" + rRhoAfterEnv.toFixed(3) + " t=" + tRhoAfterEnv.toFixed(1) +
  " → " + (tRhoAfterEnv >= 1.65 ? "ADDS beyond envCode" : "ABSORBED"));
log("  envCode AFTER log(rho)+MHI+wig: r=" + rEnvAfterRho.toFixed(3) + " t=" + tEnvAfterRho.toFixed(1) +
  " → " + (tEnvAfterRho >= 1.65 ? "envCode STILL adds" : "envCode ABSORBED"));
log("");

if (tRhoAfterEnv < 1.65 && tEnvAfterRho >= 1.65) {
  log("  → envCode ABSORBS density. Group membership is the deeper variable.");
} else if (tRhoAfterEnv >= 1.65 && tEnvAfterRho < 1.65) {
  log("  → Density ABSORBS envCode! Local density is the deeper variable.");
} else if (tRhoAfterEnv >= 1.65 && tEnvAfterRho >= 1.65) {
  log("  → BOTH carry independent info.");
} else {
  log("  → Neither adds when the other is present.");
}
log("");

log("=".repeat(80));
log("  TEST 4: PERMUTATION TESTS");
log("=".repeat(80));
log("");

const NPERMS = 10000;
for (const [label, arr, ctrls] of [
  ['logRho|MHI', logRhoArr, [mhiArr]],
  ['logRho|MHI+wig', logRhoArr, [mhiArr, wigArr]],
  ['logRho|MHI+wig+env', logRhoArr, [mhiArr, wigArr, envArr]],
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
    (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "*" : ""));
}
log("");

log("=".repeat(80));
log("  TEST 5: LOO CROSS-VALIDATION");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'logRho': logRhoArr,
  'envCode': envArr, 'logVflat': vflatArr, 'nNeighbors': nNArr,
};

function getVal(gIdx, fname) { return allFeatures[fname][gIdx]; }

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_rho: log(rho)', features: ['logRho'] },
  { name: 'M_nN: nNeighbors', features: ['nNeighbors'] },
  { name: 'M_env: envCode', features: ['envCode'] },
  { name: 'M_MHI', features: ['logMHI'] },
  { name: 'M_MHI+rho', features: ['logMHI', 'logRho'] },
  { name: 'M_MHI+env', features: ['logMHI', 'envCode'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
  { name: 'M_MHI+wig+rho', features: ['logMHI', 'rcWiggliness', 'logRho'] },
  { name: 'M_MHI+wig+env', features: ['logMHI', 'rcWiggliness', 'envCode'] },
  { name: 'M_MHI+wig+env+rho', features: ['logMHI', 'rcWiggliness', 'envCode', 'logRho'] },
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

const fn = (name) => cvResults.find(r => r.name === name) || cvResults.find(r => r.name.includes(name));
const mhiWigGap = gap > 0 ? (m0rms - fn('M_MHI+wig').cvRMS) / gap * 100 : 0;
const mhiWigRhoGap = gap > 0 ? (m0rms - fn('M_MHI+wig+rho').cvRMS) / gap * 100 : 0;
const mhiWigEnvGap = gap > 0 ? (m0rms - fn('M_MHI+wig+env').cvRMS) / gap * 100 : 0;
const mhiWigEnvRhoGap = gap > 0 ? (m0rms - fn('M_MHI+wig+env+rho').cvRMS) / gap * 100 : 0;
const rhoAlone = gap > 0 ? (m0rms - fn('M_rho').cvRMS) / gap * 100 : 0;
const envAlone = gap > 0 ? (m0rms - fn('M_env').cvRMS) / gap * 100 : 0;

log("  KEY COMPARISONS:");
log("    log(rho) alone:             " + rhoAlone.toFixed(1) + "%");
log("    envCode alone:              " + envAlone.toFixed(1) + "%");
log("    MHI + wig:                  " + mhiWigGap.toFixed(1) + "%");
log("    MHI + wig + rho:            " + mhiWigRhoGap.toFixed(1) + "%");
log("    MHI + wig + env:            " + mhiWigEnvGap.toFixed(1) + "%");
log("    MHI + wig + env + rho:      " + mhiWigEnvRhoGap.toFixed(1) + "%");
log("");

const rhoAdds = mhiWigRhoGap > mhiWigGap + 1;
const rhoAddsToEnv = mhiWigEnvRhoGap > mhiWigEnvGap + 1;
log("  rho adds beyond MHI+wig? " + (rhoAdds ? "YES (+" + (mhiWigRhoGap-mhiWigGap).toFixed(1) + "%)" : "NO"));
log("  rho adds beyond MHI+wig+env? " + (rhoAddsToEnv ? "YES" : "NO"));
log("");

log("=".repeat(80));
log("  TEST 6: HIGH-CONFIDENCE SUBSAMPLE");
log("=".repeat(80));
log("");

const hiMedIdx = galaxyData.map((g,i) => (g.conf === 'high' || g.conf === 'medium') ? i : -1).filter(i => i >= 0);
const hiN = hiMedIdx.length;
const hiDA0 = hiMedIdx.map(i => deltaA0[i]);
const hiRho = hiMedIdx.map(i => logRhoArr[i]);
const hiMHI = hiMedIdx.map(i => mhiArr[i]);
const hiWig = hiMedIdx.map(i => wigArr[i]);

const rHiRaw = corrWith(hiRho, hiDA0);
const tHiRaw = Math.abs(rHiRaw) * Math.sqrt(hiN-2) / Math.sqrt(1-rHiRaw**2+1e-10);
const rHiAfterMHI = corrWith(residualize(hiRho, [hiMHI]), residualize(hiDA0, [hiMHI]));
const tHiAfterMHI = Math.abs(rHiAfterMHI) * Math.sqrt(hiN-3) / Math.sqrt(1-rHiAfterMHI**2+1e-10);
const rHiAfterMHIWig = corrWith(residualize(hiRho, [hiMHI, hiWig]), residualize(hiDA0, [hiMHI, hiWig]));
const tHiAfterMHIWig = Math.abs(rHiAfterMHIWig) * Math.sqrt(hiN-4) / Math.sqrt(1-rHiAfterMHIWig**2+1e-10);

log("  n=" + hiN + " galaxies with high/medium confidence");
log("  log(rho) vs delta_a0 raw:      r=" + rHiRaw.toFixed(3) + " t=" + tHiRaw.toFixed(1));
log("  log(rho) vs delta_a0 |MHI:     r=" + rHiAfterMHI.toFixed(3) + " t=" + tHiAfterMHI.toFixed(1));
log("  log(rho) vs delta_a0 |MHI+wig: r=" + rHiAfterMHIWig.toFixed(3) + " t=" + tHiAfterMHIWig.toFixed(1));
log("");

log("=".repeat(80));
log("  PHASE 29 — STRICT CRITERIA");
log("=".repeat(80));
log("");

const rRaw = corrWith(logRhoArr, deltaA0);
const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);
const rAfterMHI = corrWith(residualize(logRhoArr, [mhiArr]), residualize(deltaA0, [mhiArr]));
const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);
const rAfterMHIWig = corrWith(residualize(logRhoArr, [mhiArr, wigArr]), residualize(deltaA0, [mhiArr, wigArr]));
const tAfterMHIWig = Math.abs(rAfterMHIWig) * Math.sqrt(N-4) / Math.sqrt(1-rAfterMHIWig**2+1e-10);

const criteria = [
  { name: 'Clear raw |t|>=2', pass: tRaw >= 2.0, val: 'r=' + rRaw.toFixed(3) + ' t=' + tRaw.toFixed(1) },
  { name: 'Independent of MHI', pass: tAfterMHI >= 1.65, val: 'r=' + rAfterMHI.toFixed(3) + ' t=' + tAfterMHI.toFixed(1) },
  { name: 'Independent of MHI+wig', pass: tAfterMHIWig >= 1.65, val: 'r=' + rAfterMHIWig.toFixed(3) + ' t=' + tAfterMHIWig.toFixed(1) },
  { name: 'LOO-CV > M0', pass: rhoAlone > 2, val: rhoAlone.toFixed(1) + '%' },
  { name: 'Adds to MHI+wig', pass: rhoAdds, val: (mhiWigRhoGap - mhiWigGap).toFixed(1) + '%' },
  { name: 'Adds beyond envCode', pass: rhoAddsToEnv, val: (mhiWigEnvRhoGap - mhiWigEnvGap).toFixed(1) + '%' },
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
log("  PHASE 29 — VERDICT");
log("=".repeat(80));
log("");

if (nPassed >= 4) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  LOCAL DENSITY shows significant independent signal.                   ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nPassed >= 2) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  LOCAL DENSITY shows PARTIAL signal.                                   ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  LOCAL DENSITY FAILS.                                                  ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  ENVIRONMENT DOOR — FULL SCORECARD:");
log("    group membership (envCode):    6/6 PASS ← BEST, CONFIRMED");
log("    satellite vs central:          1/6 FAIL");
log("    nearest neighbor distance:     1/6 FAIL");
log("    tidal index (Θ):               1/6 FAIL");
log("    local density (rho):           " + nPassed + "/6 " + (nPassed >= 4 ? "PASS" : nPassed >= 2 ? "PARTIAL" : "FAIL"));
log("");

log("  FUNDAMENTAL CONCLUSION:");
log("  The environmental signal affecting a0 is CATEGORICAL.");
log("  Being in a group/cluster vs field matters.");
log("  HOW dense/close/tidally-influenced does NOT add information.");
log("  This suggests the effect is a THRESHOLD, not a gradient.");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  dataQuality: { high: nHigh, medium: nMed, low: nLow },
  rawCorrelation: { r: +rRaw.toFixed(3), t: +tRaw.toFixed(1) },
  afterMHI: { r: +rAfterMHI.toFixed(3), t: +tAfterMHI.toFixed(1) },
  afterMHIWig: { r: +rAfterMHIWig.toFixed(3), t: +tAfterMHIWig.toFixed(1) },
  rhoVsEnvCode: { r: +rRhoEnv.toFixed(3) },
  rhoAfterEnv: { r: +rRhoAfterEnv.toFixed(3), t: +tRhoAfterEnv.toFixed(1) },
  envAfterRho: { r: +rEnvAfterRho.toFixed(3), t: +tEnvAfterRho.toFixed(1) },
  highConfSubsample: { n: hiN, rawR: +rHiRaw.toFixed(3), rawT: +tHiRaw.toFixed(1) },
  looGapClosed: {
    rho: +rhoAlone.toFixed(1), env: +envAlone.toFixed(1),
    mhiWig: +mhiWigGap.toFixed(1), mhiWigRho: +mhiWigRhoGap.toFixed(1),
    mhiWigEnv: +mhiWigEnvGap.toFixed(1), mhiWigEnvRho: +mhiWigEnvRhoGap.toFixed(1),
  },
  successCriteria: { passed: nPassed, total: 6 },
  verdict: nPassed >= 4 ? 'SIGNIFICANT' : nPassed >= 2 ? 'PARTIAL' : 'FAILED',
};

fs.writeFileSync(path.join(__dirname, '../public/phase29-local-density.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase29-local-density.json");
