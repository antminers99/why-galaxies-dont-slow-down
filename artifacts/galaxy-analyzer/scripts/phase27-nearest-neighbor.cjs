#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "27.0.0";
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

const neighborData = {
  'NGC2403':     { dNN: 0.60, neighbor: 'M81', confidence: 'high', note: 'Well-studied M81 group, 3D sep from Karachentsev+2002' },
  'NGC0024':     { dNN: 0.35, neighbor: 'NGC253', confidence: 'high', note: 'Sculptor group, NGC253 dominant' },
  'NGC0891':     { dNN: 1.50, neighbor: 'NGC1023', confidence: 'medium', note: 'NGC1023 group, loose association' },
  'NGC5005':     { dNN: 0.30, neighbor: 'NGC5033', confidence: 'high', note: 'Paired with NGC5033 in CVn I' },
  'NGC5033':     { dNN: 0.30, neighbor: 'NGC5005', confidence: 'high', note: 'Paired with NGC5005 in CVn I' },
  'NGC5055':     { dNN: 0.50, neighbor: 'M51/NGC5194', confidence: 'high', note: 'M51 subgroup in CVn I' },
  'NGC5371':     { dNN: 1.00, neighbor: 'NGC5353/4', confidence: 'medium', note: 'CVn I cloud, near NGC5353/4 group' },
  'NGC5907':     { dNN: 0.50, neighbor: 'NGC5866', confidence: 'high', note: 'NGC5866/Draco group' },
  'NGC7331':     { dNN: 2.00, neighbor: 'NGC7217', confidence: 'medium', note: 'Dominant in small group, nearest massive ~2 Mpc' },

  'NGC3726':     { dNN: 0.40, neighbor: 'NGC3992', confidence: 'high', note: 'UMa cluster, NGC3992 nearby' },
  'NGC3769':     { dNN: 0.15, neighbor: 'NGC3893', confidence: 'high', note: 'UMa pair with NGC3893' },
  'NGC3893':     { dNN: 0.15, neighbor: 'NGC3769', confidence: 'high', note: 'UMa pair with NGC3769' },
  'NGC4013':     { dNN: 0.30, neighbor: 'NGC3992', confidence: 'medium', note: 'UMa cluster member' },
  'NGC4100':     { dNN: 0.35, neighbor: 'NGC4157', confidence: 'medium', note: 'UMa cluster, near NGC4157' },
  'NGC4138':     { dNN: 0.40, neighbor: 'NGC4051', confidence: 'medium', note: 'UMa cluster, near NGC4051' },
  'NGC4157':     { dNN: 0.30, neighbor: 'NGC4100', confidence: 'medium', note: 'UMa cluster, near NGC4100' },
  'NGC4217':     { dNN: 0.20, neighbor: 'NGC4226', confidence: 'medium', note: 'UMa cluster, near NGC4226' },
  'NGC4559':     { dNN: 0.80, neighbor: 'NGC4565', confidence: 'medium', note: 'Coma I cloud, NGC4565 nearby' },
  'UGC06973':    { dNN: 0.40, neighbor: 'NGC3992', confidence: 'medium', note: 'UMa periphery' },

  'NGC6503':     { dNN: 5.00, neighbor: 'NGC5907?', confidence: 'low', note: 'Local Void edge, extremely isolated' },
  'NGC2915':     { dNN: 4.00, neighbor: 'NGC3109?', confidence: 'low', note: 'Isolated BCD, southern sky' },
  'NGC3741':     { dNN: 2.50, neighbor: 'NGC4214?', confidence: 'low', note: 'Isolated dwarf irregular' },
  'NGC1705':     { dNN: 3.00, neighbor: 'NGC1559?', confidence: 'low', note: 'Isolated dwarf starburst' },
  'UGC01281':    { dNN: 2.00, neighbor: 'NGC672?', confidence: 'low', note: 'Nearby edge-on dwarf, sparse region' },
  'UGC05721':    { dNN: 2.50, neighbor: 'NGC3274?', confidence: 'low', note: 'Dwarf field' },
  'UGC08490':    { dNN: 2.50, neighbor: 'NGC4214?', confidence: 'low', note: 'Dwarf (IZw36), CVn void region' },
  'NGC1003':     { dNN: 3.00, neighbor: 'NGC891?', confidence: 'low', note: 'Perseus field, sparse' },

  'NGC2841':     { dNN: 3.00, neighbor: 'NGC2683?', confidence: 'low', note: 'Relatively isolated at D=14' },
  'NGC2903':     { dNN: 3.50, neighbor: 'NGC2683?', confidence: 'low', note: 'Leo field, isolated' },
  'NGC3198':     { dNN: 2.00, neighbor: 'NGC3319?', confidence: 'low', note: 'UMa periphery but not member' },
  'NGC3521':     { dNN: 3.50, neighbor: 'NGC3627?', confidence: 'low', note: 'Leo field, Leo Triplet ~3.5 Mpc away' },
  'NGC2683':     { dNN: 3.00, neighbor: 'NGC2841?', confidence: 'low', note: 'Lynx field' },
  'NGC0289':     { dNN: 3.50, neighbor: 'NGC253?', confidence: 'low', note: 'Sculptor region but far from NGC253' },
  'NGC1090':     { dNN: 4.00, neighbor: 'NGC1068?', confidence: 'low', note: 'Cetus field' },
  'NGC7814':     { dNN: 2.50, neighbor: 'NGC7331', confidence: 'medium', note: 'Near NGC7331 on sky, similar D' },
  'NGC6015':     { dNN: 3.00, neighbor: 'NGC5907?', confidence: 'low', note: 'Relatively isolated' },
  'NGC6674':     { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Hercules field D=51 Mpc' },

  'ESO563-G021': { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Southern field D=61 Mpc' },
  'F571-8':      { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'LSB field D=53 Mpc' },
  'IC4202':      { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Distant field D=100 Mpc' },
  'NGC0801':     { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Field D=81 Mpc' },
  'NGC2955':     { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Distant field D=98 Mpc' },
  'NGC2998':     { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Field D=68 Mpc' },
  'NGC6195':     { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Distant field D=128 Mpc' },
  'UGC00128':    { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'LSB field D=65 Mpc' },
  'UGC02885':    { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Giant LSB, very isolated D=81 Mpc' },
  'UGC02916':    { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Field D=65 Mpc' },
  'UGC02953':    { dNN: 2.50, neighbor: 'NGC2768?', confidence: 'low', note: 'Field D=16 Mpc, near NGC2768 group?' },
  'UGC03205':    { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Field D=50 Mpc' },
  'UGC03546':    { dNN: 4.00, neighbor: 'unknown', confidence: 'low', note: 'Field D=29 Mpc' },
  'UGC03580':    { dNN: 3.50, neighbor: 'NGC2768?', confidence: 'low', note: 'Field D=21 Mpc' },
  'UGC06786':    { dNN: 4.00, neighbor: 'unknown', confidence: 'low', note: 'Field D=29 Mpc' },
  'UGC06787':    { dNN: 3.50, neighbor: 'unknown', confidence: 'low', note: 'Field D=21 Mpc' },
  'UGC08699':    { dNN: 4.00, neighbor: 'unknown', confidence: 'low', note: 'Field D=39 Mpc' },
  'UGC09037':    { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Field D=84 Mpc' },
  'UGC09133':    { dNN: 5.00, neighbor: 'unknown', confidence: 'low', note: 'Field D=57 Mpc' },
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

  const nd = neighborData[gal.name];
  const dNN = nd ? nd.dNN : 5.0;
  const logDNN = Math.log10(dNN);
  const confidence = nd ? nd.confidence : 'unknown';

  const ev = envLookup[gal.name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    n, logMHI, logVflat, rcWiggliness,
    dNN, logDNN, confidence, envCode,
    neighbor: nd ? nd.neighbor : 'unknown',
    note: nd ? nd.note : '',
    D: sp.D,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const nHigh = galaxyData.filter(g => g.confidence === 'high').length;
const nMed = galaxyData.filter(g => g.confidence === 'medium').length;
const nLow = galaxyData.filter(g => g.confidence === 'low').length;

log("");
log("=".repeat(80));
log("  PHASE 27: NEAREST MASSIVE NEIGHBOR DISTANCE");
log("  Does proximity to a massive galaxy predict a0 variation?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  ENVIRONMENT DOOR SCORECARD:");
sep();
log("  group membership (envCode): 6/6 PASS ← best environmental variable");
log("  satellite vs central:       1/6 FAIL");
log("  nearest neighbor distance:  ← TESTING NOW");
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
log("  CAVEAT: Without RA/Dec in SPARC, dNN for field galaxies relies on");
log("  astronomical knowledge and literature. Many distant field galaxies");
log("  receive a default dNN=5.0 Mpc (effective ceiling).");
log("");

log("  PER-GALAXY:");
const sorted = [...galaxyData].sort((a,b) => a.dNN - b.dNN);
for (const g of sorted) {
  log("    " + g.name.padEnd(16) + "dNN=" + g.dNN.toFixed(2).padStart(5) + " Mpc  " +
    g.confidence.padEnd(7) + " " + (g.neighbor||'').padEnd(14) + " logA0=" + g.logA0.toFixed(3));
}
log("");

log("=".repeat(80));
log("  TEST 1: CORRELATIONS — dNN vs delta_a0");
log("=".repeat(80));
log("");

const dNNArr = galaxyData.map(g => g.dNN);
const logDNNArr = galaxyData.map(g => g.logDNN);
const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const vflatArr = galaxyData.map(g => g.logVflat);
const nArr = galaxyData.map(g => g.n);

for (const [label, arr] of [['dNN', dNNArr], ['log(dNN)', logDNNArr]]) {
  const r = corrWith(arr, deltaA0);
  const t = Math.abs(r) * Math.sqrt(N-2) / Math.sqrt(1-r**2+1e-10);
  log("  " + label.padEnd(12) + " vs delta_a0: r=" + r.toFixed(3) + " t=" + t.toFixed(1) +
    " → " + (r > 0 ? "farther→higher a0" : "closer→higher a0"));
}
log("");

const rDNNEnv = corrWith(logDNNArr, envArr);
log("  log(dNN) vs envCode: r=" + rDNNEnv.toFixed(3) + " (collinearity check)");
log("  " + (Math.abs(rDNNEnv) > 0.6 ? "HIGH collinearity — dNN largely redundant with envCode" : "Moderate collinearity — some independent info"));
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
  { name: '|MHI+Vflat+wig+env', vals: [mhiArr, vflatArr, wigArr, envArr] },
];

for (const [label, arr] of [['log(dNN)', logDNNArr], ['dNN', dNNArr]]) {
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
log("  TEST 3: dNN AFTER envCode — DOES IT ADD?");
log("=".repeat(80));
log("");

const residAfterEnv = residualize(deltaA0, [mhiArr, wigArr, envArr]);
const rDNNafterEnv = corrWith(residualize(logDNNArr, [mhiArr, wigArr, envArr]), residAfterEnv);
const tDNNafterEnv = Math.abs(rDNNafterEnv) * Math.sqrt(N-5) / Math.sqrt(1-rDNNafterEnv**2+1e-10);

const residAfterDNN = residualize(deltaA0, [mhiArr, wigArr, logDNNArr]);
const rEnvAfterDNN = corrWith(residualize(envArr, [mhiArr, wigArr, logDNNArr]), residAfterDNN);
const tEnvAfterDNN = Math.abs(rEnvAfterDNN) * Math.sqrt(N-5) / Math.sqrt(1-rEnvAfterDNN**2+1e-10);

log("  log(dNN) AFTER envCode+MHI+wig: r=" + rDNNafterEnv.toFixed(3) + " t=" + tDNNafterEnv.toFixed(1) +
  " → " + (tDNNafterEnv >= 1.65 ? "ADDS beyond envCode" : "ABSORBED by envCode"));
log("  envCode AFTER log(dNN)+MHI+wig: r=" + rEnvAfterDNN.toFixed(3) + " t=" + tEnvAfterDNN.toFixed(1) +
  " → " + (tEnvAfterDNN >= 1.65 ? "envCode STILL adds" : "envCode ABSORBED"));
log("");

if (tDNNafterEnv < 1.65 && tEnvAfterDNN >= 1.65) {
  log("  → envCode ABSORBS dNN. Group membership is the deeper variable.");
} else if (tDNNafterEnv >= 1.65 && tEnvAfterDNN < 1.65) {
  log("  → dNN ABSORBS envCode. Physical proximity is the deeper variable.");
} else if (tDNNafterEnv >= 1.65 && tEnvAfterDNN >= 1.65) {
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
  ['log(dNN)|MHI', logDNNArr, [mhiArr]],
  ['log(dNN)|MHI+wig', logDNNArr, [mhiArr, wigArr]],
  ['log(dNN)|MHI+wig+env', logDNNArr, [mhiArr, wigArr, envArr]],
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
  log("  " + label.padEnd(24) + " |r|=" + obsR.toFixed(3) + " perm_p=" + (cnt/NPERMS).toFixed(4) +
    (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "*" : ""));
}
log("");

log("=".repeat(80));
log("  TEST 5: LOO CROSS-VALIDATION");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'logDNN': logDNNArr,
  'envCode': envArr, 'logVflat': vflatArr,
};

function getVal(gIdx, fname) { return allFeatures[fname][gIdx]; }

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_dNN: log(dNN)', features: ['logDNN'] },
  { name: 'M_env: envCode', features: ['envCode'] },
  { name: 'M_MHI', features: ['logMHI'] },
  { name: 'M_MHI+dNN', features: ['logMHI', 'logDNN'] },
  { name: 'M_MHI+env', features: ['logMHI', 'envCode'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
  { name: 'M_MHI+wig+dNN', features: ['logMHI', 'rcWiggliness', 'logDNN'] },
  { name: 'M_MHI+wig+env', features: ['logMHI', 'rcWiggliness', 'envCode'] },
  { name: 'M_MHI+wig+env+dNN', features: ['logMHI', 'rcWiggliness', 'envCode', 'logDNN'] },
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
log("  │  Model                   k    CV-RMS   vs M0   gap-closed              │");
log("  ├───────────────────────────────────────────────────────────────────────────┤");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const gc = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  │  " + r.name.padEnd(23) + r.k.toString().padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(8) + (gc + "%").padStart(13) + "             │");
}
log("  └───────────────────────────────────────────────────────────────────────────┘");
log("");

const mhiWigGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig').cvRMS) / gap * 100 : 0;
const mhiWigDNNGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig+dNN').cvRMS) / gap * 100 : 0;
const mhiWigEnvGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig+env').cvRMS) / gap * 100 : 0;
const mhiWigEnvDNNGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig+env+dNN').cvRMS) / gap * 100 : 0;
const dNNalone = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_dNN: log(dNN)').cvRMS) / gap * 100 : 0;
const envAlone = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_env: envCode').cvRMS) / gap * 100 : 0;

log("  KEY COMPARISONS:");
log("    dNN alone:                " + dNNalone.toFixed(1) + "%");
log("    envCode alone:            " + envAlone.toFixed(1) + "%");
log("    MHI + wig:                " + mhiWigGap.toFixed(1) + "%");
log("    MHI + wig + dNN:          " + mhiWigDNNGap.toFixed(1) + "%");
log("    MHI + wig + env:          " + mhiWigEnvGap.toFixed(1) + "%");
log("    MHI + wig + env + dNN:    " + mhiWigEnvDNNGap.toFixed(1) + "%");
log("");

const dNNadds = mhiWigDNNGap > mhiWigGap + 1;
const dNNaddsToEnv = mhiWigEnvDNNGap > mhiWigEnvGap + 1;
log("  dNN adds beyond MHI+wig? " + (dNNadds ? "YES (+" + (mhiWigDNNGap-mhiWigGap).toFixed(1) + "%)" : "NO"));
log("  dNN adds beyond MHI+wig+env? " + (dNNaddsToEnv ? "YES (+" + (mhiWigEnvDNNGap-mhiWigEnvGap).toFixed(1) + "%)" : "NO"));
log("");

log("=".repeat(80));
log("  TEST 6: HIGH-CONFIDENCE SUBSAMPLE");
log("  (Only galaxies with high/medium confidence dNN estimates)");
log("=".repeat(80));
log("");

const hiMedIdx = galaxyData.map((g,i) => (g.confidence === 'high' || g.confidence === 'medium') ? i : -1).filter(i => i >= 0);
const hiN = hiMedIdx.length;
const hiDeltaA0 = hiMedIdx.map(i => deltaA0[i]);
const hiLogDNN = hiMedIdx.map(i => logDNNArr[i]);
const hiMHI = hiMedIdx.map(i => mhiArr[i]);
const hiWig = hiMedIdx.map(i => wigArr[i]);
const hiEnv = hiMedIdx.map(i => envArr[i]);

const rHiRaw = corrWith(hiLogDNN, hiDeltaA0);
const tHiRaw = Math.abs(rHiRaw) * Math.sqrt(hiN-2) / Math.sqrt(1-rHiRaw**2+1e-10);

const rHiAfterMHI = corrWith(residualize(hiLogDNN, [hiMHI]), residualize(hiDeltaA0, [hiMHI]));
const tHiAfterMHI = Math.abs(rHiAfterMHI) * Math.sqrt(hiN-3) / Math.sqrt(1-rHiAfterMHI**2+1e-10);

log("  n=" + hiN + " galaxies with high/medium confidence");
log("  log(dNN) vs delta_a0 raw: r=" + rHiRaw.toFixed(3) + " t=" + tHiRaw.toFixed(1));
log("  log(dNN) vs delta_a0 |MHI: r=" + rHiAfterMHI.toFixed(3) + " t=" + tHiAfterMHI.toFixed(1));
log("");

log("=".repeat(80));
log("  PHASE 27 — STRICT CRITERIA");
log("=".repeat(80));
log("");

const rRaw = corrWith(logDNNArr, deltaA0);
const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);

const rAfterMHI = corrWith(residualize(logDNNArr, [mhiArr]), residualize(deltaA0, [mhiArr]));
const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);

const rAfterMHIWig = corrWith(residualize(logDNNArr, [mhiArr, wigArr]), residualize(deltaA0, [mhiArr, wigArr]));
const tAfterMHIWig = Math.abs(rAfterMHIWig) * Math.sqrt(N-4) / Math.sqrt(1-rAfterMHIWig**2+1e-10);

const criteria = [
  { name: 'Clear raw |t|>=2', pass: tRaw >= 2.0, val: 'r=' + rRaw.toFixed(3) + ' t=' + tRaw.toFixed(1) },
  { name: 'Independent of MHI', pass: tAfterMHI >= 1.65, val: 'r=' + rAfterMHI.toFixed(3) + ' t=' + tAfterMHI.toFixed(1) },
  { name: 'Independent of MHI+wig', pass: tAfterMHIWig >= 1.65, val: 'r=' + rAfterMHIWig.toFixed(3) + ' t=' + tAfterMHIWig.toFixed(1) },
  { name: 'LOO-CV > M0', pass: dNNalone > 2, val: dNNalone.toFixed(1) + '%' },
  { name: 'Adds to MHI+wig', pass: dNNadds, val: (mhiWigDNNGap - mhiWigGap).toFixed(1) + '%' },
  { name: 'Adds beyond envCode', pass: dNNaddsToEnv, val: (mhiWigEnvDNNGap - mhiWigEnvGap).toFixed(1) + '%' },
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
log("  PHASE 27 — VERDICT");
log("=".repeat(80));
log("");

if (nPassed >= 4) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  NEAREST NEIGHBOR DISTANCE shows significant independent signal.       ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nPassed >= 2) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  NEAREST NEIGHBOR DISTANCE shows PARTIAL signal but does not clearly   ║");
  log("  ║  improve on envCode (group membership).                                ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  NEAREST NEIGHBOR DISTANCE FAILS.                                      ║");
  log("  ║  Physical proximity does not predict a0 variation independently.       ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  ENVIRONMENT DOOR — UPDATED SCORECARD:");
log("    group membership (envCode):    6/6 PASS ← best");
log("    satellite vs central:          1/6 FAIL");
log("    nearest neighbor distance:     " + nPassed + "/6 " + (nPassed >= 4 ? "PASS" : nPassed >= 2 ? "PARTIAL" : "FAIL"));
log("");

log("  DATA QUALITY NOTE:");
log("  " + nLow + "/" + N + " galaxies have low-confidence dNN estimates");
log("  (mostly distant field galaxies assigned default 5.0 Mpc ceiling).");
log("  This truncation weakens any real signal and inflates collinearity");
log("  with envCode. A definitive test requires proper 3D positions.");
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
  dNNvsEnvCode: { r: +rDNNEnv.toFixed(3) },
  dNNafterEnv: { r: +rDNNafterEnv.toFixed(3), t: +tDNNafterEnv.toFixed(1) },
  envAfterDNN: { r: +rEnvAfterDNN.toFixed(3), t: +tEnvAfterDNN.toFixed(1) },
  looGapClosed: {
    dNN: +dNNalone.toFixed(1), env: +envAlone.toFixed(1),
    mhiWig: +mhiWigGap.toFixed(1), mhiWigDNN: +mhiWigDNNGap.toFixed(1),
    mhiWigEnv: +mhiWigEnvGap.toFixed(1), mhiWigEnvDNN: +mhiWigEnvDNNGap.toFixed(1),
  },
  successCriteria: { passed: nPassed, total: 6 },
  verdict: nPassed >= 4 ? 'SIGNIFICANT' : nPassed >= 2 ? 'PARTIAL' : 'FAILED',
  galaxyNeighbors: galaxyData.map(g => ({
    name: g.name, dNN: g.dNN, logDNN: +g.logDNN.toFixed(3),
    neighbor: g.neighbor, confidence: g.confidence
  }))
};

fs.writeFileSync(path.join(__dirname, '../public/phase27-nearest-neighbor.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase27-nearest-neighbor.json");
