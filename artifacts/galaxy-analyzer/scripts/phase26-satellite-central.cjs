#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "26.0.0";
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

function dlMeta(logA0s, ses) {
  const n = logA0s.length;
  const w = ses.map(s => 1 / (s * s + 1e-10));
  const wS = w.reduce((a, b) => a + b, 0);
  const muFE = logA0s.reduce((s, v, i) => s + w[i] * v, 0) / wS;
  let Q = 0;
  for (let i = 0; i < n; i++) Q += w[i] * (logA0s[i] - muFE) ** 2;
  const S1 = wS, S2 = w.reduce((s, v) => s + v * v, 0);
  const tau2 = Math.max(0, (Q - (n - 1)) / (S1 - S2 / S1));
  const tau = Math.sqrt(tau2);
  const W = ses.map(s => 1 / (s * s + tau2 + 1e-10));
  const WS = W.reduce((a, b) => a + b, 0);
  const mu = logA0s.reduce((s, v, i) => s + W[i] * v, 0) / WS;
  return { mu, tau, tau2, a0: Math.pow(10, mu) };
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

const satCentralMap = {
  'NGC2403': { role: 'satellite', host: 'M81/M82', note: 'Satellite of M81 group, D~3.2 Mpc' },
  'NGC0024': { role: 'satellite', host: 'NGC253', note: 'Sculptor group, satellite of NGC253' },
  'NGC0891': { role: 'central', host: 'self', note: 'Brightest in NGC1023 group sub-condensation' },
  'NGC5005': { role: 'central', host: 'self', note: 'Dominant in NGC5005/5033 pair, CVn I' },
  'NGC5033': { role: 'satellite', host: 'NGC5005', note: 'Secondary in NGC5005/5033 pair' },
  'NGC5055': { role: 'central', host: 'self', note: 'Central in M51 subgroup, CVn I cloud' },
  'NGC5371': { role: 'central', host: 'self', note: 'Bright member of CVn I cloud' },
  'NGC5907': { role: 'satellite', host: 'NGC5866', note: 'NGC5866/Draco group, not the brightest' },
  'NGC7331': { role: 'central', host: 'self', note: 'Central of NGC7331 (Deer Lick) group' },

  'NGC3726': { role: 'satellite', host: 'UMa-center', note: 'Ursa Major cluster member' },
  'NGC3769': { role: 'satellite', host: 'UMa-center', note: 'Ursa Major cluster, near NGC3893' },
  'NGC3893': { role: 'satellite', host: 'UMa-center', note: 'Ursa Major cluster member' },
  'NGC4013': { role: 'satellite', host: 'UMa-center', note: 'Ursa Major cluster member' },
  'NGC4100': { role: 'satellite', host: 'UMa-center', note: 'Ursa Major cluster member' },
  'NGC4138': { role: 'satellite', host: 'UMa-center', note: 'Ursa Major / CVn border' },
  'NGC4157': { role: 'satellite', host: 'UMa-center', note: 'Ursa Major cluster member' },
  'NGC4217': { role: 'satellite', host: 'UMa-center', note: 'Ursa Major cluster, near NGC4226' },
  'NGC4559': { role: 'satellite', host: 'ComaI', note: 'Coma I cloud member' },
  'UGC06973': { role: 'satellite', host: 'UMa-center', note: 'Ursa Major cluster periphery' },

  'ESO563-G021': { role: 'isolated', host: 'none', note: 'Field galaxy D=61 Mpc' },
  'F571-8': { role: 'isolated', host: 'none', note: 'LSB field galaxy D=53 Mpc' },
  'IC4202': { role: 'isolated', host: 'none', note: 'Field galaxy D=100 Mpc' },
  'NGC0289': { role: 'isolated', host: 'none', note: 'Field in Sculptor' },
  'NGC0801': { role: 'isolated', host: 'none', note: 'Field galaxy D=81 Mpc' },
  'NGC1003': { role: 'isolated', host: 'none', note: 'Field in Perseus' },
  'NGC1090': { role: 'isolated', host: 'none', note: 'Field in Cetus' },
  'NGC1705': { role: 'isolated', host: 'none', note: 'Isolated dwarf starburst' },
  'NGC2683': { role: 'isolated', host: 'none', note: 'Field in Lynx' },
  'NGC2841': { role: 'isolated', host: 'none', note: 'Relatively isolated D=14 Mpc' },
  'NGC2903': { role: 'isolated', host: 'none', note: 'Leo field, isolated' },
  'NGC2915': { role: 'isolated', host: 'none', note: 'Isolated blue compact dwarf' },
  'NGC2955': { role: 'isolated', host: 'none', note: 'Distant field D=98 Mpc' },
  'NGC2998': { role: 'isolated', host: 'none', note: 'Field D=68 Mpc' },
  'NGC3198': { role: 'isolated', host: 'none', note: 'UMa periphery, not member' },
  'NGC3521': { role: 'isolated', host: 'none', note: 'Leo field, isolated' },
  'NGC3741': { role: 'isolated', host: 'none', note: 'Isolated dwarf irregular' },
  'NGC6015': { role: 'isolated', host: 'none', note: 'Relatively isolated' },
  'NGC6195': { role: 'isolated', host: 'none', note: 'Distant field D=128 Mpc' },
  'NGC6503': { role: 'isolated', host: 'none', note: 'Local Void edge, very isolated' },
  'NGC6674': { role: 'isolated', host: 'none', note: 'Field in Hercules' },
  'NGC7814': { role: 'isolated', host: 'none', note: 'Relatively isolated edge-on' },
  'UGC00128': { role: 'isolated', host: 'none', note: 'LSB field' },
  'UGC01281': { role: 'isolated', host: 'none', note: 'Nearby edge-on dwarf' },
  'UGC02885': { role: 'isolated', host: 'none', note: 'Giant LSB, very isolated' },
  'UGC02916': { role: 'isolated', host: 'none', note: 'Field' },
  'UGC02953': { role: 'isolated', host: 'none', note: 'Field' },
  'UGC03205': { role: 'isolated', host: 'none', note: 'Field' },
  'UGC03546': { role: 'isolated', host: 'none', note: 'Field' },
  'UGC03580': { role: 'isolated', host: 'none', note: 'Field' },
  'UGC05721': { role: 'isolated', host: 'none', note: 'Dwarf field' },
  'UGC06786': { role: 'isolated', host: 'none', note: 'Field' },
  'UGC06787': { role: 'isolated', host: 'none', note: 'Field' },
  'UGC08490': { role: 'isolated', host: 'none', note: 'Dwarf field (IZw36)' },
  'UGC08699': { role: 'isolated', host: 'none', note: 'Field' },
  'UGC09037': { role: 'isolated', host: 'none', note: 'Field' },
  'UGC09133': { role: 'isolated', host: 'none', note: 'Field' },
};

const p11 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase11-sensitivity-lab.json'), 'utf8'));
const sparcAll = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));

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

  const sc = satCentralMap[gal.name];
  const role = sc ? sc.role : 'unknown';
  const isSatellite = role === 'satellite' ? 1 : 0;
  const isCentral = role === 'central' ? 1 : 0;
  const isIsolated = role === 'isolated' ? 1 : 0;
  const roleCode = role === 'isolated' ? 0 : role === 'central' ? 1 : role === 'satellite' ? 2 : -1;

  const envCode = (sc && (role === 'satellite' || role === 'central')) ? 1 : 0;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q,
    n, logMHI, logVflat, rcWiggliness,
    role, isSatellite, isCentral, isIsolated, roleCode, envCode,
    host: sc ? sc.host : 'unknown', note: sc ? sc.note : '',
    Vflat: sp.Vflat,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const nIso = galaxyData.filter(g => g.role === 'isolated').length;
const nCen = galaxyData.filter(g => g.role === 'central').length;
const nSat = galaxyData.filter(g => g.role === 'satellite').length;
const nUnk = galaxyData.filter(g => g.role === 'unknown').length;

log("");
log("=".repeat(80));
log("  PHASE 26: ENVIRONMENT — SATELLITE vs CENTRAL");
log("  Is the environmental signal from satellite status specifically?");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  CLASSIFICATION:");
sep();
log("");
log("  Criteria:");
log("    ISOLATED:  field galaxy, no known group membership");
log("    CENTRAL:   brightest/dominant galaxy in its group");
log("    SATELLITE: non-dominant member of group or cluster");
log("");
log("  ┌────────────────────────────────────────────────────────┐");
log("  │  Role          N     fraction                        │");
log("  ├────────────────────────────────────────────────────────┤");
log("  │  isolated      " + nIso.toString().padStart(2) + "    " + (nIso/N*100).toFixed(0) + "%".padEnd(32) + "│");
log("  │  central        " + nCen.toString().padStart(2) + "    " + (nCen/N*100).toFixed(0) + "%".padEnd(32) + "│");
log("  │  satellite     " + nSat.toString().padStart(2) + "    " + (nSat/N*100).toFixed(0) + "%".padEnd(32) + "│");
if (nUnk > 0) log("  │  unknown        " + nUnk.toString().padStart(2) + "    " + (nUnk/N*100).toFixed(0) + "%".padEnd(32) + "│");
log("  └────────────────────────────────────────────────────────┘");
log("");

log("  PER-GALAXY:");
for (const g of galaxyData) {
  log("    " + g.name.padEnd(16) + g.role.padEnd(12) + (g.host||'').padEnd(14) + "D=" + g.D.toFixed(1));
}
log("");

log("=".repeat(80));
log("  TEST 1: THREE-WAY COMPARISON — isolated / central / satellite");
log("=".repeat(80));
log("");

const roles = ['isolated', 'central', 'satellite'];
const fullDL = dlMeta(allLogA0, galaxyData.map(g => g.se));

for (const r of roles) {
  const idx = galaxyData.map((g,i) => g.role === r ? i : -1).filter(i => i >= 0);
  if (idx.length < 2) { log("  " + r + ": n=" + idx.length + " (too few)"); continue; }
  const la = idx.map(i => allLogA0[i]);
  const m = la.reduce((s,v)=>s+v,0)/la.length;
  const sd = Math.sqrt(la.reduce((s,v)=>s+(v-m)**2,0)/Math.max(la.length-1,1));
  const dl = dlMeta(la, idx.map(i => galaxyData[i].se));
  const wigM = idx.map(i => galaxyData[i].rcWiggliness).reduce((s,v)=>s+v,0)/idx.length;
  const mhiM = idx.map(i => galaxyData[i].logMHI).reduce((s,v)=>s+v,0)/idx.length;
  log("  " + r.padEnd(12) + " n=" + idx.length.toString().padStart(2) +
    "  a0=" + Math.round(Math.pow(10, m)).toString().padStart(5) +
    "  logA0=" + m.toFixed(3) + "  tau=" + dl.tau.toFixed(3) +
    "  meanWig=" + wigM.toFixed(4) + "  meanMHI=" + mhiM.toFixed(2));
}
log("");

const isoIdx = galaxyData.map((g,i) => g.role === 'isolated' ? i : -1).filter(i => i >= 0);
const satIdx = galaxyData.map((g,i) => g.role === 'satellite' ? i : -1).filter(i => i >= 0);
const cenIdx = galaxyData.map((g,i) => g.role === 'central' ? i : -1).filter(i => i >= 0);

const isoA0 = isoIdx.map(i => allLogA0[i]);
const satA0 = satIdx.map(i => allLogA0[i]);
const cenA0 = cenIdx.length >= 2 ? cenIdx.map(i => allLogA0[i]) : [];

const isoMean = isoA0.reduce((s,v)=>s+v,0)/isoA0.length;
const satMean = satA0.reduce((s,v)=>s+v,0)/satA0.length;
const isoSD = Math.sqrt(isoA0.reduce((s,v)=>s+(v-isoMean)**2,0)/(isoA0.length-1));
const satSD = Math.sqrt(satA0.reduce((s,v)=>s+(v-satMean)**2,0)/(satA0.length-1));
const isoSatSplit = Math.abs(isoMean - satMean);
const isoSatT = isoSatSplit / Math.sqrt(isoSD**2/isoIdx.length + satSD**2/satIdx.length);

log("  ISOLATED vs SATELLITE:");
log("    Split: " + isoSatSplit.toFixed(3) + " dex (t=" + isoSatT.toFixed(2) + ")");
log("    Direction: " + (isoMean > satMean ? "isolated > satellite" : "satellite > isolated"));
log("");

if (cenIdx.length >= 2) {
  const cenMean = cenA0.reduce((s,v)=>s+v,0)/cenA0.length;
  const cenSD = Math.sqrt(cenA0.reduce((s,v)=>s+(v-cenMean)**2,0)/Math.max(cenA0.length-1,1));
  const cenSatSplit = Math.abs(cenMean - satMean);
  const cenIsoSplit = Math.abs(cenMean - isoMean);
  log("  CENTRAL vs SATELLITE: " + cenSatSplit.toFixed(3) + " dex");
  log("  CENTRAL vs ISOLATED:  " + cenIsoSplit.toFixed(3) + " dex");
  log("  Central closer to: " + (Math.abs(cenMean - isoMean) < Math.abs(cenMean - satMean) ? "isolated" : "satellite"));
  log("");
}

log("=".repeat(80));
log("  TEST 2: CORRELATIONS AND CONFOUNDER STRIPPING");
log("=".repeat(80));
log("");

const satArr = galaxyData.map(g => g.isSatellite);
const roleArr = galaxyData.map(g => g.roleCode);
const envArr = galaxyData.map(g => g.envCode);
const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const vflatArr = galaxyData.map(g => g.logVflat);
const nArr = galaxyData.map(g => g.n);

const testVars = [
  { name: 'isSatellite', values: satArr },
  { name: 'roleCode', values: roleArr },
  { name: 'envCode(Phase25)', values: envArr },
];

const controlSets = [
  { name: 'raw', vals: [] },
  { name: '|MHI', vals: [mhiArr] },
  { name: '|MHI+Vflat', vals: [mhiArr, vflatArr] },
  { name: '|MHI+wig', vals: [mhiArr, wigArr] },
  { name: '|MHI+Vflat+wig', vals: [mhiArr, vflatArr, wigArr] },
  { name: '|MHI+Vflat+wig+n', vals: [mhiArr, vflatArr, wigArr, nArr] },
];

for (const tv of testVars) {
  log("  " + tv.name + " vs delta_a0:");
  log("  ┌──────────────────────────────────────────────────────────────┐");
  log("  │  Controls                 r_partial   |t|    status       │");
  log("  ├──────────────────────────────────────────────────────────────┤");
  for (const cs of controlSets) {
    const residY = residualize(deltaA0, cs.vals);
    const residX = residualize(tv.values, cs.vals);
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
log("  TEST 3: PERMUTATION TESTS (after MHI and after MHI+wig)");
log("=".repeat(80));
log("");

const NPERMS = 10000;

for (const tv of testVars) {
  for (const [label, controlVals] of [['|MHI', [mhiArr]], ['|MHI+wig', [mhiArr, wigArr]]]) {
    const residY = residualize(deltaA0, controlVals);
    const residX = residualize(tv.values, controlVals);
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
    log("  " + tv.name.padEnd(18) + label.padEnd(12) + " |r|=" + obsR.toFixed(3) +
      " perm_p=" + (cnt/NPERMS).toFixed(4) + (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "*" : ""));
  }
  log("");
}

log("=".repeat(80));
log("  TEST 4: LOO CROSS-VALIDATION");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr,
  'rcWiggliness': wigArr,
  'isSatellite': satArr,
  'roleCode': roleArr,
  'envCode': envArr,
  'logVflat': vflatArr,
};

function getVal(gIdx, fname) { return allFeatures[fname][gIdx]; }

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_sat: isSatellite', features: ['isSatellite'] },
  { name: 'M_role: roleCode', features: ['roleCode'] },
  { name: 'M_env: envCode(P25)', features: ['envCode'] },
  { name: 'M_MHI', features: ['logMHI'] },
  { name: 'M_MHI+sat', features: ['logMHI', 'isSatellite'] },
  { name: 'M_MHI+role', features: ['logMHI', 'roleCode'] },
  { name: 'M_MHI+env', features: ['logMHI', 'envCode'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
  { name: 'M_MHI+wig+sat', features: ['logMHI', 'rcWiggliness', 'isSatellite'] },
  { name: 'M_MHI+wig+role', features: ['logMHI', 'rcWiggliness', 'roleCode'] },
  { name: 'M_MHI+wig+env', features: ['logMHI', 'rcWiggliness', 'envCode'] },
  { name: 'M_full4', features: ['logMHI', 'rcWiggliness', 'envCode', 'isSatellite'] },
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
log("  │  Model                  k    CV-RMS    vs M0   gap-closed              │");
log("  ├───────────────────────────────────────────────────────────────────────────┤");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const gc = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  │  " + r.name.padEnd(22) + r.k.toString().padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(9) + (gc + "%").padStart(13) + "             │");
}
log("  └───────────────────────────────────────────────────────────────────────────┘");
log("");

const mhiGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI').cvRMS) / gap * 100 : 0;
const mhiWigGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig').cvRMS) / gap * 100 : 0;
const mhiSatGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+sat').cvRMS) / gap * 100 : 0;
const mhiWigSatGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig+sat').cvRMS) / gap * 100 : 0;
const mhiWigEnvGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig+env').cvRMS) / gap * 100 : 0;
const mhiWigRoleGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig+role').cvRMS) / gap * 100 : 0;

log("  KEY COMPARISONS:");
log("    MHI alone:              " + mhiGap.toFixed(1) + "%");
log("    MHI + wig:              " + mhiWigGap.toFixed(1) + "%");
log("    MHI + sat:              " + mhiSatGap.toFixed(1) + "%");
log("    MHI + wig + sat:        " + mhiWigSatGap.toFixed(1) + "%");
log("    MHI + wig + env(P25):   " + mhiWigEnvGap.toFixed(1) + "%");
log("    MHI + wig + role:       " + mhiWigRoleGap.toFixed(1) + "%");
log("");

log("  Does satellite/central add beyond MHI+wig?");
const satAdds = mhiWigSatGap > mhiWigGap + 1;
const roleAdds = mhiWigRoleGap > mhiWigGap + 1;
const envAdds = mhiWigEnvGap > mhiWigGap + 1;
log("    isSatellite: " + (satAdds ? "YES (+" + (mhiWigSatGap - mhiWigGap).toFixed(1) + "%)" : "NO"));
log("    roleCode:    " + (roleAdds ? "YES (+" + (mhiWigRoleGap - mhiWigGap).toFixed(1) + "%)" : "NO"));
log("    envCode(P25): " + (envAdds ? "YES (+" + (mhiWigEnvGap - mhiWigGap).toFixed(1) + "%)" : "NO"));
log("");

log("=".repeat(80));
log("  TEST 5: DOES satellite/central ABSORB envCode OR VICE VERSA?");
log("=".repeat(80));
log("");

const residAfterEnv = residualize(deltaA0, [mhiArr, wigArr, envArr]);
const rSatAfterEnv = corrWith(residualize(satArr, [mhiArr, wigArr, envArr]), residAfterEnv);
const tSatAfterEnv = Math.abs(rSatAfterEnv) * Math.sqrt(N-5) / Math.sqrt(1-rSatAfterEnv**2+1e-10);

const residAfterSat = residualize(deltaA0, [mhiArr, wigArr, satArr]);
const rEnvAfterSat = corrWith(residualize(envArr, [mhiArr, wigArr, satArr]), residAfterSat);
const tEnvAfterSat = Math.abs(rEnvAfterSat) * Math.sqrt(N-5) / Math.sqrt(1-rEnvAfterSat**2+1e-10);

log("  isSatellite AFTER envCode+MHI+wig: r=" + rSatAfterEnv.toFixed(3) + " t=" + tSatAfterEnv.toFixed(1) +
  " → " + (tSatAfterEnv >= 1.65 ? "INDEPENDENT" : "ABSORBED"));
log("  envCode AFTER isSatellite+MHI+wig: r=" + rEnvAfterSat.toFixed(3) + " t=" + tEnvAfterSat.toFixed(1) +
  " → " + (tEnvAfterSat >= 1.65 ? "INDEPENDENT" : "ABSORBED"));
log("");

if (tSatAfterEnv < 1.65 && tEnvAfterSat >= 1.65) {
  log("  → envCode ABSORBS isSatellite. Group membership matters more than sat/central distinction.");
} else if (tSatAfterEnv >= 1.65 && tEnvAfterSat < 1.65) {
  log("  → isSatellite ABSORBS envCode. Satellite status is the deeper variable.");
} else if (tSatAfterEnv >= 1.65 && tEnvAfterSat >= 1.65) {
  log("  → BOTH carry independent information. Environment has multiple channels.");
} else {
  log("  → BOTH absorbed when combined. They measure roughly the same thing.");
}
log("");

log("=".repeat(80));
log("  TEST 6: DOES ENVIRONMENT EXPLAIN rcWiggliness?");
log("=".repeat(80));
log("");

const rSatWig = corrWith(satArr, wigArr);
log("  isSatellite vs rcWiggliness: r=" + rSatWig.toFixed(3));
log("  " + (rSatWig < 0 ? "Satellites LESS wiggly" : "Satellites MORE wiggly"));
log("");

const isoWigMean = isoIdx.map(i => wigArr[i]).reduce((s,v)=>s+v,0)/isoIdx.length;
const satWigMean = satIdx.map(i => wigArr[i]).reduce((s,v)=>s+v,0)/satIdx.length;
const cenWigMean = cenIdx.length > 0 ? cenIdx.map(i => wigArr[i]).reduce((s,v)=>s+v,0)/cenIdx.length : 0;

log("  Mean rcWiggliness by role:");
log("    isolated:    " + isoWigMean.toFixed(4));
if (cenIdx.length > 0) log("    central:     " + cenWigMean.toFixed(4));
log("    satellite:   " + satWigMean.toFixed(4));
log("");

const residWigAfterMHI = residualize(deltaA0, [mhiArr, satArr]);
const rWigAfterSat = corrWith(wigArr, residWigAfterMHI);
const tWigAfterSat = Math.abs(rWigAfterSat) * Math.sqrt(N-4) / Math.sqrt(1-rWigAfterSat**2+1e-10);
log("  rcWiggliness AFTER MHI+isSatellite: r=" + rWigAfterSat.toFixed(3) + " t=" + tWigAfterSat.toFixed(1) +
  " → " + (tWigAfterSat >= 1.65 ? "STILL INDEPENDENT" : "ABSORBED"));
log("");

log("=".repeat(80));
log("  PHASE 26 — STRICT CRITERIA");
log("=".repeat(80));
log("");

const residYmhi = residualize(deltaA0, [mhiArr]);
const rSatAfterMHI = corrWith(residualize(satArr, [mhiArr]), residYmhi);
const tSatAfterMHI = Math.abs(rSatAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rSatAfterMHI**2+1e-10);

const residYmhiWig = residualize(deltaA0, [mhiArr, wigArr]);
const rSatAfterMHIWig = corrWith(residualize(satArr, [mhiArr, wigArr]), residYmhiWig);
const tSatAfterMHIWig = Math.abs(rSatAfterMHIWig) * Math.sqrt(N-4) / Math.sqrt(1-rSatAfterMHIWig**2+1e-10);

const rRaw = corrWith(satArr, deltaA0);
const tRaw = Math.abs(rRaw) * Math.sqrt(N-2) / Math.sqrt(1-rRaw**2+1e-10);

const isoDL = dlMeta(isoA0, isoIdx.map(i => galaxyData[i].se));
const satDL = dlMeta(satA0, satIdx.map(i => galaxyData[i].se));
const tauRed = Math.min(isoDL.tau, satDL.tau) / fullDL.tau;

const criteria = [
  { name: 'Clear raw |t|>=2', pass: tRaw >= 2.0, val: 't=' + tRaw.toFixed(1) },
  { name: 'Independent of MHI', pass: tSatAfterMHI >= 1.65, val: 'r=' + rSatAfterMHI.toFixed(3) + ' t=' + tSatAfterMHI.toFixed(1) },
  { name: 'Independent of MHI+wig', pass: tSatAfterMHIWig >= 1.65, val: 'r=' + rSatAfterMHIWig.toFixed(3) + ' t=' + tSatAfterMHIWig.toFixed(1) },
  { name: 'LOO-CV > M0', pass: (gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_sat')).cvRMS) / gap * 100 : 0) > 2, val: ((gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_sat')).cvRMS) / gap * 100 : 0)).toFixed(1) + '%' },
  { name: 'Adds to MHI+wig', pass: satAdds, val: (mhiWigSatGap - mhiWigGap).toFixed(1) + '%' },
  { name: 'tau reduction > 5%', pass: (1-tauRed)*100 > 5, val: ((1-tauRed)*100).toFixed(1) + '%' },
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
log("  PHASE 26 — VERDICT");
log("=".repeat(80));
log("");

if (nPassed >= 4) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  SATELLITE vs CENTRAL is a significant independent variable.           ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nPassed >= 2) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  SATELLITE vs CENTRAL shows PARTIAL signal.                            ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  SATELLITE vs CENTRAL FAILS.                                           ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  KEY QUESTION ANSWERED:");
if (tSatAfterEnv < 1.65) {
  log("  satellite/central does NOT add beyond group membership.");
  log("  The Phase 25 envCode is the better environmental variable.");
} else {
  log("  satellite/central ADDS beyond group membership.");
  log("  There are multiple environmental channels.");
}
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  roleCounts: { isolated: nIso, central: nCen, satellite: nSat, unknown: nUnk },
  isoVsSat: { splitDex: +isoSatSplit.toFixed(3), splitT: +isoSatT.toFixed(2) },
  satAfterMHI: { r: +rSatAfterMHI.toFixed(3), t: +tSatAfterMHI.toFixed(1) },
  satAfterMHIWig: { r: +rSatAfterMHIWig.toFixed(3), t: +tSatAfterMHIWig.toFixed(1) },
  satAbsorbedByEnv: tSatAfterEnv < 1.65,
  envAbsorbedBySat: tEnvAfterSat < 1.65,
  looGapClosed: {
    mhi: +mhiGap.toFixed(1), mhiWig: +mhiWigGap.toFixed(1),
    mhiSat: +mhiSatGap.toFixed(1), mhiWigSat: +mhiWigSatGap.toFixed(1),
    mhiWigEnv: +mhiWigEnvGap.toFixed(1), mhiWigRole: +mhiWigRoleGap.toFixed(1)
  },
  wigAfterSat: { r: +rWigAfterSat.toFixed(3), t: +tWigAfterSat.toFixed(1) },
  successCriteria: { passed: nPassed, total: 6 },
  verdict: nPassed >= 4 ? 'SIGNIFICANT' : nPassed >= 2 ? 'PARTIAL' : 'FAILED',
  galaxyRoles: galaxyData.map(g => ({ name: g.name, role: g.role, host: g.host }))
};

fs.writeFileSync(path.join(__dirname, '../public/phase26-satellite-central.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase26-satellite-central.json");
