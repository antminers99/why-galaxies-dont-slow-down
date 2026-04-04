#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "28.0.0";
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

const tidalData = {
  'NGC3769':     { dNN: 0.15, Mneighbor: 2.0e10, neighbor: 'NGC3893', conf: 'high',
                   note: 'UMa pair, NGC3893 is comparable mass at ~0.15 Mpc' },
  'NGC3893':     { dNN: 0.15, Mneighbor: 1.5e10, neighbor: 'NGC3769', conf: 'high',
                   note: 'UMa pair, NGC3769 smaller companion' },
  'NGC4217':     { dNN: 0.20, Mneighbor: 5.0e10, neighbor: 'NGC4226/cluster', conf: 'medium',
                   note: 'UMa cluster, multiple nearby massive galaxies' },
  'NGC4013':     { dNN: 0.30, Mneighbor: 8.0e10, neighbor: 'NGC3992(M109)', conf: 'medium',
                   note: 'UMa cluster, NGC3992 is massive spiral nearby' },
  'NGC4157':     { dNN: 0.30, Mneighbor: 3.0e10, neighbor: 'NGC4100', conf: 'medium',
                   note: 'UMa cluster members, similar mass pair' },
  'NGC5005':     { dNN: 0.30, Mneighbor: 5.0e10, neighbor: 'NGC5033', conf: 'high',
                   note: 'CVn I, paired with NGC5033' },
  'NGC5033':     { dNN: 0.30, Mneighbor: 8.0e10, neighbor: 'NGC5005', conf: 'high',
                   note: 'CVn I, NGC5005 is more massive partner' },
  'NGC0024':     { dNN: 0.35, Mneighbor: 5.0e11, neighbor: 'NGC253', conf: 'high',
                   note: 'Sculptor group, NGC253 is very massive starburst' },
  'NGC4100':     { dNN: 0.35, Mneighbor: 4.0e10, neighbor: 'NGC4157', conf: 'medium',
                   note: 'UMa cluster, near NGC4157' },
  'NGC3726':     { dNN: 0.40, Mneighbor: 8.0e10, neighbor: 'NGC3992(M109)', conf: 'high',
                   note: 'UMa cluster, NGC3992 dominant nearby' },
  'NGC4138':     { dNN: 0.40, Mneighbor: 3.0e10, neighbor: 'NGC4051', conf: 'medium',
                   note: 'UMa cluster, near Seyfert NGC4051' },
  'UGC06973':    { dNN: 0.40, Mneighbor: 8.0e10, neighbor: 'NGC3992(M109)', conf: 'medium',
                   note: 'UMa periphery, NGC3992 dominant' },
  'NGC5055':     { dNN: 0.50, Mneighbor: 3.0e10, neighbor: 'M51(NGC5194)', conf: 'high',
                   note: 'M51 subgroup, interacting pair NGC5194/5195 nearby' },
  'NGC5907':     { dNN: 0.50, Mneighbor: 2.0e10, neighbor: 'NGC5866', conf: 'high',
                   note: 'NGC5866/Draco group, NGC5866 is S0' },
  'NGC2403':     { dNN: 0.60, Mneighbor: 3.0e11, neighbor: 'M81', conf: 'high',
                   note: 'M81 group, M81 is very massive' },
  'NGC4559':     { dNN: 0.80, Mneighbor: 1.0e11, neighbor: 'NGC4565', conf: 'medium',
                   note: 'Coma I cloud, NGC4565 massive edge-on' },
  'NGC5371':     { dNN: 1.00, Mneighbor: 5.0e10, neighbor: 'NGC5353/4', conf: 'medium',
                   note: 'CVn I cloud, NGC5353/4 group nearby' },
  'NGC0891':     { dNN: 1.50, Mneighbor: 5.0e10, neighbor: 'NGC1023', conf: 'medium',
                   note: 'NGC1023 group, loose association' },
  'NGC7331':     { dNN: 2.00, Mneighbor: 3.0e10, neighbor: 'NGC7217', conf: 'medium',
                   note: 'Small group, relatively isolated dominant' },
  'NGC7814':     { dNN: 2.50, Mneighbor: 1.0e11, neighbor: 'NGC7331', conf: 'medium',
                   note: 'Near NGC7331, similar distance' },

  'NGC3198':     { dNN: 2.00, Mneighbor: 1.0e10, neighbor: 'NGC3319', conf: 'low', note: 'UMa periphery' },
  'UGC01281':    { dNN: 2.00, Mneighbor: 5.0e9,  neighbor: 'NGC672', conf: 'low', note: 'Sparse region' },
  'NGC3741':     { dNN: 2.50, Mneighbor: 1.0e10, neighbor: 'NGC4214', conf: 'low', note: 'Isolated dwarf' },
  'UGC02953':    { dNN: 2.50, Mneighbor: 5.0e10, neighbor: 'NGC2768', conf: 'low', note: 'Near NGC2768 group?' },
  'UGC05721':    { dNN: 2.50, Mneighbor: 5.0e9,  neighbor: 'NGC3274', conf: 'low', note: 'Dwarf field' },
  'UGC08490':    { dNN: 2.50, Mneighbor: 1.0e10, neighbor: 'NGC4214', conf: 'low', note: 'Dwarf CVn void' },
  'NGC1003':     { dNN: 3.00, Mneighbor: 5.0e10, neighbor: 'NGC891', conf: 'low', note: 'Perseus field' },
  'NGC1705':     { dNN: 3.00, Mneighbor: 2.0e10, neighbor: 'NGC1559', conf: 'low', note: 'Isolated starburst' },
  'NGC2683':     { dNN: 3.00, Mneighbor: 5.0e10, neighbor: 'NGC2841', conf: 'low', note: 'Lynx field' },
  'NGC2841':     { dNN: 3.00, Mneighbor: 3.0e10, neighbor: 'NGC2683', conf: 'low', note: 'Isolated' },
  'NGC6015':     { dNN: 3.00, Mneighbor: 5.0e10, neighbor: 'NGC5907', conf: 'low', note: 'Relatively isolated' },
  'NGC0289':     { dNN: 3.50, Mneighbor: 5.0e11, neighbor: 'NGC253', conf: 'low', note: 'Sculptor but far' },
  'NGC2903':     { dNN: 3.50, Mneighbor: 3.0e10, neighbor: 'NGC2683', conf: 'low', note: 'Leo field' },
  'NGC3521':     { dNN: 3.50, Mneighbor: 5.0e10, neighbor: 'NGC3627', conf: 'low', note: 'Leo Triplet far' },
  'UGC03580':    { dNN: 3.50, Mneighbor: 5.0e10, neighbor: 'NGC2768', conf: 'low', note: 'Field' },
  'UGC06787':    { dNN: 3.50, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=21' },
  'NGC1090':     { dNN: 4.00, Mneighbor: 1.0e11, neighbor: 'NGC1068', conf: 'low', note: 'Near M77' },
  'NGC2915':     { dNN: 4.00, Mneighbor: 5.0e9,  neighbor: 'NGC3109', conf: 'low', note: 'Isolated BCD' },
  'UGC03546':    { dNN: 4.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=29' },
  'UGC06786':    { dNN: 4.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=29' },
  'NGC6503':     { dNN: 5.00, Mneighbor: 5.0e10, neighbor: 'NGC5907', conf: 'low', note: 'Local Void edge' },
  'NGC6674':     { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Hercules field' },
  'ESO563-G021': { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Southern field D=61' },
  'F571-8':      { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'LSB field D=53' },
  'IC4202':      { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Distant D=100' },
  'NGC0801':     { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=81' },
  'NGC2955':     { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Distant D=98' },
  'NGC2998':     { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=68' },
  'NGC6195':     { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Distant D=128' },
  'UGC00128':    { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'LSB field D=65' },
  'UGC02885':    { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Giant LSB D=81' },
  'UGC02916':    { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=65' },
  'UGC03205':    { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=50' },
  'UGC08699':    { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=39' },
  'UGC09037':    { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=84' },
  'UGC09133':    { dNN: 5.00, Mneighbor: 3.0e10, neighbor: 'unknown', conf: 'low', note: 'Field D=57' },
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

  const td = tidalData[gal.name];
  const dNN = td ? td.dNN : 5.0;
  const Mneighbor = td ? td.Mneighbor : 3.0e10;
  const conf = td ? td.conf : 'low';

  const tidalStrength = Math.log10(Mneighbor / (dNN ** 3));
  const tidalIndex = Math.log10(Mneighbor / 1e10) - 3 * Math.log10(dNN);

  const ev = envLookup[gal.name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    n, logMHI, logVflat, rcWiggliness,
    dNN, Mneighbor, tidalStrength, tidalIndex, conf, envCode,
    neighbor: td ? td.neighbor : 'unknown',
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
log("  PHASE 28: TIDAL INDEX / TIDAL STRENGTH");
log("  Θ = log10(M_neighbor / d³) — combines mass AND distance");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  DEFINITION:");
sep();
log("  Tidal strength Θ = log10(M_neighbor / d_NN³)");
log("  where M_neighbor = estimated stellar mass of nearest massive neighbor");
log("  and d_NN = 3D separation in Mpc");
log("  Higher Θ = stronger tidal influence from environment");
log("  Following Karachentsev & Makarov (1999), Karachentsev+2004");
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

log("  PER-GALAXY (sorted by tidal strength):");
const sorted = [...galaxyData].sort((a,b) => b.tidalStrength - a.tidalStrength);
for (const g of sorted) {
  log("    " + g.name.padEnd(16) + "Θ=" + g.tidalStrength.toFixed(2).padStart(6) +
    "  dNN=" + g.dNN.toFixed(2).padStart(5) + "  logM=" + Math.log10(g.Mneighbor).toFixed(1).padStart(5) +
    "  " + g.conf.padEnd(7) + " logA0=" + g.logA0.toFixed(3));
}
log("");

const tidalArr = galaxyData.map(g => g.tidalStrength);
const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const envArr = galaxyData.map(g => g.envCode);
const vflatArr = galaxyData.map(g => g.logVflat);
const nArr = galaxyData.map(g => g.n);

log("=".repeat(80));
log("  TEST 1: RAW CORRELATIONS");
log("=".repeat(80));
log("");

const rTidRaw = corrWith(tidalArr, deltaA0);
const tTidRaw = Math.abs(rTidRaw) * Math.sqrt(N-2) / Math.sqrt(1-rTidRaw**2+1e-10);

log("  Θ vs delta_a0: r=" + rTidRaw.toFixed(3) + " t=" + tTidRaw.toFixed(1) +
  " → " + (rTidRaw < 0 ? "stronger tidal→lower a0" : "stronger tidal→higher a0"));
log("");

const rTidEnv = corrWith(tidalArr, envArr);
log("  Θ vs envCode: r=" + rTidEnv.toFixed(3) + " (collinearity with group membership)");
log("  " + (Math.abs(rTidEnv) > 0.6 ? "HIGH collinearity" : Math.abs(rTidEnv) > 0.4 ? "MODERATE collinearity" : "LOW collinearity"));
log("");

const rTidWig = corrWith(tidalArr, wigArr);
log("  Θ vs rcWiggliness: r=" + rTidWig.toFixed(3));
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

log("  Θ (tidal strength) vs delta_a0:");
log("  ┌──────────────────────────────────────────────────────────────┐");
log("  │  Controls                 r_partial   |t|    status       │");
log("  ├──────────────────────────────────────────────────────────────┤");
for (const cs of controlSets) {
  const residY = residualize(deltaA0, cs.vals);
  const residX = residualize(tidalArr, cs.vals);
  const r = corrWith(residX, residY);
  const df = N - 2 - cs.vals.length;
  const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
  const status = t >= 2.0 ? "SURVIVES" : t >= 1.65 ? "MARGINAL" : "FAILS";
  log("  │  " + cs.name.padEnd(24) + r.toFixed(3).padStart(8) + t.toFixed(1).padStart(7) + "    " + status.padEnd(10) + "│");
}
log("  └──────────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  TEST 3: Θ AFTER envCode — DOES IT ADD?");
log("=".repeat(80));
log("");

const residAfterEnv = residualize(deltaA0, [mhiArr, wigArr, envArr]);
const rTidAfterEnv = corrWith(residualize(tidalArr, [mhiArr, wigArr, envArr]), residAfterEnv);
const tTidAfterEnv = Math.abs(rTidAfterEnv) * Math.sqrt(N-5) / Math.sqrt(1-rTidAfterEnv**2+1e-10);

const residAfterTid = residualize(deltaA0, [mhiArr, wigArr, tidalArr]);
const rEnvAfterTid = corrWith(residualize(envArr, [mhiArr, wigArr, tidalArr]), residAfterTid);
const tEnvAfterTid = Math.abs(rEnvAfterTid) * Math.sqrt(N-5) / Math.sqrt(1-rEnvAfterTid**2+1e-10);

log("  Θ AFTER envCode+MHI+wig: r=" + rTidAfterEnv.toFixed(3) + " t=" + tTidAfterEnv.toFixed(1) +
  " → " + (tTidAfterEnv >= 1.65 ? "ADDS beyond envCode" : "ABSORBED by envCode"));
log("  envCode AFTER Θ+MHI+wig: r=" + rEnvAfterTid.toFixed(3) + " t=" + tEnvAfterTid.toFixed(1) +
  " → " + (tEnvAfterTid >= 1.65 ? "envCode STILL adds" : "envCode ABSORBED by Θ"));
log("");

if (tTidAfterEnv < 1.65 && tEnvAfterTid >= 1.65) {
  log("  → envCode ABSORBS Θ. Group membership is the deeper variable.");
} else if (tTidAfterEnv >= 1.65 && tEnvAfterTid < 1.65) {
  log("  → Θ ABSORBS envCode! Tidal strength is the deeper physical variable.");
} else if (tTidAfterEnv >= 1.65 && tEnvAfterTid >= 1.65) {
  log("  → BOTH carry independent info.");
} else {
  log("  → Neither adds when the other is present — they measure the same thing.");
}
log("");

log("=".repeat(80));
log("  TEST 4: PERMUTATION TESTS");
log("=".repeat(80));
log("");

const NPERMS = 10000;

for (const [label, ctrls] of [
  ['Θ|MHI', [mhiArr]],
  ['Θ|MHI+wig', [mhiArr, wigArr]],
  ['Θ|MHI+wig+env', [mhiArr, wigArr, envArr]],
]) {
  const residY = residualize(deltaA0, ctrls);
  const residX = residualize(tidalArr, ctrls);
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
  log("  " + label.padEnd(20) + " |r|=" + obsR.toFixed(3) + " perm_p=" + (cnt/NPERMS).toFixed(4) +
    (cnt/NPERMS < 0.05 ? " *" : "") + (cnt/NPERMS < 0.01 ? "*" : ""));
}
log("");

log("=".repeat(80));
log("  TEST 5: LOO CROSS-VALIDATION");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr, 'rcWiggliness': wigArr, 'tidalStrength': tidalArr,
  'envCode': envArr, 'logVflat': vflatArr,
};

function getVal(gIdx, fname) { return allFeatures[fname][gIdx]; }

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_tid: Θ', features: ['tidalStrength'] },
  { name: 'M_env: envCode', features: ['envCode'] },
  { name: 'M_MHI', features: ['logMHI'] },
  { name: 'M_MHI+tid', features: ['logMHI', 'tidalStrength'] },
  { name: 'M_MHI+env', features: ['logMHI', 'envCode'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
  { name: 'M_MHI+wig+tid', features: ['logMHI', 'rcWiggliness', 'tidalStrength'] },
  { name: 'M_MHI+wig+env', features: ['logMHI', 'rcWiggliness', 'envCode'] },
  { name: 'M_MHI+wig+env+tid', features: ['logMHI', 'rcWiggliness', 'envCode', 'tidalStrength'] },
  { name: 'M_MHI+wig+tid(repl)', features: ['logMHI', 'rcWiggliness', 'tidalStrength'] },
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
const mhiWigTidGap = gap > 0 ? (m0rms - fn('M_MHI+wig+tid').cvRMS) / gap * 100 : 0;
const mhiWigEnvGap = gap > 0 ? (m0rms - fn('M_MHI+wig+env').cvRMS) / gap * 100 : 0;
const mhiWigEnvTidGap = gap > 0 ? (m0rms - fn('M_MHI+wig+env+tid').cvRMS) / gap * 100 : 0;
const tidAlone = gap > 0 ? (m0rms - fn('M_tid').cvRMS) / gap * 100 : 0;
const envAlone = gap > 0 ? (m0rms - fn('M_env').cvRMS) / gap * 100 : 0;

log("  KEY COMPARISONS:");
log("    Θ alone:                  " + tidAlone.toFixed(1) + "%");
log("    envCode alone:            " + envAlone.toFixed(1) + "%");
log("    MHI + wig:                " + mhiWigGap.toFixed(1) + "%");
log("    MHI + wig + Θ:            " + mhiWigTidGap.toFixed(1) + "%");
log("    MHI + wig + env:          " + mhiWigEnvGap.toFixed(1) + "%");
log("    MHI + wig + env + Θ:      " + mhiWigEnvTidGap.toFixed(1) + "%");
log("");

const tidAdds = mhiWigTidGap > mhiWigGap + 1;
const tidAddsToEnv = mhiWigEnvTidGap > mhiWigEnvGap + 1;
const tidReplacesEnv = mhiWigTidGap > mhiWigEnvGap + 1;
log("  Θ adds beyond MHI+wig? " + (tidAdds ? "YES (+" + (mhiWigTidGap-mhiWigGap).toFixed(1) + "%)" : "NO"));
log("  Θ adds beyond MHI+wig+env? " + (tidAddsToEnv ? "YES" : "NO"));
log("  Θ replaces envCode? " + (tidReplacesEnv ? "YES (Θ alone > env alone)" : "NO"));
log("");

log("=".repeat(80));
log("  TEST 6: HIGH-CONFIDENCE SUBSAMPLE (n=" + (nHigh + nMed) + ")");
log("=".repeat(80));
log("");

const hiMedIdx = galaxyData.map((g,i) => (g.conf === 'high' || g.conf === 'medium') ? i : -1).filter(i => i >= 0);
const hiN = hiMedIdx.length;
const hiDA0 = hiMedIdx.map(i => deltaA0[i]);
const hiTid = hiMedIdx.map(i => tidalArr[i]);
const hiMHI = hiMedIdx.map(i => mhiArr[i]);
const hiWig = hiMedIdx.map(i => wigArr[i]);

const rHiRaw = corrWith(hiTid, hiDA0);
const tHiRaw = Math.abs(rHiRaw) * Math.sqrt(hiN-2) / Math.sqrt(1-rHiRaw**2+1e-10);

const rHiAfterMHI = corrWith(residualize(hiTid, [hiMHI]), residualize(hiDA0, [hiMHI]));
const tHiAfterMHI = Math.abs(rHiAfterMHI) * Math.sqrt(hiN-3) / Math.sqrt(1-rHiAfterMHI**2+1e-10);

const rHiAfterMHIWig = corrWith(residualize(hiTid, [hiMHI, hiWig]), residualize(hiDA0, [hiMHI, hiWig]));
const tHiAfterMHIWig = Math.abs(rHiAfterMHIWig) * Math.sqrt(hiN-4) / Math.sqrt(1-rHiAfterMHIWig**2+1e-10);

log("  n=" + hiN + " galaxies with high/medium confidence tidal estimates");
log("  Θ vs delta_a0 raw:       r=" + rHiRaw.toFixed(3) + " t=" + tHiRaw.toFixed(1));
log("  Θ vs delta_a0 |MHI:      r=" + rHiAfterMHI.toFixed(3) + " t=" + tHiAfterMHI.toFixed(1));
log("  Θ vs delta_a0 |MHI+wig:  r=" + rHiAfterMHIWig.toFixed(3) + " t=" + tHiAfterMHIWig.toFixed(1));
log("");

log("=".repeat(80));
log("  PHASE 28 — STRICT CRITERIA");
log("=".repeat(80));
log("");

const rAfterMHI = corrWith(residualize(tidalArr, [mhiArr]), residualize(deltaA0, [mhiArr]));
const tAfterMHI = Math.abs(rAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rAfterMHI**2+1e-10);
const rAfterMHIWig = corrWith(residualize(tidalArr, [mhiArr, wigArr]), residualize(deltaA0, [mhiArr, wigArr]));
const tAfterMHIWig = Math.abs(rAfterMHIWig) * Math.sqrt(N-4) / Math.sqrt(1-rAfterMHIWig**2+1e-10);

const criteria = [
  { name: 'Clear raw |t|>=2', pass: tTidRaw >= 2.0, val: 'r=' + rTidRaw.toFixed(3) + ' t=' + tTidRaw.toFixed(1) },
  { name: 'Independent of MHI', pass: tAfterMHI >= 1.65, val: 'r=' + rAfterMHI.toFixed(3) + ' t=' + tAfterMHI.toFixed(1) },
  { name: 'Independent of MHI+wig', pass: tAfterMHIWig >= 1.65, val: 'r=' + rAfterMHIWig.toFixed(3) + ' t=' + tAfterMHIWig.toFixed(1) },
  { name: 'LOO-CV > M0', pass: tidAlone > 2, val: tidAlone.toFixed(1) + '%' },
  { name: 'Adds to MHI+wig', pass: tidAdds, val: (mhiWigTidGap - mhiWigGap).toFixed(1) + '%' },
  { name: 'Adds beyond envCode', pass: tidAddsToEnv, val: (mhiWigEnvTidGap - mhiWigEnvGap).toFixed(1) + '%' },
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
log("  PHASE 28 — VERDICT");
log("=".repeat(80));
log("");

if (nPassed >= 4) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  TIDAL INDEX shows significant independent signal.                     ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nPassed >= 2) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  TIDAL INDEX shows PARTIAL signal.                                     ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  TIDAL INDEX FAILS.                                                    ║");
  log("  ║  Tidal strength does not predict a0 variation independently.           ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  ENVIRONMENT DOOR — UPDATED SCORECARD:");
log("    group membership (envCode):    6/6 PASS ← best");
log("    satellite vs central:          1/6 FAIL");
log("    nearest neighbor distance:     1/6 FAIL");
log("    tidal index (Θ):               " + nPassed + "/6 " + (nPassed >= 4 ? "PASS" : nPassed >= 2 ? "PARTIAL" : "FAIL"));
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  dataQuality: { high: nHigh, medium: nMed, low: nLow },
  rawCorrelation: { r: +rTidRaw.toFixed(3), t: +tTidRaw.toFixed(1) },
  afterMHI: { r: +rAfterMHI.toFixed(3), t: +tAfterMHI.toFixed(1) },
  afterMHIWig: { r: +rAfterMHIWig.toFixed(3), t: +tAfterMHIWig.toFixed(1) },
  tidalVsEnvCode: { r: +rTidEnv.toFixed(3) },
  tidalAfterEnv: { r: +rTidAfterEnv.toFixed(3), t: +tTidAfterEnv.toFixed(1) },
  envAfterTidal: { r: +rEnvAfterTid.toFixed(3), t: +tEnvAfterTid.toFixed(1) },
  highConfSubsample: { n: hiN, rawR: +rHiRaw.toFixed(3), rawT: +tHiRaw.toFixed(1) },
  looGapClosed: {
    tidal: +tidAlone.toFixed(1), env: +envAlone.toFixed(1),
    mhiWig: +mhiWigGap.toFixed(1), mhiWigTid: +mhiWigTidGap.toFixed(1),
    mhiWigEnv: +mhiWigEnvGap.toFixed(1), mhiWigEnvTid: +mhiWigEnvTidGap.toFixed(1),
  },
  successCriteria: { passed: nPassed, total: 6 },
  verdict: nPassed >= 4 ? 'SIGNIFICANT' : nPassed >= 2 ? 'PARTIAL' : 'FAILED',
};

fs.writeFileSync(path.join(__dirname, '../public/phase28-tidal-index.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase28-tidal-index.json");
