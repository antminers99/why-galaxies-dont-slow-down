#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "25.0.0";
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

const knownGroups = {
  'NGC2403': { group: 'M81', env: 'group', note: 'M81 group member, D~3.2 Mpc' },
  'NGC3726': { group: 'UMa', env: 'cluster', note: 'Ursa Major cluster' },
  'NGC3769': { group: 'UMa', env: 'cluster', note: 'Ursa Major cluster, near NGC3893' },
  'NGC3893': { group: 'UMa', env: 'cluster', note: 'Ursa Major cluster, NGC3893 subgroup' },
  'NGC4013': { group: 'UMa', env: 'cluster', note: 'Ursa Major cluster' },
  'NGC4100': { group: 'UMa', env: 'cluster', note: 'Ursa Major cluster' },
  'NGC4138': { group: 'UMa', env: 'cluster', note: 'Ursa Major / CVn border' },
  'NGC4157': { group: 'UMa', env: 'cluster', note: 'Ursa Major cluster' },
  'NGC4217': { group: 'UMa', env: 'cluster', note: 'Ursa Major cluster, near NGC4226' },
  'NGC4559': { group: 'ComaI', env: 'cluster', note: 'Coma I cloud / Virgo periphery' },
  'NGC5005': { group: 'CVnI', env: 'group', note: 'CVn I cloud, paired with NGC5033' },
  'NGC5033': { group: 'CVnI', env: 'group', note: 'CVn I cloud, paired with NGC5005' },
  'NGC5055': { group: 'M51', env: 'group', note: 'M51 group / CVn I cloud' },
  'NGC7331': { group: 'N7331', env: 'group', note: 'NGC7331 group (Deer Lick group)' },
  'NGC0891': { group: 'N1023', env: 'group', note: 'NGC1023 group / loose association' },
  'NGC2841': { group: 'field', env: 'field', note: 'Relatively isolated, D~14.1 Mpc' },
  'NGC2903': { group: 'field', env: 'field', note: 'Leo field, relatively isolated' },
  'NGC3198': { group: 'field', env: 'field', note: 'UMa periphery but not cluster member' },
  'NGC3521': { group: 'field', env: 'field', note: 'Leo field, relatively isolated' },
  'NGC5371': { group: 'CVnI', env: 'group', note: 'CVn I cloud' },
  'NGC5907': { group: 'N5866', env: 'group', note: 'NGC5866/Draco group' },
  'NGC6503': { group: 'field', env: 'field', note: 'Local Void edge, very isolated' },
  'NGC6015': { group: 'field', env: 'field', note: 'Relatively isolated' },
  'NGC2683': { group: 'field', env: 'field', note: 'Field galaxy in Lynx' },
  'NGC1003': { group: 'field', env: 'field', note: 'Field galaxy in Perseus' },
  'NGC1090': { group: 'field', env: 'field', note: 'Field galaxy in Cetus' },
  'NGC0289': { group: 'field', env: 'field', note: 'Field galaxy in Sculptor' },
  'NGC0024': { group: 'Sculptor', env: 'group', note: 'Sculptor group' },
  'NGC2915': { group: 'field', env: 'field', note: 'Isolated blue compact dwarf' },
  'NGC3741': { group: 'field', env: 'field', note: 'Isolated dwarf irregular' },
  'NGC1705': { group: 'field', env: 'field', note: 'Isolated dwarf starburst' },
  'NGC2955': { group: 'field', env: 'field', note: 'Distant field galaxy D~96 Mpc' },
  'NGC2998': { group: 'field', env: 'field', note: 'Field galaxy D~67 Mpc' },
  'NGC6674': { group: 'field', env: 'field', note: 'Field galaxy in Hercules' },
  'NGC6195': { group: 'field', env: 'field', note: 'Distant field galaxy' },
  'NGC7814': { group: 'field', env: 'field', note: 'Relatively isolated, edge-on' },
  'ESO563-G021': { group: 'field', env: 'field', note: 'Southern field galaxy D~61 Mpc' },
  'F571-8': { group: 'field', env: 'field', note: 'LSB field galaxy D~53 Mpc' },
  'IC4202': { group: 'field', env: 'field', note: 'Distant field galaxy D~100 Mpc' },
  'UGC00128': { group: 'field', env: 'field', note: 'LSB field galaxy' },
  'UGC01281': { group: 'field', env: 'field', note: 'Nearby edge-on dwarf' },
  'UGC02885': { group: 'field', env: 'field', note: 'Giant LSB, very isolated' },
  'UGC02916': { group: 'field', env: 'field', note: 'Field galaxy' },
  'UGC02953': { group: 'field', env: 'field', note: 'Field galaxy' },
  'UGC03205': { group: 'field', env: 'field', note: 'Field galaxy' },
  'UGC03546': { group: 'field', env: 'field', note: 'Field galaxy' },
  'UGC03580': { group: 'field', env: 'field', note: 'Field galaxy' },
  'UGC05721': { group: 'field', env: 'field', note: 'Dwarf field galaxy' },
  'UGC06786': { group: 'field', env: 'field', note: 'Field galaxy' },
  'UGC06787': { group: 'field', env: 'field', note: 'Field galaxy' },
  'UGC06973': { group: 'UMa', env: 'cluster', note: 'Ursa Major cluster periphery' },
  'UGC08490': { group: 'field', env: 'field', note: 'Dwarf field galaxy (IZw36)' },
  'UGC08699': { group: 'field', env: 'field', note: 'Field galaxy' },
  'UGC09037': { group: 'field', env: 'field', note: 'Field galaxy' },
  'UGC09133': { group: 'field', env: 'field', note: 'Field galaxy' },
};

const p11 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase11-sensitivity-lab.json'), 'utf8'));
const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));

const galaxyData = [];
for (const gal of p11.galaxies) {
  const pts = gal.localProfile.filter(p =>
    isFinite(p.log_g_bar) && isFinite(p.log_g_obs) && p.log_g_obs > p.log_g_bar * 0.5 && p.r > 0
  );
  if (pts.length < 5) continue;
  const sp = sparc.find(s => s.name === gal.name);
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

  const kg = knownGroups[gal.name];
  const env = kg ? kg.env : 'unknown';
  const groupName = kg ? kg.group : 'unknown';
  const envNote = kg ? kg.note : '';

  const isGroupMember = env === 'group' || env === 'cluster' ? 1 : 0;
  const isCluster = env === 'cluster' ? 1 : 0;
  const envCode = env === 'cluster' ? 2 : env === 'group' ? 1 : 0;

  galaxyData.push({
    name: gal.name, pts, fit,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: fit.rms / Math.sqrt(fit.n),
    inc: sp.inc, D: sp.D, T: sp.T, Q: sp.Q,
    n, logMHI, logVflat, rcWiggliness,
    env, groupName, envNote,
    isGroupMember, isCluster, envCode,
  });
}

const N = galaxyData.length;
const allLogA0 = galaxyData.map(g => g.logA0);
const meanLogA0 = allLogA0.reduce((s,v)=>s+v,0)/N;
const deltaA0 = allLogA0.map(v => v - meanLogA0);

const nField = galaxyData.filter(g => g.env === 'field').length;
const nGroup = galaxyData.filter(g => g.env === 'group').length;
const nCluster = galaxyData.filter(g => g.env === 'cluster').length;
const nUnknown = galaxyData.filter(g => g.env === 'unknown').length;

log("");
log("=".repeat(80));
log("  PHASE 25: ENVIRONMENT & NEIGHBORS — GROUP MEMBERSHIP");
log("  First variable in the 'Environment' door");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  GALAXY HISTORY DOOR — OFFICIALLY CLOSED:");
sep();
log("  sSFR: FAILED | age: FAILED | metallicity: FAILED | SFH: FAILED");
log("  merger/disturbance: photometric FAILED, RC irregularity CONFIRMED (within SPARC)");
log("  rcWiggliness = first confirmed second parameter (4/5 validation tests)");
log("");
log("  NOW OPENING: ENVIRONMENT & NEIGHBORS DOOR");
log("  First variable: group membership");
log("");

log("  GROUP MEMBERSHIP CLASSIFICATION:");
sep();
log("");
log("  Sources: Known galaxy group catalogs (Tully 2015, Karachentsev 2004,");
log("  Makarov & Karachentsev 2011, Kourkchi & Tully 2017).");
log("  Applied to our " + N + " SPARC galaxies by name matching.");
log("");
log("  ┌────────────────────────────────────────────────────────────────┐");
log("  │  Environment     N     fraction                              │");
log("  ├────────────────────────────────────────────────────────────────┤");
log("  │  field            " + nField.toString().padStart(2) + "    " + (nField/N*100).toFixed(0) + "%".padEnd(38) + "│");
log("  │  group            " + nGroup.toString().padStart(2) + "    " + (nGroup/N*100).toFixed(0) + "%".padEnd(38) + "│");
log("  │  cluster           " + nCluster.toString().padStart(2) + "    " + (nCluster/N*100).toFixed(0) + "%".padEnd(37) + "│");
if (nUnknown > 0) log("  │  unknown           " + nUnknown.toString().padStart(2) + "    " + (nUnknown/N*100).toFixed(0) + "%".padEnd(37) + "│");
log("  └────────────────────────────────────────────────────────────────┘");
log("");

log("  PER-GALAXY ASSIGNMENTS:");
for (const g of galaxyData) {
  log("    " + g.name.padEnd(16) + g.env.padEnd(10) + g.groupName.padEnd(10) + "D=" + g.D.toFixed(1) + " Mpc");
}
log("");

log("=".repeat(80));
log("  TEST 1: GROUP vs FIELD — a0 DIFFERENCE");
log("=".repeat(80));
log("");

const fieldIdx = galaxyData.map((g,i) => g.env === 'field' ? i : -1).filter(i => i >= 0);
const groupIdx = galaxyData.map((g,i) => g.env === 'group' || g.env === 'cluster' ? i : -1).filter(i => i >= 0);
const clusterIdx = galaxyData.map((g,i) => g.env === 'cluster' ? i : -1).filter(i => i >= 0);

const fieldLogA0 = fieldIdx.map(i => allLogA0[i]);
const groupLogA0 = groupIdx.map(i => allLogA0[i]);

const fieldMean = fieldLogA0.reduce((s,v)=>s+v,0)/fieldLogA0.length;
const groupMean = groupLogA0.reduce((s,v)=>s+v,0)/groupLogA0.length;
const fieldSD = Math.sqrt(fieldLogA0.reduce((s,v)=>s+(v-fieldMean)**2,0)/(fieldLogA0.length-1));
const groupSD = Math.sqrt(groupLogA0.reduce((s,v)=>s+(v-groupMean)**2,0)/(groupLogA0.length-1));
const splitDex = Math.abs(fieldMean - groupMean);
const splitT = splitDex / Math.sqrt(fieldSD**2/fieldIdx.length + groupSD**2/groupIdx.length);

const fieldDL = dlMeta(fieldLogA0, fieldIdx.map(i => galaxyData[i].se));
const groupDL = dlMeta(groupLogA0, groupIdx.map(i => galaxyData[i].se));
const fullDL = dlMeta(allLogA0, galaxyData.map(g => g.se));

log("  ┌────────────────────────────────────────────────────────────────────┐");
log("  │  Sample        N    mean a0    mean logA0   SD       tau        │");
log("  ├────────────────────────────────────────────────────────────────────┤");
log("  │  FULL          " + N.toString().padStart(2) + "   " + Math.round(Math.pow(10, meanLogA0)).toString().padStart(5) + "     " + meanLogA0.toFixed(3) + "     " + Math.sqrt(allLogA0.reduce((s,v)=>s+(v-meanLogA0)**2,0)/(N-1)).toFixed(3) + "   " + fullDL.tau.toFixed(3) + "  │");
log("  │  field         " + fieldIdx.length.toString().padStart(2) + "   " + Math.round(Math.pow(10, fieldMean)).toString().padStart(5) + "     " + fieldMean.toFixed(3) + "     " + fieldSD.toFixed(3) + "   " + fieldDL.tau.toFixed(3) + "  │");
log("  │  group+cluster " + groupIdx.length.toString().padStart(2) + "   " + Math.round(Math.pow(10, groupMean)).toString().padStart(5) + "     " + groupMean.toFixed(3) + "     " + groupSD.toFixed(3) + "   " + groupDL.tau.toFixed(3) + "  │");
log("  └────────────────────────────────────────────────────────────────────┘");
log("");
log("  Split: " + splitDex.toFixed(3) + " dex (t=" + splitT.toFixed(2) + ", p~" + (splitT > 2 ? "<0.05" : ">0.05") + ")");
log("  Direction: " + (groupMean > fieldMean ? "group > field" : "field > group"));
log("  Within-group tau reduction: " + ((1 - Math.min(fieldDL.tau, groupDL.tau) / fullDL.tau) * 100).toFixed(1) + "%");
log("");

if (clusterIdx.length >= 5) {
  const pureGroupIdx = galaxyData.map((g,i) => g.env === 'group' ? i : -1).filter(i => i >= 0);
  log("  THREE-WAY SPLIT (field / group / cluster):");
  for (const [label, idx] of [['field', fieldIdx], ['group', pureGroupIdx], ['cluster', clusterIdx]]) {
    if (idx.length < 3) continue;
    const la = idx.map(i => allLogA0[i]);
    const m = la.reduce((s,v)=>s+v,0)/la.length;
    const sd = Math.sqrt(la.reduce((s,v)=>s+(v-m)**2,0)/(la.length-1));
    const dl = dlMeta(la, idx.map(i => galaxyData[i].se));
    log("    " + label.padEnd(10) + " n=" + idx.length.toString().padStart(2) + "  a0=" + Math.round(Math.pow(10, m)).toString().padStart(5) + "  logA0=" + m.toFixed(3) + "  tau=" + dl.tau.toFixed(3));
  }
  log("");
}

log("=".repeat(80));
log("  TEST 2: CORRELATION — envCode vs delta_a0");
log("=".repeat(80));
log("");

const envCodeArr = galaxyData.map(g => g.envCode);
const isGroupArr = galaxyData.map(g => g.isGroupMember);
const mhiArr = galaxyData.map(g => g.logMHI);
const wigArr = galaxyData.map(g => g.rcWiggliness);
const vflatArr = galaxyData.map(g => g.logVflat);
const nArr = galaxyData.map(g => g.n);

const rEnvRaw = corrWith(envCodeArr, deltaA0);
const tEnvRaw = Math.abs(rEnvRaw) * Math.sqrt(N-2) / Math.sqrt(1-rEnvRaw**2+1e-10);

const rGroupRaw = corrWith(isGroupArr, deltaA0);
const tGroupRaw = Math.abs(rGroupRaw) * Math.sqrt(N-2) / Math.sqrt(1-rGroupRaw**2+1e-10);

log("  envCode (0=field, 1=group, 2=cluster) vs delta_a0:");
log("    r = " + rEnvRaw.toFixed(3) + " (t=" + tEnvRaw.toFixed(1) + ")");
log("");
log("  isGroupMember (0/1) vs delta_a0:");
log("    r = " + rGroupRaw.toFixed(3) + " (t=" + tGroupRaw.toFixed(1) + ")");
log("");

log("=".repeat(80));
log("  TEST 3: INDEPENDENCE — AFTER MHI, AFTER rcWiggliness");
log("=".repeat(80));
log("");

const controls = [
  { name: 'raw', vals: [] },
  { name: '|MHI', vals: [mhiArr] },
  { name: '|MHI+Vflat', vals: [mhiArr, vflatArr] },
  { name: '|MHI+Vflat+n', vals: [mhiArr, vflatArr, nArr] },
  { name: '|MHI+wig', vals: [mhiArr, wigArr] },
  { name: '|MHI+Vflat+wig', vals: [mhiArr, vflatArr, wigArr] },
  { name: '|MHI+Vflat+n+wig', vals: [mhiArr, vflatArr, nArr, wigArr] },
];

log("  envCode vs delta_a0, sequential confounder stripping:");
log("  ┌──────────────────────────────────────────────────────────┐");
log("  │  Controls               r_partial   |t|    status     │");
log("  ├──────────────────────────────────────────────────────────┤");

for (const c of controls) {
  const residY = residualize(deltaA0, c.vals);
  const residX = residualize(envCodeArr, c.vals);
  const r = corrWith(residX, residY);
  const df = N - 2 - c.vals.length;
  const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
  const status = t >= 2.0 ? "STRONG" : t >= 1.65 ? "MARGINAL" : t >= 1.0 ? "WEAK" : "NONE";
  log("  │  " + c.name.padEnd(22) + r.toFixed(3).padStart(8) + t.toFixed(1).padStart(7) + "    " + status.padEnd(10) + "│");
}
log("  └──────────────────────────────────────────────────────────┘");
log("");

log("  isGroupMember vs delta_a0, same controls:");
log("  ┌──────────────────────────────────────────────────────────┐");
log("  │  Controls               r_partial   |t|    status     │");
log("  ├──────────────────────────────────────────────────────────┤");

for (const c of controls) {
  const residY = residualize(deltaA0, c.vals);
  const residX = residualize(isGroupArr, c.vals);
  const r = corrWith(residX, residY);
  const df = N - 2 - c.vals.length;
  const t = Math.abs(r) * Math.sqrt(df) / Math.sqrt(1-r**2+1e-10);
  const status = t >= 2.0 ? "STRONG" : t >= 1.65 ? "MARGINAL" : t >= 1.0 ? "WEAK" : "NONE";
  log("  │  " + c.name.padEnd(22) + r.toFixed(3).padStart(8) + t.toFixed(1).padStart(7) + "    " + status.padEnd(10) + "│");
}
log("  └──────────────────────────────────────────────────────────┘");
log("");

log("=".repeat(80));
log("  TEST 4: PERMUTATION TEST");
log("=".repeat(80));
log("");

const NPERMS = 10000;
for (const [label, arr] of [['envCode', envCodeArr], ['isGroup', isGroupArr]]) {
  const residY = residualize(deltaA0, [mhiArr]);
  const obsR = Math.abs(corrWith(arr, residY));
  let cnt = 0;
  const sh = [...residY];
  for (let p = 0; p < NPERMS; p++) {
    for (let i = sh.length - 1; i > 0; i--) {
      const j = Math.floor(Math.random() * (i + 1));
      [sh[i], sh[j]] = [sh[j], sh[i]];
    }
    if (Math.abs(corrWith(arr, sh)) >= obsR) cnt++;
  }
  log("  " + label.padEnd(12) + " |r|=" + obsR.toFixed(3) + " perm_p=" + (cnt/NPERMS).toFixed(4) +
    (cnt/NPERMS < 0.05 ? " *" : ""));
}
log("");

log("=".repeat(80));
log("  TEST 5: LOO CROSS-VALIDATION");
log("=".repeat(80));
log("");

const allFeatures = {
  'logMHI': mhiArr,
  'rcWiggliness': wigArr,
  'envCode': envCodeArr,
  'isGroup': isGroupArr,
  'logVflat': vflatArr,
};

function getVal(gIdx, fname) { return allFeatures[fname][gIdx]; }

const modelDefs = [
  { name: 'M0: Universal', features: [] },
  { name: 'M_env: envCode', features: ['envCode'] },
  { name: 'M_grp: isGroup', features: ['isGroup'] },
  { name: 'M_MHI', features: ['logMHI'] },
  { name: 'M_wig', features: ['rcWiggliness'] },
  { name: 'M_MHI+env', features: ['logMHI', 'envCode'] },
  { name: 'M_MHI+grp', features: ['logMHI', 'isGroup'] },
  { name: 'M_MHI+wig', features: ['logMHI', 'rcWiggliness'] },
  { name: 'M_MHI+wig+env', features: ['logMHI', 'rcWiggliness', 'envCode'] },
  { name: 'M_MHI+wig+grp', features: ['logMHI', 'rcWiggliness', 'isGroup'] },
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

log("  ┌────────────────────────────────────────────────────────────────────────────┐");
log("  │  Model                k    CV-RMS    vs M0    gap-closed                 │");
log("  ├────────────────────────────────────────────────────────────────────────────┤");
for (const r of cvResults) {
  const vsM0 = ((1 - r.cvRMS / m0rms) * 100).toFixed(1);
  const gc = gap > 0 ? ((m0rms - r.cvRMS) / gap * 100).toFixed(1) : '0.0';
  log("  │  " + r.name.padEnd(21) + r.k.toString().padStart(3) + r.cvRMS.toFixed(5).padStart(10) +
    (vsM0 + "%").padStart(9) + (gc + "%").padStart(14) + "                │");
}
log("  └────────────────────────────────────────────────────────────────────────────┘");
log("");

const envGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_env')).cvRMS) / gap * 100 : 0;
const grpGap = gap > 0 ? (m0rms - cvResults.find(r => r.name.includes('M_grp:')).cvRMS) / gap * 100 : 0;
const mhiGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI').cvRMS) / gap * 100 : 0;
const mhiWigGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig').cvRMS) / gap * 100 : 0;
const mhiWigEnvGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig+env').cvRMS) / gap * 100 : 0;
const mhiWigGrpGap = gap > 0 ? (m0rms - cvResults.find(r => r.name === 'M_MHI+wig+grp').cvRMS) / gap * 100 : 0;

log("  KEY LOO NUMBERS:");
log("    envCode alone:         " + envGap.toFixed(1) + "%");
log("    isGroup alone:         " + grpGap.toFixed(1) + "%");
log("    MHI alone:             " + mhiGap.toFixed(1) + "%");
log("    MHI + wig:             " + mhiWigGap.toFixed(1) + "%");
log("    MHI + wig + env:       " + mhiWigEnvGap.toFixed(1) + "%");
log("    MHI + wig + grp:       " + mhiWigGrpGap.toFixed(1) + "%");
log("");

const envAddsToMHIWig = mhiWigEnvGap > mhiWigGap + 1;
log("  Does environment add beyond MHI + rcWiggliness? " +
  (envAddsToMHIWig ? "YES (+" + (mhiWigEnvGap - mhiWigGap).toFixed(1) + "%)" : "NO"));
log("");

log("=".repeat(80));
log("  TEST 6: DOES ENVIRONMENT EXPLAIN rcWiggliness?");
log("  (Is RC irregularity caused by group environment?)");
log("=".repeat(80));
log("");

const rEnvWig = corrWith(envCodeArr, wigArr);
const rGroupWig = corrWith(isGroupArr, wigArr);
log("  envCode vs rcWiggliness: r=" + rEnvWig.toFixed(3));
log("  isGroup vs rcWiggliness: r=" + rGroupWig.toFixed(3));
log("");

const fieldWig = fieldIdx.map(i => wigArr[i]);
const groupWig = groupIdx.map(i => wigArr[i]);
const fieldWigMean = fieldWig.reduce((s,v)=>s+v,0)/fieldWig.length;
const groupWigMean = groupWig.reduce((s,v)=>s+v,0)/groupWig.length;
log("  Mean rcWiggliness:");
log("    field:         " + fieldWigMean.toFixed(4));
log("    group+cluster: " + groupWigMean.toFixed(4));
log("    " + (groupWigMean > fieldWigMean ? "group galaxies MORE wiggly" : "field galaxies MORE wiggly"));
log("");

log("=".repeat(80));
log("  PHASE 25 — STRICT SUCCESS CRITERIA");
log("=".repeat(80));
log("");

const residYmhi = residualize(deltaA0, [mhiArr]);
const rEnvAfterMHI = corrWith(residualize(envCodeArr, [mhiArr]), residYmhi);
const tEnvAfterMHI = Math.abs(rEnvAfterMHI) * Math.sqrt(N-3) / Math.sqrt(1-rEnvAfterMHI**2+1e-10);

const residYmhiWig = residualize(deltaA0, [mhiArr, wigArr]);
const rEnvAfterMHIWig = corrWith(residualize(envCodeArr, [mhiArr, wigArr]), residYmhiWig);
const tEnvAfterMHIWig = Math.abs(rEnvAfterMHIWig) * Math.sqrt(N-4) / Math.sqrt(1-rEnvAfterMHIWig**2+1e-10);

const tauReduction = Math.min(fieldDL.tau, groupDL.tau) / fullDL.tau;

const criteria = [
  { name: 'Clear correlation |t|>=2', pass: tEnvRaw >= 2.0, val: 't=' + tEnvRaw.toFixed(1) },
  { name: 'Independent of MHI', pass: tEnvAfterMHI >= 1.65, val: 'r=' + rEnvAfterMHI.toFixed(3) + ' t=' + tEnvAfterMHI.toFixed(1) },
  { name: 'Independent of MHI+wig', pass: tEnvAfterMHIWig >= 1.65, val: 'r=' + rEnvAfterMHIWig.toFixed(3) + ' t=' + tEnvAfterMHIWig.toFixed(1) },
  { name: 'LOO-CV > M0', pass: envGap > 2, val: envGap.toFixed(1) + '%' },
  { name: 'Adds to MHI+wig', pass: envAddsToMHIWig, val: (mhiWigEnvGap - mhiWigGap).toFixed(1) + '%' },
  { name: 'tau reduction > 5%', pass: (1-tauReduction)*100 > 5, val: ((1-tauReduction)*100).toFixed(1) + '%' },
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
log("  PHASE 25 — VERDICT");
log("=".repeat(80));
log("");

if (nPassed >= 4) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  GROUP MEMBERSHIP shows significant independent signal.                ║");
  log("  ║  Environment matters for a0 beyond MHI and rcWiggliness.              ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else if (nPassed >= 2) {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  GROUP MEMBERSHIP shows PARTIAL signal. Some criteria met.             ║");
  log("  ║  May warrant further investigation with finer environment measures.    ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════════╗");
  log("  ║  GROUP MEMBERSHIP FAILS. No significant independent signal.            ║");
  log("  ║  Being in a group vs field does not predict a0 variation.              ║");
  log("  ║  Proceed to next environment variable: nearest neighbor / density.     ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════════╝");
}
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  environmentCounts: { field: nField, group: nGroup, cluster: nCluster, unknown: nUnknown },
  fieldVsGroup: {
    fieldMeanA0: Math.round(Math.pow(10, fieldMean)),
    groupMeanA0: Math.round(Math.pow(10, groupMean)),
    splitDex: +splitDex.toFixed(3),
    splitT: +splitT.toFixed(2),
    fieldTau: +fieldDL.tau.toFixed(3),
    groupTau: +groupDL.tau.toFixed(3),
  },
  correlations: {
    envCodeRaw: { r: +rEnvRaw.toFixed(3), t: +tEnvRaw.toFixed(1) },
    envCodeAfterMHI: { r: +rEnvAfterMHI.toFixed(3), t: +tEnvAfterMHI.toFixed(1) },
    envCodeAfterMHIWig: { r: +rEnvAfterMHIWig.toFixed(3), t: +tEnvAfterMHIWig.toFixed(1) },
  },
  envVsWiggliness: { r: +rEnvWig.toFixed(3) },
  looGapClosed: {
    env: +envGap.toFixed(1), group: +grpGap.toFixed(1),
    mhi: +mhiGap.toFixed(1), mhiWig: +mhiWigGap.toFixed(1),
    mhiWigEnv: +mhiWigEnvGap.toFixed(1), mhiWigGrp: +mhiWigGrpGap.toFixed(1)
  },
  successCriteria: { passed: nPassed, total: 6 },
  verdict: nPassed >= 4 ? 'SIGNIFICANT' : nPassed >= 2 ? 'PARTIAL' : 'FAILED',
  galaxyAssignments: galaxyData.map(g => ({ name: g.name, env: g.env, group: g.groupName }))
};

fs.writeFileSync(path.join(__dirname, '../public/phase25-group-membership.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase25-group-membership.json");
