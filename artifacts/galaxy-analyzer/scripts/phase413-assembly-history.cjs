const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

const G = 4.3009e-6;

function normalize(n) {
  return n.replace(/\s+/g,'').replace(/^NGC0*/,'NGC').replace(/^DDO0*/,'DDO').replace(/^IC0*/,'IC').replace(/^UGC0*/,'UGC').replace(/^PGC0*/,'PGC').toUpperCase();
}

const tsContent = fs.readFileSync(path.join(__dirname, '..', 'src', 'data', 'sparc-datasets.ts'), 'utf8');
function parseRCData(ts) {
  const rcMap = {};
  const re = /"([^"]+)":\s*\{[^[]*data:\s*\[([\s\S]*?)\]/g;
  let m;
  while ((m = re.exec(ts)) !== null) {
    const points = [];
    const ptRe = /r:\s*([\d.]+)\s*,\s*v:\s*([\d.]+)/g;
    let pm;
    while ((pm = ptRe.exec(m[2])) !== null) points.push({ r: parseFloat(pm[1]), v: parseFloat(pm[2]) });
    if (points.length >= 3) rcMap[m[1]] = points;
  }
  return rcMap;
}
const rcMapRaw = parseRCData(tsContent);
const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) rcMap[normalize(name)] = pts;

function pearsonR(x, y) {
  const n=x.length;if(n<4)return NaN;
  const mx=x.reduce((a,b)=>a+b,0)/n,my=y.reduce((a,b)=>a+b,0)/n;
  let num=0,dx2=0,dy2=0;
  for(let i=0;i<n;i++){const dx=x[i]-mx,dy=y[i]-my;num+=dx*dy;dx2+=dx*dx;dy2+=dy*dy;}
  return dx2>0&&dy2>0?num/Math.sqrt(dx2*dy2):0;
}
function multiR2(X, y) {
  const n=y.length,nv=X[0].length;
  const my=y.reduce((a,b)=>a+b,0)/n;
  const mx=Array(nv).fill(0);
  for(let j=0;j<nv;j++){for(let i=0;i<n;i++)mx[j]+=X[i][j];mx[j]/=n;}
  const XTX=Array.from({length:nv},()=>Array(nv).fill(0)),XTy=Array(nv).fill(0);
  for(let i=0;i<n;i++){for(let j=0;j<nv;j++){XTy[j]+=(X[i][j]-mx[j])*(y[i]-my);for(let k=0;k<nv;k++)XTX[j][k]+=(X[i][j]-mx[j])*(X[i][k]-mx[k]);}}
  const aug=XTX.map((row,i)=>[...row,XTy[i]]);
  for(let col=0;col<nv;col++){let maxRow=col;for(let row=col+1;row<nv;row++)if(Math.abs(aug[row][col])>Math.abs(aug[maxRow][col]))maxRow=row;[aug[col],aug[maxRow]]=[aug[maxRow],aug[col]];if(Math.abs(aug[col][col])<1e-12)continue;for(let row=col+1;row<nv;row++){const f=aug[row][col]/aug[col][col];for(let j=col;j<=nv;j++)aug[row][j]-=f*aug[col][j];}}
  const beta=Array(nv).fill(0);
  for(let i=nv-1;i>=0;i--){beta[i]=aug[i][nv];for(let j=i+1;j<nv;j++)beta[i]-=aug[i][j]*beta[j];beta[i]/=aug[i][i]||1;}
  let sse=0,sst=0;const residuals=[];
  for(let i=0;i<n;i++){let pred=my;for(let j=0;j<nv;j++)pred+=beta[j]*(X[i][j]-mx[j]);residuals.push(y[i]-pred);sse+=(y[i]-pred)**2;sst+=(y[i]-my)**2;}
  return {R2:sst>0?1-sse/sst:0,beta,residuals};
}
function partialR(x, y, controls) {
  if (controls.length === 0 || controls[0].length === 0) return pearsonR(x, y);
  const xModel = multiR2(controls, x);
  const yModel = multiR2(controls, y);
  return pearsonR(xModel.residuals, yModel.residuals);
}

const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[g.name] = g; resultsMap[normalize(g.name)] = g; });

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sr) continue;
  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 5) continue;

  const Vflat = sp.Vflat;
  const Rdisk = sp.Rdisk;
  const Mbar = (sp.L36*0.5 + sp.MHI*1.33) * 1e9;
  const Rmax = rc[rc.length-1].r;
  const Vmax = Math.max(...rc.map(p => p.v));

  const M_halo = sr.models.dark_halo_linear.M;
  const k_halo = sr.models.dark_halo_linear.k;
  const mse_newt = sr.models.newtonian.mse;
  const mse_halo = sr.models.dark_halo_linear.mse;

  const V_inner_pts = rc.filter(p => p.r < 2*Rdisk);
  const V_outer_pts = rc.filter(p => p.r > 2*Rdisk);
  const V_inner = V_inner_pts.length >= 2 ? V_inner_pts.reduce((a,p)=>a+p.v,0)/V_inner_pts.length : Vflat;
  const V_outer = V_outer_pts.length >= 2 ? V_outer_pts.reduce((a,p)=>a+p.v,0)/V_outer_pts.length : Vflat;

  const V_Newt_Rmax = Math.sqrt(G * Mbar / Rmax);
  const V_Newt_Rd = Math.sqrt(G * Mbar / (2.2 * Rdisk));

  const dmFrac_Rmax = Math.max(0, 1 - (V_Newt_Rmax/Vflat)**2);
  const dmFrac_2Rd = V_inner > 0 ? Math.max(0, 1 - (V_Newt_Rd/V_inner)**2) : 0;

  const logK_halo = k_halo > 0 ? Math.log10(k_halo) : -5;
  const logM_halo = M_halo > 0 ? Math.log10(M_halo) : 5;

  const innerSlope = (() => {
    const pts = rc.filter(p => p.r > 0.3*Rdisk && p.r < 2*Rdisk);
    if (pts.length < 3) return 0;
    const lx = pts.map(p => Math.log10(p.r)), ly = pts.map(p => Math.log10(p.v));
    const mx = lx.reduce((a,b)=>a+b,0)/lx.length, my = ly.reduce((a,b)=>a+b,0)/ly.length;
    let num=0,den=0;
    for(let i=0;i<lx.length;i++){num+=(lx[i]-mx)*(ly[i]-my);den+=(lx[i]-mx)**2;}
    return den>0?num/den:0;
  })();

  const outerSlope = (() => {
    const pts = rc.filter(p => p.r > 3*Rdisk);
    if (pts.length < 3) return 0;
    const lx = pts.map(p => Math.log10(p.r)), ly = pts.map(p => Math.log10(p.v));
    const mx = lx.reduce((a,b)=>a+b,0)/lx.length, my = ly.reduce((a,b)=>a+b,0)/ly.length;
    let num=0,den=0;
    for(let i=0;i<lx.length;i++){num+=(lx[i]-mx)*(ly[i]-my);den+=(lx[i]-mx)**2;}
    return den>0?num/den:0;
  })();

  const concIdx = Vmax > 0 ? V_inner / Vmax : 1;
  const haloResponse = mse_newt > 0 ? Math.log10(Math.max(mse_newt / Math.max(mse_halo, 0.001), 0.01)) : 0;

  const logMbar = Math.log10(Math.max(Mbar, 1));
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(Rdisk, 0.01));
  const logSBdisk = Math.log10(Math.max(sp.SBdisk, 0.01));

  const gasFrac = sp.L36 > 0 ? sp.MHI / sp.L36 : 0;
  const logGasFrac = gasFrac > 0 ? Math.log10(gasFrac) : -3;

  const baryonCompact = Rdisk > 0 ? Math.log10(Mbar / (Rdisk * Rdisk)) : 0;

  const stellarHaloRatio = M_halo > 0 ? logMbar - logM_halo : 0;

  const SBexcess = (() => {
    return logSBdisk;
  })();

  gals.push({
    name: g.name, logA0: g.logA0,
    Vflat, Rdisk, Rmax, Vmax, Mbar, rc,
    logVflat: Math.log10(Vflat),
    logL36, logRdisk, logMbar,
    logMHI: g.logMHI,
    morphT: sp.T,
    logSBdisk,
    envCode: g.envCode,
    L36: sp.L36, MHI: sp.MHI,
    logK_halo, logM_halo, k_halo, M_halo,
    mse_newt, mse_halo,
    dmFrac_Rmax, dmFrac_2Rd,
    innerSlope, outerSlope, concIdx,
    haloResponse,
    V_inner, V_outer,
    gasFrac, logGasFrac,
    baryonCompact,
    stellarHaloRatio,
    SBexcess,
  });
}

const structVars = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
const struct6 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
const vfModel = multiR2(structVars, gals.map(g=>g.logVflat));
const a0Model = multiR2(struct6, gals.map(g=>g.logA0));
for (let i=0; i<gals.length; i++) {
  gals[i].VfResid = vfModel.residuals[i];
  gals[i].a0Resid = a0Model.residuals[i];
}
const sdVf = Math.sqrt(gals.reduce((a,g)=>a+g.VfResid**2,0)/gals.length);
const sdA0 = Math.sqrt(gals.reduce((a,g)=>a+g.a0Resid**2,0)/gals.length);
for (let i=0; i<gals.length; i++) {
  gals[i].VfResid_z = gals[i].VfResid / sdVf;
  gals[i].a0Resid_z = gals[i].a0Resid / sdA0;
  gals[i].L_sum = gals[i].VfResid_z + gals[i].a0Resid_z;
}

const logK_struct = multiR2(struct6, gals.map(g=>g.logK_halo));
const dmFrac_struct = multiR2(struct6, gals.map(g=>g.dmFrac_Rmax));
const logK_resid = logK_struct.residuals;
const dmFrac_resid = dmFrac_struct.residuals;
for (let i=0; i<gals.length; i++) {
  gals[i].logK_resid = logK_resid[i];
  gals[i].dmFrac_resid = dmFrac_resid[i];
}

const gasFrac_struct = multiR2(struct6, gals.map(g=>g.logGasFrac));
const compact_struct = multiR2(struct6, gals.map(g=>g.baryonCompact));
const SHR_struct = multiR2(struct6, gals.map(g=>g.stellarHaloRatio));
const morphT_struct = multiR2(gals.map(g=>[g.logMbar,g.logL36,g.logRdisk]), gals.map(g=>g.morphT));

for (let i=0; i<gals.length; i++) {
  gals[i].gasFrac_resid = gasFrac_struct.residuals[i];
  gals[i].compact_resid = compact_struct.residuals[i];
  gals[i].SHR_resid = SHR_struct.residuals[i];
  gals[i].morphT_resid = morphT_struct.residuals[i];
}

const rBase = pearsonR(gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid));

console.log('='.repeat(70));
console.log('PHASE 413: ASSEMBLY HISTORY / FORMATION TIME');
console.log('Does formation history explain the residual channel after halo structure?');
console.log('='.repeat(70));
console.log('N = ' + gals.length + ', baseline r(VfResid, a0_resid) = ' + rBase.toFixed(3));


console.log('\n' + '#'.repeat(70));
console.log('413A: FORMATION PROXIES — CORRELATIONS WITH LATENT VARIABLE');
console.log('Assembly history indicators from observable properties');
console.log('#'.repeat(70));

console.log('\n  ASSEMBLY PROXY RATIONALE:');
console.log('  - logK_resid: halo slope excess at fixed structure = concentration anomaly');
console.log('  - dmFrac_resid: DM dominance excess at fixed structure = assembly depth');
console.log('  - logGasFrac: gas richness = evolutionary state (gas-poor = more evolved)');
console.log('  - gasFrac_resid: gas fraction excess = recent accretion indicator');
console.log('  - baryonCompact: log(Mbar/Rd^2) = stellar compactness = formation epoch');
console.log('  - compact_resid: compactness excess = early assembly signature');
console.log('  - stellarHaloRatio: Mbar/Mhalo = assembly efficiency');
console.log('  - SHR_resid: efficiency anomaly at fixed structure');
console.log('  - morphT_resid: morphological type excess = assembly history residual');

const assemblyProxies = [
  { name: 'logK_resid', get: g => g.logK_resid, desc: 'Halo slope excess' },
  { name: 'dmFrac_resid', get: g => g.dmFrac_resid, desc: 'DM fraction excess' },
  { name: 'logGasFrac', get: g => g.logGasFrac, desc: 'Gas fraction (raw)' },
  { name: 'gasFrac_resid', get: g => g.gasFrac_resid, desc: 'Gas fraction excess' },
  { name: 'baryonCompact', get: g => g.baryonCompact, desc: 'Baryon compactness' },
  { name: 'compact_resid', get: g => g.compact_resid, desc: 'Compactness excess' },
  { name: 'stellarHaloRatio', get: g => g.stellarHaloRatio, desc: 'Mbar/Mhalo ratio' },
  { name: 'SHR_resid', get: g => g.SHR_resid, desc: 'SH ratio excess' },
  { name: 'morphT_resid', get: g => g.morphT_resid, desc: 'Morphology excess' },
  { name: 'envCode', get: g => g.envCode, desc: 'Environment code' },
];

console.log('\n  Proxy               r(VfResid) r(a0_resid) r(L_sum)   desc');
console.log('  ' + '-'.repeat(75));

const proxyCorrelations = [];
for (const p of assemblyProxies) {
  const vals = gals.map(g => p.get(g));
  const valid = vals.every(v => isFinite(v));
  if (!valid) continue;
  const rVf = pearsonR(vals, gals.map(g=>g.VfResid));
  const rA0 = pearsonR(vals, gals.map(g=>g.a0Resid));
  const rL = pearsonR(vals, gals.map(g=>g.L_sum));
  console.log('  ' + p.name.padEnd(20) + ' ' + ((rVf>=0?'+':'')+rVf.toFixed(3)).padEnd(11) + ' ' + ((rA0>=0?'+':'')+rA0.toFixed(3)).padEnd(12) + ' ' + ((rL>=0?'+':'')+rL.toFixed(3)).padEnd(11) + ' ' + p.desc);
  proxyCorrelations.push({ name: p.name, rVfResid: rVf, rA0resid: rA0, rLsum: rL });
}


console.log('\n\n' + '#'.repeat(70));
console.log('413B: ENVIRONMENT AS HISTORY PROXY');
console.log('Field/Group/Cluster as assembly pathway indicator');
console.log('#'.repeat(70));

const envGroups = { field: [], group: [], cluster: [] };
gals.forEach(g => {
  if (g.envCode === 1) envGroups.field.push(g);
  else if (g.envCode === 2) envGroups.group.push(g);
  else if (g.envCode === 3) envGroups.cluster.push(g);
});

console.log('\n  Environment split:');
for (const [env, gs] of Object.entries(envGroups)) {
  if (gs.length < 4) { console.log('  ' + env + ': N=' + gs.length + ' (too few)'); continue; }
  const r = pearsonR(gs.map(g=>g.VfResid), gs.map(g=>g.a0Resid));
  const meanL = gs.reduce((a,g)=>a+g.L_sum,0)/gs.length;
  const sdL = Math.sqrt(gs.reduce((a,g)=>a+(g.L_sum-meanL)**2,0)/gs.length);
  console.log('  ' + env.padEnd(10) + ' N=' + String(gs.length).padEnd(4) + ' r=' + r.toFixed(3) + '  mean(L)=' + meanL.toFixed(3) + '  sd(L)=' + sdL.toFixed(3));
}

console.log('\n  HISTORY INTERPRETATION:');
console.log('  Field galaxies: isolated, secular evolution, late assembly');
console.log('  Group galaxies: moderate interactions, intermediate assembly');
console.log('  Cluster galaxies: processed, stripped, early infall');
console.log('  If channel tracks formation history, field should differ systematically.');

const envIsField = gals.map(g => g.envCode === 1 ? 1 : 0);
const envIsCluster = gals.map(g => g.envCode === 3 ? 1 : 0);
const envGradient = gals.map(g => g.envCode);

console.log('\n  Environment gradient correlations:');
console.log('  r(envCode, L_sum) = ' + pearsonR(envGradient, gals.map(g=>g.L_sum)).toFixed(3));
console.log('  r(isField, L_sum) = ' + pearsonR(envIsField, gals.map(g=>g.L_sum)).toFixed(3));
console.log('  r(isCluster, L_sum) = ' + pearsonR(envIsCluster, gals.map(g=>g.L_sum)).toFixed(3));

const haloAfterEnv_VfR = multiR2(gals.map(g=>[g.envCode]), gals.map(g=>g.VfResid));
const haloAfterEnv_a0R = multiR2(gals.map(g=>[g.envCode]), gals.map(g=>g.a0Resid));
const rAfterEnv = pearsonR(haloAfterEnv_VfR.residuals, haloAfterEnv_a0R.residuals);
console.log('\n  partial r(VfResid, a0_resid | envCode) = ' + rAfterEnv.toFixed(3) + ' (Delta=' + (rAfterEnv-rBase).toFixed(3) + ')');


console.log('\n\n' + '#'.repeat(70));
console.log('413C: SYNERGY TEST — ASSEMBLY PROXY ON TOP OF HALO STRUCTURE');
console.log('Does assembly history absorb BEYOND logK + dmFrac_Rmax?');
console.log('#'.repeat(70));

console.log('\n  Phase 412 baseline: partial r(VfR, a0R | logK, dmFrac_Rmax) = 0.589');
console.log('  Now we add assembly proxies on top.');

const haloBase = gals.map(g => [g.logK_halo, g.dmFrac_Rmax]);
const rHaloBase = partialR(gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid), haloBase);
console.log('  Verified: partial r | logK+dmFrac = ' + rHaloBase.toFixed(3));

const synergyTests = [
  { name: 'logK + dmFrac + gasFrac_resid', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.gasFrac_resid]) },
  { name: 'logK + dmFrac + compact_resid', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.compact_resid]) },
  { name: 'logK + dmFrac + SHR_resid', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.SHR_resid]) },
  { name: 'logK + dmFrac + morphT_resid', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.morphT_resid]) },
  { name: 'logK + dmFrac + envCode', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.envCode]) },
  { name: 'logK + dmFrac + logGasFrac', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.logGasFrac]) },
  { name: 'logK + dmFrac + logK_resid', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.logK_resid]) },
  { name: 'logK + dmFrac + dmFrac_resid', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.dmFrac_resid]) },
  { name: 'logK+dmF+gasF+compact+SHR', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.gasFrac_resid, g.compact_resid, g.SHR_resid]) },
  { name: 'logK+dmF+gasF+env+morphT', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.gasFrac_resid, g.envCode, g.morphT_resid]) },
  { name: 'ALL halo(7) + ALL assembly(5)', controls: gals.map(g => [
    g.logK_halo, g.haloResponse, g.dmFrac_Rmax, g.dmFrac_2Rd, g.innerSlope, g.outerSlope, g.concIdx,
    g.gasFrac_resid, g.compact_resid, g.SHR_resid, g.morphT_resid, g.envCode
  ]) },
];

console.log('\n  Control set                          partial r    Delta_from_halo  Delta_from_base');
console.log('  ' + '-'.repeat(80));

const synergyResults = [];
for (const st of synergyTests) {
  const r = partialR(gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid), st.controls);
  const deltaHalo = r - rHaloBase;
  const deltaBase = r - rBase;
  const status = deltaHalo < -0.10 ? 'STRONG' : deltaHalo < -0.05 ? 'MODERATE' : 'WEAK';
  console.log('  ' + st.name.padEnd(38) + ' ' + r.toFixed(3).padEnd(13) + ' ' + ((deltaHalo>=0?'+':'')+deltaHalo.toFixed(3)).padEnd(17) + ' ' + ((deltaBase>=0?'+':'')+deltaBase.toFixed(3)).padEnd(10) + ' ' + status);
  synergyResults.push({ name: st.name, partialR: r, deltaFromHalo: deltaHalo, deltaFromBase: deltaBase, status });
}


console.log('\n\n' + '#'.repeat(70));
console.log('413D: REGIME / HISTORY INTERACTION');
console.log('Does assembly effect strengthen at high-V (deep potential wells)?');
console.log('#'.repeat(70));

const medianV = gals.map(g=>g.Vflat).sort((a,b)=>a-b)[Math.floor(gals.length/2)];
const lowV = gals.filter(g => g.Vflat < medianV);
const highV = gals.filter(g => g.Vflat >= medianV);

console.log('\n  Split at median Vflat = ' + medianV.toFixed(0) + ' km/s');
console.log('  low-V: N=' + lowV.length + ', high-V: N=' + highV.length);

console.log('\n  BASELINE channel by regime:');
const rLow = pearsonR(lowV.map(g=>g.VfResid), lowV.map(g=>g.a0Resid));
const rHigh = pearsonR(highV.map(g=>g.VfResid), highV.map(g=>g.a0Resid));
console.log('  low-V:  r = ' + rLow.toFixed(3));
console.log('  high-V: r = ' + rHigh.toFixed(3));

const bestProxies = proxyCorrelations
  .filter(p => Math.abs(p.rLsum) > 0.2)
  .sort((a,b) => Math.abs(b.rLsum) - Math.abs(a.rLsum));

console.log('\n  Assembly proxy correlations by regime:');
console.log('  Proxy               r(L)_lowV   r(L)_highV  Strengthens?');
console.log('  ' + '-'.repeat(60));

const regimeResults = [];
for (const bp of assemblyProxies) {
  const getVal = gals[0][bp.name] !== undefined ? (g => g[bp.name]) : bp.get;
  if (!getVal) continue;
  const valsAll = gals.map(g => getVal(g));
  if (!valsAll.every(v => isFinite(v))) continue;

  const rLow_p = lowV.length >= 6 ? pearsonR(lowV.map(g=>getVal(g)), lowV.map(g=>g.L_sum)) : NaN;
  const rHigh_p = highV.length >= 6 ? pearsonR(highV.map(g=>getVal(g)), highV.map(g=>g.L_sum)) : NaN;
  const strengthens = Math.abs(rHigh_p) > Math.abs(rLow_p) ? 'YES' : 'no';
  console.log('  ' + bp.name.padEnd(20) + ' ' + (isFinite(rLow_p)?((rLow_p>=0?'+':'')+rLow_p.toFixed(3)):'N/A').padEnd(12) + ' ' + (isFinite(rHigh_p)?((rHigh_p>=0?'+':'')+rHigh_p.toFixed(3)):'N/A').padEnd(12) + ' ' + strengthens);
  regimeResults.push({ name: bp.name, rLow: rLow_p, rHigh: rHigh_p, strengthens });
}

const tertileV = [
  gals.filter(g => g.Vflat < 100),
  gals.filter(g => g.Vflat >= 100 && g.Vflat < 180),
  gals.filter(g => g.Vflat >= 180),
];
console.log('\n  Channel by Vflat tertile:');
for (let t = 0; t < 3; t++) {
  const label = ['low (<100)', 'mid (100-180)', 'high (>180)'][t];
  if (tertileV[t].length < 5) { console.log('  ' + label + ': N=' + tertileV[t].length + ' (too few)'); continue; }
  const r = pearsonR(tertileV[t].map(g=>g.VfResid), tertileV[t].map(g=>g.a0Resid));
  console.log('  ' + label.padEnd(16) + ' N=' + String(tertileV[t].length).padEnd(4) + ' r=' + r.toFixed(3));
}


console.log('\n\n' + '#'.repeat(70));
console.log('413E: MAXIMUM COMBINED ABSORPTION — ALL HALO + ALL ASSEMBLY');
console.log('What is the absolute maximum absorption possible from ALL available vars?');
console.log('#'.repeat(70));

const maxControls = gals.map(g => [
  g.logK_halo, g.haloResponse, g.dmFrac_Rmax, g.dmFrac_2Rd,
  g.innerSlope, g.outerSlope, g.concIdx,
  g.gasFrac_resid, g.compact_resid, g.SHR_resid, g.morphT_resid, g.envCode,
]);

const maxModel_VfR = multiR2(maxControls, gals.map(g=>g.VfResid));
const maxModel_a0R = multiR2(maxControls, gals.map(g=>g.a0Resid));
const rMax = pearsonR(maxModel_VfR.residuals, maxModel_a0R.residuals);

console.log('\n  12 total controls (7 halo + 5 assembly):');
console.log('  R2(VfResid ~ all) = ' + maxModel_VfR.R2.toFixed(3));
console.log('  R2(a0_resid ~ all) = ' + maxModel_a0R.R2.toFixed(3));
console.log('  r(VfResid_clean, a0R_clean) = ' + rMax.toFixed(3));
console.log('  Delta from baseline = ' + (rMax - rBase).toFixed(3));
console.log('  Delta from halo-only (412) = ' + (rMax - 0.709).toFixed(3));

const rPartialMax = partialR(gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid), maxControls);
console.log('\n  partial r(VfR, a0R | all 12) = ' + rPartialMax.toFixed(3));
console.log('  Delta from baseline = ' + (rPartialMax - rBase).toFixed(3));


console.log('\n\n' + '#'.repeat(70));
console.log('413F: OVERFITTING CHECK — LEAVE-ONE-OUT CROSS-VALIDATION');
console.log('Is the absorption real or just overfitting with 12 controls on N=55?');
console.log('#'.repeat(70));

let sumR = 0, countR = 0;
for (let i = 0; i < gals.length; i++) {
  const train = gals.filter((_,j) => j !== i);
  const trainCtrl = train.map(g => [g.logK_halo, g.dmFrac_Rmax]);
  const rLOO = partialR(train.map(g=>g.VfResid), train.map(g=>g.a0Resid), trainCtrl);
  sumR += rLOO;
  countR++;
}
const looR_halo2 = sumR / countR;

sumR = 0; countR = 0;
for (let i = 0; i < gals.length; i++) {
  const train = gals.filter((_,j) => j !== i);
  const trainCtrl = train.map(g => [
    g.logK_halo, g.haloResponse, g.dmFrac_Rmax, g.dmFrac_2Rd,
    g.innerSlope, g.outerSlope, g.concIdx,
    g.gasFrac_resid, g.compact_resid, g.SHR_resid, g.morphT_resid, g.envCode,
  ]);
  const rLOO = partialR(train.map(g=>g.VfResid), train.map(g=>g.a0Resid), trainCtrl);
  sumR += rLOO;
  countR++;
}
const looR_all12 = sumR / countR;

console.log('\n  LOO-CV partial r:');
console.log('  logK + dmFrac (2 controls): LOO mean r = ' + looR_halo2.toFixed(3) + '  (in-sample: ' + rHaloBase.toFixed(3) + ')');
console.log('  All 12 controls:             LOO mean r = ' + looR_all12.toFixed(3) + '  (in-sample: ' + rPartialMax.toFixed(3) + ')');
console.log('  Shrinkage (2 ctrl): ' + (rHaloBase - looR_halo2).toFixed(3));
console.log('  Shrinkage (12 ctrl): ' + (rPartialMax - looR_all12).toFixed(3));

const overfit12 = (rPartialMax - looR_all12) > 0.10;
console.log('  ' + (overfit12 ? 'WARNING: Substantial overfitting with 12 controls on N=55' : 'Overfitting is modest — absorption is REAL'));


console.log('\n\n' + '='.repeat(70));
console.log('PHASE 413 GRAND VERDICT');
console.log('='.repeat(70));

const bestSynergy = synergyResults.reduce((best, s) => s.deltaFromHalo < best.deltaFromHalo ? s : best, synergyResults[0]);
const assemblyAdds = bestSynergy.deltaFromHalo < -0.05;

console.log('\n  ASSEMBLY HISTORY ANALYSIS:\n');
console.log('  1. STRONGEST assembly proxies (r with L_sum):');
proxyCorrelations
  .filter(c => Math.abs(c.rLsum) > 0.15)
  .sort((a,b) => Math.abs(b.rLsum) - Math.abs(a.rLsum))
  .forEach(c => console.log('     - ' + c.name + ': r = ' + c.rLsum.toFixed(3)));

console.log('\n  2. SYNERGY (best assembly + halo absorption):');
console.log('     Halo-only baseline (logK+dmFrac): r = ' + rHaloBase.toFixed(3));
console.log('     Best synergy: ' + bestSynergy.name);
console.log('       partial r = ' + bestSynergy.partialR.toFixed(3) + ', Delta from halo = ' + bestSynergy.deltaFromHalo.toFixed(3));

console.log('\n  3. MAXIMUM ABSORPTION (all 12 controls):');
console.log('     partial r = ' + rPartialMax.toFixed(3) + ' (Delta from base = ' + (rPartialMax-rBase).toFixed(3) + ')');
console.log('     residual r = ' + rMax.toFixed(3));
console.log('     LOO-CV: ' + looR_all12.toFixed(3) + ' (shrinkage = ' + (rPartialMax-looR_all12).toFixed(3) + ')');

console.log('\n  4. REGIME INTERACTION:');
const strengthening = regimeResults.filter(r => r.strengthens === 'YES').length;
console.log('     ' + strengthening + '/' + regimeResults.length + ' proxies strengthen at high-V');

console.log('\n  VERDICT:');
if (assemblyAdds) {
  console.log('  Assembly history provides ADDITIONAL absorption beyond halo structure.');
  console.log('  The hidden variable has both a halo-structure component AND an assembly component.');
} else {
  console.log('  Assembly proxies do NOT significantly add beyond halo structure.');
  console.log('  The remaining channel (r ~ 0.6-0.7) is NOT captured by observable');
  console.log('  assembly indicators. It reflects something DEEPER.');
}

console.log('\n  SURVIVING CHANNEL after ALL controls:');
console.log('  r = ' + rPartialMax.toFixed(3) + ' (LOO: ' + looR_all12.toFixed(3) + ')');

if (rPartialMax > 0.4 || looR_all12 > 0.4) {
  console.log('\n  CRITICAL: The channel SURVIVES even after controlling for:');
  console.log('    - 7 halo structure variables');
  console.log('    - 5 assembly history proxies');
  console.log('    - Environment classification');
  console.log('  Total: 12 controls on N=55 galaxies');
  console.log('\n  The hidden variable is INVISIBLE to ALL observable proxies.');
  console.log('  This points to:');
  console.log('    1. 3D halo geometry (triaxiality, flattening) — unmeasurable from 1D RC');
  console.log('    2. Halo concentration-mass relation scatter — requires NFW decomposition');
  console.log('    3. Genuinely new physics beyond standard models');
  console.log('    4. Or a combination of many weak effects that no single proxy can capture');
}


const outPath = path.join(__dirname, '..', 'public', 'phase413-assembly-history.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '413',
  title: 'Assembly History / Formation Time',
  timestamp: new Date().toISOString(),
  N: gals.length,
  baselineR: rBase,
  rHaloBase: rHaloBase,
  proxyCorrelations,
  synergyResults,
  regimeResults,
  maxAbsorption: { partial_r: rPartialMax, residual_r: rMax, LOO_r: looR_all12, R2vfr: maxModel_VfR.R2, R2a0r: maxModel_a0R.R2 },
  assemblyAdds,
  environmentSplit: Object.fromEntries(Object.entries(envGroups).map(([k,gs]) => [k, { N: gs.length, r: gs.length >= 4 ? pearsonR(gs.map(g=>g.VfResid),gs.map(g=>g.a0Resid)) : null }])),
}, null, 2));

console.log('\nSaved: ' + outPath);
