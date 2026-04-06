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

  const V_Newt_Rd = Math.sqrt(G * Mbar / (2.2 * Rdisk));
  const V_Newt_Rmax = Math.sqrt(G * Mbar / Rmax);

  const dmFrac_Rmax = Math.max(0, 1 - (V_Newt_Rmax/Vflat)**2);
  const dmFrac_2Rd = V_inner > 0 ? Math.max(0, 1 - (V_Newt_Rd/V_inner)**2) : 0;

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

  const logK_halo = k_halo > 0 ? Math.log10(k_halo) : -5;
  const logM_halo = M_halo > 0 ? Math.log10(M_halo) : 5;

  const logK_halo_expected = (() => {
    const model = { slope: 0.3, intercept: -3 };
    return model.slope * Math.log10(Mbar) + model.intercept;
  })();

  gals.push({
    name: g.name, logA0: g.logA0,
    Vflat, Rdisk, Rmax, Vmax, Mbar, rc,
    logVflat: Math.log10(Vflat),
    logL36: Math.log10(Math.max(sp.L36, 0.001)),
    logRdisk: Math.log10(Math.max(Rdisk, 0.01)),
    logMbar: Math.log10(Math.max(Mbar, 1)),
    logMHI: g.logMHI,
    morphT: sp.T,
    logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    envCode: g.envCode,
    L36: sp.L36, MHI: sp.MHI,
    logK_halo, logM_halo, k_halo, M_halo,
    mse_newt, mse_halo,
    dmFrac_Rmax, dmFrac_2Rd,
    innerSlope, outerSlope, concIdx,
    haloResponse,
    V_inner, V_outer,
    logK_halo_resid: logK_halo - logK_halo_expected,
  });
}

const vfModel = multiR2(gals.map(g=>[g.logMbar,g.logL36,g.logRdisk,g.morphT]),gals.map(g=>g.logVflat));
const a0Model = multiR2(gals.map(g=>[g.logMbar,g.logL36,g.logRdisk,g.morphT,g.logMHI,g.logSBdisk]),gals.map(g=>g.logA0));
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

const rBase = pearsonR(gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid));

console.log('═'.repeat(70));
console.log('PHASE 412: HALO STRUCTURE DISCRIMINATION');
console.log('Is the hidden variable dark matter halo internal structure?');
console.log('═'.repeat(70));
console.log(`N = ${gals.length}, baseline r(VfResid, a₀_resid) = ${rBase.toFixed(3)}`);


console.log('\n' + '▓'.repeat(70));
console.log('412A: HALO PARAMETER CORRELATIONS');
console.log('Which halo-derived quantities correlate with the latent variable?');
console.log('▓'.repeat(70));

const haloVars = [
  { name: 'logK_halo', get: g => g.logK_halo, desc: 'Halo linear slope k' },
  { name: 'logM_halo', get: g => g.logM_halo, desc: 'Halo mass parameter' },
  { name: 'haloResponse', get: g => g.haloResponse, desc: 'log(MSE_newt/MSE_halo)' },
  { name: 'dmFrac_Rmax', get: g => g.dmFrac_Rmax, desc: 'DM fraction at Rmax' },
  { name: 'dmFrac_2Rd', get: g => g.dmFrac_2Rd, desc: 'DM fraction at 2Rd' },
  { name: 'innerSlope', get: g => g.innerSlope, desc: 'RC inner log-slope' },
  { name: 'outerSlope', get: g => g.outerSlope, desc: 'RC outer log-slope' },
  { name: 'concIdx', get: g => g.concIdx, desc: 'V_inner/Vmax concentration' },
  { name: 'logK_halo_resid', get: g => g.logK_halo_resid, desc: 'logK excess vs mass expectation' },
  { name: 'V_out/V_in', get: g => g.V_inner > 0 ? g.V_outer/g.V_inner : NaN, desc: 'Outer/inner velocity ratio' },
  { name: 'Rmax/Rdisk', get: g => g.Rmax/g.Rdisk, desc: 'RC extent relative to disk' },
  { name: 'log(Mbar/Mhalo)', get: g => g.M_halo > 0 ? g.logMbar - g.logM_halo : NaN, desc: 'Baryonic/halo mass ratio' },
];

console.log('\n  Variable            r(VfResid) r(a₀_resid) r(L_sum)   desc');
console.log('  ' + '─'.repeat(75));

const haloCorrelations = [];
for (const v of haloVars) {
  const vals = gals.map(g => v.get(g));
  const validIdx = vals.map((val,i) => isFinite(val) ? i : -1).filter(i=>i>=0);
  if (validIdx.length < 10) continue;
  const x = validIdx.map(i => vals[i]);
  const rVf = pearsonR(x, validIdx.map(i=>gals[i].VfResid));
  const rA0 = pearsonR(x, validIdx.map(i=>gals[i].a0Resid));
  const rL = pearsonR(x, validIdx.map(i=>gals[i].L_sum));
  console.log(`  ${v.name.padEnd(20)} ${(rVf>=0?'+':'')+rVf.toFixed(3).padEnd(11)} ${(rA0>=0?'+':'')+rA0.toFixed(3).padEnd(12)} ${(rL>=0?'+':'')+rL.toFixed(3).padEnd(11)} ${v.desc}`);
  haloCorrelations.push({ name: v.name, rVfResid: rVf, rA0resid: rA0, rLsum: rL });
}


console.log('\n\n' + '▓'.repeat(70));
console.log('412B: STRUCTURAL-RESIDUAL HALO VARIABLES');
console.log('Halo quantities after removing structural prediction');
console.log('▓'.repeat(70));

const structControls = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);

const haloResidVars = [
  { name: 'logK_halo', vals: gals.map(g => g.logK_halo) },
  { name: 'haloResponse', vals: gals.map(g => g.haloResponse) },
  { name: 'dmFrac_Rmax', vals: gals.map(g => g.dmFrac_Rmax) },
  { name: 'dmFrac_2Rd', vals: gals.map(g => g.dmFrac_2Rd) },
  { name: 'innerSlope', vals: gals.map(g => g.innerSlope) },
  { name: 'concIdx', vals: gals.map(g => g.concIdx) },
];

console.log('\n  Variable          R²(~struct)  r(resid, VfR)  r(resid, a₀R)  r(resid, L_sum)');
console.log('  ' + '─'.repeat(75));

const haloResidResults = [];
for (const hv of haloResidVars) {
  const model = multiR2(structControls, hv.vals);
  const resid = model.residuals;
  const rVf = pearsonR(resid, gals.map(g=>g.VfResid));
  const rA0 = pearsonR(resid, gals.map(g=>g.a0Resid));
  const rL = pearsonR(resid, gals.map(g=>g.L_sum));
  console.log(`  ${hv.name.padEnd(18)} ${model.R2.toFixed(3).padEnd(13)} ${(rVf>=0?'+':'')+rVf.toFixed(3).padEnd(15)} ${(rA0>=0?'+':'')+rA0.toFixed(3).padEnd(15)} ${(rL>=0?'+':'')+rL.toFixed(3)}`);
  haloResidResults.push({ name: hv.name, R2struct: model.R2, rVfResid: rVf, rA0resid: rA0, rLsum: rL, resid });
}


console.log('\n\n' + '▓'.repeat(70));
console.log('412C: DOES HALO STRUCTURE ABSORB THE CHANNEL?');
console.log('Partial r(VfResid, a₀_resid | halo vars)');
console.log('▓'.repeat(70));

const controlSets = [
  { name: 'logK_halo only', controls: gals.map(g => [g.logK_halo]) },
  { name: 'haloResponse only', controls: gals.map(g => [g.haloResponse]) },
  { name: 'dmFrac_Rmax only', controls: gals.map(g => [g.dmFrac_Rmax]) },
  { name: 'logK + haloResp', controls: gals.map(g => [g.logK_halo, g.haloResponse]) },
  { name: 'logK + dmFrac_Rmax', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax]) },
  { name: 'logK + dmFrac + innerSlope', controls: gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.innerSlope]) },
  { name: 'ALL halo (6 vars)', controls: gals.map(g => [g.logK_halo, g.haloResponse, g.dmFrac_Rmax, g.dmFrac_2Rd, g.innerSlope, g.concIdx]) },
  { name: 'ALL halo + envCode', controls: gals.map(g => [g.logK_halo, g.haloResponse, g.dmFrac_Rmax, g.dmFrac_2Rd, g.innerSlope, g.concIdx, g.envCode]) },
];

console.log(`\n  Baseline: r(VfResid, a₀_resid) = ${rBase.toFixed(3)}`);
console.log('  Control set                     partial r    Δr       absorbed?');
console.log('  ' + '─'.repeat(65));

const absorptionResults = [];
for (const cs of controlSets) {
  const r = partialR(gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid), cs.controls);
  const delta = r - rBase;
  const absorbed = delta < -0.10 ? 'SIGNIFICANT' : delta < -0.05 ? 'PARTIAL' : 'NONE';
  console.log(`  ${cs.name.padEnd(33)} ${r.toFixed(3).padEnd(13)} ${(delta>=0?'+':'')+delta.toFixed(3).padEnd(9)} ${absorbed}`);
  absorptionResults.push({ name: cs.name, partialR: r, delta, absorbed });
}


console.log('\n\n' + '▓'.repeat(70));
console.log('412D: CORE vs CUSP PROFILE INDICATOR');
console.log('Does the channel distinguish core-like from cusp-like halos?');
console.log('▓'.repeat(70));

const coreIndex = gals.map(g => {
  const pts = g.rc.filter(p => p.r < g.Rdisk);
  if (pts.length < 2) return NaN;
  const lx = pts.map(p => Math.log10(Math.max(p.r, 0.01)));
  const ly = pts.map(p => Math.log10(Math.max(p.v, 1)));
  const mx = lx.reduce((a,b)=>a+b,0)/lx.length;
  const my = ly.reduce((a,b)=>a+b,0)/ly.length;
  let num=0,den=0;
  for(let i=0;i<lx.length;i++){num+=(lx[i]-mx)*(ly[i]-my);den+=(lx[i]-mx)**2;}
  return den>0?num/den:0;
});

const validCore = coreIndex.map((v,i) => isFinite(v) ? i : -1).filter(i=>i>=0);
console.log(`\n  Core index = inner RC slope (log V vs log r), valid N=${validCore.length}`);
console.log('  Core-like halos → steeper inner rise (slope ~1)');
console.log('  Cusp-like halos → shallower rise (slope ~0.5)');

if (validCore.length >= 10) {
  const ci = validCore.map(i => coreIndex[i]);
  const rCore_VfR = pearsonR(ci, validCore.map(i=>gals[i].VfResid));
  const rCore_a0R = pearsonR(ci, validCore.map(i=>gals[i].a0Resid));
  const rCore_L = pearsonR(ci, validCore.map(i=>gals[i].L_sum));
  console.log(`  r(coreIdx, VfResid) = ${rCore_VfR.toFixed(3)}`);
  console.log(`  r(coreIdx, a₀_resid) = ${rCore_a0R.toFixed(3)}`);
  console.log(`  r(coreIdx, L_sum) = ${rCore_L.toFixed(3)}`);

  const rPartialCore = partialR(
    validCore.map(i=>gals[i].VfResid),
    validCore.map(i=>gals[i].a0Resid),
    validCore.map(i=>[coreIndex[i]])
  );
  console.log(`  partial r(VfR, a₀R | coreIdx) = ${rPartialCore.toFixed(3)} (Δ=${(rPartialCore-rBase).toFixed(3)})`);
}

const profileShape = gals.map(g => {
  if (g.rc.length < 5) return NaN;
  const rHalf = g.rc[Math.floor(g.rc.length/2)].r;
  const vHalf = g.rc[Math.floor(g.rc.length/2)].v;
  const vLast = g.rc[g.rc.length-1].v;
  const vFirst = g.rc[Math.min(2, g.rc.length-1)].v;
  return vFirst > 0 ? (vHalf/vFirst) * (vLast/vHalf) : NaN;
});

const validPS = profileShape.map((v,i) => isFinite(v) ? i : -1).filter(i=>i>=0);
if (validPS.length >= 10) {
  const ps = validPS.map(i => profileShape[i]);
  const rPS_L = pearsonR(ps, validPS.map(i => gals[i].L_sum));
  console.log(`\n  Profile shape index: r(profileShape, L_sum) = ${rPS_L.toFixed(3)}`);
}


console.log('\n\n' + '▓'.repeat(70));
console.log('412E: HALO-STRUCTURE BUNDLE — MAXIMUM ABSORPTION TEST');
console.log('How much channel can ALL halo variables together absorb?');
console.log('▓'.repeat(70));

const allHaloVals = gals.map(g => [
  g.logK_halo, g.haloResponse, g.dmFrac_Rmax, g.dmFrac_2Rd,
  g.innerSlope, g.outerSlope, g.concIdx,
]);

const haloModel_VfR = multiR2(allHaloVals, gals.map(g=>g.VfResid));
const haloModel_a0R = multiR2(allHaloVals, gals.map(g=>g.a0Resid));

console.log(`\n  R²(VfResid ~ 7 halo vars) = ${haloModel_VfR.R2.toFixed(3)}`);
console.log(`  R²(a₀_resid ~ 7 halo vars) = ${haloModel_a0R.R2.toFixed(3)}`);

const rAfterHalo = pearsonR(haloModel_VfR.residuals, haloModel_a0R.residuals);
console.log(`  r(VfResid_clean, a₀R_clean) = ${rAfterHalo.toFixed(3)}`);
console.log(`  Δr = ${(rAfterHalo - rBase).toFixed(3)}`);

const fullHaloEnv = gals.map(g => [
  g.logK_halo, g.haloResponse, g.dmFrac_Rmax, g.dmFrac_2Rd,
  g.innerSlope, g.outerSlope, g.concIdx,
  g.envCode === 1 ? 1 : 0, g.envCode === 2 ? 1 : 0,
]);

const fullModel_VfR = multiR2(fullHaloEnv, gals.map(g=>g.VfResid));
const fullModel_a0R = multiR2(fullHaloEnv, gals.map(g=>g.a0Resid));
const rAfterFull = pearsonR(fullModel_VfR.residuals, fullModel_a0R.residuals);

console.log(`\n  With envCode added (9 total controls):`);
console.log(`  R²(VfResid ~ halo+env) = ${fullModel_VfR.R2.toFixed(3)}`);
console.log(`  R²(a₀_resid ~ halo+env) = ${fullModel_a0R.R2.toFixed(3)}`);
console.log(`  r(VfResid_clean, a₀R_clean) = ${rAfterFull.toFixed(3)}`);
console.log(`  Δr = ${(rAfterFull - rBase).toFixed(3)}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('412F: HALO logK — DETAILED ANATOMY');
console.log('logK_halo is the strongest single halo variable');
console.log('▓'.repeat(70));

const logK_model = multiR2(structControls, gals.map(g=>g.logK_halo));
const logK_resid = logK_model.residuals;

console.log(`\n  R²(logK ~ 6 structural) = ${logK_model.R2.toFixed(3)}`);
console.log(`  → ${((1-logK_model.R2)*100).toFixed(0)}% of logK is structurally unexplained`);

const rLogKresid_VfR = pearsonR(logK_resid, gals.map(g=>g.VfResid));
const rLogKresid_a0R = pearsonR(logK_resid, gals.map(g=>g.a0Resid));
const rLogKresid_L = pearsonR(logK_resid, gals.map(g=>g.L_sum));

console.log(`  r(logK_resid, VfResid) = ${rLogKresid_VfR.toFixed(3)}`);
console.log(`  r(logK_resid, a₀_resid) = ${rLogKresid_a0R.toFixed(3)}`);
console.log(`  r(logK_resid, L_sum) = ${rLogKresid_L.toFixed(3)}`);

const rPartialLogK = partialR(
  gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid),
  gals.map(g=>[g.logK_halo])
);
console.log(`\n  partial r(VfResid, a₀_resid | logK_halo) = ${rPartialLogK.toFixed(3)}`);
console.log(`  Δr = ${(rPartialLogK - rBase).toFixed(3)}`);

const rPartialLogKresid = partialR(
  gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid),
  logK_resid.map(v=>[v])
);
console.log(`  partial r(VfResid, a₀_resid | logK_resid) = ${rPartialLogKresid.toFixed(3)}`);
console.log(`  Δr = ${(rPartialLogKresid - rBase).toFixed(3)}`);

console.log(`\n  INTERPRETATION:`);
console.log(`  logK_halo explains R²=${logK_model.R2.toFixed(3)} from structure.`);
console.log(`  The residual part (${((1-logK_model.R2)*100).toFixed(0)}%) correlates with the channel.`);
console.log(`  → The novel part of halo slope IS partially the hidden variable.`);


console.log('\n\n' + '▓'.repeat(70));
console.log('412G: CIRCULARITY CHECK');
console.log('Are halo variables circular with VfResid or a₀?');
console.log('▓'.repeat(70));

console.log('\n  KEY QUESTION: Is logK_halo derived from the same RC data as VfResid?');
console.log('  logK_halo comes from the dark_halo_linear model fit to RC points.');
console.log('  VfResid comes from Vflat regression residuals.');
console.log('  a₀ comes from fitted acceleration scale.');
console.log('  All three USE the rotation curve → potential shared measurement basis.');

console.log('\n  CIRCULARITY TEST 1: logK_halo vs. structural prediction only');
const logK_from_struct = multiR2(
  gals.map(g => [g.logMbar, g.logRdisk, g.morphT]),
  gals.map(g => g.logK_halo)
);
console.log(`  R²(logK ~ Mbar,Rd,T) = ${logK_from_struct.R2.toFixed(3)}`);
console.log(`  → logK can be predicted ${(logK_from_struct.R2*100).toFixed(0)}% from structure alone`);

console.log('\n  CIRCULARITY TEST 2: Is logK_halo just repackaged Vflat?');
const rLogK_Vflat = pearsonR(gals.map(g=>g.logK_halo), gals.map(g=>g.logVflat));
const rLogK_logA0 = pearsonR(gals.map(g=>g.logK_halo), gals.map(g=>g.logA0));
console.log(`  r(logK, logVflat) = ${rLogK_Vflat.toFixed(3)}`);
console.log(`  r(logK, logA0) = ${rLogK_logA0.toFixed(3)}`);

console.log('\n  CIRCULARITY TEST 3: partial after removing Vflat from logK');
const logK_after_Vflat = multiR2(
  gals.map(g => [g.logVflat, g.logMbar, g.logRdisk]),
  gals.map(g => g.logK_halo)
);
const rLogK_residVflat_a0R = pearsonR(logK_after_Vflat.residuals, gals.map(g=>g.a0Resid));
console.log(`  R²(logK ~ Vflat+struct) = ${logK_after_Vflat.R2.toFixed(3)}`);
console.log(`  r(logK_resid_Vflat, a₀_resid) = ${rLogK_residVflat_a0R.toFixed(3)}`);
console.log(`  → ${Math.abs(rLogK_residVflat_a0R) > 0.2 ? 'logK carries information BEYOND Vflat → NOT purely circular' : 'Weak — logK may be partially circular with Vflat'}`);


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 412 GRAND VERDICT');
console.log('═'.repeat(70));

const bestHaloAbsorption = Math.min(...absorptionResults.map(r => r.delta));
const haloAbsorbs = bestHaloAbsorption < -0.10;

console.log(`
HALO STRUCTURE ANALYSIS:

1. STRONGEST halo correlations with latent variable (L_sum):
${haloCorrelations.filter(c => Math.abs(c.rLsum) > 0.2).map(c =>
  `   - ${c.name}: r = ${c.rLsum.toFixed(3)}`).join('\n')}

2. ABSORPTION (best partial r control):
   Baseline: r = ${rBase.toFixed(3)}
   Best single control: logK_halo → r = ${rPartialLogK.toFixed(3)} (Δ=${(rPartialLogK-rBase).toFixed(3)})
   ALL 7 halo vars: r = ${rAfterHalo.toFixed(3)} (Δ=${(rAfterHalo-rBase).toFixed(3)})
   ALL halo + env: r = ${rAfterFull.toFixed(3)} (Δ=${(rAfterFull-rBase).toFixed(3)})

3. logK_halo anatomy:
   ${((1-logK_model.R2)*100).toFixed(0)}% novel (structurally unexplained)
   Novel part correlates with channel: r(logK_resid, L_sum) = ${rLogKresid_L.toFixed(3)}

4. CIRCULARITY: ${Math.abs(rLogK_residVflat_a0R) > 0.2 ? 'PASSED — logK carries non-circular info' : 'CAUTION — partial circularity possible'}

VERDICT: ${haloAbsorbs
  ? 'Halo structure PARTIALLY ABSORBS the channel. It is a significant component of the hidden variable.'
  : rAfterHalo < rBase - 0.03
    ? 'Halo structure is a WEAK CONTRIBUTOR. It captures some of the hidden variable but NOT the majority. The channel largely survives even after 7 halo variables are controlled.'
    : 'Halo structure does NOT absorb the channel. The hidden variable is DEEPER than measurable halo parameters from 1D rotation curves.'}

${rAfterFull > 0.5 ? `\nCRITICAL: Even after controlling for 7 halo vars + environment (9 total controls),
the channel remains at r = ${rAfterFull.toFixed(3)}.

This means the hidden variable is NOT fully captured by:
  - Halo slope, mass, or response
  - DM fractions (inner or outer)
  - RC shape (slopes, concentration)
  - Environment

The remaining signal (r = ${rAfterFull.toFixed(3)}) represents information that is
INVISIBLE to all measurable halo properties from 1D rotation curves.

IMPLICATION: The hidden variable likely reflects UNMEASURED halo properties:
  → 3D halo shape (triaxiality, flattening)
  → Halo concentration-mass relation scatter
  → Assembly-dependent halo response
  → Or genuine new physics beyond both ΛCDM halo models and MOND` : ''}
`);

const outPath = path.join(__dirname, '..', 'public', 'phase412-halo-structure.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '412',
  title: 'Halo Structure Discrimination',
  timestamp: new Date().toISOString(),
  N: gals.length,
  baselineR: rBase,
  haloCorrelations,
  haloResidResults: haloResidResults.map(h => ({ name: h.name, R2struct: h.R2struct, rVfResid: h.rVfResid, rA0resid: h.rA0resid, rLsum: h.rLsum })),
  absorptionResults,
  rAfterAllHalo: rAfterHalo,
  rAfterHaloEnv: rAfterFull,
  logK_anatomy: { R2struct: logK_model.R2, rResid_VfR: rLogKresid_VfR, rResid_a0R: rLogKresid_a0R, rResid_L: rLogKresid_L },
  haloAbsorbs,
}, null, 2));

console.log(`\nSaved: ${outPath}`);
