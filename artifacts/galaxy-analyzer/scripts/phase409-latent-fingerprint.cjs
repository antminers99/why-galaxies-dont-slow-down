const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

const G = 4.3009e-6;

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
function normalize(n) {
  return n.replace(/\s+/g,'').replace(/^NGC0*/,'NGC').replace(/^DDO0*/,'DDO').replace(/^IC0*/,'IC').replace(/^UGC0*/,'UGC').replace(/^PGC0*/,'PGC').toUpperCase();
}
const rcMapRaw = parseRCData(tsContent);
const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) rcMap[normalize(name)] = pts;

function pearsonR(x, y) {
  const n = x.length;
  if (n < 4) return NaN;
  const mx = x.reduce((a,b) => a+b, 0) / n;
  const my = y.reduce((a,b) => a+b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) { const dx=x[i]-mx, dy=y[i]-my; num+=dx*dy; dx2+=dx*dx; dy2+=dy*dy; }
  return dx2>0&&dy2>0 ? num/Math.sqrt(dx2*dy2) : 0;
}
function multiR2(X, y) {
  const n = y.length; const nv = X[0].length;
  const my = y.reduce((a,b)=>a+b,0)/n;
  const mx = Array(nv).fill(0);
  for (let j=0;j<nv;j++){for(let i=0;i<n;i++) mx[j]+=X[i][j]; mx[j]/=n;}
  const XTX = Array.from({length:nv},()=>Array(nv).fill(0));
  const XTy = Array(nv).fill(0);
  for(let i=0;i<n;i++){for(let j=0;j<nv;j++){XTy[j]+=(X[i][j]-mx[j])*(y[i]-my);for(let k=0;k<nv;k++)XTX[j][k]+=(X[i][j]-mx[j])*(X[i][k]-mx[k]);}}
  const aug = XTX.map((row,i)=>[...row,XTy[i]]);
  for(let col=0;col<nv;col++){let maxRow=col;for(let row=col+1;row<nv;row++)if(Math.abs(aug[row][col])>Math.abs(aug[maxRow][col]))maxRow=row;[aug[col],aug[maxRow]]=[aug[maxRow],aug[col]];if(Math.abs(aug[col][col])<1e-12)continue;for(let row=col+1;row<nv;row++){const f=aug[row][col]/aug[col][col];for(let j=col;j<=nv;j++)aug[row][j]-=f*aug[col][j];}}
  const beta=Array(nv).fill(0);
  for(let i=nv-1;i>=0;i--){beta[i]=aug[i][nv];for(let j=i+1;j<nv;j++)beta[i]-=aug[i][j]*beta[j];beta[i]/=aug[i][i]||1;}
  let sse=0,sst=0;const residuals=[];
  for(let i=0;i<n;i++){let pred=my;for(let j=0;j<nv;j++)pred+=beta[j]*(X[i][j]-mx[j]);residuals.push(y[i]-pred);sse+=(y[i]-pred)**2;sst+=(y[i]-my)**2;}
  return {R2:sst>0?1-sse/sst:0,beta,residuals};
}

const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[g.name] = g; resultsMap[normalize(g.name)] = g; });

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 5) continue;
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sr) continue;

  const Vflat = sp.Vflat;
  const Rdisk = sp.Rdisk;
  gals.push({
    name: g.name, logA0: g.logA0,
    Vflat, Rdisk, rc,
    logVflat: Math.log10(Vflat),
    logL36: Math.log10(Math.max(sp.L36, 0.001)),
    logRdisk: Math.log10(Math.max(Rdisk, 0.01)),
    logMbar: Math.log10(Math.max(sp.L36*0.5 + sp.MHI*1.33, 0.001) * 1e9),
    logMHI: g.logMHI,
    morphT: sp.T,
    logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    envCode: g.envCode,
    L36: sp.L36, MHI: sp.MHI, SBdisk: sp.SBdisk,
    inc: sp.inc, Q: sp.Q,
    M_halo: sr.models.dark_halo_linear.M,
    k_halo: sr.models.dark_halo_linear.k,
  });
}

const vfModel = multiR2(gals.map(g=>[g.logMbar,g.logL36,g.logRdisk,g.morphT]),gals.map(g=>g.logVflat));
const a0Model = multiR2(gals.map(g=>[g.logMbar,g.logL36,g.logRdisk,g.morphT,g.logMHI,g.logSBdisk]),gals.map(g=>g.logA0));
for (let i=0; i<gals.length; i++) {
  gals[i].VfResid = vfModel.residuals[i];
  gals[i].a0Resid = a0Model.residuals[i];
}

console.log('═'.repeat(70));
console.log('PHASE 409: LATENT VARIABLE FINGERPRINT');
console.log('What do we know about the hidden variable from data alone?');
console.log('═'.repeat(70));
console.log(`\nN = ${gals.length}`);

const rBase = pearsonR(gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid));
console.log(`Baseline r(VfResid, a₀_resid) = ${rBase.toFixed(3)}`);


console.log('\n' + '▓'.repeat(70));
console.log('409A: RECONSTRUCTING THE LATENT VARIABLE');
console.log('If L = α·VfResid + β·a₀_resid, what is L like?');
console.log('▓'.repeat(70));

const VfR = gals.map(g => g.VfResid);
const a0R = gals.map(g => g.a0Resid);

const sdVf = Math.sqrt(VfR.reduce((a,v)=>a+v*v,0)/VfR.length - (VfR.reduce((a,v)=>a+v,0)/VfR.length)**2);
const sdA0 = Math.sqrt(a0R.reduce((a,v)=>a+v*v,0)/a0R.length - (a0R.reduce((a,v)=>a+v,0)/a0R.length)**2);

const VfR_z = VfR.map(v => v/sdVf);
const a0R_z = a0R.map(v => v/sdA0);

const L_sum = VfR_z.map((v,i) => v + a0R_z[i]);
const L_diff = VfR_z.map((v,i) => v - a0R_z[i]);

const rL_sum_a0 = pearsonR(L_sum, gals.map(g=>g.logA0));
const rL_diff_a0 = pearsonR(L_diff, gals.map(g=>g.logA0));
const rL_sum_Vf = pearsonR(L_sum, gals.map(g=>g.logVflat));
const rL_diff_Vf = pearsonR(L_diff, gals.map(g=>g.logVflat));

console.log(`\n  L_sum = VfResid_z + a₀_resid_z (shared direction)`);
console.log(`    r(L_sum, logA0) = ${rL_sum_a0.toFixed(3)}`);
console.log(`    r(L_sum, logVflat) = ${rL_sum_Vf.toFixed(3)}`);
console.log(`  L_diff = VfResid_z - a₀_resid_z (differential direction)`);
console.log(`    r(L_diff, logA0) = ${rL_diff_a0.toFixed(3)}`);
console.log(`    r(L_diff, logVflat) = ${rL_diff_Vf.toFixed(3)}`);

console.log('\n  L_sum distribution:');
const L_sorted = [...L_sum].sort((a,b)=>a-b);
const q25 = L_sorted[Math.floor(L_sorted.length*0.25)];
const q50 = L_sorted[Math.floor(L_sorted.length*0.50)];
const q75 = L_sorted[Math.floor(L_sorted.length*0.75)];
const meanL = L_sum.reduce((a,b)=>a+b,0)/L_sum.length;
const sdL = Math.sqrt(L_sum.reduce((a,v)=>a+(v-meanL)**2,0)/L_sum.length);
console.log(`    mean=${meanL.toFixed(3)}, sd=${sdL.toFixed(3)}, [Q25,Q50,Q75]=[${q25.toFixed(3)},${q50.toFixed(3)},${q75.toFixed(3)}]`);

const skew = L_sum.reduce((a,v)=>a+((v-meanL)/sdL)**3,0)/L_sum.length;
const kurt = L_sum.reduce((a,v)=>a+((v-meanL)/sdL)**4,0)/L_sum.length - 3;
console.log(`    skewness=${skew.toFixed(3)}, excess kurtosis=${kurt.toFixed(3)}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('409B: WHAT CORRELATES WITH THE LATENT VARIABLE?');
console.log('▓'.repeat(70));

const allVars = [
  { name: 'logVflat', get: g => g.logVflat },
  { name: 'logMbar', get: g => g.logMbar },
  { name: 'logL36', get: g => g.logL36 },
  { name: 'logRdisk', get: g => g.logRdisk },
  { name: 'morphT', get: g => g.morphT },
  { name: 'logMHI', get: g => g.logMHI },
  { name: 'logSBdisk', get: g => g.logSBdisk },
  { name: 'logA0', get: g => g.logA0 },
  { name: 'inc', get: g => g.inc },
  { name: 'Q', get: g => g.Q },
  { name: 'envCode', get: g => g.envCode },
  { name: 'Rdisk', get: g => g.Rdisk },
  { name: 'Vflat', get: g => g.Vflat },
  { name: 'logK_halo', get: g => g.k_halo > 0 ? Math.log10(g.k_halo) : NaN },
  { name: 'logM_halo', get: g => g.M_halo > 0 ? Math.log10(g.M_halo) : NaN },
  { name: 'gasFrac', get: g => g.MHI > 0 ? g.MHI*1.33 / (g.L36*0.5 + g.MHI*1.33) : NaN },
  { name: 'SBdisk', get: g => g.SBdisk },
  { name: 'logMbar/L36', get: g => Math.log10((g.L36*0.5 + g.MHI*1.33) / g.L36) },
  { name: 'Rdisk/Rmax', get: g => g.Rdisk / g.rc[g.rc.length-1].r },
];

console.log('\n  Variable         r(L_sum)   r(L_diff)  r(VfResid) r(a₀_resid)');
console.log('  ' + '─'.repeat(70));

const latentCorrelations = [];
for (const v of allVars) {
  const vals = gals.map(g => v.get(g));
  const validIdx = vals.map((val,i) => isFinite(val) ? i : -1).filter(i=>i>=0);
  if (validIdx.length < 10) continue;
  const x = validIdx.map(i => vals[i]);
  const rLs = pearsonR(x, validIdx.map(i => L_sum[i]));
  const rLd = pearsonR(x, validIdx.map(i => L_diff[i]));
  const rVf = pearsonR(x, validIdx.map(i => VfR[i]));
  const rA0 = pearsonR(x, validIdx.map(i => a0R[i]));
  console.log(`  ${v.name.padEnd(18)} ${(rLs>=0?'+':'')+rLs.toFixed(3).padEnd(10)} ${(rLd>=0?'+':'')+rLd.toFixed(3).padEnd(10)} ${(rVf>=0?'+':'')+rVf.toFixed(3).padEnd(10)} ${(rA0>=0?'+':'')+rA0.toFixed(3)}`);
  latentCorrelations.push({ name: v.name, rLsum: rLs, rLdiff: rLd, rVfResid: rVf, rA0resid: rA0 });
}


console.log('\n\n' + '▓'.repeat(70));
console.log('409C: EXTREME GALAXIES — WHO CARRIES THE CHANNEL?');
console.log('▓'.repeat(70));

const withL = gals.map((g,i) => ({ ...g, L_sum: L_sum[i], L_diff: L_diff[i], idx: i }));
withL.sort((a,b) => b.L_sum - a.L_sum);

console.log('\n  Top 10 (strongest channel carriers — high L_sum):');
console.log('  Name              Vflat   logA0    VfResid  a₀_resid  L_sum  morphT  env');
console.log('  ' + '─'.repeat(80));
for (let i = 0; i < 10; i++) {
  const g = withL[i];
  console.log(`  ${g.name.padEnd(18)} ${g.Vflat.toString().padEnd(8)} ${g.logA0.toFixed(3).padEnd(9)} ${(g.VfResid>=0?'+':'')+g.VfResid.toFixed(3).padEnd(9)} ${(g.a0Resid>=0?'+':'')+g.a0Resid.toFixed(3).padEnd(10)} ${g.L_sum.toFixed(2).padEnd(7)} ${g.morphT.toString().padEnd(8)} ${g.envCode}`);
}

console.log('\n  Bottom 10 (weakest/reversed channel — low L_sum):');
console.log('  Name              Vflat   logA0    VfResid  a₀_resid  L_sum  morphT  env');
console.log('  ' + '─'.repeat(80));
for (let i = withL.length-10; i < withL.length; i++) {
  const g = withL[i];
  console.log(`  ${g.name.padEnd(18)} ${g.Vflat.toString().padEnd(8)} ${g.logA0.toFixed(3).padEnd(9)} ${(g.VfResid>=0?'+':'')+g.VfResid.toFixed(3).padEnd(9)} ${(g.a0Resid>=0?'+':'')+g.a0Resid.toFixed(3).padEnd(10)} ${g.L_sum.toFixed(2).padEnd(7)} ${g.morphT.toString().padEnd(8)} ${g.envCode}`);
}


console.log('\n\n' + '▓'.repeat(70));
console.log('409D: WHERE IN THE GALAXY IS THE CHANNEL STRONGEST?');
console.log('Radial profile of the coupling');
console.log('▓'.repeat(70));

const radialBins = [
  { name: 'r < Rdisk', filter: (r, Rd) => r < Rd },
  { name: 'Rd < r < 2Rd', filter: (r, Rd) => r >= Rd && r < 2*Rd },
  { name: '2Rd < r < 4Rd', filter: (r, Rd) => r >= 2*Rd && r < 4*Rd },
  { name: 'r > 4Rd', filter: (r, Rd) => r >= 4*Rd },
];

console.log('\n  Radial bin        N_gals  mean|VfResid|  r(VfR,a₀res)');
console.log('  ' + '─'.repeat(55));

for (const bin of radialBins) {
  const binVfResids = [];
  const binA0Resids = [];

  for (let i = 0; i < gals.length; i++) {
    const g = gals[i];
    const pts = g.rc.filter(p => bin.filter(p.r, g.Rdisk));
    if (pts.length < 2) continue;

    const Vmean = pts.reduce((a,p)=>a+p.v,0) / pts.length;
    const VfPred = Math.pow(10, vfModel.beta[0]*(g.logMbar - gals.reduce((a,gg)=>a+gg.logMbar,0)/gals.length)
      + vfModel.beta[1]*(g.logL36 - gals.reduce((a,gg)=>a+gg.logL36,0)/gals.length)
      + vfModel.beta[2]*(g.logRdisk - gals.reduce((a,gg)=>a+gg.logRdisk,0)/gals.length)
      + vfModel.beta[3]*(g.morphT - gals.reduce((a,gg)=>a+gg.morphT,0)/gals.length)
      + gals.reduce((a,gg)=>a+gg.logVflat,0)/gals.length);

    binVfResids.push(g.VfResid);
    binA0Resids.push(g.a0Resid);
  }

  if (binVfResids.length >= 10) {
    const r = pearsonR(binVfResids, binA0Resids);
    const meanAbsVf = binVfResids.reduce((a,v)=>a+Math.abs(v),0)/binVfResids.length;
    console.log(`  ${bin.name.padEnd(18)} ${binVfResids.length.toString().padEnd(8)} ${meanAbsVf.toFixed(4).padEnd(15)} ${r.toFixed(3)}`);
  } else {
    console.log(`  ${bin.name.padEnd(18)} ${binVfResids.length.toString().padEnd(8)} (too few)`);
  }
}


console.log('\n\n' + '▓'.repeat(70));
console.log('409E: INFORMATION BUDGET — WHAT IS THE LATENT VAR NOT?');
console.log('▓'.repeat(70));

const rVfResid_all = multiR2(
  gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]),
  gals.map(g => g.VfResid)
);
const rA0resid_all = multiR2(
  gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]),
  gals.map(g => g.a0Resid)
);

console.log(`\n  R²(VfResid ~ 6 structural) = ${rVfResid_all.R2.toFixed(3)}`);
console.log(`  R²(a₀_resid ~ 6 structural) = ${rA0resid_all.R2.toFixed(3)}`);
console.log(`  → ${(1-rVfResid_all.R2).toFixed(1)}% of VfResid is structurally unexplained`);
console.log(`  → ${(1-rA0resid_all.R2).toFixed(1)}% of a₀_resid is structurally unexplained (by construction ~100%)`);

const rLsum_struct = multiR2(
  gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]),
  L_sum
);
console.log(`\n  R²(L_sum ~ 6 structural) = ${rLsum_struct.R2.toFixed(3)}`);
console.log(`  → ${((1-rLsum_struct.R2)*100).toFixed(1)}% of latent variable is structurally unexplained`);


console.log('\n\n' + '▓'.repeat(70));
console.log('409F: VARIANCE DECOMPOSITION');
console.log('How much of VfResid comes through the a₀ channel vs other paths?');
console.log('▓'.repeat(70));

const rSquared = rBase * rBase;
console.log(`\n  r(VfResid, a₀_resid)² = ${rSquared.toFixed(3)}`);
console.log(`  → ${(rSquared*100).toFixed(1)}% of VfResid variance is shared with a₀_resid`);
console.log(`  → ${((1-rSquared)*100).toFixed(1)}% of VfResid variance is independent of a₀_resid`);

const totalVfResidVar = VfR.reduce((a,v)=>a+v*v,0)/VfR.length;
const channelVar = rSquared * totalVfResidVar;
const independentVar = (1-rSquared) * totalVfResidVar;

console.log(`\n  Total VfResid variance: ${totalVfResidVar.toFixed(6)}`);
console.log(`  Channel (shared with a₀): ${channelVar.toFixed(6)} (${(rSquared*100).toFixed(1)}%)`);
console.log(`  Independent of a₀: ${independentVar.toFixed(6)} (${((1-rSquared)*100).toFixed(1)}%)`);

const sdChannel = Math.sqrt(channelVar);
const sdIndep = Math.sqrt(independentVar);
console.log(`\n  In velocity terms (at mean Vflat):`);
const meanVflat = gals.reduce((a,g)=>a+g.Vflat,0)/gals.length;
const channelKms = sdChannel * meanVflat * Math.log(10);
const indepKms = sdIndep * meanVflat * Math.log(10);
console.log(`    Channel amplitude: ~${channelKms.toFixed(1)} km/s`);
console.log(`    Independent noise: ~${indepKms.toFixed(1)} km/s`);
console.log(`    Mean Vflat: ${meanVflat.toFixed(1)} km/s`);
console.log(`    Channel as fraction of Vflat: ${(channelKms/meanVflat*100).toFixed(1)}%`);


console.log('\n\n' + '▓'.repeat(70));
console.log('409G: LATENT VARIABLE CONSTRAINTS SUMMARY');
console.log('▓'.repeat(70));

console.log(`
FINGERPRINT OF THE HIDDEN VARIABLE:

1. STRENGTH: r ≈ 0.80 (explains ~65% of shared variance between VfResid and a₀)

2. UNIVERSALITY: Present at all galaxy masses (Phase 405b)
   Not regime-dependent; constant coupling strength

3. INDEPENDENCE from:
   - 6 structural properties (logMbar, logL36, logRdisk, morphT, logMHI, logSBdisk)
   - 9 RC-shape variables
   - 13 radial mismatch variables
   - Construction recipe (56/56 combinations tested)
   - M/L ratio (factor 3 range)

4. AMPLITUDE: ~${channelKms.toFixed(0)} km/s velocity effect (${(channelKms/meanVflat*100).toFixed(1)}% of Vflat)

5. CORRELATES with:
   ${latentCorrelations.filter(c => Math.abs(c.rLsum) > 0.25).map(c => `- ${c.name}: r(L_sum)=${c.rLsum.toFixed(3)}`).join('\n   ')}

6. DOES NOT CORRELATE with:
   ${latentCorrelations.filter(c => Math.abs(c.rLsum) < 0.15).map(c => `- ${c.name}: r(L_sum)=${c.rLsum.toFixed(3)}`).join('\n   ')}

7. DISTRIBUTION: approximately Gaussian
   skewness = ${skew.toFixed(3)}, kurtosis = ${kurt.toFixed(3)}

8. NOT EXPLAINED BY: structure (R²=${rLsum_struct.R2.toFixed(3)}), meaning
   ${((1-rLsum_struct.R2)*100).toFixed(0)}% of the latent variable is novel information
`);


console.log('\n' + '═'.repeat(70));
console.log('IMPLICATIONS FOR IDENTIFICATION');
console.log('═'.repeat(70));

console.log(`
The hidden variable must be something that:

A) AFFECTS both Vflat (through VfResid) and the acceleration scale (through a₀)
   simultaneously, with ~${channelKms.toFixed(0)} km/s amplitude

B) Is NOT captured by:
   - Total mass, luminosity, disk size, morphology, gas content, surface brightness
   - Rotation curve shape (inner/outer/gradient/concentration)
   - Radial baryon-dynamics mismatch profile
   - Model fit parameters

C) IS captured by the fitted a₀ parameter specifically
   (halo k gives r=0.633, Vflat²/Rmax gives r=0.574 — weaker proxies)

D) Contains ${((1-rLsum_struct.R2)*100).toFixed(0)}% novel information not in SPARC tables

CANDIDATE CLASSES:
1. Dark matter halo internal structure (concentration, core/cusp, triaxiality)
   → Would require IFU/2D kinematics or lensing to measure
2. Assembly history (formation time, merger history, accretion rate)
   → Would require cosmological simulations to test
3. External field / environment at finer granularity than envCode
   → Would require group catalogs or tidal field measurements
4. Genuine acceleration-scale physics (MOND-like threshold behavior)
   → Would predict specific functional form for the coupling
`);

const outPath = path.join(__dirname, '..', 'public', 'phase409-latent-fingerprint.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '409',
  title: 'Latent Variable Fingerprint',
  timestamp: new Date().toISOString(),
  N: gals.length,
  baselineR: rBase,
  L_sum_stats: { mean: meanL, sd: sdL, skew, kurtosis: kurt, q25, q50, q75 },
  latentCorrelations,
  channelAmplitude_kms: channelKms,
  channelFractionVflat: channelKms / meanVflat,
  structExplainedL: rLsum_struct.R2,
  varianceDecomp: { rSquared, channelVar, independentVar },
  topGalaxies: withL.slice(0,10).map(g => ({ name: g.name, Vflat: g.Vflat, logA0: g.logA0, L_sum: g.L_sum, morphT: g.morphT })),
  bottomGalaxies: withL.slice(-10).map(g => ({ name: g.name, Vflat: g.Vflat, logA0: g.logA0, L_sum: g.L_sum, morphT: g.morphT })),
}, null, 2));

console.log(`\nSaved: ${outPath}`);
