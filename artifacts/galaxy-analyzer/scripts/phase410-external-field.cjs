const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

function normalize(n) {
  return n.replace(/\s+/g,'').replace(/^NGC0*/,'NGC').replace(/^DDO0*/,'DDO').replace(/^IC0*/,'IC').replace(/^UGC0*/,'UGC').replace(/^PGC0*/,'PGC').toUpperCase();
}
function pearsonR(x, y) {
  const n = x.length; if (n < 4) return NaN;
  const mx = x.reduce((a,b)=>a+b,0)/n, my = y.reduce((a,b)=>a+b,0)/n;
  let num=0, dx2=0, dy2=0;
  for (let i=0;i<n;i++){const dx=x[i]-mx,dy=y[i]-my;num+=dx*dy;dx2+=dx*dx;dy2+=dy*dy;}
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
  return {R2:sst>0?1-sse/sst:0,beta,residuals,coeffs:beta,intercept:my};
}
function partialR(x, y, controls) {
  if (controls.length === 0 || controls[0].length === 0) return pearsonR(x, y);
  const n = x.length;
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
  gals.push({
    name: g.name, logA0: g.logA0,
    Vflat: sp.Vflat, Rdisk: sp.Rdisk,
    logVflat: Math.log10(sp.Vflat),
    logL36: Math.log10(Math.max(sp.L36, 0.001)),
    logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)),
    logMbar: Math.log10(Math.max(sp.L36*0.5 + sp.MHI*1.33, 0.001) * 1e9),
    logMHI: g.logMHI,
    morphT: sp.T,
    logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    envCode: g.envCode,
    L36: sp.L36, MHI: sp.MHI, SBdisk: sp.SBdisk,
    inc: sp.inc, Q: sp.Q, D: sp.D,
    k_halo: sr.models.dark_halo_linear.k,
    M_halo: sr.models.dark_halo_linear.M,
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
console.log('PHASE 410: EXTERNAL-FIELD DISCRIMINATION');
console.log('Is the hidden variable an external gravitational field?');
console.log('═'.repeat(70));
console.log(`N = ${gals.length}, baseline r(VfResid, a₀_resid) = ${rBase.toFixed(3)}`);


console.log('\n' + '▓'.repeat(70));
console.log('410A: ENVIRONMENT STRATIFICATION');
console.log('How does the channel behave within each environment class?');
console.log('▓'.repeat(70));

const envGroups = { 0: 'Field', 1: 'Group', 2: 'Cluster' };
for (const [code, label] of Object.entries(envGroups)) {
  const sub = gals.filter(g => g.envCode === parseInt(code));
  if (sub.length < 5) { console.log(`\n  ${label} (envCode=${code}): N=${sub.length} — too few`); continue; }
  const r = pearsonR(sub.map(g=>g.VfResid), sub.map(g=>g.a0Resid));
  const meanVfR = sub.reduce((a,g)=>a+g.VfResid,0)/sub.length;
  const meanA0R = sub.reduce((a,g)=>a+g.a0Resid,0)/sub.length;
  const meanL = sub.reduce((a,g)=>a+g.L_sum,0)/sub.length;
  const meanA0 = sub.reduce((a,g)=>a+g.logA0,0)/sub.length;
  console.log(`\n  ${label} (envCode=${code}): N=${sub.length}`);
  console.log(`    r(VfResid, a₀_resid) = ${r.toFixed(3)}`);
  console.log(`    mean VfResid = ${meanVfR>=0?'+':''}${meanVfR.toFixed(4)}`);
  console.log(`    mean a₀_resid = ${meanA0R>=0?'+':''}${meanA0R.toFixed(4)}`);
  console.log(`    mean L_sum = ${meanL>=0?'+':''}${meanL.toFixed(3)}`);
  console.log(`    mean logA0 = ${meanA0.toFixed(3)}`);
}

const fieldGals = gals.filter(g => g.envCode === 0);
const nonFieldGals = gals.filter(g => g.envCode > 0);
const rField = pearsonR(fieldGals.map(g=>g.VfResid), fieldGals.map(g=>g.a0Resid));
const rNonField = pearsonR(nonFieldGals.map(g=>g.VfResid), nonFieldGals.map(g=>g.a0Resid));
console.log(`\n  CRITICAL: Within-environment channel:`);
console.log(`    Field only (N=${fieldGals.length}): r = ${rField.toFixed(3)}`);
console.log(`    Non-field only (N=${nonFieldGals.length}): r = ${rNonField.toFixed(3)}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('410B: EFE FUNCTIONAL FORM TEST');
console.log('Does the channel follow MOND-EFE predictions?');
console.log('▓'.repeat(70));

console.log('\n  In MOND, the external field effect predicts:');
console.log('  - Galaxies in strong external fields → higher effective a₀');
console.log('  - The effect: a₀_eff ≈ a₀_int × (1 + g_ext/a₀_int)');
console.log('  - So cluster gals should show HIGHER a₀, NOT lower');
console.log('  - But we see envCode NEGATIVELY correlated with L_sum');
console.log('  - Meaning cluster gals have LOWER L_sum (weaker channel)');

const meanA0_field = fieldGals.reduce((a,g)=>a+g.logA0,0)/fieldGals.length;
const meanA0_nonfield = nonFieldGals.reduce((a,g)=>a+g.logA0,0)/nonFieldGals.length;
const meanA0resid_field = fieldGals.reduce((a,g)=>a+g.a0Resid,0)/fieldGals.length;
const meanA0resid_nonfield = nonFieldGals.reduce((a,g)=>a+g.a0Resid,0)/nonFieldGals.length;

console.log(`\n  logA0 by environment:`);
console.log(`    Field:     mean logA0 = ${meanA0_field.toFixed(3)}, mean a₀_resid = ${meanA0resid_field>=0?'+':''}${meanA0resid_field.toFixed(4)}`);
console.log(`    Non-field: mean logA0 = ${meanA0_nonfield.toFixed(3)}, mean a₀_resid = ${meanA0resid_nonfield>=0?'+':''}${meanA0resid_nonfield.toFixed(4)}`);

const a0_diff = meanA0_field - meanA0_nonfield;
const a0R_diff = meanA0resid_field - meanA0resid_nonfield;
console.log(`    Δ(logA0) = ${a0_diff>=0?'+':''}${a0_diff.toFixed(4)} (field − non-field)`);
console.log(`    Δ(a₀_resid) = ${a0R_diff>=0?'+':''}${a0R_diff.toFixed(4)}`);

const efeSign = a0R_diff > 0 ? 'FIELD higher' : 'NON-FIELD higher';
console.log(`\n  EFE SIGN TEST:`);
console.log(`    Observation: ${efeSign}`);
console.log(`    MOND-EFE expects: cluster/group galaxies → BOOSTED a₀_eff`);
console.log(`    But standard EFE with interpolation can LOWER a₀ in phantom regime`);
console.log(`    → Sign alone is inconclusive. Need functional form.`);


console.log('\n\n' + '▓'.repeat(70));
console.log('410C: DOES ENVIRONMENT ABSORB THE CHANNEL?');
console.log('Partial correlation after controlling for envCode');
console.log('▓'.repeat(70));

const envDummies = gals.map(g => [g.envCode === 1 ? 1 : 0, g.envCode === 2 ? 1 : 0]);
const r_partial_env = partialR(
  gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid),
  envDummies
);
console.log(`\n  r(VfResid, a₀_resid) baseline = ${rBase.toFixed(3)}`);
console.log(`  r(VfResid, a₀_resid | envCode) = ${r_partial_env.toFixed(3)}`);
console.log(`  Δr = ${(r_partial_env - rBase).toFixed(3)}`);
console.log(`  Absorption: ${Math.abs(r_partial_env - rBase) > 0.05 ? 'PARTIAL' : 'NONE'}`);

const envPlusHaloControls = gals.map(g => [
  g.envCode === 1 ? 1 : 0,
  g.envCode === 2 ? 1 : 0,
  g.k_halo > 0 ? Math.log10(g.k_halo) : 0,
]);
const r_partial_envHalo = partialR(
  gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid),
  envPlusHaloControls
);
console.log(`\n  r(VfResid, a₀_resid | envCode + logK_halo) = ${r_partial_envHalo.toFixed(3)}`);
console.log(`  Δr = ${(r_partial_envHalo - rBase).toFixed(3)}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('410D: SYNTHETIC EXTERNAL-FIELD PROXIES');
console.log('Building quantitative EFE-like variables');
console.log('▓'.repeat(70));

const a0_MOND = 1.2e-10;

const efeProxies = [];

const g_ext_proxy1 = gals.map(g => {
  const envStrength = [0.0, 0.5, 1.0][g.envCode];
  return envStrength;
});
efeProxies.push({ name: 'env_linear', desc: 'envCode as 0/0.5/1', vals: g_ext_proxy1 });

const g_ext_proxy2 = gals.map(g => {
  const envStrength = [0.01, 0.1, 1.0][g.envCode];
  return Math.log10(envStrength);
});
efeProxies.push({ name: 'env_log', desc: 'log(envStrength) [0.01,0.1,1]', vals: g_ext_proxy2 });

const g_ext_proxy3 = gals.map(g => {
  const a_int = Math.pow(10, g.logA0);
  const g_ext = [0, 0.1*a_int, 0.3*a_int][g.envCode];
  return g_ext > 0 ? Math.log10(g_ext) : -15;
});
efeProxies.push({ name: 'g_ext_frac_a0', desc: 'g_ext as fraction of a₀_int', vals: g_ext_proxy3 });

const g_ext_proxy4 = gals.map(g => {
  const a_int = Math.pow(10, g.logA0);
  const g_ext = [0, 0.1*a_int, 0.3*a_int][g.envCode];
  const ratio = g_ext / (a_int + g_ext + 1e-20);
  return ratio;
});
efeProxies.push({ name: 'g_ext_ratio', desc: 'g_ext/(a₀_int + g_ext)', vals: g_ext_proxy4 });

const efeCorrection = gals.map(g => {
  const a_int = Math.pow(10, g.logA0);
  const g_ext_vals = [0, 0.1, 0.3];
  const g_ext = g_ext_vals[g.envCode] * a_int;
  const a_eff = Math.sqrt(a_int*a_int + g_ext*g_ext);
  return Math.log10(a_eff) - g.logA0;
});
efeProxies.push({ name: 'efe_correction', desc: 'log(a_eff/a_int) quadrature model', vals: efeCorrection });

const efeInterp = gals.map(g => {
  const a_int = Math.pow(10, g.logA0);
  const g_ext_vals = [0, 0.1, 0.3];
  const g_ext = g_ext_vals[g.envCode] * a_int;
  const nu = 1 / (1 - Math.exp(-Math.sqrt(a_int / (a0_MOND + g_ext))));
  return Math.log10(Math.max(nu, 0.01));
});
efeProxies.push({ name: 'efe_interp', desc: 'MOND ν(a_int/(a₀+g_ext))', vals: efeInterp });

console.log('\n  Proxy               r(VfResid) r(a₀_resid) r(L_sum)   absorb?');
console.log('  ' + '─'.repeat(65));

for (const proxy of efeProxies) {
  const valid = proxy.vals.map((v,i) => isFinite(v) ? i : -1).filter(i=>i>=0);
  if (valid.length < 10) continue;
  const rVf = pearsonR(valid.map(i=>gals[i].VfResid), valid.map(i=>proxy.vals[i]));
  const rA0 = pearsonR(valid.map(i=>gals[i].a0Resid), valid.map(i=>proxy.vals[i]));
  const rL = pearsonR(valid.map(i=>gals[i].L_sum), valid.map(i=>proxy.vals[i]));

  const rPartial = partialR(
    valid.map(i=>gals[i].VfResid),
    valid.map(i=>gals[i].a0Resid),
    valid.map(i=>[proxy.vals[i]])
  );
  const absorbed = rBase - rPartial > 0.05;

  console.log(`  ${proxy.name.padEnd(20)} ${rVf.toFixed(3).padEnd(11)} ${rA0.toFixed(3).padEnd(12)} ${rL.toFixed(3).padEnd(11)} ${absorbed ? 'PARTIAL':'NONE'} (r_part=${rPartial.toFixed(3)})`);
}


console.log('\n\n' + '▓'.repeat(70));
console.log('410E: ACCELERATION-SCALE DEPENDENCE OF ENVIRONMENT EFFECT');
console.log('Does the environment effect scale with internal acceleration?');
console.log('▓'.repeat(70));

const lowA0 = gals.filter(g => g.logA0 < 3.6);
const highA0 = gals.filter(g => g.logA0 >= 3.6);

console.log(`\n  Split at logA0 = 3.6:`);
console.log(`    Low a₀ (N=${lowA0.length}): mean envCode = ${(lowA0.reduce((a,g)=>a+g.envCode,0)/lowA0.length).toFixed(2)}`);
console.log(`    High a₀ (N=${highA0.length}): mean envCode = ${(highA0.reduce((a,g)=>a+g.envCode,0)/highA0.length).toFixed(2)}`);

if (lowA0.length >= 8 && highA0.length >= 8) {
  const rLow = pearsonR(lowA0.map(g=>g.envCode), lowA0.map(g=>g.L_sum));
  const rHigh = pearsonR(highA0.map(g=>g.envCode), highA0.map(g=>g.L_sum));
  console.log(`    r(envCode, L_sum) in low a₀: ${rLow.toFixed(3)}`);
  console.log(`    r(envCode, L_sum) in high a₀: ${rHigh.toFixed(3)}`);
  console.log(`\n  EFE prediction: effect should be STRONGER at low acceleration`);
  console.log(`    (where g_ext/a_int ratio is larger)`);
  console.log(`    Observed: ${Math.abs(rLow) > Math.abs(rHigh) ? 'LOW a₀ stronger ✓ (EFE-consistent)' : 'HIGH a₀ stronger ✗ (EFE-inconsistent)'}`);
}


console.log('\n\n' + '▓'.repeat(70));
console.log('410F: ENVIRONMENT × CHANNEL INTERACTION');
console.log('Does environment modulate the channel slope, or just shift it?');
console.log('▓'.repeat(70));

for (const [code, label] of Object.entries(envGroups)) {
  const sub = gals.filter(g => g.envCode === parseInt(code));
  if (sub.length < 8) continue;

  const slope = (() => {
    const x = sub.map(g => g.a0Resid);
    const y = sub.map(g => g.VfResid);
    const mx = x.reduce((a,b)=>a+b,0)/x.length;
    const my = y.reduce((a,b)=>a+b,0)/y.length;
    let num = 0, den = 0;
    for (let i = 0; i < x.length; i++) { num += (x[i]-mx)*(y[i]-my); den += (x[i]-mx)**2; }
    return den > 0 ? num/den : 0;
  })();

  const r = pearsonR(sub.map(g=>g.VfResid), sub.map(g=>g.a0Resid));
  console.log(`  ${label.padEnd(10)} N=${sub.length.toString().padEnd(4)} slope(VfResid~a₀resid) = ${slope.toFixed(4)}, r = ${r.toFixed(3)}`);
}

console.log('\n  EFE prediction: slope should DECREASE in stronger environments');
console.log('  (external field suppresses the internal coupling)');


console.log('\n\n' + '▓'.repeat(70));
console.log('410G: COMBINED EXTERNAL-FIELD MODEL');
console.log('Best possible absorption using all environmental info');
console.log('▓'.repeat(70));

const allEnvControls = gals.map(g => [
  g.envCode === 1 ? 1 : 0,
  g.envCode === 2 ? 1 : 0,
  g.envCode * g.logA0,
  g.envCode * g.logVflat,
  g.k_halo > 0 ? Math.log10(g.k_halo) : 0,
]);

const r_full_env = partialR(
  gals.map(g=>g.VfResid), gals.map(g=>g.a0Resid),
  allEnvControls
);

console.log(`\n  Controls: envCode dummies + env×logA0 + env×logVflat + logK_halo`);
console.log(`  r(VfResid, a₀_resid | all env) = ${r_full_env.toFixed(3)}`);
console.log(`  Δr = ${(r_full_env - rBase).toFixed(3)}`);
console.log(`  Absorption: ${Math.abs(r_full_env - rBase) > 0.10 ? 'SIGNIFICANT' : Math.abs(r_full_env - rBase) > 0.05 ? 'PARTIAL' : 'NEGLIGIBLE'}`);

const rVfResid_envModel = multiR2(allEnvControls, gals.map(g=>g.VfResid));
const rA0Resid_envModel = multiR2(allEnvControls, gals.map(g=>g.a0Resid));
console.log(`\n  R²(VfResid ~ all env) = ${rVfResid_envModel.R2.toFixed(3)}`);
console.log(`  R²(a₀_resid ~ all env) = ${rA0Resid_envModel.R2.toFixed(3)}`);

const r_afterEnv = pearsonR(rVfResid_envModel.residuals, rA0Resid_envModel.residuals);
console.log(`  r(VfResid_clean, a₀resid_clean) = ${r_afterEnv.toFixed(3)}`);
console.log(`  → Channel after removing ALL environmental info: ${r_afterEnv.toFixed(3)}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('410H: EFE vs NON-EFE DISCRIMINATION');
console.log('▓'.repeat(70));

const efeScorecard = [];

const test1_signCorrect = meanA0resid_field > meanA0resid_nonfield;
efeScorecard.push({ test: 'Sign: field gals have higher a₀_resid', result: test1_signCorrect });

const test2_absorbed = Math.abs(r_partial_env - rBase) > 0.10;
efeScorecard.push({ test: 'envCode absorbs >0.10 of channel', result: test2_absorbed });

const test3_withinField = rField > 0.5;
efeScorecard.push({ test: 'Channel persists within field-only (r>0.5)', result: test3_withinField });

const test4_withinNonField = rNonField > 0.3;
efeScorecard.push({ test: 'Channel persists within non-field (r>0.3)', result: test4_withinNonField });

let test5_lowStronger = false;
if (lowA0.length >= 8 && highA0.length >= 8) {
  const rLow = pearsonR(lowA0.map(g=>g.envCode), lowA0.map(g=>g.L_sum));
  const rHigh = pearsonR(highA0.map(g=>g.envCode), highA0.map(g=>g.L_sum));
  test5_lowStronger = Math.abs(rLow) > Math.abs(rHigh);
}
efeScorecard.push({ test: 'EFE stronger at low acceleration', result: test5_lowStronger });

const test6_afterAllEnv = r_afterEnv > 0.65;
efeScorecard.push({ test: 'Channel survives after ALL env removal (r>0.65)', result: test6_afterAllEnv });

console.log('\n  EFE SCORECARD:');
console.log('  ' + '─'.repeat(60));
for (const t of efeScorecard) {
  const icon = t.result ? 'PASS' : 'FAIL';
  console.log(`  [${icon}] ${t.test}`);
}
const passCount = efeScorecard.filter(t => t.result).length;
console.log('  ' + '─'.repeat(60));
console.log(`  Score: ${passCount}/${efeScorecard.length}`);


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 410 GRAND VERDICT');
console.log('═'.repeat(70));

const channelAfterEnv = r_afterEnv;
const envExplainsChannel = Math.abs(rBase - channelAfterEnv) > 0.15;

console.log(`
External-Field Discrimination Summary:

1. Environment CORRELATES with latent variable: r(envCode, L_sum) = -0.42
   → Field galaxies carry the channel more strongly

2. Within-environment channel:
   Field only: r = ${rField.toFixed(3)}
   Non-field: r = ${rNonField.toFixed(3)}
   → Channel persists WITHIN environment classes

3. After removing ALL environmental information:
   r(VfResid, a₀_resid) = ${channelAfterEnv.toFixed(3)} (from ${rBase.toFixed(3)})
   → ${envExplainsChannel ? 'SIGNIFICANT reduction — environment matters' : 'Channel largely SURVIVES — environment is not the primary driver'}

4. EFE scorecard: ${passCount}/${efeScorecard.length} tests passed

VERDICT: ${envExplainsChannel
  ? 'Environment explains a significant portion of the channel. EFE remains a viable explanation.'
  : 'Environment is a CORRELATE but not the PRIMARY DRIVER. The hidden variable has an environmental footprint but is not reducible to external-field effects. The channel persists strongly after complete environmental decontamination.'}
`);

if (!envExplainsChannel) {
  console.log('IMPLICATION: The hidden variable is NOT primarily an external field.');
  console.log('The environment correlation (r=-0.42) may arise because:');
  console.log('  a) The true hidden variable correlates with environment (e.g., halo formation history)');
  console.log('  b) Environment selection effects create apparent L_sum-env correlation');
  console.log('  c) The hidden variable has both internal and environmental components');
  console.log('\n→ NEXT: Test MOND-like functional form (Phase 411)');
  console.log('→ Then: Halo concentration/structure (Phase 412)');
}

const outPath = path.join(__dirname, '..', 'public', 'phase410-external-field.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '410',
  title: 'External-Field Discrimination',
  timestamp: new Date().toISOString(),
  N: gals.length,
  baselineR: rBase,
  withinEnvR: { field: rField, nonField: rNonField },
  r_partial_env,
  r_afterAllEnv: channelAfterEnv,
  efeScorecard,
  envExplainsChannel,
  meanA0_field, meanA0_nonfield,
  meanA0resid_field, meanA0resid_nonfield,
}, null, 2));

console.log(`\nSaved: ${outPath}`);
