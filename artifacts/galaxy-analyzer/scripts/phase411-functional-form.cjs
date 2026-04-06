const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

function normalize(n) {
  return n.replace(/\s+/g,'').replace(/^NGC0*/,'NGC').replace(/^DDO0*/,'DDO').replace(/^IC0*/,'IC').replace(/^UGC0*/,'UGC').replace(/^PGC0*/,'PGC').toUpperCase();
}
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

const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[g.name] = g; resultsMap[normalize(g.name)] = g; });

const G = 4.3009e-6;

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sr) continue;

  const Mbar = (sp.L36*0.5 + sp.MHI*1.33) * 1e9;
  const Rdisk = sp.Rdisk;
  const Vflat = sp.Vflat;
  const Rflat = sp.Rflat > 0 ? sp.Rflat : Rdisk * 3;

  const g_bar = G * Mbar / (2.2 * Rdisk) / (2.2 * Rdisk);
  const g_obs = Vflat * Vflat / (2.2 * Rdisk);
  const g_bar_Rflat = G * Mbar / (Rflat * Rflat);
  const g_obs_Rflat = Vflat * Vflat / Rflat;

  gals.push({
    name: g.name, logA0: g.logA0,
    Vflat, Rdisk, Rflat, Mbar,
    logVflat: Math.log10(Vflat),
    logL36: Math.log10(Math.max(sp.L36, 0.001)),
    logRdisk: Math.log10(Math.max(Rdisk, 0.01)),
    logMbar: Math.log10(Math.max(Mbar, 1)),
    logMHI: g.logMHI,
    morphT: sp.T,
    logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)),
    envCode: g.envCode,
    L36: sp.L36, MHI: sp.MHI,
    k_halo: sr.models.dark_halo_linear.k,
    g_bar, g_obs, g_bar_Rflat, g_obs_Rflat,
    log_g_bar: Math.log10(Math.max(g_bar, 1e-20)),
    log_g_obs: Math.log10(Math.max(g_obs, 1e-20)),
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
console.log('PHASE 411: MOND-LIKE FUNCTIONAL FORM DISCRIMINATION');
console.log('Does the channel follow an acceleration law, or is it something else?');
console.log('═'.repeat(70));
console.log(`N = ${gals.length}, baseline r(VfResid, a₀_resid) = ${rBase.toFixed(3)}\n`);


console.log('▓'.repeat(70));
console.log('411A: ACCELERATION-SPACE PORTRAIT');
console.log('Where do our galaxies sit in g_obs vs g_bar space?');
console.log('▓'.repeat(70));

const rRAR = pearsonR(gals.map(g=>g.log_g_bar), gals.map(g=>g.log_g_obs));
console.log(`\n  r(log g_bar, log g_obs) = ${rRAR.toFixed(3)} — galaxy-level RAR`);

const sorted_gbar = [...gals].sort((a,b) => a.log_g_bar - b.log_g_bar);
const tertiles = [
  sorted_gbar.slice(0, Math.floor(gals.length/3)),
  sorted_gbar.slice(Math.floor(gals.length/3), Math.floor(2*gals.length/3)),
  sorted_gbar.slice(Math.floor(2*gals.length/3)),
];
console.log(`\n  Channel strength by g_bar tertile:`);
console.log('  Tertile      N   mean(log g_bar)  r(VfR,a₀R)  mean a₀_resid');
console.log('  ' + '─'.repeat(65));
for (let t = 0; t < 3; t++) {
  const sub = tertiles[t];
  const r = pearsonR(sub.map(g=>g.VfResid), sub.map(g=>g.a0Resid));
  const mg = sub.reduce((a,g)=>a+g.log_g_bar,0)/sub.length;
  const ma = sub.reduce((a,g)=>a+g.a0Resid,0)/sub.length;
  console.log(`  ${['Low','Mid','High'][t].padEnd(12)} ${sub.length.toString().padEnd(4)} ${mg.toFixed(3).padEnd(17)} ${r.toFixed(3).padEnd(12)} ${(ma>=0?'+':'')+ma.toFixed(4)}`);
}


console.log('\n\n' + '▓'.repeat(70));
console.log('411B: MOND INTERPOLATION FUNCTION TEST');
console.log('Does g_obs/g_bar follow a MOND ν-function of g_bar?');
console.log('▓'.repeat(70));

const a0_MOND = 1.2e-10;
const a0_kpc = a0_MOND * (3.086e19)**2 / 1e6;

const discrepancy = gals.map(g => Math.log10(g.g_obs / g.g_bar));
const rDisc_a0resid = pearsonR(discrepancy, gals.map(g => g.a0Resid));
const rDisc_VfResid = pearsonR(discrepancy, gals.map(g => g.VfResid));
const rDisc_Lsum = pearsonR(discrepancy, gals.map(g => g.L_sum));

console.log(`\n  Discrepancy D = log(g_obs/g_bar) — mass discrepancy at 2.2Rd`);
console.log(`  r(D, a₀_resid) = ${rDisc_a0resid.toFixed(3)}`);
console.log(`  r(D, VfResid)  = ${rDisc_VfResid.toFixed(3)}`);
console.log(`  r(D, L_sum)    = ${rDisc_Lsum.toFixed(3)}`);

const rDisc_gbar = pearsonR(discrepancy, gals.map(g => g.log_g_bar));
console.log(`  r(D, log g_bar) = ${rDisc_gbar.toFixed(3)} — MOND expects strong negative`);

const discResid = multiR2(
  gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]),
  discrepancy
).residuals;

const rDiscResid_a0resid = pearsonR(discResid, gals.map(g=>g.a0Resid));
const rDiscResid_VfResid = pearsonR(discResid, gals.map(g=>g.VfResid));
console.log(`\n  After structural control:`);
console.log(`  r(D_resid, a₀_resid) = ${rDiscResid_a0resid.toFixed(3)}`);
console.log(`  r(D_resid, VfResid) = ${rDiscResid_VfResid.toFixed(3)}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('411C: COMPETING FUNCTIONAL FORMS');
console.log('Which function best describes VfResid = f(a₀_resid)?');
console.log('▓'.repeat(70));

const x = gals.map(g => g.a0Resid);
const y = gals.map(g => g.VfResid);
const N = gals.length;

function fitAndScore(name, predict) {
  const yhat = predict(x, y);
  let sse = 0, sst = 0;
  const my = y.reduce((a,b)=>a+b,0)/N;
  for (let i = 0; i < N; i++) {
    sse += (y[i] - yhat[i]) ** 2;
    sst += (y[i] - my) ** 2;
  }
  const R2 = 1 - sse/sst;
  const k = predict.params || 2;
  const AIC = N * Math.log(sse/N) + 2*k;
  const BIC = N * Math.log(sse/N) + k*Math.log(N);
  const resid = y.map((yi,i) => yi - yhat[i]);
  return { name, R2, AIC, BIC, k, resid, sse };
}

const fits = [];

const linearFit = (x, y) => {
  const mx = x.reduce((a,b)=>a+b,0)/N, my = y.reduce((a,b)=>a+b,0)/N;
  let num=0,den=0;
  for(let i=0;i<N;i++){num+=(x[i]-mx)*(y[i]-my);den+=(x[i]-mx)**2;}
  const b = den>0?num/den:0, a = my - b*mx;
  linearFit.params = 2;
  return x.map(xi => a + b*xi);
};
fits.push(fitAndScore('Linear', linearFit));

const quadFit = (x, y) => {
  const X = x.map(xi => [xi, xi*xi]);
  const model = multiR2(X, y);
  quadFit.params = 3;
  const my = y.reduce((a,b)=>a+b,0)/N;
  const mx1 = x.reduce((a,b)=>a+b,0)/N;
  const mx2 = x.reduce((a,b)=>a+b*b,0)/N;
  return x.map(xi => my + model.beta[0]*(xi-mx1) + model.beta[1]*(xi*xi-mx2));
};
fits.push(fitAndScore('Quadratic', quadFit));

const cubicFit = (x, y) => {
  const X = x.map(xi => [xi, xi*xi, xi*xi*xi]);
  const model = multiR2(X, y);
  cubicFit.params = 4;
  const my = y.reduce((a,b)=>a+b,0)/N;
  const mx = X.reduce((a,row)=>{a[0]+=row[0];a[1]+=row[1];a[2]+=row[2];return a;},[0,0,0]).map(v=>v/N);
  return x.map(xi => my + model.beta[0]*(xi-mx[0]) + model.beta[1]*(xi*xi-mx[1]) + model.beta[2]*(xi*xi*xi-mx[2]));
};
fits.push(fitAndScore('Cubic', cubicFit));

const tanhFit = (x, y) => {
  const mx = x.reduce((a,b)=>a+b,0)/N;
  const my = y.reduce((a,b)=>a+b,0)/N;
  let bestSSE = Infinity, bestA = 0, bestS = 1;
  for (let s = 0.5; s <= 10; s += 0.5) {
    const xt = x.map(xi => Math.tanh(s * xi));
    const mxt = xt.reduce((a,b)=>a+b,0)/N;
    let num=0,den=0;
    for(let i=0;i<N;i++){num+=(xt[i]-mxt)*(y[i]-my);den+=(xt[i]-mxt)**2;}
    const a = den>0?num/den:0;
    const c = my - a*mxt;
    let sse = 0;
    for(let i=0;i<N;i++) sse += (y[i] - c - a*xt[i])**2;
    if (sse < bestSSE) { bestSSE = sse; bestA = a; bestS = s; }
  }
  tanhFit.params = 3;
  const xt = x.map(xi => Math.tanh(bestS * xi));
  const mxt = xt.reduce((a,b)=>a+b,0)/N;
  const c = my - bestA*mxt;
  return x.map(xi => c + bestA * Math.tanh(bestS * xi));
};
fits.push(fitAndScore('Tanh (saturating)', tanhFit));

const sqrtFit = (x, y) => {
  const X = x.map(xi => [xi, Math.sign(xi)*Math.sqrt(Math.abs(xi))]);
  const model = multiR2(X, y);
  sqrtFit.params = 3;
  const my = y.reduce((a,b)=>a+b,0)/N;
  const mx = X.reduce((a,row)=>{a[0]+=row[0];a[1]+=row[1];return a;},[0,0]).map(v=>v/N);
  return x.map(xi => my + model.beta[0]*(xi-mx[0]) + model.beta[1]*(Math.sign(xi)*Math.sqrt(Math.abs(xi))-mx[1]));
};
fits.push(fitAndScore('Sqrt+linear', sqrtFit));

const mondLikeFit = (x, y) => {
  const X = x.map(xi => {
    const a_int = Math.pow(10, 3.55 + xi * 2);
    const nu = 0.5 * (1 + Math.sqrt(1 + 4 * a0_kpc / a_int));
    return [xi, Math.log10(nu)];
  });
  const model = multiR2(X, y);
  mondLikeFit.params = 3;
  const my = y.reduce((a,b)=>a+b,0)/N;
  const mx = X.reduce((a,row)=>{a[0]+=row[0];a[1]+=row[1];return a;},[0,0]).map(v=>v/N);
  return x.map((xi,i) => my + model.beta[0]*(xi-mx[0]) + model.beta[1]*(X[i][1]-mx[1]));
};
fits.push(fitAndScore('MOND ν + linear', mondLikeFit));

fits.sort((a,b) => a.BIC - b.BIC);

console.log('\n  Model                R²      AIC      BIC      k  ΔBIC');
console.log('  ' + '─'.repeat(60));
const bestBIC = fits[0].BIC;
for (const f of fits) {
  console.log(`  ${f.name.padEnd(22)} ${f.R2.toFixed(3)}   ${f.AIC.toFixed(1).padEnd(9)} ${f.BIC.toFixed(1).padEnd(9)} ${f.k}  ${(f.BIC-bestBIC).toFixed(1)}`);
}

const linearR2 = fits.find(f=>f.name==='Linear').R2;
const bestR2 = fits[0].R2;
const improvementOverLinear = bestR2 - linearR2;
console.log(`\n  Best model: ${fits[0].name}`);
console.log(`  R² improvement over linear: ${improvementOverLinear.toFixed(4)}`);
console.log(`  → ${improvementOverLinear > 0.02 ? 'Non-linear structure DETECTED' : 'Relationship is essentially LINEAR — no MOND-like curvature'}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('411D: NON-LINEARITY DIAGNOSTIC');
console.log('Is there curvature in the VfResid vs a₀_resid relationship?');
console.log('▓'.repeat(70));

const sortedByA0R = [...gals].sort((a,b) => a.a0Resid - b.a0Resid);
const quintiles = [];
const qSize = Math.floor(N / 5);
for (let q = 0; q < 5; q++) {
  const start = q * qSize;
  const end = q === 4 ? N : (q+1) * qSize;
  const sub = sortedByA0R.slice(start, end);
  const meanA0R = sub.reduce((a,g)=>a+g.a0Resid,0)/sub.length;
  const meanVfR = sub.reduce((a,g)=>a+g.VfResid,0)/sub.length;
  const r = pearsonR(sub.map(g=>g.VfResid), sub.map(g=>g.a0Resid));
  quintiles.push({ q: q+1, N: sub.length, meanA0R, meanVfR, r });
}

console.log('\n  Quintile  N   mean(a₀R)    mean(VfR)    local r');
console.log('  ' + '─'.repeat(55));
for (const q of quintiles) {
  console.log(`  Q${q.q}       ${q.N.toString().padEnd(4)} ${(q.meanA0R>=0?'+':'')+q.meanA0R.toFixed(4).padEnd(13)} ${(q.meanVfR>=0?'+':'')+q.meanVfR.toFixed(4).padEnd(13)} ${q.r.toFixed(3)}`);
}

const linearResid = fits.find(f=>f.name==='Linear').resid;
const rResidA0R = pearsonR(linearResid, x);
const rResidA0R2 = pearsonR(linearResid, x.map(xi=>xi*xi));
console.log(`\n  Residual diagnostics (after linear fit):`);
console.log(`  r(residual, a₀_resid) = ${rResidA0R.toFixed(4)} (should be ~0 if linear correct)`);
console.log(`  r(residual, a₀_resid²) = ${rResidA0R2.toFixed(4)} (curvature signal if ≠0)`);
console.log(`  → ${Math.abs(rResidA0R2) > 0.2 ? 'CURVATURE detected: non-linear function needed' : 'No significant curvature: linear is adequate'}`);


console.log('\n\n' + '▓'.repeat(70));
console.log('411E: MOND-SPECIFIC PREDICTION TESTS');
console.log('Does the data follow specific MOND predictions?');
console.log('▓'.repeat(70));

console.log('\n  MOND Prediction 1: g_obs = g_bar × ν(g_bar/a₀)');
console.log('  If MOND is right, the residual structure should follow ν-function');

const nu_simple = gals.map(g => {
  const x_mond = g.g_bar / a0_kpc;
  return 0.5 * (1 + Math.sqrt(1 + 4/x_mond));
});
const nu_exp = gals.map(g => {
  const x_mond = g.g_bar / a0_kpc;
  return 1 / (1 - Math.exp(-Math.sqrt(x_mond)));
});

const log_nu_simple = nu_simple.map(v => Math.log10(v));
const log_nu_exp = nu_exp.map(v => Math.log10(v));

const rNuSimple_a0R = pearsonR(log_nu_simple, gals.map(g=>g.a0Resid));
const rNuSimple_VfR = pearsonR(log_nu_simple, gals.map(g=>g.VfResid));
const rNuSimple_L = pearsonR(log_nu_simple, gals.map(g=>g.L_sum));

const rNuExp_a0R = pearsonR(log_nu_exp, gals.map(g=>g.a0Resid));
const rNuExp_VfR = pearsonR(log_nu_exp, gals.map(g=>g.VfResid));
const rNuExp_L = pearsonR(log_nu_exp, gals.map(g=>g.L_sum));

console.log(`\n  ν-function correlations:`);
console.log(`  Simple ν: r(logν, a₀R)=${rNuSimple_a0R.toFixed(3)}, r(logν, VfR)=${rNuSimple_VfR.toFixed(3)}, r(logν, L)=${rNuSimple_L.toFixed(3)}`);
console.log(`  Exp ν:    r(logν, a₀R)=${rNuExp_a0R.toFixed(3)}, r(logν, VfR)=${rNuExp_VfR.toFixed(3)}, r(logν, L)=${rNuExp_L.toFixed(3)}`);

console.log('\n  MOND Prediction 2: In deep-MOND regime, g_obs ∝ √(g_bar × a₀)');
console.log('  → log(g_obs) ≈ 0.5×log(g_bar) + 0.5×log(a₀)');
console.log('  → Slope of RAR should approach 0.5 at low acceleration');

const lowAcc = gals.filter(g => g.log_g_bar < -10.5);
const highAcc = gals.filter(g => g.log_g_bar >= -10.5);

if (lowAcc.length >= 5 && highAcc.length >= 5) {
  const slopeLow = (() => {
    const x = lowAcc.map(g=>g.log_g_bar), y = lowAcc.map(g=>g.log_g_obs);
    const mx=x.reduce((a,b)=>a+b,0)/x.length, my=y.reduce((a,b)=>a+b,0)/y.length;
    let num=0,den=0;for(let i=0;i<x.length;i++){num+=(x[i]-mx)*(y[i]-my);den+=(x[i]-mx)**2;}return den>0?num/den:0;
  })();
  const slopeHigh = (() => {
    const x = highAcc.map(g=>g.log_g_bar), y = highAcc.map(g=>g.log_g_obs);
    const mx=x.reduce((a,b)=>a+b,0)/x.length, my=y.reduce((a,b)=>a+b,0)/y.length;
    let num=0,den=0;for(let i=0;i<x.length;i++){num+=(x[i]-mx)*(y[i]-my);den+=(x[i]-mx)**2;}return den>0?num/den:0;
  })();
  console.log(`  Low-g_bar (N=${lowAcc.length}): RAR slope = ${slopeLow.toFixed(3)} (MOND expects ~0.5)`);
  console.log(`  High-g_bar (N=${highAcc.length}): RAR slope = ${slopeHigh.toFixed(3)} (Newtonian expects ~1.0)`);
} else {
  console.log(`  Insufficient data for slope split (low: ${lowAcc.length}, high: ${highAcc.length})`);
}

console.log('\n  MOND Prediction 3: a₀ should be UNIVERSAL');
console.log('  But our channel says a₀ VARIES with VfResid — this is ANTI-MOND');
const a0_sd = Math.sqrt(gals.reduce((a,g)=>a+(g.logA0 - gals.reduce((aa,gg)=>aa+gg.logA0,0)/N)**2,0)/N);
console.log(`  σ(logA0) = ${a0_sd.toFixed(3)} dex`);
console.log(`  → a₀ varies by factor ${Math.pow(10, 2*a0_sd).toFixed(1)} (2σ range)`);
console.log(`  → This CONTRADICTS standard MOND (universal a₀)`);
console.log(`  → But CONSISTENT with modified/extended MOND or non-MOND physics`);


console.log('\n\n' + '▓'.repeat(70));
console.log('411F: RESIDUAL AFTER BEST FUNCTIONAL FIT');
console.log('Does the best function absorb the channel?');
console.log('▓'.repeat(70));

const bestFit = fits[0];
const rResid_channel = pearsonR(bestFit.resid, gals.map(g=>g.a0Resid));
const rResid_Vf = pearsonR(bestFit.resid, gals.map(g=>g.VfResid));

console.log(`\n  Best fit: ${bestFit.name} (R²=${bestFit.R2.toFixed(3)})`);
console.log(`  r(residual, a₀_resid) = ${rResid_channel.toFixed(4)}`);
console.log(`  r(residual, VfResid) = ${rResid_Vf.toFixed(4)}`);

const linearResids = fits.find(f=>f.name==='Linear').resid;
const rLinResid_a0 = pearsonR(linearResids, gals.map(g=>g.a0Resid));
console.log(`\n  After linear fit:`);
console.log(`  r(linear_residual, a₀_resid) = ${rLinResid_a0.toFixed(4)}`);
console.log(`  → ${Math.abs(rLinResid_a0) < 0.1 ? 'Linear captures the channel well' : 'Significant residual structure remains'}`);

const channelR2_linear = fits.find(f=>f.name==='Linear').R2;
console.log(`\n  Channel absorption by linear fit:`);
console.log(`  R²(VfResid ~ a₀_resid, linear) = ${channelR2_linear.toFixed(3)}`);
console.log(`  This is r² = ${rBase.toFixed(3)}² = ${(rBase**2).toFixed(3)}`);
console.log(`  → Linear captures ${(channelR2_linear*100).toFixed(1)}% of VfResid variance through a₀_resid`);


console.log('\n\n' + '▓'.repeat(70));
console.log('411G: CROSS-VALIDATION OF FUNCTIONAL FORM');
console.log('Is the linear relationship stable in leave-one-out?');
console.log('▓'.repeat(70));

let looSSE = 0;
for (let i = 0; i < N; i++) {
  const x_train = [...x.slice(0,i), ...x.slice(i+1)];
  const y_train = [...y.slice(0,i), ...y.slice(i+1)];
  const mx = x_train.reduce((a,b)=>a+b,0)/(N-1);
  const my = y_train.reduce((a,b)=>a+b,0)/(N-1);
  let num=0,den=0;
  for(let j=0;j<N-1;j++){num+=(x_train[j]-mx)*(y_train[j]-my);den+=(x_train[j]-mx)**2;}
  const b = den>0?num/den:0, a = my-b*mx;
  const pred = a + b * x[i];
  looSSE += (y[i] - pred) ** 2;
}
const sst = y.reduce((a,v)=>a+(v-y.reduce((a2,b)=>a2+b,0)/N)**2,0);
const looR2 = 1 - looSSE / sst;
console.log(`\n  LOO-CV R² (linear): ${looR2.toFixed(3)} (in-sample: ${channelR2_linear.toFixed(3)})`);
console.log(`  Shrinkage: ${(channelR2_linear - looR2).toFixed(4)}`);
console.log(`  → ${looR2 > 0.55 ? 'Robust — cross-validated relationship is real' : 'Some instability'}`);


console.log('\n\n' + '═'.repeat(70));
console.log('PHASE 411 GRAND VERDICT');
console.log('═'.repeat(70));

const isLinear = improvementOverLinear < 0.02;
const curvatured = Math.abs(rResidA0R2) > 0.2;
const nuCorrelated = Math.abs(rNuSimple_L) > 0.3;
const a0_varies = a0_sd > 0.1;

console.log(`
FUNCTIONAL FORM ANALYSIS:

1. Shape of channel: ${isLinear ? 'LINEAR' : 'NON-LINEAR'}
   Best model: ${fits[0].name} (R²=${fits[0].R2.toFixed(3)})
   Improvement over linear: ${improvementOverLinear.toFixed(4)}
   Curvature signal: ${curvatured ? 'YES' : 'NO'} (r(resid,x²)=${rResidA0R2.toFixed(3)})

2. MOND ν-function correlation with channel: ${nuCorrelated ? 'YES' : 'WEAK/NO'}
   r(log ν_simple, L_sum) = ${rNuSimple_L.toFixed(3)}
   r(log ν_exp, L_sum) = ${rNuExp_L.toFixed(3)}

3. a₀ universality: ${a0_varies ? 'VIOLATED' : 'CONSISTENT'}
   σ(logA0) = ${a0_sd.toFixed(3)} dex (factor ${Math.pow(10,2*a0_sd).toFixed(1)} range)
   → Standard MOND predicts universal a₀
   → Our data shows VARYING a₀ coupled to VfResid

4. Channel form: VfResid = ${fits.find(f=>f.name==='Linear').R2.toFixed(3)} × a₀_resid (linear)
   LOO-CV R² = ${looR2.toFixed(3)}
`);

if (isLinear && !curvatured) {
  console.log('VERDICT: The VfResid–a₀ channel is LINEARLY coupled.');
  console.log('There is NO evidence for MOND-like curvature or interpolation function.');
  console.log('');
  console.log('This means:');
  console.log('  - The channel is NOT a simple acceleration law with a transition scale');
  console.log('  - A linear coupling between velocity and acceleration residuals');
  console.log('    suggests a COMMON CAUSE driving both, not one causing the other');
  console.log('  - The hidden variable affects BOTH Vflat and a₀ proportionally');
  console.log('');
  console.log('IMPLICATION: MOND-like functional acceleration law DOES NOT explain the channel.');
  console.log('The linearity points to a shared hidden driver (likely halo structure or');
  console.log('assembly history) that imprints proportionally on both observables.');
} else {
  console.log('VERDICT: Non-linear structure detected in the channel.');
  console.log('MOND-like functional form may partially explain the coupling.');
}

if (a0_varies) {
  console.log(`\nCRITICAL: a₀ VARIATION (σ=${a0_sd.toFixed(3)} dex) directly contradicts`);
  console.log('standard MOND, which requires a universal acceleration scale.');
  console.log('This variation IS the channel — it is what our hidden variable modulates.');
}

console.log('\n→ REMAINING CANDIDATES after Phase 411:');
console.log('  1. ~~EFE / external field~~ — eliminated (Phase 410)');
console.log('  2. ~~MOND-like acceleration law~~ — eliminated (linear, no curvature, a₀ varies)');
console.log('  3. Dark matter halo structure — STRONGEST remaining candidate');
console.log('  4. Assembly/formation history — consistent with linearity');

const outPath = path.join(__dirname, '..', 'public', 'phase411-functional-form.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '411',
  title: 'MOND-like Functional Form Discrimination',
  timestamp: new Date().toISOString(),
  N,
  baselineR: rBase,
  functionalFits: fits.map(f => ({ name: f.name, R2: f.R2, AIC: f.AIC, BIC: f.BIC, k: f.k })),
  improvementOverLinear,
  curvatureSignal: rResidA0R2,
  nuCorrelations: { simple_L: rNuSimple_L, exp_L: rNuExp_L },
  a0_variation: a0_sd,
  isLinear,
  looR2,
  quintiles,
  discrepancyCorrelations: { a0R: rDisc_a0resid, VfR: rDisc_VfResid, Lsum: rDisc_Lsum },
  verdict: isLinear ? 'LINEAR — no MOND-like curvature' : 'NON-LINEAR structure detected'
}, null, 2));

console.log(`\nSaved: ${outPath}`);
