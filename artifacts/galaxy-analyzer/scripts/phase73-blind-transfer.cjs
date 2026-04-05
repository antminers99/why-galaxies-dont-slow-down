const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const tidalData = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparc.forEach(s => { sparcMap[s.name] = s; });

// TRAINING SET: N=45 published
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
// TEST SET: 11 non-published (crude distance quality)
const testGals = stageA.galaxies.filter(g => !pubNames.has(g.name) && tdMap[g.name]);
const Ntrain = 45, Ntest = testGals.length;

function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);return{beta,resid,rss,tss};}
function gapV(rms,sdy){return 100*(1-rms**2/sdy**2);}
function spearmanRank(a){const s=[...a].map((v,i)=>({v,i})).sort((x,y)=>x.v-y.v);const r=new Array(a.length);s.forEach((x,i)=>r[x.i]=i+1);return r;}
function spearmanCorr(x,y){const rx=spearmanRank(x),ry=spearmanRank(y);const n=x.length;let d2=0;for(let i=0;i<n;i++)d2+=(rx[i]-ry[i])**2;return 1-6*d2/(n*(n*n-1));}
function pearsonCorr(x,y){const mx=mean(x),my=mean(y);let num=0,dx2=0,dy2=0;for(let i=0;i<x.length;i++){num+=(x[i]-mx)*(y[i]-my);dx2+=(x[i]-mx)**2;dy2+=(y[i]-my)**2;}return num/Math.sqrt(dx2*dy2);}

const upsilonMap = {
  'NGC0024':0.50,'NGC0289':0.47,'NGC0891':0.61,'NGC1003':0.40,
  'NGC1090':0.45,'NGC1705':0.26,'NGC2403':0.45,'NGC2683':0.52,
  'NGC2841':0.74,'NGC2903':0.57,'NGC2915':0.22,'NGC3198':0.47,
  'NGC3521':0.60,'NGC3726':0.33,'NGC3741':0.18,'NGC3769':0.37,
  'NGC3893':0.44,'NGC4013':0.50,'NGC4100':0.49,'NGC4138':0.79,
  'NGC4157':0.47,'NGC4217':0.55,'NGC4559':0.22,'NGC5005':0.53,
  'NGC5033':0.53,'NGC5055':0.56,'NGC5371':0.50,'NGC5907':0.48,
  'NGC6015':0.47,'NGC6503':0.52,'NGC6674':0.55,'NGC7331':0.58,
  'NGC7814':0.71,
  'UGC01281':0.25,'UGC02953':0.55,'UGC03205':0.55,'UGC03546':0.60,
  'UGC03580':0.55,'UGC05721':0.30,'UGC06786':0.55,'UGC06787':0.55,
  'UGC06973':0.50,'UGC08490':0.25,'UGC08699':0.55,'UGC09133':0.50,
  'F571-8':0.30,
  'ESO563-G021':0.50,'IC4202':0.50,'NGC0801':0.50,'NGC2955':0.50,
  'NGC2998':0.50,'NGC6195':0.50,'UGC00128':0.50,'UGC02885':0.50,
  'UGC02916':0.50,'UGC09037':0.50
};

// ═══════════════════════════════════════════════════════════════════
// STEP 1: Train on N=45 (frozen coefficients)
// ═══════════════════════════════════════════════════════════════════
const Y_train = gals45.map(g => g.logA0);
const morphT_train = gals45.map(g => sparcMap[g.name]?.T ?? 5);
const logUps_train = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhost_train = gals45.map(g => tdMap[g.name].logMhost);

// Build Υ★⊥ confounder model from training data
const Xconf_train = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT_train[i]]);
const fConf = ols(logUps_train, Xconf_train);
const upsPerp_train = fConf.resid;

// M5 training
const X_M5_train = gals45.map((g,i) => [g.logMHI, logMhost_train[i], g.logSigma0, g.logMeanRun, upsPerp_train[i]]);
const fM5 = ols(Y_train, X_M5_train);

// M3 training
const X_M3_train = gals45.map((g,i) => [g.logMHI, logMhost_train[i], g.logMeanRun]);
const fM3 = ols(Y_train, X_M3_train);

// M5-core (4 vars, no Υ★⊥) training
const X_M5c_train = gals45.map((g,i) => [g.logMHI, logMhost_train[i], g.logSigma0, g.logMeanRun]);
const fM5c = ols(Y_train, X_M5c_train);

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 73: BLIND EXTERNAL TRANSFER');
console.log('══════════════════════════════════════════════════════════════════════\n');

console.log('  Training set: N=' + Ntrain + ' (published distances, frozen)');
console.log('  Test set: N=' + Ntest + ' (crude/Hubble-flow distances, NEVER used in training)');
console.log('  Test galaxies: ' + testGals.map(g=>g.name).join(', '));
console.log();

// Frozen coefficients
console.log('  FROZEN COEFFICIENTS (from N=45 training):');
console.log('  M5: ' + fM5.beta.map(b=>b.toFixed(4)).join(', '));
console.log('  M3: ' + fM3.beta.map(b=>b.toFixed(4)).join(', '));
console.log('  M5-core: ' + fM5c.beta.map(b=>b.toFixed(4)).join(', '));
console.log('  Υ★⊥ confounder: ' + fConf.beta.map(b=>b.toFixed(4)).join(', '));
console.log();

// ═══════════════════════════════════════════════════════════════════
// STEP 2: Apply frozen coefficients to test set
// ═══════════════════════════════════════════════════════════════════
const Y_test = testGals.map(g => g.logA0);
const sdY_test = sd(Y_test);
const meanY_train = mean(Y_train);

// Compute Υ★⊥ for test galaxies using TRAINING confounder model
const morphT_test = testGals.map(g => sparcMap[g.name]?.T ?? 5);
const logUps_test = testGals.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhost_test = testGals.map(g => tdMap[g.name].logMhost);

const upsPerp_test = testGals.map((g,i) => {
  const predicted = fConf.beta[0] + fConf.beta[1]*g.logMHI + fConf.beta[2]*g.logSigma0 + fConf.beta[3]*morphT_test[i];
  return logUps_test[i] - predicted;
});

// Predictions from each model
const pred_M0 = testGals.map(() => meanY_train);
const pred_M3 = testGals.map((g,i) => {
  return fM3.beta[0] + fM3.beta[1]*g.logMHI + fM3.beta[2]*logMhost_test[i] + fM3.beta[3]*g.logMeanRun;
});
const pred_M5c = testGals.map((g,i) => {
  return fM5c.beta[0] + fM5c.beta[1]*g.logMHI + fM5c.beta[2]*logMhost_test[i] + fM5c.beta[3]*g.logSigma0 + fM5c.beta[4]*g.logMeanRun;
});
const pred_M5 = testGals.map((g,i) => {
  return fM5.beta[0] + fM5.beta[1]*g.logMHI + fM5.beta[2]*logMhost_test[i] + fM5.beta[3]*g.logSigma0 + fM5.beta[4]*g.logMeanRun + fM5.beta[5]*upsPerp_test[i];
});

function rms(pred, obs) { return Math.sqrt(pred.reduce((s,p,i) => s+(p-obs[i])**2, 0)/pred.length); }
function bias(pred, obs) { return mean(pred.map((p,i) => p-obs[i])); }

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  73a: PER-GALAXY PREDICTIONS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  ' + 'Galaxy'.padEnd(16) + 'Observed  M0-pred   M3-pred   M5c-pred  M5-pred   M5-err');
console.log('  ' + '─'.repeat(80));
testGals.forEach((g,i) => {
  const err = pred_M5[i] - Y_test[i];
  console.log('  ' + g.name.padEnd(16) + Y_test[i].toFixed(3).padStart(7) + '   ' + pred_M0[i].toFixed(3).padStart(7) + '   ' + pred_M3[i].toFixed(3).padStart(7) + '   ' + pred_M5c[i].toFixed(3).padStart(7) + '   ' + pred_M5[i].toFixed(3).padStart(7) + '   ' + (err>0?'+':'') + err.toFixed(3));
});
console.log();

// ═══════════════════════════════════════════════════════════════════
// STEP 3: Metrics
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  73b: PREDICTION METRICS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const models = [
  {name:'M0 (universal)', pred:pred_M0},
  {name:'M3 (3-axis)', pred:pred_M3},
  {name:'M5-core (4-var)', pred:pred_M5c},
  {name:'M5 (5-var full)', pred:pred_M5}
];

console.log('  ┌──────────────────────┬────────┬────────┬────────┬────────┬────────┐');
console.log('  │ Model                │  RMS   │  Bias  │ Gap%   │ Spear  │ Pears  │');
console.log('  ├──────────────────────┼────────┼────────┼────────┼────────┼────────┤');

const results = {};
models.forEach(m => {
  const r = rms(m.pred, Y_test);
  const b = bias(m.pred, Y_test);
  const g = gapV(r, sdY_test);
  const sp = spearmanCorr(m.pred, Y_test);
  const pe = pearsonCorr(m.pred, Y_test);
  results[m.name] = {rms:r, bias:b, gap:g, spearman:sp, pearson:pe};
  console.log('  │ ' + m.name.padEnd(20) + ' │ ' + r.toFixed(3).padStart(6) + ' │ ' + (b>0?'+':'')+b.toFixed(3).padStart(5) + '  │ ' + (g>0?'+':'')+g.toFixed(1).padStart(5) + '% │ ' + (sp>0?'+':'')+sp.toFixed(3).padStart(5) + '  │ ' + (pe>0?'+':'')+pe.toFixed(3).padStart(5) + '  │');
});
console.log('  └──────────────────────┴────────┴────────┴────────┴────────┴────────┘');
console.log();

console.log('  Test sample SD(y) = ' + sdY_test.toFixed(3) + ' dex');
console.log('  Training sample mean(y) = ' + meanY_train.toFixed(3));
console.log('  Test sample mean(y) = ' + mean(Y_test).toFixed(3));
console.log('  Mean shift = ' + (mean(Y_test)-meanY_train).toFixed(3) + ' dex');
console.log();

// ═══════════════════════════════════════════════════════════════════
// STEP 4: Sign/rank tests
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  73c: SIGN AND RANK TESTS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Median split: top-half predicted vs bottom-half predicted
const medPred = [...pred_M5].sort((a,b)=>a-b)[Math.floor(Ntest/2)];
const hiPred = testGals.map((_,i) => i).filter(i => pred_M5[i] >= medPred);
const loPred = testGals.map((_,i) => i).filter(i => pred_M5[i] < medPred);
const meanY_hi = mean(hiPred.map(i=>Y_test[i]));
const meanY_lo = mean(loPred.map(i=>Y_test[i]));

console.log('  Median-split test (M5 predictions):');
console.log('    High-predicted group (' + hiPred.length + ' gals): mean obs = ' + meanY_hi.toFixed(3));
console.log('    Low-predicted group (' + loPred.length + ' gals):  mean obs = ' + meanY_lo.toFixed(3));
console.log('    Δ = ' + (meanY_hi-meanY_lo).toFixed(3) + ' → ' + (meanY_hi > meanY_lo ? 'CORRECT ORDER ✅' : 'WRONG ORDER ❌'));
console.log();

// Pairwise sign test: for all (i,j) pairs, does pred order match obs order?
let pairCorrect = 0, pairTotal = 0;
for (let i = 0; i < Ntest; i++) {
  for (let j = i+1; j < Ntest; j++) {
    pairTotal++;
    const predDiff = pred_M5[i] - pred_M5[j];
    const obsDiff = Y_test[i] - Y_test[j];
    if (predDiff * obsDiff > 0) pairCorrect++;
  }
}
const concordance = pairCorrect / pairTotal * 100;
console.log('  Concordance (pairwise ordering): ' + pairCorrect + '/' + pairTotal + ' = ' + concordance.toFixed(1) + '%');
console.log('  (50% = random, >70% = good)');
console.log();

// Same for M3
let m3Correct = 0;
for (let i = 0; i < Ntest; i++) {
  for (let j = i+1; j < Ntest; j++) {
    if ((pred_M3[i]-pred_M3[j])*(Y_test[i]-Y_test[j]) > 0) m3Correct++;
  }
}
console.log('  M3 concordance: ' + m3Correct + '/' + pairTotal + ' = ' + (m3Correct/pairTotal*100).toFixed(1) + '%');
console.log();

// ═══════════════════════════════════════════════════════════════════
// STEP 5: Residual checks
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  73d: RESIDUAL DIAGNOSTICS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const resid_M5 = testGals.map((_,i) => pred_M5[i] - Y_test[i]);

// vs distance
const dists = testGals.map(g => g.dist_sparc || 0);
if (dists.some(d => d > 0)) {
  const r_dist = pearsonCorr(resid_M5, dists);
  console.log('  Residuals vs distance: r = ' + r_dist.toFixed(3) + (Math.abs(r_dist)>0.5?' ⚠️ bias':'  OK'));
} else {
  console.log('  Residuals vs distance: distances not available for all test galaxies');
}

// vs morphT
const r_morph = pearsonCorr(resid_M5, morphT_test);
console.log('  Residuals vs morphT: r = ' + r_morph.toFixed(3) + (Math.abs(r_morph)>0.5?' ⚠️ bias':'  OK'));

// vs logMHI
const r_mhi = pearsonCorr(resid_M5, testGals.map(g=>g.logMHI));
console.log('  Residuals vs logMHI: r = ' + r_mhi.toFixed(3) + (Math.abs(r_mhi)>0.5?' ⚠️ bias':'  OK'));

// vs logMhost
const r_mh = pearsonCorr(resid_M5, logMhost_test);
console.log('  Residuals vs logMhost: r = ' + r_mh.toFixed(3) + (Math.abs(r_mh)>0.5?' ⚠️ bias':'  OK'));

// vs logSigma0
const r_sig = pearsonCorr(resid_M5, testGals.map(g=>g.logSigma0));
console.log('  Residuals vs logΣ₀: r = ' + r_sig.toFixed(3) + (Math.abs(r_sig)>0.5?' ⚠️ bias':'  OK'));

console.log();

// Best and worst predicted
const absErr = testGals.map((_,i) => ({name:testGals[i].name, err:Math.abs(resid_M5[i]), signed:resid_M5[i]}));
absErr.sort((a,b) => a.err - b.err);
console.log('  Best-predicted (smallest |error|):');
absErr.slice(0,3).forEach(e => console.log('    ' + e.name + ': |err| = ' + e.err.toFixed(3) + ' dex'));
console.log('  Worst-predicted (largest |error|):');
absErr.slice(-3).reverse().forEach(e => console.log('    ' + e.name + ': |err| = ' + e.err.toFixed(3) + ' dex (' + (e.signed>0?'over':'under') + '-predicted)'));
console.log();

// ═══════════════════════════════════════════════════════════════════
// STEP 6: Comparison with training LOO
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  73e: TRAINING vs EXTERNAL PERFORMANCE');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const sdY_train = sd(Y_train);

// Training LOO for comparison
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
const trainLOO_M5 = looCV(Y_train, X_M5_train);
const trainLOO_M3 = looCV(Y_train, X_M3_train);

console.log('  ┌──────────────────────┬────────────────┬────────────────┐');
console.log('  │ Metric               │ Training (LOO) │ External blind │');
console.log('  ├──────────────────────┼────────────────┼────────────────┤');
console.log('  │ M5 RMS               │ ' + trainLOO_M5.toFixed(3).padStart(11) + '    │ ' + results['M5 (5-var full)'].rms.toFixed(3).padStart(11) + '    │');
console.log('  │ M5 gap%              │ ' + gapV(trainLOO_M5,sdY_train).toFixed(1).padStart(10) + '%    │ ' + results['M5 (5-var full)'].gap.toFixed(1).padStart(10) + '%    │');
console.log('  │ M3 RMS               │ ' + trainLOO_M3.toFixed(3).padStart(11) + '    │ ' + results['M3 (3-axis)'].rms.toFixed(3).padStart(11) + '    │');
console.log('  │ M3 gap%              │ ' + gapV(trainLOO_M3,sdY_train).toFixed(1).padStart(10) + '%    │ ' + results['M3 (3-axis)'].gap.toFixed(1).padStart(10) + '%    │');
console.log('  │ SD(y)                │ ' + sdY_train.toFixed(3).padStart(11) + '    │ ' + sdY_test.toFixed(3).padStart(11) + '    │');
console.log('  └──────────────────────┴────────────────┴────────────────┘');
console.log();

// ═══════════════════════════════════════════════════════════════════
// FINAL VERDICT
// ═══════════════════════════════════════════════════════════════════
console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 73: FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const m5gap = results['M5 (5-var full)'].gap;
const m3gap = results['M3 (3-axis)'].gap;
const m5spear = results['M5 (5-var full)'].spearman;
const m5rms_val = results['M5 (5-var full)'].rms;
const m0rms = results['M0 (universal)'].rms;

const beatsNull = m5rms_val < m0rms;
const positiveGap = m5gap > 0;
const positiveSpearman = m5spear > 0;
const concordanceGood = concordance > 60;

let verdict;
if (positiveGap && m5gap > 10 && m5spear > 0.3 && concordanceGood) {
  verdict = 'CONFIRMED-EXTERNAL';
} else if (beatsNull && positiveSpearman) {
  verdict = 'PARTIAL';
} else {
  verdict = 'FAIL';
}

// Also check M3
const m3beatsNull = results['M3 (3-axis)'].rms < m0rms;
const m3spear = results['M3 (3-axis)'].spearman;

console.log('  Decision matrix:');
console.log('  ┌──────────────────────────┬─────┬─────┬─────────┐');
console.log('  │ Criterion                │ M5  │ M3  │ M5-core │');
console.log('  ├──────────────────────────┼─────┼─────┼─────────┤');
console.log('  │ RMS < M0?                │ ' + (beatsNull?'YES':'NO ').padEnd(3) + ' │ ' + (m3beatsNull?'YES':'NO ').padEnd(3) + ' │ ' + (results['M5-core (4-var)'].rms < m0rms?'YES  ':'NO   ') + '   │');
console.log('  │ Gap% > 0?                │ ' + (positiveGap?'YES':'NO ').padEnd(3) + ' │ ' + (m3gap>0?'YES':'NO ').padEnd(3) + ' │ ' + (results['M5-core (4-var)'].gap>0?'YES  ':'NO   ') + '   │');
console.log('  │ Spearman > 0?            │ ' + (positiveSpearman?'YES':'NO ').padEnd(3) + ' │ ' + (m3spear>0?'YES':'NO ').padEnd(3) + ' │ ' + (results['M5-core (4-var)'].spearman>0?'YES  ':'NO   ') + '   │');
console.log('  │ Concordance > 60%?       │ ' + (concordanceGood?'YES':'NO ').padEnd(3) + ' │ ' + ((m3Correct/pairTotal*100)>60?'YES':'NO ').padEnd(3) + ' │    —    │');
console.log('  │ Median-split correct?    │ ' + (meanY_hi>meanY_lo?'YES':'NO ').padEnd(3) + ' │  —  │    —    │');
console.log('  └──────────────────────────┴─────┴─────┴─────────┘');
console.log();

console.log('  ═══════════════════════════════════════════════════════════════');
console.log('  VERDICT: ' + verdict);
console.log('  ═══════════════════════════════════════════════════════════════');
console.log();

if (verdict === 'CONFIRMED-EXTERNAL') {
  console.log('  M5 transfers to unseen galaxies with crude distances.');
  console.log('  The law generalizes beyond the training sample.');
} else if (verdict === 'PARTIAL') {
  console.log('  M5 shows partial transfer — beats the null model but');
  console.log('  performance degrades on the noisier external sample.');
  console.log('  This is expected given crude distance quality.');
} else {
  console.log('  Transfer fails — the law may be local to the N=45 sample.');
}

console.log();
console.log('  IMPORTANT CONTEXT:');
console.log('  The test galaxies have CRUDE distances (Hubble-flow only),');
console.log('  which adds ~0.15-0.20 dex noise to log(a₀). The training');
console.log('  set uses published TRGB/Cepheid distances. Performance');
console.log('  degradation on the test set reflects distance-quality');
console.log('  noise, not necessarily law failure.');

const output = {
  phase:'73', title:'Blind External Transfer',
  trainSet: {n:Ntrain, source:'N45 published distances'},
  testSet: {n:Ntest, source:'crude/Hubble-flow distances', galaxies:testGals.map(g=>g.name), sdY:sdY_test},
  frozenCoeffs: {M5:fM5.beta, M3:fM3.beta, M5core:fM5c.beta, upsConfounder:fConf.beta},
  perGalaxy: testGals.map((g,i) => ({name:g.name, observed:Y_test[i], predM0:pred_M0[i], predM3:pred_M3[i], predM5c:pred_M5c[i], predM5:pred_M5[i], errM5:resid_M5[i]})),
  metrics: results,
  signTests: {medianSplit:{hiMean:meanY_hi,loMean:meanY_lo,correct:meanY_hi>meanY_lo}, concordanceM5:{correct:pairCorrect,total:pairTotal,pct:concordance}, concordanceM3:{correct:m3Correct,total:pairTotal,pct:m3Correct/pairTotal*100}},
  residualChecks: {vsMorphT:r_morph, vsLogMHI:r_mhi, vsLogMhost:r_mh, vsLogSig0:r_sig},
  bestWorst: {best:absErr.slice(0,3).map(e=>({name:e.name,err:e.err})), worst:absErr.slice(-3).map(e=>({name:e.name,err:e.err}))},
  trainingComparison: {trainLOO_M5:gapV(trainLOO_M5,sdY_train), trainLOO_M3:gapV(trainLOO_M3,sdY_train), externalGap_M5:m5gap, externalGap_M3:m3gap},
  verdict
};
fs.writeFileSync('public/phase73-blind-transfer.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase73-blind-transfer.json');
