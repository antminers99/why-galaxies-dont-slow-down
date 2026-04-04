const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const tidalData = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const tdMap = {};
tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {};
sparc.forEach(s => { sparcMap[s.name] = s; });

// All 56 galaxies from Stage A
const all56 = stageA.galaxies;
const gals45 = all56.filter(g => pubNames.has(g.name));
const gals11 = all56.filter(g => !pubNames.has(g.name));

function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function pearson(x,y){const n=x.length;const mx=mean(x),my=mean(y);let sxy=0,sxx=0,syy=0;for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;syy+=(y[i]-my)**2;}if(sxx<1e-15||syy<1e-15)return{r:NaN,t:NaN};const r=sxy/Math.sqrt(sxx*syy);return{r,t:r*Math.sqrt((n-2)/(1-r*r)),n};}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se:Math.sqrt(rss/(n-p)),resid,rss,tss,n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}

// Upsilon data
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
  'F571-8':0.30
};

// Also need Upsilon for the 11 holdout galaxies
// These don't have logMhost, but we can test B' (5-var model without Υ★⊥)
// For the 11 holdout: assign estimated Upsilon by morphology
const holdoutUps = {};
gals11.forEach(g => {
  const T = sparcMap[g.name]?.T || 5;
  if (T <= 2) holdoutUps[g.name] = 0.65;
  else if (T <= 4) holdoutUps[g.name] = 0.50;
  else if (T <= 6) holdoutUps[g.name] = 0.40;
  else if (T <= 8) holdoutUps[g.name] = 0.30;
  else holdoutUps[g.name] = 0.25;
  if (upsilonMap[g.name]) holdoutUps[g.name] = upsilonMap[g.name];
});

const morphT45 = gals45.map(g => sparcMap[g.name]?.T || 5);
const logUps45 = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhost45 = gals45.map(g => tdMap[g.name].logMhost);
const Y45 = gals45.map(g => g.logA0);
const sdY45 = sd(Y45);
function gapV(rms, sdy) { return 100*(1-rms**2/sdy**2); }

// Compute Υ★⊥ confounding model on training (N=45)
const Xconf = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT45[i]]);
const fConf = ols(logUps45, Xconf);
const upsPerp45 = fConf.resid;

console.log('======================================================================');
console.log('  PHASE 61: EXTERNAL TRANSFER TEST FOR B″');
console.log('  Training: N=45 (logMhost published subsample)');
console.log('  Holdout: N=11 (remaining SPARC Stage-A galaxies)');
console.log('======================================================================\n');

// ═══════════════════════════════════════════════
// TEST 1: Train on N=45, predict holdout N=11
// ═══════════════════════════════════════════════
// Problem: the 11 holdout galaxies DON'T have logMhost
// So we test the REDUCED models that don't require logMhost:
// - M_reduced: logMHI + wig + logSigma0 + logMeanRun (4 vars, no logMhost, no Υ★⊥)
// This tests whether the CORE pattern transfers

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1a: Train/Test split — Core model (no logMhost)');
console.log('  Train on N=45, Test on N=11 holdout');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Train core model on N=45
const X_core_train = gals45.map(g => [g.logMHI, g.rcWiggliness, g.logSigma0, g.logMeanRun]);
const fCore = ols(Y45, X_core_train);

console.log('  Core model trained on N=45:');
['intercept','logMHI','rcWiggliness','logSigma0','logMeanRun'].forEach((v,j) => {
  console.log('    ' + v.padEnd(16) + ' β=' + fCore.beta[j].toFixed(4));
});
console.log('    R²adj = ' + fCore.r2adj.toFixed(4));
console.log();

// Predict holdout
const Y11 = gals11.map(g => g.logA0);
const sdY11 = sd(Y11);
const X_core_test = gals11.map(g => [g.logMHI, g.rcWiggliness, g.logSigma0, g.logMeanRun]);
const pred11 = X_core_test.map(x => {
  const xi = [1, ...x];
  return xi.reduce((s,v,j) => s + v*fCore.beta[j], 0);
});
const errors11 = Y11.map((y,i) => y - pred11[i]);
const rmsTest = Math.sqrt(errors11.reduce((s,e)=>s+e*e,0)/errors11.length);
const biasTest = mean(errors11);

console.log('  Holdout prediction (N=11):');
console.log('    RMS error = ' + rmsTest.toFixed(4) + ' dex');
console.log('    Bias = ' + biasTest.toFixed(4) + ' dex');
console.log('    Holdout sd(Y) = ' + sdY11.toFixed(4) + ' dex');
console.log('    Gap% on holdout = ' + gapV(rmsTest, sdY11).toFixed(1) + '%');
console.log('    Training LOO = ' + gapV(looCV(Y45, X_core_train), sdY45).toFixed(1) + '%');
console.log();

console.log('  Per-galaxy predictions:');
gals11.forEach((g,i) => {
  const err = errors11[i];
  const dir = err > 0 ? '+' : '';
  console.log('    ' + g.name.padEnd(15) + ' obs=' + g.logA0.toFixed(3) + ' pred=' + pred11[i].toFixed(3) + ' err=' + dir + err.toFixed(3));
});
console.log();

// ═══════════════════════════════════════════════
// TEST 1b: Also train old Baseline A (with envCode instead of logMhost)
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1b: Train N=45, Test N=11 — with envCode (available for all)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// envCode is available for all 56 galaxies from the original baseline
const phase56 = JSON.parse(fs.readFileSync('public/phase56-frozen-baselines.json','utf8'));
const envMap = {};
if (phase56.baselineA && phase56.baselineA.galaxyData) {
  phase56.baselineA.galaxyData.forEach(g => { envMap[g.name] = g.envCode; });
}
// Fallback: compute from sparc
all56.forEach(g => {
  if (envMap[g.name] === undefined) {
    envMap[g.name] = (sparcMap[g.name]?.env === 'F' || sparcMap[g.name]?.env === 'field') ? 0 : 1;
  }
});

// Train with envCode on N=45
const X_env_train = gals45.map(g => [g.logMHI, g.rcWiggliness, envMap[g.name], g.logSigma0, g.logMeanRun]);
const fEnv = ols(Y45, X_env_train);

const X_env_test = gals11.map(g => [g.logMHI, g.rcWiggliness, envMap[g.name], g.logSigma0, g.logMeanRun]);
const predEnv11 = X_env_test.map(x => {
  const xi = [1, ...x];
  return xi.reduce((s,v,j) => s + v*fEnv.beta[j], 0);
});
const errorsEnv11 = Y11.map((y,i) => y - predEnv11[i]);
const rmsEnvTest = Math.sqrt(errorsEnv11.reduce((s,e)=>s+e*e,0)/errorsEnv11.length);

console.log('  With envCode, RMS holdout = ' + rmsEnvTest.toFixed(4) + ' dex');
console.log('  Gap% on holdout = ' + gapV(rmsEnvTest, sdY11).toFixed(1) + '%');
console.log();

// ═══════════════════════════════════════════════
// TEST 2: Strict train/test split WITHIN N=45
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: 50/50 train/test split WITHIN N=45 (100 reps)');
console.log('  Testing B″ (full 6-var model)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const X_Bdp = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost45[i], g.logSigma0, g.logMeanRun, upsPerp45[i]]);
const X_Bp = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost45[i], g.logSigma0, g.logMeanRun]);
const X_Ap = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost45[i], g.logSigma0]);

const nReps = 200;
const splitResults = {M1c:[], M2:[], M3:[]};

for (let rep = 0; rep < nReps; rep++) {
  // Random 50/50 split
  const idx = Array.from({length:45},(_,i)=>i);
  for(let i=44;i>0;i--){const j=Math.floor(Math.random()*(i+1));[idx[i],idx[j]]=[idx[j],idx[i]];}
  const trainIdx = idx.slice(0,23);
  const testIdx = idx.slice(23);

  const Ytrain = trainIdx.map(i=>Y45[i]);
  const Ytest = testIdx.map(i=>Y45[i]);
  const sdTest = sd(Ytest);

  // M1c (A')
  try {
    const Xt = trainIdx.map(i=>X_Ap[i]);
    const f = ols(Ytrain, Xt);
    const errs = testIdx.map(i => {
      const xi = [1,...X_Ap[i]];
      return Y45[i] - xi.reduce((s,x,j)=>s+x*f.beta[j],0);
    });
    const rms = Math.sqrt(errs.reduce((s,e)=>s+e*e,0)/errs.length);
    splitResults.M1c.push(gapV(rms, sdTest));
  } catch(e) { splitResults.M1c.push(NaN); }

  // M2 (B')
  try {
    const Xt = trainIdx.map(i=>X_Bp[i]);
    const f = ols(Ytrain, Xt);
    const errs = testIdx.map(i => {
      const xi = [1,...X_Bp[i]];
      return Y45[i] - xi.reduce((s,x,j)=>s+x*f.beta[j],0);
    });
    const rms = Math.sqrt(errs.reduce((s,e)=>s+e*e,0)/errs.length);
    splitResults.M2.push(gapV(rms, sdTest));
  } catch(e) { splitResults.M2.push(NaN); }

  // M3 (B″)
  try {
    const Xt = trainIdx.map(i=>X_Bdp[i]);
    const f = ols(Ytrain, Xt);
    const errs = testIdx.map(i => {
      const xi = [1,...X_Bdp[i]];
      return Y45[i] - xi.reduce((s,x,j)=>s+x*f.beta[j],0);
    });
    const rms = Math.sqrt(errs.reduce((s,e)=>s+e*e,0)/errs.length);
    splitResults.M3.push(gapV(rms, sdTest));
  } catch(e) { splitResults.M3.push(NaN); }
}

function validMean(a){const v=a.filter(x=>!isNaN(x));return v.length?mean(v):NaN;}
function validSD(a){const v=a.filter(x=>!isNaN(x));return v.length>1?sd(v):NaN;}
function pctile(a,p){const v=a.filter(x=>!isNaN(x)).sort((a,b)=>a-b);return v[Math.floor(p*v.length)];}

console.log('  50/50 split test gap% (200 reps):');
console.log('  ┌──────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Model    │ Mean     │ SD       │ p25      │ p75      │');
console.log('  ├──────────┼──────────┼──────────┼──────────┼──────────┤');
['M1c','M2','M3'].forEach(m => {
  const a = splitResults[m];
  const marker = m === 'M3' ? ' ◄' : '';
  console.log('  │ ' + (m+marker).padEnd(8) + ' │ ' + validMean(a).toFixed(1).padStart(5) + '%   │ ' + validSD(a).toFixed(1).padStart(5) + '%   │ ' + pctile(a,0.25).toFixed(1).padStart(5) + '%   │ ' + pctile(a,0.75).toFixed(1).padStart(5) + '%   │');
});
console.log('  └──────────┴──────────┴──────────┴──────────┴──────────┘');
console.log();

// How often does M3 beat M2?
let m3beats = 0, valid = 0;
for (let i = 0; i < nReps; i++) {
  if (!isNaN(splitResults.M3[i]) && !isNaN(splitResults.M2[i])) {
    valid++;
    if (splitResults.M3[i] > splitResults.M2[i]) m3beats++;
  }
}
console.log('  M3 beats M2 in ' + m3beats + '/' + valid + ' splits (' + (100*m3beats/valid).toFixed(1) + '%)');
console.log();

// ═══════════════════════════════════════════════
// TEST 3: Morphology-stratified transfer
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: Morphology-stratified LOO');
console.log('  Train on one morph type range, predict another');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Split: early (T<=4) vs late (T>4)
const earlyIdx = [], lateIdx = [];
gals45.forEach((g,i) => {
  const T = sparcMap[g.name]?.T || 5;
  if (T <= 4) earlyIdx.push(i); else lateIdx.push(i);
});

console.log('  Early-type (T<=4): N=' + earlyIdx.length);
console.log('  Late-type (T>4): N=' + lateIdx.length);
console.log();

// Train on early, predict late (using B' since it needs less data)
function trainPredict(trainI, testI, X, Y, label) {
  const Yt = trainI.map(i=>Y[i]), Xt = trainI.map(i=>X[i]);
  const f = ols(Yt, Xt);
  const errs = testI.map(i => {
    const xi = [1,...X[i]];
    return Y[i] - xi.reduce((s,x,j)=>s+x*f.beta[j],0);
  });
  const rms = Math.sqrt(errs.reduce((s,e)=>s+e*e,0)/errs.length);
  const sdTest = sd(testI.map(i=>Y[i]));
  console.log('  ' + label + ':');
  console.log('    RMS = ' + rms.toFixed(4) + ' dex, test sd=' + sdTest.toFixed(4) + ', gap%=' + gapV(rms, sdTest).toFixed(1) + '%');
  return {rms, gap: gapV(rms, sdTest)};
}

trainPredict(earlyIdx, lateIdx, X_Bp, Y45, 'Train early → predict late (B\x27)');
trainPredict(lateIdx, earlyIdx, X_Bp, Y45, 'Train late → predict early (B\x27)');
trainPredict(earlyIdx, lateIdx, X_Bdp, Y45, 'Train early → predict late (B″)');
trainPredict(lateIdx, earlyIdx, X_Bdp, Y45, 'Train late → predict early (B″)');
console.log();

// ═══════════════════════════════════════════════
// TEST 4: Coefficient stability across subsamples
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: Coefficient stability (bootstrap from 200 50/50 splits)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const varNames = ['intercept','logMHI','rcWig','logMhost','logΣ₀','logMR','Υ★⊥'];
const coefSamples = varNames.map(()=>[]);

for (let rep = 0; rep < 200; rep++) {
  const idx = Array.from({length:45},(_,i)=>i);
  for(let i=44;i>0;i--){const j=Math.floor(Math.random()*(i+1));[idx[i],idx[j]]=[idx[j],idx[i]];}
  const trainIdx = idx.slice(0,23);
  try {
    const f = ols(trainIdx.map(i=>Y45[i]), trainIdx.map(i=>X_Bdp[i]));
    f.beta.forEach((b,j) => coefSamples[j].push(b));
  } catch(e) {}
}

const fullFit = ols(Y45, X_Bdp);
console.log('  ┌──────────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Variable     │ Full β   │ Split μ  │ Split SD │ Stable?  │');
console.log('  ├──────────────┼──────────┼──────────┼──────────┼──────────┤');
varNames.forEach((v,j) => {
  const fullB = fullFit.beta[j];
  const splitMu = validMean(coefSamples[j]);
  const splitSd = validSD(coefSamples[j]);
  const stable = Math.abs(fullB - splitMu) < 2*splitSd ? '✅' : '⚠️';
  console.log('  │ ' + v.padEnd(12) + ' │ ' + fullB.toFixed(4).padStart(8) + ' │ ' + splitMu.toFixed(4).padStart(8) + ' │ ' + splitSd.toFixed(4).padStart(8) + ' │ ' + stable + '       │');
});
console.log('  └──────────────┴──────────┴──────────┴──────────┴──────────┘');
console.log();

// ═══════════════════════════════════════════════
// VERDICT
// ═══════════════════════════════════════════════
console.log('======================================================================');
console.log('  VERDICT');
console.log('======================================================================\n');

const holdoutGap = gapV(rmsTest, sdY11);
const splitMean3 = validMean(splitResults.M3);
const m3beatsPct = 100*m3beats/valid;

const test1pass = holdoutGap > 0;
const test2pass = splitMean3 > 30;
const test3pass = m3beatsPct > 50;

console.log('  Test 1 (holdout N=11 gap% > 0): ' + (test1pass ? 'PASS' : 'FAIL') + ' (' + holdoutGap.toFixed(1) + '%)');
console.log('  Test 2 (50/50 split mean gap% > 30): ' + (test2pass ? 'PASS' : 'FAIL') + ' (' + splitMean3.toFixed(1) + '%)');
console.log('  Test 3 (M3 beats M2 > 50% of splits): ' + (test3pass ? 'PASS' : 'FAIL') + ' (' + m3beatsPct.toFixed(1) + '%)');
console.log();

const passCount = [test1pass, test2pass, test3pass].filter(Boolean).length;
let verdict;
if (passCount === 3) verdict = 'TRANSFER CONFIRMED';
else if (passCount >= 2) verdict = 'TRANSFER PARTIAL';
else verdict = 'TRANSFER FAILS';

console.log('  ═══════════════════════════════════════════════════════════');
console.log('  ' + verdict + ' (' + passCount + '/3)');
console.log('  ═══════════════════════════════════════════════════════════');

const output = {
  phase: '61',
  title: 'External Transfer Test for B″',
  holdout: {
    n: 11, rms: rmsTest, bias: biasTest, sdY: sdY11, gap: holdoutGap,
    galaxies: gals11.map((g,i) => ({name:g.name, obs:g.logA0, pred:pred11[i], err:errors11[i]}))
  },
  splitTest: {
    nReps, splitSize: '23/22',
    M1c: {mean:validMean(splitResults.M1c), sd:validSD(splitResults.M1c)},
    M2: {mean:validMean(splitResults.M2), sd:validSD(splitResults.M2)},
    M3: {mean:validMean(splitResults.M3), sd:validSD(splitResults.M3)},
    m3BeatsM2Pct: m3beatsPct
  },
  coeffStability: varNames.map((v,j) => ({
    name:v, fullBeta:fullFit.beta[j],
    splitMean:validMean(coefSamples[j]), splitSD:validSD(coefSamples[j])
  })),
  verdict, passCount
};
fs.writeFileSync('public/phase61-external-transfer.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase61-external-transfer.json');
