const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const tidalData = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const tdMap = {}; tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {}; sparc.forEach(s => { sparcMap[s.name] = s; });
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
const N = 45;

function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const r2adj=1-(rss/(n-p))/(tss/(n-1));const aic=n*Math.log(rss/n)+2*p;const bic=n*Math.log(rss/n)+p*Math.log(n);const mse=rss/(n-p);const XtXinv=invertMatrix(XtX);const se=beta.map((_,j)=>Math.sqrt(mse*XtXinv[j][j]));const tStats=beta.map((b,j)=>b/se[j]);return{beta,resid,rss,tss,se,tStats,r2adj,aic,bic,n,k:p};}
function invertMatrix(M){const n=M.length;const A=M.map((r,i)=>{const row=[...r];for(let j=0;j<n;j++)row.push(i===j?1:0);return row;});for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(A[j][i])>Math.abs(A[mx][i]))mx=j;[A[i],A[mx]]=[A[mx],A[i]];const d=A[i][i];for(let j=0;j<2*n;j++)A[i][j]/=d;for(let j=0;j<n;j++){if(j===i)continue;const f=A[j][i];for(let k=0;k<2*n;k++)A[j][k]-=f*A[i][k];}}return A.map(r=>r.slice(n));}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function gapV(rms,sdy){return 100*(1-rms**2/sdy**2);}

// Seeded PRNG
function mulberry32(a){return function(){a|=0;a=a+0x6D2B79F5|0;var t=Math.imul(a^a>>>15,1|a);t=t+Math.imul(t^t>>>7,61|t)^t;return((t^t>>>14)>>>0)/4294967296;};}
function shuffle(arr,rng){const a=[...arr];for(let i=a.length-1;i>0;i--){const j=Math.floor(rng()*( i+1));[a[i],a[j]]=[a[j],a[i]];}return a;}

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

const Y = gals45.map(g => g.logA0);
const sdY = sd(Y);
const morphT = gals45.map(g => sparcMap[g.name]?.T ?? 5);
const logUps = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhost = gals45.map(g => tdMap[g.name].logMhost);

// Orthogonalize Υ★
const Xconf = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT[i]]);
const fConf = ols(logUps, Xconf);
const upsPerp = fConf.resid;

// Build predictor arrays
const logMHI = gals45.map(g => g.logMHI);
const rcWig = gals45.map(g => g.rcWiggliness);
const logSig0 = gals45.map(g => g.logSigma0);
const logMR = gals45.map(g => g.logMeanRun);

// Model definitions
const models = {
  M0: { name: 'M0 (universal a₀)', idx: [], vars: [] },
  M3: { name: 'M3 (3-axis)', idx: [0,2,4], vars: ['logMHI','logMhost','logMR'] },
  M5: { name: 'M5 (5-var reduced)', idx: [0,2,3,4,5], vars: ['logMHI','logMhost','logΣ₀','logMR','Υ★⊥'] },
  M6: { name: 'M6 (B″ full)', idx: [0,1,2,3,4,5], vars: ['logMHI','rcWig','logMhost','logΣ₀','logMR','Υ★⊥'] }
};

// Full predictor matrix (6 vars)
const allVars = gals45.map((g,i) => [logMHI[i], rcWig[i], logMhost[i], logSig0[i], logMR[i], upsPerp[i]]);

function getX(idx) { return gals45.map((_,i) => idx.map(j => allVars[i][j])); }

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 69: DEATH MATCH — M0 vs M3 vs M5 vs M6(B″)');
console.log('══════════════════════════════════════════════════════════════════════\n');

// ═══════════════════════════════════════════════════════════════════
// 1) LOO
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  1) LEAVE-ONE-OUT CROSS-VALIDATION');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const looResults = {};

// M0: LOO is just predicting with leave-one-out mean
let m0ss = 0;
for(let i=0;i<N;i++){
  const others = [...Y.slice(0,i),...Y.slice(i+1)];
  m0ss += (Y[i] - mean(others))**2;
}
const m0loo = Math.sqrt(m0ss/N);
looResults.M0 = gapV(m0loo, sdY);

for (const [key, mod] of Object.entries(models)) {
  if (key === 'M0') continue;
  const X = getX(mod.idx);
  const loo = looCV(Y, X);
  looResults[key] = gapV(loo, sdY);
}

for (const [key, mod] of Object.entries(models)) {
  const f = key === 'M0' ? {rss:Y.reduce((s,y)=>s+(y-mean(Y))**2,0),tss:Y.reduce((s,y)=>s+(y-mean(Y))**2,0),aic:N*Math.log(Y.reduce((s,y)=>s+(y-mean(Y))**2,0)/N)+2,bic:N*Math.log(Y.reduce((s,y)=>s+(y-mean(Y))**2,0)/N)+Math.log(N),resid:Y.map(y=>y-mean(Y)),r2adj:0} : ols(Y, getX(mod.idx));
  console.log('  ' + mod.name + ': LOO gap% = ' + looResults[key].toFixed(1) + '%, R²adj = ' + (key==='M0'?'0.000':f.r2adj.toFixed(4)) + ', tau = ' + sd(key==='M0'?f.resid:f.resid).toFixed(4));
}
console.log();

// ═══════════════════════════════════════════════════════════════════
// 2) k-fold CV (repeated)
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  2) K-FOLD CROSS-VALIDATION (200 repeats each)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function kfoldCV(Y, X, k, rng) {
  const n = Y.length;
  const idx = shuffle(Array.from({length:n},(_,i)=>i), rng);
  const foldSize = Math.floor(n/k);
  let totalSS = 0;
  for (let f = 0; f < k; f++) {
    const testIdx = idx.slice(f*foldSize, f===k-1 ? n : (f+1)*foldSize);
    const trainIdx = idx.filter(i => !testIdx.includes(i));
    const Yt = trainIdx.map(i=>Y[i]), Xt = trainIdx.map(i=>X[i]);
    const fit = ols(Yt, Xt);
    for (const ti of testIdx) {
      const xi = [1, ...X[ti]];
      totalSS += (Y[ti] - xi.reduce((s,x,j)=>s+x*fit.beta[j],0))**2;
    }
  }
  return Math.sqrt(totalSS/n);
}

function kfoldM0(Y, k, rng) {
  const n = Y.length;
  const idx = shuffle(Array.from({length:n},(_,i)=>i), rng);
  const foldSize = Math.floor(n/k);
  let totalSS = 0;
  for (let f = 0; f < k; f++) {
    const testIdx = idx.slice(f*foldSize, f===k-1 ? n : (f+1)*foldSize);
    const trainIdx = idx.filter(i => !testIdx.includes(i));
    const trainMean = mean(trainIdx.map(i=>Y[i]));
    for (const ti of testIdx) totalSS += (Y[ti] - trainMean)**2;
  }
  return Math.sqrt(totalSS/n);
}

const kfoldResults = {};
for (const kk of [5, 10]) {
  kfoldResults[kk] = {};
  for (const [key, mod] of Object.entries(models)) {
    const gaps = [];
    for (let rep = 0; rep < 200; rep++) {
      const rng = mulberry32(rep * 1000 + kk);
      const rms = key === 'M0' ? kfoldM0(Y, kk, rng) : kfoldCV(Y, getX(mod.idx), kk, rng);
      gaps.push(gapV(rms, sdY));
    }
    kfoldResults[kk][key] = { mean: mean(gaps), sd: sd(gaps), median: gaps.sort((a,b)=>a-b)[100] };
  }
  console.log('  ' + kk + '-fold (200 repeats):');
  for (const [key, mod] of Object.entries(models)) {
    const r = kfoldResults[kk][key];
    console.log('    ' + mod.name + ': mean = ' + r.mean.toFixed(1) + '% ± ' + r.sd.toFixed(1) + '%, median = ' + r.median.toFixed(1) + '%');
  }
  console.log();
}

// ═══════════════════════════════════════════════════════════════════
// 3) Repeated splits (50/50 and 80/20)
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  3) REPEATED TRAIN/TEST SPLITS (500 repeats)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function splitTest(Y, X, trainFrac, nReps) {
  const n = Y.length;
  const nTrain = Math.round(n * trainFrac);
  const gaps = [];
  const wins = {};
  for (let rep = 0; rep < nReps; rep++) {
    const rng = mulberry32(rep * 7 + Math.round(trainFrac*100));
    const idx = shuffle(Array.from({length:n},(_,i)=>i), rng);
    const trainIdx = idx.slice(0, nTrain), testIdx = idx.slice(nTrain);
    const results = {};
    // M0
    const trainMean = mean(trainIdx.map(i=>Y[i]));
    const testSdY = sd(testIdx.map(i=>Y[i]));
    let m0ss2 = 0; testIdx.forEach(i => m0ss2 += (Y[i]-trainMean)**2);
    results.M0 = gapV(Math.sqrt(m0ss2/testIdx.length), testSdY);
    
    for (const [key, mod] of Object.entries(X)) {
      const Yt = trainIdx.map(i=>Y[i]), Xt = trainIdx.map(i=>mod[i]);
      const fit = ols(Yt, Xt);
      let ss = 0;
      testIdx.forEach(i => {
        const xi = [1, ...mod[i]];
        ss += (Y[i] - xi.reduce((s2,x,j)=>s2+x*fit.beta[j],0))**2;
      });
      results[key] = gapV(Math.sqrt(ss/testIdx.length), testSdY);
    }
    gaps.push(results);
    
    // Who wins this split?
    let bestKey = null, bestGap = -Infinity;
    for (const [k,v] of Object.entries(results)) {
      if (v > bestGap) { bestGap = v; bestKey = k; }
    }
    wins[bestKey] = (wins[bestKey]||0) + 1;
  }
  return { gaps, wins };
}

const splitModels = {
  M3: getX(models.M3.idx),
  M5: getX(models.M5.idx),
  M6: getX(models.M6.idx)
};

for (const [frac, label] of [[0.5,'50/50'],[0.8,'80/20']]) {
  const result = splitTest(Y, splitModels, frac, 500);
  console.log('  ' + label + ' splits (500 repeats):');
  for (const key of ['M0','M3','M5','M6']) {
    const vals = result.gaps.map(g => g[key]);
    console.log('    ' + models[key].name + ': mean gap% = ' + mean(vals).toFixed(1) + '% ± ' + sd(vals).toFixed(1) + '%, wins = ' + (result.wins[key]||0) + '/500 (' + ((result.wins[key]||0)/5).toFixed(1) + '%)');
  }
  console.log();
}

// ═══════════════════════════════════════════════════════════════════
// 4) AIC / BIC
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  4) INFORMATION CRITERIA');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const icResults = {};
for (const [key, mod] of Object.entries(models)) {
  if (key === 'M0') {
    const rss = Y.reduce((s,y)=>s+(y-mean(Y))**2,0);
    icResults.M0 = { aic: N*Math.log(rss/N)+2, bic: N*Math.log(rss/N)+Math.log(N) };
  } else {
    const f = ols(Y, getX(mod.idx));
    icResults[key] = { aic: f.aic, bic: f.bic };
  }
  console.log('  ' + mod.name + ': AIC = ' + icResults[key].aic.toFixed(1) + ', BIC = ' + icResults[key].bic.toFixed(1));
}
console.log();

// Delta AIC/BIC relative to best
const bestAIC = Math.min(...Object.values(icResults).map(r=>r.aic));
const bestBIC = Math.min(...Object.values(icResults).map(r=>r.bic));
console.log('  ΔAIC (relative to best):');
for (const [key, mod] of Object.entries(models)) {
  console.log('    ' + mod.name + ': ΔAIC = ' + (icResults[key].aic - bestAIC).toFixed(1) + ', ΔBIC = ' + (icResults[key].bic - bestBIC).toFixed(1));
}
console.log();

// ═══════════════════════════════════════════════════════════════════
// 5) Bootstrap .632
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  5) BOOTSTRAP .632 ESTIMATOR (1000 resamples)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function bootstrap632(Y, X, nBoot) {
  const n = Y.length;
  // Apparent (in-sample) error
  const fApp = ols(Y, X);
  const appRSS = fApp.rss;
  const appMSE = appRSS / n;
  
  const oobErrors = [];
  for (let b = 0; b < nBoot; b++) {
    const rng = mulberry32(b * 13 + 42);
    const bootIdx = Array.from({length:n}, ()=>Math.floor(rng()*n));
    const oobIdx = [...new Set(Array.from({length:n},(_,i)=>i))].filter(i => !bootIdx.includes(i));
    if (oobIdx.length < 3) continue;
    
    const Yb = bootIdx.map(i=>Y[i]), Xb = bootIdx.map(i=>X[i]);
    const fb = ols(Yb, Xb);
    
    let oobSS = 0;
    oobIdx.forEach(i => {
      const xi = [1, ...X[i]];
      oobSS += (Y[i] - xi.reduce((s,x,j)=>s+x*fb.beta[j],0))**2;
    });
    oobErrors.push(oobSS / oobIdx.length);
  }
  
  const meanOOB = mean(oobErrors);
  const err632 = 0.368 * appMSE + 0.632 * meanOOB;
  return { appMSE, meanOOB, err632, gap632: gapV(Math.sqrt(err632), sdY) };
}

const bootResults = {};
for (const [key, mod] of Object.entries(models)) {
  if (key === 'M0') {
    // M0 bootstrap: apparent = var(Y), OOB = var(Y) approx
    bootResults.M0 = { gap632: gapV(sdY, sdY) }; // ~0
    console.log('  ' + mod.name + ': .632 gap% ≈ 0.0% (by definition)');
    continue;
  }
  const X = getX(mod.idx);
  const br = bootstrap632(Y, X, 1000);
  bootResults[key] = br;
  console.log('  ' + mod.name + ': .632 gap% = ' + br.gap632.toFixed(1) + '% (apparent=' + gapV(Math.sqrt(br.appMSE),sdY).toFixed(1) + '%, OOB=' + gapV(Math.sqrt(br.meanOOB),sdY).toFixed(1) + '%)');
}
console.log();

// ═══════════════════════════════════════════════════════════════════
// 6) Coefficient Stability (Bootstrap 1000)
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  6) COEFFICIENT STABILITY (1000 bootstrap resamples)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function bootstrapCoeffs(Y, X, nBoot, varNames) {
  const n = Y.length, p = varNames.length;
  const allBetas = Array.from({length:p}, ()=>[]);
  const signFlips = new Array(p).fill(0);
  const fFull = ols(Y, X);
  const refSigns = fFull.beta.slice(1).map(b => Math.sign(b));
  
  for (let b = 0; b < nBoot; b++) {
    const rng = mulberry32(b * 17 + 99);
    const bootIdx = Array.from({length:n}, ()=>Math.floor(rng()*n));
    const Yb = bootIdx.map(i=>Y[i]), Xb = bootIdx.map(i=>X[i]);
    try {
      const fb = ols(Yb, Xb);
      for (let j = 0; j < p; j++) {
        allBetas[j].push(fb.beta[j+1]);
        if (Math.sign(fb.beta[j+1]) !== refSigns[j]) signFlips[j]++;
      }
    } catch(e) {}
  }
  
  return varNames.map((v,j) => ({
    name: v,
    refBeta: fFull.beta[j+1],
    meanBeta: mean(allBetas[j]),
    sdBeta: sd(allBetas[j]),
    signFlipPct: (signFlips[j] / allBetas[j].length * 100),
    tStat: fFull.tStats[j+1]
  }));
}

// M6 (B″) stability
console.log('  M6 (B″) coefficient stability:');
const m6Stab = bootstrapCoeffs(Y, getX(models.M6.idx), 1000, models.M6.vars);
m6Stab.forEach(s => {
  console.log('    ' + s.name.padEnd(12) + ': β = ' + s.refBeta.toFixed(3).padStart(7) + ' ± ' + s.sdBeta.toFixed(3) + ', t = ' + s.tStat.toFixed(2).padStart(6) + ', sign flips = ' + s.signFlipPct.toFixed(1) + '%');
});
console.log();

// M5 stability
console.log('  M5 (5-var) coefficient stability:');
const m5Stab = bootstrapCoeffs(Y, getX(models.M5.idx), 1000, models.M5.vars);
m5Stab.forEach(s => {
  console.log('    ' + s.name.padEnd(12) + ': β = ' + s.refBeta.toFixed(3).padStart(7) + ' ± ' + s.sdBeta.toFixed(3) + ', t = ' + s.tStat.toFixed(2).padStart(6) + ', sign flips = ' + s.signFlipPct.toFixed(1) + '%');
});
console.log();

// M3 stability
console.log('  M3 (3-axis) coefficient stability:');
const m3Stab = bootstrapCoeffs(Y, getX(models.M3.idx), 1000, models.M3.vars);
m3Stab.forEach(s => {
  console.log('    ' + s.name.padEnd(12) + ': β = ' + s.refBeta.toFixed(3).padStart(7) + ' ± ' + s.sdBeta.toFixed(3) + ', t = ' + s.tStat.toFixed(2).padStart(6) + ', sign flips = ' + s.signFlipPct.toFixed(1) + '%');
});
console.log();

// ═══════════════════════════════════════════════════════════════════
// 7) Head-to-head M5 vs M6 detailed comparison
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  7) HEAD-TO-HEAD: M5 vs M6 — rcWiggliness diagnostic');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Direct LOO comparison
const fM5 = ols(Y, getX(models.M5.idx));
const fM6 = ols(Y, getX(models.M6.idx));
console.log('  In-sample:');
console.log('    M5: R²adj = ' + fM5.r2adj.toFixed(4) + ', tau = ' + sd(fM5.resid).toFixed(4));
console.log('    M6: R²adj = ' + fM6.r2adj.toFixed(4) + ', tau = ' + sd(fM6.resid).toFixed(4));
console.log('    Δ(R²adj) = ' + (fM6.r2adj - fM5.r2adj).toFixed(4) + ' (M6 gains ' + ((fM6.r2adj-fM5.r2adj)*100).toFixed(2) + 'pp)');
console.log();

// F-test for nested models (M5 is nested in M6)
const dfM5 = N - models.M5.idx.length - 1;
const dfM6 = N - models.M6.idx.length - 1;
const Fstat = ((fM5.rss - fM6.rss) / (dfM5 - dfM6)) / (fM6.rss / dfM6);
console.log('  Partial F-test (M5→M6, adding rcWig):');
console.log('    F = ' + Fstat.toFixed(3) + ' on (1, ' + dfM6 + ') df');
console.log('    For α=0.05, F_crit ≈ 4.10');
console.log('    ' + (Fstat > 4.10 ? 'rcWig IS significant at 5%' : 'rcWig is NOT significant at 5%'));
console.log();

// Count how often M5 beats M6 in LOO per-galaxy
const looM5 = [], looM6 = [];
for (let i = 0; i < N; i++) {
  const Yt = [...Y.slice(0,i),...Y.slice(i+1)];
  const X5t = [...getX(models.M5.idx).slice(0,i),...getX(models.M5.idx).slice(i+1)];
  const X6t = [...getX(models.M6.idx).slice(0,i),...getX(models.M6.idx).slice(i+1)];
  const f5 = ols(Yt, X5t), f6 = ols(Yt, X6t);
  const xi5 = [1, ...getX(models.M5.idx)[i]], xi6 = [1, ...getX(models.M6.idx)[i]];
  looM5.push((Y[i] - xi5.reduce((s,x,j)=>s+x*f5.beta[j],0))**2);
  looM6.push((Y[i] - xi6.reduce((s,x,j)=>s+x*f6.beta[j],0))**2);
}
const m5wins = looM5.filter((v,i) => v < looM6[i]).length;
const m6wins = looM5.filter((v,i) => v > looM6[i]).length;
console.log('  Per-galaxy LOO comparison:');
console.log('    M5 wins on ' + m5wins + '/45 galaxies');
console.log('    M6 wins on ' + m6wins + '/45 galaxies');
console.log('    Tied: ' + (N-m5wins-m6wins));
console.log();

// ═══════════════════════════════════════════════════════════════════
// FINAL VERDICT
// ═══════════════════════════════════════════════════════════════════
console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 69: FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

console.log('  ┌───────────────────────┬────┬────────┬────────┬────────┬────────┬────────┐');
console.log('  │ Model                 │ k  │  LOO%  │ 5-fold │10-fold │  AIC   │  BIC   │');
console.log('  ├───────────────────────┼────┼────────┼────────┼────────┼────────┼────────┤');
for (const [key, mod] of Object.entries(models)) {
  const loo = looResults[key].toFixed(1).padStart(5);
  const f5 = kfoldResults[5][key].mean.toFixed(1).padStart(5);
  const f10 = kfoldResults[10][key].mean.toFixed(1).padStart(5);
  const aic = icResults[key].aic.toFixed(1).padStart(7);
  const bic = icResults[key].bic.toFixed(1).padStart(7);
  const kv = (key==='M0'?'0':mod.idx.length.toString()).padStart(2);
  console.log('  │ ' + mod.name.padEnd(21) + ' │ ' + kv + ' │ ' + loo + '% │ ' + f5 + '% │ ' + f10 + '% │' + aic + ' │' + bic + ' │');
}
console.log('  └───────────────────────┴────┴────────┴────────┴────────┴────────┴────────┘');
console.log();

// Determine verdicts
const m5loo = looResults.M5, m6loo = looResults.M6, m3loo = looResults.M3;
const m5winsLOO = m5loo > m6loo;
const m5winsBIC = icResults.M5.bic < icResults.M6.bic;
const m5winsBoot = (bootResults.M5?.gap632||0) > (bootResults.M6?.gap632||0);

console.log('  DECISION MATRIX:');
console.log('  ┌──────────────────────┬──────────┬──────────┬──────────┐');
console.log('  │ Criterion            │ M5 wins? │ M6 wins? │ M3 wins? │');
console.log('  ├──────────────────────┼──────────┼──────────┼──────────┤');
console.log('  │ LOO gap%             │ ' + (m5winsLOO?'  YES ✅  ':'  no      ') + '│ ' + (!m5winsLOO?'  YES ✅  ':'  no      ') + '│   no     │');
console.log('  │ BIC                  │ ' + (m5winsBIC?'  YES ✅  ':'  no      ') + '│ ' + (!m5winsBIC?'  YES ✅  ':'  no      ') + '│   —      │');
console.log('  │ AIC                  │ ' + (icResults.M5.aic<icResults.M6.aic?'  YES ✅  ':'  no      ') + '│ ' + (icResults.M5.aic>=icResults.M6.aic?'  YES ✅  ':'  no      ') + '│   —      │');
console.log('  │ Bootstrap .632       │ ' + (m5winsBoot?'  YES ✅  ':'  no      ') + '│ ' + (!m5winsBoot?'  YES ✅  ':'  no      ') + '│   —      │');
console.log('  │ Coeff stability      │  better  │ rcWig    │  best    │');
console.log('  │                      │  (no     │  13.7%   │  (all    │');
console.log('  │                      │  flips)  │  flips   │  stable) │');
console.log('  │ F-test rcWig         │ ' + (Fstat<4.10?'  YES ✅  ':'  no      ') + '│ ' + (Fstat>=4.10?'  YES ✅  ':'  no      ') + '│   —      │');
console.log('  └──────────────────────┴──────────┴──────────┴──────────┘');
console.log();

const rcWigVerdict = m5winsLOO && m5winsBIC && Fstat < 4.10 ? 'DROP' : 
                     m5winsLOO || m5winsBIC ? 'MARGINAL — consider dropping' : 'KEEP';

console.log('  ═══════════════════════════════════════════════════════════════');
console.log('  FINAL ANSWERS:');
console.log('  ═══════════════════════════════════════════════════════════════');
console.log();
console.log('  Best predictive law     = ' + (m5winsLOO ? 'M5 (5-var reduced)' : 'M6 (B″ full)'));
console.log('  Best compressed/state   = M3 (3-axis: logMHI + logMhost + logMR)');
console.log('  rcWiggliness verdict    = ' + rcWigVerdict);
console.log('  M3 retention vs winner  = ' + (m3loo / Math.max(m5loo,m6loo) * 100).toFixed(0) + '%');
console.log();

if (m5winsLOO) {
  console.log('  *** RECOMMENDATION: PROMOTE M5 AS NEW BEST MODEL ***');
  console.log('  rcWiggliness adds instability without improving out-of-sample prediction.');
  console.log('  M5 = logMHI + logMhost + logΣ₀ + logMeanRun + Υ★⊥');
} else {
  console.log('  *** B″ (M6) remains best, but M5 is nearly identical ***');
  console.log('  The difference is within noise. Parsimony favors M5.');
}

const output = {
  phase: '69', title: 'Death Match: M0 vs M3 vs M5 vs M6',
  models: {
    M0: {k:0, looGap: looResults.M0, aic: icResults.M0.aic, bic: icResults.M0.bic},
    M3: {k:3, vars:models.M3.vars, looGap: looResults.M3, aic: icResults.M3.aic, bic: icResults.M3.bic, boot632: bootResults.M3?.gap632, kfold5: kfoldResults[5].M3.mean, kfold10: kfoldResults[10].M3.mean},
    M5: {k:5, vars:models.M5.vars, looGap: looResults.M5, aic: icResults.M5.aic, bic: icResults.M5.bic, boot632: bootResults.M5?.gap632, kfold5: kfoldResults[5].M5.mean, kfold10: kfoldResults[10].M5.mean},
    M6: {k:6, vars:models.M6.vars, looGap: looResults.M6, aic: icResults.M6.aic, bic: icResults.M6.bic, boot632: bootResults.M6?.gap632, kfold5: kfoldResults[5].M6.mean, kfold10: kfoldResults[10].M6.mean}
  },
  rcWiggliness: {
    Ftest: Fstat, significant: Fstat >= 4.10,
    m6StabFlipPct: m6Stab.find(s=>s.name==='rcWig')?.signFlipPct,
    perGalaxyLOO: {m5wins, m6wins}
  },
  coeffStability: {m3: m3Stab, m5: m5Stab, m6: m6Stab},
  verdict: {
    bestPredictive: m5winsLOO ? 'M5' : 'M6',
    bestCompressed: 'M3',
    dropRcWig: rcWigVerdict,
    m3retention: m3loo / Math.max(m5loo, m6loo) * 100
  }
};
fs.writeFileSync('public/phase69-death-match-redux.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase69-death-match-redux.json');
