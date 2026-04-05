const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const tidalData = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const tdMap = {};
tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {};
sparc.forEach(s => { sparcMap[s.name] = s; });
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
const N = 45;

function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);return{beta,resid,rss,tss,se:Math.sqrt(rss/(n-p)),n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}

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
function gapV(rms){return 100*(1-rms**2/sdY**2);}

const Xconf = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT[i]]);
const fConf = ols(logUps, Xconf);
const upsPerp = fConf.resid;

// Raw 6 variables
const logMHI = gals45.map(g => g.logMHI);
const rcWig = gals45.map(g => g.rcWiggliness);
const logSig0 = gals45.map(g => g.logSigma0);
const logMR = gals45.map(g => g.logMeanRun);

const X_Bdp = gals45.map((g,i) => [logMHI[i], rcWig[i], logMhost[i], logSig0[i], logMR[i], upsPerp[i]]);
const fBdp = ols(Y, X_Bdp);
const looBase = looCV(Y, X_Bdp);

console.log('======================================================================');
console.log('  PHASE 68: STATE-AXIS COMPRESSION');
console.log('  Compress B″ (6 vars) → 3 physical state axes');
console.log('======================================================================\n');

console.log('  B″ baseline: LOO gap% = ' + gapV(looBase).toFixed(1) + '%, tau = ' + sd(fBdp.resid).toFixed(4) + ' dex\n');

// ═══════════════════════════════════════════════════════════════════
// Define 3 physical state axes
// ═══════════════════════════════════════════════════════════════════

// AXIS 1: BARYONIC STATE — combines gas mass + surface density + stellar M/L
// Physical meaning: total baryonic configuration
// Components: logMHI, logΣ₀, Υ★⊥

// AXIS 2: DYNAMICAL STATE — combines RC irregularity + kinematic coherence
// Physical meaning: how far the galaxy is from equilibrium / how structured its dynamics are
// Components: rcWiggliness, logMeanRun

// AXIS 3: ENVIRONMENTAL STATE — host halo mass
// Physical meaning: depth of gravitational environment
// Components: logMhost (single variable — already a clean axis)

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  68a: SUPERVISED COMPRESSION — optimal linear combinations');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Method: for each group, find optimal linear combination by regressing 
// logA0 on group members only, then use fitted values as the composite axis

// Standardize all variables first
function standardize(arr) {
  const m = mean(arr), s = sd(arr);
  return { vals: arr.map(v => (v-m)/s), m, s };
}

const sLogMHI = standardize(logMHI);
const sRcWig = standardize(rcWig);
const sLogMhost = standardize(logMhost);
const sLogSig0 = standardize(logSig0);
const sLogMR = standardize(logMR);
const sUpsPerp = standardize(upsPerp);

// Baryonic axis: supervised combination of (logMHI, logΣ₀, Υ★⊥)
const X_bary = gals45.map((_,i) => [sLogMHI.vals[i], sLogSig0.vals[i], sUpsPerp.vals[i]]);
const f_bary = ols(Y, X_bary);
const A_bary = gals45.map((_,i) => [1, ...X_bary[i]].reduce((s,x,j) => s+x*f_bary.beta[j], 0));

console.log('  Baryonic axis (logMHI + logΣ₀ + Υ★⊥):');
console.log('    Weights (standardized): logMHI=' + f_bary.beta[1].toFixed(3) + ', logΣ₀=' + f_bary.beta[2].toFixed(3) + ', Υ★⊥=' + f_bary.beta[3].toFixed(3));
console.log('    Alone R² = ' + (1-f_bary.rss/f_bary.tss).toFixed(4));
const looBary = looCV(Y, X_bary);
console.log('    Alone LOO gap% = ' + gapV(looBary).toFixed(1) + '%');
console.log();

// Dynamical axis: supervised combination of (rcWig, logMeanRun)
const X_dyn = gals45.map((_,i) => [sRcWig.vals[i], sLogMR.vals[i]]);
const f_dyn = ols(Y, X_dyn);
const A_dyn = gals45.map((_,i) => [1, ...X_dyn[i]].reduce((s,x,j) => s+x*f_dyn.beta[j], 0));

console.log('  Dynamical axis (rcWig + logMeanRun):');
console.log('    Weights (standardized): rcWig=' + f_dyn.beta[1].toFixed(3) + ', logMeanRun=' + f_dyn.beta[2].toFixed(3));
console.log('    Alone R² = ' + (1-f_dyn.rss/f_dyn.tss).toFixed(4));
const looDyn = looCV(Y, X_dyn);
console.log('    Alone LOO gap% = ' + gapV(looDyn).toFixed(1) + '%');
console.log();

// Environmental axis: logMhost (single variable, no compression needed)
const X_env = gals45.map((_,i) => [sLogMhost.vals[i]]);
const f_env = ols(Y, X_env);

console.log('  Environmental axis (logMhost):');
console.log('    Weight (standardized): logMhost=' + f_env.beta[1].toFixed(3));
console.log('    Alone R² = ' + (1-f_env.rss/f_env.tss).toFixed(4));
const looEnv = looCV(Y, X_env);
console.log('    Alone LOO gap% = ' + gapV(looEnv).toFixed(1) + '%');
console.log();

// ═══════════════════════════════════════════════════════════════════
// 68b: 3-axis compressed model
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  68b: 3-AXIS COMPRESSED MODEL');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Method 1: Use raw group members as separate predictors but organized
// This is the "natural 3-group" model
const X_3axis = gals45.map((_,i) => [
  // Baryonic: logMHI, logΣ₀, Υ★⊥ → combined as one predictor per group
  // But we need to be careful: compression must happen INSIDE the LOO loop
  // to avoid data leakage.
  // 
  // Approach A: Pre-fixed weights (from full sample)
  // Approach B: Re-estimate weights in each LOO fold
  // 
  // For comparison, let's first test with the 3 group-composite scores
  // computed from full-sample weights (Approach A — slight optimism)
  // then verify with Approach B
  
  // For now: use the SIMPLEST member from each group
  logMHI[i],        // strongest baryonic
  logMR[i],         // strongest dynamical  
  logMhost[i]       // environmental
]);

const f3simple = ols(Y, X_3axis);
const loo3simple = looCV(Y, X_3axis);

console.log('  Model C3a: simplest-member-per-group (logMHI + logMR + logMhost):');
console.log('    R² = ' + (1-f3simple.rss/f3simple.tss).toFixed(4));
console.log('    LOO gap% = ' + gapV(loo3simple).toFixed(1) + '%');
console.log('    tau = ' + sd(f3simple.resid).toFixed(4) + ' dex');
console.log('    vs B″: Δgap = ' + (gapV(loo3simple)-gapV(looBase)).toFixed(1) + 'pp');
console.log();

// Method 2: best single from each group (by t-stat in B″)
// logMHI (t=-4.78), logMeanRun (t=+3.80), logMhost (t=-3.08)
// → same as above. Good.

// Method 3: all 3 baryonic together as group, 2 dyn together, 1 env
// = full B″ but with group labels → this IS B″

// Method 4: PCA within each group → 1 PC per group
console.log('  PCA within each group:');

// Baryonic PCA (3 vars → 1 PC)
function pcaGroup(vars) {
  const n = vars[0].length, p = vars.length;
  // Standardize
  const Z = vars.map(v => { const m = mean(v), s = sd(v); return v.map(x => (x-m)/s); });
  // Correlation matrix
  const C = Array.from({length:p}, () => new Array(p).fill(0));
  for(let i=0;i<p;i++) for(let j=0;j<p;j++) {
    let s=0; for(let k=0;k<n;k++) s+=Z[i][k]*Z[j][k];
    C[i][j]=s/(n-1);
  }
  // Power iteration
  let v = Array.from({length:p}, ()=>1/Math.sqrt(p));
  for(let iter=0;iter<200;iter++){
    const Cv=new Array(p).fill(0);
    for(let i=0;i<p;i++) for(let j=0;j<p;j++) Cv[i]+=C[i][j]*v[j];
    const norm=Math.sqrt(Cv.reduce((s,x)=>s+x*x,0));
    v=Cv.map(x=>x/norm);
  }
  // Project
  const scores = Array.from({length:n}, (_,k) => v.reduce((s,w,j) => s+w*Z[j][k], 0));
  const eigenval = v.reduce((s,w,j) => {
    let cv=0; for(let l=0;l<p;l++) cv+=C[j][l]*v[l];
    return s+w*cv;
  },0);
  return {scores, loadings:v, varExplained: eigenval/p};
}

const baryPCA = pcaGroup([logMHI, logSig0, upsPerp]);
console.log('    Baryonic PC1: var explained = ' + (baryPCA.varExplained*100).toFixed(1) + '%');
console.log('      Loadings: logMHI=' + baryPCA.loadings[0].toFixed(3) + ', logΣ₀=' + baryPCA.loadings[1].toFixed(3) + ', Υ★⊥=' + baryPCA.loadings[2].toFixed(3));

const dynPCA = pcaGroup([rcWig, logMR]);
console.log('    Dynamical PC1: var explained = ' + (dynPCA.varExplained*100).toFixed(1) + '%');
console.log('      Loadings: rcWig=' + dynPCA.loadings[0].toFixed(3) + ', logMR=' + dynPCA.loadings[1].toFixed(3));

// Env is single var, no PCA needed
console.log('    Environmental: logMhost (single var)');
console.log();

// 3-PC model
const X_3pc = gals45.map((_,i) => [baryPCA.scores[i], dynPCA.scores[i], logMhost[i]]);
const f3pc = ols(Y, X_3pc);
const loo3pc = looCV(Y, X_3pc);

console.log('  Model C3b: PCA-compressed (baryPC1 + dynPC1 + logMhost):');
console.log('    R² = ' + (1-f3pc.rss/f3pc.tss).toFixed(4));
console.log('    LOO gap% = ' + gapV(loo3pc).toFixed(1) + '%');
console.log('    tau = ' + sd(f3pc.resid).toFixed(4) + ' dex');
console.log('    vs B″: Δgap = ' + (gapV(loo3pc)-gapV(looBase)).toFixed(1) + 'pp');
console.log();

// Method 5: Supervised composites (fitted values from within-group regressions)
// Use within-group regression weights to build composites, then regress logA0 on composites
// BUT: to avoid leakage, just use the raw variables directly

// Method 6: All 2-axis combinations
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  68c: ALL 2-AXIS AND 3-AXIS COMBINATIONS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const axisGroups = {
  bary: {name:'Baryonic', vars:['logMHI','logΣ₀','Υ★⊥'], idx:[0,3,5]},
  dyn:  {name:'Dynamical', vars:['rcWig','logMR'], idx:[1,4]},
  env:  {name:'Environmental', vars:['logMhost'], idx:[2]}
};

// Test all pair combos
const pairs = [['bary','dyn'], ['bary','env'], ['dyn','env']];
const tripleVars = ['bary','dyn','env'];

pairs.forEach(([a,b]) => {
  const idxs = [...axisGroups[a].idx, ...axisGroups[b].idx];
  const X_pair = gals45.map((_,i) => idxs.map(j => X_Bdp[i][j]));
  const fPair = ols(Y, X_pair);
  const looPair = looCV(Y, X_pair);
  console.log('  ' + axisGroups[a].name + ' + ' + axisGroups[b].name + ' (' + idxs.length + ' vars):');
  console.log('    R² = ' + (1-fPair.rss/fPair.tss).toFixed(4) + ', LOO gap% = ' + gapV(looPair).toFixed(1) + '%, tau = ' + sd(fPair.resid).toFixed(4));
});
console.log();

// Full triple = B″
console.log('  All 3 groups (= B″, 6 vars):');
console.log('    R² = ' + (1-fBdp.rss/fBdp.tss).toFixed(4) + ', LOO gap% = ' + gapV(looBase).toFixed(1) + '%, tau = ' + sd(fBdp.resid).toFixed(4));
console.log();

// ═══════════════════════════════════════════════════════════════════
// 68d: Best 3-variable and 4-variable models
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  68d: BEST k-VARIABLE MODELS (exhaustive search)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const varNames = ['logMHI','rcWig','logMhost','logΣ₀','logMR','Υ★⊥'];

// All choose-k combinations
function combos(arr, k) {
  if(k===0)return[[]];if(k>arr.length)return[];
  const r=[];for(let i=0;i<=arr.length-k;i++){
    const sub=combos(arr.slice(i+1),k-1);
    sub.forEach(s=>r.push([arr[i],...s]));
  }return r;
}

for (let k = 1; k <= 5; k++) {
  const allCombos = combos([0,1,2,3,4,5], k);
  let best = null, bestGap = -Infinity;
  allCombos.forEach(combo => {
    const X_k = gals45.map((_,i) => combo.map(j => X_Bdp[i][j]));
    const loo_k = looCV(Y, X_k);
    const gap_k = gapV(loo_k);
    if (gap_k > bestGap) { bestGap = gap_k; best = combo; }
  });
  const X_best = gals45.map((_,i) => best.map(j => X_Bdp[i][j]));
  const f_best = ols(Y, X_best);
  console.log('  Best ' + k + '-var: [' + best.map(j=>varNames[j]).join(', ') + ']');
  console.log('    LOO gap% = ' + bestGap.toFixed(1) + '%, tau = ' + sd(f_best.resid).toFixed(4) + ', R² = ' + (1-f_best.rss/f_best.tss).toFixed(4));
}
console.log();
console.log('  Full 6-var (B″): LOO gap% = ' + gapV(looBase).toFixed(1) + '%');
console.log();

// ═══════════════════════════════════════════════════════════════════
// FINAL VERDICT
// ═══════════════════════════════════════════════════════════════════
console.log('======================================================================');
console.log('  PHASE 68: COMPRESSION VERDICT');
console.log('======================================================================\n');

const best3gap = gapV(loo3simple);
const loss3 = gapV(looBase) - best3gap;

console.log('  ┌──────────────────────┬──────────┬──────────┬──────────┐');
console.log('  │ Model                │ k vars   │ LOO gap% │ Loss pp  │');
console.log('  ├──────────────────────┼──────────┼──────────┼──────────┤');
console.log('  │ B″ (full)            │    6     │ ' + gapV(looBase).toFixed(1).padStart(6) + '   │   0.0    │');

// Re-compute best 3-var
const allC3 = combos([0,1,2,3,4,5],3);
let best3combo=null, best3gapV=-Infinity;
allC3.forEach(c => {
  const X3 = gals45.map((_,i) => c.map(j => X_Bdp[i][j]));
  const l3 = looCV(Y, X3);
  const g3 = gapV(l3);
  if(g3>best3gapV){best3gapV=g3;best3combo=c;}
});
console.log('  │ Best 3-var           │    3     │ ' + best3gapV.toFixed(1).padStart(6) + '   │ ' + (gapV(looBase)-best3gapV).toFixed(1).padStart(5) + '    │');
console.log('  │ C3a (MHI+MR+Mhost)   │    3     │ ' + gapV(loo3simple).toFixed(1).padStart(6) + '   │ ' + (gapV(looBase)-gapV(loo3simple)).toFixed(1).padStart(5) + '    │');
console.log('  │ C3b (PCA composites) │    3     │ ' + gapV(loo3pc).toFixed(1).padStart(6) + '   │ ' + (gapV(looBase)-gapV(loo3pc)).toFixed(1).padStart(5) + '    │');
console.log('  │ Null (M0)            │    0     │  ' + '-2.3'.padStart(5) + '   │ ' + (gapV(looBase)+2.3).toFixed(1).padStart(5) + '    │');
console.log('  └──────────────────────┴──────────┴──────────┴──────────┘');
console.log();

const compressionWorks = best3gapV > gapV(looBase) * 0.80;
console.log('  3-var retains ≥80% of B″ performance? ' + (compressionWorks ? 'YES ✅' : 'NO ❌'));
console.log('  Best 3-var: [' + best3combo.map(j=>varNames[j]).join(', ') + ']');
console.log('  Retention: ' + (best3gapV/gapV(looBase)*100).toFixed(0) + '%');
console.log();

if (compressionWorks) {
  console.log('  ═══════════════════════════════════════════════════════════');
  console.log('  COMPRESSION SUCCESSFUL');
  console.log('  B″ can be approximated with 3 state variables');
  console.log('  losing only ' + (gapV(looBase)-best3gapV).toFixed(1) + ' percentage points');
  console.log('  ═══════════════════════════════════════════════════════════');
} else {
  console.log('  ═══════════════════════════════════════════════════════════');
  console.log('  COMPRESSION PARTIAL — all 6 axes needed for full performance');
  console.log('  but 3 axes capture ' + (best3gapV/gapV(looBase)*100).toFixed(0) + '% of the signal');
  console.log('  ═══════════════════════════════════════════════════════════');
}

const output = {
  phase: '68', title: 'State-Axis Compression',
  baseline: {model:'B″', k:6, looGap: gapV(looBase), tau: sd(fBdp.resid)},
  physicalAxes: {
    baryonic: {vars:['logMHI','logΣ₀','Υ★⊥'], pcaVarExplained: baryPCA.varExplained, pcaLoadings: baryPCA.loadings},
    dynamical: {vars:['rcWig','logMeanRun'], pcaVarExplained: dynPCA.varExplained, pcaLoadings: dynPCA.loadings},
    environmental: {vars:['logMhost']}
  },
  compressed: {
    c3a: {vars:['logMHI','logMR','logMhost'], looGap: gapV(loo3simple), tau: sd(f3simple.resid)},
    c3b_pca: {vars:['baryPC1','dynPC1','logMhost'], looGap: gapV(loo3pc), tau: sd(f3pc.resid)},
    best3: {vars: best3combo.map(j=>varNames[j]), looGap: best3gapV}
  },
  retention: best3gapV/gapV(looBase)*100,
  compressionWorks
};
fs.writeFileSync('public/phase68-state-compression.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase68-state-compression.json');
