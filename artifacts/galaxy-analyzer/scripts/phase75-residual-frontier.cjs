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
function median(a){const s=[...a].sort((x,y)=>x-y);return s.length%2?s[(s.length-1)/2]:(s[s.length/2-1]+s[s.length/2])/2;}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));return{beta,resid};}
function pearsonCorr(x,y){if(x.length<3)return 0;const mx=mean(x),my=mean(y);let num=0,dx2=0,dy2=0;for(let i=0;i<x.length;i++){num+=(x[i]-mx)*(y[i]-my);dx2+=(x[i]-mx)**2;dy2+=(y[i]-my)**2;}if(dx2<1e-30||dy2<1e-30)return 0;return num/Math.sqrt(dx2*dy2);}
function spearmanRank(a){const s=[...a].map((v,i)=>({v,i})).sort((x,y)=>x.v-y.v);const r=new Array(a.length);s.forEach((x,i)=>r[x.i]=i+1);return r;}
function spearmanCorr(x,y){if(x.length<3)return 0;const rx=spearmanRank(x),ry=spearmanRank(y);const n=x.length;let d2=0;for(let i=0;i<n;i++)d2+=(rx[i]-ry[i])**2;return 1-6*d2/(n*(n*n-1));}
function tTest2(a,b){const na=a.length,nb=b.length;const ma=mean(a),mb=mean(b);const va=a.reduce((s,v)=>s+(v-ma)**2,0)/(na-1);const vb=b.reduce((s,v)=>s+(v-mb)**2,0)/(nb-1);const se=Math.sqrt(va/na+vb/nb);return{diff:ma-mb,t:(ma-mb)/(se||1e-15),sdA:Math.sqrt(va),sdB:Math.sqrt(vb)};}

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
const morphT = gals45.map(g => sparcMap[g.name]?.T || 5);
const logUps = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const Xconf = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT[i]]);
const fConf = ols(logUps, Xconf);
const upsPerp = fConf.resid;

const logMHI = gals45.map(g => g.logMHI);
const logSig0 = gals45.map(g => g.logSigma0);
const logMR = gals45.map(g => g.logMeanRun);
const envCode = gals45.map(g => g.envCode || 0);
const dist = gals45.map(g => g.dist_sparc || 0);

const X_M5 = gals45.map((_,i) => [logMHI[i], logMhost[i], logSig0[i], logMR[i], upsPerp[i]]);
const fM5 = ols(Y, X_M5);
const resid = fM5.resid;
const predicted = Y.map((y,i) => y - resid[i]);
const absResid = resid.map(r => Math.abs(r));

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 75: RESIDUAL PHYSICS FRONTIER');
console.log('══════════════════════════════════════════════════════════════════════\n');

const residRMS = Math.sqrt(resid.reduce((s,r)=>s+r*r,0)/N);
const residSD = sd(resid);
console.log('  M5 residual: RMS=' + residRMS.toFixed(3) + ' dex, SD=' + residSD.toFixed(3) + ' dex');
console.log('  Mean residual: ' + mean(resid).toFixed(4) + ' (should be ~0)');
console.log();

// ═════════ 75a: RESIDUAL STRATIFICATION ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  75a: RESIDUAL STRATIFICATION');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function stratify(name, values, splitFn) {
  const lo = [], hi = [];
  values.forEach((v, i) => { if (splitFn(v, i)) hi.push(i); else lo.push(i); });
  const rLo = lo.map(i => resid[i]), rHi = hi.map(i => resid[i]);
  const arLo = lo.map(i => absResid[i]), arHi = hi.map(i => absResid[i]);
  const tt = tTest2(rHi, rLo);
  const flag = Math.abs(tt.t) > 2 ? ' *** SIGNIFICANT' : Math.abs(tt.t) > 1.5 ? ' ** marginal' : '';
  console.log('  ' + name + ':');
  console.log('    High group (n=' + hi.length + '): mean resid = ' + (mean(rHi)>0?'+':'') + mean(rHi).toFixed(3) + ', SD = ' + sd(rHi).toFixed(3) + ', mean |resid| = ' + mean(arHi).toFixed(3));
  console.log('    Low group  (n=' + lo.length + '): mean resid = ' + (mean(rLo)>0?'+':'') + mean(rLo).toFixed(3) + ', SD = ' + sd(rLo).toFixed(3) + ', mean |resid| = ' + mean(arLo).toFixed(3));
  console.log('    Δmean = ' + tt.diff.toFixed(3) + ', t = ' + tt.t.toFixed(2) + flag);
  console.log('    SD ratio (high/low) = ' + (sd(rHi)/sd(rLo)).toFixed(2));
  console.log();
  return {name, nHi:hi.length, nLo:lo.length, meanHi:mean(rHi), meanLo:mean(rLo), sdHi:sd(rHi), sdLo:sd(rLo), t:tt.t, sig:Math.abs(tt.t)>2};
}

const medMhost = median(logMhost);
const medMHI = median(logMHI);
const medMR = median(logMR);
const medMorphT = median(morphT);

const strats = [];
strats.push(stratify('logMhost (median split)', logMhost, v => v >= medMhost));
strats.push(stratify('logMHI (median split)', logMHI, v => v >= medMHI));
strats.push(stratify('logMeanRun (median split)', logMR, v => v >= medMR));
strats.push(stratify('morphT (early<=4 vs late>4)', morphT, v => v <= 4));
strats.push(stratify('environment (field=1 vs group/cluster)', envCode, v => v > 1));
strats.push(stratify('logSigma0 (median split)', logSig0, v => v >= median(logSig0)));
strats.push(stratify('Υ★⊥ (above/below median)', upsPerp, v => v >= median(upsPerp)));

const sigStrats = strats.filter(s => s.sig);
console.log('  SUMMARY: ' + sigStrats.length + '/' + strats.length + ' splits show significant (|t|>2) bias');
if (sigStrats.length > 0) console.log('  Significant: ' + sigStrats.map(s=>s.name).join(', '));
console.log();

// ═════════ 75b: RESIDUAL CLUSTERING ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  75b: RESIDUAL CLUSTERING');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Simple k-means in standardized M5 space
const stdize = (arr) => { const m=mean(arr), s=sd(arr); return arr.map(v=>(v-m)/(s||1)); };
const feats = [stdize(logMHI), stdize(logMhost), stdize(logSig0), stdize(logMR), stdize(upsPerp)];
const points = gals45.map((_,i) => feats.map(f => f[i]));

function kmeans(pts, k, maxIter=100) {
  const n=pts.length, d=pts[0].length;
  let centers = pts.slice(0,k).map(p=>[...p]);
  let labels = new Array(n).fill(0);
  for (let iter=0; iter<maxIter; iter++) {
    let changed = false;
    for (let i=0; i<n; i++) {
      let best=0, bestD=Infinity;
      for (let c=0; c<k; c++) {
        let dist=0;
        for (let j=0; j<d; j++) dist += (pts[i][j]-centers[c][j])**2;
        if (dist < bestD) { bestD=dist; best=c; }
      }
      if (labels[i] !== best) { labels[i]=best; changed=true; }
    }
    if (!changed) break;
    for (let c=0; c<k; c++) {
      const members = pts.filter((_,i)=>labels[i]===c);
      if (members.length === 0) continue;
      for (let j=0; j<d; j++) centers[c][j] = mean(members.map(m=>m[j]));
    }
  }
  return labels;
}

for (const k of [2, 3]) {
  const labels = kmeans(points, k);
  console.log('  k=' + k + ' clustering:');
  for (let c = 0; c < k; c++) {
    const members = gals45.map((_,i)=>i).filter(i=>labels[i]===c);
    const clResid = members.map(i=>resid[i]);
    const clAbsResid = members.map(i=>absResid[i]);
    console.log('    Cluster ' + c + ' (n=' + members.length + '): mean resid = ' + (mean(clResid)>0?'+':'') + mean(clResid).toFixed(3) + ', SD = ' + sd(clResid).toFixed(3) + ', mean |resid| = ' + mean(clAbsResid).toFixed(3));
  }
  // F-test: between-cluster vs within-cluster variance of residuals
  const grandMean = mean(resid);
  let ssB = 0, ssW = 0;
  for (let c = 0; c < k; c++) {
    const members = gals45.map((_,i)=>i).filter(i=>labels[i]===c);
    const clMean = mean(members.map(i=>resid[i]));
    ssB += members.length * (clMean - grandMean)**2;
    ssW += members.map(i=>(resid[i]-clMean)**2).reduce((s,v)=>s+v,0);
  }
  const F = (ssB/(k-1)) / (ssW/(N-k));
  console.log('    F-stat (resid by cluster) = ' + F.toFixed(2) + (F>3.2?' *** significant':' (not significant)'));
  console.log();
}

// ═════════ 75c: HETEROSCEDASTICITY ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  75c: HETEROSCEDASTICITY TEST');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const residSq = resid.map(r => r*r);
const hetTests = [
  {name:'predicted a₀', vals:predicted},
  {name:'logMHI', vals:logMHI},
  {name:'logMhost', vals:logMhost},
  {name:'logΣ₀', vals:logSig0},
  {name:'logMeanRun', vals:logMR},
  {name:'Υ★⊥', vals:upsPerp},
  {name:'distance', vals:dist},
  {name:'morphT', vals:morphT},
];

console.log('  |resid|² vs predictor (Breusch-Pagan proxy):');
console.log('  ┌──────────────────┬────────┬────────┬──────────┐');
console.log('  │ Predictor        │ Pearson│ Spear  │ Flag     │');
console.log('  ├──────────────────┼────────┼────────┼──────────┤');
const hetResults = {};
hetTests.forEach(h => {
  const rP = pearsonCorr(residSq, h.vals);
  const rS = spearmanCorr(absResid, h.vals);
  const flag = Math.abs(rP) > 0.3 ? '⚠️ hetero' : 'OK';
  hetResults[h.name] = {pearson:rP, spearman:rS, flag};
  console.log('  │ ' + h.name.padEnd(16) + ' │ ' + (rP>0?'+':'')+rP.toFixed(3).padStart(5) + '  │ ' + (rS>0?'+':'')+rS.toFixed(3).padStart(5) + '  │ ' + flag.padEnd(8) + ' │');
});
console.log('  └──────────────────┴────────┴────────┴──────────┘');
console.log();

// Simple Breusch-Pagan: regress resid² on all M5 predictors
const fBP = ols(residSq, X_M5);
const bpR2 = 1 - fBP.resid.reduce((s,r)=>s+r*r,0) / residSq.reduce((s,v)=>s+(v-mean(residSq))**2,0);
const bpStat = N * bpR2;
console.log('  Breusch-Pagan test: nR² = ' + bpStat.toFixed(2) + ' (chi²(5), critical ~11.1 at 5%)');
console.log('  ' + (bpStat > 11.1 ? '⚠️ SIGNIFICANT heteroscedasticity' : 'No significant heteroscedasticity'));
console.log();

// ═════════ 75d: LOCAL FAILURE MAP ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  75d: LOCAL FAILURE MAP');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const sorted = gals45.map((g,i) => ({name:g.name, resid:resid[i], absResid:absResid[i], morphT:morphT[i], logMhost:logMhost[i], logMHI:logMHI[i], dist:dist[i]}));
sorted.sort((a,b) => b.absResid - a.absResid);

console.log('  TOP 10 LARGEST |RESIDUALS| (where M5 struggles most):');
console.log('  ' + 'Rank'.padEnd(5) + 'Galaxy'.padEnd(14) + 'resid    |resid|   morphT  logMhost  dist');
console.log('  ' + '-'.repeat(75));
sorted.slice(0,10).forEach((g,i) => {
  console.log('  ' + (i+1+'').padEnd(5) + g.name.padEnd(14) + (g.resid>0?'+':'') + g.resid.toFixed(3).padStart(6) + '   ' + g.absResid.toFixed(3).padStart(6) + '    ' + g.morphT.toString().padStart(4) + '    ' + g.logMhost.toFixed(2).padStart(6) + '    ' + g.dist.toFixed(1).padStart(5));
});
console.log();

console.log('  TOP 10 SMALLEST |RESIDUALS| (where M5 works best):');
sorted.sort((a,b) => a.absResid - b.absResid);
sorted.slice(0,10).forEach((g,i) => {
  console.log('  ' + (i+1+'').padEnd(5) + g.name.padEnd(14) + (g.resid>0?'+':'') + g.resid.toFixed(3).padStart(6) + '   ' + g.absResid.toFixed(3).padStart(6) + '    ' + g.morphT.toString().padStart(4) + '    ' + g.logMhost.toFixed(2).padStart(6) + '    ' + g.dist.toFixed(1).padStart(5));
});
console.log();

// Analyze top/bottom 10 patterns
sorted.sort((a,b) => b.absResid - a.absResid);
const worst10 = sorted.slice(0,10);
const best10 = sorted.slice(-10);
console.log('  Pattern analysis:');
console.log('    Worst-10 mean morphT = ' + mean(worst10.map(g=>g.morphT)).toFixed(1) + ', Best-10 = ' + mean(best10.map(g=>g.morphT)).toFixed(1));
console.log('    Worst-10 mean logMhost = ' + mean(worst10.map(g=>g.logMhost)).toFixed(2) + ', Best-10 = ' + mean(best10.map(g=>g.logMhost)).toFixed(2));
console.log('    Worst-10 mean dist = ' + mean(worst10.map(g=>g.dist)).toFixed(1) + ', Best-10 = ' + mean(best10.map(g=>g.dist)).toFixed(1));
console.log('    Worst-10 sign balance: ' + worst10.filter(g=>g.resid>0).length + ' positive, ' + worst10.filter(g=>g.resid<0).length + ' negative');
console.log();

// ═════════ 75e: LATENT-STATE SHAPE ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  75e: LATENT-STATE SHAPE (Distribution Analysis)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Normality tests
const sortedResid = [...resid].sort((a,b)=>a-b);
const n = N;

// Skewness
const m3 = resid.reduce((s,r) => s + ((r-mean(resid))/residSD)**3, 0) / n;
// Kurtosis (excess)
const m4 = resid.reduce((s,r) => s + ((r-mean(resid))/residSD)**4, 0) / n - 3;

console.log('  Distribution moments:');
console.log('    Skewness = ' + m3.toFixed(3) + (Math.abs(m3) > 0.5 ? ' ⚠️ asymmetric' : ' (symmetric)'));
console.log('    Excess kurtosis = ' + m4.toFixed(3) + (m4 > 1 ? ' ⚠️ heavy tails' : m4 < -1 ? ' ⚠️ light tails' : ' (near-Gaussian)'));
console.log();

// Jarque-Bera
const JB = n/6 * (m3**2 + m4**2/4);
console.log('    Jarque-Bera = ' + JB.toFixed(2) + ' (critical ~5.99 at 5%)');
console.log('    ' + (JB > 5.99 ? '⚠️ Rejects normality' : 'Cannot reject normality'));
console.log();

// Histogram
const bins = 8;
const rMin = Math.min(...resid), rMax = Math.max(...resid);
const binW = (rMax - rMin) / bins;
const hist = new Array(bins).fill(0);
resid.forEach(r => { const b = Math.min(bins-1, Math.floor((r-rMin)/binW)); hist[b]++; });
console.log('  Histogram of residuals:');
const maxH = Math.max(...hist);
for (let b = 0; b < bins; b++) {
  const lo = rMin + b * binW;
  const hi = lo + binW;
  const bar = '#'.repeat(Math.round(hist[b] / maxH * 30));
  console.log('    [' + (lo>0?'+':'') + lo.toFixed(2) + ', ' + (hi>0?'+':'') + hi.toFixed(2) + ') ' + hist[b].toString().padStart(2) + ' ' + bar);
}
console.log();

// QQ-like: compare sorted residuals vs expected normal quantiles
console.log('  QQ-check (sorted residuals vs expected Gaussian):');
const quantiles = [0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95];
const zScores = [-1.645, -1.282, -0.674, 0, 0.674, 1.282, 1.645];
console.log('    Quantile  Expected   Observed   Ratio');
quantiles.forEach((q,qi) => {
  const idx = Math.min(n-1, Math.max(0, Math.round(q*(n-1))));
  const expected = zScores[qi] * residSD;
  const observed = sortedResid[idx];
  const ratio = expected !== 0 ? observed/expected : (observed===0?1:Infinity);
  console.log('    ' + (q*100).toFixed(0).padStart(4) + '%     ' + (expected>0?'+':'') + expected.toFixed(3).padStart(6) + '    ' + (observed>0?'+':'') + observed.toFixed(3).padStart(6) + '    ' + ratio.toFixed(2));
});
console.log();

// Bimodality: Hartigan's dip-like test (simplified)
// Check if distribution has two modes
const halfN = Math.floor(n/2);
const loHalf = sortedResid.slice(0, halfN);
const hiHalf = sortedResid.slice(halfN);
const gapBetween = hiHalf[0] - loHalf[loHalf.length-1];
const loSpread = loHalf[loHalf.length-1] - loHalf[0];
const hiSpread = hiHalf[hiHalf.length-1] - hiHalf[0];
console.log('  Bimodality check:');
console.log('    Gap between halves: ' + gapBetween.toFixed(3) + ' dex');
console.log('    Lower half spread: ' + loSpread.toFixed(3) + ', Upper half spread: ' + hiSpread.toFixed(3));
console.log('    ' + (gapBetween > 0.5*residSD ? '⚠️ Possible bimodality' : 'No bimodality detected'));
console.log();

// ═════════ FINAL VERDICT ═════════
console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 75: FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const anyStratSig = sigStrats.length > 0;
const anyHetero = Object.values(hetResults).some(h => h.flag.includes('hetero'));
const bpSig = bpStat > 11.1;
const normalResid = JB < 5.99;
const noBimodal = gapBetween <= 0.5 * residSD;

let verdict;
if (anyStratSig || bpSig) {
  verdict = 'PARTIAL';
} else if (!anyStratSig && !bpSig && normalResid && noBimodal) {
  verdict = 'CONFIRMED-FRONTIER';
} else {
  verdict = 'PARTIAL';
}

console.log('  ┌─────────────────────────────────────┬────────┐');
console.log('  │ Criterion                           │ Result │');
console.log('  ├─────────────────────────────────────┼────────┤');
console.log('  │ Any stratification bias (|t|>2)?    │ ' + (anyStratSig?'YES':'NO').padEnd(6) + ' │');
console.log('  │ Any heteroscedasticity (|r|>0.3)?   │ ' + (anyHetero?'YES':'NO').padEnd(6) + ' │');
console.log('  │ Breusch-Pagan significant?          │ ' + (bpSig?'YES':'NO').padEnd(6) + ' │');
console.log('  │ Residuals Gaussian-like?            │ ' + (normalResid?'YES':'NO').padEnd(6) + ' │');
console.log('  │ No bimodality?                      │ ' + (noBimodal?'YES':'NO').padEnd(6) + ' │');
console.log('  └─────────────────────────────────────┴────────┘');
console.log();
console.log('  VERDICT: ' + verdict);
console.log();

if (verdict === 'CONFIRMED-FRONTIER') {
  console.log('  M5 residuals show no structure: no subpopulation bias,');
  console.log('  no heteroscedasticity, Gaussian-shaped, no bimodality.');
  console.log('  The remaining ~' + residRMS.toFixed(2) + ' dex appears to be a genuine');
  console.log('  irreducible frontier — possibly intrinsic scatter plus');
  console.log('  measurement noise.');
} else if (verdict === 'PARTIAL') {
  console.log('  Some residual structure remains:');
  if (anyStratSig) console.log('  - Stratification bias in: ' + sigStrats.map(s=>s.name).join(', '));
  if (anyHetero) console.log('  - Heteroscedasticity detected');
  if (bpSig) console.log('  - Breusch-Pagan test significant');
  console.log('  The residual is not yet fully featureless, suggesting');
  console.log('  potential for further refinement or that specific');
  console.log('  subpopulations drive the remaining scatter.');
}

console.log();
console.log('  Residual RMS: ' + residRMS.toFixed(3) + ' dex');
console.log('  This is the current irreducible frontier of M5.');

const output = {
  phase:'75', title:'Residual Physics Frontier',
  residualStats: {rms:residRMS, sd:residSD, skewness:m3, kurtosis:m4, JB, normal:normalResid},
  stratification: strats.map(s=>({name:s.name,nHi:s.nHi,nLo:s.nLo,meanHi:s.meanHi,meanLo:s.meanLo,sdHi:s.sdHi,sdLo:s.sdLo,t:s.t,sig:s.sig})),
  heteroscedasticity: hetResults,
  breuscPagan: {nR2:bpStat, significant:bpSig},
  localFailure: {worst10:sorted.slice(0,10).map(g=>({name:g.name,resid:g.resid})), best10:sorted.slice(-10).map(g=>({name:g.name,resid:g.resid}))},
  distributionShape: {gaussian:normalResid, bimodal:!noBimodal, skewness:m3, kurtosis:m4},
  verdict
};
fs.writeFileSync('public/phase75-residual-frontier.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase75-residual-frontier.json');
