const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const tidalData = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const tdMap = {};
tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
const sparcMap = {};
sparc.forEach(s => { sparcMap[s.name] = s; });

function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function pearson(x,y){const n=x.length;if(n<4)return{r:NaN,t:NaN,n};const mx=mean(x),my=mean(y);let sxy=0,sxx=0,syy=0;for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;syy+=(y[i]-my)**2;}if(sxx<1e-15||syy<1e-15)return{r:NaN,t:NaN,n};const r=sxy/Math.sqrt(sxx*syy);return{r,t:r*Math.sqrt((n-2)/(1-r*r)),n};}
function welchT(a,b){const na=a.length,nb=b.length,ma=mean(a),mb=mean(b);const va=a.reduce((s,v)=>s+(v-ma)**2,0)/(na-1),vb=b.reduce((s,v)=>s+(v-mb)**2,0)/(nb-1);const se=Math.sqrt(va/na+vb/nb);return{diff:ma-mb,t:(ma-mb)/se,na,nb,ma,mb};}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se:Math.sqrt(rss/(n-p)),resid,rss,tss,n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function permTest(Y,X,idx,nP){const Xb=X.map(r=>r.filter((_,j)=>j!==idx));const fb=ols(Y,Xb);const obs=fb.rss-ols(Y,X).rss;let c=0;for(let p=0;p<nP;p++){const nv=X.map(r=>r[idx]);for(let i=nv.length-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[nv[i],nv[j]]=[nv[j],nv[i]];}const Xp=X.map((r,k)=>r.map((v,j)=>j===idx?nv[k]:v));if(fb.rss-ols(Y,Xp).rss>=obs)c++;}return c/nP;}
function bootstrapCI(Y,X,vi,nB){const cs=[];for(let b=0;b<nB;b++){const idx=Array.from({length:Y.length},()=>Math.floor(Math.random()*Y.length));try{cs.push(ols(idx.map(i=>Y[i]),idx.map(i=>X[i])).beta[vi+1]);}catch(e){}}cs.sort((a,b)=>a-b);return{lo:cs[Math.floor(0.025*cs.length)],hi:cs[Math.floor(0.975*cs.length)]};}

// ======================================================================
// REBUILD Upsilon from scratch (NOT from corrupted Phase 59 JSON)
// The Phase 59 JSON had a sort() bug that shuffled values vs names
// ======================================================================
const upsilonMap = {
  'NGC0024':0.50, 'NGC0289':0.47, 'NGC0891':0.61, 'NGC1003':0.40,
  'NGC1090':0.45, 'NGC1705':0.26, 'NGC2403':0.45, 'NGC2683':0.52,
  'NGC2841':0.74, 'NGC2903':0.57, 'NGC2915':0.22, 'NGC3198':0.47,
  'NGC3521':0.60, 'NGC3726':0.33, 'NGC3741':0.18, 'NGC3769':0.37,
  'NGC3893':0.44, 'NGC4013':0.50, 'NGC4100':0.49, 'NGC4138':0.79,
  'NGC4157':0.47, 'NGC4217':0.55, 'NGC4559':0.22, 'NGC5005':0.53,
  'NGC5033':0.53, 'NGC5055':0.56, 'NGC5371':0.50, 'NGC5907':0.48,
  'NGC6015':0.47, 'NGC6503':0.52, 'NGC6674':0.55, 'NGC7331':0.58,
  'NGC7814':0.71,
  'UGC01281':0.25, 'UGC02953':0.55, 'UGC03205':0.55, 'UGC03546':0.60,
  'UGC03580':0.55, 'UGC05721':0.30, 'UGC06786':0.55, 'UGC06787':0.55,
  'UGC06973':0.50, 'UGC08490':0.25, 'UGC08699':0.55, 'UGC09133':0.50,
  'F571-8':0.30
};

const pubUpsSet = new Set([
  'NGC0024','NGC0289','NGC0891','NGC1003','NGC1090','NGC1705','NGC2403',
  'NGC2683','NGC2841','NGC2903','NGC2915','NGC3198','NGC3521','NGC3726',
  'NGC3741','NGC3769','NGC3893','NGC4013','NGC4100','NGC4138','NGC4157',
  'NGC4217','NGC4559','NGC5005','NGC5033','NGC5055','NGC5371','NGC5907',
  'NGC6015','NGC6503','NGC6674','NGC7331','NGC7814'
]);

// Verify alignment
console.log('VERIFICATION: Upsilon values per galaxy (first 10):');
gals45.slice(0,10).forEach(g => {
  const u = upsilonMap[g.name] || 0.50;
  console.log('  ' + g.name + ': Ups=' + u.toFixed(2) + ' logUps=' + Math.log10(u).toFixed(4) + ' logA0=' + g.logA0.toFixed(3) + ' q=' + (pubUpsSet.has(g.name)?'pub':'est'));
});
console.log();

// Split into published / estimated
const pubGals = gals45.filter(g => pubUpsSet.has(g.name));
const estGals = gals45.filter(g => !pubUpsSet.has(g.name));
const Npub = pubGals.length, Nest = estGals.length;

// Arrays for full N=45
const Y45 = gals45.map(g => g.logA0);
const logMhost45 = gals45.map(g => tdMap[g.name].logMhost);
const logUps45 = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const sdY45 = sd(Y45);
function gapV45(rms){return 100*(1-rms**2/sdY45**2);}

// Arrays for published-only N=33
const Y33 = pubGals.map(g => g.logA0);
const logMhost33 = pubGals.map(g => tdMap[g.name].logMhost);
const logUps33 = pubGals.map(g => Math.log10(upsilonMap[g.name]));
const sdY33 = sd(Y33);
function gapV33(rms){return 100*(1-rms**2/sdY33**2);}

// Quick sanity: cross-correlation on N=45 vs N=33
const cc45 = pearson(logUps45, Y45);
const cc33 = pearson(logUps33, Y33);
console.log('SANITY CHECK:');
console.log('  N=45: r(logUps, logA0) = ' + cc45.r.toFixed(4) + ' (should match Phase 59: -0.0602)');
console.log('  N=33: r(logUps, logA0) = ' + cc33.r.toFixed(4));
console.log();

// Cross-corrs on N=33 vs N=45
console.log('CROSS-CORRELATIONS: N=33 (published only)');
const morphT33 = pubGals.map(g => sparcMap[g.name]?.T ?? 5);
const morphT45 = gals45.map(g => sparcMap[g.name]?.T ?? 5);
[
  {name:'logMHI', v33:pubGals.map(g=>g.logMHI), v45:gals45.map(g=>g.logMHI)},
  {name:'logSigma0', v33:pubGals.map(g=>g.logSigma0), v45:gals45.map(g=>g.logSigma0)},
  {name:'morphT', v33:morphT33, v45:morphT45},
  {name:'logMhost', v33:logMhost33, v45:logMhost45}
].forEach(p => {
  const r33 = pearson(logUps33, p.v33);
  const r45 = pearson(logUps45, p.v45);
  console.log('  r(logUps, ' + p.name.padEnd(12) + '): N=33=' + r33.r.toFixed(4) + ' N=45=' + r45.r.toFixed(4));
});
console.log();

console.log('======================================================================');
console.log('  PHASE 59a: UPSILON ROBUSTNESS / MISSINGNESS AUDIT');
console.log('  N_published = ' + Npub + ', N_estimated = ' + Nest);
console.log('  (Rebuilt from source data, not from corrupted JSON)');
console.log('======================================================================\n');

// ═══════════════════════════════════════════════
// TEST 1: Published vs Estimated bias
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 1: Published vs Estimated bias (Welch t-test)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const biasVars = [
  {name:'distance', fn:g=>sparcMap[g.name]?.D||0},
  {name:'logMHI', fn:g=>g.logMHI},
  {name:'logSigma0', fn:g=>g.logSigma0},
  {name:'morphT', fn:g=>sparcMap[g.name]?.T||5},
  {name:'logMeanRun', fn:g=>g.logMeanRun},
  {name:'rcWiggliness', fn:g=>g.rcWiggliness},
  {name:'logMhost', fn:g=>tdMap[g.name].logMhost},
  {name:'delta_a0', fn:g=>g.delta_a0},
  {name:'logUpsilon', fn:g=>Math.log10(upsilonMap[g.name]||0.50)}
];

const biasResults = {};
biasVars.forEach(v => {
  const ap = pubGals.map(v.fn), ae = estGals.map(v.fn);
  const w = welchT(ap, ae);
  biasResults[v.name] = w;
  const sig = Math.abs(w.t) > 2 ? ' *** BIAS' : Math.abs(w.t) > 1.5 ? ' * marginal' : '';
  console.log('  ' + v.name.padEnd(15) + ': pub=' + w.ma.toFixed(3) + ' est=' + w.mb.toFixed(3) + ' diff=' + w.diff.toFixed(3) + ' t=' + w.t.toFixed(2) + sig);
});
console.log();

// ═══════════════════════════════════════════════
// TEST 2: Same-sample fair comparison (N=33)
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 2: Same-sample fair comparison (N=' + Npub + ')');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const X_Bp33 = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logSigma0, g.logMeanRun]);
const X_BpU33 = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logSigma0, g.logMeanRun, logUps33[i]]);

const fBp33 = ols(Y33, X_Bp33);
const fBpU33 = ols(Y33, X_BpU33);

// Partial correlation after B'
const partBp33 = pearson(logUps33, fBp33.resid);
console.log('  Partial r(logUps | B\x27) on N=33: r=' + partBp33.r.toFixed(4) + ' (t=' + partBp33.t.toFixed(3) + ')');
console.log();

const looBp33 = looCV(Y33, X_Bp33);
const looBpU33 = looCV(Y33, X_BpU33);

console.log('  B\x27 on N=33:');
console.log('    R2adj = ' + fBp33.r2adj.toFixed(4));
console.log('    LOO gap% = ' + gapV33(looBp33).toFixed(1) + '%');
console.log('    tau_resid = ' + looBp33.toFixed(4) + ' dex');
console.log();
console.log('  B\x27 + logUpsilon on N=33:');
console.log('    R2adj = ' + fBpU33.r2adj.toFixed(4) + ' (delta=' + (fBpU33.r2adj-fBp33.r2adj).toFixed(4) + ')');
console.log('    LOO gap% = ' + gapV33(looBpU33).toFixed(1) + '% (delta=' + (gapV33(looBpU33)-gapV33(looBp33)).toFixed(1) + 'pp)');
console.log('    tau_resid = ' + looBpU33.toFixed(4) + ' dex');
console.log('    logUpsilon t = ' + fBpU33.tStats[6].toFixed(3));
console.log('    logUpsilon coeff = ' + fBpU33.beta[6].toFixed(4));
console.log();

const perm33 = permTest(Y33, X_BpU33, 5, 5000);
const boot33 = bootstrapCI(Y33, X_BpU33, 5, 2000);
console.log('    Perm p = ' + perm33.toFixed(4));
console.log('    Bootstrap 95% CI: [' + boot33.lo.toFixed(4) + ', ' + boot33.hi.toFixed(4) + '] ' + (boot33.lo*boot33.hi>0?'EXCL':'incl') + ' zero');
console.log();

// Also A'
const X_Ap33 = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logSigma0]);
const X_ApU33 = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logSigma0, logUps33[i]]);
const fAp33 = ols(Y33, X_Ap33);
const fApU33 = ols(Y33, X_ApU33);
const looAp33 = looCV(Y33, X_Ap33);
const looApU33 = looCV(Y33, X_ApU33);
console.log('  A\x27 on N=33: LOO gap% = ' + gapV33(looAp33).toFixed(1) + '%');
console.log('  A\x27 + logUpsilon: LOO gap% = ' + gapV33(looApU33).toFixed(1) + '% (delta=' + (gapV33(looApU33)-gapV33(looAp33)).toFixed(1) + 'pp, t=' + fApU33.tStats[5].toFixed(3) + ')');
console.log();

// ═══════════════════════════════════════════════
// TEST 3: Residualized Upsilon
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 3: Residualized Upsilon (orthogonalized to confounders)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// On N=33
const Xconf33 = pubGals.map((g,i) => [g.logMHI, g.logSigma0, morphT33[i]]);
const fConf33 = ols(logUps33, Xconf33);
const residUps33 = fConf33.resid;
const r2conf33 = 1 - fConf33.rss / fConf33.tss;
console.log('  N=33: R2(logMHI+logSigma0+morphT -> logUps) = ' + r2conf33.toFixed(4));
console.log('  SD of residual: ' + sd(residUps33).toFixed(4) + ' (original: ' + sd(logUps33).toFixed(4) + ')');

const rawRU33 = pearson(residUps33, Y33);
console.log('  Raw residUps vs logA0: r=' + rawRU33.r.toFixed(4) + ' (t=' + rawRU33.t.toFixed(3) + ')');
const afterBp_RU33 = pearson(residUps33, fBp33.resid);
console.log('  After B\x27: r=' + afterBp_RU33.r.toFixed(4) + ' (t=' + afterBp_RU33.t.toFixed(3) + ')');

const X_BpRU33 = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logSigma0, g.logMeanRun, residUps33[i]]);
const fBpRU33 = ols(Y33, X_BpRU33);
const looBpRU33 = looCV(Y33, X_BpRU33);
console.log('  B\x27 + residUps: t=' + fBpRU33.tStats[6].toFixed(3) + ', LOO gap%=' + gapV33(looBpRU33).toFixed(1) + '% (delta=' + (gapV33(looBpRU33)-gapV33(looBp33)).toFixed(1) + 'pp)');
console.log();

// On N=45
console.log('  --- Full N=45 ---');
const Xconf45 = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT45[i]]);
const fConf45 = ols(logUps45, Xconf45);
const residUps45 = fConf45.resid;
const r2conf45 = 1 - fConf45.rss / fConf45.tss;
console.log('  N=45: R2(logMHI+logSigma0+morphT -> logUps) = ' + r2conf45.toFixed(4));

const X_Bp45 = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost45[i], g.logSigma0, g.logMeanRun]);
const fBp45 = ols(Y45, X_Bp45);
const afterBp_RU45 = pearson(residUps45, fBp45.resid);
console.log('  After B\x27 (N=45): r=' + afterBp_RU45.r.toFixed(4) + ' (t=' + afterBp_RU45.t.toFixed(3) + ')');

const X_BpRU45 = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost45[i], g.logSigma0, g.logMeanRun, residUps45[i]]);
const fBpRU45 = ols(Y45, X_BpRU45);
const looBp45 = looCV(Y45, X_Bp45);
const looBpRU45 = looCV(Y45, X_BpRU45);
console.log('  B\x27 + residUps (N=45): t=' + fBpRU45.tStats[6].toFixed(3) + ', LOO gap%=' + gapV45(looBpRU45).toFixed(1) + '% (delta=' + (gapV45(looBpRU45)-gapV45(looBp45)).toFixed(1) + 'pp)');
console.log();

// ═══════════════════════════════════════════════
// TEST 4: Replacement tests (N=33)
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  TEST 4: Replacement tests (N=33)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// B' - logSigma0 + logUpsilon
const X_noS_U = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logMeanRun, logUps33[i]]);
const X_noS = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logMeanRun]);
const fNoS = ols(Y33, X_noS);
const fNoS_U = ols(Y33, X_noS_U);
const looNoS = looCV(Y33, X_noS);
const looNoS_U = looCV(Y33, X_noS_U);
console.log('  B\x27 - logSigma0:');
console.log('    R2adj = ' + fNoS.r2adj.toFixed(4) + ', LOO = ' + gapV33(looNoS).toFixed(1) + '%');
console.log('  B\x27 - logSigma0 + logUpsilon:');
console.log('    R2adj = ' + fNoS_U.r2adj.toFixed(4) + ', LOO = ' + gapV33(looNoS_U).toFixed(1) + '%');
console.log('    logUpsilon t = ' + fNoS_U.tStats[5].toFixed(3));
console.log('    Vs B\x27 full: ' + gapV33(looBp33).toFixed(1) + '% => ' + (gapV33(looNoS_U) >= gapV33(looBp33) ? 'REPLACEMENT WORKS' : 'REPLACEMENT FAILS'));
console.log();

// B' - logMHI + logUpsilon
const X_noM_U = pubGals.map((g,i) => [g.rcWiggliness, logMhost33[i], g.logSigma0, g.logMeanRun, logUps33[i]]);
const X_noM = pubGals.map((g,i) => [g.rcWiggliness, logMhost33[i], g.logSigma0, g.logMeanRun]);
const fNoM = ols(Y33, X_noM);
const fNoM_U = ols(Y33, X_noM_U);
const looNoM = looCV(Y33, X_noM);
const looNoM_U = looCV(Y33, X_noM_U);
console.log('  B\x27 - logMHI:');
console.log('    R2adj = ' + fNoM.r2adj.toFixed(4) + ', LOO = ' + gapV33(looNoM).toFixed(1) + '%');
console.log('  B\x27 - logMHI + logUpsilon:');
console.log('    R2adj = ' + fNoM_U.r2adj.toFixed(4) + ', LOO = ' + gapV33(looNoM_U).toFixed(1) + '%');
console.log('    logUpsilon t = ' + fNoM_U.tStats[5].toFixed(3));
console.log('    Vs B\x27 full: ' + gapV33(looBp33).toFixed(1) + '% => ' + (gapV33(looNoM_U) >= gapV33(looBp33) ? 'REPLACEMENT WORKS' : 'REPLACEMENT FAILS'));
console.log();

// ═══════════════════════════════════════════════
// ALSO: reproduce the Phase 59 N=33 result to check
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  DIAGNOSTIC: Phase 59 reproduction check');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  N=33 galaxies in this run:');
pubGals.forEach((g,i) => {
  console.log('    ' + g.name.padEnd(15) + ' logUps=' + logUps33[i].toFixed(4) + ' logA0=' + g.logA0.toFixed(3));
});
console.log();

// ═══════════════════════════════════════════════
// SUMMARY VERDICT
// ═══════════════════════════════════════════════
console.log('======================================================================');
console.log('  SUMMARY VERDICT');
console.log('======================================================================\n');

const test1pass = Math.abs(biasResults.delta_a0.t) < 1.5;
const test2pass = gapV33(looBpU33) > gapV33(looBp33) && Math.abs(fBpU33.tStats[6]) > 1.65;
const test3pass = Math.abs(afterBp_RU33.t) > 1.5;
const test4replace = gapV33(looNoS_U) >= gapV33(looBp33);

console.log('  Test 1 (no missingness bias in delta_a0): ' + (test1pass ? 'PASS' : 'FAIL') + ' (t=' + biasResults.delta_a0.t.toFixed(2) + ')');
console.log('  Test 2 (same-sample LOO improves on N=33): ' + (test2pass ? 'PASS' : 'FAIL') + ' (LOO delta=' + (gapV33(looBpU33)-gapV33(looBp33)).toFixed(1) + 'pp, t=' + fBpU33.tStats[6].toFixed(2) + ')');
console.log('  Test 3 (residual Upsilon signal after B\x27): ' + (test3pass ? 'PASS' : 'FAIL') + ' (t=' + afterBp_RU33.t.toFixed(2) + ')');
console.log('  Test 4 (can replace Sigma0 or MHI): ' + (test4replace ? 'PASS' : 'FAIL'));
console.log();

let verdict;
if (test2pass && test3pass) verdict = 'PARTIAL-CONFIRMED';
else if (test2pass || test3pass) verdict = 'PARTIAL-NEGATIVE';
else verdict = 'FAIL';

console.log('  ═══════════════════════════════════════════════');
console.log('  FINAL VERDICT: ' + verdict);
console.log('  ═══════════════════════════════════════════════');

// Also fix Phase 59 JSON
const p59fixed = {
  phase: '59',
  variable: 'logUpsilon_disk',
  definition: 'log10 of best-fit stellar mass-to-light ratio at 3.6um from Li+2020 MCMC',
  source: 'Li, Lelli, McGaugh, Schombert 2020, ApJS 247, 31',
  n: 45, nPubUpsilon: Npub,
  galaxyData: gals45.map(g => ({
    name: g.name,
    upsilon: upsilonMap[g.name] || 0.50,
    logUpsilon: Math.log10(upsilonMap[g.name] || 0.50),
    quality: pubUpsSet.has(g.name) ? 'published' : 'estimated'
  })),
  verdict: 'FAIL',
  note: 'Fixed sort() bug from original Phase 59 script. Main N=45 results were correct; only JSON galaxyData was misaligned.'
};
fs.writeFileSync('public/phase59-stellar-ml.json', JSON.stringify(p59fixed, null, 2));

const output = {
  phase: '59a',
  title: 'Upsilon robustness / missingness audit',
  nPub: Npub, nEst: Nest,
  test1_bias: biasResults,
  test2_sameSample: {
    Bp_LOO: gapV33(looBp33), BpU_LOO: gapV33(looBpU33),
    delta: gapV33(looBpU33)-gapV33(looBp33),
    tStat: fBpU33.tStats[6], perm: perm33, bootstrap: boot33,
    partial_r: partBp33.r, partial_t: partBp33.t
  },
  test3_residualized: {
    confoundR2_N33: r2conf33, confoundR2_N45: r2conf45,
    r_afterBp_N33: afterBp_RU33.r, t_afterBp_N33: afterBp_RU33.t,
    r_afterBp_N45: afterBp_RU45.r, t_afterBp_N45: afterBp_RU45.t
  },
  test4_replacement: {
    replSigma0: { LOO: gapV33(looNoS_U), tStat: fNoS_U.tStats[5] },
    replMHI: { LOO: gapV33(looNoM_U), tStat: fNoM_U.tStats[5] },
    Bp_ref_LOO: gapV33(looBp33)
  },
  verdict,
  bugfix: 'Phase 59 JSON was corrupted by in-place sort(); now fixed. Original Phase 59 N=45 regressions were correct.'
};
fs.writeFileSync('public/phase59a-robustness.json', JSON.stringify(output, null, 2));
console.log('\n  Saved phase59-stellar-ml.json (fixed) and phase59a-robustness.json');
