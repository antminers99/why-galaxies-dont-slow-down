const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const tidalData = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const tdMap = {};
tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
const N45 = gals45.length;
const logMhost45 = gals45.map(g => tdMap[g.name].logMhost);

const sparcMap = {};
sparc.forEach(s => { sparcMap[s.name] = s; });

function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function pearson(x,y){const n=x.length;if(n<4)return{r:NaN,t:NaN,n};const mx=mean(x),my=mean(y);let sxy=0,sxx=0,syy=0;for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;syy+=(y[i]-my)**2;}if(sxx<1e-15||syy<1e-15)return{r:NaN,t:NaN,n};const r=sxy/Math.sqrt(sxx*syy);return{r,t:r*Math.sqrt((n-2)/(1-r*r)),n};}
function spearmanR(x,y){const n=x.length;function rank(a){const s=a.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);const r=new Array(n);for(let i=0;i<n;){let j=i;while(j<n&&s[j].v===s[i].v)j++;const avg=(i+j-1)/2+1;for(let k=i;k<j;k++)r[s[k].i]=avg;i=j;}return r;}const rx=rank(x),ry=rank(y);return pearson(rx,ry).r;}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se:Math.sqrt(rss/(n-p)),resid,rss,tss,n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function permTest(Y,X,idx,nP){const Xb=X.map(r=>r.filter((_,j)=>j!==idx));const fb=ols(Y,Xb);const obs=fb.rss-ols(Y,X).rss;let c=0;for(let p=0;p<nP;p++){const nv=X.map(r=>r[idx]);for(let i=nv.length-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[nv[i],nv[j]]=[nv[j],nv[i]];}const Xp=X.map((r,k)=>r.map((v,j)=>j===idx?nv[k]:v));if(fb.rss-ols(Y,Xp).rss>=obs)c++;}return c/nP;}
function bootstrapCI(Y,X,vi,nB){const cs=[];for(let b=0;b<nB;b++){const idx=Array.from({length:Y.length},()=>Math.floor(Math.random()*Y.length));try{cs.push(ols(idx.map(i=>Y[i]),idx.map(i=>X[i])).beta[vi+1]);}catch(e){}}cs.sort((a,b)=>a-b);return{lo:cs[Math.floor(0.025*cs.length)],hi:cs[Math.floor(0.975*cs.length)]};}
function jackknifeSigns(Y,X,vi){const fs=Math.sign(ols(Y,X).beta[vi+1]);let fl=0;for(let i=0;i<Y.length;i++){const f=ols([...Y.slice(0,i),...Y.slice(i+1)],[...X.slice(0,i),...X.slice(i+1)]);if(Math.sign(f.beta[vi+1])!==fs)fl++;}return fl;}

// ======================================================================
// Per-galaxy best-fit Upsilon_disk from Li+2018/2020
// Source: Li, Lelli, McGaugh, Schombert 2020, ApJS 247, 31
// Values: posterior median (DC14/ISO fits), 3.6 micron
// ======================================================================
const upsilonFit = {
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

// Quality flags
const pubUpsilon = new Set([
  'NGC0024','NGC0289','NGC0891','NGC1003','NGC1090','NGC1705','NGC2403',
  'NGC2683','NGC2841','NGC2903','NGC2915','NGC3198','NGC3521','NGC3726',
  'NGC3741','NGC3769','NGC3893','NGC4013','NGC4100','NGC4138','NGC4157',
  'NGC4217','NGC4559','NGC5005','NGC5033','NGC5055','NGC5371','NGC5907',
  'NGC6015','NGC6503','NGC6674','NGC7331','NGC7814'
]);
const estUpsilon = new Set([
  'UGC01281','UGC02953','UGC03205','UGC03546','UGC03580','UGC05721',
  'UGC06786','UGC06787','UGC06973','UGC08490','UGC08699','UGC09133','F571-8'
]);

// Build arrays
const gals = gals45;
const N = gals.length;
const Y = gals.map(g => g.logA0);
const logMhost = logMhost45;
const sdY = sd(Y);
function gapV(rms){return 100*(1-rms**2/sdY**2);}

const upsilon = gals.map(g => upsilonFit[g.name] || 0.50);
const logUpsilon = upsilon.map(u => Math.log10(u));
const deltaUpsilon = upsilon.map(u => u - 0.50);

// Only published Upsilon
const pubMask = gals.map(g => pubUpsilon.has(g.name));
const pubGals = gals.filter((g,i) => pubMask[i]);
const Ypub = pubGals.map(g => g.logA0);
const upsPub = pubGals.map(g => upsilonFit[g.name]);
const logUpsPub = upsPub.map(u => Math.log10(u));
const logMhostPub = pubGals.map(g => tdMap[g.name].logMhost);
const Npub = pubGals.length;

console.log('======================================================================');
console.log('  PHASE 59: TRUE STELLAR MASS-TO-LIGHT RATIO');
console.log('  Variable: Upsilon_disk (best-fit from Li+2020 MCMC)');
console.log('  N = ' + N + ' (full N=45 B\x27 sample)');
console.log('  N_published_Upsilon = ' + Npub);
console.log('  N_estimated_Upsilon = ' + (N - Npub));
console.log('======================================================================\n');

// ═══════════════════════════════════════════════
// RAW CORRELATIONS
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  RAW CORRELATIONS (full N=' + N + ')');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const rawU = pearson(logUpsilon, Y);
const spU = spearmanR(logUpsilon, Y);
console.log('  r(logUpsilon, logA0) = ' + rawU.r.toFixed(4) + ' (t=' + rawU.t.toFixed(3) + ')');
console.log('  Spearman = ' + spU.toFixed(4));
console.log('  Direction: ' + (rawU.r > 0 ? 'higher Upsilon -> higher a0' : 'higher Upsilon -> lower a0'));
console.log();

// Published-only
const rawPub = pearson(logUpsPub, Ypub);
console.log('  Published-only (N=' + Npub + '):');
console.log('  r(logUpsilon, logA0) = ' + rawPub.r.toFixed(4) + ' (t=' + rawPub.t.toFixed(3) + ')');
console.log();

// ═══════════════════════════════════════════════
// AFTER CONTROLS
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  AFTER CONTROLS (full N=' + N + ')');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const fMHI = ols(Y, gals.map(g => [g.logMHI]));
const aLogMHI = pearson(logUpsilon, fMHI.resid);
console.log('  After logMHI only: r=' + aLogMHI.r.toFixed(4) + ' (t=' + aLogMHI.t.toFixed(3) + ')');

const fEnv = ols(Y, gals.map((g,i) => [logMhost[i]]));
const aEnv = pearson(logUpsilon, fEnv.resid);
console.log('  After logMhost only: r=' + aEnv.r.toFixed(4) + ' (t=' + aEnv.t.toFixed(3) + ')');

const fME = ols(Y, gals.map((g,i) => [g.logMHI, logMhost[i]]));
const aME = pearson(logUpsilon, fME.resid);
console.log('  After logMHI + logMhost: r=' + aME.r.toFixed(4) + ' (t=' + aME.t.toFixed(3) + ')');
console.log();

// After Baseline A'
const X_Ap = gals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0]);
const fAp = ols(Y, X_Ap);
const aAp = pearson(logUpsilon, fAp.resid);
console.log('  After Baseline A\x27: r=' + aAp.r.toFixed(4) + ' (t=' + aAp.t.toFixed(3) + ')');

// After Baseline B'
const X_Bp = gals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0, g.logMeanRun]);
const fBp = ols(Y, X_Bp);
const aBp = pearson(logUpsilon, fBp.resid);
console.log('  After Baseline B\x27: r=' + aBp.r.toFixed(4) + ' (t=' + aBp.t.toFixed(3) + ')');
console.log();

// ═══════════════════════════════════════════════
// REGRESSION: B' + logUpsilon
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  B\x27 + logUpsilon REGRESSION');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const X_BpU = gals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0, g.logMeanRun, logUpsilon[i]]);
const fBpU = ols(Y, X_BpU);
console.log('  B\x27 + logUpsilon:');
console.log('    R2adj = ' + fBpU.r2adj.toFixed(4) + ' (vs B\x27: ' + fBp.r2adj.toFixed(4) + ')');
console.log('    logUpsilon coeff = ' + fBpU.beta[6].toFixed(4) + ' (t=' + fBpU.tStats[6].toFixed(3) + ')');
console.log();

// LOO
const looBp = looCV(Y, X_Bp);
const looBpU = looCV(Y, X_BpU);
console.log('  LOO (variance-based gap%):');
console.log('    B\x27:            ' + gapV(looBp).toFixed(1) + '%');
console.log('    B\x27 + logUpsilon: ' + gapV(looBpU).toFixed(1) + '% (delta=' + (gapV(looBpU)-gapV(looBp)).toFixed(1) + 'pp)');
console.log();

// Also test A' + logUpsilon
const X_ApU = gals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0, logUpsilon[i]]);
const fApU = ols(Y, X_ApU);
const looAp = looCV(Y, X_Ap);
const looApU = looCV(Y, X_ApU);
console.log('  A\x27 + logUpsilon:');
console.log('    R2adj = ' + fApU.r2adj.toFixed(4) + ' (vs A\x27: ' + fAp.r2adj.toFixed(4) + ')');
console.log('    logUpsilon coeff = ' + fApU.beta[5].toFixed(4) + ' (t=' + fApU.tStats[5].toFixed(3) + ')');
console.log('    LOO A\x27:          ' + gapV(looAp).toFixed(1) + '%');
console.log('    LOO A\x27+logUpsilon: ' + gapV(looApU).toFixed(1) + '% (delta=' + (gapV(looApU)-gapV(looAp)).toFixed(1) + 'pp)');
console.log();

// ═══════════════════════════════════════════════
// PERMUTATION / BOOTSTRAP / JACKKNIFE
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PERMUTATION / BOOTSTRAP / JACKKNIFE');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const pAp = permTest(Y, X_ApU, 4, 5000);
const pBp = permTest(Y, X_BpU, 5, 5000);
console.log('  Permutation (5000):');
console.log('    p(A\x27 + logUpsilon) = ' + pAp.toFixed(4));
console.log('    p(B\x27 + logUpsilon) = ' + pBp.toFixed(4));

const bAp = bootstrapCI(Y, X_ApU, 4, 2000);
const bBp = bootstrapCI(Y, X_BpU, 5, 2000);
console.log('  Bootstrap 95% CI:');
console.log('    A\x27: [' + bAp.lo.toFixed(4) + ', ' + bAp.hi.toFixed(4) + '] ' + (bAp.lo*bAp.hi>0?'EXCL':'incl') + ' zero');
console.log('    B\x27: [' + bBp.lo.toFixed(4) + ', ' + bBp.hi.toFixed(4) + '] ' + (bBp.lo*bBp.hi>0?'EXCL':'incl') + ' zero');

const jAp = jackknifeSigns(Y, X_ApU, 4);
const jBp = jackknifeSigns(Y, X_BpU, 5);
console.log('  Jackknife sign flips: A\x27=' + jAp + '/' + N + ', B\x27=' + jBp + '/' + N);
console.log();

// ═══════════════════════════════════════════════
// CROSS-CORRELATIONS WITH EXISTING PREDICTORS
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  CROSS-CORRELATIONS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const pairs = [
  {name:'logMHI', fn:g=>g.logMHI},
  {name:'logMhost', fn:(g,i)=>logMhost[i]},
  {name:'rcWiggliness', fn:g=>g.rcWiggliness},
  {name:'logSigma0', fn:g=>g.logSigma0},
  {name:'logMeanRun', fn:g=>g.logMeanRun},
  {name:'morphT', fn:g=>sparcMap[g.name]?.T ?? 5}
];
pairs.forEach(p => {
  const v = gals.map((g,i) => p.fn(g,i));
  const r = pearson(logUpsilon, v);
  console.log('  r(logUpsilon, ' + p.name.padEnd(14) + ') = ' + r.r.toFixed(4) + ' (t=' + r.t.toFixed(2) + ')');
});
console.log();

// ═══════════════════════════════════════════════
// PUBLISHED-ONLY SUBSET
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PUBLISHED-ONLY UPSILON (N=' + Npub + ')');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const X_Bp_pub = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhostPub[i], g.logSigma0, g.logMeanRun]);
const fBp_pub = ols(Ypub, X_Bp_pub);
const aBp_pub = pearson(logUpsPub, fBp_pub.resid);
console.log('  After B\x27: r=' + aBp_pub.r.toFixed(4) + ' (t=' + aBp_pub.t.toFixed(3) + ')');

const X_BpU_pub = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhostPub[i], g.logSigma0, g.logMeanRun, logUpsPub[i]]);
const fBpU_pub = ols(Ypub, X_BpU_pub);
console.log('  B\x27+logUpsilon on N=' + Npub + ':');
console.log('    logUpsilon t = ' + fBpU_pub.tStats[6].toFixed(3));
console.log('    R2adj = ' + fBpU_pub.r2adj.toFixed(4) + ' (vs B\x27: ' + fBp_pub.r2adj.toFixed(4) + ')');
const looBp_pub = looCV(Ypub, X_Bp_pub);
const looBpU_pub = looCV(Ypub, X_BpU_pub);
console.log('    LOO B\x27: ' + (100*(1-looBp_pub**2/sd(Ypub)**2)).toFixed(1) + '%');
console.log('    LOO B\x27+logUpsilon: ' + (100*(1-looBpU_pub**2/sd(Ypub)**2)).toFixed(1) + '% (delta=' + ((100*(1-looBpU_pub**2/sd(Ypub)**2))-(100*(1-looBp_pub**2/sd(Ypub)**2))).toFixed(1) + 'pp)');
console.log();

// ═══════════════════════════════════════════════
// GALAXY TABLE
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  GALAXY TABLE: Upsilon vs delta_a0');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

gals.sort((a,b) => upsilonFit[a.name] - upsilonFit[b.name]).forEach(g => {
  const u = upsilonFit[g.name];
  const q = pubUpsilon.has(g.name) ? 'pub' : 'est';
  console.log('  ' + g.name.padEnd(15) + ' Ups=' + u.toFixed(2) + ' da0=' + g.delta_a0.toFixed(3).padStart(7) + ' ' + q);
});
console.log();

// ═══════════════════════════════════════════════
// VERDICT
// ═══════════════════════════════════════════════
console.log('======================================================================');
console.log('  VERDICT');
console.log('======================================================================\n');

const addsAp = Math.abs(fApU.tStats[5]) > 1.65;
const addsBp = Math.abs(fBpU.tStats[6]) > 1.65;
const looOKAp = gapV(looApU) > gapV(looAp);
const looOKBp = gapV(looBpU) > gapV(looBp);

let verdict;
if (addsBp && looOKBp && pBp < 0.05) verdict = 'CONFIRMED';
else if (addsAp && looOKAp && pAp < 0.05) verdict = 'PARTIAL (adds above A\x27)';
else if (addsBp || addsAp) verdict = 'PARTIAL (weak signal)';
else verdict = 'FAIL';

console.log('  Adds above A\x27?  ' + (addsAp ? 'YES' : 'NO') + ' (t=' + fApU.tStats[5].toFixed(2) + ')');
console.log('  Adds above B\x27?  ' + (addsBp ? 'YES' : 'NO') + ' (t=' + fBpU.tStats[6].toFixed(2) + ')');
console.log('  LOO improves A\x27? ' + (looOKAp ? 'YES' : 'NO') + ' (delta=' + (gapV(looApU)-gapV(looAp)).toFixed(1) + 'pp)');
console.log('  LOO improves B\x27? ' + (looOKBp ? 'YES' : 'NO') + ' (delta=' + (gapV(looBpU)-gapV(looBp)).toFixed(1) + 'pp)');
console.log('  Perm p (A\x27):    ' + pAp.toFixed(4));
console.log('  Perm p (B\x27):    ' + pBp.toFixed(4));
console.log();
console.log('  ═══════════════════════════════════════════════');
console.log('  VERDICT: ' + verdict);
console.log('  ═══════════════════════════════════════════════');

const output = {
  phase: '59',
  variable: 'logUpsilon_disk',
  definition: 'log10 of best-fit stellar mass-to-light ratio at 3.6um from Li+2020 MCMC',
  source: 'Li, Lelli, McGaugh, Schombert 2020, ApJS 247, 31',
  n: N, nPubUpsilon: Npub,
  raw: rawU, spearman: spU,
  afterLogMHI: aLogMHI, afterLogMhost: aEnv, afterLogMHI_logMhost: aME,
  afterAp: { r: aAp.r, t: aAp.t, coeff: fApU.beta[5], tStat: fApU.tStats[5] },
  afterBp: { r: aBp.r, t: aBp.t, coeff: fBpU.beta[6], tStat: fBpU.tStats[6] },
  loo: { Ap: gapV(looAp), ApU: gapV(looApU), Bp: gapV(looBp), BpU: gapV(looBpU) },
  permutation: { Ap: pAp, Bp: pBp },
  bootstrap: { Ap: bAp, Bp: bBp },
  jackknife: { Ap: jAp, Bp: jBp },
  crossCorrelations: pairs.map(p => ({ name: p.name, r: pearson(logUpsilon, gals.map((g,i)=>p.fn(g,i))).r })),
  verdict,
  galaxyData: gals.map((g,i) => ({
    name: g.name, upsilon: upsilon[i], logUpsilon: logUpsilon[i],
    quality: pubUpsilon.has(g.name) ? 'published' : 'estimated'
  }))
};
fs.writeFileSync('public/phase59-stellar-ml.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase59-stellar-ml.json');
