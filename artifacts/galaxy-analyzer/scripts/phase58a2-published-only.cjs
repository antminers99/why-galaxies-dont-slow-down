const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const result = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8'));
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const allGals = stageA.galaxies;
const pubNames = new Set(result.tidalData.filter(t => t.quality === 'published').map(t => t.name));
const gals = allGals.filter(g => pubNames.has(g.name));
const N = gals.length;
const logMhost = gals.map(g => result.tidalData.find(t => t.name === g.name).logMhost);

function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function pearson(x,y){const n=x.length,mx=mean(x),my=mean(y);let sxy=0,sxx=0,syy=0;for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;syy+=(y[i]-my)**2;}const r=sxy/Math.sqrt(sxx*syy);return{r,t:r*Math.sqrt((n-2)/(1-r*r)),n};}
function spearmanR(x,y){const n=x.length;function rank(a){const s=a.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);const r=new Array(n);for(let i=0;i<n;){let j=i;while(j<n&&s[j].v===s[i].v)j++;const avg=(i+j-1)/2+1;for(let k=i;k<j;k++)r[s[k].i]=avg;i=j;}return r;}const rx=rank(x),ry=rank(y),mx=mean(rx),my=mean(ry);let sxy=0,sxx=0,syy=0;for(let i=0;i<n;i++){sxy+=(rx[i]-mx)*(ry[i]-my);sxx+=(rx[i]-mx)**2;syy+=(ry[i]-my)**2;}return sxy/Math.sqrt(sxx*syy);}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const se=Math.sqrt(rss/(n-p));const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se,resid,rss,tss,n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function permTest(Y,X,idx,nP){const Xb=X.map(r=>r.filter((_,j)=>j!==idx));const fb=ols(Y,Xb),ff=ols(Y,X);const obs=fb.rss-ff.rss;let c=0;for(let p=0;p<nP;p++){const nv=X.map(r=>r[idx]);for(let i=nv.length-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[nv[i],nv[j]]=[nv[j],nv[i]];}const Xp=X.map((r,k)=>r.map((v,j)=>j===idx?nv[k]:v));if(fb.rss-ols(Y,Xp).rss>=obs)c++;}return c/nP;}
function bootstrapCI(Y,X,vi,nB){const cs=[];for(let b=0;b<nB;b++){const idx=Array.from({length:Y.length},()=>Math.floor(Math.random()*Y.length));try{cs.push(ols(idx.map(i=>Y[i]),idx.map(i=>X[i])).beta[vi+1]);}catch(e){}}cs.sort((a,b)=>a-b);return{lo:cs[Math.floor(0.025*cs.length)],hi:cs[Math.floor(0.975*cs.length)],n:cs.length};}
function jackknifeSigns(Y,X,vi){const fs=Math.sign(ols(Y,X).beta[vi+1]);let fl=0;for(let i=0;i<Y.length;i++){const f=ols([...Y.slice(0,i),...Y.slice(i+1)],[...X.slice(0,i),...X.slice(i+1)]);if(Math.sign(f.beta[vi+1])!==fs)fl++;}return fl;}
const sdY = sd(gals.map(g=>g.logA0));
function gapV(rms){return 100*(1-rms**2/sdY**2);}

const Y = gals.map(g => g.logA0);

console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  PHASE 58a2: PUBLISHED-ONLY logM_host (N=' + N + ')                              ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

// Raw
const raw = pearson(logMhost, Y);
console.log('  Raw:');
console.log('    r = ' + raw.r.toFixed(4));
console.log('    t = ' + raw.t.toFixed(3));
console.log('    Spearman = ' + spearmanR(logMhost, Y).toFixed(4));
console.log();

// After envCode
const fEnv = ols(Y, gals.map(g=>[g.envCode]));
const aEnv = pearson(logMhost, fEnv.resid);
console.log('  After envCode only:');
console.log('    r = ' + aEnv.r.toFixed(4) + ', t = ' + aEnv.t.toFixed(3));

// After logMHI
const fMHI = ols(Y, gals.map(g=>[g.logMHI]));
const aMHI = pearson(logMhost, fMHI.resid);
console.log('  After logMHI only:');
console.log('    r = ' + aMHI.r.toFixed(4) + ', t = ' + aMHI.t.toFixed(3));

// After logMHI + envCode
const fME = ols(Y, gals.map(g=>[g.logMHI, g.envCode]));
const aME = pearson(logMhost, fME.resid);
console.log('  After logMHI + envCode:');
console.log('    r = ' + aME.r.toFixed(4) + ', t = ' + aME.t.toFixed(3));
console.log();

// After Baseline A
const X_A = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0]);
const fA = ols(Y, X_A);
const aA = pearson(logMhost, fA.resid);
const X_A_td = gals.map((g,i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, logMhost[i]]);
const fA_td = ols(Y, X_A_td);

console.log('  After Baseline A:');
console.log('    r = ' + aA.r.toFixed(4) + ', t = ' + aA.t.toFixed(3));
console.log('    A + logMhost: R2adj=' + fA_td.r2adj.toFixed(4) + ' (vs A: ' + fA.r2adj.toFixed(4) + ')');
console.log('    coeff(logMhost)=' + fA_td.beta[5].toFixed(4) + ' t=' + fA_td.tStats[5].toFixed(3));
console.log();

// After Baseline B  
const X_B = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun]);
const fB = ols(Y, X_B);
const aB = pearson(logMhost, fB.resid);
const X_B_td = gals.map((g,i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun, logMhost[i]]);
const fB_td = ols(Y, X_B_td);

console.log('  After Baseline B:');
console.log('    r = ' + aB.r.toFixed(4) + ', t = ' + aB.t.toFixed(3));
console.log('    B + logMhost: R2adj=' + fB_td.r2adj.toFixed(4) + ' (vs B: ' + fB.r2adj.toFixed(4) + ')');
console.log('    coeff(logMhost)=' + fB_td.beta[6].toFixed(4) + ' t=' + fB_td.tStats[6].toFixed(3));
console.log();

// LOO
const looA = looCV(Y, X_A);
const looAtd = looCV(Y, X_A_td);
const looB = looCV(Y, X_B);
const looBtd = looCV(Y, X_B_td);
console.log('  LOO (variance-based gap%):');
console.log('    A:             ' + gapV(looA).toFixed(1) + '%');
console.log('    A + logMhost:  ' + gapV(looAtd).toFixed(1) + '% (delta=' + (gapV(looAtd)-gapV(looA)).toFixed(1) + 'pp)');
console.log('    B:             ' + gapV(looB).toFixed(1) + '%');
console.log('    B + logMhost:  ' + gapV(looBtd).toFixed(1) + '% (delta=' + (gapV(looBtd)-gapV(looB)).toFixed(1) + 'pp)');
console.log();

// Perm
console.log('  Permutation (5000):');
const pA = permTest(Y, X_A_td, 4, 5000);
const pB = permTest(Y, X_B_td, 5, 5000);
console.log('    p(A + logMhost) = ' + pA.toFixed(4));
console.log('    p(B + logMhost) = ' + pB.toFixed(4));
console.log();

// Bootstrap
const bA = bootstrapCI(Y, X_A_td, 4, 2000);
const bB = bootstrapCI(Y, X_B_td, 5, 2000);
console.log('  Bootstrap 95% CI:');
console.log('    A: [' + bA.lo.toFixed(4) + ', ' + bA.hi.toFixed(4) + '] ' + (bA.lo*bA.hi>0?'EXCLUDES zero':'includes zero'));
console.log('    B: [' + bB.lo.toFixed(4) + ', ' + bB.hi.toFixed(4) + '] ' + (bB.lo*bB.hi>0?'EXCLUDES zero':'includes zero'));
console.log();

// Jackknife
const jA = jackknifeSigns(Y, X_A_td, 4);
const jB = jackknifeSigns(Y, X_B_td, 5);
console.log('  Jackknife sign flips: A=' + jA + '/' + N + ', B=' + jB + '/' + N);
console.log();

// envCode after adding logMhost
console.log('  envCode after adding logMhost:');
console.log('    t(envCode) in A+logMhost: ' + fA_td.tStats[3].toFixed(3));
console.log();

// Correlation with envCode
const corrEnvMhost = pearson(gals.map(g=>g.envCode), logMhost);
console.log('  r(envCode, logMhost) = ' + corrEnvMhost.r.toFixed(4));
console.log('  r(logMHI, logMhost) = ' + pearson(gals.map(g=>g.logMHI), logMhost).r.toFixed(4));
console.log();

// ======================================================================
// CRITICAL CHECK: Does logMhost add above envCode because it's just
// a better envCode, or because it captures something new?
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  DIAGNOSTIC: Is logMhost just a better-coded envCode?');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Replace envCode with logMhost in A
const X_A_replace = gals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0]);
const fA_rep = ols(Y, X_A_replace);
const looA_rep = looCV(Y, X_A_replace);
console.log('  A (with envCode):    R2adj=' + fA.r2adj.toFixed(4) + ', LOO gap=' + gapV(looA).toFixed(1) + '%');
console.log('  A (envCode→logMhost): R2adj=' + fA_rep.r2adj.toFixed(4) + ', LOO gap=' + gapV(looA_rep).toFixed(1) + '%');
console.log('  Difference: R2adj ' + (fA_rep.r2adj-fA.r2adj>0?'+':'') + (fA_rep.r2adj-fA.r2adj).toFixed(4) + ', LOO ' + (gapV(looA_rep)-gapV(looA)>0?'+':'') + (gapV(looA_rep)-gapV(looA)).toFixed(1) + 'pp');
console.log();

// Replace envCode with logMhost in B
const X_B_replace = gals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0, g.logMeanRun]);
const fB_rep = ols(Y, X_B_replace);
const looB_rep = looCV(Y, X_B_replace);
console.log('  B (with envCode):    R2adj=' + fB.r2adj.toFixed(4) + ', LOO gap=' + gapV(looB).toFixed(1) + '%');
console.log('  B (envCode→logMhost): R2adj=' + fB_rep.r2adj.toFixed(4) + ', LOO gap=' + gapV(looB_rep).toFixed(1) + '%');
console.log('  Difference: R2adj ' + (fB_rep.r2adj-fB.r2adj>0?'+':'') + (fB_rep.r2adj-fB.r2adj).toFixed(4) + ', LOO ' + (gapV(looB_rep)-gapV(looB)>0?'+':'') + (gapV(looB_rep)-gapV(looB)).toFixed(1) + 'pp');
console.log();

console.log('  INTERPRETATION:');
if (gapV(looA_rep) > gapV(looA) + 1) {
  console.log('  logMhost is a BETTER environmental predictor than envCode');
  console.log('  → continuous tidal depth > binary classification');
} else if (gapV(looA_rep) > gapV(looA) - 1) {
  console.log('  logMhost is COMPARABLE to envCode');
  console.log('  → continuous version carries same info, no meaningful gain');
} else {
  console.log('  logMhost is WORSE than envCode');
  console.log('  → binary flag is cleaner for this sample');
}
" 2>&1
