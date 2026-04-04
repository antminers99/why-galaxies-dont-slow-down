const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const result = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8'));

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
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se:Math.sqrt(rss/(n-p)),resid,rss,tss,n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function permTest(Y,X,idx,nP){const Xb=X.map(r=>r.filter((_,j)=>j!==idx));const fb=ols(Y,Xb),ff=ols(Y,X);const obs=fb.rss-ff.rss;let c=0;for(let p=0;p<nP;p++){const nv=X.map(r=>r[idx]);for(let i=nv.length-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[nv[i],nv[j]]=[nv[j],nv[i]];}const Xp=X.map((r,k)=>r.map((v,j)=>j===idx?nv[k]:v));if(fb.rss-ols(Y,Xp).rss>=obs)c++;}return c/nP;}
function bootstrapCI(Y,X,vi,nB){const cs=[];for(let b=0;b<nB;b++){const idx=Array.from({length:Y.length},()=>Math.floor(Math.random()*Y.length));try{cs.push(ols(idx.map(i=>Y[i]),idx.map(i=>X[i])).beta[vi+1]);}catch(e){}}cs.sort((a,b)=>a-b);return{lo:cs[Math.floor(0.025*cs.length)],hi:cs[Math.floor(0.975*cs.length)],n:cs.length};}
function jackknifeSigns(Y,X,vi){const fs=Math.sign(ols(Y,X).beta[vi+1]);let fl=0;for(let i=0;i<Y.length;i++){const f=ols([...Y.slice(0,i),...Y.slice(i+1)],[...X.slice(0,i),...X.slice(i+1)]);if(Math.sign(f.beta[vi+1])!==fs)fl++;}return fl;}
const sdY=sd(gals.map(g=>g.logA0));
function gapV(rms){return 100*(1-rms**2/sdY**2);}

const Y = gals.map(g => g.logA0);

console.log('N = ' + N + ' (published-only)\n');

// Raw
const raw = pearson(logMhost, Y);
console.log('Raw: r=' + raw.r.toFixed(4) + ' t=' + raw.t.toFixed(3) + ' Spearman=' + spearmanR(logMhost, Y).toFixed(4));

// After envCode
const fEnv=ols(Y,gals.map(g=>[g.envCode]));
const aEnv=pearson(logMhost,fEnv.resid);
console.log('After envCode: r=' + aEnv.r.toFixed(4) + ' t=' + aEnv.t.toFixed(3));

// After logMHI
const fMHI=ols(Y,gals.map(g=>[g.logMHI]));
const aMHI=pearson(logMhost,fMHI.resid);
console.log('After logMHI: r=' + aMHI.r.toFixed(4) + ' t=' + aMHI.t.toFixed(3));

// After logMHI+envCode
const fME=ols(Y,gals.map(g=>[g.logMHI,g.envCode]));
const aME=pearson(logMhost,fME.resid);
console.log('After logMHI+envCode: r=' + aME.r.toFixed(4) + ' t=' + aME.t.toFixed(3));
console.log();

// Baseline A
const XA=gals.map(g=>[g.logMHI,g.rcWiggliness,g.envCode,g.logSigma0]);
const fA=ols(Y,XA);
const aA=pearson(logMhost,fA.resid);
const XAtd=gals.map((g,i)=>[g.logMHI,g.rcWiggliness,g.envCode,g.logSigma0,logMhost[i]]);
const fAtd=ols(Y,XAtd);
console.log('After A: r=' + aA.r.toFixed(4) + ' t=' + aA.t.toFixed(3));
console.log('A+td: R2adj=' + fAtd.r2adj.toFixed(4) + ' (A=' + fA.r2adj.toFixed(4) + ') coeff=' + fAtd.beta[5].toFixed(4) + ' t=' + fAtd.tStats[5].toFixed(3));

// Baseline B
const XB=gals.map(g=>[g.logMHI,g.rcWiggliness,g.envCode,g.logSigma0,g.logMeanRun]);
const fB=ols(Y,XB);
const aB=pearson(logMhost,fB.resid);
const XBtd=gals.map((g,i)=>[g.logMHI,g.rcWiggliness,g.envCode,g.logSigma0,g.logMeanRun,logMhost[i]]);
const fBtd=ols(Y,XBtd);
console.log('After B: r=' + aB.r.toFixed(4) + ' t=' + aB.t.toFixed(3));
console.log('B+td: R2adj=' + fBtd.r2adj.toFixed(4) + ' (B=' + fB.r2adj.toFixed(4) + ') coeff=' + fBtd.beta[6].toFixed(4) + ' t=' + fBtd.tStats[6].toFixed(3));
console.log();

// LOO
const looA=looCV(Y,XA),looAtd=looCV(Y,XAtd),looB=looCV(Y,XB),looBtd=looCV(Y,XBtd);
console.log('LOO gap%: A=' + gapV(looA).toFixed(1) + ' A+td=' + gapV(looAtd).toFixed(1) + ' (d=' + (gapV(looAtd)-gapV(looA)).toFixed(1) + ')');
console.log('LOO gap%: B=' + gapV(looB).toFixed(1) + ' B+td=' + gapV(looBtd).toFixed(1) + ' (d=' + (gapV(looBtd)-gapV(looB)).toFixed(1) + ')');
console.log();

// Perm
const pA=permTest(Y,XAtd,4,5000);
const pB=permTest(Y,XBtd,5,5000);
console.log('Perm: A=' + pA.toFixed(4) + ' B=' + pB.toFixed(4));

// Boot
const bA=bootstrapCI(Y,XAtd,4,2000);
const bB=bootstrapCI(Y,XBtd,5,2000);
console.log('Boot A: [' + bA.lo.toFixed(4) + ', ' + bA.hi.toFixed(4) + '] ' + (bA.lo*bA.hi>0?'EXCL':'incl') + ' zero');
console.log('Boot B: [' + bB.lo.toFixed(4) + ', ' + bB.hi.toFixed(4) + '] ' + (bB.lo*bB.hi>0?'EXCL':'incl') + ' zero');

// Jack
const jA=jackknifeSigns(Y,XAtd,4);
const jB=jackknifeSigns(Y,XBtd,5);
console.log('Jack flips: A=' + jA + '/' + N + ' B=' + jB + '/' + N);
console.log();

// envCode after
console.log('envCode t in A+td: ' + fAtd.tStats[3].toFixed(3));
console.log('r(envCode,logMhost): ' + pearson(gals.map(g=>g.envCode),logMhost).r.toFixed(4));
console.log('r(logMHI,logMhost): ' + pearson(gals.map(g=>g.logMHI),logMhost).r.toFixed(4));
console.log();

// REPLACEMENT TEST
console.log('=== REPLACEMENT: envCode -> logMhost in A ===');
const XArep=gals.map((g,i)=>[g.logMHI,g.rcWiggliness,logMhost[i],g.logSigma0]);
const fArep=ols(Y,XArep);
const looArep=looCV(Y,XArep);
console.log('A(envCode): R2adj=' + fA.r2adj.toFixed(4) + ' LOO=' + gapV(looA).toFixed(1) + '%');
console.log('A(logMhost): R2adj=' + fArep.r2adj.toFixed(4) + ' LOO=' + gapV(looArep).toFixed(1) + '%');
console.log('Delta R2adj: ' + (fArep.r2adj-fA.r2adj).toFixed(4));
console.log('Delta LOO: ' + (gapV(looArep)-gapV(looA)).toFixed(1) + 'pp');
console.log();
console.log('=== REPLACEMENT: envCode -> logMhost in B ===');
const XBrep=gals.map((g,i)=>[g.logMHI,g.rcWiggliness,logMhost[i],g.logSigma0,g.logMeanRun]);
const fBrep=ols(Y,XBrep);
const looBrep=looCV(Y,XBrep);
console.log('B(envCode): R2adj=' + fB.r2adj.toFixed(4) + ' LOO=' + gapV(looB).toFixed(1) + '%');
console.log('B(logMhost): R2adj=' + fBrep.r2adj.toFixed(4) + ' LOO=' + gapV(looBrep).toFixed(1) + '%');
console.log('Delta R2adj: ' + (fBrep.r2adj-fB.r2adj).toFixed(4));
console.log('Delta LOO: ' + (gapV(looBrep)-gapV(looB)).toFixed(1) + 'pp');
