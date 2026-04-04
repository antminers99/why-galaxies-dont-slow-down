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
function pearson(x,y){const n=x.length;const mx=mean(x),my=mean(y);let sxy=0,sxx=0,syy=0;for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;syy+=(y[i]-my)**2;}const r=sxy/Math.sqrt(sxx*syy);return{r,t:r*Math.sqrt((n-2)/(1-r*r)),n};}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se:Math.sqrt(rss/(n-p)),resid,rss,tss,n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function permTest(Y,X,idx,nP){const Xb=X.map(r=>r.filter((_,j)=>j!==idx));const fb=ols(Y,Xb);const obs=fb.rss-ols(Y,X).rss;let c=0;for(let p=0;p<nP;p++){const nv=X.map(r=>r[idx]);for(let i=nv.length-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[nv[i],nv[j]]=[nv[j],nv[i]];}const Xp=X.map((r,k)=>r.map((v,j)=>j===idx?nv[k]:v));if(fb.rss-ols(Y,Xp).rss>=obs)c++;}return c/nP;}
function bootstrapCI(Y,X,vi,nB){const cs=[];for(let b=0;b<nB;b++){const idx=Array.from({length:Y.length},()=>Math.floor(Math.random()*Y.length));try{cs.push(ols(idx.map(i=>Y[i]),idx.map(i=>X[i])).beta[vi+1]);}catch(e){}}cs.sort((a,b)=>a-b);return{lo:cs[Math.floor(0.025*cs.length)],hi:cs[Math.floor(0.975*cs.length)],median:cs[Math.floor(0.5*cs.length)]};}
function jackknifeSigns(Y,X,vi){const fs=Math.sign(ols(Y,X).beta[vi+1]);let fl=0;for(let i=0;i<Y.length;i++){const f=ols([...Y.slice(0,i),...Y.slice(i+1)],[...X.slice(0,i),...X.slice(i+1)]);if(Math.sign(f.beta[vi+1])!==fs)fl++;}return fl;}

// ======================================================================
// Upsilon source data
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

const N = gals45.length;
const Y = gals45.map(g => g.logA0);
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const logUps = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const morphT = gals45.map(g => sparcMap[g.name]?.T || 5);
const sdY = sd(Y);
function gapV(rms){return 100*(1-rms**2/sdY**2);}

// ======================================================================
// STEP 1: Compute Υ★⊥ = residualized logUpsilon
// ======================================================================
console.log('======================================================================');
console.log('  PHASE 59b: FREEZE B" WITH RESIDUALIZED Υ★');
console.log('  N = ' + N);
console.log('======================================================================\n');

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 1: Compute Υ★⊥');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const Xconf = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT[i]]);
const fConf = ols(logUps, Xconf);
const upsPerp = fConf.resid;
const r2conf = 1 - fConf.rss / fConf.tss;

console.log('  Υ★⊥ = residual of logΥ after removing (logMHI, logΣ₀, morphT)');
console.log('  Confounders: logMHI coeff=' + fConf.beta[1].toFixed(4) + ', logΣ₀ coeff=' + fConf.beta[2].toFixed(4) + ', morphT coeff=' + fConf.beta[3].toFixed(4));
console.log('  R² of confounders on logΥ: ' + r2conf.toFixed(4) + ' (' + (r2conf*100).toFixed(1) + '%)');
console.log('  SD(Υ★⊥) = ' + sd(upsPerp).toFixed(4) + ' dex (original logΥ: ' + sd(logUps).toFixed(4) + ')');
console.log('  By construction: r(Υ★⊥, logMHI) = ' + pearson(upsPerp, gals45.map(g=>g.logMHI)).r.toFixed(6));
console.log('  By construction: r(Υ★⊥, logΣ₀) = ' + pearson(upsPerp, gals45.map(g=>g.logSigma0)).r.toFixed(6));
console.log('  By construction: r(Υ★⊥, morphT) = ' + pearson(upsPerp, morphT).r.toFixed(6));
console.log();

// Per-galaxy values
console.log('  Per-galaxy Υ★⊥:');
gals45.forEach((g,i) => {
  console.log('    ' + g.name.padEnd(15) + ' logΥ=' + logUps[i].toFixed(4) + ' Υ★⊥=' + upsPerp[i].toFixed(4));
});
console.log();

// ======================================================================
// STEP 2: Fit B' and B" on N=45
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 2: B\x27 vs B" regression on N=' + N);
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const X_Bp = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0, g.logMeanRun]);
const X_Bdp = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0, g.logMeanRun, upsPerp[i]]);

const fBp = ols(Y, X_Bp);
const fBdp = ols(Y, X_Bdp);

const looBp = looCV(Y, X_Bp);
const looBdp = looCV(Y, X_Bdp);

const varNames_Bp = ['intercept','logMHI','rcWiggliness','logMhost','logSigma0','logMeanRun'];
const varNames_Bdp = [...varNames_Bp, 'Υ★⊥'];

console.log('  ╔═══════════════════════════════════════════════════════════╗');
console.log('  ║  BASELINE B\x27 (FROZEN, reference)                       ║');
console.log('  ╚═══════════════════════════════════════════════════════════╝');
console.log();
varNames_Bp.forEach((v,j) => {
  const ci_lo = fBp.beta[j] - 1.96*fBp.seBeta[j];
  const ci_hi = fBp.beta[j] + 1.96*fBp.seBeta[j];
  console.log('    ' + v.padEnd(16) + ' β=' + fBp.beta[j].toFixed(4).padStart(8) + ' t=' + fBp.tStats[j].toFixed(3).padStart(7) + '  95% CI [' + ci_lo.toFixed(3) + ', ' + ci_hi.toFixed(3) + ']');
});
console.log('    R²adj = ' + fBp.r2adj.toFixed(4));
console.log('    LOO gap% = ' + gapV(looBp).toFixed(1) + '%');
console.log('    τ_resid = ' + looBp.toFixed(4) + ' dex');
console.log();

console.log('  ╔═══════════════════════════════════════════════════════════╗');
console.log('  ║  BASELINE B" (CANDIDATE)                                ║');
console.log('  ║  B\x27 + Υ★⊥                                              ║');
console.log('  ╚═══════════════════════════════════════════════════════════╝');
console.log();
varNames_Bdp.forEach((v,j) => {
  const ci_lo = fBdp.beta[j] - 1.96*fBdp.seBeta[j];
  const ci_hi = fBdp.beta[j] + 1.96*fBdp.seBeta[j];
  const flag = j === 6 ? ' ◄ NEW' : '';
  console.log('    ' + v.padEnd(16) + ' β=' + fBdp.beta[j].toFixed(4).padStart(8) + ' t=' + fBdp.tStats[j].toFixed(3).padStart(7) + '  95% CI [' + ci_lo.toFixed(3) + ', ' + ci_hi.toFixed(3) + ']' + flag);
});
console.log('    R²adj = ' + fBdp.r2adj.toFixed(4));
console.log('    LOO gap% = ' + gapV(looBdp).toFixed(1) + '%');
console.log('    τ_resid = ' + looBdp.toFixed(4) + ' dex');
console.log();

console.log('  COMPARISON:');
console.log('    R²adj: ' + fBp.r2adj.toFixed(4) + ' → ' + fBdp.r2adj.toFixed(4) + ' (Δ=' + (fBdp.r2adj-fBp.r2adj).toFixed(4) + ')');
console.log('    LOO:   ' + gapV(looBp).toFixed(1) + '% → ' + gapV(looBdp).toFixed(1) + '% (Δ=' + (gapV(looBdp)-gapV(looBp)).toFixed(1) + 'pp)');
console.log('    τ:     ' + looBp.toFixed(4) + ' → ' + looBdp.toFixed(4) + ' dex');
console.log('    Unexplained: ' + (100-gapV(looBdp)).toFixed(1) + '% (was ' + (100-gapV(looBp)).toFixed(1) + '%)');
console.log();

// ======================================================================
// STEP 3: Robustness — Permutation, Bootstrap, Jackknife
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 3: Robustness tests for Υ★⊥ in B"');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const permP = permTest(Y, X_Bdp, 5, 10000);
const bootCI = bootstrapCI(Y, X_Bdp, 5, 5000);
const jkFlips = jackknifeSigns(Y, X_Bdp, 5);

console.log('  Permutation test (10000 iterations):');
console.log('    p = ' + permP.toFixed(4) + (permP < 0.05 ? ' ✅ significant' : permP < 0.10 ? ' ⚠️ marginal' : ' ❌ not significant'));
console.log();
console.log('  Bootstrap 95% CI (5000 iterations):');
console.log('    Υ★⊥ coeff: [' + bootCI.lo.toFixed(4) + ', ' + bootCI.hi.toFixed(4) + ']');
console.log('    Median: ' + bootCI.median.toFixed(4));
console.log('    ' + (bootCI.lo*bootCI.hi > 0 ? '✅ Excludes zero' : '❌ Includes zero'));
console.log();
console.log('  Jackknife sign flips: ' + jkFlips + '/' + N + ' (' + (jkFlips === 0 ? '✅ perfectly stable' : jkFlips <= 2 ? '✅ stable' : '⚠️ ' + jkFlips + ' flips'));
console.log();

// Bootstrap all coefficients
console.log('  Bootstrap stability of ALL B" coefficients:');
for (let vi = 0; vi < varNames_Bdp.length; vi++) {
  if (vi === 0) continue;
  const bc = bootstrapCI(Y, X_Bdp, vi-1, 2000);
  const exc = bc.lo*bc.hi > 0 ? 'EXCL' : 'incl';
  console.log('    ' + varNames_Bdp[vi].padEnd(16) + ' [' + bc.lo.toFixed(4) + ', ' + bc.hi.toFixed(4) + '] ' + exc + ' zero');
}
console.log();

// ======================================================================
// STEP 4: Comparison table B' vs B" vs alternative
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  STEP 4: Model comparison');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Also A' + Υ★⊥
const X_Ap = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0]);
const X_Adp = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0, upsPerp[i]]);
const fAp = ols(Y, X_Ap);
const fAdp = ols(Y, X_Adp);
const looAp = looCV(Y, X_Ap);
const looAdp = looCV(Y, X_Adp);

// N=33 published-only B" side analysis
const pubUpsSet = new Set(['NGC0024','NGC0289','NGC0891','NGC1003','NGC1090','NGC1705','NGC2403',
  'NGC2683','NGC2841','NGC2903','NGC2915','NGC3198','NGC3521','NGC3726',
  'NGC3741','NGC3769','NGC3893','NGC4013','NGC4100','NGC4138','NGC4157',
  'NGC4217','NGC4559','NGC5005','NGC5033','NGC5055','NGC5371','NGC5907',
  'NGC6015','NGC6503','NGC6674','NGC7331','NGC7814']);
const pubGals = gals45.filter(g => pubUpsSet.has(g.name));
const Y33 = pubGals.map(g => g.logA0);
const sdY33 = sd(Y33);
function gapV33(rms){return 100*(1-rms**2/sdY33**2);}

// Recompute Υ★⊥ for N=33 using SAME confounding model from N=45
const upsPerp33 = pubGals.map(g => {
  const i45 = gals45.indexOf(g);
  return upsPerp[i45];
});
const logMhost33 = pubGals.map(g => tdMap[g.name].logMhost);
const X_Bp33 = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logSigma0, g.logMeanRun]);
const X_Bdp33 = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logSigma0, g.logMeanRun, upsPerp33[i]]);
const looBp33 = looCV(Y33, X_Bp33);
const looBdp33 = looCV(Y33, X_Bdp33);
const fBdp33 = ols(Y33, X_Bdp33);

// BΣ→Υ alternative (replace logSigma0 with logUpsilon raw) on N=33
const logUps33 = pubGals.map(g => Math.log10(upsilonMap[g.name]));
const X_BsU33 = pubGals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost33[i], g.logMeanRun, logUps33[i]]);
const looBsU33 = looCV(Y33, X_BsU33);

console.log('  ┌─────────────────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Model               │ R²adj    │ LOO gap% │ τ_resid  │ k (pred) │');
console.log('  ├─────────────────────┼──────────┼──────────┼──────────┼──────────┤');
console.log('  │ A\x27 (N=45)           │ ' + fAp.r2adj.toFixed(4) + '   │ ' + gapV(looAp).toFixed(1).padStart(5) + '%   │ ' + looAp.toFixed(4) + '   │ 4        │');
console.log('  │ A\x27+Υ★⊥ (N=45)       │ ' + fAdp.r2adj.toFixed(4) + '   │ ' + gapV(looAdp).toFixed(1).padStart(5) + '%   │ ' + looAdp.toFixed(4) + '   │ 5        │');
console.log('  │ B\x27 (N=45)           │ ' + fBp.r2adj.toFixed(4) + '   │ ' + gapV(looBp).toFixed(1).padStart(5) + '%   │ ' + looBp.toFixed(4) + '   │ 5        │');
console.log('  │ B" (N=45) ◄         │ ' + fBdp.r2adj.toFixed(4) + '   │ ' + gapV(looBdp).toFixed(1).padStart(5) + '%   │ ' + looBdp.toFixed(4) + '   │ 6        │');
console.log('  ├─────────────────────┼──────────┼──────────┼──────────┼──────────┤');
console.log('  │ B\x27 (N=33 pub)       │    —     │ ' + gapV33(looBp33).toFixed(1).padStart(5) + '%   │ ' + looBp33.toFixed(4) + '   │ 5        │');
console.log('  │ B" (N=33 pub)       │ ' + fBdp33.r2adj.toFixed(4) + '   │ ' + gapV33(looBdp33).toFixed(1).padStart(5) + '%   │ ' + looBdp33.toFixed(4) + '   │ 6        │');
console.log('  │ BΣ→Υ (N=33, side)   │    —     │ ' + gapV33(looBsU33).toFixed(1).padStart(5) + '%   │ ' + looBsU33.toFixed(4) + '   │ 5        │');
console.log('  └─────────────────────┴──────────┴──────────┴──────────┴──────────┘');
console.log();

// ======================================================================
// VERDICT
// ======================================================================
console.log('======================================================================');
console.log('  VERDICT');
console.log('======================================================================\n');

const looPass = gapV(looBdp) > gapV(looBp);
const tPass = Math.abs(fBdp.tStats[6]) > 1.65;
const permPass = permP < 0.10;
const bootPass = bootCI.lo * bootCI.hi > 0;
const jkPass = jkFlips <= 2;

console.log('  LOO improves:     ' + (looPass ? 'YES' : 'NO') + ' (Δ=' + (gapV(looBdp)-gapV(looBp)).toFixed(1) + 'pp)');
console.log('  t-stat > 1.65:    ' + (tPass ? 'YES' : 'NO') + ' (t=' + fBdp.tStats[6].toFixed(3) + ')');
console.log('  Perm p < 0.10:    ' + (permPass ? 'YES' : 'NO') + ' (p=' + permP.toFixed(4) + ')');
console.log('  Bootstrap excl 0: ' + (bootPass ? 'YES' : 'NO') + ' [' + bootCI.lo.toFixed(4) + ', ' + bootCI.hi.toFixed(4) + ']');
console.log('  Jackknife stable: ' + (jkPass ? 'YES' : 'NO') + ' (' + jkFlips + '/' + N + ')');
console.log();

const passCount = [looPass, tPass, permPass, bootPass, jkPass].filter(Boolean).length;
let verdict;
if (passCount >= 4) verdict = 'FREEZE B" — CONFIRMED';
else if (passCount >= 3) verdict = 'FREEZE B" — PARTIAL (proceed with caution)';
else verdict = 'DO NOT FREEZE — insufficient robustness';

console.log('  Tests passed: ' + passCount + '/5');
console.log();
console.log('  ═══════════════════════════════════════════════════════════');
console.log('  ' + verdict);
console.log('  ═══════════════════════════════════════════════════════════');

// Save output
const output = {
  phase: '59b',
  title: 'Freeze B" with residualized Υ★',
  n: N,
  upsilonPerp: {
    definition: 'Υ★⊥ = residual of log(Υ_disk) after OLS on (logMHI, logΣ₀, morphT)',
    confoundModel: {
      intercept: fConf.beta[0], logMHI: fConf.beta[1],
      logSigma0: fConf.beta[2], morphT: fConf.beta[3],
      r2: r2conf
    },
    perGalaxy: gals45.map((g,i) => ({
      name: g.name,
      logUpsilon: logUps[i],
      upsilonPerp: upsPerp[i]
    }))
  },
  baselineBp: {
    coefficients: varNames_Bp.map((v,j) => ({name:v, beta:fBp.beta[j], se:fBp.seBeta[j], t:fBp.tStats[j]})),
    r2adj: fBp.r2adj, looGap: gapV(looBp), tauResid: looBp
  },
  baselineBdp: {
    coefficients: varNames_Bdp.map((v,j) => ({name:v, beta:fBdp.beta[j], se:fBdp.seBeta[j], t:fBdp.tStats[j]})),
    r2adj: fBdp.r2adj, looGap: gapV(looBdp), tauResid: looBdp
  },
  robustness: {
    permutation: { p: permP, n: 10000 },
    bootstrap: bootCI,
    jackknife: { flips: jkFlips, n: N }
  },
  comparison: {
    Ap: { loo: gapV(looAp), tau: looAp },
    Adp: { loo: gapV(looAdp), tau: looAdp },
    Bp: { loo: gapV(looBp), tau: looBp },
    Bdp: { loo: gapV(looBdp), tau: looBdp },
    Bp_N33: { loo: gapV33(looBp33), tau: looBp33 },
    Bdp_N33: { loo: gapV33(looBdp33), tau: looBdp33 },
    BsU_N33: { loo: gapV33(looBsU33), tau: looBsU33 }
  },
  verdict,
  passCount
};
fs.writeFileSync('public/phase59b-freeze-Bdp.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase59b-freeze-Bdp.json');
