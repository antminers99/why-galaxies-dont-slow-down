const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const tidalData = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const tdMap = {};
tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const gals = stageA.galaxies.filter(g => pubNames.has(g.name));
const N = gals.length;
const logMhost = gals.map(g => tdMap[g.name].logMhost);
const distMap = {};
sparc.forEach(s => { distMap[s.name] = s.D; });

function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function pearson(x,y){const n=x.length,mx=mean(x),my=mean(y);let sxy=0,sxx=0,syy=0;for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;syy+=(y[i]-my)**2;}const r=sxy/Math.sqrt(sxx*syy);return{r,t:r*Math.sqrt((n-2)/(1-r*r)),n};}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se:Math.sqrt(rss/(n-p)),resid,rss,tss,n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function bootstrapCI(Y,X,nB){const p=X[0].length;const results=[];for(let b=0;b<nB;b++){const idx=Array.from({length:Y.length},()=>Math.floor(Math.random()*Y.length));try{results.push(ols(idx.map(i=>Y[i]),idx.map(i=>X[i])).beta);}catch(e){}}const cis=[];for(let j=0;j<=p;j++){const vals=results.map(r=>r[j]).sort((a,b)=>a-b);cis.push({lo:vals[Math.floor(0.025*vals.length)],hi:vals[Math.floor(0.975*vals.length)]});}return cis;}

const sdY = sd(gals.map(g=>g.logA0));
function gapV(rms){return 100*(1-rms**2/sdY**2);}

const Y = gals.map(g => g.logA0);
const varNames_Ap = ['logMHI','rcWiggliness','logMhost','logSigma0'];
const varNames_Bp = ['logMHI','rcWiggliness','logMhost','logSigma0','logMeanRun'];

const X_Ap = gals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0]);
const X_Bp = gals.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0, g.logMeanRun]);

const fitAp = ols(Y, X_Ap);
const fitBp = ols(Y, X_Bp);
const looAp = looCV(Y, X_Ap);
const looBp = looCV(Y, X_Bp);

// Also compute OLD baselines on same N=45 for comparison
const X_A_old = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0]);
const X_B_old = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun]);
const fitA_old = ols(Y, X_A_old);
const fitB_old = ols(Y, X_B_old);
const looA_old = looCV(Y, X_A_old);
const looB_old = looCV(Y, X_B_old);

console.log('======================================================================');
console.log('  PHASE 58a4: FREEZE UPDATED BASELINES');
console.log('  N = ' + N + ' (published logMhost subsample)');
console.log('  sdY = ' + sdY.toFixed(4) + ' dex');
console.log('======================================================================\n');

// ═══════════════════════════════════════════════
// BASELINE A' (conservative)
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  BASELINE A\' (conservative)');
console.log('  log(a0) = b0 + b1*logMHI + b2*wig + b3*logMhost + b4*logSigma0');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const bootAp = bootstrapCI(Y, X_Ap, 2000);
console.log('  Coefficients:');
console.log('    Intercept:   ' + fitAp.beta[0].toFixed(4) + ' (t=' + fitAp.tStats[0].toFixed(2) + ') CI=[' + bootAp[0].lo.toFixed(3) + ',' + bootAp[0].hi.toFixed(3) + ']');
for(let j=0;j<varNames_Ap.length;j++){
  console.log('    ' + varNames_Ap[j].padEnd(12) + ' ' + fitAp.beta[j+1].toFixed(4) + ' (t=' + fitAp.tStats[j+1].toFixed(2) + ') CI=[' + bootAp[j+1].lo.toFixed(3) + ',' + bootAp[j+1].hi.toFixed(3) + ']');
}
console.log();
console.log('  R2adj = ' + fitAp.r2adj.toFixed(4));
console.log('  SE = ' + fitAp.se.toFixed(4) + ' dex');
console.log('  LOO-RMS = ' + looAp.toFixed(4) + ' dex');
console.log('  LOO gap% = ' + gapV(looAp).toFixed(1) + '%');
console.log('  Residual tau = ' + fitAp.se.toFixed(3) + ' dex');
console.log();

console.log('  vs OLD Baseline A (envCode, same N=45):');
console.log('    R2adj: ' + fitA_old.r2adj.toFixed(4) + ' -> ' + fitAp.r2adj.toFixed(4) + ' (delta=' + (fitAp.r2adj-fitA_old.r2adj).toFixed(4) + ')');
console.log('    LOO:   ' + gapV(looA_old).toFixed(1) + '% -> ' + gapV(looAp).toFixed(1) + '% (delta=' + (gapV(looAp)-gapV(looA_old)).toFixed(1) + 'pp)');
console.log();

// ═══════════════════════════════════════════════
// BASELINE B' (extended)
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  BASELINE B\' (extended)');
console.log('  log(a0) = b0 + b1*logMHI + b2*wig + b3*logMhost + b4*logSigma0 + b5*logMeanRun');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const bootBp = bootstrapCI(Y, X_Bp, 2000);
console.log('  Coefficients:');
console.log('    Intercept:   ' + fitBp.beta[0].toFixed(4) + ' (t=' + fitBp.tStats[0].toFixed(2) + ') CI=[' + bootBp[0].lo.toFixed(3) + ',' + bootBp[0].hi.toFixed(3) + ']');
for(let j=0;j<varNames_Bp.length;j++){
  console.log('    ' + varNames_Bp[j].padEnd(12) + ' ' + fitBp.beta[j+1].toFixed(4) + ' (t=' + fitBp.tStats[j+1].toFixed(2) + ') CI=[' + bootBp[j+1].lo.toFixed(3) + ',' + bootBp[j+1].hi.toFixed(3) + ']');
}
console.log();
console.log('  R2adj = ' + fitBp.r2adj.toFixed(4));
console.log('  SE = ' + fitBp.se.toFixed(4) + ' dex');
console.log('  LOO-RMS = ' + looBp.toFixed(4) + ' dex');
console.log('  LOO gap% = ' + gapV(looBp).toFixed(1) + '%');
console.log('  Residual tau = ' + fitBp.se.toFixed(3) + ' dex');
console.log();

console.log('  vs OLD Baseline B (envCode, same N=45):');
console.log('    R2adj: ' + fitB_old.r2adj.toFixed(4) + ' -> ' + fitBp.r2adj.toFixed(4) + ' (delta=' + (fitBp.r2adj-fitB_old.r2adj).toFixed(4) + ')');
console.log('    LOO:   ' + gapV(looB_old).toFixed(1) + '% -> ' + gapV(looBp).toFixed(1) + '% (delta=' + (gapV(looBp)-gapV(looB_old)).toFixed(1) + 'pp)');
console.log();

// ═══════════════════════════════════════════════
// CORRELATION MATRIX
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PREDICTOR CORRELATION MATRIX (N=' + N + ')');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const allVars = [
  {name:'logMHI', fn:g=>g.logMHI},
  {name:'wig', fn:g=>g.rcWiggliness},
  {name:'logMhost', fn:(g,i)=>logMhost[i]},
  {name:'logSigma0', fn:g=>g.logSigma0},
  {name:'logMeanRun', fn:g=>g.logMeanRun},
  {name:'logA0', fn:g=>g.logA0}
];
const vals = allVars.map(v => gals.map((g,i)=>v.fn(g,i)));

let header = '              ';
allVars.forEach(v => { header += v.name.padStart(10); });
console.log(header);
for(let i=0;i<allVars.length;i++){
  let row = '  ' + allVars[i].name.padEnd(12);
  for(let j=0;j<allVars.length;j++){
    if(j<=i) row += pearson(vals[i],vals[j]).r.toFixed(3).padStart(10);
    else row += ''.padStart(10);
  }
  console.log(row);
}
console.log();

// ═══════════════════════════════════════════════
// PARTIAL EFFECTS
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  PARTIAL EFFECTS (B\' model)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

for(let j=0;j<varNames_Bp.length;j++){
  const v = vals[j];
  const range = Math.max(...v) - Math.min(...v);
  const effect = fitBp.beta[j+1] * range;
  const iqr = (() => { const s=[...v].sort((a,b)=>a-b); return s[Math.floor(0.75*s.length)]-s[Math.floor(0.25*s.length)]; })();
  const iqrEffect = fitBp.beta[j+1] * iqr;
  console.log('  ' + varNames_Bp[j].padEnd(12) + ': coeff=' + fitBp.beta[j+1].toFixed(4) +
    ', range=' + range.toFixed(2) + ', full effect=' + effect.toFixed(3) + ' dex' +
    ', IQR effect=' + iqrEffect.toFixed(3) + ' dex');
}
console.log();

// ═══════════════════════════════════════════════
// RESIDUAL ANALYSIS
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  RESIDUAL ANALYSIS (B\' model)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const residBp = fitBp.resid;
const absResid = residBp.map(Math.abs).sort((a,b)=>b-a);
console.log('  Residual distribution:');
console.log('    Mean = ' + mean(residBp).toFixed(4));
console.log('    SD = ' + sd(residBp).toFixed(4) + ' dex');
console.log('    Min = ' + Math.min(...residBp).toFixed(4));
console.log('    Max = ' + Math.max(...residBp).toFixed(4));
console.log('    |resid| > 0.3 dex: ' + residBp.filter(r=>Math.abs(r)>0.3).length + '/' + N);
console.log();

console.log('  Top 10 residuals:');
const sorted = gals.map((g,i) => ({name:g.name, resid:residBp[i], logA0:g.logA0}))
  .sort((a,b) => Math.abs(b.resid)-Math.abs(a.resid)).slice(0,10);
sorted.forEach(s => {
  console.log('    ' + s.name.padEnd(15) + ' resid=' + s.resid.toFixed(3) + 
    ' logA0=' + s.logA0.toFixed(3) + ' D=' + (distMap[s.name]||0).toFixed(1) + ' Mpc');
});
console.log();

// ═══════════════════════════════════════════════
// FULL COMPARISON TABLE
// ═══════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  COMPARISON: OLD vs NEW BASELINES');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('                   OLD (envCode, N=45)    NEW (logMhost, N=45)');
console.log('  ────────────────────────────────────────────────────────────');
console.log('  A R2adj:         ' + fitA_old.r2adj.toFixed(4).padStart(8) + '                   ' + fitAp.r2adj.toFixed(4).padStart(8));
console.log('  A LOO gap%:      ' + gapV(looA_old).toFixed(1).padStart(6) + '%                   ' + gapV(looAp).toFixed(1).padStart(6) + '%');
console.log('  B R2adj:         ' + fitB_old.r2adj.toFixed(4).padStart(8) + '                   ' + fitBp.r2adj.toFixed(4).padStart(8));
console.log('  B LOO gap%:      ' + gapV(looB_old).toFixed(1).padStart(6) + '%                   ' + gapV(looBp).toFixed(1).padStart(6) + '%');
console.log('  B SE (tau):      ' + fitB_old.se.toFixed(4).padStart(8) + '                   ' + fitBp.se.toFixed(4).padStart(8));
console.log();

// ═══════════════════════════════════════════════
// SAVE FROZEN BASELINES
// ═══════════════════════════════════════════════
const output = {
  phase: '58a4',
  frozenDate: '2026-04-04',
  sample: { n: N, description: 'Published logMhost subsample (UNGC+KT2017)' },
  sdY,
  baselineAp: {
    formula: 'log(a0) = b0 + b1*logMHI + b2*rcWiggliness + b3*logMhost + b4*logSigma0',
    variables: ['logMHI','rcWiggliness','logMhost','logSigma0'],
    coefficients: fitAp.beta,
    se: fitAp.seBeta,
    tStats: fitAp.tStats,
    bootstrapCI: bootAp,
    r2adj: fitAp.r2adj,
    se_model: fitAp.se,
    looRMS: looAp,
    looGap: gapV(looAp)
  },
  baselineBp: {
    formula: 'log(a0) = b0 + b1*logMHI + b2*rcWiggliness + b3*logMhost + b4*logSigma0 + b5*logMeanRun',
    variables: ['logMHI','rcWiggliness','logMhost','logSigma0','logMeanRun'],
    coefficients: fitBp.beta,
    se: fitBp.seBeta,
    tStats: fitBp.tStats,
    bootstrapCI: bootBp,
    r2adj: fitBp.r2adj,
    se_model: fitBp.se,
    looRMS: looBp,
    looGap: gapV(looBp)
  },
  oldBaselines: {
    A: { r2adj: fitA_old.r2adj, looGap: gapV(looA_old) },
    B: { r2adj: fitB_old.r2adj, looGap: gapV(looB_old) }
  },
  galaxyResiduals: gals.map((g,i) => ({
    name: g.name, logA0: g.logA0, logMhost: logMhost[i],
    residAp: fitAp.resid[i], residBp: fitBp.resid[i]
  }))
};
fs.writeFileSync('public/phase58a4-frozen-baselines.json', JSON.stringify(output, null, 2));
console.log('  Frozen baselines saved to public/phase58a4-frozen-baselines.json');
