const fs = require('fs');
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const tidalData = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8')).tidalData;

const gals = stageA.galaxies;
const N = gals.length;
const tdMap = {};
tidalData.forEach(t => { tdMap[t.name] = t; });

function mean(a){if(!a.length)return NaN;return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){if(a.length<2)return NaN;const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function pearson(x,y){const n=x.length;if(n<4)return{r:NaN,t:NaN,n};const mx=mean(x),my=mean(y);let sxy=0,sxx=0,syy=0;for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;syy+=(y[i]-my)**2;}if(sxx<1e-15||syy<1e-15)return{r:NaN,t:NaN,n};const r=sxy/Math.sqrt(sxx*syy);return{r,t:r*Math.sqrt((n-2)/(1-r*r)),n};}
function spearmanR(x,y){const n=x.length;function rank(a){const s=a.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);const r=new Array(n);for(let i=0;i<n;){let j=i;while(j<n&&s[j].v===s[i].v)j++;const avg=(i+j-1)/2+1;for(let k=i;k<j;k++)r[s[k].i]=avg;i=j;}return r;}const rx=rank(x),ry=rank(y);return pearson(rx,ry).r;}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se:Math.sqrt(rss/(n-p)),resid,rss,tss,n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function tTest2(a,b){const na=a.length,nb=b.length,ma=mean(a),mb=mean(b);const sa=sd(a),sb=sd(b);const sp=Math.sqrt(((na-1)*sa*sa+(nb-1)*sb*sb)/(na+nb-2));return{diff:ma-mb,t:(ma-mb)/(sp*Math.sqrt(1/na+1/nb)),df:na+nb-2};}

const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));
const distMap = {};
sparc.forEach(s => { distMap[s.name] = s.D; });

const pub = gals.filter(g => tdMap[g.name].quality === 'published');
const est = gals.filter(g => tdMap[g.name].quality === 'crude');
const logMhost_pub = pub.map(g => tdMap[g.name].logMhost);
const logMhost_all = gals.map(g => tdMap[g.name].logMhost);

function gapV(rms, sdv) { return 100*(1-rms**2/sdv**2); }

console.log('======================================================================');
console.log('  PHASE 58a3: logMhost ROBUSTNESS / MISSINGNESS AUDIT');
console.log('======================================================================\n');

// ═══════════════════════════════════════════════════════════════════
// 1) MISSINGNESS BIAS TEST
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  1) MISSINGNESS BIAS TEST');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
console.log('  N_published = ' + pub.length);
console.log('  N_estimated = ' + est.length + '\n');

const vars = [
  { name: 'distance', fn: g => distMap[g.name] || 0 },
  { name: 'envCode', fn: g => g.envCode },
  { name: 'logMHI', fn: g => g.logMHI },
  { name: 'rcWiggliness', fn: g => g.rcWiggliness },
  { name: 'logMeanRun', fn: g => g.logMeanRun },
  { name: 'logSigma0', fn: g => g.logSigma0 },
  { name: 'delta_a0', fn: g => g.delta_a0 },
  { name: 'logA0', fn: g => g.logA0 },
];

vars.forEach(v => {
  const a = pub.map(v.fn), b = est.map(v.fn);
  const tt = tTest2(a, b);
  const sig = Math.abs(tt.t) > 2.0 ? ' ***' : Math.abs(tt.t) > 1.65 ? ' *' : '';
  console.log('  ' + v.name.padEnd(15) + 
    ' pub=' + mean(a).toFixed(3).padStart(7) + 
    ' est=' + mean(b).toFixed(3).padStart(7) + 
    ' diff=' + tt.diff.toFixed(3).padStart(7) + 
    ' t=' + tt.t.toFixed(2).padStart(6) + sig);
});

console.log('\n  envCode distribution:');
[0,1,2].forEach(e => {
  const np = pub.filter(g=>g.envCode===e).length;
  const ne = est.filter(g=>g.envCode===e).length;
  console.log('    envCode=' + e + ': pub=' + np + ', est=' + ne);
});
console.log();

// ═══════════════════════════════════════════════════════════════════
// 2) SAME-SAMPLE COMPARISON (published-only N=45)
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  2) SAME-SAMPLE COMPARISON (N=' + pub.length + ')');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const Yp = pub.map(g => g.logA0);
const sdYp = sd(Yp);
const XAp = pub.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0]);
const XBp = pub.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun]);
const XAtdp = pub.map((g,i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, logMhost_pub[i]]);
const XBtdp = pub.map((g,i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun, logMhost_pub[i]]);

const fAp = ols(Yp, XAp);
const fAtdp = ols(Yp, XAtdp);
const fBp = ols(Yp, XBp);
const fBtdp = ols(Yp, XBtdp);

const looAp = looCV(Yp, XAp);
const looAtdp = looCV(Yp, XAtdp);
const looBp = looCV(Yp, XBp);
const looBtdp = looCV(Yp, XBtdp);

console.log('  Baseline A on N=45:');
console.log('    R2adj = ' + fAp.r2adj.toFixed(4));
console.log('    LOO gap% = ' + gapV(looAp, sdYp).toFixed(1) + '%');
console.log('  A + logMhost on N=45:');
console.log('    R2adj = ' + fAtdp.r2adj.toFixed(4));
console.log('    LOO gap% = ' + gapV(looAtdp, sdYp).toFixed(1) + '%');
console.log('    delta R2adj = ' + (fAtdp.r2adj - fAp.r2adj).toFixed(4));
console.log('    delta LOO = ' + (gapV(looAtdp, sdYp) - gapV(looAp, sdYp)).toFixed(1) + 'pp');
console.log('    logMhost t = ' + fAtdp.tStats[5].toFixed(3));
console.log();
console.log('  Baseline B on N=45:');
console.log('    R2adj = ' + fBp.r2adj.toFixed(4));
console.log('    LOO gap% = ' + gapV(looBp, sdYp).toFixed(1) + '%');
console.log('  B + logMhost on N=45:');
console.log('    R2adj = ' + fBtdp.r2adj.toFixed(4));
console.log('    LOO gap% = ' + gapV(looBtdp, sdYp).toFixed(1) + '%');
console.log('    delta R2adj = ' + (fBtdp.r2adj - fBp.r2adj).toFixed(4));
console.log('    delta LOO = ' + (gapV(looBtdp, sdYp) - gapV(looBp, sdYp)).toFixed(1) + 'pp');
console.log('    logMhost t = ' + fBtdp.tStats[6].toFixed(3));
console.log();

// ═══════════════════════════════════════════════════════════════════
// 3) WITHIN NON-FIELD ONLY (envCode > 0)
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  3) WITHIN-ENVIRONMENT STRIPS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Non-field only (envCode > 0)
const nonField = pub.filter(g => g.envCode > 0);
const nfY = nonField.map(g => g.logA0);
const nfM = nonField.map(g => tdMap[g.name].logMhost);
const nfRaw = pearson(nfM, nfY);
console.log('  Non-field only (envCode > 0, N=' + nonField.length + '):');
console.log('    r(logMhost, logA0) = ' + nfRaw.r.toFixed(4) + ' (t=' + nfRaw.t.toFixed(2) + ')');
console.log('    Spearman = ' + spearmanR(nfM, nfY).toFixed(4));

// After logMHI within non-field
if (nonField.length >= 6) {
  const nfFitMHI = ols(nfY, nonField.map(g => [g.logMHI]));
  const nfAfterMHI = pearson(nfM, nfFitMHI.resid);
  console.log('    After logMHI: r=' + nfAfterMHI.r.toFixed(4) + ' (t=' + nfAfterMHI.t.toFixed(2) + ')');
}
console.log();

// Field only (envCode = 0)
const field = pub.filter(g => g.envCode === 0);
const fY = field.map(g => g.logA0);
const fM = field.map(g => tdMap[g.name].logMhost);
const fRaw = pearson(fM, fY);
console.log('  Field only (envCode = 0, N=' + field.length + '):');
console.log('    r(logMhost, logA0) = ' + fRaw.r.toFixed(4) + ' (t=' + fRaw.t.toFixed(2) + ')');
console.log('    Spearman = ' + spearmanR(fM, fY).toFixed(4));
if (field.length >= 6) {
  const ffMHI = ols(fY, field.map(g => [g.logMHI]));
  const fAfterMHI = pearson(fM, ffMHI.resid);
  console.log('    After logMHI: r=' + fAfterMHI.r.toFixed(4) + ' (t=' + fAfterMHI.t.toFixed(2) + ')');
}
console.log();

// Group/cluster (envCode >= 1)
const grpClust = pub.filter(g => g.envCode >= 1);
const gcY = grpClust.map(g => g.logA0);
const gcM = grpClust.map(g => tdMap[g.name].logMhost);
const gcRaw = pearson(gcM, gcY);
console.log('  Group+Cluster (envCode >= 1, N=' + grpClust.length + '):');
console.log('    r(logMhost, logA0) = ' + gcRaw.r.toFixed(4) + ' (t=' + gcRaw.t.toFixed(2) + ')');
console.log();

// Cluster only (envCode = 2)
const cluster = pub.filter(g => g.envCode === 2);
const cY = cluster.map(g => g.logA0);
const cM = cluster.map(g => tdMap[g.name].logMhost);
const cRaw = pearson(cM, cY);
console.log('  Cluster only (envCode = 2, N=' + cluster.length + '):');
console.log('    r(logMhost, logA0) = ' + cRaw.r.toFixed(4) + ' (t=' + cRaw.t.toFixed(2) + ')');
console.log();

// ═══════════════════════════════════════════════════════════════════
// 4) DISTANCE-QUALITY CONTROL
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  4) DISTANCE-QUALITY CONTROL');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const dist_pub = pub.map(g => distMap[g.name] || 0);

// r(logMhost, distance)
const rMD = pearson(logMhost_pub, dist_pub);
console.log('  r(logMhost, distance) = ' + rMD.r.toFixed(4) + ' (t=' + rMD.t.toFixed(2) + ')');

// After distance
const fDist = ols(Yp, pub.map((g,i) => [dist_pub[i]]));
const aDist = pearson(logMhost_pub, fDist.resid);
console.log('  After distance: r=' + aDist.r.toFixed(4) + ' (t=' + aDist.t.toFixed(2) + ')');

// After envCode + distance
const fED = ols(Yp, pub.map((g,i) => [g.envCode, dist_pub[i]]));
const aED = pearson(logMhost_pub, fED.resid);
console.log('  After envCode+distance: r=' + aED.r.toFixed(4) + ' (t=' + aED.t.toFixed(2) + ')');

// After Baseline B + distance
const XBD = pub.map((g,i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun, dist_pub[i]]);
const fBD = ols(Yp, XBD);
const aBD = pearson(logMhost_pub, fBD.resid);
console.log('  After Baseline B + distance: r=' + aBD.r.toFixed(4) + ' (t=' + aBD.t.toFixed(2) + ')');
console.log();

// Full model: B + distance + logMhost
const XBDtd = pub.map((g,i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun, dist_pub[i], logMhost_pub[i]]);
const fBDtd = ols(Yp, XBDtd);
console.log('  B + distance + logMhost:');
console.log('    logMhost t = ' + fBDtd.tStats[7].toFixed(3));
console.log('    distance t = ' + fBDtd.tStats[6].toFixed(3));
console.log('    R2adj = ' + fBDtd.r2adj.toFixed(4) + ' (vs B: ' + fBp.r2adj.toFixed(4) + ')');
console.log();

// Near vs far split within published
const nearCut = 20;
const near = pub.filter((g,i) => dist_pub[i] <= nearCut);
const far = pub.filter((g,i) => dist_pub[i] > nearCut);
console.log('  Near (D<=' + nearCut + ' Mpc, N=' + near.length + '):');
const nearM = near.map(g => tdMap[g.name].logMhost);
const nearY = near.map(g => g.logA0);
const nearR = pearson(nearM, nearY);
console.log('    r(logMhost, logA0) = ' + nearR.r.toFixed(4) + ' (t=' + nearR.t.toFixed(2) + ')');

console.log('  Far (D>' + nearCut + ' Mpc, N=' + far.length + '):');
const farM = far.map(g => tdMap[g.name].logMhost);
const farY = far.map(g => g.logA0);
const farR = pearson(farM, farY);
console.log('    r(logMhost, logA0) = ' + farR.r.toFixed(4) + ' (t=' + farR.t.toFixed(2) + ')');
console.log();

// ═══════════════════════════════════════════════════════════════════
// 5) ADDITIONAL: Does logMhost REPLACE envCode or ADD to it?
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  5) REPLACE vs ADD diagnostic');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// B with both envCode AND logMhost
console.log('  B + logMhost (with envCode kept):');
console.log('    envCode t = ' + fBtdp.tStats[3].toFixed(3));
console.log('    logMhost t = ' + fBtdp.tStats[6].toFixed(3));
console.log('    R2adj = ' + fBtdp.r2adj.toFixed(4));
const looBoth = looCV(Yp, XBtdp);
console.log('    LOO gap% = ' + gapV(looBoth, sdYp).toFixed(1) + '%');
console.log();

// B with logMhost REPLACING envCode
const XBrep = pub.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost_pub[i], g.logSigma0, g.logMeanRun]);
const fBrep = ols(Yp, XBrep);
const looBrep = looCV(Yp, XBrep);
console.log('  B with envCode REPLACED by logMhost:');
console.log('    logMhost t = ' + fBrep.tStats[3].toFixed(3));
console.log('    R2adj = ' + fBrep.r2adj.toFixed(4));
console.log('    LOO gap% = ' + gapV(looBrep, sdYp).toFixed(1) + '%');
console.log();

console.log('  Comparison (same N=45):');
console.log('    B (envCode):         LOO = ' + gapV(looBp, sdYp).toFixed(1) + '%');
console.log('    B (logMhost replace): LOO = ' + gapV(looBrep, sdYp).toFixed(1) + '%');
console.log('    B + logMhost (both):  LOO = ' + gapV(looBoth, sdYp).toFixed(1) + '%');
console.log('    delta (replace):  ' + (gapV(looBrep, sdYp) - gapV(looBp, sdYp)).toFixed(1) + 'pp');
console.log('    delta (add):      ' + (gapV(looBoth, sdYp) - gapV(looBp, sdYp)).toFixed(1) + 'pp');
console.log();

if (gapV(looBoth, sdYp) > gapV(looBrep, sdYp) + 1) {
  console.log('  -> BOTH together better: envCode adds info beyond logMhost');
} else if (gapV(looBrep, sdYp) > gapV(looBoth, sdYp) + 1) {
  console.log('  -> REPLACEMENT better: envCode adds noise when logMhost present');
} else {
  console.log('  -> COMPARABLE: envCode and logMhost capture overlapping info');
}
console.log();

// ═══════════════════════════════════════════════════════════════════
// FINAL VERDICT
// ═══════════════════════════════════════════════════════════════════
console.log('======================================================================');
console.log('  FINAL ROBUSTNESS VERDICT');
console.log('======================================================================\n');

const missOK = Math.abs(tTest2(pub.map(g=>g.delta_a0), est.map(g=>g.delta_a0)).t) < 2.0;
const sameSampleOK = fBtdp.tStats[6] < -1.65 && gapV(looBtdp, sdYp) > gapV(looBp, sdYp);
const nonFieldOK = nonField.length >= 10 && Math.abs(nfRaw.t) > 1.0;
const distOK = Math.abs(aBD.t) > 1.0;

console.log('  1) Missingness bias: ' + (missOK ? 'PASS' : 'CONCERN') + 
  ' (delta_a0 diff t=' + tTest2(pub.map(g=>g.delta_a0), est.map(g=>g.delta_a0)).t.toFixed(2) + ')');
console.log('  2) Same-sample:     ' + (sameSampleOK ? 'PASS' : 'FAIL') + 
  ' (B+td t=' + fBtdp.tStats[6].toFixed(2) + ', LOO delta=' + 
  (gapV(looBtdp, sdYp) - gapV(looBp, sdYp)).toFixed(1) + 'pp)');
console.log('  3) Non-field strip: ' + (nonFieldOK ? 'PASS' : 'WEAK') + 
  ' (N=' + nonField.length + ', r=' + nfRaw.r.toFixed(3) + ', t=' + nfRaw.t.toFixed(2) + ')');
console.log('  4) Distance control: ' + (distOK ? 'PASS' : 'FAIL') + 
  ' (after B+dist: r=' + aBD.r.toFixed(3) + ', t=' + aBD.t.toFixed(2) + ')');
console.log();

let passed = [missOK, sameSampleOK, nonFieldOK, distOK].filter(Boolean).length;
let verdict;
if (passed === 4) verdict = 'CONFIRMED-ROBUST';
else if (passed >= 3) verdict = 'PARTIAL-ROBUST';
else if (passed >= 2) verdict = 'PARTIAL';
else verdict = 'FAIL';

console.log('  ═════════════════════════════════════════════════');
console.log('  VERDICT: ' + verdict + ' (' + passed + '/4 tests passed)');
console.log('  ═════════════════════════════════════════════════');

// Save
const output = {
  phase: '58a3',
  missingness: {
    nPub: pub.length, nEst: est.length,
    tests: vars.map(v => {
      const a = pub.map(v.fn), b = est.map(v.fn);
      const tt = tTest2(a, b);
      return { name: v.name, pubMean: mean(a), estMean: mean(b), t: tt.t };
    })
  },
  sameSample: {
    n: pub.length,
    A: { r2adj: fAp.r2adj, loo: gapV(looAp, sdYp) },
    Atd: { r2adj: fAtdp.r2adj, loo: gapV(looAtdp, sdYp), t: fAtdp.tStats[5] },
    B: { r2adj: fBp.r2adj, loo: gapV(looBp, sdYp) },
    Btd: { r2adj: fBtdp.r2adj, loo: gapV(looBtdp, sdYp), t: fBtdp.tStats[6] }
  },
  strips: {
    nonField: { n: nonField.length, r: nfRaw.r, t: nfRaw.t },
    field: { n: field.length, r: fRaw.r, t: fRaw.t },
    cluster: { n: cluster.length, r: cRaw.r, t: cRaw.t }
  },
  distanceControl: {
    afterDist: aDist, afterEnvDist: aED, afterBDist: aBD,
    bDistTd: { logMhost_t: fBDtd.tStats[7], dist_t: fBDtd.tStats[6] }
  },
  replaceVsAdd: {
    B_envCode_loo: gapV(looBp, sdYp),
    B_replace_loo: gapV(looBrep, sdYp),
    B_both_loo: gapV(looBoth, sdYp)
  },
  verdict
};
fs.writeFileSync('public/phase58a3-robustness.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase58a3-robustness.json');
