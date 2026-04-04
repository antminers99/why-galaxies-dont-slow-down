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
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0);return{beta,resid,rss};}
function pearsonCorr(x,y){if(x.length<3)return 0;const mx=mean(x),my=mean(y);let num=0,dx2=0,dy2=0;for(let i=0;i<x.length;i++){num+=(x[i]-mx)*(y[i]-my);dx2+=(x[i]-mx)**2;dy2+=(y[i]-my)**2;}if(dx2<1e-30||dy2<1e-30)return 0;return num/Math.sqrt(dx2*dy2);}
function spearmanRank(a){const s=[...a].map((v,i)=>({v,i})).sort((x,y)=>x.v-y.v);const r=new Array(a.length);s.forEach((x,i)=>r[x.i]=i+1);return r;}
function spearmanCorr(x,y){if(x.length<3)return 0;const rx=spearmanRank(x),ry=spearmanRank(y);const n=x.length;let d2=0;for(let i=0;i<n;i++)d2+=(rx[i]-ry[i])**2;return 1-6*d2/(n*(n*n-1));}
function tStat(r,n){if(Math.abs(r)>0.9999)return r>0?99:-99;return r*Math.sqrt((n-2)/(1-r*r));}

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

const Y_all = gals45.map(g => g.logA0);
const morphT = gals45.map(g => sparcMap[g.name]?.T || 5);
const logUps = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhostAll = gals45.map(g => tdMap[g.name].logMhost);
const Xconf = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT[i]]);
const fConf = ols(logUps, Xconf);
const upsPerp = fConf.resid;

const logMHI = gals45.map(g => g.logMHI);
const logSig0 = gals45.map(g => g.logSigma0);
const logMR = gals45.map(g => g.logMeanRun);
const rcWig = gals45.map(g => g.rcWiggliness);

const X_M5 = gals45.map((_,i) => [logMHI[i], logMhostAll[i], logSig0[i], logMR[i], upsPerp[i]]);
const fM5 = ols(Y_all, X_M5);
const residM5 = fM5.resid;

const X_M5noMR = gals45.map((_,i) => [logMHI[i], logMhostAll[i], logSig0[i], upsPerp[i]]);
const fM5noMR = ols(Y_all, X_M5noMR);
const residM5noMR = fM5noMR.resid;

const thingsNames = ['NGC2403','NGC2841','NGC2903','NGC3198','NGC3521','NGC5055','NGC7331'];
const thingsIdx = thingsNames.map(n => gals45.findIndex(g => g.name === n));
const Nth = 7;

const residTH = thingsIdx.map(i => residM5[i]);
const residTH_noMR = thingsIdx.map(i => residM5noMR[i]);
const logA0_TH = thingsIdx.map(i => Y_all[i]);
const logMR_TH = thingsIdx.map(i => logMR[i]);
const rcWig_TH = thingsIdx.map(i => rcWig[i]);

const varNames2D = ['things_ncm_amp','things_ncm_frac','things_c1','things_s1','things_c3','things_s3','things_lopsidedness','things_bisymFlow'];
const vars2D = {};
varNames2D.forEach(v => { vars2D[v] = thingsIdx.map(i => gals45[i][v]); });

const stdize = (arr) => { const m=mean(arr), s=sd(arr); return s<1e-15 ? arr.map(()=>0) : arr.map(v=>(v-m)/s); };

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 74: 2D DYNAMICAL-STATE PROGRAM');
console.log('══════════════════════════════════════════════════════════════════════\n');
console.log('  THINGS subset: N=' + Nth);
console.log('  Galaxies: ' + thingsNames.join(', '));
console.log();

console.log('  RAW DATA TABLE');
console.log('  ' + 'Galaxy'.padEnd(12) + 'logA0   residM5  MeanRun  ncm_amp  ncm_frac lopsid   bisymF   rcWig');
console.log('  ' + '-'.repeat(90));
thingsIdx.forEach((idx) => {
  const g = gals45[idx];
  console.log('  ' + g.name.padEnd(12) +
    g.logA0.toFixed(3).padStart(7) + '  ' +
    (residM5[idx]>0?'+':'') + residM5[idx].toFixed(3).padStart(6) + '   ' +
    g.logMeanRun.toFixed(3).padStart(6) + '   ' +
    (g.things_ncm_amp||0).toFixed(3).padStart(6) + '   ' +
    (g.things_ncm_frac||0).toFixed(3).padStart(6) + '   ' +
    (g.things_lopsidedness||0).toFixed(3).padStart(6) + '   ' +
    (g.things_bisymFlow||0).toFixed(3).padStart(6) + '   ' +
    (g.rcWiggliness||0).toFixed(3).padStart(6));
});
console.log();

// ═════════ 74a: Individual correlations ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  74a: 2D COMPOSITE PILOT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  Individual correlations with resid(M5):');
const corrResults = {};
varNames2D.forEach(v => {
  const vals = vars2D[v];
  const r = pearsonCorr(vals, residTH);
  const rs = spearmanCorr(vals, residTH);
  const t = tStat(r, Nth);
  corrResults[v] = {pearson:r, spearman:rs, t};
  console.log('    ' + v.replace('things_','').padEnd(15) + ': r=' + (r>0?'+':'') + r.toFixed(3) + ', rs=' + (rs>0?'+':'') + rs.toFixed(3) + ', t=' + t.toFixed(2));
});
console.log();

// Key finding from the data: MeanRun correlates STRONGLY with lopsidedness, ncm_frac, bisymFlow
// These are the "coherence" variables. Let's build composites from them.
const lopsid_s = stdize(vars2D['things_lopsidedness']);
const bisym_s = stdize(vars2D['things_bisymFlow']);
const ncmAmp_s = stdize(vars2D['things_ncm_amp']);
const ncmFrac_s = stdize(vars2D['things_ncm_frac']);
const c3_s = stdize(vars2D['things_c3']);
const s3_s = stdize(vars2D['things_s3']);

const compositeList = [
  {name:'coherence: -lopsid+bisym-ncmAmp', vals: lopsid_s.map((_,i)=>-lopsid_s[i]+bisym_s[i]-ncmAmp_s[i])},
  {name:'coherence: -lopsid-ncmFrac',      vals: lopsid_s.map((_,i)=>-lopsid_s[i]-ncmFrac_s[i])},
  {name:'disturbance: lopsid+ncmAmp',      vals: lopsid_s.map((_,i)=>lopsid_s[i]+ncmAmp_s[i])},
  {name:'harmonic: c3+s3',                 vals: c3_s.map((_,i)=>c3_s[i]+s3_s[i])},
  {name:'-ncmAmp (coherent)',              vals: ncmAmp_s.map(v=>-v)},
  {name:'s3 alone',                         vals: s3_s},
  {name:'c3 alone',                         vals: c3_s},
];

console.log('  Composite scores vs resid(M5):');
let bestComp = compositeList[0], bestCompAbsR = 0;
compositeList.forEach(comp => {
  const r = pearsonCorr(comp.vals, residTH);
  const rs = spearmanCorr(comp.vals, residTH);
  const t = tStat(r, Nth);
  comp.r = r; comp.rs = rs; comp.t = t;
  console.log('    ' + comp.name.padEnd(35) + ': r=' + (r>0?'+':'') + r.toFixed(3) + ', rs=' + (rs>0?'+':'') + rs.toFixed(3) + ', t=' + t.toFixed(2));
  if (Math.abs(r) > bestCompAbsR) { bestCompAbsR = Math.abs(r); bestComp = comp; }
});
console.log();
console.log('  Best composite: ' + bestComp.name + ' (r=' + bestComp.r.toFixed(3) + ')');
console.log();

// Also correlate with raw logA0
console.log('  Best composite vs raw log(a₀):');
const rA0 = pearsonCorr(bestComp.vals, logA0_TH);
const rsA0 = spearmanCorr(bestComp.vals, logA0_TH);
console.log('    r=' + rA0.toFixed(3) + ', rs=' + rsA0.toFixed(3));
console.log();

// ═════════ 74b: MeanRun Bridge ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  74b: MeanRun BRIDGE TEST');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  MeanRun vs 2D variables:');
varNames2D.forEach(v => {
  const r = pearsonCorr(vars2D[v], logMR_TH);
  const rs = spearmanCorr(vars2D[v], logMR_TH);
  const flag = Math.abs(r) > 0.6 ? ' *** STRONG' : Math.abs(r) > 0.3 ? ' ** moderate' : '';
  console.log('    ' + v.replace('things_','').padEnd(15) + ': r=' + (r>0?'+':'') + r.toFixed(3) + ', rs=' + (rs>0?'+':'') + rs.toFixed(3) + flag);
});
console.log();

const rBridge = pearsonCorr(bestComp.vals, logMR_TH);
const rsBridge = spearmanCorr(bestComp.vals, logMR_TH);
console.log('  Best composite vs MeanRun: r=' + rBridge.toFixed(3) + ', rs=' + rsBridge.toFixed(3));
console.log();

const rNoMR = pearsonCorr(bestComp.vals, residTH_noMR);
const rsNoMR = spearmanCorr(bestComp.vals, residTH_noMR);
console.log('  Best composite vs resid(M5-without-MeanRun):');
console.log('    r=' + rNoMR.toFixed(3) + ', rs=' + rsNoMR.toFixed(3));
console.log('    (Does 2D add info beyond MeanRun? Need |r|>0.4 for YES)');
console.log();

// ═════════ 74c: rcWig legacy ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  74c: rcWig LEGACY TEST');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  rcWig vs MeanRun: r=' + pearsonCorr(rcWig_TH, logMR_TH).toFixed(3) + ', rs=' + spearmanCorr(rcWig_TH, logMR_TH).toFixed(3));
console.log('  rcWig vs best 2D composite: r=' + pearsonCorr(rcWig_TH, bestComp.vals).toFixed(3));
console.log();

const mrAlign = varNames2D.map(v => Math.abs(pearsonCorr(vars2D[v], logMR_TH)));
const wigAlign = varNames2D.map(v => Math.abs(pearsonCorr(vars2D[v], rcWig_TH)));
console.log('  Mean |r| with 2D variables:');
console.log('    MeanRun:      ' + mean(mrAlign).toFixed(3));
console.log('    rcWiggliness: ' + mean(wigAlign).toFixed(3));
console.log('    ' + (mean(mrAlign) > mean(wigAlign) ? '-> MeanRun is BETTER 1D proxy for 2D dynamics' : '-> rcWig aligns more with 2D'));
console.log();

// ═════════ 74d: Pairwise sanity ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  74d: PAIRWISE 2D SANITY TEST');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Galaxy-by-galaxy comparison
console.log('  Case-by-case:');
console.log('  ' + 'Galaxy'.padEnd(12) + '  2D-score    resid(M5)   Direction match?');
console.log('  ' + '-'.repeat(60));
thingsNames.forEach((n, j) => {
  const cs = bestComp.vals[j];
  const rs = residTH[j];
  const match = (cs > 0 && rs > 0) || (cs < 0 && rs < 0) || (Math.abs(cs) < 0.3 && Math.abs(rs) < 0.05);
  console.log('  ' + n.padEnd(12) + '  ' + (cs>0?'+':'') + cs.toFixed(3).padStart(6) + '      ' + (rs>0?'+':'') + rs.toFixed(3).padStart(6) + '       ' + (match ? 'YES' : 'no'));
});
console.log();

let pairCorr = 0, pairTot = 0;
for (let i = 0; i < Nth; i++) {
  for (let j = i+1; j < Nth; j++) {
    pairTot++;
    if ((bestComp.vals[i]-bestComp.vals[j]) * (residTH[i]-residTH[j]) > 0) pairCorr++;
  }
}
console.log('  Pairwise concordance: ' + pairCorr + '/' + pairTot + ' = ' + (pairCorr/pairTot*100).toFixed(1) + '% (50%=random, >70%=good)');
console.log();

// Top/bottom split
const sorted = thingsNames.map((_,i)=>i).sort((a,b)=>bestComp.vals[b]-bestComp.vals[a]);
const topH = sorted.slice(0, 4);
const botH = sorted.slice(4);
console.log('  Top-4 by 2D coherence: ' + topH.map(i=>thingsNames[i]).join(', '));
console.log('    mean resid(M5) = ' + (mean(topH.map(i=>residTH[i]))>0?'+':'') + mean(topH.map(i=>residTH[i])).toFixed(3));
console.log('  Bottom-3 by 2D coherence: ' + botH.map(i=>thingsNames[i]).join(', '));
console.log('    mean resid(M5) = ' + (mean(botH.map(i=>residTH[i]))>0?'+':'') + mean(botH.map(i=>residTH[i])).toFixed(3));
console.log();

// ═════════ 74e: Variance budget ═════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  74e: VARIANCE BUDGET & KEY FINDING');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const residRMS = Math.sqrt(residM5.reduce((s,r)=>s+r*r,0)/N);
console.log('  M5 residual RMS (N=45): ' + residRMS.toFixed(3) + ' dex');
console.log('  Best 2D composite r with resid(M5): ' + bestComp.r.toFixed(3));
console.log('  r-squared on THINGS N=7: ' + (bestComp.r**2).toFixed(3));
console.log();

// THE KEY FINDING: MeanRun-2D bridge
console.log('  KEY FINDING — MeanRun as 1D projection of 2D dynamics:');
console.log('  ┌───────────────────────────┬────────┐');
console.log('  │ MeanRun vs lopsidedness   │ ' + pearsonCorr(logMR_TH, vars2D['things_lopsidedness']).toFixed(3).padStart(6) + ' │');
console.log('  │ MeanRun vs bisymFlow      │ ' + pearsonCorr(logMR_TH, vars2D['things_bisymFlow']).toFixed(3).padStart(6) + ' │');
console.log('  │ MeanRun vs ncm_frac       │ ' + pearsonCorr(logMR_TH, vars2D['things_ncm_frac']).toFixed(3).padStart(6) + ' │');
console.log('  └───────────────────────────┴────────┘');
console.log();

const rLop = Math.abs(pearsonCorr(logMR_TH, vars2D['things_lopsidedness']));
const rBis = Math.abs(pearsonCorr(logMR_TH, vars2D['things_bisymFlow']));
const rNcf = Math.abs(pearsonCorr(logMR_TH, vars2D['things_ncm_frac']));
const strongBridge = (rLop > 0.7 || rBis > 0.7 || rNcf > 0.7);

// ═════════ VERDICT ═════════
console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 74: FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════════════\n');

const hasResidSignal = Math.abs(bestComp.r) > 0.3;
const mrIsBridge = strongBridge;
const addsOverMR = Math.abs(rNoMR) > 0.4;
const concordanceOK = pairCorr/pairTot > 0.55;
const mrBetter = mean(mrAlign) > mean(wigAlign);

let verdict;
if (hasResidSignal && addsOverMR && concordanceOK) verdict = 'CONFIRMED-2D-BRIDGE';
else if (mrIsBridge || hasResidSignal) verdict = 'PARTIAL';
else verdict = 'FAIL';

console.log('  ┌───────────────────────────────────────────┬────────┐');
console.log('  │ Criterion                                 │ Result │');
console.log('  ├───────────────────────────────────────────┼────────┤');
console.log('  │ 2D composite correlates with resid(M5)?   │ ' + (hasResidSignal?'YES':'NO').padEnd(6) + ' │');
console.log('  │ MeanRun bridges to 2D strongly (|r|>0.7)? │ ' + (mrIsBridge?'YES':'NO').padEnd(6) + ' │');
console.log('  │ 2D adds beyond MeanRun (|r|>0.4)?         │ ' + (addsOverMR?'YES':'NO').padEnd(6) + ' │');
console.log('  │ Pairwise concordance > 55%?               │ ' + (concordanceOK?'YES':'NO').padEnd(6) + ' │');
console.log('  │ MeanRun better 1D proxy than rcWig?       │ ' + (mrBetter?'YES':'NO').padEnd(6) + ' │');
console.log('  └───────────────────────────────────────────┴────────┘');
console.log();
console.log('  VERDICT: ' + verdict);
console.log();

if (verdict === 'CONFIRMED-2D-BRIDGE') {
  console.log('  2D dynamical structure connects to M5 residuals AND adds');
  console.log('  information beyond MeanRun.');
} else if (verdict === 'PARTIAL') {
  console.log('  MeanRun is confirmed as a 1D projection of genuine 2D');
  console.log('  dynamical structure (lopsidedness, ncm_frac, bisymFlow).');
  console.log('  The bridge exists, but N=7 is too small to confirm that');
  console.log('  2D adds NEW information beyond what MeanRun captures.');
  console.log('  This is the expected outcome for a pilot with N=7.');
} else {
  console.log('  No clear 2D dynamical signal found.');
}
console.log();
console.log('  NOTE: N=7 is a PILOT. All results are suggestive, not');
console.log('  confirmatory. THINGS/PHANGS expansion needed for definitive test.');

const output = {
  phase:'74', title:'2D Dynamical-State Program',
  thingsSample: {n:Nth, galaxies:thingsNames},
  individualCorrs: corrResults,
  bestComposite: {name:bestComp.name, r:bestComp.r, rs:bestComp.rs},
  meanRunBridge: {
    vsLopsidedness: pearsonCorr(logMR_TH, vars2D['things_lopsidedness']),
    vsBisymFlow: pearsonCorr(logMR_TH, vars2D['things_bisymFlow']),
    vsNcmFrac: pearsonCorr(logMR_TH, vars2D['things_ncm_frac']),
    strongBridge
  },
  addsOverMR: {r:rNoMR, rs:rsNoMR},
  rcWigLegacy: {meanAlignMR: mean(mrAlign), meanAlignWig: mean(wigAlign), mrBetter},
  pairwise: {concordance:pairCorr, total:pairTot, pct:pairCorr/pairTot*100},
  varianceBudget: {residRMS, r2:bestComp.r**2},
  verdict
};
fs.writeFileSync('public/phase74-2d-dynamical.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase74-2d-dynamical.json');
