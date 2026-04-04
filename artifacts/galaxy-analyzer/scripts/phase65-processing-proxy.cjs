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
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se:Math.sqrt(rss/(n-p)),resid,rss,tss,n,k:p};}
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
const morphT = gals45.map(g => sparcMap[g.name]?.T || 5);
const logUps = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
function gapV(rms){return 100*(1-rms**2/sdY**2);}

const Xconf = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT[i]]);
const fConf = ols(logUps, Xconf);
const upsPerp = fConf.resid;

const X_Bdp = gals45.map((g,i) => [g.logMHI, g.rcWiggliness, logMhost[i], g.logSigma0, g.logMeanRun, upsPerp[i]]);
const fBdp = ols(Y, X_Bdp);
const residBdp = fBdp.resid;

console.log('======================================================================');
console.log('  PHASE 65: PROCESSING / INFALL PROXY HUNT');
console.log('  Residual source: B″ on N=45');
console.log('  Residual SD = ' + sd(residBdp).toFixed(4) + ' dex');
console.log('======================================================================\n');

// ═══════════════════════════════════════════════════════════════════
// Build candidate processing proxies
// ═══════════════════════════════════════════════════════════════════

// 1. HI truncation ratio: log(RHI/Rdisk)
// Stripped galaxies have smaller RHI/Rdisk
const logHIratio = gals45.map(g => {
  const rhi = sparcMap[g.name]?.RHI;
  const rd = sparcMap[g.name]?.Rdisk;
  return rhi && rhi > 0 && rd && rd > 0 ? Math.log10(rhi/rd) : null;
});

// 2. HI deficiency (already available)
const hiDef = gals45.map(g => g.hi_deficiency);

// 3. HI deficiency⊥ — residualized vs B″ confounders (logMHI, logMhost, logΣ₀)
// This removes what B″ already knows
const hiDefValid = hiDef.every(v => v !== null);

// 4. Gas extent normalized: RHI / (some size scale independent of B″ vars)
const logRHInorm = gals45.map(g => {
  const rhi = sparcMap[g.name]?.RHI;
  const D = sparcMap[g.name]?.D;
  return rhi && rhi > 0 && D && D > 0 ? Math.log10(rhi) : null;
});

// 5. Processing score composite: HI_def + log(RHI/Rd) combined
// Higher = more processed (more deficient, more truncated)
const procScore = gals45.map((g,i) => {
  const hd = hiDef[i];
  const lr = logHIratio[i];
  if (hd === null || lr === null) return null;
  // HI_def positive = deficient; logRHI/Rd smaller = more truncated
  // Combine: processing = HI_def - log(RHI/Rd) (both pointing same way when processed)
  return hd - lr;
});

// 6. Morphological disturbance proxy: |T - T_expected_from_luminosity|
// Disturbed galaxies might have unusual morphology for their luminosity

// 7. Ram pressure susceptibility: ~MHI / (Vflat² * Rdisk)
const logRPsusc = gals45.map(g => {
  const mhi = sparcMap[g.name]?.MHI;
  const vf = sparcMap[g.name]?.Vflat;
  const rd = sparcMap[g.name]?.Rdisk;
  return mhi && mhi > 0 && vf && vf > 0 && rd && rd > 0 
    ? Math.log10(mhi / (vf*vf * rd)) : null;
});

// Collect all candidates
const candidates = [
  {name: 'HI_deficiency', desc: 'HI deficiency (obs-expected)', vals: hiDef, source: 'Stage-A'},
  {name: 'logRHI/Rd', desc: 'log(RHI/Rdisk) — gas truncation', vals: logHIratio, source: 'SPARC'},
  {name: 'procScore', desc: 'HI_def - log(RHI/Rd) — processing composite', vals: procScore, source: 'derived'},
  {name: 'logRPsusc', desc: 'log(MHI/(Vflat²·Rd)) — RP susceptibility', vals: logRPsusc, source: 'derived'},
];

function pearsonCorr(x, y) {
  const n = x.length;
  const mx=mean(x),my=mean(y);
  let sxy=0,sxx=0,syy=0;
  for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;syy+=(y[i]-my)**2;}
  if(sxx<1e-15||syy<1e-15)return{r:NaN,t:NaN};
  const r=sxy/Math.sqrt(sxx*syy);
  return{r,t:r*Math.sqrt((n-2)/(1-r*r))};
}

// ═══════════════════════════════════════════════════════════════════
// TEST EACH CANDIDATE
// ═══════════════════════════════════════════════════════════════════

const results = [];

candidates.forEach(c => {
  const validIdx = [];
  for (let i = 0; i < N; i++) {
    if (c.vals[i] !== null && !isNaN(c.vals[i])) validIdx.push(i);
  }
  if (validIdx.length < 20) return;
  
  const residV = validIdx.map(i => residBdp[i]);
  const valsV = validIdx.map(i => c.vals[i]);
  const nV = validIdx.length;
  
  console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
  console.log('  CANDIDATE: ' + c.name);
  console.log('  Definition: ' + c.desc);
  console.log('  Source: ' + c.source);
  console.log('  N = ' + nV);
  console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
  
  // 1. Raw correlation with B″ residuals
  const rawCorr = pearsonCorr(residV, valsV);
  console.log('  Raw with resid(B″):');
  console.log('    r = ' + rawCorr.r.toFixed(4) + ', t = ' + rawCorr.t.toFixed(3));
  console.log();
  
  // 2. After controlling for distance
  const distV = validIdx.map(i => Math.log10(sparcMap[gals45[i].name]?.D || 10));
  const fDist = ols(valsV, distV.map(d => [d]));
  const valsDistResid = fDist.resid;
  const distCorr = pearsonCorr(residV, valsDistResid);
  console.log('  After distance control:');
  console.log('    r = ' + distCorr.r.toFixed(4) + ', t = ' + distCorr.t.toFixed(3));
  console.log();
  
  // 3. After controlling for logMhost
  const mhostV = validIdx.map(i => logMhost[i]);
  const fMhost = ols(valsV, mhostV.map(m => [m]));
  const valsMhostResid = fMhost.resid;
  const mhostCorr = pearsonCorr(residV, valsMhostResid);
  console.log('  After logMhost control:');
  console.log('    r = ' + mhostCorr.r.toFixed(4) + ', t = ' + mhostCorr.t.toFixed(3));
  console.log();
  
  // 4. Full orthogonalization: residualize vs ALL B″ predictors
  const XbdpV = validIdx.map(i => X_Bdp[i]);
  const fOrth = ols(valsV, XbdpV);
  const valsPerp = fOrth.resid;
  const perpCorr = pearsonCorr(residV, valsPerp);
  console.log('  After full B″ orthogonalization:');
  console.log('    r = ' + perpCorr.r.toFixed(4) + ', t = ' + perpCorr.t.toFixed(3));
  console.log('    B″ explains ' + ((1-fOrth.rss/fOrth.tss)*100).toFixed(1) + '% of ' + c.name);
  console.log();
  
  // 5. LOO test: B″ vs B″ + candidate
  if (nV === N) {
    const X_aug = gals45.map((g,i) => [...X_Bdp[i], c.vals[i]]);
    const looBase = looCV(Y, X_Bdp);
    const looAug = looCV(Y, X_aug);
    const fAug = ols(Y, X_aug);
    const tNew = fAug.tStats[7];
    
    console.log('  LOO test (full N=' + N + '):');
    console.log('    B″ LOO gap% = ' + gapV(looBase).toFixed(1) + '%');
    console.log('    B″ + ' + c.name + ' LOO gap% = ' + gapV(looAug).toFixed(1) + '%');
    console.log('    Δgap = ' + (gapV(looAug)-gapV(looBase)).toFixed(1) + 'pp');
    console.log('    t-stat in augmented model = ' + tNew.toFixed(3));
    console.log();
    
    // Permutation test
    let permWins = 0;
    const nPerm = 5000;
    for (let p = 0; p < nPerm; p++) {
      const pVals = [...c.vals];
      for (let i = N-1; i > 0; i--) {
        const j = Math.floor(Math.random()*(i+1));
        [pVals[i], pVals[j]] = [pVals[j], pVals[i]];
      }
      const Xp = gals45.map((g,i) => [...X_Bdp[i], pVals[i]]);
      try {
        const fp = ols(Y, Xp);
        if (Math.abs(fp.tStats[7]) >= Math.abs(tNew)) permWins++;
      } catch(e) {}
    }
    console.log('    Permutation p = ' + (permWins/nPerm).toFixed(4) + ' (' + nPerm + ' perms)');
    
    // Bootstrap CI
    const bootBetas = [];
    for (let b = 0; b < 2000; b++) {
      const idx = Array.from({length:N}, () => Math.floor(Math.random()*N));
      try {
        const fb = ols(idx.map(i=>Y[i]), idx.map(i=>X_aug[i]));
        bootBetas.push(fb.beta[7]);
      } catch(e) {}
    }
    bootBetas.sort((a,b) => a-b);
    const ci95lo = bootBetas[Math.floor(0.025*bootBetas.length)];
    const ci95hi = bootBetas[Math.floor(0.975*bootBetas.length)];
    const exclZero = ci95lo > 0 || ci95hi < 0;
    console.log('    Bootstrap 95% CI = [' + ci95lo.toFixed(4) + ', ' + ci95hi.toFixed(4) + ']');
    console.log('    Excludes zero? ' + (exclZero ? 'YES' : 'NO'));
    
    // Jackknife sign flips
    let signFlips = 0;
    const fullSign = Math.sign(fAug.beta[7]);
    for (let i = 0; i < N; i++) {
      const Yj = [...Y.slice(0,i),...Y.slice(i+1)];
      const Xj = [...X_aug.slice(0,i),...X_aug.slice(i+1)];
      try {
        const fj = ols(Yj, Xj);
        if (Math.sign(fj.beta[7]) !== fullSign) signFlips++;
      } catch(e) {}
    }
    console.log('    Jackknife sign flips = ' + signFlips + '/' + N);
    console.log();
    
    // Summary
    const pass = [
      gapV(looAug) > gapV(looBase),
      Math.abs(tNew) > 1.65,
      permWins/nPerm < 0.05,
      exclZero,
      signFlips <= 2
    ];
    const nPass = pass.filter(Boolean).length;
    
    console.log('  Robustness: ' + nPass + '/5');
    console.log('    LOO improves:     ' + (pass[0]?'✅':'❌'));
    console.log('    |t| > 1.65:       ' + (pass[1]?'✅':'❌'));
    console.log('    perm p < 0.05:    ' + (pass[2]?'✅':'❌'));
    console.log('    CI excludes 0:    ' + (pass[3]?'✅':'❌'));
    console.log('    Jackknife stable: ' + (pass[4]?'✅':'❌'));
    console.log();
    
    let verdict;
    if (nPass >= 4) verdict = 'CONFIRMED-AS-AXIS-7-CANDIDATE';
    else if (nPass >= 3) verdict = 'PARTIAL';
    else verdict = 'FAIL';
    
    console.log('  Verdict: ' + verdict);
    console.log();
    
    results.push({
      name: c.name, desc: c.desc, n: nV,
      rawR: rawCorr.r, rawT: rawCorr.t,
      mhostCtrlR: mhostCorr.r, mhostCtrlT: mhostCorr.t,
      perpR: perpCorr.r, perpT: perpCorr.t,
      looBase: gapV(looBase), looAug: gapV(looAug), looDelta: gapV(looAug)-gapV(looBase),
      tInModel: tNew,
      permP: permWins/nPerm,
      bootCI: [ci95lo, ci95hi], exclZero,
      jackFlips: signFlips,
      robustness: nPass, verdict
    });
  } else {
    results.push({
      name: c.name, desc: c.desc, n: nV,
      rawR: rawCorr.r, rawT: rawCorr.t,
      perpR: perpCorr.r, perpT: perpCorr.t,
      verdict: Math.abs(perpCorr.t) > 1.5 ? 'WORTH_EXPANDING' : 'FAIL'
    });
  }
});

// ═══════════════════════════════════════════════════════════════════
// Also test orthogonalized versions
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ORTHOGONALIZED PROCESSING PROXIES');
console.log('  (residualized vs B″ confounders before testing)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// HI_def⊥: residualize HI_def vs (logMHI, logMhost, logΣ₀, morphT)
const Xhiconf = gals45.map((g,i) => [g.logMHI, logMhost[i], g.logSigma0, morphT[i]]);
const fHiConf = ols(hiDef, Xhiconf);
const hiDefPerp = fHiConf.resid;
console.log('  HI_deficiency⊥:');
console.log('    Confounders explain ' + ((1-fHiConf.rss/fHiConf.tss)*100).toFixed(1) + '% of HI_def');

const hiPerpCorr = pearsonCorr(residBdp, hiDefPerp);
console.log('    r(resid_B″, HI_def⊥) = ' + hiPerpCorr.r.toFixed(4) + ', t = ' + hiPerpCorr.t.toFixed(3));

if (Math.abs(hiPerpCorr.t) > 1.0) {
  const X_hiAug = gals45.map((g,i) => [...X_Bdp[i], hiDefPerp[i]]);
  const looHi = looCV(Y, X_hiAug);
  const looBase = looCV(Y, X_Bdp);
  const fHiAug = ols(Y, X_hiAug);
  console.log('    LOO: B″ = ' + gapV(looBase).toFixed(1) + '%, B″+HI_def⊥ = ' + gapV(looHi).toFixed(1) + '%, Δ = ' + (gapV(looHi)-gapV(looBase)).toFixed(1) + 'pp');
  console.log('    t-stat = ' + fHiAug.tStats[7].toFixed(3));
}
console.log();

// logRHI/Rd⊥
const logHIratioClean = gals45.map((g,i) => {
  const v = logHIratio[i];
  return v !== null ? v : 0;
});
const hasRHI = logHIratio.filter(v => v !== null).length;

if (hasRHI >= 40) {
  const validHR = gals45.map((_,i) => i).filter(i => logHIratio[i] !== null);
  const Xhrconf = validHR.map(i => [gals45[i].logMHI, logMhost[i], gals45[i].logSigma0]);
  const Yhr = validHR.map(i => logHIratio[i]);
  const fHrConf = ols(Yhr, Xhrconf);
  const hrPerp = fHrConf.resid;
  
  const residHR = validHR.map(i => residBdp[i]);
  const hrPerpCorr = pearsonCorr(residHR, hrPerp);
  console.log('  logRHI/Rd⊥ (N=' + validHR.length + '):');
  console.log('    r(resid_B″, logRHI/Rd⊥) = ' + hrPerpCorr.r.toFixed(4) + ', t = ' + hrPerpCorr.t.toFixed(3));
  console.log();
}

// procScore⊥
const validPS = gals45.map((_,i) => i).filter(i => procScore[i] !== null);
if (validPS.length >= 40) {
  const XpsConf = validPS.map(i => [gals45[i].logMHI, logMhost[i], gals45[i].logSigma0]);
  const Yps = validPS.map(i => procScore[i]);
  const fPsConf = ols(Yps, XpsConf);
  const psPerp = fPsConf.resid;
  
  const residPS = validPS.map(i => residBdp[i]);
  const psPerpCorr = pearsonCorr(residPS, psPerp);
  console.log('  procScore⊥ (N=' + validPS.length + '):');
  console.log('    r(resid_B″, procScore⊥) = ' + psPerpCorr.r.toFixed(4) + ', t = ' + psPerpCorr.t.toFixed(3));
  console.log();
}

// ═══════════════════════════════════════════════════════════════════
// FINAL VERDICT
// ═══════════════════════════════════════════════════════════════════
console.log('======================================================================');
console.log('  PHASE 65: FINAL VERDICT');
console.log('======================================================================\n');

const bestResult = results.sort((a,b) => (b.robustness||0) - (a.robustness||0))[0];
const anyConfirmed = results.some(r => r.verdict === 'CONFIRMED-AS-AXIS-7-CANDIDATE');
const anyPartial = results.some(r => r.verdict === 'PARTIAL');

console.log('  ┌──────────────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Candidate        │ raw t    │ ⊥ t      │ LOO Δpp  │ Verdict  │');
console.log('  ├──────────────────┼──────────┼──────────┼──────────┼──────────┤');
results.forEach(r => {
  const lt = r.looDelta !== undefined ? r.looDelta.toFixed(1).padStart(5) : '  —  ';
  console.log('  │ ' + r.name.padEnd(16) + ' │ ' + r.rawT.toFixed(2).padStart(6) + '   │ ' + (r.perpT||0).toFixed(2).padStart(6) + '   │ ' + lt + '   │ ' + (r.verdict||'').padEnd(8) + ' │');
});
console.log('  └──────────────────┴──────────┴──────────┴──────────┴──────────┘');
console.log();

if (anyConfirmed) {
  const best = results.find(r => r.verdict === 'CONFIRMED-AS-AXIS-7-CANDIDATE');
  console.log('  ═══════════════════════════════════════════════════════════');
  console.log('  AXIS 7 CANDIDATE: ' + best.name);
  console.log('  ═══════════════════════════════════════════════════════════');
} else if (anyPartial) {
  const best = results.find(r => r.verdict === 'PARTIAL');
  console.log('  PARTIAL SIGNAL: ' + best.name + ' (' + best.robustness + '/5)');
  console.log('  Worth monitoring with better data.');
} else {
  console.log('  ═══════════════════════════════════════════════════════════');
  console.log('  NO PROCESSING PROXY CAPTURES RESIDUAL STRUCTURE');
  console.log('  ═══════════════════════════════════════════════════════════');
  console.log();
  console.log('  The ~0.10 dex intrinsic scatter after B″ is NOT explained');
  console.log('  by any available processing/infall proxy.');
  console.log('  This scatter may require:');
  console.log('    • Resolved 2D kinematic data (tilted-ring models)');
  console.log('    • Direct infall/orbit measurements');
  console.log('    • Higher-resolution mass models');
}

const output = {
  phase: '65', title: 'Processing / Infall Proxy Hunt',
  residualSource: 'B″ on N=45',
  candidates: results,
  verdict: anyConfirmed ? 'AXIS_7_CANDIDATE' : (anyPartial ? 'PARTIAL' : 'NO_PROCESSING_SIGNAL'),
  conclusion: anyConfirmed
    ? 'Processing proxy captures additional a₀ variation'
    : 'No available processing proxy explains B″ residual structure'
};
fs.writeFileSync('public/phase65-processing-proxy.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase65-processing-proxy.json');
