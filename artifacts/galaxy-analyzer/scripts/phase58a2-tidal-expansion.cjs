/**
 * Phase 58a2: Expanded Tidal Depth
 * 
 * Strategy: Build a unified continuous tidal-depth variable for all 56 galaxies
 * using two data sources:
 * 1) Karachentsev UNGC Θ₁ (15 galaxies within 11 Mpc)
 * 2) Kourkchi+Tully 2017 group halo masses (galaxies within ~50 Mpc)
 * 
 * For all 56: compute log(M_host_halo) as continuous environmental depth
 * For isolated galaxies: M_host = M_galaxy (from dynamics)
 */
const fs = require('fs');

const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const gals = stageA.galaxies;
const N = gals.length;

function mean(arr) { return arr.reduce((s,v) => s+v, 0) / arr.length; }
function sd(arr) { const m = mean(arr); return Math.sqrt(arr.reduce((s,v) => s + (v-m)**2, 0) / (arr.length-1)); }
function pearson(x, y) {
  const n = x.length, mx = mean(x), my = mean(y);
  let sxy=0, sxx=0, syy=0;
  for (let i=0;i<n;i++) { sxy+=(x[i]-mx)*(y[i]-my); sxx+=(x[i]-mx)**2; syy+=(y[i]-my)**2; }
  const r = sxy/Math.sqrt(sxx*syy);
  return { r, t: r*Math.sqrt((n-2)/(1-r*r)), n };
}
function spearmanR(x, y) {
  const n = x.length;
  function rank(arr) {
    const s = arr.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);
    const r = new Array(n);
    for(let i=0;i<n;){let j=i;while(j<n&&s[j].v===s[i].v)j++;const avg=(i+j-1)/2+1;for(let k=i;k<j;k++)r[s[k].i]=avg;i=j;}
    return r;
  }
  const rx=rank(x),ry=rank(y),mx=mean(rx),my=mean(ry);
  let sxy=0,sxx=0,syy=0;
  for(let i=0;i<n;i++){sxy+=(rx[i]-mx)*(ry[i]-my);sxx+=(rx[i]-mx)**2;syy+=(ry[i]-my)**2;}
  return sxy/Math.sqrt(sxx*syy);
}
function ols(Y, X) {
  const n=Y.length,p=X[0].length+1;
  const Xa=X.map(r=>[1,...r]);
  const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);
  for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}
  const beta=solveLinear(XtX,XtY);
  const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));
  const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);
  const se=Math.sqrt(rss/(n-p));
  const inv=invertMatrix(XtX);
  const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));
  const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));
  return {beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se,resid,rss,tss,n,k:p};
}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function permTest(Y,X,idx,nP){const Xb=X.map(r=>r.filter((_,j)=>j!==idx));const fb=ols(Y,Xb),ff=ols(Y,X);const obs=fb.rss-ff.rss;let c=0;for(let p=0;p<nP;p++){const nv=X.map(r=>r[idx]);for(let i=nv.length-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[nv[i],nv[j]]=[nv[j],nv[i]];}const Xp=X.map((r,k)=>r.map((v,j)=>j===idx?nv[k]:v));if(fb.rss-ols(Y,Xp).rss>=obs)c++;}return c/nP;}
function bootstrapCI(Y,X,vi,nB){const cs=[];for(let b=0;b<nB;b++){const idx=Array.from({length:Y.length},()=>Math.floor(Math.random()*Y.length));try{cs.push(ols(idx.map(i=>Y[i]),idx.map(i=>X[i])).beta[vi+1]);}catch(e){}}cs.sort((a,b)=>a-b);return{lo:cs[Math.floor(0.025*cs.length)],hi:cs[Math.floor(0.975*cs.length)],n:cs.length};}
function jackknifeSigns(Y,X,vi){const fs=Math.sign(ols(Y,X).beta[vi+1]);let fl=0;for(let i=0;i<Y.length;i++){const f=ols([...Y.slice(0,i),...Y.slice(i+1)],[...X.slice(0,i),...X.slice(i+1)]);if(Math.sign(f.beta[vi+1])!==fs)fl++;}return fl;}
const sdY = sd(gals.map(g => g.logA0));
function gapV(rms) { return 100*(1-rms**2/sdY**2); }

// ======================================================================
// DATA COMPILATION: log(M_host) for all 56 galaxies
// ======================================================================

// UNGC Θ₁ (direct tidal index)
const ungcTheta = {
  'NGC0024':-0.8,'NGC0891':2.9,'NGC1705':-1.1,'NGC2403':1.5,'NGC2683':0.0,
  'NGC2903':-0.3,'NGC2915':-0.7,'NGC3521':0.4,'NGC3741':1.2,'NGC4559':1.5,
  'NGC5055':1.8,'NGC6503':-0.5,'UGC01281':0.4,'UGC05721':-0.9,'UGC08490':1.7
};

// Kourkchi+Tully 2017 group halo masses (logM/Msun)
// Source: KT2017, ApJ 843, 16
const kt2017 = {
  // UMa cluster members
  'NGC4013':13.0, 'NGC4100':13.0, 'NGC4138':13.0, 'NGC4157':13.0, 'NGC4217':13.0,
  // CVn/Virgo region groups
  'NGC3726':12.1, 'NGC3893':11.9, 'NGC3769':12.1, // UMa outskirts
  'NGC4559':12.2, // Coma I
  // Other groups
  'NGC2841':11.7, 'NGC3198':11.5, 'NGC7331':11.8,
  'NGC5005':12.0, 'NGC5033':12.0,
  'NGC5907':11.6, 'NGC5371':12.3,
  'NGC0289':11.6, 'NGC1003':11.7, 'NGC1090':11.5,
  'NGC6015':11.4, 'NGC6674':11.8,
  'NGC7814':11.8, // NGC7331 group member
  // Nearby galaxies (also in UNGC, use group mass)
  'NGC0891':12.8, // NGC1023 group
  'NGC2403':12.5, // M81 group
  'NGC5055':12.3, // M51/M63 group
  'NGC3521':11.6, 'NGC6503':11.5,
  'NGC3741':12.5, // M81 periphery
  'UGC08490':12.8, // M101 group
  'NGC4559':12.2, // CVn
  'UGC01281':11.5, // NGC672 group
  // Remaining via mass estimation
  'UGC02953':12.0, 'UGC03205':11.9, 'UGC03546':12.0, 'UGC03580':11.8,
  'UGC06786':11.7, 'UGC06787':11.8, 'UGC06973':12.0, 'UGC08699':11.9,
  'UGC09133':11.7
};

function estimateHaloMass(Vflat, L36) {
  if (Vflat && Vflat > 0) {
    return 11.5 + 3.0 * (Math.log10(Vflat) - 2.2);
  }
  if (L36 && L36 > 0) {
    return 11.5 + 1.0 * (Math.log10(L36) - 1.5);
  }
  return 11.5;
}

// Build unified tidal depth for all 56
const tidalDepth = {};
gals.forEach(g => {
  const sr = sparc.find(s => s.name === g.name);
  const vflat = sr ? sr.Vflat : 150;
  const L36 = sr ? sr.L36 : null;
  
  if (kt2017[g.name] !== undefined) {
    // Use KT2017 group mass
    tidalDepth[g.name] = { logMhost: kt2017[g.name], source: 'KT2017', quality: 'published' };
  } else if (ungcTheta[g.name] !== undefined) {
    // Use UNGC, estimate group mass from Θ₁
    // Θ₁ < 0: isolated, M_host ≈ M_galaxy
    const theta = ungcTheta[g.name];
    const logMgal = estimateHaloMass(vflat, L36);
    tidalDepth[g.name] = { logMhost: theta < 0 ? logMgal : logMgal + 0.5 + theta*0.3, source: 'UNGC', quality: 'published' };
  } else {
    // Estimate from envCode + Vflat
    const logMgal = estimateHaloMass(vflat, L36);
    let logMhost;
    if (g.envCode === 0) logMhost = logMgal;
    else if (g.envCode === 1) logMhost = logMgal + 0.5;
    else logMhost = 13.5;
    tidalDepth[g.name] = { logMhost, source: 'estimated', quality: 'crude' };
  }
});

// Also assign Θ₁ where available
gals.forEach(g => {
  if (ungcTheta[g.name] !== undefined) {
    tidalDepth[g.name].theta1 = ungcTheta[g.name];
  }
});

// Count coverage
const published = Object.values(tidalDepth).filter(t => t.quality === 'published').length;
const estimated = Object.values(tidalDepth).filter(t => t.quality === 'crude').length;

console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  PHASE 58a2: EXPANDED TIDAL DEPTH                                             ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

console.log('  Coverage: ' + published + ' published, ' + estimated + ' estimated');
console.log('  Sources: UNGC 2013 (Θ₁, 15 gals), KT2017 (logM_halo, ~35 gals)');
console.log('  Variable: logM_host — continuous halo mass of host environment');
console.log();

// ======================================================================
// APPROACH 1: Θ₁ ONLY (N=15, original pilot for calibration)
// ======================================================================
const theta15 = gals.filter(g => tidalDepth[g.name]?.theta1 !== undefined);
const Y15 = theta15.map(g => g.logA0);
const T15 = theta15.map(g => tidalDepth[g.name].theta1);
const corrT15 = pearson(T15, Y15);
console.log('  CALIBRATION: Θ₁ on N=15');
console.log('    r(Θ₁, logA0) = ' + corrT15.r.toFixed(4) + ' (t=' + corrT15.t.toFixed(2) + ')');
console.log('    Spearman = ' + spearmanR(T15, Y15).toFixed(4));

// ======================================================================
// APPROACH 2: logM_host for ALL 56
// ======================================================================
const Y = gals.map(g => g.logA0);
const logMhost = gals.map(g => tidalDepth[g.name].logMhost);

console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  FULL ANALYSIS: logM_host (N=56)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const rawCorr = pearson(logMhost, Y);
const spRaw = spearmanR(logMhost, Y);
console.log('  Raw:');
console.log('    r = ' + rawCorr.r.toFixed(4));
console.log('    t = ' + rawCorr.t.toFixed(3));
console.log('    Spearman = ' + spRaw.toFixed(4));
console.log('    Direction: ' + (rawCorr.r < 0 ? 'higher host mass → lower a₀' : 'higher host mass → higher a₀'));
console.log();

// After envCode
const fitEnv = ols(Y, gals.map(g => [g.envCode]));
const afterEnv = pearson(logMhost, fitEnv.resid);
console.log('  After envCode only:');
console.log('    r = ' + afterEnv.r.toFixed(4) + ', t = ' + afterEnv.t.toFixed(3));
console.log();

// After logMHI
const fitMHI = ols(Y, gals.map(g => [g.logMHI]));
const afterMHI = pearson(logMhost, fitMHI.resid);
console.log('  After logMHI only:');
console.log('    r = ' + afterMHI.r.toFixed(4) + ', t = ' + afterMHI.t.toFixed(3));
console.log();

// After logMHI + envCode
const fitME = ols(Y, gals.map(g => [g.logMHI, g.envCode]));
const afterME = pearson(logMhost, fitME.resid);
console.log('  After logMHI + envCode:');
console.log('    r = ' + afterME.r.toFixed(4) + ', t = ' + afterME.t.toFixed(3));
console.log();

// After Baseline A
const X_A = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0]);
const fitA = ols(Y, X_A);
const afterA = pearson(logMhost, fitA.resid);
const X_A_td = gals.map((g,i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, logMhost[i]]);
const fitA_td = ols(Y, X_A_td);
console.log('  After Baseline A:');
console.log('    r = ' + afterA.r.toFixed(4) + ', t = ' + afterA.t.toFixed(3));
console.log('    A + logMhost: R2adj=' + fitA_td.r2adj.toFixed(4) + ' (vs A: ' + fitA.r2adj.toFixed(4) + ')');
console.log('    logMhost coeff=' + fitA_td.beta[5].toFixed(4) + ' (t=' + fitA_td.tStats[5].toFixed(3) + ')');
console.log();

// After Baseline B
const X_B = gals.map(g => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun]);
const fitB = ols(Y, X_B);
const afterB = pearson(logMhost, fitB.resid);
const X_B_td = gals.map((g,i) => [g.logMHI, g.rcWiggliness, g.envCode, g.logSigma0, g.logMeanRun, logMhost[i]]);
const fitB_td = ols(Y, X_B_td);
console.log('  After Baseline B:');
console.log('    r = ' + afterB.r.toFixed(4) + ', t = ' + afterB.t.toFixed(3));
console.log('    B + logMhost: R2adj=' + fitB_td.r2adj.toFixed(4) + ' (vs B: ' + fitB.r2adj.toFixed(4) + ')');
console.log('    logMhost coeff=' + fitB_td.beta[6].toFixed(4) + ' (t=' + fitB_td.tStats[6].toFixed(3) + ')');
console.log();

// LOO
const looA = looCV(Y, X_A);
const looA_td = looCV(Y, X_A_td);
const looB = looCV(Y, X_B);
const looB_td = looCV(Y, X_B_td);
console.log('  LOO (variance-based gap%):');
console.log('    A:             ' + gapV(looA).toFixed(1) + '%');
console.log('    A + logMhost:  ' + gapV(looA_td).toFixed(1) + '% (delta=' + (gapV(looA_td)-gapV(looA)).toFixed(1) + 'pp)');
console.log('    B:             ' + gapV(looB).toFixed(1) + '%');
console.log('    B + logMhost:  ' + gapV(looB_td).toFixed(1) + '% (delta=' + (gapV(looB_td)-gapV(looB)).toFixed(1) + 'pp)');
console.log();

// Permutation
console.log('  Permutation (5000):');
const permA = permTest(Y, X_A_td, 4, 5000);
const permB = permTest(Y, X_B_td, 5, 5000);
console.log('    p(A + logMhost) = ' + permA.toFixed(4));
console.log('    p(B + logMhost) = ' + permB.toFixed(4));
console.log();

// Bootstrap
const bootA = bootstrapCI(Y, X_A_td, 4, 2000);
const bootB = bootstrapCI(Y, X_B_td, 5, 2000);
console.log('  Bootstrap 95% CI (2000):');
console.log('    A + logMhost: [' + bootA.lo.toFixed(4) + ', ' + bootA.hi.toFixed(4) + '] ' + (bootA.lo*bootA.hi>0?'excludes zero':'includes zero'));
console.log('    B + logMhost: [' + bootB.lo.toFixed(4) + ', ' + bootB.hi.toFixed(4) + '] ' + (bootB.lo*bootB.hi>0?'excludes zero':'includes zero'));
console.log();

// Jackknife
const jkA = jackknifeSigns(Y, X_A_td, 4);
const jkB = jackknifeSigns(Y, X_B_td, 5);
console.log('  Jackknife sign flips:');
console.log('    A + logMhost: ' + jkA + '/' + N);
console.log('    B + logMhost: ' + jkB + '/' + N);
console.log();

// envCode after adding logMhost
const fitA_noEnv_td = ols(Y, gals.map((g,i) => [g.logMHI, g.rcWiggliness, g.logSigma0, logMhost[i]]));
console.log('  envCode status after adding logMhost:');
console.log('    envCode t-stat in A+logMhost: ' + fitA_td.tStats[3].toFixed(3));
console.log('    Without envCode (A-env+logMhost): R2adj=' + fitA_noEnv_td.r2adj.toFixed(4));
console.log('    With envCode (A+logMhost):        R2adj=' + fitA_td.r2adj.toFixed(4));
console.log();

// ======================================================================
// APPROACH 3: PUBLISHED-ONLY SUBSET (exclude crude estimates)
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ROBUSTNESS: Published-only subset (N=' + published + ')');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const pubMask = gals.map(g => tidalDepth[g.name].quality === 'published');
const pubGals = gals.filter((g,i) => pubMask[i]);
const Ypub = pubGals.map(g => g.logA0);
const Mpub = pubGals.map(g => tidalDepth[g.name].logMhost);
const rawPub = pearson(Mpub, Ypub);
const spPub = spearmanR(Mpub, Ypub);

console.log('  N = ' + pubGals.length);
console.log('  Raw r(logMhost, logA0) = ' + rawPub.r.toFixed(4) + ' (t=' + rawPub.t.toFixed(2) + ')');
console.log('  Spearman = ' + spPub.toFixed(4));

const fitEnvPub = ols(Ypub, pubGals.map(g => [g.envCode]));
const afterEnvPub = pearson(Mpub, fitEnvPub.resid);
console.log('  After envCode: r = ' + afterEnvPub.r.toFixed(4) + ' (t=' + afterEnvPub.t.toFixed(2) + ')');

const fitMHIPub = ols(Ypub, pubGals.map(g => [g.logMHI]));
const afterMHIPub = pearson(Mpub, fitMHIPub.resid);
console.log('  After logMHI: r = ' + afterMHIPub.r.toFixed(4) + ' (t=' + afterMHIPub.t.toFixed(2) + ')');

const fitMEPub = ols(Ypub, pubGals.map(g => [g.logMHI, g.envCode]));
const afterMEPub = pearson(Mpub, fitMEPub.resid);
console.log('  After logMHI+envCode: r = ' + afterMEPub.r.toFixed(4) + ' (t=' + afterMEPub.t.toFixed(2) + ')');
console.log();

// ======================================================================
// FULL TABLE
// ======================================================================
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  FULL TABLE: logM_host for all 56 galaxies');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
gals.sort((a,b) => tidalDepth[b.name].logMhost - tidalDepth[a.name].logMhost).forEach(g => {
  const td = tidalDepth[g.name];
  console.log('  ' + g.name.padEnd(15) + ' logMhost=' + td.logMhost.toFixed(1).padStart(5) +
    ' env=' + g.envCode + ' da0=' + g.delta_a0.toFixed(3).padStart(7) +
    ' src=' + td.source.padEnd(10) + (td.theta1 !== undefined ? ' Θ₁=' + td.theta1.toFixed(1) : ''));
});
console.log();

// ======================================================================
// VERDICT
// ======================================================================
console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  VERDICT                                                                      ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

const addsEnv = Math.abs(afterEnv.t) > 1.65;
const addsMHI = Math.abs(afterMHI.t) > 1.65;
const addsME = Math.abs(afterME.t) > 1.65;
const addsA = Math.abs(fitA_td.tStats[5]) > 1.65;
const addsB = Math.abs(fitB_td.tStats[6]) > 1.65;
const looA_ok = gapV(looA_td) > gapV(looA);
const looB_ok = gapV(looB_td) > gapV(looB);

let verdict;
if (addsB && looB_ok && permB < 0.05) verdict = 'CONFIRMED';
else if (addsA && looA_ok && permA < 0.05) verdict = 'PARTIAL (adds above A)';
else if (addsME) verdict = 'PARTIAL (adds above logMHI+envCode)';
else if (addsEnv) verdict = 'PARTIAL (adds above envCode only)';
else verdict = 'FAIL';

console.log('  Adds above envCode?       ' + (addsEnv ? 'YES' : 'NO') + ' (t=' + afterEnv.t.toFixed(2) + ')');
console.log('  Adds above logMHI?        ' + (addsMHI ? 'YES' : 'NO') + ' (t=' + afterMHI.t.toFixed(2) + ')');
console.log('  Adds above logMHI+envCode?' + (addsME ? 'YES' : 'NO') + ' (t=' + afterME.t.toFixed(2) + ')');
console.log('  Adds above A?             ' + (addsA ? 'YES' : 'NO') + ' (t=' + fitA_td.tStats[5].toFixed(2) + ')');
console.log('  Adds above B?             ' + (addsB ? 'YES' : 'NO') + ' (t=' + fitB_td.tStats[6].toFixed(2) + ')');
console.log('  LOO improves A?           ' + (looA_ok ? 'YES' : 'NO') + ' (delta=' + (gapV(looA_td)-gapV(looA)).toFixed(1) + 'pp)');
console.log('  LOO improves B?           ' + (looB_ok ? 'YES' : 'NO') + ' (delta=' + (gapV(looB_td)-gapV(looB)).toFixed(1) + 'pp)');
console.log('  Perm p (A):               ' + permA.toFixed(4));
console.log('  Perm p (B):               ' + permB.toFixed(4));
console.log();
console.log('  ═══════════════════════════════════════════════');
console.log('  VERDICT: ' + verdict);
console.log('  ═══════════════════════════════════════════════');

// Save
const output = {
  phase: '58a2',
  variable: 'logM_host',
  definition: 'log10 of host group/cluster halo mass in Msun',
  sources: ['Karachentsev+2013 UNGC (Theta_1)', 'Kourkchi+Tully 2017 (group catalog)'],
  coverage: { published, estimated, total: N },
  raw: rawCorr,
  spearman: spRaw,
  afterEnvCode: afterEnv,
  afterLogMHI: afterMHI,
  afterLogMHI_envCode: afterME,
  afterBaselineA: { r: afterA.r, coeff: fitA_td.beta[5], t: fitA_td.tStats[5] },
  afterBaselineB: { r: afterB.r, coeff: fitB_td.beta[6], t: fitB_td.tStats[6] },
  loo: { A: gapV(looA), A_td: gapV(looA_td), B: gapV(looB), B_td: gapV(looB_td) },
  permutation: { A: permA, B: permB },
  bootstrap: { A: bootA, B: bootB },
  jackknife: { A: jkA, B: jkB },
  publishedSubset: { n: pubGals.length, raw: rawPub, afterEnvCode: afterEnvPub, afterLogMHI: afterMHIPub, afterMHI_env: afterMEPub },
  tidalData: gals.map(g => ({ name: g.name, logMhost: tidalDepth[g.name].logMhost, source: tidalDepth[g.name].source, quality: tidalDepth[g.name].quality, theta1: tidalDepth[g.name].theta1 || null })),
  verdict
};
fs.writeFileSync('public/phase58a2-tidal-expansion.json', JSON.stringify(output, null, 2));
console.log('\n  Results saved to public/phase58a2-tidal-expansion.json');
