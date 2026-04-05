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
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);return{beta,resid,rss,tss,se:Math.sqrt(rss/(n-p)),n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}
function pearsonCorr(x,y){const n=x.length,mx=mean(x),my=mean(y);let sxy=0,sxx=0,syy=0;for(let i=0;i<n;i++){sxy+=(x[i]-mx)*(y[i]-my);sxx+=(x[i]-mx)**2;syy+=(y[i]-my)**2;}if(sxx<1e-15||syy<1e-15)return{r:NaN,t:NaN};const r=sxy/Math.sqrt(sxx*syy);return{r,t:r*Math.sqrt((n-2)/(1-r*r))};}

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
const morphT = gals45.map(g => sparcMap[g.name]?.T ?? 5);
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
console.log('  PHASE 66: WEAK-BUNDLE / LATENT RESIDUAL TEST');
console.log('  Residual source: B″ on N=45');
console.log('  Residual SD = ' + sd(residBdp).toFixed(4) + ' dex');
console.log('======================================================================\n');

// ═══════════════════════════════════════════════════════════════════
// Build bundle variables (non-circular, non-B″)
// ═══════════════════════════════════════════════════════════════════

const bundleVars = [];
const bundleNames = [];

// 1. log(Reff/Rdisk) — structural compactness (NOT circular like Vflat²/R)
gals45.forEach((g,i) => {
  const sp = sparcMap[g.name]||{};
  const v = (sp.Reff && sp.Rdisk && sp.Reff > 0 && sp.Rdisk > 0) ? Math.log10(sp.Reff/sp.Rdisk) : null;
  if (!bundleVars[0]) bundleVars[0] = [];
  bundleVars[0][i] = v;
});
bundleNames.push('logReff/Rd');

// 2. log(RHI/Rdisk) — gas extent
gals45.forEach((g,i) => {
  const sp = sparcMap[g.name]||{};
  const v = (sp.RHI && sp.Rdisk && sp.RHI > 0 && sp.Rdisk > 0) ? Math.log10(sp.RHI/sp.Rdisk) : null;
  if (!bundleVars[1]) bundleVars[1] = [];
  bundleVars[1][i] = v;
});
bundleNames.push('logRHI/Rd');

// 3. log(SBeff) — effective surface brightness
gals45.forEach((g,i) => {
  const sp = sparcMap[g.name]||{};
  const v = (sp.SBeff && sp.SBeff > 0) ? Math.log10(sp.SBeff) : null;
  if (!bundleVars[2]) bundleVars[2] = [];
  bundleVars[2][i] = v;
});
bundleNames.push('logSBeff');

// 4. log(SBdisk) — disk surface brightness
gals45.forEach((g,i) => {
  const sp = sparcMap[g.name]||{};
  const v = (sp.SBdisk && sp.SBdisk > 0) ? Math.log10(sp.SBdisk) : null;
  if (!bundleVars[3]) bundleVars[3] = [];
  bundleVars[3][i] = v;
});
bundleNames.push('logSBdisk');

// 5. morphT — morphological type (used in Υ★ confounders but NOT directly in B″)
gals45.forEach((g,i) => {
  if (!bundleVars[4]) bundleVars[4] = [];
  bundleVars[4][i] = morphT[i];
});
bundleNames.push('morphT');

// 6. inclination
gals45.forEach((g,i) => {
  const sp = sparcMap[g.name]||{};
  if (!bundleVars[5]) bundleVars[5] = [];
  bundleVars[5][i] = sp.inc || null;
});
bundleNames.push('inc');

// 7. log(D) — distance
gals45.forEach((g,i) => {
  const sp = sparcMap[g.name]||{};
  const v = (sp.D && sp.D > 0) ? Math.log10(sp.D) : null;
  if (!bundleVars[6]) bundleVars[6] = [];
  bundleVars[6][i] = v;
});
bundleNames.push('logD');

// 8. HI deficiency
gals45.forEach((g,i) => {
  if (!bundleVars[7]) bundleVars[7] = [];
  bundleVars[7][i] = g.hi_deficiency;
});
bundleNames.push('HI_def');

// 9. log(L36) — 3.6µm luminosity
gals45.forEach((g,i) => {
  const sp = sparcMap[g.name]||{};
  const v = (sp.L36 && sp.L36 > 0) ? Math.log10(sp.L36) : null;
  if (!bundleVars[8]) bundleVars[8] = [];
  bundleVars[8][i] = v;
});
bundleNames.push('logL36');

// 10. log(eVflat/Vflat) — fractional velocity error
gals45.forEach((g,i) => {
  const sp = sparcMap[g.name]||{};
  const v = (sp.eVflat && sp.Vflat && sp.Vflat > 0) ? Math.log10(sp.eVflat/sp.Vflat) : null;
  if (!bundleVars[9]) bundleVars[9] = [];
  bundleVars[9][i] = v;
});
bundleNames.push('logFracErr');

// 11. envCode
gals45.forEach((g,i) => {
  if (!bundleVars[10]) bundleVars[10] = [];
  bundleVars[10][i] = g.envCode;
});
bundleNames.push('envCode');

// Filter to variables with complete data
const validBundle = [];
const validNames = [];
for (let j = 0; j < bundleVars.length; j++) {
  const allValid = bundleVars[j].every(v => v !== null && v !== undefined && !isNaN(v));
  const nValid = bundleVars[j].filter(v => v !== null && v !== undefined && !isNaN(v)).length;
  if (nValid >= 42) {
    // Fill missing with mean
    const vals = bundleVars[j].filter(v => v !== null && !isNaN(v));
    const m = mean(vals);
    const filled = bundleVars[j].map(v => (v !== null && !isNaN(v)) ? v : m);
    validBundle.push(filled);
    validNames.push(bundleNames[j] + ' (N=' + nValid + ')');
  }
}

console.log('  Bundle variables (' + validNames.length + ' with N≥42):');
validNames.forEach((n,i) => {
  const rc = pearsonCorr(residBdp, validBundle[i]);
  console.log('    ' + (i+1) + '. ' + n.padEnd(22) + ' r=' + rc.r.toFixed(3).padStart(7) + '  t=' + rc.t.toFixed(2).padStart(6));
});
console.log();

const nBundle = validBundle.length;

// ═══════════════════════════════════════════════════════════════════
// First orthogonalize each bundle var against B″ predictors
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  ORTHOGONALIZED BUNDLE (after removing B″ predictor content)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const orthBundle = [];
for (let j = 0; j < nBundle; j++) {
  const fOrth = ols(validBundle[j], X_Bdp);
  orthBundle.push(fOrth.resid);
  const r2 = 1 - fOrth.rss/fOrth.tss;
  const rc = pearsonCorr(residBdp, fOrth.resid);
  console.log('  ' + validNames[j].padEnd(22) + ': B″ explains ' + (r2*100).toFixed(0) + '%, ⊥ r=' + rc.r.toFixed(3) + ' t=' + rc.t.toFixed(2));
}
console.log();

// ═══════════════════════════════════════════════════════════════════
// 66a — Ridge regression on residual B″
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  66a: RIDGE / OLS BUNDLE ON RESIDUAL B″');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Use orthogonalized bundle to predict residBdp
const X_orth = gals45.map((_,i) => orthBundle.map(ob => ob[i]));
const fOLS = ols(residBdp, X_orth);
const r2_bundle = 1 - fOLS.rss/fOLS.tss;
console.log('  OLS bundle (p=' + nBundle + '):');
console.log('    R² = ' + r2_bundle.toFixed(4) + ' (' + (r2_bundle*100).toFixed(1) + '%)');
console.log('    Adj R² = ' + (1 - (1-r2_bundle)*(N-1)/(N-nBundle-1)).toFixed(4));

// Ridge: manual implementation with lambda search
function ridgeCV(Y, X, lambdas) {
  const n = Y.length, p = X[0].length;
  const results = [];
  for (const lam of lambdas) {
    let cvSS = 0;
    for (let i = 0; i < n; i++) {
      const Yt = [...Y.slice(0,i),...Y.slice(i+1)];
      const Xt = [...X.slice(0,i),...X.slice(i+1)];
      const nt = n-1;
      // Ridge: (X'X + λI)^{-1} X'Y
      const XtX = Array.from({length:p},()=>new Array(p).fill(0));
      const XtYr = new Array(p).fill(0);
      for(let ii=0;ii<nt;ii++)for(let j=0;j<p;j++){XtYr[j]+=Xt[ii][j]*Yt[ii];for(let l=0;l<p;l++)XtX[j][l]+=Xt[ii][j]*Xt[ii][l];}
      for(let j=0;j<p;j++) XtX[j][j] += lam;
      const beta = solveLinear(XtX, XtYr);
      const pred = X[i].reduce((s,x,j)=>s+x*beta[j],0);
      cvSS += (Y[i]-pred)**2;
    }
    results.push({lam, cvRMS: Math.sqrt(cvSS/n)});
  }
  return results;
}

// Standardize orthBundle for Ridge
const Xstd = gals45.map((_,i) => {
  return orthBundle.map(ob => {
    const m = mean(ob), s = sd(ob);
    return s > 1e-10 ? (ob[i]-m)/s : 0;
  });
});

const lambdas = [0.001, 0.01, 0.1, 0.5, 1, 2, 5, 10, 20, 50];
const ridgeResults = ridgeCV(residBdp, Xstd, lambdas);
const bestRidge = ridgeResults.sort((a,b) => a.cvRMS - b.cvRMS)[0];

// Null model: predict 0 for all (mean of residuals ≈ 0)
const nullRMS = sd(residBdp) * Math.sqrt((N-1)/N);
const ridgeGap = 100*(1 - bestRidge.cvRMS**2 / nullRMS**2);

console.log();
console.log('  Ridge regression (standardized ⊥ bundle):');
console.log('    Best λ = ' + bestRidge.lam);
console.log('    LOO CV RMS = ' + bestRidge.cvRMS.toFixed(4));
console.log('    Null RMS = ' + nullRMS.toFixed(4));
console.log('    Ridge gap vs null = ' + ridgeGap.toFixed(1) + '%');
console.log();

// Also LOO for raw augmented model (B″ + bundle)
const X_aug = gals45.map((g,i) => [...X_Bdp[i], ...validBundle.map(ob => ob[i])]);
const looBase = looCV(Y, X_Bdp);
const looAug = looCV(Y, X_aug);
console.log('  Full model LOO:');
console.log('    B″ alone gap% = ' + gapV(looBase).toFixed(1) + '%');
console.log('    B″ + bundle gap% = ' + gapV(looAug).toFixed(1) + '%');
console.log('    Δgap = ' + (gapV(looAug)-gapV(looBase)).toFixed(1) + 'pp');
console.log();

// Permutation test on bundle
const nPerm = 2000;
let permBetter = 0;
for (let p = 0; p < nPerm; p++) {
  // Permute residuals
  const permResid = [...residBdp];
  for (let i = N-1; i > 0; i--) {
    const j = Math.floor(Math.random()*(i+1));
    [permResid[i], permResid[j]] = [permResid[j], permResid[i]];
  }
  const fPerm = ols(permResid, X_orth);
  const r2p = 1 - fPerm.rss/fPerm.tss;
  if (r2p >= r2_bundle) permBetter++;
}
const permP = permBetter / nPerm;
console.log('  Permutation test (bundle R² ≥ observed):');
console.log('    perm p = ' + permP.toFixed(4) + ' (' + nPerm + ' perms)');
console.log();

// ═══════════════════════════════════════════════════════════════════
// 66b — PCA on bundle → latent factors
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  66b: PCA / LATENT FACTORS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Standardize orthBundle columns
const Zstd = orthBundle.map(col => {
  const m = mean(col), s = sd(col);
  return col.map(v => s > 1e-10 ? (v-m)/s : 0);
});

// Compute correlation matrix of standardized ⊥ bundle
const corrMat = Array.from({length:nBundle}, () => new Array(nBundle).fill(0));
for (let i = 0; i < nBundle; i++)
  for (let j = 0; j < nBundle; j++) {
    let s = 0;
    for (let k = 0; k < N; k++) s += Zstd[i][k] * Zstd[j][k];
    corrMat[i][j] = s / (N-1);
  }

// Power iteration for top eigenvalues
function powerIteration(M, nIter=500) {
  const n = M.length;
  let v = Array.from({length:n}, () => Math.random()-0.5);
  let norm = Math.sqrt(v.reduce((s,x)=>s+x*x,0));
  v = v.map(x=>x/norm);
  let eigenvalue = 0;
  for (let iter = 0; iter < nIter; iter++) {
    const Mv = new Array(n).fill(0);
    for(let i=0;i<n;i++)for(let j=0;j<n;j++)Mv[i]+=M[i][j]*v[j];
    eigenvalue = Mv.reduce((s,x,i)=>s+x*v[i],0);
    norm = Math.sqrt(Mv.reduce((s,x)=>s+x*x,0));
    v = Mv.map(x=>x/norm);
  }
  return {eigenvalue, eigenvector: v};
}

function deflate(M, eigenvalue, eigenvector) {
  const n = M.length;
  const Md = M.map(r => [...r]);
  for(let i=0;i<n;i++)for(let j=0;j<n;j++) Md[i][j] -= eigenvalue*eigenvector[i]*eigenvector[j];
  return Md;
}

// Extract top PCs
const pcs = [];
let Mcur = corrMat.map(r=>[...r]);
const totalVar = nBundle; // sum of eigenvalues of corr matrix = p

for (let pc = 0; pc < Math.min(5, nBundle); pc++) {
  const {eigenvalue, eigenvector} = powerIteration(Mcur);
  // Project data onto PC
  const scores = gals45.map((_,i) => {
    return eigenvector.reduce((s, w, j) => s + w * Zstd[j][i], 0);
  });
  const rc = pearsonCorr(residBdp, scores);
  pcs.push({pc: pc+1, eigenvalue, varExplained: eigenvalue/totalVar, eigenvector, scores, r: rc.r, t: rc.t});
  Mcur = deflate(Mcur, eigenvalue, eigenvector);
}

console.log('  ┌─────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ PC  │ Eigenval │ Var%     │ r(resid) │ t        │');
console.log('  ├─────┼──────────┼──────────┼──────────┼──────────┤');
pcs.forEach(p => {
  console.log('  │ PC' + p.pc + ' │ ' + p.eigenvalue.toFixed(3).padStart(6) + '   │ ' + (p.varExplained*100).toFixed(1).padStart(5) + '%   │ ' + p.r.toFixed(3).padStart(7) + '  │ ' + p.t.toFixed(2).padStart(6) + '   │');
});
console.log('  └─────┴──────────┴──────────┴──────────┴──────────┘');
console.log();

// Top loadings for PC1
console.log('  PC1 loadings:');
const pc1 = pcs[0];
const loadPairs = pc1.eigenvector.map((w,i) => ({name: validNames[i], w}));
loadPairs.sort((a,b) => Math.abs(b.w) - Math.abs(a.w));
loadPairs.forEach(lp => console.log('    ' + lp.name.padEnd(22) + ': ' + lp.w.toFixed(3)));
console.log();

// LOO with PC1, PC1+PC2
const X_pc1 = gals45.map((_,i) => [...X_Bdp[i], pcs[0].scores[i]]);
const X_pc12 = gals45.map((_,i) => [...X_Bdp[i], pcs[0].scores[i], pcs[1].scores[i]]);
const looPC1 = looCV(Y, X_pc1);
const looPC12 = looCV(Y, X_pc12);

console.log('  LOO with PCs added to B″:');
console.log('    B″ alone:     gap% = ' + gapV(looBase).toFixed(1) + '%');
console.log('    B″ + PC1:     gap% = ' + gapV(looPC1).toFixed(1) + '% (Δ = ' + (gapV(looPC1)-gapV(looBase)).toFixed(1) + 'pp)');
console.log('    B″ + PC1+PC2: gap% = ' + gapV(looPC12).toFixed(1) + '% (Δ = ' + (gapV(looPC12)-gapV(looBase)).toFixed(1) + 'pp)');
console.log();

// ═══════════════════════════════════════════════════════════════════
// 66c — Stability tests
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  66c: STABILITY — BOOTSTRAP & REPEATED SPLITS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Bootstrap: does bundle R² hold up?
const nBoot = 2000;
const bootR2 = [];
for (let b = 0; b < nBoot; b++) {
  const idx = Array.from({length:N}, () => Math.floor(Math.random()*N));
  try {
    const Yb = idx.map(i => residBdp[i]);
    const Xb = idx.map(i => X_orth[i]);
    const fb = ols(Yb, Xb);
    bootR2.push(1 - fb.rss/fb.tss);
  } catch(e) {}
}
bootR2.sort((a,b) => a-b);
console.log('  Bootstrap R² of bundle on residuals:');
console.log('    Point estimate = ' + r2_bundle.toFixed(4));
console.log('    95% CI = [' + bootR2[Math.floor(0.025*bootR2.length)].toFixed(4) + ', ' + bootR2[Math.floor(0.975*bootR2.length)].toFixed(4) + ']');
console.log('    CI includes 0? ' + (bootR2[Math.floor(0.025*bootR2.length)] <= 0 ? 'YES' : 'NO'));
console.log();

// 50/50 split test
const nSplits = 500;
let splitWins = 0;
const splitDeltaGaps = [];
for (let s = 0; s < nSplits; s++) {
  const perm = gals45.map((_,i) => i);
  for (let i = N-1; i > 0; i--) {
    const j = Math.floor(Math.random()*(i+1));
    [perm[i], perm[j]] = [perm[j], perm[i]];
  }
  const train = perm.slice(0, Math.floor(N/2));
  const test = perm.slice(Math.floor(N/2));
  
  try {
    // Train B″ + bundle
    const Ytr = train.map(i => Y[i]);
    const Xtr_base = train.map(i => X_Bdp[i]);
    const Xtr_aug = train.map(i => X_aug[i]);
    
    const fBase = ols(Ytr, Xtr_base);
    const fAugTr = ols(Ytr, Xtr_aug);
    
    // Test
    let ssBase = 0, ssAug = 0;
    test.forEach(i => {
      const xb = [1,...X_Bdp[i]];
      const xa = [1,...X_aug[i]];
      const predB = xb.reduce((s,x,j)=>s+x*fBase.beta[j],0);
      const predA = xa.reduce((s,x,j)=>s+x*fAugTr.beta[j],0);
      ssBase += (Y[i]-predB)**2;
      ssAug += (Y[i]-predA)**2;
    });
    
    if (ssAug < ssBase) splitWins++;
    const sdTest = sd(test.map(i=>Y[i]));
    const rmsB = Math.sqrt(ssBase/test.length);
    const rmsA = Math.sqrt(ssAug/test.length);
    splitDeltaGaps.push(100*(1-rmsA**2/sdTest**2) - 100*(1-rmsB**2/sdTest**2));
  } catch(e) {}
}

console.log('  50/50 split test (' + nSplits + ' random splits):');
console.log('    B″+bundle beats B″ alone: ' + splitWins + '/' + nSplits + ' (' + (splitWins/nSplits*100).toFixed(1) + '%)');
console.log('    Mean Δgap = ' + mean(splitDeltaGaps).toFixed(1) + 'pp');
console.log('    Median Δgap = ' + splitDeltaGaps.sort((a,b)=>a-b)[Math.floor(splitDeltaGaps.length/2)].toFixed(1) + 'pp');
console.log();

// ═══════════════════════════════════════════════════════════════════
// FINAL VERDICT
// ═══════════════════════════════════════════════════════════════════
console.log('======================================================================');
console.log('  PHASE 66: FINAL VERDICT');
console.log('======================================================================\n');

const crit1_improves = gapV(looAug) > gapV(looBase);
const crit2_cvStable = splitWins/nSplits > 0.55;
const crit3_permSig = permP < 0.05;
const crit4_pcSignal = pcs.some(p => Math.abs(p.t) >= 1.65);

console.log('  Criterion 1 (bundle improves LOO):           ' + (crit1_improves ? 'PASS ✅' : 'FAIL ❌') + ' (Δ = ' + (gapV(looAug)-gapV(looBase)).toFixed(1) + 'pp)');
console.log('  Criterion 2 (CV stable >55% split wins):     ' + (crit2_cvStable ? 'PASS ✅' : 'FAIL ❌') + ' (' + (splitWins/nSplits*100).toFixed(0) + '%)');
console.log('  Criterion 3 (permutation p < 0.05):          ' + (crit3_permSig ? 'PASS ✅' : 'FAIL ❌') + ' (p=' + permP.toFixed(3) + ')');
console.log('  Criterion 4 (latent PC |t| ≥ 1.65):         ' + (crit4_pcSignal ? 'PASS ✅' : 'FAIL ❌'));
console.log();

const nPass = [crit1_improves, crit2_cvStable, crit3_permSig, crit4_pcSignal].filter(Boolean).length;
let verdict;
if (nPass >= 3) verdict = 'CONFIRMED-DISTRIBUTED-RESIDUAL';
else if (nPass >= 2) verdict = 'PARTIAL';
else verdict = 'FAIL';

console.log('  ═══════════════════════════════════════════════════════════');
console.log('  ' + verdict + ' (' + nPass + '/4)');
console.log('  ═══════════════════════════════════════════════════════════');
console.log();

if (verdict === 'FAIL') {
  console.log('  The ~0.10 dex intrinsic scatter after B″ is NOT explained');
  console.log('  by ANY available scalar bundle:');
  console.log('    • Not a single strong axis (Phase 62)');
  console.log('    • Not an environmental processing proxy (Phase 65)');
  console.log('    • Not a distributed weak-signal bundle (Phase 66)');
  console.log();
  console.log('  CONCLUSION:');
  console.log('  B″ is SATURATED on available scalar data.');
  console.log('  The remaining ~0.10 dex is either:');
  console.log('    1. Dynamical state information not captured by 1D scalars');
  console.log('    2. Irreducible intrinsic scatter at current data quality');
  console.log('    3. Both — some combination of unresolved physics + noise floor');
}

const output = {
  phase: '66', title: 'Weak-Bundle / Latent Residual Test',
  residualSource: 'B″ on N=45',
  bundleVars: validNames,
  olsBundleR2: r2_bundle,
  ridgeBest: {lambda: bestRidge.lam, cvRMS: bestRidge.cvRMS, gapVsNull: ridgeGap},
  looBase: gapV(looBase), looAug: gapV(looAug), looDelta: gapV(looAug)-gapV(looBase),
  permP,
  pca: pcs.map(p => ({pc:p.pc, eigenvalue:p.eigenvalue, varPct:p.varExplained*100, r:p.r, t:p.t})),
  looPC1: gapV(looPC1), looPC12: gapV(looPC12),
  splitWinRate: splitWins/nSplits,
  criteria: {crit1_improves, crit2_cvStable, crit3_permSig, crit4_pcSignal, nPass},
  verdict
};
fs.writeFileSync('public/phase66-weak-bundle.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase66-weak-bundle.json');
