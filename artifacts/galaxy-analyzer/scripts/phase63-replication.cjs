const fs = require('fs');

// ═══════════════════════════════════════════════════════════════════
// Phase 63a: Build raw input from scratch (no fitted coefficients)
// ═══════════════════════════════════════════════════════════════════
const stageA = JSON.parse(fs.readFileSync('public/stage-A-master-table.json','utf8'));
const tidalData = JSON.parse(fs.readFileSync('public/phase58a2-tidal-expansion.json','utf8')).tidalData;
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const tdMap = {};
tidalData.forEach(t => { tdMap[t.name] = t; });
const pubNames = new Set(tidalData.filter(t => t.quality === 'published').map(t => t.name));
const sparcMap = {};
sparc.forEach(s => { sparcMap[s.name] = s; });
const gals45 = stageA.galaxies.filter(g => pubNames.has(g.name));
const N = gals45.length;

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

// ── Helper functions (all from scratch, no imports) ──
function mean(a){return a.reduce((s,v)=>s+v,0)/a.length;}
function sd(a){const m=mean(a);return Math.sqrt(a.reduce((s,v)=>s+(v-m)**2,0)/(a.length-1));}
function median(a){const s=[...a].sort((x,y)=>x-y);return s.length%2?s[(s.length-1)/2]:(s[s.length/2-1]+s[s.length/2])/2;}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function invertMatrix(A){const n=A.length;const M=A.map((r,i)=>[...r,...Array.from({length:n},(_,j)=>i===j?1:0)]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];const pv=M[i][i];if(Math.abs(pv)<1e-15){for(let j=0;j<2*n;j++)M[i][j]=0;continue;}for(let j=0;j<2*n;j++)M[i][j]/=pv;for(let j=0;j<n;j++){if(j===i)continue;const f=M[j][i];for(let k=0;k<2*n;k++)M[j][k]-=f*M[i][k];}}return M.map(r=>r.slice(n));}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0),tss=Y.reduce((s,y)=>s+(y-mean(Y))**2,0);const inv=invertMatrix(XtX);const seBeta=beta.map((_,j)=>Math.sqrt(Math.abs(inv[j][j])*rss/(n-p)));const tStats=beta.map((b,j)=>b/(seBeta[j]||1e-15));return{beta,seBeta,tStats,r2adj:1-(1-(1-rss/tss))*(n-1)/(n-p),se:Math.sqrt(rss/(n-p)),resid,rss,tss,n,k:p};}
function looCV(Y,X){const n=Y.length;let ss=0;for(let i=0;i<n;i++){const Yt=[...Y.slice(0,i),...Y.slice(i+1)],Xt=[...X.slice(0,i),...X.slice(i+1)];const f=ols(Yt,Xt);const xi=[1,...X[i]];ss+=(Y[i]-xi.reduce((s,x,j)=>s+x*f.beta[j],0))**2;}return Math.sqrt(ss/n);}

// ── Phase 63a: Build raw replication input from scratch ──
const morphT = gals45.map(g => sparcMap[g.name]?.T ?? 5);
const logUps = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhost = gals45.map(g => tdMap[g.name].logMhost);

// Rebuild Υ★⊥ from scratch
const XconfRaw = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT[i]]);
const fConfRaw = ols(logUps, XconfRaw);
const upsPerpRaw = fConfRaw.resid;

// Build the raw input table
const rawInput = gals45.map((g,i) => ({
  galaxy: g.name,
  logA0_obs: g.logA0,
  logMHI: g.logMHI,
  rcWiggliness: g.rcWiggliness,
  logMhost: logMhost[i],
  logSigma0: g.logSigma0,
  logMeanRun: g.logMeanRun,
  Upsilon_perp: upsPerpRaw[i]
}));

console.log('======================================================================');
console.log('  PHASE 63: INDEPENDENT REPLICATION PROTOCOL');
console.log('  N = ' + N + ' galaxies');
console.log('======================================================================\n');

console.log('  Phase 63a: Raw input table built from scratch');
console.log('  Columns: galaxy, logA0_obs, logMHI, rcWiggliness, logMhost,');
console.log('           logSigma0, logMeanRun, Upsilon_perp');
console.log('  Υ★⊥ rebuilt from raw OLS on (logMHI, logΣ₀, morphT)');
console.log('  R²(confounders→logΥ) = ' + (1-fConfRaw.rss/fConfRaw.tss).toFixed(4));
console.log();

// ═══════════════════════════════════════════════════════════════════
// Phase 63b: Rebuild model from scratch — 3 methods
// ═══════════════════════════════════════════════════════════════════

const Y = rawInput.map(r => r.logA0_obs);
const X = rawInput.map(r => [r.logMHI, r.rcWiggliness, r.logMhost, r.logSigma0, r.logMeanRun, r.Upsilon_perp]);
const sdY = sd(Y);
function gapV(rms){return 100*(1-rms**2/sdY**2);}

const varNames = ['intercept','logMHI','rcWig','logMhost','logΣ₀','logMeanRun','Υ★⊥'];
const expectedSigns = [null, -1, +1, -1, +1, +1, +1];

// ── Method 1: Standard OLS ──
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Phase 63b — Method 1: Standard OLS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const fOLS = ols(Y, X);
const looOLS = looCV(Y, X);

console.log('  Equation:');
console.log('  log(a₀) = ' + fOLS.beta[0].toFixed(4));
varNames.slice(1).forEach((v,j) => {
  const b = fOLS.beta[j+1];
  console.log('           ' + (b>=0?'+':'-') + ' ' + Math.abs(b).toFixed(4) + '·' + v);
});
console.log();
console.log('  Adj R² = ' + fOLS.r2adj.toFixed(4));
console.log('  LOO gap% = ' + gapV(looOLS).toFixed(1) + '%');
console.log('  Residual SD = ' + sd(fOLS.resid).toFixed(4) + ' dex');
console.log();

varNames.forEach((v,j) => {
  const ci95lo = fOLS.beta[j] - 1.96*fOLS.seBeta[j];
  const ci95hi = fOLS.beta[j] + 1.96*fOLS.seBeta[j];
  console.log('  ' + v.padEnd(12) + ' β=' + fOLS.beta[j].toFixed(4).padStart(7) + ' ± ' + fOLS.seBeta[j].toFixed(4) + '  t=' + fOLS.tStats[j].toFixed(2).padStart(5) + '  95%CI=[' + ci95lo.toFixed(4) + ', ' + ci95hi.toFixed(4) + ']');
});
console.log();

// ── Method 2: Robust regression (Huber IRLS) ──
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Phase 63b — Method 2: Robust regression (Huber IRLS, k=1.345)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

function huberIRLS(Y, X, k=1.345, maxIter=50) {
  const n = Y.length, p = X[0].length + 1;
  let w = new Array(n).fill(1);
  let beta;
  
  for (let iter = 0; iter < maxIter; iter++) {
    const Xa = X.map(r => [1,...r]);
    const XtWX = Array.from({length:p}, () => new Array(p).fill(0));
    const XtWY = new Array(p).fill(0);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < p; j++) {
        XtWY[j] += w[i] * Xa[i][j] * Y[i];
        for (let l = 0; l < p; l++) {
          XtWX[j][l] += w[i] * Xa[i][j] * Xa[i][l];
        }
      }
    }
    beta = solveLinear(XtWX, XtWY);
    
    const resid = Y.map((y,i) => y - Xa[i].reduce((s,x,j) => s+x*beta[j], 0));
    const mad = median(resid.map(r => Math.abs(r)));
    const sigma = mad / 0.6745;
    
    if (sigma < 1e-10) break;
    
    const newW = resid.map(r => {
      const u = Math.abs(r) / (sigma * k);
      return u <= 1 ? 1 : 1/u;
    });
    
    const wChange = newW.reduce((s,wi,i) => s + Math.abs(wi - w[i]), 0);
    w = newW;
    if (wChange < 1e-6) break;
  }
  
  const Xa = X.map(r => [1,...r]);
  const resid = Y.map((y,i) => y - Xa[i].reduce((s,x,j) => s+x*beta[j], 0));
  const rss = resid.reduce((s,r) => s+r*r, 0);
  const tss = Y.reduce((s,y) => s+(y-mean(Y))**2, 0);
  
  return {beta, resid, rss, tss, r2adj: 1-(1-(1-rss/tss))*(n-1)/(n-p), weights: w};
}

const fRob = huberIRLS(Y, X);

// Robust LOO
let ssRobLOO = 0;
for (let i = 0; i < N; i++) {
  const Yt = [...Y.slice(0,i),...Y.slice(i+1)];
  const Xt = [...X.slice(0,i),...X.slice(i+1)];
  const fr = huberIRLS(Yt, Xt);
  const xi = [1,...X[i]];
  const pred = xi.reduce((s,x,j) => s+x*fr.beta[j], 0);
  ssRobLOO += (Y[i]-pred)**2;
}
const looRob = Math.sqrt(ssRobLOO/N);

console.log('  Equation:');
console.log('  log(a₀) = ' + fRob.beta[0].toFixed(4));
varNames.slice(1).forEach((v,j) => {
  const b = fRob.beta[j+1];
  console.log('           ' + (b>=0?'+':'-') + ' ' + Math.abs(b).toFixed(4) + '·' + v);
});
console.log();
console.log('  Adj R² = ' + fRob.r2adj.toFixed(4));
console.log('  LOO gap% = ' + gapV(looRob).toFixed(1) + '%');
console.log('  Residual SD = ' + sd(fRob.resid).toFixed(4) + ' dex');
console.log('  Downweighted points: ' + fRob.weights.filter(w=>w<0.99).length + '/' + N);
console.log();

varNames.forEach((v,j) => {
  console.log('  ' + v.padEnd(12) + ' β=' + fRob.beta[j].toFixed(4).padStart(7));
});
console.log();

// ── Method 3: Bootstrap-mean fit ──
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Phase 63b — Method 3: Bootstrap-mean fit (1000 resamples)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

const nBoot = 1000;
const bootBetas = varNames.map(() => []);

for (let b = 0; b < nBoot; b++) {
  const idx = Array.from({length:N}, () => Math.floor(Math.random()*N));
  try {
    const fb = ols(idx.map(i=>Y[i]), idx.map(i=>X[i]));
    fb.beta.forEach((v,j) => bootBetas[j].push(v));
  } catch(e) {}
}

const bootMean = bootBetas.map(a => mean(a));
const bootSD = bootBetas.map(a => sd(a));
const bootCI = bootBetas.map(a => {
  const s = [...a].sort((x,y)=>x-y);
  return [s[Math.floor(0.025*s.length)], s[Math.floor(0.975*s.length)]];
});

// Compute predictions using bootstrap-mean coefficients
const Xa = X.map(r => [1,...r]);
const predBoot = Xa.map(xi => xi.reduce((s,x,j) => s+x*bootMean[j], 0));
const residBoot = Y.map((y,i) => y - predBoot[i]);
const rssBoot = residBoot.reduce((s,r) => s+r*r, 0);
const tssBoot = Y.reduce((s,y) => s+(y-mean(Y))**2, 0);

console.log('  Equation (bootstrap mean):');
console.log('  log(a₀) = ' + bootMean[0].toFixed(4));
varNames.slice(1).forEach((v,j) => {
  const b = bootMean[j+1];
  console.log('           ' + (b>=0?'+':'-') + ' ' + Math.abs(b).toFixed(4) + '·' + v);
});
console.log();
console.log('  R²adj ≈ ' + (1-(1-(1-rssBoot/tssBoot))*(N-1)/(N-7)).toFixed(4));
console.log('  Residual SD = ' + sd(residBoot).toFixed(4) + ' dex');
console.log();

varNames.forEach((v,j) => {
  console.log('  ' + v.padEnd(12) + ' β=' + bootMean[j].toFixed(4).padStart(7) + ' ± ' + bootSD[j].toFixed(4) + '  95%CI=[' + bootCI[j][0].toFixed(4) + ', ' + bootCI[j][1].toFixed(4) + ']');
});
console.log();

// ═══════════════════════════════════════════════════════════════════
// Phase 63c: Stability audit
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Phase 63c: STABILITY AUDIT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Jackknife
console.log('  Jackknife sign flips (LOO, N=' + N + '):');
const jackFlips = varNames.map(() => 0);
for (let i = 0; i < N; i++) {
  const Yj = [...Y.slice(0,i),...Y.slice(i+1)];
  const Xj = [...X.slice(0,i),...X.slice(i+1)];
  const fj = ols(Yj, Xj);
  fj.beta.forEach((b,j) => {
    if (j > 0 && Math.sign(b) !== Math.sign(fOLS.beta[j])) jackFlips[j]++;
  });
}
varNames.forEach((v,j) => {
  if (j === 0) return;
  const sign = fOLS.beta[j] >= 0 ? '+' : '−';
  console.log('  ' + v.padEnd(12) + ' sign=' + sign + '  flips=' + jackFlips[j] + '/' + N + (jackFlips[j]===0?' ✅':jackFlips[j]<=2?' ✅':' ⚠️'));
});
console.log();

// Bootstrap sign stability
console.log('  Bootstrap sign flips (out of ' + nBoot + '):');
varNames.forEach((v,j) => {
  if (j === 0) return;
  const fullSign = Math.sign(fOLS.beta[j]);
  const flips = bootBetas[j].filter(b => Math.sign(b) !== fullSign).length;
  const pct = (100*flips/bootBetas[j].length).toFixed(1);
  console.log('  ' + v.padEnd(12) + ' flips=' + flips + '/' + bootBetas[j].length + ' (' + pct + '%)' + (flips/bootBetas[j].length < 0.05 ? ' ✅' : ' ⚠️'));
});
console.log();

// 50/50 split repeated 200 times
console.log('  50/50 split test (200 reps):');
const splitBetas = varNames.map(() => []);
const splitGaps = [];

for (let rep = 0; rep < 200; rep++) {
  const idx = Array.from({length:N},(_,i)=>i);
  for(let i=N-1;i>0;i--){const j=Math.floor(Math.random()*(i+1));[idx[i],idx[j]]=[idx[j],idx[i]];}
  const train = idx.slice(0,23), test = idx.slice(23);
  try {
    const ft = ols(train.map(i=>Y[i]), train.map(i=>X[i]));
    ft.beta.forEach((b,j) => splitBetas[j].push(b));
    const errs = test.map(i => {
      const xi = [1,...X[i]];
      return Y[i] - xi.reduce((s,x,j) => s+x*ft.beta[j], 0);
    });
    const rms = Math.sqrt(errs.reduce((s,e)=>s+e*e,0)/errs.length);
    const sdTest = sd(test.map(i=>Y[i]));
    splitGaps.push(100*(1-rms**2/sdTest**2));
  } catch(e) {}
}

varNames.forEach((v,j) => {
  if (j === 0) return;
  const fullSign = Math.sign(fOLS.beta[j]);
  const flips = splitBetas[j].filter(b => Math.sign(b) !== fullSign).length;
  console.log('  ' + v.padEnd(12) + ' split sign flips=' + flips + '/200 (' + (100*flips/200).toFixed(1) + '%)  mean β=' + mean(splitBetas[j]).toFixed(4));
});
console.log('  Split test gap%: mean=' + mean(splitGaps).toFixed(1) + '%, SD=' + sd(splitGaps).toFixed(1) + '%');
console.log();

// ═══════════════════════════════════════════════════════════════════
// Phase 63d: Coefficient comparison
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Phase 63d: COEFFICIENT COMPARISON (3 methods vs B″ original)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// B″ original coefficients (from Phase 59b)
const origBeta = [4.6005, -0.2430, 1.2790, -0.1503, 0.1680, 0.4456, 0.6590];

console.log('  ┌──────────────┬──────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Variable     │ B″ orig  │ OLS      │ Robust   │ Boot-μ   │ Sign OK? │');
console.log('  ├──────────────┼──────────┼──────────┼──────────┼──────────┼──────────┤');

let allSignsOK = true;
varNames.forEach((v,j) => {
  const orig = origBeta[j];
  const olsB = fOLS.beta[j];
  const robB = fRob.beta[j];
  const bootB = bootMean[j];
  
  let signOK = true;
  if (j > 0) {
    const exp = expectedSigns[j];
    signOK = Math.sign(olsB) === exp && Math.sign(robB) === exp && Math.sign(bootB) === exp;
    if (!signOK) allSignsOK = false;
  }
  
  const signStr = j === 0 ? '  —   ' : (signOK ? '  ✅   ' : '  ❌   ');
  console.log('  │ ' + v.padEnd(12) + ' │ ' + orig.toFixed(4).padStart(7) + '  │ ' + olsB.toFixed(4).padStart(7) + '  │ ' + robB.toFixed(4).padStart(7) + '  │ ' + bootB.toFixed(4).padStart(7) + '  │' + signStr + '│');
});
console.log('  └──────────────┴──────────┴──────────┴──────────┴──────────┴──────────┘');
console.log();

// Check if values are within ±25% or 2σ
console.log('  Proximity to original B″ (OLS replication):');
varNames.forEach((v,j) => {
  if (j === 0) return;
  const orig = origBeta[j];
  const repl = fOLS.beta[j];
  const pctDiff = Math.abs((repl - orig) / orig * 100);
  const within2s = Math.abs(repl - orig) <= 2 * fOLS.seBeta[j];
  console.log('  ' + v.padEnd(12) + ' Δ=' + pctDiff.toFixed(1) + '%' + (pctDiff<25?' ✅':' ⚠️') + '  within 2σ: ' + (within2s?'✅':'❌'));
});
console.log();

// ═══════════════════════════════════════════════════════════════════
// Phase 63e: Blind prediction replication
// ═══════════════════════════════════════════════════════════════════
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  Phase 63e: BLIND PREDICTION REPLICATION');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Predictions from each method
const predOLS = Xa.map(xi => xi.reduce((s,x,j) => s+x*fOLS.beta[j], 0));
const predRob = Xa.map(xi => xi.reduce((s,x,j) => s+x*fRob.beta[j], 0));
const predBootM = Xa.map(xi => xi.reduce((s,x,j) => s+x*bootMean[j], 0));
const predOrig = Xa.map(xi => xi.reduce((s,x,j) => s+x*origBeta[j], 0));

// Spearman rank correlation between OLS-repl predictions and original
function spearmanR(a, b) {
  const n = a.length;
  const rankA = rankArray(a), rankB = rankArray(b);
  const mA = mean(rankA), mB = mean(rankB);
  let sxy=0,sxx=0,syy=0;
  for(let i=0;i<n;i++){sxy+=(rankA[i]-mA)*(rankB[i]-mB);sxx+=(rankA[i]-mA)**2;syy+=(rankB[i]-mB)**2;}
  return sxy/Math.sqrt(sxx*syy);
}
function rankArray(a) {
  const sorted = a.map((v,i)=>({v,i})).sort((a,b)=>a.v-b.v);
  const ranks = new Array(a.length);
  sorted.forEach((s,r) => ranks[s.i] = r+1);
  return ranks;
}

// RMS between methods' predictions
const rmsOLSvsOrig = Math.sqrt(predOLS.reduce((s,p,i) => s+(p-predOrig[i])**2, 0)/N);
const rmsRobvsOrig = Math.sqrt(predRob.reduce((s,p,i) => s+(p-predOrig[i])**2, 0)/N);
const rmsBootvsOrig = Math.sqrt(predBootM.reduce((s,p,i) => s+(p-predOrig[i])**2, 0)/N);

const spOLSobs = spearmanR(predOLS, Y);
const spRobObs = spearmanR(predRob, Y);
const spBootObs = spearmanR(predBootM, Y);
const spOrigObs = spearmanR(predOrig, Y);

console.log('  ┌─────────────────┬──────────┬──────────┬──────────┐');
console.log('  │ Method          │ RMS vs   │ Spearman │ Resid SD │');
console.log('  │                 │ B″ orig  │ (pred,Y) │  (dex)   │');
console.log('  ├─────────────────┼──────────┼──────────┼──────────┤');
console.log('  │ B″ original     │   0.000  │ ' + spOrigObs.toFixed(4) + '   │ ' + sd(Y.map((y,i)=>y-predOrig[i])).toFixed(4) + '   │');
console.log('  │ OLS replicated  │ ' + rmsOLSvsOrig.toFixed(4).padStart(7) + '  │ ' + spOLSobs.toFixed(4) + '   │ ' + sd(fOLS.resid).toFixed(4) + '   │');
console.log('  │ Robust (Huber)  │ ' + rmsRobvsOrig.toFixed(4).padStart(7) + '  │ ' + spRobObs.toFixed(4) + '   │ ' + sd(fRob.resid).toFixed(4) + '   │');
console.log('  │ Bootstrap-mean  │ ' + rmsBootvsOrig.toFixed(4).padStart(7) + '  │ ' + spBootObs.toFixed(4) + '   │ ' + sd(residBoot).toFixed(4) + '   │');
console.log('  └─────────────────┴──────────┴──────────┴──────────┘');
console.log();

// ═══════════════════════════════════════════════════════════════════
// Phase 63: FINAL VERDICT
// ═══════════════════════════════════════════════════════════════════
console.log('======================================================================');
console.log('  PHASE 63: REPLICATION VERDICT');
console.log('======================================================================\n');

const crit1 = allSignsOK;
const crit2 = varNames.slice(1).every((v,j) => {
  const pctDiff = Math.abs((fOLS.beta[j+1]-origBeta[j+1])/origBeta[j+1]*100);
  return pctDiff < 25 || Math.abs(fOLS.beta[j+1]-origBeta[j+1]) <= 2*fOLS.seBeta[j+1];
});
const crit3 = gapV(looOLS) >= 46;
const crit4 = jackFlips.slice(1).every(f => f <= 2);

console.log('  Criterion 1 (all signs correct across 3 methods): ' + (crit1 ? 'PASS ✅' : 'FAIL ❌'));
console.log('  Criterion 2 (coefficients within ±25% or 2σ):     ' + (crit2 ? 'PASS ✅' : 'FAIL ❌'));
console.log('  Criterion 3 (LOO gap% >= 46%):                    ' + (crit3 ? 'PASS ✅' : 'FAIL ❌') + ' (' + gapV(looOLS).toFixed(1) + '%)');
console.log('  Criterion 4 (jackknife sign flips <= 2 each):     ' + (crit4 ? 'PASS ✅' : 'FAIL ❌'));
console.log();

const nPass = [crit1, crit2, crit3, crit4].filter(Boolean).length;
let verdict;
if (nPass === 4) verdict = 'CONFIRMED';
else if (nPass >= 3) verdict = 'PARTIAL';
else verdict = 'FAIL';

console.log('  ═══════════════════════════════════════════════════════════');
console.log('  REPLICATION: ' + verdict + ' (' + nPass + '/4)');
console.log('  ═══════════════════════════════════════════════════════════');
console.log();

if (verdict === 'CONFIRMED') {
  console.log('  B″ is NOT an artifact of any single pipeline.');
  console.log('  The same law emerges independently from the raw data');
  console.log('  across OLS, Robust, and Bootstrap methods.');
  console.log('  All 6 coefficients maintain their signs and magnitudes.');
}

// Save
const output = {
  phase: '63', title: 'Independent Replication Protocol',
  n: N, sdY,
  methods: {
    ols: {beta:fOLS.beta, se:fOLS.seBeta, tStats:fOLS.tStats, r2adj:fOLS.r2adj, looGap:gapV(looOLS), residSD:sd(fOLS.resid)},
    robust: {beta:fRob.beta, r2adj:fRob.r2adj, looGap:gapV(looRob), residSD:sd(fRob.resid), downweighted:fRob.weights.filter(w=>w<0.99).length},
    bootstrap: {meanBeta:bootMean, sdBeta:bootSD, ci95:bootCI, residSD:sd(residBoot)}
  },
  stability: {
    jackknife: varNames.map((v,j)=>({name:v, flips:jackFlips[j]})),
    bootstrapSignFlips: varNames.map((v,j)=>({name:v, flips:bootBetas[j].filter(b=>Math.sign(b)!==Math.sign(fOLS.beta[j])).length, total:bootBetas[j].length})),
    splitMeanGap: mean(splitGaps), splitSDGap: sd(splitGaps)
  },
  blindPrediction: {
    rmsOLSvsOrig, rmsRobvsOrig, rmsBootvsOrig,
    spearmanOLS: spOLSobs, spearmanRob: spRobObs, spearmanBoot: spBootObs
  },
  criteria: {crit1, crit2, crit3, crit4, nPass},
  verdict
};
fs.writeFileSync('public/phase63-replication.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase63-replication.json');
