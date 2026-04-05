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
function median(a){const s=[...a].sort((x,y)=>x-y);return s.length%2?s[(s.length-1)/2]:(s[s.length/2-1]+s[s.length/2])/2;}
function solveLinear(A,b){const n=A.length;const M=A.map((r,i)=>[...r,b[i]]);for(let i=0;i<n;i++){let mx=i;for(let j=i+1;j<n;j++)if(Math.abs(M[j][i])>Math.abs(M[mx][i]))mx=j;[M[i],M[mx]]=[M[mx],M[i]];if(Math.abs(M[i][i])<1e-15)continue;for(let j=i+1;j<n;j++){const f=M[j][i]/M[i][i];for(let k=i;k<=n;k++)M[j][k]-=f*M[i][k];}}const x=new Array(n);for(let i=n-1;i>=0;i--){x[i]=M[i][n];for(let j=i+1;j<n;j++)x[i]-=M[i][j]*x[j];x[i]/=M[i][i];}return x;}
function ols(Y,X){const n=Y.length,p=X[0].length+1;const Xa=X.map(r=>[1,...r]);const XtX=Array.from({length:p},()=>new Array(p).fill(0)),XtY=new Array(p).fill(0);for(let i=0;i<n;i++)for(let j=0;j<p;j++){XtY[j]+=Xa[i][j]*Y[i];for(let l=0;l<p;l++)XtX[j][l]+=Xa[i][j]*Xa[i][l];}const beta=solveLinear(XtX,XtY);const resid=Y.map((y,i)=>y-Xa[i].reduce((s,x,j)=>s+x*beta[j],0));const rss=resid.reduce((s,r)=>s+r*r,0);return{beta,resid,rss};}

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
const morphT = gals45.map(g => sparcMap[g.name]?.T ?? 5);
const logUps = gals45.map(g => Math.log10(upsilonMap[g.name] || 0.50));
const logMhost = gals45.map(g => tdMap[g.name].logMhost);
const Xconf = gals45.map((g,i) => [g.logMHI, g.logSigma0, morphT[i]]);
const fConf = ols(logUps, Xconf);
const upsPerp = fConf.resid;

const logMHI = gals45.map(g => g.logMHI);
const logSig0 = gals45.map(g => g.logSigma0);
const logMR = gals45.map(g => g.logMeanRun);
const names = gals45.map(g => g.name);

const varArrays = [logMHI, logMhost, logSig0, logMR, upsPerp];
const varNames = ['logMHI','logMhost','logSig0','logMeanRun','Ups_perp'];
const axisLabels = ['Gas mass','Host mass','Baryon density','Kinematic coherence','Stellar structure'];
const sdVars = varArrays.map(v => sd(v));

// M5 fit
const X_M5 = gals45.map((_,i) => [logMHI[i], logMhost[i], logSig0[i], logMR[i], upsPerp[i]]);
const fM5 = ols(Y, X_M5);
const coeffs = fM5.beta.slice(1); // 5 coefficients

console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 72: MATCHED-PAIR TEST');
console.log('  Controlled pairwise predictions from M5');
console.log('══════════════════════════════════════════════════════════════════════\n');

// For each target axis, find top 5 pairs ranked by quality score
// Quality = targetDiff_sigma / confoundDist_sigma
// where confoundDist = sqrt(sum of squared standardized diffs on non-target axes)

function findTopPairs(targetIdx, nPairs) {
  const pairs = [];
  for (let i = 0; i < N; i++) {
    for (let j = i+1; j < N; j++) {
      const targetDiff = Math.abs(varArrays[targetIdx][i] - varArrays[targetIdx][j]) / sdVars[targetIdx];
      let confound2 = 0;
      for (let k = 0; k < 5; k++) {
        if (k === targetIdx) continue;
        confound2 += ((varArrays[k][i] - varArrays[k][j]) / sdVars[k])**2;
      }
      const confoundDist = Math.sqrt(confound2);
      const quality = targetDiff / (confoundDist + 0.1); // +0.1 to avoid div by 0
      pairs.push({i, j, targetDiff, confoundDist, quality});
    }
  }
  pairs.sort((a,b) => b.quality - a.quality);
  return pairs.slice(0, nPairs);
}

const allFamilyResults = [];
const allPairDetails = [];

for (let t = 0; t < 5; t++) {
  console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
  console.log('  72' + 'abcde'[t] + ': ' + axisLabels[t].toUpperCase() + ' PAIRS (vary ' + varNames[t] + ')');
  console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
  
  const topPairs = findTopPairs(t, 5);
  let signPasses = 0;
  const absErrors = [];
  const familyPairs = [];
  
  for (let p = 0; p < topPairs.length; p++) {
    const {i, j, targetDiff, confoundDist, quality} = topPairs[p];
    
    // Orient: hiIdx has higher target value
    const hiIdx = varArrays[t][i] > varArrays[t][j] ? i : j;
    const loIdx = hiIdx === i ? j : i;
    
    const rawDelta = varArrays[t][hiIdx] - varArrays[t][loIdx];
    
    // Partial prediction (from target axis alone)
    const partialPred = coeffs[t] * rawDelta;
    
    // Full model prediction
    const predHi = [1, ...X_M5[hiIdx]].reduce((s,x,k) => s+x*fM5.beta[k], 0);
    const predLo = [1, ...X_M5[loIdx]].reduce((s,x,k) => s+x*fM5.beta[k], 0);
    const fullPred = predHi - predLo;
    
    // Observed
    const observed = Y[hiIdx] - Y[loIdx];
    
    const predictedSign = coeffs[t] > 0 ? '+' : '-';
    const observedSign = observed > 0 ? '+' : '-';
    const signMatch = predictedSign === observedSign;
    if (signMatch) signPasses++;
    
    const absErr = Math.abs(fullPred - observed);
    absErrors.push(absErr);
    
    // Confound details
    const confDetails = [];
    for (let k = 0; k < 5; k++) {
      if (k === t) continue;
      const d = Math.abs(varArrays[k][hiIdx] - varArrays[k][loIdx]) / sdVars[k];
      confDetails.push(varNames[k].replace('log','').substring(0,6) + '=' + d.toFixed(2) + 'σ');
    }
    
    const pairResult = {
      rank: p+1, axis: varNames[t], axisLabel: axisLabels[t],
      galA: names[hiIdx], galB: names[loIdx],
      targetDelta: rawDelta, targetSigma: targetDiff,
      confoundDist, quality,
      partialPred, fullPred, observed,
      signMatch, absErr,
      confDetails
    };
    familyPairs.push(pairResult);
    allPairDetails.push(pairResult);
    
    console.log('  Pair ' + (p+1) + ': ' + names[hiIdx] + ' vs ' + names[loIdx]);
    console.log('    Quality = ' + quality.toFixed(2) + '  |  Δ(' + varNames[t] + ') = ' + rawDelta.toFixed(3) + ' (' + targetDiff.toFixed(1) + 'σ)  |  Confound = ' + confoundDist.toFixed(2) + 'σ');
    console.log('    Confounds: ' + confDetails.join(', '));
    console.log('    Predicted Δlog(a₀): partial=' + (partialPred>0?'+':'') + partialPred.toFixed(3) + ', full=' + (fullPred>0?'+':'') + fullPred.toFixed(3));
    console.log('    Observed  Δlog(a₀): ' + (observed>0?'+':'') + observed.toFixed(3));
    console.log('    Sign: ' + (signMatch ? 'PASS ✅' : 'FAIL ❌') + '  |  |pred−obs| = ' + absErr.toFixed(3) + ' dex');
    console.log();
  }
  
  const medAbsErr = median(absErrors);
  allFamilyResults.push({
    axis: varNames[t], label: axisLabels[t],
    nPairs: topPairs.length, signPasses,
    medianAbsErr: medAbsErr,
    bestPair: familyPairs[0],
    pairs: familyPairs
  });
  
  console.log('  FAMILY SUMMARY: ' + signPasses + '/' + topPairs.length + ' sign passes, median |pred−obs| = ' + medAbsErr.toFixed(3) + ' dex');
  console.log();
}

// ═══════════════════════════════════════════════════════════════════
// OVERALL SUMMARY
// ═══════════════════════════════════════════════════════════════════
console.log('══════════════════════════════════════════════════════════════════════');
console.log('  PHASE 72: OVERALL SUMMARY');
console.log('══════════════════════════════════════════════════════════════════════\n');

const totalSignPasses = allFamilyResults.reduce((s,f) => s+f.signPasses, 0);
const totalPairs = allFamilyResults.reduce((s,f) => s+f.nPairs, 0);

console.log('  ┌────────────────────────┬─────────┬────────────┬───────────────┐');
console.log('  │ Family                 │  Sign   │ Median err │ Best quality  │');
console.log('  ├────────────────────────┼─────────┼────────────┼───────────────┤');
allFamilyResults.forEach(f => {
  console.log('  │ ' + f.label.padEnd(22) + ' │  ' + f.signPasses + '/' + f.nPairs + '  │   ' + f.medianAbsErr.toFixed(3) + '    │    ' + f.bestPair.quality.toFixed(2) + '       │');
});
console.log('  ├────────────────────────┼─────────┼────────────┼───────────────┤');
const overallMedianErr = median(allPairDetails.map(p => p.absErr));
console.log('  │ TOTAL                  │ ' + totalSignPasses + '/' + totalPairs + '  │   ' + overallMedianErr.toFixed(3) + '    │               │');
console.log('  └────────────────────────┴─────────┴────────────┴───────────────┘');
console.log();

// Best-controlled family
const bestFamily = allFamilyResults.reduce((best,f) => f.bestPair.quality > best.bestPair.quality ? f : best);
const worstFamily = allFamilyResults.reduce((worst,f) => f.signPasses < worst.signPasses || (f.signPasses === worst.signPasses && f.medianAbsErr > worst.medianAbsErr) ? f : worst);

console.log('  Best-controlled family: ' + bestFamily.label + ' (quality=' + bestFamily.bestPair.quality.toFixed(2) + ')');
console.log('  Weakest family: ' + worstFamily.label + ' (' + worstFamily.signPasses + '/' + worstFamily.nPairs + ' passes)');
console.log();

// Top 3 pairs overall (by quality, among those that pass)
const passingPairs = allPairDetails.filter(p => p.signMatch).sort((a,b) => b.quality - a.quality);
console.log('  ═══════════════════════════════════════════════════════════════');
console.log('  TOP 3 SHOWCASE PAIRS (highest quality, sign-correct):');
console.log('  ═══════════════════════════════════════════════════════════════');
console.log();

for (let i = 0; i < Math.min(3, passingPairs.length); i++) {
  const p = passingPairs[i];
  console.log('  #' + (i+1) + ': ' + p.galA + ' vs ' + p.galB);
  console.log('      Axis: ' + p.axisLabel + ' (' + p.axis + ')');
  console.log('      Δtarget = ' + p.targetDelta.toFixed(3) + ' (' + p.targetSigma.toFixed(1) + 'σ), confound = ' + p.confoundDist.toFixed(2) + 'σ, quality = ' + p.quality.toFixed(2));
  console.log('      Full-model predicted Δlog(a₀) = ' + (p.fullPred>0?'+':'') + p.fullPred.toFixed(3));
  console.log('      Observed Δlog(a₀) = ' + (p.observed>0?'+':'') + p.observed.toFixed(3));
  console.log('      |error| = ' + p.absErr.toFixed(3) + ' dex');
  console.log();
}

// Quality-weighted analysis: do higher-quality pairs predict better?
const sortedByQ = [...allPairDetails].sort((a,b) => b.quality - a.quality);
const topHalf = sortedByQ.slice(0, Math.floor(sortedByQ.length/2));
const botHalf = sortedByQ.slice(Math.floor(sortedByQ.length/2));
const topPassRate = topHalf.filter(p=>p.signMatch).length / topHalf.length * 100;
const botPassRate = botHalf.filter(p=>p.signMatch).length / botHalf.length * 100;
const topMedianErr = median(topHalf.map(p=>p.absErr));
const botMedianErr = median(botHalf.map(p=>p.absErr));

console.log('  Quality-stratified analysis:');
console.log('    Top-half quality pairs: sign pass = ' + topPassRate.toFixed(0) + '%, median |err| = ' + topMedianErr.toFixed(3));
console.log('    Bottom-half quality:    sign pass = ' + botPassRate.toFixed(0) + '%, median |err| = ' + botMedianErr.toFixed(3));
console.log('    Better-controlled pairs predict better? ' + (topPassRate >= botPassRate && topMedianErr <= botMedianErr ? 'YES ✅' : (topPassRate >= botPassRate ? 'PARTIALLY (sign yes, magnitude no)' : 'NO ❌')));
console.log();

// Final verdict
const passRate = totalSignPasses / totalPairs * 100;
let verdict;
if (passRate >= 80 && overallMedianErr < 0.20) verdict = 'CONFIRMED-PAIRWISE';
else if (passRate >= 60) verdict = 'PARTIAL';
else verdict = 'FAIL';

console.log('  ═══════════════════════════════════════════════════════════════');
console.log('  FINAL VERDICT: ' + verdict);
console.log('  ═══════════════════════════════════════════════════════════════');
console.log();
console.log('  Total sign passes: ' + totalSignPasses + '/' + totalPairs + ' (' + passRate.toFixed(0) + '%)');
console.log('  Median |pred−obs|: ' + overallMedianErr.toFixed(3) + ' dex');
console.log();

if (verdict === 'CONFIRMED-PAIRWISE') {
  console.log('  M5 succeeds not only as a global regression but also in');
  console.log('  controlled galaxy-to-galaxy pairwise comparisons.');
  console.log('  The law predicts the DIRECTION and approximate MAGNITUDE');
  console.log('  of a₀ differences between matched galaxies.');
}

const output = {
  phase:'72', title:'Matched-Pair Test',
  families: allFamilyResults.map(f => ({
    axis:f.axis, label:f.label, nPairs:f.nPairs, signPasses:f.signPasses,
    medianAbsErr:f.medianAbsErr,
    pairs:f.pairs.map(p => ({
      galA:p.galA, galB:p.galB, targetDelta:p.targetDelta, targetSigma:p.targetSigma,
      confoundDist:p.confoundDist, quality:p.quality,
      partialPred:p.partialPred, fullPred:p.fullPred, observed:p.observed,
      signMatch:p.signMatch, absErr:p.absErr
    }))
  })),
  topShowcasePairs: passingPairs.slice(0,3).map(p => ({galA:p.galA, galB:p.galB, axis:p.axisLabel, quality:p.quality, fullPred:p.fullPred, observed:p.observed, absErr:p.absErr})),
  qualityStratified: {topPassRate, botPassRate, topMedianErr, botMedianErr},
  summary: {totalSignPasses, totalPairs, passRate, overallMedianErr, verdict}
};
fs.writeFileSync('public/phase72-matched-pairs.json', JSON.stringify(output, null, 2));
console.log('\n  Saved to public/phase72-matched-pairs.json');
