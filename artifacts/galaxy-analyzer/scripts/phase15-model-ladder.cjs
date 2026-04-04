#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "15.0.0";

function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(80)); }

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function randn() {
  let u=0,v=0; while(!u) u=Math.random(); while(!v) v=Math.random();
  return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);
}

function fitA0(pts) {
  let lo=2.0,hi=5.0;
  for(let s=0;s<150;s++){
    const m1=lo+(hi-lo)*0.382,m2=lo+(hi-lo)*0.618;
    let c1=0,c2=0;
    for(const p of pts){
      const gb=Math.pow(10,p.log_g_bar);
      c1+=(p.log_g_obs-Math.log10(mcgaughRAR(gb,Math.pow(10,m1))))**2;
      c2+=(p.log_g_obs-Math.log10(mcgaughRAR(gb,Math.pow(10,m2))))**2;
    }
    if(c1<c2) hi=m2; else lo=m1;
  }
  const logA0=(lo+hi)/2, a0=Math.pow(10,logA0);
  let ss=0;
  for(const p of pts){const gb=Math.pow(10,p.log_g_bar); ss+=(p.log_g_obs-Math.log10(mcgaughRAR(gb,a0)))**2;}
  return {a0,logA0,rms:Math.sqrt(ss/pts.length),ss,n:pts.length};
}

function predictSS(a0,pts){
  let ss=0;
  for(const p of pts){const gb=Math.pow(10,p.log_g_bar);const pred=mcgaughRAR(gb,a0);ss+=(p.log_g_obs-Math.log10(pred>0?pred:1e-10))**2;}
  return ss;
}

function dlMeta(logA0s,ses){
  const n=logA0s.length,w=ses.map(s=>1/(s*s)),wS=w.reduce((a,b)=>a+b,0);
  const muFE=logA0s.reduce((s,v,i)=>s+w[i]*v,0)/wS;
  let Q=0; for(let i=0;i<n;i++) Q+=w[i]*(logA0s[i]-muFE)**2;
  const S1=wS,S2=w.reduce((s,v)=>s+v*v,0);
  const tau2=Math.max(0,(Q-(n-1))/(S1-S2/S1)),tau=Math.sqrt(tau2);
  const W=ses.map(s=>1/(s*s+tau2)),WS=W.reduce((a,b)=>a+b,0);
  const mu=logA0s.reduce((s,v,i)=>s+W[i]*v,0)/WS;
  return {mu,tau,tau2,a0:Math.pow(10,mu)};
}

const p11=JSON.parse(fs.readFileSync(path.join(__dirname,'../public/phase11-sensitivity-lab.json'),'utf8'));
const sparc=JSON.parse(fs.readFileSync(path.join(__dirname,'../public/sparc-table.json'),'utf8'));

const galaxyData=[];
for(const gal of p11.galaxies){
  const pts=gal.localProfile.filter(p=>isFinite(p.log_g_bar)&&isFinite(p.log_g_obs)&&p.log_g_obs>p.log_g_bar*0.5);
  if(pts.length<5) continue;
  const sp=sparc.find(s=>(s.name||s.Galaxy)===gal.name);
  let etaRot=0;
  if(pts.length>=4){
    const oh=pts.slice(Math.floor(pts.length/2)),ih=pts.slice(0,Math.floor(pts.length/2));
    if(oh.length>=2&&ih.length>=2){
      etaRot=(oh[oh.length-1].log_g_obs-oh[0].log_g_obs)/(oh[oh.length-1].log_g_bar-oh[0].log_g_bar+1e-10)-(ih[ih.length-1].log_g_obs-ih[0].log_g_obs)/(ih[ih.length-1].log_g_bar-ih[0].log_g_bar+1e-10);
    }
  }
  galaxyData.push({
    name:gal.name, pts, inc:gal.inc, D:gal.D, Vmax:gal.Vmax,
    T:gal.T||(sp?sp.T:5), etaRot, n:pts.length,
    logMHI:sp?Math.log10(sp.MHI||1e9):9,
    Rdisk:sp?(sp.Rdisk||3):3,
    SBdisk:sp?(sp.SBdisk||100):100,
    logL36:sp?Math.log10(sp.L||1e9):9,
  });
}

const N=galaxyData.length;
const perGalFits=galaxyData.map(g=>fitA0(g.pts));
const allLogA0=perGalFits.map(f=>f.logA0);
const allSE=perGalFits.map(f=>f.rms/Math.sqrt(f.n));
const baseline=dlMeta(allLogA0,allSE);
const globalA0=baseline.a0;
const globalLogA0=baseline.mu;

log("");
log("=".repeat(80));
log("  PHASE 15: MODEL LADDER — FROM M0 TO M3");
log("  a0,i = a0_global + beta1*logMHI + beta2*X2 + beta3*X3");
log("  Version "+VERSION+", "+N+" galaxies");
log("=".repeat(80));
log("");
log("  Baseline: a0 = "+Math.round(globalA0)+" (log="+globalLogA0.toFixed(4)+")");
log("");

const featureNames=['logMHI','Rdisk','logD','inc','T','etaRot','logVmax','SBdisk','logL36','n'];
function getFeature(g,fname){
  if(fname==='logMHI') return g.logMHI;
  if(fname==='Rdisk') return g.Rdisk;
  if(fname==='logD') return Math.log10(g.D);
  if(fname==='inc') return g.inc;
  if(fname==='T') return g.T;
  if(fname==='etaRot') return g.etaRot;
  if(fname==='logVmax') return Math.log10(g.Vmax);
  if(fname==='SBdisk') return g.SBdisk;
  if(fname==='logL36') return g.logL36;
  if(fname==='n') return g.n;
  return 0;
}

function linReg(X,y){
  const n=y.length,p=X[0]?X[0].length:0;
  if(p===0){const m=y.reduce((s,v)=>s+v,0)/n; return {coefs:[],intercept:m};}
  const means=[];
  for(let j=0;j<p;j++){let s=0;for(let i=0;i<n;i++)s+=X[i][j];means.push(s/n);}
  const my=y.reduce((s,v)=>s+v,0)/n;
  const Xc=X.map(row=>row.map((v,j)=>v-means[j]));
  const yc=y.map(v=>v-my);
  const XtX=Array.from({length:p},()=>Array(p).fill(0));
  const Xty=Array(p).fill(0);
  for(let i=0;i<n;i++){for(let j=0;j<p;j++){Xty[j]+=Xc[i][j]*yc[i];for(let k=0;k<p;k++)XtX[j][k]+=Xc[i][j]*Xc[i][k];}}
  for(let j=0;j<p;j++) XtX[j][j]+=1e-10;
  const A=XtX.map((r,i)=>[...r,Xty[i]]);
  for(let j=0;j<p;j++){let mx=j;for(let i=j+1;i<p;i++)if(Math.abs(A[i][j])>Math.abs(A[mx][j]))mx=i;[A[j],A[mx]]=[A[mx],A[j]];for(let i=j+1;i<p;i++){const f=A[i][j]/A[j][j];for(let k=j;k<=p;k++)A[i][k]-=f*A[j][k];}}
  const b=Array(p).fill(0);
  for(let j=p-1;j>=0;j--){b[j]=A[j][p];for(let k=j+1;k<p;k++)b[j]-=A[j][k]*b[k];b[j]/=A[j][j];}
  const intercept=my-b.reduce((s,v,j)=>s+v*means[j],0);
  return {coefs:b,intercept};
}

function predictLogA0(model,g){
  const features=model.features;
  const X=features.map(f=>getFeature(g,f));
  return model.intercept+model.coefs.reduce((s,c,j)=>s+c*X[j],0);
}

const models=[
  {name:'M0: Universal a0',features:[],desc:'One a0 for all'},
  {name:'M1: MHI only',features:['logMHI'],desc:'a0 + beta*logMHI'},
  {name:'M2: Systematics only',features:['inc','logD'],desc:'a0 + inc + logD corrections'},
  {name:'M3: MHI+sys',features:['logMHI','inc','logD'],desc:'a0 + MHI + inc + D'},
  {name:'M4: MHI+sys+morph',features:['logMHI','inc','logD','T','Rdisk'],desc:'a0 + MHI+inc+D+T+Rdisk'},
  {name:'M5: Kitchen sink',features:['logMHI','inc','logD','T','Rdisk','etaRot','logVmax','SBdisk'],desc:'Everything'},
  {name:'M6: Per-galaxy (M3-free)',features:['__pergalaxy__'],desc:'Free a0 per galaxy'},
];

log("  MODEL DEFINITIONS:");
sep();
for(const m of models) log("  "+m.name.padEnd(25)+" k="+((m.features[0]==='__pergalaxy__'?N:m.features.length+1).toString()).padStart(3)+"  "+m.desc);
log("");

log("  LOO CROSS-VALIDATION (galaxy-level)");
sep();
log("");

const cvResults=[];

for(const model of models){
  let totalSS=0,totalN=0;

  if(model.features[0]==='__pergalaxy__'){
    for(let i=0;i<N;i++){
      const fit=fitA0(galaxyData[i].pts);
      totalSS+=predictSS(fit.a0,galaxyData[i].pts);
      totalN+=galaxyData[i].pts.length;
    }
    cvResults.push({name:model.name,cvRMS:Math.sqrt(totalSS/totalN),totalSS,totalN,k:N,features:model.features});
    continue;
  }

  for(let i=0;i<N;i++){
    const trainIdx=[...Array(N).keys()].filter(j=>j!==i);
    const trainY=trainIdx.map(j=>allLogA0[j]);

    if(model.features.length===0){
      const mu=trainY.reduce((s,v)=>s+v,0)/trainY.length;
      totalSS+=predictSS(Math.pow(10,mu),galaxyData[i].pts);
    } else {
      const trainX=trainIdx.map(j=>model.features.map(f=>getFeature(galaxyData[j],f)));
      const reg=linReg(trainX,trainY);
      const testX=model.features.map(f=>getFeature(galaxyData[i],f));
      const predLogA0=reg.intercept+reg.coefs.reduce((s,c,j)=>s+c*testX[j],0);
      totalSS+=predictSS(Math.pow(10,predLogA0),galaxyData[i].pts);
    }
    totalN+=galaxyData[i].pts.length;
  }
  cvResults.push({name:model.name,cvRMS:Math.sqrt(totalSS/totalN),totalSS,totalN,k:model.features.length+1,features:model.features});
}

const m0rms=cvResults[0].cvRMS;
const m6rms=cvResults[cvResults.length-1].cvRMS;
const gap=m0rms-m6rms;

log("  ┌──────────────────────────────────────────────────────────────────────────┐");
log("  │  Model                    k    CV-RMS    vs M0     vs M6    gap-closed  │");
log("  ├──────────────────────────────────────────────────────────────────────────┤");
for(const r of cvResults){
  const vsM0=((1-r.cvRMS/m0rms)*100).toFixed(1);
  const vsM6=((r.cvRMS/m6rms-1)*100).toFixed(1);
  const gapClosed=gap>0?((m0rms-r.cvRMS)/gap*100).toFixed(1):'0.0';
  log("  │  "+r.name.padEnd(25)+r.k.toString().padStart(3)+r.cvRMS.toFixed(5).padStart(10)+
    (vsM0+"%").padStart(9)+("+"+vsM6+"%").padStart(9)+
    (gapClosed+"%").padStart(12)+"  │");
}
log("  └──────────────────────────────────────────────────────────────────────────┘");
log("");

log("  INFORMATION CRITERIA (full-sample fit)");
sep();
log("");

const icResults=[];
const totalPts=galaxyData.reduce((s,g)=>s+g.pts.length,0);

for(const model of models){
  let totalSS=0,k;
  if(model.features[0]==='__pergalaxy__'){
    for(let i=0;i<N;i++) totalSS+=predictSS(perGalFits[i].a0,galaxyData[i].pts);
    k=N;
  } else {
    k=model.features.length+1;
    if(model.features.length===0){
      for(const g of galaxyData) totalSS+=predictSS(globalA0,g.pts);
    } else {
      const X=galaxyData.map(g=>model.features.map(f=>getFeature(g,f)));
      const reg=linReg(X,allLogA0);
      for(let i=0;i<N;i++){
        const predLogA0=reg.intercept+reg.coefs.reduce((s,c,j)=>s+c*X[i][j],0);
        totalSS+=predictSS(Math.pow(10,predLogA0),galaxyData[i].pts);
      }
    }
  }
  const sigma2=totalSS/totalPts;
  const logL=-totalPts/2*Math.log(2*Math.PI*sigma2)-totalSS/(2*sigma2);
  const AIC=-2*logL+2*k;
  const BIC=-2*logL+k*Math.log(totalPts);
  icResults.push({name:model.name,k,totalSS,AIC,BIC,logL});
}

const minAIC=Math.min(...icResults.map(r=>r.AIC));
const minBIC=Math.min(...icResults.map(r=>r.BIC));

log("  ┌────────────────────────────────────────────────────────────────────┐");
log("  │  Model                    k     AIC      dAIC     BIC      dBIC  │");
log("  ├────────────────────────────────────────────────────────────────────┤");
for(const r of icResults){
  log("  │  "+r.name.padEnd(25)+r.k.toString().padStart(3)+
    r.AIC.toFixed(0).padStart(8)+(r.AIC-minAIC).toFixed(0).padStart(8)+
    r.BIC.toFixed(0).padStart(9)+(r.BIC-minBIC).toFixed(0).padStart(8)+"  │");
}
log("  └────────────────────────────────────────────────────────────────────┘");
log("");

log("  FULL-SAMPLE REGRESSION COEFFICIENTS");
sep();
log("");

const bestParsimonious=models.find(m=>m.features.length>0&&m.features[0]!=='__pergalaxy__'&&m.features.length<=5);
if(bestParsimonious){
  const feats=['logMHI','inc','logD','T','Rdisk'];
  const X=galaxyData.map(g=>feats.map(f=>getFeature(g,f)));
  const reg=linReg(X,allLogA0);
  log("  Full model: log(a0,i) = "+reg.intercept.toFixed(4));
  for(let j=0;j<feats.length;j++){
    log("    "+((reg.coefs[j]>=0?"+":"")+reg.coefs[j].toFixed(4)).padStart(10)+" × "+feats[j]);
  }
  log("");

  let resSS=0,totSS=0;
  const my=allLogA0.reduce((s,v)=>s+v,0)/N;
  for(let i=0;i<N;i++){
    const pred=reg.intercept+reg.coefs.reduce((s,c,j)=>s+c*X[i][j],0);
    resSS+=(allLogA0[i]-pred)**2;
    totSS+=(allLogA0[i]-my)**2;
  }
  const R2=1-resSS/totSS;
  const adjR2=1-(resSS/(N-feats.length-1))/(totSS/(N-1));
  const residTau=Math.sqrt(resSS/(N-feats.length-1));

  log("  R² = "+R2.toFixed(3)+", adj-R² = "+adjR2.toFixed(3));
  log("  Residual SD = "+residTau.toFixed(4)+" dex");
  log("  Original tau = "+baseline.tau.toFixed(4)+" dex");
  log("  Tau reduction = "+((1-residTau/baseline.tau)*100).toFixed(1)+"%");
  log("");
}

log("  DECISION MAP");
sep();
log("");
log("  ┌──────────────────────────────────────────────────────────────────────────┐");
log("  │                                                                        │");
log("  │  WHEN DO WE SAY a0 IS ~CONSTANT?                                       │");
log("  │  ─────────────────────────────────                                      │");
log("  │  IF a simple model (M1-M4) closes >80% of the M0→M6 gap,              │");
log("  │  THEN a0 is effectively constant + known corrections.                   │");
log("  │                                                                        │");
log("  │  WHEN DO WE SAY a0 GENUINELY VARIES?                                   │");
log("  │  ────────────────────────────────────                                   │");
log("  │  IF even the best parametric model (M5) leaves >50% of the gap open,   │");
log("  │  THEN galaxies have genuine, unpredictable a0 variation.               │");
log("  │                                                                        │");
log("  └──────────────────────────────────────────────────────────────────────────┘");
log("");

const bestParam=cvResults.filter(r=>r.features[0]!=='__pergalaxy__').sort((a,b)=>a.cvRMS-b.cvRMS)[0];
const bestParamGap=((m0rms-bestParam.cvRMS)/gap*100);
const kitchenSink=cvResults.find(r=>r.name.includes('Kitchen'));
const kitchenGap=kitchenSink?((m0rms-kitchenSink.cvRMS)/gap*100):0;

log("  RESULTS:");
log("  Best parametric model: "+bestParam.name);
log("  Gap closed by best parametric: "+bestParamGap.toFixed(1)+"%");
log("  Gap closed by kitchen sink:    "+kitchenGap.toFixed(1)+"%");
log("  Gap remaining (unexplained):   "+(100-kitchenGap).toFixed(1)+"%");
log("");

if(kitchenGap>80){
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  VERDICT: a0 IS EFFECTIVELY CONSTANT.                               ║");
  log("  ║  Known galaxy properties explain >80% of the apparent variation.     ║");
  log("  ║  The \"variation\" is mostly correctable systematics.                  ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else if(kitchenGap<50){
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  VERDICT: a0 GENUINELY VARIES BETWEEN GALAXIES.                     ║");
  log("  ║  Even with all measured properties, >50% of variation remains.      ║");
  log("  ║  No known galaxy property can predict a galaxy's a0 adequately.     ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  VERDICT: INTERMEDIATE — partially predictable variation.           ║");
  log("  ║  Known properties explain "+kitchenGap.toFixed(0)+"% of the gap.                       ║");
  log("  ║  Significant unexplained variation remains.                         ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  HIERARCHY OF MODELS (from simplest to most complex):");
log("");
for(const r of cvResults){
  const gapPct=gap>0?((m0rms-r.cvRMS)/gap*100):0;
  const barLen=Math.round(gapPct/2);
  const bar="\u2588".repeat(Math.max(0,barLen));
  log("  "+r.name.padEnd(26)+(gapPct.toFixed(1)+"%").padStart(6)+"  "+bar);
}
log("");
log("  "+"\u2500".repeat(35));
log("  0%        25%       50%       75%       100%");
log("  M0                                      M6(free)");
log("");

log("=".repeat(80));
log("  PHASE 15 FINAL ANSWER");
log("=".repeat(80));
log("");
log("  The question was: a0,i = a0 + f(Xi) — can f(Xi) replace M3?");
log("");
log("  f(Xi) with best predictors closes "+kitchenGap.toFixed(1)+"% of the gap to M3.");
log("  Remaining unexplained: "+(100-kitchenGap).toFixed(1)+"%.");
log("");
log("  This means:");
if(kitchenGap>80){
  log("  → a0 is ~constant. The apparent variation is mostly explained");
  log("    by MHI, distance, inclination, morphology, and disk size.");
  log("  → MOND's strict universality is SUPPORTED (with corrections).");
} else if(kitchenGap<50){
  log("  → a0 GENUINELY VARIES. No combination of measured properties");
  log("    can predict a galaxy's a0. The variation is either:");
  log("    (a) driven by an unmeasured variable, or");
  log("    (b) physically real — a0 differs between galaxies.");
  log("  → Strict universality is DISFAVORED by these data.");
} else {
  log("  → PARTIAL: some variation is predictable (mainly via MHI),");
  log("    but a large fraction remains unexplained.");
  log("  → The question is not fully resolved with current data.");
}
log("");
log("=".repeat(80));

const output={
  version:VERSION,
  timestamp:new Date().toISOString(),
  description:"Phase 15: Model Ladder — from M0 to M3",
  baseline:{a0:Math.round(globalA0),logA0:+globalLogA0.toFixed(4),tau:+baseline.tau.toFixed(4)},
  nGalaxies:N,
  cvResults:cvResults.map(r=>({name:r.name,k:r.k,cvRMS:+r.cvRMS.toFixed(5),gapClosed:+(gap>0?((m0rms-r.cvRMS)/gap*100):0).toFixed(1)})),
  icResults:icResults.map(r=>({name:r.name,k:r.k,AIC:+r.AIC.toFixed(1),BIC:+r.BIC.toFixed(1),dAIC:+(r.AIC-minAIC).toFixed(1),dBIC:+(r.BIC-minBIC).toFixed(1)})),
  gapClosedByBestParametric:+bestParamGap.toFixed(1),
  gapClosedByKitchenSink:+kitchenGap.toFixed(1),
  verdict:kitchenGap>80?"CONSTANT":kitchenGap<50?"VARIES":"INTERMEDIATE"
};

fs.writeFileSync(path.join(__dirname,'../public/phase15-model-ladder.json'),JSON.stringify(output,null,2));
log("\n  Results saved to public/phase15-model-ladder.json");
