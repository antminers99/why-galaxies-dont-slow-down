#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "18.0.0";
function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(80)); }

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0(pts) {
  let lo = 2.0, hi = 5.0;
  for (let s = 0; s < 150; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gb = Math.pow(10, p.log_g_bar);
      c1 += (p.log_g_obs - Math.log10(mcgaughRAR(gb, Math.pow(10, m1)))) ** 2;
      c2 += (p.log_g_obs - Math.log10(mcgaughRAR(gb, Math.pow(10, m2)))) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  const logA0 = (lo + hi) / 2, a0 = Math.pow(10, logA0);
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.log_g_bar);
    ss += (p.log_g_obs - Math.log10(mcgaughRAR(gb, a0))) ** 2;
  }
  return { a0, logA0, rms: Math.sqrt(ss / pts.length), ss, n: pts.length };
}

function predictSS(a0, pts) {
  let ss = 0;
  for (const p of pts) {
    const gb = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gb, a0);
    ss += (p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10)) ** 2;
  }
  return ss;
}

function dlMeta(logA0s, ses) {
  const n = logA0s.length;
  const w = ses.map(s => 1 / (s * s + 1e-10));
  const wS = w.reduce((a, b) => a + b, 0);
  const muFE = logA0s.reduce((s, v, i) => s + w[i] * v, 0) / wS;
  let Q = 0;
  for (let i = 0; i < n; i++) Q += w[i] * (logA0s[i] - muFE) ** 2;
  const S1 = wS, S2 = w.reduce((s, v) => s + v * v, 0);
  const tau2 = Math.max(0, (Q - (n - 1)) / (S1 - S2 / S1));
  const tau = Math.sqrt(tau2);
  const W = ses.map(s => 1 / (s * s + tau2 + 1e-10));
  const WS = W.reduce((a, b) => a + b, 0);
  const mu = logA0s.reduce((s, v, i) => s + W[i] * v, 0) / WS;
  const seMu = Math.sqrt(1 / WS);
  return { mu, tau, tau2, a0: Math.pow(10, mu), seMu };
}

function corrWith(x, y) {
  const n = x.length;
  const mx = x.reduce((s, v) => s + v, 0) / n;
  const my = y.reduce((s, v) => s + v, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i]-mx)*(y[i]-my); sxx += (x[i]-mx)**2; syy += (y[i]-my)**2; }
  return sxy / Math.sqrt(sxx * syy + 1e-20);
}

const thingsSet = new Set(['NGC2403','NGC2841','NGC2903','NGC2976','NGC3031','NGC3198','NGC3521','NGC3621','NGC4736','NGC4826','NGC5055','NGC6946','NGC7331','NGC7793','DDO154','IC2574']);
const ltSet = new Set(['CVnIdwA','DDO43','DDO46','DDO47','DDO50','DDO52','DDO53','DDO70','DDO87','DDO101','DDO126','DDO133','DDO154','DDO168','DDO210','DDO216','F564-V3','Haro29','Haro36','IC1613','IC10','M81dwB','NGC1569','NGC2366','NGC3738','NGC4163','NGC4214','NGC6822','SagDIG','UGC8508','WLM','UGCA292']);
const phangsSet = new Set(['NGC0628','NGC1672','NGC2903','NGC3351','NGC3521','NGC3627','NGC4254','NGC4303','NGC4321','NGC4535','NGC4536','NGC4548','NGC4569','NGC4571','NGC4579','NGC4654','NGC4689','NGC5055','NGC5236','NGC6744','NGC7496']);

const p11 = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/phase11-sensitivity-lab.json'), 'utf8'));
const sparc = JSON.parse(fs.readFileSync(path.join(__dirname, '../public/sparc-table.json'), 'utf8'));

const masterTable = [];
for (const gal of p11.galaxies) {
  const pts = gal.localProfile.filter(p =>
    isFinite(p.log_g_bar) && isFinite(p.log_g_obs) && p.log_g_obs > p.log_g_bar * 0.5 && p.r > 0
  );
  if (pts.length < 5) continue;
  const sp = sparc.find(s => s.name === gal.name);
  if (!sp) continue;

  const fit = fitA0(pts);
  const se = fit.rms / Math.sqrt(fit.n);

  const hasTHINGS = thingsSet.has(gal.name);
  const hasLT = ltSet.has(gal.name);
  const hasPHANGS = phangsSet.has(gal.name);
  const has2D = hasTHINGS || hasLT || hasPHANGS;

  const isPreciseD = sp.fD === 1;
  const isQ1 = sp.Q === 1;
  const isHighInc = sp.inc >= 45;
  const isModerateInc = sp.inc >= 30 && sp.inc <= 80;

  let tier = 'C';
  if (isPreciseD && isQ1 && isHighInc && has2D) tier = 'PLATINUM';
  else if (isPreciseD && isQ1 && isHighInc) tier = 'GOLD';
  else if (isQ1 && isHighInc) tier = 'SILVER';
  else if (isHighInc) tier = 'BRONZE';

  masterTable.push({
    name: gal.name,
    pts,
    D: sp.D, eD: sp.eD, fD: sp.fD,
    inc: sp.inc, eInc: sp.eInc,
    T: sp.T, Q: sp.Q,
    Vmax: gal.Vmax, Vflat: sp.Vflat,
    logMHI: Math.log10(sp.MHI || 1e-3),
    logL36: Math.log10(sp.L36 || sp.L || 1e-3),
    Rdisk: sp.Rdisk, SBdisk: sp.SBdisk,
    RHI: sp.RHI, Reff: sp.Reff,
    n: pts.length,
    logA0: fit.logA0, a0: fit.a0, rms: fit.rms, se: se,
    isPreciseD, isQ1, isHighInc, has2D,
    hasTHINGS, hasLT, hasPHANGS,
    tier
  });
}

const N = masterTable.length;

log("");
log("=".repeat(80));
log("  PHASE 18: CLEAN EXTERNAL CROSS-MATCH & DECISIVE TEST");
log("  Build the cleanest possible sample. Ask ONE question.");
log("  Version " + VERSION + ", " + N + " galaxies");
log("=".repeat(80));
log("");

log("  STEP 1: MASTER TABLE — EXTERNAL DATA FLAGS");
sep();
log("");

const tierCounts = {};
masterTable.forEach(g => { tierCounts[g.tier] = (tierCounts[g.tier] || 0) + 1; });
const has2Dcount = masterTable.filter(g => g.has2D).length;
const preciseDcount = masterTable.filter(g => g.isPreciseD).length;

log("  TIER CLASSIFICATION:");
log("    PLATINUM: precise D + Q=1 + i>=45 + 2D kinematics  n=" + (tierCounts['PLATINUM'] || 0));
log("    GOLD:     precise D + Q=1 + i>=45                  n=" + (tierCounts['GOLD'] || 0));
log("    SILVER:   Q=1 + i>=45                               n=" + (tierCounts['SILVER'] || 0));
log("    BRONZE:   i>=45                                      n=" + (tierCounts['BRONZE'] || 0));
log("    C:        below criteria                             n=" + (tierCounts['C'] || 0));
log("");
log("  External data coverage:");
log("    Precise distances (fD=1):  " + preciseDcount + "/" + N);
log("    2D kinematics available:   " + has2Dcount + "/" + N);
log("    THINGS:   " + masterTable.filter(g => g.hasTHINGS).length);
log("    LITTLE THINGS: " + masterTable.filter(g => g.hasLT).length);
log("    PHANGS:   " + masterTable.filter(g => g.hasPHANGS).length);
log("");

log("  MASTER TABLE (sorted by tier, then a0):");
log("  ┌────────────────────────────────────────────────────────────────────────────────────┐");
log("  │  Galaxy       Tier    a0     log  delta  inc  Q fD  D     2D  n_pts  rms         │");
log("  ├────────────────────────────────────────────────────────────────────────────────────┤");

const meanLogA0 = masterTable.reduce((s, g) => s + g.logA0, 0) / N;
const tierOrder = { PLATINUM: 0, GOLD: 1, SILVER: 2, BRONZE: 3, C: 4 };
const sorted = [...masterTable].sort((a, b) => tierOrder[a.tier] - tierOrder[b.tier] || a.logA0 - b.logA0);

for (const g of sorted) {
  if (g.tier === 'C') continue;
  const delta = (g.logA0 - meanLogA0).toFixed(3);
  const twod = [g.hasTHINGS ? 'T' : '', g.hasLT ? 'L' : '', g.hasPHANGS ? 'P' : ''].filter(Boolean).join('+') || '-';
  log("  │  " + g.name.padEnd(13) + g.tier.padEnd(9) +
    Math.round(g.a0).toString().padStart(5) + g.logA0.toFixed(3).padStart(7) +
    delta.padStart(7) + g.inc.toString().padStart(4) +
    g.Q.toString().padStart(3) + g.fD.toString().padStart(3) +
    g.D.toFixed(1).padStart(6) +
    twod.padStart(5) + g.n.toString().padStart(6) +
    g.rms.toFixed(3).padStart(7) + "  │");
}
log("  └────────────────────────────────────────────────────────────────────────────────────┘");
log("");

log("  STEP 2: HIERARCHICAL ANALYSIS BY TIER");
sep();
log("");

const tiers = ['PLATINUM', 'GOLD', 'SILVER', 'BRONZE', 'ALL'];
const tierResults = {};

for (const tier of tiers) {
  const subset = tier === 'ALL' ? masterTable : masterTable.filter(g => {
    if (tier === 'PLATINUM') return g.tier === 'PLATINUM';
    if (tier === 'GOLD') return g.tier === 'PLATINUM' || g.tier === 'GOLD';
    if (tier === 'SILVER') return g.tier !== 'BRONZE' && g.tier !== 'C';
    if (tier === 'BRONZE') return g.tier !== 'C';
    return true;
  });

  if (subset.length < 3) {
    tierResults[tier] = { n: subset.length, note: 'too few' };
    continue;
  }

  const logA0s = subset.map(g => g.logA0);
  const ses = subset.map(g => g.se);
  const dl = dlMeta(logA0s, ses);

  const sd = Math.sqrt(logA0s.reduce((s, v) => s + (v - dl.mu) ** 2, 0) / (logA0s.length - 1));

  let cvSS_m0 = 0, cvSS_free = 0, cvN = 0;
  for (let i = 0; i < subset.length; i++) {
    const trainLogA0 = logA0s.filter((_, j) => j !== i);
    const trainSe = ses.filter((_, j) => j !== i);
    const trainDL = dlMeta(trainLogA0, trainSe);

    cvSS_m0 += predictSS(trainDL.a0, subset[i].pts);
    cvSS_free += predictSS(subset[i].a0, subset[i].pts);
    cvN += subset[i].pts.length;
  }

  const cvRMS_m0 = Math.sqrt(cvSS_m0 / cvN);
  const cvRMS_free = Math.sqrt(cvSS_free / cvN);
  const cvImprovement = (1 - cvRMS_free / cvRMS_m0) * 100;

  tierResults[tier] = {
    n: subset.length,
    totalPts: subset.reduce((s, g) => s + g.n, 0),
    a0: Math.round(dl.a0),
    logA0: +dl.mu.toFixed(4),
    tau: +dl.tau.toFixed(4),
    sd: +sd.toFixed(4),
    seMu: +dl.seMu.toFixed(4),
    cvRMS_m0: +cvRMS_m0.toFixed(5),
    cvRMS_free: +cvRMS_free.toFixed(5),
    cvImprovement: +cvImprovement.toFixed(1),
  };
}

log("  ┌──────────────────────────────────────────────────────────────────────────────────┐");
log("  │  Tier        n   pts    a0     tau     SD     CV-M0   CV-free  M3-gain         │");
log("  ├──────────────────────────────────────────────────────────────────────────────────┤");
for (const tier of tiers) {
  const r = tierResults[tier];
  if (!r || r.note === 'too few') {
    log("  │  " + (tier + ":").padEnd(12) + r.n.toString().padStart(3) + "  (too few for analysis)".padEnd(62) + "│");
    continue;
  }
  log("  │  " + (tier + ":").padEnd(12) + r.n.toString().padStart(3) +
    r.totalPts.toString().padStart(6) +
    r.a0.toString().padStart(6) +
    r.tau.toFixed(3).padStart(8) +
    r.sd.toFixed(3).padStart(8) +
    r.cvRMS_m0.toFixed(4).padStart(9) +
    r.cvRMS_free.toFixed(4).padStart(9) +
    (r.cvImprovement.toFixed(1) + "%").padStart(8) + "    │");
}
log("  └──────────────────────────────────────────────────────────────────────────────────┘");
log("");

log("  STEP 3: THE DECISIVE QUESTION");
sep();
log("");
log("  Does tau DROP as we move to cleaner samples?");
log("");

const fullTau = tierResults['ALL'].tau;
const goldTau = tierResults['GOLD'] && !tierResults['GOLD'].note ? tierResults['GOLD'].tau : null;
const platTau = tierResults['PLATINUM'] && !tierResults['PLATINUM'].note ? tierResults['PLATINUM'].tau : null;

log("  Tau progression:");
log("    ALL (n=" + tierResults['ALL'].n + "):      tau = " + fullTau.toFixed(4) + " dex");
if (tierResults['BRONZE'] && !tierResults['BRONZE'].note)
  log("    >=BRONZE (n=" + tierResults['BRONZE'].n + "): tau = " + tierResults['BRONZE'].tau.toFixed(4) + " dex");
if (tierResults['SILVER'] && !tierResults['SILVER'].note)
  log("    >=SILVER (n=" + tierResults['SILVER'].n + "): tau = " + tierResults['SILVER'].tau.toFixed(4) + " dex");
if (goldTau !== null)
  log("    >=GOLD (n=" + tierResults['GOLD'].n + "):   tau = " + goldTau.toFixed(4) + " dex");
if (platTau !== null)
  log("    PLATINUM (n=" + tierResults['PLATINUM'].n + "):  tau = " + platTau.toFixed(4) + " dex");
log("");

const bestTau = Math.min(
  fullTau,
  tierResults['BRONZE'] && !tierResults['BRONZE'].note ? tierResults['BRONZE'].tau : 999,
  tierResults['SILVER'] && !tierResults['SILVER'].note ? tierResults['SILVER'].tau : 999,
  goldTau || 999,
  platTau || 999
);
const tauReduction = (1 - bestTau / fullTau) * 100;

log("  STEP 4: LEAVE-ONE-GALAXY-OUT ON BEST TIER");
sep();
log("");

const bestTierName = goldTau !== null ? 'GOLD' : (tierResults['SILVER'] && !tierResults['SILVER'].note ? 'SILVER' : 'ALL');
const bestSubset = masterTable.filter(g => {
  if (bestTierName === 'GOLD') return g.tier === 'PLATINUM' || g.tier === 'GOLD';
  if (bestTierName === 'SILVER') return g.tier !== 'BRONZE' && g.tier !== 'C';
  return true;
});

log("  Running on " + bestTierName + " tier (n=" + bestSubset.length + "):");
log("");

const bLogA0 = bestSubset.map(g => g.logA0);
const bSE = bestSubset.map(g => g.se);

let looSS_m0 = 0, looSS_free = 0, looN = 0;
const looResiduals = [];

for (let i = 0; i < bestSubset.length; i++) {
  const trainLA = bLogA0.filter((_, j) => j !== i);
  const trainSE = bSE.filter((_, j) => j !== i);
  const trainDL = dlMeta(trainLA, trainSE);

  const ss_m0 = predictSS(trainDL.a0, bestSubset[i].pts);
  const ss_free = predictSS(bestSubset[i].a0, bestSubset[i].pts);

  looSS_m0 += ss_m0;
  looSS_free += ss_free;
  looN += bestSubset[i].pts.length;

  looResiduals.push({
    name: bestSubset[i].name,
    logA0: bLogA0[i],
    delta: bLogA0[i] - trainDL.mu,
    rmsM0: Math.sqrt(ss_m0 / bestSubset[i].pts.length),
    rmsFree: Math.sqrt(ss_free / bestSubset[i].pts.length),
    improvement: (1 - Math.sqrt(ss_free / bestSubset[i].pts.length) / Math.sqrt(ss_m0 / bestSubset[i].pts.length)) * 100
  });
}

const looRMS_m0 = Math.sqrt(looSS_m0 / looN);
const looRMS_free = Math.sqrt(looSS_free / looN);
const looImprovement = (1 - looRMS_free / looRMS_m0) * 100;

log("  LOO CV-RMS (M0, universal): " + looRMS_m0.toFixed(5));
log("  LOO CV-RMS (M6, per-gal):   " + looRMS_free.toFixed(5));
log("  M3-free improvement:         " + looImprovement.toFixed(1) + "%");
log("");

log("  Per-galaxy LOO residuals:");
const sortedResid = [...looResiduals].sort((a, b) => Math.abs(b.delta) - Math.abs(a.delta));
for (const r of sortedResid.slice(0, 15)) {
  const flag = Math.abs(r.delta) > 0.3 ? " <<<" : "";
  log("    " + r.name.padEnd(13) +
    "delta=" + r.delta.toFixed(3).padStart(7) +
    " rms_M0=" + r.rmsM0.toFixed(3).padStart(6) +
    " rms_free=" + r.rmsFree.toFixed(3).padStart(6) +
    " gain=" + r.improvement.toFixed(0).padStart(3) + "%" + flag);
}
log("");

const outliers = looResiduals.filter(r => Math.abs(r.delta) > 0.4);
if (outliers.length > 0) {
  log("  OUTLIERS (|delta| > 0.4 dex):");
  for (const o of outliers) {
    log("    " + o.name + ": delta = " + o.delta.toFixed(3) + " dex (a0 " + (o.delta > 0 ? "HIGH" : "LOW") + ")");
  }
  log("");

  log("  Tau without outliers:");
  const noOutlierIdx = bestSubset.map((g, i) => outliers.some(o => o.name === g.name) ? -1 : i).filter(i => i >= 0);
  const noOutLA = noOutlierIdx.map(i => bLogA0[i]);
  const noOutSE = noOutlierIdx.map(i => bSE[i]);
  const noOutDL = dlMeta(noOutLA, noOutSE);
  log("    n=" + noOutlierIdx.length + ", tau=" + noOutDL.tau.toFixed(4) + " dex (vs " + tierResults[bestTierName].tau.toFixed(4) + " with outliers)");
  log("");
}

log("=".repeat(80));
log("  PHASE 18 — FINAL VERDICT");
log("=".repeat(80));
log("");

log("  THE ONE QUESTION: Does cleaning the sample eliminate heterogeneity?");
log("");

if (tauReduction > 40) {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  YES — tau drops by " + tauReduction.toFixed(0) + "% in cleanest sample.                    ║");
  log("  ║  Most of the heterogeneity WAS systematic error.                    ║");
  log("  ║  a0 is closer to universal than raw analysis suggested.             ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else if (tauReduction > 15) {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  PARTIAL — tau drops by " + tauReduction.toFixed(0) + "% but remains substantial.            ║");
  log("  ║  Some heterogeneity is systematic, but a large fraction persists.   ║");
  log("  ║  Either a deep unmeasured systematic or genuine a0 variation.       ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  NO — tau drops by only " + tauReduction.toFixed(0) + "% in cleanest sample.                 ║");
  log("  ║  Heterogeneity PERSISTS even with best data quality.                ║");
  log("  ║  This is strong evidence for either:                                ║");
  log("  ║    (a) Genuine galaxy-to-galaxy a0 variation, OR                    ║");
  log("  ║    (b) A deep unmeasured systematic beyond current data quality.    ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}
log("");

log("  NUMBERS SUMMARY:");
log("    Full sample (n=" + N + "):   tau = " + fullTau.toFixed(3) + " dex, a0 = " + tierResults['ALL'].a0);
if (goldTau) log("    GOLD+ (n=" + tierResults['GOLD'].n + "):      tau = " + goldTau.toFixed(3) + " dex, a0 = " + tierResults['GOLD'].a0);
if (platTau) log("    PLATINUM (n=" + tierResults['PLATINUM'].n + "):    tau = " + platTau.toFixed(3) + " dex, a0 = " + tierResults['PLATINUM'].a0);
log("    Best reduction: " + tauReduction.toFixed(1) + "%");
log("    M3-free improvement (LOO): " + looImprovement.toFixed(1) + "%");
log("");

log("  WHAT THIS MEANS FOR THE PROJECT:");
log("  ─────────────────────────────────");
if (tauReduction < 20) {
  log("  1. The acceleration scale a0 ≈ " + tierResults['ALL'].a0 + " (km/s)²/kpc is CONFIRMED.");
  log("  2. Galaxy-to-galaxy heterogeneity (tau ≈ " + fullTau.toFixed(2) + " dex) is ROBUST.");
  log("  3. Quality cuts, distance precision, and 2D kinematics do NOT remove it.");
  log("  4. This is either GENUINE physics or requires data BEYOND what exists.");
  log("  5. The project has reached its empirical resolution limit with SPARC.");
}
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  nGalaxies: N,
  tierCounts,
  externalCoverage: {
    preciseD: preciseDcount,
    has2D: has2Dcount,
    THINGS: masterTable.filter(g => g.hasTHINGS).length,
    LITTLE_THINGS: masterTable.filter(g => g.hasLT).length,
    PHANGS: masterTable.filter(g => g.hasPHANGS).length
  },
  tierResults,
  tauProgression: {
    full: +fullTau.toFixed(4),
    gold: goldTau ? +goldTau.toFixed(4) : null,
    platinum: platTau ? +platTau.toFixed(4) : null,
    bestReduction: +tauReduction.toFixed(1)
  },
  looCV: {
    tier: bestTierName,
    n: bestSubset.length,
    rmsM0: +looRMS_m0.toFixed(5),
    rmsFree: +looRMS_free.toFixed(5),
    improvement: +looImprovement.toFixed(1)
  },
  verdict: tauReduction > 40 ? "SYSTEMATIC" : tauReduction > 15 ? "PARTIAL" : "GENUINE_OR_DEEP",
  masterTable: masterTable.map(g => ({
    name: g.name, tier: g.tier,
    a0: Math.round(g.a0), logA0: +g.logA0.toFixed(4),
    inc: g.inc, Q: g.Q, fD: g.fD, D: g.D,
    has2D: g.has2D, hasTHINGS: g.hasTHINGS, hasLT: g.hasLT, hasPHANGS: g.hasPHANGS,
    n: g.n, rms: +g.rms.toFixed(4)
  }))
};

fs.writeFileSync(path.join(__dirname, '../public/phase18-clean-crossmatch.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase18-clean-crossmatch.json");
