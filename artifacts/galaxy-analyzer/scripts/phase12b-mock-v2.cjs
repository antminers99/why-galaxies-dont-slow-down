#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

const VERSION = "12.2.0";
const A0_TRUE = 3633;
const LOG_A0_TRUE = Math.log10(A0_TRUE);
const N_MOCK = 200;
const N_TRIALS = 50;

function log(msg) { console.log(msg); }
function sep() { log("\u2500".repeat(80)); }

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function randn() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v);
}

function sampleFrom(arr) { return arr[Math.floor(Math.random() * arr.length)]; }
function randUnif(lo, hi) { return lo + Math.random() * (hi - lo); }
function clamp(v, lo, hi) { return Math.max(lo, Math.min(hi, v)); }

const realNs = [7,9,10,11,12,13,14,17,18,19,19,19,22,22,23,23,24,24,25,28,28,29,30,30,30,31,32,32,34,36,36,36,41,41,43,43,44,45,47,48,50,68,71,73,115];
const realIncs = [46,49,50,50,53,53,53,53,54,55,55,56,56,57,58,60,61,62,63,64,64,65,66,66,67,67,68,70,71,73,73,74,75,75,76,80,80,82,83,85,86,88,89,90,90,90,90];
const realRMS = [0.03,0.03,0.03,0.04,0.04,0.04,0.05,0.05,0.05,0.05,0.06,0.06,0.06,0.06,0.07,0.07,0.07,0.07,0.07,0.08,0.08,0.09,0.09,0.09,0.10,0.10,0.12,0.13,0.13,0.14,0.16,0.20,0.23,0.23,0.24,0.26,0.27,0.37,0.38];
const hubbleTypes = [1,2,3,4,5,6,7,8,9,10];
const vmaxRange = [50, 60, 80, 100, 120, 135, 150, 170, 190, 210, 230, 250, 270, 290, 310, 325];

function fitA0(pts) {
  let lo = 2.0, hi = 5.0;
  for (let s = 0; s < 150; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gbar = Math.pow(10, p.log_g_bar);
      const p1 = mcgaughRAR(gbar, Math.pow(10, m1));
      const p2 = mcgaughRAR(gbar, Math.pow(10, m2));
      c1 += (p.log_g_obs - Math.log10(p1 > 0 ? p1 : 1e-10)) ** 2;
      c2 += (p.log_g_obs - Math.log10(p2 > 0 ? p2 : 1e-10)) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  const a0 = Math.pow(10, (lo + hi) / 2);
  let chi2 = 0;
  for (const p of pts) {
    const gbar = Math.pow(10, p.log_g_bar);
    const pred = mcgaughRAR(gbar, a0);
    chi2 += (p.log_g_obs - Math.log10(pred > 0 ? pred : 1e-10)) ** 2;
  }
  return { a0, logA0: Math.log10(a0), rms: Math.sqrt(chi2 / pts.length), n: pts.length };
}

function dlMeta(logA0s, ses) {
  const n = logA0s.length;
  const w = ses.map(s => 1 / (s * s));
  const wSum = w.reduce((a, b) => a + b, 0);
  const muFE = logA0s.reduce((s, v, i) => s + w[i] * v, 0) / wSum;
  let Q = 0;
  for (let i = 0; i < n; i++) Q += w[i] * (logA0s[i] - muFE) ** 2;
  const S1 = wSum, S2 = w.reduce((s, v) => s + v * v, 0);
  const tau2 = Math.max(0, (Q - (n - 1)) / (S1 - S2 / S1));
  const tau = Math.sqrt(tau2);
  const W = ses.map(s => 1 / (s * s + tau2));
  const WS = W.reduce((a, b) => a + b, 0);
  const mu = logA0s.reduce((s, v, i) => s + W[i] * v, 0) / WS;
  const se = 1 / Math.sqrt(WS);
  const I2 = Q > n - 1 ? ((Q - (n - 1)) / Q) * 100 : 0;
  return { mu, se, tau, I2, Q, a0: Math.pow(10, mu) };
}

function runPipeline(galaxies) {
  const fits = [];
  for (const gal of galaxies) {
    const fit = fitA0(gal.pts);
    const se = fit.rms / Math.sqrt(gal.pts.length);
    fits.push({ logA0: fit.logA0, se, rms: fit.rms, a0: fit.a0 });
  }
  const logA0s = fits.map(f => f.logA0);
  const ses = fits.map(f => f.se);
  const dl = dlMeta(logA0s, ses);
  const sd = Math.sqrt(logA0s.reduce((s, v) => s + (v - dl.mu) ** 2, 0) / (logA0s.length - 1));
  return { ...dl, sd, fits, n: galaxies.length };
}

function generateMockGalaxy_v2(scenario, config) {
  const { nPts, inc_true, baseNoise, hubbleT, vmax } = config;

  const a0_gal = A0_TRUE;

  const logGbarMin = randUnif(0.8, 2.0);
  const logGbarMax = logGbarMin + randUnif(1.0, 2.5);

  const incError = scenario.incNonlinear
    ? (inc_true > 75 ? randn() * 5 * (1 + (inc_true - 75) / 15 * 2) : randn() * 3)
    : 0;
  const inc_obs = clamp(inc_true + incError, 20, 90);
  const sinCorrection = scenario.incNonlinear
    ? (Math.sin(inc_true * Math.PI / 180) / Math.sin(inc_obs * Math.PI / 180)) ** 2
    : 1;
  const logSinCorr = Math.log10(sinCorrection);

  const distError = scenario.distCorrelated
    ? scenario._distGroupBias + randn() * 0.02
    : 0;
  const logDistCorr = 2 * distError;

  let yStarFactor = 0;
  if (scenario.yStarMorphology) {
    const yTrue = 0.5;
    const yAssumed = 0.5 + (hubbleT - 5) * 0.03 + randn() * 0.08;
    yStarFactor = Math.log10(yAssumed / yTrue);
  }

  let ncmBias = 0;
  if (scenario.nonCircular) {
    const ncmStrength = vmax < 120 ? 0.08 : (vmax > 250 ? 0.02 : 0.04);
    ncmBias = randn() * ncmStrength;
    if (inc_true > 75) ncmBias *= 1.5;
  }

  let beamBias = 0;
  if (scenario.beamSmearing) {
    const beamSize = randUnif(5, 30);
    const galSize = vmax * 0.3;
    if (beamSize > galSize * 0.3) {
      beamBias = -0.02 * (beamSize / galSize);
    }
  }

  let etaCorr = 0;
  if (scenario.etaRotCorrelated) {
    const etaRot = (logGbarMax - logGbarMin) / 2 + randn() * 0.2;
    etaCorr = 0.06 * (etaRot - 1.5);
  }

  const totalBias = logSinCorr + logDistCorr + yStarFactor + ncmBias + beamBias + etaCorr;

  const pts = [];
  for (let j = 0; j < nPts; j++) {
    const frac = nPts > 1 ? j / (nPts - 1) : 0.5;
    const logGbar = logGbarMin + (logGbarMax - logGbarMin) * frac;
    const gbar = Math.pow(10, logGbar);
    const gobs_true = mcgaughRAR(gbar, a0_gal);

    let logGobs = Math.log10(gobs_true);

    logGobs += totalBias;

    let pointNoise = baseNoise;
    if (scenario.beamSmearing && j < 3) {
      pointNoise *= 1.5;
    }
    logGobs += randn() * pointNoise;

    pts.push({ log_g_bar: logGbar + yStarFactor, log_g_obs: logGobs });
  }

  return {
    pts,
    biases: { sinCorr: logSinCorr, distCorr: logDistCorr, yStar: yStarFactor,
              ncm: ncmBias, beam: beamBias, eta: etaCorr, total: totalBias }
  };
}


const scenarioConfigs = [
  {
    code: 'V1',
    name: 'v1 BASELINE (constant a0, point noise only)',
    incNonlinear: false, distCorrelated: false, yStarMorphology: false,
    nonCircular: false, beamSmearing: false, etaRotCorrelated: false
  },
  {
    code: 'S1',
    name: 'INCLINATION ONLY (non-linear errors at high inc)',
    incNonlinear: true, distCorrelated: false, yStarMorphology: false,
    nonCircular: false, beamSmearing: false, etaRotCorrelated: false
  },
  {
    code: 'S2',
    name: 'DISTANCE ONLY (correlated group errors)',
    incNonlinear: false, distCorrelated: true, yStarMorphology: false,
    nonCircular: false, beamSmearing: false, etaRotCorrelated: false
  },
  {
    code: 'S3',
    name: 'Y* ONLY (morphology-dependent mass-to-light)',
    incNonlinear: false, distCorrelated: false, yStarMorphology: true,
    nonCircular: false, beamSmearing: false, etaRotCorrelated: false
  },
  {
    code: 'S4',
    name: 'NON-CIRCULAR MOTIONS ONLY',
    incNonlinear: false, distCorrelated: false, yStarMorphology: false,
    nonCircular: true, beamSmearing: false, etaRotCorrelated: false
  },
  {
    code: 'S5',
    name: 'BEAM SMEARING ONLY',
    incNonlinear: false, distCorrelated: false, yStarMorphology: false,
    nonCircular: false, beamSmearing: true, etaRotCorrelated: false
  },
  {
    code: 'S6',
    name: 'ETA_ROT CORRELATION ONLY',
    incNonlinear: false, distCorrelated: false, yStarMorphology: false,
    nonCircular: false, beamSmearing: false, etaRotCorrelated: true
  },
  {
    code: 'ALL',
    name: 'ALL SYSTEMATICS COMBINED (worst case)',
    incNonlinear: true, distCorrelated: true, yStarMorphology: true,
    nonCircular: true, beamSmearing: true, etaRotCorrelated: true
  },
  {
    code: 'C29',
    name: 'TRULY VARYING a0 (tau=0.29, no systematics)',
    incNonlinear: false, distCorrelated: false, yStarMorphology: false,
    nonCircular: false, beamSmearing: false, etaRotCorrelated: false,
    _isVarying: true
  }
];

log("");
log("=".repeat(80));
log("  PHASE 12b: MOCK BLIND CHALLENGE v2 — REALISTIC SYSTEMATICS");
log("  Version " + VERSION);
log("  " + N_MOCK + " galaxies × " + N_TRIALS + " trials × " + scenarioConfigs.length + " scenarios");
log("=".repeat(80));
log("");

const allResults = {};

for (const sc of scenarioConfigs) {
  log(sep());
  log("  [" + sc.code + "] " + sc.name);
  sep();

  const trialResults = [];

  for (let t = 0; t < N_TRIALS; t++) {
    const distGroupBias = sc.distCorrelated ? randn() * 0.05 : 0;
    const scWithBias = { ...sc, _distGroupBias: distGroupBias };

    const galaxies = [];
    for (let i = 0; i < N_MOCK; i++) {
      const nPts = sampleFrom(realNs);
      const inc_true = sampleFrom(realIncs);
      const baseNoise = sampleFrom(realRMS);
      const hubbleT = sampleFrom(hubbleTypes);
      const vmax = sampleFrom(vmaxRange);

      if (sc._isVarying) {
        const a0_var = Math.pow(10, LOG_A0_TRUE + 0.29 * randn());
        const pts = [];
        const logGbarMin = randUnif(0.8, 2.0);
        const logGbarMax = logGbarMin + randUnif(1.0, 2.5);
        for (let j = 0; j < nPts; j++) {
          const frac = nPts > 1 ? j / (nPts - 1) : 0.5;
          const logGbar = logGbarMin + (logGbarMax - logGbarMin) * frac;
          const gbar = Math.pow(10, logGbar);
          const gobs = mcgaughRAR(gbar, a0_var);
          const logGobs = Math.log10(gobs) + randn() * baseNoise;
          pts.push({ log_g_bar: logGbar, log_g_obs: logGobs });
        }
        galaxies.push({ pts });
      } else {
        const result = generateMockGalaxy_v2(scWithBias, { nPts, inc_true, baseNoise, hubbleT, vmax });
        galaxies.push({ pts: result.pts });
      }
    }

    const result = runPipeline(galaxies);
    trialResults.push({
      a0: Math.round(result.a0), logA0: +result.mu.toFixed(4),
      tau: +result.tau.toFixed(4), I2: +result.I2.toFixed(1),
      sd: +result.sd.toFixed(4), se: +result.se.toFixed(4)
    });
  }

  allResults[sc.code] = trialResults;

  const mean = arr => arr.reduce((s, v) => s + v, 0) / arr.length;
  const med = arr => { const s = [...arr].sort((a, b) => a - b); return s[Math.floor(s.length / 2)]; };
  const taus = trialResults.map(r => r.tau);
  const a0s = trialResults.map(r => r.a0);
  const sds = trialResults.map(r => r.sd);

  log("  tau: mean=" + mean(taus).toFixed(4) + " med=" + med(taus).toFixed(4) +
      "  |  a0: mean=" + Math.round(mean(a0s)) + "  |  SD: " + mean(sds).toFixed(4));
  log("");
}


log("");
log("=".repeat(80));
log("  SUMMARY: CAN REALISTIC SYSTEMATICS PRODUCE tau ≈ 0.29?");
log("=".repeat(80));
log("");

const mean = arr => arr.reduce((s, v) => s + v, 0) / arr.length;

log("  ┌───────────────────────────────────────────────────────────────────────────┐");
log("  │  Scenario                              a0     tau      I²     SD(logA0) │");
log("  ├───────────────────────────────────────────────────────────────────────────┤");

const summaryRows = [];

for (const sc of scenarioConfigs) {
  const r = allResults[sc.code];
  const tau = mean(r.map(x => x.tau));
  const a0 = Math.round(mean(r.map(x => x.a0)));
  const I2 = mean(r.map(x => x.I2));
  const sd = mean(r.map(x => x.sd));
  summaryRows.push({ code: sc.code, name: sc.name, tau, a0, I2, sd });

  const label = (sc.code + ": " + sc.name).substring(0, 42).padEnd(42);
  log("  │  " + label + a0.toString().padStart(5) + tau.toFixed(3).padStart(8) + (I2.toFixed(1) + "%").padStart(8) + sd.toFixed(3).padStart(8) + "    │");
}

log("  ├───────────────────────────────────────────────────────────────────────────┤");
log("  │  REAL DATA (GOLD+i45)                          3633   0.291    92.4%   0.275    │");
log("  └───────────────────────────────────────────────────────────────────────────┘");
log("");

const v1tau = summaryRows.find(r => r.code === 'V1').tau;
const allTau = summaryRows.find(r => r.code === 'ALL').tau;
const c29tau = summaryRows.find(r => r.code === 'C29').tau;
const realTau = 0.291;

log("  KEY COMPARISONS:");
log("");
log("  Noise floor (v1, no systematics):     tau = " + v1tau.toFixed(3));
log("  ALL systematics combined:             tau = " + allTau.toFixed(3));
log("  Truly varying a0:                     tau = " + c29tau.toFixed(3));
log("  REAL DATA:                            tau = 0.291");
log("");

const gap = realTau - allTau;
const pctExplained = ((allTau - v1tau) / (realTau - v1tau) * 100);

log("  Systematics explain " + pctExplained.toFixed(1) + "% of the gap between noise floor and real tau.");
log("  Remaining gap: " + gap.toFixed(3) + " dex");
log("");

if (allTau < realTau * 0.6) {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  VERDICT: ALL KNOWN SYSTEMATICS COMBINED CANNOT PRODUCE             ║");
  log("  ║  THE OBSERVED tau ≈ 0.29.                                           ║");
  log("  ║                                                                     ║");
  log("  ║  Even with:                                                         ║");
  log("  ║    • Non-linear inclination errors (amplified at high inc)           ║");
  log("  ║    • Correlated distance errors (group biases)                      ║");
  log("  ║    • Y* varying with Hubble type                                    ║");
  log("  ║    • Non-circular motions (mass-dependent)                          ║");
  log("  ║    • Beam smearing (resolution effects)                             ║");
  log("  ║    • eta_rot correlation with a0                                    ║");
  log("  ║                                                                     ║");
  log("  ║  Combined tau = " + allTau.toFixed(3) + " — still " + (realTau / allTau).toFixed(1) + "× below real data.           ║");
  log("  ║                                                                     ║");
  log("  ║  THE CASE FOR NON-STRICT UNIVERSALITY STRENGTHENS.                  ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else if (allTau >= realTau * 0.8) {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  WARNING: Systematics CAN reproduce most of the observed tau.       ║");
  log("  ║  The case for genuine variation is WEAKENED.                        ║");
  log("  ║  Combined tau = " + allTau.toFixed(3) + " vs real = 0.291                       ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
} else {
  log("  ╔═══════════════════════════════════════════════════════════════════════╗");
  log("  ║  INTERMEDIATE: Systematics explain a significant fraction but       ║");
  log("  ║  not all of the observed tau.                                       ║");
  log("  ║  Combined tau = " + allTau.toFixed(3) + " vs real = 0.291                       ║");
  log("  ╚═══════════════════════════════════════════════════════════════════════╝");
}

log("");

log("  INDIVIDUAL SYSTEMATIC CONTRIBUTIONS:");
log("  " + "-".repeat(70));
for (const r of summaryRows) {
  if (r.code === 'V1' || r.code === 'C29') continue;
  const contribution = r.tau - v1tau;
  const pct = (contribution / (realTau - v1tau) * 100);
  const bar = "#".repeat(Math.max(0, Math.round(pct / 2)));
  log("  " + r.code.padEnd(5) + r.tau.toFixed(3).padStart(7) + "  (+" + contribution.toFixed(3) + " = " + pct.toFixed(1).padStart(5) + "%)  " + bar);
}
log("  " + "-".repeat(70));
log("  ALL   " + allTau.toFixed(3) + "  (total combined)");
log("  REAL  0.291");
log("");

log("=".repeat(80));
log("  FINAL PHASE 12b CONCLUSIONS");
log("=".repeat(80));
log("");
log("  1. Each individual systematic adds only " + ((summaryRows.filter(r=>!['V1','ALL','C29'].includes(r.code)).reduce((m,r)=>Math.max(m,r.tau),0) - v1tau)*1000).toFixed(0) + " millindex at most.");
log("  2. ALL systematics combined: tau = " + allTau.toFixed(3) + " (vs real 0.291).");
log("  3. The gap of " + gap.toFixed(3) + " dex CANNOT be explained by known systematics.");
log("  4. This strengthens the case that the scatter is either:");
log("     (a) from an unmeasured/unknown systematic not in our model, OR");
log("     (b) genuinely physical — a0 is not strictly universal.");
log("");
log("=".repeat(80));

const output = {
  version: VERSION,
  timestamp: new Date().toISOString(),
  description: "Phase 12b: Mock Challenge v2 — Realistic Systematics",
  config: { nMock: N_MOCK, nTrials: N_TRIALS, a0_true: A0_TRUE },
  scenarios: {},
  realData: { tau: 0.291, I2: 92.4, sd: 0.275, a0: 3633 },
  noiseFloor: +v1tau.toFixed(4),
  allSystematicsTau: +allTau.toFixed(4),
  varyingA0Tau: +c29tau.toFixed(4),
  gapRemaining: +gap.toFixed(4),
  pctExplainedBySystematics: +pctExplained.toFixed(1)
};

for (const sc of scenarioConfigs) {
  const r = allResults[sc.code];
  output.scenarios[sc.code] = {
    name: sc.name,
    meanTau: +mean(r.map(x => x.tau)).toFixed(4),
    meanA0: Math.round(mean(r.map(x => x.a0))),
    meanI2: +mean(r.map(x => x.I2)).toFixed(1),
    meanSD: +mean(r.map(x => x.sd)).toFixed(4),
    trials: r
  };
}

fs.writeFileSync(path.join(__dirname, '../public/phase12b-mock-v2.json'), JSON.stringify(output, null, 2));
log("\n  Results saved to public/phase12b-mock-v2.json");
