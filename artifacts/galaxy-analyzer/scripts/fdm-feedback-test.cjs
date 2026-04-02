const fs = require('fs');
const path = require('path');

const G = 4.3009e-6;
const UPSILON_D = 0.5, UPSILON_B = 0.7;
const N_BOOT = 5000;
const N_REALIZATIONS = 50;
const N_PER_REAL = 300;

function seededRNG(seed) {
  let s = seed;
  return () => { s = (s * 1664525 + 1013904223) & 0x7fffffff; return s / 0x7fffffff; };
}

function gaussRNG(rng, mu, sig) {
  const u1 = rng(), u2 = rng();
  return mu + sig * Math.sqrt(-2 * Math.log(u1 + 1e-10)) * Math.cos(2 * Math.PI * u2);
}

function linearRegression(x, y) {
  const n = x.length;
  if (n < 3) return { slope: NaN, intercept: NaN, r: NaN, r2: NaN, n };
  const mx = x.reduce((a, b) => a + b, 0) / n;
  const my = y.reduce((a, b) => a + b, 0) / n;
  let sxx = 0, sxy = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxx += (x[i] - mx) ** 2; sxy += (x[i] - mx) * (y[i] - my); syy += (y[i] - my) ** 2; }
  const slope = sxy / sxx;
  const intercept = my - slope * mx;
  const r2 = sxx > 0 && syy > 0 ? (sxy * sxy) / (sxx * syy) : 0;
  const r = Math.sign(slope) * Math.sqrt(Math.abs(r2));
  return { slope, intercept, r, r2, n };
}

function pearsonR(x, y) {
  const n = x.length;
  if (n < 3) return NaN;
  const mx = x.reduce((a, b) => a + b) / n;
  const my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxy / Math.sqrt(sxx * syy);
}

function bootstrapSlope(xs, ys, nBoot, rng) {
  const n = xs.length;
  const slopes = [];
  for (let i = 0; i < nBoot; i++) {
    const bx = [], by = [];
    for (let j = 0; j < n; j++) { const idx = Math.floor(rng() * n); bx.push(xs[idx]); by.push(ys[idx]); }
    const reg = linearRegression(bx, by);
    if (!isNaN(reg.slope)) slopes.push(reg.slope);
  }
  slopes.sort((a, b) => a - b);
  const mean = slopes.reduce((a, b) => a + b, 0) / slopes.length;
  const sd = Math.sqrt(slopes.reduce((s, v) => s + (v - mean) ** 2, 0) / (slopes.length - 1));
  return { mean, sd, ci95: [slopes[Math.floor(slopes.length * 0.025)], slopes[Math.floor(slopes.length * 0.975)]], slopes };
}

function nfwMenc(r, M200, c) {
  const rho_crit = 136.18;
  const R200 = Math.cbrt(M200 / (4 / 3 * Math.PI * 200 * rho_crit));
  const rs = R200 / c;
  const x = r / rs;
  const gc = Math.log(1 + c) - c / (1 + c);
  return M200 * (Math.log(1 + x) - x / (1 + x)) / gc;
}

function expDiskMenc(r, Mstar, Rd) {
  const x = r / Rd;
  return Mstar * (1 - (1 + x) * Math.exp(-x));
}

function gasDiskMenc(r, Mgas, Rg) {
  const x = r / Rg;
  return Mgas * (1 - (1 + x) * Math.exp(-x));
}

function diCintioAlpha(logMstarMhalo) {
  const x = logMstarMhalo;
  if (x < -4.1) return -1.0;
  if (x > -1.3) return -1.0;
  const peak = -2.7;
  const width = 0.7;
  const maxFlattening = 0.8;
  const gaussian = maxFlattening * Math.exp(-0.5 * ((x - peak) / width) ** 2);
  return -1.0 + gaussian;
}

function adiabaticContraction(r, M200, c, Mbar_enc) {
  const Mdm_init = nfwMenc(r, M200, c);
  const rho_crit = 136.18;
  const R200 = Math.cbrt(M200 / (4 / 3 * Math.PI * 200 * rho_crit));
  const fb = 0.157;
  const M_total_init = Mdm_init / (1 - fb);
  const ri = r;
  const rf_ratio = (M_total_init) / (Mdm_init + Mbar_enc);
  const rf = ri / Math.max(0.5, Math.min(2.0, rf_ratio));
  const Mdm_contracted = nfwMenc(rf, M200, c);
  return Mdm_contracted;
}

function coredNFWMenc(r, M200, c, coreAlpha) {
  const rho_crit = 136.18;
  const R200 = Math.cbrt(M200 / (4 / 3 * Math.PI * 200 * rho_crit));
  const rs = R200 / c;
  const rc = rs * 0.5;
  const x = r / rs;
  const gc = Math.log(1 + c) - c / (1 + c);
  const nfwM = M200 * (Math.log(1 + x) - x / (1 + x)) / gc;
  const coreFactor = Math.pow(1 + (rc / r) ** 2, coreAlpha / 2);
  return nfwM * Math.max(0.1, Math.min(1.5, coreFactor));
}

function generateGalaxy_noFeedback(rng) {
  const logM200 = 10 + rng() * 2.5;
  const M200 = Math.pow(10, logM200);
  const c = 10 * Math.pow(M200 / 1e12, -0.1) * (0.8 + 0.4 * rng());
  const fb = 0.05 * Math.pow(M200 / 1e12, 0.3);
  const fstar = 0.6 + 0.2 * rng();
  const Mstar = M200 * fb * fstar;
  const Mgas = M200 * fb * (1 - fstar) * 1.33;
  const Rd = 0.015 * Math.cbrt(M200 / 1e12) * (0.7 + 0.6 * rng()) * 1000;
  const Rg = Rd * (1.5 + rng());

  const radii = [];
  for (let j = 1; j <= 25; j++) radii.push(Rd * 0.2 * j);
  const rmax = radii[radii.length - 1];

  const pts = [];
  for (const r of radii) {
    const Mbar_enc = expDiskMenc(r, Mstar, Rd) + gasDiskMenc(r, Mgas, Rg);
    const Mdm_enc = nfwMenc(r, M200, c);
    const vObs = Math.sqrt(G * (Mbar_enc + Mdm_enc) / r);
    const vBar = Math.sqrt(G * Mbar_enc / r);
    const vDMsq = vObs * vObs - vBar * vBar;
    const fDM = vDMsq / (vObs * vObs);
    if (fDM < 0 || fDM > 1 || vObs <= 0) continue;
    pts.push({ r, vObs, vBar, vDMsq, vDM: vDMsq > 0 ? Math.sqrt(vDMsq) : 0 });
  }
  if (pts.length < 5) return null;

  const vmax = pts.reduce((mx, p) => Math.max(mx, p.vObs), 0);
  const sigBarGal = (Mstar + Mgas) / (Math.PI * Rd * Rd);
  return extractMetrics(pts, vmax, sigBarGal, rmax);
}

function generateGalaxy_withFeedback(rng) {
  const logM200 = 10 + rng() * 2.5;
  const M200 = Math.pow(10, logM200);
  const c_base = 10 * Math.pow(M200 / 1e12, -0.1);
  const c = c_base * (0.8 + 0.4 * rng());

  const fb_peak = 0.023;
  const logMpeak = 12.0;
  const fb = fb_peak * Math.exp(-0.5 * ((logM200 - logMpeak) / 0.6) ** 2) + 0.005;
  const fstar_base = 0.3 + 0.4 * Math.min(1, Math.pow(M200 / 1e11, 0.3));
  const fstar = Math.min(0.9, fstar_base * (0.8 + 0.4 * rng()));
  const Mstar = M200 * fb * fstar;
  const Mgas = M200 * fb * (1 - fstar) * 1.33;

  const lambda = 0.03 + 0.02 * gaussRNG(rng, 0, 1);
  const rho_crit = 136.18;
  const R200 = Math.cbrt(M200 / (4 / 3 * Math.PI * 200 * rho_crit));
  const Rd = Math.max(0.3, Math.abs(lambda) * R200 / Math.sqrt(2) * (0.8 + 0.4 * rng()));
  const Rg = Rd * (1.5 + rng());

  const logMstarMhalo = Math.log10(Mstar / M200);
  const innerAlpha = diCintioAlpha(logMstarMhalo);

  const radii = [];
  for (let j = 1; j <= 25; j++) radii.push(Rd * 0.2 * j);
  const rmax = radii[radii.length - 1];

  const pts = [];
  for (const r of radii) {
    const Mbar_enc = expDiskMenc(r, Mstar, Rd) + gasDiskMenc(r, Mgas, Rg);

    let Mdm_enc;
    if (innerAlpha > -0.8) {
      Mdm_enc = coredNFWMenc(r, M200, c, innerAlpha + 1);
    } else {
      Mdm_enc = adiabaticContraction(r, M200, c, Mbar_enc);
    }

    if (Mdm_enc <= 0) Mdm_enc = nfwMenc(r, M200, c) * 0.1;

    const vObs = Math.sqrt(G * (Mbar_enc + Mdm_enc) / r);
    const vBar = Math.sqrt(G * Mbar_enc / r);
    const vDMsq = vObs * vObs - vBar * vBar;
    const fDM = vDMsq / (vObs * vObs);
    if (fDM < 0 || fDM > 1 || vObs <= 0) continue;
    pts.push({ r, vObs, vBar, vDMsq, vDM: vDMsq > 0 ? Math.sqrt(vDMsq) : 0 });
  }
  if (pts.length < 5) return null;

  const vmax = pts.reduce((mx, p) => Math.max(mx, p.vObs), 0);
  const sigBarGal = (Mstar + Mgas) / (Math.PI * Rd * Rd);
  return extractMetrics(pts, vmax, sigBarGal, rmax);
}

function extractMetrics(pts, vmax, sigBarGal, rmax) {
  const logVmax = Math.log10(vmax);
  const logSigBar = Math.log10(sigBarGal);

  let rDMdom = NaN;
  for (const pt of pts) { if (pt.vDMsq > pt.vBar * pt.vBar && isNaN(rDMdom)) rDMdom = pt.r; }

  const innerPts = pts.filter(p => p.r / rmax < 0.3 && p.vDMsq > 0);
  let alpha = NaN;
  if (innerPts.length >= 3) {
    const reg = linearRegression(innerPts.map(p => Math.log10(p.r)), innerPts.map(p => Math.log10(Math.sqrt(p.vDMsq))));
    alpha = reg.slope;
  }

  const vdmPts = pts.filter(p => p.vDMsq > 0);
  const meanVDM = vdmPts.length > 0 ? vdmPts.reduce((s, p) => s + p.vDM, 0) / vdmPts.length : 0;

  return {
    vmax, logVmax, logSigBar,
    logMeanVDM: meanVDM > 0 ? Math.log10(meanVDM) : NaN,
    rDMdomNorm: isNaN(rDMdom) ? null : rDMdom / rmax,
    alpha: isNaN(alpha) ? null : alpha,
  };
}

const rotmodDir = '/tmp/rotmod';
const sparcTablePath = '/tmp/sparc_table.mrt';
const sparcTable = {};
const tableLines = fs.readFileSync(sparcTablePath, 'utf8').split('\n');
for (const line of tableLines) {
  const parts = line.trim().split(/\s+/);
  if (parts.length < 7) continue;
  const name = parts[0];
  if (name === 'Galaxy' || name.startsWith('#') || name.startsWith('-')) continue;
  sparcTable[name] = { dist: parseFloat(parts[2]), inc: parseFloat(parts[5]) };
}

const sparcReal = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const obsGalaxies = [];
const rotmodFiles = fs.readdirSync(rotmodDir).filter(f => f.endsWith('_rotmod.dat'));

for (const file of rotmodFiles) {
  const name = file.replace('_rotmod.dat', '');
  const info = sparcTable[name];
  if (!info) continue;
  const perGalaxy = sparcReal.perGalaxy?.find(g => g.name === name);
  if (!perGalaxy || perGalaxy.sigma_bar <= 0) continue;

  const lines = fs.readFileSync(path.join(rotmodDir, file), 'utf8').trim().split('\n');
  const dataLines = lines.filter(l => !l.startsWith('#') && l.trim());

  const pts = [];
  for (const line of dataLines) {
    const p = line.trim().split(/\s+/).map(Number);
    if (p.length < 7) continue;
    const [r, vObs, eV, vGas, vDisk, vBulge] = [p[0], p[1], p[2], p[3], p[4], p[5]];
    if (r <= 0 || vObs <= 0) continue;
    const vBarSq = vGas * Math.abs(vGas) + UPSILON_D * vDisk * Math.abs(vDisk) + UPSILON_B * vBulge * Math.abs(vBulge);
    if (vBarSq <= 0) continue;
    const fDM = (vObs * vObs - vBarSq) / (vObs * vObs);
    const vDMsq = vObs * vObs - vBarSq;
    if (fDM < 0 || fDM > 1) continue;
    pts.push({ r, vObs, vBar: Math.sqrt(vBarSq), vDMsq, vDM: vDMsq > 0 ? Math.sqrt(vDMsq) : 0 });
  }
  if (pts.length < 5) continue;

  const rmax = pts[pts.length - 1].r;
  const vmax = pts.reduce((mx, p) => Math.max(mx, p.vObs), 0);
  const logSigBar = Math.log10(perGalaxy.sigma_bar);

  let rDMdom = NaN;
  for (const pt of pts) { if (pt.vDMsq > pt.vBar * pt.vBar && isNaN(rDMdom)) rDMdom = pt.r; }

  const innerPts = pts.filter(p => p.r / rmax < 0.3 && p.vDMsq > 0);
  let alpha = NaN;
  if (innerPts.length >= 3) {
    const reg = linearRegression(innerPts.map(p => Math.log10(p.r)), innerPts.map(p => Math.log10(Math.sqrt(p.vDMsq))));
    alpha = reg.slope;
  }

  const vdmPts = pts.filter(p => p.vDMsq > 0);
  const meanVDM = vdmPts.length > 0 ? vdmPts.reduce((s, p) => s + p.vDM, 0) / vdmPts.length : 0;

  obsGalaxies.push({
    name, vmax, logVmax: Math.log10(vmax), logSigBar,
    logMeanVDM: meanVDM > 0 ? Math.log10(meanVDM) : NaN,
    rDMdomNorm: isNaN(rDMdom) ? null : rDMdom / rmax,
    alpha: isNaN(alpha) ? null : alpha,
  });
}

console.log(`\n${'═'.repeat(70)}`);
console.log(`  ΛCDM + BARYONIC FEEDBACK TEST`);
console.log(`  Does feedback explain the observed baryon-halo coupling?`);
console.log(`${'═'.repeat(70)}\n`);
console.log(`  Observed galaxies: ${obsGalaxies.length}`);

console.log(`\n  Generating ΛCDM WITHOUT feedback...`);
const noFeedbackGals = [];
for (let r = 0; r < N_REALIZATIONS; r++) {
  const rng = seededRNG(42 + r * 7919);
  for (let i = 0; i < N_PER_REAL; i++) {
    const g = generateGalaxy_noFeedback(rng);
    if (g) noFeedbackGals.push(g);
  }
}
console.log(`  No-feedback galaxies: ${noFeedbackGals.length}`);

console.log(`  Generating ΛCDM WITH feedback (core formation + adiabatic contraction)...`);
const feedbackGals = [];
for (let r = 0; r < N_REALIZATIONS; r++) {
  const rng = seededRNG(42 + r * 7919);
  for (let i = 0; i < N_PER_REAL; i++) {
    const g = generateGalaxy_withFeedback(rng);
    if (g) feedbackGals.push(g);
  }
}
console.log(`  Feedback galaxies: ${feedbackGals.length}`);

const METRICS = [
  { name: 'V_DM', xKey: 'logSigBar', yKey: 'logMeanVDM', filter: g => !isNaN(g.logMeanVDM) },
  { name: 'r_DMdom', xKey: 'logSigBar', yKey: 'rDMdomNorm', filter: g => g.rDMdomNorm !== null },
  { name: 'alpha', xKey: 'logSigBar', yKey: 'alpha', filter: g => g.alpha !== null },
];

let seedCounter = 200000;
const results = { metrics: [], feedbackModel: {}, summary: {} };

console.log(`\n${'═'.repeat(70)}`);
console.log(`  COMPARISON: Observed vs No-Feedback vs With-Feedback`);
console.log(`${'═'.repeat(70)}\n`);

for (const metric of METRICS) {
  const obsF = obsGalaxies.filter(metric.filter);
  const noFbF = noFeedbackGals.filter(metric.filter);
  const fbF = feedbackGals.filter(metric.filter);

  const obsXs = obsF.map(g => g[metric.xKey]);
  const obsYs = obsF.map(g => g[metric.yKey]);
  const noFbXs = noFbF.map(g => g[metric.xKey]);
  const noFbYs = noFbF.map(g => g[metric.yKey]);
  const fbXs = fbF.map(g => g[metric.xKey]);
  const fbYs = fbF.map(g => g[metric.yKey]);

  const obsReg = linearRegression(obsXs, obsYs);
  const noFbReg = linearRegression(noFbXs, noFbYs);
  const fbReg = linearRegression(fbXs, fbYs);

  const obsR = pearsonR(obsXs, obsYs);
  const noFbR = pearsonR(noFbXs, noFbYs);
  const fbR = pearsonR(fbXs, fbYs);

  const obsBoot = bootstrapSlope(obsXs, obsYs, N_BOOT, seededRNG(seedCounter++));
  const noFbBoot = bootstrapSlope(noFbXs, noFbYs, N_BOOT, seededRNG(seedCounter++));
  const fbBoot = bootstrapSlope(fbXs, fbYs, N_BOOT, seededRNG(seedCounter++));

  const deltaNoFb = obsBoot.mean - noFbBoot.mean;
  const seNoFb = Math.sqrt(obsBoot.sd ** 2 + noFbBoot.sd ** 2);
  const sigmaNoFb = seNoFb > 0 ? Math.abs(deltaNoFb) / seNoFb : 0;

  const deltaFb = obsBoot.mean - fbBoot.mean;
  const seFb = Math.sqrt(obsBoot.sd ** 2 + fbBoot.sd ** 2);
  const sigmaFb = seFb > 0 ? Math.abs(deltaFb) / seFb : 0;

  const feedbackExplainsPercent = deltaNoFb !== 0 ? ((sigmaNoFb - sigmaFb) / sigmaNoFb * 100) : 0;

  console.log(`  ${metric.name}:`);
  console.log(`    Observed:     slope = ${obsReg.slope.toFixed(5)}, r = ${obsR.toFixed(3)}, n = ${obsF.length}`);
  console.log(`    No feedback:  slope = ${noFbReg.slope.toFixed(5)}, r = ${noFbR.toFixed(3)}, n = ${noFbF.length}`);
  console.log(`    With feedback: slope = ${fbReg.slope.toFixed(5)}, r = ${fbR.toFixed(3)}, n = ${fbF.length}`);
  console.log(`    Excess (no fb):   Δ = ${deltaNoFb.toFixed(5)}, σ = ${sigmaNoFb.toFixed(1)}`);
  console.log(`    Excess (with fb): Δ = ${deltaFb.toFixed(5)}, σ = ${sigmaFb.toFixed(1)}`);
  console.log(`    Feedback explains: ${feedbackExplainsPercent.toFixed(1)}% of the tension`);
  console.log();

  results.metrics.push({
    name: metric.name,
    observed: { slope: +obsReg.slope.toFixed(6), r: +obsR.toFixed(4), n: obsF.length, bootMean: +obsBoot.mean.toFixed(6), bootSD: +obsBoot.sd.toFixed(6), ci95: obsBoot.ci95.map(v => +v.toFixed(6)) },
    noFeedback: { slope: +noFbReg.slope.toFixed(6), r: +noFbR.toFixed(4), n: noFbF.length, bootMean: +noFbBoot.mean.toFixed(6), bootSD: +noFbBoot.sd.toFixed(6), ci95: noFbBoot.ci95.map(v => +v.toFixed(6)) },
    withFeedback: { slope: +fbReg.slope.toFixed(6), r: +fbR.toFixed(4), n: fbF.length, bootMean: +fbBoot.mean.toFixed(6), bootSD: +fbBoot.sd.toFixed(6), ci95: fbBoot.ci95.map(v => +v.toFixed(6)) },
    excessNoFeedback: { delta: +deltaNoFb.toFixed(6), se: +seNoFb.toFixed(6), sigma: +sigmaNoFb.toFixed(2) },
    excessWithFeedback: { delta: +deltaFb.toFixed(6), se: +seFb.toFixed(6), sigma: +sigmaFb.toFixed(2) },
    feedbackExplainsPercent: +feedbackExplainsPercent.toFixed(1),
    feedbackFullyExplains: sigmaFb < 2.0
  });
}

results.feedbackModel = {
  description: "ΛCDM + baryonic feedback following Di Cintio et al. (2014) and Blumenthal et al. (1986)",
  components: [
    {
      name: "Supernova-driven core formation",
      reference: "Di Cintio et al. (2014), MNRAS 437, 415",
      description: "Inner DM profile slope α depends on M_star/M_halo ratio. Maximum core formation at M_star/M_halo ~ 10^-2.7. Transforms NFW cusps into shallow cores for dwarf galaxies.",
      effect: "Creates baryon-dependent DM profile → direct coupling between Σ_bar and inner halo structure"
    },
    {
      name: "Adiabatic contraction",
      reference: "Blumenthal et al. (1986), ApJ 301, 27",
      description: "Baryonic infall deepens DM potential in massive galaxies. Halo contracts when baryons condense at center.",
      effect: "Higher Σ_bar → more contraction → steeper DM profile in massive systems"
    },
    {
      name: "Abundance matching (SMHM relation)",
      reference: "Behroozi et al. (2013), ApJ 770, 57",
      description: "Star formation efficiency peaks at M_halo ~ 10^12 M☉, creating non-linear M_star-M_halo relation.",
      effect: "Creates mass-dependent baryon fractions that can produce apparent coupling"
    }
  ],
  nGalaxies: feedbackGals.length,
  nRealizations: N_REALIZATIONS,
  perRealization: N_PER_REAL
};

const maxSigmaNoFb = Math.max(...results.metrics.map(m => m.excessNoFeedback.sigma));
const maxSigmaFb = Math.max(...results.metrics.map(m => m.excessWithFeedback.sigma));
const avgReduction = results.metrics.reduce((s, m) => s + m.feedbackExplainsPercent, 0) / results.metrics.length;
const anyFullyExplained = results.metrics.some(m => m.feedbackFullyExplains);
const allFullyExplained = results.metrics.every(m => m.feedbackFullyExplains);

results.summary = {
  maxSigmaWithoutFeedback: +maxSigmaNoFb.toFixed(1),
  maxSigmaWithFeedback: +maxSigmaFb.toFixed(1),
  averageReductionPercent: +avgReduction.toFixed(1),
  anyMetricFullyExplained: anyFullyExplained,
  allMetricsFullyExplained: allFullyExplained,
  verdict: allFullyExplained
    ? "FEEDBACK FULLY EXPLAINS the coupling. The 'excess' was due to incomplete simulation. Claim retracted."
    : anyFullyExplained
    ? `FEEDBACK PARTIALLY EXPLAINS the coupling. Some metrics drop below 2σ but ${results.metrics.filter(m => !m.feedbackFullyExplains).map(m => m.name).join(', ')} remain significant.`
    : `FEEDBACK DOES NOT EXPLAIN the coupling. Even with core formation + adiabatic contraction, excess persists at ${maxSigmaFb.toFixed(1)}σ. Original claim stands but with reduced significance.`,
  honestAssessment: allFullyExplained
    ? "The baryon-halo coupling is an expected consequence of baryonic feedback, not new physics."
    : `Peak significance drops from ${maxSigmaNoFb.toFixed(1)}σ to ${maxSigmaFb.toFixed(1)}σ when feedback is included. Feedback explains ~${avgReduction.toFixed(0)}% of the tension. ${maxSigmaFb >= 3 ? 'The remaining excess is still significant (>3σ) and may indicate physics beyond standard feedback models.' : maxSigmaFb >= 2 ? 'The remaining excess is suggestive (2-3σ) but no longer at discovery level.' : 'The remaining excess is not statistically significant (<2σ). The claim should be downgraded.'}`,
  updatedClaim: maxSigmaFb >= 5
    ? `${maxSigmaFb.toFixed(1)}σ discovery-level excess PERSISTS even with feedback`
    : maxSigmaFb >= 3
    ? `${maxSigmaFb.toFixed(1)}σ significant excess persists — strong evidence but below discovery threshold`
    : maxSigmaFb >= 2
    ? `${maxSigmaFb.toFixed(1)}σ suggestive excess — interesting but needs more data`
    : `Excess drops below 2σ — original claim not supported when feedback is included`
};

console.log(`${'═'.repeat(70)}`);
console.log(`  FINAL VERDICT`);
console.log(`${'═'.repeat(70)}`);
console.log(`  ${results.summary.verdict}`);
console.log(`  ${results.summary.honestAssessment}`);
console.log(`  Updated claim: ${results.summary.updatedClaim}`);
console.log(`${'═'.repeat(70)}\n`);

const outPath = path.join(__dirname, '..', 'public', 'feedback-test.json');
fs.writeFileSync(outPath, JSON.stringify(results, null, 2));
console.log(`Results written to ${outPath}`);
