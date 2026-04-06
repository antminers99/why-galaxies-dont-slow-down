const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');
const dq415 = require('../public/phase415-dark-quarter.json');

const G = 4.3009e-6;

function normalize(n) {
  return n.replace(/\s+/g,'').replace(/^NGC0*/,'NGC').replace(/^DDO0*/,'DDO').replace(/^IC0*/,'IC').replace(/^UGC0*/,'UGC').replace(/^PGC0*/,'PGC').toUpperCase();
}

const tsContent = fs.readFileSync(path.join(__dirname, '..', 'src', 'data', 'sparc-datasets.ts'), 'utf8');
function parseRCData(ts) {
  const rcMap = {};
  const re = /"([^"]+)":\s*\{[^[]*data:\s*\[([\s\S]*?)\]/g;
  let m;
  while ((m = re.exec(ts)) !== null) {
    const points = [];
    const ptRe = /r:\s*([\d.]+)\s*,\s*v:\s*([\d.]+)/g;
    let pm;
    while ((pm = ptRe.exec(m[2])) !== null) points.push({ r: parseFloat(pm[1]), v: parseFloat(pm[2]) });
    if (points.length >= 3) rcMap[m[1]] = points;
  }
  return rcMap;
}
const rcMapRaw = parseRCData(tsContent);
const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) rcMap[normalize(name)] = pts;

function pearsonR(x, y) {
  const n = x.length; if (n < 4) return NaN;
  const mx = x.reduce((a, b) => a + b, 0) / n, my = y.reduce((a, b) => a + b, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) { const dx = x[i] - mx, dy = y[i] - my; num += dx * dy; dx2 += dx * dx; dy2 += dy * dy; }
  return dx2 > 0 && dy2 > 0 ? num / Math.sqrt(dx2 * dy2) : 0;
}
function multiR2(X, y) {
  const n = y.length, nv = X[0].length;
  const my = y.reduce((a, b) => a + b, 0) / n;
  const mx = Array(nv).fill(0);
  for (let j = 0; j < nv; j++) { for (let i = 0; i < n; i++) mx[j] += X[i][j]; mx[j] /= n; }
  const XTX = Array.from({ length: nv }, () => Array(nv).fill(0)), XTy = Array(nv).fill(0);
  for (let i = 0; i < n; i++) { for (let j = 0; j < nv; j++) { XTy[j] += (X[i][j] - mx[j]) * (y[i] - my); for (let k = 0; k < nv; k++) XTX[j][k] += (X[i][j] - mx[j]) * (X[i][k] - mx[k]); } }
  const aug = XTX.map((row, i) => [...row, XTy[i]]);
  for (let col = 0; col < nv; col++) { let maxRow = col; for (let row = col + 1; row < nv; row++) if (Math.abs(aug[row][col]) > Math.abs(aug[maxRow][col])) maxRow = row; [aug[col], aug[maxRow]] = [aug[maxRow], aug[col]]; if (Math.abs(aug[col][col]) < 1e-12) continue; for (let row = col + 1; row < nv; row++) { const f = aug[row][col] / aug[col][col]; for (let j = col; j <= nv; j++) aug[row][j] -= f * aug[col][j]; } }
  const beta = Array(nv).fill(0);
  for (let i = nv - 1; i >= 0; i--) { beta[i] = aug[i][nv]; for (let j = i + 1; j < nv; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= aug[i][i] || 1; }
  let sse = 0, sst = 0; const residuals = [];
  for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < nv; j++) pred += beta[j] * (X[i][j] - mx[j]); residuals.push(y[i] - pred); sse += (y[i] - pred) ** 2; sst += (y[i] - my) ** 2; }
  return { R2: sst > 0 ? 1 - sse / sst : 0, beta, residuals };
}

function gaussRand() {
  let u = 0, v = 0;
  while (u === 0) u = Math.random();
  while (v === 0) v = Math.random();
  return Math.sqrt(-2.0 * Math.log(u)) * Math.cos(2.0 * Math.PI * v);
}

const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[g.name] = g; resultsMap[normalize(g.name)] = g; });

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sr) continue;
  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 5) continue;

  const Vflat = sp.Vflat;
  const Rdisk = sp.Rdisk;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  const Rmax = rc[rc.length - 1].r;

  const k_halo = sr.models.dark_halo_linear.k;
  const mse_newt = sr.models.newtonian.mse;
  const mse_halo = sr.models.dark_halo_linear.mse;

  const V_Newt_Rmax = Math.sqrt(G * Mbar / Rmax);
  const dmFrac_Rmax = Math.max(0, 1 - (V_Newt_Rmax / Vflat) ** 2);
  const logK_halo = k_halo > 0 ? Math.log10(k_halo) : -5;
  const haloResponse = mse_newt > 0 ? Math.log10(Math.max(mse_newt / Math.max(mse_halo, 0.001), 0.01)) : 0;
  const logMbar = Math.log10(Math.max(Mbar, 1));
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(Rdisk, 0.01));
  const logSBdisk = Math.log10(Math.max(sp.SBdisk, 0.01));

  const outerPts = rc.filter(p => p.r > 3 * Rdisk);
  let outerSlope = 0;
  if (outerPts.length >= 3) {
    const lx = outerPts.map(p => Math.log10(p.r)), ly = outerPts.map(p => Math.log10(p.v));
    const mx2 = lx.reduce((a, b) => a + b, 0) / lx.length, my2 = ly.reduce((a, b) => a + b, 0) / ly.length;
    let num = 0, den = 0;
    for (let i = 0; i < lx.length; i++) { num += (lx[i] - mx2) * (ly[i] - my2); den += (lx[i] - mx2) ** 2; }
    outerSlope = den > 0 ? num / den : 0;
  }

  gals.push({
    name: g.name, logA0: g.logA0,
    Vflat, Rdisk, Rmax, Mbar, rc,
    logVflat: Math.log10(Vflat),
    logL36, logRdisk, logMbar,
    logMHI: g.logMHI,
    morphT: sp.T,
    logSBdisk,
    envCode: g.envCode,
    logK_halo, dmFrac_Rmax,
    haloResponse, outerSlope,
  });
}

const N = gals.length;

const struct4 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
const struct6 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
const vfModel = multiR2(struct4, gals.map(g => g.logVflat));
const a0Model = multiR2(struct6, gals.map(g => g.logA0));
for (let i = 0; i < N; i++) {
  gals[i].VfResid = vfModel.residuals[i];
  gals[i].a0Resid = a0Model.residuals[i];
}
const sdVf = Math.sqrt(gals.reduce((a, g) => a + g.VfResid ** 2, 0) / N);
const sdA0 = Math.sqrt(gals.reduce((a, g) => a + g.a0Resid ** 2, 0) / N);
for (let i = 0; i < N; i++) {
  gals[i].VfResid_z = gals[i].VfResid / sdVf;
  gals[i].a0Resid_z = gals[i].a0Resid / sdA0;
  gals[i].L_sum = gals[i].VfResid_z + gals[i].a0Resid_z;
}

const bestControls = gals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.envCode]);
const Lsum_from_best = multiR2(bestControls, gals.map(g => g.L_sum));
const trueDQ = Lsum_from_best.residuals;

const trueFingerprint = {
  rVf: pearsonR(trueDQ, gals.map(g => g.VfResid)),
  rA0: pearsonR(trueDQ, gals.map(g => g.a0Resid)),
  rHaloResp: pearsonR(trueDQ, gals.map(g => g.haloResponse)),
  rOuterSlope: pearsonR(trueDQ, gals.map(g => g.outerSlope)),
  rVflat: pearsonR(trueDQ, gals.map(g => g.Vflat)),
};

console.log('='.repeat(70));
console.log('PHASE 417: RESOLVE THE DARK QUARTER');
console.log('='.repeat(70));
console.log('N = ' + N);
console.log('True DQ fingerprint:');
console.log('  r(DQ, VfResid) = ' + trueFingerprint.rVf.toFixed(3));
console.log('  r(DQ, a0_resid) = ' + trueFingerprint.rA0.toFixed(3));
console.log('  r(DQ, haloResponse) = ' + trueFingerprint.rHaloResp.toFixed(3));
console.log('  r(DQ, outerSlope) = ' + trueFingerprint.rOuterSlope.toFixed(3));
console.log('  r(DQ, Vflat) = ' + trueFingerprint.rVflat.toFixed(3));


console.log('\n\n' + '#'.repeat(70));
console.log('417B: 3D HALO SIMULATION TRACK');
console.log('Can triaxiality/disequilibrium REPRODUCE the DQ fingerprint?');
console.log('#'.repeat(70));


console.log('\n\n' + '='.repeat(50));
console.log('417B.1: NFW TRIAXIAL HALO MODEL');
console.log('Full 3D NFW halo with axis ratios and random orientation');
console.log('='.repeat(50));

function nfwVcirc(r, rs, rho0) {
  const x = r / rs;
  const mass = 4 * Math.PI * rho0 * rs * rs * rs * (Math.log(1 + x) - x / (1 + x));
  return Math.sqrt(G * mass / r);
}

function triaxialNFW(r_obs, rs, rho0, axisRatios, viewAngles) {
  const { b_a, c_a } = axisRatios;
  const { theta, phi_v } = viewAngles;

  const ct = Math.cos(theta), st = Math.sin(theta);
  const cp = Math.cos(phi_v), sp = Math.sin(phi_v);

  const q_eff_sq = (ct * ct) * (cp * cp + sp * sp / (b_a * b_a)) +
    (st * st) / (c_a * c_a);
  const q_eff = Math.sqrt(q_eff_sq);

  const r_3d = r_obs * q_eff;
  return nfwVcirc(r_3d, rs, rho0);
}

function generateTriaxialGalaxy(g, triaxModel) {
  const Mbar = g.Mbar;
  const rc = g.rc;
  const Rdisk = g.Rdisk;

  const c = 5 + 15 * Math.random();
  const Mvir = Mbar * (10 + 40 * Math.random());
  const Rvir = Math.pow(3 * Mvir / (4 * Math.PI * 200 * 0.127), 1 / 3);
  const rs = Rvir / c;
  const rho0 = Mvir / (4 * Math.PI * rs * rs * rs * (Math.log(1 + c) - c / (1 + c)));

  let b_a, c_a;
  if (triaxModel.correlateWithConc) {
    const cNorm = (c - 10) / 5;
    b_a = Math.max(0.4, Math.min(1.0, triaxModel.b_a_mean + triaxModel.b_a_scatter * gaussRand() - triaxModel.concCorr * cNorm * 0.1));
    c_a = Math.max(0.3, Math.min(b_a, triaxModel.c_a_mean + triaxModel.c_a_scatter * gaussRand() - triaxModel.concCorr * cNorm * 0.15));
  } else {
    b_a = Math.max(0.4, Math.min(1.0, triaxModel.b_a_mean + triaxModel.b_a_scatter * gaussRand()));
    c_a = Math.max(0.3, Math.min(b_a, triaxModel.c_a_mean + triaxModel.c_a_scatter * gaussRand()));
  }

  const theta = Math.acos(2 * Math.random() - 1);
  const phi_v = 2 * Math.PI * Math.random();

  const newRc = rc.map(p => {
    const V_bar = Math.sqrt(G * Mbar / p.r) * Math.min(1, Math.sqrt(p.r / (2.2 * Rdisk)));
    const V_halo = triaxialNFW(p.r, rs, rho0, { b_a, c_a }, { theta, phi_v });
    const V_tot = Math.sqrt(V_bar * V_bar + V_halo * V_halo);
    return { r: p.r, v: Math.max(V_tot, 5) };
  });

  const Vflat = newRc.reduce((mx, p) => Math.max(mx, p.v), 0);
  const V_Newt_Rmax = Math.sqrt(G * Mbar / g.Rmax);
  const dmFrac = Math.max(0, 1 - (V_Newt_Rmax / Vflat) ** 2);

  const logK_halo_new = Math.log10(Math.max(Vflat * Vflat / (G * Mbar / (2.2 * Rdisk)), 0.01));

  let outerSl = 0;
  const outerPts = newRc.filter(p => p.r > 3 * Rdisk);
  if (outerPts.length >= 3) {
    const lx = outerPts.map(p => Math.log10(p.r)), ly = outerPts.map(p => Math.log10(p.v));
    const mx2 = lx.reduce((a, b) => a + b, 0) / lx.length, my2 = ly.reduce((a, b) => a + b, 0) / ly.length;
    let n2 = 0, d2 = 0;
    for (let i = 0; i < lx.length; i++) { n2 += (lx[i] - mx2) * (ly[i] - my2); d2 += (lx[i] - mx2) ** 2; }
    outerSl = d2 > 0 ? n2 / d2 : 0;
  }

  const mse_newt_new = newRc.reduce((a, p) => a + (p.v - Math.sqrt(G * Mbar / p.r)) ** 2, 0) / newRc.length;
  const mse_halo_new = newRc.reduce((a, p) => {
    const vModel = Math.sqrt(G * Mbar / p.r + Vflat * Vflat * 0.5);
    return a + (p.v - vModel) ** 2;
  }, 0) / newRc.length;
  const haloResp = mse_newt_new > 0 ? Math.log10(Math.max(mse_newt_new / Math.max(mse_halo_new, 0.01), 0.01)) : 0;

  return {
    ...g,
    Vflat,
    logVflat: Math.log10(Math.max(Vflat, 1)),
    rc: newRc,
    dmFrac_Rmax: dmFrac,
    logK_halo: logK_halo_new,
    haloResponse: haloResp,
    outerSlope: outerSl,
    triax: { b_a, c_a, theta, phi_v, c, rs },
  };
}

function runFullPipeline(mockGals) {
  const n = mockGals.length;
  const s4 = mockGals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
  const s6 = mockGals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
  const vf = multiR2(s4, mockGals.map(g => g.logVflat));
  const a0 = multiR2(s6, mockGals.map(g => g.logA0));

  const sdV = Math.sqrt(vf.residuals.reduce((a, v) => a + v * v, 0) / n) || 1;
  const sdA = Math.sqrt(a0.residuals.reduce((a, v) => a + v * v, 0) / n) || 1;
  const Lsum = vf.residuals.map((v, i) => v / sdV + a0.residuals[i] / sdA);

  const ctrl = mockGals.map(g => [g.logK_halo, g.dmFrac_Rmax, g.envCode]);
  const Lmod = multiR2(ctrl, Lsum);
  const dq = Lmod.residuals;

  const rChannel = pearsonR(vf.residuals, a0.residuals);
  const rDQ_VfR = pearsonR(dq, vf.residuals);
  const rDQ_a0R = pearsonR(dq, a0.residuals);
  const rDQ_haloResp = pearsonR(dq, mockGals.map(g => g.haloResponse));
  const rDQ_outerSlope = pearsonR(dq, mockGals.map(g => g.outerSlope));
  const rDQ_Vflat = pearsonR(dq, mockGals.map(g => g.Vflat));

  return {
    rChannel, rDQ_VfR, rDQ_a0R, rDQ_haloResp, rDQ_outerSlope, rDQ_Vflat,
    bilateral: rDQ_VfR > 0.3 && rDQ_a0R > 0.3,
    dq,
  };
}

function fingerprintMatch(result) {
  const targets = [
    { key: 'rDQ_VfR', target: trueFingerprint.rVf, weight: 2 },
    { key: 'rDQ_a0R', target: trueFingerprint.rA0, weight: 2 },
    { key: 'rDQ_haloResp', target: trueFingerprint.rHaloResp, weight: 1 },
    { key: 'rDQ_outerSlope', target: trueFingerprint.rOuterSlope, weight: 1 },
    { key: 'rDQ_Vflat', target: trueFingerprint.rVflat, weight: 1 },
  ];
  let score = 0, total = 0;
  for (const t of targets) {
    const diff = Math.abs(result[t.key] - t.target);
    score += t.weight * Math.max(0, 1 - diff / 0.3);
    total += t.weight;
  }
  return score / total;
}

const triaxModels = [
  { name: 'Spherical control', b_a_mean: 1.0, c_a_mean: 1.0, b_a_scatter: 0.0, c_a_scatter: 0.0, correlateWithConc: false, concCorr: 0 },
  { name: 'Mildly triaxial', b_a_mean: 0.85, c_a_mean: 0.75, b_a_scatter: 0.08, c_a_scatter: 0.10, correlateWithConc: false, concCorr: 0 },
  { name: 'Moderately triaxial', b_a_mean: 0.75, c_a_mean: 0.60, b_a_scatter: 0.10, c_a_scatter: 0.12, correlateWithConc: false, concCorr: 0 },
  { name: 'Strongly triaxial', b_a_mean: 0.65, c_a_mean: 0.50, b_a_scatter: 0.12, c_a_scatter: 0.15, correlateWithConc: false, concCorr: 0 },
  { name: 'Conc-correlated mild', b_a_mean: 0.85, c_a_mean: 0.75, b_a_scatter: 0.08, c_a_scatter: 0.10, correlateWithConc: true, concCorr: 0.3 },
  { name: 'Conc-correlated moderate', b_a_mean: 0.75, c_a_mean: 0.60, b_a_scatter: 0.10, c_a_scatter: 0.12, correlateWithConc: true, concCorr: 0.5 },
  { name: 'Conc-correlated strong', b_a_mean: 0.65, c_a_mean: 0.50, b_a_scatter: 0.12, c_a_scatter: 0.15, correlateWithConc: true, concCorr: 0.7 },
];

const nTrials = 200;

console.log('\n--- 417B.1: Triaxial NFW halo models ---');
console.log('  Each model: ' + nTrials + ' realizations of N=' + N + ' mock galaxies');
console.log('  True fingerprint: rVf=' + trueFingerprint.rVf.toFixed(3) + ' rA0=' + trueFingerprint.rA0.toFixed(3) + ' rHR=' + trueFingerprint.rHaloResp.toFixed(3) + ' rOS=' + trueFingerprint.rOuterSlope.toFixed(3) + ' rVflat=' + trueFingerprint.rVflat.toFixed(3));
console.log('');

const modelResults = {};

for (const model of triaxModels) {
  const trials = [];
  for (let t = 0; t < nTrials; t++) {
    const mockGals = gals.map(g => generateTriaxialGalaxy(g, model));
    const result = runFullPipeline(mockGals);
    trials.push({ ...result, fmatch: fingerprintMatch(result) });
  }

  const med = (arr) => { const s = [...arr].sort((a, b) => a - b); return s[Math.floor(s.length / 2)]; };
  const p95 = (arr) => { const s = [...arr].sort((a, b) => a - b); return s[Math.floor(s.length * 0.95)]; };

  const summary = {
    rChannel_med: med(trials.map(t => t.rChannel)),
    rDQ_VfR_med: med(trials.map(t => t.rDQ_VfR)),
    rDQ_a0R_med: med(trials.map(t => t.rDQ_a0R)),
    rDQ_haloResp_med: med(trials.map(t => t.rDQ_haloResp)),
    rDQ_outerSlope_med: med(trials.map(t => t.rDQ_outerSlope)),
    rDQ_Vflat_med: med(trials.map(t => t.rDQ_Vflat)),
    bilateralRate: trials.filter(t => t.bilateral).length / nTrials,
    fmatch_med: med(trials.map(t => t.fmatch)),
    fmatch_p95: p95(trials.map(t => t.fmatch)),
    fmatch_max: Math.max(...trials.map(t => t.fmatch)),
  };

  modelResults[model.name] = summary;

  console.log('  ' + model.name);
  console.log('    Channel r: ' + summary.rChannel_med.toFixed(3));
  console.log('    DQ: rVfR=' + summary.rDQ_VfR_med.toFixed(3) + ' rA0R=' + summary.rDQ_a0R_med.toFixed(3) + ' bilateral=' + (summary.bilateralRate * 100).toFixed(0) + '%');
  console.log('    Fingerprint: rHR=' + summary.rDQ_haloResp_med.toFixed(3) + ' rOS=' + summary.rDQ_outerSlope_med.toFixed(3) + ' rVflat=' + summary.rDQ_Vflat_med.toFixed(3));
  console.log('    Match score: med=' + summary.fmatch_med.toFixed(3) + ' 95th=' + summary.fmatch_p95.toFixed(3) + ' max=' + summary.fmatch_max.toFixed(3));
  console.log('');
}


console.log('\n' + '='.repeat(50));
console.log('417B.2: DISEQUILIBRIUM + TRIAXIALITY COMBINED');
console.log('Physically motivated: merging/accreting halos');
console.log('='.repeat(50));

function generateDisequilGalaxy(g, diseqModel) {
  const mockG = generateTriaxialGalaxy(g, diseqModel.triaxParams);

  const amplitude = diseqModel.amplitude;
  const massSlope = diseqModel.massSlope;
  const outerBoost = diseqModel.outerBoost;
  const responseCorr = diseqModel.responseCorr;

  const massFactor = massSlope * (g.logMbar - 10);
  const perturbation = amplitude * gaussRand() + massFactor;

  const newRc = mockG.rc.map(p => {
    const rNorm = p.r / g.Rdisk;
    const radialFactor = rNorm > 3 ? (1 + outerBoost * (rNorm - 3) / 10) : 1;
    const v = p.v * (1 + perturbation * radialFactor);
    return { r: p.r, v: Math.max(v, 5) };
  });

  const Vflat = newRc.reduce((mx, p) => Math.max(mx, p.v), 0);
  const V_Newt = Math.sqrt(G * g.Mbar / g.Rmax);
  const dmFrac = Math.max(0, 1 - (V_Newt / Vflat) ** 2);

  let outerSl = 0;
  const outerPts = newRc.filter(p => p.r > 3 * g.Rdisk);
  if (outerPts.length >= 3) {
    const lx = outerPts.map(p => Math.log10(p.r)), ly = outerPts.map(p => Math.log10(p.v));
    const mx2 = lx.reduce((a, b) => a + b, 0) / lx.length, my2 = ly.reduce((a, b) => a + b, 0) / ly.length;
    let n2 = 0, d2 = 0;
    for (let i = 0; i < lx.length; i++) { n2 += (lx[i] - mx2) * (ly[i] - my2); d2 += (lx[i] - mx2) ** 2; }
    outerSl = d2 > 0 ? n2 / d2 : 0;
  }

  const logK_new = Math.log10(Math.max(Vflat * Vflat / (G * g.Mbar / (2.2 * g.Rdisk)), 0.01));
  const mse_newt_new = newRc.reduce((a, p) => a + (p.v - Math.sqrt(G * g.Mbar / p.r)) ** 2, 0) / newRc.length;
  const mse_halo_new = newRc.reduce((a, p) => {
    const vModel = Math.sqrt(G * g.Mbar / p.r + Vflat * Vflat * 0.3);
    return a + (p.v - vModel) ** 2;
  }, 0) / newRc.length;
  const haloResp = Math.log10(Math.max(mse_newt_new / Math.max(mse_halo_new, 0.01), 0.01));

  return {
    ...mockG,
    Vflat, logVflat: Math.log10(Math.max(Vflat, 1)),
    rc: newRc, dmFrac_Rmax: dmFrac,
    logK_halo: logK_new, haloResponse: haloResp, outerSlope: outerSl,
  };
}

const diseqModels = [
  {
    name: 'Mild triax + mild diseq',
    triaxParams: { b_a_mean: 0.85, c_a_mean: 0.75, b_a_scatter: 0.08, c_a_scatter: 0.10, correlateWithConc: false, concCorr: 0 },
    amplitude: 0.03, massSlope: 0, outerBoost: 0.1, responseCorr: 0,
  },
  {
    name: 'Mod triax + outer-boosted diseq',
    triaxParams: { b_a_mean: 0.75, c_a_mean: 0.60, b_a_scatter: 0.10, c_a_scatter: 0.12, correlateWithConc: true, concCorr: 0.3 },
    amplitude: 0.05, massSlope: 0.01, outerBoost: 0.3, responseCorr: 0.2,
  },
  {
    name: 'Strong triax + mass-dep diseq',
    triaxParams: { b_a_mean: 0.65, c_a_mean: 0.50, b_a_scatter: 0.12, c_a_scatter: 0.15, correlateWithConc: true, concCorr: 0.5 },
    amplitude: 0.08, massSlope: 0.02, outerBoost: 0.5, responseCorr: 0.3,
  },
  {
    name: 'DQ-optimized (tuned to fingerprint)',
    triaxParams: { b_a_mean: 0.70, c_a_mean: 0.55, b_a_scatter: 0.15, c_a_scatter: 0.15, correlateWithConc: true, concCorr: 0.6 },
    amplitude: 0.06, massSlope: 0.015, outerBoost: 0.4, responseCorr: 0.4,
  },
];

const diseqResults = {};

for (const model of diseqModels) {
  const trials = [];
  for (let t = 0; t < nTrials; t++) {
    const mockGals = gals.map(g => generateDisequilGalaxy(g, model));
    const result = runFullPipeline(mockGals);
    trials.push({ ...result, fmatch: fingerprintMatch(result) });
  }

  const med = (arr) => { const s = [...arr].sort((a, b) => a - b); return s[Math.floor(s.length / 2)]; };
  const p95 = (arr) => { const s = [...arr].sort((a, b) => a - b); return s[Math.floor(s.length * 0.95)]; };

  const summary = {
    rChannel_med: med(trials.map(t => t.rChannel)),
    rDQ_VfR_med: med(trials.map(t => t.rDQ_VfR)),
    rDQ_a0R_med: med(trials.map(t => t.rDQ_a0R)),
    rDQ_haloResp_med: med(trials.map(t => t.rDQ_haloResp)),
    rDQ_outerSlope_med: med(trials.map(t => t.rDQ_outerSlope)),
    rDQ_Vflat_med: med(trials.map(t => t.rDQ_Vflat)),
    bilateralRate: trials.filter(t => t.bilateral).length / nTrials,
    fmatch_med: med(trials.map(t => t.fmatch)),
    fmatch_p95: p95(trials.map(t => t.fmatch)),
    fmatch_max: Math.max(...trials.map(t => t.fmatch)),
  };

  diseqResults[model.name] = summary;

  console.log('\n  ' + model.name);
  console.log('    Channel r: ' + summary.rChannel_med.toFixed(3));
  console.log('    DQ: rVfR=' + summary.rDQ_VfR_med.toFixed(3) + ' rA0R=' + summary.rDQ_a0R_med.toFixed(3) + ' bilateral=' + (summary.bilateralRate * 100).toFixed(0) + '%');
  console.log('    Fingerprint: rHR=' + summary.rDQ_haloResp_med.toFixed(3) + ' rOS=' + summary.rDQ_outerSlope_med.toFixed(3) + ' rVflat=' + summary.rDQ_Vflat_med.toFixed(3));
  console.log('    Match score: med=' + summary.fmatch_med.toFixed(3) + ' 95th=' + summary.fmatch_p95.toFixed(3) + ' max=' + summary.fmatch_max.toFixed(3));
}


console.log('\n\n' + '='.repeat(50));
console.log('417A: IFU / 2D OBSERVATIONAL TRACK');
console.log('Target selection for follow-up observations');
console.log('='.repeat(50));

const dqRanked = gals.map((g, i) => ({
  name: g.name, dq: trueDQ[i], Vflat: g.Vflat, logMbar: g.logMbar,
  morphT: g.morphT, env: g.envCode, haloResp: g.haloResponse,
})).sort((a, b) => b.dq - a.dq);

console.log('\n--- 417A.1: IFU target selection ---');
console.log('  Criterion: high |DQ|, matched in mass/Vflat for control');

const targets = dqRanked.slice(0, 5);
const controls = [];
for (const tgt of targets) {
  const match = dqRanked.filter(g =>
    Math.abs(g.dq) < 0.5 &&
    Math.abs(g.logMbar - tgt.logMbar) < 0.3 &&
    Math.abs(g.Vflat - tgt.Vflat) < 40 &&
    g.name !== tgt.name
  ).sort((a, b) => Math.abs(a.dq) - Math.abs(b.dq));
  if (match.length > 0) controls.push({ target: tgt.name, control: match[0].name, targetDQ: tgt.dq, controlDQ: match[0].dq });
}

console.log('\n  Target          DQ     Vflat  logMbar  ->  Control           DQ');
console.log('  ' + '-'.repeat(75));
for (const pair of controls) {
  const tgt = dqRanked.find(g => g.name === pair.target);
  const ctrl = dqRanked.find(g => g.name === pair.control);
  console.log('  ' + pair.target.padEnd(18) + (pair.targetDQ >= 0 ? '+' : '') + pair.targetDQ.toFixed(2).padEnd(7) + ' ' + (tgt ? tgt.Vflat.toFixed(0) : '?').padEnd(7) + ' ' + (tgt ? tgt.logMbar.toFixed(2) : '?').padEnd(9) + '->  ' + pair.control.padEnd(18) + (pair.controlDQ >= 0 ? '+' : '') + pair.controlDQ.toFixed(2));
}

console.log('\n--- 417A.2: What to look for in IFU data ---');
console.log('  High DQ galaxies should show (if DQ = 3D halo projection):');
console.log('    1. Non-circular motions: kinematic position angle twists');
console.log('    2. Asymmetric velocity fields: lopsided rotation');
console.log('    3. Radial flows: oval/bar-driven inflows');
console.log('    4. Outer kinematic warps: misaligned outer disk');
console.log('    5. Velocity dispersion anomalies: unrelaxed components');
console.log('  Low DQ controls should show NONE of these at the same level.');

const ifuSurveys = [
  { name: 'THINGS HI', desc: 'High-res HI velocity fields for NGC2841, NGC5005, NGC3521, NGC2903, NGC5055' },
  { name: 'HERACLES CO', desc: 'CO velocity fields for NGC2841, NGC5005, NGC3521' },
  { name: 'DiskMass Survey', desc: 'Stellar velocity dispersions for edge-on galaxies' },
  { name: 'SPARC extended', desc: 'Some galaxies have multi-slit data beyond published RC' },
];

console.log('\n--- 417A.3: Available 2D/IFU surveys ---');
for (const s of ifuSurveys) {
  console.log('  ' + s.name.padEnd(20) + ' ' + s.desc);
}


console.log('\n\n' + '='.repeat(70));
console.log('PHASE 417 GRAND VERDICT');
console.log('='.repeat(70));

let bestModel = '';
let bestScore = 0;
for (const [name, res] of [...Object.entries(modelResults), ...Object.entries(diseqResults)]) {
  if (res.fmatch_med > bestScore) { bestScore = res.fmatch_med; bestModel = name; }
}

console.log('\n  417B — 3D HALO SIMULATION:');
console.log('    Best model: ' + bestModel + ' (match score = ' + bestScore.toFixed(3) + ')');
console.log('    Target fingerprint: rVfR=' + trueFingerprint.rVf.toFixed(3) + ' rA0R=' + trueFingerprint.rA0.toFixed(3) + ' rHR=' + trueFingerprint.rHaloResp.toFixed(3));

const bestRes = modelResults[bestModel] || diseqResults[bestModel];
if (bestRes) {
  console.log('    Best achieved: rVfR=' + bestRes.rDQ_VfR_med.toFixed(3) + ' rA0R=' + bestRes.rDQ_a0R_med.toFixed(3) + ' rHR=' + bestRes.rDQ_haloResp_med.toFixed(3));
}

if (bestScore > 0.7) {
  console.log('    -> STRONG MATCH: 3D halo geometry CAN reproduce the DQ fingerprint.');
  console.log('       The Dark Quarter is likely the projection of triaxial/non-equilibrium halos.');
} else if (bestScore > 0.4) {
  console.log('    -> PARTIAL MATCH: 3D halo geometry reproduces SOME features of the DQ.');
  console.log('       Additional physics may be needed for the complete picture.');
} else {
  console.log('    -> POOR MATCH: 3D halo geometry alone CANNOT reproduce the DQ fingerprint.');
  console.log('       The Dark Quarter may require fundamentally different physics.');
}

console.log('\n  417A — IFU FOLLOW-UP:');
console.log('    ' + controls.length + ' target-control pairs identified for IFU observation.');
console.log('    Key prediction: high-DQ galaxies show non-circular motions, low-DQ controls do not.');

console.log('\n  OVERALL:');
console.log('    The VfResid-a0 channel decomposes into:');
console.log('    1. 75% = concentration-mass scatter (IDENTIFIED, Phase 414)');
console.log('    2. 25% = Dark Quarter (CONFIRMED REAL, Phase 415-416)');
console.log('    3. DQ fingerprint: fast rotators, strong halo response, rising outer RCs');
console.log('    4. 3D halo model match: ' + bestScore.toFixed(3));
console.log('    5. Resolution requires IFU/2D kinematic data or cosmological simulations.');


const outPath = path.join(__dirname, '..', 'public', 'phase417-resolve-dq.json');
fs.writeFileSync(outPath, JSON.stringify({
  phase: '417',
  title: 'Resolve the Dark Quarter',
  timestamp: new Date().toISOString(),
  N,
  trueFingerprint,
  triaxModels: modelResults,
  diseqModels: diseqResults,
  bestModel, bestScore,
  ifuTargets: controls,
}, null, 2));
console.log('\nSaved: ' + outPath);
