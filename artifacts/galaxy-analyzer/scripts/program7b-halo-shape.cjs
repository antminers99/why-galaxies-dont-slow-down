const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparc = require('../public/sparc-table.json');
const sparcResults = require('../public/sparc-results.json');

const G = 4.3009e-6;

function normalize(n) {
  return n.replace(/\s+/g, '').replace(/^NGC0*/, 'NGC').replace(/^DDO0*/, 'DDO').replace(/^IC0*/, 'IC').replace(/^UGC0*/, 'UGC').replace(/^PGC0*/, 'PGC').toUpperCase();
}

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
  const resid = [];
  for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < nv; j++) pred += beta[j] * (X[i][j] - mx[j]); resid.push(y[i] - pred); }
  return { residuals: resid, beta };
}

const sparcMap = {};
sparc.forEach(g => { sparcMap[g.name] = g; sparcMap[normalize(g.name)] = g; });
const resultsMap = {};
sparcResults.perGalaxy.forEach(g => { resultsMap[g.name] = g; resultsMap[normalize(g.name)] = g; });

const tsContent = fs.readFileSync(path.join(__dirname, '..', 'src', 'data', 'sparc-datasets.ts'), 'utf8');
const rcMapRaw = {};
const re = /"([^"]+)":\s*\{[^[]*data:\s*\[([\s\S]*?)\]/g;
let m;
while ((m = re.exec(tsContent)) !== null) {
  const pts = [];
  const ptRe = /r:\s*([\d.]+)\s*,\s*v:\s*([\d.]+)/g;
  let pm;
  while ((pm = ptRe.exec(m[2])) !== null) pts.push({ r: parseFloat(pm[1]), v: parseFloat(pm[2]) });
  if (pts.length >= 3) rcMapRaw[m[1]] = pts;
}
const rcMap = {};
for (const [name, pts] of Object.entries(rcMapRaw)) rcMap[normalize(name)] = pts;

function fitHaloProfile(rc, Mbar, Rdisk, Rmax) {
  const haloVelocities = [];
  for (const p of rc) {
    if (p.r < 0.3) continue;
    const encMass = Mbar * Math.min(p.r / Rmax, 1);
    const V_bar = Math.sqrt(G * encMass / p.r);
    const V_halo_sq = Math.max(p.v * p.v - V_bar * V_bar, 0);
    const V_halo = Math.sqrt(V_halo_sq);
    haloVelocities.push({ r: p.r, V_halo, V_obs: p.v, V_bar });
  }
  if (haloVelocities.length < 5) return null;

  function nfwVelocity(r, V200, c, Rs) {
    const x = r / Rs;
    const gx = Math.log(1 + x) - x / (1 + x);
    const gc = Math.log(1 + c) - c / (1 + c);
    return V200 * Math.sqrt(gx / (x * gc));
  }

  function burkertVelocity(r, rho0, r0) {
    const x = r / r0;
    const term1 = Math.log(1 + x);
    const term2 = 0.5 * Math.log(1 + x * x);
    const term3 = Math.atan(x);
    const Mx = rho0 * r0 * r0 * r0 * Math.PI * (term1 + term2 - term3);
    return Math.sqrt(G * Math.abs(Mx) / r);
  }

  function pIsoVelocity(r, rho0, rc) {
    const Vsq = 4 * Math.PI * G * rho0 * rc * rc * (1 - (rc / r) * Math.atan(r / rc));
    return Math.sqrt(Math.max(Vsq, 0));
  }

  function fitMSE(predictFn, haloV) {
    let sse = 0;
    for (const h of haloV) {
      const pred = predictFn(h.r);
      sse += (h.V_halo - pred) ** 2;
    }
    return sse / haloV.length;
  }

  let bestNFW = { mse: Infinity, c: 0, V200: 0 };
  for (let c = 2; c <= 40; c += 2) {
    for (let V200 = 20; V200 <= 400; V200 += 20) {
      const Rs = Rmax / c;
      const mse = fitMSE(r => nfwVelocity(r, V200, c, Rs), haloVelocities);
      if (mse < bestNFW.mse) bestNFW = { mse, c, V200, Rs };
    }
  }
  for (let dc = -3; dc <= 3; dc += 0.5) {
    for (let dV = -15; dV <= 15; dV += 5) {
      const c = bestNFW.c + dc;
      const V200 = bestNFW.V200 + dV;
      if (c < 1 || V200 < 5) continue;
      const Rs = Rmax / c;
      const mse = fitMSE(r => nfwVelocity(r, V200, c, Rs), haloVelocities);
      if (mse < bestNFW.mse) bestNFW = { mse, c, V200, Rs };
    }
  }

  let bestBurkert = { mse: Infinity, rho0: 0, r0: 0 };
  for (let r0 = 0.5; r0 <= 20; r0 += 1) {
    for (let lrho = -4; lrho <= 0; lrho += 0.5) {
      const rho0 = Math.pow(10, lrho);
      const mse = fitMSE(r => burkertVelocity(r, rho0, r0), haloVelocities);
      if (mse < bestBurkert.mse) bestBurkert = { mse, rho0, r0 };
    }
  }
  for (let dr = -0.8; dr <= 0.8; dr += 0.2) {
    for (let dlrho = -0.4; dlrho <= 0.4; dlrho += 0.1) {
      const r0 = bestBurkert.r0 + dr;
      const rho0 = bestBurkert.rho0 * Math.pow(10, dlrho);
      if (r0 < 0.1 || rho0 < 1e-6) continue;
      const mse = fitMSE(r => burkertVelocity(r, rho0, r0), haloVelocities);
      if (mse < bestBurkert.mse) bestBurkert = { mse, rho0, r0 };
    }
  }

  let bestPiso = { mse: Infinity, rho0: 0, rc: 0 };
  for (let rc = 0.5; rc <= 20; rc += 1) {
    for (let lrho = -4; lrho <= 0; lrho += 0.5) {
      const rho0 = Math.pow(10, lrho);
      const mse = fitMSE(r => pIsoVelocity(r, rho0, rc), haloVelocities);
      if (mse < bestPiso.mse) bestPiso = { mse, rho0, rc };
    }
  }
  for (let dr = -0.8; dr <= 0.8; dr += 0.2) {
    for (let dlrho = -0.4; dlrho <= 0.4; dlrho += 0.1) {
      const rc = bestPiso.rc + dr;
      const rho0 = bestPiso.rho0 * Math.pow(10, dlrho);
      if (rc < 0.1 || rho0 < 1e-6) continue;
      const mse = fitMSE(r => pIsoVelocity(r, rho0, rc), haloVelocities);
      if (mse < bestPiso.mse) bestPiso = { mse, rho0, rc };
    }
  }

  const innerSlope = computeInnerSlope(haloVelocities);

  const radialProfile = computeRadialProfile(haloVelocities, Rdisk, Rmax);

  const profiles = [
    { name: 'NFW', mse: bestNFW.mse, params: { c: bestNFW.c, V200: bestNFW.V200 } },
    { name: 'Burkert', mse: bestBurkert.mse, params: { rho0: bestBurkert.rho0, r0: bestBurkert.r0 } },
    { name: 'pIso', mse: bestPiso.mse, params: { rho0: bestPiso.rho0, rc: bestPiso.rc } },
  ];
  profiles.sort((a, b) => a.mse - b.mse);
  const bestFamily = profiles[0].name;
  const familyAdvantage = profiles.length > 1 ? Math.log10(Math.max(profiles[1].mse / profiles[0].mse, 0.01)) : 0;

  return {
    nfw: bestNFW,
    burkert: bestBurkert,
    pIso: bestPiso,
    bestFamily,
    familyAdvantage,
    innerSlope,
    radialProfile,
    concentration: bestNFW.c,
    nfwMSE: bestNFW.mse,
    burkertMSE: bestBurkert.mse,
    pisoMSE: bestPiso.mse,
  };
}

function computeInnerSlope(haloV) {
  const inner = haloV.filter(h => h.r > 0.3 && h.r < 3);
  if (inner.length < 3) return NaN;
  const logR = inner.map(h => Math.log10(h.r));
  const logVh = inner.map(h => Math.log10(Math.max(h.V_halo, 0.1)));
  const mr = logR.reduce((a, b) => a + b, 0) / logR.length;
  const mv = logVh.reduce((a, b) => a + b, 0) / logVh.length;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mr) * (logVh[i] - mv); sxx += (logR[i] - mr) ** 2; }
  return sxx > 0 ? sxy / sxx : 0;
}

function computeRadialProfile(haloV, Rdisk, Rmax) {
  const zones = [
    { name: 'inner', min: 0, max: Rdisk },
    { name: 'disk', min: Rdisk, max: 3 * Rdisk },
    { name: 'outer', min: 3 * Rdisk, max: Infinity },
  ];
  const result = {};
  let totalSupport = 0;
  for (const z of zones) {
    const pts = haloV.filter(h => h.r >= z.min && h.r < z.max);
    const support = pts.reduce((s, h) => s + h.V_halo * h.V_halo, 0);
    result[z.name] = { n: pts.length, support };
    totalSupport += support;
  }
  for (const z of zones) {
    result[z.name].fraction = totalSupport > 0 ? result[z.name].support / totalSupport : 0;
  }

  const innerFrac = result.inner.fraction;
  const outerFrac = result.outer.fraction;
  const diskFrac = result.disk.fraction;

  let transWidth = 0;
  if (haloV.length >= 5) {
    const fracs = haloV.map(h => h.V_obs > 0 ? h.V_halo / h.V_obs : 0);
    let r10 = haloV[0].r, r90 = haloV[haloV.length - 1].r;
    for (let i = 0; i < fracs.length; i++) {
      if (fracs[i] >= 0.1) { r10 = haloV[i].r; break; }
    }
    for (let i = 0; i < fracs.length; i++) {
      if (fracs[i] >= 0.9) { r90 = haloV[i].r; break; }
    }
    transWidth = Rdisk > 0 ? (r90 - r10) / Rdisk : 0;
  }

  return { innerFrac, diskFrac, outerFrac, transWidth };
}


const gals = [];
for (const g of d56.perGalaxy) {
  const sp = sparcMap[g.name] || sparcMap[normalize(g.name)];
  if (!sp || sp.Vflat <= 0) continue;
  const sr = resultsMap[g.name] || resultsMap[normalize(g.name)];
  if (!sr) continue;
  const rc = rcMap[normalize(g.name)];
  if (!rc || rc.length < 5) continue;
  const Vflat = sp.Vflat, Rdisk = sp.Rdisk;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  const Rmax = rc[rc.length - 1].r;
  const logMbar = Math.log10(Math.max(Mbar, 1));
  const logL36 = Math.log10(Math.max(sp.L36, 0.001));
  const logRdisk = Math.log10(Math.max(Rdisk, 0.01));
  const logSBdisk = Math.log10(Math.max(sp.SBdisk, 0.01));
  const logK = sr.models.dark_halo_linear.k > 0 ? Math.log10(sr.models.dark_halo_linear.k) : -5;
  const hR = sr.models.newtonian.mse > 0 ? Math.log10(Math.max(sr.models.newtonian.mse / Math.max(sr.models.dark_halo_linear.mse, 0.001), 0.01)) : 0;
  const V_Newt_Rmax = Math.sqrt(G * Mbar / Rmax);
  const dmFrac = Math.max(0, 1 - (V_Newt_Rmax / Vflat) ** 2);

  let rcSmooth = 0;
  if (rc.length >= 5) {
    const sm = [];
    for (let i = 2; i < rc.length - 2; i++) sm.push((rc[i - 2].v + rc[i - 1].v + rc[i].v + rc[i + 1].v + rc[i + 2].v) / 5);
    let ss = 0;
    for (let i = 0; i < sm.length; i++) ss += (rc[i + 2].v - sm[i]) ** 2;
    rcSmooth = 1 - Math.sqrt(ss / sm.length) / Vflat;
  }

  const haloFit = fitHaloProfile(rc, Mbar, Rdisk, Rmax);
  if (!haloFit) continue;

  gals.push({
    name: g.name, Vflat, Rdisk, Rmax, Mbar, rc,
    logVflat: Math.log10(Vflat), logL36, logRdisk, logMbar,
    logMHI: g.logMHI, morphT: sp.T, logA0: g.logA0,
    logSBdisk, envCode: g.envCode,
    logK, dmFrac, hR, rcSmooth,
    haloFit,
    concentration: haloFit.concentration,
    bestFamily: haloFit.bestFamily,
    familyAdvantage: haloFit.familyAdvantage,
    innerSlope: haloFit.innerSlope,
    innerFrac: haloFit.radialProfile.innerFrac,
    diskFrac: haloFit.radialProfile.diskFrac,
    outerFrac: haloFit.radialProfile.outerFrac,
    transWidth: haloFit.radialProfile.transWidth,
    nfwMSE: haloFit.nfwMSE,
    burkertMSE: haloFit.burkertMSE,
    pisoMSE: haloFit.pisoMSE,
  });
}

const N = gals.length;
const struct4 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]);
const struct6 = gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]);
const vfM = multiR2(struct4, gals.map(g => g.logVflat));
const a0M = multiR2(struct6, gals.map(g => g.logA0));
for (let i = 0; i < N; i++) { gals[i].VfR = vfM.residuals[i]; gals[i].a0R = a0M.residuals[i]; }
const sdVf = Math.sqrt(gals.reduce((a, g) => a + g.VfR ** 2, 0) / N);
const sdA0 = Math.sqrt(gals.reduce((a, g) => a + g.a0R ** 2, 0) / N);
for (let i = 0; i < N; i++) {
  gals[i].VfR_z = gals[i].VfR / sdVf;
  gals[i].a0R_z = gals[i].a0R / sdA0;
  gals[i].L_sum = gals[i].VfR_z + gals[i].a0R_z;
}
const ctrlV = gals.map(g => [g.logK, g.dmFrac, g.envCode]);
const Lresid = multiR2(ctrlV, gals.map(g => g.L_sum));
for (let i = 0; i < N; i++) gals[i].dq = Lresid.residuals[i];
gals.sort((a, b) => b.dq - a.dq);

const concRegress = multiR2(gals.map(g => [g.logMbar, g.logVflat]), gals.map(g => g.concentration));
for (let i = 0; i < N; i++) gals[i].concResid = concRegress.residuals[i];

console.log('='.repeat(70));
console.log('PROGRAM 7B: HALO SHAPE, NOT HALO STRENGTH');
console.log('Is H about halo SHAPE (concentration, core/cusp, radial');
console.log('redistribution) rather than halo AMPLITUDE?');
console.log('='.repeat(70));
console.log('\nN = ' + N + ' galaxies with halo profile fits');


console.log('\n\n' + '#'.repeat(70));
console.log('7B-1: CONCENTRATION AT FIXED MASS');
console.log('Are high-H galaxies over- or under-concentrated relative to');
console.log('mass-matched expectations?');
console.log('#'.repeat(70));

const targetNames = ['NGC2841', 'NGC3741', 'ESO563-G021'];
function matchScore(t, c) {
  const dVf = Math.abs(Math.log10(t.Vflat) - Math.log10(c.Vflat)) / 0.3;
  const dMb = Math.abs(t.logMbar - c.logMbar) / 0.5;
  const dRd = Math.abs(Math.log10(t.Rdisk) - Math.log10(c.Rdisk)) / 0.3;
  const dT = Math.abs(t.morphT - c.morphT) / 3;
  return Math.sqrt(dVf ** 2 + dMb ** 2 + dRd ** 2 + dT ** 2);
}

const targets = targetNames.map(n => gals.find(g => g.name === n)).filter(Boolean);
const controls = [];
const used = new Set(targetNames);
for (const target of targets) {
  const lowDQ = gals.filter(g => g.dq < 0 && !used.has(g.name));
  const scored = lowDQ.map(g => ({ g, score: matchScore(target, g) })).sort((a, b) => a.score - b.score);
  const best = scored[0].g;
  controls.push(best);
  used.add(best.name);
}
const pairs = targets.map((t, i) => ({ target: t, control: controls[i] }));

console.log('\n  RESIDUAL CONCENTRATION (c_NFW after controlling for mass and Vflat):');
console.log('  ' + '-'.repeat(90));
console.log('  ' + 'Pair'.padEnd(6) + 'Target'.padEnd(18) + 'c_NFW'.padEnd(8) + 'c_resid'.padEnd(10) + 'Control'.padEnd(18) + 'c_NFW'.padEnd(8) + 'c_resid'.padEnd(10) + 'Different?');
console.log('  ' + '-'.repeat(90));

let test1Pass = 0;
let test1Direction = [];
for (let i = 0; i < pairs.length; i++) {
  const t = pairs[i].target, c = pairs[i].control;
  const diff = t.concResid - c.concResid;
  const diffAbs = Math.abs(diff);
  const isDifferent = diffAbs > 2;
  if (isDifferent) test1Pass++;
  test1Direction.push(diff > 0 ? 'over' : 'under');
  console.log('  ' + (i + 1 + '').padEnd(6) + t.name.padEnd(18) + (t.concentration + '').padEnd(8) + t.concResid.toFixed(2).padEnd(10) + c.name.padEnd(18) + (c.concentration + '').padEnd(8) + c.concResid.toFixed(2).padEnd(10) + (isDifferent ? 'YES (' + (diff > 0 ? '+' : '') + diff.toFixed(1) + ')' : 'no'));
}

const r_dq_concResid = pearsonR(gals.map(g => g.dq), gals.map(g => g.concResid));
const t_stat_conc = r_dq_concResid * Math.sqrt((N - 2) / (1 - r_dq_concResid * r_dq_concResid));
console.log('\n  Full-sample: r(DQ, concResid) = ' + r_dq_concResid.toFixed(3) + ' (t=' + t_stat_conc.toFixed(2) + ', p' + (Math.abs(t_stat_conc) > 2 ? '<' : '>') + '0.05)');
console.log('  TEST 7B-1: ' + test1Pass + '/3 pairs show systematic concentration difference');
const consistentDirection = test1Direction.every(d => d === test1Direction[0]);
console.log('  Direction: ' + (consistentDirection ? 'CONSISTENT (' + test1Direction[0] + '-concentrated)' : 'MIXED'));
console.log('  ' + (test1Pass >= 2 && consistentDirection ? 'PASS' : 'FAIL'));


console.log('\n\n' + '#'.repeat(70));
console.log('7B-2: CORE/CUSP FAMILY TEST');
console.log('Do high-H galaxies prefer a different halo family?');
console.log('NFW (cuspy) vs Burkert (cored) vs pseudo-isothermal (cored)');
console.log('#'.repeat(70));

console.log('\n  BEST-FIT FAMILY BY GALAXY:');
console.log('  ' + '-'.repeat(95));
console.log('  ' + 'Pair'.padEnd(6) + 'Target'.padEnd(18) + 'Family'.padEnd(10) + 'NFW_MSE'.padEnd(12) + 'Burk_MSE'.padEnd(12) + 'pIso_MSE'.padEnd(12) + 'Control'.padEnd(18) + 'Family');
console.log('  ' + '-'.repeat(95));

let test2FamilyDiff = 0;
for (let i = 0; i < pairs.length; i++) {
  const t = pairs[i].target, c = pairs[i].control;
  const familyDiff = t.bestFamily !== c.bestFamily;
  if (familyDiff) test2FamilyDiff++;
  console.log('  ' + (i + 1 + '').padEnd(6) + t.name.padEnd(18) + t.bestFamily.padEnd(10) + t.nfwMSE.toFixed(2).padEnd(12) + t.burkertMSE.toFixed(2).padEnd(12) + t.pisoMSE.toFixed(2).padEnd(12) + c.name.padEnd(18) + c.bestFamily);
}

const familyCounts = { highH: {}, lowH: {} };
const quintileSize = Math.floor(N / 5);
const Q1 = gals.slice(0, quintileSize);
const Q5 = gals.slice(N - quintileSize);
for (const g of Q1) familyCounts.highH[g.bestFamily] = (familyCounts.highH[g.bestFamily] || 0) + 1;
for (const g of Q5) familyCounts.lowH[g.bestFamily] = (familyCounts.lowH[g.bestFamily] || 0) + 1;

console.log('\n  QUINTILE FAMILY DISTRIBUTION:');
console.log('  Q1 (high-H, n=' + Q1.length + '): ' + Object.entries(familyCounts.highH).map(([k, v]) => k + '=' + v).join(', '));
console.log('  Q5 (low-H,  n=' + Q5.length + '): ' + Object.entries(familyCounts.lowH).map(([k, v]) => k + '=' + v).join(', '));

const q1NFWfrac = (familyCounts.highH['NFW'] || 0) / Q1.length;
const q5NFWfrac = (familyCounts.lowH['NFW'] || 0) / Q5.length;
const q1CoredFrac = ((familyCounts.highH['Burkert'] || 0) + (familyCounts.highH['pIso'] || 0)) / Q1.length;
const q5CoredFrac = ((familyCounts.lowH['Burkert'] || 0) + (familyCounts.lowH['pIso'] || 0)) / Q5.length;

console.log('  NFW fraction:   Q1=' + (q1NFWfrac * 100).toFixed(0) + '% vs Q5=' + (q5NFWfrac * 100).toFixed(0) + '%');
console.log('  Cored fraction: Q1=' + (q1CoredFrac * 100).toFixed(0) + '% vs Q5=' + (q5CoredFrac * 100).toFixed(0) + '%');

const familyShift = Math.abs(q1NFWfrac - q5NFWfrac) > 0.15;
console.log('\n  TEST 7B-2: ' + test2FamilyDiff + '/3 pairs have different family preference');
console.log('  Quintile family shift (>15%): ' + (familyShift ? 'YES' : 'NO'));
console.log('  ' + (test2FamilyDiff >= 2 || familyShift ? 'PASS' : 'FAIL'));


console.log('\n\n' + '#'.repeat(70));
console.log('7B-3: RADIAL SUPPORT REDISTRIBUTION');
console.log('How is halo support distributed radially?');
console.log('(inner vs disk vs outer zone fractions)');
console.log('#'.repeat(70));

console.log('\n  RADIAL HALO SUPPORT FRACTIONS (V_halo^2 share by zone):');
console.log('  ' + '-'.repeat(110));
console.log('  ' + 'Pair'.padEnd(6) + 'Target'.padEnd(18) + 'inner%'.padEnd(10) + 'disk%'.padEnd(10) + 'outer%'.padEnd(10) + 'transW'.padEnd(10) + 'Control'.padEnd(18) + 'inner%'.padEnd(10) + 'disk%'.padEnd(10) + 'outer%'.padEnd(10) + 'transW');
console.log('  ' + '-'.repeat(110));

let test3Consistent = 0;
const redistributionPattern = [];
for (let i = 0; i < pairs.length; i++) {
  const t = pairs[i].target, c = pairs[i].control;
  const dInner = t.innerFrac - c.innerFrac;
  const dOuter = t.outerFrac - c.outerFrac;
  const pattern = dInner > 0.05 ? 'inner-heavy' : dOuter > 0.05 ? 'outer-heavy' : 'similar';
  redistributionPattern.push(pattern);
  if (pattern !== 'similar') test3Consistent++;
  console.log('  ' + (i + 1 + '').padEnd(6) + t.name.padEnd(18) + (t.innerFrac * 100).toFixed(1).padEnd(10) + (t.diskFrac * 100).toFixed(1).padEnd(10) + (t.outerFrac * 100).toFixed(1).padEnd(10) + t.transWidth.toFixed(2).padEnd(10) + c.name.padEnd(18) + (c.innerFrac * 100).toFixed(1).padEnd(10) + (c.diskFrac * 100).toFixed(1).padEnd(10) + (c.outerFrac * 100).toFixed(1).padEnd(10) + c.transWidth.toFixed(2));
}

console.log('\n  QUINTILE RADIAL PROFILE:');
const q1Inner = Q1.reduce((s, g) => s + g.innerFrac, 0) / Q1.length;
const q5Inner = Q5.reduce((s, g) => s + g.innerFrac, 0) / Q5.length;
const q1Outer = Q1.reduce((s, g) => s + g.outerFrac, 0) / Q1.length;
const q5Outer = Q5.reduce((s, g) => s + g.outerFrac, 0) / Q5.length;
const q1Disk = Q1.reduce((s, g) => s + g.diskFrac, 0) / Q1.length;
const q5Disk = Q5.reduce((s, g) => s + g.diskFrac, 0) / Q5.length;
const q1Trans = Q1.reduce((s, g) => s + g.transWidth, 0) / Q1.length;
const q5Trans = Q5.reduce((s, g) => s + g.transWidth, 0) / Q5.length;

console.log('  ' + 'Zone'.padEnd(15) + 'Q1 (high-H)'.padEnd(15) + 'Q5 (low-H)'.padEnd(15) + 'Diff');
console.log('  ' + '-'.repeat(55));
console.log('  ' + 'inner'.padEnd(15) + (q1Inner * 100).toFixed(1).padEnd(15) + (q5Inner * 100).toFixed(1).padEnd(15) + ((q1Inner - q5Inner > 0 ? '+' : '') + ((q1Inner - q5Inner) * 100).toFixed(1) + '%'));
console.log('  ' + 'disk'.padEnd(15) + (q1Disk * 100).toFixed(1).padEnd(15) + (q5Disk * 100).toFixed(1).padEnd(15) + ((q1Disk - q5Disk > 0 ? '+' : '') + ((q1Disk - q5Disk) * 100).toFixed(1) + '%'));
console.log('  ' + 'outer'.padEnd(15) + (q1Outer * 100).toFixed(1).padEnd(15) + (q5Outer * 100).toFixed(1).padEnd(15) + ((q1Outer - q5Outer > 0 ? '+' : '') + ((q1Outer - q5Outer) * 100).toFixed(1) + '%'));
console.log('  ' + 'transWidth'.padEnd(15) + q1Trans.toFixed(2).padEnd(15) + q5Trans.toFixed(2).padEnd(15) + ((q1Trans - q5Trans > 0 ? '+' : '') + (q1Trans - q5Trans).toFixed(2)));

const r_dq_innerFrac = pearsonR(gals.map(g => g.dq), gals.map(g => g.innerFrac));
const r_dq_outerFrac = pearsonR(gals.map(g => g.dq), gals.map(g => g.outerFrac));
const r_dq_transWidth = pearsonR(gals.map(g => g.dq), gals.map(g => g.transWidth));
const r_dq_innerSlope = pearsonR(gals.map(g => g.dq), gals.map(g => g.innerSlope));

console.log('\n  FULL-SAMPLE CORRELATIONS:');
console.log('  r(DQ, innerFrac)  = ' + r_dq_innerFrac.toFixed(3));
console.log('  r(DQ, outerFrac)  = ' + r_dq_outerFrac.toFixed(3));
console.log('  r(DQ, transWidth) = ' + r_dq_transWidth.toFixed(3));
console.log('  r(DQ, innerSlope) = ' + r_dq_innerSlope.toFixed(3));

const redistribConsistent = redistributionPattern.filter(p => p !== 'similar').length >= 2 &&
  redistributionPattern.filter(p => p === redistributionPattern.find(pp => pp !== 'similar')).length >= 2;
console.log('\n  TEST 7B-3: ' + test3Consistent + '/3 pairs show radial redistribution');
console.log('  Pattern consistent: ' + (redistribConsistent ? 'YES' : 'NO'));
console.log('  ' + (test3Consistent >= 2 && redistribConsistent ? 'PASS' : 'FAIL'));


console.log('\n\n' + '#'.repeat(70));
console.log('7B-4: SHAPE + QUIETNESS COUPLING');
console.log('Is H the combination of halo-shape + kinematic quietness?');
console.log('#'.repeat(70));

const r_dq_smooth = pearsonR(gals.map(g => g.dq), gals.map(g => g.rcSmooth));
const r_dq_hR = pearsonR(gals.map(g => g.dq), gals.map(g => g.hR));

const shapeQuietProduct = gals.map(g => g.concResid * g.rcSmooth);
const r_dq_sqProduct = pearsonR(gals.map(g => g.dq), shapeQuietProduct);

const innerSlopeQuiet = gals.map(g => g.innerSlope * g.rcSmooth);
const r_dq_isqProduct = pearsonR(gals.map(g => g.dq), innerSlopeQuiet);

const familyQuiet = gals.map(g => (g.bestFamily === 'NFW' ? 1 : 0) * g.rcSmooth);
const r_dq_fqProduct = pearsonR(gals.map(g => g.dq), familyQuiet);

console.log('\n  INDIVIDUAL vs PRODUCT CORRELATIONS WITH DQ:');
console.log('  ' + '-'.repeat(60));
console.log('  ' + 'Variable'.padEnd(30) + 'r(DQ,X)'.padEnd(12) + 'Significant?');
console.log('  ' + '-'.repeat(60));

const test4Metrics = [
  { name: 'rcSmoothness (alone)', r: r_dq_smooth },
  { name: 'haloResponse (alone)', r: r_dq_hR },
  { name: 'concResid (alone)', r: r_dq_concResid },
  { name: 'innerSlope (alone)', r: r_dq_innerSlope },
  { name: 'concResid x smooth', r: r_dq_sqProduct },
  { name: 'innerSlope x smooth', r: r_dq_isqProduct },
  { name: 'NFW_flag x smooth', r: r_dq_fqProduct },
];

let bestProduct = null;
let bestProductR = 0;
for (const met of test4Metrics) {
  const t = met.r * Math.sqrt((N - 2) / (1 - met.r * met.r));
  const sig = Math.abs(t) > 2;
  console.log('  ' + met.name.padEnd(30) + ((met.r >= 0 ? '+' : '') + met.r.toFixed(3)).padEnd(12) + (sig ? 'yes (t=' + t.toFixed(2) + ')' : 'no'));
  if (met.name.includes('x') && Math.abs(met.r) > Math.abs(bestProductR)) {
    bestProduct = met.name;
    bestProductR = met.r;
  }
}

const productBeatsSingle = Math.abs(bestProductR) > Math.abs(r_dq_smooth) && Math.abs(bestProductR) > Math.abs(r_dq_concResid);
console.log('\n  Best product: ' + bestProduct + ' (r=' + bestProductR.toFixed(3) + ')');
console.log('  Product beats singles: ' + (productBeatsSingle ? 'YES — shape+quiet coupling exists' : 'NO — no synergy'));

console.log('\n  MATCHED PAIR SHAPE+QUIET:');
for (let i = 0; i < pairs.length; i++) {
  const t = pairs[i].target, c = pairs[i].control;
  console.log('  Pair ' + (i + 1) + ': ' + t.name + ' [smooth=' + t.rcSmooth.toFixed(3) + ', innerSlope=' + t.innerSlope.toFixed(3) + ', family=' + t.bestFamily + ']');
  console.log('         ' + c.name + ' [smooth=' + c.rcSmooth.toFixed(3) + ', innerSlope=' + c.innerSlope.toFixed(3) + ', family=' + c.bestFamily + ']');
  const bothDiff = Math.abs(t.innerSlope - c.innerSlope) > 0.1 && t.rcSmooth > 0.97;
  console.log('         Shape diff + quiet: ' + (bothDiff ? 'YES' : 'no'));
}

console.log('\n  TEST 7B-4: ' + (productBeatsSingle ? 'PASS' : 'FAIL'));


console.log('\n\n' + '='.repeat(70));
console.log('PROGRAM 7B GRAND VERDICT');
console.log('='.repeat(70));

const allTests = [
  { name: '7B-1: Concentration at fixed mass', pass: test1Pass >= 2 && consistentDirection },
  { name: '7B-2: Core/cusp family preference', pass: test2FamilyDiff >= 2 || familyShift },
  { name: '7B-3: Radial support redistribution', pass: test3Consistent >= 2 && redistribConsistent },
  { name: '7B-4: Shape + quietness coupling', pass: productBeatsSingle },
];

let totalPass = 0;
for (const t of allTests) {
  if (t.pass) totalPass++;
  console.log('  ' + (t.pass ? 'PASS' : 'FAIL').padEnd(6) + t.name);
}
console.log('  Total: ' + totalPass + '/4');

if (totalPass >= 3) {
  console.log('\n  VERDICT: HALO SHAPE CONFIRMED');
  console.log('  H is about halo SHAPE, not halo STRENGTH.');
  console.log('  The downstream coupling operates through concentration,');
  console.log('  profile family, or radial redistribution — not amplitude.');
} else if (totalPass >= 2) {
  console.log('\n  VERDICT: PARTIAL SUPPORT FOR HALO SHAPE');
  console.log('  Some evidence for shape-based coupling, but not conclusive.');
  console.log('  M2 needs refinement: gamma_hR should be reinterpreted.');
} else if (totalPass >= 1) {
  console.log('\n  VERDICT: WEAK / INCONCLUSIVE');
  console.log('  Halo shape hypothesis has minimal support.');
  console.log('  The gamma_hR coupling may be entirely epiphenomenal.');
} else {
  console.log('\n  VERDICT: HALO SHAPE ALSO FALSIFIED');
  console.log('  Neither halo amplitude (7A) nor halo shape (7B) explains');
  console.log('  the DQ-haloResponse coupling. gamma_hR is epiphenomenal.');
}


const outPath = path.join(__dirname, '..', 'public', 'program7b-halo-shape.json');
fs.writeFileSync(outPath, JSON.stringify({
  program: '7B',
  title: 'Halo Shape, Not Halo Strength',
  timestamp: new Date().toISOString(),
  N,
  pairs: pairs.map(p => ({
    target: { name: p.target.name, dq: p.target.dq, concentration: p.target.concentration, concResid: p.target.concResid, bestFamily: p.target.bestFamily, innerSlope: p.target.innerSlope, innerFrac: p.target.innerFrac, outerFrac: p.target.outerFrac, transWidth: p.target.transWidth, rcSmooth: p.target.rcSmooth },
    control: { name: p.control.name, dq: p.control.dq, concentration: p.control.concentration, concResid: p.control.concResid, bestFamily: p.control.bestFamily, innerSlope: p.control.innerSlope, innerFrac: p.control.innerFrac, outerFrac: p.control.outerFrac, transWidth: p.control.transWidth, rcSmooth: p.control.rcSmooth },
  })),
  tests: {
    t1_concentration: { pairPass: test1Pass, direction: test1Direction, consistent: consistentDirection, r_dq_concResid, pass: test1Pass >= 2 && consistentDirection },
    t2_family: { pairDiff: test2FamilyDiff, quintileFamilyCounts: familyCounts, familyShift, pass: test2FamilyDiff >= 2 || familyShift },
    t3_redistribution: { pairConsistent: test3Consistent, patterns: redistributionPattern, pass: test3Consistent >= 2 && redistribConsistent },
    t4_shapeQuiet: { productBeatsSingle, bestProduct, bestProductR, pass: productBeatsSingle },
  },
  correlations: { r_dq_concResid, r_dq_innerFrac, r_dq_outerFrac, r_dq_transWidth, r_dq_innerSlope, r_dq_smooth, r_dq_hR },
  totalPass,
  verdict: totalPass >= 3 ? 'HALO SHAPE CONFIRMED' : totalPass >= 2 ? 'PARTIAL SUPPORT' : totalPass >= 1 ? 'WEAK/INCONCLUSIVE' : 'HALO SHAPE ALSO FALSIFIED',
}, null, 2));
console.log('\nSaved: ' + outPath);
