#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 94: CIRCULARITY vs GENUINE EXTRA DYNAMICS');
console.log('  Is r(a₀, slope | structural) = -0.294 real or artifact?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0fromPts(pts) {
  let lo = 0.5, hi = 5.0;
  for (let s = 0; s < 200; s++) {
    const m1 = lo + (hi - lo) * 0.382, m2 = lo + (hi - lo) * 0.618;
    let c1 = 0, c2 = 0;
    for (const p of pts) {
      const gb = Math.pow(10, p.logGbar);
      c1 += (p.logGobs - Math.log10(mcgaughRAR(gb, Math.pow(10, m1)))) ** 2;
      c2 += (p.logGobs - Math.log10(mcgaughRAR(gb, Math.pow(10, m2)))) ** 2;
    }
    if (c1 < c2) hi = m2; else lo = m1;
  }
  return { logA0: (lo + hi) / 2, a0: Math.pow(10, (lo + hi) / 2) };
}

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function sd(a) { if (a.length < 2) return 0; const m = mean(a); return Math.sqrt(a.reduce((s, v) => s + (v - m) ** 2, 0) / (a.length - 1)); }
function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function olsRegress(X, y) {
  const n = y.length, p = X[0].length;
  if (n <= p + 1) return null;
  const Xt = Array.from({ length: p }, (_, i) => X.map(row => row[i]));
  const XtX = Array.from({ length: p }, (_, i) =>
    Array.from({ length: p }, (_, j) =>
      Xt[i].reduce((s, _, k) => s + Xt[i][k] * Xt[j][k], 0)));
  const Xty = Xt.map(col => col.reduce((s, v, k) => s + v * y[k], 0));
  const aug = XtX.map((row, i) => [...row, Xty[i]]);
  for (let i = 0; i < p; i++) {
    let maxRow = i;
    for (let k = i + 1; k < p; k++) if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) maxRow = k;
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];
    if (Math.abs(aug[i][i]) < 1e-12) return null;
    for (let k = 0; k < p; k++) {
      if (k === i) continue;
      const f = aug[k][i] / aug[i][i];
      for (let j = i; j <= p; j++) aug[k][j] -= f * aug[i][j];
    }
  }
  const beta = aug.map((row, i) => row[p] / row[i]);
  const yhat = X.map(row => row.reduce((s, v, j) => s + v * beta[j], 0));
  const residuals = y.map((v, i) => v - yhat[i]);
  const ss_res = residuals.reduce((s, r) => s + r * r, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, residuals, yhat };
}

function partialCorrMulti(xFn, yFn, zFns, data) {
  const xv = data.map(xFn), yv = data.map(yFn);
  const Z = data.map(g => [1, ...zFns.map(fn => fn(g))]);
  const fitXZ = olsRegress(Z, xv);
  const fitYZ = olsRegress(Z, yv);
  if (!fitXZ || !fitYZ) return NaN;
  return pearsonR(fitXZ.residuals, fitYZ.residuals);
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const T = parseInt(line.substring(12, 14).trim());
  const L36 = parseFloat(line.substring(40, 47).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const Q = parseInt(line.substring(112, 115).trim());
  table1[name] = { T, L36, Rdisk, MHI, Vflat, Q };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, vgas, vdisk, vbul });
}

const p25 = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase25-group-membership.json'), 'utf-8'));
const envLookup = {};
for (const g of p25.galaxyAssignments) envLookup[g.name] = g;

const allGalaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    const gObs = pt.vobs * pt.vobs / pt.rad;
    const gBar = Math.abs(vBarSq) / pt.rad;
    if (gBar <= 0 || gObs <= 0 || !isFinite(gBar) || !isFinite(gObs)) continue;
    pts.push({ r: pt.rad, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), vobs: pt.vobs, vgas: pt.vgas, vdisk: pt.vdisk });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);

  const fitAll = fitA0fromPts(sorted);
  if (!isFinite(fitAll.logA0) || fitAll.logA0 < 0.5 || fitAll.logA0 > 5.5) continue;

  const half = Math.floor(sorted.length / 2);
  const innerPts = sorted.slice(0, half);
  const outerPts = sorted.slice(half);

  const fitInner = fitA0fromPts(innerPts);
  const fitOuter = fitA0fromPts(outerPts);

  const third = Math.floor(sorted.length / 3);
  const firstThird = sorted.slice(0, third);
  const lastThird = sorted.slice(2 * third);
  const fitFirst = fitA0fromPts(firstThird);
  const fitLast = fitA0fromPts(lastThird);

  const logR = outerPts.map(p => Math.log10(p.r));
  const logV = outerPts.map(p => Math.log10(p.vobs));
  const mrr = mean(logR), mvv = mean(logV);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR.length; i++) { sxy += (logR[i] - mrr) * (logV[i] - mvv); sxx += (logR[i] - mrr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const innerLogR = innerPts.map(p => Math.log10(p.r));
  const innerLogV = innerPts.map(p => Math.log10(p.vobs));
  const imr = mean(innerLogR), imv = mean(innerLogV);
  let isxy = 0, isxx = 0;
  for (let i = 0; i < innerLogR.length; i++) { isxy += (innerLogR[i] - imr) * (innerLogV[i] - imv); isxx += (innerLogR[i] - imr) ** 2; }
  const innerSlope = isxx > 0 ? isxy / isxx : 0;

  const rarResid = sorted.map(p => {
    const gb = Math.pow(10, p.logGbar);
    const pred = mcgaughRAR(gb, fitAll.a0);
    return p.logGobs - Math.log10(pred > 0 ? pred : 1e-20);
  });
  let currentSign = Math.sign(rarResid[0]);
  let runLen = 1;
  const residRun = [];
  for (let i = 1; i < rarResid.length; i++) {
    const s = Math.sign(rarResid[i]);
    if (s === currentSign) runLen++;
    else { residRun.push(runLen); runLen = 1; currentSign = s; }
  }
  residRun.push(runLen);
  const meanRunLen = residRun.reduce((s, v) => s + v, 0) / residRun.length;

  const ev = envLookup[name];
  const envCode = ev ? (ev.env === 'cluster' ? 2 : ev.env === 'group' ? 1 : 0) : 0;
  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);

  const vMax = Math.max(...sorted.map(p => p.vobs));
  const vLast3 = mean(sorted.slice(-3).map(p => p.vobs));
  const flatness = vMax > 0 ? vLast3 / vMax : 1;

  const outerCV = sd(outerPts.map(p => p.vobs)) / mean(outerPts.map(p => p.vobs));

  const rarRMS = Math.sqrt(rarResid.reduce((s, r) => s + r * r, 0) / rarResid.length);

  allGalaxies.push({
    name, nPts: pts.length,
    logA0all: fitAll.logA0,
    logA0inner: fitInner.logA0,
    logA0outer: fitOuter.logA0,
    logA0first: fitFirst.logA0,
    logA0last: fitLast.logA0,
    logMHI: Math.log10(t1.MHI),
    logMeanRun: Math.log10(meanRunLen > 0.1 ? meanRunLen : 0.1),
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    Vflat: t1.Vflat, Q: t1.Q, T: t1.T,
    envCode, fgas,
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    outerSlope, innerSlope, flatness, outerCV, rarRMS,
  });
}

const hi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK);
console.log('  High-regime: N=' + hi.length + ' (≥8 pts)\n');

const structControls = [g => g.logSigma0, g => g.logVflat, g => g.fgas, g => g.logMHI];

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 1: BASELINE — reproduce Phase 93 residual');
console.log('══════════════════════════════════════════════════════════════\n');

const baseResid = partialCorrMulti(g => g.logA0all, g => g.outerSlope, structControls, hi);
console.log('  r(logA0_all, outerSlope | Σ₀,Vflat,fgas,MHI) = ' + baseResid.toFixed(3));
console.log('  (Phase 93 value: -0.294)\n');

console.log('══════════════════════════════════════════════════════════════');
console.log('  TEST 2: SPLIT-HALF — a₀ from inner, slope from outer');
console.log('  Eliminates shared-data circularity');
console.log('══════════════════════════════════════════════════════════════\n');

const rInnerA0_outerSlope_raw = pearsonR(hi.map(g => g.logA0inner), hi.map(g => g.outerSlope));
const rInnerA0_outerSlope_ctrl = partialCorrMulti(g => g.logA0inner, g => g.outerSlope, structControls, hi);

console.log('  a₀(inner half) vs outerSlope:');
console.log('    Raw r = ' + rInnerA0_outerSlope_raw.toFixed(3));
console.log('    Partial r (| structural) = ' + rInnerA0_outerSlope_ctrl.toFixed(3));

const rOuterA0_innerSlope_raw = pearsonR(hi.map(g => g.logA0outer), hi.map(g => g.innerSlope));
const rOuterA0_innerSlope_ctrl = partialCorrMulti(g => g.logA0outer, g => g.innerSlope, structControls, hi);

console.log('\n  a₀(outer half) vs innerSlope:');
console.log('    Raw r = ' + rOuterA0_innerSlope_raw.toFixed(3));
console.log('    Partial r (| structural) = ' + rOuterA0_innerSlope_ctrl.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 3: THIRDS — a₀ from first third, slope from last third');
console.log('  Maximum separation of data');
console.log('══════════════════════════════════════════════════════════════\n');

const logR_last = hi.map(g => {
  return g.logA0last;
});

const rFirstA0_outerSlope_raw = pearsonR(hi.map(g => g.logA0first), hi.map(g => g.outerSlope));
const rFirstA0_outerSlope_ctrl = partialCorrMulti(g => g.logA0first, g => g.outerSlope, structControls, hi);

console.log('  a₀(first 1/3) vs outerSlope:');
console.log('    Raw r = ' + rFirstA0_outerSlope_raw.toFixed(3));
console.log('    Partial r (| structural) = ' + rFirstA0_outerSlope_ctrl.toFixed(3));

const rLastA0_innerSlope_raw = pearsonR(hi.map(g => g.logA0last), hi.map(g => g.innerSlope));
const rLastA0_innerSlope_ctrl = partialCorrMulti(g => g.logA0last, g => g.innerSlope, structControls, hi);

console.log('\n  a₀(last 1/3) vs innerSlope:');
console.log('    Raw r = ' + rLastA0_innerSlope_raw.toFixed(3));
console.log('    Partial r (| structural) = ' + rLastA0_innerSlope_ctrl.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 4: CONTROL FOR RC-SHAPE SUMMARIES');
console.log('  Add MeanRun, rarRMS, flatness, outerCV to controls');
console.log('══════════════════════════════════════════════════════════════\n');

const extControls = [
  ...structControls,
  g => g.logMeanRun,
  g => g.rarRMS,
  g => g.flatness,
  g => g.outerCV,
];

const extResid = partialCorrMulti(g => g.logA0all, g => g.outerSlope, extControls, hi);
console.log('  r(logA0, outerSlope | structural + RC-shape vars) = ' + extResid.toFixed(3));

const extResidNoCV = partialCorrMulti(g => g.logA0all, g => g.outerSlope,
  [...structControls, g => g.logMeanRun, g => g.rarRMS], hi);
console.log('  r(logA0, outerSlope | structural + MeanRun + rarRMS) = ' + extResidNoCV.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 5: INNER a₀ vs OUTER SLOPE with EXTENDED CONTROLS');
console.log('══════════════════════════════════════════════════════════════\n');

const splitExtCtrl = partialCorrMulti(g => g.logA0inner, g => g.outerSlope, extControls, hi);
console.log('  r(a₀_inner, outerSlope | all controls) = ' + splitExtCtrl.toFixed(3));

const splitExtCtrl2 = partialCorrMulti(g => g.logA0inner, g => g.outerSlope,
  [...structControls, g => g.logMeanRun], hi);
console.log('  r(a₀_inner, outerSlope | structural + MeanRun) = ' + splitExtCtrl2.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 6: PERMUTATION TEST — Is split-half residual significant?');
console.log('══════════════════════════════════════════════════════════════\n');

function permPartialCorr(xFn, yFn, zFns, data, nPerm) {
  const real = partialCorrMulti(xFn, yFn, zFns, data);
  if (isNaN(real)) return { real: NaN, p: 1 };
  let nBetter = 0;
  const yVals = data.map(yFn);
  for (let i = 0; i < nPerm; i++) {
    const shuffled = [...yVals].sort(() => Math.random() - 0.5);
    const shuffData = data.map((g, j) => ({ ...g, _shuffY: shuffled[j] }));
    const pr = partialCorrMulti(xFn, g => g._shuffY, zFns, shuffData);
    if (Math.abs(pr) >= Math.abs(real)) nBetter++;
  }
  return { real, p: nBetter / nPerm };
}

const permSplitHalf = permPartialCorr(g => g.logA0inner, g => g.outerSlope, structControls, hi, 2000);
console.log('  Split-half: r(a₀_inner, outerSlope | structural) = ' + permSplitHalf.real.toFixed(3) +
  ', perm p = ' + permSplitHalf.p.toFixed(3));

const permFull = permPartialCorr(g => g.logA0all, g => g.outerSlope, structControls, hi, 2000);
console.log('  Full:       r(a₀_all, outerSlope | structural) = ' + permFull.real.toFixed(3) +
  ', perm p = ' + permFull.p.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 7: CONSISTENCY — a₀(inner) vs a₀(outer)');
console.log('══════════════════════════════════════════════════════════════\n');

const rInOuter = pearsonR(hi.map(g => g.logA0inner), hi.map(g => g.logA0outer));
console.log('  r(a₀_inner, a₀_outer) = ' + rInOuter.toFixed(3));
console.log('  (High → a₀ is a stable galaxy property, not noise)');

const sdInner = sd(hi.map(g => g.logA0inner));
const sdOuter = sd(hi.map(g => g.logA0outer));
const sdAll = sd(hi.map(g => g.logA0all));
console.log('  sd(a₀_inner) = ' + sdInner.toFixed(3) + ', sd(a₀_outer) = ' + sdOuter.toFixed(3) + ', sd(a₀_all) = ' + sdAll.toFixed(3));

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  TEST 8: CORRELATION DECOMPOSITION');
console.log('  How much of the baseline -0.294 is from shared RC?');
console.log('══════════════════════════════════════════════════════════════\n');

console.log('  Baseline (same RC):      r(a₀_all, slope | struct) = ' + baseResid.toFixed(3));
console.log('  Split-half (diff RC):    r(a₀_inner, slope | struct) = ' + rInnerA0_outerSlope_ctrl.toFixed(3));
console.log();

const ratio = Math.abs(rInnerA0_outerSlope_ctrl) / Math.abs(baseResid);
console.log('  Ratio (split/full) = ' + ratio.toFixed(2));
if (ratio > 0.7) {
  console.log('  → Most of the residual SURVIVES split-half');
  console.log('  → This is GENUINE extra dynamics, NOT circularity');
} else if (ratio > 0.3) {
  console.log('  → PARTIAL circularity: some residual survives, some is shared-RC');
} else {
  console.log('  → The residual is MOSTLY circularity from shared RC data');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

let verdict;
if (ratio > 0.6 && permSplitHalf.p < 0.05) {
  verdict = 'GENUINE-DYNAMICS';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  GENUINE EXTRA DYNAMICS: a₀ carries real dynamical info     ║');
  console.log('  ║  beyond baryonic structure (survives split-half + perm)     ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else if (ratio > 0.3 && permSplitHalf.p < 0.1) {
  verdict = 'PARTIAL-GENUINE';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  PARTIALLY GENUINE: some real dynamical content survives    ║');
  console.log('  ║  split-half, but weaker than full-RC estimate              ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else if (ratio < 0.3) {
  verdict = 'CIRCULARITY';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  CIRCULARITY: the residual is mostly shared-RC artifact    ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
} else {
  verdict = 'INCONCLUSIVE';
  console.log('  ╔══════════════════════════════════════════════════════════════╗');
  console.log('  ║  INCONCLUSIVE: signal present but not statistically robust ║');
  console.log('  ╚══════════════════════════════════════════════════════════════╝');
}

const output = {
  phase: '94',
  title: 'Circularity vs Genuine Extra Dynamics',
  hiN: hi.length,
  baseline: { r: +baseResid.toFixed(3) },
  splitHalf: {
    innerA0_outerSlope: { raw: +rInnerA0_outerSlope_raw.toFixed(3), partial: +rInnerA0_outerSlope_ctrl.toFixed(3) },
    outerA0_innerSlope: { raw: +rOuterA0_innerSlope_raw.toFixed(3), partial: +rOuterA0_innerSlope_ctrl.toFixed(3) },
  },
  thirds: {
    firstA0_outerSlope: { raw: +rFirstA0_outerSlope_raw.toFixed(3), partial: +rFirstA0_outerSlope_ctrl.toFixed(3) },
    lastA0_innerSlope: { raw: +rLastA0_innerSlope_raw.toFixed(3), partial: +rLastA0_innerSlope_ctrl.toFixed(3) },
  },
  extendedControls: { r: +extResid.toFixed(3), rNoCV: +extResidNoCV.toFixed(3) },
  splitHalfExtended: { r: +splitExtCtrl.toFixed(3) },
  permutation: {
    splitHalf: { r: +permSplitHalf.real.toFixed(3), p: +permSplitHalf.p.toFixed(3) },
    full: { r: +permFull.real.toFixed(3), p: +permFull.p.toFixed(3) },
  },
  innerOuterConsistency: +rInOuter.toFixed(3),
  ratio: +ratio.toFixed(2),
  verdict,
};

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase94-circularity-test.json'), JSON.stringify(output, null, 2));
console.log('\n  Saved: public/phase94-circularity-test.json');
