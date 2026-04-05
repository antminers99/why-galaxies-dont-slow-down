#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 95: WHAT IS THE HIDDEN INNER QUANTITY?');
console.log('  Which inner-derived RC property explains a₀(inner)→outerSlope?');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const VFLAT_BREAK = 70;

function mcgaughRAR(gbar, a0) {
  const y = gbar / a0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function fitA0fromPts(pts) {
  if (pts.length < 3) return { logA0: NaN, a0: NaN };
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
    pts.push({ r: pt.rad, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), vobs: pt.vobs, vgas: pt.vgas, vdisk: pt.vdisk, vbul: pt.vbul || 0 });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const fitAll = fitA0fromPts(sorted);
  if (!isFinite(fitAll.logA0) || fitAll.logA0 < 0.5 || fitAll.logA0 > 5.5) continue;

  const half = Math.floor(sorted.length / 2);
  const innerPts = sorted.slice(0, half);
  const outerPts = sorted.slice(half);

  const fitInner = fitA0fromPts(innerPts);

  const logR_outer = outerPts.map(p => Math.log10(p.r));
  const logV_outer = outerPts.map(p => Math.log10(p.vobs));
  const mr = mean(logR_outer), mv = mean(logV_outer);
  let sxy = 0, sxx = 0;
  for (let i = 0; i < logR_outer.length; i++) { sxy += (logR_outer[i] - mr) * (logV_outer[i] - mv); sxx += (logR_outer[i] - mr) ** 2; }
  const outerSlope = sxx > 0 ? sxy / sxx : 0;

  const innerLogR = innerPts.map(p => Math.log10(p.r));
  const innerLogV = innerPts.map(p => Math.log10(p.vobs));
  const imr = mean(innerLogR), imv = mean(innerLogV);
  let isxy = 0, isxx = 0;
  for (let i = 0; i < innerLogR.length; i++) { isxy += (innerLogR[i] - imr) * (innerLogV[i] - imv); isxx += (innerLogR[i] - imr) ** 2; }
  const innerSlope = isxx > 0 ? isxy / isxx : 0;

  const innerVobs = innerPts.map(p => p.vobs);
  const innerVmax = Math.max(...innerVobs);
  const innerVmean = mean(innerVobs);
  const innerVrise = innerPts.length > 1 ? (innerVobs[innerVobs.length - 1] - innerVobs[0]) / (innerPts[innerPts.length - 1].r - innerPts[0].r) : 0;

  const innerGbar = innerPts.map(p => Math.pow(10, p.logGbar));
  const innerGobs = innerPts.map(p => Math.pow(10, p.logGobs));
  const innerGbarRange = Math.max(...innerPts.map(p => p.logGbar)) - Math.min(...innerPts.map(p => p.logGbar));
  const innerGobsRange = Math.max(...innerPts.map(p => p.logGobs)) - Math.min(...innerPts.map(p => p.logGobs));

  const innerMassConc = innerPts.length > 2 ? (() => {
    const vSq = innerPts.map(p => p.vobs * p.vobs);
    const rArr = innerPts.map(p => p.r);
    const enc = vSq.map((v, i) => v * rArr[i]);
    const midIdx = Math.floor(innerPts.length / 2);
    return enc[midIdx] > 0 ? enc[enc.length - 1] / enc[midIdx] : NaN;
  })() : NaN;

  const innerBaryonFrac = (() => {
    const bf = [];
    for (const p of innerPts) {
      const vBarSq = UPSILON_DISK * p.vdisk * Math.abs(p.vdisk) + UPSILON_BULGE * p.vbul * Math.abs(p.vbul) + p.vgas * Math.abs(p.vgas);
      if (p.vobs > 0) bf.push(Math.abs(vBarSq) / (p.vobs * p.vobs));
    }
    return mean(bf);
  })();

  let transitionRadius = NaN;
  const a0_cgs = 1.2e-10;
  for (let i = 0; i < innerPts.length; i++) {
    const gb = Math.pow(10, innerPts[i].logGbar);
    if (gb <= a0_cgs && i > 0) {
      transitionRadius = innerPts[i].r;
      break;
    }
  }
  const logTransR = isFinite(transitionRadius) && transitionRadius > 0 ? Math.log10(transitionRadius) : NaN;

  const innerCurvature = (() => {
    if (innerPts.length < 5) return NaN;
    const n = innerPts.length;
    const x = innerPts.map(p => Math.log10(p.r));
    const y = innerPts.map(p => Math.log10(p.vobs));
    const mx2 = mean(x);
    const X = x.map(xi => [1, xi - mx2, (xi - mx2) ** 2]);
    const fit2 = olsRegress(X, y);
    return fit2 ? fit2.beta[2] : NaN;
  })();

  const innerRARresid = (() => {
    const resid = innerPts.map(p => {
      const gb = Math.pow(10, p.logGbar);
      const pred = mcgaughRAR(gb, fitInner.a0);
      return p.logGobs - Math.log10(pred > 0 ? pred : 1e-20);
    });
    return Math.sqrt(resid.reduce((s, r) => s + r * r, 0) / resid.length);
  })();

  const innerRARslope = (() => {
    const x = innerPts.map(p => p.logGbar);
    const y = innerPts.map(p => p.logGobs);
    const mx3 = mean(x), my3 = mean(y);
    let sxy3 = 0, sxx3 = 0;
    for (let i = 0; i < x.length; i++) { sxy3 += (x[i] - mx3) * (y[i] - my3); sxx3 += (x[i] - mx3) ** 2; }
    return sxx3 > 0 ? sxy3 / sxx3 : NaN;
  })();

  const innerDiscrepancy = (() => {
    const disc = innerPts.map(p => p.logGobs - p.logGbar);
    return mean(disc);
  })();

  const innerDiscSlope = (() => {
    const x = innerPts.map(p => p.logGbar);
    const y = innerPts.map(p => p.logGobs - p.logGbar);
    const mx4 = mean(x), my4 = mean(y);
    let sxy4 = 0, sxx4 = 0;
    for (let i = 0; i < x.length; i++) { sxy4 += (x[i] - mx4) * (y[i] - my4); sxx4 += (x[i] - mx4) ** 2; }
    return sxx4 > 0 ? sxy4 / sxx4 : NaN;
  })();

  const Sigma0 = t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk * t1.Rdisk : 1);
  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);

  allGalaxies.push({
    name, nPts: pts.length,
    logA0inner: fitInner.logA0,
    logA0all: fitAll.logA0,
    outerSlope, innerSlope,
    innerVmax, innerVmean, innerVrise,
    innerGbarRange, innerGobsRange,
    innerMassConc: isFinite(innerMassConc) ? innerMassConc : NaN,
    innerBaryonFrac,
    logTransR,
    innerCurvature: isFinite(innerCurvature) ? innerCurvature : NaN,
    innerRARresid, innerRARslope,
    innerDiscrepancy, innerDiscSlope,
    logSigma0: Math.log10(Sigma0 > 0 ? Sigma0 : 1e-3),
    logVflat: Math.log10(t1.Vflat > 0 ? t1.Vflat : 1),
    Vflat: t1.Vflat, fgas,
    logMHI: Math.log10(t1.MHI),
  });
}

const hi = allGalaxies.filter(g => g.Vflat >= VFLAT_BREAK && isFinite(g.logA0inner));
console.log('  High-regime with valid inner a₀: N=' + hi.length + '\n');

const structControls = [g => g.logSigma0, g => g.logVflat, g => g.fgas, g => g.logMHI];

const baselineRef = partialCorrMulti(g => g.logA0inner, g => g.outerSlope, structControls, hi);
console.log('  BASELINE: r(a₀_inner, outerSlope | structural) = ' + baselineRef.toFixed(3));
console.log('  This is what we need to explain.\n');

console.log('══════════════════════════════════════════════════════════════');
console.log('  PART A: WHICH INNER QUANTITY CORRELATES WITH outerSlope?');
console.log('══════════════════════════════════════════════════════════════\n');

const candidates = [
  { name: 'innerSlope', fn: g => g.innerSlope, label: 'Inner RC slope (logV vs logR)' },
  { name: 'innerCurvature', fn: g => g.innerCurvature, label: 'Inner RC curvature (quadratic term)' },
  { name: 'innerVrise', fn: g => g.innerVrise, label: 'Inner velocity rise (dV/dR)' },
  { name: 'logInnerVmax', fn: g => Math.log10(g.innerVmax > 0 ? g.innerVmax : 1), label: 'Inner Vmax' },
  { name: 'innerBaryonFrac', fn: g => g.innerBaryonFrac, label: 'Inner baryon fraction (Vbar²/Vobs²)' },
  { name: 'innerGbarRange', fn: g => g.innerGbarRange, label: 'Inner gbar range (dex)' },
  { name: 'innerRARslope', fn: g => g.innerRARslope, label: 'Inner RAR slope (logGobs vs logGbar)' },
  { name: 'innerRARresid', fn: g => g.innerRARresid, label: 'Inner RAR RMS residual' },
  { name: 'innerDiscrepancy', fn: g => g.innerDiscrepancy, label: 'Inner mass discrepancy (mean logGobs-logGbar)' },
  { name: 'innerDiscSlope', fn: g => g.innerDiscSlope, label: 'Inner discrepancy slope (d(Gobs-Gbar)/dGbar)' },
  { name: 'innerMassConc', fn: g => g.innerMassConc, label: 'Inner mass concentration (enc_last/enc_mid)' },
];

const validCandidates = [];
for (const c of candidates) {
  const vals = hi.map(c.fn).filter(v => isFinite(v));
  if (vals.length < hi.length * 0.7) {
    console.log('  ' + c.name + ': SKIPPED (too many NaN: ' + vals.length + '/' + hi.length + ')');
    continue;
  }
  const hiValid = hi.filter(g => isFinite(c.fn(g)));

  const rRawOuter = pearsonR(hiValid.map(c.fn), hiValid.map(g => g.outerSlope));
  const rRawA0 = pearsonR(hiValid.map(c.fn), hiValid.map(g => g.logA0inner));
  const rPartOuter = partialCorrMulti(c.fn, g => g.outerSlope, structControls, hiValid);
  const rPartA0 = partialCorrMulti(c.fn, g => g.logA0inner, structControls, hiValid);

  console.log('  ' + c.name + ' (' + c.label + ')');
  console.log('    raw r(X, outerSlope) = ' + rRawOuter.toFixed(3) + '   partial = ' + rPartOuter.toFixed(3));
  console.log('    raw r(X, a₀_inner)   = ' + rRawA0.toFixed(3) + '   partial = ' + rPartA0.toFixed(3));
  console.log();

  validCandidates.push({
    ...c,
    n: hiValid.length,
    rRawOuter, rRawA0, rPartOuter, rPartA0,
    absPartOuter: Math.abs(rPartOuter),
    absPartA0: Math.abs(rPartA0),
  });
}

validCandidates.sort((a, b) => b.absPartOuter - a.absPartOuter);
console.log('\n  RANKING by |partial r(X, outerSlope)|:');
for (let i = 0; i < validCandidates.length; i++) {
  const c = validCandidates[i];
  console.log('    ' + (i + 1) + '. ' + c.name + ': ' + c.rPartOuter.toFixed(3));
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART B: DOES ANY INNER QUANTITY ABSORB THE a₀→outerSlope');
console.log('  RESIDUAL? (The key mediation test)');
console.log('══════════════════════════════════════════════════════════════\n');

const mediationResults = [];
for (const c of validCandidates) {
  const hiValid = hi.filter(g => isFinite(c.fn(g)));
  const extCtrl = [...structControls, c.fn];
  const rAfter = partialCorrMulti(g => g.logA0inner, g => g.outerSlope, extCtrl, hiValid);
  const reduction = 1 - Math.abs(rAfter) / Math.abs(baselineRef);
  console.log('  + ' + c.name + ':');
  console.log('    r(a₀_inner, outerSlope | struct + ' + c.name + ') = ' + rAfter.toFixed(3));
  console.log('    Reduction: ' + (reduction * 100).toFixed(1) + '%');
  console.log();
  mediationResults.push({ name: c.name, rAfter: +rAfter.toFixed(3), reduction: +(reduction * 100).toFixed(1) });
}

mediationResults.sort((a, b) => b.reduction - a.reduction);
console.log('\n  RANKING by reduction of a₀→outerSlope residual:');
for (let i = 0; i < mediationResults.length; i++) {
  const m = mediationResults[i];
  console.log('    ' + (i + 1) + '. ' + m.name + ': r=' + m.rAfter + ' (' + m.reduction + '% reduction)');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART C: BEST INNER MODEL — combine top candidates');
console.log('══════════════════════════════════════════════════════════════\n');

const top3 = mediationResults.slice(0, 3).map(m => validCandidates.find(c => c.name === m.name));
if (top3.length >= 3 && top3[0] && top3[1] && top3[2]) {
  const hiAll3 = hi.filter(g => top3.every(c => isFinite(c.fn(g))));
  const allInnerCtrl = [...structControls, top3[0].fn, top3[1].fn, top3[2].fn];
  const rAllInner = partialCorrMulti(g => g.logA0inner, g => g.outerSlope, allInnerCtrl, hiAll3);
  console.log('  Top 3 inner candidates: ' + top3.map(c => c.name).join(', '));
  console.log('  r(a₀_inner, outerSlope | struct + top3) = ' + rAllInner.toFixed(3));
  console.log('  N = ' + hiAll3.length);

  const totalReduction = 1 - Math.abs(rAllInner) / Math.abs(baselineRef);
  console.log('  Total reduction: ' + (totalReduction * 100).toFixed(1) + '%\n');

  if (Math.abs(rAllInner) < 0.1) {
    console.log('  → Top 3 inner quantities FULLY ABSORB the a₀ residual');
    console.log('  → a₀(inner) is a PROXY for these inner RC properties');
  } else if (Math.abs(rAllInner) < 0.2) {
    console.log('  → Top 3 inner quantities absorb MOST of the a₀ residual');
    console.log('  → a₀(inner) is partially a proxy, partially something else');
  } else {
    console.log('  → Top 3 inner quantities only partially absorb the a₀ residual');
    console.log('  → a₀(inner) carries information BEYOND these inner summaries');
  }
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART D: DOES THE INNER QUANTITY ITSELF PREDICT outerSlope');
console.log('  BETTER THAN a₀(inner)?');
console.log('══════════════════════════════════════════════════════════════\n');

const bestMediator = mediationResults[0];
const bestC = validCandidates.find(c => c.name === bestMediator.name);
if (bestC) {
  const hiV = hi.filter(g => isFinite(bestC.fn(g)));
  const rBestDirect = partialCorrMulti(bestC.fn, g => g.outerSlope, structControls, hiV);
  const rA0Direct = partialCorrMulti(g => g.logA0inner, g => g.outerSlope, structControls, hiV);

  console.log('  Best mediator: ' + bestC.name);
  console.log('  r(' + bestC.name + ', outerSlope | struct) = ' + rBestDirect.toFixed(3));
  console.log('  r(a₀_inner, outerSlope | struct) = ' + rA0Direct.toFixed(3));

  if (Math.abs(rBestDirect) > Math.abs(rA0Direct)) {
    console.log('  → ' + bestC.name + ' is a BETTER predictor than a₀(inner) itself!');
    console.log('  → a₀(inner) may be carrying this information indirectly');
  } else {
    console.log('  → a₀(inner) remains the better predictor');
    console.log('  → a₀(inner) integrates multiple inner RC features');
  }
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  PART E: INNER MASS DISCREPANCY as halo proxy');
console.log('  Does inner (Gobs-Gbar) encode halo contribution?');
console.log('══════════════════════════════════════════════════════════════\n');

const discC = validCandidates.find(c => c.name === 'innerDiscrepancy');
const discSC = validCandidates.find(c => c.name === 'innerDiscSlope');
if (discC && discSC) {
  const hiV = hi.filter(g => isFinite(discC.fn(g)) && isFinite(discSC.fn(g)));
  const rDisc = partialCorrMulti(discC.fn, g => g.outerSlope, structControls, hiV);
  const rDiscS = partialCorrMulti(discSC.fn, g => g.outerSlope, structControls, hiV);

  console.log('  Inner mass discrepancy (mean) → outerSlope: partial r = ' + rDisc.toFixed(3));
  console.log('  Inner discrepancy slope → outerSlope: partial r = ' + rDiscS.toFixed(3));

  const bothCtrl = [...structControls, discC.fn, discSC.fn];
  const rAfterDisc = partialCorrMulti(g => g.logA0inner, g => g.outerSlope, bothCtrl, hiV);
  console.log('  r(a₀_inner, outerSlope | struct + disc + discSlope) = ' + rAfterDisc.toFixed(3));
  console.log('  (baseline: ' + baselineRef.toFixed(3) + ')');

  const haloReduction = 1 - Math.abs(rAfterDisc) / Math.abs(baselineRef);
  console.log('  Halo-proxy reduction: ' + (haloReduction * 100).toFixed(1) + '%');
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  FINAL VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const output = {
  phase: '95',
  title: 'What Is the Hidden Inner Quantity?',
  hiN: hi.length,
  baseline: +baselineRef.toFixed(3),
  candidateRankByOuterSlope: validCandidates.map(c => ({
    name: c.name, label: c.label,
    rPartOuter: +c.rPartOuter.toFixed(3),
    rPartA0: +c.rPartA0.toFixed(3),
  })),
  mediationRank: mediationResults,
};

if (top3.length >= 3 && top3[0] && top3[1] && top3[2]) {
  const hiAll3 = hi.filter(g => top3.every(c => isFinite(c.fn(g))));
  const allInnerCtrl = [...structControls, top3[0].fn, top3[1].fn, top3[2].fn];
  const rAllInner = partialCorrMulti(g => g.logA0inner, g => g.outerSlope, allInnerCtrl, hiAll3);
  output.top3combined = {
    names: top3.map(c => c.name),
    rAfterAll: +rAllInner.toFixed(3),
    reduction: +((1 - Math.abs(rAllInner) / Math.abs(baselineRef)) * 100).toFixed(1),
  };
}

fs.writeFileSync(path.join(__dirname, '..', 'public', 'phase95-inner-quantity.json'), JSON.stringify(output, null, 2));
console.log('  Saved: public/phase95-inner-quantity.json');
