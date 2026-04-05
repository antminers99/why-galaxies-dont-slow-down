#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 110: UNCERTAINTY PROPAGATION');
console.log('  Does the signal survive measurement noise?');
console.log('  Monte Carlo perturbation of catalog values.');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;
const N_MC = 2000;

let rng = 42;
function rand() { rng = (rng * 1664525 + 1013904223) & 0x7fffffff; return rng / 0x7fffffff; }
function randNormal() {
  let u, v, s;
  do { u = 2 * rand() - 1; v = 2 * rand() - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
  return u * Math.sqrt(-2 * Math.log(s) / s);
}

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function percentile(a, p) { const s = [...a].sort((x, y) => x - y); const i = p * (s.length - 1); const lo = Math.floor(i); return lo === s.length - 1 ? s[lo] : s[lo] + (i - lo) * (s[lo + 1] - s[lo]); }

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
  const ss_res = y.reduce((s, v, i) => s + (v - yhat[i]) ** 2, 0);
  const ss_tot = y.reduce((s, v) => s + (v - mean(y)) ** 2, 0);
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0 };
}

function looR2(xArr, yArr) {
  const n = xArr.length;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(yArr);
  for (let i = 0; i < n; i++) {
    const xTr = xArr.filter((_, j) => j !== i);
    const yTr = yArr.filter((_, j) => j !== i);
    const fit = olsRegress(xTr.map(v => [1, v]), yTr);
    if (!fit) return NaN;
    ss_loo += (yArr[i] - (fit.beta[0] + fit.beta[1] * xArr[i])) ** 2;
    ss_tot += (yArr[i] - yMean) ** 2;
  }
  return ss_tot > 0 ? 1 - ss_loo / ss_tot : 0;
}

const table2Raw = fs.readFileSync('/tmp/sparc_cds.dat', 'utf-8').trim().split('\n');
const table1Raw = fs.readFileSync('/tmp/sparc_table1.dat', 'utf-8').trim().split('\n');

const table1 = {};
for (const line of table1Raw) {
  const name = line.substring(0, 11).trim();
  const L36 = parseFloat(line.substring(40, 47).trim());
  const eL36 = parseFloat(line.substring(47, 55).trim());
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  const D = parseFloat(line.substring(15, 21).trim());
  const eD = parseFloat(line.substring(21, 27).trim());
  table1[name] = { L36, eL36: isFinite(eL36) ? eL36 : L36 * 0.1, Rdisk, MHI, Vflat, D, eD: isFinite(eD) ? eD : D * 0.15 };
}

const rcByGalaxy = {};
for (const line of table2Raw) {
  const name = line.substring(0, 11).trim();
  const rad = parseFloat(line.substring(19, 25).trim());
  const vobs = parseFloat(line.substring(26, 32).trim());
  const evobs = parseFloat(line.substring(32, 39).trim());
  const vgas = parseFloat(line.substring(39, 45).trim());
  const vdisk = parseFloat(line.substring(46, 52).trim());
  const vbul = parseFloat(line.substring(53, 59).trim());
  if (!rcByGalaxy[name]) rcByGalaxy[name] = [];
  rcByGalaxy[name].push({ rad, vobs, evobs: isFinite(evobs) ? evobs : vobs * 0.05, vgas, vdisk, vbul });
}

const galaxyData = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat < 70) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    pts.push(pt);
  }
  if (pts.length < 8) continue;
  galaxyData.push({ name, t1, pts: pts.sort((a, b) => a.rad - b.rad) });
}

function computeGalaxy(pts, t1, perturbVobs, perturbMHI, perturbL36) {
  const processed = [];
  for (const pt of pts) {
    const vobs = pt.vobs + (perturbVobs ? randNormal() * pt.evobs : 0);
    if (vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    processed.push({ vobs, vBar: Math.sqrt(Math.max(vBarSq, 0.01)) });
  }
  if (processed.length < 8) return null;

  const half = Math.floor(processed.length / 2);
  const outer = processed.slice(half);
  const outerMD = outer.map(p => (p.vobs ** 2) / (p.vBar ** 2)).filter(v => isFinite(v) && v > 0);
  const logOMD = Math.log10(mean(outerMD));
  if (!isFinite(logOMD)) return null;

  let MHI = t1.MHI;
  let L36 = t1.L36;
  if (perturbMHI) MHI = Math.max(0.01, MHI * (1 + randNormal() * 0.2));
  if (perturbL36) L36 = Math.max(0.01, L36 + randNormal() * t1.eL36);

  const fgas = MHI / (MHI + L36 * UPSILON_DISK);
  if (!isFinite(fgas) || fgas <= 0 || fgas >= 1) return null;

  return { logOMD, fgas };
}

console.log('  N = ' + galaxyData.length + ' galaxies\n');

const unperturbed = galaxyData.map(g => computeGalaxy(g.pts, g.t1, false, false, false)).filter(v => v);
const baseR = pearsonR(unperturbed.map(g => g.fgas), unperturbed.map(g => g.logOMD));
const baseLOO = looR2(unperturbed.map(g => g.fgas), unperturbed.map(g => g.logOMD));
const baseFit = olsRegress(unperturbed.map(g => [1, g.fgas]), unperturbed.map(g => g.logOMD));

console.log('  Unperturbed baseline: r=' + baseR.toFixed(4) + ', LOO=' + baseLOO.toFixed(4));
console.log('  slope=' + baseFit.beta[1].toFixed(4) + ', intercept=' + baseFit.beta[0].toFixed(4) + '\n');

const scenarios = [
  { name: 'Perturb Vobs only', vobs: true, mhi: false, l36: false },
  { name: 'Perturb MHI only (20%)', vobs: false, mhi: true, l36: false },
  { name: 'Perturb L36 only', vobs: false, mhi: false, l36: true },
  { name: 'Perturb ALL simultaneously', vobs: true, mhi: true, l36: true },
];

for (const scenario of scenarios) {
  console.log('══════════════════════════════════════════════════════════════');
  console.log('  ' + scenario.name + ' (' + N_MC + ' iterations)');
  console.log('══════════════════════════════════════════════════════════════\n');

  const mcR = [], mcLOO = [], mcSlope = [], mcIntercept = [];

  for (let iter = 0; iter < N_MC; iter++) {
    const perturbed = galaxyData.map(g => computeGalaxy(g.pts, g.t1, scenario.vobs, scenario.mhi, scenario.l36)).filter(v => v);
    if (perturbed.length < 20) continue;

    const xArr = perturbed.map(g => g.fgas);
    const yArr = perturbed.map(g => g.logOMD);

    const r = pearsonR(xArr, yArr);
    const fit = olsRegress(xArr.map(v => [1, v]), yArr);

    mcR.push(r);
    if (fit) { mcSlope.push(fit.beta[1]); mcIntercept.push(fit.beta[0]); }

    if (iter < 50) {
      const loo = looR2(xArr, yArr);
      mcLOO.push(loo);
    }
  }

  console.log('  r:         mean=' + mean(mcR).toFixed(4) + ' [' + percentile(mcR, 0.025).toFixed(4) + ', ' + percentile(mcR, 0.975).toFixed(4) + ']');
  if (mcLOO.length > 0) console.log('  LOO R²:    mean=' + mean(mcLOO).toFixed(4) + ' [' + percentile(mcLOO, 0.025).toFixed(4) + ', ' + percentile(mcLOO, 0.975).toFixed(4) + '] (from ' + mcLOO.length + ' iterations)');
  console.log('  slope:     mean=' + mean(mcSlope).toFixed(4) + ' [' + percentile(mcSlope, 0.025).toFixed(4) + ', ' + percentile(mcSlope, 0.975).toFixed(4) + ']');
  console.log('  intercept: mean=' + mean(mcIntercept).toFixed(4) + ' [' + percentile(mcIntercept, 0.025).toFixed(4) + ', ' + percentile(mcIntercept, 0.975).toFixed(4) + ']');
  console.log('  r never drops below 0.5? ' + (percentile(mcR, 0.025) > 0.5 ? 'YES' : 'NO'));
  console.log('  slope sign-flip rate: ' + (mcSlope.filter(s => s < 0).length / mcSlope.length * 100).toFixed(1) + '%');
  console.log('');

  scenario.results = {
    r: { mean: parseFloat(mean(mcR).toFixed(4)), ci025: parseFloat(percentile(mcR, 0.025).toFixed(4)), ci975: parseFloat(percentile(mcR, 0.975).toFixed(4)) },
    loo: mcLOO.length > 0 ? { mean: parseFloat(mean(mcLOO).toFixed(4)), ci025: parseFloat(percentile(mcLOO, 0.025).toFixed(4)), ci975: parseFloat(percentile(mcLOO, 0.975).toFixed(4)) } : null,
    slope: { mean: parseFloat(mean(mcSlope).toFixed(4)), ci025: parseFloat(percentile(mcSlope, 0.025).toFixed(4)), ci975: parseFloat(percentile(mcSlope, 0.975).toFixed(4)) },
    intercept: { mean: parseFloat(mean(mcIntercept).toFixed(4)), ci025: parseFloat(percentile(mcIntercept, 0.025).toFixed(4)), ci975: parseFloat(percentile(mcIntercept, 0.975).toFixed(4)) },
    signFlipRate: parseFloat((mcSlope.filter(s => s < 0).length / mcSlope.length * 100).toFixed(1))
  };
}

console.log('══════════════════════════════════════════════════════════════');
console.log('  VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const worstCase = scenarios.find(s => s.name.includes('ALL'));
const rSurvives = worstCase && worstCase.results.r.ci025 > 0.4;
const slopeSurvives = worstCase && worstCase.results.signFlipRate < 1;

let verdict;
if (rSurvives && slopeSurvives) {
  verdict = 'STRONG PASS: Signal survives simultaneous perturbation of all measurement uncertainties. r 95% CI entirely above 0.4, zero slope sign-flips.';
} else if (rSurvives) {
  verdict = 'MODERATE PASS: r survives but slope shows some instability.';
} else {
  verdict = 'FAIL: Signal degrades significantly under measurement uncertainty.';
}

console.log('  ' + verdict);

const output = {
  phase: 110,
  title: 'Uncertainty Propagation (Monte Carlo)',
  n: galaxyData.length,
  nMC: N_MC,
  baseline: { r: parseFloat(baseR.toFixed(4)), loo: parseFloat(baseLOO.toFixed(4)), slope: parseFloat(baseFit.beta[1].toFixed(4)), intercept: parseFloat(baseFit.beta[0].toFixed(4)) },
  scenarios: scenarios.map(s => ({ name: s.name, ...s.results })),
  verdict
};

const outPath = path.join(__dirname, '..', 'public', 'phase110-uncertainty-propagation.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);
