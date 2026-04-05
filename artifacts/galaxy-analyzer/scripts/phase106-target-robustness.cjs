#!/usr/bin/env node
const fs = require('fs');
const path = require('path');

console.log('======================================================================');
console.log('  PHASE 106: TARGET DEFINITION ROBUSTNESS');
console.log('  Does fgas survive when we change HOW we measure outer support?');
console.log('  6 target definitions. Same predictor. Same galaxies.');
console.log('======================================================================\n');

const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function mean(a) { return a.length ? a.reduce((s, v) => s + v, 0) / a.length : NaN; }
function median(a) { const s = [...a].sort((x, y) => x - y); const m = Math.floor(s.length / 2); return s.length % 2 ? s[m] : (s[m - 1] + s[m]) / 2; }

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
  return { beta, r2: ss_tot > 0 ? 1 - ss_res / ss_tot : 0, ss_res };
}

function pearsonR(x, y) {
  const n = x.length; if (n < 5) return NaN;
  const mx = mean(x), my = mean(y);
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (x[i] - mx) * (y[i] - my); sxx += (x[i] - mx) ** 2; syy += (y[i] - my) ** 2; }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

function looR2(xArr, yArr) {
  const n = xArr.length;
  let ss_loo = 0, ss_tot = 0;
  const yMean = mean(yArr);
  for (let i = 0; i < n; i++) {
    const xTrain = xArr.filter((_, j) => j !== i);
    const yTrain = yArr.filter((_, j) => j !== i);
    const fit = olsRegress(xTrain.map(v => [1, v]), yTrain);
    if (!fit) return NaN;
    const pred = fit.beta[0] + fit.beta[1] * xArr[i];
    ss_loo += (yArr[i] - pred) ** 2;
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
  const Rdisk = parseFloat(line.substring(71, 76).trim());
  const MHI = parseFloat(line.substring(86, 93).trim());
  const Vflat = parseFloat(line.substring(100, 105).trim());
  table1[name] = { L36, Rdisk, MHI, Vflat };
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

const galaxies = [];
for (const [name, rcPoints] of Object.entries(rcByGalaxy)) {
  const t1 = table1[name];
  if (!t1 || !t1.MHI || t1.MHI <= 0 || !t1.L36 || t1.L36 <= 0 || t1.Vflat < 70) continue;

  const pts = [];
  for (const pt of rcPoints) {
    if (pt.rad <= 0 || pt.vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * pt.vdisk * Math.abs(pt.vdisk) +
                   UPSILON_BULGE * (pt.vbul || 0) * Math.abs(pt.vbul || 0) +
                   pt.vgas * Math.abs(pt.vgas);
    pts.push({ r: pt.rad, vobs: pt.vobs, vBar: Math.sqrt(Math.max(vBarSq, 0.01)),
               vobsSq: pt.vobs ** 2, vBarSq: Math.max(vBarSq, 0.01) });
  }
  if (pts.length < 8) continue;

  const sorted = [...pts].sort((a, b) => a.r - b.r);
  const half = Math.floor(sorted.length / 2);
  const outerPts = sorted.slice(half);
  const rMax = sorted[sorted.length - 1].r;
  const rdisk = t1.Rdisk > 0 ? t1.Rdisk : 1;

  const fgas = t1.MHI / (t1.MHI + t1.L36 * UPSILON_DISK);
  if (!isFinite(fgas)) continue;

  const ratios = outerPts.map(p => p.vobsSq / p.vBarSq).filter(v => isFinite(v) && v > 0);
  const diffs = outerPts.map(p => p.vobsSq - p.vBarSq).filter(v => isFinite(v));
  if (ratios.length < 3) continue;

  const targets = {};

  targets.logMeanRatio = Math.log10(mean(ratios));

  targets.logMedianRatio = Math.log10(median(ratios));

  const lastPt = outerPts[outerPts.length - 1];
  const lastRatio = lastPt.vobsSq / lastPt.vBarSq;
  targets.logLastPoint = isFinite(lastRatio) && lastRatio > 0 ? Math.log10(lastRatio) : NaN;

  const target2Rd = sorted.filter(p => Math.abs(p.r - 2 * rdisk) / (2 * rdisk) < 0.3);
  if (target2Rd.length > 0) {
    const closest = target2Rd.sort((a, b) => Math.abs(a.r - 2 * rdisk) - Math.abs(b.r - 2 * rdisk))[0];
    const r2rd = closest.vobsSq / closest.vBarSq;
    targets.logAt2Rd = isFinite(r2rd) && r2rd > 0 ? Math.log10(r2rd) : NaN;
  } else {
    targets.logAt2Rd = NaN;
  }

  targets.meanDiff = mean(diffs);

  targets.linearMeanRatio = mean(ratios);

  const valid = Object.values(targets).some(v => isFinite(v));
  if (!valid) continue;

  galaxies.push({ name, fgas, targets });
}

console.log('  N = ' + galaxies.length + ' galaxies\n');

const targetDefs = [
  { name: 'log10(mean Vobs2/Vbar2) outer [BASELINE]', key: 'logMeanRatio' },
  { name: 'log10(median Vobs2/Vbar2) outer', key: 'logMedianRatio' },
  { name: 'log10(Vobs2/Vbar2) at last point', key: 'logLastPoint' },
  { name: 'log10(Vobs2/Vbar2) at 2*Rdisk', key: 'logAt2Rd' },
  { name: 'mean(Vobs2 - Vbar2) outer [linear diff]', key: 'meanDiff' },
  { name: 'mean(Vobs2/Vbar2) outer [linear ratio]', key: 'linearMeanRatio' },
];

console.log('══════════════════════════════════════════════════════════════');
console.log('  RESULTS BY TARGET DEFINITION');
console.log('══════════════════════════════════════════════════════════════\n');

const targetResults = [];

for (const td of targetDefs) {
  const valid = galaxies.filter(g => isFinite(g.targets[td.key]));
  if (valid.length < 20) {
    console.log('  ' + td.name);
    console.log('    SKIPPED — only ' + valid.length + ' valid galaxies\n');
    targetResults.push({ name: td.name, key: td.key, n: valid.length, skipped: true });
    continue;
  }

  const xArr = valid.map(g => g.fgas);
  const yArr = valid.map(g => g.targets[td.key]);

  const r = pearsonR(xArr, yArr);
  const fit = olsRegress(xArr.map(v => [1, v]), yArr);
  const r2 = fit ? fit.r2 : NaN;
  const loo = looR2(xArr, yArr);
  const slope = fit ? fit.beta[1] : NaN;
  const intercept = fit ? fit.beta[0] : NaN;

  console.log('  ' + td.name);
  console.log('    N=' + valid.length + '  r=' + r.toFixed(3) + '  R2=' + r2.toFixed(3) + '  LOO R2=' + loo.toFixed(3));
  console.log('    slope=' + slope.toFixed(4) + '  intercept=' + intercept.toFixed(4) + '\n');

  targetResults.push({
    name: td.name, key: td.key, n: valid.length,
    r: parseFloat(r.toFixed(3)), r2: parseFloat(r2.toFixed(3)),
    loo: parseFloat(loo.toFixed(3)),
    slope: parseFloat(slope.toFixed(4)), intercept: parseFloat(intercept.toFixed(4))
  });
}

console.log('══════════════════════════════════════════════════════════════');
console.log('  COMPETITOR CHECK: logSigma0 across same targets');
console.log('══════════════════════════════════════════════════════════════\n');

const Sigma0Arr = galaxies.map(g => {
  const t1 = table1[g.name];
  return Math.log10(t1.L36 / (t1.Rdisk > 0 ? t1.Rdisk ** 2 : 1));
});

for (const td of targetDefs) {
  const valid = galaxies.map((g, i) => ({ ...g, logSigma0: Sigma0Arr[i] })).filter(g => isFinite(g.targets[td.key]));
  if (valid.length < 20) continue;

  const xArr = valid.map(g => g.logSigma0);
  const yArr = valid.map(g => g.targets[td.key]);
  const r = pearsonR(xArr, yArr);
  const loo = looR2(xArr, yArr);

  const fgasResult = targetResults.find(t => t.key === td.key);
  const fgasR = fgasResult && !fgasResult.skipped ? fgasResult.r : NaN;
  const fgasLOO = fgasResult && !fgasResult.skipped ? fgasResult.loo : NaN;
  const winner = Math.abs(fgasLOO) > Math.abs(loo) ? 'fgas' : 'logSigma0';

  console.log('  ' + td.name);
  console.log('    logSigma0: r=' + r.toFixed(3) + '  LOO=' + loo.toFixed(3) + '  | fgas: r=' + fgasR.toFixed(3) + '  LOO=' + fgasLOO.toFixed(3) + '  -> winner: ' + winner);
}

console.log('\n══════════════════════════════════════════════════════════════');
console.log('  VERDICT');
console.log('══════════════════════════════════════════════════════════════\n');

const logTargets = targetResults.filter(t => !t.skipped && t.key.startsWith('log'));
const stableCount = logTargets.filter(t => t.loo > 0.3).length;
const totalLog = logTargets.length;

let verdict;
if (stableCount === totalLog && totalLog >= 3) {
  verdict = 'STRONG PASS: fgas is robust across all ' + totalLog + ' log-based target definitions (all LOO > 0.3).';
} else if (stableCount >= Math.ceil(totalLog * 0.6)) {
  verdict = 'MODERATE PASS: fgas survives ' + stableCount + '/' + totalLog + ' log-based targets.';
} else {
  verdict = 'FAIL: fgas collapses when target definition changes (' + stableCount + '/' + totalLog + ' survive).';
}

const linearTargets = targetResults.filter(t => !t.skipped && !t.key.startsWith('log'));
if (linearTargets.length > 0) {
  const linSurvive = linearTargets.filter(t => Math.abs(t.r) > 0.4).length;
  verdict += ' Linear targets: ' + linSurvive + '/' + linearTargets.length + ' survive (r>0.4).';
}

console.log('  ' + verdict);

const output = {
  phase: 106,
  title: 'Target Definition Robustness',
  n: galaxies.length,
  targetResults,
  verdict
};

const outPath = path.join(__dirname, '..', 'public', 'phase106-target-robustness.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log('\n  Results saved to: ' + outPath);
