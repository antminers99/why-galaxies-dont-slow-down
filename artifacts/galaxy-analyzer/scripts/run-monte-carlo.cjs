const fs = require('fs');
const path = require('path');

const G = 4.3009e-6;
const g_dag_kpc = 3.7032;
const Y_disk = 0.5, Y_bulge = 0.7;
const rotmodDir = '/tmp/rotmod';
const N_MC = 1000;

const raw = fs.readFileSync('/tmp/sparc_table.mrt', 'utf8');
const tlines = raw.split('\n');
let lastSep = -1;
for (let i = 0; i < tlines.length; i++) { if (tlines[i].startsWith('---')) lastSep = i; }
const galaxyProps = {};
for (let i = lastSep + 1; i < tlines.length; i++) {
  const parts = tlines[i].trim().split(/\s+/);
  if (parts.length < 18) continue;
  galaxyProps[parts[0]] = {
    dist: parseFloat(parts[2]), eDist: parseFloat(parts[3]),
    inc: parseFloat(parts[5]), eInc: parseFloat(parts[6]),
    Q: parseInt(parts[17])
  };
}

const files = fs.readdirSync(rotmodDir).filter(f => f.endsWith('_rotmod.dat'));
const rawData = {};
for (const file of files) {
  const name = file.replace('_rotmod.dat', '');
  const content = fs.readFileSync(path.join(rotmodDir, file), 'utf8');
  const dataLines = content.split('\n').filter(l => l.trim() && !l.startsWith('#'));
  const points = [];
  for (const line of dataLines) {
    const p = line.trim().split(/\s+/).map(Number);
    if (p.length >= 6 && !isNaN(p[0]) && p[0] > 0 && p[1] > 0) {
      points.push({ r: p[0], Vobs: p[1], errV: p[2], Vgas: p[3], Vdisk: p[4], Vbul: p[5] });
    }
  }
  if (points.length >= 3) rawData[name] = points;
}

function gaussRandom() {
  let u, v, s;
  do { u = Math.random() * 2 - 1; v = Math.random() * 2 - 1; s = u * u + v * v; } while (s >= 1 || s === 0);
  return u * Math.sqrt(-2 * Math.log(s) / s);
}

function linreg(xs, ys) {
  const n = xs.length; if (n < 3) return { b: 0, r: 0, R2: 0, rmse: 0, a: 0 };
  const mx = xs.reduce((s, v) => s + v, 0) / n, my = ys.reduce((s, v) => s + v, 0) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) { sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; syy += (ys[i] - my) ** 2; }
  const b = sxx > 0 ? sxy / sxx : 0, a = my - b * mx, r = (sxx > 0 && syy > 0) ? sxy / Math.sqrt(sxx * syy) : 0;
  const rmse = Math.sqrt(ys.reduce((s, y, i) => s + (y - (a + b * xs[i])) ** 2, 0) / n);
  return { a, b, r, R2: r * r, rmse };
}

function multireg(x1s, x2s, ys) {
  const n = x1s.length;
  const mx1 = x1s.reduce((s, v) => s + v, 0) / n, mx2 = x2s.reduce((s, v) => s + v, 0) / n, my = ys.reduce((s, v) => s + v, 0) / n;
  let s11 = 0, s12 = 0, s22 = 0, s1y = 0, s2y = 0, syy = 0;
  for (let i = 0; i < n; i++) { const d1 = x1s[i] - mx1, d2 = x2s[i] - mx2, dy = ys[i] - my; s11 += d1 * d1; s12 += d1 * d2; s22 += d2 * d2; s1y += d1 * dy; s2y += d2 * dy; syy += dy * dy; }
  const det = s11 * s22 - s12 * s12; if (Math.abs(det) < 1e-20) return { a: 0, b1: 0, b2: 0, R2: 0, rmse: Infinity };
  const b1 = (s22 * s1y - s12 * s2y) / det, b2 = (s11 * s2y - s12 * s1y) / det, a = my - b1 * mx1 - b2 * mx2;
  const ssr = ys.reduce((s, y, i) => s + (y - (a + b1 * x1s[i] + b2 * x2s[i])) ** 2, 0);
  return { a, b1, b2, R2: 1 - ssr / (syy > 0 ? syy : 1), rmse: Math.sqrt(ssr / n) };
}

function partialCorr(xs, ys, zs) { const rxy = linreg(xs, ys).r, rxz = linreg(xs, zs).r, ryz = linreg(ys, zs).r; const d = Math.sqrt((1 - rxz ** 2) * (1 - ryz ** 2)); return d > 0 ? (rxy - rxz * ryz) / d : 0; }

function processGalaxiesForMC(incPerts, distPerts) {
  const galaxies = [];
  for (const [name, points] of Object.entries(rawData)) {
    const props = galaxyProps[name];
    if (!props) continue;
    const incOrig = props.inc * Math.PI / 180;
    const incNew = (props.inc + (incPerts?.[name] || 0)) * Math.PI / 180;
    const distRatio = (props.dist + (distPerts?.[name] || 0)) / props.dist;
    const sinRatio = Math.sin(incOrig) / Math.sin(incNew);
    const processed = [];
    for (const p of points) {
      const Vobs_corr = p.Vobs * sinRatio;
      const r_corr = p.r * distRatio;
      const Vd2 = p.Vdisk < 0 ? -p.Vdisk * p.Vdisk : p.Vdisk * p.Vdisk;
      const Vb2 = p.Vbul < 0 ? -p.Vbul * p.Vbul : p.Vbul * p.Vbul;
      const Vg2 = p.Vgas < 0 ? -p.Vgas * p.Vgas : p.Vgas * p.Vgas;
      const Vbar2 = Y_disk * Vd2 + Y_bulge * Vb2 + Vg2;
      const g_obs = Vobs_corr * Vobs_corr / r_corr;
      const g_bar = Math.abs(Vbar2) / r_corr;
      if (g_bar > 0 && g_obs > 0) {
        const x = Math.sqrt(g_bar / g_dag_kpc), g_rar = g_bar / (1 - Math.exp(-x));
        processed.push({ r: r_corr, Vobs: Vobs_corr, Vbar2, g_obs, g_bar, g_rar,
          delta_rar: Math.log10(g_obs) - Math.log10(g_rar) });
      }
    }
    if (processed.length < 3) continue;
    const Vmax = Math.max(...processed.map(p => p.Vobs));
    const r_fid = 2 * (Vmax / 70);
    const ptsInFid = processed.filter(p => p.r <= r_fid);
    let M_bar_fid;
    if (ptsInFid.length > 0) { const lf = ptsInFid[ptsInFid.length - 1]; M_bar_fid = Math.abs(lf.Vbar2) * lf.r / G; }
    else { M_bar_fid = Math.abs(processed[0].Vbar2) * processed[0].r / G; }
    const sigma_bar = M_bar_fid / (Math.PI * r_fid * r_fid);
    const outerPts = processed.filter(p => p.r > r_fid);
    const meanDeltaOuter = outerPts.length > 0 ? outerPts.reduce((s, p) => s + Math.abs(p.delta_rar), 0) / outerPts.length : 0;
    if (sigma_bar > 0 && isFinite(meanDeltaOuter) && meanDeltaOuter !== 0) {
      galaxies.push({ name, Vmax, sigma_bar, meanDeltaOuter });
    }
  }
  return galaxies;
}

const mcSlopes = [], mcPartialR = [], mcDeltaAIC = [], mcRawR = [];
console.log('Running ' + N_MC + ' Monte Carlo iterations...');
for (let t = 0; t < N_MC; t++) {
  const incPerts = {}, distPerts = {};
  for (const [name, props] of Object.entries(galaxyProps)) {
    let newInc = props.inc + gaussRandom() * props.eInc;
    newInc = Math.max(10, Math.min(89, newInc));
    incPerts[name] = newInc - props.inc;
    let newDist = props.dist + gaussRandom() * props.eDist;
    newDist = Math.max(0.1, newDist);
    distPerts[name] = newDist - props.dist;
  }
  const galaxies = processGalaxiesForMC(incPerts, distPerts);
  const logSig = galaxies.map(g => Math.log10(g.sigma_bar));
  const deltas = galaxies.map(g => g.meanDeltaOuter);
  const logVmax = galaxies.map(g => Math.log10(g.Vmax));
  const fit = linreg(logSig, deltas);
  mcSlopes.push(fit.b);
  mcRawR.push(fit.r);
  const rPart = partialCorr(logSig, deltas, logVmax);
  mcPartialR.push(rPart);
  const mA = linreg(logVmax, deltas);
  const mB = multireg(logVmax, logSig, deltas);
  const n = galaxies.length;
  const aicA = n * Math.log(mA.rmse ** 2) + 4;
  const aicB = n * Math.log(mB.rmse ** 2) + 6;
  mcDeltaAIC.push(aicB - aicA);
  if ((t + 1) % 100 === 0) console.log('  Iteration ' + (t + 1) + '/' + N_MC);
}

const sorted_b = [...mcSlopes].sort((a, b) => a - b);
const sorted_pr = [...mcPartialR].sort((a, b) => a - b);
const sorted_aic = [...mcDeltaAIC].sort((a, b) => a - b);

const ci95_b = [sorted_b[Math.floor(0.025 * N_MC)], sorted_b[Math.floor(0.975 * N_MC)]];
const ci99_b = [sorted_b[Math.floor(0.005 * N_MC)], sorted_b[Math.floor(0.995 * N_MC)]];
const fracNeg = mcSlopes.filter(b => b < 0).length / N_MC;
const meanB = mcSlopes.reduce((s, v) => s + v, 0) / N_MC;
const stdB = Math.sqrt(mcSlopes.reduce((s, v) => s + (v - meanB) ** 2, 0) / N_MC);
const ci95_pr = [sorted_pr[Math.floor(0.025 * N_MC)], sorted_pr[Math.floor(0.975 * N_MC)]];
const meanPR = mcPartialR.reduce((s, v) => s + v, 0) / N_MC;
const stdPR = Math.sqrt(mcPartialR.reduce((s, v) => s + (v - meanPR) ** 2, 0) / N_MC);
const ci95_aic = [sorted_aic[Math.floor(0.025 * N_MC)], sorted_aic[Math.floor(0.975 * N_MC)]];
const meanAIC = mcDeltaAIC.reduce((s, v) => s + v, 0) / N_MC;
const stdAIC = Math.sqrt(mcDeltaAIC.reduce((s, v) => s + (v - meanAIC) ** 2, 0) / N_MC);
const fracAICneg = mcDeltaAIC.filter(v => v < 0).length / N_MC;
const meanRawR = mcRawR.reduce((s, v) => s + v, 0) / N_MC;
const pass = ci95_b[1] < 0;

console.log('\nSlope b: mean=' + meanB.toFixed(4) + ' std=' + stdB.toFixed(4));
console.log('95% CI: [' + ci95_b[0].toFixed(4) + ', ' + ci95_b[1].toFixed(4) + ']');
console.log('Fraction negative: ' + fracNeg.toFixed(3));
console.log('Partial r|Vmax: mean=' + meanPR.toFixed(4) + ' std=' + stdPR.toFixed(4));
console.log('ΔAIC: mean=' + meanAIC.toFixed(1) + ' std=' + stdAIC.toFixed(1));
console.log('PASS:', pass);

const binCount = 30;
const bMin = Math.min(...mcSlopes), bMax = Math.max(...mcSlopes);
const binWidth = (bMax - bMin) / binCount;
const histogram = Array.from({ length: binCount }, (_, i) => {
  const lo = bMin + i * binWidth;
  const hi = lo + binWidth;
  const count = mcSlopes.filter(b => b >= lo && (i === binCount - 1 ? b <= hi : b < hi)).length;
  return { binCenter: +((lo + hi) / 2).toFixed(4), count, lo: +lo.toFixed(4), hi: +hi.toFixed(4) };
});

const mcData = {
  nIterations: N_MC,
  perturbed: ['inclination', 'distance'],
  slope: { mean: meanB, std: stdB, ci95: ci95_b, ci99: ci99_b, fracNegative: fracNeg, values: mcSlopes.map(v => +v.toFixed(5)) },
  partialR: { mean: meanPR, std: stdPR, ci95: ci95_pr },
  deltaAIC: { mean: meanAIC, std: stdAIC, ci95: ci95_aic, fracNegative: fracAICneg },
  rawR: { mean: meanRawR },
  histogram,
  pass,
};

const existing = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
existing.monteCarlo = mcData;
fs.writeFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), JSON.stringify(existing, null, 2));
console.log('Saved to rar-analysis-real.json');
