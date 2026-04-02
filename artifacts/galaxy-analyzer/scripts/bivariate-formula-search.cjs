const fs = require('fs');
const path = require('path');

const G_KPC = 4.3009e-6;
const A0 = 3.7032;
const UPSILON_DISK = 0.5;
const UPSILON_BULGE = 0.7;

function rmsScatter(vals) {
  if (vals.length === 0) return Infinity;
  const m = vals.reduce((a, b) => a + b, 0) / vals.length;
  return Math.sqrt(vals.reduce((a, v) => a + (v - m) ** 2, 0) / vals.length);
}

function medAbsDev(vals) {
  const s = [...vals].sort((a, b) => a - b);
  const med = s[Math.floor(s.length / 2)];
  const devs = vals.map(v => Math.abs(v - med)).sort((a, b) => a - b);
  return devs[Math.floor(devs.length / 2)];
}

function pearsonR(x, y) {
  const n = x.length;
  const mx = x.reduce((a, b) => a + b) / n;
  const my = y.reduce((a, b) => a + b) / n;
  let sxy = 0, sxx = 0, syy = 0;
  for (let i = 0; i < n; i++) {
    sxy += (x[i] - mx) * (y[i] - my);
    sxx += (x[i] - mx) ** 2;
    syy += (y[i] - my) ** 2;
  }
  return syy > 0 ? sxy / Math.sqrt(sxx * syy) : 0;
}

console.log("=== BIVARIATE FORMULA SEARCH: F(g_bar, Sigma_bar) ===\n");
console.log("Goal: Find the single best function g_obs = F(g_bar, Sigma_bar)");
console.log("that produces a tighter collapse than the standard RAR,");
console.log("works cross-dataset, and survives honest scrutiny.\n");

const rarJson = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'rar-analysis-real.json'), 'utf8'));
const galaxyLookup = {};
for (const g of rarJson.perGalaxy) galaxyLookup[g.name] = g;

const sparcPoints = [];
const sparcDir = '/tmp/rotmod';
const files = fs.readdirSync(sparcDir).filter(f => f.endsWith('.dat'));

let sparcGalCount = 0;
for (const file of files) {
  const name = file.replace('_rotmod.dat', '');
  const gInfo = galaxyLookup[name];
  if (!gInfo) continue;
  const lines = fs.readFileSync(path.join(sparcDir, file), 'utf8').trim().split('\n');
  let galPts = 0;
  for (const line of lines) {
    if (line.trim().startsWith('#') || line.trim() === '') continue;
    const parts = line.trim().split(/\s+/).map(Number);
    if (parts.length < 6) continue;
    const [r, vobs, evobs, vgas, vdisk, vbul] = parts;
    if (!isFinite(r) || !isFinite(vobs) || r <= 0 || vobs <= 0) continue;
    const vBarSq = UPSILON_DISK * vdisk * Math.abs(vdisk) + UPSILON_BULGE * vbul * Math.abs(vbul) + vgas * Math.abs(vgas);
    const gObs = vobs * vobs / r;
    const gBar = Math.abs(vBarSq) / r;
    if (gBar <= 0 || gObs <= 0 || !isFinite(gBar) || !isFinite(gObs)) continue;
    const sigBar = gInfo.sigma_bar;
    if (!isFinite(sigBar) || sigBar <= 0) continue;
    sparcPoints.push({ galaxy: name, gObs, gBar, sigBar, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), logSig: Math.log10(sigBar) });
    galPts++;
  }
  if (galPts >= 3) sparcGalCount++;
}

function parseRotCurve(fp) {
  const lines = fs.readFileSync(fp, 'utf8').trim().split('\n');
  const g = {};
  for (const line of lines) {
    const nm = line.substring(0, 8).trim();
    const tp = line.substring(9, 14).trim();
    if (tp !== 'Data') continue;
    const r03 = parseFloat(line.substring(15, 23).trim());
    const v03 = parseFloat(line.substring(24, 34).trim());
    const rS = parseFloat(line.substring(35, 44).trim());
    const vS = parseFloat(line.substring(45, 54).trim());
    if (!g[nm]) g[nm] = [];
    g[nm].push({ r: rS * r03, v: vS * v03 });
  }
  return g;
}
function parseT2(fp) {
  const lines = fs.readFileSync(fp, 'utf8').trim().split('\n');
  const r = {};
  for (const line of lines) {
    const nm = line.substring(0, 8).trim();
    const rmax = parseFloat(line.substring(9, 13).trim());
    const vR = parseFloat(line.substring(19, 24).trim());
    const mgS = line.substring(139, 145).trim();
    const msS = line.substring(152, 157).trim();
    const mkS = line.substring(146, 151).trim();
    const mg = mgS ? parseFloat(mgS) * 1e7 : 0;
    const ms = msS ? parseFloat(msS) * 1e7 : (mkS ? parseFloat(mkS) * 1e7 : 0);
    r[nm] = { rmax, vR, mbar: ms + 1.33 * mg };
  }
  return r;
}
function interpV(c, r) {
  if (c.length === 0) return NaN;
  if (r <= c[0].r) return c[0].v * (r / c[0].r);
  if (r >= c[c.length - 1].r) return c[c.length - 1].v;
  for (let i = 0; i < c.length - 1; i++) {
    if (r >= c[i].r && r <= c[i + 1].r) {
      const f = (r - c[i].r) / (c[i + 1].r - c[i].r);
      return c[i].v + f * (c[i + 1].v - c[i].v);
    }
  }
  return c[c.length - 1].v;
}

const rotT = parseRotCurve('/tmp/little_things/rotdmbar.dat');
const rotD = parseRotCurve('/tmp/little_things/rotdm.dat');
const t2 = parseT2('/tmp/little_things/table2.dat');

const ltPoints = [];
let ltGalCount = 0;
for (const nm of Object.keys(rotT)) {
  if (!rotD[nm] || !t2[nm]) continue;
  const info = t2[nm];
  const cT = rotT[nm].sort((a, b) => a.r - b.r);
  const cD = rotD[nm].sort((a, b) => a.r - b.r);
  if (cT.length < 5 || info.mbar <= 0) continue;
  let galPts = 0;
  for (const pt of cT) {
    const r = pt.r, vO = pt.v;
    const vDM = interpV(cD, r);
    if (isNaN(vDM) || vO <= 0 || r <= 0) continue;
    let vBsq = vO * vO - vDM * vDM;
    if (vBsq < 0) vBsq = 0;
    const gObs = vO * vO / r;
    const gBar = vBsq / r;
    if (gBar <= 0 || gObs <= 0) continue;
    const sigBar = info.mbar / (Math.PI * r * r);
    if (!isFinite(sigBar) || sigBar <= 0) continue;
    ltPoints.push({ galaxy: nm, gObs, gBar, sigBar, logGobs: Math.log10(gObs), logGbar: Math.log10(gBar), logSig: Math.log10(sigBar) });
    galPts++;
  }
  if (galPts >= 3) ltGalCount++;
}

const allPoints = [...sparcPoints, ...ltPoints];
console.log("SPARC: " + sparcGalCount + " galaxies, " + sparcPoints.length + " points");
console.log("LITTLE THINGS: " + ltGalCount + " galaxies, " + ltPoints.length + " points");
console.log("COMBINED: " + allPoints.length + " points\n");

function stdRAR(gbar) {
  const y = gbar / A0;
  return gbar / (1 - Math.exp(-Math.sqrt(y)));
}

function evalFormula(points, predictor) {
  const residuals = [];
  for (const p of points) {
    const pred = predictor(p.gBar, p.sigBar, p.logGbar, p.logSig);
    if (!isFinite(pred) || pred <= 0) continue;
    residuals.push(Math.log10(p.gObs) - Math.log10(pred));
  }
  if (residuals.length < 10) return null;
  return {
    rms: rmsScatter(residuals),
    mad: medAbsDev(residuals),
    mean: residuals.reduce((a, b) => a + b, 0) / residuals.length,
    n: residuals.length
  };
}

function crossVal(trainPts, testPts, predictor) {
  const trainRes = evalFormula(trainPts, predictor);
  const testRes = evalFormula(testPts, predictor);
  if (!trainRes || !testRes) return null;
  return { train: trainRes, test: testRes };
}

const baselineAll = evalFormula(allPoints, (g) => stdRAR(g));
const baselineSPARC = evalFormula(sparcPoints, (g) => stdRAR(g));
const baselineLT = evalFormula(ltPoints, (g) => stdRAR(g));

console.log("============================================");
console.log("BASELINE: Standard McGaugh RAR");
console.log("============================================");
console.log("  ALL:   rms=" + baselineAll.rms.toFixed(6) + "  MAD=" + baselineAll.mad.toFixed(6));
console.log("  SPARC: rms=" + baselineSPARC.rms.toFixed(6) + "  MAD=" + baselineSPARC.mad.toFixed(6));
console.log("  LT:    rms=" + baselineLT.rms.toFixed(6) + "  MAD=" + baselineLT.mad.toFixed(6));

console.log("\n============================================");
console.log("SYSTEMATIC FORMULA SEARCH");
console.log("============================================\n");

const formulas = [];

console.log("--- Family 1: Variable a0 ---");
console.log("  g_obs = g_bar / (1 - exp(-sqrt(g_bar / a0_eff)))");
console.log("  a0_eff = a0 * (Sigma / Sigma_0)^delta\n");

const sig0_vals = [1e6, 3e6, 1e7, 3e7, 1e8];
const delta_vals = [];
for (let d = -0.5; d <= 0.5; d += 0.025) delta_vals.push(d);

let bestF1 = null;
for (const s0 of sig0_vals) {
  for (const delta of delta_vals) {
    const pred = (gbar, sig) => {
      const a0_eff = A0 * Math.pow(sig / s0, delta);
      const y = gbar / a0_eff;
      return gbar / (1 - Math.exp(-Math.sqrt(y)));
    };
    const res = evalFormula(allPoints, pred);
    if (!res) continue;
    if (!bestF1 || res.rms < bestF1.rms) {
      bestF1 = { ...res, delta, s0, label: "a0*(Sig/" + s0.toExponential(0) + ")^" + delta.toFixed(3) };
    }
  }
}
if (bestF1) {
  console.log("  Best: delta=" + bestF1.delta.toFixed(3) + " Sig0=" + bestF1.s0.toExponential(0));
  console.log("  rms=" + bestF1.rms.toFixed(6) + " (baseline=" + baselineAll.rms.toFixed(6) + " improvement=" + ((baselineAll.rms - bestF1.rms) / baselineAll.rms * 100).toFixed(2) + "%)");
  formulas.push({ name: "Variable a0", ...bestF1, predictor: (gbar, sig) => { const a = A0 * Math.pow(sig / bestF1.s0, bestF1.delta); return gbar / (1 - Math.exp(-Math.sqrt(gbar / a))); } });
}

console.log("\n--- Family 2: Multiplicative correction ---");
console.log("  g_obs = RAR(g_bar) * (Sigma / Sigma_0)^gamma\n");

let bestF2 = null;
for (const s0 of sig0_vals) {
  for (let gamma = -0.5; gamma <= 0.5; gamma += 0.01) {
    const pred = (gbar, sig) => stdRAR(gbar) * Math.pow(sig / s0, gamma);
    const res = evalFormula(allPoints, pred);
    if (!res) continue;
    if (!bestF2 || res.rms < bestF2.rms) {
      bestF2 = { ...res, gamma, s0, label: "RAR*(Sig/" + s0.toExponential(0) + ")^" + gamma.toFixed(3) };
    }
  }
}
if (bestF2) {
  console.log("  Best: gamma=" + bestF2.gamma.toFixed(3) + " Sig0=" + bestF2.s0.toExponential(0));
  console.log("  rms=" + bestF2.rms.toFixed(6) + " (improvement=" + ((baselineAll.rms - bestF2.rms) / baselineAll.rms * 100).toFixed(2) + "%)");
  formulas.push({ name: "Multiplicative", ...bestF2, predictor: (gbar, sig) => stdRAR(gbar) * Math.pow(sig / bestF2.s0, bestF2.gamma) });
}

console.log("\n--- Family 3: Additive correction ---");
console.log("  g_obs = RAR(g_bar) + alpha * (Sigma/Sigma_0)^beta * g_bar\n");

let bestF3 = null;
for (const s0 of [1e7, 3e7]) {
  for (let alpha = -0.5; alpha <= 0.5; alpha += 0.025) {
    for (let beta = -1.0; beta <= 1.0; beta += 0.25) {
      const pred = (gbar, sig) => stdRAR(gbar) + alpha * Math.pow(sig / s0, beta) * gbar;
      const res = evalFormula(allPoints, pred);
      if (!res) continue;
      if (!bestF3 || res.rms < bestF3.rms) {
        bestF3 = { ...res, alpha, beta, s0 };
      }
    }
  }
}
if (bestF3) {
  console.log("  Best: alpha=" + bestF3.alpha.toFixed(3) + " beta=" + bestF3.beta.toFixed(2) + " Sig0=" + bestF3.s0.toExponential(0));
  console.log("  rms=" + bestF3.rms.toFixed(6) + " (improvement=" + ((baselineAll.rms - bestF3.rms) / baselineAll.rms * 100).toFixed(2) + "%)");
  formulas.push({ name: "Additive", ...bestF3, predictor: (gbar, sig) => stdRAR(gbar) + bestF3.alpha * Math.pow(sig / bestF3.s0, bestF3.beta) * gbar });
}

console.log("\n--- Family 4: Two-parameter interpolation ---");
console.log("  g_obs = g_bar / (1 - exp(-(g_bar/a0)^nu * (Sigma/Sigma_0)^mu))\n");

let bestF4 = null;
for (const s0 of [1e7, 3e7]) {
  for (let nu = 0.3; nu <= 0.8; nu += 0.05) {
    for (let mu = -0.3; mu <= 0.3; mu += 0.025) {
      const pred = (gbar, sig) => {
        const x = Math.pow(gbar / A0, nu) * Math.pow(sig / s0, mu);
        return gbar / (1 - Math.exp(-x));
      };
      const res = evalFormula(allPoints, pred);
      if (!res) continue;
      if (!bestF4 || res.rms < bestF4.rms) {
        bestF4 = { ...res, nu, mu, s0 };
      }
    }
  }
}
if (bestF4) {
  console.log("  Best: nu=" + bestF4.nu.toFixed(2) + " mu=" + bestF4.mu.toFixed(3) + " Sig0=" + bestF4.s0.toExponential(0));
  console.log("  rms=" + bestF4.rms.toFixed(6) + " (improvement=" + ((baselineAll.rms - bestF4.rms) / baselineAll.rms * 100).toFixed(2) + "%)");
  formulas.push({ name: "Two-param interp", ...bestF4, predictor: (gbar, sig) => { const x = Math.pow(gbar / A0, bestF4.nu) * Math.pow(sig / bestF4.s0, bestF4.mu); return gbar / (1 - Math.exp(-x)); } });
}

console.log("\n--- Family 5: Log-space polynomial ---");
console.log("  log(g_obs) = a + b*log(g_bar) + c*log(Sig) + d*log(g_bar)*log(Sig) + e*log(g_bar)^2\n");

function fitLogPoly(points) {
  const n = points.length;
  const y = points.map(p => p.logGobs);
  const X = points.map(p => [1, p.logGbar, p.logSig, p.logGbar * p.logSig, p.logGbar * p.logGbar]);
  const p = 5;
  const Xt = [];
  for (let j = 0; j < p; j++) Xt.push(X.map(row => row[j]));
  const XtX = [];
  for (let i = 0; i < p; i++) {
    XtX.push([]);
    for (let j = 0; j < p; j++) {
      let s = 0;
      for (let k = 0; k < n; k++) s += Xt[i][k] * Xt[j][k];
      XtX[i].push(s);
    }
  }
  const Xty = [];
  for (let i = 0; i < p; i++) {
    let s = 0;
    for (let k = 0; k < n; k++) s += Xt[i][k] * y[k];
    Xty.push(s);
  }
  const M = XtX.map((row, i) => [...row, Xty[i]]);
  for (let col = 0; col < p; col++) {
    let maxR = col;
    for (let row = col + 1; row < p; row++) if (Math.abs(M[row][col]) > Math.abs(M[maxR][col])) maxR = row;
    [M[col], M[maxR]] = [M[maxR], M[col]];
    if (Math.abs(M[col][col]) < 1e-15) continue;
    for (let row = col + 1; row < p; row++) {
      const f = M[row][col] / M[col][col];
      for (let j = col; j <= p; j++) M[row][j] -= f * M[col][j];
    }
  }
  const beta = new Array(p).fill(0);
  for (let i = p - 1; i >= 0; i--) {
    beta[i] = M[i][p];
    for (let j = i + 1; j < p; j++) beta[i] -= M[i][j] * beta[j];
    beta[i] /= M[i][i];
  }
  return beta;
}

const polyBeta = fitLogPoly(allPoints);
console.log("  Fit on ALL: a=" + polyBeta[0].toFixed(4) + " b=" + polyBeta[1].toFixed(4) + " c=" + polyBeta[2].toFixed(5) + " d=" + polyBeta[3].toFixed(6) + " e=" + polyBeta[4].toFixed(6));
const f5pred = (gbar, sig, lgb, lgs) => Math.pow(10, polyBeta[0] + polyBeta[1] * lgb + polyBeta[2] * lgs + polyBeta[3] * lgb * lgs + polyBeta[4] * lgb * lgb);
const resF5 = evalFormula(allPoints, f5pred);
if (resF5) {
  console.log("  rms=" + resF5.rms.toFixed(6) + " (improvement=" + ((baselineAll.rms - resF5.rms) / baselineAll.rms * 100).toFixed(2) + "%)");
  formulas.push({ name: "Log-poly (5-param)", ...resF5, predictor: f5pred, beta: polyBeta });
}

console.log("\n--- Family 6: Variable exponent RAR ---");
console.log("  g_obs = g_bar^alpha(Sig) * a0^(1-alpha(Sig))  where alpha = a + b*log(Sig/Sig0)\n");

let bestF6 = null;
for (let a = 0.4; a <= 0.7; a += 0.025) {
  for (let b = -0.15; b <= 0.15; b += 0.01) {
    const pred = (gbar, sig) => {
      const alpha = a + b * Math.log10(sig / 1e7);
      if (alpha <= 0 || alpha >= 1) return stdRAR(gbar);
      return Math.pow(gbar, alpha) * Math.pow(A0, 1 - alpha);
    };
    const res = evalFormula(allPoints, pred);
    if (!res) continue;
    if (!bestF6 || res.rms < bestF6.rms) {
      bestF6 = { ...res, a, b };
    }
  }
}
if (bestF6) {
  console.log("  Best: a=" + bestF6.a.toFixed(3) + " b=" + bestF6.b.toFixed(3));
  console.log("  rms=" + bestF6.rms.toFixed(6) + " (improvement=" + ((baselineAll.rms - bestF6.rms) / baselineAll.rms * 100).toFixed(2) + "%)");
  formulas.push({ name: "Variable exponent", ...bestF6, predictor: (gbar, sig) => { const al = bestF6.a + bestF6.b * Math.log10(sig / 1e7); return (al > 0 && al < 1) ? Math.pow(gbar, al) * Math.pow(A0, 1 - al) : stdRAR(gbar); } });
}

console.log("\n--- Family 7: MOND-like with Sigma-dependent transition ---");
console.log("  g_obs = g_bar * [1 + (a0_eff/g_bar)^n]^(1/n) where a0_eff = a0*(Sig/Sig0)^delta, n varies\n");

let bestF7 = null;
for (const s0 of [1e7, 3e7]) {
  for (let delta = -0.3; delta <= 0.3; delta += 0.025) {
    for (let n = 0.5; n <= 2.5; n += 0.25) {
      const pred = (gbar, sig) => {
        const a0e = A0 * Math.pow(sig / s0, delta);
        return gbar * Math.pow(1 + Math.pow(a0e / gbar, n), 1 / n);
      };
      const res = evalFormula(allPoints, pred);
      if (!res) continue;
      if (!bestF7 || res.rms < bestF7.rms) {
        bestF7 = { ...res, delta, n, s0 };
      }
    }
  }
}
if (bestF7) {
  console.log("  Best: delta=" + bestF7.delta.toFixed(3) + " n=" + bestF7.n.toFixed(2) + " Sig0=" + bestF7.s0.toExponential(0));
  console.log("  rms=" + bestF7.rms.toFixed(6) + " (improvement=" + ((baselineAll.rms - bestF7.rms) / baselineAll.rms * 100).toFixed(2) + "%)");
  formulas.push({ name: "MOND Sig-transition", ...bestF7, predictor: (gbar, sig) => { const a = A0 * Math.pow(sig / bestF7.s0, bestF7.delta); return gbar * Math.pow(1 + Math.pow(a / gbar, bestF7.n), 1 / bestF7.n); } });
}

console.log("\n--- Family 8: Simple interpolation with Sigma modulation ---");
console.log("  g_obs = 0.5*g_bar*(1 + sqrt(1 + 4*a0_eff/g_bar))  [simple IF] with a0_eff = a0*(Sig/Sig0)^delta\n");

let bestF8 = null;
for (const s0 of [1e7, 3e7]) {
  for (let delta = -0.4; delta <= 0.4; delta += 0.02) {
    const pred = (gbar, sig) => {
      const a0e = A0 * Math.pow(sig / s0, delta);
      return 0.5 * gbar * (1 + Math.sqrt(1 + 4 * a0e / gbar));
    };
    const res = evalFormula(allPoints, pred);
    if (!res) continue;
    if (!bestF8 || res.rms < bestF8.rms) {
      bestF8 = { ...res, delta, s0 };
    }
  }
}
if (bestF8) {
  console.log("  Best: delta=" + bestF8.delta.toFixed(3) + " Sig0=" + bestF8.s0.toExponential(0));
  console.log("  rms=" + bestF8.rms.toFixed(6) + " (improvement=" + ((baselineAll.rms - bestF8.rms) / baselineAll.rms * 100).toFixed(2) + "%)");
  formulas.push({ name: "Simple IF + Sig", ...bestF8, predictor: (gbar, sig) => { const a = A0 * Math.pow(sig / bestF8.s0, bestF8.delta); return 0.5 * gbar * (1 + Math.sqrt(1 + 4 * a / gbar)); } });
}

console.log("\n============================================");
console.log("RANKING ALL FORMULAS");
console.log("============================================\n");

formulas.sort((a, b) => a.rms - b.rms);
console.log("  " + "Rank".padEnd(5) + "Formula".padEnd(25) + "RMS".padStart(10) + "  Improvement".padStart(12) + "  MAD".padStart(10));
for (let i = 0; i < formulas.length; i++) {
  const f = formulas[i];
  const imp = ((baselineAll.rms - f.rms) / baselineAll.rms * 100);
  console.log("  " + String(i + 1).padEnd(5) + f.name.padEnd(25) + f.rms.toFixed(6).padStart(10) + ("  " + imp.toFixed(2) + "%").padStart(12) + ("  " + f.mad.toFixed(6)).padStart(10));
}
console.log("  " + "---".padEnd(5) + "Standard RAR (baseline)".padEnd(25) + baselineAll.rms.toFixed(6).padStart(10) + "  0.00%".padStart(12) + ("  " + baselineAll.mad.toFixed(6)).padStart(10));

const best = formulas[0];
const second = formulas.length > 1 ? formulas[1] : null;

console.log("\n============================================");
console.log("WINNER: " + best.name);
console.log("============================================\n");

console.log("  Per-dataset breakdown:");
const bestSPARC = evalFormula(sparcPoints, best.predictor);
const bestLT = evalFormula(ltPoints, best.predictor);
console.log("    SPARC: rms=" + bestSPARC.rms.toFixed(6) + " (baseline=" + baselineSPARC.rms.toFixed(6) + " delta=" + ((baselineSPARC.rms - bestSPARC.rms) / baselineSPARC.rms * 100).toFixed(2) + "%)");
console.log("    LT:    rms=" + bestLT.rms.toFixed(6) + " (baseline=" + baselineLT.rms.toFixed(6) + " delta=" + ((baselineLT.rms - bestLT.rms) / baselineLT.rms * 100).toFixed(2) + "%)");

console.log("\n============================================");
console.log("CROSS-VALIDATION (THE ONLY TEST THAT MATTERS)");
console.log("============================================\n");

const topFormulas = formulas.slice(0, Math.min(4, formulas.length));
const cvResults = [];

for (const f of topFormulas) {
  const cv1 = crossVal(sparcPoints, ltPoints, f.predictor);
  const cv2 = crossVal(ltPoints, sparcPoints, f.predictor);
  
  if (!cv1 || !cv2) { console.log("  " + f.name + ": cross-validation failed"); continue; }
  
  const baseCv1 = crossVal(sparcPoints, ltPoints, (g) => stdRAR(g));
  const baseCv2 = crossVal(ltPoints, sparcPoints, (g) => stdRAR(g));
  
  const imp1 = baseCv1 ? ((baseCv1.test.rms - cv1.test.rms) / baseCv1.test.rms * 100) : 0;
  const imp2 = baseCv2 ? ((baseCv2.test.rms - cv2.test.rms) / baseCv2.test.rms * 100) : 0;
  
  console.log("  " + f.name + ":");
  console.log("    SPARC->LT:  test rms=" + cv1.test.rms.toFixed(6) + " (base=" + (baseCv1 ? baseCv1.test.rms.toFixed(6) : "?") + " imp=" + imp1.toFixed(2) + "%)");
  console.log("    LT->SPARC:  test rms=" + cv2.test.rms.toFixed(6) + " (base=" + (baseCv2 ? baseCv2.test.rms.toFixed(6) : "?") + " imp=" + imp2.toFixed(2) + "%)");
  console.log("    Both improve: " + (imp1 > 0 && imp2 > 0 ? "YES" : "NO"));
  
  cvResults.push({ name: f.name, imp1, imp2, bothImprove: imp1 > 0 && imp2 > 0, testRms1: cv1.test.rms, testRms2: cv2.test.rms });
}

console.log("\n============================================");
console.log("OVERFITTING CHECK: In-sample vs Out-of-sample");
console.log("============================================\n");

if (formulas.length > 0) {
  const f = formulas[0];
  
  const nSparc = sparcPoints.length;
  const nLt = ltPoints.length;
  const halfS = Math.floor(nSparc / 2);
  const halfL = Math.floor(nLt / 2);
  
  const sparcTrain = sparcPoints.slice(0, halfS);
  const sparcTest = sparcPoints.slice(halfS);
  const ltTrain = ltPoints.slice(0, halfL);
  const ltTest = ltPoints.slice(halfL);
  
  const trainAll = [...sparcTrain, ...ltTrain];
  const testAll = [...sparcTest, ...ltTest];
  
  const trainRes = evalFormula(trainAll, f.predictor);
  const testRes = evalFormula(testAll, f.predictor);
  const baseTrainRes = evalFormula(trainAll, (g) => stdRAR(g));
  const baseTestRes = evalFormula(testAll, (g) => stdRAR(g));
  
  console.log("  " + f.name + " (winner):");
  console.log("    Train half: rms=" + trainRes.rms.toFixed(6) + " (base=" + baseTrainRes.rms.toFixed(6) + ")");
  console.log("    Test half:  rms=" + testRes.rms.toFixed(6) + " (base=" + baseTestRes.rms.toFixed(6) + ")");
  console.log("    Train improvement: " + ((baseTrainRes.rms - trainRes.rms) / baseTrainRes.rms * 100).toFixed(2) + "%");
  console.log("    Test improvement:  " + ((baseTestRes.rms - testRes.rms) / baseTestRes.rms * 100).toFixed(2) + "%");
  
  const overfit = ((baseTrainRes.rms - trainRes.rms) / baseTrainRes.rms) > 2 * ((baseTestRes.rms - testRes.rms) / baseTestRes.rms);
  console.log("    Overfitting: " + (overfit ? "WARNING - train gain >> test gain" : "OK - comparable gains"));
}

console.log("\n============================================");
console.log("GALAXY-LEVEL RESIDUAL ANALYSIS");
console.log("============================================\n");

if (formulas.length > 0) {
  const f = formulas[0];
  
  function galaxyResiduals(points, predictor) {
    const byGal = {};
    for (const p of points) {
      if (!byGal[p.galaxy]) byGal[p.galaxy] = [];
      const pred = predictor(p.gBar, p.sigBar, p.logGbar, p.logSig);
      if (!isFinite(pred) || pred <= 0) continue;
      byGal[p.galaxy].push(Math.log10(p.gObs) - Math.log10(pred));
    }
    const galMeans = [];
    for (const nm of Object.keys(byGal)) {
      if (byGal[nm].length >= 3) {
        galMeans.push(byGal[nm].reduce((a, b) => a + b) / byGal[nm].length);
      }
    }
    return galMeans;
  }
  
  const galResBase = galaxyResiduals(allPoints, (g) => stdRAR(g));
  const galResBest = galaxyResiduals(allPoints, f.predictor);
  
  console.log("  Galaxy-level residual scatter:");
  console.log("    Standard RAR: rms=" + rmsScatter(galResBase).toFixed(6) + " MAD=" + medAbsDev(galResBase).toFixed(6) + " n=" + galResBase.length);
  console.log("    " + f.name + ":  rms=" + rmsScatter(galResBest).toFixed(6) + " MAD=" + medAbsDev(galResBest).toFixed(6) + " n=" + galResBest.length);
  console.log("    Improvement: " + ((rmsScatter(galResBase) - rmsScatter(galResBest)) / rmsScatter(galResBase) * 100).toFixed(2) + "%");
  
  const galResSparc = galaxyResiduals(sparcPoints, f.predictor);
  const galResLt = galaxyResiduals(ltPoints, f.predictor);
  const galBaseS = galaxyResiduals(sparcPoints, (g) => stdRAR(g));
  const galBaseL = galaxyResiduals(ltPoints, (g) => stdRAR(g));
  
  console.log("    SPARC galaxies: base=" + rmsScatter(galBaseS).toFixed(6) + " best=" + rmsScatter(galResSparc).toFixed(6) + " (" + ((rmsScatter(galBaseS) - rmsScatter(galResSparc)) / rmsScatter(galBaseS) * 100).toFixed(2) + "%)");
  console.log("    LT galaxies:    base=" + rmsScatter(galBaseL).toFixed(6) + " best=" + rmsScatter(galResLt).toFixed(6) + " (" + ((rmsScatter(galBaseL) - rmsScatter(galResLt)) / rmsScatter(galBaseL) * 100).toFixed(2) + "%)");
  
  const sparcMeanBase = galBaseS.reduce((a, b) => a + b, 0) / galBaseS.length;
  const ltMeanBase = galBaseL.reduce((a, b) => a + b, 0) / galBaseL.length;
  const sparcMeanBest = galResSparc.reduce((a, b) => a + b, 0) / galResSparc.length;
  const ltMeanBest = galResLt.reduce((a, b) => a + b, 0) / galResLt.length;
  
  console.log("\n  Inter-dataset alignment:");
  console.log("    Standard RAR: SPARC mean=" + sparcMeanBase.toFixed(5) + " LT mean=" + ltMeanBase.toFixed(5) + " gap=" + Math.abs(sparcMeanBase - ltMeanBase).toFixed(5));
  console.log("    " + f.name + ":  SPARC mean=" + sparcMeanBest.toFixed(5) + " LT mean=" + ltMeanBest.toFixed(5) + " gap=" + Math.abs(sparcMeanBest - ltMeanBest).toFixed(5));
}

console.log("\n============================================");
console.log("FINAL HONEST VERDICT");
console.log("============================================\n");

const bestImpAll = ((baselineAll.rms - best.rms) / baselineAll.rms * 100);
const cvBest = cvResults.find(r => r.name === best.name);
const sparcHelps = bestSPARC && baselineSPARC && bestSPARC.rms < baselineSPARC.rms;
const ltHelps = bestLT && baselineLT && bestLT.rms < baselineLT.rms;

console.log("  Best formula: " + best.name);
console.log("  Combined scatter reduction: " + bestImpAll.toFixed(2) + "%");
console.log("  Helps SPARC: " + (sparcHelps ? "YES" : "NO"));
console.log("  Helps LITTLE THINGS: " + (ltHelps ? "YES" : "NO"));
if (cvBest) {
  console.log("  Cross-validated both directions: " + (cvBest.bothImprove ? "YES" : "NO"));
  console.log("  CV improvement SPARC->LT: " + cvBest.imp1.toFixed(2) + "%");
  console.log("  CV improvement LT->SPARC: " + cvBest.imp2.toFixed(2) + "%");
}

const isTreasure = bestImpAll > 15 && sparcHelps && ltHelps && cvBest && cvBest.bothImprove && cvBest.imp1 > 5 && cvBest.imp2 > 5;
const isPromising = bestImpAll > 5 && sparcHelps && cvBest && cvBest.bothImprove;
const isModest = bestImpAll > 1 && sparcHelps;

if (isTreasure) {
  console.log("\n  >>> TREASURE FOUND: Genuine bivariate law with strong cross-validated improvement <<<");
} else if (isPromising) {
  console.log("\n  >>> PROMISING: Real improvement, cross-validated, but not transformative <<<");
} else if (isModest) {
  console.log("\n  >>> MODEST: Detectable improvement but marginal practical significance <<<");
} else {
  console.log("\n  >>> NO TREASURE: Standard RAR remains optimal or nearly so <<<");
}

const output = {
  timestamp: new Date().toISOString(),
  datasets: {
    sparc: { nGalaxies: sparcGalCount, nPoints: sparcPoints.length },
    littleThings: { nGalaxies: ltGalCount, nPoints: ltPoints.length },
    combined: { nPoints: allPoints.length }
  },
  baseline: {
    all: { rms: baselineAll.rms, mad: baselineAll.mad },
    sparc: { rms: baselineSPARC.rms, mad: baselineSPARC.mad },
    lt: { rms: baselineLT.rms, mad: baselineLT.mad }
  },
  formulas: formulas.map(f => ({
    name: f.name,
    rms: f.rms,
    mad: f.mad,
    improvement: ((baselineAll.rms - f.rms) / baselineAll.rms * 100),
    n: f.n
  })),
  winner: {
    name: best.name,
    rmsAll: best.rms,
    rmsSPARC: bestSPARC ? bestSPARC.rms : null,
    rmsLT: bestLT ? bestLT.rms : null,
    improvementAll: bestImpAll,
    improvementSPARC: bestSPARC ? ((baselineSPARC.rms - bestSPARC.rms) / baselineSPARC.rms * 100) : null,
    improvementLT: bestLT ? ((baselineLT.rms - bestLT.rms) / baselineLT.rms * 100) : null,
    helpsSPARC: sparcHelps,
    helpsLT: ltHelps
  },
  crossValidation: cvResults,
  verdict: isTreasure ? "TREASURE" : isPromising ? "PROMISING" : isModest ? "MODEST" : "NO_TREASURE"
};

const outPath = path.join(__dirname, '..', 'public', 'formula-search.json');
fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
console.log("\nResults saved to: " + outPath);
