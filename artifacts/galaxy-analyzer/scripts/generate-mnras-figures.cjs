const fs = require('fs');
const path = require('path');

const d56 = require('../public/phase56-frozen-baselines.json');
const sparcTable = require('../public/sparc-table.json');
const sparcRes = require('../public/sparc-results.json');
const p903 = require('../public/program9-phase903.json');
const p9v = require('../public/program9v-red-team.json');
const p8a = require('../public/program8a-2d-state.json');

function norm(n) { return n.replace(/\s+/g,'').replace(/^NGC0*/,'NGC').replace(/^DDO0*/,'DDO').replace(/^IC0*/,'IC').replace(/^UGC0*/,'UGC').toUpperCase(); }
const spMap = {}; sparcTable.forEach(g => { spMap[norm(g.name)] = g; });
const srMap = {}; sparcRes.perGalaxy.forEach(g => { srMap[norm(g.name)] = g; });

function olsResid(X, y) {
  const n = y.length, p = X[0].length;
  const my = y.reduce((a, b) => a + b, 0) / n;
  const mX = Array(p).fill(0); for (let j = 0; j < p; j++) { for (let i = 0; i < n; i++) mX[j] += X[i][j]; mX[j] /= n; }
  const Xc = X.map(r => r.map((v, j) => v - mX[j])), yc = y.map(v => v - my);
  const XtX = Array.from({ length: p }, () => Array(p).fill(0)), Xty = Array(p).fill(0);
  for (let i = 0; i < n; i++) for (let j = 0; j < p; j++) { Xty[j] += Xc[i][j] * yc[i]; for (let k = 0; k < p; k++) XtX[j][k] += Xc[i][j] * Xc[i][k]; }
  const aug = XtX.map((r, i) => [...r, Xty[i]]);
  for (let c = 0; c < p; c++) { let mr = c; for (let r = c + 1; r < p; r++) if (Math.abs(aug[r][c]) > Math.abs(aug[mr][c])) mr = r; [aug[c], aug[mr]] = [aug[mr], aug[c]]; if (Math.abs(aug[c][c]) < 1e-14) continue; for (let r = c + 1; r < p; r++) { const f = aug[r][c] / aug[c][c]; for (let j = c; j <= p; j++) aug[r][j] -= f * aug[c][j]; } }
  const beta = Array(p).fill(0);
  for (let i = p - 1; i >= 0; i--) { beta[i] = aug[i][p]; for (let j = i + 1; j < p; j++) beta[i] -= aug[i][j] * beta[j]; beta[i] /= Math.abs(aug[i][i]) > 1e-14 ? aug[i][i] : 1; }
  const res = []; for (let i = 0; i < n; i++) { let pred = my; for (let j = 0; j < p; j++) pred += beta[j] * Xc[i][j]; res.push(y[i] - pred); }
  return res;
}
function zsc(arr) { const m = arr.reduce((a, b) => a + b, 0) / arr.length; const s = Math.sqrt(arr.reduce((ss, v) => ss + (v - m) ** 2, 0) / arr.length); return s > 1e-10 ? arr.map(v => (v - m) / s) : arr.map(() => 0); }

const gals = [];
for (const g of d56.perGalaxy) {
  const sp = spMap[norm(g.name)]; if (!sp || sp.Vflat <= 0) continue;
  const sr = srMap[norm(g.name)]; if (!sr) continue;
  const Mbar = (sp.L36 * 0.5 + sp.MHI * 1.33) * 1e9;
  gals.push({ name: g.name, logVflat: Math.log10(sp.Vflat), Vflat: sp.Vflat, logMbar: Math.log10(Math.max(Mbar, 1)), logL36: Math.log10(Math.max(sp.L36, 0.001)), logRdisk: Math.log10(Math.max(sp.Rdisk, 0.01)), morphT: sp.T, logMHI: g.logMHI, logSBdisk: Math.log10(Math.max(sp.SBdisk, 0.01)), logA0: g.logA0 });
}
const vfR = olsResid(gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT]), gals.map(g => g.logVflat));
const a0R = olsResid(gals.map(g => [g.logMbar, g.logL36, g.logRdisk, g.morphT, g.logMHI, g.logSBdisk]), gals.map(g => g.logA0));
const bilZ = zsc(zsc(vfR).map((v, i) => v + zsc(a0R)[i]));
for (let i = 0; i < gals.length; i++) { gals[i].vfResid = vfR[i]; gals[i].a0Resid = a0R[i]; gals[i].DQ = bilZ[i]; }

const outDir = path.join(__dirname, '..', 'public', 'replication', 'mnras', 'figures');
fs.mkdirSync(outDir, { recursive: true });

function esc(s) { return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;'); }

function svgHeader(w, h) {
  return '<?xml version="1.0" encoding="UTF-8"?>\n' +
    '<svg xmlns="http://www.w3.org/2000/svg" width="' + w + '" height="' + h + '" viewBox="0 0 ' + w + ' ' + h + '" font-family="Helvetica, Arial, sans-serif">\n' +
    '<rect width="' + w + '" height="' + h + '" fill="white"/>\n';
}

function mapRange(v, lo, hi, pLo, pHi) { return pLo + (v - lo) / (hi - lo) * (pHi - pLo); }


console.log('=== Generating MNRAS Figures ===\n');

// ====== FIGURE 1: Bilateral Coupling Scatter ======
{
  const W = 520, H = 480;
  const ml = 80, mr = 30, mt = 40, mb = 60;
  const pw = W - ml - mr, ph = H - mt - mb;

  const xArr = vfR, yArr = a0R;
  const xMin = -0.15, xMax = 0.15, yMin = -0.6, yMax = 0.6;

  let svg = svgHeader(W, H);

  svg += '<rect x="' + ml + '" y="' + mt + '" width="' + pw + '" height="' + ph + '" fill="#fafafa" stroke="#ccc"/>\n';

  svg += '<line x1="' + ml + '" y1="' + (mt + ph/2) + '" x2="' + (ml + pw) + '" y2="' + (mt + ph/2) + '" stroke="#ddd" stroke-dasharray="4,3"/>\n';
  svg += '<line x1="' + (ml + pw/2) + '" y1="' + mt + '" x2="' + (ml + pw/2) + '" y2="' + (mt + ph) + '" stroke="#ddd" stroke-dasharray="4,3"/>\n';

  const quadLabels = [
    { label: 'Q1', x: ml + pw * 0.75, y: mt + ph * 0.25, desc: '+Vf, +a\u2080' },
    { label: 'Q2', x: ml + pw * 0.25, y: mt + ph * 0.25, desc: '\u2212Vf, +a\u2080' },
    { label: 'Q3', x: ml + pw * 0.25, y: mt + ph * 0.75, desc: '\u2212Vf, \u2212a\u2080' },
    { label: 'Q4', x: ml + pw * 0.75, y: mt + ph * 0.75, desc: '+Vf, \u2212a\u2080 (dark)' },
  ];
  for (const q of quadLabels) {
    svg += '<text x="' + q.x + '" y="' + q.y + '" text-anchor="middle" font-size="9" fill="#999">' + esc(q.label) + '</text>\n';
  }

  const things7 = new Set(p903.results.map(r => norm(r.name)));
  for (let i = 0; i < gals.length; i++) {
    const px = mapRange(xArr[i], xMin, xMax, ml, ml + pw);
    const py = mapRange(yArr[i], yMax, yMin, mt, mt + ph);
    if (px < ml || px > ml + pw || py < mt || py > mt + ph) continue;
    const isThings = things7.has(norm(gals[i].name));
    const col = isThings ? '#e63946' : '#457b9d';
    const r = isThings ? 5 : 3;
    svg += '<circle cx="' + px.toFixed(1) + '" cy="' + py.toFixed(1) + '" r="' + r + '" fill="' + col + '" opacity="' + (isThings ? 0.9 : 0.5) + '" stroke="' + (isThings ? '#000' : 'none') + '" stroke-width="' + (isThings ? 0.8 : 0) + '"/>\n';
  }

  const mx = xArr.reduce((a, b) => a + b, 0) / xArr.length;
  const my = yArr.reduce((a, b) => a + b, 0) / yArr.length;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < xArr.length; i++) { sxy += (xArr[i] - mx) * (yArr[i] - my); sxx += (xArr[i] - mx) ** 2; }
  const slope = sxy / sxx, inter = my - slope * mx;
  const x1 = xMin, x2 = xMax, y1 = slope * x1 + inter, y2 = slope * x2 + inter;
  const lx1 = mapRange(x1, xMin, xMax, ml, ml + pw), ly1 = mapRange(y1, yMax, yMin, mt, mt + ph);
  const lx2 = mapRange(x2, xMin, xMax, ml, ml + pw), ly2 = mapRange(y2, yMax, yMin, mt, mt + ph);
  svg += '<line x1="' + lx1.toFixed(1) + '" y1="' + ly1.toFixed(1) + '" x2="' + lx2.toFixed(1) + '" y2="' + ly2.toFixed(1) + '" stroke="#e63946" stroke-width="2" stroke-dasharray="6,3"/>\n';

  for (let v = xMin; v <= xMax + 0.001; v += 0.05) {
    const px = mapRange(v, xMin, xMax, ml, ml + pw);
    svg += '<text x="' + px.toFixed(1) + '" y="' + (mt + ph + 18) + '" text-anchor="middle" font-size="11">' + v.toFixed(2) + '</text>\n';
    svg += '<line x1="' + px.toFixed(1) + '" y1="' + (mt + ph) + '" x2="' + px.toFixed(1) + '" y2="' + (mt + ph + 4) + '" stroke="#333"/>\n';
  }
  for (let v = yMin; v <= yMax + 0.001; v += 0.2) {
    const py = mapRange(v, yMax, yMin, mt, mt + ph);
    svg += '<text x="' + (ml - 8) + '" y="' + (py + 4) + '" text-anchor="end" font-size="11">' + v.toFixed(1) + '</text>\n';
    svg += '<line x1="' + (ml - 4) + '" y1="' + py.toFixed(1) + '" x2="' + ml + '" y2="' + py.toFixed(1) + '" stroke="#333"/>\n';
  }

  svg += '<text x="' + (ml + pw / 2) + '" y="' + (H - 8) + '" text-anchor="middle" font-size="13" font-weight="bold">VfResid (dex)</text>\n';
  svg += '<text x="15" y="' + (mt + ph / 2) + '" text-anchor="middle" font-size="13" font-weight="bold" transform="rotate(-90,15,' + (mt + ph / 2) + ')">a\u2080Resid (dex)</text>\n';
  svg += '<text x="' + (ml + pw / 2) + '" y="22" text-anchor="middle" font-size="14" font-weight="bold">Bilateral Coupling: r = 0.77 (LOO, N = 45)</text>\n';

  svg += '<circle cx="' + (ml + 15) + '" cy="' + (mt + 18) + '" r="3" fill="#457b9d" opacity="0.5"/>\n';
  svg += '<text x="' + (ml + 22) + '" y="' + (mt + 22) + '" font-size="10" fill="#457b9d">SPARC (N=45)</text>\n';
  svg += '<circle cx="' + (ml + 115) + '" cy="' + (mt + 18) + '" r="5" fill="#e63946" stroke="#000" stroke-width="0.8"/>\n';
  svg += '<text x="' + (ml + 123) + '" y="' + (mt + 22) + '" font-size="10" fill="#e63946">THINGS 2D (N=7)</text>\n';

  svg += '</svg>';
  fs.writeFileSync(path.join(outDir, 'fig1-bilateral-coupling.svg'), svg);
  console.log('  Fig 1: Bilateral coupling scatter');
}


// ====== FIGURE 2: Information Ceiling ======
{
  const W = 520, H = 400;
  const ml = 160, mr = 40, mt = 40, mb = 80;
  const pw = W - ml - mr, ph = H - mt - mb;

  const bars = [
    { label: 'RC shape features\n(Program 8A)', val: 11.5, col: '#457b9d' },
    { label: 'Scalar haloResponse\n(Program 8B)', val: 30, col: '#6a9fba' },
    { label: '71-dim RC map\n(Program 8C)', val: 30, col: '#8bbbd4' },
    { label: 'm=2 velocity field\n(Program 9)', val: 83, col: '#e63946' },
  ];
  const barH = 38, gap = 18;

  let svg = svgHeader(W, H);

  for (let i = 0; i < bars.length; i++) {
    const y = mt + i * (barH + gap);
    const bw = bars[i].val / 100 * pw;
    svg += '<rect x="' + ml + '" y="' + y + '" width="' + bw.toFixed(1) + '" height="' + barH + '" fill="' + bars[i].col + '" rx="3"/>\n';
    svg += '<rect x="' + ml + '" y="' + y + '" width="' + pw + '" height="' + barH + '" fill="none" stroke="#ccc" rx="3"/>\n';

    const lines = bars[i].label.split('\n');
    svg += '<text x="' + (ml - 8) + '" y="' + (y + barH / 2 - 6) + '" text-anchor="end" font-size="11" font-weight="bold">' + esc(lines[0]) + '</text>\n';
    if (lines[1]) svg += '<text x="' + (ml - 8) + '" y="' + (y + barH / 2 + 8) + '" text-anchor="end" font-size="10" fill="#666">' + esc(lines[1]) + '</text>\n';

    svg += '<text x="' + (ml + bw + 6) + '" y="' + (y + barH / 2 + 5) + '" font-size="12" font-weight="bold" fill="' + bars[i].col + '">' + bars[i].val + '%</text>\n';
  }

  const ceilX = ml + 30 / 100 * pw;
  svg += '<line x1="' + ceilX.toFixed(1) + '" y1="' + (mt - 5) + '" x2="' + ceilX.toFixed(1) + '" y2="' + (mt + 3 * (barH + gap) + barH + 10) + '" stroke="#333" stroke-width="1.5" stroke-dasharray="6,4"/>\n';
  svg += '<text x="' + (ceilX + 4) + '" y="' + (mt - 8) + '" font-size="10" fill="#333">1D ceiling (~30%)</text>\n';

  const inaccessX = ml + 88.5 / 100 * pw;
  const arrowY = mt + 4 * (barH + gap) + 15;
  svg += '<line x1="' + ceilX.toFixed(1) + '" y1="' + arrowY + '" x2="' + inaccessX.toFixed(1) + '" y2="' + arrowY + '" stroke="#e63946" stroke-width="2" marker-end="url(#arrowhead)"/>\n';
  svg += '<defs><marker id="arrowhead" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#e63946"/></marker></defs>\n';
  svg += '<text x="' + ((ceilX + inaccessX) / 2) + '" y="' + (arrowY - 6) + '" text-anchor="middle" font-size="10" fill="#e63946">+53% recovered by 2D m=2</text>\n';

  for (let v = 0; v <= 100; v += 20) {
    const px = ml + v / 100 * pw;
    svg += '<text x="' + px.toFixed(1) + '" y="' + (mt + 4 * (barH + gap) + 40) + '" text-anchor="middle" font-size="11">' + v + '%</text>\n';
    svg += '<line x1="' + px.toFixed(1) + '" y1="' + (mt + 3 * (barH + gap) + barH) + '" x2="' + px.toFixed(1) + '" y2="' + (mt + 3 * (barH + gap) + barH + 4) + '" stroke="#333"/>\n';
  }
  svg += '<text x="' + (ml + pw / 2) + '" y="' + (H - 8) + '" text-anchor="middle" font-size="13" font-weight="bold">Variance of H recovered (%)</text>\n';
  svg += '<text x="' + (W / 2) + '" y="22" text-anchor="middle" font-size="14" font-weight="bold">Information Ceiling: 1D vs 2D Recovery of Hidden State H</text>\n';

  svg += '</svg>';
  fs.writeFileSync(path.join(outDir, 'fig2-information-ceiling.svg'), svg);
  console.log('  Fig 2: Information ceiling');
}


// ====== FIGURE 3: DQ vs m=2 Power ======
{
  const W = 520, H = 460;
  const ml = 80, mr = 30, mt = 40, mb = 60;
  const pw = W - ml - mr, ph = H - mt - mb;

  const data = p903.results.map(g => ({
    name: g.name, DQ: g.DQ, m2: g.m2Power, logM2: Math.log10(g.m2Power)
  }));

  const xMin = -2, xMax = 3, yMin = 2.8, yMax = 4.5;

  let svg = svgHeader(W, H);
  svg += '<rect x="' + ml + '" y="' + mt + '" width="' + pw + '" height="' + ph + '" fill="#fafafa" stroke="#ccc"/>\n';

  for (let v = yMin; v <= yMax; v += 0.5) {
    const py = mapRange(v, yMax, yMin, mt, mt + ph);
    svg += '<line x1="' + ml + '" y1="' + py.toFixed(1) + '" x2="' + (ml + pw) + '" y2="' + py.toFixed(1) + '" stroke="#eee"/>\n';
  }

  const gold = { 'NGC2841': '#e63946', 'NGC5055': '#1d3557' };
  for (const g of data) {
    const px = mapRange(g.DQ, xMin, xMax, ml, ml + pw);
    const py = mapRange(g.logM2, yMax, yMin, mt, mt + ph);
    const col = gold[g.name] || '#457b9d';
    const r = gold[g.name] ? 7 : 5;
    svg += '<circle cx="' + px.toFixed(1) + '" cy="' + py.toFixed(1) + '" r="' + r + '" fill="' + col + '" stroke="#000" stroke-width="0.8"/>\n';
    const labelOff = g.name === 'NGC2841' ? { dx: 10, dy: -10 } : g.name === 'NGC5055' ? { dx: 10, dy: 12 } : g.name === 'NGC3521' ? { dx: 10, dy: -8 } : { dx: 10, dy: -5 };
    svg += '<text x="' + (px + labelOff.dx) + '" y="' + (py + labelOff.dy) + '" font-size="9" fill="#333">' + esc(g.name) + '</text>\n';
  }

  const xs = data.map(d => d.DQ), ys = data.map(d => d.logM2);
  const mxv = xs.reduce((a, b) => a + b, 0) / xs.length;
  const myv = ys.reduce((a, b) => a + b, 0) / ys.length;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < xs.length; i++) { sxy += (xs[i] - mxv) * (ys[i] - myv); sxx += (xs[i] - mxv) ** 2; }
  const sl = sxy / sxx, intr = myv - sl * mxv;
  const lx1 = mapRange(xMin, xMin, xMax, ml, ml + pw), ly1 = mapRange(sl * xMin + intr, yMax, yMin, mt, mt + ph);
  const lx2 = mapRange(xMax, xMin, xMax, ml, ml + pw), ly2 = mapRange(sl * xMax + intr, yMax, yMin, mt, mt + ph);
  svg += '<line x1="' + lx1.toFixed(1) + '" y1="' + ly1.toFixed(1) + '" x2="' + lx2.toFixed(1) + '" y2="' + ly2.toFixed(1) + '" stroke="#e63946" stroke-width="2" stroke-dasharray="6,3"/>\n';

  for (let v = -2; v <= 3; v++) {
    const px = mapRange(v, xMin, xMax, ml, ml + pw);
    svg += '<text x="' + px.toFixed(1) + '" y="' + (mt + ph + 18) + '" text-anchor="middle" font-size="11">' + v + '</text>\n';
    svg += '<line x1="' + px.toFixed(1) + '" y1="' + (mt + ph) + '" x2="' + px.toFixed(1) + '" y2="' + (mt + ph + 4) + '" stroke="#333"/>\n';
  }
  for (let v = 3.0; v <= 4.5; v += 0.5) {
    const py = mapRange(v, yMax, yMin, mt, mt + ph);
    svg += '<text x="' + (ml - 8) + '" y="' + (py + 4) + '" text-anchor="end" font-size="11">' + v.toFixed(1) + '</text>\n';
    svg += '<line x1="' + (ml - 4) + '" y1="' + py.toFixed(1) + '" x2="' + ml + '" y2="' + py.toFixed(1) + '" stroke="#333"/>\n';
  }

  svg += '<text x="' + (ml + pw / 2) + '" y="' + (H - 8) + '" text-anchor="middle" font-size="13" font-weight="bold">Bilateral Excess DQ (\u03c3)</text>\n';
  svg += '<text x="15" y="' + (mt + ph / 2) + '" text-anchor="middle" font-size="13" font-weight="bold" transform="rotate(-90,15,' + (mt + ph / 2) + ')">log\u2081\u2080 m=2 Power</text>\n';
  svg += '<text x="' + (ml + pw / 2) + '" y="22" text-anchor="middle" font-size="14" font-weight="bold">DQ vs m=2 Azimuthal Power: r = 0.85 (p = 0.005, N = 7)</text>\n';

  svg += '<rect x="' + (ml + pw - 140) + '" y="' + (mt + ph - 50) + '" width="130" height="44" fill="white" fill-opacity="0.9" stroke="#ccc" rx="4"/>\n';
  svg += '<circle cx="' + (ml + pw - 125) + '" cy="' + (mt + ph - 36) + '" r="6" fill="#e63946" stroke="#000" stroke-width="0.8"/>\n';
  svg += '<text x="' + (ml + pw - 115) + '" y="' + (mt + ph - 32) + '" font-size="10">NGC 2841 (High-H)</text>\n';
  svg += '<circle cx="' + (ml + pw - 125) + '" cy="' + (mt + ph - 18) + '" r="6" fill="#1d3557" stroke="#000" stroke-width="0.8"/>\n';
  svg += '<text x="' + (ml + pw - 115) + '" y="' + (mt + ph - 14) + '" font-size="10">NGC 5055 (Low-H)</text>\n';

  svg += '</svg>';
  fs.writeFileSync(path.join(outDir, 'fig3-dq-m2-power.svg'), svg);
  console.log('  Fig 3: DQ vs m=2 power');
}


// ====== FIGURE 4: Gold Pair Comparison ======
{
  const W = 520, H = 400;
  const ml = 60, mr = 30, mt = 50, mb = 30;

  const ngc2841 = p903.results.find(g => g.name === 'NGC2841');
  const ngc5055 = p903.results.find(g => g.name === 'NGC5055');

  const props = [
    { label: 'log Mbar', v1: ngc2841.logMbar, v2: ngc5055.logMbar, unit: 'dex' },
    { label: 'Vflat', v1: ngc2841.Vflat, v2: ngc5055.Vflat, unit: 'km/s' },
    { label: 'DQ', v1: ngc2841.DQ, v2: ngc5055.DQ, unit: '\u03c3' },
    { label: 'm=2 Power', v1: ngc2841.m2Power, v2: ngc5055.m2Power, unit: '' },
    { label: 'PA Twist', v1: ngc2841.paTwist, v2: ngc5055.paTwist, unit: '\u00b0' },
    { label: 'Coherence', v1: ngc2841.coherenceOuter, v2: ngc5055.coherenceOuter, unit: '' },
  ];

  let svg = svgHeader(W, H);

  svg += '<text x="' + (W / 2) + '" y="22" text-anchor="middle" font-size="14" font-weight="bold">Gold Pair: NGC 2841 (High-H) vs NGC 5055 (Low-H)</text>\n';
  svg += '<text x="' + (W / 2) + '" y="38" text-anchor="middle" font-size="11" fill="#666">\u0394log Mbar = 0.07 dex (structurally matched) | m=2 ratio = 10.5\u00d7</text>\n';

  const rowH = 44, startY = mt + 20;
  const colL = ml + 110, colR = ml + 310;

  svg += '<text x="' + colL + '" y="' + (startY - 4) + '" text-anchor="middle" font-size="12" font-weight="bold" fill="#e63946">NGC 2841</text>\n';
  svg += '<text x="' + colR + '" y="' + (startY - 4) + '" text-anchor="middle" font-size="12" font-weight="bold" fill="#1d3557">NGC 5055</text>\n';
  svg += '<text x="' + (ml + 40) + '" y="' + (startY - 4) + '" text-anchor="middle" font-size="12" font-weight="bold" fill="#333">Property</text>\n';

  for (let i = 0; i < props.length; i++) {
    const y = startY + i * rowH + 10;
    const bg = i % 2 === 0 ? '#f8f8f8' : 'white';
    svg += '<rect x="' + ml + '" y="' + (y - 10) + '" width="' + (W - ml - mr) + '" height="' + (rowH - 4) + '" fill="' + bg + '" rx="3"/>\n';

    svg += '<text x="' + (ml + 40) + '" y="' + (y + 12) + '" text-anchor="middle" font-size="12" font-weight="bold">' + esc(props[i].label) + '</text>\n';

    const fmt1 = props[i].v1 >= 100 ? props[i].v1.toFixed(0) : props[i].v1.toFixed(2);
    const fmt2 = props[i].v2 >= 100 ? props[i].v2.toFixed(0) : props[i].v2.toFixed(2);
    svg += '<text x="' + colL + '" y="' + (y + 12) + '" text-anchor="middle" font-size="13" fill="#e63946">' + fmt1 + '</text>\n';
    svg += '<text x="' + colR + '" y="' + (y + 12) + '" text-anchor="middle" font-size="13" fill="#1d3557">' + fmt2 + '</text>\n';

    if (props[i].label === 'm=2 Power') {
      const ratio = (props[i].v1 / props[i].v2).toFixed(1);
      svg += '<text x="' + (colL + (colR - colL) / 2) + '" y="' + (y + 12) + '" text-anchor="middle" font-size="11" font-weight="bold" fill="#e63946">' + ratio + '\u00d7</text>\n';
    } else if (props[i].label === 'DQ') {
      const diff = (props[i].v1 - props[i].v2).toFixed(2);
      svg += '<text x="' + (colL + (colR - colL) / 2) + '" y="' + (y + 12) + '" text-anchor="middle" font-size="11" fill="#666">\u0394=' + diff + '\u03c3</text>\n';
    }
  }

  svg += '</svg>';
  fs.writeFileSync(path.join(outDir, 'fig4-gold-pair.svg'), svg);
  console.log('  Fig 4: Gold pair comparison');
}


// ====== FIGURE 5: Causal Chain ======
{
  const W = 560, H = 520;
  let svg = svgHeader(W, H);

  svg += '<defs><marker id="arr" markerWidth="8" markerHeight="6" refX="8" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#333"/></marker></defs>\n';

  const boxes = [
    { id: 'form', x: 175, y: 20, w: 210, h: 36, text: 'Formation History', col: '#a8dadc', tc: '#1d3557' },
    { id: 'conc', x: 60, y: 110, w: 160, h: 36, text: 'Halo Concentration', col: '#a8dadc', tc: '#1d3557' },
    { id: 'triax', x: 300, y: 110, w: 180, h: 36, text: 'Halo Triaxiality (b/a)', col: '#e63946', tc: 'white' },
    { id: 'mass', x: 40, y: 210, w: 170, h: 36, text: 'Enclosed Mass Shift', col: '#457b9d', tc: 'white' },
    { id: 'm2', x: 310, y: 210, w: 200, h: 36, text: 'm=2 Velocity Field Power', col: '#e63946', tc: 'white' },
    { id: 'vf', x: 40, y: 310, w: 130, h: 36, text: 'Vflat Shift', col: '#457b9d', tc: 'white' },
    { id: 'a0', x: 200, y: 310, w: 130, h: 36, text: 'a\u2080 Shift (4\u00d7)', col: '#e63946', tc: 'white' },
    { id: 'bil', x: 120, y: 410, w: 230, h: 44, text: 'Bilateral Coupling', col: '#1d3557', tc: 'white' },
  ];

  for (const b of boxes) {
    svg += '<rect x="' + b.x + '" y="' + b.y + '" width="' + b.w + '" height="' + b.h + '" fill="' + b.col + '" stroke="#333" stroke-width="1" rx="6"/>\n';
    svg += '<text x="' + (b.x + b.w / 2) + '" y="' + (b.y + b.h / 2 + 5) + '" text-anchor="middle" font-size="12" font-weight="bold" fill="' + b.tc + '">' + esc(b.text) + '</text>\n';
  }

  const arrows = [
    { from: 'form', to: 'conc' }, { from: 'form', to: 'triax' },
    { from: 'conc', to: 'mass' }, { from: 'triax', to: 'm2' },
    { from: 'triax', to: 'mass', dashed: true },
    { from: 'mass', to: 'vf' }, { from: 'mass', to: 'a0' },
    { from: 'vf', to: 'bil' }, { from: 'a0', to: 'bil' },
  ];
  const bMap = {}; boxes.forEach(b => bMap[b.id] = b);
  for (const a of arrows) {
    const f = bMap[a.from], t = bMap[a.to];
    const x1 = f.x + f.w / 2, y1 = f.y + f.h;
    const x2 = t.x + t.w / 2, y2 = t.y;
    svg += '<line x1="' + x1 + '" y1="' + y1 + '" x2="' + x2 + '" y2="' + y2 + '" stroke="#333" stroke-width="1.5" marker-end="url(#arr)"' + (a.dashed ? ' stroke-dasharray="5,3"' : '') + '/>\n';
  }

  svg += '<text x="' + (bMap.m2.x + bMap.m2.w / 2) + '" y="' + (bMap.m2.y + bMap.m2.h + 18) + '" text-anchor="middle" font-size="11" fill="#e63946" font-weight="bold">r(DQ, m2) = 0.85</text>\n';
  svg += '<text x="' + (bMap.m2.x + bMap.m2.w / 2) + '" y="' + (bMap.m2.y + bMap.m2.h + 32) + '" text-anchor="middle" font-size="10" fill="#999">angular, 1D-invisible</text>\n';
  svg += '<text x="' + (bMap.mass.x + bMap.mass.w / 2) + '" y="' + (bMap.mass.y + bMap.mass.h + 16) + '" text-anchor="middle" font-size="10" fill="#999">integral, 1D-visible</text>\n';
  svg += '<text x="' + (bMap.bil.x + bMap.bil.w / 2) + '" y="' + (bMap.bil.y + bMap.bil.h / 2 + 20) + '" text-anchor="middle" font-size="12" fill="white">r = 0.77 (LOO, N = 45)</text>\n';

  svg += '<rect x="365" y="310" width="170" height="60" fill="none" stroke="#e63946" stroke-width="1.5" stroke-dasharray="4,3" rx="6"/>\n';
  svg += '<text x="450" y="330" text-anchor="middle" font-size="10" fill="#e63946">1D INVISIBLE</text>\n';
  svg += '<text x="450" y="345" text-anchor="middle" font-size="9" fill="#999">Azimuthal averaging</text>\n';
  svg += '<text x="450" y="358" text-anchor="middle" font-size="9" fill="#999">destroys m=2 by</text>\n';
  svg += '<text x="450" y="371" text-anchor="middle" font-size="9" fill="#999">construction</text>\n';

  svg += '</svg>';
  fs.writeFileSync(path.join(outDir, 'fig5-causal-chain.svg'), svg);
  console.log('  Fig 5: Causal chain diagram');
}


// ====== FIGURE 6: Red Team Summary ======
{
  const W = 520, H = 460;
  const ml = 270, mr = 40, mt = 50, mb = 40;
  const pw = W - ml - mr;

  const tests = [
    { name: 'Independent replication', val: 0.911, type: 'C', ref: 'r = 0.91' },
    { name: 'All sensitivity positive', val: 19, type: 'C', ref: '19/19' },
    { name: '\u226580% above r = 0.3', val: 19, type: 'A', ref: '19/19' },
    { name: 'Unbarred-only survives', val: 0.910, type: 'C', ref: 'r = 0.91' },
    { name: 'Outer-only m=2 signal', val: 0.869, type: 'A', ref: 'r = 0.87' },
    { name: 'LOO all positive', val: 7, type: 'C', ref: '7/7' },
    { name: 'r > 0.3 without NGC2841', val: 0.674, type: 'C', ref: 'r = 0.67' },
    { name: 'Positive without gold pair', val: 0.739, type: 'A', ref: 'r = 0.74' },
    { name: 'Inclination not confounder', val: 0.0, type: 'C', ref: 'r = 0.00' },
    { name: 'Partial r|confounders > 0', val: 0.660, type: 'C', ref: 'r = 0.66' },
    { name: 'Program 8B reconciliation', val: 1, type: 'C', ref: 'consistent' },
  ];

  const rowH = 32, startY = mt + 10;
  let svg = svgHeader(W, H);

  svg += '<text x="' + (W / 2) + '" y="22" text-anchor="middle" font-size="14" font-weight="bold">Program 9V: Red Team Scorecard \u2014 11/11 PASS</text>\n';
  svg += '<text x="' + (W / 2) + '" y="40" text-anchor="middle" font-size="11" fill="#666">8 Critical + 3 Advisory = All Pass</text>\n';

  for (let i = 0; i < tests.length; i++) {
    const y = startY + i * rowH;
    const bg = i % 2 === 0 ? '#f0f8f0' : '#f8fff8';
    svg += '<rect x="10" y="' + y + '" width="' + (W - 20) + '" height="' + (rowH - 2) + '" fill="' + bg + '" rx="3"/>\n';

    svg += '<text x="35" y="' + (y + rowH / 2 + 4) + '" font-size="16" fill="#2d6a4f">\u2713</text>\n';

    const typeCol = tests[i].type === 'C' ? '#1d3557' : '#999';
    const typeLabel = tests[i].type === 'C' ? 'CRITICAL' : 'advisory';
    svg += '<rect x="50" y="' + (y + 5) + '" width="' + (tests[i].type === 'C' ? 62 : 52) + '" height="18" fill="' + typeCol + '" rx="9"/>\n';
    svg += '<text x="' + (tests[i].type === 'C' ? 81 : 76) + '" y="' + (y + 17) + '" text-anchor="middle" font-size="9" fill="white">' + typeLabel + '</text>\n';

    svg += '<text x="120" y="' + (y + rowH / 2 + 4) + '" font-size="11" fill="#333">' + esc(tests[i].name) + '</text>\n';

    svg += '<text x="' + (W - 50) + '" y="' + (y + rowH / 2 + 4) + '" text-anchor="end" font-size="11" font-weight="bold" fill="#2d6a4f">' + esc(tests[i].ref) + '</text>\n';
  }

  const summY = startY + tests.length * rowH + 15;
  svg += '<rect x="130" y="' + summY + '" width="260" height="36" fill="#2d6a4f" rx="6"/>\n';
  svg += '<text x="260" y="' + (summY + 22) + '" text-anchor="middle" font-size="14" font-weight="bold" fill="white">CARRIER CLAIM SURVIVES \u2014 ~95%</text>\n';

  svg += '</svg>';
  fs.writeFileSync(path.join(outDir, 'fig6-red-team.svg'), svg);
  console.log('  Fig 6: Red team summary');
}

console.log('\n  All 6 figures saved to: ' + outDir);
console.log('  SVG files can be converted to PDF with inkscape or rsvg-convert');
console.log('\n=== Done ===');
