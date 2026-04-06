const fs = require('fs');
const path = require('path');

const outDir = path.join(__dirname, '..', 'paper', 'figures');
if (!fs.existsSync(outDir)) fs.mkdirSync(outDir, { recursive: true });

const phase302 = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase302-regime-law.json'), 'utf8'));
const phase303 = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase303-physical-interpretation.json'), 'utf8'));
const phase301 = JSON.parse(fs.readFileSync(path.join(__dirname, '..', 'public', 'phase301-vfresid-drivers.json'), 'utf8'));

const W = 500, H = 380;
const M = { top: 50, right: 30, bottom: 60, left: 70 };
const PW = W - M.left - M.right;
const PH = H - M.top - M.bottom;

function svgHeader(width, height, title) {
  return `<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" width="${width}" height="${height}" font-family="serif">
<rect width="${width}" height="${height}" fill="white"/>
<text x="${width/2}" y="25" text-anchor="middle" font-size="13" font-weight="bold">${title}</text>
`;
}

function svgFooter() { return '</svg>\n'; }

function escXml(s) { return String(s).replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;'); }

function makeAxis(xLabel, yLabel, xTicks, yTicks, xMin, xMax, yMin, yMax) {
  let s = '';
  s += `<line x1="${M.left}" y1="${H-M.bottom}" x2="${M.left+PW}" y2="${H-M.bottom}" stroke="black" stroke-width="1"/>`;
  s += `<line x1="${M.left}" y1="${M.top}" x2="${M.left}" y2="${H-M.bottom}" stroke="black" stroke-width="1"/>`;
  s += `<text x="${M.left+PW/2}" y="${H-10}" text-anchor="middle" font-size="12">${escXml(xLabel)}</text>`;
  s += `<text x="15" y="${M.top+PH/2}" text-anchor="middle" font-size="12" transform="rotate(-90,15,${M.top+PH/2})">${escXml(yLabel)}</text>`;
  for (const t of xTicks) {
    const x = M.left + (t - xMin) / (xMax - xMin) * PW;
    s += `<line x1="${x}" y1="${H-M.bottom}" x2="${x}" y2="${H-M.bottom+5}" stroke="black"/>`;
    s += `<text x="${x}" y="${H-M.bottom+18}" text-anchor="middle" font-size="10">${t}</text>`;
  }
  for (const t of yTicks) {
    const y = H - M.bottom - (t - yMin) / (yMax - yMin) * PH;
    s += `<line x1="${M.left-5}" y1="${y}" x2="${M.left}" y2="${y}" stroke="black"/>`;
    s += `<text x="${M.left-8}" y="${y+4}" text-anchor="end" font-size="10">${t}</text>`;
  }
  return s;
}

function px(val, min, max) { return M.left + (val - min) / (max - min) * PW; }
function py(val, min, max) { return H - M.bottom - (val - min) / (max - min) * PH; }

function fig1_hierarchy() {
  const data = [
    { label: 'Core\n(3-axis)', gap: 44.1, color: '#2563eb' },
    { label: 'Core +\nVfResid', gap: 61.1, color: '#7c3aed' },
    { label: 'Full\n(5-axis)', gap: 65.4, color: '#059669' }
  ];
  const w = 500, h = 350;
  const ml = 80, mr = 30, mt = 50, mb = 70;
  const pw = w - ml - mr, ph = h - mt - mb;
  let s = svgHeader(w, h, 'Model Hierarchy: LOO Cross-Validation Performance');

  const barW = pw / data.length * 0.6;
  const gap = pw / data.length;
  for (let i = 0; i < data.length; i++) {
    const cx = ml + gap * (i + 0.5);
    const bh = (data[i].gap / 80) * ph;
    const by = mt + ph - bh;
    s += `<rect x="${cx - barW/2}" y="${by}" width="${barW}" height="${bh}" fill="${data[i].color}" opacity="0.8" rx="3"/>`;
    s += `<text x="${cx}" y="${by - 8}" text-anchor="middle" font-size="12" font-weight="bold">${data[i].gap}%</text>`;
    const lines = data[i].label.split('\n');
    lines.forEach((line, j) => {
      s += `<text x="${cx}" y="${mt + ph + 18 + j * 14}" text-anchor="middle" font-size="10">${escXml(line)}</text>`;
    });
  }

  s += `<line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="black" stroke-width="1"/>`;
  s += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="black" stroke-width="1"/>`;
  for (let t = 0; t <= 80; t += 20) {
    const y = mt + ph - (t / 80) * ph;
    s += `<line x1="${ml-5}" y1="${y}" x2="${ml}" y2="${y}" stroke="black"/>`;
    s += `<text x="${ml-8}" y="${y+4}" text-anchor="end" font-size="10">${t}%</text>`;
    if (t > 0) s += `<line x1="${ml}" y1="${y}" x2="${ml+pw}" y2="${y}" stroke="#ddd" stroke-width="0.5" stroke-dasharray="4,3"/>`;
  }
  s += `<text x="18" y="${mt+ph/2}" text-anchor="middle" font-size="11" transform="rotate(-90,18,${mt+ph/2})">LOO gap (% variance explained)</text>`;

  s += `<text x="${ml + pw/2}" y="${h - 5}" text-anchor="middle" font-size="10" fill="#666">N = 45 published-quality SPARC galaxies</text>`;

  const annotations = [
    { from: 0, to: 1, label: '+17.1 pp', y: mt + ph - (61.1/80)*ph - 30 },
    { from: 1, to: 2, label: '+4.3 pp', y: mt + ph - (65.4/80)*ph - 30 }
  ];
  for (const a of annotations) {
    const x1 = ml + gap * (a.from + 0.5);
    const x2 = ml + gap * (a.to + 0.5);
    const xm = (x1 + x2) / 2;
    s += `<line x1="${x1}" y1="${a.y+12}" x2="${x2}" y2="${a.y+12}" stroke="#333" stroke-width="1" marker-start="url(#arrowL)" marker-end="url(#arrowR)"/>`;
    s += `<text x="${xm}" y="${a.y+6}" text-anchor="middle" font-size="9" fill="#333">${a.label}</text>`;
  }
  s += `<defs>
    <marker id="arrowR" markerWidth="6" markerHeight="6" refX="5" refY="3" orient="auto"><path d="M0,0 L6,3 L0,6" fill="none" stroke="#333" stroke-width="0.8"/></marker>
    <marker id="arrowL" markerWidth="6" markerHeight="6" refX="1" refY="3" orient="auto"><path d="M6,0 L0,3 L6,6" fill="none" stroke="#333" stroke-width="0.8"/></marker>
  </defs>`;

  s += svgFooter();
  return s;
}

function fig2_external() {
  const data = [
    { label: 'Full\n(N=59)', coreGap: -8.1, fullGap: 8.2 },
    { label: 'High-V\n(N=16)', coreGap: 16.7, fullGap: 34.3 },
    { label: 'V-High-V\n(N=8)', coreGap: -16.7, fullGap: 48.7 },
    { label: 'Q1+HV\n(N=11)', coreGap: 25.4, fullGap: 59.0 }
  ];
  const w = 520, h = 380;
  const ml = 80, mr = 30, mt = 50, mb = 85;
  const pw = w - ml - mr, ph = h - mt - mb;
  let s = svgHeader(w, h, 'External Validation: Frozen-Coefficient Transfer');

  const yMin = -30, yMax = 70;
  s += `<line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="black" stroke-width="1"/>`;
  s += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="black" stroke-width="1"/>`;
  for (let t = -20; t <= 60; t += 20) {
    const yp = mt + ph - ((t - yMin) / (yMax - yMin)) * ph;
    s += `<line x1="${ml-5}" y1="${yp}" x2="${ml}" y2="${yp}" stroke="black"/>`;
    s += `<text x="${ml-8}" y="${yp+4}" text-anchor="end" font-size="10">${t}%</text>`;
    s += `<line x1="${ml}" y1="${yp}" x2="${ml+pw}" y2="${yp}" stroke="#eee" stroke-width="0.5" stroke-dasharray="4,3"/>`;
  }
  const zeroY = mt + ph - ((0 - yMin) / (yMax - yMin)) * ph;
  s += `<line x1="${ml}" y1="${zeroY}" x2="${ml+pw}" y2="${zeroY}" stroke="#999" stroke-width="1" stroke-dasharray="6,3"/>`;

  const groupW = pw / data.length;
  const barW = groupW * 0.3;
  for (let i = 0; i < data.length; i++) {
    const cx = ml + groupW * (i + 0.5);
    const coreH = ((data[i].coreGap - yMin) / (yMax - yMin)) * ph;
    const fullH = ((data[i].fullGap - yMin) / (yMax - yMin)) * ph;
    const coreY = mt + ph - coreH;
    const fullY = mt + ph - fullH;

    s += `<rect x="${cx - barW - 2}" y="${Math.min(coreY, zeroY)}" width="${barW}" height="${Math.abs(coreY - zeroY)}" fill="#f59e0b" opacity="0.8" rx="2"/>`;
    s += `<rect x="${cx + 2}" y="${Math.min(fullY, zeroY)}" width="${barW}" height="${Math.abs(fullY - zeroY)}" fill="#2563eb" opacity="0.8" rx="2"/>`;

    s += `<text x="${cx - barW/2 - 2}" y="${coreY + (data[i].coreGap >= 0 ? -5 : 14)}" text-anchor="middle" font-size="9">${data[i].coreGap}%</text>`;
    s += `<text x="${cx + barW/2 + 2}" y="${fullY - 5}" text-anchor="middle" font-size="9" font-weight="bold">${data[i].fullGap}%</text>`;

    const lines = data[i].label.split('\n');
    lines.forEach((l, j) => {
      s += `<text x="${cx}" y="${mt + ph + 18 + j * 13}" text-anchor="middle" font-size="10">${escXml(l)}</text>`;
    });
  }

  s += `<text x="18" y="${mt+ph/2}" text-anchor="middle" font-size="11" transform="rotate(-90,18,${mt+ph/2})">Transfer gap (%)</text>`;

  const ly = h - 10;
  s += `<rect x="${ml + 10}" y="${ly - 9}" width="12" height="10" fill="#f59e0b" opacity="0.8" rx="2"/>`;
  s += `<text x="${ml + 26}" y="${ly}" font-size="10">Core (improved Mhost)</text>`;
  s += `<rect x="${ml + 170}" y="${ly - 9}" width="12" height="10" fill="#2563eb" opacity="0.8" rx="2"/>`;
  s += `<text x="${ml + 186}" y="${ly}" font-size="10">Core + VfResid</text>`;

  s += svgFooter();
  return s;
}

function fig3_regime() {
  const cumData = phase302.test1_activationShape.internal.cumulative;
  let extCum = [];
  if (phase302.test1_activationShape.external) {
    extCum = phase302.test1_activationShape.external.cumulative || [];
  }

  const w = 520, h = 380;
  let s = svgHeader(w, h, 'Regime Dependence: VfResid-a\u2080 Coupling vs Vflat Threshold');

  const xMin = 40, xMax = 210, yMin = 0.3, yMax = 1.0;
  s += makeAxis('Vflat threshold (km/s)', 'r(VfResid, log a\u2080)', [50, 80, 100, 120, 150, 180, 200], [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], xMin, xMax, yMin, yMax);

  let points = cumData.map(d => ({ x: px(d.threshold, xMin, xMax), y: py(d.r, yMin, yMax) }));
  if (points.length > 1) {
    let pathD = `M${points[0].x},${points[0].y}`;
    for (let i = 1; i < points.length; i++) pathD += ` L${points[i].x},${points[i].y}`;
    s += `<path d="${pathD}" fill="none" stroke="#2563eb" stroke-width="2"/>`;
  }
  for (const p of points) {
    s += `<circle cx="${p.x}" cy="${p.y}" r="4" fill="#2563eb"/>`;
  }

  if (extCum.length > 0) {
    let ePts = extCum.map(d => ({ x: px(d.threshold, xMin, xMax), y: py(d.r, yMin, yMax) }));
    if (ePts.length > 1) {
      let pathD = `M${ePts[0].x},${ePts[0].y}`;
      for (let i = 1; i < ePts.length; i++) pathD += ` L${ePts[i].x},${ePts[i].y}`;
      s += `<path d="${pathD}" fill="none" stroke="#dc2626" stroke-width="2" stroke-dasharray="6,3"/>`;
    }
    for (const p of ePts) {
      s += `<circle cx="${p.x}" cy="${p.y}" r="4" fill="#dc2626"/>`;
    }
  }

  const threshX = px(120, xMin, xMax);
  s += `<line x1="${threshX}" y1="${M.top}" x2="${threshX}" y2="${H-M.bottom}" stroke="#999" stroke-width="1" stroke-dasharray="4,3"/>`;
  s += `<text x="${threshX+5}" y="${M.top+15}" font-size="9" fill="#666">120 km/s</text>`;

  s += `<circle cx="${M.left+PW-120}" cy="${M.top+15}" r="4" fill="#2563eb"/>`;
  s += `<text x="${M.left+PW-112}" y="${M.top+19}" font-size="10">Internal (N=45)</text>`;
  if (extCum.length > 0) {
    s += `<circle cx="${M.left+PW-120}" cy="${M.top+32}" r="4" fill="#dc2626"/>`;
    s += `<text x="${M.left+PW-112}" y="${M.top+36}" font-size="10">External (N=59)</text>`;
  }

  s += svgFooter();
  return s;
}

function fig4_vfresid_scatter() {
  const w = 500, h = 400;
  let s = svgHeader(w, h, 'VfResid as Baryon-Halo Coupling Proxy');

  const drivers = phase301.singlePredictorRankings;
  const intData = drivers.internal;

  const barData = intData.slice(0, 8).map(d => ({
    name: d.name === 'mondImprove' ? 'MOND Improv.' :
          d.name === 'logMeanRun' ? 'log MR' :
          d.name === 'envCode' ? 'Environ.' :
          d.name === 'lhOuter' ? 'lhOuter' :
          d.name === 'logMhost' ? 'log Mhost' :
          d.name === 'dhlMSE' ? 'DHL MSE' :
          d.name,
    r: Math.abs(d.r),
    sign: d.r >= 0 ? '+' : '-',
    color: d.r >= 0 ? '#2563eb' : '#dc2626'
  }));

  const ml = 110, mr = 30, mt = 50, mb = 50;
  const pw = w - ml - mr, ph = h - mt - mb;
  const rowH = ph / barData.length;

  s += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="black" stroke-width="1"/>`;
  s += `<line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="black" stroke-width="1"/>`;
  for (let t = 0; t <= 0.6; t += 0.2) {
    const x = ml + (t / 0.7) * pw;
    s += `<line x1="${x}" y1="${mt+ph}" x2="${x}" y2="${mt+ph+5}" stroke="black"/>`;
    s += `<text x="${x}" y="${mt+ph+18}" text-anchor="middle" font-size="10">${t.toFixed(1)}</text>`;
    s += `<line x1="${x}" y1="${mt}" x2="${x}" y2="${mt+ph}" stroke="#eee" stroke-width="0.5" stroke-dasharray="4,3"/>`;
  }
  s += `<text x="${ml+pw/2}" y="${h-5}" text-anchor="middle" font-size="11">|r| with VfResid</text>`;

  for (let i = 0; i < barData.length; i++) {
    const y = mt + rowH * (i + 0.5);
    const bw = (barData[i].r / 0.7) * pw;
    s += `<rect x="${ml}" y="${y - rowH * 0.35}" width="${bw}" height="${rowH * 0.7}" fill="${barData[i].color}" opacity="0.7" rx="2"/>`;
    s += `<text x="${ml - 5}" y="${y + 4}" text-anchor="end" font-size="10">${barData[i].name}</text>`;
    s += `<text x="${ml + bw + 5}" y="${y + 4}" font-size="9" fill="#333">${barData[i].sign}${barData[i].r.toFixed(3)}</text>`;
  }

  s += svgFooter();
  return s;
}

function fig5_hypothesis() {
  const scores = phase303.hypothesisScores;
  const w = 500, h = 350;
  const ml = 140, mr = 30, mt = 50, mb = 50;
  const pw = w - ml - mr, ph = h - mt - mb;

  let s = svgHeader(w, h, 'Physical Interpretation: Hypothesis Scores');

  const colors = {
    H4_dynIntegration: '#059669',
    H1_haloResponse: '#2563eb',
    H3_feedbackImprint: '#7c3aed',
    H2_assemblyHistory: '#f59e0b'
  };
  const labels = {
    H4_dynIntegration: 'H4: Dynamical Integration',
    H1_haloResponse: 'H1: Halo Response',
    H3_feedbackImprint: 'H3: Feedback Imprint',
    H2_assemblyHistory: 'H2: Assembly History'
  };

  const rowH = ph / scores.length;

  s += `<line x1="${ml}" y1="${mt}" x2="${ml}" y2="${mt+ph}" stroke="black" stroke-width="1"/>`;
  s += `<line x1="${ml}" y1="${mt+ph}" x2="${ml+pw}" y2="${mt+ph}" stroke="black" stroke-width="1"/>`;
  for (let t = 0; t <= 10; t += 2) {
    const x = ml + (t / 12) * pw;
    s += `<line x1="${x}" y1="${mt+ph}" x2="${x}" y2="${mt+ph+5}" stroke="black"/>`;
    s += `<text x="${x}" y="${mt+ph+18}" text-anchor="middle" font-size="10">${t}</text>`;
    s += `<line x1="${x}" y1="${mt}" x2="${x}" y2="${mt+ph}" stroke="#eee" stroke-width="0.5" stroke-dasharray="4,3"/>`;
  }
  s += `<text x="${ml+pw/2}" y="${h-5}" text-anchor="middle" font-size="11">Evidence score</text>`;

  for (let i = 0; i < scores.length; i++) {
    const y = mt + rowH * (i + 0.5);
    const bw = (scores[i].score / 12) * pw;
    const c = colors[scores[i].name] || '#999';
    const label = labels[scores[i].name] || scores[i].name;
    s += `<rect x="${ml}" y="${y - rowH * 0.35}" width="${bw}" height="${rowH * 0.7}" fill="${c}" opacity="0.8" rx="3"/>`;
    s += `<text x="${ml - 5}" y="${y + 4}" text-anchor="end" font-size="11">${escXml(label)}</text>`;
    s += `<text x="${ml + bw + 8}" y="${y + 5}" font-size="13" font-weight="bold" fill="${c}">${scores[i].score}</text>`;
  }

  s += svgFooter();
  return s;
}

function fig6_synthesis() {
  const w = 560, h = 440;
  let s = svgHeader(w, h, 'Synthesis: Hierarchical Coupling Law');

  const boxW = 160, boxH = 55;
  const cx = w / 2;

  const boxes = [
    { x: cx - boxW/2, y: 50, label: 'Structural Core', sub: 'logMHI + logMhost + MR', sub2: 'LOO = 44.1%', color: '#dbeafe', border: '#2563eb' },
    { x: cx - boxW/2, y: 145, label: 'VfResid Channel', sub: 'Kinematic residual', sub2: '+17.1 pp \u2192 61.1%', color: '#ede9fe', border: '#7c3aed' },
    { x: cx - boxW/2, y: 240, label: 'lhOuter Channel', sub: 'Outer halo quality', sub2: '+4.3 pp \u2192 65.4%', color: '#d1fae5', border: '#059669' },
    { x: cx - boxW/2, y: 340, label: 'Regime Gate', sub: 'Vflat \u2265 120 km/s activates', sub2: 'External: P > 99.6%', color: '#fef3c7', border: '#f59e0b' }
  ];

  for (const b of boxes) {
    s += `<rect x="${b.x}" y="${b.y}" width="${boxW}" height="${boxH+10}" rx="6" fill="${b.color}" stroke="${b.border}" stroke-width="2"/>`;
    s += `<text x="${b.x+boxW/2}" y="${b.y+18}" text-anchor="middle" font-size="12" font-weight="bold">${escXml(b.label)}</text>`;
    s += `<text x="${b.x+boxW/2}" y="${b.y+34}" text-anchor="middle" font-size="10" fill="#555">${escXml(b.sub)}</text>`;
    s += `<text x="${b.x+boxW/2}" y="${b.y+50}" text-anchor="middle" font-size="10" fill="#333" font-weight="bold">${escXml(b.sub2)}</text>`;
  }

  for (let i = 0; i < boxes.length - 1; i++) {
    const fromY = boxes[i].y + boxH + 10;
    const toY = boxes[i+1].y;
    s += `<line x1="${cx}" y1="${fromY}" x2="${cx}" y2="${toY}" stroke="#333" stroke-width="1.5" marker-end="url(#arrowDown)"/>`;
  }

  s += `<defs><marker id="arrowDown" markerWidth="8" markerHeight="8" refX="4" refY="7" orient="auto"><path d="M0,0 L4,7 L8,0" fill="none" stroke="#333" stroke-width="1"/></marker></defs>`;

  const sideX = cx + boxW/2 + 30;
  s += `<text x="${sideX}" y="${boxes[1].y + 20}" font-size="10" fill="#7c3aed">r(haloK) = 0.60</text>`;
  s += `<text x="${sideX}" y="${boxes[1].y + 35}" font-size="10" fill="#7c3aed">36\u201369% irreducible</text>`;
  s += `<text x="${sideX}" y="${boxes[1].y + 50}" font-size="10" fill="#7c3aed">Mediates all halo proxies</text>`;
  s += `<line x1="${cx+boxW/2}" y1="${boxes[1].y+boxH/2+5}" x2="${sideX-5}" y2="${boxes[1].y+35}" stroke="#7c3aed" stroke-width="1" stroke-dasharray="3,2"/>`;

  const sideX2 = cx - boxW/2 - 140;
  s += `<text x="${sideX2}" y="${boxes[3].y + 20}" font-size="10" fill="#b45309">H4: Dynamical Integration</text>`;
  s += `<text x="${sideX2}" y="${boxes[3].y + 35}" font-size="10" fill="#b45309">Score = 10 (best)</text>`;
  s += `<text x="${sideX2}" y="${boxes[3].y + 50}" font-size="10" fill="#b45309">slope ratio 2.6:1</text>`;
  s += `<line x1="${cx-boxW/2}" y1="${boxes[3].y+boxH/2+5}" x2="${sideX2+140}" y2="${boxes[3].y+35}" stroke="#b45309" stroke-width="1" stroke-dasharray="3,2"/>`;

  s += svgFooter();
  return s;
}

console.log('Generating 6 publication figures...');

const figures = [
  { name: 'fig1_hierarchy.svg', fn: fig1_hierarchy, desc: 'Model hierarchy progression' },
  { name: 'fig2_external.svg', fn: fig2_external, desc: 'External validation transfer' },
  { name: 'fig3_regime.svg', fn: fig3_regime, desc: 'Regime dependence' },
  { name: 'fig4_vfresid.svg', fn: fig4_vfresid_scatter, desc: 'VfResid driver analysis' },
  { name: 'fig5_hypothesis.svg', fn: fig5_hypothesis, desc: 'Physical interpretation scores' },
  { name: 'fig6_synthesis.svg', fn: fig6_synthesis, desc: 'Synthesis flow diagram' }
];

for (const fig of figures) {
  const svg = fig.fn();
  fs.writeFileSync(path.join(outDir, fig.name), svg);
  console.log(`  [OK] ${fig.name} — ${fig.desc}`);
}

console.log(`\nAll figures written to ${outDir}`);
console.log('\nConverting SVGs to PDFs for LaTeX inclusion...');

const { execSync } = require('child_process');
let convertOK = true;
for (const fig of figures) {
  const svgPath = path.join(outDir, fig.name);
  const pdfName = fig.name.replace('.svg', '.pdf');
  const pdfPath = path.join(outDir, pdfName);
  try {
    execSync(`rsvg-convert -f pdf -o "${pdfPath}" "${svgPath}"`, { timeout: 10000 });
    console.log(`  [OK] ${pdfName}`);
  } catch (e) {
    console.log(`  [SKIP] ${pdfName} — rsvg-convert not available, using SVG`);
    convertOK = false;
    break;
  }
}

if (!convertOK) {
  try {
    execSync('which inkscape', { timeout: 5000 });
    for (const fig of figures) {
      const svgPath = path.join(outDir, fig.name);
      const pdfName = fig.name.replace('.svg', '.pdf');
      const pdfPath = path.join(outDir, pdfName);
      try {
        execSync(`inkscape "${svgPath}" --export-filename="${pdfPath}" --export-type=pdf`, { timeout: 30000 });
        console.log(`  [OK] ${pdfName} (inkscape)`);
      } catch(e2) {
        console.log(`  [SKIP] ${pdfName}`);
      }
    }
  } catch(e) {
    console.log('  No PDF converter available. SVG files created successfully.');
    console.log('  To include in LaTeX, convert SVGs to PDF externally or use svg package.');
  }
}

console.log('\nDone!');
