const fs = require('fs');
const mathjs = require('/home/runner/workspace/artifacts/galaxy-analyzer/node_modules/mathjs');

const galaxies = JSON.parse(fs.readFileSync('/tmp/sparc_parsed.json', 'utf8'));
const G = 4.3009e-6;

const MODELS = {
  newtonian: {
    name: 'Newtonian',
    eval: (r, M, k, a) => { const v2 = (G * M) / r; return v2 > 0 ? Math.sqrt(v2) : null; },
    hasK: false, hasA: false,
    description: 'v = sqrt(GM/r)'
  },
  dark_halo_linear: {
    name: 'Dark Halo (Linear)',
    eval: (r, M, k, a) => { const v2 = (G * M) / r + k * r; return v2 > 0 ? Math.sqrt(v2) : null; },
    hasK: true, hasA: false,
    description: 'v = sqrt(GM/r + kr)'
  },
  dark_halo_flat: {
    name: 'Dark Halo (Flat/Isothermal)',
    eval: (r, M, k, a) => { const v2 = (G * M) / r + k; return v2 > 0 ? Math.sqrt(v2) : null; },
    hasK: true, hasA: false,
    description: 'v = sqrt(GM/r + k)'
  },
  modified_gravity_halo: {
    name: 'Modified Gravity + Halo',
    eval: (r, M, k, a) => { const v2 = (G * M) / (r + a) + k * r; return v2 > 0 ? Math.sqrt(v2) : null; },
    hasK: true, hasA: true,
    description: 'v = sqrt(GM/(r+a) + kr)'
  },
  mond: {
    name: 'MOND-inspired',
    eval: (r, M, k, a) => { const gN = (G * M) / (r * r); const v2 = Math.sqrt(gN * a) * r + (G * M) / r; return v2 > 0 ? Math.sqrt(v2) : null; },
    hasK: false, hasA: true,
    description: 'v = sqrt(sqrt(gN*a₀)*r + GM/r)'
  },
  log_halo: {
    name: 'Logarithmic Halo (NFW-like)',
    eval: (r, M, k, a) => { const v2 = (G * M) / r + k * Math.log(1 + r / a); return v2 > 0 ? Math.sqrt(v2) : null; },
    hasK: true, hasA: true,
    description: 'v = sqrt(GM/r + k*ln(1+r/a))'
  },
  transition: {
    name: 'Transition Model',
    eval: (r, M, k, a) => { const vN2 = (G * M) / r; if (vN2 <= 0) return null; const vN = Math.sqrt(vN2); return vN * (1 + k * r / 100); },
    hasK: true, hasA: false,
    description: 'v = sqrt(GM/r)*(1+kr/100)'
  }
};

function computeMSE(modelEval, data, M, k, a) {
  let sum = 0, count = 0;
  for (let i = 0; i < data.length; i++) {
    const p = data[i];
    const pred = modelEval(p.r, M, k, a);
    if (pred !== null && isFinite(pred)) { sum += (p.v - pred) ** 2; count++; }
  }
  return count > 0 ? sum / count : Infinity;
}

function optimizeModel(modelKey, data) {
  const model = MODELS[modelKey];
  const evalFn = model.eval;
  let bestM = 1e11, bestK = 0, bestA = 1, bestMSE = Infinity;

  const mVals = [1e7, 5e7, 1e8, 5e8, 1e9, 5e9, 1e10, 5e10, 1e11, 5e11, 1e12];
  const kVals = model.hasK ? [0, 1, 5, 10, 25, 50, 100, 200, 500, 1000, 3000] : [0];
  const aVals = model.hasA ? [0.05, 0.2, 0.5, 1, 3, 5, 10, 20, 50] : [1];

  for (const tM of mVals) {
    for (const tK of kVals) {
      for (const tA of aVals) {
        const mse = computeMSE(evalFn, data, tM, tK, tA);
        if (mse < bestMSE) { bestMSE = mse; bestM = tM; bestK = tK; bestA = tA; }
      }
    }
  }

  for (let pass = 0; pass < 2; pass++) {
    const fineM = [];
    for (let i = -3; i <= 3; i++) fineM.push(bestM * (1 + i * 0.2));
    const fineK = model.hasK ? [] : [0];
    if (model.hasK) {
      for (let i = -5; i <= 5; i++) fineK.push(Math.max(0, bestK + i * Math.max(bestK * 0.1, 1)));
    }
    const fineA = model.hasA ? [] : [1];
    if (model.hasA) {
      for (let i = -4; i <= 4; i++) fineA.push(Math.max(0.001, bestA + i * Math.max(bestA * 0.12, 0.1)));
    }

    for (const tM of fineM) {
      for (const tK of fineK) {
        for (const tA of fineA) {
          const mse = computeMSE(evalFn, data, tM, tK, tA);
          if (mse < bestMSE) { bestMSE = mse; bestM = tM; bestK = tK; bestA = tA; }
        }
      }
    }
  }

  return { M: bestM, k: bestK, a: bestA, mse: bestMSE };
}

function displayName(raw) {
  const ngcMatch = raw.match(/^NGC(\d+)$/);
  if (ngcMatch) return 'NGC ' + parseInt(ngcMatch[1]);
  const ugcMatch = raw.match(/^UGC(\d+)$/);
  if (ugcMatch) return 'UGC ' + parseInt(ugcMatch[1]);
  const icMatch = raw.match(/^IC(\d+)$/);
  if (icMatch) return 'IC ' + parseInt(icMatch[1]);
  const ddoMatch = raw.match(/^DDO(\d+)$/);
  if (ddoMatch) return 'DDO ' + parseInt(ddoMatch[1]);
  return raw;
}

console.log('=== FULL SPARC ANALYSIS ===');
console.log(`Testing ${Object.keys(MODELS).length} models on ${galaxies.length} galaxies`);
console.log(`G = ${G} kpc·(km/s)²/M☉`);

const modelKeys = Object.keys(MODELS);
const allResults = [];
const modelStats = {};
for (const mk of modelKeys) {
  modelStats[mk] = { kValues: [], aValues: [], mValues: [], mseValues: [], improvements: [], wins: 0, innerImpr: [], outerImpr: [] };
}

const startTime = Date.now();
let processed = 0;

for (const galaxy of galaxies) {
  processed++;
  const name = displayName(galaxy.name);
  const data = galaxy.data;

  if (processed % 10 === 0 || processed === 1 || processed === galaxies.length) {
    const elapsed = ((Date.now() - startTime) / 1000).toFixed(0);
    process.stdout.write(`[${elapsed}s] ${processed}/${galaxies.length}: ${name.padEnd(20)}\r`);
  }

  const results = {};
  for (const mk of modelKeys) {
    results[mk] = optimizeModel(mk, data);
  }

  const newtMSE = results.newtonian.mse;

  const mid = Math.floor(data.length / 2);
  const innerData = data.slice(0, mid);
  const outerData = data.slice(mid);

  const galaxyRow = {
    name,
    rawName: galaxy.name,
    distance: galaxy.distance,
    pointCount: data.length,
    maxR: galaxy.maxR,
    maxV: galaxy.maxV,
    models: {}
  };

  for (const mk of modelKeys) {
    const r = results[mk];
    const improvPct = newtMSE > 0 && isFinite(newtMSE) ? ((newtMSE - r.mse) / newtMSE * 100) : 0;

    let mseInner = 0, mseOuter = 0, cIn = 0, cOut = 0;
    let mseInnerN = 0, mseOuterN = 0;

    innerData.forEach(p => {
      const pc = MODELS[mk].eval(p.r, r.M, r.k, r.a);
      const pn = MODELS.newtonian.eval(p.r, results.newtonian.M, 0, 1);
      if (pc !== null) { mseInner += (p.v - pc) ** 2; cIn++; }
      if (pn !== null) { mseInnerN += (p.v - pn) ** 2; }
    });
    outerData.forEach(p => {
      const pc = MODELS[mk].eval(p.r, r.M, r.k, r.a);
      const pn = MODELS.newtonian.eval(p.r, results.newtonian.M, 0, 1);
      if (pc !== null) { mseOuter += (p.v - pc) ** 2; cOut++; }
      if (pn !== null) { mseOuterN += (p.v - pn) ** 2; }
    });
    mseInner = cIn > 0 ? mseInner / cIn : 0;
    mseOuter = cOut > 0 ? mseOuter / cOut : 0;
    mseInnerN = cIn > 0 ? mseInnerN / cIn : 0;
    mseOuterN = cOut > 0 ? mseOuterN / cOut : 0;

    const innerImpr = mseInnerN > 0 ? ((mseInnerN - mseInner) / mseInnerN * 100) : 0;
    const outerImpr = mseOuterN > 0 ? ((mseOuterN - mseOuter) / mseOuterN * 100) : 0;

    galaxyRow.models[mk] = {
      M: r.M, k: r.k, a: r.a, mse: r.mse,
      improvementVsNewton: improvPct,
      mseInner, mseOuter, innerImprovement: innerImpr, outerImprovement: outerImpr
    };

    if (mk !== 'newtonian') {
      modelStats[mk].kValues.push(r.k);
      modelStats[mk].aValues.push(r.a);
      modelStats[mk].mValues.push(r.M);
      modelStats[mk].mseValues.push(r.mse);
      modelStats[mk].improvements.push(improvPct);
      if (improvPct > 5) modelStats[mk].wins++;
      modelStats[mk].innerImpr.push(innerImpr);
      modelStats[mk].outerImpr.push(outerImpr);
    }
  }

  allResults.push(galaxyRow);
}

const totalTime = ((Date.now() - startTime) / 1000).toFixed(1);
console.log(`\nDone in ${totalTime}s. Processing ${galaxies.length} galaxies.\n`);

function mean(arr) { return arr.length > 0 ? arr.reduce((s, v) => s + v, 0) / arr.length : 0; }
function stddev(arr) { const m = mean(arr); return Math.sqrt(arr.reduce((s, v) => s + (v - m) ** 2, 0) / arr.length); }
function median(arr) { const s = [...arr].sort((a, b) => a - b); const mid = Math.floor(s.length / 2); return s.length % 2 ? s[mid] : (s[mid - 1] + s[mid]) / 2; }
function cv(arr) { const m = mean(arr); return m > 0 ? (stddev(arr) / m * 100) : Infinity; }

let report = '';
report += '╔══════════════════════════════════════════════════════════════════════════════════╗\n';
report += '║           تحليل منحنيات الدوران — قاعدة بيانات SPARC الكاملة                    ║\n';
report += '║           GALAXY ROTATION CURVE ANALYSIS — FULL SPARC DATABASE                  ║\n';
report += '║           175 Galaxies • 3,391 Data Points • 7 Models                           ║\n';
report += '╚══════════════════════════════════════════════════════════════════════════════════╝\n\n';

report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n';
report += '1. METHODOLOGY — المنهجية\n';
report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n';
report += 'Data Source: SPARC (Spitzer Photometry and Accurate Rotation Curves)\n';
report += '             Lelli, McGaugh & Schombert (2016)\n';
report += `             ${galaxies.length} late-type galaxies, ${galaxies.reduce((s, g) => s + g.pointCount, 0)} total data points\n\n`;
report += `Gravitational Constant: G = ${G} kpc·(km/s)²/M☉ (physical value)\n\n`;
report += 'Units:\n';
report += '  r — kiloparsecs (kpc), 1 kpc = 3,260 light-years\n';
report += '  v — kilometers per second (km/s)\n';
report += '  M — solar masses (M☉)\n';
report += '  k — fitting parameter (units depend on formula)\n';
report += '  a — core radius / acceleration scale\n\n';
report += 'Optimization: 2-phase grid search per galaxy per model\n';
report += '  Newtonian baseline optimized independently per galaxy (fair comparison)\n\n';

report += 'الادعاء المُختبَر — Claim Under Test:\n';
report += '  "إضافة حد يعتمد على المسافة يحسّن fit على عدد كبير من المجرات،\n';
report += '   وبمعامل k شبه ثابت — مما يشير إلى تأثير فيزيائي حقيقي."\n\n';
report += '  "Adding a distance-dependent term improves fit across many galaxies,\n';
report += '   with a roughly consistent parameter k."\n\n';

report += 'Models Tested — النماذج المُختبرة:\n';
for (const mk of modelKeys) {
  const m = MODELS[mk];
  report += `  ${m.name.padEnd(32)} ${m.description}\n`;
}
report += '\n';

report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n';
report += '2. MODEL COMPARISON SUMMARY — ملخص مقارنة النماذج\n';
report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n';

const summaryHeader = 'Model'.padEnd(32) + 'Avg MSE'.padStart(12) + 'Med MSE'.padStart(12) + 'Avg Impr%'.padStart(11) + 'Wins/175'.padStart(10) + 'Avg k'.padStart(10) + 'k CV%'.padStart(8) + 'Stable?'.padStart(10) + '\n';
report += summaryHeader;
report += '─'.repeat(105) + '\n';

const rankings = [];
for (const mk of modelKeys) {
  if (mk === 'newtonian') {
    report += MODELS[mk].name.padEnd(32) + '(baseline)'.padStart(12) + ''.padStart(12) + '0.0%'.padStart(11) + '—'.padStart(10) + '—'.padStart(10) + '—'.padStart(8) + '—'.padStart(10) + '\n';
    continue;
  }
  const s = modelStats[mk];
  const avgMSE = mean(s.mseValues);
  const medMSE = median(s.mseValues);
  const avgImpr = mean(s.improvements);
  const kCV = model_hasK(mk) ? cv(s.kValues) : 0;
  const kStable = !model_hasK(mk) ? 'N/A' : kCV < 20 ? 'STABLE' : kCV < 50 ? 'moderate' : kCV < 80 ? 'variable' : 'UNSTABLE';

  rankings.push({ mk, avgImpr, wins: s.wins, kCV, kStable, avgMSE });

  report += MODELS[mk].name.padEnd(32) +
    avgMSE.toFixed(1).padStart(12) +
    medMSE.toFixed(1).padStart(12) +
    (avgImpr.toFixed(1) + '%').padStart(11) +
    `${s.wins}/175`.padStart(10) +
    (model_hasK(mk) ? mean(s.kValues).toFixed(1) : '—').padStart(10) +
    (model_hasK(mk) ? (kCV.toFixed(0) + '%') : '—').padStart(8) +
    kStable.padStart(10) + '\n';
}

function model_hasK(mk) { return MODELS[mk].hasK; }

rankings.sort((a, b) => b.avgImpr - a.avgImpr);
report += '\n';
report += `🏆 أفضل نموذج — Best Model: ${MODELS[rankings[0].mk].name}\n`;
report += `   Average Improvement: ${rankings[0].avgImpr.toFixed(1)}%\n`;
report += `   Wins: ${rankings[0].wins}/175 galaxies\n\n`;

report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n';
report += '3. REGIONAL ANALYSIS — التحليل الإقليمي (داخلي مقابل خارجي)\n';
report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n';
report += '"إذا كان التحسن يأتي من الأطراف (المناطق الخارجية)، فهذا يتوافق مع مشكلة\n';
report += ' الدوران المجري حيث تفشل جاذبية نيوتن عند المسافات الكبيرة."\n\n';

const regHeader = 'Model'.padEnd(32) + 'Inner Impr%'.padStart(14) + 'Outer Impr%'.padStart(14) + 'Outer>Inner?'.padStart(14) + '\n';
report += regHeader;
report += '─'.repeat(74) + '\n';

for (const mk of modelKeys) {
  if (mk === 'newtonian') continue;
  const s = modelStats[mk];
  const avgIn = mean(s.innerImpr);
  const avgOut = mean(s.outerImpr);
  const verdict = avgOut > avgIn * 1.2 ? 'YES ✓' : avgIn > avgOut * 1.2 ? 'NO (inner)' : 'BALANCED';

  report += MODELS[mk].name.padEnd(32) +
    (avgIn.toFixed(1) + '%').padStart(14) +
    (avgOut.toFixed(1) + '%').padStart(14) +
    verdict.padStart(14) + '\n';
}
report += '\n';

report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n';
report += '4. k PARAMETER STABILITY — استقرار معامل k (Dark Halo Linear)\n';
report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n';

const dhlK = modelStats.dark_halo_linear.kValues;
report += `Mean k:     ${mean(dhlK).toFixed(2)}\n`;
report += `Median k:   ${median(dhlK).toFixed(2)}\n`;
report += `Std Dev:    ${stddev(dhlK).toFixed(2)}\n`;
report += `CV:         ${cv(dhlK).toFixed(1)}%\n`;
report += `Min k:      ${Math.min(...dhlK).toFixed(2)}\n`;
report += `Max k:      ${Math.max(...dhlK).toFixed(2)}\n\n`;

const kBins = [0, 5, 10, 25, 50, 100, 200, 500, 1000, Infinity];
report += 'توزيع k — k Distribution:\n';
for (let i = 0; i < kBins.length - 1; i++) {
  const lo = kBins[i], hi = kBins[i + 1];
  const count = dhlK.filter(v => v >= lo && v < hi).length;
  const bar = '█'.repeat(Math.round(count / galaxies.length * 60));
  const label = hi === Infinity ? `${lo}+` : `${lo}-${hi}`;
  report += `  ${label.padEnd(10)} ${bar} ${count} (${(count / galaxies.length * 100).toFixed(0)}%)\n`;
}
report += '\n';

const kVerdict = cv(dhlK) < 20
  ? '✓ k ثابت — يشير إلى معامل فيزيائي حقيقي وليس مجرد fitting.\n  ✓ k IS CONSISTENT — suggests a real physical parameter.'
  : cv(dhlK) < 50
    ? '~ k متوسط الاستقرار — يوجد نمط لكن مع تباين ملحوظ.\n  ~ k is MODERATELY CONSISTENT — pattern exists but with variation.'
    : '✗ k غير مستقر — يتغير كثيرًا بين المجرات؛ قد يكون يمتص خصائص خاصة بكل مجرة.\n  ✗ k is NOT CONSISTENT — varies too much; may absorb galaxy-specific structure.';
report += `الحكم — Verdict: ${kVerdict}\n\n`;

report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n';
report += '5. HEAD-TO-HEAD — مقارنة مباشرة: Halo vs MOND vs Log Halo\n';
report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n';

const contenders = ['dark_halo_linear', 'mond', 'log_halo', 'dark_halo_flat', 'modified_gravity_halo'];
let h2hHeader = 'Galaxy'.padEnd(20);
for (const c of contenders) h2hHeader += MODELS[c].name.substring(0, 14).padStart(16);
h2hHeader += '  Winner\n';
report += h2hHeader;
report += '─'.repeat(h2hHeader.length) + '\n';

const h2hWins = {};
for (const c of contenders) h2hWins[c] = 0;

for (const row of allResults) {
  let line = row.name.substring(0, 19).padEnd(20);
  let bestC = contenders[0], bestMSE = Infinity;
  for (const c of contenders) {
    const mse = row.models[c].mse;
    line += mse.toFixed(1).padStart(16);
    if (mse < bestMSE) { bestMSE = mse; bestC = c; }
  }
  h2hWins[bestC]++;
  line += '  ' + MODELS[bestC].name.substring(0, 16);
  report += line + '\n';
}

report += '\nعدد الانتصارات — Win Counts:\n';
for (const c of contenders) {
  const pct = (h2hWins[c] / galaxies.length * 100).toFixed(1);
  report += `  ${MODELS[c].name.padEnd(32)} ${h2hWins[c].toString().padStart(4)} wins (${pct}%)\n`;
}
report += '\n';

report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n';
report += '6. FIVE SIGNS — العلامات الخمس لنظرية واعدة\n';
report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n';

const bestModel = rankings[0];
const bStats = modelStats[bestModel.mk];

const sign1 = bestModel.wins > galaxies.length * 0.6;
const sign2 = model_hasK(bestModel.mk) ? cv(bStats.kValues) < 50 : true;
const avgOuter = mean(bStats.outerImpr);
const avgInner = mean(bStats.innerImpr);
const sign3 = avgOuter > avgInner;
const numParams = (MODELS[bestModel.mk].hasK ? 1 : 0) + (MODELS[bestModel.mk].hasA ? 1 : 0) + 1;
const sign4 = numParams <= 3;

let competitiveCount = 0;
for (const row of allResults) {
  const bestMSE = row.models[bestModel.mk].mse;
  const mondMSE = row.models.mond.mse;
  const logMSE = row.models.log_halo.mse;
  if (bestMSE <= mondMSE * 1.1 && bestMSE <= logMSE * 1.1) competitiveCount++;
}
const sign5 = competitiveCount > galaxies.length * 0.5;

report += `النموذج الأفضل — Best model: ${MODELS[bestModel.mk].name}\n\n`;
report += `${sign1 ? '✓' : '✗'} 1. التحسن على عدد كبير من المجرات — Improvement on many galaxies\n`;
report += `     ${bestModel.wins}/175 (${(bestModel.wins / 175 * 100).toFixed(0)}%) — Threshold: >60% → ${sign1 ? 'PASS' : 'FAIL'}\n\n`;
report += `${sign2 ? '✓' : '✗'} 2. k مستقر — k is stable\n`;
report += `     CV = ${model_hasK(bestModel.mk) ? cv(bStats.kValues).toFixed(1) + '%' : 'N/A'} — Threshold: <50% → ${sign2 ? 'PASS' : 'FAIL'}\n\n`;
report += `${sign3 ? '✓' : '✗'} 3. التحسن عند الأطراف — Improvement at outer radii\n`;
report += `     Inner ${avgInner.toFixed(1)}% vs Outer ${avgOuter.toFixed(1)}% → ${sign3 ? 'PASS' : 'FAIL'}\n\n`;
report += `${sign4 ? '✓' : '✗'} 4. عدد قليل من المعاملات — Few parameters\n`;
report += `     ${numParams} parameters → ${sign4 ? 'PASS' : 'FAIL'}\n\n`;
report += `${sign5 ? '✓' : '✗'} 5. منافس لـ MOND/Halo — Competitive with MOND/Log Halo\n`;
report += `     Wins or ties on ${competitiveCount}/175 → ${sign5 ? 'PASS' : 'FAIL'}\n\n`;

const passCount = [sign1, sign2, sign3, sign4, sign5].filter(Boolean).length;
report += `━━ النتيجة — OVERALL: ${passCount}/5 ━━\n`;
if (passCount >= 4) report += '→ النظرية واعدة جدًا. تستحق تحقيقًا أعمق.\n→ Theory shows STRONG promise. Warrants deeper investigation.\n';
else if (passCount >= 3) report += '→ النظرية واعدة بشكل معتدل. بعض المخاوف تبقى.\n→ Theory shows MODERATE promise. Some concerns remain.\n';
else if (passCount >= 2) report += '→ النظرية ضعيفة الوعد. مشاكل كبيرة يجب معالجتها.\n→ Theory shows WEAK promise. Significant issues to address.\n';
else report += '→ النظرية لا تحقق المعايير الأساسية. غالبًا نتيجة مرونة التكييف فقط.\n→ Theory does NOT meet basic criteria. Likely curve-fitting artifact.\n';
report += '\n';

report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n';
report += '7. DETAILED PER-GALAXY RESULTS — نتائج كل مجرة\n';
report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n';

const detailHeader = 'Galaxy'.padEnd(18) + 'Dist'.padStart(7) + 'Pts'.padStart(5) + 'maxR'.padStart(7) + 'maxV'.padStart(7) +
  ' | Newt'.padStart(9) + 'Halo'.padStart(9) + 'MOND'.padStart(9) + 'LogH'.padStart(9) +
  ' | k'.padStart(8) + 'M'.padStart(10) + 'Impr%'.padStart(8) + '\n';
report += detailHeader;
report += '─'.repeat(106) + '\n';

for (const row of allResults) {
  const n = row.models.newtonian;
  const h = row.models.dark_halo_linear;
  const m = row.models.mond;
  const l = row.models.log_halo;

  report += row.name.substring(0, 17).padEnd(18) +
    (row.distance ? row.distance.toFixed(1) : '—').padStart(7) +
    row.pointCount.toString().padStart(5) +
    row.maxR.toFixed(1).padStart(7) +
    row.maxV.toFixed(0).padStart(7) +
    ' |' +
    n.mse.toFixed(0).padStart(8) +
    h.mse.toFixed(0).padStart(9) +
    m.mse.toFixed(0).padStart(9) +
    l.mse.toFixed(0).padStart(9) +
    ' |' +
    h.k.toFixed(1).padStart(7) +
    h.M.toExponential(1).padStart(10) +
    (h.improvementVsNewton.toFixed(1) + '%').padStart(8) + '\n';
}
report += '\n';

report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n';
report += '8. CONCLUSION — الاستنتاج\n';
report += '━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n';

report += `Analysis Date: ${new Date().toISOString()}\n`;
report += `Computation Time: ${totalTime}s\n`;
report += `Galaxies: ${galaxies.length}, Points: ${galaxies.reduce((s, g) => s + g.pointCount, 0)}\n\n`;

report += `Best: ${MODELS[rankings[0].mk].name}\n`;
report += `  Avg Improvement: ${rankings[0].avgImpr.toFixed(1)}%\n`;
report += `  Wins: ${rankings[0].wins}/175\n\n`;

report += `Runner-up: ${MODELS[rankings[1].mk].name}\n`;
report += `  Avg Improvement: ${rankings[1].avgImpr.toFixed(1)}%\n`;
report += `  Wins: ${rankings[1].wins}/175\n\n`;

report += `MOND:\n`;
report += `  Avg Improvement: ${mean(modelStats.mond.improvements).toFixed(1)}%\n`;
report += `  Wins: ${modelStats.mond.wins}/175\n\n`;

if (passCount >= 3) {
  report += 'التقييم — ASSESSMENT:\n';
  report += '  التصحيح المعتمد على المسافة يُظهر تحسنًا ذا معنى عبر جزء كبير من عينة SPARC.\n';
  report += '  The distance-dependent correction shows meaningful improvement.\n\n';
  report += '  الخطوات التالية:\n';
  report += '  1. هل يرتبط k بخصائص المجرة (اللمعان، الكتلة، المورفولوجيا)?\n';
  report += '  2. اختبار على بيانات مستقلة (THINGS, LITTLE THINGS)\n';
  report += '  3. التحقق من القيود الكونية (CMB, BAO, LSS)\n';
  report += '  4. اشتقاق الأصل النظري للحد الإضافي\n';
} else {
  report += 'التقييم — ASSESSMENT:\n';
  report += '  النموذج يُظهر بعض التحسن لكنه لا يتفوق بوضوح على البدائل الموجودة.\n';
  report += '  The model shows improvement but does not clearly outperform existing alternatives.\n';
}

report += '\n═══════════════════════════════════════════════════════════════════════════════════\n';
report += '                              END OF REPORT — نهاية التقرير\n';
report += '═══════════════════════════════════════════════════════════════════════════════════\n';

const outputPath = '/home/runner/workspace/artifacts/galaxy-analyzer/public/sparc-full-analysis.txt';
fs.writeFileSync(outputPath, report);
console.log(`Report: ${outputPath} (${(report.length / 1024).toFixed(1)} KB)`);

const jsonPath = '/home/runner/workspace/artifacts/galaxy-analyzer/public/sparc-results.json';
fs.writeFileSync(jsonPath, JSON.stringify({
  metadata: { date: new Date().toISOString(), galaxyCount: galaxies.length, totalPoints: galaxies.reduce((s, g) => s + g.pointCount, 0), computeTime: totalTime },
  summary: Object.fromEntries(modelKeys.filter(mk => mk !== 'newtonian').map(mk => {
    const s = modelStats[mk];
    return [mk, { name: MODELS[mk].name, avgMSE: mean(s.mseValues), avgImprovement: mean(s.improvements), wins: s.wins, avgK: mean(s.kValues), kCV: cv(s.kValues) }];
  })),
  perGalaxy: allResults,
  rankings: rankings.map(r => ({ model: MODELS[r.mk].name, avgImprovement: r.avgImpr, wins: r.wins })),
  fiveSigns: { sign1, sign2, sign3, sign4, sign5, passCount }
}, null, 2));
console.log(`JSON: ${jsonPath}`);

console.log(`\n=== QUICK SUMMARY ===`);
console.log(`Best: ${MODELS[rankings[0].mk].name} — ${rankings[0].avgImpr.toFixed(1)}% avg, ${rankings[0].wins}/175 wins`);
console.log(`2nd:  ${MODELS[rankings[1].mk].name} — ${rankings[1].avgImpr.toFixed(1)}% avg, ${rankings[1].wins}/175 wins`);
console.log(`MOND: ${mean(modelStats.mond.improvements).toFixed(1)}% avg, ${modelStats.mond.wins}/175 wins`);
console.log(`Criteria: ${passCount}/5`);
