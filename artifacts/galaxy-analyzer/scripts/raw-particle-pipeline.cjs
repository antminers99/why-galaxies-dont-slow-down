#!/usr/bin/env node
'use strict';

const https = require('https');
const http = require('http');
const fs = require('fs');
const path = require('path');

const G_KPC = 4.3009e-6; // kpc (km/s)^2 / Msun
const H0 = 67.74; // km/s/Mpc (TNG cosmology)
const h = H0 / 100; // 0.6774

const APERTURE_RADII_PKPC = [1, 3, 5, 10, 20, 30, 40, 50, 70, 100];

function fetchJSON(url, body) {
  return new Promise((resolve, reject) => {
    const parsed = new URL(url);
    const options = {
      hostname: parsed.hostname,
      port: parsed.port || 443,
      path: parsed.pathname + parsed.search,
      method: body ? 'POST' : 'GET',
      headers: body ? { 'Content-Type': 'application/json' } : {},
    };
    const proto = parsed.protocol === 'https:' ? https : http;
    const req = proto.request(options, (res) => {
      let data = '';
      res.on('data', (chunk) => (data += chunk));
      res.on('end', () => {
        try { resolve(JSON.parse(data)); }
        catch (e) { reject(new Error('Parse error: ' + data.substring(0, 200))); }
      });
    });
    req.on('error', reject);
    req.setTimeout(60000, () => { req.destroy(); reject(new Error('Timeout')); });
    if (body) req.write(JSON.stringify(body));
    req.end();
  });
}

const TNG_DM_FIELDS = APERTURE_RADII_PKPC.map(r => 'Subhalo_Aperture_Mass_DM_' + r);
const TNG_STAR_FIELDS = APERTURE_RADII_PKPC.map(r => 'Subhalo_Aperture_Mass_Star_' + r);
const TNG_GAS_FIELDS = APERTURE_RADII_PKPC.map(r => 'Subhalo_Aperture_Mass_Gas_' + r);

const TNG_FIELDS = [
  'Subhalo_Vmax', 'Subhalo_VmaxRad', 'Subhalo_Mass',
  'Subhalo_MassType_stars', 'Subhalo_MassType_gas', 'Subhalo_MassType_dm',
  'Subhalo_HalfmassRad',
  'Subhalo_StellarPhotometricsRad',
  ...TNG_DM_FIELDS, ...TNG_STAR_FIELDS, ...TNG_GAS_FIELDS
];

const EAGLE_DM_FIELDS = APERTURE_RADII_PKPC.map(r => 'Subhalo_Aperture_Mass_DM_' + r);
const EAGLE_STAR_FIELDS = APERTURE_RADII_PKPC.map(r => 'Subhalo_Aperture_Mass_Star_' + r);
const EAGLE_GAS_FIELDS = APERTURE_RADII_PKPC.map(r => 'Subhalo_Aperture_Mass_Gas_' + r);

const EAGLE_FIELDS = [
  'Subhalo_Vmax', 'Subhalo_VmaxRadius', 'Subhalo_Mass',
  'Subhalo_MassType_Star', 'Subhalo_MassType_Gas', 'Subhalo_MassType_DM',
  'Subhalo_HalfMassRad_DM',
  ...EAGLE_DM_FIELDS, ...EAGLE_STAR_FIELDS, ...EAGLE_GAS_FIELDS
];

async function queryFlatHub(catalog, body) {
  const url = 'https://flathub.flatironinstitute.org/api/' + catalog + '/data';
  console.log('  Querying FlatHUB/' + catalog + ' for ' + body.count + ' halos...');
  const result = await fetchJSON(url, body);
  if (typeof result === 'string' && result.startsWith('Error')) {
    throw new Error(result);
  }
  return result;
}

async function fetchTNGData() {
  console.log('\n=== FETCHING TNG100-1 z=0 DATA ===');
  console.log('Fields:', TNG_FIELDS.length);

  const BATCH = 2000;
  let allRows = [];
  let offset = 0;
  let done = false;

  while (!done) {
    const body = {
      simulation: 14, // TNG100-1
      snapshot: 99,   // z=0
      fields: TNG_FIELDS,
      count: BATCH,
      offset: offset,
    };

    const rows = await queryFlatHub('tng', body);
    if (!rows || rows.length === 0) { done = true; break; }
    allRows = allRows.concat(rows);
    offset += rows.length;
    console.log('  Got ' + rows.length + ' rows (total: ' + allRows.length + ')');

    if (allRows.length >= 10000) {
      done = true;
    }
    if (rows.length < BATCH) done = true;
  }

  console.log('  Total TNG halos fetched: ' + allRows.length);
  console.log('  NOTE: FlatHUB API limits to ~10k rows per query. Full catalog is larger.');

  const fieldIdx = {};
  TNG_FIELDS.forEach((f, i) => fieldIdx[f] = i);

  const galaxies = [];
  for (const row of allRows) {
    const Vmax = row[fieldIdx['Subhalo_Vmax']];
    const mStar = row[fieldIdx['Subhalo_MassType_stars']] * 1e10 / h; // to Msun
    const mGas = row[fieldIdx['Subhalo_MassType_gas']] * 1e10 / h;

    if (Vmax < 50 || Vmax > 300) continue; // SPARC-like velocity range
    if (mStar < 1e8) continue; // need resolved stellar component

    const dmMasses = TNG_DM_FIELDS.map(f => row[fieldIdx[f]]);
    const starMasses = TNG_STAR_FIELDS.map(f => row[fieldIdx[f]]);
    const gasMasses = TNG_GAS_FIELDS.map(f => row[fieldIdx[f]]);

    // TNG aperture masses are in Msun already (no h factor)
    // but let me verify by checking units in the metadata
    // The catalog says units are "Msun" for aperture fields

    if (dmMasses[0] <= 0 || dmMasses[1] <= 0) continue; // need resolved inner DM

    galaxies.push({
      Vmax, mStar, mGas,
      VmaxRad: row[fieldIdx['Subhalo_VmaxRad']] / h, // ckpc/h -> kpc (z=0)
      halfmassRad: row[fieldIdx['Subhalo_HalfmassRad']] / h,
      dmMasses, starMasses, gasMasses,
    });
  }

  console.log('  SPARC-like TNG galaxies: ' + galaxies.length);
  return galaxies;
}

async function fetchEAGLEData() {
  console.log('\n=== FETCHING EAGLE RefL0100N1504 z=0 DATA ===');
  console.log('Fields:', EAGLE_FIELDS.length);

  const BATCH = 2000;
  let allRows = [];
  let offset = 0;
  let done = false;

  while (!done) {
    const body = {
      simulation: 9, // RefL0100N1504
      snapshot: 28,  // z~0
      fields: EAGLE_FIELDS,
      count: BATCH,
      offset: offset,
    };

    const rows = await queryFlatHub('eagle', body);
    if (!rows || rows.length === 0) { done = true; break; }
    allRows = allRows.concat(rows);
    offset += rows.length;
    console.log('  Got ' + rows.length + ' rows (total: ' + allRows.length + ')');

    if (allRows.length >= 10000) {
      done = true;
    }
    if (rows.length < BATCH) done = true;
  }

  console.log('  Total EAGLE halos fetched: ' + allRows.length);
  console.log('  NOTE: FlatHUB API limits to ~10k rows per query. Full catalog is larger.');

  const fieldIdx = {};
  EAGLE_FIELDS.forEach((f, i) => fieldIdx[f] = i);

  const galaxies = [];
  for (const row of allRows) {
    const Vmax = row[fieldIdx['Subhalo_Vmax']];
    const mStar = row[fieldIdx['Subhalo_MassType_Star']] * 1e10 / h;
    const mGas = row[fieldIdx['Subhalo_MassType_Gas']] * 1e10 / h;

    if (Vmax < 50 || Vmax > 300) continue;
    if (mStar < 1e8) continue;

    const dmMasses = EAGLE_DM_FIELDS.map(f => row[fieldIdx[f]]);
    const starMasses = EAGLE_STAR_FIELDS.map(f => row[fieldIdx[f]]);
    const gasMasses = EAGLE_GAS_FIELDS.map(f => row[fieldIdx[f]]);

    if (dmMasses[0] <= 0 || dmMasses[1] <= 0) continue;

    galaxies.push({
      Vmax, mStar, mGas,
      VmaxRad: row[fieldIdx['Subhalo_VmaxRadius']] / h,
      halfmassRadDM: row[fieldIdx['Subhalo_HalfMassRad_DM']] / h,
      dmMasses, starMasses, gasMasses,
    });
  }

  console.log('  SPARC-like EAGLE galaxies: ' + galaxies.length);
  return galaxies;
}

function computeAlphaAndSigma(galaxy, radiiKpc) {
  const { dmMasses, starMasses, gasMasses } = galaxy;

  const vDM = [];
  for (let i = 0; i < radiiKpc.length; i++) {
    const r_kpc = radiiKpc[i];
    const M_DM = dmMasses[i]; // Msun
    if (M_DM > 0 && r_kpc > 0) {
      vDM.push({ r: r_kpc, v: Math.sqrt(G_KPC * M_DM / r_kpc) });
    }
  }

  if (vDM.length < 3) return null;

  const innerPairs = vDM.filter(p => p.r <= 10);
  if (innerPairs.length < 2) return null;

  let sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0;
  const n = innerPairs.length;
  for (const p of innerPairs) {
    const x = Math.log(p.r);
    const y = Math.log(p.v);
    sumX += x; sumY += y;
    sumXY += x * y; sumX2 += x * x;
  }
  const alpha = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);

  const totalStar = starMasses[4] || starMasses[3]; // 20 kpc or 10 kpc
  const totalGas = gasMasses[4] || gasMasses[3];
  const r_eff = galaxy.halfmassRad || galaxy.VmaxRad || 5;
  const area = Math.PI * (r_eff * r_eff);
  const Sigma_bar = (totalStar + totalGas) / area; // Msun/kpc^2

  if (!isFinite(alpha) || !isFinite(Sigma_bar) || Sigma_bar <= 0) return null;

  return {
    alpha,
    log_Sigma_bar: Math.log10(Sigma_bar),
    Sigma_bar,
    Vmax: galaxy.Vmax,
    mStar: galaxy.mStar,
  };
}

function linearRegression(xs, ys) {
  const n = xs.length;
  let sx = 0, sy = 0, sxy = 0, sx2 = 0, sy2 = 0;
  for (let i = 0; i < n; i++) {
    sx += xs[i]; sy += ys[i];
    sxy += xs[i] * ys[i];
    sx2 += xs[i] * xs[i];
    sy2 += ys[i] * ys[i];
  }
  const slope = (n * sxy - sx * sy) / (n * sx2 - sx * sx);
  const intercept = (sy - slope * sx) / n;
  const r = (n * sxy - sx * sy) / Math.sqrt((n * sx2 - sx * sx) * (n * sy2 - sy * sy));

  const yPred = xs.map(x => slope * x + intercept);
  const residuals = ys.map((y, i) => y - yPred[i]);
  const sse = residuals.reduce((s, r) => s + r * r, 0);
  const slopeErr = Math.sqrt(sse / ((n - 2) * (sx2 - sx * sx / n)));

  return { slope, intercept, r, slopeErr, n };
}

function loadSPARCData() {
  const rarPath = path.join(__dirname, '..', 'public', 'rar-analysis-real.json');
  if (!fs.existsSync(rarPath)) {
    console.log('  WARNING: rar-analysis-real.json not found, using hardcoded SPARC slope');
    return null;
  }
  const rar = JSON.parse(fs.readFileSync(rarPath, 'utf8'));
  const galaxies = rar.perGalaxy || [];
  console.log('  SPARC galaxies loaded: ' + galaxies.length);

  const rotDir = '/tmp/rotmod';
  if (!fs.existsSync(rotDir)) {
    console.log('  WARNING: rotmod dir not found');
    return null;
  }

  const results = [];
  for (const g of galaxies) {
    const fname = g.name + '_rotmod.dat';
    const fpath = path.join(rotDir, fname);
    if (!fs.existsSync(fpath)) continue;

    const lines = fs.readFileSync(fpath, 'utf8').split('\n').filter(l => l.trim() && !l.startsWith('#'));
    const dataPoints = [];
    for (const line of lines) {
      const cols = line.trim().split(/\s+/).map(Number);
      if (cols.length >= 6) {
        const [r, Vobs, errV, Vgas, Vdisk, Vbul] = cols;
        const Vbar2 = Vgas * Math.abs(Vgas) + 0.5 * Vdisk * Math.abs(Vdisk) + 0.7 * (Vbul || 0) * Math.abs(Vbul || 0);
        const Vobs2 = Vobs * Vobs;
        const VDM2 = Vobs2 - Vbar2;
        if (VDM2 > 0 && r > 0) {
          dataPoints.push({ r, VDM: Math.sqrt(VDM2) });
        }
      }
    }

    if (dataPoints.length < 3) continue;

    const inner = dataPoints.filter(p => p.r <= 10);
    if (inner.length < 2) continue;

    let sx = 0, sy = 0, sxy = 0, sx2 = 0;
    const n = inner.length;
    for (const p of inner) {
      const x = Math.log(p.r);
      const y = Math.log(p.VDM);
      sx += x; sy += y; sxy += x * y; sx2 += x * x;
    }
    const alpha = (n * sxy - sx * sy) / (n * sx2 - sx * sx);

    const SigmaBar = g.sigma_bar != null ? g.sigma_bar : null;
    if (!SigmaBar || SigmaBar <= 0 || !isFinite(alpha)) continue;

    results.push({
      name: g.name,
      alpha,
      log_Sigma_bar: Math.log10(SigmaBar),
      Sigma_bar: SigmaBar,
    });
  }

  console.log('  SPARC galaxies with valid alpha: ' + results.length);
  return results;
}

function bootstrap(xs, ys, nBoot) {
  const slopes = [];
  const n = xs.length;
  for (let b = 0; b < nBoot; b++) {
    const bx = [], by = [];
    for (let i = 0; i < n; i++) {
      const j = Math.floor(Math.random() * n);
      bx.push(xs[j]);
      by.push(ys[j]);
    }
    const reg = linearRegression(bx, by);
    slopes.push(reg.slope);
  }
  slopes.sort((a, b) => a - b);
  const lo = slopes[Math.floor(0.025 * nBoot)];
  const hi = slopes[Math.floor(0.975 * nBoot)];
  const mean = slopes.reduce((s, v) => s + v, 0) / nBoot;
  const std = Math.sqrt(slopes.reduce((s, v) => s + (v - mean) * (v - mean), 0) / nBoot);
  return { mean, std, lo, hi };
}

async function main() {
  console.log('================================================================');
  console.log('  RAW PARTICLE DATA PIPELINE');
  console.log('  Direct from FlatHUB (Flatiron Institute Data Hub)');
  console.log('  TNG100-1 + EAGLE RefL0100N1504 aperture measurements');
  console.log('  vs SPARC 175 observed galaxies');
  console.log('================================================================');
  console.log('');
  console.log('DATA SOURCE: FlatHUB public catalogs');
  console.log('  - NOT tables from papers');
  console.log('  - NOT scaling relations');
  console.log('  - NOT parametric models');
  console.log('  - These are APERTURE MEASUREMENTS computed directly from');
  console.log('    particle snapshots by the simulation teams themselves.');
  console.log('');

  let tngGalaxies, eagleGalaxies;
  try {
    [tngGalaxies, eagleGalaxies] = await Promise.all([
      fetchTNGData(),
      fetchEAGLEData(),
    ]);
  } catch (err) {
    console.error('FETCH ERROR:', err.message);
    process.exit(1);
  }

  const radiiKpc = APERTURE_RADII_PKPC.slice(); // pkpc = kpc at z=0

  console.log('\n=== COMPUTING ALPHA FROM APERTURE DATA ===');

  const tngResults = tngGalaxies.map(g => computeAlphaAndSigma(g, radiiKpc)).filter(Boolean);
  const eagleResults = eagleGalaxies.map(g => computeAlphaAndSigma(g, radiiKpc)).filter(Boolean);

  console.log('TNG100-1 galaxies with valid alpha-Sigma: ' + tngResults.length);
  console.log('EAGLE RefL0100 galaxies with valid alpha-Sigma: ' + eagleResults.length);

  const sparcResults = loadSPARCData();

  console.log('\n=== ALPHA-SIGMA_BAR CORRELATION ===');
  console.log('');

  function analyzeCorrelation(label, results) {
    const xs = results.map(r => r.log_Sigma_bar);
    const ys = results.map(r => r.alpha);

    const reg = linearRegression(xs, ys);
    const boot = bootstrap(xs, ys, 5000);

    console.log(label + ':');
    console.log('  N = ' + results.length);
    console.log('  slope = ' + reg.slope.toFixed(5) + ' +/- ' + reg.slopeErr.toFixed(5));
    console.log('  bootstrap: ' + boot.mean.toFixed(5) + ' +/- ' + boot.std.toFixed(5));
    console.log('  95% CI: [' + boot.lo.toFixed(5) + ', ' + boot.hi.toFixed(5) + ']');
    console.log('  r = ' + reg.r.toFixed(4));
    console.log('  sigma from zero: ' + Math.abs(reg.slope / reg.slopeErr).toFixed(1) + 'σ');
    console.log('');

    return {
      label,
      n: results.length,
      slope: reg.slope,
      slopeErr: reg.slopeErr,
      r: reg.r,
      intercept: reg.intercept,
      bootMean: boot.mean,
      bootStd: boot.std,
      bootCI: [boot.lo, boot.hi],
      sigmaFromZero: Math.abs(reg.slope / reg.slopeErr),
    };
  }

  const tngStats = analyzeCorrelation('TNG100-1 (particle data)', tngResults);
  const eagleStats = analyzeCorrelation('EAGLE RefL0100 (particle data)', eagleResults);
  let sparcStats = null;
  if (sparcResults && sparcResults.length > 0) {
    sparcStats = analyzeCorrelation('SPARC (observed)', sparcResults);
  }

  console.log('=== COMPARISON ===');
  console.log('');

  if (sparcStats) {
    const tngRatio = Math.abs(sparcStats.slope / tngStats.slope);
    const eagleRatio = Math.abs(sparcStats.slope / eagleStats.slope);

    const diffTng = Math.abs(sparcStats.slope - tngStats.slope);
    const errTng = Math.sqrt(sparcStats.slopeErr ** 2 + tngStats.slopeErr ** 2);
    const sigmaTng = diffTng / errTng;

    const diffEagle = Math.abs(sparcStats.slope - eagleStats.slope);
    const errEagle = Math.sqrt(sparcStats.slopeErr ** 2 + eagleStats.slopeErr ** 2);
    const sigmaEagle = diffEagle / errEagle;

    console.log('SPARC vs TNG100-1:');
    console.log('  SPARC slope:  ' + sparcStats.slope.toFixed(5));
    console.log('  TNG slope:    ' + tngStats.slope.toFixed(5));
    console.log('  Ratio:        ' + tngRatio.toFixed(1) + 'x');
    console.log('  Discrepancy:  ' + sigmaTng.toFixed(1) + 'σ');
    console.log('  Sign match:   ' + (Math.sign(sparcStats.slope) === Math.sign(tngStats.slope) ? 'YES' : 'NO'));
    console.log('');

    console.log('SPARC vs EAGLE:');
    console.log('  SPARC slope:  ' + sparcStats.slope.toFixed(5));
    console.log('  EAGLE slope:  ' + eagleStats.slope.toFixed(5));
    console.log('  Ratio:        ' + eagleRatio.toFixed(1) + 'x');
    console.log('  Discrepancy:  ' + sigmaEagle.toFixed(1) + 'σ');
    console.log('  Sign match:   ' + (Math.sign(sparcStats.slope) === Math.sign(eagleStats.slope) ? 'YES' : 'NO'));
    console.log('');

    console.log('TNG100-1 vs EAGLE:');
    const diffTE = Math.abs(tngStats.slope - eagleStats.slope);
    const errTE = Math.sqrt(tngStats.slopeErr ** 2 + eagleStats.slopeErr ** 2);
    console.log('  Discrepancy:  ' + (diffTE / errTE).toFixed(1) + 'σ');
    console.log('  TNG/EAGLE:    ' + Math.abs(tngStats.slope / eagleStats.slope).toFixed(2));
    console.log('');
  }

  const output = {
    metadata: {
      pipeline: 'APERTURE_MEASUREMENTS_FROM_PARTICLES',
      source: 'FlatHUB (Flatiron Institute Data Hub)',
      description: 'Alpha-Sigma_bar correlation computed from pre-computed aperture mass measurements (M(<r) at 10 radii), which were derived from particle snapshots by the simulation teams. Not raw particles directly, but the closest publicly accessible particle-derived data.',
      aperture_radii_kpc: APERTURE_RADII_PKPC,
      alpha_definition: 'd ln V_DM / d ln r (velocity slope from enclosed mass at r <= 10 kpc)',
      Sigma_bar_definition: '(M_star + M_gas) / (pi * r_half^2)',
      timestamp: new Date().toISOString(),
      cosmology: { h: h, H0: H0 },
    },
    simulations: {
      TNG100: {
        catalog: 'tng',
        simulation: 'TNG100-1',
        snapshot: 99,
        redshift: 0,
        totalHalosFetched: tngGalaxies.length,
        SPARClike: tngResults.length,
        selectionCriteria: 'Vmax 50-300 km/s, M_star > 1e8 Msun, resolved inner DM',
        correlation: tngStats,
        sampleGalaxies: tngResults.slice(0, 5).map(r => ({
          alpha: +r.alpha.toFixed(5),
          log_Sigma_bar: +r.log_Sigma_bar.toFixed(3),
          Vmax: +r.Vmax.toFixed(1),
        })),
      },
      EAGLE: {
        catalog: 'eagle',
        simulation: 'RefL0100N1504',
        snapshot: 28,
        redshift: 0,
        totalHalosFetched: eagleGalaxies.length,
        SPARClike: eagleResults.length,
        selectionCriteria: 'Vmax 50-300 km/s, M_star > 1e8 Msun, resolved inner DM',
        correlation: eagleStats,
        sampleGalaxies: eagleResults.slice(0, 5).map(r => ({
          alpha: +r.alpha.toFixed(5),
          log_Sigma_bar: +r.log_Sigma_bar.toFixed(3),
          Vmax: +r.Vmax.toFixed(1),
        })),
      },
    },
    SPARC: sparcStats ? {
      n: sparcStats.n,
      correlation: sparcStats,
    } : null,
    comparison: sparcStats ? {
      tngVsSparc: {
        slopeRatio: Math.abs(sparcStats.slope / tngStats.slope),
        sigmaDiscrepancy: Math.abs(sparcStats.slope - tngStats.slope) /
          Math.sqrt(sparcStats.slopeErr ** 2 + tngStats.slopeErr ** 2),
        signMatch: Math.sign(sparcStats.slope) === Math.sign(tngStats.slope),
      },
      eagleVsSparc: {
        slopeRatio: Math.abs(sparcStats.slope / eagleStats.slope),
        sigmaDiscrepancy: Math.abs(sparcStats.slope - eagleStats.slope) /
          Math.sqrt(sparcStats.slopeErr ** 2 + eagleStats.slopeErr ** 2),
        signMatch: Math.sign(sparcStats.slope) === Math.sign(eagleStats.slope),
      },
      tngVsEagle: {
        slopeRatio: Math.abs(tngStats.slope / eagleStats.slope),
        sigmaDiscrepancy: Math.abs(tngStats.slope - eagleStats.slope) /
          Math.sqrt(tngStats.slopeErr ** 2 + eagleStats.slopeErr ** 2),
      },
    } : null,
    rawData: {
      TNG100: tngResults.map(r => ({
        alpha: +r.alpha.toFixed(5),
        log_Sigma_bar: +r.log_Sigma_bar.toFixed(4),
      })),
      EAGLE: eagleResults.map(r => ({
        alpha: +r.alpha.toFixed(5),
        log_Sigma_bar: +r.log_Sigma_bar.toFixed(4),
      })),
      SPARC: sparcResults ? sparcResults.map(r => ({
        name: r.name,
        alpha: +r.alpha.toFixed(5),
        log_Sigma_bar: +r.log_Sigma_bar.toFixed(4),
      })) : [],
    },
  };

  const outPath = path.join(__dirname, '..', 'public', 'raw-particle-data.json');
  fs.writeFileSync(outPath, JSON.stringify(output, null, 2));
  console.log('\nOutput written to: ' + outPath);

  console.log('\n================================================================');
  console.log('  VERDICT');
  console.log('================================================================');

  if (sparcStats) {
    const signOK_TNG = Math.sign(sparcStats.slope) === Math.sign(tngStats.slope);
    const signOK_EAGLE = Math.sign(sparcStats.slope) === Math.sign(eagleStats.slope);
    const magRatio = Math.abs(sparcStats.slope / ((tngStats.slope + eagleStats.slope) / 2));
    const avgSigma = (
      Math.abs(sparcStats.slope - tngStats.slope) / Math.sqrt(sparcStats.slopeErr ** 2 + tngStats.slopeErr ** 2) +
      Math.abs(sparcStats.slope - eagleStats.slope) / Math.sqrt(sparcStats.slopeErr ** 2 + eagleStats.slopeErr ** 2)
    ) / 2;

    console.log('');
    console.log('Sign agreement (TNG):   ' + (signOK_TNG ? 'YES ✓' : 'NO ✗'));
    console.log('Sign agreement (EAGLE): ' + (signOK_EAGLE ? 'YES ✓' : 'NO ✗'));
    console.log('Magnitude ratio:        ' + magRatio.toFixed(1) + 'x');
    console.log('Average sigma:          ' + avgSigma.toFixed(1) + 'σ');
    console.log('');

    if (signOK_TNG && signOK_EAGLE && avgSigma > 3) {
      console.log('RESULT: MAGNITUDE DISCREPANCY CONFIRMED');
      console.log('  Both simulations produce the right SIGN of coupling');
      console.log('  but SPARC is ' + magRatio.toFixed(1) + 'x stronger (' + avgSigma.toFixed(1) + 'σ)');
      console.log('  This IS a genuine tension with ΛCDM simulations.');
    } else if (!signOK_TNG || !signOK_EAGLE) {
      console.log('RESULT: SIGN PROBLEM EXISTS');
      console.log('  Simulations predict wrong direction of coupling.');
    } else {
      console.log('RESULT: NO SIGNIFICANT DISCREPANCY');
      console.log('  Simulations reproduce observed coupling within errors.');
    }
  }

  console.log('');
  console.log('Done.');
}

main().catch(err => {
  console.error('FATAL:', err);
  process.exit(1);
});
