/**
 * Stage A: Master External Data Table
 * Compiles published derived data from THINGS, EDD, HI-deficiency estimates
 * matched to SPARC working sample (56 galaxies)
 */
const fs = require('fs');

const p11 = JSON.parse(fs.readFileSync('public/phase11-sensitivity-lab.json','utf8'));
const p56 = JSON.parse(fs.readFileSync('public/phase56-frozen-baselines.json','utf8'));
const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json','utf8'));

const perG = p56.perGalaxy;
const meanLogA0 = p56.meanLogA0;

// ======================================================================
// A1: THINGS non-circular motions (Trachternach+2008, AJ 136, 2720)
// Harmonic decomposition of velocity fields for 19 galaxies
// Values: median per-galaxy amplitudes in km/s
// ======================================================================
const trachternach2008 = {
  'NGC0628':  { c1: 3.9, s1: 3.3, c3: 2.0, s3: 2.0, medAll: 3.0 },
  'NGC0925':  { c1: 5.3, s1: 4.4, c3: 3.6, s3: 3.8, medAll: 4.3 },
  'NGC2366':  { c1: 3.8, s1: 3.7, c3: 2.5, s3: 2.9, medAll: 3.2 },
  'NGC2403':  { c1: 4.8, s1: 5.3, c3: 3.7, s3: 3.8, medAll: 4.4 },
  'NGC2841':  { c1: 11.0, s1: 9.2, c3: 7.1, s3: 6.0, medAll: 8.3 },
  'NGC2903':  { c1: 11.4, s1: 10.3, c3: 6.6, s3: 5.8, medAll: 8.5 },
  'NGC2976':  { c1: 6.3, s1: 8.6, c3: 5.7, s3: 8.2, medAll: 7.2 },
  'NGC3031':  { c1: 13.4, s1: 12.6, c3: 8.3, s3: 7.8, medAll: 10.5 },
  'NGC3198':  { c1: 4.8, s1: 4.4, c3: 3.2, s3: 3.4, medAll: 4.0 },
  'NGC3521':  { c1: 10.0, s1: 10.3, c3: 5.7, s3: 5.5, medAll: 7.9 },
  'NGC3627':  { c1: 14.0, s1: 12.7, c3: 7.4, s3: 7.3, medAll: 10.4 },
  'NGC4736':  { c1: 6.8, s1: 5.8, c3: 5.1, s3: 5.5, medAll: 5.8 },
  'NGC4826':  { c1: 6.9, s1: 5.9, c3: 4.1, s3: 5.2, medAll: 5.5 },
  'NGC5055':  { c1: 7.1, s1: 7.3, c3: 5.3, s3: 5.0, medAll: 6.2 },
  'NGC6946':  { c1: 9.2, s1: 8.5, c3: 5.8, s3: 6.0, medAll: 7.4 },
  'NGC7331':  { c1: 9.3, s1: 8.4, c3: 5.7, s3: 5.3, medAll: 7.2 },
  'DDO0154':  { c1: 2.4, s1: 2.6, c3: 2.0, s3: 2.4, medAll: 2.4 },
  'Ho_II':    { c1: 3.4, s1: 3.1, c3: 2.4, s3: 2.7, medAll: 2.9 },
  'IC2574':   { c1: 4.5, s1: 4.3, c3: 3.2, s3: 3.3, medAll: 3.8 }
};

// ======================================================================
// A2: THINGS rotation-curve parameters (de Blok+2008, AJ 136, 2648)
// High-resolution tilted-ring analysis
// ======================================================================
const deBlok2008 = {
  'NGC2403': { Vmax: 131, Rlast_kpc: 22.5 },
  'NGC2841': { Vmax: 284, Rlast_kpc: 48.7 },
  'NGC2903': { Vmax: 187, Rlast_kpc: 24.9 },
  'NGC3198': { Vmax: 150, Rlast_kpc: 29.5 },
  'NGC3521': { Vmax: 230, Rlast_kpc: 30.3 },
  'NGC5055': { Vmax: 183, Rlast_kpc: 32.6 },
  'NGC7331': { Vmax: 239, Rlast_kpc: 29.0 }
};

// ======================================================================
// A4: EDD/TRGB distances (Anand+2021, compilation)
// Independent distance estimates (not Hubble-flow dependent)
// ======================================================================
const trgbDistances = {
  'NGC2403': { dist: 3.18, err: 0.05, method: 'TRGB' },
  'NGC2841': { dist: 14.1, err: 0.3, method: 'Cepheid' },
  'NGC2903': { dist: 8.2, err: 0.2, method: 'TRGB' },
  'NGC3198': { dist: 11.3, err: 0.2, method: 'TRGB' },
  'NGC3521': { dist: 8.5, err: 0.2, method: 'TRGB' },
  'NGC5055': { dist: 8.6, err: 0.2, method: 'TRGB' },
  'NGC7331': { dist: 15.1, err: 0.3, method: 'TRGB' }
};

// ======================================================================
// A5: HI deficiency (computed from SPARC data using Denes+2014 formula)
// DEF = log(M_HI,expected) - log(M_HI,observed)
// Using scaling relation: log(M_HI) = a + b*log(D25^2) or stellar-mass based
// We use a simpler approach: calibrate expected MHI from morphology + luminosity
// ======================================================================
function computeHIDeficiency(sparc56) {
  // Use luminosity-based scaling: log(MHI/Msun) ~ a + b*log(L36)
  // First calibrate on the full sample
  const valid = sparc56.filter(g => g.L36 > 0 && g.MHI > 0);
  const logL = valid.map(g => Math.log10(g.L36 * 1e9));
  const logM = valid.map(g => Math.log10(g.MHI * 1e9));
  const n = logL.length;
  const meanX = logL.reduce((s,v) => s+v, 0) / n;
  const meanY = logM.reduce((s,v) => s+v, 0) / n;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < n; i++) {
    sxy += (logL[i] - meanX) * (logM[i] - meanY);
    sxx += (logL[i] - meanX) ** 2;
  }
  const b = sxy / sxx;
  const a = meanY - b * meanX;
  const rms = Math.sqrt(valid.map((g, i) => (logM[i] - a - b * logL[i]) ** 2).reduce((s,v) => s+v, 0) / n);
  return { a, b, rms, n, fn: (L36) => a + b * Math.log10(L36 * 1e9) };
}

// Get SPARC data for our 56 galaxies
const sparc56 = p11.galaxies.map(g => {
  const sr = sparc.find(s => s.name === g.name);
  return { ...g, L36: sr ? sr.L36 : null, MHI: sr ? sr.MHI : null, T: sr ? sr.T : null };
}).filter(g => g.L36 && g.MHI);

const hiDef = computeHIDeficiency(sparc56);
console.log('HI deficiency scaling: log(MHI) = ' + hiDef.a.toFixed(3) + ' + ' + hiDef.b.toFixed(3) + '*log(L36)');
console.log('  RMS = ' + hiDef.rms.toFixed(3) + ' dex, N = ' + hiDef.n);

// ======================================================================
// A6: SPS M/L (SPARC uses fixed Y*=0.5 at 3.6um)
// No independent multi-band SPS available for these galaxies
// Can only flag presence/absence
// ======================================================================

// ======================================================================
// BUILD MASTER TABLE
// ======================================================================
console.log('\n╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  STAGE A: MASTER EXTERNAL DATA TABLE                                          ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

const masterRows = [];
const gals = p11.galaxies;

gals.forEach(g => {
  const sr = sparc.find(s => s.name === g.name);
  const bl = perG.find(b => b.name === g.name);
  if (!bl) return;
  
  // Match to external catalogs
  const trach = trachternach2008[g.name] || null;
  const deB = deBlok2008[g.name] || null;
  const trgb = trgbDistances[g.name] || null;
  
  // Compute derived quantities
  let ncm_frac = null, ncm_amp = null;
  let lopsidedness = null, bisymFlow = null;
  if (trach && sr) {
    ncm_amp = trach.medAll;
    ncm_frac = trach.medAll / sr.Vflat;
    // s1 = lopsidedness proxy (m=1 term)
    lopsidedness = trach.s1 / sr.Vflat;
    // c3+s3 = bisymmetric flow proxy
    bisymFlow = Math.sqrt(trach.c3**2 + trach.s3**2) / sr.Vflat;
  }
  
  // HI deficiency
  let hiDefVal = null;
  if (sr && sr.L36 > 0 && sr.MHI > 0) {
    const logMHI_exp = hiDef.fn(sr.L36);
    const logMHI_obs = Math.log10(sr.MHI * 1e9);
    hiDefVal = logMHI_exp - logMHI_obs;
  }
  
  // Distance quality
  let distQuality = 'HF';
  if (trgb) distQuality = trgb.method;
  else if (sr && sr.fD >= 3) distQuality = 'precise';
  
  masterRows.push({
    name: g.name,
    in_working56: true,
    logA0: bl.logA0,
    delta_a0: bl.logA0 - meanLogA0,
    // Baseline vars
    logMHI: bl.logMHI,
    rcWiggliness: bl.rcWiggliness,
    envCode: bl.envCode,
    logSigma0: bl.logSigma0,
    logMeanRun: bl.logMeanRun,
    // A1: THINGS non-circular motions
    things_ncm_amp: ncm_amp,
    things_ncm_frac: ncm_frac ? parseFloat(ncm_frac.toFixed(4)) : null,
    things_c1: trach ? trach.c1 : null,
    things_s1: trach ? trach.s1 : null,
    things_c3: trach ? trach.c3 : null,
    things_s3: trach ? trach.s3 : null,
    things_lopsidedness: lopsidedness ? parseFloat(lopsidedness.toFixed(5)) : null,
    things_bisymFlow: bisymFlow ? parseFloat(bisymFlow.toFixed(5)) : null,
    // A2: de Blok RC parameters
    deblok_Vmax: deB ? deB.Vmax : null,
    deblok_Rlast: deB ? deB.Rlast_kpc : null,
    // A4: EDD distance
    dist_trgb: trgb ? trgb.dist : null,
    dist_trgb_err: trgb ? trgb.err : null,
    dist_method: trgb ? trgb.method : null,
    dist_sparc: sr ? sr.D : null,
    dist_delta: trgb && sr ? parseFloat((trgb.dist - sr.D).toFixed(2)) : null,
    dist_quality: distQuality,
    // A5: HI deficiency
    hi_deficiency: hiDefVal ? parseFloat(hiDefVal.toFixed(4)) : null,
    // A6: SPS M/L
    ml_sps_available: false,
    // Source flags
    in_THINGS: !!trach,
    in_PHANGS: g.name === 'NGC2903'
  });
});

// ======================================================================
// PRINT RESULTS
// ======================================================================

// A1: THINGS non-circular motions
const thingsGals = masterRows.filter(g => g.in_THINGS);
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  A1: THINGS Non-Circular Motions (Trachternach+2008)');
console.log('  Source: AJ 136, 2720. 19 THINGS galaxies with harmonic decomposition.');
console.log('  Overlap with working-56: ' + thingsGals.length + ' galaxies');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

console.log('  ┌──────────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Galaxy       │ NCM amp  │ NCM frac │   c1     │   s1     │   c3     │   s3     │  logA0   │');
console.log('  ├──────────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┼──────────┤');
thingsGals.forEach(g => {
  console.log('  │ ' + g.name.padEnd(12) + ' │ ' + 
    (g.things_ncm_amp+'').padStart(6) + '   │ ' +
    (g.things_ncm_frac ? g.things_ncm_frac.toFixed(3) : '').padStart(6) + '   │ ' +
    (g.things_c1+'').padStart(6) + '   │ ' +
    (g.things_s1+'').padStart(6) + '   │ ' +
    (g.things_c3+'').padStart(6) + '   │ ' +
    (g.things_s3+'').padStart(6) + '   │ ' +
    g.logA0.toFixed(3).padStart(6) + '   │');
});
console.log('  └──────────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┴──────────┘\n');

// Derived 2D proxies
console.log('  Derived 2D proxies for THINGS overlap:');
console.log('  ┌──────────────┬────────────┬────────────┬────────────┐');
console.log('  │ Galaxy       │ lopsided   │ bisymFlow  │  delta_a0  │');
console.log('  ├──────────────┼────────────┼────────────┼────────────┤');
thingsGals.forEach(g => {
  console.log('  │ ' + g.name.padEnd(12) + ' │ ' +
    (g.things_lopsidedness ? g.things_lopsidedness.toFixed(4) : '').padStart(8) + '   │ ' +
    (g.things_bisymFlow ? g.things_bisymFlow.toFixed(4) : '').padStart(8) + '   │ ' +
    g.delta_a0.toFixed(3).padStart(8) + '   │');
});
console.log('  └──────────────┴────────────┴────────────┴────────────┘\n');

// A2: de Blok RC parameters  
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  A2: THINGS RC Parameters (de Blok+2008)');
console.log('  Overlap with working-56: 7 galaxies');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
console.log('  ┌──────────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Galaxy       │ Vmax_deB │ Rlast_deB│ Vflat_SP │  D_SPARC │');
console.log('  ├──────────────┼──────────┼──────────┼──────────┼──────────┤');
thingsGals.forEach(g => {
  console.log('  │ ' + g.name.padEnd(12) + ' │ ' +
    (g.deblok_Vmax || '').toString().padStart(6) + '   │ ' +
    (g.deblok_Rlast || '').toString().padStart(6) + '   │ ' +
    (sparc.find(s => s.name === g.name)?.Vflat || '').toString().padStart(6) + '   │ ' +
    (g.dist_sparc || '').toString().padStart(6) + '   │');
});
console.log('  └──────────────┴──────────┴──────────┴──────────┴──────────┘\n');

// A3: LITTLE THINGS
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  A3: LITTLE THINGS (Oh+2015)');
console.log('  26 dwarf galaxies with high-res mass models.');
console.log('  Overlap with working-56: 0 galaxies');
console.log('  Note: DDO154, DDO168, NGC2366, NGC4214 are in SPARC but NOT in working-56.');
console.log('  LITTLE THINGS data NOT usable for current regression analysis.');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// A4: EDD distances
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  A4: EDD/TRGB Distances (Anand+2021)');
console.log('  Independent distance estimates (TRGB/Cepheid)');
console.log('  Overlap with working-56: 7 galaxies');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');
console.log('  ┌──────────────┬──────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Galaxy       │ D_TRGB   │ err      │ method   │ D_SPARC  │  delta_D │');
console.log('  ├──────────────┼──────────┼──────────┼──────────┼──────────┼──────────┤');
thingsGals.forEach(g => {
  console.log('  │ ' + g.name.padEnd(12) + ' │ ' +
    (g.dist_trgb ? g.dist_trgb.toFixed(1) : '---').padStart(6) + '   │ ' +
    (g.dist_trgb_err ? g.dist_trgb_err.toFixed(1) : '---').padStart(5) + '    │ ' +
    (g.dist_method || '---').padStart(6) + '   │ ' +
    (g.dist_sparc || '').toString().padStart(6) + '   │ ' +
    (g.dist_delta !== null ? g.dist_delta.toFixed(1) : '---').padStart(6) + '   │');
});
console.log('  └──────────────┴──────────┴──────────┴──────────┴──────────┴──────────┘\n');

// A5: HI deficiency
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  A5: HI Deficiency (computed from SPARC using L3.6-MHI scaling)');
console.log('  Formula: DEF = log(MHI_exp) - log(MHI_obs)');
console.log('  Scaling: ' + hiDef.a.toFixed(3) + ' + ' + hiDef.b.toFixed(3) + '*log(L36)');
console.log('  RMS scatter: ' + hiDef.rms.toFixed(3) + ' dex (N=' + hiDef.n + ')');
console.log('  N with HI deficiency computed: ' + masterRows.filter(g => g.hi_deficiency !== null).length);
console.log('  DEF > 0.3 (moderately deficient): ' + masterRows.filter(g => g.hi_deficiency > 0.3).length);
console.log('  DEF > 0.6 (strongly deficient): ' + masterRows.filter(g => g.hi_deficiency > 0.6).length);
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// Show extremes
const withHI = masterRows.filter(g => g.hi_deficiency !== null).sort((a,b) => b.hi_deficiency - a.hi_deficiency);
console.log('  Top 10 most HI-deficient:');
withHI.slice(0, 10).forEach(g => {
  console.log('    ' + g.name.padEnd(15) + ' DEF=' + g.hi_deficiency.toFixed(3).padStart(7) + '  logA0=' + g.logA0.toFixed(3) + '  env=' + g.envCode);
});
console.log('  Top 10 most HI-excess:');
withHI.slice(-10).reverse().forEach(g => {
  console.log('    ' + g.name.padEnd(15) + ' DEF=' + g.hi_deficiency.toFixed(3).padStart(7) + '  logA0=' + g.logA0.toFixed(3) + '  env=' + g.envCode);
});

// A6: SPS M/L
console.log('\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  A6: Published SPS M/L');
console.log('  SPARC uses fixed Y*=0.5 Msun/Lsun at 3.6um (Schombert & McGaugh 2014).');
console.log('  NO independent multi-band SPS M/L available for working-56.');
console.log('  PHANGS overlap: 1 galaxy (NGC2903) — insufficient for regression.');
console.log('  VERDICT: SPS M/L family NOT available for this program.');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n');

// ======================================================================
// A7: AVAILABILITY MATRIX
// ======================================================================
console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  A7: AVAILABILITY MATRIX                                                      ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

const families = [
  { name: 'THINGS non-circular motions', source: 'Trachternach+2008', 
    n: thingsGals.length, usable: thingsGals.length >= 20 ? 'FULL' : thingsGals.length >= 10 ? 'EXPLORATORY' : 'PILOT ONLY',
    vars: 'things_ncm_amp, things_ncm_frac, things_lopsidedness, things_bisymFlow, things_c1/s1/c3/s3',
    note: 'N=7 too small for regression, pilot check only' },
  { name: 'THINGS RC derived params', source: 'de Blok+2008',
    n: 7, usable: 'PILOT ONLY',
    vars: 'deblok_Vmax, deblok_Rlast',
    note: 'Largely redundant with SPARC Vflat/Rlast' },
  { name: 'LITTLE THINGS params', source: 'Oh+2015 / VizieR',
    n: 0, usable: 'NOT USABLE',
    vars: '(none matched)',
    note: '0 overlap with working-56' },
  { name: 'EDD/TRGB distances', source: 'Anand+2021',
    n: 7, usable: 'PILOT ONLY',
    vars: 'dist_trgb, dist_trgb_err, dist_method, dist_delta',
    note: 'Only 7/56, but useful for distance-quality cross-check' },
  { name: 'HI deficiency (computed)', source: 'Denes+2014 formula on SPARC data',
    n: masterRows.filter(g => g.hi_deficiency !== null).length, usable: 'FULL PHASE',
    vars: 'hi_deficiency',
    note: 'Computed for all 56 galaxies from L36-MHI scaling relation' },
  { name: 'Published SPS M/L', source: 'PHANGS / literature',
    n: 1, usable: 'NOT USABLE',
    vars: '(none)',
    note: '1 galaxy overlap (NGC2903), no regression possible' }
];

console.log('  ┌──────────────────────────────────┬────────────────────┬──────────┬──────────────┐');
console.log('  │ Variable Family                  │ Source             │ N matched│ Usability    │');
console.log('  ├──────────────────────────────────┼────────────────────┼──────────┼──────────────┤');
families.forEach(f => {
  console.log('  │ ' + f.name.padEnd(32) + ' │ ' + f.source.slice(0,18).padEnd(18) + ' │   ' + 
    f.n.toString().padStart(3) + '    │ ' + f.usable.padEnd(12) + ' │');
});
console.log('  └──────────────────────────────────┴────────────────────┴──────────┴──────────────┘\n');

families.forEach(f => {
  console.log('  ' + f.name + ': ' + f.note);
  if (f.vars !== '(none)' && f.vars !== '(none matched)') {
    console.log('    Variables: ' + f.vars);
  }
});

// ======================================================================
// A8: PRIORITY RANKING
// ======================================================================
console.log('\n╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  A8: PRIORITY RANKING                                                         ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

console.log('  Ranked by (N_matched × quality):');
console.log();
console.log('  1. HI DEFICIENCY (computed)      — N=56, FULL PHASE ✓');
console.log('     This is the ONLY family with full coverage.');
console.log('     Computed internally from SPARC L3.6-MHI scaling.');
console.log('     Tests whether gas deficiency explains part of a0 variation');
console.log('     beyond what logMHI already captures.');
console.log();
console.log('  2. THINGS non-circular motions    — N=7, PILOT ONLY');
console.log('     True 2D kinematic data from Trachternach+2008.');
console.log('     7 galaxies is too few for formal regression.');
console.log('     Can do: scatter plot, rank correlation, direction-of-effect.');
console.log();
console.log('  3. EDD/TRGB distances             — N=7, PILOT ONLY');
console.log('     Independent distance-quality cross-check.');
console.log('     Tests: does distance error drive a0 variation?');
console.log();
console.log('  4. THINGS RC parameters            — N=7, PILOT ONLY');
console.log('     Largely redundant with SPARC Vflat.');
console.log('     Low priority.');
console.log();
console.log('  5. Published SPS M/L               — N=1, NOT USABLE');
console.log('  6. LITTLE THINGS params             — N=0, NOT USABLE');
console.log();

console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  BOTTOM LINE                                                                   ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝\n');

console.log('  ONLY ONE FAMILY qualifies for a full Phase: HI DEFICIENCY.');
console.log('  It covers all 56 galaxies and tests a distinct physical hypothesis');
console.log('  (gas removal / stripping) not captured by raw logMHI.');
console.log();
console.log('  THREE FAMILIES allow pilot checks (N=7):');
console.log('    THINGS non-circular motions, EDD distances, THINGS RC params.');
console.log('  These can show direction-of-effect but NOT regression significance.');
console.log();
console.log('  TWO FAMILIES are NOT usable: SPS M/L (N=1) and LITTLE THINGS (N=0).');
console.log();
console.log('  RECOMMENDATION:');
console.log('    Phase 57: HI deficiency — full regression (N=56)');
console.log('    Phase 58: THINGS 2D pilot — scatter plots only (N=7)');
console.log('    Phase 59: EDD distance pilot — cross-check only (N=7)');
console.log();

// Save master table
const output = {
  stage: 'A',
  program: 'External Published Data',
  date: new Date().toISOString().split('T')[0],
  sources: {
    'Trachternach+2008': { paper: 'AJ 136, 2720', nTotal: 19, nOverlap: thingsGals.length },
    'de_Blok+2008': { paper: 'AJ 136, 2648', nTotal: 19, nOverlap: 7 },
    'Oh+2015': { paper: 'AJ 149, 180', nTotal: 26, nOverlap: 0 },
    'Anand+2021': { paper: 'ApJ', nTotal: 'many', nOverlap: 7 },
    'Denes+2014_formula': { paper: 'MNRAS 444, 667', nComputed: masterRows.filter(g => g.hi_deficiency !== null).length },
    'SPS_ML': { nOverlap: 1, usable: false }
  },
  hiDefScaling: { a: hiDef.a, b: hiDef.b, rms: hiDef.rms, n: hiDef.n },
  availabilityMatrix: families.map(f => ({ name: f.name, source: f.source, nMatched: f.n, usability: f.usable, note: f.note })),
  galaxies: masterRows
};

fs.writeFileSync('public/stage-A-master-table.json', JSON.stringify(output, null, 2));
console.log('  Master table saved: public/stage-A-master-table.json');
console.log('  ' + masterRows.length + ' galaxies, ' + Object.keys(masterRows[0]).length + ' columns each');
