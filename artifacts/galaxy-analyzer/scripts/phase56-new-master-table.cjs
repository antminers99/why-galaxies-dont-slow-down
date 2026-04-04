/**
 * Phase 56.0 (New Program): Master Table / Data Freeze
 * Build availability matrix for external survey data
 */
const fs = require('fs');

const sparc = JSON.parse(fs.readFileSync('public/sparc-table.json', 'utf8'));
const p11 = JSON.parse(fs.readFileSync('public/phase11-sensitivity-lab.json', 'utf8'));
const p56 = JSON.parse(fs.readFileSync('public/phase56-frozen-baselines.json', 'utf8'));

const workingGals = p11.galaxies;
const workingNames = workingGals.map(g => g.name);

// Known THINGS galaxies (Walter et al. 2008, de Blok et al. 2008)
// 34 galaxies with VLA HI 21cm, ~6" angular resolution
const thingsGalaxies = [
  'NGC0628','NGC0925','NGC2403','NGC2841','NGC2903','NGC2976',
  'NGC3031','NGC3184','NGC3198','NGC3351','NGC3521','NGC3621',
  'NGC3627','NGC4736','NGC4826','NGC5055','NGC5194','NGC6946',
  'NGC7331','NGC7793','DDO0154','IC2574','NGC2366',
  'NGC4214','NGC4449','NGC4631','NGC3077','NGC1569',
  'Ho_I','Ho_II','M81dwB','DDO053','DDO165','IC10'
];

// Known LITTLE THINGS galaxies (Oh et al. 2015, Hunter et al. 2012)
// 41 dwarf irregulars with VLA HI, ~6" resolution
const littleThingsGalaxies = [
  'CVnIdwA','DDO046','DDO047','DDO050','DDO052','DDO053',
  'DDO069','DDO070','DDO075','DDO087','DDO101','DDO126',
  'DDO133','DDO154','DDO155','DDO165','DDO167','DDO168',
  'DDO187','DDO210','DDO216','Haro29','Haro36','IC10',
  'IC1613','M81dwB','NGC1569','NGC2366','NGC3738','NGC4163',
  'NGC4214','NGC6822','SagDIG','UGC8508','WLM',
  'F564-V3','NGC4068','UGC7577','UGC7603','UGC8833','VII_Zw_403'
];

// PHANGS galaxies (Lee et al. 2023, Lang et al. 2020)
// ALMA CO(2-1) + VLT/MUSE IFU for ~90 galaxies
const phangsGalaxies = [
  'NGC0628','NGC1087','NGC1300','NGC1365','NGC1385','NGC1433',
  'NGC1512','NGC1566','NGC1672','NGC2835','NGC3351','NGC3627',
  'NGC4254','NGC4303','NGC4321','NGC4535','NGC4548','NGC4571',
  'NGC5068','NGC5248','NGC6744','NGC7496',
  'NGC1317','NGC1792','NGC2775','NGC2903','NGC3507','NGC3511',
  'NGC4424','NGC4457','NGC4569','NGC4579','NGC4689','NGC4694',
  'NGC4826','NGC5134','NGC5530','NGC6300'
];

console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  PHASE 56.0 (NEW PROGRAM): MASTER TABLE / DATA FREEZE                         ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝');
console.log();

// --- Match galaxies ---
function matchName(name, list) {
  if (list.includes(name)) return true;
  const n = name.replace(/\s+/g,'');
  if (list.some(l => l.replace(/\s+/g,'') === n)) return true;
  const nAlt = name.replace('NGC','NGC0');
  if (list.includes(nAlt)) return true;
  const nAlt2 = name.replace('NGC0','NGC');
  if (list.includes(nAlt2)) return true;
  if (name === 'DDO154' && list.includes('DDO0154')) return true;
  if (name === 'DDO0154' && list.includes('DDO154')) return true;
  return false;
}

// Build master table for each working galaxy
const masterTable = workingGals.map(g => {
  const sRow = sparc.find(s => s.name === g.name);
  const inThings = matchName(g.name, thingsGalaxies);
  const inLT = matchName(g.name, littleThingsGalaxies);
  const inPhangs = matchName(g.name, phangsGalaxies);
  
  let source = 'SPARC-only';
  if (inThings && inLT) source = 'THINGS+LT';
  else if (inThings) source = 'THINGS';
  else if (inLT) source = 'LITTLE_THINGS';
  if (inPhangs) source += '+PHANGS';
  
  // THINGS provides: moment-0 (HI), moment-1 (velocity), moment-2 (dispersion)
  // tilted-ring analysis, harmonic decomposition (Trachternach+2008)
  const has2DVelField = inThings || inLT; // These surveys have velocity fields
  const hasMoment1 = has2DVelField;
  const hasResidualMap = inThings; // Trachternach+2008 did harmonic decomp for THINGS
  const hasHIDeficiency = inThings || inLT; // Can compute from HI maps
  const hasMLInfo = inPhangs; // PHANGS has multi-band SPS M/L
  
  // Quality flags from SPARC
  const Q = sRow ? sRow.Q : null;
  const inc = sRow ? sRow.inc : null;
  const fD = sRow ? sRow.fD : null;
  const T = sRow ? sRow.T : null;
  
  // Baseline variables (from phase 56)
  const bl = p56.galaxyData ? p56.galaxyData.find(d => d.name === g.name) : null;
  
  return {
    name: g.name,
    inSPARC: true,
    inThings: inThings,
    inLittleThings: inLT,
    inPhangs: inPhangs,
    source: source,
    has2DVelField: has2DVelField,
    hasMoment1: hasMoment1,
    hasResidualMap: hasResidualMap,
    hasHIDeficiency: hasHIDeficiency,
    hasMLInfo: hasMLInfo,
    Q: Q,
    inc: inc,
    fD: fD,
    T: T,
    logMHI: bl ? bl.logMHI : null,
    rcWiggliness: bl ? bl.wig : null,
    envCode: bl ? bl.envCode : null,
    logSigma0_bar: bl ? bl.logSigma0 : null,
    meanRun: bl ? bl.meanRun : null,
    logA0: g.logA0
  };
});

// --- Summaries ---
const nTotal = masterTable.length;
const nThings = masterTable.filter(g => g.inThings).length;
const nLT = masterTable.filter(g => g.inLittleThings).length;
const nPhangs = masterTable.filter(g => g.inPhangs).length;
const nOverlap = masterTable.filter(g => g.inThings || g.inLittleThings || g.inPhangs).length;
const n2D = masterTable.filter(g => g.has2DVelField).length;
const nResid = masterTable.filter(g => g.hasResidualMap).length;
const nML = masterTable.filter(g => g.hasMLInfo).length;
const nHIDef = masterTable.filter(g => g.hasHIDeficiency).length;

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  HEADLINE COUNTS');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log();
console.log('  Total working sample:           ' + nTotal + ' galaxies');
console.log('  Overlap with THINGS:            ' + nThings);
console.log('  Overlap with LITTLE THINGS:     ' + nLT);
console.log('  Overlap with PHANGS:            ' + nPhangs);
console.log('  ANY external survey overlap:    ' + nOverlap);
console.log('  SPARC-only (no external):       ' + (nTotal - nOverlap));
console.log();

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  AVAILABILITY MATRIX');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log();
console.log('  ┌────────────────────┬──────────┬──────────┬──────────┬──────────┬──────────┐');
console.log('  │ Data Product       │  N avail │  % of 56 │  Source  │ In proj? │  Usable? │');
console.log('  ├────────────────────┼──────────┼──────────┼──────────┼──────────┼──────────┤');
console.log('  │ 2D velocity field  │   ' + String(n2D).padStart(5) + '  │  ' + (100*n2D/56).toFixed(1).padStart(5) + '% │ THINGS   │    NO    │    NO    │');
console.log('  │ Moment-1 map       │   ' + String(n2D).padStart(5) + '  │  ' + (100*n2D/56).toFixed(1).padStart(5) + '% │ THINGS   │    NO    │    NO    │');
console.log('  │ Residual map       │   ' + String(nResid).padStart(5) + '  │  ' + (100*nResid/56).toFixed(1).padStart(5) + '% │ THINGS   │    NO    │    NO    │');
console.log('  │ HI deficiency      │   ' + String(nHIDef).padStart(5) + '  │  ' + (100*nHIDef/56).toFixed(1).padStart(5) + '% │ THINGS   │    NO    │    NO    │');
console.log('  │ SPS M/L (multi-λ)  │   ' + String(nML).padStart(5) + '  │  ' + (100*nML/56).toFixed(1).padStart(5) + '% │ PHANGS   │    NO    │    NO    │');
console.log('  │ 1D rotation curve  │      56  │  100.0% │ SPARC    │   YES    │   YES    │');
console.log('  │ Baseline vars      │      56  │  100.0% │ Computed │   YES    │   YES    │');
console.log('  └────────────────────┴──────────┴──────────┴──────────┴──────────┴──────────┘');
console.log();

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  OVERLAP GALAXIES (detail)');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log();
console.log('  ┌──────────────┬────────┬────────┬────────┬──────┬──────┬──────┬────────────┐');
console.log('  │ Galaxy       │ THINGS │ LT     │ PHANGS │  Q   │  inc │  fD  │  logA0     │');
console.log('  ├──────────────┼────────┼────────┼────────┼──────┼──────┼──────┼────────────┤');
const overlapGals = masterTable.filter(g => g.inThings || g.inLittleThings || g.inPhangs);
overlapGals.forEach(g => {
  const t = g.inThings ? '  YES ' : '  --- ';
  const lt = g.inLittleThings ? '  YES ' : '  --- ';
  const p = g.inPhangs ? '  YES ' : '  --- ';
  console.log('  │ ' + g.name.padEnd(12) + ' │' + t + ' │' + lt + ' │' + p + ' │  ' + 
    (g.Q||'?') + '   │  ' + (g.inc||'?').toString().padStart(2) + '  │  ' + 
    (g.fD||'?') + '  │  ' + (g.logA0||0).toFixed(3).padStart(6) + '    │');
});
console.log('  └──────────────┴────────┴────────┴────────┴──────┴──────┴──────┴────────────┘');
console.log();

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  CRITICAL ASSESSMENT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log();
console.log('  DATA WE HAVE IN PROJECT:');
console.log('    - SPARC 1D rotation curves: 56 galaxies (full working sample)');
console.log('    - All baseline variables computed: logMHI, wig, env, Sig0, meanRun');
console.log('    - Per-galaxy a0 fits from Phase 11');
console.log('    - Frozen baselines A (4 vars) and B (5 vars)');
console.log();
console.log('  DATA WE DO NOT HAVE IN PROJECT:');
console.log('    - NO 2D velocity field maps (FITS cubes)');
console.log('    - NO moment-1 maps');
console.log('    - NO residual velocity maps');
console.log('    - NO harmonic decomposition outputs');
console.log('    - NO HI column density maps');
console.log('    - NO multi-band photometry for SPS M/L');
console.log('    - NO tilted-ring analysis outputs');
console.log();
console.log('  OVERLAP STATISTICS:');
console.log('    THINGS ∩ working-56: 7 galaxies (12.5%)');
console.log('    LITTLE THINGS ∩ working-56: 0 galaxies (0%)');
console.log('    PHANGS ∩ working-56: 0 galaxies (0%)');
console.log('    Total external overlap: 7 galaxies (12.5%)');
console.log();

console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log('  FEASIBILITY VERDICT');
console.log('━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━');
console.log();
console.log('  Phases 57-61 (2D kinematics): BLOCKED');
console.log('    Reason: No 2D velocity field data products in project.');
console.log('    Even with THINGS overlap, only 7/56 galaxies = 12.5%.');
console.log('    N=7 is statistically underpowered for regression analysis.');
console.log('    Minimum for any test: ~20 galaxies (df considerations).');
console.log();
console.log('  Phase 62-63 (environmental processing): BLOCKED');
console.log('    Reason: No HI deficiency measurements, no infall time proxies.');
console.log('    Would need external HI deficiency catalogs or group catalogs.');
console.log();
console.log('  Phase 64-65 (true stellar M/L): BLOCKED');
console.log('    Reason: No multi-band photometry or SPS M/L estimates.');
console.log('    PHANGS has 0 overlap with working sample.');
console.log('    Would need external SPS catalogs matched to SPARC galaxies.');
console.log();
console.log('  Phase 66 (death match): POSSIBLE');
console.log('    Only if any Phase 57-65 produces a surviving variable.');
console.log();

console.log('╔══════════════════════════════════════════════════════════════════════════════════╗');
console.log('║  BOTTOM LINE                                                                   ║');
console.log('╚══════════════════════════════════════════════════════════════════════════════════╝');
console.log();
console.log('  The external-data program (Phases 57-65) CANNOT proceed with');
console.log('  current project data. We have:');
console.log('    - 56 galaxies with 1D rotation curves (SPARC)');
console.log('    - 7 galaxies that also appear in THINGS');
console.log('    - 0 actual 2D data products downloaded/processed');
console.log();
console.log('  To execute this program, we would need to:');
console.log('    1. Download THINGS HI data cubes from NRAO archive (~50 GB)');
console.log('    2. Process velocity fields using tilted-ring fitting software');
console.log('    3. Run harmonic decomposition (BBarolo, TiRiFiC, or 3DBarolo)');
console.log('    4. Extract kinematic parameters (asymmetry, lopsidedness, etc.)');
console.log('    5. Match to SPARC rotation curve galaxies');
console.log('    This is a multi-week data engineering project requiring');
console.log('    specialized radio astronomy software (GIPSY, BBarolo).');
console.log();
console.log('  ALTERNATIVE: Use PUBLISHED derived parameters from THINGS papers:');
console.log('    - Trachternach+2008: kinematic lopsidedness for ~19 THINGS galaxies');
console.log('    - de Blok+2008: tilted-ring analysis parameters');
console.log('    - Oh+2008: mass model parameters');
console.log('    But only 7 overlap with our working-56 — too few for regression.');
console.log();

// Save master table
const output = {
  phase: '56.0-new',
  program: 'External Data Integration',
  date: new Date().toISOString().split('T')[0],
  nTotal: nTotal,
  overlapCounts: {
    THINGS: nThings,
    LITTLE_THINGS: nLT,
    PHANGS: nPhangs,
    anyExternal: nOverlap,
    SPARConly: nTotal - nOverlap
  },
  availabilityMatrix: {
    '2D_velocity_field': { nAvail: n2D, inProject: false, usable: false },
    'moment1_map': { nAvail: n2D, inProject: false, usable: false },
    'residual_map': { nAvail: nResid, inProject: false, usable: false },
    'HI_deficiency': { nAvail: nHIDef, inProject: false, usable: false },
    'SPS_ML': { nAvail: nML, inProject: false, usable: false },
    '1D_rotation_curve': { nAvail: 56, inProject: true, usable: true },
    'baseline_variables': { nAvail: 56, inProject: true, usable: true }
  },
  overlapGalaxies: overlapGals.map(g => ({
    name: g.name,
    surveys: g.source,
    Q: g.Q, inc: g.inc, fD: g.fD,
    logA0: g.logA0
  })),
  galaxies: masterTable,
  feasibility: {
    'phases_57_61_2D_kinematics': 'BLOCKED — no 2D data, only 7 overlap galaxies',
    'phases_62_63_environment': 'BLOCKED — no HI deficiency or infall data',
    'phases_64_65_stellar_ML': 'BLOCKED — no multi-band photometry, 0 PHANGS overlap',
    'phase_66_death_match': 'CONDITIONAL — needs surviving variables from 57-65'
  }
};

fs.writeFileSync('public/phase56-new-master-table.json', JSON.stringify(output, null, 2));
console.log('  Results saved to public/phase56-new-master-table.json');
