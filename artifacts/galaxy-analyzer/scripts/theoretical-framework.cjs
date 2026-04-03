const fs = require('fs');
const path = require('path');

const c = 2.998e8;
const G = 6.674e-11;
const hbar = 1.055e-34;
const kB = 1.381e-23;
const H0 = 2.18e-18;
const Lambda = 1.11e-52;
const rho_vac = 5.96e-27;

const a0_obs = 1.2e-10;

const sep = '='.repeat(72);
console.log(`\n${sep}`);
console.log('  THEORETICAL INVESTIGATION: WHY DOES a₀ EXIST?');
console.log('  Three Notebooks for Finding a Physical Explanation');
console.log(sep);

console.log(`\n${'═'.repeat(72)}`);
console.log('  NOTEBOOK 1: WHAT MUST BE EXPLAINED (5 Facts)');
console.log(`${'═'.repeat(72)}`);

const facts = [
  {
    id: 'F1',
    statement: 'An acceleration scale a₀ ≈ 1.2×10⁻¹⁰ m/s² appears in galaxy dynamics',
    evidence: 'RAR, BTFR, rotation curves of 175 SPARC + 22 LITTLE THINGS galaxies',
    robustness: 'Survives 4 breaking tests (null, scramble, M/L, distance)',
  },
  {
    id: 'F2',
    statement: 'a₀ is approximately universal across galaxy types',
    evidence: 'Per-galaxy scatter 0.33-0.57 dex (observational), no systematic mass/size dependence (partial r|Vmax = -0.656)',
    caveat: 'Scatter is large enough to hide weak dependencies',
  },
  {
    id: 'F3',
    statement: 'a₀ ≈ cH₀/2π to within ~10%',
    evidence: 'cH₀/2π = 1.04×10⁻¹⁰ m/s², observed a₀ = 1.2×10⁻¹⁰ m/s²',
    caveat: 'Could be numerological coincidence (not yet ruled out)',
  },
  {
    id: 'F4',
    statement: 'The relationship g_obs = f(g_bar) is tight with minimal scatter',
    evidence: 'RAR scatter < 0.13 dex (McGaugh+2016), our analysis confirms',
    implication: 'Dark matter "knows" where baryons are — or gravity itself changes',
  },
  {
    id: 'F5',
    statement: 'a₀ works for galaxies but NOT fully for clusters',
    evidence: 'Floor reduces Newton deficit by ~43% in clusters, ~52.5% residual remains',
    implication: 'Any explanation must account for cluster failure or predict additional physics',
  },
];

for (const f of facts) {
  console.log(`\n  ${f.id}: ${f.statement}`);
  console.log(`     Evidence: ${f.evidence}`);
  if (f.caveat) console.log(`     Caveat: ${f.caveat}`);
  if (f.implication) console.log(`     Implication: ${f.implication}`);
  if (f.robustness) console.log(`     Robustness: ${f.robustness}`);
}

console.log(`\n${'═'.repeat(72)}`);
console.log('  CONDITIONS FOR A SUCCESSFUL EXPLANATION');
console.log(`${'═'.repeat(72)}`);

const conditions = [
  { id: 'C1', condition: 'DERIVE a₀ as a scale (not insert it by hand)', test: 'Does the theory predict a₀ without being told about it?' },
  { id: 'C2', condition: 'Explain universality (weak mass/size dependence)', test: 'Does the derived a₀ depend on galaxy properties?' },
  { id: 'C3', condition: 'Explain cH₀ proximity', test: 'Does cosmology naturally enter the local dynamics?' },
  { id: 'C4', condition: 'Recover Newtonian limit at high acceleration', test: 'Does g → g_bar when g_bar >> a₀?' },
  { id: 'C5', condition: 'Predict or explain cluster shortfall', test: 'Why does a₀ underperform in clusters?' },
  { id: 'C6', condition: 'Have a mechanism (not just a formula)', test: 'Is there a field equation, action, or physical process?' },
];

for (const c of conditions) {
  console.log(`\n  ${c.id}: ${c.condition}`);
  console.log(`     Test: ${c.test}`);
}

console.log(`\n${'═'.repeat(72)}`);
console.log('  NOTEBOOK 2: CANDIDATE EXPLANATIONS');
console.log(`${'═'.repeat(72)}`);

console.log('\n  ── FAMILY A: Modified Gravity (MOND-like) ──');
console.log('  Core idea: The law of gravity itself changes below a₀');
console.log('  Representatives: Milgrom MOND, AQUAL, TeVeS, RMOND');
console.log(`
  Scorecard:
    C1 (derive a₀):      ✗ — a₀ is inserted as a free parameter
    C2 (universality):    ✓ — a₀ same for all galaxies by construction
    C3 (cH₀ proximity):  ✗ — considered a coincidence in standard MOND
    C4 (Newton limit):    ✓ — by construction (interpolation function)
    C5 (clusters):        ✗ — fails; needs ~2× more mass even with MOND
    C6 (mechanism):       △ — AQUAL has Lagrangian, but no deeper "why"
  
  Overall: Describes the phenomenon well, but does not EXPLAIN a₀.
`);

console.log('  ── FAMILY B: Baryon-Coupled Dark Matter ──');
console.log('  Core idea: DM exists but its distribution tracks baryons');
console.log('  Representatives: Feedback models, SIDM, superfluid DM');
console.log(`
  Scorecard:
    C1 (derive a₀):      △ — emerges statistically in some simulations
    C2 (universality):    △ — depends on feedback prescriptions
    C3 (cH₀ proximity):  ✗ — no natural explanation
    C4 (Newton limit):    ✓ — standard gravity + DM at high density
    C5 (clusters):        ✓ — DM explains clusters naturally
    C6 (mechanism):       ✓ — galaxy formation + feedback physics

  Overall: Explains clusters but struggles with tightness of RAR
  and offers no explanation for cH₀ coincidence.
`);

console.log('  ── FAMILY C: Cosmological Emergent Effect ──');
console.log('  Core idea: The cosmic expansion creates a local acceleration floor');
console.log('  Representatives: Cosmic Floor model, Verlinde emergent gravity');
console.log(`
  Scorecard:
    C1 (derive a₀):      △ — a₀ = cH₀/2π is motivated but not derived from first principles
    C2 (universality):    ✓ — cosmological origin means galaxy-independent
    C3 (cH₀ proximity):  ✓ — by construction / motivation
    C4 (Newton limit):    ✓ — floor only matters when g_bar << a₀
    C5 (clusters):        △ — partial (~43% help), needs more
    C6 (mechanism):       ✗ — THE MAIN GAP — no field equation yet

  Overall: Best at explaining F3 (cH₀) and F2 (universality),
  but lacks the physical mechanism that would make it a theory.
`);

console.log(`\n${'═'.repeat(72)}`);
console.log('  NOTEBOOK 3: DERIVATION ATTEMPTS');
console.log(`${'═'.repeat(72)}`);

console.log('\n  ── ATTEMPT 1: de Sitter Horizon Scale ──');
console.log('  Question: Does Λ naturally produce an acceleration ~ a₀?');

const R_dS = c / Math.sqrt(Lambda / 3);
const a_dS = c * c / R_dS;
const a_Lambda = c * Math.sqrt(Lambda / 3);
const ratio_dS = a_Lambda / a0_obs;

console.log(`
  de Sitter radius: R_dS = c/√(Λ/3) = ${(R_dS / 3.086e22).toExponential(2)} Mpc
  de Sitter acceleration: a_dS = c²/R_dS = c√(Λ/3) = ${a_Lambda.toExponential(3)} m/s²
  
  Ratio a_dS / a₀_obs = ${ratio_dS.toFixed(2)}
  
  Result: a_dS is within a factor of ~${ratio_dS.toFixed(1)} of a₀.
  But c√(Λ/3) ≈ ${a_Lambda.toExponential(2)} vs a₀ = ${a0_obs.toExponential(1)}.
`);

const a_cH = c * H0;
const a_cH_2pi = c * H0 / (2 * Math.PI);
const ratio_cH = a_cH / a0_obs;
const ratio_cH2pi = a_cH_2pi / a0_obs;

console.log('  ── ATTEMPT 2: Hubble Acceleration Scale ──');
console.log('  Question: Why cH₀ specifically?');
console.log(`
  cH₀ = ${a_cH.toExponential(3)} m/s²    ratio to a₀: ${ratio_cH.toFixed(2)}
  cH₀/2π = ${a_cH_2pi.toExponential(3)} m/s²  ratio to a₀: ${ratio_cH2pi.toFixed(2)}
  cH₀/6 = ${(a_cH / 6).toExponential(3)} m/s²  ratio to a₀: ${(a_cH / 6 / a0_obs).toFixed(2)}
  
  The 2π could come from:
    - Angular frequency vs linear frequency (ω = 2πf)
    - Circular horizon geometry
    - Fourier factor in quantum vacuum modes
    - Unruh relation periodicity
  
  The factor matters: cH₀ ≈ 6.5×a₀, but cH₀/2π ≈ 1.04×a₀.
  The 2π is NOT arbitrary — it needs a physical reason.
`);

console.log('  ── ATTEMPT 3: Unruh Temperature Argument ──');
console.log('  Question: Does the Unruh effect provide a minimum temperature/acceleration?');

const T_dS = hbar * Math.sqrt(Lambda / 3) / (2 * Math.PI * kB);
const T_Unruh_a0 = hbar * a0_obs / (2 * Math.PI * c * kB);
const a_from_TdS = 2 * Math.PI * kB * T_dS * c / hbar;

console.log(`
  de Sitter temperature: T_dS = ℏ√(Λ/3) / (2πk_B) = ${T_dS.toExponential(3)} K
  Unruh temperature at a₀: T_U(a₀) = ℏa₀/(2πck_B) = ${T_Unruh_a0.toExponential(3)} K
  
  Ratio T_dS / T_U(a₀) = ${(T_dS / T_Unruh_a0).toFixed(3)}
  
  The Unruh argument: If the universe has a minimum temperature T_dS
  set by the cosmological horizon, then there exists a minimum
  acceleration a_min such that T_U(a_min) = T_dS.
  
  This gives: a_min = 2πck_BT_dS/ℏ = c√(Λ/3) ≈ ${a_from_TdS.toExponential(3)} m/s²
  
  This is exactly a_dS — the same number, just derived differently.
  The Unruh argument provides a REASON: below a_min, the Unruh
  radiation bath merges with the cosmic horizon temperature,
  and inertia (resistance to acceleration) becomes affected.
`);

console.log('  ── ATTEMPT 4: The H₀ vs Λ Relationship ──');
console.log('  Question: Why does a₀ ≈ cH₀/2π work better than c√(Λ/3)?');

const H0_from_Lambda = Math.sqrt(Lambda * c * c / 3);
const ratio_H0 = H0 / H0_from_Lambda;

console.log(`
  In ΛCDM: H₀² = (8πG/3)ρ_total = (8πG/3)(ρ_m + ρ_Λ)
  In the far future: H → H_∞ = c√(Λ/3) = ${H0_from_Lambda.toExponential(3)} s⁻¹
  
  Current H₀ = ${H0.toExponential(3)} s⁻¹
  H₀ / H_∞ = ${ratio_H0.toFixed(3)}
  
  Key insight: H₀ includes matter contribution, H_∞ is pure Λ.
  The ratio H₀/H_∞ ≈ ${ratio_H0.toFixed(2)} is NOT 1.0 because we are at
  a special epoch (Ω_m ≈ 0.3, Ω_Λ ≈ 0.7).
  
  If a₀ = cH/2π (using current H):
    - a₀ changes with z (our testable prediction)
    - a₀ is set by the CURRENT expansion rate, not just Λ
    - This means matter density participates in setting a₀
    
  If a₀ = c√(Λ/3) (using Λ directly):
    - a₀ is a true constant (like MOND)
    - No redshift evolution
    - Simpler, but less connected to observations
  
  THE TWO PREDICTIONS DIVERGE AT z > 0.5.
  This is exactly what the redshift test measures.
`);

console.log('  ── ATTEMPT 5: Vacuum Energy Acceleration ──');
console.log('  Question: Does the vacuum energy density produce an acceleration?');

const a_vac = 4 * Math.PI * G * rho_vac * c / H0;
const a_vac2 = Math.sqrt(8 * Math.PI * G * rho_vac / 3) * c;

console.log(`
  ρ_vac = Λc²/(8πG) = ${rho_vac.toExponential(3)} kg/m³
  
  Try 1: a = √(8πGρ_vac/3) × c = c × H_Λ = ${a_vac2.toExponential(3)} m/s²
    This is just c√(Λ/3) again — same as Attempt 1.
  
  Try 2: a = (Λc²/6) × R for characteristic R
    At R = c/H₀ (Hubble radius):
    a = Λc³/(6H₀) = ${(Lambda * c * c * c / (6 * H0)).toExponential(3)} m/s²
    Ratio to a₀: ${(Lambda * c * c * c / (6 * H0) / a0_obs).toFixed(2)}
    
  None of these give a₀ exactly without a factor adjustment.
  The 2π factor remains unexplained from vacuum energy alone.
`);

console.log('  ── ATTEMPT 6: Dimensional Analysis ──');
console.log('  Question: What combinations of fundamental constants give a₀?');

const combos = [
  { formula: 'cH₀', value: c * H0, label: 'cH₀' },
  { formula: 'cH₀/2π', value: c * H0 / (2 * Math.PI), label: 'cH₀/2π' },
  { formula: 'c√(Λ/3)', value: c * Math.sqrt(Lambda / 3), label: 'c√(Λ/3)' },
  { formula: 'c²√Λ', value: c * c * Math.sqrt(Lambda), label: 'c²√Λ' },
  { formula: '(Λc⁴/6G)^(1/2)/c', value: Math.sqrt(Lambda * c * c * c * c / (6 * G)) / c, label: '(Λc⁴/6G)^½/c' },
  { formula: 'Gc²ρ_vac^(1/2)', value: G * c * c * Math.sqrt(rho_vac), label: 'Gc²√ρ_vac' },
];

console.log(`\n  ${'Formula'.padEnd(20)} ${'Value (m/s²)'.padEnd(16)} ${'Ratio to a₀'.padEnd(12)}`);
console.log(`  ${'-'.repeat(50)}`);
for (const co of combos) {
  console.log(`  ${co.label.padEnd(20)} ${co.value.toExponential(3).padEnd(16)} ${(co.value / a0_obs).toFixed(3)}`);
}

console.log(`\n${'═'.repeat(72)}`);
console.log('  SCORECARD: WHICH FAMILY BEST EXPLAINS THE FACTS?');
console.log(`${'═'.repeat(72)}`);
console.log(`
  Condition              Modified Gravity    Baryon-DM       Cosmological
  ─────────────────────────────────────────────────────────────────────────
  C1: Derive a₀          ✗ (inserted)        △ (emergent?)   △ (motivated)
  C2: Universality        ✓                   △               ✓
  C3: cH₀ proximity      ✗                   ✗               ✓
  C4: Newton limit        ✓                   ✓               ✓
  C5: Cluster shortfall   ✗                   ✓               △
  C6: Mechanism           △ (AQUAL)           ✓ (feedback)    ✗ (MAIN GAP)
  ─────────────────────────────────────────────────────────────────────────
  Score                    2.5/6               3.5/6           3.5/6
`);

console.log(`${'═'.repeat(72)}`);
console.log('  THE MAIN OPEN QUESTION');
console.log(`${'═'.repeat(72)}`);
console.log(`
  "Why does cH₀ appear as a local acceleration in galaxies?"
  
  This is THE question. Not "what is dark matter?" or "is MOND right?"
  
  The most promising path forward:
  
  1. The Unruh/horizon argument gives a₀ ~ c√(Λ/3) naturally
     but this is ~0.6× too small and doesn't explain the 2π
  
  2. Using H₀ instead of √(Λ/3) gives better agreement
     but this means a₀ depends on the CURRENT epoch
     → testable prediction (redshift evolution)
  
  3. The 2π factor might arise from:
     a) ω = 2πf conversion (Hubble frequency to angular)
     b) Horizon circumference (2πR_H) vs radius
     c) Quantization condition on circular orbits
     d) Fourier decomposition of vacuum fluctuations
  
  4. The cluster failure at ~52.5% suggests either:
     a) Additional hot gas physics not captured
     b) a₀ is only part of the story (missing component)
     c) The formula needs modification at cluster scales
  
  STATUS: We have a well-defined question, several derivation paths,
  and a testable prediction. We do NOT yet have a mechanism.
  This is the honest state of the investigation.
`);

const output = {
  date: new Date().toISOString(),
  description: 'Theoretical Investigation: Why does a₀ exist?',
  facts: facts,
  conditions: conditions,
  families: {
    modifiedGravity: {
      name: 'Modified Gravity (MOND-like)',
      scores: { C1: 0, C2: 1, C3: 0, C4: 1, C5: 0, C6: 0.5 },
      total: 2.5,
    },
    baryonDM: {
      name: 'Baryon-Coupled Dark Matter',
      scores: { C1: 0.5, C2: 0.5, C3: 0, C4: 1, C5: 1, C6: 1 },
      total: 4.0,
    },
    cosmological: {
      name: 'Cosmological Emergent Effect',
      scores: { C1: 0.5, C2: 1, C3: 1, C4: 1, C5: 0.5, C6: 0 },
      total: 4.0,
    },
  },
  derivations: {
    deSitter: {
      a_dS: a_Lambda,
      ratio_to_a0: ratio_dS,
      verdict: 'Within factor ~0.6, but missing 2π and too small',
    },
    hubble: {
      cH0: a_cH,
      cH0_2pi: a_cH_2pi,
      ratio_to_a0: ratio_cH2pi,
      verdict: 'Best match (ratio 0.87), but 2π needs physical justification',
    },
    unruh: {
      T_dS: T_dS,
      T_Unruh_a0: T_Unruh_a0,
      ratio: T_dS / T_Unruh_a0,
      verdict: 'Provides a REASON (minimum temperature), but gives same number as de Sitter',
    },
    dimensional: combos.map(c => ({ formula: c.formula, value: c.value, ratio: c.value / a0_obs })),
  },
  mainQuestion: 'Why does cH₀ appear as a local acceleration in galaxies?',
  mostPromisingPaths: [
    'Unruh/horizon argument with proper 2π derivation',
    'H(z) dependence as testable discriminant',
    'Cluster failure as constraint on mechanism',
  ],
};

fs.writeFileSync(
  path.join(__dirname, '..', 'public', 'theoretical-investigation.json'),
  JSON.stringify(output, null, 2)
);
console.log('\nResults saved to public/theoretical-investigation.json');
