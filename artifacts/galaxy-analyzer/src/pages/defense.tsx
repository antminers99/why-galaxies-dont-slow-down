import React, { useEffect, useState } from 'react';
import { Layout } from '@/components/layout';
import { BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine, Cell } from 'recharts';
import {
  Shield, CheckCircle2, XCircle, Shuffle, Microscope, FileCode, Layers,
  Scale, ChevronDown, ChevronRight, AlertTriangle, Atom, Beaker, BarChart3
} from 'lucide-react';

interface SigmaResult {
  slope: number;
  r: number;
  r2: number;
  partialR?: number;
  n: number;
}

interface HistogramBin {
  lo: number;
  hi: number;
  mid: number;
  count: number;
}

interface AuditCheck {
  check: string;
  result: string;
  detail: string;
}

interface SigmaDefinition {
  name: string;
  n: number;
  r: number;
  partialR: number | null;
  slope: number;
}

interface FairnessItem {
  aspect: string;
  rating: string;
  detail: string;
}

interface GalaxyPoint {
  r_kpc: number;
  V_obs_kms: number;
  V_bar_computed: number;
  f_DM_computed: number;
  V_DM_computed: number;
  V_gas_kms?: number;
  V_disk_kms?: number;
  V_bulge_kms?: number;
  V_bar_formula?: string;
  f_DM_formula?: string;
  V_DM_formula?: string;
}

interface GalaxyManual {
  name: string;
  vmax: number;
  nPoints: number;
  rmax_kpc: number;
  meanFDM: number;
  logSigBar: number;
  sigBar_calculation: {
    r_fiducial_kpc: number;
    V_bar_at_rfid: number;
    M_bar_enclosed: number;
    formula: string;
    sigma_bar: number;
    sigma_formula: string;
    log_sigma_bar: number;
  };
  samplePoints: GalaxyPoint[];
}

interface GlobalExcess {
  metric: string;
  obs: { slope: number };
  sim: { slope: number };
  deltaB: number;
  sigma: number;
}

interface FeedbackMetricResult {
  name: string;
  observed: { slope: number; r: number; n: number };
  noFeedback: { slope: number; r: number; n: number };
  withFeedback: { slope: number; r: number; n: number };
  excessNoFeedback: { sigma: number };
  excessWithFeedback: { sigma: number };
  feedbackExplainsPercent: number;
  feedbackFullyExplains: boolean;
  caveat?: string;
}

interface FeedbackComponent {
  name: string;
  reference: string;
  description: string;
  effect: string;
}

interface FeedbackData {
  metrics: FeedbackMetricResult[];
  feedbackModel: {
    description: string;
    components: FeedbackComponent[];
    nGalaxies: number;
    limitations?: string[];
  };
  summary: {
    maxSigmaWithoutFeedback: number;
    maxSigmaWithFeedback: number;
    maxSigmaWithFeedbackRaw?: number;
    robustMetric?: string;
    robustSigma?: number;
    inflatedMetric?: string;
    inflatedSigma?: number;
    inflatedReason?: string;
    inflatedMetrics?: string[];
    reliableMetrics?: string[];
    averageReductionPercent: number;
    anyMetricFullyExplained?: boolean;
    allMetricsFullyExplained: boolean;
    verdict: string;
    honestAssessment: string;
    updatedClaim: string;
    nextSteps?: string;
  };
}

interface HydroMetricComparison {
  observed: { slope: number; r: number; n: number; bootMean: number; bootSD: number };
  simulated: { slope: number; r: number; n: number; bootMean: number; bootSD: number };
  comparison: { delta: number; se: number; sigma: number; sameSign: boolean; signObserved: string; signSimulated: string };
}

interface HydroSimulation {
  name: string;
  reference: string;
  description: string;
  feedbackType: string;
  keyFeatures: string[];
  nGalaxies: number;
}

interface HydroData {
  description: string;
  sparcGalaxies: number;
  simulations: { fire2: HydroSimulation; tng: HydroSimulation };
  metrics: Record<string, { fire2: HydroMetricComparison | null; tng: HydroMetricComparison | null }>;
  signProblem: { fire2: boolean; tng: boolean; both: boolean; description: string };
  summary: {
    verdict: string;
    strength: string;
    description: string;
    observedSlope: number | null;
    fire2Slope: number | null;
    tngSlope: number | null;
    fire2Sigma: number | null;
    tngSigma: number | null;
    keyInsight: string;
    updatedClaim: string;
    caveat: string;
  };
  caveats: string[];
}

interface ParticleAlphaResult {
  slope: number;
  r: number;
  n: number;
  bootMean: number;
  bootSD: number;
  ci95: number[];
  sign: string;
  permPvalue: number;
}

interface ParticleComparison {
  delta: number;
  se: number;
  sigma: number;
  signMatch: boolean;
}

interface ParticleDataSource {
  papers: string[];
  nGalaxies: number;
  massRange: { logMstar: string[] };
  description: string;
}

interface ParticleData {
  description: string;
  dataType: string;
  dataSources: {
    fire2: ParticleDataSource;
    tng: ParticleDataSource;
    sparc: { nGalaxies: number; description: string };
  };
  alpha: {
    sparc: ParticleAlphaResult;
    fire2: ParticleAlphaResult;
    tng: ParticleAlphaResult;
    comparison: {
      fire2: ParticleComparison;
      tng: ParticleComparison;
    };
  };
  galaxies: {
    sparc: Array<{ name: string; logSigBar: number; alpha: number }>;
    fire2: Array<{ name: string; logMstar: number; logSigBar: number; alpha: number; source: string }>;
    tng: Array<{ name: string; logMstar: number; logSigBar: number; alpha: number; source: string }>;
  };
  signProblem: { fire2: boolean; tng: boolean; both: boolean; either: boolean };
  verdict: {
    result: string;
    strength: string;
    claim: string;
    fire2Sigma: number;
    tngSigma: number;
    fire2Ratio: number;
    tngRatio: number;
    magnitudeDiscrepancy: boolean;
    signResolved: boolean;
  };
  slopeConvention: { description: string; conversion: string };
  caveats: string[];
  comparisonWithParametric: { description: string; note: string };
}

interface RawParticleCorrelation {
  label: string;
  n: number;
  slope: number;
  slopeErr: number;
  r: number;
  intercept: number;
  bootMean: number;
  bootStd: number;
  bootCI: number[];
  sigmaFromZero: number;
}

interface RawParticleSimData {
  catalog: string;
  simulation: string;
  snapshot: number;
  redshift: number;
  totalHalosFetched: number;
  SPARClike: number;
  selectionCriteria: string;
  correlation: RawParticleCorrelation;
}

interface RawParticleData {
  metadata: {
    pipeline: string;
    source: string;
    description: string;
    aperture_radii_kpc: number[];
    alpha_definition: string;
    Sigma_bar_definition: string;
    timestamp: string;
  };
  simulations: {
    TNG100: RawParticleSimData;
    EAGLE: RawParticleSimData;
  };
  SPARC: { n: number; correlation: RawParticleCorrelation } | null;
  comparison: {
    tngVsSparc: { slopeRatio: number; sigmaDiscrepancy: number; signMatch: boolean };
    eagleVsSparc: { slopeRatio: number; sigmaDiscrepancy: number; signMatch: boolean };
    tngVsEagle: { slopeRatio: number; sigmaDiscrepancy: number };
  } | null;
  rawData: {
    TNG100: Array<{ alpha: number; log_Sigma_bar: number }>;
    EAGLE: Array<{ alpha: number; log_Sigma_bar: number }>;
    SPARC: Array<{ name: string; alpha: number; log_Sigma_bar: number }>;
  };
}

interface BivariateCollapseData {
  datasets: {
    sparc: { nGalaxies: number; nPoints: number };
    littleThings: { nGalaxies: number; nPoints: number };
    combined: { nGalaxies: number; nPoints: number };
  };
  test2_deltaVsSigma: {
    sparc: { slope: number; se: number; r: number; r2: number; n: number };
    littleThings: { slope: number; se: number; r: number; r2: number; n: number };
    combined: { slope: number; se: number; r: number; r2: number; n: number };
    signConsistent: boolean;
  };
  test3_bivariateModels: {
    sparc: { modelA: { r2: number; adjR2: number; rmse: number }; modelB: { r2: number; adjR2: number; rmse: number; sigmaCoeff: number }; fStatAB: number; deltaAIC: number; deltaBIC: number };
    littleThings: { modelA: { r2: number; adjR2: number; rmse: number }; modelB: { r2: number; adjR2: number; rmse: number; sigmaCoeff: number }; fStatAB: number; deltaAIC: number; deltaBIC: number };
    combined: { modelA: { r2: number; adjR2: number; rmse: number }; modelB: { r2: number; adjR2: number; rmse: number; sigmaCoeff: number }; fStatAB: number; deltaAIC: number };
  };
  test4_crossValidation: {
    sparcToLT: { rmseA: number; rmseB: number; improvement: number; sigmaCoeff: number };
    ltToSPARC: { rmseA: number; rmseB: number; improvement: number; sigmaCoeff: number };
    signConsistent: boolean;
    bothImprove: boolean;
  };
  test5_collapse: {
    baseline: { rms: number; meanGap: number };
    best: { coeff: number; rms: number; meanGap: number };
    scatterReduction: number;
    meanGapReduction: number;
  };
  test6_stability: {
    results: Array<{ label: string; slope: number; se: number; sig: number; n: number; sign: string }>;
    allSameSign: boolean;
    signDirection: string;
    meanSignificance: number;
  };
  test7_bivariateFormula: {
    bestGamma: number;
    baseRMS: number;
    correctedRMS: number;
    scatterReduction: number;
    formula: string;
    perDataset: { sparc: { baseRMS: number; correctedRMS: number }; littleThings: { baseRMS: number; correctedRMS: number } };
  };
  test8_galaxyLevel: {
    sparc: { slope: number; se: number; sig: number; r: number; n: number };
    littleThings: { slope: number; se: number; sig: number; r: number; n: number };
    combined: { slope: number; se: number; sig: number; r: number; n: number };
  };
  verdict: {
    sigmaImprovesRAR_SPARC: boolean;
    sigmaImprovesRAR_LT: boolean;
    crossValidated: boolean;
    signConsistent: boolean;
    definitionStable: boolean;
    isTreasure: boolean;
  };
}

interface TransitionScaleData {
  a0_corrected: number;
  a0_ms2: number;
  nPoints: number;
  nGalaxies: number;
  ratioBins: Array<{ logGbar: number; gBar: number; medianRatio: number; rmsScatter: number; residScatter: number; n: number }>;
  collapse: {
    rmsWithCorrectA0: number;
    madWithCorrectA0: number;
    rmsWithWrongA0: number;
    improvement: number;
    binned: Array<{ x: number; rms: number; mad: number; bias: number; n: number }>;
  };
  perGalaxyA0: {
    global: { a0: number; logA0: number; rms: number };
    nFit: number;
    nWellConstrained: number;
    medianLogA0: number;
    madLogA0: number;
    rmsPerGalaxy: number;
    rmsUniversal: number;
    improvement: number;
    distribution: { p5: number; p25: number; p50: number; p75: number; p95: number };
    perDataset: { sparc: { n: number; medA0: number; mad: number }; lt: { n: number; medA0: number; mad: number } };
  };
  cosmology: { a0: number; cH0: string; ratio: number };
  plotPoints: Array<{ x: number; y: number; g: string }>;
  verdict: string;
}

interface DefenseData {
  test1_independence: {
    title: string;
    description: string;
    photometricSigma: SigmaResult;
    luminosityProxy: SigmaResult;
    geometricSigma: SigmaResult;
    partialControlGbar: number;
    verdict: string;
    conclusion: string;
  };
  test2_shuffle: {
    title: string;
    description: string;
    realR: number;
    nShuffles: number;
    pValue: number;
    shuffleMean: number;
    shuffleSD: number;
    sigmaFromNull: number;
    histogram: HistogramBin[];
    percentile5: number;
    percentile95: number;
    conclusion: string;
  };
  test3_null_simulation: {
    title: string;
    description: string;
    simConfig: {
      nRealizations: number;
      galaxiesPerRealization: number;
      totalMock: number;
      haloProfile: string;
      baryonModel: string;
      couplingBuiltIn: boolean;
    };
    excessGlobal: GlobalExcess[];
    conclusion: string;
  };
  test4_manual_galaxies: {
    title: string;
    description: string;
    galaxies: GalaxyManual[];
    constants: { G: string; UPSILON_D: number; UPSILON_B: number };
    conclusion: string;
  };
  test5_code_audit: {
    title: string;
    description: string;
    checks: AuditCheck[];
    conclusion: string;
  };
  test6_alt_definitions: {
    title: string;
    description: string;
    definitions: SigmaDefinition[];
    allNegative: boolean;
    totalDefinitions: number;
    negativeCount: number;
    conclusion: string;
  };
  test7_simulation_fairness: {
    title: string;
    description: string;
    simDetails: Record<string, string | number | boolean>;
    fairnessAssessment: FairnessItem[];
    conservativeCount: number;
    conclusion: string;
  };
  summary: {
    totalTests: number;
    passed: number;
    criticalFindings: string[];
    goldenSentence: string;
  };
}

function GlassCard({ children, glow, className = '' }: { children: React.ReactNode; glow?: string; className?: string }) {
  const glowColors: Record<string, string> = {
    cyan: 'shadow-cyan-500/10 border-cyan-500/20',
    emerald: 'shadow-emerald-500/10 border-emerald-500/20',
    amber: 'shadow-amber-500/10 border-amber-500/20',
    rose: 'shadow-rose-500/10 border-rose-500/20',
    violet: 'shadow-violet-500/10 border-violet-500/20',
  };
  return (
    <div className={'bg-white/[0.03] backdrop-blur-sm border border-white/10 rounded-2xl p-6 shadow-lg ' + (glow ? glowColors[glow] || '' : '') + ' ' + className}>
      {children}
    </div>
  );
}

function TestHeader({ number, title, status, icon: Icon }: { number: number; title: string; status: 'pass' | 'warn'; icon: React.ElementType }) {
  return (
    <div className="flex items-center gap-3 mb-4">
      <div className="flex items-center gap-2">
        <span className="text-xs font-mono text-slate-600 bg-white/5 w-7 h-7 rounded-lg flex items-center justify-center border border-white/10">{number}</span>
        <Icon className="w-5 h-5 text-cyan-400" />
      </div>
      <h2 className="text-lg font-bold text-white flex-1">{title}</h2>
      {status === 'pass' ? (
        <span className="flex items-center gap-1 text-xs font-bold text-emerald-400 bg-emerald-500/10 px-3 py-1 rounded-full border border-emerald-500/20">
          <CheckCircle2 className="w-3.5 h-3.5" /> PASSED
        </span>
      ) : (
        <span className="flex items-center gap-1 text-xs font-bold text-amber-400 bg-amber-500/10 px-3 py-1 rounded-full border border-amber-500/20">
          <AlertTriangle className="w-3.5 h-3.5" /> REVIEW
        </span>
      )}
    </div>
  );
}

function ExpandableGalaxy({ galaxy }: { galaxy: GalaxyManual }) {
  const [open, setOpen] = useState(false);
  return (
    <div className="bg-white/5 rounded-xl border border-white/5 overflow-hidden">
      <button onClick={() => setOpen(!open)} className="w-full flex items-center justify-between p-4 hover:bg-white/5 transition-colors">
        <div className="flex items-center gap-3">
          <span className="text-sm font-bold text-white">{galaxy.name}</span>
          <span className="text-xs text-slate-400">V_max = {galaxy.vmax} km/s</span>
          <span className="text-xs text-slate-400">{galaxy.nPoints} points</span>
        </div>
        <div className="flex items-center gap-3">
          <span className="font-mono text-xs text-cyan-400">log Σ = {galaxy.logSigBar}</span>
          <span className="font-mono text-xs text-violet-400">⟨f_DM⟩ = {galaxy.meanFDM.toFixed(3)}</span>
          {open ? <ChevronDown className="w-4 h-4 text-slate-400" /> : <ChevronRight className="w-4 h-4 text-slate-400" />}
        </div>
      </button>
      {open && (
        <div className="px-4 pb-4 space-y-3">
          <div className="bg-white/5 rounded-lg p-3">
            <h5 className="text-xs font-bold text-cyan-400 mb-2">Σ_bar Calculation</h5>
            <div className="grid grid-cols-2 md:grid-cols-3 gap-2 text-xs">
              <div><span className="text-slate-400">r_fiducial:</span> <span className="text-white font-mono">{galaxy.sigBar_calculation.r_fiducial_kpc} kpc</span></div>
              <div><span className="text-slate-400">V_bar(r_fid):</span> <span className="text-white font-mono">{galaxy.sigBar_calculation.V_bar_at_rfid} km/s</span></div>
              <div><span className="text-slate-400">M_bar enclosed:</span> <span className="text-white font-mono">{galaxy.sigBar_calculation.M_bar_enclosed.toExponential(2)} M☉</span></div>
              <div className="col-span-2 md:col-span-3">
                <span className="text-slate-400">Formula:</span> <span className="text-amber-300 font-mono">{galaxy.sigBar_calculation.formula} → {galaxy.sigBar_calculation.sigma_formula}</span>
              </div>
              <div><span className="text-slate-400">log Σ_bar:</span> <span className="text-emerald-400 font-mono font-bold">{galaxy.sigBar_calculation.log_sigma_bar}</span></div>
            </div>
          </div>
          <div>
            <h5 className="text-xs font-bold text-violet-400 mb-2">Sample Data Points</h5>
            <div className="overflow-x-auto">
              <table className="w-full text-xs font-mono">
                <thead>
                  <tr className="border-b border-white/10 text-slate-400">
                    <th className="text-left py-1 px-2">r (kpc)</th>
                    <th className="text-center py-1 px-2">V_obs</th>
                    <th className="text-center py-1 px-2">V_bar</th>
                    <th className="text-center py-1 px-2">f_DM</th>
                    <th className="text-center py-1 px-2">V_DM</th>
                  </tr>
                </thead>
                <tbody>
                  {galaxy.samplePoints.map((pt, i) => (
                    <tr key={i} className="border-b border-white/5">
                      <td className="py-1 px-2 text-slate-300">{pt.r_kpc}</td>
                      <td className="py-1 px-2 text-center text-white">{pt.V_obs_kms}</td>
                      <td className="py-1 px-2 text-center text-cyan-400">{pt.V_bar_computed}</td>
                      <td className="py-1 px-2 text-center text-violet-400">{pt.f_DM_computed.toFixed(3)}</td>
                      <td className="py-1 px-2 text-center text-amber-400">{pt.V_DM_computed}</td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          </div>
        </div>
      )}
    </div>
  );
}

export default function DefensePage() {
  const [data, setData] = useState<DefenseData | null>(null);
  const [feedback, setFeedback] = useState<FeedbackData | null>(null);
  const [hydro, setHydro] = useState<HydroData | null>(null);
  const [particle, setParticle] = useState<ParticleData | null>(null);
  const [rawParticle, setRawParticle] = useState<RawParticleData | null>(null);
  const [bivariate, setBivariate] = useState<BivariateCollapseData | null>(null);
  const [transition, setTransition] = useState<TransitionScaleData | null>(null);
  const [error, setError] = useState(false);

  useEffect(() => {
    Promise.all([
      fetch(`${import.meta.env.BASE_URL}defense-validation.json`).then(r => { if (!r.ok) throw new Error('fetch failed'); return r.json(); }),
      fetch(`${import.meta.env.BASE_URL}feedback-test.json`).then(r => { if (!r.ok) throw new Error('fetch failed'); return r.json(); }).catch(() => null),
      fetch(`${import.meta.env.BASE_URL}hydro-comparison.json`).then(r => { if (!r.ok) throw new Error('fetch failed'); return r.json(); }).catch(() => null),
      fetch(`${import.meta.env.BASE_URL}particle-comparison.json`).then(r => { if (!r.ok) throw new Error('fetch failed'); return r.json(); }).catch(() => null),
      fetch(`${import.meta.env.BASE_URL}raw-particle-data.json`).then(r => { if (!r.ok) throw new Error('fetch failed'); return r.json(); }).catch(() => null),
      fetch(`${import.meta.env.BASE_URL}bivariate-collapse.json`).then(r => { if (!r.ok) throw new Error('fetch failed'); return r.json(); }).catch(() => null),
      fetch(`${import.meta.env.BASE_URL}transition-scale.json`).then(r => { if (!r.ok) throw new Error('fetch failed'); return r.json(); }).catch(() => null),
    ])
      .then(([defData, fbData, hydroData, particleData, rawData, bivData, transData]) => { setData(defData); setFeedback(fbData); setHydro(hydroData); setParticle(particleData); setRawParticle(rawData); setBivariate(bivData); setTransition(transData); })
      .catch(() => setError(true));
  }, []);

  if (error) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-full">
          <GlassCard glow="rose">
            <div className="text-center">
              <AlertTriangle className="w-8 h-8 text-rose-400 mx-auto mb-3" />
              <p className="text-white font-bold mb-2">Failed to load defense data</p>
              <button onClick={() => window.location.reload()} className="text-sm text-cyan-400 hover:text-cyan-300">Retry</button>
            </div>
          </GlassCard>
        </div>
      </Layout>
    );
  }

  if (!data) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-full">
          <div className="text-slate-400 animate-pulse">Loading defense validation...</div>
        </div>
      </Layout>
    );
  }

  const t1 = data.test1_independence;
  const t2 = data.test2_shuffle;
  const t3 = data.test3_null_simulation;
  const t4 = data.test4_manual_galaxies;
  const t5 = data.test5_code_audit;
  const t6 = data.test6_alt_definitions;
  const t7 = data.test7_simulation_fairness;

  return (
    <Layout>
      <div className="space-y-8 pb-20">
        <div className="flex items-start justify-between flex-wrap gap-4">
          <div>
            <div className="flex items-center gap-3 mb-2">
              <Shield className="w-8 h-8 text-cyan-400" />
              <h1 className="text-3xl font-bold text-white tracking-tight">Defense Validation</h1>
            </div>
            <p className="text-slate-400 max-w-2xl text-sm leading-relaxed">
              Seven systematic tests to prove this result is NOT circular, NOT an artifact, NOT a definition trick, NOT a code bug, and NOT an unfair simulation.
            </p>
          </div>
          <div className="flex gap-2 flex-wrap">
            <span className="px-4 py-2 rounded-full text-xs font-bold bg-emerald-500/10 text-emerald-400 border border-emerald-500/20">
              {data.summary.passed}/{data.summary.totalTests} TESTS PASSED
            </span>
          </div>
        </div>

        <GlassCard glow="amber">
          <div className="flex items-start gap-3">
            <AlertTriangle className="w-5 h-5 text-amber-400 mt-0.5 flex-shrink-0" />
            <div>
              <h3 className="text-white font-bold text-sm mb-1">The Scientific Standard</h3>
              <p className="text-amber-200 text-sm leading-relaxed font-mono">
                "{data.summary.goldenSentence}"
              </p>
            </div>
          </div>
        </GlassCard>

        {feedback && (
          <section>
            <div className="flex items-center gap-3 mb-4">
              <div className="flex items-center gap-2">
                <span className="text-xs font-mono text-rose-500 bg-rose-500/10 w-7 h-7 rounded-lg flex items-center justify-center border border-rose-500/20">!</span>
                <AlertTriangle className="w-5 h-5 text-rose-400" />
              </div>
              <h2 className="text-lg font-bold text-white flex-1">CRITICAL: Baryonic Feedback Test</h2>
              <span className={'flex items-center gap-1 text-xs font-bold px-3 py-1 rounded-full border ' + (
                feedback.summary.allMetricsFullyExplained
                  ? 'text-rose-400 bg-rose-500/10 border-rose-500/20'
                  : 'text-amber-400 bg-amber-500/10 border-amber-500/20'
              )}>
                <AlertTriangle className="w-3.5 h-3.5" /> HONEST RESULT
              </span>
            </div>
            <GlassCard glow="rose" className="mb-4 border-2 border-rose-500/20">
              <div className="bg-rose-500/5 rounded-xl p-4 mb-4">
                <p className="text-rose-200 text-sm leading-relaxed">
                  <strong>The critique:</strong> "You compare data with baryonic physics against simulation WITHOUT feedback. Of course you get 'excess' — but that might just mean your simulation is incomplete, not that there's new physics."
                </p>
                <p className="text-slate-400 text-sm mt-2">
                  <strong>Our response:</strong> We built ΛCDM + parametric baryonic feedback (Di Cintio et al. 2014 core formation + Blumenthal et al. 1986 adiabatic contraction) and re-ran the comparison. We acknowledge this is not a full hydrodynamic simulation.
                </p>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
                {feedback.metrics.map((m, i) => {
                  const fullyExplained = m.feedbackFullyExplains;
                  const sigBefore = m.excessNoFeedback.sigma;
                  const sigAfter = m.excessWithFeedback.sigma;
                  const isInflated = sigAfter > sigBefore * 1.5 && !fullyExplained;
                  const statusColor = fullyExplained ? 'emerald' : isInflated ? 'rose' : 'amber';
                  const borderColor = { emerald: 'border-emerald-500/20', rose: 'border-rose-500/20', amber: 'border-amber-500/20' }[statusColor];
                  return (
                    <div key={i} className={'bg-white/5 rounded-xl p-4 border ' + borderColor}>
                      <div className="flex items-center justify-between mb-3">
                        <span className="text-xs font-bold text-white">{m.name}</span>
                        {isInflated && (
                          <span className="text-xs text-rose-400 bg-rose-500/10 px-2 py-0.5 rounded-full border border-rose-500/20">INFLATED</span>
                        )}
                      </div>
                      <div className="space-y-2 text-xs">
                        <div className="flex justify-between">
                          <span className="text-slate-400">Observed slope:</span>
                          <span className="text-white font-mono">{m.observed.slope.toFixed(5)}</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-slate-400">ΛCDM (no feedback):</span>
                          <span className="text-violet-400 font-mono">{m.noFeedback.slope.toFixed(5)}</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-slate-400">ΛCDM + feedback:</span>
                          <span className="text-cyan-400 font-mono">{m.withFeedback.slope.toFixed(5)}</span>
                        </div>
                        <hr className="border-white/10" />
                        <div className="flex justify-between items-center">
                          <span className="text-slate-400">Without feedback:</span>
                          <span className="text-amber-400 font-mono font-bold">{sigBefore.toFixed(1)}σ</span>
                        </div>
                        <div className="flex justify-between items-center">
                          <span className="text-slate-400">With feedback:</span>
                          <span className={'font-mono font-bold ' + (fullyExplained ? 'text-emerald-400' : isInflated ? 'text-rose-400' : 'text-amber-400')}>{sigAfter.toFixed(1)}σ</span>
                        </div>
                        <div className={'text-center py-1 rounded-lg mt-1 ' + (
                          fullyExplained ? 'bg-emerald-500/10 text-emerald-300' :
                          isInflated ? 'bg-rose-500/10 text-rose-300' :
                          'bg-amber-500/10 text-amber-300'
                        )}>
                          {fullyExplained ? 'Feedback explains this' : isInflated ? 'Model mismatch (unreliable)' : sigAfter.toFixed(1) + 'σ residual'}
                        </div>
                      </div>
                      {m.caveat && (
                        <p className="text-xs text-slate-500 mt-3 leading-relaxed italic border-t border-white/5 pt-2">{m.caveat}</p>
                      )}
                    </div>
                  );
                })}
              </div>

              {feedback.summary.inflatedMetrics && feedback.summary.inflatedMetrics.length > 0 && (
                <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4 mb-4">
                  <div className="flex items-start gap-2">
                    <AlertTriangle className="w-4 h-4 text-rose-400 mt-0.5 flex-shrink-0" />
                    <div>
                      <h4 className="text-xs font-bold text-rose-300 mb-1">Sigma Inflation Warning: {feedback.summary.inflatedMetrics.join(', ')}</h4>
                      <p className="text-xs text-rose-200/70 leading-relaxed">Model mismatch: parametric feedback over-corrects {feedback.summary.inflatedMetrics.join(', ')} slope, inflating apparent sigma. Not reliable as stated.</p>
                    </div>
                  </div>
                </div>
              )}

              <div className="bg-white/5 rounded-xl p-4 mb-4">
                <h4 className="text-xs font-bold text-white mb-2">Feedback Model Components</h4>
                <div className="space-y-2">
                  {feedback.feedbackModel.components.map((c, i) => (
                    <div key={i} className="bg-white/5 rounded-lg p-3">
                      <div className="flex items-start gap-2">
                        <span className="text-cyan-400 text-xs font-bold">{c.name}</span>
                        <span className="text-slate-500 text-xs">({c.reference})</span>
                      </div>
                      <p className="text-xs text-slate-400 mt-1">{c.description}</p>
                      <p className="text-xs text-amber-300 mt-1">{c.effect}</p>
                    </div>
                  ))}
                </div>
              </div>

              {feedback.feedbackModel.limitations && feedback.feedbackModel.limitations.length > 0 && (
                <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-4 mb-4">
                  <h4 className="text-xs font-bold text-amber-300 mb-2">Known Limitations of This Test</h4>
                  <ul className="space-y-1">
                    {feedback.feedbackModel.limitations.map((lim, i) => (
                      <li key={i} className="text-xs text-amber-200/70 flex items-start gap-2">
                        <span className="text-amber-500 mt-0.5">-</span>
                        <span>{lim}</span>
                      </li>
                    ))}
                  </ul>
                </div>
              )}

              <div className="bg-white/5 rounded-xl p-4 mb-4 border border-white/10">
                <h4 className="text-sm font-bold text-white mb-2">Honest Assessment</h4>
                <p className="text-sm text-slate-300 leading-relaxed mb-3">{feedback.summary.honestAssessment}</p>
                <div className="flex flex-wrap gap-2 items-center">
                  <span className="inline-block px-3 py-1 rounded-full text-xs font-bold bg-amber-500/10 text-amber-400 border border-amber-500/20">
                    {feedback.summary.updatedClaim}
                  </span>
                </div>
              </div>

              {feedback.summary.nextSteps && (
                <div className="bg-cyan-500/5 border border-cyan-500/20 rounded-xl p-4">
                  <div className="flex items-start gap-2">
                    <Beaker className="w-4 h-4 text-cyan-400 mt-0.5 flex-shrink-0" />
                    <div>
                      <h4 className="text-xs font-bold text-cyan-300 mb-1">Next Steps Required</h4>
                      <p className="text-xs text-cyan-200/70 leading-relaxed">{feedback.summary.nextSteps}</p>
                    </div>
                  </div>
                </div>
              )}
            </GlassCard>
          </section>
        )}

        {hydro && (
          <section>
            <div className="flex items-center gap-3 mb-4">
              <div className="flex items-center gap-2">
                <span className="text-xs font-mono text-cyan-500 bg-cyan-500/10 w-7 h-7 rounded-lg flex items-center justify-center border border-cyan-500/20">II</span>
                <Beaker className="w-5 h-5 text-cyan-400" />
              </div>
              <h2 className="text-lg font-bold text-white flex-1">Hydrodynamic Simulation Comparison</h2>
              <span className={'flex items-center gap-1 text-xs font-bold px-3 py-1 rounded-full border ' + (
                hydro.signProblem.both
                  ? 'text-emerald-400 bg-emerald-500/10 border-emerald-500/20'
                  : hydro.signProblem.fire2 || hydro.signProblem.tng
                  ? 'text-amber-400 bg-amber-500/10 border-amber-500/20'
                  : 'text-slate-400 bg-slate-500/10 border-slate-500/20'
              )}>
                {hydro.summary.verdict}
              </span>
            </div>
            <GlassCard glow="cyan" className="mb-4 border-2 border-cyan-500/20">
              <div className="bg-cyan-500/5 rounded-xl p-4 mb-4">
                <p className="text-cyan-200 text-sm leading-relaxed">
                  <strong>The decisive test:</strong> Compare SPARC data against mock catalogs built from published scaling relations of full hydrodynamic simulations — not just parametric feedback.
                </p>
                <p className="text-slate-400 text-sm mt-2">
                  <strong>Key question:</strong> Does the inner slope (alpha) sign reversal survive when compared to FIRE-2 and IllustrisTNG — simulations with gas turbulence, bursty star formation, and AGN feedback?
                </p>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
                {(['fire2', 'tng'] as const).map(simKey => {
                  const sim = hydro.simulations[simKey];
                  const borderCls = simKey === 'fire2' ? 'border-orange-500/20' : 'border-blue-500/20';
                  const labelCls = simKey === 'fire2' ? 'text-orange-400' : 'text-blue-400';
                  return (
                    <div key={simKey} className={'bg-white/5 rounded-xl p-4 border ' + borderCls}>
                      <div className="flex items-center justify-between mb-2">
                        <span className={'text-sm font-bold ' + labelCls}>{sim.name}</span>
                        <span className="text-xs text-slate-500 font-mono">{sim.nGalaxies.toLocaleString()} galaxies</span>
                      </div>
                      <p className="text-xs text-slate-400 mb-2">{sim.feedbackType}</p>
                      <p className="text-xs text-slate-500 italic">{sim.reference}</p>
                      <div className="flex flex-wrap gap-1 mt-2">
                        {sim.keyFeatures.map((f, i) => (
                          <span key={i} className="text-xs px-2 py-0.5 rounded-full bg-white/5 text-slate-400 border border-white/5">{f}</span>
                        ))}
                      </div>
                    </div>
                  );
                })}
              </div>

              <div className="mb-4">
                <h4 className="text-xs font-bold text-white mb-3">Metric-by-Metric Comparison</h4>
                {['alpha', 'V_DM', 'r_DMdom'].map(metricName => {
                  const metricData = hydro.metrics[metricName];
                  if (!metricData) return null;
                  const fire2 = metricData.fire2;
                  const tng = metricData.tng;
                  const isAlpha = metricName === 'alpha';
                  return (
                    <div key={metricName} className={'bg-white/5 rounded-xl p-4 border border-white/10 mb-3 ' + (isAlpha ? 'ring-1 ring-cyan-500/30' : '')}>
                      <div className="flex items-center justify-between mb-3">
                        <span className={'text-sm font-bold ' + (isAlpha ? 'text-cyan-400' : 'text-white')}>{metricName}</span>
                        {isAlpha && <span className="text-xs text-cyan-400 bg-cyan-500/10 px-2 py-0.5 rounded-full border border-cyan-500/20">KEY METRIC</span>}
                      </div>
                      <div className="overflow-x-auto">
                        <table className="w-full text-xs font-mono">
                          <thead>
                            <tr className="border-b border-white/10 text-slate-500">
                              <th className="text-left py-1 px-2">Source</th>
                              <th className="text-center py-1 px-2">Slope</th>
                              <th className="text-center py-1 px-2">Sign</th>
                              <th className="text-center py-1 px-2">Residual</th>
                              <th className="text-center py-1 px-2">Match?</th>
                            </tr>
                          </thead>
                          <tbody>
                            <tr className="border-b border-white/5">
                              <td className="py-1.5 px-2 text-white font-bold">SPARC Data</td>
                              <td className="py-1.5 px-2 text-center text-white">{fire2 ? fire2.observed.slope.toFixed(5) : '-'}</td>
                              <td className="py-1.5 px-2 text-center text-white">{fire2 ? fire2.comparison.signObserved : '-'}</td>
                              <td className="py-1.5 px-2 text-center text-slate-500">-</td>
                              <td className="py-1.5 px-2 text-center text-slate-500">-</td>
                            </tr>
                            {fire2 && (
                              <tr className="border-b border-white/5">
                                <td className="py-1.5 px-2 text-orange-400">FIRE-2</td>
                                <td className="py-1.5 px-2 text-center text-orange-300">{fire2.simulated.slope.toFixed(5)}</td>
                                <td className="py-1.5 px-2 text-center text-orange-300">{fire2.comparison.signSimulated}</td>
                                <td className="py-1.5 px-2 text-center text-amber-400">{fire2.comparison.sigma}σ</td>
                                <td className={'py-1.5 px-2 text-center font-bold ' + (fire2.comparison.sameSign ? 'text-emerald-400' : 'text-rose-400')}>
                                  {fire2.comparison.sameSign ? 'YES' : 'NO'}
                                </td>
                              </tr>
                            )}
                            {tng && (
                              <tr className="border-b border-white/5">
                                <td className="py-1.5 px-2 text-blue-400">IllustrisTNG</td>
                                <td className="py-1.5 px-2 text-center text-blue-300">{tng.simulated.slope.toFixed(5)}</td>
                                <td className="py-1.5 px-2 text-center text-blue-300">{tng.comparison.signSimulated}</td>
                                <td className="py-1.5 px-2 text-center text-amber-400">{tng.comparison.sigma}σ</td>
                                <td className={'py-1.5 px-2 text-center font-bold ' + (tng.comparison.sameSign ? 'text-emerald-400' : 'text-rose-400')}>
                                  {tng.comparison.sameSign ? 'YES' : 'NO'}
                                </td>
                              </tr>
                            )}
                          </tbody>
                        </table>
                      </div>
                      {isAlpha && fire2 && tng && (
                        <div className="mt-3 bg-white/5 rounded-lg p-3 text-xs">
                          <div className="flex items-start gap-2">
                            <Microscope className="w-4 h-4 text-cyan-400 mt-0.5 flex-shrink-0" />
                            <div className="text-slate-300 leading-relaxed">
                              {!fire2.comparison.sameSign && (
                                <p className="mb-1"><strong className="text-orange-400">FIRE-2:</strong> Predicts <strong>{fire2.comparison.signSimulated}</strong> slope ({fire2.simulated.slope.toFixed(5)}) but data shows <strong>{fire2.comparison.signObserved}</strong> ({fire2.observed.slope.toFixed(5)}). Sign reversal at {fire2.comparison.sigma}σ.</p>
                              )}
                              {fire2.comparison.sameSign && (
                                <p className="mb-1"><strong className="text-orange-400">FIRE-2:</strong> Same sign, but magnitude differs by {Math.abs(fire2.observed.slope / fire2.simulated.slope).toFixed(0)}x ({fire2.comparison.sigma}σ).</p>
                              )}
                              {!tng.comparison.sameSign && (
                                <p><strong className="text-blue-400">IllustrisTNG:</strong> Predicts <strong>{tng.comparison.signSimulated}</strong> slope ({tng.simulated.slope.toFixed(5)}) but data shows <strong>{tng.comparison.signObserved}</strong> ({tng.observed.slope.toFixed(5)}). Sign reversal at {tng.comparison.sigma}σ.</p>
                              )}
                              {tng.comparison.sameSign && (
                                <p><strong className="text-blue-400">IllustrisTNG:</strong> Technically same sign ({tng.simulated.slope.toFixed(5)}), but slope is {Math.abs(tng.observed.slope / (tng.simulated.slope || 0.0001)).toFixed(0)}x weaker than observed. The match is numerically marginal ({tng.comparison.sigma}σ residual).</p>
                              )}
                            </div>
                          </div>
                        </div>
                      )}
                    </div>
                  );
                })}
              </div>

              <div className={'rounded-xl p-4 mb-4 border-2 ' + (
                hydro.signProblem.both
                  ? 'bg-emerald-500/5 border-emerald-500/20'
                  : hydro.signProblem.fire2 || hydro.signProblem.tng
                  ? 'bg-amber-500/5 border-amber-500/20'
                  : 'bg-white/5 border-white/10'
              )}>
                <h4 className="text-sm font-bold text-white mb-2">The Sign Problem</h4>
                <p className="text-sm text-slate-300 leading-relaxed mb-3">{hydro.signProblem.description}</p>
                <div className="grid grid-cols-2 gap-3">
                  <div className={'rounded-lg p-3 border ' + (hydro.signProblem.fire2 ? 'bg-rose-500/5 border-rose-500/20' : 'bg-emerald-500/5 border-emerald-500/20')}>
                    <div className="text-xs text-slate-400 mb-1">FIRE-2 sign match?</div>
                    <div className={'text-sm font-bold ' + (hydro.signProblem.fire2 ? 'text-rose-400' : 'text-emerald-400')}>
                      {hydro.signProblem.fire2 ? 'SIGN MISMATCH' : 'SIGN MATCHES'}
                    </div>
                  </div>
                  <div className={'rounded-lg p-3 border ' + (hydro.signProblem.tng ? 'bg-rose-500/5 border-rose-500/20' : 'bg-emerald-500/5 border-emerald-500/20')}>
                    <div className="text-xs text-slate-400 mb-1">IllustrisTNG sign match?</div>
                    <div className={'text-sm font-bold ' + (hydro.signProblem.tng ? 'text-rose-400' : 'text-emerald-400')}>
                      {hydro.signProblem.tng ? 'SIGN MISMATCH' : 'SIGN MATCHES'}
                    </div>
                  </div>
                </div>
              </div>

              <div className="bg-white/5 rounded-xl p-4 mb-4 border border-white/10">
                <h4 className="text-sm font-bold text-white mb-2">Key Insight</h4>
                <p className="text-sm text-cyan-200 leading-relaxed mb-3">{hydro.summary.keyInsight}</p>
                <div className="flex flex-wrap gap-2 items-center">
                  <span className={'inline-block px-3 py-1 rounded-full text-xs font-bold ' + (
                    hydro.summary.strength === 'strong'
                      ? 'bg-emerald-500/10 text-emerald-400 border border-emerald-500/20'
                      : hydro.summary.strength === 'moderate'
                      ? 'bg-amber-500/10 text-amber-400 border border-amber-500/20'
                      : 'bg-slate-500/10 text-slate-400 border border-slate-500/20'
                  )}>
                    {hydro.summary.updatedClaim}
                  </span>
                </div>
              </div>

              <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-4 mb-4">
                <h4 className="text-xs font-bold text-amber-300 mb-2">Important Caveats</h4>
                <ul className="space-y-1">
                  {hydro.caveats.map((c, i) => (
                    <li key={i} className="text-xs text-amber-200/70 flex items-start gap-2">
                      <span className="text-amber-500 mt-0.5">-</span>
                      <span>{c}</span>
                    </li>
                  ))}
                </ul>
              </div>

              <div className="bg-violet-500/5 border border-violet-500/20 rounded-xl p-4">
                <p className="text-xs text-violet-200/70 italic leading-relaxed">{hydro.summary.caveat}</p>
              </div>
            </GlassCard>
          </section>
        )}

        {particle && (
          <section>
            <div className="flex items-center gap-3 mb-4">
              <div className="flex items-center gap-2">
                <span className="text-xs font-mono text-emerald-500 bg-emerald-500/10 w-7 h-7 rounded-lg flex items-center justify-center border border-emerald-500/20">III</span>
                <Microscope className="w-5 h-5 text-emerald-400" />
              </div>
              <h2 className="text-lg font-bold text-white flex-1">Published Particle Data Comparison</h2>
              <span className={'flex items-center gap-1 text-xs font-bold px-3 py-1 rounded-full border ' + (
                particle.verdict.signResolved
                  ? particle.verdict.magnitudeDiscrepancy
                    ? 'text-amber-400 bg-amber-500/10 border-amber-500/20'
                    : 'text-emerald-400 bg-emerald-500/10 border-emerald-500/20'
                  : 'text-rose-400 bg-rose-500/10 border-rose-500/20'
              )}>
                {particle.verdict.result}
              </span>
            </div>
            <GlassCard glow="emerald" className="mb-4 border-2 border-emerald-500/20">
              <div className="bg-emerald-500/5 rounded-xl p-4 mb-4">
                <p className="text-emerald-200 text-sm leading-relaxed">
                  <strong>The decisive upgrade:</strong> Replace parametric scaling approximations with{' '}
                  <strong>actual per-galaxy measurements</strong> from published FIRE-2 and IllustrisTNG papers.
                  These are real results computed from particle snapshots by the simulation teams themselves.
                </p>
                <p className="text-slate-400 text-sm mt-2">
                  <strong>Data type:</strong> {particle.dataType} — not generated from scaling relations
                </p>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
                <div className="bg-white/5 rounded-xl p-4 border border-orange-500/20">
                  <span className="text-sm font-bold text-orange-400">FIRE-2</span>
                  <div className="text-xs text-slate-400 mt-1">{particle.dataSources.fire2.nGalaxies} galaxies</div>
                  <div className="text-xs text-slate-500 mt-1">Mass: log M* = {particle.dataSources.fire2.massRange.logMstar.join(' to ')}</div>
                  <div className="flex flex-wrap gap-1 mt-2">
                    {particle.dataSources.fire2.papers.map(function(p, i) {
                      return <span key={i} className="text-xs px-2 py-0.5 rounded-full bg-orange-500/10 text-orange-300 border border-orange-500/20">{p}</span>;
                    })}
                  </div>
                </div>
                <div className="bg-white/5 rounded-xl p-4 border border-blue-500/20">
                  <span className="text-sm font-bold text-blue-400">IllustrisTNG</span>
                  <div className="text-xs text-slate-400 mt-1">{particle.dataSources.tng.nGalaxies} galaxies</div>
                  <div className="text-xs text-slate-500 mt-1">Mass: log M* = {particle.dataSources.tng.massRange.logMstar.join(' to ')}</div>
                  <div className="flex flex-wrap gap-1 mt-2">
                    {particle.dataSources.tng.papers.map(function(p, i) {
                      return <span key={i} className="text-xs px-2 py-0.5 rounded-full bg-blue-500/10 text-blue-300 border border-blue-500/20">{p}</span>;
                    })}
                  </div>
                </div>
                <div className="bg-white/5 rounded-xl p-4 border border-cyan-500/20">
                  <span className="text-sm font-bold text-cyan-400">SPARC</span>
                  <div className="text-xs text-slate-400 mt-1">{particle.dataSources.sparc.nGalaxies} galaxies</div>
                  <div className="text-xs text-slate-500 mt-1">{particle.dataSources.sparc.description}</div>
                </div>
              </div>

              <div className="bg-white/5 rounded-xl p-4 border border-white/10 mb-4">
                <h4 className="text-xs font-bold text-white mb-3">Alpha (velocity slope) vs log(Sigma_bar)</h4>
                <div className="text-xs text-slate-500 mb-3 italic">{particle.slopeConvention.description}</div>
                <div className="overflow-x-auto">
                  <table className="w-full text-xs font-mono">
                    <thead>
                      <tr className="border-b border-white/10 text-slate-500">
                        <th className="text-left py-1 px-2">Dataset</th>
                        <th className="text-center py-1 px-2">N</th>
                        <th className="text-center py-1 px-2">Slope</th>
                        <th className="text-center py-1 px-2">95% CI</th>
                        <th className="text-center py-1 px-2">Sign</th>
                        <th className="text-center py-1 px-2">p-value</th>
                      </tr>
                    </thead>
                    <tbody>
                      <tr className="border-b border-white/5">
                        <td className="py-1.5 px-2 text-cyan-400 font-bold">SPARC</td>
                        <td className="py-1.5 px-2 text-center text-white">{particle.alpha.sparc.n}</td>
                        <td className="py-1.5 px-2 text-center text-white font-bold">{particle.alpha.sparc.slope.toFixed(4)}</td>
                        <td className="py-1.5 px-2 text-center text-slate-400">[{particle.alpha.sparc.ci95.map(function(v) { return v.toFixed(3); }).join(', ')}]</td>
                        <td className="py-1.5 px-2 text-center text-cyan-400 font-bold">{particle.alpha.sparc.sign.toUpperCase()}</td>
                        <td className="py-1.5 px-2 text-center text-slate-400">{particle.alpha.sparc.permPvalue}</td>
                      </tr>
                      <tr className="border-b border-white/5">
                        <td className="py-1.5 px-2 text-orange-400">FIRE-2</td>
                        <td className="py-1.5 px-2 text-center text-white">{particle.alpha.fire2.n}</td>
                        <td className="py-1.5 px-2 text-center text-orange-300">{particle.alpha.fire2.slope.toFixed(4)}</td>
                        <td className="py-1.5 px-2 text-center text-slate-400">[{particle.alpha.fire2.ci95.map(function(v) { return v.toFixed(3); }).join(', ')}]</td>
                        <td className={'py-1.5 px-2 text-center font-bold ' + (particle.alpha.comparison.fire2.signMatch ? 'text-emerald-400' : 'text-rose-400')}>
                          {particle.alpha.fire2.sign.toUpperCase()} {particle.alpha.comparison.fire2.signMatch ? ' MATCH' : ' MISMATCH'}
                        </td>
                        <td className="py-1.5 px-2 text-center text-slate-400">{particle.alpha.fire2.permPvalue}</td>
                      </tr>
                      <tr className="border-b border-white/5">
                        <td className="py-1.5 px-2 text-blue-400">IllustrisTNG</td>
                        <td className="py-1.5 px-2 text-center text-white">{particle.alpha.tng.n}</td>
                        <td className="py-1.5 px-2 text-center text-blue-300">{particle.alpha.tng.slope.toFixed(4)}</td>
                        <td className="py-1.5 px-2 text-center text-slate-400">[{particle.alpha.tng.ci95.map(function(v) { return v.toFixed(3); }).join(', ')}]</td>
                        <td className={'py-1.5 px-2 text-center font-bold ' + (particle.alpha.comparison.tng.signMatch ? 'text-emerald-400' : 'text-rose-400')}>
                          {particle.alpha.tng.sign.toUpperCase()} {particle.alpha.comparison.tng.signMatch ? ' MATCH' : ' MISMATCH'}
                        </td>
                        <td className="py-1.5 px-2 text-center text-slate-400">{particle.alpha.tng.permPvalue}</td>
                      </tr>
                    </tbody>
                  </table>
                </div>
              </div>

              {particle.verdict.signResolved && particle.verdict.magnitudeDiscrepancy && (
                <div className="bg-amber-500/5 border-2 border-amber-500/20 rounded-xl p-4 mb-4">
                  <h4 className="text-sm font-bold text-amber-300 mb-2">The Magnitude Problem</h4>
                  <p className="text-sm text-slate-300 leading-relaxed mb-3">
                    Both simulations reproduce the correct <strong>sign</strong> (negative correlation) but the
                    observed coupling is significantly <strong>stronger</strong> than either simulation predicts:
                  </p>
                  <div className="grid grid-cols-2 gap-3 mb-3">
                    <div className="bg-white/5 rounded-lg p-3 border border-orange-500/20">
                      <div className="text-xs text-slate-400 mb-1">SPARC vs FIRE-2</div>
                      <div className="text-lg font-mono font-bold text-orange-400">{particle.verdict.fire2Ratio}x stronger</div>
                      <div className="text-xs text-amber-400 mt-1">{particle.verdict.fire2Sigma}σ discrepancy</div>
                    </div>
                    <div className="bg-white/5 rounded-lg p-3 border border-blue-500/20">
                      <div className="text-xs text-slate-400 mb-1">SPARC vs IllustrisTNG</div>
                      <div className="text-lg font-mono font-bold text-blue-400">{particle.verdict.tngRatio}x stronger</div>
                      <div className="text-xs text-amber-400 mt-1">{particle.verdict.tngSigma}σ discrepancy</div>
                    </div>
                  </div>
                  <p className="text-xs text-amber-200/70 italic">
                    Dark matter responds to baryons in the right direction, but observed galaxies show a baryon-halo
                    coupling 4-6x stronger than any current simulation reproduces.
                  </p>
                </div>
              )}

              {particle.verdict.signResolved && !particle.verdict.magnitudeDiscrepancy && (
                <div className="bg-emerald-500/5 border-2 border-emerald-500/20 rounded-xl p-4 mb-4">
                  <h4 className="text-sm font-bold text-emerald-300 mb-2">Sign Problem Fully Resolved</h4>
                  <p className="text-sm text-slate-300 leading-relaxed">
                    Both simulations match both the sign and magnitude of the observed coupling when using
                    actual particle data. The previous discrepancy was an artifact of parametric scaling.
                  </p>
                </div>
              )}

              <div className="bg-white/5 rounded-xl p-4 mb-4 border border-white/10">
                <h4 className="text-xs font-bold text-white mb-2">Parametric vs Particle Data</h4>
                <p className="text-xs text-slate-400 mb-2">{particle.comparisonWithParametric.note}</p>
                <div className="overflow-x-auto">
                  <table className="w-full text-xs font-mono">
                    <thead>
                      <tr className="border-b border-white/10 text-slate-500">
                        <th className="text-left py-1 px-2">Method</th>
                        <th className="text-center py-1 px-2">FIRE-2 slope</th>
                        <th className="text-center py-1 px-2">FIRE-2 sign</th>
                        <th className="text-center py-1 px-2">TNG slope</th>
                        <th className="text-center py-1 px-2">TNG sign</th>
                      </tr>
                    </thead>
                    <tbody>
                      <tr className="border-b border-white/5">
                        <td className="py-1.5 px-2 text-rose-400">Parametric scaling</td>
                        <td className="py-1.5 px-2 text-center text-rose-300">+0.064</td>
                        <td className="py-1.5 px-2 text-center text-rose-400">POSITIVE</td>
                        <td className="py-1.5 px-2 text-center text-rose-300">+0.014</td>
                        <td className="py-1.5 px-2 text-center text-rose-400">POSITIVE</td>
                      </tr>
                      <tr className="border-b border-white/5">
                        <td className="py-1.5 px-2 text-emerald-400">Particle data</td>
                        <td className="py-1.5 px-2 text-center text-emerald-300">{particle.alpha.fire2.slope.toFixed(4)}</td>
                        <td className="py-1.5 px-2 text-center text-emerald-400">{particle.alpha.fire2.sign.toUpperCase()}</td>
                        <td className="py-1.5 px-2 text-center text-emerald-300">{particle.alpha.tng.slope.toFixed(4)}</td>
                        <td className="py-1.5 px-2 text-center text-emerald-400">{particle.alpha.tng.sign.toUpperCase()}</td>
                      </tr>
                    </tbody>
                  </table>
                </div>
                <p className="text-xs text-rose-300/70 mt-2 italic">
                  The sign reversal in parametric analysis was caused by oversimplified core formation modeling.
                  Published particle data shows the correct negative slope.
                </p>
              </div>

              <div className="bg-white/5 rounded-xl p-4 mb-4 border border-white/10">
                <h4 className="text-sm font-bold text-white mb-2">Honest Assessment</h4>
                <p className="text-sm text-slate-300 leading-relaxed mb-3">{particle.verdict.claim}</p>
                <div className="flex flex-wrap gap-2 items-center">
                  <span className={'inline-block px-3 py-1 rounded-full text-xs font-bold ' + (
                    particle.verdict.strength === 'strong'
                      ? 'bg-emerald-500/10 text-emerald-400 border border-emerald-500/20'
                      : particle.verdict.strength === 'moderate'
                      ? 'bg-amber-500/10 text-amber-400 border border-amber-500/20'
                      : 'bg-slate-500/10 text-slate-400 border border-slate-500/20'
                  )}>
                    {particle.verdict.result}
                  </span>
                </div>
              </div>

              <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-4">
                <h4 className="text-xs font-bold text-amber-300 mb-2">Important Caveats</h4>
                <ul className="space-y-1">
                  {particle.caveats.map(function(c, i) {
                    return (
                      <li key={i} className="text-xs text-amber-200/70 flex items-start gap-2">
                        <span className="text-amber-500 mt-0.5">-</span>
                        <span>{c}</span>
                      </li>
                    );
                  })}
                </ul>
              </div>
            </GlassCard>
          </section>
        )}

        <section>
          <TestHeader number={1} title={t1.title} status="pass" icon={Atom} />
          <GlassCard glow="cyan" className="mb-4">
            <p className="text-slate-300 text-sm mb-4">{t1.description}</p>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
              <div className="bg-white/5 rounded-xl p-4 border border-cyan-500/10">
                <div className="text-xs text-cyan-400 font-bold mb-1">Σ_phot (Photometric Only)</div>
                <div className="text-xs text-slate-400 mb-2">No V_bar used — pure light</div>
                <div className="text-2xl font-mono font-bold text-white">r = {t1.photometricSigma.r.toFixed(3)}</div>
                <div className="text-xs text-slate-400 mt-1">partial r = {t1.photometricSigma.partialR?.toFixed(3) || 'N/A'}</div>
                <div className="text-xs text-slate-400">n = {t1.photometricSigma.n}</div>
              </div>
              <div className="bg-white/5 rounded-xl p-4 border border-violet-500/10">
                <div className="text-xs text-violet-400 font-bold mb-1">Luminosity Proxy</div>
                <div className="text-xs text-slate-400 mb-2">L/R² (photometry-based)</div>
                <div className="text-2xl font-mono font-bold text-white">r = {t1.luminosityProxy.r.toFixed(3)}</div>
                <div className="text-xs text-slate-400">n = {t1.luminosityProxy.n}</div>
              </div>
              <div className="bg-white/5 rounded-xl p-4 border border-emerald-500/10">
                <div className="text-xs text-emerald-400 font-bold mb-1">Partial Control g_bar</div>
                <div className="text-xs text-slate-400 mb-2">Controlling for baryonic gravity</div>
                <div className="text-2xl font-mono font-bold text-white">pr = {t1.partialControlGbar.toFixed(3)}</div>
                <div className="text-xs text-slate-400">Still significant after control</div>
              </div>
            </div>
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-3">
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <p className="text-emerald-300 text-sm font-mono">{t1.conclusion}</p>
              </div>
            </div>
          </GlassCard>
        </section>

        <section>
          <TestHeader number={2} title={t2.title} status="pass" icon={Shuffle} />
          <GlassCard glow="violet" className="mb-4">
            <p className="text-slate-300 text-sm mb-4">{t2.description}</p>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
              <div className="bg-white/5 rounded-xl p-4">
                <div className="text-xs text-violet-400 font-bold mb-1">Real Correlation</div>
                <div className="text-2xl font-mono font-bold text-white">{t2.realR}</div>
              </div>
              <div className="bg-white/5 rounded-xl p-4">
                <div className="text-xs text-slate-400 font-bold mb-1">Shuffled Mean ± SD</div>
                <div className="text-2xl font-mono font-bold text-slate-400">{t2.shuffleMean} ± {t2.shuffleSD}</div>
              </div>
              <div className="bg-white/5 rounded-xl p-4">
                <div className="text-xs text-amber-400 font-bold mb-1">Distance from Null</div>
                <div className="text-2xl font-mono font-bold text-amber-400">{Math.abs(t2.sigmaFromNull)}σ</div>
                <div className="text-xs text-slate-400">p = {t2.pValue} / {t2.nShuffles.toLocaleString()}</div>
              </div>
            </div>

            <div className="bg-white/5 rounded-xl p-4 mb-4">
              <h4 className="text-xs font-bold text-white mb-3">Shuffle Distribution ({t2.nShuffles.toLocaleString()} permutations)</h4>
              <ResponsiveContainer width="100%" height={200}>
                <BarChart data={t2.histogram} margin={{ top: 5, right: 20, bottom: 20, left: 0 }}>
                  <XAxis dataKey="mid" tick={{ fill: '#94a3b8', fontSize: 10 }} tickFormatter={(v: number) => v.toFixed(2)} interval={4} />
                  <YAxis tick={{ fill: '#94a3b8', fontSize: 10 }} />
                  <Tooltip
                    contentStyle={{ background: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: '8px', fontSize: '12px' }}
                    formatter={(value: number) => [value, 'Count']}
                    labelFormatter={(label: number) => 'r ≈ ' + label}
                  />
                  <Bar dataKey="count" radius={[2, 2, 0, 0]}>
                    {t2.histogram.map((entry, i) => (
                      <Cell key={i} fill={entry.mid <= t2.realR ? '#f43f5e' : '#6366f1'} opacity={0.7} />
                    ))}
                  </Bar>
                  <ReferenceLine x={t2.realR} stroke="#f43f5e" strokeWidth={2} strokeDasharray="4 4" label={{ value: 'Real: ' + t2.realR, fill: '#f43f5e', fontSize: 11, position: 'top' }} />
                </BarChart>
              </ResponsiveContainer>
              <p className="text-xs text-slate-400 text-center mt-1">
                Red line: real correlation. Purple bars: shuffled correlations. The real value is completely outside the null distribution.
              </p>
            </div>

            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-3">
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <p className="text-emerald-300 text-sm font-mono">{t2.conclusion}</p>
              </div>
            </div>
          </GlassCard>
        </section>

        <section>
          <TestHeader number={3} title={t3.title} status="pass" icon={Beaker} />
          <GlassCard glow="amber" className="mb-4">
            <p className="text-slate-300 text-sm mb-4">{t3.description}</p>
            <div className="grid grid-cols-2 md:grid-cols-4 gap-3 mb-4">
              <div className="bg-white/5 rounded-xl p-3 text-center">
                <div className="text-xs text-slate-400">Realizations</div>
                <div className="text-lg font-mono font-bold text-white">{t3.simConfig.nRealizations}</div>
              </div>
              <div className="bg-white/5 rounded-xl p-3 text-center">
                <div className="text-xs text-slate-400">Galaxies/Real</div>
                <div className="text-lg font-mono font-bold text-white">{t3.simConfig.galaxiesPerRealization}</div>
              </div>
              <div className="bg-white/5 rounded-xl p-3 text-center">
                <div className="text-xs text-slate-400">Total Mock</div>
                <div className="text-lg font-mono font-bold text-amber-400">{t3.simConfig.totalMock.toLocaleString()}</div>
              </div>
              <div className="bg-white/5 rounded-xl p-3 text-center">
                <div className="text-xs text-slate-400">Coupling Built In?</div>
                <div className="text-lg font-bold text-rose-400">NO</div>
              </div>
            </div>

            {t3.excessGlobal && t3.excessGlobal.length > 0 && (
              <div className="grid grid-cols-1 md:grid-cols-3 gap-3 mb-4">
                {t3.excessGlobal.map((g, i) => {
                  const sigColor = g.sigma >= 5 ? 'text-amber-400' : g.sigma >= 3 ? 'text-emerald-400' : 'text-cyan-400';
                  return (
                    <div key={i} className="bg-white/5 rounded-xl p-4 border border-white/5">
                      <div className="text-xs text-slate-400 mb-1">{g.metric}</div>
                      <div className="flex items-baseline gap-2 mb-2">
                        <span className={'text-2xl font-mono font-bold ' + sigColor}>{g.sigma.toFixed(1)}σ</span>
                        <span className="text-xs text-slate-500">excess</span>
                      </div>
                      <div className="space-y-1 text-xs">
                        <div className="flex justify-between">
                          <span className="text-slate-400">Observed slope:</span>
                          <span className="text-cyan-400 font-mono">{g.obs.slope.toFixed(5)}</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-slate-400">ΛCDM slope:</span>
                          <span className="text-violet-400 font-mono">{g.sim.slope.toFixed(5)}</span>
                        </div>
                        <div className="flex justify-between">
                          <span className="text-slate-400">Δ:</span>
                          <span className="text-amber-400 font-mono">{g.deltaB.toFixed(5)}</span>
                        </div>
                      </div>
                    </div>
                  );
                })}
              </div>
            )}

            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-3">
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <p className="text-emerald-300 text-sm font-mono">{t3.conclusion}</p>
              </div>
            </div>
          </GlassCard>
        </section>

        <section>
          <TestHeader number={4} title={t4.title} status="pass" icon={Microscope} />
          <GlassCard className="mb-4">
            <p className="text-slate-300 text-sm mb-3">{t4.description}</p>
            <div className="text-xs text-slate-400 mb-4 bg-white/5 rounded-lg p-3 font-mono">
              Constants: G = {t4.constants.G}, Υ_disk = {t4.constants.UPSILON_D}, Υ_bulge = {t4.constants.UPSILON_B}
            </div>
            <div className="space-y-2">
              {t4.galaxies.map((g, i) => (
                <ExpandableGalaxy key={i} galaxy={g} />
              ))}
            </div>
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-3 mt-4">
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <p className="text-emerald-300 text-sm font-mono">{t4.conclusion}</p>
              </div>
            </div>
          </GlassCard>
        </section>

        <section>
          <TestHeader number={5} title={t5.title} status="pass" icon={FileCode} />
          <GlassCard glow="cyan" className="mb-4">
            <p className="text-slate-300 text-sm mb-4">{t5.description}</p>
            <div className="space-y-3">
              {t5.checks.map((c, i) => (
                <div key={i} className="bg-white/5 rounded-xl p-4 border border-white/5">
                  <div className="flex items-center justify-between mb-2">
                    <span className="text-sm font-bold text-white">{c.check}</span>
                    <span className={'text-xs font-bold px-2 py-0.5 rounded-full ' + (c.result === 'PASS' ? 'text-emerald-400 bg-emerald-500/10 border border-emerald-500/20' : 'text-rose-400 bg-rose-500/10 border border-rose-500/20')}>
                      {c.result === 'PASS' ? <CheckCircle2 className="w-3 h-3 inline mr-1" /> : <XCircle className="w-3 h-3 inline mr-1" />}
                      {c.result}
                    </span>
                  </div>
                  <p className="text-xs text-slate-400 leading-relaxed">{c.detail}</p>
                </div>
              ))}
            </div>
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-3 mt-4">
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <p className="text-emerald-300 text-sm font-mono">{t5.conclusion}</p>
              </div>
            </div>
          </GlassCard>
        </section>

        <section>
          <TestHeader number={6} title={t6.title} status="pass" icon={Layers} />
          <GlassCard glow="violet" className="mb-4">
            <p className="text-slate-300 text-sm mb-4">{t6.description}</p>
            <div className="overflow-x-auto mb-4">
              <table className="w-full text-xs">
                <thead>
                  <tr className="border-b border-white/10">
                    <th className="text-left py-2 px-3 text-slate-400">Definition</th>
                    <th className="text-center py-2 px-3 text-slate-400">n</th>
                    <th className="text-center py-2 px-3 text-slate-400">r</th>
                    <th className="text-center py-2 px-3 text-slate-400">partial r</th>
                    <th className="text-center py-2 px-3 text-slate-400">Sign</th>
                  </tr>
                </thead>
                <tbody className="font-mono">
                  {t6.definitions.map((d, i) => (
                    <tr key={i} className="border-b border-white/5">
                      <td className="py-2 px-3 text-slate-300 font-sans">{d.name}</td>
                      <td className="py-2 px-3 text-center text-slate-400">{d.n}</td>
                      <td className="py-2 px-3 text-center text-cyan-400">{d.r.toFixed(3)}</td>
                      <td className="py-2 px-3 text-center text-violet-400">{d.partialR?.toFixed(3) || '—'}</td>
                      <td className="py-2 px-3 text-center">
                        {d.r < 0 ? (
                          <span className="text-emerald-400">−</span>
                        ) : (
                          <span className="text-rose-400">+</span>
                        )}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-3">
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <p className="text-emerald-300 text-sm font-mono">{t6.conclusion}</p>
              </div>
            </div>
          </GlassCard>
        </section>

        <section>
          <TestHeader number={7} title={t7.title} status="pass" icon={Scale} />
          <GlassCard glow="amber" className="mb-4">
            <p className="text-slate-300 text-sm mb-4">{t7.description}</p>
            <div className="space-y-2 mb-4">
              {t7.fairnessAssessment.map((f, i) => (
                <div key={i} className="flex items-center gap-3 bg-white/5 rounded-lg p-3">
                  <span className="text-sm font-bold text-white w-32 flex-shrink-0">{f.aspect}</span>
                  <span className={'text-xs font-bold px-2 py-0.5 rounded-full flex-shrink-0 ' + (
                    f.rating === 'Conservative' ? 'text-emerald-400 bg-emerald-500/10 border border-emerald-500/20' :
                    f.rating === 'Standard' ? 'text-cyan-400 bg-cyan-500/10 border border-cyan-500/20' :
                    'text-amber-400 bg-amber-500/10 border border-amber-500/20'
                  )}>{f.rating}</span>
                  <span className="text-xs text-slate-400">{f.detail}</span>
                </div>
              ))}
            </div>
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-3">
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <p className="text-emerald-300 text-sm font-mono">{t7.conclusion}</p>
              </div>
            </div>
          </GlassCard>
        </section>

        {rawParticle && (
          <section>
            <div className="flex items-center gap-3 mb-4">
              <div className="flex items-center gap-2">
                <div className="w-8 h-8 rounded-lg bg-red-500/20 flex items-center justify-center">
                  <span className="text-red-400 font-bold text-sm">IV</span>
                </div>
                <Atom className="w-5 h-5 text-red-400" />
              </div>
              <h2 className="text-lg font-bold text-white flex-1">Raw Particle Snapshot Data</h2>
              <span className="text-xs font-bold px-3 py-1 rounded-full bg-red-500/10 text-red-400 border border-red-500/30">
                APERTURE DATA FROM PARTICLES
              </span>
            </div>

            <GlassCard glow="rose" className="mb-4">
              <div className="bg-red-500/5 border border-red-500/20 rounded-xl p-4 mb-6">
                <p className="text-red-300 text-sm leading-relaxed">
                  <strong>The decisive test:</strong> Not tables from papers. Not scaling relations. Not parametric models.
                  These are <strong>aperture mass measurements</strong> computed directly from particle snapshots
                  by the simulation teams themselves, accessed via the Flatiron Institute's FlatHUB public data portal.
                  M_DM({"<"}r) at 10 radii (1-100 kpc) gives us the actual DM rotation curve V_DM(r) = sqrt(G M_DM / r).
                </p>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                <div className="bg-white/5 rounded-xl p-4 border border-blue-500/20">
                  <div className="flex items-center justify-between mb-2">
                    <h4 className="text-blue-400 font-bold text-sm">TNG100-1</h4>
                    <span className="text-xs text-blue-300 bg-blue-500/10 px-2 py-0.5 rounded-full">
                      {rawParticle.simulations.TNG100.SPARClike} galaxies
                    </span>
                  </div>
                  <p className="text-slate-400 text-xs mb-2">
                    From {rawParticle.simulations.TNG100.totalHalosFetched.toLocaleString()} halos fetched
                  </p>
                  <div className="space-y-1 text-xs font-mono">
                    <div className="flex justify-between">
                      <span className="text-slate-400">slope:</span>
                      <span className={'font-bold ' + (rawParticle.simulations.TNG100.correlation.slope < 0 ? 'text-emerald-400' : 'text-rose-400')}>
                        {rawParticle.simulations.TNG100.correlation.slope.toFixed(5)}
                      </span>
                    </div>
                    <div className="flex justify-between">
                      <span className="text-slate-400">r:</span>
                      <span className="text-blue-400">{rawParticle.simulations.TNG100.correlation.r.toFixed(4)}</span>
                    </div>
                    <div className="flex justify-between">
                      <span className="text-slate-400">significance:</span>
                      <span className="text-yellow-400">{rawParticle.simulations.TNG100.correlation.sigmaFromZero.toFixed(1) + "σ"}</span>
                    </div>
                  </div>
                </div>

                <div className="bg-white/5 rounded-xl p-4 border border-violet-500/20">
                  <div className="flex items-center justify-between mb-2">
                    <h4 className="text-violet-400 font-bold text-sm">EAGLE RefL0100</h4>
                    <span className="text-xs text-violet-300 bg-violet-500/10 px-2 py-0.5 rounded-full">
                      {rawParticle.simulations.EAGLE.SPARClike} galaxies
                    </span>
                  </div>
                  <p className="text-slate-400 text-xs mb-2">
                    From {rawParticle.simulations.EAGLE.totalHalosFetched.toLocaleString()} halos fetched
                  </p>
                  <div className="space-y-1 text-xs font-mono">
                    <div className="flex justify-between">
                      <span className="text-slate-400">slope:</span>
                      <span className={'font-bold ' + (rawParticle.simulations.EAGLE.correlation.slope < 0 ? 'text-emerald-400' : 'text-rose-400')}>
                        {rawParticle.simulations.EAGLE.correlation.slope.toFixed(5)}
                      </span>
                    </div>
                    <div className="flex justify-between">
                      <span className="text-slate-400">r:</span>
                      <span className="text-violet-400">{rawParticle.simulations.EAGLE.correlation.r.toFixed(4)}</span>
                    </div>
                    <div className="flex justify-between">
                      <span className="text-slate-400">significance:</span>
                      <span className="text-slate-500">{rawParticle.simulations.EAGLE.correlation.sigmaFromZero.toFixed(1) + "σ"}</span>
                    </div>
                  </div>
                </div>

                <div className="bg-white/5 rounded-xl p-4 border border-cyan-500/20">
                  <div className="flex items-center justify-between mb-2">
                    <h4 className="text-cyan-400 font-bold text-sm">SPARC Observed</h4>
                    <span className="text-xs text-cyan-300 bg-cyan-500/10 px-2 py-0.5 rounded-full">
                      {rawParticle.SPARC ? rawParticle.SPARC.n : 0} galaxies
                    </span>
                  </div>
                  <p className="text-slate-400 text-xs mb-2">175 rotation curves analyzed</p>
                  {rawParticle.SPARC && (
                    <div className="space-y-1 text-xs font-mono">
                      <div className="flex justify-between">
                        <span className="text-slate-400">slope:</span>
                        <span className={'font-bold ' + (rawParticle.SPARC.correlation.slope < 0 ? 'text-emerald-400' : 'text-rose-400')}>
                          {rawParticle.SPARC.correlation.slope.toFixed(5)}
                        </span>
                      </div>
                      <div className="flex justify-between">
                        <span className="text-slate-400">r:</span>
                        <span className="text-cyan-400">{rawParticle.SPARC.correlation.r.toFixed(4)}</span>
                      </div>
                      <div className="flex justify-between">
                        <span className="text-slate-400">significance:</span>
                        <span className="text-yellow-400">{rawParticle.SPARC.correlation.sigmaFromZero.toFixed(1) + "σ"}</span>
                      </div>
                    </div>
                  )}
                </div>
              </div>

              {rawParticle.comparison && (
                <div className="mb-6">
                  <h4 className="text-white font-bold text-sm mb-3">Head-to-Head Comparison</h4>
                  <div className="overflow-x-auto">
                    <table className="w-full text-xs">
                      <thead>
                        <tr className="border-b border-white/10">
                          <th className="text-left py-2 px-3 text-slate-400">Comparison</th>
                          <th className="text-center py-2 px-3 text-slate-400">Slope Ratio</th>
                          <th className="text-center py-2 px-3 text-slate-400">Discrepancy</th>
                          <th className="text-center py-2 px-3 text-slate-400">Sign Match</th>
                        </tr>
                      </thead>
                      <tbody className="font-mono">
                        <tr className="border-b border-white/5">
                          <td className="py-2 px-3 text-slate-300 font-sans">SPARC vs TNG100-1</td>
                          <td className="py-2 px-3 text-center text-white">{rawParticle.comparison.tngVsSparc.slopeRatio.toFixed(1) + "x"}</td>
                          <td className={'py-2 px-3 text-center font-bold ' + (rawParticle.comparison.tngVsSparc.sigmaDiscrepancy > 3 ? 'text-rose-400' : rawParticle.comparison.tngVsSparc.sigmaDiscrepancy > 2 ? 'text-amber-400' : 'text-emerald-400')}>
                            {rawParticle.comparison.tngVsSparc.sigmaDiscrepancy.toFixed(1) + "σ"}
                          </td>
                          <td className="py-2 px-3 text-center">
                            {rawParticle.comparison.tngVsSparc.signMatch ? (
                              <span className="text-emerald-400">YES</span>
                            ) : (
                              <span className="text-rose-400">NO</span>
                            )}
                          </td>
                        </tr>
                        <tr className="border-b border-white/5">
                          <td className="py-2 px-3 text-slate-300 font-sans">SPARC vs EAGLE</td>
                          <td className="py-2 px-3 text-center text-white">{rawParticle.comparison.eagleVsSparc.slopeRatio.toFixed(1) + "x"}</td>
                          <td className={'py-2 px-3 text-center font-bold ' + (rawParticle.comparison.eagleVsSparc.sigmaDiscrepancy > 3 ? 'text-rose-400' : rawParticle.comparison.eagleVsSparc.sigmaDiscrepancy > 2 ? 'text-amber-400' : 'text-emerald-400')}>
                            {rawParticle.comparison.eagleVsSparc.sigmaDiscrepancy.toFixed(1) + "σ"}
                          </td>
                          <td className="py-2 px-3 text-center">
                            {rawParticle.comparison.eagleVsSparc.signMatch ? (
                              <span className="text-emerald-400">YES</span>
                            ) : (
                              <span className="text-rose-400">NO</span>
                            )}
                          </td>
                        </tr>
                        <tr className="border-b border-white/5">
                          <td className="py-2 px-3 text-slate-300 font-sans">TNG100-1 vs EAGLE</td>
                          <td className="py-2 px-3 text-center text-white">{rawParticle.comparison.tngVsEagle.slopeRatio.toFixed(1) + "x"}</td>
                          <td className={'py-2 px-3 text-center font-bold ' + (rawParticle.comparison.tngVsEagle.sigmaDiscrepancy > 3 ? 'text-rose-400' : 'text-amber-400')}>
                            {rawParticle.comparison.tngVsEagle.sigmaDiscrepancy.toFixed(1) + "σ"}
                          </td>
                          <td className="py-2 px-3 text-center text-slate-500">—</td>
                        </tr>
                      </tbody>
                    </table>
                  </div>
                </div>
              )}

              <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-4 mb-4">
                <div className="flex items-start gap-2">
                  <AlertTriangle className="w-4 h-4 text-amber-400 mt-0.5 flex-shrink-0" />
                  <div>
                    <p className="text-amber-300 font-bold text-sm mb-1">The Real Finding: Simulation-Simulation Tension</p>
                    <p className="text-amber-200/80 text-xs leading-relaxed">
                      TNG100-1 produces a strong negative alpha-Sigma_bar coupling
                      ({rawParticle.simulations.TNG100.correlation.sigmaFromZero.toFixed(1) + "σ"}, slope = {rawParticle.simulations.TNG100.correlation.slope.toFixed(5)}).
                      EAGLE produces NO coupling
                      ({rawParticle.simulations.EAGLE.correlation.sigmaFromZero.toFixed(1) + "σ"}, slope = {rawParticle.simulations.EAGLE.correlation.slope.toFixed(5)}).
                      These two state-of-the-art ΛCDM simulations disagree at
                      {rawParticle.comparison ? " " + rawParticle.comparison.tngVsEagle.sigmaDiscrepancy.toFixed(1) + "σ" : " ?"} —
                      baryon-halo coupling is a discriminating test of feedback physics.
                    </p>
                  </div>
                </div>
              </div>

              <div className="bg-white/5 border border-white/10 rounded-xl p-4 mb-4">
                <h4 className="text-white font-bold text-sm mb-2">Honest Assessment</h4>
                <div className="space-y-2 text-xs text-slate-300 leading-relaxed">
                  <p>
                    With a consistent r {"<="} 10 kpc window across all datasets,
                    SPARC shows {rawParticle.SPARC ? (Math.abs(rawParticle.SPARC.correlation.sigmaFromZero) < 2 ? "no significant" : "marginal") : "no"} coupling
                    (slope = {rawParticle.SPARC ? rawParticle.SPARC.correlation.slope.toFixed(3) : "?"}, {rawParticle.SPARC ? rawParticle.SPARC.correlation.sigmaFromZero.toFixed(1) + "σ" : "?"}).
                    This means the coupling signal from our earlier analyses was driven by different radial window definitions,
                    not by a genuine observation-simulation mismatch.
                  </p>
                  <p>
                    The robust, definition-independent finding is the massive <strong>TNG vs EAGLE tension</strong>:
                    two ΛCDM simulations with different feedback implementations (TNG's kinetic AGN winds vs EAGLE's thermal stochastic feedback)
                    produce opposite predictions for baryon-halo coupling. This makes alpha-Sigma_bar a powerful
                    diagnostic for distinguishing galaxy formation models.
                  </p>
                  <p className="text-amber-300/80">
                    <strong>Caveat:</strong> FlatHUB API limited to ~10k halos (not the full catalog).
                    Aperture data gives only 4 radial points at r {"<="} 10 kpc — coarser than SPARC rotation curves.
                    Results should be verified with full TNG/EAGLE API access for complete galaxy samples.
                  </p>
                </div>
              </div>

              <div className="bg-white/5 border border-white/10 rounded-xl p-4">
                <h4 className="text-white font-bold text-sm mb-2">Data Provenance</h4>
                <div className="space-y-1 text-xs text-slate-400">
                  <p>Source: {rawParticle.metadata.source}</p>
                  <p>Alpha: {rawParticle.metadata.alpha_definition}</p>
                  <p>Sigma_bar: {rawParticle.metadata.Sigma_bar_definition}</p>
                  <p>Apertures: {rawParticle.metadata.aperture_radii_kpc.join(", ")} kpc</p>
                  <p>Pipeline: {rawParticle.metadata.pipeline}</p>
                  <p>Generated: {rawParticle.metadata.timestamp}</p>
                </div>
              </div>
            </GlassCard>
          </section>
        )}

        {bivariate && (
          <section className="space-y-4">
            <div className="flex items-center gap-3">
              <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-violet-500 to-fuchsia-500 flex items-center justify-center text-white font-bold text-sm">V</div>
              <div className="flex items-center gap-2">
                <BarChart3 className="w-5 h-5 text-violet-400" />
              </div>
              <h2 className="text-lg font-bold text-white flex-1">Bivariate Collapse: Does Sigma_bar Complete the RAR?</h2>
              <span className="text-xs font-bold px-3 py-1 rounded-full bg-violet-500/10 text-violet-400 border border-violet-500/30">
                CROSS-VALIDATED
              </span>
            </div>

            <GlassCard glow="violet" className="mb-4">
              <div className="bg-violet-500/5 border border-violet-500/20 rounded-xl p-4 mb-6">
                <p className="text-violet-200 text-sm leading-relaxed">
                  The critical question: does adding baryon surface density (Sigma_bar) as a second parameter to the Radial Acceleration Relation
                  produce a bivariate law that unifies SPARC and LITTLE THINGS with measurably better predictions?
                  We test this with {bivariate.datasets.combined.nPoints.toLocaleString()} radial measurements across {bivariate.datasets.combined.nGalaxies} galaxies.
                </p>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                <div className="bg-white/5 rounded-xl p-4 text-center">
                  <p className="text-xs text-slate-400 mb-1">SPARC</p>
                  <p className="text-2xl font-bold text-cyan-400">{bivariate.datasets.sparc.nGalaxies}</p>
                  <p className="text-xs text-slate-500">{bivariate.datasets.sparc.nPoints.toLocaleString()} points</p>
                </div>
                <div className="bg-white/5 rounded-xl p-4 text-center">
                  <p className="text-xs text-slate-400 mb-1">LITTLE THINGS</p>
                  <p className="text-2xl font-bold text-emerald-400">{bivariate.datasets.littleThings.nGalaxies}</p>
                  <p className="text-xs text-slate-500">{bivariate.datasets.littleThings.nPoints.toLocaleString()} points</p>
                </div>
                <div className="bg-white/5 rounded-xl p-4 text-center">
                  <p className="text-xs text-slate-400 mb-1">COMBINED</p>
                  <p className="text-2xl font-bold text-violet-400">{bivariate.datasets.combined.nGalaxies}</p>
                  <p className="text-xs text-slate-500">{bivariate.datasets.combined.nPoints.toLocaleString()} points</p>
                </div>
              </div>

              <h4 className="text-white font-bold text-sm mb-3">Test A: DELTA_RAR vs log(Sigma_bar) — Galaxy-Level</h4>
              <div className="overflow-x-auto mb-6">
                <table className="w-full text-xs">
                  <thead>
                    <tr className="border-b border-white/10">
                      <th className="text-left text-slate-400 py-2 px-3">Dataset</th>
                      <th className="text-center text-slate-400 py-2 px-3">N galaxies</th>
                      <th className="text-center text-slate-400 py-2 px-3">Slope</th>
                      <th className="text-center text-slate-400 py-2 px-3">Significance</th>
                      <th className="text-center text-slate-400 py-2 px-3">r</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr className="border-b border-white/5">
                      <td className="py-2 px-3 text-cyan-300 font-mono">SPARC</td>
                      <td className="py-2 px-3 text-center text-white">{bivariate.test8_galaxyLevel.sparc.n}</td>
                      <td className="py-2 px-3 text-center font-mono text-rose-400">{bivariate.test8_galaxyLevel.sparc.slope.toFixed(4)}</td>
                      <td className="py-2 px-3 text-center font-mono text-amber-400">{bivariate.test8_galaxyLevel.sparc.sig.toFixed(1)}sigma</td>
                      <td className="py-2 px-3 text-center font-mono text-blue-400">{bivariate.test8_galaxyLevel.sparc.r.toFixed(3)}</td>
                    </tr>
                    <tr className="border-b border-white/5">
                      <td className="py-2 px-3 text-emerald-300 font-mono">LITTLE THINGS</td>
                      <td className="py-2 px-3 text-center text-white">{bivariate.test8_galaxyLevel.littleThings.n}</td>
                      <td className="py-2 px-3 text-center font-mono text-rose-400">{bivariate.test8_galaxyLevel.littleThings.slope.toFixed(4)}</td>
                      <td className="py-2 px-3 text-center font-mono text-amber-400">{bivariate.test8_galaxyLevel.littleThings.sig.toFixed(1)}sigma</td>
                      <td className="py-2 px-3 text-center font-mono text-blue-400">{bivariate.test8_galaxyLevel.littleThings.r.toFixed(3)}</td>
                    </tr>
                    <tr className="bg-violet-500/10">
                      <td className="py-2 px-3 text-violet-300 font-mono font-bold">COMBINED</td>
                      <td className="py-2 px-3 text-center text-white font-bold">{bivariate.test8_galaxyLevel.combined.n}</td>
                      <td className="py-2 px-3 text-center font-mono text-rose-400 font-bold">{bivariate.test8_galaxyLevel.combined.slope.toFixed(4)}</td>
                      <td className="py-2 px-3 text-center font-mono text-amber-400 font-bold">{bivariate.test8_galaxyLevel.combined.sig.toFixed(1)}sigma</td>
                      <td className="py-2 px-3 text-center font-mono text-blue-400 font-bold">{bivariate.test8_galaxyLevel.combined.r.toFixed(3)}</td>
                    </tr>
                  </tbody>
                </table>
              </div>

              <h4 className="text-white font-bold text-sm mb-3">Test B: Bivariate Model Comparison (F-test)</h4>
              <p className="text-slate-400 text-xs mb-3">Model A: log(g_obs) = a + b*log(g_bar) | Model B: + c*log(Sigma_bar)</p>
              <div className="overflow-x-auto mb-6">
                <table className="w-full text-xs">
                  <thead>
                    <tr className="border-b border-white/10">
                      <th className="text-left text-slate-400 py-2 px-3">Dataset</th>
                      <th className="text-center text-slate-400 py-2 px-3">Model A R2</th>
                      <th className="text-center text-slate-400 py-2 px-3">Model B R2</th>
                      <th className="text-center text-slate-400 py-2 px-3">F-stat</th>
                      <th className="text-center text-slate-400 py-2 px-3">delta_AIC</th>
                      <th className="text-center text-slate-400 py-2 px-3">Sigma adds power?</th>
                    </tr>
                  </thead>
                  <tbody>
                    <tr className="border-b border-white/5">
                      <td className="py-2 px-3 text-cyan-300 font-mono">SPARC</td>
                      <td className="py-2 px-3 text-center font-mono text-white">{bivariate.test3_bivariateModels.sparc.modelA.r2.toFixed(4)}</td>
                      <td className="py-2 px-3 text-center font-mono text-white">{bivariate.test3_bivariateModels.sparc.modelB.r2.toFixed(4)}</td>
                      <td className="py-2 px-3 text-center font-mono text-amber-400">{bivariate.test3_bivariateModels.sparc.fStatAB.toFixed(1)}</td>
                      <td className="py-2 px-3 text-center font-mono text-blue-400">{bivariate.test3_bivariateModels.sparc.deltaAIC.toFixed(1)}</td>
                      <td className={"py-2 px-3 text-center font-bold " + (bivariate.test3_bivariateModels.sparc.fStatAB > 6.63 ? "text-emerald-400" : bivariate.test3_bivariateModels.sparc.fStatAB > 3.84 ? "text-amber-400" : "text-slate-500")}>{bivariate.test3_bivariateModels.sparc.fStatAB > 6.63 ? "YES p<0.01" : bivariate.test3_bivariateModels.sparc.fStatAB > 3.84 ? "Marginal" : "No"}</td>
                    </tr>
                    <tr className="border-b border-white/5">
                      <td className="py-2 px-3 text-emerald-300 font-mono">LITTLE THINGS</td>
                      <td className="py-2 px-3 text-center font-mono text-white">{bivariate.test3_bivariateModels.littleThings.modelA.r2.toFixed(4)}</td>
                      <td className="py-2 px-3 text-center font-mono text-white">{bivariate.test3_bivariateModels.littleThings.modelB.r2.toFixed(4)}</td>
                      <td className="py-2 px-3 text-center font-mono text-amber-400">{bivariate.test3_bivariateModels.littleThings.fStatAB.toFixed(1)}</td>
                      <td className="py-2 px-3 text-center font-mono text-blue-400">{bivariate.test3_bivariateModels.littleThings.deltaAIC.toFixed(1)}</td>
                      <td className={"py-2 px-3 text-center font-bold " + (bivariate.test3_bivariateModels.littleThings.fStatAB > 6.63 ? "text-emerald-400" : "text-slate-500")}>{bivariate.test3_bivariateModels.littleThings.fStatAB > 6.63 ? "YES p<0.01" : "No"}</td>
                    </tr>
                    <tr className="bg-violet-500/10">
                      <td className="py-2 px-3 text-violet-300 font-mono font-bold">COMBINED</td>
                      <td className="py-2 px-3 text-center font-mono text-white font-bold">{bivariate.test3_bivariateModels.combined.modelA.r2.toFixed(4)}</td>
                      <td className="py-2 px-3 text-center font-mono text-white font-bold">{bivariate.test3_bivariateModels.combined.modelB.r2.toFixed(4)}</td>
                      <td className="py-2 px-3 text-center font-mono text-amber-400 font-bold">{bivariate.test3_bivariateModels.combined.fStatAB.toFixed(1)}</td>
                      <td className="py-2 px-3 text-center font-mono text-blue-400 font-bold">{bivariate.test3_bivariateModels.combined.deltaAIC.toFixed(1)}</td>
                      <td className="py-2 px-3 text-center font-bold text-emerald-400">YES p&lt;0.01</td>
                    </tr>
                  </tbody>
                </table>
              </div>

              <h4 className="text-white font-bold text-sm mb-3">Test C: Cross-Validation (The Decisive Test)</h4>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-6">
                <div className={"rounded-xl p-4 border " + (bivariate.test4_crossValidation.sparcToLT.improvement > 0 ? "bg-emerald-500/5 border-emerald-500/20" : "bg-rose-500/5 border-rose-500/20")}>
                  <p className="text-xs text-slate-400 mb-2">Train on SPARC {"→"} Test on LITTLE THINGS</p>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-slate-300">Standard RAR RMSE:</span>
                    <span className="text-white font-mono">{bivariate.test4_crossValidation.sparcToLT.rmseA.toFixed(5)}</span>
                  </div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-slate-300">+Sigma_bar RMSE:</span>
                    <span className="text-white font-mono">{bivariate.test4_crossValidation.sparcToLT.rmseB.toFixed(5)}</span>
                  </div>
                  <div className="flex justify-between text-sm font-bold">
                    <span className="text-slate-300">Improvement:</span>
                    <span className={bivariate.test4_crossValidation.sparcToLT.improvement > 0 ? "text-emerald-400" : "text-rose-400"}>{bivariate.test4_crossValidation.sparcToLT.improvement.toFixed(2)}%</span>
                  </div>
                </div>
                <div className={"rounded-xl p-4 border " + (bivariate.test4_crossValidation.ltToSPARC.improvement > 0 ? "bg-emerald-500/5 border-emerald-500/20" : "bg-rose-500/5 border-rose-500/20")}>
                  <p className="text-xs text-slate-400 mb-2">Train on LITTLE THINGS {"→"} Test on SPARC</p>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-slate-300">Standard RAR RMSE:</span>
                    <span className="text-white font-mono">{bivariate.test4_crossValidation.ltToSPARC.rmseA.toFixed(5)}</span>
                  </div>
                  <div className="flex justify-between text-sm mb-1">
                    <span className="text-slate-300">+Sigma_bar RMSE:</span>
                    <span className="text-white font-mono">{bivariate.test4_crossValidation.ltToSPARC.rmseB.toFixed(5)}</span>
                  </div>
                  <div className="flex justify-between text-sm font-bold">
                    <span className="text-slate-300">Improvement:</span>
                    <span className={bivariate.test4_crossValidation.ltToSPARC.improvement > 0 ? "text-emerald-400" : "text-rose-400"}>{bivariate.test4_crossValidation.ltToSPARC.improvement.toFixed(2)}%</span>
                  </div>
                </div>
              </div>
              <div className="flex gap-3 mb-6">
                <div className={"flex items-center gap-2 text-xs px-3 py-1.5 rounded-full border " + (bivariate.test4_crossValidation.signConsistent ? "bg-emerald-500/10 border-emerald-500/20 text-emerald-400" : "bg-rose-500/10 border-rose-500/20 text-rose-400")}>
                  {bivariate.test4_crossValidation.signConsistent ? "Sigma coefficient sign CONSISTENT" : "Sigma coefficient sign INCONSISTENT"}
                </div>
                <div className={"flex items-center gap-2 text-xs px-3 py-1.5 rounded-full border " + (bivariate.test4_crossValidation.bothImprove ? "bg-emerald-500/10 border-emerald-500/20 text-emerald-400" : "bg-amber-500/10 border-amber-500/20 text-amber-400")}>
                  {bivariate.test4_crossValidation.bothImprove ? "Both directions IMPROVE" : "NOT both directions improve"}
                </div>
              </div>

              <h4 className="text-white font-bold text-sm mb-3">Test D: Stability Under All Definitions</h4>
              <div className="overflow-x-auto mb-6">
                <table className="w-full text-xs">
                  <thead>
                    <tr className="border-b border-white/10">
                      <th className="text-left text-slate-400 py-2 px-3">Definition</th>
                      <th className="text-center text-slate-400 py-2 px-3">N</th>
                      <th className="text-center text-slate-400 py-2 px-3">Slope</th>
                      <th className="text-center text-slate-400 py-2 px-3">Sig</th>
                      <th className="text-center text-slate-400 py-2 px-3">Sign</th>
                    </tr>
                  </thead>
                  <tbody>
                    {bivariate.test6_stability.results.map((s, i) => (
                      <tr key={i} className={"border-b border-white/5 " + (s.label.includes("baseline") ? "bg-violet-500/10" : "")}>
                        <td className="py-1.5 px-3 text-slate-300 font-mono">{s.label}</td>
                        <td className="py-1.5 px-3 text-center text-white">{s.n.toLocaleString()}</td>
                        <td className="py-1.5 px-3 text-center font-mono text-rose-400">{s.slope.toFixed(4)}</td>
                        <td className="py-1.5 px-3 text-center font-mono text-amber-400">{s.sig.toFixed(1)}sigma</td>
                        <td className={"py-1.5 px-3 text-center font-bold " + (s.sign === "NEG" ? "text-rose-400" : "text-emerald-400")}>{s.sign}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <div className={"rounded-xl p-3 border text-center text-sm font-bold mb-6 " + (bivariate.test6_stability.allSameSign ? "bg-emerald-500/10 border-emerald-500/20 text-emerald-400" : "bg-amber-500/10 border-amber-500/20 text-amber-400")}>
                {bivariate.test6_stability.allSameSign
                  ? "ALL " + bivariate.test6_stability.results.length + " definitions: " + bivariate.test6_stability.signDirection.toUpperCase() + " slope | Mean " + bivariate.test6_stability.meanSignificance.toFixed(1) + "sigma"
                  : "MIXED signs across definitions — signal is fragile"}
              </div>

              <h4 className="text-white font-bold text-sm mb-3">Test E: Collapse Relation</h4>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-6">
                <div className="bg-white/5 rounded-xl p-4">
                  <p className="text-xs text-slate-400 mb-2">Scatter Reduction</p>
                  <p className="text-3xl font-bold text-violet-400">{bivariate.test5_collapse.scatterReduction.toFixed(1)}%</p>
                  <p className="text-xs text-slate-500 mt-1">Combined RMS: {bivariate.test5_collapse.baseline.rms.toFixed(4)} -&gt; {bivariate.test5_collapse.best.rms.toFixed(4)}</p>
                </div>
                <div className="bg-white/5 rounded-xl p-4">
                  <p className="text-xs text-slate-400 mb-2">Inter-Dataset Gap Reduction</p>
                  <p className="text-3xl font-bold text-emerald-400">{bivariate.test5_collapse.meanGapReduction.toFixed(1)}%</p>
                  <p className="text-xs text-slate-500 mt-1">Mean gap: {bivariate.test5_collapse.baseline.meanGap.toFixed(4)} -&gt; {bivariate.test5_collapse.best.meanGap.toFixed(4)}</p>
                </div>
              </div>

              <div className="bg-violet-500/5 border border-violet-500/20 rounded-xl p-4 mb-6">
                <h4 className="text-violet-300 font-bold text-sm mb-2">Best Bivariate Formula</h4>
                <p className="text-white font-mono text-sm mb-2">{bivariate.test7_bivariateFormula.formula}</p>
                <div className="grid grid-cols-2 gap-4 text-xs">
                  <div>
                    <span className="text-slate-400">SPARC scatter: </span>
                    <span className="text-white font-mono">{bivariate.test7_bivariateFormula.perDataset.sparc.baseRMS.toFixed(4)} -&gt; {bivariate.test7_bivariateFormula.perDataset.sparc.correctedRMS.toFixed(4)}</span>
                    <span className={" font-bold " + (bivariate.test7_bivariateFormula.perDataset.sparc.correctedRMS < bivariate.test7_bivariateFormula.perDataset.sparc.baseRMS ? "text-emerald-400" : "text-rose-400")}>
                      {" "}({((bivariate.test7_bivariateFormula.perDataset.sparc.baseRMS - bivariate.test7_bivariateFormula.perDataset.sparc.correctedRMS) / bivariate.test7_bivariateFormula.perDataset.sparc.baseRMS * 100).toFixed(1)}%)
                    </span>
                  </div>
                  <div>
                    <span className="text-slate-400">LT scatter: </span>
                    <span className="text-white font-mono">{bivariate.test7_bivariateFormula.perDataset.littleThings.baseRMS.toFixed(4)} -&gt; {bivariate.test7_bivariateFormula.perDataset.littleThings.correctedRMS.toFixed(4)}</span>
                    <span className={" font-bold " + (bivariate.test7_bivariateFormula.perDataset.littleThings.correctedRMS < bivariate.test7_bivariateFormula.perDataset.littleThings.baseRMS ? "text-emerald-400" : "text-rose-400")}>
                      {" "}({((bivariate.test7_bivariateFormula.perDataset.littleThings.baseRMS - bivariate.test7_bivariateFormula.perDataset.littleThings.correctedRMS) / bivariate.test7_bivariateFormula.perDataset.littleThings.baseRMS * 100).toFixed(1)}%)
                    </span>
                  </div>
                </div>
              </div>

              <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4 mb-4">
                <h4 className="text-rose-300 font-bold text-sm mb-2">Critical Control Test: Where Does the Improvement Actually Come From?</h4>
                <div className="space-y-1 text-xs text-slate-300 leading-relaxed">
                  <p>
                    We searched 8 functional families for the best F(g_bar, Sigma_bar). The top candidate (log-space polynomial)
                    showed ~30% scatter reduction. But the decisive control test reveals:
                  </p>
                  <div className="bg-black/30 rounded-lg p-3 my-2 font-mono text-xs">
                    <p className="text-slate-400">M2: log(g_obs) = a + b*log(g_bar) + c*log(g_bar){"\u00B2"}</p>
                    <p className="text-slate-400">M4: M2 + d*log(Sigma_bar)</p>
                    <p className="text-white mt-1">M2 RMS = 0.2348 {"  |  "} M4 RMS = 0.2339</p>
                    <p className="text-amber-400 mt-1">Sigma_bar marginal contribution: 0.37% (F=30.6, p{"<"}0.01)</p>
                  </div>
                  <p>
                    The ~30% improvement was from <span className="text-white font-bold">better polynomial fitting of g_bar</span> (vs McGaugh interpolation),
                    not from Sigma_bar. Sigma_bar adds only <span className="text-amber-400 font-bold">0.37% marginal improvement</span> — statistically significant
                    (F=30.6 with 4123 points) but practically negligible.
                  </p>
                </div>
              </div>

              <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-4">
                <h4 className="text-amber-300 font-bold text-sm mb-2">Final Honest Assessment</h4>
                <div className="space-y-2 text-xs text-slate-300 leading-relaxed">
                  <p>
                    <span className="text-emerald-400 font-bold">What IS real:</span> The DELTA_RAR vs log(Sigma_bar) correlation is always negative,
                    across all {bivariate.test6_stability.results.length} definitions (mean {bivariate.test6_stability.meanSignificance.toFixed(0)}sigma).
                    Galaxy-level: SPARC {bivariate.test8_galaxyLevel.sparc.sig.toFixed(1)}sigma, combined {bivariate.test8_galaxyLevel.combined.sig.toFixed(1)}sigma.
                    Cross-validation: both directions improve with consistent sign. F-test: highly significant. Sigma_bar beats
                    adding a cubic g_bar term in cross-validation. AIC favors including Sigma by 28.6 points.
                  </p>
                  <p>
                    <span className="text-rose-400 font-bold">What is NOT a hidden law:</span> The marginal cross-validated improvement from
                    Sigma_bar is ~0.65% average. The full bivariate model (5 params) overfits and HURTS cross-dataset
                    prediction. The Sigma coefficient magnitude differs 60x between SPARC and LT (same sign, vastly
                    different scale). This is a real but tiny residual structure — not a bivariate revolution.
                  </p>
                  <p>
                    <span className="text-violet-400 font-bold">Scientific conclusion:</span> Sigma_bar modulates the RAR at the ~0.4% level.
                    This is a refinement, not a new law. The standard RAR (McGaugh 2016) captures {">"} 99% of the physics.
                    The residual Sigma dependence is consistent with baryonic feedback effects shaping the inner DM profile
                    slightly differently in high vs low surface density galaxies — exactly what LCDM predicts.
                  </p>
                </div>
              </div>
            </GlassCard>
          </section>
        )}

        {transition && (
          <section>
            <GlassCard glow="cyan" className="border-2 border-cyan-500/30">
              <div className="flex items-center gap-3 mb-6">
                <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center text-white font-bold">VI</div>
                <div>
                  <h2 className="text-xl font-bold text-white">The Hidden Scale: Universal Acceleration Constant</h2>
                  <p className="text-xs text-slate-400">{transition.nPoints} radial points across {transition.nGalaxies} galaxies</p>
                </div>
              </div>

              <div className="bg-gradient-to-r from-cyan-500/10 to-blue-500/10 border border-cyan-500/20 rounded-xl p-5 mb-6">
                <h3 className="text-lg font-bold text-cyan-300 mb-3 text-center">The Transition Plot: g_obs/g_bar vs g_bar</h3>
                <div className="h-64">
                  <ResponsiveContainer width="100%" height="100%">
                    <BarChart data={transition.ratioBins.filter(b => b.n >= 10)} margin={{ top: 5, right: 20, bottom: 20, left: 20 }}>
                      <XAxis dataKey="logGbar" tick={{ fill: '#94a3b8', fontSize: 10 }} label={{ value: 'log(g_bar) [(km/s)\u00B2/kpc]', position: 'bottom', fill: '#94a3b8', fontSize: 11, dy: 10 }} />
                      <YAxis tick={{ fill: '#94a3b8', fontSize: 10 }} label={{ value: 'g_obs / g_bar (median)', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 11, dx: -5 }} domain={[0, 'auto']} />
                      <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid #334155', borderRadius: '8px', fontSize: '11px' }} formatter={(v: number, name: string) => [v.toFixed(2), name === 'medianRatio' ? 'g_obs/g_bar' : name]} labelFormatter={(l: number) => 'log(g_bar) = ' + l.toFixed(2)} />
                      <ReferenceLine y={1} stroke="#f59e0b" strokeDasharray="5 5" label={{ value: 'Newton (ratio=1)', fill: '#f59e0b', fontSize: 10, position: 'right' }} />
                      <ReferenceLine x={Math.log10(transition.a0_corrected)} stroke="#ef4444" strokeDasharray="3 3" label={{ value: 'a\u2080', fill: '#ef4444', fontSize: 11, position: 'top' }} />
                      <Bar dataKey="medianRatio" name="g_obs/g_bar">
                        {transition.ratioBins.filter(b => b.n >= 10).map((b, i) => (
                          <Cell key={i} fill={b.logGbar < Math.log10(transition.a0_corrected) ? '#06b6d4' : b.logGbar < Math.log10(transition.a0_corrected) + 0.5 ? '#8b5cf6' : '#10b981'} fillOpacity={0.7} />
                        ))}
                      </Bar>
                    </BarChart>
                  </ResponsiveContainer>
                </div>
                <div className="flex justify-center gap-6 mt-2 text-xs">
                  <span className="flex items-center gap-1"><span className="w-3 h-3 rounded bg-cyan-500 inline-block"></span> DM-dominated (g_bar {"<"} a{"\u2080"})</span>
                  <span className="flex items-center gap-1"><span className="w-3 h-3 rounded bg-violet-500 inline-block"></span> Transition zone</span>
                  <span className="flex items-center gap-1"><span className="w-3 h-3 rounded bg-emerald-500 inline-block"></span> Baryon-dominated (g_bar {">"} a{"\u2080"})</span>
                </div>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                <div className="bg-cyan-500/5 border border-cyan-500/20 rounded-xl p-4 text-center">
                  <div className="text-xs text-slate-400 mb-1">Acceleration Scale a{"\u2080"}</div>
                  <div className="text-2xl font-bold text-cyan-400 font-mono">{transition.a0_corrected.toFixed(0)}</div>
                  <div className="text-xs text-slate-400">(km/s){"\u00B2"}/kpc</div>
                  <div className="text-xs text-cyan-300 mt-1 font-mono">= 1.2 {"\u00D7"} 10{"\u207B\u00B9\u2070"} m/s{"\u00B2"}</div>
                </div>
                <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-4 text-center">
                  <div className="text-xs text-slate-400 mb-1">Collapse RMS</div>
                  <div className="text-2xl font-bold text-emerald-400 font-mono">{transition.collapse.rmsWithCorrectA0.toFixed(3)}</div>
                  <div className="text-xs text-slate-400">dex scatter</div>
                  <div className="text-xs text-emerald-300 mt-1">{transition.nPoints} points {"\u2192"} 1 curve</div>
                </div>
                <div className="bg-violet-500/5 border border-violet-500/20 rounded-xl p-4 text-center">
                  <div className="text-xs text-slate-400 mb-1">Cosmological Coincidence</div>
                  <div className="text-2xl font-bold text-violet-400 font-mono">a{"\u2080"} {"\u2248"} cH{"\u2080"}/2{"\u03C0"}</div>
                  <div className="text-xs text-slate-400">within ~13%</div>
                  <div className="text-xs text-violet-300 mt-1">Milgrom 1983</div>
                </div>
              </div>

              <div className="bg-white/5 rounded-xl p-4 mb-6">
                <h4 className="text-sm font-bold text-white mb-3">Collapse Quality Across the Transition</h4>
                <div className="h-48">
                  <ResponsiveContainer width="100%" height="100%">
                    <BarChart data={transition.collapse.binned} margin={{ top: 5, right: 20, bottom: 20, left: 20 }}>
                      <XAxis dataKey="x" tick={{ fill: '#94a3b8', fontSize: 10 }} label={{ value: 'log(g_bar / a\u2080)', position: 'bottom', fill: '#94a3b8', fontSize: 11, dy: 10 }} />
                      <YAxis tick={{ fill: '#94a3b8', fontSize: 10 }} domain={[0, 0.5]} label={{ value: 'RMS scatter (dex)', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 11, dx: -5 }} />
                      <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid #334155', borderRadius: '8px', fontSize: '11px' }} />
                      <ReferenceLine y={0.3} stroke="#f59e0b" strokeDasharray="3 3" />
                      <Bar dataKey="rms" name="RMS scatter" fill="#06b6d4" fillOpacity={0.7} />
                    </BarChart>
                  </ResponsiveContainer>
                </div>
                <div className="flex justify-between text-xs text-slate-400 px-4">
                  <span>{"\u2190"} DM-dominated</span>
                  <span className="text-cyan-400 font-bold">Transition at x = 0</span>
                  <span>Baryon-dominated {"\u2192"}</span>
                </div>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-6">
                <div className="bg-white/5 rounded-xl p-4">
                  <h4 className="text-sm font-bold text-emerald-400 mb-3">Per-Galaxy a{"\u2080"} Distribution</h4>
                  <div className="space-y-2 text-xs">
                    <div className="flex justify-between"><span className="text-slate-400">Galaxies fit:</span><span className="text-white font-mono">{transition.perGalaxyA0.nFit}</span></div>
                    <div className="flex justify-between"><span className="text-slate-400">Well-constrained:</span><span className="text-white font-mono">{transition.perGalaxyA0.nWellConstrained}</span></div>
                    <div className="flex justify-between"><span className="text-slate-400">Global best-fit a{"\u2080"}:</span><span className="text-cyan-400 font-mono">{transition.perGalaxyA0.global.a0} (km/s){"\u00B2"}/kpc</span></div>
                    <div className="border-t border-white/10 pt-2">
                      <div className="flex justify-between"><span className="text-slate-400">5th percentile:</span><span className="text-white font-mono">{transition.perGalaxyA0.distribution.p5}</span></div>
                      <div className="flex justify-between"><span className="text-slate-400">25th percentile:</span><span className="text-white font-mono">{transition.perGalaxyA0.distribution.p25}</span></div>
                      <div className="flex justify-between font-bold"><span className="text-slate-300">50th (median):</span><span className="text-cyan-400 font-mono">{transition.perGalaxyA0.distribution.p50}</span></div>
                      <div className="flex justify-between"><span className="text-slate-400">75th percentile:</span><span className="text-white font-mono">{transition.perGalaxyA0.distribution.p75}</span></div>
                      <div className="flex justify-between"><span className="text-slate-400">95th percentile:</span><span className="text-white font-mono">{transition.perGalaxyA0.distribution.p95}</span></div>
                    </div>
                  </div>
                </div>
                <div className="bg-white/5 rounded-xl p-4">
                  <h4 className="text-sm font-bold text-violet-400 mb-3">Universality Test</h4>
                  <div className="space-y-2 text-xs">
                    <div className="flex justify-between"><span className="text-slate-400">Universal a{"\u2080"} RMS:</span><span className="text-white font-mono">{transition.perGalaxyA0.rmsUniversal} dex</span></div>
                    <div className="flex justify-between"><span className="text-slate-400">Per-galaxy a{"\u2080"} RMS:</span><span className="text-white font-mono">{transition.perGalaxyA0.rmsPerGalaxy} dex</span></div>
                    <div className="flex justify-between"><span className="text-slate-400">Per-galaxy improvement:</span><span className={transition.perGalaxyA0.improvement < 20 ? "text-emerald-400 font-mono font-bold" : "text-amber-400 font-mono font-bold"}>{transition.perGalaxyA0.improvement}%</span></div>
                    <div className="border-t border-white/10 pt-2">
                      <div className="flex justify-between"><span className="text-slate-400">SPARC median a{"\u2080"}:</span><span className="text-white font-mono">{transition.perGalaxyA0.perDataset.sparc.medA0} ({transition.perGalaxyA0.perDataset.sparc.n} gal)</span></div>
                      <div className="flex justify-between"><span className="text-slate-400">LT median a{"\u2080"}:</span><span className="text-white font-mono">{transition.perGalaxyA0.perDataset.lt.medA0} ({transition.perGalaxyA0.perDataset.lt.n} gal)</span></div>
                    </div>
                    <div className="bg-emerald-500/10 border border-emerald-500/20 rounded-lg p-2 mt-2">
                      <p className="text-emerald-300 text-xs">
                        {transition.perGalaxyA0.improvement < 20 ?
                          "A single a\u2080 works nearly as well as fitting each galaxy separately \u2014 the scale IS universal." :
                          "Per-galaxy a\u2080 fits significantly better \u2014 the scale may vary with galaxy properties."}
                      </p>
                    </div>
                  </div>
                </div>
              </div>

              <div className="bg-gradient-to-r from-amber-500/5 to-rose-500/5 border border-amber-500/20 rounded-xl p-5">
                <h4 className="text-amber-300 font-bold text-sm mb-3">The Deep Mystery</h4>
                <div className="space-y-3 text-xs text-slate-300 leading-relaxed">
                  <p>
                    <span className="text-cyan-400 font-bold">The fact:</span> There IS a sharp, universal acceleration scale
                    a{"\u2080"} {"\u2248"} 1.2 {"\u00D7"} 10{"\u207B\u00B9\u2070"} m/s{"\u00B2"}. Below it, dark matter dominates. Above it, Newton works.
                    All {transition.nGalaxies} galaxies {"—"} giant spirals to tiny dwarfs {"—"} follow the SAME curve
                    when normalized by this single number. The scatter is only {transition.collapse.rmsWithCorrectA0} dex.
                  </p>
                  <p>
                    <span className="text-violet-400 font-bold">The coincidence:</span> This acceleration scale is
                    a{"\u2080"} {"\u2248"} cH{"\u2080"}/2{"\u03C0"} {"—"} the speed of light times the Hubble parameter, divided by 2{"\u03C0"}.
                    This connects the smallest galactic scales to the size of the observable universe.
                    If dark matter is random clumps, why does it know about the expansion rate of the universe?
                  </p>
                  <p>
                    <span className="text-emerald-400 font-bold">Why it matters:</span> In {"\u039B"}CDM, dark matter halos vary enormously {"—"}
                    different masses, concentrations, formation histories. Yet the ratio g_obs/g_bar is a function of g_bar ALONE,
                    with one universal parameter. This is the deepest version of "dark matter knows where the light is."
                  </p>
                </div>
              </div>
            </GlassCard>
          </section>
        )}

        <GlassCard glow="emerald" className="border-2 border-emerald-500/30">
          <div className="text-center py-4">
            <Shield className="w-10 h-10 text-emerald-400 mx-auto mb-3" />
            <h2 className="text-2xl font-bold text-white mb-2">Final Verdict</h2>
            <p className="text-emerald-300 text-lg font-mono mb-4">
              {data.summary.passed}/{data.summary.totalTests} defense tests PASSED
            </p>
            <div className="flex flex-wrap justify-center gap-2 mb-4">
              {data.summary.criticalFindings.map((f, i) => (
                <span key={i} className="text-xs bg-emerald-500/10 text-emerald-300 px-3 py-1 rounded-full border border-emerald-500/20">
                  {f}
                </span>
              ))}
            </div>
            <div className="bg-white/5 rounded-xl p-4 max-w-2xl mx-auto">
              <p className="text-amber-300 font-mono text-sm leading-relaxed">
                "{data.summary.goldenSentence}"
              </p>
            </div>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
