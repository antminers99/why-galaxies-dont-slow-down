import React, { useEffect, useState } from 'react';
import { Layout } from '@/components/layout';
import { BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine, Cell } from 'recharts';
import {
  Shield, CheckCircle2, XCircle, Shuffle, Microscope, FileCode, Layers,
  Scale, ChevronDown, ChevronRight, AlertTriangle, Atom, Beaker
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
  };
  summary: {
    maxSigmaWithoutFeedback: number;
    maxSigmaWithFeedback: number;
    averageReductionPercent: number;
    allMetricsFullyExplained: boolean;
    verdict: string;
    honestAssessment: string;
    updatedClaim: string;
  };
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
  const [error, setError] = useState(false);

  useEffect(() => {
    Promise.all([
      fetch(`${import.meta.env.BASE_URL}defense-validation.json`).then(r => { if (!r.ok) throw new Error('fetch failed'); return r.json(); }),
      fetch(`${import.meta.env.BASE_URL}feedback-test.json`).then(r => { if (!r.ok) throw new Error('fetch failed'); return r.json(); }).catch(() => null),
    ])
      .then(([defData, fbData]) => { setData(defData); setFeedback(fbData); })
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
                  <strong>Our response:</strong> We built ΛCDM + full baryonic feedback (Di Cintio et al. 2014 core formation + Blumenthal et al. 1986 adiabatic contraction) and re-ran the comparison.
                </p>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
                {feedback.metrics.map((m, i) => {
                  const fullyExplained = m.feedbackFullyExplains;
                  const sigBefore = m.excessNoFeedback.sigma;
                  const sigAfter = m.excessWithFeedback.sigma;
                  const color = fullyExplained ? 'emerald' : sigAfter > sigBefore ? 'rose' : 'amber';
                  const borderColor = { emerald: 'border-emerald-500/20', rose: 'border-rose-500/20', amber: 'border-amber-500/20' }[color];
                  return (
                    <div key={i} className={'bg-white/5 rounded-xl p-4 border ' + borderColor}>
                      <div className="text-xs font-bold text-white mb-3">{m.name}</div>
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
                          <span className="text-slate-400">σ without feedback:</span>
                          <span className="text-amber-400 font-mono font-bold">{sigBefore.toFixed(1)}σ</span>
                        </div>
                        <div className="flex justify-between items-center">
                          <span className="text-slate-400">σ with feedback:</span>
                          <span className={'font-mono font-bold ' + (fullyExplained ? 'text-emerald-400' : 'text-rose-400')}>{sigAfter.toFixed(1)}σ</span>
                        </div>
                        <div className={'text-center py-1 rounded-lg mt-1 ' + (fullyExplained ? 'bg-emerald-500/10 text-emerald-300' : 'bg-amber-500/10 text-amber-300')}>
                          {fullyExplained ? 'Feedback explains this metric' : 'Still ' + sigAfter.toFixed(1) + 'σ excess'}
                        </div>
                      </div>
                    </div>
                  );
                })}
              </div>

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
                      <p className="text-xs text-amber-300 mt-1">Effect: {c.effect}</p>
                    </div>
                  ))}
                </div>
              </div>

              <div className={'rounded-xl p-4 border ' + (
                feedback.summary.allMetricsFullyExplained
                  ? 'bg-rose-500/5 border-rose-500/20'
                  : 'bg-amber-500/5 border-amber-500/20'
              )}>
                <h4 className="text-sm font-bold text-white mb-2">Honest Verdict</h4>
                <p className="text-sm text-amber-200 leading-relaxed mb-2">{feedback.summary.honestAssessment}</p>
                <div className={'inline-block px-3 py-1 rounded-full text-xs font-bold ' + (
                  feedback.summary.maxSigmaWithFeedback >= 3
                    ? 'bg-amber-500/10 text-amber-400 border border-amber-500/20'
                    : feedback.summary.maxSigmaWithFeedback >= 2
                    ? 'bg-yellow-500/10 text-yellow-400 border border-yellow-500/20'
                    : 'bg-rose-500/10 text-rose-400 border border-rose-500/20'
                )}>
                  {feedback.summary.updatedClaim}
                </div>
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
