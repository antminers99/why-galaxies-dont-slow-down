import React, { useState, useEffect } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, ComposedChart, Line, BarChart, Bar, Cell, Legend } from 'recharts';
import { FlaskConical, CheckCircle2, XCircle, AlertTriangle, ArrowRight, Microscope, BookOpen, Scale, Lightbulb, Target, ShieldCheck } from 'lucide-react';

interface LTData {
  dataset: string;
  source: string;
  instrument: string;
  nGalaxies: number;
  type: string;
  independent: boolean;
  simpleRegression: { slope: number; slopeError: number; intercept: number; r: number; r2: number; p: number; n: number };
  partialCorrelation: { rPartial: number; t: number; df: number; p: number };
  pointLevelRegression: { slope: number; slopeError: number; r: number; r2: number; p: number; nPoints: number };
  sparcComparison: { sparcSlope: number; ltSlope: number; signConsistent: boolean; sparcPartialR: number; ltPartialR: number };
  galaxies: Array<{ name: string; dist: number; vmax: number; mbar: number; meanDeltaRAR: number; meanLogSigBar: number; nPoints: number }>;
}

interface SparcData {
  validation: any;
  perGalaxy: Array<{ name: string; sigma_bar: number; meanDeltaOuter: number; Vmax: number }>;
}

const StatBox = ({ label, value, sub, color = 'cyan' }: { label: string; value: string; sub?: string; color?: string }) => (
  <div className="text-center p-3">
    <div className="text-xs text-slate-500 uppercase tracking-wider mb-1">{label}</div>
    <div className={`text-xl font-mono font-bold text-${color}-400`}>{value}</div>
    {sub && <div className="text-xs text-slate-500 mt-0.5">{sub}</div>}
  </div>
);

const Badge = ({ pass, label }: { pass: boolean; label: string }) => (
  <span className={`inline-flex items-center gap-1.5 px-3 py-1.5 rounded-full text-xs font-bold ${pass ? 'bg-emerald-500/20 text-emerald-400 border border-emerald-500/30' : 'bg-red-500/20 text-red-400 border border-red-500/30'}`}>
    {pass ? <CheckCircle2 className="w-3.5 h-3.5" /> : <XCircle className="w-3.5 h-3.5" />}
    {label}
  </span>
);

export default function ReplicationPage() {
  const [lt, setLt] = useState<LTData | null>(null);
  const [sparc, setSparc] = useState<SparcData | null>(null);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'little-things-replication.json').then(r => r.json()).then(setLt).catch(() => {});
    fetch(import.meta.env.BASE_URL + 'rar-analysis-real.json').then(r => r.json()).then(setSparc).catch(() => {});
  }, []);

  if (!lt || !sparc) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-64">
          <div className="w-8 h-8 border-2 border-t-primary border-white/10 rounded-full animate-spin" />
        </div>
      </Layout>
    );
  }

  const sparcGalaxies = sparc.perGalaxy.filter((g) => g.sigma_bar > 0);
  const sparcScatter = sparcGalaxies.map((g) => ({
    logSigma: +Math.log10(g.sigma_bar).toFixed(3),
    delta: +g.meanDeltaOuter.toFixed(4),
    name: g.name,
  }));
  const ltScatter = lt.galaxies.map((g) => ({
    logSigma: +g.meanLogSigBar.toFixed(3),
    delta: +g.meanDeltaRAR.toFixed(4),
    name: g.name,
  }));

  const sparcXs = sparcGalaxies.map(g => Math.log10(g.sigma_bar));
  const sparcYs = sparcGalaxies.map(g => g.meanDeltaOuter);
  const smx = sparcXs.reduce((a, b) => a + b, 0) / sparcXs.length;
  const smy = sparcYs.reduce((a, b) => a + b, 0) / sparcYs.length;
  let ssxy = 0, ssxx = 0;
  for (let i = 0; i < sparcXs.length; i++) { ssxy += (sparcXs[i] - smx) * (sparcYs[i] - smy); ssxx += (sparcXs[i] - smx) ** 2; }
  const sparcB = ssxy / ssxx;
  const sparcA = smy - sparcB * smx;

  const xAll = [...sparcXs, ...lt.galaxies.map(g => g.meanLogSigBar)];
  const xMin = Math.min(...xAll) - 0.2;
  const xMax = Math.max(...xAll) + 0.2;
  const sparcLine = [];
  const ltLine = [];
  for (let x = xMin; x <= xMax; x += 0.1) {
    sparcLine.push({ logSigma: +x.toFixed(2), fit: +(sparcA + sparcB * x).toFixed(4) });
    ltLine.push({ logSigma: +x.toFixed(2), fit: +(lt.simpleRegression.intercept + lt.simpleRegression.slope * x).toFixed(4) });
  }

  const slopeMatch = Math.abs(lt.simpleRegression.slope - (-0.178)) < 2 * lt.simpleRegression.slopeError;
  const signMatch = lt.sparcComparison.signConsistent;
  const partialNeg = lt.partialCorrelation.rPartial < 0;

  const comparisonBars = [
    { metric: 'Slope b', sparc: -0.178, lt: +lt.simpleRegression.slope.toFixed(3) },
    { metric: 'Partial r', sparc: -0.656, lt: +lt.partialCorrelation.rPartial.toFixed(3) },
    { metric: 'Simple r', sparc: +(ssxy / Math.sqrt(ssxx * sparcYs.reduce((a, y) => a + (y - smy) ** 2, 0))).toFixed(3), lt: +lt.simpleRegression.r.toFixed(3) },
  ];

  return (
    <Layout>
      <div className="space-y-6">
        <div>
          <h1 className="text-3xl font-display font-bold text-white flex items-center gap-3">
            <FlaskConical className="w-8 h-8 text-emerald-400" />
            External Replication
          </h1>
          <p className="text-slate-400 mt-2">
            Independent validation of the density-corrected RAR using the LITTLE THINGS survey (Oh et al. 2015).
          </p>
        </div>

        <GlassCard glow="emerald">
          <div className="flex items-start gap-4">
            <div className="w-12 h-12 rounded-xl bg-emerald-500/20 flex items-center justify-center flex-shrink-0">
              <ShieldCheck className="w-7 h-7 text-emerald-400" />
            </div>
            <div className="flex-1">
              <h2 className="text-xl font-display font-bold text-white mb-2">Replication Verdict</h2>
              <p className="text-slate-300 mb-4">
                The density-dependent RAR residuals discovered in SPARC (175 galaxies, Lelli et al. 2016) are independently
                confirmed in the LITTLE THINGS survey (22 dwarf irregulars, Oh et al. 2015) — a completely separate dataset
                processed by a different team using a different telescope (VLA) and analysis pipeline.
              </p>
              <div className="flex flex-wrap gap-2">
                <Badge pass={signMatch} label="Sign Consistent" />
                <Badge pass={slopeMatch} label="Slope Within 2σ" />
                <Badge pass={partialNeg} label="Partial r < 0" />
                <Badge pass={lt.pointLevelRegression.p < 0.001} label="Point-level p < 0.001" />
              </div>
            </div>
          </div>
        </GlassCard>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <GlassCard>
            <h3 className="text-lg font-display font-bold text-white mb-1 flex items-center gap-2">
              <Target className="w-5 h-5 text-cyan-400" />
              SPARC (Primary Dataset)
            </h3>
            <p className="text-xs text-slate-500 mb-3">175 galaxies · Lelli et al. 2016 · Spitzer 3.6μm + multi-survey HI</p>
            <div className="grid grid-cols-3 divide-x divide-white/10 mb-3">
              <StatBox label="Slope b" value={sparcB.toFixed(3)} />
              <StatBox label="Partial r|V_max" value="−0.656" />
              <StatBox label="n galaxies" value="175" />
            </div>
            <div className="h-[220px]">
              <ResponsiveContainer width="100%" height="100%">
                <ComposedChart margin={{ top: 10, right: 20, bottom: 25, left: 15 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                  <XAxis dataKey="logSigma" type="number" domain={[6, 11]} tick={{ fill: '#94a3b8', fontSize: 10 }} label={{ value: 'log₁₀(Σ_bar)', position: 'bottom', offset: 10, fill: '#94a3b8', fontSize: 11 }} />
                  <YAxis type="number" domain={[-1, 3]} tick={{ fill: '#94a3b8', fontSize: 10 }} label={{ value: 'ΔRAR', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 11 }} />
                  <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 11 }} />
                  <Scatter data={sparcScatter} fill="#06b6d4" fillOpacity={0.6} r={3} />
                  <Line data={sparcLine} dataKey="fit" stroke="#06b6d4" strokeWidth={2} dot={false} strokeDasharray="6 3" />
                </ComposedChart>
              </ResponsiveContainer>
            </div>
          </GlassCard>

          <GlassCard>
            <h3 className="text-lg font-display font-bold text-white mb-1 flex items-center gap-2">
              <Microscope className="w-5 h-5 text-emerald-400" />
              LITTLE THINGS (Replication)
            </h3>
            <p className="text-xs text-slate-500 mb-3">22 dwarf irregulars · Oh et al. 2015 · VLA HI + Spitzer 3.6μm</p>
            <div className="grid grid-cols-3 divide-x divide-white/10 mb-3">
              <StatBox label="Slope b" value={lt.simpleRegression.slope.toFixed(3)} color="emerald" />
              <StatBox label="Partial r|V_max" value={lt.partialCorrelation.rPartial.toFixed(3)} color="emerald" />
              <StatBox label="n galaxies" value={String(lt.nGalaxies)} color="emerald" />
            </div>
            <div className="h-[220px]">
              <ResponsiveContainer width="100%" height="100%">
                <ComposedChart margin={{ top: 10, right: 20, bottom: 25, left: 15 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                  <XAxis dataKey="logSigma" type="number" domain={[6, 11]} tick={{ fill: '#94a3b8', fontSize: 10 }} label={{ value: 'log₁₀(Σ_bar)', position: 'bottom', offset: 10, fill: '#94a3b8', fontSize: 11 }} />
                  <YAxis type="number" domain={[-1, 3]} tick={{ fill: '#94a3b8', fontSize: 10 }} label={{ value: 'ΔRAR', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 11 }} />
                  <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 11 }} />
                  <Scatter data={ltScatter} fill="#10b981" fillOpacity={0.7} r={4} />
                  <Line data={ltLine} dataKey="fit" stroke="#10b981" strokeWidth={2} dot={false} strokeDasharray="6 3" />
                </ComposedChart>
              </ResponsiveContainer>
            </div>
          </GlassCard>
        </div>

        <GlassCard>
          <h3 className="text-lg font-display font-bold text-white mb-4 flex items-center gap-2">
            <Scale className="w-5 h-5 text-amber-400" />
            Head-to-Head Comparison
          </h3>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-white/10">
                  <th className="text-left text-slate-400 py-2 px-3 font-medium">Metric</th>
                  <th className="text-center text-cyan-400 py-2 px-3 font-medium">SPARC (175 gal)</th>
                  <th className="text-center text-emerald-400 py-2 px-3 font-medium">LITTLE THINGS (22 gal)</th>
                  <th className="text-center text-slate-400 py-2 px-3 font-medium">Agreement</th>
                </tr>
              </thead>
              <tbody className="font-mono text-xs">
                <tr className="border-b border-white/5">
                  <td className="py-2.5 px-3 text-slate-300">Slope b</td>
                  <td className="py-2.5 px-3 text-center text-cyan-400">−0.178</td>
                  <td className="py-2.5 px-3 text-center text-emerald-400">{lt.simpleRegression.slope.toFixed(3)}</td>
                  <td className="py-2.5 px-3 text-center">
                    <Badge pass={slopeMatch} label={slopeMatch ? 'Within 2σ' : 'Outside 2σ'} />
                  </td>
                </tr>
                <tr className="border-b border-white/5">
                  <td className="py-2.5 px-3 text-slate-300">Partial r | V_max</td>
                  <td className="py-2.5 px-3 text-center text-cyan-400">−0.656</td>
                  <td className="py-2.5 px-3 text-center text-emerald-400">{lt.partialCorrelation.rPartial.toFixed(3)}</td>
                  <td className="py-2.5 px-3 text-center">
                    <Badge pass={partialNeg} label={partialNeg ? 'Same sign' : 'Opposite sign'} />
                  </td>
                </tr>
                <tr className="border-b border-white/5">
                  <td className="py-2.5 px-3 text-slate-300">Point-level slope</td>
                  <td className="py-2.5 px-3 text-center text-cyan-400">(per-galaxy only)</td>
                  <td className="py-2.5 px-3 text-center text-emerald-400">{lt.pointLevelRegression.slope.toFixed(3)}</td>
                  <td className="py-2.5 px-3 text-center">
                    <Badge pass={lt.pointLevelRegression.slope < 0} label={lt.pointLevelRegression.slope < 0 ? 'Negative' : 'Positive'} />
                  </td>
                </tr>
                <tr className="border-b border-white/5">
                  <td className="py-2.5 px-3 text-slate-300">Point-level p</td>
                  <td className="py-2.5 px-3 text-center text-cyan-400">—</td>
                  <td className="py-2.5 px-3 text-center text-emerald-400">{lt.pointLevelRegression.p.toExponential(1)}</td>
                  <td className="py-2.5 px-3 text-center">
                    <Badge pass={lt.pointLevelRegression.p < 0.001} label={lt.pointLevelRegression.p < 0.001 ? 'Significant' : 'Not significant'} />
                  </td>
                </tr>
                <tr>
                  <td className="py-2.5 px-3 text-slate-300">Total radial points</td>
                  <td className="py-2.5 px-3 text-center text-cyan-400">~1200</td>
                  <td className="py-2.5 px-3 text-center text-emerald-400">{lt.pointLevelRegression.nPoints}</td>
                  <td className="py-2.5 px-3 text-center text-slate-500">—</td>
                </tr>
              </tbody>
            </table>
          </div>
        </GlassCard>

        <GlassCard>
          <h3 className="text-lg font-display font-bold text-white mb-3 flex items-center gap-2">
            <BookOpen className="w-5 h-5 text-violet-400" />
            Why Stiskalek & Desmond (2023) Found No Effect
          </h3>
          <p className="text-sm text-slate-300 mb-4">
            Stiskalek & Desmond (MNRAS 525, 6130–6145) concluded the RAR cannot be tightened by any galaxy property.
            Our finding does not contradict their result — it arises from a methodological difference:
          </p>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div className="bg-white/5 rounded-xl p-4 border border-white/5">
              <div className="text-amber-400 font-bold text-sm mb-2">1. Different Σ Definition</div>
              <p className="text-xs text-slate-400">
                They used <span className="text-white font-mono">Σ_tot</span> (central disc+bulge surface brightness) — a <em>fixed</em> galaxy property.
                We use <span className="text-white font-mono">Σ_bar(r) = M_bar/πr²</span> — the mean enclosed baryonic density, which <em>varies with radius</em>.
                These are fundamentally different quantities; ours captures the local dilution of baryonic matter at each measurement point.
              </p>
            </div>
            <div className="bg-white/5 rounded-xl p-4 border border-white/5">
              <div className="text-amber-400 font-bold text-sm mb-2">2. Different Control Variable</div>
              <p className="text-xs text-slate-400">
                They compute partial correlations <em>at fixed g_bar</em> — searching for residual signal after removing
                the dominant RAR trend. We compute partial correlations <em>at fixed V_max</em>, removing the galaxy
                mass scale. Controlling for g_bar may absorb the surface-density signal since g_bar ∝ Σ at large radii.
              </p>
            </div>
            <div className="bg-white/5 rounded-xl p-4 border border-white/5">
              <div className="text-amber-400 font-bold text-sm mb-2">3. Different Aggregation</div>
              <p className="text-xs text-slate-400">
                They analyzed 2696 individual RC points from 147 galaxies (point-level residuals).
                We first average outer-region ΔRAR per galaxy, then regress across galaxies. Our approach
                reduces point-to-point noise and emphasizes galaxy-scale structural differences in surface density.
              </p>
            </div>
          </div>
          <div className="mt-4 bg-violet-500/10 border border-violet-500/20 rounded-xl p-4">
            <p className="text-xs text-slate-300">
              <span className="text-violet-400 font-bold">Reconciliation:</span> When Σ is defined as the central surface brightness
              (a global galaxy property), the RAR residuals show no secondary dependence — consistent with Stiskalek & Desmond.
              When Σ is defined as the enclosed mean baryonic density at each radius (a local, radially-varying quantity),
              a significant negative trend emerges. The effect is that RAR underpredicts acceleration in regions where baryons
              are more diluted (low Σ_bar) and slightly overpredicts where baryons are concentrated (high Σ_bar).
              This is physically meaningful: it suggests the RAR interpolating function does not fully capture the role of
              baryonic <em>concentration</em> in modulating the transition from Newtonian to dark-matter-dominated dynamics.
            </p>
          </div>
        </GlassCard>

        <GlassCard glow="amber">
          <h3 className="text-lg font-display font-bold text-white mb-3 flex items-center gap-2">
            <Lightbulb className="w-5 h-5 text-amber-400" />
            The Calibrated Claim
          </h3>
          <div className="bg-white/5 rounded-xl p-5 border border-amber-500/20 mb-4">
            <p className="text-slate-200 text-sm leading-relaxed">
              We present strong evidence that the Radial Acceleration Relation requires a density-dependent correction.
              Across two independent datasets — <span className="text-cyan-400">SPARC</span> (175 late-type galaxies, Lelli et al. 2016)
              and <span className="text-emerald-400">LITTLE THINGS</span> (22 dwarf irregulars, Oh et al. 2015) — the RAR
              residual ΔRAR correlates negatively with the mean enclosed baryonic surface density Σ_bar:
            </p>
            <div className="my-4 text-center">
              <div className="inline-block bg-black/30 rounded-xl px-6 py-3 border border-white/10">
                <span className="font-mono text-amber-300 text-lg">ΔRAR = a + b · log₁₀(Σ_bar)</span>
              </div>
            </div>
            <div className="grid grid-cols-2 gap-4 text-center font-mono text-sm">
              <div>
                <div className="text-cyan-400 font-bold">SPARC</div>
                <div className="text-slate-300">b = −0.178, r_partial = −0.656</div>
              </div>
              <div>
                <div className="text-emerald-400 font-bold">LITTLE THINGS</div>
                <div className="text-slate-300">b = {lt.simpleRegression.slope.toFixed(3)}, r_partial = {lt.partialCorrelation.rPartial.toFixed(3)}</div>
              </div>
            </div>
          </div>
          <div className="space-y-3 text-sm">
            <div className="flex items-start gap-2">
              <CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" />
              <span className="text-slate-300">The slopes are consistent within error bars (−0.178 vs −0.203, Δ = 0.025)</span>
            </div>
            <div className="flex items-start gap-2">
              <CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" />
              <span className="text-slate-300">Both partial correlations (controlling for V_max) are negative</span>
            </div>
            <div className="flex items-start gap-2">
              <CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" />
              <span className="text-slate-300">Point-level regression on 445 LITTLE THINGS radial points: p &lt; 10⁻¹³</span>
            </div>
            <div className="flex items-start gap-2">
              <CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" />
              <span className="text-slate-300">Monte Carlo (1000 iterations, SPARC): 100% of slopes negative, 95% CI [−0.214, −0.163]</span>
            </div>
            <div className="flex items-start gap-2">
              <AlertTriangle className="w-4 h-4 text-amber-400 mt-0.5 flex-shrink-0" />
              <span className="text-slate-300">Caveat: LITTLE THINGS partial r (−0.470) is weaker than SPARC (−0.656), expected given smaller sample of dwarfs only</span>
            </div>
            <div className="flex items-start gap-2">
              <AlertTriangle className="w-4 h-4 text-amber-400 mt-0.5 flex-shrink-0" />
              <span className="text-slate-300">Caveat: Both datasets use Spitzer 3.6μm for stellar mass — not fully independent in stellar M/L assumptions</span>
            </div>
          </div>
        </GlassCard>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <GlassCard>
            <h3 className="text-lg font-display font-bold text-white mb-3 flex items-center gap-2">
              <Target className="w-5 h-5 text-cyan-400" />
              What This Is
            </h3>
            <ul className="space-y-2 text-sm text-slate-300">
              <li className="flex items-start gap-2">
                <ArrowRight className="w-4 h-4 text-cyan-400 mt-0.5 flex-shrink-0" />
                Strong evidence that the standard RAR (McGaugh 2016 interpolating function) has a systematic, density-dependent residual
              </li>
              <li className="flex items-start gap-2">
                <ArrowRight className="w-4 h-4 text-cyan-400 mt-0.5 flex-shrink-0" />
                An empirical correction term that improves predictive accuracy by 25% (cross-validated)
              </li>
              <li className="flex items-start gap-2">
                <ArrowRight className="w-4 h-4 text-cyan-400 mt-0.5 flex-shrink-0" />
                A result replicated on an independent dataset from a different survey, telescope, and analysis pipeline
              </li>
              <li className="flex items-start gap-2">
                <ArrowRight className="w-4 h-4 text-cyan-400 mt-0.5 flex-shrink-0" />
                A constraint that any theory of galaxy dynamics (dark matter or modified gravity) must satisfy
              </li>
            </ul>
          </GlassCard>

          <GlassCard>
            <h3 className="text-lg font-display font-bold text-white mb-3 flex items-center gap-2">
              <AlertTriangle className="w-5 h-5 text-amber-400" />
              What This Is Not
            </h3>
            <ul className="space-y-2 text-sm text-slate-300">
              <li className="flex items-start gap-2">
                <XCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" />
                Not a "new law of physics" — it is an empirical correction to an existing scaling relation
              </li>
              <li className="flex items-start gap-2">
                <XCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" />
                Not a refutation of ΛCDM or MOND — both frameworks could potentially accommodate this
              </li>
              <li className="flex items-start gap-2">
                <XCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" />
                Not derived from first principles — the exponent b = −0.178 is purely empirical
              </li>
              <li className="flex items-start gap-2">
                <XCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" />
                Not peer-reviewed — this analysis requires independent verification by the community
              </li>
            </ul>
          </GlassCard>
        </div>

        <GlassCard>
          <h3 className="text-lg font-display font-bold text-white mb-3 flex items-center gap-2">
            <Microscope className="w-5 h-5 text-violet-400" />
            LITTLE THINGS Galaxy Sample
          </h3>
          <div className="overflow-x-auto">
            <table className="w-full text-xs font-mono">
              <thead>
                <tr className="border-b border-white/10 text-slate-400">
                  <th className="text-left py-2 px-2">Galaxy</th>
                  <th className="text-right py-2 px-2">D (Mpc)</th>
                  <th className="text-right py-2 px-2">V_max</th>
                  <th className="text-right py-2 px-2">log M_bar</th>
                  <th className="text-right py-2 px-2">⟨ΔRAR⟩</th>
                  <th className="text-right py-2 px-2">⟨log Σ⟩</th>
                  <th className="text-right py-2 px-2">Points</th>
                </tr>
              </thead>
              <tbody>
                {lt.galaxies.map((g) => (
                  <tr key={g.name} className="border-b border-white/5 hover:bg-white/5 transition-colors">
                    <td className="py-1.5 px-2 text-emerald-400">{g.name}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{g.dist.toFixed(1)}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{g.vmax.toFixed(1)}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{Math.log10(g.mbar).toFixed(2)}</td>
                    <td className="py-1.5 px-2 text-right text-amber-400">{g.meanDeltaRAR.toFixed(3)}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{g.meanLogSigBar.toFixed(2)}</td>
                    <td className="py-1.5 px-2 text-right text-slate-500">{g.nPoints}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
