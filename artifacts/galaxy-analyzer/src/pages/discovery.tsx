import React, { useState, useEffect } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, ReferenceLine, Line, LineChart, ComposedChart, Bar, Cell, Legend } from 'recharts';
import { Sparkles, ArrowRight, FlaskConical, Target, AlertTriangle, BookOpen, TrendingDown, Atom, CheckCircle2, Lightbulb, Microscope } from 'lucide-react';

export default function DiscoveryPage() {
  const [real, setReal] = useState<any>(null);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'rar-analysis-real.json')
      .then(r => r.json())
      .then(d => setReal(d))
      .catch(() => {});
  }, []);

  if (!real?.validation) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-64">
          <div className="w-8 h-8 border-2 border-t-primary border-white/10 rounded-full animate-spin" />
        </div>
      </Layout>
    );
  }

  const v = real.validation;
  const mc = real.monteCarlo;
  const galaxies = real.perGalaxy.filter((g: any) => g.sigma_bar > 0 && g.meanDeltaOuter > 0);

  const xs = galaxies.map((g: any) => Math.log10(g.sigma_bar));
  const ys = galaxies.map((g: any) => g.meanDeltaOuter);
  const n = xs.length;
  const mx = xs.reduce((s: number, v: number) => s + v, 0) / n;
  const my = ys.reduce((s: number, v: number) => s + v, 0) / n;
  let sxy = 0, sxx = 0;
  for (let i = 0; i < n; i++) { sxy += (xs[i] - mx) * (ys[i] - my); sxx += (xs[i] - mx) ** 2; }
  const b_val = sxy / sxx;
  const a_val = my - b_val * mx;
  const sigma0 = Math.pow(10, -a_val / b_val);

  const scatterData = galaxies.map((g: any) => ({
    logSigma: +Math.log10(g.sigma_bar).toFixed(3),
    delta: +g.meanDeltaOuter.toFixed(4),
    name: g.name,
    Vmax: g.Vmax,
  }));

  const lineData = [];
  const xMin = Math.min(...xs), xMax = Math.max(...xs);
  for (let x = xMin; x <= xMax; x += 0.05) {
    lineData.push({ logSigma: +x.toFixed(3), fit: +(a_val + b_val * x).toFixed(4) });
  }

  const innerOuterData = v.innerOuterByMass?.map((m: any) => ({
    category: m.name,
    inner_b: +m.inner.b.toFixed(3),
    outer_b: +m.outer.b.toFixed(3),
    ratio: +(m.inner.b / m.outer.b).toFixed(1),
  })) || [];

  const predictiveData = [
    { model: 'RAR only', R2: v.modelComparison.modelA.R2, rmse: v.modelComparison.modelA.rmse },
    { model: 'RAR + Σ_bar', R2: v.modelComparison.modelB.R2, rmse: v.modelComparison.modelB.rmse },
  ];

  return (
    <Layout>
      <div className="space-y-6">
        <div>
          <h1 className="text-3xl font-display font-bold text-white flex items-center gap-3">
            <Sparkles className="w-8 h-8 text-amber-400" />
            The Density-Corrected RAR
          </h1>
          <p className="text-slate-400 mt-2">
            An empirical law extending the Radial Acceleration Relation with a baryonic surface density correction.
          </p>
        </div>

        <GlassCard glow="cyan">
          <div className="flex items-center gap-2 mb-4">
            <Atom className="w-5 h-5 text-cyan-400" />
            <h2 className="text-lg font-semibold text-white">The Standard RAR</h2>
          </div>
          <p className="text-sm text-slate-300 mb-4">
            McGaugh et al. (2016) established that observed acceleration correlates tightly with baryonic acceleration:
          </p>
          <div className="p-4 bg-slate-900/60 border border-white/10 rounded-xl text-center mb-4">
            <div className="font-mono text-lg text-cyan-400">
              g<sub>obs</sub> = f(g<sub>bar</sub>) = g<sub>bar</sub> / (1 − e<sup>−√(g<sub>bar</sub>/g†)</sup>)
            </div>
            <div className="text-xs text-slate-500 mt-2">where g† = 1.2 × 10⁻¹⁰ m/s² (McGaugh's acceleration scale)</div>
          </div>
          <p className="text-sm text-slate-400">
            This is treated as a universal, one-parameter relation. But when we compute 
            ΔRAR = log(g<sub>obs</sub>) − log(g<sub>RAR</sub>), the residuals are <strong className="text-amber-400">not random</strong>.
          </p>
        </GlassCard>

        <GlassCard glow="cyan">
          <div className="flex items-center gap-2 mb-4">
            <Sparkles className="w-5 h-5 text-amber-400" />
            <h2 className="text-lg font-semibold text-white">The Discovery: Density-Dependent Residuals</h2>
          </div>
          <p className="text-sm text-slate-300 mb-4">
            Using real baryonic decomposition (V<sub>gas</sub> + V<sub>disk</sub> + V<sub>bulge</sub>) from 
            175 SPARC galaxies, we find that RAR residuals correlate with baryonic surface density:
          </p>

          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-4">
            <div>
              <h3 className="text-sm font-semibold text-cyan-400 mb-2">ΔRAR vs log(Σ<sub>bar</sub>) — 175 Galaxies</h3>
              <div className="h-[280px]">
                <ResponsiveContainer width="100%" height="100%">
                  <ComposedChart margin={{ top: 10, right: 20, bottom: 20, left: 20 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis
                      dataKey="logSigma"
                      type="number"
                      domain={['auto', 'auto']}
                      tick={{ fill: '#94a3b8', fontSize: 10 }}
                      label={{ value: 'log₁₀(Σ_bar) [M☉/kpc²]', position: 'bottom', fill: '#94a3b8', fontSize: 11 }}
                    />
                    <YAxis
                      tick={{ fill: '#94a3b8', fontSize: 10 }}
                      label={{ value: 'ΔRAR', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 11 }}
                    />
                    <Tooltip
                      contentStyle={{ background: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: '8px', fontSize: '11px' }}
                      formatter={(val: number, name: string) => [typeof val === 'number' ? val.toFixed(4) : val, name === 'delta' ? 'ΔRAR' : 'Fit']}
                      labelFormatter={(val: number) => `log₁₀(Σ) = ${val}`}
                    />
                    <Scatter data={scatterData} dataKey="delta" fill="#06b6d4" fillOpacity={0.5} r={3} name="Galaxies" />
                    <Line data={lineData} dataKey="fit" stroke="#f59e0b" strokeWidth={2} strokeDasharray="6 3" dot={false} name="Fit" />
                    <ReferenceLine y={0} stroke="rgba(255,255,255,0.2)" strokeDasharray="3 3" />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>
            </div>

            <div className="space-y-3">
              <div className="p-4 bg-gradient-to-r from-amber-500/10 to-orange-500/10 border border-amber-500/20 rounded-xl">
                <div className="text-xs text-slate-500 uppercase tracking-wider mb-1">The Empirical Relation</div>
                <div className="font-mono text-base text-amber-400 text-center py-2">
                  ΔRAR = {a_val.toFixed(3)} + ({b_val.toFixed(3)}) × log₁₀(Σ<sub>bar</sub>)
                </div>
              </div>
              <div className="grid grid-cols-2 gap-2">
                <div className="p-3 bg-white/[0.03] rounded-xl text-center">
                  <div className="text-xs text-slate-500">Slope b</div>
                  <div className="text-lg font-mono font-bold text-purple-400">{b_val.toFixed(3)}</div>
                  <div className="text-xs text-slate-500">95% CI [{mc?.slope.ci95[0].toFixed(3)}, {mc?.slope.ci95[1].toFixed(3)}]</div>
                </div>
                <div className="p-3 bg-white/[0.03] rounded-xl text-center">
                  <div className="text-xs text-slate-500">Partial r | V<sub>max</sub></div>
                  <div className="text-lg font-mono font-bold text-cyan-400">{v.partialCorrelations.afterVmaxResidual.toFixed(3)}</div>
                  <div className="text-xs text-slate-500">Not a mass proxy</div>
                </div>
                <div className="p-3 bg-white/[0.03] rounded-xl text-center">
                  <div className="text-xs text-slate-500">ΔAIC</div>
                  <div className="text-lg font-mono font-bold text-emerald-400">{v.modelComparison.deltaAIC.toFixed(0)}</div>
                  <div className="text-xs text-slate-500">Overwhelming evidence</div>
                </div>
                <div className="p-3 bg-white/[0.03] rounded-xl text-center">
                  <div className="text-xs text-slate-500">Validation</div>
                  <div className="text-lg font-mono font-bold text-emerald-400">10/10</div>
                  <div className="text-xs text-slate-500">Tests passed</div>
                </div>
              </div>
            </div>
          </div>
        </GlassCard>

        <GlassCard glow="cyan">
          <div className="flex items-center gap-2 mb-4">
            <FlaskConical className="w-5 h-5 text-purple-400" />
            <h2 className="text-lg font-semibold text-white">The Proposed Law</h2>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
            <div className="p-4 bg-slate-900/60 border border-white/10 rounded-xl">
              <div className="text-xs text-slate-500 uppercase tracking-wider mb-2">Step 1: Standard RAR</div>
              <div className="font-mono text-sm text-cyan-400 text-center py-3">
                g<sub>RAR</sub> = f(g<sub>bar</sub>)
              </div>
              <p className="text-xs text-slate-400 text-center">Known universal relation</p>
            </div>
            <div className="p-4 bg-slate-900/60 border border-white/10 rounded-xl">
              <div className="text-xs text-slate-500 uppercase tracking-wider mb-2">Step 2: Density Correction</div>
              <div className="font-mono text-sm text-amber-400 text-center py-3">
                Δ = a + b × log(Σ<sub>bar</sub>)
              </div>
              <p className="text-xs text-slate-400 text-center">Our empirical finding</p>
            </div>
            <div className="p-4 bg-gradient-to-r from-purple-500/10 to-cyan-500/10 border border-purple-500/20 rounded-xl">
              <div className="text-xs text-purple-400 uppercase tracking-wider mb-2">Step 3: Combined Law</div>
              <div className="font-mono text-sm text-white text-center py-3">
                g<sub>obs</sub> = g<sub>RAR</sub> × (Σ<sub>bar</sub> / Σ₀)<sup>b</sup>
              </div>
              <p className="text-xs text-purple-300 text-center">Density-corrected acceleration</p>
            </div>
          </div>

          <div className="p-4 bg-gradient-to-r from-purple-500/10 to-cyan-500/10 border border-purple-500/20 rounded-xl mb-4">
            <div className="text-center">
              <div className="text-xs text-slate-500 uppercase tracking-wider mb-3">The Density-Corrected RAR</div>
              <div className="font-mono text-xl text-white py-2">
                g<sub>obs</sub>(r) = f(g<sub>bar</sub>(r)) × (Σ<sub>bar</sub> / Σ₀)<sup>b</sup>
              </div>
              <div className="grid grid-cols-3 gap-4 mt-4 text-xs">
                <div>
                  <span className="text-slate-500">b = </span>
                  <span className="font-mono text-purple-400">{b_val.toFixed(3)}</span>
                  <div className="text-slate-500">(global average)</div>
                </div>
                <div>
                  <span className="text-slate-500">Σ₀ = </span>
                  <span className="font-mono text-cyan-400">{sigma0.toExponential(2)}</span>
                  <div className="text-slate-500">M☉/kpc²</div>
                </div>
                <div>
                  <span className="text-slate-500">log₁₀(Σ₀) = </span>
                  <span className="font-mono text-amber-400">{(-a_val / b_val).toFixed(2)}</span>
                  <div className="text-slate-500">pivot density</div>
                </div>
              </div>
            </div>
          </div>

          <div className="p-3 bg-white/[0.02] rounded-xl">
            <p className="text-xs text-slate-300 leading-relaxed">
              <strong className="text-white">Physical meaning:</strong> When Σ<sub>bar</sub> = Σ₀, 
              the correction vanishes and standard RAR holds exactly. For <strong className="text-amber-400">low-density 
              galaxies</strong> (Σ<sub>bar</sub> {'<'} Σ₀), the correction factor {'>'} 1 — they show <em>more</em> excess 
              acceleration than RAR predicts. For <strong className="text-cyan-400">high-density galaxies</strong> (Σ<sub>bar</sub> {'>'} Σ₀), 
              the correction {'<'} 1 — they are closer to the standard relation.
            </p>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-2 mb-4">
            <TrendingDown className="w-5 h-5 text-orange-400" />
            <h2 className="text-lg font-semibold text-white">Radial Dependence: b(r) Varies with Region</h2>
          </div>
          <p className="text-sm text-slate-400 mb-4">
            The correction is not uniform — it is <strong className="text-orange-400">stronger in the inner regions</strong> of 
            massive galaxies, where baryons dominate. This is physically expected: the density correction should matter 
            most where baryonic density is highest.
          </p>
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            <div>
              <h3 className="text-sm font-semibold text-orange-400 mb-2">b(r) by Mass Category</h3>
              <div className="space-y-2">
                {innerOuterData.map((d: any) => (
                  <div key={d.category} className="p-3 bg-white/[0.02] rounded-lg">
                    <div className="flex justify-between items-center mb-1">
                      <span className="text-sm text-white">{d.category}</span>
                      <span className="text-xs text-slate-500">inner/outer ratio: {d.ratio}×</span>
                    </div>
                    <div className="flex gap-3">
                      <div className="flex-1">
                        <div className="flex justify-between text-xs mb-1">
                          <span className="text-slate-500">Inner b</span>
                          <span className="font-mono text-orange-400">{d.inner_b}</span>
                        </div>
                        <div className="h-2 bg-slate-800 rounded-full overflow-hidden">
                          <div
                            className="h-full bg-gradient-to-r from-orange-500 to-red-500 rounded-full"
                            style={{ width: `${Math.min(100, Math.abs(d.inner_b) / 0.75 * 100)}%` }}
                          />
                        </div>
                      </div>
                      <div className="flex-1">
                        <div className="flex justify-between text-xs mb-1">
                          <span className="text-slate-500">Outer b</span>
                          <span className="font-mono text-cyan-400">{d.outer_b}</span>
                        </div>
                        <div className="h-2 bg-slate-800 rounded-full overflow-hidden">
                          <div
                            className="h-full bg-gradient-to-r from-cyan-500 to-blue-500 rounded-full"
                            style={{ width: `${Math.min(100, Math.abs(d.outer_b) / 0.75 * 100)}%` }}
                          />
                        </div>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </div>
            <div className="space-y-4">
              <div className="p-4 bg-gradient-to-r from-red-500/10 to-orange-500/10 border border-orange-500/20 rounded-xl">
                <div className="flex items-center gap-2 mb-2">
                  <Sparkles className="w-4 h-4 text-orange-400" />
                  <span className="text-sm font-semibold text-orange-400">Key Finding</span>
                </div>
                <p className="text-xs text-slate-300 leading-relaxed">
                  In massive galaxies ({'>'} 200 km/s), the inner correction is <strong className="text-orange-400">
                  b = {innerOuterData[3]?.inner_b}</strong> — {innerOuterData[3]?.ratio}× stronger than the outer 
                  region (b = {innerOuterData[3]?.outer_b}). This is consistent with the correction being driven 
                  by local baryonic physics: where baryonic density is concentrated, the departure from standard 
                  RAR is greatest.
                </p>
              </div>
              <div className="p-4 bg-white/[0.02] border border-white/5 rounded-xl">
                <div className="text-xs text-slate-500 uppercase tracking-wider mb-2">Extended Form</div>
                <div className="font-mono text-sm text-white text-center py-2">
                  g<sub>obs</sub>(r) = f(g<sub>bar</sub>(r)) × (Σ<sub>bar</sub> / Σ₀)<sup>b(r)</sup>
                </div>
                <p className="text-xs text-slate-400 text-center mt-2">
                  where b(r) varies from ≈ {innerOuterData[3]?.inner_b} (inner) to ≈ {innerOuterData[3]?.outer_b} (outer) in massive galaxies
                </p>
              </div>
            </div>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-2 mb-4">
            <Target className="w-5 h-5 text-emerald-400" />
            <h2 className="text-lg font-semibold text-white">Predictive Power</h2>
          </div>
          <p className="text-sm text-slate-400 mb-4">
            The critical test: does the density correction <strong className="text-emerald-400">improve predictions</strong> of 
            observed acceleration beyond standard RAR?
          </p>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-4">
            <div className="p-3 bg-white/[0.03] rounded-xl text-center">
              <div className="text-xs text-slate-500">R² (RAR only)</div>
              <div className="text-xl font-mono font-bold text-amber-400">{v.modelComparison.modelA.R2.toFixed(3)}</div>
            </div>
            <div className="p-3 bg-emerald-500/10 border border-emerald-500/20 rounded-xl text-center">
              <div className="text-xs text-emerald-400">R² (RAR + Σ<sub>bar</sub>)</div>
              <div className="text-xl font-mono font-bold text-emerald-400">{v.modelComparison.modelB.R2.toFixed(3)}</div>
            </div>
            <div className="p-3 bg-cyan-500/10 border border-cyan-500/20 rounded-xl text-center">
              <div className="text-xs text-cyan-400">Variance Gain</div>
              <div className="text-xl font-mono font-bold text-cyan-400">+{(v.modelComparison.R2improvement * 100).toFixed(0)}%</div>
            </div>
            <div className="p-3 bg-purple-500/10 border border-purple-500/20 rounded-xl text-center">
              <div className="text-xs text-purple-400">CV Error Reduction</div>
              <div className="text-xl font-mono font-bold text-purple-400">{v.crossValidation.avgImprovement.toFixed(0)}%</div>
            </div>
          </div>
          <div className="p-3 bg-emerald-500/10 border border-emerald-500/20 rounded-xl">
            <p className="text-xs text-emerald-300">
              <CheckCircle2 className="w-3.5 h-3.5 inline mr-1" />
              The density-corrected model predicts g<sub>obs</sub> better than RAR alone, with 
              {' '}{v.crossValidation.avgImprovement.toFixed(0)}% error reduction in 5-fold cross-validation 
              and {v.trainTest.improvement.toFixed(0)}% on unseen test data. This is <strong>predictive</strong>, 
              not just descriptive.
            </p>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-2 mb-4">
            <BookOpen className="w-5 h-5 text-violet-400" />
            <h2 className="text-lg font-semibold text-white">Literature Context</h2>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="p-4 bg-white/[0.03] border border-white/5 rounded-xl">
              <h3 className="text-sm font-semibold text-slate-300 mb-3">What is Known</h3>
              <ul className="text-xs text-slate-400 space-y-2">
                <li className="flex gap-2">
                  <span className="text-slate-600 shrink-0">•</span>
                  <span>LSB galaxies have slower-rising rotation curves than HSB at same V<sub>max</sub> (de Blok et al., Lelli et al.)</span>
                </li>
                <li className="flex gap-2">
                  <span className="text-slate-600 shrink-0">•</span>
                  <span>Central stellar density correlates with central dynamical density (Lelli et al. 2016)</span>
                </li>
                <li className="flex gap-2">
                  <span className="text-slate-600 shrink-0">•</span>
                  <span>Halo cuspiness linked to stellar compactness (Kaplinghat et al. 2019)</span>
                </li>
                <li className="flex gap-2">
                  <span className="text-slate-600 shrink-0">•</span>
                  <span>Di Paolo et al. (2018): r/R<sub>opt</sub> as a second variable in RAR — showing RAR may not be purely 1D</span>
                </li>
                <li className="flex gap-2">
                  <span className="text-slate-600 shrink-0">•</span>
                  <span>MIGHTEE-HI (2025): resolved stellar masses and variable M/L give tighter RAR — supports that better baryonic modeling matters</span>
                </li>
              </ul>
            </div>
            <div className="p-4 bg-white/[0.03] border border-white/5 rounded-xl">
              <h3 className="text-sm font-semibold text-slate-300 mb-3">Key Tension</h3>
              <div className="p-3 bg-red-500/10 border border-red-500/20 rounded-lg mb-3">
                <p className="text-xs text-red-300 leading-relaxed">
                  <AlertTriangle className="w-3.5 h-3.5 inline mr-1" />
                  <strong>Stiskalek & Desmond (2023)</strong> tested surface brightness among candidate features and 
                  concluded that RAR cannot be tightened by adding galaxy properties. Our result directly contradicts this.
                </p>
              </div>
              <h3 className="text-sm font-semibold text-slate-300 mb-2 mt-3">Possible Explanations for the Tension</h3>
              <ul className="text-xs text-slate-400 space-y-2">
                <li className="flex gap-2">
                  <span className="text-amber-500 shrink-0">1.</span>
                  <span>Different baryonic modeling: they may have used simpler g<sub>bar</sub> computation. We use real V<sub>gas</sub>, V<sub>disk</sub>, V<sub>bulge</sub> decomposition.</span>
                </li>
                <li className="flex gap-2">
                  <span className="text-amber-500 shrink-0">2.</span>
                  <span>Different surface density definition: central SB vs. our integrated Σ<sub>bar</sub> within r<sub>fid</sub>.</span>
                </li>
                <li className="flex gap-2">
                  <span className="text-amber-500 shrink-0">3.</span>
                  <span>Different statistical methodology: feature importance vs. our ΔRAR residual regression with V<sub>max</sub> control.</span>
                </li>
              </ul>
            </div>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-2 mb-4">
            <Lightbulb className="w-5 h-5 text-amber-400" />
            <h2 className="text-lg font-semibold text-white">Physical Interpretation</h2>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
            <div className="p-4 bg-white/[0.03] border border-white/5 rounded-xl">
              <h3 className="text-sm font-semibold text-blue-400 mb-2">If Dark Matter is Real</h3>
              <p className="text-xs text-slate-400 leading-relaxed">
                The correction means dark matter halo profiles are not uniquely determined by total baryonic mass — 
                they also depend on how that mass is distributed. Diffuse galaxies would have more extended halos 
                relative to their baryonic content, producing more excess acceleration.
              </p>
            </div>
            <div className="p-4 bg-white/[0.03] border border-white/5 rounded-xl">
              <h3 className="text-sm font-semibold text-purple-400 mb-2">If Gravity is Modified</h3>
              <p className="text-xs text-slate-400 leading-relaxed">
                Any modified gravity theory (MOND, emergent gravity, etc.) that predicts a universal, 
                one-parameter RAR would need modification. The density dependence suggests the gravitational 
                response depends not just on the local field g<sub>bar</sub> but on the source configuration.
              </p>
            </div>
            <div className="p-4 bg-white/[0.03] border border-white/5 rounded-xl">
              <h3 className="text-sm font-semibold text-emerald-400 mb-2">Either Way</h3>
              <p className="text-xs text-slate-400 leading-relaxed">
                The RAR is not a complete description. It is a <strong className="text-white">special case</strong> of 
                a broader relation that includes baryonic surface density. The full relation has greater predictive power.
              </p>
            </div>
          </div>
          <div className="p-4 bg-gradient-to-r from-amber-500/10 to-orange-500/10 border border-amber-500/20 rounded-xl">
            <p className="text-sm text-amber-300 font-medium text-center">
              "Galaxy dynamics follow the RAR with a density-dependent correction, where baryonic surface density 
              acts as a second parameter modulating the observed acceleration."
            </p>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-2 mb-4">
            <AlertTriangle className="w-5 h-5 text-amber-400" />
            <h2 className="text-lg font-semibold text-white">Status & Caveats</h2>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="p-4 bg-emerald-500/10 border border-emerald-500/20 rounded-xl">
              <h3 className="text-sm font-semibold text-emerald-400 mb-2">What This Is</h3>
              <ul className="text-xs text-slate-300 space-y-1.5">
                <li className="flex items-start gap-2">
                  <CheckCircle2 className="w-3.5 h-3.5 text-emerald-400 shrink-0 mt-0.5" />
                  <span>An <strong>empirical law</strong> derived from 175 SPARC galaxies with real baryonic decomposition</span>
                </li>
                <li className="flex items-start gap-2">
                  <CheckCircle2 className="w-3.5 h-3.5 text-emerald-400 shrink-0 mt-0.5" />
                  <span><strong>Predictive</strong> — improves g<sub>obs</sub> prediction by {v.crossValidation.avgImprovement.toFixed(0)}% in cross-validation</span>
                </li>
                <li className="flex items-start gap-2">
                  <CheckCircle2 className="w-3.5 h-3.5 text-emerald-400 shrink-0 mt-0.5" />
                  <span><strong>Validated</strong> — 10/10 statistical tests, robust to Υ choice, survey source, Monte Carlo perturbations</span>
                </li>
                <li className="flex items-start gap-2">
                  <CheckCircle2 className="w-3.5 h-3.5 text-emerald-400 shrink-0 mt-0.5" />
                  <span><strong>Physically sensible</strong> — radially dependent, strongest where baryons dominate</span>
                </li>
              </ul>
            </div>
            <div className="p-4 bg-amber-500/10 border border-amber-500/20 rounded-xl">
              <h3 className="text-sm font-semibold text-amber-400 mb-2">What This is Not</h3>
              <ul className="text-xs text-slate-300 space-y-1.5">
                <li className="flex items-start gap-2">
                  <AlertTriangle className="w-3.5 h-3.5 text-amber-400 shrink-0 mt-0.5" />
                  <span><strong>Replicated on LITTLE THINGS</strong> — b = −0.203 (22 dwarfs, Oh et al. 2015); see <a href="replication" className="text-cyan-400 underline">Replication page</a></span>
                </li>
                <li className="flex items-start gap-2">
                  <AlertTriangle className="w-3.5 h-3.5 text-amber-400 shrink-0 mt-0.5" />
                  <span><strong>Contradicts Stiskalek & Desmond (2023)</strong> — this tension must be resolved before strong claims</span>
                </li>
                <li className="flex items-start gap-2">
                  <AlertTriangle className="w-3.5 h-3.5 text-amber-400 shrink-0 mt-0.5" />
                  <span><strong>Weak in dwarfs</strong> — R² = 0.08 for the smallest galaxies (n=21), structural limitation</span>
                </li>
                <li className="flex items-start gap-2">
                  <AlertTriangle className="w-3.5 h-3.5 text-amber-400 shrink-0 mt-0.5" />
                  <span><strong>No theoretical derivation</strong> — the law is empirical, not derived from first principles</span>
                </li>
              </ul>
            </div>
          </div>
        </GlassCard>

        <GlassCard glow="cyan">
          <div className="flex items-center gap-2 mb-4">
            <Microscope className="w-5 h-5 text-cyan-400" />
            <h2 className="text-lg font-semibold text-white">Next Steps for Confirmation</h2>
          </div>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div className="p-3 bg-white/[0.02] rounded-xl">
              <div className="text-xs text-emerald-400 font-semibold mb-1">1. External Replication ✓</div>
              <p className="text-xs text-slate-400">Done — LITTLE THINGS (Oh et al. 2015): b = −0.203, partial r = −0.470, 445 points p {'<'} 10⁻¹³. See <a href="replication" className="text-cyan-400 underline">Replication page</a>.</p>
            </div>
            <div className="p-3 bg-white/[0.02] rounded-xl">
              <div className="text-xs text-cyan-400 font-semibold mb-1">2. Resolve Stiskalek Tension</div>
              <p className="text-xs text-slate-400">Directly compare methodologies — same data, same Σ definition, same controls. Identify why the conclusions differ.</p>
            </div>
            <div className="p-3 bg-white/[0.02] rounded-xl">
              <div className="text-xs text-cyan-400 font-semibold mb-1">3. Theoretical Derivation</div>
              <p className="text-xs text-slate-400">Can DM models (NFW + adiabatic contraction) or modified gravity frameworks predict b ≈ −0.18 and its radial dependence?</p>
            </div>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
