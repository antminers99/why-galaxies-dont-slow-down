import React, { useState, useEffect } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
  ResponsiveContainer, ScatterChart, Scatter, Cell, ReferenceLine, Label,
  LineChart, Line
} from 'recharts';
import {
  Shield, CheckCircle2, XCircle, AlertTriangle,
  Shuffle, FlipVertical, Ruler, Target, Zap
} from 'lucide-react';

interface BreakResults {
  baseline: { nGalaxies: number; slope: number; r: number; partialR_Vmax: number };
  test1_null: { pValueSlope: number; pValueR: number; sigma: number; nullMeanSlope: number; nullStdSlope: number; pass: boolean };
  test2_scramble: { pValueSlope: number; pValuePartialR: number; sigma: number; scrambleMeanSlope: number; scrambleStdSlope: number; pass: boolean };
  test3_ml: { results: Array<{ Yd: number; Yb: number; slope: number; r: number; partialR: number }>; allNegative: boolean; slopeRange: number[]; partialRRange: number[]; pass: boolean };
  test4_distance: { meanSlope: number; stdSlope: number; ci95: number[]; fracNegative: number; fracPartialRneg: number; pass: boolean };
  a0_stability: { originalA0: number; perturbedMean: number; perturbedStd: number; ci95: number[]; cv_percent: number };
  overall: { passed: number; total: number; verdict: string };
}

export default function BreakTestPage() {
  const [data, setData] = useState<BreakResults | null>(null);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'break-test-results.json')
      .then(r => r.json())
      .then(d => setData(d))
      .catch(() => {});
  }, []);

  if (!data) {
    return (
      <Layout>
        <div className="text-center text-slate-400 mt-20">Loading Phase C results...</div>
      </Layout>
    );
  }

  const tests = [
    { name: 'Null Test', icon: Shuffle, pass: data.test1_null.pass, sigma: data.test1_null.sigma },
    { name: '\u03A3_bar Scramble', icon: FlipVertical, pass: data.test2_scramble.pass, sigma: data.test2_scramble.sigma },
    { name: 'M/L Stress', icon: Ruler, pass: data.test3_ml.pass, sigma: null },
    { name: 'Distance MC', icon: Target, pass: data.test4_distance.pass, sigma: null },
  ];

  const mlData = data.test3_ml.results.map(r => ({
    label: `\u03A5d=${r.Yd} \u03A5b=${r.Yb}`,
    slope: r.slope,
    partialR: r.partialR,
    Yd: r.Yd,
  }));

  return (
    <Layout>
      <div className="space-y-8 max-w-7xl mx-auto">
        <div className="text-center space-y-3">
          <div className="flex items-center justify-center gap-3">
            <Shield className="w-8 h-8 text-amber-400" />
            <h1 className="text-3xl font-bold text-white font-['Space_Grotesk']">
              Phase C: Breaking the Idea
            </h1>
          </div>
          <p className="text-slate-400 max-w-2xl mx-auto text-sm leading-relaxed">
            If you can{"'"}t break your own idea, it starts to deserve a physical explanation.
            If you can break it, you{"'"}ve saved yourself years.
          </p>
        </div>

        <GlassCard glow="amber">
          <div className="flex items-center gap-3 mb-4">
            <Zap className="w-5 h-5 text-amber-400" />
            <h2 className="text-lg font-semibold text-white">Overall Verdict</h2>
          </div>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
            {tests.map((t, i) => (
              <div key={i} className={`p-4 rounded-xl border ${t.pass
                ? 'border-emerald-500/30 bg-emerald-500/5'
                : 'border-red-500/30 bg-red-500/5'
              }`}>
                <div className="flex items-center gap-2 mb-2">
                  {t.pass
                    ? <CheckCircle2 className="w-5 h-5 text-emerald-400" />
                    : <XCircle className="w-5 h-5 text-red-400" />
                  }
                  <span className="text-sm font-medium text-white">{t.name}</span>
                </div>
                <p className="text-xs text-slate-400">
                  {t.pass ? 'Signal survives' : 'Signal vulnerable'}
                  {t.sigma !== null && ` (${Math.abs(t.sigma).toFixed(1)}\u03C3)`}
                </p>
              </div>
            ))}
          </div>
          <div className={`p-4 rounded-xl text-center ${
            data.overall.passed === 4
              ? 'bg-emerald-500/10 border border-emerald-500/30'
              : 'bg-amber-500/10 border border-amber-500/30'
          }`}>
            <p className="text-lg font-bold text-white">
              {data.overall.passed}/{data.overall.total} tests passed
            </p>
            <p className="text-sm text-slate-300 mt-1">{data.overall.verdict}</p>
          </div>
        </GlassCard>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <GlassCard glow="cyan">
            <div className="flex items-center gap-2 mb-1">
              <Shuffle className="w-5 h-5 text-cyan-400" />
              <h3 className="text-base font-semibold text-white">Test 1: Null Test</h3>
              {data.test1_null.pass
                ? <CheckCircle2 className="w-4 h-4 text-emerald-400 ml-auto" />
                : <XCircle className="w-4 h-4 text-red-400 ml-auto" />
              }
            </div>
            <p className="text-xs text-slate-400 mb-4">
              Shuffle {"\u0394"}RAR across galaxies 1000 times. If the signal is real, shuffling kills it.
            </p>
            <div className="grid grid-cols-2 gap-3 mb-4">
              <StatBox label="Null mean slope" value={data.test1_null.nullMeanSlope.toFixed(4)} />
              <StatBox label="Real slope" value={data.baseline.slope.toFixed(4)} highlight />
              <StatBox label="p-value" value={data.test1_null.pValueSlope < 0.001 ? '< 0.001' : data.test1_null.pValueSlope.toFixed(4)} />
              <StatBox label="Significance" value={`${Math.abs(data.test1_null.sigma).toFixed(1)}\u03C3`} highlight />
            </div>
            <NullDistChart
              nullMean={data.test1_null.nullMeanSlope}
              nullStd={data.test1_null.nullStdSlope}
              realValue={data.baseline.slope}
              label="Slope"
            />
          </GlassCard>

          <GlassCard glow="purple">
            <div className="flex items-center gap-2 mb-1">
              <FlipVertical className="w-5 h-5 text-purple-400" />
              <h3 className="text-base font-semibold text-white">Test 2: {"\u03A3"}_bar Scramble</h3>
              {data.test2_scramble.pass
                ? <CheckCircle2 className="w-4 h-4 text-emerald-400 ml-auto" />
                : <XCircle className="w-4 h-4 text-red-400 ml-auto" />
              }
            </div>
            <p className="text-xs text-slate-400 mb-4">
              Randomly reassign {"\u03A3"}_bar between galaxies. Tests whether coupling is specific.
            </p>
            <div className="grid grid-cols-2 gap-3 mb-4">
              <StatBox label="Scrambled mean slope" value={data.test2_scramble.scrambleMeanSlope.toFixed(4)} />
              <StatBox label="Real slope" value={data.baseline.slope.toFixed(4)} highlight />
              <StatBox label="p-value (partial r)" value={data.test2_scramble.pValuePartialR < 0.001 ? '< 0.001' : data.test2_scramble.pValuePartialR.toFixed(4)} />
              <StatBox label="Significance" value={`${Math.abs(data.test2_scramble.sigma).toFixed(1)}\u03C3`} highlight />
            </div>
            <NullDistChart
              nullMean={data.test2_scramble.scrambleMeanSlope}
              nullStd={data.test2_scramble.scrambleStdSlope}
              realValue={data.baseline.slope}
              label="Slope"
            />
          </GlassCard>
        </div>

        <GlassCard glow="cyan">
          <div className="flex items-center gap-2 mb-1">
            <Ruler className="w-5 h-5 text-cyan-400" />
            <h3 className="text-base font-semibold text-white">Test 3: Extreme M/L Stress Test</h3>
            {data.test3_ml.pass
              ? <CheckCircle2 className="w-4 h-4 text-emerald-400 ml-auto" />
              : <XCircle className="w-4 h-4 text-red-400 ml-auto" />
            }
          </div>
          <p className="text-xs text-slate-400 mb-4">
            21 combinations of {"\u03A5"}_disk (0.2{"\u2013"}0.8) and {"\u03A5"}_bulge (0.5{"\u2013"}1.0), each fully recomputed.
            All slopes must remain negative for the signal to survive.
          </p>
          <div className="grid grid-cols-3 gap-3 mb-6">
            <StatBox label="Slope range" value={`${data.test3_ml.slopeRange[0].toFixed(3)} to ${data.test3_ml.slopeRange[1].toFixed(3)}`} />
            <StatBox label="Partial r range" value={`${data.test3_ml.partialRRange[0].toFixed(3)} to ${data.test3_ml.partialRRange[1].toFixed(3)}`} />
            <StatBox label="All negative?" value={data.test3_ml.allNegative ? 'YES' : 'NO'} highlight />
          </div>
          <div className="h-64">
            <ResponsiveContainer width="100%" height="100%">
              <BarChart data={mlData} margin={{ top: 5, right: 20, bottom: 30, left: 20 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
                <XAxis dataKey="label" tick={{ fill: '#94a3b8', fontSize: 9 }} angle={-45} textAnchor="end" height={60} />
                <YAxis tick={{ fill: '#94a3b8', fontSize: 11 }}>
                  <Label value="Slope b" angle={-90} position="insideLeft" style={{ fill: '#94a3b8', fontSize: 12 }} />
                </YAxis>
                <Tooltip
                  contentStyle={{ backgroundColor: '#1e293b', border: '1px solid #334155', borderRadius: 8 }}
                  labelStyle={{ color: '#e2e8f0' }}
                  formatter={(v: number) => [v.toFixed(4), 'Slope']}
                />
                <ReferenceLine y={0} stroke="#ef4444" strokeDasharray="4 4" />
                <Bar dataKey="slope" radius={[4, 4, 0, 0]}>
                  {mlData.map((_, i) => (
                    <Cell key={i} fill={mlData[i].slope < 0 ? '#22d3ee' : '#ef4444'} fillOpacity={0.7} />
                  ))}
                </Bar>
              </BarChart>
            </ResponsiveContainer>
          </div>
        </GlassCard>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <GlassCard glow="purple">
            <div className="flex items-center gap-2 mb-1">
              <Target className="w-5 h-5 text-purple-400" />
              <h3 className="text-base font-semibold text-white">Test 4: Distance Monte Carlo</h3>
              {data.test4_distance.pass
                ? <CheckCircle2 className="w-4 h-4 text-emerald-400 ml-auto" />
                : <XCircle className="w-4 h-4 text-red-400 ml-auto" />
              }
            </div>
            <p className="text-xs text-slate-400 mb-4">
              Perturb every galaxy{"'"}s distance by {"\u00B1"}30% (Gaussian), 1000 iterations.
              Propagate through {"\u03A3"}_bar, V_max, and {"\u0394"}RAR.
            </p>
            <div className="grid grid-cols-2 gap-3 mb-4">
              <StatBox label="Mean perturbed slope" value={data.test4_distance.meanSlope.toFixed(4)} />
              <StatBox label={"\u00B1 1\u03C3"} value={data.test4_distance.stdSlope.toFixed(4)} />
              <StatBox label="95% CI" value={`[${data.test4_distance.ci95[0].toFixed(3)}, ${data.test4_distance.ci95[1].toFixed(3)}]`} />
              <StatBox label="Fraction negative" value={`${(data.test4_distance.fracNegative * 100).toFixed(1)}%`} highlight />
            </div>
            <div className="grid grid-cols-2 gap-3">
              <StatBox label="Partial r fraction neg" value={`${(data.test4_distance.fracPartialRneg * 100).toFixed(1)}%`} />
              <StatBox label="Real slope" value={data.baseline.slope.toFixed(4)} />
            </div>
          </GlassCard>

          <GlassCard glow="amber">
            <div className="flex items-center gap-2 mb-1">
              <AlertTriangle className="w-5 h-5 text-amber-400" />
              <h3 className="text-base font-semibold text-white">a{"\u2080"} Stability</h3>
            </div>
            <p className="text-xs text-slate-400 mb-4">
              How stable is the acceleration scale a{"\u2080"} when distances are perturbed by {"\u00B1"}30%?
            </p>
            <div className="grid grid-cols-2 gap-3 mb-4">
              <StatBox label="Original a\u2080" value={`${data.a0_stability.originalA0} (km/s)\u00B2/kpc`} />
              <StatBox label="Perturbed mean" value={`${data.a0_stability.perturbedMean.toFixed(1)} \u00B1 ${data.a0_stability.perturbedStd.toFixed(1)}`} />
              <StatBox label="95% CI" value={`[${data.a0_stability.ci95[0].toFixed(0)}, ${data.a0_stability.ci95[1].toFixed(0)}]`} />
              <StatBox label="CV" value={`${data.a0_stability.cv_percent.toFixed(1)}%`} highlight />
            </div>
            <p className="text-xs text-slate-500 mt-2">
              Note: Perturbed mean uses V{"\u00B2"}/R proxy. The low CV ({data.a0_stability.cv_percent.toFixed(1)}%)
              confirms a{"\u2080"} is well-constrained even under extreme distance noise.
            </p>
          </GlassCard>
        </div>

        <GlassCard glow="none">
          <h3 className="text-base font-semibold text-white mb-3">What This Means</h3>
          <div className="space-y-3 text-sm text-slate-300 leading-relaxed">
            <p>
              The signal {"\u2014"} the anti-correlation between dark matter fraction and baryonic surface density,
              and the universal acceleration scale a{"\u2080"} {"\u2014"} survives all four breaking attempts:
            </p>
            <ul className="space-y-2 ml-4">
              <li className="flex gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 shrink-0" />
                <span>It is NOT an artifact of random pairing (null test p {"<"} 0.001)</span>
              </li>
              <li className="flex gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 shrink-0" />
                <span>It is NOT an artifact of {"\u03A3"}_bar labeling (scramble p {"<"} 0.001)</span>
              </li>
              <li className="flex gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 shrink-0" />
                <span>It is NOT an artifact of mass-to-light ratio assumptions (21 {"\u03A5"} combinations, all negative)</span>
              </li>
              <li className="flex gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 shrink-0" />
                <span>It is NOT an artifact of distance errors ({"\u00B1"}30% {"\u2192"} 100% still negative)</span>
              </li>
            </ul>
            <p className="text-slate-400 text-xs mt-4 italic">
              This does NOT prove the Cosmic Floor hypothesis. It proves the PHENOMENON is real.
              The interpretation (modified gravity, feedback, or cosmological coupling) remains open.
            </p>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}

function StatBox({ label, value, highlight }: { label: string; value: string; highlight?: boolean }) {
  return (
    <div className="p-3 rounded-lg bg-slate-800/50 border border-slate-700/50">
      <p className="text-[10px] text-slate-500 uppercase tracking-wider mb-1">{label}</p>
      <p className={`text-sm font-mono font-semibold ${highlight ? 'text-cyan-400' : 'text-slate-200'}`}>
        {value}
      </p>
    </div>
  );
}

function NullDistChart({ nullMean, nullStd, realValue, label }: {
  nullMean: number; nullStd: number; realValue: number; label: string;
}) {
  const points = [];
  for (let i = -4; i <= 4; i += 0.2) {
    const x = nullMean + i * nullStd;
    const y = Math.exp(-0.5 * i * i) / (nullStd * Math.sqrt(2 * Math.PI));
    points.push({ x: parseFloat(x.toFixed(4)), y: parseFloat(y.toFixed(2)) });
  }

  return (
    <div className="h-40">
      <ResponsiveContainer width="100%" height="100%">
        <LineChart data={points} margin={{ top: 5, right: 20, bottom: 20, left: 20 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
          <XAxis dataKey="x" tick={{ fill: '#94a3b8', fontSize: 10 }} tickFormatter={(v: number) => v.toFixed(2)}>
            <Label value={label} position="bottom" offset={0} style={{ fill: '#94a3b8', fontSize: 11 }} />
          </XAxis>
          <YAxis tick={false} />
          <Line type="monotone" dataKey="y" stroke="#64748b" strokeWidth={2} dot={false} />
          <ReferenceLine x={parseFloat(realValue.toFixed(4))} stroke="#ef4444" strokeWidth={2} strokeDasharray="4 4">
            <Label value={`Real: ${realValue.toFixed(3)}`} position="top" style={{ fill: '#ef4444', fontSize: 10 }} />
          </ReferenceLine>
          <ReferenceLine x={parseFloat(nullMean.toFixed(4))} stroke="#22d3ee" strokeWidth={1} strokeDasharray="2 2" />
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
}
