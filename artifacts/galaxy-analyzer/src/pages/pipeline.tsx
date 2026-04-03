import React, { useState, useEffect } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
  ResponsiveContainer, Cell, ReferenceLine, Label,
  LineChart, Line, Legend
} from 'recharts';
import {
  FlaskConical, Trophy, AlertTriangle, CheckCircle2,
  XCircle, ArrowRight, Target, Layers, GitCompare
} from 'lucide-react';

interface PipelineData {
  parameters: { true_shift_pct: number; sigmaV: number; sigmaM_dex: number; sigmaInc_deg: number; N_MC: number };
  comparison: Array<{ name: string; N20: number; N50: number; N100: number }>;
  methods: {
    naive: Array<{ N: number; meanShift: number; stdShift: number; detection_rate: number }>;
    btfr_matched: Array<{ N: number; meanShift: number; stdShift: number; detection_rate: number }>;
    ratio: Array<{ N: number; meanShift: number; stdShift: number; detection_rate: number }>;
    global_fit: Array<{ N: number; meanA0: number; stdA0: number; trueA0: number; bias: number; detection_rate: number }>;
    stacked_diff: Array<{ N: number; meanSlope: number; stdSlope: number; sigma: number; detection_rate: number }>;
    optimal: Array<{ N: number; meanRatio: number; stdRatio: number; detection_rate: number }>;
  };
}

const METHOD_COLORS: Record<string, string> = {
  'Naive V_flat': '#ef4444',
  'Mass-matched BTFR': '#22d3ee',
  'Ratio V/V_bar': '#a78bfa',
  'Global a\u2080 fit': '#f59e0b',
  'Stacked differential': '#34d399',
  'Optimal combined': '#f472b6',
};

export default function PipelinePage() {
  const [data, setData] = useState<PipelineData | null>(null);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'measurement-pipeline.json')
      .then(r => r.json())
      .then(d => setData(d))
      .catch(() => {});
  }, []);

  if (!data) {
    return (
      <Layout>
        <div className="text-center text-slate-400 mt-20">Loading pipeline data...</div>
      </Layout>
    );
  }

  const comparisonN50 = data.comparison.map(c => ({
    name: c.name,
    rate: Math.round(c.N50 * 100),
    color: METHOD_COLORS[c.name] || '#64748b',
  }));

  const evolutionData = [10, 20, 30, 50, 80, 100].map(N => {
    const row: any = { N };
    for (const c of data.comparison) {
      const methods: Record<string, any[]> = {
        'Naive V_flat': data.methods.naive,
        'Mass-matched BTFR': data.methods.btfr_matched,
        'Ratio V/V_bar': data.methods.ratio,
        'Global a\u2080 fit': data.methods.global_fit,
        'Stacked differential': data.methods.stacked_diff,
        'Optimal combined': data.methods.optimal,
      };
      const arr = methods[c.name];
      if (arr) {
        const entry = arr.find((r: any) => r.N === N);
        if (entry) row[c.name] = Math.round(entry.detection_rate * 100);
      }
    }
    return row;
  });

  const pipelineSteps = [
    { step: 1, title: 'Sample Selection', detail: '30-50 galaxies at z=0.5-1.5, M=10\u2079-10\u00B9\u00B9 M\u2609, inc 30-75\u00B0', icon: Target },
    { step: 2, title: 'Measurement', detail: 'V_flat from tilted-ring fit + M_bar from SED + gas', icon: FlaskConical },
    { step: 3, title: 'Mass-Matched Comparison', detail: 'Bin by mass, compare V within bins against SPARC z=0', icon: Layers },
    { step: 4, title: 'Fit Evolution', detail: 'V_flat(z)/V_flat(0) = E(z)^(\u03B1/4). Floor: \u03B1=1, MOND: \u03B1=0', icon: GitCompare },
    { step: 5, title: 'Significance', detail: 'Likelihood ratio test, bootstrap CI, require \u22653\u03C3', icon: CheckCircle2 },
  ];

  return (
    <Layout>
      <div className="space-y-8 max-w-7xl mx-auto">
        <div className="text-center space-y-3">
          <div className="flex items-center justify-center gap-3">
            <FlaskConical className="w-8 h-8 text-cyan-400" />
            <h1 className="text-3xl font-bold text-white font-['Space_Grotesk']">
              Measurement Pipeline
            </h1>
          </div>
          <p className="text-slate-400 max-w-2xl mx-auto text-sm leading-relaxed">
            We don{"'"}t make the signal bigger {"\u2014"} we make the noise smaller until the signal shows through.
          </p>
          <div className="flex items-center justify-center gap-4 text-xs text-slate-500">
            <span>Signal: +{data.parameters.true_shift_pct.toFixed(1)}% V shift</span>
            <span>{"\u2022"}</span>
            <span>{"\u03C3"}(V)={data.parameters.sigmaV * 100}%</span>
            <span>{"\u2022"}</span>
            <span>{"\u03C3"}(logM)={data.parameters.sigmaM_dex} dex</span>
            <span>{"\u2022"}</span>
            <span>{data.parameters.N_MC} MC trials</span>
          </div>
        </div>

        <GlassCard glow="cyan">
          <div className="flex items-center gap-3 mb-4">
            <Trophy className="w-6 h-6 text-amber-400" />
            <h2 className="text-lg font-semibold text-white">Key Finding</h2>
          </div>
          <div className="bg-cyan-500/10 border border-cyan-500/30 rounded-xl p-5 text-center">
            <p className="text-2xl font-bold text-cyan-400 mb-2">Mass-matched BTFR wins overwhelmingly</p>
            <p className="text-sm text-slate-300">
              93% detection rate with just 20 galaxies {"\u2014"} vs 20% for naive comparison.
              Mass matching removes the dominant systematic error.
            </p>
            <div className="mt-3 grid grid-cols-3 gap-3 max-w-md mx-auto">
              <div className="text-center">
                <div className="text-xl font-bold font-mono text-amber-400">93%</div>
                <div className="text-[10px] text-slate-500 uppercase">N=20</div>
              </div>
              <div className="text-center">
                <div className="text-xl font-bold font-mono text-emerald-400">99%</div>
                <div className="text-[10px] text-slate-500 uppercase">N=50</div>
              </div>
              <div className="text-center">
                <div className="text-xl font-bold font-mono text-cyan-400">100%</div>
                <div className="text-[10px] text-slate-500 uppercase">N=100</div>
              </div>
            </div>
          </div>
        </GlassCard>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <GlassCard glow="purple">
            <h3 className="text-base font-semibold text-white mb-4">Detection Rate at N=50 (6 Methods)</h3>
            <div className="h-72">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={comparisonN50} layout="vertical" margin={{ top: 5, right: 30, bottom: 5, left: 10 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
                  <XAxis type="number" domain={[0, 100]} tick={{ fill: '#94a3b8', fontSize: 11 }} tickFormatter={v => `${v}%`} />
                  <YAxis type="category" dataKey="name" tick={{ fill: '#94a3b8', fontSize: 10 }} width={130} />
                  <Tooltip
                    contentStyle={{ backgroundColor: '#1e293b', border: '1px solid #334155', borderRadius: 8 }}
                    formatter={(v: number) => [`${v}%`, 'Detection rate']}
                  />
                  <ReferenceLine x={50} stroke="#f59e0b" strokeDasharray="4 4">
                    <Label value="50%" position="top" style={{ fill: '#f59e0b', fontSize: 10 }} />
                  </ReferenceLine>
                  <Bar dataKey="rate" radius={[0, 4, 4, 0]}>
                    {comparisonN50.map((entry, i) => (
                      <Cell key={i} fill={entry.color} fillOpacity={0.8} />
                    ))}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
            </div>
          </GlassCard>

          <GlassCard glow="cyan">
            <h3 className="text-base font-semibold text-white mb-4">Detection Rate vs Sample Size</h3>
            <div className="h-72">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={evolutionData} margin={{ top: 5, right: 20, bottom: 20, left: 10 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
                  <XAxis dataKey="N" tick={{ fill: '#94a3b8', fontSize: 11 }}>
                    <Label value="N galaxies" position="bottom" offset={0} style={{ fill: '#94a3b8', fontSize: 11 }} />
                  </XAxis>
                  <YAxis tick={{ fill: '#94a3b8', fontSize: 11 }} domain={[0, 100]} tickFormatter={v => `${v}%`} />
                  <Tooltip
                    contentStyle={{ backgroundColor: '#1e293b', border: '1px solid #334155', borderRadius: 8 }}
                    formatter={(v: number, name: string) => [`${v}%`, name]}
                  />
                  <ReferenceLine y={50} stroke="#64748b" strokeDasharray="4 4" />
                  {Object.entries(METHOD_COLORS).map(([name, color]) => (
                    <Line
                      key={name}
                      type="monotone"
                      dataKey={name}
                      stroke={color}
                      strokeWidth={name === 'Mass-matched BTFR' ? 3 : 1.5}
                      dot={name === 'Mass-matched BTFR'}
                      strokeDasharray={name === 'Mass-matched BTFR' ? undefined : '4 4'}
                    />
                  ))}
                  <Legend wrapperStyle={{ fontSize: '10px' }} />
                </LineChart>
              </ResponsiveContainer>
            </div>
          </GlassCard>
        </div>

        <GlassCard>
          <h3 className="text-base font-semibold text-white mb-2">Why Mass-Matching Wins</h3>
          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div className="bg-red-500/5 border border-red-500/20 rounded-xl p-4">
              <h4 className="text-sm font-bold text-red-400 mb-2">The Problem</h4>
              <p className="text-xs text-slate-300">
                At high z, mass uncertainty ({"\u03C3"}=0.35 dex) is 2.5{"\u00D7"} larger
                than the velocity signal (15.7%). Different mass distributions between
                z=0 and z=1 samples create a false offset that drowns the real shift.
              </p>
            </div>
            <div className="bg-cyan-500/5 border border-cyan-500/20 rounded-xl p-4">
              <h4 className="text-sm font-bold text-cyan-400 mb-2">The Solution</h4>
              <p className="text-xs text-slate-300">
                Bin galaxies by mass. Within each bin, mass scatter is controlled.
                The only variable left is V_flat{"\u2014"}exactly what we want to measure.
                Compare each bin{"'"}s V_flat(z{">"}0) to SPARC V_flat(z=0) at the same mass.
              </p>
            </div>
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-4">
              <h4 className="text-sm font-bold text-emerald-400 mb-2">The Result</h4>
              <p className="text-xs text-slate-300">
                20 mass-matched galaxies {">"} 100 unmatched galaxies.
                Mass matching provides a 4.5{"\u00D7"} boost in detection power,
                equivalent to increasing the sample by 20{"\u00D7"}.
              </p>
            </div>
          </div>
        </GlassCard>

        <GlassCard glow="amber">
          <h3 className="text-base font-semibold text-white mb-5">Recommended 5-Step Pipeline</h3>
          <div className="space-y-4">
            {pipelineSteps.map((s, i) => (
              <div key={i} className="flex gap-4 items-start">
                <div className="flex flex-col items-center">
                  <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-amber-500 to-orange-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">
                    {s.step}
                  </div>
                  {i < pipelineSteps.length - 1 && (
                    <div className="w-0.5 h-6 bg-amber-500/20 mt-1" />
                  )}
                </div>
                <div className="flex-1 pb-2">
                  <div className="flex items-center gap-2">
                    <s.icon className="w-4 h-4 text-amber-400" />
                    <h4 className="text-sm font-bold text-white">{s.title}</h4>
                  </div>
                  <p className="text-xs text-slate-400 mt-1">{s.detail}</p>
                </div>
              </div>
            ))}
          </div>
        </GlassCard>

        <GlassCard glow="none">
          <h3 className="text-base font-semibold text-white mb-3">The Golden Rule</h3>
          <div className="bg-gradient-to-r from-amber-500/10 to-cyan-500/10 border border-amber-500/20 rounded-xl p-5 text-center">
            <p className="text-lg font-bold text-amber-400 italic">
              {"\u201C"}Compare like with like. Control the mass, and the velocity shift reveals itself.{"\u201D"}
            </p>
          </div>
          <div className="mt-4 grid grid-cols-1 md:grid-cols-2 gap-4">
            <div>
              <h4 className="text-sm font-bold text-emerald-400 mb-2">Lessons Learned</h4>
              <ul className="text-xs text-slate-300 space-y-1.5">
                <li className="flex gap-2"><CheckCircle2 className="w-3.5 h-3.5 text-emerald-400 mt-0.5 shrink-0" /><span>Mass matching is the single most important strategy</span></li>
                <li className="flex gap-2"><CheckCircle2 className="w-3.5 h-3.5 text-emerald-400 mt-0.5 shrink-0" /><span>V_flat is better than full rotation curves (less noise)</span></li>
                <li className="flex gap-2"><CheckCircle2 className="w-3.5 h-3.5 text-emerald-400 mt-0.5 shrink-0" /><span>20 matched galaxies {">"} 100 unmatched galaxies</span></li>
              </ul>
            </div>
            <div>
              <h4 className="text-sm font-bold text-red-400 mb-2">What Didn{"'"}t Work</h4>
              <ul className="text-xs text-slate-300 space-y-1.5">
                <li className="flex gap-2"><XCircle className="w-3.5 h-3.5 text-red-400 mt-0.5 shrink-0" /><span>Direct a{"\u2080"} fitting has {"\u223C"}37% bias from mass errors</span></li>
                <li className="flex gap-2"><XCircle className="w-3.5 h-3.5 text-red-400 mt-0.5 shrink-0" /><span>V/V_bar ratio doesn{"'"}t cancel enough systematics</span></li>
                <li className="flex gap-2"><AlertTriangle className="w-3.5 h-3.5 text-amber-400 mt-0.5 shrink-0" /><span>Combined method needs bias calibration work</span></li>
              </ul>
            </div>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
