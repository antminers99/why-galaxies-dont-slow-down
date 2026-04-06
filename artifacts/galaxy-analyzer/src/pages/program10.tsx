import React, { useState, useEffect } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Cell, ReferenceLine } from 'recharts';
import { Telescope, CheckCircle2, XCircle, AlertTriangle, Globe, Search, BarChart3, Crosshair } from 'lucide-react';

interface P10Data {
  program: number;
  phase: string;
  title: string;
  timestamp: string;
  keyFinding: string;
  currentSample: {
    N: number;
    survey: string;
    r: number;
    pValue: number;
    spearmanRho: number;
    partialR_givenVflat: number;
    bootstrap95CI: [number, number];
  };
  surveyOverlaps: Record<string, any>;
  expansionPath: {
    maxCurrentPublic: number;
    nearTerm: string[];
    proposalTargets: Array<{ name: string; DQ: number; Vflat: number; dist: number }>;
  };
  results: Array<{
    name: string;
    DQ: number;
    Vflat: number;
    dist: number;
    m2Power: number;
    m1Power: number;
    m0Power: number;
    m2m0: number;
    paTwist: number;
    coherenceOuter: number;
    nValid: number;
    survey: string;
    source: string;
  }>;
}

const StatBox = ({ label, value, sub, color = 'cyan' }: { label: string; value: string; sub?: string; color?: string }) => (
  <div className="text-center p-3">
    <div className="text-xs text-slate-500 uppercase tracking-wider mb-1">{label}</div>
    <div className={"text-xl font-mono font-bold " + (color === 'cyan' ? 'text-cyan-400' : color === 'amber' ? 'text-amber-400' : color === 'emerald' ? 'text-emerald-400' : color === 'violet' ? 'text-violet-400' : color === 'red' ? 'text-red-400' : 'text-slate-300')}>{value}</div>
    {sub && <div className="text-xs text-slate-500 mt-0.5">{sub}</div>}
  </div>
);

export default function Program10Page() {
  const [data, setData] = useState<P10Data | null>(null);
  const [loading, setLoading] = useState(true);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'program10-multi-survey.json')
      .then(r => r.json())
      .then(d => { setData(d); setLoading(false); })
      .catch(() => setLoading(false));
  }, []);

  if (loading) return <Layout title="Program 10"><div className="text-center py-20 text-slate-400">Loading...</div></Layout>;
  if (!data) return <Layout title="Program 10"><div className="text-center py-20 text-red-400">No Program 10 data available</div></Layout>;

  const scatterData = data.results.map(r => ({
    name: r.name,
    DQ: r.DQ,
    logM2: Math.log10(r.m2Power),
    Vflat: r.Vflat,
    survey: r.survey,
  }));

  const surveys = [
    { name: 'MaNGA', confirmed: data.surveyOverlaps.MaNGA?.confirmed || 0, status: 'zero', reason: data.surveyOverlaps.MaNGA?.reason || '' },
    { name: 'THINGS', confirmed: data.surveyOverlaps.THINGS?.confirmed || 0, status: 'complete', reason: 'Fully processed in Program 9' },
    { name: 'PHANGS', confirmed: data.surveyOverlaps.PHANGS?.confirmed || 0, status: 'overlap', reason: 'All overlap with THINGS sample' },
    { name: 'CALIFA', confirmed: data.surveyOverlaps.CALIFA?.confirmed || 0, status: 'partial', reason: 'Needs verification; new galaxies lack quality fits' },
    { name: 'LITTLE THINGS', confirmed: data.surveyOverlaps.LITTLE_THINGS?.confirmed || 0, status: 'no-fit', reason: 'Galaxies lack quality RAR a0 fits' },
  ];

  return (
    <Layout title="Program 10 — Survey Census & Expansion Path">
      <div className="space-y-6">

        <GlassCard>
          <div className="flex items-center gap-3 mb-4">
            <Telescope className="w-6 h-6 text-violet-400" />
            <h2 className="text-xl font-display font-bold text-white">Program 10.1: Definitive Survey Census & Expansion Path</h2>
          </div>
          <div className="bg-amber-500/10 border border-amber-500/30 rounded-lg p-4 mb-4">
            <div className="flex items-start gap-2">
              <AlertTriangle className="w-5 h-5 text-amber-400 mt-0.5 flex-shrink-0" />
              <div>
                <div className="text-amber-300 font-semibold text-sm">Key Finding</div>
                <div className="text-amber-200/80 text-sm mt-1">{data.keyFinding}</div>
              </div>
            </div>
          </div>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-2 bg-slate-800/50 rounded-lg">
            <StatBox label="Sample Size" value={"N = " + data.currentSample.N} sub="THINGS only" color="cyan" />
            <StatBox label="r(DQ, log m2)" value={data.currentSample.r.toFixed(3)} sub={"p = " + data.currentSample.pValue.toFixed(4)} color="emerald" />
            <StatBox label="Partial r | Vflat" value={data.currentSample.partialR_givenVflat.toFixed(3)} sub="Not mass-driven" color="violet" />
            <StatBox label="Bootstrap 95% CI" value={"[" + data.currentSample.bootstrap95CI[0].toFixed(2) + ", " + data.currentSample.bootstrap95CI[1].toFixed(2) + "]"} sub="Excludes zero" color="amber" />
          </div>
        </GlassCard>

        <div className="grid md:grid-cols-2 gap-6">
          <GlassCard>
            <h3 className="text-lg font-display font-bold text-white mb-4 flex items-center gap-2">
              <BarChart3 className="w-5 h-5 text-cyan-400" />
              DQ vs log(m=2 Power)
            </h3>
            <ResponsiveContainer width="100%" height={300}>
              <ScatterChart margin={{ top: 10, right: 20, bottom: 30, left: 10 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="#334155" />
                <XAxis type="number" dataKey="DQ" name="DQ" stroke="#94a3b8" fontSize={11}
                  label={{ value: 'DQ (bilateral residual)', position: 'bottom', offset: 15, fill: '#94a3b8', fontSize: 11 }} />
                <YAxis type="number" dataKey="logM2" name="log(m2)" stroke="#94a3b8" fontSize={11}
                  label={{ value: 'log(m=2 power)', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 11 }} />
                <ReferenceLine x={0} stroke="#475569" strokeDasharray="4 4" />
                <Tooltip
                  content={({ active, payload }) => {
                    if (!active || !payload?.length) return null;
                    const d = payload[0].payload;
                    return (
                      <div className="bg-slate-800 border border-slate-600 rounded p-2 text-xs">
                        <div className="text-white font-bold">{d.name}</div>
                        <div className="text-slate-300">DQ: {d.DQ.toFixed(2)}</div>
                        <div className="text-cyan-400">log(m2): {d.logM2.toFixed(3)}</div>
                        <div className="text-slate-400">Vflat: {d.Vflat} km/s</div>
                      </div>
                    );
                  }}
                />
                <Scatter data={scatterData} fill="#22d3ee">
                  {scatterData.map((entry, i) => (
                    <Cell key={i} fill={entry.DQ > 0 ? '#f59e0b' : '#6366f1'} />
                  ))}
                </Scatter>
              </ScatterChart>
            </ResponsiveContainer>
            <div className="text-xs text-slate-500 mt-2 text-center">
              r = {data.currentSample.r.toFixed(3)}, p = {data.currentSample.pValue.toFixed(4)} | Gold = high-H, Purple = low-H
            </div>
          </GlassCard>

          <GlassCard>
            <h3 className="text-lg font-display font-bold text-white mb-4 flex items-center gap-2">
              <Globe className="w-5 h-5 text-emerald-400" />
              Survey Overlap Verification
            </h3>
            <div className="space-y-3">
              {surveys.map(s => (
                <div key={s.name} className="flex items-start gap-3 p-2 bg-slate-800/30 rounded">
                  {s.status === 'complete' ? <CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" /> :
                   s.status === 'zero' ? <XCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" /> :
                   <AlertTriangle className="w-4 h-4 text-amber-400 mt-0.5 flex-shrink-0" />}
                  <div className="flex-1 min-w-0">
                    <div className="flex items-center gap-2">
                      <span className="text-white font-semibold text-sm">{s.name}</span>
                      <span className="text-xs font-mono text-slate-400">{s.confirmed} galaxies</span>
                    </div>
                    <div className="text-xs text-slate-500 mt-0.5">{s.reason}</div>
                  </div>
                </div>
              ))}
            </div>
            <div className="mt-3 p-2 bg-red-500/10 border border-red-500/20 rounded text-xs text-red-300">
              MaNGA correction: Previous census used target catalog (candidates), not observed catalog. Verified against SDSS DR17 mangaDrpAll.
            </div>
          </GlassCard>
        </div>

        <GlassCard>
          <h3 className="text-lg font-display font-bold text-white mb-4 flex items-center gap-2">
            <Search className="w-5 h-5 text-amber-400" />
            Processed Galaxies
          </h3>
          <div className="overflow-x-auto">
            <table className="w-full text-xs font-mono">
              <thead>
                <tr className="border-b border-white/10 text-slate-400">
                  <th className="text-left py-2 px-2">Galaxy</th>
                  <th className="text-right py-2 px-2">DQ</th>
                  <th className="text-right py-2 px-2">log(m2)</th>
                  <th className="text-right py-2 px-2">m2/m0</th>
                  <th className="text-right py-2 px-2">PA twist</th>
                  <th className="text-right py-2 px-2">Coherence</th>
                  <th className="text-right py-2 px-2">Vflat</th>
                  <th className="text-right py-2 px-2">D (Mpc)</th>
                  <th className="text-right py-2 px-2">Survey</th>
                </tr>
              </thead>
              <tbody>
                {data.results.sort((a, b) => b.DQ - a.DQ).map(r => (
                  <tr key={r.name} className="border-b border-white/5 hover:bg-white/5 transition-colors">
                    <td className={"py-1.5 px-2 " + (r.DQ > 0 ? 'text-amber-400' : 'text-indigo-400')}>{r.name}</td>
                    <td className="py-1.5 px-2 text-right text-white font-bold">{r.DQ.toFixed(2)}</td>
                    <td className="py-1.5 px-2 text-right text-cyan-400">{Math.log10(r.m2Power).toFixed(3)}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{r.m2m0.toFixed(3)}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{r.paTwist.toFixed(1) + "\u00B0"}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{(r.coherenceOuter * 100).toFixed(0) + "%"}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{Math.round(r.Vflat)}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{r.dist.toFixed(1)}</td>
                    <td className="py-1.5 px-2 text-right text-emerald-400">{r.survey}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </GlassCard>

        <GlassCard>
          <h3 className="text-lg font-display font-bold text-white mb-4 flex items-center gap-2">
            <Crosshair className="w-5 h-5 text-red-400" />
            Priority Targets for Future Observations
          </h3>
          <p className="text-sm text-slate-400 mb-3">
            Top SPARC galaxies without 2D kinematic data, ranked by prediction power (|DQ|).
            These are the highest-value targets for dedicated VLA/ASKAP proposals.
          </p>
          <div className="overflow-x-auto">
            <table className="w-full text-xs font-mono">
              <thead>
                <tr className="border-b border-white/10 text-slate-400">
                  <th className="text-left py-2 px-2">Rank</th>
                  <th className="text-left py-2 px-2">Galaxy</th>
                  <th className="text-right py-2 px-2">DQ</th>
                  <th className="text-right py-2 px-2">Vflat</th>
                  <th className="text-right py-2 px-2">D (Mpc)</th>
                  <th className="text-right py-2 px-2">Prediction</th>
                </tr>
              </thead>
              <tbody>
                {data.expansionPath.proposalTargets.slice(0, 10).map((t, i) => (
                  <tr key={t.name} className="border-b border-white/5 hover:bg-white/5 transition-colors">
                    <td className="py-1.5 px-2 text-slate-500">{i + 1}</td>
                    <td className={"py-1.5 px-2 " + (t.DQ > 0 ? 'text-amber-400' : 'text-indigo-400')}>{t.name}</td>
                    <td className="py-1.5 px-2 text-right text-white font-bold">{t.DQ.toFixed(2)}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{Math.round(t.Vflat)}</td>
                    <td className="py-1.5 px-2 text-right text-slate-300">{t.dist.toFixed(1)}</td>
                    <td className={"py-1.5 px-2 text-right font-semibold " + (t.DQ > 0 ? 'text-amber-400' : 'text-indigo-400')}>
                      {t.DQ > 0 ? 'Strong m=2' : 'Weak m=2'}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
          <div className="mt-3 grid md:grid-cols-3 gap-3">
            <div className="p-2 bg-emerald-500/10 border border-emerald-500/20 rounded text-xs text-emerald-300">
              <strong>N needed for alpha=0.01:</strong> 11 galaxies
            </div>
            <div className="p-2 bg-amber-500/10 border border-amber-500/20 rounded text-xs text-amber-300">
              <strong>N needed for alpha=0.001:</strong> 15 galaxies
            </div>
            <div className="p-2 bg-violet-500/10 border border-violet-500/20 rounded text-xs text-violet-300">
              <strong>Near-term:</strong> WALLABY (ASKAP), PHANGS expansion
            </div>
          </div>
        </GlassCard>

        <GlassCard>
          <h3 className="text-lg font-display font-bold text-white mb-3">Robustness Summary</h3>
          <div className="grid md:grid-cols-2 gap-4">
            <div className="space-y-2 text-sm">
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <span className="text-slate-300">Permutation p {"<"} 0.005 (50k replicates)</span>
              </div>
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <span className="text-slate-300">LOO: 7/7 positive, min r = 0.691</span>
              </div>
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <span className="text-slate-300">Bootstrap 95% CI excludes zero</span>
              </div>
              <div className="flex items-center gap-2">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <span className="text-slate-300">Partial r(DQ, m2 | Vflat) = {data.currentSample.partialR_givenVflat.toFixed(3)}</span>
              </div>
            </div>
            <div className="space-y-2 text-sm">
              <div className="flex items-center gap-2">
                <AlertTriangle className="w-4 h-4 text-amber-400" />
                <span className="text-slate-300">Small N = 7 limits statistical power</span>
              </div>
              <div className="flex items-center gap-2">
                <AlertTriangle className="w-4 h-4 text-amber-400" />
                <span className="text-slate-300">Single survey (THINGS) only</span>
              </div>
              <div className="flex items-center gap-2">
                <AlertTriangle className="w-4 h-4 text-amber-400" />
                <span className="text-slate-300">Cannot reach N=20 with public data</span>
              </div>
              <div className="flex items-center gap-2">
                <AlertTriangle className="w-4 h-4 text-amber-400" />
                <span className="text-slate-300">Independent survey verification pending</span>
              </div>
            </div>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
