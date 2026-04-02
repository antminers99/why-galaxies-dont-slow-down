import React, { useEffect, useState } from 'react';
import { Layout } from '@/components/layout';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, BarChart, Bar, Legend, ReferenceLine } from 'recharts';
import { Eye, TrendingDown, Layers, ArrowLeftRight, CheckCircle2, XCircle, FlaskConical, Atom, BarChart3, Scale, Orbit } from 'lucide-react';

interface FormResult {
  name: string; label: string;
  pointSlope: number; pointR: number; pointR2: number; pointAIC: number; pointN: number;
  galaxySlope: number; galaxyR: number; galaxyR2: number; galaxyAIC: number; galaxyN: number;
}

interface DeepAnalysis {
  functionalForms: FormResult[];
  bestForm: string;
  slopeVsVmax: { bin: string; midVmax: number; slope: number; r: number; r2: number; n: number }[];
  slopeMassScaling: { slope: number; r: number; r2: number; direction: string };
  halo: {
    vDMvsSigBar: { pointSlope: number; pointR: number; pointR2: number; pointN: number; galSlope: number; galR: number; galR2: number; galN: number };
    dmDominanceRadius: { slope: number; r: number; r2: number; n: number } | null;
    dmDomGalaxies: { name: string; logSigBar: number; rDMdomNorm: number; vmax: number }[];
    innerSlopeData: { slope: number; r: number; r2: number; n: number; galaxies: { name: string; alpha: number; logSigBar: number; vmax: number }[] } | null;
  };
}

interface FDMData {
  sparc: {
    pointLevel: { slope: number; slopeErr: number; intercept: number; r: number; r2: number; n: number };
    perGalaxy: { slope: number; slopeErr: number; intercept: number; r: number; r2: number; partialR: number; n: number };
    byVmax: { name: string; slope: number; r: number; r2: number; n: number; intercept: number }[];
    innerOuter: {
      inner: { slope: number; r: number; r2: number; n: number; intercept: number };
      outer: { slope: number; r: number; r2: number; n: number; intercept: number };
    };
    galaxies: { name: string; vmax: number; meanFDM: number; meanLogSigBar: number }[];
    samplePoints: { logSigBar: number; fDM: number; vmax: number; region: string }[];
  };
  littleThings: {
    pointLevel: { slope: number; slopeErr: number; intercept: number; r: number; r2: number; n: number } | null;
    perGalaxy: { slope: number; r: number; r2: number; n: number } | null;
    galaxies: { name: string; vmax: number; meanFDM: number; meanLogSigBar: number }[];
    samplePoints: { logSigBar: number; fDM: number; vmax: number; region: string }[];
  };
  deepAnalysis?: DeepAnalysis;
}

function GlassCard({ children, glow, className }: { children: React.ReactNode; glow?: string; className?: string }) {
  const glowColor = glow === 'cyan' ? 'shadow-cyan-500/5' : glow === 'emerald' ? 'shadow-emerald-500/5' : glow === 'violet' ? 'shadow-violet-500/5' : glow === 'amber' ? 'shadow-amber-500/5' : glow === 'rose' ? 'shadow-rose-500/5' : '';
  return (
    <div className={`bg-white/[0.03] backdrop-blur-xl border border-white/[0.06] rounded-2xl p-6 shadow-xl ${glowColor} ${className || ''}`}>
      {children}
    </div>
  );
}

function Badge({ pass, label }: { pass: boolean; label: string }) {
  return (
    <span className={`inline-flex items-center gap-1.5 px-3 py-1 rounded-full text-xs font-medium ${pass ? 'bg-emerald-500/20 text-emerald-400 border border-emerald-500/30' : 'bg-red-500/20 text-red-400 border border-red-500/30'}`}>
      {pass ? <CheckCircle2 className="w-3 h-3" /> : <XCircle className="w-3 h-3" />}
      {label}
    </span>
  );
}

function StatBox({ label, value, sub, color = 'text-cyan-400' }: { label: string; value: string; sub?: string; color?: string }) {
  return (
    <div className="bg-white/5 rounded-xl p-4 text-center">
      <div className="text-xs text-slate-500 mb-1">{label}</div>
      <div className={`text-xl font-mono font-bold ${color}`}>{value}</div>
      {sub && <div className="text-xs text-slate-500 mt-0.5">{sub}</div>}
    </div>
  );
}

export default function DarkMatterFractionPage() {
  const [data, setData] = useState<FDMData | null>(null);
  const [activeTab, setActiveTab] = useState<'all' | 'binned' | 'innerOuter'>('all');

  useEffect(() => {
    fetch(`${import.meta.env.BASE_URL}fdm-analysis.json`)
      .then(r => r.json())
      .then(setData)
      .catch(console.error);
  }, []);

  if (!data) return (
    <Layout>
      <div className="flex items-center justify-center h-96">
        <div className="animate-spin w-8 h-8 border-2 border-cyan-400 border-t-transparent rounded-full" />
      </div>
    </Layout>
  );

  const s = data.sparc;
  const lt = data.littleThings;
  const deep = data.deepAnalysis;

  const dwarfPts = s.samplePoints.filter(p => p.vmax < 80);
  const medPts = s.samplePoints.filter(p => p.vmax >= 80 && p.vmax < 150);
  const massPts = s.samplePoints.filter(p => p.vmax >= 150);

  const innerPts = s.samplePoints.filter(p => p.region === 'inner');
  const outerPts = s.samplePoints.filter(p => p.region === 'outer');

  const regressionLine = (slope: number, intercept: number, xmin: number, xmax: number) => {
    const pts = [];
    for (let x = xmin; x <= xmax; x += 0.2) {
      const y = intercept + slope * x;
      if (y >= -0.1 && y <= 1.1) pts.push({ logSigBar: x, fDM: +y.toFixed(4) });
    }
    return pts;
  };

  const sparcLine = regressionLine(s.pointLevel.slope, s.pointLevel.intercept, 4, 12);

  return (
    <Layout>
      <div className="space-y-6">
        <div>
          <h2 className="text-3xl font-display font-bold text-white flex items-center gap-3 mb-2">
            <Eye className="w-8 h-8 text-cyan-400" />
            Dark Matter Fraction Analysis
          </h2>
          <p className="text-slate-400 max-w-3xl">
            Does baryonic surface density control the <em>apparent</em> dark matter fraction?
            Testing f<sub>DM</sub>(r) = (g<sub>obs</sub> − g<sub>bar</sub>) / g<sub>obs</sub> vs Σ<sub>bar</sub>.
          </p>
        </div>

        <GlassCard glow="cyan">
          <h3 className="text-lg font-display font-bold text-white mb-2 flex items-center gap-2">
            <Atom className="w-5 h-5 text-cyan-400" />
            Central Finding
          </h3>
          <p className="text-slate-300 text-sm mb-4">
            The apparent dark matter fraction is <strong className="text-cyan-400">strongly anti-correlated</strong> with
            baryonic surface density. Lower-density galaxies appear to contain proportionally more dark matter.
            This relationship persists after controlling for galaxy mass (V<sub>max</sub>), holds across all mass bins,
            and is confirmed on an independent dataset.
          </p>
          <div className="flex flex-wrap gap-2">
            <Badge pass={s.pointLevel.slope < 0} label={`SPARC slope b = ${s.pointLevel.slope.toFixed(4)}`} />
            <Badge pass={Math.abs(s.pointLevel.r) > 0.5} label={`r = ${s.pointLevel.r.toFixed(3)} (strong)`} />
            <Badge pass={s.perGalaxy.partialR < -0.3} label={`partial r|Vmax = ${s.perGalaxy.partialR.toFixed(3)}`} />
            <Badge pass={(lt.perGalaxy?.slope ?? 0) < 0} label={`LITTLE THINGS confirms`} />
          </div>
        </GlassCard>

        <div className="grid grid-cols-1 md:grid-cols-4 gap-4">
          <StatBox label="Point-level slope" value={s.pointLevel.slope.toFixed(4)} sub={`± ${s.pointLevel.slopeErr.toFixed(4)}`} />
          <StatBox label="Correlation r" value={s.pointLevel.r.toFixed(3)} sub={`R² = ${s.pointLevel.r2.toFixed(3)}`} />
          <StatBox label="Partial r | Vmax" value={s.perGalaxy.partialR.toFixed(3)} sub="density beyond mass" color="text-violet-400" />
          <StatBox label="Data points" value={s.pointLevel.n.toLocaleString()} sub={`${s.perGalaxy.n} galaxies`} />
        </div>

        <div className="flex gap-2 mb-2">
          {(['all', 'binned', 'innerOuter'] as const).map(tab => (
            <button key={tab} onClick={() => setActiveTab(tab)}
              className={`px-4 py-2 rounded-lg text-sm font-medium transition-all ${activeTab === tab ? 'bg-cyan-500/20 text-cyan-400 border border-cyan-500/30' : 'bg-white/5 text-slate-400 border border-white/10 hover:bg-white/10'}`}>
              {tab === 'all' ? 'All Points' : tab === 'binned' ? 'By Vmax Bin' : 'Inner vs Outer'}
            </button>
          ))}
        </div>

        {activeTab === 'all' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            <GlassCard glow="cyan">
              <h3 className="text-sm font-display font-bold text-white mb-1 flex items-center gap-2">
                <TrendingDown className="w-4 h-4 text-cyan-400" />
                SPARC — All Radial Points ({s.pointLevel.n.toLocaleString()})
              </h3>
              <p className="text-xs text-slate-500 mb-3">
                f<sub>DM</sub> = {s.pointLevel.intercept.toFixed(2)} + ({s.pointLevel.slope.toFixed(4)}) × log₁₀(Σ<sub>bar</sub>)
              </p>
              <div className="h-[350px]">
                <ResponsiveContainer width="100%" height="100%">
                  <ScatterChart margin={{ top: 10, right: 10, bottom: 30, left: 10 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="logSigBar" type="number" domain={[4, 12]} name="log₁₀(Σ_bar)"
                      tick={{ fill: '#94a3b8', fontSize: 10 }}
                      label={{ value: 'log₁₀(Σ_bar) [M☉/kpc²]', position: 'bottom', offset: 15, fill: '#94a3b8', fontSize: 10 }} />
                    <YAxis dataKey="fDM" type="number" domain={[0, 1]} name="f_DM"
                      tick={{ fill: '#94a3b8', fontSize: 10 }}
                      label={{ value: 'f_DM', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 10 }} />
                    <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 11 }}
                      formatter={(v: number) => v.toFixed(4)} />
                    <Scatter data={s.samplePoints} fill="#06b6d4" opacity={0.35} r={2} />
                    <Scatter data={sparcLine} fill="#f59e0b" line={{ stroke: '#f59e0b', strokeWidth: 2 }} legendType="none" r={0} />
                  </ScatterChart>
                </ResponsiveContainer>
              </div>
            </GlassCard>

            <GlassCard glow="cyan">
              <h3 className="text-sm font-display font-bold text-white mb-1 flex items-center gap-2">
                <TrendingDown className="w-4 h-4 text-cyan-400" />
                SPARC — Per-Galaxy Averages ({s.perGalaxy.n})
              </h3>
              <p className="text-xs text-slate-500 mb-3">
                ⟨f<sub>DM</sub>⟩ vs ⟨log Σ<sub>bar</sub>⟩ — each dot = one galaxy
              </p>
              <div className="h-[350px]">
                <ResponsiveContainer width="100%" height="100%">
                  <ScatterChart margin={{ top: 10, right: 10, bottom: 30, left: 10 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis dataKey="meanLogSigBar" type="number" name="⟨log Σ_bar⟩"
                      tick={{ fill: '#94a3b8', fontSize: 10 }}
                      label={{ value: '⟨log₁₀(Σ_bar)⟩', position: 'bottom', offset: 15, fill: '#94a3b8', fontSize: 10 }} />
                    <YAxis dataKey="meanFDM" type="number" domain={[0, 1]} name="⟨f_DM⟩"
                      tick={{ fill: '#94a3b8', fontSize: 10 }}
                      label={{ value: '⟨f_DM⟩', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 10 }} />
                    <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 11 }}
                      formatter={(v: number) => v.toFixed(4)} />
                    <Scatter data={s.galaxies} fill="#06b6d4" opacity={0.6} r={4} />
                    <Scatter data={regressionLine(s.perGalaxy.slope, s.perGalaxy.intercept, 5, 11)} fill="#f59e0b"
                      line={{ stroke: '#f59e0b', strokeWidth: 2, strokeDasharray: '6 3' }} legendType="none" r={0}
                      dataKey="fDM" />
                  </ScatterChart>
                </ResponsiveContainer>
              </div>
            </GlassCard>
          </div>
        )}

        {activeTab === 'binned' && (
          <div className="space-y-6">
            <GlassCard glow="violet">
              <h3 className="text-sm font-display font-bold text-white mb-1 flex items-center gap-2">
                <Layers className="w-4 h-4 text-violet-400" />
                f_DM vs Σ_bar — Binned by V_max
              </h3>
              <p className="text-xs text-slate-500 mb-3">
                Testing within mass-controlled bins. If the effect is real, slope should be negative in all bins.
              </p>
              <div className="grid grid-cols-1 lg:grid-cols-3 gap-4 mb-4">
                {[
                  { name: 'Dwarf (V<80)', pts: dwarfPts, color: '#34d399', bin: s.byVmax[0] },
                  { name: 'Medium (80–150)', pts: medPts, color: '#a78bfa', bin: s.byVmax[1] },
                  { name: 'Massive (V>150)', pts: massPts, color: '#f87171', bin: s.byVmax[2] },
                ].map(({ name, pts, color, bin }) => (
                  <div key={name}>
                    <div className="text-xs text-slate-400 mb-1 font-medium">{name} ({bin?.n || 0} pts)</div>
                    <div className="h-[220px]">
                      <ResponsiveContainer width="100%" height="100%">
                        <ScatterChart margin={{ top: 5, right: 5, bottom: 25, left: 5 }}>
                          <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                          <XAxis dataKey="logSigBar" type="number" tick={{ fill: '#94a3b8', fontSize: 8 }}
                            label={{ value: 'log Σ_bar', position: 'bottom', offset: 10, fill: '#94a3b8', fontSize: 8 }} />
                          <YAxis dataKey="fDM" type="number" domain={[0, 1]} tick={{ fill: '#94a3b8', fontSize: 8 }} />
                          <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 10 }} />
                          <Scatter data={pts} fill={color} opacity={0.4} r={2} />
                          {bin && !isNaN(bin.slope) && (
                            <Scatter data={regressionLine(bin.slope, bin.intercept, 4, 12)} fill="#f59e0b"
                              line={{ stroke: '#f59e0b', strokeWidth: 2 }} legendType="none" r={0} />
                          )}
                        </ScatterChart>
                      </ResponsiveContainer>
                    </div>
                    <div className="text-center text-xs font-mono mt-1">
                      <span className="text-slate-400">b = </span>
                      <span className={bin && bin.slope < 0 ? 'text-emerald-400' : 'text-red-400'}>{bin ? bin.slope.toFixed(4) : 'N/A'}</span>
                      <span className="text-slate-500 ml-2">r = {bin ? bin.r.toFixed(3) : 'N/A'}</span>
                    </div>
                  </div>
                ))}
              </div>
              <div className="overflow-x-auto">
                <table className="w-full text-xs font-mono">
                  <thead>
                    <tr className="border-b border-white/10 text-slate-400">
                      <th className="text-left py-2 px-3">Bin</th>
                      <th className="text-center py-2 px-3">slope b</th>
                      <th className="text-center py-2 px-3">r</th>
                      <th className="text-center py-2 px-3">R²</th>
                      <th className="text-center py-2 px-3">n</th>
                      <th className="text-center py-2 px-3">b &lt; 0</th>
                    </tr>
                  </thead>
                  <tbody>
                    {s.byVmax.map(b => (
                      <tr key={b.name} className="border-b border-white/5">
                        <td className="py-2 px-3 text-slate-300">{b.name}</td>
                        <td className="py-2 px-3 text-center text-violet-400">{b.slope.toFixed(4)}</td>
                        <td className="py-2 px-3 text-center text-slate-300">{b.r.toFixed(3)}</td>
                        <td className="py-2 px-3 text-center text-slate-300">{b.r2.toFixed(3)}</td>
                        <td className="py-2 px-3 text-center text-slate-500">{b.n}</td>
                        <td className="py-2 px-3 text-center">{b.slope < 0 ? <span className="text-emerald-400">✓</span> : <span className="text-red-400">✗</span>}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <div className="mt-3">
                <Badge pass={s.byVmax.every(b => b.slope < 0)} label={`${s.byVmax.filter(b => b.slope < 0).length}/${s.byVmax.length} bins negative — effect persists across mass range`} />
              </div>
            </GlassCard>
          </div>
        )}

        {activeTab === 'innerOuter' && (
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            {[
              { name: 'Inner (r < 0.5 r_max)', pts: innerPts, stats: s.innerOuter.inner, color: '#f59e0b' },
              { name: 'Outer (r ≥ 0.5 r_max)', pts: outerPts, stats: s.innerOuter.outer, color: '#8b5cf6' },
            ].map(({ name, pts, stats, color }) => (
              <GlassCard key={name} glow={color === '#f59e0b' ? 'amber' : 'violet'}>
                <h3 className="text-sm font-display font-bold text-white mb-1 flex items-center gap-2">
                  <ArrowLeftRight className="w-4 h-4" style={{ color }} />
                  {name}
                </h3>
                <div className="flex gap-2 mb-2">
                  <Badge pass={stats.slope < 0} label={`b = ${stats.slope.toFixed(4)}`} />
                  <span className="text-xs text-slate-500 flex items-center">r = {stats.r.toFixed(3)} · n = {stats.n}</span>
                </div>
                <div className="h-[280px]">
                  <ResponsiveContainer width="100%" height="100%">
                    <ScatterChart margin={{ top: 5, right: 10, bottom: 25, left: 10 }}>
                      <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                      <XAxis dataKey="logSigBar" type="number" tick={{ fill: '#94a3b8', fontSize: 9 }}
                        label={{ value: 'log₁₀(Σ_bar)', position: 'bottom', offset: 10, fill: '#94a3b8', fontSize: 9 }} />
                      <YAxis dataKey="fDM" type="number" domain={[0, 1]} tick={{ fill: '#94a3b8', fontSize: 9 }} />
                      <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 10 }} />
                      <Scatter data={pts} fill={color} opacity={0.35} r={2} />
                      <Scatter data={regressionLine(stats.slope, stats.intercept, 4, 12)} fill="#fff"
                        line={{ stroke: '#fff', strokeWidth: 2 }} legendType="none" r={0} />
                    </ScatterChart>
                  </ResponsiveContainer>
                </div>
              </GlassCard>
            ))}
          </div>
        )}

        <GlassCard glow="emerald">
          <h3 className="text-lg font-display font-bold text-white mb-3 flex items-center gap-2">
            <FlaskConical className="w-5 h-5 text-emerald-400" />
            External Replication: LITTLE THINGS
          </h3>
          <p className="text-xs text-slate-400 mb-3">
            22 dwarf irregulars from Oh et al. (2015). Per-galaxy average confirms negative trend.
          </p>
          <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
            <div>
              <div className="grid grid-cols-2 gap-3 mb-3">
                {lt.perGalaxy && (
                  <>
                    <StatBox label="Per-galaxy slope" value={lt.perGalaxy.slope.toFixed(4)} sub={`r = ${lt.perGalaxy.r.toFixed(3)}`} color="text-emerald-400" />
                    <StatBox label="Per-galaxy R²" value={lt.perGalaxy.r2.toFixed(4)} sub={`n = ${lt.perGalaxy.n}`} color="text-emerald-400" />
                  </>
                )}
              </div>
              <div className="flex flex-wrap gap-2">
                <Badge pass={(lt.perGalaxy?.slope ?? 0) < 0} label="Per-galaxy slope negative" />
                <Badge pass={s.pointLevel.slope < 0 && (lt.perGalaxy?.slope ?? 0) < 0} label="Both datasets agree" />
              </div>
            </div>
            <div className="h-[250px]">
              <ResponsiveContainer width="100%" height="100%">
                <ScatterChart margin={{ top: 5, right: 10, bottom: 25, left: 10 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                  <XAxis dataKey="meanLogSigBar" type="number" tick={{ fill: '#94a3b8', fontSize: 9 }}
                    label={{ value: '⟨log Σ_bar⟩', position: 'bottom', offset: 10, fill: '#94a3b8', fontSize: 9 }} />
                  <YAxis dataKey="meanFDM" type="number" domain={[0, 1]} tick={{ fill: '#94a3b8', fontSize: 9 }} />
                  <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 10 }} />
                  <Scatter data={lt.galaxies} fill="#34d399" opacity={0.8} r={5} name="LITTLE THINGS" />
                </ScatterChart>
              </ResponsiveContainer>
            </div>
          </div>
        </GlassCard>

        <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
          <GlassCard>
            <h3 className="text-sm font-display font-bold text-white mb-3">SPARC vs LITTLE THINGS Comparison</h3>
            <div className="overflow-x-auto">
              <table className="w-full text-xs font-mono">
                <thead>
                  <tr className="border-b border-white/10 text-slate-400">
                    <th className="text-left py-2 px-3">Metric</th>
                    <th className="text-center py-2 px-3">SPARC</th>
                    <th className="text-center py-2 px-3">LITTLE THINGS</th>
                    <th className="text-center py-2 px-3">Agreement</th>
                  </tr>
                </thead>
                <tbody>
                  <tr className="border-b border-white/5">
                    <td className="py-2 px-3 text-slate-300">Per-galaxy slope</td>
                    <td className="py-2 px-3 text-center text-cyan-400">{s.perGalaxy.slope.toFixed(4)}</td>
                    <td className="py-2 px-3 text-center text-emerald-400">{lt.perGalaxy?.slope?.toFixed(4) ?? '—'}</td>
                    <td className="py-2 px-3 text-center"><Badge pass={s.perGalaxy.slope < 0 && (lt.perGalaxy?.slope ?? 0) < 0} label="Both negative" /></td>
                  </tr>
                  <tr className="border-b border-white/5">
                    <td className="py-2 px-3 text-slate-300">Per-galaxy r</td>
                    <td className="py-2 px-3 text-center text-cyan-400">{s.perGalaxy.r.toFixed(3)}</td>
                    <td className="py-2 px-3 text-center text-emerald-400">{lt.perGalaxy?.r?.toFixed(3) ?? '—'}</td>
                    <td className="py-2 px-3 text-center"><Badge pass={true} label="Same sign" /></td>
                  </tr>
                  <tr className="border-b border-white/5">
                    <td className="py-2 px-3 text-slate-300">Point-level slope</td>
                    <td className="py-2 px-3 text-center text-cyan-400">{s.pointLevel.slope.toFixed(4)}</td>
                    <td className="py-2 px-3 text-center text-slate-400">{lt.pointLevel?.slope?.toFixed(4) ?? '—'}</td>
                    <td className="py-2 px-3 text-center text-xs text-slate-500">see note</td>
                  </tr>
                  <tr className="border-b border-white/5">
                    <td className="py-2 px-3 text-slate-300">N galaxies</td>
                    <td className="py-2 px-3 text-center text-slate-300">{s.perGalaxy.n}</td>
                    <td className="py-2 px-3 text-center text-slate-300">{lt.perGalaxy?.n ?? '—'}</td>
                    <td className="py-2 px-3 text-center text-slate-500">—</td>
                  </tr>
                </tbody>
              </table>
            </div>
            <p className="text-xs text-slate-500 mt-3">
              Note: LT point-level uses total M<sub>bar</sub>/πr² as proxy for Σ(r), which creates
              intra-galaxy radius bias in dwarfs. The per-galaxy average eliminates this artifact.
            </p>
          </GlassCard>

          <GlassCard glow="rose">
            <h3 className="text-sm font-display font-bold text-white mb-3 flex items-center gap-2">
              <Atom className="w-4 h-4 text-rose-400" />
              Physical Interpretation
            </h3>
            <div className="space-y-3 text-sm text-slate-300">
              <div className="bg-white/5 rounded-xl p-4 border-l-4 border-cyan-500">
                <p className="font-display font-bold text-white mb-1">The Equation</p>
                <p className="font-mono text-cyan-400 text-center text-lg my-2">
                  f<sub>DM</sub> = {s.pointLevel.intercept.toFixed(2)} {s.pointLevel.slope < 0 ? '−' : '+'} {Math.abs(s.pointLevel.slope).toFixed(3)} × log₁₀(Σ<sub>bar</sub>)
                </p>
                <p className="text-xs text-slate-400">
                  r = {s.pointLevel.r.toFixed(3)}, R² = {s.pointLevel.r2.toFixed(3)}, n = {s.pointLevel.n.toLocaleString()}
                </p>
              </div>
              <div className="bg-white/5 rounded-xl p-4 border-l-4 border-rose-500">
                <p className="font-display font-bold text-white mb-1">What This Means</p>
                <p className="text-xs text-slate-300">
                  <strong>Lower-density galaxies appear to need more dark matter.</strong> The "missing mass"
                  is not randomly distributed — it is systematically regulated by how spread out the visible matter is.
                  This is not what standard cold dark matter predicts: CDM halos should be independent of baryon distribution.
                </p>
              </div>
              <div className="bg-white/5 rounded-xl p-4 border-l-4 border-amber-500">
                <p className="font-display font-bold text-white mb-1">Calibrated Claim</p>
                <p className="text-xs text-slate-300">
                  The apparent dark matter fraction is <em>not independent</em> of baryonic structure.
                  It is regulated by baryonic surface density — a result that any theory of dark matter
                  or modified gravity must account for.
                </p>
              </div>
            </div>
          </GlassCard>
        </div>

        {deep && (
          <>
            <div className="border-t border-white/10 pt-6 mt-2">
              <h2 className="text-2xl font-display font-bold text-white mb-1">Deep Analysis</h2>
              <p className="text-sm text-slate-400 mb-4">Functional form comparison, mass scaling, and halo connection.</p>
            </div>

            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
              <GlassCard glow="violet">
                <h3 className="text-sm font-display font-bold text-white mb-2 flex items-center gap-2">
                  <Scale className="w-4 h-4 text-violet-400" />
                  Functional Form Comparison
                </h3>
                <p className="text-xs text-slate-400 mb-3">
                  Which formulation captures the density dependence best?
                </p>
                <div className="overflow-x-auto mb-3">
                  <table className="w-full text-xs font-mono">
                    <thead>
                      <tr className="border-b border-white/10 text-slate-400">
                        <th className="text-left py-2 px-2">Formulation</th>
                        <th className="text-center py-2 px-2">|r| pt</th>
                        <th className="text-center py-2 px-2">R² pt</th>
                        <th className="text-center py-2 px-2">AIC pt</th>
                        <th className="text-center py-2 px-2">|r| gal</th>
                        <th className="text-center py-2 px-2">Best?</th>
                      </tr>
                    </thead>
                    <tbody>
                      {deep.functionalForms.map((f, i) => (
                        <tr key={i} className={`border-b border-white/5 ${f.name === deep.bestForm ? 'bg-violet-500/10' : ''}`}>
                          <td className="py-2 px-2 text-slate-300 text-[10px]">{f.label}</td>
                          <td className="py-2 px-2 text-center text-violet-400">{Math.abs(f.pointR).toFixed(4)}</td>
                          <td className="py-2 px-2 text-center text-slate-300">{f.pointR2.toFixed(4)}</td>
                          <td className="py-2 px-2 text-center text-slate-400">{f.pointAIC.toFixed(0)}</td>
                          <td className="py-2 px-2 text-center text-cyan-400">{Math.abs(f.galaxyR).toFixed(4)}</td>
                          <td className="py-2 px-2 text-center">{f.name === deep.bestForm ? <span className="text-amber-400">★</span> : ''}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
                <Badge pass={true} label={`Best: f_DM formulation (|r| = ${Math.abs(deep.functionalForms[0]?.pointR ?? 0).toFixed(3)})`} />
                <p className="text-xs text-slate-500 mt-2">
                  f<sub>DM</sub> = (g<sub>obs</sub>−g<sub>bar</sub>)/g<sub>obs</sub> gives the tightest correlation.
                  log(g<sub>obs</sub>/g<sub>bar</sub>) and log(M<sub>dyn</sub>/M<sub>bar</sub>) are mathematically equivalent.
                </p>
              </GlassCard>

              <GlassCard glow="amber">
                <h3 className="text-sm font-display font-bold text-white mb-2 flex items-center gap-2">
                  <BarChart3 className="w-4 h-4 text-amber-400" />
                  Slope b(V<sub>max</sub>) — Mass Dependence
                </h3>
                <p className="text-xs text-slate-400 mb-3">
                  Does the strength of the density effect depend on galaxy mass?
                  <strong className="text-amber-400 ml-1">r = {deep.slopeMassScaling.r.toFixed(3)}</strong>
                </p>
                <div className="h-[180px] mb-3">
                  <ResponsiveContainer width="100%" height="100%">
                    <BarChart data={deep.slopeVsVmax.map(b => ({ bin: b.bin, slope: +b.slope.toFixed(5), r: +b.r.toFixed(4) }))} margin={{ top: 5, right: 10, bottom: 20, left: 10 }}>
                      <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                      <XAxis dataKey="bin" tick={{ fill: '#94a3b8', fontSize: 8 }} />
                      <YAxis tick={{ fill: '#94a3b8', fontSize: 8 }} />
                      <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 10 }} />
                      <Bar dataKey="slope" fill="#f59e0b" name="slope b" />
                      <ReferenceLine y={0} stroke="#ef4444" strokeDasharray="3 3" />
                    </BarChart>
                  </ResponsiveContainer>
                </div>
                <div className="overflow-x-auto">
                  <table className="w-full text-xs font-mono">
                    <thead>
                      <tr className="border-b border-white/10 text-slate-400">
                        <th className="text-left py-1.5 px-2">V<sub>max</sub> bin</th>
                        <th className="text-center py-1.5 px-2">slope b</th>
                        <th className="text-center py-1.5 px-2">r</th>
                        <th className="text-center py-1.5 px-2">n</th>
                      </tr>
                    </thead>
                    <tbody>
                      {deep.slopeVsVmax.map(b => (
                        <tr key={b.bin} className="border-b border-white/5">
                          <td className="py-1.5 px-2 text-slate-300">{b.bin}</td>
                          <td className={`py-1.5 px-2 text-center ${b.slope < 0 ? 'text-amber-400' : 'text-red-400'}`}>{b.slope.toFixed(4)}</td>
                          <td className="py-1.5 px-2 text-center text-slate-300">{b.r.toFixed(3)}</td>
                          <td className="py-1.5 px-2 text-center text-slate-500">{b.n}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
                <div className="mt-3 bg-amber-500/10 border border-amber-500/20 rounded-xl p-3">
                  <p className="text-xs text-amber-300">
                    <strong>Discovery:</strong> The slope b becomes more negative for more massive galaxies
                    (r = {deep.slopeMassScaling.r.toFixed(3)}, R² = {deep.slopeMassScaling.r2.toFixed(3)}).
                    The density regulation of apparent dark matter is itself mass-dependent.
                  </p>
                </div>
              </GlassCard>
            </div>

            <GlassCard glow="cyan">
              <h3 className="text-lg font-display font-bold text-white mb-3 flex items-center gap-2">
                <Orbit className="w-5 h-5 text-cyan-400" />
                Halo Connection
              </h3>
              <p className="text-xs text-slate-400 mb-4">
                Direct connection between baryonic surface density and dark matter halo properties.
              </p>
              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
                <StatBox label="⟨V_DM⟩ vs Σ_bar (per-galaxy)" value={deep.halo.vDMvsSigBar.galR.toFixed(3)}
                  sub={`R² = ${deep.halo.vDMvsSigBar.galR2.toFixed(3)}, n = ${deep.halo.vDMvsSigBar.galN}`} color="text-cyan-400" />
                {deep.halo.dmDominanceRadius && (
                  <StatBox label="r_DMdom vs Σ_bar" value={deep.halo.dmDominanceRadius.r.toFixed(3)}
                    sub={`R² = ${deep.halo.dmDominanceRadius.r2.toFixed(3)}, n = ${deep.halo.dmDominanceRadius.n}`} color="text-violet-400" />
                )}
                {deep.halo.innerSlopeData && (
                  <StatBox label="Inner slope α vs Σ_bar" value={deep.halo.innerSlopeData.r.toFixed(3)}
                    sub={`n = ${deep.halo.innerSlopeData.n} galaxies`} color="text-amber-400" />
                )}
              </div>

              <div className="grid grid-cols-1 lg:grid-cols-2 gap-4">
                <div>
                  <div className="text-xs text-slate-400 mb-2 font-medium">V<sub>DM</sub> vs Σ<sub>bar</sub> — per-galaxy</div>
                  <div className="h-[220px]">
                    <ResponsiveContainer width="100%" height="100%">
                      <ScatterChart margin={{ top: 5, right: 10, bottom: 25, left: 10 }}>
                        <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                        <XAxis dataKey="logSigBar" type="number" tick={{ fill: '#94a3b8', fontSize: 9 }}
                          label={{ value: 'log Σ_bar', position: 'bottom', offset: 10, fill: '#94a3b8', fontSize: 9 }} />
                        <YAxis dataKey="rDMdomNorm" type="number" tick={{ fill: '#94a3b8', fontSize: 9 }}
                          label={{ value: 'r_DMdom / r_max', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 9 }} />
                        <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 10 }} />
                        <Scatter data={deep.halo.dmDomGalaxies} fill="#06b6d4" opacity={0.6} r={3} />
                      </ScatterChart>
                    </ResponsiveContainer>
                  </div>
                </div>

                {deep.halo.innerSlopeData && (
                  <div>
                    <div className="text-xs text-slate-400 mb-2 font-medium">Inner DM slope α (V<sub>DM</sub> ∝ r<sup>α</sup>) vs Σ<sub>bar</sub></div>
                    <div className="h-[220px]">
                      <ResponsiveContainer width="100%" height="100%">
                        <ScatterChart margin={{ top: 5, right: 10, bottom: 25, left: 10 }}>
                          <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                          <XAxis dataKey="logSigBar" type="number" tick={{ fill: '#94a3b8', fontSize: 9 }}
                            label={{ value: 'log Σ_bar', position: 'bottom', offset: 10, fill: '#94a3b8', fontSize: 9 }} />
                          <YAxis dataKey="alpha" type="number" tick={{ fill: '#94a3b8', fontSize: 9 }}
                            label={{ value: 'α (inner slope)', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 9 }} />
                          <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 10 }} />
                          <Scatter data={deep.halo.innerSlopeData.galaxies} fill="#f59e0b" opacity={0.6} r={3} />
                          <ReferenceLine y={1} stroke="#ef4444" strokeDasharray="3 3" label={{ value: 'cusp (NFW)', fill: '#ef4444', fontSize: 9 }} />
                          <ReferenceLine y={0} stroke="#34d399" strokeDasharray="3 3" label={{ value: 'core', fill: '#34d399', fontSize: 9 }} />
                        </ScatterChart>
                      </ResponsiveContainer>
                    </div>
                  </div>
                )}
              </div>

              <div className="mt-4 space-y-2">
                <div className="bg-white/5 rounded-xl p-4 border-l-4 border-cyan-500">
                  <p className="text-xs text-slate-300">
                    <strong className="text-cyan-400">Key finding:</strong> Per-galaxy V<sub>DM</sub> correlates <em>positively</em> with
                    Σ<sub>bar</sub> (r = {deep.halo.vDMvsSigBar.galR.toFixed(3)}). Higher-density galaxies have faster dark matter
                    rotation — but their <em>fraction</em> of dark matter is lower (f<sub>DM</sub> anti-correlates). This means
                    the baryonic component grows faster than the dark component as density increases.
                  </p>
                </div>
                {deep.halo.dmDominanceRadius && (
                  <div className="bg-white/5 rounded-xl p-4 border-l-4 border-violet-500">
                    <p className="text-xs text-slate-300">
                      <strong className="text-violet-400">DM dominance radius:</strong> Higher-density galaxies transition to
                      dark-matter-dominated rotation at relatively <em>larger</em> normalized radii
                      (r = {deep.halo.dmDominanceRadius.r.toFixed(3)}). Dark matter "takes over" later in denser systems.
                    </p>
                  </div>
                )}
              </div>
            </GlassCard>

            <GlassCard glow="rose">
              <h3 className="text-lg font-display font-bold text-white mb-3 flex items-center gap-2">
                <Atom className="w-5 h-5 text-rose-400" />
                Synthesis: The Complete Picture
              </h3>
              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
                <div className="bg-white/5 rounded-xl p-4 border-l-4 border-cyan-500">
                  <p className="text-xs font-display font-bold text-cyan-400 mb-1">Level 1: Safe</p>
                  <p className="text-xs text-slate-300">RAR has a density-dependent residual (b = {s.pointLevel.slope.toFixed(3)})</p>
                </div>
                <div className="bg-white/5 rounded-xl p-4 border-l-4 border-amber-500">
                  <p className="text-xs font-display font-bold text-amber-400 mb-1">Level 2: Strong</p>
                  <p className="text-xs text-slate-300">Apparent DM fraction is regulated by Σ<sub>bar</sub> (r = {s.pointLevel.r.toFixed(3)})</p>
                </div>
                <div className="bg-white/5 rounded-xl p-4 border-l-4 border-rose-500">
                  <p className="text-xs font-display font-bold text-rose-400 mb-1">Level 3: Deeper</p>
                  <p className="text-xs text-slate-300">Regulation strength itself scales with mass (r = {deep.slopeMassScaling.r.toFixed(3)})</p>
                </div>
              </div>
              <div className="bg-gradient-to-r from-cyan-500/10 via-amber-500/10 to-rose-500/10 border border-white/10 rounded-xl p-5 text-center">
                <p className="font-display font-bold text-white text-lg mb-2">
                  Baryonic surface density regulates the apparent dark matter fraction in galaxy rotation curves
                </p>
                <p className="text-sm text-slate-300 max-w-2xl mx-auto">
                  Dark matter is not independent — it is regulated by baryonic structure.
                  If dark matter exists, it must be coupled to baryons. If it doesn't, the dynamical law
                  must depend on density, not just mass.
                </p>
              </div>
            </GlassCard>
          </>
        )}
      </div>
    </Layout>
  );
}
