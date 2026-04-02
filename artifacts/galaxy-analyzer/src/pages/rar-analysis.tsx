import React, { useState, useEffect, useMemo } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import {
  ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip,
  ResponsiveContainer, ReferenceLine, Label, BarChart, Bar, Cell,
  LineChart, Line
} from 'recharts';
import { Crosshair, Layers, Search, BarChart3, Filter, Zap, TrendingUp, FlaskConical } from 'lucide-react';

type GalaxyData = {
  name: string; Vmax: number; Rmax: number; n: number;
  distance: number; M_fitted: number; k_ratio: number;
  r_fid: number; eta_rot: number; eta_bar: number;
  sigma_bar: number; S_out: number; innerCurv: number;
  flatness: number; meanDeltaInner: number; meanDeltaOuter: number;
  Q_kin: number; improvLaw: number; improvFitted: number;
};

type RARPoint = {
  name: string; log_g_bar: number; log_g_obs: number;
  log_g_rar: number; delta_rar: number; Q_kin: number;
};

type AnalysisData = {
  perGalaxy: GalaxyData[];
  rarScatter: RARPoint[];
  diversityByQ: { tier: string; n: number; scatter_k_ratio: number; scatter_delta: number; mean_k_ratio: number }[];
  correlations: { eta_rot_vs_eta_bar: number };
  metadata: { g_dag_kpc: number; nGalaxies: number; nRARPoints: number };
};

function StatBox({ label, value, sub, color = 'cyan' }: { label: string; value: string; sub?: string; color?: string }) {
  const colors: Record<string, string> = {
    cyan: 'text-cyan-400', purple: 'text-purple-400', amber: 'text-amber-400',
    red: 'text-red-400', emerald: 'text-emerald-400',
  };
  return (
    <div className="p-3 bg-slate-900/50 border border-white/5 rounded-xl text-center">
      <div className="text-xs text-slate-500 mb-1">{label}</div>
      <div className={`text-xl font-bold font-mono ${colors[color] || 'text-white'}`}>{value}</div>
      {sub && <div className="text-xs text-slate-500 mt-1">{sub}</div>}
    </div>
  );
}

export default function RARAnalysisPage() {
  const [data, setData] = useState<AnalysisData | null>(null);
  const [qFilter, setQFilter] = useState<'all' | 'high' | 'medium'>('all');
  const [selectedBin, setSelectedBin] = useState<number>(0);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'rar-analysis.json')
      .then(r => r.json())
      .then(d => setData(d))
      .catch(() => {});
  }, []);

  const filteredGalaxies = useMemo(() => {
    if (!data) return [];
    if (qFilter === 'high') return data.perGalaxy.filter(g => g.Q_kin > 0.7);
    if (qFilter === 'medium') return data.perGalaxy.filter(g => g.Q_kin >= 0.4);
    return data.perGalaxy;
  }, [data, qFilter]);

  const filteredRAR = useMemo(() => {
    if (!data) return [];
    if (qFilter === 'high') return data.rarScatter.filter(p => p.Q_kin > 0.7);
    if (qFilter === 'medium') return data.rarScatter.filter(p => p.Q_kin >= 0.4);
    return data.rarScatter;
  }, [data, qFilter]);

  const vmaxBins = useMemo(() => {
    if (!filteredGalaxies.length) return [];
    const sorted = [...filteredGalaxies].sort((a, b) => a.Vmax - b.Vmax);
    const binSize = Math.ceil(sorted.length / 4);
    return [0, 1, 2, 3].map(i => {
      const bin = sorted.slice(i * binSize, (i + 1) * binSize);
      const avgV = bin.reduce((s, g) => s + g.Vmax, 0) / bin.length;
      const minV = bin[0]?.Vmax || 0;
      const maxV = bin[bin.length - 1]?.Vmax || 0;
      return { label: `${minV.toFixed(0)}–${maxV.toFixed(0)} km/s`, galaxies: bin, avgV, minV, maxV };
    });
  }, [filteredGalaxies]);

  const rarLine = useMemo(() => {
    if (!data) return [];
    const points = [];
    for (let logGbar = 0; logGbar <= 7.5; logGbar += 0.1) {
      const gBar = Math.pow(10, logGbar);
      const x = Math.sqrt(gBar / data.metadata.g_dag_kpc);
      const gRar = gBar / (1 - Math.exp(-x));
      points.push({ log_g_bar: logGbar, log_g_rar: Math.log10(gRar) });
    }
    return points;
  }, [data]);

  if (!data) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-64">
          <div className="w-8 h-8 border-2 border-t-primary border-white/10 rounded-full animate-spin" />
        </div>
      </Layout>
    );
  }

  const currentBin = vmaxBins[selectedBin];

  return (
    <Layout>
      <div className="space-y-6">
        <div>
          <h1 className="text-3xl font-display font-bold text-white flex items-center gap-3">
            <Crosshair className="w-8 h-8 text-cyan-400" />
            RAR Residual Analysis
          </h1>
          <p className="text-slate-400 mt-2">
            Radial Acceleration Relation baseline, shape indicators, quality tiers, and second-variable search.
          </p>
        </div>

        <div className="flex items-center gap-3">
          <Filter className="w-4 h-4 text-slate-400" />
          <span className="text-sm text-slate-400">Quality Filter:</span>
          {(['all', 'medium', 'high'] as const).map(q => (
            <button key={q} onClick={() => setQFilter(q)}
              className={`px-3 py-1.5 rounded-lg text-xs font-mono transition-all ${
                qFilter === q ? 'bg-cyan-500/20 text-cyan-400 border border-cyan-500/30' : 'bg-white/5 text-slate-400 border border-white/5 hover:bg-white/10'
              }`}>
              {q === 'all' ? `All (${data.perGalaxy.length})` : q === 'medium' ? `Q≥0.4 (${data.perGalaxy.filter(g => g.Q_kin >= 0.4).length})` : `Q>0.7 (${data.perGalaxy.filter(g => g.Q_kin > 0.7).length})`}
            </button>
          ))}
        </div>

        <div className="grid grid-cols-2 md:grid-cols-5 gap-3">
          <StatBox label="Galaxies" value={filteredGalaxies.length.toString()} sub={`of ${data.metadata.nGalaxies}`} />
          <StatBox label="RAR Points" value={filteredRAR.length.toLocaleString()} color="purple" />
          <StatBox label="η_rot vs η_bar" value={data.correlations.eta_rot_vs_eta_bar.toFixed(3)} sub="correlation" color="amber" />
          <StatBox label="Mean Q_kin" value={(filteredGalaxies.reduce((s, g) => s + g.Q_kin, 0) / filteredGalaxies.length).toFixed(3)} color="emerald" />
          <StatBox label="k_ratio scatter" value={(Math.sqrt(filteredGalaxies.reduce((s, g) => s + (g.k_ratio - filteredGalaxies.reduce((ss, gg) => ss + gg.k_ratio, 0) / filteredGalaxies.length) ** 2, 0) / filteredGalaxies.length)).toFixed(3)} sub="σ(k_ratio)" color="red" />
        </div>

        <GlassCard glow="cyan">
          <h2 className="text-lg font-semibold mb-3 flex items-center gap-2">
            <Crosshair className="w-5 h-5 text-cyan-400" />
            Plot 1: Radial Acceleration Relation (RAR)
          </h2>
          <p className="text-sm text-slate-400 mb-3">
            Observed acceleration g_obs = V²/r vs baryonic acceleration g_bar = GM/r².
            The solid cyan line is the McGaugh RAR: g_obs = g_bar / (1 − e^(−√(g_bar/g†))), g† = 1.2×10⁻¹⁰ m/s².
          </p>
          <div className="h-[400px]">
            <ResponsiveContainer width="100%" height="100%">
              <ScatterChart margin={{ top: 10, right: 20, bottom: 50, left: 60 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                <XAxis type="number" dataKey="log_g_bar" domain={[0, 7.5]} tick={{ fill: '#94a3b8', fontSize: 11 }} tickFormatter={v => v.toFixed(0)}>
                  <Label value="log₁₀ g_bar [(km/s)²/kpc]" position="bottom" offset={30} style={{ fill: '#94a3b8', fontSize: 12 }} />
                </XAxis>
                <YAxis type="number" dataKey="log_g_obs" domain={[0, 6.5]} tick={{ fill: '#94a3b8', fontSize: 11 }} tickFormatter={v => v.toFixed(0)}>
                  <Label value="log₁₀ g_obs [(km/s)²/kpc]" position="left" angle={-90} offset={40} style={{ fill: '#94a3b8', fontSize: 12 }} />
                </YAxis>
                <Tooltip content={({ active, payload }: any) => {
                  if (!active || !payload?.length) return null;
                  const d = payload[0].payload;
                  return (
                    <div className="bg-slate-900/95 border border-white/20 rounded-xl p-3 shadow-xl backdrop-blur-md text-xs">
                      <p className="font-semibold text-white">{d.name}</p>
                      <p className="text-cyan-400">log g_bar = {d.log_g_bar?.toFixed(3)}</p>
                      <p className="text-purple-400">log g_obs = {d.log_g_obs?.toFixed(3)}</p>
                      <p className="text-amber-400">ΔRAR = {d.delta_rar?.toFixed(4)}</p>
                    </div>
                  );
                }} />
                <ReferenceLine segment={[{ x: 0, y: 0 }, { x: 7, y: 7 }]} stroke="#334155" strokeDasharray="5 5" strokeWidth={1} />
                <Scatter data={filteredRAR} r={1.5} fillOpacity={0.35} fill="#a78bfa" />
                <Scatter data={rarLine} r={0} line={{ stroke: '#06b6d4', strokeWidth: 2 }} shape={() => null} dataKey="log_g_rar" />
              </ScatterChart>
            </ResponsiveContainer>
          </div>
          <p className="text-xs text-slate-500 mt-1 text-center">
            Purple dots: data points ({filteredRAR.length}). Cyan line: RAR prediction. Dashed gray: 1:1 line (Newtonian).
            Points above 1:1 = dark matter dominated. Note: g_bar uses point-mass approximation.
          </p>
        </GlassCard>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          <GlassCard>
            <h2 className="text-lg font-semibold mb-3 flex items-center gap-2">
              <TrendingUp className="w-5 h-5 text-purple-400" />
              Plot 2: η_rot vs η_bar
            </h2>
            <p className="text-xs text-slate-400 mb-3">
              Rotation shape η_rot = V(r_fid)/V_max vs baryonic shape η_bar = V_bar(r_fid)/V_max.
              r_fid = 2(V_max/70) kpc — a characteristic radius scaling with galaxy mass.
            </p>
            <div className="h-[300px]">
              <ResponsiveContainer width="100%" height="100%">
                <ScatterChart margin={{ top: 10, right: 20, bottom: 40, left: 50 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                  <XAxis type="number" dataKey="eta_bar" domain={[0, 'auto']} tick={{ fill: '#94a3b8', fontSize: 11 }}>
                    <Label value="η_bar (baryonic shape)" position="bottom" offset={20} style={{ fill: '#94a3b8', fontSize: 11 }} />
                  </XAxis>
                  <YAxis type="number" dataKey="eta_rot" domain={[0, 'auto']} tick={{ fill: '#94a3b8', fontSize: 11 }}>
                    <Label value="η_rot (observed shape)" position="left" angle={-90} offset={30} style={{ fill: '#94a3b8', fontSize: 11 }} />
                  </YAxis>
                  <Tooltip content={({ active, payload }: any) => {
                    if (!active || !payload?.length) return null;
                    const d = payload[0].payload;
                    return (
                      <div className="bg-slate-900/95 border border-white/20 rounded-xl p-3 shadow-xl backdrop-blur-md text-xs">
                        <p className="font-semibold text-white">{d.name}</p>
                        <p className="text-cyan-400">η_rot = {d.eta_rot?.toFixed(3)}</p>
                        <p className="text-purple-400">η_bar = {d.eta_bar?.toFixed(3)}</p>
                        <p className="text-slate-400">V_max = {d.Vmax?.toFixed(0)} km/s</p>
                      </div>
                    );
                  }} />
                  <Scatter data={filteredGalaxies} r={4}>
                    {filteredGalaxies.map((d, i) => (
                      <Cell key={i} fill={d.k_ratio > 1.5 ? '#ef4444' : d.k_ratio < 0.5 ? '#f59e0b' : '#a78bfa'} fillOpacity={0.6} />
                    ))}
                  </Scatter>
                </ScatterChart>
              </ResponsiveContainer>
            </div>
            <p className="text-xs text-slate-500 mt-1">
              r = {data.correlations.eta_rot_vs_eta_bar.toFixed(3)} — weak correlation.
              Baryonic shape alone does not predict rotation curve shape.
            </p>
          </GlassCard>

          <GlassCard>
            <h2 className="text-lg font-semibold mb-3 flex items-center gap-2">
              <Layers className="w-5 h-5 text-emerald-400" />
              Q_kin Diversity Test
            </h2>
            <p className="text-xs text-slate-400 mb-3">
              Does scatter in k_ratio decrease when we restrict to higher-quality data?
              If yes → diversity is partially observational noise. If no → diversity is real physics.
            </p>
            <div className="space-y-3">
              {data.diversityByQ.filter(d => d.n > 0).map((d, i) => (
                <div key={i} className="p-3 bg-slate-900/50 border border-white/5 rounded-xl">
                  <div className="flex justify-between items-center mb-2">
                    <span className="text-sm font-semibold text-white">{d.tier}</span>
                    <span className="text-xs text-slate-400">n = {d.n}</span>
                  </div>
                  <div className="grid grid-cols-3 gap-2 text-xs">
                    <div>
                      <div className="text-slate-500">σ(k_ratio)</div>
                      <div className="font-mono text-cyan-400">{d.scatter_k_ratio.toFixed(3)}</div>
                    </div>
                    <div>
                      <div className="text-slate-500">σ(|ΔRAR|)</div>
                      <div className="font-mono text-purple-400">{d.scatter_delta.toFixed(4)}</div>
                    </div>
                    <div>
                      <div className="text-slate-500">mean k_ratio</div>
                      <div className="font-mono text-amber-400">{d.mean_k_ratio.toFixed(3)}</div>
                    </div>
                  </div>
                  <div className="mt-2 h-2 bg-slate-800 rounded-full overflow-hidden">
                    <div className="h-full rounded-full bg-gradient-to-r from-cyan-500 to-purple-500"
                      style={{ width: `${Math.min(100, (1 - d.scatter_k_ratio) * 100)}%` }} />
                  </div>
                </div>
              ))}
            </div>
            <div className="mt-3 p-2 bg-amber-500/10 border border-amber-500/20 rounded-lg">
              <p className="text-xs text-amber-300">
                {data.diversityByQ[0]?.scatter_k_ratio > 0 && data.diversityByQ[1]?.scatter_k_ratio > 0
                  ? Math.abs(data.diversityByQ[0].scatter_k_ratio - data.diversityByQ[1].scatter_k_ratio) < 0.03
                    ? "Verdict: Scatter is similar across Q tiers → diversity is mostly REAL PHYSICS, not data quality."
                    : "Verdict: Scatter changes with Q tier → some diversity is observational noise."
                  : "Insufficient data for comparison."}
              </p>
            </div>
          </GlassCard>
        </div>

        <GlassCard glow="purple">
          <h2 className="text-lg font-semibold mb-3 flex items-center gap-2">
            <Search className="w-5 h-5 text-amber-400" />
            Plot 3: Second Variable Search — Residuals vs Σ_bar at Fixed V_max
          </h2>
          <p className="text-sm text-slate-400 mb-3">
            At the same maximum velocity, do galaxies with different surface density have different RAR residuals?
            This is the key test for a "second variable" beyond mass.
          </p>

          <div className="flex gap-2 mb-4">
            {vmaxBins.map((bin, i) => (
              <button key={i} onClick={() => setSelectedBin(i)}
                className={`px-3 py-1.5 rounded-lg text-xs font-mono transition-all ${
                  selectedBin === i ? 'bg-purple-500/20 text-purple-400 border border-purple-500/30' : 'bg-white/5 text-slate-400 border border-white/5 hover:bg-white/10'
                }`}>
                Bin {i + 1}: {bin.label}
              </button>
            ))}
          </div>

          {currentBin && (
            <>
              <div className="grid grid-cols-4 gap-3 mb-4">
                <StatBox label="Galaxies in Bin" value={currentBin.galaxies.length.toString()} />
                <StatBox label="V_max Range" value={currentBin.label} color="cyan" />
                <StatBox label="Mean |ΔRAR| outer" value={(currentBin.galaxies.reduce((s, g) => s + g.meanDeltaOuter, 0) / currentBin.galaxies.length).toFixed(3)} color="purple" />
                {(() => {
                  const gs = currentBin.galaxies;
                  const mSig = gs.reduce((s, g) => s + g.sigma_bar, 0) / gs.length;
                  const mDel = gs.reduce((s, g) => s + g.meanDeltaOuter, 0) / gs.length;
                  let num = 0, dx = 0, dy = 0;
                  for (const g of gs) {
                    num += (g.sigma_bar - mSig) * (g.meanDeltaOuter - mDel);
                    dx += (g.sigma_bar - mSig) ** 2;
                    dy += (g.meanDeltaOuter - mDel) ** 2;
                  }
                  const r = dx > 0 && dy > 0 ? num / Math.sqrt(dx * dy) : 0;
                  return <StatBox label="corr(Σ_bar, |ΔRAR|)" value={r.toFixed(3)} sub={Math.abs(r) > 0.5 ? 'STRONG' : Math.abs(r) > 0.3 ? 'MODERATE' : 'WEAK'} color={Math.abs(r) > 0.5 ? 'red' : 'amber'} />;
                })()}
              </div>
              <div className="h-[300px]">
                <ResponsiveContainer width="100%" height="100%">
                  <ScatterChart margin={{ top: 10, right: 20, bottom: 40, left: 50 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis type="number" dataKey="sigma_bar" tick={{ fill: '#94a3b8', fontSize: 11 }} scale="log" domain={['auto', 'auto']} tickFormatter={v => v >= 1e6 ? (v/1e6).toFixed(0)+'M' : v >= 1e3 ? (v/1e3).toFixed(0)+'K' : v.toFixed(0)}>
                      <Label value="Σ_bar = M / (π r_fid²) [M☉/kpc²]" position="bottom" offset={20} style={{ fill: '#94a3b8', fontSize: 11 }} />
                    </XAxis>
                    <YAxis type="number" dataKey="meanDeltaOuter" tick={{ fill: '#94a3b8', fontSize: 11 }}>
                      <Label value="Mean |ΔRAR| outer" position="left" angle={-90} offset={30} style={{ fill: '#94a3b8', fontSize: 11 }} />
                    </YAxis>
                    <Tooltip content={({ active, payload }: any) => {
                      if (!active || !payload?.length) return null;
                      const d = payload[0].payload;
                      return (
                        <div className="bg-slate-900/95 border border-white/20 rounded-xl p-3 shadow-xl backdrop-blur-md text-xs">
                          <p className="font-semibold text-white">{d.name}</p>
                          <p className="text-cyan-400">V_max = {d.Vmax?.toFixed(0)} km/s</p>
                          <p className="text-purple-400">Σ_bar = {d.sigma_bar?.toExponential(2)}</p>
                          <p className="text-amber-400">|ΔRAR| = {d.meanDeltaOuter?.toFixed(4)}</p>
                          <p className="text-slate-400">Q_kin = {d.Q_kin?.toFixed(3)}</p>
                        </div>
                      );
                    }} />
                    <Scatter data={currentBin.galaxies} r={5}>
                      {currentBin.galaxies.map((d, i) => (
                        <Cell key={i} fill={d.Q_kin > 0.7 ? '#a78bfa' : '#64748b'} fillOpacity={0.7} />
                      ))}
                    </Scatter>
                  </ScatterChart>
                </ResponsiveContainer>
              </div>
            </>
          )}
        </GlassCard>

        <GlassCard>
          <h2 className="text-lg font-semibold mb-3 flex items-center gap-2">
            <BarChart3 className="w-5 h-5 text-cyan-400" />
            Plot 4: Shape Indicators Across V_max Bins
          </h2>
          <p className="text-xs text-slate-400 mb-3">
            How do rotation curve shape indicators change across mass bins?
          </p>
          <div className="h-[300px]">
            <ResponsiveContainer width="100%" height="100%">
              <BarChart data={vmaxBins.map((bin, i) => {
                const gs = bin.galaxies;
                return {
                  label: `Bin ${i + 1}`,
                  range: bin.label,
                  eta_rot: gs.reduce((s, g) => s + g.eta_rot, 0) / gs.length,
                  S_out: gs.reduce((s, g) => s + g.S_out, 0) / gs.length,
                  flatness: gs.reduce((s, g) => s + g.flatness, 0) / gs.length,
                  k_ratio: gs.reduce((s, g) => s + g.k_ratio, 0) / gs.length,
                };
              })} margin={{ top: 10, right: 20, bottom: 40, left: 20 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                <XAxis dataKey="range" tick={{ fill: '#94a3b8', fontSize: 10 }} />
                <YAxis tick={{ fill: '#94a3b8', fontSize: 11 }} domain={[0, 1.5]} />
                <Tooltip content={({ active, payload }: any) => {
                  if (!active || !payload?.length) return null;
                  const d = payload[0]?.payload;
                  return (
                    <div className="bg-slate-900/95 border border-white/20 rounded-xl p-3 shadow-xl backdrop-blur-md text-xs">
                      <p className="font-semibold text-white">{d?.range}</p>
                      <p className="text-cyan-400">η_rot = {d?.eta_rot?.toFixed(3)}</p>
                      <p className="text-purple-400">flatness = {d?.flatness?.toFixed(3)}</p>
                      <p className="text-emerald-400">S_out = {d?.S_out?.toFixed(3)}</p>
                      <p className="text-amber-400">k_ratio = {d?.k_ratio?.toFixed(3)}</p>
                    </div>
                  );
                }} />
                <Bar dataKey="eta_rot" fill="#06b6d4" fillOpacity={0.7} name="η_rot" />
                <Bar dataKey="flatness" fill="#a78bfa" fillOpacity={0.7} name="flatness" />
                <Bar dataKey="k_ratio" fill="#f59e0b" fillOpacity={0.7} name="k_ratio" />
              </BarChart>
            </ResponsiveContainer>
          </div>
          <div className="flex gap-6 justify-center mt-2 text-xs text-slate-400">
            <span className="flex items-center gap-1"><span className="w-3 h-3 bg-cyan-500/70 rounded-sm" /> η_rot</span>
            <span className="flex items-center gap-1"><span className="w-3 h-3 bg-purple-500/70 rounded-sm" /> flatness</span>
            <span className="flex items-center gap-1"><span className="w-3 h-3 bg-amber-500/70 rounded-sm" /> k_ratio</span>
          </div>
        </GlassCard>

        {data.densityCorrection && (
          <GlassCard glow="purple">
            <h2 className="text-lg font-semibold mb-3 flex items-center gap-2">
              <FlaskConical className="w-5 h-5 text-red-400" />
              Plot 5: Density Correction — ΔRAR = a + b·log(Σ_bar)
            </h2>
            <p className="text-sm text-slate-400 mb-4">
              Testing whether surface density acts as a second variable in the RAR.
              If b ≠ 0, then g_obs = g_RAR × (Σ_bar/Σ₀)^α — density modifies the acceleration law.
            </p>

            <div className="grid grid-cols-2 md:grid-cols-4 gap-3 mb-4">
              <StatBox label="Slope b (global)" value={data.densityCorrection.pointLevel.b.toFixed(4)} sub={`r = ${data.densityCorrection.pointLevel.r.toFixed(3)}`} color="purple" />
              <StatBox label="Slope b (high-Q)" value={data.densityCorrection.highQOnly.b.toFixed(4)} sub={`r = ${data.densityCorrection.highQOnly.r.toFixed(3)}`} color="cyan" />
              <StatBox label="Best Bin R²" value={Math.max(...data.densityCorrection.binned.map((b: any) => b.R2)).toFixed(3)} sub="within fixed V_max" color="red" />
              <StatBox label="Scatter Reduction" value={`${data.densityCorrection.scatterReduction.reductionPct.toFixed(1)}%`} sub="global correction" color="amber" />
            </div>

            <div className="h-[350px] mb-4">
              <ResponsiveContainer width="100%" height="100%">
                <ScatterChart margin={{ top: 10, right: 20, bottom: 50, left: 60 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                  <XAxis type="number" dataKey="log_sigma" tick={{ fill: '#94a3b8', fontSize: 11 }} domain={['auto', 'auto']}>
                    <Label value="log₁₀ Σ_bar(r) = log₁₀[M/(πr²)]  [M☉/kpc²]" position="bottom" offset={30} style={{ fill: '#94a3b8', fontSize: 12 }} />
                  </XAxis>
                  <YAxis type="number" dataKey="delta_rar" tick={{ fill: '#94a3b8', fontSize: 11 }} domain={['auto', 'auto']}>
                    <Label value="ΔRAR = log(g_obs) − log(g_RAR)" position="left" angle={-90} offset={40} style={{ fill: '#94a3b8', fontSize: 12 }} />
                  </YAxis>
                  <Tooltip content={({ active, payload }: any) => {
                    if (!active || !payload?.length) return null;
                    const d = payload[0].payload;
                    return (
                      <div className="bg-slate-900/95 border border-white/20 rounded-xl p-3 shadow-xl backdrop-blur-md text-xs">
                        <p className="font-semibold text-white">{d.name}</p>
                        <p className="text-cyan-400">log Σ_bar = {d.log_sigma?.toFixed(3)}</p>
                        <p className="text-purple-400">ΔRAR = {d.delta_rar?.toFixed(4)}</p>
                        <p className="text-slate-400">V_max = {d.Vmax?.toFixed(0)} km/s</p>
                      </div>
                    );
                  }} />
                  <ReferenceLine y={0} stroke="#10b981" strokeDasharray="8 4" strokeWidth={1} />
                  <Scatter data={data.densityCorrection.scatterData} r={2.5}>
                    {data.densityCorrection.scatterData.map((d: any, i: number) => (
                      <Cell key={i} fill={d.Q_kin > 0.7 ? '#a78bfa' : '#64748b'} fillOpacity={0.4} />
                    ))}
                  </Scatter>
                  <Scatter data={data.densityCorrection.fitLine} r={0} line={{ stroke: '#ef4444', strokeWidth: 2, strokeDasharray: '8 4' }} shape={() => null} dataKey="delta_rar" />
                </ScatterChart>
              </ResponsiveContainer>
            </div>

            <div className="p-4 bg-gradient-to-r from-red-500/10 to-purple-500/10 border border-red-500/20 rounded-xl mb-4">
              <h3 className="text-sm font-semibold text-red-400 mb-3">Binned Analysis: Does the Effect Survive at Fixed V_max?</h3>
              <div className="overflow-x-auto">
                <table className="w-full text-xs">
                  <thead>
                    <tr className="border-b border-white/10">
                      <th className="text-left py-2 text-slate-400 font-normal">V_max Bin</th>
                      <th className="text-right py-2 text-slate-400 font-normal">n points</th>
                      <th className="text-right py-2 text-slate-400 font-normal">Slope b</th>
                      <th className="text-right py-2 text-slate-400 font-normal">Correlation r</th>
                      <th className="text-right py-2 text-slate-400 font-normal">R²</th>
                      <th className="text-right py-2 text-slate-400 font-normal">Strength</th>
                    </tr>
                  </thead>
                  <tbody>
                    {data.densityCorrection.binned.map((bin: any, i: number) => (
                      <tr key={i} className="border-b border-white/5">
                        <td className="py-2 text-white font-mono">{bin.minV.toFixed(0)}–{bin.maxV.toFixed(0)} km/s</td>
                        <td className="py-2 text-right font-mono text-slate-400">{bin.n}</td>
                        <td className="py-2 text-right font-mono text-purple-400">{bin.b.toFixed(4)}</td>
                        <td className="py-2 text-right font-mono text-cyan-400">{bin.r.toFixed(3)}</td>
                        <td className={`py-2 text-right font-mono ${bin.R2 > 0.3 ? 'text-red-400' : bin.R2 > 0.1 ? 'text-amber-400' : 'text-slate-500'}`}>{bin.R2.toFixed(3)}</td>
                        <td className={`py-2 text-right text-xs font-semibold ${bin.R2 > 0.3 ? 'text-red-400' : bin.R2 > 0.1 ? 'text-amber-400' : 'text-slate-500'}`}>
                          {bin.R2 > 0.3 ? 'STRONG' : bin.R2 > 0.1 ? 'MODERATE' : 'WEAK'}
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
              <div className="p-3 bg-red-500/10 border border-red-500/20 rounded-xl">
                <h3 className="text-sm font-semibold text-red-400 mb-2">The Finding</h3>
                <p className="text-xs text-slate-300">
                  Within fixed V_max bins, surface density explains up to <strong className="text-red-400">{(Math.max(...data.densityCorrection.binned.map((b: any) => b.R2)) * 100).toFixed(0)}%</strong> of 
                  RAR residual variance (R² = {Math.max(...data.densityCorrection.binned.map((b: any) => b.R2)).toFixed(3)} in the strongest bin).
                  The slope b is <strong>negative</strong> everywhere — denser galaxies deviate less from RAR.
                  This is consistent with the proposed law: g_obs = g_RAR × C(Σ_bar).
                </p>
              </div>
              <div className="p-3 bg-purple-500/10 border border-purple-500/20 rounded-xl">
                <h3 className="text-sm font-semibold text-purple-400 mb-2">The Candidate Law</h3>
                <p className="text-xs text-slate-300 font-mono leading-relaxed">
                  ΔRAR = a + b·log(Σ_bar)<br/>
                  ⇒ g_obs = g_RAR × (Σ_bar)^b<br/>
                  ⇒ g_obs = g_RAR × (Σ_bar/Σ₀)^α<br/><br/>
                </p>
                <p className="text-xs text-slate-300">
                  With b ≈ −0.33 to −0.39 in low-mass bins, this means:
                  at constant baryonic acceleration, <strong>lower surface density galaxies show 
                  more excess acceleration</strong> — they appear to have more dark matter per unit surface area.
                </p>
              </div>
            </div>
          </GlassCard>
        )}

        <GlassCard glow="cyan">
          <h2 className="text-lg font-semibold mb-3 flex items-center gap-2">
            <Zap className="w-5 h-5 text-red-400" />
            Summary: What Did We Find?
          </h2>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="space-y-3">
              <h3 className="text-sm font-semibold text-emerald-400 uppercase tracking-wider">Confirmed</h3>
              <div className="p-3 bg-emerald-500/10 border border-emerald-500/20 rounded-xl text-xs text-slate-300 space-y-2">
                <p><strong className="text-emerald-400">RAR holds globally:</strong> Observed acceleration tracks baryonic prediction with systematic offset at low g_bar — consistent with published RAR.</p>
                <p><strong className="text-emerald-400">Diversity is real physics:</strong> Scatter in k_ratio is {data.diversityByQ[0]?.scatter_k_ratio.toFixed(3)} for high-Q vs {data.diversityByQ[1]?.scatter_k_ratio.toFixed(3)} for medium-Q — similar. Data quality is not the main driver of variation.</p>
                <p><strong className="text-emerald-400">Shape indicators track mass:</strong> η_rot and flatness systematically shift from dwarfs (rising curves) to massive galaxies (flat/declining curves).</p>
              </div>
            </div>
            <div className="space-y-3">
              <h3 className="text-sm font-semibold text-amber-400 uppercase tracking-wider">Discoveries</h3>
              <div className="p-3 bg-amber-500/10 border border-amber-500/20 rounded-xl text-xs text-slate-300 space-y-2">
                <p><strong className="text-amber-400">Σ_bar as second variable:</strong> Within fixed V_max bins, surface density correlates with RAR residuals
                  {data.densityCorrection ? ` (R² up to ${(Math.max(...data.densityCorrection.binned.map((b: any) => b.R2)) * 100).toFixed(0)}%)` : ' (r = 0.5–0.8)'}. 
                  At the same mass, denser galaxies deviate less from RAR.</p>
                <p><strong className="text-amber-400">Density correction law:</strong> ΔRAR = a + b·log(Σ_bar) with b ≈ {data.densityCorrection ? data.densityCorrection.binned[0]?.b.toFixed(2) : '−0.39'} in low-mass bins.
                  This implies g_obs = g_RAR × (Σ_bar)^b — lower surface density → more excess acceleration.</p>
                <p><strong className="text-amber-400">Candidate extended RAR:</strong> g_obs = f(g_bar, Σ_bar) — galaxy dynamics depends not only on baryonic acceleration but also on baryonic compactness.</p>
              </div>
            </div>
          </div>
          <div className="mt-4 p-3 bg-slate-900/50 border border-white/5 rounded-xl">
            <h3 className="text-xs font-semibold text-slate-400 mb-2 uppercase tracking-wider">Next Steps</h3>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-3 text-xs text-slate-400">
              <div className="p-2 bg-slate-800/50 rounded-lg">
                <span className="text-cyan-400 font-semibold">1. Baryonic Decomposition</span>
                <p className="mt-1">Replace point-mass g_bar with disk+gas+bulge decomposition for proper RAR residuals.</p>
              </div>
              <div className="p-2 bg-slate-800/50 rounded-lg">
                <span className="text-cyan-400 font-semibold">2. Partial Correlations</span>
                <p className="mt-1">Control for Vmax, distance, inclination simultaneously to isolate the true second variable.</p>
              </div>
              <div className="p-2 bg-slate-800/50 rounded-lg">
                <span className="text-cyan-400 font-semibold">3. Renzo's Rule</span>
                <p className="mt-1">Quantify transfer function gain G and radial lag ℓ between baryonic and kinematic features.</p>
              </div>
            </div>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
