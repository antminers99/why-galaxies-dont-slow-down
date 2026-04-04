import React, { useState, useEffect } from 'react';
import { GlassCard } from '@/components/ui/glass-card';
import {
  ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip,
  ResponsiveContainer, BarChart, Bar, Cell, ReferenceLine,
  Label
} from 'recharts';

type FigureData = {
  nGalaxies: number;
  nPoints: number;
  a0Global: number;
  axisRange: { xMin: number; xMax: number; yMin: number; yMax: number };
  rarPoints: { x: number; y: number }[];
  rarCurve: { x: number; y: number }[];
  histogram: {
    bins: { binCenter: number; count: number }[];
    meanLogA0: number;
    sdLogA0: number;
  };
  varianceDecomp: {
    total: { variance: number; rms: number; pct: number };
    within: { variance: number; rms: number; pct: number };
    between: { variance: number; rms: number; pct: number };
    sdLogA0: number;
    tauDL: number;
  };
};

function DecompBar({ label, value, maxVal, color, rms }: {
  label: string; value: number; maxVal: number; color: string; rms: string;
}) {
  const pct = Math.max(1, (value / maxVal) * 100);
  return (
    <div className="mb-3">
      <div className="flex justify-between text-xs mb-1">
        <span className="text-slate-300">{label}</span>
        <span className="font-mono text-slate-400">{rms} dex</span>
      </div>
      <div className="w-full h-5 bg-slate-800/80 rounded-full overflow-hidden">
        <div
          className="h-full rounded-full transition-all duration-700"
          style={{ width: pct + '%', backgroundColor: color }}
        />
      </div>
    </div>
  );
}

export default function Figure1Panel() {
  const [data, setData] = useState<FigureData | null>(null);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'figure1-data.json')
      .then(r => r.json())
      .then(d => setData(d))
      .catch(() => {});
  }, []);

  if (!data) {
    return (
      <GlassCard glow="cyan">
        <div className="text-center text-slate-400 py-8">Loading Figure 1 data...</div>
      </GlassCard>
    );
  }

  const { rarPoints, rarCurve, histogram, varianceDecomp } = data;
  const decomp = varianceDecomp;

  return (
    <GlassCard glow="cyan">
      <div className="flex items-center gap-2 mb-4">
        <h2 className="text-lg font-bold text-white">
          Figure 1: RAR Scatter Decomposition
        </h2>
        <span className="text-xs text-slate-500 font-mono">
          {data.nGalaxies} GOLD+i45 galaxies (benchmark) | {data.nPoints} points
        </span>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-4">

        <div className="bg-slate-900/60 border border-white/5 rounded-xl p-3">
          <h3 className="text-xs font-semibold text-cyan-400 mb-2 uppercase tracking-wider">
            A. RAR Point Cloud
          </h3>
          <p className="text-[10px] text-slate-500 mb-2">
            Within-galaxy RMS = {decomp.within.rms} dex (cf. McGaugh 0.13)
          </p>
          <ResponsiveContainer width="100%" height={280}>
            <ScatterChart margin={{ top: 5, right: 10, left: 0, bottom: 20 }}>
              <CartesianGrid strokeDasharray="3 3" stroke="#334155" opacity={0.3} />
              <XAxis
                dataKey="x"
                type="number"
                domain={[data.axisRange.xMin, data.axisRange.xMax]}
                name="log g_bar"
                tick={{ fill: '#94a3b8', fontSize: 10 }}
                tickCount={5}
              >
                <Label value="log g_bar [(km/s)^2/kpc]" position="bottom" offset={5} style={{ fill: '#64748b', fontSize: 10 }} />
              </XAxis>
              <YAxis
                dataKey="y"
                type="number"
                domain={[data.axisRange.yMin, data.axisRange.yMax]}
                name="log g_obs"
                tick={{ fill: '#94a3b8', fontSize: 10 }}
                tickCount={5}
              >
                <Label value="log g_obs" angle={-90} position="insideLeft" offset={10} style={{ fill: '#64748b', fontSize: 10 }} />
              </YAxis>
              <Tooltip
                contentStyle={{ backgroundColor: '#1e293b', border: '1px solid #334155', borderRadius: 8, fontSize: 11 }}
              />
              <Scatter
                data={rarPoints}
                fill="#22d3ee"
                fillOpacity={0.2}
                r={1.5}
                isAnimationActive={false}
              />
              <Scatter
                data={rarCurve}
                fill="none"
                stroke="#f59e0b"
                strokeWidth={2}
                line={{ stroke: '#f59e0b', strokeWidth: 2 }}
                r={0}
                isAnimationActive={false}
                name="McGaugh RAR"
                legendType="none"
              />
              <Scatter
                data={[{ x: data.axisRange.xMin, y: data.axisRange.yMin }, { x: data.axisRange.xMax, y: data.axisRange.yMax }]}
                fill="none"
                stroke="#64748b"
                line={{ stroke: '#64748b', strokeWidth: 1, strokeDasharray: '6 4' }}
                r={0}
                isAnimationActive={false}
                name="1:1"
                legendType="none"
              />
            </ScatterChart>
          </ResponsiveContainer>
          <div className="flex gap-4 justify-center text-[10px] mt-1">
            <span className="text-cyan-400">&#9679; Data points</span>
            <span className="text-amber-400">&#9473; McGaugh RAR</span>
            <span className="text-slate-400">- - 1:1 (Newtonian)</span>
          </div>
        </div>

        <div className="bg-slate-900/60 border border-white/5 rounded-xl p-3">
          <h3 className="text-xs font-semibold text-purple-400 mb-2 uppercase tracking-wider">
            B. Per-Galaxy log(a&#8320;) Distribution
          </h3>
          <p className="text-[10px] text-slate-500 mb-2">
            SD = {histogram.sdLogA0} dex | tau(DL) = {decomp.tauDL} dex
          </p>
          <ResponsiveContainer width="100%" height={280}>
            <BarChart
              data={histogram.bins.filter(b => b.count > 0 || (b.binCenter >= 2.8 && b.binCenter <= 4.3))}
              margin={{ top: 5, right: 10, left: 0, bottom: 20 }}
            >
              <CartesianGrid strokeDasharray="3 3" stroke="#334155" opacity={0.3} />
              <XAxis
                dataKey="binCenter"
                type="number"
                domain={[2.8, 4.3]}
                tick={{ fill: '#94a3b8', fontSize: 10 }}
                tickCount={6}
              >
                <Label value="log a_0 [(km/s)^2/kpc]" position="bottom" offset={5} style={{ fill: '#64748b', fontSize: 10 }} />
              </XAxis>
              <YAxis
                tick={{ fill: '#94a3b8', fontSize: 10 }}
              >
                <Label value="Count" angle={-90} position="insideLeft" offset={10} style={{ fill: '#64748b', fontSize: 10 }} />
              </YAxis>
              <Tooltip
                contentStyle={{ backgroundColor: '#1e293b', border: '1px solid #334155', borderRadius: 8, fontSize: 11 }}
              />
              <Bar dataKey="count" isAnimationActive={false}>
                {histogram.bins.filter(b => b.count > 0 || (b.binCenter >= 2.8 && b.binCenter <= 4.3)).map((entry, i) => (
                  <Cell key={i} fill={
                    Math.abs(entry.binCenter - histogram.meanLogA0) < 0.15 ? '#a78bfa' : '#6366f1'
                  } fillOpacity={0.7} />
                ))}
              </Bar>
              <ReferenceLine
                x={histogram.meanLogA0}
                stroke="#f59e0b"
                strokeWidth={2}
                strokeDasharray="6 3"
              >
                <Label value="mean" position="top" style={{ fill: '#f59e0b', fontSize: 10 }} />
              </ReferenceLine>
              <ReferenceLine
                x={Math.log10(3703)}
                stroke="#22d3ee"
                strokeWidth={1.5}
                strokeDasharray="3 3"
              >
                <Label value="McGaugh" position="top" style={{ fill: '#22d3ee', fontSize: 9 }} offset={12} />
              </ReferenceLine>
            </BarChart>
          </ResponsiveContainer>
        </div>

        <div className="bg-slate-900/60 border border-white/5 rounded-xl p-3">
          <h3 className="text-xs font-semibold text-amber-400 mb-2 uppercase tracking-wider">
            C. Variance Decomposition
          </h3>
          <p className="text-[10px] text-slate-500 mb-2">
            ANOVA of point residuals (global a&#8320;)
          </p>

          <div className="mt-4 px-2">
            <DecompBar
              label={"Total (global a\u2080)"}
              value={decomp.total.variance}
              maxVal={decomp.total.variance}
              color="#22d3ee"
              rms={decomp.total.rms}
            />
            <DecompBar
              label="Within-galaxy"
              value={decomp.within.variance}
              maxVal={decomp.total.variance}
              color="#34d399"
              rms={decomp.within.rms}
            />
            <DecompBar
              label="Between-galaxy"
              value={decomp.between.variance}
              maxVal={decomp.total.variance}
              color="#f59e0b"
              rms={decomp.between.rms}
            />
          </div>

          <div className="mt-4 space-y-2 text-xs">
            <div className="flex justify-between p-2 bg-slate-800/60 rounded-lg">
              <span className="text-slate-400">Within-galaxy fraction</span>
              <span className="font-mono text-emerald-400">{decomp.within.pct}%</span>
            </div>
            <div className="flex justify-between p-2 bg-slate-800/60 rounded-lg">
              <span className="text-slate-400">Between-galaxy fraction</span>
              <span className="font-mono text-amber-400">{decomp.between.pct}%</span>
            </div>
          </div>

          <div className="mt-4 p-3 bg-slate-800/40 border border-white/5 rounded-xl">
            <h4 className="text-[10px] font-semibold text-slate-400 uppercase mb-2">
              Four Scatter Measures
            </h4>
            <div className="space-y-1.5 text-[11px]">
              <div className="flex justify-between">
                <span className="text-cyan-400">1. Total point RMS</span>
                <span className="font-mono text-slate-300">{decomp.total.rms} dex</span>
              </div>
              <div className="flex justify-between">
                <span className="text-emerald-400">2. Within-galaxy RMS</span>
                <span className="font-mono text-slate-300">{decomp.within.rms} dex</span>
              </div>
              <div className="flex justify-between">
                <span className="text-amber-400">3. Between-galaxy RMS</span>
                <span className="font-mono text-slate-300">{decomp.between.rms} dex</span>
              </div>
              <div className="flex justify-between">
                <span className="text-purple-400">4. SD per-galaxy log(a&#8320;)</span>
                <span className="font-mono text-slate-300">{decomp.sdLogA0} dex</span>
              </div>
              <div className="flex justify-between border-t border-white/10 pt-1.5 mt-1.5">
                <span className="text-red-400">tau (DL hierarchical)</span>
                <span className="font-mono text-slate-300">{decomp.tauDL} dex</span>
              </div>
            </div>
          </div>

          <div className="mt-3 p-2 bg-emerald-950/30 border border-emerald-500/20 rounded-lg">
            <p className="text-[10px] text-emerald-300 leading-relaxed">
              Literature reports #1 or #2 (point-level).
              Our tau is closest to #4 (galaxy-level).
              No tension: different quantities.
            </p>
          </div>
        </div>
      </div>
    </GlassCard>
  );
}
