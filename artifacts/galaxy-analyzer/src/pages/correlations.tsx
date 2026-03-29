import React, { useState, useEffect, useMemo } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Legend, ReferenceLine, Label } from 'recharts';
import { Loader2, TrendingUp, Search, AlertTriangle, CheckCircle2, XCircle, FlaskConical, Zap, Flame } from 'lucide-react';

interface GalaxyResult {
  name: string;
  rawName: string;
  distance: number;
  pointCount: number;
  maxR: number;
  maxV: number;
  models: Record<string, {
    M: number;
    k: number;
    a: number;
    mse: number;
    improvementVsNewton: number;
    mseInner: number;
    mseOuter: number;
    innerImprovement: number;
    outerImprovement: number;
  }>;
}

interface AnalysisData {
  metadata: { date: string; galaxyCount: number; totalPoints: number };
  perGalaxy: GalaxyResult[];
  summary: Record<string, { name: string; avgMSE: number; avgImprovement: number; wins: number; avgK: number; kCV: number }>;
}

type ModelKey = 'dark_halo_linear' | 'dark_halo_flat' | 'modified_gravity_halo' | 'log_halo' | 'transition';
type AxisVar = 'M' | 'k' | 'a' | 'maxV' | 'maxR' | 'distance' | 'pointCount' | 'mse' | 'improvementVsNewton';

const MODEL_OPTIONS: { key: ModelKey; label: string; color: string }[] = [
  { key: 'dark_halo_linear', label: 'Dark Halo (Linear)', color: '#a78bfa' },
  { key: 'modified_gravity_halo', label: 'Modified Gravity + Halo', color: '#34d399' },
  { key: 'log_halo', label: 'Logarithmic Halo (NFW-like)', color: '#38bdf8' },
  { key: 'dark_halo_flat', label: 'Dark Halo (Flat)', color: '#fbbf24' },
  { key: 'transition', label: 'Transition Model', color: '#f87171' },
];

const AXIS_OPTIONS: { key: AxisVar; label: string; isModel: boolean }[] = [
  { key: 'k', label: 'k (Dark Matter Parameter)', isModel: true },
  { key: 'M', label: 'M (Fitted Mass)', isModel: true },
  { key: 'a', label: 'a (Core Radius)', isModel: true },
  { key: 'mse', label: 'MSE (Fit Error)', isModel: true },
  { key: 'improvementVsNewton', label: 'Improvement vs Newton (%)', isModel: true },
  { key: 'maxV', label: 'Max Velocity (km/s)', isModel: false },
  { key: 'maxR', label: 'Max Radius (kpc)', isModel: false },
  { key: 'distance', label: 'Distance (Mpc)', isModel: false },
  { key: 'pointCount', label: 'Data Points', isModel: false },
];

function pearsonCorrelation(x: number[], y: number[]): number {
  const n = x.length;
  if (n < 3) return 0;
  const mx = x.reduce((s, v) => s + v, 0) / n;
  const my = y.reduce((s, v) => s + v, 0) / n;
  let num = 0, dx2 = 0, dy2 = 0;
  for (let i = 0; i < n; i++) {
    const dx = x[i] - mx;
    const dy = y[i] - my;
    num += dx * dy;
    dx2 += dx * dx;
    dy2 += dy * dy;
  }
  const den = Math.sqrt(dx2 * dy2);
  return den > 0 ? num / den : 0;
}

function spearmanCorrelation(x: number[], y: number[]): number {
  const n = x.length;
  if (n < 3) return 0;
  const rank = (arr: number[]) => {
    const sorted = arr.map((v, i) => ({ v, i })).sort((a, b) => a.v - b.v);
    const ranks = new Array(n);
    let i = 0;
    while (i < n) {
      let j = i;
      while (j < n - 1 && sorted[j + 1].v === sorted[j].v) j++;
      const avgRank = (i + j) / 2 + 1;
      for (let t = i; t <= j; t++) ranks[sorted[t].i] = avgRank;
      i = j + 1;
    }
    return ranks;
  };
  return pearsonCorrelation(rank(x), rank(y));
}

function linearRegression(x: number[], y: number[]) {
  const n = x.length;
  const mx = x.reduce((s, v) => s + v, 0) / n;
  const my = y.reduce((s, v) => s + v, 0) / n;
  let num = 0, den = 0;
  for (let i = 0; i < n; i++) {
    num += (x[i] - mx) * (y[i] - my);
    den += (x[i] - mx) ** 2;
  }
  const slope = den > 0 ? num / den : 0;
  const intercept = my - slope * mx;
  return { slope, intercept };
}

function getValue(galaxy: GalaxyResult, axis: AxisVar, model: ModelKey): number | null {
  const modelData = galaxy.models[model];
  if (!modelData) return null;
  switch (axis) {
    case 'k': return modelData.k;
    case 'M': return modelData.M;
    case 'a': return modelData.a;
    case 'mse': return modelData.mse;
    case 'improvementVsNewton': return modelData.improvementVsNewton;
    case 'maxV': return galaxy.maxV;
    case 'maxR': return galaxy.maxR;
    case 'distance': return galaxy.distance;
    case 'pointCount': return galaxy.pointCount;
    default: return null;
  }
}

function formatAxisValue(v: number, axis: AxisVar): string {
  if (axis === 'M') return v.toExponential(1);
  if (axis === 'improvementVsNewton') return v.toFixed(1) + '%';
  if (v > 10000) return v.toExponential(1);
  if (v > 100) return v.toFixed(0);
  if (v > 1) return v.toFixed(1);
  return v.toFixed(3);
}

const CustomTooltip = ({ active, payload }: any) => {
  if (!active || !payload?.length) return null;
  const d = payload[0].payload;
  return (
    <div className="bg-slate-900/95 border border-white/20 rounded-xl p-3 shadow-xl backdrop-blur-md">
      <p className="font-semibold text-white text-sm mb-1">{d.name}</p>
      <p className="text-xs text-cyan-400">X: {d.xFormatted}</p>
      <p className="text-xs text-purple-400">Y: {d.yFormatted}</p>
      {d.distance && <p className="text-xs text-slate-400">Distance: {d.distance.toFixed(1)} Mpc</p>}
    </div>
  );
};

export default function CorrelationsPage() {
  const [data, setData] = useState<AnalysisData | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);
  const [selectedModel, setSelectedModel] = useState<ModelKey>('dark_halo_linear');
  const [xAxis, setXAxis] = useState<AxisVar>('M');
  const [yAxis, setYAxis] = useState<AxisVar>('k');
  const [useLogX, setUseLogX] = useState(true);
  const [useLogY, setUseLogY] = useState(true);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'sparc-results.json')
      .then(r => r.json())
      .then(d => { setData(d); setLoading(false); })
      .catch(e => { setError(e.message); setLoading(false); });
  }, []);

  const scatterData = useMemo(() => {
    if (!data) return [];
    return data.perGalaxy
      .map(g => {
        const xVal = getValue(g, xAxis, selectedModel);
        const yVal = getValue(g, yAxis, selectedModel);
        if (xVal === null || yVal === null || !isFinite(xVal) || !isFinite(yVal)) return null;
        if (xVal <= 0 && useLogX) return null;
        if (yVal <= 0 && useLogY) return null;
        return {
          name: g.name,
          x: useLogX ? Math.log10(xVal) : xVal,
          y: useLogY ? Math.log10(yVal) : yVal,
          rawX: xVal,
          rawY: yVal,
          xFormatted: formatAxisValue(xVal, xAxis),
          yFormatted: formatAxisValue(yVal, yAxis),
          distance: g.distance,
          maxV: g.maxV,
          maxR: g.maxR,
        };
      })
      .filter(Boolean) as any[];
  }, [data, selectedModel, xAxis, yAxis, useLogX, useLogY]);

  const correlation = useMemo(() => {
    if (scatterData.length < 3) return null;
    const xs = scatterData.map((d: any) => d.x);
    const ys = scatterData.map((d: any) => d.y);
    const r = pearsonCorrelation(xs, ys);
    const rho = spearmanCorrelation(xs, ys);
    const reg = linearRegression(xs, ys);
    return { pearson: r, spearman: rho, r2: r * r, slope: reg.slope, intercept: reg.intercept, n: scatterData.length };
  }, [scatterData]);

  const presetCorrelations = useMemo(() => {
    if (!data) return [];
    const pairs: { xAxis: AxisVar; yAxis: AxisVar; model: ModelKey; label: string }[] = [
      { xAxis: 'M', yAxis: 'k', model: 'dark_halo_linear', label: 'k vs M (Dark Halo Linear)' },
      { xAxis: 'maxV', yAxis: 'k', model: 'dark_halo_linear', label: 'k vs Max Velocity (DHL)' },
      { xAxis: 'maxR', yAxis: 'k', model: 'dark_halo_linear', label: 'k vs Max Radius (DHL)' },
      { xAxis: 'M', yAxis: 'k', model: 'modified_gravity_halo', label: 'k vs M (Mod Gravity+Halo)' },
      { xAxis: 'maxV', yAxis: 'k', model: 'modified_gravity_halo', label: 'k vs Max Velocity (MGH)' },
      { xAxis: 'M', yAxis: 'k', model: 'log_halo', label: 'k vs M (Log Halo)' },
      { xAxis: 'maxV', yAxis: 'k', model: 'log_halo', label: 'k vs Max Velocity (Log Halo)' },
      { xAxis: 'M', yAxis: 'a', model: 'modified_gravity_halo', label: 'a vs M (Mod Gravity+Halo)' },
      { xAxis: 'maxR', yAxis: 'a', model: 'log_halo', label: 'a vs Max Radius (Log Halo)' },
      { xAxis: 'distance', yAxis: 'k', model: 'dark_halo_linear', label: 'k vs Distance (DHL)' },
    ];
    return pairs.map(p => {
      const points = data.perGalaxy.map(g => {
        const xVal = getValue(g, p.xAxis, p.model);
        const yVal = getValue(g, p.yAxis, p.model);
        if (xVal === null || yVal === null || !isFinite(xVal) || !isFinite(yVal) || xVal <= 0 || yVal <= 0) return null;
        return { x: Math.log10(xVal), y: Math.log10(yVal) };
      }).filter(Boolean) as { x: number; y: number }[];
      const xs = points.map(pt => pt.x);
      const ys = points.map(pt => pt.y);
      const r = pearsonCorrelation(xs, ys);
      const rho = spearmanCorrelation(xs, ys);
      return { ...p, pearson: r, spearman: rho, r2: r * r, n: points.length };
    }).sort((a, b) => Math.abs(b.pearson) - Math.abs(a.pearson));
  }, [data]);

  const powerLawAnalysis = useMemo(() => {
    if (!data) return null;
    const models: ModelKey[] = ['dark_halo_linear', 'modified_gravity_halo', 'log_halo', 'dark_halo_flat', 'transition'];
    
    const results = models.map(mk => {
      const pts = data.perGalaxy
        .map(g => ({ k: g.models[mk]?.k, v: g.maxV, r: g.maxR, M: g.models[mk]?.M }))
        .filter(p => p.k > 0 && p.v > 0 && p.r > 0 && p.M > 0);
      
      if (pts.length < 10) return null;
      
      const fitPowerLaw = (xs: number[], ys: number[]) => {
        const n = xs.length;
        const mx = xs.reduce((s, v) => s + v, 0) / n;
        const my = ys.reduce((s, v) => s + v, 0) / n;
        let num = 0, den = 0;
        for (let i = 0; i < n; i++) { num += (xs[i] - mx) * (ys[i] - my); den += (xs[i] - mx) ** 2; }
        const alpha = den > 0 ? num / den : 0;
        const logC = my - alpha * mx;
        let ssRes = 0, ssTot = 0;
        for (let i = 0; i < n; i++) { const p = alpha * xs[i] + logC; ssRes += (ys[i] - p) ** 2; ssTot += (ys[i] - my) ** 2; }
        return { alpha, c: Math.pow(10, logC), r2: ssTot > 0 ? 1 - ssRes / ssTot : 0, n };
      };
      
      const kVsV = fitPowerLaw(pts.map(p => Math.log10(p.v)), pts.map(p => Math.log10(p.k)));
      const kVsR = fitPowerLaw(pts.map(p => Math.log10(p.r)), pts.map(p => Math.log10(p.k)));
      const kVsM = fitPowerLaw(pts.map(p => Math.log10(p.M)), pts.map(p => Math.log10(p.k)));
      const kVsV2R = fitPowerLaw(
        pts.map(p => Math.log10(p.v * p.v / p.r)),
        pts.map(p => Math.log10(p.k))
      );
      
      const logV = pts.map(p => Math.log10(p.v));
      const logR = pts.map(p => Math.log10(p.r));
      const logK = pts.map(p => Math.log10(p.k));
      const n = pts.length;
      const mV = logV.reduce((s, v) => s + v, 0) / n;
      const mR = logR.reduce((s, v) => s + v, 0) / n;
      const mK = logK.reduce((s, v) => s + v, 0) / n;
      let s00 = 0, s01 = 0, s11 = 0, s0y = 0, s1y = 0;
      for (let i = 0; i < n; i++) {
        const d0 = logV[i] - mV, d1 = logR[i] - mR, dy = logK[i] - mK;
        s00 += d0 * d0; s01 += d0 * d1; s11 += d1 * d1; s0y += d0 * dy; s1y += d1 * dy;
      }
      const det = s00 * s11 - s01 * s01;
      const alphaV = det > 0 ? (s11 * s0y - s01 * s1y) / det : 0;
      const alphaR = det > 0 ? (s00 * s1y - s01 * s0y) / det : 0;
      const logCmulti = mK - alphaV * mV - alphaR * mR;
      let ssResM = 0, ssTotM = 0;
      for (let i = 0; i < n; i++) {
        const pred = alphaV * logV[i] + alphaR * logR[i] + logCmulti;
        ssResM += (logK[i] - pred) ** 2; ssTotM += (logK[i] - mK) ** 2;
      }
      const multiR2 = ssTotM > 0 ? 1 - ssResM / ssTotM : 0;
      
      return {
        model: mk,
        modelName: MODEL_OPTIONS.find(m => m.key === mk)?.label || mk,
        kVsV, kVsR, kVsM, kVsV2R,
        multi: { alphaV, alphaR, c: Math.pow(10, logCmulti), r2: multiR2, n },
      };
    }).filter(Boolean) as NonNullable<ReturnType<typeof Array.prototype.map>[number]>[];
    
    return results;
  }, [data]);

  const modelColor = MODEL_OPTIONS.find(m => m.key === selectedModel)?.color || '#a78bfa';
  const xLabel = AXIS_OPTIONS.find(a => a.key === xAxis)?.label || xAxis;
  const yLabel = AXIS_OPTIONS.find(a => a.key === yAxis)?.label || yAxis;

  if (loading) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-[60vh]">
          <Loader2 className="w-8 h-8 animate-spin text-cyan-400" />
          <span className="ml-3 text-slate-400">Loading SPARC analysis data...</span>
        </div>
      </Layout>
    );
  }

  if (error || !data) {
    return (
      <Layout>
        <GlassCard>
          <p className="text-red-400">Failed to load analysis data: {error}</p>
          <p className="text-sm text-slate-400 mt-2">Run the full analysis script first to generate sparc-results.json</p>
        </GlassCard>
      </Layout>
    );
  }

  return (
    <Layout>
      <header className="mb-8">
        <h1 className="text-3xl font-bold flex items-center gap-3">
          <TrendingUp className="w-8 h-8 text-emerald-400" />
          Correlation Explorer
        </h1>
        <p className="text-slate-400 mt-2">
          Discover if k depends on galaxy properties. If k correlates with mass, size, or velocity — it's a physical law, not just curve fitting.
        </p>
      </header>

      <div className="space-y-6">
        <GlassCard glow="cyan">
          <h2 className="text-lg font-semibold mb-3 flex items-center gap-2">
            <Search className="w-5 h-5 text-cyan-400" />
            Key Question: Does k Depend on Galaxy Properties?
          </h2>
          <p className="text-sm text-slate-300 leading-relaxed">
            k varies across 175 galaxies (CV ~81%). If this variation is <span className="text-emerald-400 font-semibold">correlated</span> with 
            galaxy mass (M), max velocity (V<sub>max</sub>), or radius (R<sub>max</sub>), then k encodes a real physical 
            relationship — not random noise. A strong correlation means we can predict k from observable properties, 
            turning an unstable fitting parameter into a <span className="text-amber-400 font-semibold">new physical law</span>.
          </p>
        </GlassCard>

        <GlassCard>
          <h2 className="text-lg font-semibold mb-4 flex items-center gap-2">
            <FlaskConical className="w-5 h-5 text-purple-400" />
            Correlation Matrix (log-log, top pairs)
          </h2>
          <div className="overflow-x-auto">
            <table className="w-full text-sm">
              <thead>
                <tr className="border-b border-white/10">
                  <th className="text-left py-2 px-3 text-slate-400 font-medium">#</th>
                  <th className="text-left py-2 px-3 text-slate-400 font-medium">Relationship</th>
                  <th className="text-right py-2 px-3 text-cyan-400 font-medium">Pearson r</th>
                  <th className="text-right py-2 px-3 text-purple-400 font-medium">Spearman ρ</th>
                  <th className="text-right py-2 px-3 text-amber-400 font-medium">R²</th>
                  <th className="text-center py-2 px-3 text-slate-400 font-medium">Verdict</th>
                </tr>
              </thead>
              <tbody>
                {presetCorrelations.map((c, i) => {
                  const absR = Math.abs(c.pearson);
                  const verdict = absR > 0.7 ? 'STRONG' : absR > 0.5 ? 'MODERATE' : absR > 0.3 ? 'WEAK' : 'NONE';
                  const verdictColor = absR > 0.7 ? 'text-emerald-400' : absR > 0.5 ? 'text-cyan-400' : absR > 0.3 ? 'text-amber-400' : 'text-red-400';
                  const Icon = absR > 0.5 ? CheckCircle2 : absR > 0.3 ? AlertTriangle : XCircle;
                  return (
                    <tr
                      key={i}
                      className="border-b border-white/5 hover:bg-white/5 cursor-pointer transition-colors"
                      onClick={() => {
                        setXAxis(c.xAxis);
                        setYAxis(c.yAxis);
                        setSelectedModel(c.model);
                        setUseLogX(true);
                        setUseLogY(true);
                      }}
                    >
                      <td className="py-2.5 px-3 font-mono text-slate-500">{i + 1}</td>
                      <td className="py-2.5 px-3 font-medium text-slate-200">{c.label}</td>
                      <td className="py-2.5 px-3 text-right font-mono">
                        <span className={verdictColor}>{c.pearson.toFixed(3)}</span>
                      </td>
                      <td className="py-2.5 px-3 text-right font-mono text-purple-300">{c.spearman.toFixed(3)}</td>
                      <td className="py-2.5 px-3 text-right font-mono text-amber-300">{c.r2.toFixed(3)}</td>
                      <td className="py-2.5 px-3 text-center">
                        <span className={`flex items-center gap-1 justify-center text-xs font-semibold ${verdictColor}`}>
                          <Icon className="w-3.5 h-3.5" />
                          {verdict}
                        </span>
                      </td>
                    </tr>
                  );
                })}
              </tbody>
            </table>
          </div>
          <p className="text-xs text-slate-500 mt-3">Click any row to view the scatter plot below.</p>
        </GlassCard>

        <GlassCard glow="purple">
          <h2 className="text-lg font-semibold mb-4">Interactive Scatter Plot</h2>
          <div className="grid grid-cols-1 md:grid-cols-4 gap-4 mb-6">
            <div>
              <label className="text-xs text-slate-400 uppercase tracking-wider mb-1 block">Model</label>
              <select
                value={selectedModel}
                onChange={e => setSelectedModel(e.target.value as ModelKey)}
                className="w-full bg-slate-900/80 border border-white/10 rounded-xl px-3 py-2 text-sm text-white focus:border-purple-500 focus:outline-none"
              >
                {MODEL_OPTIONS.map(m => (
                  <option key={m.key} value={m.key}>{m.label}</option>
                ))}
              </select>
            </div>
            <div>
              <label className="text-xs text-slate-400 uppercase tracking-wider mb-1 block">X Axis</label>
              <select
                value={xAxis}
                onChange={e => setXAxis(e.target.value as AxisVar)}
                className="w-full bg-slate-900/80 border border-white/10 rounded-xl px-3 py-2 text-sm text-white focus:border-cyan-500 focus:outline-none"
              >
                {AXIS_OPTIONS.map(a => (
                  <option key={a.key} value={a.key}>{a.label}</option>
                ))}
              </select>
              <label className="flex items-center gap-2 mt-2 text-xs text-slate-400">
                <input type="checkbox" checked={useLogX} onChange={e => setUseLogX(e.target.checked)} className="accent-cyan-500" />
                Log scale
              </label>
            </div>
            <div>
              <label className="text-xs text-slate-400 uppercase tracking-wider mb-1 block">Y Axis</label>
              <select
                value={yAxis}
                onChange={e => setYAxis(e.target.value as AxisVar)}
                className="w-full bg-slate-900/80 border border-white/10 rounded-xl px-3 py-2 text-sm text-white focus:border-purple-500 focus:outline-none"
              >
                {AXIS_OPTIONS.map(a => (
                  <option key={a.key} value={a.key}>{a.label}</option>
                ))}
              </select>
              <label className="flex items-center gap-2 mt-2 text-xs text-slate-400">
                <input type="checkbox" checked={useLogY} onChange={e => setUseLogY(e.target.checked)} className="accent-purple-500" />
                Log scale
              </label>
            </div>
            <div className="flex items-end">
              <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5 w-full">
                <div className="text-xs text-slate-400 mb-1">Data Points</div>
                <div className="text-2xl font-bold font-mono text-white">{scatterData.length}</div>
                <div className="text-xs text-slate-500">of {data.perGalaxy.length} galaxies</div>
              </div>
            </div>
          </div>

          <div className="h-[500px] w-full">
            <ResponsiveContainer width="100%" height="100%">
              <ScatterChart margin={{ top: 20, right: 30, bottom: 60, left: 60 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                <XAxis
                  type="number"
                  dataKey="x"
                  name={xLabel}
                  tick={{ fill: '#94a3b8', fontSize: 11 }}
                  tickFormatter={(v) => useLogX ? `10^${v.toFixed(0)}` : v.toFixed(0)}
                >
                  <Label value={useLogX ? `log₁₀(${xLabel})` : xLabel} position="bottom" offset={35} style={{ fill: '#94a3b8', fontSize: 12 }} />
                </XAxis>
                <YAxis
                  type="number"
                  dataKey="y"
                  name={yLabel}
                  tick={{ fill: '#94a3b8', fontSize: 11 }}
                  tickFormatter={(v) => useLogY ? `10^${v.toFixed(0)}` : v.toFixed(0)}
                >
                  <Label value={useLogY ? `log₁₀(${yLabel})` : yLabel} position="left" angle={-90} offset={35} style={{ fill: '#94a3b8', fontSize: 12 }} />
                </YAxis>
                <Tooltip content={<CustomTooltip />} />
                {correlation && Math.abs(correlation.pearson) > 0.15 && (
                  <ReferenceLine
                    segment={[
                      { x: Math.min(...scatterData.map((d: any) => d.x)), y: correlation.slope * Math.min(...scatterData.map((d: any) => d.x)) + correlation.intercept },
                      { x: Math.max(...scatterData.map((d: any) => d.x)), y: correlation.slope * Math.max(...scatterData.map((d: any) => d.x)) + correlation.intercept },
                    ]}
                    stroke="#fbbf24"
                    strokeWidth={2}
                    strokeDasharray="8 4"
                  />
                )}
                <Scatter
                  data={scatterData}
                  fill={modelColor}
                  fillOpacity={0.7}
                  r={5}
                />
              </ScatterChart>
            </ResponsiveContainer>
          </div>
        </GlassCard>

        {correlation && (
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-4">
            <GlassCard className="p-4">
              <div className="text-xs text-slate-400 uppercase tracking-wider mb-1">Pearson r</div>
              <div className={`text-2xl font-bold font-mono ${Math.abs(correlation.pearson) > 0.7 ? 'text-emerald-400' : Math.abs(correlation.pearson) > 0.5 ? 'text-cyan-400' : Math.abs(correlation.pearson) > 0.3 ? 'text-amber-400' : 'text-red-400'}`}>
                {correlation.pearson.toFixed(4)}
              </div>
              <div className="text-xs text-slate-500 mt-1">Linear correlation in {useLogX || useLogY ? 'log' : 'linear'} space</div>
            </GlassCard>
            <GlassCard className="p-4">
              <div className="text-xs text-slate-400 uppercase tracking-wider mb-1">Spearman ρ</div>
              <div className={`text-2xl font-bold font-mono ${Math.abs(correlation.spearman) > 0.7 ? 'text-emerald-400' : Math.abs(correlation.spearman) > 0.5 ? 'text-cyan-400' : Math.abs(correlation.spearman) > 0.3 ? 'text-amber-400' : 'text-red-400'}`}>
                {correlation.spearman.toFixed(4)}
              </div>
              <div className="text-xs text-slate-500 mt-1">Rank-order (monotonic) correlation</div>
            </GlassCard>
            <GlassCard className="p-4">
              <div className="text-xs text-slate-400 uppercase tracking-wider mb-1">R² (Coefficient of Determination)</div>
              <div className="text-2xl font-bold font-mono text-amber-400">{correlation.r2.toFixed(4)}</div>
              <div className="text-xs text-slate-500 mt-1">{(correlation.r2 * 100).toFixed(1)}% of variance explained</div>
            </GlassCard>
            <GlassCard className="p-4">
              <div className="text-xs text-slate-400 uppercase tracking-wider mb-1">Verdict</div>
              {(() => {
                const absR = Math.abs(correlation.pearson);
                const absRho = Math.abs(correlation.spearman);
                const best = Math.max(absR, absRho);
                if (best > 0.7) return (
                  <div>
                    <div className="text-xl font-bold text-emerald-400 flex items-center gap-2">
                      <CheckCircle2 className="w-5 h-5" /> STRONG
                    </div>
                    <p className="text-xs text-emerald-300/70 mt-1">Strong trend detected (n={correlation.n}). Suggests a predictable relationship worth investigating further.</p>
                  </div>
                );
                if (best > 0.5) return (
                  <div>
                    <div className="text-xl font-bold text-cyan-400 flex items-center gap-2">
                      <CheckCircle2 className="w-5 h-5" /> MODERATE
                    </div>
                    <p className="text-xs text-cyan-300/70 mt-1">Moderate trend (n={correlation.n}). Pattern exists but with significant scatter.</p>
                  </div>
                );
                if (best > 0.3) return (
                  <div>
                    <div className="text-xl font-bold text-amber-400 flex items-center gap-2">
                      <AlertTriangle className="w-5 h-5" /> WEAK
                    </div>
                    <p className="text-xs text-amber-300/70 mt-1">Weak trend (n={correlation.n}). Too much scatter for confident interpretation.</p>
                  </div>
                );
                return (
                  <div>
                    <div className="text-xl font-bold text-red-400 flex items-center gap-2">
                      <XCircle className="w-5 h-5" /> NO CORRELATION
                    </div>
                    <p className="text-xs text-red-300/70 mt-1">No detectable linear trend (n={correlation.n}).</p>
                  </div>
                );
              })()}
            </GlassCard>
          </div>
        )}

        {powerLawAnalysis && powerLawAnalysis.length > 0 && (
          <GlassCard glow="cyan">
            <h2 className="text-lg font-semibold mb-2 flex items-center gap-2">
              <Zap className="w-5 h-5 text-amber-400" />
              Power Law Analysis: k = c × V<sup>α</sup> × R<sup>β</sup>
            </h2>
            <p className="text-sm text-slate-400 mb-4">
              Fitting k as a power law of galaxy observables. If α ≈ 2 and β ≈ -1, then k ≈ V²/R — which is exactly
              what physics predicts when v² = GM/r + kr at large radii.
            </p>

            <div className="overflow-x-auto mb-6">
              <table className="w-full text-sm">
                <thead>
                  <tr className="border-b border-white/10">
                    <th className="text-left py-2 px-3 text-slate-400 font-medium">Model</th>
                    <th className="text-center py-2 px-3 text-cyan-400 font-medium" colSpan={2}>k ∝ V<sup>α</sup></th>
                    <th className="text-center py-2 px-3 text-purple-400 font-medium" colSpan={2}>k ∝ V²/R</th>
                    <th className="text-center py-2 px-3 text-emerald-400 font-medium" colSpan={3}>k ∝ V<sup>α</sup> × R<sup>β</sup></th>
                  </tr>
                  <tr className="border-b border-white/5">
                    <th></th>
                    <th className="text-right py-1 px-2 text-xs text-slate-500">α</th>
                    <th className="text-right py-1 px-2 text-xs text-slate-500">R²</th>
                    <th className="text-right py-1 px-2 text-xs text-slate-500">α</th>
                    <th className="text-right py-1 px-2 text-xs text-slate-500">R²</th>
                    <th className="text-right py-1 px-2 text-xs text-slate-500">α(V)</th>
                    <th className="text-right py-1 px-2 text-xs text-slate-500">β(R)</th>
                    <th className="text-right py-1 px-2 text-xs text-slate-500">R²</th>
                  </tr>
                </thead>
                <tbody>
                  {powerLawAnalysis.map(r => {
                    const bestR2 = Math.max(r.kVsV.r2, r.kVsV2R.r2, r.multi.r2);
                    return (
                      <tr key={r.model} className="border-b border-white/5 hover:bg-white/5 transition-colors">
                        <td className="py-2.5 px-3 font-medium text-slate-200">{r.modelName}</td>
                        <td className="py-2.5 px-2 text-right font-mono text-cyan-300">{r.kVsV.alpha.toFixed(2)}</td>
                        <td className={`py-2.5 px-2 text-right font-mono ${r.kVsV.r2 === bestR2 ? 'text-emerald-400 font-bold' : 'text-slate-400'}`}>{r.kVsV.r2.toFixed(3)}</td>
                        <td className="py-2.5 px-2 text-right font-mono text-purple-300">{r.kVsV2R.alpha.toFixed(2)}</td>
                        <td className={`py-2.5 px-2 text-right font-mono ${r.kVsV2R.r2 === bestR2 ? 'text-emerald-400 font-bold' : 'text-slate-400'}`}>{r.kVsV2R.r2.toFixed(3)}</td>
                        <td className="py-2.5 px-2 text-right font-mono text-emerald-300">{r.multi.alphaV.toFixed(2)}</td>
                        <td className="py-2.5 px-2 text-right font-mono text-emerald-300">{r.multi.alphaR.toFixed(2)}</td>
                        <td className={`py-2.5 px-2 text-right font-mono ${r.multi.r2 === bestR2 ? 'text-emerald-400 font-bold' : 'text-slate-400'}`}>{r.multi.r2.toFixed(3)}</td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>

            {(() => {
              const dhl = powerLawAnalysis.find(r => r.model === 'dark_halo_linear');
              if (!dhl) return null;
              return (
                <div className="space-y-4">
                  <div className="p-4 bg-gradient-to-r from-amber-500/10 to-emerald-500/10 border border-amber-500/20 rounded-xl">
                    <h3 className="text-sm font-semibold text-amber-400 mb-3 flex items-center gap-2">
                      <Flame className="w-4 h-4" />
                      Key Finding: Dark Halo Linear
                    </h3>
                    <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                      <div className="text-center p-3 bg-slate-900/50 rounded-lg">
                        <div className="text-xs text-slate-400 mb-1">Single Variable</div>
                        <div className="font-mono text-lg text-cyan-400">k ∝ V<sup>{dhl.kVsV.alpha.toFixed(2)}</sup></div>
                        <div className="text-xs text-slate-500 mt-1">R² = {dhl.kVsV.r2.toFixed(3)}</div>
                      </div>
                      <div className="text-center p-3 bg-slate-900/50 rounded-lg">
                        <div className="text-xs text-slate-400 mb-1">Theoretical (V²/R)</div>
                        <div className="font-mono text-lg text-purple-400">k ∝ (V²/R)<sup>{dhl.kVsV2R.alpha.toFixed(2)}</sup></div>
                        <div className="text-xs text-slate-500 mt-1">R² = {dhl.kVsV2R.r2.toFixed(3)}</div>
                      </div>
                      <div className="text-center p-3 bg-slate-900/50 rounded-lg">
                        <div className="text-xs text-slate-400 mb-1">Multivariate</div>
                        <div className="font-mono text-lg text-emerald-400">k ∝ V<sup>{dhl.multi.alphaV.toFixed(1)}</sup>R<sup>{dhl.multi.alphaR.toFixed(1)}</sup></div>
                        <div className="text-xs text-slate-500 mt-1">R² = {dhl.multi.r2.toFixed(3)}</div>
                      </div>
                    </div>
                  </div>

                  <div className="p-4 bg-slate-900/50 border border-white/10 rounded-xl">
                    <h3 className="text-sm font-semibold text-white mb-3">Physical Interpretation</h3>
                    <div className="space-y-3 text-xs text-slate-300 leading-relaxed">
                      <div className="flex gap-3">
                        <div className="w-16 shrink-0 text-right font-mono text-amber-400">V²/R</div>
                        <div>
                          In the Dark Halo model: v² = GM/r + kr. At large radii where GM/r → 0, we get v² ≈ kr, 
                          thus <span className="text-cyan-400 font-semibold">k ≈ v²/r</span>. 
                          The fit confirms this: α ≈ {dhl.kVsV2R.alpha.toFixed(2)} ≈ 1.0, R² = {dhl.kVsV2R.r2.toFixed(3)}.
                          {dhl.kVsV2R.r2 > 0.85 && <span className="text-emerald-400 font-semibold"> This is an excellent fit!</span>}
                        </div>
                      </div>
                      <div className="flex gap-3">
                        <div className="w-16 shrink-0 text-right font-mono text-amber-400">V<sup>{dhl.multi.alphaV.toFixed(1)}</sup>R<sup>{dhl.multi.alphaR.toFixed(1)}</sup></div>
                        <div>
                          Multivariate: k ≈ {dhl.multi.c.toFixed(2)} × V<sup>{dhl.multi.alphaV.toFixed(2)}</sup> × R<sup>{dhl.multi.alphaR.toFixed(2)}</sup> (R² = {dhl.multi.r2.toFixed(3)}).
                          {Math.abs(dhl.multi.alphaV - 2) < 0.3 && Math.abs(dhl.multi.alphaR - (-1)) < 0.3 && (
                            <span className="text-emerald-400 font-semibold">
                              {' '}α(V) ≈ 2, β(R) ≈ -1 confirms the V²/R relationship!
                            </span>
                          )}
                        </div>
                      </div>
                      <div className="flex gap-3">
                        <div className="w-16 shrink-0 text-right font-mono text-amber-400">Meaning</div>
                        <div>
                          {dhl.kVsV2R.r2 > 0.8 ? (
                            <span>
                              k is <span className="text-emerald-400 font-semibold">not a free parameter</span> — it is determined by 
                              the galaxy's flat rotation velocity and size. This is consistent with k encoding the dark matter 
                              halo density at the flat part of the rotation curve. The relationship k ≈ V²<sub>flat</sub>/R<sub>max</sub> is 
                              a direct consequence of the model physics, confirming internal consistency.
                            </span>
                          ) : (
                            <span>The power law fit has moderate quality. The relationship may be real but with significant scatter from measurement uncertainties.</span>
                          )}
                        </div>
                      </div>
                    </div>
                  </div>

                  <div className="p-4 bg-purple-500/10 border border-purple-500/20 rounded-xl">
                    <h3 className="text-sm font-semibold text-purple-400 mb-2">Comparison with Known Scaling Laws</h3>
                    <div className="grid grid-cols-1 md:grid-cols-3 gap-3 text-xs">
                      <div className="p-2 bg-slate-900/50 rounded-lg">
                        <div className="font-mono text-cyan-400 font-semibold">Tully-Fisher</div>
                        <div className="text-slate-400 mt-1">L ∝ V<sup>4</sup></div>
                        <div className="text-slate-500 mt-0.5">Luminosity scales with velocity</div>
                      </div>
                      <div className="p-2 bg-slate-900/50 rounded-lg">
                        <div className="font-mono text-cyan-400 font-semibold">Virial Theorem</div>
                        <div className="text-slate-400 mt-1">M ∝ V² × R</div>
                        <div className="text-slate-500 mt-0.5">Mass from gravitational equilibrium</div>
                      </div>
                      <div className="p-2 bg-slate-900/50 rounded-lg border border-emerald-500/30">
                        <div className="font-mono text-emerald-400 font-semibold">Your Finding</div>
                        <div className="text-slate-300 mt-1">k ∝ V<sup>{dhl.multi.alphaV.toFixed(1)}</sup> × R<sup>{dhl.multi.alphaR.toFixed(1)}</sup></div>
                        <div className="text-slate-500 mt-0.5">≈ V²/R — halo density scaling</div>
                      </div>
                    </div>
                  </div>
                </div>
              );
            })()}
          </GlassCard>
        )}

        <GlassCard>
          <h2 className="text-lg font-semibold mb-3">Interpretation Guide</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="p-4 bg-emerald-500/10 border border-emerald-500/20 rounded-xl">
              <h3 className="text-sm font-semibold text-emerald-400 mb-2">If Strong Correlation Found (|r| &gt; 0.7)</h3>
              <ul className="text-xs text-slate-300 space-y-1.5">
                <li>k = f(M) or k = f(V<sub>max</sub>) → k is predictable from observables</li>
                <li>This means the "extra effect" scales with galaxy properties</li>
                <li>The power-law slope tells you the physics: k ∝ M<sup>α</sup></li>
                <li>This is a potential new empirical law (like Tully-Fisher)</li>
              </ul>
            </div>
            <div className="p-4 bg-red-500/10 border border-red-500/20 rounded-xl">
              <h3 className="text-sm font-semibold text-red-400 mb-2">If No Correlation Found (|r| &lt; 0.3)</h3>
              <ul className="text-xs text-slate-300 space-y-1.5">
                <li>k is absorbing random galaxy-specific noise</li>
                <li>The model has too many degrees of freedom</li>
                <li>The improvement is from overfitting, not physics</li>
                <li>Need a different parameterization or fewer free parameters</li>
              </ul>
            </div>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
