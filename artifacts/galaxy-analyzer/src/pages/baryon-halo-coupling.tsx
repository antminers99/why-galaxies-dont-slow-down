import React, { useEffect, useState } from 'react';
import { Layout } from '@/components/layout';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, BarChart, Bar, Legend, ReferenceLine, Cell, Label } from 'recharts';
import { Link2, CheckCircle2, XCircle, Layers, ArrowLeftRight, TrendingUp, Activity, Orbit, Target, GitBranch } from 'lucide-react';

interface FingerprintItem {
  metric: string;
  n: number;
  slope: number;
  slopeErr: number;
  r: number;
  r2: number;
  partialR: number;
  bootMean: number;
  bootCI95: number[];
  expectedSign: string;
  matchesExpected: boolean;
}

interface InnerOuterData {
  innerVDM: { slope: number; r: number; partialR: number; n: number };
  outerVDM: { slope: number; r: number; partialR: number; n: number };
  innerFDM: { slope: number; r: number; n: number };
  outerFDM: { slope: number; r: number; n: number };
  innerStronger: { vdm: boolean; fdm: boolean };
}

interface ResidualTest {
  baselineR2: number;
  residualSlope: number;
  residualIntercept: number;
  residualR: number;
  residualR2: number;
  n: number;
  sigBarExplains: boolean;
}

interface CouplingLaw {
  intercept: number;
  slope: number;
  r2: number;
  r: number;
  formula: string;
}

interface MassBin {
  label: string;
  n: number;
  slope?: number;
  r?: number;
  skipped: boolean;
}

interface PlotPoint {
  logSigBar: number;
  logVDM?: number;
  rDMdomNorm?: number;
  alpha?: number;
  residual?: number;
  vmax?: number;
  name?: string;
}

interface CouplingAnalysis {
  fingerprint: FingerprintItem[];
  innerOuter: InnerOuterData;
  residualTests: { vdm: ResidualTest; rDMdom: ResidualTest; alpha: ResidualTest };
  couplingLaws: { vdm: CouplingLaw; rDMdom: CouplingLaw; alpha: CouplingLaw };
  massBins: MassBin[];
  tests: { name: string; pass: boolean }[];
  passCount: number;
  couplingWins: boolean;
  verdict: string;
  plotData: {
    vdm: PlotPoint[];
    rDMdom: PlotPoint[];
    alpha: PlotPoint[];
    residuals: { vdm: PlotPoint[]; rDMdom: PlotPoint[]; alpha: PlotPoint[] };
  };
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

function SectionHeader({ icon: Icon, title, subtitle }: { icon: any; title: string; subtitle: string }) {
  return (
    <div className="flex items-center gap-3 mb-6">
      <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500/20 to-violet-500/20 flex items-center justify-center border border-cyan-500/20">
        <Icon className="w-5 h-5 text-cyan-400" />
      </div>
      <div>
        <h2 className="text-lg font-display font-bold text-white">{title}</h2>
        <p className="text-xs text-slate-500">{subtitle}</p>
      </div>
    </div>
  );
}

const vmaxColor = (vmax: number) => {
  if (vmax < 60) return '#22d3ee';
  if (vmax < 100) return '#06b6d4';
  if (vmax < 150) return '#8b5cf6';
  if (vmax < 200) return '#a78bfa';
  return '#f472b6';
};

function CouplingScatterPlot({ data, xKey, yKey, xLabel, yLabel, reg, color = '#06b6d4' }: {
  data: PlotPoint[];
  xKey: string;
  yKey: string;
  xLabel: string;
  yLabel: string;
  reg?: { slope: number; intercept: number; r: number };
  color?: string;
}) {
  const chartData = data.filter(d => d[xKey as keyof PlotPoint] != null && d[yKey as keyof PlotPoint] != null);

  const trendData = reg ? (() => {
    const xs = chartData.map(d => d[xKey as keyof PlotPoint] as number);
    const xMin = Math.min(...xs);
    const xMax = Math.max(...xs);
    return [
      { x: xMin, y: reg.intercept + reg.slope * xMin },
      { x: xMax, y: reg.intercept + reg.slope * xMax },
    ];
  })() : [];

  return (
    <ResponsiveContainer width="100%" height={300}>
      <ScatterChart margin={{ top: 10, right: 30, bottom: 25, left: 20 }}>
        <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
        <XAxis type="number" dataKey="x" name={xLabel} stroke="#64748b" tick={{ fill: '#94a3b8', fontSize: 11 }}>
          <Label value={xLabel} position="bottom" offset={10} style={{ fill: '#94a3b8', fontSize: 12 }} />
        </XAxis>
        <YAxis type="number" dataKey="y" name={yLabel} stroke="#64748b" tick={{ fill: '#94a3b8', fontSize: 11 }}>
          <Label value={yLabel} angle={-90} position="insideLeft" offset={-5} style={{ fill: '#94a3b8', fontSize: 12 }} />
        </YAxis>
        <Tooltip
          content={({ active, payload }) => {
            if (!active || !payload?.length) return null;
            const d = payload[0].payload;
            return (
              <div className="bg-slate-900/95 border border-white/10 rounded-lg p-3 text-xs">
                {d.name && <div className="text-white font-medium mb-1">{d.name}</div>}
                <div className="text-cyan-400">{xLabel}: {d.x?.toFixed(3)}</div>
                <div className="text-emerald-400">{yLabel}: {d.y?.toFixed(4)}</div>
                {d.vmax && <div className="text-slate-400">V_max: {d.vmax} km/s</div>}
              </div>
            );
          }}
        />
        <Scatter
          data={chartData.map(d => ({
            x: d[xKey as keyof PlotPoint],
            y: d[yKey as keyof PlotPoint],
            name: d.name,
            vmax: d.vmax,
          }))}
          fill={color}
          fillOpacity={0.6}
          r={4}
        >
          {chartData.map((d, i) => (
            <Cell key={i} fill={d.vmax ? vmaxColor(d.vmax) : color} />
          ))}
        </Scatter>
        {reg && (
          <Scatter
            data={trendData}
            fill="none"
            line={{ stroke: '#f59e0b', strokeWidth: 2, strokeDasharray: '6 3' }}
            legendType="none"
            r={0}
          />
        )}
      </ScatterChart>
    </ResponsiveContainer>
  );
}

export default function BaryonHaloCouplingPage() {
  const [coupling, setCoupling] = useState<CouplingAnalysis | null>(null);
  const [activePhase, setActivePhase] = useState<1 | 2 | 3>(1);

  useEffect(() => {
    fetch(`${import.meta.env.BASE_URL}fdm-analysis.json`)
      .then(r => r.json())
      .then(d => { if (d.couplingAnalysis) setCoupling(d.couplingAnalysis); })
      .catch(() => {});
  }, []);

  if (!coupling) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-96">
          <div className="animate-spin w-8 h-8 border-2 border-cyan-500 border-t-transparent rounded-full" />
        </div>
      </Layout>
    );
  }

  const fp = coupling.fingerprint;
  const io = coupling.innerOuter;
  const rt = coupling.residualTests;
  const laws = coupling.couplingLaws;
  const pd = coupling.plotData;

  return (
    <Layout>
      <div className="space-y-8">
        <div className="flex items-start justify-between">
          <div>
            <h1 className="text-3xl font-display font-bold text-white flex items-center gap-3">
              <div className="w-12 h-12 rounded-2xl bg-gradient-to-br from-violet-500 to-cyan-500 flex items-center justify-center shadow-lg shadow-violet-500/20">
                <Link2 className="w-7 h-7 text-white" />
              </div>
              Baryon–Halo Coupling
            </h1>
            <p className="text-slate-400 mt-2 max-w-2xl">
              Does baryonic surface density track dark matter halo structure itself — not just f_DM?
              If V_DM, r_DMdom, and inner slope all correlate with Sigma_bar after removing Vmax dependence,
              baryons couple directly to halo properties.
            </p>
          </div>
          <div className={`px-5 py-3 rounded-2xl text-sm font-bold ${coupling.couplingWins ? 'bg-emerald-500/20 text-emerald-400 border border-emerald-500/30' : 'bg-amber-500/20 text-amber-400 border border-amber-500/30'}`}>
            {coupling.passCount}/6 TESTS PASS
          </div>
        </div>

        <GlassCard glow="emerald" className="border-emerald-500/20">
          <div className="flex items-center gap-2 mb-3">
            <Target className="w-5 h-5 text-emerald-400" />
            <h3 className="text-white font-bold">Verdict</h3>
          </div>
          <p className="text-emerald-300 text-sm leading-relaxed font-mono">{coupling.verdict}</p>
          <div className="flex flex-wrap gap-2 mt-4">
            {coupling.tests.map((t, i) => (
              <Badge key={i} pass={t.pass} label={t.name} />
            ))}
          </div>
        </GlassCard>

        <div className="flex gap-2">
          {[1, 2, 3].map(phase => (
            <button
              key={phase}
              onClick={() => setActivePhase(phase as 1 | 2 | 3)}
              className={`px-4 py-2 rounded-xl text-sm font-medium transition-all ${activePhase === phase ? 'bg-cyan-500/20 text-cyan-400 border border-cyan-500/30' : 'bg-white/5 text-slate-400 border border-white/10 hover:text-white'}`}
            >
              {phase === 1 ? 'Coupling Fingerprint' : phase === 2 ? 'Inner vs Outer' : 'Residual Tests'}
            </button>
          ))}
        </div>

        {activePhase === 1 && (
          <div className="space-y-6">
            <SectionHeader icon={Activity} title="Phase 1: Coupling Fingerprint" subtitle="Raw correlations + partial r controlling for Vmax + bootstrap CI" />

            <div className="grid grid-cols-3 gap-4">
              {fp.map((f, i) => (
                <GlassCard key={i} glow={f.matchesExpected ? 'cyan' : 'rose'}>
                  <div className="flex items-center justify-between mb-3">
                    <span className="text-sm font-bold text-white">{f.metric} vs log(Sigma_bar)</span>
                    <Badge pass={f.matchesExpected} label={f.matchesExpected ? 'Expected sign' : 'Wrong sign'} />
                  </div>
                  <div className="grid grid-cols-2 gap-3 mb-3">
                    <StatBox label="Raw r" value={f.r.toFixed(4)} color={Math.abs(f.r) > 0.5 ? 'text-emerald-400' : Math.abs(f.r) > 0.3 ? 'text-cyan-400' : 'text-amber-400'} />
                    <StatBox label="Partial r|Vmax" value={f.partialR.toFixed(4)} color={Math.abs(f.partialR) > 0.3 ? 'text-emerald-400' : Math.abs(f.partialR) > 0.15 ? 'text-cyan-400' : 'text-slate-400'} />
                    <StatBox label="Slope" value={f.slope.toFixed(4)} sub={`± ${f.slopeErr.toFixed(4)}`} />
                    <StatBox label="Bootstrap" value={f.bootMean.toFixed(4)} sub={`[${f.bootCI95[0].toFixed(3)}, ${f.bootCI95[1].toFixed(3)}]`} />
                  </div>
                  <div className="text-xs text-slate-500 text-center">n = {f.n} galaxies</div>
                </GlassCard>
              ))}
            </div>

            <div className="grid grid-cols-3 gap-6">
              <GlassCard>
                <h4 className="text-sm font-bold text-white mb-2">V_DM vs Sigma_bar</h4>
                <CouplingScatterPlot
                  data={pd.vdm}
                  xKey="logSigBar"
                  yKey="logVDM"
                  xLabel="log(Sigma_bar)"
                  yLabel="log(V_DM)"
                  reg={laws.vdm}
                />
                <p className="text-xs text-slate-500 mt-1 text-center font-mono">
                  R² = {laws.vdm.r2.toFixed(3)} | {laws.vdm.formula}
                </p>
              </GlassCard>

              <GlassCard>
                <h4 className="text-sm font-bold text-white mb-2">r_DMdom vs Sigma_bar</h4>
                <CouplingScatterPlot
                  data={pd.rDMdom}
                  xKey="logSigBar"
                  yKey="rDMdomNorm"
                  xLabel="log(Sigma_bar)"
                  yLabel="r_DMdom / r_max"
                  reg={laws.rDMdom}
                />
                <p className="text-xs text-slate-500 mt-1 text-center font-mono">
                  R² = {laws.rDMdom.r2.toFixed(3)} | {laws.rDMdom.formula}
                </p>
              </GlassCard>

              <GlassCard>
                <h4 className="text-sm font-bold text-white mb-2">Inner Slope vs Sigma_bar</h4>
                <CouplingScatterPlot
                  data={pd.alpha}
                  xKey="logSigBar"
                  yKey="alpha"
                  xLabel="log(Sigma_bar)"
                  yLabel="alpha (V_DM ~ r^alpha)"
                  reg={laws.alpha}
                  color="#a78bfa"
                />
                <p className="text-xs text-slate-500 mt-1 text-center font-mono">
                  R² = {laws.alpha.r2.toFixed(3)} | {laws.alpha.formula}
                </p>
              </GlassCard>
            </div>

            <GlassCard glow="violet">
              <SectionHeader icon={Orbit} title="Coupling Laws" subtitle="Empirical relationships between baryon surface density and halo properties" />
              <div className="space-y-3">
                {[
                  { law: laws.vdm, label: 'V_DM' },
                  { law: laws.rDMdom, label: 'r_DMdom' },
                  { law: laws.alpha, label: 'alpha' },
                ].map((item, i) => (
                  <div key={i} className="flex items-center gap-4 bg-white/5 rounded-xl p-4">
                    <div className="flex-1">
                      <span className="font-mono text-sm text-cyan-400">{item.law.formula}</span>
                    </div>
                    <div className="flex gap-6 text-xs">
                      <div>
                        <span className="text-slate-500">Slope: </span>
                        <span className="text-white font-mono">{item.law.slope.toFixed(4)}</span>
                      </div>
                      <div>
                        <span className="text-slate-500">r: </span>
                        <span className="text-white font-mono">{item.law.r.toFixed(4)}</span>
                      </div>
                      <div>
                        <span className="text-slate-500">R²: </span>
                        <span className="text-white font-mono">{item.law.r2.toFixed(3)}</span>
                      </div>
                    </div>
                  </div>
                ))}
              </div>
            </GlassCard>

            <GlassCard>
              <SectionHeader icon={Layers} title="Mass-Binned Coupling" subtitle="Does V_DM coupling persist across mass regimes?" />
              <div className="grid grid-cols-2 gap-4">
                {coupling.massBins.filter(b => !b.skipped).map((b, i) => (
                  <div key={i} className="bg-white/5 rounded-xl p-4">
                    <div className="text-sm font-bold text-white mb-2">{b.label}</div>
                    <div className="grid grid-cols-3 gap-3">
                      <StatBox label="n" value={String(b.n)} />
                      <StatBox label="slope" value={b.slope?.toFixed(4) || '-'} />
                      <StatBox label="r" value={b.r?.toFixed(4) || '-'} color={Math.abs(b.r || 0) > 0.3 ? 'text-cyan-400' : 'text-slate-400'} />
                    </div>
                  </div>
                ))}
              </div>
            </GlassCard>
          </div>
        )}

        {activePhase === 2 && (
          <div className="space-y-6">
            <SectionHeader icon={ArrowLeftRight} title="Phase 2: Inner vs Outer Coupling" subtitle="Where in the galaxy is baryon-halo coupling strongest?" />

            <div className="grid grid-cols-2 gap-6">
              <GlassCard glow={io.innerStronger.vdm ? 'emerald' : 'amber'}>
                <h4 className="text-sm font-bold text-white mb-4 flex items-center gap-2">
                  <Activity className="w-4 h-4 text-cyan-400" />
                  V_DM Coupling by Region
                </h4>
                <div className="space-y-4">
                  <div className="bg-white/5 rounded-xl p-4">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-sm text-slate-300">Inner (r {'<'} 0.5 r_max)</span>
                      {io.innerStronger.vdm && <span className="text-xs text-emerald-400 bg-emerald-500/10 px-2 py-0.5 rounded">Stronger</span>}
                    </div>
                    <div className="grid grid-cols-3 gap-2">
                      <StatBox label="r" value={io.innerVDM.r.toFixed(4)} color="text-cyan-400" />
                      <StatBox label="partial r" value={io.innerVDM.partialR.toFixed(4)} />
                      <StatBox label="n" value={String(io.innerVDM.n)} />
                    </div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-sm text-slate-300">Outer (r {'≥'} 0.5 r_max)</span>
                      {!io.innerStronger.vdm && <span className="text-xs text-emerald-400 bg-emerald-500/10 px-2 py-0.5 rounded">Stronger</span>}
                    </div>
                    <div className="grid grid-cols-3 gap-2">
                      <StatBox label="r" value={io.outerVDM.r.toFixed(4)} color="text-cyan-400" />
                      <StatBox label="partial r" value={io.outerVDM.partialR.toFixed(4)} />
                      <StatBox label="n" value={String(io.outerVDM.n)} />
                    </div>
                  </div>
                </div>
              </GlassCard>

              <GlassCard glow={io.innerStronger.fdm ? 'emerald' : 'amber'}>
                <h4 className="text-sm font-bold text-white mb-4 flex items-center gap-2">
                  <TrendingUp className="w-4 h-4 text-violet-400" />
                  f_DM Coupling by Region
                </h4>
                <div className="space-y-4">
                  <div className="bg-white/5 rounded-xl p-4">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-sm text-slate-300">Inner (r {'<'} 0.5 r_max)</span>
                      {io.innerStronger.fdm && <span className="text-xs text-emerald-400 bg-emerald-500/10 px-2 py-0.5 rounded">Stronger</span>}
                    </div>
                    <div className="grid grid-cols-2 gap-2">
                      <StatBox label="r" value={io.innerFDM.r.toFixed(4)} color="text-violet-400" />
                      <StatBox label="n" value={String(io.innerFDM.n)} />
                    </div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-sm text-slate-300">Outer (r {'≥'} 0.5 r_max)</span>
                      {!io.innerStronger.fdm && <span className="text-xs text-emerald-400 bg-emerald-500/10 px-2 py-0.5 rounded">Stronger</span>}
                    </div>
                    <div className="grid grid-cols-2 gap-2">
                      <StatBox label="r" value={io.outerFDM.r.toFixed(4)} color="text-violet-400" />
                      <StatBox label="n" value={String(io.outerFDM.n)} />
                    </div>
                  </div>
                </div>
              </GlassCard>
            </div>

            <GlassCard glow="cyan">
              <h4 className="text-sm font-bold text-white mb-3">Interpretation</h4>
              <div className="text-sm text-slate-300 space-y-2">
                <p>
                  {io.innerStronger.fdm
                    ? 'Inner f_DM coupling is stronger (|r| = ' + Math.abs(io.innerFDM.r).toFixed(3) + ' vs ' + Math.abs(io.outerFDM.r).toFixed(3) + '), consistent with baryonic surface density directly shaping the inner halo profile where baryon concentration is highest.'
                    : 'Outer f_DM coupling is stronger, suggesting extended halo properties are more responsive to baryonic surface density.'}
                </p>
                <p>
                  {io.innerStronger.vdm
                    ? 'Inner V_DM raw correlation is slightly stronger (|r| = ' + Math.abs(io.innerVDM.r).toFixed(3) + ' vs ' + Math.abs(io.outerVDM.r).toFixed(3) + '), but both regions show significant coupling.'
                    : 'Both regions show significant V_DM coupling with comparable strength.'}
                </p>
                <p className="text-cyan-400 font-medium">
                  This radial gradient in coupling strength is a key prediction of baryon-halo coupling models and is not expected from pure feedback scenarios.
                </p>
              </div>
            </GlassCard>
          </div>
        )}

        {activePhase === 3 && (
          <div className="space-y-6">
            <SectionHeader icon={GitBranch} title="Phase 3: Residual Tests" subtitle="After removing Vmax dependence, does Sigma_bar still explain halo properties?" />

            <p className="text-sm text-slate-400">
              We regress each halo property on log(Vmax) alone, then test whether the residuals still correlate with log(Sigma_bar).
              A significant residual correlation proves Sigma_bar carries independent information about halo structure beyond what galaxy mass provides.
            </p>

            <div className="grid grid-cols-3 gap-6">
              {[
                { key: 'vdm' as const, label: 'Delta V_DM', color: 'cyan' },
                { key: 'rDMdom' as const, label: 'Delta r_DMdom', color: 'emerald' },
                { key: 'alpha' as const, label: 'Delta alpha', color: 'violet' },
              ].map((item) => {
                const test = rt[item.key];
                return (
                  <GlassCard key={item.key} glow={test.sigBarExplains ? 'emerald' : 'amber'}>
                    <div className="flex items-center justify-between mb-3">
                      <span className="text-sm font-bold text-white">{item.label}</span>
                      <Badge pass={test.sigBarExplains} label={test.sigBarExplains ? 'Sigma explains' : 'Weak signal'} />
                    </div>
                    <div className="grid grid-cols-2 gap-3 mb-4">
                      <StatBox label="Baseline R² (Vmax)" value={test.baselineR2.toFixed(4)} sub="Vmax alone" color="text-slate-400" />
                      <StatBox label="Residual r" value={test.residualR.toFixed(4)} sub={`slope: ${test.residualSlope.toFixed(4)}`} color={test.sigBarExplains ? 'text-emerald-400' : 'text-amber-400'} />
                    </div>
                    <CouplingScatterPlot
                      data={pd.residuals[item.key]}
                      xKey="logSigBar"
                      yKey="residual"
                      xLabel="log(Sigma_bar)"
                      yLabel={item.label}
                      reg={{ slope: test.residualSlope, intercept: test.residualIntercept, r: test.residualR }}
                      color={item.color === 'cyan' ? '#06b6d4' : item.color === 'emerald' ? '#10b981' : '#8b5cf6'}
                    />
                    <p className="text-xs text-slate-500 mt-2 text-center">
                      n = {test.n} | residual R² = {test.residualR2.toFixed(4)}
                    </p>
                  </GlassCard>
                );
              })}
            </div>

            <GlassCard glow="emerald">
              <h4 className="text-sm font-bold text-white mb-3">Residual Test Summary</h4>
              <div className="text-sm text-slate-300 space-y-2">
                <p>
                  <strong className="text-white">V_DM residuals:</strong> After removing the Vmax trend (R² = {rt.vdm.baselineR2.toFixed(3)} baseline),
                  Sigma_bar still explains residual variation with r = {rt.vdm.residualR.toFixed(3)}.
                  {rt.vdm.sigBarExplains ? ' This confirms Sigma_bar carries independent structural information.' : ' The signal is weak.'}
                </p>
                <p>
                  <strong className="text-white">r_DMdom residuals:</strong> Vmax explains little of r_DMdom variance (R² = {rt.rDMdom.baselineR2.toFixed(3)}),
                  and Sigma_bar independently predicts the DM dominance radius with r = {rt.rDMdom.residualR.toFixed(3)}.
                </p>
                <p>
                  <strong className="text-white">alpha residuals:</strong> Inner slope shows {rt.alpha.sigBarExplains ? 'significant' : 'weak'} residual coupling
                  (r = {rt.alpha.residualR.toFixed(3)}).
                </p>
                <p className="text-emerald-400 font-medium mt-3">
                  {Object.values(rt).filter(r => r.sigBarExplains).length >= 2
                    ? 'At least two halo properties show residual coupling to Sigma_bar after removing mass dependence. This is the strongest evidence that baryonic surface density directly tracks halo structure.'
                    : 'Residual coupling is present but not uniformly strong across all halo properties.'}
                </p>
              </div>
            </GlassCard>
          </div>
        )}
      </div>
    </Layout>
  );
}
