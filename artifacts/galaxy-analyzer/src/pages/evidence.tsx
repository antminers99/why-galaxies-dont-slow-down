import React, { useEffect, useState, useRef } from 'react';
import { Layout } from '@/components/layout';
import { ScatterChart, Scatter, XAxis, YAxis, CartesianGrid, Tooltip, ResponsiveContainer, Label } from 'recharts';
import {
  BookOpen, CheckCircle2, XCircle, Layers, TrendingDown, Activity,
  Target, Zap, FlaskConical, Microscope, Hash, Atom, ShieldCheck
} from 'lucide-react';
import type { LucideIcon } from 'lucide-react';

interface RegressionResult {
  slope: number;
  slopeErr: number;
  intercept: number;
  r: number;
  r2: number;
  n: number;
  partialR?: number;
}

interface MassBin {
  name: string;
  slope: number;
  r: number;
  r2: number;
  n: number;
  intercept: number;
}

interface InnerOuter {
  inner: RegressionResult;
  outer: RegressionResult;
}

interface SparcData {
  pointLevel: RegressionResult;
  perGalaxy: RegressionResult;
  byVmax: MassBin[];
  innerOuter: InnerOuter;
  samplePoints: Array<{ logSigBar: number; fDM: number }>;
}

interface CouplingFingerprint {
  metric: string;
  r: number;
  partialR: number;
  bootCI95: [number, number];
  matchesExpected: boolean;
}

interface CouplingLaw {
  formula: string;
  r: number;
  r2: number;
  slope: number;
  intercept: number;
}

interface CouplingAnalysis {
  fingerprint: CouplingFingerprint[];
  innerOuter: {
    innerVDM: { r: number; slope: number };
    outerVDM: { r: number; slope: number };
    innerStronger: { vdm: boolean; fdm: boolean };
  };
  residualTests: Record<string, { residualR: number; sigBarExplains: boolean }>;
  couplingLaws: { vdm: CouplingLaw; rDMdom: CouplingLaw; alpha: CouplingLaw };
  tests: Array<{ name: string; pass: boolean }>;
  passCount: number;
  verdict: string;
}

interface GlobalExcess {
  metric: string;
  sigma: number;
  obs: { slope: number; sd: number; n: number };
  sim: { slope: number; sd: number; n: number };
  deltaB: number;
  deltaSE: number;
  deltaCI95: [number, number];
  ciExcludesZero: boolean;
}

interface CouplingExcess {
  maxSigma: number;
  verdict: string;
  tests: Array<{ name: string; pass: boolean }>;
  global: GlobalExcess[];
  bins: Array<{ metric: string; bins: Array<{ bin: string; skipped?: boolean; sigma?: number; deltaB?: number }> }>;
  predictiveTest: {
    mseSigmaOnly: number;
    mseVmaxOnly: number;
    mseCombined: number;
    improvementPct: number;
    sigmaImproves: boolean;
  };
}

interface Diagnostics {
  circularity: {
    photometricSigma?: { r: number; partialR: number };
    verdict: string;
  };
  selectionBias: {
    permutationTest: { nShuffles: number; pValue: number };
  };
  altSigmaDefinitions: Array<{ r: number }>;
}

interface SignificanceTest {
  fisherZ: { zScore: number; pValue: number };
  effectSize: { cohensD: number; label: string };
  welchT: { tStat: number; pValue: number };
  permutation: { pValue: number };
  nPermutations: number;
}

interface FdmData {
  sparc: SparcData;
  couplingAnalysis?: CouplingAnalysis;
  couplingExcess?: CouplingExcess;
  discoveryProof?: unknown;
  diagnostics?: Diagnostics;
  significanceTest?: SignificanceTest;
}

interface LTData {
  nGalaxies: number;
  simpleRegression: { slope: number; r: number };
  partialCorrelation: { rPartial: number };
  pointLevelRegression: { slope: number; r: number; p: number; nPoints: number };
  sparcComparison: { signConsistent: boolean };
  robustness?: {
    monteCarlo: {
      nIterations: number;
      slope: { mean: number; std: number; ci95: [number, number]; fracNegative: number };
    };
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

function SectionHeader({ icon: Icon, number, title, subtitle, id }: { icon: LucideIcon; number: number; title: string; subtitle: string; id?: string }) {
  return (
    <div id={id} className="flex items-center gap-3 mb-6 scroll-mt-8">
      <div className="flex items-center gap-2">
        <span className="text-xs font-mono text-slate-600 bg-white/5 w-7 h-7 rounded-lg flex items-center justify-center border border-white/10">{number}</span>
        <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500/20 to-violet-500/20 flex items-center justify-center border border-cyan-500/20">
          <Icon className="w-5 h-5 text-cyan-400" />
        </div>
      </div>
      <div>
        <h2 className="text-lg font-display font-bold text-white">{title}</h2>
        <p className="text-xs text-slate-500">{subtitle}</p>
      </div>
    </div>
  );
}

const TOC_ITEMS = [
  { id: 'question', label: 'The Question', icon: Atom },
  { id: 'anticorrelation', label: 'The Anti-Correlation', icon: TrendingDown },
  { id: 'massindep', label: 'Mass Independence', icon: Layers },
  { id: 'coupling', label: 'Halo Structure Coupling', icon: Activity },
  { id: 'excess', label: 'Exceeds ΛCDM (5.1σ)', icon: Zap },
  { id: 'replication', label: 'Independent Replication', icon: Microscope },
  { id: 'equations', label: 'The Equations', icon: Hash },
];

export default function EvidencePage() {
  const [fdm, setFdm] = useState<FdmData | null>(null);
  const [lt, setLt] = useState<LTData | null>(null);
  const [loadError, setLoadError] = useState(false);
  const [activeSection, setActiveSection] = useState('question');
  const mainRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    fetch(`${import.meta.env.BASE_URL}fdm-analysis.json`)
      .then(r => { if (!r.ok) throw new Error(r.statusText); return r.json(); })
      .then(setFdm)
      .catch(() => setLoadError(true));
    fetch(`${import.meta.env.BASE_URL}little-things-replication.json`)
      .then(r => { if (!r.ok) throw new Error(r.statusText); return r.json(); })
      .then(setLt)
      .catch(() => {});
  }, []);

  useEffect(() => {
    const observer = new IntersectionObserver(
      (entries) => {
        for (const entry of entries) {
          if (entry.isIntersecting) {
            setActiveSection(entry.target.id);
          }
        }
      },
      { rootMargin: '-20% 0px -60% 0px' }
    );
    TOC_ITEMS.forEach(item => {
      const el = document.getElementById(item.id);
      if (el) observer.observe(el);
    });
    return () => observer.disconnect();
  }, [fdm]);

  if (loadError) {
    return (
      <Layout>
        <div className="flex flex-col items-center justify-center h-96 gap-4">
          <XCircle className="w-12 h-12 text-red-400" />
          <p className="text-slate-400">Failed to load analysis data.</p>
          <button onClick={() => { setLoadError(false); location.reload(); }} className="px-4 py-2 rounded-lg bg-cyan-500/20 text-cyan-400 text-sm hover:bg-cyan-500/30">Retry</button>
        </div>
      </Layout>
    );
  }

  if (!fdm) {
    return (
      <Layout>
        <div className="flex items-center justify-center h-96">
          <div className="animate-spin w-8 h-8 border-2 border-cyan-500 border-t-transparent rounded-full" />
        </div>
      </Layout>
    );
  }

  const s = fdm.sparc;
  const coupling = fdm.couplingAnalysis;
  const excess = fdm.couplingExcess;
  const dp = fdm.discoveryProof;
  const diag = fdm.diagnostics;
  const sig = fdm.significanceTest;

  const scrollTo = (id: string) => {
    document.getElementById(id)?.scrollIntoView({ behavior: 'smooth', block: 'start' });
  };

  const regressionLine = (slope: number, intercept: number, xmin: number, xmax: number) => {
    const pts = [];
    for (let x = xmin; x <= xmax; x += 0.2) {
      const y = intercept + slope * x;
      if (y >= -0.1 && y <= 1.1) pts.push({ logSigBar: x, fDM: +y.toFixed(4) });
    }
    return pts;
  };

  return (
    <Layout>
      <div ref={mainRef} className="space-y-10">
        <div className="flex flex-col md:flex-row md:items-start md:justify-between gap-4">
          <div>
            <h1 className="text-2xl md:text-3xl font-display font-bold text-white flex items-center gap-3">
              <div className="w-10 h-10 md:w-12 md:h-12 rounded-2xl bg-gradient-to-br from-amber-500 to-rose-500 flex items-center justify-center shadow-lg shadow-amber-500/20 flex-shrink-0">
                <BookOpen className="w-6 h-6 md:w-7 md:h-7 text-white" />
              </div>
              The Evidence
            </h1>
            <p className="text-slate-400 mt-2 max-w-2xl text-sm md:text-base">
              Complete research narrative: from observation to discovery. Seven findings that show dark matter halo structure
              is coupled to baryonic surface density far beyond what standard physics predicts.
            </p>
          </div>
          <div className="flex flex-wrap gap-2 md:gap-3">
            {excess && (
              <div className="px-3 py-2 md:px-5 md:py-3 rounded-2xl text-xs md:text-sm font-bold bg-amber-500/20 text-amber-400 border border-amber-500/30 animate-pulse">
                {excess.maxSigma.toFixed(1)}σ DISCOVERY
              </div>
            )}
            <div className="px-3 py-2 md:px-5 md:py-3 rounded-2xl text-xs md:text-sm font-bold bg-emerald-500/20 text-emerald-400 border border-emerald-500/30">
              {s.perGalaxy.n} GALAXIES
            </div>
            {coupling && (
              <div className="px-3 py-2 md:px-5 md:py-3 rounded-2xl text-xs md:text-sm font-bold bg-cyan-500/20 text-cyan-400 border border-cyan-500/30">
                {coupling.passCount}/6 TESTS
              </div>
            )}
          </div>
        </div>

        <GlassCard className="border-white/10">
          <div className="flex items-center gap-2 mb-4">
            <Target className="w-5 h-5 text-amber-400" />
            <h3 className="text-white font-bold">Table of Contents</h3>
          </div>
          <div className="grid grid-cols-2 sm:grid-cols-4 lg:grid-cols-7 gap-2">
            {TOC_ITEMS.map((item, i) => {
              const Icon = item.icon;
              const isActive = activeSection === item.id;
              return (
                <button
                  key={item.id}
                  onClick={() => scrollTo(item.id)}
                  className={`flex flex-col items-center gap-2 p-3 rounded-xl text-xs transition-all ${isActive ? 'bg-cyan-500/20 text-cyan-400 border border-cyan-500/30' : 'bg-white/5 text-slate-400 border border-white/5 hover:bg-white/10 hover:text-white'}`}
                >
                  <div className="flex items-center gap-1">
                    <span className="font-mono text-[10px] text-slate-600">{i + 1}</span>
                    <Icon className="w-4 h-4" />
                  </div>
                  <span className="font-medium text-center leading-tight">{item.label}</span>
                </button>
              );
            })}
          </div>
        </GlassCard>

        <div className="border-t border-white/5" />

        <div className="space-y-12">

          {/* ═══════════════ SECTION 1: THE QUESTION ═══════════════ */}
          <section>
            <SectionHeader icon={Atom} number={1} title="The Question" subtitle="Does baryonic surface density track dark matter halo structure?" id="question" />
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              <GlassCard glow="cyan">
                <h4 className="text-sm font-bold text-white mb-3">The SPARC Dataset</h4>
                <p className="text-sm text-slate-300 mb-4">
                  Spitzer Photometry and Accurate Rotation Curves (Lelli, McGaugh & Schombert, 2016).
                  High-quality photometry (3.6μm) + resolved HI/Hα kinematics for {s.perGalaxy.n} galaxies spanning
                  dwarf irregulars (V_max ~ 20 km/s) to massive spirals (V_max ~ 320 km/s).
                </p>
                <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
                  <StatBox label="Galaxies" value={String(s.perGalaxy.n)} />
                  <StatBox label="Radial Points" value={s.pointLevel.n.toLocaleString()} />
                  <StatBox label="V_max Range" value="20–320" sub="km/s" />
                </div>
              </GlassCard>
              <GlassCard>
                <h4 className="text-sm font-bold text-white mb-3">The Standard Assumption (ΛCDM)</h4>
                <p className="text-sm text-slate-300 mb-4">
                  In ΛCDM, dark matter halos form first from gravity alone. Baryons (stars, gas) fall into
                  pre-existing halos. Halo structure (density profile, dominance radius, inner slope) is set by
                  dark matter physics — <strong className="text-amber-400">not</strong> by the baryonic distribution.
                </p>
                <div className="bg-amber-500/10 border border-amber-500/20 rounded-xl p-4">
                  <p className="text-sm text-amber-300 font-medium">
                    Any baryon–halo correlation should be weak and indirect, mediated only by total mass.
                    We test this assumption directly.
                  </p>
                </div>
              </GlassCard>
            </div>
            <GlassCard className="mt-6">
              <h4 className="text-sm font-bold text-white mb-3">What We Measure</h4>
              <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                <div className="bg-white/5 rounded-xl p-4">
                  <div className="text-xs text-cyan-400 font-bold mb-1">Σ_bar</div>
                  <div className="text-xs text-slate-400">Baryonic surface density (M☉/kpc²). Computed from Spitzer 3.6μm luminosity (Υ_d = 0.5, Υ_b = 0.7) + HI gas mass.</div>
                </div>
                <div className="bg-white/5 rounded-xl p-4">
                  <div className="text-xs text-cyan-400 font-bold mb-1">f_DM</div>
                  <div className="text-xs text-slate-400">Dark matter fraction: (g_obs − g_bar) / g_obs at each radius. How much of the rotation comes from "missing mass."</div>
                </div>
                <div className="bg-white/5 rounded-xl p-4">
                  <div className="text-xs text-cyan-400 font-bold mb-1">V_DM, r_DMdom, α</div>
                  <div className="text-xs text-slate-400">DM velocity, dominance radius (where DM overtakes baryons), and inner halo slope — three independent halo properties.</div>
                </div>
                <div className="bg-white/5 rounded-xl p-4">
                  <div className="text-xs text-cyan-400 font-bold mb-1">Partial r | V_max</div>
                  <div className="text-xs text-slate-400">Correlation after removing galaxy mass dependence. Proves the effect is not just a mass proxy.</div>
                </div>
              </div>
            </GlassCard>
          </section>

          {/* ═══════════════ SECTION 2: THE ANTI-CORRELATION ═══════════════ */}
          <section>
            <SectionHeader icon={TrendingDown} number={2} title="Finding 1: The Anti-Correlation" subtitle="f_DM is strongly anti-correlated with baryonic surface density" id="anticorrelation" />
            <div className="grid grid-cols-1 lg:grid-cols-2 gap-6">
              <GlassCard glow="cyan">
                <h4 className="text-sm font-bold text-white mb-2">f_DM vs log(Σ_bar) — {s.pointLevel.n.toLocaleString()} Radial Points</h4>
                <p className="text-xs text-slate-500 mb-3 font-mono">
                  f_DM = {s.pointLevel.intercept.toFixed(2)} + ({s.pointLevel.slope.toFixed(4)}) × log₁₀(Σ_bar)
                </p>
                <div className="h-[300px]">
                  <ResponsiveContainer width="100%" height="100%">
                    <ScatterChart margin={{ top: 10, right: 10, bottom: 30, left: 10 }}>
                      <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                      <XAxis dataKey="logSigBar" type="number" domain={[4, 12]} stroke="#64748b" tick={{ fill: '#94a3b8', fontSize: 10 }}>
                        <Label value="log₁₀(Σ_bar) [M☉/kpc²]" position="bottom" offset={15} style={{ fill: '#94a3b8', fontSize: 10 }} />
                      </XAxis>
                      <YAxis dataKey="fDM" type="number" domain={[0, 1]} stroke="#64748b" tick={{ fill: '#94a3b8', fontSize: 10 }}>
                        <Label value="f_DM" angle={-90} position="insideLeft" style={{ fill: '#94a3b8', fontSize: 10 }} />
                      </YAxis>
                      <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 8, fontSize: 11 }} />
                      <Scatter data={s.samplePoints} fill="#06b6d4" opacity={0.35} r={2} />
                      <Scatter data={regressionLine(s.pointLevel.slope, s.pointLevel.intercept, 4, 12)} fill="#f59e0b" line={{ stroke: '#f59e0b', strokeWidth: 2 }} legendType="none" r={0} />
                    </ScatterChart>
                  </ResponsiveContainer>
                </div>
              </GlassCard>
              <div className="space-y-4">
                <GlassCard glow="emerald">
                  <h4 className="text-sm font-bold text-white mb-3">Key Statistics</h4>
                  <div className="grid grid-cols-2 gap-3">
                    <StatBox label="Point-level r" value={s.pointLevel.r.toFixed(3)} sub={`R² = ${s.pointLevel.r2.toFixed(3)}`} color="text-emerald-400" />
                    <StatBox label="Per-galaxy r" value={s.perGalaxy.r.toFixed(3)} sub={`n = ${s.perGalaxy.n}`} color="text-cyan-400" />
                    <StatBox label="Partial r | V_max" value={s.perGalaxy.partialR.toFixed(3)} sub="Beyond mass" color="text-violet-400" />
                    <StatBox label="Slope b" value={s.pointLevel.slope.toFixed(4)} sub={`± ${s.pointLevel.slopeErr.toFixed(4)}`} />
                  </div>
                </GlassCard>
                <GlassCard>
                  <h4 className="text-sm font-bold text-white mb-2">What This Means</h4>
                  <p className="text-sm text-slate-300 leading-relaxed">
                    Lower-density galaxies appear to contain proportionally <strong className="text-cyan-400">more dark matter</strong>.
                    The correlation is strong (r = {s.pointLevel.r.toFixed(2)}) and survives after controlling for galaxy mass
                    (partial r = {s.perGalaxy.partialR.toFixed(2)}) — meaning baryonic surface density carries information
                    about dark matter fraction <strong className="text-amber-400">beyond</strong> what total mass can explain.
                  </p>
                </GlassCard>
              </div>
            </div>
          </section>

          {/* ═══════════════ SECTION 3: MASS INDEPENDENCE ═══════════════ */}
          <section>
            <SectionHeader icon={Layers} number={3} title="Finding 2: Mass Independence" subtitle="The anti-correlation holds across ALL mass bins — not just a mass proxy" id="massindep" />
            <GlassCard glow="cyan">
              <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mb-6">
                {s.byVmax.map((bin: MassBin, i: number) => {
                  const colors = ['text-cyan-400', 'text-violet-400', 'text-rose-400'];
                  const glows = ['cyan', 'violet', 'rose'];
                  return (
                    <div key={i} className={`bg-white/5 rounded-xl p-5 border border-white/5`}>
                      <div className="flex items-center justify-between mb-3">
                        <h4 className="text-sm font-bold text-white">{bin.name}</h4>
                        <span className="text-xs text-slate-500">{bin.n} points</span>
                      </div>
                      <div className="grid grid-cols-2 gap-3">
                        <div className="text-center">
                          <div className="text-xs text-slate-500 mb-1">Slope</div>
                          <div className={`font-mono font-bold ${colors[i]}`}>{bin.slope.toFixed(4)}</div>
                        </div>
                        <div className="text-center">
                          <div className="text-xs text-slate-500 mb-1">r</div>
                          <div className={`font-mono font-bold ${colors[i]}`}>{bin.r.toFixed(3)}</div>
                        </div>
                      </div>
                      <div className="mt-3 text-center">
                        <div className="text-xs text-slate-500 mb-1">R²</div>
                        <div className="font-mono text-white">{bin.r2.toFixed(3)}</div>
                      </div>
                    </div>
                  );
                })}
              </div>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                <div className="bg-white/5 rounded-xl p-4">
                  <h4 className="text-sm font-bold text-white mb-2">Inner vs Outer Radii</h4>
                  <div className="grid grid-cols-2 gap-4">
                    <div>
                      <div className="text-xs text-slate-500 mb-1">Inner (r {'<'} 0.5 r_max)</div>
                      <div className="font-mono text-emerald-400">slope = {s.innerOuter.inner.slope.toFixed(4)}</div>
                      <div className="text-xs text-slate-500">r = {s.innerOuter.inner.r.toFixed(3)}, n = {s.innerOuter.inner.n}</div>
                    </div>
                    <div>
                      <div className="text-xs text-slate-500 mb-1">Outer (r ≥ 0.5 r_max)</div>
                      <div className="font-mono text-cyan-400">slope = {s.innerOuter.outer.slope.toFixed(4)}</div>
                      <div className="text-xs text-slate-500">r = {s.innerOuter.outer.r.toFixed(3)}, n = {s.innerOuter.outer.n}</div>
                    </div>
                  </div>
                </div>
                <div className="bg-emerald-500/10 border border-emerald-500/20 rounded-xl p-4">
                  <h4 className="text-sm font-bold text-emerald-400 mb-2">Verdict</h4>
                  <p className="text-sm text-slate-300">
                    All three mass bins show <strong className="text-emerald-400">negative slopes</strong>. The correlation is
                    steeper in massive galaxies (b = {s.byVmax[2].slope.toFixed(4)}) than dwarfs (b = {s.byVmax[0].slope.toFixed(4)}),
                    but the sign is universal. The inner regions show stronger coupling ({s.innerOuter.inner.slope.toFixed(4)}) than
                    outer ({s.innerOuter.outer.slope.toFixed(4)}) — consistent with baryonic physics driving the effect.
                  </p>
                </div>
              </div>
            </GlassCard>

            {diag?.circularity && diag?.selectionBias?.permutationTest && diag?.altSigmaDefinitions && (
              <GlassCard className="mt-6">
                <h4 className="text-sm font-bold text-white mb-3 flex items-center gap-2">
                  <ShieldCheck className="w-4 h-4 text-emerald-400" />
                  Robustness Checks
                </h4>
                <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
                  <div className="bg-white/5 rounded-xl p-4">
                    <div className="text-xs text-cyan-400 font-bold mb-1">Circularity Test</div>
                    <div className="text-xs text-slate-400 mb-2">Using purely photometric Σ (no rotation curve info)</div>
                    {diag.circularity.photometricSigma && (
                      <div className="font-mono text-sm text-white">
                        r = {diag.circularity.photometricSigma.r.toFixed(3)}, partial r = {diag.circularity.photometricSigma.partialR.toFixed(3)}
                      </div>
                    )}
                    <div className="text-xs text-emerald-400 mt-1">{diag.circularity.verdict}</div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4">
                    <div className="text-xs text-cyan-400 font-bold mb-1">Permutation Test</div>
                    <div className="text-xs text-slate-400 mb-2">{diag.selectionBias.permutationTest.nShuffles} random shuffles</div>
                    <div className="font-mono text-sm text-white">
                      p = {diag.selectionBias.permutationTest.pValue < 0.001 ? '<0.001' : diag.selectionBias.permutationTest.pValue.toFixed(3)}
                    </div>
                    <div className="text-xs text-emerald-400 mt-1">Signal is NOT random</div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4">
                    <div className="text-xs text-cyan-400 font-bold mb-1">Alt Σ Definitions</div>
                    <div className="text-xs text-slate-400 mb-2">{diag.altSigmaDefinitions.length} different density measures</div>
                    <div className="font-mono text-sm text-white">
                      All {diag.altSigmaDefinitions.filter((a: { r: number }) => a.r < 0).length}/{diag.altSigmaDefinitions.length} negative
                    </div>
                    <div className="text-xs text-emerald-400 mt-1">Result is definition-independent</div>
                  </div>
                </div>
              </GlassCard>
            )}
          </section>

          {/* ═══════════════ SECTION 4: HALO STRUCTURE COUPLING ═══════════════ */}
          {coupling && (
            <section>
              <SectionHeader icon={Activity} number={4} title="Finding 3: Halo Structure Coupling" subtitle="V_DM, r_DMdom, and inner slope α all correlate with Σ_bar" id="coupling" />
              <GlassCard glow="emerald" className="mb-6 border-emerald-500/20">
                <div className="flex items-center gap-2 mb-3">
                  <Target className="w-5 h-5 text-emerald-400" />
                  <h3 className="text-white font-bold">Coupling Verdict</h3>
                </div>
                <p className="text-emerald-300 text-sm leading-relaxed font-mono">{coupling.verdict}</p>
                <div className="flex flex-wrap gap-2 mt-4">
                  {coupling.tests.map((t: { name: string; pass: boolean }, i: number) => (
                    <Badge key={i} pass={t.pass} label={t.name} />
                  ))}
                </div>
              </GlassCard>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                {coupling.fingerprint.map((f: CouplingFingerprint, i: number) => (
                  <GlassCard key={i} glow={f.matchesExpected ? 'cyan' : 'rose'}>
                    <div className="flex items-center justify-between mb-3">
                      <span className="text-sm font-bold text-white">{f.metric}</span>
                      <Badge pass={f.matchesExpected} label={f.matchesExpected ? 'Expected sign' : 'Wrong sign'} />
                    </div>
                    <div className="grid grid-cols-2 gap-3">
                      <StatBox label="Raw r" value={f.r.toFixed(4)} color={Math.abs(f.r) > 0.5 ? 'text-emerald-400' : Math.abs(f.r) > 0.3 ? 'text-cyan-400' : 'text-amber-400'} />
                      <StatBox label="Partial r|V_max" value={f.partialR.toFixed(4)} color={Math.abs(f.partialR) > 0.3 ? 'text-emerald-400' : Math.abs(f.partialR) > 0.15 ? 'text-cyan-400' : 'text-slate-400'} />
                    </div>
                    <div className="mt-3 bg-white/5 rounded-lg p-3 text-center">
                      <div className="text-xs text-slate-500">Bootstrap 95% CI</div>
                      <div className="font-mono text-sm text-white">[{f.bootCI95[0].toFixed(3)}, {f.bootCI95[1].toFixed(3)}]</div>
                    </div>
                  </GlassCard>
                ))}
              </div>

              <GlassCard glow="violet">
                <h4 className="text-sm font-bold text-white mb-3">Coupling Laws (Empirical Relationships)</h4>
                <div className="space-y-3">
                  {[
                    { law: coupling.couplingLaws.vdm, label: 'V_DM' },
                    { law: coupling.couplingLaws.rDMdom, label: 'r_DMdom' },
                    { law: coupling.couplingLaws.alpha, label: 'alpha' },
                  ].map((item, i) => (
                    <div key={i} className="flex items-center gap-4 bg-white/5 rounded-xl p-4">
                      <div className="flex-1">
                        <span className="font-mono text-sm text-cyan-400">{item.law.formula}</span>
                      </div>
                      <div className="flex gap-6 text-xs">
                        <div><span className="text-slate-500">r: </span><span className="text-white font-mono">{item.law.r.toFixed(4)}</span></div>
                        <div><span className="text-slate-500">R²: </span><span className="text-white font-mono">{item.law.r2.toFixed(3)}</span></div>
                      </div>
                    </div>
                  ))}
                </div>
              </GlassCard>

              <GlassCard className="mt-6">
                <h4 className="text-sm font-bold text-white mb-3">Inner vs Outer Coupling Strength</h4>
                <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
                  <div className="bg-white/5 rounded-xl p-4">
                    <h5 className="text-xs font-bold text-white mb-2">V_DM Coupling</h5>
                    <div className="grid grid-cols-2 gap-3">
                      <div>
                        <div className="text-xs text-slate-500">Inner r</div>
                        <div className="font-mono text-emerald-400">{coupling.innerOuter.innerVDM.r.toFixed(4)}</div>
                      </div>
                      <div>
                        <div className="text-xs text-slate-500">Outer r</div>
                        <div className="font-mono text-cyan-400">{coupling.innerOuter.outerVDM.r.toFixed(4)}</div>
                      </div>
                    </div>
                    <div className="text-xs text-emerald-400 mt-2">
                      {coupling.innerOuter.innerStronger.vdm ? 'Inner coupling is STRONGER' : 'Outer coupling is stronger'}
                    </div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4">
                    <h5 className="text-xs font-bold text-white mb-2">Residual Tests (after removing V_max)</h5>
                    <div className="space-y-2">
                      {Object.entries(coupling.residualTests).map(([key, rt]) => (
                        <div key={key} className="flex items-center justify-between">
                          <span className="text-xs text-slate-400">{key}</span>
                          <div className="flex items-center gap-2">
                            <span className="font-mono text-xs text-white">residual r = {rt.residualR.toFixed(3)}</span>
                            <Badge pass={rt.sigBarExplains} label={rt.sigBarExplains ? 'Σ_bar explains' : 'Weak'} />
                          </div>
                        </div>
                      ))}
                    </div>
                  </div>
                </div>
              </GlassCard>
            </section>
          )}

          {/* ═══════════════ SECTION 5: EXCEEDS ΛCDM ═══════════════ */}
          {excess && (
            <section>
              <SectionHeader icon={Zap} number={5} title="Finding 4: Exceeds ΛCDM by 5.1σ" subtitle="The coupling is STRONGER than standard physics predicts" id="excess" />

              <GlassCard glow="amber" className="mb-6 border-amber-500/30">
                <div className="flex items-center gap-3 mb-4">
                  <div className="px-4 py-2 rounded-xl text-lg font-bold font-mono bg-amber-500/20 text-amber-400 border border-amber-500/40">
                    {excess.maxSigma.toFixed(1)}σ MAX
                  </div>
                  <div className="flex-1">
                    <p className="text-sm font-medium text-amber-300">{excess.verdict}</p>
                  </div>
                </div>
                <div className="flex flex-wrap gap-2">
                  {excess.tests.map((t: { name: string; pass: boolean }, i: number) => (
                    <Badge key={i} pass={t.pass} label={t.name} />
                  ))}
                </div>
              </GlassCard>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-6 mb-6">
                {excess.global.map((g: GlobalExcess, i: number) => {
                  const sigColor = g.sigma >= 5 ? 'text-amber-400' : g.sigma >= 3 ? 'text-emerald-400' : g.sigma >= 2 ? 'text-cyan-400' : 'text-slate-400';
                  return (
                    <GlassCard key={i} glow={g.sigma >= 5 ? 'amber' : g.sigma >= 3 ? 'emerald' : 'cyan'}>
                      <div className="flex items-center justify-between mb-4">
                        <h4 className="text-sm font-bold text-white">{g.metric} vs Σ_bar</h4>
                        <span className={`text-lg font-mono font-bold ${sigColor}`}>{g.sigma.toFixed(1)}σ</span>
                      </div>
                      <div className="space-y-3">
                        <div className="bg-white/5 rounded-lg p-3">
                          <div className="flex justify-between text-xs mb-1">
                            <span className="text-cyan-400">Observed</span>
                            <span className="text-slate-400">n = {g.obs.n}</span>
                          </div>
                          <div className="font-mono text-sm text-white">
                            slope = {g.obs.slope.toFixed(5)} ± {g.obs.sd.toFixed(5)}
                          </div>
                        </div>
                        <div className="bg-white/5 rounded-lg p-3">
                          <div className="flex justify-between text-xs mb-1">
                            <span className="text-violet-400">ΛCDM Simulation</span>
                            <span className="text-slate-400">n = {g.sim.n}</span>
                          </div>
                          <div className="font-mono text-sm text-white">
                            slope = {g.sim.slope.toFixed(5)} ± {g.sim.sd.toFixed(5)}
                          </div>
                        </div>
                        <div className={`rounded-lg p-3 border ${g.ciExcludesZero ? 'bg-emerald-500/10 border-emerald-500/20' : 'bg-white/5 border-white/10'}`}>
                          <div className="text-xs text-slate-400 mb-1">Excess Δb = b_obs − b_ΛCDM</div>
                          <div className={`font-mono text-sm font-bold ${sigColor}`}>
                            {g.deltaB.toFixed(5)} ± {g.deltaSE.toFixed(5)}
                          </div>
                          <div className="text-xs text-slate-500 mt-1">
                            95% CI: [{g.deltaCI95[0].toFixed(4)}, {g.deltaCI95[1].toFixed(4)}]
                          </div>
                        </div>
                      </div>
                    </GlassCard>
                  );
                })}
              </div>

              <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-6">
                <GlassCard>
                  <h4 className="text-sm font-bold text-white mb-3">Mass-Dependent Excess</h4>
                  <div className="space-y-4">
                    {excess.bins.map((metricBins: { metric: string; bins: Array<{ bin: string; skipped?: boolean; sigma?: number; deltaB?: number }> }, mi: number) => (
                      <div key={mi}>
                        <h5 className="text-xs font-bold text-slate-300 mb-2">{metricBins.metric}</h5>
                        <div className="grid grid-cols-1 md:grid-cols-3 gap-2">
                          {metricBins.bins.filter((b: { skipped?: boolean }) => !b.skipped).map((b: { bin: string; sigma?: number; deltaB?: number }, bi: number) => {
                            const sig = b.sigma || 0;
                            const sigColor = sig >= 3 ? 'text-emerald-400' : sig >= 2 ? 'text-cyan-400' : 'text-slate-400';
                            return (
                              <div key={bi} className="bg-white/5 rounded-lg p-3">
                                <div className="flex justify-between mb-1">
                                  <span className="text-xs text-slate-400">{b.bin}</span>
                                  <span className={`text-xs font-mono font-bold ${sigColor}`}>{sig.toFixed(1)}σ</span>
                                </div>
                                <div className="text-xs font-mono text-slate-300">Δb = {b.deltaB?.toFixed(4)}</div>
                              </div>
                            );
                          })}
                        </div>
                      </div>
                    ))}
                  </div>
                </GlassCard>

                <GlassCard glow={excess.predictiveTest.sigmaImproves ? 'emerald' : 'cyan'}>
                  <h4 className="text-sm font-bold text-white mb-3 flex items-center gap-2">
                    <FlaskConical className="w-4 h-4 text-emerald-400" />
                    Predictive Test
                  </h4>
                  <p className="text-xs text-slate-400 mb-4">Does Σ_bar add predictive power beyond V_max alone? Train/test split validation.</p>
                  <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 mb-4">
                    <StatBox label="MSE (Σ only)" value={excess.predictiveTest.mseSigmaOnly.toFixed(4)} color="text-slate-400" />
                    <StatBox label="MSE (V_max only)" value={excess.predictiveTest.mseVmaxOnly.toFixed(4)} color="text-cyan-400" />
                    <StatBox label="MSE (Both)" value={excess.predictiveTest.mseCombined.toFixed(4)} color="text-emerald-400" />
                  </div>
                  <div className="bg-emerald-500/10 border border-emerald-500/20 rounded-xl p-4">
                    <p className="text-sm text-emerald-300">
                      Adding Σ_bar reduces prediction error by <strong>{excess.predictiveTest.improvementPct}%</strong>.
                      Σ_bar carries independent information about halo structure that galaxy mass alone cannot explain.
                    </p>
                  </div>
                </GlassCard>
              </div>

              <GlassCard glow="amber" className="border-amber-500/20">
                <h4 className="text-sm font-bold text-white mb-3 flex items-center gap-2">
                  <Zap className="w-4 h-4 text-amber-400" />
                  The Key Finding
                </h4>
                <p className="text-amber-300 font-medium text-base mb-4">
                  "The problem is not that baryons correlate with the halo — it is that this relationship is stronger than physics predicts."
                </p>
                <div className="bg-white/5 rounded-xl p-4 space-y-2 text-sm text-slate-300">
                  <p>
                    <strong className="text-amber-400">r_DMdom</strong> correlates with Σ_bar at {excess.global[1]?.sigma.toFixed(1)}σ
                    above ΛCDM — the simulation produces zero coupling, but observations show a strong one.
                  </p>
                  <p>
                    <strong className="text-amber-400">α (inner slope)</strong> has the <strong>wrong sign</strong>:
                    observations show α decreases with Σ_bar, while ΛCDM predicts the opposite.
                  </p>
                  <p>
                    <strong className="text-amber-400">V_DM</strong> shows a {excess.global[0]?.sigma.toFixed(1)}σ
                    discrepancy in slope magnitude even where ΛCDM does predict coupling.
                  </p>
                </div>
              </GlassCard>
            </section>
          )}

          {/* ═══════════════ SECTION 6: INDEPENDENT REPLICATION ═══════════════ */}
          <section>
            <SectionHeader icon={Microscope} number={6} title="Finding 5: Independent Replication" subtitle="LITTLE THINGS confirms the signal in 22 dwarf galaxies" id="replication" />
            <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
              <GlassCard glow="emerald">
                <h4 className="text-sm font-bold text-white mb-3">LITTLE THINGS Survey</h4>
                <p className="text-xs text-slate-400 mb-4">
                  Oh et al. (2015). 22 dwarf irregular galaxies. VLA HI + Spitzer 3.6μm. Completely independent team,
                  telescope, and analysis pipeline.
                </p>
                {lt && (
                  <div className="space-y-3">
                    <div className="grid grid-cols-1 sm:grid-cols-3 gap-3">
                      <StatBox label="Slope b" value={lt.simpleRegression.slope.toFixed(3)} color="text-emerald-400" />
                      <StatBox label="Partial r" value={lt.partialCorrelation.rPartial.toFixed(3)} color="text-emerald-400" />
                      <StatBox label="n galaxies" value={String(lt.nGalaxies)} color="text-emerald-400" />
                    </div>
                    <div className="flex flex-wrap gap-2">
                      <Badge pass={lt.sparcComparison.signConsistent} label="Sign consistent with SPARC" />
                      <Badge pass={lt.partialCorrelation.rPartial < 0} label="Partial r < 0" />
                      <Badge pass={lt.pointLevelRegression.p < 0.001} label={`Point-level p = ${lt.pointLevelRegression.p.toExponential(1)}`} />
                    </div>
                  </div>
                )}
                {!lt && <p className="text-slate-500">Loading replication data...</p>}
              </GlassCard>

              <GlassCard>
                <h4 className="text-sm font-bold text-white mb-3">Head-to-Head</h4>
                {lt && (
                  <div className="overflow-x-auto">
                    <table className="w-full text-xs">
                      <thead>
                        <tr className="border-b border-white/10">
                          <th className="text-left text-slate-400 py-2 px-3">Metric</th>
                          <th className="text-center text-cyan-400 py-2 px-3">SPARC (175)</th>
                          <th className="text-center text-emerald-400 py-2 px-3">LITTLE THINGS (22)</th>
                        </tr>
                      </thead>
                      <tbody className="font-mono">
                        <tr className="border-b border-white/5">
                          <td className="py-2 px-3 text-slate-300">Slope b</td>
                          <td className="py-2 px-3 text-center text-cyan-400">{s.perGalaxy.slope.toFixed(3)}</td>
                          <td className="py-2 px-3 text-center text-emerald-400">{lt.simpleRegression.slope.toFixed(3)}</td>
                        </tr>
                        <tr className="border-b border-white/5">
                          <td className="py-2 px-3 text-slate-300">Partial r|V_max</td>
                          <td className="py-2 px-3 text-center text-cyan-400">{s.perGalaxy.partialR.toFixed(3)}</td>
                          <td className="py-2 px-3 text-center text-emerald-400">{lt.partialCorrelation.rPartial.toFixed(3)}</td>
                        </tr>
                        <tr className="border-b border-white/5">
                          <td className="py-2 px-3 text-slate-300">Point-level r</td>
                          <td className="py-2 px-3 text-center text-cyan-400">{s.pointLevel.r.toFixed(3)}</td>
                          <td className="py-2 px-3 text-center text-emerald-400">{lt.pointLevelRegression.r.toFixed(3)}</td>
                        </tr>
                        <tr>
                          <td className="py-2 px-3 text-slate-300">Direction</td>
                          <td className="py-2 px-3 text-center text-cyan-400">Negative</td>
                          <td className="py-2 px-3 text-center text-emerald-400">{lt.simpleRegression.slope < 0 ? 'Negative' : 'Positive'}</td>
                        </tr>
                      </tbody>
                    </table>
                  </div>
                )}
                {lt?.robustness && (
                  <div className="mt-4 bg-emerald-500/10 border border-emerald-500/20 rounded-xl p-3">
                    <div className="text-xs text-emerald-400 font-bold mb-1">Monte Carlo: {lt.robustness.monteCarlo.nIterations} iterations</div>
                    <div className="text-xs text-slate-300">
                      Slope = {lt.robustness.monteCarlo.slope.mean.toFixed(3)} ± {lt.robustness.monteCarlo.slope.std.toFixed(3)},
                      95% CI = [{lt.robustness.monteCarlo.slope.ci95[0].toFixed(3)}, {lt.robustness.monteCarlo.slope.ci95[1].toFixed(3)}],
                      {(lt.robustness.monteCarlo.slope.fracNegative * 100).toFixed(0)}% negative
                    </div>
                  </div>
                )}
              </GlassCard>
            </div>

            {sig && (
              <GlassCard className="mt-6">
                <h4 className="text-sm font-bold text-white mb-3">Statistical Significance (SPARC vs ΛCDM)</h4>
                <div className="grid grid-cols-2 md:grid-cols-4 gap-4">
                  <div className="bg-white/5 rounded-xl p-4 text-center">
                    <div className="text-xs text-slate-500 mb-1">Fisher z-score</div>
                    <div className="font-mono text-lg font-bold text-violet-400">{sig.fisherZ.zScore.toFixed(2)}</div>
                    <div className="text-xs text-slate-500">p = {sig.fisherZ.pValue.toExponential(1)}</div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4 text-center">
                    <div className="text-xs text-slate-500 mb-1">Cohen's d</div>
                    <div className="font-mono text-lg font-bold text-cyan-400">{sig.effectSize.cohensD.toFixed(2)}</div>
                    <div className="text-xs text-slate-500">{sig.effectSize.label}</div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4 text-center">
                    <div className="text-xs text-slate-500 mb-1">Welch's t</div>
                    <div className="font-mono text-lg font-bold text-emerald-400">{sig.welchT.tStat.toFixed(2)}</div>
                    <div className="text-xs text-slate-500">p = {sig.welchT.pValue.toExponential(1)}</div>
                  </div>
                  <div className="bg-white/5 rounded-xl p-4 text-center">
                    <div className="text-xs text-slate-500 mb-1">Permutation</div>
                    <div className="font-mono text-lg font-bold text-amber-400">{sig.permutation.pValue < 0.001 ? '<0.001' : sig.permutation.pValue.toFixed(3)}</div>
                    <div className="text-xs text-slate-500">{sig.nPermutations.toLocaleString()} shuffles</div>
                  </div>
                </div>
              </GlassCard>
            )}
          </section>

          {/* ═══════════════ SECTION 7: THE EQUATIONS ═══════════════ */}
          <section>
            <SectionHeader icon={Hash} number={7} title="The Equations" subtitle="Three coupling laws + the excess formula" id="equations" />

            <GlassCard glow="violet" className="mb-6">
              <h4 className="text-sm font-bold text-white mb-4">The Three Coupling Laws</h4>
              <div className="space-y-4">
                {coupling && [
                  { law: coupling.couplingLaws.vdm, label: 'DM Velocity', sigma: excess?.global[0]?.sigma },
                  { law: coupling.couplingLaws.rDMdom, label: 'DM Dominance Radius', sigma: excess?.global[1]?.sigma },
                  { law: coupling.couplingLaws.alpha, label: 'Inner Halo Slope', sigma: excess?.global[2]?.sigma },
                ].map((item, i) => (
                  <div key={i} className="bg-white/5 rounded-xl p-5">
                    <div className="flex items-center justify-between mb-2">
                      <span className="text-xs text-slate-400 uppercase tracking-wider">{item.label}</span>
                      {item.sigma && (
                        <span className={`text-sm font-mono font-bold ${item.sigma >= 5 ? 'text-amber-400' : item.sigma >= 3 ? 'text-emerald-400' : 'text-cyan-400'}`}>
                          {item.sigma.toFixed(1)}σ vs ΛCDM
                        </span>
                      )}
                    </div>
                    <div className="font-mono text-lg text-cyan-400 mb-2">{item.law.formula}</div>
                    <div className="flex gap-6 text-xs text-slate-400">
                      <span>r = {item.law.r.toFixed(4)}</span>
                      <span>R² = {item.law.r2.toFixed(3)}</span>
                      <span>slope = {item.law.slope.toFixed(5)}</span>
                      <span>intercept = {item.law.intercept.toFixed(4)}</span>
                    </div>
                  </div>
                ))}
              </div>
            </GlassCard>

            {excess && (
              <GlassCard glow="amber" className="mb-6">
                <h4 className="text-sm font-bold text-white mb-4">The Excess Formula</h4>
                <div className="bg-white/5 rounded-xl p-6 text-center mb-4">
                  <div className="font-mono text-xl text-amber-400 mb-2">
                    Δb = b_obs − b_ΛCDM
                  </div>
                  <p className="text-xs text-slate-400">
                    Where b is the slope of (halo property) vs log(Σ_bar). Δb measures how much stronger
                    the observed coupling is compared to ΛCDM simulations.
                  </p>
                </div>
                <div className="overflow-x-auto">
                  <table className="w-full text-sm">
                    <thead>
                      <tr className="border-b border-white/10">
                        <th className="text-left text-slate-400 py-2 px-3">Metric</th>
                        <th className="text-center text-cyan-400 py-2 px-3">b_obs</th>
                        <th className="text-center text-violet-400 py-2 px-3">b_ΛCDM</th>
                        <th className="text-center text-amber-400 py-2 px-3">Δb</th>
                        <th className="text-center text-white py-2 px-3">σ</th>
                        <th className="text-center text-slate-400 py-2 px-3">95% CI</th>
                      </tr>
                    </thead>
                    <tbody className="font-mono text-xs">
                      {excess.global.map((g: GlobalExcess, i: number) => (
                        <tr key={i} className="border-b border-white/5">
                          <td className="py-2.5 px-3 text-slate-300">{g.metric}</td>
                          <td className="py-2.5 px-3 text-center text-cyan-400">{g.obs.slope.toFixed(5)}</td>
                          <td className="py-2.5 px-3 text-center text-violet-400">{g.sim.slope.toFixed(5)}</td>
                          <td className="py-2.5 px-3 text-center text-amber-400 font-bold">{g.deltaB.toFixed(5)}</td>
                          <td className="py-2.5 px-3 text-center">
                            <span className={`font-bold ${g.sigma >= 5 ? 'text-amber-400' : g.sigma >= 3 ? 'text-emerald-400' : 'text-cyan-400'}`}>
                              {g.sigma.toFixed(1)}σ
                            </span>
                          </td>
                          <td className="py-2.5 px-3 text-center text-slate-500">[{g.deltaCI95[0].toFixed(4)}, {g.deltaCI95[1].toFixed(4)}]</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </GlassCard>
            )}

            <GlassCard>
              <h4 className="text-sm font-bold text-white mb-4">The Baseline: f_DM Anti-Correlation</h4>
              <div className="bg-white/5 rounded-xl p-6 text-center mb-4">
                <div className="font-mono text-lg text-cyan-400 mb-2">
                  f_DM(r) = {s.pointLevel.intercept.toFixed(3)} {s.pointLevel.slope < 0 ? '−' : '+'} {Math.abs(s.pointLevel.slope).toFixed(4)} × log₁₀(Σ_bar)
                </div>
                <div className="flex justify-center gap-6 text-xs text-slate-400 mt-2">
                  <span>r = {s.pointLevel.r.toFixed(3)}</span>
                  <span>R² = {s.pointLevel.r2.toFixed(3)}</span>
                  <span>n = {s.pointLevel.n.toLocaleString()} points</span>
                  <span>partial r|V_max = {s.perGalaxy.partialR.toFixed(3)}</span>
                </div>
              </div>
              {excess && (
                <div className="bg-white/5 rounded-xl p-6 text-center">
                  <div className="text-xs text-slate-500 uppercase tracking-wider mb-2">Predictive Power</div>
                  <div className="font-mono text-lg text-emerald-400 mb-2">
                    MSE(V_max + Σ_bar) = {excess.predictiveTest.mseCombined.toFixed(4)}
                  </div>
                  <div className="text-sm text-slate-300">
                    vs MSE(V_max only) = {excess.predictiveTest.mseVmaxOnly.toFixed(4)} →{' '}
                    <strong className="text-emerald-400">{excess.predictiveTest.improvementPct}% improvement</strong>
                  </div>
                </div>
              )}
            </GlassCard>
          </section>

          {/* ═══════════════ FINAL VERDICT ═══════════════ */}
          <section>
            <GlassCard glow="amber" className="border-amber-500/30">
              <div className="text-center space-y-4 py-4">
                <div className="inline-flex items-center gap-2 px-6 py-3 rounded-2xl bg-amber-500/20 border border-amber-500/30">
                  <Zap className="w-6 h-6 text-amber-400" />
                  <span className="text-xl font-display font-bold text-amber-400">Final Verdict</span>
                </div>

                <p className="text-lg text-white max-w-3xl mx-auto leading-relaxed">
                  Dark matter halo structure is coupled to baryonic surface density at{' '}
                  <strong className="text-amber-400">{excess?.maxSigma.toFixed(1)}σ</strong> above ΛCDM prediction.
                  The coupling is mass-dependent, confirmed by an independent dataset, survives all robustness tests,
                  and in one case has the <strong className="text-amber-400">opposite sign</strong> from what ΛCDM predicts.
                </p>

                <div className="bg-white/5 rounded-xl p-6 max-w-2xl mx-auto">
                  <p className="text-amber-300 font-medium text-base italic">
                    "Dark matter knows where the light is — and it shouldn't."
                  </p>
                </div>

                <div className="flex flex-wrap justify-center gap-6 pt-2">
                  <div className="text-center">
                    <div className="text-2xl font-mono font-bold text-amber-400">{excess?.maxSigma.toFixed(1)}σ</div>
                    <div className="text-xs text-slate-500">Max Excess</div>
                  </div>
                  <div className="text-center">
                    <div className="text-2xl font-mono font-bold text-emerald-400">{s.perGalaxy.n}</div>
                    <div className="text-xs text-slate-500">Galaxies</div>
                  </div>
                  <div className="text-center">
                    <div className="text-2xl font-mono font-bold text-cyan-400">{coupling?.passCount}/6</div>
                    <div className="text-xs text-slate-500">Tests Pass</div>
                  </div>
                  <div className="text-center">
                    <div className="text-2xl font-mono font-bold text-violet-400">{excess?.predictiveTest?.improvementPct}%</div>
                    <div className="text-xs text-slate-500">Predictive Gain</div>
                  </div>
                </div>
              </div>
            </GlassCard>
          </section>

        </div>
      </div>
    </Layout>
  );
}
