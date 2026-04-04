import React, { useState, useEffect } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip,
  ResponsiveContainer, Cell, ReferenceLine, Label,
  ScatterChart, Scatter
} from 'recharts';
import {
  ScrollText, CheckCircle2, XCircle, AlertTriangle,
  Target, Layers, FlaskConical, Telescope, ArrowRight,
  ChevronDown, ChevronRight, BookOpen, Sigma, Atom,
  ShieldCheck, BarChart3, GitCompare, Search
} from 'lucide-react';

type Phase10Data = {
  nGalaxies: number;
  baselineTau: number;
  candidates: Array<{
    name: string;
    nValid: number;
    r: number;
    r2: number;
    pPerm: number;
    tauNew: number;
    deltaTauPct: number;
    deltaAIC: number;
    significant: boolean;
    tStat: number;
    beta: number;
    seBeta: number;
    pointRMS_M0: number;
    pointRMS_M1: number;
    pointImprPct: number;
  }>;
  verdict: string;
};

type FigureData = {
  nGalaxies: number;
  nPoints: number;
  rarPoints: { x: number; y: number }[];
  rarCurve: { x: number; y: number }[];
  axisRange: { xMin: number; xMax: number; yMin: number; yMax: number };
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

function Section({ id, title, icon: Icon, children, defaultOpen = false }: {
  id: string; title: string; icon: any; children: React.ReactNode; defaultOpen?: boolean;
}) {
  const [open, setOpen] = useState(defaultOpen);
  return (
    <div id={id} className="mb-4">
      <button
        onClick={() => setOpen(!open)}
        className="w-full flex items-center gap-3 p-4 bg-slate-800/50 hover:bg-slate-800/80 border border-white/5 rounded-xl transition-all text-left"
      >
        <Icon className="w-5 h-5 text-cyan-400 flex-shrink-0" />
        <span className="text-sm font-semibold text-white flex-1">{title}</span>
        {open ? <ChevronDown className="w-4 h-4 text-slate-400" /> : <ChevronRight className="w-4 h-4 text-slate-400" />}
      </button>
      {open && (
        <div className="mt-2 p-4 bg-slate-900/40 border border-white/5 rounded-xl text-sm text-slate-300 leading-relaxed space-y-3">
          {children}
        </div>
      )}
    </div>
  );
}

function Stat({ label, value, sub, color = 'cyan' }: { label: string; value: string; sub?: string; color?: string }) {
  const colors: Record<string, string> = {
    cyan: 'text-cyan-400', purple: 'text-purple-400', amber: 'text-amber-400',
    red: 'text-red-400', emerald: 'text-emerald-400',
  };
  return (
    <div className="p-3 bg-slate-800/60 border border-white/5 rounded-xl text-center">
      <div className="text-[10px] text-slate-500 uppercase tracking-wider mb-1">{label}</div>
      <div className={'text-lg font-bold font-mono ' + (colors[color] || 'text-white')}>{value}</div>
      {sub && <div className="text-[10px] text-slate-500 mt-1">{sub}</div>}
    </div>
  );
}

function Tag({ children, color = 'cyan' }: { children: React.ReactNode; color?: string }) {
  const styles: Record<string, string> = {
    cyan: 'bg-cyan-950/50 text-cyan-300 border-cyan-500/20',
    emerald: 'bg-emerald-950/50 text-emerald-300 border-emerald-500/20',
    amber: 'bg-amber-950/50 text-amber-300 border-amber-500/20',
    red: 'bg-red-950/50 text-red-300 border-red-500/20',
    purple: 'bg-purple-950/50 text-purple-300 border-purple-500/20',
  };
  return (
    <span className={'inline-block px-2 py-0.5 text-[10px] font-mono border rounded ' + (styles[color] || styles.cyan)}>
      {children}
    </span>
  );
}

export default function CanonicalSummaryPage() {
  const [p10, setP10] = useState<Phase10Data | null>(null);
  const [fig, setFig] = useState<FigureData | null>(null);

  useEffect(() => {
    const base = import.meta.env.BASE_URL;
    Promise.all([
      fetch(base + 'phase10-second-param-results.json').then(r => r.ok ? r.json() : null).catch(() => null),
      fetch(base + 'figure1-data.json').then(r => r.ok ? r.json() : null).catch(() => null),
    ]).then(([a, b]) => { setP10(a); setFig(b); });
  }, []);

  const candsSorted = p10 ? [...p10.candidates].sort((a, b) => a.deltaTauPct - b.deltaTauPct) : [];

  return (
    <Layout>
      <div className="flex items-center gap-3 mb-2">
        <ScrollText className="w-7 h-7 text-cyan-400" />
        <div>
          <h1 className="text-2xl font-display font-bold text-white">Canonical Summary</h1>
          <p className="text-sm text-slate-400">
            Definitive findings, methodology, and open questions
          </p>
        </div>
      </div>

      <GlassCard glow="cyan">
        <h2 className="text-base font-bold text-white mb-3">The Central Question</h2>
        <p className="text-slate-300 text-sm mb-4">
          Does a real acceleration transition scale exist in galaxy rotation curves?
          If so, can it be represented by a single characteristic value a&#8320; after sample
          cleaning and systematic treatment?
        </p>

        <div className="p-4 bg-gradient-to-r from-cyan-950/40 to-purple-950/40 border border-cyan-500/20 rounded-xl">
          <h3 className="text-sm font-bold text-cyan-400 mb-2">The Short Answer</h3>
          <p className="text-slate-200 text-sm">
            <strong>Yes</strong>, a real acceleration transition scale exists. Our best current estimate:
          </p>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3 mt-3">
            <Stat label="a&#8320;" value="3633" sub="(km/s)&#178;/kpc" />
            <Stat label="a&#8320;" value="1.18" sub="&#215;10&#8315;&#185;&#8304; m/s&#178;" color="purple" />
            <Stat label="tau" value="0.291" sub="dex" color="amber" />
            <Stat label="I&#178;" value="92.4%" sub="heterogeneity" color="red" />
          </div>
          <p className="text-[11px] text-slate-400 mt-3">
            This is <strong>not</strong> proof that a&#8320; is an exact universal constant.
            The literature value ~1.20&#215;10&#8315;&#185;&#8304; m/s&#178; falls within our range.
          </p>
        </div>
      </GlassCard>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-4 mt-4">

        <GlassCard glow="none">
          <h3 className="text-sm font-bold text-emerald-400 mb-3 flex items-center gap-2">
            <CheckCircle2 className="w-4 h-4" /> What Is Established
          </h3>
          <ul className="space-y-2 text-sm text-slate-300">
            <li className="flex items-start gap-2">
              <CheckCircle2 className="w-3.5 h-3.5 text-emerald-500 mt-0.5 flex-shrink-0" />
              <span>A real, tight RAR exists in the data</span>
            </li>
            <li className="flex items-start gap-2">
              <CheckCircle2 className="w-3.5 h-3.5 text-emerald-500 mt-0.5 flex-shrink-0" />
              <span>A transition scale of order ~10&#8315;&#185;&#8304; m/s&#178; exists</span>
            </li>
            <li className="flex items-start gap-2">
              <CheckCircle2 className="w-3.5 h-3.5 text-emerald-500 mt-0.5 flex-shrink-0" />
              <span>Result survives sample cleaning, Y&#8902; marginalization, distance marginalization, hierarchical modeling</span>
            </li>
            <li className="flex items-start gap-2">
              <CheckCircle2 className="w-3.5 h-3.5 text-emerald-500 mt-0.5 flex-shrink-0" />
              <span>Result survives kinematic contamination tests</span>
            </li>
            <li className="flex items-start gap-2">
              <CheckCircle2 className="w-3.5 h-3.5 text-emerald-500 mt-0.5 flex-shrink-0" />
              <span>Within-galaxy scatter (0.122 dex) matches literature (0.13 dex)</span>
            </li>
          </ul>
        </GlassCard>

        <GlassCard glow="none">
          <h3 className="text-sm font-bold text-red-400 mb-3 flex items-center gap-2">
            <XCircle className="w-4 h-4" /> What Is NOT Established
          </h3>
          <ul className="space-y-2 text-sm text-slate-300">
            <li className="flex items-start gap-2">
              <XCircle className="w-3.5 h-3.5 text-red-500 mt-0.5 flex-shrink-0" />
              <span>That a&#8320; is an exact universal constant across all galaxies</span>
            </li>
            <li className="flex items-start gap-2">
              <XCircle className="w-3.5 h-3.5 text-red-500 mt-0.5 flex-shrink-0" />
              <span>That MOND is proven</span>
            </li>
            <li className="flex items-start gap-2">
              <XCircle className="w-3.5 h-3.5 text-red-500 mt-0.5 flex-shrink-0" />
              <span>That dark matter is ruled out</span>
            </li>
            <li className="flex items-start gap-2">
              <XCircle className="w-3.5 h-3.5 text-red-500 mt-0.5 flex-shrink-0" />
              <span>That the cosmological relation a&#8320; = cH&#8320;/2&#960; is physically confirmed</span>
            </li>
            <li className="flex items-start gap-2">
              <XCircle className="w-3.5 h-3.5 text-red-500 mt-0.5 flex-shrink-0" />
              <span>That galaxy-level heterogeneity is fully explained</span>
            </li>
          </ul>
        </GlassCard>
      </div>

      <GlassCard glow="purple" className="mt-4">
        <h2 className="text-base font-bold text-white mb-3 flex items-center gap-2">
          <Layers className="w-5 h-5 text-purple-400" /> The Critical Distinction
        </h2>
        <p className="text-sm text-slate-300 mb-4">
          Two types of scatter must never be conflated:
        </p>
        <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
          <div className="p-4 bg-emerald-950/30 border border-emerald-500/20 rounded-xl">
            <h4 className="text-xs font-bold text-emerald-400 uppercase tracking-wider mb-2">
              Point-Level Scatter
            </h4>
            <p className="text-sm text-slate-300 mb-2">
              How tightly do individual data points within each galaxy follow
              the global RAR curve?
            </p>
            <div className="flex items-center gap-2">
              <span className="text-2xl font-bold font-mono text-emerald-400">0.122</span>
              <span className="text-xs text-slate-400">dex (within-galaxy)</span>
            </div>
            <p className="text-[10px] text-emerald-400/70 mt-1">
              Matches McGaugh+2016 (0.13 dex)
            </p>
          </div>
          <div className="p-4 bg-amber-950/30 border border-amber-500/20 rounded-xl">
            <h4 className="text-xs font-bold text-amber-400 uppercase tracking-wider mb-2">
              Galaxy-Level Heterogeneity
            </h4>
            <p className="text-sm text-slate-300 mb-2">
              How much do galaxies differ from each other in their best-fit a&#8320;?
            </p>
            <div className="flex items-center gap-2">
              <span className="text-2xl font-bold font-mono text-amber-400">0.291</span>
              <span className="text-xs text-slate-400">dex (between-galaxy tau)</span>
            </div>
            <p className="text-[10px] text-amber-400/70 mt-1">
              Different quantity from literature scatter
            </p>
          </div>
        </div>
        <p className="text-xs text-slate-400 mt-3">
          These are <strong>different quantities</strong>. The RAR is tight at the point level
          but heterogeneous at the galaxy level. Both findings are correct and compatible.
        </p>
      </GlassCard>

      <div className="mt-4 space-y-0">

        <Section id="data" title="Data & Sample" icon={BarChart3} defaultOpen={false}>
          <div className="grid grid-cols-2 md:grid-cols-4 gap-3 mb-3">
            <Stat label="Original" value="175+22" sub="SPARC + LT galaxies" />
            <Stat label="Total points" value="3,755" sub="v4.0 radial" color="purple" />
            <Stat label="GOLD+i45" value="59" sub="clean sample" color="emerald" />
            <Stat label="Clean points" value="1,789" sub="GOLD+i45" color="amber" />
          </div>
          <p>
            <strong>GOLD+i45 criteria:</strong> V&#8344;&#8336;&#8339; &#8805; 50 km/s,
            inclination &#8805; 45&#176;, n &#8805; 5 points, g&#8346;&#8336;&#8339; dynamic range &#8805; 1.0 dex.
          </p>
        </Section>

        <Section id="model" title="Model: McGaugh RAR Only" icon={Sigma} defaultOpen={false}>
          <div className="p-3 bg-slate-800/60 border border-white/5 rounded-xl font-mono text-sm text-center mb-3">
            g&#8348;&#8346;&#8347; = g&#8346;&#8336;&#8339; / (1 - exp(-&#8730;(g&#8346;&#8336;&#8339; / a&#8320;)))
          </div>
          <p>
            A single estimator was used throughout to prevent mixing different functional forms
            within the same inference. a&#8320; is a phenomenological parameter extracted from data,
            not derived from first principles.
          </p>
        </Section>

        <Section id="pipeline" title="v4.0 Pipeline Upgrades" icon={GitCompare} defaultOpen={false}>
          <div className="space-y-2">
            <div className="flex items-start gap-2">
              <ArrowRight className="w-3.5 h-3.5 text-cyan-400 mt-1 flex-shrink-0" />
              <span>Full sample instead of subsampling</span>
            </div>
            <div className="flex items-start gap-2">
              <ArrowRight className="w-3.5 h-3.5 text-cyan-400 mt-1 flex-shrink-0" />
              <span>Per-galaxy Y&#8902; marginalization over [0.3, 0.8]</span>
            </div>
            <div className="flex items-start gap-2">
              <ArrowRight className="w-3.5 h-3.5 text-cyan-400 mt-1 flex-shrink-0" />
              <span>Distance marginalization within reported uncertainties</span>
            </div>
            <div className="flex items-start gap-2">
              <ArrowRight className="w-3.5 h-3.5 text-cyan-400 mt-1 flex-shrink-0" />
              <span>DerSimonian-Laird hierarchical random-effects model</span>
            </div>
          </div>
          <p className="mt-3 text-slate-400">
            This moved a&#8320; from 3374 (v3.0) to 3633 (v4.0) &mdash; toward the literature value.
            However, tau did not decrease substantially: the galaxy-level heterogeneity is
            not merely a Y&#8902; artifact.
          </p>
        </Section>

        <Section id="systematics" title="Systematic Tests Passed" icon={ShieldCheck} defaultOpen={false}>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
            <div className="p-3 bg-emerald-950/20 border border-emerald-500/10 rounded-lg">
              <h5 className="text-xs font-bold text-emerald-400 mb-1">Non-circular motions</h5>
              <p className="text-xs text-slate-400">Inner vs outer split was small. Not the source of a&#8320;.</p>
            </div>
            <div className="p-3 bg-emerald-950/20 border border-emerald-500/10 rounded-lg">
              <h5 className="text-xs font-bold text-emerald-400 mb-1">Pressure support</h5>
              <p className="text-xs text-slate-400">Small correction. Did not destroy the signal or the scale.</p>
            </div>
            <div className="p-3 bg-emerald-950/20 border border-emerald-500/10 rounded-lg">
              <h5 className="text-xs font-bold text-emerald-400 mb-1">Quality / Morphology</h5>
              <p className="text-xs text-slate-400">Effects small to moderate. Not destructive to the result.</p>
            </div>
            <div className="p-3 bg-emerald-950/20 border border-emerald-500/10 rounded-lg">
              <h5 className="text-xs font-bold text-emerald-400 mb-1">Distance method</h5>
              <p className="text-xs text-slate-400">b=0.061, t=0.56, p=0.55. NOT significant after covariate control.</p>
            </div>
          </div>
          <div className="mt-3 p-3 bg-amber-950/20 border border-amber-500/10 rounded-lg">
            <h5 className="text-xs font-bold text-amber-400 mb-1">Remaining concern</h5>
            <p className="text-xs text-slate-400">
              Distance and inclination splits show large effect sizes but
              fail to reach statistical significance (jackknife t &lt; 2).
              Status: <Tag color="amber">unresolved, not decisive</Tag>
            </p>
          </div>
        </Section>

        <Section id="tau" title="What tau = 0.291 Means" icon={AlertTriangle} defaultOpen={true}>
          <p>
            This is arguably the most important finding in the entire project.
          </p>
          <p className="mt-2">
            tau means there is <strong>real scatter between galaxies</strong> in their
            best-fit a&#8320; values that cannot be explained by within-galaxy fitting noise alone.
          </p>
          <p className="mt-2">
            However, this does <strong>not</strong> automatically mean a&#8320; itself varies.
            It could also mean:
          </p>
          <ul className="mt-2 space-y-1 ml-4">
            <li className="flex items-start gap-2">
              <span className="text-cyan-400 font-bold">1.</span>
              <span>a&#8320; is universal, but there is an additional galaxy-dependent correction term</span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-cyan-400 font-bold">2.</span>
              <span>A hidden systematic variable biases the inference per galaxy</span>
            </li>
            <li className="flex items-start gap-2">
              <span className="text-cyan-400 font-bold">3.</span>
              <span>A coupling with some galaxy property not explicitly in the model</span>
            </li>
          </ul>
        </Section>

        <Section id="phase10" title="Phase 10: Second-Parameter Search" icon={Search} defaultOpen={true}>
          <p className="mb-3">
            Tested 10 candidate second parameters to see if any explains the galaxy-level heterogeneity.
            The model: log(a&#8320;)&#8337;&#8342;&#8342; = &#956; + &#946; &#215; X&#8345;&#8338;&#8339;&#8344;.
          </p>
          {p10 && (
            <>
              <div className="overflow-x-auto">
                <table className="w-full text-xs font-mono">
                  <thead>
                    <tr className="text-slate-500 border-b border-white/5">
                      <th className="text-left py-2 pr-2">#</th>
                      <th className="text-left py-2 pr-4">Candidate</th>
                      <th className="text-right py-2 pr-4">r</th>
                      <th className="text-right py-2 pr-4">p(perm)</th>
                      <th className="text-right py-2 pr-4">dAIC</th>
                      <th className="text-right py-2 pr-4">tau drop</th>
                      <th className="text-center py-2">Sig?</th>
                    </tr>
                  </thead>
                  <tbody>
                    {candsSorted.map((c, i) => (
                      <tr key={c.name} className={c.significant ? 'text-emerald-300' : 'text-slate-400'}>
                        <td className="py-1.5 pr-2">{i + 1}</td>
                        <td className="py-1.5 pr-4">{c.name}</td>
                        <td className="py-1.5 pr-4 text-right">{c.r > 0 ? '+' : ''}{c.r.toFixed(3)}</td>
                        <td className="py-1.5 pr-4 text-right">{c.pPerm.toFixed(3)}</td>
                        <td className="py-1.5 pr-4 text-right">{c.deltaAIC > 0 ? '+' : ''}{c.deltaAIC.toFixed(1)}</td>
                        <td className="py-1.5 pr-4 text-right">{c.deltaTauPct.toFixed(1)}%</td>
                        <td className="py-1.5 text-center">{c.significant ? '***' : ''}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>

              <div className="mt-4 grid grid-cols-1 md:grid-cols-3 gap-3">
                {candsSorted.filter(c => c.significant).map(c => (
                  <div key={c.name} className="p-3 bg-emerald-950/20 border border-emerald-500/15 rounded-xl">
                    <div className="text-xs font-bold text-emerald-400 mb-1">{c.name}</div>
                    <div className="text-xs text-slate-400 space-y-0.5">
                      <div>r = {c.r > 0 ? '+' : ''}{c.r.toFixed(3)}, t = {c.tStat.toFixed(2)}</div>
                      <div>p(perm) = {c.pPerm.toFixed(4)}</div>
                      <div>tau: {p10.baselineTau.toFixed(3)} &#8594; {c.tauNew.toFixed(3)} ({c.deltaTauPct.toFixed(1)}%)</div>
                      <div>Point RMS: {c.pointImprPct.toFixed(1)}% improvement</div>
                    </div>
                  </div>
                ))}
              </div>

              <div className="mt-4 p-3 bg-amber-950/20 border border-amber-500/15 rounded-xl">
                <p className="text-xs text-amber-300">
                  <strong>Key result:</strong> Even the best candidate (eta_rot) reduces tau by only 5.7%.
                  Over 94% of galaxy-level heterogeneity remains unexplained by any measured property.
                </p>
              </div>
            </>
          )}
        </Section>

        <Section id="hypothesis" title="The Alternative Hypothesis" icon={Atom} defaultOpen={true}>
          <div className="p-4 bg-gradient-to-r from-purple-950/30 to-cyan-950/30 border border-purple-500/15 rounded-xl">
            <p className="text-sm text-slate-200 italic">
              &ldquo;Perhaps a&#8320; is approximately universal, but there exists a
              galaxy-dependent addition/shift on top of the simple RAR model.&rdquo;
            </p>
          </div>
          <p className="mt-3">
            This hypothesis remains <strong>open and unresolved</strong>.
            Phase 10 found suggestive candidates (rotation curve shape, Hubble type)
            but none explains the bulk of heterogeneity.
          </p>
          <p className="mt-2">
            The stronger framing is not &ldquo;a&#8320; varies from galaxy to galaxy&rdquo; but rather:
          </p>
          <div className="mt-2 p-3 bg-cyan-950/20 border border-cyan-500/15 rounded-xl text-sm text-cyan-200">
            a&#8320; may be truly universal, but the effective RAR may depend on a second
            variable not represented in the single-parameter model.
          </div>
        </Section>

        <Section id="literature" title="Position in the Literature" icon={BookOpen} defaultOpen={false}>
          <div className="space-y-3">
            <div className="p-3 bg-slate-800/40 rounded-lg">
              <h5 className="text-xs font-bold text-emerald-400 mb-1">Compatible with:</h5>
              <ul className="text-xs text-slate-400 space-y-1 ml-2">
                <li>Tight RAR at point level</li>
                <li>a&#8320; of order 1.2&#215;10&#8315;&#185;&#8304; m/s&#178;</li>
                <li>Priors and nuisance treatment affect the inferred value</li>
                <li>Universality is not a settled question</li>
              </ul>
            </div>
            <div className="p-3 bg-slate-800/40 rounded-lg">
              <h5 className="text-xs font-bold text-cyan-400 mb-1">Our addition:</h5>
              <ul className="text-xs text-slate-400 space-y-1 ml-2">
                <li>Point-level tightness &#8800; galaxy-level homogeneity &mdash; must be separated</li>
                <li>Clean pipeline that explicitly decomposes the two</li>
                <li>Systematic search over 10 second parameters &mdash; none decisive</li>
              </ul>
            </div>
          </div>
        </Section>

        <Section id="dontsay" title="What We Must NOT Claim" icon={XCircle} defaultOpen={false}>
          <div className="grid grid-cols-1 md:grid-cols-2 gap-2">
            {[
              'We discovered that a\u2080 is definitively non-universal',
              'We solved the dark matter puzzle',
              'We proved MOND',
              'We confirmed the cosmological origin of a\u2080',
            ].map((s, i) => (
              <div key={i} className="flex items-center gap-2 p-2 bg-red-950/20 border border-red-500/10 rounded-lg">
                <XCircle className="w-3.5 h-3.5 text-red-500 flex-shrink-0" />
                <span className="text-xs text-red-300">{s}</span>
              </div>
            ))}
          </div>
          <p className="mt-3 text-xs text-slate-400">
            These claims exceed the current evidence.
          </p>
        </Section>

      </div>

      <GlassCard glow="cyan" className="mt-4">
        <h2 className="text-base font-bold text-white mb-3 flex items-center gap-2">
          <Target className="w-5 h-5 text-cyan-400" /> The Correct Claim
        </h2>
        <div className="p-4 bg-gradient-to-r from-cyan-950/40 to-emerald-950/40 border border-cyan-500/20 rounded-xl">
          <p className="text-sm text-slate-200 leading-relaxed">
            A real acceleration transition scale exists in galaxy rotation curves.
            After sample cleaning and treatment of Y&#8902; and distance via hierarchical
            modeling, we obtain:
          </p>
          <p className="text-lg font-bold font-mono text-cyan-400 my-3 text-center">
            a&#8320; = 3633 (km/s)&#178;/kpc = 1.18&#215;10&#8315;&#185;&#8304; m/s&#178;
          </p>
          <p className="text-sm text-slate-200 leading-relaxed">
            This is consistent with the literature, but non-negligible heterogeneity
            persists between galaxies (tau = 0.291 dex, I&#178; = 92.4%).
            Therefore, it is <strong>not yet possible to claim that a&#8320; is an exact
            universal constant</strong>.
          </p>
        </div>
      </GlassCard>

      <GlassCard glow="amber" className="mt-4">
        <h2 className="text-base font-bold text-white mb-3 flex items-center gap-2">
          <Telescope className="w-5 h-5 text-amber-400" /> The Open Question
        </h2>
        <div className="p-4 bg-amber-950/30 border border-amber-500/15 rounded-xl">
          <p className="text-sm text-slate-200 leading-relaxed mb-3">
            The question is no longer <em>&ldquo;Does a scale exist?&rdquo;</em>
          </p>
          <p className="text-sm text-amber-200 font-medium leading-relaxed">
            It is now: <em>&ldquo;Does the galaxy-level variation mean a&#8320; itself differs,
            or is there an additional term on top of the simple RAR?&rdquo;</em>
          </p>
        </div>
        <div className="mt-3 grid grid-cols-1 md:grid-cols-2 gap-3">
          <div className="p-3 bg-slate-800/40 border border-white/5 rounded-lg">
            <h5 className="text-xs font-bold text-purple-400 mb-2">Model A</h5>
            <p className="text-xs text-slate-400">
              a&#8320; universal + scatter only
            </p>
          </div>
          <div className="p-3 bg-slate-800/40 border border-white/5 rounded-lg">
            <h5 className="text-xs font-bold text-cyan-400 mb-2">Model B</h5>
            <p className="text-xs text-slate-400">
              a&#8320; universal + galaxy-dependent correction term (EFE, morphology,
              kinematic bias, gas fraction, environment)
            </p>
          </div>
        </div>
        <p className="text-xs text-slate-400 mt-3">
          Phase 10 tested Model B with available data. eta_rot is suggestive but explains only ~6%.
          The bulk of heterogeneity requires unmeasured variables or is genuinely intrinsic.
        </p>
      </GlassCard>

      <div className="mt-4 p-4 bg-slate-800/30 border border-white/5 rounded-xl">
        <h3 className="text-xs font-bold text-slate-400 uppercase tracking-wider mb-2">
          Quick Reference for AI Systems
        </h3>
        <div className="font-mono text-[11px] text-slate-500 space-y-1">
          <div>Problem: existence and universality of galaxy acceleration scale</div>
          <div>Dataset: 197 galaxies, 3755 radial points</div>
          <div>Clean sample: GOLD+i45, 59 galaxies, 1789 points</div>
          <div>Model: McGaugh RAR only</div>
          <div>Estimator: hierarchical DL with per-galaxy Y* and distance marginalization</div>
          <div>Best value: a0 = 3633 (km/s)^2/kpc = 1.18e-10 m/s^2</div>
          <div>Statistical: log a0 = 3.560 +/- 0.041 dex</div>
          <div>Heterogeneity: tau = 0.291 dex, I^2 = 92.4%</div>
          <div>Within-galaxy scatter: 0.122 dex (matches McGaugh 0.13)</div>
          <div>Literature: consistent | Universality: not established</div>
          <div>Best 2nd param: eta_rot (r=0.334, p=0.010, -5.7% tau)</div>
          <div>Cosmological: suggestive only | a0/(cH0/2pi) = 1.130</div>
        </div>
      </div>

    </Layout>
  );
}
