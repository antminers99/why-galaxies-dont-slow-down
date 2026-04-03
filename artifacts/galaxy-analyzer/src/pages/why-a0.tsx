import React, { useState, useEffect } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import {
  RadarChart, Radar, PolarGrid, PolarAngleAxis, PolarRadiusAxis,
  ResponsiveContainer, BarChart, Bar, XAxis, YAxis, CartesianGrid, Tooltip, Cell,
  ScatterChart, Scatter, ReferenceLine, Label
} from 'recharts';
import {
  BookOpenCheck, AlertTriangle, CheckCircle2, XCircle, HelpCircle,
  Atom, Orbit, Globe2, Sigma, FlaskConical, Lightbulb, Search
} from 'lucide-react';

interface TheoryData {
  facts: Array<{ id: string; statement: string; evidence: string; caveat?: string; implication?: string; robustness?: string }>;
  conditions: Array<{ id: string; condition: string; test: string }>;
  families: Record<string, { name: string; scores: Record<string, number>; total: number }>;
  derivations: {
    deSitter: { a_dS: number; ratio_to_a0: number; verdict: string };
    hubble: { cH0: number; cH0_2pi: number; ratio_to_a0: number; verdict: string };
    unruh: { T_dS: number; T_Unruh_a0: number; ratio: number; verdict: string };
    dimensional: Array<{ formula: string; value: number; ratio: number }>;
  };
  mainQuestion: string;
  mostPromisingPaths: string[];
}

const ScoreIcon = ({ score }: { score: number }) => {
  if (score >= 1) return <CheckCircle2 className="w-4 h-4 text-emerald-400" />;
  if (score >= 0.5) return <HelpCircle className="w-4 h-4 text-amber-400" />;
  return <XCircle className="w-4 h-4 text-red-400" />;
};

const conditionLabels: Record<string, string> = {
  C1: 'Derive a\u2080',
  C2: 'Universality',
  C3: 'cH\u2080 link',
  C4: 'Newton limit',
  C5: 'Clusters',
  C6: 'Mechanism',
};

export default function WhyA0Page() {
  const [data, setData] = useState<TheoryData | null>(null);
  const [activeNotebook, setActiveNotebook] = useState(0);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'theoretical-investigation.json')
      .then(r => r.json())
      .then(d => setData(d))
      .catch(() => {});
  }, []);

  if (!data) {
    return (
      <Layout>
        <div className="text-center text-slate-400 mt-20">Loading theoretical investigation...</div>
      </Layout>
    );
  }

  const radarData = Object.keys(conditionLabels).map(key => ({
    condition: conditionLabels[key],
    modGrav: (data.families.modifiedGravity?.scores[key] ?? 0) * 100,
    baryonDM: (data.families.baryonDM?.scores[key] ?? 0) * 100,
    cosmological: (data.families.cosmological?.scores[key] ?? 0) * 100,
  }));

  const barData = [
    { name: 'Modified\nGravity', score: data.families.modifiedGravity?.total ?? 0, fill: '#f97316' },
    { name: 'Baryon-DM\nCoupling', score: data.families.baryonDM?.total ?? 0, fill: '#a78bfa' },
    { name: 'Cosmological\nEmergent', score: data.families.cosmological?.total ?? 0, fill: '#06b6d4' },
  ];

  const dimData = data.derivations.dimensional.map(d => ({
    formula: d.formula,
    ratio: d.ratio,
    value: d.value,
    fill: Math.abs(d.ratio - 1) < 0.2 ? '#06b6d4' : Math.abs(d.ratio - 1) < 0.5 ? '#f59e0b' : '#ef4444',
  }));

  const notebooks = [
    { label: 'Facts', icon: BookOpenCheck },
    { label: 'Candidates', icon: FlaskConical },
    { label: 'Derivations', icon: Sigma },
  ];

  return (
    <Layout>
      <div className="space-y-8 max-w-7xl mx-auto">
        <div className="text-center space-y-3">
          <div className="flex items-center justify-center gap-3">
            <Search className="w-8 h-8 text-cyan-400" />
            <h1 className="text-3xl md:text-4xl font-display font-bold text-white">
              Why Does a{'\u2080'} Exist?
            </h1>
          </div>
          <p className="text-slate-400 max-w-2xl mx-auto">
            Phase D: From phenomenon to explanation. Three notebooks for finding the physical origin of the acceleration floor.
          </p>
          <div className="flex items-center justify-center gap-6 text-xs font-mono text-slate-500">
            <span>5 facts</span>
            <span className="text-slate-600">{'•'}</span>
            <span>6 conditions</span>
            <span className="text-slate-600">{'•'}</span>
            <span>3 families</span>
            <span className="text-slate-600">{'•'}</span>
            <span>6 derivation attempts</span>
          </div>
        </div>

        <GlassCard glow="cyan">
          <div className="flex items-start gap-3">
            <Lightbulb className="w-6 h-6 text-cyan-400 mt-1 flex-shrink-0" />
            <div>
              <h3 className="text-lg font-display font-bold text-white mb-2">The Central Question</h3>
              <p className="text-xl text-cyan-300 font-display italic">
                {'"'}Why does cH{'\u2080'} appear as a local acceleration in galaxies?{'"'}
              </p>
              <p className="text-sm text-slate-400 mt-2">
                Not {'"'}what is dark matter?{'"'} Not {'"'}is MOND right?{'"'} The question is more precise:
                an acceleration scale a{'\u2080'} {'\u2248'} cH{'\u2080'}/2{'\u03C0'} appears universally in galaxy dynamics.
                Any successful theory must explain why.
              </p>
            </div>
          </div>
        </GlassCard>

        <div className="flex gap-2 justify-center">
          {notebooks.map((nb, i) => (
            <button
              key={i}
              onClick={() => setActiveNotebook(i)}
              className={`flex items-center gap-2 px-5 py-2.5 rounded-xl transition-all font-medium text-sm ${
                activeNotebook === i
                  ? 'bg-cyan-500/20 text-cyan-300 border border-cyan-500/40'
                  : 'bg-white/5 text-slate-400 border border-white/10 hover:bg-white/10'
              }`}
            >
              <nb.icon className="w-4 h-4" />
              {nb.label}
            </button>
          ))}
        </div>

        {activeNotebook === 0 && (
          <div className="space-y-6">
            <h2 className="text-xl font-display font-bold text-white flex items-center gap-2">
              <BookOpenCheck className="w-5 h-5 text-cyan-400" />
              Notebook 1: What Must Be Explained
            </h2>
            <p className="text-sm text-slate-400">Five empirical facts. No philosophy, no interpretation. Just what the data says.</p>

            <div className="grid gap-4">
              {data.facts.map((f, i) => (
                <GlassCard key={f.id} glow="none">
                  <div className="flex gap-4">
                    <div className="w-10 h-10 rounded-lg bg-cyan-500/20 flex items-center justify-center flex-shrink-0">
                      <span className="text-cyan-300 font-mono font-bold text-sm">{f.id}</span>
                    </div>
                    <div className="space-y-1 flex-1">
                      <p className="text-white font-medium">{f.statement}</p>
                      <p className="text-sm text-slate-400"><span className="text-slate-500">Evidence:</span> {f.evidence}</p>
                      {f.caveat && (
                        <p className="text-sm text-amber-400/80 flex items-center gap-1">
                          <AlertTriangle className="w-3 h-3" /> Caveat: {f.caveat}
                        </p>
                      )}
                      {f.implication && (
                        <p className="text-sm text-cyan-400/80">{'\u2192'} {f.implication}</p>
                      )}
                      {f.robustness && (
                        <p className="text-sm text-emerald-400/80 flex items-center gap-1">
                          <CheckCircle2 className="w-3 h-3" /> {f.robustness}
                        </p>
                      )}
                    </div>
                  </div>
                </GlassCard>
              ))}
            </div>

            <GlassCard glow="purple">
              <h3 className="text-lg font-display font-bold text-white mb-4">6 Conditions for a Successful Explanation</h3>
              <div className="grid md:grid-cols-2 gap-3">
                {data.conditions.map(c => (
                  <div key={c.id} className="flex gap-3 items-start">
                    <span className="text-xs font-mono text-purple-400 bg-purple-500/10 px-2 py-1 rounded flex-shrink-0">{c.id}</span>
                    <div>
                      <p className="text-sm text-white font-medium">{c.condition}</p>
                      <p className="text-xs text-slate-500">{c.test}</p>
                    </div>
                  </div>
                ))}
              </div>
            </GlassCard>
          </div>
        )}

        {activeNotebook === 1 && (
          <div className="space-y-6">
            <h2 className="text-xl font-display font-bold text-white flex items-center gap-2">
              <FlaskConical className="w-5 h-5 text-purple-400" />
              Notebook 2: Candidate Explanations
            </h2>

            <div className="grid lg:grid-cols-3 gap-4">
              {[
                {
                  key: 'modifiedGravity',
                  icon: Atom,
                  color: 'orange',
                  borderColor: 'border-orange-500/30',
                  bgColor: 'bg-orange-500/10',
                  textColor: 'text-orange-400',
                  title: 'Modified Gravity',
                  subtitle: 'MOND, AQUAL, TeVeS',
                  idea: 'The law of gravity itself changes below a\u2080',
                  strengths: ['Universal by construction', 'Newton limit built-in'],
                  weaknesses: ['a\u2080 inserted, not derived', 'Fails in clusters', 'cH\u2080 unexplained'],
                },
                {
                  key: 'baryonDM',
                  icon: Orbit,
                  color: 'purple',
                  borderColor: 'border-purple-500/30',
                  bgColor: 'bg-purple-500/10',
                  textColor: 'text-purple-400',
                  title: 'Baryon-DM Coupling',
                  subtitle: 'Feedback, SIDM, Superfluid',
                  idea: 'DM exists but its distribution tracks baryons',
                  strengths: ['Explains clusters', 'Has feedback mechanism', 'Newton limit natural'],
                  weaknesses: ['RAR too tight for stochastic process', 'No cH\u2080 explanation', 'a\u2080 not derived'],
                },
                {
                  key: 'cosmological',
                  icon: Globe2,
                  color: 'cyan',
                  borderColor: 'border-cyan-500/30',
                  bgColor: 'bg-cyan-500/10',
                  textColor: 'text-cyan-400',
                  title: 'Cosmological Emergent',
                  subtitle: 'Cosmic Floor, Verlinde',
                  idea: 'Cosmic expansion creates a local acceleration floor',
                  strengths: ['Best at cH\u2080', 'Natural universality', 'Testable via z-evolution'],
                  weaknesses: ['NO MECHANISM (main gap)', 'Partial cluster help only', 'a\u2080 motivated, not derived'],
                },
              ].map(fam => {
                const scores = data.families[fam.key]?.scores ?? {};
                return (
                  <GlassCard key={fam.key} glow="none">
                    <div className={`border-l-2 ${fam.borderColor} pl-4 space-y-3`}>
                      <div className="flex items-center gap-2">
                        <div className={`w-8 h-8 rounded-lg ${fam.bgColor} flex items-center justify-center`}>
                          <fam.icon className={`w-4 h-4 ${fam.textColor}`} />
                        </div>
                        <div>
                          <h3 className="font-display font-bold text-white text-sm">{fam.title}</h3>
                          <p className="text-xs text-slate-500">{fam.subtitle}</p>
                        </div>
                      </div>

                      <p className="text-xs text-slate-400 italic">{fam.idea}</p>

                      <div className="space-y-1">
                        {Object.entries(conditionLabels).map(([key, label]) => (
                          <div key={key} className="flex items-center gap-2 text-xs">
                            <ScoreIcon score={scores[key] ?? 0} />
                            <span className="text-slate-400">{label}</span>
                          </div>
                        ))}
                      </div>

                      <div className="pt-2 border-t border-white/5">
                        <div className={`text-lg font-mono font-bold ${fam.textColor}`}>
                          {(data.families[fam.key]?.total ?? 0).toFixed(1)}/6
                        </div>
                      </div>

                      <div className="space-y-1">
                        <p className="text-xs text-emerald-400 font-medium">Strengths:</p>
                        {fam.strengths.map((s, i) => (
                          <p key={i} className="text-xs text-slate-400 pl-2">+ {s}</p>
                        ))}
                        <p className="text-xs text-red-400 font-medium mt-1">Weaknesses:</p>
                        {fam.weaknesses.map((w, i) => (
                          <p key={i} className="text-xs text-slate-400 pl-2">- {w}</p>
                        ))}
                      </div>
                    </div>
                  </GlassCard>
                );
              })}
            </div>

            <div className="grid md:grid-cols-2 gap-6">
              <GlassCard glow="none">
                <h3 className="text-base font-display font-bold text-white mb-4">Radar Comparison</h3>
                <ResponsiveContainer width="100%" height={300}>
                  <RadarChart data={radarData}>
                    <PolarGrid stroke="#334155" />
                    <PolarAngleAxis dataKey="condition" tick={{ fill: '#94a3b8', fontSize: 11 }} />
                    <PolarRadiusAxis angle={90} domain={[0, 100]} tick={false} axisLine={false} />
                    <Radar name="Modified Gravity" dataKey="modGrav" stroke="#f97316" fill="#f97316" fillOpacity={0.15} strokeWidth={2} />
                    <Radar name="Baryon-DM" dataKey="baryonDM" stroke="#a78bfa" fill="#a78bfa" fillOpacity={0.15} strokeWidth={2} />
                    <Radar name="Cosmological" dataKey="cosmological" stroke="#06b6d4" fill="#06b6d4" fillOpacity={0.15} strokeWidth={2} />
                    <Tooltip contentStyle={{ backgroundColor: '#0f172a', border: '1px solid #334155', borderRadius: '8px', fontSize: 12 }} />
                  </RadarChart>
                </ResponsiveContainer>
                <div className="flex justify-center gap-4 mt-2 text-xs">
                  <span className="flex items-center gap-1"><span className="w-3 h-0.5 bg-orange-500 inline-block" /> Mod Gravity</span>
                  <span className="flex items-center gap-1"><span className="w-3 h-0.5 bg-purple-400 inline-block" /> Baryon-DM</span>
                  <span className="flex items-center gap-1"><span className="w-3 h-0.5 bg-cyan-400 inline-block" /> Cosmological</span>
                </div>
              </GlassCard>

              <GlassCard glow="none">
                <h3 className="text-base font-display font-bold text-white mb-4">Overall Scores</h3>
                <ResponsiveContainer width="100%" height={300}>
                  <BarChart data={barData} layout="vertical">
                    <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
                    <XAxis type="number" domain={[0, 6]} tick={{ fill: '#94a3b8', fontSize: 11 }} />
                    <YAxis type="category" dataKey="name" tick={{ fill: '#94a3b8', fontSize: 11 }} width={100} />
                    <Tooltip contentStyle={{ backgroundColor: '#0f172a', border: '1px solid #334155', borderRadius: '8px', fontSize: 12 }} />
                    <Bar dataKey="score" radius={[0, 6, 6, 0]}>
                      {barData.map((d, i) => <Cell key={i} fill={d.fill} />)}
                    </Bar>
                  </BarChart>
                </ResponsiveContainer>
              </GlassCard>
            </div>
          </div>
        )}

        {activeNotebook === 2 && (
          <div className="space-y-6">
            <h2 className="text-xl font-display font-bold text-white flex items-center gap-2">
              <Sigma className="w-5 h-5 text-cyan-400" />
              Notebook 3: Derivation Attempts
            </h2>
            <p className="text-sm text-slate-400">Can we derive a{'\u2080'} from deeper principles? Six attempts, scored by how close they get.</p>

            <div className="grid md:grid-cols-2 gap-4">
              <GlassCard glow="none">
                <div className="space-y-3">
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-mono text-cyan-400 bg-cyan-500/10 px-2 py-1 rounded">#1</span>
                    <h3 className="font-display font-bold text-white text-sm">de Sitter Horizon</h3>
                  </div>
                  <p className="text-xs text-slate-400">
                    a_dS = c{'\u221A'}({'\u039B'}/3) = 7.07{'\u00D7'}10{'\u207B\u00B9\u00B9'} m/s{'\u00B2'}
                  </p>
                  <div className="flex items-center gap-2">
                    <div className="flex-1 bg-white/5 rounded-full h-3 overflow-hidden">
                      <div className="h-full bg-amber-500 rounded-full" style={{ width: `${(data.derivations.deSitter.ratio_to_a0) * 100}%` }} />
                    </div>
                    <span className="text-xs font-mono text-amber-400">{data.derivations.deSitter.ratio_to_a0.toFixed(2)}{'\u00D7'} a{'\u2080'}</span>
                  </div>
                  <p className="text-xs text-slate-500">{data.derivations.deSitter.verdict}</p>
                </div>
              </GlassCard>

              <GlassCard glow="cyan">
                <div className="space-y-3">
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-mono text-cyan-400 bg-cyan-500/10 px-2 py-1 rounded">#2</span>
                    <h3 className="font-display font-bold text-white text-sm">Hubble Scale (Best Match)</h3>
                    <span className="text-xs bg-cyan-500/20 text-cyan-300 px-2 py-0.5 rounded-full">CLOSEST</span>
                  </div>
                  <p className="text-xs text-slate-400">
                    cH{'\u2080'}/2{'\u03C0'} = 1.04{'\u00D7'}10{'\u207B\u00B9\u2070'} m/s{'\u00B2'}
                  </p>
                  <div className="flex items-center gap-2">
                    <div className="flex-1 bg-white/5 rounded-full h-3 overflow-hidden">
                      <div className="h-full bg-cyan-500 rounded-full" style={{ width: `${(data.derivations.hubble.ratio_to_a0) * 100}%` }} />
                    </div>
                    <span className="text-xs font-mono text-cyan-400">{data.derivations.hubble.ratio_to_a0.toFixed(2)}{'\u00D7'} a{'\u2080'}</span>
                  </div>
                  <p className="text-xs text-slate-500">{data.derivations.hubble.verdict}</p>
                  <div className="bg-cyan-500/5 border border-cyan-500/20 rounded-lg p-2">
                    <p className="text-xs text-cyan-300">The 2{'\u03C0'} could arise from: angular vs linear frequency, horizon circumference, circular orbit quantization, or vacuum mode Fourier structure.</p>
                  </div>
                </div>
              </GlassCard>

              <GlassCard glow="none">
                <div className="space-y-3">
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-mono text-cyan-400 bg-cyan-500/10 px-2 py-1 rounded">#3</span>
                    <h3 className="font-display font-bold text-white text-sm">Unruh Temperature</h3>
                  </div>
                  <p className="text-xs text-slate-400">
                    If universe has T_min = T_dS, then a_min = 2{'\u03C0'}ck_BT_dS/{'\u210F'}
                  </p>
                  <p className="text-xs text-slate-400">
                    T_dS / T_Unruh(a{'\u2080'}) = {data.derivations.unruh.ratio.toFixed(3)}
                  </p>
                  <p className="text-xs text-amber-400">Provides a REASON (minimum temperature) but gives same number as de Sitter</p>
                </div>
              </GlassCard>

              <GlassCard glow="none">
                <div className="space-y-3">
                  <div className="flex items-center gap-2">
                    <span className="text-xs font-mono text-cyan-400 bg-cyan-500/10 px-2 py-1 rounded">#4</span>
                    <h3 className="font-display font-bold text-white text-sm">H{'\u2080'} vs {'\u039B'} Discriminant</h3>
                  </div>
                  <p className="text-xs text-slate-400">
                    If a{'\u2080'} = cH/2{'\u03C0'}: changes with redshift (testable!)<br/>
                    If a{'\u2080'} = c{'\u221A'}({'\u039B'}/3): constant forever
                  </p>
                  <div className="bg-emerald-500/10 border border-emerald-500/20 rounded-lg p-2">
                    <p className="text-xs text-emerald-300 font-medium">These predictions DIVERGE at z {'>'} 0.5 — this is our measurable test</p>
                  </div>
                </div>
              </GlassCard>
            </div>

            <GlassCard glow="none">
              <h3 className="text-base font-display font-bold text-white mb-4">Dimensional Analysis: How Close Can We Get?</h3>
              <ResponsiveContainer width="100%" height={250}>
                <ScatterChart>
                  <CartesianGrid strokeDasharray="3 3" stroke="#1e293b" />
                  <XAxis
                    type="category"
                    dataKey="formula"
                    tick={{ fill: '#94a3b8', fontSize: 10 }}
                    allowDuplicatedCategory={false}
                  />
                  <YAxis
                    domain={[0, 6]}
                    tick={{ fill: '#94a3b8', fontSize: 11 }}
                    label={{ value: 'Ratio to a\u2080', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 11 }}
                  />
                  <ReferenceLine y={1} stroke="#06b6d4" strokeDasharray="5 5" strokeWidth={2}>
                    <Label value="a\u2080 = 1.0" fill="#06b6d4" fontSize={11} position="right" />
                  </ReferenceLine>
                  <Tooltip
                    contentStyle={{ backgroundColor: '#0f172a', border: '1px solid #334155', borderRadius: '8px', fontSize: 12 }}
                    formatter={(v: number) => [v.toFixed(3), 'Ratio']}
                  />
                  <Scatter data={dimData} dataKey="ratio">
                    {dimData.map((d, i) => <Cell key={i} fill={d.fill} r={8} />)}
                  </Scatter>
                </ScatterChart>
              </ResponsiveContainer>
              <p className="text-xs text-slate-500 text-center mt-2">
                <span className="text-cyan-400">{'\u25CF'}</span> Within 20% of a{'\u2080'}
                {' '}<span className="text-amber-400">{'\u25CF'}</span> Within 50%
                {' '}<span className="text-red-400">{'\u25CF'}</span> More than 50% off
              </p>
            </GlassCard>
          </div>
        )}

        <GlassCard glow="amber">
          <h3 className="text-lg font-display font-bold text-white mb-3 flex items-center gap-2">
            <AlertTriangle className="w-5 h-5 text-amber-400" />
            Honest Status
          </h3>
          <div className="grid md:grid-cols-2 gap-6">
            <div>
              <h4 className="text-sm font-medium text-emerald-400 mb-2">What We Have</h4>
              <ul className="space-y-1.5 text-sm text-slate-300">
                <li className="flex items-start gap-2"><CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" /> A well-defined question</li>
                <li className="flex items-start gap-2"><CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" /> 5 empirical facts surviving robustness tests</li>
                <li className="flex items-start gap-2"><CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" /> 6 clear conditions any explanation must meet</li>
                <li className="flex items-start gap-2"><CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" /> A testable prediction (a{'\u2080'}(z) divergence at z {'>'} 0.5)</li>
                <li className="flex items-start gap-2"><CheckCircle2 className="w-4 h-4 text-emerald-400 mt-0.5 flex-shrink-0" /> A measurement pipeline ready for JWST data</li>
              </ul>
            </div>
            <div>
              <h4 className="text-sm font-medium text-red-400 mb-2">What We Don{"'"}t Have</h4>
              <ul className="space-y-1.5 text-sm text-slate-300">
                <li className="flex items-start gap-2"><XCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" /> A mechanism (field equation, action, or Lagrangian)</li>
                <li className="flex items-start gap-2"><XCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" /> Physical origin of the 2{'\u03C0'} factor</li>
                <li className="flex items-start gap-2"><XCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" /> Resolution of the cluster problem (~52.5% residual)</li>
                <li className="flex items-start gap-2"><XCircle className="w-4 h-4 text-red-400 mt-0.5 flex-shrink-0" /> Proof that cH{'\u2080'} link is not coincidence</li>
              </ul>
            </div>
          </div>
          <div className="mt-4 pt-4 border-t border-white/10">
            <p className="text-sm text-amber-300 italic text-center">
              {'"'}The gap between a motivated formula and a derived theory is where real physics lives.{'"'}
            </p>
          </div>
        </GlassCard>

        <GlassCard glow="none">
          <h3 className="text-lg font-display font-bold text-white mb-3">Most Promising Paths Forward</h3>
          <div className="grid md:grid-cols-3 gap-4">
            {data.mostPromisingPaths.map((p, i) => (
              <div key={i} className="flex items-start gap-3">
                <div className="w-7 h-7 rounded-full bg-cyan-500/20 flex items-center justify-center flex-shrink-0 mt-0.5">
                  <span className="text-cyan-300 font-mono text-xs font-bold">{i + 1}</span>
                </div>
                <p className="text-sm text-slate-300">{p}</p>
              </div>
            ))}
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
