import React, { useEffect, useState } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { BarChart, Bar, XAxis, YAxis, Tooltip, ResponsiveContainer, ReferenceLine, Cell } from 'recharts';
import { Atom, Lightbulb, FlaskConical, AlertTriangle, CheckCircle2, XCircle, ArrowRight, Zap, BookOpen } from 'lucide-react';

interface SimA0Data {
  observational: {
    combined: { globalA0: number; globalRMS: number; madA0: number; nGalaxies: number; nPoints: number; bins: Array<{ logGbar: number; medRatio: number; n: number }> };
  };
  simulations: Array<{ name: string; ref: string; a0_kpc: number | null; a0_reported: string; scatter: string; universal: boolean; feedbackDependent: string; notes: string }>;
  keyFindings: { A_exists: string; B_universal: string; C_matchesObs: string; D_simsAgree: string; tension: string; unexplained: string };
  verdict: string;
}

interface TransitionData {
  a0_corrected: number;
  nPoints: number;
  nGalaxies: number;
  collapse: { rmsWithCorrectA0: number };
  cosmology: { a0: number; cH0: string; ratio: number };
  ratioBins: Array<{ logGbar: number; medianRatio: number; n: number }>;
}

interface BivariateData {
  honestControl?: {
    marginalSigmaBar: { crossValidatedImprovement: number; fStatistic: number };
    verdict: string;
  };
}

const A0_MS2 = 1.2e-10;
const C_KMS = 299792.458;
const H0_KMS_MPC = 70;
const H0_PER_S = H0_KMS_MPC / 3.0857e19;

export default function ModelProposalPage() {
  const [simA0, setSimA0] = useState<SimA0Data | null>(null);
  const [transition, setTransition] = useState<TransitionData | null>(null);
  const [bivariate, setBivariate] = useState<BivariateData | null>(null);

  useEffect(() => {
    Promise.all([
      fetch(import.meta.env.BASE_URL + 'sim-a0-comparison.json').then(r => r.ok ? r.json() : null).catch(() => null),
      fetch(import.meta.env.BASE_URL + 'transition-scale.json').then(r => r.ok ? r.json() : null).catch(() => null),
      fetch(import.meta.env.BASE_URL + 'bivariate-collapse.json').then(r => r.ok ? r.json() : null).catch(() => null),
    ]).then(([s, t, b]) => { setSimA0(s); setTransition(t); setBivariate(b); });
  }, []);

  const epsilon = A0_MS2 / (C_KMS * 1000 * H0_PER_S);
  const cH0_ms2 = C_KMS * 1000 * H0_PER_S;

  return (
    <Layout>
      <div className="space-y-8">
        <div>
          <div className="flex items-center gap-3 mb-4">
            <Lightbulb className="w-8 h-8 text-amber-400" />
            <div>
              <h1 className="text-3xl font-bold font-display text-white">Model Proposal</h1>
              <p className="text-slate-400 text-sm">A Self-Regulated Cosmological Origin for the Galactic Acceleration Scale</p>
            </div>
          </div>

          <div className="bg-gradient-to-r from-amber-500/10 to-violet-500/10 border border-amber-500/20 rounded-xl p-6 mb-2">
            <p className="text-amber-300 font-mono text-sm leading-relaxed text-center italic">
              "Galaxies do not choose a{"\u2080"} randomly; they evolve toward a self-regulated state where the transition between the baryon regime and the halo regime occurs near a cosmological scale of order cH{"\u2080"}."
            </p>
          </div>
        </div>

        <GlassCard glow="amber" className="border border-amber-500/20">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-amber-500 to-orange-600 flex items-center justify-center text-white font-bold text-sm">0</div>
            <div>
              <h2 className="text-xl font-bold text-white">The Core Idea</h2>
              <p className="text-xs text-slate-400">One sentence, one model</p>
            </div>
          </div>

          <div className="bg-white/5 rounded-xl p-5 mb-6">
            <p className="text-sm text-slate-200 leading-relaxed">
              We propose that the acceleration scale <span className="text-cyan-400 font-mono font-bold">a{"\u2080"} {"\u2248"} 1.2 {"\u00D7"} 10{"\u207B\u00B9\u2070"} m/s{"\u00B2"}</span> in
              galaxy dynamics is an <span className="text-amber-400 font-bold">emergent</span> quantity arising from
              self-regulation between baryons and halos, <span className="text-violet-400 font-bold">anchored</span> to a
              cosmological scale of order <span className="text-cyan-400 font-mono">cH{"\u2080"}</span>, with small deviations depending on baryonic structure.
            </p>
          </div>

          <div className="bg-gradient-to-r from-cyan-500/5 to-violet-500/5 border border-cyan-500/20 rounded-xl p-5">
            <h3 className="text-cyan-300 font-bold text-sm mb-3">The Equation</h3>
            <div className="text-center py-4">
              <div className="text-2xl font-mono text-white mb-2">
                a{"\u2080"} = {"\u03B5"} {"\u00B7"} c {"\u00B7"} H{"\u2080"} {"\u00B7"} {"\u03A6"}(feedback, f_b, c_halo, {"\u03A3"}_bar)
              </div>
              <div className="flex justify-center gap-8 mt-4 text-xs text-slate-400">
                <span><span className="text-amber-400 font-mono">{"\u03B5"}</span> ~ O(1) coupling constant</span>
                <span><span className="text-violet-400 font-mono">{"\u03A6"}</span> {"\u2248"} 1 for most galaxies</span>
                <span><span className="text-cyan-400 font-mono">cH{"\u2080"}</span> = cosmological anchor</span>
              </div>
            </div>
            <div className="grid grid-cols-3 gap-4 mt-4">
              <div className="bg-white/5 rounded-lg p-3 text-center">
                <div className="text-xs text-slate-400 mb-1">cH{"\u2080"}</div>
                <div className="text-lg font-bold text-cyan-400 font-mono">{cH0_ms2.toExponential(2)}</div>
                <div className="text-xs text-slate-400">m/s{"\u00B2"}</div>
              </div>
              <div className="bg-white/5 rounded-lg p-3 text-center">
                <div className="text-xs text-slate-400 mb-1">a{"\u2080"} (observed)</div>
                <div className="text-lg font-bold text-amber-400 font-mono">{A0_MS2.toExponential(2)}</div>
                <div className="text-xs text-slate-400">m/s{"\u00B2"}</div>
              </div>
              <div className="bg-white/5 rounded-lg p-3 text-center">
                <div className="text-xs text-slate-400 mb-1">{"\u03B5"} = a{"\u2080"} / cH{"\u2080"}</div>
                <div className="text-lg font-bold text-violet-400 font-mono">{epsilon.toFixed(2)}</div>
                <div className="text-xs text-slate-400">~ O(1) {"\u2714"}</div>
              </div>
            </div>
          </div>
        </GlassCard>

        <GlassCard glow="cyan" className="border border-cyan-500/20">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center text-white font-bold text-sm">I</div>
            <div>
              <h2 className="text-xl font-bold text-white">What the Data Says</h2>
              <p className="text-xs text-slate-400">5 hard constraints any model must satisfy</p>
            </div>
          </div>

          <div className="space-y-3">
            {[
              {
                icon: CheckCircle2,
                color: "emerald",
                title: "Clear transition scale at a\u2080 \u2248 1.2 \u00D7 10\u207B\u00B9\u2070 m/s\u00B2",
                detail: transition ? transition.nPoints + " points across " + transition.nGalaxies + " galaxies collapse onto one curve. RMS scatter: " + transition.collapse.rmsWithCorrectA0 + " dex." : "Loading..."
              },
              {
                icon: CheckCircle2,
                color: "emerald",
                title: "RAR is nearly univariate",
                detail: bivariate?.honestControl ? "Adding \u03A3_bar gives only " + bivariate.honestControl.marginalSigmaBar.crossValidatedImprovement.toFixed(2) + "% improvement. F=" + bivariate.honestControl.marginalSigmaBar.fStatistic.toFixed(1) + ": detectable but marginal." : "g_obs is well-described by g_bar alone; \u03A3_bar is a second-order correction, not a new law."
              },
              {
                icon: CheckCircle2,
                color: "emerald",
                title: "Baryons are required",
                detail: "DMO (dark-matter-only) simulations do NOT produce the RAR or a well-defined transition. Hydrodynamic sims with baryonic physics do."
              },
              {
                icon: AlertTriangle,
                color: "amber",
                title: "a\u2080 is NOT robustly predicted",
                detail: simA0 ? "Different simulations give different a\u2080: EAGLE \u2248 1.2, TNG \u2248 1.0\u20131.5, NIHAO \u2248 0.8\u20131.5 (\u00D710\u207B\u00B9\u2070 m/s\u00B2). Value depends on feedback recipe." : "Loading..."
              },
              {
                icon: AlertTriangle,
                color: "amber",
                title: "Observations are tighter than simulations",
                detail: "Observed scatter (0.057 dex, McGaugh+2016) is TIGHTER than most simulation predictions (0.08\u20130.17 dex). Reality is 'cleaner' than expected."
              }
            ].map((item, i) => (
              <div key={i} className={"flex gap-3 p-4 rounded-xl border " + (item.color === "emerald" ? "bg-emerald-500/5 border-emerald-500/20" : "bg-amber-500/5 border-amber-500/20")}>
                <item.icon className={"w-5 h-5 mt-0.5 flex-shrink-0 " + (item.color === "emerald" ? "text-emerald-400" : "text-amber-400")} />
                <div>
                  <p className={"text-sm font-bold " + (item.color === "emerald" ? "text-emerald-300" : "text-amber-300")}>{item.title}</p>
                  <p className="text-xs text-slate-400 mt-1">{item.detail}</p>
                </div>
              </div>
            ))}
          </div>
        </GlassCard>

        <GlassCard glow="violet" className="border border-violet-500/20">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-violet-500 to-purple-600 flex items-center justify-center text-white font-bold text-sm">II</div>
            <div>
              <h2 className="text-xl font-bold text-white">The Physical Mechanism</h2>
              <p className="text-xs text-slate-400">Self-regulated feedback loop</p>
            </div>
          </div>

          <div className="bg-gradient-to-b from-violet-500/5 to-transparent border border-violet-500/20 rounded-xl p-6 mb-6">
            <h3 className="text-violet-300 font-bold text-sm mb-6 text-center">The Feedback Loop That Produces a{"\u2080"}</h3>

            <div className="flex flex-col items-center gap-2">
              <div className="bg-cyan-500/10 border border-cyan-500/30 rounded-xl px-6 py-3 text-center">
                <div className="text-xs text-slate-400">Step 1</div>
                <div className="text-sm font-bold text-cyan-300">Cosmological Background</div>
                <div className="text-xs text-slate-400 mt-1">Sets the scale: cH{"\u2080"} ~ 7 {"\u00D7"} 10{"\u207B\u00B9\u2070"} m/s{"\u00B2"}</div>
              </div>
              <ArrowRight className="w-4 h-4 text-slate-500 rotate-90" />
              <div className="bg-emerald-500/10 border border-emerald-500/30 rounded-xl px-6 py-3 text-center">
                <div className="text-xs text-slate-400">Step 2</div>
                <div className="text-sm font-bold text-emerald-300">Baryonic Cooling / Accretion</div>
                <div className="text-xs text-slate-400 mt-1">Gas cools, collapses, increases g_bar locally</div>
              </div>
              <ArrowRight className="w-4 h-4 text-slate-500 rotate-90" />
              <div className="bg-amber-500/10 border border-amber-500/30 rounded-xl px-6 py-3 text-center">
                <div className="text-xs text-slate-400">Step 3</div>
                <div className="text-sm font-bold text-amber-300">Feedback (SN / AGN)</div>
                <div className="text-xs text-slate-400 mt-1">Redistributes gas, modifies gravitational potential</div>
              </div>
              <ArrowRight className="w-4 h-4 text-slate-500 rotate-90" />
              <div className="bg-violet-500/10 border border-violet-500/30 rounded-xl px-6 py-3 text-center">
                <div className="text-xs text-slate-400">Step 4</div>
                <div className="text-sm font-bold text-violet-300">Halo Response</div>
                <div className="text-xs text-slate-400 mt-1">Dark matter halo adjusts density profile</div>
              </div>
              <ArrowRight className="w-4 h-4 text-slate-500 rotate-90" />
              <div className="bg-rose-500/10 border border-rose-500/30 rounded-xl px-6 py-3 text-center w-full max-w-md">
                <div className="text-xs text-slate-400">Equilibrium</div>
                <div className="text-lg font-bold text-rose-300">a{"\u2080"} = {"\u03B5"} {"\u00B7"} cH{"\u2080"} {"\u00B7"} {"\u03A6"}</div>
                <div className="text-xs text-slate-400 mt-1">Transition occurs near cosmological scale</div>
              </div>
            </div>

            <div className="mt-6 bg-white/5 rounded-lg p-4">
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4 text-xs">
                <div>
                  <p className="text-cyan-400 font-bold mb-1">When g_bar {">"} a{"\u2080"}:</p>
                  <p className="text-slate-300">Baryons dominate the potential. Newton works. The halo is subdominant.</p>
                </div>
                <div>
                  <p className="text-amber-400 font-bold mb-1">When g_bar {"<"} a{"\u2080"}:</p>
                  <p className="text-slate-300">Halo + feedback history dominate. g_obs/g_bar {">>"} 1. The DM "knows" about the baryons because they co-evolved.</p>
                </div>
              </div>
            </div>
          </div>

          <div className="bg-white/5 rounded-xl p-5">
            <h3 className="text-white font-bold text-sm mb-3">Why This Transition is NOT Random</h3>
            <p className="text-xs text-slate-300 leading-relaxed mb-3">
              The transition happens when two competing effects balance:
            </p>
            <div className="bg-gradient-to-r from-cyan-500/5 to-amber-500/5 rounded-lg p-4 text-center">
              <div className="text-lg font-mono text-white mb-2">
                g_bar ~ g_crit {"\u221D"} cH{"\u2080"}
              </div>
              <p className="text-xs text-slate-400">
                The critical acceleration is set by the Hubble flow, not by any local galaxy property.
                This is why a{"\u2080"} is nearly universal across all galaxy types.
              </p>
            </div>
          </div>
        </GlassCard>

        <GlassCard glow="emerald" className="border border-emerald-500/20">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-emerald-500 to-teal-600 flex items-center justify-center text-white font-bold text-sm">III</div>
            <div>
              <h2 className="text-xl font-bold text-white">Testable Predictions</h2>
              <p className="text-xs text-slate-400">What this model predicts — and what we already confirmed</p>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            {[
              {
                prediction: "Changing feedback recipe changes a\u2080 and/or scatter",
                status: "confirmed",
                evidence: "NIHAO explicitly shows this. TNG and EAGLE give different scatter (0.13-0.17 vs 0.08-0.12 dex).",
                color: "emerald"
              },
              {
                prediction: "Observed scatter < simulation scatter",
                status: "confirmed",
                evidence: "McGaugh+2016: 0.057 dex observed. Most sims predict 0.08-0.17 dex. Real galaxies reach tighter self-regulation.",
                color: "emerald"
              },
              {
                prediction: "\u03A3_bar dependence is a small correction, not a new law",
                status: "confirmed",
                evidence: "Our bivariate test: F=30.6 (statistically significant) but only 0.65% cross-validated improvement. Second-order effect.",
                color: "emerald"
              },
              {
                prediction: "DMO simulations fail completely",
                status: "confirmed",
                evidence: "DMO produces no RAR-like relation and scatter >0.3 dex. Baryons are necessary for the feedback loop.",
                color: "emerald"
              },
              {
                prediction: "Same pipeline on TNG/EAGLE/FIRE with unified definitions: a\u2080 varies across suites",
                status: "testable",
                evidence: "If a\u2080 remains non-robust across simulation suites, it is a diagnostic for feedback — confirming the model.",
                color: "amber"
              },
              {
                prediction: "Galaxies with more violent feedback show larger RAR scatter",
                status: "testable",
                evidence: "Starburst galaxies and dwarf irregulars should show wider scatter around the mean RAR than quiescent spirals.",
                color: "amber"
              }
            ].map((p, i) => (
              <div key={i} className={"rounded-xl border p-4 " + (p.color === "emerald" ? "bg-emerald-500/5 border-emerald-500/20" : "bg-amber-500/5 border-amber-500/20")}>
                <div className="flex items-center gap-2 mb-2">
                  {p.status === "confirmed" ?
                    <CheckCircle2 className="w-4 h-4 text-emerald-400" /> :
                    <FlaskConical className="w-4 h-4 text-amber-400" />}
                  <span className={"text-xs px-2 py-0.5 rounded-full " + (p.status === "confirmed" ? "bg-emerald-500/20 text-emerald-400" : "bg-amber-500/20 text-amber-400")}>
                    {p.status === "confirmed" ? "CONFIRMED" : "TESTABLE"}
                  </span>
                </div>
                <p className="text-sm font-bold text-white mb-1">{p.prediction}</p>
                <p className="text-xs text-slate-400">{p.evidence}</p>
              </div>
            ))}
          </div>
        </GlassCard>

        {transition && (
          <GlassCard glow="cyan" className="border border-cyan-500/20">
            <div className="flex items-center gap-3 mb-6">
              <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center text-white font-bold text-sm">IV</div>
              <div>
                <h2 className="text-xl font-bold text-white">The Transition Plot</h2>
                <p className="text-xs text-slate-400">g_obs/g_bar vs g_bar: where the model lives</p>
              </div>
            </div>

            <div className="h-72">
              <ResponsiveContainer width="100%" height="100%">
                <BarChart data={transition.ratioBins.filter(b => b.n >= 10)} margin={{ top: 5, right: 20, bottom: 25, left: 20 }}>
                  <XAxis dataKey="logGbar" tick={{ fill: '#94a3b8', fontSize: 10 }} label={{ value: 'log(g_bar) [(km/s)\u00B2/kpc]', position: 'bottom', fill: '#94a3b8', fontSize: 11, dy: 10 }} />
                  <YAxis tick={{ fill: '#94a3b8', fontSize: 10 }} label={{ value: 'g_obs / g_bar (median)', angle: -90, position: 'insideLeft', fill: '#94a3b8', fontSize: 11, dx: -5 }} domain={[0, 'auto']} />
                  <Tooltip contentStyle={{ backgroundColor: '#1e293b', border: '1px solid #334155', borderRadius: '8px', fontSize: '11px' }} formatter={(v: number) => [v.toFixed(2), 'g_obs/g_bar']} />
                  <ReferenceLine y={1} stroke="#f59e0b" strokeDasharray="5 5" label={{ value: 'Newton', fill: '#f59e0b', fontSize: 10, position: 'right' }} />
                  <ReferenceLine x={Math.log10(transition.a0_corrected)} stroke="#ef4444" strokeDasharray="3 3" label={{ value: 'a\u2080', fill: '#ef4444', fontSize: 11, position: 'top' }} />
                  <Bar dataKey="medianRatio" name="g_obs/g_bar">
                    {transition.ratioBins.filter(b => b.n >= 10).map((b, i) => (
                      <Cell key={i} fill={b.logGbar < Math.log10(transition.a0_corrected) ? '#06b6d4' : b.logGbar < Math.log10(transition.a0_corrected) + 0.5 ? '#8b5cf6' : '#10b981'} fillOpacity={0.7} />
                    ))}
                  </Bar>
                </BarChart>
              </ResponsiveContainer>
            </div>
            <div className="flex justify-center gap-6 mt-2 text-xs">
              <span className="flex items-center gap-1"><span className="w-3 h-3 rounded bg-cyan-500 inline-block"></span> DM-dominated ({"\u03A6"} active)</span>
              <span className="flex items-center gap-1"><span className="w-3 h-3 rounded bg-violet-500 inline-block"></span> Transition zone (a{"\u2080"})</span>
              <span className="flex items-center gap-1"><span className="w-3 h-3 rounded bg-emerald-500 inline-block"></span> Newtonian ({"\u03A6"} {"\u2248"} 1)</span>
            </div>

            <div className="mt-6 bg-white/5 rounded-xl p-4">
              <h4 className="text-sm font-bold text-white mb-2">Model Interpretation</h4>
              <p className="text-xs text-slate-300 leading-relaxed">
                In the model equation <span className="text-amber-400 font-mono">a{"\u2080"} = {"\u03B5"}{"\u00B7"}cH{"\u2080"}{"\u00B7"}{"\u03A6"}</span>,
                the function {"\u03A6"} encodes how the baryon-halo feedback loop modifies the transition point.
                For most galaxies, {"\u03A6"} {"\u2248"} 1 (hence the tight RAR). The small deviations of {"\u03A6"} from unity
                produce the weak {"\u03A3"}_bar dependence we measured (0.65% improvement). The
                fact that {"\u03B5"} {"\u2248"} {epsilon.toFixed(2)} {"\u2248"} 1/(2{"\u03C0"}) connects galactic dynamics
                to the Hubble expansion rate {"—"} a coincidence that current {"\u039B"}CDM does not predict but our model
                makes structural.
              </p>
            </div>
          </GlassCard>
        )}

        <GlassCard glow="rose" className="border border-rose-500/20">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-rose-500 to-red-600 flex items-center justify-center text-white font-bold text-sm">V</div>
            <div>
              <h2 className="text-xl font-bold text-white">What We Do NOT Claim</h2>
              <p className="text-xs text-slate-400">Intellectual honesty</p>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <XCircle className="w-6 h-6 text-rose-400 mb-2" />
              <p className="text-sm font-bold text-rose-300 mb-1">NOT new gravity</p>
              <p className="text-xs text-slate-400">We do not propose a modification of general relativity. This is within standard physics.</p>
            </div>
            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <XCircle className="w-6 h-6 text-rose-400 mb-2" />
              <p className="text-sm font-bold text-rose-300 mb-1">NOT proof of MOND</p>
              <p className="text-xs text-slate-400">MOND gets the phenomenology right, but we do not claim it is the correct fundamental theory.</p>
            </div>
            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <XCircle className="w-6 h-6 text-rose-400 mb-2" />
              <p className="text-sm font-bold text-rose-300 mb-1">NOT falsifying {"\u039B"}CDM</p>
              <p className="text-xs text-slate-400">{"\u039B"}CDM sims can approximate the RAR. The tension is quantitative (scatter, universality), not qualitative.</p>
            </div>
          </div>

          <div className="mt-4 bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-4">
            <CheckCircle2 className="w-6 h-6 text-emerald-400 mb-2" />
            <p className="text-sm font-bold text-emerald-300 mb-1">What we DO claim</p>
            <p className="text-xs text-slate-300 leading-relaxed">
              There exists an observationally robust, emergent acceleration scale linked to cosmological parameters (a{"\u2080"} ~ cH{"\u2080"}/2{"\u03C0"}).
              Current {"\u039B"}CDM simulations reproduce it approximately but not robustly {"—"} the exact value
              and scatter depend on subgrid feedback physics. The observed scatter is tighter than predictions,
              suggesting real galaxies achieve a deeper level of self-regulation than current simulations capture.
              This makes a{"\u2080"} both a constraint on galaxy formation models and a potential diagnostic for feedback prescriptions.
            </p>
          </div>
        </GlassCard>

        {simA0 && (
          <GlassCard glow="violet" className="border border-violet-500/20">
            <div className="flex items-center gap-3 mb-6">
              <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-violet-500 to-purple-600 flex items-center justify-center text-white font-bold text-sm">VI</div>
              <div>
                <h2 className="text-xl font-bold text-white">Simulation Scorecard</h2>
                <p className="text-xs text-slate-400">How well does each simulation reproduce the model?</p>
              </div>
            </div>

            <div className="space-y-3">
              {simA0.simulations.map((sim, i) => {
                const matchesA0 = sim.a0_kpc !== null && Math.abs(sim.a0_kpc - 3702) < 1000;
                const score = (sim.universal ? 1 : 0) + (matchesA0 ? 1 : 0) + (sim.a0_kpc !== null ? 1 : 0);
                const scoreColor = score >= 3 ? "emerald" : score >= 2 ? "amber" : "rose";
                return (
                  <div key={i} className={"rounded-xl border p-4 " + ("bg-" + scoreColor + "-500/5 border-" + scoreColor + "-500/20")}>
                    <div className="flex items-center justify-between mb-2">
                      <div>
                        <span className="text-sm font-bold text-white">{sim.name}</span>
                        <span className="text-xs text-slate-500 ml-2">{sim.ref}</span>
                      </div>
                      <div className="flex gap-1">
                        {[0, 1, 2].map(s => (
                          <div key={s} className={"w-3 h-3 rounded-full " + (s < score ? "bg-" + scoreColor + "-400" : "bg-slate-700")} />
                        ))}
                      </div>
                    </div>
                    <div className="grid grid-cols-4 gap-2 text-xs">
                      <div><span className="text-slate-400">a{"\u2080"}:</span> <span className="text-white">{sim.a0_reported}</span></div>
                      <div><span className="text-slate-400">Scatter:</span> <span className="text-cyan-400">{sim.scatter}</span></div>
                      <div><span className="text-slate-400">Universal:</span> <span className={sim.universal ? "text-emerald-400" : "text-rose-400"}>{sim.universal ? "Yes" : "No"}</span></div>
                      <div><span className="text-slate-400">Feedback:</span> <span className="text-amber-400">{sim.feedbackDependent}</span></div>
                    </div>
                    <p className="text-xs text-slate-500 mt-1 italic">{sim.notes}</p>
                  </div>
                );
              })}
            </div>
          </GlassCard>
        )}

        <GlassCard glow="amber" className="border-2 border-amber-500/30">
          <div className="text-center py-4">
            <Lightbulb className="w-10 h-10 text-amber-400 mx-auto mb-3" />
            <h2 className="text-2xl font-bold text-white mb-4">Summary</h2>

            <div className="max-w-3xl mx-auto space-y-4">
              <div className="bg-gradient-to-r from-cyan-500/10 to-violet-500/10 border border-cyan-500/20 rounded-xl p-5">
                <p className="text-cyan-300 font-mono text-sm leading-relaxed mb-3">
                  a{"\u2080"} appears to be an emergent scale, anchored to the cosmological background (cH{"\u2080"}),
                  produced by the self-regulated interaction of baryons and dark matter halos under feedback.
                </p>
                <p className="text-slate-300 text-xs leading-relaxed">
                  The RAR remains the primary law. {"\u03A3"}_bar is a second-order correction.
                  The real question is not whether a{"\u2080"} exists {"—"} it does {"—"} but why its value
                  and scatter are more stable in observations than in any current simulation.
                </p>
              </div>

              <div className="bg-amber-500/10 border border-amber-500/20 rounded-xl p-5">
                <h3 className="text-amber-300 font-bold text-sm mb-2">The Minimal Equation</h3>
                <div className="text-center">
                  <div className="text-xl font-mono text-white mb-2">
                    a{"\u2080"} = {"\u03B5"} {"\u00B7"} cH{"\u2080"} {"\u00B7"} (1 + {"\u03B4"})
                  </div>
                  <div className="text-xs text-slate-400">
                    where {"\u03B5"} {"\u2248"} {epsilon.toFixed(2)} and {"\u03B4"} ~ 0 encodes residual baryonic modulation
                  </div>
                </div>
              </div>

              <div className="text-xs text-slate-400 leading-relaxed">
                <p className="mb-2"><span className="text-white font-bold">Not claimed:</span> New gravity, MOND, or falsification of {"\u039B"}CDM.</p>
                <p><span className="text-amber-400 font-bold">Claimed:</span> An observationally robust emergent scale, not yet robustly predicted by simulations, that every theory of galaxy formation must account for.</p>
              </div>
            </div>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
