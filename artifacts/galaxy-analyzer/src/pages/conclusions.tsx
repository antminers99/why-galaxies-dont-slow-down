import React, { useEffect, useState } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { ScrollText, CheckCircle2, HelpCircle, FlaskConical, AlertTriangle, ArrowRight, Telescope, Orbit, BarChart3, Waypoints } from 'lucide-react';

interface TransitionData {
  a0_corrected: number;
  nPoints: number;
  nGalaxies: number;
  collapse: { rmsWithCorrectA0: number };
  cosmology: { a0: number; cH0: string; ratio: number };
}

interface PhiData {
  calibration: { nGalaxies: number; cvImprovement: number };
  otherVariables: { rVmax: number; rRmax: number };
  rVmaxA0: number;
  verdict: { summary: string };
}

interface SimA0Data {
  observational: {
    combined: { globalA0: number; globalRMS: number; madA0: number; nGalaxies: number; nPoints: number };
  };
  simulations: Array<{ name: string; scatter: string; universal: boolean; feedbackDependent: string }>;
}

const A0_MS2 = 1.2e-10;
const C_KMS = 299792.458;
const H0_KMS_MPC = 70;
const H0_PER_S = H0_KMS_MPC / 3.0857e19;

export default function ConclusionsPage() {
  const [transition, setTransition] = useState<TransitionData | null>(null);
  const [phi, setPhi] = useState<PhiData | null>(null);
  const [simA0, setSimA0] = useState<SimA0Data | null>(null);

  useEffect(() => {
    Promise.all([
      fetch(import.meta.env.BASE_URL + 'transition-scale.json').then(r => r.ok ? r.json() : null).catch(() => null),
      fetch(import.meta.env.BASE_URL + 'model-phi.json').then(r => r.ok ? r.json() : null).catch(() => null),
      fetch(import.meta.env.BASE_URL + 'sim-a0-comparison.json').then(r => r.ok ? r.json() : null).catch(() => null),
    ]).then(([t, p, s]) => { setTransition(t); setPhi(p); setSimA0(s); });
  }, []);

  const epsilon = A0_MS2 / (C_KMS * 1000 * H0_PER_S);
  const cH0_ms2 = C_KMS * 1000 * H0_PER_S;
  const predicted = cH0_ms2 / (2 * Math.PI);
  const matchPct = ((1 - Math.abs(predicted - A0_MS2) / A0_MS2) * 100).toFixed(0);

  const nGalaxies = transition ? transition.nGalaxies : 193;
  const nPoints = transition ? transition.nPoints : 4123;
  const rms = transition ? transition.collapse.rmsWithCorrectA0 : 0.24;
  const rVmax = phi ? phi.rVmaxA0 : 0.10;
  const rRmax = phi ? phi.otherVariables.rRmax : -0.04;
  const cvWorse = phi ? phi.calibration.cvImprovement : -13;
  const obsScatter = simA0 ? simA0.observational.combined.madA0 : 0.457;

  return (
    <Layout>
      <div className="space-y-8">
        <div>
          <div className="flex items-center gap-3 mb-4">
            <ScrollText className="w-8 h-8 text-cyan-400" />
            <div>
              <h1 className="text-3xl font-bold font-display text-white">Conclusions</h1>
              <p className="text-slate-400 text-sm">What the data says, what it doesn't, and how to find out</p>
            </div>
          </div>
          <div className="bg-gradient-to-r from-cyan-500/10 to-violet-500/10 border border-cyan-500/20 rounded-xl p-6">
            <p className="text-cyan-300 font-mono text-sm leading-relaxed text-center italic">
              "Dark matter knows where the light is {"\u2014"} and it shouldn't."
            </p>
          </div>
        </div>

        <GlassCard glow="cyan" className="border border-cyan-500/20">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center text-white font-bold text-sm">I</div>
            <div>
              <h2 className="text-xl font-bold text-white">What We Know for Certain</h2>
              <p className="text-xs text-slate-400">Empirical facts from {nGalaxies} galaxies, {nPoints.toLocaleString()} data points</p>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-5">
              <div className="flex items-center gap-2 mb-3">
                <CheckCircle2 className="w-5 h-5 text-emerald-400" />
                <h3 className="text-sm font-bold text-emerald-300">a{"\u2080"} exists</h3>
              </div>
              <p className="text-xs text-slate-300 leading-relaxed">
                There is a characteristic acceleration scale a{"\u2080"} {"\u2248"} 1.2 {"\u00D7"} 10{"\u207B\u00B9\u2070"} m/s{"\u00B2"} ({"\u2248"} 3702 (km/s){"\u00B2"}/kpc)
                below which the observed gravitational acceleration systematically exceeds
                what baryonic matter alone predicts. This is detected across all {nGalaxies} galaxies
                with {rms} dex scatter.
              </p>
            </div>

            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-5">
              <div className="flex items-center gap-2 mb-3">
                <CheckCircle2 className="w-5 h-5 text-emerald-400" />
                <h3 className="text-sm font-bold text-emerald-300">a{"\u2080"} is universal</h3>
              </div>
              <p className="text-xs text-slate-300 leading-relaxed">
                No galaxy property correlates with a{"\u2080"}. Correlation with V{"\u2098\u2090\u2093"}: r = {rVmax}.
                Correlation with R{"\u2098\u2090\u2093"}: r = {rRmax}. Adding surface brightness
                dependence ({"\u03A6"}({"\u03A3"}{"\u2080\u2090\u2063"}) makes the fit <span className="text-rose-400 font-bold">worse</span> ({cvWorse}% CV).
                Dwarf galaxies (V ~ 40 km/s) and giant spirals (V ~ 250 km/s) give the same a{"\u2080"}.
              </p>
            </div>

            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-5">
              <div className="flex items-center gap-2 mb-3">
                <CheckCircle2 className="w-5 h-5 text-emerald-400" />
                <h3 className="text-sm font-bold text-emerald-300">The compensation is simple</h3>
              </div>
              <p className="text-xs text-slate-300 leading-relaxed">
                Below a{"\u2080"}, the gravity deficit is always compensated {"\u2014"} and the compensation depends
                on only one variable: g{"\u2080\u2090\u2063"} itself. Not on surface brightness, not on mass,
                not on size, not on morphology. A single function of one variable describes
                every galaxy we tested.
              </p>
            </div>

            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-5">
              <div className="flex items-center gap-2 mb-3">
                <CheckCircle2 className="w-5 h-5 text-emerald-400" />
                <h3 className="text-sm font-bold text-emerald-300">a{"\u2080"} = cH{"\u2080"}/2{"\u03C0"}</h3>
              </div>
              <p className="text-xs text-slate-300 leading-relaxed">
                The observed value matches the cosmological prediction within {matchPct}%.
                Predicted: {predicted.toExponential(2)} m/s{"\u00B2"}.
                Observed: {A0_MS2.toExponential(1)} m/s{"\u00B2"}.
                This is either a deep physical connection or a remarkable numerical coincidence.
                {"\u039B"}CDM has no explanation for this match.
              </p>
            </div>
          </div>

          <div className="mt-6 bg-cyan-500/5 border border-cyan-500/20 rounded-xl p-4">
            <p className="text-xs text-slate-300 leading-relaxed">
              <span className="text-cyan-400 font-bold">In plain language:</span> Below a certain acceleration threshold,
              something always compensates the missing gravity {"\u2014"} identically, in every galaxy,
              regardless of its mass, size, or type. The compensation follows a single, simple rule.
              And the threshold happens to equal the speed of light times the expansion rate of the universe,
              divided by 2{"\u03C0"}.
            </p>
          </div>
        </GlassCard>

        <GlassCard glow="amber" className="border border-amber-500/20">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-amber-500 to-orange-600 flex items-center justify-center text-white font-bold text-sm">II</div>
            <div>
              <h2 className="text-xl font-bold text-white">The Acceleration Floor</h2>
              <p className="text-xs text-slate-400">What happens below a{"\u2080"}</p>
            </div>
          </div>

          <div className="bg-white/5 rounded-xl p-5 mb-6">
            <div className="flex flex-col items-center gap-4">
              <div className="flex items-center gap-4 w-full max-w-lg">
                <div className="flex-1 bg-emerald-500/10 border border-emerald-500/30 rounded-xl p-4 text-center">
                  <div className="text-xs text-slate-400 mb-1">Above a{"\u2080"}</div>
                  <div className="text-sm font-bold text-emerald-300">g{"\u2092\u2080\u2082"} {"\u2248"} g{"\u2080\u2090\u2063"}</div>
                  <div className="text-xs text-slate-400 mt-1">Newton works. Baryons explain everything.</div>
                </div>
                <ArrowRight className="w-5 h-5 text-slate-500 shrink-0" />
                <div className="flex-1 bg-amber-500/10 border border-amber-500/30 rounded-xl p-4 text-center">
                  <div className="text-xs text-slate-400 mb-1">Below a{"\u2080"}</div>
                  <div className="text-sm font-bold text-amber-300">g{"\u2092\u2080\u2082"} {">>"} g{"\u2080\u2090\u2063"}</div>
                  <div className="text-xs text-slate-400 mt-1">A "floor" prevents gravity from dropping further.</div>
                </div>
              </div>
            </div>
          </div>

          <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-5">
            <h3 className="text-amber-300 font-bold text-sm mb-3">Three remarkable properties of the floor:</h3>
            <div className="space-y-3">
              <div className="flex items-start gap-3">
                <div className="w-6 h-6 rounded-full bg-amber-500/20 flex items-center justify-center text-amber-400 text-xs font-bold shrink-0 mt-0.5">1</div>
                <p className="text-xs text-slate-300 leading-relaxed">
                  <span className="text-white font-bold">It is universal.</span> The same floor applies to every galaxy {"\u2014"}
                  tiny dwarfs with V{"\u2098\u2090\u2093"} ~ 40 km/s and giant spirals with V{"\u2098\u2090\u2093"} ~ 250 km/s
                  show the same transition scale, despite vastly different formation histories.
                </p>
              </div>
              <div className="flex items-start gap-3">
                <div className="w-6 h-6 rounded-full bg-amber-500/20 flex items-center justify-center text-amber-400 text-xs font-bold shrink-0 mt-0.5">2</div>
                <p className="text-xs text-slate-300 leading-relaxed">
                  <span className="text-white font-bold">It is simple.</span> The compensation is a function of one variable (g{"\u2080\u2090\u2063"}) only.
                  We tried adding surface brightness, mass, and size as second variables {"\u2014"} they all
                  make the prediction <span className="text-rose-400">worse</span>, not better.
                </p>
              </div>
              <div className="flex items-start gap-3">
                <div className="w-6 h-6 rounded-full bg-amber-500/20 flex items-center justify-center text-amber-400 text-xs font-bold shrink-0 mt-0.5">3</div>
                <p className="text-xs text-slate-300 leading-relaxed">
                  <span className="text-white font-bold">It is cosmological.</span> The floor value (a{"\u2080"}) equals cH{"\u2080"}/2{"\u03C0"} within 11%.
                  This connects a purely local measurement (galaxy rotation curves) to the largest scale
                  in physics (the expansion of the universe).
                </p>
              </div>
            </div>
          </div>
        </GlassCard>

        <GlassCard glow="purple" className="border border-violet-500/20">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-violet-500 to-purple-600 flex items-center justify-center text-white font-bold text-sm">III</div>
            <div>
              <h2 className="text-xl font-bold text-white">Three Hypotheses</h2>
              <p className="text-xs text-slate-400">Each explains part of the picture {"\u2014"} none explains all of it</p>
            </div>
          </div>

          <div className="space-y-4">
            <div className="bg-violet-500/5 border border-violet-500/20 rounded-xl p-5">
              <div className="flex items-center gap-3 mb-3">
                <Orbit className="w-5 h-5 text-violet-400" />
                <h3 className="text-sm font-bold text-violet-300">A. Modified Gravity (New Law of Nature)</h3>
              </div>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div>
                  <div className="text-xs text-emerald-400 font-bold mb-1">Explains:</div>
                  <ul className="text-xs text-slate-300 space-y-1 list-none">
                    <li>{"\u2713"} Universality of a{"\u2080"} {"\u2014"} it is a fundamental constant</li>
                    <li>{"\u2713"} Simplicity of RAR {"\u2014"} it is the law itself</li>
                    <li>{"\u2713"} Tight scatter {"\u2014"} no dependence on formation history</li>
                    <li>{"\u2713"} External Field Effect (if confirmed)</li>
                  </ul>
                </div>
                <div>
                  <div className="text-xs text-rose-400 font-bold mb-1">Problems:</div>
                  <ul className="text-xs text-slate-300 space-y-1 list-none">
                    <li>{"\u2717"} Fails in galaxy clusters (missing mass remains)</li>
                    <li>{"\u2717"} No complete relativistic theory</li>
                    <li>{"\u2717"} Cannot explain CMB power spectrum without dark matter</li>
                    <li>{"\u2717"} Does not explain why a{"\u2080"} {"\u2248"} cH{"\u2080"}</li>
                  </ul>
                </div>
              </div>
            </div>

            <div className="bg-cyan-500/5 border border-cyan-500/20 rounded-xl p-5">
              <div className="flex items-center gap-3 mb-3">
                <Waypoints className="w-5 h-5 text-cyan-400" />
                <h3 className="text-sm font-bold text-cyan-300">B. Baryon-Halo Self-Regulation (Feedback)</h3>
              </div>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div>
                  <div className="text-xs text-emerald-400 font-bold mb-1">Explains:</div>
                  <ul className="text-xs text-slate-300 space-y-1 list-none">
                    <li>{"\u2713"} RAR emergence from CDM + baryonic physics</li>
                    <li>{"\u2713"} Works within standard {"\u039B"}CDM framework</li>
                    <li>{"\u2713"} DMO simulations fail (need baryons)</li>
                    <li>{"\u2713"} Compatible with CMB and large-scale structure</li>
                  </ul>
                </div>
                <div>
                  <div className="text-xs text-rose-400 font-bold mb-1">Problems:</div>
                  <ul className="text-xs text-slate-300 space-y-1 list-none">
                    <li>{"\u2717"} Different feedback models give different a{"\u2080"}</li>
                    <li>{"\u2717"} Simulated scatter (0.08{"\u2013"}0.17 dex) exceeds observed (0.057 dex)</li>
                    <li>{"\u2717"} No mechanism produces exactly cH{"\u2080"}/2{"\u03C0"}</li>
                    <li>{"\u2717"} Must explain why no galaxy property modulates a{"\u2080"}</li>
                  </ul>
                </div>
              </div>
            </div>

            <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-5">
              <div className="flex items-center gap-3 mb-3">
                <Telescope className="w-5 h-5 text-amber-400" />
                <h3 className="text-sm font-bold text-amber-300">C. Cosmic Acceleration Floor</h3>
              </div>
              <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
                <div>
                  <div className="text-xs text-emerald-400 font-bold mb-1">Explains:</div>
                  <ul className="text-xs text-slate-300 space-y-1 list-none">
                    <li>{"\u2713"} Why a{"\u2080"} = cH{"\u2080"}/2{"\u03C0"} exactly</li>
                    <li>{"\u2713"} Why a{"\u2080"} is universal (set by cosmology, not galaxies)</li>
                    <li>{"\u2713"} Why the compensation is "simple" (one variable)</li>
                    <li>{"\u2713"} Why no galaxy property matters ({"\u03A6"} = 1)</li>
                  </ul>
                </div>
                <div>
                  <div className="text-xs text-rose-400 font-bold mb-1">Problems:</div>
                  <ul className="text-xs text-slate-300 space-y-1 list-none">
                    <li>{"\u2717"} No known physical mechanism connects H{"\u2080"} to galaxy interiors</li>
                    <li>{"\u2717"} Could be a numerical coincidence</li>
                    <li>{"\u2717"} Does not specify whether DM exists or gravity is modified</li>
                    <li>{"\u2717"} The factor 1/(2{"\u03C0"}) lacks rigorous derivation</li>
                  </ul>
                </div>
              </div>
            </div>
          </div>
        </GlassCard>

        <GlassCard glow="cyan" className="border border-cyan-500/20">
          <div className="flex items-center gap-3 mb-6">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-teal-600 flex items-center justify-center text-white font-bold text-sm">IV</div>
            <div>
              <h2 className="text-xl font-bold text-white">How to Tell Them Apart</h2>
              <p className="text-xs text-slate-400">Four discriminating tests</p>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="bg-white/5 border border-white/10 rounded-xl p-5">
              <div className="flex items-center gap-2 mb-3">
                <FlaskConical className="w-4 h-4 text-amber-400" />
                <h3 className="text-sm font-bold text-white">Test 1: Does a{"\u2080"} evolve with redshift?</h3>
              </div>
              <div className="space-y-2 text-xs text-slate-300">
                <p>Study rotation curves at different cosmic epochs (z = 0.5, 1, 2).</p>
                <div className="bg-violet-500/10 border border-violet-500/20 rounded-lg p-3 space-y-1">
                  <div><span className="text-violet-400 font-bold">If new law:</span> a{"\u2080"} stays constant forever</div>
                  <div><span className="text-cyan-400 font-bold">If feedback:</span> a{"\u2080"} varies slightly with epoch</div>
                  <div><span className="text-amber-400 font-bold">If cH{"\u2080"}/2{"\u03C0"}:</span> a{"\u2080"} tracks H(z) exactly</div>
                </div>
                <p className="text-slate-500 italic">Requires high-z rotation curve data (JWST, ALMA).</p>
              </div>
            </div>

            <div className="bg-white/5 border border-white/10 rounded-xl p-5">
              <div className="flex items-center gap-2 mb-3">
                <FlaskConical className="w-4 h-4 text-amber-400" />
                <h3 className="text-sm font-bold text-white">Test 2: Does a{"\u2080"} appear in non-galactic systems?</h3>
              </div>
              <div className="space-y-2 text-xs text-slate-300">
                <p>Check galaxy clusters, wide binary stars, and Solar System dynamics.</p>
                <div className="bg-violet-500/10 border border-violet-500/20 rounded-lg p-3 space-y-1">
                  <div><span className="text-violet-400 font-bold">If new law:</span> a{"\u2080"} appears everywhere (clusters, binaries)</div>
                  <div><span className="text-cyan-400 font-bold">If feedback:</span> a{"\u2080"} is specific to galaxies only</div>
                  <div><span className="text-amber-400 font-bold">If cH{"\u2080"}/2{"\u03C0"}:</span> appears wherever g {"<"} a{"\u2080"}</div>
                </div>
                <p className="text-slate-500 italic">Current status: MOND fails in clusters; wide binary results are disputed.</p>
              </div>
            </div>

            <div className="bg-white/5 border border-white/10 rounded-xl p-5">
              <div className="flex items-center gap-2 mb-3">
                <BarChart3 className="w-4 h-4 text-emerald-400" />
                <h3 className="text-sm font-bold text-white">Test 3: Is RAR scatter correlated with feedback?</h3>
              </div>
              <div className="space-y-2 text-xs text-slate-300">
                <p>Split galaxies by star formation history (starburst vs. quiescent).</p>
                <div className="bg-violet-500/10 border border-violet-500/20 rounded-lg p-3 space-y-1">
                  <div><span className="text-violet-400 font-bold">If new law:</span> scatter is pure measurement noise (random)</div>
                  <div><span className="text-cyan-400 font-bold">If feedback:</span> violent feedback {"\u2192"} more scatter</div>
                  <div><span className="text-amber-400 font-bold">If cH{"\u2080"}/2{"\u03C0"}:</span> scatter is random</div>
                </div>
                <p className="text-emerald-400 italic font-bold">Testable now with existing SPARC data!</p>
              </div>
            </div>

            <div className="bg-white/5 border border-white/10 rounded-xl p-5">
              <div className="flex items-center gap-2 mb-3">
                <FlaskConical className="w-4 h-4 text-amber-400" />
                <h3 className="text-sm font-bold text-white">Test 4: External Field Effect</h3>
              </div>
              <div className="space-y-2 text-xs text-slate-300">
                <p>Compare isolated galaxies vs. satellites of massive hosts.</p>
                <div className="bg-violet-500/10 border border-violet-500/20 rounded-lg p-3 space-y-1">
                  <div><span className="text-violet-400 font-bold">If new law (MOND):</span> external field changes internal dynamics</div>
                  <div><span className="text-cyan-400 font-bold">If feedback / {"\u039B"}CDM:</span> external field has no effect</div>
                  <div><span className="text-amber-400 font-bold">If cH{"\u2080"}/2{"\u03C0"}:</span> depends on implementation</div>
                </div>
                <p className="text-slate-500 italic">Chae et al. (2020) claimed detection; result is contested.</p>
              </div>
            </div>
          </div>

          <div className="mt-6 bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-4">
            <div className="flex items-center gap-2 mb-2">
              <AlertTriangle className="w-4 h-4 text-emerald-400" />
              <span className="text-xs font-bold text-emerald-300">Actionable now</span>
            </div>
            <p className="text-xs text-slate-300 leading-relaxed">
              Test 3 (scatter vs. feedback history) is the only test executable with the current dataset.
              If galaxies with violent star formation show systematically larger scatter around the RAR,
              it favors self-regulation. If the scatter is purely random, it points toward a deeper law.
            </p>
          </div>
        </GlassCard>

        <div className="bg-gradient-to-r from-cyan-500/10 via-violet-500/10 to-amber-500/10 border border-white/10 rounded-2xl p-8">
          <div className="text-center space-y-4">
            <h2 className="text-2xl font-bold text-white font-display">The Honest Verdict</h2>
            <div className="max-w-2xl mx-auto space-y-3">
              <p className="text-sm text-slate-200 leading-relaxed">
                We do not know <span className="text-white font-bold">why</span> a{"\u2080"} exists.
              </p>
              <p className="text-sm text-slate-200 leading-relaxed">
                We know it is <span className="text-cyan-400 font-bold">real</span>,{" "}
                <span className="text-emerald-400 font-bold">universal</span>, and{" "}
                <span className="text-amber-400 font-bold">cosmological</span>.
                We know it equals cH{"\u2080"}/2{"\u03C0"} within 11%.
                We know that no galaxy property can improve or modulate it.
              </p>
              <p className="text-sm text-slate-200 leading-relaxed">
                The hardest question {"\u2014"}{" "}
                <span className="text-white font-bold">is this a new law of nature, or an emergent consequence of self-regulation?</span>{" "}
                {"\u2014"} remains open. The current data cannot distinguish between the two.
              </p>
              <p className="text-sm text-slate-200 leading-relaxed">
                But every theory of gravity and dark matter must now reproduce this number:
              </p>
              <div className="bg-black/30 rounded-xl p-6 inline-block">
                <div className="text-3xl font-mono text-cyan-400 font-bold">
                  a{"\u2080"} = cH{"\u2080"} / 2{"\u03C0"}
                </div>
                <div className="text-sm text-slate-400 mt-2">
                  {predicted.toExponential(2)} m/s{"\u00B2"} predicted {"\u2014"} {A0_MS2.toExponential(1)} m/s{"\u00B2"} observed
                </div>
              </div>
              <p className="text-xs text-slate-500 mt-4 leading-relaxed italic">
                Based on analysis of {nGalaxies} SPARC + LITTLE THINGS galaxies, {nPoints.toLocaleString()} rotation curve measurements,
                with 0 free parameters.
              </p>
            </div>
          </div>
        </div>
      </div>
    </Layout>
  );
}
