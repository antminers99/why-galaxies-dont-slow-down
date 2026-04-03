import React, { useState, useEffect, useMemo } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { Sigma, ArrowRight, CheckCircle2, Sparkles, Telescope, FlaskConical, Atom, TrendingUp } from 'lucide-react';
import {
  ScatterChart, Scatter, LineChart, Line, XAxis, YAxis, CartesianGrid,
  Tooltip, ResponsiveContainer, ReferenceLine, ComposedChart
} from 'recharts';

interface TransitionData {
  a0_ms2: number;
  a0_corrected: number;
  cosmology: { a0: number; cH0: number; ratio: number };
  plotPoints: { x: number; y: number; g: string }[];
}

interface PhiData {
  model: { epsilon: number; oneOver2Pi: number; cH0_ms2: number; cH0_kpc: number };
}

function MathBlock({ children, label }: { children: React.ReactNode; label?: string }) {
  return (
    <div className="my-4">
      {label && <div className="text-xs text-slate-500 mb-1">{label}</div>}
      <div className="bg-slate-900/80 border border-white/10 rounded-xl px-4 sm:px-6 py-4 font-mono text-sm sm:text-base text-cyan-300 overflow-x-auto">
        {children}
      </div>
    </div>
  );
}

function StepCard({ step, title, icon: Icon, children, color = "cyan" }: {
  step: number; title: string; icon: any; children: React.ReactNode; color?: string
}) {
  const colors: Record<string, string> = {
    cyan: "from-cyan-500 to-blue-600",
    purple: "from-violet-500 to-purple-600",
    amber: "from-amber-500 to-orange-600",
    emerald: "from-emerald-500 to-teal-600",
    rose: "from-rose-500 to-pink-600"
  };
  return (
    <GlassCard>
      <div className="flex items-center gap-3 mb-5">
        <div className={"w-10 h-10 rounded-xl bg-gradient-to-br " + (colors[color] || colors.cyan) + " flex items-center justify-center text-white font-bold text-sm flex-shrink-0"}>
          {step}
        </div>
        <div className="flex items-center gap-2 min-w-0">
          <Icon className="w-5 h-5 text-slate-400 flex-shrink-0" />
          <h3 className="text-lg font-bold text-white truncate">{title}</h3>
        </div>
      </div>
      {children}
    </GlassCard>
  );
}

function LimitCheck({ regime, condition, result, icon }: {
  regime: string; condition: string; result: string; icon: string
}) {
  return (
    <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-4">
      <div className="flex items-center gap-2 mb-2">
        <CheckCircle2 className="w-4 h-4 text-emerald-400" />
        <span className="text-sm font-bold text-emerald-400">{regime}</span>
        <span className="text-xs text-slate-500 ml-auto">{icon}</span>
      </div>
      <div className="text-xs text-slate-400 mb-1">
        {"When: "}{condition}
      </div>
      <div className="font-mono text-sm text-white">
        {result}
      </div>
    </div>
  );
}

export default function EquationPage() {
  const [transition, setTransition] = useState<TransitionData | null>(null);
  const [phi, setPhi] = useState<PhiData | null>(null);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + "transition-scale.json")
      .then(r => r.json()).then(setTransition).catch(() => {});
    fetch(import.meta.env.BASE_URL + "model-phi.json")
      .then(r => r.json()).then(setPhi).catch(() => {});
  }, []);

  const a0_ms2 = transition ? transition.a0_ms2 : 1.2e-10;
  const a0_kpc = transition ? transition.a0_corrected : 3702;
  const cH0 = phi ? phi.model.cH0_ms2 : 6.801e-10;
  const epsilon = phi ? phi.model.epsilon : 0.1764;
  const oneOver2Pi = phi ? phi.model.oneOver2Pi : 0.1592;

  const modelCurve = useMemo(() => {
    const points = [];
    for (let logGbar = -0.5; logGbar <= 5.5; logGbar += 0.05) {
      const gBar = Math.pow(10, logGbar);
      const logA0 = Math.log10(a0_kpc);
      const gObsFloor = Math.sqrt(gBar * gBar + a0_kpc * a0_kpc * (gBar / (gBar + a0_kpc)));
      const logGobs = Math.log10(gObsFloor);
      const gObsNewton = gBar;
      const logNewton = logGbar;
      points.push({
        logGbar,
        logGobs,
        logNewton: logGbar,
        logMOND: logGbar < logA0 ? 0.5 * logGbar + 0.5 * logA0 : logGbar
      });
    }
    return points;
  }, [a0_kpc]);

  const rarPoints = useMemo(() => {
    if (!transition?.plotPoints) return [];
    return transition.plotPoints
      .filter(p => p.x > 0 && p.y > 0)
      .map(p => ({ logGbar: p.x, logGobs: p.y }));
  }, [transition]);

  const deviationCurve = useMemo(() => {
    const points = [];
    for (let logGbar = -0.5; logGbar <= 5.0; logGbar += 0.1) {
      const gBar = Math.pow(10, logGbar);
      const ratio = a0_kpc / gBar;
      const deviation = Math.sqrt(1 + ratio * ratio * (1 / (1 + ratio))) - 1;
      points.push({ logGbar, deviation: deviation * 100 });
    }
    return points;
  }, [a0_kpc]);

  const tfPoints = useMemo(() => {
    const points = [];
    for (let logM = 7; logM <= 12; logM += 0.2) {
      const M = Math.pow(10, logM);
      const G = 4.3009e-6;
      const Vflat = Math.pow(G * M * a0_kpc, 0.25);
      const logV = Math.log10(Vflat);
      points.push({ logM, logV });
    }
    return points;
  }, [a0_kpc]);

  return (
    <Layout>
      <div className="space-y-8">
        <div>
          <h1 className="text-2xl md:text-3xl font-display font-bold text-white flex items-center gap-3">
            <div className="w-10 h-10 md:w-12 md:h-12 rounded-2xl bg-gradient-to-br from-amber-500 to-orange-600 flex items-center justify-center shadow-lg shadow-amber-500/20 flex-shrink-0">
              <Sigma className="w-5 h-5 md:w-7 md:h-7 text-white" />
            </div>
            The Equation
          </h1>
          <p className="text-slate-400 mt-2 text-sm md:text-base max-w-3xl">
            Building a testable equation from observation. Not a fit {"\u2014"} a derivation from the data itself.
          </p>
        </div>

        <div className="bg-gradient-to-r from-amber-500/10 to-orange-500/10 border border-amber-500/20 rounded-2xl p-4 sm:p-6 text-center">
          <p className="font-mono text-amber-300 text-sm sm:text-base italic leading-relaxed">
            {'"'}Goal: one equation that recovers Newton at high acceleration,{<br/>}
            gives RAR at low acceleration, and predicts Tully-Fisher {"\u2014"}{<br/>}
            with zero free parameters.{'"'}
          </p>
        </div>

        <StepCard step={1} title="Starting Point" icon={Telescope} color="cyan">
          <p className="text-sm text-slate-300 mb-4">
            We have two measurable accelerations at every point in every galaxy:
          </p>
          <div className="grid grid-cols-1 sm:grid-cols-2 gap-4 mb-4">
            <div className="bg-cyan-500/5 border border-cyan-500/20 rounded-xl p-4">
              <div className="font-mono text-cyan-400 text-lg mb-1">g<sub>bar</sub></div>
              <p className="text-xs text-slate-400">Acceleration predicted from visible (baryonic) matter alone: stars + gas, using Newtonian gravity.</p>
              <div className="font-mono text-xs text-slate-500 mt-2">g<sub>bar</sub> = GM<sub>bar</sub>(r) / r{"\u00B2"}</div>
            </div>
            <div className="bg-purple-500/5 border border-purple-500/20 rounded-xl p-4">
              <div className="font-mono text-purple-400 text-lg mb-1">g<sub>obs</sub></div>
              <p className="text-xs text-slate-400">Acceleration actually observed from rotation curves: V{"\u00B2"}(r)/r. Always exceeds g<sub>bar</sub> in the outer parts.</p>
              <div className="font-mono text-xs text-slate-500 mt-2">g<sub>obs</sub> = V{"\u00B2"}(r) / r</div>
            </div>
          </div>
          <p className="text-sm text-slate-300">
            We want a function <span className="font-mono text-white">g<sub>obs</sub> = f(g<sub>bar</sub>)</span> that captures the observed relationship across all galaxies with a single universal scale.
          </p>
        </StepCard>

        <StepCard step={2} title="The Acceleration Floor" icon={Sparkles} color="amber">
          <p className="text-sm text-slate-300 mb-4">
            The data shows a characteristic acceleration scale <span className="font-mono text-amber-400">a{"\u2080"}</span> below which a
            {" "}{"\""}compensation{"\"" } appears. The simplest interpretation: there is a {" "}
            <span className="text-amber-400 font-bold">minimum acceleration floor</span> that the universe enforces.
          </p>

          <MathBlock label="Core Equation (Acceleration Floor Form)">
            <div className="text-center">
              <span className="text-white">g<sub>obs</sub>{"\u00B2"}</span>
              <span className="text-slate-500"> = </span>
              <span className="text-cyan-400">g<sub>bar</sub>{"\u00B2"}</span>
              <span className="text-slate-500"> + </span>
              <span className="text-amber-400">a{"\u2080"}{"\u00B2"}</span>
              <span className="text-slate-500"> {"\u00D7"} </span>
              <span className="text-purple-400">{"\u03BD"}(g<sub>bar</sub>/a{"\u2080"})</span>
            </div>
          </MathBlock>

          <p className="text-xs text-slate-400 mb-4">
            where <span className="font-mono text-purple-400">{"\u03BD"}(x) = x/(1+x)</span> is a smooth interpolation that turns on the floor only when needed.
          </p>

          <MathBlock label="Expanded form">
            <div className="text-center">
              <span className="text-white">g<sub>obs</sub></span>
              <span className="text-slate-500"> = </span>
              <span className="text-slate-400">{"\u221A"}(</span>
              <span className="text-cyan-400">g<sub>bar</sub>{"\u00B2"}</span>
              <span className="text-slate-500"> + </span>
              <span className="text-amber-400">a{"\u2080"}{"\u00B2"}</span>
              <span className="text-slate-500"> {"\u00D7"} </span>
              <span className="text-cyan-400">g<sub>bar</sub></span>
              <span className="text-slate-500"> / (</span>
              <span className="text-cyan-400">g<sub>bar</sub></span>
              <span className="text-slate-500"> + </span>
              <span className="text-amber-400">a{"\u2080"}</span>
              <span className="text-slate-500">)</span>
              <span className="text-slate-400">)</span>
            </div>
          </MathBlock>

          <p className="text-xs text-slate-500 mt-3 italic">
            This is not a fit to data {"\u2014"} it is the simplest algebraic form that satisfies both required limits.
          </p>
        </StepCard>

        <StepCard step={3} title="Limit Verification" icon={CheckCircle2} color="emerald">
          <p className="text-sm text-slate-300 mb-4">
            Any valid equation must recover known physics in both limits. Let us verify:
          </p>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-6">
            <LimitCheck
              regime="High Acceleration (Newtonian)"
              condition={"g_bar \u226B a\u2080"}
              result={"g_obs\u00B2 \u2248 g_bar\u00B2 + a\u2080\u00B2 \u2248 g_bar\u00B2 \u2192 g_obs = g_bar \u2714"}
              icon={"\uD83C\uDF1F"}
            />
            <LimitCheck
              regime="Low Acceleration (Deep MOND)"
              condition={"g_bar \u226A a\u2080"}
              result={"\u03BD \u2192 g_bar/a\u2080, so g_obs\u00B2 \u2248 g_bar\u00B7a\u2080 \u2192 g_obs = \u221A(g_bar\u00B7a\u2080) \u2714"}
              icon={"\uD83C\uDF0C"}
            />
          </div>

          <div className="bg-slate-800/50 border border-white/5 rounded-xl p-4">
            <div className="text-xs text-slate-500 mb-2">Why the deep-MOND limit matters:</div>
            <p className="text-xs text-slate-300">
              When g<sub>obs</sub> = {"\u221A"}(g<sub>bar</sub>{"\u00B7"}a{"\u2080"}), the rotation velocity becomes:
            </p>
            <div className="font-mono text-sm text-amber-400 text-center my-2">
              V{"\u00B2"}/r = {"\u221A"}(GM/r{"\u00B2"} {"\u00D7"} a{"\u2080"})
            </div>
            <p className="text-xs text-slate-300">
              which gives <span className="font-mono text-white">V{"\u2074"} = G{"\u00B7"}M{"\u00B7"}a{"\u2080"}</span> {"\u2014"} a flat rotation curve
              that depends only on total mass. This is the <span className="text-amber-400 font-bold">Baryonic Tully-Fisher Relation</span>.
            </p>
          </div>
        </StepCard>

        <StepCard step={4} title="The Cosmological Connection" icon={Atom} color="purple">
          <p className="text-sm text-slate-300 mb-4">
            Now we make the leap. Instead of treating a{"\u2080"} as a fitted parameter, we {" "}
            <span className="text-purple-400 font-bold">derive it from cosmology</span>:
          </p>

          <MathBlock label="The Cosmological Prediction">
            <div className="text-center text-lg">
              <span className="text-amber-400">a{"\u2080"}</span>
              <span className="text-slate-500"> = </span>
              <span className="text-purple-400">c {"\u00D7"} H{"\u2080"}</span>
              <span className="text-slate-500"> / </span>
              <span className="text-cyan-400">2{"\u03C0"}</span>
            </div>
          </MathBlock>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 mb-4">
            <div className="bg-slate-800/50 rounded-lg p-3 text-center">
              <div className="text-xs text-slate-500">Predicted</div>
              <div className="font-mono text-purple-400">{(cH0 / (2 * Math.PI)).toExponential(2)} m/s{"\u00B2"}</div>
            </div>
            <div className="bg-slate-800/50 rounded-lg p-3 text-center">
              <div className="text-xs text-slate-500">Observed</div>
              <div className="font-mono text-cyan-400">{a0_ms2.toExponential(1)} m/s{"\u00B2"}</div>
            </div>
            <div className="bg-slate-800/50 rounded-lg p-3 text-center">
              <div className="text-xs text-slate-500">Mismatch</div>
              <div className="font-mono text-amber-400">~{(Math.abs((cH0 / (2 * Math.PI) - a0_ms2) / a0_ms2) * 100).toFixed(0)}%</div>
            </div>
          </div>

          <MathBlock label="The Complete Equation (Zero Free Parameters)">
            <div className="text-center text-lg">
              <span className="text-white">g<sub>obs</sub></span>
              <span className="text-slate-500"> = </span>
              <span className="text-slate-400">{"\u221A"}(</span>
              <span className="text-cyan-400">g<sub>bar</sub>{"\u00B2"}</span>
              <span className="text-slate-500"> + </span>
              <span className="text-slate-400">(</span>
              <span className="text-purple-400">cH{"\u2080"}/2{"\u03C0"}</span>
              <span className="text-slate-400">){"\u00B2"}</span>
              <span className="text-slate-500"> {"\u00D7"} </span>
              <span className="text-cyan-400">g<sub>bar</sub></span>
              <span className="text-slate-500"> / (</span>
              <span className="text-cyan-400">g<sub>bar</sub></span>
              <span className="text-slate-500"> + </span>
              <span className="text-purple-400">cH{"\u2080"}/2{"\u03C0"}</span>
              <span className="text-slate-500">)</span>
              <span className="text-slate-400">)</span>
            </div>
          </MathBlock>

          <p className="text-xs text-slate-400 mt-3">
            This equation has <span className="text-white font-bold">zero adjustable parameters</span>. Every quantity is either
            measured locally (g<sub>bar</sub>) or fixed by cosmology (c, H{"\u2080"}).
          </p>
        </StepCard>

        <StepCard step={5} title="Confronting the Data" icon={TrendingUp} color="cyan">
          <p className="text-sm text-slate-300 mb-4">
            The equation plotted against {rarPoints.length > 0 ? rarPoints.length.toLocaleString() : "4,123"} observed data points
            from SPARC + LITTLE THINGS galaxies:
          </p>

          <div className="h-[350px] sm:h-[400px] mb-4">
            <ResponsiveContainer width="100%" height="100%">
              <ComposedChart margin={{ top: 10, right: 20, bottom: 40, left: 50 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                <XAxis
                  dataKey="logGbar"
                  type="number"
                  domain={[-0.5, 5.5]}
                  label={{ value: "log g_bar [(km/s)\u00B2/kpc]", position: "bottom", offset: 20, style: { fill: '#94a3b8', fontSize: 11 } }}
                  tick={{ fill: '#64748b', fontSize: 10 }}
                />
                <YAxis
                  type="number"
                  domain={[-0.5, 5.5]}
                  label={{ value: "log g_obs [(km/s)\u00B2/kpc]", angle: -90, position: "insideLeft", offset: -35, style: { fill: '#94a3b8', fontSize: 11 } }}
                  tick={{ fill: '#64748b', fontSize: 10 }}
                />
                <Tooltip
                  contentStyle={{ backgroundColor: '#0f172a', border: '1px solid rgba(255,255,255,0.1)', borderRadius: 12, fontSize: 11 }}
                  itemStyle={{ color: '#e2e8f0' }}
                />
                {rarPoints.length > 0 && (
                  <Scatter data={rarPoints} dataKey="logGobs" fill="rgba(6,182,212,0.15)" r={1.5} name="Observed" />
                )}
                <Line data={modelCurve} dataKey="logGobs" stroke="#f59e0b" strokeWidth={2.5} dot={false} name="Our Model" type="monotone" />
                <Line data={modelCurve} dataKey="logNewton" stroke="#64748b" strokeWidth={1} strokeDasharray="6 4" dot={false} name="Newton (1:1)" type="monotone" />
                <Line data={modelCurve} dataKey="logMOND" stroke="#a78bfa" strokeWidth={1} strokeDasharray="3 3" dot={false} name="Simple MOND" type="monotone" />
              </ComposedChart>
            </ResponsiveContainer>
          </div>

          <div className="flex flex-wrap gap-4 text-xs justify-center">
            <span className="flex items-center gap-1"><span className="w-3 h-3 rounded-full bg-cyan-500/30 inline-block" /> Data</span>
            <span className="flex items-center gap-1"><span className="w-6 h-0.5 bg-amber-500 inline-block" /> Our Model</span>
            <span className="flex items-center gap-1"><span className="w-6 h-0.5 bg-slate-500 inline-block border-dashed" /> Newton</span>
            <span className="flex items-center gap-1"><span className="w-6 h-0.5 bg-violet-400 inline-block" /> MOND</span>
          </div>
        </StepCard>

        <StepCard step={6} title="Tully-Fisher Derivation" icon={TrendingUp} color="emerald">
          <p className="text-sm text-slate-300 mb-4">
            In the deep low-acceleration regime (outer galaxy), our equation predicts:
          </p>

          <div className="space-y-3 mb-4">
            <MathBlock label="Step 1: Deep limit">
              <span className="text-white">g<sub>obs</sub></span>
              <span className="text-slate-500"> {"\u2248"} </span>
              <span className="text-amber-400">{"\u221A"}(g<sub>bar</sub> {"\u00D7"} a{"\u2080"})</span>
            </MathBlock>

            <MathBlock label="Step 2: Substitute V\u00B2/r and GM/r\u00B2">
              <span className="text-white">V{"\u00B2"}/r</span>
              <span className="text-slate-500"> = </span>
              <span className="text-amber-400">{"\u221A"}(GM/r{"\u00B2"} {"\u00D7"} a{"\u2080"})</span>
            </MathBlock>

            <MathBlock label="Step 3: Solve for V">
              <div className="text-center text-lg">
                <span className="text-emerald-400 font-bold">V<sub>flat</sub>{"\u2074"}</span>
                <span className="text-slate-500"> = </span>
                <span className="text-white">G {"\u00D7"} M<sub>bar</sub> {"\u00D7"} a{"\u2080"}</span>
              </div>
            </MathBlock>
          </div>

          <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-4">
            <div className="flex items-center gap-2 mb-2">
              <CheckCircle2 className="w-4 h-4 text-emerald-400" />
              <span className="text-sm font-bold text-emerald-400">Baryonic Tully-Fisher Relation {"\u2014"} Derived, Not Assumed</span>
            </div>
            <p className="text-xs text-slate-300 leading-relaxed">
              The flat rotation velocity scales as M<sup>1/4</sup> {"\u2014"} exactly the observed
              Baryonic Tully-Fisher Relation. This is not inserted by hand; it {" "}
              <span className="text-emerald-400 font-bold">emerges automatically</span> from the equation in the low-acceleration limit.
              This is a strong consistency check.
            </p>
          </div>

          <div className="h-[250px] mt-4">
            <ResponsiveContainer width="100%" height="100%">
              <LineChart data={tfPoints} margin={{ top: 10, right: 20, bottom: 35, left: 50 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                <XAxis
                  dataKey="logM"
                  type="number"
                  domain={[7, 12]}
                  label={{ value: "log M_bar [M\u2609]", position: "bottom", offset: 15, style: { fill: '#94a3b8', fontSize: 11 } }}
                  tick={{ fill: '#64748b', fontSize: 10 }}
                />
                <YAxis
                  dataKey="logV"
                  type="number"
                  domain={[1, 3]}
                  label={{ value: "log V_flat [km/s]", angle: -90, position: "insideLeft", offset: -35, style: { fill: '#94a3b8', fontSize: 11 } }}
                  tick={{ fill: '#64748b', fontSize: 10 }}
                />
                <Line dataKey="logV" stroke="#10b981" strokeWidth={2} dot={false} name="Predicted V_flat" />
              </LineChart>
            </ResponsiveContainer>
          </div>
          <p className="text-xs text-slate-500 text-center mt-1">
            Predicted Tully-Fisher relation: slope = 1/4 in log-log space (V {"\u221D"} M<sup>0.25</sup>)
          </p>
        </StepCard>

        <StepCard step={7} title="Three Physical Readings" icon={FlaskConical} color="purple">
          <p className="text-sm text-slate-300 mb-4">
            The same equation admits three different physical interpretations.
            Distinguishing between them requires new data {"\u2014"} the equation alone cannot decide.
          </p>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div className="bg-cyan-500/5 border border-cyan-500/20 rounded-xl p-4">
              <div className="text-sm font-bold text-cyan-400 mb-2">Reading A: Modified Gravity</div>
              <p className="text-xs text-slate-300 mb-3">
                The equation IS the fundamental law. Gravity itself changes below a{"\u2080"}.
                There is no dark matter.
              </p>
              <div className="text-xs text-slate-500">
                <span className="text-emerald-400">{"\u2713"}</span> Simple, one law for all galaxies{<br/>}
                <span className="text-rose-400">{"\u2717"}</span> Struggles with galaxy clusters{<br/>}
                <span className="text-rose-400">{"\u2717"}</span> No relativistic extension proven
              </div>
            </div>

            <div className="bg-purple-500/5 border border-purple-500/20 rounded-xl p-4">
              <div className="text-sm font-bold text-purple-400 mb-2">Reading B: Effective Dark Matter</div>
              <p className="text-xs text-slate-300 mb-3">
                Dark matter exists, but its local density is locked to baryons:
                g<sub>DM</sub> = g<sub>obs</sub> - g<sub>bar</sub>, fully determined by g<sub>bar</sub>.
              </p>
              <div className="text-xs text-slate-500">
                <span className="text-emerald-400">{"\u2713"}</span> Consistent with CMB + lensing{<br/>}
                <span className="text-rose-400">{"\u2717"}</span> No mechanism for tight coupling{<br/>}
                <span className="text-rose-400">{"\u2717"}</span> Why does it know about a{"\u2080"}?
              </div>
            </div>

            <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-4">
              <div className="text-sm font-bold text-amber-400 mb-2">Reading C: Emergent Cosmic Effect</div>
              <p className="text-xs text-slate-300 mb-3">
                a{"\u2080"} = cH{"\u2080"}/2{"\u03C0"} is not coincidence {"\u2014"} the cosmic expansion creates a
                minimum acceleration floor that couples to local gravity.
              </p>
              <div className="text-xs text-slate-500">
                <span className="text-emerald-400">{"\u2713"}</span> Explains the cH{"\u2080"} coincidence{<br/>}
                <span className="text-emerald-400">{"\u2713"}</span> Predicts a{"\u2080"} evolves with H(z){<br/>}
                <span className="text-rose-400">{"\u2717"}</span> No mechanism yet
              </div>
            </div>
          </div>
        </StepCard>

        <StepCard step={8} title="Unique Predictions" icon={Telescope} color="rose">
          <p className="text-sm text-slate-300 mb-4">
            To elevate this equation from {'"'}observation summary{'"'} to {'"'}testable theory,{'"'} it must make {" "}
            <span className="text-rose-400 font-bold">predictions that can be falsified</span>:
          </p>

          <div className="space-y-4">
            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <div className="flex items-center gap-2 mb-2">
                <div className="w-6 h-6 rounded-full bg-rose-500/20 flex items-center justify-center text-rose-400 text-xs font-bold">1</div>
                <span className="text-sm font-bold text-rose-400">Redshift Evolution</span>
              </div>
              <p className="text-xs text-slate-300">
                If a{"\u2080"} = cH{"\u2080"}/2{"\u03C0"}, then at redshift z the scale becomes a{"\u2080"}(z) = cH(z)/2{"\u03C0"}.
                Since H(z) {">"} H{"\u2080"} at high z, rotation curves of distant galaxies should show a {" "}
                <span className="text-rose-400 font-bold">larger</span> acceleration floor.
              </p>
              <MathBlock>
                <span className="text-rose-400">a{"\u2080"}(z)</span>
                <span className="text-slate-500"> = </span>
                <span className="text-white">c {"\u00D7"} H(z) / 2{"\u03C0"}</span>
                <span className="text-slate-500"> = </span>
                <span className="text-amber-400">a{"\u2080"} {"\u00D7"} E(z)</span>
              </MathBlock>
              <p className="text-xs text-slate-500 italic">
                where E(z) = H(z)/H{"\u2080"}. At z=1, E {"\u2248"} 1.7, so a{"\u2080"} should be ~70% larger.
                JWST and ALMA can test this with high-z rotation curves.
              </p>
            </div>

            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <div className="flex items-center gap-2 mb-2">
                <div className="w-6 h-6 rounded-full bg-rose-500/20 flex items-center justify-center text-rose-400 text-xs font-bold">2</div>
                <span className="text-sm font-bold text-rose-400">Galaxy Clusters</span>
              </div>
              <p className="text-xs text-slate-300">
                Apply the same equation to galaxy clusters. If the floor model is universal,
                it must reduce the missing mass problem in clusters {"\u2014"} but likely not eliminate it entirely.
                The <span className="text-rose-400 font-bold">residual deficit</span> in clusters would quantify
                how much {'"'}real{'"'} dark matter is needed beyond the floor effect.
              </p>
            </div>

            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <div className="flex items-center gap-2 mb-2">
                <div className="w-6 h-6 rounded-full bg-rose-500/20 flex items-center justify-center text-rose-400 text-xs font-bold">3</div>
                <span className="text-sm font-bold text-rose-400">Wide Binaries</span>
              </div>
              <p className="text-xs text-slate-300">
                At separations {">"} 7000 AU, binary star orbits enter the regime g {"<"} a{"\u2080"}.
                Our equation predicts their orbital velocities should deviate from Kepler in a
                {" "}specific, calculable way. Gaia DR4 data can test this directly.
              </p>
            </div>

            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <div className="flex items-center gap-2 mb-2">
                <div className="w-6 h-6 rounded-full bg-rose-500/20 flex items-center justify-center text-rose-400 text-xs font-bold">4</div>
                <span className="text-sm font-bold text-rose-400">Scatter vs. Feedback</span>
              </div>
              <p className="text-xs text-slate-300">
                If the scatter around the RAR is random (not correlated with star formation rate,
                morphology, or environment), this supports a fundamental law.
                If scatter correlates with feedback strength, it supports {"\u039B"}CDM + baryonic processes.
                Our data shows <span className="text-amber-400 font-bold">no correlation</span> {"\u2014"} but more samples are needed.
              </p>
            </div>
          </div>
        </StepCard>

        <GlassCard glow="amber">
          <div className="text-center space-y-4 py-4">
            <div className="text-xs text-amber-400 font-bold tracking-widest">THE EQUATION</div>
            <div className="font-mono text-lg sm:text-xl text-white leading-relaxed">
              g<sub>obs</sub> = {"\u221A"}(g<sub>bar</sub>{"\u00B2"} + (cH{"\u2080"}/2{"\u03C0"}){"\u00B2"} {"\u00D7"} g<sub>bar</sub>/(g<sub>bar</sub> + cH{"\u2080"}/2{"\u03C0"}))
            </div>
            <div className="h-px bg-gradient-to-r from-transparent via-amber-500/30 to-transparent" />
            <div className="grid grid-cols-1 sm:grid-cols-4 gap-3 text-xs">
              <div className="bg-slate-800/50 rounded-lg p-3">
                <div className="text-slate-500 mb-1">Free Parameters</div>
                <div className="font-mono text-2xl text-amber-400 font-bold">0</div>
              </div>
              <div className="bg-slate-800/50 rounded-lg p-3">
                <div className="text-slate-500 mb-1">Recovers Newton</div>
                <div className="font-mono text-emerald-400">{"\u2714"} at high g</div>
              </div>
              <div className="bg-slate-800/50 rounded-lg p-3">
                <div className="text-slate-500 mb-1">Gives RAR</div>
                <div className="font-mono text-emerald-400">{"\u2714"} at low g</div>
              </div>
              <div className="bg-slate-800/50 rounded-lg p-3">
                <div className="text-slate-500 mb-1">Predicts Tully-Fisher</div>
                <div className="font-mono text-emerald-400">{"\u2714"} V{"\u2074"} = GMa{"\u2080"}</div>
              </div>
            </div>

            <div className="text-xs text-slate-500 italic max-w-xl mx-auto">
              Proposed name: <span className="text-amber-400 font-bold">Cosmic Acceleration Floor Model</span>{<br/>}
              Status: Empirically motivated, not yet a complete theory. Missing: physical mechanism,
              relativistic formulation, cluster-scale validation.
            </div>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
