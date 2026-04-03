import React, { useState, useMemo, useEffect } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { Telescope, CheckCircle2, XCircle, Sparkles, Info, Clock, AlertTriangle, BarChart3 } from 'lucide-react';
import {
  LineChart, Line, XAxis, YAxis, CartesianGrid,
  Tooltip, ResponsiveContainer, ReferenceLine, Legend,
  ErrorBar, ComposedChart, Scatter, Area,
  BarChart, Bar, Cell
} from 'recharts';

const G_KPC = 4.3009e-6;
const A0_LOCAL = 3702;
const A0_MS2 = 1.2e-10;
const H0_KMS_MPC = 67.4;

function Ez(z: number, omM: number, omL: number): number {
  return Math.sqrt(omM * Math.pow(1 + z, 3) + omL);
}

function Hz(z: number, omM: number, omL: number): number {
  return H0_KMS_MPC * Ez(z, omM, omL);
}

function a0AtZ(z: number, omM: number, omL: number): number {
  return A0_LOCAL * Ez(z, omM, omL);
}

function a0AtZ_ms2(z: number, omM: number, omL: number): number {
  return A0_MS2 * Ez(z, omM, omL);
}

function cosmicFloorVelocity(r: number, M: number, a0: number): number {
  if (r <= 0) return 0;
  const gbar = G_KPC * M / (r * r);
  const gobs2 = gbar * gbar + a0 * a0 * gbar / (gbar + a0);
  return Math.sqrt(Math.sqrt(gobs2) * r);
}

function newtonianVelocity(r: number, M: number): number {
  if (r <= 0) return 0;
  return Math.sqrt(G_KPC * M / r);
}

function btfrVflat(M: number, a0: number): number {
  return Math.pow(G_KPC * M * a0, 0.25);
}

function generateMockJWSTData(z: number, omM: number, omL: number) {
  const zValues = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0];
  const rng = (seed: number) => {
    let s = seed;
    return () => { s = (s * 16807 + 0) % 2147483647; return s / 2147483647; };
  };

  return zValues.filter(zv => zv <= z + 0.1).map(zv => {
    const rand = rng(Math.round(zv * 1000));
    const e = Ez(zv, omM, omL);
    const trueA0 = A0_MS2 * e;
    const noise = (rand() - 0.5) * 0.15;
    const measured = trueA0 * (1 + noise);
    const errSize = trueA0 * (0.11 + 0.045 * zv);
    return {
      z: zv,
      measured: measured * 1e10,
      errorBar: errSize * 1e10,
      cosmicFloor: (A0_MS2 * e) * 1e10,
      mond: A0_MS2 * 1e10,
    };
  });
}

function MathBlock({ children }: { children: React.ReactNode }) {
  return (
    <div className="my-3 bg-slate-900/80 border border-white/10 rounded-xl px-4 sm:px-6 py-3 font-mono text-sm sm:text-base text-cyan-300 overflow-x-auto">
      {children}
    </div>
  );
}

interface FeasibilityData {
  signalTable: Array<{ z: number; Ez: number; delta_a0_pct: number; delta_Vflat_pct: number }>;
  existingSurveys: Array<{ name: string; zRange: string; nGalaxies: string | number; method: string; resolution: string; vPrecision: string; a0Constraint: string; usability: string }>;
  powerAnalysis: { signal_deltaV_frac: number; results: Array<{ sigmaV_frac: number; 'N_3\u03C3': number; 'N_5\u03C3': number }> };
  timeline: { now: { verdict: string }; near: { verdict: string; when: string }; future: { verdict: string; when: string } };
  honestAssessment: { canTestNow: boolean; canTestSoon: string; canTestDefinitively: string; mainObstacle: string; alternativeTests: string[] };
}

export default function RedshiftLabPage() {
  const [z, setZ] = useState(1.0);
  const [omM, setOmM] = useState(0.315);
  const [omL, setOmL] = useState(0.685);
  const [flatUniverse, setFlatUniverse] = useState(true);
  const [refMass] = useState(1e11);
  const [feasibility, setFeasibility] = useState<FeasibilityData | null>(null);

  useEffect(() => {
    fetch(import.meta.env.BASE_URL + 'redshift-feasibility.json')
      .then(r => r.json())
      .then(d => setFeasibility(d))
      .catch(() => {});
  }, []);

  const e = useMemo(() => Ez(z, omM, omL), [z, omM, omL]);
  const hZ = useMemo(() => Hz(z, omM, omL), [z, omM, omL]);
  const a0z = useMemo(() => a0AtZ(z, omM, omL), [z, omM, omL]);
  const a0z_ms2 = useMemo(() => a0AtZ_ms2(z, omM, omL), [z, omM, omL]);

  const evolutionData = useMemo(() => {
    const points = [];
    const seed = 42;
    for (let zi = 0; zi <= 3.0; zi += 0.05) {
      const ei = Ez(zi, omM, omL);
      const idx = Math.round(zi * 20);
      const scatter = Math.sin(idx * 7.3 + seed) * 0.1 + Math.cos(idx * 3.1) * 0.05;
      points.push({
        z: parseFloat(zi.toFixed(2)),
        cosmicFloor: parseFloat((A0_MS2 * ei * 1e10).toFixed(4)),
        mond: parseFloat((A0_MS2 * 1e10).toFixed(4)),
        lcdmUpper: parseFloat((A0_MS2 * (1.15 + scatter) * 1e10).toFixed(4)),
        lcdmLower: parseFloat((A0_MS2 * (0.85 + scatter) * 1e10).toFixed(4)),
        lcdmMid: parseFloat((A0_MS2 * (1.0 + scatter) * 1e10).toFixed(4)),
      });
    }
    return points;
  }, [omM, omL]);

  const mockData = useMemo(() => generateMockJWSTData(z, omM, omL), [z, omM, omL]);

  const rotationCurveZ0 = useMemo(() => {
    const points = [];
    for (let r = 0.5; r <= 50; r += 0.5) {
      points.push({
        r: parseFloat(r.toFixed(1)),
        newton: parseFloat(newtonianVelocity(r, refMass).toFixed(2)),
        floor_z0: parseFloat(cosmicFloorVelocity(r, refMass, A0_LOCAL).toFixed(2)),
        floor_zN: parseFloat(cosmicFloorVelocity(r, refMass, a0z).toFixed(2)),
      });
    }
    return points;
  }, [refMass, a0z]);

  const btfrData = useMemo(() => {
    const masses = [];
    for (let logM = 7; logM <= 12; logM += 0.2) {
      const M = Math.pow(10, logM);
      masses.push({
        logM: parseFloat(logM.toFixed(1)),
        vflat_z0: parseFloat(btfrVflat(M, A0_LOCAL).toFixed(2)),
        vflat_zN: parseFloat(btfrVflat(M, a0z).toFixed(2)),
      });
    }
    return masses;
  }, [a0z]);

  const pctIncrease = ((e - 1) * 100).toFixed(0);

  return (
    <Layout>
      <div className="space-y-8">
        <div>
          <div className="flex items-center gap-3 mb-4">
            <Telescope className="w-8 h-8 text-violet-400" />
            <div>
              <h1 className="text-3xl font-bold font-display text-white">Redshift Evolution Lab</h1>
              <p className="text-slate-400 text-sm">The killer test: does a{"\u2080"} evolve with cosmic expansion?</p>
            </div>
          </div>
          <div className="bg-gradient-to-r from-violet-500/10 to-cyan-500/10 border border-violet-500/20 rounded-xl p-4 sm:p-6">
            <p className="text-violet-300 font-mono text-sm leading-relaxed text-center italic">
              "If a{"\u2080"} = cH{"\u2080"}/2{"\u03C0"}, then a{"\u2080"}(z) = cH(z)/2{"\u03C0"} {"\u2014"} and distant galaxies should show a larger acceleration floor."
            </p>
          </div>
        </div>

        <GlassCard glow="purple">
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-violet-500 to-purple-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">1</div>
            <div>
              <h3 className="text-lg font-bold text-white">The Core Prediction</h3>
              <p className="text-xs text-slate-400">What distinguishes the Cosmic Floor from MOND and {"\u039B"}CDM</p>
            </div>
          </div>

          <p className="text-sm text-slate-300 mb-4">
            If a{"\u2080"} is truly set by cosmology (a{"\u2080"} = cH/2{"\u03C0"}), then it must change
            with the Hubble parameter H(z). Since H was larger in the past, galaxies at
            high redshift should have a <span className="text-violet-400 font-bold">larger</span> acceleration floor:
          </p>

          <MathBlock>
            <span className="text-violet-400">a{"\u2080"}(z)</span>
            <span className="text-slate-500"> = </span>
            <span className="text-white">cH(z)/2{"\u03C0"}</span>
            <span className="text-slate-500"> = </span>
            <span className="text-amber-400">a{"\u2080"} {"\u00D7"} E(z)</span>
          </MathBlock>

          <MathBlock>
            <span className="text-slate-400">E(z) = H(z)/H{"\u2080"} = </span>
            <span className="text-cyan-400">{"\u221A"}({"\u03A9"}{"\u2098"}(1+z){"\u00B3"} + {"\u03A9"}{"\u039B"})</span>
          </MathBlock>

          <div className="grid grid-cols-1 sm:grid-cols-3 gap-3 mt-4">
            <div className="bg-violet-500/10 border border-violet-500/20 rounded-xl p-3 text-center">
              <div className="text-xs text-slate-500 mb-1">MOND</div>
              <div className="text-sm font-bold text-violet-400">a{"\u2080"} = const</div>
              <div className="text-xs text-slate-400 mt-1">Never changes</div>
            </div>
            <div className="bg-slate-500/10 border border-slate-500/20 rounded-xl p-3 text-center">
              <div className="text-xs text-slate-500 mb-1">{"\u039B"}CDM / Feedback</div>
              <div className="text-sm font-bold text-slate-400">No prediction</div>
              <div className="text-xs text-slate-400 mt-1">a{"\u2080"} is emergent</div>
            </div>
            <div className="bg-amber-500/10 border border-amber-500/20 rounded-xl p-3 text-center">
              <div className="text-xs text-slate-500 mb-1">Cosmic Floor</div>
              <div className="text-sm font-bold text-amber-400">a{"\u2080"} {"\u00D7"} E(z)</div>
              <div className="text-xs text-slate-400 mt-1">Scales with H(z)</div>
            </div>
          </div>
        </GlassCard>

        <GlassCard glow="cyan">
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">2</div>
            <div>
              <h3 className="text-lg font-bold text-white">Interactive Controls</h3>
              <p className="text-xs text-slate-400">Adjust redshift and cosmological parameters</p>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-6">
            <div>
              <div className="flex justify-between items-end mb-2">
                <label className="font-mono text-violet-400 font-medium text-sm">Redshift z</label>
                <span className="font-mono text-xs text-white">{z.toFixed(2)}</span>
              </div>
              <input
                type="range" min={0} max={3} step={0.05}
                value={z}
                onChange={(e) => setZ(Number(e.target.value))}
                className="w-full accent-violet-500"
              />
              <div className="flex justify-between text-xs text-slate-500 mt-1">
                <span>z=0 (now)</span>
                <span>z=3 (~11.5 Gyr ago)</span>
              </div>
            </div>

            <div>
              <div className="flex justify-between items-end mb-2">
                <label className="font-mono text-cyan-400 font-medium text-sm">{"\u03A9"}{"\u2098"}</label>
                <span className="font-mono text-xs text-white">{omM.toFixed(3)}</span>
              </div>
              <input
                type="range" min={0.1} max={0.5} step={0.005}
                value={omM}
                onChange={(e) => { const v = Number(e.target.value); setOmM(v); if (flatUniverse) setOmL(1 - v); }}
                className="w-full accent-cyan-500"
              />
              <div className="flex justify-between text-xs text-slate-500 mt-1">
                <span>0.1</span>
                <span>Planck: 0.315</span>
                <span>0.5</span>
              </div>
            </div>

            <div>
              <div className="flex justify-between items-end mb-2">
                <label className="font-mono text-emerald-400 font-medium text-sm">{"\u03A9"}{"\u039B"}</label>
                <span className="font-mono text-xs text-white">{omL.toFixed(3)}</span>
              </div>
              <input
                type="range" min={0.0} max={1.0} step={0.005}
                value={omL}
                onChange={(e) => { const v = Number(e.target.value); setOmL(v); if (flatUniverse) setOmM(1 - v); }}
                className="w-full accent-emerald-500"
                disabled={flatUniverse}
              />
              <div className="flex items-center gap-2 mt-1">
                <label className="flex items-center gap-1.5 text-xs text-slate-400 cursor-pointer">
                  <input
                    type="checkbox"
                    checked={flatUniverse}
                    onChange={(e) => { setFlatUniverse(e.target.checked); if (e.target.checked) setOmL(1 - omM); }}
                    className="accent-emerald-500"
                  />
                  Flat universe ({"\u03A9"}{"\u2098"}+{"\u03A9"}{"\u039B"}=1)
                </label>
              </div>
            </div>

            <div className="space-y-3">
              <div className="bg-slate-800/50 border border-white/10 rounded-xl p-3">
                <div className="grid grid-cols-2 gap-2 text-xs">
                  <div>
                    <span className="text-slate-500">E(z):</span>
                    <span className="text-white font-mono ml-1">{e.toFixed(3)}</span>
                  </div>
                  <div>
                    <span className="text-slate-500">H(z):</span>
                    <span className="text-white font-mono ml-1">{hZ.toFixed(1)} km/s/Mpc</span>
                  </div>
                  <div>
                    <span className="text-slate-500">a{"\u2080"}(z):</span>
                    <span className="text-amber-400 font-mono ml-1">{a0z.toFixed(0)} (km/s){"\u00B2"}/kpc</span>
                  </div>
                  <div>
                    <span className="text-slate-500">Increase:</span>
                    <span className="text-emerald-400 font-mono ml-1">+{pctIncrease}%</span>
                  </div>
                </div>
              </div>
              {!flatUniverse && (omM + omL) !== 1 && (
                <div className="text-xs text-amber-400">
                  {"\u26A0"} Non-flat: {"\u03A9"}{"\u2098"}+{"\u03A9"}{"\u039B"} = {(omM + omL).toFixed(3)} {(omM + omL) > 1 ? "(closed)" : "(open)"}
                </div>
              )}
            </div>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-amber-500 to-orange-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">3</div>
            <div>
              <h3 className="text-lg font-bold text-white">a{"\u2080"} Evolution: Three Models vs. Mock JWST Data</h3>
              <p className="text-xs text-slate-400">How the acceleration scale changes with redshift under each framework</p>
            </div>
          </div>

          <div className="h-[350px] sm:h-[400px]">
            <ResponsiveContainer width="100%" height="100%">
              <ComposedChart data={evolutionData} margin={{ top: 10, right: 30, left: 10, bottom: 10 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                <XAxis
                  dataKey="z" type="number" domain={[0, 3]}
                  stroke="rgba(255,255,255,0.3)"
                  tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                  label={{ value: 'Redshift z', position: 'insideBottom', offset: -5, fill: 'rgba(255,255,255,0.5)', fontSize: 12 }}
                />
                <YAxis
                  stroke="rgba(255,255,255,0.3)"
                  tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                  label={{ value: 'a\u2080 (\u00D710\u207B\u00B9\u2070 m/s\u00B2)', angle: -90, position: 'insideLeft', fill: 'rgba(255,255,255,0.5)', fontSize: 12, dx: -5 }}
                  domain={[0, 'auto']}
                />
                <Tooltip
                  contentStyle={{ backgroundColor: 'rgba(15,23,42,0.95)', border: '1px solid rgba(255,255,255,0.1)', borderRadius: '12px', fontSize: '12px' }}
                  formatter={(value: number, name: string) => {
                    const labels: Record<string, string> = { cosmicFloor: 'Cosmic Floor', mond: 'MOND (constant)', lcdmMid: '\u039BCDM (emergent, scattered)' };
                    if (name === 'lcdmUpper' || name === 'lcdmLower') return [null, null];
                    return [value.toFixed(3) + ' \u00D710\u207B\u00B9\u2070 m/s\u00B2', labels[name] || name];
                  }}
                  labelFormatter={(v) => 'z = ' + Number(v).toFixed(2)}
                />
                <Area type="monotone" dataKey="lcdmUpper" stroke="none" fill="rgba(239,68,68,0.1)" legendType="none" />
                <Area type="monotone" dataKey="lcdmLower" stroke="none" fill="rgba(10,14,26,1)" legendType="none" />
                <Line type="monotone" dataKey="lcdmMid" stroke="#ef4444" strokeWidth={1.5} strokeDasharray="4 4" dot={false} name={"\u039BCDM (no prediction)"} />
                <Line type="monotone" dataKey="cosmicFloor" stroke="#f59e0b" strokeWidth={3} dot={false} name="Cosmic Floor" />
                <Line type="monotone" dataKey="mond" stroke="#8b5cf6" strokeWidth={2} strokeDasharray="8 4" dot={false} name="MOND (constant)" />
                <ReferenceLine x={z} stroke="rgba(255,255,255,0.3)" strokeDasharray="4 4" label={{ value: 'z=' + z.toFixed(1), fill: 'rgba(255,255,255,0.5)', fontSize: 11 }} />
                <Scatter data={mockData} dataKey="measured" fill="#22d3ee" name="Mock JWST">
                  <ErrorBar dataKey="errorBar" width={4} strokeWidth={1.5} stroke="#22d3ee" direction="y" />
                </Scatter>
                <Legend
                  wrapperStyle={{ fontSize: '11px', color: 'rgba(255,255,255,0.6)' }}
                />
              </ComposedChart>
            </ResponsiveContainer>
          </div>

          <div className="mt-4 bg-cyan-500/5 border border-cyan-500/20 rounded-xl p-4">
            <div className="flex items-center gap-2 mb-2">
              <Info className="w-4 h-4 text-cyan-400" />
              <span className="text-xs font-bold text-cyan-400">Mock Data Points</span>
            </div>
            <p className="text-xs text-slate-400">
              Cyan dots show simulated JWST/ALMA measurements assuming the Cosmic Floor is correct.
              Error bars grow with redshift (harder to measure distant galaxies).
              Real data from JWST rotation curve programs will replace these.
            </p>
          </div>
        </GlassCard>

        <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
          <GlassCard>
            <div className="flex items-center gap-3 mb-4">
              <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-emerald-500 to-teal-600 flex items-center justify-center text-white font-bold text-xs flex-shrink-0">4a</div>
              <div>
                <h3 className="text-base font-bold text-white">Rotation Curves: z=0 vs z={z.toFixed(1)}</h3>
                <p className="text-xs text-slate-400">Same galaxy (M=10{"\u00B9\u00B9"} M{"\u2609"}), different epochs</p>
              </div>
            </div>

            <div className="h-[280px]">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={rotationCurveZ0} margin={{ top: 5, right: 15, left: 5, bottom: 5 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                  <XAxis
                    dataKey="r" type="number"
                    stroke="rgba(255,255,255,0.3)"
                    tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 10 }}
                    label={{ value: 'r (kpc)', position: 'insideBottom', offset: -3, fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                  />
                  <YAxis
                    stroke="rgba(255,255,255,0.3)"
                    tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 10 }}
                    label={{ value: 'V (km/s)', angle: -90, position: 'insideLeft', fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                  />
                  <Tooltip
                    contentStyle={{ backgroundColor: 'rgba(15,23,42,0.95)', border: '1px solid rgba(255,255,255,0.1)', borderRadius: '12px', fontSize: '11px' }}
                    formatter={(v: number, name: string) => {
                      const labels: Record<string, string> = { newton: 'Newtonian', floor_z0: 'Floor (z=0)', floor_zN: 'Floor (z=' + z.toFixed(1) + ')' };
                      return [v.toFixed(1) + ' km/s', labels[name] || name];
                    }}
                  />
                  <Line type="monotone" dataKey="newton" stroke="#f97316" strokeWidth={1.5} strokeDasharray="5 3" dot={false} name="Newton" />
                  <Line type="monotone" dataKey="floor_z0" stroke="#22d3ee" strokeWidth={2} dot={false} name="Floor z=0" />
                  <Line type="monotone" dataKey="floor_zN" stroke="#f59e0b" strokeWidth={2.5} dot={false} name={"Floor z=" + z.toFixed(1)} />
                  <Legend wrapperStyle={{ fontSize: '10px' }} />
                </LineChart>
              </ResponsiveContainer>
            </div>

            <p className="text-xs text-slate-400 mt-3">
              At z={z.toFixed(1)}, the acceleration floor is {pctIncrease}% higher {"\u2192"} the flat
              part of the rotation curve sits <span className="text-amber-400 font-bold">higher</span>.
              Newtonian curve is identical (gravity doesn't change with expansion).
            </p>
          </GlassCard>

          <GlassCard>
            <div className="flex items-center gap-3 mb-4">
              <div className="w-8 h-8 rounded-lg bg-gradient-to-br from-rose-500 to-pink-600 flex items-center justify-center text-white font-bold text-xs flex-shrink-0">4b</div>
              <div>
                <h3 className="text-base font-bold text-white">Tully-Fisher at z=0 vs z={z.toFixed(1)}</h3>
                <p className="text-xs text-slate-400">V{"\u2074"} = GMa{"\u2080"}(z) shifts the BTFR</p>
              </div>
            </div>

            <div className="h-[280px]">
              <ResponsiveContainer width="100%" height="100%">
                <LineChart data={btfrData} margin={{ top: 5, right: 15, left: 5, bottom: 5 }}>
                  <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                  <XAxis
                    dataKey="logM" type="number" domain={[7, 12]}
                    stroke="rgba(255,255,255,0.3)"
                    tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 10 }}
                    label={{ value: 'log\u2081\u2080(M/M\u2609)', position: 'insideBottom', offset: -3, fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                  />
                  <YAxis
                    stroke="rgba(255,255,255,0.3)"
                    tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 10 }}
                    label={{ value: 'V_flat (km/s)', angle: -90, position: 'insideLeft', fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                  />
                  <Tooltip
                    contentStyle={{ backgroundColor: 'rgba(15,23,42,0.95)', border: '1px solid rgba(255,255,255,0.1)', borderRadius: '12px', fontSize: '11px' }}
                    formatter={(v: number, name: string) => {
                      const labels: Record<string, string> = { vflat_z0: 'BTFR (z=0)', vflat_zN: 'BTFR (z=' + z.toFixed(1) + ')' };
                      return [v.toFixed(1) + ' km/s', labels[name] || name];
                    }}
                  />
                  <Line type="monotone" dataKey="vflat_z0" stroke="#22d3ee" strokeWidth={2} dot={false} name="z=0" />
                  <Line type="monotone" dataKey="vflat_zN" stroke="#f59e0b" strokeWidth={2.5} dot={false} name={"z=" + z.toFixed(1)} />
                  <Legend wrapperStyle={{ fontSize: '10px' }} />
                </LineChart>
              </ResponsiveContainer>
            </div>

            <p className="text-xs text-slate-400 mt-3">
              The Baryonic Tully-Fisher Relation shifts upward at high z because V{"\u2074"} = GMa{"\u2080"}(z).
              At z={z.toFixed(1)}, V{"\u209B\u2090\u209C"} is {((Math.pow(e, 0.25) - 1) * 100).toFixed(1)}% higher at fixed mass.
              This is a <span className="text-rose-400 font-bold">unique prediction</span> {"\u2014"} MOND predicts no shift.
            </p>
          </GlassCard>
        </div>

        <GlassCard glow="amber">
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-rose-500 to-red-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">5</div>
            <div>
              <h3 className="text-lg font-bold text-white">Falsifiability Decision Table</h3>
              <p className="text-xs text-slate-400">What future observations would confirm or kill each model</p>
            </div>
          </div>

          <div className="overflow-x-auto">
            <table className="w-full text-xs sm:text-sm">
              <thead>
                <tr className="border-b border-white/10">
                  <th className="text-left py-3 px-2 text-slate-400 font-medium">Observation</th>
                  <th className="text-center py-3 px-2 text-violet-400 font-medium">MOND</th>
                  <th className="text-center py-3 px-2 text-slate-400 font-medium">{"\u039B"}CDM</th>
                  <th className="text-center py-3 px-2 text-amber-400 font-medium">Cosmic Floor</th>
                </tr>
              </thead>
              <tbody>
                <tr className="border-b border-white/5">
                  <td className="py-3 px-2 text-slate-300 font-medium">a{"\u2080"} constant at all z</td>
                  <td className="py-3 px-2 text-center"><CheckCircle2 className="w-4 h-4 text-emerald-400 mx-auto" /></td>
                  <td className="py-3 px-2 text-center text-slate-500">{"\u2014"}</td>
                  <td className="py-3 px-2 text-center"><XCircle className="w-4 h-4 text-rose-400 mx-auto" /></td>
                </tr>
                <tr className="border-b border-white/5">
                  <td className="py-3 px-2 text-slate-300 font-medium">a{"\u2080"}(z) = a{"\u2080"}{"\u00D7"}E(z)</td>
                  <td className="py-3 px-2 text-center"><XCircle className="w-4 h-4 text-rose-400 mx-auto" /></td>
                  <td className="py-3 px-2 text-center text-slate-500">{"\u2014"}</td>
                  <td className="py-3 px-2 text-center"><CheckCircle2 className="w-4 h-4 text-emerald-400 mx-auto" /></td>
                </tr>
                <tr className="border-b border-white/5">
                  <td className="py-3 px-2 text-slate-300 font-medium">a{"\u2080"} varies randomly with z</td>
                  <td className="py-3 px-2 text-center"><XCircle className="w-4 h-4 text-rose-400 mx-auto" /></td>
                  <td className="py-3 px-2 text-center"><CheckCircle2 className="w-4 h-4 text-emerald-400 mx-auto" /></td>
                  <td className="py-3 px-2 text-center"><XCircle className="w-4 h-4 text-rose-400 mx-auto" /></td>
                </tr>
                <tr className="border-b border-white/5">
                  <td className="py-3 px-2 text-slate-300 font-medium">BTFR shifts at high z</td>
                  <td className="py-3 px-2 text-center"><XCircle className="w-4 h-4 text-rose-400 mx-auto" /></td>
                  <td className="py-3 px-2 text-center text-slate-500">unclear</td>
                  <td className="py-3 px-2 text-center"><CheckCircle2 className="w-4 h-4 text-emerald-400 mx-auto" /></td>
                </tr>
                <tr className="border-b border-white/5">
                  <td className="py-3 px-2 text-slate-300 font-medium">No a{"\u2080"} scale at all</td>
                  <td className="py-3 px-2 text-center"><XCircle className="w-4 h-4 text-rose-400 mx-auto" /></td>
                  <td className="py-3 px-2 text-center"><CheckCircle2 className="w-4 h-4 text-emerald-400 mx-auto" /></td>
                  <td className="py-3 px-2 text-center"><XCircle className="w-4 h-4 text-rose-400 mx-auto" /></td>
                </tr>
              </tbody>
            </table>
          </div>

          <div className="mt-4 grid grid-cols-1 sm:grid-cols-3 gap-3">
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-3">
              <div className="flex items-center gap-2 mb-1">
                <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                <span className="text-xs font-bold text-emerald-400">Confirms Model</span>
              </div>
              <p className="text-xs text-slate-400">
                If a{"\u2080"}(z=1) {"\u2248"} {(A0_MS2 * Ez(1, omM, omL) * 1e10).toFixed(2)}{"\u00D7"}10{"\u207B\u00B9\u2070"} m/s{"\u00B2"} ({"\u00B1"}20%)
              </p>
            </div>
            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-3">
              <div className="flex items-center gap-2 mb-1">
                <XCircle className="w-4 h-4 text-rose-400" />
                <span className="text-xs font-bold text-rose-400">Kills Model</span>
              </div>
              <p className="text-xs text-slate-400">
                If a{"\u2080"}(z=1) = {(A0_MS2 * 1e10).toFixed(2)}{"\u00D7"}10{"\u207B\u00B9\u2070"} m/s{"\u00B2"} (no change)
              </p>
            </div>
            <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-3">
              <div className="flex items-center gap-2 mb-1">
                <Sparkles className="w-4 h-4 text-amber-400" />
                <span className="text-xs font-bold text-amber-400">When?</span>
              </div>
              <p className="text-xs text-slate-400">
                JWST rotation curve programs (Cycle 3+), ALMA CO kinematics at z{">"}1
              </p>
            </div>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-emerald-500 to-teal-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">6</div>
            <div>
              <h3 className="text-lg font-bold text-white">Quantitative Predictions at z = {z.toFixed(1)}</h3>
              <p className="text-xs text-slate-400">Numbers for the proposal / paper</p>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div className="bg-slate-800/50 border border-white/10 rounded-xl p-4 text-center">
              <div className="text-xs text-slate-500 mb-1">a{"\u2080"}(z={z.toFixed(1)})</div>
              <div className="text-xl font-bold font-mono text-amber-400">{a0z_ms2.toExponential(2)}</div>
              <div className="text-xs text-slate-500">m/s{"\u00B2"}</div>
              <div className="text-xs text-slate-400 mt-1">{a0z.toFixed(0)} (km/s){"\u00B2"}/kpc</div>
            </div>
            <div className="bg-slate-800/50 border border-white/10 rounded-xl p-4 text-center">
              <div className="text-xs text-slate-500 mb-1">E(z={z.toFixed(1)})</div>
              <div className="text-xl font-bold font-mono text-cyan-400">{e.toFixed(4)}</div>
              <div className="text-xs text-slate-500">H(z)/H{"\u2080"}</div>
              <div className="text-xs text-slate-400 mt-1">+{pctIncrease}% above local</div>
            </div>
            <div className="bg-slate-800/50 border border-white/10 rounded-xl p-4 text-center">
              <div className="text-xs text-slate-500 mb-1">BTFR shift</div>
              <div className="text-xl font-bold font-mono text-violet-400">{((Math.pow(e, 0.25) - 1) * 100).toFixed(1)}%</div>
              <div className="text-xs text-slate-500">V{"\u209B\u2090\u209C"} increase at fixed M</div>
              <div className="text-xs text-slate-400 mt-1">V{"\u2074"} {"\u221D"} a{"\u2080"}(z)</div>
            </div>
          </div>
        </GlassCard>

        {feasibility && (
          <>
            <GlassCard glow="amber">
              <div className="flex items-center gap-3 mb-5">
                <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-amber-500 to-red-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">7</div>
                <div>
                  <h3 className="text-lg font-bold text-white">Feasibility Analysis: Can We Actually Test This?</h3>
                  <p className="text-xs text-slate-400">Honest power analysis with real observational constraints</p>
                </div>
              </div>

              <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-6">
                <div className="bg-red-500/10 border border-red-500/20 rounded-xl p-4 text-center">
                  <div className="flex items-center justify-center gap-2 mb-2">
                    <XCircle className="w-4 h-4 text-red-400" />
                    <span className="text-xs font-bold text-red-400">NOW (2024-25)</span>
                  </div>
                  <p className="text-sm font-bold text-white">{feasibility.timeline.now.verdict}</p>
                  <p className="text-xs text-slate-400 mt-1">Errors too large vs signal</p>
                </div>
                <div className="bg-amber-500/10 border border-amber-500/20 rounded-xl p-4 text-center">
                  <div className="flex items-center justify-center gap-2 mb-2">
                    <AlertTriangle className="w-4 h-4 text-amber-400" />
                    <span className="text-xs font-bold text-amber-400">NEAR ({feasibility.timeline.near.when})</span>
                  </div>
                  <p className="text-sm font-bold text-white">{feasibility.timeline.near.verdict}</p>
                  <p className="text-xs text-slate-400 mt-1">JWST NIRSpec IFS</p>
                </div>
                <div className="bg-emerald-500/10 border border-emerald-500/20 rounded-xl p-4 text-center">
                  <div className="flex items-center justify-center gap-2 mb-2">
                    <CheckCircle2 className="w-4 h-4 text-emerald-400" />
                    <span className="text-xs font-bold text-emerald-400">FUTURE ({feasibility.timeline.future.when})</span>
                  </div>
                  <p className="text-sm font-bold text-white">{feasibility.timeline.future.verdict}</p>
                  <p className="text-xs text-slate-400 mt-1">ELT/TMT + ALMA</p>
                </div>
              </div>

              <h4 className="text-sm font-semibold text-white mb-3 flex items-center gap-2">
                <BarChart3 className="w-4 h-4 text-amber-400" />
                Power Analysis: How Many Galaxies at z{"\u223C"}1?
              </h4>
              <p className="text-xs text-slate-400 mb-3">
                Signal: {"\u0394"}V/V = {(feasibility.powerAnalysis.signal_deltaV_frac * 100).toFixed(1)}% at z=1.
                How many galaxies needed to detect this shift?
              </p>
              <div className="overflow-x-auto">
                <table className="w-full text-xs">
                  <thead>
                    <tr className="border-b border-white/10">
                      <th className="text-left py-2 px-3 text-slate-400">{"\u03C3"}(V)/V per galaxy</th>
                      <th className="text-center py-2 px-3 text-slate-400">N for 3{"\u03C3"}</th>
                      <th className="text-center py-2 px-3 text-slate-400">N for 5{"\u03C3"}</th>
                      <th className="text-center py-2 px-3 text-slate-400">Instrument</th>
                    </tr>
                  </thead>
                  <tbody>
                    {feasibility.powerAnalysis.results.map((r, i) => {
                      const instruments = ['ELT (ideal)', 'JWST (best)', 'JWST (typical)', 'KMOS3D', 'Ground AO', 'Seeing-limited'];
                      return (
                        <tr key={i} className="border-b border-white/5">
                          <td className="py-2 px-3 text-slate-300 font-mono">{(r.sigmaV_frac * 100).toFixed(0)}%</td>
                          <td className={`py-2 px-3 text-center font-mono font-bold ${r['N_3\u03C3'] <= 30 ? 'text-emerald-400' : r['N_3\u03C3'] <= 80 ? 'text-amber-400' : 'text-red-400'}`}>{r['N_3\u03C3']}</td>
                          <td className={`py-2 px-3 text-center font-mono font-bold ${r['N_5\u03C3'] <= 50 ? 'text-emerald-400' : r['N_5\u03C3'] <= 150 ? 'text-amber-400' : 'text-red-400'}`}>{r['N_5\u03C3']}</td>
                          <td className="py-2 px-3 text-center text-slate-500">{instruments[i] || ''}</td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>
              <p className="text-xs text-slate-500 mt-3 italic">
                Green = feasible with existing/planned programs. Amber = possible with dedicated survey. Red = impractical.
              </p>
            </GlassCard>

            <GlassCard glow="cyan">
              <div className="flex items-center gap-3 mb-5">
                <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">8</div>
                <div>
                  <h3 className="text-lg font-bold text-white">Existing High-z Surveys</h3>
                  <p className="text-xs text-slate-400">What rotation curve data exists at z {">"} 0.5?</p>
                </div>
              </div>
              <div className="space-y-3">
                {feasibility.existingSurveys.map((s, i) => {
                  const isBest = s.name.includes('JWST');
                  return (
                    <div key={i} className={`p-3 rounded-xl border ${isBest ? 'border-cyan-500/30 bg-cyan-500/5' : 'border-white/5 bg-slate-800/30'}`}>
                      <div className="flex items-center justify-between mb-1">
                        <span className={`text-sm font-bold ${isBest ? 'text-cyan-400' : 'text-white'}`}>{s.name}</span>
                        <span className="text-xs text-slate-500">z = {s.zRange}, N = {s.nGalaxies}</span>
                      </div>
                      <div className="grid grid-cols-2 gap-2 text-xs mb-2">
                        <div><span className="text-slate-500">Method:</span> <span className="text-slate-300">{s.method}</span></div>
                        <div><span className="text-slate-500">Resolution:</span> <span className="text-slate-300">{s.resolution}</span></div>
                      </div>
                      <p className="text-xs text-slate-400">{s.usability}</p>
                    </div>
                  );
                })}
              </div>
            </GlassCard>

            <GlassCard glow="purple">
              <div className="flex items-center gap-3 mb-4">
                <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-purple-500 to-violet-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">9</div>
                <div>
                  <h3 className="text-lg font-bold text-white">Alternative External Tests</h3>
                  <p className="text-xs text-slate-400">If direct RAR at high z is too hard, what else?</p>
                </div>
              </div>
              <div className="space-y-3">
                {feasibility.honestAssessment.alternativeTests.map((t, i) => (
                  <div key={i} className="flex gap-3 items-start">
                    <div className="w-6 h-6 rounded-full bg-purple-500/20 flex items-center justify-center text-purple-400 text-xs font-bold flex-shrink-0 mt-0.5">{i + 1}</div>
                    <p className="text-sm text-slate-300">{t}</p>
                  </div>
                ))}
              </div>
              <div className="mt-4 bg-purple-500/5 border border-purple-500/20 rounded-xl p-3">
                <p className="text-xs text-slate-400">
                  <span className="text-purple-400 font-bold">Main obstacle:</span> {feasibility.honestAssessment.mainObstacle}
                </p>
              </div>
            </GlassCard>
          </>
        )}

        <GlassCard>
          <div className="flex items-center gap-3 mb-4">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-slate-500 to-gray-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">{"\u2139"}</div>
            <div>
              <h3 className="text-lg font-bold text-white">Honest Assessment</h3>
              <p className="text-xs text-slate-400">What this page shows vs. what it doesn't</p>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4">
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-4">
              <h4 className="text-sm font-bold text-emerald-400 mb-2">What this shows</h4>
              <ul className="text-xs text-slate-300 space-y-1">
                <li>{"\u2022"} The Cosmic Floor makes a unique, quantitative, falsifiable prediction</li>
                <li>{"\u2022"} This prediction separates it from both MOND and {"\u039B"}CDM</li>
                <li>{"\u2022"} The numbers are concrete enough for an observing proposal</li>
                <li>{"\u2022"} BTFR evolution is a second, independent test</li>
              </ul>
            </div>
            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <h4 className="text-sm font-bold text-rose-400 mb-2">What this doesn't prove</h4>
              <ul className="text-xs text-slate-300 space-y-1">
                <li>{"\u2022"} Mock data are simulations, not real observations</li>
                <li>{"\u2022"} Selection effects at high z may complicate analysis</li>
                <li>{"\u2022"} The Cosmic Floor may be correct AND still not the full story</li>
                <li>{"\u2022"} No physical mechanism is proposed here {"\u2014"} only the prediction</li>
              </ul>
            </div>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
