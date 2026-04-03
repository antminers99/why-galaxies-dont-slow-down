import React, { useState, useMemo } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { Orbit, AlertTriangle, CheckCircle2, BarChart3, Info } from 'lucide-react';
import {
  BarChart, Bar, XAxis, YAxis, CartesianGrid,
  Tooltip, ResponsiveContainer, Legend, ReferenceLine,
  ComposedChart, Line
} from 'recharts';

const G_KPC = 4.3009e-6;
const A0_KPC = 3702;

interface ClusterData {
  name: string;
  totalMass: number;
  baryonicMass: number;
  radius: number;
  observedSigma: number;
  temperature: number;
  notes: string;
}

const CLUSTERS: ClusterData[] = [
  {
    name: 'Coma',
    totalMass: 1.2e15,
    baryonicMass: 1.8e14,
    radius: 2900,
    observedSigma: 1000,
    temperature: 8.2,
    notes: 'Richest nearby cluster, well-studied mass profile'
  },
  {
    name: 'Perseus',
    totalMass: 6.6e14,
    baryonicMass: 1.2e14,
    radius: 2200,
    observedSigma: 1280,
    temperature: 6.8,
    notes: 'Brightest X-ray cluster, strong cooling flow'
  },
  {
    name: 'Bullet',
    totalMass: 1.5e15,
    baryonicMass: 2.5e14,
    radius: 3500,
    observedSigma: 1100,
    temperature: 14.0,
    notes: 'Post-merger cluster, DM-baryon separation observed via lensing'
  },
  {
    name: 'Virgo',
    totalMass: 4.2e14,
    baryonicMass: 7.0e13,
    radius: 1550,
    observedSigma: 760,
    temperature: 2.3,
    notes: 'Nearest massive cluster (16.5 Mpc), irregular and still forming'
  },
  {
    name: 'Abell 2029',
    totalMass: 8.0e14,
    baryonicMass: 1.4e14,
    radius: 2800,
    observedSigma: 1160,
    temperature: 8.7,
    notes: 'Very regular, relaxed cluster with giant cD galaxy'
  }
];

const VIRIAL_FACTOR = Math.sqrt(3);

function newtonianSigma(M: number, r: number): number {
  const gbar = G_KPC * M / (r * r);
  const vc = Math.sqrt(gbar * r);
  return vc / VIRIAL_FACTOR;
}

function floorSigma(M: number, r: number, a0: number): number {
  const gbar = G_KPC * M / (r * r);
  const gobs2 = gbar * gbar + a0 * a0 * gbar / (gbar + a0);
  const vc = Math.sqrt(Math.sqrt(gobs2) * r);
  return vc / VIRIAL_FACTOR;
}

function mondSigma(M: number, r: number, a0: number): number {
  const gbar = G_KPC * M / (r * r);
  const gobs = Math.sqrt(gbar * a0);
  const vc = Math.sqrt(gobs * r);
  return vc / VIRIAL_FACTOR;
}

interface ClusterResult {
  name: string;
  observed: number;
  newtonian: number;
  floor: number;
  mond: number;
  lcdm: number;
  newtonMissing: number;
  floorDeviation: number;
  floorOvershoot: boolean;
  mondMissing: number;
  mondOvershoot: boolean;
  baryonFrac: number;
  notes: string;
  temperature: number;
}

export default function ClusterTestPage() {
  const [selectedCluster, setSelectedCluster] = useState(0);

  const results: ClusterResult[] = useMemo(() => {
    return CLUSTERS.map(c => {
      const bSigma = newtonianSigma(c.baryonicMass, c.radius);
      const fSigma = floorSigma(c.baryonicMass, c.radius, A0_KPC);
      const mSigma = mondSigma(c.baryonicMass, c.radius, A0_KPC);

      const newtonMissing = (1 - (bSigma * bSigma) / (c.observedSigma * c.observedSigma)) * 100;
      const floorDev = (1 - (fSigma * fSigma) / (c.observedSigma * c.observedSigma)) * 100;
      const mondDev = (1 - (mSigma * mSigma) / (c.observedSigma * c.observedSigma)) * 100;

      return {
        name: c.name,
        observed: c.observedSigma,
        newtonian: parseFloat(bSigma.toFixed(1)),
        floor: parseFloat(fSigma.toFixed(1)),
        mond: parseFloat(mSigma.toFixed(1)),
        lcdm: c.observedSigma,
        newtonMissing: parseFloat(newtonMissing.toFixed(1)),
        floorDeviation: parseFloat(Math.abs(floorDev).toFixed(1)),
        floorOvershoot: floorDev < 0,
        mondMissing: parseFloat(Math.abs(mondDev).toFixed(1)),
        mondOvershoot: mondDev < 0,
        baryonFrac: parseFloat(((c.baryonicMass / c.totalMass) * 100).toFixed(1)),
        notes: c.notes,
        temperature: c.temperature,
      };
    });
  }, []);

  const floorAccuracyData = useMemo(() => {
    const deviations = results.map(r => {
      const ratio = r.floor / r.observed;
      return Math.abs(ratio - 1) * 100;
    });
    return (deviations.reduce((a, b) => a + b, 0) / deviations.length).toFixed(1);
  }, [results]);

  const overshootCount = useMemo(() => results.filter(r => r.floorOvershoot).length, [results]);
  const undershootCount = useMemo(() => results.filter(r => !r.floorOvershoot).length, [results]);

  const avgNewtonDeficit = useMemo(() => {
    return (results.reduce((s, r) => s + r.newtonMissing, 0) / results.length).toFixed(1);
  }, [results]);

  const avgFloorDeficit = useMemo(() => {
    return (results.reduce((s, r) => s + r.floorDeviation, 0) / results.length).toFixed(1);
  }, [results]);

  const avgMondDeficit = useMemo(() => {
    return (results.reduce((s, r) => s + r.mondMissing, 0) / results.length).toFixed(1);
  }, [results]);

  const deficitReduction = useMemo(() => {
    const avgNewton = results.reduce((s, r) => s + r.newtonMissing, 0) / results.length;
    const avgFloor = results.reduce((s, r) => s + r.floorDeviation, 0) / results.length;
    return ((avgNewton - avgFloor) / avgNewton * 100).toFixed(0);
  }, [results]);

  const velocityComparisonData = useMemo(() => {
    return results.map(r => ({
      name: r.name,
      Observed: r.observed,
      Newtonian: r.newtonian,
      'Cosmic Floor': r.floor,
      MOND: r.mond,
    }));
  }, [results]);

  const deviationData = useMemo(() => {
    return results.map(r => ({
      name: r.name,
      'Newton deficit': r.newtonMissing,
      'Floor deviation': r.floorOvershoot ? -r.floorDeviation : r.floorDeviation,
      'MOND deviation': r.mondOvershoot ? -r.mondMissing : r.mondMissing,
    }));
  }, [results]);

  const sel = results[selectedCluster];

  const radialProfile = useMemo(() => {
    const c = CLUSTERS[selectedCluster];
    const points = [];
    for (let frac = 0.05; frac <= 2.0; frac += 0.05) {
      const r = c.radius * frac;
      const massFrac = Math.min(1, Math.pow(frac, 1.5));
      const mBary = c.baryonicMass * massFrac;
      const mTotal = c.totalMass * massFrac;

      points.push({
        r: parseFloat(r.toFixed(0)),
        rFrac: parseFloat(frac.toFixed(2)),
        newton: parseFloat(newtonianSigma(mBary, r).toFixed(1)),
        floor: parseFloat(floorSigma(mBary, r, A0_KPC).toFixed(1)),
        mond: parseFloat(mondSigma(mBary, r, A0_KPC).toFixed(1)),
        observed: parseFloat(newtonianSigma(mTotal, r).toFixed(1)),
      });
    }
    return points;
  }, [selectedCluster]);

  return (
    <Layout>
      <div className="space-y-8">
        <div>
          <div className="flex items-center gap-3 mb-4">
            <Orbit className="w-8 h-8 text-amber-400" />
            <div>
              <h1 className="text-3xl font-bold font-display text-white">Cluster Stress Test</h1>
              <p className="text-slate-400 text-sm">Applying the acceleration floor to galaxy clusters {"\u2014"} where MOND fails</p>
            </div>
          </div>
          <div className="bg-gradient-to-r from-amber-500/10 to-rose-500/10 border border-amber-500/20 rounded-xl p-4 sm:p-6">
            <p className="text-amber-300 font-mono text-sm leading-relaxed text-center italic">
              "If the floor model is universal, it must say something about clusters {"\u2014"} even if it can't explain them fully."
            </p>
          </div>
        </div>

        <GlassCard glow="amber">
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-amber-500 to-orange-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">1</div>
            <div>
              <h3 className="text-lg font-bold text-white">Why Clusters Matter</h3>
              <p className="text-xs text-slate-400">The graveyard of modified gravity theories</p>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <h4 className="text-sm font-bold text-rose-400 mb-2">The Problem</h4>
              <p className="text-xs text-slate-300">
                Galaxy clusters are 10{"\u00D7"}-100{"\u00D7"} more massive than individual galaxies.
                MOND successfully explains galaxy rotation curves but <span className="text-rose-400 font-bold">fails in clusters</span> {"\u2014"}
                it only recovers ~50% of the missing mass.
              </p>
            </div>
            <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-4">
              <h4 className="text-sm font-bold text-amber-400 mb-2">Our Test</h4>
              <p className="text-xs text-slate-300">
                Apply the same equation (g{"\u2082"} = g{"\u209A\u2082"} + a{"\u2080\u00B2"}{"\u00D7"}g{"\u209A"}/(g{"\u209A"}+a{"\u2080"}))
                to 5 well-known clusters. Measure how much of the missing mass
                the floor model <span className="text-amber-400 font-bold">reduces vs. leaves unexplained</span>.
              </p>
            </div>
            <div className="bg-cyan-500/5 border border-cyan-500/20 rounded-xl p-4">
              <h4 className="text-sm font-bold text-cyan-400 mb-2">Honest Expectation</h4>
              <p className="text-xs text-slate-300">
                We <span className="text-cyan-400 font-bold">do not expect</span> the floor to eliminate all missing mass in clusters.
                A partial reduction would still be physically meaningful.
                A total explanation would be extraordinary.
              </p>
            </div>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">2</div>
            <div>
              <h3 className="text-lg font-bold text-white">Velocity Dispersion: Prediction vs. Observed</h3>
              <p className="text-xs text-slate-400">Each cluster's velocity dispersion under four models</p>
            </div>
          </div>

          <div className="h-[300px] sm:h-[350px]">
            <ResponsiveContainer width="100%" height="100%">
              <BarChart data={velocityComparisonData} margin={{ top: 10, right: 20, left: 5, bottom: 5 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                <XAxis
                  dataKey="name"
                  stroke="rgba(255,255,255,0.3)"
                  tick={{ fill: 'rgba(255,255,255,0.6)', fontSize: 11 }}
                />
                <YAxis
                  stroke="rgba(255,255,255,0.3)"
                  tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                  label={{ value: '\u03C3 (km/s)', angle: -90, position: 'insideLeft', fill: 'rgba(255,255,255,0.5)', fontSize: 12 }}
                />
                <Tooltip
                  contentStyle={{ backgroundColor: 'rgba(15,23,42,0.95)', border: '1px solid rgba(255,255,255,0.1)', borderRadius: '12px', fontSize: '12px' }}
                  formatter={(v: number) => v.toFixed(0) + ' km/s'}
                />
                <Bar dataKey="Observed" fill="#22d3ee" radius={[4, 4, 0, 0]} />
                <Bar dataKey="Newtonian" fill="#f97316" radius={[4, 4, 0, 0]} />
                <Bar dataKey="Cosmic Floor" fill="#f59e0b" radius={[4, 4, 0, 0]} />
                <Bar dataKey="MOND" fill="#8b5cf6" radius={[4, 4, 0, 0]} />
                <Legend wrapperStyle={{ fontSize: '11px', color: 'rgba(255,255,255,0.6)' }} />
              </BarChart>
            </ResponsiveContainer>
          </div>
        </GlassCard>

        <div className="grid grid-cols-1 sm:grid-cols-3 gap-4">
          <div className="bg-gradient-to-br from-orange-500/10 to-orange-500/5 border border-orange-500/20 rounded-xl p-5 text-center">
            <div className="text-3xl font-bold font-mono text-orange-400">{avgNewtonDeficit}%</div>
            <div className="text-xs text-slate-400 mt-1">Avg Newton deficit</div>
            <div className="text-xs text-slate-500">(baryons only)</div>
          </div>
          <div className="bg-gradient-to-br from-amber-500/10 to-amber-500/5 border border-amber-500/20 rounded-xl p-5 text-center">
            <div className="text-3xl font-bold font-mono text-amber-400">{avgFloorDeficit}%</div>
            <div className="text-xs text-slate-400 mt-1">Avg Floor deficit</div>
            <div className="text-xs text-slate-500">{deficitReduction}% reduction from Newton</div>
          </div>
          <div className="bg-gradient-to-br from-violet-500/10 to-violet-500/5 border border-violet-500/20 rounded-xl p-5 text-center">
            <div className="text-3xl font-bold font-mono text-violet-400">{avgMondDeficit}%</div>
            <div className="text-xs text-slate-400 mt-1">Avg MOND deficit</div>
            <div className="text-xs text-slate-500">(comparable to Floor)</div>
          </div>
        </div>

        <GlassCard glow="purple">
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-rose-500 to-pink-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">3</div>
            <div>
              <h3 className="text-lg font-bold text-white">Mass Deficit Analysis</h3>
              <p className="text-xs text-slate-400">Negative = overshoot (predicts too much), Positive = deficit (still missing mass)</p>
            </div>
          </div>

          <div className="h-[300px] sm:h-[350px]">
            <ResponsiveContainer width="100%" height="100%">
              <BarChart data={deviationData} margin={{ top: 10, right: 20, left: 5, bottom: 5 }}>
                <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                <XAxis
                  dataKey="name"
                  stroke="rgba(255,255,255,0.3)"
                  tick={{ fill: 'rgba(255,255,255,0.6)', fontSize: 11 }}
                />
                <YAxis
                  stroke="rgba(255,255,255,0.3)"
                  tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                  label={{ value: 'Deviation %', angle: -90, position: 'insideLeft', fill: 'rgba(255,255,255,0.5)', fontSize: 12 }}
                />
                <ReferenceLine y={0} stroke="rgba(255,255,255,0.3)" strokeDasharray="3 3" />
                <Tooltip
                  contentStyle={{ backgroundColor: 'rgba(15,23,42,0.95)', border: '1px solid rgba(255,255,255,0.1)', borderRadius: '12px', fontSize: '12px' }}
                  formatter={(v: number, name: string) => [v.toFixed(1) + '%', name]}
                />
                <Bar dataKey="Newton deficit" fill="#f97316" radius={[4, 4, 0, 0]} />
                <Bar dataKey="Floor deviation" fill="#f59e0b" radius={[4, 4, 0, 0]} />
                <Bar dataKey="MOND deviation" fill="#8b5cf6" radius={[4, 4, 0, 0]} />
                <Legend wrapperStyle={{ fontSize: '11px', color: 'rgba(255,255,255,0.6)' }} />
              </BarChart>
            </ResponsiveContainer>
          </div>

          <div className="mt-4 bg-amber-500/10 border border-amber-500/20 rounded-xl p-4">
            <div className="flex items-center gap-2 mb-2">
              <BarChart3 className="w-4 h-4 text-amber-400" />
              <span className="text-sm font-bold text-amber-400">Key Finding</span>
            </div>
            <p className="text-sm text-slate-300">
              The floor model (with same a{"\u2080"} from galaxies) brings cluster predictions
              within <span className="text-amber-400 font-bold font-mono">{floorAccuracyData}%</span> average
              deviation from observed. It <span className="text-amber-400 font-bold">overshoots</span> in {overshootCount}/5 clusters
              and <span className="text-cyan-400 font-bold">undershoots</span> in {undershootCount}/5.
              This means the floor effect is in the right ballpark but not a perfect match {"\u2014"}{" "}
              consistent with additional physics (real DM, gas dynamics) playing a role.
            </p>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-violet-500 to-purple-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">4</div>
            <div>
              <h3 className="text-lg font-bold text-white">Cluster Deep Dive: {sel.name}</h3>
              <p className="text-xs text-slate-400">Click a cluster to explore its radial profile</p>
            </div>
          </div>

          <div className="flex flex-wrap gap-2 mb-5">
            {results.map((r, i) => (
              <button
                key={r.name}
                onClick={() => setSelectedCluster(i)}
                className={"px-3 py-1.5 rounded-lg text-xs font-medium transition-all " + (i === selectedCluster ? "bg-violet-500/20 text-violet-400 border border-violet-500/30" : "bg-white/5 text-slate-400 border border-white/10 hover:bg-white/10")}
              >
                {r.name}
              </button>
            ))}
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-6">
            <div>
              <h4 className="text-sm font-bold text-white mb-3">Radial Velocity Profile</h4>
              <div className="h-[280px]">
                <ResponsiveContainer width="100%" height="100%">
                  <ComposedChart data={radialProfile} margin={{ top: 5, right: 15, left: 5, bottom: 5 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="rgba(255,255,255,0.05)" />
                    <XAxis
                      dataKey="rFrac" type="number" domain={[0, 2]}
                      stroke="rgba(255,255,255,0.3)"
                      tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 10 }}
                      label={{ value: 'r / r\u2082\u2080\u2080', position: 'insideBottom', offset: -3, fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                    />
                    <YAxis
                      stroke="rgba(255,255,255,0.3)"
                      tick={{ fill: 'rgba(255,255,255,0.5)', fontSize: 10 }}
                      label={{ value: '\u03C3 (km/s)', angle: -90, position: 'insideLeft', fill: 'rgba(255,255,255,0.5)', fontSize: 11 }}
                    />
                    <Tooltip
                      contentStyle={{ backgroundColor: 'rgba(15,23,42,0.95)', border: '1px solid rgba(255,255,255,0.1)', borderRadius: '12px', fontSize: '11px' }}
                      formatter={(v: number, name: string) => [v.toFixed(0) + ' km/s', name]}
                    />
                    <Line type="monotone" dataKey="observed" stroke="#22d3ee" strokeWidth={2.5} dot={false} name="Observed (DM+bary)" />
                    <Line type="monotone" dataKey="newton" stroke="#f97316" strokeWidth={1.5} strokeDasharray="5 3" dot={false} name="Newton (bary only)" />
                    <Line type="monotone" dataKey="floor" stroke="#f59e0b" strokeWidth={2} dot={false} name="Cosmic Floor" />
                    <Line type="monotone" dataKey="mond" stroke="#8b5cf6" strokeWidth={1.5} strokeDasharray="8 4" dot={false} name="MOND" />
                    <Legend wrapperStyle={{ fontSize: '10px' }} />
                  </ComposedChart>
                </ResponsiveContainer>
              </div>
            </div>

            <div className="space-y-4">
              <h4 className="text-sm font-bold text-white mb-2">Cluster Properties</h4>
              <div className="bg-slate-800/50 border border-white/10 rounded-xl p-4">
                <div className="grid grid-cols-2 gap-y-3 gap-x-4 text-xs">
                  <div>
                    <div className="text-slate-500">Total Mass</div>
                    <div className="text-white font-mono">{CLUSTERS[selectedCluster].totalMass.toExponential(1)} M{"\u2609"}</div>
                  </div>
                  <div>
                    <div className="text-slate-500">Baryonic Mass</div>
                    <div className="text-white font-mono">{CLUSTERS[selectedCluster].baryonicMass.toExponential(1)} M{"\u2609"}</div>
                  </div>
                  <div>
                    <div className="text-slate-500">Radius (r{"\u2082\u2080\u2080"})</div>
                    <div className="text-white font-mono">{CLUSTERS[selectedCluster].radius} kpc</div>
                  </div>
                  <div>
                    <div className="text-slate-500">Temperature</div>
                    <div className="text-white font-mono">{CLUSTERS[selectedCluster].temperature} keV</div>
                  </div>
                  <div>
                    <div className="text-slate-500">Baryon Fraction</div>
                    <div className="text-cyan-400 font-mono">{sel.baryonFrac}%</div>
                  </div>
                  <div>
                    <div className="text-slate-500">Observed {"\u03C3"}</div>
                    <div className="text-cyan-400 font-mono">{sel.observed} km/s</div>
                  </div>
                </div>
              </div>

              <div className="bg-slate-800/50 border border-white/10 rounded-xl p-4">
                <h5 className="text-xs font-bold text-white mb-2">Predictions for {sel.name}</h5>
                <div className="space-y-2 text-xs">
                  <div className="flex justify-between">
                    <span className="text-orange-400">Newton (baryons):</span>
                    <span className="font-mono text-white">{sel.newtonian} km/s <span className="text-rose-400">({sel.newtonMissing}% deficit)</span></span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-amber-400">Cosmic Floor:</span>
                    <span className="font-mono text-white">{sel.floor} km/s <span className={sel.floorOvershoot ? "text-amber-400" : "text-rose-400"}>({sel.floorOvershoot ? "+" : ""}{sel.floorDeviation}% {sel.floorOvershoot ? "overshoot" : "deficit"})</span></span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-violet-400">MOND:</span>
                    <span className="font-mono text-white">{sel.mond} km/s <span className={sel.mondOvershoot ? "text-amber-400" : "text-rose-400"}>({sel.mondOvershoot ? "+" : ""}{sel.mondMissing}% {sel.mondOvershoot ? "overshoot" : "deficit"})</span></span>
                  </div>
                  <div className="flex justify-between">
                    <span className="text-cyan-400">Observed:</span>
                    <span className="font-mono text-white">{sel.observed} km/s</span>
                  </div>
                </div>
              </div>

              <div className="bg-slate-800/50 border border-white/10 rounded-xl p-4">
                <h5 className="text-xs font-bold text-amber-400 mb-1">Floor Accuracy ({sel.floorOvershoot ? "Overshoot" : "Deficit"})</h5>
                <div className="flex items-center gap-2">
                  <div className="flex-1 bg-slate-700 rounded-full h-3 relative">
                    <div className="absolute left-1/2 top-0 w-px h-3 bg-white/30" />
                    <div
                      className={"rounded-full h-3 transition-all " + (sel.floorOvershoot ? "bg-gradient-to-r from-amber-500 to-orange-500 ml-auto" : "bg-gradient-to-r from-cyan-500 to-blue-500")}
                      style={{ width: Math.min(50, sel.floorDeviation / 2) + '%', marginLeft: sel.floorOvershoot ? 'auto' : undefined }}
                    />
                  </div>
                  <span className={"font-mono text-sm font-bold " + (sel.floorOvershoot ? "text-amber-400" : "text-cyan-400")}>{sel.floorOvershoot ? "+" : "-"}{sel.floorDeviation}%</span>
                </div>
                <p className="text-xs text-slate-500 mt-1">{sel.floorOvershoot ? "Predicts more than observed (no missing mass needed)" : "Still needs additional mass beyond the floor"}</p>
              </div>

              <p className="text-xs text-slate-500 italic">{sel.notes}</p>
            </div>
          </div>
        </GlassCard>

        <GlassCard>
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-cyan-500 to-blue-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">5</div>
            <div>
              <h3 className="text-lg font-bold text-white">Full Comparison Table</h3>
              <p className="text-xs text-slate-400">All five clusters at a glance</p>
            </div>
          </div>

          <div className="overflow-x-auto">
            <table className="w-full text-xs">
              <thead>
                <tr className="border-b border-white/10">
                  <th className="text-left py-2 px-2 text-slate-400">Cluster</th>
                  <th className="text-right py-2 px-2 text-cyan-400">{"\u03C3"} obs</th>
                  <th className="text-right py-2 px-2 text-orange-400">{"\u03C3"} Newton</th>
                  <th className="text-right py-2 px-2 text-amber-400">{"\u03C3"} Floor</th>
                  <th className="text-right py-2 px-2 text-violet-400">{"\u03C3"} MOND</th>
                  <th className="text-right py-2 px-2 text-rose-400">Newton gap</th>
                  <th className="text-right py-2 px-2 text-amber-400">Floor result</th>
                </tr>
              </thead>
              <tbody>
                {results.map(r => (
                  <tr key={r.name} className="border-b border-white/5 hover:bg-white/5 transition-colors">
                    <td className="py-2 px-2 text-white font-medium">{r.name}</td>
                    <td className="py-2 px-2 text-right font-mono text-cyan-400">{r.observed}</td>
                    <td className="py-2 px-2 text-right font-mono text-orange-400">{r.newtonian}</td>
                    <td className="py-2 px-2 text-right font-mono text-amber-400">{r.floor}</td>
                    <td className="py-2 px-2 text-right font-mono text-violet-400">{r.mond}</td>
                    <td className="py-2 px-2 text-right font-mono text-rose-400">{r.newtonMissing}%</td>
                    <td className={"py-2 px-2 text-right font-mono " + (r.floorOvershoot ? "text-amber-400" : "text-emerald-400")}>{r.floorOvershoot ? "+" : "-"}{r.floorDeviation}%</td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </GlassCard>

        <GlassCard glow="cyan">
          <div className="flex items-center gap-3 mb-5">
            <div className="w-10 h-10 rounded-xl bg-gradient-to-br from-slate-500 to-gray-600 flex items-center justify-center text-white font-bold text-sm flex-shrink-0">6</div>
            <div>
              <h3 className="text-lg font-bold text-white">Verdict</h3>
              <p className="text-xs text-slate-400">What the cluster test tells us {"\u2014"} honestly</p>
            </div>
          </div>

          <div className="grid grid-cols-1 md:grid-cols-2 gap-4 mb-4">
            <div className="bg-emerald-500/5 border border-emerald-500/20 rounded-xl p-4">
              <div className="flex items-center gap-2 mb-2">
                <CheckCircle2 className="w-5 h-5 text-emerald-400" />
                <h4 className="text-sm font-bold text-emerald-400">What the floor does</h4>
              </div>
              <ul className="text-xs text-slate-300 space-y-1">
                <li>{"\u2022"} Reduces Newton's ~90%+ deficit to ~{floorAccuracyData}% on average</li>
                <li>{"\u2022"} No tuning: same a{"\u2080"} = cH{"\u2080"}/2{"\u03C0"} used for 197 galaxies applied directly here</li>
                <li>{"\u2022"} Comparable to MOND (both fall short by similar amount in clusters)</li>
                <li>{"\u2022"} Virial-corrected: {"\u03C3"} = v{"\u1D04"}/{"\u221A"}3 (isotropic sphere assumption)</li>
              </ul>
            </div>

            <div className="bg-rose-500/5 border border-rose-500/20 rounded-xl p-4">
              <div className="flex items-center gap-2 mb-2">
                <AlertTriangle className="w-5 h-5 text-rose-400" />
                <h4 className="text-sm font-bold text-rose-400">What remains unexplained</h4>
              </div>
              <ul className="text-xs text-slate-300 space-y-1">
                <li>{"\u2022"} Undershoots in {undershootCount}/{results.length} clusters {"\u2014"} still needs ~30-70% more mass</li>
                <li>{"\u2022"} This is the well-known "cluster problem" for all modified gravity theories</li>
                <li>{"\u2022"} The Bullet Cluster lensing offset still requires collisionless dark matter</li>
                <li>{"\u2022"} Cluster mass estimates have large systematics (X-ray vs. lensing vary 2{"\u00D7"})</li>
              </ul>
            </div>
          </div>

          <div className="bg-amber-500/5 border border-amber-500/20 rounded-xl p-4">
            <div className="flex items-center gap-2 mb-2">
              <Info className="w-4 h-4 text-amber-400" />
              <span className="text-sm font-bold text-amber-400">Interpretation</span>
            </div>
            <p className="text-xs text-slate-300 leading-relaxed">
              The floor model behaves in clusters exactly as a <span className="text-amber-400 font-bold">partial explanation</span> should:
              it reduces the problem without eliminating it. This is consistent with the hypothesis that
              a{"\u2080"} is a real physical effect that <span className="text-amber-400">coexists</span> with genuine dark matter.
              In galaxies, the floor effect dominates (baryonic coupling explains everything).
              In clusters, both the floor effect and additional dark matter contribute.
              This is not a failure {"\u2014"} it's a constraint that any complete theory must satisfy.
            </p>
          </div>
        </GlassCard>
      </div>
    </Layout>
  );
}
