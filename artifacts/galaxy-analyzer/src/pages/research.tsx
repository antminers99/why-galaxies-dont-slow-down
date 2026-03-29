import React, { useState } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { useGalaxy } from '@/hooks/use-galaxy';
import { Loader2, Microscope, Trophy, BarChart3, AlertTriangle, Info, CheckCircle2, XCircle, Minus, FileJson, ArrowUpDown, MapPin } from 'lucide-react';
import { formatScientific } from '@/lib/utils';

export default function ResearchPage() {
  const { 
    datasets, activeDatasetIds, loadAllSamples,
    benchmarkResult, isBenchmarking, runFullBenchmark, modelParams
  } = useGalaxy();

  const [expandedFormula, setExpandedFormula] = useState<string | null>(null);
  const [sortBy, setSortBy] = useState<'improvement' | 'consistency'>('improvement');

  const activeCount = datasets.filter(d => activeDatasetIds.includes(d.id)).length;

  const sortedFormulas = benchmarkResult?.formulas
    ? [...benchmarkResult.formulas].sort((a, b) => 
        sortBy === 'improvement' ? b.avgImprovement - a.avgImprovement : a.kStdDev - b.kStdDev
      )
    : [];

  const exportBenchmark = () => {
    if (!benchmarkResult) return;
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(benchmarkResult, null, 2));
    const a = document.createElement('a');
    a.setAttribute("href", dataStr);
    a.setAttribute("download", "benchmark-results.json");
    document.body.appendChild(a);
    a.click();
    a.remove();
  };

  return (
    <Layout>
      <header className="flex justify-between items-end mb-8">
        <div>
          <h1 className="text-3xl font-bold">Research Lab</h1>
          <p className="text-slate-400 mt-2">Systematic benchmark, regional analysis, and parameter consistency testing.</p>
        </div>
        {benchmarkResult && (
          <button 
            onClick={exportBenchmark}
            className="flex items-center gap-2 px-4 py-2 bg-slate-800 hover:bg-slate-700 text-white rounded-xl border border-white/10 transition-colors"
          >
            <FileJson className="w-4 h-4" /> Export Results
          </button>
        )}
      </header>

      <div className="space-y-6">

        <GlassCard>
          <h2 className="text-lg font-semibold mb-3 flex items-center gap-2">
            <Info className="w-5 h-5 text-cyan-400" />
            Unit System & Physical Constants
          </h2>
          <p className="text-sm text-slate-400 mb-4">
            All calculations use consistent astrophysical units. Understanding these is essential before interpreting any fitted parameter as physically meaningful.
          </p>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-3">
            <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
              <div className="font-mono text-cyan-400 text-sm font-bold">r — Radius</div>
              <div className="text-xs text-slate-300 mt-1">kiloparsecs (kpc)</div>
              <div className="text-xs text-slate-500 mt-0.5">1 kpc = 3,260 light-years = 3.086 x 10^19 m</div>
            </div>
            <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
              <div className="font-mono text-cyan-400 text-sm font-bold">v — Velocity</div>
              <div className="text-xs text-slate-300 mt-1">kilometers per second (km/s)</div>
              <div className="text-xs text-slate-500 mt-0.5">Circular orbital velocity at radius r</div>
            </div>
            <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
              <div className="font-mono text-amber-400 text-sm font-bold">G — Gravitational Constant</div>
              <div className="text-xs text-slate-300 mt-1">{modelParams.G.toExponential(4)} kpc (km/s)^2 / M_sun</div>
              <div className="text-xs text-slate-500 mt-0.5">Real physical value in these units. Not arbitrary.</div>
            </div>
            <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
              <div className="font-mono text-green-400 text-sm font-bold">M — Galaxy Mass</div>
              <div className="text-xs text-slate-300 mt-1">Solar masses (M_sun)</div>
              <div className="text-xs text-slate-500 mt-0.5">Fitted parameter — represents total baryonic mass enclosed.</div>
            </div>
            <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
              <div className="font-mono text-purple-400 text-sm font-bold">k — Dark Matter Parameter</div>
              <div className="text-xs text-slate-300 mt-1">Units depend on formula</div>
              <div className="text-xs text-slate-500 mt-0.5">In sqrt(GM/r + kr): k has units (km/s)^2/kpc. Fitted, not yet physical.</div>
            </div>
            <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
              <div className="font-mono text-rose-400 text-sm font-bold">a — Core Radius</div>
              <div className="text-xs text-slate-300 mt-1">kiloparsecs (kpc)</div>
              <div className="text-xs text-slate-500 mt-0.5">Scale length where gravity softening or halo transition occurs.</div>
            </div>
          </div>
          <div className="mt-4 p-3 bg-amber-500/10 border border-amber-500/20 rounded-lg">
            <div className="flex items-start gap-2">
              <AlertTriangle className="w-4 h-4 text-amber-400 mt-0.5 shrink-0" />
              <div className="text-xs text-amber-200">
                <strong>Important:</strong> G is the real gravitational constant expressed in kpc-(km/s)^2/M_sun units. M is a physical mass estimate. However, k and a are fitting parameters — their values are meaningful <em>within this unit system</em> but should not be treated as universal physical constants until validated against independent measurements.
                If k is consistent across many galaxies, it suggests a real physical effect. If it varies wildly, it may be absorbing galaxy-specific structure.
              </div>
            </div>
          </div>
        </GlassCard>

        <GlassCard glow="purple">
          <div className="flex items-center justify-between mb-4">
            <div>
              <h2 className="text-lg font-semibold flex items-center gap-2">
                <Microscope className="w-5 h-5 text-purple-400" />
                Full Benchmark
              </h2>
              <p className="text-sm text-slate-400 mt-1">
                Test all 7 formula presets across your {activeCount} active galaxies. Each formula is independently optimized per galaxy. Newtonian baseline optimized independently for fair comparison.
              </p>
            </div>
            <div className="flex items-center gap-3">
              <button
                onClick={runFullBenchmark}
                disabled={isBenchmarking || activeCount === 0}
                className="flex items-center gap-2 px-5 py-2.5 bg-gradient-to-r from-purple-600 to-cyan-600 hover:from-purple-500 hover:to-cyan-500 disabled:opacity-50 disabled:cursor-not-allowed text-white font-semibold rounded-xl transition-all shadow-lg shadow-purple-500/20"
              >
                {isBenchmarking ? (
                  <><Loader2 className="w-5 h-5 animate-spin" /> Running Benchmark...</>
                ) : (
                  <><Microscope className="w-5 h-5" /> Run Benchmark</>
                )}
              </button>
            </div>
          </div>
          {activeCount === 0 && (
            <p className="text-amber-400 text-sm">No active galaxies. Load sample datasets first.</p>
          )}
        </GlassCard>

        {benchmarkResult && (
          <>
            <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
              <GlassCard className="p-4" glow="cyan">
                <div className="text-xs text-slate-400 uppercase tracking-wider mb-2 font-semibold flex items-center gap-1">
                  <Trophy className="w-4 h-4" /> Best Formula
                </div>
                <div className="text-xl font-bold text-cyan-400">{benchmarkResult.bestFormula}</div>
                <div className="text-sm text-slate-400 mt-1">
                  {sortedFormulas[0]?.avgImprovement.toFixed(1)}% avg improvement over Newtonian
                </div>
              </GlassCard>
              <GlassCard className="p-4">
                <div className="text-xs text-slate-400 uppercase tracking-wider mb-2 font-semibold flex items-center gap-1">
                  <ArrowUpDown className="w-4 h-4" /> k Consistency
                </div>
                <div className="text-sm text-slate-200 leading-relaxed">{benchmarkResult.kConsistency}</div>
              </GlassCard>
              <GlassCard className="p-4">
                <div className="text-xs text-slate-400 uppercase tracking-wider mb-2 font-semibold flex items-center gap-1">
                  <BarChart3 className="w-4 h-4" /> Galaxies Tested
                </div>
                <div className="text-xl font-bold font-mono">{activeCount}</div>
                <div className="text-sm text-slate-400 mt-1">
                  Across {benchmarkResult.formulas.length} formula variants
                </div>
              </GlassCard>
            </div>

            <GlassCard>
              <div className="flex items-center justify-between mb-4">
                <h2 className="text-lg font-semibold">Formula Rankings</h2>
                <div className="flex gap-2">
                  <button
                    onClick={() => setSortBy('improvement')}
                    className={`px-3 py-1.5 rounded-lg text-xs font-medium border transition-colors ${
                      sortBy === 'improvement' ? 'text-white bg-white/10 border-white/20' : 'text-slate-400 border-transparent hover:bg-white/5'
                    }`}
                  >
                    By Improvement
                  </button>
                  <button
                    onClick={() => setSortBy('consistency')}
                    className={`px-3 py-1.5 rounded-lg text-xs font-medium border transition-colors ${
                      sortBy === 'consistency' ? 'text-white bg-white/10 border-white/20' : 'text-slate-400 border-transparent hover:bg-white/5'
                    }`}
                  >
                    By k Consistency
                  </button>
                </div>
              </div>
              <div className="overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b border-white/10">
                      <th className="text-left py-2 px-3 text-slate-400 font-medium">#</th>
                      <th className="text-left py-2 px-3 text-slate-400 font-medium">Formula</th>
                      <th className="text-right py-2 px-3 text-slate-400 font-medium">Avg MSE</th>
                      <th className="text-right py-2 px-3 text-slate-400 font-medium">Avg Improvement</th>
                      <th className="text-right py-2 px-3 text-slate-400 font-medium">Wins</th>
                      <th className="text-right py-2 px-3 text-purple-400 font-medium">Avg k</th>
                      <th className="text-right py-2 px-3 text-purple-400 font-medium">k StdDev</th>
                      <th className="text-right py-2 px-3 text-rose-400 font-medium">Avg a</th>
                    </tr>
                  </thead>
                  <tbody>
                    {sortedFormulas.map((f, i) => (
                      <tr 
                        key={f.formulaId} 
                        className={`border-b border-white/5 cursor-pointer transition-colors ${
                          expandedFormula === f.formulaId ? 'bg-purple-500/10' : 'hover:bg-white/5'
                        }`}
                        onClick={() => setExpandedFormula(expandedFormula === f.formulaId ? null : f.formulaId)}
                      >
                        <td className="py-2.5 px-3 font-mono text-slate-500">{i + 1}</td>
                        <td className="py-2.5 px-3">
                          <div className="font-medium">{f.formulaName}</div>
                          <div className="text-xs font-mono text-cyan-400/70">{f.formula}</div>
                        </td>
                        <td className="py-2.5 px-3 text-right font-mono text-sm">{formatScientific(f.avgMSE)}</td>
                        <td className="py-2.5 px-3 text-right">
                          <span className={`font-mono font-medium ${f.avgImprovement > 50 ? 'text-emerald-400' : f.avgImprovement > 20 ? 'text-cyan-400' : f.avgImprovement > 0 ? 'text-amber-400' : 'text-red-400'}`}>
                            {f.avgImprovement.toFixed(1)}%
                          </span>
                        </td>
                        <td className="py-2.5 px-3 text-right font-mono">{f.winsCount}/{activeCount}</td>
                        <td className="py-2.5 px-3 text-right font-mono text-purple-300">{f.avgK.toFixed(1)}</td>
                        <td className="py-2.5 px-3 text-right font-mono text-purple-300/60">{f.kStdDev.toFixed(1)}</td>
                        <td className="py-2.5 px-3 text-right font-mono text-rose-300">{f.avgA.toFixed(1)}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </GlassCard>

            {expandedFormula && (() => {
              const formula = benchmarkResult.formulas.find(f => f.formulaId === expandedFormula);
              if (!formula) return null;
              return (
                <GlassCard>
                  <div className="flex items-center justify-between mb-4">
                    <h2 className="text-lg font-semibold">
                      Per-Galaxy Results: <span className="text-purple-400">{formula.formulaName}</span>
                    </h2>
                    <div className="text-sm text-slate-400">
                      Avg k = <span className="text-purple-400 font-mono">{formula.avgK.toFixed(1)}</span> 
                      {' '}(StdDev: <span className="font-mono">{formula.kStdDev.toFixed(1)}</span>)
                    </div>
                  </div>
                  <div className="overflow-x-auto">
                    <table className="w-full text-sm">
                      <thead>
                        <tr className="border-b border-white/10">
                          <th className="text-left py-2 px-3 text-slate-400 font-medium">Galaxy</th>
                          <th className="text-right py-2 px-3 text-orange-400 font-medium">Newton MSE</th>
                          <th className="text-right py-2 px-3 text-purple-400 font-medium">Custom MSE</th>
                          <th className="text-right py-2 px-3 text-slate-400 font-medium">Improvement</th>
                          <th className="text-right py-2 px-3 text-purple-400 font-medium">k</th>
                          <th className="text-right py-2 px-3 text-green-400 font-medium">M</th>
                          <th className="text-right py-2 px-3 text-rose-400 font-medium">a</th>
                          <th className="text-center py-2 px-3 text-slate-400 font-medium">
                            <span className="flex items-center gap-1 justify-center"><MapPin className="w-3 h-3" /> Region</span>
                          </th>
                          <th className="text-right py-2 px-3 text-slate-400 font-medium">Pts</th>
                        </tr>
                      </thead>
                      <tbody>
                        {formula.galaxyResults.map(g => (
                          <tr key={g.galaxyId} className="border-b border-white/5 hover:bg-white/5 transition-colors">
                            <td className="py-2.5 px-3 font-medium">{g.galaxyName}</td>
                            <td className="py-2.5 px-3 text-right font-mono text-xs">{formatScientific(g.mseNewton)}</td>
                            <td className="py-2.5 px-3 text-right font-mono text-xs">{formatScientific(g.mseCustom)}</td>
                            <td className="py-2.5 px-3 text-right">
                              <span className={`font-mono font-medium ${g.improvementPct > 50 ? 'text-emerald-400' : g.improvementPct > 20 ? 'text-cyan-400' : g.improvementPct > 0 ? 'text-amber-400' : 'text-red-400'}`}>
                                {g.improvementPct.toFixed(1)}%
                              </span>
                            </td>
                            <td className="py-2.5 px-3 text-right font-mono text-purple-300">{g.bestK}</td>
                            <td className="py-2.5 px-3 text-right font-mono text-green-300 text-xs">{g.bestM.toExponential(1)}</td>
                            <td className="py-2.5 px-3 text-right font-mono text-rose-300">{g.bestA.toFixed(1)}</td>
                            <td className="py-2.5 px-3 text-center">
                              <div className="flex items-center gap-1 justify-center text-xs">
                                <span className={g.innerImprovement > 20 ? 'text-emerald-400' : g.innerImprovement > 0 ? 'text-slate-300' : 'text-red-400'}>
                                  Inner {g.innerImprovement.toFixed(0)}%
                                </span>
                                <span className="text-slate-600">|</span>
                                <span className={g.outerImprovement > 20 ? 'text-emerald-400' : g.outerImprovement > 0 ? 'text-slate-300' : 'text-red-400'}>
                                  Outer {g.outerImprovement.toFixed(0)}%
                                </span>
                              </div>
                            </td>
                            <td className="py-2.5 px-3 text-right text-slate-500 font-mono text-xs">{g.pointCount}</td>
                          </tr>
                        ))}
                      </tbody>
                    </table>
                  </div>

                  <div className="mt-4 grid grid-cols-1 md:grid-cols-2 gap-4">
                    <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
                      <h4 className="text-xs text-slate-400 uppercase tracking-wider mb-2 font-semibold">k Value Distribution</h4>
                      <div className="space-y-1">
                        {formula.galaxyResults.map(g => {
                          const maxK = Math.max(...formula.galaxyResults.map(x => x.bestK), 1);
                          const widthPct = (g.bestK / maxK) * 100;
                          return (
                            <div key={g.galaxyId} className="flex items-center gap-2 text-xs">
                              <span className="w-28 text-slate-400 truncate">{g.galaxyName}</span>
                              <div className="flex-1 bg-slate-800 rounded-full h-4 overflow-hidden">
                                <div 
                                  className="h-full bg-gradient-to-r from-purple-600 to-purple-400 rounded-full flex items-center justify-end pr-1"
                                  style={{ width: `${Math.max(widthPct, 5)}%` }}
                                >
                                  <span className="text-[10px] font-mono text-white">{g.bestK}</span>
                                </div>
                              </div>
                            </div>
                          );
                        })}
                      </div>
                    </div>
                    <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
                      <h4 className="text-xs text-slate-400 uppercase tracking-wider mb-2 font-semibold">Regional Improvement</h4>
                      <p className="text-xs text-slate-500 mb-3">Does the improvement come from inner regions (near center) or outer regions (far from center)?</p>
                      {(() => {
                        const avgInner = formula.galaxyResults.reduce((s, g) => s + g.innerImprovement, 0) / formula.galaxyResults.length;
                        const avgOuter = formula.galaxyResults.reduce((s, g) => s + g.outerImprovement, 0) / formula.galaxyResults.length;
                        return (
                          <div className="space-y-3">
                            <div>
                              <div className="flex justify-between text-xs mb-1">
                                <span className="text-slate-300">Inner Region (near center)</span>
                                <span className={`font-mono ${avgInner > 20 ? 'text-emerald-400' : 'text-slate-400'}`}>{avgInner.toFixed(1)}%</span>
                              </div>
                              <div className="bg-slate-800 rounded-full h-3">
                                <div className="h-full bg-gradient-to-r from-cyan-600 to-cyan-400 rounded-full" style={{ width: `${Math.min(Math.max(avgInner, 0), 100)}%` }} />
                              </div>
                            </div>
                            <div>
                              <div className="flex justify-between text-xs mb-1">
                                <span className="text-slate-300">Outer Region (far from center)</span>
                                <span className={`font-mono ${avgOuter > 20 ? 'text-emerald-400' : 'text-slate-400'}`}>{avgOuter.toFixed(1)}%</span>
                              </div>
                              <div className="bg-slate-800 rounded-full h-3">
                                <div className="h-full bg-gradient-to-r from-amber-600 to-amber-400 rounded-full" style={{ width: `${Math.min(Math.max(avgOuter, 0), 100)}%` }} />
                              </div>
                            </div>
                            <p className="text-xs text-slate-500 mt-2">
                              {avgOuter > avgInner * 1.5 
                                ? "Improvement concentrated at outer radii — consistent with the galaxy rotation problem where Newtonian gravity fails at large distances."
                                : avgInner > avgOuter * 1.5
                                  ? "Improvement concentrated at inner radii — the formula corrects core dynamics more than outer halo."
                                  : "Improvement distributed across the full curve — the formula captures both inner and outer dynamics."
                              }
                            </p>
                          </div>
                        );
                      })()}
                    </div>
                  </div>
                </GlassCard>
              );
            })()}
          </>
        )}
      </div>
    </Layout>
  );
}
