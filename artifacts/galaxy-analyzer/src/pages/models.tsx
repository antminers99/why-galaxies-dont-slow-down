import React, { useState } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { useGalaxy, FORMULA_PRESETS } from '@/hooks/use-galaxy';
import { AlertCircle, Play, FileJson, Zap, FlaskConical, CheckCircle2, XCircle, Minus, Loader2 } from 'lucide-react';
import { formatScientific } from '@/lib/utils';
import * as math from 'mathjs';

export default function ModelsPage() {
  const { 
    modelParams, updateModelParams, getInsights, evaluateModel,
    applyPreset, autoOptimize, isOptimizing, optimizationLog,
    datasets, activeDatasetIds
  } = useGalaxy();
  const [formulaInput, setFormulaInput] = useState(modelParams.formula);
  const [error, setError] = useState<string | null>(null);
  const [activeCategory, setActiveCategory] = useState<string | null>(null);
  
  const insights = getInsights();

  const testFormula = (formula: string) => {
    try {
      const testScope = { r: 10, G: modelParams.G, M: modelParams.M, k: modelParams.k, a: modelParams.a };
      const result = math.evaluate(formula, testScope);
      if (typeof result !== 'number' || isNaN(result) || !isFinite(result)) {
        throw new Error("Formula yields invalid result (NaN or Infinity)");
      }
      setError(null);
      updateModelParams({ formula });
    } catch (err: any) {
      setError(err.message || "Invalid mathematical expression");
    }
  };

  const handleFormulaKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      testFormula(formulaInput);
    }
  };

  const handlePresetClick = (preset: typeof FORMULA_PRESETS[0]) => {
    setFormulaInput(preset.formula);
    applyPreset(preset);
    setError(null);
  };

  const exportParams = () => {
    const exportData = {
      formula: modelParams.formula,
      parameters: { G: modelParams.G, M: modelParams.M, k: modelParams.k, a: modelParams.a },
      insights: {
        mseCustom: insights.mseCustom,
        mseNewtonian: insights.mseNewton,
        betterModel: insights.betterModel,
        generalizationScore: insights.generalizationScore,
        perGalaxy: insights.perGalaxy
      },
      timestamp: new Date().toISOString()
    };
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(exportData, null, 2));
    const a = document.createElement('a');
    a.setAttribute("href", dataStr);
    a.setAttribute("download", "model-parameters.json");
    document.body.appendChild(a);
    a.click();
    a.remove();
  };

  const categories = [
    { id: 'additive', label: 'Additive Term', color: 'text-cyan-400 bg-cyan-500/10 border-cyan-500/30' },
    { id: 'modified_gravity', label: 'Modified Gravity', color: 'text-amber-400 bg-amber-500/10 border-amber-500/30' },
    { id: 'transition', label: 'Transition', color: 'text-purple-400 bg-purple-500/10 border-purple-500/30' }
  ];

  const filteredPresets = activeCategory 
    ? FORMULA_PRESETS.filter(p => p.category === activeCategory)
    : FORMULA_PRESETS;

  return (
    <Layout>
      <header className="flex justify-between items-end mb-8">
        <div>
          <h1 className="text-3xl font-bold">Model Builder</h1>
          <p className="text-slate-400 mt-2">Define hypotheses, test formulas, and auto-optimize parameters.</p>
        </div>
        <button 
          onClick={exportParams}
          className="flex items-center gap-2 px-4 py-2 bg-slate-800 hover:bg-slate-700 text-white rounded-xl border border-white/10 transition-colors"
        >
          <FileJson className="w-4 h-4" /> Export JSON
        </button>
      </header>

      <div className="grid grid-cols-1 xl:grid-cols-3 gap-6">
        
        <div className="xl:col-span-1 space-y-6">
          <GlassCard glow="purple" className="flex flex-col">
            <h2 className="text-lg font-semibold mb-4 flex items-center gap-2">
              <FlaskConical className="w-5 h-5 text-purple-400" />
              Custom Formula <span className="text-purple-400 font-mono text-sm">v(r)</span>
            </h2>
            
            <div className="bg-slate-900/80 rounded-xl p-4 border border-white/10 mb-3 focus-within:border-purple-500/50 transition-colors">
              <div className="flex items-center text-slate-400 font-mono mb-2 text-xs">
                <span className="text-purple-400 mr-2">let</span> v = 
              </div>
              <input 
                type="text"
                value={formulaInput}
                onChange={(e) => setFormulaInput(e.target.value)}
                onBlur={() => testFormula(formulaInput)}
                onKeyDown={handleFormulaKeyDown}
                className="w-full bg-transparent outline-none text-lg font-mono text-white placeholder-slate-600"
                spellCheck={false}
              />
            </div>

            {error ? (
              <div className="flex items-start gap-2 text-red-400 text-sm p-3 bg-red-500/10 rounded-lg border border-red-500/20">
                <AlertCircle className="w-4 h-4 shrink-0 mt-0.5" />
                <p>{error}</p>
              </div>
            ) : (
              <div className="flex items-start gap-2 text-emerald-400 text-sm p-3 bg-emerald-500/10 rounded-lg border border-emerald-500/20">
                <Play className="w-4 h-4 shrink-0 mt-0.5" />
                <p>Formula compiles successfully.</p>
              </div>
            )}

            <div className="mt-6 grid grid-cols-2 gap-2 text-xs font-mono">
              <div className="bg-slate-800/50 p-2 rounded border border-white/5"><span className="text-cyan-400">r</span> : radius</div>
              <div className="bg-slate-800/50 p-2 rounded border border-white/5"><span className="text-amber-400">G</span> : gravity</div>
              <div className="bg-slate-800/50 p-2 rounded border border-white/5"><span className="text-green-400">M</span> : mass</div>
              <div className="bg-slate-800/50 p-2 rounded border border-white/5"><span className="text-purple-400">k</span> : custom</div>
              <div className="bg-slate-800/50 p-2 rounded border border-white/5"><span className="text-rose-400">a</span> : core radius</div>
              <div className="bg-slate-800/50 p-2 rounded border border-white/5 text-slate-500">sqrt, log, ^</div>
            </div>
          </GlassCard>

          <GlassCard>
            <h2 className="text-lg font-semibold mb-4">Parameters</h2>
            
            <div className="space-y-6">
              <div>
                <div className="flex justify-between items-end mb-2">
                  <label className="font-mono text-green-400 font-medium text-sm">M <span className="text-slate-500 text-xs font-sans">(Galaxy Mass)</span></label>
                  <span className="font-mono text-xs">{formatScientific(modelParams.M)} M☉</span>
                </div>
                <input 
                  type="range" min={1e9} max={1e12} step={1e9}
                  value={modelParams.M}
                  onChange={(e) => updateModelParams({ M: Number(e.target.value) })}
                  className="w-full accent-green-500"
                />
              </div>

              <div>
                <div className="flex justify-between items-end mb-2">
                  <label className="font-mono text-purple-400 font-medium text-sm">k <span className="text-slate-500 text-xs font-sans">(Dark Matter)</span></label>
                  <span className="font-mono text-xs">{modelParams.k}</span>
                </div>
                <input 
                  type="range" min={0} max={200} step={1}
                  value={modelParams.k}
                  onChange={(e) => updateModelParams({ k: Number(e.target.value) })}
                  className="w-full accent-purple-500"
                />
              </div>

              <div>
                <div className="flex justify-between items-end mb-2">
                  <label className="font-mono text-rose-400 font-medium text-sm">a <span className="text-slate-500 text-xs font-sans">(Core Radius)</span></label>
                  <span className="font-mono text-xs">{modelParams.a.toFixed(1)} kpc</span>
                </div>
                <input 
                  type="range" min={0.1} max={30} step={0.1}
                  value={modelParams.a}
                  onChange={(e) => updateModelParams({ a: Number(e.target.value) })}
                  className="w-full accent-rose-500"
                />
              </div>

              <div className="pt-4 border-t border-white/10 opacity-70 hover:opacity-100 transition-opacity">
                <div className="flex justify-between items-end mb-2">
                  <label className="font-mono text-amber-400 font-medium text-sm">G <span className="text-slate-500 text-xs font-sans">(Grav. Const)</span></label>
                  <span className="font-mono text-[10px]">{modelParams.G.toExponential(4)}</span>
                </div>
                <input 
                  type="range" min={1e-6} max={1e-5} step={1e-7}
                  value={modelParams.G}
                  onChange={(e) => updateModelParams({ G: Number(e.target.value) })}
                  className="w-full accent-amber-500"
                />
              </div>
            </div>

            <button
              onClick={autoOptimize}
              disabled={isOptimizing || activeDatasetIds.length === 0}
              className="mt-6 w-full flex items-center justify-center gap-2 px-4 py-3 bg-gradient-to-r from-purple-600 to-cyan-600 hover:from-purple-500 hover:to-cyan-500 disabled:opacity-50 disabled:cursor-not-allowed text-white font-semibold rounded-xl transition-all shadow-lg shadow-purple-500/20"
            >
              {isOptimizing ? (
                <><Loader2 className="w-5 h-5 animate-spin" /> Optimizing...</>
              ) : (
                <><Zap className="w-5 h-5" /> Auto-Optimize Parameters</>
              )}
            </button>
          </GlassCard>
        </div>

        <div className="xl:col-span-2 space-y-6">
          <GlassCard>
            <h2 className="text-lg font-semibold mb-4">Formula Library</h2>
            <p className="text-sm text-slate-400 mb-4">Pre-built hypothesis formulas. Click to apply.</p>
            
            <div className="flex gap-2 mb-4 flex-wrap">
              <button
                onClick={() => setActiveCategory(null)}
                className={`px-3 py-1.5 rounded-lg text-xs font-medium border transition-colors ${
                  !activeCategory ? 'text-white bg-white/10 border-white/20' : 'text-slate-400 border-transparent hover:bg-white/5'
                }`}
              >
                All
              </button>
              {categories.map(cat => (
                <button
                  key={cat.id}
                  onClick={() => setActiveCategory(activeCategory === cat.id ? null : cat.id)}
                  className={`px-3 py-1.5 rounded-lg text-xs font-medium border transition-colors ${
                    activeCategory === cat.id ? cat.color : 'text-slate-400 border-transparent hover:bg-white/5'
                  }`}
                >
                  {cat.label}
                </button>
              ))}
            </div>

            <div className="grid grid-cols-1 md:grid-cols-2 gap-3">
              {filteredPresets.map(preset => (
                <button
                  key={preset.id}
                  onClick={() => handlePresetClick(preset)}
                  className={`text-left p-4 rounded-xl border transition-all hover:scale-[1.01] ${
                    modelParams.formula === preset.formula 
                      ? 'bg-purple-500/15 border-purple-500/40 shadow-lg shadow-purple-500/10' 
                      : 'bg-slate-800/30 border-white/5 hover:border-white/15 hover:bg-slate-800/50'
                  }`}
                >
                  <div className="flex items-center justify-between mb-1">
                    <span className="font-semibold text-sm">{preset.name}</span>
                    {modelParams.formula === preset.formula && (
                      <span className="text-xs text-purple-400 font-mono">ACTIVE</span>
                    )}
                  </div>
                  <div className="font-mono text-xs text-cyan-400 mb-2">{preset.description}</div>
                  <p className="text-xs text-slate-400 leading-relaxed">{preset.physicalMeaning}</p>
                </button>
              ))}
            </div>
          </GlassCard>

          {insights.perGalaxy.length > 0 && (
            <GlassCard>
              <div className="flex items-center justify-between mb-4">
                <h2 className="text-lg font-semibold">Per-Galaxy Comparison</h2>
                <div className="text-sm">
                  <span className="text-slate-400 mr-2">Generalization:</span>
                  <span className={`font-bold font-mono ${
                    insights.generalizationScore >= 80 ? 'text-emerald-400' :
                    insights.generalizationScore >= 50 ? 'text-amber-400' : 'text-red-400'
                  }`}>
                    {insights.generalizationScore.toFixed(0)}%
                  </span>
                </div>
              </div>
              
              <div className="overflow-x-auto">
                <table className="w-full text-sm">
                  <thead>
                    <tr className="border-b border-white/10">
                      <th className="text-left py-2 px-3 text-slate-400 font-medium">Galaxy</th>
                      <th className="text-right py-2 px-3 text-orange-400 font-medium font-mono">Newton MSE</th>
                      <th className="text-right py-2 px-3 text-purple-400 font-medium font-mono">Custom MSE</th>
                      <th className="text-center py-2 px-3 text-slate-400 font-medium">Winner</th>
                      <th className="text-right py-2 px-3 text-slate-400 font-medium">Points</th>
                    </tr>
                  </thead>
                  <tbody>
                    {insights.perGalaxy.map(g => (
                      <tr key={g.galaxyId} className="border-b border-white/5 hover:bg-white/5 transition-colors">
                        <td className="py-2.5 px-3 font-medium">{g.galaxyName}</td>
                        <td className="py-2.5 px-3 text-right font-mono text-sm">{formatScientific(g.mseNewton)}</td>
                        <td className="py-2.5 px-3 text-right font-mono text-sm">{formatScientific(g.mseCustom)}</td>
                        <td className="py-2.5 px-3 text-center">
                          {g.winner === 'Custom' ? (
                            <span className="inline-flex items-center gap-1 text-emerald-400 text-xs"><CheckCircle2 className="w-3.5 h-3.5" /> Custom</span>
                          ) : g.winner === 'Newtonian' ? (
                            <span className="inline-flex items-center gap-1 text-orange-400 text-xs"><XCircle className="w-3.5 h-3.5" /> Newton</span>
                          ) : (
                            <span className="inline-flex items-center gap-1 text-slate-400 text-xs"><Minus className="w-3.5 h-3.5" /> Tie</span>
                          )}
                        </td>
                        <td className="py-2.5 px-3 text-right text-slate-500 font-mono text-xs">{g.pointCount}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              
              {insights.generalizationScore < 50 && insights.perGalaxy.length > 1 && (
                <div className="mt-4 p-3 bg-amber-500/10 border border-amber-500/20 rounded-lg text-xs text-amber-300">
                  Low generalization score. The model may be overfitting to specific galaxies. Try a different formula or load more datasets.
                </div>
              )}
              {insights.generalizationScore >= 80 && insights.perGalaxy.length >= 3 && (
                <div className="mt-4 p-3 bg-emerald-500/10 border border-emerald-500/20 rounded-lg text-xs text-emerald-300">
                  High generalization! The custom model consistently outperforms Newtonian across galaxies. This formula captures real physics.
                </div>
              )}
            </GlassCard>
          )}

          {optimizationLog.length > 0 && (
            <GlassCard>
              <h2 className="text-lg font-semibold mb-4 flex items-center gap-2">
                <Zap className="w-5 h-5 text-amber-400" />
                Optimization Log
              </h2>
              <div className="bg-slate-900/80 rounded-xl p-4 border border-white/10 max-h-[300px] overflow-y-auto font-mono text-xs leading-relaxed">
                {optimizationLog.map((line, i) => (
                  <div key={i} className={`${
                    line.startsWith('RESULT:') ? 'text-cyan-400 font-bold mt-2' :
                    line.includes('Improvement:') ? 'text-emerald-400' :
                    line.includes('Phase') ? 'text-amber-400 mt-1' :
                    line === '' ? 'h-2' : 'text-slate-300'
                  }`}>
                    {line || '\u00A0'}
                  </div>
                ))}
              </div>
            </GlassCard>
          )}
        </div>
      </div>
    </Layout>
  );
}
