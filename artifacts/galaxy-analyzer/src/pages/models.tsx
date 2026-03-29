import React, { useState, useEffect } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { useGalaxy } from '@/hooks/use-galaxy';
import { Save, AlertCircle, Play, FileJson } from 'lucide-react';
import { formatScientific } from '@/lib/utils';

export default function ModelsPage() {
  const { modelParams, updateModelParams, getInsights, evaluateModel } = useGalaxy();
  const [formulaInput, setFormulaInput] = useState(modelParams.formula);
  const [error, setError] = useState<string | null>(null);
  
  const insights = getInsights();

  // Test formula validity
  const testFormula = (formula: string) => {
    try {
      // Create a test scope
      const scope = { r: 10, ...modelParams, formula: undefined };
      const testVal = evaluateModel(10, 'custom');
      if (testVal === null) throw new Error("Invalid formula syntax");
      if (isNaN(testVal)) throw new Error("Formula yields NaN");
      setError(null);
      updateModelParams({ formula });
    } catch (err: any) {
      setError(err.message || "Invalid mathematical expression");
    }
  };

  const handleFormulaBlur = () => {
    testFormula(formulaInput);
  };

  const handleFormulaKeyDown = (e: React.KeyboardEvent) => {
    if (e.key === 'Enter') {
      testFormula(formulaInput);
    }
  };

  const exportParams = () => {
    const dataStr = "data:text/json;charset=utf-8," + encodeURIComponent(JSON.stringify(modelParams, null, 2));
    const downloadAnchorNode = document.createElement('a');
    downloadAnchorNode.setAttribute("href",     dataStr);
    downloadAnchorNode.setAttribute("download", "model-parameters.json");
    document.body.appendChild(downloadAnchorNode);
    downloadAnchorNode.click();
    downloadAnchorNode.remove();
  };

  return (
    <Layout>
      <header className="flex justify-between items-end mb-8">
        <div>
          <h1 className="text-3xl font-bold">Model Builder</h1>
          <p className="text-slate-400 mt-2">Define custom physics and tune parameters.</p>
        </div>
        <button 
          onClick={exportParams}
          className="flex items-center gap-2 px-4 py-2 bg-slate-800 hover:bg-slate-700 text-white rounded-xl border border-white/10 transition-colors"
        >
          <FileJson className="w-4 h-4" /> Export JSON
        </button>
      </header>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-8">
        
        {/* Formula Builder */}
        <GlassCard glow="purple" className="flex flex-col">
          <h2 className="text-xl font-semibold mb-6 flex items-center gap-2">
            Velocity Formula <span className="text-purple-400 font-mono">v(r)</span>
          </h2>
          
          <div className="bg-slate-900/80 rounded-xl p-4 border border-white/10 mb-4 focus-within:border-purple-500/50 transition-colors">
            <div className="flex items-center text-slate-400 font-mono mb-2 text-sm">
              <span className="text-purple-400 mr-2">let</span> v = 
            </div>
            <input 
              type="text"
              value={formulaInput}
              onChange={(e) => setFormulaInput(e.target.value)}
              onBlur={handleFormulaBlur}
              onKeyDown={handleFormulaKeyDown}
              className="w-full bg-transparent outline-none text-xl font-mono text-white placeholder-slate-600"
              spellCheck={false}
            />
          </div>

          {error ? (
            <div className="flex items-start gap-2 text-red-400 text-sm mt-2 p-3 bg-red-500/10 rounded-lg border border-red-500/20">
              <AlertCircle className="w-5 h-5 shrink-0" />
              <p>{error}</p>
            </div>
          ) : (
            <div className="flex items-start gap-2 text-emerald-400 text-sm mt-2 p-3 bg-emerald-500/10 rounded-lg border border-emerald-500/20">
              <Play className="w-5 h-5 shrink-0" />
              <p>Formula compiles successfully.</p>
            </div>
          )}

          <div className="mt-8 space-y-3">
            <h3 className="text-sm font-semibold text-slate-400 uppercase tracking-wider">Available Variables</h3>
            <div className="grid grid-cols-2 gap-2 text-sm font-mono">
              <div className="bg-slate-800/50 p-2 rounded border border-white/5"><span className="text-cyan-400">r</span> : radius</div>
              <div className="bg-slate-800/50 p-2 rounded border border-white/5"><span className="text-amber-400">G</span> : gravity</div>
              <div className="bg-slate-800/50 p-2 rounded border border-white/5"><span className="text-green-400">M</span> : mass</div>
              <div className="bg-slate-800/50 p-2 rounded border border-white/5"><span className="text-purple-400">k</span> : custom</div>
            </div>
            <p className="text-xs text-slate-500 mt-2">Uses mathjs syntax. Supports sqrt(), ^, log(), etc.</p>
          </div>
        </GlassCard>

        {/* Parameters */}
        <div className="space-y-6">
          <GlassCard>
            <h2 className="text-xl font-semibold mb-6">Model Parameters</h2>
            
            <div className="space-y-8">
              {/* Parameter M */}
              <div>
                <div className="flex justify-between items-end mb-2">
                  <label className="font-mono text-green-400 font-medium">M <span className="text-slate-500 text-sm font-sans">(Galaxy Mass)</span></label>
                  <span className="font-mono text-sm">{formatScientific(modelParams.M)} M☉</span>
                </div>
                <input 
                  type="range" 
                  min={1e9} 
                  max={1e12} 
                  step={1e9}
                  value={modelParams.M}
                  onChange={(e) => updateModelParams({ M: Number(e.target.value) })}
                  className="w-full accent-green-500"
                />
              </div>

              {/* Parameter k */}
              <div>
                <div className="flex justify-between items-end mb-2">
                  <label className="font-mono text-purple-400 font-medium">k <span className="text-slate-500 text-sm font-sans">(Dark Matter factor)</span></label>
                  <span className="font-mono text-sm">{modelParams.k}</span>
                </div>
                <input 
                  type="range" 
                  min={0} 
                  max={200} 
                  step={1}
                  value={modelParams.k}
                  onChange={(e) => updateModelParams({ k: Number(e.target.value) })}
                  className="w-full accent-purple-500"
                />
              </div>

              {/* Parameter G (Usually constant, but allowed to tweak for fun) */}
              <div className="pt-6 border-t border-white/10 opacity-70 hover:opacity-100 transition-opacity">
                <div className="flex justify-between items-end mb-2">
                  <label className="font-mono text-amber-400 font-medium">G <span className="text-slate-500 text-sm font-sans">(Gravitational Const)</span></label>
                  <span className="font-mono text-xs">{modelParams.G.toExponential(4)}</span>
                </div>
                <input 
                  type="range" 
                  min={1e-6} 
                  max={1e-5} 
                  step={1e-7}
                  value={modelParams.G}
                  onChange={(e) => updateModelParams({ G: Number(e.target.value) })}
                  className="w-full accent-amber-500"
                />
              </div>
            </div>
          </GlassCard>

          {/* Real-time Insights snippet */}
          <GlassCard className="bg-slate-900/60 border-purple-500/20">
             <h3 className="text-sm font-semibold text-slate-400 mb-3 uppercase tracking-wider">Live Fit Analysis</h3>
             <div className="flex justify-between items-center">
               <span className="text-slate-300">Custom Model MSE:</span>
               <span className="font-mono text-lg text-purple-400">{formatScientific(insights.mseCustom)}</span>
             </div>
             <div className="flex justify-between items-center mt-2">
               <span className="text-slate-300">Vs Newtonian:</span>
               <span className={`font-medium ${insights.betterModel === 'Custom' ? 'text-emerald-400' : 'text-amber-400'}`}>
                 {insights.betterModel === 'Custom' ? 'Performing Better' : 'Performing Worse'}
               </span>
             </div>
          </GlassCard>
        </div>

      </div>
    </Layout>
  );
}
