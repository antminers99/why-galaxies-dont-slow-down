import React from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { useGalaxy } from '@/hooks/use-galaxy';
import { Link } from 'wouter';
import { ArrowRight, Database, TrendingUp, Sparkles, Activity, Zap } from 'lucide-react';
import { formatScientific } from '@/lib/utils';

export default function Dashboard() {
  const { datasets, activeDatasetIds, getInsights, modelParams, loadAllSamples } = useGalaxy();
  const insights = getInsights();

  const activeCount = datasets.filter(d => activeDatasetIds.includes(d.id)).length;

  return (
    <Layout>
      <header className="mb-10">
        <h1 className="text-4xl md:text-5xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-white to-slate-400">
          Galaxy Rotation Analysis
        </h1>
        <p className="mt-3 text-slate-400 max-w-2xl text-lg">
          Investigate the missing mass problem. Upload rotation curves, build mathematical models, and discover where standard physics breaks down.
        </p>
      </header>

      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-4 gap-5 mb-10">
        <GlassCard glow="cyan" className="flex flex-col">
          <div className="flex items-center gap-3 mb-4">
            <div className="p-3 bg-cyan-500/10 rounded-lg">
              <Database className="w-6 h-6 text-cyan-400" />
            </div>
            <h3 className="text-lg font-semibold">Datasets</h3>
          </div>
          <p className="text-4xl font-mono mb-1">{datasets.length}</p>
          <p className="text-slate-400 text-sm mb-4 flex-1">{activeCount} active rotation curves</p>
          <Link href="/upload" className="inline-flex items-center text-cyan-400 hover:text-cyan-300 font-medium transition-colors text-sm">
            Manage Data <ArrowRight className="w-4 h-4 ml-2" />
          </Link>
        </GlassCard>

        <GlassCard glow="amber" className="flex flex-col">
          <div className="flex items-center gap-3 mb-4">
            <div className="p-3 bg-amber-500/10 rounded-lg">
              <Sparkles className="w-6 h-6 text-amber-400" />
            </div>
            <h3 className="text-lg font-semibold">Anomalies</h3>
          </div>
          <p className="text-4xl font-mono text-amber-400 mb-1">{insights.anomalies.length}</p>
          <p className="text-slate-400 text-sm mb-4 flex-1">Points deviating {'>'}15% from model</p>
          <Link href="/analysis" className="inline-flex items-center text-amber-400 hover:text-amber-300 font-medium transition-colors text-sm">
            View Analysis <ArrowRight className="w-4 h-4 ml-2" />
          </Link>
        </GlassCard>

        <GlassCard glow="purple" className="flex flex-col">
          <div className="flex items-center gap-3 mb-4">
            <div className="p-3 bg-purple-500/10 rounded-lg">
              <TrendingUp className="w-6 h-6 text-purple-400" />
            </div>
            <h3 className="text-lg font-semibold">Model Status</h3>
          </div>
          <div className="mb-1">
            <span className="text-xs text-slate-400">Winner:</span>
            <p className="text-2xl font-bold text-purple-400">{insights.betterModel}</p>
          </div>
          <div className="flex gap-3 text-xs font-mono text-slate-300 mb-4 flex-1">
            <div>N: {formatScientific(insights.mseNewton)}</div>
            <div>C: {formatScientific(insights.mseCustom)}</div>
          </div>
          <Link href="/models" className="inline-flex items-center text-purple-400 hover:text-purple-300 font-medium transition-colors text-sm">
            Tune Model <ArrowRight className="w-4 h-4 ml-2" />
          </Link>
        </GlassCard>

        <GlassCard className="flex flex-col" glow={insights.generalizationScore >= 80 ? 'cyan' : 'none'}>
          <div className="flex items-center gap-3 mb-4">
            <div className="p-3 bg-emerald-500/10 rounded-lg">
              <Activity className="w-6 h-6 text-emerald-400" />
            </div>
            <h3 className="text-lg font-semibold">Generalization</h3>
          </div>
          <p className={`text-4xl font-mono mb-1 ${
            insights.generalizationScore >= 80 ? 'text-emerald-400' :
            insights.generalizationScore >= 50 ? 'text-amber-400' : 'text-red-400'
          }`}>
            {insights.generalizationScore.toFixed(0)}%
          </p>
          <p className="text-slate-400 text-sm mb-4 flex-1">
            Custom wins on {insights.perGalaxy.filter(g => g.winner === 'Custom').length}/{insights.perGalaxy.length} galaxies
          </p>
          <Link href="/models" className="inline-flex items-center text-emerald-400 hover:text-emerald-300 font-medium transition-colors text-sm">
            Compare <ArrowRight className="w-4 h-4 ml-2" />
          </Link>
        </GlassCard>
      </div>

      <div className="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-8">
        <div>
          <h2 className="text-xl font-semibold mb-4">Current Model</h2>
          <GlassCard>
            <div className="grid grid-cols-2 gap-4 font-mono text-sm mb-4">
              <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
                <div className="text-slate-400 text-xs mb-1">G (Gravity)</div>
                <div className="text-white">{formatScientific(modelParams.G)}</div>
              </div>
              <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
                <div className="text-slate-400 text-xs mb-1">M (Mass)</div>
                <div className="text-white">{formatScientific(modelParams.M)}</div>
              </div>
              <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
                <div className="text-slate-400 text-xs mb-1">k (Dark Matter)</div>
                <div className="text-white">{modelParams.k}</div>
              </div>
              <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
                <div className="text-slate-400 text-xs mb-1">a (Core Radius)</div>
                <div className="text-white">{modelParams.a.toFixed(1)}</div>
              </div>
            </div>
            <div className="p-3 bg-slate-900/50 rounded-xl border border-white/5">
              <div className="text-slate-400 text-xs mb-1 font-sans">Active Formula</div>
              <div className="text-lg text-purple-400 font-mono">{modelParams.formula}</div>
            </div>
          </GlassCard>
        </div>

        <div>
          <h2 className="text-xl font-semibold mb-4">Quick Actions</h2>
          <div className="space-y-3">
            <button
              onClick={loadAllSamples}
              className="w-full p-4 bg-slate-800/40 border border-white/5 rounded-xl hover:bg-slate-800/60 transition-colors text-left flex items-center gap-4"
            >
              <div className="p-3 bg-cyan-500/10 rounded-lg">
                <Database className="w-5 h-5 text-cyan-400" />
              </div>
              <div>
                <div className="font-medium">Load All 5 Sample Galaxies</div>
                <div className="text-xs text-slate-400 mt-0.5">M31, NGC 3198, Milky Way, NGC 6503, UGC 2885</div>
              </div>
            </button>
            <Link href="/models" className="block w-full p-4 bg-slate-800/40 border border-white/5 rounded-xl hover:bg-slate-800/60 transition-colors text-left">
              <div className="flex items-center gap-4">
                <div className="p-3 bg-purple-500/10 rounded-lg">
                  <Zap className="w-5 h-5 text-purple-400" />
                </div>
                <div>
                  <div className="font-medium">Auto-Optimize Parameters</div>
                  <div className="text-xs text-slate-400 mt-0.5">Let the optimizer find the best k, a, M values</div>
                </div>
              </div>
            </Link>
            <Link href="/analysis" className="block w-full p-4 bg-slate-800/40 border border-white/5 rounded-xl hover:bg-slate-800/60 transition-colors text-left">
              <div className="flex items-center gap-4">
                <div className="p-3 bg-amber-500/10 rounded-lg">
                  <Sparkles className="w-5 h-5 text-amber-400" />
                </div>
                <div>
                  <div className="font-medium">Discovery Mode</div>
                  <div className="text-xs text-slate-400 mt-0.5">Find where standard physics breaks down</div>
                </div>
              </div>
            </Link>
          </div>
        </div>
      </div>
    </Layout>
  );
}
