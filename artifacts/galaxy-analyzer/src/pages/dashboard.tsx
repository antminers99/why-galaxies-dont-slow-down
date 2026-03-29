import React from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { useGalaxy } from '@/hooks/use-galaxy';
import { Link } from 'wouter';
import { ArrowRight, Database, TrendingUp, Sparkles } from 'lucide-react';
import { formatScientific } from '@/lib/utils';
import { motion } from 'framer-motion';

export default function Dashboard() {
  const { datasets, getInsights, modelParams } = useGalaxy();
  const insights = getInsights();

  return (
    <Layout>
      <header className="mb-10">
        <h1 className="text-4xl md:text-5xl font-bold bg-clip-text text-transparent bg-gradient-to-r from-white to-slate-400">
          Galaxy Rotation Analysis
        </h1>
        <p className="mt-3 text-slate-400 max-w-2xl text-lg">
          Investigate the missing mass problem. Upload rotation curves, build mathematical models, and visualize dark matter halos.
        </p>
      </header>

      <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-6 mb-10">
        <GlassCard glow="cyan" className="flex flex-col">
          <div className="flex items-center gap-3 mb-4">
            <div className="p-3 bg-cyan-500/10 rounded-lg">
              <Database className="w-6 h-6 text-cyan-400" />
            </div>
            <h3 className="text-xl font-semibold">Datasets</h3>
          </div>
          <p className="text-4xl font-mono mb-2">{datasets.length}</p>
          <p className="text-slate-400 mb-6 flex-1">Active rotation curves loaded</p>
          <Link href="/upload" className="inline-flex items-center text-cyan-400 hover:text-cyan-300 font-medium transition-colors">
            Manage Data <ArrowRight className="w-4 h-4 ml-2" />
          </Link>
        </GlassCard>

        <GlassCard glow="amber" className="flex flex-col">
          <div className="flex items-center gap-3 mb-4">
            <div className="p-3 bg-amber-500/10 rounded-lg">
              <Sparkles className="w-6 h-6 text-amber-400" />
            </div>
            <h3 className="text-xl font-semibold">Anomalies</h3>
          </div>
          <p className="text-4xl font-mono text-amber-400 mb-2">{insights.anomalies.length}</p>
          <p className="text-slate-400 mb-6 flex-1">Data points deviating {'>'}15% from model</p>
          <Link href="/analysis" className="inline-flex items-center text-amber-400 hover:text-amber-300 font-medium transition-colors">
            View Analysis <ArrowRight className="w-4 h-4 ml-2" />
          </Link>
        </GlassCard>

        <GlassCard glow="purple" className="flex flex-col">
          <div className="flex items-center gap-3 mb-4">
            <div className="p-3 bg-purple-500/10 rounded-lg">
              <TrendingUp className="w-6 h-6 text-purple-400" />
            </div>
            <h3 className="text-xl font-semibold">Model Status</h3>
          </div>
          <div className="mb-2">
            <span className="text-sm text-slate-400">Winning Model:</span>
            <p className="text-2xl font-bold text-purple-400">{insights.betterModel}</p>
          </div>
          <div className="flex gap-4 text-sm font-mono text-slate-300 mb-6 flex-1">
            <div>N_MSE: {formatScientific(insights.mseNewton)}</div>
            <div>C_MSE: {formatScientific(insights.mseCustom)}</div>
          </div>
          <Link href="/models" className="inline-flex items-center text-purple-400 hover:text-purple-300 font-medium transition-colors">
            Tune Parameters <ArrowRight className="w-4 h-4 ml-2" />
          </Link>
        </GlassCard>
      </div>

      <h2 className="text-2xl font-semibold mb-6">Current Model Parameters</h2>
      <GlassCard>
        <div className="grid grid-cols-1 sm:grid-cols-3 gap-6 font-mono text-sm">
          <div className="p-4 bg-slate-900/50 rounded-xl border border-white/5">
            <div className="text-slate-400 mb-1">G (Gravity Const)</div>
            <div className="text-lg text-white">{formatScientific(modelParams.G)}</div>
          </div>
          <div className="p-4 bg-slate-900/50 rounded-xl border border-white/5">
            <div className="text-slate-400 mb-1">M (Mass Est.)</div>
            <div className="text-lg text-white">{formatScientific(modelParams.M)}</div>
          </div>
          <div className="p-4 bg-slate-900/50 rounded-xl border border-white/5">
            <div className="text-slate-400 mb-1">k (Dark Matter)</div>
            <div className="text-lg text-white">{modelParams.k}</div>
          </div>
        </div>
        <div className="mt-6 p-4 bg-slate-900/50 rounded-xl border border-white/5">
          <div className="text-slate-400 mb-1 font-sans">Active Custom Formula</div>
          <div className="text-xl text-purple-400 font-mono tracking-wider">{modelParams.formula}</div>
        </div>
      </GlassCard>
    </Layout>
  );
}
