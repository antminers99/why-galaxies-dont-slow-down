import React, { useRef, useMemo } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { useGalaxy } from '@/hooks/use-galaxy';
import { Switch } from '@/components/ui/switch';
import { Download, AlertTriangle, TrendingUp, BarChart3 } from 'lucide-react';
import { formatScientific } from '@/lib/utils';
import { 
  ComposedChart, Scatter, Line, XAxis, YAxis, CartesianGrid, 
  Tooltip, Legend, ResponsiveContainer
} from 'recharts';
import * as htmlToImage from 'html-to-image';

export default function AnalysisPage() {
  const { 
    datasets, activeDatasetIds, 
    showObserved, showNewtonian, showCustom, discoveryMode,
    toggleLayer, generateChartData, getInsights
  } = useGalaxy();

  const chartRef = useRef<HTMLDivElement>(null);
  const data = useMemo(() => generateChartData(), [generateChartData]);
  const insights = useMemo(() => getInsights(), [getInsights]);

  const activeDatasets = datasets.filter(d => activeDatasetIds.includes(d.id));

  const yDomain = useMemo(() => {
    if (activeDatasets.length === 0) return [0, 300];
    let maxV = 0;
    activeDatasets.forEach(ds => {
      ds.data.forEach(p => { if (p.v > maxV) maxV = p.v; });
    });
    return [0, Math.ceil(maxV * 1.4 / 50) * 50];
  }, [activeDatasets]);

  const exportGraph = async () => {
    if (chartRef.current) {
      try {
        const dataUrl = await htmlToImage.toPng(chartRef.current, { backgroundColor: '#0f172a' });
        const link = document.createElement('a');
        link.download = 'rotation-curve.png';
        link.href = dataUrl;
        link.click();
      } catch (err) {
        console.error('Failed to export graph', err);
      }
    }
  };

  const fitQuality = (mse: number): { label: string; color: string } => {
    if (mse < 500) return { label: 'Excellent', color: 'text-emerald-400' };
    if (mse < 2000) return { label: 'Good', color: 'text-cyan-400' };
    if (mse < 10000) return { label: 'Moderate', color: 'text-amber-400' };
    return { label: 'Poor', color: 'text-red-400' };
  };

  const customFit = fitQuality(insights.mseCustom);
  const newtonFit = fitQuality(insights.mseNewton);

  return (
    <Layout>
      <header className="flex justify-between items-end mb-6">
        <div>
          <h1 className="text-3xl font-bold">Rotation Curve Analysis</h1>
          <p className="text-slate-400 mt-1">Visualize and compare predictive models against observed data.</p>
        </div>
        <button 
          onClick={exportGraph}
          className="flex items-center gap-2 px-4 py-2 bg-slate-800 hover:bg-slate-700 text-white rounded-xl border border-white/10 transition-colors"
        >
          <Download className="w-4 h-4" /> Export PNG
        </button>
      </header>

      <div className="grid grid-cols-1 xl:grid-cols-4 gap-6">
        <div className="xl:col-span-1 space-y-6 flex flex-col">
          <GlassCard className="flex-1">
            <h3 className="font-semibold mb-6 border-b border-white/10 pb-4">Visualization Layers</h3>
            
            <div className="space-y-6">
              <div className="flex items-center justify-between">
                <div>
                  <div className="flex items-center gap-2">
                    <div className="w-3 h-3 rounded-full bg-cyan-500" />
                    <span className="font-medium">Observed Data</span>
                  </div>
                  <p className="text-xs text-slate-500 mt-1">Actual measured velocities</p>
                </div>
                <Switch checked={showObserved} onCheckedChange={() => toggleLayer('observed')} accent="cyan" />
              </div>

              <div className="flex items-center justify-between">
                <div>
                  <div className="flex items-center gap-2">
                    <div className="w-3 h-1 bg-orange-500" />
                    <span className="font-medium">Newtonian Model</span>
                  </div>
                  <p className="text-xs text-slate-500 mt-1">Expected without dark matter</p>
                </div>
                <Switch checked={showNewtonian} onCheckedChange={() => toggleLayer('newtonian')} accent="amber" />
              </div>

              <div className="flex items-center justify-between">
                <div>
                  <div className="flex items-center gap-2">
                    <div className="w-3 h-1 bg-purple-500" />
                    <span className="font-medium">Custom Model</span>
                  </div>
                  <p className="text-xs text-slate-500 mt-1">User defined formula</p>
                </div>
                <Switch checked={showCustom} onCheckedChange={() => toggleLayer('custom')} accent="purple" />
              </div>

              <div className="pt-6 border-t border-white/10">
                <div className="flex items-center justify-between">
                  <div>
                    <div className="flex items-center gap-2 text-amber-400">
                      <AlertTriangle className="w-4 h-4" />
                      <span className="font-medium">Discovery Mode</span>
                    </div>
                    <p className="text-xs text-slate-500 mt-1">Highlight large deviations</p>
                  </div>
                  <Switch checked={discoveryMode} onCheckedChange={() => toggleLayer('discovery')} accent="amber" />
                </div>
              </div>
            </div>
          </GlassCard>
        </div>

        <div className="xl:col-span-3 space-y-6">
          <GlassCard className="p-2 relative flex flex-col" glow={discoveryMode ? 'amber' : 'none'}>
            {activeDatasets.length === 0 ? (
              <div className="h-[500px] flex items-center justify-center text-slate-500">
                Please activate a dataset in the Datasets tab to view analysis.
              </div>
            ) : (
              <div className="w-full h-[500px] p-4" ref={chartRef}>
                <ResponsiveContainer width="100%" height="100%">
                  <ComposedChart data={data} margin={{ top: 20, right: 30, left: 20, bottom: 20 }}>
                    <CartesianGrid strokeDasharray="3 3" stroke="#334155" opacity={0.5} />
                    <XAxis 
                      dataKey="r" 
                      type="number" 
                      domain={['dataMin', 'dataMax']} 
                      name="Radius" 
                      unit=" kpc" 
                      stroke="#94a3b8"
                      tick={{fill: '#94a3b8', fontSize: 12}}
                      label={{ value: 'Radius (kpc)', position: 'insideBottom', offset: -10, fill: '#94a3b8' }}
                    />
                    <YAxis 
                      name="Velocity" 
                      unit=" km/s" 
                      stroke="#94a3b8"
                      tick={{fill: '#94a3b8', fontSize: 12}}
                      domain={yDomain}
                      label={{ value: 'Velocity (km/s)', angle: -90, position: 'insideLeft', offset: 10, fill: '#94a3b8' }}
                    />
                    <Tooltip 
                      contentStyle={{ backgroundColor: '#1e293b', borderColor: '#334155', color: '#f8fafc', borderRadius: '12px' }}
                      itemStyle={{ color: '#f8fafc' }}
                      labelFormatter={(val) => `Radius: ${Number(val).toFixed(2)} kpc`}
                    />
                    <Legend wrapperStyle={{ paddingTop: '20px' }}/>

                    {showNewtonian && (
                      <Line 
                        type="monotone" 
                        dataKey="vNewtonian" 
                        name="Newtonian Model" 
                        stroke="#f97316" 
                        strokeWidth={2} 
                        strokeDasharray="6 4"
                        dot={false}
                        isAnimationActive={false}
                        connectNulls={false}
                      />
                    )}
                    {showCustom && (
                      <Line 
                        type="monotone" 
                        dataKey="vCustom" 
                        name="Custom Model" 
                        stroke="#a855f7" 
                        strokeWidth={3} 
                        dot={false}
                        isAnimationActive={false}
                        connectNulls={false}
                      />
                    )}

                    {(showObserved || discoveryMode) && activeDatasets.map(ds => (
                      <Scatter 
                        key={ds.id} 
                        name={ds.name} 
                        dataKey={`vObs_${ds.id}`} 
                        fill={ds.color || "#06b6d4"} 
                        opacity={showObserved ? 0.9 : 0.2}
                      />
                    ))}

                    {discoveryMode && (
                      <Scatter 
                        name="Anomalies (>15% dev)" 
                        dataKey="anomaly" 
                        fill="#f59e0b" 
                        shape="cross"
                        isAnimationActive={true}
                      />
                    )}
                  </ComposedChart>
                </ResponsiveContainer>
              </div>
            )}
          </GlassCard>

          <div className="grid grid-cols-1 md:grid-cols-3 gap-4">
            <GlassCard className="p-4">
              <div className="flex items-center gap-2 mb-3 text-sm font-semibold text-slate-400 uppercase tracking-wider">
                <BarChart3 className="w-4 h-4" />
                Newtonian MSE
              </div>
              <div className="text-2xl font-mono text-orange-400 mb-1">{formatScientific(insights.mseNewton)}</div>
              <div className={`text-sm font-medium ${newtonFit.color}`}>{newtonFit.label} Fit</div>
            </GlassCard>

            <GlassCard className="p-4">
              <div className="flex items-center gap-2 mb-3 text-sm font-semibold text-slate-400 uppercase tracking-wider">
                <BarChart3 className="w-4 h-4" />
                Custom MSE
              </div>
              <div className="text-2xl font-mono text-purple-400 mb-1">{formatScientific(insights.mseCustom)}</div>
              <div className={`text-sm font-medium ${customFit.color}`}>{customFit.label} Fit</div>
            </GlassCard>

            <GlassCard className="p-4" glow={insights.betterModel === 'Custom' ? 'purple' : 'none'}>
              <div className="flex items-center gap-2 mb-3 text-sm font-semibold text-slate-400 uppercase tracking-wider">
                <TrendingUp className="w-4 h-4" />
                Winner
              </div>
              <div className={`text-2xl font-bold mb-1 ${
                insights.betterModel === 'Custom' ? 'text-purple-400' : 
                insights.betterModel === 'Newtonian' ? 'text-orange-400' : 'text-slate-300'
              }`}>
                {insights.betterModel}
              </div>
              <div className="text-sm text-slate-400">
                {insights.anomalies.length} anomalous point{insights.anomalies.length !== 1 ? 's' : ''} detected
              </div>
            </GlassCard>
          </div>
        </div>
      </div>
    </Layout>
  );
}
