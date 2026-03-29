import React, { useCallback, useState } from 'react';
import { Layout } from '@/components/layout';
import { GlassCard } from '@/components/ui/glass-card';
import { useGalaxy } from '@/hooks/use-galaxy';
import { UploadCloud, X, Check, Search, Eye, Layers } from 'lucide-react';
import { useDropzone } from 'react-dropzone';
import { motion, AnimatePresence } from 'framer-motion';

const SAMPLE_INFO: Record<string, { desc: string; color: string }> = {
  "M31 (Andromeda)": { desc: "High mass spiral galaxy", color: "hover:bg-cyan-500/10 hover:text-cyan-400 hover:border-cyan-500/30" },
  "NGC 3198": { desc: "Extended flat rotation curve", color: "hover:bg-purple-500/10 hover:text-purple-400 hover:border-purple-500/30" },
  "Milky Way": { desc: "Our home galaxy", color: "hover:bg-green-500/10 hover:text-green-400 hover:border-green-500/30" },
  "NGC 6503": { desc: "Dwarf spiral with dark halo", color: "hover:bg-orange-500/10 hover:text-orange-400 hover:border-orange-500/30" },
  "UGC 2885": { desc: "Giant galaxy, flat at 300 km/s", color: "hover:bg-rose-500/10 hover:text-rose-400 hover:border-rose-500/30" },
};

export default function UploadPage() {
  const { 
    datasets, uploadDataset, loadSampleDataset, loadAllSamples,
    removeDataset, toggleDatasetActive, activeDatasetIds, sampleDatasetNames
  } = useGalaxy();
  const [error, setError] = useState<string | null>(null);
  const [previewId, setPreviewId] = useState<string | null>(null);

  const onDrop = useCallback(async (acceptedFiles: File[]) => {
    setError(null);
    try {
      for (const file of acceptedFiles) {
        await uploadDataset(file);
      }
    } catch (err: any) {
      setError(err.message || "Failed to parse file.");
    }
  }, [uploadDataset]);

  const { getRootProps, getInputProps, isDragActive } = useDropzone({ 
    onDrop,
    accept: {
      'text/csv': ['.csv'],
      'text/plain': ['.dat', '.txt']
    }
  });

  return (
    <Layout>
      <header className="mb-8">
        <h1 className="text-3xl font-bold">Data Management</h1>
        <p className="text-slate-400 mt-2">Upload and manage galaxy rotation curve datasets.</p>
      </header>

      <div className="grid grid-cols-1 lg:grid-cols-3 gap-8">
        
        <div className="lg:col-span-1 space-y-6">
          <GlassCard className="p-0 overflow-hidden">
            <div 
              {...getRootProps()} 
              className={`p-10 border-2 border-dashed flex flex-col items-center justify-center text-center cursor-pointer transition-all duration-300 ${
                isDragActive ? 'border-cyan-400 bg-cyan-500/10' : 'border-slate-600 bg-slate-800/20 hover:bg-slate-800/40 hover:border-cyan-500/50'
              }`}
            >
              <input {...getInputProps()} />
              <div className={`p-4 rounded-full mb-4 transition-colors ${isDragActive ? 'bg-cyan-500/20 text-cyan-400' : 'bg-slate-700/50 text-slate-400'}`}>
                <UploadCloud className="w-8 h-8" />
              </div>
              <h3 className="font-semibold text-lg mb-1">Upload Data</h3>
              <p className="text-sm text-slate-400">Drag & drop .csv or .dat files here</p>
              <p className="text-xs text-slate-500 mt-4">Needs columns: radius (r), velocity (v)</p>
            </div>
            
            {error && (
              <div className="p-4 bg-red-500/10 border-t border-red-500/20 text-red-400 text-sm">
                Error: {error}
              </div>
            )}
          </GlassCard>

          <GlassCard>
            <div className="flex items-center justify-between mb-4">
              <h3 className="font-semibold text-slate-200">Sample Datasets</h3>
              <button
                onClick={loadAllSamples}
                className="flex items-center gap-1.5 px-3 py-1.5 text-xs font-medium bg-gradient-to-r from-cyan-600 to-purple-600 hover:from-cyan-500 hover:to-purple-500 text-white rounded-lg transition-all"
              >
                <Layers className="w-3.5 h-3.5" /> Load All
              </button>
            </div>
            <div className="space-y-2">
              {sampleDatasetNames.map(name => {
                const info = SAMPLE_INFO[name] || { desc: "Sample galaxy", color: "hover:bg-slate-500/10 hover:text-slate-300" };
                const isLoaded = datasets.some(d => d.name === name);
                return (
                  <button 
                    key={name}
                    onClick={() => loadSampleDataset(name)}
                    className={`w-full text-left px-4 py-3 rounded-xl bg-slate-800/50 border border-transparent transition-all ${
                      isLoaded ? 'opacity-60' : info.color
                    }`}
                    disabled={isLoaded}
                  >
                    <div className="flex items-center justify-between">
                      <div className="font-medium text-sm">{name}</div>
                      {isLoaded && <span className="text-xs text-slate-500">Loaded</span>}
                    </div>
                    <div className="text-xs text-slate-500 mt-0.5">{info.desc}</div>
                  </button>
                );
              })}
            </div>
          </GlassCard>
        </div>

        <div className="lg:col-span-2 space-y-6">
          <h2 className="text-xl font-semibold mb-2">Loaded Datasets ({datasets.length})</h2>
          
          <AnimatePresence>
            {datasets.length === 0 && (
              <motion.div 
                initial={{ opacity: 0 }} animate={{ opacity: 1 }}
                className="p-12 text-center text-slate-500 border border-dashed border-slate-700 rounded-2xl"
              >
                <Search className="w-12 h-12 mx-auto mb-4 opacity-20" />
                <p>No datasets loaded. Upload or select a sample to begin.</p>
              </motion.div>
            )}

            {datasets.map((ds) => (
              <motion.div
                key={ds.id}
                initial={{ opacity: 0, x: -20 }}
                animate={{ opacity: 1, x: 0 }}
                exit={{ opacity: 0, scale: 0.95 }}
              >
                <div className="flex items-center justify-between p-4 bg-slate-800/40 border border-white/5 rounded-t-2xl backdrop-blur-sm group hover:bg-slate-800/60 transition-colors">
                  <div className="flex items-center gap-4">
                    <button 
                      onClick={() => toggleDatasetActive(ds.id)}
                      className={`w-6 h-6 rounded border flex items-center justify-center transition-colors ${
                        activeDatasetIds.includes(ds.id) 
                          ? 'bg-cyan-500 border-cyan-400 text-white' 
                          : 'border-slate-500 text-transparent hover:border-cyan-400'
                      }`}
                    >
                      <Check className="w-4 h-4" />
                    </button>
                    <div 
                      className="w-8 h-8 rounded-full shrink-0"
                      style={{ backgroundColor: ds.color || '#06b6d4' }}
                    />
                    <div>
                      <h4 className="font-medium text-slate-100">{ds.name}</h4>
                      <p className="text-xs font-mono text-slate-400">{ds.data.length} points • max r: {Math.max(...ds.data.map(d=>d.r)).toFixed(1)} kpc</p>
                    </div>
                  </div>
                  
                  <div className="flex items-center gap-2">
                    <button
                      onClick={() => setPreviewId(previewId === ds.id ? null : ds.id)}
                      className="p-2 text-slate-400 hover:text-cyan-400 hover:bg-cyan-400/10 rounded-lg transition-all"
                      title="Preview Data"
                    >
                      <Eye className="w-5 h-5" />
                    </button>
                    <button 
                      onClick={() => removeDataset(ds.id)}
                      className="p-2 text-slate-500 hover:text-red-400 hover:bg-red-400/10 rounded-lg transition-all"
                      title="Remove Dataset"
                    >
                      <X className="w-5 h-5" />
                    </button>
                  </div>
                </div>
                
                <AnimatePresence>
                  {previewId === ds.id && (
                    <motion.div
                      initial={{ height: 0, opacity: 0 }}
                      animate={{ height: 'auto', opacity: 1 }}
                      exit={{ height: 0, opacity: 0 }}
                      className="overflow-hidden bg-slate-900/60 border border-t-0 border-white/5 rounded-b-2xl"
                    >
                      <div className="max-h-[300px] overflow-auto">
                        <table className="w-full text-sm font-mono">
                          <thead className="sticky top-0 bg-slate-800/90 backdrop-blur-sm">
                            <tr>
                              <th className="px-4 py-2 text-left text-slate-400 font-medium">#</th>
                              <th className="px-4 py-2 text-left text-cyan-400 font-medium">Radius (kpc)</th>
                              <th className="px-4 py-2 text-left text-cyan-400 font-medium">Velocity (km/s)</th>
                            </tr>
                          </thead>
                          <tbody>
                            {ds.data.map((point, i) => (
                              <tr key={i} className="border-t border-white/5 hover:bg-white/5 transition-colors">
                                <td className="px-4 py-2 text-slate-500">{i + 1}</td>
                                <td className="px-4 py-2 text-slate-200">{point.r.toFixed(2)}</td>
                                <td className="px-4 py-2 text-slate-200">{point.v.toFixed(2)}</td>
                              </tr>
                            ))}
                          </tbody>
                        </table>
                      </div>
                    </motion.div>
                  )}
                </AnimatePresence>
              </motion.div>
            ))}
          </AnimatePresence>
        </div>
      </div>
    </Layout>
  );
}
