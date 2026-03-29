import React, { createContext, useContext, useState, useCallback, ReactNode, useEffect } from 'react';
import Papa from 'papaparse';
import * as math from 'mathjs';

export interface DataPoint {
  r: number; // Radius (kpc)
  v: number; // Velocity (km/s)
}

export interface Dataset {
  id: string;
  name: string;
  data: DataPoint[];
  color?: string;
}

export interface ModelParams {
  G: number;
  M: number;
  k: number;
  formula: string;
  [key: string]: string | number; // For dynamic custom params
}

export interface InsightStats {
  mseNewton: number;
  mseCustom: number;
  betterModel: 'Newtonian' | 'Custom' | 'Tie' | 'N/A';
  anomalies: DataPoint[];
}

interface GalaxyContextType {
  datasets: Dataset[];
  activeDatasetIds: string[];
  modelParams: ModelParams;
  showObserved: boolean;
  showNewtonian: boolean;
  showCustom: boolean;
  discoveryMode: boolean;
  
  // Actions
  uploadDataset: (file: File) => Promise<void>;
  loadSampleDataset: (name: string) => void;
  removeDataset: (id: string) => void;
  toggleDatasetActive: (id: string) => void;
  updateModelParams: (params: Partial<ModelParams>) => void;
  toggleLayer: (layer: 'observed' | 'newtonian' | 'custom' | 'discovery') => void;
  
  // Computed
  evaluateModel: (r: number, type: 'newtonian' | 'custom') => number | null;
  getInsights: () => InsightStats;
  generateChartData: () => any[];
}

const defaultParams: ModelParams = {
  G: 4.3009e-6, // Gravitational constant approx in kpc (km/s)^2 / M_sun
  M: 1e11,      // 100 billion solar masses
  k: 50,        // Custom dark matter halo parameter
  formula: "sqrt((G * M) / r + k * r)",
};

const SAMPLE_DATASETS: Record<string, Dataset> = {
  "M31 (Andromeda)": {
    id: "sample-m31",
    name: "M31 (Andromeda)",
    data: [
      { r: 2, v: 100 }, { r: 4, v: 210 }, { r: 6, v: 235 }, { r: 8, v: 245 },
      { r: 10, v: 250 }, { r: 15, v: 255 }, { r: 20, v: 250 }, { r: 25, v: 245 },
      { r: 30, v: 240 }, { r: 35, v: 238 }
    ],
    color: "hsl(189, 94%, 43%)"
  },
  "NGC 3198": {
    id: "sample-ngc3198",
    name: "NGC 3198",
    data: [
      { r: 1, v: 50 }, { r: 3, v: 120 }, { r: 5, v: 140 }, { r: 10, v: 150 },
      { r: 15, v: 155 }, { r: 20, v: 152 }, { r: 25, v: 150 }, { r: 30, v: 148 }
    ],
    color: "hsl(280, 80%, 60%)"
  },
  "Milky Way": {
    id: "sample-mw",
    name: "Milky Way",
    data: [
      { r: 1, v: 200 }, { r: 2, v: 215 }, { r: 3, v: 220 }, { r: 4, v: 225 },
      { r: 5, v: 225 }, { r: 6, v: 220 }, { r: 8, v: 220 }, { r: 10, v: 210 },
      { r: 12, v: 215 }, { r: 14, v: 220 }, { r: 16, v: 230 }, { r: 18, v: 235 },
      { r: 20, v: 240 }, { r: 25, v: 235 }, { r: 30, v: 230 }
    ],
    color: "hsl(142, 76%, 45%)"
  }
};

const GalaxyContext = createContext<GalaxyContextType | undefined>(undefined);

export const GalaxyProvider = ({ children }: { children: ReactNode }) => {
  const [datasets, setDatasets] = useState<Dataset[]>([SAMPLE_DATASETS["M31 (Andromeda)"]]);
  const [activeDatasetIds, setActiveDatasetIds] = useState<string[]>(["sample-m31"]);
  const [modelParams, setModelParams] = useState<ModelParams>(defaultParams);
  
  const [showObserved, setShowObserved] = useState(true);
  const [showNewtonian, setShowNewtonian] = useState(true);
  const [showCustom, setShowCustom] = useState(true);
  const [discoveryMode, setDiscoveryMode] = useState(false);

  const uploadDataset = async (file: File) => {
    return new Promise<void>((resolve, reject) => {
      Papa.parse(file, {
        header: true,
        dynamicTyping: true,
        skipEmptyLines: true,
        complete: (results) => {
          try {
            // Find columns that look like radius and velocity
            const headers = results.meta.fields || [];
            let rCol = headers.find(h => h.toLowerCase().includes('r') || h.toLowerCase() === 'rad' || h.toLowerCase() === 'radius');
            let vCol = headers.find(h => h.toLowerCase().includes('v') || h.toLowerCase() === 'vel' || h.toLowerCase() === 'velocity');
            
            // Fallback to first two columns if no obvious names
            if (!rCol || !vCol) {
              if (headers.length >= 2) {
                rCol = headers[0];
                vCol = headers[1];
              } else {
                throw new Error("Could not identify radius and velocity columns.");
              }
            }

            const data: DataPoint[] = results.data
              .map((row: any) => ({
                r: Number(row[rCol]),
                v: Number(row[vCol])
              }))
              .filter(p => !isNaN(p.r) && !isNaN(p.v) && p.r > 0)
              .sort((a, b) => a.r - b.r);

            if (data.length === 0) throw new Error("No valid numeric data found.");

            const newDataset: Dataset = {
              id: `ds-${Date.now()}`,
              name: file.name,
              data,
              color: `hsl(${Math.floor(Math.random() * 360)}, 80%, 60%)`
            };

            setDatasets(prev => [...prev, newDataset]);
            setActiveDatasetIds(prev => [...new Set([...prev, newDataset.id])]);
            resolve();
          } catch (e) {
            reject(e);
          }
        },
        error: (error) => reject(error)
      });
    });
  };

  const loadSampleDataset = (name: string) => {
    const ds = SAMPLE_DATASETS[name];
    if (ds && !datasets.find(d => d.id === ds.id)) {
      setDatasets(prev => [...prev, ds]);
      setActiveDatasetIds(prev => [...new Set([...prev, ds.id])]);
    } else if (ds) {
      // Already loaded, just make it active
      setActiveDatasetIds(prev => [...new Set([...prev, ds.id])]);
    }
  };

  const removeDataset = (id: string) => {
    setDatasets(prev => prev.filter(d => d.id !== id));
    setActiveDatasetIds(prev => prev.filter(a => a !== id));
  };

  const toggleDatasetActive = (id: string) => {
    setActiveDatasetIds(prev => 
      prev.includes(id) ? prev.filter(a => a !== id) : [...prev, id]
    );
  };

  const updateModelParams = (params: Partial<ModelParams>) => {
    setModelParams(prev => ({ ...prev, ...params }));
  };

  const toggleLayer = (layer: 'observed' | 'newtonian' | 'custom' | 'discovery') => {
    if (layer === 'observed') setShowObserved(p => !p);
    if (layer === 'newtonian') setShowNewtonian(p => !p);
    if (layer === 'custom') setShowCustom(p => !p);
    if (layer === 'discovery') setDiscoveryMode(p => !p);
  };

  const evaluateModel = useCallback((r: number, type: 'newtonian' | 'custom'): number | null => {
    if (r <= 0) return 0;
    try {
      if (type === 'newtonian') {
        return Math.sqrt((modelParams.G * modelParams.M) / r);
      } else {
        const scope = { r, ...modelParams };
        return math.evaluate(modelParams.formula, scope);
      }
    } catch (e) {
      return null; // Invalid formula
    }
  }, [modelParams]);

  const generateChartData = useCallback(() => {
    const activeData = datasets.filter(d => activeDatasetIds.includes(d.id));
    if (activeData.length === 0) return [];

    let maxR = 0;
    let maxV = 0;
    activeData.forEach(ds => {
      ds.data.forEach(p => { 
        if (p.r > maxR) maxR = p.r; 
        if (p.v > maxV) maxV = p.v;
      });
    });
    
    const yClamp = maxV * 1.5;
    
    const points = [];
    const steps = 100;
    const stepSize = Math.max(maxR / steps, 0.1);

    for (let r = 0.5; r <= maxR * 1.1; r += stepSize) {
      const p: any = { r: Number(r.toFixed(2)) };
      if (showNewtonian) {
        const vN = evaluateModel(r, 'newtonian');
        if (vN !== null && vN <= yClamp && vN >= 0) p.vNewtonian = Number(vN.toFixed(2));
      }
      if (showCustom) {
        const vC = evaluateModel(r, 'custom');
        if (vC !== null && vC <= yClamp && vC >= 0) p.vCustom = Number(vC.toFixed(2));
      }
      points.push(p);
    }

    // Merge observed data onto the closest generated points or just return separately
    // For Recharts ComposedChart, it's easier to have a continuous domain.
    // We will let Recharts handle the interpolation for lines, and scatter for points.
    
    // Actually, returning a unified array sorted by R is best for Recharts continuous XAxis
    let combinedData = [...points];
    
    if (showObserved || discoveryMode) {
      activeData.forEach(ds => {
        ds.data.forEach(obsPoint => {
          let dp: any = { r: obsPoint.r };
          dp[`vObs_${ds.id}`] = obsPoint.v;
          
          if (discoveryMode) {
            const vMod = evaluateModel(obsPoint.r, 'custom') || 0;
            const diff = Math.abs(obsPoint.v - vMod);
            // Threshold for anomaly (e.g. 15% deviation)
            if (diff > obsPoint.v * 0.15) {
              dp.anomaly = obsPoint.v;
            }
          }
          combinedData.push(dp);
        });
      });
    }

    return combinedData.sort((a, b) => a.r - b.r);
  }, [datasets, activeDatasetIds, showObserved, showNewtonian, showCustom, discoveryMode, evaluateModel]);

  const getInsights = useCallback((): InsightStats => {
    const activeData = datasets.filter(d => activeDatasetIds.includes(d.id));
    let mseNewt = 0;
    let mseCust = 0;
    let count = 0;
    const anomalies: DataPoint[] = [];

    activeData.forEach(ds => {
      ds.data.forEach(p => {
        const vNewt = evaluateModel(p.r, 'newtonian') || 0;
        const vCust = evaluateModel(p.r, 'custom') || 0;
        
        mseNewt += Math.pow(p.v - vNewt, 2);
        mseCust += Math.pow(p.v - vCust, 2);
        count++;

        if (Math.abs(p.v - vCust) > p.v * 0.15) {
          anomalies.push(p);
        }
      });
    });

    if (count > 0) {
      mseNewt /= count;
      mseCust /= count;
    }

    let better: InsightStats['betterModel'] = 'N/A';
    if (count > 0) {
      if (mseCust < mseNewt * 0.95) better = 'Custom';
      else if (mseNewt < mseCust * 0.95) better = 'Newtonian';
      else better = 'Tie';
    }

    return {
      mseNewton: mseNewt,
      mseCustom: mseCust,
      betterModel: better,
      anomalies
    };
  }, [datasets, activeDatasetIds, evaluateModel]);


  return (
    <GalaxyContext.Provider value={{
      datasets, activeDatasetIds, modelParams, 
      showObserved, showNewtonian, showCustom, discoveryMode,
      uploadDataset, loadSampleDataset, removeDataset, toggleDatasetActive,
      updateModelParams, toggleLayer,
      evaluateModel, getInsights, generateChartData
    }}>
      {children}
    </GalaxyContext.Provider>
  );
};

export const useGalaxy = () => {
  const context = useContext(GalaxyContext);
  if (context === undefined) {
    throw new Error('useGalaxy must be used within a GalaxyProvider');
  }
  return context;
};
