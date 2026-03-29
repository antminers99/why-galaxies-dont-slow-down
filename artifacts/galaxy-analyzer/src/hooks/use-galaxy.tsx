import React, { createContext, useContext, useState, useCallback, ReactNode } from 'react';
import Papa from 'papaparse';
import * as math from 'mathjs';

export interface DataPoint {
  r: number;
  v: number;
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
  a: number;
  formula: string;
}

export interface PerGalaxyMSE {
  galaxyName: string;
  galaxyId: string;
  mseNewton: number;
  mseCustom: number;
  winner: 'Newtonian' | 'Custom' | 'Tie';
  pointCount: number;
}

export interface ResidualPoint {
  r: number;
  galaxyName: string;
  observed: number;
  newtonian: number;
  custom: number;
  residualNewton: number;
  residualCustom: number;
  absResidualCustom: number;
  percentDeviation: number;
  isAnomaly: boolean;
}

export interface InsightStats {
  mseNewton: number;
  mseCustom: number;
  betterModel: 'Newtonian' | 'Custom' | 'Tie' | 'N/A';
  anomalies: DataPoint[];
  perGalaxy: PerGalaxyMSE[];
  residuals: ResidualPoint[];
  generalizationScore: number;
}

export interface FormulaPreset {
  id: string;
  name: string;
  formula: string;
  description: string;
  physicalMeaning: string;
  category: 'additive' | 'modified_gravity' | 'transition';
}

export const FORMULA_PRESETS: FormulaPreset[] = [
  {
    id: 'newtonian',
    name: 'Pure Newtonian',
    formula: 'sqrt((G * M) / r)',
    description: 'v = √(GM/r)',
    physicalMeaning: 'Standard Keplerian — predicts declining curve at large r. Fails for galaxies.',
    category: 'modified_gravity'
  },
  {
    id: 'additive_linear',
    name: 'Dark Halo (Linear)',
    formula: 'sqrt((G * M) / r + k * r)',
    description: 'v = √(GM/r + kr)',
    physicalMeaning: 'Adds a force that grows with distance — models a dark matter halo with linear density.',
    category: 'additive'
  },
  {
    id: 'additive_constant',
    name: 'Dark Halo (Flat)',
    formula: 'sqrt((G * M) / r + k)',
    description: 'v = √(GM/r + k)',
    physicalMeaning: 'Adds a constant velocity floor — models isothermal dark matter halo.',
    category: 'additive'
  },
  {
    id: 'modified_gravity',
    name: 'Modified Gravity',
    formula: 'sqrt((G * M) / (r + a))',
    description: 'v = √(GM/(r+a))',
    physicalMeaning: 'Gravity "softens" at a core radius a — prevents singularity and flattens inner curve.',
    category: 'modified_gravity'
  },
  {
    id: 'modified_gravity_halo',
    name: 'Modified Gravity + Halo',
    formula: 'sqrt((G * M) / (r + a) + k * r)',
    description: 'v = √(GM/(r+a) + kr)',
    physicalMeaning: 'Combines softened gravity with dark matter halo — two free parameters for fitting.',
    category: 'modified_gravity'
  },
  {
    id: 'transition',
    name: 'Transition Model',
    formula: 'sqrt((G * M) / r) * (1 + k * r / 100)',
    description: 'v = √(GM/r) × (1 + kr/100)',
    physicalMeaning: 'Newtonian at small r, gradually boosted at large r — models transition to dark matter dominated regime.',
    category: 'transition'
  },
  {
    id: 'mond_like',
    name: 'MOND-inspired',
    formula: 'sqrt(sqrt((G * M) / r) * a + (G * M) / r)',
    description: 'v = √(√(GM/r)·a + GM/r)',
    physicalMeaning: 'Inspired by Modified Newtonian Dynamics — gravity transitions at a critical acceleration scale.',
    category: 'modified_gravity'
  },
  {
    id: 'log_halo',
    name: 'Logarithmic Halo',
    formula: 'sqrt((G * M) / r + k * log(1 + r / a))',
    description: 'v = √(GM/r + k·ln(1+r/a))',
    physicalMeaning: 'NFW-inspired logarithmic dark matter halo profile — widely used in astrophysics.',
    category: 'additive'
  }
];

interface GalaxyContextType {
  datasets: Dataset[];
  activeDatasetIds: string[];
  modelParams: ModelParams;
  showObserved: boolean;
  showNewtonian: boolean;
  showCustom: boolean;
  discoveryMode: boolean;
  isOptimizing: boolean;
  optimizationLog: string[];

  uploadDataset: (file: File) => Promise<void>;
  loadSampleDataset: (name: string) => void;
  loadAllSamples: () => void;
  removeDataset: (id: string) => void;
  toggleDatasetActive: (id: string) => void;
  updateModelParams: (params: Partial<ModelParams>) => void;
  toggleLayer: (layer: 'observed' | 'newtonian' | 'custom' | 'discovery') => void;
  applyPreset: (preset: FormulaPreset) => void;
  autoOptimize: () => void;

  evaluateModel: (r: number, type: 'newtonian' | 'custom') => number | null;
  evaluateFormulaWithParams: (formula: string, params: ModelParams, r: number) => number | null;
  getInsights: () => InsightStats;
  generateChartData: () => any[];
  generateResidualData: () => ResidualPoint[];
  sampleDatasetNames: string[];
}

const defaultParams: ModelParams = {
  G: 4.3009e-6,
  M: 1e11,
  k: 50,
  a: 5,
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
  },
  "NGC 6503": {
    id: "sample-ngc6503",
    name: "NGC 6503",
    data: [
      { r: 0.5, v: 40 }, { r: 1, v: 70 }, { r: 2, v: 95 }, { r: 3, v: 108 },
      { r: 4, v: 115 }, { r: 5, v: 118 }, { r: 7, v: 120 }, { r: 9, v: 121 },
      { r: 12, v: 120 }, { r: 15, v: 119 }, { r: 18, v: 118 }, { r: 20, v: 117 }
    ],
    color: "hsl(32, 95%, 54%)"
  },
  "UGC 2885": {
    id: "sample-ugc2885",
    name: "UGC 2885",
    data: [
      { r: 5, v: 200 }, { r: 10, v: 280 }, { r: 15, v: 295 }, { r: 20, v: 300 },
      { r: 25, v: 298 }, { r: 30, v: 300 }, { r: 40, v: 305 }, { r: 50, v: 300 },
      { r: 60, v: 298 }, { r: 70, v: 295 }, { r: 80, v: 300 }
    ],
    color: "hsl(350, 80%, 55%)"
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
  const [isOptimizing, setIsOptimizing] = useState(false);
  const [optimizationLog, setOptimizationLog] = useState<string[]>([]);

  const uploadDataset = async (file: File) => {
    return new Promise<void>((resolve, reject) => {
      Papa.parse(file, {
        header: true,
        dynamicTyping: true,
        skipEmptyLines: true,
        complete: (results) => {
          try {
            const headers = results.meta.fields || [];
            let rCol = headers.find(h => /^r$|radius|rad/i.test(h));
            let vCol = headers.find(h => /^v$|velocity|vel/i.test(h));
            
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
              name: file.name.replace(/\.(csv|dat|txt)$/i, ''),
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
      setActiveDatasetIds(prev => [...new Set([...prev, ds.id])]);
    }
  };

  const loadAllSamples = () => {
    const allSamples = Object.values(SAMPLE_DATASETS);
    setDatasets(prev => {
      const existing = new Set(prev.map(d => d.id));
      const toAdd = allSamples.filter(s => !existing.has(s.id));
      return [...prev, ...toAdd];
    });
    setActiveDatasetIds(prev => [...new Set([...prev, ...allSamples.map(s => s.id)])]);
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

  const applyPreset = (preset: FormulaPreset) => {
    setModelParams(prev => ({ ...prev, formula: preset.formula }));
  };

  const evaluateFormulaWithParams = useCallback((formula: string, params: ModelParams, r: number): number | null => {
    if (r <= 0) return 0;
    try {
      const scope = { r, G: params.G, M: params.M, k: params.k, a: params.a };
      const result = math.evaluate(formula, scope);
      if (typeof result !== 'number' || isNaN(result) || !isFinite(result)) return null;
      return result;
    } catch {
      return null;
    }
  }, []);

  const evaluateModel = useCallback((r: number, type: 'newtonian' | 'custom'): number | null => {
    if (r <= 0) return 0;
    try {
      if (type === 'newtonian') {
        const val = Math.sqrt((modelParams.G * modelParams.M) / r);
        return isFinite(val) ? val : null;
      } else {
        return evaluateFormulaWithParams(modelParams.formula, modelParams, r);
      }
    } catch {
      return null;
    }
  }, [modelParams, evaluateFormulaWithParams]);

  const computeMSEForParams = useCallback((formula: string, params: ModelParams): number => {
    const activeData = datasets.filter(d => activeDatasetIds.includes(d.id));
    let sum = 0;
    let count = 0;
    activeData.forEach(ds => {
      ds.data.forEach(p => {
        const predicted = evaluateFormulaWithParams(formula, params, p.r);
        if (predicted !== null) {
          sum += Math.pow(p.v - predicted, 2);
          count++;
        }
      });
    });
    return count > 0 ? sum / count : Infinity;
  }, [datasets, activeDatasetIds, evaluateFormulaWithParams]);

  const autoOptimize = useCallback(() => {
    setIsOptimizing(true);
    setOptimizationLog([]);
    const log: string[] = [];
    
    log.push("Starting parameter optimization...");
    log.push(`Formula: ${modelParams.formula}`);
    log.push(`Active galaxies: ${activeDatasetIds.length}`);
    
    let bestK = modelParams.k;
    let bestA = modelParams.a;
    let bestM = modelParams.M;
    let bestMSE = computeMSEForParams(modelParams.formula, modelParams);
    
    log.push(`Initial MSE: ${bestMSE.toFixed(2)}`);
    log.push("");
    log.push("Phase 1: Coarse grid search...");

    const kValues = Array.from({ length: 40 }, (_, i) => i * 5);
    const aValues = [0.1, 0.5, 1, 2, 3, 5, 8, 10, 15, 20];
    const mValues = [1e10, 5e10, 1e11, 2e11, 5e11, 1e12];
    
    let iterations = 0;
    
    for (const testM of mValues) {
      for (const testK of kValues) {
        for (const testA of aValues) {
          const testParams = { ...modelParams, k: testK, a: testA, M: testM };
          const mse = computeMSEForParams(modelParams.formula, testParams);
          iterations++;
          if (mse < bestMSE) {
            bestMSE = mse;
            bestK = testK;
            bestA = testA;
            bestM = testM;
          }
        }
      }
    }
    
    log.push(`  Tested ${iterations} combinations`);
    log.push(`  Best coarse: k=${bestK}, a=${bestA.toFixed(1)}, M=${bestM.toExponential(1)}`);
    log.push(`  MSE: ${bestMSE.toFixed(2)}`);
    log.push("");
    log.push("Phase 2: Fine-tuning around best values...");
    
    const fineK = Array.from({ length: 21 }, (_, i) => Math.max(0, bestK - 10 + i));
    const fineA = Array.from({ length: 21 }, (_, i) => Math.max(0.01, bestA - 2 + i * 0.2));
    const fineM = Array.from({ length: 11 }, (_, i) => bestM * (0.5 + i * 0.1));
    
    iterations = 0;
    for (const testM of fineM) {
      for (const testK of fineK) {
        for (const testA of fineA) {
          const testParams = { ...modelParams, k: testK, a: testA, M: testM };
          const mse = computeMSEForParams(modelParams.formula, testParams);
          iterations++;
          if (mse < bestMSE) {
            bestMSE = mse;
            bestK = testK;
            bestA = testA;
            bestM = testM;
          }
        }
      }
    }
    
    log.push(`  Tested ${iterations} fine combinations`);
    log.push(`  Optimized: k=${bestK}, a=${bestA.toFixed(2)}, M=${bestM.toExponential(2)}`);
    log.push(`  Final MSE: ${bestMSE.toFixed(2)}`);
    log.push("");
    
    const newtonMSE = computeMSEForParams('sqrt((G * M) / r)', { ...modelParams, M: bestM });
    const improvement = ((newtonMSE - bestMSE) / newtonMSE * 100);
    
    log.push(`Newtonian MSE: ${newtonMSE.toFixed(2)}`);
    log.push(`Custom MSE: ${bestMSE.toFixed(2)}`);
    log.push(`Improvement: ${improvement.toFixed(1)}%`);
    log.push("");
    
    if (improvement > 50) {
      log.push("RESULT: Significant improvement over Newtonian model!");
    } else if (improvement > 20) {
      log.push("RESULT: Moderate improvement. Model shows promise.");
    } else if (improvement > 0) {
      log.push("RESULT: Slight improvement. Consider trying different formulas.");
    } else {
      log.push("RESULT: No improvement over Newtonian. Try a different formula structure.");
    }
    
    setModelParams(prev => ({ ...prev, k: bestK, a: bestA, M: bestM }));
    setOptimizationLog(log);
    setIsOptimizing(false);
  }, [modelParams, activeDatasetIds, computeMSEForParams]);

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
    const steps = 120;
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

    let combinedData = [...points];
    
    if (showObserved || discoveryMode) {
      activeData.forEach(ds => {
        ds.data.forEach(obsPoint => {
          let dp: any = { r: obsPoint.r };
          dp[`vObs_${ds.id}`] = obsPoint.v;
          
          if (discoveryMode) {
            const vMod = evaluateModel(obsPoint.r, 'custom') || 0;
            const diff = Math.abs(obsPoint.v - vMod);
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

  const generateResidualData = useCallback((): ResidualPoint[] => {
    const activeData = datasets.filter(d => activeDatasetIds.includes(d.id));
    const residuals: ResidualPoint[] = [];

    activeData.forEach(ds => {
      ds.data.forEach(p => {
        const vNewt = evaluateModel(p.r, 'newtonian') || 0;
        const vCust = evaluateModel(p.r, 'custom') || 0;
        const residualN = p.v - vNewt;
        const residualC = p.v - vCust;
        const percentDev = vCust !== 0 ? Math.abs(residualC / p.v * 100) : 0;
        
        residuals.push({
          r: p.r,
          galaxyName: ds.name,
          observed: p.v,
          newtonian: vNewt,
          custom: vCust,
          residualNewton: Number(residualN.toFixed(2)),
          residualCustom: Number(residualC.toFixed(2)),
          absResidualCustom: Number(Math.abs(residualC).toFixed(2)),
          percentDeviation: Number(percentDev.toFixed(1)),
          isAnomaly: percentDev > 15
        });
      });
    });

    return residuals.sort((a, b) => a.r - b.r);
  }, [datasets, activeDatasetIds, evaluateModel]);

  const getInsights = useCallback((): InsightStats => {
    const activeData = datasets.filter(d => activeDatasetIds.includes(d.id));
    let mseNewt = 0;
    let mseCust = 0;
    let count = 0;
    const anomalies: DataPoint[] = [];
    const perGalaxy: PerGalaxyMSE[] = [];

    activeData.forEach(ds => {
      let gMseN = 0, gMseC = 0, gCount = 0;
      
      ds.data.forEach(p => {
        const vNewt = evaluateModel(p.r, 'newtonian') || 0;
        const vCust = evaluateModel(p.r, 'custom') || 0;
        
        const errN = Math.pow(p.v - vNewt, 2);
        const errC = Math.pow(p.v - vCust, 2);
        
        mseNewt += errN;
        mseCust += errC;
        gMseN += errN;
        gMseC += errC;
        count++;
        gCount++;

        if (Math.abs(p.v - vCust) > p.v * 0.15) {
          anomalies.push(p);
        }
      });

      if (gCount > 0) {
        gMseN /= gCount;
        gMseC /= gCount;
        perGalaxy.push({
          galaxyName: ds.name,
          galaxyId: ds.id,
          mseNewton: gMseN,
          mseCustom: gMseC,
          winner: gMseC < gMseN * 0.95 ? 'Custom' : gMseN < gMseC * 0.95 ? 'Newtonian' : 'Tie',
          pointCount: gCount
        });
      }
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

    const customWins = perGalaxy.filter(g => g.winner === 'Custom').length;
    const generalizationScore = perGalaxy.length > 0 ? (customWins / perGalaxy.length) * 100 : 0;

    return {
      mseNewton: mseNewt,
      mseCustom: mseCust,
      betterModel: better,
      anomalies,
      perGalaxy,
      residuals: generateResidualData(),
      generalizationScore
    };
  }, [datasets, activeDatasetIds, evaluateModel, generateResidualData]);


  return (
    <GalaxyContext.Provider value={{
      datasets, activeDatasetIds, modelParams, 
      showObserved, showNewtonian, showCustom, discoveryMode,
      isOptimizing, optimizationLog,
      uploadDataset, loadSampleDataset, loadAllSamples, removeDataset, toggleDatasetActive,
      updateModelParams, toggleLayer, applyPreset, autoOptimize,
      evaluateModel, evaluateFormulaWithParams, getInsights, generateChartData, generateResidualData,
      sampleDatasetNames: Object.keys(SAMPLE_DATASETS)
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
