import React, { createContext, useContext, useState, useCallback, ReactNode } from 'react';
import Papa from 'papaparse';
import * as math from 'mathjs';
import { SPARC_DATASETS } from '@/data/sparc-datasets';

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

export interface BenchmarkGalaxyResult {
  galaxyName: string;
  galaxyId: string;
  mseNewton: number;
  mseCustom: number;
  improvementPct: number;
  bestK: number;
  bestA: number;
  bestM: number;
  mseInner: number;
  mseOuter: number;
  innerImprovement: number;
  outerImprovement: number;
  pointCount: number;
}

export interface FormulaBenchmark {
  formulaId: string;
  formulaName: string;
  formula: string;
  avgMSE: number;
  avgImprovement: number;
  winsCount: number;
  galaxyResults: BenchmarkGalaxyResult[];
  avgK: number;
  kStdDev: number;
  avgA: number;
  avgM: number;
}

export interface BenchmarkResult {
  formulas: FormulaBenchmark[];
  bestFormula: string;
  kConsistency: string;
  timestamp: string;
}

export const FORMULA_PRESETS: FormulaPreset[] = [
  {
    id: 'newtonian',
    name: 'Pure Newtonian',
    formula: 'sqrt((G * M) / r)',
    description: 'v = sqrt(GM/r)',
    physicalMeaning: 'Standard Keplerian — predicts declining curve at large r. Fails for galaxies.',
    category: 'modified_gravity'
  },
  {
    id: 'additive_linear',
    name: 'Dark Halo (Linear)',
    formula: 'sqrt((G * M) / r + k * r)',
    description: 'v = sqrt(GM/r + kr)',
    physicalMeaning: 'Adds a force that grows with distance — models a dark matter halo with linear density.',
    category: 'additive'
  },
  {
    id: 'additive_constant',
    name: 'Dark Halo (Flat)',
    formula: 'sqrt((G * M) / r + k)',
    description: 'v = sqrt(GM/r + k)',
    physicalMeaning: 'Adds a constant velocity floor — models isothermal dark matter halo.',
    category: 'additive'
  },
  {
    id: 'modified_gravity',
    name: 'Modified Gravity',
    formula: 'sqrt((G * M) / (r + a))',
    description: 'v = sqrt(GM/(r+a))',
    physicalMeaning: 'Gravity "softens" at a core radius a — prevents singularity and flattens inner curve.',
    category: 'modified_gravity'
  },
  {
    id: 'modified_gravity_halo',
    name: 'Modified Gravity + Halo',
    formula: 'sqrt((G * M) / (r + a) + k * r)',
    description: 'v = sqrt(GM/(r+a) + kr)',
    physicalMeaning: 'Combines softened gravity with dark matter halo — two free parameters for fitting.',
    category: 'modified_gravity'
  },
  {
    id: 'transition',
    name: 'Transition Model',
    formula: 'sqrt((G * M) / r) * (1 + k * r / 100)',
    description: 'v = sqrt(GM/r) * (1 + kr/100)',
    physicalMeaning: 'Newtonian at small r, gradually boosted at large r — models transition to dark matter dominated regime.',
    category: 'transition'
  },
  {
    id: 'mond_like',
    name: 'MOND-inspired',
    formula: 'sqrt(sqrt((G * M) / r) * a + (G * M) / r)',
    description: 'v = sqrt(sqrt(GM/r)*a + GM/r)',
    physicalMeaning: 'Inspired by Modified Newtonian Dynamics — gravity transitions at a critical acceleration scale.',
    category: 'modified_gravity'
  },
  {
    id: 'log_halo',
    name: 'Logarithmic Halo',
    formula: 'sqrt((G * M) / r + k * log(1 + r / a))',
    description: 'v = sqrt(GM/r + k*ln(1+r/a))',
    physicalMeaning: 'NFW-inspired logarithmic dark matter halo profile — widely used in astrophysics.',
    category: 'additive'
  },
  {
    id: 'accel_floor',
    name: 'Acceleration Floor (a\u2080)',
    formula: 'sqrt(r * sqrt((G * M / r^2)^2 + a^2 * (G * M / r^2) / (G * M / r^2 + a)))',
    description: 'g_obs = \u221A(g_bar\u00B2 + a\u2080\u00B2 \u00D7 g_bar/(g_bar+a\u2080)), V = \u221A(g_obs \u00D7 r)',
    physicalMeaning: 'Acceleration-based model: the universe enforces a minimum acceleration floor a\u2080. Recovers Newton at high g, gives RAR + flat rotation at low g. Use a = a\u2080 in (km/s)\u00B2/kpc.',
    category: 'acceleration'
  },
  {
    id: 'accel_floor_cosmo',
    name: 'Cosmic Floor (cH\u2080/2\u03C0)',
    formula: 'sqrt(r * sqrt((G * M / r^2)^2 + 3702^2 * (G * M / r^2) / (G * M / r^2 + 3702)))',
    description: 'a\u2080 = cH\u2080/2\u03C0 \u2248 3702 (km/s)\u00B2/kpc — zero free parameters',
    physicalMeaning: 'The acceleration floor is fixed by cosmology: a\u2080 = cH\u2080/2\u03C0. No adjustable parameters beyond G and M. Predicts flat rotation, RAR, and Tully-Fisher simultaneously.',
    category: 'acceleration'
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
  benchmarkResult: BenchmarkResult | null;
  isBenchmarking: boolean;

  uploadDataset: (file: File) => Promise<void>;
  loadSampleDataset: (name: string) => void;
  loadAllSamples: () => void;
  removeDataset: (id: string) => void;
  toggleDatasetActive: (id: string) => void;
  updateModelParams: (params: Partial<ModelParams>) => void;
  toggleLayer: (layer: 'observed' | 'newtonian' | 'custom' | 'discovery') => void;
  applyPreset: (preset: FormulaPreset) => void;
  autoOptimize: () => void;
  runFullBenchmark: () => void;

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
  a: 3702,
  formula: "sqrt(r * sqrt((G * M / r^2)^2 + a^2 * (G * M / r^2) / (G * M / r^2 + a)))",
};

const SAMPLE_DATASETS: Record<string, Dataset> = SPARC_DATASETS;

const GalaxyContext = createContext<GalaxyContextType | undefined>(undefined);

export const GalaxyProvider = ({ children }: { children: ReactNode }) => {
  const firstGalaxy = SAMPLE_DATASETS["NGC 3198"] || Object.values(SAMPLE_DATASETS)[0];
  const [datasets, setDatasets] = useState<Dataset[]>([firstGalaxy]);
  const [activeDatasetIds, setActiveDatasetIds] = useState<string[]>([firstGalaxy.id]);
  const [modelParams, setModelParams] = useState<ModelParams>(defaultParams);
  
  const [showObserved, setShowObserved] = useState(true);
  const [showNewtonian, setShowNewtonian] = useState(true);
  const [showCustom, setShowCustom] = useState(true);
  const [discoveryMode, setDiscoveryMode] = useState(false);
  const [isOptimizing, setIsOptimizing] = useState(false);
  const [optimizationLog, setOptimizationLog] = useState<string[]>([]);
  const [benchmarkResult, setBenchmarkResult] = useState<BenchmarkResult | null>(null);
  const [isBenchmarking, setIsBenchmarking] = useState(false);

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

  const optimizeFormulaOnDataset = useCallback((formula: string, ds: Dataset, baseG: number): { k: number; a: number; M: number; mse: number } => {
    let bestK = 0, bestA = 1, bestM = 1e11, bestMSE = Infinity;

    const computeDS_MSE = (f: string, p: ModelParams, data: DataPoint[]): number => {
      let sum = 0, count = 0;
      data.forEach(pt => {
        const pred = evaluateFormulaWithParams(f, p, pt.r);
        if (pred !== null) { sum += (pt.v - pred) ** 2; count++; }
      });
      return count > 0 ? sum / count : Infinity;
    };

    const kVals = Array.from({ length: 25 }, (_, i) => i * 8);
    const aVals = [0.1, 0.5, 1, 2, 3, 5, 8, 10, 15, 20];
    const mVals = [1e9, 5e9, 1e10, 5e10, 1e11, 2e11, 5e11, 1e12];

    for (const tM of mVals) {
      for (const tK of kVals) {
        for (const tA of aVals) {
          const params: ModelParams = { G: baseG, M: tM, k: tK, a: tA, formula: formula };
          const mse = computeDS_MSE(formula, params, ds.data);
          if (mse < bestMSE) { bestMSE = mse; bestK = tK; bestA = tA; bestM = tM; }
        }
      }
    }

    const fineK = Array.from({ length: 17 }, (_, i) => Math.max(0, bestK - 8 + i));
    const fineA = Array.from({ length: 11 }, (_, i) => Math.max(0.01, bestA - 2.5 + i * 0.5));
    const fineM = Array.from({ length: 9 }, (_, i) => bestM * (0.6 + i * 0.1));

    for (const tM of fineM) {
      for (const tK of fineK) {
        for (const tA of fineA) {
          const params: ModelParams = { G: baseG, M: tM, k: tK, a: tA, formula: formula };
          const mse = computeDS_MSE(formula, params, ds.data);
          if (mse < bestMSE) { bestMSE = mse; bestK = tK; bestA = tA; bestM = tM; }
        }
      }
    }

    return { k: bestK, a: bestA, M: bestM, mse: bestMSE };
  }, [evaluateFormulaWithParams]);

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
          const testParams: ModelParams = { ...modelParams, k: testK, a: testA, M: testM };
          const mse = computeMSEForParams(modelParams.formula, testParams);
          iterations++;
          if (mse < bestMSE) { bestMSE = mse; bestK = testK; bestA = testA; bestM = testM; }
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
          const testParams: ModelParams = { ...modelParams, k: testK, a: testA, M: testM };
          const mse = computeMSEForParams(modelParams.formula, testParams);
          iterations++;
          if (mse < bestMSE) { bestMSE = mse; bestK = testK; bestA = testA; bestM = testM; }
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

  const runFullBenchmark = useCallback(() => {
    setIsBenchmarking(true);

    const benchmarkData = datasets.filter(d => activeDatasetIds.includes(d.id));
    if (benchmarkData.length === 0) { setIsBenchmarking(false); return; }

    const baseG = modelParams.G;
    const formulasToTest = FORMULA_PRESETS.filter(p => p.id !== 'newtonian');

    const newtonianBaselineCache = new Map<string, { M: number; mse: number }>();
    benchmarkData.forEach(ds => {
      const newtOpt = optimizeFormulaOnDataset('sqrt((G * M) / r)', ds, baseG);
      newtonianBaselineCache.set(ds.id, { M: newtOpt.M, mse: newtOpt.mse });
    });

    const formulaResults: FormulaBenchmark[] = formulasToTest.map(preset => {
      const galaxyResults: BenchmarkGalaxyResult[] = benchmarkData.map(ds => {
        const opt = optimizeFormulaOnDataset(preset.formula, ds, baseG);

        const newtBaseline = newtonianBaselineCache.get(ds.id)!;
        const newtMSE = newtBaseline.mse;
        const newtParams: ModelParams = { G: baseG, M: newtBaseline.M, k: 0, a: 1, formula: 'sqrt((G * M) / r)' };

        const mid = Math.floor(ds.data.length / 2);
        const innerData = ds.data.slice(0, mid);
        const outerData = ds.data.slice(mid);

        const optParams: ModelParams = { G: baseG, M: opt.M, k: opt.k, a: opt.a, formula: preset.formula };
        let mseIn = 0, mseOut = 0, mseInN = 0, mseOutN = 0;
        let cIn = 0, cOut = 0;

        innerData.forEach(p => {
          const pc = evaluateFormulaWithParams(preset.formula, optParams, p.r);
          const pn = evaluateFormulaWithParams('sqrt((G * M) / r)', newtParams, p.r);
          if (pc !== null) { mseIn += (p.v - pc) ** 2; cIn++; }
          if (pn !== null) { mseInN += (p.v - pn) ** 2; }
        });
        outerData.forEach(p => {
          const pc = evaluateFormulaWithParams(preset.formula, optParams, p.r);
          const pn = evaluateFormulaWithParams('sqrt((G * M) / r)', newtParams, p.r);
          if (pc !== null) { mseOut += (p.v - pc) ** 2; cOut++; }
          if (pn !== null) { mseOutN += (p.v - pn) ** 2; }
        });
        mseIn = cIn > 0 ? mseIn / cIn : 0;
        mseOut = cOut > 0 ? mseOut / cOut : 0;
        mseInN = cIn > 0 ? mseInN / cIn : 0;
        mseOutN = cOut > 0 ? mseOutN / cOut : 0;

        const improvPct = newtMSE > 0 ? ((newtMSE - opt.mse) / newtMSE * 100) : 0;
        const innerImpr = mseInN > 0 ? ((mseInN - mseIn) / mseInN * 100) : 0;
        const outerImpr = mseOutN > 0 ? ((mseOutN - mseOut) / mseOutN * 100) : 0;

        return {
          galaxyName: ds.name,
          galaxyId: ds.id,
          mseNewton: newtMSE,
          mseCustom: opt.mse,
          improvementPct: improvPct,
          bestK: opt.k,
          bestA: opt.a,
          bestM: opt.M,
          mseInner: mseIn,
          mseOuter: mseOut,
          innerImprovement: innerImpr,
          outerImprovement: outerImpr,
          pointCount: ds.data.length
        };
      });

      const kValues = galaxyResults.map(g => g.bestK);
      const avgK = kValues.reduce((s, v) => s + v, 0) / kValues.length;
      const kVariance = kValues.reduce((s, v) => s + (v - avgK) ** 2, 0) / kValues.length;
      const kStdDev = Math.sqrt(kVariance);

      const avgMSE = galaxyResults.reduce((s, g) => s + g.mseCustom, 0) / galaxyResults.length;
      const avgImpr = galaxyResults.reduce((s, g) => s + g.improvementPct, 0) / galaxyResults.length;
      const wins = galaxyResults.filter(g => g.improvementPct > 5).length;
      const avgA = galaxyResults.reduce((s, g) => s + g.bestA, 0) / galaxyResults.length;
      const avgM = galaxyResults.reduce((s, g) => s + g.bestM, 0) / galaxyResults.length;

      return {
        formulaId: preset.id,
        formulaName: preset.name,
        formula: preset.formula,
        avgMSE,
        avgImprovement: avgImpr,
        winsCount: wins,
        galaxyResults,
        avgK,
        kStdDev,
        avgA,
        avgM
      };
    });

    formulaResults.sort((a, b) => b.avgImprovement - a.avgImprovement);
    const best = formulaResults[0];

    let kConsistency = 'N/A';
    if (best) {
      const cv = best.avgK > 0 ? (best.kStdDev / best.avgK * 100) : 100;
      if (cv < 20) kConsistency = `Highly consistent (CV=${cv.toFixed(0)}%). k is stable across galaxies.`;
      else if (cv < 50) kConsistency = `Moderately consistent (CV=${cv.toFixed(0)}%). Some variation across galaxies.`;
      else kConsistency = `Inconsistent (CV=${cv.toFixed(0)}%). k varies significantly — may be galaxy-dependent.`;
    }

    setBenchmarkResult({
      formulas: formulaResults,
      bestFormula: best?.formulaName || 'N/A',
      kConsistency,
      timestamp: new Date().toISOString()
    });
    setIsBenchmarking(false);
  }, [datasets, activeDatasetIds, modelParams.G, evaluateFormulaWithParams, optimizeFormulaOnDataset]);

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
      isOptimizing, optimizationLog, benchmarkResult, isBenchmarking,
      uploadDataset, loadSampleDataset, loadAllSamples, removeDataset, toggleDatasetActive,
      updateModelParams, toggleLayer, applyPreset, autoOptimize, runFullBenchmark,
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
