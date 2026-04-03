# Workspace

## Overview

This pnpm workspace monorepo, built with TypeScript, serves as a professional dark matter research tool. Its core purpose is to analyze galaxy rotation curves using real-world data and advanced modeling. The project aims to provide a comprehensive platform for researchers to:

1.  **Analyze Galaxy Data**: Ingest and visualize astronomical datasets, specifically focusing on galaxy rotation curves.
2.  **Model Dark Matter Hypotheses**: Develop, test, and optimize various dark matter and modified gravity models against observed data.
3.  **Validate Theories**: Conduct rigorous statistical validation of models, perform replication studies, and explore fundamental physics questions related to the baryon-dark matter connection.

The project addresses a critical need in astrophysics for tools that facilitate the testing and refinement of dark matter theories, offering a powerful platform for scientific discovery.

## User Preferences

I prefer detailed explanations and clear communication. Please ask before making major changes. I appreciate an iterative development approach.

**CRITICAL**: The PRIMARY report is now `artifacts/galaxy-analyzer/diagnostic-report.txt`. All future analysis work and edits go here. The old theoretical report (`galaxy-rotation-curve-analysis-full-report.txt`) is ARCHIVED — read-only reference.

**Skills**: Three custom skills store all project data and context:
- `galaxy-data` — Physical constants, dataset summaries, computed test results
- `galaxy-analysis` — How to run scripts, edit reports, code patterns
- `galaxy-project` — Project context, current status, strategic direction

## System Architecture

The project is structured as a pnpm monorepo, facilitating shared libraries and independent application development.

**Monorepo Structure:**

*   `artifacts/`: Deployable applications (`api-server`, `galaxy-analyzer`).
*   `lib/`: Shared libraries (`api-spec`, `api-client-react`, `api-zod`, `db`).
*   `scripts/`: Utility scripts.

**Galaxy Rotation Curve Analyzer (Frontend - `galaxy-analyzer`):**

*   **UI/UX:** Dark space theme with glassmorphism cards. Uses `cyan/teal` for observed data, `orange dashed` for Newtonian models, `purple solid` for custom models, and `amber` for anomalies. Fonts include Space Grotesk (display), Inter (sans), and JetBrains Mono (monospace). **Fully responsive**: sidebar collapses to hamburger menu on screens < 1024px; grids use responsive breakpoints (grid-cols-1 → md:grid-cols-2/3); GlassCard padding adapts to screen size.
*   **Technical Stack:** React 19, Vite, Tailwind CSS v4, Recharts for visualization, mathjs for formula parsing, PapaParse for CSV, Framer Motion for animations. It is a fully client-side application.
*   **Key Features:**
    *   **Data Handling:** Dataset upload (.csv/.dat), parsing, and preview. Integration of 175 SPARC galaxies with 3,391 data points.
    *   **Visualization:** Interactive scatter/line plots with observed data, Newtonian, and custom model overlays.
    *   **Model Building:** Custom math formula input, parameter sliders (G, M, k, a), 8 formula presets.
    *   **Optimization:** 2-phase grid search optimizer for model parameters.
    *   **Analysis Tools:** Discovery Mode for anomalies, per-galaxy MSE table, residual analysis, generalization score.
    *   **Research Lab:** Full benchmark across all formulas and galaxies, regional analysis, k-consistency scoring, units normalization.
    *   **Advanced Research Pages:** Correlations Explorer, Theoretical Framework (derivation, Law Test, density profiles), Stress Test (outlier analysis, categories), RAR Residual Analysis (density correction fit, validation pipeline), Validation Pipeline (statistical tests, replication, Monte Carlo), Discovery (density-corrected RAR law), External Replication (LITTLE THINGS data), Dark Matter Fraction, Baryon-Halo Coupling, Redshift Lab (+ feasibility analysis + high-z survey catalog), Cluster Test, Break Test (Phase C: 4 robustness tests), Pipeline (6-method measurement pipeline simulation for JWST a₀(z) extraction).
*   **Design System:** Constants like `G = 4.3009e-6 kpc·(km/s)²/M_sun` are used. Formula validation is performed using `mathjs` before persistence.

**API Server (Backend - `api-server`):**

*   **Technical Stack:** Express 5.
*   **Structure:** Routes are defined in `src/routes/` and leverage shared Zod schemas for validation and Drizzle ORM for persistence.
*   **Core Functionality:** Provides API endpoints, including a health check.

**Database (`lib/db`):**

*   **Technology:** PostgreSQL with Drizzle ORM.
*   **Functionality:** Manages database connection and schema definitions.

**API Specification & Codegen (`lib/api-spec`, `lib/api-zod`, `lib/api-client-react`):**

*   **Specification:** OpenAPI 3.1 (`openapi.yaml`) defines the API contract.
*   **Code Generation:** Orval generates:
    *   React Query hooks and a fetch client (`lib/api-client-react`).
    *   Zod schemas for request/response validation (`lib/api-zod`).

**TypeScript Configuration:**

*   Uses TypeScript 5.9.
*   Employs composite projects (`composite: true`) and project references for efficient type-checking across packages. `tsc --build --emitDeclarationOnly` is used for type declaration generation, while esbuild/Vite handles actual JavaScript bundling.

## External Dependencies

*   **Monorepo Tool:** pnpm workspaces
*   **Package Manager:** pnpm
*   **Node.js:** Version 24
*   **TypeScript:** Version 5.9
*   **API Framework:** Express 5
*   **Database:** PostgreSQL
*   **ORM:** Drizzle ORM
*   **Validation:** Zod (`zod/v4`), `drizzle-zod`
*   **API Codegen:** Orval
*   **Build Tools:** esbuild (CJS bundle), Vite
*   **Frontend Libraries:** React 19, Tailwind CSS v4, Recharts, mathjs, PapaParse, Framer Motion
*   **Data Sources:** SPARC database (Lelli, McGaugh & Schombert 2016), LITTLE THINGS (Oh et al. 2015 via VizieR)

### Key Analysis Scripts
*   `scripts/definitive-v4.cjs` — **CURRENT DEFINITIVE PIPELINE (v4.0)**: Full sample (3755 pts), per-galaxy Y* + distance marginalization (65 grid evals/galaxy), DerSimonian-Laird hierarchical model. HEADLINE: a₀=3633, tau=0.291, I²=92.4%, ratio=1.130. Includes 11-split ratio stability test (CV=14.1%).
*   `scripts/final-measurement-v3.cjs` — v3.0 pipeline (archived baseline): 2062 subsampled pts, fixed Y*, hier. a₀=3374.
*   `scripts/fdm-coupling.cjs` — Coupling analysis: V_DM, r_DMdom, α vs Σ_bar with partial correlations, bootstrap CIs, residual tests, inner/outer splits. 6/6 tests pass.
*   `scripts/fdm-coupling-excess.cjs` — Coupling excess vs ΛCDM: 50×300 mock galaxies, Δb comparison with 10,000 bootstrap. Results: r_DMdom 5.1σ, V_DM 3.8σ, α 2.7σ. Predictive test: Σ_bar improves V_DM prediction by 16.6% over Vmax alone.
*   `scripts/fdm-excess.cjs` — Discovery proof: Δb per Vmax bin (observed − ΛCDM), 10,000 bootstrap. Overall Δb = −0.034 (2.1σ).
*   `scripts/fdm-significance.cjs` — Statistical significance: bootstrap CIs, permutation test, Fisher z-test, Cohen's d, Welch's t.
*   `scripts/fdm-simulation.cjs` — NFW+exponential disk mock galaxy simulation (3 scenarios)
*   `scripts/fdm-deep-analysis.cjs` — Deep analysis: functional forms, mass scaling, halo connection

### Key Pages (14 total)
*   The Evidence (`/evidence`) — Comprehensive 7-section research narrative consolidating all findings. Sections: The Question, Anti-Correlation, Mass Independence, Halo Structure Coupling, Exceeds ΛCDM (5.1σ), Independent Replication, The Equations. Includes TOC with IntersectionObserver scroll tracking and Final Verdict banner.
*   Baryon–Halo Coupling (`/baryon-halo-coupling`) — 4 phases: Fingerprint, Inner/Outer, Residuals, Excess vs ΛCDM. 5.1σ discovery-level excess.