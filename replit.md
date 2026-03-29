# Workspace

## Overview

pnpm workspace monorepo using TypeScript. Each package manages its own dependencies.

## Stack

- **Monorepo tool**: pnpm workspaces
- **Node.js version**: 24
- **Package manager**: pnpm
- **TypeScript version**: 5.9
- **API framework**: Express 5
- **Database**: PostgreSQL + Drizzle ORM
- **Validation**: Zod (`zod/v4`), `drizzle-zod`
- **API codegen**: Orval (from OpenAPI spec)
- **Build**: esbuild (CJS bundle)

## Structure

```text
artifacts-monorepo/
├── artifacts/              # Deployable applications
│   ├── api-server/         # Express API server
│   └── galaxy-analyzer/    # Galaxy Rotation Curve Analyzer (React + Vite)
├── lib/                    # Shared libraries
│   ├── api-spec/           # OpenAPI spec + Orval codegen config
│   ├── api-client-react/   # Generated React Query hooks
│   ├── api-zod/            # Generated Zod schemas from OpenAPI
│   └── db/                 # Drizzle ORM schema + DB connection
├── scripts/                # Utility scripts (single workspace package)
│   └── src/                # Individual .ts scripts, run via `pnpm --filter @workspace/scripts run <script>`
├── pnpm-workspace.yaml     # pnpm workspace (artifacts/*, lib/*, lib/integrations/*, scripts)
├── tsconfig.base.json      # Shared TS options (composite, bundler resolution, es2022)
├── tsconfig.json           # Root TS project references
└── package.json            # Root package with hoisted devDeps
```

## Galaxy Rotation Curve Analyzer

A professional dark matter research tool built with React + Vite + Tailwind CSS.

### Features
- **Dataset Upload**: Drag & drop .csv/.dat files, parse radius/velocity columns, data preview table
- **175 SPARC Galaxies**: Real rotation curves from Lelli, McGaugh & Schombert (2016) SPARC database. 3,391 data points across 175 nearby galaxies with searchable catalog.
- **Visualization**: Interactive Recharts scatter + line plot with observed data, Newtonian model, custom model overlays
- **Model Builder**: Custom math formula input using mathjs, parameter sliders (G, M, k, a), 8 formula presets (Dark Halo Linear/Flat, Modified Gravity, MOND-inspired, Logarithmic Halo, etc.)
- **Auto-Parameter Optimizer**: 2-phase grid search (coarse + fine-tune) for k, a, M
- **Discovery Mode**: Highlights anomalies where model deviates >15% from observed data
- **Per-Galaxy MSE Table**: Individual MSE comparison per galaxy with winner indicator
- **Residual Analysis**: Residual chart showing model fit quality across radii
- **Generalization Score**: Percentage of galaxies where custom model beats Newtonian
- **Research Lab**: Full benchmark (all 7 formulas × all 15 galaxies), regional analysis (inner vs outer), k-consistency scoring (CV metric), units normalization panel
- **Export**: PNG graph export (html-to-image), JSON model/benchmark export

### Tech
- React 19, Vite, Tailwind CSS v4, Recharts, mathjs, PapaParse, Framer Motion
- Dark space theme with glassmorphism cards
- Fully client-side (no backend needed)

### Key Files
- `artifacts/galaxy-analyzer/src/hooks/use-galaxy.tsx` — Core state management (datasets, model params, chart data, benchmark engine, 15 sample galaxies, formula presets)
- `artifacts/galaxy-analyzer/src/pages/dashboard.tsx` — Dashboard overview
- `artifacts/galaxy-analyzer/src/pages/upload.tsx` — Dataset management with preview table (15 sample galaxies)
- `artifacts/galaxy-analyzer/src/pages/analysis.tsx` — Visualization with layer toggles, residual chart, per-galaxy MSE, discovery mode
- `artifacts/galaxy-analyzer/src/pages/models.tsx` — Formula builder, parameter tuning, auto-optimizer, formula library
- `artifacts/galaxy-analyzer/src/pages/research.tsx` — Research Lab: units panel, full benchmark, formula rankings, regional analysis, k-consistency
- `artifacts/galaxy-analyzer/src/components/ui/glass-card.tsx` — Glassmorphism card component
- `artifacts/galaxy-analyzer/src/components/layout.tsx` — Sidebar navigation layout

### Design System
- Background: #0a0e1a (deep space), glassmorphism cards (bg-white/5, backdrop-blur)
- Colors: cyan/teal (observed data), orange dashed (Newtonian), purple solid (custom model), amber (anomalies)
- Fonts: Space Grotesk (display), Inter (sans), JetBrains Mono (monospace)
- G = 4.3009e-6 kpc·(km/s)²/M_sun — real physical constant
- ModelParams interface: explicit fields only (G, M, k, a, formula) — no index signature
- Formula validation uses real mathjs evaluation at r=10 before persisting

## TypeScript & Composite Projects

Every package extends `tsconfig.base.json` which sets `composite: true`. The root `tsconfig.json` lists all packages as project references. This means:

- **Always typecheck from the root** — run `pnpm run typecheck` (which runs `tsc --build --emitDeclarationOnly`). This builds the full dependency graph so that cross-package imports resolve correctly. Running `tsc` inside a single package will fail if its dependencies haven't been built yet.
- **`emitDeclarationOnly`** — we only emit `.d.ts` files during typecheck; actual JS bundling is handled by esbuild/tsx/vite...etc, not `tsc`.
- **Project references** — when package A depends on package B, A's `tsconfig.json` must list B in its `references` array. `tsc --build` uses this to determine build order and skip up-to-date packages.

## Root Scripts

- `pnpm run build` — runs `typecheck` first, then recursively runs `build` in all packages that define it
- `pnpm run typecheck` — runs `tsc --build --emitDeclarationOnly` using project references

## Packages

### `artifacts/api-server` (`@workspace/api-server`)

Express 5 API server. Routes live in `src/routes/` and use `@workspace/api-zod` for request and response validation and `@workspace/db` for persistence.

- Entry: `src/index.ts` — reads `PORT`, starts Express
- App setup: `src/app.ts` — mounts CORS, JSON/urlencoded parsing, routes at `/api`
- Routes: `src/routes/index.ts` mounts sub-routers; `src/routes/health.ts` exposes `GET /health` (full path: `/api/health`)
- Depends on: `@workspace/db`, `@workspace/api-zod`
- `pnpm --filter @workspace/api-server run dev` — run the dev server
- `pnpm --filter @workspace/api-server run build` — production esbuild bundle (`dist/index.cjs`)
- Build bundles an allowlist of deps (express, cors, pg, drizzle-orm, zod, etc.) and externalizes the rest

### `lib/db` (`@workspace/db`)

Database layer using Drizzle ORM with PostgreSQL. Exports a Drizzle client instance and schema models.

- `src/index.ts` — creates a `Pool` + Drizzle instance, exports schema
- `src/schema/index.ts` — barrel re-export of all models
- `src/schema/<modelname>.ts` — table definitions with `drizzle-zod` insert schemas (no models definitions exist right now)
- `drizzle.config.ts` — Drizzle Kit config (requires `DATABASE_URL`, automatically provided by Replit)
- Exports: `.` (pool, db, schema), `./schema` (schema only)

Production migrations are handled by Replit when publishing. In development, we just use `pnpm --filter @workspace/db run push`, and we fallback to `pnpm --filter @workspace/db run push-force`.

### `lib/api-spec` (`@workspace/api-spec`)

Owns the OpenAPI 3.1 spec (`openapi.yaml`) and the Orval config (`orval.config.ts`). Running codegen produces output into two sibling packages:

1. `lib/api-client-react/src/generated/` — React Query hooks + fetch client
2. `lib/api-zod/src/generated/` — Zod schemas

Run codegen: `pnpm --filter @workspace/api-spec run codegen`

### `lib/api-zod` (`@workspace/api-zod`)

Generated Zod schemas from the OpenAPI spec (e.g. `HealthCheckResponse`). Used by `api-server` for response validation.

### `lib/api-client-react` (`@workspace/api-client-react`)

Generated React Query hooks and fetch client from the OpenAPI spec (e.g. `useHealthCheck`, `healthCheck`).

### `scripts` (`@workspace/scripts`)

Utility scripts package. Each script is a `.ts` file in `src/` with a corresponding npm script in `package.json`. Run scripts via `pnpm --filter @workspace/scripts run <script>`. Scripts can import any workspace package (e.g., `@workspace/db`) by adding it as a dependency in `scripts/package.json`.
