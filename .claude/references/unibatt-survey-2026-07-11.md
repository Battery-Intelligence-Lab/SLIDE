# unibatt / glide-wasm — survey digest for SLIDE M10 (WASM) + M11 (Studio) (scout, 2026-07-11)

> Source repo: `C:\D\git\unibatt` (Rust/WASM, Volkan). Read the named files before designing M10/M11.

## 1. Colour palette (exact hex + paths)

**CSS theme tokens** — `glide-wasm/web/src/styles/main.css` (`:root`, lines 5-27):
- Dark backgrounds: `--bg-body: #0f0f23`, `--bg-panel: #1a1a2e`, `--bg-input`/`--bg-plot: #16213e`
- Text: `--text-primary: #e0e0e0`, `--text-secondary: #a0a0b8`, `--text-muted: #8888a8`
- Borders: `--border: #2a2a4a`, `--border-focus: #0072B2`
- Accents (Wong colourblind-safe): `--accent-green: #009E73`, `--accent-blue: #0072B2`,
  `--accent-orange: #D55E00`, `--accent-cyan: #00b4d8`
- Oxford identity (chrome only): `--ox-blue: #002147`, `--ox-blue-600: #122f53`, `--ox-blue-300: #49B6FF`
- `--radius: 8px`; fonts `'JetBrains Mono'/'Fira Code'/'Cascadia Code'` mono, `'Segoe UI'/system-ui` sans

**Light theme** (`[data-theme="light"]`, lines 30-40): `--bg-body: #f5f5f5`, panels/input `#ffffff`,
`--bg-plot: #f8f8f8`, `--text-primary: #1a1a2e`. Secondary button `#005f99`, orange hover `#e66b00`,
blue hover `#0088d4`.

**Plotly data-series palette (Wong)** — `glide-wasm/web/src/components/chart.js:6`:
`COLORS = ['#0072B2', '#D55E00', '#009E73', '#F0E442', '#CC79A7', '#56B4E9', '#E69F00']`
Chart layout (chart.js:10-16 dark / 21-27 light): `paper_bgcolor:#1a1a2e`, `plot_bgcolor:#16213e`,
gridcolor `#2a2a4a`, zeroline/linecolor `#3a3a5a`, font `#e0e0e0`; light = `#ffffff`/`#f8f8f8`/`#e0e0e0`/`#cccccc`.

Design rule stated in the CSS: **Oxford blue for chrome only, never for data**; Wong palette reserved for
data series.

## 2. WASM architecture

- Toolchain: `wasm-pack` + `wasm-bindgen` (NOT trunk/emscripten/yew/leptos). Build:
  `wasm-pack build crates/glide-wasm-bindgen --target web --out-dir ../../web/pkg` (`scripts/build-wasm.sh`;
  echoed in `glide-wasm/CLAUDE.md`).
- Bindgen crate `glide-wasm/crates/glide-wasm-bindgen/Cargo.toml`: `crate-type = ["cdylib","rlib"]`,
  edition 2024; deps wasm-bindgen 0.2, serde-wasm-bindgen 0.6, js-sys, web-sys, getrandom (`wasm_js`),
  console_error_panic_hook, rand/rand_chacha/rand_distr, faer 0.22 (no_std linalg); core `glide` crate
  with `default-features = false`.
- `glide-wasm/Cargo.toml` is its own nested workspace (`[profile.release] opt-level="s", lto=true`) —
  deliberately NOT a member of the root workspace so WASM/no_std deps don't pollute the native build.
- Exposed API (`crates/glide-wasm-bindgen/src/lib.rs`, JsValue via serde): `init`, `generate_time_series`,
  `extract_segments`, `detect_sign_convention`, `optimize`, `evaluate_params`, `evaluate_params_diag`,
  `forward_backward_pass`, `get_default_ocv`, `set_ocv_curve`, `has_custom_ocv`. Internal modules:
  `battery_sim.rs` (SimpleECMBattery, UDDS/pulse/CC/calendar), `optimizer.rs` (DifferentialEvolution,
  LbfgsB), `schimpe.rs` (Schimpe LFP aging), `segment_extract.rs`.
- Frontend: **vanilla JS ES modules**, no framework. Entry `web/src/main.js` hand-builds tabs/panels via
  `createElement`. `wasm-loader.js` dynamically imports `../pkg/glide_wasm_bindgen.js` with graceful mock
  fallback if pkg not built.
- Charting: Plotly.js (`plotly-basic-2.35.2` via CDN, `web/index.html:9`); YAML via js-yaml CDN; parquet
  via hyparquet npm.
- **Web workers** offload heavy WASM calls: `web/src/worker.js`, `eval-worker.js`, `multi-worker.js`
  (spawned by fitting.js). Cancel = `worker.terminate()` because COOP/COEP headers (needed for
  SharedArrayBuffer atomics) break CDN scripts — noted in `vite.config.js:14-16`.

## 3. What the app does today

"GLiDE-WASM — GP-ECM Battery Diagnostics": 7-tab wizard (`main.js:17-25`): Data input (synthetic
calendar/cc_cycling/UDDS generation or CSV/parquet/YAML upload) → OCV (view/edit, default LFP,
`data/ocv-lfp.json`) → Segments (extraction + sign-convention auto-detect with modal flip prompt) →
Fitting (DE + L-BFGS-B in workers, progress bar) → Estimation (GP-ECM forward/backward, NLL) → Sweep →
Documentation. Forms + sliders (`components/slider.js`) + Plotly + tables + drag-drop file zones +
dark/light toggle persisted in localStorage. All simulation client-side in WASM; no backend.

## 4. Build / deploy patterns worth imitating

- Two deliverables: (a) Vite dev/prod site — `web/package.json` scripts; Vite 6 + `vite-plugin-wasm` +
  `vite-plugin-top-level-await`, `build.target:'esnext'`, `worker.format:'es'` (`web/vite.config.js`).
  (b) **Single self-contained offline HTML** — `scripts/build-single-html.js` (esbuild-bundles JS/CSS,
  base64-embeds the .wasm, inlines Plotly+js-yaml) → `dist/glide.html` (~2.3 MB, double-click offline;
  OCV CSVs decimated; optional `--recipient` watermark). The most reproducible pattern for a shareable
  design tool.
- Workspace layout: root Cargo workspace (`glide`, `pyglide`, resolver 3, shared workspace deps);
  glide-wasm nested independent workspace.
- CI (`.github/workflows/CI.yml`): maturin-generated, covers ONLY Python wheels + Rust tests; **no CI for
  glide-wasm/web** — wasm + site built manually via the two scripts.
- Docs to mine: `glide-wasm/CLAUDE.md`, `glide-wasm/CHECKLIST.md` (living invariants),
  `glide-wasm/docs/superpowers/specs/2026-03-31-glide-wasm-migration-design.md`.
- Not present: trunk, emscripten, yew/leptos/react, Tailwind, WebGPU/three.js, custom canvas rendering.
