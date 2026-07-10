# PAY-4 cross-tool benchmark

PAY-4 is a reproducible positioning measurement, not a correctness gate. It compares:

- one Chen2020 SPM cell over a 1C discharge followed by CCCV charge against PyBaMM 26.x using `IDAKLUSolver` (`rtol=1e-6`, `atol=1e-8`); and
- Chen2020 SPM 16s4p and 1s100p packs over 60 ten-second 5 A/cell discharge steps against liionpack 0.3.12 using `CasadiManager`.

The pack cell resistance is 100 micro-ohm in both tools. liionpack also requires 0.1 nano-ohm connector resistors; SLIDE's corresponding connections are ideal. Final voltage is included in every result so a fast but different workload cannot pass unnoticed.

## Run

Build the two Release executables, create a separate Python environment for liionpack, then run the driver from the repository root:

```powershell
cmake --build build-release --target benchmark_PAY4_slide_single benchmark_PAY4_slide_pack
py -3.11 -m venv $env:TEMP\slide-liionpack-bench
& $env:TEMP\slide-liionpack-bench\Scripts\pip.exe install liionpack==0.3.12 pandas==2.2.3
.venv\Scripts\python.exe benchmark\benchmark_PAY4_cross_tool.py `
  --liionpack-python $env:TEMP\slide-liionpack-bench\Scripts\python.exe `
  --output PAY4-result.json
```

liionpack 0.3.12 currently installs PyBaMM 24.9, whereas the single-cell reference uses the Phase-7 target PyBaMM 26.6.2.0. This unavoidable version split is emitted in the JSON. pandas is capped below 3 because liionpack 0.3.12 mutates a mask through an API that became read-only in pandas 3.

## Timing boundary

SLIDE setup covers parameter absorption, spectral compilation, batch allocation, topology compilation, and solver configuration. Its solve measurement restores a preallocated checkpoint before each repetition. PyBaMM setup covers model/parameter/experiment construction plus `build_for_experiment`; the first solve is reported separately from warm repeated solves because IDAKLU initializes more cached state on that call. liionpack setup uses its `setup_only=True` path; solve times its own protocol loop and output materialization.

liionpack necessarily records its output arrays while the SLIDE pack harness advances state without a recorder. Both setup and solve are therefore reported, but the pack comparison represents end-to-end tool positioning rather than a kernel-only ratio. Rerun on a quiet machine before quoting results.

## 2026-07-10 development-machine run

Windows 11, Python 3.13.9, PyBaMM 26.6.2.0, liionpack 0.3.12/PyBaMM 24.9.0, five SLIDE repetitions:

| Case | SLIDE solve median | Comparator solve | Speedup | Conservative speedup |
|---|---:|---:|---:|---:|
| single cell vs warm IDAKLU | 0.682 ms | 12.373 ms | 18.13x | n/a |
| 16s4p vs liionpack | 0.751 ms | 1.024 s | 1,363x | 1,165x |
| 1s100p vs liionpack | 1.470 ms | 1.705 s | 1,160x | 1,111x |

Cold setup-plus-first-solve is 71.4x faster than PyBaMM. SLIDE/liionpack final pack voltages differ by 14.38 mV for 16 series cells (0.90 mV/cell) and 0.90 mV for the 100p pack.
