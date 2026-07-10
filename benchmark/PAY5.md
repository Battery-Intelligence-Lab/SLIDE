# PAY-5: CUDA batch positioning

P8-G2 uses `benchmark_PAY5_cuda` to compare the same base-isothermal `nch=12`
exact-modal SPM composition on the CPU SoA batch and the optional CUDA backend.
The CUDA kernel maps one lane to one thread. Setup, state upload, warmed solve,
and state download are timed separately; current H2D is warmed before the solve.

The workload is 100,000 lanes and 360 × 10 s 1C steps. It uses a symmetric
180-step discharge / 180-step charge so every lane stays inside the OCV domain,
with SOC 0.65–0.90, ±5% solid diffusivity, and ±5% contact resistance. Three
repetitions alternate CPU/GPU ordering and report medians. Timing is admitted
only when the final full arena satisfies the registered state band and terminal
voltage differs by no more than 2 µV.

Qualified development-machine result (2026-07-10):

- GPU: NVIDIA RTX 4000 Ada Generation, compute capability 8.9, 20,475 MiB
- CUDA toolkit/compiler: 13.0.48; driver 596.72
- CPU setup: 0.03491 s; GPU setup: 0.20710 s
- H2D: 0.00317 s; D2H: 0.01278 s
- CPU warmed solve median: 8.60211 s
- GPU warmed solve median: 0.36075 s
- solve speedup: **23.85×** (target ≥5×, useful floor ≥2×)
- maximum arena error: 1.78e-14; maximum normalized gate usage: 9.00e-6
- maximum terminal-voltage error: 2.66e-15 V
- one 70,492,416-byte device allocation; zero device-wide synchronizations

The result qualifies consumer-Ada f64 performance rather than promising the same
ratio on every GPU. Reproduce with a Release CUDA build and run
`benchmark_PAY5_cuda`; the executable emits the complete machine-readable JSON.
