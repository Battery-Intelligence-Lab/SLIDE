# Handoff — Phase 9B parser fuzz gate complete (2026-07-10)

## Outcome

P9-G2 is complete. Experiment, BPX, and liionpack CSV now have bounded
raw-byte libFuzzer targets with committed corpora/dictionaries, poisoned
failure-atomicity sentinels, independent success invariants, Debug/Release
ASan+UBSan coverage, explicit size-limit replays, and a matching GitHub
workflow. P9-G1, P9-G3, and P9-G4 remain open.

## Architecture

- `cmake/SlideCoreTarget.cmake` is the single production source-list helper.
  Fuzzers compile `slide_core_fuzz` separately, so sanitizer and Windows
  static-CRT/iterator ABI settings cannot contaminate `slide_core` consumers.
- `FuzzSupport.hpp` deeply compares every public output field. Failure must
  preserve a populated poison Experiment/ParameterSet/topology; success must
  remove all poison and satisfy parser-output invariants. Named callback poison
  targets carry exact markers, while an opaque reference/volatile IEEE-bit
  classifier keeps finiteness checks valid under Release fast-math.
- The netlist oracle does not call `validateElectricalNetlist`. It independently
  proves cell/path/location bijection, endpoints/resistors, connectivity,
  sparsity, ladder classification/orientation, and imported thermal metadata.
- Fuzz targets disable C++ module scanning because no target contains modules;
  this avoids a false Ubuntu dependency on the separately packaged
  `clang-scan-deps` executable.

## Adversarial evidence

- Experiment append mutation: exit 77 on the documented seed.
- BPX merge mutation: exit 77 on the valid SPM seed.
- Netlist stray thermal-scratch mutation: exit 77 on the parallel seed.
- Release Experiment NaN-publication mutation: exit 77 on the documented seed;
  the mutation used a quiet-NaN bit pattern and the reverted seed replay passes.
- A normal trailing newline originally prevented the Experiment corpus from
  reaching success. The adapter now suppresses only that synthetic final empty
  line; interior empty steps remain parser input.
- A new `Hold at 4.2 V until 3.8 V` seed then found P9-B29: parsing published a
  segment rejected by the runner's own semantic validator. The ordinary test
  was red with 7 failures before parsing adopted `validSegment`; Debug and
  Release now pass all 203 Experiment assertions.
- Linux Debug and fast-math Release both replay a 65,537-byte Experiment step
  and 4,194,305-byte BPX/netlist documents under ASan+UBSan+LSan.

Exact 60-second Linux campaign counts:

| Build | Experiment | BPX | NetlistCsv |
|-------|-----------:|----:|-----------:|
| Debug | 149,399 | 73,981 | 94,933 |
| Release | 291,223 | 546,881 | 556,946 |

No crash, timeout, UB, or leak was reported; the slowest input was below the
one-second reporting resolution. Windows Clang 21 Debug/Release ASan+UBSan
campaigns and limit replays also pass. Windows alone uses
`detect_container_overflow=0` because the prebuilt libFuzzer runtime and MSVC
STL annotation ABI disagree; Linux retains full container/leak checking.

Ordinary core-only Debug/Release smoke tests pass 1/1, full-tree Debug/Release
`slide_core` builds pass, and a combined Windows build links/runs the ordinary
core smoke beside the instrumented fuzzer. Workflow YAML is syntax-checked.
The hosted Debug/Release jobs are committed but cannot run until pushed.

## Next

Start with the ThreadPool audit already identified by review: production has no
parallel integration, failure selection is scheduling-dependent rather than
lowest-index deterministic, const-lvalue callbacks fail to compile, and the
core target does not explicitly link `Threads::Threads`. Add red tests first,
fix ownership/determinism/linkage, then run the P9-G1 TSan lane. Continue to
P9-G3 measured Status-branch coverage and the remaining CUDA/binding/
sensitivity subsystem passes before closing P9-G4 or beginning 9C.
