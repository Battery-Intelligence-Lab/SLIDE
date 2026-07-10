# Phase 9B adversarial bug ledger — 2026-07-10

This is the live evidence ledger required by P9-G4. Every confirmed defect must
have a short regression written before its production fix. Refuted candidates
remain recorded so a later pass does not revive them without new evidence.

## Confirmed defects

| ID | Severity | Subsystem | Pre-fix red evidence | Fix | Post-fix evidence | State |
|----|----------|-----------|----------------------|-----|-------------------|-------|
| P9-B01 | High | `CyclerV2` transaction | The registered throwing post-advance event callback returned `Invalid_parameters` after mutating the arena; the initial `[P9]` run contributed state/output mismatches to 13 failed assertions. | Restore the full arena on integrator, observable, control, indicator, bisection-trial, and final-event failures before any event-aware step is committed. | Debug and fast-math Release `[P9]`: 42/42 assertions, including byte-exact rollback; full Experiment binary 183/183 in each build. | FIXED |
| P9-B02 | Medium | Experiment validation/parser | Direct segments accepted invalid enum and NaN/Inf metadata; scaled durations could overflow after the only finiteness check. Tests were registered before the validator was added. | Exhaustively validate public descriptors before output/state initialization; reject invalid integrator values; recheck the duration after unit scaling. | Debug and Release direct-validation/event suite 42/42; documented grammar/atomicity suite 56/56; full subsystem 183/183 in each build. | FIXED |
| P9-B03 | High | `Recorder` thin ordering | With capacity one, steps 0, 2, duplicate 2, and backward 1 all returned success; the thinning count reached 3 instead of remaining 1. Three registered assertions failed. | Track the last accepted cadence point independently of stored slots, update it before both store and thin outcomes, and reset it on `clear()`. | Focused Debug/Release 20/20; full Recorder 98/98 in each build. | FIXED |
| P9-B04 | Medium | `Recorder` derived metadata | Bit-injected NaN elapsed time and `1e300 A / 1e-300 m²` both returned success and consumed a slot; four registered assertions failed. | Validate elapsed time and each derived density with the integer-barrier `is_finite` guard before advancing the watermark or writing storage. | Focused fast-math Release plus Debug 20/20; full Recorder 98/98 each. | FIXED |
| P9-B05 | Low | `RecorderConfig` enum | `static_cast<BackpressurePolicy>(255)` returned success and marked the recorder configured; two registered assertions failed after B03/B04 were already green. | Exhaustively accept only `stop` or `thin` before allocations/commit. | Focused Debug/Release 20/20; configuration remains atomic. | FIXED |
| P9-B06 | High | `CompressedRecording` allocation bound | A valid 64-byte header claiming `UINT64_MAX` one-value snapshots threw `length_error`; `CHECK_NOTHROW` and status assertions failed. | Bound snapshot count by `(file_size - file_header)/block_header` before any vector allocation and map `length_error` to `Invalid_parameters`. | Focused Debug/Release 67/67; no exception and failed open leaves the reader invalid. | FIXED |
| P9-B07 | Medium | compressed codec wire value | CRC-valid file and block codec fields of `256` narrowed to `CompressionCodec::none`; open returned success and two assertions failed. | Accept only the exact 32-bit wire values for `none`/`zstd` before casting; availability remains a second check. | Focused Debug/Release 67/67; full async 233/233 each. | FIXED |
| P9-B08 | Medium | `AsyncRecorder` derived density | Finite current/area whose quotient overflowed returned success and produced one written snapshot; two assertions failed. | Validate the quotient with `is_finite` before locking or publishing the slot. | Invalid enqueue, successful empty finish, zero snapshots in both builds. | FIXED |
| P9-B09 | Low | `AsyncRecorderConfig` enum | Policy value 255 configured successfully and left the object configured; two assertions failed. | Exhaustively accept `block` or `thin` before allocations/file creation. | Focused Debug/Release 67/67; configuration is unchanged. | FIXED |
| P9-B10 | High | `AsyncRecorder::configure` rollback | Faulting the measured long-path copy returned `Numerical_failure` but left `configured()==true`; one registered assertion failed. The member publish also preceded context ownership transfer. | Copy the path into a local before all resource creation; swap only during the no-throw member-commit phase. | Allocation executable Debug/Release 15/15; failed object remains unconfigured and then configures/finishes successfully. | FIXED |

## Refuted candidates

| Candidate | Adversarial evidence | Result |
|-----------|----------------------|--------|
| Cycler chooses a later event when two roots fall in one accepted step. | Independent indicators at 0.25 s and 0.75 s with a 1 s trial select the named 0.25 s event. | REFUTED |
| Cycler misses an event exactly on a step boundary. | Indicator root exactly at 1 s with a 1 s step terminates at exactly 1 s. | REFUTED |
| Pack topology accepts empty/zero groups or duplicate thermal endpoints. | Compile-time validation paths reject the malformed descriptors before workspace mutation. | REFUTED (read-only audit; mechanised tests pending ledger expansion) |
| Ordinary finite thermal spans and one-cell/zero-thermal-edge packs are intrinsically invalid. | Compile validation and the scalar topology algebra support these cases; no failing invariant was found. | REFUTED (read-only audit) |

## Validation protocol

- All scenarios are parser-only or at most one accepted model step; no long
  trajectory was used.
- Debug and Release are separate binaries. Release uses the repository's
  fast-math flags, and NaN/Inf tests construct IEEE bit patterns directly.
- Full sanitizer, fuzz, branch-coverage, and subsystem-suite evidence will be
  appended before P9-G1 through P9-G4 can close.
