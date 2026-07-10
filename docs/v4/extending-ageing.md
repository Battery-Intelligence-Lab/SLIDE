---
layout: default
title: Add an ageing mechanism
nav_order: 7
---

# Add a v4 ageing mechanism

Ageing mechanisms are compile-time stages inside an SPM batch. They may be selected once per archetype through cold options, but they are not per-cell virtual objects.

1. **Derive the mechanism.** State the rate law, sign convention, units, required observables, conserved inventory, valid parameter domain, and a result that would falsify the implementation.
2. **Create a named parameter pack.** Follow `SeiParams`, `SurfaceCrackParams`, `LamParams`, or `LithiumPlatingParams`. Validate all non-finite and unphysical values before replacing a valid batch.
3. **Declare every persistent state.** Extend `SpmStateLayout`/`declareSpmState`, which `SpmPipeline::declareLayout` invokes for the shared SPM arena. Hidden “previous” values are forbidden because arena snapshots must be complete.
4. **Write a scalar-generic physics body.** Keep the equation separate from lane/mask scaffolding and support the scalar types needed by forward sensitivities. Reuse `SpmObservables`; never recompute a private, inconsistent voltage or stress path.
5. **Add an allocation-free lane stage.** Preallocate scratch at construction, skip a disabled mask before entering lane loops, and add contributions to the already-zeroed derivative. Do not overwrite another mechanism's contribution.
6. **Compose it in `SpmPipeline`.** Place it after the shared observable/stress stage, add its layout and parameters to `SpmPipelineParams`, and include it in the exact restart/rollback path.
7. **Compile it in `SpmFactory`.** Add a cold option/mask, validate incompatible combinations, and map Python/MATLAB option names to that one registry selection.
8. **Protect every bridge.** Update recording for any new observable, dual/sensitivity propagation, and CUDA only when the mechanism is genuinely supported; otherwise make device preflight reject it explicitly.

Required tests are: direct rate-law points against an independent calculation, both current signs where relevant, invalid parameters and states, disabled-stage zero effect, additive composition with another mechanism, heterogeneous/non-SIMD lane counts, zero accepted-step allocations, and bitwise checkpoint/restart. A lower voltage RMSE does not validate an ageing mechanism if lithium inventory, active fraction, surface area, or local concentration becomes unphysical.

Adding only an `AgingMechanismSpec` entry is insufficient; that type is cold metadata and does not connect a rate law to the production pipeline.
