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
4. **Write a scalar-generic physics body.** Keep the equation in its named mechanism header and separate from lane/mask scaffolding. Support the scalar types needed by forward sensitivities, retain the selected scalar-generic operand in piecewise branches, and test both smooth sides against finite differences. Reuse `SpmObservables`; never recompute a private, inconsistent voltage or stress path.
5. **Use the common hot-path-allocation-free scaffold.** Give the body a named `Basic...Output<Real>` view and a construction-time scratch owner backed by `detail::AgeingScratchStorage`; scratch allocation happens only during construction. Clear accumulated or multi-model fields with `detail::clear_ageing_fields` before traversal; a field that is unconditionally overwritten once per lane need not be cleared. Use `detail::for_each_ageing_lane_while_success` and, for numbered alternatives, `detail::for_each_enabled_ageing_model_lane`. The numerical order is part of the contract: models and lanes are both traversed in ascending order.
6. **Compose it transactionally in `SpmPipeline`.** Place it after the shared observable/stress stage and invoke it through `detail::evaluate_ageing_stage`: disabled is a no-op, a failed compute returns without applying scratch, and only a successful compute adds to the already-zeroed derivative. Add its layout and parameters to `SpmPipelineParams` and include it in the exact restart/rollback path.
7. **Compile it in `SpmFactory`.** Add a cold option/mask, validate incompatible combinations, and map Python/MATLAB option names to that one registry selection.
8. **Protect every bridge.** Update recording for any new observable, dual/sensitivity propagation, and CUDA only when the mechanism is genuinely supported; otherwise make device preflight reject it explicitly.

Required tests are: direct nonzero rate-law points against an independent calculation, both current signs where relevant, every smooth piecewise branch checked by finite differences for supported scalar types, invalid parameters and states, disabled-stage zero effect, failed-stage no-publication through the full pipeline, additive composition with another mechanism, heterogeneous/non-SIMD lane counts, batch-versus-isolated lane equality, A/B/A scratch reset, zero accepted-step allocations, and bitwise checkpoint/restart. Update the 9C-2 structural gate when adding a field, model mask, or stage so field mappings and pipeline order cannot drift on toolchains where recorded bits are disabled. A lower voltage RMSE does not validate an ageing mechanism if lithium inventory, active fraction, surface area, or local concentration becomes unphysical.

Adding only an `AgingMechanismSpec` entry is insufficient; that type is cold metadata and does not connect a rate law to the production pipeline.
