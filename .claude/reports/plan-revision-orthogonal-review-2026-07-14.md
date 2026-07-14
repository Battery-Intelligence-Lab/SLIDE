# Orthogonal review of the 2026-07-14 PLAN.md revision (breadth wave, goals 15–23)

Reviewer: independent Fable critic (Agent, adversarial brief; agentId aa8e7a15d582f6e32).
Subject: the uncommitted PLAN.md diff absorbing Volkan's 16-point directive + the integrator-research
follow-up (goals 15–23, §1.4 VC table, §3.23–3.25, D-43..D-49, ladder M1.0/M12.0/M13–M19).
Verdict: **SOUND-WITH-FIXES**. All 12 findings applied before commit; fixes listed per finding.

| # | Sev | Finding (compressed) | Fix applied |
|---|-----|----------------------|-------------|
| 1 | CRITICAL | The new §8 row recorded "orthogonal review applied; fixes in the revision commit" BEFORE the review ran and before any commit existed — a fabricated verification record, in a document whose own §0 forbids exactly that | Row rewritten AFTER the review with the true finding counts and this artifact's path |
| 2 | MAJOR | Status header (line 4) said "goals 15–22, M13–M18", contradicting §1/§6/§8 (15–23, M13–M19) — the first line every session reads | Header corrected |
| 3 | MAJOR | Goal 17 promised `device="cuda"` on EVERY batch-compiled tier while CUDA SPMe/DFN/FVM are *(opt)* boxes and no box implemented ECM/Schiffer CUDA; VC-3 bound only "claiming" tiers → M16.2 could pass vacuously | Goal 17 + M16.1/M16.2 rebound to a RECORDED GPU-tier list (minimum SPM + ECM, ECM implemented in M16.1); non-listed tiers must report the gap explicitly; M16.2 names the list so it cannot pass vacuously |
| 4 | MAJOR | VC-2 ("every fit-relevant surface through PyBOP") had no PyBOP box for ECM — the most common PyBOP target — nor a Schiffer decision | New M13.6 (ECM PyBOP fit); M14.2 must decide + record Schiffer fit-relevance; VC-2 enforcement column updated |
| 5 | MAJOR | M12.3 release notes still said "lead-acid out of scope Q16" (now overturned, D-43) and pointed "v6" at items the defined v6 ladder does not carry | M12.3 rewritten: Schiffer arrives in v6 at M14; D-36 closure/`.yp`/FMU marked unscheduled beyond v6.0.0; adjoint = M16.5 evaluation only |
| 6 | MAJOR | "PyBaMM ships no Schiffer counterpart, so VC-1 does not bind" asserted from memory with no verify marker, and conflated the ageing half with the electrical half (PyBaMM's origins include lead-acid porous-electrode models) | §3.23 + M14.2: VC-1 binding decided at M14.1 SEPARATELY per half; the no-counterpart claim marked WEB-VERIFY at M14.1 |
| 7 | MINOR | G-SHIM "within registered bands" named no compared quantity — a script that merely completes could tick it | M15.3: named `Solution` variables per script vs committed pinned-PyBaMM outputs, band per script; completion alone never ticks |
| 8 | MINOR | M1.0's "per-binary floors must not drop" is ill-defined once a binary is split (successors have no frozen floor) | M1.0: sum of successor counts ≥ predecessor floor; floors re-frozen per successor |
| 9 | MINOR | Device-token drift (§3.8 `{cpu,gpu}` vs D-45 `{"cpu","cuda"}`); M16.1 tier list drifted from goal 17 | §3.8 token unified to `{"cpu","cuda"}`; tier list reconciled under finding 3's fix |
| 10 | MINOR | §3.25(1)/M17.3 registered a years-long hypothesis while the oracle (full run) is necessarily SHORT — the tickable claim was weaker than the registered one by construction | Hypothesis re-registered at the scale the SHORT oracle can judge; long-horizon claim demoted to a separate qualified, non-blocking record |
| 11 | MINOR | Directive items 8 ("different ageing formulas") and 9 ("easy battery construction") missing from the already-covered pointer map; M18.1 implied an impossible full compiler×OS cross-product (MSVC is Windows-only) | Two pointer clauses added; M18.1 restated over APPLICABLE pairs with a committed matrix |
| 12 | MINOR | Smaller unmarked memory claims: D-43's "SimSES ships it"; VC-1's tier/counterpart names | Both marked verify-at-box |

Coverage check (by the critic, against the directive verbatim): all 16 items + the integrator follow-up
map to a contract row, section, or box; no contradiction found with D-06, D-13, D-24, non-negotiable #4,
§0.3(5), or PC-4; M12.0 and M1.0 placements coherent; §5.2/§5.6/§5.7 references match file convention.

Residual risks the critic accepted (recorded, not hidden): the M13–M19 wave sits AFTER v5.0.0, so
directive items 5/6 (shim, new tiers) do not gate v5 — sequencing chosen to avoid churning the frozen
v5 scope; and every from-memory citation in §3.23/§3.25 remains unverified until its box (M13.3, M14.1,
M17.1), by design.
