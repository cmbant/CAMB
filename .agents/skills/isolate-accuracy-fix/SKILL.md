---
name: isolate-accuracy-fix
description: "Isolate minimal efficient code changes to improve numerical accuracy or fix regression mismatches in CAMB. Use for run_tests.py failures, check_accuracy instability, step-size or switch sensitivity, low-k or recombination-era regressions, and when accuracy fixes must preserve runtime."
argument-hint: "Describe the failing test, ini file, output mismatch, or accuracy-sensitive regime to isolate."
user-invocable: true
---

# Isolate Minimal Accuracy Fixes

Use this skill when CAMB has a numerical mismatch, stability failure, or accuracy regression and the goal is to find the smallest code change that fixes it without paying unnecessary runtime cost.

This workflow is for targeted isolation, not broad tuning. Prefer a local fix to a global increase in tolerances, sampling, or step control.

## When to Use

- `fortran/tests/run_tests.py` fails on one or a few parameter files.
- `camb/check_accuracy.py` shows unstable outputs or larger-than-expected boosted-reference differences.
- A change in step size, switch timing, approximation handoff, or tolerance appears to affect results.
- A proposed fix improves accuracy but may slow unrelated runs.
- You need to identify which switch, redshift range, `k` range, or output family is actually sensitive.

## Main Goal

Find the narrowest code path and the narrowest phase-space region that controls the failing behavior, then implement the smallest performant fix that removes the mismatch.

## Core Rules

1. Start from a concrete failing anchor.
   Use a failing ini file, output stem, test command, or changed routine.

2. Reproduce in a clean comparison environment.
   If the workspace is dirty, use a disposable copy under `/tmp` so A/B tests do not disturb user edits.

3. Isolate the list of physically distinct failing outputs and address them one at a time.
   E.g. low-l EE, mpk, and high-L TT failures are likely to be different issues. Focus separately on each in turn.

4. Do not infer causality from a command option just because it is present.
   Vary or remove the suspected trigger and check whether the same physical output still fails. E.g. a command with high lensing accuracy can still expose an unrelated low-l EE issue.

5. a) Probe broadly first. Sweep AccuracyBoost / lAccuracyBoost / lSampleBoost
   to confirm the issue is accuracy-controlled, not a bug.
   The camb/check_accuracy.py script command can do this for stability under changes in accuracy paremeters, as well as providing info on which other less broad boost parameters affect the result. You may need to wait many minutes while it runs.
   b) Then change one lever at a time when narrowing the cause, more locally that just a boost parameter.

6. Prefer local overrides over global defaults.
   If a change is only needed in a specific regime (a redshift window, k-range, or output type), test a localized change before touching global defaults. A broad change is acceptable when the precision gain is broad and the runtime cost is negligible; narrow the target when the fix is expensive.

7. Treat runtime as a first-class constraint.
   A fix is incomplete if it passes but slows transfer or matter-power runs for no precision benefit.

8. Your changes should make physical sense.
   Target changes to relevant physical time/k/other scales, rather than be ad-hoc patches that treat symptom rather than accuracy of general underlying calculation. Give physical explanation for each physically distinct change that you make, and if it makes no sense, investigate further.

9. When changing a sampling grid, identify all consumers of that grid.
   CMB source grids, transfer grids, and lensing grids can be shared with matter-power or source-window paths. Run an adjacent shared-path check before accepting the change.

10. Treat index-based regression tests as grid-fragile.
    If a test checks `array[i]`, decide whether `i` represents a physical point or just the old sampling. For grid changes, validate by interpolating to fixed physical coordinates before updating expected values.

11. Separate "more accurate" from "different grid".
    A passing boosted comparison is not enough if unrelated outputs move. Check a like-for-like physical comparison or a higher-accuracy reference before updating tests.

12. Prefer direct control variables over proxy boosts.
    If a broad boost works, identify which concrete derived quantity changed
    (e.g. `dlnk`, `dk`, `lognum`, sample number, switch location, integration limit). Prefer changing that quantity directly rather than keeping a boost that also
    changes unrelated ranges.

## Acceptance Matrix

For each proposed local fix, test:

- The original failing command, suppressing unrelated output families only when needed to keep focus.
- The same physical output with the suspected trigger removed or varied.
- An adjacent path that shares the changed grid, switch, or tolerance.
- The default unit/regression test covering that shared path.
- Runtime before and after for the target run.
- For any accepted boost-derived fix, document the concrete internal quantity
  being changed and confirm unrelated quantities controlled by the same boost
  are not part of the final patch.
- Never use a local scaled boost fix where a more specific/physically-targeted fix to specific k/time/output/etc might work available.
- No temporary or unjustified code changes remain in the working tree

## Deliverable

Leave the final, minimal fix applied to the working tree as unstaged changes so the user can review the diff directly.

Clean up any A/B scaffolding before finishing where not likely to be useful again:

- Remove disposable `/tmp` copies used for comparison.
- Revert any exploratory edits, debug prints, or temporary parameter bumps that are not part of the final fix.

Report should include (for each separate physical issue):

- The minimal failing reproducer (ini + command).
- The localized change(s) (file + lines) that resolves it
- Before/after numbers for the failing diagnostic.
- Runtime impact on a representative non-failing run (cpu times).
- Physical justification for why the change targets the right scale/domain.
