---
name: isolate-accuracy-fix
description: "Isolate minimal, physically-scoped accuracy fixes in CAMB. Triages distinct failures, addresses the one the user asked for (or asks which, if unspecified and several exist), delegates diagnosis to an executing sensitivity probe, runs the acceptance matrix on rebuilt code, then delegates a final audit to a read-only checker, and writes a changelog entry per fix. Retains sole patch authority. Use for run_tests.py failures, check_accuracy instability, step-size/switch sensitivity, low-k or recombination-era regressions, and accuracy fixes that must preserve runtime."
argument-hint: "Describe the failing test, ini file, output mismatch, or accuracy-sensitive regime to isolate."
user-invocable: true
---

# Isolate Minimal Accuracy Fixes (Orchestrator)

You own the failing anchor, the working tree, and the final report.
You are the ONLY actor allowed to edit the working tree.
You do NOT run broad sweeps yourself (delegate to the probe). You DO run the
acceptance matrix yourself, on rebuilt code. You do NOT grade your own reasoning
(delegate the final audit to the checker).

Work strictly ONE physically distinct issue at a time, each as its own unstaged
diff plus its own changelog entry. Never combine fixes for different issues.

## Disposable tracked copy

For every `/tmp` reproduction or A/B copy, use this skill's
`isolate-accuracy-fix/scripts/make_tmp_tracked_copy.py`. It creates a detached git worktree at HEAD,
initializes submodules, applies staged tracked changes, then applies unstaged
tracked changes, including tracked submodule diffs, then builds CAMB in the
copy. It intentionally ignores untracked files. Record the printed copy path,
source HEAD, and build summary in the evidence bundle.

## Step 1 — Triage (you)

Reproduce the failing command once in a helper-created `/tmp` copy unless the
command is demonstrably read-only and emits no build or output artifacts. Group
failures by LIKELY DISTINCT PHYSICAL CAUSE, not by file or line.

- Distinct (non-exhaustive examples): low-l EE vs mpk vs high-L TT vs
  recombination-era TT vs lensing, etc.
- NOT distinct: many adjacent l-bins of the same spectrum failing together.

Emit a candidate list, one row per distinct issue:

  id | physical_output | l/k/z region | one-line evidence | likely_shared_cause?(y/n)

Triage is a HYPOTHESIS. The probe may collapse two apparent issues into one
controlling quantity, or split one. Report that back rather than assuming.

## Step 2 — Choose (HARD GATE)

- If the user's invocation already named the aspect/output/test to focus on:
  select that issue and proceed. Do NOT ask.
- Else if exactly one distinct issue: state the selected issue and proceed.
- Else (several distinct issues AND the prompt did not disambiguate):
  PRESENT THE LIST AND ASK THE USER WHICH ONE TO ADDRESS.
  Do NOT auto-pick. Do NOT begin a fix for more than one.

## Step 3 — Diagnose (delegate to `sensitivity-probe` subagent)

Spawn and wait for the probe for the ONE chosen issue. Pass it the ini, exact
command, the single failing output, and the suspected trigger. The probe may need many minutes for rebuilds,
sensitivity runs, or `check_accuracy.py`; wait for its report rather than
timing out early.

Require back the full Sensitivity Report defined in that agent. It must end the
progressive drill-down at CODE LINES, not at a boost:
- broad_axis: the first boost in the fixed grid that restores the output;
- sub_parameter: the specific sub-accuracy parameter under that axis;
- functions: where that sub-parameter enters the calculation;
- controlling_quantity: concrete var + file:line (the drill-down endpoint);
- coordinate_locality_table: whether the failure localizes in available
  physical/numerical coordinates, with code vars and measured evidence;
- candidate_physical_levers: narrower direct levers tested or rejected, with
  accuracy and timing evidence;
- full_physical_extent: an EVIDENCE TABLE, one row per co-affected output
  family, each with a code link AND a measured or shared-grid justification.

If accuracy_controlled is false, STOP and report a likely bug. Do not patch.

## Step 4 — Implement (you)

Change ONE direct control variable named in the report, scoped to the FULL
physical extent the probe evidenced — no narrower (do not gate on a flag that
excludes a co-affected path).

Hard locality gate:
- Reject a constant floor, broad boost, or component-wide multiplier unless the
  probe's locality audit shows why narrower physical/numerical coordinates are
  unavailable, ineffective, unstable, or more expensive than the broad change.
- If the affected code is evaluated over a coordinate such as k, l, time,
  redshift, mass, sample/momentum index, or output family, either gate/scale the
  change in that coordinate or include measured evidence explaining why not.
- Do not use a proxy gate when a more direct computed quantity is available and
  untested.

On breadth, balance scope against COST:
- Prefer a local override to a global default WHEN the broad change carries real
  runtime cost.
- If the final change has non-negligible runtime cost, test smaller floors and
  coordinate-local gates before accepting the chosen value/scope. Record the
  accuracy and timing reason they were rejected.
- BUT if a slightly broader / more conservative change (e.g. a bit more accuracy
  than the minimum that passes) has NEGLIGIBLE numerical cost, prefer the safe
  conservative choice over a brittle just-barely-passing threshold.
- Narrow aggressively only where the fix is expensive.

Give a physical reason. If there is none, return to Step 3.

## Step 5 — Acceptance matrix (you run it, on rebuilt code)

REBUILD the working tree first, then run ALL rows against the rebuilt binary —
never against a stale build. Record results; the checker cannot run them:
- original failing command now passes;
- same physical output with the suspected trigger removed/varied;
- adjacent path sharing the changed grid/switch/tolerance;
- EVERY co-affected output family the probe flagged (under-scope check);
- continuity/edge checks for every new or decisive hard gate, branch predicate,
  or threshold in the final fix. For any gate over an otherwise continuous
  physical or numerical control, test representative points on both sides of the
  boundary, and record why any remaining discontinuity is physically/code-path
  justified. Prefer physically continuous accuracy changes where sensible
  (with k, q, time, momentum, etc.).
- the default regression test covering the shared path;
- cpu time before/after on a representative NON-failing run, DISCARDING the
  first CAMB call;
- broad-vs-local comparison: if the final patch uses a broad boost/floor,
  include the Step 4 narrower-lever evidence (accuracy + timing), or note why no
  such lever exists in the affected code path;
- reference independence: the final default must not be validated only against a
  reference sharing the same boosted controlling quantity. Include an independent
  convergence check, e.g. a still-higher value of the same control, a different
  physically justified reference, or a fixed-coordinate convergence sequence;
- like-for-like check: a passing BOOSTED comparison is not sufficient —
  interpolate to fixed physical coords before blessing index-based values.

If any row fails, return to Step 4 (or Step 3) before continuing.

## Step 6 — Final Checker Review (delegate to `accuracy-fix-checker` subagent)

Once the acceptance matrix passes for the current candidate fix, spawn and wait
for one checker round. It runs nothing; it audits reasoning and the evidence you
already produced. Send the full bundle so it does NOT redo discovery:

- the proposed final diff (or key changed files);
- the physical mechanism explanation;
- the failing reproducer and before/after numbers;
- the Step 5 acceptance-matrix results (confirming tests ran on a rebuild);
- reference-convergence evidence (boosted/higher-accuracy comparison);
- the timing summary (first CAMB call discarded);
- the probe's complete Sensitivity Report (esp. full_physical_extent table);
- the coordinate_locality_table and candidate_physical_levers, including any
  narrower alternatives rejected before accepting a broad boost/floor;
- known residual risks or deliberately ignored failures.

Act on the verdict:
- PASS: proceed.
- PASS_WITH_NOTES: proceed only when the notes are clearly optional extensions,
  residual caveats, or documentation notes that do not require changing the code,
  scope, evidence bundle, acceptance matrix, or final claim. Record them in the
  changelog/final report.
- BLOCK, or any PASS_WITH_NOTES item that actually requires refinement: continue
  the fix loop. Fix the issue or gather the missing evidence LOCALLY yourself,
  rebuild, re-run the affected acceptance-matrix rows, then run one new checker
  round on the updated evidence. Do not finalize after a BLOCK without a
  subsequent PASS or acceptable PASS_WITH_NOTES. Only downgrade the claim instead
  when the checker identified a clear extension outside this selected physical
  issue.

## Step 7 — Finalize this issue (you)

- Leave ONLY this issue's minimal fix and its changelog as unstaged changes.
  Nothing else from another issue or from exploration may remain in the tree.
- Write docs/changelog/<YYYY-MM-DD>-<slug>.md using the template below.
- Clean up: delete /tmp build copies, revert debug prints and temporary boosts.

Changelog template:

  # <slug>: <one-line summary>
  - Date / issue id
  - Symptom & minimal reproducer (ini + command)
  - Root cause: drill-down broad_axis -> sub_parameter -> function ->
    controlling quantity (var + file:line)
  - Full physical extent: output families / k / z / l affected, with the evidence
    that the fix scope matches it (not too narrow; not broader than cost justifies)
  - Locality / rejected alternatives: narrower coordinates or physical levers
    tested, with accuracy/timing reason for accepting the final scope
  - Change: file + lines
  - Before/after numbers for the failing diagnostic
  - Acceptance-matrix summary (note tests ran on a rebuild)
  - Runtime impact (cpu before/after, first call discarded) on a non-failing run
  - Physical justification for the targeted scale/domain
  - Checker verdict (one line)

## Step 8 — Next issue (you)

If untreated distinct issues remain, surface the list again and ASK whether to
continue. Each accepted issue starts a fresh probe -> implement -> acceptance
matrix -> checker -> diff -> changelog cycle. Never carry an in-progress change
between issues.

## Hard rules

Cross-cutting invariants not fully covered by a single step above (the breadth /
locality / scope doctrine lives in Step 4; do not restate it here):

1. Reproduce A/B only in a helper-created clean /tmp copy of the CODE (not
   output/plot/data directories); never disturb user edits.
2. Always rebuild before running any test or timing; never test a stale binary.
3. Do not infer causality from a command option just because it is present.
4. Treat index-based regression tests as grid-fragile.
5. Runtime is a first-class constraint; discard the first CAMB call when timing.
