# Night-mode operating rules

Authorized 2026-07-22, modeled directly on LOFT's own `NIGHT_MODE.md`
(`hari251086/LOFT`) at the user's explicit request ("switch on night time
mode" / "refer LOFT project for correct directions"). Not a permanent
project convention — a dated authorization for continuing the OREM global
RPE investigation unattended. Safe to delete once the user is back and has
reviewed the results.

## Rule 1 — Continue the RPE investigation autonomously, track everything via GitHub issues

Keep working through the plan (`E:\claude\plans\lucky-juggling-dongarra.md`,
reachable via the `C:\Users\hari2\.claude\plans` junction) without asking for
confirmation between steps. For every finding, bug, or decision: file a
GitHub issue or comment on an existing one, not just a commit message —
issue #32 is the parent tracking issue for this investigation.

- `hari251086/OREM`, branch `HS-dev` (permanent working branch; `main` is
  the stable branch, receives merges weekly/on-demand only — do NOT merge
  to main on every commit). See [[feedback_orem_git_workflow]].
- P1-P4 priority labels, same scheme as KSROP/OREM's existing labels — see
  [[reference_gh_issue_priority_scheme]].
- No `Co-Authored-By: Claude` trailer in commit messages — see
  [[feedback_no_claude_coauthor_ksrop]] (same convention extended to OREM).
- Update README and post an issue comment on every code/tests/summary
  change — see [[feedback_update_readme_issue]].
- Never use C: drive for scratch/temp — `TMP`/`TEMP` already point to
  `E:\claude`; plan files now live on `E:\claude\plans` via the junction at
  `C:\Users\hari2\.claude\plans` — see [[feedback_no_c_drive_scratch]].
- Digging into a suspicious result (an unexpectedly clean RPE, a wildly-off
  number, a paper's claim that doesn't match OREM's own code) before
  reporting it as-is is this investigation's whole methodology — keep doing
  that rather than taking aggregate results or a single paper at face value.

**Current state (2026-07-22)**: Phase 0 literature reading is the active
task — user explicitly said to finish "reading and compiling and
remembering everything before we do trade off and implementation," so no
Phase 1+ code changes until that's done and reported. See
[[project_orem_rpe_investigation_plan]] and
[[reference_orem_reentry_literature]] for exactly what's read vs.
remaining. Once reading is complete: synthesize findings into the plan file
and memory, report back to the user, and only then move to trade-off
discussion / Phase 0's synthetic harness build — do not skip ahead on the
assumption that "more reading" and "start implementing" are the same
decision.

**Never suggest or cite Gondelach et al. 2017** — see
[[feedback_avoid_gondelach_paper]]. This exclusion survives any autonomous
continuation; don't reconsider it without the user raising it again.

## Rule 2 — Token/usage budget discipline: keep activities continuous, don't burn the session out early

- **Prefer direct tool use over spawning a subagent** whenever the task
  fits in a normal turn — a targeted `Read`/`Grep`/`Glob`, a `Bash` command.
  Only reach for a subagent when the task is genuinely broad/parallelizable
  research that would otherwise cost many sequential turns here.
- **When subagents are the right tool, keep the batch small and scoped
  tightly** — 2-3 well-scoped agents, not 5+, sequential or in small waves
  rather than one big parallel burst unless genuinely time-critical.
- **Long-running compute (test suites, campaign runs, builds) belongs in a
  background `Bash`/`PowerShell` task, not a subagent.** A backgrounded
  command costs one tool call to start and a notification to read later —
  the compute time itself is free with respect to token/usage budget.
- **Keep turns concise.** Don't re-read large files/outputs already known
  from earlier in the conversation; don't re-print unchanged status; act,
  then report the result plainly.
- **Batch related decisions instead of trickling them out** — e.g. decide
  a whole reading list or test plan up front rather than a long back-and-
  forth of single-item checks each requiring a fresh round of analysis.

## Rule 3 — Consciously avoid reaching the session usage limit

- Before launching something expensive (a batch of subagents, a large
  build/test run), consider whether it can be done more cheaply first
  (Rule 2), and whether it can be staged/sequenced rather than fired all
  at once.
- If several expensive operations are queued, spread them rather than
  bursting all in the same turn — finish and commit/checkpoint one piece of
  work before starting the next, so a limit hit mid-stream leaves a clean,
  resumable state (a closed-out issue comment, a pushed commit, an updated
  plan file) rather than a half-finished operation.
- This is a *soft, conscious* constraint, not a hard token budget to track
  number-by-number — default to the cheaper mechanism (Rule 2) and notice
  when a plan is unusually expensive before committing to it.

## Rule 4 — Resume automatically through session/usage limits

If a usage or session limit interrupts work mid-task despite Rules 2-3, the
same handling applies: wait for the limit to reset, then resume the exact
task that was interrupted, without waiting for the user to say so. Don't
restart from scratch — pick up where it left off using whatever partial
results/state already exist (the plan file's current content, an open
issue's existing comments, memory files already updated this session).
