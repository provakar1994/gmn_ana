# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Workflow rules

### 1. Plan first, get consent, then act
- Do **not** make substantial changes to any file without first presenting a plan and getting explicit approval.
- Break every task into small, self-contained steps. Execute them **one at a time**, confirming with the user before moving to the next step.
- Small, obviously-requested edits (typo fixes, a one-line change the user asked for) can proceed directly. If in doubt, ask.

### 2. Commit regularly, but only after review
- Commit after each completed step, so every change is small and easy to review or revert.
- **Before every commit, show the user the changes (diff summary + proposed commit message) and wait for their review and explicit approval.** Never commit without it.
- Write clear, descriptive commit messages explaining *what* changed and *why*.

### 3. Keep a change log
- Record every change in [`CHANGELOG.md`](CHANGELOG.md) at the repository root (create it if it doesn't exist).
- One entry per change, newest at the bottom, in this format:

  ```markdown
  ## YYYY-MM-DD HH:MM
  - **Files:** `path/to/file.cpp`, ...
  - **Change:** short description of what changed and why.
  ```
- Update the log in the same commit as the change it describes.

## Coding style
- Keep code as **simple and human-readable** as possible. Prefer clarity over cleverness.
- Add **ample comments**: explain the purpose of each function/block, the physics meaning of cuts and variables, and any non-obvious logic.
- Follow coding best practices: descriptive names, small focused functions, no magic numbers (use named constants or config entries), no duplicated code, and consistency with the style of surrounding code.

## Repository layout
- `gmn_ana/` is the repository base; all paths below are relative to it.
- `include/`, `src/` — shared library code (headers and implementations).
- `scripts/` — scripts that analyze **real data**, and scripts that do **data/MC comparisons**.
- `scripts/simulation/` — scripts that analyze **simulated (MC) data**.
- `tests/` — unit tests for the library (`test_*.cpp`, run with `ctest`); small test inputs in `tests/data/`. Add or update a test with every library change that affects physics results.

Place new scripts in the matching directory above.

## Coordination with the MC consistency audit (../MC)
Data/MC consistency work is coordinated from the MC project, which owns the decisions.
- Task specs for this repo: ../MC/handoff/  (start here: read the spec of the task you're on)
- Open items and findings: ../MC/CONSISTENCY_TODO.md  (items M2, M5, ...)
- Decisions, accepted discrepancies, evidence: ../MC/AUDIT_REPORT.md
Rules:
- Treat "Decisions already made" in a spec as fixed; raise disagreements with the user instead of changing them.
- When a handed-off task is done, update only its status line in ../MC/CONSISTENCY_TODO.md
  (gmn_ana commit hash + one line) and add a dated Findings entry with the requested results.
  Do not edit other parts of ../MC files.
- Log changes in this repo's CHANGELOG.md as usual.
