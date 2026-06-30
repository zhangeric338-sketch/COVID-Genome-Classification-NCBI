# tools/CLAUDE.md

Directory-specific guidance for scripts under `tools/`.

## Scope
- Keep scripts task-focused (download, utility, maintenance).
- Prefer idempotent behavior where practical.

## Download and I/O Rules
- Default dataset size: `0.05 GB`.
- Never enable full dataset unless explicitly requested.
- Respect API rate limits and retry policies.
- Check for existing outputs before writing new files.
- Validate outputs are non-empty and structurally valid.

## Script Design
- Keep CLI options explicit and documented.
- Fail with clear error messages when prerequisites are missing.
- Avoid hidden side effects outside expected output directories.

## Quick Validation
- Exercise changed CLI paths with minimal-size inputs.
- Confirm command help text remains accurate after argument changes.
