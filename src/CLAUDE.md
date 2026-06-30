# src/CLAUDE.md

Directory-specific guidance for code under `src/`.

## Scope
- Keep business logic in `src/`.
- Avoid mixing script-style side effects into importable modules.

## Implementation Guidance
- Prefer small, composable functions with clear inputs/outputs.
- Preserve existing public interfaces unless explicitly asked to change them.
- Keep preprocessing and model code deterministic where possible.

## Safety and Compatibility
- Maintain compatibility with `main.py` usage patterns.
- Do not silently change expected data formats.
- If changing assumptions, update caller code in the same task.

## Quick Validation
- Run targeted module-level checks or existing execution paths touched by edits.
- Verify that imports from `main.py` still resolve cleanly.
