# AGENTS.md

Instructions for AI coding agents (Codex, Cursor, Copilot, …) are in
[CLAUDE.md](course_starter/CLAUDE.md). They apply to every agent. The most important points:

- Any AI tool may be used on the course data. **Never commit the raw data** in `data/`;
  code and rendered reports are committed and pushed to the group's own GitHub repo.
- Use **R** / tidyverse and the helpers in `R/` (`source("R/setup.R")`).
- `A_vs_B` means `diff = log2(A) - log2(B)`. Always say which group is higher.
- Plot raw per-sample values before trusting a fold change; flag `lod_imputed` estimates.
- Cite checkable sources; keep a lab notebook in `notebook/`.
