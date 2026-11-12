# Work log

## 2026-11-11 — Step 1 · Design verification

### What I checked

I ran `1_CheckDesign.R` in RStudio and checked the contrast
definitions in the `contrasts` sheet of `DE_total.xlsx`.

### What I observed

- The sample sheet contains 27 unique samples, with three biological replicates in each of the nine groups.
- For `TSB_vs_SN_at_anc`, a negative value means that abundance is higher in SN than in TSB.
- An earlier automated run stopped because the course setup loaded an unavailable `dplyr` package. The saved script reads the sample sheet directly with Base R.

### AI assistance

Codex helped locate the files, draft the script and verify the saved check. I verified the contrast names and signs against the `contrasts` sheet in `data/raw/total/DE_total.xlsx`.

## Step 7 — Phospho-STYH interactions

[2026-11-12 · Analysis, table and volcano plots](2026-11-12_step07.md)
