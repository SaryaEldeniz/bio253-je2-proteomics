# Step 1 — Experimental design and contrasts

## Experimental design

| Genotype | TSB | SNpulse | SN |
| --- | ---: | ---: | ---: |
| `anc` — JE2 wild type | 3 | 3 | 3 |
| `stk` — `stk1::Tn` mutant | 3 | 3 | 3 |
| `stp` — `stp1::Tn` mutant | 3 | 3 | 3 |

The experiment comprises 27 biological samples. Each sample was measured
in the total-proteome and phospho-enriched workflows.

- **TSB** is the control medium.
- **SNpulse** is growth in TSB followed by 15 min of PASN exposure.
- **SN** is approximately 15 h in 70% TSB and 30% PASN.
- **PASN** is cell-free *Pseudomonas aeruginosa* supernatant.
- Raw filenames use `pulseSN`. The sample table and contrasts use `SNpulse`.

## Contrast checklist

| Question | Contrast in the current DEA |
| --- | --- |
| What changes in the wild type after longer exposure? | `TSB_vs_SN_at_anc` |
| What changes in the wild type after the pulse? | `TSB_vs_SNpulse_at_anc` |
| How does `stk` differ from the wild type in the control medium? | `stk_vs_anc_at_TSB` |
| How does `stp` differ from the wild type in the control medium? | `stp_vs_anc_at_TSB` |
| Does the longer exposure response differ in `stk`? | `interaction_stk_vs_anc_at_TSB_vs_SN` |
| Does the longer exposure response differ in `stp`? | `interaction_stp_vs_anc_at_TSB_vs_SN` |
| Does the pulse response differ in `stk`? | `interaction_stk_vs_anc_at_TSB_vs_SNpulse` |
| Does the pulse response differ in `stp`? | `interaction_stp_vs_anc_at_TSB_vs_SNpulse` |

## Reading the signs

- For `A_vs_B`, `diff` is log2(A) minus log2(B). A negative value
  for `TSB_vs_SN_at_anc` means abundance is higher in SN.
- `interaction_stk_vs_anc_at_TSB_vs_SN` calculates
  `(stk_TSB - anc_TSB) - (stk_SN - anc_SN)` on the log scale.
  It compares the `stk`–wild-type difference between the two conditions.

## Verification

- [x] Ran `1_CheckDesign.R` and confirmed 27 unique samples.
- [x] Confirmed three biological replicates in each of the nine groups.
- [x] Checked the contrast names and the signs against the `contrasts`
      sheet in `data/raw/total/DE_total.xlsx`.

## Sources

- [Course README](../../README.md)
- Local sample sheet `data/derived/samples.tsv`
- Local workbook `data/raw/total/DE_total.xlsx`, sheet `contrasts`
