# Reconciling the lecturer screenshot and Step 4 counts

The screenshot is reproducible from the August 2026 total-proteome export.
The project uses the October 2026 total-proteome export. Both use source
FDR <0.1 and |diff| >1 for the counts below; FDR alone gives different counts.

| Stage | stk vs anc interaction | stp vs anc interaction | stp vs stk interaction |
|---|---:|---:|---:|
| August export, all entries passing the rule (screenshot) | 120 | 104 | 130 |
| October export, all entries passing the same rule | 102 | 102 | 131 |
| October, keep ordinary CDS entries (our selected set) | 95 | 88 | 118 |
| Also require source estimate_type == observed | 85 | 78 | 102 |
| Also require >=2/3 measured values in each relevant group | 67 | 63 | 86 |

All interactions in this report compare TSB versus SN. The August export has
2,348 entries per contrast; October has 2,324 before our classification filter,
and 2,229 ordinary CDS entries after that filter.

## Dataset-version difference

The August source model names include Linear_Model_moderated and
Imputed_Mean_moderated. The current October source uses lm_impute. This comparison
establishes which file reproduces the screenshot; it does not isolate which
upstream processing/model change caused each individual statistical difference.
The extra classification filter is not the whole explanation for the count gap.

Verified sources:

- August: `data/DEAs/DEAs/DEA_20260802_WUtotal_JE2_2f_vsn/Results_WU_total_JE2_2f/DE_WUtotal_JE2_2f.xlsx`
- October: `data/raw/total/DE_total.xlsx`, represented by `data/derived/dea_total.tsv.gz`.
- Step 4: `results/step04/total/annotated_all.tsv.gz`.

The first two paths are relative to the locally configured course DATA_ROOT;
the Step 4 path is relative to the code project.

## Ordinary-protein filter

Among October rows passing the FDR/effect rule, we excluded:

| Classification | stk vs anc | stp vs anc | stp vs stk |
|---|---:|---:|---:|
| Metadata contaminants | 6 | 12 | 11 |
| Pseudogene/frame translations | 1 | 2 | 2 |
| Total excluded selected entries | 7 | 14 | 13 |

These exclusions use the supplied richer metadata classification. They are not
removals because a protein lacks a good UniProt mapping. All source rows remain
in the audit; source FDRs and effects are not recalculated.

## Conservative sensitivity subset

The conservative subset is an additional project sensitivity view, not the
lecturer's significance definition and not a newly fitted statistical test.

It requires all three:

1. Ordinary CDS, source FDR <0.1 and |diff| >1.
2. Source model label observed rather than lod_imputed.
3. At least 2 of 3 actual measurements in EACH of the four groups used by the interaction.

For stk versus ancestor, those groups are anc_TSB, anc_SN, stk_TSB and stk_SN.
For stp versus ancestor, replace stk with stp. For stp versus stk, use stk_TSB,
stk_SN, stp_TSB and stp_SN. Eight measurements concentrated in only three groups
do not meet the rule; every relevant group needs at least two.

Sequential removals from the selected set are:

| Removal | stk vs anc | stp vs anc | stp vs stk |
|---|---:|---:|---:|
| Remove source-model LOD-flagged rows first | 10 | 10 | 16 |
| Then remove remaining rows with <2/3 in a relevant group | 18 | 15 | 16 |
| Retained conservative hits | 67 | 63 | 86 |

The source flag is model-wide in this dataset. It can reflect an unobserved
group outside the displayed TSB/SN comparison. Some real, well-measured effects
are therefore excluded by this deliberately cautious source-label restriction.
Likewise, reproducible detection-versus-missing patterns can be biologically
interesting. Exclusion from the conservative subset does not mean a hit is false
or unimportant. Keep the main selected list and use this subset as a robustness
check, with actual group counts available for review.

## Reproduce and inspect

Run `source("scripts/04_analysis/04_reconcile_counts.R")` from this project's root. The
script calculates both export counts and exact filtering steps. Detailed tables
and identities of the excluded/nonconservative hits are saved locally under
`results/step04/reconciliation/`. No source analysis or dataset choice was changed.
