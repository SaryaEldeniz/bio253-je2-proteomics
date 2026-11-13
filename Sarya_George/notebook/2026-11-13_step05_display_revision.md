# 2026-11-13: observed-only complete enrichment figures

## Request

Remove SNpulse from the two sample-evidence plots; show every significant GO/KEGG
ORA/GSEA result excluding source LOD-imputed entries; sort GSEA rows by descending
signed NES and ORA rows by descending fold enrichment. Explain GO_BP and the
GO/KEGG versus ORA/GSEA distinction.

## Analysis policy

Added `observed_only` rather than substituting the earlier `observed_replicated`
policy: that earlier sensitivity would also remove observed proteins with only
1/3 measured values in a relevant group. The requested universe is the same
high-quality singleton UniProt mapping as before, with finite source statistics
and `estimate_type == "observed"`. It contains 1,858 proteins in each comparison;
200 source-LOD entries are excluded relative to the 2,058-protein primary universe.
Sparse observed proteins remain eligible and their support is auditable.

The source protein effects, FDRs and ranking statistics are unchanged. ORA lists,
GSEA ranks and the eligible background all use the observed-only universe before
recalculating enrichment. Kept the same size10-500 sets, protein hit thresholds,
seeded multilevel fgsea engine and study-wide BH family across ten comparisons
and five original collections, separately by method. UniProt keywords remain
in that correction family; the requested displays include KEGG and GO BP/MF/CC.

Completed 15,720 directional ORA tests and 7,860 GSEA tests. Across the four
displayed GO/KEGG collections there are 158 significant ORA term-comparisons
and 1,374 significant GSEA term-comparisons (study FDR<0.1). Shared memberships
mean these are not counts of independent biological modules.

## Figures and code

The two existing sample PNG/PDF pairs keep the same proteins and now use only
TSB and SN with compact positions. These plots use measured sample values;
their original primary-based protein selection is separate from the new
observed-only enrichment policy.

The new all-results renderer writes one complete PNG/PDF/orderedTSV per
method/collection/comparison: 80 figure pairs, including explicit empty results.
There is no Top-N display cap. GSEA is ordered largest positive NES to most
negative; ORA directions remain separate and each sorts highest fold enrichment
to lowest. Long plots expand vertically and vector PDFs support zooming. The
original primary Top-N figures remain a distinct initial view, clearly linked
from the updated guide to the new observed-only index.

New entry points: `run_project_step("enrichment_observed")` and
`run_project_step("enrichment_all_figures")`. `read_enrichment()` accepts
`policy="observed_only"`. All analysis and figure code is saved in the group
folders; raw data and full caches remain ignored.

## Verification

The original testthat suite and seven ORA checks passed. A separate result
validator independently reconstructs observed-only eligibility from the source,
checks no extra replicate filter, verifies the source statistics, foregrounds,
backgrounds, overlaps and leading edges, independently recalculates both BH
correction levels, and checks a real ORA result against one-sided Fisher testing.
Its final coverage checks compare every plotted TSV with every significant
GO/KEGG test, check all80 outputs and verify descending panel row order.
See reports/step05/observed_only/validation.txt for completed checks.
Representative long/short/empty enrichment plots and both revised sample PNGs
were visually checked; this is not a claim of manual inspection of all80plots.

## Interpretation

BP=Biological Process; MF=Molecular Function; CC=Cellular Component. GO defines
linked functional/process/location categories; KEGG defines pathway maps.
ORA tests concentration among selected hits; GSEA tests concentration toward
the ends of the full ranking. Enrichment does not directly measure pathway
activity or establish Stk1/Stp1 phosphorylation targets.

## AI assistance

Codex implemented the observed-only analysis and documentation. Separate agents
handled the two-media sample plots, the uncapped sorted figure renderer, and
independent eligibility/statistics/display validation. The user specified the
exclusion, display scope and ordering. No new source differential models were fit.
