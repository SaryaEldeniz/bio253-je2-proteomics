# 2026-11-17 — Step 6 parts 5–7: STRING networks and contrast halos

Question: how do the 47 user-selected proteins and Stk1/Stp1 connect in STRING,
which additional partners are useful, and how do their three protein interactions
compare in a concentric-halo network?

The student supplied seven STRING files. An asynchronous question confirmed a
functional network with >=0.4 confidence for the MAIN analysis. Kept >=0.7 as a
supplementary view, without changing the student's requested main threshold.

Ran from the group R project:

```r
source("scripts/06_network/06_import_string_exports.R")
source("scripts/06_network/06_expand_string_network.R")
source("scripts/06_network/06_review_network_proteins.R")
source("scripts/06_network/06_plot_string_networks.R")
source("scripts/06_network/06_write_network_guide.R")
```

Import: all 49 UniProt accessions matched exact alias tokens uniquely to taxon
93061 STRING IDs. All coordinate/degree/FASTA node sets agree. Recomputed degrees
match the degree export. Corrected swapped description/SMART-URL headings in the
protein annotation TSV. Saved raw downloads under ignored project data/raw/.
The 3,725 functional-annotation rows are memberships, not enrichment results.

Baseline: 49 nodes, 204 unique undirected edges, 10 isolates at >=0.4; 104 edges
and 14 isolates at >=0.7. Live version-specific STRING API reproduced every seed
edge and score exactly (zero differences). Version is based on student context
and explicit 12.5 API address; XML version2.5 is the PSI-MI schema.

Expansion: requested up to5 partners per seed; 242 rows (2–5 each) yielded126
distinct external candidates. Queried the175-node induced discovered network:
1,289 edges. Resolved via exact UniProt STRING cross-references:115 unique JE2
maps,6 exact reference maps without JE2 match,4 ambiguous,1 unresolved. Unmapped
measurement status remains unknown, not falsely unmeasured. The first shell
request failed DNS; obtained turn-scoped network permission and succeeded.
API responses and query provenance were cached; no experimental effects uploaded.

Compact selection: first preserve a contextual partner in each selected module,
then fill the12-node display budget by anchor/bridge connections, measurement
availability and confidence. No NES, effect magnitude or protein FDR entered this
selection. All126 candidates remain in full tables with decisions. Fixed a
join-suffix collision during development before results were produced.

Final compact network: all49 original nodes +12 partners =61 nodes,363 edges.
176 edges meet>=0.7;187 are .4–.7. Added PheT,RpsG,B7H15_06340,RibBA,Ldh2,Ldh1,
GlmU,GatA,Fmt,GuaA,AtpC,B7H15_06820; all12 uniquely mapped and measured. GuaA/PheT
connect26/25 original seeds; this is selected-network degree, not global hubness.

Three ring colours: innerstk/anc,middlestp/anc,outerstp/stk; shared signed-effect
scale±2 with display-only clipping, greyNA distinct fromzero, FDR<.1 dots. Constant
ring widths, module centre colours, dark anchors, dashed added centres, daggers
for retained annotation review. Original values preserved. Source-observed and
>=2/3 support required for quantitative halos; deletion-affected anchors can be
grey. Eight network views inPNG/PDF/SVG: original49,strict49,expanded61,anchor
detail,four module zooms. All nodes including isolates retained in completeviews;
zoom membership manifest saved. Deterministic base-R force layout required no
extra packages; an initial stretched layout was fixed after visual inspection.

Experimental review: original47 candidates include9 large/38 smaller supported
effects; seven existing membership flags retained. Among12 added proteins,Ldh1
has stk/anc−1.158,FDR.000322,stp/stk+.758,FDR.00904; Ldh2 and B7H15_06340 have
smaller supported effects. Nine additions are context only under current rule.
Ldh1 SN−TSB responses anc+2.048,stk+3.207,stp+2.449: increase strongest instk.
Generated measured normalized TSB/SN sample panels for all12 (two pages). RibBA
has2/3 instk_TSB and appreciable spread; no measurements removed or model refit.

Evidence caution: Stk1–Stp1 .998 and Stp1–Ndk .796 are original. New high-score
Fmt/B7H15_06340 anchor links are strongly genomic-neighbourhood/coexpression
supported. Ldh1/2–Stp1 .442 links have experimental-channel support but moderate
combined confidence. GlmU–Stk1 .513 is text-mining-only. No direct-substrate or
physical-binding claims. Three interaction effects are algebraically dependent;
stp reducedStk1/polar-effect caveat remains. Protein abundance is not flux or
phosphorylation. No source FDRs/estimates changed; no new enrichment test run.

Checks: independent agent audited mappings, all49 seed preservation, all204
edge scores, every61×3 original effect/FDR, every finite halo, context ranking,
four module coverage slots, and the exact12 selected added IDs. Raw/cache data
are git-ignored. Final root visual checks inspected expanded and anchor PNGs
and added sample page1; plotting agent inspected all network views. No push.

AI assistance: Codex wrote/ran R code and report; one read-only agent independently
audited evidence and one agent implemented/visually reviewed halo plots. Next:
student inspects module zooms/evidence tables, chooses hypotheses and compares
with phosphosite data before asserting Stk1/Stp1 targeting.
