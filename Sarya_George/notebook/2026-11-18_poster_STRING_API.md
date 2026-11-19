# 2026-11-18 — Fresh API networks and poster redesign

User requested live STRING API retrieval for the exact49 selected proteins and
the original67, then compact poster figures with central anchors, colours suited
to green colour-vision weakness, labelled pathway encirculations and overlap-aware
membership representation. Earlier figures and historical enrichment unchanged.

Read the user-linked STRING guide as a starting reference and official API docs
for current endpoint/parameter semantics. Kept the project R workflow. Obtained
turn-scoped network permission; no credentials or experimental measurements sent.

New scripts in scripts/06_network/:

1. 06_fetch_poster_networks.R — exact-input mapping/network/annotation/FASTA/links.
2. 06_fetch_string_reference_images.R — native STRING reference PNGs (evidence
   view, Leiden resolution0.7, no added colour/white nodes).
3. 06_compare_poster_annotations.R — live-versus-course membership audit.
4. 06_plot_poster_networks.R — compact, reproducible poster layout and halos.
5. 06_check_poster_outputs.R — independent source/output checks.
6. 06_write_poster_guide.R — results, provenance and reproduction guide.

Exact sets are not nested:42 shared,7 only49,25 only67,union74. All74 mapped
uniquely to reference taxon93061 and verified exact UniProt STRING crossrefs.
Curated manifests prevent five ineligible alternate JE2 rows from multiplying
the join. Used STRING12.5,functional,required_score400,add_nodes0.

Results:49nodes204edges10isolates (104edges>=.7);67nodes390edges4isolates
(234edges>=.7). Fresh49 reproduces all204 original scores exactly. Functional
annotations contain1082/1235 term rows with inputGenes lists; every gene scope
and term count reconciled. Reference FASTAs match exact49/67 IDs. All retrieval
records HTTP200 on18Nov2026; cached reruns preserve those retrieval times.
RawTSVs,FASTAs,nativePNGs and link records archived under ignored
data/raw/step06_string/poster_api/; public-resource-export ZIP prepared locally.

Important annotation finding: current STRING KEGG purine membership10/24 differs
from course scope24/33 (49/67). Other3pathway memberships agree exactly. Asked
student optional preference via asynchronous question, continued independent
work, then after reasonable opportunity stated default to fresh KEGG groups.
The student subsequently explicitly confirmed fresh STRING KEGG assignments;
the completed live-mode figures already match that choice.
All proteins remain;14/9 not assigned to these four fresh terms are shown as
additional selected context. No historical GSEA membership/NES/FDR recalculated.
Missing current annotation is not proof of an incorrect original assignment.
There are no actual multi-memberships across the four selected terms; graphics
support overlaps and multiple badges without fabricating them or duplicating nodes.

Poster output: exact49/exact67 each PNG300dpi,vectorPDF/SVG. Stk1/Stp1 pinned
symmetrically centrally for prominence, not network-derived centrality. Neutral
node centres; blue/orange/charcoal/magenta outlines + distinct line styles,
P/R/N/T badges and neutralCcontext. Ring orderinnerstk/anc,middlestp/anc,
outerstp/stk; original signed effects unchanged, common+/-2 colour scale with
display-only clipping, NAgrey, sourceFDR<.1 dots. No green-vs-red identification
required for pathway membership. Palette is supplemented by labels/styles;
no claim to simulate every individual's colour perception.

Compact canvas:49=16x9.86in;67=16x10.52in. Footer11.2%/10.5%. Isolates
packed with groups, not a distant extra row. Minimum centre distance1.72 vs
node diameter1.42. Neither anchor lies in a pathway/context outline; minimum
anchor-ring-to-outline clearance1.15/0.59. Coordinates and boundaries exported.

Development checks corrected native-pipe syntax and summary variable shadowing
before completion. First plot pass exposed locale glyph artifacts and clipped
footer; switched plot-visible labels to ASCII and two short footer lines. Root
visually inspected both final live figures; plotting agent inspected and checked
geometry. An independent read-only agent audited exact sets, API source/mappings,
all effects/FDRs, annotation gene counts and the course/live discrepancy. R checks
reconcile inputs, edges, all source values, halos and six nonempty poster exports.

AI assistance: Codex wrote retrieval/comparison/review scripts and guide; one
agent independently audited the data and another implemented/verified layout.
