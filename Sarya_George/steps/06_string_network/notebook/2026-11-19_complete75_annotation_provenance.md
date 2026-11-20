# Restore eight candidates and clarify the origin of annotation questions

The student asked for the source of the questioned assignments and requested
all eight omitted proteins in the poster network. Clarified that the review
flags were authored by Codex after comparing supplied pathway memberships with
JE2/reference descriptions. They were not warnings from the lecturer or STRING.
The earliest list appears in 05_extract_module_drivers.R; step6 reads the
expanded module_annotation_review.csv. The quantitative criteria were met by
each of the eight in at least one interaction; the additional manual annotation
filter alone excluded them from the original67. This was too restrictive for
exploratory network analysis. No claim that a course annotation is wrong follows
from an unfamiliar product name or an absent current reference assignment.

06_prepare_complete75.R restores Q2G010,Q2FYK7,Q2FYG2,Q2FUT0,Q2FV52,Q2FXM9,
Q2FWF8,Q2G0W1 to the67. Exact75 protein_keys match sample_plot_index.csv.
Fresh75 mapping/API network/functional annotation/reference FASTA/link exports
were retrieved from STRING12.5 (functional,score400,no added nodes) and UniProt.
All75 mappings match exact UniProt STRING cross-references. Earlier49/67 requests
reuse cached responses with their original retrieval timestamps.

75nodes,412associations>=.4,238>=.7,10isolates. Isolates retained. Current four
KEGG groups contain P24,R8,N7,T17,plus17context and2anchors. All eight additions
fall outside these four current terms; other KEGG terms, if returned, are listed
in restored_eight_annotation_provenance.csv. This later source disagreement is
separate from the earlier manual judgment. Course terms and review notes remain
auditable in course fields. Historical ORA/GSEA and primary selection tables
are unchanged; no phosphosite overlap analysis performed.

06_plot_poster_networks.R now supports JE2_POSTER_SETS=75 and merges geometry
records while preserving the previous49/67 outputs. NewPNG4800x3331,16x11.10in
at300dpi plusPDF/SVG. All eight restored in context, centralStk1/Stp1,three
original interaction halos,accessible outline colours/styles. No node overlap,
clipping,or anchor containment. Mincentre1.72>diameter1.42;clearance0.59.
Root and plot agent visually inspected PNG. Native STRING75 image also archived
with HTTP200. Raw responses stay ignored.

06_write_complete75_guide.R writes the exact eight-protein audit with manual
CSV line numbers, source UniProt links, current full KEGG terms and network
degree, plus a narrative explaining the distinction between authored judgment
and database disagreement. README points to complete75_README.md.

Independent 06_check_poster_outputs.R run passed: exact lists75 and8added, all
source effects/FDR unchanged, observed>=2/3 halo eligibility, raw current KEGG
membership, node/edge/isolate reconciliation,49 original204edge-score regression,
nine nonempty figure exports and three geometry rows. Only existingRlocale and
package-build warnings. No commit/push requested in this turn.

AI assistance: Codex implemented retrieval/provenance work, a second agent
updated and inspected the plot, a third independently extended validation,
and a read-only agent traced original flags and primary function sources.
