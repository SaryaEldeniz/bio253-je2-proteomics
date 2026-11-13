# Step 5: ORA and GSEA methods, in plain language

This guide describes the declared Step 5 protocol. The generated results report
contains the actual run counts, findings, warnings and figures. Enrichment groups
protein changes into functions; it does not measure pathway activity directly.

## 1. Start from the existing October protein results

The input is the Step 4 total-proteome comparison table, not the original screenshot
and not the August export. Source measurements and statistical models are not
refitted here. The ten comparisons are:

- PASN responses within the ancestor, kinase mutant and phosphatase mutant.
- Each mutant versus the ancestor in TSB and in SN.
- Three genotype-by-medium interactions: stk versus ancestor, stp versus ancestor,
  and stp versus stk.

The 15-minute pulse comparisons are outside this declared primary analysis.
The supplied protein-level effects, FDRs, statistics and estimate-type labels stay
unchanged. Step 2 sample validation remains skipped at the user's request.

## 2. Understand what is being counted

An experimental protein entry has a JE2 `protein_key`. Its mapped UniProt accession
is the reference identifier used for the primary enrichment. A pathway or term is
a set of those identifiers. Some sets are metabolic pathways, while GO functions,
GO locations and UniProt keywords describe other biological properties.

The primary mapping is to the supplied NCTC8325 reference-proteome orthologs. It
does not imply that the reference sequence is identical to JE2. Phosphosite
coordinates remain attached to JE2 and are not transferred by this analysis.

## 3. Define a primary universe before seeing pathway results

The primary universe retains ordinary protein-coding entries with finite source
effect, FDR and ranking statistic, a quality-checked UniProt mapping, and one
eligible JE2 entry per UniProt target within the total-proteome DEA universe.

The Step 3 map contains 2,229 ordinary total-DEA protein entries. Of these, 2,107
have good reference mappings, representing 2,080 distinct target accessions. There
are 49 source rows involved in many-to-one target collisions. Excluding all those
collision rows leaves 2,058 singleton entries before any contrast-specific finite
value checks. The generated coverage table records the final universe actually
used in each run.

Good mapping uses the supplied policy: BLAST E-value <=1e-5 and reference-length
ratio between 0.9 and 1.1. Weak mappings and target collisions remain in the audit;
they are not assigned whichever statistic would yield the strongest enrichment.

The same primary namespace is used for ORA and GSEA. This makes their findings
easier to compare. It also avoids counting one reference protein several times.

## 4. Keep the annotation snapshot fixed

The memberships come from the supplied local JE2/STRING terms, read from the
Step 3 `terms_cds.tsv.gz` table. For the primary analysis, these memberships are
transferred to the eligible UniProt IDs and duplicate target-term pairs are
removed. This is not a fresh download of all UniProt or STRING annotations.

Converting a one-to-one identifier does not itself improve an enrichment result;
it improves interoperability and provides a controlled connection to reference
resources. Annotation source and version can matter more than the spelling of
the identifier. Input hashes and the code record which snapshot was used.

The supplied local resource describes STRING v12 and the custom organism
`STRG0A35SHG`. Its exact minor release and GO/KEGG annotation dates are unknown.
The teaching slides also show STRING 12.5 examples for a different custom organism
(`STRG0B09NBK`); that does not establish the version of this analysis snapshot.
Live API checks use the mapped reference identifiers and reference taxon 93061,
because custom organisms are unsupported by the public STRING API. Those checks
are separate from the fixed memberships used for enrichment.

The declared collections are:

| Collection | What it groups |
|---|---|
| KEGG | Pathway memberships. |
| GO Biological Process | Biological processes in which proteins participate. |
| GO Molecular Function | Molecular functions associated with proteins. |
| GO Cellular Component | Cellular structures, complexes or locations. |
| UniProt keywords | Broad curated/assigned annotation categories. |

The other supplied collections are not added after inspecting favourable results.
Some broad categories may have little experimental specificity; use the individual
contributing proteins and literature to interpret them.

## 5. Use measured/tested backgrounds, not the whole genome

The background is every eligible tested identifier for the specific policy and
contrast. It includes identifiers with no annotation in a particular collection.
For example, the primary universe is not restricted to only KEGG-annotated
proteins when testing a KEGG term. This is a declared shared-background choice.

The foreground, ranks and term memberships must use the same identifier system.
The code reports collection-specific annotation coverage so lack of annotation
is not mistaken for complete knowledge of nonmembership. Another defensible
question could use collection-annotated-only backgrounds, but that would be a
different analysis with a different denominator.

For each run, intersect every term with the actual universe first. Retain terms
with effective size 10 through 500, excluding a set that spans the entire
universe. A term's effective size is the number of eligible measured/tested
members, not its total database or genome-wide size.

## 6. ORA: ask whether selected changes concentrate in a term

Within each contrast, a protein enters the selected foreground when the source
FDR is <0.1 and the absolute source log2 effect is >1. These are the same Step 4
thresholds, applied after the declared universe restrictions. This is a combined
significance/effect-size selection, not a newly fitted test against a twofold
threshold.

Positive and negative source effects form separate foregrounds. The code does
not combine opposite directions into one discovery list. It tests every eligible
term, including zero-overlap terms, rather than first discarding terms with poor
overlap and thereby shrinking the multiple-testing correction.

For one term, the table contains:

- N: number of proteins in the eligible background.
- n: number of selected proteins in the directional foreground.
- K: number of background proteins belonging to the term.
- k: number of foreground proteins belonging to the term.
- Expected overlap: n*K/N.
- Fold enrichment: (k/n)/(K/N), when n is nonzero.
- Upper-tail hypergeometric p-value: probability of at least k hits under this
  sampling model, calculated as `phyper(k-1, K, N-K, n, lower.tail=FALSE)`.

This upper-tail test is equivalent to the corresponding one-sided Fisher exact
test. A zero-overlap term has p=1. An empty directional foreground also has p=1;
fold enrichment is unavailable rather than divided by zero. Keep the member IDs
and readable labels so every result can be traced back to proteins.

## 7. GSEA: ask whether an entire term tends toward one end of the ranking

GSEA uses all eligible finite signed DEA statistics, including proteins that fail
the Step 4 significance threshold. A statistic combines direction with the
estimated uncertainty; it is not simply a sorted list of significant fold changes.
No top-hit cutoff is applied to this ranked input.

The rank vector has one value per unique identifier. It is sorted in descending
statistic order. Equal values are ordered deterministically by identifier, without
adding artificial noise. The tied-score fraction is recorded.

The engine is the official `fgsea::fgseaMultilevel` adaptive multilevel Monte Carlo
implementation. The declared settings are a fixed recorded seed, weighted score
parameter 1, standard signed scoring, sample size 101, initial simple-permutation
count at least 1,000, effective set sizes 10-500, and a single-process serial
backend. The run's exact settings and package version are saved with its outputs.

The p-value lower estimation bound is 1e-10. Results at that bound are labelled;
they are not reported as exact zeros. `log2err` and engine warnings remain visible.
An unestimable term stays missing with a review status, rather than receiving a
fabricated p-value. A fixed seed makes the software result reproducible; it does
not remove biological uncertainty.

GSEA outputs include the enrichment score (ES), normalized enrichment score
(NES), estimated p-value, adjusted values, effective size and leading-edge
members. The leading edge is the subset driving the ranked enrichment pattern.
It is a useful list for measurement plots, functional reading and later network
annotation, not a list of experimentally validated direct targets.

## 8. Read positive and negative results in the original contrast direction

For `TSB_vs_SN`, a positive effect/statistic is higher in TSB and a negative one
is higher in SN. Therefore positive NES concentrates toward the higher-TSB side;
negative NES concentrates toward the higher-SN side.

For interactions, define R_g=log2(SN_g)-log2(TSB_g). The three source interaction
effects are R_anc-R_stk, R_anc-R_stp and R_stk-R_stp. A positive NES on the first
interaction therefore indicates concentration toward proteins with more positive
PASN responses in the ancestor relative to stk. This can represent larger
increases or smaller decreases. It does not mean higher absolute abundance in
the ancestor, and it does not by itself establish pathway activation.

The three interactions are related comparisons, not independent experiments.
Their effect identity does not permit subtraction of their enrichment p-values,
NES values or adjusted significance levels.

## 9. Separate discovery from declared sensitivity analyses

| Policy | Identifier space and eligibility |
|---|---|
| `primary` | Quality-checked, unique UniProt targets; all source-model labels retained. |
| `replicated` | Primary mapping plus >=2/3 observed values in every relevant group. |
| `observed_replicated` | Replicated policy plus source model label `observed`. |
| `local_JE2` | All ordinary eligible JE2 protein entries using their supplied local terms; no UniProt-quality or collision exclusion. |

Each sensitivity policy defines its own eligible background, directional
foregrounds and complete ranked vector. An interaction uses four relevant groups;
a simple comparison uses two. Restricting only the foreground while leaving an
unmatched background would answer a different question.

The source `lod_imputed` label is model-wide. It can arise from a group outside
the displayed contrast. The primary policy therefore retains those rows and
reports sensitivity instead of assuming every flagged contrast lacks measured
support. Conversely, an `observed` label does not guarantee dense replication.

Local-JE2 sensitivity changes both the identifier unit and mapping eligibility;
differences are not evidence that one naming system is intrinsically better.
The policy with the smallest p-value is not selected retrospectively as the main
result. A primary result that persists in sensible sensitivities is easier to
interpret; a sensitivity-only result is labelled accordingly.

## 10. Account for the number of questions tested

ORA and GSEA are adjusted separately. For each method and policy:

- `FDR_local` applies BH correction across the tested terms/collections for one
  contrast. ORA includes both directional foregrounds in that local family.
- `FDR_study` applies BH correction across all ten contrasts and tested
  collections in that policy. ORA again includes both directions.
- The original fgsea engine-adjusted value is retained separately, rather than
  silently overwritten by this project's broader correction.

The main discovery threshold is study-adjusted FDR <0.1 in the predeclared
primary policy. Local-adjusted and sensitivity results remain available for
transparency. These are enrichment-level adjusted p-values, distinct from the
protein-level source FDR used to select the ORA foreground.

Eligible failed GSEA tests remain missing but contribute to the declared family
size when adjusting the finite p-values. Testing more sensitivity policies does
not justify choosing whichever produces the most favourable answer.

## 11. Reduce presentation redundancy without hiding the full table

GO terms overlap heavily and can describe the same proteins at different levels.
Ten enriched terms need not mean ten independent biological mechanisms. Figures
show readable representative subsets with a documented selection rule, while the
full tables retain every tested term and its members. A pathway heatmap shows
patterns across contrasts; it does not replace an interaction test.

Use the strongest supported terms, contributing proteins, their sample values
and primary literature together when naming modules. Report results with no
study-level significance honestly rather than forcing a module name.

## 12. What this step does not prove

ORA is a selected-list overrepresentation test, and preranked GSEA is a
rank/gene-set enrichment test. GSEA does not perform sample-label permutation
here and does not add replicate samples to the experiment. Both methods inherit
the source DEA and annotation limitations, including correlated proteins and
overlapping sets.

Neither proves metabolic flux, direct protein binding, kinase-substrate
relationships or causal regulation. The documented reduction of Stk1 in the stp
strain complicates pure phosphatase-loss interpretations. PASN is mixed spent
medium, so these comparisons do not isolate secreted toxins from nutrient effects.

Phosphosite-level enrichment and motif analysis remain outside this protein-level
Step 5 protocol until localization and parent-protein abundance have been examined.

## Requested observed-only display, 2026-11-13

The [all-results index](../reports/step05/observed_only/README.md) is a separate
analysis with `policy="observed_only"`. It retains quality-mapped singleton CDS
entries with finite source statistics and source model label `observed` (1,858
proteins per comparison), without adding a replicate-count filter. Source
`lod_imputed` rows are excluded from both foregrounds and the background/ranking.
This is different from the earlier `observed_replicated` sensitivity.

ORA and GSEA are rerun for this universe, keeping all five original collections
in the study-level BH families across the ten comparisons. The displayed/exported
GO BP/MF/CC and KEGG findings use FDR_study<0.1 and the unchanged effective set
size filter of 10-500. No further Top-N limit is applied. GSEA panels sort by
descending signed NES, including negative scores at the bottom. ORA directions
remain separate; each panel sorts by descending fold enrichment. Per-plot TSVs
record the exact displayed order. Comparisons with no passing terms get explicit
empty plots; an empty result is not evidence of biological absence.

The two individual-protein sample plots retain the selected proteins but now
show only TSB and SN. Their values are measured normalized abundances, not
imputed values. Their original primary-based protein selection is retained and
is distinct from the observed-only enrichment analysis.

## Primary method references

- [Official fgsea tutorial](https://bioconductor.org/packages/release/bioc/vignettes/fgsea/inst/doc/fgsea-tutorial.html): preranked input, multilevel estimation, size limits, leading edges, enrichment plots and overrepresentation universes.
- [R hypergeometric documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Hypergeometric.html): upper-tail probability definition and parameters.
- [R multiple-testing documentation](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html): BH adjustment and declared family size.
