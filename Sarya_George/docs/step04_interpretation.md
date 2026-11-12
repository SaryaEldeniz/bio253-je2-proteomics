# Step 4: understanding the comparison tables

This step selects the relevant results from the lecturer's existing differential
analysis and attaches the Step 3 protein annotations. It does not refit the models.
Step 2 remains skipped at the user's request. Individual measurements are shown
for selected candidates to help interpret these particular comparisons.

## What the ten comparisons ask

The first group describes the PASN response within each genotype:

| Contrast | Question | Meaning of a positive source `diff` |
|---|---|---|
| `TSB_vs_SN_at_anc` | What changes after long PASN exposure in the ancestor? | Higher in ancestor TSB, therefore lower in ancestor SN. |
| `TSB_vs_SN_at_stk` | What changes after long PASN exposure in the kinase mutant? | Higher in kinase-mutant TSB, therefore lower in kinase-mutant SN. |
| `TSB_vs_SN_at_stp` | What changes after long PASN exposure in the phosphatase mutant? | Higher in phosphatase-mutant TSB, therefore lower in phosphatase-mutant SN. |

For all three, a negative source `diff` means higher in SN. The source direction
is retained rather than silently reversed. TSB is the control and SN is the long
exposure to cell-free PASN.

The second group compares a mutant with the ancestor in one medium:

| Contrast | Question | Meaning of a positive source `diff` |
|---|---|---|
| `stk_vs_anc_at_TSB` | How does the kinase mutant differ in control medium? | Higher in the kinase mutant in TSB. |
| `stp_vs_anc_at_TSB` | How does the phosphatase mutant differ in control medium? | Higher in the phosphatase mutant in TSB. |
| `stk_vs_anc_at_SN` | How does the kinase mutant differ under long PASN exposure? | Higher in the kinase mutant in SN. |
| `stp_vs_anc_at_SN` | How does the phosphatase mutant differ under long PASN exposure? | Higher in the phosphatase mutant in SN. |

The final group asks whether the PASN response differs between genotypes. These
are statistical interactions, distinct from the protein-network edges in STRING:

| Contrast | Question |
|---|---|
| `interaction_stk_vs_anc_at_TSB_vs_SN` | Does the kinase mutant's PASN response differ from the ancestor's? |
| `interaction_stp_vs_anc_at_TSB_vs_SN` | Does the phosphatase mutant's PASN response differ from the ancestor's? |
| `interaction_stp_vs_stk_at_TSB_vs_SN` | Do the two mutants differ in their PASN responses? |

The first two interactions compare each mutant with the ancestor. The last is a
follow-up comparison between the mutants. The 15-minute SNpulse comparisons remain
outside this main Step 4 set.

## Read an interaction as a difference between responses

Define a protein's PASN response in genotype g as:

`R_g = log2(SN_g) - log2(TSB_g)`.

Then the three interaction effects are:

- I1, kinase versus ancestor: `R_anc - R_stk`.
- I2, phosphatase versus ancestor: `R_anc - R_stp`.
- I3, phosphatase versus kinase: `R_stk - R_stp`.

For I1, a positive value means the ancestor response is more positive than the
kinase-mutant response. A negative value means the kinase-mutant response is more
positive. The same principle applies to I2 and I3 using their listed references.
"More positive" could mean a larger increase or a smaller decrease; inspect the
within-genotype responses before describing the biological pattern.

For example, suppose an illustrative protein rises fourfold in the ancestor but
twofold in the kinase mutant after PASN exposure. Their log2 responses are 2 and
1. I1 is therefore `2 - 1 = 1`. The ancestor's response is more positive. This
example alone says nothing about statistical significance, and an effect exactly
equal to 1 does not pass a strict `abs(diff) > 1` selection rule.

The three effect estimates satisfy `I3 = I2 - I1`. They are therefore not three
independent pieces of evidence. This identity concerns the effect estimates;
standard errors, test statistics, p-values and FDR values must not be subtracted.
An interaction magnitude above 1 represents a response-ratio difference above
twofold, rather than a simple twofold protein abundance difference.

## Statistical selection and supporting measurements

The course's descriptive hit selection is `FDR < 0.1` and `abs(diff) > 1`.
Both conditions must be satisfied. This combines statistical support with an
effect-size filter; it is not automatically a formal test against the twofold
threshold. A row that fails selection is not proof that its effect is zero.

The supplied FDR values are preserved unchanged. Annotation filtering does not
trigger a new p-value adjustment in Step 4. Do not claim the source FDR controls
errors jointly across all ten comparisons unless the source analysis establishes
that family of tests. Likewise, a significant response in one genotype and a
nonsignificant response in another do not by themselves establish a difference
between genotypes; use the interaction comparison.

The source `estimate_type` is preserved as well. Its labels require care in this
dataset: all selected total-proteome contrasts have the same source-label counts,
and some `lod_imputed` rows have three measured values in both groups of a displayed
comparison. The label can reflect missingness elsewhere in the model, rather than
only in the groups displayed in that particular comparison.

Therefore read the source label alongside the contrast-specific observation
counts. Those counts answer how many actual measurements support each group used
by that contrast. A simple comparison uses two groups; a genotype-by-medium
interaction uses four. A source label of `observed` does not guarantee three
measurements per relevant group, and a label of `lod_imputed` does not identify
which relevant group, if any, is missing all measurements.

The conservative sensitivity view requires both a source label of `observed` and
at least two measured replicates in every relevant group. This is a transparent
restriction of the existing results, not a new fitted model or a new FDR
calculation. Preserve and inspect results outside that subset rather than treating
them as biologically absent. Never replace missing log2 measurements by zero.

## Total-proteome results and preliminary phosphosite results

The primary Step 4 ordinary-protein analysis uses metadata `entry_class == "cds"`.
The full annotated audit retains the other source entries. Some contaminants were
not flagged in the DEA's original contaminant column; Step 3's richer annotation
classification supplies the ordinary-protein filter. UniProt mapping quality and
many-to-one mappings remain visible. A weak or absent reference mapping does not
make a JE2 measurement disappear from the ordinary-protein table.

The total-proteome table describes protein abundance changes. It supplies the
main Step 4 plots and the later functional-enrichment inputs. The preliminary
phosphosite table attaches source site effects to their parent-protein annotation.
Its rows remain sites, not unique proteins. Site-level rows must not be counted as
independent proteins in a protein-enrichment background.

Preliminary phosphosite selection is not a final list of regulated phosphorylation
sites. Localization confidence, site detection, ambiguous assignments and parent
protein abundance require the planned Step 7 evaluation. Keep exact JE2 sequence
coordinates even when functional annotation uses a mapped UniProt ortholog. A
phosphosite signal change can partly reflect a change in the amount of its parent
protein; subtraction or a visual comparison alone does not supply an adjusted
significance test.

## What can be concluded at this stage?

Step 4 establishes which existing comparison rows meet the declared selection,
their directions, annotation and supporting measurement patterns. It prepares
auditable inputs for pathway analysis. It does not yet establish enriched modules,
direct kinase substrates, phosphorylation motifs or causal regulation.

Keep the documented reduced Stk1 abundance in the stp strain in mind when
interpreting mutant comparisons. A stp effect is not automatically a pure
phosphatase-loss effect. PASN is mixed spent medium; these comparisons alone do
not separate secreted-molecule effects from nutrient depletion.
