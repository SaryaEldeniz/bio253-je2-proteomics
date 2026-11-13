# Step 5: which modules are worth following up?

The clearest main story is a coordinated decrease in purine-biosynthesis protein
abundance during long PASN exposure, with a stronger decrease in the kinase
mutant. A useful second story is a different nitrate/nitrogen response in the
phosphatase mutant. These are abundance-response modules to investigate in the
later phosphosite analysis; they are not established direct Stk1/Stp1 substrates.

All enrichment FDRs below are the study-wide BH values for the stated method and
policy, across ten contrasts and five annotation collections. The main threshold
is FDR <0.1. Sensitivity ranges include all four predeclared policies: `primary`,
`replicated`, `observed_replicated` and `local_JE2`. Each policy has its own tested
background and correction family. Agreement between policies is a robustness
check, not an independent replication experiment.

The complete evidence is in the [primary GSEA table](gsea_primary_significant.tsv),
[primary ORA table](ora_primary_significant.tsv),
[universe/foreground audit](universe_and_hit_counts.tsv).
Full all-tested tables are generated locally as `results/step05/gsea_all.tsv.gz`
and `results/step05/ora_all.tsv.gz`; those detailed caches are excluded from Git.
The [methods guide](../../docs/step05_methods.md) gives the exact selection rules.

## Read the interaction signs first

Define `R_g = log2(SN_g) - log2(TSB_g)`, the PASN response in genotype `g`.

| Native interaction | Source effect | Positive end of its ranking |
|---|---|---|
| `interaction_stk_vs_anc_at_TSB_vs_SN` | `R_anc - R_stk` | More positive PASN response in the ancestor. |
| `interaction_stp_vs_anc_at_TSB_vs_SN` | `R_anc - R_stp` | More positive PASN response in the ancestor. |
| `interaction_stp_vs_stk_at_TSB_vs_SN` | `R_stk - R_stp` | More positive PASN response in the kinase mutant. |

"More positive" can mean a larger increase or a smaller decrease. A positive NES
on the first interaction does not mean that pathway proteins increase in SN.
Use the within-genotype changes and measured sample values to distinguish those
possibilities. The native formulas and sign labels are also preserved in the
[comparison plan](comparison_plan.tsv).

## 1. Main module: purine biosynthesis, especially de novo IMP synthesis

The broad [KEGG purine-metabolism map](https://rest.kegg.jp/get/map00230) is the
strongest coherent kinase-dependent response among the inspected KEGG modules.
The narrower GO term "de novo IMP biosynthetic process" points to the same
biosynthetic machinery, rather than treating every member of the broad purine
map as equivalent.

| Evidence | Primary result | Result across all four policies |
|---|---|---|
| Purine metabolism, kinase-versus-ancestor interaction, GSEA | NES +1.932; FDR 0.000696 | Positive NES throughout; FDR 0.000286-0.000854. |
| Purine metabolism, phosphatase-versus-kinase interaction, GSEA | NES -1.786; FDR 0.00512 | Negative NES throughout; FDR 0.00336-0.00783. |
| Purine metabolism, kinase-versus-ancestor positive ORA | 8 hits among 57 selected proteins; 73 tested pathway members; 3.96-fold enrichment; FDR 0.0723 | FDR 0.0119-0.0723. |
| De novo IMP biosynthetic process, kinase-versus-ancestor positive ORA | 4/57 selected proteins, 11 tested term members; 13.13-fold enrichment; FDR 0.0211 | This row is supporting primary evidence; the broad purine sensitivity results above are the declared robustness comparison. |

The narrow ORA hits are `purM`, `purF`, `purS` and `purC`. The GSEA leading edges
also include other purine enzymes; GSEA is allowed to use their ranked evidence
even when their individual effects fail the ORA cutoff.

The direction is a **stronger fall in purine-enzyme abundance in long SN in the
kinase mutant**. Purine metabolism has positive environmental NES in all three
`TSB_vs_SN` comparisons, meaning higher abundance toward the TSB end of those
rankings: ancestor +1.917 (FDR 0.000348), kinase mutant +2.052 (FDR 0.0000608),
and phosphatase mutant +1.863 (FDR 0.00112). The positive kinase-versus-ancestor
interaction then means that the SN-minus-TSB response is more negative in the
kinase mutant. The negative phosphatase-versus-kinase interaction is consistent
with a larger fall in the kinase mutant than in the phosphatase mutant.

The [measured sample-value plot](../../figures/step05/step05_candidate_sample_values.png)
shows six purine enzymes with 3/3 measured values in each displayed TSB/SN group.
Their abundances fall visibly in long SN. The 15-minute pulse is omitted from
the revised sample plots and is outside this ten-contrast enrichment analysis.
These six identifiers were also checked against current reference entries:
[PurD Q2FZI5](https://rest.uniprot.org/uniprotkb/Q2FZI5.json),
[PurH Q2FZI6](https://rest.uniprot.org/uniprotkb/Q2FZI6.json),
[PurM Q2FZI8](https://rest.uniprot.org/uniprotkb/Q2FZI8.json),
[PurF Q2FZI9](https://rest.uniprot.org/uniprotkb/Q2FZI9.json),
[PurL Q2FZJ0](https://rest.uniprot.org/uniprotkb/Q2FZJ0.json), and
[PurQ Q2FZJ1](https://rest.uniprot.org/uniprotkb/Q2FZJ1.json).
These are mapped NCTC8325 reference identifiers, not a claim of identical JE2
phosphosite coordinates.

One concrete example is PurF (`B7H15_05545`, Q2FZI9). Its source
kinase-versus-ancestor interaction effect is +1.314, with protein-level FDR
0.0000197 and 3/3 measurements in every relevant group. Other plotted members
have smaller effects, which helps explain why the ORA foreground and GSEA
leading edge are different lists. The
[candidate measurement-support table](selected_candidate_measurement_support.tsv)
records each source effect, FDR and detection pattern.

**The phosphatase-versus-ancestor purine interaction is not significant:** primary
NES +0.996, FDR 0.725; all-policy FDR 0.636-0.725. Thus the evidence supports a
kinase-associated difference and a difference between the mutants, without
establishing a detectable whole-purine-module difference between stp and ancestor.
An individual purine protein can still differ in that contrast.

A defensible presentation sentence is: "Long PASN exposure is associated with
lower abundance of purine-biosynthesis proteins, and this decrease is stronger
in the kinase mutant than in the ancestor and phosphatase mutant." Avoid calling
this activation of purine metabolism or a measured change in purine flux.

## 2. Second module: nitrate/nitrogen response in the phosphatase mutant

The [KEGG nitrogen-metabolism map](https://rest.kegg.jp/get/map00910) is enriched
at the positive end of the phosphatase-versus-ancestor interaction: primary
NES +1.838, FDR 0.0257, with 12 tested pathway members. It remains positive and
significant in all four policies, with FDR 0.0115-0.0257. This is GSEA-supported
module evidence; it is not an additional significant primary KEGG ORA result.

Its leading edge contains respiratory nitrate-reductase subunits
`B7H15_13270`, `B7H15_13280` and `B7H15_13285`, the annotated nitrite-reductase
small subunit `B7H15_13295`, and glutamate-synthase large subunit `B7H15_02470`.
All five have the source-model label `observed` and 3/3 measured values in the
ancestor TSB/SN and phosphatase-mutant TSB/SN groups.

The nitrate-reductase component gives a concrete interpretation of the positive
interaction: **its SN induction is larger in the ancestor than in the stp
mutant**. For the alpha subunit (`B7H15_13285`), source `TSB-SN` effects are
-2.813 in the ancestor and -1.908 in stp. Therefore `R_anc-R_stp` is +0.905,
with protein-level FDR 0.000476. The gamma and beta subunits show the same
direction, with interaction effects +1.495 and +0.583, respectively. The
nitrite-reductase and glutamate-synthase proteins broaden the module beyond a
single nitrate-reductase complex; do not assume identical absolute changes
for every pathway member.

The useful hypothesis is that the phosphatase-mutant background alters the
nitrate/nitrogen response to long PASN exposure. Protein abundance does not
establish nitrate consumption, respiratory flux, oxygen availability, or a direct
phosphorylation mechanism. Several contributing subunits may be coregulated;
they are not independent biological replicates.

## 3. Additional robust lead: oxidative phosphorylation in the kinase mutant

For the kinase-versus-ancestor interaction, oxidative phosphorylation has
primary NES -1.800 and FDR 0.0156. Its direction remains negative in every
policy, with FDR 0.0154-0.0182. The leading edge includes ATP-synthase subunits
and respiratory-chain-associated proteins, including `atpD`, `atpH`, `atpG`,
`atpA`, `qoxB` and `qoxA`.

The negative interaction means a more positive SN-minus-TSB abundance response
in stk than in the ancestor. Separately, the `stk_vs_anc_at_SN` GSEA result is
positive (NES +1.985, FDR 0.00106), placing these proteins toward higher abundance
in stk during SN exposure. These are related abundance observations. They do
not establish greater ATP production or respiratory activity. This is a useful
third module if there is time for sample-value and phosphosite follow-up.

## 4. Strong ORA candidate with an important limit: branched-chain amino acids

"Valine, leucine and isoleucine biosynthesis" is strongly overrepresented in the
positive kinase-versus-ancestor ORA foreground: 5 of 57 selected proteins belong
to the 12-member tested set, giving 15.04-fold enrichment and FDR 0.00248.
The five hits are `ilvA`, `leuB`, `leuA`, `B7H15_11390` and `ilvD`.
The full local-JE2 sensitivity also supports this ORA result (FDR 0.00219).

However, the corresponding whole-set primary GSEA is not significant at the
study threshold (NES +1.552, FDR 0.168). The replicated and
observed-replicated policies leave only seven eligible set members, below the
predeclared minimum size of ten, so that term is **not tested** in those policies.
Its missing sensitivity rows are not negative results. Present this as a
focused ORA-supported candidate, rather than claiming agreement between both
methods or all sensitivity analyses.

## Place Stk1/Stp1 in the story without overclaiming

The current reference identities of
[Stk1/PrkC Q2FZ64](https://rest.uniprot.org/uniprotkb/Q2FZ64.json) and
[Stp1 Q2FZ65](https://rest.uniprot.org/uniprotkb/Q2FZ65.json) were verified in the
[UniProt checks](live_uniprot_checks.tsv). The
[resource provenance note](resource_checks.md) distinguishes the frozen course
enrichment memberships from the small live STRING context network. Functional
network edges do not demonstrate direct binding or kinase-substrate relationships.

The sample-value plot also makes the existing mutant caveat visible. Stp1 values
in the stp background are highly variable, with a high sample in the SN
group; its TSB group has only 1/3 measured values. Stk1 has 0/3 measured values
in stp SN. A missing detection is not a measured zero. Consequently, the stp
phenotype cannot automatically be interpreted as an isolated phosphatase effect
with normal kinase abundance. This step preserves the supplied DEA rather than
revalidating the mutants or refitting its model.

For the next phase, prioritize localized phosphosites on the supported module
members, examine parent-protein abundance, and then inspect motif context.
The enrichment and leading-edge lists nominate candidates for those checks;
they do not by themselves establish where phosphorylation explains the response.
