# Functional summary of the three total-proteome interactions

This figure connects the protein-level interaction volcanoes to their functional annotations. It combines all significant observed-only KEGG GSEA findings across the three interactions with a compact, explicitly representative GO Biological Process ORA panel. It does not refit the protein models, rerun enrichment, update annotation memberships or reinterpret nonsignificant results as zero effects.

## Figure and reproduction

- [Poster SVG](../../figures/poster_functional_summary/functional_summary_interactions.svg).
- [PNG preview](../../figures/poster_functional_summary/functional_summary_interactions.png).
- [PDF](../../figures/poster_functional_summary/functional_summary_interactions.pdf).
- [R exporter](../../scripts/05_poster_functional_summary.R).

Open `Sarya_George.Rproj` and run:

```r
source("steps/05_enrichment/scripts/05_poster_functional_summary.R",
       encoding = "UTF-8")
```

The exporter reads the existing observed-only Step 5 results. Full enrichment caches remain local under ignored `results/step05/observed_only/`. Original tables and figures remain available in the [complete observed-only results](../observed_only/README.md). The [Step 5 methods](../../docs/step05_methods.md) document eligibility, ranking, background, annotation snapshot and multiplicity choices.

## What the KEGG matrix includes

There are **13 significant pathway–interaction pairs across 12 distinct KEGG pathways**: 7 in stk versus ancestor, 3 in stp versus ancestor, and 3 in stp versus stk. Purine metabolism occurs in two comparisons. The display includes their full union, without a Top-N cutoff or a restriction to the four historically highlighted modules.

Columns follow the same genotype-comparison order as the three volcano panels. Rows are grouped by the first comparison in which they are significant, then ordered by descending signed NES within that block. This preserves one common row order across the three columns. A numeric/coloured cell denotes study-wide FDR <0.10; a dash denotes a tested pathway that did not pass that threshold. A dash is **not** a zero NES, proof of absence or proof that two genotypes respond identically.

| Interaction | KEGG pathway | NES | Study FDR | Effective set size | Leading-edge proteins | PASN response more positive toward |
| --- | --- | ---: | ---: | ---: | ---: | --- |
| stk versus ancestor | Purine metabolism | +1.960 | 0.00142 | 73 | 21 | ancestor |
| stk versus ancestor | Lysine biosynthesis | +1.885 | 0.0215 | 22 | 11 | ancestor |
| stk versus ancestor | Two-component system | +1.532 | 0.0388 | 122 | 42 | ancestor |
| stk versus ancestor | Biosynthesis of amino acids | +1.502 | 0.0667 | 96 | 38 | ancestor |
| stk versus ancestor | RNA degradation | −1.684 | 0.0542 | 20 | 12 | stk |
| stk versus ancestor | Arginine biosynthesis | −1.766 | 0.0414 | 16 | 8 | stk |
| stk versus ancestor | Oxidative phosphorylation | −1.891 | 0.0152 | 24 | 12 | stk |
| stp versus ancestor | Nitrogen metabolism | +1.865 | 0.0198 | 12 | 5 | ancestor |
| stp versus ancestor | Ribosome | −1.582 | 0.0579 | 56 | 25 | stp |
| stp versus ancestor | Fructose and mannose metabolism | −1.858 | 0.0364 | 17 | 7 | stp |
| stp versus stk | Aminoacyl-tRNA biosynthesis | +1.870 | 0.00540 | 27 | 12 | stk |
| stp versus stk | Glycerolipid metabolism | +1.833 | 0.0181 | 23 | 10 | stk |
| stp versus stk | Purine metabolism | −1.829 | 0.00614 | 73 | 23 | stp |

Values are rounded for reading; exact estimates are in the three source tables: [stk–ancestor](../observed_only/step05_KEGG__interaction_stk_vs_anc_at_TSB_vs_SN__GSEA__all_observed.tsv), [stp–ancestor](../observed_only/step05_KEGG__interaction_stp_vs_anc_at_TSB_vs_SN__GSEA__all_observed.tsv), and [stp–stk](../observed_only/step05_KEGG__interaction_stp_vs_stk_at_TSB_vs_SN__GSEA__all_observed.tsv).

Leading edges are the proteins contributing to the ranked enrichment pattern, not proven causal regulators or an additive decomposition of NES. The effective set size is the number of eligible members in the tested background, rather than all genome/database members.

## Read the interaction and NES directions correctly

Let `R_g = log2(SN_g) − log2(TSB_g)`. The native interaction for a panel labelled first genotype versus second genotype is:

```text
I = (first_TSB − first_SN) − (second_TSB − second_SN)
  = R_second − R_first
```

GSEA ranks the supplied signed protein DEA statistic in descending order. Positive NES means pathway members concentrate toward positive interaction statistics; negative NES means concentration toward negative statistics. The statistic incorporates the source estimate's uncertainty; the ranked input is not simply a list of fold changes or significant proteins.

- In **stk versus ancestor**, positive NES favours a more positive SN-minus-TSB response in the ancestor; negative NES favours stk.
- In **stp versus ancestor**, positive NES favours a more positive response in the ancestor; negative NES favours stp.
- In **stp versus stk**, positive NES favours a more positive response in stk; negative NES favours stp.

“More positive” can mean a larger increase or a smaller decrease. A pathway's sign does not establish that it increases within either genotype, which genotype has the higher absolute abundance, or whether its metabolic flux/activity rises. NES is a normalized enrichment score, not a fold change or an estimate in log2 abundance units. The three interaction effects are algebraically related comparisons, not three independent experiments; their NES or FDR values cannot be subtracted to reconstruct one another.

## GO-BP ORA was not empty

There is no statistical reason to describe GO-BP ORA as having no results across this project. For the requested observed-only interactions it found **8 significant term–direction results for stk versus ancestor, and 0 for each of the other two comparisons**. The eight results consist of three distinct overlap groups:

| Representative term | Sign | Observed overlap | Foreground | Effective set size | Fold enrichment | Study FDR | Terms sharing exactly these hits |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Branched-chain amino acid biosynthetic process | positive | 5 | 51 | 10 | 18.216 | 0.000993 | 2 |
| De novo IMP biosynthetic process | positive | 4 | 51 | 11 | 13.248 | 0.0189 | 5 |
| Cytolysis | negative | 3 | 27 | 11 | 18.768 | 0.0463 | 1 |

The representative for each group is the term with the smallest study FDR among terms with exactly the same directional overlap IDs. This reduces presentation redundancy without retesting or discarding the full results. The [full eight-term source table](../observed_only/step05_GO_BP__interaction_stk_vs_anc_at_TSB_vs_SN__ORA__all_observed.tsv) remains available. The [stp–ancestor](../observed_only/step05_GO_BP__interaction_stp_vs_anc_at_TSB_vs_SN__ORA__all_observed.tsv) and [stp–stk](../observed_only/step05_GO_BP__interaction_stp_vs_stk_at_TSB_vs_SN__ORA__all_observed.tsv) tables contain headers and no passing rows.

The plotted representatives are then sorted by descending fold enrichment: cytolysis, branched-chain amino-acid biosynthesis, and de novo IMP biosynthesis. Bar lengths use a shared zero baseline and scale; endpoint labels are fold enrichment, not protein fold changes. Positive/negative describes the direction of the selected protein interactions; all bars show overrepresentation.

- The two branched-chain amino-acid terms share **ilvA, leuB, leuA, B7H15_11390 and ilvD**.
- The five IMP/purine-monophosphate terms share **purM, purF, purS and purC**.
- Cytolysis contains **hlgB, B7H15_12080 and hly** in the supplied mapping.

The IMP result supports the purine GSEA theme using the subset of proteins that individually pass the volcano selection rule. Branched-chain amino-acid biosynthesis adds detail to the broader amino-acid response. Cytolysis adds a distinct annotated theme that is not contained in the significant KEGG list; it should not be dismissed as wholly redundant. It does not show that bacterial lysis, haemolysis or toxin secretion was measured or increased in this experiment.

Leaving the entire GO-BP ORA plot off the poster can be an **editorial space/redundancy choice**, with these results retained in the supplementary material. It would be inaccurate to justify that choice by saying GO-BP ORA was negative everywhere. The compact three-representative panel lets the poster retain its relevant support and the cytolysis result without displaying eight overlapping GO rows. No passing terms in the other comparisons means no detected overrepresentation at the stated threshold and background; it does not establish biological absence.

## Why ORA and GSEA need not find the same functions

ORA starts from selected positive/negative protein lists. GSEA uses all eligible ranked proteins, including individually nonsignificant estimates. A coherent rank shift across many modest effects can therefore produce GSEA enrichment without an ORA result. For the historically highlighted significant oxidative-phosphorylation and aminoacyl-tRNA results, each leading edge contains 12 proteins but **none individually passes both protein-selection thresholds**. That is a concrete reason to show GSEA alongside volcanoes rather than treating enrichment as merely a re-labelling of their coloured points.

Conversely, a small selected group can yield ORA enrichment without establishing a consistent shift throughout the entire ranked gene set. Neither method is a validation of the other: they ask related but different questions on the same underlying experiment, with overlapping annotations.

## Population, methods and why the counts differ from the volcanoes

The total-proteome volcanoes show 2,229 eligible JE2 protein entries per interaction, including observed and source LOD-imputed models. This functional summary uses the user's later **observed-only** enrichment policy: **1,858 uniquely quality-mapped UniProt entries per interaction**, excluding source LOD-imputed models from both foregrounds and backgrounds/ranks. Reference-target collisions, weak/absent mappings and finite-statistic eligibility also affect entry retention. Renaming JE2 IDs to UniProt alone does not create biological evidence or guarantee identical JE2/reference proteins.

This is not the stricter `observed_replicated` policy. No extra ≥2/3 measurement filter was added. The observed-only universes retain 102, 71 and 101 proteins, respectively, failing ≥2/3 measured support in at least one relevant group. The global observed flag does not guarantee three measured replicates in every group. [Universe and selection audit](../observed_only/universe_audit.tsv).

The ORA directional foregrounds are **51 positive/27 negative**, **29/32**, and **25/60**, respectively, using source protein FDR <0.10 and absolute log2 effect >1. The universe contains all eligible tested identifiers, including ones lacking an annotation in a particular collection; it is not the whole genome or only the KEGG-annotated proteins.

Memberships come from the frozen course/STRING terms transferred through the supplied mapping, rather than fresh live STRING network memberships. Sets are intersected with the eligible universe and retained at effective size 10–500. For every comparison, 71 KEGG and 419 GO-BP sets are eligible. Positive and negative ORA directions are tested separately, including zero-overlap terms. The hypergeometric test is an upper-tail overrepresentation test.

Preranked GSEA used `fgsea::fgseaMultilevel`, signed standard scoring, weighted parameter 1, deterministic descending score/identifier tie ordering, recorded seeds, sample size 101, initial simple-permutation count 1,000, and estimation floor 1e−10. No failed observed-only GSEA estimates occur in the saved run.

Significance here is **study-wide BH FDR <0.10**, adjusted separately for each method across all ten declared contrasts and all five original collections (including UniProt keywords, even though this figure shows only GO-BP/KEGG). The observed-only families have **7,860 GSEA tests** and **15,720 directional ORA tests**. These enrichment FDRs differ from the original protein FDR used for ORA hit selection. [Settings](../observed_only/settings.txt), [test counts](../observed_only/counts.tsv) and [run log](../observed_only/gsea_run_log.tsv).

Earlier primary-policy figures retained source LOD labels and used a different 2,058-protein background. For example, their KEGG GSEA counts for these interactions were 9/3/3 and GO-BP ORA counts were 6/0/0. The observed-only counts are 7/3/3 and 8/0/0. These are distinct declared universes, not interchangeable alternative numbers for the same plot. The older four-module overview was also a selected subset: five significant module–interaction pairs across four pathways, with seven nonsignificant pairs supplied separately. The present KEGG panel restores all 13 significant pairs across 12 pathways.

## Supported conclusions and practical limits

The interaction volcanoes identify individual protein changes whose TSB-to-PASN response differs between genotypes. The functional summary shows that those response differences also have organized pathway/process structure. Relative patterns involve nucleotide and amino-acid biosynthesis, energy-related protein groups, and translation-associated functions, with distinct signatures in the kinase and phosphatase mutant comparisons.

The stk comparison has opposing purine/amino-acid and oxidative-phosphorylation rank patterns. The stp comparison identifies nitrogen metabolism and opposing ribosomal/carbohydrate-related patterns. Direct mutant comparison distinguishes aminoacyl-tRNA/glycerolipid and purine patterns. The positive and negative purine results in two columns are consistent with the signed source contrasts; they do not make the two columns independent replications.

These are **associations of protein abundance responses with genotype and environment**. They do not prove direct Stk1/Stp1 substrates, protein binding, pathway activation, metabolic flux or a purely phosphatase-specific effect. The course notes document reduced Stk1 in the stp strain, and PASN is mixed spent medium with potential nutrient as well as secreted-product effects. Reference mapping and overlapping pathway memberships remain limitations; historically flagged unexpected memberships are retained to reproduce the reported NES rather than silently corrected after seeing the result. Live STRING assignment checks belong to the separate Step 6 network analysis and do not retrospectively redefine this enrichment snapshot.

Suggested poster conclusion:

> Stk1/Stp1 disruption is associated with distinct functional patterns of protein responses to PASN, involving nucleotide and amino-acid metabolism, energy-related proteins and translation-associated functions. GO-BP ORA supports selected biosynthetic changes in the stk comparison, while ranked GSEA additionally identifies coordinated shifts among proteins that do not individually pass the volcano thresholds.

AI assistance: Codex audited source tables, eligibility and score direction and prepared this interpretation. All numerical claims above are traceable to the saved Step 5 outputs and methods; no new biological measurement or enrichment analysis was performed for the summary figure.

## Figure caption, data and rendering checks

Suggested caption: **Functional patterns of genotype-dependent protein responses to PASN.** (A) KEGG preranked GSEA across the three genotype-by-medium interactions. All 12 pathways significant in at least one interaction are included; the 13 coloured numeric cells pass study-wide BH FDR <0.10. Values are normalized enrichment scores, not fold changes. Dashes mark tested nonsignificant cells. (B) Representative GO-BP ORA terms for the stk–ancestor interaction, selected by the smallest study FDR per exact hit-overlap group and displayed by descending fold enrichment. Bars use a common scale; annotations give study FDR and overlapping selected-protein counts. Eight significant GO terms share these three hit groups; the other two comparisons have no significant GO-BP ORA terms. Both methods use 1,858 quality-mapped observed-model proteins and the original course annotation snapshot. Enrichment describes relative abundance-response patterns, not pathway activity or flux.

- [All 36 KEGG matrix estimates](KEGG_matrix_all_36_tests.tsv).
- [All 13 significant KEGG pairs](KEGG_all_13_significant_pairs.tsv).
- [All eight GO-BP ORA terms](GO_BP_ORA_all_8_terms.tsv).
- [The three GO-BP representatives](GO_BP_ORA_3_representatives.tsv).
- [Counts by interaction](enrichment_counts.tsv).
- [Rendering checks](render_checks.txt).

At its intended 529.333×250-mm placed size the figure uses Arial 24-pt text, 26-pt panel headings and a 32-pt title. Blue/orange sign colours avoid the previously problematic green/light-blue combinations. The SVG consists of vector paths and shapes with explicit physical dimensions. Independent QA matched all 36 matrix estimates/FDRs to the full source cache, confirmed all 13 passing and 23 nonsignificant cells, independently checked the ORA overlaps/ratios/signs/representative rule, parsed the SVG with no raster images, and checked the PNG for clipping and overlap. No tests were added to mirror rendering implementation; source-consistency assertions run in the exporter.
