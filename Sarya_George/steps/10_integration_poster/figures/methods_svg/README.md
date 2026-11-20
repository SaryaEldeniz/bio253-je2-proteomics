# Poster methods: editable vector illustrations

Original scientific illustrations in a blue, teal and mint palette, with navy
labels and outlines. All text and drawings remain SVG objects; there are no
embedded raster images or external assets. The individual panels have a
transparent background; the combined layout has a removable white background.
These are original drawings, not BioRender exports or licensed BioRender assets.

## Files

- `00_complete_methods.svg`: all three panels in a single vertical layout.
- `01_experimental_design.svg`: PAO1 spent medium, sterile filtration, the
  JE2 ancestor/kinase/phosphatase mutants, and the nine-group sample matrix.
- `02_proteomics_layers.svg`: protein processing, total-proteome and Fe-IMAC
  branches, LC-MS/MS and the two measurement outputs.
- `03_analysis_workflow.svg`: identifier mapping, protein interactions,
  enrichment, contrast-based STRING selection, phosphosite comparisons,
  Venn sets, JE2 sequence logos and parent-protein/site integration.
- Matching PNGs are previews. Use the SVGs in the poster for scalable artwork.
- `JE2_methods_SVGs.zip` bundles the four SVGs and this guide.

Insert SVGs as pictures in PowerPoint, Inkscape or Illustrator. In applications
supporting SVG-to-shape conversion, ungroup to edit individual elements; the
SVG source also retains live text. Preserve the aspect ratio when resizing.
Use the separate panels when the combined vertical version does not suit the
poster layout. Labels and symbols supplement colour, so colours alone do not
identify genotype, condition or measurement layer.

## Source-verified design

The primary design source is the original course `Proteomics-Data/README.md`,
Biology in10lines and Design sections, cross-checked with the configured
`data/derived/samples.tsv`. Local project records are
`../../reports/step01_design.md` and `../../notebook/1_ExperimentalDesign.md`.

- Three genotypes: JE2 ancestor (`anc`), `stk1::Tn` (`stk`), `stp1::Tn` (`stp`).
- Three conditions: TSB control; TSB growth then15minPASN (`SNpulse`); approximately
  15h in70%TSB+30%PASN (`SN`). The pulse-specific final mixture is not independently
  detailed, so no extra pulse concentration is invented.
- Each genotype/condition has three biological replicates:27samples total.
  Each of those samples has total-protein and phospho-enriched measurement
  layers; these are not54biological replicates.
- PASN is cell-free, sterile-filtered PAO1 spent medium. PAO1 is illustrated
  upstream of filtration and is not mixed alive with JE2.
- Current project methods are based on October analysis exports. The focal
  protein interactions compare TSB with long SN; the pulse remains visible in
  the full experimental design.
- The stp strain also has reduced Stk1, as documented in course guidance.
  The gene diagram marks the transposon insertion, without implying an isolated
  phosphatase-only physiological effect.

Generic protein-processing and LC-MS/MS icons illustrate the verified assay
class. No unverified extraction chemistry, digestion enzyme, instrument model
or preparation duration is specified. The older Day2 ancestor/evolved lecture
and the Schwyter publication describe a different experiment and were not used
to import their sample counts, exposure times or instrument details.

## Analysis distinctions

The top analysis branch concerns total-protein genotype-by-medium interactions.
The lower branch concerns phosphosite mutant-versus-ancestor comparisons within
TSB and SN, with separate phosphosite interaction tests. Venn sets are therefore
not depicted as derived from interaction-selected hits.

GO/KEGG ORA and ranked GSEA use the observed-only analysis. Network candidates
come from protein contrasts and measurement support, independently of NES or
enrichment significance; the direct arrow bypasses enrichment. STRING12.5 is
functional at confidence>=0.4. Its edges do not establish direct phosphorylation.
The integration endpoint is comparison and candidate prioritisation, not already
completed occupancy modelling or functional validation.

Detailed criteria remain in `../../docs/step05_methods.md`,
`../../reports/step06/poster_networks/complete75_README.md`,
`../../reports/step09_exclusive_phosphosite_proteins.md` and
`../../reports/step07_phospho_parent_comparison.md`. Enrichment and network
selection do not share one identical filter. STYH denotes a search setting,
not proof that every displayed protein has all four residue modifications.

Icons resembling charts, networks and sequence logos are schematic method
symbols; their illustrative colours and shapes are not experimental results.

## Reproduce

From `Sarya_George.Rproj`:

```r
source("steps/10_integration_poster/scripts/methods_illustrations/01_render_design_analysis_svgs.R")
source("steps/10_integration_poster/scripts/methods_illustrations/02_render_omics_methods_svg.R")
source("steps/10_integration_poster/scripts/methods_illustrations/03_assemble_methods_svgs.R")
```

These R scripts author the SVGs using native vector elements and `xml2`.
Optional PNG previews use `steps/10_integration_poster/scripts/methods_illustrations/render_methods_previews.cjs`
with Node.js and sharp. The vector files themselves need neither Node nor
a browser to remain editable. Structural checks are saved in `validation.txt`.
