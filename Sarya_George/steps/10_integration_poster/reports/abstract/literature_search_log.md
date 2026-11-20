# Literature research for the abstract

The student explicitly requires primary references only. All three core and
two supporting scientific articles report original experiments; none is a
review. Original experimental designs/results and MEDLINE publication types
were checked. The lecture PDFs provide writing/search guidance only and are
excluded from the scientific bibliography. `abstract_core_references.ris`
contains the three core papers; `references.ris` contains all five primary
research papers. Article links point to the original publications via DOI.

Search date:19November2026. Applied the supplied bibliography lecture's
database/author-search approach,reference-management advice and emphasis on
scientific merit rather than journal impact factor. The lecture lists options;
it does not require searching every named platform. Used PubMed and primary
article records,with existing course papers as starting points. No impact-factor
ranking or unverified AI-generated citation was used.

## Bounded discovery queries

- PubMed-domain DOI searches: `10.1074/mcp.RA120.002232`,
  `10.1128/msystems.00262-25`, `10.1128/mbio.03153-22`,
  `10.1128/jb.00101-25`.
- Author/title searches: `"Serine-threonine phosphoregulation by PknB and Stp"`,
  `Huemer PknB Stp quiescence 2022`.
- Gap check: `Staphylococcus Stk1 Stp1 Pseudomonas supernatant`.

Screened primary evidence for species/strain,perturbation,experimental setting,
and whether results concern protein abundance,phosphorylation or a functional
phenotype. This was a focused background/gap check,not a systematic review or
proof that no similar study exists. The gap is therefore cautiously phrased as
incompletely resolved,without a first-discovery claim.

## Core background references

1. [Schwyter etal.,2025](https://pubmed.ncbi.nlm.nih.gov/40741853/): JE2 proteome
   responses to Pseudomonas supernatant. Closest exposure background; ancestor/
   evolved comparison differs from the current kinase/phosphatase mutants.
2. [Prust etal.,2021](https://pubmed.ncbi.nlm.nih.gov/33444734/): extensive
   staphylococcal phosphorylation and downstream Stk1-associated candidates.
   It does not validate every mutant-associated site as a direct substrate.
3. [Huemer etal.,2023](https://pubmed.ncbi.nlm.nih.gov/36595572/): phosphoregulation,
   stress-related metabolic/translation changes and quiescence/tolerance.
   Mainly acid-stressed CowanI with someJE2 follow-up; current experiment does
   not demonstrate antibiotic tolerance. Bibliographic year is2023 despite
   copyright2022 in the supplied PDF.

## Supporting context

4. [Niggli etal.,2023](https://pubmed.ncbi.nlm.nih.gov/37646506/): strain-specific
   adaptation to secreted inhibitory compounds. Evolutionary resistance is
   distinct from the current transposon-mutant protein comparison.
5. [Roman-Rodriguez etal.,2025](https://pubmed.ncbi.nlm.nih.gov/41114581/): direct
   evidence for HQNO-mediated respiratory inhibition and fermentative fitness.
   Uses USA300LAC/PA14 and pure-compound experiments,not the exact JE2/PAO1
   mixture. The peer-reviewed final article is used instead of its preprint.

Exact metadata was then retrieved through NCBI's PubMed EFetch API,checked
against all5expectedDOIs and written as a TSV and RIS for Zotero/EndNote.
`bibliography_API_provenance.tsv` preserves the request and retrieval timestamp.
The public MEDLINE response remains under ignored `data/raw/abstract_literature/`.
The first restricted-network attempt failed; after network permission was
granted,the same R script completed and verified all5records. No experimental
data were sent to PubMed.
