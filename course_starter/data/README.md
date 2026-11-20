# data/ – local only, never committed

The data are **unpublished**. Everything in this folder except this README, the
checksum file and the two check scripts is git-ignored. Do not remove those lines from
`.gitignore`: the data files never go to GitHub (they are unpublished, and too large
anyway). Your code and rendered reports do. Any AI tool may be used to analyse them.

## How to get the data

The data are **not** in this repository. The lecturer distributes them as **one zip file**.

1. **Extract** the zip into this `data/` folder. Afterwards you have `data/raw/`,
   `data/derived/` and `data/reference/` next to this README (see *What you get* below).
2. **Check the files** with the MD5 checksums in [`checksums.md5`](checksums.md5)
   (run from the repository root):

   ```bash
   bash data/check_data.sh          # macOS / Linux
   Rscript data/check_data.R        # Windows, or anywhere R is installed
   ```

   Every file should say `OK`, and the last line `Data check passed.`
   * `MISSING`: a file is not where it should be. Did you extract into `data/`, and not into
     `data/<zip name>/`?
   * `CORRUPT`: the file differs from the lecturer's version (incomplete download, broken
     extraction, or the file was edited). Download and extract again.

   Experts can also run `cd data && md5sum -c checksums.md5` directly.

That's it: the tidy tables in `data/derived/` are already included, and you can start with
`analyses/01_explore_data.qmd`. If something looks odd later, run the check again: if it
passes, the problem is in the code, not in the data. Never edit files in `data/`; write
your own results elsewhere (e.g. `results/`).

*For the lecturer:* the zip is built from this folder with
`zip -r bio253_je2_data.zip raw derived reference -x '*.DS_Store'` (run inside `data/`).
If any data file changes, regenerate the checksums with
`find raw derived reference -type f ! -name .DS_Store | LC_ALL=C sort | xargs md5sum > checksums.md5`
and commit the new `checksums.md5`. `scripts/00_fetch_data.R` rebuilds `raw/` and `derived/`
from the original course folder.

## What you get

```
data/raw/
  MANIFEST.tsv                       source path, size and md5 of every copied file
  total/DE_total.xlsx                prolfquapp DEA, total proteome (Oct 2026 run)
  phospho_STYH/DE_phospho_STYH.xlsx  prolfquapp DEA, phospho-sites, S/T/Y/H search (Oct 2026)
  phospho_STY/DE_phospho_STY.xlsx    prolfquapp DEA, phospho-sites, S/T/Y search (Aug 2026, older model)
  sites/combined_site_STY.tsv        FragPipe site tables (raw, wide)
  sites/combined_site_STYH.tsv
  resources/JE2_proteins.fasta       search database (JE2 CP020619.1 + pseudogene translations + contaminants)
  resources/JE2_genome.gff           genome annotation
  resources/string_*.txt             STRING v12 aliases, enrichment terms, physical links for JE2
data/derived/                        tidy long tables (tab-separated, gzipped)
  protein_annotation.tsv.gz          locus_tag, product, protein_id, location, NCTC8325 ortholog
  samples.tsv                        27 samples: genotype, medium, replicate, group
  dea_<dataset>.tsv.gz               DEA results, one row per protein/site x contrast
  abundances_<dataset>.tsv.gz        normalised abundances, one row per protein/site x sample
  sites_STY.tsv.gz, sites_STYH.tsv.gz  FragPipe sites, one row per site x sample (0 -> NA)
data/reference/                      lecturer reference tables, if available (e.g. pHis tiers)
```

`<dataset>` is `total`, `phospho_STYH` or `phospho_STY`. Load the tables with the
helpers (`read_dea()`, `read_abundances()`, `read_sites()`, …) rather than by hand.
The helpers attach the protein annotation and set the factor levels.

## Things to know

* Contrast `A_vs_B`: `diff = log2(A) − log2(B)`.
* `estimate_type = lod_imputed`: a group mean was imputed, not measured.
* Truncated total-proteome sample names (`2anc_pulseSN_1`) are repaired to `JE2anc_pulseSN_1`.
* The phospho DEAs are site-level. `site` = `<protein>_<residue><position>~<peptide>`.
* More material in the course folder, not copied here: the HTML QC/DEA reports,
  `SummarizedExperiment.rds`, `AnnData.h5ad` and `*.parquet` in each `Results_*` folder,
  plus ORA/GSEA input files.
