# Historical exploratory script; not an official final analysis.
# Read scripts/90_exploratory_archive/README.md before running.
source(if (file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
if(normalizePath(getwd(),winslash="/") != normalizePath(REPO_ROOT,winslash="/")) stop("Run archived scripts from the Sarya_George project folder.")
dir.create(file.path(REPO_ROOT,"results","exploratory"),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(REPO_ROOT,"figures","exploratory"),recursive=TRUE,showWarnings=FALSE)

source('R/setup.R')
out <- 'results/exploratory/G_inspect.txt'
sink(out, split = TRUE)
annot <- read_annotation()
cat('ANNOTATION SCHEMA\n'); print(names(annot))
targets <- find_protein(c('50S ribosomal protein L9','HtrA','trypsin-like serine protease'))
cat('CANDIDATE ANNOTATION\n'); print(targets, width=250)
sites <- read_sites('STYH')
cat('SITE SCHEMA\n'); print(names(sites))
cat('RplI candidate observed site identities\n')
print(sites |> filter(protein_key %in% targets$protein_key) |> distinct(site_id,locus_tag,product,residue,position,peptide,best_loc_prob), n=50, width=200)
dea <- read_dea('phospho_STYH', 'stp_vs_anc_at_TSB')
cat('DEA SCHEMA\n'); print(names(dea))
cat('CANDIDATE EXACT PUBLISHED POSITION\n')
preview <- sites |> filter(protein_key %in% targets$protein_key, (residue=='S' & position==92) | (residue=='T' & position==12))
print(preview |> distinct(site_id,locus_tag,product,residue,position,peptide,best_loc_prob), width=200)
if(nrow(preview)>0) {
 cat('OBSERVED REPLICATES AT EXACT PUBLISHED POSITION\n')
 print(preview |> group_by(site_id,genotype,medium) |> summarise(n_observed=sum(!is.na(intensity)),.groups='drop'),n=30)
 p <- plot_sites(unique(preview$site_id),sites)
 ggsave('figures/exploratory/G_known_site_preview.png',p,width=10,height=4.5,dpi=170)
}
cat('MEASURED PROTEIN CANDIDATE DETAILS\n'); print(targets |> select(protein_key,locus_tag,protein_id,product,staa8_uniprot,staa8_gene), width=250)
sink()
