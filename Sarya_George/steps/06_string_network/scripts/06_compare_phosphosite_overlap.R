# Exact protein overlap; keep total-protein interactions and within-medium site
# effects separate. Neither STRING edges nor overlap establish direct substrates.
source(if(file.exists("Sarya_George.Rproj")) "R/setup.R" else "Sarya_George/R/setup.R")
out <- file.path(REPO_ROOT,"steps/06_string_network/reports/phosphosite_overlap")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
key <- tribble(
  ~uniprot_acc,~locus_tag,~protein,
  "Q2FV52","B7H15_14290","IsaA",
  "Q2FZY9","B7H15_04690","CsbD family",
  "Q2G227","B7H15_00775","Phosphopentomutase",
  "Q2FZK0","B7H15_05485","QoxB",
  "Q2FXM9","B7H15_09250","Pyk",
  "Q2FXL6","B7H15_09315","UspA",
  "Q2FWE5","B7H15_11735","GlyA",
  "Q2FZH6","B7H15_05615","HPr",
  "Q2FX13","B7H15_10475","Hypothetical protein",
  "Q2FZD2","B7H15_05950","Thioredoxin",
  "Q2FYI5","B7H15_07625","GpsB",
  "Q2FZ28","B7H15_06530","HslU")
# Transcribed exactly from the student's table; CsbD S5 occurs in both media.
sites <- tribble(
  ~uniprot_acc,~residue,~position,~medium,~mutant,
  "Q2FV52","S",194,"TSB","stk",
  "Q2FZY9","S",5,"TSB","stp",
  "Q2FZY9","S",5,"SN","stp",
  "Q2FZY9","T",52,"SN","stp",
  "Q2G227","T",85,"TSB","stk",
  "Q2FZK0","T",552,"TSB","stp",
  "Q2FXM9","T",537,"TSB","stp",
  "Q2FXL6","T",157,"TSB","stp",
  "Q2FWE5","Y",51,"TSB","stk",
  "Q2FZH6","S",52,"SN","stk",
  "Q2FX13","S",51,"SN","stk",
  "Q2FZD2","S",71,"SN","stp",
  "Q2FYI5","T",93,"SN","stk",
  "Q2FZ28","T",217,"SN","stp")
annotation <- read_tsv(file.path(DATA_ROOT,"derived/protein_annotation.tsv.gz"),show_col_types=FALSE)
identity <- key |> left_join(annotation |> select(protein_key,locus_tag,staa8_uniprot),
  by="locus_tag",relationship="one-to-one")
stopifnot(nrow(identity)==12L,!anyNA(identity$protein_key),
  all(identity$uniprot_acc==identity$staa8_uniprot))
nodes_dir <- file.path(REPO_ROOT,"steps/06_string_network/reports/poster_networks")
nodes <- read_csv(file.path(nodes_dir,"nodes_75_live.csv"),show_col_types=FALSE)
comparison <- identity
for(set in c("49","67","75")) {
  n <- read_csv(file.path(nodes_dir,paste0("nodes_",set,"_live.csv")),show_col_types=FALSE)
  comparison[[paste0("in_",set)]] <- comparison$uniprot_acc %in% n$uniprot_acc
}
stopifnot(sum(comparison$in_49)==2L,sum(comparison$in_67)==1L,sum(comparison$in_75)==3L,
  setequal(comparison$uniprot_acc[comparison$in_75],c("Q2FV52","Q2FXM9","Q2FZK0")))
comparison <- comparison |> left_join(nodes |> select(uniprot_acc,network_protein_key=protein_key,
  network_locus=locus_tag,string_id,degree,module_group),by="uniprot_acc",relationship="one-to-one")
stopifnot(all(comparison$protein_key[comparison$in_75]==comparison$network_protein_key[comparison$in_75]),
  all(comparison$locus_tag[comparison$in_75]==comparison$network_locus[comparison$in_75]))
write_csv(comparison,file.path(out,"all12_protein_overlap.csv"))
overlap <- filter(comparison,in_75)

phospho <- read_tsv(file.path(DATA_ROOT,"derived/dea_phospho_STYH.tsv.gz"),show_col_types=FALSE)
total <- read_tsv(file.path(DATA_ROOT,"derived/dea_total.tsv.gz"),show_col_types=FALSE)
site_hits <- sites |> left_join(identity |> select(uniprot_acc,protein_key,protein,locus_tag),
  by="uniprot_acc",relationship="many-to-one") |>
  mutate(contrast=paste0(mutant,"_vs_anc_at_",medium)) |>
  left_join(phospho |> select(protein_key,residue,position,contrast,site_id,
      site_diff=diff,site_FDR=FDR,site_estimate_type=estimate_type,Peptide,SequenceWindow),
    by=c("protein_key","residue","position","contrast"),relationship="one-to-one")
stopifnot(nrow(site_hits)==14L,all(site_hits$site_estimate_type=="observed"),
  all(site_hits$site_FDR<.1),all(abs(site_hits$site_diff)>1))
# 'Exclusive' describes one thresholded mutant-vs-ancestor set, not detection
# exclusively in that mutant and not proof that two mutant effects differ.
other_effect <- phospho |> transmute(protein_key,residue,position,
  other_contrast=contrast,other_diff=diff,other_FDR=FDR,other_estimate_type=estimate_type)
site_hits <- site_hits |> mutate(other_contrast=paste0(if_else(mutant=="stk","stp","stk"),"_vs_anc_at_",medium)) |>
  left_join(other_effect,by=c("protein_key","residue","position","other_contrast"),relationship="one-to-one") |>
  mutate(other_meets_selection=other_estimate_type=="observed" & other_FDR<.1 & abs(other_diff)>1,
    in_75=uniprot_acc %in% overlap$uniprot_acc)
stopifnot(!any(site_hits$other_meets_selection))
write_csv(site_hits,file.path(out,"all14_site_hits_checked.csv"))
joined <- site_hits |> filter(in_75) |>
  left_join(total |> select(protein_key,contrast,parent_diff=diff,parent_FDR=FDR,
    parent_estimate_type=estimate_type),by=c("protein_key","contrast"),relationship="one-to-one") |>
  mutate(site_direction=if_else(site_diff>0,"Higher phosphosite signal in mutant than ancestor",
      "Lower phosphosite signal in mutant than ancestor"),
    parent_direction=if_else(parent_diff>0,"Higher parent abundance in mutant than ancestor",
      "Lower parent abundance in mutant than ancestor"))
write_csv(joined,file.path(out,"three_overlapping_sites_and_parent_TSB.csv"))

codes <- c("interaction_stk_vs_anc_at_TSB_vs_SN","interaction_stp_vs_anc_at_TSB_vs_SN",
  "interaction_stp_vs_stk_at_TSB_vs_SN")
site_interactions <- phospho |> semi_join(joined,by=c("protein_key","residue","position")) |>
  filter(contrast %in% codes) |> left_join(overlap |> select(protein_key,protein,uniprot_acc),
    by="protein_key",relationship="many-to-one")
write_csv(site_interactions,file.path(out,"three_sites_PASN_interaction_estimates.csv"))
protein_interactions <- nodes |> semi_join(overlap,by="protein_key")
write_csv(protein_interactions,file.path(out,"three_proteins_original_interaction_estimates.csv"))
responses <- total |> semi_join(overlap,by="protein_key") |>
  filter(contrast %in% c("TSB_vs_SN_at_anc","TSB_vs_SN_at_stk","TSB_vs_SN_at_stp")) |>
  left_join(overlap |> select(protein_key,protein,uniprot_acc),by="protein_key",relationship="many-to-one") |>
  mutate(SN_minus_TSB=-diff,SN_response_direction=if_else(SN_minus_TSB>0,"Higher in SN than TSB","Lower in SN than TSB"))
stopifnot(nrow(responses)==9L,all(responses$estimate_type=="observed"))
write_csv(responses,file.path(out,"three_protein_SN_responses.csv"))

edges <- read_csv(file.path(nodes_dir,"edges_75.csv"),show_col_types=FALSE) |>
  filter(node_a %in% overlap$string_id | node_b %in% overlap$string_id) |>
  left_join(nodes |> select(node_a=string_id,name_a=display_name,uniprot_a=uniprot_acc),by="node_a",relationship="many-to-one") |>
  left_join(nodes |> select(node_b=string_id,name_b=display_name,uniprot_b=uniprot_acc),by="node_b",relationship="many-to-one")
write_csv(edges,file.path(out,"overlap_protein_network_edges.csv"))

# Verify JE2 site coordinates against JE2 FASTA; a reference accession alone
# does not validate site numbering across strains.
fasta <- readLines(file.path(DATA_ROOT,"raw/resources/JE2_proteins.fasta"),warn=FALSE)
h <- which(startsWith(fasta,">")); ends <- c(h[-1]-1L,length(fasta))
keys <- sub("^>tr\\|([^|]+)\\|.*","\\1",fasta[h])
seqs <- vapply(match(joined$protein_key,keys),function(i) paste0(fasta[(h[i]+1L):ends[i]],collapse=""),character(1))
stopifnot(!anyNA(seqs),all(substr(seqs,joined$position,joined$position)==joined$residue))
sequence_check <- joined |> select(uniprot_acc,protein,protein_key,residue,position) |>
  mutate(JE2_length=nchar(seqs),JE2_window=vapply(seq_along(seqs),function(i)
    substr(seqs[i],max(1,position[i]-11),min(nchar(seqs[i]),position[i]+10)),character(1)),
    residue_verified=TRUE)
write_csv(sequence_check,file.path(out,"JE2_site_coordinate_checks.csv"))

abund_total <- read_tsv(file.path(DATA_ROOT,"derived/abundances_total.tsv.gz"),show_col_types=FALSE) |>
  semi_join(overlap,by="protein_key") |> filter(medium %in% c("TSB","SN")) |>
  left_join(overlap |> select(protein_key,protein),by="protein_key",relationship="many-to-one") |>
  mutate(measurement="Total protein",site_label="parent")
abund_sites <- read_tsv(file.path(DATA_ROOT,"derived/abundances_phospho_STYH.tsv.gz"),show_col_types=FALSE) |>
  semi_join(joined,by="site_id") |> filter(medium %in% c("TSB","SN")) |>
  left_join(overlap |> select(protein_key,protein),by="protein_key",relationship="many-to-one") |>
  mutate(measurement="Phosphosite",site_label=paste0(residue,position))
sample_values <- bind_rows(abund_total,abund_sites) |>
  select(protein,protein_key,measurement,site_label,sample,genotype,medium,replicate,normalized_abundance,
    any_of("Localization_Probability"))
# Save only small summaries to reports, not the individual source measurements.
detection <- sample_values |> group_by(protein,measurement,site_label,genotype,medium) |>
  summarise(n_observed=sum(is.finite(normalized_abundance)),
    mean_log2=if(any(is.finite(normalized_abundance))) mean(normalized_abundance,na.rm=TRUE) else NA_real_,
    .groups="drop")
write_csv(detection,file.path(out,"sample_detection_and_means.csv"))
loc <- abund_sites |> group_by(protein,site_label) |>
  summarise(min_exported_localization=min(Localization_Probability,na.rm=TRUE),
    max_exported_localization=max(Localization_Probability,na.rm=TRUE),
    min_localization_observed=min(Localization_Probability[is.finite(normalized_abundance)],na.rm=TRUE),
    max_localization_observed=max(Localization_Probability[is.finite(normalized_abundance)],na.rm=TRUE),.groups="drop")
raw_site_summary <- read_tsv(file.path(DATA_ROOT,"derived/sites_STYH.tsv.gz"),show_col_types=FALSE) |>
  semi_join(joined,by="site_id") |>
  left_join(overlap |> select(protein_key,protein),by="protein_key",relationship="many-to-one") |>
  group_by(protein) |> summarise(best_loc_prob=min(best_loc_prob,na.rm=TRUE),.groups="drop")
loc <- loc |> left_join(raw_site_summary,by="protein",relationship="one-to-one")
write_csv(loc,file.path(out,"site_localization_summary.csv"))
loc_groups <- abund_sites |> group_by(protein,site_label,genotype,medium) |>
  summarise(n_observed=sum(is.finite(normalized_abundance)),
    min_localization_observed=min(Localization_Probability[is.finite(normalized_abundance)],na.rm=TRUE),
    max_localization_observed=max(Localization_Probability[is.finite(normalized_abundance)],na.rm=TRUE),.groups="drop")
write_csv(loc_groups,file.path(out,"site_localization_by_group.csv"))

# Plot the actual observed values, retaining empty groups as missing. No imputation.
sample_values <- sample_values |> mutate(genotype=factor(genotype,levels=c("anc","stk","stp")),
  medium=factor(medium,levels=c("TSB","SN")),
  panel=paste0(protein," | ",if_else(measurement=="Total protein","total protein",site_label)))
panels <- c(paste0(c("IsaA","Pyk","QoxB")," | total protein"),"IsaA | S194","Pyk | T537","QoxB | T552")
sample_values$panel <- factor(sample_values$panel,levels=panels)
detection$panel <- factor(paste0(detection$protein," | ",if_else(detection$measurement=="Total protein","total protein",detection$site_label)),levels=panels)
plot <- ggplot(sample_values,aes(medium,normalized_abundance,colour=genotype)) +
  geom_point(position=position_jitterdodge(jitter.width=.08,dodge.width=.65,seed=19),size=2.1,na.rm=TRUE) +
  stat_summary(fun=mean,geom="point",shape=95,size=7,position=position_dodge(width=.65),na.rm=TRUE) +
  geom_text(data=detection,aes(y=-Inf,label=paste0(n_observed,"/3")),
    position=position_dodge(width=.65),vjust=-.4,size=3,show.legend=FALSE) +
  facet_wrap(~panel,ncol=3,scales="free_y",drop=FALSE) +
  scale_colour_manual(values=c(anc="#333333",stk="#0072B2",stp="#D55E00"),drop=FALSE) +
  scale_y_continuous(expand=expansion(mult=c(.18,.08))) +
  labs(x=NULL,y="Source normalized log2 abundance",colour="Genotype",
    title="Proteins shared by the 75-node network and phosphosite list",
    subtitle="Top: total protein. Bottom: phosphosite signal; panels have separate y scales.",
    caption="Points are observed sample values; bars are means; n/3 labels show detection. Missing values remain NA.\nPhosphosite signal is not phosphorylation occupancy; these assays are not on a common absolute scale.") +
  theme_bw(base_size=11) + theme(legend.position="top",strip.background=element_rect(fill="#F3F5F7"),
    plot.caption=element_text(hjust=0),panel.grid.minor=element_blank())
fig <- file.path(REPO_ROOT,"steps/06_string_network/figures/phosphosite_overlap")
dir.create(fig,recursive=TRUE,showWarnings=FALSE)
ggsave(file.path(fig,"three_overlap_sample_values.png"),plot,width=11,height=7,dpi=240,bg="white")
writeLines(trimws(capture.output(sessionInfo()),which="right"),file.path(out,"sessionInfo.txt"))
writeLines(c("Exact-UniProt and JE2-locus matches:49=2,67=1,75=3.",
  "All14 user-listed site hits reproduce the original observed FDR<0.1 and |effect|>1 exclusive selection.",
  "Allthree overlapping site coordinates verified against JE2 FASTA.",
  "Source effects and FDR retained; no new model or overlap-enrichment test.",
  "No direct Stk1/Stp1 substrate claim follows from protein overlap."),file.path(out,"checks.txt"))
print(joined |> select(protein,residue,position,contrast,site_diff,site_FDR,site_estimate_type,
  parent_diff,parent_FDR,parent_estimate_type),width=160)
print(site_interactions |> select(protein,contrast,diff,FDR,estimate_type),n=9,width=160)
print(loc)
message("Overlap audit written:3/12 proteins in75;49hadIsaA/Pyk;67hadQoxB.")
