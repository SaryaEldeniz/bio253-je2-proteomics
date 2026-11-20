# Tests for the parsing helpers. These do NOT need the data.

test_that("sample names are parsed, truncated names are repaired", {
  x <- parse_sample_name(c("2anc_pulseSN_1",
                           "JE2stk_TSB_3",
                           "20260612_053_rep_C41881_S1143466_JE2stk_pulseSN_1_Enriched",
                           "20260610_027_C41881_S1143474_JE2stp_SN_3_Enriched Intensity"))
  expect_equal(x$sample, c("JE2anc_pulseSN_1", "JE2stk_TSB_3", "JE2stk_pulseSN_1", "JE2stp_SN_3"))
  expect_equal(as.character(x$genotype), c("anc", "stk", "stk", "stp"))
  expect_equal(as.character(x$medium), c("SNpulse", "TSB", "SNpulse", "SN"))
  expect_equal(x$replicate, c(1L, 3L, 1L, 3L))
  expect_equal(as.character(x$group), c("anc_SNpulse", "stk_TSB", "stk_SNpulse", "stp_SN"))
  expect_equal(x$reinjected, c(FALSE, FALSE, TRUE, FALSE))
})

test_that("unparseable names give NA instead of nonsense", {
  x <- parse_sample_name("QC_standard_01")
  expect_true(is.na(x$sample))
  expect_true(is.na(x$genotype))
})

test_that("FASTA descriptions are parsed", {
  d <- paste("FGCZ [locus_tag=B7H15_13375] [protein=amino acid ABC transporter",
             "substrate-binding protein] [protein_id=ARG46976.1]",
             "[location=complement(2535968..2536747)] [gbkey=CDS]",
             "BLASTORTHO tr_Q2FVL2_Q2FVL2_STAA8 evalue=0")
  p <- parse_description(d)
  expect_equal(p$locus_tag, "B7H15_13375")
  expect_equal(p$protein_id, "ARG46976.1")
  expect_equal(p$product, "amino acid ABC transporter substrate-binding protein")
  expect_equal(p$staa8_uniprot, "Q2FVL2")
  expect_equal(p$strand, "-")
  expect_true(is.na(p$staa8_gene))           # TrEMBL entry: no gene name

  sp <- parse_description("FGCZ [locus_tag=B7H15_13880] [protein=surface protein] [pseudo=true] BLASTORTHO sp_Q2G2B2_SASG_STAA8 evalue=0")
  expect_equal(sp$staa8_gene, "SASG")
  expect_true(sp$pseudo)

  fs <- parse_description("FGCZ ID=cds-B7H15_10525; ID=cds-B7H15_10525;Parent=gene-B7H15_10525;Note=frameshifted")
  expect_equal(fs$locus_tag, "B7H15_10525")
  expect_true(fs$pseudo)
})

test_that("protein keys are extracted from tr|KEY|KEY ids", {
  expect_equal(protein_key(c("tr|lcl_CP020619.1_prot_ARG46976.1_2593|lcl_CP020619.1_prot_ARG46976.1_2593",
                             "lcl_CP020619.1_prot_2694")),
               c("lcl_CP020619.1_prot_ARG46976.1_2593", "lcl_CP020619.1_prot_2694"))
})

test_that("explain_contrast states the A - B rule", {
  expect_match(explain_contrast("TSB_vs_SN_at_anc"), "log2\\(TSB\\) - log2\\(SN\\), within anc")
  expect_match(explain_contrast("TSB_vs_SN_at_anc"), "diff > 0 means HIGHER in TSB")
  expect_match(explain_contrast("stk_vs_anc"), "averaged over the other factor")
  expect_match(explain_contrast("interaction_stk_vs_anc_at_TSB_vs_SN"),
               "\\(stk - anc in TSB\\) - \\(stk - anc in SN\\)")
})

test_that("pHis peptide helpers handle ';' variants and lower case", {
  expect_equal(n_sty_acceptors("ANhEQIIEAIK;MANhEQIIEAIK"), 0)
  expect_equal(n_sty_acceptors("MEQNSYVIIDETGIhAR"), 3)
  expect_equal(n_sty_acceptors("AShK;HK"), 0)          # fewest over variants
  expect_true(co_phospho_sty("MEQNsYVIIDETGIhAR"))
  set <- c("MEQNSYVIIDETGIHAR")
  expect_true(all_variants_in("MEQNSYVIIDETGIhAR", set))
  expect_false(all_variants_in("MEQNSYVIIDETGIhAR;MEQNSYVIIDETGIHARPATMLVQTASK", set))
  expect_true(any_variant_in("MEQNSYVIIDETGIhAR;MEQNSYVIIDETGIHARPATMLVQTASK", set))
})
