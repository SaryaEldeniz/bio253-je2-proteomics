# Tests on the real data. Skipped if scripts/00_fetch_data.R was not run yet.

has_data <- function() file.exists(derived_path("dea_total.tsv.gz"))

test_that("the sign convention holds: diff of A_vs_B = mean(A) - mean(B)", {
  skip_if_not(has_data(), "data not fetched")
  tcyA <- "B7H15_13375"
  means <- read_abundances("total") |>
    filter(locus_tag == tcyA, group %in% c("anc_TSB", "anc_SN")) |>
    group_by(group) |>
    summarise(m = mean(normalized_abundance), n = sum(!is.na(normalized_abundance)))
  expect_equal(means$n, c(3L, 3L))                    # all observed, no imputation
  own_diff <- means$m[means$group == "anc_TSB"] - means$m[means$group == "anc_SN"]

  excel_diff <- read_dea("total", contrasts = "TSB_vs_SN_at_anc") |>
    filter(locus_tag == tcyA) |>
    pull(diff)
  expect_equal(own_diff, excel_diff, tolerance = 0.01)
  expect_lt(excel_diff, 0)                            # TcyA is higher in SN in the ancestor
})

test_that("abundance tables have 27 samples with the expected design", {
  skip_if_not(has_data(), "data not fetched")
  samples <- read_samples()
  expect_equal(nrow(samples), 27)
  expect_equal(as.vector(table(samples$group)), rep(3L, 9))
  ab <- read_abundances("total")
  expect_false(any(str_detect(ab$sample, "^2anc|^2stk|^2stp")))   # truncated names repaired
  expect_setequal(unique(ab$sample), samples$sample)
})

test_that("site tables: zeros became NA, residues and positions parsed", {
  skip_if_not(has_data(), "data not fetched")
  styh <- read_sites("STYH")
  expect_false(any(styh$intensity == 0, na.rm = TRUE))
  expect_setequal(unique(styh$residue), c("S", "T", "Y", "H"))
  expect_false(anyNA(styh$position))
  expect_true(all(is.na(styh$loc_prob[is.na(styh$intensity)])))
  sty <- read_sites("STY")
  expect_false("H" %in% sty$residue)
  expect_equal(n_distinct(styh$sample), 27)
})

test_that("find_protein finds the Tcy proteins by locus tag and protein id", {
  skip_if_not(has_data(), "data not fetched")
  expect_equal(nrow(find_protein("B7H15_13375")), 1)
  expect_equal(find_protein("ARG46976")$locus_tag, "B7H15_13375")
  expect_equal(nrow(find_protein(c("B7H15_13375", "B7H15_13370", "B7H15_13365"))), 3)
})
