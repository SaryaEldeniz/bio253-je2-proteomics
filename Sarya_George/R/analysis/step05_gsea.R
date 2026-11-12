# Preranked GSEA via the official fgsea adaptive multilevel implementation.
# Main Step 5 code supplies unique target IDs and defines multiplicity families.

run_step05_gsea <- function(stats, pathways, seed = 20261112L,
                            min_size = 10L, max_size = 500L, eps = 1e-10,
                            sample_size = 101L, n_perm_simple = 1000L) {
  if (!is.numeric(stats) || length(stats) < 2L || is.null(names(stats)) ||
      anyNA(names(stats)) || any(!nzchar(names(stats))) || anyDuplicated(names(stats)) ||
      any(!is.finite(stats))) {
    stop("GSEA requires a finite numeric score vector with unique, nonempty feature-ID names.")
  }
  if (all(stats == 0)) stop("GSEA requires at least one nonzero ranking score.")
  if (!is.list(pathways) || is.null(names(pathways)) || anyNA(names(pathways)) ||
      any(!nzchar(names(pathways))) || anyDuplicated(names(pathways)) ||
      any(!vapply(pathways, function(x) is.character(x) && !anyNA(x) && all(nzchar(x)), logical(1)))) {
    stop("Supply a named list of character gene sets with unique, nonempty pathway names.")
  }
  is_integer <- function(x) length(x) == 1L && is.finite(x) && x == floor(x)
  if (!is_integer(min_size) || min_size < 1L || !is_integer(max_size) || max_size < min_size) {
    stop("Pathway size limits must be positive integers with max_size >= min_size.")
  }
  if (!is_integer(seed) || seed < 0 || seed > .Machine$integer.max) stop("Use one valid integer seed.")
  if (length(eps) != 1L || !is.finite(eps) || eps < 0 || eps >= 1) stop("eps must be in [0, 1).")
  if (!is_integer(sample_size) || sample_size < 3L || sample_size %% 2L != 1L) {
    stop("sample_size must be an odd integer >=3.")
  }
  if (!is_integer(n_perm_simple) || n_perm_simple < 1000L) stop("n_perm_simple must be an integer >=1000.")
  # Tied ranks are ordered deterministically by identifier, with no added jitter.
  stats <- stats[order(-stats, names(stats), method = "radix")]
  intersected <- lapply(pathways, function(x) intersect(unique(x), names(stats)))
  sizes <- lengths(intersected)
  eligible <- sizes >= min_size & sizes <= max_size & sizes < length(stats)
  prepared <- intersected[eligible]
  method <- list(engine = "fgsea::fgseaMultilevel", seed = as.integer(seed),
                 min_size = as.integer(min_size), max_size = as.integer(max_size),
                 eps = eps, sample_size = as.integer(sample_size),
                 n_perm_simple = as.integer(n_perm_simple), gseaParam = 1,
                 scoreType = "std", backend = "BiocParallel::SerialParam",
                 n_ranked_features = length(stats), n_pathways_requested = length(pathways),
                 n_pathways_eligible = length(prepared),
                 tied_score_fraction = sum(duplicated(stats)) / length(stats),
                 tie_policy = "descending score; ID order for ties; no jitter",
                 score_universe = "all supplied eligible ranked features")
  if (!length(prepared)) {
    result <- tibble::tibble(pathway = character(), ES = double(), NES = double(),
                             pval = double(), padj_engine = double(), size = integer(),
                             leadingEdge = list(), log2err = double(),
                             status = character(), pval_at_floor = logical())
    attr(result, "gsea_method") <- method
    attr(result, "gsea_warnings") <- "No pathways met the declared effective-size limits."
    attr(result, "elapsed_seconds") <- 0
    return(result)
  }
  if (!requireNamespace("fgsea", quietly = TRUE) || !requireNamespace("BiocParallel", quietly = TRUE)) {
    stop("Install compatible Bioconductor packages fgsea and BiocParallel before running GSEA; no substitute test is used.")
  }
  method$fgsea_version <- tryCatch(as.character(utils::packageVersion("fgsea")), error = function(e) NA_character_)
  saved_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  saved_rng <- RNGkind()
  on.exit({
    do.call(RNGkind, as.list(saved_rng))
    if (is.null(saved_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", saved_seed, envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(as.integer(seed))
  messages <- character()
  start <- proc.time()[["elapsed"]]
  raw <- withCallingHandlers(
    fgsea::fgseaMultilevel(pathways = prepared, stats = stats,
                           minSize = as.integer(min_size), maxSize = as.integer(max_size),
                           eps = eps, sampleSize = as.integer(sample_size),
                           nPermSimple = as.integer(n_perm_simple), scoreType = "std",
                           gseaParam = 1, BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  elapsed <- proc.time()[["elapsed"]] - start
  result <- tibble::as_tibble(raw)
  attr(result, ".internal.selfref") <- NULL
  if (!"log2err" %in% names(result)) result$log2err <- NA_real_
  result <- result |>
    dplyr::rename(padj_engine = padj) |>
    dplyr::mutate(pval_at_floor = is.finite(pval) & eps > 0 & pval <= eps,
                  status = dplyr::case_when(!is.finite(pval) | !is.finite(NES) ~ "not_estimable_review_engine_warning",
                                            pval_at_floor ~ "estimated_p_at_declared_floor",
                                            !is.finite(log2err) ~ "estimated_uncertainty_not_available_review_warning",
                                            TRUE ~ "estimated")) |>
    dplyr::select(pathway, ES, NES, pval, padj_engine, size, leadingEdge, log2err, status, pval_at_floor) |>
    dplyr::arrange(pathway)
  # Do not replace failed/NA estimates with invented P-values.
  if (length(messages)) for (item in unique(messages)) message("fgsea warning retained: ", item)
  attr(result, "gsea_method") <- method
  attr(result, "gsea_warnings") <- unique(messages)
  attr(result, "elapsed_seconds") <- unname(elapsed)
  result
}
