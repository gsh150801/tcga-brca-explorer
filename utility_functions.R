split_query_str <- function(query_str) {
  ids <- unlist(strsplit(query_str, split = " "))
}

retrieve_tcga_data <- function(conn, ids) {
  mutations_df <- getProfileData(conn,
    genes = ids,
    geneticProfiles = "brca_tcga_mutations",
    caseList = "brca_tcga_all"
  )
  names(mutations_df) <- paste0(names(mutations_df), "_mutations")
  mutations_df[] <- lapply(mutations_df, as.character)
  mutations_df[is.na(mutations_df)] <- "(germline)"
  mutations_df[] <- lapply(
    mutations_df,
    function(x) {
      x[x %in% "NaN"] <- NA
      x
    }
  )
  mutations_df[] <- lapply(
    mutations_df,
    function(x) relevel(factor(x), ref = "(germline)")
  )
  mutations_df$subjid <- row.names(mutations_df)

  gistic_df <- getProfileData(conn,
    genes = ids,
    geneticProfiles = "brca_tcga_gistic",
    caseList = "brca_tcga_all"
  )
  names(gistic_df) <- paste0(names(gistic_df), "_gistic")
  gistic_df[] <- lapply(gistic_df, function(x) factor(x,
      levels = -2:2,
      labels = c("Deep depletion", "Shallow depletion", "Diploid", "Gain", "Amplification")
    ))
  gistic_df$subjid <- row.names(gistic_df)

  rna_df <- getProfileData(conn,
    genes = ids,
    geneticProfiles = "brca_tcga_rna_seq_v2_mrna",
    caseList = "brca_tcga_all"
  )

  retrieved_ids <- names(rna_df)
  retrieved_ids <- ids[ids %in% retrieved_ids] # keep original order

  names(rna_df) <- paste0(names(rna_df), "_rna")
  rna_df[] <- lapply(rna_df, function(x) log2(x + 1))
  rna_df$subjid <- row.names(rna_df)

  # clin_df <- getClinicalData(conn,
  #   caseList = "brca_tcga_all")
  # clin_df$subjid <- row.names(clin_df)
  # names(clin_df) <- tolower(names(clin_df))

  tcga_df <- Reduce(
    function(x, y) merge(x, y, by = "subjid", all = TRUE),
    list(mutations_df, gistic_df, rna_df)
  )

  return(
    list(
      "ids" = retrieved_ids,
      "data" = tcga_df
    )
  )
}

perform_subtype_classification <- function(conn, pam50centroids) {
  exprs <- getProfileData(conn,
    genes = row.names(pam50centroids),
    geneticProfiles = "brca_tcga_rna_seq_v2_mrna",
    caseList = "brca_tcga_rna_seq_v2_mrna"
  )
  exprs[] <- lapply(exprs, function(x) log2(x + 1))
  exprs <- t(as.matrix(exprs))

  ids <- intersect(row.names(pam50centroids), row.names(exprs))
  rel_exprs <- sweep(exprs[ids, ], 1, apply(exprs[ids, ], 1, "median", na.rm = TRUE))
  correlations <- cor(rel_exprs, pam50centroids[ids, ], method = "spearman")
  subtype_idx <- unlist(apply(correlations[colnames(exprs), ], 1, which.max))
  subtypecd <- colnames(pam50centroids)[subtype_idx]

  subtype_lkup <- c(
    "LA" = "Luminal A",
    "LB" = "Luminal B",
    "H2" = "HER2-enriched",
    "BL" = "Basal-like",
    "NBL" = "Normal breast-like"
  )

  data.frame(
    row.names = colnames(exprs),
    subjid = colnames(exprs),
    subtypecd = factor(subtypecd,
      levels = names(subtype_lkup)
    ),
    subtype = factor(subtype_lkup[subtypecd],
      levels = subtype_lkup
    ),
    subtype2 = factor(subtype_lkup[subtypecd],
      levels = subtype_lkup[c("H2", "BL", "NBL", "LA", "LB")]
    ),
    stringsAsFactors = FALSE
  )
}
