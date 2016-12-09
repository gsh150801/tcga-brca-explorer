retreive_data <- function(ids) {
  conn <- CGDS("http://www.cbioportal.org/public-portal/")
  
  mutations_df <- getProfileData(conn,
    genes = ids,
    geneticProfiles = "brca_tcga_mutations",
    caseList = "brca_tcga_all")
  names(mutations_df) <- paste0(names(mutations_df), "_mutations")
  mutations_df[] <- lapply(mutations_df, as.character)
  mutations_df[is.na(mutations_df)] <- "(none)"
  mutations_df[] <- lapply(mutations_df, 
    function(x) {
      x[x %in% "NaN"] <- NA
      x})
  mutations_df[] <- lapply(mutations_df, 
    function(x) relevel(factor(x), ref = "(none)"))

  gistic_df <- getProfileData(conn,
    genes = ids,
    geneticProfiles = "brca_tcga_gistic",
    caseList = "brca_tcga_all")
  names(gistic_df) <- paste0(names(gistic_df), "_gistic")
  gistic_df[] <- lapply(gistic_df, function(x) factor(x, 
    levels = -2:2, 
    labels = c("Deep depletion", "Shallow depletion", "Diploid", "Gain", "Amplification")))
  
  rna_df <- getProfileData(conn,
    genes = ids,
    geneticProfiles = "brca_tcga_rna_seq_v2_mrna",
    caseList = "brca_tcga_all")
  names(rna_df) <- paste0(names(rna_df), "_rna")
  rna_df[] <- lapply(rna_df, function(x) log2(x + 1))
  
  tcga_df <- cbind(mutations_df, gistic_df, rna_df) %>%
    identity()
  
  return(tcga_df)
}
