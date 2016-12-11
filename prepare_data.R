library(Biobase)
library(org.Hs.eg.db)

data("pam50", package = "genefu")
pam50centroids <- pam50$centroids

rownames(pam50centroids) <- unlist(mget(
  as.character(pam50$centroids.map$EntrezGene.ID), org.Hs.egSYMBOL, 
  ifnotfound = NA))

lkup <- c(
  "LumA" = "LA", 
  "LumB" = "LB", 
  "Her2" = "H2", 
  "Basal" = "BL", 
  "Normal" = "NBL")
colnames(pam50centroids) <- lkup[colnames(pam50$centroids)]

save(pam50centroids, file = file.path("data", "pam50centroids.rda"))
