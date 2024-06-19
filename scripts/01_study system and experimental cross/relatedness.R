# script written by andy lee on 06/19/2024; contact at andymuanlee@gmail.com
# Rversion:

library(reshape2)
library(snpR)

setwd("~/KW/09_angsd/relatedness/relatedness_individuals/")


library(related)


### calculate eculidan distances between samples 
sample_meta <- read.csv("~/KW/09_angsd/neutral/sample_info.csv")
all_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_allsnps.vcf", sample.meta = sample_meta)

# all_monnap <- filter_snps(all_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")
all_mon <- filter_snps(all_vcf[site=c("MON")], maf = 0.05, maf_facets = "site")
all_nap <- filter_snps(all_vcf[site=c("NAP")], maf = 0.05, maf_facets = "site")
pca_mon <- plot_clusters(all_mon, "site", simplify_output = TRUE,)$pca$data
pca_nap <- plot_clusters(all_nap, "site", simplify_output = TRUE,)$pca$data

## write a function to compare within and among family euclidan distance on a PCA 
euclid.dist <- function(pca){
  dists <- dist(method = "euclidian", pca[,c("PC1", "PC2")])
  dists <- as.data.frame(as.matrix(dists))
  dists2 <- cbind(dists, sample_ID = pca$.sample.id, family = pca$family)
  dists.m <- reshape2::melt(dists2, id.vars = c("sample_ID", "family"))
  dists.m$variable <- as.numeric(as.character(dists.m$variable))
  
  # dists.m2 <- cbind(dists.m, sample_2_ID = dists.m$variable)
  
  dists.m2 <- merge(dists.m, pca, by.x = "variable", by.y = ".sample.id")
  dists.m2 <- dists.m2[-which(dists.m2$sample_ID.x == dists.m2$variable), ]
  dists.m2$dup_family <- dists.m2$family.x==dists.m2$family.y

  return(t.test(dists.m2$value[dists.m2$dup_family], dists.m2$value[!dists.m2$dup_family]))
}

euclid.dist(pca_mon)
euclid.dist(pca_nap)

## double check the dist calculations are correct
s1 <- pca[1,c('PC1','PC2')]
s2 <- pca[2,c('PC1','PC2')]
sqrt((s1$PC1 - s2$PC1)^2 + (s1$PC2 -s2$PC2)^2)

all_monnap <- filter_snps(all_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")
pca <- plot_clusters(all_monnap, "site", simplify_output = TRUE,)$pca$data