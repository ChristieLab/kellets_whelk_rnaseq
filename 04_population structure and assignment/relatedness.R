# script written by andy lee on 06/19/2024; contact at andymuanlee@gmail.com
# Rversion:

library(reshape2)
library(snpR)
library(ggplot2)
setwd("~/KW/09_angsd/relatedness/relatedness_individuals/")

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
  dists.m2 <- merge(dists.m, pca, by.x = "variable", by.y = ".sample.id")
  dists.m2 <- dists.m2[-which(dists.m2$sample_ID.x == dists.m2$variable), ]
  dists.m2 <- dists.m2[,1:11]
  dists.m2$dup_family <- dists.m2$family.x==dists.m2$family.y
  return(dists.m2)
  # return(t.test(dists.m2$value[dists.m2$dup_family], dists.m2$value[!dists.m2$dup_family]))
}

mon_dist <- euclid.dist(pca_mon)
nap_dist <- euclid.dist(pca_nap)

dist_dat <- mon_dist

## plot distributions 
plot_dist_distribution <- function(dist_dat){
  # split up related and unrelated samples
  unrelated <- dist_dat[!dist_dat$dup_family, ]
  related <- dist_dat[dist_dat$dup_family,]
  
  # make the plot
  p <- ggplot(data = unrelated) +
    geom_histogram(aes(x=value), binwidth = 0.01) + # plot euclid dist distribution
    geom_vline(data=unrelated, aes(xintercept=mean(value)), size=0.3, linetype="dotdash") + # add in unrelated samples 
    
    geom_histogram(data= related, aes(x=value), binwidth = 0.01, color="darkorange") + 
    geom_vline(data=related, aes(xintercept=mean(value)), size=0.3, linetype="dotdash", color="darkorange") +
    
    coord_cartesian(xlim=c(0,1)) +
    xlab(label = "Pairwise Euclidean Distance") 
  
  return(p)
}

cowplot::plot_grid(
  plot_dist_distribution(mon_dist) + ggtitle("Monterey"),
  plot_dist_distribution(nap_dist) + ggtitle("Naples"), 
  ncol=1
)

