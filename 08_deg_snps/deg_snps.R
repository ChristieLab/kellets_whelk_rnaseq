#### deg snps
# install.packages("remotes")

vignette("snpR_introduction")
library("snpR")
library("ape")
library("data.table")

setwd("~/KW/08_deg_snps")

#### read in data 
sample_meta <- read.csv("sampleinfo.csv")
my.dat <- read_vcf("rnaspades_annotated_degs.vcf",
                           sample.meta=sample_meta)
head(sample.meta(my.dat))
head(snp.meta(my.dat))

#### filter SNPs 
dat <- filter_snps(my.dat, maf = 0.05, min_loci=0.8)
format_snps(dat, output="vcf", chr="CHROM", outfile="filt_maf05_minloci80_rnaspades_annotated_degs.vcf")

??filter_snps



# export filtered SNP vcf 
tree <- calc_tree(dat)
class(tree$.base$.base)
plot.phylo(tree$.base$.base)


#### get pairwise FST 

dat <- calc_pairwise_fst(dat, facets="site")
get.snpR.stats(dat, "site", "fst")
?calc_pairwise_fst
#### plot PCA
my.colors <- c("goldenrod1","darkred")
plot_clusters(dat, facets="site", alt.palette = my.colors)
?plot_clusters

#### Structure plot using admixture and all SNPs
plot_structure(dat, method="admixture", k=2, admixture_path = "./admixture_macosx-1.3.0/admixture", facet = "site", iterations = 100000, burnin = 20000)

#### Structure plot using STRUCTURE #### 
#### subset SNPs to run structure 
subsamp <- sample(nrow(dat), 5000, FALSE)
filt.dat <- filter_snps(dat[subsamp])
genos <- genotypes(filt.dat)
genos[genos == "."] <- "NN"
filt.dat <- import.snpR.data(genos, snp.meta(filt.dat), sample.meta(filt.dat))

### plot structure 
## plot structure for k= 1-4, try 10 runs per K and use clumpp 
p   <- plot_structure(filt.dat, method="structure", k=2:4, structure_path = "./STRUCTURE/structure",facet = "site", iteration = 100000, burnin = 10000)

pk2 <- readRDS("from_Will/k2_plot.RDS")
p_k1_4 <- readRDS("from_Will//k1_4_plot.RDS")

ggplot(p_k1_4$plot_data, (aes(x=p_k1_4$plot_data$ID, y=p_k1_4$plot_data$Percentage))) + 
         geom_bar(aes(x=p_k1_4$plot_data$ID,  y=p_k1_4$plot_data$Percentage, fill=p_k1_4$plot_data$Cluster))

## plot structure for k=2 
my.colors <- c("goldenrod1","darkred")
p_structure_k2 <- plot_structure(filt.dat, method="structure", k=2, structure_path = "./STRUCTURE/structure",facet = "site", iteration = 100, burnin = 100, alt.palette = my.colors, cleanup = TRUE)

## Evanno Method 
## Code modified from Will Hemstrom's script 

deg_structure <- read.csv(file="structrure_ln.txt", sep=" ")
deg_structure <- as.data.table(deg_structure)

evanno <- deg_structure[, mean(est_ln_prob), by = K]
colnames(evanno)[2] <- "mean_est_ln_prob"
evanno$lnpK <- NA
evanno$lnppK <- NA
evanno$deltaK <- NA
evanno <- evanno[order(evanno$K),]
evanno$sd_est_ln_prob <- deg_structure[, sqrt(var(est_ln_prob)), by = K][[2]]
evanno$lnpK[-1] <- evanno$mean_est_ln_prob[-1] - evanno$mean_est_ln_prob[-nrow(evanno)]
evanno$lnppK[-nrow(evanno)] <- abs(evanno$lnpK[-nrow(evanno)] - evanno$lnpK[-1])
evanno$deltaK[-c(1, nrow(evanno))] <- abs(evanno$lnppK)[-c(1, nrow(evanno))]/evanno$sd_est_ln_prob[-c(1, nrow(evanno))] # no reason to resolve for ln''(K)
infs <- which(is.infinite(evanno$deltaK))
if(length(infs) > 0){
  evanno$deltaK[infs] <- NA
}
evanno_m <- data.table::melt(evanno, id.vars = c("K"))
evanno_m <- evanno_m[which(evanno_m$variable %in% c("deltaK", "mean_est_ln_prob")),]
evanno_m[,variable := as.character(variable)]
evanno_m$value[evanno_m$variable == "deltaK"] <- log10(evanno_m$value[evanno_m$variable == "deltaK"])
evanno_m$variable[evanno_m$variable == "deltaK"] <- "log[10](Delta*K)"
evanno_m$variable[evanno_m$variable == "mean_est_ln_prob"] <- "bar(ln(Prob))"


