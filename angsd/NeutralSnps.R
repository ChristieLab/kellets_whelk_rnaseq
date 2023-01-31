# Neutral SNPs 
library("snpR")
library("ape")
library("data.table")

setwd("~/KW/09_neutral_snps")

#### read in data 
sample_meta <- read.csv("sample_info.csv")

snp.dat <- read.table("test.geno")
genos <- snp.dat[,-c(1:2)]
snp_meta <- snp.dat[,1:2]

full.dat <- import.snpR.data(genos ,snp.meta = snp_meta, sample.meta = sample_meta)
head(sample.meta(full.dat))
head(snp.meta(full.dat))

#### filter SNPs 
dat <- filter_snps(full.dat, maf = 0.05, min_loci=0.8)
dat2 <- filter_snps(full.dat, min_loci=0.8)
#### # export filtered SNP vcf 
format_snps(dat, output="vcf", chr="CHROM", outfile="filt_maf05_minloci80_neutral.vcf")

### get pairwise FST 
dat <- calc_pairwise_fst(dat, facets="site")
get.snpR.stats(dat, "site", "fst")

# without MAF filter, FST value is much lower 
dat2 <- calc_pairwise_fst(dat2, facets="site")
get.snpR.stats(dat2, "site", "fst") 


#### plot PCA
my.colors <- c("goldenrod1","darkred")
plot_clusters(dat, facets="site")
plot_clusters(dat2, facets="site")

?plot_clusters

#### Structure plot using admixture and all SNPs
