#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###
### analyzing Kellet's whelks RNA-Seq SNPs  ###
### written by Andy Lee April 2023, contact at andymuanlee@gmail.com, all rights reserved 
### 
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++###

# load packages and set working diretory 
library("snpR")
library("ggplot2")
library("ggtree")
library("fastreeR")
library(tidyverse)
library(gridExtra)
library(dplyr)
library(adegenet)
library(pegas)
library("poppr")
library("ape") 

setwd("~/KW/09_angsd/")

### load raw vcf files 
sample_meta <- read.csv("~/KW/09_angsd/neutral/sample_info.csv")
deg_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_degsnps.vcf", sample.meta = sample_meta)
all_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_allsnps.vcf", sample.meta = sample_meta)

### filter by maf within each site 
all_filt <- filter_snps(all_vcf, maf = 0.05, maf_facets = "site") 
snp.meta(all_filt)$ID <- paste0(snp.meta(all_filt)$CHROM, "_",snp.meta(all_filt)$position)
deg_filt <- filter_snps(deg_vcf, maf = 0.05, maf_facets = "site")
snp.meta(deg_filt)$ID <- paste0(snp.meta(deg_filt)$CHROM, "_",snp.meta(deg_filt)$position)


plot_clusters(deg_filt, facets="site")

####format different populations============
## without maternal cross 
deg_nocross <- filter_snps(deg_vcf[site=-"NAPxMON"], maf = 0.05, maf_facets = "site")
snp.meta(deg_nocross)$ID <- paste0(snp.meta(deg_nocross)$CHROM, "_",snp.meta(deg_nocross)$position)
all_nocross <- filter_snps(all_vcf[site=-"NAPxMON"], maf = 0.05, maf_facets = "site")
snp.meta(all_nocross)$ID <- paste0(snp.meta(all_nocross)$CHROM, "_",snp.meta(all_nocross)$position)

## MON and NAP only 
deg_monnap <- filter_snps(deg_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")
all_monnap <- filter_snps(all_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")

# Filter out "high" FST snps
all_nocross <- calc_pairwise_fst(all_nocross, facets="site")
fst <- get.snpR.stats(all_nocross, "site", "fst")
fsts <- fst$pairwise

keep <- unique(fsts$ID[which(fsts$fst <= 0.01)])
keep <- which(snp.meta(all_nocross)$ID %in% keep)
lowfst <- all_nocross[keep,]

nsnps(lowfst)/nsnps(all_nocross)


# Filter out high importance snps in random forest 
#rfdf <- readRDS("rfdf.rds")
#rfdf[(nrow(sorted.rfdf) * 0.9),]$Importance
## keep <- unique(rfdf$snp_pos[which(rfdf$Importance <= rfdf[(nrow(sorted.rfdf) * 0.9#),]$Importance)])
#keep <- unique(rfdf$snp_pos[which(rfdf$Importance <= 0)]) 
#keep <- which(snp.meta(all_nocross)$ID %in% keep)
#low_import <- all_nocross[keep,] # 60941 snps 
#plot_clusters(low_import, facets="site", alt.palette=manual_col3)
#nsnps(low_import)/nsnps(all_nocross)
#colnames(rfdf)


### PCAs ============================================================
# set manual colors 
manual_col4 <- c("#000004FF", "#D1426FFF", "#972C80FF", "#FEB77EFF")
manual_col3 <- c("#000004FF", "#D1426FFF", "#FEB77EFF")
# manual_col2 <- c("#FEB77EFF", "#972C80FF")





deg_pca <- plot_clusters(deg_nocross, facets="site", alt.palette = manual_col3)
all_pca <- plot_clusters(all_nocross, facets="site", alt.palette = manual_col3)
low_fst_pca <- plot_clusters(lowfst, facets="site", alt.palette=manual_col3)

plot_clusters(lowfst, facets="site", alt.palette=manual_col3) 


p_pca_deg <- deg_pca$plots$pca + theme(legend.position = "none")
p_pca_all <- all_pca$plots$pca + theme(legend.position = "none")
p_pca_low <- low_fst_pca$plots$pca + theme(legend.position = "none")

## microsat PCA ====
## data preparation 
usat <- read.csv("kw_usat.csv")
num_rows <- nrow(usat) / 2  # data contains an even number of rows
num_cols <- ncol(usat) - 1  # Subtract 1 to exclude the first column (individual names)
microsat_data <- matrix(NA, nrow = num_rows, ncol = num_cols)

i=2
row_index <- 1
for (i in seq(2, nrow(usat), by = 2)){
  combined_row <- paste0(usat[i-1, -1], "/", usat[i, -1])  # Combine alleles with "/"
  if (length(combined_row) == num_cols) {
    microsat_data[row_index, ] <- combined_row
    row_index <- row_index + 1
  }
}
names <- usat[,1]
names <- names[nzchar(names)]
colnames(microsat_data) <- colnames(usat)[-1]
rownames(microsat_data) <- names

genind <- df2genind(as.data.frame(microsat_data), ncode=9, sep="/") # create s4 object 
pops <- c(rep("MON", 78), rep("DIC",32), rep("POL", 64), rep("NAP",47)) 
genind$pop <- as.factor(pops)

genind_main <- genind[pop=c("MON","NAP","POL")] #remove DIC from analysis 
x.gen <- tab(genind_main)
pca.x <- dudi.pca(x.gen, scannf = FALSE, nf = 2)
pca.x$eig
PC1var <- round(pca.x$eig[1]/sum(pca.x$eig)*100, 2)
PC2var <- round(pca.x$eig[2]/sum(pca.x$eig)*100, 2)


pca.x$li$pop <- pop(genind_main)
p_pca_usat <- ggplot(pca.x$li)  +
  geom_point(aes(x=Axis1, y=Axis2, col=pop)) +
  scale_color_manual(values=c("#FEB77EFF", "#972C80FF","#000004FF")) + 
  labs(x=paste0("PC1 (", PC1var,"%)"),
       y=paste0("PC2 (", PC2var,"%)")) +
  # ggtitle("PCA of 189 KW individuals at 9 microsatellite loci") +
  theme_bw() + 
  theme(legend.position = "none") 


###STRUCTURE plots ===============================================
## wild pops only (no cross)
#deg_nocross 
#all_nocross 
#subsamp.nocross <- sample(nrow(all_nocross), 5000, FALSE)
#all_nocross_subsamp <- filter_snps(all_nocross[subsamp.nocross])

#saveRDS(deg_nocross, file = "deg_nocross.filt.rds") #1644 SNPs
#saveRDS(all_nocross_subsamp, file = "all_nocross.filt.rds") #5000 SNPs

#subsamp.lowfst <- sample(nrow(lowfst), 5000, FALSE)
#saveRDS(lowfst, file = "low_fst_nocross.filt.rds")

## plot structure for k= 1-4, try 10 runs per K and use clumpp, run by Will Hemstrom 
structure_plots <- readRDS(file = "all_and_degs_nocross_structure_plots.RDS")

p_str_usat <- structure_plots$microsats +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.text.y = element_blank()
  ) 

p_str_all <- structure_plots$all +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text.y = element_blank())

p_str_deg <- structure_plots$degs +
  theme(legend.position = "none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
) 


### Phylogenetic Tree=========
## fastreeR 
## Note: cannot get bootstrap values
### Sample stastistic of the VCF file
# myVcfIstats <- fastreeR::vcf2istats(inputFile = "/Users/andy/KW/09_angsd/kw_all_filt_maf05_nocross.vcf")
# plot(myVcfIstats[,7:9])
# 
### Calculate distance from vcf  
# myVcfDist <- fastreeR::vcf2dist(inputFile = "/Users/andy/KW/09_angsd/kw_all_filt_maf05_nocross.vcf", threads = 2)
#       
### histogram of distance
# graphics::hist(myVcfDist, breaks = 100, main=NULL, xlab = "Distance", xlim = c(0,max(myVcfDist)))

### myVcfTree <- fastreeR::dist2tree(inputDist = myVcfDist)

## fastreeR uses vcf files as input 
format_snps(all_nocross, output="vcf", chr="CHROM", outfile="all_nocross_filt.vcf")
format_snps(deg_nocross, output="vcf", chr="CHROM", outfile="deg_nocross_filt.vcf")
format_snps(lowfst , output="vcf", chr="CHROM", outfile="low_fst_nocross.vcf")

### get trees directly
allsnps_tree <- fastreeR::vcf2tree(inputFile = "/Users/andy/KW/09_angsd/all_nocross_filt.vcf", threads = 2)
allsnps_tree <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/04_deg_pipeline/data/", "", allsnps_tree) 
#allsnps_tree <- gsub("*.sorted.bam", "", simple_tree)


allsnptree <- ape::read.tree(text= allsnps_tree) #use ggtree to customize tree plot
all_tree <- ggtree(allsnptree, layout="daylight") %<+% sample_meta + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label"), size=1.5) + scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF")) +
  theme(legend.position = "none") 
  # ggtitle("fastreeR using all SNPs") 

degsnps_tree <- fastreeR::vcf2tree(inputFile = "/Users/andy/KW/09_angsd/deg_nocross_filt.vcf", threads = 2)
degsnps_tree <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/a05_angsd/degs/data/", "", degsnps_tree) 
deg_sample_meta <- dplyr::filter(sample_meta, site!="NAPxMON") #dplyr
deg_sample_meta$sample_ID <-  gsub("sorted","deg", deg_sample_meta$sample_ID)
deg_sample_meta <- as.data.frame(deg_sample_meta)

degsnps_tree <- ape::read.tree(text= degsnps_tree)
deg_tree <- ggtree(degsnps_tree, layout="daylight") %<+% deg_sample_meta + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label"), size=1.5) +
  scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF")) +
  theme(legend.position = "none")

# low importance tree
low_fst_nocross <- fastreeR::vcf2tree(inputFile = "/Users/andy/KW/09_angsd/low_fst_nocross.vcf", threads = 2)
low_fst_tree <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/04_deg_pipeline/data/", "", low_fst_nocross) 

#use ggtree to customize tree plot
low_fst_tree <- ape::read.tree(text= low_fst_tree)
p_low_fst_tree <- ggtree(low_fst_tree, layout="daylight") %<+% sample_meta + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label"), size=1.5) + scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF")) +
  theme(legend.position = "none") 

# ggtitle("fastreeR using all SNPs") 



### microsat NJ tree
## Plot NJ tree ====
set.seed(42)
genind_main # genind obj with all individuals from MON, NAP and POL 
sub_genind <- genind_main[sample(nInd(genind_main), 50)] #subset for a visually simpler tree
# d <- dist(genind_main) #full tree
d <- dist(sub_genind)
njtree <- d %>%
  nj() %>%    # calculate neighbor-joining tree
  ladderize() # organize branches by clade

nInd(genind_main)
pop_info <- as.data.frame(cbind(rownames(sub_genind$tab), as.character(pop(sub_genind))))
# pop_info <- as.data.frame(cbind(rownames(genind_main$tab), as.character(pop(genind_main)))) # full tree
colnames(pop_info) <- c("sample","site")

usat_sub_tree <- ggtree(njtree, layout="daylight") %<+% pop_info + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label"), size=1.5) +   
  scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF")) + 
  theme(legend.position = "none")
  # ggtitle("NJ tree using 9 microsatellite loci and 189 KW inviduals from 3 populations") 


### plot figure 2 ====
lay <- cbind(c(1,2,3),
             c(4,5,6),
             c(7,8,9),
             c(10,11,12))

pdf("~/KW/figures/figure 2/fig2.pdf", width = 11, height = 8.5)
grid.arrange(p_pca_usat, p_str_usat, usat_sub_tree,
             p_pca_low, p_str_low_fst, p_low_fst_tree,
             p_pca_all,  p_str_all,  all_tree, 
             p_pca_deg,  p_str_deg,  deg_tree, 
             layout_matrix= lay)
dev.off()

#pdf("~/KW/figures/supplemental/usat_full_tree.pdf", width = 8.5, height = 11)
#usat_tree
#dev.off()


# 1 snp per contig -- supplemental ====
set.seed(42)
keep <- snp.meta(all_nocross) %>% group_by(CHROM) %>% sample_n(1)
keep <- which(snp.meta(all_nocross)$ID %in% keep$ID)

keep_deg <- snp.meta(deg_nocross) %>% group_by(CHROM) %>% sample_n(1)
keep_deg <- which(snp.meta(deg_nocross)$ID %in% keep_deg$ID)

all_no_ld <- all_nocross[keep,]
deg_no_ld <- deg_nocross[keep_deg,]

# PCAs
p_all_no_ld_pca <- plot_clusters(all_no_ld, facets="site", alt.palette = manual_col3)
p_deg_no_ld_pca <- plot_clusters(deg_no_ld, facets="site", alt.palette = manual_col3)

# NJ tree

# transcriptome wide tree, 1 snp per contig 
format_snps(all_no_ld, output="vcf", chr="CHROM", outfile="all_no_ld.vcf")
all_no_ld_tree <- fastreeR::vcf2tree(inputFile = "/Users/andy/KW/09_angsd/all_no_ld.vcf", threads = 2)
all_no_ld_tree <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/04_deg_pipeline/data/", "", all_no_ld_tree) 
all_no_ld_tree <- ape::read.tree(text= all_no_ld_tree)
p_all_no_ld_tree <- ggtree(all_no_ld_tree, layout="daylight") %<+% sample_meta + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label"), size=1.5) + scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF")) +
  theme(legend.position = "none") + 
  ggtitle("transcriptome-wide SNPs NJ tree, 1 SNP per contig")

### DEG Tree, 1 snp per contig 
format_snps(deg_no_ld, output="vcf", chr="CHROM", outfile="deg_no_ld.vcf")
deg_no_ld_tree <- fastreeR::vcf2tree(inputFile = "/Users/andy/KW/09_angsd/deg_no_ld.vcf", threads = 2)
deg_no_ld_tree <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/a05_angsd/degs/data/", "", deg_no_ld_tree) 
deg_no_ld_tree <- ape::read.tree(text= deg_no_ld_tree)
deg_sample_meta <- dplyr::filter(sample_meta, site!="NAPxMON") #dplyr
deg_sample_meta$sample_ID <-  gsub("sorted","deg", deg_sample_meta$sample_ID)
deg_sample_meta <- as.data.frame(deg_sample_meta)

p_deg_no_ld_tree <- ggtree(deg_no_ld_tree, layout="daylight") %<+% deg_sample_meta + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label"), size=1.5) + scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF")) +
  theme(legend.position = "none") + 
  ggtitle("DEGs only NJ tree, 1 SNP per contig") 


# Structure from Will 
saveRDS(deg_no_ld, file = "deg_no_ld.rds") 
saveRDS(all_no_ld, file = "all_no_ld.rds")

lay <- rbind(c(1,4),
             c(2,5),
             c(3,6))


pdf("~/KW/figures/supplemental /onesnppercontig.pdf", width = 11, height = 8.5)
grid.arrange(p_all_no_ld_pca$plots$pca, p_all_no_ld_tree, 
             p_deg_no_ld_pca$plots$pca, p_deg_no_ld_tree)
dev.off()








####### Figure 3 =====
# FST distribution ==============
### look at FST distribution 
# get pairwise FST values to plot against importance 
deg_monnap <- calc_pairwise_fst(deg_monnap, facets="site")
all_monnap <- calc_pairwise_fst(all_monnap, facets="site")
all_monnap_fst <-  get.snpR.stats(all_monnap, "site", "fst")
deg_monnap_fst <-  get.snpR.stats(deg_monnap, "site", "fst")

all_fst_sorted <- all_monnap_fst$pairwise[order(-all_monnap_fst$pairwise$fst), ]
str(all_fst_sorted)

bins <- seq(-0.1, 0.5, by=0.01)
deg_fst <- deg_monnap_fst$pairwise

p <- ggplot(deg_fst) +
  geom_histogram(aes(x=fst), breaks=bins) + 
  ggtitle("deg snps FST distribution") +
  theme_bw()

all_fst <- all_monnap_fst$pairwise
p1 <- ggplot(all_fst) +
  geom_histogram(aes(x=fst), breaks=bins) + 
  ggtitle("transcriptome-wide snps FST distribution") + 
  theme_bw()

#pdf("~/KW/figures/supplemental /fst_distribution_monnap.pdf", width = 8.5, height = 8.5)
grid.arrange(p, p1, ncol=1)
#dev.off()

# FST vs. LFC plot ==========
## get FST values and LFC for each contig 
all_monnap <- calc_pairwise_fst(all_monnap, facets="site")
deg_monnap <- calc_pairwise_fst(deg_monnap, facets="site")
all_monnap_fst <- get.snpR.stats(all_monnap, facets="site", stats="fst")
deg_monnap_fst <- get.snpR.stats(deg_monnap, facets="site", stats="fst")

x <- all_monnap_fst$pairwise
z <- deg_monnap_fst$pairwise
res <- readRDS(file = "deseq_monnap.rds")
res$CHROM <- rownames(res)
  temp <- do.call(rbind, strsplit(res$CHROM, "\\|"))
  res$CHROM <- temp[,1]
y <- res

df <- merge(x, y, by="CHROM")
View(df)

df$fst
df$log2FoldChange

p <- ggplot(df) +
  geom_point(aes(x=fst, y=log2FoldChange), size = 0.3) + 
  geom_smooth(method=lm, aes(x=fst, y=log2FoldChange)) +
  ggtitle("Transcriptome-wide SNPs")
p

## deg snps only 
deg_df <- merge(z, y, by="CHROM")
p2 <- ggplot(deg_df) +
  geom_point(aes(x=fst, y=log2FoldChange), size = 0.3) + 
  geom_smooth(method=lm, aes(x=fst, y=log2FoldChange)) +
  ggtitle("DEG SNPs")

## test for interaction 
p_lfc_fst <- ggplot() + 
  geom_point(data=df, aes(x=fst, y=log2FoldChange), size = 0.3) + 
  geom_smooth(data=df, method=lm, aes(x=fst, y=log2FoldChange), col="black") +
  geom_point(data=deg_df, aes(x=fst, y=log2FoldChange), size = 0.3, col="darkred") + 
  geom_smooth(data=deg_df, method=lm, aes(x=fst, y=log2FoldChange), col="darkred") + 
  theme_bw()
  
# with text
p_lfc_fst_text <- p_lfc_fst + 
  geom_text(aes(x = .29, y = -1.70, label= paste0("Transcriptome-wide SNPs, Slope: ", round(coef(lm(log2FoldChange ~ fst, data = df))["fst"], 2),", p < ", round(summary(lm(log2FoldChange ~ fst, data = df))$coefficients[2, 4], 2)))) + 
  geom_text(aes(x = .34, y = 1.65, label= paste0("DEG SNPs, Slope: ", round(coef(lm(log2FoldChange ~ fst, data = deg_df))["fst"], 2),", p < ", round(summary(lm(log2FoldChange ~ fst, data = deg_df))$coefficients[2, 4], 2))), color="darkred") +
  theme(legend.position = "none")

## only one FST per contig 
set.seed(42)
sampled_df <- df %>%
  group_by(CHROM) %>%
  sample_n(1)

p1 <- ggplot(sampled_df, aes(x=fst, y=log2FoldChange)) +
  geom_point(size=0.3) + 
  geom_smooth(method=lm, aes(x=fst, y=log2FoldChange)) +
  ggtitle("Transcriptome-wide SNPs, only one SNP per contig")

#pdf("~/KW/figures/supplemental /fstvslfc.pdf", width = 11, height = 8.5)
#grid.arrange(p, p1, p2, ncol=3)
#dev.off()

##### Supplemental ==== 
# SNP Verification =====

### make sure that all contigs are within the DEG transcriptome 
deg_angsd <- filter_snps(deg_monnap, maf = 0.05, maf_facets = "site")
deg_angsd_meta <- snp.meta(deg_angsd)
deg_angsd_snp_CHROM <- as.data.frame(deg_angsd_meta$CHROM)

length(unique(deg_angsd_meta$CHROM)) # number of contigs in the deg set with SNPs
deg_contig_names <- read.csv("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/05_angsd/degs/rnaspades_degs_contignames.txt", header = FALSE)
table(deg_angsd_snp_CHROM$`deg_angsd_meta$CHROM` %in% deg_contig_names$V1) # should all return true 

## Figure out how many snps are on each contig 
hist(table(deg_angsd_snp_CHROM), breaks=30)

sort(table(deg_angsd_snp_CHROM))
# highest SNPs per contig is on NODE_16533_length_6329_cov_5591.178389_g366_i14
# NODE_87022_length_3018_cov_8698.421732_g5957_i36
# NODE_12782_length_6875_cov_7230.180388_g3001_i3


##  higher than expected SNPs per contigs, checking to see if the SNPs on the same contigs are high in LD 
snp.meta(deg_monnap)$CHROM <- gsub("\\.", "__", snp.meta(deg_monnap)$CHROM)
deg_monnap <- calc_pairwise_ld(deg_monnap, facets = "site")
deg_monnap_ld <- get.snpR.stats(deg_monnap, "site", "LD")
View(deg_monnap_ld$prox)

deg_monnap <- readRDS("~/../Downloads/deg_monnap_LD.RDS")
snp.meta(deg_monnap)$CHROM <- gsub("\\.", "__", snp.meta(deg_monnap)$CHROM)
deg_monnap <- calc_pairwise_ld(deg_monnap, facets = "site.CHROM", verbose = TRUE, CLD = FALSE)
deg_monnap_ld <- get.snpR.stats(deg_monnap, "site.CHROM", "LD")
snp.tab <- snp.meta(deg_monnap)
snp.tab <- table(snp.tab$CHROM)
deg_monnap_ld$prox[,meanDprime := mean(Dprime, na.rm = TRUE), by = .(s1_CHROM, sample.subfacet)]
r <- deg_monnap_ld$prox
r <- unique(r[,c("s1_CHROM", "sample.subfacet", "meanDprime")])
r$nsnps <- snp.tab[match(r$s1_CHROM, names(snp.tab))]
library(ggplot2)
ggplot(r, aes(x = nsnps, y = meanDprime, color = sample.subfacet)) +
  geom_point() +
  theme_bw() +
  khroma::scale_color_highcontrast()

# Examine duplicates ====
dup_deg_monnap <- deg_monnap[family=c("T1E1", "T1E18", "T1E2", "T1E4", "T1E5", "T1E6", "T6E18", "T6E21", "T6E23", "T6E24", "T6E44", "T6E45")] 

p2 <- plot_clusters(dup_deg_monnap, facets="family")
p <- plot_clusters(deg_monnap, , alt.palette = manual_col2)
df <- as.data.frame(p$data)
df.dup <- subset(df, pca.family %in% c("T1E1", "T1E18", "T1E2", "T1E4", "T1E5", "T1E6", "T6E18", "T6E21", "T6E23", "T6E24", "T6E44", "T6E45"))

dupplot <- p$plots$pca + 
  geom_line(data=  df.dup, aes(pca.PC1, pca.PC2, group=pca.family), color="darkgray") + 
  geom_point(data = df.dup, aes(pca.PC1, pca.PC2, color=pca.site)) +
  theme_test()


### export x,y coordinates for mark to calculate distances 
# xy.coord <- dupplot$data[ , c(1,3,4,9,10)]
# xy.coord$duplicate <- xy.coord$family %in% c("T1E1", "T1E18", "T1E2", "T1E4", "T1E5", "T1E6", "T6E18", # "T6E21", "T6E23", "T6E24", "T6E44", "T6E45")
# write.csv(xy.coord, "kw_SNP_xycoord.csv")

#### p-values from Mark 
rnaseq_mon_pval <- read.table("./duplicate_comparison/p-vals_rnaseq_Monterey.txt")
rnaseq_nap_pval <- read.table("./duplicate_comparison/p-vals_rnaseq_Naples.txt")
snps_mon_pval   <- read.table("./duplicate_comparison/p-vals_snps_Monterey.txt")
snps_nap_pval   <- read.table("./duplicate_comparison/p-vals_snps_Naples.txt")

# Code Archive ====
### filter SNPs 
# removing snps most out of HWE 
# look within site and do multiple testing correction 
# remove ones most out of HWE -- maybe where pools are heterogeneity 
# tried maf = 0.01 and 0.03, does not change patterns much 

#all_filt <- filter_snps(all_vcf, hwe = .05, hwe_facets ="site", fwe_method = "holm", maf = 0.05, maf_facets = "site") #63031
#deg_filt <- filter_snps(deg_vcf, hwe = .05, hwe_facets ="site", fwe_method = "holm", maf = 0.05, maf_facets = "site")



####### test angsd vs. gatk ###### use %in%

monnap_sample_meta <- sample_meta[sample_meta$site==c("MON") | sample_meta$site==c("NAP") , ]

deg_gatk <- import.snpR.data("../5_GATK/rnaspades_annotated_degs.vcf.gz", sample.meta = monnap_sample_meta)

deg_gatk <- filter_snps(deg_gatk, hwe = .05, hwe_facets ="site", fwe_method = "holm", maf = 0.05, maf_facets = "site")
deg_angsd_meta <- snp.meta(deg_angsd)


deg_angsd_snp_pos <- paste0(deg_angsd_meta$CHROM, "_", deg_angsd_meta$position)
deg_gatk_snp_pos <- paste0(deg_gatk_meta$CHROM, "_", deg_gatk_meta$position)


#### trees stuff 
### export for raxml (fasta or phylip format)
nocross.filt <- filter_snps(all_vcf[site=-"NAPxMON"], hwe = .05, hwe_facets ="site", fwe_method = "holm", maf = 0.05, maf_facets = "site")
subsamp.nocross.filt <- sample(nrow(nocross.filt), 5000, FALSE)
all_nocross.subsamp.filt <- filter_snps(nocross.filt[subsamp.nocross.filt])
format_snps(all_nocross.subsamp.filt, output="fasta", chr="CHROM", outfile="/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/raxml/all_noscross_subsamp_filt.fasta") 

all.monnap.filt <- filter_snps(all_vcf[site=-c("NAPxMON","POL")], hwe = .05, hwe_facets ="site", fwe_method = "holm", maf = 0.05, maf_facets = "site")
subsamp.monnap <- sample(nrow(all.monnap.filt), 5000, FALSE)
all_monnap_subsamp_filt  <- filter_snps(all.monnap.filt[subsamp.monnap])
sample.meta(all_monnap_subsamp_filt)
#format_snps(all_monnap.subsamp, output="fasta", chr="CHROM", outfile="/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/raxml/all_monnap_subsamp.fasta") 

#format_snps(deg_monnap, output="fasta", chr="CHROM", outfile="/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/raxml/deg_monnap_maf05.fasta") 

#format_snps(deg_nocross, output="fasta", chr="CHROM", outfile="/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/raxml/deg_nocross_maf05.fasta") 

raxml_tree <- read.tree("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/raxml/RAxML_bestTree.all_nocross_subsamp.tre") 
# read in the tree output from SNPhylo 
raxml_tree$tip.label <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/04_deg_pipeline/data/","",raxml_tree$tip.label)

ggtree(raxml_tree, layout="daylight") %<+% sample_meta + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label")) + 
  scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF", NAPxMON="#d1426fff")) +
  theme(legend.position = "none") +
  ggtitle("RAXML using 5k neutral SNPs") 


raxml_deg_tree <- read.tree("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/raxml/deg_tree/RAxML_bestTree.deg_nocross.tre") 
raxml_deg_tree$tip.label <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/a05_angsd/degs/data/","",raxml_deg_tree$tip.label)
# deg_sample_meta <- filter(sample_meta, site==c("MON","NAP","POL")) #dplyr
deg_sample_meta <- sample_meta
deg_sample_meta$sample_ID <-  gsub("sorted","deg", deg_sample_meta$sample_ID)
deg_sample_meta <- as.data.frame(deg_sample_meta)

ggtree(raxml_deg_tree, layout="daylight") %<+% deg_sample_meta + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label")) + 
  scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF", NAPxMON="#d1426fff")) +
  theme(legend.position = "none") +
  ggtitle("RAXML using deg SNPs (nocross)") 

# export monnap vcf for snphylo 
snp.meta(all_monnap)$CHROM <- 1
format_snps(all_nocross, output="vcf", chr="CHROM", outfile="/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/monnap_maf05.vcf") 
format_snps(all_monnap, output="vcf", chr="CHROM", outfile="/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/monnap_maf05.vcf") 


# plot snphylo trees 
kw_phyDat <- phangorn::read.phyDat("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/snphylo/kw_deg.fasta", format = "fasta") # read in aligned fasta file, an output from SNPhylo
kw_tree <- read.tree("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/01_tree/snphylo/kw_deg.ml.tree")# read in the tree output from SNPhylo 
sample_meta$sample_ID <- gsub(".sorted.bam","",sample_meta$sample_ID)
ggtree(kw_tree, layout="daylight") %<+% sample_meta + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label")) + 
  scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF", NAPxMON="#d1426fff")) +
  theme(legend.position = "none") +
  ggtitle("snphylo using deg SNPs and all populations") 



## evanno 
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


#### export filtered SNP vcf 
#format_snps(deg_dat, output="vcf", chr="CHROM", outfile="kw_deg_filt_maf05.vcf")
#format_snps(all_dat, output="vcf", chr="CHROM", outfile="kw_all_filt_maf05.vcf")

## explort snp position for more SNP calling of TAE samples 
#all_snp_rf <- paste0(snp.meta(all_filt)$CHROM,":",snp.meta(all_filt)$position)
#write.csv(all_snp_rf, "/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelk_tae/scripts/kw_all_snp_rf.txt", row.names=FALSE, quote=FALSE, col.names = FALSE )
