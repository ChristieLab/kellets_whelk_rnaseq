# script written by andy lee on 02/15/2023; contact at andymuanlee@gmail.com
# Rversion:
# Neutral SNPs 
library("snpR")
library("ggplot2")
library("ggtree")
library("fastreeR")

setwd("~/KW/09_angsd/")

sample_meta <- read.csv("~/KW/09_angsd//neutral/sample_info.csv")

deg_vcf <- import.snpR.data("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/05_angsd/degs/scripts/kw_deg.vcf", sample.meta = sample_meta)
all_vcf <- import.snpR.data("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/05_angsd/neutral_all/kw_neutral.vcf", sample.meta = sample_meta)

# nodeg_vcf <-import.snpR.data("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes/Bell Scratch/kellets_whelks_rnaseq/a05_angsd/neutral_nodegs/scripts/kw_nodeg.vcf", sample.meta = sample_meta)

#### filter SNPs 
deg_dat <- filter_snps(deg_vcf, maf = 0.05)
all_dat <- filter_snps(all_vcf, maf = 0.05)



# nodeg_dat <- filter_snps(neutral_nodeg_vcf, maf = 0.05, min_loci=0.8)

#### export filtered SNP vcf 
#format_snps(deg_dat, output="vcf", chr="CHROM", outfile="kw_deg_filt_maf05.vcf")
#format_snps(all_dat, output="vcf", chr="CHROM", outfile="kw_all_filt_maf05.vcf")

## calculate statistics 
deg_dat <- calc_pairwise_fst(deg_dat, facets="site")
deg_dat <- calc_het_hom_ratio(deg_dat, facets="site")
deg_dat <- calc_fis(deg_dat)
deg_dat <- calc_ho(deg_dat)
deg_dat <- calc_he(deg_dat)

all_dat <- calc_pairwise_fst(all_dat, facets="site")
all_dat <- calc_het_hom_ratio(all_dat, facets="site")
all_dat <- calc_fis(all_dat)
all_dat <- calc_ho(all_dat)
all_dat <- calc_he(all_dat)

# nodeg_dat <- calc_pairwise_fst(nodeg_dat, facets="site")
# nodeg_dat <- calc_het_hom_ratio(nodeg_dat, facets="site")
# nodeg_dat <- calc_fis(nodeg_dat)
# nodeg_dat <- calc_ho(nodeg_dat)
# nodeg_dat <- calc_he(nodeg_dat)

##### view statistics 

## pairwise FST values 
get.snpR.stats(deg_dat, "site", "fst")
get.snpR.stats(all_dat, "site", "fst")
# get.snpR.stats(nodeg_dat, "site", "fst")


get.snpR.stats(deg_dat, "site", stats="het_hom_ratio",)
get.snpR.stats(deg_dat, stats="fis")
get.snpR.stats(deg_dat, stats="ho")
get.snpR.stats(deg_dat, stats="he")


get.snpR.stats(all_dat, "site", stats="het_hom_ratio",)
get.snpR.stats(all_dat, stats="fis")
get.snpR.stats(all_dat, stats="ho")
get.snpR.stats(all_dat, stats="he")

# get.snpR.stats(nodeg_dat, "site", stats="het_hom_ratio",)
# get.snpR.stats(nodeg_dat, stats="fis")
# get.snpR.stats(nodeg_dat, stats="ho")
# get.snpR.stats(nodeg_dat, stats="he")

############################################### Format different populations ##################################

### without maternal cross 
deg_nocross <- filter_snps(deg_dat[site=-"NAPxMON"], maf = 0.05)
all_nocross <- filter_snps(all_dat[site=-"NAPxMON"], maf = 0.05)
# nodeg_nocross <- nodeg_dat[site=-"NAPxMON"]

### MON and NAP only 
deg_monnap <- filter_snps(deg_dat[site=-c("NAPxMON","POL")], maf = 0.05)
all_monnap <- filter_snps(all_dat[site=-c("NAPxMON","POL")], maf = 0.05)

################################################ plot PCA ################################################################
manual_col4 <- c("#FEB77EFF", "#972C80FF","#d1426fff","#000004FF")
manual_col3 <- c("#FEB77EFF", "#972C80FF","#000004FF")
manual_col2 <- c("#FEB77EFF", "#972C80FF")


plot_clusters(deg_dat, facets="site", alt.palette = manual_col4)   # GO into SI 

deg_monnap <- deg_dat[site=-c("NAPxMON", "POL")]
dup_deg_monnap <- deg_dat[family=c("T1E1", "T1E18", "T1E2", "T1E4", "T1E5", "T1E6", "T6E18", "T6E21", "T6E23", "T6E24", "T6E44", "T6E45")] 

p2 <- plot_clusters(dup_deg_monnap, facets="family")
p <- plot_clusters(deg_monnap, facets="site", alt.palette = manual_col2)
p$plots$pca + 
  ggtitle("A. PCA of DEG SNPs ")

df <- as.data.frame(p$data)
df.dup <- subset(df, pca.family %in% c("T1E1", "T1E18", "T1E2", "T1E4", "T1E5", "T1E6", "T6E18", "T6E21", "T6E23", "T6E24", "T6E44", "T6E45"))

### connect the duplicates 
dupplot <- p$plots$pca + 
  geom_line(data=  df.dup, aes(pca.PC1, pca.PC2, group=pca.family), color="darkgray") + 
  geom_point(data = df.dup, aes(pca.PC1, pca.PC2, color=pca.site)) +
  theme_test()


### export x,y coordinates for mark to calculate distances 
# xy.coord <- dupplot$data[ , c(1,3,4,9,10)]
# xy.coord$duplicate <- xy.coord$family %in% c("T1E1", "T1E18", "T1E2", "T1E4", "T1E5", "T1E6", "T6E18", # "T6E21", "T6E23", "T6E24", "T6E44", "T6E45")
# write.csv(xy.coord, "kw_SNP_xycoord.csv")

###### Random Forest to see how many SNPs exlpain the vairance among samples ###############
rf_deg <- run_random_forest(deg_monnap, response = "site", importance = "permutation")
rf_all <- run_random_forest(all_monnap, response = "site", importance = "permutation")

par(mar = c(1, 1, 1, 1))
plot(density(rf_deg$models$.base_.base$model$variable.importance))
plot(density(rf_all$models$.base_.base$model$variable.importance))

rf_deg$models$.base_.base$predictions
rf_all$models$.base_.base$predictions

deg_monnap_meta <- snp.meta(deg_monnap)
deg_monnap_snp_pos <- paste0(deg_monnap_meta$CHROM, "_", deg_monnap_meta$position)
all_monnap_meta <- snp.meta(all_monnap)
all_monnap_snp_pos <- paste0(all_monnap_meta$CHROM, "_", all_monnap_meta$position)


deg_snps_in_transcriptome <- all_monnap_snp_pos %in% deg_monnap_snp_pos 

# double check that %in% works as expected 
deg_monnap_snp_pos[101]
grep(deg_monnap_snp_pos[101], all_monnap_snp_pos)
deg_snps_in_transcriptome[1741]

rfdf <- as.data.frame(cbind(rf_all$models$.base_.base$model$variable.importance,deg_snps_in_transcriptome))
plot(rfdf$V1, col= factor(rfdf$deg_snps_in_transcriptome), ylab="Importance in random forest", main="Importance of each SNP in random forest model (MON and NAP only)")
 
################################################### Phylogenetic Tree ####################################################
###################################################### fastreeR ##########################################################
##### Note: cannot get bootstrap values


# https://4va.github.io/biodatasci/r-ggtree.html ggtree tutorial 
# create tree 

format_snps(deg_nocross, output="vcf", chr="CHROM", outfile="kw_deg_filt_maf05_nocross.vcf" )####
format_snps(all_nocross, output="vcf", chr="CHROM", outfile="kw_all_filt_maf05_nocross.vcf" )####

### Sample stastistic of the VCF file
### myVcfIstats <- fastreeR::vcf2istats(inputFile = "/Users/andy/KW/09_angsd/kw_all_filt_maf05_nocross.vcf")
### plot(myVcfIstats[,7:9])
### 
### ### Calculate distance from vcf  
### myVcfDist <- fastreeR::vcf2dist(inputFile = "/Users/andy/KW/09_angsd/kw_all_filt_maf05_nocross.vcf", threads = 2)
###       
### ### histogram of distance
### graphics::hist(myVcfDist, breaks = 100, main=NULL, 
###                     xlab = "Distance", xlim = c(0,max(myVcfDist)))
###       
### myVcfTree <- fastreeR::dist2tree(inputDist = myVcfDist)

allsnps_tree <- fastreeR::vcf2tree(inputFile = "/Users/andy/KW/09_angsd/kw_all_filt_maf05_nocross.vcf", threads = 2)
allsnps_tree <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/04_deg_pipeline/data/", "", allsnps_tree) 
# allsnps_tree <- gsub("*.sorted.bam", "", simple_tree)

#use ggtree to customize tree plot
allsnptree <- ape::read.tree(text= allsnps_tree)

p <- ggtree(allsnptree, layout="daylight") %<+% sample_meta + 
 geom_tiplab(aes(label = factor(site), color = site, geom = "label")) + scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF")) +
 theme(legend.position = "none") +
 ggtitle("fastreeR using all SNPs") 

######## plot deg tree #########

degsnps_tree <-  fastreeR::vcf2tree(inputFile = "/Users/andy/KW/09_angsd/kw_deg_filt_maf05_nocross.vcf", threads = 2)

#sub sample names 

#deg_simple_tree <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/a05_angsd/degs/data/","",degsnps_tree)
#deg_simple_tree <- gsub("deg", "sorted", deg_simple_tree) 

deg_tree <- ape::read.tree(text= degsnps_tree)


# Get the dna format 

#format_snps(deg_dat, output="vcf", chr="CHROM", outfile="kw_deg_filt_maf05.vcf")

p2 <- ggtree(deg_tree, layout="daylight") %<+% sample_meta + 
  geom_tiplab(aes(label = factor(site), color = site, geom = "label")) + scale_color_manual(values=c(MON = "#FEB77EFF", NAP = "#972C80FF", POL="#000004FF")) +
  theme(legend.position = "none") +
  ggtitle("fastreeR using deg SNPs") 
 # geom_text(aes(label=node), hjust=-.3) + 

### code to try to bootstrap, doesn't work 
#dna <- format_snps(deg_nocross, output="sn")
# dna <- gsub("/scratch/bell/lee3617/kellets_whelks_rnaseq/a05_angsd/degs/data/","",dna)#
# dna <- gsub("deg", "sorted", dna) 
# ape::boot.phylo(deg_tree, dna, FUN = fastreeR::vcf2dist(inputFile = "/Users/andy/KW/09_angsd/kw_deg_filt_maf05_nocross.vcf"))


################################################### STRUCTURE Plots########################################################


## all samples 
# saveRDS(deg_dat, file = "deg.rds") #1596 SNPs
# saveRDS(all_dat, file = "all_dat.rds") #5000 SNPs

## wild pops only (no cross)
deg_nocross <- deg_vcf[site=-"NAPxMON"]
deg_nocross.filt <- filter_snps(deg_nocross, maf = 0.05) #1644

all_nocross <- all_vcf[site=-"NAPxMON"]
all_nocross.filt <- filter_snps(all_nocross, maf = 0.05)
subsamp <- sample(nrow(all_nocross.filt), 5000, FALSE)
all_nocross.filt <- filter_snps(all_nocross.filt[subsamp])

deg_monnap <- deg_vcf[site=c("MON","NAP")] 
deg_monnap.filt <- filter_snps(deg_monnap, maf = 0.05) 

all_monnap <- all_vcf[site=c("MON","NAP")]
all_monnap.filt <- filter_snps(all_monnap, maf = 0.05)
subsamp.monnap <- sample(nrow(all_monnap.filt), 5000, FALSE)
all_monnap.filt <- filter_snps(all_monnap.filt[subsamp.monnap])


saveRDS(deg_nocross.filt, file = "deg_nocross.filt.rds") #1644 SNPs
saveRDS(all_nocross.filt, file = "all_nocross.filt.rds") #5000 SNPs

saveRDS(deg_monnap.filt, file = "deg_monnap.filt.rds") #1582
saveRDS(all_monnap.filt, file = "all_monnap.filt.rds") #5000 SNPs

## plot structure for k= 1-4, try 10 runs per K and use clumpp 
plot_structure(dat, method="structure", k=4:6, structure_path = "./STRUCTURE/structure",facet = "site", iteration = 100000, burnin = 20000, alt.palette = manual_col)


p_structure <- plot_structure(deg_monnap, method = "structure", structure_path = "../08_deg_snps/STRUCTURE/structure", facet = "site", iterations = 100000, burnin = 10000, k = 2, cleanup = TRUE, alt.palette = manual_col2)

p_structure$plot + 
  ggtitle("B. STRUCTURE analysis using DEG SNPs") +
  xlab("Location")

# plot_structure(nodeg_dat,method = "structure", structure_path = "../08_deg_snps/STRUCTURE/structure", facet = "site", iterations = 10000, iteration = 10000, burnin = 100, k = 2:4, cleanup = TRUE)


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





####### test angsd vs. gatk ###### use %in%
monnap_sample_meta <- sample_meta[sample_meta$site==c("MON") | sample_meta$site==c("NAP") , ]
deg_gatk <- import.snpR.data("../5_GATK/rnaspades_annotated_degs.vcf.gz", sample.meta = monnap_sample_meta)

deg_gatk <- filter_snps(deg_gatk, maf = 0.05)

deg_monnap_meta <- snp.meta(deg_monnap)
deg_gatk_meta <- snp.meta(deg_gatk)

deg_monnap_snp_pos <- paste0(deg_monnap_meta$CHROM, "_", deg_monnap_meta$position)
deg_gatk_snp_pos <- paste0(deg_gatk_meta$CHROM, "_", deg_gatk_meta$position)

table(deg_monnap_snp_pos %in% deg_gatk_snp_pos) #661 shared SNPs, 835 identified by ANGSD but not in gatk 
table(deg_gatk_snp_pos%in% deg_monnap_snp_pos) #661 shared SNPs, 40361 identified by gatk but not in angsd 

# can merge SNPs from gatk AND angsd if desired 
m <- merge_snpRdata(deg_monnap, deg_gatk)
?merge_snpRdata
