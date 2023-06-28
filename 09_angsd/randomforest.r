#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++ Random Forest to see how many SNPs explain the vairance among samples +++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(snpR)
library(patchwork)
library(ggplot2)
library(ggblend)

sample_meta <- read.csv("~/KW/09_angsd/neutral/sample_info.csv")
deg_vcf <- import.snpR.data("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/05_angsd/degs/scripts/kw_deg.vcf", sample.meta = sample_meta)
all_vcf <- import.snpR.data("/Users/andy/Library/Group Containers/G69SCX94XU.duck/Library/Application Support/duck/Volumes.noindex/Bell Scratch/kellets_whelks_rnaseq/05_angsd/neutral_all/kw_neutral.vcf", sample.meta = sample_meta)

deg_monnap <- filter_snps(deg_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")
all_monnap <- filter_snps(all_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")

###### Random Forest to see how many SNPs exlpain the vairance among samples ###############
rf_deg <- run_random_forest(deg_monnap, response = "site", importance = "permutation", num.trees =1000000)
rf_all <- run_random_forest(all_monnap, response = "site", importance = "permutation", num.trees = 1000000)

par(mar = c(1, 1, 1, 1))
plot(density(rf_deg$models$.base_.base$model$variable.importance))
plot(density(rf_all$models$.base_.base$model$variable.importance))

# See how well random forest models perfromed 
rf_deg$models$.base_.base$predictions
rf_all$models$.base_.base$predictions

# use metadata to find DEGs 
deg_monnap_meta <- snp.meta(deg_monnap)
deg_monnap_snp_pos <- paste0(deg_monnap_meta$CHROM, "_", deg_monnap_meta$position)
all_monnap_meta <- snp.meta(all_monnap)
all_monnap_snp_pos <- paste0(all_monnap_meta$CHROM, "_", all_monnap_meta$position)
deg_snps_in_transcriptome <- all_monnap_snp_pos %in% deg_monnap_snp_pos 

# get pairwise FST values to plot against importance 
deg_monnap <- calc_pairwise_fst(deg_monnap, facets="site")
all_monnap <- calc_pairwise_fst(all_monnap, facets="site")
deg_monnap_fst <-  get.snpR.stats(deg_monnap, "site", "fst")
all_monnap_fst <-  get.snpR.stats(all_monnap, "site", "fst")
all_monnap_fst$pairwise

bins <- seq(-0.1, 0.5, by=0.01)

deg_fst <- deg_monnap_fst$pairwise
p_degsnps <- ggplot(deg_fst) +
  geom_histogram(aes(x=fst), breaks=bins) + 
  ggtitle("deg snps FST distribution")

all_fst <- all_monnap_fst$pairwise

p_allsnps <- ggplot(all_fst) +
  geom_histogram(aes(x=fst), breaks=bins) + 
  ggtitle("transcriptome-wide snps FST distribution")

grid.arrange(p, p1, ncol=1)


# double checking the 
rf_all$models$.base_.base$model$variable.importance
all_monnap_fst$pairwise$snp_pos <- paste0(all_monnap_fst$pairwise$CHROM, "_",all_monnap_fst$pairwise$position)
identical(all_monnap_fst$pairwise$snp_pos, all_monnap_snp_pos)

# create importance plot 
rfdf <- as.data.frame(cbind(all_monnap_fst$pairwise$snp_pos,
                          all_monnap_fst$pairwise$fst,
                            deg_snps_in_transcriptome, 
                            rf_all$models$.base_.base$model$variable.importance
                            ))


colnames(rfdf) <- c("snp_pos", "FST", "DEG", "Importance")
rfdf$Importance <- as.numeric(rfdf$Importance)
rfdf$FST <- as.numeric(rfdf$FST)
rfdf$DEG  <- as.factor(rfdf$DEG)

rfdf <- rfdf[order(rfdf$Importance), ] # sort by importance 
saveRDS(rfdf, "rfdf.rds")

plot(rfdf$Importance, col= factor(rfdf$DEG), ylab="Importance in random forest", main="Importance of each SNP in random forest model (MON and NAP only)")

# plot importance value next to box plot of all snps and degs only 
p_rf_index <- ggplot(rfdf) + 
  geom_point(aes(x= 1:nrow(rfdf), y=Importance, color= DEG, shape=DEG)) + #* (blend("lighten") + blend("multiply", alpha = 0.5)) +
  scale_color_manual(values = c("black","darkred")) +
  theme_test() + 
  theme(legend.position="none") + 
  xlab('index')

p_rf_fst <- ggplot(rfdf) +
  geom_point(aes(x = FST, y = Importance, color=DEG), size=0.3) * (blend("lighten") + blend("multiply", alpha = 0.5)) +
  scale_color_manual(values = c("black", "red3"), labels = c(
    "Non-DEG SNPs", "DEG SNPs")) + 
  theme_test() +
  theme(legend.position = c(0.25, 0.75),
        legend.title = element_blank(), 
        text = element_text(size=16))

tiff("~/Downloads/", width=12, height=5, units = "in", res=300)
p
dev.off()
p_rf_fst + inset_element(p_rf_index, 
                                 left= 0.15, 
                                 bottom= 0.6, 
                                 right= 0.7,
                                 top= 0.9,
                                 align_to = "plot")

### from angsdsnp.r LFC VS FST plots
p_lfc_fst <- ggplot() + 
  geom_point(data=df, aes(x=fst, y=log2FoldChange), size = 0.3) + #* (blend("lighten") + blend("multiply", alpha = 0.5))+ 
  geom_smooth(data=df, method=lm, aes(x=fst, y=log2FoldChange), col="black") +
  geom_point(data=deg_df, aes(x=fst, y=log2FoldChange), size = 0.3, col="red3") +#* (blend("lighten") + blend("multiply", alpha = 0.5)) + 
  geom_smooth(data=deg_df, method=lm, aes(x=fst, y=log2FoldChange), col="red3") +
  theme_test() +  
  theme(text = element_text(size=16)) + labs(x="FST")

p_lfc_fst_text <- p_lfc_fst + 
  geom_text(aes(x = .29, y = -1.70, label= paste0("Transcriptome-wide SNPs, Slope: ", round(coef(lm(log2FoldChange ~ fst, data = df))["fst"], 2),", p < ", round(summary(lm(log2FoldChange ~ fst, data = df))$coefficients[2, 4], 2)))) + 
  geom_text(aes(x = .34, y = 1.65, label= paste0("DEG SNPs, Slope: ", round(coef(lm(log2FoldChange ~ fst, data = deg_df))["fst"], 2),", p < ", round(summary(lm(log2FoldChange ~ fst, data = deg_df))$coefficients[2, 4], 2))), color="red3") + theme_test() +  
  theme(text = element_text(size=16))







##### filter out Importance value < 0 

all_monnap <- calc_pairwise_fst(all_monnap, facets="site")
fst <- get.snpR.stats(all_monnap, "site", "fst")
fsts <- fst$pairwise

unique(rfdf$snp_pos[which(rfdf$Importance <= 0)])
keep <- unique(rfdf$snp_pos[which(rfdf$Importance <= 0)])
keep <- which(snp.meta(all_nocross)$ID %in% keep)
low_import <- all_nocross[keep,]
plot_clusters(low_import, facets="site", alt.palette=manual_col3) 

View(rfdf)




# violin plot below 
p1 <- ggplot(rfdf) +
  geom_violin(aes(x = DEG, y = Importance)) +
  geom_point(aes(x = DEG, y = Importance, color=DEG), position = position_jitter(seed = 1, width = 0.125), size=0.1) +
  scale_color_manual(values = c("black","red")) +
  scale_x_discrete(labels=c("0" = "neutral SNPs", "1" = "DEG SNPs")) + 
  theme_test()
