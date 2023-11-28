#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++ Random Forest to see how many SNPs explain the vairance among samples +++#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(snpR)
library(grid)
library(gridExtra)
library(ggplot2)
library(ggblend)
library(cowplot)
library(ggplotify)
setwd("~/KW/09_angsd")

sample_meta <- read.csv("~/KW/09_angsd/neutral/sample_info.csv")
deg_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_deg.vcf", sample.meta = sample_meta)
all_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_neutral.vcf", sample.meta = sample_meta)

deg_monnap <- filter_snps(deg_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")
all_monnap <- filter_snps(all_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")

###### Random Forest to see how many SNPs exlpain the vairance among samples ###############
rf_deg <- run_random_forest(deg_monnap, response = "site", importance = "permutation", num.trees =1000000)
rf_all <- run_random_forest(all_monnap, response = "site", importance = "permutation", num.trees = 1000000)

par(mar = c(1, 1, 1, 1))
plot(density(rf_deg$models$.base_.base$model$variable.importance))
plot(density(rf_all$models$.base_.base$model$variable.importance))

# See how well random forest models performed 
table(rf_deg$models$.base_.base$predictions)
table(rf_all$models$.base_.base$predictions)


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

median(all_monnap_fst$pairwise$fst) #0.02762004
median(deg_monnap_fst$pairwise$fst) #0.02382336

# double checking the dataframe 
rf_all$models$.base_.base$model$variable.importance
all_monnap_fst$pairwise$snp_pos <- paste0(all_monnap_fst$pairwise$CHROM, "_",all_monnap_fst$pairwise$position)
identical(all_monnap_fst$pairwise$snp_pos, all_monnap_snp_pos)


### create importance plot ====
rfdf <- as.data.frame(cbind(all_monnap_fst$pairwise$CHROM, 
                            all_monnap_fst$pairwise$snp_pos,
                            all_monnap_fst$pairwise$fst,
                            deg_snps_in_transcriptome, 
                            rf_all$models$.base_.base$model$variable.importance
                            ))

colnames(rfdf) <- c("contig", "snp_pos", "FST", "DEG", "Importance")
rfdf$Importance <- as.numeric(rfdf$Importance)
rfdf$FST <- as.numeric(rfdf$FST)
rfdf$DEG  <- as.factor(rfdf$DEG)

## Investigating the outlier 

outlier_contig <- snp.meta(all_monnap)[snp.meta(all_monnap)$CHROM == "NODE_170811_length_1937_cov_17418.954399_g51510_i0", ]

outlier_contig_fst <- rfdf[rfdf$contig == "NODE_170811_length_1937_cov_17418.954399_g51510_i0", ]




#saveRDS(rfdf, "rfdf.rds")
# rfdf <- readRDS("rfdf.rds")

# saveRDS(outlier_contig_fst, "/Users/andy/KW/14_tpi/tpi_rf.rds")

## COMMAND to get sequence to blast to find gene name of outlier 
# grep -A1 ">NODE_170811_length_1937_cov_17418.954399_g51510_i0" rnaspades_assembly_calpoly_annotated_swissprot_nr.fasta

# triosephosphate isomerase


# plot(rfdf$Importance, col= factor(rfdf$DEG), ylab="Importance in random forest", main="Importance of each SNP in random forest model (MON and NAP only)")

# plot importance value next to box plot of all snps and degs only 
rfdf <- rfdf[order(rfdf$Importance),]

p_rf_index <- ggplot(rfdf) + 
  geom_point(aes(x= 1:nrow(rfdf), y=Importance, color= DEG, shape=DEG), size=0.6) + #* (blend("lighten") + blend("multiply", alpha = 0.4))+
  scale_color_manual(values = c("black","darkorange2")) +
  theme_test() + 
  theme(legend.position="none",
        text = element_text(size=13)) + 
  xlab('Rank order')

p_rf_fst <- ggplot(rfdf) +
  geom_point(data=subset(rfdf, DEG %in% FALSE), aes(x = FST, y = Importance, color=DEG), size=0.5, color = "black") + #* (blend("lighten") + blend("multiply", alpha = 0.4)) +
  geom_point(data=subset(rfdf, DEG %in% TRUE), aes(x = FST, y = Importance, color=DEG), size=0.5, color = "darkorange2") +#  *(blend("lighten") + blend("multiply", alpha = 0.4))+
  #scale_color_manual(values = c("black", "darkorange2"), labels = c("Non-DEG SNPs", "DEG SNPs")) + 
  theme_test() +
  theme(legend.position = c(0.25, .35),
        legend.title = element_blank(), 
        text = element_text(size=18)) + 
  guides(colour = guide_legend(override.aes = list(size=5))) + # change guide size
  xlab(expression(italic(F[ST]))) + 
  geom_segment(x=median(all_monnap_fst$pairwise$fst), # add a segment for median fst
               xend=median(all_monnap_fst$pairwise$fst), 
               y=0, yend=0.00045, 
               color="black", linetype="dashed", linewidth=0.2) + 
  geom_segment(x=median(deg_monnap_fst$pairwise$fst), 
               xend=median(deg_monnap_fst$pairwise$fst), 
               y=0, yend=0.00045, 
               color="darkorange2", linetype="dashed", linewidth=0.2) 

p_rf <- p_rf_fst + 
  annotation_custom(
    grob = ggplotGrob(p_rf_index),
    xmin = -Inf,
    xmax = 0.235,
    ymin = .00095,
    ymax = .0024
  )


  

# t_grob <- grobTree(textGrob("Triosephosphate Isomerase", x = 0.7, y= 0.95, hjust=0, gp=gpar(col="black", fontsize=8, fontface="italic")))
#p_rf + annotation_custom(t_grob)

### from angsdsnp.r LFC VS FST plots
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

# p <- ggplot(df) +
#    geom_point(aes(x=fst, y=log2FoldChange), size = 0.4) + 
#    geom_smooth(method=lm, aes(x=fst, y=log2FoldChange)) +
#    ggtitle("Transcriptome-wide SNPs")


# deg snps only 
deg_df <- merge(z, y, by="CHROM")
# p2 <- ggplot(deg_df) +
#    geom_point(aes(x=fst, y=log2FoldChange), size = 0.4) + 
#    geom_smooth(method=lm, aes(x=fst, y=log2FoldChange)) +
#    ggtitle("DEG SNPs")

### both together
p_lfc_fst <- ggplot() + 
  geom_point(data=df, aes(x=fst, y=log2FoldChange), size = 0.4) + #* (blend("lighten") + blend("multiply", alpha = 0.4))+ 
  geom_smooth(data=df, method=lm, aes(x=fst, y=log2FoldChange), col="black") +
  geom_point(data=deg_df, aes(x=fst, y=log2FoldChange), size = 0.4, col="darkorange2") +#  * (blend("lighten") + blend("multiply", alpha = 0.4)) + 
  geom_smooth(data=deg_df, method=lm, aes(x=fst, y=log2FoldChange), col="darkorange2") +
  theme_test() +  
  theme(text = element_text(size=18)) + xlab(expression(italic(F[ST])))

p_lfc_fst_text <- p_lfc_fst + 
  geom_text(aes(x = .125, y = 2.75, label= paste0("Non-DEG SNPs: Slope = ", round(coef(lm(log2FoldChange ~ fst, data = df))["fst"], 5),", p < 0.00001 ")), hjust=0) + 
  geom_text(aes(x = .125, y = 3, label= paste0("DEG SNPs: Slope = ", round(coef(lm(log2FoldChange ~ fst, data = deg_df))["fst"], 2),", p < 0.00001 "), hjust=0), color="darkorange2") + theme_test() +  
  theme(text = element_text(size=18))

# calculating p-value for each regression line 
summary(lm(log2FoldChange ~ fst, data = df))$coefficients[2, 4]
summary(lm(log2FoldChange ~ fst, data = deg_df))$coefficients[2, 4]

gap::chow.test(y1 = deg_df$log2FoldChange, x1 = deg_df$fst, y2=df$log2FoldChange, x2 = df$fst)


### Putting it all together :) ====
#or 
cairo_pdf("~/KW/figures/figure 3/rf.pdf", width= 7.25 , height= 8)
#svg("~/KW/figures/figure 3/rf.svg", width= 7.25 , height= 8)
plot_grid(
  p_lfc_fst_text,
  p_rf,
  ncol = 1, 
  rel_heights= c(3,3), 
  align = "v"
)
dev.off()


### FST Distribution/Density/Box plots for supplemental ====
# Create a grouping variable for each dataset
deg_fst <- deg_monnap_fst$pairwise
all_fst <- all_monnap_fst$pairwise

all_fst$loci <- "Non-DEG SNPs"
deg_fst$loci <- "DEG SNPs"

# Combine the datasets
combined_data <- rbind(all_fst, deg_fst)

# box plot
p_fst_box <- ggplot(data = combined_data, aes(y = loci, x = fst, color = loci)) +
  geom_boxplot() +
  scale_color_manual(values=c("darkorange2", "black")) +
  theme_test() + 
  xlab(expression(italic(F[ST]))) + 
  theme(axis.title.y = element_blank(), legend.position = "bottom") + 
  guides(color=guide_legend(title="Loci Type"))

# density plot 
p_fst_den <- ggplot() + 
  geom_density(data=all_fst, aes(x=fst), adjust = 5) + 
  geom_density(data=deg_fst, aes(x=fst), col="darkorange2", adjust = 5) +
  theme_test() + 
  xlab(expression(italic(F[ST]))) +
  ylab('Density') +
  theme(text = element_text(size=16))


pdf("~/KW/figures/supplemental /fst_distribution.pdf", 7.25 , 11)
plot_grid(
  p_fst_box,
  p_fst_den, 
  ncol = 1, 
  rel_heights= c(1,.50), 
  align = "v", 
  labels = "AUTO"
)
dev.off()

# violin plot ====
p1 <- ggplot(rfdf) +
  geom_violin(aes(x = DEG, y = Importance)) +
  geom_point(aes(x = DEG, y = Importance, color=DEG), position = position_jitter(seed = 1, width = 0.125), size=0.1) +
  scale_color_manual(values = c("black","red")) +
  scale_x_discrete(labels=c("0" = "neutral SNPs", "1" = "DEG SNPs")) + 
  theme_test()
 # dist plot ====

# 
bins <- seq(-0.1, 0.5, by=0.01)
deg_fst <- deg_monnap_fst$pairwise
p_degsnps <- ggplot(deg_fst) +
  geom_histogram(aes(x=fst), breaks=bins) + 
  ggtitle("DEG snps FST distribution")

all_fst <- all_monnap_fst$pairwise

p_allsnps <- ggplot(all_fst) +
  geom_histogram(aes(x=fst), breaks=bins) + 
  ggtitle("Transcriptome-wide snps FST distribution")