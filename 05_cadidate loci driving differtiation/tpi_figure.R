#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#+++ Explore TPI, a putatively adaptive gene   +++#
#+++ Code written by Andy Lee
#+++ Last Updated 7/10/2024
#+++ All rights reserved, contact andymuanlee@gmail.com  
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Load the packages
library(GenomeInfoDb)
library(Biostrings)
library(BSgenome)
library(ggplot2)
library(seqinr)
library(snpR)
library('DESeq2')
library("rtracklayer")
library("dplyr")

### ORF figure ====
# Load fasta and ORF sequence\
setwd("/Users/andy/KW/14_tpi")
contig  <- read.fasta("tpi_contig.fasta")

# write.csv(snp_data, "tpi_snp_SI.csv")
### Get FST and importance values for each snp 
tpi_rf <- readRDS("/Users/andy/KW/14_tpi/tpi_rf.rds") # from RandomForest Script 
tpi_rf$Position <- as.numeric(gsub(".+_(\\d+)$", "\\1", tpi_rf$snp_pos))
tpi_rf$Region <- NA
tpi_rf$Region[1:4] <- "Downstream"
tpi_rf$Region[5:8] <- "ORF"

tpi_rf_startend <- rbind(tpi_rf, matrix(NA, nrow = 2, ncol = NCOL(tpi_rf), 
                 dimnames = list(NULL, colnames(tpi_rf))))
tpi_rf_startend <- tpi_rf_startend[,c(6,7)]
tpi_rf_startend$Position[9:10] <- c(0, length(seq))


write.csv(tpi_rf, "tpi_snp_info.csv")

contig_sequence <- readDNAStringSet("tpi_contig.fasta")
seq <- DNAString(contig_sequence[1]$NODE_170811_length_1937_cov_17418.954399_g51510_i0)
orf <- contig[921:1664]

p_orf <- ggplot() + 
  geom_rect(aes(xmin = 0, xmax = length(seq), ymin = 0.45, ymax = 0.55), fill = "gray", alpha = 0.3) + 
  geom_rect(aes(xmin = 921, xmax = 1664,ymin = 0.45, ymax = 0.55), fill = "darkorange2", alpha = 0.5) +
  ## snp label 
  ggrepel::geom_text_repel(aes(x = Position, y = 0.551, label = 1:8), 
                           data = tpi_rf, 
                           direction = "x", 
                           hjust = 0.5,
                           vjust = 0,
                           # alpha=0.5, 
                           #nudge_x = .0001, 
                           nudge_y = .025, 
                           segment.alpha = 1, 
                           size =4) +
  ## position number of each snp 
  ggrepel::geom_text_repel(aes(x = Position, y = 0.45, label = Position), 
                           data = tpi_rf_startend, 
                           direction = "x", 
                           hjust = 0.5,
                           vjust = 1,
                           alpha=0.5, 
                           segment.alpha = 0.5,
                           nudge_y = -0.02,
                           size =3) +
  ## FST values 
  ggrepel::geom_text_repel(aes(x = Position, y = 0.41, label = round(FST, 2)), 
                           data = tpi_rf, 
                           direction = "x", 
                           hjust = 0.5,
                           alpha=0.5, 
                           nudge_y = -.02, 
                           segment.alpha = 0.5, 
                           size =3) +
  ## draw lines 
  geom_segment(aes(x = Position, xend = Position, y = 0.45, yend = 0.55, color = Region, linetype=Region), data = tpi_rf, size = .5) + 
  scale_color_manual(values=c("black", "darkorange2")) +
  scale_linetype_manual(values=c("dotted", "solid" )) +
  ## RF Importance 
  geom_point(data = tpi_rf, aes(x = Position, y = 0.41, size = Importance)) +
  theme_minimal() +
  scale_x_continuous(position="top") +
  xlab("Position") +
  theme(axis.text.y = element_blank(),    
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank()) 

pdf("~/KW/figures/figure 4 /p_orf.pdf", width = 7, height = 2.8)
p_orf
dev.off()



## Expression (normalized gene count)  ====
sample_info <- read.csv("~/KW/4_DESeq2/annotated_rnaspades/kw_sample_info.csv", row.names = 1)

dds.MON.NAP <- readRDS("~/KW/11_wgcna/dds.RDS")
assembly <- readGFF("~/KW/4_DESeq2/annotated_rnaspades/stringtie_all_merged_kw.gtf") 

gene_idx.MON.NAP      <- match(dds.MON.NAP@rowRanges@partitioning@NAMES, assembly$gene_id)
seq.name.MON.NAP      <- as.character(assembly$seqid[gene_idx.MON.NAP])
transcript.id.MON.NAP <- assembly$transcript_id[gene_idx.MON.NAP]
gene.id.MON.NAP       <- assembly$gene_id[gene_idx.MON.NAP]
xloc.MON.NAP          <- assembly$exon_number[gene_idx.MON.NAP]
contig_names.MON.NAP  <- as.data.frame(cbind(seq.name.MON.NAP, gene.id.MON.NAP, transcript.id.MON.NAP, xloc.MON.NAP))
dds.MON.NAP@rowRanges@partitioning@NAMES <- paste(contig_names.MON.NAP[,1])

normalized_counts <- counts(dds.MON.NAP, normalized = TRUE)
normalized_counts["NODE_170811_length_1937_cov_17418.954399_g51510_i0", ]

tpi_counts <- normalized_counts["NODE_170811_length_1937_cov_17418.954399_g51510_i0", ]
tpi_counts <- as.data.frame(tpi_counts)

tpi_counts$site <- NA
tpi_counts$site[1:30] <- "MON"
tpi_counts$site[31:54] <- "NAP"

summary_data <- tpi_counts %>%
  group_by(site) %>%
  summarize(mean_count = mean(tpi_counts),
            std_error = sd(tpi_counts) / sqrt(n()))

summary_data$site <- factor(summary_data$site, levels = c("NAP", "MON"))

p_normcount <- ggplot(data = summary_data, aes(x = site, y = mean_count)) +
  geom_bar(stat = "identity", aes(fill = site), alpha = 0.85) + 
  scale_fill_manual(values=c("#D1426FFF", "#000004FF")) +
  geom_errorbar(aes(ymin = mean_count - std_error, ymax = mean_count + std_error), width = 0.2, position = position_dodge(0.9)) +
  labs(x = "Site", y = "Mean Normalized Count") +
  theme_bw() + 
  theme(legend.position = "none")

pdf("~/KW/figures/figure 4 /p_nromcount.pdf", width = 4, height = 4)
p_normcount
dev.off()


### Allele freq of tpi ====
sample_meta <- read.csv("~/KW/09_angsd/neutral/sample_info.csv")
deg_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_degsnps.vcf", sample.meta = sample_meta)
deg_monnap <- filter_snps(deg_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")
snp.meta(deg_monnap)$CHROM <- gsub("\\.", "__", snp.meta(deg_monnap)$CHROM)

tpi_snps <- deg_monnap[CHROM = c("NODE_170811_length_1937_cov_17418__954399_g51510_i0")]

af <- tabulate_allele_frequency_matrix(tpi_snps, "site")
af_matrix <- get.snpR.stats(af, "site", "allele_frequency_matrix")

af_m <- t(as.data.frame(af_matrix$site$.base))
Gene = rownames(af_m)
Site = colnames(af_m)

af_df <- melt(af_m, id.vars = "Gene", variable.name = "Site")
colnames(af_df) <- c("Gene", "Site", "Frequency")
af_df$Site <- factor(af_df$Site, levels = c("NAP", "MON"))


nap_af <- af_df[af_df$Site=="NAP", ]
mon_df <- af_df[af_df$Site=="MON", ]

# find the minor allele in NAP 
# Initialize an empty dataframe to store the results for minor allele
# Iterate through the data frame
nap_maf <- data.frame(Gene = character(0), Site = character(0), Frequency = numeric(0))

for (i in seq(1, nrow(nap_af), by = 2)) {
  gene1 <- nap_af$Gene[i]
  gene2 <- nap_af$Gene[i + 1]
  freq1 <- nap_af$Frequency[i]
  freq2 <- nap_af$Frequency[i + 1]
  
  # Determine which gene has the smaller frequency
  if (freq1 <= freq2) {
    result_row <- nap_af[i, ]
  } else {
    result_row <- nap_af[i + 1, ]
  }
  
  # Append the result to the result dataframe
  nap_maf <- rbind(nap_maf, result_row)
}

# reconfigure the df
maf_comb <- merge(nap_maf, mon_df, by="Gene", all.x=TRUE, all.y=FALSE)
comb_mon <- maf_comb[, c(1,4,5)]
comb_nap <- maf_comb[, -c(4,5)]

colnames(comb_mon) <- c("Gene","Site", "Frequency")
colnames(comb_nap) <- c("Gene","Site", "Frequency")

comb <- rbind(comb_mon, comb_nap)
comb$gene_id <- as.numeric(gsub("[^0-9]", "", comb$Gene))
comb <- comb[order(comb$gene_id), ]
comb$id <- rep(8:1, each = 2 )

# add ORF info 
comb$orf <- NA
comb$orf[1:8] <- "ORF"
comb$orf[9:16]  <- "Downstream"

# make the plot: 
p_af <- ggplot(data = comb, aes(x = Site, y = Frequency, group = Gene, color=orf))  +  
  geom_point() +
  geom_line(aes(linetype=orf)) + #position = ggstance::position_dodgev(height = 0.02)) +
  ggrepel::geom_text_repel(data= subset(comb, Site %in% "MON"), 
                           aes(label=id),
                           direction = "y", 
                           hjust = 0, 
                           nudge_x = 0.2, 
                           color = "black", 
                           alpha=.5, size=3.5) +
  ylim (0, 1) +
  labs(x = "Site", y = "Allele Frequency") + 
  scale_color_manual(values=c("darkorange2","black")) +
  theme_bw() + 
  theme(legend.position = "none")

pdf("~/KW/figures/figure 4 /p_af.pdf", width = 4, height = 4)
p_af
dev.off()

### double check with some random alleles ====
sample <- deg_monnap[sample(nrow(deg_monnap), 1000, FALSE)]

af <- tabulate_allele_frequency_matrix(sample, "site")
af_matrix <- get.snpR.stats(af, "site", "allele_frequency_matrix")

data <- t(as.data.frame(af_matrix$site$.base))
Gene = rownames(data)
Site = colnames(data)

df <- melt(data, id.vars = "Gene", variable.name = "Site")
colnames(df) <- c("Gene", "Site", "Frequency")
df$Site <- factor(df$Site, levels = c("NAP", "MON"))

# find the minor allele in NAP 
nap_af <- df[df$Site=="NAP", ]
mon_df <- df[df$Site=="MON", ]


# Iterate through the data frame
nap_maf <- data.frame(Gene = character(0), Site = character(0), Frequency = numeric(0))
for (i in seq(1, nrow(nap_af), by = 2)) {
  gene1 <- nap_af$Gene[i]
  gene2 <- nap_af$Gene[i + 1]
  freq1 <- nap_af$Frequency[i]
  freq2 <- nap_af$Frequency[i + 1]
  
  # Determine which gene has the smaller frequency
  if (freq1 <= freq2) {
    result_row <- nap_af[i, ]
  } else {
    result_row <- nap_af[i + 1, ]
  }
  
  # Append the result to the result dataframe
  nap_maf <- rbind(nap_maf, result_row)
}

maf_comb <- merge(nap_maf, mon_df, by="Gene", all.x=TRUE, all.y=FALSE)

# reconfigure the df
comb_mon <- maf_comb[, c(1,4,5)]
comb_nap <- maf_comb[, -c(4,5)]
colnames(comb_mon) <- c("Gene","Site", "Frequency")
colnames(comb_nap) <- c("Gene","Site", "Frequency")
# comb <- rbind(comb_mon, comb_nap)
comb2 <- rbind(comb_mon, comb_nap)

ggplot()  +  
  geom_point(data = comb, aes(x = Site, y = Frequency, group = Gene)) +
  geom_line(data = comb, aes(x = Site, y = Frequency, group = Gene), color="darkblue", alpha=.4)  +
  geom_point(data = comb2, aes(x = Site, y = Frequency, group = Gene)) +
  geom_line(data = comb2, aes(x = Site, y = Frequency, group = Gene), color="darkred", alpha=.3) +
  ylim (0, 1.01) +
  labs(x = "Site", y = "Allele Frequency") + 
  theme_bw() + 
  theme(legend.position = "none") 

dev.off()

median(comb$Frequency)
median(comb2$Frequency)

hist(comb$Frequency, breaks=30)
hist(comb2$Frequency, breaks=30)

t.test(comb$Frequency, comb2$Frequency, conf.level = 0.95)
## calculate LD ====
### load raw vcf files 
sample_meta <- read.csv("~/KW/09_angsd/neutral/sample_info.csv")
all_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_allsnps.vcf", sample.meta = sample_meta)
all_monnap <- filter_snps(all_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")
snp.meta(all_monnap)$CHROM <- gsub("\\.", "__", snp.meta(all_monnap)$CHROM)
tpi_contig <- all_monnap[CHROM = "NODE_170811_length_1937_cov_17418__954399_g51510_i0"]

tpi_contig <- calc_pairwise_ld(tpi_contig, "CHROM", CLD = T, use.ME = T)
ld <- get.snpR.stats(tpi_contig, "CHROM", "LD")
View(ld$prox)

mean(ld$prox$CLD) #0.09308462
range(ld$prox$CLD) # 0.002262443 1.000000000

ld_tab <- ld$prox[,c(2,10,17,20:23)]
write.csv(ld_tab, "tpi_ld.csv")
