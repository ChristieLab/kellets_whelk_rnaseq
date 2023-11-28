

# Install other required packages
# install.packages("ggplot2")
# install.packages("seqinr")
# BiocManager::install("Biostrings")
# BiocManager::install("BSgenome")
# BiocManager::install("GenomeInfoDb")

# Load the packages
library(GenomeInfoDb)
library(Biostrings)
library(BSgenome)
library(ggplot2)
library(seqinr)

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

write.csv(tpi_rf, "tpi_snp_info.csv")

contig_sequence <- readDNAStringSet("tpi_contig.fasta")
seq <- DNAString(contig_sequence[1]$NODE_170811_length_1937_cov_17418.954399_g51510_i0)
orf <- contig[921:1664]

p_orf <- ggplot() + 
  geom_rect(aes(xmin = 1, xmax = length(seq), ymin = 0.25, ymax = 0.75), fill = "gray", alpha = 0.3) + 
  geom_rect(aes(xmin = 537, xmax = 1583, ymin = 0.25, ymax = 0.75), fill = "darkorange2", alpha = 0.5) +
  geom_text(aes(x = 770, y = 0.5, label = "TPI ORF"), color = "black", size = 3) +
  ggrepel::geom_text_repel(aes(x = Position, y = 0.751, label = 1:8), 
                            data = tpi_rf, 
                            direction = "x", 
                            hjust = 0.5,
                            alpha=0.5, 
                            #nudge_x = .0001, 
                            nudge_y = 0.085, 
                            segment.alpha = 0.5, 
                            size =3) +
  geom_segment(aes(x = Position, xend = Position, y = 0.25, yend = 0.75, color = Region, linetype=Region), data = tpi_rf, size = .5) + 
  scale_color_manual(values=c("black", "darkorange2")) +
  scale_linetype_manual(values=c("dotted", "solid" )) +
  geom_point(data = tpi_rf, aes(x = Position, y = 0.20, size = Importance)) +
  theme_minimal() +
  ggrepel::geom_text_repel(aes(x = Position, y = 0.20, label = round(FST, 2)), 
                           data = tpi_rf, 
                           direction = "x", 
                           hjust = 0.5,
                           alpha=0.5, 
                           nudge_y = -0.085, 
                           segment.alpha = 0.5, 
                           size =3) +
  xlab("Position") +
  theme(axis.text.y = element_blank(),    
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank(),
        panel.grid = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank()) 

pdf("~/KW/figures/figure 4 /p_orf.pdf", width = 7, height = 3)
p_orf
dev.off()


## Expression (normalized gene count)  ====
library('DESeq2')
library("rtracklayer")
library("dplyr")

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
deg_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_deg.vcf", sample.meta = sample_meta)
deg_monnap <- filter_snps(deg_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")
snp.meta(deg_monnap)$CHROM <- gsub("\\.", "__", snp.meta(deg_monnap)$CHROM)

tpi_snps <- deg_monnap[CHROM = c("NODE_170811_length_1937_cov_17418__954399_g51510_i0")]


genotypes <- format_snps(tpi_snps, output="NN")
geno_table <- t(genotypes[,c(1,2,4,5,8:ncol(genotypes))])
#write.csv(geno_table, "~/Downloads/kw_tpi_geno.csv")


af <- tabulate_allele_frequency_matrix(tpi_snps, "site")
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
comb$id <- rep(1:8, each = 2 )

# add ORF info 
comb$orf <- NA
comb$orf[9:16] <- "ORF"
comb$orf[1:8]  <- "Downstream"

comb
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
  ylim (0, 0.6) +
  labs(x = "Site", y = "Allele Frequency") + 
  scale_color_manual(values=c("darkorange2","black")) +
  theme_bw() + 
  theme(legend.position = "none")

pdf("~/KW/figures/figure 4 /p_af.pdf", width = 4, height = 4)
p_af
dev.off()

# ggplot(data = subset_df, aes(x = Site, y = Frequency, group = Gene, color = Gene)) +
#   geom_line() +
#   geom_point() + 
#   ylim (0, 1) +
#   labs(x = "Site", y = "Frequency") +
#   theme_minimal()

# plot_clusters(tpi_snps, facets="site")

pdf("~/KW/figures/figure 4 /p_orf.pdf", width = 6, height = 4)
p_orf
p_af
p_normcount

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

### metacoder plot ====
### get top snps to make metacoder plot
all_monnap <- calc_pairwise_fst(all_monnap, facets="site")
snp.meta(all_monnap)$ID <- paste0(snp.meta(all_monnap)$CHROM, "_",snp.meta(all_monnap)$position)

all_monnap <- calc_pairwise_fst(all_monnap, facets="site")
all_monnap_fst <-  get.snpR.stats(all_monnap, "site", "fst")
fsts <- all_monnap_fst$pairwise

keep <- unique(fsts$ID[which(fsts$fst >= 0.3)])
keep <- which(snp.meta(all_monnap)$ID %in% keep)
high_fst  <- all_monnap[keep,]
unique(snp.meta(high_fst)$CHROM) # one is not unique 

snp.meta(high_fst)$CHROM <- gsub("\\__", "\\.", snp.meta(high_fst)$CHROM)
out <- paste0(">", snp.meta(high_fst)$CHROM)

write.table(out, "~/Downloads/highfst03_contignames.txt",row.names=FALSE,sep="\t", quote = FALSE) #output contigs for eggnog annotation 

  ### run enrichment using eggnog output and 
library(topGO)

allgeneID2GO <- readMappings(file = "/Users/andy/KW/14_tpi/all_ref_genes_to_go.txt")

# interestID2GO <- readMappings(file = "/Users/andy/KW/14_tpi/genes_to_go.txt")
interestID2GO <- readMappings(file = "/Users/andy/KW/14_tpi/genes_03_to_go.txt")
str(geneID2GO)

geneNames <- names(allgeneID2GO)
myInterestingGenes <- names(interestID2GO)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

GOdata <- new("topGOdata", 
              description = "Kellet's Whelks most differentiated genes",
              ontology = "BP",
              allGenes = geneList, 
              annot = annFUN.gene2GO, 
              gene2GO = allgeneID2GO)

genes(GOdata) ## list of genes 
numGenes(GOdata) # number of genes with GO Terms

num.ann.genes <- countGenesInTerm(GOdata) ## the number of annotated genes
num.ann.genes # get counts of go terms 
ann.genes <- genesInTerm(GOdata) ## get the annotations
head(ann.genes) # genes with each go term 

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultweight01 <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")

allRes <- GenTable(GOdata, classicFisher = resultFisher, weight01Fisher = resultweight01, 
                   orderBy = "weight01Fisher", ranksOf = "classicFisher", topNodes = 20)

allRes$weight01Fisher <-as.numeric(allRes$weight01Fisher)
allRes <- allRes[allRes$weight01Fisher < 0.001, ]

allRes
write.csv(allRes, "GO_FST_03.csv")
### Visualize in metacoder
library(GO.db)
library(metacoder)
terms <- AnnotationDbi::select(GO.db, allRes$GO.ID, columns = c('TERM','ONTOLOGY'), keytype = "GOID") 
bp.terms <- terms[terms$ONTOLOGY=="BP", ] # grab only BP terms 
bp.terms <- bp.terms[!is.na(bp.terms), ]

gobpparents <- as.list(GOBPPARENTS)
gobpparents <- gobpparents[!is.na(gobpparents)]

get_parents = function(x, current=x, all_paths = FALSE, verbose = TRUE, valid_relationships = c("isa")) {
  # Get immediate children of current taxon
  parents = tryCatch({
    possible_parents <- as.list(gobpparents[x[1]])[[1]] #have to manually change this line 
    if (! is.null(valid_relationships)) {
      possible_parents <- possible_parents[names(possible_parents) %in% valid_relationships]
    }
    names(AnnotationDbi::Term(possible_parents))
  }, error = function(e) {
    c()
  })
  # only go down one path if desired
  if (! all_paths) {
    parents <- parents[1]
  }
  parents <- parents[parents != "all"]
  if (is.null(parents)) {
    return(c())
  } else if (length(parents) == 0) {
    cat(length(x))
    return(paste0(collapse = ";", AnnotationDbi::Term(x)))
  } else {
    next_x <- lapply(parents, function(y) c(y, x))
    
    # Run this function on them to get their output
    child_output <- lapply(next_x, get_parents, all_paths = all_paths)
    output <- unlist(child_output, recursive = FALSE)
    
    return(output)
  }
}
allRes$weight01Fisher <-as.numeric(allRes$weight01Fisher)
  
sig.GOID <- allRes$GO.ID[allRes$weight01Fisher < 0.001]

bpterms = lapply(sig.GOID, get_parents, all_paths = FALSE)
bpres   = data.frame(class=unlist(bpterms))

write.table(bpres, "/Users/andy/KW/14_tpi/temp.csv", sep=",", col.names=TRUE, row.names=FALSE)
bpres=read.table("/Users/andy/KW/14_tpi/temp.csv", header=TRUE, sep=",")

taxdata <- parse_tax_data(bpres, class_sep = ";")

set.seed(420) ## seeds: red = 50, thistle = 100, 105 for DEGs 

heat_tree(taxdata, 
          node_label = taxon_names, ### modified by AL 
          # node_size = colsandsize$ngenes,
          # node_size_trans = "log10",
          node_size_range = c(0.01, 0.01),
          # node_label_size_trans = "log10",
          node_label_size_range = c(0.01, 0.01),
          # edge_size_trans = "log10",
          edge_size_range = c(0.004, 0.004),
          node_color = "royalblue",
          # node_color_trans = "linear",
          # node_color_range = diverging_palette(),
          # node_color_interval = c(-4, 4),
          # edge_color_trans = "linear",
          # edge_color_range = diverging_palette(),
          # edge_color_interval =  c(-4, 4),
          node_label_max = 500,
          # node_color_axis_label = "Factor change",
          # node_size_axis_label = "Number of genes",
          layout = "da", initial_layout = "re"
)

