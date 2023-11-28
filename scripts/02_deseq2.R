#---------------- Data Description -------------------------------------------------------------------------------------#
# script written by mark christie on 11/25/2019; contact at markchristie1500@gmail.com
# Rversion:
# DESeq2 version:
# script modified from Avril Harder, available at: https://github.com/ChristieLab/Salmo_salar_RNAseq
# data consists of hood river steelhead reads aligned with hisat2 > stringtie > htseq
# covariates include sex, cross date, and sequencing machine
# relevant sample data can be found in "sample_info_all.csv"

#set working directory, load libraries, list files
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# BiocManager::install("rtracklayer")
# BiocManager::install("gplots")

library('DESeq2')
library("rtracklayer")
library("ggplot2")
library("gplots")


setwd("~/KW/4_DESeq2/annotated_rnaspades")
list.files()

#---------------- Data loading and prep --------------------------------------------------------------------------------#
# read in count data generated with featureCounts
# single end and paired end reads aligned to same merged gtf file "stringtie_all_merged.gtf"
# gtfs used to create merged file found in "mergelist.txt"

all.gene.counts <- read.table("featureCount_counts.txt", sep="\t", header=T, row.names=1, check.names=F)

#may not want to drop below to find gene locations
drop.vars   <- names(all.gene.counts) %in% c("Chr","Start", "End", "Strand", "Length")
all.gene.counts <- all.gene.counts[!drop.vars]
head(all.gene.counts)

samples <- read.csv("kw_sample_info.csv", row.names = 1)
samples


#---------------- Sample processing ------------------------------------------------------------------------------------#
# subset by "population". we are mainly interested in MON x NAP, but will make pairwise comparisions between all combinations

POL <- which(samples[, 2] == "POL") # might want to combine POL with NAP later
NAPxMON <- which(samples[, 2] == "NAPxMON")
MON <- which(samples[, 2] == "MON")
NAP <- which(samples[ ,2] == "NAP")


# gene.counts.MON.NAP <- all.gene.counts[, c(MON, NAP)]
gene.counts.MON.POL <- all.gene.counts[, c(MON, POL)]
# gene.counts.NAP.POL <- all.gene.counts[, c(NAP, POL)]
# gene.counts.NAPxMON.MON <- all.gene.counts[, c(NAPxMON, MON)]
# gene.counts.NAPxMON.NAP <- all.gene.counts[, c(NAPxMON, NAP)]

# read in CSV with sample information
# col.data.MON.NAP <- samples[c(MON, NAP),]
col.data.MON.POL <- samples[c(MON, POL),]
# col.data.NAP.POL <- samples[c(NAP, POL),]
# col.data.NAPxMON.MON <- samples[c(NAPxMON, MON),]
# col.data.NAPxMON.NAP <- samples[c(NAPxMON, NAP),]

# check to make sure that sample names are in the same order in the gene count table and the sample info table; must be TRUE or do not proceed
# all(rownames(col.data.MON.NAP) %in% colnames(gene.counts.MON.NAP))
# all(rownames(col.data.MON.NAP) == colnames(gene.counts.MON.NAP))

all(rownames(col.data.MON.POL) %in% colnames(gene.counts.MON.POL))
all(rownames(col.data.MON.POL) == colnames(gene.counts.MON.POL))
# 
# all(rownames(col.data.NAP.POL) %in% colnames(gene.counts.NAP.POL))
# all(rownames(col.data.NAP.POL) == colnames(gene.counts.NAP.POL))
# 
# all(rownames(col.data.NAPxMON.MON) %in% colnames(gene.counts.NAPxMON.MON))
# all(rownames(col.data.NAPxMON.MON) == colnames(gene.counts.NAPxMON.MON))
# 
# all(rownames(col.data.NAPxMON.NAP) %in% colnames(gene.counts.NAPxMON.NAP))
# all(rownames(col.data.NAPxMON.NAP) == colnames(gene.counts.NAPxMON.NAP))


#----------------------- Construction of DESeqDataSet object and DE analysis ---------------------------------------------------
# test for effects of treatment while controlling for the effect of family; dds = DESeqDataSet
# variable order matters! Covariates go first, variable of interest goes last in design 

# dds.MON.NAP <- DESeqDataSetFromMatrix(countData = gene.counts.MON.NAP,
#                               colData = col.data.MON.NAP,
#                               design = ~ as.factor(site))
# dds.MON.NAP$site <- relevel(dds.MON.NAP$site, "NAP") ## set the untreated group as the control
# dds.MON.NAP <- DESeq(dds.MON.NAP) # DESeq() = differential expression analysis based on the negative binomial distribution
# plotDispEsts(dds.MON.NAP, ylim=c(1e-6, 1e1)) # examine distribution of dispersion values

# # MON x POL 
dds.MON.POL <- DESeqDataSetFromMatrix(countData = gene.counts.MON.POL,
                                      colData = col.data.MON.POL,
                                      design = ~ as.factor(site))
dds.MON.POL$site <- relevel(dds.MON.POL$site, "POL") 
dds.MON.POL <- DESeq(dds.MON.POL) 
plotDispEsts(dds.MON.POL, ylim=c(1e-6, 1e1))   

# dds.NAP.POL <- DESeqDataSetFromMatrix(countData = gene.counts.NAP.POL,
#                                       colData = col.data.NAP.POL,
#                                       design = ~ as.factor(site))
# dds.NAP.POL$site <- relevel(dds.NAP.POL$site, "POL") 
# dds.NAP.POL <- DESeq(dds.NAP.POL) 
# plotDispEsts(dds.NAP.POL, ylim=c(1e-6, 1e1))   
# 
# dds.NAPxMON.MON <- DESeqDataSetFromMatrix(countData = gene.counts.NAPxMON.MON,
#                                       colData = col.data.NAPxMON.MON,
#                                       design = ~ as.factor(site))
# dds.NAPxMON.MON$site <- relevel(dds.NAPxMON.MON$site, "MON") 
# dds.NAPxMON.MON <- DESeq(dds.NAPxMON.MON) 
# plotDispEsts(dds.NAPxMON.MON, ylim=c(1e-6, 1e1))   
# 
# dds.NAPxMON.NAP <- DESeqDataSetFromMatrix(countData = gene.counts.NAPxMON.NAP,
#                                       colData = col.data.NAPxMON.NAP,
#                                       design = ~ as.factor(site))
# dds.NAPxMON.NAP$site <- relevel(dds.NAPxMON.NAP$site, "NAP") 
# dds.NAPxMON.NAP <- DESeq(dds.NAPxMON.NAP) 
# plotDispEsts(dds.NAPxMON.NAP, ylim=c(1e-6, 1e1))  

# dds.MON.NAP
dds.MON.POL
# dds.NAP.POL
# dds.NAPxMON.MON
# dds.NAPxMON.NAP

#-----------------------------------------------------------------------------------------------------------------------#
##----------------------- match gene names  -----------------------------_# 
## renames genes in a DESeqDataSet to match gene names in the reference annotation (i.e., replaces MSTRG herever gene name is known)
## create congig_names, which is a table linking gene name, transcript ID, and xloc information for later istinguishing different isoforms of the same gene (same gene start and end) and transcripts generated from ifferent paralogs of the same gene (by xloc)


assembly <- readGFF("stringtie_all_merged_kw.gtf")             ## read in merged GTF

# gene_idx.MON.NAP      <- match(dds.MON.NAP@rowRanges@partitioning@NAMES, assembly$gene_id)
# seq.name.MON.NAP      <- as.character(assembly$seqid[gene_idx.MON.NAP])
# transcript.id.MON.NAP <- assembly$transcript_id[gene_idx.MON.NAP]
# gene.id.MON.NAP       <- assembly$gene_id[gene_idx.MON.NAP]
# xloc.MON.NAP          <- assembly$exon_number[gene_idx.MON.NAP]
# contig_names.MON.NAP  <- as.data.frame(cbind(seq.name.MON.NAP, gene.id.MON.NAP, transcript.id.MON.NAP, xloc.MON.NAP))
# dds.MON.NAP@rowRanges@partitioning@NAMES <- paste(contig_names.MON.NAP[,1],contig_names.MON.NAP[,2],contig_names.MON.NAP[ ,3], sep="|")

gene_idx.MON.POL      <- match(dds.MON.POL@rowRanges@partitioning@NAMES, assembly$gene_id)
seq.name.MON.POL      <- as.character(assembly$seqid[gene_idx.MON.POL])
transcript.id.MON.POL <- assembly$transcript_id[gene_idx.MON.POL]
gene.id.MON.POL       <- assembly$gene_id[gene_idx.MON.POL]
xloc.MON.POL          <- assembly$exon_number[gene_idx.MON.POL]
contig_names.MON.POL  <- as.data.frame(cbind(seq.name.MON.POL, gene.id.MON.POL, transcript.id.MON.POL, xloc.MON.POL))
dds.MON.POL@rowRanges@partitioning@NAMES <- paste(contig_names.MON.POL[,1],contig_names.MON.POL[,2],contig_names.MON.POL[,3], sep="|")

# gene_idx.NAP.POL      <- match(dds.NAP.POL@rowRanges@partitioning@NAMES, assembly$gene_id)
# seq.name.NAP.POL      <- as.character(assembly$seqid[gene_idx.NAP.POL])
# transcript.id.NAP.POL <- assembly$transcript_id[gene_idx.NAP.POL]
# gene.id.NAP.POL       <- assembly$gene_id[gene_idx.NAP.POL]
# xloc.NAP.POL          <- assembly$exon_number[gene_idx.NAP.POL]
# contig_names.NAP.POL  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# dds.NAP.POL@rowRanges@partitioning@NAMES <- paste(contig_names.NAP.POL[,1],contig_names.NAP.POL[,2],contig_names.NAP.POL[,3], sep="|")

# gene_idx.NAPxMON.MON      <- match(dds.NAPxMON.MON@rowRanges@partitioning@NAMES, assembly$gene_id)
# seq.name.NAPxMON.MON      <- as.character(assembly$seqid[gene_idx.NAPxMON.MON])
# transcript.id.NAPxMON.MON <- assembly$transcript_id[gene_idx.NAPxMON.MON]
# gene.id.NAPxMON.MON       <- assembly$gene_id[gene_idx.NAPxMON.MON]
# xloc.NAPxMON.MON          <- assembly$exon_number[gene_idx.NAPxMON.MON]
# contig_names.NAPxMON.MON  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# dds.NAPxMON.MON@rowRanges@partitioning@NAMES <- paste(contig_names.NAPxMON.MON[,1],contig_names.NAPxMON.MON[,2],contig_names.NAPxMON.MON[,3], sep="|")

# gene_idx.NAPxMON.NAP      <- match(dds.NAPxMON.NAP@rowRanges@partitioning@NAMES, assembly$gene_id)
# seq.name.NAPxMON.NAP      <- as.character(assembly$seqid[gene_idx.NAPxMON.NAP])
# transcript.id.NAPxMON.NAP <- assembly$transcript_id[gene_idx.NAPxMON.NAP]
# gene.id.NAPxMON.NAP       <- assembly$gene_id[gene_idx.NAPxMON.NAP]
# xloc.NAPxMON.NAP          <- assembly$exon_number[gene_idx.NAPxMON.NAP]
# contig_names.NAPxMON.NAP  <- as.data.frame(cbind(seq.name, gene.id, transcript.id, xloc))
# dds.NAPxMON.NAP@rowRanges@partitioning@NAMES <- paste(contig_names.NAPxMON.NAP[,1],contig_names.NAPxMON.NAP[,2],contig_names.NAPxMON.NAP[,3], sep="|")

# which(duplicated(dds.MON.NAP@rowRanges@partitioning@NAMES)) # check to make sure that no gene names in dds are duplicates
which(duplicated(dds.MON.POL@rowRanges@partitioning@NAMES)) 
# which(duplicated(dds.NAP.POL@rowRanges@partitioning@NAMES)) 
# which(duplicated(dds.NAPxMON.MON@rowRanges@partitioning@NAMES))
# which(duplicated(dds.NAPxMON.NAP@rowRanges@partitioning@NAMES))

#----------------- Extract and examine results from DESeq analysis -----------------------------------------------------#
# res.MON.NAP     <- results(dds.MON.NAP    , alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
res.MON.POL     <- results(dds.MON.POL    , alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
# res.NAP.POL     <- results(dds.NAP.POL    , alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
# res.NAPxMON.MON <- results(dds.NAPxMON.MON, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)
# res.NAPxMON.NAP <- results(dds.NAPxMON.NAP, alpha=0.05, pAdjustMethod = "BH", independentFiltering = TRUE)

# length(which(res.MON.NAP[, 6]     < 0.05)) #2770
length(which(res.MON.POL[, 6]     < 0.05)) #39522
# length(which(res.NAP.POL[, 6]     < 0.05)) 
# length(which(res.NAPxMON.MON[, 6] < 0.05)) 
# length(which(res.NAPxMON.NAP[, 6] < 0.05)) 

# summary(res.MON.NAP    )
summary(res.MON.POL    )
# summary(res.NAP.POL    )
# summary(res.NAPxMON.MON)
# summary(res.NAPxMON.NAP)

## examine MA plot (normalized counts vs. log2fold changes)
# plotMA(res.MON.NAP    , ylim = c(-6, 6), xlim=c(1e-2,1e6)) 
plotMA(res.MON.POL    , ylim = c(-6, 6), xlim=c(1e-2,1e6))
# plotMA(res.NAP.POL    , ylim = c(-6, 6), xlim=c(1e-2,1e6))
# plotMA(res.NAPxMON.MON, ylim = c(-6, 6), xlim=c(1e-2,1e6))
# plotMA(res.NAPxMON.NAP, ylim = c(-6, 6), xlim=c(1e-2,1e6))

## examine histograms of p-values and adjusted p-values
# hist(res.MON.NAP$pvalue, breaks=20, col="grey")
# hist(res.MON.NAP$padj, breaks=20, col="grey") 

hist(res.MON.POL$pvalue, breaks=20, col="grey")
hist(res.MON.POL$padj, breaks=20, col="grey") 
# 
# hist(res.NAP.POL$pvalue, breaks=20, col="grey")
# hist(res.NAP.POL$padj, breaks=20, col="grey") 
# 
# hist(res.NAPxMON.MON$pvalue, breaks=20, col="grey")
# hist(res.NAPxMON.MON$padj, breaks=20, col="grey")
# 
# hist(res.NAPxMON.NAP$pvalue, breaks=20, col="grey")
# hist(res.NAPxMON.NAP$padj, breaks=20, col="grey")

#-----------------------------------------------------------------------------------------------------------------------#

################### DATA MANIPULATION AND EXPLORATION ############################
#--<>--<>--<>--<>-- Log transformation and distance calculation --<>--<>--<>--<>--
## rlog() = transforms count data to the log2 scale in a way that minimizes differences between
## samples for genes with small counts, and which normalizes with respect to library size
#rld <- rlog(dds)
# vst is faster alternative
# head(assay(rld.MON.NAP))
# rld.MON.NAP     <- vst(dds.MON.NAP)
rld.MON.POL     <- vst(dds.MON.POL)
# rld.NAP.POL     <- vst(dds.NAP.POL)
# rld.NAPxMON.MON <- vst(dds.NAPxMON.MON)
# rld.NAPxMON.NAP <- vst(dds.NAPxMON.NAP)

## calculate Euclidean distances between all samples to examine overall similarity
# sampleDists.MON.NAP     <- dist(t(assay(rld.MON.NAP)))
sampleDists.MON.POL     <- dist(t(assay(rld.MON.POL)))
# sampleDists.NAP.POL     <- dist(t(assay(rld.NAP.POL)))
# sampleDists.NAPxMON.MON <- dist(t(assay(rld.NAPxMON.MON)))
# sampleDists.NAPxMON.NAP <- dist(t(assay(rld.NAPxMON.NAP)))

# sampleDistMatrix.MON.NAP     <- as.matrix(sampleDists.MON.NAP)    
sampleDistMatrix.MON.POL     <- as.matrix(sampleDists.MON.POL)    
# sampleDistMatrix.NAP.POL     <- as.matrix(sampleDists.NAP.POL)    
# sampleDistMatrix.NAPxMON.MON <- as.matrix(sampleDists.NAPxMON.MON)
# sampleDistMatrix.NAPxMON.NAP <- as.matrix(sampleDists.NAPxMON.NAP)

#-----------------------------------------------------------------------------------------------------------------------#
## for plotPCA(), >getMethod("plotPCA","DESeqTransform") to see source code

# --<>--<>--<>--<>-- PCA to examine effect of treatment --<>--<>--<>--<>--<>--<>--
## uses n most variable genes -- not necessarily sig DEGs!
# plot.treat.data.MON.NAP     <- plotPCA(rld.MON.NAP    , intgroup = c("site"), returnData=TRUE, n=1500)
plot.treat.data.MON.POL     <- plotPCA(rld.MON.POL    , intgroup = c("site"), returnData=TRUE, n=1500)
# plot.treat.data.NAP.POL     <- plotPCA(rld.NAP.POL    , intgroup = c("site"), returnData=TRUE, n=1500)
# plot.treat.data.NAPxMON.MON <- plotPCA(rld.NAPxMON.MON, intgroup = c("site"), returnData=TRUE, n=1500)
# plot.treat.data.NAPxMON.NAP <- plotPCA(rld.NAPxMON.NAP, intgroup = c("site"), returnData=TRUE, n=1500)

my.colors <- c("dodgerblue", "orange")

# percentVar.MON.NAP     <- round(100*attr(plot.treat.data.MON.NAP, "percentVar"))
percentVar.MON.POL     <- round(100*attr(plot.treat.data.MON.POL, "percentVar"))
# percentVar.NAP.POL     <- round(100*attr(plot.treat.data.NAP.POL, "percentVar"))
# percentVar.NAPxMON.MON <- round(100*attr(plot.treat.data.NAPxMON.MON, "percentVar"))
# percentVar.NAPxMON.NAP <- round(100*attr(plot.treat.data.NAPxMON.NAP, "percentVar"))

# ggplot(plot.treat.data.MON.NAP, aes(PC1, PC2, color=site)) +
#   scale_color_manual(values=c(my.colors)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar.MON.NAP[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar.MON.NAP[2],"% variance")) +
#   coord_fixed()

 ggplot(plot.treat.data.MON.POL, aes(PC1, PC2, color=site)) +
   scale_color_manual(values=c(my.colors)) +
   geom_point(size=3) +
   xlab(paste0("PC1: ",percentVar.MON.POL[1],"% variance")) +
   ylab(paste0("PC2: ",percentVar.MON.POL[2],"% variance")) +
   coord_fixed()
# 
# ggplot(plot.treat.data.NAP.POL, aes(PC1, PC2, color=site)) +
#   scale_color_manual(values=c(my.colors)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar.NAP.POL[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar.NAP.POL[2],"% variance")) +
#   coord_fixed()
# 
# ggplot(plot.treat.data.NAPxMON.MON, aes(PC1, PC2, color=site)) +
#   scale_color_manual(values=c(my.colors)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar.NAPxMON.MON[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar.NAPxMON.MON[2],"% variance")) +
#   coord_fixed()
# 
# ggplot(plot.treat.data.NAPxMON.NAP, aes(PC1, PC2, color=site)) +
#   scale_color_manual(values=c(my.colors)) +
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar.NAPxMON.NAP[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar.NAPxMON.NAP[2],"% variance")) +
#   coord_fixed()
#-----------------------------------------------------------------------------------------------------------------------#

################### ID OF DEGs ############################################
#--<>--<>--<>--<>-- How many and which genes are sig DE? --<>--<>--<>--<>--
       
#sum(res.MON.NAP$padj < 0.05, na.rm = T) 
#sum(res.MON.NAP$padj < 0.05 & abs(res.MON.NAP$log2FoldChange) >= 1, na.rm = T)    
sum(res.MON.POL$padj < 0.05, na.rm = T)
sum(res.MON.POL$padj < 0.05 & abs(res.MON.POL$log2FoldChange) >= 1, na.rm = T) #17226 
# sum(res.NAP.POL$padj < 0.05, na.rm = T) 
# sum(res.NAP.POL$padj < 0.05 & abs(res$log2FoldChange) >= 0, na.rm = T)    
# sum(res.NAPxMON.MON$padj < 0.05, na.rm = T)  
# sum(res.NAPxMON.MON$padj < 0.05 & abs(res$log2FoldChange) >= 0, na.rm = T)  
# sum(res.NAPxMON.NAP$padj < 0.05, na.rm = T)  
# sum(res.NAPxMON.NAP$padj < 0.05 & abs(res$log2FoldChange) >= 0, na.rm = T)

# resSig       <- res[which(res$padj < 0.05),]   ## put all genes w/ padj < 0.05 in resSig

# resSig.MON.NAP     <- res.MON.NAP[which(res.MON.NAP$padj < 0.05), ]    
resSig.MON.POL     <- res.MON.POL[which(res.MON.POL$padj < 0.05), ]    
# resSig.NAP.POL     <- res.NAP.POL[which(res.NAP.POL$padj < 0.05), ]    
# resSig.NAPxMON.MON <- res.NAPxMON.MON[which(res.NAPxMON.MON$padj < 0.05), ]
# resSig.NAPxMON.NAP <- res.NAPxMON.NAP[which(res.NAPxMON.NAP$padj < 0.05), ]

## examine MA plot
# plotMA(resSig.MON.NAP, ylim=c(-6,6), xlim=c(1e-2,1e6))  
plotMA(resSig.MON.POL, ylim=c(-6,6), xlim=c(1e-2,1e6))
# plotMA(resSig.NAP.POL, ylim=c(-6,6), xlim=c(1e-2,1e6))
# plotMA(resSig.NAPxMON.MON, ylim=c(-6,6), xlim=c(1e-2,1e6))
# plotMA(resSig.NAPxMON.NAP, ylim=c(-6,6), xlim=c(1e-2,1e6))

## sort by log2foldchange and show most down-regulated
# head(resSig.MON.NAP[order(resSig.MON.NAP$log2FoldChange),])    
head(resSig.MON.POL[order(resSig.MON.POL$log2FoldChange),]) 
# head(resSig.NAP.POL[order(resSig.NAP.POL$log2FoldChange),]) 
# head(resSig.NAPxMON.MON[order(resSig.NAPxMON.MON$log2FoldChange),]) 
# head(resSig.NAPxMON.NAP[order(resSig.NAPxMON.NAP$log2FoldChange),]) 

## sort by log2foldchange and show most up-regulated
tail(resSig.MON.NAP[order(resSig.MON.NAP$log2FoldChange),])    
# tail(resSig.MON.POL[order(resSig.MON.POL$log2FoldChange),])
# tail(resSig.NAP.POL[order(resSig.NAP.POL$log2FoldChange),])
# tail(resSig.NAPxMON.MON[order(resSig.NAPxMON.MON$log2FoldChange),])
# tail(resSig.NAPxMON.NAP[order(resSig.NAPxMON.NAP$log2FoldChange),])

## get trinity contig name for significant genes 
sigGene_id.MON.NAP <- match(rownames(resSig.MON.NAP), contig_names.MON.NAP$gene.id) 
sigGene.MON.NAP <- contig_names.MON.NAP[sigGene_id.MON.NAP,]
resSig.MON.NAP$seq_name <- sigGene.MON.NAP$seq.name ## add to resSig df 

# sigGene_id.MON.POL <- match(rownames(resSig.MON.POL), contig_names.MON.POL$gene.id) 
# sigGene.MON.POL <- contig_names.MON.POL[sigGene_id.MON.POL,]
# resSig.MON.POL$seq_name <- sigGene.MON.POL$seq.name ## add to resSig df 
# 
# sigGene_id.NAP.POL <- match(rownames(resSig.NAP.POL), contig_names.NAP.POL$gene.id) 
# sigGene.NAP.POL <- contig_names.NAP.POL[sigGene_id.NAP.POL,]
# resSig.NAP.POL$seq_name <- sigGene.NAP.POL$seq.name ## add to resSig df 
# 
# sigGene_id.NAPxMON.MON <- match(rownames(resSig.NAPxMON.MON), contig_names.NAPxMON.MON$gene.id) 
# sigGene.NAPxMON.MON <- contig_names.NAPxMON.MON[sigGene_id.NAPxMON.MON,]
# resSig.NAPxMON.MON$seq_name <- sigGene.NAPxMON.MON$seq.name ## add to resSig df 
# 
# sigGene_id.NAPxMON.NAP <- match(rownames(resSig.NAPxMON.NAP), contig_names.NAPxMON.NAP$gene.id) 
# sigGene.NAPxMON.NAP <- contig_names.NAPxMON.NAP[sigGene_id.NAPxMON.NAP,]
# resSig.NAPxMON.NAP$seq_name <- sigGene.NAPxMON.NAP$seq.name ## add to resSig df 

## write CSV with all DEGs, padj < 0.05 
write.csv(as.data.frame(resSig.MON.NAP), file="degs_padj05.MON.NAP.csv") 
# write.csv(as.data.frame(resSig.MON.POL), file="degs_padj05.MON.POL.csv")
# write.csv(as.data.frame(resSig.NAP.POL), file="degs_padj05.NAP.POL.csv")
# write.csv(as.data.frame(resSig.NAPxMON.MON), file="degs_padj05.NAPxMON.MON.csv")
# write.csv(as.data.frame(resSig.NAPxMON.NAP), file="degs_padj05.NAPxMON.NAP.csv")

####### get readnames for piccard FilterSamReads 
readnames <- resSig.MON.NAP[,7]
readnames_MSTRG <- rownames(resSig.MON.NAP) 
write.csv(readnames, file="deg_contignames_noquote.txt", row.names = FALSE, quote=FALSE) ## write read_names.txt for piccard FilterSamReads to subset DEGs reads
#write.csv(readnames_MSTRG, file="read_names_MSTRG_noquote.txt", row.names = FALSE, quote=FALSE) ## write read_names.txt for piccard FilterSamReads to subset DEGs reads

#-----------------------------------------------------------------------------------------------------------------------#

################### VISUALIZING DEGs #######################################

#-----------------------------------------------------------------------------------------------------------------------#

#--<>--<>--<>--<>-- PCA plots - DEGs by cut-offs --<>--<>--<>--<>--
p.cutoff <- 0.05  # fdr p-value
fc.cutoff <- 1    # log fold change 

# topSigGenes.MON.NAP <- res.MON.NAP[which(res.MON.NAP$padj < p.cutoff & abs(res.MON.NAP$log2FoldChange) >= fc.cutoff),]
# num.genes.MON.NAP   <- topSigGenes.MON.NAP@nrows
# rld.Sig.MON.NAP     <- rld.MON.NAP[rownames(rld.MON.NAP) %in% rownames(topSigGenes.MON.NAP)]

topSigGenes.MON.POL <- res.MON.POL[which(res.MON.POL$padj < p.cutoff & abs(res.MON.POL$log2FoldChange) >= fc.cutoff),]
num.genes.MON.POL   <- topSigGenes.MON.POL@nrows
rld.Sig.MON.POL     <- rld.MON.POL[rownames(rld.MON.POL) %in% rownames(topSigGenes.MON.POL)]
# 
# topSigGenes.NAP.POL <- res.NAP.POL[which(res.NAP.POL$padj < p.cutoff & abs(res.NAP.POL$log2FoldChange) >= fc.cutoff),]
# num.genes.NAP.POL   <- topSigGenes.NAP.POL@nrows
# rld.Sig.NAP.POL     <- rld.NAP.POL[rownames(rld.NAP.POL) %in% rownames(topSigGenes.NAP.POL)]

# topSigGenes.NAPxMON.MON <- res.NAPxMON.MON[which(res.NAPxMON.MON$padj < p.cutoff & abs(res.NAPxMON.MON$log2FoldChange) >= fc.cutoff),]
# num.genes.NAPxMON.MON   <- topSigGenes.NAPxMON.MON@nrows
# rld.Sig.NAPxMON.MON     <- rld.NAPxMON.MON[rownames(rld.NAPxMON.MON) %in% rownames(topSigGenes.NAPxMON.MON)]

# topSigGenes.NAPxMON.NAP <- res.NAPxMON.NAP[which(res.NAPxMON.NAP$padj < p.cutoff & abs(res.NAPxMON.NAP$log2FoldChange) >= fc.cutoff),]
# num.genes.NAPxMON.NAP   <- topSigGenes.NAPxMON.NAP@nrows
# rld.Sig.NAPxMON.NAP     <- rld.NAPxMON.NAP[rownames(rld.NAPxMON.NAP) %in% rownames(topSigGenes.NAPxMON.NAP)]



## color PCA by treatment

# pdf(file="~/KW/4_DESeq2/annotated_rnaspades/plots/AllPlots.pdf",
#     width = 8.5,
#     height = 11)

# plot.all.data.MON.NAP <- plotPCA(rld.Sig.MON.NAP, intgroup = c("family", "site"), returnData=TRUE, ntop=num.genes.MON.NAP)
# percentVar.MON.NAP <- round(100*attr(plot.all.data.MON.NAP, "percentVar"))
# p <- ggplot(plot.all.data.MON.NAP, aes(PC1, PC2, color=site, shape=site)) +
#   ggtitle(paste("Top",num.genes.MON.NAP,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +
#   scale_color_manual(values=c("#972C80FF","#FEB77EFF")) + 
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar.MON.NAP[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar.MON.NAP[2],"% variance")) +
#   coord_fixed() +
#   theme_bw() +
#   theme(legend.position = "none") 
# 


plot.all.data.MON.POL <- plotPCA(rld.Sig.MON.POL, intgroup = c("family", "site"), returnData=TRUE, ntop=num.genes.MON.POL)
percentVar.MON.POL <- round(100*attr(plot.all.data.MON.POL, "percentVar"))
ggplot(plot.all.data.MON.POL, aes(PC1, PC2, color=site, shape=site)) +
  ggtitle(paste("Top",num.genes.MON.POL,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +
  scale_color_manual(values=c("firebrick4","darkgoldenrod3")) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar.MON.POL[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.MON.POL[2],"% variance")) +
  coord_fixed()

# plot.all.data.NAP.POL <- plotPCA(rld.Sig.NAP.POL, intgroup = c("family", "site"), returnData=TRUE, ntop=num.genes.NAP.POL)
# percentVar.NAP.POL <- round(100*attr(plot.all.data.NAP.POL, "percentVar"))
# my.colors <- beyonce_palette(66,9,type=c("continuous"))
# ggplot(plot.all.data.NAP.POL, aes(PC1, PC2, color=site, shape=site)) +
#   ggtitle(paste("Top",num.genes.NAP.POL,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +   
#   scale_color_manual(values=c("firebrick4","darkgoldenrod3")) + 
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar.NAP.POL[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar.NAP.POL[2],"% variance")) +
#   coord_fixed()
# 
# plot.all.data.NAPxMON.MON <- plotPCA(rld.Sig.NAPxMON.MON, intgroup = c("family", "site"), returnData=TRUE, ntop=num.genes.NAPxMON.MON)
# percentVar.NAPxMON.MON <- round(100*attr(plot.all.data.NAPxMON.MON, "percentVar"))
# my.colors <- beyonce_palette(66,9,type=c("continuous"))
# ggplot(plot.all.data.NAPxMON.MON, aes(PC1, PC2, color=site, shape=site)) +
#   ggtitle(paste("Top",num.genes.NAPxMON.MON,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +   
#   scale_color_manual(values=c("firebrick4","darkgoldenrod3")) + 
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar.NAPxMON.MON[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar.NAPxMON.MON[2],"% variance")) +
#   coord_fixed()
# 
# plot.all.data.NAPxMON.NAP <- plotPCA(rld.Sig.NAPxMON.NAP, intgroup = c("family", "site"), returnData=TRUE, ntop=num.genes.NAPxMON.NAP)
# percentVar.NAPxMON.NAP <- round(100*attr(plot.all.data.NAPxMON.NAP, "percentVar"))
# my.colors <- beyonce_palette(66,9,type=c("continuous"))
# ggplot(plot.all.data.NAPxMON.NAP, aes(PC1, PC2, color=site, shape=site)) +
#   ggtitle(paste("Top",num.genes.NAPxMON.NAP,"DEGs, padj < ",p.cutoff,"and log2FC >",fc.cutoff)) +   
#   scale_color_manual(values=c("firebrick4","darkgoldenrod3")) + 
#   geom_point(size=3) +
#   xlab(paste0("PC1: ",percentVar.NAPxMON.NAP[1],"% variance")) +
#   ylab(paste0("PC2: ",percentVar.NAPxMON.NAP[2],"% variance")) +
#   coord_fixed()

## plot all genes

plot.all.genes.MON.NAP <- plotPCA(rld.MON.NAP, intgroup = c("family", "site"), returnData=TRUE)
percentVar.allgenes.MON.NAP <- round(100*attr(plot.all.genes.MON.NAP, "percentVar"))

## color all genes PCA by treatment
ggplot(plot.all.genes.MON.NAP, aes(PC1, PC2, color=site, shape=site)) +
  ggtitle(paste("All genes with nonzero total read count")) +   
  scale_color_manual(values=c("#972C80FF","#FEB77EFF")) + 
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar.allgenes.MON.NAP[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.allgenes.MON.NAP[2],"% variance")) + 
  coord_fixed() 

# dev.off()


#-----color duplicates-----------------------------------------------------------------#
n_occur <- data.frame(table(plot.all.data.MON.NAP$family))
n_occur2 <- n_occur[n_occur$Freq>1,]
duplicate.plot <- plot.all.data.MON.NAP[plot.all.data.MON.NAP$family %in% n_occur2$Var1,]
nonduplicate.plot <- plot.all.data.MON.NAP[!plot.all.data.MON.NAP$family %in% n_occur2$Var1,]

ggplot(nonduplicate.plot, aes(PC1, PC2, shape=site))+
  geom_point(size=2)+
  xlab(paste0("PC1: ",percentVar.MON.NAP[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.MON.NAP[2],"% variance")) +
  coord_fixed() +
  geom_line(data=duplicate.plot, aes(group=family), color="grey") + 
  geom_point(data=duplicate.plot, aes(PC1, PC2, shape=site, color=family), size=2) +
  scale_fill_viridis()
  theme_test()
  

##------- OTHER PLOTS ------------##

## construct heatmaps to examine sample patterns with different clustering methods
rownames(sampleDistMatrix) <- paste(rld$family, rld$treatment, sep="-")
colnames(sampleDistMatrix) <- paste(rld$family, rld$treatment, sep="-")
dists.ward <- hclust(sampleDists, method="ward.D2")
heatmap.2(sampleDistMatrix, trace="none", Rowv=as.dendrogram(dists.ward),
          margins=c(6,7))
dists.avg <- hclust(sampleDists, method="average")
heatmap.2(sampleDistMatrix, trace="none", Rowv=as.dendrogram(dists.avg),
          margins=c(6,7))
heatmap.2(sampleDistMatrix, trace="none", dendrogram="column",
          scale="row", margins=c(6,7))

#--<>--<>--<>--<>-- Heatmap w/ top significant DEGs --<>--<>--<>--<>--

topSigGenes <- head(resSig[order(resSig$padj),],n=10)      ## get 10 genes with lowest padj values
colours <- beyonce_palette(64,300,type=c("continuous"))

#!# genes on y, samples on x
heatmap.2(assay(rld)[rownames(topSigGenes),], 
          scale="row", ## plot heatmap using these genes
          trace="none", dendrogram="column",margins=c(4,10), key.title=NA,
          cexRow=1, cexCol=1,
          col = colours)
