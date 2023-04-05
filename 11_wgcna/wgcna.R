### Weighted gene co-expression network analysis (WGCNA) for RNA-seq data ###


### written by Andy Lee
### modified from Avril Harder and Mark Christie 
### Largely based on WGCNA tutorials at: 
### https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/ 

setwd("~/KW/11_wgcna")
#BiocManager::install("WGCNA")
#BiocManager::install("impute")
#BiocManager::install("preprocessCore")

library("DESeq2")
library("WGCNA")
library("dplyr")
library(GO.db)

##### DATA LOADING AND PREP #####
# load("through_logxform.RData") ## saved during DESeq2 analysis

#### DESeq 2 
## all.gene.counts <- read.table("~/KW/4_DESeq2/annotated_rnaspades/featureCount_counts## .txt", sep="\t", header=T, row.names=1, check.names=F)
## drop.vars   <- names(all.gene.counts) %in% c("Chr","Start", "End", "Strand", ## "Length")
## all.gene.counts <- all.gene.counts[!drop.vars]
## samples <- read.csv("~/KW/4_DESeq2/annotated_rnaspades/kw_sample_info.csv")
## 
## m1 <- which(samples[, 3] == "POL") # might want to combine POL with NAP later
## m2 <- which(samples[, 3] == "NAPxMON")
## m3 <- which(samples[, 3] == "MON")
## m4 <- which(samples[ ,3] == "NAP")
## 
## crosses <- c(m1, m2)
## 
## # should all be sites we want to exclude (for now)!
## samples[-crosses, ]
## # with sample 22 excluded (outlier) and all HxW and WxH offspring excluded
## gene.counts <- all.gene.counts[, -crosses]
## gene.counts
## # read in CSV with sample information
## samples <- read.csv("~/KW/4_DESeq2/annotated_rnaspades/kw_sample_info.csv", row## .names =1)
## col.data <- samples[-crosses, ]
## 
## 
## # check to make sure that sample names are in the same order in the gene count table ## and the sample info table; must be TRUE or do not proceed
## all(rownames(col.data) %in% colnames(gene.counts))
## all(rownames(col.data) == colnames(gene.counts))
## 
## dds <- DESeqDataSetFromMatrix(countData = gene.counts,
##                               colData = col.data,
##                               design = ~ as.factor(site))
## dds$site <- relevel(dds$site, "NAP") ## set the untreated group as the control
## #$site <- relevel(dds$site, "MON")
## dds

traitData <- readRDS("mon_nap_trait.rds")
dds <- readRDS("dds.rds")

#### load data for WCGNA 
new.dds <- estimateSizeFactors(dds)
counts <- counts(new.dds, normalized=TRUE)

options(stringsAsFactors = FALSE) # iimportant setting, don't delete 
dim(counts) ## 54 samples,  167,051 transcripts
rownames(counts) ## sample names

## check for genes and samples with too many missing values
gsg <- goodSamplesGenes(counts, verbose = 3)
gsg$allOK ## <-- this should return "TRUE" ==> if it doesn't, run below loop

# if (!gsg$allOK)
# {
#   ## this is where genes that weren't "good" were filtered out.
#   ## optionally, print the gene and sample names that were removed:
#   if (sum(!gsg$goodGenes)>0)
#     printFlush(paste("Removing genes:", paste(names(counts)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0)
#     printFlush(paste("Removing samples:", paste(rownames(counts)[!gsg$goodSamples], collapse = ", ")));
#   # Remove the offending genes and samples from the data:
#   counts = counts[gsg$goodSamples, gsg$goodGenes]
# }

## how many genes have at least 30 samples with non-zero counts? all 
dim(counts[rowSums(counts==0) <= 30,])

## make this the new criteria for counts to speed things up for outlier check step
counts <- counts[rowSums(counts==0) <= 30,]

# ## cluster samples to see if there are any obvious outliers
# sampleTree = hclust(dist(counts), method = "average")
# par(cex = 0.6)
# par(mar = c(2,4.5,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)
# 
# ## Static tree cutting
# # Plot a line to show the cut
# abline(h = 80000, col = "red", lty=3);
# # Determine cluster under the line
# clust = cutreeStatic(sampleTree, cutHeight = 80000, minSize = 10)
# table(clust) ## how many will be cut [0] and how many will be kept [1]
# # clust 1 contains the 35 samples we want to keep.
# keepSamples = (clust==1)
# # datExpr = counts[keepSamples, ]
datExpr <- t(counts)
# sampleTree = hclust(dist(datExpr), method = "average")
# par(cex = 0.6)
# par(mar = c(2,4.5,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 1.5, cex.main = 2)

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
## gene expression data now ready for expression analysis

## read in trait data and match info. to expression samples
## stack and treatment ID recoded as #s

# traitData <- col.data
dim(traitData)
names(traitData)

## make data frame to match expression data
sampnames <- rownames(datExpr)
traitrows <- match(sampnames, rownames(traitData))
datTraits <- traitData[traitrows, 2, drop=FALSE]

datTraits$site <- as.numeric(as.factor(datTraits$site))

# rownames(datTraits) <- rownames(traitData[traitrows, ])
collectGarbage();

## important data are now in datExpr and datTraits
## want to see how the traits relate to the sample dendrogram
sampleTree2 <- hclust(dist(datExpr), method="average")

## convert traits to a colr representation; white = low values, red = high value, grey = missing value
traitColors <- numbers2colors(datTraits, signed=FALSE)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels=names(datTraits),
                    main="Sample dendrogram and trait heatmap")

save(datExpr, datTraits, file = "kw-WGCNA-dataInput.RData")

##### Network Construction #####
# choose a set of soft-thresholding powers to start with
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
## call the network topology analysis function
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
## plot the results
par(mfrow = c(1,2))
cex1 = 0.9
## scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

## this line corresponds to using an R^2 cut-off of h
abline(h=0.975,col="darkred")
abline(h=0.90,col="red")

## mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
## going with 6 

## construct gene network and ID modules
## power determined above
## TOM = topological overlap matrix
## blockwiseModules() == automatic module detection via dynamic tree cutting

?blockwiseModules()

net <- blockwiseModules(datExpr, power = 6,
                        TOMType = "signed", minModuleSize = 30,
                        reassignThreshold = 0, mergeCutHeight = 0.10, # manual calls for 0.25, Avril uses 0.1
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = TRUE,
                        networkType = "signed",
                        saveTOMFileBase = "kw_blockwise_TOM", 
                        verbose = 3)

saveRDS(net, "wgcna_net.rds")
net <- readRDS("wgcna_net.rds")

## check to see how many modules were ID'ed and the size of the modules
## IF THIS ERRORS OUT: run search(); ".GlobalEnv" and "package:WGCNA" should be the first 2 in the list
## if they're not, restart R and only load WGCNA

# save.image("rep4_allRData.RData")

##### STARTING POINT FOR AUTO MODULE ID W/ POWER=6, MINMODSIZE=30 , MERGECUTHEIGHT=0.10 #####

# load("rep4_allRData.RData")

table(net$colors) ## corresponds to modules and numbers of genes in each module. genes in module "0" weren't assigned to any module


## convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)
## plot the dendrogram and the module colors underneath


plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

pdf(file="Cluster Dendrograms.pdf", 
    width = 5,
    height = 5 )

for (i in 1:33) {
  plotDendroAndColors(net$dendrograms[[i]], mergedColors[net$blockGenes[[i]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}

dev.off()

## save the module assignment and module eigengene information necessary for subsequent analysis
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs;
geneTree <- net$dendrograms[[1]];

##### Relate modules to external info and ID important genes #####

## define numbers of genes and samples
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)

## recalculate MEs with color labels
MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
# var.explain <- moduleEigengenes(datExpr, moduleColors)$varExplained
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, datTraits, use = "p")  #Mark added this line
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)


# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "  (", signif(moduleTraitPvalue, 1), ")", sep = ""); #modified by andy
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               #  =ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

##### GETTING LIST OF GENE IDS #####
sig.modules <- moduleTraitPvalue[apply(moduleTraitPvalue, 1, function(row) {any(row < 0.05)}),]

list.sig.modules <- names(sig.modules)
for (i in 1:length(list.sig.modules)) {
  list.sig.modules[i] <- substring(list.sig.modules[i], 3) ## remove "ME" from the beginnings of all module names
}

module.gene.IDs <- data.frame()
i <- 1
for (module in list.sig.modules) {
  tmp.gene.IDs <- data.frame()
  mod.genes <- dimnames(datExpr)[[2]][moduleColors==module]
  tmp.gene.IDs[i:(length(mod.genes)),1] <- module
  tmp.gene.IDs[i:(length(mod.genes)),2] <- mod.genes
  module.gene.IDs <- rbind(module.gene.IDs, tmp.gene.IDs)
}
module.gene.IDs
sum(table(module.gene.IDs))

#### added by Andy 
#### create gene_names, which is a table linking gene name, transcript ID, and xloc information for later distinguishing different isoforms of the same gene (same gene start and end) and transcripts generated from different paralogs of the same gene (by xloc)
assembly <- rtracklayer::readGFF("~/KW/4_DESeq2/annotated_rnaspades/stringtie_all_merged_kw.gtf")             ## read in merged GTF
gene_idx <- match(dds@rowRanges@partitioning@NAMES, assembly$gene_id)
seq.name      <- as.character(assembly$seqid[gene_idx])
gene.id <- assembly$gene_id[gene_idx]
# xloc           <- assembly$exon_number[gene_idx]
# transcript.id  <- assembly$transcript_id[gene_idx]
contig_names     <- as.data.frame(cbind(seq.name, gene.id))

sigGene_id <- match(module.gene.IDs$V2, contig_names$gene.id) ## get trinity contig name for significant genes 
sigGene <- contig_names[sigGene_id,]
module.gene.IDs$contigname <- sigGene$seq.name

write.csv(module.gene.IDs, "all_module_gene_IDs.csv")


##### BUILD LIST OF MODULES SIG. ASSOC. WITH VARS, GET CORS AND P-VALUES #####
head(moduleTraitCor)
head(moduleTraitPvalue)

all(rownames(moduleTraitCor) == rownames(moduleTraitPvalue))
all(colnames(moduleTraitCor) == colnames(moduleTraitPvalue))

p.vals <- data.frame(module=character(),
                     n.genes=numeric(),
                     pred.var=character(),
                     cor=numeric(),
                     p.value=numeric(),
                     stringsAsFactors=FALSE)
temp <- data.frame(module=character(),
                   n.genes=numeric(),
                   pred.var=character(),
                   cor=numeric(),
                   p.value=numeric(),
                   stringsAsFactors=FALSE)

for(row in 1:nrow(moduleTraitPvalue)) {
  for (col in 1:ncol(moduleTraitPvalue)) {
    if (moduleTraitPvalue[row,col]<0.05) {
      temp.mod <- substring(rownames(moduleTraitPvalue)[row],3)
      temp[1,1] <- rownames(moduleTraitPvalue)[row]
      temp[1,3] <- colnames(moduleTraitPvalue)[col]
      temp[1,4] <- moduleTraitCor[row,col]
      temp[1,5] <- moduleTraitPvalue[row,col]
      temp[1,2] <- length(colnames(datExpr)[moduleColors==temp.mod])
      p.vals <- rbind(p.vals, temp)
    }
  }
}

p.vals
plot(abs(p.vals$cor), p.vals$p.value)  # relationship between correlation and p-value

sequential.bonferroni   <- p.adjust(p.vals[, 5], method = "holm")
pvals.out   <- cbind(p.vals, sequential.bonferroni)
# write.csv(pvals.out,"significant_modules.csv", row.names=F)


# pvals <- sequential.bonferroni  # did not use Bonferroni here
pvals <- pvals.out$p.value
top.modules <- p.vals[which(pvals <= 0.05),]
rownames(top.modules) <- top.modules[,1]
top.modules <- top.modules[,-1]
top.modules[order(top.modules$p.value), ]

top.modules.seqbon <- p.vals[which(sequential.bonferroni <= 0.05),]
rownames(top.modules.seqbon) <- top.modules.seqbon[,1]
top.modules.seqbon <- top.modules.seqbon[,-1]
top.modules.seqbon[order(top.modules.seqbon$p.value), ] 

# read avril's paper 
# table S4, MIGHT HAVE TO DO SOMETHING ABOUT SHARED GENES 


####gene module membership and gene module significance============================================#
  
m1 <- match(p.vals[, 1], names(MEs))

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, datTraits, use = "p")); # modified 
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(datTraits), sep="");
names(GSPvalue) = paste("p.GS.", names(datTraits), sep="");

OUT <- NULL
# top.names <- substring(rownames(top.modules), 3)
top.names <- substring(rownames(top.modules.seqbon), 3)
for(n in 1:length(top.names)){
  module = top.names[n]
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  
  #sizeGrWindow(7, 7);
  #par(mfrow = c(1,1));
  #verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), abs(geneTraitSignificance[moduleGenes, 1]),
  #xlab = paste("Module Membership in", module, "module"),
  #ylab = "Gene significance for site",
  #main = paste("Module membership vs. gene significance\n"),
  #cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
  gene.module.membership  <-  abs(geneModuleMembership[moduleGenes, column])  # abs implies direction of gene specific association removed here
  gene.trait.significance <-  abs(geneTraitSignificance[moduleGenes, 1])
  tp <- GSPvalue[moduleGenes, 1]
  gene.id <- rownames(GSPvalue)[moduleGenes]
  sigGene_id <- match(gene.id, contig_names$gene.id) ## get trinity contig name for significant genes 
  sigGene <- contig_names[sigGene_id,]
  genename <- sigGene$seq.name # prev 3 lines added by andy
  output <- cbind(module, top.modules[n, -2], gene.id, genename, gene.module.membership, gene.trait.significance, tp)
  OUT <- rbind(OUT, output)
}

OUT #warning message here is okay

saveRDS(OUT, "wgcna_output.rds")
write.csv(OUT, "wgcna_output.csv")

#OUT <- readRDS("wgcna_output.rds")
unique(OUT$module)





###### match genename with GO Terms from eggNOG
refGenes2Go <- read.table("/Users/andy/KW/07_geneontology/eggnog_results/all_ref_genes_to_go.txt") 
refGenes2Go$V1 <- substr(refGenes2Go$V1, 1, nchar(refGenes2Go$V1)-2)

OUT$genename %in% refGenes2Go$V1 # not all genes ID'd by WGCNA have GO's associated with them. Some of them are not even included in the eggNOG output 


WGCNA_result_with_GO <- OUT %>% left_join(refGenes2Go, by = join_by(genename == V1), relationship = "many-to-many")
colnames(WGCNA_result_with_GO)[10] <- "ID"
WGCNA_keep_GO <- filter(WGCNA_result_with_GO, ID != "NA" & ID != "-")
  
thistle <- WGCNA_keep_GO[WGCNA_keep_GO$module=="thistle", ]
thistle.go <- strsplit(thistle$ID, split = ",")
thistle.go <- unlist(thistle.go)
thistle.go.count <- as.data.frame(table(thistle.go))

red <- WGCNA_keep_GO[WGCNA_keep_GO$module=="red", ]
red.go <- strsplit(red$ID, split = ",")
red.go <- unlist(red.go)
red.go.count <- as.data.frame(table(red.go))

goterms <- Term(GOTERM)
goterms <- strsplit(goterms, split = "\t")
goterms <- as.data.frame(as.matrix(goterms))
goterms$ID <- rownames(goterms)

thistle.go.count <- thistle.go.count %>% left_join(goterms, by = join_by(thistle.go == ID))
thistle.go.count <- thistle.go.count[,c(1,3,2)]
colnames(thistle.go.count) <- c("ID", "Term", "Freq")
thistle.go.count <- apply(thistle.go.count, 2, as.character)

red.go.count <- red.go.count %>% left_join(goterms, by = join_by(red.go == ID))
red.go.count <- red.go.count[,c(1,3,2)]
colnames(red.go.count) <- c("ID", "Term", "Freq")
red.go.count <- apply(red.go.count, 2, as.character)

write.csv(thistle.go.count, "thistle_unique_GO_terms.csv")
write.csv(red.go.count, "red_unique_GO_terms.csv")
