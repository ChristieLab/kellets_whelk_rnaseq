# Top Go 

#BiocManager::install("topGO")

library(Rgraphviz)
library(topGO)
# browseVignettes("topGO")
setwd("/Users/andy/KW/07_geneontology")

############ Comparing DEGs against reference transcriptome 

geneID2GO <- readMappings(file = "eggnog_results/all_ref_genes_to_go.txt")
degID2GO <- readMappings(file = "eggnog_results/all_deg_genes_to_go.txt")

str(head(geneID2GO))

geneNames <- names(geneID2GO)
myInterestingGenes <- names(degID2GO)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
str(geneList)

GOdata <- new("topGOdata", 
              description = "Kellet's Whelks DEGs to Ref",
              ontology = "BP",
              allGenes = geneList, 
              annot = annFUN.gene2GO, 
              gene2GO = geneID2GO)

# working with the topGO object
description(GOdata) # dscription 
genes(GOdata) ## list of genes 
numGenes(GOdata) # number of genes with GO Terms n = 106581


sg <- sigGenes(GOdata) # get significant genes 
str(sg)
numSigGenes(GOdata) # n = 124 in DEG set, n = 1225 compared to ref

graph(GOdata) ## returns the GO graph 
ug <- usedGO(GOdata) ##GOs used in analysis 

num.ann.genes <- countGenesInTerm(GOdata) ## the number of annotated genes
num.ann.genes # get counts of go terms 
ann.genes <- genesInTerm(GOdata) ## get the annotations
head(ann.genes) # genes with each go term 

# get stats for some terms 
# termStat(GOdata, sel.terms) #sel terms is a list of GO terms of interest

# running GO enrichment getSigGroups(), can run on specific go terms or group of go terms 
### get statistics of enrichment tests
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisher.elim <- runTest(GOdata, algorithm = "elim", statistic = "fisher")
# resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks") # doesn't apply to our study, uses expression data to find DEGs
# more flexible way of running enrichment tests, can define own statistics by using getSigGroups(), see manual 


## combine results of different tests 
# allRes <- GenTable(GOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "classicFisher", ranksOf = "elimFisher", topNodes = 10)

allRes <- GenTable(GOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim,
         orderBy = "elimFisher", ranksOf = "classicFisher", topNodes = 1225)

write.csv(allRes, "topGO_deg2ref_results.csv")
go_withPval <- cbind(allRes$GO.ID, allRes$elimFisher)
colnames(go_withPval) <- c("goID","elimFisher p-value")
write.csv(go_withPval, "deg2ref_enrichedGOid.csv")
# write.csv(ug, "allGOused.csv")


## get p-value 
pValue.classic <- score(resultFisher)
pValue.elim <- score(resultFisher.elim)[names(pValue.classic)]
gstat <- termStat(GOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4

# classic vs. elim (elim more conservative)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize)
hist(pValue.classic, 50, xlab ="p-values")
hist(pValue.elim, 50, xlab ="p-values")

# show top 5 nodes GO map 
showSigOfNodes(GOdata, score(resultFisher.elim), firstSigNodes = 1, useInfo = 'all')






############ find significant genes within the DEG set 
deggeneNames <- myInterestingGenes
randomInterestingGenes <- sample(geneNames, length(geneNames) / 10)
randomgeneList <- factor(as.integer(geneNames %in% randomInterestingGenes))
names(randomgeneList) <- deggeneNames
str(geneList)

degonlyGOdata <-new("topGOdata", 
                    description = "Kellet's Whelks DEGs only",
                    ontology = "BP",
                    allGenes = geneList, 
                    annot = annFUN.gene2GO, 
                    gene2GO = geneID2GO)

description(degonlyGOdata) # dscription 
genes(degonlyGOdata) ## list of genes 
numGenes(degonlyGOdata) # number of genes with GO Terms n = 106581


sg <- sigGenes(degonlyGOdata) # get significant genes 
str(sg)
numSigGenes(degonlyGOdata) # n = 124 in DEG set, n = 1225 compared to ref

graph(degonlyGOdata) ## returns the GO graph 
ug <- usedGO(degonlyGOdata) ##GOs used in analysis 

num.ann.genes <- countGenesInTerm(degonlyGOdata) ## the number of annotated genes
num.ann.genes # get counts of go terms 
ann.genes <- genesInTerm(degonlyGOdata) ## get the annotations
head(ann.genes) # genes with each go term 

# get stats for some terms 
# termStat(degonlyGOdata, sel.terms) #sel terms is a list of GO terms of interest

# running GO enrichment getSigGroups(), can run on specific go terms or group of go terms 
### get statistics of enrichment tests
resultFisher <- runTest(degonlyGOdata, algorithm = "classic", statistic = "fisher")
resultFisher.elim <- runTest(degonlyGOdata, algorithm = "elim", statistic = "fisher")
# resultKS <- runTest(degonlyGOdata, algorithm = "classic", statistic = "ks") # doesn't apply to our study, uses expression data to find DEGs
# more flexible way of running enrichment tests, can define own statistics by using getSigGroups(), see manual 


## combine results of different tests 
allRes <- GenTable(degonlyGOdata, classicFisher = resultFisher, elimFisher = resultFisher.elim, orderBy = "classicFisher", ranksOf = "elimFisher", topNodes = 10)


################## explore enriched genes within the DEG set ###################
