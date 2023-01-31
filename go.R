# Gene Ontology R 
library(tidyr)

setwd("~/KW/07_geneontology")

dat<-read.csv("uniprot_alldegs_go.tsv", sep = "\t")



dat <- separate_rows(dat, Gene.Ontology.IDs, sep = "; ")

table <- table(dat$Gene.Ontology.IDs)

tail(table)
table <- table[order(table, decreasing = TRUE)]

# get gene names to input into Quick GO to figure out which ones in what category 
table.df <- data.frame(table)
table.reduced <- subset(table.df, Freq>15)
table.reduced$Var1
write.csv(table.reduced$Var1,file = "gonames.csv", quote = FALSE, col.names = NULL)


## 
top_bp_gonames <-  c("GO:0007283", "GO:0007165", "GO:0007018", "GO:0060271", "GO:0045893", "GO:0006357", "GO:0045944", "GO:0030317", "GO:0003341", "GO:0015031", "GO:0015074", "GO:0030154", "GO:0006915", "GO:0002181", "GO:0006468", "GO:0006412", "GO:0006508") 



go_terms <- read.csv("goterms.txt", sep="\t")
top_bp_go <- table.df[table.df$Var1 %in% go_terms$term, ]
names(top_bp_go)=c("term","Freq")

combined <- inner_join(go_terms, top_bp_go)
combined <- combined[order(combined$Freq, decreasing = TRUE),]
