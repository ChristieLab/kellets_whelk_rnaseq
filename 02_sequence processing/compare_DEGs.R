setwd("~/KW/4_DESeq2/annotated_rnaspades/")

mon.pol <- read.csv("degs_padj05.MON.POL.csv")
nap.pol <- read.csv("degs_padj05.NAP.POL.csv")
mon.nap <- read.csv("degs_padj05_rnaspades_annotated.csv")
napxmon.mon <- read.csv("degs_padj05.NAPxMON.MON.csv")
napxmon.nap <- read.csv("degs_padj05.NAPxMON.NAP.csv")


nrow(mon.nap)
nrow(napxmon.mon)
nrow(napxmon.nap)

## try collapsing POL isoforms 
mon.pol$gene_names <- gsub("_i[0-9]+$", "", mon.pol$seq_name)
length(unique(mon.pol$seq_name))
length(unique(mon.pol$gene_name))

nap.pol$gene_names <- gsub("_i[0-9]+$", "", nap.pol$seq_name)
length(unique(nap.pol$seq_name))
length(unique(nap.pol$gene_name))

length(intersect(mon.nap$seq_name, napxmon.mon$seq_name)) / 2770
length(intersect(mon.nap$seq_name, napxmon.nap$seq_name)) / 2770

