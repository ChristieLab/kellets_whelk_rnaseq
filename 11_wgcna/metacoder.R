#===========================================================================#
# Script written by Andy Lee on 4/5/2023
# Script modified from Janna Willoughby and Avril Harder
#
# This script generates metacoder trees for a given list of GO terms;
# See code_notes.pdf for additional info. on functionality
#
# Update March 17, 2021: script was previously only run in R v3 and some 
# required packages may not be compatible with R v4.
#===========================================================================#

##### AMH 4/14/23: this still works for me; I'm running R v4.0.3 #####

library(GO.db)
library(metacoder)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GOfuncR)

# set module of interest
module <- "thistle" # darkred

setwd("/Users/andy/KW/11_wgcna/")

##### READ IN ALL MODULE GENE:GO INFORMATION #####
## get list of all GO IDs and corresponding terms
terms <- Term(GOTERM)
terms <- as.data.frame(terms)
terms <- cbind(rownames(terms), terms)
rownames(terms) <- c()
colnames(terms) <- c("go.id","term")

# ## get list of gene:GO term relationships
# ##
# ## gene.name       go.id
# ##
# rep.mods <- read.csv("/Users/Avril/Documents/rna_seq/seq_data_processing_notes_and_analyses/wgcna/go_for_janna/rep4_module_gene_go_info.csv")
# ## keep info for module of interest
# rep.mods <- rep.mods[which(rep.mods$module==module),]
# rep.mods$go.id <- as.character(rep.mods$go.id)
# rep.mods <- rep.mods[,c(2,4)]


##### GET GO TERMS UNIQUE TO MODULE OF INTEREST #####
##
## ID       term       count
##
uni.terms <- read.csv(paste0(module,"_unique_GO_terms.csv"))


##### GET LIST OF INTERESTING GO TERMS (filt.go) THAT ARE ALSO UNIQUE TO MODULE #####
n <- 30 ## can limit the # of GO terms you want to include, just to see how it runs
final.go <- head(uni.terms, n=n)
# final.go <- uni.terms[which(uni.terms$ID %in% filt.go),]
# rm(list=ls()[! ls() %in% c("final.go","module")])

##### PASS FINAL GO LIST TO METACODER TO PLOT TREES #####
setwd("~/Desktop/")

bpdata <- as.data.frame(final.go$ID)
colnames(bpdata) = c("GO.ID")
bpdata$GO.ID <- as.character(bpdata$GO.ID)

# function needed to parse GO IDs
##### AL 4/1/7/23: the function does not work, writing a new one. 
##### below are test codes 
##### need to get only biological processes? 


get_parent_nodes("GO:0000315")

##### 4/17/2023 AL: Code below works fine

# modify, pull go term relationships
bpterms = lapply(bpdata$GO.ID, term_class, all_paths = FALSE, type = GOBPPARENTS) #this line can take forever
bpres   = data.frame(class=unlist(bpterms))

# write/write data (annoyingly the only way I can get it to work, but at least it works)
write.table(bpres, "data/temp.csv", sep=",", col.names=TRUE, row.names=FALSE)
bpres=read.table("data/temp.csv", header=TRUE, sep=",")

##### AMH 4/14/23: the function above does not work for me. I found an old copy of temp.csv, though, so you can see how the information was formatted. looks like you need to figure out how to get all the ancestor terms for each term you want to display, then collapse that ordered list using '|' as a separator #####

##### AMH 4/14/23: reading in the old file I found to see if everything else still works #####

##### AL 4/17/23: The following code works with metacoder 0.3.6



## bpdata contains the GO IDs for the terminal nodes
## bpres and bpdata contain the same information, different formats = path information for each terminal node (GO Terms)

#parse GO data
##### bpres=read.table("temp.csv", header=TRUE, sep=",") # test with Avril's file 

bpres2 <- data.frame(lapply(bpres, function(x) {
  gsub("\\|", ";", x)
})) # change format for parse_tax_data()

testdata <- parse_tax_data(bpres2, class_sep = ";")

# data = separate(bpres, c("class" = -1) , c("1","2","3","4","5","6","7","8","9","10","11") , sep = '\\|')

##### AMH 4/14/23: big warning about deprecated command here; looks like the old version of metacoder uses a deprecated dplyr command inside the parse_taxonomy_table command -- maybe why my output figure looks so shitty #####

# parse GO data (parse_taxonomy_table replaced by taxa::parse_tax_data)
# data <- parse_tax_data(bpres,
#                        class_sep = "\\|")


#create figure
tempdata = filter_taxa(testdata, n_supertaxa <= 500) #filters terms that are WAAAAY out from the middle
## n_supertaxa sets # of nodes that can be in a single path from center to terminal nodes (cuts from terminal end, not internal nodes)
##### AMH 4/14/23: had to update lazyeval package here #####

pdf(paste("",module,"meta_labels.pdf", sep=""), width=5, height=5, useDingbats=FALSE, onefile=FALSE)
#grey70      grey40      dodgerblue2 tomato2  chartreuse3  darkorchid2 yellow3 orange2 magenta1
par(bg=NA)
set.seed(100) ## seeds: orange = 35       lightsteelblue1 = 75       paleturquoise = 100
heat_tree(tempdata, 
          node_label = taxon_names, ### modified by AL 
          # node_size = colsandsize$ngenes,
          # node_size_trans = "log10",
          node_size_range = c(0.01, 0.01),
          # node_label_size_trans = "log10",
          node_label_size_range = c(0.01, 0.01),
          # edge_size_trans = "log10",
          edge_size_range = c(0.004, 0.004),
          node_color = module,
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
dev.off()



## run this version to find best seed #s for each module's tree - have to highlight and run manually :/
# par(bg=NA)
# sede <- sample(1:100,1)
# pdf(paste0("/Users/Avril/Desktop/",module,"/",sede,".pdf"), width=5, height=5, useDingbats=FALSE, onefile=F)
# set.seed(sede) #grey70      grey40      dodgerblue2 tomato2  chartreuse3  darkorchid2 yellow3 orange2 magenta1
# heat_tree(tempdata, #node_label = tempdata$taxon_data$name,
#           # node_size = colsandsize$ngenes,
#           # node_size_trans = "log10",
#           node_size_range = c(0.01, 0.01),
#           # node_label_size_trans = "log10",
#           node_label_size_range = c(0.01, 0.01),
#           # edge_size_trans = "log10",
#           edge_size_range = c(0.004, 0.004),
#           node_color = module,
#           # node_color_trans = "linear",
#           # node_color_range = diverging_palette(),
#           # node_color_interval = c(-4, 4),
#           # edge_color_trans = "linear",
#           # edge_color_range = diverging_palette(),
#           # edge_color_interval =  c(-4, 4),
#           node_label_max = 500,
#           # node_color_axis_label = "Factor change",
#           # node_size_axis_label = "Number of genes",
#           layout = "da", initial_layout = "re"
# )
# dev.off()

