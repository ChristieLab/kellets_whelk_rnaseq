#===========================================================================#
# Script written by Andy Lee on 4/5/2023
# Script modified from Janna Willoughby and Avril Harder
#
# This script generates metacoder trees for a given list of GO terms;
# See code_notes.pdf for additional info. on functionality
#
#===========================================================================#


### some notes about working with GO.db 
### https://biology.stackexchange.com/questions/44887/get-parent-go-terms-of-go-term-vector



# BiocManager::install("metacoder")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("AnnotationDbi")
# BiocManager::install("gtable")
library(GO.db)
library(metacoder)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)

# set module of interest
setwd("/Users/andy/KW/11_wgcna/enriched_go_bp/")
setwd("/Users/andy/KW/11_wgcna/enriched_go_MF/")

##### GET GO TERMS UNIQUE TO MODULE OF INTEREST #####
module <-"red" # significant after correction 
module <-"thistle" # significant after correction 
module <-"darkorange"
module <-"floralwhite"
module <-"darkviolet"
module <-"palevioletred3"
module <-"brown2"
module <-"lightcyan"
module <-"yellowgreen"
module <- "darkgreen"
module <- "green"
module <- "grey" #errored out for some reason for MF
module <- "darkorange2"
module <- "purple"
module <- "darkmagenta"
module <- "cyan"
module <- "plum1"
module <- "grey60"
module <- "darkslateblue"
module <- "darkgrey"
module <- "darkseagreen4"

uni.terms <- read.csv(paste0("/Users/andy/KW/11_wgcna/enriched_go_bp/", module, "_module_top20_enrichedterms.csv")) # bp

uni.terms <- uni.terms[1:19,]

# uni.terms <- read.csv(paste0("/Users/andy/KW/11_wgcna/enriched_go_MF/", module, "_module_top20_MF_enrichedterms.csv")) 
#uni.terms <- read.csv(paste0("/Users/andy/KW/07_geneontology/topGO_deg2ref_results.csv")) # bp
#uni.terms <- uni.terms[1:20,]
##### GET LIST OF INTERESTING GO TERMS (filt.go) THAT ARE ALSO UNIQUE TO MODULE #####

terms <- AnnotationDbi::select(GO.db, uni.terms$GO.ID, columns = c('TERM','ONTOLOGY'), keytype = "GOID") 
bp.terms <- terms[terms$ONTOLOGY=="BP", ] # grab only BP terms 
bp.terms <- bp.terms[!is.na(bp.terms), ]

#mf.terms <- terms[terms$ONTOLOGY=="MF", ] # grab only mf terms 
#mf.terms <- mf.terms[!is.na(mf.terms), ]

# n <- 30 ## can limit the # of GO terms you want to include, just to see how it runs
# final.go <- head(bp.terms, n=n)

# final.go <- uni.terms[which(uni.terms$ID %in% filt.go),]
# rm(list=ls()[! ls() %in% c("final.go","module")])

##### PASS FINAL GO LIST TO METACODER TO PLOT TREES #####
# setwd("~/Desktop/")


# function needed to parse GO IDs
##### AL 4/17/23: the function does not work, writing a new one. 

#function needed to parse GO IDs

gobpparents <- as.list(GOBPPARENTS)
gobpparents <- gobpparents[!is.na(gobpparents)] # Remove GO IDs that do not have any parent
# gomfparents <- as.list(GOMFPARENTS)
# gomfparents <- gomfparents[!is.na(gomfparents)]


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


##### 4/17/2023 AL: Code below works fine
# modify, pull go term relationships


  
bpterms = lapply(uni.terms$GO.ID, get_parents, all_paths = FALSE)

#uni.terms.temp <- c("GO:1990091", "GO:0014718", "GO:0070315", "GO:1990092")

bpres   = data.frame(class=unlist(bpterms))

###write/write data (annoyingly the only way I can get it to work, but at least it works)

write.table(bpres, "/Users/andy/KW/11_wgcna/data/temp.csv", sep=",", col.names=TRUE, row.names=FALSE)
bpres=read.table("/Users/andy/KW/11_wgcna/data/temp.csv", header=TRUE, sep=",")


##### AL 4/17/23: The following code works with metacoder 0.3.6

## bpdata contains the GO IDs for the terminal nodes
## bpres and bpdata contain the same information, different formats = path information for each terminal node (GO Terms)

#parse GO data

testdata <- parse_tax_data(bpres, class_sep = ";")

#  mfterms = lapply(uni.terms$GO.ID, get_parents, all_paths = FALSE)
#  mfres   <- data.frame(class=unlist(mfterms))
#  write.table(mfres, "temp.csv", sep=",", col.names=TRUE, row.names=FALSE)
#  mfres=read.table("temp.csv", header=TRUE, sep=",")
#  testdata <- parse_tax_data(mfres, class_sep = ";")

# data = separate(bpres, c("class" = -1) , c("1","2","3","4","5","6","7","8","9","10","11") , sep = '\\|')

##### AMH 4/14/23: big warning about deprecated command here; looks like the old version of metacoder uses a deprecated dplyr command inside the parse_taxonomy_table command -- maybe why my output figure looks so shitty #####

# parse GO data (parse_taxonomy_table replaced by taxa::parse_tax_data)
# data <- parse_tax_data(bpres,
#                        class_sep = "\\|")


#create figure
# tempdata = filter_taxa(testdata, n_supertaxa <= 5) #filters terms that are WAAAAY out from the middle
## n_supertaxa sets # of nodes that can be in a single path from center to terminal nodes (cuts from terminal end, not internal nodes)
##### AMH 4/14/23: had to update lazyeval package here #####

#pdf(paste("",module,"meta_labels.pdf", sep=""), width=5, height=5, useDingbats=FALSE, onefile=FALSE)
#grey70      grey40      dodgerblue2 tomato2  chartreuse3  darkorchid2 yellow3 orange2 magenta1
#par(bg=NA)
set.seed(100) ## seeds: red = 50, thistle = 100, 105 for DEGs 

heat_tree(testdata, 
          node_label = taxon_names, ### modified by AL 
          # node_size = colsandsize$ngenes,
          # node_size_trans = "log10",
          node_size_range = c(0.01, 0.01),
          # node_label_size_trans = "log10",
          node_label_size_range = c(0.01, 0.01),
          # edge_size_trans = "log10",
          edge_size_range = c(0.004, 0.004),
          node_color = "blue",
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
