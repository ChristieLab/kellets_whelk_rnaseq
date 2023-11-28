# Get all shared module terms 
library(dplyr)
library(stringr)
file_dir <- "/Users/andy/KW/11_wgcna/enriched_go_bp/"
file_names <- list.files(file_dir, pattern = "_module_all_enrichedterms.csv", full.names = TRUE)
enriched_list <- lapply(file_names, read.csv)
go_id_list <- lapply(enriched_list, function(df) df$GO.ID)
common_characters <- Reduce(intersect, go_id_list)

# significant modules after boniferrei correction = red and thistle 
sig_file_names <- c("/Users/andy/KW/11_wgcna/enriched_go_bp/red_module_all_enrichedterms.csv", "/Users/andy/KW/11_wgcna/enriched_go_bp/thistle_module_all_enrichedterms.csv")

enriched_list <- lapply(sig_file_names, read.csv)
go_id_list <- lapply(enriched_list, function(df) df$GO.ID)
shared_terms <- Reduce(intersect, go_id_list)
shared_dat <- filter(enriched_list[[1]], GO.ID %in% shared_terms )
colnames(shared_dat)[colnames(shared_dat) == "result1"] <- "P-value"



write.csv(shared_dat, "shared_terms_thistleandred.csv")

shared_dat
sig_shared_dat <- shared_dat[shared_dat$result1 < 0.001, ]

### find unique terms in each 
red.go <- read.csv("/Users/andy/KW/11_wgcna/enriched_go_bp/red_module_all_enrichedterms.csv")
thistle.go <- read.csv("/Users/andy/KW/11_wgcna/enriched_go_bp/thistle_module_all_enrichedterms.csv")

unique_red <- red.go[!(red.go$GO.ID %in% thistle.go$GO.ID), ]
sig_unique_red <- unique_red[unique_red$result1 < 0.001, ]
nrow(sig_unique_red)

unique_thistle <- thistle.go[!(thistle.go$GO.ID %in% red.go$GO.ID), ]
sig_unique_thistle <- unique_thistle[thistle.go$result1 < 0.001, ]
nrow(sig_unique_thistle)

unique_red <- anti_join(red.go, thistle.go, by = "GO.ID")
unique_thistle <- anti_join(red.go, thistle.go, by = "GO.ID")