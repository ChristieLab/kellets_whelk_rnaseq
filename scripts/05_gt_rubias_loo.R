# Script to explore gtseq data 
library(snpR); library(ggplot2); library(ranger); library(forestError); library(rubias); library(dplyr)
setwd("/Users/andy/gtseq")

#### Load in Data ====
## Get sample info 
sample_meta <- read.csv("/Users/andy/gtseq/gtseek_sampleinfo.csv")
site_info <- read.csv("/Users/andy/gtseq/gtseek_site_info.csv")
adult_or_recruit <- read.csv("/Users/andy/gtseq/gtseek_sampleinfo_age.csv")

sample_meta$cap_label <-  as.character(gsub(".*([A-Za-z]\\d+)$", "\\1", sample_meta$Sample))
sample_meta$Plate <- as.numeric(gsub(".*KW[_-]([0-9]+).*", "\\1", sample_meta$Sample))
sample_meta$Site_Letter  <- as.character((gsub(".*([A-Za-z])\\d*$", "\\1", sample_meta$Sample)))
sample_meta_merged <- merge(sample_meta, site_info, by.x="Site_Letter", by.y="Site_Letter", all.x=TRUE)
sample_meta_merged_age <- merge(sample_meta_merged, adult_or_recruit, by.x="cap_label", by.y="cap_label", all.x=TRUE)

## Load genotypes 
run1_genotypes  <- read.csv("/Users/andy/gtseq/Whelk_Prod1_Genotypes.csv", row.names=1)
rerun_genotypes <- read.csv("/Users/andy/gtseq/Whelk_P14-rerun_Genotypes.csv", row.names=1)

all_genotypes <- rbind(run1_genotypes,rerun_genotypes)
genotypes <- t(all_genotypes[, -(1:5)])

## Create snpR object 
allpop    <- import.snpR.data(genotypes, mDat = "00", rows_per_individual = 1, sample.meta = sample_meta_merged_age)
snp.meta(allpop)$SNP    <- row.names(t(all_genotypes[, -(1:5)]))
snp.meta(allpop)$Contig  <- gsub("(_i\\d+)_\\d+$", "\\1", snp.meta(allpop)$SNP)
adult_allpop <- allpop[tissue_type = "Adult"]

saveRDS(adult_allpop, "adult_allpop.rds")

### get snps from maincross 
### load raw vcf files 
sample_meta <- read.csv("~/KW/09_angsd/neutral/sample_info.csv")
deg_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_deg.vcf", sample.meta = sample_meta)
deg_monnap <- filter_snps(deg_vcf[site=c("MON","NAP")], maf = 0.05, maf_facets = "site") ### filter by maf within each site 
snp.meta(deg_monnap)$ID <- paste0(snp.meta(deg_monnap)$CHROM, "_",snp.meta(deg_monnap)$position)

saveRDS(deg_monnap, "deg_monnap.rds")
## Crossing study pops only 
gt_monnap  <- adult_allpop[Site_Code=c("MON","NAR")]
saveRDS(gt_monnap, "gt_monnap.rds")

gtsnps    <- monnap[snp.meta(monnap)$SNP %in% snp.meta(deg_monnap)$ID] ### get shared snps 
gtsnps <- calc_pairwise_fst(gtsnps, facets="Site_Code")

gtsnps_fst <- get.snpR.stats(gtsnps, "Site_Code", "fst")
gtsnps_fst$pairwise[,c(2,3,4,6)]

### Compare with FST of the same SNPs in the RNA dataset
deg_monnap <-  calc_pairwise_fst(deg_monnap, facets="site")
deg_monnap_fst <- get.snpR.stats(deg_monnap, "site", "fst")
shared_snps_fst <- deg_monnap_fst$pairwise[deg_monnap_fst$pairwise$ID %in% gtsnps_fst$pairwise$SNP, ]
shared_snps_fst_monnap <- shared_snps_fst[1:6,c(2,5,10)]
shared_snps_fst_monnap[c(4,5,3,1,6,2),]

### Format Shared SNPs between GT data and main cross
sample.meta(gtsnps)$indiv          <- sample.meta(gtsnps)$Sample
gt_ref                      <- format_snps(gtsnps, facets="indiv", "rafm") 
gt_ref$sample_type          <- "reference"
gt_ref$collection           <- as.character(sample.meta(gtsnps)$Site_Code)
gt_ref$indiv                <- as.character(sample.meta(gtsnps)$indiv)
gt_ref$repunit              <- as.character(sample.meta(gtsnps)$Site_Code) # get the site to mutate 
gt_ref <- gt_ref %>% # use site to create range
  mutate(repunit = case_when(
    sample.meta(gtsnps)$Site_Code %in% "MON" ~ "Expanded Range",
    sample.meta(gtsnps)$Site_Code %in% c("NAR","POL") ~ "Historical Range" 
  ))
gt_ref <- gt_ref %>%
  select(indiv, sample_type, repunit, collection, everything())
gt_ref[gt_ref == ""] <- NA

### Format all 305 shared SNPs ====
crossing <- monnap
sample.meta(crossing)$indiv          <- sample.meta(crossing)$Sample
gt_all                      <- format_snps(crossing, facets="indiv", "rafm") 
gt_all$sample_type          <- "reference"
gt_all$collection           <- as.character(sample.meta(crossing)$Site_Code)
gt_all$indiv                <- sample.meta(crossing)$indiv 
gt_all$repunit              <- as.character(sample.meta(crossing)$Site_Code)
gt_all <- gt_all %>%
  mutate(repunit = case_when(
    sample.meta(crossing)$Site_Code %in% "MON" ~ "Expanded Range",
    sample.meta(crossing)$Site_Code %in% c("NAR","POL") ~ "Historical Range" 
  ))
gt_all <- gt_all %>%
  select(indiv, sample_type, repunit, collection, everything())
gt_all[gt_all == ""] <- NA

## Self Assignment ====
# sa_gt_snps   <- self_assign(reference = gt_ref, gen_start_col = 5)
# sa_gtcontig  <- self_assign(reference = gtcontig_ref, gen_start_col = 5)

sa_allgtsnps <- self_assign(reference = gt_all, gen_start_col = 5)

## Explore results with just shared SNPs 
dat <- sa_allgtsnps %>%
  select(indiv, collection, inferred_collection, log_likelihood, scaled_likelihood) %>%
  group_by(indiv) %>%
  filter(scaled_likelihood==max(scaled_likelihood)) %>%
  # filter(scaled_likelihood >= 0.45) %>%
  group_by(indiv) # %>%

result <- addmargins(table(dat$collection, dat$inferred_collection))[-3,]
percent <- rbind(result[1,1]/result[1,3], result[2,2]/result[2,3])
result
cbind(result, percent)

hist(dat$scaled_likelihood)

filtdat <- sa_allgtsnps %>%
  select(indiv, collection, inferred_collection, log_likelihood, scaled_likelihood) %>%
  group_by(indiv) %>%
  # filter(scaled_likelihood==max(scaled_likelihood)) %>%
  filter(scaled_likelihood >= 0.5) %>%
  group_by(indiv) # %>%

addmargins(table(filtdat$collection, filtdat$inferred_collection))[-3,]


### Bootstrap, using



### Using gt samples with shared snps between gtseq and maincross to make predictions on maincross samples ====
ref_gt_snps                      <- format_snps(gtsnps, facets="indiv", "rafm") 
ref_gt_snps$sample_type          <- "reference"
ref_gt_snps$repunit              <- as.character(sample.meta(gtsnps)$Site_Code)
ref_gt_snps <- ref_gt_snps %>% # use site to create range
  mutate(repunit = case_when(
    sample.meta(gtsnps)$Site_Code %in% "MON" ~ "Expanded Range",
    sample.meta(gtsnps)$Site_Code %in% c("NAR","POL") ~ "Historical Range" 
  ))

ref_gt_snps$collection           <- as.character(sample.meta(gtsnps)$Site_Code)
ref_gt_snps$indiv                <- as.character(sample.meta(gtsnps)$indiv)
ref_gt_snps <- ref_gt_snps %>%
  select(indiv, sample_type, repunit, collection, everything())
ref_gt_snps[ref_gt_snps == ""] <- NA

deg_shared_snps    <- deg_monnap[snp.meta(deg_monnap)$ID %in% snp.meta(monnap)$SNP ]
sample.meta(deg_shared_snps)$indiv              <- sub("^(\\w+)\\..*", "\\1",  sample.meta(deg_shared_snps)$sample_ID)
deg_mix                      <- format_snps(deg_shared_snps, facets="indiv", "rafm") 
deg_mix$sample_type          <- "mixture"
deg_mix$collection           <- paste0("mc_",as.character(sample.meta(deg_shared_snps)$site))
deg_mix$indiv                <- as.character(sample.meta(deg_shared_snps)$indiv)
deg_mix$repunit              <- NA # get the site to mutate 

deg_mix <- deg_mix %>%
  select(indiv, sample_type, repunit, collection, everything())
deg_mix[deg_mix == ""] <- NA


mix_mc_gt <- infer_mixture(reference = gt_ref, mixture = deg_mix, 
                            gen_start_col = 5)   

mixprop <- mix_mc_gt$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi)) %>%
  #left_join(mix_kw_gtseq$bootstrapped_proportions) %>%
  #ungroup() %>%
  arrange(desc(repprop))
mixprop

#Individual assignment probabilities to each collection site 
mix_main_gt_ests_site <- mix_mc_gt$indiv_posterior %>%
  # filter(mixture_collection != "NAPxMON") %>%
  group_by(mixture_collection, indiv, collection) %>%
  summarise(rep_pofz = sum(PofZ)) 

mix_main_gt_ests_site_topassign_site <- mix_main_gt_ests_site  %>% #Return top assignment for each indiv
  group_by(indiv) %>%
  filter(rep_pofz==max(rep_pofz)) #%>%

# Calculate the percentage of matches between "inferred" and "site" within each site
mix_main_gt_ests_site_topassign_site %>%
  group_by(mixture_collection) %>%
  mutate(site = sub("^mc_", "", mixture_collection)) %>%
  mutate(site = sub("NAR", "NAP", site)) %>%
  summarise(percentage_match = sum(site == collection) / n())


### Using maincrosswith shared snps between gtseq and maincross to make predictions on gtseq samples ====

deg_shared_snps    <- deg_monnap[snp.meta(deg_monnap)$ID %in% snp.meta(monnap)$SNP ]
sample.meta(deg_shared_snps)$indiv              <- sub("^(\\w+)\\..*", "\\1",  sample.meta(deg_shared_snps)$sample_ID)
deg_ref                      <- format_snps(deg_shared_snps, facets="indiv", "rafm") 
deg_ref$sample_type          <- "reference"
deg_ref$collection           <- as.character(sample.meta(deg_shared_snps)$site)
deg_ref$indiv                <- as.character(sample.meta(deg_shared_snps)$indiv)
deg_ref$repunit              <- as.character(sample.meta(deg_shared_snps)$site)
deg_ref              <- deg_ref %>% # use site to create range
  mutate(repunit = case_when(
    sample.meta(deg_shared_snps)$site %in% "MON" ~ "Expanded Range",
    sample.meta(deg_shared_snps)$site %in% "NAP" ~ "Historical Range" 
  ))
deg_ref <- deg_ref %>%
  select(indiv, sample_type, repunit, collection, everything())
deg_ref[deg_ref == ""] <- NA


gt_mix                      <- format_snps(gtsnps, facets="indiv", "rafm") 
gt_mix$sample_type          <- "mixture"
gt_mix$repunit              <- NA
gt_mix$collection           <- paste0("gt_",as.character(sample.meta(gtsnps)$Site_Code))
gt_mix$indiv                <- as.character(sample.meta(gtsnps)$indiv)
gt_mix <- gt_mix %>%
  select(indiv, sample_type, repunit, collection, everything())
gt_mix[gt_mix == ""] <- NA

mix_mc_gt<- infer_mixture(reference = deg_ref, mixture = gt_mix, 
                             gen_start_col = 5)   

mixprop <- mix_mc_gt$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi)) %>%
  #left_join(mix_kw_gtseq$bootstrapped_proportions) %>%
  #ungroup() %>%
  arrange(desc(repprop))
mixprop

#Individual assignment probabilities to each collection site 
mix_main_gt_ests_site <- mix_mc_gt$indiv_posterior %>%
  # filter(mixture_collection != "NAPxMON") %>%
  group_by(mixture_collection, indiv, collection) %>%
  summarise(rep_pofz = sum(PofZ)) 

mix_main_gt_ests_site_topassign_site <- mix_main_gt_ests_site  %>% #Return top assignment for each indiv
  group_by(indiv) %>%
  filter(rep_pofz==max(rep_pofz)) #%>%

# Calculate the percentage of matches between "inferred" and "site" within each site
mix_main_gt_ests_site_topassign_site %>%
  group_by(mixture_collection) %>%
  mutate(site = sub("^gt_", "", mixture_collection)) %>%
  mutate(site = sub("NAR", "NAP", site)) %>%
  summarise(percentage_match = sum(site == collection) / n())
