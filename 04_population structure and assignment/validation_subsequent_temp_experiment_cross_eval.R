### This is a script written by Andy Lee; last modified Aug 25, 2023
### Conduct population assignments on TAE and Maternal Cross data 
### 

### Load packages and set wd 
library("snpR")
library("ranger")
library(forestError)
library(dplyr)
setwd("~/KW/12_validation")

### Load in datasets 
tae_sample_meta <- read.csv("~/KW/12_validation/tae_sampleinfo.csv", na.strings = "#N/A")
tae_sample_meta <- tae_sample_meta[!is.na(tae_sample_meta$sample_ID), ]
tae_dat <- import.snpR.data("~/KW/12_validation/kw_tae_raw.vcf.gz", sample.meta = tae_sample_meta)

sample_meta <- read.csv("~/KW/09_angsd/neutral/sample_info.csv")
all_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_neutral.vcf", sample.meta = sample_meta)
### filter by maf within each site 
all_filt <- filter_snps(all_vcf, maf = 0.05, maf_facets = "site") 
deg_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_deg.vcf", sample.meta = sample_meta)
deg_filt <- filter_snps(deg_vcf, maf = 0.05, maf_facets = "site") 

sample.meta(tae_dat)$experiment <- "tae"
sample.meta(all_filt)$experiment <- "maincross"
sample.meta(deg_filt)$experiment <- "maincross"


# comb <- merge_snpRdata(tae_dat, all_filt, by.snp=c("CHROM","position"), all.x.snps = FALSE, all.y.snps = FALSE)
# comb <- merge_snpRdata(tae_dat, deg_filt, by.snp=c("CHROM","position"), all.x.snps = FALSE, all.y.snps = FALSE)

comb.nocross <- filter_snps(comb[site=-c("NAPxMON")])
comb.tae <-   comb.nocross[experiment="tae"]
comb.maincross <- comb.nocross[experiment="maincross"]

### Run with only MON/NAP 
#comb.monnap <- comb[site=c("MON","NAP")]
#comb.tae <-   comb.monnap[experiment="TAE"]
#comb.maincross <- comb.monnap[experiment="cross"]

### Running random forest in ranger 
comb <- merge_snpRdata(tae_dat, deg_filt, by.snp=c("CHROM","position"), all.x.snps = FALSE, all.y.snps = FALSE)
comb.nocross <- filter_snps(comb[site=-c("NAPxMON")])
comb.tae <-   comb.nocross[experiment="tae"]
# comb.tae.main <- comb.tae[site=c("MON","NAP")]
comb.maincross <- comb.nocross[experiment="maincross"]

maincross_sn <- format_snps(comb.maincross, output="sn")
maincross_train   <- t(maincross_sn[,-c(1:12)])
maincross_sites   <- as.factor(sample.meta(comb.maincross)$site)

tae_sn <- format_snps(comb.tae.main, output="sn")
tae_test    <- t(tae_sn[,-c(1:12)])
tae_sites   <- as.factor(sample.meta(comb.tae.main)$site)


rf <- ranger(x=maincross_train, importance="permutation", keep.inbag = TRUE, y=maincross_sites, num.trees = 1000000)

eval <- predict(rf, tae_test, num.trees = 1000000)
table(obs=tae_sites, pred=eval$predictions)

### Calculate error rate 
err <- forestError::quantForestError(rf, # the forest 
                                      X.train = maincross_train, # the training data
                                      X.test =  tae_test, # the test data
                                      Y.train = maincross_sites)  # the classifications of the training data
cbind(tae_sites,err)

### Summarize results 
eval.dat <- as.data.frame(table(obs=tae_sites, pred=eval$predictions))
eval.dat$obs <- as.character(eval.dat$obs)
eval.dat$pred <- as.character(eval.dat$pred)

result <- eval.dat %>%
  group_by(obs) %>%
  summarize(correct_percent = sum(Freq[pred == obs]) / sum(Freq))
result

# Calculate the correct assignment for each range
eval.dat <- eval.dat %>% 
  mutate(
    PredRange = case_when(
      pred %in% c("MON", "DIC") ~ "Expanded Range", 
      pred %in% c("NAP", "POL") ~ "Historical Range", 
      ), 
    ObsRange = case_when(
      obs %in% c("MON", "DIC") ~ "Expanded Range", 
      obs %in% c("NAP", "POL") ~ "Historical Range", 
    )
  )

# Calculate correct assignments for each range
ranger_correct_assignments <- eval.dat %>%
  group_by(obs) %>%
  summarize(total_correct = sum(Freq[ObsRange == PredRange]),
            total_freq = sum(Freq)) %>%
mutate(percentage_correct = total_correct / total_freq * 100)
ranger_correct_assignments

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
###     train with maincross and predict with NAPxMON           
###     calculate error rates
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

all.maincross <- deg_filt[site=-c("NAPxMON")]
all.testcross <- deg_filt[site="NAPxMON"]

maincross_sn <- format_snps(all.maincross, output="sn")
maincross_train <- t(maincross_sn[,-c(1:7)])
maincross_sites <- as.factor(sample.meta(all.maincross)$site)

testcross_sn <- format_snps(all.testcross, output="sn")
testcross_test <- t(testcross_sn[,-c(1:7)])
testcross_sites = as.factor(sample.meta(all.testcross)$site)

rf_all <- ranger(x=maincross_train, importance="permutation", keep.inbag = TRUE, y=maincross_sites, num.trees = 1000000)

eval <- predict(rf_all, testcross_test, num.trees = 1000000)
table(obs=testcross_sites, pred=eval$predictions)

err <- forestError::quantForestError(rf_all, # the forest 
                                     X.train = maincross_train, # the training data
                                     X.test =  testcross_test, # the test data
                                     Y.train = maincross_sites)  # the classifications of the training data

# compare err rate with maincross predictions to TAE ======================================
sample.meta(tae_dat)$experiment <- "tae"
sample.meta(deg_filt)$experiment <- "maincross"

comb <- merge_snpRdata(tae_dat, deg_filt, by.snp=c("CHROM","position"), all.x.snps = FALSE, all.y.snps = FALSE)
comb.nocross <- filter_snps(comb[site=-("NAPxMON")])
comb.tae <-   comb.nocross[experiment="tae"]
comb.tae.main <- comb.tae[site=c("MON","NAP")]
comb.maincross <- comb.nocross[experiment="maincross"]

maincross_sn <- format_snps(comb.maincross, output="sn")
maincross_train   <- t(maincross_sn[,-c(1:12)])
maincross_sites   <- as.factor(sample.meta(comb.maincross)$site)

tae_sn <- format_snps(comb.tae.main, output="sn")
tae_test    <- t(tae_sn[,-c(1:12)])
tae_sites   <- as.factor(sample.meta(comb.tae.main)$site)

rf_tae <- ranger(x=maincross_train, importance="permutation", keep.inbag = TRUE, y=maincross_sites, num.trees = 100000)

err1 <- forestError::quantForestError(rf_tae, # the forest 
                                     X.train = maincross_train, # the training data
                                     X.test =  tae_test, # the test data
                                     Y.train = maincross_sites)  # the classifications of the training data
cbind(tae_sites,err1)



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### train with TAE and predict with cross, can we predict POL?  ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
rf_tae <- run_random_forest(comb.tae, response = "site", importance = "permutation")
rf_tae$models$.base_.base$predictions

sn_cross <- format_snps(comb.maincross, output="sn")
predict_data_cross <- cbind.data.frame(site=sample.meta(comb.maincross)$site, t(sn_cross[, 13:ncol(sn_cross)]))

eval_tae <- predict(rf_tae$models$.base_.base$model, predict_data_cross)
table(obs=predict_data_cross$site, pred=eval_tae$predictions) 




### Rubias cross validation 
### RUBIAS ====
library(rubias)
library(tidyverse)
library(snpR)
library(ggpattern)

setwd("~/KW/13_assignment/rubias")
### format input data 
### the data frame must have these 4 columns: sample_type, repunit, collection, indiv 

## Transcriptome-wide snps: 
kw_sample_meta <- read.csv("~/KW/09_angsd/neutral/sample_info.csv")
all_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_neutral.vcf", sample.meta = kw_sample_meta)
all_nocross <- filter_snps(all_vcf[site=-"NAPxMON"], maf=0.05, maf_facets = "site")
# snp.meta(all_nocross)$ID <- paste0(snp.meta(all_nocross)$CHROM, "_",snp.meta(all_nocross)$position)

# all_monnap <- filter_snps(all_vcf[site=-c("NAPxMON","POL")], maf = 0.05, maf_facets = "site")

## DEG SNPs 
deg_vcf <- import.snpR.data("~/KW/09_angsd/vcfs/kw_deg.vcf", sample.meta = kw_sample_meta)
deg_nocross <- filter_snps(deg_vcf[site=-"NAPxMON"], maf = 0.05, maf_facets = "site")


## load TAE data (to be assigned)
tae_sample_meta <- read.csv("~/KW/12_validation/tae_sampleinfo.csv", na.strings = "#N/A")
tae_sample_meta <- tae_sample_meta[!is.na(tae_sample_meta$sample_ID), ]
tae_sample_meta_short <- tae_sample_meta[, 1:4]
tae_dat <- import.snpR.data("~/KW/12_validation/kw_tae_raw.vcf.gz", sample.meta = tae_sample_meta_short)

## Set reference 
# ref <- all_nocross 
kw <- deg_nocross 

## Add sample type before merging 
sample.meta(kw)$indiv              <- sub("^(\\w+)\\..*", "\\1",  sample.meta(kw)$sample_ID)
sample.meta(kw)$sample_type        <- "mixture"
sample.meta(tae_dat)$sample_type    <- "reference"
sample.meta(tae_dat)$indiv          <- sample.meta(tae_dat)$sample_ID

### merge data 
merge <- merge_snpRdata(kw, tae_dat, all.x.snps=TRUE, all.y.snps=FALSE)
#sample.meta(merge)

### Format the reference and mixture for rubias
tae_reference <- filter_snps(merge[sample_type="reference"], non_poly=FALSE, bi_al=FALSE)
tae_reference                      <- format_snps(tae_reference, facets="indiv", "rafm") 
tae_reference$sample_type          <- "reference"
tae_reference$collection           <- as.character(sample.meta(tae_dat)$site)
tae_reference$indiv                <- as.character(tae_reference$indiv)
# tae_reference$repunit              <- as.character(sample.meta(ref)$site)
tae_reference <- tae_reference %>%
  mutate(repunit = case_when(
    sample.meta(tae_dat)$site %in% c("MON","DIC") ~ "Expanded Range",
    sample.meta(tae_dat)$site %in% c("NAP","POL") ~ "Historical Range" 
  ))
tae_reference <- tae_reference %>%
  select(indiv, sample_type, repunit, collection, everything())

kw_mixture  <- filter_snps(merge[sample_type="mixture"], non_poly=FALSE, bi_al=FALSE)
kw                     <- format_snps(kw_mixture, facets="indiv", "rafm") 
kw$sample_type          <- "mixture"
kw$repunit              <- NA
kw$collection           <- paste0("kw_", as.character(sample.meta(kw)$site))
kw$indiv                <- as.character(kw$indiv)

kw <- kw %>%
  select(indiv, sample_type, repunit, collection, everything())


###### mixture analyses code below
mix_kw_tae <- infer_mixture(reference = tae_reference, 
                            mixture = kw, 
                            gen_start_col = 5 )   

## Mixture among tae individuals 
mixprop <- mix_kw_tae$mixing_proportions %>%
  group_by(mixture_collection, repunit) %>%
  summarise(repprop = sum(pi)) %>%
  #left_join(mix_kw_tae$bootstrapped_proportions) %>%
  #ungroup() %>%
  arrange(desc(repprop))
mixprop

#Individual assignment probabilities to each collection site 
kw_tae_rep_indiv_ests_site <- mix_kw_tae$indiv_posterior %>%
  filter(mixture_collection != "NAPxMON") %>%
  group_by(mixture_collection, indiv, collection) %>%
  summarise(rep_pofz = sum(PofZ)) 

kw_tae_rep_indiv_topassign_site <- kw_tae_rep_indiv_ests_site  %>% #Return top assignment for each indiv
  group_by(indiv) %>%
  filter(rep_pofz==max(rep_pofz)) #%>%

#Individual assignment probabilities to each reporting unit
kw_tae_rep_indiv_ests_range <- mix_kw_tae$indiv_posterior %>%
  group_by(mixture_collection, indiv, repunit) %>%
  summarise(rep_pofz = sum(PofZ))

kw_tae_rep_indiv_topassign_range <- kw_tae_rep_indiv_ests_range  %>% #Return top assignment to range for each indiv
  group_by(indiv) %>%
  filter(rep_pofz==max(rep_pofz)) 

# Calculate the percentage of matches between "inferred" and "site" within each site
grouped_data <- kw_tae_rep_indiv_topassign_site %>%
  group_by(mixture_collection)

percent_matches_to_site <- grouped_data %>%
  mutate(site = sub("^kw_", "", mixture_collection)) %>%
  summarise(percentage_match = sum(site == collection) / n())
percent_matches_to_site

## calculation for rep unit 
tae_dic_mon_expanded <- kw_tae_rep_indiv_topassign_range %>%
  filter(mixture_collection %in% c("kw_MON")) %>%
  group_by(mixture_collection) %>%
  summarize(percentage_match_to_range = sum(repunit == "Expanded Range") / n())

tae_pol_nap_historical <- kw_tae_rep_indiv_topassign_range %>%
  filter(mixture_collection %in% c("kw_POL", "kw_NAP")) %>%
  group_by(mixture_collection) %>%
  summarize(percentage_match_to_range = sum(repunit == "Historical Range") / n())

percent_match_to_range <- rbind(tae_dic_mon_expanded,tae_pol_nap_historical)

### result tables: 
site_dat  <- percent_matches_to_site %>% mutate(assignment_type = "site")
range_dat <- percent_match_to_range %>% mutate(assignment_type = "range")

site_dat
range_dat


### Archive of older code: ====
#### Running random forest in snpR written with Will in June 2023
rf <- run_random_forest(comb.maincross, response = "site", importance = "permutation", num.trees = 1000000)
rf$models$.base_.base$predictions
table(rf$models$.base_.base$predictions)

sn_tae <- format_snps(comb.tae, output="sn")
predict_data <- cbind.data.frame(site=sample.meta(comb.tae)$site, t(sn_tae[, 13:ncol(sn_tae)]))

rf$models$.base_.base$model$forest$independent.variable.names               
colnames(predict_data)

eval <- predict(rf$models$.base_.base$model, predict_data)
table(obs=predict_data$site, pred=eval$predictions)