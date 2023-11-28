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
#ref <- all_nocross 
ref <- deg_nocross 

## Add sample type before merging 
sample.meta(ref)$indiv              <- sub("^(\\w+)\\..*", "\\1",  sample.meta(ref)$sample_ID)
sample.meta(ref)$sample_type        <- "reference"
sample.meta(tae_dat)$sample_type    <- "mixture"
sample.meta(tae_dat)$indiv          <- sample.meta(tae_dat)$sample_ID

### merge data 
merge <- merge_snpRdata(ref, tae_dat, all.x.snps=TRUE, all.y.snps=FALSE)
#sample.meta(merge)

### Format the reference and mixture for rubias
kw_reference <- filter_snps(merge[sample_type="reference"], non_poly=FALSE, bi_al=FALSE)
kw_reference                      <- format_snps(kw_reference, facets="indiv", "rafm") 
kw_reference$sample_type          <- "reference"
kw_reference$collection           <- as.character(sample.meta(ref)$site)
kw_reference$indiv                <- as.character(kw_reference$indiv)
# kw_reference$repunit              <- as.character(sample.meta(ref)$site)
kw_reference <- kw_reference %>%
  mutate(repunit = case_when(
    sample.meta(ref)$site %in% "MON" ~ "Expanded Range",
    sample.meta(ref)$site %in% c("NAP","POL") ~ "Historical Range" 
  ))
kw_reference <- kw_reference %>%
  select(indiv, sample_type, repunit, collection, everything())

tae_mixture  <- filter_snps(merge[sample_type="mixture"], non_poly=FALSE, bi_al=FALSE)
tae                      <- format_snps(tae_mixture, facets="indiv", "rafm") 
tae$sample_type          <- "mixture"
tae$repunit              <- NA
tae$collection           <- paste0("TAE_", as.character(sample.meta(tae_dat)$site))


tae$indiv                <- as.character(tae$indiv)
tae <- tae %>%
  select(indiv, sample_type, repunit, collection, everything())

# write.csv(tae, "tae_rubias_genotypes.txt")
# write.csv(kw_reference, "kw_ref_rubias_genotypes.txt")
# kw_tae_merged <- rbind(kw_reference, tae)

# kw_reference <- read.csv("kw_ref_rubias_genotypes.txt")
# tae_mixture  <- read.csv("tae_rubias_genotypes.txt")

########################################################################
# perform genetic mixture analysis (following parameters by D'Aloia et al 2021)
########################################################################
###### test mixture analyses code below
mix_kw_tae <- infer_mixture(reference = kw_reference, 
                            mixture = tae, 
                            gen_start_col = 5 )   

saveRDS(mix_kw_tae, "kw_allsnps_rubias_result.rds")
# saveRDS(mix_kw_tae, "kw_degsnps_rubias_result.rds")

mix_kw_tae <- readRDS("kw_allsnps_rubias_result.rds")
# mix_kw_tae <- readRDS("kw_degsnps_rubias_result.rds")

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
  mutate(site = sub("^TAE_", "", mixture_collection)) %>%
  summarise(percentage_match = sum(site == collection) / n())
percent_matches_to_site

## calculation for rep unit 
tae_dic_mon_expanded <- kw_tae_rep_indiv_topassign_range %>%
filter(mixture_collection %in% c("TAE_DIC", "TAE_MON")) %>%
group_by(mixture_collection) %>%
summarize(percentage_match_to_range = sum(repunit == "Expanded Range") / n())

tae_pol_nap_historical <- kw_tae_rep_indiv_topassign_range %>%
  filter(mixture_collection %in% c("TAE_POL", "TAE_NAP")) %>%
  group_by(mixture_collection) %>%
  summarize(percentage_match_to_range = sum(repunit == "Historical Range") / n())

percent_match_to_range <- rbind(tae_dic_mon_expanded,tae_pol_nap_historical)

### result tables: 
site_dat  <- percent_matches_to_site %>% mutate(assignment_type = "site")
range_dat <- percent_match_to_range %>% mutate(assignment_type = "range")

allsnp_dat  <- bind_rows(site_dat, range_dat) %>% mutate(loci = "All SNPs")
#degsnps_dat <- bind_rows(site_dat, range_dat) %>% mutate(loci = "DEG SNPs")
merged_data <- as.data.frame(bind_rows(allsnp_dat, degsnps_dat))

merged_data

### make bar plots 
merged_data_comb <- merged_data %>%
  group_by(mixture_collection, assignment_type, loci) %>%
  summarize(percentage = ifelse(all(is.na(percentage_match)), 
                                first(percentage_match_to_range), 
                                first(percentage_match)), .groups = "drop")
merged_data_comb <- merged_data_comb %>% mutate(incorrect_percentage=1-percentage)

## Make bar plots for each location 


mon <-ggplot(merged_data_comb %>% filter(mixture_collection == "TAE_MON"), aes(x=loci, y= percentage, fill = loci, pattern=assignment_type)) + 
  geom_bar_pattern(stat = "identity", position = "identity") +
  geom_text(aes(label=round(percentage, digits=2)), vjust=-.3, size=7) +
  scale_y_continuous(limits=c(0,1)) + 
  scale_fill_manual(values = c("DEG SNPs" = "darkorange2", "All SNPs" = "black")) + 
  scale_pattern_manual(values=c(range="stripe", site="none")) +
  facet_wrap(~ factor(assignment_type, levels=c("range", "site"))) +
  theme_void() + 
  theme(text = element_text(size = 16), 
        legend.position = "none", 
        strip.text=element_blank(),
        # panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))

dic <- ggplot(merged_data_comb %>% filter(mixture_collection == "TAE_DIC"), aes(x=loci, y= percentage, fill = loci, pattern=assignment_type)) + 
  geom_bar_pattern(stat = "identity", position = "identity") +
  geom_text(aes(label=round(percentage, digits=2)),vjust=-.04, size=7) +
  scale_y_continuous(limits=c(0,1)) + 
  scale_fill_manual(values = c("DEG SNPs" = "darkorange2", "All SNPs" = "black")) + 
  scale_pattern_manual(values=c(range="stripe", site="none")) +
  facet_wrap(~ factor(assignment_type, levels=c("range", "site"))) +
  theme_void() + 
  theme(text = element_text(size = 16), 
        legend.position = "none", 
        strip.text=element_blank(), 
        # panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))

nap <- ggplot(merged_data_comb %>% filter(mixture_collection == "TAE_NAP"), aes(x=loci, y= percentage, fill = loci, pattern=assignment_type)) + 
  geom_bar_pattern(stat = "identity", position = "identity") +
  geom_text(aes(label=round(percentage, digits=2)), vjust=-.3, size=7) +
  scale_y_continuous(limits=c(0,1)) + 
  scale_fill_manual(values = c("DEG SNPs" = "darkorange2", "All SNPs" = "black")) + 
  scale_pattern_manual(values=c(range="stripe", site="none")) +
  facet_wrap(~ factor(assignment_type, levels=c("range", "site"))) +
  theme_void() + 
  theme(text = element_text(size = 16), 
        legend.position = "none", 
        strip.text=element_blank(), 
        # panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))
        
pol <- ggplot(merged_data_comb %>% filter(mixture_collection == "TAE_POL"), aes(x=loci, y= percentage, fill = loci, pattern=assignment_type)) + 
  geom_bar_pattern(stat = "identity", position = "identity") +
  geom_text(aes(label=round(percentage, digits=2)), vjust=-.3, size=7) +
  scale_y_continuous(limits=c(0,1)) + 
  scale_fill_manual(values = c("DEG SNPs" = "darkorange2", "All SNPs" = "black")) + 
  scale_pattern_manual(values=c(range="stripe", site="none")) +
  facet_wrap(~ factor(assignment_type, levels=c("range", "site"))) +
  theme_void() + 
  theme(text = element_text(size = 16), 
        legend.position = "none", 
        strip.text=element_blank(), 
        # panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))

ggsave("mon.png", mon, bg = "transparent", width=4, height=4)
ggsave("nap.png", nap, bg = "transparent", width=4, height=4)
ggsave("pol.png", pol, bg = "transparent", width=4, height=4)
ggsave("dic.png", dic, bg = "transparent", width=4, height=4)


## addtional pngs for presentation
## mon and nap with site only: 

mon_site <- ggplot(merged_data_comb %>% filter(mixture_collection == "TAE_MON", assignment_type == "site"), aes(x=loci, y= percentage, fill = loci, pattern=assignment_type)) + 
  geom_bar_pattern(stat = "identity", position = "identity") +
  geom_text(aes(label=round(percentage, digits=2)), vjust=-.3, size=7) +
  scale_y_continuous(limits=c(0,1)) + 
  scale_fill_manual(values = c("DEG SNPs" = "darkorange2", "All SNPs" = "black")) + 
  scale_pattern_manual(values=c(range="stripe", site="none")) +
  # facet_wrap(~ factor(assignment_type, levels=c("range", "site"))) +
  theme_void() + 
  theme(text = element_text(size = 16), 
        legend.position = "none", 
        strip.text=element_blank(),
        # panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))

nap_site <- ggplot(merged_data_comb %>% filter(mixture_collection == "TAE_NAP", assignment_type == "site"), aes(x=loci, y= percentage, fill = loci, pattern=assignment_type)) + 
  geom_bar_pattern(stat = "identity", position = "identity") +
  geom_text(aes(label=round(percentage, digits=2)), vjust=-.3, size=7) +
  scale_y_continuous(limits=c(0,1)) + 
  scale_fill_manual(values = c("DEG SNPs" = "darkorange2", "All SNPs" = "black")) + 
  scale_pattern_manual(values=c(range="stripe", site="none")) +
  facet_wrap(~ factor(assignment_type, levels=c("range", "site"))) +
  theme_void() + 
  theme(text = element_text(size = 16), 
        legend.position = "none", 
        strip.text=element_blank(), 
        # panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))

dic <- ggplot(merged_data_comb %>% filter(mixture_collection == "TAE_DIC"), aes(x=loci, y= percentage, fill = loci, pattern=assignment_type)) + 
  geom_bar_pattern(stat = "identity", position = "identity") +
  geom_text(aes(label=round(percentage, digits=2)),vjust=-.04, size=7) +
  scale_y_continuous(limits=c(0,1)) + 
  scale_fill_manual(values = c("DEG SNPs" = "darkorange2", "All SNPs" = "black")) + 
  scale_pattern_manual(values=c(range="stripe", site="none")) +
  facet_wrap(~ factor(assignment_type, levels=c("range", "site"))) +
  theme_void() + 
  theme(text = element_text(size = 16), 
        legend.position = "none", 
        strip.text=element_blank(), 
        # panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))
pol <- ggplot(merged_data_comb %>% filter(mixture_collection == "TAE_POL"), aes(x=loci, y= percentage, fill = loci, pattern=assignment_type)) + 
  geom_bar_pattern(stat = "identity", position = "identity") +
  geom_text(aes(label=round(percentage, digits=2)), vjust=-.3, size=7) +
  scale_y_continuous(limits=c(0,1)) + 
  scale_fill_manual(values = c("DEG SNPs" = "darkorange2", "All SNPs" = "black")) + 
  scale_pattern_manual(values=c(range="stripe", site="none")) +
  facet_wrap(~ factor(assignment_type, levels=c("range", "site"))) +
  theme_void() + 
  theme(text = element_text(size = 16), 
        legend.position = "none", 
        strip.text=element_blank(), 
        # panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA))












###   Explore accuracy of marker panels ### ====
## Self Assignment ====
sa_kw <- self_assign(reference = kw_reference, gen_start_col = 5)
head(sa_kw, n = 190)

sa_to_site <- sa_chinook %>%
  group_by(indiv, collection, repunit, inferred_collection) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood)) %>% 
  group_by(indiv) %>%
  filter(repu_scaled_like==max(repu_scaled_like)) 

sa_to_range <- sa_chinook %>%
  group_by(indiv, collection, repunit, inferred_repunit) %>%
  summarise(repu_scaled_like = sum(scaled_likelihood)) %>% 
  group_by(indiv) %>%
  filter(repu_scaled_like==max(repu_scaled_like)) 



### Run LOO models ====
rep_prop <- kw_reference  %>%
  group_by( repunit)  %>%
  summarize(repprop =n() ) %>%
  mutate(proportion = repprop/sum(repprop))

#Use these proportions to supply Dirichlet random variable
arep <- rep_prop %>%
  ungroup() %>%
  mutate(dirichlet = 10 * proportion) %>%
  select(repunit, dirichlet)
arep

set.seed(100)
kw_loo_sims <- assess_reference_loo(reference = kw_reference, 
                                    gen_start_col = 5, 
                                    reps = 100,            #How many simulated mixtures to create
                                    mixsize = 500,         #Num individuals in each simulated mixture 
                                    alpha_repunit = arep, 
                                    return_indiv_posteriors = TRUE)

#saveRDS(kw_loo_sims, "deg_kw_loo_sims.RDS")
kw_loo_sims <- readRDS('deg_kw_loo_sims.RDS')

### Retrieving the individual simulated fish posteriors
# print out the indiv posteriors
kw_loo_sims$indiv_posteriors


# summarise things by rep unit
collection_pof <- kw_loo_sims$indiv_posteriors %>%
  filter(collection == simulated_collection) %>%
  group_by(iter, indiv, simulated_collection, collection) %>%  # first aggregate over reporting units
  summarise(collection_PofZ = sum(PofZ)) %>%
  ungroup() %>%
  arrange(collection, simulated_collection) %>%
  mutate(simulated_collection = factor(simulated_collection, levels = unique(simulated_collection)))

# summarise things by collection
repu_pofzs <- kw_loo_sims$indiv_posteriors %>%
  filter(repunit == simulated_repunit) %>%
  group_by(iter, indiv, simulated_collection, repunit) %>%  # first aggregate over reporting units
  summarise(repu_PofZ = sum(PofZ)) %>%
  ungroup() %>%
  arrange(repunit, simulated_collection) %>%
  mutate(simulated_collection = factor(simulated_collection, levels = unique(simulated_collection)))

# also get the number of simulated individuals from each collection
num_simmed <- kw_loo_sims$indiv_posteriors %>%
  group_by(iter, indiv) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::count(simulated_collection)

# note, the last few steps make simulated collection a factor so that collections within
# the same repunit are grouped together in the plot.

# now, plot it
p <- ggplot(repu_pofzs, aes(x = simulated_collection, y = repu_PofZ)) +
  geom_boxplot(aes(colour = repunit)) +
  geom_text(data = num_simmed, mapping = aes(y = 1.025, label = n), angle = 90, hjust = 0, vjust = 0.5, size = 3) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 9, vjust = 0.5)) +
  ylim(c(NA, 1.05)) +
  scale_color_viridis(discrete=TRUE, name="Range") +
  theme_bw()+
  labs(x = "Simulated Reporting Unit", y = "Probability of assignment to colletion site")

ggsave("~/KW/figures/supplemental /LOO_result.png", p, width=5, height=5, units="in" )

#Show  values below 0.9 for each simulatedcollection 
num_below09 <- kw_loo_sims$indiv_posteriors %>%
  filter(collection == simulated_collection) %>%
  group_by(iter, indiv, simulated_collection, collection) %>%  # first aggregate over collections
  summarise(collection_PofZ = sum(PofZ))  %>%
  filter(collection_PofZ < 0.95) %>%
  arrange(collection, collection_PofZ) %>%
  group_by(collection) %>%
  tally 

