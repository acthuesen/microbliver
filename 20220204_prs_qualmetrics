# DEPENDENCIES AND WD ----
library(readxl)
library(tidyverse)
setwd('workdir') # where genetic data are and go

# init an empty dataframe ----
qualmetrics <- data.frame(matrix(ncol=8,nrow=0))
colnames(qualmetrics) <- c('pgsid','trait','tot_nsnps','tot_wtsum','mis_nsnps','mis_wtsum','prop_mis_nsnps','prop_mis_wtsum')

# prs example ----
# load in pgs data (from the pgs catalog)
pgs <- read.table('/PGS0000_trait.txt.annotated', header=T)

# # NB if missing rsid, get pgs info from .txt file, not .annotated
# pgs <- read.table('/PGS000_trait.txt', header=T, skip = 11, fill=T, na.strings = "", sep='\t')
# pgs <- pgs %>%
#   mutate(pos_key = case_when(is.na(rsID) ~ paste(chr_position, effect_allele, sep=':'),
#                              TRUE ~ paste(chr_name, chr_position, sep=':'))) %>%
#   select(pos_key, effect_allele, effect_weight)

# create summary (no snps and weightsum)
pgs_sum <- pgs %>%
  summarise(tot_wtsum = sum(abs(effect_weight)),
            tot_nsnps = length(pos_key)) %>%
  mutate(trait = 'trait')

# load in the calculated nopred data. 
# Add in trait id (not necessary, just for good measure)
# Filter out nopred comments that aren't NOSNP and NOALLELE
npr <- read.table('/trait.nopred', header=F, col.names=c('nopred','pos_key'))
npr <- npr %>%
  mutate(trait = 'trait') %>%
  filter(nopred == 'NOSNP' | nopred == 'NOALLELE')

trait <- left_join(npr,pgs, by='pos_key')

# check leftover (should be 0)
anti_trait <- trait %>%
  filter(is.na(effect_weight))
rm(anti_t2d)

# calculate missing wt sum and missing
trait_sum <- trait %>%
  summarise(mis_wtsum = sum(abs(effect_weight)),
            mis_nsnps = length(pos_key)) %>%
  mutate(trait = 'trait',
         pgsid = 'pgsid')

trait_sum <- full_join(trait_sum, pgs_trait_sum, by='trait')
trait_sum <- trait_sum %>%
  mutate(prop_mis_nosnp = mis_nsnps/tot_nsnps,
         prop_mis_wtsum = mis_wtsum/tot_wtsum)

qualmetrics <- rbind(trait_sum, qualmetrics)
rm(pgs,npr,trait,pgs_sum,triat_sum)


# grs's ----
# grs's don't (at present) have a log file/nopred output
# they do have a weight matrix
# if snps are missing, this shows up in an error message
# missing weights requires manually putting these in
grs <- read.table('/grs_weighted.txt', header=T)
grs <- grs %>%
  mutate(chrposall = paste(POSID,ALLELE,sep=':'))
grs_sum <- grs %>% 
  summarise(tot_nsnps = length(RSID),
         tot_wtsum = sum(BETA))
chrposall <- c('1:161643560:T','6:32585055:A') # put missing snps here
grs_missing <- as.data.frame(chrposall)
grs_missing <- left_join(grs_missing,grs, by='chrposall')
grs_missing_sum <- grs_missing %>%
  summarise(mis_nsnps = length(RSID),
            mis_wtsum = sum(BETA))
grs_sum <- cbind(grs_sum,grs_missing_sum)
grs_sum <- grs_sum %>%
  mutate(prop_mis_nosnp = mis_nsnps/tot_nsnps,
         prop_mis_wtsum = mis_wtsum/tot_wtsum,
         trait = 'trait',
         pgsid = 'grsid')
rm(grs,chrposall,grs_missing,grs_missing_sum)

qualmetrics <- rbind(qualmetrics,grs_sum)


# QC flags ----
qualmetrics <- qualmetrics %>%
  mutate(flag = case_when(prop_mis_nosnp  >= 0.1 & prop_mis_wtsum < 0.1 ~ 'high_missing_snps',
                          prop_mis_nosnp  < 0.1 & prop_mis_wtsum >= 0.1 ~ 'high_missing_wtsm',
                          prop_mis_nosnp  >= 0.1 & prop_mis_wtsum >= 0.1 ~ 'high_missing_snps_wtsum',
                          TRUE ~ 'no_flag'))

write.table(qualmetrics, 'qualmetrics.txt', row.names = F)
