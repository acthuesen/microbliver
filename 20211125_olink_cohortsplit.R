# WHAT THIS DOES ----
# - Split o-link data by cohort
# - convert to wide format
# - correct CBMRID
# - Add in timepoint/sample variables

# Note that the o-link qc has already been done before this
# Make sure you're working on the "Global QC" version of the data

# DEPENDENCIES ----
library(tidyverse)
library(readxl)

# GET DATA ----
# here's data 
olinkdat <- read.csv('Q:/Projects/MicrobLiver/Results/Olink/Global_QC/QC_olink.csv')

# CLEAN UP COHORTS ----
# HOL-F, HOL-N, Control, and CONTROL are all "artefacts" from QC. Remove them.
# Bridge are the samples used for the bridging procedure. Remove them. 
olinkdat <- olinkdat %>%
  filter(cohort != 'Control',
         cohort != 'CONTROL',
         cohort != 'HOL-N',
         cohort != 'HOL-F', 
         cohort != 'Bridge')

# TIEP and TIPSKC/ are spelling errors. Correct them. 
olinkdat <- olinkdat %>%
  mutate(cohort = case_when(cohort == 'TIEP' ~ 'TIPS',
                            cohort == 'TIPSKC/' ~ 'TIPS',
                            TRUE ~ cohort))

# IDs ---- 
# atm, id is a concatenated string of cohort and id (and sometimes sample info)
# cohort already exists as a variable, so take that out
# in the end we want a "CBMRID" which is cohort and padded id -sample -timepoint

# following matches non-capturing group of letters at beginning of string
# gsub with nothing to make newid
olinkdat <- olinkdat %>%
  mutate(newid = gsub("^(?:[A-Z]+)", "", SampleID))

# for HOL samples, reverse the newid order, since they're reversed
olinkdat <- olinkdat %>%
  mutate(newid = case_when(cohort == 'HOL' ~ str_replace(newid,"([^-]+)-([^-]+)","\\2-\\1"),
                           TRUE ~ newid))
  
# separate newid to get sample/time info
olinkdat <- olinkdat %>%
  separate(newid, into=c('newid','temp_sampleinfo'), extra='merge', sep='-')

# and because TIPS isn't separated by a delimiter, extract with regex for TIPS
olinkdat <- olinkdat %>%
  mutate(temp_sampleinfo = case_when(cohort == 'TIPS' ~ gsub("[0-9]+", "", newid),
                                     TRUE ~ temp_sampleinfo),
         newid = case_when(cohort == 'TIPS' ~ gsub('[A-Z]+', '', newid),
                           TRUE ~ newid))

# and because HOL has some duplicated samples, drop those with the suffix _b
# then drop the suffix 
olinkdat <- olinkdat %>% 
  filter(!(str_detect(newid, '_b') & cohort == 'HOL')) %>%
  mutate(newid = case_when(cohort == 'HOL' & str_detect(newid, '_a') ~ str_replace(newid, '_a', ''),
                           TRUE ~ newid))

# pad the id str so that it's 3 digits (4 for ALD)
olinkdat <- olinkdat %>%
  mutate(newid = case_when(cohort == 'ALD' ~ str_pad(newid, width=4, side='left', '0'),
                           cohort != 'ALD' ~ str_pad(newid, width=3, side='left', '0')))

# and finally concat the newid and cohort
olinkdat <- olinkdat %>%
  unite('newid', c(cohort, newid), remove=FALSE, sep='')

# SAMPLE TYPE AND TIMEPOINT ----
# As the only cohort, TIPS only has sampletype 
# As the only cohort, ALCO has BOTH timepoints and sampletype
# The following have only timepoints: HCOT, HOL, MLGB, PRF, ROC, RFX
# The following have no sample variables: ALD, HP
olinkdat <- olinkdat %>%
  mutate(sampletype = NA_character_,
         timepoint = NA_character_) %>%
  mutate(sampletype = case_when(cohort == 'TIPS' ~ temp_sampleinfo,
                                TRUE ~ sampletype),
         sampletype = case_when(cohort == 'ALCO' ~ str_extract(temp_sampleinfo, "S[0-9]+"),
                                TRUE ~ sampletype)) %>%
  mutate(timepoint = case_when(cohort == 'ALCO' ~ str_extract(temp_sampleinfo, 'T[0-9]+'),
                               TRUE ~ timepoint),
         timepoint = case_when(cohort == 'HCOT' | 
                                 cohort == 'HOL' | 
                                 cohort == 'MLGB' | 
                                 cohort == 'PRF' | 
                                 cohort == 'RDC' | 
                                 cohort == 'RFX' ~ temp_sampleinfo,
                               TRUE ~ timepoint)) %>%
  select(-temp_sampleinfo)

# cleanup timepoint notation for HCOT samples
# add in baseline "timepoint" for HOL where visit number refers to visit after baseline
olinkdat <- olinkdat %>%
  mutate(timepoint = case_when(cohort == 'HCOT' & timepoint == 'A' ~ 'BL',
                               cohort == 'HCOT' & timepoint == 'B' ~ 'OF',
                               cohort == 'HCOT' & timepoint == 'C' ~ 'TR',
                               cohort == 'HOL' & is.na(timepoint) ~ 'BL',
                               TRUE ~ timepoint))

# CBMR ID ----
olinkdat <- olinkdat %>%
  unite('CBMRID', c(newid,sampletype,timepoint), sep='-', na.rm=T, remove = F)

# WIDE FORMAT ----
olinkdat_wide <- olinkdat %>%
  select(-c(OlinkID, UniProt, MissingFreq, Panel, Panel_Version, PlateID, QC_Warning, LOD, NPX, Normalization, Index, no_missing_obs, total_obs, propmissing, keep_assay)) %>%
  pivot_wider(names_from = Assay, values_from = corrected_NPX)

# SPLIT BY COHORT ----
 setwd('Q:/Projects/MicrobLiver/Results/Olink/Global_QC/split_cohorts/') # this is where your output goes

olinkdat_wide %>%
  group_split(cohort) %>%
  sapply(., function (x) 
    write.csv(x, file=paste0("QC_Olink", unique(x$cohort), "_wide.csv"), row.names=F))

olinkdat %>%
  group_split(cohort) %>%
  sapply(., function (x) 
    write.csv(x, file=paste0("QC_Olink", unique(x$cohort), "_long.csv"), row.names=F))
