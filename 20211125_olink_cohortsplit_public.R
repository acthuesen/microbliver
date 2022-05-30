# Updated 05 May 2022
# WHAT THIS DOES ----
# - Split o-link data by cohort
# - convert to wide format
# - correct CBMRID
# - Add in timepoint/sample variables

# Note that the o-link qc has already been done before this

# This script is updated for batch 3


# DEPENDENCIES ----
library(tidyverse)
library(readxl)

# GET DATA ----
# here's data 
olinkdat <- read.csv('PATH/Olink/Global_QC/QC_olink.csv')

# CLEAN UP COHORTS ----
# HOL-F, HOL-N, Control, and CONTROL are all "artefacts" from QC. Remove them.
# Bridge are the samples used for the bridging procedure. Remove them. 
olinkdat <- olinkdat %>%
  filter(cohort != 'Control',
         cohort != 'CONTROL',
         cohort != 'HOL-N',
         cohort != 'HOL-F', 
         cohort != 'Bridge')

# IDs ---- 
# atm, id is a concatenated string of cohort and id (and sometimes sample info)
# cohort already exists as a variable, so take that out
# in the end we want a "CBMRID" which is cohort and padded id -sample -timepoint

# following matches non-capturing group of letters at beginning of string
# gsub with nothing to make newid
olinkdat <- olinkdat %>%
  mutate(newid = gsub("^(?:[a-zA-Z]+)", "", SampleID))

## ADDED IN THIS VERSION: added removal of leading _ (PK and AHEP samples)
olinkdat <- olinkdat %>%
  mutate(newid = gsub("^_", "", newid))

# for HOL samples, reverse the newid order, since they're reversed
olinkdat <- olinkdat %>%
  mutate(newid = case_when(cohort == 'HOL' ~ str_replace(newid,"([^-]+)-([^-]+)","\\2-\\1"),
                           TRUE ~ newid))

# separate newid to get sample/time info
## IN THIS VERSION: nb, this also separates the IID from the FID-IID in PK
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
                           cohort == 'SIP' ~ str_pad(newid, width=5, side='left', '0'),
                           cohort == 'PK' | cohort == 'AHEP' ~ newid,
                           cohort != 'ALD' | cohort != 'SIP' | cohort != 'PK' | cohort != 'AHEP' ~ 
                             str_pad(newid, width=3, side='left', '0')))

# and finally concat the newid and cohort
olinkdat <- olinkdat %>%
  unite('newid', c(cohort, newid), remove=FALSE, sep='')

# SAMPLE TYPE AND TIMEPOINT ----
# As the only cohort, TIPS only has sampletype 
# As the only cohort, ALCO has BOTH timepoints and sampletype
# The following have only timepoints: HCOT, HOL, MLGB, PRF, ROC, RFX
# The following have no sample variables: ALD, HP
## CHANGED IN THIS VERSION: AHEP, PK, ProD, SIP all also have no sample variables
## However, for both PK and AHEP something is delineated by a hyphen which needs to be added back in
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
  mutate(otherinfo = case_when(cohort == 'AHEP' | 
                                       cohort == 'PK' ~ temp_sampleinfo)) %>%
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
  unite('CBMRID', c(newid,sampletype,timepoint,otherinfo), sep='-', na.rm=T, remove = F) %>%
  select(-otherinfo)

# WIDE FORMAT ----
olinkdat_wide <- olinkdat %>%
  select(-c(OlinkID, UniProt, MissingFreq, Panel, Panel_Version, PlateID, QC_Warning, LOD, corrected_NPX, Normalization, Index, no_missing_obs, total_obs, MissingProp)) %>%
  pivot_wider(names_from = Assay, values_from = NPX)

# SPLIT BY COHORT ----
setwd('PATH/Olink/Global_QC/split_cohorts/') # this is where your output goes

olinkdat_wide %>%
  group_split(cohort) %>%
  sapply(., function (x) 
    write.table(x, file=paste0("QC_Olink", unique(x$cohort), "_wide.tsv"), row.names=F, sep='\t'))

olinkdat %>%
  group_split(cohort) %>%
  sapply(., function (x) 
    write.table(x, file=paste0("QC_Olink", unique(x$cohort), "_long.tsv"), row.names=F, sep='\t'))

