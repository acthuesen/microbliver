# DEPENDENCIES AND WD ----
library(readxl)
library(tidyverse)
setwd('Q:/Projects/MicrobLiver/Results/Genetics/SNPs_PRSs_for_datahub/') # where genetic data are and go

# ID MASTER FILE ----
# Some ID's belong to multiple cohorts. These affiliations are encoded in the ID MASTER FILE
idcompare <- read_excel('preprocessing/IDcompare_ALD_SIP_ALCO_RFX_PRF_HP.xlsx', sheet='clean', na="NA")

# Four ID's from ALCO cohort are duplicated, remove from the master
idcompare <- idcompare %>%
  mutate(Alcochallenge_id = na_if(Alcochallenge_id, 211),
         Alcochallenge_id = na_if(Alcochallenge_id, 214),
         Alcochallenge_id = na_if(Alcochallenge_id, 215), 
         Alcochallenge_id = na_if(Alcochallenge_id, 302))


# SINGLE SNPS / SPLIT BY COHORT ----
# load snp data
genedata <- read.csv('preprocessing/snps/microbliver_snps.geno.csv', check.names = F) #check names ensures that chr isn't added to numeric vars (i.e. optional)

# fix id col which doesn't have a varname with Carsten's SNP extractor
colnames(genedata)[1] <- 'id' 

# drop numeric cohort identifiers from the ID name
# also drop the pnpla3 snp because it's added from the previous imputation batch (see below)
# also drop the id col again, because we create a better id
genedata <- genedata %>%
  mutate(CBMR_ID = str_replace_all(id, "[0-9]+x", "")) %>%
  select(-`22:44324727:C:G`) %>%
  select(-id)

# merge ids with data
genedata <- left_join(genedata, idcompare, by='CBMR_ID')

# get genotyped pnpla3
pnpla3 <- read.csv('preprocessing/snps/pnpla3snp_1kg.geno.csv', header=T, check.names = F)

# fix id as above and then drop
# correct spelling mistake in cohort which only exists in this dataset
colnames(pnpla3)[1] <- 'id' 
pnpla3 <- pnpla3 %>%
  mutate(CBMR_ID = str_replace_all(id, "[0-9]+x", "")) %>%
  mutate(CBMR_ID = gsub('AlCO', 'ALCO', CBMR_ID)) %>%
  select(-id)

# merge with previous data
genedata <- full_join(genedata,pnpla3, by='CBMR_ID')

# unite cohorts into one "multiple cohorts" column and then separate into multiple rows
# remove the sample from SIP that was supposed to be removed during imputation (SIPremove)
genedata <- genedata %>%
  mutate(cohort1 = str_replace_all(CBMR_ID, "[0-9]+", "")) %>%
  mutate(cohort2 = case_when(is.na(Rifaximin_id)==FALSE ~ "RFX")) %>%
  mutate(cohort3 = case_when(is.na(Alcochallenge_id)==FALSE ~ "ALCO")) %>%
  mutate(cohort4 = case_when(is.na(Profermin_id)==FALSE ~ "PRF")) %>%
  unite("cohort", cohort1:cohort4, na.rm=T, remove=FALSE) %>%
  select(-Rifaximin_id, -Profermin_id, -Alcochallenge_id, -(cohort1:cohort4)) %>%
  separate_rows(cohort, sep="_") %>%
  filter(cohort != 'SIPremove')

# split by cohort
genedata %>% 
  group_split(cohort) %>% 
    sapply(., function (x) 
    write.csv(x, file=paste0("preprocessing/split/extractedsnps_", unique(x$cohort), ".csv"), row.names=F))


# PRS/GRS / SPLIT BY COHORT ----
# Unlike SNP data where multiple SNPs are in one initial file, each PRS/GRS is its own file
# therefore, load multiple files and merge into one before correcting ids etc and splitting
# load all prs files
prss <- list.files("preprocessing/prs", pattern=".prs", ignore.case=T)

# init a new dataframe (so they can merge)
allprs <- data.frame(IID=character())
# merge all prs files 
for(prs in prss){
  tempdata <- read.table(paste0("preprocessing/prs/", unique(prs)), header=T) #load
  allprs <- full_join(allprs, tempdata, by="IID") # merge
  rm(tempdata) # cleanup
}

# load all grs files
grss <- list.files("preprocessing/prs", pattern=".txt", ignore.case=T)

# init a new dataframe (so they can merge)
allgrs <- data.frame(IID=character())
# merge all prs files 
for(grs in grss){
  tempdata <- read.table(paste0("preprocessing/prs/", unique(grs)), header=T, col.names=c('IID',grs)) #load
  allgrs <- full_join(allgrs, tempdata, by="IID") # merge
  rm(tempdata) # cleanup
}

# join prs's and grs's
allrs <- full_join(allprs, allgrs, by='IID')

# cleanup id
allrs <- allrs %>%
  mutate(CBMR_ID = str_replace_all(IID, "[0-9]+x", ""))

# merge ids with data
allrs <- left_join(allrs, idcompare, by='CBMR_ID')

# unite cohorts and separate
# also cleanup column names
allrs <- allrs %>%
  mutate(cohort1 = str_replace_all(CBMR_ID, "[0-9]+", "")) %>%
  mutate(cohort2 = case_when(is.na(Rifaximin_id)==FALSE ~ "RFX")) %>%
  mutate(cohort3 = case_when(is.na(Alcochallenge_id)==FALSE ~ "ALCO")) %>%
  mutate(cohort4 = case_when(is.na(Profermin_id)==FALSE ~ "PRF")) %>%
  unite("cohort", cohort1:cohort4, na.rm=T, remove=FALSE) %>%
  select(-Rifaximin_id, -Profermin_id, -Alcochallenge_id, -IID, -(cohort1:cohort4)) %>%
  separate_rows(cohort, sep="_") %>%
  filter(cohort != 'SIPremove') %>%
  rename_with(~str_remove(., '.txt')) %>%
  rename_with(~str_remove(., '_MBL_correctX'))

# SPLIT BY COHORT ----
allrs %>% 
  group_split(cohort) %>% 
  sapply(., function (x) 
    write.csv(x, file=paste0("preprocessing/split/riskscores_", unique(x$cohort), ".csv"), row.names=F))
