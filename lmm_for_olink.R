# --------------------------------------------------------------------------- #
# First version: 24 May 2022
# Latest revision: 30 May 2022
# 
# Analysis of Olink data with arbitrary number of parallel groups/timepoints. 
# 
# Data should be long format and contain: 
#   1. A grouping/intervention variable ('group')
#   2. A time variable coded as integer and factor ('visit_int'+'visit_fac')
#   3. An ID variable
# --------------------------------------------------------------------------- #

# -------------------------- #
# Load data and dependencies #
# -------------------------- #
library(tidyverse)
library(nlme)
library(broom.mixed)

setwd('wdir')

df <- read.table('data', header=T)

# --------------------------- #
# Fit LMM: gls (nlme) version #
# --------------------------- # 

# Nest assays
# Tidy model output -> unnest
# Adjust for multiple comparisons

lmm <- df %>%
  group_by(Assay) %>%
  nest() %>%
  mutate(lmm = map(data, ~ gls(
    int_NPX~visit_fac * group,
    data=.x,
    correlation=corSymm(form=~visit_int|newid),
    weights=varIdent(form=~1|visit_int))),
    tidied = map(lmm, tidy)) %>%
  unnest(tidied) %>%
  ungroup() %>%
  select(-c(data, lmm, statistic)) %>%
  filter(str_detect(term, 'group|visit')) %>%
  group_by(term) %>%
  mutate(adjp = p.adjust(p.value, method='fdr'))


# ------------- #
# Volcano plots #
# ------------- #

# Extract factor of interest from the term
lmm_forjoin <- lmm %>%
  filter(!str_detect(term, 'group')) %>%
  mutate(visit_fac = gsub('visit_fac',"",term),
         visit_fac = as.factor(visit_fac))

# Calculate changes compared to the reference
df_wfc <- df %>%
  filter(visit_fac != '1') %>%
  inner_join(df %>% 
               filter(visit_fac =='1'), by=c('newid','Assay'), suffix=c('','1')) %>%
  select(-c(group1, visit_int1,visit_fac1)) %>%
  mutate(fc = NPX-NPX1) %>%
  group_by(Assay,visit_fac) %>%
  summarise(mean_fc = mean(fc, na.rm=T)) 

# Join your model output with your data
df_wfc_wlmm <- df_wfc %>%
  left_join(., lmm_forjoin, by=c('Assay','visit_fac')) %>%
  mutate(nom_sigupdown = case_when(mean_fc<0 & p.value<0.05 ~ 'sigdown',
                                   mean_fc>0 & p.value<0.05  ~ 'sigup',
                                   TRUE ~ 'insig')) %>%
  mutate(act_sigupdown = case_when(mean_fc<0 & adjp<0.05 ~ 'sigdown',
                                   mean_fc>0 & adjp<0.05 ~ 'sigup',
                                   TRUE ~ 'insig'))
rm(lmm_forjoin,df_wfc)

#Pick some good colors for fill and color
vpc <- c('#8D3344','#1d738b','#cccccc')
names(vpc) <- c('sigup','sigdown','insig')
vpf <- c('#8D3344','#1d738b','#ffffff')
names(vpf) <- c('sigup','sigdown','insig')

# Create your volcano plot
df_wfc_wlmm %>% {
    ggplot(., aes(x=mean_fc, 
                  y=-log10(p.value), 
                  group=Assay, 
                  color=nom_sigupdown,
                  fill=act_sigupdown)) + 
      geom_vline(xintercept = 0, linetype='dashed', color='gray50') +
      geom_point(shape=21) + 
      geom_text(data=filter(., p.value<0.05),
                aes(label=Assay),
                nudge_y= 0.1, size=2) + 
      scale_color_manual(values = vpc) +
      scale_fill_manual(values = vpf) +
      labs(title='THIS IS A GOOD TITLE FOR A PLOT',
           subtitle = 'Fold change of cytokine level at each visit',
           x='log2(Fold Change)',
           y='-log10(p-value)') +
      facet_wrap(~visit_fac, scales='free') +
      theme_bw() + 
      theme(
        legend.position = 'none',
        strip.background = element_rect(fill='white')
      )
}  
