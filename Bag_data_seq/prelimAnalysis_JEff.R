library(tidyverse)

dat <- read_tsv('Raw_data/Jeff_Prelim/SMM/ITS_output_clean.tsv')

# looking at guilds
dat %>% 
  mutate(guild2 = case_when(trophicMode == 'Saprotroph' ~ 'sap', 
                            guild == 'Arbuscular Mycorrhizal' ~ 'amf', 
                            guild == 'Ectomycorrhizal' ~ 'ecm', 
                            TRUE ~ 'other_unknown')) %>% 
  group_by(barcode, guild2) %>% 
  summarise(count=sum(count)) %>% 
  ungroup() %>% 
  pivot_wider(names_from=guild2, values_from=count, values_fill=0) %>% 
  mutate(lrr_ecm_to_sap = log((ecm+0.5)/(sap+0.5)), 
         lrr_ecm_to_amf = log((ecm+0.5)/(amf+0.5)), 
         ra_sap = sap/(amf+ecm+other_unknown+sap), 
         ra_ecm = ecm/(amf+ecm+other_unknown+sap), 
         ra_amf = amf/(amf+ecm+other_unknown+sap), 
         # following creates rough categories of mycorrhizal fungal relative abundance
         # tries to account for low affinity of primers to AM fungal DNA
         # very arbitrary, interpret with caution
         lev_sap = case_when(ra_sap < 0.1 ~ 'vlow', 
                             ra_sap >= 0.1 & ra_sap < 0.3 ~ 'low', 
                             ra_sap >= 0.3 & ra_sap < 0.6 ~ 'med', 
                             ra_sap >= 0.6 ~ 'high'), 
         lev_ecm = case_when(ra_ecm < 0.1 ~ 'vlow', 
                             ra_ecm >= 0.1 & ra_ecm < 0.3 ~ 'low', 
                             ra_ecm >= 0.3 & ra_ecm < 0.6 ~ 'med', 
                             ra_ecm >= 0.6 ~ 'high'), 
         lev_amf = case_when(ra_amf < 0.0001 ~ 'vlow', 
                             ra_amf >= 0.0001 & ra_amf < 0.01 ~ 'low', 
                             ra_amf >= 0.01 & ra_amf < 0.03 ~ 'med', 
                             ra_amf >= 0.03 ~ 'high'), 
         lev_mf = case_when(lev_ecm == 'high' | lev_amf == 'high' ~ 'high', 
                            lev_ecm == 'med' | lev_amf == 'med' ~ 'med', 
                            lev_ecm == 'low' | lev_amf == 'low' ~ 'low', 
                            lev_ecm == 'vlow' & lev_amf == 'vlow' ~ 'vlow'), 
         lev_sap = factor(lev_sap, levels=c('vlow', 'low', 'med', 'high')), 
         lev_amf = factor(lev_amf, levels=c('vlow', 'low', 'med', 'high')), 
         lev_ecm = factor(lev_ecm, levels=c('vlow', 'low', 'med', 'high')), 
         lev_mf = factor(lev_mf, levels=c('vlow', 'low', 'med', 'high'))) %>% 
  left_join(dat %>% 
              group_by(barcode) %>% 
              summarise(mass_extracted = first(mass_extracted))) -> tmp

# rough grouping of samples based on relative abundance of sequence reads for guilds
# SAP fungi
table(tmp$lev_sap)
# AM fungi
table(tmp$lev_amf)
# EcM fungi
table(tmp$lev_ecm)
# AM or EcM fungi
table(tmp$lev_mf)

# relative abundance of SAP fungal reads
ggplot(tmp, aes(x=fct_reorder(barcode, ra_sap), y=ra_sap)) + 
  geom_bar(stat='identity') + 
  theme(axis.text.x=element_text(angle=90, size=5))

# relative abundance of EcM fungal reads
ggplot(tmp, aes(x=fct_reorder(barcode, ra_ecm), y=ra_ecm)) + 
  geom_bar(stat='identity') + 
  theme(axis.text.x=element_text(angle=90, size=5))

# relative abundance of AM fungal reads
ggplot(tmp, aes(x=fct_reorder(barcode, ra_amf), y=ra_amf)) + 
  geom_bar(stat='identity') + 
  theme(axis.text.x=element_text(angle=90, size=5))

# log response ratios of EcM fungal reads to saprotroph reads
ggplot(tmp, aes(x=fct_reorder(barcode, lrr_ecm_to_sap), y=lrr_ecm_to_sap)) + 
  geom_bar(stat='identity') + 
  theme(axis.text.x=element_text(angle=90, size=5))

# log response ratios of EcM fungal reads to AM fungal reads
ggplot(tmp, aes(x=fct_reorder(barcode, lrr_ecm_to_amf), y=lrr_ecm_to_amf)) + 
  geom_bar(stat='identity') + 
  theme(axis.text.x=element_text(angle=90, size=5))

# amount of biomass extracted for different levels of mycorrhizal fungal reads
ggplot(tmp, aes(x=lev_mf, y=mass_extracted)) + 
  geom_boxplot()

# amount of biomass extracted for different levels of saprotrophic fungal reads
ggplot(tmp, aes(x=lev_sap, y=mass_extracted)) + 
  geom_boxplot()

