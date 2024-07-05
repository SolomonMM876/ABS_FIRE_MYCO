install.packages("devtools")
devtools::install_github("ropenscilabs/datastorr")
devtools::install_github("traitecoevo/fungaltraits")
library(fungaltraits)
fungal_traits<-fungal_traits()

na.omit(fungal_traits$em_expl)


left_join(out%>% rename(species=Species),fungal_traits%>%select(species,Genus,em_expl,em_text), by= 'species')
out$Genus
fungal_traits$Genus
