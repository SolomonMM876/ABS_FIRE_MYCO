library(readr)
library(tidyverse)
library(Hmsc)
library(ape) # we need this to construct a taxonomic tree
library(dplyr)


#set wd
localDir = "HMSC_ABS"
dataDir = file.path(localDir, "data")
modelDir = file.path(localDir,"models")

data = read.csv(file.path(dataDir, "Bag_data.csv"))

PhyData=read.csv(file.path(dataDir,'Trait_Phylo_data.csv'))%>%
  mutate(across(everything(), as.factor))

#OTUs with both genome size and explo type
PhyData_filtered<-PhyData%>%filter(!is.na(exploration_type))%>%filter(!is.na(mean_gs))

PhyData<-PhyData%>%
  select(-exploration_type,-Ecm_lineage,-mean_gs)

TrData=read.csv(file.path(dataDir,'Trait_Phylo_data.csv'))%>%
  mutate(across(everything(), as.factor))%>%
  select(OTU:last_col())

#OTUs with both genome size and explo type
TrData_filtered<-TrData%>%filter(!is.na(exploration_type))%>%filter(!is.na(mean_gs))




#order TrData and Y data the same
otu_order <- PhyData_filtered$OTU

# READ AND MODIFY ENVIRONMENTAL DATA (BEGINNING)####
#Transform variables to more normal distributions and set factors
XData <-data%>%
  select(Fire.Interval,Fire.Severity,Ortho_P_mg_kg,Nitrate_mg_kg,Ammonia_mg_kg,pH,readcount,
         log10_biomass_day,perc_myco_host_freq)%>%
  mutate(Fire.Interval=as.factor(Fire.Interval),
         Fire.Severity=as.factor(Fire.Severity),
         Ortho_P_mg_kg=log10(Ortho_P_mg_kg),
         Nitrate_mg_kg=log10(Nitrate_mg_kg),
         Ammonia_mg_kg=log10(Ammonia_mg_kg),
         pH=log10(pH)
  )
#set factor order
XData$Fire.Severity <- factor(XData$Fire.Severity, levels = c("Low", "High"))
XData$Fire.Interval <- factor(XData$Fire.Interval, levels = c("Long", "Short"))
##########

#Select Community Data
YData<-data%>%
  select(starts_with('ITSall'))%>%
  select(all_of(otu_order))


# I could reduce the number of species here depending if I want to focus on more or less common
#but considering how few IDs I actually have I will leave it as is
# I would use something with sel.sp = colSums(YData>0)>=10 and then filter for only those cols...

P = colMeans(YData>0, na.rm = TRUE)
A = colSums(YData, na.rm = TRUE)/sum(YData, na.rm = TRUE)

par(mfrow=c(1,2))
hist(P,xlim=c(0,1),breaks = seq(from=0,to=1,by=0.1), col = "grey", xlab = "Prevalence")
hist(log(A+0.001,base=10),breaks=10, col = "grey", xlab = "Abundance")

#look at relative abundances of species
library(ggplot2)
data.frame(
  Species = names(YData),
  Relative_Abundance = colSums(YData, na.rm = TRUE)/sum(YData, na.rm = TRUE))%>%
  ggplot( aes(x = reorder(Species, Relative_Abundance), y = Relative_Abundance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip to make labels more readable
  labs(title = "Relative Abundance of Fungal Species",
       x = "Species",
       y = "Relative Abundance") +
  theme_minimal() 

# with these data the simplest starting point
# would be to fit a presence-absence model with decay class and log-transformed read-count as explanatory variables,
# and the random effect of the log as a random effect. This model can be defined as follows:

Y = as.matrix(YData)

# For this model We truncate the data to presence-absence
#Y = Y>0 #leaving full species occurrences 


# convert Y to numerical either as Y = 1*Y or as

Y = apply(Y,MARGIN = 2,FUN = as.numeric)

##################################################################################################
# SET UP THE MODEL (BEGINNING)
##################################################################################################
# STUDY DESIGN

studyDesign = data.frame(
  Location = data$Location,
  Transect = data$Transect,
  Site = data$Site
) %>%
  mutate(across(everything(), as.factor))


rL_location = HmscRandomLevel(units = studyDesign$Location)
rL_transect = HmscRandomLevel(units = studyDesign$Transect)
rL_site = HmscRandomLevel(units = studyDesign$Site)
#rL_spatial= HmscRandomLevel(sData=xy, longlat=TRUE)


# REGRESSION MODEL FOR ENVIRONMENTAL COVARIATES.
XFormula = ~Fire.Severity + Fire.Interval + log(readcount)+ 
  log10_biomass_day+ perc_myco_host_freq +
  Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg+ pH #Soil Properties

# CONSTRUCT TAXONOMICAL TREE TO BE USED AS PROXY FOR PHYLOGENETIC TREE
taxonomicTree <- as.phylo(~phylum/class/order/family/genus/OTU, data = PhyData_filtered, collapse = FALSE)
taxonomicTree$edge.length = rep(1,length(taxonomicTree$edge))
plot(taxonomicTree,cex=0.5) 


#Trait matrix, which I dont have because I dont have species traits here
TrFormula = ~exploration_type +  Ecm_lineage + mean_gs

#make sure TrData and Y matrix have some species in same order

# Get column names from Y and row names from Trdata
Y_cols <- colnames(Y)

TrData<-TrData_filtered%>%
  mutate(OTU = factor(OTU, levels = Y_cols)) %>%
  arrange(OTU)%>%
  column_to_rownames('OTU')

Y = as.matrix(YData)
Y = 1*(Y>0)
#build model
m = Hmsc(Y = Y, XData = XData, XFormula = XFormula, 
         phyloTree = taxonomicTree,
         TrData = TrData,TrFormula = TrFormula ,
         studyDesign = studyDesign, 
         ranLevels = list(Location = rL_location,
                          Transect = rL_transect, 
                          Site = rL_site),
         distr="probit")#poisson

#Above I have set up the model

#Below I will double check the structure and design
######################################################

# To understand what was done with the above script, let us look at the objects:
# The study design simply lists the id of each log.
head(studyDesign)

#The right-hand side (rL) describes the structure of 
# random effect. To see the structure, evaluate
rL_site

# It is always a good idea to look at how the XData and XFormula are translate to the design matrix X. 
# These can be seen as e.g. by evaluating the following

m$XFormula
head(m$XData)
head(m$X)
head(m$TrData)

# Note further that there is "1 trait" even if we did not define any traits. To see why this is the case, evaluate
head(m$Tr)

#good idea to look at the model object as well:
m

# COMBINING AND SAVING MODELS (START)
##################################################################################################
models = list(m)
names(models) = c("presence-absence model_w_Traits")
save(models, file = file.path(modelDir, "unfitted_models.RData"))
# TESTING THAT MODELS FIT WITHOUT ERRORS (START)
##################################################################################################
for(i in 1:length(models)){
  print(i)
  sampleMcmc(models[[i]],samples=2)
}
##################################################################################################
# TESTING THAT MODELS FIT WITHOUT ERRORS (END)


load(file=file.path(modelDir,"unfitted_models.RData"))

nm = length(models)
samples_list = c(250,250,250)
thin_list = c(250,500,1000)
nChains = 2
nParallel = nChains
Lst = 1
while(Lst <= length(samples_list)){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  print(paste0("thin = ",as.character(thin),"; samples = ",as.character(samples)))
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){
    print("model had been fitted already")
  } else {
    print(date())
    for (mi in 1:nm) {
      print(paste0("model = ",names(models)[mi]))
      m = models[[mi]]
      m = sampleMcmc(m, samples = samples, thin=thin,
                     adaptNf=rep(ceiling(0.4*samples*thin),m$nr),
                     transient = ceiling(0.5*samples*thin),
                     nChains = nChains,
                     nParallel = nParallel)
      models[[mi]] = m
    }
    save(models,file=filename)
  }
  Lst = Lst + 1
}
#this is if we want to test a set thin, whereas below is if we want to let it run for a while


# nChains = 2
# samples = 10
# thin = 1 # try with thin = 1, thin = 10, thin = 100, etc.
# transient<-0.3*thin*samples
# nParallel = nChains
# 
# model<-sampleMcmc(m, samples = samples, thin=thin,
#                   adaptNf=rep(ceiling(0.4*samples*thin),m$nr),
#                   transient = ceiling(0.5*samples*thin),
#                   nChains = nChains,
#                   nParallel = nParallel)
# filename=file.path(modelDir, paste0("models_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
# save(model,file=filename)


##################################################################################################

# hurdle model
# models = list()
# for (i in 1){
#   Y = as.matrix(YData)
#   {Y = 1*(Y>0)}
#   XFormula =     ~~Fire.Severity + Fire.Interval + log(readcount)+ 
#     log10_biomass_day+ perc_myco_host_freq +
#     Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg+ pH
#   
#   m = Hmsc(Y = Y, XData = XData, XFormula = XFormula,
#            phyloTree = taxonomicTree,
#            TrData = TrData,TrFormula = TrFormula ,
#            studyDesign = studyDesign, ranLevels = list(Location = rL_location,
#                                                        Transect = rL_transect,
#                                                        Site = rL_site),
#            distr=switch(i,"probit"),
#            YScale = TRUE)
#   models[[i]] = m
# }
# names(models) = c("probit")

#All of this that is below is if we want to run different iterations of the model in a loop 
##################################################################
#################################################################

# In all models, we also add a random effect at the sampling unit level. The random effect models associations among the species,
# which is what we are primarily interested about.
#i=1 is a lognormal Poisson mode
#i=2 and i=3 form together a hurdle model that separates presence-absence variation from abundance variation
#i=2 is a probit model that is fitted to sequences countstruncated to presence-absence
#i=3 is a normal model that is fitted to log-transformed sequence counts conditional on presence. 
#This means that for i=3, sampling units where the species is not present are shown in the Y matrix as missing data (NA) rather than as zeros

# this script is if I want to test different distributions
# models = list()
# for (i in 1:3){
#   Y = as.matrix(YData)
#   if (i==2) {Y = 1*(Y>0)}
#   if (i==3) {
#     Y[Y==0] = NA
#     Y = log(Y)
#   }
#       XFormula =     ~~Fire.Severity + Fire.Interval + log(readcount)+ 
#                         log10_biomass_day+ perc_myco_host_freq +
#                         Ortho_P_mg_kg+ Nitrate_mg_kg+ Ammonia_mg_kg+ pH
#       
#     m = Hmsc(Y = Y, XData = XData, XFormula = XFormula,
#              phyloTree = taxonomicTree,
#              TrData = TrData,TrFormula = TrFormula ,
#              studyDesign = studyDesign, ranLevels = list(Location = rL_location,
#                                                          Transect = rL_transect,
#                                                          Site = rL_site),
#              distr=switch(i,"lognormal poisson","probit","normal"),
#              YScale = TRUE)
#   models[[i]] = m
# }
# names(models) = c("lognormal_poisson","probit","normal")


# In the above script, we have used the option Yscale = TRUE to scale the response data to zero mean and unit variance.
# is discussed in more detail in Section 8.3 of the book, this scaling influences only the normal model,
# and it is done to make the default priors of Hmsc compatible with the data.


#############################################################
# # Thus, in summary, running the model fitting for thin = 1, 10, 100, 1000 typically saves a lot of time,
# 
# samples = 250
# nChains = 3
# 
# for (thin in c(250,500,1000)){
#   transient = 0.3*thin*samples
#   for (i in 1:3){
#     cat("Running model =", i, "with thin =", thin, "\n")
#     models[[i]] = sampleMcmc(models[[i]], thin = thin, samples = samples, transient = transient,
#                                     nChains = nChains, nParallel = nChains )
#   }
#   filename=file.path(modelDir, paste0("models_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
#   save(models,file=filename)
# }
# 

