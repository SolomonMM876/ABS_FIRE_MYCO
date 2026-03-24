
### load libraries
library(Hmsc)


### set up workspace (only reference) ###

# get selected predictors from ordination
load('output/workspace/ordistep_reference.RData')

# response and predictor tables
Y <- otus %>% filter(preds$Treat == 'Reference') %>% 
  select_if(function(x)sum(x)>0) %>% 
  as.matrix()
X <- preds %>% filter(Treat == 'Reference') %>% 
  select(gsub('`', '', attr(formula(m1.ord), 'term.labels'))) %>% 
  mutate(across(where(is.character), as_factor))

# study design and coords
XY <- xy 

# remove NA rows
x <- complete.cases(X)
X <- X[x, ]
Y <- Y[x, ]
XY <- XY[x, ]

# random effect structure -- study design
rL1 <- HmscRandomLevel(units = XY$Plot_ID)
rL2 <- HmscRandomLevel(units = XY$Subplot)

# random effect structure -- spatial latent variable
rL0 <- HmscRandomLevel(sData=XY %>% 
                         group_by(Plot_ID) %>% 
                         slice(1) %>% 
                         column_to_rownames('Plot_ID') %>% 
                         select(X, Y), longlat=TRUE)

# fixed effects
XFormula <- ~ .

## prepare models
models <- list()

# presence-absence model
models[[1]] <- Hmsc(Y = 1*(Y>0), XData = X, XFormula = XFormula, 
                    studyDesign = XY %>% select(Plot_ID, Subplot) %>% mutate(across(everything(), as_factor)), 
                    ranLevels = list(Plot_ID=rL1, Subplot=rL2, Plot_ID=rL0),
                    distr="probit",
                    YScale = TRUE)


# abundance model
Y2 <- Y
Y2[Y2==0] <- NA
models[[2]] <- Hmsc(Y = log(Y2), XData = X, XFormula = XFormula, 
                    studyDesign = XY %>% select(Plot_ID, Subplot) %>% mutate(across(everything(), as_factor)), 
                    ranLevels = list(Plot_ID=rL1, Subplot=rL2, Plot_ID=rL0),
                    distr="normal",
                    YScale = TRUE)

# model parameters for testing
# thin <- 1
# samples <- 2
# nChains <- 2
# transient <- 1

# model parameters
thin <- 20
samples <- 500
nChains <- 2
transient <- 0.3*thin*samples

# export workspace
 save(models, thin, samples, nChains, transient, file='output/workspace/hmsc_reference/mcmc_input.RData')


## fit models on server

# ssh jeff@203.101.228.37
# cd /home/jeff/s01212ss-hie6-vm-space/RProjects/snowies

# start R
# /opt/R/4.0.3/bin/R

# load library
library(Hmsc)

# load workspace
load('output/workspace/hmsc_reference/mcmc_input.RData')

# fit models -- one in each screen

models1 <- sampleMcmc(models[[1]], samples = samples, thin = thin, transient = transient, nChains = nChains, 
                      nParallel = 2)
save(models1, file='output/workspace/hmsc_reference/mcmc_output1.RData')

models2 <- sampleMcmc(models[[2]], samples = samples, thin = thin, transient = transient, nChains = nChains, 
                      nParallel = 2)
save(models2, file='output/workspace/hmsc_reference/mcmc_output2.RData')


# move to one screen, load output and summarise posterior estimates

load('output/workspace/hmsc_reference/mcmc_output1.RData')
load('output/workspace/hmsc_reference/mcmc_output2.RData')

models[[1]] <- models1
models[[2]] <- models2
rm(models1, models2)

## check model convergence

mpost1 <- convertToCodaObject(models[[1]])
psrf.beta1 <- gelman.diag(mpost1$Beta, multivariate=FALSE)$psrf
summary(psrf.beta1)
quantile(psrf.beta1[, 1], seq(0, 1, 0.05))

mpost2 <- convertToCodaObject(models[[2]])
psrf.beta2 <- gelman.diag(mpost2$Beta, multivariate=FALSE)$psrf
summary(psrf.beta2)
quantile(psrf.beta2[, 1], seq(0, 1, 0.05))

predY1 <- computePredictedValues(models[[1]], expected=FALSE, nParallel=2)
MF1 <- evaluateModelFit(hM=models[[1]], predY=predY1)
mean(MF1$TjurR2, na.rm=T)  # 0.38
mean(MF1$AUC, na.rm=T)  # 0.92
#WAIC <- computeWAIC(models[[1]])

predY2 <- computePredictedValues(models[[2]], expected=FALSE, nParallel=2)
MF2 <- evaluateModelFit(hM=models[[2]], predY=predY2)
mean(MF2$R2, na.rm=T)  # 0.36

# # two-fold cross validation.
 partition <- createPartition(models[[2]], nfolds=2)
 preds_new <- computePredictedValues(models[[2]], partition=partition, nParallel=2)
 MF_pred <- evaluateModelFit(hM=models[[2]], predY=preds_new)
# 
 save.image(file='HMSC_ABS/models_updated/mcmc_convergence.RData')


## get posterior estimates

# beta
postBeta1 <- getPostEstimate(models[[1]], parName='Beta', start=1)
postBeta2 <- getPostEstimate(models[[2]], parName='Beta', start=1)

# omega
OmegaCor1 <- computeAssociations(models[[1]], start=1)
OmegaCor2 <- computeAssociations(models[[2]], start=1)


## save results
# save(models, mpost1, mpost2, postBeta1, postBeta2, OmegaCor1, OmegaCor2, file='output/workspace/hmsc_reference/mcmc_output.RData')

