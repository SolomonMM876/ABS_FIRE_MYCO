# THIS SCRIPT EXPLORES THE PARAMETER ESTIMATES OF HMSC MODELS OF THE FUNGAL EXAMPLE (SECTION 7.9) OF THE BOOK

# The preliminaries are as in script S1

# Change the working directory if necessary
# wd = "C:/HMSC_book/R-scripts/Section_7_9_fungi"
# setwd(wd)

# localDir = "HMSC"
#data.directory = file.path(localDir, "data")
model.directory = file.path("HMSC/models")
library(Hmsc)
set.seed(1)


# We first read in the object models that includes all six fitted models
# You may change the thin parameter to the highest thin for which you have fitted the model, to get as reliable results as possible

nChains = 2
samples = 100
thin = 1
filename=file.path(model.directory, paste0("models_chains_",as.character(nChains),"_samples_",as.character(samples),"_thin_",as.character(thin)))
load(filename)

#We first plot the responses of the species to the environmental covariates in the models
# that include both sequencing depth and Fire Regime (the second index is 2 in m = models[[i]][[2]]).

# By setting the value of i the lognormal Poisson model
# on sequence counts (i=1), the probit model on presence-absence data (i=2), and normal model on log-transformed 
# sequence count conditional on presence (i=3). Note that only the last plot may be shown in your screen,
# so you need to scroll back with the images to see the other ones
i=2
m = models[[i]][[2]]
postBeta = getPostEstimate(m, parName="Beta")
plotBeta(m, post=postBeta, param="Support", supportLevel = 0.95, spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))

# In these heatmaps of estimated species niches. Red (respectively, blue) colour shows parameters that are estimated
# to be positive (respectively, negative) with at least 0.95 posterior probability. The intensity of the colour refers
# to the posterior mean estimate of the parameter. 


# interpret from the beta-plot and best visualised through gradient plots.
# We next construct a gradient plot over decay classes to the presence-absence model (i=2) 
# For simplicity, we call this model(models[[2]][[2]]) as m

m = models[[2]][[2]]


# We explore the gradient over decay class, and thus we set focalVariable = "DC"
# We wish to normalize readcount to its mean over all the data, and thus we set non.focalVariables = list("readcount"=list(1))
# See help for the function constructGradient to see the options for non-focal variables

Gradient = constructGradient(m,focalVariable = "Fire.Severity",#Fire.Interval
                             non.focalVariables = list("readcount"=list(1),
                                                       "Fire.Interval"=list(1)))

# It is usually a good idea to look at the Gradient object that has been constructed, to ensure that it is as it should be
# Thus, you may evaluate

Gradient$XDataNew

# We next make predicted species communities over this gradient

predY = predict(m,Gradient = Gradient, expected = TRUE)

# Let's explore this object:

class(predY)

# It is a list...

length(predY)

# ..of length 200. This is because the predictions are done for each posterior sample, and we stored 
# 100 samples for each of the two chains, thus 200 in total.

dim(predY[[1]])

# each prediction has a dimension of 2 x 19, as the predictions are made for 4 sampling units (as explained above)
# for each of the 31 species.

head(predY[[1]])


# We next visualize the predictions with the plotGradient function, either for individual species
# (measure="Y", index selecting the species) or species richness (measure="S"). With the option
# showData, one can decide whether the raw data is plotted as well. With the option jigger,
# one can avoid overlapping points. With the option prob, one can choose the credible intervals to be plotted 

prob = c(0.25,0.5,0.75)
# summary(m$Y) #select species of interest 
plotGradient(m, Gradient, pred=predY, measure="Y", index=23, showData = TRUE, jigger = 0.15, prob = prob)
plotGradient(m, Gradient, pred=predY, measure="Y", index=16, showData = TRUE, jigger = 0.15, prob = prob)
plotGradient(m, Gradient, pred=predY, measure="S", showData = TRUE, jigger = 0.15, prob = prob)

#no differences in species richness between the Fire Severities or Fire Frequencies 

# We next plot the species association networks. To do so, we need the corrplot package

library(corrplot)

# We plot the associations here only for the presence-absence model (i=2) for the version
# that does not account for Fire Regime in the fixed effects (j=1)

m = models[[2]][[1]]

# To plot the association network for any of the other models, simply pick another part of the models object
# We first use computeAssociations to derive the correlations from the object m

OmegaCor = computeAssociations(m)

# Let's explore this object

class(OmegaCor)

# It is a list..

length(OmegaCor)

# ...of length 1. This is because we included only one random effect to the model, i.e. that at the sampling unit level

class(OmegaCor[[1]])

# Also the first element of this list is a list. It is a named list with posterior mean (mean) values and the 
# posterior supports for the correlations being positive (support). These can be explore numerically as follows:

head(OmegaCor[[1]]$mean)
head(OmegaCor[[1]]$support)

# We next plot the network of assocications, colouring only those that are either positive or negative
# with at least 95% posterior probability

supportLevel = 0.95
toPlot = ((OmegaCor[[1]]$support>supportLevel) + (OmegaCor[[1]]$support<(1-supportLevel))>0)* OmegaCor[[1]]$mean
dev.off
corrplot(toPlot, method = "color", col=colorRampPalette(c("blue","white","red"))(200))


# With m = models[[2]][[1]] you will visualize the "raw associations" of the presence-absence model
# Repeating the above with m = models[[2]][[2]] you will visualize the "residual associations" of the presence-absence model
# There are many more raw association than residual associations. 
# This is typically to be expected, as the raw associations are also generated by
# differential responses of the species to the environmental conditions, here the decay class.

# We may also visualize the association network as a model-based ordination using the biPlot function
# Let us first do that e.g. for the raw associations of the presence-absence model 

m = models[[2]][[1]]
biPlot(m,etaPost = getPostEstimate(m,"Eta"),lambdaPost = getPostEstimate(m,"Lambda"), colVar = 1)

# The circles are the sampling units and the triangles the dots. The sampling units are coloured
# according the Fire Severity because we selected colVar = 1 and Fire Severity class is the first column of XData

head(m$XData)

# Note that even if decay class is part of XData, in this model it is not part of XFormula,
# and thus not part of the X-matrix:

head(m$X)

# Let us do a similar plot for residual associations:

m = models[[2]][[2]]
biPlot(m,etaPost = getPostEstimate(m,"Eta"),lambdaPost = getPostEstimate(m,"Lambda"), colVar = 1)










