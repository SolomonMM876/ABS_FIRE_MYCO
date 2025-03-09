# THIS SCRIPT EVALUATES MCMC CONVERGENCE OF HMSC MODELS OF THE FUNGAL EXAMPLE (SECTION 7.9) OF THE BOOK

# The preliminaries are as in script S1

# Change the working directory if necessary

modelDir = file.path("HMSC_ABS/models")
library(Hmsc)
set.seed(1)

# We first read in the object models that includes all six fitted models
# We recommed running the script first with setting thin = 1
# Then rerun it with thin = 10, thin = 100, ... to examine how the convergence statistics improve



nChains = 2
samples = 10
thin = 10
filename=file.path(modelDir, paste0("models_thin_",as.character(thin),"_samples_",as.character(samples),"_chains_",as.character(nChains),'.Rdata'))
load(filename)

nChains = 4
samples = 250
thin = 100
filename=file.path(modelDir, paste0("model_chains_",as.character(nChains),"_samples_",
                                    as.character(samples),"_thin_",as.character(thin)))
load(filename)

# We evaluate here MCMC convergence only for the $\beta$ and $\Omega$ parameters, and only in terms of the potential scale reduction factor.
# While we could evaluate MCMC convergence also in terms of effective sample size, we note that usually
# examining the potential scale reduction factor is more critical: if convergence is satisfactory in terms of it,
# then it is usually satisfactory also in terms of effective sample size

par(mfrow=c(3,2))
for (i in 1:3){
  mpost = convertToCodaObject(models[[i]][[1]])#remove the 1 later when you run it
  psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  psrf.omega = gelman.diag(mpost$Omega[[1]],multivariate=FALSE)$psrf
  mymain = switch (i, "lognormal poisson", "Probit", "Normal")
  myxlab1 = switch (i, "", "", "psrf (beta)")
  myxlab2 = switch (i, "", "", "psrf (Omega)")
  hist(psrf.beta, main=mymain, xlab = myxlab1)
  hist(psrf.omega, main=mymain, xlab = myxlab2)
}
dev.off()
#this is the alternative way to look at the model outputs:
par(mfrow=c(3,2))

  mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
  psrf.beta = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
  psrf.omega = gelman.diag(mpost$Omega[[1]],multivariate=FALSE)$psrf
  mymain = "Probit"
  myxlab1 = "psrf (beta)"
  myxlab2 ="psrf (Omega)"
  hist(psrf.beta, main=mymain, xlab = myxlab1)
  hist(psrf.omega, main=mymain, xlab = myxlab2)

dev.off()

##################################################################################################
# MAKE THE SCRIPT REPRODUCIBLE (BEGINNING)
##################################################################################################
set.seed(1)
##################################################################################################
## MAKE THE SCRIPT REPRODUCIBLE (END)
##################################################################################################


##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (BEGINNING)
##################################################################################################
showBeta = NULL #Default: showBeta = TRUE, convergence shown for beta-parameters
showGamma = NULL #Default: showGamma = FALSE, convergence not shown for gamma-parameters
showOmega = NULL #Default: showOmega = FALSE, convergence not shown for Omega-parameters
maxOmega = NULL #Default: convergence of Omega shown for 50 randomly selected species pairs
showRho = NULL #Default: showRho = FALSE, convergence not shown for rho-parameters
showAlpha = NULL #Default: showAlpha = FALSE, convergence not shown for alpha-parameters
##################################################################################################
# SETTING COMMONLY ADJUSTED PARAMETERS TO NULL WHICH CORRESPONDS TO DEFAULT CHOICE (END)
##################################################################################################

##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (BEGINNING)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################
showBeta = TRUE
showGamma = TRUE
showOmega = TRUE
maxOmega = 100
showRho = TRUE
showAlpha = TRUE
##################################################################################################
# CHANGE DEFAULT OPTIONS BY REMOVING COMMENT AND SETTING VALUE (END)
# NOTE THAT THIS IS THE ONLY SECTION OF THE SCRIPT THAT YOU TYPICALLY NEED TO MODIFY
##################################################################################################

##################################################################################################
# SET DIRECTORIES (BEGINNING)
##################################################################################################
localDir="HMSC_ABS"
modelDir = file.path(localDir, "models")
resultDir = file.path(localDir, "results")
##################################################################################################
# SET DIRECTORIES (END)
##################################################################################################

if(is.null(showBeta)) showBeta = TRUE
if(is.null(showGamma)) showGamma = FALSE
if(is.null(showOmega)) showOmega = FALSE
if(is.null(maxOmega)) maxOmega = 50
if(is.null(showRho)) showRho = FALSE
if(is.null(showAlpha)) showAlpha = FALSE

library(Hmsc)
library(colorspace)
library(vioplot)

samples_list = c(10)
  #c(250,250,250,250,250,250)
thin_list = c(10)#c(1,10,100,500,1000)
nst = length(thin_list)
nChains = 2

text.file = file.path(resultDir,"/MCMC_convergence.txt")
cat("MCMC Convergennce statistics\n\n",file=text.file,sep="")

ma.beta = NULL
na.beta = NULL
ma.gamma = NULL
na.gamma = NULL
ma.omega= NULL
na.omega = NULL
ma.alpha = NULL
na.alpha = NULL  
ma.rho = NULL
na.rho = NULL
Lst = 1


while(Lst <= nst){
  thin = thin_list[Lst]
  samples = samples_list[Lst]
  
  
  filename = file.path(modelDir,paste("models_thin_", as.character(thin),
                                      "_samples_", as.character(samples),
                                      "_chains_",as.character(nChains),
                                      ".Rdata",sep = ""))
  if(file.exists(filename)){
    load(filename)
    cat(c("\n",filename,"\n\n"),file=text.file,sep="",append=TRUE)
    nm = length(models)
    for(j in 1:nm){
      mpost = convertToCodaObject(models[[j]], spNamesNumbers = c(T,F), covNamesNumbers = c(T,F))
      nr = models[[j]]$nr
      cat(c("\n",names(models)[j],"\n\n"),file=text.file,sep="",append=TRUE)
      if(showBeta){
        psrf = gelman.diag(mpost$Beta,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\nbeta\n\n",file=text.file,sep="",append=TRUE)
        cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
        if(is.null(ma.beta)){
          ma.beta = psrf[,1]
          na.beta = paste0(as.character(thin),",",as.character(samples))
        } else {
          ma.beta = cbind(ma.beta,psrf[,1])
          if(j==1){
            na.beta = c(na.beta,paste0(as.character(thin),",",as.character(samples)))
          } else {
            na.beta = c(na.beta,"")
          }
        }
      }
      if(showGamma){
        psrf = gelman.diag(mpost$Gamma,multivariate=FALSE)$psrf
        tmp = summary(psrf)
        cat("\ngamma\n\n",file=text.file,sep="",append=TRUE)
        cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
        if(is.null(ma.gamma)){
          ma.gamma = psrf[,1]
          na.gamma = paste0(as.character(thin),",",as.character(samples))
        } else {
          ma.gamma = cbind(ma.gamma,psrf[,1])
          if(j==1){
            na.gamma = c(na.gamma,paste0(as.character(thin),",",as.character(samples)))
          } else {
            na.gamma = c(na.gamma,"")
          }
        }
      }
      if(showRho & !is.null(mpost$Rho)){
        psrf = gelman.diag(mpost$Rho,multivariate=FALSE)$psrf
        cat("\nrho\n\n",file=text.file,sep="",append=TRUE)
        cat(psrf[1],file=text.file,sep="\n",append=TRUE)
      }
      if(showOmega & nr>0){
        cat("\nomega\n\n",file=text.file,sep="",append=TRUE)
        for(k in 1:nr){
          cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
          tmp = mpost$Omega[[k]]
          z = dim(tmp[[1]])[2]
          if(z > maxOmega){
            sel = sample(1:z, size = maxOmega)
            for(i in 1:length(tmp)){
              tmp[[i]] = tmp[[i]][,sel]
            }
          }
          psrf = gelman.diag(tmp, multivariate = FALSE)$psrf
          tmp = summary(psrf)
          cat(tmp[,1],file=text.file,sep="\n",append=TRUE)
          if(is.null(ma.omega)){
            ma.omega = psrf[,1]
            na.omega = paste0(as.character(thin),",",as.character(samples))
          } else {
            ma.omega = cbind(ma.omega,psrf[,1])
            if(j==1){
              na.omega = c(na.omega,paste0(as.character(thin),",",as.character(samples)))
            } else {
              na.omega = c(na.omega,"")
            }
          }
        }
      }
      if(showAlpha & nr>0){
        for(k in 1:nr){
          if(models[[j]]$ranLevels[[k]]$sDim>0){
            cat("\nalpha\n\n",file=text.file,sep="\n",append=TRUE)
            cat(c("\n",names(models[[j]]$ranLevels)[k],"\n\n"),file=text.file,sep="",append=TRUE)
            psrf = gelman.diag(mpost$Alpha[[k]],multivariate = FALSE)$psrf
            cat(psrf[,1],file=text.file,sep="\n",append=TRUE)            
          }
        }
      }
    }
  }
  Lst = Lst + 1
}

pdf(file= file.path(resultDir,"/MCMC_convergence.pdf"))
if(showBeta){
  par(mfrow=c(2,1))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0,max(ma.beta)),main="psrf(beta)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.beta,col=rainbow_hcl(nm),names=na.beta,ylim=c(0.9,1.1),main="psrf(beta)")
}
if(showGamma){
  par(mfrow=c(2,1))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0,max(ma.gamma)),main="psrf(gamma)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.gamma,col=rainbow_hcl(nm),names=na.gamma,ylim=c(0.9,1.1),main="psrf(gamma)")
}
if(showOmega & !is.null(ma.omega)){
  par(mfrow=c(2,1))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0,max(ma.omega)),main="psrf(omega)")
  legend("topright",legend = names(models), fill=rainbow_hcl(nm))
  vioplot(ma.omega,col=rainbow_hcl(nm),names=na.omega,ylim=c(0.9,1.1),main="psrf(omega)")
}
dev.off()