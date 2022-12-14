### LOO Comparison among models 
#Note: Analysing remotly

####  Calling packages ####
#Analysis
library(ape)
library(brms)
library(posterior)

library(tidybayes)

#Setting pathway to files
setwd("Files pathway")

load("FULLinterBayesianModel.RData")

load("BayesianInteractionsFULLModel4.RData")

load("PAIREDinterBayesianModel.RData")

load("ADDIinterBayesianModel2.RData.RData")

load("BayesianSST22Symbionts2.RData")

load("BayesianSimplerModel0.RData")

# Renaming objects from the analyses (as many as models generated)

log0 <- logliktable1[!rowSums(!is.finite(logliktable1)),]
log22 <- logliktable22[!rowSums(!is.finite(logliktable22)),]
logADDI2 <- logliktableADDI2[!rowSums(!is.finite(logliktableADDI2)),]
logPAIRED <- logliktablePAIRED[!rowSums(!is.finite(logliktablePAIRED)),]
log4 <- logliktable4[!rowSums(!is.finite(logliktable4)),]
logI <- logliktableI[!rowSums(!is.finite(logliktableI)),]

loolog1 <- loo(log0)
loolog2 <- loo(log22) 
loolog3 <- loo(logADDI2)
loolog4 <- loo(logPAIRED)
loolog5 <- loo(log4) 
loolog6 <- loo(logI) 
#Note: I could get loo() of each object and then calculate manually the elpd_diff or...
#Compare all models 
loo_models <- loo_compare(loolog1,loolog2,loolog3,loolog4,loolog5,loolog6)


save.image('Save RData in specified pathway')