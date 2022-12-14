# R version 4.1.1 (R Core Development Team 2021)


####  Calling packages ####
#Analysis
library(ape)
library(brms)
library(posterior)

library(tidybayes)

#Setting pathway to files
setwd("Files pathway")


###### Temperature analysis across species

###### Calling the database ####
numeggs <- read.csv("MeanLifehistoryCompleteDatabase.csv", header =T) 
###### TREES

###### Calling the multiPhylo (1000 phylogenetic trees)
Supertrees<-read.tree("supertreeMultiPhylo.tre")   # 1000 phylogenetic trees, Object from the class multiPhylo


##### Prune supertree for analysis
####   Matching database with tree tips #### 
pruned.Supertrees20 <-lapply(Supertrees,drop.tip,Supertrees[[1]]$tip.label[-na.omit(match(numeggs$spp_phylo,Supertrees[[1]]$tip.label))])    # R base function lapply
class(pruned.Supertrees20) <-"multiPhylo"   # assigning class is necessary because lapply returns a generic list with no class assigned

numeggs <- droplevels(numeggs[numeggs$spp_phylo %in% pruned.Supertrees20[[1]]$tip.label, ])   # remove species that don't match

### Defining cofactors (just spec_mean_cf)
numeggs$spec_mean_cf <- 
  with(numeggs, sapply(split(sstmean, spp_phylo), mean)[spp_phylo])

numeggs$spec_mean_cfVAR <- 
  with(numeggs, sapply(split(sstvarall, spp_phylo), mean)[spp_phylo])

numeggs$spec_mean_cfChla <- 
  with(numeggs, sapply(split(chlamean, spp_phylo), mean)[spp_phylo])  

randomtrees <-sample(pruned.Supertrees20, size=100, replace = FALSE)

#######################################################
#   FULL  Model (model_repeat4)
#######################################################

### Work with a randomly selected sub-group of trees  #####
spp_phylotable4 <- matrix(NA,ncol=7) 
colnames(spp_phylotable4) <- c('Estimate', 'Est.Error', 'l-95% CI', 'u-95% CI', 'Rhat','Bulk_ESS','Tail_ESS')

specie_nametable4 <- matrix(NA,ncol=7)
colnames(specie_nametable4) <- c('Estimate', 'Est.Error', 'l-95% CI', 'u-95% CI', 'Rhat','Bulk_ESS','Tail_ESS')

intertable4 <- data.frame(matrix(ncol = 4))
colnames(intertable4) <- c('.chain','.iteration','.draw','b_Intercept')

spec_meantable4 <- data.frame(matrix(ncol = 4))
colnames(spec_meantable4) <- c('.chain','.iteration','.draw','b_polyspec_mean_cf21')

spec_meantable42 <- data.frame(matrix(ncol = 4))
colnames(spec_meantable42) <- c('.chain','.iteration','.draw','b_polyspec_mean_cf22')

spec_mean_cfVARtable4 <- data.frame(matrix(ncol = 4))
colnames(spec_mean_cfVARtable4) <- c('.chain','.iteration','.draw','b_spec_mean_cfVAR')

spec_mean_cfChlatable4 <- data.frame(matrix(ncol = 4))
colnames(spec_mean_cfChlatable4) <- c('.chain','.iteration','.draw','b_spec_mean_cfChla')

SymbiontsYestable4 <- data.frame(matrix(ncol = 4))
colnames(SymbiontsYestable4) <- c('.chain','.iteration','.draw','b_SymbiontsYes')

ModeReproductionHermaphroditetable4 <- data.frame(matrix(ncol = 4))
colnames(ModeReproductionHermaphroditetable4) <- c('.chain','.iteration','.draw','b_ModeReproductionHermaphrodite')

logliktable4 <- matrix(NA, ncol=236)

seed_brms = 123

hyp.list4 <-list()

models.list4 <-list()
A <-list()

for (i in 1:length(randomtrees)) {
  A[[i]] <- ape::vcv.phylo(randomtrees[[i]])
  
  model_repeat4 <- brm(
    lg.eggsize ~ poly(spec_mean_cf, 2)*spec_mean_cfVAR*spec_mean_cfChla + ModeReproduction*Symbionts +(1|gr(spp_phylo, cov = A)) + (1|species_name)+(1|Location_name), 
    data = numeggs, 
    family = gaussian(), 
    data2 = list(A = A[[i]]),
    prior = c(
      prior(normal(0,10), "b"),
      prior(normal(0,50), "Intercept"),
      prior(student_t(3,0,20), "sd"),
      prior(student_t(3,0,20), "sigma")
    ),
    sample_prior = TRUE, chains = 3, cores = 3, 
    iter = 10000, warmup = 1000, thin = 5, control =list(adapt_delta=0.9999, max_treedepth = 15), seed = seed_brms
  )
  
  
  randomsp4 <- summary(model_repeat4)$random$species_name
  randomphylo4 <- summary(model_repeat4)$random$spp_phylo
  specie_nametable4 <- rbind(specie_nametable4, randomsp4)
  spp_phylotable4 <- rbind(spp_phylotable4, randomphylo4)
  
  intercept4 <- as.data.frame(spread_draws(model_repeat4, b_Intercept))
  meantemp4<- as.data.frame(spread_draws(model_repeat4, b_polyspec_mean_cf21))
  meantemp42<- as.data.frame(spread_draws(model_repeat4, b_polyspec_mean_cf22))
  tempvar4<- as.data.frame(spread_draws(model_repeat4, b_spec_mean_cfVAR))
  chloro4<- as.data.frame(spread_draws(model_repeat4, b_spec_mean_cfChla))
  symbiont4<- as.data.frame(spread_draws(model_repeat4, b_SymbiontsYes))
  modereprod4<- as.data.frame(spread_draws(model_repeat4, b_ModeReproductionHermaphrodite))
  
  intertable4 <- rbind(intertable4, intercept4)
  spec_meantable4 <- rbind(spec_meantable4, meantemp4)
  spec_meantable42 <- rbind(spec_meantable42, meantemp42)
  spec_mean_cfVARtable4 <- rbind(spec_mean_cfVARtable4, tempvar4)
  spec_mean_cfChlatable4 <- rbind(spec_mean_cfChlatable4, chloro4)
  SymbiontsYestable4 <- rbind(SymbiontsYestable4, symbiont4)
  ModeReproductionHermaphroditetable4 <- rbind(ModeReproductionHermaphroditetable4, modereprod4)
  
  loglik4 <- log_lik(model_repeat4)
  logliktable4 <- rbind(logliktable4, loglik4)
  
  models.list4<-c(models.list4,list(summary(model_repeat4)))
  
  hyp4 <- paste(
    "sd_spp_phylo__Intercept^2 /", 
    "(sd_spp_phylo__Intercept^2 + sd_species_name__Intercept^2 + sd_Location_name__Intercept^2 + sigma^2) = 0"
  )
  hyptable4 <- hypothesis(model_repeat4, hyp4, class = NULL)
  hyp.list4 <- c(hyp.list4, list(hyptable4))
  
}


save.image('Save RData in specified pathway')