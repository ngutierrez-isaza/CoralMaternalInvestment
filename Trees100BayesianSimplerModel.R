# R version 4.1.1 (R Core Development Team 2021)


####  Calling packages ####
#Analysis
library(ape)
library(brms)
library(posterior)

library(tidybayes)

#Setting pathway to files
setwd("Files pathway")

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
#   FULL  Model (model_repeat1)
#######################################################

### Work with a randomly selected sub-group of trees  #####
spp_phylotable1 <- matrix(NA,ncol=7) 
colnames(spp_phylotable1) <- c('Estimate', 'Est.Error', 'l-95% CI', 'u-95% CI', 'Rhat','Bulk_ESS','Tail_ESS')

specie_nametable1 <- matrix(NA,ncol=7)
colnames(specie_nametable1) <- c('Estimate', 'Est.Error', 'l-95% CI', 'u-95% CI', 'Rhat','Bulk_ESS','Tail_ESS')

intertable1 <- data.frame(matrix(ncol = 4))
colnames(intertable1) <- c('.chain','.iteration','.draw','b_Intercept')

SymbiontsYestable1 <- data.frame(matrix(ncol = 4))
colnames(SymbiontsYestable1) <- c('.chain','.iteration','.draw','b_SymbiontsYes')

logliktable1 <- matrix(NA, ncol=236)

seed_brms = 123

hyp.list1 <-list()

models.list1 <-list()
A <-list()

for (i in 1:length(randomtrees)) {
  A[[i]] <- ape::vcv.phylo(randomtrees[[i]])
  
  model_repeat1 <- brm(
    lg.eggsize ~ Symbionts +(1|gr(spp_phylo, cov = A)) + (1|species_name)+(1|Location_name), 
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
  
  
  randomsp1 <- summary(model_repeat1)$random$species_name
  randomphylo1 <- summary(model_repeat1)$random$spp_phylo
  specie_nametable1 <- rbind(specie_nametable1, randomsp1)
  spp_phylotable1 <- rbind(spp_phylotable1, randomphylo1)
  
  intercept1 <- as.data.frame(spread_draws(model_repeat1, b_Intercept))
  symbiont1<- as.data.frame(spread_draws(model_repeat1, b_SymbiontsYes))
  
  intertable1 <- rbind(intertable1, intercept1)
  SymbiontsYestable1 <- rbind(SymbiontsYestable1, symbiont1)
  
  loglik1 <- log_lik(model_repeat1)
  logliktable1 <- rbind(logliktable1, loglik1)
  
  models.list1<-c(models.list1,list(summary(model_repeat1)))
  
  hyp1 <- paste(
    "sd_spp_phylo__Intercept^2 /", 
    "(sd_spp_phylo__Intercept^2 + sd_species_name__Intercept^2 + sd_Location_name__Intercept^2 + sigma^2) = 0"
  )
  hyptable1 <- hypothesis(model_repeat1, hyp1, class = NULL)
  hyp.list1 <- c(hyp.list1, list(hyptable1))
  
}


save.image('Save RData in specified pathway')