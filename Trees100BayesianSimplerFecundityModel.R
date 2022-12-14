# R version 4.1.1 (R Core Development Team 2021)


####  Calling packages ####
#Analysis
library(ape)
library(brms)
library(posterior)

# Packages required to follow vignette tidybayes-brms

library(tidybayes)

 

#Setting pathway to files
setwd("Files pathway")


###### Calling the database ####
numeggs <- read.csv("FecundityHHT.csv", header =T) 
###### TREES

###### Calling the multiPhylo (1000 phylogenetic trees)
Supertrees<-read.tree("supertreeMultiPhylo.tre")   # 1000 phylogenetic trees, Object from the class multiPhylo


##### Prune supertree for analysis
####   Matching database with tree tips #### 
pruned.Supertrees20 <-lapply(Supertrees,drop.tip,Supertrees[[1]]$tip.label[-na.omit(match(numeggs$spp_phylo,Supertrees[[1]]$tip.label))])    # R base function lapply
class(pruned.Supertrees20) <-"multiPhylo"   # assigning class is necessary because lapply returns a generic list with no class assigned

numeggs <- droplevels(numeggs[numeggs$spp_phylo %in% pruned.Supertrees20[[1]]$tip.label, ])   # remove species that don't match

### Defining cofactors (just spec_mean_cf)
numeggs$lg_fecundity <-log(numeggs$Eggs_per_area)

numeggs$spec_mean_cfFec <- 
  with(numeggs, sapply(split(lg_fecundity, spp_phylo), mean)[spp_phylo])

numeggs$within_spec_cfFec <- numeggs$lg_fecundity - numeggs$spec_mean_cfFec

randomtrees <-sample(pruned.Supertrees20, size=100, replace = FALSE)

#######################################################
#   FULL  Model Gaussian (model_repeat2) Intercept model
#######################################################

### Work with a randomly selected sub-group of trees  #####
spp_phylotable2 <- matrix(NA,ncol=7) 
colnames(spp_phylotable2) <- c('Estimate', 'Est.Error', 'l-95% CI', 'u-95% CI', 'Rhat','Bulk_ESS','Tail_ESS')

specie_nametable2 <- matrix(NA,ncol=7)
colnames(specie_nametable2) <- c('Estimate', 'Est.Error', 'l-95% CI', 'u-95% CI', 'Rhat','Bulk_ESS','Tail_ESS')

intertable2 <- data.frame(matrix(ncol = 4))
colnames(intertable2) <- c('.chain','.iteration','.draw','b_Intercept')

logliktable2 <- matrix(NA, ncol=35)

seed_brms = 123

hyp.list2 <-list()

models.list2 <-list()
A <-list()

for (i in 1:length(randomtrees)) {
  A[[i]] <- ape::vcv.phylo(randomtrees[[i]])
  
  model_repeat2 <- brm(
    lg.eggsize ~ 1 +(1|gr(spp_phylo, cov = A)) + (1|species_name)+(1|Location_name), 
    data = numeggs, 
    family = gaussian(), 
    data2 = list(A = A[[i]]),
    prior = c(
      prior(normal(0,50), "Intercept"),
      prior(student_t(3,0,20), "sd"),
      prior(student_t(3,0,20), "sigma")
    ),
    sample_prior = TRUE, chains = 3, cores = 3, 
    iter = 10000, warmup = 1000, thin = 5, control =list(adapt_delta=0.9999, max_treedepth = 15), seed = seed_brms
  )
  
  
  randomsp2 <- summary(model_repeat2)$random$species_name
  randomphylo2 <- summary(model_repeat2)$random$spp_phylo
  specie_nametable2 <- rbind(specie_nametable2, randomsp2)
  spp_phylotable2 <- rbind(spp_phylotable2, randomphylo2)
  
  intercept2 <- as.data.frame(spread_draws(model_repeat2, b_Intercept))
  
  intertable2 <- rbind(intertable2, intercept2)
  
  loglik2 <- log_lik(model_repeat2)
  logliktable2 <- rbind(logliktable2, loglik2)
  
  models.list2<-c(models.list2,list(summary(model_repeat2)))
  
  hyp2 <- paste(
    "sd_spp_phylo__Intercept^2 /", 
    "(sd_spp_phylo__Intercept^2 + sd_species_name__Intercept^2 + sd_Location_name__Intercept^2 + sigma^2) = 0"
  )
  hyptable2 <- hypothesis(model_repeat2, hyp2, class = NULL)
  hyp.list2 <- c(hyp.list2, list(hyptable2))
  
}


save.image('Save RData in specified pathway')