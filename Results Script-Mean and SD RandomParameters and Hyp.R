# Calculating mean and SD values for random parameters ####
library(dplyr)

# Means
colMeans(location_nametable0, na.rm=T)
colMeans(specie_nametable0, na.rm=T)
colMeans(spp_phylotable0, na.rm=T)

# SD
spp_phylotable2 %>% summarise_if(na.rm=T, is.numeric, sd)
# Calculating sd one column at the time
sd(specie_nametable1[,3], na.rm=T)

# Calculating mean and SD values for the hypothesis list ####
library(dplyr)
hyp2 <- sapply(hyp.list2, '[', 1)

# Means
bind_rows(hyp2, .id = 'id') %>% 
  mutate(Estimate   = as.numeric(Estimate  ))  %>%
  group_by(id) %>% 
  summarise(Estimate = mean(Estimate, na.rm = TRUE))

bind_rows(hyp2, .id = 'id') %>% 
  mutate(CI.Lower   = as.numeric(CI.Lower  ))  %>%
  group_by(id) %>% 
  summarise(CI.Lower = mean(CI.Lower, na.rm = TRUE))

bind_rows(hyp2, .id = 'id') %>% 
  mutate(CI.Upper   = as.numeric(CI.Upper  ))  %>%
  group_by(id) %>% 
  summarise(CI.Upper = mean(CI.Upper, na.rm = TRUE))

# SD
bind_rows(hyp2, .id = 'id') %>% 
  mutate(Estimate   = as.numeric(Estimate  ))  %>%
  group_by(id) %>% 
  summarise(Estimate = sd(Estimate, na.rm = TRUE))

bind_rows(hyp2, .id = 'id') %>% 
  mutate(CI.Lower   = as.numeric(CI.Lower  ))  %>%
  group_by(id) %>% 
  summarise(CI.Lower = sd(CI.Lower, na.rm = TRUE))

bind_rows(hyp2, .id = 'id') %>% 
  mutate(CI.Upper   = as.numeric(CI.Upper  ))  %>%
  group_by(id) %>% 
  summarise(CI.Upper = sd(CI.Upper, na.rm = TRUE))


### Figures ####
plot(model_repeat2, N=2, ask=FALSE) # trace plots used in Appendix S3, Fig S3.4a

pp_check(model_repeat2, nsamples = 100) # Figure S3.4b

plot(conditional_effects(model_repeat2), points = F) %>% par(bg="white") # Figure S3.5a
plot(marginal_effects(model_repeat1)) # Figure S3.5b

