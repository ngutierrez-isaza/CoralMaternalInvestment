########################################################################
#    Results
#######################################################################


# Model validation
#####  trace plots 
plot(model_repeat1, N=2, ask=FALSE)
##### posterior predictive checks
pp_check(model_repeat1, nsamples = 100)


#### plot conditional effects
plot(conditional_effects(model_repeat1), points = T)

#Note: cofactor_name can be replaced by intercept (intertable), spec_mean_cf(spec_meantable), within_spec_cf (within_spectable) or any parameter from the model that want to be checked

