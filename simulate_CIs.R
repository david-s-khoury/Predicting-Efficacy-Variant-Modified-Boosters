# 1) Generate the parameter set

# 2) for each set, generate the eight curves at 1% efficacy intervals - this is what takes a long time

# Now - each 1% efficacy is associated with a neut Ab, adn that decays over time, 
# so we can pre-calculate the neuts at whcih we need to calculate teh efficacies, 
# assuming we know the decay rates

# 3) for each of the 8 curves, extract the (0.025,.5,.975) quantiles at each efficacy interval


calculate_eff_from_neuts_and_params<- function(log10neutRs.matrix, curve_boost_params, min_eff){
  log10_boost = curve_boost_params$log10_boost
  log10_VSVimprovement = curve_boost_params$log10_VSVimprovement
  log10_VSVimprovement_matched = curve_boost_params$log10_VSVimprovement_matched
  log10_VSVimprovement_nonmatched = curve_boost_params$log10_VSVimprovement_nonmatched
  
  # Number of time columns input
  n_time_cols = ncol(log10neutRs.matrix)
  
  # Basic neuts
  entire.log10neutR.matrix = cbind(log10neutRs.matrix,
                               log10neutRs.matrix+log10_boost,
                               log10neutRs.matrix+log10_boost+log10_VSVimprovement,
                               log10neutRs.matrix+log10_boost+log10_VSVimprovement_matched,
                               log10neutRs.matrix+log10_boost+log10_VSVimprovement_nonmatched)
  
  # Convert the matrix to a vector (by columns)
  entire.log10neutR.vector = as.vector(entire.log10neutR.matrix)
  
  sympt_eff = 100*LogisticModel_PercentUninfected(entire.log10neutR.vector,
                                                  curve_boost_params$sigma, curve_boost_params$hill, curve_boost_params$IC50)
  severe_eff = 100*LogisticModel_PercentUninfected(entire.log10neutR.vector,
                                                   curve_boost_params$sigma, curve_boost_params$hill_severe, curve_boost_params$IC50_severe)
  
  #unboosted_eff = 100*LogisticModel_PercentUninfected(log10neutR,curve_boost_params$sigma,curve_boost_params$hill,curve_boost_params$IC50)
  #unboosted_severe_eff = 100*LogisticModel_PercentUninfected(log10neutR,curve_boost_params$sigma,curve_boost_params$hill_severe,curve_boost_params$IC50_severe)
  
  #boosted_eff = 100*LogisticModel_PercentUninfected(log10neutR+log10_boost,curve_boost_params$sigma,curve_boost_params$hill,curve_boost_params$IC50)
  #boosted_severe_eff = 100*LogisticModel_PercentUninfected(log10neutR+log10_boost,curve_boost_params$sigma,curve_boost_params$hill_severe,curve_boost_params$IC50_severe)
  
  #boosted_improved_eff = 100*LogisticModel_PercentUninfected(log10neutR+log10_boost+log10_VSVimprovement,curve_boost_params$sigma,curve_boost_params$hill,curve_boost_params$IC50)
  #boosted_improved_severe_eff = 100*LogisticModel_PercentUninfected(log10neutR+log10_boost+log10_VSVimprovement,curve_boost_params$sigma,curve_boost_params$hill_severe,curve_boost_params$IC50_severe)
  #effs = data.frame(unboosted_eff,boosted_eff,boosted_improved_eff,unboosted_severe_eff,boosted_severe_eff,boosted_improved_severe_eff)
  #colnames(effs) = c('unboosted_sympt','boosted_sympt','boosted_improved_sympt','unboosted_severe','boosted_severe','boosted_improved_severe')
  
  # Now convert back to matrix
  entire.sympt_eff.matrix = matrix(sympt_eff,nrow=nrow(entire.log10neutR.matrix))
  entire.severe_eff.matrix = matrix(sympt_eff,nrow=nrow(entire.log10neutR.matrix))
  
  effs = list()
  i = 0
  effs$sympt.unboosted = pmax(entire.sympt_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  effs$sympt.boosted = pmax(entire.sympt_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  effs$sympt.boosted.improved = pmax(entire.sympt_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  effs$sympt.boosted.improved.matched = pmax(entire.sympt_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  effs$sympt.boosted.improved.nonmatched = pmax(entire.sympt_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  
  i = 0
  effs$severe.unboosted = pmax(entire.severe_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  effs$severe.boosted = pmax(entire.severe_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  effs$severe.boosted.improved = pmax(entire.severe_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  effs$severe.boosted.improved.matched = pmax(entire.severe_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  effs$severe.boosted.improved.nonmatched = pmax(entire.severe_eff.matrix[,c( (i*n_time_cols+1) : ((i+1)*n_time_cols) )],min_eff)
  i=i+1
  
  
  effs
}

calculate_eff_after_time_from_neuts_and_params = function(log10neutRs, curve_boost_params, time_period = 364/2, min_eff = 0){
  
  log10_e_factor = log10(exp(1))
  
  # Extract the decay rate we are using
  decayrate = curve_boost_params$decay
  #decay the neuts - this now needs to be a matrix
  
  # Matrix of neuts - each column is the neuts we have entered
  log10neutRs.matrix = matrix(rep(log10neutRs,length(time_period)),nrow=length(log10neutRs))
  
  # Matrix of times - each row is the times
  time_period.matrix = matrix(rep(time_period,length(log10neutRs)),nrow=length(log10neutRs), byrow = T)
  log10decayed_neut.matrix = log10neutRs.matrix-time_period.matrix*decayrate*log10_e_factor 
  
  # This is where we need to set to a minimum of the pre-boost levels. 
  # Probably need another input variable.
  
  # Get the new efficacy
  decayed_eff = calculate_eff_from_neuts_and_params(log10decayed_neut.matrix, curve_boost_params, min_eff)
  decayed_eff
}

# This calculates the average efficacy over time based on starting (mean) efficacy (for symptomatic)
calculate_eff_avg_over_time_from_effs_and_params = function(effs, curve_boost_params, time_period = 364/2, min_eff_vect=0){
  log10neutRs = log10(get_neut_from_efficacy(effs))
  eff_avg = calculate_eff_avg_over_time_from_neuts_and_params(log10neutRs, curve_boost_params, time_period, min_eff_vect)
  eff_avg
}

# This calculates the average efficacy over time based on starting log10 of neut Ratios (for symptomatic)
calculate_eff_avg_over_time_from_neuts_and_params = function(log10neutRs, curve_boost_params, time_period = 364/2, min_eff_vect=0){
  time_series = seq(0,time_period)
  
  input.effs = get_efficacy_from_neut(10^log10neutRs)
  effs = calculate_eff_after_time_from_neuts_and_params(log10neutRs, curve_boost_params, time_series, min_eff_vect)
  
  
  ############
  # Check from here
  ##########
  
  # Now need to integrate for each list entry and each row of that matrix.
  
  eff_avg.matrix = matrix(0,ncol=length(effs),nrow=length(log10neutRs))
  colnames(eff_avg.matrix) = names(effs)
  rownames(eff_avg.matrix) = 100*round(input.effs,3)
  for (eff.type.num in c(1:length(effs))){
    eff.type = names(effs)[eff.type.num]
    for (input.eff.num in c(1:length(input.effs))){
      input.eff = input.effs[input.eff.num]
      efficacy.series = effs[[eff.type.num]][input.eff.num,]
      eff_avg.matrix[input.eff.num,eff.type.num] = integrate( splinefun(time_series,efficacy.series), 0, time_period)$value/time_period
    }
    
  }
  eff_avg.matrix
}




if(T){
mean_log10boost = mean(anc_boost_log10fc)
sd_log10boost = sd(anc_boost_log10fc)

mean_log10VSVimprovement = mean(log10(vsv_improvement_variant))
sd_log10VSVimprovement = sd(log10(vsv_improvement_variant))
mean_log10VSVimprovement_matched = mean(log10(vsv_improvement_variant_matched))
sd_log10VSVimprovement_matched = sd(log10(vsv_improvement_variant_matched))
mean_log10VSVimprovement_nonmatched = mean(log10(vsv_improvement_variant_nonmatched))
sd_log10VSVimprovement_nonmatched = sd(log10(vsv_improvement_variant_nonmatched))


boost.improv.params = c(mean_log10boost,sd_log10boost,
                        mean_log10VSVimprovement,sd_log10VSVimprovement,
                        mean_log10VSVimprovement_matched, sd_log10VSVimprovement_matched,
                        mean_log10VSVimprovement_nonmatched, sd_log10VSVimprovement_nonmatched)

names(boost.improv.params) = c('log10_boost', 'log10_boostSD', 'log10_vsvimprovement', 'log10_vsvimprovementSD',
                               'log10_vsvimprovement_matched', 'log10_vsvimprovement_matchedSD',
                               'log10_vsvimprovement_nonmatched', 'log10_vsvimprovement_nonmatchedSD')
nruns=10
params = select_parameters_from_distributions(nruns, boost.improv.params)

gap = .1
start_effs = c(0.01,seq(gap,.99,gap),.99)

sim_starttime=Sys.time()

effs_avged=list()
for (i in c(1:nrow(params))){
  effs_avged[[i]] = calculate_eff_avg_over_time_from_effs_and_params(start_effs,params[i,])
  print(Sys.time())
}
sim_runtime=Sys.time()-sim_starttime
## Now we need to combine

print('Up to here')
#sympt_effs = effs[,c(1:3)]
#severe_effs = effs[,c(4:6)]
#colnames(sympt_effs) = c('unboosted','boosted','boosted_improved')
#colnames(severe_effs) = c('unboosted','boosted','boosted_improved')

colMeans(sympt_effs)
colMeans(severe_effs)

sim_endtime=Sys.time()
sim_runtime=sim_endtime-sim_starttime
print(paste0('Ran ',nruns,' iterations at 50% eff: RunTime = ',format(round(sim_runtime,2))))

}