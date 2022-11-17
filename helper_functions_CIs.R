source('./select_bootstrap_params.R')
# This calculates the average efficacy over time based on starting (mean) efficacy (for symptomatic)
# It is a wrapper function for if the neuts are input.
calculate_eff_avg_over_time_from_effs_and_params = function(effs, curve_boost_params, time_period = 364/2, min_eff_vect=0, method = 'standard'){
  log10neutRs = log10(get_neut_from_efficacy(effs))
  eff_avg = calculate_eff_avg_over_time_from_neuts_and_params(log10neutRs, curve_boost_params, time_period, min_eff_vect, method)
  eff_avg
}

# This calculates the average efficacy over time based on starting log10 of neut Ratios (for symptomatic)
calculate_eff_avg_over_time_from_neuts_and_params = function(log10neutRs, curve_boost_params, time_period = 364/2, min_eff_vect=0, method='standard'){
  
  input.effs = get_efficacy_from_neut(10^log10neutRs)
  
  eff_avg.matrix = matrix(0,ncol=16,nrow=length(log10neutRs))
  if(method=='standard'){
    time_series = seq(0,time_period)
    effs = calculate_eff_after_time_from_neuts_and_params(log10neutRs, curve_boost_params, time_series, min_eff_vect)
  } else if (method=='combined'){
    #effs = calculate_eff_after_time_from_neuts_and_params(log10neutRs, curve_boost_params, time_period, min_eff_vect, method)
    effs = calculate_eff_from_neuts_and_params(log10neutRs, curve_boost_params, time_period, min_eff_vect,method) 
  }
  
  # Add all the results that were in a list into an output matrix
  for (eff.type.num in c(1:length(effs))){
    eff.type = names(effs)[eff.type.num]
    for (input.eff.num in c(1:length(input.effs))){
      input.eff = input.effs[input.eff.num]
      
      if (method=='standard') {
        if (is.null(nrow(effs[[eff.type.num]]))){
          efficacy.series = effs[[eff.type.num]]
        } else {
          efficacy.series = effs[[eff.type.num]][input.eff.num,]
        }
        eff_avg.matrix[input.eff.num,eff.type.num] = integrate( splinefun(time_series,efficacy.series), 0, time_period)$value/time_period
      } else if (method == 'combined'){
        eff_avg.matrix[input.eff.num,eff.type.num] = effs[[eff.type.num]][input.eff.num]
      }
    }
    
  }
  # Name the rows and columns of the matrix
  colnames(eff_avg.matrix) = names(effs)
  rownames(eff_avg.matrix) = 100*round(input.effs,3)
  
  eff_avg.matrix
}

calculate_eff_from_neuts_and_params<- function(log10neutRs.matrix, curve_boost_params,time_period=NA, min_eff, method='standard'){
  
  # Extract the boost / fold improvement params
  log10_boost = curve_boost_params$log10_boost
  log10_VSVimprovement = curve_boost_params$log10_VSVimprovement
  log10_VSVimprovement_matched = curve_boost_params$log10_VSVimprovement_matched
  log10_VSVimprovement_nonmatched = curve_boost_params$log10_VSVimprovement_nonmatched
  if(log10_VSVimprovement<0){
    print(paste0('log10_VSVimprovement < 0 (',round(log10_VSVimprovement,3),')'))
  }
  # Number of time columns input (this will be 1 if the method is combined, and the whole decay interval if the method is standard)
  n_time_cols = ncol(log10neutRs.matrix)
  if(is.null(n_time_cols)){
    n_time_cols=1
  }
  # Starting neuts
  entire.log10neutR.matrix = cbind(log10neutRs.matrix,
                                   log10neutRs.matrix+log10_boost,
                                   log10neutRs.matrix+log10_boost+log10_VSVimprovement,
                                   log10neutRs.matrix+log10_boost+log10_VSVimprovement_matched,
                                   log10neutRs.matrix+log10_boost+log10_VSVimprovement_nonmatched)
  
  # Convert the matrix to a vector (by columns)
  entire.log10neutR.vector = as.vector(entire.log10neutR.matrix)
  
  if(method == 'standard') {
    sympt_eff = 100*LogisticModel_PercentUninfected(entire.log10neutR.vector,
                                                    curve_boost_params$sigma, curve_boost_params$hill, curve_boost_params$IC50)
    severe_eff = 100*LogisticModel_PercentUninfected(entire.log10neutR.vector,
                                                     curve_boost_params$sigma, curve_boost_params$hill_severe, curve_boost_params$IC50_severe)
  }
  else if (method == 'combined') {
    # This is a new method that calculates the average efficacy over a time period, given a decay rate.
    sympt_eff = 100*LogisticModel_PercentUninfected_OverTime(entire.log10neutR.vector,
                                                             curve_boost_params$sigma, curve_boost_params$hill, curve_boost_params$IC50,
                                                             curve_boost_params$decay, time_period)
    severe_eff = 100*LogisticModel_PercentUninfected_OverTime(entire.log10neutR.vector,
                                                              curve_boost_params$sigma, curve_boost_params$hill_severe, curve_boost_params$IC50_severe,
                                                              curve_boost_params$decay, time_period)
    
  }
  # Now convert back to a matrix
  entire.sympt_eff.matrix = matrix(sympt_eff,nrow=nrow(entire.log10neutR.matrix))
  entire.severe_eff.matrix = matrix(severe_eff,nrow=nrow(entire.log10neutR.matrix))
  
  # Add everything to a list
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
  
  # This is the difference in paired efficacies
  effs$sympt.improved.diff  = effs$sympt.boosted.improved-effs$sympt.boosted
  effs$sympt.improved.matched.diff  = effs$sympt.boosted.improved.matched-effs$sympt.boosted
  effs$sympt.improved.nonmatched.diff  = effs$sympt.boosted.improved.nonmatched-effs$sympt.boosted
  effs$severe.improved.diff  = effs$severe.boosted.improved-effs$severe.boosted
  effs$severe.improved.matched.diff  = effs$severe.boosted.improved.matched-effs$severe.boosted
  effs$severe.improved.nonmatched.diff  = effs$severe.boosted.improved.nonmatched-effs$severe.boosted
  
  # Return the list of efficacies
  effs
}

calculate_eff_after_time_from_neuts_and_params = function(log10neutRs, curve_boost_params, time_period = 364/2, min_eff = 0, method='deborah'){
  
  if (method == 'combined'){
    # If we are using the combined function this is an unnecessary step and should never get called
    decayed_eff = calculate_eff_from_neuts_and_params(log10neutRs, curve_boost_params, time_period, min_eff,method)  
  } else if (method=='standard'){
    # Matrix of times - each row is the times
    time_period.matrix = matrix(rep(time_period,length(log10neutRs)),nrow=length(log10neutRs), byrow = T)
    
    log10_e_factor = log10(exp(1))
    # Extract the decay rate we are using
    decayrate = curve_boost_params$decay
    # Matrix of neuts - each column is the neuts we have entered
    log10neutRs.matrix = matrix(rep(log10neutRs,length(time_period)),nrow=length(log10neutRs))
    #decay the neuts - this now needs to be a matrix
    log10decayed_neut.matrix = log10neutRs.matrix-time_period.matrix*decayrate*log10_e_factor # Dont do this for david
    # Get the new efficacy
    decayed_eff = calculate_eff_from_neuts_and_params(log10decayed_neut.matrix, curve_boost_params, time_period=NA, min_eff)
  } 
  
  decayed_eff
}

load('./natmed_parameters.RData')

sig=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1]
logk = tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2)[1]

# Fitted IC50 values
C50=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,1)
C50_Severe=tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Severe$estimate,1)

decay = log(2)/108

ProbRemainUninfected_AverageTime=function(log10Titre,logk,C50_log10){1/(1+exp(-exp(logk)*(log10Titre-C50_log10)))}

LogisticModel_Symptomatic_OverTime=function(mu_titre_init_log10,TimeInterval=182){
  eff = mu_titre_init_log10
  # Only use the non-infinite indices
  use.inds = which(!is.infinite(mu_titre_init_log10))
  eff[use.inds] = LogisticModel_PercentUninfected_OverTime(mu_titre_init_log10[use.inds], sig, logk, C50, decay, TimeInterval)
  eff[is.infinite(eff) & eff < 0] = 0
  eff[is.infinite(eff) & eff > 0] = 1
  eff
}

LogisticModel_Severe_OverTime=function(mu_titre_init_log10,TimeInterval=182){
  eff = mu_titre_init_log10
  # Only use the non-infinite indices
  use.inds = which(!is.infinite(mu_titre_init_log10))
  eff[use.inds] = LogisticModel_PercentUninfected_OverTime(mu_titre_init_log10[use.inds], sig, logk, C50_Severe, decay, TimeInterval)
  eff[is.infinite(eff) & eff < 0] = 0
  eff[is.infinite(eff) & eff > 0] = 1
  eff
}

LogisticModel_PercentUninfected_OverTime=function(mu_titre_init_log10,sig_titre_log10,logk,C50_log10,DecayRate,TimeInterval){
  NumInteration<-max(length(mu_titre_init_log10),length(sig_titre_log10),length(logk),length(C50_log10),length(DecayRate),length(TimeInterval))
  Output<-NULL
  
  decayrate_log10scale=DecayRate*log10(exp(1))
  
  if (length(C50_log10)==1) {
    C50_log10=rep(C50_log10,NumInteration)
  }
  
  if (length(logk)==1) {
    logk=rep(logk,NumInteration)
  }
  
  if (length(sig_titre_log10)==1) {
    sig_titre_log10=rep(sig_titre_log10,NumInteration)
  }
  
  if (length(mu_titre_init_log10)==1) {
    mu_titre_init_log10=rep(mu_titre_init_log10,NumInteration)
  }
  
  if (length(decayrate_log10scale)==1) {
    decayrate_log10scale=rep(decayrate_log10scale,NumInteration)
  }
  
  if (length(TimeInterval)==1) {
    TimeInterval=rep(TimeInterval,NumInteration)
  }
  
  for (i in 1:NumInteration) {
    Step=sig_titre_log10[i]*0.001
    
    IntegralVector=seq(mu_titre_init_log10[i]-5*sig_titre_log10[i],mu_titre_init_log10[i]+5*sig_titre_log10[i],by=Step)
    
    Output[i]=(1/(decayrate_log10scale[i]*TimeInterval[i]))*sum(ProbRemainUninfected(IntegralVector,logk[i],C50_log10[i])*(pnorm(IntegralVector-mu_titre_init_log10[i]+decayrate_log10scale[i]*TimeInterval[i],0,sig_titre_log10[i])-pnorm(IntegralVector-mu_titre_init_log10[i],0,sig_titre_log10[i])))*Step
  }
  Output
}