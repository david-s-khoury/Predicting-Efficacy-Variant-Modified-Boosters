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