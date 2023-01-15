# Simulation parameters

# What are the varying parameters
# 1) Sigma
mean_sigma=SummaryTable_Efficacy_NeutRatio_SD_SEM$PooledSD[1]
sd_sigma=SummaryTable_Efficacy_NeutRatio_SD_SEM$SE_PooledSD[1]

# 2) Hill Coeff
# 3) IC50
cov_hill_IC50<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPool$hessian)[9:10,9:10]
mean_hill_IC50=c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool$estimate,2))
sd_hill_IC50=sqrt(diag(cov_hill_IC50))

# 4) Efficacy -> Starting Neuts assume normal dist of starting log10 neut ratios
# Use pfizer as a baseline
#study = 'Pfizer'
#eff_init= SummaryTable_Efficacy_NeutRatio_SD_SEM$Efficacy[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==study]
#eff_lower95 = SummaryTable_Efficacy_NeutRatio_SD_SEM$Lower[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==study]/100
#eff_upper95 = SummaryTable_Efficacy_NeutRatio_SD_SEM$Upper[SummaryTable_Efficacy_NeutRatio_SD_SEM$Study==study]/100

#sd_log10neutR = estimate_sd_from_CIs(mean_log10neutR,log10neutR_low,log10neutR_high)

# 6) Shift to severe
#mean_shift = 20/3
#sd_shift = .05*mean_shift
#Note that  think the sd ofo this shift is too small becuase the off diagonal elements for the hill coefficient are not considered. 
#Really we should look at a 3x3 matrix that includes the last 3 elements and accounts for hill, IC50 and shift
indicies<-c(-1,0)+length(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate)
CovS<-solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$hessian)[indicies,indicies]
mean_shift=sum(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate[indicies]*c(-1,1))
sd_shift=sqrt(sum(c(diag(CovS)^2,2*CovS[1,2]*CovS[2,1])))

# 7) hill, IC50 for severe
cov_hill_IC50_severe = solve(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$hessian)[15:16,15:16]
mean_hill_IC50_severe=c(tail(FittedLogistic_RawEfficacy_MeanRept_SDPool_Diff_EC50_Only$estimate,3))[1:2]
sd_hill_IC50_severe=sqrt(diag(cov_hill_IC50_severe))

# 8) Ab decay rate 
mean_decay = log(2)/108
sd_decay = estimate_sd_from_CIs(mean_decay,log(2)/159,log(2)/82)



select_parameters_from_distributions = function(nruns, boostimprovementparams){
  
  bip = boostimprovementparams
  paramnames = c('sigma','hill','IC50','decay','severe_shift','hill_severe','IC50_severe',
                 'log10_boost','log10_VSVimprovement','log10_VSVimprovement_matched','log10_VSVimprovement_nonmatched')
  means = c(mean_sigma, mean_hill_IC50, mean_decay,mean_shift,mean_hill_IC50_severe, 
            bip['log10_boost'], bip['log10_vsvimprovement'], bip['log10_vsvimprovement_matched'], bip['log10_vsvimprovement_nonmatched'])
  sds = c(sd_sigma, sd_hill_IC50, sd_decay, sd_shift,sd_hill_IC50_severe,
          bip['log10_boostSD'], bip['log10_vsvimprovementSD'], bip['log10_vsvimprovement_matchedSD'], bip['log10_vsvimprovement_nonmatchedSD'])
  names(means)=paramnames
  names(sds)=paramnames
  paramseeds = randomLHS(nruns,length(means))
  params = data.frame(paramseeds)
  colnames(params)<-paramnames
  tstart=Sys.time()
  for (i in c(1:(ncol(paramseeds)))){
    params[,i] = qnorm(paramseeds[,i], mean=means[i],sd = sds[i])   
  }
  params[colnames(params) %in% c('hill','IC50')]= rmvnorm(nruns,mean=mean_hill_IC50, sigma = cov_hill_IC50)
  params[colnames(params) %in% c('hill_severe','IC50_severe')]= rmvnorm(nruns,mean=mean_hill_IC50_severe, sigma = cov_hill_IC50_severe)
  
  params  
}

# select_base_parameters_from_distributions = function(withBoosting=F){
#   paramnames = c('sigma','hill','IC50','decay','severe_shift','hill_severe','IC50_severe')
#   means = c(mean_sigma, mean_hill_IC50, mean_decay,mean_shift,mean_hill_IC50_severe)
#   sds = c(sd_sigma, sd_hill_IC50, sd_decay, sd_shift,sd_hill_IC50_severe)
#   names(means)=paramnames
#   names(sds)=paramnames
#   paramseeds = rep(0.5,length(means))
#   params = paramseeds
#   names(params)<-paramnames
#   for (i in c(1:(length(paramseeds)-2))){
#     params[i] = qnorm(paramseeds[i], mean=means[i],sd = sds[i])   
#   }
#   
#   if(withBoosting){
#     nvax = length(vaccines)
#   } else{
#     nvax = sum(vaccines!='Boosted')
#   }
#   nvariants = length(variants_to_use)
#   paramstruct = data.frame(matrix(rep(params,each=nvariants*nvax),nrow=nvariants*nvax))
#   colnames(paramstruct) = paramnames
#   paramstruct$variant = rep(c(0:(nvariants-1))/nvariants,each=nvax)
#   paramstruct$serum = rep(c(0:(nvax-1))/length(vaccines),nvariants)
#   paramstruct
# }