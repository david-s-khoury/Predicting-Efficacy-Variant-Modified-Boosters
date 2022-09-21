# geometric mean
geomean = function(x){
  10^(mean(log10(x),na.rm=T))
}

cis <- function (x,upperlower='l',takeLog=T, rounding = 1){
  if (takeLog){
    x=log10(x)
  }
  a=c()
  if (upperlower=='l' | upperlower==.25){
    a=c(a,mean(x)-1.96*sd(x)/sqrt(length(x)))
  }
  if(upperlower=='u'| upperlower==.75){
    a=c(a,mean(x)+1.96*sd(x)/sqrt(length(x)))
  }
  if (takeLog){
    a=10^a
  }
  round(a, rounding)
}


get_eff_after_over_var = function (str){
  variable = as.vector(str_split(str_split(str, 'after', simplify=T),'over',simplify=T))
  var1 = str_trunc(variable[1], nchar(variable[1])-1,'right',ellipsis = '')
  var2 = str_trunc(variable[2], nchar(variable[2])-1,'left',ellipsis = '')
  time = as.vector(str_split(var2,'_',simplify=T))[1]
  time_end = str_locate(var2,'_')[1]
  if (!is.na(time_end)){
    var2 = str_trunc(var2, nchar(var2)-time_end+1,'left',ellipsis = '')
  } else {
    var2=''
  }
  return_var = paste0(var1, var2)
  return_var
  
}
get_timing = function(str){
  variable = as.vector(str_split(str_split(str, 'after', simplify=T),'over',simplify=T))
  var1 = str_trunc(variable[1], nchar(variable[1])-1,'right',ellipsis = '')
  var2 = str_trunc(variable[2], nchar(variable[2])-1,'left',ellipsis = '')
  time = as.vector(str_split(var2,'_',simplify=T))[1]
  time
}

# Calculate the efficacy after a certain amount of time (based on the current efficacy)
# default time period is 6 months
# Note: efficacy_type = 3 is symptomatic disease
# efficacy_type = 4 is severe disease
eff_after_time = function(eff, efficacy_type=3, time_period = 364/2, min_eff = 0){
  
  log10_e_factor = log10(exp(1))
  
  decayrate = log(2)/108
  # Get the current neuts
  neut = get_neut_from_efficacy(eff, efficacy_type)
  #decay them
  log10decayed_neut = log10(neut)-time_period*decayrate*log10_e_factor 
  
  # This is where we need to set to a minimum of the pre-boost levels. 
  # Probably need another input variable.
  
  # Get the new efficacy
  decayed_eff = pmax(get_efficacy_from_neut(10^log10decayed_neut,efficacy_type), min_eff)
  decayed_eff
}

# Calculate the difference between two starting efficacies after a certain period of time
eff_diff_after_time = function(eff1, eff2=0, efficacy_type=3, time_period = 364/2, min_eff1 = 0, min_eff2 = 0){
  eff1 = eff_after_time(eff1, efficacy_type, time_period, min_eff1)
  eff2 = eff_after_time(eff2, efficacy_type, time_period, min_eff2)
  return_val = eff1-eff2
  return_val
}

# Calculate the average difference in efficacies over a certain period of time based on the two starting efficacies
eff_diff_over_time = function(eff_vect1, eff_vect2=0, efficacy_type=3, time_period = 364/2, min_eff_vect1=0, min_eff_vect2=0){
  time_series = seq(0,time_period)
  eff_avg_vect = eff_vect1
  if(length(min_eff_vect1)==1){
    min_eff_vect1 = rep(min_eff_vect1, length(eff_vect1))
  }
  if(length(min_eff_vect2)==1){
    min_eff_vect2 = rep(min_eff_vect2, length(eff_vect2))
  }
  for (i in c(1:length(eff_vect1))){
    effs1 = eff_after_time(eff_vect1[i], efficacy_type, time_series, min_eff_vect1[i])
    effs2 = eff_after_time(eff_vect2[i], efficacy_type, time_series, min_eff_vect2[i])
    eff_avg_vect[i] = integrate(splinefun(time_series,effs1-effs2), 0, time_period)$value/time_period
  }
  eff_avg_vect
}

# Aferage efficacy over a period of time based ont eh startign efficacies
eff_avg_over_time = function(eff_vect, efficacy_type=3, time_period = 364/2, min_eff_vect=0){
  time_series = seq(0,time_period)
  eff_avg_vect = eff_vect
  i=1
  if(length(min_eff_vect)==1){
    min_eff_vect = rep(min_eff_vect, length(eff_vect))
  }
  # This should be done as an apply function for efficiency
  for (eff in eff_vect){
    effs = eff_after_time(eff, efficacy_type, time_series, min_eff_vect[i])
    eff_avg_vect[i] = integrate(splinefun(time_series,effs), 0, time_period)$value/time_period
    i=i+1
  }
  eff_avg_vect
}

# case_diff_over_time = function(eff1, eff2=0, efficacy_type=3, time_period = 364/2){
#   decayrate = log(2)/108
#   time_series = seq(0,time_period)
#   
#   if (time_period == 0){
#     return (1000*abs(eff1-eff2))
#   }
#   diffs=c()
#   neut1 = get_neut_from_efficacy(eff1, efficacy_type)
#   if (eff2!=0){
#     neut2 = get_neut_from_efficacy(eff2, efficacy_type)
#   }
#   for (i in c(1:length(neut1))){
#     neut_series1 = log10(neut1[i])-time_series*decayrate
#     eff_series1 = get_efficacy_from_neut(10^neut_series1,efficacy_type)
#     saved1 = integrate(splinefun(time_series,eff_series1), 0, time_period)$value
#     if (eff2!=0){ 
#       neut_series2 = log10(neut2[i])-time_series*decayrate 
#       eff_series2 = get_efficacy_from_neut(10^neut_series2,efficacy_type)
#       saved2 = integrate(splinefun(time_series,eff_series2), 0, time_period)$value
#     } else {
#       saved2 = rep(0,length(saved1))
#     }
#     diffs[i] = abs(saved2-saved1)/time_period
#   }
#   diffs
# }
# 
# 
# case_diff_after_time = function(eff1, eff2, efficacy_type=3, time_period = 364/2){
#   neut1 = get_neut_from_efficacy(eff1, efficacy_type)
#   neut2 = get_neut_from_efficacy(eff2, efficacy_type)
#   
#   decayrate = log(2)/108
#   neut1 = log10(neut1)-time_period*decayrate
#   neut2 = log10(neut2)-time_period*decayrate
#   eff1 = get_efficacy_from_neut(10^neut1,efficacy_type)
#   eff2 = get_efficacy_from_neut(10^neut2,efficacy_type)
#   return_val = case_diff_over_time(eff1, eff2, efficacy_type, 1)
#   return_val
# }
# 
