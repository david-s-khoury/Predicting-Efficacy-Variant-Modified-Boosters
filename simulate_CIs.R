# 1) Generate the parameter set

# 2) for each set, generate the eight curves at 1% efficacy intervals - this is what takes a long time

# Now - each 1% efficacy is associated with a neut Ab, adn that decays over time, 
# so we can pre-calculate the neuts at whcih we need to calculate teh efficacies, 
# assuming we know the decay rates

# 3) for each of the 8 curves, extract the (0.025,.5,.975) quantiles at each efficacy interval


# Which method are we using to calculate the CIs? 
# Combined combines two intergrals and is MUCH faster

use.method = 'combined'

make.eff_avged = F
make.new.params = F

# Set up the starting efficacies
lower_sympt_effs = seq(1e-5,9e-4,1e-5)
#lower_sympt_effs = c(1e-5,5e-5,seq(1e-4,9e-4,1e-4))
upper_sympt_effs = seq(.001,.999,.001)
#upper_sympt_effs = seq(.01,.99,.01)

baseline_sympt_effs = c(lower_sympt_effs,upper_sympt_effs)
start_effs = baseline_sympt_effs
gap = .25
#start_effs = c(0.01,seq(gap,.99,gap),.99)
#start_effs = c(1e-5,seq(2e-5,1e-3,2.5e-5),seq(2e-3,1e-2,5e-4),seq(2e-2,.98,2e-2),seq(.99,.999,.001))

nruns=10000
savenameCIs =  paste0('./simulatedCIoutput/simulated_CI_results_n',nruns,'_neff',length(start_effs),'_',paste('params',paste(c(ancestral_boost,vsv_rises),collapse ='_'),sep='_'))

if(make.eff_avged){
  
  if(make.new.params){
    # Are we making a new parameter set?
    
    # depermine the parameters
    mean_log10boost = mean(anc_boost_log10fc)
    se_log10boost = sd(anc_boost_log10fc)/sqrt(length(anc_boost_log10fc))
    
    mean_log10VSVimprovement = mean(log10(vsv_improvement_variant))
    se_log10VSVimprovement = sd(log10(vsv_improvement_variant))/sqrt(length(vsv_improvement_variant))
    mean_log10VSVimprovement_matched = mean(log10(vsv_improvement_variant_matched))
    se_log10VSVimprovement_matched = sd(log10(vsv_improvement_variant_matched))/sqrt(length(vsv_improvement_variant_matched))
    mean_log10VSVimprovement_nonmatched = mean(log10(vsv_improvement_variant_nonmatched))
    se_log10VSVimprovement_nonmatched = sd(log10(vsv_improvement_variant_nonmatched))/sqrt(length(vsv_improvement_variant_nonmatched))
    
    
    boost.improv.params = c(mean_log10boost,se_log10boost,
                            mean_log10VSVimprovement,se_log10VSVimprovement,
                            mean_log10VSVimprovement_matched, se_log10VSVimprovement_matched,
                            mean_log10VSVimprovement_nonmatched, se_log10VSVimprovement_nonmatched)
    
    names(boost.improv.params) = c('log10_boost', 'log10_boostSD', 'log10_vsvimprovement', 'log10_vsvimprovementSD',
                                   'log10_vsvimprovement_matched', 'log10_vsvimprovement_matchedSD',
                                   'log10_vsvimprovement_nonmatched', 'log10_vsvimprovement_nonmatchedSD')
    params = select_parameters_from_distributions(nruns, boost.improv.params)
  }
  
  sim_starttime=Sys.time()
  
  # now start the simulation
  #These are what we will output
  col.names = c("sympt.unboosted","sympt.boosted",  "sympt.boosted.improved",
                "sympt.boosted.improved.matched","sympt.boosted.improved.nonmatched",  
                "severe.unboosted", "severe.boosted", "severe.boosted.improved",
                "severe.boosted.improved.matched", "severe.boosted.improved.nonmatched",
                "sympt.improved.diff",
                "sympt.improved.matched.diff",
                "sympt.improved.nonmatched.diff",
                "severe.improved.diff",
                "severe.improved.matched.diff",
                "severe.improved.nonmatched.diff")
  # Name teh rows by the baseline starting efficacy
  row.names = 100*round(start_effs,3)
  matrix.names = paste0('param.set.',c(1:nrow(params)))
  
  #This is an array of the averaged effs that we will fill
  effs_avged=array(NA,dim=c(length(start_effs), length(col.names),nrow(params)),
                   dimnames = list(row.names, col.names, matrix.names))
  for (i in c(1:nrow(params))){
    effs_avged[,,i] = calculate_eff_avg_over_time_from_effs_and_params(start_effs,params[i,], method=use.method)
    if((i%%100)==0){
      sim_midtime=Sys.time()
      sim_midruntime=sim_midtime-sim_starttime
      print(paste0('Completed ',i,' iterations, (t=',format(round(sim_midruntime,2)),')'))
      save(effs_avged, file = paste0(savenameCIs,'_intermediate_i=',i,'.RData'))
    }
  }
  
  effs_avged_matchedvsnonmatched=array(NA,dim=c(length(start_effs), 2,nruns),
                   dimnames = list(row.names, c('sympt.matchedvsnonmatched.boost.diff','severe.matchedvsnonmatched.boost.diff'), matrix.names))
  effs_avged_matchedvsnonmatched[,'sympt.matchedvsnonmatched.boost.diff',] = effs_avged[,'sympt.boosted.improved.matched',]-effs_avged[,'sympt.boosted.improved.nonmatched',]
  effs_avged_matchedvsnonmatched[,'severe.matchedvsnonmatched.boost.diff',] = effs_avged[,'severe.boosted.improved.matched',]-effs_avged[,'severe.boosted.improved.nonmatched',]
  
  # Combine
  effs_avged2 = array(NA,dim=c(length(start_effs), 18,nruns),
                      dimnames = list(row.names, c(col.names,'sympt.matchedvsnonmatched.boost.diff','severe.matchedvsnonmatched.boost.diff'), matrix.names))
  effs_avged2[,1:16,]=effs_avged
  effs_avged2[,17:18,]=effs_avged_matchedvsnonmatched
  effs_avged = effs_avged2
  
  # Now generate the quantiles
  quantile.results = list()
  quantile.results$lower95.effs_avged = effs_avged[,,1]
  quantile.results$upper95.effs_avged = effs_avged[,,1]
  quantile.results$median.effs_avged = effs_avged[,,1]
  use.quantiles = c(.025,.5,.975)


  for(m in c(1:nrow(effs_avged[,,1]))){
    for (n in c(1:ncol(effs_avged[,,1]))){
      quantile.values = quantile(effs_avged[m,n,],use.quantiles)
      quantile.results$lower95.effs_avged[m,n] =  quantile.values[1] 
      quantile.results$median.effs_avged[m,n] =  quantile.values[2] 
      quantile.results$upper95.effs_avged[m,n] =  quantile.values[3] 
    }
  }
  
  sim_endtime=Sys.time()
  sim_runtime=sim_endtime-sim_starttime
  print(paste0('method = ',use.method))
  print(paste0('Ran ',nruns,' iterations with vector of length ',length(start_effs),': RunTime = ',format(round(sim_runtime,2))))
  
  ## Now we need to combine
  save(effs_avged, quantile.results, start_effs, file=paste0(savenameCIs,'_final.RData'))
} else {
  load(file=paste0(savenameCIs,'_final.RData'))
}


# Now combined the quantiles into one dataframe that includes variables for lower, upper and median
quantile.results.combined = data.frame()
for (i in c(1:length(quantile.results))){
  quantile.results.combined.temp = reshape2::melt(quantile.results[[i]]) %>% 
    rename(baseline_eff = Var1, eff_type = Var2, new_eff = value) 
  # Name the column
  names(quantile.results.combined.temp)[3] = strsplit(names(quantile.results)[i],'\\.')[[1]][1]
  
  if (nrow(quantile.results.combined)==0) {
    quantile.results.combined = quantile.results.combined.temp
  } else {
    quantile.results.combined = bind_cols(quantile.results.combined,quantile.results.combined.temp[,3])
    names(quantile.results.combined)[ncol(quantile.results.combined)] = names(quantile.results.combined.temp)[3]
  } 
}

# Add these things so that we can do the plotting
quantile.results.combined = quantile.results.combined %>%
  #rename(baseline_sympt_eff = start_eff, avg_eff_over_6m=median) %>%
  mutate(GMR=eff_type,
         GMRdiff=eff_type,
         outcome = ifelse(str_split(as.character(quantile.results.combined$eff_type),'\\.', simplify = T)[,1]=='sympt','Symptomatic','Severe'),
         eff_type_group = ifelse(outcome=='Symptomatic', str_remove(eff_type,'sympt.'),str_remove(eff_type,'severe.')),
         lower95adj = lower95,
         upper95adj = upper95
  )
quantile.results.combined$outcome = factor(quantile.results.combined$outcome, levels = c('Symptomatic','Severe'))

# Fix the CIs so there is no overlap on the difference plot
quantile.results.combined$lower95adj[quantile.results.combined$eff_type_group=='improved.matched.diff']=quantile.results.combined$upper95[quantile.results.combined$eff_type_group=='improved.diff']
quantile.results.combined$upper95adj[quantile.results.combined$eff_type_group=='improved.nonmatched.diff']=quantile.results.combined$lower95[quantile.results.combined$eff_type_group=='improved.diff']


# Add the colour scales for plotting
newefficacycolourscale = gmr_colours
names(newefficacycolourscale) = c('unboosted','boosted', 'boosted.improved','boosted.improved.matched','boosted.improved.nonmatched')
newimprovementcolourscale = improvement_cols
names(newimprovementcolourscale) = paste0(c('improved','improved.matched', 'improved.nonmatched'),'.diff')
newalphascale = c(.15,.05,.09)
names(newalphascale) = paste0(c('improved','improved.matched', 'improved.nonmatched'),'.diff')

# Add in the confidence intervals to the plots
sub_plot1_avg_over_6m_CI=sub_plot1_avg_over_6m+
  ggnewscale::new_scale_fill()+
  geom_ribbon(data=filter(quantile.results.combined, eff_type_group=='boosted' | eff_type_group=='boosted.improved' ),aes(x=baseline_eff/100,y=median,ymin=lower95/100, ymax=upper95/100, fill=eff_type_group),colour=NA, alpha=.15)+
  scale_fill_manual(values=newefficacycolourscale,name=fold_improvement_legend_name, labels = gmr_labels, guide='none')
# print(sub_plot1_avg_over_6m_CI)
ggsave(paste0(dir$plots,'Average_Saved_with_CIs.pdf'),sub_plot1_avg_over_6m_CI, width=12,height=6)

sub_plot_improvement_avg_over_6m_by_outcome_CI = sub_plot_improvement_avg_over_6m_by_outcome +
  ggnewscale::new_scale_fill()+
  geom_ribbon(data=filter(quantile.results.combined, startsWith(eff_type_group,'improved')),aes(x=baseline_eff/100,y=median,ymin=lower95adj, ymax=upper95adj, fill=eff_type_group, alpha=eff_type_group),colour=NA)+
  scale_fill_manual(values=newimprovementcolourscale,name=fold_improvement_legend_name, labels = gmr_labels, guide='none')+
  coord_cartesian(xlim = x_expand, ylim = y_improvement_expand, expand=F)+
  scale_alpha_manual(values=newalphascale)
ggsave(paste0(dir$plots,'Average_Improvement_with_CIs.pdf'),sub_plot_improvement_avg_over_6m_by_outcome_CI, width=12,height=6)

  
