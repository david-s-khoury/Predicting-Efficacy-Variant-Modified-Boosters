library(stringr)
library(ggplot2)
library(dplyr)
library(reshape2)
library(RColorBrewer)

source('./WP4_modelImplementationFunctions.R')
source('./helper_functions.R')
max_neuts = 10^5
min_neuts = 10^-5
# population size
pop = 1e6
# force of infection (unused)
foi_use=.1
foi_gap=.1
multiple_fois = FALSE
fois = ifelse (multiple_fois, seq(foi_gap,1,foi_gap), foi_use)

# potential disease outcomes
outcomes = c('Symptomatic', 'Severe')

# This just tells us what to plot / what structures to make
remake_plot_neut_eff_struct = TRUE # This needs to be set to true for the first run, 
                                  # but it takes 2-3 mins to run (on my slower machine), 
                                  #so only set to TRUE when you have changed to boosts 
                                  #and / or the fold changes you want to include in the plot
save_plot_neut_eff_struct = TRUE
reload_plot_neut_eff_struct = FALSE

do_after_sixm = TRUE
include_average_over_6m = TRUE
#This is to check results over 4, 8, 12 week period - takes much longer to run 
do_after_other_weeks = FALSE
include_average_over_other_weeks = FALSE


# Unused
# Probability of severe disease
psevere = .01

# Make sure you have run importVSVdata first 

# extract boosts with ancestral and variant vaccines
anc_boost_log10fc = filter(vsvdata,BoosterType=='Ancestral', !is.na(log10foldchange))$log10foldchange
vsv_boost_log10fc = filter(vsvdata,BoosterGroup1=='Variant',!is.na(log10foldchange))$log10foldchange

# Extract improvement with ancestral and variant vaccines as well as for matched and un-matched
vsv_improvement_ancestral = filter(vsv_improvement_data_new,VariantGroup=='Ancestral')$VSVimprovement
vsv_improvement_variant = filter(vsv_improvement_data_new,VariantGroup=='Variant')$VSVimprovement
vsv_improvement_variant_matched = filter(vsv_improvement_data_new,VariantGroup=='Variant', Matching == matching_text[1])$VSVimprovement
vsv_improvement_variant_nonmatched = filter(vsv_improvement_data_new,VariantGroup=='Variant', Matching == matching_text[2])$VSVimprovement

# This is where we should extract the matched / unmatched rises if we want to use them instead of the general variant rise

# Define the ancestral boost and set up the list of possible boosts
ancestral_boost = round(10^mean(anc_boost_log10fc),1)
boosts = c(1,seq(2,20,2), ancestral_boost)

# Define the vsv rise and set up the list of rises
vsv_rise = round(geomean(vsv_improvement_variant),2)
vsv_rise_matched = round(geomean(vsv_improvement_variant_matched),2)
vsv_rise_nonmatched = round(geomean(vsv_improvement_variant_nonmatched),2)
vsv_rises = c(vsv_rise,vsv_rise_matched,vsv_rise_nonmatched)
rises = c(seq(1,2,.25),vsv_rises) # add the matched / unmatched rises here

# Set up the baseline neuts and efficacies for indexing. 
# This is for equally spaced neuts. 
baseline_neuts = 10^seq(log10(min_neuts),log10(max_neuts),.1)
baseline_sympt_effs = get_efficacy_from_neut(baseline_neuts)
baseline_severe_effs = get_efficacy_from_neut(baseline_neuts,4)

# This is for equally spaced sympt effs either close together (better) or further apart (quicker)
lower_sympt_effs = seq(1e-5,9e-4,1e-5)
#lower_sympt_effs = c(1e-5,5e-5,seq(1e-4,9e-4,1e-4))
upper_sympt_effs = seq(.001,.999,.001)
#upper_sympt_effs = seq(.01,.99,.01)

baseline_sympt_effs = c(lower_sympt_effs,upper_sympt_effs)

baseline_neuts = get_neut_from_efficacy(baseline_sympt_effs)
baseline_severe_effs = get_efficacy_from_neut(baseline_neuts,4)

# Set up the neut and efficacy structure to use going forward
neut_eff_struct = data.frame(baseline_neut = baseline_neuts, 
                             new_neut = baseline_neuts, 
                             GMR = 1, 
                             outcome=outcomes[1],
                             baseline_eff = baseline_sympt_effs,
                             baseline_sympt_eff = baseline_sympt_effs,
                             new_eff = baseline_sympt_effs, 
                             incremental_eff=1-(1-baseline_sympt_effs)/(1-baseline_sympt_effs),
                             NNTTfactor = NaN,
                             eff_diff = 0,
                             extraAverted = 0) %>%
  rbind(
    data.frame(baseline_neut = baseline_neuts, 
               new_neut = baseline_neuts, 
               GMR = 1, 
               outcome=outcomes[2],
               baseline_eff = baseline_severe_effs, 
               baseline_sympt_eff = baseline_sympt_effs,
               new_eff = baseline_severe_effs, 
               incremental_eff=1-(1-baseline_severe_effs)/(1-baseline_severe_effs),
               NNTTfactor = NaN,
               eff_diff = 0,
               extraAverted = 0)
  )

#Multiply all the different boosting options by all the different VSV fold rises.
gmrs = sort(unique(as.vector(outer(boosts,rises))))

# Now add all these potential rises to the neut and efficacy structure.
for (gmr in tail(gmrs,-1)){
  new_neuts=baseline_neuts*gmr
  new_neuts[new_neuts>max_neuts]=NaN
  
  new_sympt_effs = get_efficacy_from_neut(new_neuts)
  new_sympt_effs[new_neuts>max_neuts]=NaN
  new_severe_effs = get_efficacy_from_neut(new_neuts,4)
  new_severe_effs[new_neuts>max_neuts]=NaN
  
  new_neut_eff_struct = data.frame(baseline_neut = baseline_neuts, 
                                   new_neut = new_neuts, 
                                   GMR = gmr,
                                   outcome=outcomes[1],
                                   baseline_eff = baseline_sympt_effs, 
                                   baseline_sympt_eff = baseline_sympt_effs,
                                   new_eff = new_sympt_effs, 
                                   incremental_eff=1-(1-new_sympt_effs)/(1-baseline_sympt_effs),
                                   NNTTfactor = 1/(new_sympt_effs-baseline_sympt_effs),
                                   eff_diff = (new_sympt_effs-baseline_sympt_effs),
                                   extraAverted = (new_sympt_effs-baseline_sympt_effs)/baseline_sympt_effs) %>%
    
    rbind(
      data.frame(baseline_neut = baseline_neuts, 
                 new_neut = new_neuts, 
                 GMR = gmr, 
                 outcome=outcomes[2],
                 baseline_eff = baseline_severe_effs, 
                 baseline_sympt_eff = baseline_sympt_effs,
                 new_eff = new_severe_effs, 
                 incremental_eff=1-(1-new_severe_effs)/(1-baseline_severe_effs),
                 NNTTfactor = 1/(new_severe_effs-baseline_severe_effs)*(1-baseline_severe_effs)/(1-baseline_sympt_effs),
                 eff_diff=(new_severe_effs-baseline_severe_effs),
                 extraAverted = (new_severe_effs-baseline_severe_effs)/baseline_severe_effs)
    )
  # Now add the latest gmr increase to the neut / efficacy structure.
  neut_eff_struct = rbindlist(list(neut_eff_struct,new_neut_eff_struct)) # rbindlist is much faster than rbind
}

# This is where we add in results for different (or only one) fois
neut_eff_struct_base = neut_eff_struct
for(foi in fois){
  neut_eff_struct_new = neut_eff_struct_base
  neut_eff_struct_new$foi = foi
  neut_eff_struct_new$NNTT = neut_eff_struct_new$NNTTfactor/foi
  neut_eff_struct_new$NNTT[neut_eff_struct_new$outcome == 'severe'] = neut_eff_struct_new$NNTT[neut_eff_struct_new$outcome == 'severe']/psevere
  if (is.null(neut_eff_struct$foi)){
    neut_eff_struct=neut_eff_struct_new
  } else{
    neut_eff_struct = rbindlist(list(neut_eff_struct,neut_eff_struct_new))
  }
}

# Ensure outcome is a factor
neut_eff_struct$outcome = factor(neut_eff_struct$outcome, levels = outcomes)

# This is where we add in the efficacies at 6 month, and the average over 6 months. 
# It takes a long time (i.e. 2-3 min), so we sometimes set remake_plot_neut_eff_struct=FALSE
gmrs_for_subplots=c(1,ancestral_boost, ancestral_boost*vsv_rises) # add in ancestral_boost*matched / unmated rises
savename =  paste0('./plot_neut_eff_struct_addweeks_',paste('GMRs',paste(gmrs_for_subplots,collapse ='_'),sep='_'),'.RDS')

if (reload_plot_neut_eff_struct){
  plot_neut_eff_struct = readRDS(savename)
} else if (remake_plot_neut_eff_struct){
  # This is the neut / efficacy structure that we use for plotting only
  plot_neut_eff_struct = filter(neut_eff_struct, foi==foi_use, GMR %in% gmrs_for_subplots)
  
  # If we want to include values after six months
  
  # Currently this calculates the baseline_eff_after 6m, but we dont need to do that if we aren't decaying the baseline 
  if(do_after_sixm){
    plot_neut_eff_struct = plot_neut_eff_struct %>% 
      mutate (new_eff_after_6m = ifelse(outcome==outcomes[1],eff_after_time(new_eff),eff_after_time(new_eff,4)),
              baseline_eff_after_sixm = ifelse(outcome==outcomes[1],eff_after_time(baseline_eff),eff_after_time(baseline_eff,4)),
      )
  }
  if(do_after_other_weeks){
    plot_neut_eff_struct = plot_neut_eff_struct %>% 
      mutate (
              new_eff_after_4w = ifelse(outcome==outcomes[1],eff_after_time(new_eff,time_period = 4*7),eff_after_time(new_eff,4,time_period = 4*7)),
              baseline_eff_after_4w = ifelse(outcome==outcomes[1],eff_after_time(baseline_eff,time_period = 4*7),eff_after_time(baseline_eff,4,time_period = 4*7)),
              
              new_eff_after_8w = ifelse(outcome==outcomes[1],eff_after_time(new_eff,time_period = 8*7),eff_after_time(new_eff,4,time_period = 8*7)),
              baseline_eff_after_8w = ifelse(outcome==outcomes[1],eff_after_time(baseline_eff,time_period = 8*7),eff_after_time(baseline_eff,4,time_period = 8*7)),
              
              new_eff_after_12w = ifelse(outcome==outcomes[1],eff_after_time(new_eff,time_period = 12*7),eff_after_time(new_eff,4,time_period = 12*7)),
              baseline_eff_after_12w = ifelse(outcome==outcomes[1],eff_after_time(baseline_eff,time_period = 12*7),eff_after_time(baseline_eff,4,time_period = 12*7))
      )
  }
  
  if (include_average_over_6m){
    plot_neut_eff_struct = plot_neut_eff_struct %>% 
      mutate (avg_eff_over_6m = ifelse(outcome ==outcomes[1],eff_avg_over_time(new_eff),eff_avg_over_time(new_eff,4)),
              avg_eff_over_6m_cf_baseline = ifelse(outcome == outcomes[1],eff_diff_over_time(new_eff, baseline_eff),
                                                   eff_diff_over_time(new_eff, baseline_eff,4))
      )
  }
  if (include_average_over_other_weeks) {
    plot_neut_eff_struct = plot_neut_eff_struct %>% 
      mutate (
              avg_eff_over_4w = ifelse(outcome ==outcomes[1],eff_avg_over_time(new_eff,time_period = 4*7),eff_avg_over_time(new_eff,4,time_period = 4*7)),
              avg_eff_over_4w_cf_baseline = ifelse(outcome == outcomes[1],eff_diff_over_time(new_eff, baseline_eff,time_period = 4*7),
                                                   eff_diff_over_time(new_eff, baseline_eff,4,time_period = 4*7)),
              
              avg_eff_over_8w = ifelse(outcome ==outcomes[1],eff_avg_over_time(new_eff,time_period = 8*7),eff_avg_over_time(new_eff,4,time_period = 8*7)),
              avg_eff_over_8w_cf_baseline = ifelse(outcome == outcomes[1],eff_diff_over_time(new_eff, baseline_eff,time_period = 8*7),
                                                   eff_diff_over_time(new_eff, baseline_eff,4,time_period = 8*7)),
              
              avg_eff_over_12w = ifelse(outcome ==outcomes[1],eff_avg_over_time(new_eff,time_period = 12*7),eff_avg_over_time(new_eff,4,time_period = 12*7)),
              avg_eff_over_12w_cf_baseline = ifelse(outcome == outcomes[1],eff_diff_over_time(new_eff, baseline_eff,time_period = 12*7),
                                                   eff_diff_over_time(new_eff, baseline_eff,4,time_period = 12*7)),
    )
  } 
  # Don't let anything decay below baseline - really should only do this if do_after_6m or include_average_over_6m are true, but not worrying about that now
  plot_neut_eff_struct = plot_neut_eff_struct %>%
    mutate(new_eff_after_6m = ifelse(new_eff_after_6m<baseline_eff, baseline_eff, new_eff_after_6m),
           baseline_eff_after_sixm = ifelse(baseline_eff_after_sixm<baseline_eff, baseline_eff, baseline_eff_after_sixm),
           avg_eff_over_6m = ifelse(avg_eff_over_6m<baseline_eff, baseline_eff, avg_eff_over_6m)
           #avg_eff_over_6m_cf_baseline = ifelse(baseline_eff_after_sixm<baseline_eff, baseline_eff, baseline_eff_after_sixm),
    )
  if(do_after_other_weeks & include_average_over_other_weeks){
    plot_neut_eff_struct = plot_neut_eff_struct %>%
      mutate(
           new_eff_after_4w = ifelse(new_eff_after_4w<baseline_eff, baseline_eff, new_eff_after_4w),
           baseline_eff_after_4w = ifelse(baseline_eff_after_4w<baseline_eff, baseline_eff, baseline_eff_after_4w),
           avg_eff_over_4w = ifelse(avg_eff_over_4w<baseline_eff, baseline_eff, avg_eff_over_4w),
           
           new_eff_after_8w = ifelse(new_eff_after_8w<baseline_eff, baseline_eff, new_eff_after_8w),
           baseline_eff_after_8w = ifelse(baseline_eff_after_8w<baseline_eff, baseline_eff, baseline_eff_after_8w),
           avg_eff_over_8w = ifelse(avg_eff_over_8w<baseline_eff, baseline_eff, avg_eff_over_8w),
           
           new_eff_after_12w = ifelse(new_eff_after_12w<baseline_eff, baseline_eff, new_eff_after_12w),
           baseline_eff_after_12w = ifelse(baseline_eff_after_12w<baseline_eff, baseline_eff, baseline_eff_after_12w),
           avg_eff_over_12w = ifelse(avg_eff_over_12w<baseline_eff, baseline_eff, avg_eff_over_12w)
    )
  }
  if(save_plot_neut_eff_struct){
    saveRDS(plot_neut_eff_struct, savename)
  }
}
# Next step is to look at the difference between ancestral and variant boosters
# We create a new structure for plotting called plot_neut_eff_struct_improvement
# This includes the improvement for a range of fold changes
# Dummy structures
plot_neut_eff_struct_improvement_sympt = data.frame(baseline_neut = baseline_neuts, 
                                     outcome = 'Symptomatic', 
                                     baseline_sympt_eff = baseline_sympt_effs, 
                                     baseline_eff = baseline_sympt_effs
                                     )
plot_neut_eff_struct_improvement_severe = data.frame(baseline_neut = baseline_neuts, 
                                     outcome = 'Severe', 
                                     baseline_sympt_eff = baseline_sympt_effs, 
                                     baseline_eff = baseline_severe_effs
)

makeNewValues = T
if(makeNewValues){
  add_cols = 3
  if(do_after_sixm){
    add_cols = add_cols+2
  } 
  if(do_after_other_weeks){
    add_cols = add_cols+6
  }
  plot_neut_eff_struct_improvement = data.frame(matrix(nrow=0,ncol = ncol(plot_neut_eff_struct_improvement_sympt)+add_cols))
  colnames(plot_neut_eff_struct_improvement)<-colnames(c(plot_neut_eff_struct_improvement_sympt,'from','GMRdiff','value'))


  use_boosts = boosts
  use_rises = rises

  use_boosts = ancestral_boost
  use_rises = vsv_rises

  for (boost in use_boosts){
    for (rise in use_rises){
      for (this_outcome in outcomes){
        # Set up a temporary structure for this boost, rise and outcome
        if (this_outcome == outcomes[1]){
          plot_neut_eff_struct_improvement_temp = plot_neut_eff_struct_improvement_sympt
        } else {
          plot_neut_eff_struct_improvement_temp = plot_neut_eff_struct_improvement_severe
        }
        plot_neut_eff_struct_improvement_temp$from = boost
        plot_neut_eff_struct_improvement_temp$GMRdiff = rise
        # Work out the difference in efficacy
        plot_neut_eff_struct_improvement_temp$value = filter(plot_neut_eff_struct, GMR == boost*rise, outcome == this_outcome)$new_eff -
                                                            filter(plot_neut_eff_struct, GMR == boost, outcome == this_outcome)$new_eff
        if(do_after_sixm){
          plot_neut_eff_struct_improvement_temp$value6m = filter(plot_neut_eff_struct, GMR == boost*rise, outcome == this_outcome)$new_eff_after_6m - 
                                                            filter(plot_neut_eff_struct, GMR == boost, outcome == this_outcome)$new_eff_after_6m
        }
        if (do_after_other_weeks){
          plot_neut_eff_struct_improvement_temp$value4w = filter(plot_neut_eff_struct, GMR == boost*rise, outcome == this_outcome)$new_eff_after_4w - 
            filter(plot_neut_eff_struct, GMR == boost, outcome == this_outcome)$new_eff_after_4w
          plot_neut_eff_struct_improvement_temp$value8w = filter(plot_neut_eff_struct, GMR == boost*rise, outcome == this_outcome)$new_eff_after_8w - 
            filter(plot_neut_eff_struct, GMR == boost, outcome == this_outcome)$new_eff_after_8w
          plot_neut_eff_struct_improvement_temp$value12w = filter(plot_neut_eff_struct, GMR == boost*rise, outcome == this_outcome)$new_eff_after_12w - 
            filter(plot_neut_eff_struct, GMR == boost, outcome == this_outcome)$new_eff_after_12w
        }
        if(include_average_over_6m){
          plot_neut_eff_struct_improvement_temp$value6mavg = filter(plot_neut_eff_struct, GMR == boost*rise, outcome == this_outcome)$avg_eff_over_6m -
                                                            filter(plot_neut_eff_struct, GMR == boost, outcome == this_outcome)$avg_eff_over_6m
        }
        if (include_average_over_other_weeks){
          plot_neut_eff_struct_improvement_temp$value4wavg = filter(plot_neut_eff_struct, GMR == boost*rise, outcome == this_outcome)$avg_eff_over_4w -
            filter(plot_neut_eff_struct, GMR == boost, outcome == this_outcome)$avg_eff_over_4w
          plot_neut_eff_struct_improvement_temp$value8wavg = filter(plot_neut_eff_struct, GMR == boost*rise, outcome == this_outcome)$avg_eff_over_8w -
            filter(plot_neut_eff_struct, GMR == boost, outcome == this_outcome)$avg_eff_over_8w
          plot_neut_eff_struct_improvement_temp$value12wavg = filter(plot_neut_eff_struct, GMR == boost*rise, outcome == this_outcome)$avg_eff_over_12w -
            filter(plot_neut_eff_struct, GMR == boost, outcome == this_outcome)$avg_eff_over_12w
          
        }
        plot_neut_eff_struct_improvement = rbindlist(list(plot_neut_eff_struct_improvement,plot_neut_eff_struct_improvement_temp))      
      }
    }
  }
  plot_neut_eff_struct_improvement$outcome = factor(plot_neut_eff_struct_improvement$outcome, levels = outcomes)
}
