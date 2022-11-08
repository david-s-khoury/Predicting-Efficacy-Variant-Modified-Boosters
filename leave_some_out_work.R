studies = unique(vsvdata$Study)
nstudies = length(studies)
n_leave_out = 1
n_to_use = nstudies-n_leave_out
n_runs = 1000

# extract boosts with ancestral and variant vaccines
anc_boost_log10fc = filter(vsvdata,BoosterType=='Ancestral', !is.na(log10foldchange)) %>%
  select(c('Study','log10foldchange'))

# Extract improvement with ancestral and variant vaccines as well as for matched and un-matched
vsv_improvement_variant = filter(vsv_improvement_data_new,VariantGroup=='Variant') %>% 
  select(c('Study','VSVimprovement'))
vsv_improvement_variant_matched = filter(vsv_improvement_data_new,VariantGroup=='Variant', Matching == matching_text[1]) %>% 
  select(c('Study','VSVimprovement'))
vsv_improvement_variant_nonmatched = filter(vsv_improvement_data_new,VariantGroup=='Variant', Matching == matching_text[2]) %>% 
  select(c('Study','VSVimprovement'))

# Define the ancestral boost and set up the list of possible boosts
ancestral_boost = round(10^mean(anc_boost_log10fc$log10foldchange),1)
# Define the vsv rise and set up the list of rises
vsv_rise = round(geomean(vsv_improvement_variant$VSVimprovement),2)
vsv_rise_matched = round(geomean(vsv_improvement_variant_matched$VSVimprovement),2)
vsv_rise_nonmatched = round(geomean(vsv_improvement_variant_nonmatched$VSVimprovement),2)

calculated_rises = data.frame(AncBoost = ancestral_boost, AvgImprov = vsv_rise, MatchedImprov = vsv_rise_matched, NonMatchedImprov = vsv_rise_nonmatched )

for (i in c(1:n_runs)){
  use_studies = studies[sample(studies, n_to_use, replace = F)]
  
  anc_boost_log10fc_use = filter(anc_boost_log10fc, Study %in% use_studies )$log10foldchange
  ancestral_boost = round(10^mean(anc_boost_log10fc_use),1)
  
  vsv_improvement_variant_use = filter(vsv_improvement_variant, Study %in% use_studies )$VSVimprovement
  vsv_rise = round(geomean(vsv_improvement_variant_use),2)
  
  vsv_improvement_variant_matched_use = filter(vsv_improvement_variant_matched, Study %in% use_studies )$VSVimprovement
  vsv_rise_matched = round(geomean(vsv_improvement_variant_matched_use),2)
  
  vsv_improvement_variant_nonmatched_use = filter(vsv_improvement_variant_nonmatched, Study %in% use_studies )$VSVimprovement
  vsv_rise_nonmatched = round(geomean(vsv_improvement_variant_nonmatched_use),2)
  calculated_rises = rbind(calculated_rises,
                           c(ancestral_boost, vsv_rise, vsv_rise_matched, vsv_rise_nonmatched ))
  
}
calculated_rises_for_plot = melt(calculated_rises[c(2:nrow(calculated_rises)),])
rise_hists = ggplot(calculated_rises_for_plot, aes(x=value))+facet_wrap(~variable, scales = 'free')+geom_histogram()
print(rise_hists)
