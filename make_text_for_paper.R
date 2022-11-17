# This needs to be run after make plots

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


#fold rise text
rounding_level = 1
rounding_level2 = 2

# anc_boost_log10fc_noanc = filter(vsvdata,BoosterType=='Ancestral', !is.na(log10foldchange), Variant != 'Ancestral')$log10foldchange
# vsv_boost_log10fc_noanc = filter(vsvdata,BoosterGroup1=='Variant',!is.na(log10foldchange), Variant != 'Ancestral')$log10foldchange
# anc_boost_log10fc = anc_boost_log10fc_noanc
# vsv_boost_log10fc = vsv_boost_log10fc_noanc

ancestral_variant_testing_data = filter(vsvdata, Variant=='Ancestral')
voc_variant_testing_data = filter(vsvdata, Variant!='Ancestral')

ancestral_variant_testing_improv = filter(vsv_improvement_data_new, VariantGroup=='Ancestral')
voc_variant_testing_improv = filter(vsv_improvement_data_new, VariantGroup=='Variant')


neut_of_ancestral_comparison=t.test(ancestral_variant_testing_data$log10foldchange[ancestral_variant_testing_data$BoosterGroup1=='Ancestral'],
                                    ancestral_variant_testing_data$log10foldchange[ancestral_variant_testing_data$BoosterGroup1=='Variant'],paired=F)

neut_of_voc_comparison=t.test(voc_variant_testing_data$log10foldchange[voc_variant_testing_data$BoosterGroup1=='Ancestral'],
                              voc_variant_testing_data$log10foldchange[voc_variant_testing_data$BoosterGroup1=='Variant'],paired=F)

neut_of_ancestral_improv_comparison=t.test(ancestral_variant_testing_improv$log10VSVimprovement)
neut_of_voc_improv_comparison=t.test(voc_variant_testing_improv$log10VSVimprovement)

vsv_rise = round(geomean(vsv_improvement_variant),2)
vsv_rise_matched = round(geomean(vsv_improvement_variant_matched),2)
vsv_rise_nonmatched = round(geomean(vsv_improvement_variant_nonmatched),2)

vsv_rises_monovalent = filter(vsv_improvement_data_new,VariantGroup=='Variant', Valency=='Monovalent')$VSVimprovement
vsv_rises_bivalent = filter(vsv_improvement_data_new,VariantGroup=='Variant', Valency=='Bivalent')$VSVimprovement
valency_comparison = t.test(log10(vsv_rises_monovalent), log10(vsv_rises_bivalent), paried = F)

vsv_rises_priorinfect = filter(vsv_improvement_data_new,VariantGroup=='Variant', PriorStatusGroup=='Infected')$VSVimprovement
vsv_rises_prioruninf = filter(vsv_improvement_data_new,VariantGroup=='Variant', PriorStatusGroup=='Uninfected')$VSVimprovement
prior_status_comparison = t.test(log10(vsv_rises_prioruninf), log10(vsv_rises_priorinfect), paried = F)

vsv_rises_priorprimary = filter(vsv_improvement_data_new,VariantGroup=='Variant', PriorDoses=='2 Doses')$VSVimprovement
vsv_rises_priorboosted = filter(vsv_improvement_data_new,VariantGroup=='Variant', PriorDoses=='3+ Doses')$VSVimprovement
prior_doses_comparison = t.test(log10(vsv_rises_priorprimary), log10(vsv_rises_priorboosted), paried = F)

para2_text=c()
i=1
para2_text[i] = paste0('Considering in vitro neutralisation for all variants reported in the studies, we found that an ancestral-based vaccine increased neutralisation titres by a mean of ',
              round(ancestral_boost,rounding_level),'-fold compared with pre-booster titres (95%CI ',
              cis(10^anc_boost_log10fc,'l', rounding = rounding_level),'-',
              cis(10^anc_boost_log10fc,'u', rounding = rounding_level),') (Fig 1A). Although variant-modified vaccines did not show an improvement in neutralisation towards ancestral SARS-CoV-2 compared to the ancestral-based vaccines (p= ',
              round(neut_of_ancestral_improv_comparison$p.value,rounding_level2),'), the variant-modified vaccines produced more potent neutralisatio of the variants tested. Considering only neutralisation of SARS-CoV-2 variant strains, we found that variant-modified vaccines on average produced ',
              vsv_rises[1],'-fold [95% CI ',
              cis(vsv_improvement_variant,'l', rounding = rounding_level),'-',cis(vsv_improvement_variant,'u', rounding = rounding_level),
              '] higher titres than the equivalent ancestral-based vaccine (p=',
              round(neut_of_voc_improv_comparison$p.value,13),', fig 1B).')
i=i+1
para2_text[i]  = paste0('Boosting was slightly higher to homologous strains (',
                            round(vsv_rise_matched, rounding_level2),' vs. ',
                            round(vsv_rise_nonmatched, rounding_level2),', p=',
                            round(t.test(log10(vsv_improvement_variant_matched),log10(vsv_improvement_variant_nonmatched))$p.value,3),').  We found no significant differences when stratifying results for monovalent versus bivalent vaccines (p = ',
                            round(valency_comparison$p.value, 2),', see Figure SX).')
i=i+1
para2_text[i]  = paste0('Boosting was higher against homologous (',
                               round(vsv_rise_matched, rounding_level2),
                        ' [95% CI ',
                        cis(vsv_improvement_variant_matched,'l', rounding = rounding_level),
                        '-',
                        cis(vsv_improvement_variant_matched,'u', rounding = rounding_level),
                        '])',
                        ' vs. non-homologous (',
                        round(vsv_rise_nonmatched, rounding_level2),
                        ' [95% CI ',
                        cis(vsv_improvement_variant_nonmatched,'l', rounding = rounding_level),
                        '-',
                        cis(vsv_improvement_variant_nonmatched,'u', rounding = rounding_level),
                        ']) strains (p=',
                        round(t.test(log10(vsv_improvement_variant_matched),log10(vsv_improvement_variant_nonmatched))$p.value,2),
                        '). We found no significant difference in the relative improvement conferred by a variant-modified vaccine over an ancestral-based vaccine after stratifying results by vaccine valency (monovalent versus bivalent, p = ',
                               round(valency_comparison$p.value, 2),'), history of prior infection (p=',
                               round(prior_status_comparison$p.value, 2),') and the number of previous vaccines subjects had received (p = ',
                               round(prior_doses_comparison$p.value, 2),') (see Figs S3,S4).')
print(para2_text)


discussion_text=c()
i=1

discussion_text[i] = paste0('we find that an ancestral-based booster giving an ',
                            ancestral_boost,
                            '-fold boost ',
                            '[95%CI ',
                            cis(10^anc_boost_log10fc,'l', rounding = rounding_level),'-',
                            cis(10^anc_boost_log10fc,'u', rounding = rounding_level),
                            '] ',
                            'in neutralising antibody titres would increase the average protection over a six month period against symptomatic infection from 50% to ',
                            round(100*filter(plot_neut_eff_struct_gmr_use, outcome=='Symptomatic', baseline_sympt_eff==.5, GMR==ancestral_boost)$avg_eff_over_6m,rounding_level),
                            '% [95% CI ',
                            round(filter(quantile.results.combined, baseline_eff==50,eff_type=='sympt.boosted')$lower95,rounding_level),
                            '-',
                            round(filter(quantile.results.combined, baseline_eff==50,eff_type=='sympt.boosted')$upper95,rounding_level),
                            ']. Directly comparing the incremental benefit of a variant-modified booster versus an ancestral-based booster, a variant-modified booster is predicted to provide an additional ',
                            round(100*filter(plot_neut_eff_struct_gmr_use, outcome=='Symptomatic', baseline_sympt_eff==.5, GMR==ancestral_boost*vsv_rise)$avg_eff_over_6m-100*filter(plot_neut_eff_struct_gmr_use, outcome=='Symptomatic', baseline_sympt_eff==.5, GMR==ancestral_boost)$avg_eff_over_6m,rounding_level),
                            ' percentage points ',
                            
                            ' [95% CI ',
                            round(filter(quantile.results.combined, baseline_eff==50,eff_type=='sympt.improved.diff')$lower95,rounding_level),
                            '-',
                            round(filter(quantile.results.combined, baseline_eff==50,eff_type=='sympt.improved.diff')$upper95,rounding_level),
                            ']',
                            ' protection against symptomatic compared to the ancestral-based vaccine, Fig 2A and C).'
                            )
i=i+1

discussion_text[i] = paste0('A population with 50% protection from symptomatic infection is predicted to have ',
    round(100*filter(plot_neut_eff_struct_gmr_use, outcome=='Severe', baseline_sympt_eff==.5, GMR==ancestral_boost)$baseline_eff,rounding_level),
    '% protection from severe COVID-19. Boosting with an ancestral-based booster',
    ' is expected to increase this to an average of ',
    round(100*filter(plot_neut_eff_struct_gmr_use, outcome=='Severe', baseline_sympt_eff==.5, GMR==ancestral_boost)$avg_eff_over_6m,rounding_level),
    
    '% [95% CI ',
    round(filter(quantile.results.combined, baseline_eff==50,eff_type=='severe.boosted')$lower95,rounding_level),
    '-',
    round(filter(quantile.results.combined, baseline_eff==50,eff_type=='severe.boosted')$upper95,rounding_level),
    ']',
    
    
    ' protection from severe COVID-19 over the six month period following boosting, and a variant-modified booster',
    #' producing ',
    #round(vsv_rise,rounding_level),
    #'-fold ',
    #'[95% CI ',
    #cis(vsv_improvement_variant,'l', rounding = rounding_level),'-',cis(vsv_improvement_variant,'u', rounding = rounding_level),
    #']',
    #' higher neutralising antibody titres',
    ' is predicted to provide an additional ',
    round(100*filter(plot_neut_eff_struct_gmr_use, outcome=='Severe', baseline_sympt_eff==.5, GMR==ancestral_boost*vsv_rise)$avg_eff_over_6m-100*filter(plot_neut_eff_struct_gmr_use, outcome=='Severe', baseline_sympt_eff==.5, GMR==ancestral_boost)$avg_eff_over_6m,rounding_level),
    ' percentage points',
    
    ' [95% CI ',
    round(filter(quantile.results.combined, baseline_eff==50,eff_type=='severe.improved.diff')$lower95,rounding_level),
    '-',
    round(filter(quantile.results.combined, baseline_eff==50,eff_type=='severe.improved.diff')$upper95,rounding_level),
    ', p= <.001]',
    
    ' of protection from severe COVID-19 compared to an ancestral-based booster, Fig 2B and D). '
    )
i=i+1

maxdiff = max(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise, outcome=='Symptomatic')$value6mavg)
maxdiff_ind = which(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise, outcome=='Symptomatic')$value6mavg==maxdiff)
maxdiff_sympt_eff = as.numeric(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise, outcome=='Symptomatic')[maxdiff_ind, 'baseline_sympt_eff'])
maxdiff_neuts = as.numeric(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise, outcome=='Symptomatic')[maxdiff_ind, 'baseline_neut'])
maxdiff_severe_eff = as.numeric(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise, outcome=='Severe', baseline_sympt_eff==maxdiff_sympt_eff)$baseline_eff)


discussion_text[i] = paste0('In general, the lower the pre-booster immunity the greater the relative benefit of a variant-modified booster compared to an ancestral-based booster, peaking with an additional ',
                            round(100*maxdiff,rounding_level),
                            ' [95% CI ',
                            round(filter(quantile.results.combined, eff_type == 'sympt.improved.diff',baseline_eff==100*maxdiff_sympt_eff)$lower95,rounding_level),
                            '-',
                            round(filter(quantile.results.combined, eff_type == 'sympt.improved.diff',baseline_eff==100*maxdiff_sympt_eff)$upper95,rounding_level),
                            '] ',
                            'percentage points additional protection from symptomatic infection when the population has only ',
                            round(100*maxdiff_sympt_eff, rounding_level),
                            '%',
                            
                            ' pre-existing protection against symptomatic infection prior to the boost (which corresponds to ',
                            round(100*maxdiff_severe_eff,rounding_level),'%',
                            
                            ' protection against severe COVID-19 (Fig 2D)). At this level of population immunity, a variant modified vaccine would provide an additional ',
                            round(100*filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise, outcome=='Severe', baseline_sympt_eff==maxdiff_sympt_eff)$value6mavg, rounding_level),
                            ' percentage points [95% CI ',
                            round(filter(quantile.results.combined, eff_type == 'severe.improved.diff',baseline_eff==100*maxdiff_sympt_eff)$lower95,rounding_level),
                            '-',
                            round(filter(quantile.results.combined, eff_type == 'severe.improved.diff',baseline_eff==100*maxdiff_sympt_eff)$upper95,rounding_level),
                            ']',
                            ' protection from severe COVID-19 compared to an ancestral-based vaccine.'
                            )
i=i+1

maxdiff.matched = max(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_matched, outcome=='Symptomatic')$value6mavg)
maxdiff_ind.matched = which(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_matched, outcome=='Symptomatic')$value6mavg==maxdiff.matched)
maxdiff_sympt_eff.matched = as.numeric(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_matched, outcome=='Symptomatic')[maxdiff_ind.matched, 'baseline_sympt_eff'])
maxdiff_neuts.matched = as.numeric(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_matched, outcome=='Symptomatic')[maxdiff_ind.matched, 'baseline_neut'])
maxdiff_severe_eff.matched = as.numeric(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_matched, outcome=='Severe', baseline_sympt_eff==maxdiff_sympt_eff.matched)$baseline_eff)

maxdiff.nonmatched = max(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_nonmatched, outcome=='Symptomatic')$value6mavg)
maxdiff_ind.nonmatched = which(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_nonmatched, outcome=='Symptomatic')$value6mavg==maxdiff.nonmatched)
maxdiff_sympt_eff.nonmatched = as.numeric(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_nonmatched, outcome=='Symptomatic')[maxdiff_ind.nonmatched, 'baseline_sympt_eff'])
maxdiff_neuts.nonmatched = as.numeric(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_nonmatched, outcome=='Symptomatic')[maxdiff_ind.nonmatched, 'baseline_neut'])
maxdiff_severe_eff.nonmatched = as.numeric(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise_nonmatched, outcome=='Severe', baseline_sympt_eff==maxdiff_sympt_eff.nonmatched)$baseline_eff)




discussion_text[i] = paste0('We found a maximum benefit of a ',
                            round(100*maxdiff.matched,rounding_level),
                            ' [',
                            round(filter(quantile.results.combined, eff_type == 'sympt.improved.matched.diff',baseline_eff==100*maxdiff_sympt_eff.matched)$lower95,rounding_level),
                            '-',
                            round(filter(quantile.results.combined, eff_type == 'sympt.improved.matched.diff',baseline_eff==100*maxdiff_sympt_eff.matched)$upper95,rounding_level),
                            '] ',
                        
                            ' percentage point improvement in protection from symptomatic infection for boosters against antigenically related strains and ',
                            round(100*maxdiff.nonmatched,rounding_level),
                            ' [',
                            round(filter(quantile.results.combined, eff_type == 'sympt.improved.nonmatched.diff',baseline_eff==100*maxdiff_sympt_eff.nonmatched)$lower95,rounding_level),
                            '-',
                            round(filter(quantile.results.combined, eff_type == 'sympt.improved.nonmatched.diff',baseline_eff==100*maxdiff_sympt_eff.nonmatched)$upper95,rounding_level),
                            '] ',
                        
                            ' percentage point improvement for antigenically distinct boosters (Fig 2C and D). At the level of population immunity corresponding to this maximum improvement in protection against symptomatic disease, protection against severe COVID-19 would be improved by ',
                            round(100*filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise, outcome=='Severe', baseline_sympt_eff==maxdiff_sympt_eff.matched)$value6mavg, rounding_level),
                            ' [',
                            round(filter(quantile.results.combined, eff_type == 'severe.improved.matched.diff',baseline_eff==100*maxdiff_sympt_eff.matched)$lower95,rounding_level),
                            '-',
                            round(filter(quantile.results.combined, eff_type == 'severe.improved.matched.diff',baseline_eff==100*maxdiff_sympt_eff.matched)$upper95,rounding_level),
                            ']',
                            ' and ',
                            round(100*filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise, outcome=='Severe', baseline_sympt_eff==maxdiff_sympt_eff.nonmatched)$value6mavg, rounding_level),
                            ' [',
                            round(filter(quantile.results.combined, eff_type == 'severe.improved.nonmatched.diff',baseline_eff==100*maxdiff_sympt_eff.nonmatched)$lower95,rounding_level),
                            '-',
                            round(filter(quantile.results.combined, eff_type == 'severe.improved.nonmatched.diff',baseline_eff==100*maxdiff_sympt_eff.nonmatched)$upper95,rounding_level),
                            ']',
                            ' percentage points for an antigenically similar / distinct variant booster respectively.')
i=i+1


discussion_text[i]=paste0('We found that for a population with 50% pre-boost protection from symptomatic SARS-CoV-2, antigenic matching of the booster is predicted to provide and additional ',
                          round(100*filter(plot_neut_eff_struct_gmr_use, outcome=='Symptomatic', baseline_sympt_eff==.5, GMR==ancestral_boost*vsv_rise_matched)$avg_eff_over_6m -
                             100*filter(plot_neut_eff_struct_gmr_use, outcome=='Symptomatic', baseline_sympt_eff==.5, GMR==ancestral_boost*vsv_rise_nonmatched)$avg_eff_over_6m,rounding_level),
                          ' [95% CI ',
                          round(filter(quantile.results.combined, baseline_eff==50,eff_type=='sympt.matchedvsnonmatched.boost.diff')$lower95,rounding_level),
                          '-',
                          round(filter(quantile.results.combined, baseline_eff==50,eff_type=='sympt.matchedvsnonmatched.boost.diff')$upper95,rounding_level),
                          '] and ',
                          round(100*filter(plot_neut_eff_struct_gmr_use, outcome=='Severe', baseline_sympt_eff==.5, GMR==ancestral_boost*vsv_rise_matched)$avg_eff_over_6m-
                              100*filter(plot_neut_eff_struct_gmr_use, outcome=='Severe', baseline_sympt_eff==.5, GMR==ancestral_boost*vsv_rise_nonmatched)$avg_eff_over_6m,rounding_level),
                          ' [95% CI ',
                          round(filter(quantile.results.combined, baseline_eff==50,eff_type=='severe.matchedvsnonmatched.boost.diff')$lower95,rounding_level),
                          '-',
                          round(filter(quantile.results.combined, baseline_eff==50,eff_type=='severe.matchedvsnonmatched.boost.diff')$upper95,rounding_level),
                          '] percentage points protection from symptomatic and severe COVID-19 respectively compared to an antigenically unmatched booster'
                          )

print(discussion_text)                    
                         
                         