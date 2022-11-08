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
vsv_rises_priorboosted = filter(vsv_improvement_data_new,VariantGroup=='Variant', PriorDoses=='3 Doses')$VSVimprovement
prior_doses_comparison = t.test(log10(vsv_rises_priorprimary), log10(vsv_rises_priorboosted), paried = F)


para2_text = paste0('Considering in vitro neutralisation for all variants reported in the studies, we found that an ancestral-based vaccine increased neutralisation titres by a mean of ',
              round(ancestral_boost,rounding_level),'-fold compared with pre-booster titres (95%CI ',
              cis(10^anc_boost_log10fc,'l', rounding = rounding_level),'-',
              cis(10^anc_boost_log10fc,'u', rounding = rounding_level),') (Fig 1A). Although variant-modified vaccines did not show an improvement in neutralisation towards ancestral SARS-CoV-2 compared to the ancestral-based vaccines (p= ',
              round(neut_of_ancestral_improv_comparison$p.value,rounding_level2),'), the variant-modified vaccines produced more potent neutralisatio of the variants tested. Considering only neutralisation of SARS-CoV-2 variant strains, we found that variant-modified vaccines on average produced ',
              vsv_rises[1],'-fold [95% CI ',
              cis(vsv_improvement_variant,'l', rounding = rounding_level),'-',cis(vsv_improvement_variant,'u', rounding = rounding_level),
              '] higher titres than the equivalent ancestral-based vaccine (p=',
              round(neut_of_voc_improv_comparison$p.value,13),', fig 1B).')

para2_bottom_text  = paste0('Boosting was slightly higher to homologous strains (',
                            round(vsv_rise_matched, rounding_level2),' vs. ',
                            round(vsv_rise_nonmatched, rounding_level2),', p=',
                            round(t.test(log10(vsv_improvement_variant_matched),log10(vsv_improvement_variant_nonmatched))$p.value,5),').  We found no significant differences when stratifying results for monovalent versus bivalent vaccines (p = ',
                            round(valency_comparison$p.value, 2),', see Figure SX).')

para2_bottom_text_v2  = paste0('Boosting was slightly higher to homologous antigens (',
                               round(vsv_rise_matched, rounding_level2),' vs. ',
                               round(vsv_rise_nonmatched, rounding_level2),', p=',
                               round(t.test(log10(vsv_improvement_variant_matched),log10(vsv_improvement_variant_nonmatched))$p.value,5),'). We found the relative improvement conferred by a variant-modified vaccine over an ancestral-based vaccine remained unchanged after stratifying results by vaccine valency (monovalent versus bivalent, p = ',
                               round(valency_comparison$p.value, 2),'), history of prior infection (p=',
                               round(prior_status_comparison$p.value, 2),') and the number of previous vaccines subjects had received (p = ',
                               round(prior_doses_comparison$p.value, 2),') (see Fig S1-3).')
print(para2_text)
print(para2_bottom_text)
print(para2_bottom_text_v2)