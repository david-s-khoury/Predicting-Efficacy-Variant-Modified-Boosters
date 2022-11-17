
# Do the fold rise / improvement plots

y_improv_label = 'Benefit of Variant Booster (fold)'
x_improv_label = 'Strain tested in-vitro'
x_matching_improv_label = 'Vaccine vs Variant Immunogen'
x_valency_improv_label = 'Vaccine composition'
improv_title = 'Improvement for Specific Vaccine'
improv_breaks = c(.5,1,2,4,8)

# First plot everything together
vsv_improvement_plot_all_data = ggplot(vsv_improvement_data_new, aes(x=1,y=VSVimprovement))+
  geom_boxplot(fill='grey95',outlier.shape = NA)+
  geom_point(aes(colour=Variant),size=3,position = position_jitter(width=.1))+
  theme_classic()+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  #stat_summary(aes(label=round(after_stat(y), 2), fontface='bold'),fun=geomean, geom="text",vjust=-1,position = position_dodge(0.9)) +
  scale_colour_manual(name='Variant Tested',values = variant_colours)+
  labs(y=y_improv_label,title=paste0('Non-Ancestral variants: ',improv_title),x='All Boosters Combined')
#print(vsv_improvement_plot_all_data)
ggsave(paste0(dir$plots,'ImprovementFromVSV_All.pdf'),vsv_improvement_plot_all_data, width=6,height=5)

# Then split by variant tested - either ancestral or other variants
vsv_improvement_plot_all_data_by_variant = ggplot(vsv_improvement_data_new, aes(x=VariantGroup,y=log10(VSVimprovement),fill=VariantGroup))+
  geom_boxplot(aes(alpha=.5),outlier.shape = NA)+
  geom_point(aes(colour=Variant, shape=FirstAuthor),size=3,position = position_jitter(width=.1))+
  theme_classic()+
  stat_summary(aes(label=round(10^after_stat(y), 2), fontface='bold'),fun=mean, geom="text",vjust=-1,position = position_dodge(.9)) +
  stat_compare_means(comparisons = list(c('Ancestral','Variant')),paired=F,label=c("p.signif"),method='t.test')+
  stat_compare_means(paired=F,label.y = .6,method='t.test')+
  scale_shape_manual(name = 'Study',values=study_shapes, labels = study_shape_labels)+
  scale_colour_manual(name='Variant Tested',values = variant_colours)+
  scale_y_continuous(breaks = log10(improv_breaks),labels =improv_breaks) +
  #scale_y_log10()+
  scale_fill_manual(values=c('Ancestral'='darkslategray1','Variant'='darkseagreen1'), guide='none')+
  labs(y=y_improv_label,title=paste0('All variants: ',improv_title),x=x_improv_label)
print(vsv_improvement_plot_all_data_by_variant)
ggsave(paste0(dir$plots,'ImprovementFromVSV_All_byVariant.pdf'),vsv_improvement_plot_all_data_by_variant, width=6,height=5)

# Make improvement plots split by matched / unmatched
vsv_improvement_plot_split_by_matching = ggplot(filter(vsv_improvement_data_new,Variant!='Ancestral'), aes(x=Matching,y=log10(VSVimprovement),fill=Matching))+
  geom_boxplot(aes(alpha=.6),outlier.shape = NA)+
  geom_point(aes(x=Matching,colour=Variant, shape=FirstAuthor),size=3,position = position_jitter(width=.1))+
  stat_summary(aes(label=round(10^after_stat(y), 2), fontface='bold'),fun=mean, geom="text",vjust=-1,position = position_dodge(0.9)) +
  stat_compare_means(comparisons = list(matching_text),paired=F,method='t.test',label=c("p.signif"))+
  stat_compare_means(paired=F,method='t.test',label.y=.63)+
  theme_classic()+
  #scale_x_discrete(labels=matching_text)+
  #scale_y_continuous(limits = c(.5,3), breaks = seq(.5,3,.5))+
  scale_fill_discrete(name = 'Variant / Booster Match', guide='none')+
  #scale_shape_discrete(name='Vaccine Immunogen vs Variant')+
  scale_y_continuous(breaks = log10(improv_breaks),labels =improv_breaks) +
  scale_shape_manual(name = 'Study',values=study_shapes, labels = study_shape_labels)+
  scale_colour_manual(name='Variant Tested in vitro',values = variant_colours)+
  labs(y=y_improv_label,title=paste0('Non-Ancestral variants: ',improv_title),x=x_matching_improv_label)
#print((vsv_improvement_plot_split_by_matching)
ggsave(paste0(dir$plots,'ImprovementFromVSV_split_by_matching.pdf'),vsv_improvement_plot_split_by_matching, width=6,height=5)

# Redo these plots with specific immunogens only
do_specific_variants = F
if (do_specific_variants){
# Omicron
vsv_improvement_plot_split_by_matching_omicron = vsv_improvement_plot_split_by_matching+
  labs(title=(paste0('Omicron: ',improv_title)))
vsv_improvement_plot_split_by_matching_omicron$data = filter(vsv_improvement_plot_split_by_matching_omicron$data, VSVBoosterPart1=='Omicron BA.1'| VSVBoosterPart2 == 'Omicron BA.1')
ggsave(paste0(dir$plots,'ImprovementFromVSV_split_by_matching_omicron.pdf'),vsv_improvement_plot_split_by_matching_omicron, width=6,height=5)

# Beta
vsv_improvement_plot_split_by_matching_beta = vsv_improvement_plot_split_by_matching+
  labs(title=(paste0('Beta: ',improv_title)))
vsv_improvement_plot_split_by_matching_beta$data = filter(vsv_improvement_plot_split_by_matching_beta$data, VSVBoosterPart1=='Beta'| VSVBoosterPart2 == 'Beta')
ggsave(paste0(dir$plots,'ImprovementFromVSV_split_by_matching_beta.pdf'),vsv_improvement_plot_split_by_matching_beta, width=6,height=5)

#Delta
vsv_improvement_plot_split_by_matching_delta = vsv_improvement_plot_split_by_matching+
  labs(title=(paste0('Delta: ',improv_title)))
vsv_improvement_plot_split_by_matching_delta$data = filter(vsv_improvement_plot_split_by_matching_delta$data, VSVBoosterPart1=='Delta'| VSVBoosterPart2 == 'Delta')
ggsave(paste0(dir$plots,'ImprovementFromVSV_split_by_matching_delta.pdf'),vsv_improvement_plot_split_by_matching_delta, width=6,height=5)
}
# Now do the matched / unmatched plots by valency and then vaccine manufacturer
vsv_improvement_plot_pure_or_mix1 = vsv_improvement_plot_split_by_matching +
  facet_wrap(~Valency)
#ggsave(paste0(dir$plots,'ImprovementFromVSV_split_by_matching_ByValency.pdf'),vsv_improvement_plot_pure_or_mix1, width=6,height=5)

vsv_improvement_plot_by_manufacturer1 = vsv_improvement_plot_split_by_matching+facet_wrap(~BoosterManufacturer)
#ggsave(paste0(dir$plots,'ImprovementFromVSV_split_by_matching_ByManufacturer.pdf'),vsv_improvement_plot_by_manufacturer1, width=8,height=5)


# Now split by valency instead of matching
vsv_improvement_plot_split_by_valency = ggplot(filter(vsv_improvement_data_new,Variant!='Ancestral'), aes(x=Valency,y=log10(VSVimprovement),fill=Valency))+
  geom_boxplot(outlier.shape = NA, alpha=.6)+
  geom_point(aes(x=Valency,colour=Variant, shape = FirstAuthor),size=3,position = position_jitter(width=.1))+
  stat_summary(aes(label=round(10^after_stat(y), 2), fontface='bold'),fun=mean, geom="text",vjust=-1,position = position_dodge(0.9)) +
  stat_compare_means(comparisons = list(c('Bivalent','Monovalent')),paired=F,method='t.test',label=c("p.signif"))+
  stat_compare_means(paired=F,method='t.test',label.y = .75)+
  theme_classic()+
  scale_shape_manual(name = 'Reference',values=study_shapes, labels = study_shape_labels_2)+
  scale_fill_manual(name = 'Vaccine composition', guide='none',values=valency_colours)+
  scale_colour_manual(name='Variant Tested',values = variant_colours)+
  scale_y_continuous(breaks = log10(improv_breaks),labels =improv_breaks) +
  labs(y=y_improv_label,x=x_valency_improv_label)
#print(vsv_improvement_plot_split_by_valency)
ggsave(paste0(dir$plots,'ImprovementFromVSV_ValencySplit.pdf'),vsv_improvement_plot_split_by_valency, width=6,height=5)

# And redo the above also dividing by matched / unmatched and by manufacturer
vsv_improvement_plot_pure_or_mix2 = vsv_improvement_plot_split_by_valency +
  facet_wrap(~Matching)
#print(vsv_improvement_plot_pure_or_mix2)
ggsave(paste0(dir$plots,'ImprovementFromVSV_ValencySplitByMatching.pdf'),vsv_improvement_plot_pure_or_mix2, width=6,height=5)

vsv_improvement_plot_by_manufacturer2 = vsv_improvement_plot_split_by_valency+facet_wrap(~BoosterManufacturer)
#print(vsv_improvement_plot_by_manufacturer2)
#ggsave(paste0(dir$plots,'ImprovementFromVSV_ValencySplitByManufacturer.pdf'),vsv_improvement_plot_by_manufacturer2, width=8,height=5)

# vsv_improvement_plot_ancestral = vsv_improvement_plot+
#   labs(y='Fold Increase in Neut Titre with Variant Specific Vaccine',title='Improvement for Variant Vaccine (Including Ancestral Variants)')
# vsv_improvement_plot_ancestral$data = vsv_improvement_data_new
# print(vsv_improvement_plot_ancestral)
# vsv_improvement_plot_ancestral_pure_or_mix = vsv_improvement_plot_ancestral +
#   facet_wrap(~BoosterGroup2)
# print(vsv_improvement_plot_ancestral_pure_or_mix)
# 
# vsv_improvement_plot_by_manufacturer = vsv_improvement_plot_ancestral+facet_wrap(~BoosterManufacturer)
# print(vsv_improvement_plot_by_manufacturer)
# 

# Now split by prior Status / Boosting
vsv_improvement_plot_split_by_priorstatus = ggplot(filter(vsv_improvement_data_new,Variant!='Ancestral',PriorStatusGroup!='Mixed'), aes(x=PriorStatusGroup,y=log10(VSVimprovement),fill=PriorStatusGroup))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(x=PriorStatusGroup,colour=Variant, shape = FirstAuthor),size=3,position = position_jitter(width=.1))+
  stat_summary(aes(label=round(10^after_stat(y), 2), fontface='bold'),fun=mean, geom="text",vjust=-1,position = position_dodge(0.9)) +
  stat_compare_means(comparisons = list(c('Infected','Uninfected')),paired=F,method='t.test',label=c("p.signif"))+
  stat_compare_means(paired=F,method='t.test',label.y = .75)+
  theme_classic()+
  scale_shape_manual(name = 'Reference',values=study_shapes, labels = study_shape_labels_2)+
  scale_fill_manual(name = 'PriorStatus', guide='none',values=c('lightpink','lightblue'))+ #'thistle',
  scale_colour_manual(name='Variant Tested',values = variant_colours)+
  scale_x_discrete(labels=c('None','Previously Infected'))+
  scale_y_continuous(breaks = log10(improv_breaks),labels =improv_breaks) +
  labs(y=y_improv_label,x='Previous SARS-COV-2 Infection')
#print(vsv_improvement_plot_split_by_priorstatus)
ggsave(paste0(dir$plots,'ImprovementFromVSV_PriorStatusSplit.pdf'),vsv_improvement_plot_split_by_priorstatus, width=6,height=5)

vsv_improvement_plot_split_by_priordoses = ggplot(filter(vsv_improvement_data_new,Variant!='Ancestral'), aes(x=PriorDoses,y=log10(VSVimprovement),fill=PriorDoses))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(x=PriorDoses,colour=Variant, shape = FirstAuthor),size=3,position = position_jitter(width=.1))+
  stat_summary(aes(label=round(10^after_stat(y), 2), fontface='bold'),fun=mean, geom="text",vjust=-1,position = position_dodge(0.9)) +
  stat_compare_means(comparisons = list(c(1,2)),paired=F,method='t.test',label=c("p.signif"))+
  stat_compare_means(paired=F,method='t.test',label.y = .75)+
  theme_classic()+
  scale_shape_manual(name = 'Reference',values=study_shapes, labels = study_shape_labels_2)+
  scale_fill_manual(name = 'PriorStatus', guide='none',values=c('thistle','bisque'))+
  scale_colour_manual(name='Variant Tested',values = variant_colours)+
  scale_y_continuous(breaks = log10(improv_breaks),labels =improv_breaks) +
  #scale_x_discrete(labels=c('None','Prior Infection'))
  labs(y=y_improv_label,x='Previous Vaccination History' )
#print((vsv_improvement_plot_split_by_priordoses)
ggsave(paste0(dir$plots,'ImprovementFromVSV_PriorDosesSplit.pdf'),vsv_improvement_plot_split_by_priordoses, width=6,height=5)



BoosterTypeNames = c('Ancestral + Omicron BA.1', 'Ancestral + Omicron BA.5','Ancestral + Beta', 'Omicron BA.1', 'Beta', 'Beta + Omicron BA.1', 'Beta + Delta', 'Delta + Omicron BA.1')
names(BoosterTypeNames) = c('Ancestral_BA1', 'Ancestral_BA5', 'Ancestral_Beta', 'BA1', 'Beta', 'Beta_BA1', 'Beta_Delta', 'Delta_BA1')
vsv_improvement_plot_all_data_by_vaccine = ggplot(vsv_improvement_data_new, aes(x=BoosterType,y=log10(VSVimprovement), fill=VariantGroup))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(colour=Variant, shape=FirstAuthor),size=3,position = position_jitter(width=.1))+
  theme_classic()+
  stat_summary(aes(label=round(10^after_stat(y), 2), fontface='bold'),fun=mean, geom="text",vjust=-1,position = position_dodge(.9)) +
  stat_compare_means(paired=F,label.y = .8,method='anova')+
  scale_shape_manual(name = 'Study',values=study_shapes, labels = study_shape_labels)+
  scale_colour_manual(name='Variant Tested',values = variant_colours)+
  scale_y_continuous(breaks = log10(improv_breaks),labels =improv_breaks) +
  scale_x_discrete(labels =BoosterTypeNames) +
  #scale_y_log10()+
  scale_fill_manual(values=c('Ancestral'='darkslategray1','Variant'='darkseagreen1'), guide='none')+
  labs(y=y_improv_label,title=paste0('Improvement By Booster Type'),x='Booster Vaccine Composition')+
  facet_wrap(~VariantGroup, ncol=1)+
  theme(axis.text.x = element_text(angle=45,hjust=1))
#print(vsv_improvement_plot_all_data_by_vaccine)
ggsave(paste0(dir$plots,'ImprovementFromVSV_All_byBoosterType.pdf'),vsv_improvement_plot_all_data_by_vaccine, width=9,height=8)


