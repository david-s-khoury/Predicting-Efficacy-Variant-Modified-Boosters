# Now plot the fold changes in neut titre after boosting
# We will do this in a number of ways, including by vaccine type

# ploting variables / constants
fc_breaks=c(1,2,5,10,20,50,100)
y_label_anc_fold_change='Increase in Neutralisation Titer'
x_label_anc_fold_change=''

# This is the basic fold change plot - which shows the fold change in neuts after boosting.
fc_plot=ggplot(vsvdata, aes(x=BoosterType,y=log10foldchange))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(mapping=aes(colour=Variant,shape=FirstAuthor),size=3,position = position_jitter(width=.1))+
  stat_summary(aes(label=round(10^after_stat(y), 1),fontface='bold'),fun=mean, geom="text", vjust=-1,position = position_dodge(0.9)) +
  #scale_y_continuous(breaks = log10(drop_breaks),labels = drop_breaks)+
  scale_y_continuous(breaks=log10(fc_breaks),labels=fc_breaks)+
  scale_shape_manual(name = 'Study',values=study_shapes, labels = study_shape_labels)+
  scale_colour_manual(name='Variant Tested in vitro',values = variant_colours, guide='none')+
  #scale_fill_manual(values = drop_colours)+
  theme_classic()+
  theme(axis.text.x = element_text(angle=30,vjust=1,hjust=1))+
  labs(title='Fold change after boosting',y='Fold Change',x='Booster Composition')
#print(fc_plot)
ggsave(paste0(dir$plots,'FoldChangeFromVaccines_ByVaccine.pdf'),fc_plot, width=8,height=5)

# This is the plot we are going to use in the manuscript - only ancestral boosting.
fc_plot_use = fc_plot
fc_plot_use$data = filter(fc_plot$data,BoosterType == 'Ancestral')
fc_plot_use$mapping$x = quo(`BoosterGroup1`)
fc_plot_use = fc_plot_use+
  scale_x_discrete(labels = 'Ancestral Booster')+
  theme(axis.text.x = element_text(angle=0,vjust=1,hjust=.5))+
  labs(y=y_label_anc_fold_change,x=x_label_anc_fold_change)
#print(fc_plot_use)
ggsave(paste0(dir$plots,'FoldChangeFromVaccines_AncOnly.pdf'),fc_plot_use, width=5,height=5)


# This is the new supplementary plot - only ancestral boosting.
fc_plot_supp = fc_plot
fc_plot_supp$data = fc_plot_supp$data %>%
  mutate(VariantGroup = ifelse(Variant=='Ancestral','Ancestral','Non-Ancestral Variant'),
         BoosterGroup1 = ifelse(BoosterGroup1=='Ancestral','Ancestral Booster','Variant Booster'))
fc_plot_supp$mapping$x = quo(`VariantGroup`)
fc_plot_supp = fc_plot_supp+
  #scale_x_discrete(labels = 'Booster Composition')+
  theme(axis.text.x = element_text(angle=30,vjust=1,hjust=1))+
  labs(y=y_label_anc_fold_change,x='Variant Tested in vitro')+
  facet_wrap(~BoosterGroup1)
#print(fc_plot_supp)
ggsave(paste0(dir$plots,'FoldChangeFromVaccines_AllVaccines.pdf'),fc_plot_supp, width=5,height=5)



vsv_new_data = vsvdata %>% 
  reshape2::dcast(formula=cast_formula,
                  value.var = 'log10neutsafterboost')
vsv_new_data$index = rownames(vsv_new_data)
vsv_new_data = vsv_new_data %>% 
  reshape2::melt(id.vars = c("Study", "ComparisonGroup",  "FirstAuthor","PriorStatusGroup","PriorDoses","Variant","index"))%>%
  filter(!is.na(value)) %>%
  mutate(BoosterType = ifelse(variable!='Ancestral', 'Variant Booster','Ancestral Booster'))

vsv_new_data_merged = merge(filter(vsv_new_data, BoosterType=='Ancestral Booster'),
                            filter(vsv_new_data, BoosterType=='Variant Booster'),
                            by=c("Study", "ComparisonGroup",  "FirstAuthor","PriorStatusGroup","PriorDoses","Variant","index"),
                            all=T) %>%
  rename(BoosterComposition1 = variable.x,
         BoosterComposition2 = variable.y,
         log10fc1 = value.x,
         log10fc2 = value.y) %>%
  #filter(!is.na(log10fc1), !is.na(log10fc2)) %>% 
  mutate(
    VariantGroup = ifelse(Variant=='Ancestral','Ancestral Variant','Non-Ancestral Variant'),
  )

vsv_new_data_merged_paired = vsv_new_data_merged %>% filter(!is.na(log10fc1), !is.na(log10fc2))
vsv_new_data_merged_all = vsv_new_data_merged

with(filter(vsv_new_data_merged_paired, Variant=='Ancestral', !is.na(log10fc1), !is.na(log10fc2)), t.test(log10fc1,log10fc2, paired=T))$p.value
with(filter(vsv_new_data_merged_paired, Variant!='Ancestral', !is.na(log10fc1), !is.na(log10fc2)), t.test(log10fc1,log10fc2, paired=T))$p.value

vsv_new_data_merged_melted_paired = vsv_new_data_merged_paired %>% 
  reshape2::melt(measure.vars = c('log10fc1','log10fc2')) %>%
  mutate(BoosterType = ifelse(variable=='log10fc1',BoosterType.x,BoosterType.y),
         log10foldchange = value) %>%
  select(-c('BoosterComposition2'))%>%
  distinct()

vsv_new_data_merged_melted_all = vsv_new_data_merged_all %>% 
  reshape2::melt(measure.vars = c('log10fc1','log10fc2')) %>%
  mutate(BoosterType = ifelse(variable=='log10fc1',BoosterType.x,BoosterType.y),
         log10foldchange = value) %>%
  select(-c('BoosterComposition2'))%>%
  distinct()

vsv_new_data_merged_melted_all_forlines = vsv_new_data_merged_all %>% 
  mutate(joiner = rownames(vsv_new_data_merged_all) )%>%
  reshape2::melt(measure.vars = c('log10fc1','log10fc2')) %>%
  mutate(BoosterType = ifelse(variable=='log10fc1',BoosterType.x,BoosterType.y),
         log10foldchange = value) %>%
  select(-c('BoosterComposition2'))%>%
  distinct()

fc_plot_supp_gmt = ggplot(filter(vsv_new_data_merged_melted_paired, !is.na(BoosterType), VariantGroup=='Ancestral Variant'),
                       aes(x=BoosterType,y=log10foldchange))+
  geom_boxplot(outlier.shape=NA)+
  geom_point(mapping=aes(colour=Variant,shape=FirstAuthor),size=3, position = position_jitter(0.2))+
  scale_shape_manual(name = 'Study',values=study_shapes, labels = study_shape_labels)+
  scale_colour_manual(name='Variant Tested in vitro',values = variant_colours, guide='none')+
  scale_y_continuous(breaks = c(3,log10(3000),4,log10(30000)), labels=c('1,000','3,000','10,000','30,000'))+
  theme_classic()+
  facet_adjust+
  theme(axis.text.x = element_text(angle=30,vjust=1,hjust=1))+
  labs(title='GMT of Neutralisation Titers Agaist the Ancestral Variant',y='GMT of Neutralisaiton Titers after Boosting',x='Booster Composition')
#print(fc_plot_supp_gmt)  
ggsave(paste0(dir$manuscript_plots,'Figure_Response_BoosterComposition.pdf'), fc_plot_supp_gmt, width=7,height=6)



fc_plot_supp2 = ggplot(filter(vsv_new_data_merged_melted_all, !is.na(BoosterType)),
                       aes(x=BoosterType,y=log10foldchange))+
  geom_boxplot(outlier.shape=NA)+
  scale_y_continuous(breaks=log10(fc_breaks),labels=fc_breaks)+
  scale_shape_manual(name = 'Study',values=study_shapes, labels = study_shape_labels)+
  scale_colour_manual(name='Variant Tested in vitro',values = variant_colours, guide='none')+
  theme_classic()+
  theme(axis.text.x = element_text(angle=30,vjust=1,hjust=1))+
  labs(title='Fold change after boosting',y='Fold Change',x='Booster Composition')+
  facet_wrap(~VariantGroup) #+
  #stat_compare_means(method='t.test',comparisons = list(c(1,2)))
#print(fc_plot_supp2)

fc_plot_supp2_alldata = fc_plot_supp2+geom_point(mapping=aes(colour=Variant,shape=FirstAuthor),size=3, position = position_dodge(0.2))
#print(fc_plot_supp2_alldata)

fc_plot_supp2_paireddata = fc_plot_supp2_alldata
fc_plot_supp2_paireddata$data = vsv_new_data_merged_melted_paired 
#print(fc_plot_supp2_paireddata)

fc_plot_supp2_alldata_withlines = fc_plot_supp2
fc_plot_supp2_alldata_withlines$data = filter(vsv_new_data_merged_melted_all_forlines, !is.na(BoosterType))
fc_plot_supp2_alldata_withlines = fc_plot_supp2_alldata_withlines +
  geom_point(mapping=aes(colour=Variant,shape=FirstAuthor, group=joiner),size=3, position = position_dodge(0.2))+
  geom_line(aes(group=joiner, colour=Variant), position = position_dodge(0.2))+
  stat_summary(aes(label=round(10^after_stat(y), 1),fontface='bold'),fun=mean, geom="text", vjust=-1,position = position_dodge(0.9))
  
#print(fc_plot_supp2_alldata_withlines)



fc_plot_1 = fc_plot
#fc_plot_1$data = filter(fc_plot$data,Variant != 'Ancestral')
fc_plot_1$mapping$x = quo(`BoosterGroup1`)
fc_plot_1 = fc_plot_1+
  stat_compare_means(comparisons = list(c('Ancestral','Variant')),paired=F,method='t.test',label=c("p.signif"))+
  stat_compare_means(paired=F,method='t.test',label.y = log10(100))+
  theme(axis.text.x = element_text(angle=0,vjust=1,hjust=.5))+
  labs(title='Fold change after boosting',y='Fold Change',x='Ancestral or Variant Vaccine')
ggsave(paste0(dir$plots,'FoldChangeFromVaccines_AncvsVar.pdf'),fc_plot_1, width=5,height=5)


fc_plot_1 = fc_plot
fc_plot_1$mapping$x = quo(`BoosterGroup1`)
fc_plot_1 = fc_plot_1+
  stat_compare_means(comparisons = list(c('Ancestral','Variant')),paired=F,method='t.test',label=c("p.signif"))+
  stat_compare_means(paired=F,method='t.test',label.y = log10(100))+
  theme(axis.text.x = element_text(angle=0,hjust=.5))+
  labs(title='Fold change after boosting',y='Fold Change',x='Ancestral or Variant Vaccine')
ggsave(paste0(dir$plots,'FoldChangeFromVaccines_AncvsVar.pdf'),fc_plot_1, width=6,height=5)

fc_plot_2 = fc_plot
fc_plot_2$mapping$x = quo(`BoosterGroup2`)
fc_plot_2 = fc_plot_2+
  stat_compare_means(comparisons = list(c('Ancestral','Bivalent'),c('Ancestral','Monovalent'),c('Bivalent','Monovalent')),paired=F,method='t.test',label=c("p.signif"))+
  stat_compare_means(paired=F,method='anova',label.y = log10(100))+
  theme(axis.text.x = element_text(angle=0,hjust=.5))+
  labs(title='Fold change after boosting',y='Fold Change',x='Booster Type')
ggsave(paste0(dir$plots,'FoldChangeFromVaccines_AncvsValency.pdf'),fc_plot_2, width=6,height=5)

fc_plot_3 = fc_plot
fc_plot_3$mapping$x = quo(`PriorStatusGroup`)
fc_plot_3 = fc_plot_3+
  stat_compare_means(comparisons = list(c('Infected','Uninfected')),paired=F,method='t.test',label=c("p.signif"))+
  stat_compare_means(paired=F,method='anova',label.y = log10(100))+
  theme(axis.text.x = element_text(angle=0,hjust=.5))+
  labs(title='Fold change after boosting',y='Fold Change',x='Booster Type')
ggsave(paste0(dir$plots,'FoldChangeFromVaccines_PriorStatus.pdf'),fc_plot_3, width=6,height=5)

fc_plot_4 = fc_plot
fc_plot_4$mapping$x = quo(`PriorDoses`)
fc_plot_4 = fc_plot_4+
  stat_compare_means(comparisons = list(c('Primary Course','Boosted')),paired=F,method='t.test',label=c("p.signif"))+
  stat_compare_means(paired=F,method='anova',label.y = log10(100))+
  theme(axis.text.x = element_text(angle=0,hjust=.5))+
  labs(title='Fold change after boosting',y='Fold Change',x='Booster Type')
ggsave(paste0(dir$plots,'FoldChangeFromVaccines_PriorDoses.pdf'),fc_plot_4, width=6,height=5)
