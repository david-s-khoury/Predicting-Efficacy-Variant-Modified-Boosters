# Now plot the fold changes in neut titre after boosting
# We will do this in a number of ways, including by vaccine type

# ploting variables / constants
fc_breaks=c(1,2,5,10,20,50,100)
y_label_anc_fold_change='Fold change in Neutralisaiton GMT after boosting'
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
print(fc_plot_supp)
ggsave(paste0(dir$plots,'FoldChangeFromVaccines_AllVaccines.pdf'),fc_plot_supp, width=5,height=5)



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
