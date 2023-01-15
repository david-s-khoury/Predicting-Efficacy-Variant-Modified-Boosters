vsv_improvement_plot_split_by_valency_omicron = ggplot(filter(vsv_improvement_data_new,
                                                      Variant!='Ancestral', 
                                                      (VSVBoosterPart1=='Omicron BA.1' & VSVBoosterPart2=='') |(VSVBoosterPart1=='Ancestral' & VSVBoosterPart2=='Omicron BA.1')
                                                      ), 
                                               aes(x=Valency,y=VSVimprovement,fill=Valency))+
  geom_boxplot()+
  geom_point(aes(x=Valency,colour=Variant, shape = FirstAuthor),size=3,position = position_jitter(width=.1))+
  stat_summary(aes(label=round(after_stat(y), 2), fontface='bold'),fun=geomean, geom="text",vjust=-1,position = position_dodge(0.9)) +
  stat_compare_means(comparisons = list(c('Bivalent','Monovalent')),paired=F,method='t.test',label=c("p.signif"))+
  #stat_compare_means(paired=F,method='t.test',label.y = 2.75)+
  theme_classic()+
  scale_shape_manual(name = 'Reference',values=study_shapes, labels = study_shape_labels_2)+
  scale_fill_manual(name = 'Vaccine composition', guide='none',values=valency_colours)+
  scale_colour_manual(name='Variant Tested',values = variant_colours)+
  facet_wrap(~Matching)
  labs(y=y_improv_label,x=paste(x_valency_improv_label, '(must include BA.1)'))
print(vsv_improvement_plot_split_by_valency_omicron)
ggsave(paste0(dir$plots,'ImprovementFromVSV_ValencySplit_Omicron.pdf'),vsv_improvement_plot_split_by_valency, width=6,height=5)
