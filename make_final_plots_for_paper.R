# Combine plots for figure in paper
no_title_no_legend = theme(plot.title = element_blank(),legend.position="none")
no_y_axis = theme(axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
facet_adjust = theme(strip.background.x = element_rect(colour=NULL,linetype = 'blank'),
                           strip.text = element_text(face='bold'))
large_font_theme = theme(text = element_text(size=16))
# Data focused plot
paper_data_plot_A = fc_plot_use
paper_data_plot_B = vsv_improvement_plot_all_data_by_variant
paper_data_plot_C = vsv_improvement_plot_split_by_matching
paper_data_plot=plot_grid(paper_data_plot_A+no_title_no_legend+large_font_theme + theme(plot.margin = unit(c(4,4,1,1), "lines")),
                          paper_data_plot_B+no_title_no_legend+large_font_theme+ theme(plot.margin = unit(c(4,3,1,0), "lines")), 
                          paper_data_plot_C+no_title_no_legend+no_y_axis+large_font_theme+ theme(plot.margin = unit(c(4,0,1,0), "lines")),
                          nrow=1,
                          rel_widths = c(2.3,2.2,1.5),
                          labels = list('A','B','C'))

# Legend for data focused plot
dummy_plot = paper_data_plot_B
other_variant_list =  c('Omicron BA.2.75','Omicron BA.2.75.2','Omicron BA.2','Omicron BA.4.6','OmicronBF.7','OmicronBQ.1.1',  
                        'OmicronXBB.1','BA.2.75.2','BQ.1.1')
dummy_plot$data = dummy_plot$data %>% 
  mutate(Variant=ifelse(Variant %in% other_variant_list, 
                        'Omicron (other)',Variant))
dummy_variant_colours = variant_colours[names(variant_colours) %in% unique(dummy_plot$data$Variant)]
dummy_plot = dummy_plot+scale_colour_manual(values=dummy_variant_colours)
paper_legend_plot = plot_grid(#as_ggplot(get_legend(paper_data_plot_A)),
  as_ggplot(get_legend(dummy_plot)$grobs[[1]])+large_font_theme,
  as_ggplot(get_legend(dummy_plot)$grobs[[3]])+large_font_theme,
                              #as_ggplot(get_legend(sub_plot1a)),
                              ncol=1)
                              #align = 'v', rel_heights=c(1,1.5))
paper_legend_plot_2col = plot_grid(as_ggplot(get_legend(dummy_plot)$grobs[[2]])+large_font_theme,
                                   as_ggplot(get_legend(dummy_plot)$grobs[[3]])+large_font_theme,
                                   ncol=2, align = 'h')
paper_data_plot_final = plot_grid(paper_data_plot, paper_legend_plot_2col,nrow=1,rel_widths = c(4.2,1))

#print(paper_data_plot_final)
ggsave(paste0(dir$plots,'Paper_Data_Plot.pdf'), paper_data_plot_final, width=12,height=6)
ggsave(paste0(dir$manuscript_plots,'Figure_1_Data_Plot.pdf'), paper_data_plot_final, width=14,height=6)



paper_sim_plot_AB_legend = sub_plot1_avg_over_6m
paper_sim_plot_AB = sub_plot1_avg_over_6m_CI+ theme(panel.spacing = unit(2, "lines"))
paper_sim_plot_CD_legend = sub_plot_improvement_avg_over_6m_by_outcome
paper_sim_plot_CD = sub_plot_improvement_avg_over_6m_by_outcome_CI+ theme(panel.spacing = unit(2, "lines"))
paper_sim_plot=plot_grid(paper_sim_plot_AB+no_title_no_legend+large_font_theme+facet_adjust+
                           theme(plot.margin = unit(c(0,0,1,1), "lines")),
                         paper_sim_plot_CD+no_title_no_legend+large_font_theme+facet_adjust +
                           theme(plot.margin = unit(c(1,0,0,1), "lines")), 
                         nrow=2,labels = list('A','C'))
paper_sim_legend_plot = plot_grid(as_ggplot(get_legend(paper_sim_plot_AB_legend)),
                                  as_ggplot(get_legend(paper_sim_plot_CD_legend)),
                                  ncol=1, 
                                  align = 'v')
paper_sim_plot_final = plot_grid(paper_sim_plot, 
                                 paper_sim_legend_plot,
                                 nrow=1,rel_widths = c(4,1.2))
ggsave(paste0(dir$plots,'Paper_Sim_Plot.pdf'), paper_sim_plot_final, width=8,height=8)
ggsave(paste0(dir$manuscript_plots,'Figure_2_Sim_Plot.pdf'), paper_sim_plot_final, width=11,height=8)


# supp_plot_1 = plot_grid(fc_plot_supp2+no_title_no_legend+facet_adjust,
#                         paper_legend_plot,nrow=1,rel_widths = c(4,1.5))
# 
# supp_plot_1alldata = plot_grid(fc_plot_supp2_alldata+no_title_no_legend+facet_adjust,
#                         paper_legend_plot,nrow=1,rel_widths = c(4,1.5))
# supp_plot_1paireddata = plot_grid(fc_plot_supp2_paireddata+no_title_no_legend+facet_adjust,
#                                paper_legend_plot,nrow=1,rel_widths = c(4,1.5))
# supp_plot_1alldata_withlines = plot_grid(fc_plot_supp2_alldata_withlines+no_title_no_legend+facet_adjust,
#                                   paper_legend_plot,nrow=1,rel_widths = c(4,1.5))
# 
# ggsave(paste0(dir$manuscript_plots,'Figure_S1_AllBoosters_Plot.pdf'), supp_plot_1, width=9,height=6)
# ggsave(paste0(dir$manuscript_plots,'Figure_S1_alldata_AllBoosters_Plot.pdf'), supp_plot_1alldata, width=9,height=6)
# ggsave(paste0(dir$manuscript_plots,'Figure_S1_paireddata_AllBoosters_Plot.pdf'), supp_plot_1paireddata, width=9,height=6)
# ggsave(paste0(dir$manuscript_plots,'Figure_S1_alldata_withlines_AllBoosters_Plot.pdf'), supp_plot_1alldata_withlines, width=9,height=6)



#margin
#c(t,r,b,l)
supp_plot_BoosterType = plot_grid(vsv_improvement_plot_all_data_by_vaccine+no_title_no_legend+facet_adjust+large_font_theme+
                                    stat_summary(aes(label = str_c('n=',round(..y..,2)), y = stage(log10(VSVimprovement), after_stat = log10(.5))),fun=length, geom="text") +
                          theme(plot.margin = unit(c(0,2,0,4), "lines")),
                        paper_legend_plot_2col,
                        nrow=1,rel_widths = c(4,2.5))
ggsave(paste0(dir$manuscript_plots,'Figure_S1_BoosterType_Plot.pdf'), supp_plot_BoosterType, width=10,height=8)


supp_plot_ValencyStatusDoses = plot_grid(vsv_improvement_plot_split_by_valency+no_title_no_legend+large_font_theme+
                                           stat_summary(aes(label = str_c('n=',round(..y..,2)), y = stage(log10(VSVimprovement), after_stat = log10(.5))),fun=length, geom="text") ,
                        vsv_improvement_plot_split_by_priorstatus+no_title_no_legend+theme(axis.title.y = element_blank())+large_font_theme+
                          stat_summary(aes(label = str_c('n=',round(..y..,2)), y = stage(log10(VSVimprovement), after_stat = log10(.5))),fun=length, geom="text") ,
                        vsv_improvement_plot_split_by_priordoses+no_title_no_legend+theme(axis.title.y = element_blank())+large_font_theme+
                          stat_summary(aes(label = str_c('n=',round(..y..,2)), y = stage(log10(VSVimprovement), after_stat = log10(.5))),fun=length, geom="text"),
                        paper_legend_plot_2col,
                        nrow=1, 
                        align = 'h',
                        rel_widths = c(1.1,1,1,1),
                        labels = list('A','B','C'))
ggsave(paste0(dir$manuscript_plots,'Figure_S2_ValencyStatusDoses_Plot.pdf'), supp_plot_ValencyStatusDoses, width=14,height=5)



supp_plot_ValencyMatching = plot_grid(vsv_improvement_plot_pure_or_mix2+no_title_no_legend+facet_adjust+large_font_theme+
                                        stat_summary(aes(label = str_c('n=',round(..y..,2)), y = stage(log10(VSVimprovement), after_stat = log10(.5))),fun=length, geom="text"),
                        paper_legend_plot_2col,nrow=1,rel_widths = c(4,2.5))
ggsave(paste0(dir$manuscript_plots,'Figure_S3_ValencyMatching_Plot.pdf'), supp_plot_ValencyMatching, width=10,height=4)




make_older_plots=F
if (make_older_plots){
  
  # Older plots
  # Simulation plots
  paper_sim_plot_D = sub_plot1_sympt_avg_over_6m
  paper_sim_plot_E = sub_plot_improvement_avg_over_6m
  paper_sim_plot=plot_grid(paper_sim_plot_D+no_title_no_legend,
                           paper_sim_plot_E+no_title_no_legend, 
                           nrow=1,labels = list('D','E'))
  paper_sim_legend_plot = plot_grid(as_ggplot(get_legend(paper_sim_plot_D)),
                                    as_ggplot(get_legend(paper_sim_plot_E)),
                                    ncol=1, 
                                    align = 'v')
  paper_sim_plot_final = plot_grid(paper_sim_plot, paper_sim_legend_plot,nrow=1,rel_widths = c(4,1))
  ggsave(paste0(dir$plots,'Paper_Sim_Plot_Old.pdf'), paper_sim_plot_final, width=12,height=6)
  
  # Final plot for paper
  paper_plot = plot_grid(paper_data_plot,paper_sim_plot, rel_widths=c(1,1),nrow=2)
  paper_legend = plot_grid(as_ggplot(get_legend(paper_data_plot_A)),
                           as_ggplot(get_legend(paper_data_plot_B)),
                           as_ggplot(get_legend(paper_sim_plot_D)),
                           as_ggplot(get_legend(paper_sim_plot_E)),
                           ncol=1, 
                           align = 'v',
                           rel_heights=c(1.5,2,.7,.7))
  print(paper_legend)
  paper_plot_total = plot_grid(paper_plot,
                               paper_legend,
                               ncol=2,
                               rel_widths=c(4,1))
  ggsave(paste0(dir$plots,'Paper_Plot_Old.pdf'), paper_plot_total, width=12,height=9)
  ggsave(paste0(dir$plots,'Paper_Plot_Old.jpg'), paper_plot_total, width=12,height=9)
  
}