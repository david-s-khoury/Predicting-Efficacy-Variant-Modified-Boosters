# Combine plots for figure in paper
no_title_no_legend = theme(plot.title = element_blank(),legend.position="none")
no_y_axis = theme(axis.title.y = element_blank(), axis.line.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

# Data focused plot
paper_data_plot_A = fc_plot_use
paper_data_plot_B = vsv_improvement_plot_all_data_by_variant
paper_data_plot_C = vsv_improvement_plot_split_by_matching
paper_data_plot=plot_grid(paper_data_plot_A+no_title_no_legend,
                          paper_data_plot_B+no_title_no_legend, 
                          paper_data_plot_C+no_title_no_legend+no_y_axis,
                          nrow=1,
                          rel_widths = c(1,1.5,1.5),
                          labels = list('A','B','C'))
# Legend for data focused plot
paper_legend_plot = plot_grid(#as_ggplot(get_legend(paper_data_plot_A)),
                              as_ggplot(get_legend(paper_data_plot_B)),
                              #as_ggplot(get_legend(sub_plot1a)),
                              ncol=1)
                              #align = 'v', rel_heights=c(1,1.5))
paper_data_plot_final = plot_grid(paper_data_plot, paper_legend_plot,nrow=1,rel_widths = c(4,1))
#print(paper_data_plot_final)
ggsave(paste0(dir$plots,'Paper_Data_Plot.pdf'), paper_data_plot_final, width=12,height=6)
ggsave(paste0(dir$manuscript_plots,'Figure_1_Data_Plot.pdf'), paper_data_plot_final, width=12,height=6)



paper_sim_plot_AB = sub_plot1_avg_over_6m
paper_sim_plot_CD = sub_plot_improvement_avg_over_6m_by_outcome
paper_sim_plot=plot_grid(paper_sim_plot_AB+no_title_no_legend,
                         paper_sim_plot_CD+no_title_no_legend, 
                         nrow=2,labels = list('A','C'))
paper_sim_legend_plot = plot_grid(as_ggplot(get_legend(paper_sim_plot_AB)),
                                  as_ggplot(get_legend(paper_sim_plot_CD)),
                                  ncol=1, 
                                  align = 'v')
paper_sim_plot_final = plot_grid(paper_sim_plot, paper_sim_legend_plot,nrow=1,rel_widths = c(4,1.2))
ggsave(paste0(dir$plots,'Paper_Sim_Plot.pdf'), paper_sim_plot_final, width=8,height=8)
ggsave(paste0(dir$manuscript_plots,'Figure_2_Sim_Plot.pdf'), paper_sim_plot_final, width=11,height=8)

supp_plot_1 = plot_grid(vsv_improvement_plot_split_by_valency+no_title_no_legend,
                        vsv_improvement_plot_split_by_priorstatus+no_title_no_legend+theme(axis.title.y = element_blank()),
                        vsv_improvement_plot_split_by_priordoses+no_title_no_legend+theme(axis.title.y = element_blank()),
                        as_ggplot(get_legend(vsv_improvement_plot_split_by_valency)),
                        nrow=1, 
                        align = 'h',
                        rel_widths = c(1.1,1,1,.5),
                        labels = list('A','B','C'))
ggsave(paste0(dir$manuscript_plots,'Figure_S1_ValencyStatusDoses_Plot.pdf'), supp_plot_1, width=14,height=5)


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