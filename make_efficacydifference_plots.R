# For plotting
percent_breaks = seq(0,1,.2)
percent_break_labels = paste0(100*percent_breaks,'%')
neut_breaks = 10^seq(-5,5,1)
neut_break_labels = as.character(neut_breaks) 
improv_labels=paste0(rises,'-fold')
improv_colours = colorRampPalette(brewer.pal(9, "RdPu"))(8)[5:8]
pre_boost_label_text='pre-boost population protection from symptomatic infection'
base_improv_label_text='Percentage point improvement over no boost'
base_improv_title_text = 'Improvement over no boost'
y_label_efficacy = 'Protection'
percentage_point_y_title = 'Percentage point improvement vs ancestral booster'
improv_legend_name = 'Additional protection\n(variant boost vs ancestral)'

theme_for_sim_plots = theme_classic()+
                        theme(legend.key.width = unit(50,'pt'),
                              legend.key.height = unit(30,'pt'))
fold_improvement_legend_name = 'Predicted Protection'
x_expand = c(-.02,1.05)
y_improvement_expand =c(0,15)

pop_multiplier=1

tall_width = 8
tall_height = 12
# Which gmrs will we plot?
gmrs_for_subplots=c(1,ancestral_boost, ancestral_boost*vsv_rises) # add in ancestral_boost*matched / unmathed rises

gmr_use_indices = c(1,2,3,4,5)
vsv_rise_use_indices = gmr_use_indices-2
vsv_rise_use_indices = vsv_rise_use_indices[vsv_rise_use_indices>0]

gmr_colours = c('1'='grey','ancestral'='#088AFF','vsv'='#FF5050', 'vsv_matched'='#AF0404', 'vsv_nonmatched'='#FF9A8F')
gmr_colours = gmr_colours[c(1:(2+length(vsv_rises)))] # use the right number of colours
names(gmr_colours) = c('1',ancestral_boost,ancestral_boost*vsv_rises)
gmr_colours = gmr_colours[gmr_use_indices]

gmr_fills = c('1'='#FFFCC5','ancestral'=NA,'vsv'=NA, 'vsv_matched'=NA, 'vsv_nonmatched'=NA)
gmr_fills = gmr_fills[c(1:(2+length(vsv_rises)))] # use the right number of fills
names(gmr_fills) = c('1',ancestral_boost,ancestral_boost*vsv_rises)
gmr_fills = gmr_fills[gmr_use_indices]
alpha_fill = .3

gmr_linetypes = c('1'='solid','ancestral'='solid','vsv'='solid', 'vsv_matched'='32', 'vsv_nonmatched'='32')
gmr_linetypes = gmr_linetypes[c(1:(2+length(vsv_rises)))] # use the right number of fills
names(gmr_linetypes) = c('1',ancestral_boost,ancestral_boost*vsv_rises)
gmr_linetypes = gmr_linetypes[gmr_use_indices]

gmr_labels = c('pre-boost protection',
               paste0('\n(',ancestral_boost,'-fold)'))
for (vsv_rise in vsv_rises){
  gmr_labels = c(gmr_labels,
                 paste0('\n(additional ',vsv_rise,'-fold)') )
  
}
gmr_label_starts = c("","Ancestral Boost","Any Variant Boost","Matched Variant","Non-Matched Variant")
gmr_labels = paste0(gmr_label_starts, gmr_labels)
gmr_labels = gmr_labels[gmr_use_indices]

gmr_sizes = c(1,1,1,.5,.5)
names(gmr_sizes) = names(gmr_colours)
gmr_sizes = gmr_sizes[gmr_use_indices]


vsv_linetypes = tail(gmr_linetypes,3)
vsv_linetypes = vsv_linetypes[c(1:(length(vsv_rises)))] # use the right number of fills
names(vsv_linetypes) = c(vsv_rises)
vsv_linetypes = vsv_linetypes[pmax(0,gmr_use_indices-2)]

vsv_labels = c('Any Variant Boost\n','Matched Variant Boost\n','Non-Matched Variant Boost\n')
for (i in c(1:length(vsv_rises))){
  vsv_rise = vsv_rises[i]
  vsv_labels[i] = paste0(vsv_labels[i],'(additional ',vsv_rise,'-fold)')
}
vsv_labels = vsv_labels[c(1:(length(vsv_rises)))] # use the right number of fills
vsv_labels = vsv_labels[pmax(0,gmr_use_indices-2)]

vsv_labels = tail(gmr_labels,3)
names(vsv_labels) = c(vsv_rises)

vsv_sizes = tail(gmr_sizes,3)
names(vsv_sizes) = c(vsv_rises)
vsv_sizes = vsv_sizes[pmax(0,gmr_use_indices-2)]

improvement_cols = c('sympt_vsv'='#009933','sympt_vsv_matched'='#006622','sympt_vsv_nonmatched'='#00e64d','severe_vsv'='#ff9900','severe_vsv_matched'='#b36b00','severe_vsv_nonmatched'='#ffc266')
names(improvement_cols) = c(paste0(outcomes[1],vsv_rises),paste0(outcomes[2],vsv_rises))
improvement_cols = tail(improvement_cols,3)
names(improvement_cols) = vsv_rises

improvement_linetypes = c(vsv_linetypes,vsv_linetypes)
improvement_linetypes = vsv_linetypes
names(improvement_linetypes) = names(improvement_cols)

improvement_labels = c(paste(outcomes[1],vsv_labels,sep='\n'),paste(outcomes[2],vsv_labels,sep='\n'))
improvement_labels = vsv_labels

plot_neut_eff_struct_gmr_use = plot_neut_eff_struct %>%
  filter(GMR %in% gmrs_for_subplots[gmr_use_indices])
plot_neut_eff_struct_gmr_use_improvement = plot_neut_eff_struct_improvement %>%
  filter(GMRdiff %in% vsv_rises[vsv_rise_use_indices])


# This is for the vertical / horizontal dashed lines at 50% pre-boost eff
midpointSymptData = filter(plot_neut_eff_struct_gmr_use, 
                           GMR==ancestral_boost, 
                           abs(baseline_sympt_eff-.5)==min(abs(baseline_sympt_eff-.5)),
                           outcome==outcomes[1])
midpointSevereData = filter(plot_neut_eff_struct_gmr_use, 
                           GMR==ancestral_boost, 
                           abs(baseline_sympt_eff-.5)==min(abs(baseline_sympt_eff-.5)),
                           outcome==outcomes[2])

# New efficacy for symptomatic only
sub_plot1_sympt = ggplot(filter(plot_neut_eff_struct_gmr_use,outcome==outcomes[1]), aes(x=baseline_sympt_eff,y=pop_multiplier*new_eff, colour=as.factor(GMR), linetype = as.factor(GMR)))+
  geom_ribbon(aes(ymax=pop_multiplier*new_eff,fill=as.factor(GMR)),ymin=0,alpha=alpha_fill,colour=NA)+
  geom_line(mapping = aes(size=as.factor(GMR))) +
  scale_size_manual(values = gmr_sizes, guide='none')+
  theme_for_sim_plots+
  scale_colour_manual(values=gmr_colours,name=fold_improvement_legend_name, labels = gmr_labels)+
  scale_fill_manual(values=gmr_fills,name=fold_improvement_legend_name, labels = gmr_labels,guide='none', na.value = NA)+
  scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
  scale_x_continuous(breaks = percent_breaks,labels = percent_break_labels)+
  scale_y_continuous(breaks = percent_breaks, labels = percent_break_labels, expand=c(0,0))+
  labs(title=y_label_efficacy, x=pre_boost_label_text,y=y_label_efficacy)+
  geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=pop_multiplier),colour='black',size=.5,linetype='dashed')+
  geom_segment(data = midpointSymptData, aes(x=0,y=pop_multiplier*new_eff, xend=baseline_sympt_eff, yend=pop_multiplier*new_eff),colour='black',size=.5,linetype='dashed')+
  coord_cartesian(xlim = x_expand, expand=F)
#print(sub_plot1_sympt)

# New efficacy for both outcomes
sub_plot1 = sub_plot1_sympt+
  facet_wrap(~outcome, scales = 'free_y')
sub_plot1$data=plot_neut_eff_struct_gmr_use

# This is protection after 6 months
if(do_after_sixm){
  sub_plot1_sympt_after_6m = ggplot(filter(plot_neut_eff_struct_gmr_use,outcome==outcomes[1]), aes(x=baseline_sympt_eff,y=pop_multiplier*new_eff_after_6m, colour=as.factor(GMR), linetype = as.factor(GMR)))+
    geom_ribbon(aes(ymax=pop_multiplier*new_eff_after_6m,fill=as.factor(GMR)),ymin=0,alpha=alpha_fill,colour=NA)+
    geom_line(mapping = aes(size=as.factor(GMR))) +
    scale_size_manual(values = gmr_sizes, guide='none')+
    theme_for_sim_plots+
    scale_colour_manual(values=gmr_colours,name=fold_improvement_legend_name, labels = gmr_labels)+
    scale_fill_manual(values=gmr_fills,name=fold_improvement_legend_name, labels = gmr_labels,guide='none', na.value = NA)+
    scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
    scale_x_continuous(breaks = percent_breaks,labels = percent_break_labels,expand=c(0.001,0.001))+
    scale_y_continuous(breaks = percent_breaks, labels = percent_break_labels, expand=c(0,0))+
    labs(title=y_label_efficacy, x=pre_boost_label_text,y=y_label_efficacy)+
    geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=pop_multiplier),colour='black',size=.5,linetype='dashed')+
    geom_segment(data = midpointSevereData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=pop_multiplier),colour='black',size=.5,linetype='dashed')+
    coord_cartesian(xlim = x_expand, expand=F)
  #print(sub_plot1_sympt_after_6m)
  
  # New efficacy for both outcomes
  sub_plot1_after_6m = sub_plot1_sympt_after_6m+
    facet_wrap(~outcome, scales = 'free_y')
  sub_plot1_after_6m$data=plot_neut_eff_struct_gmr_use
  
  ggsave(paste0(dir$plots,'Now_Later_Saved.pdf'),plot_grid(sub_plot1,sub_plot1_after_6m,nrow=2), width=8,height=8)
}


# This is average protection over next 6 months
if (include_average_over_6m){
  sub_plot1_sympt_avg_over_6m = ggplot(filter(plot_neut_eff_struct_gmr_use,outcome==outcomes[1]), aes(x=baseline_sympt_eff,y=pop_multiplier*avg_eff_over_6m, colour=as.factor(GMR), linetype = as.factor(GMR)))+
    geom_ribbon(aes(ymax=pop_multiplier*avg_eff_over_6m,fill=as.factor(GMR)),ymin=0,colour=NA, alpha=alpha_fill)+
    geom_line(mapping = aes(size=as.factor(GMR))) +
    scale_size_manual(values = gmr_sizes, guide='none')+
    theme_for_sim_plots+
    scale_colour_manual(values=gmr_colours,name=fold_improvement_legend_name, labels = gmr_labels)+
    scale_fill_manual(values=gmr_fills,name=fold_improvement_legend_name, labels = gmr_labels,guide='none', na.value = NA)+
    scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
    scale_x_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'))+
    scale_y_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'))+
    labs(title=y_label_efficacy, x=pre_boost_label_text,y=y_label_efficacy)+
    geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=pop_multiplier),colour='black',size=.5,linetype='dashed')+
    geom_segment(data = midpointSevereData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=pop_multiplier),colour='black',size=.5,linetype='dashed')+
    coord_cartesian(xlim = x_expand, expand=F)
  #print(sub_plot1_sympt_avg_over_6m)
  # New efficacy for both outcomes
  sub_plot1_avg_over_6m = sub_plot1_sympt_avg_over_6m+
    facet_wrap(~outcome, scales = 'free_y')
  sub_plot1_avg_over_6m$data=plot_neut_eff_struct_gmr_use

  ggsave(paste0(dir$plots,'Now_Later_Average_Saved.pdf'),plot_grid(sub_plot1,sub_plot1_after_6m, sub_plot1_avg_over_6m,nrow=3), width=tall_width,height=tall_height)
}


# Now we look at difference between boosting and not boosting
line_max = max(with(plot_neut_eff_struct_gmr_use,pop_multiplier*(new_eff-baseline_eff)),na.rm=T)
sub_plot2 = ggplot(plot_neut_eff_struct_gmr_use, aes(x=baseline_sympt_eff,y=pop_multiplier*(new_eff-baseline_eff), colour=as.factor(GMR), linetype = as.factor(GMR)))+
  geom_line(mapping = aes(size=as.factor(GMR))) +
  scale_size_manual(values = gmr_sizes, guide='none')+
  theme_for_sim_plots+
  scale_colour_manual(values=gmr_colours,name=fold_improvement_legend_name, labels = gmr_labels)+
  scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
  scale_x_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'))+
  scale_y_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'),expand=c(0.1,0))+
  labs(title=base_improv_title_text, x=pre_boost_label_text,y=base_improv_label_text)+
  facet_wrap(~outcome, scales = 'free_y')+
  geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=line_max),colour='black',size=.5,linetype='dashed')+
  geom_segment(data = midpointSymptData, aes(x=0,y=pop_multiplier*(new_eff-baseline_eff), xend=baseline_sympt_eff, yend=pop_multiplier*(new_eff-baseline_eff)),colour='black',size=.5,linetype='dashed')+
  coord_cartesian(xlim = x_expand, expand=F)
#print(sub_plot2)

if(do_after_sixm){
  line_max = max(with(plot_neut_eff_struct_gmr_use,pop_multiplier*(new_eff_after_6m-baseline_eff_after_sixm)),na.rm=T)
  sub_plot2_after_6m = ggplot(plot_neut_eff_struct_gmr_use, aes(x=baseline_sympt_eff,y=pop_multiplier*(new_eff_after_6m-baseline_eff_after_sixm), colour=as.factor(GMR), linetype = as.factor(GMR)))+
    geom_line(mapping = aes(size=as.factor(GMR))) +
    scale_size_manual(values = gmr_sizes, guide='none')+
    theme_for_sim_plots+
    scale_colour_manual(values=gmr_colours,name=fold_improvement_legend_name, labels = gmr_labels)+
    scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
    scale_x_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'),expand=c(0.1,0))+
    scale_y_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'),expand=c(0.1,0))+
    labs(title=paste0(base_improv_title_text, ' (After 6m)'), x=pre_boost_label_text,y=base_improv_label_text)+
    facet_wrap(~outcome, scales = 'free_y')+
    geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=line_max),colour='black',size=.5,linetype='dashed')+
    geom_segment(data = midpointSymptData, aes(x=0,y=pop_multiplier*(new_eff_after_6m-baseline_eff_after_sixm), xend=baseline_sympt_eff, yend=pop_multiplier*(new_eff_after_6m-baseline_eff_after_sixm)),colour='black',size=.5,linetype='dashed')+
    coord_cartesian(xlim = x_expand, expand=F)
  #print(sub_plot2_after_6m)
  ggsave(paste0(dir$plots,'Now_Later_Difference_Saved.pdf'),plot_grid(sub_plot2,sub_plot2_after_6m,nrow=2), width=8,height=8)
}
if (include_average_over_6m){
  line_max = max(with(plot_neut_eff_struct_gmr_use,pop_multiplier*avg_eff_over_6m_cf_baseline),na.rm=T)
  sub_plot2_avg_over_6m = ggplot(plot_neut_eff_struct_gmr_use, aes(x=baseline_sympt_eff,y=pop_multiplier*avg_eff_over_6m_cf_baseline, colour=as.factor(GMR), linetype = as.factor(GMR)))+
    geom_line(mapping = aes(size=as.factor(GMR))) +
    scale_size_manual(values = gmr_sizes, guide='none')+
    theme_for_sim_plots+
    scale_colour_manual(values=gmr_colours,name=fold_improvement_legend_name, labels = gmr_labels)+
    scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
    scale_x_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'),expand=c(0,0))+
    scale_y_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'),expand=c(0.1,0))+
    labs(title=paste0(base_improv_title_text, ' (Avg over 6m)'), x=pre_boost_label_text, y=base_improv_label_text)+
    facet_wrap(~outcome, scales = 'free_y')+
    geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=line_max),colour='black',size=.5,linetype='dashed')+
    geom_segment(data = midpointSymptData, aes(x=0,y=pop_multiplier*avg_eff_over_6m_cf_baseline, xend=baseline_sympt_eff, yend=pop_multiplier*avg_eff_over_6m_cf_baseline),colour='black',size=.5,linetype='dashed')+
    coord_cartesian(xlim = x_expand, expand=F)
  #print(sub_plot2b)
  ggsave(paste0(dir$plots,'Now_Later_Average_Difference_Saved.pdf'),plot_grid(sub_plot2,sub_plot2_after_6m, sub_plot2_avg_over_6m,nrow=3), width=tall_width,height=tall_height)
}

# Save instantaneous eff and diff after boosting plots as PaperPlot_1
ggsave(paste0(dir$plots,'PaperPlot_1.pdf'),plot_grid(sub_plot1,sub_plot2,nrow=2), width=8,height=8)

midpointSymptData = filter(plot_neut_eff_struct_gmr_use_improvement, 
                           GMRdiff==vsv_rise, 
                           from==ancestral_boost,
                           abs(baseline_sympt_eff-.5)==min(abs(baseline_sympt_eff-.5)))
line_max = max(with(plot_neut_eff_struct_gmr_use_improvement,pop_multiplier*value),na.rm=T)

sub_plot3 = ggplot(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise), aes(x=baseline_sympt_eff,y=pop_multiplier*value, colour=as.character(GMRdiff)))+
  geom_line(mapping = aes(size=as.character(GMRdiff))) +
  scale_size_manual(values = vsv_sizes)+
  theme_for_sim_plots+
  scale_colour_manual(values=improv_colours,name='Improvement of\nNew Vaccine', labels = improv_labels, guide='none')+
  scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
  scale_x_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'))+
  scale_y_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'),expand=c(0.1,0))+
  labs(title='Additional Cases Averted per 1,000 Cases in an Unvaccinated Population', x='pre-Boost Population Symptomatic Efficacy',y='Additional Cases Averted / 1,000 Original Cases')+
  #facet_grid(from~outcome, scales = 'free_y',switch = 'y',labeller = as_labeller(c(`1` = "Vs No Boost", `2` = "Vs Two Fold Boost", '4'='Vs Four Fold Boost','6'='Vs Six Fold Boost', 'symptomatic'='Symptomatic Cases','severe'='Severe Cases')))+
  facet_wrap(~outcome)+
  geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=line_max),colour='black',size=.5,linetype='dashed')+
  geom_segment(data = midpointSymptData, aes(x=0,y=pop_multiplier*value, xend=baseline_sympt_eff, yend=pop_multiplier*value),colour='black',size=.5,linetype='dashed')+
  theme(strip.placement = 'outside')+
  geom_text(data = midpointSymptData,aes(x=.25,y=pop_multiplier*value,label=paste0(round(pop_multiplier*value,0),' cases/1,000')))+
  coord_cartesian(xlim = x_expand, expand=F)

#print(sub_plot3)
ggsave(paste0(dir$plots,'PaperPlot_2.pdf'),sub_plot3, width=8,height=5)

if(do_after_sixm){
  sub_plot3_after_6m = ggplot(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise), aes(x=baseline_sympt_eff,y=pop_multiplier*value6m, colour=as.character(GMRdiff)))+
    geom_line(mapping = aes(size=as.character(GMRdiff))) +
    scale_size_manual(values = vsv_sizes)+
    theme_for_sim_plots+
    scale_colour_manual(values=improv_colours,name='Improvement of\nNew Vaccine', labels = improv_labels, guide='none')+
    scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
    scale_x_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'))+
    scale_y_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'),expand=c(0.1,0))+
    labs(title='Additional Cases Averted after 6 months per 1,000 Cases in an Unvaccinated Population', x='pre-Boost Population Symptomatic Efficacy',y='Additional Cases Averted / 1,000 Original Cases')+
    #facet_grid(from~outcome, scales = 'free_y',switch = 'y',labeller = as_labeller(c(`1` = "Vs No Boost", `2` = "Vs Two Fold Boost", '4'='Vs Four Fold Boost','6'='Vs Six Fold Boost', 'symptomatic'='Symptomatic Cases','severe'='Severe Cases')))+
    facet_wrap(~outcome)+
    geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=line_max),colour='black',size=.5,linetype='dashed')+
    geom_segment(data = midpointSymptData, aes(x=0,y=pop_multiplier*value6m, xend=baseline_sympt_eff, yend=pop_multiplier*value6m),colour='black',size=.5,linetype='dashed')+
    theme(strip.placement = 'outside')+
    geom_text(data = midpointSymptData,aes(x=.25,y=pop_multiplier*value6m,label=paste0(round(pop_multiplier*value6m,0),' cases/1,000')))+
    coord_cartesian(xlim = x_expand, expand=F)
  #print(sub_plot3a)
}
if(include_average_over_6m){
  sub_plot3_avg_over_6m = ggplot(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rise), aes(x=baseline_sympt_eff,y=pop_multiplier*value6mavg, colour=as.character(GMRdiff)))+
    geom_line(mapping = aes(size=as.character(GMRdiff))) +
    scale_size_manual(values = vsv_sizes)+
    theme_for_sim_plots+
    scale_colour_manual(values=improv_colours,name='Improvement of\nNew Vaccine', labels = improv_labels, guide='none')+
    scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
    scale_x_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'))+
    scale_y_continuous(breaks = percent_breaks,labels = paste0(100*percent_breaks,'%'),expand=c(0.1,0))+
    labs(title='Additional Cases Averted after 6 months per 1,000 Cases in an Unvaccinated Population', x='pre-Boost Population Symptomatic Efficacy',y='Additional Cases Averted / 1,000 Original Cases')+
    #facet_grid(from~outcome, scales = 'free_y',switch = 'y',labeller = as_labeller(c(`1` = "Vs No Boost", `2` = "Vs Two Fold Boost", '4'='Vs Four Fold Boost','6'='Vs Six Fold Boost', 'symptomatic'='Symptomatic Cases','severe'='Severe Cases')))+
    facet_wrap(~outcome)+
    geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=line_max),colour='black',size=.5,linetype='dashed')+
    geom_segment(data = midpointSymptData, aes(x=0,y=pop_multiplier*value6mavg, xend=baseline_sympt_eff, yend=pop_multiplier*value6mavg),colour='black',size=.5,linetype='dashed')+
    theme(strip.placement = 'outside')+
    geom_text(data = midpointSymptData,aes(x=.25,y=pop_multiplier*value6mavg,label=paste0(round(pop_multiplier*value6mavg,0),' cases/1,000')))+
    coord_cartesian(xlim = x_expand, expand=F)
  #print(sub_plot3b)
  ggsave(paste0(dir$plots,'OverSixm_Average_Difference_Saved.pdf'),plot_grid(sub_plot3,sub_plot3_after_6m, sub_plot3_avg_over_6m,nrow=3), width=tall_width,height=tall_height)
}
#Save the difference as Paper_Plot_2
ggsave(paste0(dir$plots,'PaperPlot_2.pdf'),sub_plot3, width=8,height=5)

scale_color_manual(values = colorRampPalette(brewer.pal(9, "YlGnBu"))(12)[6:12])

#print 2x2 plot
# grid_plot_for_paper = plot_grid(fc_plot_1+theme(plot.title = element_blank(),legend.position="none"),
#                                 vsv_improvement_plot_all_data_by_variant+labs(y='Fold improvement for Variant Vaccine')+theme(plot.title = element_blank(),legend.position="none"),
#                                 vsv_improvement_plot_HHsplit+labs(y='Fold improvement for Variant Vaccine')+theme(plot.title = element_blank(),legend.position="none"),
#                                 sub_plot1_sympt+theme(plot.title = element_blank(), legend.position="none"),
#                                 ncol=2,nrow=2)
# legends_for_paper = plot_grid(as_ggplot(get_legend(fc_plot_2)),
#                               as_ggplot(get_legend(vsv_improvement_plot_HHsplit)),
#                               as_ggplot(get_legend(sub_plot1_sympt)),
#                               ncol=1,nrow=3)
# grid_plot_for_paper_with_legend = plot_grid(grid_plot_for_paper, 
#                                             legends_for_paper,
#                                             nrow=1,
#                                             rel_widths = c(5,1))
# ggsave(paste0(dir$plots,'PaperPlot_3.pdf'),grid_plot_for_paper_with_legend,width=10,height=8)


# Set up dashed lines
midpointSymptImprovementNumbers = filter(filter(plot_neut_eff_struct_gmr_use_improvement, outcome==outcomes[1]),
                           GMRdiff %in% vsv_rises, 
                           from==ancestral_boost,
                           abs(baseline_sympt_eff-.5)==min(abs(baseline_sympt_eff-.5)))
midpointSevereImprovementNumbers = filter(filter(plot_neut_eff_struct_gmr_use_improvement, outcome=='Severe'), 
                            GMRdiff %in% vsv_rises, 
                            from==ancestral_boost,
                            abs(baseline_sympt_eff-.5)==min(abs(baseline_sympt_eff-.5)))

midpointSymptImprovementLine = midpointSymptImprovementNumbers %>% filter(GMRdiff == vsv_rises[1])
midpointSevereImprovementLine = midpointSevereImprovementNumbers %>% filter(GMRdiff == vsv_rises[1])
line_max_improvement = max(with(plot_neut_eff_struct_gmr_use_improvement,pop_multiplier*value),na.rm=T)

# Improvement subplot
multiplier = 100
shift = .35
sub_plot_improvement = ggplot(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rises), 
                              aes(x=baseline_sympt_eff,y=multiplier*value, colour=paste0(outcome,GMRdiff), linetype = paste0(outcome,GMRdiff)))+
  geom_line(mapping = aes(size=as.character(GMRdiff))) +
  scale_size_manual(values = vsv_sizes, guide='none')+
  theme_for_sim_plots+
  scale_colour_manual(values=improvement_cols,name=improv_legend_name,labels = improvement_labels)+
  scale_linetype_manual(values=improvement_linetypes,name=improv_legend_name, labels = improvement_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
  #scale_linetype_manual(values=vsv_linetypes,name=fold_improvement_legend_name, labels = vsv_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
  scale_x_continuous(breaks = percent_breaks,labels = percent_break_labels, expand = c(0.001,.001))+
  labs(title='Improvement From Variant Vaccine', x=pre_boost_label_text,y=percentage_point_y_title)+
  theme(strip.placement = 'outside')+
  # Vertical dashed lines
  geom_segment(data = midpointSymptImprovementLine, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=multiplier*line_max_improvement),colour='black',size=.5,linetype='dashed')+
  geom_segment(data = midpointSevereImprovementLine, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=multiplier*line_max_improvement),colour='black',size=.5,linetype='dashed')+
  # Text
  geom_text(data = midpointSymptImprovementNumbers,aes(x=GMRdiff*shift,y=multiplier*value,label=round(multiplier*value,1)), show.legend = F)+
  geom_text(data = midpointSevereImprovementNumbers,aes(x=GMRdiff*shift,y=multiplier*value,label=round(multiplier*value,1)), show.legend = F)+
  coord_cartesian(xlim = x_expand, ylim = y_improvement_expand, expand=F)
#print(sub_plot_improvement)

sub_plot_improvement_after_6m = ggplot(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rises), 
                                       aes(x=baseline_sympt_eff,y=multiplier*value6m, colour=paste0(outcome,GMRdiff), linetype = paste0(outcome,GMRdiff)))+
  geom_line(mapping = aes(size=as.character(GMRdiff))) +
  scale_size_manual(values = vsv_sizes, guide='none')+
  theme_for_sim_plots+
  scale_colour_manual(values=improvement_cols,name=improv_legend_name,labels = improvement_labels)+
  scale_linetype_manual(values=improvement_linetypes,name=improv_legend_name, labels = improvement_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
  #scale_colour_manual(values=c('darkgreen','orange'),name='Disease Outcome')+
  #scale_linetype_manual(values=vsv_linetypes,name=fold_improvement_legend_name, labels = vsv_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
  scale_x_continuous(breaks = percent_breaks,labels = percent_break_labels, expand = c(0.001,.001))+
  labs(title='Improvement From Variant Vaccine', x=pre_boost_label_text,y=percentage_point_y_title)+
  theme(strip.placement = 'outside')+
  # Vertical dashed lines
  geom_segment(data = midpointSymptImprovementLine, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=multiplier*line_max_improvement),colour='black',size=.5,linetype='dashed')+
  geom_segment(data = midpointSevereImprovementLine, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=multiplier*line_max_improvement),colour='black',size=.5,linetype='dashed')+
  # Text
  geom_text(data = midpointSymptImprovementNumbers,aes(x=GMRdiff*shift,y=multiplier*value6m,label=round(multiplier*value6m,1)), show.legend = F)+
  geom_text(data = midpointSevereImprovementNumbers,aes(x=GMRdiff*shift,y=multiplier*value6m,label=round(multiplier*value6m,1)), show.legend = F)+
  coord_cartesian(xlim = x_expand, ylim = y_improvement_expand, expand=F)

sub_plot_improvement_avg_over_6m = ggplot(filter(plot_neut_eff_struct_gmr_use_improvement,from==ancestral_boost,GMRdiff==vsv_rises), 
                                          aes(x=baseline_sympt_eff,y=multiplier*value6mavg, colour=as.character(GMRdiff), linetype = as.character(GMRdiff)))+
  geom_line(mapping = aes(size=as.character(GMRdiff))) +
  scale_size_manual(values = vsv_sizes, guide='none')+
  theme_for_sim_plots+
  scale_colour_manual(values=improvement_cols,name=improv_legend_name,labels = improvement_labels)+
  scale_linetype_manual(values=improvement_linetypes,name=improv_legend_name, labels = improvement_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
  #scale_colour_manual(values=c('darkgreen','orange'),name='Disease Outcome')+
  #scale_linetype_manual(values=vsv_linetypes,name=fold_improvement_legend_name, labels = vsv_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
  scale_x_continuous(breaks = percent_breaks,labels = percent_break_labels, expand = c(0.001,.001))+
  labs(title='Improvement From Variant Vaccine', x=pre_boost_label_text,y=percentage_point_y_title)+
  # Vertical dashed lines
  geom_segment(data = midpointSymptImprovementLine, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=multiplier*line_max_improvement),colour='black',size=.5,linetype='dashed')+
  geom_segment(data = midpointSevereImprovementLine, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=multiplier*line_max_improvement),colour='black',size=.5,linetype='dashed')+
  # Text
  geom_text(data = midpointSymptImprovementNumbers,aes(x=GMRdiff*shift,y=multiplier*value6mavg,label=round(multiplier*value6mavg,1)), show.legend = F)+
  geom_text(data = midpointSevereImprovementNumbers,aes(x=GMRdiff*shift,y=multiplier*value6mavg,label=round(multiplier*value6mavg,1)), show.legend = F)+
  theme(strip.placement = 'outside')+
  coord_cartesian(xlim = x_expand, ylim = y_improvement_expand, expand=F)
#print(sub_plot_improvement_avg_over_6m)
ggsave(paste0(dir$plots,'Now_Later_Average_Improvement.pdf'),plot_grid(sub_plot_improvement,sub_plot_improvement_after_6m, sub_plot_improvement_avg_over_6m,nrow=3), width=tall_width,height=tall_height)

sub_plot_improvement_avg_over_6m_by_outcome = sub_plot_improvement_avg_over_6m+facet_wrap(~outcome)

#print(sub_plot_improvement)


#### 4, 8, 12 week plots
new_eff_over_time_cols = which(str_detect(colnames(plot_neut_eff_struct_gmr_use),'after') | str_detect(colnames(plot_neut_eff_struct_gmr_use),'over'))
plot_neut_eff_struct_weekly = plot_neut_eff_struct_gmr_use %>%
                                melt(measure.vars = new_eff_over_time_cols,na.rm=T, value.name = 'calc_eff_value') %>%
                                mutate(eff_after_over_var = get_eff_after_over_var(variable),
                                       after_over_time = get_timing(variable)
                                       )

plot_neut_eff_struct_weekly$variable = as.character(plot_neut_eff_struct_weekly$variable)
these_inds = startsWith(plot_neut_eff_struct_weekly$variable,'baseline_eff_after')
plot_neut_eff_struct_weekly$eff_after_over_var[these_inds] = 'baseline_eff_after'
a=str_split(plot_neut_eff_struct_weekly$variable[these_inds],'_',simplify='T')
plot_neut_eff_struct_weekly$after_over_time[these_inds] = a[,ncol(a)]

these_inds = startsWith(plot_neut_eff_struct_weekly$variable,'new_eff_after')
plot_neut_eff_struct_weekly$eff_after_over_var[these_inds] = 'new_eff_after'
a=str_split(plot_neut_eff_struct_weekly$variable[these_inds],'_',simplify='T')
plot_neut_eff_struct_weekly$after_over_time[these_inds] = a[,ncol(a)]

these_inds = startsWith(plot_neut_eff_struct_weekly$variable,'avg_eff_over')&!endsWith(plot_neut_eff_struct_weekly$variable,'baseline')
plot_neut_eff_struct_weekly$eff_after_over_var[these_inds] = 'avg_eff_over'
a=str_split(plot_neut_eff_struct_weekly$variable[these_inds],'_',simplify='T')
plot_neut_eff_struct_weekly$after_over_time[these_inds] = a[,ncol(a)]

these_inds = startsWith(plot_neut_eff_struct_weekly$variable,'avg_eff_over')&endsWith(plot_neut_eff_struct_weekly$variable,'baseline')
plot_neut_eff_struct_weekly$eff_after_over_var[these_inds] = 'avg_eff_over_cf_baseline'
a=str_split(plot_neut_eff_struct_weekly$variable[these_inds],'_',simplify='T')
plot_neut_eff_struct_weekly$after_over_time[these_inds] = a[,4]

plot_neut_eff_struct_weekly$after_over_time = factor(plot_neut_eff_struct_weekly$after_over_time, labels = c('4w','8w','12w','6m','sixm'))
  
plots_over_time = ggplot(filter(plot_neut_eff_struct_weekly, eff_after_over_var %in% c('avg_eff_over')), aes(x=baseline_sympt_eff,y=pop_multiplier*calc_eff_value, colour=as.factor(GMR), linetype = as.factor(GMR)))+
    geom_ribbon(aes(ymax=pop_multiplier*calc_eff_value,fill=as.factor(GMR)),ymin=0,alpha=alpha_fill,colour=NA)+
    geom_line(mapping = aes(size=as.factor(GMR))) +
    scale_size_manual(values = gmr_sizes, guide='none')+
    theme_for_sim_plots+
    scale_colour_manual(values=gmr_colours,name=fold_improvement_legend_name, labels = gmr_labels)+
    scale_fill_manual(values=gmr_fills,name=fold_improvement_legend_name, labels = gmr_labels,guide='none', na.value = NA)+
    scale_linetype_manual(values=gmr_linetypes,name=fold_improvement_legend_name, labels = gmr_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
    scale_x_continuous(breaks = percent_breaks,labels = percent_break_labels,expand=c(0.001,0.001))+
    scale_y_continuous(breaks = percent_breaks, labels = percent_break_labels, expand=c(0,0))+
    labs(title=y_label_efficacy, x=pre_boost_label_text,y=y_label_efficacy)+
    #geom_segment(data = midpointSymptData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=pop_multiplier),colour='black',size=.5,linetype='dashed')+
    #geom_segment(data = midpointSevereData, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=pop_multiplier),colour='black',size=.5,linetype='dashed')+
    coord_cartesian(xlim = x_expand, expand=F,y=c(-.01,1.01))+
    facet_grid(after_over_time~outcome, scales = 'free_y')
  print(plots_over_time)
  ggsave(paste0(dir$plots,'ImprovedProtectionOverTimes.pdf'),plots_over_time, width=7,height=8)

  
  
  new_eff_over_time_cols = which(str_detect(colnames(plot_neut_eff_struct_improvement),'value'))
  plot_neut_eff_struct_improvement_weekly = plot_neut_eff_struct_improvement %>%
    melt(measure.vars = new_eff_over_time_cols,na.rm=T, value.name = 'calc_improvement')

  
  plot_neut_eff_struct_improvement_weekly$variable = as.character(plot_neut_eff_struct_improvement_weekly$variable)
  plot_neut_eff_struct_improvement_weekly$after_over_time = ''
  
  these_inds = which(endsWith(plot_neut_eff_struct_improvement_weekly$variable,'avg'))
  plot_neut_eff_struct_improvement_weekly$eff_after_over_var[these_inds] = 'value_avg'
  these_inds = which(plot_neut_eff_struct_improvement_weekly$variable == 'value4wavg')
  plot_neut_eff_struct_improvement_weekly$after_over_time[these_inds]='4w'
  these_inds = which(plot_neut_eff_struct_improvement_weekly$variable == 'value8wavg')
  plot_neut_eff_struct_improvement_weekly$after_over_time[these_inds]='8w'
  these_inds = which(plot_neut_eff_struct_improvement_weekly$variable == 'value12wavg')
  plot_neut_eff_struct_improvement_weekly$after_over_time[these_inds]='12w'
  these_inds = which(plot_neut_eff_struct_improvement_weekly$variable == 'value6mavg')
  plot_neut_eff_struct_improvement_weekly$after_over_time[these_inds]='6m'

  
  # This bit isnt right - only right for averages
  these_inds = which(!endsWith(plot_neut_eff_struct_improvement_weekly$variable,'avg'))
  plot_neut_eff_struct_improvement_weekly$eff_after_over_var[these_inds] = 'value'
  a=str_trunc(plot_neut_eff_struct_improvement_weekly$variable[these_inds],nchar(plot_neut_eff_struct_improvement_weekly$variable[these_inds])-3,'right',ellipsis='')
  plot_neut_eff_struct_improvement_weekly$after_over_time[these_inds] = str_trunc(a,nchar(a)-2,'left',ellipsis = '')
  
  plot_neut_eff_struct_improvement_weekly = filter(plot_neut_eff_struct_improvement_weekly,eff_after_over_var %in% c('value_avg'))
  plot_neut_eff_struct_improvement_weekly$after_over_time = factor(plot_neut_eff_struct_improvement_weekly$after_over_time, labels = c('4w','8w','12w','6m'))
  
  
# This is improvement
  improvement_over_time = ggplot(filter(plot_neut_eff_struct_improvement_weekly,
                                        from==ancestral_boost,
                                        GMRdiff==vsv_rises,
                                       eff_after_over_var %in% c('value_avg')), 
                                            aes(x=baseline_sympt_eff,y=multiplier*calc_improvement, colour=as.character(GMRdiff), linetype = as.character(GMRdiff)))+
    geom_line(mapping = aes(size=as.character(GMRdiff))) +
    scale_size_manual(values = vsv_sizes, guide='none')+
    theme_for_sim_plots+
    scale_colour_manual(values=improvement_cols,name=improv_legend_name,labels = improvement_labels)+
    scale_linetype_manual(values=improvement_linetypes,name=improv_legend_name, labels = improvement_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
    #scale_colour_manual(values=c('darkgreen','orange'),name='Disease Outcome')+
    #scale_linetype_manual(values=vsv_linetypes,name=fold_improvement_legend_name, labels = vsv_labels, guide = guide_legend(override.aes = list(fill = "#FFFFFF")))+
    scale_x_continuous(breaks = percent_breaks,labels = percent_break_labels, expand = c(0.001,.001))+
    labs(title='Improvement From Variant Vaccine', x=pre_boost_label_text,y=percentage_point_y_title)+
    # Vertical dashed lines
    #geom_segment(data = midpointSymptImprovementLine, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=multiplier*line_max_improvement),colour='black',size=.5,linetype='dashed')+
    #geom_segment(data = midpointSevereImprovementLine, aes(x=baseline_sympt_eff,y=0, xend=baseline_sympt_eff, yend=multiplier*line_max_improvement),colour='black',size=.5,linetype='dashed')+
    # Text
    #geom_text(data = midpointSymptImprovementNumbers,aes(x=GMRdiff*shift,y=multiplier*value6mavg,label=round(multiplier*value6mavg,1)), show.legend = F)+
    #geom_text(data = midpointSevereImprovementNumbers,aes(x=GMRdiff*shift,y=multiplier*value6mavg,label=round(multiplier*value6mavg,1)), show.legend = F)+
    theme(strip.placement = 'outside')+
    coord_cartesian(xlim = x_expand, ylim = y_improvement_expand, expand=F)+
    facet_grid(after_over_time~outcome, scales = 'free_y')
  print(improvement_over_time)
  ggsave(paste0(dir$plots,'Improvement_Avg_By_Time.pdf'),improvement_over_time, width=7,height=8)




