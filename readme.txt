Analysis and plots for this project can be performed by running the file 
make_all_plots_for_paper.R

This will: 
1) import the data, ('./importVSVdata.R')
2) make the fold change ('./make_foldchange_plots.R') and improvements plots ('./make_improvement_plots.R')
3) generate the simulated data ('./make_simulated_dataset.R') - Note this can take up to 5 mins on a standard computer
4) make the simulation plots './make_efficacydifference_plots.R') 
5) finalise the plots for the manuscript ('./make_final_plots_for_paper.R')
and
6) run a file to generate some of the output text for the manuscript ('./make_text_for_paper.R') 

The following libraries will need to be installed: readxl, stringr, ggplot2, dplyr, cowplot, stats, reshape2, ggpubr, RColorBrewer
