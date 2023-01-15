Analysis and plots for this project can be performed by running the file 
make_all_plots_for_paper.R

This will: 
1) import the data, ('./importVSVdata.R')
2) make the fold change ('./make_foldchange_plots.R') and improvements plots ('./make_improvement_plots.R')
3) generate the simulated data ('./make_simulated_dataset.R') - Note this can take up to 5 mins on a standard computer
4) make the simulation plots './make_efficacydifference_plots.R') 
5) Generate the confidence interval data - note this can take 1-2 days on standard computer ('./simulate_CIs.R')
At the moment the iterations are set not to run, and they would be loaded from a file (if the file existed - which it doesn't online as it is too large). 
If the variables 
make.eff_avged and make.new.params on lines 17 and 18 are set to true, then the bootstrapping will be done again.
6) finalise the plots for the manuscript ('./make_final_plots_for_paper.R')
and
7) run a file to generate some of the output text for the manuscript ('./make_text_for_paper.R') 
8) run a file to check the p-values for the manuscript ('./check_p_values.R') 

The following libraries will need to be installed: readxl, stringr, ggplot2, dplyr, cowplot, stats, reshape2, ggpubr, RColorBrewer, data.table
