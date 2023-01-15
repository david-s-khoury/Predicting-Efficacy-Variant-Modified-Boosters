library(readxl)
library(stringr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(stats)
library(reshape2)
library(ggpubr)
library(RColorBrewer)

total_starttime = Sys.time()

source('./importVSVdata.R')
source('./make_foldchange_plots.R')
source('./make_improvement_plots.R')

sim_starttime=Sys.time()
source('./make_simulated_dataset.R')
sim_endtime=Sys.time()
sim_runtime=sim_endtime-sim_starttime
print(paste0('make_simulated_dataset.R: RunTime = ',format(round(sim_runtime,2))))

source('./make_efficacydifference_plots.R')
source('./simulate_CIs.R')
source('./make_final_plots_for_paper.R')
source('./make_text_for_paper.R')
source('./check_p_values')

total_endtime=Sys.time()
total_runtime=total_endtime-total_starttime
print(paste0('make_all_plots_for_paper.R: RunTime = ',format(round(total_runtime,2))))
