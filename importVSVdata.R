library(lifecycle)
library(data.table)
library(readxl)
library(stringr)
library(ggplot2)
library(dplyr)
library(cowplot)
library(stats)
library(reshape2)
library(ggpubr)
library(lhs)
library(mvtnorm)


# Details for importing / saving data
sheets = list(vsvdata = 'Sheet1'
              )
ranges = list(vsvdata='A1:T1000') #(changed from A1:T91 to A1:T1000 after DK check)
# Maybe should add following columns
ranges = list(vsvdata='A1:AH1000') #(changed from A1:T91 to A1:T1000 after DK check)
dir = list(base = '../',
           code = './',
           data = paste0('./'),
           plots = paste0('./Figures/'),
           manuscript_plots = './Figures/')
files = list(vsvdata = paste0(dir$data,'MasterSheet_FoldChanges_github_resubmission.xlsx'))

load(str_c(dir$data,'natmed_parameters_sympt_and_severe.RData'))

source('./WP4_modelImplementationFunctions.R')
source('./helper_functions.R')
source('./helper_functions_CIs.R')


# Plotting constants
plot_type='pdf'
plain_omicron_colour = 'tan'
variant_colours=c('Ancestral'='grey', 
                  'Omicron BA.1'='orange',
                  'Omicron BA.4/5'='orangered2', 
                  'Omicron BA.2' = plain_omicron_colour,
                  'Omicron BA.2.75' = plain_omicron_colour,
                  'Omicron BA.2.75.2' = plain_omicron_colour,
                  'Omicron BA.4.6' = plain_omicron_colour,
                  'Omicron (other)' = plain_omicron_colour,
                  'OmicronBF.7' = plain_omicron_colour,
                  'OmicronBQ.1.1' = plain_omicron_colour,
                  'OmicronXBB.1' = plain_omicron_colour,
                  'Beta'='olivedrab',
                  'Delta'='Blue','Gamma'='purple')
variant_colour_labels=variant_colours
  
valency_colours=c('Bivalent'='#7CAE00','Monovalent'='#C77CFF')
matching_text = c('Matched','Non-Matched')


read.data=T
if(read.data){
  vsvdata_read = data.frame(read_xlsx(path = str_c(files$vsvdata),sheet = sheets$vsvdata, range = ranges$vsvdata))
  vsvdata_read$Variant[vsvdata_read$Variant=='OmicronBA.1']='Omicron BA.1'
  vsvdata_read$Variant[vsvdata_read$Variant=='OmicronBA.4']='Omicron BA.4/5'
  vsvdata_read$Variant[vsvdata_read$Variant=='OmicronBA.4/5']='Omicron BA.4/5'
  vsvdata_read$Variant[vsvdata_read$Variant=='OmicronBA.4/5']='Omicron BA.4/5'
  vsvdata_read$Variant[vsvdata_read$Variant=='OmicronBA.5']='Omicron BA.4/5'
  vsvdata_read$Variant[vsvdata_read$Variant=='OmicronBA.2']='Omicron BA.2'
  vsvdata_read$Variant[vsvdata_read$Variant=='OmicronBA.2.75']='Omicron BA.2.75'
  vsvdata_read$Variant[vsvdata_read$Variant=='OmicronBA.2.75.2']='Omicron BA.2.75.2'
  vsvdata_read$Variant[vsvdata_read$Variant=='OmicronBA.4.6']='Omicron BA.4.6'
}

# Add fields to imported data
vsvdata = vsvdata_read %>% 
  filter(!is.na(Study)) %>% # remove the rows that have NA for a study as they correspond to blank rows (added after DK check)
  filter(Booster!='infection') %>% # remove infection rows
  #filter(PriorStatusGroup=='Infected') %>% #Dont include new infected data
  select(which(!startsWith(colnames(vsvdata_read),'...'))) %>%
  mutate(foldchange =	NeutsAfterBoost/NeutsBeforeBoost,
         log10foldchange = log10(foldchange),
         log10foldchangesurrogate = ifelse(is.na(log10foldchange),log10(NeutsAfterBoost),log10foldchange),
         #log10foldchangesurrogate = log10(NeutsAfterBoost),
         log10neutsafterboost = log10(NeutsAfterBoost),
         BoosterGroup1 = if_else(BoosterType=='Ancestral', 'Ancestral','Variant'),
         BoosterGroup2 = case_when(str_detect(BoosterType,'_')~'Bivalent',
                                   BoosterType =='Ancestral'~'Ancestral',
                                   T~'Monovalent'),
         PriorDoses = case_when(nPriorDoses==2~'2 Doses',
                                nPriorDoses>=3~'3+ Doses',
                                T~'Other')
  )
vsvdata$PriorDoses = factor(vsvdata$PriorDoses, levels = c('2 Doses','3+ Doses'))
vsvdata$PriorStatusGroup = factor(vsvdata$PriorStatusGroup, levels = c('Uninfected','Mixed','Infected'))
  
# The code above introduces a surrogate fold change of the neut value only IF the pre boost neuts weren't recorded.
# That is fine, as long as this is the same for ALL vaccines in that group. In the case of Pfizer 30mg 
# from their presentation it doesnt work as there is fold change data for the BA1 and Ancestral_BA1 vaccines
# but only neut data for teh Ancestral vaccine.
# for now I need to fix this manually.
fix_surrogates=with(vsvdata,which(Study==7&Variant=='Omicron BA.1'&Age=='>56y.o.'))
vsvdata$log10foldchangesurrogate[fix_surrogates]=log10(vsvdata$NeutsAfterBoost[fix_surrogates])


# Now calculate a new data frame (vvs_improvement_data) that works out the improvement for using an ancestral vaccine
# choose which improvement measure to use
#improvement_measure = 'surrogate'
improvement_measure = 'neutsafterboost'
#improvement_measure = 'foldchangecfpreboost'


#cast_formula_string = 'Study+FirstAuthor+Journal+StudyType+BoosterManufacturer+BoosterDose+Assay+PreviousVaccine+TimeSinceLastVaccine+Infection.Status+Age+Variant+DaysAfterBoost~BoosterType'
cast_formula_string = paste0('Study+FirstAuthor+Journal+StudyType',
                        '+Assay+PreviousVaccine+PriorStatusGroup+PriorDoses+AgeGroup',
                        '+Variant',
                        '+TimeSinceBoostGp+BoosterManufacturer+BoosterDose',
                              '~BoosterType')
cast_formula_string = paste0('Study+ComparisonGroup+FirstAuthor+PriorStatusGroup+PriorDoses+Variant',
                             '~BoosterType')
cast_formula = as.formula(cast_formula_string)
if (improvement_measure == 'surrogate'){
  #This uses surrogate fold change - make it explicit - are we using fold change cf pre boost, or only neuts after boosting
  vsv_improvement_data = vsvdata %>% 
    reshape2::dcast(formula=cast_formula,
          value.var = 'log10foldchangesurrogate') %>%
    measure_name = 'log10foldchangesurrogate'
    filter(!is.na(Ancestral))
} else if (improvement_measure == 'neutsafterboost'){
  # This uses only neuts after boosting
  vsv_improvement_data = vsvdata %>% 
    reshape2::dcast(formula=cast_formula,
        value.var = 'log10neutsafterboost') %>%
    filter(!is.na(Ancestral))
  measure_name = 'log10neutsafterboost'
} else if(improvement_measure == 'foldchangecfpreboost'){
  # This uses fold change in neuts cf pre boost
  vsv_improvement_data = vsvdata %>% 
    reshape2::dcast(formula=cast_formula,
          value.var = 'log10foldchange') %>%
    filter(!is.na(Ancestral))
  measure_name = 'log10foldchange'
}

# Columns for ancestral boost in new data struct
col_ancestral = which(colnames(vsv_improvement_data)=='Ancestral')
# columns for vsv boost
col_vsvs = col_ancestral+c(1:(ncol(vsv_improvement_data)-col_ancestral))

# Now calculate the improvement for the VSV vaccine
vsv_improvement_data_new = vsv_improvement_data
vsv_improvement_data_new[,col_vsvs]=vsv_improvement_data[,col_vsvs]-vsv_improvement_data[,col_ancestral]
vsv_improvement_data_new = vsv_improvement_data_new %>%
  reshape2::melt(measure.vars = col_vsvs,na.rm=T, value.name = 'log10VSVimprovement') %>%
  rename(BoosterType = variable, log10_srgt_fc_anc=Ancestral) %>%
  merge(
    reshape2::melt(vsv_improvement_data, measure.vars = col_vsvs,na.rm=T, value.name = measure_name) %>%
      rename(BoosterType = variable, srgt_fc_anc=Ancestral),
    by=c(1:(col_ancestral+1))
  ) %>%
  mutate(VSVimprovement = 10^log10VSVimprovement,
         Valency = case_when(str_detect(BoosterType,'_')~'Bivalent',
                                   T~'Monovalent'),
         VSVBoosterPart1 = str_split(BoosterType,'_',simplify=T)[,1],
         VSVBoosterPart1 = case_when(VSVBoosterPart1=='BA1' ~'Omicron BA.1',
                                     VSVBoosterPart1=='BA5' ~'Omicron BA.4/5',
                                     T~VSVBoosterPart1),
         VSVBoosterPart2 = str_split(BoosterType,'_',simplify=T)[,2],
         VSVBoosterPart2 = case_when(VSVBoosterPart2=='BA1' ~'Omicron BA.1',
                                     VSVBoosterPart2=='BA5' ~'Omicron BA.4/5',
                                     T~VSVBoosterPart2),
         Matching = ifelse((Variant==VSVBoosterPart1 | Variant ==VSVBoosterPart2),matching_text[1],matching_text[2]),
         VariantGroup = ifelse(Variant == 'Ancestral','Ancestral','Variant')
         #srgt_fc_anc=round(10^log10_srgt_fc_anc,1),
         #srgt_fc_vsv=round(10^log10_srgt_fc_vsv,1)
  )

#drop_colours=c('TRUE'='forestgreen','FALSE'='orange','NA'='grey')
study_shapes = vsvdata$Study
names(study_shapes) = vsvdata$FirstAuthor
study_shape_labels =  paste0(vsvdata$FirstAuthor,' et. al.')
study_shape_labels[study_shapes==7] = 'Pfizer (FDA powerpoint)'
study_shape_labels[study_shapes==3] = paste0(study_shape_labels[study_shapes==3],' (1)')
study_shape_labels[study_shapes==4] = paste0(study_shape_labels[study_shapes==4],' (2)')
study_shape_labels[study_shapes==5] = 'Moderna (press release)'
study_shape_labels[study_shapes==14] = 'Davis Gardiner et. al.'
study_shape_labels[study_shapes==15] = 'Pfizer (BA.5 press release)'
study_shapes[study_shapes==15]=35

vsvdata$Reference =  case_when((vsvdata$Study<6)~ vsvdata$Study,
                               (vsvdata$Study >6)~(vsvdata$Study-1),
                               T~0)
study_shape_labels_2 = paste0('(',vsvdata$Reference,')')
