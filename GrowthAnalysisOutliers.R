#packages
library(lme4)
library(readr)
library (car)
library(emmeans)
library(boot)
library(esquisse)
library(dplyr)
library(nlme)
library(car)
library(ggplot2)
library(plotrix) #for standard error

#data
Growth_Data_forR_1 <- read_csv("/Users/sophiemontague/Desktop/MontagueORCC/Oyster Weight Data/growth_phase2.1_test.csv", 
                             col_types = cols(Phase_1_temp = col_factor(), 
                                              Phase_1_DO =col_factor(), 
                                              Phase_2.1_temp =col_factor(), 
                                              Phase_2.1_DO=col_factor(), 
                                              Phase_1_rep =col_factor(), 
                                              Phase_2_rep =col_factor(),
                                              Ratio_tissue_shell_mg = col_double(),
                                              Phase1_Phase2_rep = col_factor(),
                                              Phase_1_treat = col_factor(),
                                              Phase_2_treat = col_factor(),
                                              Phase_2_rep_R = col_factor(),
                                              Phase_1_rep_R = col_factor()))


str(Growth_Data_forR_1) #check data was read in correctly

#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups, contr.poly
getOption("contrasts") 

#filter out rows with NA in Actual_tissue_growth_mg only to not include data from dead oysters in analysis
#also filter out one replicate for each
Growth_Data_forR_clean_1 <- Growth_Data_forR[!is.na(Growth_Data_forR_1$Actual_tissue_growth_mg), ]
View(Growth_Data_forR_clean_1)

#visualize data with outliers removed
ggplot(Growth_Data_forR_clean_1) +
  aes(x = Actual_shell_growth_mg, y = Actual_tissue_growth_mg) +
  geom_point(colour = "#112446") +
  theme_classic()

##Subset for phase 1 treatments and find means for shell pre
control <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_1_treat == "Cont", ]
View(control)
mean(control$Actual_shell_pre_mg) #415.483
    #to see how the means change from removing outliers
control_1 <- Growth_Data_forR_clean_1[Growth_Data_forR_clean_1$Phase_1_treat == "Cont", ]
View(control_1)
mean(control_1$Actual_shell_pre_mg) #416.2745

ggplot(control_1) +
  aes(x = Actual_shell_pre_mg, y = Phase_1_rep_R) +
  geom_boxplot(fill = "#112446") +
  theme_minimal() #exclude cont 06

both <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_1_treat == "Both", ]
View(both)
mean(both$Actual_shell_pre_mg) #352.8925

both_1 <- Growth_Data_forR_clean_1[Growth_Data_forR_clean_1$Phase_1_treat == "Both", ]
View(both)
mean(both_1$Actual_shell_pre_mg) #352.3444

ggplot(both_1) +
  aes(x = Actual_shell_pre_mg, y = Phase_1_rep_R) +
  geom_boxplot(fill = "#112446") +
  theme_minimal() # exclude both 02

hyp <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_1_treat == "Hyp", ]
View(hyp)
mean(hyp$Actual_shell_pre_mg) #355.7689

hyp_1 <- Growth_Data_forR_clean_1[Growth_Data_forR_clean_1$Phase_1_treat == "Hyp", ]
View(hyp)
mean(hyp_1$Actual_shell_pre_mg) #355.524

ggplot(hyp_1) +
  aes(x = Actual_shell_pre_mg, y = Phase_1_rep_R) +
  geom_boxplot(fill = "#112446") +
  theme_minimal() # exclude hyp 01

warm <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_1_treat == "Warm", ]
View(warm)
mean(warm$Actual_shell_pre_mg) #385.5639

warm_1 <- Growth_Data_forR_clean_1[Growth_Data_forR_clean_1$Phase_1_treat == "Warm", ]
View(warm)
mean(warm$Actual_shell_pre_mg) #385.5639

ggplot(warm_1) +
  aes(x = Actual_shell_pre_mg, y = Phase_1_rep_R) +
  geom_boxplot(fill = "#112446") +
  theme_minimal() # exclude warm 02

