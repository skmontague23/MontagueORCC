
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
Growth_Data_forR <- read_csv("/Users/sophiemontague/Desktop/MontagueORCC/Oyster Weight Data/growth_phase2.1_weightsSKM.csv", 
                             col_types = cols(Phase_1_temp = col_factor(), 
                                              Phase_1_DO =col_factor(), 
                                              Phase_2.1_temp =col_factor(), 
                                              Phase_2.1_DO=col_factor(), 
                                              Phase_1_rep =col_factor(), 
                                              Phase_2_rep =col_factor(),
                                              Ratio_tissue_shell_mg = col_double(),
                                              Phase1_Phase2_rep = col_factor(),
                                              Phase_1_treat = col_factor(),
                                              Phase_2_treat = col_factor()))

str(Growth_Data_forR) #check data was read in correctly

#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups
getOption("contrasts") 

#weights
vf_1 <- varIdent(form = ~1 | Phase_1_treat) #test running separately
vf <- varIdent(form = ~1 | Phase1_Phase2_treat)

#filter out rows with NA in Actual_tissue_growth_mg only to not include data from dead oysters in analysis
Growth_Data_forR_clean <- Growth_Data_forR[!is.na(Growth_Data_forR$Actual_tissue_growth_mg), ]
View(Growth_Data_forR_clean)

#subset for phase 2 treatments
control_df <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_2_treat == "Cont", ]
View(control_df)

both_df <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_2_treat == "Both", ]
View(both_df)

hyp_df <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_2_treat == "Hyp", ]
View(hyp_df)

warm_df <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_2_treat == "Warm", ]
View(warm_df)

#Phase 2 control
#tissue
Clme1 <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = control_df)
anova(Clme1, type = "m")

leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, control_df) #passes
CT <- residuals(Clme1) 
qqnorm(CT) #looks good

#shell
Clme2 <- lme(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = control_df)
anova(Clme2, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, control_df) #passes
CS <- residuals(Clme2) 
qqnorm(CS) #looks good

#tissue:shell
Clme3 <- lme(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = control_df)
anova(Clme3, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, control_df) #passes
CTS <- residuals(Clme3) 
qqnorm(CTS) #looks good

#Phase 2 Both
#tissue
Clme4 <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = both_df)
anova(Clme4, type = "m")

leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, both_df) #passes
BT <- residuals(Clme4) 
qqnorm(BT) #looks good

#shell
Clme5 <- lme(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             weights = varIdent(vf_1), 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = both_df)
anova(Clme5, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, both_df) #passes
BS <- residuals(Clme5) 
qqnorm(BS) #looks good

#tissue:shell
Clme6 <- lme(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = both_df)
anova(Clme6, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, control_df) #passes
BTS <- residuals(Clme6) 
qqnorm(BTS) #looks good

#Phase 2 Hyp

#tissue
Clme7 <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = hyp_df)
anova(Clme7, type = "m")

leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, hyp_df) #passes
HT <- residuals(Clme7) 
qqnorm(HT) #looks good

#shell
Clme8 <- lme(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = hyp_df)
anova(Clme8, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, hyp_df) #passes
HS <- residuals(Clme8) 
qqnorm(HS) #looks good

#tissue:shell
Clme9 <- lme(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = hyp_df)
anova(Clme9, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, control_df) #passes
CTS <- residuals(Clme3) 
qqnorm(CTS) #looks good

#Phase 2 Warm

#tissue
Clme10 <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = warm_df)
anova(Clme10, type = "m")

leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, warm_df) #passes
WT <- residuals(Clme10) 
qqnorm(WT) #looks good

#shell
Clme11 <- lme(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = warm_df)
anova(Clme11, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, warm_df) #passes
WS <- residuals(Clme11) 
qqnorm(WS) #looks good

#tissue:shell
Clme12 <- lme(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = warm_df)
anova(Clme12, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, warm_df) #passes
WTS <- residuals(Clme12) 
qqnorm(WTS) #looks good






##Subset for phase 1 treatments
control <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_1_treat == "Cont", ]
View(control)

  #etc.... more in "GrowthAnalysisSizeClass.R"

#Phase 1 control
#tissue
ClmeA <- lme(Actual_tissue_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = control)
anova(ClmeA, type = "m")

leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, control) #passes
CTA <- residuals(ClmeA) 
qqnorm(CTA) #looks good

#shell
ClmeB <- lme(Actual_shell_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = control)
anova(ClmeB, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, control) #passes
CSB <- residuals(ClmeB) 
qqnorm(CSB) #looks good

#tissue:shell
ClmeC <- lme(Ratio_tissue_shell_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = control)
anova(ClmeC, type = "m")




