
#packages
library(lme4) #for linear mixed effects model
library(readr) #to read in data
library(emmeans) #for post hocs
library(esquisse) #interface for building plots
library(dplyr) #for data wrangling
library(nlme) #for linear mixed effects model
library(car) #for the Levene test
library(ggplot2) #for graphs
library(plotrix) #for standard error

#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups
getOption("contrasts") 


#subset for phase 2 treatments: Effect of phase 1 treatment on oysters in each phase 2 treatment:
  #Subset into phase 2 treatments, look at phase 1 DO and temp as fixed variables

control_df <- Growth_Data_forR_full[Growth_Data_forR_full$Phase_2_treat == "Cont", ]
View(control_df)

both2_df <- Growth_Data_forR_full[Growth_Data_forR_full$Phase_2_treat == "Both", ]
View(both_df)

hyp2_df <- Growth_Data_forR_full[Growth_Data_forR_full$Phase_2_treat == "Hyp", ]
View(hyp_df)

warm2_df <- Growth_Data_forR_full[Growth_Data_forR_full$Phase_2_treat == "Warm", ]
View(warm_df)

#Phase 2 control
#tissue
Clme1 <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_tissue_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = control2_df)
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
Clme4 <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_tissue_pre_mg,
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
Clme7 <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_tissue_pre_mg,
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
Clme10 <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_tissue_pre_mg,
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






##Subset for phase 1 treatments, How did the oysters in each phase 1 treatment respond to phase 2.
control <- Growth_Data_forR_clean[Growth_Data_forR_clean$Phase_1_treat == "Cont", ]
View(control)

  #etc.... more in "GrowthAnalysisSizeClass.R"


#Phase 1 control
#tissue
ClmeA <- lme(Actual_tissue_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_tissue_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = control)
anova(ClmeA, type = "m")


leveneTest(Actual_tissue_pre_mg~Phase1_Phase2_treat, control) #passes
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


#Phase 1 both
#tissue
ClmeD <- lme(Actual_tissue_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_tissue_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = both)
anova(ClmeD, type = "m")

leveneTest(Actual_tissue_pre_mg~Phase1_Phase2_treat, both) #passes
CTA <- residuals(ClmeA) 
qqnorm(CTA) #looks good

#shell
ClmeE <- lme(Actual_shell_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = both)
anova(ClmeE, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, both) #passes
CSB <- residuals(ClmeB) 
qqnorm(CSB) #looks good

#tissue:shell
ClmeF <- lme(Ratio_tissue_shell_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = both)
anova(ClmeF, type = "m")

#Phase 1 hyp
#tissue
ClmeG <- lme(Actual_tissue_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_tissue_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = hyp)
anova(ClmeG, type = "m")

leveneTest(Actual_tissue_pre_mg~Phase1_Phase2_treat, hyp) #passes
CTA <- residuals(ClmeA) 
qqnorm(CTA) #looks good

#shell
ClmeH <- lme(Actual_shell_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = hyp)
anova(ClmeH, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, hyp) #passes
CSB <- residuals(ClmeB) 
qqnorm(CSB) #looks good

#tissue:shell
ClmeI <- lme(Ratio_tissue_shell_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = hyp)
anova(ClmeI, type = "m")

leveneTest(Ratio_tissue_shell_mg~Phase1_Phase2_treat, hyp) #passes
CSB <- residuals(ClmeB) 
qqnorm(CSB) #looks good


#Phase 1 warm
#tissue
ClmeJ <- lme(Actual_tissue_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_tissue_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = warm)
anova(ClmeJ, type = "m")

leveneTest(Actual_tissue_pre_mg~Phase1_Phase2_treat, warm) #passes
CTA <- residuals(ClmeA) 
qqnorm(CTA) #looks good

#shell
ClmeK <- lme(Actual_shell_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = warm)
anova(ClmeK, type = "m")

leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, warm) #passes
CSB <- residuals(ClmeB) 
qqnorm(CSB) #looks good

#tissue:shell
ClmeL <- lme(Ratio_tissue_shell_mg ~ Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg,
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = warm)
anova(ClmeL, type = "m")

leveneTest(Ratio_tissue_shell_mg~Phase1_Phase2_treat, warm) #passes
CSB <- residuals(ClmeB) 
qqnorm(CSB) #looks good
