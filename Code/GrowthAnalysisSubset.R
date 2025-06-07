
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




##From Size Class R script

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
str(Growth_Data_forR)

#set contrasts ALWAYS RUN
options(contrasts = c("contr.treatment","contr.poly")) #could also be contr.sum for equal groups
getOption("contrasts")

##Subset for phase 1 treatments
control <- Growth_Data_forR_full[Growth_Data_forR_full$Phase_1_treat == "Cont", ]
View(control)

both <- Growth_Data_forR_full[Growth_Data_forR_full$Phase_1_treat == "Both", ]
View(both)

hyp <- Growth_Data_forR_full[Growth_Data_forR_full$Phase_1_treat == "Hyp", ]
View(hyp)

warm <- Growth_Data_forR_full[Growth_Data_forR_full$Phase_1_treat == "Warm", ]
View(warm)

#control Phase 1, initial phase 2 stats
min(control$Actual_shell_pre_mg) #121.6519
max(control$Actual_shell_pre_mg) #1443.633
mean(control$Actual_shell_pre_mg) #415.483
sd(control$Actual_shell_pre_mg)/sqrt(length((control$Actual_shell_pre_mg))) #se #12.7764
ggplot(control) +
  aes(x = Actual_shell_pre_mg, fill = Phase_1_rep_R) +
  geom_histogram(bins = 30L) +
  scale_fill_hue(direction = 1) +
  theme_classic()
#which rep has the highest mean
Control1 <-control[control$Phase_1_rep == "01", ] #second highest
mean(Control1$Actual_shell_pre_mg) #446.6151
Control2 <-control[control$Phase_1_rep == "02", ]
mean(Control2$Actual_shell_pre_mg) #365.2395
Control3 <-control[control$Phase_1_rep == "03", ]
mean(Control3$Actual_shell_pre_mg) #392.7069
Control4 <-control[control$Phase_1_rep == "04", ]
mean(Control4$Actual_shell_pre_mg) #358.0041
Control5 <-control[control$Phase_1_rep == "05", ]
mean(Control5$Actual_shell_pre_mg) #432.2688
Control6 <-control[control$Phase_1_rep == "06", ]
mean(Control6$Actual_shell_pre_mg) #497.0167 EXCLUDE CONT06

# Calculate mean and standard error for each replicate
colnames(control)
summary_stats_cont <- control %>%
  filter(Phase_1_rep %in% c("01", "02", "03", "04", "05", "06")) %>%
  group_by(Phase_1_rep) %>%
  summarise(
    mean_dryweight_pre = mean(Dry_weight_pre, na.rm = TRUE),
    se_dryweight_pre = sd(Dry_weight_pre, na.rm = TRUE) / sqrt(n()))

# Rename replicates for plotting (e.g., "Control1", "Control2", ...)
summary_stats_cont$Replicate <- factor(summary_stats_cont$Phase_1_rep,
                                       levels = c("01", "02", "03", "04", "05", "06"),
                                       labels = c("Control1", "Control2", "Control3", "Control4", "Control5", "Control6"))

# Plot means with standard error bars
ggplot(summary_stats_cont, aes(x = Replicate, y = mean_dryweight_pre)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_dryweight_pre - se_dryweight_pre, 
                    ymax = mean_dryweight_pre + se_dryweight_pre), 
                width = 0.2) +
  labs(x = "Tank Replicate",
       y = "Mean Shell Mass (mg)") +
  theme_classic()



#both Phase 1, initial phase 2 stats
min(both$Actual_shell_pre_mg) #96.6342
max(both$Actual_shell_pre_mg) #1263.758
mean(both$Actual_shell_pre_mg) #352.8925
sd(both$Actual_shell_pre_mg)/sqrt(length((both$Actual_shell_pre_mg))) #se #12.20958
ggplot(both) +
  aes(x = Actual_shell_pre_mg, fill = Phase_1_rep_R) +
  geom_histogram(bins = 30L) +
  scale_fill_hue(direction = 1) +
  theme_classic()
#which rep has the lowest mean
Both1 <-both[both$Phase_1_rep == "01", ]
mean(Both1$Actual_shell_pre_mg) #365.007
Both2 <-both[both$Phase_1_rep == "02", ]
mean(Both2$Actual_shell_pre_mg) #319.7146 EXCLUDE BOTH02
Both3 <-both[both$Phase_1_rep == "03", ]
mean(Both3$Actual_shell_pre_mg) #386.9968
Both4 <-both[both$Phase_1_rep == "04", ]
mean(Both4$Actual_shell_pre_mg) #343.8161
Both5 <-both[both$Phase_1_rep == "05", ]
mean(Both5$Actual_shell_pre_mg) #346.6168
Both6 <-both[both$Phase_1_rep == "06", ]
mean(Both6$Actual_shell_pre_mg) #347.4679

# Calculate mean and standard error for each replicate
summary_stats_both <- both %>%
  filter(Phase_1_rep %in% c("01", "02", "03", "04", "05", "06")) %>%
  group_by(Phase_1_rep) %>%
  summarise(
    mean_dryweight_pre = mean(Dry_weight_pre, na.rm = TRUE),
    se_dryweight_pre = sd(Dry_weight_pre, na.rm = TRUE) / sqrt(n()))

# Rename replicates for plotting (e.g., "Control1", "Control2", ...)
summary_stats_both$Replicate <- factor(summary_stats_both$Phase_1_rep,
                                       levels = c("01", "02", "03", "04", "05", "06"),
                                       labels = c("Both1", "Both2", "Both3", "Both4", "Both5", "Both6"))

# Plot means with standard error bars
ggplot(summary_stats_both, aes(x = Replicate, y = mean_dryweight_pre)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_dryweight_pre - se_dryweight_pre, 
                    ymax = mean_dryweight_pre + se_dryweight_pre), 
                width = 0.2) +
  labs(x = "Tank Replicate",
       y = "Mean Shell Mass (mg)") +
  theme_classic()



#hyp Phase 1, initial phase 2 stats
min(hyp$Actual_shell_pre_mg) #71.77482
max(hyp$Actual_shell_pre_mg) #1300.652
mean(hyp$Actual_shell_pre_mg) #355.7689
sd(hyp$Actual_shell_pre_mg)/sqrt(length((hyp$Actual_shell_pre_mg))) #se #11.33287
ggplot(hyp) +
  aes(x = Actual_shell_pre_mg, fill = Phase_1_rep_R) +
  geom_histogram(bins = 30L) +
  scale_fill_hue(direction = 1) +
  theme_classic()

#which rep has the lowest mean
Hyp1 <-hyp[hyp$Phase_1_rep == "01", ]
mean(Hyp1$Actual_shell_pre_mg) #348.7181
Hyp2 <-hyp[hyp$Phase_1_rep == "02", ]
mean(Hyp2$Actual_shell_pre_mg) #373.0695 
Hyp3 <-hyp[hyp$Phase_1_rep == "03", ]
mean(Hyp3$Actual_shell_pre_mg) #393.3238
Hyp4 <-hyp[hyp$Phase_1_rep == "04", ]
mean(Hyp4$Actual_shell_pre_mg) #353.2899
Hyp5 <-hyp[hyp$Phase_1_rep == "05", ]
mean(Hyp5$Actual_shell_pre_mg) #328.0052 EXCLUDE HYP05
Hyp6 <-hyp[hyp$Phase_1_rep == "06", ]
mean(Hyp6$Actual_shell_pre_mg) #336.2902

# Calculate mean and standard error for each replicate
summary_stats_hyp <- hyp %>%
  filter(Phase_1_rep %in% c("01", "02", "03", "04", "05", "06")) %>%
  group_by(Phase_1_rep) %>%
  summarise(
    mean_dryweight_pre = mean(Dry_weight_pre, na.rm = TRUE),
    se_dryweight_pre = sd(Dry_weight_pre, na.rm = TRUE) / sqrt(n()))

# Rename replicates for plotting (e.g., "Control1", "Control2", ...)
summary_stats_hyp$Replicate <- factor(summary_stats_hyp$Phase_1_rep,
                                      levels = c("01", "02", "03", "04", "05", "06"),
                                      labels = c("Hyp1", "Hyp2", "Hyp3", "Hyp4", "Hyp5", "Hyp6"))

# Plot means with standard error bars
ggplot(summary_stats_hyp, aes(x = Replicate, y = mean_dryweight_pre)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_dryweight_pre - se_dryweight_pre, 
                    ymax = mean_dryweight_pre + se_dryweight_pre), 
                width = 0.2) +
  labs(x = "Tank Replicate",
       y = "Mean Shell Mass (mg)") +
  theme_classic()


#warm Phase 1, initial phase 2 stats
min(warm$Actual_shell_pre_mg) #93.15072
max(warm$Actual_shell_pre_mg) #1531.67
mean(warm$Actual_shell_pre_mg) #385.5639
sd(warm$Actual_shell_pre_mg)/sqrt(length((warm$Actual_shell_pre_mg))) #se #11.57248
ggplot(warm) +
  aes(x = Actual_shell_pre_mg, fill = Phase_1_rep_R) +
  geom_histogram(bins = 30L) +
  scale_fill_hue(direction = 1) +
  theme_classic()

#which rep has the highest mean
Warm1 <-warm[warm$Phase_1_rep == "01", ]
mean(Warm1$Actual_shell_pre_mg) #436.8607 EXLUDE WARM01
Warm2 <-warm[warm$Phase_1_rep == "02", ]
mean(Warm2$Actual_shell_pre_mg) #416.0262 
Warm3 <-warm[warm$Phase_1_rep == "03", ]
mean(Warm3$Actual_shell_pre_mg) #345.7839
Warm4 <-warm[warm$Phase_1_rep == "04", ]
mean(Warm4$Actual_shell_pre_mg) #369.299
Warm5 <-warm[warm$Phase_1_rep == "05", ]
mean(Warm5$Actual_shell_pre_mg) #370.842 
Warm6 <-warm[warm$Phase_1_rep == "06", ]
mean(Warm6$Actual_shell_pre_mg) #372.667

# Calculate mean and standard error for each replicate
summary_stats_warm <- warm %>%
  filter(Phase_1_rep %in% c("01", "02", "03", "04", "05", "06")) %>%
  group_by(Phase_1_rep) %>%
  summarise(
    mean_dryweight_pre = mean(Dry_weight_pre, na.rm = TRUE),
    se_dryweight_pre = sd(Dry_weight_pre, na.rm = TRUE) / sqrt(n()))

# Rename replicates for plotting (e.g., "Control1", "Control2", ...)
summary_stats_warm$Replicate <- factor(summary_stats_warm$Phase_1_rep,
                                       levels = c("01", "02", "03", "04", "05", "06"),
                                       labels = c("Warm1", "Warm2", "Warm3", "Warm4", "Warm5", "Warm6"))

# Plot means with standard error bars
ggplot(summary_stats_warm, aes(x = Replicate, y = mean_dryweight_pre)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_dryweight_pre - se_dryweight_pre, 
                    ymax = mean_dryweight_pre + se_dryweight_pre), 
                width = 0.2) +
  labs(x = "Tank Replicate",
       y = "Mean Shell Mass (mg)") +
  theme_classic()



##Only using the mean +/- one standard deviation from shell mass pre.

mean(Growth_Data_forR_clean$Actual_shell_pre_mg) #mean
sd(Growth_Data_forR_clean$Actual_shell_pre_mg)/sqrt(length((Growth_Data_forR_clean$Actual_shell_pre_mg))) #se
sd(Growth_Data_forR_clean$Actual_shell_pre_mg)

Growth_Data_forR_test <- Growth_Data_forR %>%
  filter(!is.na(Actual_tissue_growth_mg) & Actual_shell_pre_mg < 581.4607 & Actual_shell_pre_mg > 172.6935)
View(Growth_Data_forR_test)
nrow(Growth_Data_forR_test)

#mean and standard error calculations
377.0771+204.3836
377.0771-204.3836

#Run model on the middle chunk of data to see what happens
leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, Growth_Data_forR_clean) #passes
T1 <- residuals(Tlme1) 
qqnorm(T1) #looks good



##From Growth ANalysis outliers

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


