#control tissue mean & stdev
#load packages
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
library(plotrix)

#Read in dataset, set column types
getwd()
setwd("/Users/sophiemontague/Desktop/MontagueORCC/Oyster Weight Data")

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

