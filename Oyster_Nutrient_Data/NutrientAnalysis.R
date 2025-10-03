#load packages
library(lme4) #for linear mixed effects model
library(readr) #to read in data
library(tidyr)
library(emmeans) #for post hocs
library(esquisse) #interface for building plots
library(dplyr) #for data wrangling
library(nlme) #for linear mixed effects model
library(car) #for the Levene test
library(ggplot2) #for graphs
library(plotrix) #for standard error
library(performance) #model testing
library(lmtest) #another option than levene test for variance

#Read in dataset, set column types
getwd()
setwd("/Users/sophiemontague/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Nutrient_Data")


#### Tissue N and C ####
Nutrient_df <- read_csv("Phase1_nutrient_working.csv", 
                        col_types = cols(wt_percent_N = col_number(), 
                                         wt_percentC = col_number()))

Nutrient_dfclean <- Nutrient_df%>%
  filter(N_exclude != "Y" | is.na(N_exclude))


#plot with mean and SD, NITROGEN
summary_stats <- Nutrient_dfclean %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_growth = mean(wt_percent_N, na.rm = TRUE),
    se_growth = std.error(wt_percent_N, na.rm = TRUE))
View(summary_stats)

unique(Nutrient_dfclean$Phase_1_treat)

#reorder Phase_1_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "% N") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend



#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 


## % N
#nothing significant
Nm1 <- lmer(wt_percent_N ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Nutrient_dfclean, REML=TRUE)
Anova(Nm1, test="F", type="III")

#post hocs
emmeans(Tm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#diagnostics
leveneTest(wt_percent_N~Phase_1_treat, Nutrient_dfclean)
m1.e <- residuals(Nm1) 
qqnorm(m1.e)
qqline(m1.e)



## % C
#nothing significant
Cm1 <- lmer(wt_percent_C ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Nutrient_dfclean, REML=TRUE)
Anova(Cm1, test="F", type="III")

#post hocs
emmeans(Cm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#diagnostics
leveneTest(wt_percent_C~Phase_1_treat, Nutrient_dfclean)
m1.e <- residuals(Cm1) 
qqnorm(m1.e)
qqline(m1.e)




#Calculate summary stats for Nitrogen to multiply with mass
summary_stats_N <- Nutrient_dfclean %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_N, na.rm = TRUE),
    se_growth = std.error(wt_percent_N, na.rm = TRUE))
View(summary_stats_N)

#Now C
summary_stats_C <- Nutrient_dfclean %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_C, na.rm = TRUE),
    se_growth = std.error(wt_percent_C, na.rm = TRUE))
View(summary_stats_C)




##Apply percentages to weights for PHASE 1
Nutrientsbyweight <- read_csv("/Users/sophiemontague/Desktop/Nutrients_byweight.csv", 
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

colnames(Nutrientsbyweight)
N_C_byweight <- Nutrientsbyweight %>%
  filter(Exclude_all != "Y" | is.na(Exclude_all)) %>%
  mutate(
    Phase2N_tissue = case_when(
      Phase_1_treat == "Cont" ~ Actual_tissue_pre_mg * 0.07741667,
      Phase_1_treat == "Warm" ~ Actual_tissue_pre_mg * 0.07421429,
      Phase_1_treat == "Hyp"  ~ Actual_tissue_pre_mg * 0.08007692,
      Phase_1_treat == "Both" ~ Actual_tissue_pre_mg * 0.06581818,
      TRUE ~ NA_real_
    ))%>%
    mutate(
      Phase2C_tissue = case_when(
        Phase_1_treat == "Cont" ~ Actual_tissue_pre_mg * 0.3853333,
        Phase_1_treat == "Warm" ~ Actual_tissue_pre_mg * 0.3741429,
        Phase_1_treat == "Hyp"  ~ Actual_tissue_pre_mg * 0.4050769,
        Phase_1_treat == "Both" ~ Actual_tissue_pre_mg * 0.3490909,
        TRUE ~ NA_real_)
  )

View(N_C_byweight)


## % N by weight
Nm1 <- lmer(log(Phase2N_tissue) ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = N_C_byweight, REML=TRUE)
Anova(Nm1, test="F", type="III")

#post hocs
emmeans(Nm1,specs = pairwise ~ Phase_1_DO, adjust = "none")
emmeans(Nm1,specs = pairwise ~ Phase_1_temp, adjust = "none")

#diagnostics
leveneTest(log(Phase2N_tissue)~Phase_1_treat, N_C_byweight)
m1.e <- residuals(Nm1) 
qqnorm(m1.e)
qqline(m1.e)


#plot with mean and SD, NITROGEN
summary_stats <- N_C_byweight %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_growth = mean(Phase2N_tissue, na.rm = TRUE),
    se_growth = std.error(Phase2N_tissue, na.rm = TRUE))
View(summary_stats)

unique(Nutrient_dfclean$Phase_1_treat)

#reorder Phase_1_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Nitrogen per weight (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend


## % C by weight
Cm1 <- lmer(log(Phase2C_tissue) ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = N_C_byweight, REML=TRUE)
Anova(Cm1, test="F", type="III")

#post hocs
emmeans(Cm1,specs = pairwise ~ Phase_1_DO, adjust = "none")
emmeans(Cm1,specs = pairwise ~ Phase_1_temp, adjust = "none")

#diagnostics
leveneTest(log(Phase2C_tissue)~Phase_1_treat, N_C_byweight)
m1.e <- residuals(Cm1) 
qqnorm(m1.e)
qqline(m1.e)


#plot with mean and SD, CARBON
summary_stats <- N_C_byweight %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_growth = mean(Phase2C_tissue, na.rm = TRUE),
    se_growth = std.error(Phase2C_tissue, na.rm = TRUE))
View(summary_stats)

unique(Nutrient_dfclean$Phase_1_treat)

#reorder Phase_1_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Carbon per weight (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend





####Shell N and C ####
Shell_nutrient_df <- read_csv("Phase1_shellnutrient_working.csv", 
                        col_types = cols(wt_percent_N = col_number(), 
                                         wt_percent_C = col_number()))

#Shell_nutrient_dfclean <- Shell_nutrient_df%>%
#  filter(N_exclude != "Y" | is.na(N_exclude))


#plot with mean and SD, NITROGEN
summary_stats_s <- Shell_nutrient_df %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_growth = mean(wt_percent_N, na.rm = TRUE),
    se_growth = std.error(wt_percent_N, na.rm = TRUE))
View(summary_stats_s)

unique(Nutrient_dfclean$Phase_1_treat)

#reorder Phase_1_treat
summary_stats_s$Phase_1_treat <- factor(summary_stats_s$Phase_1_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "% N") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend



#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 


## % N
#nothing significant
Nm1 <- lmer(wt_percent_N ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Shell_nutrient_df, REML=TRUE)
Anova(Nm1, test="F", type="III")

#post hocs
emmeans(Nm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#diagnostics
leveneTest(wt_percent_N~Phase_1_treat, Shell_nutrient_df)
m1.e <- residuals(Nm1) 
qqnorm(m1.e)
qqline(m1.e)



## % C
#nothing significant
Cm1 <- lmer(wt_percent_C ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Shell_nutrient_df, REML=TRUE)
Anova(Cm1, test="F", type="III")

#post hocs
emmeans(Cm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#diagnostics
leveneTest(wt_percent_C~Phase_1_treat, Shell_nutrient_df)
m1.e <- residuals(Cm1) 
qqnorm(m1.e)
qqline(m1.e)




#Calculate summary stats for Nitrogen to multiply with mass
summary_stats_N <- Shell_nutrient_df %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_N, na.rm = TRUE),
    se_growth = std.error(wt_percent_N, na.rm = TRUE))
View(summary_stats_N)

#Now C
summary_stats_C <- Nutrient_dfclean %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_C, na.rm = TRUE),
    se_growth = std.error(wt_percent_C, na.rm = TRUE))
View(summary_stats_C)




##Apply percentages to weights for PHASE 1
Shell_Nutrientsbyweight <- read_csv("/Users/sophiemontague/Desktop/Nutrients_byweight.csv", 
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

colnames(Nutrientsbyweight)
N_C_byweight <- Shell_Nutrientsbyweight %>%
  filter(Exclude_all != "Y" | is.na(Exclude_all)) %>%
  mutate(
    Phase2N_shell = case_when(
      Phase_1_treat == "Cont" ~ Actual_shell_pre_mg * 0.07741667,
      Phase_1_treat == "Warm" ~ Actual_shell_pre_mg * 0.07421429,
      Phase_1_treat == "Hyp"  ~ Actual_shell_pre_mg * 0.08007692,
      Phase_1_treat == "Both" ~ Actual_shell_pre_mg * 0.06581818,
      TRUE ~ NA_real_
    ))%>%
  mutate(
    Phase2C_shell = case_when(
      Phase_1_treat == "Cont" ~ Actual_shell_pre_mg * 0.3853333,
      Phase_1_treat == "Warm" ~ Actual_shell_pre_mg * 0.3741429,
      Phase_1_treat == "Hyp"  ~ Actual_shell_pre_mg * 0.4050769,
      Phase_1_treat == "Both" ~ Actual_shell_pre_mg * 0.3490909,
      TRUE ~ NA_real_)
  )

View(N_C_byweight)

####
                                          