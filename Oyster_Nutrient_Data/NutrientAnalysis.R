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
#read in data
Tissue_Nutrient_df <- read_csv("Phase1_tissuenutrient_working.csv", 
                        col_types = cols(wt_percent_N = col_number(), 
                                         wt_percent_C = col_number()))
Tissue_Nutrient_df
#filter out data points (errors/replicates)
N_T_df <- Tissue_Nutrient_df%>%
  filter(N_exclude != "Y" | is.na(N_exclude))

C_T_df <- Tissue_Nutrient_df%>%
  filter(C_exclude != "Y" | is.na(C_exclude))
View(C_T_df)

#plot with mean and SD, TISSUE, NITROGEN
summary_stats <- N_T_df %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_N, na.rm = TRUE),
    se_growth = std.error(wt_percent_N, na.rm = TRUE))
View(summary_stats)

#reorder Phase_1_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "% N in Tissue") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  ylim(6,9)+
  theme(legend.position = "none") # Remove legend


#plot with mean and SD, TISSUE, CARBON
summary_stats <- C_T_df %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_C, na.rm = TRUE),
    se_growth = std.error(wt_percent_C, na.rm = TRUE))
View(summary_stats)

#reorder Phase_1_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "% C in Tissue") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  ylim(35,45)+
  theme(legend.position = "none") # Remove legend





#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 

N_T_df
## % N
#nothing significant with or without log transformation
Nm1 <- lmer(log(wt_percent_N) ~ Phase_1_DO*Phase_1_temp + 
              (1|Phase_1_rep_R), data = N_T_df, REML=TRUE)
Anova(Nm1, test="F", type="III")
N_T_df
#post hocs
#emmeans(Nm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#diagnostics
leveneTest(wt_percent_N~Phase_1_DO*Phase_1_temp, N_T_df)
m1.e <- residuals(Nm1) 
qqnorm(m1.e)
qqline(m1.e)



## % C
#nothing significant
Cm1 <- lmer(log(wt_percent_C) ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = C_T_df, REML=TRUE)
Anova(Cm1, test="F", type="III")

#post hocs
#emmeans(Cm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#diagnostics
leveneTest(log(wt_percent_C)~Phase_1_DO*Phase_1_temp, C_T_df)
m1.e <- residuals(Cm1) 
qqnorm(m1.e)
qqline(m1.e)




#Calculate summary stats for Nitrogen to multiply with TISSUE mass
summary_stats_TN <- N_T_df %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_N, na.rm = TRUE),
    se_growth = std.error(wt_percent_N, na.rm = TRUE))
View(summary_stats_TN)

#Now C
summary_stats_TC <- C_T_df %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_C, na.rm = TRUE),
    se_growth = std.error(wt_percent_C, na.rm = TRUE))
View(summary_stats_TC)
View(C_T_df)



##multiply percentages by weights for PHASE 1
#read in the data
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

#get mg nutrients in TISSUE
N_C_bytissueweight <- Nutrientsbyweight %>%
  filter(Exclude_pre_analysis != "Y" | is.na(Exclude_pre_analysis)) %>%
  mutate(
    Phase1N_tissue = case_when(
      Phase_1_treat == "Cont" ~ ((Actual_tissue_pre_mg * 7.741667)/100),
      Phase_1_treat == "Warm" ~ ((Actual_tissue_pre_mg * 7.641667)/100),
      Phase_1_treat == "Hyp"  ~ ((Actual_tissue_pre_mg * 7.890909)/100),
      Phase_1_treat == "Both" ~ ((Actual_tissue_pre_mg * 7.077778)/100),
      TRUE ~ NA_real_
    ))%>%
    mutate(
      Phase1C_tissue = case_when(
        Phase_1_treat == "Cont" ~ ((Actual_tissue_pre_mg * 38.53333)/100),
        Phase_1_treat == "Warm" ~ ((Actual_tissue_pre_mg * 38.47500)/100),
        Phase_1_treat == "Hyp"  ~ ((Actual_tissue_pre_mg * 39.71667)/100),
        Phase_1_treat == "Both" ~ ((Actual_tissue_pre_mg * 39.78182)/100),
        TRUE ~ NA_real_)
  )

View(N_C_bytissueweight)



#read in nutrient data
N_C_bytissueweight <- read_csv("Phase1_nutrient_merged.csv")

## mg N in tissue
Nm1 <- lmer(log(Phase1N_tissue) ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = N_C_bytissueweight, REML=TRUE)
Anova(Nm1, test="F", type="III")

#post hocs
emmeans(Nm1,specs = pairwise ~ Phase_1_DO, adjust = "none")
emmeans(Nm1,specs = pairwise ~ Phase_1_temp, adjust = "none")

exp(2.66) #hyp
exp(2.81) #norm
(16.60992-14.29629)/16.60992 #14 percent less nitrogen stored in oyster tissue under hypoxia

#diagnostics
leveneTest(log(Phase1N_tissue)~Phase_1_DO*Phase_1_temp, N_C_bytissueweight)
m1.e <- residuals(Nm1) 
qqnorm(m1.e)
qqline(m1.e)


#plot with mean and SD, TISSUE NITROGEN
summary_stats <- N_C_bytissueweight %>%
  group_by(Phase_1_DO) %>%
  summarize(
    mean_growth = mean(Phase1N_tissue, na.rm = TRUE),
    se_growth = std.error(Phase1N_tissue, na.rm = TRUE))

View(summary_stats)

#reorder Phase_1_treat
summary_stats$Phase_1_DO <- factor(summary_stats$Phase_1_DO, 
                                      levels = c("Norm", "Hyp"))

ggplot(summary_stats, aes(x = Phase_1_DO, y = mean_growth, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.08, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 Treatment", y = "Nitrogen in Tissue (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  ylim(16,20)+
  theme(legend.position = "none") # Remove legend


## % C by weight
Cm1 <- lmer(log(Phase1C_tissue) ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = N_C_bytissueweight, REML=TRUE)
Anova(Cm1, test="F", type="III")

#post hocs
emmeans(Cm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#diagnostics
leveneTest(log(Phase1C_tissue)~Phase_1_DO*Phase_1_temp, N_C_bytissueweight)
m1.e <- residuals(Cm1) 
qqnorm(m1.e)
qqline(m1.e)

#plot with mean and SD, HYPOXIA ON TISSUE NITROGEN
summary_stats <- N_C_bytissueweight %>%
  group_by(Phase_1_DO) %>%
  summarize(
    mean_growth = mean(Phase1C_tissue, na.rm = TRUE),
    se_growth = std.error(Phase1C_tissue, na.rm = TRUE))

View(summary_stats)

#reorder Phase_1_treat
summary_stats$Phase_1_DO <- factor(summary_stats$Phase_1_DO, 
                                   levels = c("Norm", "Hyp"))

ggplot(summary_stats, aes(x = Phase_1_DO, y = mean_growth, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 Treatment", y = "Carbon in Tissue (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  theme(legend.position = "none") # Remove legend

#plot with mean and SD, WARMING ON TISSUE NITROGEN
summary_stats <- N_C_bytissueweight %>%
  group_by(Phase_1_temp) %>%
  summarize(
    mean_growth = mean(Phase1C_tissue, na.rm = TRUE),
    se_growth = std.error(Phase1C_tissue, na.rm = TRUE))

View(summary_stats)

#reorder Phase_1_treat
summary_stats$Phase_1_temp <- factor(summary_stats$Phase_1_temp, 
                                   levels = c("Ambient", "Warm"))

ggplot(summary_stats, aes(x = Phase_1_temp, y = mean_growth, color = Phase_1_temp)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Ambient" = "darkblue", "Warm" = "#B00149")) +
  labs(x = "Phase 1 Treatment", y = "Carbon in Tissue (%mg)") +
  theme(legend.position = "none") # Remove legend


#plot with mean and SD, CARBON
summary_stats <- N_C_bytissueweight %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_growth = mean(Phase2C_tissue, na.rm = TRUE),
    se_growth = std.error(Phase2C_tissue, na.rm = TRUE))
View(summary_stats)


#reorder Phase_1_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Carbon in Tissue (%mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend





####Shell N and C ####
Shell_Nutrient_df <- read_csv("Phase1_shellnutrient_working.csv", 
                        col_types = cols(wt_percent_N = col_number(), 
                                         wt_percent_C = col_number()))

N_S_df <- Shell_Nutrient_df%>%
  filter(N_exclude != "Y" | is.na(N_exclude))

C_S_df <- Shell_Nutrient_df%>%
  filter(C_exclude != "Y" | is.na(C_exclude))


#plot with mean and SD, shell, NITROGEN
summary_stats <- N_S_df %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_N, na.rm = TRUE),
    se_growth = std.error(wt_percent_N, na.rm = TRUE))

View(summary_stats)

#reorder Phase_1_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "% N in Shell") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  ylim(0.25, 0.32)+
  theme(legend.position = "none") # Remove legend


#plot with mean and SD, TISSUE, CARBON
summary_stats <- C_S_df %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_C, na.rm = TRUE),
    se_growth = std.error(wt_percent_C, na.rm = TRUE))
View(summary_stats)

#reorder Phase_1_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "% C in Shell") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  ylim(15, 18) +
  theme(legend.position = "none") # Remove legend

## % N analysis
Nm1 <- lmer(wt_percent_N ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = N_S_df, REML=TRUE)
Anova(Nm1, test="F", type="III")

#post hocs
emmeans(Cm1,specs = pairwise ~ Phase_1_DO, adjust = "none")
emmeans(Cm1,specs = pairwise ~ Phase_1_temp, adjust = "none")

#diagnostics
leveneTest(wt_percent_N~Phase_1_DO*Phase_1_temp, N_S_df)
m1.e <- residuals(Nm1) 
qqnorm(m1.e)
qqline(m1.e)

## % C analysis
Cm1 <- lmer(wt_percent_C ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = C_S_df, REML=TRUE)
Anova(Cm1, test="F", type="III")

#post hocs
emmeans(Cm1,specs = pairwise ~ Phase_1_DO, adjust = "none")
emmeans(Cm1,specs = pairwise ~ Phase_1_temp, adjust = "none")

#diagnostics
leveneTest(wt_percent_C~Phase_1_treat, C_S_df)
m1.e <- residuals(Cm1) 
qqnorm(m1.e)
qqline(m1.e)



##nutrients by weight shell
##using this as a df to add shell N and C on to tissue N and C
N_C_bytissueweight

#nitrogen summmary stats
summary_stats_SN <- N_S_df %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_N, na.rm = TRUE),
    se_growth = std.error(wt_percent_N, na.rm = TRUE))
View(summary_stats_SN)

#carbon summary stats
summary_stats_SC <- C_S_df %>%
  group_by(Phase_1_treat) %>%
  summarize(
    mean_growth = mean(wt_percent_C, na.rm = TRUE),
    se_growth = std.error(wt_percent_C, na.rm = TRUE))
View(summary_stats_SC)

#update equations with summary stats

N_C_byshellweight <- N_C_bytissueweight %>%
  filter(Exclude_pre_analysis != "Y" | is.na(Exclude_pre_analysis)) %>%
  mutate(
    Phase1N_shell = case_when(
      Phase_1_treat == "Cont" ~ ((Actual_shell_pre_mg * 0.2833333)/100),
      Phase_1_treat == "Warm" ~ ((Actual_shell_pre_mg * 0.2916667)/100),
      Phase_1_treat == "Hyp"  ~ ((Actual_shell_pre_mg * 0.2916667)/100),
      Phase_1_treat == "Both" ~ ((Actual_shell_pre_mg * 0.3000000)/100),
      TRUE ~ NA_real_
    ))%>%
  mutate(
    Phase1C_shell = case_when(
      Phase_1_treat == "Cont" ~ ((Actual_shell_pre_mg * 16.20000)/100),
      Phase_1_treat == "Warm" ~ ((Actual_shell_pre_mg * 16.23333)/100),
      Phase_1_treat == "Hyp"  ~ ((Actual_shell_pre_mg * 16.24167)/100),
      Phase_1_treat == "Both" ~ ((Actual_shell_pre_mg * 16.71667)/100),
      TRUE ~ NA_real_)
  )

View(N_C_byshellweight)

as.numeric(N_C_byshellweight$Phase1N_tissue);as.numeric(N_C_byshellweight$Phase1C_tissue);as.numeric(N_C_byshellweight$Phase1N_shell);as.numeric(N_C_byshellweight$Phase1C_shell)

#write CSV with nitrogen and carbon in mg for each oyster
#write_csv(N_C_byshellweight, "~/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Nutrient_Data/Phase1_nutrient_merged.csv")


#read in nutrient data
N_C_byshellweight <- read_csv("Phase1_nutrient_merged.csv")

## % N by weight in shell
Nm1 <- lmer(Phase1N_shell ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = N_C_byshellweight, REML=TRUE)

Anova(Nm1, test="F", type="III")

#post hocs
emmeans(Nm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#calculate
(1.15-1.04)/1.15

#diagnostics
leveneTest(Phase1N_shell~Phase_1_DO*Phase_1_temp, N_C_byshellweight)
m1.e <- residuals(Nm1) 
qqnorm(m1.e)
qqline(m1.e)


#plot with mean and SD, HYPOXIA ON SHELL NITROGEN
summary_stats <- N_C_byshellweight %>%
  group_by(Phase_1_DO) %>%
  summarize(
    mean_growth = mean(Phase1N_shell, na.rm = TRUE),
    se_growth = std.error(Phase1N_shell, na.rm = TRUE))

View(summary_stats)

#reorder Phase_1_treat
summary_stats$Phase_1_DO <- factor(summary_stats$Phase_1_DO, 
                                   levels = c("Norm", "Hyp"))

ggplot(summary_stats, aes(x = Phase_1_DO, y = mean_growth, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.08, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 Treatment", y = "Nitrogen in Shell (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  ylim(1, 1.2)+
  theme(legend.position = "none") # Remove legend



## % C by weight in shell
Cm1 <- lmer(Phase1C_shell ~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = N_C_byshellweight, REML=TRUE)
Anova(Cm1, test="F", type="III")

#post hocs
emmeans(Cm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#calculate
(64.8-57.9)/64.8
#diagnostics
leveneTest(Phase1C_shell~Phase_1_treat, N_C_byshellweight)
m1.e <- residuals(Cm1) 
qqnorm(m1.e)
qqline(m1.e)

#plot with mean and SD, HYPOXIA ON SHELL CARBON
summary_stats <- N_C_byshellweight %>%
  group_by(Phase_1_DO) %>%
  summarize(
    mean_growth = mean(Phase1C_shell, na.rm = TRUE),
    se_growth = std.error(Phase1C_shell, na.rm = TRUE))

View(summary_stats)

#reorder Phase_1_treat
summary_stats$Phase_1_DO <- factor(summary_stats$Phase_1_DO, 
                                   levels = c("Norm", "Hyp"))

ggplot(summary_stats, aes(x = Phase_1_DO, y = mean_growth, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.08, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  theme(panel.background = element_rect(fill = "#E5E5E5")) + #fill background light grey
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 Treatment", y = "Carbon in Shell (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  theme(legend.position = "none") # Remove legend
