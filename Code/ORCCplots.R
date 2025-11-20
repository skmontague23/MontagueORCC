
#create mean and standard error plots

#load packages
library(ggplot2)
library(dplyr)
library(plotrix) #for standard error

#### TISSUE GROWTH PLOTS ####
#All possible treatment combinations

#effect of phase 1, TISSUE
pre_summary_stats_t <- Growth_Data_forR_pre %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_pre = mean(Actual_tissue_pre_mg, na.rm = TRUE),
    se_pre = std.error(Actual_tissue_pre_mg, na.rm = TRUE))

  #reorder effect of phase 1 Phase_1_treat
  pre_summary_stats_t$Phase_1_treat <- factor(pre_summary_stats_t$Phase_1_treat, 
                                            levels = c("Cont", "Warm","Hyp", "Both"))
  #plot TISSUE
  ggplot(pre_summary_stats_t, aes(x = Phase_1_treat, y = mean_pre, color = Phase_1_treat)) +
    geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
    geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                  width = 0.1, position = position_dodge(0.9)) + # Error bars for SD
    theme_classic(base_size = 18)+
    scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
    labs(x = "Phase 1 Treatment", y = "Tissue Mass (mg)") +
    scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
    theme(legend.position = "none") # Remove legend

  
#all, TISSUE
summary_stats_t <- Growth_Data_forR_full %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Actual_tissue_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_mg, na.rm = TRUE))

  #reorder Phase_1_treat and Phase_2_treat
  summary_stats_t$Phase_1_treat <- factor(summary_stats_t$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))
  summary_stats_t$Phase_2_treat <- factor(summary_stats_t$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))
  
  #plot with mean and SD, TISSUE
  ggplot(summary_stats_t, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
    geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
    geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                  width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
    theme_classic(base_size = 20) +
    guides(color = "none") + # Remove legend for color
    facet_wrap(vars(Phase_2_treat), 
               labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")), 
               scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
    scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
    labs(x = "Phase 1 Treatment", y = "Mean Tissue Growth (mg)") +
    scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
    theme(legend.position = "none") # Remove legend

  
#standardize initial size, TISSUE
filt_summary_stats_t <- Growth_Data_filtered %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Actual_tissue_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_mg, na.rm = TRUE))

  #reorder Phase_1_treat and Phase_2_treat for standardize initial size
  filt_summary_stats_t$Phase_1_treat <- factor(filt_summary_stats_t$Phase_1_treat, 
                                             levels = c("Cont", "Warm","Hyp", "Both"))
  filt_summary_stats_t$Phase_2_treat <- factor(filt_summary_stats_t$Phase_2_treat, 
                                             levels = c("Cont", "Warm","Hyp", "Both"))
  #plot, TISSUE
  ggplot(filt_summary_stats_t, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
    geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
    geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                  width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
    theme_classic(base_size = 20) +
    guides(color = "none") + # Remove legend for color
    facet_wrap(vars(Phase_2_treat), 
               labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")), 
               scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
    scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
    labs(x = "Phase 1 Treatment", y = "Mean Tissue Growth (mg)") +
    scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
    theme(legend.position = "none") # Remove legend

  
  

#### SHELL GROWTH ####
#All possible treatment combinations

#effect of phase 1 data, SHELL
pre_summary_stats_s <- Growth_Data_forR_pre %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_pre = mean(Actual_shell_pre_mg, na.rm = TRUE),
    se_pre = std.error(Actual_shell_pre_mg, na.rm = TRUE))

  #reorder effect of phase 1 Phase_1_treat
  pre_summary_stats_s$Phase_1_treat <- factor(pre_summary_stats_s$Phase_1_treat, 
                                            levels = c("Cont", "Warm","Hyp",  "Both"))
  
  #plot SHELL
  ggplot(pre_summary_stats_s, aes(x = Phase_1_treat, y = mean_pre, color = Phase_1_treat)) +
    geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
    geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                  width = 0.1, position = position_dodge(0.9)) + # Error bars for SD
    theme_classic(base_size = 18)+
    scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
    labs(x = "Phase 1 Treatment", y = "Shell Mass (mg)") +
    scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
    ylim(330, 450)+
    theme(legend.position = "none") # Remove legend

#all, SHELL
summary_stats_s <- Growth_Data_forR_full%>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Actual_shell_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_shell_growth_mg, na.rm = TRUE))

  #reorder Phase_1_treat and Phase_2_treat
  summary_stats_s$Phase_1_treat <- factor(summary_stats_s$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp",  "Both"))
  summary_stats_s$Phase_2_treat <- factor(summary_stats_s$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp",  "Both"))
  
  #plot with mean and SD, SHELL
  ggplot(summary_stats_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
    geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
    geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                  width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
    theme_classic(base_size = 18) +
    guides(color = "none") + # Remove legend for color
    facet_wrap(vars(Phase_2_treat), labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")),
               scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
    scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
    labs(x = "Phase 1 Treatment", y = "Mean Shell Growth (mg)") +
    scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
    theme(legend.position = "none") # Remove legend

#standardize initial size, SHELL
filt_summary_stats_s <- Growth_Data_filtered %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Actual_shell_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_shell_growth_mg, na.rm = TRUE))

  #reorder Phase_1_treat and Phase_2_treat for standardize initial size
  filt_summary_stats_s$Phase_1_treat <- factor(filt_summary_stats_s$Phase_1_treat, 
                                             levels = c("Cont", "Warm","Hyp",  "Both"))
  filt_summary_stats_s$Phase_2_treat <- factor(filt_summary_stats_s$Phase_2_treat, 
                                             levels = c("Cont", "Warm","Hyp",  "Both"))


  

#### TISSUE:SHELL PLOTS ####
#All possible treatment combinations

#effect of phase 1 data, T:S
pre_summary_stats_t_s <- Growth_Data_forR_pre %>%
  group_by(Phase_1_temp) %>%
  mutate(
    mean_pre = mean(Ratio_tissue_shell_pre_mg, na.rm = TRUE),
    se_pre = std.error(Ratio_tissue_shell_pre_mg, na.rm = TRUE))

  #reorder effect of phase 1 Phase_1_treat
  pre_summary_stats_t_s$Phase_1_treat <- factor(pre_summary_stats_t_s$Phase_1_treat, 
                                              levels = c("Cont", "Warm","Hyp",  "Both"))
  
  #plot with mean and SD, T:S
  ggplot(pre_summary_stats_t_s, aes(x = Phase_1_treat, y = mean_pre, color = Phase_1_treat)) +
    geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
    geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                  width = 0.1, position = position_dodge(0.9)) + # Error bars for SD
    theme_classic(base_size = 18) +
    guides(color = "none") + # Remove legend for color
    scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
    labs(x = "Phase 1 Treatment", y = "Tissue:Shell Mass (mg)") +
    scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
    ylim(0.58, 0.66)+
    theme(legend.position = "none") # Remove legend

  
#all T:S
summary_stats_t_s <-Growth_Data_forR_full %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Ratio_tissue_shell_mg, na.rm = TRUE),
    se_growth = std.error(Ratio_tissue_shell_mg, na.rm = TRUE))

  #reorder Phase_1_treat and Phase_2_treat
  summary_stats_t_s$Phase_1_treat <- factor(summary_stats_t_s$Phase_1_treat, 
                                          levels = c("Cont", "Warm","Hyp",  "Both"))
  summary_stats_t_s$Phase_2_treat <- factor(summary_stats_t_s$Phase_2_treat, 
                                          levels = c("Cont", "Warm","Hyp",  "Both"))
  
  #plot with mean and SD, T:S
  ggplot(summary_stats_t_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
    geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
    geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                  width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
    theme_classic(base_size = 18) +
    guides(color = "none") + # Remove legend for color
    facet_wrap(vars(Phase_2_treat), 
               labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")), 
               scales = "fixed", nrow = 1) +
    scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
    labs(title = "Phase 2 Treatment", x = "Phase 1 Treatment", y = "Mean Tissue:Shell Growth (mg)") +
    scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
    ylim(0.35,0.55)+
    theme(
      plot.title = element_text(size = 18, hjust = 0.5)  # <-- title font size here
    )

  
#standardize initial size, T:S
filt_summary_stats_t_s <- Growth_Data_filtered %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Ratio_tissue_shell_mg, na.rm = TRUE),
    se_growth = std.error(Ratio_tissue_shell_mg, na.rm = TRUE))



#### MEAT YIELD PLOTS ####
##wet tissue weight / total weight = meat yield

#effects of phase 1, MEAT YIELD
head(Growth_Data_forR_pre_my$meat_yield)
View(Growth_Data_forR_pre_my)
head(pre_summary_stats_my)

pre_summary_stats_my <- Growth_Data_forR_pre_my %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_pre = mean(meat_yield, na.rm = TRUE),
    se_pre = std.error(meat_yield, na.rm = TRUE))
colnames(pre_summary_stats_my)

#order
pre_summary_stats_my$Phase_1_treat <- factor(pre_summary_stats_my$Phase_1_treat, 
                                             levels = c("Cont", "Warm","Hyp", "Both"))
#Meat yield plot with mean and SD
ggplot(pre_summary_stats_my, aes(x = Phase_1_treat, y = mean_pre, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.1, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Meat Yield") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend


#all, MEAT YEILD
colnames(Growth_Data_forR_full)

Growth_Data_forR_full%>%
  mutate(meat_yield = (Actual_tissue_growth_mg/(Actual_tissue_growth_mg+Actual_shell_growth_mg))) -> Growth_Data_forR_full_my

summary_stats_my <- Growth_Data_forR_full_my %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(meat_yield, na.rm = TRUE),
    se_growth = std.error(meat_yield, na.rm = TRUE))

summary_stats_my$Phase_2_treat <- factor(summary_stats_my$Phase_2_treat, 
                                         levels = c("Cont", "Warm","Hyp", "Both"))
summary_stats_my$Phase_1_treat <- factor(summary_stats_my$Phase_1_treat, 
                                         levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats_my, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), 
             labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")), 
             scales = "fixed", nrow = 1) +
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(title = "Phase 2 Treatment", x = "Phase 1 Treatment", y = "Meat Yield") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none", plot.title = element_text(size = 18, hjust = 0.5)) # Remove legend







#### SIGNIFICANT INTERACTIONS PLOTS

#TISSUE PLOTS, significant interactions

##effects of phase 1, PHASE 1 DO & TISSUE
pre_summary_stats_t <- Growth_Data_forR_pre %>%
  group_by(Phase_1_DO) %>%
  mutate(
    mean_pre = mean(Actual_tissue_pre_mg, na.rm = TRUE),
    se_pre = std.error(Actual_tissue_pre_mg, na.rm = TRUE))
#reorder
pre_summary_stats_t$Phase_1_DO <- factor(pre_summary_stats_t$Phase_1_DO, levels = c("Norm", "Hyp"))
#Tissue plot with mean and SD
ggplot(pre_summary_stats_t, aes(x = Phase_1_DO, y = mean_pre, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  labs(x = "Phase 1 DO", y = "Tissue Mass (mg)") +
  ylim(210,260)+
  theme(legend.position = "none") # Remove legend





#SHELL PLOTS, significant interactions

#shell growth and PHASE 2 DO
stats_s <- Growth_Data_forR_full %>%
  group_by(Phase_2.1_DO) %>%
  mutate(
    mean_growth = mean(Actual_shell_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_shell_growth_mg, na.rm = TRUE))
#reorder
stats_s$Phase_2.1_DO <- factor(stats_s$Phase_2.1_DO, 
                               levels = c("Norm","Hyp"))
#Shell plot with mean and SD
ggplot(stats_s, aes(x = Phase_2.1_DO, y = mean_growth, color = Phase_2.1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  labs(x = "Phase 2 DO", y = "Shell Growth (mg)") +
  ylim(200, 300)+
  theme(legend.position = "none") # Remove legend


#shell growth and Phase_1_DO:Phase_1_temp
stats_s <- Growth_Data_forR_full %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_growth = mean(Actual_shell_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_shell_growth_mg, na.rm = TRUE))
#reorder
stats_s$Phase_1_treat <- factor(stats_s$Phase_1_treat, 
                                levels = c("Cont", "Warm","Hyp", "Both"))
#Shell plot with mean and SD
ggplot(stats_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.1, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 16) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  labs(x = "Phase 1 Treatment", y = "Shell Growth in Phase 2 (mg)") +
  ylim(230,280) +
  theme(legend.position = "none") # Remove legend


#shell growth and Phase_1_temp:Phase_2.1_DO
stats_s <- Growth_Data_forR_full %>%
  group_by(Phase_1_temp, Phase_2.1_DO) %>%
  mutate(
    mean_growth = mean(Actual_shell_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_shell_growth_mg, na.rm = TRUE))

#reorder Phase_1_temp:Phase_2.1_DO
stats_s$Phase_1_temp <- factor(stats_s$Phase_1_temp, 
                                        levels = c("Ambient", "Warm"))
stats_s$Phase_2.1_DO <- factor(stats_s$Phase_2.1_DO, 
                                        levels = c("Norm","Hyp"))

#plot with mean and SD, SHELL
ggplot(stats_s, aes(x = Phase_1_temp, y = mean_growth, color = Phase_1_temp)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.1, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2.1_DO), labeller = as_labeller(c("Norm" = "Normoxic", "Hyp" = "Hypoxic")),
             scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_manual(values = c("Warm" = "#B00149", "Ambient" = "darkblue")) +
  labs(x = "Phase 1 Temperature", y = "Mean Shell Growth (mg)") +
  scale_x_discrete(labels = c("Norm" = "Normoxic", "Hyp" = "Hypoxic")) +
  theme(legend.position = "none") # Remove legend





#TISSUE:SHELL PLOTS, significant interactions

##effects of phase 1, PHASE 1 DO & TISSUE:SHELL
pre_summary_stats_t_s <- Growth_Data_forR_pre %>%
  group_by(Phase_1_DO) %>%
  mutate(
    mean_pre = mean(Ratio_tissue_shell_pre_mg, na.rm = TRUE),
    se_pre = std.error(Ratio_tissue_shell_pre_mg, na.rm = TRUE))

pre_summary_stats_t_s$Phase_1_DO <- factor(pre_summary_stats_t_s$Phase_1_DO, 
                                         levels = c("Norm", "Hyp"))
#Tissue:Shell plot with mean and SD
ggplot(pre_summary_stats_t_s, aes(x = Phase_1_DO, y = mean_pre, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 DO", y = "Tissue:Shell Mass (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  ylim(0.61, 0.65)+
  theme(legend.position = "none") # Remove legend


##effects of phase 1, PHASE 1 TEMP & TISSUE:SHELL
pre_summary_stats_t_s <- Growth_Data_forR_pre %>%
  group_by(Phase_1_temp) %>%
  mutate(
    mean_pre = mean(Ratio_tissue_shell_pre_mg, na.rm = TRUE),
    se_pre = std.error(Ratio_tissue_shell_pre_mg, na.rm = TRUE))

View(pre_summary_stats_t_s)

pre_summary_stats_t_s$Phase_1_temp <- factor(pre_summary_stats_t_s$Phase_1_temp, 
                                           levels = c("Ambient", "Warm"))
#Tissue:Shell plot with mean and SD
ggplot(pre_summary_stats_t_s, aes(x = Phase_1_temp, y = mean_pre, color = Phase_1_temp)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Warm" = "#B00149", "Ambient" = "darkblue")) +
  labs(x = "Phase 1 Temp", y = "Tissue:Shell Mass (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  ylim(0.61, 0.65)+
  theme(legend.position = "none") # Remove legend








##STandardized initial size
#Tissue plot with mean and SD
ggplot(filt_summary_stats_t, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic() +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_brewer(palette = "Set2") + # Use color palette for points
  labs(x = "Phase 1 Treatment", y = "Mean Tissue Growth (mg)") +
  theme(legend.position = "none") # Remove legend

#Shell plot with mean and SD
ggplot(filt_summary_stats_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic() +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_brewer(palette = "Set2") + # Use color palette for points
  labs(x = "Phase 1 Treatment", y = "Mean Shell Growth (mg)") +
  theme(legend.position = "none") # Remove legend

#Normoxic vs hypoxic phase 2
colnames(filt_summary_stats_s)
ggplot(filt_summary_stats_s, aes(x = Phase_2.1_DO, y = mean_pre, color = Phase_2.1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  labs(x = "Phase 1 DO", y = "Tissue Mass (mg)") +
  theme(legend.position = "none") # Remove legend

#Tissue:Shell plot with mean and SD
ggplot(filt_summary_stats_t_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic() +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_brewer(palette = "Set2") + # Use color palette for points
  labs(x = "Phase 1 Treatment", y = "Mean Tissue:Shell Growth (mg)") +
  theme(legend.position = "none") # Remove legend



#error bars are crazy on T:S plot, visualize data
#histogram
ggplot(Growth_Data_forR_clean) +
  aes(x = Ratio_tissue_shell_mg) +
  geom_histogram(bins = 30L, fill = "#112446") +
  theme_classic()
#scatterplot
ggplot(Growth_Data_forR_clean) +
  aes(x = Actual_shell_growth_mg, y = Actual_tissue_growth_mg) +
  geom_point(colour = "#112446") +
  theme_classic()



ggplot(Growth_Data_forR_clean) +
  aes(x = Phase_1_treat, y = Actual_tissue_growth_mg, fill = Phase_1_treat) + # Map fill to Phase_1_treat
  geom_boxplot() +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(vars(Phase_2_treat), scales = "free", nrow = 1) + # Arrange facets side by side
  scale_fill_brewer(palette = "Set2") # Use a color palette for distinct colors

ggplot(Growth_Data_forR_clean) +
  aes(x = Phase_1_treat, y = Actual_shell_growth_mg, fill = Phase_1_treat) + # Map fill to Phase_1_treat
  geom_boxplot() +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(vars(Phase_2_treat), scales = "free", nrow = 1) + # Arrange facets side by side
  scale_fill_brewer(palette = "Set2") # Use a color palette for distinct colors

ggplot(Growth_Data_forR_clean) +
  aes(x = Phase_1_treat, y = Ratio_tissue_shell_mg, fill = Phase_1_treat) + # Map fill to Phase_1_treat
  geom_boxplot() +
  theme_classic() +
  theme(legend.position="none") +
  facet_wrap(vars(Phase_2_treat), scales = "free", nrow = 1) + # Arrange facets side by side
  scale_fill_brewer(palette = "Set2") # Use a color palette for distinct colors



# Convert emmeans and contrasts data to data frames
emmeans_df <- as.data.frame(emm2$emmeans)
contrasts_df <- as.data.frame(emm2$contrasts)

# Basic plot for emmeans
ggplot(emmeans_df, aes(x = interaction(Phase_1_DO, Phase_2.1_DO), 
                       y = emmean, fill = Phase_1_temp)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Phase_2.1_temp ~ Phase_1_temp) +
  labs(x = "DO Levels (Phase 1 and Phase 2.1)", 
       y = "Estimated Marginal Means", 
       title = "Effects of Phase-wise DO and Temperature Interactions") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Overlay contrasts as annotations (optional based on contrast significance)
# Example: Adding contrast annotations if significant contrasts found
contrasts_df <- contrasts_df %>% 
  mutate(significant = ifelse(p.value < 0.05, "Significant", "Not Significant"))

ggplot(contrasts_df, aes(x = interaction(Phase_1_DO, Phase_2.1_DO), 
                         y = emmean, fill = Phase_1_temp)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(Phase_2.1_temp ~ Phase_1_temp) +
  labs(x = "DO Levels (Phase 1 and Phase 2.1)", 
       y = "Estimated Marginal Means", 
       title = "Effects of Phase-wise DO and Temperature Interactions") +
  theme_minimal() +
  theme(legend.position = "bottom")

#looking at shell growth data by treatment
ggplot(Growth_Data_forR) +
  aes(x = Phase1_Phase2_treat, y = Actual_shell_growth_mg) +
  geom_point(colour = "#112446") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 27L, vjust = 0.6))

#looking at shell growth by initial shell
ggplot(Growth_Data_forR) +
  aes(x = Actual_shell_pre_mg, y = Actual_shell_growth_mg) +
  geom_point(shape = "bullet", colour = "#112446") +
  geom_smooth(method = "lm", colour = "red", se = FALSE) + 
  theme_classic()

#looking at distribution of tissue and shell growth proportions

ggplot(Growth_Data_forR_full) +
  aes(x = prop_shell_growth, y = prop_tissue_growth) +
  geom_point(colour = "#112446") +
  theme_classic()

ggplot(Growth_Data_forR_full) +
  aes(x = prop_tissue_growth) +
  geom_histogram(bins = 30L, fill = "#112446") +
  theme_classic()

ggplot(Growth_Data_forR_full) +
  aes(x = prop_shell_growth) +
  geom_histogram(bins = 30L, fill = "#112446") +
  theme_classic()

ggplot(Growth_Data_forR_full) +
  aes(x= Actual_shell_pre_mg, y= Actual_tissue_pre_mg, colour = Phase_1_treat) +
  geom_point() +
  theme_minimal()


colnames(Growth_Data_forR_full)

##for formatting tables into overleaf
install.packages("stargazer")
library(stargazer)
stargazer(m1, type = "latex")

#normal table
install.packages("flextable")
library(flextable)
flextable(Warm1)





##graphs for conference presentation

#tissue growth
stats_t <- Growth_Data_forR_full %>%
  group_by(Phase_2.1_DO) %>%
  mutate(
    mean_growth = mean(Actual_tissue_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_mg, na.rm = TRUE))
#reorder
stats_t$Phase_2.1_DO <- factor(stats_t$Phase_2.1_DO, 
                               levels = c("Norm","Hyp"))
#Tissue plot with mean and SD
ggplot(stats_t, aes(x = Phase_2.1_DO, y = mean_growth, color = Phase_2.1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  labs(x = "Phase 2 DO", y = "Tissue Growth (mg)") +
  ylim(90, 130)+
  theme(legend.position = "none") # Remove legend



#tissue:shell growth
stats_ts <- Growth_Data_forR_full %>%
  group_by(Phase_1_temp, Phase_2.1_DO) %>%
  mutate(
    mean_growth = mean(Ratio_tissue_shell_mg, na.rm = TRUE),
    se_growth = std.error(Ratio_tissue_shell_mg, na.rm = TRUE))
#order
stats_ts$Phase_2.1_DO <- factor(stats_ts$Phase_2.1_DO, 
                                levels = c("Norm", "Hyp"))
stats_ts$Phase_1_temp <- factor(stats_ts$Phase_1_temp, 
                                levels = c("Ambient", "Warm"))
#Tissue:Shell plot with mean and SD
ggplot(stats_ts, aes(x = Phase_1_temp, y = mean_growth, color = Phase_1_temp)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  facet_wrap(~ Phase_2.1_DO, labeller = as_labeller(c("Hyp" = "Hypoxic", "Norm" = "Normoxic"))) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Warm" = "#B00149", "Ambient" = "darkblue")) +
  labs(title = "Phase 2 DO", x = "Phase 1 Temperature", y = "Mean Tissue:Shell Growth (mg)") +
  theme(legend.position = "none",
        plot.title = element_text(size = 18, hjust = 0.5))

#meat yield
stats_ts_my <- Growth_Data_forR_full_my %>%
  group_by(Phase_1_temp, Phase_2.1_DO) %>%
  mutate(
    mean_growth = mean(meat_yield, na.rm = TRUE),
    se_growth = std.error(meat_yield, na.rm = TRUE))
#order
stats_ts$Phase_2.1_DO <- factor(stats_ts$Phase_2.1_DO, 
                                levels = c("Norm", "Hyp"))
stats_ts$Phase_1_temp <- factor(stats_ts$Phase_1_temp, 
                                levels = c("Ambient", "Warm"))
#Tissue:Shell plot with mean and SD
ggplot(stats_ts, aes(x = Phase_1_temp, y = mean_growth, color = Phase_1_temp)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  facet_wrap(~ Phase_2.1_DO, labeller = as_labeller(c("Hyp" = "Hypoxic", "Norm" = "Normoxic"))) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Warm" = "#B00149", "Ambient" = "darkblue")) +
  labs(title = "Phase 2 DO", x = "Phase 1 Temperature", y = "Mean Tissue:Shell Growth (mg)") +
  theme(legend.position = "none",
        plot.title = element_text(size = 18, hjust = 0.5))




##STANDARDIZE SIZE AT START
#Tissue plot with mean and SD
ggplot(filt_summary_stats_t, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), 
             labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")), 
             scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Mean Tissue Growth (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend

#Shell plot with mean and SD
ggplot(filt_summary_stats_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")),
             scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Mean Shell Growth (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend

#Tissue:Shell plot with mean and SD
ggplot(filt_summary_stats_t_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), 
             labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")), 
             scales = "fixed", nrow = 1) +
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(title = "Phase 2 Treatment", x = "Phase 1 Treatment", y = "Mean Tissue:Shell Growth (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none", plot.title = element_text(size = 18, hjust = 0.5)) # Remove legend

#smaller graphs for subset to standardize initial size
Growth_Data_filtered <- Growth_Data_forR_full %>%
  filter(Phase_1_rep_R != "Cont01") %>%
  filter(Phase_1_rep_R != "Both02") %>%
  filter(Phase_1_rep_R != "Hyp06") %>%
  filter(Phase_1_rep_R != "Warm02") #works

#standardize initial size TISSUE, P2 DO
filt_summary_stats_t <- Growth_Data_filtered %>%
  group_by(Phase_2.1_DO) %>%
  mutate(
    mean_growth = mean(Actual_tissue_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_mg, na.rm = TRUE))
#order
filt_summary_stats_t$Phase_2.1_DO <- factor(filt_summary_stats_t$Phase_2.1_DO, 
                                            levels = c("Norm", "Hyp"))
#tissue phase 2 DO
ggplot(filt_summary_stats_t, aes(x = Phase_2.1_DO, y = mean_growth, color = Phase_2.1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  labs(x = "Phase 2 DO", y = "Mean Tissue Growth (mg)") +
  ylim(90,130)+
  theme(legend.position = "none") # Remove legend


#standardize initial size TISSUE, P2 TEMP & DO
filt_summary_stats_t <- Growth_Data_filtered %>%
  group_by(Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Actual_tissue_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_mg, na.rm = TRUE))

# Reorder Phase_1_treat and Phase_2_treat for standardize initial size
filt_summary_stats_t$Phase_1_treat <- factor(filt_summary_stats_t$Phase_1_treat, 
                                             levels = c("Cont", "Warm","Hyp", "Both"))
filt_summary_stats_t$Phase_2_treat <- factor(filt_summary_stats_t$Phase_2_treat, 
                                             levels = c("Cont", "Warm","Hyp", "Both"))
filt_summary_stats_s$Phase_1_treat <- factor(filt_summary_stats_s$Phase_1_treat, 
                                             levels = c("Cont", "Warm","Hyp",  "Both"))
filt_summary_stats_s$Phase_2_treat <- factor(filt_summary_stats_s$Phase_2_treat, 
                                             levels = c("Cont", "Warm","Hyp",  "Both"))
filt_summary_stats_t_s$Phase_1_treat <- factor(filt_summary_stats_t_s$Phase_1_treat, 
                                               levels = c("Cont", "Warm","Hyp",  "Both"))
filt_summary_stats_t_s$Phase_2_treat <- factor(filt_summary_stats_t_s$Phase_2_treat, 
                                               levels = c("Cont", "Warm","Hyp",  "Both"))
#tissue phase 2 DO
ggplot(filt_summary_stats_t, aes(x = Phase_2_treat, y = mean_growth, color = Phase_2_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.1, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control")) +
  labs(x = "Phase 2 Treatment", y = "Mean Tissue Growth (mg)") +
  ylim(80, 140)+
  theme(legend.position = "none") # Remove legend

#standardize initial size TISSUE, P1 TEMP & DO AND P2 TEMP
filt_summary_stats_t <- Growth_Data_filtered %>%
  group_by(Phase_1_treat, Phase_2.1_temp) %>%
  mutate(
    mean_growth = mean(Actual_tissue_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_mg, na.rm = TRUE))
#order
filt_summary_stats_t$Phase_2.1_temp <- factor(filt_summary_stats_t$Phase_2.1_temp, 
                                              levels = c("Ambient", "Warm"))
filt_summary_stats_t$Phase_1_treat <- factor(filt_summary_stats_t$Phase_1_treat, 
                                             levels = c("Cont", "Warm","Hyp", "Both"))
#t
ggplot(filt_summary_stats_t, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control")) +
  facet_wrap(vars(Phase_2.1_temp), 
             scales = "fixed", nrow = 1) +
  labs(title = "Phase 2 Temperature", x = "Phase 1 Treatment", y = "Mean Tissue Growth (mg)") +
  theme(legend.position = "none", plot.title = element_text(size = 18, hjust = 0.5)) # Remove legend

#standardize initial size SHELL, P2 DO
filt_summary_stats_s <- Growth_Data_filtered %>%
  group_by(Phase_2.1_DO) %>%
  mutate(
    mean_growth = mean(Actual_shell_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_shell_growth_mg, na.rm = TRUE))
#order
filt_summary_stats_s$Phase_2.1_DO <- factor(filt_summary_stats_s$Phase_2.1_DO, 
                                            levels = c("Norm", "Hyp"))
#shell phase 2 DO
ggplot(filt_summary_stats_s, aes(x = Phase_2.1_DO, y = mean_growth, color = Phase_2.1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  labs(x = "Phase 2 DO", y = "Mean Shell Growth (mg)") +
  ylim(200,300)+
  theme(legend.position = "none") # Remove legend

View(filt_summary_stats_s)




#MAP of CBL and Chesapeake
install.packages(c("ggplot2", "sf", "ggspatial", "rnaturalearth", "rnaturalearthdata"))
install.packages("devtools")
library(ggplot2); library(sf); library(rnaturalearth); library(ggspatial); 

devtools::install_github("ropensci/rnaturalearthhires")
states <- ne_states(country = "United States of America", returnclass = "sf")

# CBL coordinates
labs <- data.frame(
  lon = c(-76.4547, -76.1238, -76.610031),
  lat = c(38.3189, 38.5926, 38.434068),
  label = c("CBL", "HPL", "Patuxent"))

bbox <- c(xmin = -77.5, xmax = -75.5, ymin = 37.5, ymax = 39.5) #background box

# Plot
ggplot(states) +
  geom_rect(aes(xmin = bbox["xmin"], xmax = bbox["xmax"],
                ymin = bbox["ymin"], ymax = bbox["ymax"]),
            fill = "lightblue") +
  geom_sf(fill = "springgreen4", color = "black") +
  coord_sf(xlim = c(-77.5, -75.5), ylim = c(37.5, 39.5), expand = FALSE) +
  geom_point(data = labs, aes(lon, lat), color = "brown4", size = 1.5) +
  geom_text(data = labs, aes(lon, lat, label = label), nudge_y = -0.05, size = 3) +
  annotation_scale(location = "bl", width_hint = 0.3) +
  annotation_north_arrow(location = "tl", which_north = "true", 
                         style = north_arrow_orienteering,
                         height = unit(0.8, "cm"),
                         width = unit(0.8, "cm")) +
  labs(title = NULL, x = "Longitude", y = "Latitude") +
  theme_minimal()



#### effects of phase 1 after standardizing initial size ####

Growth_Data_filtered_pre <- Growth_Data_forR_pre %>%
  filter(Phase_1_rep_R != "Cont01") %>%
  filter(Phase_1_rep_R != "Both02") %>%
  filter(Phase_1_rep_R != "Hyp06") %>%
  filter(Phase_1_rep_R != "Warm02") #works

Growth_Data_filtered_pre <- Growth_Data_forR_pre %>%
  filter(Phase_1_rep_R != "Cont01") %>%
  filter(Phase_1_rep_R != "Both02") %>%
  filter(Phase_1_rep_R != "Hyp05") %>%
  filter(Phase_1_rep_R != "Warm02") #better

#effect of phase 1, after standardizing initial size, TISSUE
filt_pre_summary_stats_t <- Growth_Data_filtered_pre %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_pre = mean(Actual_tissue_pre_mg, na.rm = TRUE),
    se_pre = std.error(Actual_tissue_pre_mg, na.rm = TRUE))

#reorder effect of phase 1 Phase_1_treat
filt_pre_summary_stats_t$Phase_1_treat <- factor(filt_pre_summary_stats_t$Phase_1_treat, 
                                            levels = c("Cont", "Warm","Hyp", "Both"))
#plot TISSUE
ggplot(filt_pre_summary_stats_t, aes(x = Phase_1_treat, y = mean_pre, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.1, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18)+
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Tissue Mass (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend


#effect of phase 1, after standardizing initial size, SHELL
filt_pre_summary_stats_s <- Growth_Data_filtered_pre %>%
  group_by(Phase_1_treat) %>%
  mutate(
    mean_pre = mean(Actual_shell_pre_mg, na.rm = TRUE),
    se_pre = std.error(Actual_shell_pre_mg, na.rm = TRUE))

#reorder effect of phase 1 Phase_1_treat
filt_pre_summary_stats_s$Phase_1_treat <- factor(filt_pre_summary_stats_s$Phase_1_treat, 
                                                 levels = c("Cont", "Warm","Hyp", "Both"))
#plot Shell
ggplot(filt_pre_summary_stats_s, aes(x = Phase_1_treat, y = mean_pre, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.1, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18)+
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Tissue Mass (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend




#Shell area plots
#read in data from final merged csv
mergedarea_df <- ...

#dataset
merged_df_cleaned_pre <- mergedarea_df %>%
  filter(Exclude_pre_analysis != "Y" | is.na(Exclude_pre_analysis)) %>%
  mutate(Area_growth_mm2 = Area_post_mm2 - Area_pre_mm2,
         Feret_growth_mm = Feret_post_mm - Feret_pre_mm)

colnames(merged_df_cleaned_pre)

#SHELL PLOTS, significant interactions

##effects of phase 1, PHASE 1 DO & SHELL AREA
#pre area
pre_summary_stats_area <- merged_df_cleaned_pre %>%
  group_by(Phase_1_DO) %>%
  mutate(
    mean_pre = mean(Area_pre_mm2, na.rm = TRUE),
    se_pre = std.error(Area_pre_mm2, na.rm = TRUE))

pre_summary_stats_area$Phase_1_DO <- factor(pre_summary_stats_area$Phase_1_DO, 
                                         levels = c("Norm", "Hyp"))
#area plot with mean and SD
ggplot(pre_summary_stats_area, aes(x = Phase_1_DO, y = mean_pre, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 DO", y = "Shell Area (mm^2)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  theme(legend.position = "none") # Remove legend


##effects of phase 1, PHASE 1 DO & FERET
#pre feret
colnames(merged_df_cleaned_pre)
pre_summary_stats_feret <- merged_df_cleaned_pre %>%
  group_by(Phase_1_DO) %>%
  mutate(
    mean_pre = mean(Feret_pre_mm, na.rm = TRUE),
    se_pre = std.error(Feret_pre_mm, na.rm = TRUE))

pre_summary_stats_feret$Phase_1_DO <- factor(pre_summary_stats_feret$Phase_1_DO, 
                                            levels = c("Norm", "Hyp"))
#area plot with mean and SD
ggplot(pre_summary_stats_feret, aes(x = Phase_1_DO, y = mean_pre, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.05, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 DO", y = "Feret (mm)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  theme(legend.position = "none") # Remove legend










## model results
library(emmeans)
emm <- emmeans(Am1, ~ Phase_1_DO)
emm_df <- as.data.frame(emm)

ggplot(emm_df, aes(x = Phase_1_DO, y = emmean, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.05, position = position_dodge(0.9)) +
  theme_classic(base_size = 18) +
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 DO", y = "Model-adjusted Shell Area (mm^2)")


#### side by side comparison of raw data vs model ####

library(dplyr)
library(ggplot2)
library(emmeans)
library(patchwork) # For side-by-side plots

# Calculate raw means
pre_summary_stats_area <- merged_df_cleaned_pre %>%
  group_by(Phase_1_DO) %>%
  summarise(
    mean_pre = mean(Area_pre_mm2, na.rm = TRUE),
    se_pre = sd(Area_pre_mm2, na.rm = TRUE) / sqrt(sum(!is.na(Area_pre_mm2))),
    .groups = "drop"
  )

# Convert to factor with custom levels
pre_summary_stats_area$Phase_1_DO <- factor(pre_summary_stats_area$Phase_1_DO, 
                                            levels = c("Norm", "Hyp"))

# Raw mean plot
p1 <- ggplot(pre_summary_stats_area, aes(x = Phase_1_DO, y = mean_pre, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean_pre - se_pre, ymax = mean_pre + se_pre), 
                width = 0.05, position = position_dodge(0.9)) +
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 DO", y = "Raw mean shell area (mm^2)", title = "Raw Means") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  theme(legend.position = "none")

# Model and emmeans
Am1 <- lmer(Area_pre_mm2 ~ Phase_1_DO*Phase_1_temp+Actual_shell_pre_mg+
              (1|Phase_1_rep_R), data = merged_df_cleaned_pre, REML=TRUE)
emm <- emmeans(Am1, ~ Phase_1_DO)
emm_df <- as.data.frame(emm)

emm_df$Phase_1_DO <- factor(emm_df$Phase_1_DO, levels = c("Norm", "Hyp"))

# Model-adjusted mean plot
p2 <- ggplot(emm_df, aes(x = Phase_1_DO, y = emmean, color = Phase_1_DO)) +
  geom_point(size = 4, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE),
                width = 0.05, position = position_dodge(0.9)) +
  theme_classic(base_size = 16) +
  scale_color_manual(values = c("Hyp" = "darkmagenta", "Norm" = "seagreen")) +
  labs(x = "Phase 1 DO", y = "Model-adjusted mean shell area (mm^2)", title = "Model-Adjusted Means") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Norm" = "Normoxic")) +
  theme(legend.position = "none")

# Show plots side by side
p1 + p2

