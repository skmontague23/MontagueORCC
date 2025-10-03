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

R.version.string
citation("lme4")

#Read in dataset, set column types
getwd()
setwd("/Users/sophiemontague/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Weight_Data")

Year2_Growth <- read_csv("/Users/sophiemontague/Desktop/MontagueORCC_repo/MontagueORCC/Year_2_Data/growth_YEAR2_weights_merged.csv", 
                             col_types = cols(Phase_1_temp = col_factor(), 
                                              Phase_1_DO =col_factor(), 
                                              Phase_2.1_temp =col_factor(), 
                                              Phase_2.1_DO=col_factor(), 
                                              Phase_1_rep =col_factor(), 
                                              Phase_2_rep =col_factor(),
                                              Phase_1_treat = col_factor(),
                                              Phase_2_treat = col_factor(),
                                              Phase_2_rep_R = col_factor(),
                                              Phase_1_rep_R = col_factor(),
                                              double_single_dead_year2 = col_factor()))
colnames(Year2_Growth)
#filter out dead and untagged
Year2_Growth_clean <- Year2_Growth%>%
  drop_na(Actual_tissue_growth_year2_mg) %>%
  filter(Exclude_year2 != "Y" | is.na(Exclude_year2)) %>%
  mutate(
    Actual_tissue_growth_year2_mg = as.numeric(Actual_tissue_growth_year2_mg),
    Actual_shell_growth_year2_mg = as.numeric(Actual_shell_growth_year2_mg),
    Ratio_TS_year2 = as.numeric(Ratio_TS_year2))

Year2_Growth_clean %>%
  count(Phase_1_treat, Phase_2_treat)
Year2_Growth_clean$double_single_dead_year2
#plot with mean and SD, 16 point TISSUE
summary_stats <- Year2_Growth_clean %>%
  drop_na(Ratio_TS_year2) %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Dry_weight_year2, na.rm = TRUE),
    se_growth = std.error(Dry_weight_year2, na.rm = TRUE))

#reorder Phase_1_treat and Phase_2_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))
summary_stats$Phase_2_treat <- factor(summary_stats$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), 
             labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")), 
             scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Whole weight (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend

View(Year2_Growth_clean)
#counting dead oysters
unique(Year2_Growth_clean$Notes_year2)

dead_counts <- Year2_Growth_clean %>%
  drop_na(double_single_dead_year2)%>%
  filter(double_single_dead_year2 == "Dead") %>%
  count(Phase_1_treat, Phase_2_treat)


# Plot as bar chart
ggplot(dead_counts, aes(x = Phase_1_treat, y = n, fill = Phase_1_treat)) +
  geom_col(position = "dodge") +
  facet_wrap(~Phase_2_treat) +
  labs(
    title = "Count of Dead Oysters by Treatments",
    x = "Phase 1 Treatment",
    y = "Number of Dead Oysters"
  ) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none")






#check levels
str(Year2_Growth_clean)
levels(Year2_Growth_clean$Phase_1_DO)
levels(Year2_Growth_clean$Phase_1_temp)
levels(Year2_Growth_clean$Phase_2.1_DO)
levels(Year2_Growth_clean$Phase_2.1_temp)
levels(Year2_Growth_clean$Phase_1_treat)
levels(Year2_Growth_clean$Phase_2_treat)

##QAQC
#look at minimum tissue growth
sorted_vector <- sort(Year2_Growth_clean$Actual_tissue_growth_year2_mg)
lowest_30_values <- sorted_vector[1:30]
lowest_30_values

#look at minimum shell growth
min(Year2_Growth_clean$Actual_shell_growth_year2_mg)
sorted_vector <- sort(Year2_Growth_clean$Actual_shell_growth_year2_mg)
lowest_30_values <- sorted_vector[1:30]
lowest_30_values



#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 


##tissue growth (mg), having the Actual_tissue_pre_mg as the covariate makes different post hocs significant
m1 <- lmer(Actual_tissue_growth_year2_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_tissue_growth_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Year2_Growth_clean, REML=TRUE)
Anova(m1, test="F", type="III")

#posthocs
emmeans(m1,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp, adjust = "none")
emmeans(m1,specs = pairwise ~ Phase_1_temp:Phase_2.1_temp:Phase_2.1_DO, adjust = "none")

#diagnostics
leveneTest(Actual_tissue_growth_year2_mg~Phase_1_treat*Phase_2_treat, Year2_Growth_clean) 
m1.e <- residuals(m1) 
qqnorm(m1.e)
qqline(m1.e)

#plot with mean and SD, 16 point TISSUE
summary_stats_t <- Year2_Growth_clean %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Actual_tissue_growth_year2_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_year2_mg, na.rm = TRUE))

#reorder Phase_1_treat and Phase_2_treat
summary_stats_t$Phase_1_treat <- factor(summary_stats_t$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))
summary_stats_t$Phase_2_treat <- factor(summary_stats_t$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))

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
  labs(x = "Phase 1 Treatment", y = "Mean Normalized Tissue Growth (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend


#plot with phase 1 and 2 temp, TISSUE
summary_stats_t <- Year2_Growth_clean %>%
  filter(!is.na(Actual_tissue_growth_year2_mg)) %>%
  group_by(Phase_1_temp, Phase_2.1_temp) %>%
  mutate(
    mean_growth = mean(Actual_tissue_growth_year2_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_year2_mg, na.rm = TRUE))

#reorder Phase_1_treat and Phase_2_treat
summary_stats_t$Phase_1_temp <- factor(summary_stats_t$Phase_1_temp, 
                                        levels = c("Ambient", "Warm"))
summary_stats_t$Phase_2.1_temp <- factor(summary_stats_t$Phase_2.1_temp, 
                                        levels = c("Ambient", "Warm"))

ggplot(summary_stats_t, aes(x = Phase_1_temp, y = mean_growth, color = Phase_1_temp)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2.1_temp), 
             labeller = as_labeller(c("Warm" = "Warm", "Ambient" = "Ambient")), 
             scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  labs(x = "Phase 1 Treatment", y = "Year 2 Tissue Growth (mg)") +
  theme(legend.position = "none") # Remove legend




## Shell growth
##tissue growth (mg)
pos_shell<- Year2_Growth_clean%>%
  filter(Actual_shell_growth_year2_mg >= 0)

m2 <- lmer(Actual_shell_growth_year2_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_tissue_post_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = pos_shell, REML=TRUE)
Anova(m2, test="F", type="III")

#posthocs
emmeans(m2,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp, adjust = "none") 

#diagnostics
leveneTest(Actual_shell_growth_year2_mg~Phase_1_treat*Phase_2_treat, Year2_Growth_clean) 
m1.e <- residuals(m2) 
qqnorm(m1.e)
qqline(m1.e)


summary_stats_s <- Year2_Growth_clean %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Actual_shell_growth_year2_mg, na.rm = TRUE),
    se_growth = std.error(Actual_shell_growth_year2_mg, na.rm = TRUE))

#reorder Phase_1_treat and Phase_2_treat
summary_stats_s$Phase_1_treat <- factor(summary_stats_s$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))
summary_stats_s$Phase_2_treat <- factor(summary_stats_s$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(summary_stats_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 20) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), 
             labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")), 
             scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Year 2 Tissue Growth (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend


#####tissue:Shell growth (mg)
colnames(Year2_Growth_clean)
m3 <- lmer(Ratio_TS_year2~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Year2_Growth_clean, REML=TRUE)
Anova(m3, test="F", type="III")

#posthocs
emmeans(m3,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp, adjust = "none") 

#diagnostics
leveneTest(Ratio_TS_year2~Phase_1_treat*Phase_2_treat, Year2_Growth_clean) 
m1.e <- residuals(m3) 
qqnorm(m1.e)
qqline(m1.e)
