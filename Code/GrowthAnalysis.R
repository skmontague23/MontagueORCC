#load packages
install.packages("plotrix")
install.packages("performance")
install.packages("lmtest")

library(lme4) #for linear mixed effects model
library(readr) #to read in data
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

Growth_Data_forR <- read_csv("/Users/sophiemontague/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Weight_Data/growth_phase2.1_weightsSKM.csv", 
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

#check levels
str(Growth_Data_forR)
levels(Growth_Data_forR$Phase_1_DO)
levels(Growth_Data_forR$Phase_1_temp)
levels(Growth_Data_forR$Phase_2.1_DO)
levels(Growth_Data_forR$Phase_2.1_temp)
levels(Growth_Data_forR$Phase_1_treat)
levels(Growth_Data_forR$Phase_2_treat)

#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 

#filter to not include data from dead oysters/ doubles in analysis
  #pre excludes doubles but keeps oysters that died in the second phase, as all oysters were alive for first measurements
  #full excludes doubles and dead oysters
Growth_Data_forR_pre <- Growth_Data_forR %>%
  filter(Exclude_pre_analysis != "Y" | is.na(Exclude_pre_analysis)) %>%
  mutate(prop_tissue_growth = Actual_tissue_growth_mg / Actual_tissue_pre_mg,
         prop_shell_growth = Actual_shell_growth_mg / Actual_shell_pre_mg,
         meat_yield = (Actual_tissue_pre_mg/(Dry_weight_pre)))
View(Growth_Data_forR_pre)
table(Growth_Data_forR$Exclude_pre_analysis, useNA = "ifany")
nrow(Growth_Data_forR_pre)

Growth_Data_forR_full <- Growth_Data_forR %>%
  filter(Exclude_all != "Y" | is.na(Exclude_all)) %>%
  mutate(prop_tissue_growth = Actual_tissue_growth_mg / Actual_tissue_pre_mg,
    prop_shell_growth = Actual_shell_growth_mg / Actual_shell_pre_mg, 
    tissuegrowthratio = Actual_tissue_post_mg/Actual_tissue_pre_mg,
    shellgrowthratio = Actual_shell_post_mg/Actual_shell_pre_mg,
    log_shell_growth_mg = log(Actual_shell_growth_mg),
    whole_growth_mg = Dry_weight_post - Dry_weight_pre,
    meat_yield_post = Actual_tissue_post_mg/Dry_weight_post,
    shell_yield_post = Actual_shell_post_mg/Dry_weight_post)

table(Growth_Data_forR$Exclude_all, useNA = "ifany")
table(Growth_Data_forR$Exclude_pre_analysis, useNA = "ifany")
View(Growth_Data_forR_full)


#### Effect of phase 1 at the start of phase 2 ####
#not excluding any replicates
##tissue growth (mg)
#Tm1, passes the levene test regardless, Q-Q plot looks better when log transformed, same ish results regardless
Tm1 <- lmer(log(Actual_tissue_pre_mg)~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Growth_Data_forR_pre, REML=TRUE)
Anova(Tm1, test="F", type="III")

  #post hocs
emmeans(Tm1,specs = pairwise ~ Phase_1_DO, adjust = "none")

  #actual value 
exp(5.25) #hyp
exp(5.37) #norm
(214.8629-190.5663)/214.8629

  #diagnostics
leveneTest(log(Actual_tissue_pre_mg)~Phase_1_treat, Growth_Data_forR_pre)
m1.e <- residuals(Tm1) 
qqnorm(m1.e)
qqline(m1.e)

plot(resid(Tm1), log(Growth_Data_forR_pre$Actual_tissue_pre_mg))

ggplot(Growth_Data_forR_pre) +
  aes(x = log(Actual_tissue_pre_mg)) +
  geom_histogram(bins = 30L, fill = "#112446") +
  theme_classic()


##shell growth (mg)
#Sm1, only passes levene's when not log transformed, not log transforming even though Q-Q plot looks exponential
Sm1 <- lmer(Actual_shell_pre_mg~ Phase_1_DO*Phase_1_temp +
                (1|Phase_1_rep_R), data = Growth_Data_forR_pre, REML=TRUE)
Anova(Sm1, test="F", type="III")

  #diagnostics
leveneTest(Actual_shell_pre_mg~Phase_1_treat, Growth_Data_forR_pre)
m1.e <- residuals(Sm1) 
qqnorm(m1.e)
qqline(m1.e)
ggplot(Growth_Data_forR_pre) +
  aes(x = log(Actual_shell_pre_mg)) +
  geom_histogram(bins = 30L, fill = "#112446") +
  theme_classic()

  #posthocs
emmeans(Sm1,specs = pairwise ~ Phase_1_DO, adjust = "none") 



##tissue:shell growth (mg)
colnames(Growth_Data_forR_pre)
TSm1 <- lmer(Ratio_tissue_shell_pre_mg~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Growth_Data_forR_pre, REML=TRUE)
Anova(TSm1, test="F", type="III")
  #diagnostics
leveneTest(Ratio_tissue_shell_pre_mg~Phase_1_treat, Growth_Data_forR_pre)
m1.e <- residuals(TSm1) 
qqnorm(m1.e)
qqline(residuals(TSm1))
AIC(TSm1)

  #posthocs
emmeans(TSm1, specs = pairwise ~ Phase_1_temp, adjust = "none") 
(0.643-0.618) /0.643

emmeans(TSm1, specs = pairwise ~ Phase_1_DO, adjust = "none")
(0.641-0.620) /0.641

##wet tissue weight / total weight = meat yield

colnames(Growth_Data_forR_pre_my)

qqnorm(Growth_Data_forR_pre$Ratio_tissue_shell_pre_mg)

TSm2 <- lmer(meat_yield ~ Phase_1_DO*Phase_1_temp +
               (1|Phase_1_rep_R), data = Growth_Data_forR_pre, REML=TRUE)
Anova(TSm2, test="F", type="III")
#diagnostics
leveneTest(meat_yield~Phase_1_treat, Growth_Data_forR_pre)
m1.e <- residuals(TSm2) 
qqnorm(m1.e)
qqline(residuals(TSm2))
AIC(TSm2)

#posthocs
emmeans(TSm2, specs = pairwise ~ Phase_1_temp, adjust = "none") 
emmeans(TSm2, specs = pairwise ~ Phase_1_DO, adjust = "none") 




#### FULL LMER MODELS ####
#phase 1 and phase 2

##tissue growth (mg)
m1 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_tissue_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_full, REML=TRUE)
Anova(m1, test="F", type="III")

#posthocs
emmeans(m1,specs = pairwise ~ Phase_2.1_DO, adjust = "none") 

(119.3-95.5)/119.3

#diagnostics
leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, Growth_Data_forR_full) 
m1.e <- residuals(m1) 
qqnorm(m1.e)
qqline(m1.e)


##shell growth (mg)

#m2 <- lmer(Actual_shell_growth_mg~ Phase_1_temp*Phase_1_DO*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
#             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
#             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_full, REML=TRUE)
#Anova(m2, test="F", type="III")
#leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, Growth_Data_forR_full)
#m3.e <- residuals(m2) 
#qqnorm(m3.e)
#qqline(m3.e)

#posthocs
#emmeans(m2,specs = pairwise ~ Phase_2.1_DO, adjust = "none") 
(289-225)/289

#log, using this model for analysis
Lm2 <- lmer(log_shell_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
              (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
              (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_full, REML=TRUE)
Anova(Lm2, test="F", type="III")

leveneTest(log_shell_growth_mg~Phase1_Phase2_treat, Growth_Data_forR_full)
m3.e <- residuals(Lm2)
qqnorm(m3.e)
qqline(m3.e)
emmeans(Lm2,specs = pairwise ~ Phase_2.1_DO, adjust = "none")
exp(5.23)
exp(5.51)
(247.1511-186.7928)/247.1511
emmeans(Lm2,specs = pairwise ~ Phase_1_DO*Phase_1_temp, adjust = "none")
exp(5.41)
exp(5.34)
(223.6316- 208.5127)/223.6316
exp(5.34)
exp(5.41)
(223.6316- 208.5127)/223.6316
emmeans(Lm2,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO, adjust = "none")
exp(5.22)
exp(5.54)


##tissue:shell growth 
m3 <- lmer(Ratio_tissue_shell_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_full, REML=TRUE)
Anova(m3, test="F", type="III")

#posthocs
emmeans(m3,specs = pairwise ~ Phase_1_temp:Phase_2.1_DO, adjust="none")
emmeans(m3,specs = pairwise ~ Phase_1_DO:Phase_1_temp:Phase_2.1_DO:Phase_2.1_temp, adjust="none")

View(Growth_Data_forR)

#diagnostics
leveneTest(Ratio_tissue_shell_mg~Phase1_Phase2_treat, Growth_Data_forR_full) 
m3.e <- residuals(m3) 
qqnorm(m3.e) #not quite normal but ancova is robust to non-normality

summary_stats_t_s <- Growth_Data_forR_full %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  summarise(
    mean_growth = mean(Ratio_tissue_shell_mg, na.rm = TRUE),
    se_growth = std.error(Ratio_tissue_shell_mg, na.rm = TRUE),
    .groups = "drop")


ggplot(summary_stats_t_s) +
  aes(
    x = Phase_1_treat,
    y = mean_growth,
    fill = Phase_2_treat
  ) +
  geom_point(position = position_dodge(width = 0.9), shape = 21, size = 3) +
  geom_errorbar(
    aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth),
    width = 0.2,
    position = position_dodge(width = 0.9)
  ) +
  scale_fill_hue(direction = 1) +
  theme_minimal()

##meat yield phase 2
m4 <- lmer(meat_yield_post~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_full, REML=TRUE)
Anova(m4, test="F", type="III")

#posthocs
emmeans(m4,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp, adjust="none")

View(Growth_Data_forR)

#diagnostics
leveneTest(meat_yield_post~Phase1_Phase2_treat, Growth_Data_forR_full) 
m3.e <- residuals(m4) 
qqnorm(m3.e) #not quite normal but ancova is robust to non-normality


##shell yield phase 2
m5 <- lmer(shell_yield_post~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_full, REML=TRUE)
Anova(m5, test="F", type="III")

#posthocs
emmeans(m5,specs = pairwise ~ Phase_2.1_DO, adjust="none")
emmeans(m5,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp, adjust="none")

View(Growth_Data_forR)

#diagnostics
leveneTest(shell_yield_post~Phase1_Phase2_treat, Growth_Data_forR_full) 
m3.e <- residuals(m5) 
qqnorm(m3.e) #not quite normal but ancova is robust to non-normality



#### Normalized by Whole Weight ####
#Phase 1
#### Effect of phase 1 at the start of phase 2 ####
#not excluding any replicates
##tissue growth (mg)
Tm1_norm <- lmer((Actual_tissue_pre_mg/Dry_weight_pre)~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Growth_Data_forR_pre, REML=TRUE)
Anova(Tm1_norm, test="F", type="III")

#post hocs
emmeans(Tm1_norm,specs = pairwise ~ Phase_1_DO, adjust = "none")
emmeans(Tm1_norm,specs = pairwise ~ Phase_1_temp, adjust = "none")

#diagnostics
leveneTest((Actual_tissue_pre_mg/Dry_weight_pre)~Phase_1_treat, Growth_Data_forR_pre)
m1.e <- residuals(Tm1_norm) 
qqnorm(m1.e)
qqline(m1.e)


## normalized shell growth (mg)
Sm1_norm <- lmer((Actual_shell_pre_mg/Dry_weight_pre)~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Growth_Data_forR_pre, REML=TRUE)
Anova(Sm1_norm, test="F", type="III")

#diagnostics
leveneTest((Actual_shell_pre_mg/Dry_weight_pre)~Phase_1_treat, Growth_Data_forR_pre)
m1.e <- residuals(Sm1_norm) 
qqnorm(m1.e)
qqline(m1.e)

#posthocs
emmeans(Sm1_norm,specs = pairwise ~ Phase_1_DO, adjust = "none") 
emmeans(Sm1_norm,specs = pairwise ~ Phase_1_temp, adjust = "none") 


##Phase 2
##normalized tissue growth (mg)
m1_norm <- lmer((Actual_tissue_growth_mg/whole_growth_mg)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), 
            data = Growth_Data_forR_full, REML=TRUE)
Anova(m1_norm, test="F", type="III")

#posthocs
emmeans(m1_norm,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO, adjust = "none")
emmeans(m1_norm,specs = pairwise ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO, adjust = "none") 

(119.3-95.5)/119.3

#diagnostics
leveneTest((Actual_tissue_growth_mg/whole_growth_mg)~Phase1_Phase2_treat, Growth_Data_forR_full)
m1.e <- residuals(m1_norm) 
qqnorm(m1.e)
qqline(m1.e)

#all, TISSUE
summary_stats_t <- Growth_Data_forR_full %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean((Actual_tissue_growth_mg/whole_growth_mg), na.rm = TRUE),
    se_growth = std.error((Actual_tissue_growth_mg/whole_growth_mg), na.rm = TRUE))

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
  labs(x = "Phase 1 Treatment", y = "Mean Normalized Tissue Growth (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend


##normalized shell growth (mg)
m2_norm <- lmer((Actual_shell_growth_mg/whole_growth_mg)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), 
           data = Growth_Data_forR_full, REML=TRUE)
Anova(m2_norm, test="F", type="III")

#posthocs
emmeans(m2_norm,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO, adjust = "none")
emmeans(m2_norm,specs = pairwise ~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp, adjust = "none")

(119.3-95.5)/119.3

#diagnostics
leveneTest((Actual_shell_growth_mg/whole_growth_mg)~Phase1_Phase2_treat, Growth_Data_forR_full)
m1.e <- residuals(m2_norm) 
qqnorm(m1.e)
qqline(m1.e)

#Phase_1_temp*Phase_2.1_DO graph
ss_normshell <- Growth_Data_forR_full %>%
  group_by(Phase_1_temp, Phase_2.1_DO) %>%
  summarise(
    mean_growth = mean((Actual_shell_growth_mg/whole_growth_mg), na.rm = TRUE),
    se_growth = std.error((Actual_shell_growth_mg/whole_growth_mg), na.rm = TRUE))

nrow(ss_normshell)


ss_normshell$Phase_1_temp <- factor(ss_normshell$Phase_1_temp, 
                                     levels = c("Ambient", "Warm"))
ss_normshell$Phase_2.1_DO <- factor(ss_normshell$Phase_2.1_DO, 
                                      levels = c("Norm","Hyp"))
View(ss_normshell)
ggplot(ss_normshell) +
  aes(x = Phase_1_temp, y = mean_growth) +
  geom_point(colour = "#112446") +
  labs(x = "Phase 1 Temp",
       y = "Norm Shell Growth (mg)",
       title = "Phase 2 DO")+
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), width = 0.2) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 18L,
                              hjust = 0.5),
    axis.title.y = element_text(size = 18L),
    axis.title.x = element_text(size = 18L)
  )+
  facet_wrap(vars(Phase_2.1_DO))+ scale_fill_hue()

#4 way interaction
#all, SHELL
summary_stats_s <- Growth_Data_forR_full %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean((Actual_shell_growth_mg/whole_growth_mg), na.rm = TRUE),
    se_growth = std.error((Actual_shell_growth_mg/whole_growth_mg), na.rm = TRUE))

#reorder Phase_1_treat and Phase_2_treat
summary_stats_s$Phase_1_treat <- factor(summary_stats_s$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))
summary_stats_s$Phase_2_treat <- factor(summary_stats_s$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))

#plot with mean and SD, 
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
  labs(x = "Phase 1 Treatment", y = "Mean Normalized Shell Growth (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend

####


##### Growth Proportion ####
##normalized tissue growth (mg)
m1_prop <- lmer(log(prop_tissue_growth +0.55)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+
                  (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                  (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), 
                data = Growth_Data_forR_full, REML=TRUE)
Anova(m1_prop, test="F", type="III")

#posthocs
emmeans(m1_prop,specs = pairwise ~ Phase_1_DO, adjust = "none")
emmeans(m1_prop,specs = pairwise ~ Phase_2.1_DO, adjust = "none")
emmeans(m1_prop,specs = pairwise ~ Phase_1_DO*Phase_1_temp, adjust = "none")

(119.3-95.5)/119.3

#diagnostics
leveneTest(log(prop_tissue_growth+0.55)~Phase1_Phase2_treat, Growth_Data_forR_full)
min(Growth_Data_forR_full$prop_tissue_growth)
m1.e <- residuals(m1_prop) 
qqnorm(m1.e)
qqline(m1.e)


##normalized shell growth (mg)
m2_prop <- lmer(log(prop_shell_growth)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+
                  (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                  (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), 
                data = Growth_Data_forR_full, REML=TRUE)
Anova(m2_prop, test="F", type="III")

#posthocs
emmeans(m2_prop,specs = pairwise ~ Phase_1_DO, adjust = "none")
emmeans(m2_prop,specs = pairwise ~ Phase_2.1_DO, adjust = "none")
emmeans(m2_prop,specs = pairwise ~ Phase_1_DO*Phase_1_temp, adjust = "none")


#diagnostics
leveneTest(log(prop_shell_growth)~Phase1_Phase2_treat, Growth_Data_forR_full)
m1.e <- residuals(m2_prop) 
qqnorm(m1.e)
qqline(m1.e)

####


####Post Phase 2 mass analysis####
##tissue mass post
m1_post <- lmer(Actual_tissue_post_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_tissue_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_full, REML=TRUE)
Anova(m1_post, test="F", type="III")

#posthocs
emmeans(m1_post,specs = pairwise ~ Phase_2.1_DO, adjust = "none")

#diagnostics
leveneTest(Actual_tissue_post_mg~Phase1_Phase2_treat, Growth_Data_forR_full) 
m1.e <- residuals(m1_post) 
qqnorm(m1.e)
qqline(m1.e)


##shell mass post
m2_post <- lmer(Actual_shell_post_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_tissue_pre_mg+
                  (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                  (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_full, REML=TRUE)
Anova(m2_post, test="F", type="III")

#posthocs
emmeans(m2_post,specs = pairwise ~ Phase_1_DO, adjust = "none")
emmeans(m2_post,specs = pairwise ~ Phase_2.1_DO, adjust = "none")
emmeans(m2_post,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO, adjust = "none")

#diagnostics
leveneTest(Actual_shell_post_mg~Phase1_Phase2_treat, Growth_Data_forR_full)
m1.e <- residuals(m2_post) 
qqnorm(m1.e)
qqline(m1.e)

#PLOT SHELL POST PHASE 2
summary_stats_s <- Growth_Data_forR_full%>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean_growth = mean(Actual_shell_post_mg, na.rm = TRUE),
    se_growth = std.error(Actual_shell_post_mg, na.rm = TRUE))

#reorder Phase_1_treat and Phase_2_treat
summary_stats_s$Phase_1_treat <- factor(summary_stats_s$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp",  "Both"))
summary_stats_s$Phase_2_treat <- factor(summary_stats_s$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp",  "Both"))

#plot with mean and SD, SHELL POST
ggplot(summary_stats_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")),
             scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Mean Shell Mass Post Phase 2 (mg)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend


##tissue:shell mass post
m3_post <- lmer((Actual_tissue_post_mg/Actual_shell_post_mg)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_tissue_pre_mg+
                  (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                  (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_full, REML=TRUE)
Anova(m3_post, test="F", type="III")

#posthocs
emmeans(m3_post,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO, adjust = "none")
emmeans(m3_post,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO, adjust = "none")

#diagnostics
leveneTest((Actual_tissue_post_mg/Actual_shell_post_mg)~Phase1_Phase2_treat, Growth_Data_forR_full)
m1.e <- residuals(m3_post) 
qqnorm(m1.e)
qqline(m1.e)

####


#### Growth Ratio ####
tissuegrowthratio
shellgrowthratio
## tissue growth (mg)
m1_ratio <- lmer((1/tissuegrowthratio)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+
                  (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                  (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), 
                data = Growth_Data_forR_full, REML=TRUE)
Anova(m1_ratio, test="F", type="III")

#posthocs
emmeans(m1_ratio,specs = pairwise ~ Phase_2.1_DO, adjust = "none")
emmeans(m1_ratio,specs = pairwise ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp, adjust = "none") 

#diagnostics
leveneTest((1/tissuegrowthratio)~Phase1_Phase2_treat, Growth_Data_forR_full)
m1.e <- residuals(m1_ratio) 
qqnorm(m1.e)
qqline(m1.e)

## shell growth (mg)
m2_ratio <- lmer((1/shellgrowthratio)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+
                   (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                   (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), 
                 data = Growth_Data_forR_full, REML=TRUE)
Anova(m2_ratio, test="F", type="III")

#posthocs
emmeans(m1_ratio,specs = pairwise ~ Phase_2.1_DO, adjust = "none")
emmeans(m1_ratio,specs = pairwise ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp, adjust = "none") 

#diagnostics
leveneTest((1/shellgrowthratio)~Phase1_Phase2_treat, Growth_Data_forR_full)
m1.e <- residuals(m2_ratio) 
qqnorm(m1.e)
qqline(m1.e)

####



####Standardize Initial Size####
#means for each rep calculated in GrowthAnalysisSizeClass.R
#second attempt works better to equalize size at the start of phase 2, this one works with no covariate
Growth_Data_filtered <- Growth_Data_forR_full %>%
  filter(Phase_1_rep_R != "Cont01") %>%
  filter(Phase_1_rep_R != "Both02") %>%
  filter(Phase_1_rep_R != "Hyp06") %>%
  filter(Phase_1_rep_R != "Warm02") #works, this is the one I've been working with

Growth_Data_filtered <- Growth_Data_forR_full %>%
  filter(Phase_1_rep_R != "Cont01") %>%
  filter(Phase_1_rep_R != "Both02") %>%
  filter(Phase_1_rep_R != "Hyp05") %>%
  filter(Phase_1_rep_R != "Warm02") #better?



#check the effect of phase 1 at the start of phase 2
# excluding replicates, there should now be no effect
##tissue growth (mg)
Tm2 <- lmer(Actual_tissue_pre_mg~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Growth_Data_filtered_pre, REML=TRUE)
Anova(Tm2, test="F", type="III")
#diagnostics
leveneTest(Actual_tissue_pre_mg~Phase1_Phase2_treat, Growth_Data_filtered_pre)
m1.e <- residuals(Tm2) 
qqnorm(m1.e)
qqline(m1.e)

##shell growth (mg)
Sm2 <- lmer(Actual_shell_pre_mg~ Phase_1_DO*Phase_1_temp + 
              (1|Phase_1_rep_R), data = Growth_Data_filtered_pre, REML=TRUE)
Anova(Sm2, test="F", type="III")
#diagnostics
leveneTest(Actual_shell_pre_mg~Phase1_Phase2_treat, Growth_Data_filtered)
m1.e <- residuals(Sm2) 
qqnorm(m1.e)
qqline(m1.e)

##tissue:shell growth (mg)
TSm2 <- lmer(Ratio_tissue_shell_pre_mg~ Phase_1_DO*Phase_1_temp +
               (1|Phase_1_rep_R), data = Growth_Data_filtered_pre, REML=TRUE)
Anova(TSm2, test="F", type="III")
#diagnostics
leveneTest(Ratio_tissue_shell_pre_mg~Phase1_Phase2_treat, Growth_Data_filtered)
m1.e <- residuals(TSm2) 
qqnorm(m1.e)
qqline(m1.e)



#run tissue model excluding replicates from phase 1 treatments
m1_rep <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_tissue_pre_mg+
                 (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                 (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_filtered, REML=TRUE)
Anova(m1_rep, test="F", type="III")

#posthocs
emmeans(m1_rep,specs = pairwise ~ Phase_2.1_DO, adjust="none")
(118.5 - 96.9)/118.5

emmeans(m1_rep,specs = pairwise ~ Phase_2.1_temp*Phase_2.1_DO, adjust="none")
(126.0 - 86.1)/126.0

emmeans(m1_rep,specs = pairwise ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp, adjust="none")

#diagnostics
leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, Growth_Data_filtered) #passes
m1.e <- residuals(m1_rep) 
qqnorm(m1.e)


#run shell model excluding replicates from phase 1 treatments
Lm2_rep <- lmer(log(Actual_shell_growth_mg)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
                  (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                  (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_filtered, REML=TRUE)
Anova(Lm2_rep, test="F", type="III")

#posthocs
emm1 <- emmeans(Lm2_rep,specs = pairwise ~ Phase_2.1_DO, adjust="none")
emm1$emmeans 
emm1$contrasts 

(286-226)/286

#diagnostics
leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, Growth_Data_filtered) 
m1.e <- residuals(Lm2_rep)
qqnorm(m1.e)
qqline(m1.e)

View(Growth_Data_filtered)

#T:S model excluding reps
Lm3_rep <- lmer(Ratio_tissue_shell_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+ 
                  (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                  (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_filtered, REML=TRUE)
Anova(Lm3_rep, test="F", type="III")


#posthocs
emmeans(Lm3_rep,specs = pairwise ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO, adjust = "none") 
#More tissue:shell growth in hypoxic conditions than in normoxic conditions


#diagnostics
leveneTest(Ratio_tissue_shell_mg~Phase1_Phase2_treat, Growth_Data_filtered) #passes
m1.e <- residuals(Lm3_rep) 
qqnorm(m1.e)





#Try standardizing by filtering out just the top 20 
Growthdata_20out_pre<- Growth_Data_forR_pre%>%
  group_by(Phase_1_treat) %>%
  arrange(if_else(Phase_1_treat %in% c("Cont", "Warm"), desc(Dry_weight_pre), Dry_weight_pre)) %>%
  mutate(row_num = row_number()) %>%
  filter(
    !(Phase_1_treat %in% c("Cont", "Warm") & row_num <= 6),
    !(Phase_1_treat %in% c("Hyp", "Both") & row_num <= 6)
  ) %>%
  select(-row_num) %>%
  ungroup()

Growthdata_20out<- Growth_Data_forR_full%>%
  group_by(Phase_1_treat) %>%
  arrange(if_else(Phase_1_treat %in% c("Cont", "Warm"), desc(Dry_weight_pre), Dry_weight_pre)) %>%
  mutate(row_num = row_number()) %>%
  filter(
    !(Phase_1_treat %in% c("Cont", "Warm") & row_num <= 6),
    !(Phase_1_treat %in% c("Hyp", "Both") & row_num <= 6)
  ) %>%
  select(-row_num) %>%
  ungroup()
View(Growthdata_20out)


#check the effect of phase 1 at the start of phase 2
# excluding INDIVIDUAL OYSTERS
##tissue growth (mg)
Tm2 <- lmer(Actual_tissue_pre_mg~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = Growthdata_20out, REML=TRUE)
Anova(Tm2, test="F", type="III")
#diagnostics
leveneTest(Actual_tissue_pre_mg~Phase1_Phase2_treat, Growthdata_20out)
m1.e <- residuals(Tm2) 
qqnorm(m1.e)
qqline(m1.e)

##shell growth (mg)
Sm2 <- lmer(Actual_shell_pre_mg~ Phase_1_DO*Phase_1_temp + 
              (1|Phase_1_rep_R), data = Growthdata_20out, REML=TRUE)
Anova(Sm2, test="F", type="III")
#diagnostics
leveneTest(log(Actual_shell_pre_mg)~Phase1_Phase2_treat, Growthdata_20out)
m1.e <- residuals(Sm2) 
qqnorm(m1.e)
qqline(m1.e)

##tissue:shell growth (mg)
TSm2 <- lmer(Ratio_tissue_shell_pre_mg~ Phase_1_DO*Phase_1_temp +
               (1|Phase_1_rep_R), data = Growthdata_20out, REML=TRUE)
Anova(TSm2, test="F", type="III")
#diagnostics
leveneTest(Ratio_tissue_shell_pre_mg~Phase1_Phase2_treat, Growthdata_20out)
m1.e <- residuals(TSm2) 
qqnorm(m1.e)
qqline(m1.e)


#run tissue model excluding replicates from phase 1 treatments
m1_rep <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_tissue_pre_mg+
                 (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                 (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growthdata_20out, REML=TRUE)
Anova(m1_rep, test="F", type="III")

#posthocs
emmeans(m1_rep,specs = pairwise ~ Phase_2.1_DO, adjust="none")
(119- 94)/119

#diagnostics
leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, Growthdata_20out) #passes
m1.e <- residuals(m1_rep) 
qqnorm(m1.e)


#run shell model excluding replicates from phase 1 treatments
Lm2_rep <- lmer(log(Actual_shell_growth_mg)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
                  (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                  (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growthdata_20out, REML=TRUE)
Anova(Lm2_rep, test="F", type="III")

#posthocs
emmeans(Lm2_rep,specs = pairwise ~ Phase_2.1_DO, adjust="none")
exp(5.50); exp(5.24)
(244.6919-188.6701)/244.6919
emmeans(Lm2_rep,specs = pairwise ~ Phase_1_DO:Phase_1_temp, adjust="none")

#diagnostics
leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, Growthdata_20out) 
m1.e <- residuals(m1_rep)
qqnorm(m1.e)

View(Growth_Data_filtered)

#T:S model
Lm3_rep <- lmer(Ratio_tissue_shell_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+ 
                  (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                  (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growthdata_20out, REML=TRUE)
Anova(Lm3_rep, test="F", type="III")


#posthocs
emmeans(Lm3_rep,specs = pairwise ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO, adjust = "none") 
#More tissue:shell growth in hypoxic conditions than in normoxic conditions


#diagnostics
leveneTest(Ratio_tissue_shell_mg~Phase1_Phase2_treat, Growth_Data_filtered) #passes
m1.e <- residuals(Lm3_rep) 
qqnorm(m1.e)


