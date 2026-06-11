library(readr)
library(tidyverse)

getwd()
setwd("~/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Nutrient_Data/Phase2_Nutrient")

nutrientdata <- read_csv("~/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Nutrient_Data/Phase2_Nutrient/Phase2_nutrient_working.csv")

preprocessing_info <- read_csv("~/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Nutrient_Data/NutrientPreprocessing/Oyster_Nutrient_Processing_Phase2.1.csv")

colnames(nutrientdata)
colnames(preprocessing_info)

nutrientdata <- nutrientdata%>%
  janitor::clean_names()%>%
  rename(Vial_num = vial_num)%>%
  mutate(Vial_num = as.character(Vial_num))

preprocessing_info <- preprocessing_info %>%
  mutate(Vial_num = as.character(Vial_num))%>%
  mutate(Tissue_type = as.character(Tissue_type))%>%
  rename(filt_out_preprocessing = filt_out)



?left_join()

joined<- full_join(
  nutrientdata,
  preprocessing_info,
  by = "Vial_num",
  copy = FALSE,
  suffix = c(".x", ".y"),
  keep = NULL
)

View(joined)

#write_csv(joined, "~/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Nutrient_Data/Phase2_Nutrient/Phase2_nutrient_merged.csv")


#read in nutrient data
joined <- read.csv("Phase2_nutrient_merged.csv")
View(joined)

#Filter samples to make a tissue and shell dataset
S1 <- joined%>%
  filter(joined$Tissue_type == "S")


T1 <- joined%>%
  filter(joined$Tissue_type == "T")

View(S1)

repack <- S1 %>%
  filter(filt_out == "Y" | n_exclude == "Y")




#load packages
library(ggplot2)
library(dplyr)
library(plotrix) #for standard error

#### NITROGEN STORAGE PLOTS #### in tissue
#All possible treatment combinations

colnames(T1)

TN <- T1 %>%
  filter(is.na(filt_out) | filt_out != "Y",
         is.na(n_exclude) | n_exclude != "Y")



summary_stats <- TN %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  mutate(
    mean = mean(wt_percent_n, na.rm = TRUE),
    se = std.error(wt_percent_n, na.rm = TRUE))

View(summary_stats)

#reorder Phase_1_treat and Phase_2_treat
summary_stats$Phase_1_treat <- factor(summary_stats$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))
summary_stats$Phase_2_treat <- factor(summary_stats$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))

#plot with mean and SD, TISSUE
ggplot(summary_stats, aes(x = Phase_1_treat, y = mean, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic(base_size = 18) +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), 
             labeller = as_labeller(c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")), 
             scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_manual(values = c("Hyp" = "steelblue3", "Warm" = "palevioletred", "Cont" = "burlywood3", "Both" = "plum3")) +
  labs(x = "Phase 1 Treatment", y = "Nitrogen in Tissue (%)") +
  scale_x_discrete(labels = c("Hyp" = "Hypoxic", "Cont" = "Control", "Warm" = "Warm", "Both" = "Both")) +
  theme(legend.position = "none") # Remove legend





#models

#phase 1 and phase 2

options(contrasts = c("contr.sum","contr.sum")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 

##tissue growth (mg)
m1 <- lmer(wt_percent_n~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_tissue_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = TN, REML=TRUE)
Anova(m1, test="F", type="III")

#posthocs
emmeans(m1,specs = pairwise ~ Phase_2.1_DO, adjust = "none") 

(119.3-95.5)/119.3

#diagnostics
leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, Growth_Data_forR_full) 
m1.e <- residuals(m1) 
qqnorm(m1.e)
qqline(m1.e)
