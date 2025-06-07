#set working directory
getwd()
setwd("/Users/sophiemontague/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Weight_Data")

library(readr)
library(dplyr)

#merge area data frame with growth data frame
shellarea_df <- read_csv("OysterGrowth2024_area.csv")
View(shellarea_df)
Growth_Data_forR <- read_csv("/Users/sophiemontague/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Weight_Data/growth_phase2.1_weightsSKM.csv", 
                             col_types = cols(Phase_1_temp = col_factor(), 
                                              Phase_1_DO =col_factor(), 
                                              Phase_2.1_temp =col_factor(), 
                                              Phase_2.1_DO=col_factor(), 
                                              Phase_1_rep =col_factor(), 
                                              Phase_2_rep =col_factor(),
                                              Phase_1_treat = col_factor(),
                                              Phase_2_treat = col_factor(),
                                              Phase_2_rep_R = col_factor(),
                                              Phase_1_rep_R = col_factor()))

?merge
merged_df <- merge(Growth_Data_forR, shellarea_df, by = "Sample_Name")
View(merged_df)
nrow(Growth_Data_forR)

#check which rows were mismatched between datasets
not_in_shellarea <- anti_join(Growth_Data_forR, shellarea_df, by = "Sample_Name")
View(not_in_shellarea) #missing 16 oysters
not_in_growth <- anti_join(shellarea_df, Growth_Data_forR, by = "Sample_Name")
View(not_in_growth) #None in shell area that are not in growth



#Edit dataset to exclude doubles and dead ones
merged_df_cleaned <- merged_df %>%
  filter(Exclude_all != "Y" | is.na(Exclude_all)) %>%
  mutate(prop_tissue_growth = Actual_tissue_growth_mg / Actual_tissue_pre_mg,
         prop_shell_growth = Actual_shell_growth_mg / Actual_shell_pre_mg, 
         tissuegrowthratio = Actual_tissue_post_mg/Actual_tissue_pre_mg,
         shellgrowthratio = Actual_shell_post_mg/Actual_shell_pre_mg,
         log_shell_growth_mg = log(Actual_shell_growth_mg))


#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 


#run model on shell area growth
Am1 <- lmer(log(`Area_mm^2`)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned, REML=TRUE)
Anova(Am1, test="F", type="III")

#posthocs
emmeans(Am1,specs = pairwise ~ Phase_1_DO, adjust = "none") 

#diagnostics
leveneTest(log(`Area_mm^2`)~Phase1_Phase2_treat, merged_df_cleaned) 
m1.e <- residuals(Am1) 
qqnorm(m1.e)
qqline(m1.e)

AIC(Am1)
