#set working directory
getwd()
setwd("/Users/sophiemontague/Desktop/MontagueORCC_repo/MontagueORCC/Oyster_Weight_Data")

library(lme4) #for linear mixed effects model
library(readr) #to read in data
library(emmeans) #for post hocs
library(esquisse) #interface for building plots
library(dplyr) #for data wrangling
library(nlme) #for linear mixed effects model
library(car) #for the Levene test
library(ggplot2) #for graphs
library(plotrix) #for standard error
library(readxl) #for excel workbooks

#load data
area_pre_df <- read_csv("growth_phase2.1_areaPRE.csv")
area_post_df <- read_csv("growth_phase2.1_areaPOST.csv")
View(area_pre_df)
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

#merge PRE AND POST PHASE 2 area data frame with growth data frame
colnames(area_pre_df)
mergedpre_df <- Growth_Data_forR %>%
  select(Sample_Name,
         Phase_1_DO,
         Phase_1_temp,
         Phase_2.1_DO,
         Phase_2.1_temp,
         Phase_1_treat,
         Phase_2_treat,
         Phase_1_rep_R,
         Phase_2_rep_R,
                 Actual_shell_pre_mg,
                 Actual_tissue_pre_mg,
                 Actual_shell_post_mg,
                 Actual_tissue_post_mg,
                 Actual_shell_growth_mg,
                 Actual_tissue_growth_mg,
                 Exclude_all,
                 Exclude_pre_analysis,
                 Notes_pre,
                 Notes_post)%>%
  left_join(
    area_pre_df %>%
      select(Sample_Name, Area_pre_mm2, Feret_pre_mm, NotesShell_Pre),
    by = "Sample_Name")%>%
  
  left_join(
    area_post_df %>%
      select(Sample_Name, Area_post_mm2, Feret_post_mm, ShellNotes_post),
    by = "Sample_Name")

View(mergedpre_df)
nrow(Growth_Data_forR)
nrow(mergedpre_df)

#check which oysters have missing area data PRE phase 2
missing_area_pre <- mergedpre_df %>%
  filter(is.na(Area_pre_mm2))
View(missing_area_pre) #8 oysters missing data pre phase 2

#check which oysters have missing area data POST phase 2
missing_area_post <- mergedpre_df %>%
  filter(is.na(Area_post_mm2))
View(missing_area_post) #36 oysters missing data pre phase 2, now 32 with matching area data to oyster tags, now down to 27 with matching more tags with the rest of the unmatched area data

#saving a csv of missing areas to annotate
write.csv(missing_area_post, "~/Desktop/missing_area_post.csv", row.names = FALSE)

#Visualize mass x area from the two datasets to make sure there aren't outiers from merging
ggplot(merged_df_cleaned) +
  aes(x = Actual_shell_pre_mg, y = Area_pre_mm2) +
  geom_point(colour = "blue4") +
  theme_classic()

ggplot(merged_df_cleaned) +
  aes(x = Area_post_mm2, y = Area_pre_mm2) +
  geom_point(colour = "blue4") +
  theme_classic()

ggplot(merged_df_cleaned) +
  aes(x = Feret_growth_mm, y = Area_growth_mm2) +
  geom_point(colour = "blue4") +
  theme_classic()

#Make a list of occurences the shell area growth is negative
negative_shellarea <- merged_df_cleaned %>%
  filter(Area_growth_mm2 < 0)
nrow(negative_shellarea) #92 is a lot, after excluding the ones with negative shell mass

#now positive
positive_shellarea <- merged_df_cleaned %>%
  filter(Area_growth_mm2 > 0)

#Make a list of occurences the Feret growth is negative
negative_feret <- merged_df_cleaned %>%
  filter(Feret_growth_mm < 0)
nrow(negative_feret) #96 is a lot

#now positive
positive_feret <- merged_df_cleaned %>%
  filter(Feret_growth_mm > 0)

#get the overlap between negative area and Feret sample names

# Get the vectors of Sample_Name from each dataframe
neg_shell_names <- negative_shellarea$Sample_Name
neg_feret_names <- negative_feret$Sample_Name

# Find the intersection and count overlap
overlap_names <- intersect(neg_shell_names, neg_feret_names)
length(overlap_names) #64


## run the Area Growth (mm^2) model on just the positive data
PAm2 <- lmer(Area_growth_mm2 ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
              (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
              (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive_shellarea, REML=TRUE)
Anova(PAm2, test="F", type="III")

View(negative_shellarea)

## run the Feret (mm) model on just the positive data
PFm2 <- lmer(Feret_growth_mm ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
               (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
               (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive_feret, REML=TRUE)
Anova(PFm2, test="F", type="III")

View(negative_shellarea)

#check which rows were mismatched between datasets, should add more onto the 8 that are missing from before
not_in_mergeddf <- anti_join(Growth_Data_forR, mergedgrowth_df, by = "Sample_Name")
View(not_in_mergeddf) #still only missing 8 oysters
not_in_growth_post <- anti_join(mergedgrowth_df, Growth_Data_forR, by = "Sample_Name")
View(not_in_growth_post) #None in shell area that are not in growth


#Edit dataset to exclude doubles and dead ones, calculate growth in area and Feret
merged_df_cleaned <- mergedpre_df %>%
  filter(Exclude_all != "Y" | is.na(Exclude_all)) %>%
  mutate(Area_growth_mm2 = Area_post_mm2 - Area_pre_mm2,
         Feret_growth_mm = Feret_post_mm - Feret_pre_mm)


#Visualize mass x area from the two datasets to make sure there aren't outiers from merging
ggplot(mergedgrowth_df) +
  aes(x = Actual_shell_post_mg, y = Area_post_mm2) +
  geom_point(colour = "blue4") +
  theme_classic()
#REMOVES 38 ROWS, WHY?

missing_area <- merged_df_cleaned %>%
  filter(is.na(Area_growth_mm2))
View(missing_area)

#Visualize GROWTH mass x area from the two datasets to make sure there aren't outiers from merging
ggplot(merged_df_cleaned) +
  aes(x = Actual_shell_growth_mg, y = Area_growth_mm2) +
  geom_point(colour = "blue4") +
  theme_classic()
#SOME OUTLIERS TO ADDRESS

#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 


#run model on shell area size at the start of phase 2
#effects of phase 1
Am1 <- lmer( Area_pre_mm2 ~ Phase_1_DO*Phase_1_temp+Actual_shell_pre_mg+
             (1|Phase_1_rep_R), data = merged_df_cleaned, REML=TRUE)
Anova(Am1, test="F", type="III")

#posthocs
emmeans(Am1,specs = pairwise ~ Phase_1_DO, adjust = "none")

#diagnostics
leveneTest(Area_pre_mm2~ Phase_1_treat*Phase_2_treat, merged_df_cleaned)
m1.e <- residuals(Am1) 
qqnorm(m1.e)
qqline(m1.e)

AIC(Am1)


#run model on feret
Fm1 <- lmer(Feret_pre_mm~ Phase_1_DO*Phase_1_temp+Actual_shell_pre_mg+
              (1|Phase_1_rep_R), data = merged_df_cleaned, REML=TRUE)
Anova(Fm1, test="F", type="III")

#posthocs
emmeans(Fm1,specs = pairwise ~ Phase_1_DO, adjust = "none") 

#diagnostics
leveneTest(log(Feret_pre_mm)~ Phase_1_treat*Phase_2_treat, merged_df_cleaned) 
m1.e <- residuals(Fm1) 
qqnorm(m1.e)
qqline(m1.e)


####FULL LMER MODEL####

## Area Growth (mm^2)
Am2 <- lmer(Area_growth_mm2 ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned, REML=TRUE)
Anova(Am2, test="F", type="III")

#posthocs
emmeans(Am2,specs = pairwise ~ XYZ, adjust = "none") 

#diagnostics
leveneTest(Area_growth_mm2~Phase_1_treat*Phase_2_treat, merged_df_cleaned) 
m1.e <- residuals(Am2) 
qqnorm(m1.e)
qqline(m1.e)


## Feret Growth (mm)
Fm2 <- lmer(Feret_growth_mm ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
              (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
              (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned, REML=TRUE)
Anova(Fm2, test="F", type="III")

#posthocs
emmeans(Am2,specs = pairwise ~ XYZ, adjust = "none") 

#diagnostics
leveneTest(Area_growth_mm2~Phase_1_treat*Phase_2_treat, merged_df_cleaned) 
m1.e <- residuals(Am2) 
qqnorm(m1.e)
qqline(m1.e)


