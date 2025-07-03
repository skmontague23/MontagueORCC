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
mergedarea_df <- Growth_Data_forR %>%
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
      select(Sample_Name, Area_pre_mm2, Feret_pre_mm, AreaNotes_Pre),
    by = "Sample_Name")%>%
  
  left_join(
    area_post_df %>%
      select(Sample_Name, Area_post_mm2, Feret_post_mm, AreaNotes_Post),
    by = "Sample_Name")

View(mergedarea_df)
nrow(Growth_Data_forR)
nrow(mergedarea_df)

#check which oysters have missing area data PRE phase 2
colnames(mergedarea_df)
cleanedarea_df <- mergedarea_df %>%
  filter(Exclude_pre_analysis != "Y" | is.na(Exclude_pre_analysis))

missing_area_pre <- cleanedarea_df %>%
  filter(is.na(Area_pre_mm2))
View(missing_area_pre) #6 oysters missing data pre phase 2

#check which oysters have missing area data POST phase 2
cleanedarea_df <- mergedarea_df %>%
  filter(Exclude_all != "Y" | is.na(Exclude_all))

missing_area_post <- cleanedarea_df %>%
  filter(is.na(Area_post_mm2))
View(missing_area_post) #36 oysters missing data pre phase 2, now 32 with matching area data to oyster tags, now down to 23 with matching more tags with the rest of the unmatched area data

#saving csv of missing areas to annotate while correcting
  #pre
write.csv(missing_area_pre, "~/Desktop/missing_area_pre.csv", row.names = FALSE)
  #post
write.csv(missing_area_post, "~/Desktop/missing_area_post.csv", row.names = FALSE)



##Working with the data
#get rid of the missing and dead 
#Edit dataset to exclude doubles and dead ones, calculate growth in area and Feret
merged_df_cleaned_all <- mergedarea_df %>%
  filter(Exclude_all != "Y" | is.na(Exclude_all)) %>%
  mutate(Area_growth_mm2 = Area_post_mm2 - Area_pre_mm2,
         Feret_growth_mm = Feret_post_mm - Feret_pre_mm)

View(mergedarea_df)
  
#Visualize mass x area from the two datasets to make sure there aren't outiers from merging
merged_df_cleaned_pre <- mergedarea_df %>%
  filter(Exclude_pre_analysis != "Y" | is.na(Exclude_pre_analysis)) %>%
  mutate(Area_growth_mm2 = Area_post_mm2 - Area_pre_mm2,
         Feret_growth_mm = Feret_post_mm - Feret_pre_mm)

ggplot(merged_df_cleaned_pre) +
  aes(x = Actual_shell_pre_mg, y = Area_pre_mm2) +
  geom_point(colour = "blue4") +
  theme_classic()

ggplot(merged_df_cleaned_all) +
  aes(x = Area_post_mm2, y = Area_pre_mm2) +
  geom_point(colour = "blue4") +
  theme_classic() #6 missing from pre, 23 missing from post, 29 missing overall with no overlap

ggplot(merged_df_cleaned_all) +
  aes(x = Feret_growth_mm, y = Area_growth_mm2) +
  geom_point(colour = "blue4") +
  theme_classic()

#Make a list of occurences the shell area growth is negative
negative_shellarea <- merged_df_cleaned_all %>%
  filter(Area_growth_mm2 < 0)
nrow(negative_shellarea) #89 is a lot, especially after excluding the ones with negative shell mass
View(negative_shellarea)

negative_shellarea %>%
  group_by(Phase_2_rep_R) %>%
  summarise(n = n())

hyp02<- negative_shellarea %>%
  filter(Phase_2_rep_R == "Hyp02")
View(hyp02)

both06<- negative_shellarea %>%
  filter(Phase_2_rep_R == "Both06")
View(both06)

#now positive
positive_shellarea <- merged_df_cleaned_all %>%
  filter(Area_growth_mm2 >= 0)

#Make a list of occurences the Feret growth is negative
negative_feret <- merged_df_cleaned_all %>%
  filter(Feret_growth_mm < 0)
nrow(negative_feret) #94 is a lot

#now positive
positive_feret <- merged_df_cleaned_all %>%
  filter(Feret_growth_mm >= 0)

#get the overlap between negative area and Feret sample names

# Get the vectors of Sample_Name from each dataframe
neg_shell_names <- negative_shellarea$Sample_Name
neg_feret_names <- negative_feret$Sample_Name

# Find the intersection and count overlap
overlap_names <- intersect(neg_shell_names, neg_feret_names)
length(overlap_names) #62 oysters lost shell area AND length out of 92


## run the Area Growth (mm^2) model on just the positive data
    #set contrasts ALWAYS RUN
    options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
    getOption("contrasts")
    
PAm2 <- lmer(Area_growth_mm2 ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
              (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
              (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive_shellarea, REML=TRUE)
Anova(PAm2, test="F", type="III")

#posthocs
emmeans(PAm2,specs = pairwise ~ Phase_1_DO, adjust = "none") #grew less area in hypoxia
emmeans(PAm2,specs = pairwise ~ Phase_2.1_DO, adjust = "none") #grew less area in hypoxia
emmeans(PAm2,specs = pairwise ~ Phase_2.1_temp*Phase_2.1_DO, adjust = "none") #warm, hypoxic, and both grew less than control

#diagnostics
leveneTest(Area_pre_mm2~ Phase_1_treat*Phase_2_treat, positive_shellarea) #passes
m1.e <- residuals(PAm2) #looks ok
qqnorm(m1.e)
qqline(m1.e)

View(negative_shellarea)

## run the Feret (mm) model on just the positive data
PFm2 <- lmer(Feret_growth_mm ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
               (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
               (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive_feret, REML=TRUE)
Anova(PFm2, test="F", type="III")

#posthocs
emmeans(PFm2,specs = pairwise ~ Phase_1_DO, adjust = "none") #grew less in hypoxia
emmeans(PFm2,specs = pairwise ~ Phase_2.1_DO, adjust = "none") # grew less in hypoxia
emmeans(PFm2,specs = pairwise ~ Phase_2.1_temp*Phase_2.1_DO, adjust = "none") #warm, hypoxic, and both grew less than control
  #not quite significant
emmeans(PFm2,specs = pairwise ~ Phase_1_DO*Phase_1_temp, adjust = "none") # P1warm grew more in P2 than P1 both, hypoxic, and control
emmeans(PFm2,specs = pairwise ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp, adjust = "none") 
  #P1 both re-exposed to p2 warming grew less than every other group in any P2 temperature (P1 cont exposed to warming camparison is not quite significant)
  #P1 warm re-exposed to p2 warming grew more than p1 cont and p1 hyp exposed to p2 warming

  #visualize this last post hoc to visually compare
  ss_positive_feret <- positive_feret %>%
  group_by(Phase_1_treat, Phase_2.1_temp) %>%
  mutate(
    mean_growth = mean(Feret_growth_mm, na.rm = TRUE),
    se_growth = std.error(Feret_growth_mm, na.rm = TRUE))

  ggplot(ss_positive_feret) +
  aes(x = Phase_1_treat, y = mean_growth) +
  geom_point(colour = "#112446") +
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), width = 0.2) +
  theme_classic() +
  facet_wrap(vars(Phase_2.1_temp))

#diagnostics
leveneTest(Area_pre_mm2~ Phase_1_treat*Phase_2_treat, merged_df_cleaned)
m1.e <- residuals(Am1) 
qqnorm(m1.e)
qqline(m1.e)

View(negative_shellarea)



#### Effects of Phase 1 ####
#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 

#copied from before, this is the cleaned pre dataset model will run on
merged_df_cleaned_pre <- mergedarea_df %>%
  filter(Exclude_pre_analysis != "Y" | is.na(Exclude_pre_analysis)) %>%
  mutate(Area_growth_mm2 = Area_post_mm2 - Area_pre_mm2,
         Feret_growth_mm = Feret_post_mm - Feret_pre_mm)

#run model on shell area size at the start of phase 2
#effects of phase 1
Am1 <- lmer(log(Area_pre_mm2) ~ Phase_1_DO*Phase_1_temp+ 
             (1|Phase_1_rep_R), data = merged_df_cleaned_pre, REML=TRUE)
Anova(Am1, test="F", type="III")

#posthocs
emmeans(Am1,specs = pairwise ~ Phase_1_DO, adjust = "none") #hyp grew less than norm
(231-223)/231

#diagnostics
leveneTest(log(Area_pre_mm2)~ Phase_1_treat*Phase_2_treat, merged_df_cleaned_pre) #passes
m1.e <- residuals(Am1) #ok
qqnorm(m1.e)
qqline(m1.e)

AIC(Am1)


#run model on feret
Fm1 <- lmer(log(Feret_pre_mm)~ Phase_1_DO*Phase_1_temp +
              (1|Phase_1_rep_R), data = merged_df_cleaned_pre, REML=TRUE)
Anova(Fm1, test="F", type="III")

#posthocs
emmeans(Fm1,specs = pairwise ~ Phase_1_DO, adjust = "none") # hyp grew less than norm

#diagnostics
leveneTest(log(Feret_pre_mm)~ Phase_1_treat*Phase_2_treat, merged_df_cleaned_pre) #passes
m1.e <- residuals(Fm1) #good
qqnorm(m1.e)
qqline(m1.e)


####FULL LMER MODEL####

#copied from before

#Edit dataset to exclude doubles and dead ones, calculate growth in area and Feret
merged_df_cleaned_all <- mergedarea_df %>%
  filter(Exclude_all != "Y" | is.na(Exclude_all)) %>%
  mutate(Area_growth_mm2 = Area_post_mm2 - Area_pre_mm2,
         Feret_growth_mm = Feret_post_mm - Feret_pre_mm) %>%
  filter( !is.na(Area_growth_mm2)) %>%
  mutate(Area_positivegrowth_mm2 = Area_growth_mm2 + 306)%>%
  mutate(Feret_positivegrowth_mm = Feret_growth_mm + 9)

min(merged_df_cleaned_all$Area_growth_mm2)
min(merged_df_cleaned_all$Feret_growth_mm)

## Area Growth (mm^2)
Am2 <- lmer(Area_growth_mm2 ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(Am2, test="F", type="III")

#with a log transformation + constant, the results are the same (slightly higher p-values) but pass the levene test
Am2_test <- lmer(log(Area_positivegrowth_mm2) ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
              (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
              (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(Am2_test, test="F", type="III")
#posthocs
emmeans(Am2_test,specs = pairwise ~ Phase_1_DO, adjust = "none") #grew less in hyp
emmeans(Am2_test,specs = pairwise ~ Phase_2.1_DO, adjust = "none") #grew less in hyp
emmeans(Am2_test,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO, adjust = "none") #grew less in both and hypoxic compared to control
#not significant
emmeans(Am2_test,specs = pairwise ~ Phase_2.1_temp*Phase_2.1_DO, adjust = "none") #grew more after warm than both or hypoxic
    #diagnostics
  leveneTest(log(Area_positivegrowth_mm2)~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #passes
  m1.e <- residuals(Am2_test) 
  qqnorm(m1.e)
  qqline(m1.e)

#posthocs
emmeans(Am2,specs = pairwise ~ Phase_1_DO, adjust = "none") #grew less in hyp
emmeans(Am2,specs = pairwise ~ Phase_2.1_DO, adjust = "none") #grew less in hyp
emmeans(Am2,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO, adjust = "none") #grew less in both and hypoxic compared to control
  #not significant
emmeans(Am2,specs = pairwise ~ Phase_2.1_temp*Phase_2.1_DO, adjust = "none") #grew more after warm than both or hypoxic
  
  #visualize phase 1 treatment effect
  ss_areagrowth <- merged_df_cleaned_all %>%
  group_by(Phase_1_treat) %>%
    summarise(
    mean_growth = mean(Area_growth_mm2, na.rm = TRUE),
    se_growth = std.error(Area_growth_mm2, na.rm = TRUE))
nrow(ss_areagrowth)

ss_areagrowth$Phase_1_treat <- factor(ss_areagrowth$Phase_1_treat, 
                                            levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(ss_areagrowth) +
  aes(x = Phase_1_treat, y = mean_growth) +
  geom_point(colour = "#112446") +
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), width = 0.2) +
  theme_classic()

#visualize phase 2 treatment effect
ss_areagrowth <- merged_df_cleaned_all %>%
  group_by(Phase_2_treat) %>%
  summarise(
    mean_growth = mean(Area_growth_mm2, na.rm = TRUE),
    se_growth = std.error(Area_growth_mm2, na.rm = TRUE))
nrow(ss_areagrowth)

ss_areagrowth$Phase_2_treat <- factor(ss_areagrowth$Phase_2_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(ss_areagrowth) +
  aes(x = Phase_2_treat, y = mean_growth) +
  geom_point(colour = "#112446") +
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), width = 0.2) +
  theme_classic()

#visualize 3 way interaction
ss_areagrowth <- merged_df_cleaned_all %>%
  group_by(Phase_1_temp, Phase_2_treat) %>%
  summarise(
    mean_growth = mean(Area_growth_mm2, na.rm = TRUE),
    se_growth = std.error(Area_growth_mm2, na.rm = TRUE))
nrow(ss_areagrowth)
levels(ss_areagrowth$Phase_1_temp)

ss_areagrowth$Phase_1_temp <- factor(ss_areagrowth$Phase_1_temp, 
                                      levels = c("Ambient", "Warm"))
ss_areagrowth$Phase_2_treat <- factor(ss_areagrowth$Phase_2_treat, 
                                      levels = c("Cont", "Warm","Hyp", "Both"))

ggplot(ss_areagrowth) +
  aes(x = Phase_2_treat, y = mean_growth) +
  geom_point(colour = "#112446") +
  labs(x = "Phase 2 Treatment",
       y = "Shell Area Growth (mm^2)",
       title = "Phase 1 Temperature")+
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), width = 0.2) +
  theme_classic()+
  theme(
    plot.title = element_text(size = 18L,
                              hjust = 0.5),
    axis.title.y = element_text(size = 18L),
    axis.title.x = element_text(size = 18L)
  )+
  facet_wrap(vars(Phase_1_temp))+ scale_fill_hue()


ggplot(Alkalinity_data) +
  aes(x = treatment, y = alkalinity_umolL) +
  geom_point(colour = "#112446") +
  labs(
    x = "Phase 2 Treatment",
    y = "Shell Area Growth (mm^2)",
    title = "Phase 1 Temperature"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(size = 18L,
                              hjust = 0.5),
    axis.title.y = element_text(size = 15L),
    axis.title.x = element_text(size = 18L),
    axis.text.x = element_text(size = 12L)
  )

#diagnostics
leveneTest(Area_growth_mm2~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #does not pass
m1.e <- residuals(Am2) 
qqnorm(m1.e)
qqline(m1.e)


## Feret Growth (mm)
Fm2 <- lmer(Feret_growth_mm ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
              (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
              (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(Fm2, test="F", type="III")

#posthocs
emmeans(Fm2,specs = pairwise ~ Phase_1_DO, adjust = "none") #hyp grew less
emmeans(Fm2,specs = pairwise ~ Phase_2.1_DO, adjust = "none") #hyp grew less
emmeans(Fm2,specs = pairwise ~ Phase_1_DO*Phase_1_temp, adjust = "none") #p1 both grew less than p1 warm and p1 control
  #not quite significant
emmeans(Fm2,specs = pairwise ~ Phase_2.1_temp*Phase_2.1_DO, adjust = "none")

#diagnostics
leveneTest(Feret_growth_mm~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #doesn't pass
m1.e <- residuals(Fm2) #ok
qqnorm(m1.e)
qqline(m1.e)

## LOG Feret Growth (mm)
Fm2_test <- lmer(log(Feret_positivegrowth_mm) ~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
              (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
              (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(Fm2_test, test="F", type="III")

#posthocs
emmeans(Fm2_test,specs = pairwise ~ , adjust = "none")

#diagnostics
leveneTest(log(Feret_positivegrowth_mm)~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #doesn't pass
m1.e <- residuals(Fm2_test) #ok
qqnorm(m1.e)
qqline(m1.e)




options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups sum
getOption("contrasts") 

#### AREA AS A COVARIATE TO MASS ANALYSIS ####

##tissue growth (mg)
m1_co_area <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Area_growth_mm2+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(m1_co_area, test="F", type="III")

#posthocs
emmeans(m1_co_area,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp, adjust = "none") 
  #not quite significant, oysters exposed to early life warming don't do as well as oysters in early life ambient water when reexposed to warming

#diagnostics
leveneTest(Actual_tissue_growth_mg~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #pass
m1.e <- residuals(m1_co_area) #ok, tails trail off
qqnorm(m1.e)
qqline(m1.e)


#now doing Feret as the tissue covariate

m1_co_feret <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Feret_growth_mm+
                     (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                     (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(m1_co_feret, test="F", type="III")

#better, p value for Phase_1_temp:Phase_2.1_temp is 0.07 instead of 0.09

#posthocs
emmeans(m1_co_feret,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp, adjust = "none") 
#not quite significant, oysters exposed to early life warming don't do as well as oysters in ambient water when exposed to warming

#diagnostics
leveneTest(Actual_tissue_growth_mg~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #pass
m1.e <- residuals(m1_co_feret) #pretty exponential, would probably need a log transformation
qqnorm(m1.e)
qqline(m1.e)

AIC(m1_co_area, m1_co_feret) #m1_co_area better than both the others
AIC(m1)




##Shell growth (mg)
m2_co_area <- lmer(log(Actual_shell_growth_mg)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Area_growth_mm2+
                     (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                     (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(m2_co_area, test="F", type="III")

#posthocs
emmeans(m2_co_area,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO, adjust = "none") 
#oysters grow less shell in hypoxia when they have previously been exposed to early life warming
#oysters exposed to early life warming grow more in p2 normoxia than hypoxia

#diagnostics
leveneTest(log(Actual_shell_growth_mg)~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #pass
m1.e <- residuals(m2_co_area) #ok, tails trail off
qqnorm(m1.e)
qqline(m1.e)


#now doing Feret as the tissue covariate

m2_co_feret <- lmer(Actual_shell_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Feret_growth_mm+
                      (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                      (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(m2_co_feret, test="F", type="III")

#better, p value for Phase_1_temp:Phase_2.1_temp is 0.07 instead of 0.09

#posthocs
emmeans(m2_co_feret,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO, adjust = "none") 
#not quite significant, oysters exposed to early life warming grow more shell in p2 normoxia than in p2 hypoxia

#diagnostics
leveneTest(Actual_tissue_growth_mg~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #pass
m1.e <- residuals(m2_co_feret) #pretty exponential, would probably need a log transformation
qqnorm(m1.e)
qqline(m1.e)




##Tissue: Shell growth (mg)
m3_co_area <- lmer((Actual_tissue_growth_mg/Actual_shell_growth_mg)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Area_growth_mm2+
                     (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                     (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(m3_co_area, test="F", type="III")

#posthocs
emmeans(m3_co_area,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO, adjust = "none") 
  #oysters exposed to early life warming don't grow as mcuh in normoxia later on

  #almost significant SKM
emmeans(m3_co_area,specs = pairwise ~ Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO, adjust = "none")
  #visualize
ss <- merged_df_cleaned_all %>%
  group_by(Phase_1_temp, Phase_2_treat) %>%
  summarise(
    mean_growth = mean(Actual_tissue_growth_mg / Actual_shell_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_mg / Actual_shell_growth_mg, na.rm = TRUE),
    .groups = "drop")

ggplot(ss) +
  aes(x = Phase_1_temp, y = mean_growth, color = Phase_2_treat) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), width = 0.2) +
  theme_classic() +
  facet_wrap(vars(Phase_2_treat))

#diagnostics
leveneTest((Actual_tissue_growth_mg/Actual_shell_growth_mg)~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #pass
m1.e <- residuals(m3_co_area) #ok, tails trail off
qqnorm(m1.e)
qqline(m1.e)


#now doing Feret as the tissue covariate

m3_co_feret <- lmer((Actual_tissue_growth_mg/Actual_shell_growth_mg)~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Feret_growth_mm+
                      (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
                      (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = merged_df_cleaned_all, REML=TRUE)
Anova(m3_co_feret, test="F", type="III")

#posthocs
emmeans(m3_co_feret,specs = pairwise ~ Phase_1_temp*Phase_2.1_DO, adjust = "none")
#oysters exposed to early life warming rather than ambient, grow less relative tissue in p2 normoxia

#diagnostics
leveneTest((Actual_tissue_growth_mg/Actual_shell_growth_mg)~Phase_1_treat*Phase_2_treat, merged_df_cleaned_all) #pass
m1.e <- residuals(m3_co_feret) #tails trail off
qqnorm(m1.e)
qqline(m1.e)
