#load packages
install.packages("plotrix")
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
library(plotrix) #for standard error

#growth and N content
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
levels(Growth_Data_forR$Phase_1_DO)
levels(Growth_Data_forR$Phase_1_temp)
levels(Growth_Data_forR$Phase_2.1_DO)
levels(Growth_Data_forR$Phase_2.1_temp)
levels(Growth_Data_forR$Phase_1_treat)
levels(Growth_Data_forR$Phase_2_treat)


#set contrasts ALWAYS RUN
options(contrasts = c("contr.sum","contr.poly")) #could also be contr.treatment for unequal groups
getOption("contrasts") 

#filter out rows with NA in Actual_tissue_growth_mg only to not include data from dead oysters in analysis
Growth_Data_forR_clean <- Growth_Data_forR[!is.na(Growth_Data_forR$Actual_tissue_growth_mg), ]
View(Growth_Data_forR_clean)
summary(Growth_Data_forR_clean) #check data

# Filter to remove NA values in Actual_tissue_growth_mg and see if results change with middle size range
mean(Growth_Data_forR_clean$Ratio_tissue_shell_mg) #mean
sd(Growth_Data_forR_clean$Ratio_tissue_shell_mg)/sqrt(length((Growth_Data_forR_clean$Ratio_tissue_shell_mg))) #se
sd(Growth_Data_forR_clean$Ratio_tissue_shell_mg)

Growth_Data_forR_test <- Growth_Data_forR %>%
  filter(!is.na(Actual_tissue_growth_mg) & Ratio_tissue_shell_mg < 1.080456 & Ratio_tissue_shell_mg > -1.01848)
View (Growth_Data_forR_test)

#mean and standard error calculations
0.03098755+1.049468
0.03098755-1.049468

#start model selection pg. 75 in Zuur et al, not used for actual running model
#tissue data
#no random effects included yet, this runs
Test1 <- lm(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data=Growth_Data_forR_clean)
anova(Test1)

#trying to figure out proportional variances to treatments, if variances are homogeneous in treatments, not working
vf1Fixed <- varFixed(~Actual_shell_pre_mg) ## have also tried this being Phase1_Phase2_treat
Test1 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, weights = vf1Fixed, data=Growth_Data_forR_clean)
plot(Test1, which = c(1), col = 1, add.smooth = FALSE, caption = "")

vf2Fixed <- varIdent(form = ~1 | Phase1_Phase2_treat)
Test1 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, weights = vf2Fixed, data=Growth_Data_forR_clean)

vf1DO=varIdent(form=~1|Phase_1_DO)#1
vf1T=varIdent(form=~1|Phase_1_temp)#2
vf2DO=varIdent(form=~1|Phase_2.1_DO)#3
vf2T=varIdent(form=~1|Phase_2.1_temp)#4
vf1DO1T=varIdent(form=~1|Phase_1_DO*Phase_1_temp)#5
vf1DO2DO=varIdent(form=~1|Phase_1_DO*Phase_2.1_DO)#6
vf1T2DO=varIdent(form=~1|Phase_1_temp*Phase_2.1_DO)#7
vf1DO2T=varIdent(form=~1|Phase_1_DO*Phase_2.1_temp)#8
vf1T2T=varIdent(form=~1|Phase_1_temp*Phase_2.1_temp)#9
vf2DO2T=varIdent(form=~1|Phase_2.1_DO*Phase_2.1_temp)#10
vf1DO1T2DO=varIdent(form=~1|Phase_1_DO*Phase_1_temp*Phase_2.1_DO)#11
vf1DO1T2T=varIdent(form=~1|Phase_1_DO*Phase_1_temp*Phase_2.1_temp)#12
vf1DO2DO2T=varIdent(form=~1|Phase_1_DO*Phase_2.1_DO*Phase_2.1_temp)#13
vf1T2DO2T=varIdent(form=~1|Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp)#14
vf4=varIdent(form=~1|Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp)#15

  #tissue weights model, 4 way interaction is the best weight for model
gls1 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO)
gls2 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T)
gls3 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2DO)
gls4 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2T)
gls5 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO1T)
gls6 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO2DO)
gls7 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T2DO)
gls8 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO2T)
gls9 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T2T)
gls10 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2DO2T)
gls11 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO1T2DO)
gls12 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO1T2T)
gls13 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO2DO2T)
gls14 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T2DO2T)
gls15 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf4)
gls16 <- gls(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2Fixed)
  
  #shell weights model, 4 way interaction is the best weight for model
gls1 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO)
gls2 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T)
gls3 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2DO)
gls4 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2T)
gls5 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO1T)
gls6 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO2DO)
gls7 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T2DO)
gls8 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO2T)
gls9 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T2T)
gls10 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2DO2T)
gls11 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO1T2DO)
gls12 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO1T2T)
gls13 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO2DO2T)
gls14 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T2DO2T)
gls15 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf4)
gls16 <- gls(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2Fixed)

#tissue:shell weights model, 4 way interaction is the best weight for model
gls1 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO)
gls2 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T)
gls3 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2DO)
gls4 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2T)
gls5 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO1T)
gls6 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO2DO)
gls7 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T2DO)
gls8 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO2T)
gls9 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T2T)
gls10 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2DO2T)
gls11 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO1T2DO)
gls12 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO1T2T)
gls13 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1DO2DO2T)
gls14 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf1T2DO2T)
gls15 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf4)
gls16 <- gls(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, data = Growth_Data_forR_clean, weights = vf2Fixed)

AIC(gls1,gls2,gls3,gls4,gls5,gls6,gls7,gls8,gls9,gls10,gls11,gls12,gls13,gls14,gls15) #lowest AIC score is best model
AIC(gls16)


#working model with variances and random effects
#save variance to use in model, model did not pass lavene test without this
vf2Fixed <- varIdent(form = ~1 | Phase1_Phase2_treat)
vf_1 <- varIdent(form = ~1 | Phase_1_treat) #test running separately
vf_2 <- varIdent(form = ~1 | Phase_2_treat) #test running separately

a <- varIdent(form = ~1 | Phase_1_DO)
Mlmea <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp+ Actual_shell_pre_mg,
             method = "REML", 
             #is there a way to make the fixed effects nested? use ML instead of REML if yes pg. 122
             weights = a, 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = Growth_Data_forR_clean)
b <- varIdent(form = ~1 | Phase_1_temp)
Mlmeb <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp+ Actual_shell_pre_mg,
             method = "REML", 
             weights = b, 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = Growth_Data_forR_clean)
c <- varIdent(form = ~1 | Phase_2.1_DO)

AIC(Mlmea, Mlmeb)

#tissue growth models

  #diagnostics
leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, Growth_Data_forR_clean) #passes but still added varIdent for unequal variances
m1.e <- residuals(Mlme1) 
qqnorm(m1.e) #fairly normal, anova is robust to non-normality

  #model using lme
Mlme1 <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp+ Actual_shell_pre_mg,
             method = "REML", 
             #is there a way to make the fixed effects nested? use ML instead of REML if yes pg. 122
             weights = varIdent(vf2Fixed), 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = Growth_Data_forR_clean)
anova(Mlme1, type = "m")
  #testing AIC of model to compare to others
AIC(Mlme1)

  #posthocs
emm1 <- emmeans(Mlme1,specs = pairwise ~ Phase_2.1_DO, adjust = "none") 
emm1$emmeans 
emm1$contrasts

  #just phase 1, TEST
MlmeA <- lme(Actual_tissue_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             weights = varIdent(vf_1), 
             random = ~1 | Phase_1_rep_R, data = Growth_Data_forR_clean)
anova(MlmeA, type = "m")
  #just phase 2, TEST
MlmeB <- lme(Actual_tissue_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp+ Actual_shell_pre_mg,
             method = "REML", 
             weights = varIdent(vf_2), 
             random = ~1 | Phase_2_rep_R, data = Growth_Data_forR_clean)
anova(MlmeB, type = "m")


#shell growth models

  #diagnostics
leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, Growth_Data_forR_clean) #doesn't pass so added varIdent for unequal variances
m2.e <- residuals(Mlme2) 
qqnorm(m2.e) #looks good

  #model using lme
Mlme2 <- lme(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg, 
             method = "REML", 
             weights = varIdent(vf2Fixed), 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = Growth_Data_forR_clean)
anova(Mlme2, type = "m")

#posthocs
emm2 <- emmeans(Mlme2,specs = pairwise ~ Phase_2.1_DO, adjust = "none") 
emm2$emmeans
emm2$contrasts

  #just phase 1
MlmeC <- lme(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg,
             method = "REML", 
             weights = varIdent(vf_1), 
             random = ~1 | Phase_1_rep_R, data = Growth_Data_forR_clean)
anova(MlmeC, type = "m")
  #just phase 2
MlmeD <- lme(Actual_shell_growth_mg ~ Phase_2.1_DO * Phase_2.1_temp+ Actual_shell_pre_mg,
             method = "REML", 
             weights = varIdent(vf_2), 
             random = ~1 | Phase_2_rep_R, data = Growth_Data_forR_clean)
anova(MlmeD, type = "m")


#Tissue:shell growth

  #diagnostics
leveneTest(Ratio_tissue_shell_mg~Phase1_Phase2_treat, Growth_Data_forR_clean) #passes but still added varIdent for unequal variances
m3.e <- residuals(Mlme3) 
qqnorm(m3.e) #looks good in the middle, trails off towards the higher end

  #model using lme
Mlme3 <- lme(Ratio_tissue_shell_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp + Actual_shell_pre_mg, 
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = Growth_Data_forR_test)
anova(Mlme3, type = "m")

AIC(Mlme3)

#test model to compare Phase_1_treat and Phase_2_treat effects
Mlme4 <- lme(Ratio_tissue_shell_mg ~ Phase_1_treat * Phase_2_treat + Actual_shell_pre_mg, 
             method = "REML", 
             random = list(~1 | Phase_2_rep_R, ~1 | Phase_1_rep_R, ~1 | Phase1_Phase2_rep), data = Growth_Data_forR_clean)
anova(Mlme4, type = "m")


#Doing analysis on initial tissue and shell mass (not growth) at start of phase 2 -- effects of phase 1 (no phase 2 in model)
  #tissue
leveneTest(Actual_tissue_pre_mg~Phase_1_treat, Growth_Data_forR_clean) #passes but still added varIdent for unequal variances
m3.e <- residuals(Mlme3) 
qqnorm(m3.e)

Mlme_i1 <- lme( Actual_tissue_pre_mg ~ Phase_1_DO * Phase_1_temp + Actual_shell_pre_mg, 
             method = "REML", 
             weights = varIdent(vf_1), 
             random = list( ~1 | Phase_1_rep_R), data = Growth_Data_forR_clean)
anova(Mlme_i1, type = "m")
emmi1 <- emmeans(Mlme_i1, specs = pairwise ~ Phase_1_DO*Phase_1_temp, adjust="none")
emmi1$emmeans
emmi1$contrasts 

  #shell
Mlme_i2 <- lme( Actual_shell_pre_mg ~ Phase_1_DO * Phase_1_temp + Actual_tissue_pre_mg, 
                method = "REML", 
                weights = varIdent(vf_1), 
                random = list( ~1 | Phase_1_rep_R), data = Growth_Data_forR_clean)
anova(Mlme_i2, type = "m")

emmi2 <- emmeans(Mlme_i2, specs = pairwise ~ Phase_1_DO*Phase_1_temp, adjust="none")
emmi2$emmeans
emmi2$contrasts 

emmi3 <- emmeans(Mlme_i2, specs = pairwise ~ Phase_1_DO, adjust="none")
emmi3$emmeans
emmi3$contrasts 
  #t:s
Mlme_i3 <- lme( (Actual_tissue_pre_mg/Actual_shell_pre_mg) ~ Phase_1_DO * Phase_1_temp + Actual_tissue_pre_mg, 
                method = "REML", 
                weights = varIdent(vf_1), 
                random = list( ~1 | Phase_1_rep_R), data = Growth_Data_forR_clean)
anova(Mlme_i3, type = "m")

emmi4 <- emmeans(Mlme_i3, specs = pairwise ~ Phase_1_DO, adjust="none")
emmi4$emmeans
emmi4$contrasts 

emmi5 <- emmeans(Mlme_i3, specs = pairwise ~ Phase_1_temp, adjust="none")
emmi5$emmeans
emmi5$contrasts 

#check variance of models
plot(Mlme_i3, which = c(1), col = 1, add.smooth = FALSE, caption = "")

#shell growth
vf2Fixed <- varIdent(form = ~1 | Phase1_Phase2_treat)
Mlme2 <- lme(Actual_shell_growth_mg ~ Phase_1_DO * Phase_1_temp * Phase_2.1_DO * Phase_2.1_temp, method = "REML",
             weights = varIdent(Actual_shell_growth_mg), #this model works because the variances are numeric
             random = ~1 | Phase_2_rep_R / Phase_1_rep_R, data = Growth_Data_forR_clean)
anova(Mlme2)
AIC(Mlme2)

#check variance of models
plot(Mlme2, which = c(1), col = 1, add.smooth = FALSE, caption = "")

####LMER MODELS####

##tissue growth (mg)
m1 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_temp*Phase_2.1_DO+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_clean, REML=TRUE)
Anova(m1, test="F", type="III")

#testing
m1 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase1_Phase2_rep), data = Growth_Data_forR_clean, REML=TRUE)
Anova(m1, test="F", type="III")

#posthocs
emm1 <- emmeans(m1,specs = pairwise ~ Phase_2.1_DO, adjust = "none") 
emm1$emmeans 
emm1$contrasts

#diagnostics
leveneTest(Actual_tissue_growth_mg~Phase1_Phase2_treat, Growth_Data_forR_clean) 
m1.e <- residuals(m1) 
qqnorm(m1.e) #not quite normal but ancova is robust to non-normality


##shell growth (mg)
m2 <- lmer(Actual_shell_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_clean, REML=TRUE)
oneway.test(Actual_tissue_growth_mg ~ Phase1_Phase2_treat, data = Growth_Data_forR_clean, var.equal = FALSE)
Anova(m2, test="F", type="III")

#posthocs
emm2 <- emmeans(m2,specs = pairwise ~ Phase_2.1_DO, adjust = "none") 
emm2$emmeans #More shell growth in normoxic conditions than in hypoxic conditions
emm2$contrasts

#diagnostics
leveneTest(Actual_shell_growth_mg~Phase1_Phase2_treat, Growth_Data_forR_clean) #failed
oneway.test(Actual_shell_growth_mg~Phase1_Phase2_treat, data = Growth_Data_forR_clean, var.equal = FALSE)
m2.e <- residuals(m2) 
qqnorm(m2.e) #not quite normal but ancova is robust to non-normality

##tissue:shell growth 
m3 <- lmer(Ratio_tissue_shell_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR_clean, REML=TRUE)
Anova(m3, test="F", type="III")

#posthocs
emm3 <- emmeans(m3,specs = pairwise ~ Phase_1_DO*Phase_2.1_DO*Phase_2.1_temp , adjust="none")
emm3$emmeans
emm3$contrasts #none are significant

emm3 <- emmeans(m3,specs = pairwise ~ Phase_1_DO, adjust="none")
emm3$emmeans 
emm3$contrasts #no significant difference between Hyp and Norm

#diagnostics
leveneTest(Ratio_tissue_shell_mg~Phase1_Phase2_treat, Growth_Data_forR_clean) 
m3.e <- residuals(m3) 
qqnorm(m3.e) #not quite normal but ancova is robust to non-normality


#try other adjustments Dunnett’s test is used when comparing each treatment level to a control group
emm3 <- emmeans(m3,specs = pairwise ~ Phase_1_DO, adjust="dunnett", Phase_1_treat = "Cont")
emm3$emmeans
emm3$contrasts #literally doesn't change anything

#Sidák Adjustment (for fewer comparisons)

emm3 <- emmeans(m3,specs = pairwise ~ Phase_1_DO, adjust="sidak")
emm3$emmeans
emm3$contrasts #literally doesn't change anything



####SIZE CLASSES####
min(Growth_Data_forR_clean$Actual_tissue_pre_mg)
max(Growth_Data_forR_clean$Actual_tissue_pre_mg)

positive1 <- Growth_Data_forR_clean[Growth_Data_forR_clean$Actual_tissue_pre_mg >= 0 & 
                                     Growth_Data_forR_clean$Actual_tissue_pre_mg <= 100, ]

nrow(positive1)

m4 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive1, REML=TRUE)
Anova(m4, test="F", type="III")


positive2 <- Growth_Data_forR_clean[Growth_Data_forR_clean$Actual_tissue_pre_mg >= 100 & 
                                     Growth_Data_forR_clean$Actual_tissue_pre_mg <= 200, ]
nrow(positive2)

m5 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive2, REML=TRUE)
Anova(m5, test="F", type="III")

positive3 <- Growth_Data_forR_clean[Growth_Data_forR_clean$Actual_tissue_pre_mg >= 200 & 
                                      Growth_Data_forR_clean$Actual_tissue_pre_mg <= 300, ]
nrow(positive3)

m6 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive3, REML=TRUE)
Anova(m6, test="F", type="III")


positive4 <- Growth_Data_forR_clean[Growth_Data_forR_clean$Actual_tissue_pre_mg >= 300 & 
                                      Growth_Data_forR_clean$Actual_tissue_pre_mg <= 400, ]
nrow(positive4)

m7 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive4, REML=TRUE)
Anova(m7, test="F", type="III")


positive5 <- Growth_Data_forR_clean[Growth_Data_forR_clean$Actual_tissue_pre_mg >= 400 & 
                                      Growth_Data_forR_clean$Actual_tissue_pre_mg <= 500, ]
nrow(positive5)

m8 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive5, REML=TRUE)
Anova(m8, test="F", type="III")



positive6 <- Growth_Data_forR_clean[Growth_Data_forR_clean$Actual_tissue_pre_mg >= 700 & 
                                      Growth_Data_forR_clean$Actual_tissue_pre_mg <= 800, ]
nrow(positive6)

m9 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep_R)+(1|Phase_1_rep_R)+
             (1|Phase_2_rep_R:Phase_1_DO)+(1|Phase_2_rep_R:Phase_1_temp)+(1|Phase_2_rep_R:Phase_1_DO:Phase_1_temp), data = positive6, REML=TRUE)
Anova(m9, test="F", type="III")

#diagnostics
leveneTest(Ratio_tissue_shell_mg~Phase1_Phase2_treat, positive1) 
m5.e <- residuals(m5) 
qqnorm(m3.e) #not quite normal but ancova is robust to non-normality

#post hoc
emm5 <- emmeans(m5,specs = pairwise ~ Phase_1_DO , adjust="none")
emm3$emmeans
emm3$contrasts #none are significant

positive2 <- Growth_Data_forR_clean[Growth_Data_forR_clean$Actual_tissue_growth_mg >= 200 & 
                                      Growth_Data_forR_clean$Actual_tissue_growth_mg <= 400, ]

m6 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep)+(1|Phase_1_rep)+
             (1|Phase_2_rep:Phase_1_DO)+(1|Phase_2_rep:Phase_1_temp)+(1|Phase_2_rep:Phase_1_DO:Phase_1_temp), data = positive2, REML=TRUE)
Anova(m6, test="F", type="III")

positive3 <- Growth_Data_forR_clean[Growth_Data_forR_clean$Actual_tissue_growth_mg >= 400 & 
                                      Growth_Data_forR_clean$Actual_tissue_growth_mg <= 600, ]

m7 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep)+(1|Phase_1_rep)+
             (1|Phase_2_rep:Phase_1_DO)+(1|Phase_2_rep:Phase_1_temp)+(1|Phase_2_rep:Phase_1_DO:Phase_1_temp), data = positive3, REML=TRUE)
Anova(m7, test="F", type="III") ##There aren't enough oysters





second_m1 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep)+(1|Phase_1_rep)+
             (1|Phase_2_rep:Phase_1_DO)+(1|Phase_2_rep:Phase_1_temp)+(1|Phase_2_rep:Phase_1_DO:Phase_1_temp), data = Growth_Data_clean, REML=TRUE)
Anova(second_m1, test="F", type="III")

 second_m2 <- lmer(Actual_shell_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep)+(1|Phase_1_rep)+
             (1|Phase_2_rep:Phase_1_DO)+(1|Phase_2_rep:Phase_1_temp)+(1|Phase_2_rep:Phase_1_DO:Phase_1_temp), data = Growth_Data_clean, REML=TRUE)
Anova(second_m2, test="F", type="III")

second_m3 <- lmer(Ratio_tissue_shell_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep)+(1|Phase_1_rep)+
             (1|Phase_2_rep:Phase_1_DO)+(1|Phase_2_rep:Phase_1_temp)+(1|Phase_2_rep:Phase_1_DO:Phase_1_temp), data = Growth_Data_clean, REML=TRUE)
Anova(second_m3, test="F", type="III")



##tissue growth (mg)
m1 <- lmer(Actual_tissue_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep)+(1|Phase_1_rep)+
             (1|Phase_2_rep:Phase_1_DO)+(1|Phase_2_rep:Phase_1_temp)+(1|Phase_2_rep:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR, REML=TRUE)
Anova(m1, test="F", type="III")

#diagnostics
leveneTest(TissueGrowthmg~OverallTreatmentCombination, Growth_Data_forR) 
m1.e <- residuals(m1) 
qqnorm(m1.e) #not quite normal but ancova is robust to non-normality

#posthocs
emm1 <- emmeans(m1,specs = pairwise ~ Phase_1_temp, adjust = "none") 
emm1$emmeans
emm1$contrasts


##shell growth (mg)
m2 <- lmer(Actual_shell_growth_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep)+(1|Phase_1_rep)+
             (1|Phase_2_rep:Phase_1_DO)+(1|Phase_2_rep:Phase_1_temp)+(1|Phase_2_rep:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR, REML=TRUE)
Anova(m2, test="F", type="III")

#diagnostics
leveneTest(ShellGrowthmg~OverallTreatmentCombination, Growth_Data_forR) 
m2.e <- residuals(m2)
qqnorm(m2.e) #not quite normal but ancova is robust to non-normality

#posthocs
emm2 <- emmeans(second_m2,specs = pairwise ~ Phase_2.1_DO, adjust = "none") 
emm2$emmeans
emm2$contrasts


##tissue:shell growth 
m3 <- lmer(na.omit(Ratio_tissue_shell_mg)~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep)+(1|Phase_1_rep)+
             (1|Phase_2_rep:Phase_1_DO)+(1|Phase_2_rep:Phase_1_temp)+(1|Phase_2_rep:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR, REML=TRUE)
Anova(m3, test="F", type="III")

#diagnostics
leveneTest(TissueShellGrowth~OverallTreatmentCombination, Growth_Data_forR) 
m3.e <- residuals(m3) 
qqnorm(m3.e) #not quite normal but ancova is robust to non-normality





#create mean and standard error plots

# Calculate means and standard deviations for each Phase_1_treat and Phase_2_treat group
  #tissue growth
summary_stats_t <- Growth_Data_forR_clean %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  summarise(mean_growth = mean(Actual_tissue_growth_mg, na.rm = TRUE),
    se_growth = std.error(Actual_tissue_growth_mg, na.rm = TRUE))
View(summary_stats_t)

  #shell growth
summary_stats_s <- Growth_Data_forR_clean %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  summarise(mean_growth = mean(Actual_shell_growth_mg, na.rm = TRUE),
            se_growth = std.error(Actual_shell_growth_mg, na.rm = TRUE))

#tissue:shell growth
summary_stats_t_s <- Growth_Data_forR_clean %>%
  group_by(Phase_1_treat, Phase_2_treat) %>%
  summarise(mean_growth = mean(Ratio_tissue_shell_mg, na.rm = TRUE),
            se_growth = std.error(Ratio_tissue_shell_mg, na.rm = TRUE))



  #11.19.24
  #plot the figure with points
# Reorder Phase_1_treat and Phase_2_treat
summary_stats_t$Phase_1_treat <- factor(summary_stats_t$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))
summary_stats_t$Phase_2_treat <- factor(summary_stats_t$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp", "Both"))
summary_stats_s$Phase_1_treat <- factor(summary_stats_s$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp",  "Both"))
summary_stats_s$Phase_2_treat <- factor(summary_stats_s$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp",  "Both"))
summary_stats_t_s$Phase_1_treat <- factor(summary_stats_t_s$Phase_1_treat, 
                                        levels = c("Cont", "Warm","Hyp",  "Both"))
summary_stats_t_s$Phase_2_treat <- factor(summary_stats_t_s$Phase_2_treat, 
                                        levels = c("Cont", "Warm","Hyp",  "Both"))

#Tissue plot with mean and SD
ggplot(summary_stats_t, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
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
ggplot(summary_stats_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
  geom_point(size = 4, position = position_dodge(0.9)) + # Plot means as points
  geom_errorbar(aes(ymin = mean_growth - se_growth, ymax = mean_growth + se_growth), 
                width = 0.2, position = position_dodge(0.9)) + # Error bars for SD
  theme_classic() +
  guides(color = "none") + # Remove legend for color
  facet_wrap(vars(Phase_2_treat), scales = "fixed", nrow = 1) + # Facet by Phase_2_treat
  scale_color_brewer(palette = "Set2") + # Use color palette for points
  labs(x = "Phase 1 Treatment", y = "Mean Shell Growth (mg)") +
  theme(legend.position = "none") # Remove legend

#Tissue:Shell plot with mean and SD
ggplot(summary_stats_t_s, aes(x = Phase_1_treat, y = mean_growth, color = Phase_1_treat)) +
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




##for formatting tables into overleaf
install.packages("stargazer")
library(stargazer)
stargazer(m1, type = "latex")

#normal table
install.packages("flextable")
library(flextable)
flextable(Warm1)

