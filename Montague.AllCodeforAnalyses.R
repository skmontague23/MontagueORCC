#code adapted from Donelan et al.
#Paper Title

#import data and load packages
library(lme4)
library(readr)
library (car)
library(emmeans)
library(boot)

#water quality
Water.Qual <- read_csv("Desktop/Manuscripts 2019-/Oyster Carryover19/Data for Figshare/1.WaterQuality.csv", col_types = cols(Replicate = col_factor(), DO = col_factor(), Temperature=col_factor(), CyclePhase=col_factor(), Year=col_factor(), StartEndHypoxicCycle=col_factor()))
str(Water.Qual)

#Alkalinity and pCO2
Alk <- read_csv("Desktop/Manuscripts 2019-/Oyster Carryover19/Data for Figshare/2.AlkalinitypCO2.csv", col_types = cols(Replicate = col_factor(), Temperature=col_factor(), Year=col_factor()))
str(Alk)


#### SKM start here ####
#growth and N content
Growth_Data_forR <- read_csv("Desktop/MontagueORCC/Oyster Weight Data/growth_phase2.1_weightsSKM.csv", 
                             col_types = cols(Phase_1_temp = col_factor(), 
                                              Phase_1_DO =col_factor(), 
                                              Phase_2.1_temp =col_factor(), 
                                              Phase_2.1_DO=col_factor(), 
                                              Phase_1_rep =col_factor(), 
                                              Phase_2_rep =col_factor(),
                                              Ratio_tissue_shell_mg = col_double()))
str(Growth_Data_forR)
levels(Growth_Data_forR$Phase_1_DO)
levels(Growth_Data_forR$Phase_1_temp)
levels(Growth_Data_forR$Phase_2.1_DO)
levels(Growth_Data_forR$Phase_2.1_temp)


#allometric equations
allo.equations <- read_csv("~/Desktop/Manuscripts 2019-/**SUBMITTED/Oyster Carryover19/Data for Figshare/4.AllometricEquationsData.csv", col_types = cols(OverallTreatmentCombination = col_factor()))

#shell lengths for bootstrapping, all oysters
bootstrap.allo <- read_csv("~/Desktop/Manuscripts 2019-/**SUBMITTED/Oyster Carryover19/Data for Figshare/5.ReefShellLengths.csv")
bootstrap.allo <- as.data.frame(bootstrap.allo)


#water quality
#LMMs and ANCOVAs for response variables
#set contrasts
options(contrasts = c("contr.sum","contr.poly"))
getOption("contrasts") 

#DO during hypoxic phase of cycle, Year 1
w1 <- lmer(DissolvedOxygenmgL ~ DO*Temperature+(1|Replicate), data = subset(Water.Qual, CyclePhase=="Hypoxia" & Year=="1"), REML=TRUE)
Anova(w1, test="F", type="III")

#DO during hypoxic phase of cycle, Year 2
w2 <- lmer(DissolvedOxygenmgL ~ DO*Temperature+(1|Replicate), data = subset(Water.Qual, CyclePhase=="Hypoxia" & Year=="2"), REML=TRUE)
Anova(w2, test="F", type="III")

#DO during normoxic phase of cycle, Year 1
w3 <- lm(DissolvedOxygenmgL ~ DO*Temperature, data = subset(Water.Qual, CyclePhase=="Normoxia" & Year=="1"))
Anova(w3, test="F", type="III")

#DO during normoxic phase of cycle, Year 2
w4 <- lm(DissolvedOxygenmgL ~ DO*Temperature, data = subset(Water.Qual, CyclePhase=="Normoxia" & Year=="2"))
Anova(w4, test="F", type="III")

#Temperature, Year 1
w5 <- lmer(TemperatureC ~ DO*Temperature+(1|Replicate), data = subset(Water.Qual, CyclePhase=="Hypoxia" & Year=="1"), REML=TRUE)
Anova(w5, test="F", type="III")

#Temperature, Year 2
w6 <- lmer(TemperatureC ~ DO*Temperature+(1|Replicate), data = subset(Water.Qual, CyclePhase=="Hypoxia" & Year=="2"), REML=TRUE)
Anova(w6, test="F", type="III")

#pH, Year 1
w7 <- lmer(pH ~ DO*Temperature+(1|Replicate), data = subset(Water.Qual, Year=="1"), REML=TRUE)
Anova(w7, test="F", type="III")

#pH, Year 2
w8 <- lmer(pH ~ DO*Temperature+(1|Replicate), data = subset(Water.Qual, Year=="2"), REML=TRUE)
Anova(w8, test="F", type="III")

#Alkalinity, Year 1
w9 <- lm(TotalAlkalinityumolkg ~ Temperature, data = subset(Alk, Year=="1"))
Anova(w9, test="F", type="III")

#Alkalinity, Year 2
w10 <- lm(TotalAlkalinityumolkg ~ Temperature, data = subset(Alk, Year=="2"))
Anova(w10, test="F", type="III")

#pCO2, Year 1
w11 <- lm(pCO2uatm ~ Temperature, data = subset(Alk, Year=="1"))
Anova(w11, test="F", type="III")

#pCO2, Year 2
w12 <- lm(pCO2uatm ~ Temperature, data = subset(Alk, Year=="2"))
Anova(w12, test="F", type="III")

#### SKM 1,  Growth ####
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
emm1 <- emmeans(m1,specs = pairwise ~ Year1DO*Year1Temp, adjust = "none") 
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
emm2 <- emmeans(m2,specs = pairwise ~ Year2DO*Year2Temp, adjust = "none") 
emm2$emmeans
emm2$contrasts


##tissue:shell growth 
m3 <- lmer(Ratio_tissue_shell_mg~ Phase_1_DO*Phase_1_temp*Phase_2.1_DO*Phase_2.1_temp+Actual_shell_pre_mg+
             (1|Phase_2_rep)+(1|Phase_1_rep)+
             (1|Phase_2_rep:Phase_1_DO)+(1|Phase_2_rep:Phase_1_temp)+(1|Phase_2_rep:Phase_1_DO:Phase_1_temp), data = Growth_Data_forR, REML=TRUE)
Anova(m3, test="F", type="III")

#diagnostics
leveneTest(TissueShellGrowth~OverallTreatmentCombination, Growth_Data_forR) 
m3.e <- residuals(m3) 
qqnorm(m3.e) #not quite normal but ancova is robust to non-normality

#posthocs
emm3 <- emmeans(m3,specs = pairwise ~ Phase_1_DO:Phase_2.1_DO:Phase_2.1_temp, adjust="none")
emm3$emmeans
emm3$contrasts


##Nitrogen in tissue growth (mg)
m4 <- lmer(TissueGrowthNmg~ Year1DO*Year1Temp*Year2DO*Year2Temp+InitialShellMassmg+
             (1|Year2Rep)+(1|Year1Rep)+
             (1|Year2Rep:Year1DO)+(1|Year2Rep:Year1Temp)+(1|Year2Rep:Year1DO:Year1Temp), data = Growth_Data_forR, REML=TRUE)
Anova(m4, test="F", type="III")

#diagnostics
leveneTest(TissueGrowthNmg~OverallTreatmentCombination, Growth_Data_forR) 
m4.e <- residuals(m4) 
qqnorm(m4.e) #not quite normal but ancova is robust to non-normality

#posthocs
emm4 <- emmeans(m4,specs = pairwise ~ Year1DO*Year1Temp*Year2Temp, adjust="none")
emm4$emmeans
emm4$contrasts


##Nitrogen in shell growth (mg)
m5 <- lmer(ShellGrowthNmg~ Year1DO*Year1Temp*Year2DO*Year2Temp+InitialShellMassmg+
             (1|Year2Rep)+(1|Year1Rep)+
             (1|Year2Rep:Year1DO)+(1|Year2Rep:Year1Temp)+(1|Year2Rep:Year1DO:Year1Temp), data = Growth_Data_forR, REML=TRUE)
Anova(m5, test="F", type="III")

#diagnostics
leveneTest(ShellGrowthNmg~OverallTreatmentCombination, Growth_Data_forR) 
m5.e <- residuals(m5) 
qqnorm(m5.e) #not quite normal but ancova is robust to non-normality

#posthocs
emm5 <- emmeans(m5,specs = pairwise ~ Year1DO*Year1Temp*Year2Temp, adjust="none")
emm5$emmeans
emm5$contrasts

emm6 <- emmeans(m5,specs = pairwise ~ Year1DO*Year1Temp*Year2DO, adjust="none")
emm6$emmeans
emm6$contrasts


#generating allometric equations for Harris Creek reef extrapolation from experimental oysters
#different relationship for each of the 16 treatment combinations incorporates size differences
#I'm sure there's a cleaner way to do this

#Norm,Amb,Norm,Amb
WBWB.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Amb,Norm,Amb"))
summary(WBWB.allot)

WBWB.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Amb,Norm,Amb"))
summary(WBWB.allos)

#Norm,Amb,Norm,Warm
WBWP.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Amb,Norm,Warm"))
summary(WBWP.allot)

WBWP.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Amb,Norm,Warm"))
summary(WBWP.allos)

#Norm,Amb,Hyp,Amb
WBBB.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Amb,Hyp,Amb"))
summary(WBBB.allot)

WBBB.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Amb,Hyp,Amb"))
summary(WBBB.allos)

#Norm,Amb,Hyp,Warm
WBBP.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Amb,Hyp,Warm"))
summary(WBBP.allot)

WBBP.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Amb,Hyp,Warm"))
summary(WBWB.allos)

#Norm,Warm,Norm,Amb
WPWB.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Warm,Norm,Amb"))
summary(WPWB.allot)

WPWB.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Warm,Norm,Amb"))
summary(WPWB.allos)

#Norm,Warm,Norm,Warm
WPWP.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Warm,Norm,Warm"))
summary(WPWP.allot)

WPWP.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Warm,Norm,Warm"))
summary(WPWP.allos)

#Norm,Warm,Hyp,Amb
WPBB.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Warm,Hyp,Amb"))
summary(WPBB.allot)

WPBB.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Warm,Hyp,Amb"))
summary(WPBB.allos)

#Norm,Warm,Hyp,Warm
WPBP.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Warm,Hyp,Warm"))
summary(WPBP.allot)

WPBP.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Norm,Warm,Hyp,Warm"))
summary(WPBP.allos)

#Hyp,Amb,Norm,Amb
BBWB.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Amb,Norm,Amb"))
summary(BBWB.allot)

BBWB.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Amb,Norm,Amb"))
summary(BBWB.allos)

#Hyp,Amb,Norm,Warm
BBWP.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Amb,Norm,Warm"))
summary(BBWP.allot)

BBWP.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Amb,Norm,Warm"))
summary(BBWP.allos)

#Hyp,Amb,Hyp,Amb
BBBB.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Amb,Hyp,Amb"))
summary(BBBB.allot)

BBBB.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Amb,Hyp,Amb"))
summary(BBBB.allos)

#Hyp,Amb,Hyp,Warm
BBBP.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Amb,Hyp,Warm"))
summary(BBBP.allot)

BBBP.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Amb,Hyp,Warm"))
summary(BBBP.allos)

#Hyp,Warm,Norm,Amb
BPWB.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Warm,Norm,Amb"))
summary(BPWB.allot)

BPWB.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Warm,Norm,Amb"))
summary(BPWB.allos)

#Hyp,Warm,Norm,Warm
BPWP.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Warm,Norm,Warm"))
summary(BPWP.allot)

BPWP.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Warm,Norm,Warm"))
summary(BPWP.allos)

#Hyp,Warm,Hyp,Amb
BPBB.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Warm,Hyp,Amb"))
summary(BPBB.allot)

BPBB.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Warm,Hyp,Amb"))
summary(BPBB.allos)

#Hyp,Warm,Hyp,Warm
BPBP.allot <- lm(lnTissuemassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Warm,Hyp,Warm"))
summary(BPBP.allot)

BPBP.allos <- lm(lnShellmassmg ~ lnShelllengthmm, data = subset(allo.equations, OverallTreatmentCombination=="Hyp,Warm,Hyp,Warm"))
summary(BPBP.allos)


##bootstrapping means of extrapolated N on restored oyster reef
#using allometric equations from above and shell lengths of ALL oysters from 2020 Maryland Oyster Monitoring Report
# to generate bootstrapped mean tissueand shell N and 95% CIs

#in a potential Normoxic,Ambient environment

#Normoxic,Ambient early life
#tissue for mean
WBWBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.456)*(sl^1.98))/1000)*0.07517*121.42*0.001*4046.86)
}
#equation notation of e : first set of parentheses = allometric equation from Table S5
#0.07517 = proportion of tissue that is N from Table S4
#121.42 density of oysters in m2 from 2020 report
#0.001 convert g N to kg N
#4046.86 = convert m2 to acre

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, WBWBt, R=1000)
myBootstrap.tiss

#shell for mean
WBWBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.158)*(sl^1.63))/1000)*0.001083*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, WBWBs, R=1000)
myBootstrap.shell

#ttl for CI
WBWBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.456)*(sl^1.98))/1000)*0.07517*121.42*0.001*4046.86)
  s <- (((exp(2.158)*(sl^1.63))/1000)*0.001083*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, WBWBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Normoxic,Warm early life
#tissue for mean
WPWBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.985)*(sl^1.76))/1000)*0.08817*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, WPWBt, R=1000)
myBootstrap.tiss

#shell for mean
WPWBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.795)*(sl^1.40))/1000)*0.0012*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, WPWBs, R=1000)
myBootstrap.shell

#ttl for CI
WPWBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.985)*(sl^1.76))/1000)*0.08817*121.42*0.001*4046.86)
  s <- (((exp(2.795)*(sl^1.40))/1000)*0.0012*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, WPWBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Ambient early life
#tissue for mean
BBWBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.487)*(sl^1.94))/1000)*0.09567*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, BBWBt, R=1000)
myBootstrap.tiss

#shell for mean
BBWBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.347)*(sl^1.55))/1000)*0.001333*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, BBWBs, R=1000)
myBootstrap.shell

#ttl for CI
BBWBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <-  (((exp(0.487)*(sl^1.94))/1000)*0.09567*121.42*0.001*4046.86)
  s <- (((exp(2.347)*(sl^1.55))/1000)*0.001333*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, BBWBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Warm early life
#tissue for mean
BPWBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.339)*(sl^1.99))/1000)*0.08533*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, BPWBt, R=1000)
myBootstrap.tiss

#shell for mean
BPWBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.22)*(sl^1.59))/1000)*0.001*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, BPWBs, R=1000)
myBootstrap.shell

#ttl for CI
BPWBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.339)*(sl^1.99))/1000)*0.08533*121.42*0.001*4046.86)
  s <- (((exp(2.22)*(sl^1.59))/1000)*0.001*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, BPWBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#in a potential Normoxic,Warm environment

#Normoxic,Ambient early life
#tissue for mean
WBWPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.73)*(sl^1.50))/1000)*0.09683*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, WBWPt, R=1000)
myBootstrap.tiss

#shell for mean
WBWPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(3.343)*(sl^1.20))/1000)*0.0018*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, WBWPs, R=1000)
myBootstrap.shell

#ttl for CI
WBWPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.73)*(sl^1.50))/1000)*0.09683*121.42*0.001*4046.86)
  s <- (((exp(3.343)*(sl^1.20))/1000)*0.0018*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, WBWPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Normoxic,Warm early life
#tissue for mean
WPWPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.851)*(sl^1.44))/1000)*0.08617*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, WPWPt, R=1000)
myBootstrap.tiss

#shell for mean
WPWPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(3.825)*(sl^1.01))/1000)*0.0017*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, WPWPs, R=1000)
myBootstrap.shell

#ttl for CI
WPWPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.851)*(sl^1.44))/1000)*0.08617*121.42*0.001*4046.86)
  s <- (((exp(3.825)*(sl^1.01))/1000)*0.0017*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, WPWPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Ambient early life
#tissue for mean
BBWPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.487)*(sl^1.55))/1000)*0.086*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, BBWPt, R=1000)
myBootstrap.tiss

#shell for mean
BBWPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(3.192)*(sl^1.23))/1000)*0.001617*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, BBWPs, R=1000)
myBootstrap.shell

#ttl for CI
BBWPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.487)*(sl^1.55))/1000)*0.086*121.42*0.001*4046.86)
  s <- (((exp(3.192)*(sl^1.23))/1000)*0.001617*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, BBWPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Warm early life
#tissue for mean
BPWPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.631)*(sl^1.52))/1000)*0.08766*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, BPWPt, R=1000)
myBootstrap.tiss

#shell for mean
BPWPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(3.005)*(sl^1.30))/1000)*0.0013*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, BPWPs, R=1000)
myBootstrap.shell

#ttl for CI
BPWPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.631)*(sl^1.52))/1000)*0.08766*121.42*0.001*4046.86)
  s <- (((exp(3.005)*(sl^1.30))/1000)*0.0013*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, BPWPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#in a potential Hypoxic,Ambient environment

#Normoxic,Ambient early life
#tissue for mean
WBBBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.34)*(sl^1.62))/1000)*0.095*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, WBBBt, R=1000)
myBootstrap.tiss

#shell for mean
WBBBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.962)*(sl^1.32))/1000)*0.001483*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, WBBBs, R=1000)
myBootstrap.shell

#ttl for CI
WBBBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.34)*(sl^1.62))/1000)*0.095*121.42*0.001*4046.86)
  s <- (((exp(2.962)*(sl^1.32))/1000)*0.001483*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, WBBBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Normoxic,Warm early life
#tissue for mean
WPBBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.763)*(sl^1.82))/1000)*0.08766*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, WPBBt, R=1000)
myBootstrap.tiss

#shell for mean
WPBBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.533)*(sl^1.47))/1000)*0.0015*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, WPBBs, R=1000)
myBootstrap.shell

#ttl for CI
WPBBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <-(((exp(0.763)*(sl^1.82))/1000)*0.08766*121.42*0.001*4046.86)
  s <- (((exp(2.533)*(sl^1.47))/1000)*0.0015*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, WPBBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Ambient early life
#tissue for mean
BBBBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.915)*(sl^1.74))/1000)*0.101*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, BBBBt, R=1000)
myBootstrap.tiss

#shell for mean
BBBBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.573)*(sl^1.44))/1000)*0.001366*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, BBBBs, R=1000)
myBootstrap.shell

#ttl for CI
BBBBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.915)*(sl^1.74))/1000)*0.101*121.42*0.001*4046.86)
  s <- (((exp(2.573)*(sl^1.44))/1000)*0.001366*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, BBBBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Warm early life
#tissue for mean
BPBBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.436)*(sl^1.55))/1000)*0.08666*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, BPBBt, R=1000)
myBootstrap.tiss

#shell for mean
BPBBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.895)*(sl^1.33))/1000)*0.001283*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, BPBBs, R=1000)
myBootstrap.shell

#ttl for CI
BPBBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.436)*(sl^1.55))/1000)*0.08666*121.42*0.001*4046.86)
  s <- (((exp(2.895)*(sl^1.33))/1000)*0.001283*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, BPBBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#in a Hypoxic,Warm potential reef environment

#Normoxic,Ambient early life
#tissue for mean
WBBPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.911)*(sl^1.83))/1000)*0.0875*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, WBBPt, R=1000)
myBootstrap.tiss

#shell for mean
WBBPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.642)*(sl^1.47))/1000)*0.0015*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, WBBPs, R=1000)
myBootstrap.shell

#ttl for CI
WBBPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.911)*(sl^1.83))/1000)*0.0875*121.42*0.001*4046.86)
  s <- (((exp(2.642)*(sl^1.47))/1000)*0.0015*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, WBBPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Normoxic,Warm early life
#tissue for mean
WPBPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.917)*(sl^1.80))/1000)*0.08567*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, WPBPt, R=1000)
myBootstrap.tiss

#shell for mean
WPBPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.586)*(sl^1.47))/1000)*0.00138*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, WPBPs, R=1000)
myBootstrap.shell

#ttl for CI
WPBPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.917)*(sl^1.80))/1000)*0.08567*121.42*0.001*4046.86)
  s <- (((exp(2.586)*(sl^1.47))/1000)*0.00138*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, WPBPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Ambient early life
#tissue for mean
BBBPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.967)*(sl^1.78))/1000)*0.084*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, BBBPt, R=1000)
myBootstrap.tiss

#shell for mean
BBBPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.755)*(sl^1.41))/1000)*0.001233*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, BBBPs, R=1000)
myBootstrap.shell

#ttl for CI
BBBPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.967)*(sl^1.78))/1000)*0.084*121.42*0.001*4046.86)
  s <- (((exp(2.755)*(sl^1.41))/1000)*0.001233*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, BBBPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Warm early life
#tissue for mean
BPBPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.765)*(sl^1.86))/1000)*0.0893*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.allo, BPBPt, R=1000)
myBootstrap.tiss

#shell for mean
BPBPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.435)*(sl^1.52))/1000)*0.001317*121.42*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.allo, BPBPs, R=1000)
myBootstrap.shell
plot(myBootstrap.shell)

#ttl for CI
BPBPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.765)*(sl^1.86))/1000)*0.0893*121.42*0.001*4046.86)
  s <- (((exp(2.435)*(sl^1.52))/1000)*0.001317*121.42*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.allo, BPBPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))


####bootstrapping means of extrapolated N on restored oyster reef
#using allometric equations from above and shell lengths OF SPAT-SIZED OYSTERS ONLY (<40mm) from 2020 Maryland Oyster Monitoring Report
#to generate bootstrapped mean tissue and shell N and 95% CIs

#the only difference between the code for all oysters above and spat-sized only oysters here is that this is run on the spat sized subset of data only and
#the density is of spat-sized oysters only (spat-sized oysters: 38.76 oysters m-2 from the 2020 report).

#subset data so using spat-sixed oysters only
bootstrap.alloSUB <- bootstrap.allo[ which(bootstrap.allo$ShellLengthmm < 40), ]
bootstrap.alloSUB <- as.data.frame(bootstrap.alloSUB)

#in a potential Normoxic,Ambient environment

#Normoxic,Ambient early life
#tissue for mean

WBWBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.456)*(sl^1.98))/1000)*0.07517*38.76*0.001*4046.86)
}
#equation notation of e : first set of parentheses = allometric equation from Table S5
#0.07517 = proportion of tissue that is N from Table S4
#spat density: 38.76 in m2 from 2020 report 
#0.001 convert g N to kg N
#4046.86 = convert m2 to acre

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, WBWBt, R=1000)
myBootstrap.tiss

#shell for mean
WBWBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.158)*(sl^1.63))/1000)*0.001083*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, WBWBs, R=1000)
myBootstrap.shell

#ttl for CI
WBWBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.456)*(sl^1.98))/1000)*0.07517*38.76*0.001*4046.86)
  s <- (((exp(2.158)*(sl^1.63))/1000)*0.001083*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, WBWBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Normoxic,Warm early life
#tissue for mean
WPWBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.985)*(sl^1.76))/1000)*0.08817*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, WPWBt, R=1000)
myBootstrap.tiss

#shell for mean
WPWBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.795)*(sl^1.40))/1000)*0.0012*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, WPWBs, R=1000)
myBootstrap.shell

#ttl for CI
WPWBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.985)*(sl^1.76))/1000)*0.08817*38.76*0.001*4046.86)
  s <- (((exp(2.795)*(sl^1.40))/1000)*0.0012*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, WPWBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Ambient early life
#tissue for mean
BBWBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.487)*(sl^1.94))/1000)*0.09567*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, BBWBt, R=1000)
myBootstrap.tiss

#shell for mean
BBWBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.347)*(sl^1.55))/1000)*0.001333*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, BBWBs, R=1000)
myBootstrap.shell

#ttl for CI
BBWBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <-  (((exp(0.487)*(sl^1.94))/1000)*0.09567*38.76*0.001*4046.86)
  s <- (((exp(2.347)*(sl^1.55))/1000)*0.001333*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, BBWBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Warm early life
#tissue for mean
BPWBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.339)*(sl^1.99))/1000)*0.08533*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, BPWBt, R=1000)
myBootstrap.tiss

#shell for mean
BPWBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.22)*(sl^1.59))/1000)*0.001*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, BPWBs, R=1000)
myBootstrap.shell

#ttl for CI
BPWBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.339)*(sl^1.99))/1000)*0.08533*38.76*0.001*4046.86)
  s <- (((exp(2.22)*(sl^1.59))/1000)*0.001*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, BPWBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#in a potential Normoxic,Warm environment

#Normoxic,Ambient early life
#tissue for mean
WBWPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.73)*(sl^1.50))/1000)*0.09683*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, WBWPt, R=1000)
myBootstrap.tiss

#shell for mean
WBWPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(3.343)*(sl^1.20))/1000)*0.0018*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, WBWPs, R=1000)
myBootstrap.shell

#ttl for CI
WBWPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.73)*(sl^1.50))/1000)*0.09683*38.76*0.001*4046.86)
  s <- (((exp(3.343)*(sl^1.20))/1000)*0.0018*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, WBWPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Normoxic,Warm early life
#tissue for mean
WPWPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.851)*(sl^1.44))/1000)*0.08617*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, WPWPt, R=1000)
myBootstrap.tiss

#shell for mean
WPWPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(3.825)*(sl^1.01))/1000)*0.0017*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, WPWPs, R=1000)
myBootstrap.shell

#ttl for CI
WPWPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.851)*(sl^1.44))/1000)*0.08617*38.76*0.001*4046.86)
  s <- (((exp(3.825)*(sl^1.01))/1000)*0.0017*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, WPWPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Ambient early life
#tissue for mean
BBWPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.487)*(sl^1.55))/1000)*0.086*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, BBWPt, R=1000)
myBootstrap.tiss

#shell for mean
BBWPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(3.192)*(sl^1.23))/1000)*0.001617*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, BBWPs, R=1000)
myBootstrap.shell

#ttl for CI
BBWPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.487)*(sl^1.55))/1000)*0.086*38.76*0.001*4046.86)
  s <- (((exp(3.192)*(sl^1.23))/1000)*0.001617*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, BBWPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Warm early life
#tissue for mean
BPWPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.631)*(sl^1.52))/1000)*0.08766*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, BPWPt, R=1000)
myBootstrap.tiss

#shell for mean
BPWPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(3.005)*(sl^1.30))/1000)*0.0013*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, BPWPs, R=1000)
myBootstrap.shell

#ttl for CI
BPWPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.631)*(sl^1.52))/1000)*0.08766*38.76*0.001*4046.86)
  s <- (((exp(3.005)*(sl^1.30))/1000)*0.0013*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, BPWPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#in a potential Hypoxic,Ambient environment

#Normoxic,Ambient early life
#tissue for mean
WBBBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.34)*(sl^1.62))/1000)*0.095*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, WBBBt, R=1000)
myBootstrap.tiss

#shell for mean
WBBBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.962)*(sl^1.32))/1000)*0.001483*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, WBBBs, R=1000)
myBootstrap.shell

#ttl for CI
WBBBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.34)*(sl^1.62))/1000)*0.095*38.76*0.001*4046.86)
  s <- (((exp(2.962)*(sl^1.32))/1000)*0.001483*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, WBBBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Normoxic,Warm early life
#tissue for mean
WPBBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.763)*(sl^1.82))/1000)*0.08766*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, WPBBt, R=1000)
myBootstrap.tiss

#shell for mean
WPBBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.533)*(sl^1.47))/1000)*0.0015*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, WPBBs, R=1000)
myBootstrap.shell

#ttl for CI
WPBBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <-(((exp(0.763)*(sl^1.82))/1000)*0.08766*38.76*0.001*4046.86)
  s <- (((exp(2.533)*(sl^1.47))/1000)*0.0015*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, WPBBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Ambient early life
#tissue for mean
BBBBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.915)*(sl^1.74))/1000)*0.101*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, BBBBt, R=1000)
myBootstrap.tiss

#shell for mean
BBBBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.573)*(sl^1.44))/1000)*0.001366*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, BBBBs, R=1000)
myBootstrap.shell

#ttl for CI
BBBBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.915)*(sl^1.74))/1000)*0.101*38.76*0.001*4046.86)
  s <- (((exp(2.573)*(sl^1.44))/1000)*0.001366*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, BBBBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Warm early life
#tissue for mean
BPBBt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(1.436)*(sl^1.55))/1000)*0.08666*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, BPBBt, R=1000)
myBootstrap.tiss

#shell for mean
BPBBs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.895)*(sl^1.33))/1000)*0.001283*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, BPBBs, R=1000)
myBootstrap.shell

#ttl for CI
BPBBttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(1.436)*(sl^1.55))/1000)*0.08666*38.76*0.001*4046.86)
  s <- (((exp(2.895)*(sl^1.33))/1000)*0.001283*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, BPBBttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#in a Hypoxic,Warm potential reef environment

#Normoxic,Ambient early life
#tissue for mean
WBBPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.911)*(sl^1.83))/1000)*0.0875*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, WBBPt, R=1000)
myBootstrap.tiss

#shell for mean
WBBPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.642)*(sl^1.47))/1000)*0.0015*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, WBBPs, R=1000)
myBootstrap.shell

#ttl for CI
WBBPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.911)*(sl^1.83))/1000)*0.0875*38.76*0.001*4046.86)
  s <- (((exp(2.642)*(sl^1.47))/1000)*0.0015*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, WBBPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Normoxic,Warm early life
#tissue for mean
WPBPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.917)*(sl^1.80))/1000)*0.08567*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, WPBPt, R=1000)
myBootstrap.tiss

#shell for mean
WPBPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.586)*(sl^1.47))/1000)*0.00138*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, WPBPs, R=1000)
myBootstrap.shell

#ttl for CI
WPBPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.917)*(sl^1.80))/1000)*0.08567*38.76*0.001*4046.86)
  s <- (((exp(2.586)*(sl^1.47))/1000)*0.00138*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, WPBPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Ambient early life
#tissue for mean
BBBPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.967)*(sl^1.78))/1000)*0.084*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, BBBPt, R=1000)
myBootstrap.tiss

#shell for mean
BBBPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.755)*(sl^1.41))/1000)*0.001233*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, BBBPs, R=1000)
myBootstrap.shell

#ttl for CI
BBBPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.967)*(sl^1.78))/1000)*0.084*38.76*0.001*4046.86)
  s <- (((exp(2.755)*(sl^1.41))/1000)*0.001233*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, BBBPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

#Hypoxic,Warm early life
#tissue for mean
BPBPt <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(0.765)*(sl^1.86))/1000)*0.0893*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.tiss <- boot(bootstrap.alloSUB, BPBPt, R=1000)
myBootstrap.tiss

#shell for mean
BPBPs <- function(data, indices, sl){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  e <- (((exp(2.435)*(sl^1.52))/1000)*0.001317*38.76*0.001*4046.86)
}

set.seed(12345)
myBootstrap.shell <- boot(bootstrap.alloSUB, BPBPs, R=1000)
myBootstrap.shell

#ttl for CI
BPBPttl <- function(data, indices, t,s){
  dt<-data[indices,]
  sl <- mean(dt[indices])
  t <- (((exp(0.765)*(sl^1.86))/1000)*0.0893*38.76*0.001*4046.86)
  s <- (((exp(2.435)*(sl^1.52))/1000)*0.001317*38.76*0.001*4046.86)
  ttl <- t+s
}

set.seed(12345)
myBootstrap.ttl <- boot(bootstrap.alloSUB, BPBPttl, R=1000)
myBootstrap.ttl
boot.ci(myBootstrap.ttl, type=c('basic','norm'))

