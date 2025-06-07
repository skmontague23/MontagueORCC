## script just to format Orion data for analysis 
## Last modified: JGM 06/24/2024

# loading in packages
library(tidyverse)
library(stringr)

##### PHASE 1 ##### code adapted from JGM
# replace with location of repo directory/file location
data <- read_csv("~/Desktop/MontagueORCC/ORION_wq_data/phase1_master_orion.csv")
head(data)
data <- janitor::clean_names(data)
colnames(data)
View(data)
# can choose here what you want - I've chosen the date, time, temperature, and DO (&pH)
data1 <- select(data, tank_num, date_time_2, temperature_5, concentration_value, p_h_value)

# adding another column to data1 for Treatment, warm treat, and DO treat
data1$treatment <- 'empty'
data1$temp_treat <- 'empty'
data1$DO_treat <- 'empty'

# now assigning treatments to tank numbers
# Treatment = common name which refers to both DO and temp treatment
# temp_treat = temperature treatment (warm or ambient)
# DO_treat = DO treatment (hypoxic or normoxic)

# for loop and ifelse statements to assign treatments
for (i in 1:length(data$tank_num)) { # goes through each row and looks at the tank number
  tank <- data$tank_num[i] # indexing to select a single tank 
  if (tank <= 6) { # if tank is less than or equal to 6, it's a control tank
    data1$treatment[i] <- 'control'
    data1$temp_treat[i] <- 'ambient'
    data1$DO_treat[i] <- 'normoxic'
  } else if (7 <= tank && tank <= 12) { # if tank is between 7-12, it's a hypoxic tank
    data1$treatment[i] <- 'hypoxic'
    data1$temp_treat[i] <- 'ambient'
    data1$DO_treat[i] <- 'hypoxic'
  } else if (13 <= tank && tank <= 18) { # if tank is between 13-18, it's a both tank
    data1$treatment[i] <- 'both'
    data1$temp_treat[i] <- 'warm'
    data1$DO_treat[i] <- 'hypoxic'
  } else if (tank >= 19) { # if tank is greater than or equal to 19, it's a warm tank
    data1$treatment[i] <- 'warm' 
    data1$temp_treat[i] <- 'warm' 
    data1$DO_treat[i] <- 'normoxic'
  }
}

# have R recognize data and time format
data1$date_time <- as.POSIXct(data1$date_time_2, format = "%m/%d/%Y %H:%M")

# but can also split date and time into their own columns as well
split_text <- strsplit(data1$date_time_2, " ") # split text by a space
split_df <- data.frame(do.call(rbind, split_text)) # create a new df of our new 2 columns
colnames(split_df) <- c('date','time') # rename the columns to date and time
data2 <- cbind(split_df, data1) # add the new columns to our OG data frame 
colnames(data2)
data2 <- select(data2, date, time, tank_num, date_time_2, temperature_5, 
                concentration_value, treatment, DO_treat, temp_treat, date_time, p_h_value)
# want to remove the ':' in time so I can manipulate later (so 1230 instead of 12:30)
data2$time <- gsub(":", "", data2$time)
# get R to recongize the date column as dates
data2$date <- as.Date(data2$date, format = "%m/%d/%Y")

# assigning the cycle at which the reading was taken:
# morning = normoxia, afternoon = hypoxia

# create new column for cycle 
data2$cycle <- 'empty'
# time is now a numeric value (and remove rows with NAs)
data2$time <- as.numeric(data2$time)
data2 <- na.omit(data2)
# for loop and ifelse statements to assign treatments
for (i in 1:length(data2$time)) {
  t <- data2$time[i]
  if (t < 1200) { # if it's the morning (before noon), HOPE is in normoxia
    data2$cycle[i] <- 'normoxic'
  } else if (t > 1200) { # if it's the afternoon (after noon), HOPE is in hypoxia
    data2$cycle[i] <- 'hypoxic'
  } 
}

# removing extra characters from temperature and DO concentration

# removing *C from temp
data2$temperature_5 <- as.numeric(str_extract(data2$temperature_5, "\\d+\\.?\\d*"))
# removing mg/L from DO concentration
data2$concentration_value <- gsub("mg/L", "", data2$concentration_value)
data2$concentration_value <- as.numeric(data2$concentration_value)

# removing "pH" from pH value
data2$p_h_value <- gsub("pH", "", data2$p_h_value)
data2$p_h_value <- as.numeric(data2$p_h_value)

#rename cols to include units
View(data2)
colnames(data2)

phase1_formatted_orion <- data2%>%
  select(-date_time_2) %>% #remove redundant column
  rename(temp_C = temperature_5)%>%
  rename(oxygen_conc_mgL = concentration_value) %>%
  rename(pH = p_h_value) #rename a few columns

View(phase1_formatted_orion)

# now formatted df is called phase1_formatted_orion
# ready for analysis

write.csv(phase1_formatted_orion, "phase1_formatted_orion.csv", row.names = FALSE)
getwd()


##### PHASE 2 ####### code adapted from JGM

# replace with location of repo directory/file location
data <- read_csv("~/Desktop/MontagueORCC/ORION_wq_data/phase2_master_orion.csv")
data <- janitor::clean_names(data)
colnames(data)
head(data)
# can choose here what you want - I've chosen the date, time, temperature, and DO 
data1 <- select(data, tank_num, date_time_2, temperature_5, concentration_value, p_h_value)

# adding another column to data1 for Treatment, warm treat, and DO treat
data1$treatment <- 'empty'
data1$temp_treat <- 'empty'
data1$DO_treat <- 'empty'

# now assigning treatments to tank numbers
# Treatment = common name which refers to both DO and temp treatment
# temp_treat = temperature treatment (warm or ambient)
# DO_treat = DO treatment (hypoxic or normoxic)

# for loop and ifelse statements to assign treatments
for (i in 1:length(data$tank_num)) { # goes through each row and looks at the tank number
  tank <- data$tank_num[i] # indexing to select a single tank 
  if (tank <= 6) { # if tank is less than or equal to 6, it's a control tank
    data1$treatment[i] <- 'control'
    data1$temp_treat[i] <- 'ambient'
    data1$DO_treat[i] <- 'normoxic'
  } else if (7 <= tank && tank <= 12) { # if tank is between 7-12, it's a hypoxic tank
    data1$treatment[i] <- 'hypoxic'
    data1$temp_treat[i] <- 'ambient'
    data1$DO_treat[i] <- 'hypoxic'
  } else if (13 <= tank && tank <= 18) { # if tank is between 13-18, it's a both tank
    data1$treatment[i] <- 'both'
    data1$temp_treat[i] <- 'warm'
    data1$DO_treat[i] <- 'hypoxic'
  } else if (tank >= 19) { # if tank is greater than or equal to 19, it's a warm tank
    data1$treatment[i] <- 'warm' 
    data1$temp_treat[i] <- 'warm' 
    data1$DO_treat[i] <- 'normoxic'
  }
}

# have R recognize data and time format
data1$date_time <- as.POSIXct(data1$date_time_2, format = "%m/%d/%Y %H:%M")

# but can also split date and time into their own columns as well
split_text <- strsplit(data1$date_time_2, " ") # split text by a space
split_df <- data.frame(do.call(rbind, split_text)) # create a new df of our new 2 columns
colnames(split_df) <- c('date','time') # rename the columns to date and time
data2 <- cbind(split_df, data1) # add the new columns to our OG data frame 
colnames(data2)
data2 <- select(data2, date, time, tank_num, date_time_2, temperature_5, concentration_value, treatment, DO_treat, temp_treat, date_time, p_h_value)
# want to remove the ':' in time so I can manipulate later (so 1230 instead of 12:30)
data2$time <- gsub(":", "", data2$time)
# get R to recongize the date column as dates
data2$date <- as.Date(data2$date, format = "%m/%d/%Y")

# assigning the cycle at which the reading was taken:
# morning = normoxia, afternoon = hypoxia

# create new column for cycle 
data2$cycle <- 'empty'
# time is now a numeric value (and remove rows with NAs)
data2$time <- as.numeric(data2$time)
data2 <- na.omit(data2)
# for loop and ifelse statements to assign treatments
for (i in 1:length(data2$time)) {
  t <- data2$time[i]
  if (t < 1200) { # if it's the morning (before noon), HOPE is in normoxia
    data2$cycle[i] <- 'normoxic'
  } else if (t > 1200) { # if it's the afternoon (after noon), HOPE is in hypoxia
    data2$cycle[i] <- 'hypoxic'
  } 
}

# removing extra characters from temperature and DO concentration

# removing *C from temp
data2$temperature_5 <- as.numeric(str_extract(data2$temperature_5, "\\d+\\.?\\d*"))
# removing mg/L from DO concentration
data2$concentration_value <- gsub("mg/L", "", data2$concentration_value)
data2$concentration_value <- as.numeric(data2$concentration_value)

# removing "pH" from pH value
data2$p_h_value <- gsub("pH", "", data2$p_h_value)
data2$p_h_value <- as.numeric(data2$p_h_value)

#rename cols to include units
View(data2)
colnames(data2)

phase2_formatted_orion <- data2%>%
  select(-date_time_2) %>% #remove redundant column
  rename(temp_C = temperature_5)%>%
  rename(oxygen_conc_mgL = concentration_value)%>%
  rename(pH = p_h_value)

View(phase2_formatted_orion)

# now formatted df is called phase2_formatted_orion
# ready for analysis

write.csv(phase2_formatted_orion, "phase2_formatted_orion.csv", row.names = FALSE)
getwd()
setwd("~/Desktop/MontagueORCC_repo/MontagueORCC/ORION_wq_data")

###SKM WQ analysis

View(phase2_formatted_orion)
colnames(phase2_formatted_orion)
library(lme4)
library(readr)
library(car)
library(emmeans)
library(boot)
library(esquisse)
library(ggplot2)
library(dplyr)

#read in data
phase1_formatted_orion <- read_csv("phase1_formatted_orion.csv")
phase2_formatted_orion <- read_csv("phase2_formatted_orion.csv")
phase2_formatted_orion <- phase2_formatted_orion%>%
  filter(measurement_error != "Y" | is.na(measurement_error))
min(phase2_formatted_orion$pH)
head(phase1_formatted_orion)
head(phase2_formatted_orion)

#need to adjust the time frame for the phase 1 temperature exposure
#oysters only experienced warming after 06/11/24
colnames(phase1_formatted_orion)
colnames(phase2_formatted_orion)
head(phase1_formatted_orion)

phase1_temp_orion <- phase1_formatted_orion%>%
  filter(date >= "2024-06-11")

#check to see what treatments were applied before the 6/11
before_phase1_temp <- phase1_formatted_orion%>%
  filter(date < "2024-06-11")

View(before_phase1_temp)
#hypoxia, DO
model1 <- lmer(oxygen_conc_mgL ~ DO_treat*temp_treat+(1|tank_num), data = subset(before_phase1_temp, cycle=="hypoxic"), REML=TRUE)
Anova(model1, test="F", type="III")
emmeans(model1, pairwise ~ DO_treat, adjust = "none")
emmeans(model1, pairwise ~ temp_treat, adjust = "none")
emmeans(model1, pairwise ~ DO_treat*temp_treat, adjust = "none")
#normoxia, DO
model2 <- lmer(oxygen_conc_mgL ~ DO_treat*temp_treat+(1|tank_num), data = subset(before_phase1_temp, cycle=="normoxic"), REML=TRUE)
Anova(model2, test="F", type="III")
emmeans(model2, pairwise ~ DO_treat, adjust = "none")
#temperature
model3 <- lmer(temp_C ~ DO_treat*temp_treat+(1|tank_num), data = before_phase1_temp, REML=TRUE)
Anova(model3, test="F", type="III")
emmeans(model3, pairwise ~ temp_treat, adjust = "none")
colnames(before_phase1_temp)
#graphs
phase1_formatted_orion$date_time <- as.POSIXct(gsub(" UTC", "", phase1_formatted_orion$date_time), format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
phase2_formatted_orion$date_time <- as.POSIXct(phase2_formatted_orion$date_time, 
                                               format = "%m/%d/%y %H:%M", 
                                               tz = "UTC")

#visualize phase 1 pH
ggplot(phase1_formatted_orion) +
  aes(x = date_time, y = pH, colour = treatment) +
  geom_point() +
  theme_classic()
#phase 2
ggplot(phase2_formatted_orion) +
  aes(x = date_time, y = pH, colour = treatment) +
  geom_point() +
  theme_classic()
sum(is.na(phase2_formatted_orion$date_time))
View(phase2_formatted_orion)

#visualize phase 1 temp
ggplot(phase1_formatted_orion) +
  aes(x = date_time, y = temp_C, colour = treatment) +
  geom_point() +
  theme_classic()

#visualize phase 2 temp
ggplot(phase2_formatted_orion) +
  aes(x = date_time, y = temp_C, colour = treatment) +
  geom_point() +
  theme_classic()

#visualize phase 1 DO
ggplot(subset(phase1_formatted_orion, cycle=="normoxic")) +
  aes(x = date_time, y = oxygen_conc_mgL, colour = treatment) +
  geom_point() +
  theme_classic()

ggplot(subset(phase1_formatted_orion, cycle=="hypoxic")) +
  aes(x = date_time, y = oxygen_conc_mgL, colour = treatment) +
  geom_point() +
  theme_classic()

#visualize phase 2 DO

ggplot(subset(phase2_formatted_orion, cycle=="normoxic")) +
  aes(x = treatment, y = oxygen_conc_mgL, colour = treatment) +
  geom_point() +
  theme_classic()

ggplot(subset(phase2_formatted_orion, cycle=="normoxic")) +
  aes(x = date_time, y = oxygen_conc_mgL, colour = treatment) +
  geom_point() +
  theme_classic()

ggplot(subset(phase2_formatted_orion, cycle=="hypoxic")) +
  aes(x = date_time, y = oxygen_conc_mgL, colour = treatment) +
  geom_point() +
  theme_classic()


#stats
#set contrasts
options(contrasts = c("contr.sum","contr.poly"))
getOption("contrasts") 

#DO during hypoxic phase of cycle, Phase 1
w1 <- lmer(oxygen_conc_mgL ~ DO_treat*temp_treat+(1|tank_num), data = subset(phase1_formatted_orion, cycle=="hypoxic"), REML=TRUE)
Anova(w1, test="F", type="III")
emmeans(w1, pairwise ~ DO_treat, adjust = "none")
emmeans(w1, pairwise ~ temp_treat, adjust = "none")

#DO during hypoxic phase of cycle, Phase 2
w2 <- lmer(oxygen_conc_mgL ~ DO_treat*temp_treat+(1|tank_num), data = subset(phase2_formatted_orion, cycle=="hypoxic"), REML=TRUE)
Anova(w2, test="F", type="III")
emmeans(w2, pairwise ~ DO_treat, adjust = "none")
emmeans(w2, pairwise ~ temp_treat, adjust = "none")
emmeans(w2, pairwise ~ DO_treat*temp_treat, adjust = "none")

#DO during normoxic phase of cycle, Phase 1
w3 <- lmer(oxygen_conc_mgL ~ DO_treat*temp_treat+(1|tank_num), data = subset(phase1_formatted_orion, cycle=="normoxic"), REML=TRUE)
Anova(w3, test="F", type="III")

#DO during normoxic phase of cycle, Phase 2
w4 <- lmer(oxygen_conc_mgL ~ DO_treat*temp_treat+(1|tank_num), data = subset(phase2_formatted_orion, cycle=="normoxic"), REML=TRUE)
Anova(w4, test="F", type="III")
emmeans(w4, pairwise ~ temp_treat, adjust = "none")

#Temperature, Phase 1
w5 <- lmer(temp_C ~ DO_treat*temp_treat+(1|tank_num), data = phase1_temp_orion, REML=TRUE)
Anova(w5, test="F", type="III")
emmeans(w5, pairwise ~ temp_treat, adjust = "none")

#the whole time frame
w5 <- lmer(temp_C ~ DO_treat*temp_treat+(1|tank_num), data = phase1_formatted_orion, REML=TRUE)
Anova(w5, test="F", type="III")
emmeans(w5, pairwise ~ temp_treat, adjust = "none")


#Temperature, Phase 2 
w6 <- lmer(temp_C ~ DO_treat*temp_treat+(1|tank_num), data = phase2_formatted_orion, REML=TRUE)
Anova(w6, test="F", type="III")
emmeans(w6, pairwise ~ DO_treat, adjust = "none")
emmeans(w6, pairwise ~ temp_treat, adjust = "none")
emmeans(w6, pairwise ~ DO_treat*temp_treat, adjust = "none")

#pH, Phase 1
w7 <- lmer(pH ~ DO_treat*temp_treat+(1|tank_num), data = phase1_formatted_orion, REML=TRUE)
Anova(w7, test="F", type="III")
emmeans(w7, pairwise ~ DO_treat, adjust = "none")
emmeans(w7, pairwise ~ temp_treat, adjust = "none")

View(phase1_formatted_orion)
controlpH <- phase1_formatted_orion%>%
  filter(treatment == "control")
View(controlpH)
mean(controlpH$pH)

View(phase1_formatted_orion)
hyppH <- phase1_formatted_orion%>%
  filter(treatment == "hypoxic")
View(hyppH)
mean(hyppH$pH)

View(phase1_formatted_orion)
warmpH <- phase1_formatted_orion%>%
  filter(treatment == "warm")
View(warmpH)
mean(warmpH$pH)

View(phase1_formatted_orion)
bothpH <- phase1_formatted_orion%>%
  filter(treatment == "both")
View(bothpH)
mean(bothpH$pH)

mean(controlpH$pH)
mean(hyppH$pH)
mean(warmpH$pH)
mean(bothpH$pH)



View(phase2_formatted_orion)
controlpH2 <- phase2_formatted_orion%>%
  filter(treatment == "control")
View(controlpH2)
mean(controlpH2$pH)

View(phase1_formatted_orion)
hyppH2 <- phase2_formatted_orion%>%
  filter(treatment == "hypoxic")
View(hyppH2)
mean(hyppH2$pH)

View(phase2_formatted_orion)
warmpH2 <- phase2_formatted_orion%>%
  filter(treatment == "warm")
View(warmpH2)
mean(warmpH2$pH)

View(phase2_formatted_orion)
bothpH2 <- phase2_formatted_orion%>%
  filter(treatment == "both")
View(bothpH2)
mean(bothpH2$pH)

mean(controlpH2$pH)
mean(hyppH2$pH)
mean(warmpH2$pH)
mean(bothpH2$pH)

#pH, Phase 2
w8 <- lmer(pH ~ DO_treat*temp_treat+(1|tank_num), data = phase2_formatted_orion, REML=TRUE)
Anova(w8, test="F", type="III")
emmeans(w8, pairwise ~ DO_treat*temp_treat, adjust = "none")


#Alk data, analyze
Alk <- read_csv("~/Desktop/MontagueORCC/ORION_wq_data/Alkalinity_cleaned_ORCC.csv")
colnames(Alk)

#Alkalinity, Phase 1
w9 <- lm(alkalinity_umolL ~ DO_treat*temp_treat, data = subset(Alk, Phase == "1"))
Anova(w9, test="F", type="III")

Alk$treatment <- factor(Alk$treatment, 
                        levels = c("Cont", "Warm","Hyp",  "Both"))
ggplot(subset(Alk, Phase == "1")) +
  aes(x = treatment, y = alkalinity_umolL, colour = treatment) +
  geom_point() +
  theme_classic()

#Alkalinity, Phase 2
w10 <- lm(alkalinity_umolL ~ DO_treat*temp_treat, data = subset(Alk, Phase == "2"))
Anova(w10, test="F", type="III")

ggplot(subset(Alk, Phase == "2")) +
  aes(x = treatment, y = alkalinity_umolL, colour = treatment) +
  geom_point() +
  theme_classic()

#pCO2
library(dplyr)
library(seacarb)

Alk_clean <- Alk %>%
  filter(
    !is.na(pH),
    !is.na(alkalinity_molkg),
    !is.na(salinity_psu),
    !is.na(temp_C),
    !is.na(pressure_barr))%>%
  rowwise() %>%
  mutate(pCO2uatm = carb(
    flag = 8,
    var1 = pH,
    var2 = (alkalinity_umolL/(10^6)),
    S = salinity_psu,
    T = temp_C,
    P = pressure_barr
  )$pCO2) %>%
  ungroup()
View(Alk_clean)

#Phase 1
w11 <- lm(pCO2uatm ~ DO_treat*temp_treat, data = subset(Alk_clean, Phase =="1"))
Anova(w11, test="F", type="III")

Alk_clean$treatment <- factor(Alk_clean$treatment, 
                        levels = c("Cont", "Warm","Hyp",  "Both"))
ggplot(subset(Alk_clean, Phase == "1")) +
  aes(x = treatment, y = pCO2uatm, colour = treatment) +
  geom_point() +
  theme_classic()

#pCO2, Phase 2
w12 <- lm(pCO2uatm ~ DO_treat*temp_treat, data = subset(Alk_clean, Phase=="2"))
Anova(w12, test="F", type="III")

ggplot(subset(Alk_clean, Phase == "2")) +
  aes(x = treatment, y = pCO2uatm, colour = treatment) +
  geom_point() +
  theme_classic()

carb(flag = 8, var1 = 7.96, var2 = 1402.67/(10^6),
              S = 9, T = 25.1, P = 1.018)




