library(esquisse)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readxl)
library(lubridate)
library(plotrix)
library(ggplot2)
install.packages("lubridate")

#Wellfleet 1
Wellfleet1 <- read_excel("Wellfleet1 2025-03-21 16_11_07 EDT (Data EDT).xlsx")
View(Wellfleet1_high)

Wellfleet1_high <- Wellfleet1 %>%
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-12-13 12:00:00', tz = "UTC") & 
           `Date-Time (EST/EDT)` <= as.POSIXct('2025-03-10 12:00:00', tz = "UTC"))%>%
  mutate(site = "Wellfleet 1, High Pit")

ggplot(Wellfleet1_high) +
  aes(x = `Date-Time (EST/EDT)`, y = `Temperature , °C`) +
  geom_line(colour = "#112446") +
  theme_classic() +
  labs(title = "Wellfleet 1, High")+
  ylim(0, 6.5)

#Wellfleet 2
Wellfleet2 <- read_excel("Wellfleet2 2025-03-21 16_10_07 EDT (Data EDT).xlsx")
View(Wellfleet2_low)

Wellfleet2_low <- Wellfleet2 %>%
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-12-13 12:00:00', tz = "UTC") & 
         `Date-Time (EST/EDT)` <= as.POSIXct('2025-03-10 12:00:00', tz = "UTC"))%>%
  mutate(site = "Wellfleet 2, Low Pit")


ggplot(Wellfleet2_low) +
  aes(x = `Date-Time (EST/EDT)`, y = `Temperature , °C`) +
  geom_line(colour = "#112446") +
  theme_classic() +
  labs(title = "Wellfleet 2, Low")+
  ylim(0, 6.5)

library(dplyr)
library(lubridate)



# View filtered data
filtered_data



range(Wellfleet2_low$`Temperature , °C`)
range(Wellfleet2_low$`Date-Time (EST/EDT)`)


#### West Island ####
library(readxl)
library(dplyr)
library(ggplot2)
library(lubridate)

# Read and filter with site label
Wisland_surface1 <- read_excel("WI_S1 W isle surf 1 2025-04-11 13_21_40 EDT.xlsx") %>%
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-11-15 12:00:00', tz = "UTC") & 
           `Date-Time (EST/EDT)` <= as.POSIXct('2025-04-09 12:00:00', tz = "UTC")) %>%
  mutate(site = "West Island Surface 1")

Wisland_surface2 <- read_excel("Wellfleet2 2025-03-21 16_10_07 EDT (Data EDT).xlsx") %>%
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-11-15 12:00:00', tz = "UTC") & 
           `Date-Time (EST/EDT)` <= as.POSIXct('2025-04-09 12:00:00', tz = "UTC")) %>%
  mutate(site = "West Island Surface 2")

Wisland_bottom1 <- read_excel("WI_B1 W isle bottom1 2025-04-11 13_22_39 EDT.xlsx") %>%
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-11-15 12:00:00', tz = "UTC") & 
           `Date-Time (EST/EDT)` <= as.POSIXct('2025-04-09 12:00:00', tz = "UTC")) %>%
  mutate(site = "West Island Bottom 1")

Wisland_bottom2 <- read_excel("WI_B2 W isle bottom2 2025-04-11 13_24_33 EDT.xlsx") %>%
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-11-15 12:00:00', tz = "UTC") & 
           `Date-Time (EST/EDT)` <= as.POSIXct('2025-04-09 12:00:00', tz = "UTC")) %>%
  mutate(site = "West Island Bottom 2")


all_data <- bind_rows(Wisland_surface1, Wisland_bottom1, Wisland_bottom2, Wellfleet1_high, Wellfleet2_low)

ggplot(all_data, aes(x = `Date-Time (EST/EDT)`, y = `Temperature , °C`, color = site)) +
  geom_line() +
  theme_classic() +
  labs(
       x = "Date-Time (EST/EDT)",
       y = "Temperature (°C)",
       color = "Site") +
  xlim(as.POSIXct('2024-11-15 12:00:00', tz = "UTC"),
       as.POSIXct('2025-04-09 12:00:00', tz = "UTC"))

WI_data <- bind_rows(Wisland_surface1, Wisland_bottom1, Wisland_bottom2)

ggplot(WI_data, aes(x = `Date-Time (EST/EDT)`, y = `Temperature , °C`, color = site)) +
  geom_line() +
  theme_classic() +
  labs(
    x = "Date-Time (EST/EDT)",
    y = "Temperature (°C)",
    color = "Site") +
  xlim(as.POSIXct('2024-12-15 12:00:00', tz = "UTC"),
       as.POSIXct('2025-04-01 12:00:00', tz = "UTC"))

Wisland_surface1 <- read_excel("WI_S1 W isle surf 1 2025-04-11 13_21_40 EDT.xlsx")
Wisland_surface1%>% 
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-11-15 12:00:00', tz = "UTC") & 
         `Date-Time (EST/EDT)` <= as.POSIXct('2025-04-09 12:00:00', tz = "UTC"))%>%
  ggplot() +
  aes(x = `Date-Time (EST/EDT)`, y = `Temperature , °C`) +
  geom_line(colour = "#112446") +
  theme_classic() +
  labs(title = "West Island Surface 1")

Wisland_surface2 <- read_excel("Wellfleet2 2025-03-21 16_10_07 EDT (Data EDT).xlsx")
str(Wisland_surface2$`Date-Time (EST/EDT)`)

Wisland_surface2%>%
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-11-15 12:00:00', tz = "UTC") & 
           `Date-Time (EST/EDT)` <= as.POSIXct('2025-04-09 12:00:00', tz = "UTC"))%>%
  ggplot() +
  aes(x = `Date-Time (EST/EDT)`, y = `Temperature , °C`) +
  geom_line(colour = "#112446") +
  theme_classic() +
  labs(title = "West Island Surface 2")

Wisland_bottom1 <- read_excel("WI_B1 W isle bottom1 2025-04-11 13_22_39 EDT.xlsx")
Wisland_bottom1%>%
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-11-15 12:00:00', tz = "UTC") & 
           `Date-Time (EST/EDT)` <= as.POSIXct('2025-04-09 12:00:00', tz = "UTC"))%>%
  ggplot() +
  aes(x = `Date-Time (EST/EDT)`, y = `Temperature , °C`) +
  geom_line(colour = "#112446") +
  theme_classic() +
  labs(title = "West Island Bottom 1")

Wisland_bottom2 <- read_excel("WI_B2 W isle bottom2 2025-04-11 13_24_33 EDT.xlsx")
Wisland_bottom2%>%
  filter(`Date-Time (EST/EDT)` >= as.POSIXct('2024-11-15 12:00:00', tz = "UTC") & 
         `Date-Time (EST/EDT)` <= as.POSIXct('2025-04-09 12:00:00', tz = "UTC"))%>%
  ggplot() +
  aes(x = `Date-Time (EST/EDT)`, y = `Temperature , °C`) +
  geom_line(colour = "#112446") +
  theme_classic() +
  labs(title = "West Island Bottom 2")

colnames(Wisland_bottom2)




#Wellfleet summary stats
colnames(Wellfleet1_high)
View(Wellfleet1_high)

  #add month column
Wellfleet1_high <- Wellfleet1_high %>%
  mutate(month = month(`Date-Time (EST/EDT)`))

summary_stats_Wellfleet <- Wellfleet1_high %>%
  group_by(month) %>%
  summarise(
    mean_temp = mean(`Temperature , °C`, na.rm = TRUE),
    se_temp = sd(`Temperature , °C`, na.rm = TRUE),
    .groups = 'drop')

Wellfleet1_high <- Wellfleet1_high %>%
  mutate(month = month(`Date-Time (EST/EDT)`, label = TRUE, abbr = TRUE))

summary_stats_Wellfleet$month <- factor(summary_stats_Wellfleet$month, levels=c("Dec", "Jan", "Feb", "Mar"))

  #plot
ggplot(summary_stats_Wellfleet, aes(x = month, y = mean_temp, color = month)) +
  geom_point(size = 4, position = position_dodge(0.9)) +
  geom_errorbar(aes(ymin = mean_temp - se_temp, ymax = mean_temp + se_temp),
                width = 0.2, position = position_dodge(0.9)) +
  theme_classic(base_size = 20) +
  labs(x = "Month", y = "Temperature (°C)") +
  ylim(0,6)+
  theme(legend.position = "none")

View(summary_stats_Wellfleet)
