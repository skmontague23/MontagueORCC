setwd("/Users/sophiemontague/Desktop/MontagueORCC/Year 2 Data")
getwd()

#read in data
library(readr)
df <- read_csv("Oysters to pull for Gonads.csv")
df <- read_csv("Oysters to pull for Spawning.csv") #Also assigning phase 1 information to this sheet

#adapted from JGM's code
#adding phase 1 treatment based on tag color

# white = control
# orange = both
# blue = hypoxic
# green = warm

for (i in 1:length(df$Tag_color)) { # loop through each row in df
  tag <- df$Tag_color[i] # read that row's tag color
  if (tag == 'B') { # if the tag color is blue, it's hypoxic
    df$Phase_1_treat[i] <- 'Hyp'
  } else if (tag == 'G') { # if the tag color is green, it's warm
    df$Phase_1_treat[i] <- 'Warm'
  } else if (tag == 'O') { # if the tag color is orange, it's both
    df$Phase_1_treat[i] <- 'Both'
  } else if (tag == 'W') { # if the tag color is white, it's control
    df$Phase_1_treat[i] <- 'Cont'
  }
}

# assigning tag numbers to phase 1 rep
rep1 <- c(1,2,13,14,25,26,37,38,49,50,61,62)
rep2 <- c(3,4,15,16,27,28,39,40,51,52,63,64)
rep3 <- c(5,6,17,18,29,30,41,42,53,54,65,66)
rep4 <- c(7,8,19,20,31,32,43,44,55,56,67,68)
rep5 <- c(9,10,21,22,33,34,45,46,57,58,69,70)
rep6 <- c(11,12,23,24,35,36,47,48,59,60,71,72)

# changing to numeric class
df$Tag_num <- as.numeric(df$Tag_num)

# initialize empty list to store unmatched tag numbers
unmatched_rows <- list()# assigning tag numbers to phase 1 rep
rep1 <- c(1,2,13,14,25,26,37,38,49,50,61,62)
rep2 <- c(3,4,15,16,27,28,39,40,51,52,63,64)
rep3 <- c(5,6,17,18,29,30,41,42,53,54,65,66)
rep4 <- c(7,8,19,20,31,32,43,44,55,56,67,68)
rep5 <- c(9,10,21,22,33,34,45,46,57,58,69,70)
rep6 <- c(11,12,23,24,35,36,47,48,59,60,71,72)

# changing to numeric class
df$Tag_num <- as.numeric(df$Tag_num)

# initialize empty list to store unmatched tag numbers
unmatched_rows <- list()

# convert list to df
unmatched_df <- do.call(rbind, unmatched_rows)
# print data frame of odd tag numbers
unmatched_df

# for loop to add phase 1 rep info
for (i in 1:length(df$Tag_num)) { # loop through each row in df
  tag <- df$Tag_num[i] # read the tag number 
  if (tag %in% rep1) { # if tag number is in the rep1 list, the rep is 01
    df$Phase_1_rep[i] <- '01'
  } else if (tag %in% rep2) {
    df$Phase_1_rep[i] <- '02'
  } else if (tag %in% rep3) {
    df$Phase_1_rep[i] <- '03'
  } else if (tag %in% rep4) {
    df$Phase_1_rep[i] <- '04'
  } else if (tag %in% rep5) {
    df$Phase_1_rep[i] <- '05'
  } else if (tag %in% rep6) {
    df$Phase_1_rep[i] <- '06'
  } else { # if there are any tag numbers that don't match 
    unmatched_rows[[length(unmatched_rows) + 1]] <- df[i, ] # add to list
  }
}

View(df)
# then write csv file
write.csv(df, '/Users/sophiemontague/Desktop/MontagueORCC/Year 2 Data/Oysters to pull for Spawning1.csv')


#check if the phase 1 treatment and reps are even

sum(df$Phase_1_treat == "Cont")
sum(df$Phase_1_treat == "Cont" & df$Phase_1_rep == "01")
sum(df$Phase_1_treat == "Cont" & df$Phase_1_rep == "02")
sum(df$Phase_1_treat == "Cont" & df$Phase_1_rep == "03")
sum(df$Phase_1_treat == "Cont" & df$Phase_1_rep == "04")
sum(df$Phase_1_treat == "Cont" & df$Phase_1_rep == "05")
sum(df$Phase_1_treat == "Cont" & df$Phase_1_rep == "06")
sum(df$Phase_1_treat == "Warm")
sum(df$Phase_1_treat == "Warm" & df$Phase_1_rep == "01")
sum(df$Phase_1_treat == "Warm" & df$Phase_1_rep == "02")
sum(df$Phase_1_treat == "Warm" & df$Phase_1_rep == "03")
sum(df$Phase_1_treat == "Warm" & df$Phase_1_rep == "04")
sum(df$Phase_1_treat == "Warm" & df$Phase_1_rep == "05")
sum(df$Phase_1_treat == "Warm" & df$Phase_1_rep == "06")
sum(df$Phase_1_treat == "Hyp")
sum(df$Phase_1_treat == "Hyp" & df$Phase_1_rep == "01")
sum(df$Phase_1_treat == "Hyp" & df$Phase_1_rep == "02")
sum(df$Phase_1_treat == "Hyp" & df$Phase_1_rep == "03")
sum(df$Phase_1_treat == "Hyp" & df$Phase_1_rep == "04")
sum(df$Phase_1_treat == "Hyp" & df$Phase_1_rep == "05")
sum(df$Phase_1_treat == "Hyp" & df$Phase_1_rep == "06")
sum(df$Phase_1_treat == "Both")

sum(df$Phase_1_rep == "01")
as.data.frame(df$Phase_1_rep == "01")
sum(df$Phase_1_rep == "02")
sum(df$Phase_1_rep == "03")
sum(df$Phase_1_rep == "04")
sum(df$Phase_1_rep == "05")
sum(df$Phase_1_rep == "06")

#one that got put in rep 1 was supposed to go in rep 5. for gonads

sum(df$Phase_2_treat == "Cont")
sum(df$Phase_2_treat == "Warm")
sum(df$Phase_2_treat == "Hyp")
sum(df$Phase_2_treat == "Both")


