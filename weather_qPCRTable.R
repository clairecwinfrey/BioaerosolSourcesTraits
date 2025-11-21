# Weather_qPCRTable
# November 25, 2024, edited November 20, 2025

# Script description: The purpose of this script is to create weather and qPCR plots for a supplemental table,
# which I assume will be the first piece of supplemental information created
##################
# I. SET-UP 
##################
# 1. LIBRARIES
library(vegan)
library(phyloseq)
library(tidyr)
library(dplyr)
library(purrr)
library(patchwork)
library(lubridate)
library(stringr)
sessionInfo()

list.files()

# 2. Bring in important data and files 
## qPCR data made in habitatAnalyses.R
I6SqPCRdata_forTable <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/I6SqPCRdata_forTable_Nov17-2025.csv") 
head(I6SqPCRdata_forTable)
ITSqPCRdata_forTable <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/ITSqPCRdata_forTable_Nov17-2025.csv") 
head(ITSqPCRdata_forTable)
## (Most) Weather and meterological info
SRS_metDatAll <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/SRS_metDatAll.csv") #bring in the SRS meterological data from Brian Viner. Note that this is 4 hours
# ahead of EDT time, as UTC
head(SRS_metDatAll)
colnames(SRS_metDatAll)
windDatLogger <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/windDataLoggerApril23_2023.csv") # data from the wind data logger
head(windDatLogger)
SRS_humidity <- read.csv("~/Desktop/CU_Research/SRS_Aeromicrobiome/BioinformaticsAndMetadata/SRS_humidity.csv") #bring in the SRS humidity data from Brian Viner. Note that this is in EDT
head(SRS_humidity)
colnames(SRS_humidity)[2] <- "relHumidity_2m" #change humidity column name to remove X at beginning
head(SRS_humidity)
## Phyloseq objects FOR JUST AIR DATA 
airOnly_16Sr_noAirSingsDoubs_phyloseq <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/air16S_noSingDoubs_ps_sept25.rds")
airOnly_ITS_noAirSingsDoubs_phyloseq <- readRDS(file="~/Desktop/CU_Research/SRS_Aeromicrobiome/scriptsDoubleCheck/RobjectsToReCheck/airOnly_ITS_noAirSingsDoubs_ps.rds")

##################
# II CLEAN UP MET DATA AND MERGE WITH WIND DATA LOGGER DATA
##################
# Tackle the SRS metadata first. 
head(SRS_metDatAll)
str(SRS_metDatAll)
colnames(SRS_metDatAll)
# Since "windDatLogger" has wind speed, gust, and direction, will only keep columns pertaining to
# rainfall and temperature at 2m.
SRS_metDatAll_smaller <- SRS_metDatAll %>% 
  select(-Climo.18m.Wind.Direction..degrees.CW.from.North.) %>% 
  select(-Climo.36m.Wind.Direction..degrees.CW.from.North.) %>% 
  select(-Climo.18m.Temperature...C.) %>% 
  select(-Climo.36m.Temperature...C.) %>% 
  select(-Climo.61m.Temperature...C.) %>% 
  select(-Climo.18m.Wind.Speed..m.s.) %>% 
  select(-Climo.36m.Wind.Speed..m.s.) %>% 
  select(-Climo.61m.Wind.Speed..m.s.) %>% 
  select(-A.Area.61m.Wind.Direction..degrees.CW.from.North.) %>% 
  select(-A.Area.61m.Temperature...C.) %>% 
  select(-A.Area.61m.Wind.Speed..m.s.) %>% 
  select(-H.Area.61m.Wind.Direction..degrees.CW.from.North.) %>% 
  select(-H.Area.61m.Temperature...C.) %>% 
  select(-H.Area.61m.Wind.Speed..m.s.) %>% 
  select(-P.Area.61m.Wind.Direction..degrees.CW.from.North.) %>% 
  select(-P.Area.61m.Temperature...C.) %>% 
  select(-P.Area.61m.Wind.Speed..m.s.) %>% 
  select(-Climo.61m.Wind.Direction..degrees.CW.from.North.) %>% 
  separate(Timestamp..UTC., into = c("Date", "time"), sep = " ") #break up Timestamp... UTC column into date and time
colnames(SRS_metDatAll_smaller) #has what we want!
# View(SRS_metDatAll_smaller)
# Remove weird extra empty rows that were imported from Excel
SRS_metDatAll_smaller <- SRS_metDatAll_smaller[-which(is.na(SRS_metDatAll_smaller$time)),]
# View(SRS_metDatAll_smaller) #looks good!

# Now need to make the date as it is in some of the other datasets. So first, check the format
head(windDatLogger) #e.g. 16-Jun-22
unique(as.data.frame(as.matrix(sample_data(airOnly_ITS_noAirSingsDoubs_phyloseq)))$DateSetOut) #e.g. 16-Jun-22

# # Edit these dates to match the formatting above
# unique(SRS_metDatAll_smaller$Date)
# unique(windDatLogger$Date)
# SRS_metDatAll_smaller <- SRS_metDatAll_smaller %>%
#   mutate(
#     date = case_when(
#       Date == "6/16/22" ~ "16-Jun-22",
#       Date == "6/17/22" ~ "17-Jun-22",
#       Date == "6/18/22" ~ "18-Jun-22",
#       Date == "6/19/22" ~ "19-Jun-22",
#       Date == "6/20/22" ~ "20-Jun-22",
#       Date == "6/21/22" ~ "21-Jun-22",
#       Date == "6/22/22" ~ "22-Jun-22",
#       Date == "6/23/22" ~ "23-Jun-22",
#       Date == "6/24/22" ~ "24-Jun-22",
#       Date == "6/25/22" ~ "25-Jun-22",
#       Date == "6/26/22" ~ "26-Jun-22",
#       Date == "6/27/22" ~ "27-Jun-22",
#       Date == "6/28/22" ~ "28-Jun-22",
#       Date == "6/29/22" ~ "29-Jun-22",
#       Date == "6/30/22" ~ "30-Jun-22",
#       Date == "7/1/22" ~ "1-Jul-22",
#       Date == "7/2/22" ~ "2-Jul-22",
#       Date == "7/3/22" ~ "3-Jul-22",
#       Date == "7/4/22" ~ "4-Jul-22",
#       Date == "7/5/22" ~ "5-Jul-22",
#       Date == "7/6/22" ~ "6-Jul-22",
#       Date == "7/7/22" ~ "7-Jul-22",
#       Date == "7/8/22" ~ "8-Jul-22" ,
#       Date == "7/9/22" ~ "9-Jul-22",
#       TRUE ~ NA_character_
#     )
#   )
# Clean up some of the column names
colnames(SRS_metDatAll_smaller)
colnames(SRS_metDatAll_smaller) <- c(colnames(SRS_metDatAll_smaller)[1], colnames(SRS_metDatAll_smaller)[2],
                                     "climo2mWind", "inchesRain15min", "climo2mTemp", "climo2mWSpeed",
                                     "PareaInchesRain15min")
colnames(SRS_metDatAll_smaller)

# Change time to be in EDT, not UTC. UTC is 4 hours ahead of EDT, so need to subtract 4 hours
# and occassionally change date too.. this is challenging since times are characters and have colons
# in them, so can't be made into numeric. 
SRS_metDatAll_smaller$time
# To avoid this issue, and since I don't feel like using complicated time/date tools in R...
# Will Separate times into hours and minutes as new columns, make hour column numeric and subtract 4
# hours when necessary, then convert back to character and re-merge with original time
colnames(SRS_metDatAll_smaller)
# i. Make hours separate and numeric
SRS_metDatAll_smaller <- SRS_metDatAll_smaller %>% 
  separate(time, into = c("hour", "minute"), sep = ":") #break up time into hours and minutes
str(SRS_metDatAll_smaller)
SRS_metDatAll_smaller$hour <- as.numeric(SRS_metDatAll_smaller$hour) #make hours numeric
str(SRS_metDatAll_smaller)

# ii. If hours are between 4 and 23, simply subtract 4 hours and can keep the date the same too
unique(SRS_metDatAll_smaller$hour) 
SRS_metDatAll_smaller <- SRS_metDatAll_smaller %>% 
  mutate(
    hour_fixed = case_when(
      hour >=4 ~ hour - 4, #if greater than or equal to 4, subtract four
      hour == 0 ~ 20, #midnight is 8p, 20:00
      hour == 1 ~ 21, #1 am is 9pm, 21:00
      hour == 2 ~ 22,  #2 am is 10pm, 22:00
      hour == 3 ~ 23, #3 am is 11pm, 23:00
    )
  )
# check out a few to confirm that this worked as expected:
cbind(SRS_metDatAll_smaller$hour, SRS_metDatAll_smaller$hour_fixed) #looks correct!

# iii. Change date to be the one prior for those with hour_fixed between 20 and 23
unique(SRS_metDatAll_smaller$Date)
# First, do a couple of checks to make sure that converting works as expected -- it does!
dateCheck <- as.Date(SRS_metDatAll_smaller$Date, "%m/%d/%y")
unique(dateCheck)[16] #get first day in July and see if substracting one yields June 30th!
unique(dateCheck)[16] - 1 #It does! "2022-06-30"
# Make new, date column
SRS_metDatAll_smaller$DateFormatted <- as.Date(SRS_metDatAll_smaller$Date, "%m/%d/%y") # make a new date column
str(SRS_metDatAll_smaller$DateFormatted)
cbind.data.frame(SRS_metDatAll_smaller$DateFormatted, SRS_metDatAll_smaller$Date) #looks correct
# Make a new column such that, if hour is 20-23, then DateFixed becomes the day prior
SRS_metDatAll_smaller <- SRS_metDatAll_smaller %>% 
  mutate(
    DateFixed = case_when(
      hour >=20 ~ DateFormatted - 1, #if hour is 20, 21, 22, or 23, take previous day as date
      hour < 20 ~ DateFormatted #if hour less than 20 keep it as it is! 
    )
  )
# Check a few to make sure that this worked as expected:
cbind.data.frame(SRS_metDatAll_smaller$hour, SRS_metDatAll_smaller$DateFormatted, SRS_metDatAll_smaller$DateFixed) #looks correct

# Finally, clean all of this up so that can move on with the rest of the script, making a new object called
# "SRS_metDatAll_fixed"
# View(SRS_metDatAll_smaller)
head(SRS_metDatAll_smaller)
# i. Clean up time
SRS_metDatAll_fixed <- SRS_metDatAll_smaller %>% 
  select(-Date) %>% #remove "Date" because it's wrong
  select(-hour) %>% #remove "hour" because it's wrong
  select(-DateFormatted) %>% #removed DateFormatted
  tidyr::unite(col=time, c("hour_fixed", "minute"), sep=":") #make one new time!
head(SRS_metDatAll_fixed)
# ii. Re-name "DateFixed" to make it "Date", as in windLoggerData
colnames(SRS_metDatAll_fixed)
colnames(SRS_metDatAll_fixed)[7] <- "Date"
colnames(SRS_metDatAll_fixed)

# Get time of day to be the same in SRS_metDatAll_fixed and windDatLogger
# These overall look the same, except windDatLogger has seconds as well, which aren't necessary. 
SRS_metDatAll_fixed$time
windDatLogger$Time
# View(windDatLogger)
# Remove the seconds from the windDatLogger dataframe to make it match SRS_metDatAll_fixed,
# and rename it "time" lowercase for merging.
windDatLogger <- windDatLogger %>% 
  mutate(time = sub(":00$", "", Time))
windDatLogger$time
class(windDatLogger$time)
SRS_metDatAll_fixed$time

# Remove date formatting in SRS_metDatAll because windDatLogger is not easily made to a formatted date
SRS_metDatAll_fixed$Date <- as.character(SRS_metDatAll_fixed$Date) #make character
class(SRS_metDatAll_fixed$Date)
unique(SRS_metDatAll_fixed$Date)

# Make Date formatting match, specifically make both look like windDatLogger$Date
unique(windDatLogger$Date) #e.g. "16-Jun-22"
unique(SRS_metDatAll_fixed$Date) #e.g. "2022-07-03"
length(unique(windDatLogger$Date)) == length(unique(SRS_metDatAll_fixed$Date)) 
# This is longer because SRS_metDatAll_fixed$Date has 2022-06-15, which it doesn't need
# Remove the 2022-06-15
SRS_metDatAll_fixed <- SRS_metDatAll_fixed[-which(SRS_metDatAll_fixed$Date %in% "2022-06-15"),]
length(unique(windDatLogger$Date)) == length(unique(SRS_metDatAll_fixed$Date)) 

# i. Make format the same, starting with SRS_metDatAll_fixed
SRS_metDatAll_fixed <- SRS_metDatAll_fixed %>% 
  mutate(
    DateFixed = ymd(Date) |> format("%d-%b-%y")
  )
unique(SRS_metDatAll_fixed$DateFixed)
# Spot check old version with new version, looks good!
cbind(SRS_metDatAll_fixed$Date, SRS_metDatAll_fixed$DateFixed)
# ii. Change windDatLogger$Date too to have 0s before single-digit day of month
windDatLogger <- windDatLogger %>% 
  mutate(
    Date = dmy(Date) |> format("%d-%b-%y")
  )
unique(windDatLogger$Date)
# View(SRS_metDatAll_fixed)
# Double check that these are now formatted the same:
unique(windDatLogger$Date) %in% unique(SRS_metDatAll_fixed$DateFixed)
# Shows no incorrectly formatted dates
unique(SRS_metDatAll_fixed$DateFixed)[which(unique(SRS_metDatAll_fixed$DateFixed) %in% unique(windDatLogger$Date) == FALSE)]
head(windDatLogger)
head(SRS_metDatAll_fixed)
# Remove "Date" from SRS_metDatAll_fixed
SRS_metDatAll_fixed <- SRS_metDatAll_fixed %>% 
  select(-Date)
head(SRS_metDatAll_fixed)
# Rename "DateFixed" as "Date" to be able to merge with windDatLogger
colnames(SRS_metDatAll_fixed)[which(colnames(SRS_metDatAll_fixed) %in% "DateFixed")] <- "Date"
head(SRS_metDatAll_fixed)

# Merge these two data frames, keeping only the relevant times from windDatLogger
head(windDatLogger)
head(SRS_metDatAll_fixed)
allMeterologicalData_1 <- left_join(windDatLogger, SRS_metDatAll_fixed, by = c("Date", "time")) #"time" in both is formatted without seconds
dim(SRS_metDatAll_fixed); dim(windDatLogger)
setdiff(allMeterologicalData_1$Date, windDatLogger$Date) # these are the same!
nrow(windDatLogger) == nrow(allMeterologicalData_1) #checks that the same number of rows was preserved from windDatLogger
allMeterologicalData_1 <- allMeterologicalData_1 %>% 
  select(-Time) #remove "Time" since it is formatted wrong and to avoid confusion later

# Check it all out!
# View(allMeterologicalData_1)

##################
# III ADD IN HUMIDITY
##################
# Lines below show how SRS_humidity needs to be fixed before merging: 1)  "Timestamp._EDT" needs to be split into "Date" and "time",
# 2) Date needs to be reformatted to be like 16-Jun-22 instead of 6/16/22. 
colnames(SRS_humidity)
head(SRS_humidity)
colnames(allMeterologicalData_1); head(allMeterologicalData_1)

# 1. Reformat SRS_humidity as specified above
unique(allMeterologicalData_1$Date)
colnames(SRS_humidity)
SRS_humidity_2 <- SRS_humidity %>% 
  separate(Timestamp._EDT, into = c("Date", "time"), sep = " ") %>% 
  mutate(
    DateFixed = mdy(Date) |> format("%d-%b-%y") #Get a new, correctly formatted date column
  ) 
# View(SRS_humidity_2) #looks good!
head(SRS_humidity_2) 
# Because DateFixed looks correct, can drop "Date" then rename "DateFixed" "Date" for merging below:
SRS_humidity_2 <- SRS_humidity_2 %>% 
  select(-Date)
head(SRS_humidity_2)
colnames(SRS_humidity_2)[colnames(SRS_humidity_2) %in% "DateFixed"] <- "Date" #rename column

# 2. Merge with allMeterologicalData_1
colnames(SRS_humidity_2)
colnames(allMeterologicalData_1)
allMeterologicalData_2 <- left_join(allMeterologicalData_1, SRS_humidity_2, by = c("Date", "time")) #"time" in both is formatted without seconds
nrow(windDatLogger) == nrow(allMeterologicalData_2) #checks that the same number of rows was preserved from windDatLogger
setdiff(allMeterologicalData_2$Date, windDatLogger$Date) # these are the same!
head(allMeterologicalData_2)

##################
# IV. GET INFO ON EACH SAMPLING PERIOD
##################
# Now, because nearly all sampling days were unevenly split between two days, need to combine all of the 15-min data for 
# each 24 hour sampling period and then get means of this.

# 1. Get metadata for 16S and ITS data
I6S_airMeta <- as.data.frame(as.matrix(sample_data(airOnly_16Sr_noAirSingsDoubs_phyloseq)))
# View(I6S_airMeta)
ITS_airMeta <- as.data.frame(as.matrix(sample_data(airOnly_ITS_noAirSingsDoubs_phyloseq)))
# View(ITS_airMeta)
# The metadata for each shows that "daysOut" is a unique identifier for each 24-hour sampling period. It also shows that the 
# dates are formatted differently from the meterological data, which will need to be updated.

# 2. Fix dates in the I6S and ITS metadata to make 2022 just 22, as in the meteorological data
I6S_airMeta <- I6S_airMeta %>%
  mutate(across(everything(), ~ gsub("2022", "22", .)))
ITS_airMeta <- ITS_airMeta %>%
  mutate(across(everything(), ~ gsub("2022", "22", .)))

# 3. Remove weird spaces in daysOut
ITS_airMeta$daysOut <- gsub(x=ITS_airMeta$daysOut, pattern = " ", replacement = "")
I6S_airMeta$daysOut <- gsub(x=I6S_airMeta$daysOut, pattern = " ", replacement = "")

# 3. For each unique sampling day "daysOut" make a vector containing "DateSetOut" and "DateRetrieved"
## i. Pre-allocate a matrix, where each row contains a 2 string vector for "DateSetOut" and "DateRetrieved" for each day
I6S_datesPerSampPeriod <- matrix(nrow=length(unique(I6S_airMeta$daysOut)), ncol = 4)
colnames(I6S_datesPerSampPeriod) <- c("daysOut", "DateSetOut", "DateRetrieved", "EU")
I6S_datesPerSampPeriod

ITS_datesPerSampPeriod <- matrix(nrow=length(unique(ITS_airMeta$daysOut)), ncol = 4)
colnames(ITS_datesPerSampPeriod) <- c("daysOut", "DateSetOut", "DateRetrieved", "EU")
ITS_datesPerSampPeriod

# 4. Fix a few of the dates
class(I6S_airMeta$DateRetrieved)
I6S_airMeta <- I6S_airMeta %>% 
mutate(
  DateRetrieved = dmy(DateRetrieved) |> format("%d-%b-%y"), #Get a new, correctly formatted date column
  DateSetOut = dmy(DateSetOut) |> format("%d-%b-%y")
) 

ITS_airMeta <- ITS_airMeta %>% 
  mutate(
    DateRetrieved = dmy(DateRetrieved) |> format("%d-%b-%y"), #Get a new, correctly formatted date column
    DateSetOut = dmy(DateSetOut) |> format("%d-%b-%y")
  ) 

# ii. Write for loops which grab each unique  "DateSetOut", "DateRetrieved", and EU for each "daysOut"
for (i in 1:length(unique(I6S_airMeta$daysOut))){
  I6S_datesPerSampPeriod[i,1] <- unique(I6S_airMeta$daysOut)[i] #in each ith row of column 1 put the ith unique daysOut
  # Pulls out which rows in the whole I6S_airMeta match the ith unique daysOut, isolating the DateSetOut. In other words, this code
  # # finds out how each unique daysOut maps to DateSetOut. Puts this in DateSetOut column (column 2)
  I6S_datesPerSampPeriod[i,2] <- unique(I6S_airMeta$DateSetOut[which(I6S_airMeta$daysOut == unique(I6S_airMeta$daysOut)[i])])
  # Does the same as previous line but pulls out DateRetrieved
  I6S_datesPerSampPeriod[i,3] <- unique(I6S_airMeta$DateRetrieved[which(I6S_airMeta$daysOut == unique(I6S_airMeta$daysOut)[i])])
  # Does the same as previous 2 lines but pulls out EU
  I6S_datesPerSampPeriod[i,4] <- unique(I6S_airMeta$EU[which(I6S_airMeta$daysOut == unique(I6S_airMeta$daysOut)[i])])
}

# I6S_datesPerSampPeriod #checked a few manually against I6S_airMeta and they look good!

for (i in 1:length(unique(ITS_airMeta$daysOut))){
  ITS_datesPerSampPeriod[i,1] <- unique(ITS_airMeta$daysOut)[i]
  ITS_datesPerSampPeriod[i,2] <- unique(ITS_airMeta$DateSetOut[which(ITS_airMeta$daysOut == unique(ITS_airMeta$daysOut)[i])])
  ITS_datesPerSampPeriod[i,3] <- unique(ITS_airMeta$DateRetrieved[which(ITS_airMeta$daysOut == unique(ITS_airMeta$daysOut)[i])])
  ITS_datesPerSampPeriod[i,4] <- unique(ITS_airMeta$EU[which(ITS_airMeta$daysOut == unique(ITS_airMeta$daysOut)[i])])
}

ITS_datesPerSampPeriod #checked a few manually against ITS_airMeta and they look good!

# iii. Check for loops
# Tests with 2nd unique day in I6S to show that it's indexing correctly.
unique(I6S_airMeta$daysOut)[2]
I6S_airMeta$daysOut[which(I6S_airMeta$daysOut == unique(I6S_airMeta$daysOut)[2])]
I6S_airMeta$daysOut[which(I6S_airMeta$daysOut == unique(I6S_airMeta$daysOut)[2])] == unique(I6S_airMeta$daysOut)[2]
unique(I6S_airMeta$DateSetOut[which(I6S_airMeta$daysOut == unique(I6S_airMeta$daysOut)[2])])
unique(I6S_airMeta$DateRetrieved[which(I6S_airMeta$daysOut == unique(I6S_airMeta$daysOut)[2])])

# 4. For 16S, make a list, named dayOut and EU, where each element has all of the data from allMeterologicalData_2 for all of the days
I6S_meteroDatByDay <- vector("list", length=nrow(I6S_datesPerSampPeriod))
# Name each element a combination of EU and daysout
for (i in 1:length(I6S_meteroDatByDay)){
  names(I6S_meteroDatByDay)[[i]] <- paste0(I6S_datesPerSampPeriod[i,4],"_",I6S_datesPerSampPeriod[i,1])
}
names(I6S_meteroDatByDay)
colnames(allMeterologicalData_2)
head(allMeterologicalData_2)
# For loop to separate out allMeterologicalData_2 by each sampling period
for (i in 1:length(I6S_meteroDatByDay)){
  I6S_meteroDatByDay[[i]] <- allMeterologicalData_2 %>%
    # Filter out allMeterologicalData2 by each ith in I6S_datesPerSampPeriod DateSetOut and DateRetrieved AND EU (because some days I sampled two EUs)
    filter(Date %in% c(I6S_datesPerSampPeriod[i,2], I6S_datesPerSampPeriod[i,3]) & EU %in% I6S_datesPerSampPeriod[i,4])
}
# QUICK CHECKS -- LOOKS GOOD!
# View(I6S_meteroDatByDay[[1]])
# View(I6S_meteroDatByDay[[16]])
# View(I6S_meteroDatByDay[[12]]) #this is last day out, day 23 (looks correct when checked against windDataLoggerApril23_2023EXCEL)

# Double check that classes are correct for everything-- they are not.
colnames(I6S_meteroDatByDay[[1]])
str(I6S_meteroDatByDay)
# Get the names of variables that I want to make numeric.
colnames(I6S_meteroDatByDay[[1]])
colsToNumeric <- c("Speed", "Gust", "avgSpeedSampPeriod", "climo2mWind",
                     "inchesRain15min", "climo2mTemp", "climo2mWSpeed", "PareaInchesRain15min", "relHumidity_2m")
# Make all these columns numeric across the list
I6S_meteroDatByDay <- map(I6S_meteroDatByDay, ~ .x %>%
                           mutate(across(all_of(colsToNumeric), ~ as.numeric(.)))
)
str(I6S_meteroDatByDay)

# 5. For ITS, make a list, named dayOut and EU, where each element has all of the data from allMeterologicalData_2 for all of the days
ITS_meteroDatByDay <- vector("list", length=nrow(ITS_datesPerSampPeriod)) 
# Name each element a combination of EU and daysout
for (i in 1:length(ITS_meteroDatByDay)){
  names(ITS_meteroDatByDay)[[i]] <- paste0(ITS_datesPerSampPeriod[i,4],"_",ITS_datesPerSampPeriod[i,1])
}
names(ITS_meteroDatByDay)
colnames(allMeterologicalData_2)
head(allMeterologicalData_2)

for (i in 1:length(ITS_meteroDatByDay)){
  ITS_meteroDatByDay[[i]] <- allMeterologicalData_2 %>% 
    filter(Date %in% c(ITS_datesPerSampPeriod[i,2], ITS_datesPerSampPeriod[i,3]) & EU %in% ITS_datesPerSampPeriod[i,4])
}
# LOOKS GOOD!
# View(ITS_meteroDatByDay[[1]])
# View(ITS_meteroDatByDay[[2]])

# Make all these columns numeric across the list
ITS_meteroDatByDay <- map(ITS_meteroDatByDay, ~ .x %>%
                            mutate(across(all_of(colsToNumeric), ~ as.numeric(.)))
)
str(ITS_meteroDatByDay)

##################
# V. WITHIN EACH OF THE DAYS, GET MEAN VALUES FOR TEMPERATURE, WIND SPEED, HIGH, LOW, AND MEAN TEMP,
# # AND SUMMED RAINFALL (USE ONLY CLIMO RAINFALL DATA)
##################
# 1. Pre-allocate a list for ITS and 16S to hold everything
I6S_meteroSummed <- vector("list", length=length(I6S_meteroDatByDay)) 
# Name each element a combination of EU and daysout
names(I6S_meteroSummed) <- names(I6S_meteroDatByDay)
colnames(I6S_meteroDatByDay[[1]])

ITS_meteroSummed <- vector("list", length=length(ITS_meteroDatByDay)) 
# Name each element a combination of EU and daysout
names(ITS_meteroSummed) <- names(ITS_meteroDatByDay)
colnames(ITS_meteroDatByDay[[1]])

# 2. 16S: For loop to get mean and summed values
for (j in 1:length(I6S_meteroSummed)){ #looping over each separate element in the list
  I6S_meteroSummed[[j]] <- I6S_meteroDatByDay[[j]] %>% 
    summarise(
      totalRain_AreaP_cm = sum(PareaInchesRain15min, na.rm = TRUE) * 2.54, #rain is total amount x 2.54 to convert inches to cm
      totalRain_climo_cm = sum(inchesRain15min, na.rm = TRUE) * 2.54,
      mean_climo2mTemp = mean(climo2mTemp, na.rm = TRUE),
      mean_climo2mWSpeed = mean(climo2mWSpeed, na.rm = TRUE),
      max_climo2mTemp = max(climo2mTemp, na.rm = TRUE),
      min_climo2mTemp = min(climo2mTemp, na.rm = TRUE),
      meanRelHumidity_2m = mean(relHumidity_2m, na.rm = TRUE),
      meanTotalRainFall_cm = mean(c(totalRain_AreaP_cm, totalRain_climo_cm)), #mean rainfall is takes from two areas
      meanGust = mean(Gust, na.rm = TRUE), #gust takes from winddata logger
      meanSpeed_dataLogger = mean(Speed, na.rm = TRUE) #speed also takes from winddata logger
      
    )
}
str(I6S_meteroSummed)
I6S_meteroSummed[[1]]

# 3. ITS: For loop to get mean and summed values
for (j in 1:length(ITS_meteroSummed)){
  ITS_meteroSummed[[j]] <- ITS_meteroDatByDay[[j]] %>% 
    summarise(
      totalRain_AreaP_cm = sum(PareaInchesRain15min, na.rm = TRUE) * 2.54,
      totalRain_climo_cm = sum(inchesRain15min, na.rm = TRUE) * 2.54,
      mean_climo2mTemp = mean(climo2mTemp, na.rm = TRUE),
      mean_climo2mWSpeed = mean(climo2mWSpeed, na.rm = TRUE),
      max_climo2mTemp = max(climo2mTemp, na.rm = TRUE),
      min_climo2mTemp = min(climo2mTemp, na.rm = TRUE),
      meanRelHumidity_2m = mean(relHumidity_2m, na.rm = TRUE),
      meanTotalRainFall_cm = mean(c(totalRain_AreaP_cm, totalRain_climo_cm)),
      meanGust = mean(Gust, na.rm = TRUE),
      meanSpeed_dataLogger = mean(Speed, na.rm = TRUE)
    )
}
# Look at some 
ITS_meteroSummed[[1]]
ITS_meteroSummed[[14]]

# 4. Add in dayOut and EU to the objects made above
for (j in 1:length(I6S_meteroSummed)){
  I6S_meteroSummed[[j]]$EU <- unique(I6S_meteroDatByDay[[j]]$EU) #add in EU
  # Add in the day sampled
  I6S_meteroSummed[[j]]$daysOut <- str_extract(names(I6S_meteroSummed)[[j]], "_\\d+$") %>% 
    str_remove(("^_"))
}
names(I6S_meteroSummed)[[1]]
I6S_meteroSummed[[8]]
names(I6S_meteroSummed)[[8]]

# Check that EU and day out are the same in name as in data-- they are!
for (j in 1:length(I6S_meteroSummed)){
  # Gets the name of each of the 16 dfs in I6S_meteroSummed and checks that they match I6S_meteroSummed data within
  print(names(I6S_meteroSummed)[[j]] == paste0(I6S_meteroSummed[[j]]$EU, "_", I6S_meteroSummed[[j]]$daysOut))
}

# ITS
for (j in 1:length(ITS_meteroSummed)){
  ITS_meteroSummed[[j]]$EU <- unique(ITS_meteroDatByDay[[j]]$EU) #add in EU
  # Add in the day sampled
  ITS_meteroSummed[[j]]$daysOut <- str_extract(names(ITS_meteroSummed)[[j]], "_\\d+$") %>% 
    str_remove(("^_"))
}
# View(ITS_meteroSummed[[12]])

# Check that EU and day out are the same in name as in data-- they are!
for (j in 1:length(ITS_meteroSummed)){
  print(names(ITS_meteroSummed)[[j]] == paste0(ITS_meteroSummed[[j]]$EU, "_", ITS_meteroSummed[[j]]$daysOut))
}

ITS_meteroSummed[[11]]

##################
# VI. MERGE RELEVANT SAMPLE METADATA WITH METEORLOGICAL DATA FROM EACH SAMPLING DAY
##################
head(ITS_airMeta)
head(I6S_airMeta)
# 1. Make a column called sampleName to make it match with the qPCR data later on
I6S_airMeta$sampleName <- str_extract(rownames(I6S_airMeta), "_\\d+$") %>% 
  str_remove(("^_"))
cbind(I6S_airMeta$sampleName, rownames(I6S_airMeta)) #looks correct

ITS_airMeta$sampleName <- str_extract(rownames(ITS_airMeta), "_\\d+$") %>% 
  str_remove(("^_"))
cbind(ITS_airMeta$sampleName, rownames(ITS_airMeta)) #looks correct

# 2. Make all of dataframes in the list above into one dataframe with an additional
# column indicating the EU and daysOut combination
# I6S
I6Sall_summedMeterol <- I6S_meteroSummed %>%
  imap_dfr(~ mutate(.x, listName = .y))
# View(I6Sall_summedMeterol) #looks correct
# ITS
ITSall_summedMeterol <- ITS_meteroSummed %>%
  imap_dfr(~ mutate(.x, listName = .y))
# View(ITSall_summedMeterol)

# 3. Merge relevant columns of the metadata with the bigger dataframes made above
# I6S
colnames(I6S_airMeta)
# i. First only keep relevant columns and try to maybe get them in the order I want
I6S_airMeta_trimmed <- I6S_airMeta %>% 
  select("sampleName", "EU","HabitatAir", "DateSetOut","DateRetrieved", "daysOut") %>% 
  tibble::rownames_to_column(var= "fullSampName") %>%  #change rownames to a new column 
  mutate(across(everything(), ~ gsub(" ", "", .))) #remove all spaces
# View(I6S_airMeta_trimmed)

# ITS
colnames(ITS_airMeta)
# First only keep relevant columns and try to maybe get them in the order I want
ITS_airMeta_trimmed <- ITS_airMeta %>% 
  select("sampleName", "EU","HabitatAir", "DateSetOut","DateRetrieved", "daysOut") %>% 
  tibble::rownames_to_column(var= "fullSampName") %>%  #change rownames to a new column 
  mutate(across(everything(), ~ gsub(" ", "", .))) #remove all spaces 
ITS_airMeta_trimmed 
# View(ITS_airMeta_trimmed)

# ii. Merge them!
colnames(I6S_airMeta_trimmed)
colnames(I6Sall_summedMeterol)
head(I6S_airMeta_trimmed)
head(I6Sall_summedMeterol)
unique(I6S_airMeta_trimmed$daysOut)
unique(I6Sall_summedMeterol$daysOut)
# I6S
I6S_metero_Meta_1 <- left_join(I6S_airMeta_trimmed, I6Sall_summedMeterol, by = "daysOut")
# View(I6S_metero_Meta_1)
I6S_metero_Meta_1$EU.x == I6S_metero_Meta_1$EU.y #since the above is true remove EU.y and jsut have EU
I6S_metero_Meta_1$EU.y <- NULL
colnames(I6S_metero_Meta_1)[colnames(I6S_metero_Meta_1) %in% "EU.x"] <- "EU" #rename EU.x EU
# View(I6S_metero_Meta_1)
# ITS
ITS_metero_Meta_1 <- left_join(ITS_airMeta_trimmed, ITSall_summedMeterol, by = "daysOut")
# View(ITS_metero_Meta_1)
ITS_metero_Meta_1$EU.x == ITS_metero_Meta_1$EU.y #since the above is true remove EU.y and jsut have EU
ITS_metero_Meta_1$EU.y <- NULL
colnames(ITS_metero_Meta_1)[colnames(ITS_metero_Meta_1) %in% "EU.x"] <- "EU" #rename EU.x EU
# View(ITS_metero_Meta_1)

##################
# VI. MERGE WITH QPCR DATA TO MAKE THE FINAL TABLES
##################
# 1. First, tidy up the qPCR dataframes a bit
str(I6SqPCRdata_forTable)
I6SqPCRdata_forTable$sampleName <- as.character(I6SqPCRdata_forTable$sampleName)
I6SqPCRdata_forTable$X <- NULL #remove weird X column that was from re-ordering originally
str(ITSqPCRdata_forTable)
ITSqPCRdata_forTable$sampleName <- as.character(ITSqPCRdata_forTable$sampleName)
ITSqPCRdata_forTable$X <- NULL #remove weird X column that was from re-ordering originally

# I6S
I6S_weather_qPCR1 <- left_join(I6SqPCRdata_forTable, I6S_metero_Meta_1, by = "sampleName")
colnames(I6S_weather_qPCR1)
I6S_weather_qPCR_final <- I6S_weather_qPCR1 %>% 
  select("sampleName","EU","HabitatAir", "DateSetOut", "meanSpeed_dataLogger", "meanGust",
         "mean_climo2mTemp", "max_climo2mTemp", "min_climo2mTemp", "meanRelHumidity_2m","totalRain_climo_cm",
         "I6Scopies")
# View(I6S_weather_qPCR_final)
colnames(I6S_weather_qPCR_final)

# ITS
ITS_weather_qPCR1 <- left_join(ITSqPCRdata_forTable, ITS_metero_Meta_1, by = "sampleName")
colnames(ITS_weather_qPCR1)
ITS_weather_qPCR_final <- ITS_weather_qPCR1 %>% 
  select("sampleName","EU","HabitatAir", "DateSetOut", "meanSpeed_dataLogger", "meanGust",
         "mean_climo2mTemp", "max_climo2mTemp", "min_climo2mTemp", "meanRelHumidity_2m","totalRain_climo_cm",
         "ITScopies")
# View(ITS_weather_qPCR_final)
colnames(ITS_weather_qPCR_final)
 
# All together for the figure
colnames(I6S_weather_qPCR_final)
setdiff(ITS_weather_qPCR_final$sampleName, I6S_weather_qPCR_final$sampleName)
setdiff(I6S_weather_qPCR_final$sampleName, ITS_weather_qPCR_final$sampleName) #this shows that there are no 16S samples that are not present 
# in the ITS samples, meaning that the joining below is fine
final_I6SITS_table <- full_join(ITS_weather_qPCR_final, I6S_weather_qPCR_final[,c(1, 12)], by = "sampleName") #just take sampleName and 16SgenomeEquiv
final_I6SITS_table
# View(final_I6SITS_table)

# Export to put in manuscript!!
# write.csv(final_I6SITS_table, file ="~/Desktop/CU_Research/SRS_Aeromicrobiome/rObjectsSaved/MS_figures/supplement_formatted/final_I6SITS_table_Nov21_2025.csv") #made Nov. 21, 2025

##################
# VII. PLOTS
##################
# First, get a summary table for meteorological conditions each day
colnames(final_I6SITS_table) #could have done this from earlier data object, but this is easier
weatherSummedByDay <- final_I6SITS_table %>% 
  group_by(DateSetOut, EU) %>% #group by EU too to keep EU
  summarise(sampDayMeanTemp = mean(mean_climo2mTemp),
            sampDayMinTemp = mean(min_climo2mTemp),
            sampDayMaxTemp = mean(max_climo2mTemp),
            sampDayTotalRainFall_cm = mean(totalRain_climo_cm), #just climo to make it match other climo variables
            sampDayMeanWindSpeed = mean(meanSpeed_dataLogger),
            sampDayMeanGust = mean(meanGust),
            sampDayMeanRelHumidity = mean(meanRelHumidity_2m))
# View(weatherSummedByDay)

# Add back in "daysOut" column for ease of plotting
I6S_metero_Meta_1 # can take from this above given that no days were dropped in the bacterial dataset
which(colnames(I6S_metero_Meta_1) %in% c("DateSetOut", "daysOut"))
# Get unique days to make the merge below go as desired.
daysSimple <- I6S_metero_Meta_1[,which(colnames(I6S_metero_Meta_1) %in% c("DateSetOut", "daysOut"))] %>% distinct()
daysSimple
colnames(weatherSummedByDay)
dim(weatherSummedByDay) #16, 9
weatherSummedByDay <- merge(weatherSummedByDay, daysSimple, by = "DateSetOut")
str(weatherSummedByDay)
# make "daysOut" numeric for plotting as x-axis below!
weatherSummedByDay$daysOut <- as.numeric(weatherSummedByDay$daysOut)
str(weatherSummedByDay)

# 1. TEMPERATURE
# Format data long for plotting 
weatherSummedByDay_long <- pivot_longer(weatherSummedByDay, cols = c("sampDayMeanTemp", "sampDayMinTemp", "sampDayMaxTemp"), 
                                        names_to = "temperatureType", 
                                        values_to = "temperature")
# View(weatherSummedByDay_long)

# Make the plot!
colnames(weatherSummedByDay_long)
str(weatherSummedByDay_long)
tempPlot_withEU <- ggplot(weatherSummedByDay_long, aes(x = daysOut, y = temperature)) +
  geom_line(aes(linetype = temperatureType), color = "black", size = 1) +           # Line for each temperature type
  geom_point(aes(color = EU), size = 3) +          # Points colored by location
  labs(title = "Temperature trends by sampling day",
       x = "Day", y = "Temperature (°C)",
       color = "Location",
       linetype = "Temperature Type") +
  theme_bw()
tempPlot_withEU

# Make a plot that does not have the EUs colored
unique(weatherSummedByDay_long$temperatureType)
unique(as.data.frame(cbind(weatherSummedByDay_long$daysOut, weatherSummedByDay_long$DateSetOut)))
# Re-label daysOut to be DateSetOut on the Plot Itself
dateLabels <-  as.data.frame(matrix(ncol=2, nrow=16))
colnames(dateLabels) <- c("breaks", "labels")
dateLabels[1,] <- c(1,"16 June")
dateLabels[2,] <- c(3,"18 June")
dateLabels[3,] <- c(5,"20 June")
dateLabels[4,] <- c(7,"22 June")
dateLabels[5,] <- c(8,"23 June")
dateLabels[6,] <- c(9,"24 June")
dateLabels[7,] <- c(11,"26 June")
dateLabels[8,] <- c(12,"27 June")
dateLabels[9,] <- c(14,"29 June")
dateLabels[10,] <- c(15,"30 June")
dateLabels[11,] <- c(16,"1 July")
dateLabels[12,] <- c(17,"2 July")
dateLabels[13,] <- c(19,"4 July")
dateLabels[14,] <- c(20,"5 July")
dateLabels[15,] <- c(21,"6 July")
dateLabels[16,] <- c(23,"8 July")
dateLabels$breaks <- as.numeric(dateLabels$breaks)
str(dateLabels)
str(weatherSummedByDay_long)

tempPlot_noEUs <- ggplot(weatherSummedByDay_long, aes(x = daysOut, y = temperature, color = temperatureType)) +
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +  
  # scale_color_manual(values = c("sampDayMeanTemp" = "black", "sampDayMinTemp" = "lightblue", "sampDayMaxTemp" = "red3")) + 
  scale_color_manual(
    values =c("sampDayMeanTemp" = "black", "sampDayMinTemp" = "blue", "sampDayMaxTemp" = "red3"),  #line colors
    labels = c("sampDayMeanTemp" = "mean temp.", 
               "sampDayMaxTemp" = "max. temp.", 
               "sampDayMinTemp" = "min. temp."),
    name = ""
  ) +
  scale_x_continuous(
    breaks = dateLabels$breaks,
    labels = dateLabels$labels) + # Change day labels 
  labs(title = "Temperature trends by sampling day",
       x = NULL, y = "Temperature (°C)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(
    axis.title.y = element_text(size = 10),
    axis.ticks.length = unit(0.2, "cm"),  #increase tick mark length
    axis.text.x = element_text(size = 10),  #increase size of x-axis labels
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10), #increase the size of the legend text 
    plot.title = element_text(size = 12),
    legend.position = "bottom"
  )
tempPlot_noEUs

# 2. WIND
# Format data long for plotting so that can have wind speed and gust on the same plot
colnames(weatherSummedByDay)
windSummedByDay_long <- pivot_longer(weatherSummedByDay, cols = c("sampDayMeanGust", "sampDayMeanWindSpeed"), 
                                     names_to = "windSpeedType", 
                                     values_to = "speed")
# View(windSummedByDay_long)
colnames(windSummedByDay_long)
unique(windSummedByDay_long$windSpeedType)
windPlot_noEUs <- ggplot(windSummedByDay_long, aes(x = daysOut, y = speed, color = windSpeedType)) +
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +  
  scale_color_manual(
    values =c("sampDayMeanGust" = "#af8dc3", "sampDayMeanWindSpeed" = "#7fbf7b"),  #line colors
    labels = c("sampDayMeanGust" = "mean wind\ngust", 
               "sampDayMeanWindSpeed" = "mean wind\nspeed"),
    name = ""
  ) +
  scale_x_continuous(
    breaks = dateLabels$breaks,
    labels = dateLabels$labels) + #change day labels 
  labs(title = "Wind trends by sampling day",
       x = NULL, y = "Wind speed (km per hour)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(
    axis.title.y = element_text(size = 10),
    axis.ticks.length = unit(0.2, "cm"),  #increase tick mark length
    axis.ticks = element_line(size = .7), #increase tick mark thickness
    axis.text.x = element_text(size = 10),  #increase size of x-axis labels
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10), #increase the size of the legend text 
    plot.title = element_text(size = 12),
    legend.position = "bottom"
  )
windPlot_noEUs

# 3. HUMIDITY
colnames(weatherSummedByDay)
humidityPlot_noEUs <- ggplot(weatherSummedByDay, aes(x = daysOut, y = sampDayMeanRelHumidity)) +
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +  
  scale_x_continuous(
    breaks = dateLabels$breaks,
    labels = dateLabels$labels) + #change day labels 
  labs(title = "Relative humidity by sampling day",
       x = NULL, y = "relative humidity (%)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(
    axis.title.y = element_text(size = 10),
    axis.ticks.length = unit(0.2, "cm"),  #increase tick mark length
    axis.ticks = element_line(size = .7), #increase tick mark thickness
    axis.text.x = element_text(size = 10),  #increase size of x-axis labels
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10), #increase the size of the legend text 
    plot.title = element_text(size = 12)
  )
humidityPlot_noEUs
min(weatherSummedByDay$sampDayMeanRelHumidity) #38.02856
max(weatherSummedByDay$sampDayMeanRelHumidity) #90.5766
mean(weatherSummedByDay$sampDayMeanRelHumidity) #70.66712

# 4. RAINFALL
colnames(weatherSummedByDay)
rainPlot_noEUs <- ggplot(weatherSummedByDay, aes(x = daysOut, y = sampDayTotalRainFall_cm)) +
  geom_line(linewidth = 1) + 
  geom_point(size = 2) +  
  scale_x_continuous(
    breaks = dateLabels$breaks,
    labels = dateLabels$labels) + #change day labels 
  labs(title = "Total rainfall by sampling day",
       x = NULL, y = "rainfall (cm)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(
    axis.title.y = element_text(size = 10),
    axis.ticks.length = unit(0.2, "cm"),  #increase tick mark length
    axis.ticks = element_line(size = .7), #increase tick mark thickness
    axis.text.x = element_text(size = 10),  #increase size of x-axis labels
    axis.text.y = element_text(size = 10),
    legend.text = element_text(size = 10), #increase the size of the legend text 
    plot.title = element_text(size = 12)
  )
rainPlot_noEUs

# Plot all of these together
grid.arrange(tempPlot_noEUs, windPlot_noEUs, humidityPlot_noEUs, rainPlot_noEUs, ncol=1)


plot_grid(tempPlot_noEUs, windPlot_noEUs, humidityPlot_noEUs, rainPlot_noEUs, 
          ncol = 1, 
          align = "h",     # vertically align plots
          axis = "l")

library(patchwork)
(tempPlot_noEUs/windPlot_noEUs/humidityPlot_noEUs/rainPlot_noEUs) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "right") #move the legends to the right, will have to re-adjust in InkScape

# Format for journal
# tiff(filename="weatherPlotsSupplemental.tiff",height=5600,width=5200,units="px",res=800,compression="lzw")
# (tempPlot_noEUs/windPlot_noEUs/humidityPlot_noEUs/rainPlot_noEUs) +
#   plot_layout(guides = "collect") &
#   theme(legend.position = "right") #move the legends to the right, will have to re-adjust in InkScape
# dev.off()
