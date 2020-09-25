rm(list = ls()); 
library(tidyverse)
library(dplyr)
library(stringr)
#----------------------------------------------------------------------------
# Texas Mobility analysis
#----------------------------------------------------------------------------
# File from: https://www.google.com/covid19/mobility/
Global_Mobility_Report <- read.csv("Global_Mobility_Report.csv")

#Tidy:
GMR_TX = Global_Mobility_Report %>% filter(country_region == "United States",
                                           sub_region_1 == "Texas",
                                           sub_region_2 != "") %>%
  select(-c("country_region_code", "country_region","sub_region_1","iso_3166_2_code","census_fips_code")) %>%
  rename(county = sub_region_2) %>%
  mutate(county =  strsplit(county, " ")[[1]][1], date = as.Date(date, "%Y-%m-%d")) 
write.csv(GMR_TX, file = "GMR_TX.csv", fileEncoding = "macroman")

d = dat_cf.first10 %>% left_join(GMR_TX, by = c("date","county" ))
View(d)
