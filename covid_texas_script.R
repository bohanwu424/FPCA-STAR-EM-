rm(list = ls()); 
#----------------------------------------------------------------------------
# Texas COVID analysis
#----------------------------------------------------------------------------
library(tidyverse)
# File from: https://dshs.texas.gov/coronavirus/AdditionalData.aspx
#----------------------------------------------------------------------------
# Hospitalization data
#----------------------------------------------------------------------------
dat_h = read_csv("Texas COVID-19 Hospitalizations by TSA.csv")

# Tidy:
dat_h = dat_h %>% pivot_longer(
  cols = colnames(dat_h)[-(1:2)], 
  names_to = "date", 
  values_to = "hospitalizations") %>%
  mutate(date = as.Date(paste(date, "/2020", sep=''), "%m/%d/%Y"))

dat_h_city = dat_h %>% 
  filter(`TSA Name` == "Statewide Total")
  #filter(`TSA Name` == "Houston")

plot(dat_h_city$date, dat_h_city$hospitalizations)
#----------------------------------------------------------------------------
# Fatalities data
#----------------------------------------------------------------------------
dat_f = read_csv(file = 'Texas COVID-19 Fatality Count Data by County.csv') 

# Tidy:
dat_f = dat_f %>% pivot_longer(
  cols = colnames(dat_f)[-1], 
  names_to = "date", 
  values_to = "total fatalities") %>%
  mutate(date = as.Date(paste(date, "-2020", sep=''), "%m-%d-%Y"))

# Add the new fatalities by county:
dat_f = dat_f %>% 
  group_by(`County Name`) %>% 
  mutate(fatalities = c(NA, diff(`total fatalities`)))

dat_f_city = dat_f %>% 
  filter(`County Name` == "Harris")


plot(dat_f_city$date, dat_f_city$`total fatalities`)
plot(dat_f_city$date, dat_f_city$fatalities)
#----------------------------------------------------------------------------
# Case count data
#----------------------------------------------------------------------------
dat_c = read_csv(file = 'Texas COVID-19 Case Count Data by County.csv') 
d0 = unlist(lapply(gregexpr(pattern = '0', colnames(dat_c)[-(1:2)]),
       function(z) z[1]))
colnames(dat_c)[-(1:2)] = substr(colnames(dat_c)[-(1:2)], d0, 500)

# Tidy:
dat_c = dat_c %>% pivot_longer(
  cols = colnames(dat_c)[-(1:2)], 
  names_to = "date", 
  values_to = "total cases") %>%
  mutate(date = as.Date(paste(date, "-2020", sep=''), "%m-%d-%Y"))

# Add the new fatalities by county:
dat_c = dat_c %>% 
  group_by(`County Name`) %>% 
  mutate(cases = c(NA, diff(`total cases`)))

dat_c_city = dat_c %>% 
  filter(`County Name` == "Harris")

plot(dat_c_city$date, dat_c_city$`total cases`)
plot(dat_c_city$date, dat_c_city$cases)
#----------------------------------------------------------------------------
# Merge cases and fatalities:
#----------------------------------------------------------------------------
dat_cf = inner_join(dat_c, dat_f)

dat_cf_city = dat_cf %>% 
  filter(`County Name` == "Harris")
plot(dat_cf_city$date, log(dat_cf_city$cases + 1))
lines(dat_cf_city$date, log(dat_cf_city$fatalities + 1), type = 'p', pch =6)
lines(dat_h_city$date, log(dat_h_city$hospitalizations + 1), type='p', pch = 5)