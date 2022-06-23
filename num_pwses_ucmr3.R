# Testing for number of PWSs in UCMR3 
# https://echo.epa.gov/tools/data-downloads
# https://echo.epa.gov/tools/data-downloads/sdwa-download-summary#sdwis
# https://www.govinfo.gov/content/pkg/FR-2012-05-02/pdf/2012-9978.pdf

library(tidyverse)

setwd('C:/Users/Aakas/Desktop/Stuff for DWJ/PFAS Project')

sdwa_intrest <- read_csv('Data/CSVs/SDWA_present/SDWA_PUB_WATER_SYSTEMS.csv')

ucmr_3 <- read_csv('Data/CSVs/UCMR_3_data.csv')

pwsids <- ucmr_3$PWSID %>% unique()

names <- ucmr_3$PWSName %>% unique()

facid <- ucmr_3$FacilityID %>% unique()

facnames <- ucmr_3$FacilityName %>% unique()


num_cws_plus_ntncws <- function(sdwa_dataframe) {
  sdwa_dataframe <- sdwa_dataframe[is.na(sdwa_dataframe$PWS_DEACTIVATION_DATE), ]
  tabled <- sdwa_dataframe$PWS_TYPE_CODE %>% table()
  if (length(tabled) == 3) {
    return((tabled[1] + tabled[2]) %>% as.numeric())
  }
  if (length(tabled) == 4) {
    return((tabled[1] + tabled[3]) %>% as.numeric())
  }
}


cws_ntncws_2022 <- num_cws_plus_ntncws(sdwa_intrest)

# We can get all CWSs active in 2014 by working backwards and filtering for years
# where the CWS was first reported

sdwa_intrest_2015 <- sdwa_intrest %>%
  mutate(year_created = FIRST_REPORTED_DATE %>%
           substr(7, 13) %>%
           as.integer()) %>%
  filter(year_created <= 2015)

sdwa_intrest_2014 <- sdwa_intrest_2015 %>% filter(year_created <= 2014)
sdwa_intrest_2013 <- sdwa_intrest_2015 %>% filter(year_created <= 2013)

cws_ntncws_2015 <- num_cws_plus_ntncws(sdwa_intrest_2015) 
cws_ntncws_2014 <- num_cws_plus_ntncws(sdwa_intrest_2014)
cws_ntncws_2013 <- num_cws_plus_ntncws(sdwa_intrest_2013) 

percent_tested <- data.frame(year = c(2013, 2014, 2015, 2022),
                             ntncws_and_cws = c(cws_ntncws_2013, cws_ntncws_2014,
                                                cws_ntncws_2015, cws_ntncws_2022)) %>%
  mutate(percent_tested = 100 * (facnames %>% length()) / ntncws_and_cws)
