library(tidyverse)
library(lubridate)

first_date = as.Date("2020-03-30")
last_date = as.Date("2021-10-31")

oc_counts_data <- httr::GET("https://services2.arcgis.com/LORzk2hk9xzHouw9/arcgis/rest/services/oc_covid_count/FeatureServer/0/query?where=1%3D1&outFields=*&outSR=4326&f=json",
                            query = list(
                              outFields = "*",
                              where = "1=1"
                            )
) %>%
  httr::content(as = "text", encoding = "UTF-8") %>%
  jsonlite::fromJSON(flatten = T) %>%
  magrittr::use_series(features) %>%
  as_tibble() %>%
  transmute(
    dates = as_date(as_datetime(attributes.date/1000)),
    I = attributes.dailycases_specimen,
    Region = "OC"
  ) %>% filter(dates >= first_date)  %>% filter(dates <= last_date) %>% arrange(dates)
