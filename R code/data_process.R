######## Last updated 09/26/23
######## by Mikyoung Jun

######## This file contains codes to process terrorism data from GTD, and covariates used from other sources in order to prepare them in a format 
########  to be used for the analysis


####### Terrorism data from GTD ########

### the following shows how one might process data for each country (just an example) ###

library(readxl)

gtd1 = read_excel("globalterrorismdb_0522dist.xlsx")
gtd2 = read_excel("globalterrorismdb_2021Jan-June_1222dist.xlsx")

cols = c("eventid", "iyear", "imonth", "iday", "country_txt", "longitude", "latitude", "specificity", "gname")
gtd1 = gtd1[, cols]
gtd2 = gtd2[, cols]

gtd = rbind(gtd1, gtd2)

gtd$time = paste(gtd$iyear, gtd$imonth, gtd$iday, sep="-")

gtd$date = as.Date(gtd$time, "%Y-%m-%d")

head(gtd)

gtd$time = NULL

colnames(gtd)[2:5] = c("year", "month", "day", "country")
head(gtd)

### Afghanistan (year 2013 as an example) ###

afg13 = gtd %>%
  filter(country == "Afghanistan", year==2013) %>%
  filter(gname=="Taliban") 

afg13 = afg13 %>%
  filter_at(vars(longitude, latitude), all_vars(!is.na(.)))


save(afg13, file="xxx.RData") ## xxx for file name


### Nigeria ###

## do similarly to above but extract for Boko Haram and Fulani Extremists 





####### population data from WorldPop #######

library(raster)
pop2013_wp = raster("ppp_2013_1km_Aggregated.tif")

xy = afg13[, c("longitude", "latitude")]

afg13$pop = extract(x=pop2013_wp, y=xy, method="bilinear")



####### elevation data ########

## elevation data was not included in the final analysis of the paper although it was explored in the preliminary analysis. 
## look up the following to get more info on high resolution elevation data used for the preliminary analysis:
## https://cran.r-project.org/web/packages/elevatr/vignettes/introduction_to_elevatr.html


