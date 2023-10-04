##### There are multiple data sets used in the work for terrorism and accompanying covariate info. #####

### Terrorism data ###

1) Source: GTD: https://www.start.umd.edu/gtd/contact/download

START (National Consortium for the Study of Terrorism and Responses to Terrorism). (2022). Global Terrorism Database 1970 - 2020 [data file]. https://www.start.umd.edu/gtd

2) The analysis in the paper considers two countries, Afghanistan and Nigeria. For Afghanistan, attacks by Taliban from Jan 1, 2002 to Dec 31, 2013 were considered. For Nigeria, attachks by 
Boko Haram and Fulani Extremists from Jan 1, 2010 to Dec 31, 2017 were considered.
3) The following variables for each attack record were used: 
- eventid: a 12-digit Event ID
- year, month, day: year/month/day of the month, in/on which the incident occurred
- country: country where the incident occurred
- longitude, latitude: longitude and latitude (based on WGS1984 standards) of the city in which the event occurred
- specificity: the geospatial resolution of the latitude and longitude fields.
    1 = event occurred in city/village/town and lat/long is for that location
    2 = event occurred in city/village/town and no lat/long could be found, so coordinates are for centroid of smallest subnational administrative region identified
    3 = event did not occur in city/village/town, so coordinates are for centroid of smallest subnational administrative region identified
    4 = no 2nd order or smaller region could be identified, so coordinates are for center of 1st order administrative region
    5 = no 1st order administrative region could be identified for the location of the attack, so latitude and longitude are unknown
- gname: the name of the group that carried out the attack
- date: year-month-day


### Covariate information ###
1) Population: from the World Pop project, https://hub.worldpop.org/geodata/summary?id=24771 (resolution of approximately 1km at the equator)
2) Elevation: from an R package, elevatr. 

Check out data_process.R in this directory for examples on how we processed the data.
