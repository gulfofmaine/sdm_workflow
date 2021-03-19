#download data from ftp://ftp.dfo-mpo.gc.ca/MarPED/
#extract it somewhere on your computer

#install the Mar.datawrangling package from github
library(devtools)
install_github("Maritimes/Mar.datawrangling")
library(Mar.datawrangling)

#set data.dir to download extraction location
data.dir = "C:/<wherever you put the data>/RVSurvey_YYYYMMDD"

#Please see wiki at https://github.com/Maritimes/Mar.datawrangling/wiki for how to 
#interrogate the data

#Here's a sample extraction of Summer 2017 thorny skate data (sp code = 201)

#step 1: load data from the files you downloaded into R
  get_data('rv', data.dir = data.dir)

#step 2: filter the data as desired (overwriting the loaded objects)
  get_survey('rv', survey = "SUMMER")
  GSSPECIES = GSSPECIES[GSSPECIES$CODE == 201, ]
  GSMISSIONS = GSMISSIONS[GSMISSIONS$YEAR == 2017,]

#step 3: run self_filter to drop all of the records that are not relevant to 
#        the filters you specified above      
  self_filter('rv',keep_nullsets = T)

#step 4 (optional): save the data locally as a csv and a shapefile.
  save_data('rv',filename = "thornyskates_summer2017_4X",formats = c('csv','shp'))

#messed up a filter or need to look for something else? re-load the data, and start again...
  get_data('rv', data.dir = data.dir)  