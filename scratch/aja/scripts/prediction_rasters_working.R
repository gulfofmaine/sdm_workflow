fishSDM.prediction.df<- function(rcp85.mu.dir, rcp85.pct05.dir, rcp85.pct95.dir, rcp45.mu.dir, rcp45.pct05.dir, rcp45.pct95.dir, oisst.dir, sp.in, dates.baseline, dates.future, seasonal.mu, season, model.dat) {
  
  library(tidyverse)
  library(maptools)
  library(raster)
  library(rgeos)
  library(geosphere)
  library(zoo)
  
  suppressWarnings(sapply(list.files(pattern = "[.]R$", path = "~/Dropbox/Andrew/Work/GMRI/AllRFunctions/", full.names = TRUE), source))
  
  if(FALSE) {
    rcp85.mu.dir<- paste(res.data.path, "CMIP5_SST/RCP85_mu.grd", sep = "")
    rcp85.pct05.dir<- paste(res.data.path, "CMIP5_SST/RCP85_5th.grd", sep = "")
    rcp85.pct95.dir<- paste(res.data.path, "CMIP5_SST/RCP85_95th.grd", sep = "")
    rcp45.mu.dir<- paste(res.data.path, "CMIP5_SST/RCP45_mu.grd", sep = "")
    rcp45.pct05.dir<- paste(res.data.path, "CMIP5_SST/RCP45_5th.grd", sep = "")
    rcp45.pct95.dir<- paste(res.data.path, "CMIP5_SST/RCP45_95th.grd", sep = "")
    oisst.dir = paste(res.data.path, "OISST/OISSTNWAtl_Projected.grd", sep = "")
    sp.in = "~/Box/RES_Data/Shapefiles/"
    dates.baseline = c("2013-10-16", "2014-10-16", "2015-10-16", "2016-10-16", "2017-10-16", "2018-10-16")
    dates.future =  c("2055-10-16", "2100-10-16")
    dates.baseline = c("2011-10-16", "2012-10-16", "2013-10-16", "2014-10-16", "2015-10-16")
    dates.future =  c("2004-10-16", "2012-10-16", "2025-10-16", "2040-10-16", "2055-10-16", "2100-10-16")
    seasonal.mu = TRUE
    season = "Fall"
    model.dat = here::here("scratch/aja/data", "model_dat.rds")
  }
  
  ## Projections
  proj.wgs84<- CRS("+init=epsg:4326") #WGS84
  proj.utm<- CRS("+init=epsg:2960") #UTM 19
  
  ##### Start
  ## Baseline SSTs
  # Empty stack
  pred.rast.stack<- stack()
  
  # Add OISST
  name.ind<- nlayers(pred.rast.stack)+1
  stack0<- stack(oisst.dir)
  
  # Move to monthly?
  oisst.min<- gsub("X", "", min(names(stack0)))
  oisst.min.date<- as.Date(gsub("[.]", "-", oisst.min))
  oisst.max<- gsub("X", "", max(names(stack0)))
  oisst.max.date<- as.Date(gsub("[.]", "-", oisst.max))
  
  # Calculate monthly mean temperature -- this would be compared to the sstclim data (monthly climate ensemble)
  oisst.dates<- seq.Date(from = oisst.min.date, to = oisst.max.date, by = "day")
  oisst.dat<- setZ(stack0, oisst.dates)
  
  # Aggregate daily to monthly data
  oisst.monthly <- zApply(oisst.dat, by = as.yearmon, mean)
  
  #### Mean seasonal temperature
  if(seasonal.mu) {
    
    # Basline
    #Baseline stack, store seasonal means and then average them all
    sst.stack<- stack()
    
    years<- seq(from = 2013, to = 2018, by = 1)
    
    for(i in seq_along(years)){
      
      dates.use<- switch(season,
                         "Fall" = paste(c("Sep", "Oct", "Nov"), rep(years[i]), sep = "."),
                         "Spring" = paste(c("Mar", "Apr", "May"), rep(years[i]), sep = "."),
                         "Summer" = paste(c("Jun", "Jul", "Aug", "Sep"), rep(years[i]), sep = "."))
      
      sst.temp<- calc(oisst.monthly[[which(names(oisst.monthly) %in% dates.use)]], mean)
      
      sst.stack<- stack(sst.stack, sst.temp)
      print(years[i])
    }
    
    names(sst.stack)<- paste(season, years, sep = ".")
    
    sst.basemeans<- calc(sst.stack[[c(1,2,3,4,5)]], mean)
    
    pred.rast.stack<- stack(pred.rast.stack, sst.basemeans)
    names(pred.rast.stack)[c(1)]<- c("Baseline")
    
    # Climate
    name.ind<- nlayers(pred.rast.stack)+1
    stack0<- projectRaster(stack(rcp85.mu.dir), crs = proj.utm)
    stack0<- resample(stack0, oisst.monthly[[1]])
    clim.dates<- seq.Date(from = as.Date("1982-01-16"), to = as.Date("2300-12-16"), by = "month")
    clim.stack<- setZ(stack0, clim.dates)
    clim.stack.zind<- getZ(clim.stack)
    
    #Baseline stack, store seasonal means and then average them all
    rcp85.mu.stack<- stack()
    
    years<- c("2055", "2100")
    
    for(i in seq_along(years)){
      
      dates.use<- switch(season,
                         "Fall" = as.Date(paste(rep(years[i]), c("09-16", "10-16", "11-16"), sep = "-")),
                         "Spring" = as.Date(paste(rep(years[i]), c("03-16", "04-16", "05-16"), sep = "-")),
                         "Summer" = as.Date(paste(rep(years[i]), c("06-16", "07-16", "08-16", "09-16"), sep = "-")))
      
      clim.mu.temp<- calc(clim.stack[[which(clim.stack.zind %in% dates.use)]], mean)
      
      rcp85.mu.stack<- stack(rcp85.mu.stack, clim.mu.temp)
      print(years[i])
    }
    
    names(rcp85.mu.stack)<- paste(season, years, "rcp85.mu", sep = ".")
    
    # Climate 5th
    name.ind<- nlayers(pred.rast.stack)+1
    stack0<- projectRaster(stack(rcp85.pct05.dir), crs = proj.utm)
    stack0<- resample(stack0, oisst.monthly[[1]])
    clim.dates<- seq.Date(from = as.Date("1982-01-16"), to = as.Date("2300-12-16"), by = "month")
    clim.stack<- setZ(stack0, clim.dates)
    clim.stack.zind<- getZ(clim.stack)
    
    #Baseline stack, store seasonal means and then average them all
    rcp85.pct05.stack<- stack()
    
    for(i in seq_along(years)){
      
      dates.use<- switch(season,
                         "Fall" = as.Date(paste(rep(years[i]), c("09-16", "10-16", "11-16"), sep = "-")),
                         "Spring" = as.Date(paste(rep(years[i]), c("03-16", "04-16", "05-16"), sep = "-")),
                         "Summer" = as.Date(paste(rep(years[i]), c("06-16", "07-16", "08-16", "09-16"), sep = "-")))
      
      clim.pct05.temp<- calc(clim.stack[[which(clim.stack.zind %in% dates.use)]], mean)
      
      rcp85.pct05.stack<- stack(rcp85.pct05.stack, clim.pct05.temp)
      print(years[i])
    }
    
    names(rcp85.pct05.stack)<- paste(season, years, "rcp85.pct05", sep = ".")
    
    # Climate -- 95th
    name.ind<- nlayers(pred.rast.stack)+1
    stack0<- projectRaster(stack(rcp85.pct95.dir), crs = proj.utm)
    stack0<- resample(stack0, oisst.monthly[[1]])
    clim.dates<- seq.Date(from = as.Date("1982-01-16"), to = as.Date("2300-12-16"), by = "month")
    clim.stack<- setZ(stack0, clim.dates)
    clim.stack.zind<- getZ(clim.stack)
    
    #Baseline stack, store seasonal means and then average them all
    rcp85.pct95.stack<- stack()
    
    for(i in seq_along(years)){
      
      dates.use<- switch(season,
                         "Fall" = as.Date(paste(rep(years[i]), c("09-16", "10-16", "11-16"), sep = "-")),
                         "Spring" = as.Date(paste(rep(years[i]), c("03-16", "04-16", "05-16"), sep = "-")),
                         "Summer" = as.Date(paste(rep(years[i]), c("06-16", "07-16", "08-16", "09-16"), sep = "-")))
      
      clim.pct95.temp<- calc(clim.stack[[which(clim.stack.zind %in% dates.use)]], mean)
      
      rcp85.pct95.stack<- stack(rcp85.pct95.stack, clim.pct95.temp)
      print(years[i])
    }
    
    names(rcp85.pct95.stack)<- paste(season, years, "rcp85.pct95", sep = ".")
    
    # RCP45
    # Climate
    name.ind<- nlayers(pred.rast.stack)+1
    stack0<- projectRaster(stack(rcp45.mu.dir), crs = proj.utm)
    stack0<- resample(stack0, oisst.monthly[[1]])
    clim.dates<- seq.Date(from = as.Date("1982-01-16"), to = as.Date("2300-12-16"), by = "month")
    clim.stack<- setZ(stack0, clim.dates)
    clim.stack.zind<- getZ(clim.stack)
    
    #Baseline stack, store seasonal means and then average them all
    rcp45.mu.stack<- stack()
    
    for(i in seq_along(years)){
      
      dates.use<- switch(season,
                         "Fall" = as.Date(paste(rep(years[i]), c("09-16", "10-16", "11-16"), sep = "-")),
                         "Spring" = as.Date(paste(rep(years[i]), c("03-16", "04-16", "05-16"), sep = "-")),
                         "Summer" = as.Date(paste(rep(years[i]), c("06-16", "07-16", "08-16", "09-16"), sep = "-")))
      
      clim.mu.temp<- calc(clim.stack[[which(clim.stack.zind %in% dates.use)]], mean)
      
      rcp45.mu.stack<- stack(rcp45.mu.stack, clim.mu.temp)
      print(years[i])
    }
    
    names(rcp45.mu.stack)<- paste(season, years, "rcp45.mu", sep = ".")
    
    # Climate 5th
    name.ind<- nlayers(pred.rast.stack)+1
    stack0<- projectRaster(stack(rcp45.pct05.dir), crs = proj.utm)
    stack0<- resample(stack0, oisst.monthly[[1]])
    clim.dates<- seq.Date(from = as.Date("1982-01-16"), to = as.Date("2300-12-16"), by = "month")
    clim.stack<- setZ(stack0, clim.dates)
    clim.stack.zind<- getZ(clim.stack)
    
    #Baseline stack, store seasonal means and then average them all
    rcp45.pct05.stack<- stack()
    
    for(i in seq_along(years)){
      
      dates.use<- switch(season,
                         "Fall" = as.Date(paste(rep(years[i]), c("09-16", "10-16", "11-16"), sep = "-")),
                         "Spring" = as.Date(paste(rep(years[i]), c("03-16", "04-16", "05-16"), sep = "-")),
                         "Summer" = as.Date(paste(rep(years[i]), c("06-16", "07-16", "08-16", "09-16"), sep = "-")))
      
      clim.pct05.temp<- calc(clim.stack[[which(clim.stack.zind %in% dates.use)]], mean)
      
      rcp45.pct05.stack<- stack(rcp45.pct05.stack, clim.pct05.temp)
      print(years[i])
    }
    
    names(rcp45.pct05.stack)<- paste(season, years, "rcp45.pct05", sep = ".")
    
    # Climate -- 95th
    name.ind<- nlayers(pred.rast.stack)+1
    stack0<- projectRaster(stack(rcp45.pct95.dir), crs = proj.utm)
    stack0<- resample(stack0, oisst.monthly[[1]])
    clim.dates<- seq.Date(from = as.Date("1982-01-16"), to = as.Date("2300-12-16"), by = "month")
    clim.stack<- setZ(stack0, clim.dates)
    clim.stack.zind<- getZ(clim.stack)
    
    #Baseline stack, store seasonal means and then average them all
    rcp45.pct95.stack<- stack()
    
    for(i in seq_along(years)){
      
      dates.use<- switch(season,
                         "Fall" = as.Date(paste(rep(years[i]), c("09-16", "10-16", "11-16"), sep = "-")),
                         "Spring" = as.Date(paste(rep(years[i]), c("03-16", "04-16", "05-16"), sep = "-")),
                         "Summer" = as.Date(paste(rep(years[i]), c("06-16", "07-16", "08-16", "09-16"), sep = "-")))
      
      clim.pct95.temp<- calc(clim.stack[[which(clim.stack.zind %in% dates.use)]], mean)
      
      rcp45.pct95.stack<- stack(rcp45.pct95.stack, clim.pct95.temp)
      print(years[i])
    }
    
    names(rcp45.pct95.stack)<- paste(season, years, "rcp45.pct95", sep = ".")
    
    # Add it to the pred rast
    pred.rast.stack<- stack(pred.rast.stack, rcp85.mu.stack, rcp85.pct05.stack, rcp85.pct95.stack, rcp45.mu.stack, rcp45.pct05.stack, rcp45.pct95.stack)
    
    # Other predictors
    # Add depth
    neshelf.bathy<- raster(paste(sp.in, "NEShelf_Etopo1_bathy.tiff", sep = ""))
    proj4string(neshelf.bathy)<- proj.wgs84
    depth.temp<- projectRaster(neshelf.bathy, crs = proj.utm)
    DEPTH<- resample(depth.temp, pred.rast.stack[[1]])
    pred.rast.stack<- stack(pred.rast.stack, DEPTH)
    
    # Add TRI
    TRI.temp<- terrain(DEPTH, opt = "TRI")
    TRI<- resample(TRI.temp, pred.rast.stack[[1]])
    pred.rast.stack<- stack(pred.rast.stack, TRI)
    names(pred.rast.stack)[26:27]<- c("DEPTH", "TRI")
    
    # Get these variables
    #Mask out points outside of NELME
    nelme.rast<- pred.rast.stack[[1]]
    nelme.rast[]<- NA
    nelme<- read_sf(paste(sp.in, "NELME_regions/NELME_sf.shp", sep = ""))
    nelme.utm<- st_transform(nelme, proj.utm)
    nelme.buff<- gBuffer(as_Spatial(nelme.utm), width = 40000)
    nelme.rast<- rasterize(nelme.buff, nelme.rast)
    pred.rast.stack.m<- raster::mask(pred.rast.stack, mask = nelme.rast, inverse = FALSE)
    
    # Now values
    pred.rast.stack.m.wgs84<- projectRaster(pred.rast.stack.m, crs = proj.wgs84)
    pred.df<- raster::as.data.frame(pred.rast.stack.m.wgs84, xy = T)
    points.wgs84<- pred.df
    coordinates(points.wgs84)<- ~x+y
    proj4string(points.wgs84)<- proj.wgs84
    
    # Getting sediment type information
    neshelf.sediment<- as_Spatial(read_sf(paste(sp.in, "Sediment/TNC_benthicsediment.shp", sep = "")))
    proj4string(neshelf.sediment)<- proj.wgs84
    
    # Use over to get polygon data at each point
    pred.df$SED.SIZE<- over(points.wgs84, neshelf.sediment)$SEDIMENT
    pred.df$SED.TYPE<- over(points.wgs84, neshelf.sediment)$GRPSED
    
    # Add in new along/cross shelf position
    pred.df$SHELF_POS<- distCosine(points.wgs84, cbind(-75, 35), r=6378137)/1000
    
    return(pred.df)
    
  } else {
    for(i in 1:length(oisst.windows)) {
      window<- round(oisst.windows[i]/30, 0)
      dates.stack.temp<- stack()
      
      for(k in 1:length(dates.baseline)) {
        date.new<- unlist(strsplit(gsub("-", ".", format(as.Date(dates.baseline[k]), "%Y-%b")), "[.]"))
        date2<- paste(date.new[2], date.new[1], sep = ".")
        
        stack.start<- which(names(oisst.monthly) == date2)
        stack.t<- oisst.monthly[[(stack.start-(window-1)):stack.start]]
        datemu.t<- calc(stack.t, mean)
        dates.stack.temp<- stack(dates.stack.temp, datemu.t)
      }
      mu.t<- calc(dates.stack.temp, mean) 
      pred.rast.stack<- stack(pred.rast.stack, mu.t)
      names(pred.rast.stack)[name.ind]<- paste("d", oisst.windows[i], "MU.OISST", sep = "")
      name.ind<- name.ind+1
      print(paste(window, " is done", sep = ""))
    }
    
    ## Climate projections SSTs
    name.ind<- nlayers(pred.rast.stack)+1
    stack0<- stack(climate.dir)
    pred.rast.stack<- raster::resample(pred.rast.stack, stack0[[1]])
    
    for(i in 1:length(dates.future)) {
      year<- format(as.Date(dates.future[[i]]), "%Y")
      stack.start<- which(gsub("[.]", "-", gsub("X", "", names(stack0))) == dates.future[i])
      
      for(k in 1:length(climate.windows)) {
        window<- round(climate.windows[k]/30, 0)
        stack.t<- stack0[[(stack.start-(window-1)):stack.start]]
        datemu.t<- calc(stack.t, mean)
        pred.rast.stack<- stack(pred.rast.stack, datemu.t)
        names(pred.rast.stack)[name.ind]<- paste(year, ".", window, "MO.Clim", sep = "")
        name.ind<- name.ind+1
        print(paste(dates.future[i], climate.windows[k], " is done", sep = " "))        
      }
    }
    
    # Add depth
    neshelf.bathy<- raster(paste(sp.in, "NEShelf_etopo1_bathy_reclass.tif", sep = ""))
    proj4string(neshelf.bathy)<- proj.wgs84
    depth.temp<- projectRaster(neshelf.bathy, crs = proj.utm)
    DEPTH<- resample(depth.temp, pred.rast.stack[[1]])
    pred.rast.stack<- stack(pred.rast.stack, DEPTH)
    
    # Add TRI
    TRI.temp<- terrain(DEPTH, opt = "TRI")
    TRI<- resample(TRI.temp, pred.rast.stack[[1]])
    pred.rast.stack<- stack(pred.rast.stack, TRI)
    
    # Get these variables
    #Mask out points outside of NELME
    nelme.rast<- pred.rast.stack[[1]]
    nelme.rast[]<- NA
    nelme<- readShapePoly(paste("~/Dropbox/Andrew/Work/GMRI/AllGIS/nelme.shp", sep = ""))
    proj4string(nelme)<- proj.wgs84
    nelme.utm<- spTransform(nelme, proj.utm)
    nelme.buff<- gBuffer(nelme.utm, width = 40000)
    nelme.rast<- rasterize(nelme.buff, nelme.rast)
    pred.rast.stack.m<- mask(pred.rast.stack, mask = nelme.rast, inverse = FALSE)
    
    # Species specific biomass
    # Load it
    temp.space <- new.env()
    temp.df <- load(model.dat, temp.space)
    dat.all <- get(temp.df, temp.space)
    rm(temp.space)
    
    # Add in new along/cross shelf position
    proj.wgs84<- CRS("+init=epsg:4326") #WGS84
    proj.utm<- CRS("+init=epsg:2960") #UTM 19
    
    pts<- data.frame("x" = dat.all$DECDEG_BEGLON, "y" = dat.all$DECDEG_BEGLAT)
    coordinates(pts)<- ~x+y
    proj4string(pts)<- proj.utm
    pts.sp<- data.frame(spTransform(pts, proj.wgs84))
    
    dat.all$SHELF_POS<- distCosine(pts.sp, cbind(-75, 35), r=6378137)/1000
    
    # Filter
    species.all<- read.csv("~/Dropbox/Andrew/Work/GMRI/AllData/Assesmentfishspecies.csv")
    dat<- filter(dat.all, COMNAME %in% species.all$COMNAME)
    
    # Formatting datasets for this specific modeling strucutre
    dat$BIOMASS.ADJ<- ifelse(dat$bio.abund.flag == "FLAG1", 1, dat$BIOMASS) # Abundance, but no biomass recorded. Set these records to biomass = 1.
    dat$BIOMASS.ADJ<- ifelse(is.na(dat$bio.abund.flag), 0, dat$BIOMASS.ADJ) # True absences, change NA to 0.
    
    # Create BIOMASS.MOD for modeling log biomass
    dat<- dat %>%
      mutate(.,
             "BIOMASS.LOG" = log(BIOMASS.ADJ),
             "BIOMASS.MOD" = ifelse(BIOMASS.LOG == -Inf | BIOMASS.LOG <=0, NA, BIOMASS.LOG),
             "RANDOM" = rnorm(nrow(.)))
    
    # Create Mean Biomass
    load("~/Dropbox/Andrew/Work/GMRI/AllData/stratum.area.Rdata")
    dat<- left_join(dat, stratum.area, by = "STRATUM")
    
    # Weighted mean biomass
    # Get STRATUM code and STRATUM_AREA and only unique observations for each combo
    t1<- dplyr::select(dat, STRATUM, STRATUM_AREA)
    t2<- t1[!duplicated(t1["STRATUM"]),]
    totstrwt<-sum(t2$STRATUM_AREA, na.rm = TRUE)
    
    # Get STRATUM proportion area relative to total area surveyed
    ratio.dat<- data.frame("STRATUM" = t2$STRATUM, "STRATUM.RATIO" = t2$STRATUM_AREA/totstrwt)
    
    # Get a count of the unique number of tows per stratum per year  
    tows_unique<- dat[!duplicated(dat["ID"]),] 
    numtows<-aggregate(ID ~ EST_YEAR + STRATUM, length, data= tows_unique)
    colnames(numtows)<-c('EST_YEAR','STRATUM','COUNT')
    
    dat<- left_join(dat, numtows, by = c("EST_YEAR", "STRATUM"))
    dat<- left_join(dat, ratio.dat, by = "STRATUM")
    
    # Calculate total species biomass per strata per year
    strat.biomass<-aggregate(BIOMASS.ADJ ~ EST_YEAR + SVSPP + STRATUM + STRATUM.RATIO, sum, data = dat)
    strat.biomass<- left_join(strat.biomass, numtows, by = c("EST_YEAR", "STRATUM"))
    
    # Mean biomass per year/species/strata 
    strat.biomass$mean.biomass<- strat.biomass$BIOMASS.ADJ/strat.biomass$COUNT
    
    # Area weighted mean biomass per year/species/strata 
    strat.biomass$BIOMASS.WMEAN<-strat.biomass$mean.biomass*strat.biomass$STRATUM.RATIO
    
    # Mean Annual Biomass
    annual.seasonal.wtmean<- strat.biomass %>%
      group_by(., EST_YEAR, SEASON, SVSPP) %>%
      summarise(., "BIOMASS.WMEAN" = sum(BIOMASS.WMEAN)) 
    
    avg.seasonal.wtmean<- annual.seasonal.wtmean %>%
      group_by(., SEASON, SVSPP) %>%
      summarize(., "AVG.BIOMASS.WMEAN" = mean(BIOMASS.WMEAN, na.rm = T))
    
    annual.seasonal.wtmean<- annual.seasonal.wtmean %>%
      left_join(., avg.seasonal.wtmean, by = c("SEASON", "SVSPP")) %>%
      mutate(., "BIOMASS.WMEAN.ANOM" = BIOMASS.WMEAN - AVG.BIOMASS.WMEAN)
    
    # Add back to full data
    dat<- left_join(dat, annual.wtmean, by = c("EST_YEAR", "SVSPP"))
    
    # Get baseline observations
    dat$TRAIN.TEST<- ifelse(as.Date(dat$DATE) >= "1986-01-01" & as.Date(dat$DATE) <= "2010-12-31", "TRAIN", 
                            ifelse(as.Date(dat$DATE) >= "2011-01-01" & as.Date(dat$DATE) <= "2016-12-31", "TEST", "Neither")) 
    
    dat.test<- dat %>%
      group_by(., COMNAME, TRAIN.TEST) %>%
      dplyr::filter(., TRAIN.TEST == "TEST") %>%
      nest(.key = "TEST.DATA") %>%
      arrange(COMNAME) 
    
    # Rasters ready
    stack.id<- nlayers(pred.rast.stack.m)
    count.rast.temp<- pred.rast.stack.m[[1]]
    count.rast.temp[]<- NA
    names(count.rast.temp)<- "Count"
    
    for(i in seq_along(as.character(unique(dat.test$COMNAME))))  {
      temp<- dat.test$TEST.DATA[[i]]
      temp.sp<- temp
      coordinates(temp.sp)<- ~DECDEG_BEGLON+DECDEG_BEGLAT
      proj4string(temp.sp)<- proj.utm
      
      # Get average weighted mean biomass in each grid cell
      temp$cell.id<- extract(count.rast.temp, temp.sp, cellnumbers=TRUE)[,1]
      temp2<- temp %>%
        group_by(cell.id) %>%
        dplyr::summarise(., mean = mean(BIOMASS.ADJ))
      
      # Save values in raster
      out.rast<- count.rast.temp
      out.rast[temp2$cell.id]<- temp2$mean
      
      # Add em to the stack
      pred.rast.stack.m<- stack(pred.rast.stack.m, out.rast)
      names(pred.rast.stack.m)[nlayers(pred.rast.stack.m)]<- paste(dat.test$COMNAME[i], ".", "Biomass", sep = "")
    }
    
    # Now values
    points<- data.frame(coordinates(pred.rast.stack.m[[1]]))
    coordinates(points)<- ~x+y
    proj4string(points)<- proj.utm
    
    pred.df<- data.frame(raster::extract(pred.rast.stack.m, points))
    names(pred.df)[21]<- "DEPTH"
    names(pred.df)[22]<- "TRI"
    
    # Add sediment
    neshelf.sediment<- readShapePoly(paste(sp.in, "TNC_benthicsediment.shp", sep = ""))
    proj4string(neshelf.sediment)<- proj.wgs84
    neshelf.sediment.utm<- spTransform(neshelf.sediment, proj.utm)
    neshelf.sediment.pts<- over(points, neshelf.sediment.utm)$GRPSED
    pred.df$SED.TYPE<- as.factor(neshelf.sediment.pts)
    
    points.wgs84<- spTransform(points, proj.wgs84)
    pred.df$x<- coordinates(points.wgs84)[,1]
    pred.df$y<- coordinates(points.wgs84)[,2]
    
    return(pred.df)
  }
}

res.data.path<- "~/Box/RES_Data/"
fall.rast.pred<- fishSDM.prediction.df(rcp45.mu.dir = paste(res.data.path, "CMIP5_SST/RCP45_mu.grd", sep = ""),  rcp45.pct05.dir = paste(res.data.path, "CMIP5_SST/RCP45_5th.grd", sep = ""), rcp45.pct95.dir = paste(res.data.path, "CMIP5_SST/RCP45_95th.grd", sep = ""), rcp85.mu.dir = paste(res.data.path, "CMIP5_SST/RCP85_mu.grd", sep = ""), rcp85.pct05.dir = paste(res.data.path, "CMIP5_SST/RCP85_5th.grd", sep = ""), rcp85.pct95.dir = paste(res.data.path, "CMIP5_SST/RCP85_95th.grd", sep = ""), oisst.dir = paste(res.data.path, "OISST/OISSTNWAtl_Projected.grd", sep = ""), sp.in = "~/Box/RES_Data/Shapefiles/", dates.baseline = c("2013-10-16", "2014-10-16", "2015-10-16", "2016-10-16", "2017-10-16", "2018-10-16"), dates.future =  c("2055-10-16", "2100-10-16"), seasonal.mu = TRUE, season = "Fall", model.dat = here::here("scratch/aja/data", "model_dat.rds"))
fall.rast.pred$SEASON<- "FALL"
saveRDS(fall.rast.pred, file = here("scratch/aja/Data", "fall.rast.preds.01122021.rds"))

spring.rast.pred<- fishSDM.prediction.df(rcp45.mu.dir = paste(res.data.path, "CMIP5_SST/RCP45_mu.grd", sep = ""),  rcp45.pct05.dir = paste(res.data.path, "CMIP5_SST/RCP45_5th.grd", sep = ""), rcp45.pct95.dir = paste(res.data.path, "CMIP5_SST/RCP45_95th.grd", sep = ""), rcp85.mu.dir = paste(res.data.path, "CMIP5_SST/RCP85_mu.grd", sep = ""),rcp85.pct05.dir = paste(res.data.path, "CMIP5_SST/RCP85_5th.grd", sep = ""), rcp85.pct95.dir = paste(res.data.path, "CMIP5_SST/RCP85_95th.grd", sep = ""), oisst.dir = "~/Box/Mills Lab/Projects/COCA15_ClimVuln/COCA-SDM/Data/1981to2015_OISST_Projected.grd", sp.in = "~/Box/Data/Shapefiles/", dates.baseline = c("2011-10-16", "2012-10-16", "2013-10-16", "2014-10-16", "2015-10-16"), dates.future =  c("2004-10-16", "2012-10-16", "2025-10-16", "2040-10-16", "2055-10-16", "2100-10-16"), seasonal.mu = TRUE, season = "Spring", model.dat = "~/GitHub/COCA/Data/dat.fitted.biomass.rds")
spring.rast.pred$SEASON<- "SPRING"
saveRDS(spring.rast.pred, file = here("scratch/aja/Data", "spring.rast.preds.01122021.rds"))







## RCP4.5
nelme<- st_read("~/GitHub/SDMSeasonality/Data/NELME_clipped.shp")
proj.wgs84<- "+init=epsg:4326"
st_crs(nelme)<- proj.wgs84

clim.proj<- readRDS("~/Box/Data/CMIP5_SST/ProcessedSSTProjections.rds")

# Mean, rcp45 and rcp85
clim.summs<- clim.proj %>%
  separate(Year, into = c("Year", "Month", "Day")) %>%
  group_by(Scenario, Year, Month, x, y) %>%
  summarize("Mean" = mean(SST, na.rm = TRUE),
            "Pct5th" = quantile(SST, probs = c(0.05), na.rm = TRUE, names = FALSE),
            "Pct95th" = quantile(SST, probs = c(0.95), na.rm = TRUE, names = FALSE))

clim.summs<- clim.summs %>%
  group_by(Scenario, Year, Month) %>%
  nest()

df_to_rast<- function(df, stat, mask = nelme) {
  if(FALSE){
    df<- clim.summs$data[[1]]
    stat<- "Mean.SST"
  }
  
  df.temp<- df %>%
    dplyr::select(x, y, stat)
  
  rast.temp<- rasterFromXYZ(df.temp)
  rast.out<- mask(rast.temp, nelme)
  return(rast.out)
}

clim.summs<- clim.summs %>%
  mutate(., "RasterStack.Mean" = map2(data, "Mean", df_to_rast),
         "RasterStack.Pct05" = map2(data, "Pct5th", df_to_rast),
         "RasterStack.Pct95" = map2(data, "Pct95th", df_to_rast))
