source("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/R/vast_functions.R")
library(VAST)
library(sf)
library(tidyverse)

data(georges_bank_haddock_fall, package = "FishStatsUtils")
dat_use<- georges_bank_haddock_fall
strata_use_orig<- load_example(data_set = "GB_fall_haddock")$strata.limits

nmfs<- st_read("~/Box/RES_Data/shapefiles/nmfs_trawl_regions/nmfs_trawl_regions_collection.geojson")
strata_data<- data.frame("Region" = "Georges_Bank")

nelme<- st_read("~/Box/RES_Data/Shapefiles/large_marine_ecosystems/northeast_us_continental_shelf_exterior.geojson")

# Region
region_keep<- nmfs %>%
  filter(., finstr_id %in% strata_use_orig[[1]][1:15]) %>%
  st_transform(., st_crs(nelme)) %>%
  st_union() %>%
  st_as_sf(., strata_data)

# Shapefiles for indices
shape_keep_1<- nmfs %>%
  filter(., finstr_id %in% strata_use_orig[[1]][1:5]) %>%
  st_transform(., st_crs(nelme)) %>%
  st_union() %>%
  st_as_sf(., strata_data)
shape_keep_1$Region[1]<- "Georges_Bank_A"
shape_keep_2<- nmfs %>%
  filter(., finstr_id %in% strata_use_orig[[1]][6:15]) %>%
  st_transform(., st_crs(nelme)) %>%
  st_union() %>%
  st_as_sf(., strata_data)
shape_keep_2$Region[1]<- "Georges_Bank_B"
shape_keep_3<- nmfs %>%
  filter(., finstr_id %in% strata_use_orig[[1]][1:15]) %>%
  st_transform(., st_crs(nelme)) %>%
  st_union() %>%
  st_as_sf(., strata_data)
shape_keep_3$Region[1]<- "All"

# shape_keep<- bind_rows(shape_keep_1, shape_keep_2, shape_keep_3)
shape_keep<- bind_rows(shape_keep_1, shape_keep_2, shape_keep_3)

# I *think* we will need a new column that basically catches the data that are inside this polygon. Check this -- I don't think we need the strata column in the data...
dat_use_sf<- dat_use %>%
  st_as_sf(., coords = c("LONGITUDE", "LATITUDE"), crs = st_crs(nelme), remove = FALSE)

dat_use2<- dat_use_sf %>%
  st_drop_geometry()

strata_use<- data.frame("STRATA" = c("Georges_Bank_A", "Georges_Bank_B", "All"))

# Transform crs of shapefile to common WGS84 lon/lat format.
cell_size = 2500
region_shapefile<- region_keep
region_wgs84<- st_transform(region_shapefile, crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ")

# Get UTM zone
lon<- sum(st_bbox(region_wgs84)[c(1,3)])/2
utm_zone<- floor((lon + 180)/6)+1

# Transform to the UTM zone
crs_utm<- st_crs(paste0("+proj=utm +zone=",utm_zone," +ellps=WGS84 +datum=WGS84 +units=m +no_defs "))
region_utm<- st_transform(region_wgs84, crs = crs_utm)

# Make extrapolation grid with sf
region_grid<- st_as_sf(st_make_grid(region_utm, cellsize = cell_size, what = "centers"), crs = crs_utm) 

# Now get only the points that fall within the shape polygon
points_keep<- data.frame("pt_row" = seq(from = 1, to = nrow(region_grid), by = 1), "in_out" = st_intersects(region_grid, region_utm, sparse = FALSE))            
region_grid<- region_grid %>%
  mutate(., "in_poly" = st_intersects(region_grid, region_utm, sparse = FALSE)) %>%
  filter(., in_poly == TRUE)

# Convert back to WGS84 lon/lat, as that is what VAST expects and add in flag for the shapefile indices...
extrap_grid<- region_grid %>%
  st_transform(., crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ") %>%
  st_join(., shape_keep, join = st_within) %>%
  mutate(., "Lon" = as.numeric(st_coordinates(.)[,1]),
         "Lat" = as.numeric(st_coordinates(.)[,2])) %>%
  st_drop_geometry() %>%
  dplyr::select(., Lon, Lat, Region) %>%
  mutate(., Area_km2=((cell_size/1000)^2),
         STRATA = factor(Region, levels = shape_keep$Region, labels = shape_keep$Region))
str(extrap_grid)

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "User", 
                          purpose = "index2", 
                          strata.limits = strata_use, 
                          bias.correct = FALSE,
                          knot_method = "grid")

# Run model
fit_g = fit_model_aja( settings = settings, 
                       input_grid = extrap_grid,
                       Lat_i = dat_use2[,'LATITUDE'], 
                       Lon_i = dat_use2[,'LONGITUDE'], 
                       t_i = dat_use2[,'YEAR'], 
                       c_i = rep(0, nrow(dat_use2)), 
                       b_i = dat_use2[,'CATCH_WT_CAL'], 
                       a_i = rep(0.0112 * 1.852^2, nrow(dat_use2)),
                       v_i = rep("missing", nrow(dat_use2)),
                       run_model = TRUE,
                       index_shapes = shape_keep)

plot(fit_g)


######
## Getting abundance index time series
######
get_vast_index_timeseries<- function(vast_fit, nice_category_names, index_scale = c("raw", "log"), out_dir){
  
  if(FALSE){
    vast_fit = fit_g
    index_scale = "raw"
  }
  
  TmbData<- vast_fit$data_list
  Sdreport<- vast_fit$parameter_estimates$SD
  
  # Time series steps
  time_ind<- 1:TmbData$n_t
  time_labels<- unique(vast_fit$data_frame$t_i)[time_ind]

  # Index regions
  index_regions_ind<- 1:TmbData$n_l
  index_regions = names(vast_fit$settings$strata.limits)[index_regions_ind]
 
  # Categories
  categories_ind<- 1:TmbData$n_c
  
  # Get the index information
  SD<- TMB::summary.sdreport(Sdreport)
  SD_stderr<- TMB:::as.list.sdreport(Sdreport, what = "Std. Error", report = TRUE)
  SD_estimate<- TMB:::as.list.sdreport(Sdreport, what = "Estimate", report = TRUE)
  if(vast_fit$settings$bias.correct == TRUE && "unbiased" %in% names(Sdreport)){
    SD_estimate_biascorrect<- TMB:::as.list.sdreport(Sdreport, what = "Est.(bias.correct)", report = TRUE)
  }
  
  # Now, populate array with values
  Index_ctl = log_Index_ctl = array(NA, dim = c(unlist(TmbData[c('n_c','n_t','n_l')]), 2), dimnames = list(categories_ind, time_labels, index_regions, c('Estimate','Std. Error')))
  
  if(index_scale == "raw"){
    if(vast_fit$settings$bias.correct == TRUE && "unbiased" %in% names(Sdreport)){
      Index_ctl[] = SD[which(rownames(SD) == "Index_ctl"),c('Est. (bias.correct)','Std. Error')]
    } else {
      Index_ctl[]<- SD[which(rownames(SD) == "Index_ctl"), c('Estimate','Std. Error')]
    }
    index_res_array<- Index_ctl
  } else {
    if(vast_fit$settings$bias.correct == TRUE && "unbiased" %in% names(Sdreport)){
      log_Index_ctl[] = SD[which(rownames(SD) == "ln_Index_ctl"),c('Est. (bias.correct)','Std. Error')]
    } else {
      log_Index_ctl[]<- SD[which(rownames(SD) == "ln_Index_ctl"), c('Estimate','Std. Error')]
    }
    index_res_array<- log_Index_ctl
  }
  
  # Data manipulation to get out out the array and to something more "plottable"
  for(i in seq_along(categories_ind)){
    index_array_temp<- index_res_array[i, , , ]
    index_res_temp_est<- data.frame("Time" = as.numeric(rownames(index_array_temp[,,1])), "Category" = categories_ind[i], index_array_temp[,,1]) %>%
      gather(., "Index_Region", "Index_Estimate", -Time, -Category)
    index_res_temp_sd<- data.frame("Time" = as.numeric(rownames(index_array_temp[,,1])), "Category" = categories_ind[i], index_array_temp[,,2]) %>%
      gather(., "Index_Region", "Index_SD", -Time, -Category)
    index_res_temp_out<- index_res_temp_est %>%
      left_join(., index_res_temp_sd)
    
    if(i == 1){
      index_res_out<- index_res_temp_out
    } else {
      index_res_out<- bind_rows(index_res_out, index_res_temp_out)
    }
  }
  
  # Save and return it
  write.csv(index_res_out, file = paste(out_dir, "Biomass_Index_", index_scale, "_", nice_category_names, ".csv", sep = ""))
  return(index_res_out)
}

plot_vast_index_timeseries<- function(index_res_df, nice_category_names, nice_xlab, nice_ylab, paneling = c("Category", "Index_Region", "None"), color_pal = c('#66c2a5','#fc8d62','#8da0cb'), out_dir){
  
  if(FALSE){
    index_res_df<- index_res_out
    nice_category_names<- "American lobster"
    nice_xlab = "Year"
    nice_ylab = "Biomass index (metric tons)"
    paneling<- "None"
  }
  
  if(paneling == "none"){
    if(!is.null(color_pal)){
      colors_use<- color_pal
    } else {
      color_pal<- c('#66c2a5','#fc8d62','#8da0cb','#e78ac3','#a6d854')
      colors_use<- color_pal[1:length(unique(index_res_df$Index_Region))]
    }
    plot_out<- ggplot() +
      geom_errorbar(data = index_res_df, aes(x = Time, ymin = Index_Estimate - Index_SD, ymax = Index_Estimate + Index_SD, color = Index_Region)) + 
      geom_point(data = index_res_df, aes(x = Time, y = Index_Estimate, color = Index_Region)) +
      scale_color_manual(values = colors_use) +
      xlab({{nice_xlab}}) +
      ylab({{nice_ylab}}) +
      ggtitle({{nice_category_names}}) + 
      theme_bw() +
      theme(legend.title = element_blank())
  }
  
  # Save and return the plot
  ggsave(plot_out, file = paste(out_dir, "Biomass_Index_", index_scale, "_", nice_category_names, ".jpg", sep = ""))
}
  
  
  
  if(is.null(year_labels)) year_labels = 1:TmbData$n_t
  if(is.null(years_to_plot)) years_to_plot = 1:TmbData$n_t
  if(is.null(strata_names)) strata_names = 1:TmbData$n_l
  if(is.null(category_names)) category_names = 1:TmbData$n_c
  
  
  Index_ctl = log_Index_ctl = array(NA, dim = c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames = list(category_names, year_labels, strata_names, c('Estimate','Std. Error')) )
  
}

TmbData
Sdreport

Index_ctl[] = SD[which(rownames(SD)==ParName),c('Est. (bias.correct)','Std. Error')]
Index_ctl[] = SD[which(rownames(SD)==ParName),c('Estimate','Std. Error')]

plot_index( Index_ctl=array(Index_ctl[,,,'Estimate'],dim(Index_ctl)[1:3]),
            sd_Index_ctl=array(log_Index_ctl[,,,'Std. Error'],dim(log_Index_ctl)[1:3]),
            year_labels=year_labels,
            years_to_plot=years_to_plot,
            strata_names=strata_names,
            category_names=category_names,
            DirName=DirName,
            PlotName=paste0(PlotName,"-",Plot_suffix[plotI],".png"),
            interval_width=interval_width,
            width=width,
            height=height,
            xlab="Year",
            ylab="Index",
            scale="log",
            plot_args=list("log"=ifelse(plot_log==TRUE,"y","")),
            Yrange=Yrange )





######
## Thinking towards editing fitted object with new data....
######
strata_use_orig<- load_example(data_set = "GB_fall_haddock")$strata.limits
print(strata_use_orig)

strata_use<- list("Georges_Bank_A" = strata_use_orig[[1]][1:5], "Georges_Bank_B" = strata_use_orig[[1]][6:15], "All" = strata_use_orig[[1]][1:15])

# We'd first fit the model, no adjustments...
# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "northwest_atlantic", 
                          purpose = "index2",  
                          bias.correct = FALSE,
                          knot_method = "grid", 
                          strata.limits = strata_use)

# Run model
fit_base = fit_model( settings = settings, 
                       Lat_i = dat_use2[,'LATITUDE'], 
                       Lon_i = dat_use2[,'LONGITUDE'], 
                       t_i = dat_use2[,'YEAR'], 
                       c_i = rep(0, nrow(dat_use2)), 
                       b_i = dat_use2[,'CATCH_WT_CAL'], 
                       a_i = rep(0.0112 * 1.852^2, nrow(dat_use2)),
                       v_i = rep("missing", nrow(dat_use2)),
                       run_model = FALSE)

names(fit_base)
names(fit_g)

str(fit_base$extrapolation_list$a_el)
str(fit_g$extrapolation_list$a_el)

names(fit_base$data_list)
names(fit_g$data_list)

str(fit_base$data_list$a_gl)
str(fit_g$data_list$a_gl)


# What will we need to edit? What if all we did was change the extrapolation object?
extrapolation_list_adjust<- Prepare_User_Extrapolation_Data_Fn_aja(input_grid = extrap_grid, strata.limits = strata_use, projargs = fit_base$extrapolation_list$projargs, flip_around_dateline =fit_base$extrapolation_list$flip_around_dateline, index_shapes = shape_keep)

# Now refit...
fit_new<- fit_model(settings = settings, 
                    input_grid = extrap_grid,
                    Lat_i = dat_use2[,'LATITUDE'], 
                    Lon_i = dat_use2[,'LONGITUDE'], 
                    t_i = dat_use2[,'YEAR'], 
                    c_i = rep(0, nrow(dat_use2)), 
                    b_i = dat_use2[,'CATCH_WT_CAL'], 
                    a_i = rep(0.0112 * 1.852^2, nrow(dat_use2)),
                    v_i = rep("missing", nrow(dat_use2)),
                    run_model = FALSE,
                    "spatial_args_input$Extrapolation_List" = extrapolation_list_adjust)

names(fit_base$data_list)
str(fit_base$data_list$Ais_ij)
