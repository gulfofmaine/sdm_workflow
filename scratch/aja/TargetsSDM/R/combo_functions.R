####  Common Resources  ####


####  Functions  ####
####
## Dynamic loading functions -- each takes in a directory path and reads in files from it. 
species_read_csv<- function(species_table_dir){
  # Read it
  species_table<- read.csv(paste(species_table_dir, "species_table.csv", sep = "/"))
  
  # Return it
  return(species_table)
}

land_read_sf<- function(land_sf_dir){
  land_sf<- st_read(land_sf_dir)
  return(land_sf)
}

## Cleaning and processing functions

# Function to convert points to sf object
points_to_sf<- function(points){
  sf_out<- st_as_sf(points, coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326, remove = FALSE)
  return(sf_out)
}

#' @title Bind NMFS and DFO tow data 
#' 
#' @description A short function that will row bind a NMFS tows dataframe and DFO tows dataframe to create one dataframe that has all tows combined, with one row per unique tow.
#'
#' @param nmfs_tows = NMFS tows dataframe, created by `nmfs_get_tows()`
#' @param dfo_tows = DFO tows dataframe, created by `dfo_get_tows()`
#' @param out_dir = Directory to save the combined dataframe as an .rds file
#' 
#' @return A datafame with information of all unique tows. This file is also saved in out_dir. 
#' 
#' @export

bind_nmfs_dfo_tows<- function(nmfs_tows, dfo_tows, out_dir){
  
  # For debugging
  if(FALSE){
    nmfs_tows = readRDS(here::here("scratch/aja/TargetsSDM/data/nmfs/clean/nmfs_tows.rds"))
    dfo_tows = readRDS(here::here("scratch/aja/TargetsSDM/data/dfo/clean/dfo_tows.rds"))
    out_dir = here::here("scratch/aja/TargetsSDM/data/combined")
  }
  
  # Load in both datasets, add survey column
  nmfs_tows<- nmfs_tows %>%
    mutate(., "SURVEY" = rep("NMFS", nrow(.))) 
  nmfs_tows$SEASON<- toupper(nmfs_tows$SEASON)
  dfo_tows<- dfo_tows %>%
    mutate(., "SURVEY" = rep("DFO", nrow(.))) %>%
    dplyr::select(., -DIST)
  
  # Bind them together
  all_tows<- bind_rows(nmfs_tows, dfo_tows) %>%
    arrange(., DATE)
  
  # Return and save
  saveRDS(all_tows, file = paste(out_dir, "all_tows.rds", sep = "/"))
  return(all_tows)
}

#' @title Bind NMFS and DFO "tidy" occupancy data
#' 
#' @description A short function that will row bind a NMFS tidy occupancy dataframe with a DFO tidy occupancy dataframe to create one data frame that has all combined observation, where each row is a unique tow - species - occupancy record for all species in `species_table`.
#'
#' @param nmfs_tidy_occu = NMFS tidy occupancy dataframe, created by `nmfs_make_tidy_occu()`
#' @param dfo_tidy_occu = DFO tidy occupancy dataframe, created by `dfo_make_tidy_occu()`
#' @param out_dir = Directory to save the combined dataframe as an .rds file
#' 
#' @return A datafame with information of all occupancy records. This file is also saved in out_dir. 
#' 
#' @export

bind_nmfs_dfo_tidy_occu<- function(nmfs_tidy_occu, dfo_tidy_occu, out_dir){
  
  # For debugging
  if(FALSE){
    nmfs_tidy_occu = readRDS(here::here("scratch/aja/TargetsSDM/data/nmfs/clean/nmfs_tidy_occu.rds"))
    dfo_tidy_occu = readRDS(here::here("scratch/aja/TargetsSDM/data/dfo/clean/dfo_tidy_occu.rds"))
    out_dir = here::here("scratch/aja/TargetsSDM/data/combined")
  }
  
  # Load in both datasets, add survey column
  nmfs_tidy_occu<- nmfs_tidy_occu %>%
    mutate(., "SURVEY" = rep("NMFS", nrow(.)))
  dfo_tidy_occu<- dfo_tidy_occu %>%
    mutate(., "SURVEY" = rep("DFO", nrow(.)))
  
  # Bind them together
  all_tidy_occu<- bind_rows(nmfs_tidy_occu, dfo_tidy_occu)
  
  # Return and save
  saveRDS(all_tidy_occu, file = paste(out_dir, "all_tidy_occu.rds", sep = "/"))
  return(all_tidy_occu)
}

#' @title Make a "tidy" model dataframe
#' 
#' @description This function brings together information on the unique tows, which will usually include measures of effort as well as habitat covariates, with the occupancy data to create a "tidy" model dataframe. The "tidy" model dataframe is set up so that each row is a unique tow - species - occupancy record with additional information for habitat covariates that we might include in a species distribution model. 
#'
#' @param all_tows = A dataframe with information on all unique tows, including any potential habitat covariates describing characteristics of these locations
#' @param all_tidy_occu = A dataframe with all occupancy records, created by `bind_nmfs_dfo_tidy_occu`
#' @param out_dir = Directory to save the tidy model dataframe as an .rds file
#' 
#' @return A tidy model datafame with all the information (tows, habitat covariates, species occurrences) needed to fit a species distribution model. This file is also saved in out_dir. 
#' 
#' @export
make_tidy_mod_data<- function(all_tows, all_tidy_occu, out_dir){
  
  # For debugging
  if(FALSE){
    tar_load(all_tows_with_all_covs)
    all_tows = all_tows_with_all_covs
    tar_load(all_tidy_occu)
    all_tows = readRDS(here::here("scratch/aja/TargetsSDM/data/combined/all_tows.rds"))
    all_tidy_occu = readRDS(here::here("scratch/aja/TargetsSDM/data/combined/all_tidy_occu.rds"))
    out_dir = here::here("scratch/aja/TargetsSDM/data/combined")
  }
  
  # Need to join up the tow info with the tidy occu
  tidy_mod_data<- all_tidy_occu %>%
    left_join(., all_tows, by = c("ID", "SURVEY"))
  
  # Keep only what we need..
  cov_names<- names(all_tows)[-which(names(all_tows) %in% c("ID", "DATE", "EST_YEAR", "SEASON", "SURVEY", "SVVESSEL", "DECDEG_BEGLAT", "DECDEG_BEGLON", "NMFS_SVSPP", "DFO_SPEC", "PRESENCE", "BIOMASS", "ABUNDANCE"))]
  tidy_mod_data_out<- tidy_mod_data %>%
    dplyr::select(., ID, DATE, EST_YEAR, SEASON, SURVEY, SVVESSEL, DECDEG_BEGLAT, DECDEG_BEGLON, NMFS_SVSPP, DFO_SPEC, PRESENCE, BIOMASS, ABUNDANCE, {{cov_names}}) 
  
  # Drop NAs, return and save
  tidy_mod_data_out<- tidy_mod_data_out %>%
    drop_na(., {{cov_names}})
  saveRDS(tidy_mod_data_out, file = paste(out_dir, "tidy_mod_data.rds", sep = "/"))
  return(tidy_mod_data_out)
}

#' @title Create suite of trawl data summaries
#' 
#' @description This function... 
#'
#' @param all_tows = A dataframe with information on all unique tows, including any potential habitat covariates describing characteristics of these locations
#' @param out_dir = Directory to save the tidy model dataframe as an .rds file
#' 
#' @return Plots
#' 
#' @export
plot_trawl_dat_summs<- function(all_tows, fit_year_min, fit_year_max, grid, land_sf, land_color = "#d9d9d9",  xlim = c(-85, -55), ylim = c(30, 50), out_dir_maps, out_dir_tables){
  
  # For debugging
  if(FALSE){
    setwd(here::here("scratch/aja/TargetsSDM/"))
    tar_load(all_tows_with_all_covs)
    all_tows<- all_tows_with_all_covs
    tar_load(all_tidy_occu)
    tar_load(land_sf)
    land_color = "#d9d9d9"
    xlim = c(-77, -58)
    ylim = c(35, 47)
    fit_year_min<- 1985
    fit_year_max<- 2017
    grid = rotate(raster(here::here("scratch/aja/TargetsSDM/data/supporting", "Rast0.25grid.grd")))
    out_dir_maps = here::here("scratch/aja/TargetsSDM/results/plots_maps")
    out_dir_tables = here::here("scratch/aja/TargetsSDM/results/tables")
    
  }
  
  # Filter to focus years
  trawl_dat_temp<- all_tows %>% 
    filter(., EST_YEAR >= fit_year_min & EST_YEAR <= fit_year_max)
  
  # Timing of seasonal surveys...
  trawl_dat_dates_table<- trawl_dat_temp %>%
    group_by(., SURVEY, SEASON) %>%
    summarize(., "Min_Date" = min(DATE),
              "Max_Date" = max(DATE))
  
  # Summary table and plot
  trawl_dat_summ_table<- trawl_dat_temp %>%
    group_by(., SURVEY, SEASON, EST_YEAR) %>%
    summarize(., "Total_Tows_Season_Year" = n())
  
  trawl_dat_summ_plot<- ggplot() +
    geom_point(data = trawl_dat_summ_table, aes(x = EST_YEAR, y = Total_Tows_Season_Year, color = SURVEY, shape = SEASON)) +
    geom_line() +
    scale_color_manual(name = "Survey", values = c("#7570b3", "#1b9e77")) +
    scale_shape_manual(name = "Season", values = c(1, 2, 3)) +
    theme_bw()
  
  # Spatial map for fall and spring to explore potential overlap...
  trawl_map_dat_fall<- trawl_dat_temp %>%
    filter(., SEASON %in% c("FALL"))
  
  trawl_map_dat_spring<- trawl_dat_temp %>%
    filter(., SEASON %in% c("SPRING"))
  
  fall_maps<- ggplot() +
    geom_point(data = trawl_map_dat_fall, aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = SURVEY, shape = SURVEY), size = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Survey", values = c("#7570b3", "#1b9e77")) +
    scale_shape_manual(name = "Survey", values = c(1, 2)) +
    geom_sf(data = land_sf, fill = land_color, lwd = 0.2) +
    coord_sf(xlim =  xlim, ylim = ylim, expand = FALSE) +
  theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05)) +
    facet_wrap(~EST_YEAR)
  
  spring_maps<-  ggplot() +
    geom_point(data = trawl_map_dat_spring, aes(x = DECDEG_BEGLON, y = DECDEG_BEGLAT, color = SURVEY, shape = SURVEY), size = 0.5, alpha = 0.5) +
    scale_color_manual(name = "Survey", values = c("#7570b3", "#1b9e77")) +
    scale_shape_manual(name = "Survey", values = c(1, 2)) +
    geom_sf(data = land_sf, fill = land_color, lwd = 0.2) +
    coord_sf(xlim =  xlim, ylim = ylim, expand = FALSE) +
    theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05)) +
    facet_wrap(~EST_YEAR)
  
  # Samples per grid cell in the spring and fall...
  samp_rasts_fall<- raster::stack()
  samp_rasts_spring<- raster::stack()

  for(i in seq_along(unique(trawl_map_dat$EST_YEAR))){
    # Get the data
    temp_dat_fall<- trawl_map_dat_fall %>%
      dplyr::filter(., EST_YEAR == unique(trawl_map_dat$EST_YEAR)[i]) %>%
      st_as_sf(., coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326, remove = FALSE)
    temp_dat_spring<- trawl_map_dat_spring %>%
      dplyr::filter(., EST_YEAR == unique(trawl_map_dat$EST_YEAR)[i]) %>%
      st_as_sf(., coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326, remove = FALSE)
    
    # Rasterize and save in stack
    samp_rast_fall_temp<- rasterize(temp_dat_fall, y = grid, field = "SURVEY", fun = function(x,...) length(unique(na.omit(x))))
    samp_rast_spring_temp<- rasterize(temp_dat_spring, y = grid, field = "SURVEY", fun = function(x,...) length(unique(na.omit(x))))
    
    samp_rasts_fall<- raster::stack(samp_rasts_fall, samp_rast_fall_temp)
    samp_rasts_spring<- raster::stack(samp_rasts_spring, samp_rast_spring_temp)
  }
  
  # Plot...
  samp_rasts_fall_df<- as.data.frame(samp_rasts_fall, xy = TRUE)
  colnames(samp_rasts_fall_df)[3:ncol(samp_rasts_fall_df)]<- paste("Surveys", unique(trawl_map_dat$EST_YEAR), sep = "_")
  samp_rasts_fall_df<- samp_rasts_fall_df %>%
    pivot_longer(., cols = !c(x, y), names_to = "Survey", values_to = "Count") %>%
    drop_na(Count) %>%
    mutate(., "Count_Fac" = factor(Count, levels = c(1, 2)))
  samp_rasts_spring_df<- as.data.frame(samp_rasts_spring, xy = TRUE)
  colnames(samp_rasts_spring_df)[3:ncol(samp_rasts_spring_df)]<- paste("Surveys", unique(trawl_map_dat$EST_YEAR), sep = "_")
  samp_rasts_spring_df<- samp_rasts_spring_df %>%
    pivot_longer(., cols = !c(x, y), names_to = "Survey", values_to = "Count") %>%
    drop_na(Count) %>%
    mutate(., "Count_Fac" = factor(Count, levels = c(1, 2)))
  
  # Plot...
  fall_plot<- ggplot() +
    geom_tile(data = samp_rasts_fall_df, aes(x = x, y = y, fill = Count_Fac, alpha = Count_Fac)) +
    scale_fill_manual(name = "Unique survey samples", values = c("#1b9e77", "#7570b3")) +
    scale_alpha_manual(name = "Unique survey samples", values = c(0.25, 1)) +
    geom_sf(data = land_sf, fill = land_color, lwd = 0.2) +
    coord_sf(xlim =  xlim, ylim = ylim, expand = FALSE) +
    theme(panel.background = element_rect(fill = "white"), strip.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05)) +
    facet_wrap(~Survey, ncol = 5)
  ggsave(paste0(out_dir_maps, "/Fall_Samples.jpg"), fall_plot, width = 10, height = 15)
  
  spring_plot<- ggplot() +
    geom_tile(data = samp_rasts_spring_df, aes(x = x, y = y, fill = Count_Fac, alpha = Count_Fac)) +
    scale_fill_manual(name = "Unique survey samples", values = c("#1b9e77", "#7570b3")) +
    scale_alpha_manual(name = "Unique survey samples", values = c(0.25, 1)) +
    geom_sf(data = land_sf, fill = land_color, lwd = 0.2) +
    coord_sf(xlim =  xlim, ylim = ylim, expand = FALSE) +
    theme(panel.background = element_rect(fill = "white"), strip.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05)) +
    facet_wrap(~Survey, ncol = 5)
  ggsave(paste0(out_dir_maps, "/Spring_Samples.jpg"), spring_plot, width = 10, height = 15)
  
  # For the spring season, can we visualize the catches from each survey in overlapping grid points?
  # First need to get the grid cell each point is within...
  spring_sf<- trawl_map_dat_spring %>%
    st_as_sf(., coords = c("DECDEG_BEGLON", "DECDEG_BEGLAT"), crs = 4326, remove = FALSE) %>%
    mutate(., "Raster_Cell" = cellFromXY_sf(grid, .))
  
  # Get duplicated "Raster_Cell" columns...
  spring_keep_NMFS<- spring_sf %>%
    filter(., SURVEY == "NMFS") %>%
    group_by(., EST_YEAR, SEASON) %>%
    mutate(., "Year_Raster_Cell" = paste(EST_YEAR, Raster_Cell, sep = "_")) %>%
    dplyr::select(., ID, EST_YEAR, SEASON, SURVEY, Year_Raster_Cell) %>%
    st_drop_geometry()
  
  spring_keep_DFO<- spring_sf %>%
    filter(., SURVEY == "DFO") %>%
    group_by(., EST_YEAR, SEASON) %>%
    mutate(., "Year_Raster_Cell" = paste(EST_YEAR, Raster_Cell, sep = "_")) %>%
    dplyr::select(., ID, EST_YEAR, SEASON, SURVEY, Year_Raster_Cell) %>%
    st_drop_geometry()
  
  spring_keep_NMFS_IDs<- spring_keep_NMFS$ID[which(spring_keep_NMFS$Year_Raster_Cell %in% spring_keep_DFO$Year_Raster_Cell)]
  spring_keep_NMFS<-spring_keep_NMFS %>%
    filter(., ID %in% spring_keep_NMFS_IDs) %>%
    ungroup() %>%
    dplyr::select(., ID, Year_Raster_Cell)
  
  spring_keep_DFO_IDs<- spring_keep_DFO$ID[which(spring_keep_DFO$Year_Raster_Cell %in% spring_keep_NMFS$Year_Raster_Cell)]
  spring_keep_DFO<- spring_keep_DFO %>%
    filter(., ID %in% spring_keep_DFO_IDs) %>%
    ungroup() %>%
    dplyr::select(., ID, Year_Raster_Cell)
  
  # Bring in catch data...
  spring_catch_comp_NMFS<- all_tidy_occu %>%
    filter(., ID %in% spring_keep_NMFS$ID) %>%
    rename(., c(NMFS_Presence = PRESENCE, NMFS_Biomass = BIOMASS)) %>%
    left_join(., spring_keep_NMFS, by = c("ID" = "ID")) %>%
    dplyr::select(Year_Raster_Cell, NMFS_SVSPP, NMFS_Presence, NMFS_Biomass)
  
  spring_catch_comp_DFO<- all_tidy_occu %>%
    filter(., ID %in% spring_keep_DFO$ID) %>%
    rename(., c(DFO_Presence = PRESENCE, DFO_Biomass = BIOMASS)) %>%
    left_join(., spring_keep_DFO, by = c("ID" = "ID")) %>%
    dplyr::select(Year_Raster_Cell, NMFS_SVSPP, DFO_Presence, DFO_Biomass)
  
  spring_catch_comp<- spring_catch_comp_NMFS %>%
    left_join(., spring_catch_comp_DFO, by = c("NMFS_SVSPP" = "NMFS_SVSPP", "Year_Raster_Cell" = "Year_Raster_Cell")) %>%
    mutate(., "Sum_Biomass" = NMFS_Biomass + DFO_Biomass) %>%
    filter(., Sum_Biomass > 0)
  
  # Which species might have multiple comparisons?
  spring_catch_comp_table<- spring_catch_comp %>%
    group_by(., NMFS_SVSPP) %>%
    summarize(., "Comparison_Samples" = n()) %>%
    filter(., Comparison_Samples > 500)
  
  spring_catch_comp_plot<- spring_catch_comp %>%
    filter(., NMFS_SVSPP %in% spring_catch_comp_table$NMFS_SVSPP)
  
  spring_catch_comp_plot_out<- ggplot() +
    geom_point(data = spring_catch_comp_plot, aes(x = NMFS_Biomass, y = DFO_Biomass)) +
    theme_bw() +
    facet_wrap(~NMFS_SVSPP, scales = "free")
}


cellFromXY_sf <- function(object, xy) {
  if (inherits(xy, "sf")) {
    xy <- spbabel::sptable(xy) %>% dplyr::select(x_, y_)
    cellFromXY(object, as.matrix(xy))
  } else {
    raster::cellFromXY(object, xy)
  }
}