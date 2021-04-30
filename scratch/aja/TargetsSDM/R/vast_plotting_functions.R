vast_plot_density<- function(vast_fit, locations, all_times = all_times, plot_times = NULL, land_sf, xlim, ylim, panel_or_gif = "gif", out_dir, ...){
  if(FALSE){
    vast_fit = fit
    locations = "grid"
    all_times = as.character(unique(full_data$year_season))
    plot_times = NULL
    land_color = "#d9d9d9"
    res_data_path = "~/Box/RES_Data/"
    land_sf = st_read(paste(res_data_path, "Shapefiles/ne_50m_land/ne_50m_land.shp", sep = ""))
    xlim = c(-85, -60)
    ylim = c(25, 48)
    panel_or_gif = "gif"
    panel_cols 
    panel_rows
  }
  
  # Get prediction dataframe together
  if(locations == "points"){
    # Plotting prediction points
  } else {
    # Plotting at spatial knots...
    # Getting prediction array
    pred_array<- log(vast_fit$Report$D_gct+1)
    
    # Getting time info
    if(!is.null(plot_times)){
      plot_times<- all_times[which(all_times) %in% plot_times]
    } else {
      plot_times<- all_times
    }
    
    # Getting spatial information
    spat_data<- vast_fit$extrapolation_list
    loc_g<- spat_data$Data_Extrap[which(spat_data$Data_Extrap[, "Include"] > 0), c("Lon", "Lat")]
    CRS_orig<- sp::CRS("+proj=longlat")
    CRS_proj<- sp::CRS(spat_data$projargs)
    land_sf_proj<- st_transform(st_crop(land_sf, xmin = xlim[1], ymin = ylim[1], xmax = xlim[2], ymax = ylim[2]),  crs = CRS_proj)
    
    # Looping through...
    rasts_out<- vector("list", dim(pred_array)[3])
    rasts_range<- pred_array
    rast_lims<- c(round(min(rasts_range)-0.000001, 2), round(max(rasts_range) + 0.0000001, 2))
    
    if(dim(pred_array)[3] == 1){
      df<- data.frame(loc_g, z = pred_array[,1,])
      points_ll = st_as_sf(data_df, coords = c("Lon", "Lat"), crs = CRS_orig)
      points_proj = points_ll %>%
        st_transform(., crs = CRS_proj)
      points_bbox<- st_bbox(points_proj)
      raster_proj<- st_rasterize(points_proj)
      
      plot_out<- ggplot() +
        geom_stars(data = raster_proj, aes(x = x, y = y, fill = z)) +
        scale_fill_viridis_c(name = "Density", option = "viridis", na.value = "transparent", limits = rast_lims) +
        geom_sf(data = land_sf_proj, fill = land_color, lwd = 0.2) +
        coord_sf(xlim = points_bbox[c(1,3)], ylim = points_bbox[c(2,4)], expand = FALSE, datum = sf::st_crs(CRS_proj)) 
        theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05))
      
      ggsave(filename = paste(out_dir, file_name, ".png", sep = ""), plot_out, width = 11, height = 8, units = "in")
    } else {
    
      for (tI in 1:dim(pred_array)[3]) {
        data_df<- data.frame(loc_g, z = pred_array[,1,tI])
        points_ll = st_as_sf(data_df, coords = c("Lon", "Lat"), crs = CRS_orig)
        points_proj = points_ll %>%
          st_transform(., crs = CRS_proj)
        points_bbox<- st_bbox(points_proj)
        n_cells<- nrow(data_df)
        cell_size = mean(diff(points_bbox[c(1,3)]), diff(points_bbox[c(2,4)]))/floor(sqrt(n_cells))
        
        raster_proj<- as.data.frame(plotKML::vect2rast(as_Spatial(points_proj), cell.size = cell_size, fun = mean))
        names(raster_proj)[2:3]<- c("x", "y")
        
        time_plot_use<- plot_times[tI]
        
        rasts_out[[tI]]<- ggplot() +
          geom_tile(data = raster_proj, aes(x = x, y = y, fill = z)) +
          scale_fill_viridis_c(name = "Log(density+1)", option = "viridis", na.value = "transparent", limits = rast_lims) +
          geom_sf(data = land_sf_proj, fill = land_color, lwd = 0.2, na.rm = TRUE) +
          coord_sf(xlim = points_bbox[c(1,3)], ylim = points_bbox[c(2,4)], expand = FALSE, datum = sf::st_crs(CRS_proj)) +
          #annotate("text", x = points_bbox[3]-175000, y = points_bbox[2] + 100000, label = tolower(time_plot_use)) +
          ggtitle(time_plot_use)
          theme(panel.background = element_rect(fill = "white"), panel.border = element_rect(fill = NA), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title = element_blank(), plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt"))
      }
      if(panel_or_gif == "panel"){
        # Panel plot
        all_plot<- wrap_plots(rasts_out, ncol = panel_cols, nrow = panel_rows, guides = "collect", theme(plot.margin = margin(t = 0.05, r = 0.05, b = 0.05, l = 0.05, unit = "pt")))
        ggsave(filename = paste(working_dir, file_name, ".png", sep = ""), all.plot, width = 11, height = 8, units = "in")
      } else {
        # Make a gif
        plot_loop_func<- function(plot_list){
          for (i in seq_along(plot_list)) {
            plot_use<- plot_list[[i]]
            print(plot_use)
          }
        }
        invisible(save_gif(plot_loop_func(rasts_out), paste0(out_dir, "LogDensity.gif"), delay = 0.5, width = 600, height = 800, progress = FALSE))
      }
    }
  }
}

