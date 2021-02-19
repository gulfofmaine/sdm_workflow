static_extract<- function(rast_stack, stack_names, sf_points, out_path, df_sf = "sf"){
  ## Details
  # This function reads a raster layer or stack and a sf spatial points object. The function then does a simple extraction by overlaying the spatial points onto the raster stack and gathering the value for each raster layer in the stack at all of the point locations. The function then returns the original sf_points data set as either an sf object or data frame with columns for the raster layers as a tibble and also saves this file as "model_covs.rds".
  
  # Args:
  # rast_stack = Static raster layer or stack of static raster layers 
  # stack_names = The names we want to use for the extracted values for each of the layers. These MUST be in the same order of the raster layers in the stack!!!
  # sf_points = SF spatial points object specifying the locations where we want to extract raster layer values
  # out_path = Path to save processed rds data file.
  # df_sf = Character string one of "sf" or "df" signaling whether the returned object should be a sf object or data frame
  
  # Returns: Either an SF object or data frame, which is also saved as an .rds file 
  
  ## Start function
  # Preliminaries -----------------------------------------------------------
  # Library check helper function -- not sure the best place to have this?
  library_check<- function(libraries) {
    ## Details
    # This function will check for and then either load or install any libraries needed to run subsequent functions
    
    # Args:
    # libraries = Vector of required library names
    
    # Returns: NA, just downloads/loads libraries into current work space.
    
    ## Start function
    lapply(libraries, FUN = function(x) {
      if (!require(x, character.only = TRUE)) {
        install.packages(x, dependencies = TRUE)
        library(x, character.only = TRUE)
      }
    })
    ## End function
  }
  
  # Load libraries, using library_check to download if it is missing
  libraries_needed<- c("tidyverse", "gmRi", "raster", "sf")
  library_check(libraries_needed)
  
  # For debugging
  if(FALSE){
    rast_stack = raster::stack(paste(shared.path(os.use = os_use, group = "root", folder = "RES Data/Shapefiles/"), "NEShelf_Etopo1_bathy.tiff", sep = ""))
    stack_names = "DEPTH"
    sf_points = trawl_sf
    out_path = here::here("/scratch/aja/data/")
    df_sf = "sf"
  }
  
  # A few checks...
  # Do the projections match?
  if(!st_crs(rast_stack) == st_crs(sf_points)){
    print("Check that projection in raster stack and spatial points match")
    stop()
  }
  
  # Are there duplicate records?
  if(any(duplicated(sf_points, by = c(EST_DATE, geometry)))){
    print("Check `sf_points` and remove duplicated observations to reduce extraction time")
    stop()
  }
 
  # Extraction -----------------------------------------------------------
  sf_extract<- data.frame(raster::extract(rast_stack, sf_points))
  names(sf_extract)<- stack_names
  
  # Bind to sf_unique 
  sf_points<- bind_cols(sf_points, sf_extract)
 
  # Write out and return processed file -----------------------------------------------------------
  if(df_sf == "sf"){
    # Keep it as sf object
    out<- sf_points
    saveRDS(out, file = paste(out_path, "model_covs.rds", sep = ""))
    return(out)
  } else {
    # Drop the geometry and save the data frame
    out<- st_drop_geometry(sf_points)
    saveRDS(out, file = paste(out_path, "model_covs.rds", sep = ""))
    return(out)
  }
  
 
  #########
  ## End
  #########
}
