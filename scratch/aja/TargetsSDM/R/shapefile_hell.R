######
## Code to figure out the different regions we will use to calculate the index of abundance...
######
nelme<- st_read("~/Box/RES_Data/Shapefiles/large_marine_ecosystems/northeast_us_continental_shelf_exterior.geojson")
plot(st_geometry(nelme))

sshelf<- st_read("~/Box/RES_Data/Shapefiles/large_marine_ecosystems/scotian_shelf_exterior.geojson")
plot(st_geometry(sshelf))

full<- bind_rows(nelme, sshelf) 
st_write(full, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_polys.shp")

full_region<- full %>%
  st_union()
st_write(full_region, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/region_shape.shp")

dfo<- st_read("~/Box/Mills Lab/Projects/DFO_survey_data/strata_shapefiles/MaritimesRegionEcosystemAssessmentBoundary.shp")
plot(st_geometry(dfo))

eez<- st_read("~/Box/RES_Data/Shapefiles/World_EEZ_v11_20191118/eez_v11.shp")
str(eez)

us_ca<- eez %>% 
  filter(., GEONAME %in% c("United States Exclusive Economic Zone", "Canadian Exclusive Economic Zone"))

plot(st_geometry(full_region))
plot(st_geometry(dfo), col = "red", add = TRUE)
plot(st_geometry(us_ca), col = "blue", add = TRUE)
plot(st_geometry(dfo), col = "red", add = TRUE)

us_eez<- us_ca %>%
  filter(., GEONAME == "United States Exclusive Economic Zone")
ca_eez<- us_ca %>%
  filter(., GEONAME == "Canadian Exclusive Economic Zone")

dfo_crop<- dfo %>%
  st_difference(., us_eez)
plot(st_geometry(dfo_crop), col = "green", add = TRUE)

dfo_crop<- st_as_sf(data.frame(Region = "Canada"), geometry=dfo_crop$geometry)

st_write(dfo_crop, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/maritimes_regions_canadian_only.shp", append = FALSE)

nmfs_us_only<- st_read("~/Box/RES_Data/Shapefiles/nmfs_trawl_regions/nmfs_trawl_regions_collection.geojson") %>%
  st_union() %>%
  st_transform(., crs = st_crs(dfo_crop))
nmfs_us_only<- nmfs_us_only %>%
  st_difference(., dfo_crop)
nmfs_us_only<- st_as_sf(data.frame(Region = "US"), geometry=nmfs_us_only)
nmfs_us_only2<- nmfs_us_only %>%
  st_difference(., ca_eez)
st_write(nmfs_us_only2, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/nmfs_us_only.shp", append = FALSE)

full_survey_region<- nmfs_us_only2 %>%
  bind_rows(., dfo_crop) %>%
  st_union()
full_survey_region<- st_as_sf(data.frame(Region = "NMFS_and_DFO"), geometry=full_survey_region)

st_write(full_survey_region, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/NMFS_DFO.shp", append = FALSE)


#####
## Status check...
#####

nmfs_us_only<- st_read("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/nmfs_us_only.shp") 
nmfs_us_only<- st_as_sf(data.frame(Region = "NMFS"), geometry=nmfs_us_only$geometry)
st_write(nmfs_us_only, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/NMFS.shp", append = FALSE)

dfo_ca_only<- st_read("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/maritimes_regions_canadian_only.shp")
dfo_ca_only<- st_as_sf(data.frame(Region = "DFO"), geometry=dfo_ca_only$geometry)
st_write(dfo_ca_only, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/DFO.shp", append = FALSE)

full_region<- st_read("~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/region_shapefile/full_survey_region.shp")
full_Region<- st_as_sf(data.frame(Region = "NMFS_and_DFO"), geometry=full_region$geometry)
st_write(full_region, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/NMFS_DFO.shp", append = FALSE)

plot(st_geometry(full_region), lwd = 4)
plot(st_geometry(nmfs_us_only), border = "red", add = TRUE, lwd = 2)
plot(st_geometry(dfo_ca_only), border = "blue", add = TRUE, lwd = 2)


## Looks good, now how about the biogeographic regions? Scotian Shelf - Gulf of Maine? - Georges Bank - SNE - MAB. For each of these, also want to "crop" them by our full region polygon...
sshelf<- st_read("~/Box/RES_Data/Shapefiles/large_marine_ecosystems/scotian_shelf_exterior.geojson") %>%
  st_transform(., crs = st_crs(dfo_crop))
plot(st_geometry(sshelf), border = "green", add = TRUE, lwd = 2)

gom<- st_read("~/Box/RES_Data/Shapefiles/nmfs_trawl_regions/nmfs_trawl_gulf_of_maine.geojson") %>% 
  st_union() %>%
  st_transform(., crs = st_crs(dfo_crop))
gom<- st_as_sf(data.frame(Region = "Gulf_of_Maine"), geometry=gom)
plot(st_geometry(gom), col = "white", border = "purple", add = TRUE, lwd = 2)
st_write(gom, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/gom.shp", append = FALSE)

sne<- st_read("~/Box/RES_Data/Shapefiles/nmfs_trawl_regions/nmfs_trawl_southern_new_england.geojson") %>% 
  st_union() %>%
  st_transform(., crs = st_crs(dfo_crop))
sne<- st_as_sf(data.frame(Region = "Southern_New_England"), geometry=sne)

plot(st_geometry(sne), fill = "white", border = "purple", add = TRUE, lwd = 2)
st_write(sne, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/sne.shp", append = FALSE)

gb<- st_read("~/Box/RES_Data/Shapefiles/nmfs_trawl_regions/nmfs_trawl_georges_bank.geojson") %>% 
  st_union() %>%
  st_transform(., crs = st_crs(dfo_crop))
gb<- st_as_sf(data.frame(Region = "Georges_Bank"), geometry=gb)
plot(st_geometry(gb), fill = "white", border = "purple", add = TRUE, lwd = 2)
st_write(gb, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/gb.shp", append = FALSE)

mab<- st_read("~/Box/RES_Data/Shapefiles/nmfs_trawl_regions/nmfs_trawl_mid_atlantic_bight.geojson") %>% 
  st_union() %>%
  st_transform(., crs = st_crs(dfo_crop))
mab<- st_as_sf(data.frame(Region = "Mid_Atlantic_Bight"), geometry=mab)

plot(st_geometry(mab), fill = "white", border = "purple", add = TRUE, lwd = 2)
st_write(mab, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/mab.shp", append = FALSE)

# Figuring out the SS area...
sshelf2<- dfo_ca_only %>%
  st_difference(., gom) %>%
  st_difference(., gb)
plot(st_geometry(sshelf2), border = "red", lwd =  2, add = TRUE)

bof<- st_read("~/Box/RES_Data/Shapefiles/GulfOfMainePhysioRegions/PhysioRegions_WGS84.shp") %>%
  filter(., Region == "Bay of Fundy")
sshelf3<- sshelf2 %>%
  st_difference(., bof)
sshelf<- st_as_sf(data.frame(Region = "Scotian_Shelf"), geometry=sshelf3$geometry) %>%
  st_zm(., drop = TRUE)
plot(st_geometry(sshelf), border = "blue", lwd = 2, add = TRUE)
st_write(sshelf, "~/GitHub/sdm_workflow/scratch/aja/TargetsSDM/data/supporting/index_shapefiles/sshelf.shp", append = FALSE)



%>%
  st_difference(., full_region)
plot(st_geometry(sshelf2), border = "green", lwd =  2, add = TRUE)








sshelf2<- sshelf %>%
  st_intersects(., dfo_ca_only)
plot(st_geometry(sshelf2), border = "green", add = TRUE, lwd = 2)



us_ca2<- us_ca %>%
  st_crop(., )


data(georges_bank_haddock_fall, package = "FishStatsUtils")
dat_use<- georges_bank_haddock_fall
strata_use<- load_example(data_set = "GB_fall_haddock")$strata.limits
print(strata_use)


# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "northwest_atlantic", 
                          purpose = "index2", 
                          strata.limits = strata_use, 
                          bias.correct = FALSE )

# Run model
fit = fit_model( settings = settings, 
                 Lat_i = dat_use[,'LATITUDE'], 
                 Lon_i = dat_use[,'LONGITUDE'], 
                 t_i = dat_use[,'YEAR'], 
                 c_i = rep(0, nrow(dat_use)), 
                 b_i = dat_use[,'CATCH_WT_CAL'], 
                 a_i = rep(0.0112 * 1.852^2, nrow(dat_use)),
                 v_i = rep("missing", nrow(dat_use)))

# Plot results
plot( fit )
colSums( fit$extrapolation_list$a_el )

#####
## Trying different strata -- easy way, change the numbers
#####
strata_use_orig<- load_example(data_set = "GB_fall_haddock")$strata.limits
print(strata_use_orig)

strata_use<- list("Georges_Bank" = strata_use_orig[[1]][1:5])


# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "northwest_atlantic", 
                          purpose = "index2", 
                          strata.limits = strata_use, 
                          bias.correct = FALSE )

# Run model
fit_b = fit_model( settings = settings, 
                 Lat_i = dat_use[,'LATITUDE'], 
                 Lon_i = dat_use[,'LONGITUDE'], 
                 t_i = dat_use[,'YEAR'], 
                 c_i = rep(0, nrow(dat_use)), 
                 b_i = dat_use[,'CATCH_WT_CAL'], 
                 a_i = rep(0.0112 * 1.852^2, nrow(dat_use)),
                 v_i = rep("missing", nrow(dat_use)))

# Plot results
plot( fit_b )
colSums( fit_b$extrapolation_list$a_el )

#####
## Try 2 -- that worked. Now how about adding them both together?
#####
strata_use<- list("Georges_Bank_A" = strata_use_orig[[1]][1:5], "Georges_Bank_B" = strata_use_orig[[1]][6:15])


# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "northwest_atlantic", 
                          purpose = "index2", 
                          strata.limits = strata_use, 
                          bias.correct = FALSE )

# Run model
fit_c = fit_model( settings = settings, 
                   Lat_i = dat_use[,'LATITUDE'], 
                   Lon_i = dat_use[,'LONGITUDE'], 
                   t_i = dat_use[,'YEAR'], 
                   c_i = rep(0, nrow(dat_use)), 
                   b_i = dat_use[,'CATCH_WT_CAL'], 
                   a_i = rep(0.0112 * 1.852^2, nrow(dat_use)),
                   v_i = rep("missing", nrow(dat_use)))

# Plot results
plot( fit_c )
colSums( fit_c$extrapolation_list$a_el )


#####
## That also worked, how about one with "all" the strata??
#####
strata_use<- list("Georges_Bank_A" = strata_use_orig[[1]][1:5], "Georges_Bank_B" = strata_use_orig[[1]][6:15], "All" = strata_use_orig[[1]][1:15])


# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "northwest_atlantic", 
                          purpose = "index2", 
                          strata.limits = strata_use, 
                          bias.correct = FALSE )

# Run model
fit_d = fit_model( settings = settings, 
                   Lat_i = dat_use[,'LATITUDE'], 
                   Lon_i = dat_use[,'LONGITUDE'], 
                   t_i = dat_use[,'YEAR'], 
                   c_i = rep(0, nrow(dat_use)), 
                   b_i = dat_use[,'CATCH_WT_CAL'], 
                   a_i = rep(0.0112 * 1.852^2, nrow(dat_use)),
                   v_i = rep("missing", nrow(dat_use)))

# Plot results
plot( fit_d )
colSums( fit_d$extrapolation_list$a_el )


#####
## Alright...all of that worked. What about polygons instead?
#####

# Read in trawl data polygon
nmfs<- st_read("~/Box/RES_Data/shapefiles/nmfs_trawl_regions/nmfs_trawl_regions_collection.geojson")
strata_data<- data.frame("Region" = "Georges_Bank")

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
# extrap_grid<- region_grid %>%
#   st_transform(., crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ") %>%
#   st_join(., shape_keep, join = st_within) %>%
#   mutate(., "Lon" = as.numeric(st_coordinates(.)[,1]),
#          "Lat" = as.numeric(st_coordinates(.)[,2])) %>%
#   st_drop_geometry() %>%
#   dplyr::select(., Lon, Lat, Region) %>%
#   mutate(., Area_km2=((cell_size/1000)^2),
#          STRATA = factor(Region, levels = c("Georges_Bank_A", "Georges_Bank_B", "All"), labels = c("Georges_Bank_A", "Georges_Bank_B", "All")))
extrap_grid<- region_grid %>%
  st_transform(., crs = "+proj=longlat +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0 ") %>%
  st_join(., shape_keep, join = st_within) %>%
  mutate(., "Lon" = as.numeric(st_coordinates(.)[,1]),
         "Lat" = as.numeric(st_coordinates(.)[,2])) %>%
  st_drop_geometry() %>%
  dplyr::select(., Lon, Lat, Region) %>%
  mutate(., Area_km2=((cell_size/1000)^2),
         STRATA = factor(Region, levels = c("Georges_Bank_A", "Georges_Bank_B", "All"), labels = c("Georges_Bank_A", "Georges_Bank_B", "All")))
str(extrap_grid)

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "User", 
                          purpose = "index2", 
                          strata.limits = strata_use, 
                          bias.correct = FALSE,
                          knot_method = "grid")

# Run model
fit_f = fit_model( settings = settings, 
                   input_grid = extrap_grid,
                   Lat_i = dat_use2[,'LATITUDE'], 
                   Lon_i = dat_use2[,'LONGITUDE'], 
                   t_i = dat_use2[,'YEAR'], 
                   c_i = rep(0, nrow(dat_use2)), 
                   b_i = dat_use2[,'CATCH_WT_CAL'], 
                   a_i = rep(0.0112 * 1.852^2, nrow(dat_use2)),
                   v_i = rep("missing", nrow(dat_use2)),
                   run_model = TRUE)

plot(fit_f)

# Not working....result is the same as the "all" (104k ish) and it should be around 26k


#####
## What happens when things are "working"
strata_use_orig<- load_example(data_set = "GB_fall_haddock")$strata.limits
print(strata_use_orig)

strata_use<- list("Georges_Bank" = strata_use_orig[[1]][1:5])

# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "northwest_atlantic", 
                          purpose = "index2", 
                          strata.limits = strata_use, 
                          bias.correct = FALSE )

Region = settings$Region
projargs = NA
zone = settings$zone
strata.limits = settings$strata.limits
create_strata_per_region = FALSE
max_cells = settings$max_cells
input_grid = extrap_grid
observations_LL = NULL
grid_dim_km = c(2, 2)
maximum_distance_from_sample = NULL
grid_in_UTM = TRUE
grid_dim_LL = c(0.1, 0.1)
region = c("south_coast", "west_coast")
strata_to_use = c("SOG", "WCVI", "QCS", "HS", "WCHG")
epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[1]
survey = "Chatham_rise"
surveyname = "propInWCGBTS"
flip_around_dateline
nstart = 100
area_tolerance = 0.05
backwards_compatible_kmeans = FALSE
DirPath = paste0(getwd(), "/") 

rI = 1
Extrapolation_List = NULL

if (missing(flip_around_dateline)) 
  flip_around_dateline = FALSE
strata.limits = strata.limits
epu_to_use = epu_to_use
projargs = projargs
zone = zone
flip_around_dateline = flip_around_dateline

message("Using strata ", strata.limits)
if (any(tolower(epu_to_use) %in% "all")) {
  epu_to_use <- c("Georges_Bank", "Mid_Atlantic_Bight", 
                  "Scotian_Shelf", "Gulf_of_Maine", "Other")
}
utils::data(northwest_atlantic_grid, package = "FishStatsUtils")
Data_Extrap <- northwest_atlantic_grid

# What is this???
Data_Extrap_sf<- Data_Extrap %>%
  st_as_sf(., coords = c("Lon", "Lat"), crs = st_crs(nmfs))
plot(st_geometry(Data_Extrap_sf)) #Each one of these knots has an area associated with it...

# How about our grid??
our_grid_sf<- extrap_grid %>%
  st_as_sf(., coords = c("Lon", "Lat"), crs = st_crs(nmfs))
plot(st_geometry(our_grid_sf), col = "red", add = TRUE)

# Other way??
plot(st_geometry(our_grid_sf), col = "red")
plot(st_geometry(Data_Extrap_sf), add = TRUE)

Tmp = cbind(BEST_DEPTH_M = 0, BEST_LAT_DD = Data_Extrap[, "Lat"], BEST_LON_DD = Data_Extrap[, "Lon"])
if (length(strata.limits) == 1 && strata.limits[1] == "EPU") {
  Data_Extrap <- Data_Extrap[Data_Extrap$EPU %in% epu_to_use, ]
  Data_Extrap$EPU <- droplevels(Data_Extrap$EPU)
  a_el = matrix(NA, nrow = nrow(Data_Extrap), ncol = length(epu_to_use), dimnames = list(NULL, epu_to_use))
  Area_km2_x = Data_Extrap[, "Area_in_survey_km2"]
  for (l in 1:ncol(a_el)) {
    a_el[, l] = ifelse(Data_Extrap[, "EPU"] %in% epu_to_use[[l]], Area_km2_x, 0)
  }
} else {
  a_el = as.data.frame(matrix(NA, nrow = nrow(Data_Extrap), ncol = length(strata.limits), dimnames = list(NULL, names(strata.limits))))
  Area_km2_x = Data_Extrap[, "Area_in_survey_km2"]
  ######
  ## THIS SEEMS TO BE THE KEY!!!
  for (l in 1:ncol(a_el)) {
    a_el[, l] = ifelse(Data_Extrap[, "stratum_number"] %in%  strata.limits[[l]], Area_km2_x, 0)
  }
}
tmpUTM = project_coordinates(X = Data_Extrap[, "Lon"], Y = Data_Extrap[, "Lat"], projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline)
Data_Extrap = cbind(Data_Extrap, Include = 1)
Data_Extrap[, c("E_km", "N_km")] = tmpUTM[, c("X", "Y")]
Return = list(a_el = a_el, Data_Extrap = Data_Extrap, zone = attr(tmpUTM, "zone"), projargs = attr(tmpUTM, "projargs"), flip_around_dateline = flip_around_dateline, Area_km2_x = Area_km2_x)


Extrapolation_List = Prepare_NWA_Extrapolation_Data_Fn(strata.limits = strata.limits, epu_to_use = epu_to_use, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline, ...)

if (max_cells < nrow(Return$Data_Extrap)) {
  message("# Reducing extrapolation-grid from ", nrow(Return$Data_Extrap), " to ", max_cells, " cells for Region(s): ", paste(Region, collapse = ", "))
  loc_orig = Return$Data_Extrap[, c("E_km", "N_km")]
  loc_orig = loc_orig[which(Return$Area_km2_x > 0), ]
  Kmeans = make_kmeans(n_x = max_cells, loc_orig = loc_orig, nstart = nstart, randomseed = 1, iter.max = 1000, DirPath = DirPath, Save_Results = TRUE, kmeans_purpose = "extrapolation", backwards_compatible_kmeans = backwards_compatible_kmeans)
  Kmeans[["cluster"]] = RANN::nn2(data = Kmeans[["centers"]], query = Return$Data_Extrap[, c("E_km", "N_km")], k = 1)$nn.idx[, 1]
  aggregate_vector = function(values_x, index_x, max_index, FUN = sum) {
    tapply(values_x, INDEX = factor(index_x, levels = 1:max_index), FUN = FUN)
  }
  a_el = matrix(NA, nrow = max_cells, ncol = ncol(Return$a_el))
  for (lI in 1:ncol(Return$a_el)) {
    a_el[, lI] = aggregate_vector(values_x = Return$a_el[, lI], index_x = Kmeans$cluster, max_index = max_cells)
  }
  Area_km2_x = aggregate_vector(values_x = Return$Area_km2_x, index_x = Kmeans$cluster, max_index = max_cells)
  Include = aggregate_vector(values_x = Return$Data_Extrap[, "Include"], index_x = Kmeans$cluster, max_index = max_cells, FUN = function(vec) { 
    any(vec > 0)
    })
  lonlat_g = project_coordinates(X = Kmeans$centers[, "E_km"], Y = Kmeans$centers[, "N_km"], projargs = "+proj=longlat +ellps=WGS84", origargs = Return$projargs)
  Data_Extrap = cbind(Lon = lonlat_g[, 1], Lat = lonlat_g[, 2], Include = Include, Kmeans$centers)
  Return = list(a_el = a_el, Data_Extrap = Data_Extrap, zone = Return$zone, projargs = Return$projargs, flip_around_dateline = Return$flip_around_dateline, Area_km2_x = Area_km2_x)
}
if (length(Region) > 1 & create_strata_per_region == TRUE) {
  Return$a_el = cbind(Total = rowSums(Return$a_el), Return$a_el)
}


# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "northwest_atlantic", 
                          purpose = "index2", 
                          strata.limits = strata_use, 
                          bias.correct = FALSE 

####
## What happens with our workflow
colSums(fit_f$extrapolation_list$a_el)

Region = settings$Region
projargs = NA
zone = settings$zone
strata.limits = settings$strata.limits
create_strata_per_region = FALSE
max_cells = settings$max_cells
input_grid = extrap_grid
observations_LL = NULL
grid_dim_km = c(2, 2)
maximum_distance_from_sample = NULL
grid_in_UTM = TRUE
grid_dim_LL = c(0.1, 0.1)
region = c("south_coast", "west_coast")
strata_to_use = c("SOG", "WCVI", "QCS", "HS", "WCHG")
epu_to_use = c("All", "Georges_Bank", "Mid_Atlantic_Bight", "Scotian_Shelf", "Gulf_of_Maine", "Other")[1]
survey = "Chatham_rise"
surveyname = "propInWCGBTS"
flip_around_dateline
nstart = 100
area_tolerance = 0.05
backwards_compatible_kmeans = FALSE
DirPath = paste0(getwd(), "/") 

rI = 1
Extrapolation_List = NULL
#strata_use<- list("Georges_Bank_A" = strata_use_orig[[1]][1:5])

if (tolower(Region[rI]) == "user") {
  if (is.null(input_grid)) {
    stop("Because you're using a user-supplied region, please provide 'input_grid' input")
  }
  if (!(all(c("Lat", "Lon", "Area_km2") %in% colnames(input_grid)))) {
    stop("'input_grid' must contain columns named 'Lat', 'Lon', and 'Area_km2'")
  }
  if (missing(flip_around_dateline)) 
    flip_around_dateline = FALSE
  Data_Extrap <- input_grid
  Area_km2_x = Data_Extrap[, "Area_km2"]
  Tmp = cbind(BEST_LAT_DD = Data_Extrap[, "Lat"], BEST_LON_DD = Data_Extrap[, "Lon"])
  if ("Depth" %in% colnames(Data_Extrap)) {
    Tmp = cbind(Tmp, BEST_DEPTH_M = Data_Extrap[, "Depth"])
  }
  a_el = as.data.frame(matrix(NA, nrow = nrow(Data_Extrap), ncol = nrow(strata.limits), dimnames = list(NULL, strata.limits[, "STRATA"])))
  for (l in 1:ncol(a_el)) {
    a_el[,l]<- match_strata_fn_aja(points = Tmp, strata_dataframe = strata.limits[l, , drop = FALSE], index_shapes = index_shapes)
    a_el[, l] = ifelse(is.na(a_el[, l]), 0, Area_km2_x)
  }
  tmpUTM = project_coordinates(X = Data_Extrap[, "Lon"], Y = Data_Extrap[, "Lat"], projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline)
  Data_Extrap = cbind(Data_Extrap, Include = 1)
  if (all(c("E_km", "N_km") %in% colnames(Data_Extrap))) {
    Data_Extrap[, c("E_km", "N_km")] = tmpUTM[, c("X", "Y")]
  } else {
    Data_Extrap = cbind(Data_Extrap, E_km = tmpUTM[, "X"], N_km = tmpUTM[, "Y"])
  }
  Extrapolation_List_Use = list(a_el = a_el, Data_Extrap = Data_Extrap, zone = attr(tmpUTM, "zone"), projargs = attr(tmpUTM, "projargs"), flip_around_dateline = flip_around_dateline, Area_km2_x = Area_km2_x)
}


# Make settings (turning off bias.correct to save time for example)
settings = make_settings( n_x = 250, 
                          Region = "User", 
                          purpose = "index2", 
                          strata.limits = strata_use, 
                          bias.correct = FALSE,
                          knot_method = "grid")

# Run model
fit_g = fit_model( settings = settings, 
                   input_grid = extrap_grid,
                   extrapolation_list = Extrapolation_List_Use,
                   Lat_i = dat_use2[,'LATITUDE'], 
                   Lon_i = dat_use2[,'LONGITUDE'], 
                   t_i = dat_use2[,'YEAR'], 
                   c_i = rep(0, nrow(dat_use2)), 
                   b_i = dat_use2[,'CATCH_WT_CAL'], 
                   a_i = rep(0.0112 * 1.852^2, nrow(dat_use2)),
                   v_i = rep("missing", nrow(dat_use2)),
                   run_model = TRUE)

plot(fit_g)


####
## Edited all functions...
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
colSums(fit_g$extrapolation_list$a_el)



























if (is.null(Extrapolation_List)) {
  if (is.null(observations_LL)) {
    stop("Because you're using a new Region[rI], please provide 'observations_LL' input with columns named `Lat` and `Lon`")
  }
  if (missing(flip_around_dateline)) 
    flip_around_dateline = FALSE
  Extrapolation_List = Prepare_Other_Extrapolation_Data_Fn(strata.limits = strata.limits, observations_LL = observations_LL, grid_dim_km = grid_dim_km, maximum_distance_from_sample = maximum_distance_from_sample, grid_in_UTM = grid_in_UTM, grid_dim_LL = grid_dim_LL, projargs = projargs, zone = zone, flip_around_dateline = flip_around_dateline)
}

