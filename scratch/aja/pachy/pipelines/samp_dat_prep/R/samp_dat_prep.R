#!/usr/bin/env Rscript
library(argparser, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(here, quietly = TRUE)

#####
## Workhorse function
samp_dat_prep<- function(in_csv_dir, species_vec, year_min, year_max){

  # Get full file path from the directory
  in_csv_path<- list.files(in_csv_dir, full.names = TRUE)

  # Read it in
  dat_temp<- read.csv(in_csv_path)

  # Filter if necessary, species first
  if(!is.null(species_vec)){
    dat_temp<- dat_temp %>%
      filter(., svspp %in% species_vec)
  } 
  
  # Now years
  if(!is.null(year_min)){
    dat_temp<- dat_temp %>%
      filter(., year_orig >= year_min & year_orig <= year_max)
  }

  # Return and save
  return(dat_temp)
}

#####
## Wrapper function for writing out
samp_dat_prep_out_rds<- function(in_csv_dir, species_vec, year_min, year_max, out_rds_dir){
  
  # Run samp_dat_prep
  samp_dat_prep_out<- samp_dat_prep(in_csv_dir, species_vec, year_min, year_max)
  
  # Save it
  saveRDS(samp_dat_prep_out, file = paste0(out_rds_dir, "samp_dat_prep_new3.rds"))
  
}

if(!interactive()){
  # Build up our command line argument parser
  parser<- arg_parser("Prep sample data")
  # Our first argument, the directory where our input file resides
  parser<- add_argument(parser, "in_csv_dir", help = "Input .csv file directory")
  # Our second argument, the species to keep
  parser<- add_argument(parser, "species_vec", type = "numeric", help = "The species svspp code", default = NULL)
  # Our third argument, the min year to keep
  parser<- add_argument(parser, "year_min", type = "numeric", help = "The minimum year to keep", default = NULL)
  # Our fourth argument, the max year to keep
  parser<- add_argument(parser, "year_max", type = "numeric", help = "The minimum year to keep", default = NULL)
  # Our fifth argument, the directory to save the prepped sample Rdata file
  parser<- add_argument(parser, "out_rds_dir", help = "Output Rdata directory")
  
  # Parse the arguments
  args<- parse_args(parser)
  print(args)
  
  print("You got the change incorporated")
  print("You did it again, nice job!")
  print("You did it again and again, finally!!!")


  # Prep sample data
  samp_dat_prep_out_rds(args$in_csv_dir, args$species_vec, args$year_min, args$year_max, args$out_rds_dir)
}