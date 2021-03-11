#!/usr/bin/env Rscript
library(argparser, quietly=TRUE)
library(tidyverse, quietly=TRUE)
library(here, quietly = TRUE)

#####
## Workhorse function
#' @title Prep sample dataframe
#' 
#' @description This function takes in a cleaned csv file from a specific directory and subsets it to include only observations for species designated by species_vec and for observations made within a certain time span
#'
#' @param in_csv_dir Directory that has the cleaned sample csv file, where each row constitutes a unique tow-species occurrence observation
#' @param species_vec A vector of svspp species identifying code. If NULL, then all species are kept.
#' @param year_min An integer value specifying the minimum year to keep observations. If NULL, then all years are kept.
#' @param year_max An integer value specifying the maximum year to keep observations. If NULL, then all years are kept.
#' @param out_rds_dir Directory to save prepped sample rdata file.
#' 
#' @return sample data frame 
#' 
samp_dat_prep<- function(in_csv_dir, species_vec, year_min, year_max, out_rds_dir){

  # For debugging and walking through the function
  if(FALSE){
    in_csv_dir = "~/GitHub/sdm_workflow/scratch/aja/pachy/pfs/samp_dat_prep"
    species_vec = 
    year_min = 
    year_max = 
    out_rds_dir = "~/GitHub/sdm_workflow/scratch/aja/pachy/pfs/samp_dat_prep"
  }

  # Get full file path from the directory
  in_csv_path<- list.files(in_csv_dir, full.names = TRUE)

  # Read it in
  dat_temp<- read.csv(in_csv_path[[1]])

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
  saveRDS(dat_temp, paste0(out_rds_dir, "samp_dat_prepped.rds"))
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
  
  # Prep sample data
  samp_dat_prep(args$in_csv_dir, args$species_vec, args$year_min, args$year_max, args$out_rds_dir)
}