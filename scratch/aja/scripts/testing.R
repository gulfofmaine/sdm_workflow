##### Andrew's code #####
library(tidyverse)
library(palmerpenguins)

# Read in the penguins data -- how would we handle file paths?
pen_data<- read_csv("~/Box/sdm_workflow/data/penguins.csv")

# Do something
pen_data %>% 
  count(species)
  
