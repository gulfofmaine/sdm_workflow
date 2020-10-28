##### Andrew's code #####
library(tidyverse)
library(palmerpenguins)
library(argparser)

parser<- arg_parser("Explore penguin data")
parser<- add_argument(parser, "penguin_path", help="Path to CSV containing penguin data")
args<- parse_args(parser)

# Read in the penguins data -- how would we handle file paths?
pen_data<- read_csv(args$penguin_path)

# Do something
pen_data %>% 
  count(species)
  
