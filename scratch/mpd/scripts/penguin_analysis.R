library(tidyverse)

####  Data  ####
pen_data <- read_csv("~/Box/Mills Lab/Projects/sdm_workflow/data/penguins.csv")

glimpse(pen_data)

pen_data %>% 
  group_by(species, sex) %>% 
  summarise(bill_ln = mean(bill_length_mm)) %>% 
  ggplot() + geom_col(aes(species, bill_ln, fill = sex), position = "dodge")
