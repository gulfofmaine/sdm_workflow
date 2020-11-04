# Penguin Flex
# 10/27/2020


####  Packages  ####
library(tidyverse)
library(gmRi)

####  Data  ####
sdm_path <- shared.path(group = "Mills Lab", folder = "Projects/sdm_workflow/data/")
peng <- read_csv(paste0(sdm_path, "penguins.csv"))


####  Cool  Analysis  ####
glimpse(peng)

peng %>% 
  pivot_longer(names_to = "metric", cols =  c(bill_length_mm, bill_depth_mm, flipper_length_mm, body_mass_g)) %>% 
  mutate(metric = str_replace_all(metric, "_", " "),
         metric = str_to_sentence(metric)) %>% 
  ggplot(aes(year, value, fill = species, color = species)) +
  geom_jitter(aes(shape = sex), size = .75, width = .1, height = 0) +
  facet_grid(metric ~ island, scales = "free") +
  scale_x_continuous(breaks = c(2007, 2008, 2009)) +
  geom_smooth(method = "lm",
              formula = y ~ x,
              se = FALSE) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "", y = "Measurement", caption = "Changes in body features within penguin pecies over time, for individual island populations.")

