# Base image, only need tidyverse
FROM rocker/tidyverse:4.0.4

# create an R user
ENV USER rstudio

# cd but for Docker
WORKDIR /home/$USER/

# Additions to base image
COPY ./pipelines/samp_dat_prep/R/DockerConfig/add_requirements.R /home/$USER/add_requirements.R
RUN Rscript /home/$USER/add_requirements.R

# Copy all R functions for this pipeline step and drop in new directory
COPY ./pipelines/samp_dat_prep/R/ /home/$USER/R/
COPY ./R/ /home/$USER/R/

# Make them executable
RUN ls -lR /home/$USER/R/
RUN chmod +x  /home/$USER/R/*.R