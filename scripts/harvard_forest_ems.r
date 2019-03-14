#=================================================================================================#
# harvard_forest_ems.r
# Model comparison at Harvard Forest (Little Prospect Hill) EMS tower site
# Adam Erickson, PhD, Washington State University
# Contact: adam.michael.erickson@gmail.com
# March 14, 2019
# License: Apache 2.0
#
# Copyright 2019 Washington State University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#=================================================================================================#

source("R/classes.r")
source("R/text.r")
source("R/utilities.r")

# Set base directory
base <- file.path()

#=================================================================================================#
# Metrics
#=================================================================================================#
metric_df <- data.frame(
  metric      = c("nee","biomass_ag","c_ag","n_ag","biomass_bg","c_bg","n_bg","c_so","n_so","r_soil",
                  "anpp","biomass_sp","abundance_sp"),
  units       = c("kg C m^2 year^-1","kg mass m^2","kg C m^2","kg N m^2","kg mass m^2","kg C m^2",
                  "kg N m^2","kg C m^2","kg N m^2","kg C m^2 year^-1","kg mass m^2 year^-1",
                  "kg mass m^2","percent"),
  description = c("Net ecosystem exchange","Aboveground biomass","Aboveground C","Aboveground N",
                  "Belowground biomass","Belowground C","Belowground N","Soil organic C",
                  "Soil organic N","Respiration","Aboveground net primary production",
                  "Species total biomass","Species relative abundance")
)
metric_df
stargazer::stargazer(metric_df, type="latex", summary=FALSE, rownames=FALSE,
                     out=file.path(base, "observations", "metrics.tex"))

#=================================================================================================#
# Model: Sortie-PPA
#=================================================================================================#
ppa_hf <- PPA_SiBGC$new(file.path(base, "models", "ppa_sibgc", "hf_ems"))

# Check pre-defined parameters
ppa_hf$name
ppa_hf$version
ppa_hf$directory
ppa_hf$files
ppa_hf$arguments
ppa_hf$command

# Site name
ppa_hf$site <- "hf_ems"
ppa_hf$site

# Required: Set a simulation timescale (years)
ppa_hf$timescale <- 2002:2012
ppa_hf$timescale

# Required: Set the plot size (m-2)
ppa_hf$area <- 10681.42
ppa_hf$area

# Update predictions if model was previously run; requires a timescale be set
ppa_hf$update_predictions()
ppa_hf$predictions

# Check the spatial bounds
ppa_hf$bounds

# Modify arguments to command; requires update_command().
ppa_hf$arguments$verbose <- ""
ppa_hf$update_command()

# Optional: Backup original parameter files before running model
#model$archive_parameters(overwrite=FALSE)

# Set observation data list and parsing methods
# Pre-processing to populate model$observations using get_observations("all")
source("R/harvard_forest_datasets.r")
ppa_hf$set_observations(observations=observation_list, method=get_observations)
ppa_hf$files$observations
ppa_hf$observations

# Check get_observations methods; precomputed lists
names(ppa_hf$observations)

#-------------------------------------------------------------------------------------------------#
# Run the model
out <- ppa_hf$run()
out
#-------------------------------------------------------------------------------------------------#

# Model predictions
names(ppa_hf$observations)
names(ppa_hf$predictions)

# Check observations and predictions
ppa_hf$observations
ppa_hf$predictions

# Remove initial year ANPP
ppa_hf$predictions$anpp <- ppa_hf$predictions$anpp[-1,]

# Plot predictions against observations; same variable name must be used!
# Legend: observations (blue) and predictions (red)
ppa_hf$plot("nee")
ppa_hf$plot("biomass_ag")
ppa_hf$plot("c_ag")
ppa_hf$plot("n_ag")
ppa_hf$plot("biomass_bg")
ppa_hf$plot("c_bg")
ppa_hf$plot("n_bg")
ppa_hf$plot("r_soil")
ppa_hf$plot("anpp")
ppa_hf$plot("c_so")
ppa_hf$plot("n_so")
ppa_hf$plot("biomass_sp")
ppa_hf$plot("abundance_sp")

# SOC constant
ppa_hf$observations$c_so$c_so
ppa_hf$predictions$c_so$c_so

# SON constant
ppa_hf$observations$n_so$n_so
ppa_hf$predictions$n_so$n_so

# Save PDF plots to figures folder in base directory
ppa_hf$save_pdf(file.path(base, "figures"), "nee")
ppa_hf$save_pdf(file.path(base, "figures"), "biomass_ag")
ppa_hf$save_pdf(file.path(base, "figures"), "c_ag")
ppa_hf$save_pdf(file.path(base, "figures"), "n_ag")
ppa_hf$save_pdf(file.path(base, "figures"), "biomass_bg")
ppa_hf$save_pdf(file.path(base, "figures"), "c_bg")
ppa_hf$save_pdf(file.path(base, "figures"), "n_bg")
ppa_hf$save_pdf(file.path(base, "figures"), "r_soil")
ppa_hf$save_pdf(file.path(base, "figures"), "anpp")
ppa_hf$save_pdf(file.path(base, "figures"), "biomass_sp",   width=10, height=5)
ppa_hf$save_pdf(file.path(base, "figures"), "abundance_sp", width=10, height=5)

# Plot observations
plot(ppa_hf$observations$nee)
plot(ppa_hf$observations$biomass_ag)
plot(ppa_hf$observations$c_ag)
plot(ppa_hf$observations$n_ag)
plot(ppa_hf$observations$biomass_bg)
plot(ppa_hf$observations$c_bg)
plot(ppa_hf$observations$n_bg)
plot(ppa_hf$observations$r_soil)
plot(ppa_hf$observations$anpp)
plot(ppa_hf$observations$biomass_sp)
plot(ppa_hf$observations$abundance_sp)
ppa_hf$observations$c_so
ppa_hf$observations$n_so

# Plot predictions
plot(ppa_hf$predictions$nee)
plot(ppa_hf$predictions$biomass_ag)
plot(ppa_hf$predictions$c_ag)
plot(ppa_hf$predictions$n_ag)
plot(ppa_hf$predictions$biomass_bg)
plot(ppa_hf$predictions$c_bg)
plot(ppa_hf$predictions$n_bg)
plot(ppa_hf$predictions$r_soil)
plot(ppa_hf$predictions$anpp)
plot(ppa_hf$predictions$biomass_sp)
plot(ppa_hf$predictions$abundance_sp)
plot(ppa_hf$predictions$c_so)
plot(ppa_hf$predictions$n_so)


#=================================================================================================#
# Model: LANDIS-2
#=================================================================================================#
landis_hf <- LANDIS_2$new(file.path(base, "models", "landis_2", "hf_ems"))

# Check pre-defined parameters
landis_hf$name
landis_hf$version
landis_hf$directory
landis_hf$files
landis_hf$arguments
landis_hf$command

# Site name
landis_hf$site <- "hf_ems"
landis_hf$site

# Required: Set a simulation timescale (years)
landis_hf$timescale <- 2002:2012
landis_hf$timescale

# Required: Set the plot size (m-2)
landis_hf$area <- 10681.42
landis_hf$area

# Update predictions if model was previously run; requires a timescale be set
landis_hf$update_predictions()
landis_hf$predictions

# Check the spatial bounds
landis_hf$bounds

# Modify arguments to command; requires update_command().
landis_hf$arguments$scenario <- "scenario.txt"
landis_hf$update_command()

# Optional: Backup original parameter files before running model
#model$archive_parameters(overwrite=FALSE)

# Set observation data list and parsing methods
# Pre-processing to populate model$observations using get_observations("all")
source("R/harvard_forest_datasets.r")
landis_hf$set_observations(observations=observation_list, method=get_observations)
landis_hf$files$observations
landis_hf$observations

# Check get_observations methods; precomputed lists
names(landis_hf$observations)

#-------------------------------------------------------------------------------------------------#
# Run the model
out <- landis_hf$run()
out
#-------------------------------------------------------------------------------------------------#

# Model predictions
names(landis_hf$observations)
names(landis_hf$predictions)

# Check observations and predictions
landis_hf$observations
landis_hf$predictions

# Modify species names to match observation data a posteriori
lut <- data.frame(
  code = c("ACRU","QURU"),
  name = c("RedMaple","RedOak")
)
landis_hf$predictions$biomass_sp$species   <- lut$code[match(landis_hf$predictions$biomass_sp$species,   lut$name)]
landis_hf$predictions$abundance_sp$species <- lut$code[match(landis_hf$predictions$abundance_sp$species, lut$name)]

# Plot predictions against observations; same variable name must be used!
# Legend: observations (blue) and predictions (red)
landis_hf$plot("nee")
landis_hf$plot("biomass_ag")
landis_hf$plot("c_ag")
landis_hf$plot("n_ag")
landis_hf$plot("biomass_bg")
landis_hf$plot("c_bg")
landis_hf$plot("n_bg")
landis_hf$plot("r_soil")
landis_hf$plot("anpp")
landis_hf$plot("c_so")
landis_hf$plot("n_so")
landis_hf$plot("biomass_sp")
landis_hf$plot("abundance_sp")

# SOC constant
landis_hf$observations$c_so$c_so
landis_hf$predictions$c_so$c_so

# SON constant
landis_hf$observations$n_so$n_so
landis_hf$predictions$n_so$n_so

# Save PDF plots to figures folder in base directory
landis_hf$save_pdf(file.path(base, "figures"), "nee")
landis_hf$save_pdf(file.path(base, "figures"), "biomass_ag")
landis_hf$save_pdf(file.path(base, "figures"), "c_ag")
landis_hf$save_pdf(file.path(base, "figures"), "n_ag")
landis_hf$save_pdf(file.path(base, "figures"), "biomass_bg")
landis_hf$save_pdf(file.path(base, "figures"), "c_bg")
landis_hf$save_pdf(file.path(base, "figures"), "n_bg")
landis_hf$save_pdf(file.path(base, "figures"), "r_soil")
landis_hf$save_pdf(file.path(base, "figures"), "anpp")
landis_hf$save_pdf(file.path(base, "figures"), "biomass_sp",   width=10, height=5)
landis_hf$save_pdf(file.path(base, "figures"), "abundance_sp", width=10, height=5)

# Plot observations
plot(landis_hf$observations$nee)
plot(landis_hf$observations$biomass_ag)
plot(landis_hf$observations$c_ag)
plot(landis_hf$observations$n_ag)
plot(landis_hf$observations$biomass_bg)
plot(landis_hf$observations$c_bg)
plot(landis_hf$observations$n_bg)
plot(landis_hf$observations$r_soil)
plot(landis_hf$observations$anpp)
plot(landis_hf$observations$biomass_sp)
plot(landis_hf$observations$abundance_sp)
landis_hf$observations$c_so
landis_hf$observations$n_so

# Plot predictions
plot(landis_hf$predictions$nee)
plot(landis_hf$predictions$biomass_ag)
plot(landis_hf$predictions$c_ag)
plot(landis_hf$predictions$n_ag)
plot(landis_hf$predictions$biomass_bg)
plot(landis_hf$predictions$c_bg)
plot(landis_hf$predictions$n_bg)
plot(landis_hf$predictions$r_soil)
plot(landis_hf$predictions$anpp)
plot(landis_hf$predictions$biomass_sp)
plot(landis_hf$predictions$abundance_sp)
plot(landis_hf$predictions$c_so)
plot(landis_hf$predictions$n_so)


#=================================================================================================#
# Intercomparison plots and LaTeX table
#=================================================================================================#
ic_hf <- Intercomparison$new(models=list(ppa_hf, landis_hf))

ic_hf$names
ic_hf$fitness

ic_hf$plot("nee")
ic_hf$plot("biomass_ag")
ic_hf$plot("c_ag")
ic_hf$plot("n_ag")
ic_hf$plot("biomass_bg")
ic_hf$plot("c_bg")
ic_hf$plot("n_bg")
ic_hf$plot("r_soil")
ic_hf$plot("anpp")
ic_hf$plot("biomass_sp")
ic_hf$plot("abundance_sp")

ic_hf$save_pdf(file.path(base, "figures"), "nee")
ic_hf$save_pdf(file.path(base, "figures"), "biomass_ag")
ic_hf$save_pdf(file.path(base, "figures"), "c_ag")
ic_hf$save_pdf(file.path(base, "figures"), "n_ag")
ic_hf$save_pdf(file.path(base, "figures"), "biomass_bg")
ic_hf$save_pdf(file.path(base, "figures"), "c_bg")
ic_hf$save_pdf(file.path(base, "figures"), "n_bg")
ic_hf$save_pdf(file.path(base, "figures"), "r_soil")
ic_hf$save_pdf(file.path(base, "figures"), "anpp")
ic_hf$save_pdf(file.path(base, "figures"), "biomass_sp",   width=9, height=3)
ic_hf$save_pdf(file.path(base, "figures"), "abundance_sp", width=9, height=3)

ic_hf$save_tex(base)
