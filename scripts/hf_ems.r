#=================================================================================================#
# hf_ems.r
# Model comparison at Harvard Forest (Little Prospect Hill) EMS tower site
# Adam Erickson, Washington State University
# June 29, 2018
#
# License:
# Copyright 2018 Washington State University
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
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
base <- file.path("C:/Users/bob/forestmodels")

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
model <- SortiePPA$new(file.path(base, "models", "sortie_ppa", "hf_ems"))

# Check pre-defined parameters
model$name
model$version
model$directory
model$files
model$arguments
model$command

# Site name
model$site <- "hf_ems"
model$site

# Required: Set a simulation timescale (years)
model$timescale <- 2002:2012
model$timescale

# Required: Set the plot size (m-2)
model$area <- 10681.42
model$area

# Update predictions if model was previously run; requires a timescale be set
model$update_predictions()
model$predictions

# Check the spatial bounds
model$bounds

# Modify arguments to command; requires update_command().
model$arguments$verbose <- ""
model$update_command()

# Optional: Backup original parameter files before running model
#model$archive_parameters(overwrite=FALSE)

# Set observation data list and parsing methods
# Pre-processing to populate model$observations using get_observations("all")
source("hf_ems_data.r")
model$set_observations(observations=observation_list, method=get_observations)
model$files$observations
model$observations

# Check get_observations methods; precomputed lists
names(model$observations)

#-------------------------------------------------------------------------------------------------#
# Run the model
out <- model$run()
out
#-------------------------------------------------------------------------------------------------#

# Model predictions
names(model$observations)
names(model$predictions)

# Check observations and predictions
model$observations
model$predictions

# Generate model fitness table
model$get_fitness()
hf_sortie_fitness <- model$fitness

# Plot predictions against observations; same variable name must be used!
# Legend: observations (blue) and predictions (red)
model$plot("nee")
model$plot("biomass_ag")
model$plot("c_ag")
model$plot("n_ag")
model$plot("biomass_bg")
model$plot("c_bg")
model$plot("n_bg")
model$plot("r_soil")
model$plot("anpp")
model$plot("c_so")
model$plot("n_so")
model$plot("biomass_sp")
model$plot("abundance_sp")

# SOC constant
model$observations$c_so$c_so
model$predictions$c_so$c_so

# SON constant
model$observations$n_so$n_so
model$predictions$n_so$n_so

# Save PDF plots to figures folder in base directory
model$save_pdf(file.path(base, "figures"), "nee")
model$save_pdf(file.path(base, "figures"), "biomass_ag")
model$save_pdf(file.path(base, "figures"), "c_ag")
model$save_pdf(file.path(base, "figures"), "n_ag")
model$save_pdf(file.path(base, "figures"), "biomass_bg")
model$save_pdf(file.path(base, "figures"), "c_bg")
model$save_pdf(file.path(base, "figures"), "n_bg")
model$save_pdf(file.path(base, "figures"), "r_soil")
model$save_pdf(file.path(base, "figures"), "anpp")
model$save_pdf(file.path(base, "figures"), "biomass_sp",   width=10, height=5)
model$save_pdf(file.path(base, "figures"), "abundance_sp", width=10, height=5)

# Plot observations
plot(model$observations$nee)
plot(model$observations$biomass_ag)
plot(model$observations$c_ag)
plot(model$observations$n_ag)
plot(model$observations$biomass_bg)
plot(model$observations$c_bg)
plot(model$observations$n_bg)
plot(model$observations$r_soil)
plot(model$observations$anpp)
plot(model$observations$biomass_sp)
plot(model$observations$abundance_sp)
model$observations$c_so
model$observations$n_so

# Plot predictions
plot(model$predictions$nee)
plot(model$predictions$biomass_ag)
plot(model$predictions$c_ag)
plot(model$predictions$n_ag)
plot(model$predictions$biomass_bg)
plot(model$predictions$c_bg)
plot(model$predictions$n_bg)
plot(model$predictions$r_soil)
plot(model$predictions$anpp)
plot(model$predictions$biomass_sp)
plot(model$predictions$abundance_sp)
plot(model$predictions$c_so)
plot(model$predictions$n_so)


#=================================================================================================#
# Model: Landis-2
#=================================================================================================#
model <- Landis2$new(file.path(base, "models", "landis_2", "hf_ems"))

# Check pre-defined parameters
model$name
model$version
model$directory
model$files
model$arguments
model$command

# Site name
model$site <- "hf_ems"
model$site

# Required: Set a simulation timescale (years)
model$timescale <- 2002:2012
model$timescale

# Required: Set the plot size (m-2)
model$area <- 10681.42
model$area

# Update predictions if model was previously run; requires a timescale be set
model$update_predictions()
model$predictions

# Check the spatial bounds
model$bounds

# Modify arguments to command; requires update_command().
model$arguments$scenario <- "scenario.txt"
model$update_command()

# Optional: Backup original parameter files before running model
#model$archive_parameters(overwrite=FALSE)

# Set observation data list and parsing methods
# Pre-processing to populate model$observations using get_observations("all")
source("hf_ems_data.r")
model$set_observations(observations=observation_list, method=get_observations)
model$files$observations
model$observations

# Check get_observations methods; precomputed lists
names(model$observations)

#-------------------------------------------------------------------------------------------------#
# Run the model
out <- model$run()
out
#-------------------------------------------------------------------------------------------------#

# Model predictions
names(model$observations)
names(model$predictions)

# Check observations and predictions
model$observations
model$predictions

# Modify species names to match observation data a posteriori
lut <- data.frame(
  code = c("ACRU","QURU"),
  name = c("RedMaple","RedOak")
)
model$predictions$biomass_sp$species   = lut$code[match(model$predictions$biomass_sp$species,   lut$name)]
model$predictions$abundance_sp$species = lut$code[match(model$predictions$abundance_sp$species, lut$name)]

# Generate model fitness table
model$get_fitness()
hf_landis_fitness <- model$fitness

# Plot predictions against observations; same variable name must be used!
# Legend: observations (blue) and predictions (red)
model$plot("nee")
model$plot("biomass_ag")
model$plot("c_ag")
model$plot("n_ag")
model$plot("biomass_bg")
model$plot("c_bg")
model$plot("n_bg")
model$plot("r_soil")
model$plot("anpp")
model$plot("c_so")
model$plot("n_so")
model$plot("biomass_sp")
model$plot("abundance_sp")

# SOC constant
model$observations$c_so$c_so
model$predictions$c_so$c_so

# SON constant
model$observations$n_so$n_so
model$predictions$n_so$n_so

# Save PDF plots to figures folder in base directory
model$save_pdf(file.path(base, "figures"), "nee")
model$save_pdf(file.path(base, "figures"), "biomass_ag")
model$save_pdf(file.path(base, "figures"), "c_ag")
model$save_pdf(file.path(base, "figures"), "n_ag")
model$save_pdf(file.path(base, "figures"), "biomass_bg")
model$save_pdf(file.path(base, "figures"), "c_bg")
model$save_pdf(file.path(base, "figures"), "n_bg")
model$save_pdf(file.path(base, "figures"), "r_soil")
model$save_pdf(file.path(base, "figures"), "anpp")
model$save_pdf(file.path(base, "figures"), "biomass_sp",   width=10, height=5)
model$save_pdf(file.path(base, "figures"), "abundance_sp", width=10, height=5)

# Plot observations
plot(model$observations$nee)
plot(model$observations$biomass_ag)
plot(model$observations$c_ag)
plot(model$observations$n_ag)
plot(model$observations$biomass_bg)
plot(model$observations$c_bg)
plot(model$observations$n_bg)
plot(model$observations$r_soil)
plot(model$observations$anpp)
plot(model$observations$biomass_sp)
plot(model$observations$abundance_sp)
model$observations$c_so
model$observations$n_so

# Plot predictions
plot(model$predictions$nee)
plot(model$predictions$biomass_ag)
plot(model$predictions$c_ag)
plot(model$predictions$n_ag)
plot(model$predictions$biomass_bg)
plot(model$predictions$c_bg)
plot(model$predictions$n_bg)
plot(model$predictions$r_soil)
plot(model$predictions$anpp)
plot(model$predictions$biomass_sp)
plot(model$predictions$abundance_sp)
plot(model$predictions$c_so)
plot(model$predictions$n_so)


#=================================================================================================#
# Calculate fitness for both models and export to a LaTeX table
#=================================================================================================#
hf_results <- cbind(hf_sortie_fitness, hf_landis_fitness[,-1])
means      <- as.numeric(colMeans(hf_results[,-1], na.rm=TRUE))
means      <- data.frame("mean", means[1], means[2], means[3], means[4], means[5], means[6])
colnames(means) <- c("metric","r2","rmse","mae","r2","rmse","mae")
hf_results <- rbind(hf_results, means)
hf_results
stargazer::stargazer(hf_results, type="latex", digits=2, summary=FALSE, rownames=FALSE, align=TRUE,
                     header=FALSE, out=file.path(base, "hf_fitness.tex"))


