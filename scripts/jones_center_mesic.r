#=================================================================================================#
# jones_center_mesic.r
# Model comparison at Jones Ecological Research Center (mesic) Red Dirt tower site
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
# Model: PPA-SiBGC
#=================================================================================================#
ppa_jerc <- PPA_SiBGC$new(file.path(base, "models", "ppa_sibgc", "jerc_rd"))

# Check pre-defined parameters
ppa_jerc$name
ppa_jerc$version
ppa_jerc$directory
ppa_jerc$files
ppa_jerc$arguments
ppa_jerc$command

# Site name
ppa_jerc$site <- "jerc_rd"
ppa_jerc$site

# Required: Set a simulation timescale (years)
ppa_jerc$timescale <- 2009:2013
ppa_jerc$timescale

# Required: Set the plot size (m-2)
ppa_jerc$area <- 10000
ppa_jerc$area

# Update predictions if model was previously run; requires a timescale be set
ppa_jerc$update_predictions()
ppa_jerc$predictions

# Check the spatial bounds
ppa_jerc$bounds

# Modify arguments to command; requires update_command().
ppa_jerc$arguments$verbose <- ""
ppa_jerc$update_command()

# Optional: Backup original parameter files before running model
#model$archive_parameters(overwrite=FALSE)

# Set observation data list and methods
# Pre-processing to populate model$observations using get_observations("all")
source("R/jones_center_datasets.r")
ppa_jerc$set_observations(observations=observation_list, method=get_observations)
ppa_jerc$files$observations
ppa_jerc$observations

# Check get_observations methods; precomputed lists
names(ppa_jerc$observations)

#-------------------------------------------------------------------------------------------------#
# Run the model
out <- ppa_jerc$run()
out
#-------------------------------------------------------------------------------------------------#

# Model predictions
names(ppa_jerc$observations)
names(ppa_jerc$predictions)

# Check observations and predictions
ppa_jerc$observations
ppa_jerc$predictions

# Remove initial year ANPP
ppa_jerc$predictions$anpp = ppa_jerc$predictions$anpp[-1,]

# Plot predictions against observations; same variable name must be used!
# Legend: observations (blue) and predictions (red)
ppa_jerc$plot("nee")
ppa_jerc$plot("biomass_ag")
ppa_jerc$plot("c_ag")
ppa_jerc$plot("n_ag")
ppa_jerc$plot("biomass_bg")
ppa_jerc$plot("c_bg")
ppa_jerc$plot("n_bg")
ppa_jerc$plot("r_soil")
ppa_jerc$plot("anpp")
ppa_jerc$plot("c_so")
ppa_jerc$plot("n_so")
ppa_jerc$plot("biomass_sp")
ppa_jerc$plot("abundance_sp")

# SOC constant
ppa_jerc$observations$c_so$c_so
unique(ppa_jerc$predictions$c_so$c_so)

# SON constant
ppa_jerc$observations$n_so$n_so
unique(ppa_jerc$predictions$n_so$n_so)

# Save PDF plots to figures folder in base directory
ppa_jerc$save_pdf(file.path(base, "figures"), "nee")
ppa_jerc$save_pdf(file.path(base, "figures"), "biomass_ag")
ppa_jerc$save_pdf(file.path(base, "figures"), "c_ag")
ppa_jerc$save_pdf(file.path(base, "figures"), "n_ag")
ppa_jerc$save_pdf(file.path(base, "figures"), "biomass_bg")
ppa_jerc$save_pdf(file.path(base, "figures"), "c_bg")
ppa_jerc$save_pdf(file.path(base, "figures"), "n_bg")
ppa_jerc$save_pdf(file.path(base, "figures"), "r_soil")
ppa_jerc$save_pdf(file.path(base, "figures"), "anpp")
ppa_jerc$save_pdf(file.path(base, "figures"), "biomass_sp",   width=10, height=5)
ppa_jerc$save_pdf(file.path(base, "figures"), "abundance_sp", width=10, height=5)

# Plot observations
plot(ppa_jerc$observations$nee)
plot(ppa_jerc$observations$biomass_ag)
plot(ppa_jerc$observations$c_ag)
plot(ppa_jerc$observations$n_ag)
plot(ppa_jerc$observations$biomass_bg)
plot(ppa_jerc$observations$c_bg)
plot(ppa_jerc$observations$n_bg)
plot(ppa_jerc$observations$r_soil)
plot(ppa_jerc$observations$anpp)
plot(ppa_jerc$observations$biomass_sp)
plot(ppa_jerc$observations$abundance_sp)
ppa_jerc$observations$c_so
ppa_jerc$observations$n_so

# Plot predictions
plot(ppa_jerc$predictions$nee)
plot(ppa_jerc$predictions$biomass_ag)
plot(ppa_jerc$predictions$c_ag)
plot(ppa_jerc$predictions$n_ag)
plot(ppa_jerc$predictions$biomass_bg)
plot(ppa_jerc$predictions$c_bg)
plot(ppa_jerc$predictions$n_bg)
plot(ppa_jerc$predictions$r_soil)
plot(ppa_jerc$predictions$anpp)
plot(ppa_jerc$predictions$biomass_sp)
plot(ppa_jerc$predictions$abundance_sp)
plot(ppa_jerc$predictions$c_so)
plot(ppa_jerc$predictions$n_so)


#=================================================================================================#
# Model: LANDIS_2
#=================================================================================================#
landis_jerc <- LANDIS_2$new(file.path(base, "models", "landis_2", "jerc_rd"))

# Check pre-defined parameters
landis_jerc$name
landis_jerc$version
landis_jerc$directory
landis_jerc$files
landis_jerc$arguments
landis_jerc$command

# Site name
landis_jerc$site <- "jerc_rd"
landis_jerc$site

# Required: Set a simulation timescale (years)
landis_jerc$timescale <- 2009:2013
landis_jerc$timescale

# Required: Set the plot size (m-2)
landis_jerc$area <- 10000
landis_jerc$area

# Update predictions if model was previously run; requires a timescale be set
landis_jerc$update_predictions()
landis_jerc$predictions

# Check the spatial bounds
landis_jerc$bounds

# Modify arguments to command; requires update_command().
landis_jerc$arguments$scenario <- "scenario.txt"
landis_jerc$update_command()

# Optional: Backup original parameter files before running model
#model$archive_parameters(overwrite=FALSE)

# Set observation data list and methods
# Pre-processing to populate model$observations using get_observations("all")
source("R/jones_center_datasets.r")
landis_jerc$set_observations(observation_list, get_observations)
landis_jerc$files$observations
landis_jerc$observations

# Check get_observations methods; precomputed lists
names(landis_jerc$observations)

#-------------------------------------------------------------------------------------------------#
# Run the model
out <- landis_jerc$run()
out
#-------------------------------------------------------------------------------------------------#

# Model predictions
names(landis_jerc$observations)
names(landis_jerc$predictions)

# Check observations and predictions
landis_jerc$observations
landis_jerc$predictions

# Modify species names to match observation data a posteriori
lut <- data.frame(
  code = c("QUIN","PIPA","QULA"),
  name = c("BlueJackOak","LongPine","TurkeyOak" )
)
landis_jerc$predictions$biomass_sp$species   = lut$code[match(landis_jerc$predictions$biomass_sp$species,   lut$name)]
landis_jerc$predictions$abundance_sp$species = lut$code[match(landis_jerc$predictions$abundance_sp$species, lut$name)]

# Plot predictions against observations; same variable name must be used!
# Legend: observations (blue) and predictions (red)
landis_jerc$plot("nee")
landis_jerc$plot("biomass_ag")
landis_jerc$plot("c_ag")
landis_jerc$plot("n_ag")
landis_jerc$plot("biomass_bg")
landis_jerc$plot("c_bg")
landis_jerc$plot("n_bg")
landis_jerc$plot("r_soil")
landis_jerc$plot("anpp")
landis_jerc$plot("c_so")
landis_jerc$plot("n_so")
landis_jerc$plot("biomass_sp")
landis_jerc$plot("abundance_sp")

# SOC constant
landis_jerc$observations$c_so$c_so
unique(landis_jerc$predictions$c_so$c_so)

# SON constant
landis_jerc$observations$n_so$n_so
unique(landis_jerc$predictions$n_so$n_so)

# Save PDF plots to figures folder in base directory
landis_jerc$save_pdf(file.path(base, "figures"), "nee")
landis_jerc$save_pdf(file.path(base, "figures"), "biomass_ag")
landis_jerc$save_pdf(file.path(base, "figures"), "c_ag")
landis_jerc$save_pdf(file.path(base, "figures"), "n_ag")
landis_jerc$save_pdf(file.path(base, "figures"), "biomass_bg")
landis_jerc$save_pdf(file.path(base, "figures"), "c_bg")
landis_jerc$save_pdf(file.path(base, "figures"), "n_bg")
landis_jerc$save_pdf(file.path(base, "figures"), "r_soil")
landis_jerc$save_pdf(file.path(base, "figures"), "anpp")
landis_jerc$save_pdf(file.path(base, "figures"), "biomass_sp",   width=10, height=5)
landis_jerc$save_pdf(file.path(base, "figures"), "abundance_sp", width=10, height=5)

# Plot observations
plot(landis_jerc$observations$nee)
plot(landis_jerc$observations$biomass_ag)
plot(landis_jerc$observations$c_ag)
plot(landis_jerc$observations$n_ag)
plot(landis_jerc$observations$biomass_bg)
plot(landis_jerc$observations$c_bg)
plot(landis_jerc$observations$n_bg)
plot(landis_jerc$observations$r_soil)
plot(landis_jerc$observations$anpp)
plot(landis_jerc$observations$biomass_sp)
plot(landis_jerc$observations$abundance_sp)
landis_jerc$observations$c_so
landis_jerc$observations$n_so

# Plot predictions
plot(landis_jerc$predictions$nee)
plot(landis_jerc$predictions$biomass_ag)
plot(landis_jerc$predictions$c_ag)
plot(landis_jerc$predictions$n_ag)
plot(landis_jerc$predictions$biomass_bg)
plot(landis_jerc$predictions$c_bg)
plot(landis_jerc$predictions$n_bg)
plot(landis_jerc$predictions$r_soil)
plot(landis_jerc$predictions$anpp)
plot(landis_jerc$predictions$biomass_sp)
plot(landis_jerc$predictions$abundance_sp)
plot(landis_jerc$predictions$c_so)
plot(landis_jerc$predictions$n_so)


#=================================================================================================#
# Intercomparison plots and LaTeX table
#=================================================================================================#
ic_jerc <- Intercomparison$new(models=list(ppa_jerc, landis_jerc))

ic_jerc$names
ic_jerc$fitness

ic_jerc$plot("nee")
ic_jerc$plot("biomass_ag")
ic_jerc$plot("c_ag")
ic_jerc$plot("n_ag")
ic_jerc$plot("biomass_bg")
ic_jerc$plot("c_bg")
ic_jerc$plot("n_bg")
ic_jerc$plot("r_soil")
ic_jerc$plot("anpp")
ic_jerc$plot("biomass_sp")
ic_jerc$plot("abundance_sp")

ic_jerc$save_pdf(file.path(base, "figures"), "nee")
ic_jerc$save_pdf(file.path(base, "figures"), "biomass_ag")
ic_jerc$save_pdf(file.path(base, "figures"), "c_ag")
ic_jerc$save_pdf(file.path(base, "figures"), "n_ag")
ic_jerc$save_pdf(file.path(base, "figures"), "biomass_bg")
ic_jerc$save_pdf(file.path(base, "figures"), "c_bg")
ic_jerc$save_pdf(file.path(base, "figures"), "n_bg")
ic_jerc$save_pdf(file.path(base, "figures"), "r_soil")
ic_jerc$save_pdf(file.path(base, "figures"), "anpp")
ic_jerc$save_pdf(file.path(base, "figures"), "biomass_sp",   width=9, height=3)
ic_jerc$save_pdf(file.path(base, "figures"), "abundance_sp", width=9, height=3)

ic_jerc$save_tex(base)


#=================================================================================================#
# Combined results for both sites
#=================================================================================================#
hf_fitness   <- ic_hf$fitness
jerc_fitness <- ic_jerc$fitness

hf_fitness
jerc_fitness

hf_fitness[hf_fitness$metric=="mean",]
jerc_fitness[jerc_fitness$metric=="mean",]

colMeans(rbind(hf_fitness[hf_fitness$metric=="mean",], jerc_fitness[jerc_fitness$metric=="mean",])[,-1])


#=================================================================================================#
# Calculate fitness for both models and export to a LaTeX table
#=================================================================================================#
#jerc_results <- cbind(jerc_sortie_fitness, jerc_landis_fitness[,-1])
#means        <- as.numeric(colMeans(jerc_results[,-1], na.rm=TRUE))
#means        <- data.frame("mean", means[1], means[2], means[3], means[4], means[5], means[6],
#                           means[7], means[8])
#colnames(means) <- c("metric","r2","rmse","mae","me","r2","rmse","mae","me")
#jerc_results <- rbind(jerc_results, means)
#jerc_results
#stargazer::stargazer(jerc_results, type="latex", digits=2, summary=FALSE, rownames=FALSE, align=TRUE,
#                     header=FALSE, out=file.path(base, "jerc_fitness.tex"))
#
## Avergate fitness between the two sites
#hf_results[hf_results$metric=="mean",]
#jerc_results[jerc_results$metric=="mean",]
#
#colMeans(rbind(hf_results[hf_results$metric=="mean",], jerc_results[jerc_results$metric=="mean",])[,-1])

