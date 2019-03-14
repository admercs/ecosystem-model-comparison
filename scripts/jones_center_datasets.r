#=================================================================================================#
# jones_center_datasets.r
# Datasets and parsing functions for Jones Ecological Research Center Red Dirt flux tower (mesic)
# Note: Mesic (Red Dirt tower) Sites: 3, 3100, 3200, and 3300
#       Mesic (Red Dirt tower) Plots: 17, 21, 25, and 29
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

#-------------------------------------------------------------------------------------------------#
# General data:
#-------------------------------------------------------------------------------------------------#
# "JC002-01-tree-survey.csv"             species, DBH, 1995
# "JC002-02-banded-tree-survey.csv"      species, DBH, height, 1995
# "JC002-03-oak-age.csv"                 oak species, DBH, height, age, DBH increment, 1989-1994
# "JC003-01-tree-BA-biomass.csv"         longleaf pine DBH, growth, basal area, biomass, 1996-2006
# "JC003-02-dbh-01-14.csv"               species, DBH increment, 2001-2013
# "JC003-03-tree-productivity-09-14.csv" genus, DBH increment, 2009-2014
# "JC003-04-anpp-tree.csv"               species, DBH, basal area, biomass compartments, ANPP, 2004-2014
# "JC003-07-standing-litter-CN.csv"      litter chemistry, 2008-2011
# "JC004-01-root-CN.csv"                 root diameter, C:N, 2011-2013
# "JC005-01-tree-productivity.csv"       species, DBH, height, crown length, age, DBH increment, 1995-1999
# "JC005-02-dbh-88-99.csv"               genus, DBH, 1988-1999
# "JC005-03-oak-harvest.csv"             oak biomass compartments, 1995-1996
# "JC005-04-oak-harvest-branch.csv"      oak biomass compartments, 1995-1996
# "JC005-05-pine-harvest.csv"            pine biomass compartments, 1995-1996
# "JC005-06-pine-harvest-branch.csv"     pine biomass compartments, 1995-1996
# "JC005-07-anpp-tree.csv"               genus, DBH, basal area, biomass compartments, ANPP, 1989-1999
# "JC008-01-litterfall-wt.csv"           genus, litterfall mass, 2002-2014
# "JC008-02-litterfall-cn.csv"           genus, litterfall mass, C, N, 2002-2012
# "JC009-01-litterfall-weight.csv"       genus, litterfall mass, 1995-2000
# "JC009-03-litterfall-SP-chn.csv"       genus, litterfall C, N, 1995-1998
# "JC010-01-exchange.csv"                NEE, GEE, CO2, 2008-2013
# "JC010-02-met.csv"                     meteorology, 2008-2014
# "JC011-01-LI6400.csv"                  r_soil, 2009-2015 (manual)
# "JC011-03-LI8100.csv"                  r_soil, 2010-2015 (automated)
# "JC012-02-Monthly.csv"                 meteorology, 2002-2015
# "JC013-01-Daily.csv"                   GAEMN meteorology (tmax, tmin, precip), 2000-2015
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
# Metrics and units
#-------------------------------------------------------------------------------------------------#
# Note: kg m-2 (10 Mg ha-1); µmol m-2 sec-1 (1e-6 * 60 * 60 * 24 * 365.24 mol m-2 year-1)
# nee ................. kg C m-2 year-1
# biomass_ag .......... kg m-2
# c_ag ................ kg C m-2
# n_ag ................ kg N m-2
# biomass_bg .......... kg m-2
# c_bg ................ kg C m-2
# n_bg ................ kg N m-2
# c_so ................ kg C m-2
# n_so ................ kg N m-2
# r_soil .............. kg C m-2 year-1
# anpp ................ kg m-2 year-1
# biomass_sp .......... kg m-2
# abundance_sp ........ % of total n-trees

#-------------------------------------------------------------------------------------------------#
# Selected data:
#-------------------------------------------------------------------------------------------------#
# "JC004-01-root-CN.csv"
#   root/stem C:N ratio
# "JC008-02-litterfall-cn.csv"
#   litter/leaf C:N ratio
# "JC003-02-dbh-01-14.csv" or...
# "JC003-04-anpp-tree.csv"
#   timescale: 2004-2014
#   variables: biomass_ag, c_ag, n_ag, biomass_bg, c_ng, n_bg, anpp, species biomass
# "JC011-01-LI6400.csv"
#   timescale: 2009-2015
#   variables: r_soil
# "JC003-07-standing-litter-CN.csv"
#   timescale: 2008-2011
#   variables: c_so, n_so (estimated from litter C:N, biomass per Harvard Forest regression model)
# "JC010-01-exchange.csv"
#   timescale: 2008-2013
#   variables: nee
# "JC010-02-met.csv" (used as forcing data)
#   timescale: 2008-2014
#   variables: monthly tmean, tmax, pmean
#   Weather station: GAEMN Newton or NOAA Newton 8

#-------------------------------------------------------------------------------------------------#
# Time period:
#-------------------------------------------------------------------------------------------------#
# Intersection: 2009-2011 (3 years)
# Union:        2004-2015 (12 years)
# Climate:      2009-2013 (5 years)

#-------------------------------------------------------------------------------------------------#
# Base directory
#-------------------------------------------------------------------------------------------------#
base_jerc_data <- file.path(base, "observations", "jerc_rd")

#-------------------------------------------------------------------------------------------------#
# Species look-up table
#-------------------------------------------------------------------------------------------------#
# Using USDA PLANTS database codes: https://plants.usda.gov
species_lut <- data.frame(
  code  = c("PIPA","PRSE","QUFA","QUHE","QUIN","QULA","QULV","QUMA","QUNI","QUST","QUVI"),
  genus = c("Pinus","Prunus","Quercus","Quercus","Quercus","Quercus","Quercus","Quercus",
            "Quercus","Quercus","Quercus"),
  name  = c("Pinus palustris","Prunus serotina","Quercus falcata","Quercus hemispherica",
            "Quercus incana","Quercus laurifolia","Quercus laevis","Quercus margaretta",
            "Quercus nigra","Quercus stellata","Quercus virginiana")
)
#message(paste("Species LUT:", jsonlite::toJSON(species_lut, pretty=TRUE)))

#-------------------------------------------------------------------------------------------------#
# Observation file list
#-------------------------------------------------------------------------------------------------#
observation_list <- list(
  climate   = file.path(base_jerc_data, "JC010-02-met.csv"),
  fluxes    = file.path(base_jerc_data, "JC010-01-exchange.csv"),
  inventory = file.path(base_jerc_data, "JC003-04-anpp-tree.csv"),
  root      = file.path(base_jerc_data, "JC004-01-root-CN.csv"),
  litter    = file.path(base_jerc_data, "JC003-07-standing-litter-CN.csv"),
  r_soil    = file.path(base_jerc_data, "JC011-01-LI6400.csv")
)

#-------------------------------------------------------------------------------------------------#
# Functions
#-------------------------------------------------------------------------------------------------#

# Tree species aboveground biomass equations; kg; DBH cm; Chojnacky, Heath, Jenkins (2014)
get_tree_biomass_ag <- function(code, dbh) {
  switch(code,                                        #  WSG
         PIPA = { exp(-3.0506 + 2.6465 * log(dbh)) }, # 0.54
         PRSE = { exp(-2.2118 + 2.4133 * log(dbh)) }, # 0.41
         QUFA = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.52
         QUHE = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.59 oaks
         QUIN = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.59 oaks
         QULA = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.59 oaks
         QULV = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.59 oaks; species name?
         QUMA = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.59 oaks
         QUNI = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.56
         QUST = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.60
         QUVI = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.59
         NA                                           # else
  )
}

# Tree species belowground biomass equations; kg; DBH cm; Chojnacky, Heath, Jenkins (2014)
get_tree_biomass_bg <- function(code, dbh) {
  biomass_ag       = get_tree_biomass_ag(code, dbh)
  root_coarse_frac = exp(-1.4485 - 0.03476 * log(dbh))
  root_fine_frac   = exp(-1.8629 - 0.77534 * log(dbh))
  biomass_bg       = (biomass_ag * root_coarse_frac) + (biomass_ag * root_fine_frac)
  return(biomass_bg)
}

# Root biomass fraction by DBH class; treated as constants per Chojnacky, Heath, Jenkins (2014)
#plot(1:100, sapply(1:100, function(dbh) exp(-1.4485 - 0.03476 * log(dbh))), type="l")
#plot(1:100, sapply(1:100, function(dbh) exp(-1.8629 - 0.77534 * log(dbh))), type="l")

# Main

# Parent function containing child functions
get_observations <- function(metric="all") {
  funs = list(get_nee(), get_biomass_ag(), get_c_ag(), get_n_ag(), get_biomass_bg(),
              get_c_bg(), get_n_bg(), get_c_so(), get_n_so(), get_r_soil(), get_anpp(),
              get_biomass_sp(), get_abundance_sp())
  switch(metric,
    nee           = { get_nee() },
    biomass_ag    = { get_biomass_ag() },
    c_ag          = { get_c_ag() },
    n_ag          = { get_n_ag() },
    biomass_bg    = { get_biomass_bg() },
    c_bg          = { get_c_bg() },
    n_bg          = { get_n_bg() },
    c_so          = { get_c_so() },
    n_so          = { get_n_so() },
    r_soil        = { get_r_soil() },
    anpp          = { get_anpp() },
    biomass_sp    = { get_biomass_sp() },
    abundance_sp  = { get_abundance_sp() },
    biomass_total = { cbind(get_biomass_ag()[,1], get_biomass_ag()[,2] + get_biomass_bg()[,2]) },
    c_total       = { cbind(get_c_ag()[,1], get_c_ag()[,2] + get_c_bg()[,2]) },
    n_total       = { cbind(get_n_ag()[,1], get_n_ag()[,2] + get_n_bg()[,2]) },
    all           = { setNames(lapply(funs, eval),
                                   c("nee","biomass_ag","c_ag","n_ag","biomass_bg",
                                     "c_bg","n_bg","c_so","n_so","r_soil","anpp",
                                     "biomass_sp","abundance_sp")) },
    { stop(paste("Variable", metric, "not found")) }
  )
}

# Test
#test = get_observations("all")

# NEE; kg C m-2 year-1; Converted from µmol CO2 m-2 sec-1 (1e-6 moles)
get_nee <- function() {
  #df = read.csv(self$files$observations$flux, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_jerc_data, "JC010-01-exchange.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  df[df == -9999] = NA
  df$c_kg = df$nee * 1e-6 * 12 / 1000       # convert µmol CO2 to mol CO2 to g C to kg C
  df$c_kg = df$c_kg * 60 * 60 * 24 * 365.24 # convert kg C m-2 second-1 to kg C m-2 year-1
  return(setNames(aggregate(df$c_kg, by=list(df$year), FUN=mean), c("year","nee")))
}

# Aboveground biomass; kg m-2; 4*2500 m2 plots (10000 m2)
get_biomass_ag <- function() {
  #df = read.csv(self$files$observations$inventory, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_jerc_data, "JC003-04-anpp-tree.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  # ^ filters down species to PIPA QUNI QUVI QUIN
  colnames(df)[colnames(df) == "sb"] = "biomass_stem_branch"
  df$genus   = species_lut$genus[match(df$species, species_lut$code)]
  df$biomass_foliage = ifelse(df$genus == "Pinus", df$ba*0.0366, df$ba*0.0211) # model from Mitchell et al. (1999)
  df$biomass_ag = df$biomass_stem_branch + df$biomass_foliage
  df$biomass_ag = df$biomass_ag / (4*2500) # convert kg to kg m-2
  return(setNames(aggregate(df$biomass_ag, by=list(df$year), FUN=sum), c("year","biomass_ag")))
}

# Aboveground carbon; kg C m-2; 4x2500 m2 plots (10000 m2)
get_c_ag <- function() {
  df = get_biomass_ag()
  df$c_ag = df$biomass_ag * 0.5
  return(df[,c("year","c_ag")])
}

# Aboveground nitrogen; kg N m-2; 4x2500 m2 plots (10000 m2)
get_n_ag <- function() {
  #df = read.csv(self$files$observations$inventory, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_jerc_data, "JC003-04-anpp-tree.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  colnames(df)[colnames(df) == "sb"] = "biomass_stem_branch"
  #df = df[df$site == c(3,3100,3200,3300) & df$year %in% self$timescale, c("year","species","ba","sb")]
  # ^ filters down species to PIPA QUNI QUVI QUIN
  #stopifnot(all.equal(df$year, self$timescale))
  df$genus   = species_lut$genus[match(df$species, species_lut$code)]
  df$biomass_foliage = ifelse(df$genus == "Pinus", df$ba*0.0366, df$ba*0.0211) # model from Mitchell et al. (1999)
  df$cn_stem = ifelse(df$genus == "Pinus", 133.72073, 96.37033)
  df$cn_leaf = ifelse(df$genus == "Pinus", 255.10303, 85.25917)
  df$n_ag    = 0.5 * ((df$biomass_stem_branch / df$cn_stem) + (df$biomass_foliage / df$cn_leaf))
  df$n_ag    = df$n_ag / (4*2500) # convert kg to kg m-2
  return(setNames(aggregate(df$n_ag, by=list(df$year), FUN=sum), c("year","n_ag")))
}

# Belowground biomass; kg m-2; 4x2500 m2 plots (10000 m2)
# Uses general AGB-to-BGB model of Chojnacky, Heath, Jenkins (2014)
get_biomass_bg <- function() {
  df = get_biomass_ag()
  df$root_coarse = df$biomass_ag * 0.22
  df$root_fine   = df$biomass_ag * 0.01
  df$biomass_bg  = df$root_coarse + df$root_fine
  return(df[,c("year","biomass_bg")])
}

# Belowground C; kg C m-2; 4x2500 m2 plots (10000 m2)
get_c_bg <- function () {
  bgb = get_biomass_bg()
  #df = read.csv(self$files$observations$root, header=TRUE, stringsAsFactors=FALSE)
  #df  = read.csv(file.path(base_jerc_data, "JC004-01-root-CN.csv"), header=TRUE, stringsAsFactors=FALSE)
  #colnames(df) = tolower(colnames(df))
  #df = df[df$site == 3,]
  #bgb$c_bg = bgb$biomass_bg * mean(df$c, na.rm=TRUE)
  bgb$c_bg = bgb$biomass_bg * 0.5
  return(bgb[c("year","c_bg")])
}

# Belowground N; kg N m-2; 4x2500 m2 plots (10000 m2)
get_n_bg <- function () {
  bgb = get_biomass_bg()
  #df = read.csv(self$files$observations$root, header=TRUE, stringsAsFactors=FALSE)
  df  = read.csv(file.path(base_jerc_data, "JC004-01-root-CN.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  bgb$n_bg = bgb$biomass_bg * mean(df$n, na.rm=TRUE)
  return(bgb[c("year","n_bg")])
}

# Soil organic C; kg C m-2; Approximated with standing litter
get_c_so <- function() {
  #df = read.csv(self$files$observations$litter, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_jerc_data, "JC003-07-standing-litter-CN.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  df$kg_m2_c = df$mg_ha_c / 10
  return(setNames(aggregate(df$kg_m2_c, by=list(df$year), FUN=sum), c("year","c_so")))
}

# Soil organic N; kg N m-2; Approximated with standing litter
get_n_so <- function() {
  #df = read.csv(self$files$observations$litter, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_jerc_data, "JC003-07-standing-litter-CN.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  df$kg_m2_n = df$kg_ha_n / 10000
  return(setNames(aggregate(df$kg_m2_n, by=list(df$year), FUN=sum), c("year","n_so")))
}

# Soil respiration; µmol CO2 m-2 sec-1 (1e-6 moles); kg C m-2 year-1
get_r_soil <- function() {
  #df = read.csv(self$files$observations$r_soil, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_jerc_data, "JC011-01-LI6400.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  df$efflux = df$efflux * 1e-6 * 12 / 1000      # convert µmol CO2 to mol CO2 to g C to kg C
  df$efflux = df$efflux * 60 * 60 * 24 * 365.24 # convert kg C m-2 second-1 to kg C m-2 year-1
  return(setNames(aggregate(df$efflux, by=list(df$year), FUN=mean), c("year","r_soil")))
}

# Annual net primary productivity (ANPP); kg year-1; kg C m-2 year-1
get_anpp <- function() {
  #df = read.csv(self$files$observations$inventory, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_jerc_data, "JC003-04-anpp-tree.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  df$anpp = df$anpp / (4*2500) * 0.5 # kg mass to kg C m-2
  return(setNames(aggregate(df$anpp, by=list(df$year), FUN=sum), c("year","anpp")))
}

# Total species aboveground biomass; kg m-2; 4x2500 m2 plots (10000 m2)
get_biomass_sp <- function() {
  #df = read.csv(self$files$observations$inventory, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_jerc_data, "JC003-04-anpp-tree.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  # ^ filters down species to PIPA QUNI QUVI QUIN
  colnames(df)[colnames(df) == "sb"] = "biomass_stem_branch"
  df$genus   = species_lut$genus[match(df$species, species_lut$code)]
  df$biomass_foliage = ifelse(df$genus == "Pinus", df$ba*0.0366, df$ba*0.0211) # model from Mitchell et al. (1999)
  df$biomass_ag = df$biomass_stem_branch + df$biomass_foliage
  df$biomass_ag = df$biomass_ag / (4*2500)
  return(setNames(aggregate(df$biomass_ag, by=list(df$year, df$species), FUN=sum), c("year","species","biomass_ag")))
}

# Species relative abundance; percent of n-trees
get_abundance_sp <- function() {
  #df = read.csv(self$files$observations$inventory, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_jerc_data, "JC003-04-anpp-tree.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site == 3,]
  ntrees_yr      = setNames(aggregate(df$species, by=list(df$year), FUN=length), c("year","ntrees"))
  df$ntrees_yr   = ntrees_yr$ntrees[match(df$year, ntrees_yr$year)]
  df$fraction_yr = 1 / df$ntrees_yr # relative abundance of each tree
  return(setNames(aggregate(df$fraction_yr, by=list(df$year, df$species), FUN=sum), c("year","species","abundance")))
}
