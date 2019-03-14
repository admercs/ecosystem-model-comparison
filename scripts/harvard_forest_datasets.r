#=================================================================================================#
# harvard_forest_datasets.r
# Datasets and parsing functions for Harvard Forest EMS flux tower (Little Prospect Hill)
# Note: site class is ems (EMS tower)
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
# Datasets specific to EMS tower:
#-------------------------------------------------------------------------------------------------#
# HF004: Canopy-Atmosphere Exchange of Carbon, Water and Energy at Harvard Forest EMS Tower since 1991
# HF006: Soil Respiration, Temperature and Moisture at Harvard Forest EMS Tower since 1995
# HF060: EMS - Automated High-Frequency Methane Data
# HF066: Concentrations and Surface Exchange of Air Pollutants at Harvard Forest EMS Tower since 1990
# HF067: CFCs and Radiatively Important Trace Species at Harvard Forest EMS Tower 1996-2005
# HF068: Soil Respiration Along a Hydrological Gradient at Harvard Forest EMS Tower 2003-2006
# HF069: Biomass Inventories at Harvard Forest EMS Tower since 1993
# HF102: Radiation Measurements at Harvard Forest EMS Tower 1991-2007
# HF145: Hydrocarbon Concentrations at Harvard Forest EMS Tower 1992-2001
# HF206: Microclimate at Harvard Forest HEM, LPH and EMS Towers since 2005
# HF209: Isotopic Composition of Net Ecosystem CO2 Exchange at Harvard Forest EMS Tower since 2011
# HF214: Fluxes of Carbonyl Sulfide at Harvard Forest EMS Tower since 2010
# HF288: Fluxes of Molecular Hydrogen (H2) at Harvard Forest EMS Tower 2010-2012
# HF306: Atmospheric Oxygen and Carbon Dioxide at Harvard Forest EMS Tower since 2006

#-------------------------------------------------------------------------------------------------#
# General data:
#-------------------------------------------------------------------------------------------------#
# "hf069-13-annual-summ.csv"          biomass_ag (woody), ANPP (woody), recruitment, mortality, 1993-2017
# "hf069-15-annual-summ.csv"          biomass_ag (woody), ANPP (woody), recruitment, mortality, 1994-2013 (old)
# "hf069-14-understory-summ.csv"      understory regenertion, biomass, 2004-2006
# "hf001-04-monthly-m.csv"            meteorology, 2001-2014
# "hf004-02-filled.csv"               NEE/GEE, respiration (ecosystem), 1991-2015
# "hf069-09-trees.csv"                species, DBH, 1993-2017
# "hf069-11-trees.csv"                species, DBH, 1998-2013
# "hf008-04-tree.csv"                 species, dbh, 1998-2008
# "hf213-01-hf-inventory.csv"         species, dbh, coordinates, 2009-2011
# "hf154-03-lph-cores.csv"            ANPP, 2005
# "hf133-03-productivity.csv"         biomass, 1996-1997
# "hf137-01-tree.csv"                 biomass, ANPP, species, 1996 & 2006
# "hf008-01-foliage.csv"              canopy chemistry, 1988-2011
# "hf008-02-litter.csv"               canopy chemistry, 1988-2002
# "hf062-01-plot.csv"                 canopy chemistry, 1992
# "hf062-02-tree.csv"                 canopy chemistry, 1992
# "hf069-06-canopy-chem.csv"          canopy chemistry, 1998-2008
# "hf008-07-decomp.csv"               litter chemistry, 1998-2008
# "hf069-05-litter-chem.csv"          litter chemistry, 1998-2000
# "hf047-01-regen.csv"                regeneration, 1990-2000
# "hf005-05-soil-respiration.csv"     r_soil, 1991-2002, C = control
# "hf006-01-soil-respiration.csv"     r_soil, 1995-2010
# "hf006-03-autochamber.csv"          r_soil, 2003
# "hf068-01-soil-resp.csv"            r_soil, 2003-2006
# "hf045-01-soil-resp.csv"            r_soil, 2006-2012
# "hf148-02-lph-soil-respiration.csv" r_soil, 2002-2007
# "hf194-01-soil-resp.csv"            r_soil, 1993-2006
# "hf194-02-synthesis-data"           r_soil, 1991-2008
# "hf243-01-soil-respiration.csv"     r_soil (Finzi), 2009-2010
# "hf040-02-species.csv"              species abundance, 1937-1991
# "hf015-06-species.csv"              species codes
# "hf040-01-species.csv"              species codes
# "hf222-01-vegetation.csv"           species codes
# "hf039-02-tree.csv"                 species DBH-height, 1937
# "hf005-01-trace-gas.csv"            co2, ch4, n20, 1991-2002
# "hf003-01-plant.csv"                tree canopy dimensions, 1991
# "hf069-16-allometries.txt"          allometry equations
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
# "hf069-13-annual-summ.csv" (annual woody) or...
# "hf069-09-trees.csv"
#   timescale: 1993-2017
#   variables: biomass_ag and modeled c_ag, n_ag, biomass_bg, c_bg, n_bg, anpp, biomass_sp, abundance_sp
# "hf006-01-soil-respiration.csv"
#   timescale: 1995-2010
#   variables: r_soil,
# "hf015-05-soils.csv"
#   timescale: 1937-2013
#   variables: c_so, n_so
# "hf004-02-filled.csv"
#   timescale: 1991-2015
#   variables: nee

#-------------------------------------------------------------------------------------------------#
# Time period:
#-------------------------------------------------------------------------------------------------#
# Intersection: 1996-2009 (14 years)
# Union:        1937-2017 (81 years)
# Climate:      2002-2012 (11 years)

# Example: download observation data
#url  <- "http://harvardforest.fas.harvard.edu/data/p06/hf069/hf069-09-trees.csv"
#dest <- file.path(base, "observations","hf_ems","hf069-09-trees.csv")
#download_file(url=url, destination=dest)

#-------------------------------------------------------------------------------------------------#
# Base directory
#-------------------------------------------------------------------------------------------------#
base_hf_data <- file.path(base, "observations", "hf_ems")

#-------------------------------------------------------------------------------------------------#
# Species look-up table
#-------------------------------------------------------------------------------------------------#
# http://harvardforest.fas.harvard.edu/sites/harvardforest.fas.harvard.edu/files/data/p06/hf069/hf069-17-allometries.txt
# http://harvardforest.fas.harvard.edu:8080/exist/apps/datasets/showData.html?id=hf003
# Using USDA PLANTS database codes: https://plants.usda.gov
species_lut <- data.frame(
  abbr  = c("AB","AC","ash","bb",
            "BB","BC","beech","bo",
            "cherry","chestnut","EH","FPC",
            "gb","haw","hbb","HBB",
            "hem","HT","nwr","NWR",
            "PB","rm","RM","ro",
            "RO","rp","sm","SM",
            "SPI","SS","wb","wh",
            "WH","wo","wp","WP",
            "ws","WS","yb","YB"),
  name  = c("Fagus grandifolia","Castenea dentate","Fraxinus americana","Betula lenta",
            "Betula lenta","Prunus serotina","Fagus grandifolia","Quercus velutina",
            "Prunus serotina","Castenea dentate","Tsuga canadensis","?",
            "Betula populifolia","Crataegus","Vaccinium","Vaccinium",
            "Tsuga canadensis","Crataegus","Viburnum","Viburnum",
            "Betula papyrifera","Acer rubrum","Acer rubrum","Quercus rubra",
            "Quercus rubra","Pinus resinosa","Acer pennsylvanicum","Acer pennsylvanicum",
            "Spiraea","?","Ilex verticillata","Hamamelis virginiania",
            "Hamamelis virginiania","Quercus alba","Pinus strobus","Pinus strobus",
            "Picea glauca","Picea glauca","Betula alleghaniensis","Betula alleghaniensis"),
  genus = c("Fagus","Castenea","Fraxinus","Betula","Betula","Prunus","Fagus","Quercus","Prunus",
            "Castenea","Tsuga","","Betula","Crataegus","Vaccinium","Vaccinium","Tsuga","Crataegus",
            "Viburnum","Viburnum","Betula","Acer","Acer","Quercus","Quercus","Pinus","Acer",
            "Acer","Spiraea","","Ilex","Hamamelis","Hamamelis","Quercus","Pinus","Pinus","Picea",
            "Picea","Betula","Betula"),
  code  = c("FAGR","CADE","FRAM","BELE",
            "BELE","PRSE","FAGR","QUVE",
            "PRSE","CADE","TSCA","",
            "BEPO","CRSP","VACO","VACO",
            "TSCA","CRSP","VICA","VICA",
            "BEPA","ACRU","ACRU","QURU",
            "QURU","PIRE","ACPE","ACPE",
            "","","ILVE","HAVI",
            "HAVI","QUAL","PIST","PIST",
            "PIGL","PIGL","BEAL","BEAL")
)
#message(paste("Species LUT:", jsonlite::toJSON(species_lut, pretty=TRUE)))

#-------------------------------------------------------------------------------------------------#
# Observation file list
#-------------------------------------------------------------------------------------------------#
observation_list <- list(
  climate      = file.path(base_hf_data, "observations", "hf_ems", "hf001-04-monthly-m.csv"),
  fluxes       = file.path(base_hf_data, "observations", "hf_ems", "hf004-02-filled.csv"),
  inventory    = file.path(base_hf_data, "observations", "hf_ems", "hf069-09-trees.csv"),
  root_cn      = file.path(base_hf_data, "observations", "hf_ems", "hf278-04-root-carbon-nitrogen.csv"),
  litter_cn    = file.path(base_hf_data, "observations", "hf_ems", "hf069-06-canopy-chem.csv"),
  soils        = file.path(base_hf_data, "observations", "hf_ems", "hf015-05-soils-cn.csv"),
  r_soil       = file.path(base_hf_data, "observations", "hf_ems", "hf006-01-soil-respiration.csv"),
  anpp         = file.path(base_hf_data, "observations", "hf_ems", "hf069-13-annual-summ.csv")
)

#-------------------------------------------------------------------------------------------------#
# Functions
#-------------------------------------------------------------------------------------------------#

# Tree species aboveground biomass equations; kg; DBH cm; Chojnacky, Heath, Jenkins (2014)
get_tree_biomass_ag <- function(code, dbh) {
  switch(code,                                   #  WSG
    ACPE = { exp(-2.0470 + 2.3852 * log(dbh)) }, # 0.44
    ACRU = { exp(-2.0470 + 2.3852 * log(dbh)) }, # 0.49
    BEAL = { exp(-1.8096 + 2.3480 * log(dbh)) }, # 0.55
    BELE = { exp(-1.8096 + 2.3480 * log(dbh)) }, # 0.60
    BEPA = { exp(-2.2271 + 2.4513 * log(dbh)) }, # 0.48
    BEPO = { exp(-2.2271 + 2.4513 * log(dbh)) }, # 0.45
    CADE = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.40
    FAGR = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.56
    FRAM = { exp(-1.8384 + 2.3524 * log(dbh)) }, # 0.55
    PIGL = { exp(-2.1364 + 2.3233 * log(dbh)) }, # 0.37
    PIRE = { exp(-2.6177 + 2.4638 * log(dbh)) }, # 0.41
    PIST = { exp(-2.6177 + 2.4638 * log(dbh)) }, # 0.34
    PRSE = { exp(-2.2118 + 2.4133 * log(dbh)) }, # 0.47
    QUAL = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.60
    QURU = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.56
    QUVE = { exp(-2.0705 + 2.4410 * log(dbh)) }, # 0.56
    TSCA = { exp(-2.3480 + 2.3876 * log(dbh)) }, # 0.38
    CRSP = { NA },                               # shrub
    HAVI = { NA },                               # shrub
    ILVE = { NA },                               # shrub
    VACO = { NA },                               # shrub
    VICA = { NA },                               # shrub
    NA                                           # else
  )
}

# Tree species belowground biomass equations; kg; DBH cm; Chojnacky, Heath, Jenkins (2014)
get_tree_biomass_bg <- function(code, dbh) {
  biomass_ag       = get_tree_biomass_ag(code, dbh)
  root_coarse_frac = exp(-1.4485 - 0.03476 * log(dbh))
  root_fine_frac   = exp(-1.8629 - 0.77534 * log(dbh))
  soil_frac        = 0.68
  biomass_bg       = (biomass_ag * root_coarse_frac) + (biomass_ag * root_fine_frac) + (biomass_ag * soil_frac)
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

# NEE; kg C m-2 year-1; Converted from µmol CO2 m-2 sec-1 (10^-6 moles)
get_nee <- function() {
  #df = read.csv(self$files$observations$flux, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf004-02-filled.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df$c_kg = df$nee * 1e-6 * 12 / 1000       # convert µmol CO2 to mol CO2 to g C to kg C
  df$c_kg = df$c_kg * 60 * 60 * 24 * 365.24 # convert kg C m-2 second-1 to kg C m-2 year-1
  return(setNames(aggregate(df$c_kg, by=list(df$year), FUN=mean), c("year","nee")))
}

# Aboveground biomass; kg m-2; 34 EMS tower plots (10m radius; (pi*10^2)*34 = 10681.42 m2)
get_biomass_ag <- function() {
  #df = read.csv(self$files$observations$inventory, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf069-09-trees.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  #df = df[df$site =="ems" & df$tree.type == "live",]
  df = df[df$site =="ems" & df$tree.type == "live" & !is.na(df$year) & df$year > 1993,]
  df$code = species_lut$code[match(df$species, species_lut$abbr)]
  df = df[df$code %in% c("ACPE","ACRU","BEAL","BELE","BEPO","FAGR","FRAM","PIGL",
                         "PIRE","PIST","PRSE","QURU","QUVE","TSCA"),] # modeled species
  # Mean annual individual tree DBH
  df = setNames(aggregate(df$dbh, by=list(df$year, df$code, df$plot, df$tag), FUN=mean, na.rm=TRUE),
                c("year","code","plot","tag","dbh"))
  df = df[,c("year","code","dbh")]
  df$biomass_ag = apply(df, 1, function(x) get_tree_biomass_ag(x[2], as.numeric(x[3])))
  agb = setNames(aggregate(df$biomass_ag, by=list(df$year), FUN=sum, na.rm=TRUE), c("year","biomass_ag"))
  agb$biomass_ag = agb$biomass_ag / 10681.42 # convert kg to kg m-2
  return(agb)
}

# Aboveground carbon; kg C m-2; 34 tower plots (10m radius; (pi*10^2)*34 = 10681.42 m2)
get_c_ag <- function() {
  df      = get_biomass_ag()
  df$c_ag = df$biomass_ag * 0.5
  return(df[,c("year","c_ag")])
}

# Aboveground carbon; kg N m-2; 34 tower plots (10m radius; (pi*10^2)*34 = 10681.42 m2)
get_n_ag <- function() {
  agc = get_c_ag()
  #df = read.csv(self$files$observations$root_cn, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf278-04-root-carbon-nitrogen.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df$cn   = df$c.n / 100
  cn_stem = mean(df$cn, na.rm=TRUE)
  #df = read.csv(self$files$observations$leaf_cn, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf069-06-canopy-chem.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df$cn_leaf = df$c.percent / df$n.percent
  cn_leaf    = mean(df$cn_leaf, na.rm=TRUE)
  agc$n_ag   = agc$c_ag / ((cn_stem * 0.80) + (cn_leaf * 0.20)) # assumes general 80/20 stem/leaf fraction
  return(agc[,c("year","n_ag")])
}

# Belowground biomass; kg m-2; 34 tower plots (10m radius; (pi*10^2)*34 = 10681.42 m2)
# Uses general AGB-to-BGB model of Chojnacky, Heath, Jenkins (2014)
get_biomass_bg <- function() {
  df = read.csv(file.path(base_hf_data, "hf069-09-trees.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site =="ems" & df$tree.type == "live",]
  df$code = species_lut$code[match(df$species, species_lut$abbr)]
  df = df[df$code %in% c("ACPE","ACRU","BEAL","BELE","BEPO","FAGR","FRAM","PIGL",
                         "PIRE","PIST","PRSE","QURU","QUVE","TSCA"),] # modeled species
  # Mean annual individual tree DBH
  df = setNames(aggregate(df$dbh, by=list(df$year, df$code, df$plot, df$tag), FUN=mean, na.rm=TRUE),
                c("year","code","plot","tag","dbh"))
  df = df[,c("year","code","dbh")]
  df$biomass_ag = apply(df, 1, function(x) get_tree_biomass_ag(x[2], as.numeric(x[3])))
  #df$root_coarse = df$biomass_ag * 0.22
  #df$root_fine   = df$biomass_ag * 0.01
  df$root_coarse = df$biomass_ag * exp(-1.4485 - 0.03476 * log(df$dbh))
  df$root_fine   = df$biomass_ag * exp(-1.8629 - 0.77534 * log(df$dbh))
  df$soil        = df$biomass_ag * 0.68
  df$biomass_bg  = df$root_coarse + df$root_fine + df$soil
  bgb = setNames(aggregate(df$biomass_bg, by=list(df$year), FUN=sum, na.rm=TRUE), c("year","biomass_bg"))
  bgb$biomass_bg = bgb$biomass_bg / 10681.42 # convert kg to kg m-2
  return(bgb[,c("year","biomass_bg")])
}

# Belowground C; kg C m-2; 34 tower plots (10m radius; (pi*10^2)*34 = 10681.42 m2)
get_c_bg <- function () {
  bgb = get_biomass_bg()
  #df = read.csv(self$files$observations$root_cn, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf278-04-root-carbon-nitrogen.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  bgb$c_bg = bgb$biomass_bg * (mean(df$c.per, na.rm=TRUE) / 100)
  return(bgb[,c("year","c_bg")])
}

# Belowground N; kg N m-2; 34 tower plots (10m radius; (pi*10^2)*34 = 10681.42 m2)
get_n_bg <- function () {
  bgb = get_biomass_bg()
  #df = read.csv(self$files$observations$root_cn, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf278-04-root-carbon-nitrogen.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  bgb$n_bg = bgb$biomass_bg * (mean(df$n.per, na.rm=TRUE) / 100)
  return(bgb[,c("year","n_bg")])
}

# Soil organic C; g C cm-3; 240 (22.5 x 22.5m) 506.25 m2 plots (121500 m2 total) of 269 original plots
get_c_so <- function() {
  #df = read.csv(self$files$observations$soils, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf015-05-soils-cn.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  bd   = df$bdgcc.corr           # bulk density or g soil per cm3
  c_gg = df$gcgsoil              # g C per g soil
  c_so = bd * c_gg / 1000 * 1e+6 # g C per cm3 soil converted to kg per m3
  c_so = mean(c_so, na.rm=TRUE)
  return(data.frame(year=NA, c_so=c_so))
}

# Soil organic N; g N cm-3; 240 (22.5 x 22.5m) 506.25 m2 plots (121500 m2 total) of 269 original plots
get_n_so <- function() {
  #df = read.csv(self$files$observations$soils, header=TRUE, stringsAsFactors=FALSE)
  df  = read.csv(file.path(base_hf_data, "hf015-05-soils-cn.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  bd   = df$bdgcc.corr           # bulk density or g soil per cm3
  n_gg = df$gngsoil              # g N per g soil
  n_so = bd * n_gg / 1000 * 1e+6 # g C per cm3 soil converted to kg per m3
  n_so = mean(n_so, na.rm=TRUE)
  return(data.frame(year=NA, n_so=n_so))
}

# Soil respiration; mg C m-2 hour-1; kg C m-2 year-1; 6 sites near EMS tower
get_r_soil <- function() {
  #df = read.csv(self$files$observations$r_soil, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf006-01-soil-respiration.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df$date = as.POSIXlt(df$date, format="%m/%d/%Y", tz="EST")
  df$year = as.numeric(format(df$date, "%Y"))
  df$flux = df$flux / 1e+6 * 24 * 365.24 # convert mg C m-2 hour-1 to kg C m-2 year-1
  return(setNames(aggregate(df$flux, by=list(df$year), FUN=mean, na.rm=TRUE), c("year","r_soil")))
}

# Annual net primary productivity (ANPP); Mg year-1; kg C m-2 year-1
get_anpp <- function() {
  #df = read.csv(self$files$observations$anpp, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf069-13-annual-summ.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site =="ems",]
  df = df[!is.na(df$anpp),]
  df$anpp = df$anpp * 1000 / 10681.42 * 0.5 # Mg mass year-1 to kg C m-2 year-1
  return(df[,c("year","anpp")])
}

# Total species aboveground biomass; kg m-2; 34 tower plots (10m radius; (pi*10^2)*34 = 10681.42 m2)
get_biomass_sp <- function() {
  #df = read.csv(self$files$observations$inventory, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf069-09-trees.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site =="ems" & df$tree.type == "live",]
  df$code  = species_lut$code[match(df$species, species_lut$abbr)]
  df       = df[df$code %in% c("ACPE","ACRU","BEAL","BELE","BEPO","FAGR","FRAM","PIGL",
                               "PIRE","PIST","PRSE","QURU","QUVE","TSCA"),] # modeled species
  # Mean annual individual tree DBH
  df = setNames(aggregate(df$dbh, by=list(df$year, df$code, df$plot, df$tag), FUN=mean, na.rm=TRUE),
                c("year","code","plot","tag","dbh"))
  df = df[,c("year","code","dbh")]
  df$biomass_ag = apply(df, 1, function(x) get_tree_biomass_ag(x[2], as.numeric(x[3])))
  biomass_sp    = aggregate(df$biomass_ag, by=list(df$year, df$code), FUN=sum, na.rm=TRUE)
  colnames(biomass_sp) = c("year","species","biomass_ag")
  biomass_sp$biomass_ag = biomass_sp$biomass_ag / 10681.42
  return(biomass_sp)
}

# Species relative abundance; percent of n-trees
get_abundance_sp <- function() {
  #df = read.csv(self$files$observations$inventory, header=TRUE, stringsAsFactors=FALSE)
  df = read.csv(file.path(base_hf_data, "hf069-09-trees.csv"), header=TRUE, stringsAsFactors=FALSE)
  colnames(df) = tolower(colnames(df))
  df = df[df$site =="ems" & df$tree.type == "live",]
  df$code        = species_lut$code[match(df$species, species_lut$abbr)]
  df             = df[df$code %in% c("ACPE","ACRU","BEAL","BELE","BEPO","FAGR","FRAM","PIGL",
                                     "PIRE","PIST","PRSE","QURU","QUVE","TSCA"),] # modeled species
  ntrees_yr      = setNames(aggregate(df$code, by=list(df$year), FUN=length), c("year","ntrees"))
  df$ntrees_yr   = ntrees_yr$ntrees[match(df$year, ntrees_yr$year)]
  df$fraction_yr = 1 / df$ntrees_yr # relative abundance of each tree
  return(setNames(aggregate(df$fraction_yr, by=list(df$year, df$code), FUN=sum), c("year","species","abundance")))
}
