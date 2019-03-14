#!/usr/bin/env Rscript
#----------------------------------------------------------------------------------------------------------------#
# Filename:
#   ppa_sibgc_v1.r
# Description:
#   R implementation of the Perfect Plasticity Approximation (PPA) with simply biogeochemistry (PPA-SiBGC).
# Version:
#   1.0
# Authors:
#   Adam Erickson,   Washington State University
#   Nikolay Strigul, Washington State University
# Modified:
#   June 30, 2018
# Procedure:
#   (1) Generate constant growth, mortality, and regeneration coefficients from plot data
#   (2) Generate above- and below-ground stoichiometric coefficients from databases
#   (3) Run the Sortie-PPA model:
#     (*) Generate cohorts from tree inventory data
#     (*) Initialize values for allometry, biomass, C, and N from species, type, DBH
#     (a) Soil respiration
#     (b) Regeneratation
#     (c) Sort by height descending
#     (d) Calculate cumulative crown area (CCA)
#     (e) Calculate Z* height and nearest cohort
#     (f) Mortality
#     (g) Growth
#     (h) Allometry
#     (i) Biomass fractions
#     (j) Carbon fractions
#     (k) Nitrogen fractions
#     (l) Save yearly results to CSV
#     (m) Save final results to CSV
#   (4) Compare model results against empirical data
# Notes:
#   Sortie-PPA is based on the assumption that tree canopies are perfectly plastic. This reduces
#   the variance in canopy join heights to zero, creating a single canopy join height known as z*.
#   The approximation allows for model simplification and mathematical tractability. The spatial
#   location of trees is discarded. Growth and mortality are modeled by their mean values for each
#   species and type (adults = above z*; saplings = below z*). Coefficients are used to modify
#   growth rates per the fraction of light received by trees above and below the z* threshold.
#   Aboveground biomass is modeled using equations from Chojnacky, Heath, and Jenkins (2014).
#   Belowground biomass as well as C,N dynamics are modeled using allometry and stoichiometry.
# Examples:
#   Rscript --vanilla ppa_v50.r
#   Rscript --vanilla ppa_v50.r --wd /Users/null/test --verbose
#----------------------------------------------------------------------------------------------------------------#
Version <- "5.0"

# Debugging
options(error=quote({dump.frames(to.file=TRUE); q()}))

# Print banner and pause
Banner <- readLines("sortie_ppa_ansi_shadow.txt")
cat(Banner, sep="\n")

# Print version
message("")
message("Authors: Adam Erickson, Nikolay Strigul")
message(paste("Version:", Version, "\n"))
Sys.sleep(0)

# Remove existing outputs
outputs_files   <- dir("outputs", recursive=FALSE, full.names=TRUE)
outputs_deleted <- file.remove(outputs_files)

# Fetch trailing command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Default verbosity setting
Verbose <- FALSE

# Set working directory
if (length(args) == 0) {
  message(paste0("Working directory: ", getwd()))
} else if (length(args) == 1) {
  Verbose <- args[1] %in% c("--v", "--verbose")
  stopifnot(Verbose)
} else if (length(args) == 2) {
  flag   <- args[1]
  stopifnot(flag %in% c("--wd", "--dir"))
  folder <- args[2]
  if (flag=="--wd" || flag=="--dir") {
    setwd(folder)
    message(paste("Working directory:", getwd()))
  } else {
    message("Flag not recognized.")
  }
} else if (length(args) == 3) {
  flag   <- args[1]
  stopifnot(flag %in% c("--wd", "--dir"))
  folder  <- args[2]
  setwd(folder)
  message(paste("Working directory:", getwd()))
  Verbose <- args[3] %in% c("--v", "--verbose")
  stopifnot(Verbose)
} else {
  stop("Error: Only 0 or 2 arguments may be passed.", call.=FALSE)
}

# Libraries #

# Functions #

# Round to the nearest integer by step
round_int <- function(x, step) { 
  return(step * round(x / step))
}

# Bin trees into cohorts by rounding DBH to the nearest cm; check n-trees
# Expects CSV file with species, type, dbh, age (optional)
GenerateCohorts <- function(Trees, interval=1) {
  n_trees <- nrow(Trees)
  Trees$dbh <- round_int(Trees$dbh, interval)
  if("age" %in% tolower(colnames(Trees))) {
    age <- aggregate(Trees$age, by=list(Trees$species, Trees$type, Trees$dbh), FUN=mean)
    colnames(age) <- c("species","type","dbh","age")
    Trees <- subset(Trees, select= -c(age))
    Trees <- merge(Trees, age, by=c("species","type","dbh"))
  }
  cohorts <- data.frame(table(Trees$species, Trees$type, Trees$dbh))
  colnames(cohorts) <- c("species","type","dbh","n_trees")
  cohorts$n_trees <- cohorts$n_trees + 1
  Trees <- merge(Trees, cohorts, by=c("species","type","dbh"))
  Trees <- Trees[!duplicated(Trees[,c("species","type","dbh")]), ]
  Trees <- Trees[Trees$dbh > 0 & !is.na(Trees$n_trees), ]
  Trees$id <- 1:nrow(Trees)
  if ("age" %in% tolower(colnames(Trees))) {
    Trees <- Trees[, c("id","species","type","age","dbh","n_trees")]
  } else {
    Trees <- Trees[, c("id","species","type","dbh","n_trees")]
  }
  message(paste(n_trees, "trees merged to", nrow(Trees), "cohorts containing",
                sum(Trees$n_trees), "trees"))
  return(Trees)
}

# Create biomass, C, and N compartments
Initialize <- function(Trees, ...) {
  # Create columns
  Trees$height              <- NA
  Trees$crown_l             <- NA
  Trees$crown_r             <- NA
  Trees$crown_a             <- NA
  Trees$ba                  <- NA
  Trees$crown_a_cumsum      <- NA
  Trees$biomass_ag          <- NA
  Trees$biomass_stem        <- NA
  Trees$biomass_branch      <- NA
  Trees$biomass_leaf        <- NA
  Trees$biomass_root_coarse <- NA
  Trees$biomass_root_fine   <- NA
  Trees$biomass_soil        <- NA
  Trees$biomass_bg          <- NA
  Trees$biomass_total       <- NA
  Trees$c_ag                <- NA
  Trees$c_bg                <- NA
  Trees$c_total             <- NA
  Trees$anpp                <- NA
  Trees$n_ag                <- NA
  Trees$n_bg                <- NA
  Trees$n_total             <- NA
  Trees$annp                <- NA
  # Update each species and type; vectorized
  for (species in unique(Trees$species)) {
    if (!species %in% AllometryLUT$species) { stop("Species not found in allometry table") }
    s <- Trees$species==species
    for (type in unique(Trees[s,]$type)) {
      st <- Trees$species==species & Trees$type==type
      index <- AllometryLUT$species==species & AllometryLUT$type==type
      allometry <- AllometryLUT[index,]
      Trees[st,]$height  <- 1.35 + (30-1.35) * (1-exp(-(allometry$h_coeff * Trees[st,]$dbh)))
      Trees[st,]$height  <- ifelse(Trees[st,]$height > 30, 30, Trees[st,]$height) # limit height to 30m
      Trees[st,]$type    <- ifelse(Trees[st,]$dbh > 10, "adult", "sapling")
      Trees[st,]$crown_l <- allometry$cd  * Trees[st,]$dbh
      Trees[st,]$crown_r <- allometry$cr1 * Trees[st,]$dbh^allometry$cr2
      Trees[st,]$crown_a <- Trees[st,]$crown_r^2 * pi
      Trees[st,]$ba      <- Trees[st,]$dbh^2 * 0.00007854 # forester's constant for DBH (cm); m2
    }
    # Biomass allocation or partitioning
    if (!species %in% BiomassLUT$species) { stop("Species not found in biomass table") }
    index <- BiomassLUT$species==species
    biomass <- BiomassLUT[index,]
    if (CohortMode) {
      Trees[s,]$biomass_ag <- exp(biomass$b0 + biomass$b1 * log(Trees[s,]$dbh)) * Trees[s,]$n_trees
    } else {
      Trees[s,]$biomass_ag <- exp(biomass$b0 + biomass$b1 * log(Trees[s,]$dbh))
    }
    Trees[s,]$biomass_stem        <- Trees[s,]$biomass_ag * biomass$fraction_stem
    Trees[s,]$biomass_branch      <- Trees[s,]$biomass_ag * biomass$fraction_branch
    Trees[s,]$biomass_leaf        <- Trees[s,]$biomass_ag * biomass$fraction_leaf
    Trees[s,]$biomass_root_coarse <- Trees[s,]$biomass_ag * exp(-1.4485 - 0.03476 * log(Trees[s,]$dbh))
    Trees[s,]$biomass_root_fine   <- Trees[s,]$biomass_ag * exp(-1.8629 - 0.77534 * log(Trees[s,]$dbh))
    Trees[s,]$biomass_soil        <- Trees[s,]$biomass_ag * biomass$fraction_soil
    Trees[s,]$biomass_bg          <- Trees[s,]$biomass_root_coarse + Trees[s,]$biomass_root_fine +
                                     Trees[s,]$biomass_soil
    Trees[s,]$biomass_total       <- Trees[s,]$biomass_ag + Trees[s,]$biomass_bg
    # C fraction
    #index <- CarbonLUT$species==species
    carbon <- CarbonLUT #[index,]
    Trees[s,]$c_ag <- Trees[s,]$biomass_ag *
      ((biomass$fraction_stem   * carbon$fraction_stem)   +
       (biomass$fraction_branch * carbon$fraction_branch) +
       (biomass$fraction_leaf   * carbon$fraction_leaf))
    Trees[s,]$c_bg <- Trees[s,]$biomass_ag *
      ((biomass$fraction_root * carbon$fraction_root) +
       (biomass$fraction_soil * carbon$fraction_soil))
    Trees[s,]$c_total <- Trees[s,]$c_ag + Trees[s,]$c_bg
    Trees[s,]$anpp    <- NA
    # C:N stoichiometry
    if (!species %in% StoichiometryLUT$species) { stop("Species not found in stoichiometry table") }
    index <- StoichiometryLUT$species==species
    stoichiometry <- StoichiometryLUT[index,]
    Trees[s,]$n_ag <- Trees[s,]$c_ag /
      ((biomass$fraction_stem   * stoichiometry$cn_stem)   +
       (biomass$fraction_branch * stoichiometry$cn_branch) +
       (biomass$fraction_leaf   * stoichiometry$cn_leaf))
    Trees[s,]$n_bg <- Trees[s,]$c_bg /
      ((biomass$fraction_root * stoichiometry$cn_root)  +
       (biomass$fraction_soil * stoichiometry$cn_soil))
    Trees[s,]$n_total <- Trees[s,]$n_ag + Trees[s,]$n_bg
    Trees[s,]$annp    <- NA
  }
  return(Trees)
}

# Ensure column names and tree type are lower case
ensure_lowercase <- function(Trees) {
  colnames(Trees) <- tolower(colnames(Trees))
  Trees$type <- tolower(Trees$type)
  return(Trees)
}

# Check if year is a leap year (366 days instead of 365)
leap_year <- function(year) {
  return(ifelse((year %%4 == 0 & year %%100 != 0) | year %%400 == 0, TRUE, FALSE))
}

# Soil respiration, monthly mean; heterotrophic + autotrophic
# g C m-2 day-1; kg C m-2 day-1
# Parameters:
#  Ta = monthly mean temperature (°C)
#  P  = monthly mean precipitation (cm)
# Raich, Potter, Bhagawati, 2002, Interannual variability in global soil respiration, 1980-94, Global Change Biology, 8, pp. 800-12.
# https://onlinelibrary.wiley.com/doi/abs/10.1046/j.1365-2486.2002.00511.x
respiration_soil <- function(Ta, P) {
  if (Ta < -13.3) {
    Rs <- 0
  } else {
    if (Ta > 33.5) { Ta <- 33.5 } # bound the temperature response
    e  <- exp(1)  # Euler's constant
    Q  <- 0.05452 # soil respiration rate change wrt temperature (°C-1)
    F  <- 1.250   # soil respiration rate at 0°C (g C m-2 day-1)
    K  <- 4.259   # hyperbolic respiration-rainfall curve half-saturation constant (mm month-1)
    Rs <- F * e^(Q * Ta) * (P / (K + P))
  }
  return(Rs / 1000)
}

# Simple SOC model
# Mg C ha-1 cm-1; kg C m-2; midpoint of profile depth; Approximates integral in 1 cm steps
# Domke et al. (2017) Toward inventory-based estimates of soil organic carbon in forests of the United States, Ecological Applications, 27(4), pp. 1223–1235.
# https://www.fs.fed.us/nrs/pubs/jrnl/2017/nrs_2017_domke_001.pdf
soc_depth = function(order, depth_cm=100) {
  soc_table = data.frame(
    order     = c("All","Alfisols","Andisols","Aridisols","Entisols","Histosols","Inceptisols",
                  "Mollisols","Spodosols","Ultisols","Vertisols"),
    intercept = c(1.1795,1.1122,1.3837,0.2065,0.9300,1.6227,1.1631,1.0163,1.4262,1.1576,0.5145),
    slope     = c(-0.8228,-0.8330,-0.8425,-0.1300,-0.7207,-1.0109,-0.7331,-0.6214,-0.9801,-0.8867,-0.2427)
  )
  coeffs = as.numeric(soc_table[soc_table$order==order,][,c("intercept","slope")])
  soc    = sapply(seq(1, depth_cm, 1), function(x) 10^(coeffs[1] + coeffs[2] * log10(x)))
  soc    = sum(soc) / 10  # convert to kg m-2 integrated over depth
  return(soc)
}

# Time #
start_time <- proc.time()
loop_times <- c()

# Data #

# Load lookup tables (LUTs) into memory; expects CSV files in lut folder
Configuration    <- read.csv("configuration.csv",     stringsAsFactors=FALSE, nrows=1)
Trees            <- read.csv("trees.csv",             stringsAsFactors=FALSE)
Climate          <- read.csv("climate.csv",           stringsAsFactors=FALSE)
GrowthLUT        <- read.csv("lut/growth.csv",        stringsAsFactors=FALSE)
MortalityLUT     <- read.csv("lut/mortality.csv",     stringsAsFactors=FALSE)
RegenerationLUT  <- read.csv("lut/regeneration.csv",  stringsAsFactors=FALSE)
AllometryLUT     <- read.csv("lut/allometry.csv",     stringsAsFactors=FALSE)
BiomassLUT       <- read.csv("lut/biomass.csv",       stringsAsFactors=FALSE)
StoichiometryLUT <- read.csv("lut/stoichiometry.csv", stringsAsFactors=FALSE)
CarbonLUT        <- read.csv("lut/carbon.csv",        stringsAsFactors=FALSE, nrows=1)

# Model #

# Optional: process in parallel
#nCores <- as.integer(Configuration$n_cores)
#if (nCores > 1) {
#  library(doParallel)
#  cl <- parallel::makeCluster(nCores)
#  doParallel::registerDoParallel(cl)
#}

# Optional: Generate cohorts from initial tree list
CohortMode <- as.logical(Configuration$cohort_mode)
if (CohortMode == TRUE) {
  message("Generating cohorts...")
  Trees <- GenerateCohorts(Trees)
}

# Initialize tree biomass
message("Initializing...")
Trees <- Initialize(Trees)

# Ensure column names and tree type are lower case
Trees <- ensure_lowercase(Trees)

# Number of years to run model; e.g., 2000-2049 (50 years)
StartYear <- as.numeric(Configuration$start_year) 
EndYear   <- as.numeric(Configuration$end_year)
Years     <- StartYear:EndYear # 200 years

# Amount of light suppressed cohorts receive
UnderstoryAdult   <- as.numeric(Configuration$understory_adult)   # 0.1
UnderstorySapling <- as.numeric(Configuration$understory_sapling) # 0.7

# id placeholder
id_holder <- nrow(Trees)

# Run the simulation
message("Starting simulation...")
message("=========================================================================")
for (year in Years) {
  
  loop_start_time <- proc.time()
  
  if (Verbose) {
    message("=========================================================================")
    message(paste("Year:", year))
    message(paste("n-Trees:", nrow(Trees)))
    message("=========================================================================")
  }

  # Append or update year; vectorized
  Trees$year <- year

  # Increment tree age; vectorized; regeneration yields age 0
  Trees$age <- Trees$age + 1

  # Calculate soil respiration
  if (Verbose) {
    message("Computing soil respiration...")
  }
  climate  <- Climate[Climate$year==year,]
  n_days <- list(
    "1"  = 31,
    "2"  = ifelse(leap_year(year), 29, 28),
    "3"  = 31,
    "4"  = 30,
    "5"  = 31,
    "6"  = 30,
    "7"  = 31,
    "8"  = 31,
    "9"  = 30,
    "10" = 31,
    "11" = 30,
    "12" = 31
  )
  Rs_month <- c()
  for (month in climate$month) {
    month_days <- as.numeric(n_days[as.character(month)])
    month_i    <- climate[climate$month==month,]
    Rs         <- respiration_soil(Ta=month_i$temperature, P=month_i$precipitation)
    Rs         <- Rs * month_days
    Rs_month   <- c(Rs_month, Rs)
  }
  Rs_year        <- sum(Rs_month) # * Configuration$field_area
  respiration_df <- data.frame(year=year, r_soil=Rs_year)

  # Calculate soil organic C and N
  c_so   <- soc_depth(order=Configuration$soil_order, depth_cm=100)
  n_so   <- c_so / mean(StoichiometryLUT$cn_soil, na.rm=TRUE)
  som_df <- data.frame(year=year, c_so=c_so, n_so=n_so)

  # Sort trees by descending height; vectorized
  if (Verbose) {
    message("Sorting...")
  }
  Trees <- Trees[order(Trees$height, decreasing=TRUE), ]

  # Calculate cumulative crown area; vectorized
  if (Verbose) {
    message("Calculating culumative crown area...")
  }
  if (CohortMode) {
    Trees$crown_a_cumsum <- cumsum(Trees$crown_a * Trees$n_trees)
  } else {
    Trees$crown_a_cumsum <- cumsum(Trees$crown_a)
  }

  # Calculate z* and closest cohort; vectorized
  if (Verbose) {
    message("Calculating Z*...")
  }
  index    <- which.min(abs(Trees$crown_a_cumsum - Configuration$field_area))
  zstar_df <- Trees[index,]
  if (Verbose) {
    message(paste("Z* height:", zstar_df$height))
  }

  # Dataframe placeholders
  mortality_df <- data.frame(matrix(ncol=nrow(Trees), nrow=0))
  colnames(mortality_df) <- colnames(Trees)
  regeneration_df <- data.frame(matrix(ncol=nrow(Trees), nrow=0))
  colnames(regeneration_df) <- colnames(Trees)

  # Loop over each species
  for (species in unique(Trees$species)) {

    # Skip species not found in lookup tables
    if (!species %in% GrowthLUT$species) {
      next
    }

    # Species index
    s <- Trees$species==species
    
    # Regeneration; species
    if (Verbose) {
      message("Applying regeneration...")
    }
    n_replicates <- length(which(s))
    index <- RegenerationLUT$species==species
    n_saplings <- RegenerationLUT[index,]$mean
    # new sapling used in regeneration
    new_sapling <- data.frame(species=species, type="sapling", age=0, dbh=1,
      height=2, crown_r=0.1, crown_l=0.845, crown_a=0.0314, ba=0.00008,
      crown_a_cumsum=NA, biomass_ag=NA, biomass_stem=NA, biomass_branch=NA,
      biomass_leaf=NA, biomass_root_coarse=NA, biomass_root_fine=NA,
      biomass_soil=NA, biomass_bg=NA, biomass_total=NA, c_ag=NA, c_bg=NA,
      c_total=NA, anpp=NA, n_ag=NA, n_bg=NA, n_total=NA, annp=NA, year=year)
    if (CohortMode & n_replicates > 0) {
      new_sapling$id <- id_holder
      new_sapling$n_trees <- n_saplings
      Trees <- rbind(Trees, new_sapling)
      regeneration_df <- rbind(regeneration_df, new_sapling)
      id_holder <- id_holder + 1 #+ n_replicates
    } else if (n_replicates > 0){
      n_trees <- nrow(trees)
      new_saplings <- do.call(rbind, replicate(n_replicates * n_saplings, new_sapling, simplify=FALSE))
      new_saplings$id <- id_holder:(id_holder + (n_replicates * n_saplings) - 1)
      Trees <- rbind(Trees, new_saplings)
      regeneration_df <- rbind(regeneration_df, new_saplings)
      id_holder <- id_holder + (n_replicates * n_saplings)
    }

    # Recalculate species index after regeneration
    s <- Trees$species==species

    # Loop over each type within each species
    for (type in unique(Trees[s,]$type)) {

      # Species and type index
      st <- Trees$species==species & Trees$type==type

      # Mortality; species & type
      if (Verbose) {
        message("Applying mortality...")
      }
      index <- MortalityLUT$species==species & MortalityLUT$type==type
      mortality <- MortalityLUT[index,]
      random_uniform <- runif(1, min=0, max=1)
      if (random_uniform < mortality$probability) {
        kill_fraction <- runif(1, min=0, max=1) # random uniform
        if (CohortMode) {
          kill_n <- round(nrow(Trees[st,]) * kill_fraction)
          if (kill_n > 0) {
            # Kill fraction of cohorts
            kill_index <- sample(which(st), kill_n, replace=FALSE)
            mortality_df <- rbind(mortality_df,  Trees[kill_index,])
            Trees <- Trees[-kill_index,] # kill trees
            # Optional: kill fraction within cohorts and then remove dead cohorts
            #Trees[st,]$n_trees <- Trees[st,]$n_trees - kill_n # kill trees
            #mortality_df <- cbind(data.frame(year=year), Trees[st,])
            #mortality_df$n_trees <- round(mortality_df$n_trees * kill_fraction)
          }
          #dead_cohorts <- Trees$n_trees < 1 | is.na(Trees$n_trees) | is.na(Trees$dbh)
          #if(any(dead_cohorts)) {
          #  Trees <- Trees[-which(dead_cohorts),] # remove dead trees/cohorts
          #  message("Removal")
          #  stopifnot(!anyNA(Trees[st,]$dbh) | !anyNA(Trees[st,]$n_trees)) # "PIAB" "sapling"
          #}
        } else {
          kill_n <- round(nrow(Trees[st,]) * kill_fraction)
          if (kill_n > 0) {
            kill_index <- sample(which(st), kill_n, replace=FALSE)
            mortality_df <- rbind(mortality_df, Trees[kill_index,])
            Trees <- Trees[-kill_index,] # kill trees
          }
        }
      }

      # Recalculate species and type index after mortality
      st <- Trees$species==species & Trees$type==type

      # Skip ahead if species or type is all dead
      if (nrow(Trees[st,]) < 1) {
        next
      }

      # Growth; species & type
      index <- GrowthLUT$species==species & GrowthLUT$type==type
      growth <- GrowthLUT[index,]
      if (type=="adult") {
        Trees[st,]$dbh <- ifelse(Trees[st,]$crown_a_cumsum < Configuration$field_area,
                                 Trees[st,]$dbh + growth$mean,
                                 Trees[st,]$dbh + growth$mean * UnderstoryAdult)
      } else if (type=="sapling") {
        Trees[st,]$dbh <- Trees[st,]$dbh + growth$mean * UnderstorySapling
      } else {
        message("Error: Type not recognized")
      }

      # Allometry; species & type
      if (Verbose) {
        message("Calculating allometry...")
      }
      index <- AllometryLUT$species==species & AllometryLUT$type==type
      allometry <- AllometryLUT[index,]
      Trees[st,]$height  <- 1.35 + (30-1.35) * (1-exp(-(allometry$h_coeff * Trees[st,]$dbh)))
      Trees[st,]$height  <- ifelse(Trees[st,]$height > 30, 30, Trees[st,]$height) # limit height to 30m
      Trees[st,]$type    <- ifelse(Trees[st,]$dbh > 10, "adult", "sapling")
      Trees[st,]$crown_l <- allometry$cd  * Trees[st,]$dbh
      Trees[st,]$crown_r <- allometry$cr1 * Trees[st,]$dbh^allometry$cr2
      Trees[st,]$crown_a <- Trees[st,]$crown_r^2 * pi
      Trees[st,]$ba      <- Trees[st,]$dbh^2 * 0.00007854 # forester's constant for DBH (cm); m2
    } # end type loop

    # Calculate species index after regeneration and mortality
    s <- Trees$species==species

    # Skip ahead if species is all dead
    if (nrow(Trees[s,]) < 1) {
      next
    }

    # Biomass expansion factors and allocation or partitioning; species
    # DBH to AGB and AGB to BGB (Chojnacky, Heath, Jenkins, 2014)
    # To add...? BGB to leaf and steam (Enquist and Niklas, 2002)
    if (Verbose) {
      message("Calculating biomass...")
    }
    index <- BiomassLUT$species==species
    biomass <- BiomassLUT[index,]
    if (CohortMode) {
      Trees[s,]$biomass_ag <- exp(biomass$b0 + biomass$b1 * log(Trees[s,]$dbh)) * Trees[s,]$n_trees
    } else {
      Trees[s,]$biomass_ag <- exp(biomass$b0 + biomass$b1 * log(Trees[s,]$dbh))
    }
    Trees[s,]$biomass_stem        <- Trees[s,]$biomass_ag * biomass$fraction_stem
    Trees[s,]$biomass_branch      <- Trees[s,]$biomass_ag * biomass$fraction_branch
    Trees[s,]$biomass_leaf        <- Trees[s,]$biomass_ag * biomass$fraction_leaf
    Trees[s,]$biomass_root_coarse <- Trees[s,]$biomass_ag * exp(-1.4485 - 0.03476 * log(Trees[s,]$dbh))
    Trees[s,]$biomass_root_fine   <- Trees[s,]$biomass_ag * exp(-1.8629 - 0.77534 * log(Trees[s,]$dbh))
    Trees[s,]$biomass_soil        <- Trees[s,]$biomass_ag * biomass$fraction_soil
    Trees[s,]$biomass_bg          <- Trees[s,]$biomass_root_coarse + Trees[s,]$biomass_root_fine +
                                     Trees[s,]$biomass_soil
    Trees[s,]$biomass_total       <- Trees[s,]$biomass_ag + Trees[s,]$biomass_bg

    # Year naught total C and N for ANPP and ANNP calculations
    c_total_previous <- Trees[s,]$c_total
    n_total_previous <- Trees[s,]$n_total

    # C fraction; species
    if (Verbose) {
      message("Calculating C fraction...")
    }
    #index <- CarbonLUT$species==species
    carbon <- CarbonLUT #[index,]
    Trees[s,]$c_ag <- Trees[s,]$biomass_ag * 
      (biomass$fraction_stem   * carbon$fraction_stem) +
      (biomass$fraction_branch * carbon$fraction_branch) +
      (biomass$fraction_leaf   * carbon$fraction_leaf)
    Trees[s,]$c_bg <- Trees[s,]$biomass_ag *
      (biomass$fraction_root * carbon$fraction_root) +
      (biomass$fraction_soil * carbon$fraction_soil)
    Trees[s,]$c_total <- Trees[s,]$c_ag + Trees[s,]$c_bg
    if (year > min(Years)) {
      Trees[s,]$anpp <- Trees[s,]$c_total - c_total_previous
    }
    # C:N stoichiometry; species
    if (Verbose) {
      message("Calculating N fraction...")
    }
    index <- StoichiometryLUT$species==species
    stoichiometry <- StoichiometryLUT[index,]
    Trees[s,]$n_ag <- Trees[s,]$c_ag /
      (biomass$fraction_stem   * stoichiometry$cn_stem)   +
      (biomass$fraction_branch * stoichiometry$cn_branch) +
      (biomass$fraction_leaf   * stoichiometry$cn_leaf)
    Trees[s,]$n_bg <- Trees[s,]$c_bg /
      (biomass$fraction_root * stoichiometry$cn_root) +
      (biomass$fraction_soil * stoichiometry$cn_soil)
    Trees[s,]$n_total <- Trees[s,]$n_ag + Trees[s,]$n_bg
    if (year > min(Years)) {
      Trees[s,]$annp <- Trees[s,]$n_total - n_total_previous
    }
  } # end species loop

  # Fix order of regeneration_df columns
  regeneration_df <- regeneration_df[,colnames(Trees)]
  
  # Calculate NEE and add to fluxes; NEE calculation based on Equation 1 in:
  #   Clark et al. (2001) Measuring Net Primary Production in Forests: Concepts and Field Methods,
  #   Ecological Applications, 11(2), pp. 356-370.
  # and
  #   Reichstein et al. (2005) On the separation of net ecosystem exchange into assimilation and
  #   ecosystem respiration: review and improved algorithm, Global Change Biology, 11, pp. 1424–1439.
  # Convert C to CO2 in ANPP
  fluxes_df <- data.frame(
    year   = year,
    nee    = (respiration_df$r_soil * Configuration$field_area) - sum(Trees$anpp, na.rm=TRUE),
    r_soil = (respiration_df$r_soil * Configuration$field_area)
  )

  # Save annual tree/cohort results
  write.csv(Trees, file=paste0("outputs/trees_", year, ".csv"), row.names=FALSE,
    fileEncoding="UTF-8")

  # Append results to CSV
  header_boolean <- year == Years[1]

  suppressWarnings(write.table(som_df, file="outputs/som.csv",
    append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  suppressWarnings(write.table(fluxes_df, file="outputs/fluxes.csv",
    append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  suppressWarnings(write.table(zstar_df, file="outputs/zstar.csv",
    append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  suppressWarnings(write.table(regeneration_df, file="outputs/regeneration.csv",
    append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  # Mortality is non-deterministic
  if (nrow(mortality_df) > 0 | header_boolean) {
    suppressWarnings(write.table(mortality_df, file="outputs/mortality.csv",
      append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  }

  # Calculate loop time
  loop_time <- proc.time() - loop_start_time
  loop_time <- as.numeric(loop_time["elapsed"])
  message(paste(paste0("Year: ", year, ","), "Duration (sec):", round(loop_time, 2)))
  loop_times <- c(loop_times, loop_time)
} # end year loop

# Calculate total time
total_time <- proc.time() - start_time
total_time <- as.numeric(total_time["elapsed"])
message("=========================================================================")
message("Simulation complete")
message(paste("Years:", length(Years)))
message(paste("Duration (sec):", round(total_time, 2)))
message("=========================================================================")

# Export loop times to CSV
loop_times_df <- data.frame(year=Years, time=loop_times)
write.table(loop_times_df, file="outputs/loop_times.csv", col.names=TRUE,
  row.names=FALSE, sep=",", fileEncoding="UTF-8")

# End #
