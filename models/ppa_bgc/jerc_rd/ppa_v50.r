#!/usr/bin/env Rscript
#----------------------------------------------------------------------------------------------------------------#
# Filename:
#   ppa_v50.R
# Description:
#   PPA-BGC model for R
#   R implementation of the Perfect Plasticity Approximation (PPA)
# Version:
#   5.0
# Authors:
#   Adam Erickson,   Washington State University
#   Nikolay Strigul, Washington State University
# Modified:
#   June 30, 2018
# Examples:
#   Rscript --vanilla ppa_v50.r
#   Rscript --vanilla ppa_v50.r --wd /folder/test
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
#----------------------------------------------------------------------------------------------------------------#
Version <- "5.0"

options(error=quote({dump.frames(to.file=TRUE); q()}))

message("")
message("Authors: Adam Erickson, Nikolay Strigul")
message(paste("Version:", Version, "\n"))
Sys.sleep(0)

outputs_files   <- dir("outputs", recursive=FALSE, full.names=TRUE)
outputs_deleted <- file.remove(outputs_files)

args <- commandArgs(trailingOnly=TRUE)

if (length(args) == 0) {
  message(paste0("Working directory: ", getwd()))
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
} else {
  stop("Error: Only 0 or 2 arguments may be passed.", call.=FALSE)
}

# Functions #

round_int <- function(x, step) { 
  return(step * round(x / step))
}

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

Initialize <- function(Trees, ...) {
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
  for (species in unique(Trees$species)) {
    if (!species %in% AllometryLUT$species) { stop("Species not found in allometry table") }
    s <- Trees$species==species
    for (type in unique(Trees[s,]$type)) {
      st <- Trees$species==species & Trees$type==type
      index <- AllometryLUT$species==species & AllometryLUT$type==type
      allometry <- AllometryLUT[index,]
      Trees[st,]$height  <- 1.35 + (30-1.35) * (1-exp(-(allometry$h_coeff * Trees[st,]$dbh)))
      Trees[st,]$height  <- ifelse(Trees[st,]$height > 30, 30, Trees[st,]$height)
      Trees[st,]$type    <- ifelse(Trees[st,]$dbh > 10, "adult", "sapling")
      Trees[st,]$crown_l <- allometry$cd  * Trees[st,]$dbh
      Trees[st,]$crown_r <- allometry$cr1 * Trees[st,]$dbh^allometry$cr2
      Trees[st,]$crown_a <- Trees[st,]$crown_r^2 * pi
      Trees[st,]$ba      <- Trees[st,]$dbh^2 * 0.00007854
    }
    if (!species %in% BiomassLUT$species) { stop("Species not found in biomass table") }
    index <- BiomassLUT$species==species
    biomass <- BiomassLUT[index,]
    if (CohortMode) {
      Trees[s,]$biomass_ag <- (biomass$b0 + biomass$b1 * log(Trees[s,]$dbh)) * Trees[s,]$n_trees
    } else {
      Trees[s,]$biomass_ag <- (biomass$b0 + biomass$b1 * log(Trees[s,]$dbh))
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
    carbon <- CarbonLUT
    Trees[s,]$c_ag <- Trees[s,]$biomass_ag *
      ((biomass$fraction_stem   * carbon$fraction_stem)   +
       (biomass$fraction_branch * carbon$fraction_branch) +
       (biomass$fraction_leaf   * carbon$fraction_leaf))
    Trees[s,]$c_bg <- Trees[s,]$biomass_ag *
      ((biomass$fraction_root * carbon$fraction_root) +
       (biomass$fraction_soil * carbon$fraction_soil))
    Trees[s,]$c_total <- Trees[s,]$c_ag + Trees[s,]$c_bg
    Trees[s,]$anpp    <- NA
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

ensure_lowercase <- function(Trees) {
  colnames(Trees) <- tolower(colnames(Trees))
  Trees$type <- tolower(Trees$type)
  return(Trees)
}

leap_year <- function(year) {
  return(ifelse((year %%4 == 0 & year %%100 != 0) | year %%400 == 0, TRUE, FALSE))
}

respiration_soil <- function(Ta, P) {
  if (Ta < -13.3) {
    Rs <- 0
  } else {
    if (Ta > 33.5) { Ta <- 33.5 }
    e  <- exp(1)
    Q  <- 0.05452
    F  <- 1.250
    K  <- 4.259
    Rs <- F * e^(Q * Ta) * (P / (K + P))
  }
  return(Rs)
}

soc_depth = function(order, depth_cm=100) {
  soc_table = data.frame(
    order     = c("All","Alfisols","Andisols","Aridisols","Entisols","Histosols","Inceptisols",
                  "Mollisols","Spodosols","Ultisols","Vertisols"),
    intercept = c(1.1795,1.1122,1.3837,0.2065,0.9300,1.6227,1.1631,1.0163,1.4262,1.1576,0.5145),
    slope     = c(-0.8228,-0.8330,-0.8425,-0.1300,-0.7207,-1.0109,-0.7331,-0.6214,-0.9801,-0.8867,-0.2427)
  )
  coeffs = as.numeric(soc_table[soc_table$order==order,][,c("intercept","slope")])
  soc    = sapply(seq(1, depth_cm, 1), function(x) 10^(coeffs[1] + coeffs[2] * log10(x)))
  soc    = sum(soc) / 10
  return(soc)
}

start_time <- proc.time()
loop_times <- c()

# Data #

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

CohortMode <- as.logical(Configuration$cohort_mode)
if (CohortMode == TRUE) {
  message("Generating cohorts...")
  Trees <- GenerateCohorts(Trees)
}

message("Initializing...")
Trees <- Initialize(Trees)

Trees <- ensure_lowercase(Trees)

StartYear <- as.numeric(Configuration$start_year) 
EndYear   <- as.numeric(Configuration$end_year)
Years     <- StartYear:EndYear # 200 years

UnderstoryAdult   <- as.numeric(Configuration$understory_adult)
UnderstorySapling <- as.numeric(Configuration$understory_sapling)

id_holder <- nrow(Trees)

message("Starting simulation...")
message("=========================================================================")
for (year in Years) {
  
  loop_start_time <- proc.time()

  Trees$year <- year

  Trees$age <- Trees$age + 1

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
    Rs         <- Rs * month_days / 1000
    Rs_month   <- c(Rs_month, Rs)
  }
  Rs_year        <- sum(Rs_month, na.rm=TRUE) * Configuration$field_area
  respiration_df <- data.frame(year=year, r_soil=Rs_year)

  c_so   <- soc_depth(order=Configuration$soil_order, depth_cm=100)
  n_so   <- c_so / mean(StoichiometryLUT$cn_soil, na.rm=TRUE)
  som_df <- data.frame(year=year, c_so=c_so, n_so=n_so)

  Trees <- Trees[order(Trees$height, decreasing=TRUE), ]

  if (CohortMode) {
    Trees$crown_a_cumsum <- cumsum(Trees$crown_a * Trees$n_trees)
  } else {
    Trees$crown_a_cumsum <- cumsum(Trees$crown_a)
  }

  index    <- which.min(abs(Trees$crown_a_cumsum - Configuration$field_area))
  zstar_df <- Trees[index,]

  mortality_df <- data.frame(matrix(ncol=nrow(Trees), nrow=0))
  colnames(mortality_df) <- colnames(Trees)
  regeneration_df <- data.frame(matrix(ncol=nrow(Trees), nrow=0))
  colnames(regeneration_df) <- colnames(Trees)

  for (species in unique(Trees$species)) {
    if (!species %in% GrowthLUT$species) {
      next
    }

    s <- Trees$species==species
    
    n_replicates <- length(which(s))
    index <- RegenerationLUT$species==species
    n_saplings <- RegenerationLUT[index,]$mean
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
      id_holder <- id_holder + 1
    } else if (n_replicates > 0) {
      n_trees <- nrow(trees)
      new_saplings <- do.call(rbind, replicate(n_replicates * n_saplings, new_sapling, simplify=FALSE))
      new_saplings$id <- id_holder:(id_holder + (n_replicates * n_saplings) - 1)
      Trees <- rbind(Trees, new_saplings)
      regeneration_df <- rbind(regeneration_df, new_saplings)
      id_holder <- id_holder + (n_replicates * n_saplings)
    }

    s <- Trees$species==species

    for (type in unique(Trees[s,]$type)) {

      st <- Trees$species==species & Trees$type==type

      index <- MortalityLUT$species==species & MortalityLUT$type==type
      mortality <- MortalityLUT[index,]
      random_uniform <- runif(1, min=0, max=1)
      if (random_uniform < mortality$probability) {
        kill_fraction <- runif(1, min=0, max=1)
        if (CohortMode) {
          kill_n <- round(nrow(Trees[st,]) * kill_fraction)
          if (kill_n > 0) {
            kill_index <- sample(which(st), kill_n, replace=FALSE)
            mortality_df <- rbind(mortality_df,  Trees[kill_index,])
            Trees <- Trees[-kill_index,]
        } else {
          kill_n <- round(nrow(Trees[st,]) * kill_fraction)
          if (kill_n > 0) {
            kill_index <- sample(which(st), kill_n, replace=FALSE)
            mortality_df <- rbind(mortality_df, Trees[kill_index,])
            Trees <- Trees[-kill_index,]
          }
        }
      }

      st <- Trees$species==species & Trees$type==type

      if (nrow(Trees[st,]) < 1) {
        next
      }

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

      index <- AllometryLUT$species==species & AllometryLUT$type==type
      allometry <- AllometryLUT[index,]
      Trees[st,]$height  <- 1.35 + (30-1.35) * (1-exp(-(allometry$h_coeff * Trees[st,]$dbh)))
      Trees[st,]$height  <- ifelse(Trees[st,]$height > 30, 30, Trees[st,]$height)
      Trees[st,]$type    <- ifelse(Trees[st,]$dbh > 10, "adult", "sapling")
      Trees[st,]$crown_l <- allometry$cd  * Trees[st,]$dbh
      Trees[st,]$crown_r <- allometry$cr1 * Trees[st,]$dbh^allometry$cr2
      Trees[st,]$crown_a <- Trees[st,]$crown_r^2 * pi
      Trees[st,]$ba      <- Trees[st,]$dbh^2 * 0.00007854
    }

    s <- Trees$species==species

    if (nrow(Trees[s,]) < 1) {
      next
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

    c_total_previous <- Trees[s,]$c_total
    n_total_previous <- Trees[s,]$n_total

    carbon <- CarbonLUT
    Trees[s,]$c_ag <- Trees[s,]$biomass_ag * 
      (biomass$fraction_stem   * carbon$fraction_stem) +
      (biomass$fraction_branch * carbon$fraction_branch) +
      (biomass$fraction_leaf   * carbon$fraction_leaf)
    Trees[s,]$c_bg <- Trees[s,]$biomass_ag *
      (biomass$fraction_root * carbon$fraction_root) +
      (biomass$fraction_soil * carbon$fraction_soil)
    Trees[s,]$c_total <- Trees[s,]$c_ag + Trees[s,]$c_bg
    if (year > Years[1]) {
      Trees[s,]$anpp <- Trees[s,]$c_total - c_total_previous
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
  }

  regeneration_df <- regeneration_df[,colnames(Trees)]
  
  fluxes_df <- data.frame(
    year   = year,
    nee    = sum(Trees$anpp, na.rm=TRUE) - respiration_df$r_soil,
    r_soil = respiration_df$r_soil
  )

  write.csv(Trees, file=paste0("outputs/trees_", year, ".csv"), row.names=FALSE,
    fileEncoding="UTF-8")

  header_boolean <- (year == Years[1])

  suppressWarnings(write.table(som_df, file="outputs/som.csv",
    append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  suppressWarnings(write.table(fluxes_df, file="outputs/fluxes.csv",
    append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  suppressWarnings(write.table(zstar_df, file="outputs/zstar.csv",
    append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  suppressWarnings(write.table(regeneration_df, file="outputs/regeneration.csv",
    append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  if (nrow(mortality_df) > 0 | header_boolean) {
    suppressWarnings(write.table(mortality_df, file="outputs/mortality.csv",
      append=TRUE, col.names=header_boolean, row.names=FALSE, sep=",", fileEncoding="UTF-8"))
  }

  loop_time <- proc.time() - loop_start_time
  loop_time <- as.numeric(loop_time["elapsed"])
  message(paste(paste0("Year: ", year, ","), "Duration (sec):", round(loop_time, 2)))
  loop_times <- c(loop_times, loop_time)
}

total_time <- proc.time() - start_time
total_time <- as.numeric(total_time["elapsed"])
message("=========================================================================")
message("Simulation complete")
message(paste("Years:", length(Years)))
message(paste("Duration (sec):", round(total_time, 2)))
message("=========================================================================")

loop_times_df <- data.frame(year=Years, time=loop_times)
write.table(loop_times_df, file="outputs/loop_times.csv", col.names=TRUE,
  row.names=FALSE, sep=",", fileEncoding="UTF-8")

# End #
