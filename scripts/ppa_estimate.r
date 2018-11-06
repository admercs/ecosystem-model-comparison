#=================================================================================================#
# Functions to generate growth and mortality parameters from inventory data
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

# Clean R environment and check session
rm(list=ls())
graphics.off()
sessionInfo()

# Libraries #
library(stats)
library(jsonlite)
library(ggplot2)

# Functions #

# Estimate growth rate (DBH increment) from inventory data
# This function assumes a dataframe with the columns { tag, plot, year, species }
# If no plots exist, set: df$plot <- 1
growth_rate <- function(df, overstory=TRUE, threshold=10) {
  growth = data.frame()
  for (year in unique(df$year)) {
    for (plot in unique(df$plot[df$year==year])) {
      for (species in unique(df$species[df$year==year & df$plot==plot])) {
        if (year == min(df$year)) {
          growth = rbind(growth, data.frame(year=year, plot=plot, species=species, increment=NA))
        } else {
          if (overstory == TRUE) {
            dbh          = mean(df$dbh[df$year==year     & df$plot==plot & df$species==species & df$dbh > threshold], na.rm=TRUE)
            dbh_previous = mean(df$dbh[df$year==(year-1) & df$plot==plot & df$species==species & df$dbh > threshold], na.rm=TRUE)
          } else {
            dbh          = mean(df$dbh[df$year==year     & df$plot==plot & df$species==species & df$dbh <= threshold], na.rm=TRUE)
            dbh_previous = mean(df$dbh[df$year==(year-1) & df$plot==plot & df$species==species & df$dbh <= threshold], na.rm=TRUE)
          }
          increment = dbh - dbh_previous
          increment = ifelse(dbh < dbh_previous, NA, increment) # remove the effect of disturbance
          growth    = rbind(growth, data.frame(year=year, plot=plot, species=species, increment=increment))
        }
      }
    }
  }
  i_rate = aggregate(growth$increment, by=list(year=growth$year, species=growth$species), FUN=mean, na.rm=TRUE)
  colnames(i_rate) = c("year","species","increment")
  return(i_rate)
}


# Estimate mortality rate from inventory data
# This function assumes a dataframe with the columns { tag, plot, year, species }
# If no plots exist, set: df$plot <- 1
mortality_rate <- function(df, overstory=TRUE, threshold=10) {
  mortality = data.frame()
  for (year in unique(df$year)) {
    for (plot in unique(df$plot[df$year==year])) {
      for (species in unique(df$species[df$year==year & df$plot==plot])) {
        if (year == min(df$year)) {
          mortality = rbind(mortality, data.frame(year=year, plot=plot, species=species, mu=NA))
        } else {
          if (overstory == TRUE) {
            tags          = unique(df$tag[df$year==year     & df$plot==plot & df$species==species & df$dbh > threshold])
            tags_previous = unique(df$tag[df$year==(year-1) & df$plot==plot & df$species==species & df$dbh > threshold])
          } else {
            tags          = unique(df$tag[df==year          & df$plot==plot & df$species==species & df$dbh <= threshold])
            tags_previous = unique(df$tag[df$year==(year-1) & df$plot==plot & df$species==species & df$dbh <= threshold])
          }
          mu = 1 - (sum(tags %in% tags_previous) / length(tags_previous))
          mortality = rbind(mortality, data.frame(year=year, plot=plot, species=species, mu=mu))
        }
      }
    }
  }
  mu_rate = aggregate(mortality$mu, by=list(year=mortality$year, species=mortality$species), FUN=mean, na.rm=TRUE)
  colnames(mu_rate) = c("year","species","mu")
  return(mu_rate)
}

# Estimate recruitment rate from inventory data
# This function assumes a dataframe with the columns { tag, plot, year, species }
# If no plots exist, set: df$plot <- 1
recruitment_rate <- function(df) {
  recruitment = data.frame()
  for (year in unique(df$year)) {
    for (plot in unique(df$plot[df$year==year])) {
      for (species in unique(df$species[df$year==year & df$plot==plot])) {
        if (year == min(df$year)) {
          recruitment = rbind(recruitment, data.frame(year=year, plot=plot, species=species, f=NA))
        } else {
          tags          = unique(df$tag[df==year          & df$plot==plot & df$species==species])
          tags_previous = unique(df$tag[df$year==(year-1) & df$plot==plot & df$species==species])
          f = sum(!tags %in% tags_previous)
          recruitment = rbind(recruitment, data.frame(year=year, plot=plot, species=species, f=f))
        }
      }
    }
  }
  f_rate = aggregate(recruitment$f, by=list(year=recruitment$year, species=recruitment$species), FUN=sum, na.rm=TRUE)
  colnames(f_rate) = c("year","species","f")
  return(f_rate)
}

# All functions assume the following columns are present:
# { tag, plot, year, species }

# HF EMS; Site: EMS
DataDir <- ""
setwd(DataDir)
trees <- read.csv(file.path(DataDir, "hf069-09-trees.csv"), header=TRUE, stringsAsFactors=FALSE)
colnames(trees) <- tolower(colnames(trees))
trees <- trees[trees$site=="ems",]
trees$species <- species_lut$code[match(trees$species, species_lut$abbr)] # use species LUT for HF
trees <- trees[trees$tree.type == "live" & !trees$species %in% c("ILVE","") & !is.na(trees$species) &
                 !is.na(trees$year) & trees$year > 1993,]
head(trees, 3)
sapply(trees, class)

# JERC RD; Sites: 3, 3100, 3200, and 3300
DataDir <- ""
list.files(DataDir)
trees <- read.csv(file.path(DataDir, "JC003-04-anpp-tree.csv"), header=TRUE, stringsAsFactors=FALSE)
colnames(trees) <- tolower(colnames(trees))
trees <- subset(trees, site==3) # get all the data we can when calculating rates by masking this
head(trees, 3)

# Tree List #

# HF EMS
# Required columns: id, species, type, age, dbh, height, crown_r, crown_l, crown_, ba
initial <- trees[trees$year == 1998,] # first year after the gap
initial <- setNames(aggregate(initial$dbh, by=list(initial$species, initial$plot, initial$tag),
                              FUN=mean, na.rm=TRUE), c("species","plot","tag","dbh"))
initial$id <- 1:nrow(initial)
initial <- initial[,c("id","species","dbh")]
initial <- cbind(initial, data.frame(type=NA, age=NA, height=NA, crown_r=NA, crown_l=NA, crown_a=NA, ba=NA))
initial <- initial[,c("id","species","type","age","dbh","height","crown_r","crown_l","crown_a","ba")]
initial$type <- ifelse(initial$dbh > 5, "adult", "sapling")
head(initial, 3)
sort(unique(initial$species))
write.csv(initial, file.path("", "trees.csv"), row.names=FALSE)
# update species list in overall file
trees <- trees[trees$species %in% as.character(sort(unique(initial$species))),]

# JERC RD
# Required columns: id, species, type, age, dbh, height, crown_r, crown_l, crown_, ba
initial <- trees[trees$year == 2004,]
any(duplicated(initial[,c("year","site","plot","tag","species")]))
initial <- setNames(aggregate(trees$dbh, by=list(trees$year, trees$species, trees$plot, trees$tag),
                              FUN=mean, na.rm=TRUE), c("year","species","plot","tag","dbh"))
initial$id <- 1:nrow(initial)
initial <- initial[,c("id","species","dbh_mean","ba")]
initial <- cbind(initial, data.frame(type=NA, age=NA, height=NA, crown_r=NA, crown_l=NA, crown_a=NA))
initial <- initial[,c("id","species","type","age","dbh","height","crown_r","crown_l","crown_a","ba")]
initial$type <- ifelse(initial$dbh > 5, "adult", "sapling")
head(initial, 3)
sort(unique(initial$species))
write.csv(initial, file.path("", "trees.csv"), row.names=FALSE)
# update species list in overall file
trees <- trees[trees$species %in% as.character(sort(unique(initial$species))),]

# Growth #

# overstory
increment_overstory <- growth_rate(trees, overstory=TRUE, threshold=5)
ggplot2::ggplot(increment_overstory, aes(x=year, y=increment, color=species)) + geom_line() + theme_minimal()
setNames(aggregate(increment_overstory$increment, by=list(increment_overstory$species), mean, na.rm=TRUE),
         c("species","increment_overstory"))

# understory
increment_understory <- growth_rate(trees, overstory=FALSE, threshold=12)
ggplot2::ggplot(increment_understory, aes(x=year, y=increment, color=species)) + geom_line() + theme_minimal()
setNames(aggregate(increment_understory$increment, by=list(increment_understory$species), mean, na.rm=TRUE),
         c("species","increment_understory"))

# Mortality #

# overstory
mu_overstory <- mortality_rate(trees, overstory=TRUE, threshold=5)
ggplot2::ggplot(mu_overstory, aes(x=year, y=mu, color=species)) + geom_line() + theme_minimal()
setNames(aggregate(mu_overstory$mu, by=list(mu_overstory$species), mean, na.rm=TRUE),
         c("species","mu_overstory"))

# understory
mu_understory <- mortality_rate(trees, overstory=FALSE, threshold=10)
ggplot2::ggplot(mu_understory, aes(x=year, y=mu, color=species)) + geom_line() + theme_minimal()
setNames(aggregate(mu_understory$mu, by=list(mu_understory$species), mean, na.rm=TRUE),
         c("species","mu_understory"))

# Recruitment #

# understory, naturally
f_rate <- recruitment_rate(trees)
ggplot2::ggplot(f_rate, aes(x=year, y=f, color=species)) + geom_line() + theme_minimal()
setNames(aggregate(f_rate$f, by=list(f_rate$species), mean, na.rm=TRUE), c("species","f"))
