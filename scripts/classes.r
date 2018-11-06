#=================================================================================================#
# classes.r
# Model classes and methods
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

#' Model R6 class
#'
#' Generic wrapper object for forest biogeochemistry models.
#' @docType class
#' @author Adam Erickson, PhD, Washington State University, \email{adam.michael.erickson@@gmail.com}
#' @export
#' @return Object of \code{\link{R6::R6Class}} with methods for running forest biogeochemistry models.
#' @format \code{\link{R6::R6Class}} object.
#' @examples
#' sortie <- Model$new(list(name="SortiePPA", version="5.0", directory="/sortie", executable="Rscript", arguments="--vanilla sortie_ppa_50.r"))
#' sortie$parameters <- c(0.7, 0.1)
#' sortie$files <- list()
#' print(sortie$command)
#' sortie$run()
#'
Model <- R6::R6Class("Model",
  public = list(
    #' @section Public:
    #' Public variables and methods for class Model.
    #' @field name Name character. Defaults to character().
    #' @field version Version character. Defaults to character().
    #' @field directory Directory character. Defaults to character().
    #' @field arguments Named list for model executable. Defaults to list().
    #' @field command Command character to run model. Defaults to character().
    #' @field files File input/output list. Defaults to list(inputs=list(), outputs=list(), observations=list()).
    #' @field results List for storing model results and parameters used. Defaults to list().
    #' Optional parameters
    #' @field timescale Numeric vector timescale parameter. Defaults to c().
    #' @field bounds Numeric vector bounding box or polygon coordinates. Defaults to c().
    #' Bounding box coordinates are generally of the form (Xmin, Xmax, Ymin, Ymax).
    name           = character(),
    version        = character(),
    directory      = character(),
    arguments      = list(),
    command        = character(),
    files          = list(inputs=list(), outputs=list(), observations=list()),
    predictions    = list(),
    observations   = list(),
    timescale      = c(),
    area           = c(),
    metrics        = c("nee","biomass_ag","c_ag","n_ag","biomass_bg","c_bg","n_bg","r_soil","anpp",
                       "c_so","n_so","biomass_sp","abundance_sp"),
    fitness        = data.frame(),
    # Optional helper parameters
    site           = character(),
    bounds         = c(),
    species        = list(),
    #' Initialize method
    #' This function is executed by default upon class initialization.
    #' @param directory Model directory path character.
    #' @return Model object.
    #' @example model$new()
    initialize = function(directory=NA) {
      if(is.na(directory)) {
        stop("Please provide a model directory path")
      } else if(!dir.exists(directory)) {
        stop("Model directory does not exist")
      } else {
        self$directory     = directory
        self$command       = trimws(paste(unlist(self$arguments), collapse=" "))
        self$loss_function = eval
        message(paste("Model initialized:", self$name, "version", self$version, "at", self$directory))
        message(paste("Using run command:", self$command))
      }
    },
    #' Run method
    #' @param ... Only use class variables.
    #' @return Command line output character.
    #' @example model$run()
    run = function(...) {
      if (!dir.exists(self$directory)) { stop("Model directory does not exist") }
      current_wd = getwd()
      setwd(self$directory)
      message(paste("Running model:", self$name))
      ptm    = proc.time()
      out    = system(self$command, intern=FALSE, wait=TRUE, timeout=0, ...)
      timing = proc.time() - ptm
      message(paste("Elapsed (sec):", round(timing["elapsed"], 2)))
      setwd(current_wd)
      stopifnot(out == 0)
      message("Model run complete")
      return(invisible(list(out=out, timing=timing)))
    },
    #' Update command method
    #' @param ... Only use class variables.
    #' @return None.
    #' @example model$update()
    update_command = function(...) {
      self$command = trimws(paste(unlist(self$arguments), collapse=" "))
      message(paste("Run command:", self$command))
      return(invisible(0))
    },
    #' Create backup of model input parameter files
    #' @note Required by optimize() and montecarlo().
    #' @param overwrite Boolean for whether to overwrite existing files.
    archive_parameters = function(overwrite=FALSE, ...) {
      archive_dir = file.path(self$directory, "archive")
      if (!dir.exists(archive_dir)) {
        message("Creating `archive` directory")
        dir.create(archive_dir)
        existing_files = NULL
      } else {
        existing_files = list.files(archive_dir, full.names=TRUE)
      }
      if (length(existing_files) > 0 & overwrite == TRUE) {
        message("Removing existing archive")
        file.remove(dir(archive_dir, all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE))
      } else if (length(existing_files) > 0 & overwrite == FALSE) {
        stop("Archive already exists. Specify `overwrite=TRUE` to replace existing files.")
      }
      message(paste("Copying parameter files to:"), archive_dir)
      input_files = as.character(unlist(self$files$inputs))
      file.copy(from=input_files, to=archive_dir, overwrite=overwrite, recursive=TRUE,
                copy.mode=TRUE, copy.date=TRUE)
      message("Parameter files successfully archived")
      return(invisible(0))
    },
    #' Update parameters
    #' @note Required by optimize() and montecarlo().
    #' @note Assumes CSV file format, override otherwise.
    #' @param ... Only use class variables.
    update_parameters = function(...) {
      for (i in seq(length(self$parameters))) {
        list_id    = names(self$parameters[i])
        list_names = unlist(strsplit(list_id, split="\\."))
        parameters = self$parameters[[list_names]]
        inputfile  = read.csv(self$files$inputs[[list_names]], header=TRUE, stringsAsFactors=FALSE)
        inputfile[names(parameters)] = unlist(parameters)
        write.csv(inputfile, self$files$inputs[[list_names]], row.names=FALSE)
      }
      return(invisible(0))
    },
    #' Optimization/Monte Carlo model prediction function placeholder.
    #' @param ... Only use class variables.
    #' @note Required by optimize().
    #' @note Pass user-defined function that computes variable of interest in temporal range from simulation outputs.
    get_predictions = function(...) {},
    #' Placeholder for function to update predictions stored in object.
    #' @param ... Only use class variables.
    #' @note Pass user-defined function that computes variable of interest in temporal range from simulation outputs.
    update_predictions = function(...) {},
    #' Optimization/Monte Carlo observation metric retrieval function.
    #' @param ... Only use class variables.
    #' @note Required by optimize(). Depends on set_parameters() function.
    get_observations = NULL,
    #' Set observations method; a hack around R6 class environment restrictions.
    #' @param observation_list List of file paths to observation files.
    #' @param get_method Function for retrieving metrics from observation list files.
    #' @note Required by optimize(). Depends on set_parameters() function.
    set_observations = function(observations, method, ...) {
      self$files$observations = observations
      environment(method)     = parent.env(environment())
      self$get_observations   = method
      self$observations       = method()
      message("Observations set")
    },
    #' Plot model predictions against observations for most recent model run
    #' @param ... Only use class variables.
    #' @note Requires run() be called and, get_predictions() and get_observations() be defined.
    #' @example model$plot()
    plot = function(variable=NA, limits=list(xlim=NA, ylim=NA), labels=list(xlab=NA, ylab=NA), ...) {
      if (length(self$observations) < 1) { self$observations = self$get_observations() }
      if (length(self$predictions)  < 1) { self$predictions  = self$get_predictions()  }
      if (anyNA(unlist(labels))) {
        if (grepl("^(biomass|abundance)_sp$", variable)) {
          ylab = tools::toTitleCase(gsub("^(biomass|abundance)_sp$", "\\1 Species", variable))
        } else if (grepl("^r_soil$", variable)) {
          ylab = gsub("^r_soil$", "Respiration Soil", variable)
        } else if (grepl("^biomass_(ag|bg)", variable)) {
          ylab = paste("Biomass", gsub("^(biomass_)|(ag|bg)$", "\\U\\2", variable, perl=TRUE))
        } else {
          ylab = gsub("(^(anpp|nee)$|^(c|n)_(ag|bg)$|^(c|n)_so$)", "\\U\\1", variable, perl=TRUE)
          ylab = gsub("_", " ", ylab)
        }
        labels = list(xlab="Year", ylab=ylab)
      }
      x = self$observations[[variable]]
      y = self$predictions[[variable]]
      years = Reduce(intersect, list(self$timescale, x$year, y$year))
      x = x[x$year %in% years,]
      x = x[order(x$year),]
      y = y[y$year %in% years,]
      y = y[order(y$year),]
      if (grepl("^(biomass|abundance)_sp$", variable)) {
        y_var = colnames(y)[grep("biomass|abundance", colnames(y))]
        p1 = ggplot2::ggplot(data=x, ggplot2::aes_string(x="year", y=y_var, color="species", label="species")) +
          ggplot2::labs(x=labels$xlab, y=labels$ylab) +
          ggplot2::guides(colour=ggplot2::guide_legend(override.aes=list(size=4))) +
          ggplot2::theme(legend.key=ggplot2::element_rect(fill=NA)) +
          ggplot2::geom_line() + ggplot2::theme_minimal()  +
          ggplot2::geom_text(data=x[x$year == max(years),],
                             ggplot2::aes_string(label="species", x="year + 0.5", y=y_var, color="species"))
        p2 = ggplot2::ggplot(data=y, ggplot2::aes_string(x="year", y=y_var, color="species", label="species")) +
          ggplot2::labs(x=labels$xlab, y="") +
          ggplot2::guides(colour=ggplot2::guide_legend(override.aes=list(size=4))) +
          ggplot2::theme(legend.key=ggplot2::element_rect(fill=NA)) +
          ggplot2::geom_line() + ggplot2::theme_minimal()  +
          ggplot2::geom_text(data=y[y$year == max(years),],
                             ggplot2::aes_string(label="species", x="year + 0.5", y=y_var, color="species"))
        print(ggpubr::ggarrange(p1, p2, widths=c(5,5), ncol=2, nrow=1, align="h", legend="right", common.legend=TRUE))
      } else {
        r2_model   = r_squared(x=x[,2], y=y[,2])
        rmse_model = rmse(x=x[,2], y=y[,2])
        ylim = c(min(x[,2], y[,2]), max(x[,2], y[,2]))
        if (anyNA(unlist(limits))) {
          plot(x=years, y=x[,2], type="l", xlab=labels$xlab, ylab=labels$ylab, ylim=ylim, col="blue", lwd=2)
          lines(x=years, y=y[,2], col="red", lwd=2)
        } else {
          plot(x=years, y=x[,2], type="l", xlab=labels$xlab, ylab=labels$ylab, xlim=limits$xlim,
               ylim=limits$ylim, col="blue", lwd=2)
          lines(x=years, y=y[,2], col="red", lwd=2)
        }
        #abline(0, 1, col="grey")
        mtext(paste0("R2: ", round(r2_model, 2), " RMSE: ", round(rmse_model, 2)), line=-1)
      }
      return(invisible(0))
    },
    #' Helper function to automatically produce PDF files of plots
    #' @param ... Only use class variables.
    save_pdf = function(directory, variable, width=5, height=5, family="Helvetica", paper="special", ...) {
      file = paste0(self$name, "_", self$site, "_", variable, ".pdf")
      pdf(file=file.path(directory, file), width=width, height=height, family=family, paper=paper)
      self$plot(variable=variable)
      graphics.off()
      message(paste("PDF file location:", file.path(directory, file)))
      return(invisible(0))
    },
    #' Helper function to generate a model fitness table
    #' @note Requires self$metrics, self$timescale, self$observations, self$predictions
    #' @param ... Only use class variables.
    get_fitness = function(...) {
      if (length(self$metrics) < 1) { stop("Please define the model metrics") }
      if (length(self$timescale) < 1) { stop("Please define the model timescale") }
      df = data.frame(metric=NULL, r2=NULL, rmse=NULL, mae=NULL)
      for (metric in self$metrics) {
        years   = Reduce(intersect, list(self$timescale,
                                         self$observations[[metric]]$year,
                                         self$predictions [[metric]]$year))
        if (any(grepl("_so", metric))) {
          x  = mean(self$observations[[metric]][,2], na.rm=TRUE)
          y  = mean(self$predictions [[metric]][,2], na.rm=TRUE)
          df = rbind(df, data.frame(metric=metric, r2=NA, rmse=rmse(x=x, y=y), mae=mae(x=x, y=y)))
        } else if (grepl("_sp", metric)) {
          x  = self$observations[[metric]][self$observations[[metric]]$year %in% years,]
          y  = self$predictions[[metric]] [self$predictions [[metric]]$year %in% years,]
          species_xy = intersect(x$species, y$species)
          #r2_scaling = length(species_xy) / length(x$species)
          if (length(species_xy) < 1) { stop("Species names do not match") }
          intermediate = data.frame(x=NULL, y=NULL)
          for (year in years) {
            for (species in species_xy) {
              x_i = x[x$year==year & x$species==species,][,3]
              y_i = y[y$year==year & y$species==species,][,3]
              if (length(x_i) != length(y_i)) {
                x_i = NA
                y_i = NA
              }
              intermediate = rbind(intermediate, data.frame(x=x_i, y=y_i))
            }
          }
          df = rbind(df, data.frame(
            metric = metric,
            r2     = r_squared(x=intermediate$x, y=intermediate$y),
            rmse   = rmse(x=intermediate$x, y=intermediate$y),
            mae    = mae(x=intermediate$x, y=intermediate$y)
            ))
        }
        else {
          x  = self$observations[[metric]][self$observations[[metric]]$year %in% years,][,2]
          y  = self$predictions[[metric]] [self$predictions [[metric]]$year %in% years,][,2]
          df = rbind(df, data.frame(metric=metric, r2=r_squared(x=x, y=y), rmse=rmse(x=x, y=y),
                                    mae=mae(x=x, y=y)))
        }
      }
      self$fitness = df
      return(self$fitness)
    }
  )
)


#' SortiePPA R6 class
#'
#' SortiePPA R6 subclass of Model R6 class.
#' @docType class
#' @author Adam Erickson, PhD, Washington State University, \email{adam.michael.erickson@@gmail.com}
#' @export
#' @return Object of \code{\link{R6::R6Class}} with methods for running forest biogeochemistry models.
#' @format \code{\link{R6::R6Class}} object.
#'
SortiePPA <- R6::R6Class("SortiePPA",
  inherit = Model,
  public = list(
    #' @section Public:
    #' Public variables and methods for class SortiePPA.
    #' Fields are inherited from Model class.
    name    = "Sortie-PPA",
    version = 5.0,
    #' Override default initialization
    initialize = function(directory=NA) {
      self$arguments = list(
        program = "Rscript",
        clean   = "--vanilla",
        file    = "ppa_v50.r",
        verbose = "--verbose"
      )
      self$files$inputs = list(
        configuration = file.path(directory, "configuration.csv"),
        climate       = file.path(directory, "climate.csv"),
        trees         = file.path(directory, "trees.csv"),
        lut = list(
          allometry     = file.path(directory, "lut", "allometry.csv"),
          biomass       = file.path(directory, "lut", "biomass.csv"),
          carbon        = file.path(directory, "lut", "carbon.csv"),
          growth        = file.path(directory, "lut", "growth.csv"),
          mortality     = file.path(directory, "lut", "mortality.csv"),
          regeneration  = file.path(directory, "lut", "regeneration.csv"),
          stoichiometry = file.path(directory, "lut", "stoichiometry.csv")
        )
      )
      self$files$outputs = list(
        fluxes       = file.path(directory, "outputs", "fluxes.csv"),
        trees        = file.path(directory, "outputs", "..."),
        mortality    = file.path(directory, "outputs", "mortality.csv"),
        regeneration = file.path(directory, "outputs", "regeneration.csv"),
        som          = file.path(directory, "outputs", "som.csv"),
        zstar        = file.path(directory, "outputs", "zstar.csv")
      )
      self$files$outputs$trees = get_files(file.path(self$directory, "outputs"), "trees", "csv")
      super$initialize(directory)
    },
    #' Override default run method to list output files
    #' @param ... Only use class variables.
    run = function(...) {
      out = super$run()
      self$files$outputs$trees = get_files(file.path(self$directory, "outputs"), "trees", "csv")
      self$update_predictions()
      return(invisible(out))
    },
    #' Override default get_prediction methods
    #' Variables and units:
    #' nee ................. kg C m-2 year-1
    #' biomass_ag .......... kg m-2
    #' c_ag ................ kg C m-2
    #' n_ag ................ kg N m-2
    #' biomass_bg .......... kg m-2
    #' c_bg ................ kg C m-2
    #' n_bg ................ kg N m-2
    #' c_so ................ kg C m-2
    #' n_so ................ kg N m-2
    #' r_soil .............. kg C m-2 year-1
    #' anpp ................ kg m-2 year-1
    #' biomass_sp .......... kg m-2
    #' abundance_sp ........ % of total n-trees
    #' @param variable Character for metric to select.
    #' @param ... Only class variables.
    get_predictions = function(variable="all", ...) {
      if (length(self$timescale) < 2) { stop("Please set a timescale before getting predictions") }
      #Fluxes
      fluxes  = read.csv(self$files$outputs$fluxes, header=TRUE, stringsAsFactors=FALSE)
      fluxes$nee    = fluxes$nee / self$area
      fluxes$r_soil = fluxes$r_soil / self$area
      # Trees
      trees   = do.call(rbind, lapply(self$files$outputs$trees, read.csv, header=TRUE, stringsAsFactors=FALSE))
      columns = c("biomass_ag","c_ag","n_ag","biomass_bg","c_bg","n_bg","anpp")
      annual  = aggregate(trees[,columns], list(trees$year), sum, na.rm=TRUE)
      colnames(annual)[1] = "year"
      annual[,columns] = annual[,columns] / self$area # adjust for site area
      # SOM
      som = read.csv(self$files$outputs$som, header=TRUE, stringsAsFactors=FALSE)
      if (!"n_trees" %in% colnames(trees)) { trees$n_trees <- 1 }
      # Species biomass
      biomass_sp = setNames(aggregate(trees$biomass_total, list(trees$year, trees$species), sum, na.rm=TRUE),
                                 c("year","species","biomass_ag"))
      biomass_sp$biomass_ag = biomass_sp$biomass_ag / self$area # adjust for site area
      # Species relative abundance
      abundance_sp = setNames(aggregate(trees$n_trees, list(trees$year, trees$species), sum, na.rm=TRUE),
                                   c("year","species","n_trees"))
      abundance_total = setNames(aggregate(trees$n_trees, list(trees$year), sum, na.rm=TRUE), c("year","n_trees"))
      abundance_sp$total     = abundance_total$n_trees[match(abundance_sp$year, abundance_total$year)]
      abundance_sp$abundance = abundance_sp$n_trees / abundance_sp$total
      switch(variable,
        nee          = { fluxes[,c("year","nee")] },
        biomass_ag   = { annual[,c("year","biomass_ag")] },
        c_ag         = { annual[,c("year","c_ag")] },
        n_ag         = { annual[,c("year","n_ag")] },
        biomass_bg   = { annual[,c("year","biomass_bg")] },
        c_bg         = { annual[,c("year","c_bg")] },
        n_bg         = { annual[,c("year","n_bg")] },
        c_so         = { som[,c("year","c_so")] },
        n_so         = { som[,c("year","n_so")] },
        r_soil       = { fluxes[,c("year","r_soil")] },
        anpp         = { annual[,c("year","anpp")] },
        biomass_sp   = { biomass_sp },
        abundance_sp = { abundance_sp[,c("year","species","abundance")] },
        all = {
          list(
            nee          = fluxes[,c("year","nee")],
            biomass_ag   = annual[,c("year","biomass_ag")],
            c_ag         = annual[,c("year","c_ag")],
            n_ag         = annual[,c("year","n_ag")],
            biomass_bg   = annual[,c("year","biomass_bg")],
            c_bg         = annual[,c("year","c_bg")],
            n_bg         = annual[,c("year","n_bg")],
            c_so         = som[,c("year","c_so")],
            n_so         = som[,c("year","n_so")],
            r_soil       = fluxes[,c("year","r_soil")],
            anpp         = annual[,c("year","anpp")],
            biomass_sp   = biomass_sp,
            abundance_sp = abundance_sp[,c("year","species","abundance")]
          )
        },
        stop("Please specify a variable name or 'all'")
      )
    },
    #' Update predictions stored in object
    #' @param ... Only use class variables.
    update_predictions = function(...) {
      if (length(self$timescale) < 1) { stop("Output files must exist and timescale must be set") }
      if (length(self$files$outputs$trees) < 1) {
        self$files$outputs$trees = get_files(file.path(self$directory, "outputs"), "trees", "csv")
      }
      self$predictions = self$get_predictions()
      return(invisible(0))
    }
  )
)


#' Landis2 R6 class
#'
#' Landis2 R6 subclass of Model R6 class.
#' @docType class
#' @author Adam Erickson, PhD, Washington State University, \email{adam.michael.erickson@@gmail.com}
#' @export
#' @return Object of \code{\link{R6::R6Class}} with methods for running forest biogeochemistry models.
#' @format \code{\link{R6::R6Class}} object.
#'
Landis2 <- R6::R6Class("Landis2",
  inherit = Model,
  public = list(
    #' @section Public:
    #' Public variables and methods for class LandisII.
    #' Fields are inherited from Model class.
    name    = "Landis2",
    version = 6.2,
    #' Override default initialization
    initialize = function(directory=NA) {
      self$arguments = list(
        program  = "landis-ii",
        scenario = "scenario.txt"
      )
      self$files$inputs = list(
        core = list(
          scenario   = file.path(directory, "scenario.txt"),
          species    = file.path(directory, "species.txt"),
          ecoregions = file.path(directory, "ecoregions.txt")
        ),
        succession = list(
          necn = list(
            main       = file.path(directory, "succession", "necn", "main.txt"),
            ic         = file.path(directory, "succession", "necn", "ic.txt"),
            reductions = file.path(directory, "succession", "necn", "reductions.txt"),
            rasters = list(
              ic             = file.path(directory, "succession", "necn", "rasters", "ic.img"),
              soil_depth     = file.path(directory, "succession", "necn", "rasters", "soil_depth.img"),
              soil_drainage  = file.path(directory, "succession", "necn", "rasters", "soil_drainage.img"),
              flow_base      = file.path(directory, "succession", "necn", "rasters", "flow_base.img"),
              flow_storm     = file.path(directory, "succession", "necn", "rasters", "flow_storm.img"),
              field_capacity = file.path(directory, "succession", "necn", "rasters", "field_capacity.img"),
              wilting_point  = file.path(directory, "succession", "necn", "rasters", "wilting_point.img"),
              percent_sand   = file.path(directory, "succession", "necn", "rasters", "percent_sand.img"),
              percent_clay   = file.path(directory, "succession", "necn", "rasters", "percent_clay.img"),
              som1c_surface  = file.path(directory, "succession", "necn", "rasters", "som1_c_surface.img"),
              som1n_surface  = file.path(directory, "succession", "necn", "rasters", "som1_n_surface.img"),
              som1c_soil     = file.path(directory, "succession", "necn", "rasters", "som1_c_soil.img"),
              som1n_soil     = file.path(directory, "succession", "necn", "rasters", "som1_n_soil.img"),
              som2c          = file.path(directory, "succession", "necn", "rasters", "som2_c.img"),
              som2n          = file.path(directory, "succession", "necn", "rasters", "som2_n.img"),
              som3c          = file.path(directory, "succession", "necn", "rasters", "som3_c.img"),
              som3n          = file.path(directory, "succession", "necn", "rasters", "som3_n.img"),
              dead_wood      = file.path(directory, "succession", "necn", "rasters", "dead_wood.img"),
              dead_roots     = file.path(directory, "succession", "necn", "rasters", "dead_roots.img")
            ),
            climate = list(
              generator = file.path(directory, "succession", "necn", "climate", "generator.txt"),
              data      = file.path(directory, "succession", "necn", "climate", "data.txt"),
              spinup    = data
            )
          ),
          age_only = list(),
          biomass  = list(),
          pnet     = list()
        ),
        disturbance = list(
          base_fire    = list(),
          base_wind    = list(),
          base_harvest = list(),
          scrapple     = list()
        ),
        statistics = list(
          cohort_statistics = list(),
          reclassification  = list(),
          biomass           = list(),
          biomass_community = list()
        )
      )
      self$files$outputs = list(
        core = list(),
        succession = list(
          necn = list(
            logs = list(
              # Time (simulation year), ClimateRegionName, NEEC, AGB, AG_NPPC, BG_NPPC, TotalN, TotalNdep
              annual    = file.path(directory, "outputs", "succession", "necn", "logs", "NECN-succession-log.csv"),
              # Time (simulation year), Month, ClimateRegionName, avgNPPtc, avgResp, avgNEE, Ndep
              monthly   = file.path(directory, "outputs", "succession", "necn", "logs", "NECN-succession-monthly-log.csv"),
              # Year (simulation year), Month, SpeciesName, CohortAge
              calibrate = file.path(directory, "outputs", "succession", "necn", "logs", "NECN-calibrate-log.csv")
            ),
            rasters = list(
              anpp          = file.path(directory, "outputs", "succession", "necn", "rasters", "..."),
              nee           = file.path(directory, "outputs", "succession", "necn", "rasters", "..."),
              h2o_budget    = file.path(directory, "outputs", "succession", "necn", "rasters", "..."),
              h2o_available = file.path(directory, "outputs", "succession", "necn", "rasters", "..."),
              lai           = file.path(directory, "outputs", "succession", "necn", "rasters", "..."),
              n_soil        = file.path(directory, "outputs", "succession", "necn", "rasters", "..."),
              c_som         = file.path(directory, "outputs", "succession", "necn", "rasters", "..."),
              c_total       = file.path(directory, "outputs", "succession", "necn", "rasters", "...")
            )
          ),
          age_only = list(),
          biomass  = list(),
          pnet     = list()
        ),
        disturbance = list(
          base_fire    = list(),
          base_wind    = list(),
          base_harvest = list(),
          scrapple     = list()
        ),
        statistics = list(
          species_age       = list(),
          cohort_statistics = list(),
          reclassification  = list(),
          biomass           = list(),
          biomass_community = list()
        )
      )
      super$initialize(directory)
    },
    #' Set LANDIS_VERSION environmental variable for executable
    #' @note Necessary if multiple LANDIS-II versions are installed.
    #' @param ... Only use class variables.
    set_version = function(...) {
      Sys.putenv(LANDIS_VERSION=self$version)
      return(invisible(self$version))
    },
    #' Method to read species file to list from txt
    #' @note Stores results in parameters$species list.
    #' @param ... Only use class variables.
    read_species = function(...) {
      txt     = readLines(self$files$input$core$species)
      header  = read_index(txt, "LandisData", col.names=c("key","value"))
      name    = c("species","longevity_age","maturity_age","tolerance_shade","tolerance_fire",
                  "dispersal_distance_effective","dispersal_distance_max","regeneration_vegetative_p",
                  "resprouting_age_min","resprouting_age_max","regeneration_postfire")
      species = read_index(txt, "LandisData|>>", col.names=name, not=TRUE)
      # Write to object: self$parameters$core$species
      return(invisible(list(header=header, species=species)))
    },
    #' Method to write species file to txt
    #' @note Expects parameters$succession$necn was first populated by read_necn().
    #' @param overwrite Boolean switch for overwriting an existing file.
    write_species = function(overwrite=TRUE, ...) {
      path = file.path(self$directory, "succession", "necn", "species.txt")
      if (file.exists(path) & overwrite==TRUE) {
        file.remove(path)
      } else if (file.exists(path) & overwrite==FALSE) {
        stop("File exists. Set overwrite=TRUE (default) to replace existing file.")
      }
      df = self$parameters$core$species
      banner = "LandisData  Species"
      header = c("species","longevity_age","maturity_age","tolerance_shade","tolerance_fire",
                 "dispersal_distance_effective","dispersal_distance_max","regeneration_vegetative_p",
                 "resprouting_age_min","resprouting_age_max","regeneration_postfire")
      write(banner, path, append=TRUE)
      write(" ",    path, append=TRUE)
      write(header, path, append=TRUE)
      write(t(df),  path, ncolumns=ncol(df), append=TRUE, sep=" ")
      write(" ",    path, append=TRUE)
      return(invisible(0))
    },
    #' Method to read necn_hydro file to list from txt
    #' @note Designed for NECN version 5.
    #' @note Parses tables within text file to dataframes within list.
    #' @param ... Only use class variables.
    read_necn_main = function(...) {
      txt  = read_clean(self$files$input$succession$necn$main)
      # Extract key-value pairs
      kv   = extract_between(txt, before=FALSE, after="MaximumLAI", col.names=c("key","value"))
      # LAI maximum by shade class table
      lai  = extract_between(txt, before="MaximumLAI", after="LightEstablishmentTable",
                             col.names=c("shade_class","lai_max"))
      # Light-establishment table
      let  = extract_between(txt, before="LightEstablishmentTable", after="SpeciesParameters",
                            col.names=paste0("actual_shade_",0:5), row.names=paste0("shade_class_",1:5))
      # Species table
      name = c("species","pft","n_fixer","gdd_min","gdd_max","t_min_january","drought_max","leaf_longevity",
               "epicormic_resprout","lignin_leaf","lignin_root_fine","lignin_wood","lignin_root_coarse",
               "cn_leaf","cn_root_fine","cn_wood","cn_root_coarse","cn_litter","anpp_max","biomass_max")
      spp  = extract_between(txt, before="SpeciesParameters", after="FunctionalGroupParameters",
                             col.names=name)
      # Functional group table
      name = c("pft","index","ppdf1_t_mean","ppdf2_t_max","ppdf3_t_shape","ppdf4_t_shape","anpp_leaf",
               "biomass_to_lai","wood_c_lai_halfmax","lai_max","pprpts2","pprpts3","decay_rate_wood",
               "mortality_wood_monthly","mortality_shape","leaf_abscission_month","root_coarse","root_fine")
      pft  = extract_between(txt, before="FunctionalGroupParameters", after="FireReductionParameters",
                             col.names=name)
      # Fire reduction table
      name = c("severity_class","reduction_wood","reduction_litter","reduction_som")
      fire = extract_between(txt, before="FireReductionParameters", after="HarvestReductionParameters",
                             col.names=name)
      # Harvest reduction table
      name = c("name","reduction_wood","reduction_litter","reduction_som","removal_wood","removal_leaf")
      harvest = extract_between(txt, before="HarvestReductionParameters", after=FALSE, col.names=name)
      # Write to object: self$parameters$succession$necn$main
      return(invisible(list(kv=kv, lai=lai, establishment=let, species=spp,
                            pft=pft, fire=fire, harvest=harvest)))
    },
    #' Method to write necn_hydro file from list to txt
    #' @note Designed for NECN version 5.
    #' @note Expects parameters$succession$necn was first populated by read_necn().
    #' @param overwrite Boolean switch for overwriting an existing file.
    write_necn_main = function(overwrite=TRUE, ...) {
      path   = file.path(self$directory, "succession", "necn", "main.txt")
      if (file.exists(path) & overwrite==TRUE) {
        file.remove(path)
      } else if (file.exists(path) & overwrite==FALSE) {
        stop("File exists. Set overwrite=TRUE (default) to replace existing file.")
      }
      for (i in seq(length(self$parameters$succession$necn$main))) {
        df = self$parameters$succession$necn$main[[i]]
        header = ifelse(i==2, "MaximumLAI",
                 ifelse(i==3, "LightEstablishmentTable",
                 ifelse(i==4, "SpeciesParameters",
                 ifelse(i==5, "FunctionalGroupParameters",
                 ifelse(i==6, "FireReductionParameters",
                 ifelse(i==7, "HarvestReductionParameters",
                 ""))))))
        if (i > 1) { write(header, path, append=TRUE) }
        write(t(df), path, ncolumns=ncol(df), append=TRUE, sep=" ")
        write(" ", path, append=TRUE)
      }
      return(invisible(0))
    },
    #' Override default self$update_parameters() method
    #' @note Selectively updates parameters using a named list.
    #' @note Landis-2 requires custom read-write functions due to the txt format used.
    #' @note Required by optimize() and montecarlo().
    #' @param ... Only use class variables.
    update_parameters = function(...) {
      if (length(self$parameters$core$species) < 1 & length(self$parameters$succession$necn$main) < 1) {
        message("No parameters to update")
        return(invisible(1))
      }
      if (length(self$parameters$core$species) > 0) {
        species_file = self$read_species()
        parameters   = self$parameters$core$species
        species_file[names(parameters)] = unlist(parameters)
        self$parameters$core$species = species_file
        self$write_species()
      }
      if (length(self$parameters$succession$necn$main) > 0) {
        necn_file  = self$read_necn_main()
        parameters = self$parameters$succession$necn$main
        necn_file[names(parameters)] = unlist(parameters)
        self$parameters$succession$necn$main = necn_file
        self$write_necn_main()
      }
      return(invisible(0))
    },
    #' Override default run method to list output files
    #' @param ... Only use class variables.
    run = function(...) {
      out = super$run()
      message("Cleaning up")
      self$get_outputs()
      self$update_predictions()
      return(invisible(out))
    },
    #' Get and arrange outputs from NECN v5 succession model
    #' @param ... Only use class variables.
    get_outputs = function(...) {
      dir_move(from=file.path(self$directory, "NECN"),     to=file.path(self$directory, "outputs", "succession", "necn", "rasters"))
      dir_move(from=file.path(self$directory, "Metadata"), to=file.path(self$directory, "outputs", "succession", "necn", "metadata"))
      outputs = get_files(self$directory, c("NECN-.*?-log","Climate-.*?-log","Landis-.*?-log"), "*")
      files_move(files=outputs, to=file.path(self$directory, "outputs", "succession", "necn", "logs"))
      # Log files
      self$files$outputs$succession$necn$logs$annual    = file.path(self$directory, "outputs", "succession", "necn", "logs", "NECN-succession-log.csv")
      self$files$outputs$succession$necn$logs$monthly   = file.path(self$directory, "outputs", "succession", "necn", "logs", "NECN-succession-monthly-log.csv")
      self$files$outputs$succession$necn$logs$calibrate = file.path(self$directory, "outputs", "succession", "necn", "logs", "NECN-calibrate-log.csv")
      # Raster images
      self$files$outputs$succession$necn$rasters$anpp          = get_files(file.path(self$directory, "outputs", "succession", "necn", "rasters"), "AG_NPP",         "img")
      self$files$outputs$succession$necn$rasters$nee           = get_files(file.path(self$directory, "outputs", "succession", "necn", "rasters"), "ANEE",           "img")
      self$files$outputs$succession$necn$rasters$h2o_budget    = get_files(file.path(self$directory, "outputs", "succession", "necn", "rasters"), "water-budget",   "img")
      self$files$outputs$succession$necn$rasters$h2o_available = get_files(file.path(self$directory, "outputs", "succession", "necn", "rasters"), "AvailableWater", "img")
      self$files$outputs$succession$necn$rasters$lai           = get_files(file.path(self$directory, "outputs", "succession", "necn", "rasters"), "LAI",            "img")
      self$files$outputs$succession$necn$rasters$n_soil        = get_files(file.path(self$directory, "outputs", "succession", "necn", "rasters"), "SoilN",          "img")
      self$files$outputs$succession$necn$rasters$c_som         = get_files(file.path(self$directory, "outputs", "succession", "necn", "rasters"), "SOMTC",          "img")
      self$files$outputs$succession$necn$rasters$c_total       = get_files(file.path(self$directory, "outputs", "succession", "necn", "rasters"), "TotalC",         "img")
      return(invisible(0))
    },
    #' Override default get_prediction methods
    #' Variables and units:
    #' nee ................. kg C m-2 year-1
    #' biomass_ag .......... kg m-2
    #' c_ag ................ kg C m-2
    #' n_ag ................ kg N m-2
    #' biomass_bg .......... kg m-2
    #' c_bg ................ kg C m-2
    #' n_bg ................ kg N m-2
    #' c_so ................ kg C m-2
    #' n_so ................ kg N m-2
    #' r_soil .............. kg C m-2 year-1
    #' anpp ................ kg m-2 year-1
    #' biomass_sp .......... kg m-2
    #' abundance_sp ........ % of total n-trees
    #' @param variable Character for metric to select.
    #' @param ... Only class variables.
    get_predictions = function(variable="all", ...) {
      if (length(self$timescale) < 2) { stop("Please set a timescale before getting predictions") }
      # Annual log
      annual       = read.csv(self$files$outputs$succession$necn$logs$annual,  header=TRUE, stringsAsFactors=FALSE)
      annual$year  = annual$Time + min(self$timescale)
      annual       = annual[annual$year %in% min(self$timescale):max(self$timescale),]
      # Monthly log
      monthly      = read.csv(self$files$outputs$succession$necn$logs$monthly, header=TRUE, stringsAsFactors=FALSE)
      monthly$year = monthly$Time + min(self$timescale)
      monthly      = monthly[monthly$year %in% min(self$timescale):max(self$timescale),]
      variables    = c("avgNPPtc","avgResp","avgNEE") # only avgResp used
      annualized   = aggregate(monthly[,variables], by=list(monthly$year), FUN=mean, na.rm=TRUE)
      colnames(annualized)[1] <- "year"
      # Calibration log
      calibrate    = read.csv(self$files$outputs$succession$necn$logs$calibrate, header=TRUE, stringsAsFactors=FALSE)
      calibrate$year = calibrate$Year + min(self$timescale)
      calibrate    = calibrate[calibrate$year %in% min(self$timescale):max(self$timescale),]
      # Species abundance
      abundance_sp = setNames(aggregate(calibrate$SpeciesName,
                                             by=list(calibrate$year, calibrate$SpeciesName), length),
                                   c("year","species","n_trees"))
      abundance_total = setNames(aggregate(calibrate$SpeciesName, by=list(calibrate$year), length),
                                   c("year","n_trees"))
      abundance_sp$total     = abundance_total$n_trees[match(abundance_sp$year, abundance_total$year)]
      abundance_sp$abundance = abundance_sp$n_trees / abundance_sp$total
      # Species biomass
      biomass_sp = setNames(aggregate(calibrate$Bcohort,
                                      by=list(calibrate$year, calibrate$SpeciesName),
                                      sum, na.rm=TRUE), c("year","species","biomass_ag"))
      biomass_sp$biomass_ag = biomass_sp$biomass_ag / 1000 / self$area
      nee        = data.frame(year=annual$year, nee        = annual$NEEC / 1000)
      biomass_ag = data.frame(year=annual$year, biomass_ag = annual$AGB / 1000)
      c_ag       = data.frame(year=annual$year, c_ag       = annual$AGB * 0.5 / 1000)
      n_ag       = data.frame(year=annual$year, n_ag       = (annual$N_Leaf + annual$N_Wood) / 1000)
      biomass_bg = data.frame(year=annual$year, biomass_bg = (((annual$C_LiveFRoot + annual$C_LiveCRoot) * 2) + (annual$AGB * 0.68)) / 1000)
      c_bg       = data.frame(year=annual$year, c_bg       = (annual$C_LiveFRoot + annual$C_LiveCRoot) / 1000)
      n_bg       = data.frame(year=annual$year, n_bg       = (annual$N_FRoot + annual$N_CRoot) / 1000)
      c_so       = data.frame(year=annual$year, c_so       = annual$SOMTC / 1000)
      n_so       = data.frame(year=annual$year, n_so       = (annual$N_SOM1surf + annual$N_SOM1soil +
                                                             annual$N_SOM2 + annual$N_SOM3) / 1000)
      r_soil     = data.frame(year=annualized$year, r_soil = annualized$avgResp / 1000)
      anpp       = data.frame(year=annual$year, anpp       = (annual$AG_NPPC + annual$BG_NPPC) * 2 / 1000)
      switch(variable,
        nee          = { nee },
        biomass_ag   = { biomass_ag },
        c_ag         = { c_ag },
        n_ag         = { n_ag },
        biomass_bg   = { biomass_bg },
        c_bg         = { c_bg },
        n_bg         = { n_bg },
        c_so         = { c_so },
        n_so         = { n_so },
        r_soil       = { r_soil },
        anpp         = { anpp },
        biomass_sp   = { biomass_sp },
        abundance_sp = { abundance_sp[,c("year","species","abundance")] },
        all = {
          list(
            nee          = nee,
            biomass_ag   = biomass_ag,
            c_ag         = c_ag,
            n_ag         = n_ag,
            biomass_bg   = biomass_bg,
            c_bg         = c_bg,
            n_bg         = n_bg,
            c_so         = c_so,
            n_so         = n_so,
            r_soil       = r_soil,
            anpp         = anpp,
            biomass_sp   = biomass_sp,
            abundance_sp = abundance_sp[,c("year","species","abundance")]
          )
        },
        stop("Please specify a variable name or 'all'")
      )
    },
    #' Update predictions stored in object
    #' @param ... Only use class variables.
    update_predictions = function(...) {
      if (length(self$files$outputs) < 1 & length(self$timescale) < 1) {
        stop("Output files must exist and timescale must be set")
      }
      self$predictions = self$get_predictions()
      return(invisible(0))
    }
  )
)
