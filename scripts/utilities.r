#=================================================================================================#
# utilities.r
# Generic utility methods
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

#' Estimate empirical PDF and sample from CDF; bw = SD of density kernel; value range = 0:1.
#' @param x Numeric vector of empirical values.
#' @param n Numeric number of samples to draw.
#' @param from Numeric minimum sample value.
#' @param to Numeric maximum sample value.
#' @param bw Numeric smoothing bandwidth or standard deviation for Gaussian kernel.
#' @return Numeric vector.
sample_empirical_cdf = function(x, n, from=0, to=1, bw=0.1, ...) {
  x = x[!is.na(x) & is.finite(x)]
  if (length(x) < 1) { return(NA) }
  x_pdf = stats::density(x=x, kernel="gaussian", from=from, to=to, bw=bw, n=512)
  return(stats::approx(cumsum(x_pdf$y)/sum(x_pdf$y), x_pdf$x, stats::runif(n=n))$y)
}

#' Helper function to calculate mean from arbitrary inputs.
#' @param ... Numeric vector or multiple numeric arguments.
#' @return Numeric value.
mean0 = function(...) {
  x = as.numeric(unlist(list(...)))
  return(sum(x, na.rm=TRUE) / length(x))
}

#' Coefficient of determination (R2) helper method.
#' @note Different methods of calculating R2 produce slightly different results.
#' @note Similar to Nash-Sutcliffe efficiency coefficient in the general case.
#' @param x Numeric observation values.
#' @param y Numeric prediction values.
#' @return Numeric value.
r_squared = function(x, y) {
  if (length(x) != length(y)) { stop("Lengths of y and y_hat differ") }
  #return(1 - sum((x-y)^2) / sum((x-mean(x))^2))
  #return(cov(x, y, method="pearson") / (sd(x) * sd(y)))
  return(cor(x, y, method="pearson")^2)
}

#' Nash-Sutcliffe efficiency coefficient (NSE) helper method.
#' @note Implementation based on R2 in Nash and Sutcliffe (1970) River flow forecasting through conceptual models.
#' @note This metric is "analogous to the coefficent of determination."
#' @param x Numeric observation values.
#' @param y Numeric prediction values.
#' @return Numeric value.
nse = function(x, y) {
  if (length(x) != length(y)) { stop("Lengths of y and y_hat differ") }
  return((sum((x-mean(x))^2) - sum((x-y)^2)) / sum((x-mean(x))^2))
}

#' Root mean squared error (RMSE) helper method.
#' @param x Numeric observation values.
#' @param y Numeric prediction values.
#' @return Numeric value.
rmse = function(x, y) {
  if (length(x) != length(y)) { stop("Lengths of y and y_hat differ") }
  return(sqrt(sum((y-x)^2) / length(x)))
}

#' Mean absolute error (MAE) helper method.
#' @param x Numeric observation values.
#' @param y Numeric prediction values.
#' @return Numeric value.
mae = function(x, y) {
  if (length(x) != length(y)) { stop("Lengths of y and y_hat differ") }
  return(mean(abs(y-x)))
}

#' Mean error (ME) helper method.
#' @note Equivalent to bias.
#' @param x Numeric observation values.
#' @param y Numeric prediction values.
#' @return Numeric value.
me = function(x, y) {
  if (length(x) != length(y)) { stop("Lengths of y and y_hat differ") }
  return(mean(y-x))
}

#' Helper method for listing files.
#' @note List files in folder with extension matching pattern
#' @param folder Character folder path containing files.
#' @param filnames Character vector of filenames to search.
#' @param extension Character file extension to search.
#' @return Character vector of files matching pattern.
get_files = function(path, filenames, extension, ...) {
  pattern = paste0("^(", paste(filenames, collapse="|"), ").*?", ".", extension, "$")
  return(list.files(path=path, pattern=pattern, full.names=TRUE))
}

#' Helper method for moving multiple files.
#' @param files Character vector of file paths
#' @param to Character directory target
#' @param recursive Boolean recursive switch
#' @param overwrite Boolean overwrite switch
#' @return 0
files_move = function(files, to, recursive=TRUE, overwrite=TRUE) {
  if (!dir.exists(file.path(to))) { dir.create(path=file.path(to), recursive=recursive) }
  file.copy(from=files, to=file.path(to), recursive=recursive, overwrite=overwrite, copy.mode=TRUE)
  file.remove(files)
  return(invisible(0))
}

#' Helper method for moving directories.
#' @param from Character directory origin
#' @param to Character directory target
#' @param recursive Boolean recursive switch
#' @param overwrite Boolean overwrite switch
#' @return 0
dir_move = function(from, to, recursive=TRUE, overwrite=TRUE) {
  if (!dir.exists(paths=file.path(to))) { dir.create(path=file.path(to), recursive=recursive) }
  files = list.files(from, full.names=TRUE, recursive=FALSE)
  file.copy(from=files, to=file.path(to), recursive=recursive, overwrite=overwrite, copy.mode=TRUE)
  unlink(file.path(from), recursive=recursive)
  return(invisible(0))
}

#' Helper method for downloading files
#' @note Download observation data and add to list
#' @param url Character url to a file.
#' @param destination Character path to save downloaded file.
download_file = function(url, destination) {
  curl::curl_download(url=url, destfile=destination, quiet=FALSE, mode="w")
  message(paste("Download successful:", destination))
  return(invisible(0))
}

#' Helper method to interpolate points to raster
#' @param x Numeric x,y coordinates or n-dimensional vector of predictor variables.
#' @param y Numeric 1-d vector for target variable.
#' @param xi Numeric x,y coordinates or n-dimensional data.frame() for imputing values.
#' @param method Character method selection.
#' @param fun Function for caret trainControl or RandomFields model.
#' @param multicore Boolean flag for multicore computation using OpenMP.
#' @param spatial Boolean flag for whether xi contains x,y coordinates.
#' @param proj4string Character Proj4 projection string for raster image.
#' @param plot Boolean flag for plotting results.
#' @examples
#' Classification
#' data(iris)
#' intrain = caret::createDataPartition(iris$Species, p=0.7, list=FALSE)
#' train   = iris[ intrain,]
#' test    = iris[-intrain,]
#' x       = train[,1:4]
#' y       = train[,5]
#' xi      = test[,1:4]
#' yi      = test[,5]
#' p = interpolate_points(x=x, y=y, xi=xi, method="knn", multicore=FALSE, spatial=FALSE, plot=TRUE)
#' sum(p == yi) / length(p)
#' p = interpolate_points(x=x, y=y, xi=xi, method="xgbTree", multicore=TRUE, spatial=FALSE, plot=TRUE)
#' sum(p == yi) / length(p)
#' Spatial regression
#' x  = cbind(x=runif(100, 1, 100), y=runif(100, 1, 100))
#' y  = rowSums(xy)
#' xi = expand.grid(x=seq(1, 100, length.out=100), y=seq(1, 100, length.out=100))
#' p = interpolate_points(x=x, y=y, xi=xi, method="knn", fun=NA, multicore=FALSE, spatial=TRUE,
#'                        proj4string="+proj=longlat +ellps=WGS84 +datum=WGS84", plot=TRUE)
#' or
#' p = interpolate_points(x=x, y=y, xi=xi, method="rf", fun=NA, multicore=TRUE, spatial=TRUE,
#'                        proj4string="+proj=longlat +ellps=WGS84 +datum=WGS84", plot=TRUE)
#' Kriging
#' xi  = data.frame(x=seq(1, 100, length.out=100), y=seq(1, 100, length.out=100))
#' fun = RandomFields::RMexp() + RandomFields::RMtrend(mean=NA)
#' fun = ~1 + RandomFields::RMwhittle(scale=NA, var=NA, nu=NA) + RandomFields::RMnugget(var=NA)
#' p = interpolate_points(x=x, y=y, xi=xi, method="kriging", fun=fun, multicore=FALSE, spatial=TRUE,
#'                        proj4string="+proj=longlat +ellps=WGS84 +datum=WGS84", plot=TRUE)
interpolate_points = function(x, y, xi=NA, method="knn", fun=NA, multicore=FALSE, spatial=TRUE,
                              proj4string="+proj=longlat +ellps=WGS84 +datum=WGS84", plot=FALSE, ...) {
  if (method != "kriging") {
    if (multicore == TRUE) { doMC::registerDoMC(cores=parallel::detectCores()-1) }
    if (!is.na(fun)) {
      ctl = fun
    } else {
      ctl  = caret::trainControl(method="repeatedcv", number=3, repeats=3, search="grid")
    }
    mod  = caret::train(x, y, method=method, trControl=ctl, preProcess=c("center","scale"), tuneLength=10)
    yhat = predict(mod, newdata=xi)
    if (spatial == TRUE) {
      coords = colnames(x) %in% c("x","y","lat","lon","latitude","longitude")
      proj4crs = sp::CRS(proj4string)
      yhat = sp::SpatialPixelsDataFrame(points=sp::SpatialPoints(xi[,coords]), data=data.frame(y=yhat),
                                        proj4string=proj4crs)
    }
  } else if (method == "kriging" & spatial == TRUE) {
    if (length(unique(round(diff(xi[,1]),6))) > 2 | length(unique(round(diff(xi[,2]),6))) > 2) {
      stop("Please enter an evenly spaced grid for xi parameter")
    }
    if (any(duplicated(xi[,1]) | duplicated(xi[,2]))) {
      warning("Parameter xi appears to be a grid, creating list of unique coordinates", immediate.=TRUE)
      xi = data.frame(x=min(xi[,1]):max(xi[,1]), y=min(xi[,2]):max(xi[,2]))
    }
    proj4crs = sp::CRS(proj4string)
    coords = colnames(x) %in% c("x","y","lat","lon","latitude","longitude")
    spdf = sp::SpatialPointsDataFrame(coords=x[,coords], data=data.frame(y=y), proj4string=proj4crs)
    fit  = RandomFields::RFfit(fun, data=spdf)
    yhat = RandomFields::RFinterpolate(fit, x=xi[,coords], grid=TRUE, data=spdf)
  } else {
    stop("Parameters incorrectly specified")
  }
  if (spatial == TRUE) {
    yhat = raster::raster(yhat)
    if (plot == TRUE) { raster::plot(yhat, col=viridis::magma(256)) }
  }
  return(yhat)
}

#' Helper method to generate raster image
#' @note Used to generate random maps, point value maps, or matrix maps.
#' @note A RandomFields model can also be passed to generator.
#' @param values Numeric matrix, single value, "random", or "scrf" for spatially correlated random fields.
#' @param path Character file path for saving the result.
#' @param format Character raster file format.
#' @param params Optional numeric vector of parameters to pass to distribution function.
#' @param generator Optional statistical generator for drawing random values. Default is runif().
#' @param size Optional numeric vector of length 2 for the n*n size of the map.
#' @param extent Optional numeric vector of length 4 for the image bounds.
#' @param crs Optional character vector proj4-compatible coordinate reference system.
#' @examples
#' f  = "/Users/julia/rasters"
#' g = RandomFields::RMexp(var=3, scale=5) + RandomFields::RMnugget(var=3)
#' r = generate_raster("randomfield", size=c(100,100), generator=g, filename=f)
generate_raster = function(values="random", size=c(10,10), bounds=c(-1,1,-1,1), generator=runif,
                           params=c(prod(size)), filename=NA, format="GTiff",
                           proj4string="+proj=longlat +ellps=WGS84 +datum=WGS84", plot=TRUE, ...) {
  if (values != "random" & class(values) == "numeric" & class(values) != "matrix") {
    if (length(values) != prod(size)) { stop("Number of values must equal the product of size") }
    values = matrix(values, nrow=size[1], ncol=size[2])
  } else if (values == "random" & class(generator) == "function") {
    values = do.call(generator, as.list(params))
    values = matrix(values, nrow=size[1], ncol=size[2])
  } else if (values == "randomfields" & class(generator)[1] == "RMmodel") {
    x = 1:size[1]
    y = 1:size[2]
    values = RandomFields::RFsimulate(generator, x, y, grid=TRUE)
  } else {
    stop("Parameters incorrectly specified")
  }
  ri = raster::raster(values)
  raster::extent(ri) = bounds
  raster::projection(ri) = raster::crs(proj4string)
  if (!is.na(filename) == TRUE) { raster::writeRaster(ri, filename=filename, format=format) }
  if (plot == TRUE) { raster::plot(ri, col=viridis::magma(256)) }
  return(ri)
}
