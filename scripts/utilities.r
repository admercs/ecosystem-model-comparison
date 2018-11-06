#=================================================================================================#
# utilities.r
# Generic utility methods
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

#' Coefficient of determination (R2) helper method.
#' @param x Numeric observation values.
#' @param y Numeric prediction values.
#' @return Numeric value.
r_squared = function(x, y) {
  if (length(x) != length(y)) { stop("Lengths of y and y_hat differ") }
  return(cor(x, y, method="pearson")^2)
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
