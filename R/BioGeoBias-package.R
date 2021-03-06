
#' @name BioGeoBias-package
#' @aliases BioGeoBias
#' @docType package
#' @title Package overview
#' @description This package provides various functions to assess and account
#' for sampling bias in species occurrence data from GBIF. Functions are
#' hierarchically structured to allow users with varying levels of expertise to
#' either use an out-of-the box solution or to apply the specific method they
#' deem most suitable to their needs. Users new to R should start off using the
#' \code{bias_correction} function.
#'
#' @section Beginners:
#'
#' \describe{
#'  \item{\strong{bias_correction}}{
#'  If you are not sure which bias correction approach is the right one
#'  for your particular case, try using the \code{bias_correction} function.
#'  It will try to automatically choose the approporate method for the species
#'  in question based on the number of records available for the target group.
#'  Pay close attention to warning messages and erros, as they will give some
#'  indication whether the approach was successful or not, and which parameters
#'  might be causing problems.}}
#'
#' @section Advanced Users:
#'
#' \describe{
#'  \item{\strong{create_target_group_background_data}}{
#'  This function uses the \code{name_backbone} and \code{occ_search}
#'  functions of the \pkg{rgbif} package to generate a set of occurrence records
#'  from the target group. This species occurrence dataset can subsequently be
#'  used as background (pseudo-absence) dataset in any presence-only species
#'  distribution modelling approach}}
#'
#' \describe{
#'  \item{\strong{create_bias_grid}}{
#'  This function uses the web tile service of the GBIF map API to
#'  generate a map of species occurrence record density. The output is a raster
#'  file that can be used as a bias grid, for example in the Maxent GUI
#'  application. This method is particularly useful when the target group is
#'  extremely data or species rich and widespread, because it circumvents the
#'  download of large quantities of data by using the aggregated map format
#'  instead.}}
#'
#' -----------------------------------------------------------------------------
#'
#' @section Issues to resolve:
#' \describe{
#'  \item{issue title}{
#'   \itemize{
#'    \item what to do
#'   }
#'  }
#'  \item{rename functions}{
#'   \itemize{
#'    \item Need some examples
#'   }
#'  }
#'  \item{add shapefile and gridcell support}{
#'   \itemize{
#'    \item Need some examples
#'   }
#'  }
#' }
#'
#' @section Plan for Future Functionality Additions:
#'  \describe{
#'   \item{Species accumulation curves}{
#'   \itemize{
#'    \item Add functionality to compute sample completeness using species
#'    accumulation curves.
#'    }
#'   }
#'   \item{Map API in rgbif}{
#'   \itemize{
#'    \item Create a pull request to include the map API call in the rgbif package
#'    }
#'   }
#'  }
#'
#' @importFrom dplyr select
#' @importFrom dplyr filter
#' @importFrom stats complete.cases
#' @importFrom magrittr "%>%"
#' @importFrom magrittr "%<>%"
#' @importFrom magrittr "%$%"
#' @importFrom rgbif taxrank
#' @importFrom rgbif name_backbone
#' @importFrom rgbif occ_count
#' @importFrom rgbif occ_search
#' @importFrom utils read.csv
#' @importFrom utils read.table
#' @importFrom utils write.table
#' @importFrom raster stack
#' @importFrom raster reclassify
#' @importFrom raster plot
#' @importClassesFrom raster RasterLayer


NULL
