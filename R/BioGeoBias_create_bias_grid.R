
#' @title Create a Bias Grid from GBIF web map API
#' @description This function takes a species name and a rank, and returns a
#' raster of number of species occurrence records per cell which can be used as
#' a sampling grid for bias correction in SDM applications.
#' @param species_name (character) a valid taxon name, presumably a species name.
#' @param target_rank (character) a taxon rank for the target group (optional).
#' @param kingdom (character) the kingdom of the species. Can be supplied to
#' avoid possible confusion when matching names. Should be one of c('animalia',
#' 'plantae','archaea','bacteria','fungi','protozoa','viruses') (optional).
#' @param lonMin minimum longitude of a rectancular bounding box restricting the
#' search for species occurrences of the target group (optional).
#' @param lonMax maximum longitude of a rectancular bounding box restricting the
#' search for species occurrences of the target group (optional).
#' @param latMin minimum latitude of a rectancular bounding box restricting the
#' search for species occurrences of the target group (optional).
#' @param latMax maximum latitude of a rectancular bounding box restricting the
#' search for species occurrences of the target group (optional).
#' @author Jan Laurens Geffert, \email{laurensgeffert@@gmail.com}
#' @details hello
#' @references \url{http://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors}
#' @seealso \url{http://www.ecography.org/accepted-article/performance-tradeoffs-target-group-bias-correction-species-distribution-models}
#' @keywords GBIF, sampling bias, web map tile, raster, grid
#' @export
create_bias_grid <- function(
  taxonkey = 1,
  target_rank = rgbif::taxrank(),
  kingdom = c(NULL,
              'animalia',
              'plantae',
              'archaea',
              'bacteria',
              'fungi',
              'protozoa',
              'viruses'),
  lonMin = NULL,
  lonMax = NULL,
  latMin = NULL,
  latMax = NULL){

  call_map_api(
    taxonkey = taxonkey,
    nbreaks = 50)
  }


