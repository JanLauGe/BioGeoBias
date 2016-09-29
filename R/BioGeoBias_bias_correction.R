
#' @title Choose bias correction approach and generate bias correction data
#' @description This function takes a species name and a rank, and returns a
#' dataset suitable for bias correction by invoking
#' \code{create_target_group_background_data} or
#' \code{create_bias_grid}, depending on the number of species occurrence
#' records available for the species and the spatial extent of the selection.
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
#' @keywords GBIF, sampling bias, SDM
#' @export
bias_correction <- function(
  species_name,
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

  # check valid variable values for the function
  if(length(target_rank) > 1){
    # If no target rank is given, use class as default
    message('\nNo valid target rank supplied, using "class"...\n')
    target_rank <- 'class'
  }
  if(class(target_rank) != 'character'){
    # If invalid target rank is given, stop with error
    stop('value supplied to argument "target_rank" is not a character. Try rgbif::taxrank() to get a summary of valid inputs for this argument')
  }
  target_rank = match.arg(target_rank)
  if(class(species_name) != 'character'){
    stop('value supplied to argument "species_name" is not a character. This should be the latin binomial of the species you want to model')
  }
  if(!is.null(kingdom)){
    kingdom = match.arg(kingdom)
  }
  for(v in c(lonMin, lonMax, latMin, latMax)){
    if(!is.numeric(v) & !is.integer(v) & !is.null(v)){
      stop('values supplied to arguments lonMin, lonMax, latMin, latMax are not numeric. These should be numerical values giving the minimum and maximum longitute and latitude for the extent of the occurrence search.')
    }
  }
  # If min & max coordinates are supplied, use spatial filter
  latFilter <- ifelse(!is.null(latMin) & !is.null(latMax), TRUE, FALSE)
  lonFilter <- ifelse(!is.null(lonMin) & !is.null(lonMax), TRUE, FALSE)

  # GBIF name query --------------------------------------------------------------

  # Get name key of the target group from gbif backbone taxonomy.
  # kingdom is used if supplied, but ignored otherwise
  message('\nChecking GBFI for target group taxon key...\n')
  NameData <- name_backbone(
    name = species_name,
    kingdom = kingdom)

  # Check if name_backbone returned a valid result
  if(NameData$matchType == 'NONE') stop('No valid taxon key for the target species. Are you sure you provied a valid latin binomial name?')
  # Check if name_backbone name matching confidence is 95%
  else if(NameData$confidence < 95){
    warning('Name matching confidence was less than 95%. Please make sure that the matched taxon is the one you want.')
    message(paste0(
      'Name matched with ', NameData$confidence, ' confidence \n'
    ))}
  # Check if name_backbone returned a valid taxon key for target group
  if(is.null(NameData[paste0(target_rank ,'Key')][[1]])){
    stop('No valid taxon key for the target group. Perhaps you should try a different rank?')
  }

  # Get the taxon key of the target group
  Key <- NameData[paste0(target_rank ,'Key')]
  # Print information about the matched species
  message(paste0(
    'Name matched to:\n  ',
    NameData$scientificName,
    ', Taxon key: ',
    NameData$usageKey))
  # Print information about the matched target group
  message(paste0(
    'Target group selected:\n  ',
    target_rank, ' ', NameData[target_rank],
    ', Taxon key: ', Key), '\n')

  # Get the number of occurrences in the target group
  Count <- occ_count(
    taxonKey = Key,
    georeferenced = TRUE)


  ifelse(Count > 200000, useMap <- TRUE, useMap <- FALSE)


  if(useMap == FALSE){
    # Use target group background ==============================================
    this <- create_target_group_background_data(
      species_name = species_name,
      target_rank = target_rank,
      lonMin = NULL,
      lonMax = NULL,
      latMin = NULL,
      latMax = NULL,
      limit = 200000)
  }else if(useMap == TRUE){
    # Use map api ==============================================================
    this <- create_bias_grid(
      taxonkey = as.numeric(Key),
      lonMin = NULL,
      lonMax = NULL,
      latMin = NULL,
      latMax = NULL)
  }
#
#   # Use accumulation curve =====================================================
#   this <- species_accumulation_curve_estimation(
#     species_name = species_name,
#     target_rank = target_rank,
#     lonMin = NULL,
#     lonMax = NULL,
#     latMin = NULL,
#     latMax = NULL)
# }
}
