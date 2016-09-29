
#' @title Get background data using the Target Group Background Approach (TGBA)
#' @description This function takes a species name and a rank, and returns a set
#' of background data for species distribution models of the species using the
#' Target Group Background Approach (TGBA).
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
#' @details With the TGBA, occurrence records of other species in the given
#' taxonomic group are considered background points. Assuming that sampling
#' intensity is spatially evenly distributed across all species in the group
#' these occurrences become a proxy of the spatial distribution of sampling
#' effort. Using them as background data results in models that explicitly
#' account for sampling bias and can create improved predictions of species
#' distributions.
#' The function connects to the GBIF API using functions from the R
#' package rgbif. First, the name of the species is matched against the GBIF
#' backbone taxonomy. Next, the name key of the taxon corresponding to the
#' parent taxon of the given rank is extracted. An occurrence search for
#' georeferenced records of that taxon returns the dataset we want.
#' @references \url{http://en.wikipedia.org/wiki/List_of_Crayola_crayon_colors}
#' @seealso \url{http://www.ecography.org/accepted-article/performance-tradeoffs-target-group-bias-correction-species-distribution-models}
#' @keywords sampling bias, TGBA, SDM
#' @export
create_target_group_background_data <- function(
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
  latMax = NULL,
  yearMin = NULL,
  yearMax = NULL,
  limit = 200000){

  #issues:
  # - include spatial filter
  # - include date filter
  # - include functionality to use custom species list instead of taxon


  # Validate variables -----------------------------------------------------------
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

  # Error handling
  if(Count == 0) break('No records found for target group')
  if(Count < 1000){
    warning(paste0('Warning: Less than 1,000 observations.
            Consider using a higher rank than ', target_rank, ' as target group'))
  }
  message(paste0('Total number of records: ', Count), '\n')

  # GBIF occurrences query -------------------------------------------------------
  occs <- occ_search(
    taxonKey = Key,
    return = 'data',
    limit = limit,
    decimalLatitude = ifelse(latFilter,
      paste(latMin,latMax,sep=','),
      paste('-90','90',sep=',')),
    decimalLongitude = ifelse(lonFilter,
      paste(lonMin,lonMax,sep=','),
      paste('-180','180',sep=',')),
    yearMin = yearMin,
    yearMax = yearMax,
    hasCoordinate = TRUE,
    hasGeospatialIssue = FALSE)

  if(is.null(dim(occs))){
    stop('No records found for the parameters specified')
  }
  else if(nrow(occs) >= 199999){
    warning('It looks like the query has reached the GBIF API hard limit of 200000 records per search. Consider selecting a lower taxonomic rank for your target group using the "target_rank" argument or restrict the spatial extent of your query using "latMin", "latMax", "lonMin", and "lonMax".')
  }

  return(occs)
}



