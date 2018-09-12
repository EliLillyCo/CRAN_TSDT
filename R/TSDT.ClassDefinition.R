#!usr/bin/dev/R
#################################################################################
# FILENAME   : TSDT.ClassDefinition.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 02/04/2013
# DESCRIPTION: Implementation of TSDT class.
#################################################################################

# Declare class TSDT_CutpointDistribution
setClass( "TSDT_CutpointDistribution", package = "TSDT" )

#' @title TSDT
#' @description TSDT is a container class for TSDT samples and metadata.
#' @seealso \link{TSDT}, \linkS4class{TSDT_Sample},
#' \linkS4class{TSDT_CutpointDistribution}
#' @slot parameters List of parameters used in construction of TSDT samples.
#' @slot samples Vector of \linkS4class{TSDT_Sample} objects.
#' @slot superior_subgroups data.frame containing summary statistics for superior
#' subgroups
#' @slot cutpoints An object of class  \linkS4class{TSDT_CutpointDistribution}.
#' @slot distributions A list of distributions of TSDT statistics.
#' @return Object of class TSDT
setClass( "TSDT",
        
         representation = representation(
                           parameters = "list",
                           samples = "list",
                           superior_subgroups = "data.frame",
                           cutpoints = "TSDT_CutpointDistribution",
                           distributions = "list" ),

            validity = function( object ){

              validity_string <- NULL 
         
              for( i in 1:length(object@samples) ){

                # Test depth of tree does not exceed maxdepth. ctree has a
                # different node-numbering convention so the logiv for
                # determining tree depth is different.
                if( object@parameters$tree_builder %in% c("rpart") ){
                  
                  if( floor( log2( max( object@samples[[i]]@subgroups$NodeID ) ) ) > object@parameters$maxdepth )
                      validity_string <- paste0( validity_string, "\nERROR: depth of tree in sample ", i, " exceeds maxdepth" )
                }else if( object@parameters$tree_builder %in% c("ctree") ){
                  
                  if( max( object@samples[[i]]@subgroups$Depth ) > object@parameters$maxdepth )
                      validity_string <- paste0( validity_string, "\nERROR: depth of tree in sample ", i, " exceeds maxdepth" )
                }
                
                if( min( object@samples[[i]]@subgroups$NodeSize ) < object@parameters$min_subgroup_n_trt )
                    validity_string <- paste0( validity_string, "\nERROR: a subgroup in sample ", i, " has fewer than min_subgroup_n_trt treatment arm subjects" )
                
                if( object@parameters$min_subgroup_n_control > 0 && NROW( object@samples[[i]]@subgroups ) > 1 ){
                  
                  for( j in 2:NROW( object@samples[[i]]@subgroups ) ){
                    
                    if( object@samples[[i]]@subgroups$Superior_Inbag_Subgroup[[j]] ){
                      
                      subgroup__ <- subgroup( splits = object@samples[[i]]@subgroups,
                                             node = object@samples[[i]]@subgroups$NodeID[[j]],
                                             xdata = object@samples[[i]]@inbag )
                      
                      control_n <- length( which( subgroup__$trt == object@parameters$trt_control ) )
                      
                      if( control_n < object@parameters$min_subgroup_n_control )
                          validity_string <- paste0( validity_string, "\nERROR: a subgroup in sample ", i, " has fewer than min_subgroup_n_control control arm subjects" )
                    }
                  }
                }
                
                if( length( which( object@samples[[i]]@subgroups$NodeID < 0 ) ) > 2 * object@parameters$rootcompete )
                    validity_string <- paste0( validity_string, "\nERROR: the number of competitor splits in sample ", i, " exceeds rootcompete" )
              }
              
              if( !is.null(validity_string) ){
                return( validity_string )
              }
              
              # else
              return( TRUE )
              
         },
         package = "TSDT" )

setMethod( f = "initialize", signature = "TSDT",
          
           definition = function( .Object, parameters, samples,
                                  superior_subgroups, cutpoints, distributions ){

            .Object@parameters = parameters
            .Object@samples = samples
            .Object@superior_subgroups = superior_subgroups
            .Object@cutpoints = cutpoints
            .Object@distributions = distributions
            
            validity_check <- validObject( .Object )
            
            if( validity_check != TRUE )
                stop( validity_check )
            
            return( .Object )
         })

setMethod( f = "show", signature = "TSDT",

           definition = function( object ){

            N_rows <- NROW( object@superior_subgroups )
            ROWS_TO_DISPLAY <- 10 # How many rows to show. Others are omitted from display.
             
            cat( "Summary Statistics for Superior Subgroups\n" )
 
            if( N_rows <= ROWS_TO_DISPLAY )
                base::print( object@superior_subgroups )
            else{
              base::print( object@superior_subgroups[1:ROWS_TO_DISPLAY,] )
              cat( paste0( "...\n...\n...\n\nNOTE: ", N_rows - ROWS_TO_DISPLAY, " rows omitted...\n" ) )
              cat( "      Use summary() method to see all subgroups.\n" )
            }

         })

#################################################################################
# Accessor methods for class data members                                       #
#################################################################################

#' @title Summary function for class TSDT.
#' @description Summary function for class TSDT.
#' @seealso \link{TSDT}
#' @param object An object of class \linkS4class{TSDT}.
#' @return A data.frame containing the superior subgroups identified by TSDT.
#' @export
#' @docType methods
#' @rdname summary-methods
setMethod( f = "summary", signature = "TSDT",
         
          definition = function( object ){
            
            return( object@superior_subgroups )

         })

parameters <- function( object ){return( object@parameters )}
samples <- function( object ){return( object@samples )}

#' @title Get distribution of cutpoints for subgroups.
#' @description Get distribution of cutpoints for subgroups.
#' @param object An object of class TSDT
#' @param subgroup A string decscription of a subgroup (optional)
#' @param subsub A string description of a sub-subgroup (optional)
#' @seealso \link{TSDT}
#' @return A vector containing the subgroup cutpoints.
#' @export
cutpoints <-function( object, subgroup = NULL, subsub = NULL ){

  if( is.null( subgroup ) )
      return( object@cutpoints )
  
  # else
  subgroup <- as.character( subgroup )
  
  if( !is.null( subsub ) )
      subsub <- as.character( subsub )
  
  return( get_cutpoints( object@cutpoints, subgroup, subsub ) )
}

## #' @title Get distribution of a specified summary statistic.
## #' @description Get distribution of a specified summary statistic.
## #' @param object An object of class TSDT
## #' @param statistic The desired summary statistic
## #' @param subgroup A string decscription of a subgroup (optional)
## #' @seealso \link{TSDT}
## #' @return A vector containing the distribution of values for the specified statistic.
## #' @export
## #' @docType methods
## #' @rdname distribution-methods
## distribution <- function( object, statistic, subgroup = NULL ){

##   if( is.null( subgroup ) )
##       return( object@distributions[[ statistic ]] )

##   #else
##   return( object@distributions[[ statistic ]][[ subgroup ]] )
## }



## END OF FILE
