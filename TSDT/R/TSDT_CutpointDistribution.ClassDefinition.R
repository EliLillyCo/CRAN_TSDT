#!usr/bin/dev/R
#################################################################################
# FILENAME   : TSDT_CutpointDistribution.ClassDefinition.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 02/04/2013
# DESCRIPTION: Implementation of TSDT_CutpointDistribution class. This class
#              stores the specific numeric cutpoints for each generic subgroup
#              defined on a continuous split variable. If the subgroup contains
#              more than one split variable a distribution of numeric cutpoints
#              is collected for each continuous split variable in the subgroup
#              definition.
#################################################################################

requireNamespace( "hash", quietly = TRUE )

#' @title TSDT_CutpointDistribution
#' @description Implementation of TSDT_CutpointDistribution class. This class
#  stores the specific numeric cutpoints for each generic subgroup defined on a
#' continuous split variable. If the subgroup contains more than one split
#' variable a distribution of numeric cutpoints is collected for each continuous
#' split variable in the subgroup definition.
#' @seealso \link{TSDT}, \link[hash]{hash}
#' @slot Cutpoints An object of class \link[hash]{hash-class}
#' @return Object of class TSDT_CutpointDistribution
setClass( "TSDT_CutpointDistribution",
          representation = representation( Cutpoints = "hash" ),
          package = "TSDT" )

setMethod( f = "initialize", signature =  "TSDT_CutpointDistribution",
          definition = function(.Object, Cutpoints = NULL ){
            requireNamespace( "hash", quietly = TRUE )
            if( is.null( Cutpoints ) || !is(Cutpoints, 'hash') ){
              .Object@Cutpoints = hash()
            }else{
              .Object@Cutpoints = Cutpoints
            }
            return( .Object )
          })


setGeneric( name = "set_cutpoints",
            def = function( .Object, subgroup ){standardGeneric("set_cutpoints")} )

setMethod( f = "set_cutpoints", signature = "TSDT_CutpointDistribution",

          definition = function( .Object,
            subgroup ){
            
            generic <- get_generic_subgroup( subgroup )
            subgroups <- strsplit( subgroup, " & " )[[1]]

            generic_subgroups <- subgroups
            
            for( k in 1:length(generic_subgroups) )
              generic_subgroups[[k]] <- get_generic_subgroup( generic_subgroups[[k]] )
            
            split_type <- ifelse( grepl( subgroups, pattern = "%in%|%nin%" ), "nominal", "continuous" )

            # Extract numeric cutpoints
            cutpoint <- rep( NA, length(subgroups) )

            for( j in 1:length(subgroups) ){
                
                if( split_type[j] == "continuous" ){
                    
                    if( grepl( subgroups[j], pattern = ">=" ) ){
                        split_value_start <- regexpr( pattern = ">=", subgroups[j], useBytes = TRUE ) + 2
                    }else if( grepl( subgroups[j], pattern = "<=" ) ){
                        split_value_start <- regexpr( pattern = "<=", subgroups[j], useBytes = TRUE ) + 2
                    }else if( grepl( subgroups[j], pattern = ">" ) ){
                        split_value_start <- regexpr( pattern = ">", subgroups[j], useBytes = TRUE ) + 1
                    }else if( grepl( subgroups[j], pattern = "<" ) ){
                        split_value_start <- regexpr( pattern = "<", subgroups[j], useBytes = TRUE ) + 1
                    }
                    
                    cutpoint[j] <- substr( subgroups[j], start = split_value_start, stop = nchar(subgroups[j] ) )

                    rm( split_value_start )
                }
            }
            
            # key already exists in hash
            if( has.key( generic, .Object@Cutpoints )[[1]] ){
                
                for( i in 1:length( subgroups ) )
                    if( split_type[i] == "continuous" )
                        eval( parse( text = paste0( ' .Object@Cutpoints[[ generic ]]$`',
                                         generic_subgroups[i], '`<- c( .Object@Cutpoints[[ generic ]]$`',
                                         generic_subgroups[i], '`,as.numeric(', cutpoint[i], '))' ) ) )
                
            }

            # key does not already exist in hash
            else{
                for( i in 1:length( subgroups ) )
                    if( split_type[i] == "continuous" )
                        eval( parse( text = paste0( '.Object@Cutpoints[[ generic ]]$`',
                                        generic_subgroups[i], '`<- as.numeric( cutpoint[i] )' ) ) )
            }
            
            invisible()
        })

#' @title get_cutpoints
#' @description Accessor method for cutpoints slot in TSDT objects.
#' @details The summary results from TSDT provide a set of 'anonymized' subgroups
#' in a form similar to 'X1<xxxxx'. The variable X1 may have been selected as a
#' splitting variable in several bootstrapped samples. The exact numerical
#' cutpoint for X1 could vary from one sample to the next. The get_cutpoints
#' method returns all the numerical cutpoints associated with this subgroup. If
#' the subgroup is a compound subgroup defined on more than one spliting variable
#' the user can specify the 'subsub' parameter to get the cutpoints associated
#' with a particular component of the subgroup.
#' @param .Object A TSDT object.
#' @param subgroup The anonymized subgroup.
#' @param subsub A particular component of the subgroup to retrieve.
#' @export
#' @docType methods
#' @rdname get_cutpoints-methods
#' @examples
#' \dontrun{
#' example( TSDT )
#' ## You can access the cutpoints slot of a TSDT object directly
#' ex2@cutpoints
#'
#' ## You can also use the accessor method
#' get_cutpoints( ex2@cutpoints, subgroup = 'X1<xxxxx' )
#' 
#' ## Retrieving a compound subgroup defined on multiple splits
#' get_cutpoints( ex2, subgroup = 'X1<xxxxx & X1>=xxxxx' )
#'
#' ## Retrieving a single component from the compound subgroup
#' get_cutpoints( ex2, subgroup = 'X1<xxxxx & X1>=xxxxx', subsub = 'X1>=xxxxx' )
#' }
setGeneric( name = "get_cutpoints",
            def = function( .Object, subgroup, subsub = NULL ){standardGeneric("get_cutpoints")} )


#' @rdname get_cutpoints-methods
#' @aliases get_cutpoints,ANY-method
setMethod( f = "get_cutpoints", signature = "TSDT_CutpointDistribution",

          definition = function( .Object,
            subgroup = character,
            subsub = NULL ){

            if( is.null( subsub ) )
              return( .Object@Cutpoints[[ subgroup ]] )

            else{
              sg <- .Object@Cutpoints[[ subgroup ]]
              eval( parse( text = paste0( "return( sg$`", subsub, "`)" ) ) )
            }
            
         })


## Define a get_cutpoints accesor methed when .Object is of class TSDT
#' @rdname get_cutpoints-methods
#' @aliases get_cutpoints,ANY-method
setMethod( f = "get_cutpoints", signature = "TSDT",

          definition = function( .Object,
            subgroup = character,
            subsub = NULL ){

            ## Call TSDT::get_cutpoints method that takes an .Object of class
            ## TSDT_CutpointDistribution.
            return( get_cutpoints( .Object@cutpoints, subgroup = subgroup, subsub = subsub ) )
            
          })


setGeneric( name = "delete_cutpoints",
            def = function( .Object ){standardGeneric("delete_cutpoints")} )

setMethod( f = "delete_cutpoints", signature = "TSDT_CutpointDistribution",

          definition = function( .Object ){

            clear( .Object@Cutpoints )
            obj <- deparse( substitute( .Object ) )
            rm( list = obj, envir = parent.frame() )
         })

setGeneric( name = "append_cutpoints",
            def = function( .Object, new__ ){standardGeneric("append_cutpoints")} )

setMethod( f = "append_cutpoints", signature = "TSDT_CutpointDistribution",

          definition = function( .Object, new__ ){

            for( key in keys( new__@Cutpoints ) ){

              values <- new__@Cutpoints[[key]]

              # A key in new__ already exists in .Object
              if( has.key( key, .Object@Cutpoints ) ){

                for( v in 1:length(values) ){

                  # Append each individual split value to the appropriate subsub
                  eval( parse( text = paste0( ".Object@Cutpoints[[key]]$`", names(values)[[v]],
                                 "`<-c(.Object@Cutpoints[[key]]$`", names(values)[[v]],
                                "`,", values[[v]], ")" ) ) )
                }
              }
              # A key in new__ does not exist in .Object
              else{
                .Object@Cutpoints[[ key ]] <- values
              }
            }
          })


## END OF FILE
