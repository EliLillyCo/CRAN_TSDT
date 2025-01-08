#!usr/bin/dev/R
#################################################################################
# FILENAME   : TSDT_Sample.ClassDefinition.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 02/04/2013
# DESCRIPTION: Implementation of TSDT_Sample class.
#################################################################################

#' @title TSDT_Sample
#' @description TSDT_Sample is a container class containing the in-bag and
#' out-of-bag data from a subsampled or bootstrapped dataset. This container
#' class also contains a data.frame containing the parsed tree that is fit on the
#' in-bag data.
#' @seealso \link{TSDT}
#' @slot inbag A data.frame containing in-bag data
#' @slot oob A data.frame containing out-of-bag data
#' @slot subgroups A data.frame containing a parsed tree
#' @return Object of class TSDT_Sample
setClass( "TSDT_Sample",
          representation = representation(
                           inbag = "data.frame",
                           oob = "data.frame",
                           subgroups = "data.frame" ),

            validity = function( object ){

              validity_string <- NULL

              # Validate values in object@subgroups
              if( NROW( object@subgroups ) > 0 ){

                if( any( object@subgroups$Superior_Inbag_Subgroup %nin% c( TRUE, FALSE ) ) ){
                  
                  validity_string <- paste0( validity_string,
                                            "\nERROR: Superior_Inbag_Subgroup must be in {FALSE,TRUE}" )
                }
              
              
                if( any( object@subgroups$Superior_OOB_Subgroup %nin% c( TRUE, FALSE ) ) ){
                  
                  validity_string <- paste0( validity_string,
                                            "\nERROR: Superior_OOB_Subgroup must be in {FALSE,TRUE}" )
                }
                
                
                if( any( object@subgroups$Superior_Inbag_Subgroup == FALSE &
                        !is.na( object@subgroups$External_Consistency ) ) )
                  validity_string <- paste0( validity_string,
                                            "\nERROR: external consistency must be NA ",
                                            "if the in-bag subgroup is not a superior ",
                                            "subgroup" )
                
                
                test_superior_oob_subgroup <- subset( object@subgroups$Superior_OOB_Subgroup, object@subgroups$Superior_Inbag_Subgroup == TRUE )
                test_external_consistency  <- subset( object@subgroups$External_Consistency, object@subgroups$Superior_Inbag_Subgroup == TRUE )
                
                if( any( test_external_consistency %nin% c( TRUE, FALSE ) ) ){
                  
                  validity_string <- paste0( validity_string,
                                            "\nERROR: external consistency must not be in {FALSE,TRUE} ",
                                            "if Superior_Inbag_Subgroup is TRUE" )
                  
                  
                }
                
                if( any( test_superior_oob_subgroup != test_external_consistency ) )
                  validity_string <- paste0( validity_string,
                                            "\nERROR: if in-bag subgroup is a superior ",
                                            "subgroup, then external consistency is ",
                                            "true if OOB subgroup is also superior, ",
                                            "false otherwise" )
              }
              
              if( !is.null(validity_string) )
                return( validity_string )
              
              # else
              return( TRUE )
              
            },
         package = "TSDT" )

setMethod( f = "initialize", signature = "TSDT_Sample",

          definition = function( .Object,
            inbag = data.frame( NULL ), 
            oob = data.frame( NULL ),
            subgroups = data.frame( NULL ) ){
           
           .Object@inbag = inbag
           .Object@oob = oob
           .Object@subgroups = subgroups

           validity_check <- validObject( .Object )

           if( validity_check != TRUE )
             stop( validity_check )
           
           return( .Object )
         })

#################################################################################
# Accessor methods for class data members                                       #
#################################################################################

# Accessor method for in-bag data
setGeneric( name = "inbag",
            def = function(object){standardGeneric("inbag")} )

setMethod( f = "inbag", signature = "TSDT_Sample",
          
          definition = function( object ){
            
            return( object@inbag )
            
          })

# Accessor method for out-of-bag data
setGeneric( name = "oob",
            def = function(object){standardGeneric("oob")} )

setMethod( f = "oob", signature = "TSDT_Sample",

           definition = function( object ){

             return( object@oob )

         })

# Accessor method for the TSDT tree
setGeneric( name = "subgroups",
            def = function(object){standardGeneric("subgroups")} )

setMethod( f = "subgroups", signature = "TSDT_Sample",

           definition = function( object ){

             return( object@subgroups )

           })

## END OF FILE
