#!usr/bin/dev/R
#################################################################################
# FILENAME   : Bootstrap.ClassDefinition.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 09/09/2012
# DESCRIPTION: Implementation of Bootstrap and BootstrapStatistic container
#              classes.
#################################################################################

#' @title Bootstrap
#' @description Bootstrap is a container class for bootstrap samples.
#' @seealso \linkS4class{BootstrapStatistic}, \link{bootstrap}
#' @slot inbag In-bag bootstrap sample.
#' @slot oob Out-of-bag bootstrap sample.
#' @return Object of class Bootstrap
setClass( "Bootstrap",
         representation( inbag = "data.frame", 
                         oob = "data.frame" ),
         package = "TSDT" )


setMethod( f = "initialize", signature = "Bootstrap",

          definition = function( .Object,
            inbag = data.frame( NULL ), 
            oob = data.frame( NULL ) ){

           .Object@inbag = inbag
           .Object@oob = oob
           
           return( .Object )
         })



#' @title BootstrapStatistic
#' @description BootstrapStatistic is a container class for bootstrap samples augmented with a computed statistic.
#' @seealso \linkS4class{Bootstrap}, \link{bootstrap}
#' @slot statname The name of a (possibly user-defined) statistic to compute on the bootstrap sample.
#' @slot arglist A list of arguments passed to the function referenced by statname.
#' @slot variable The name of the variable on which to compute statname.
#' @slot inbag_stat The value of statname for the in-bag bootstrapped sample. 
#' @slot oob_stat The value of statname for the out-of-bag bootstrapped sample.
#' @return Object of class BootstrapStatistic
setClass( Class = "BootstrapStatistic", contains = "Bootstrap",
         representation( statname = "character",
                         arglist = "list",
                         variable = "character",
                         inbag_stat = "numeric",
                         oob_stat = "numeric" ) )


setMethod( f = "initialize", signature = "BootstrapStatistic",

          definition = function( .Object, inbag, oob, statname, arglist,
            variable, inbag_stat, oob_stat ){

            .Object@inbag = inbag
            .Object@oob = oob
            .Object@statname = statname
            .Object@arglist = arglist
            .Object@variable = variable
            .Object@inbag_stat = inbag_stat
            .Object@oob_stat = oob_stat
            
            return( .Object )
          })



## END OF FILE
