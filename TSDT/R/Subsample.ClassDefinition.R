#!usr/bin/dev/R
#################################################################################
# FILENAME   : Subsample.ClassDefinition.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 02/27/2013
# DESCRIPTION: Implementation of Subsample container class
#################################################################################

#' @title Subsample
#' @description Subsmaple is a container class for subsamples.
#' @seealso \link{subsample}
#' @slot training Training data.
#' @slot validation Validation data.
#' @slot test Test data.
#' @return Object of class
setClass( "Subsample",
         
         representation( training = "data.frame",
                         validation = "data.frame",
                         test = "data.frame" ),
         package = "TSDT" )

## END OF FILE
