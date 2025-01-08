#!usr/bin/dev/R
#################################################################################
# FILENAME   : subsample.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 09/09/2012
# DESCRIPTION: 
#################################################################################

get_subsamples <- function( x, trt, trt_control, training_fraction, validation_fraction, test_fraction ){
  
  if( !is.data.frame(x) )
      x <- data.frame( x )
  
  control      <- subset( x, trt == trt_control )
  experimental <- subset( x, trt != trt_control )

  n_control      <- NROW( control )
  n_experimental <- NROW( experimental )
  
  # Compute sample size for each dataset
  n_experimental_validation <- round(validation_fraction * n_experimental)
  n_experimental_test <- round(test_fraction * n_experimental)
  n_experimental_training <- n_experimental - n_experimental_test - n_experimental_validation
  
  experimental_training   <- NULL
  experimental_validation <- NULL
  experimental_test       <- NULL

  # Randomly allocate observation to test and validation datasets, and
  # allocate remaining observations to training.
  if( n_experimental > 0 ){
    
    experimental_rows <- 1:n_experimental

    experimental_training_rows <- NULL
    experimental_validation_rows <- NULL
    experimental_test_rows <- NULL
    
    if( n_experimental_validation > 0 ){
      experimental_validation_rows <- sample( experimental_rows, size = n_experimental_validation, replace = FALSE )
      experimental_rows <- experimental_rows[ -experimental_validation_rows ]
    }
    
    if( n_experimental_test > 0 ){
      experimental_test_rows <- sample( experimental_rows, size = n_experimental_test, replace = FALSE )
      experimental_rows <- experimental_rows[ -experimental_test_rows ]
    }
    
    if( n_experimental_training > 0 )
        experimental_training_rows <- experimental_rows
    

    experimental_training   <- experimental[ experimental_training_rows, ]
    experimental_validation <- experimental[ experimental_validation_rows, ]
    experimental_test       <- experimental[ experimental_test_rows, ]
    
  }

  # Repeat for control arm
  n_control_validation <- round(validation_fraction * n_control)
  n_control_test <- round(test_fraction * n_control)
  n_control_training <- n_control - n_control_test - n_control_validation

  control_rows <- 1:n_control

  control_training_rows <- NULL
  control_validation_rows <- NULL
  control_test_rows <- NULL

  if( n_control_validation > 0 ){
    control_validation_rows <- sample( control_rows, size = n_control_validation, replace = FALSE )
    control_rows <- control_rows[ -control_validation_rows ]
  }

  if( n_control_test > 0 ){
    control_test_rows <- sample( control_rows, size = n_control_test, replace = FALSE )
    control_rows <- control_rows[ -control_test_rows ]
  }

  if( n_control_training > 0 )
    control_training_rows <- control_rows

  
  control_training   <- control[ control_training_rows, ]
  control_validation <- control[ control_validation_rows, ]
  control_test       <- control[ control_test_rows, ]
  
  training   <- rbind( control_training, experimental_training )
  validation <- rbind( control_validation, experimental_validation )
  test       <- rbind( control_test, experimental_test )
    
  rownames( training ) <- NULL
  rownames( validation ) <- NULL
  rownames( test ) <- NULL
  
  return( new( "Subsample", training = training, validation = validation, test = test ) )
}

#' @title subsample
#' @description Generate a vector of subsamples.
#' @details Each subsample will contain training, validation, and test data in
#' proportions specified by the user. If a treatment variable is supplied the
#' ratio of treatments will be preserved as closely as possible.
#' @seealso \linkS4class{Subsample}
#' @param x <Source data to subsample.
#' @param trt Treatment variable. (optional)
#' @param trt_control Value for treatment control arm. Defaulte value is
#' 'Control'.
#' @param training_fraction Fraction of source data to include in training
#' subsample.
#' @param validation_fraction Fraction of source data to include in validation
#' subsample.
#' @param test_fraction Fraction of source data to include in test subsample.
#' @param n_samples Number of subsamples to generate.
#' @return Vector of objects of class Subsample.
#' @examples
#' ## Generate example data frame containing response and treatment
#' N <- 50
#' x <- data.frame( runif( N ) )
#' names( x ) <- "response"
#' x$treatment <- factor( sample( c("Control","Experimental"), size = N,
#'                        prob = c(0.8,0.2), replace = TRUE ) )
#' 
#' ## Generate two subsamples
#' ex1 <- subsample( x,
#'                   training_fraction = 0.9,
#'                   test_fraction = 0.1,
#'                   n_samples = 2 )
#' 
#' ## Generate two subsamples preserving treatment ratio
#' ex2 <- subsample( x,
#'                   trt = x$treatment,
#'                   trt_control = "Control",
#'                   training_fraction = 0.7,
#'                   validation_fraction = 0.2,
#'                   test_fraction = 0.1,
#'                   n_samples = 2 )
#' @export
subsample <- function( x,
                       trt = NULL,
                       trt_control = "Control",
                       training_fraction = NULL,
                       validation_fraction = NULL,
                       test_fraction = NULL,
                       n_samples = 1 ){
  
  if( !is.data.frame( x ) )
      x <- as.data.frame( x )
  
  if( is.null( trt ) ){
    trt <- rep( "SINGLE_ARM", NROW(x) )
    trt_control <- "SINGLE_ARM"
  }
  
  # Validate selections for sampling fractions
  if( !is.null(training_fraction) && !is.null(validation_fraction) && is.null(test_fraction) ){
    test_fraction <- round(1 - training_fraction - validation_fraction, digits = 6)
  }else if( !is.null(training_fraction) && is.null(validation_fraction) && !is.null(test_fraction) ){
    validation_fraction <-round( 1 - training_fraction - test_fraction, digits = 6)
  }else if( is.null(training_fraction) && !is.null(validation_fraction) && !is.null(test_fraction) ){
    training_fraction <- round(1 - test_fraction - validation_fraction, digits = 6)
  }else if( is.null(training_fraction) && is.null(validation_fraction) && is.null(test_fraction) ){
    stop( "ERROR: At least two of {training_fraction,validation_fraction,test_fraction} must be non-null" )
  }
  
  if( training_fraction < 0 || training_fraction > 1 ||
     validation_fraction < 0 || validation_fraction > 1 ||
     test_fraction < 0 || test_fraction > 1 ){
    stop( "ERROR: training_fraction, validation_fraction, and test_fraction must all be in [0,1]" )
  }
  
  if( training_fraction + validation_fraction + test_fraction - 1 > 1E-06 ){
    stop( "ERROR: training_fraction + validation_fraction + test_fraction must equal 1" )
  }
  
  ###Vector of objects of class Subsample.
  samples <- lapply( rep( "Subsample", n_samples ), new )
  
  for( i in 1:n_samples ) 
      samples[[i]] <- get_subsamples( x, trt, trt_control, training_fraction, validation_fraction, test_fraction )
  
  return( samples )
}
    
## END OF FILE
