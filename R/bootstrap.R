#!usr/bin/dev/R
#################################################################################
# FILENAME   : bootstrap.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 09/09/2012
# DESCRIPTION: Functions to generate bootstrap samples and return in-bag and
#              out-of-bag data. Optionally, a statistic may be computed on each
#              bootstrap sample. Results are returned in a vector of objects of
#              class Bootstrap or BootstrapStatistic.
#################################################################################

bootstrap_bags <- function( x, trt, trt_control ){
  
  if( !is.data.frame(x) )
    x <- data.frame( x )
  
  if( trt_control == "SINGLE_ARM" ){

    n_control <- 0
    n_experimental <- NROW( x )

    control <- NULL
    experimental <- x
    
    control_inbag_rows <- NULL
    control_oob_rows   <- NULL
  }else{
    
    n_control      <- length( which(trt == trt_control) )
    n_experimental <- length( which(trt != trt_control) )

    control      <- subset( x, trt == trt_control )
    experimental <- subset( x, trt != trt_control )
  
    control_inbag_rows <- sample( 1:n_control, size = NROW( control ), replace = TRUE )
    control_oob_rows   <- which( 1:n_control %nin% control_inbag_rows )
  }
  
  experimental_inbag_rows <- sample( 1:n_experimental, size = NROW( experimental ), replace = TRUE )
  experimental_oob_rows    <- which( 1:n_experimental %nin% experimental_inbag_rows )

  experimental_inbag <- experimental[experimental_inbag_rows,]
  experimental_oob <- experimental[experimental_oob_rows,]
  
  names( experimental_inbag ) <- names( x )
  names( experimental_oob )   <- names( x )
  
  if( is.null( control ) ){
    inbag <- experimental_inbag
    oob   <- experimental_oob
  }else{

    control_inbag <- control[control_inbag_rows,]
    control_oob <- control[control_oob_rows,]
    
    names( control_inbag ) <- names( x )
    names( control_oob )   <- names( x )
    
    inbag <- rbind( experimental_inbag, control_inbag )
    
    oob   <- rbind( experimental_oob, control_oob )
  }
  
  if( !is.data.frame( inbag ) )
    inbag <- data.frame( inbag )

  if( !is.data.frame( oob ) )
    oob <- data.frame( oob )

  if( NCOL( inbag ) != NCOL( x ) )
    stop( "ERROR: wrong number of columns in in-bag sample" )

  if( NCOL( oob ) != NCOL( x ) )
    stop( "ERROR: wrong number of columns in out-of-bag sample" )


  if( trt_control != 'SINGLE_ARM' ){
    if( NROW( control_inbag ) != n_control )
        stop( 'ERROR: incorrect number of control samples in in-bags data' )
    
    if( NROW( experimental_inbag ) != n_experimental )
        stop( 'ERROR: incorrect number of experimental samples in in-bags data' )
  }

  rownames( inbag ) <- NULL
  rownames( oob ) <- NULL
  names( inbag ) <- names(x)
  names( oob ) <- names(x)
  
  return( list( inbag = inbag, oob = oob ) )
  
}


#' @title bootstrap
#' @description Generate a vector of bootstrap samples.
#' @details Each bootstrap sample will retain the in-bag and out-of-bag data.
#' Optionally, the user may specify a function to compute a statistic for each
#' in-bag and out-of-bag sample. This function may be a built-in R function
#' (e.g. mean, median, etc.) or a user-defined function (see Examples). If no
#' statistic function is provided bootstrap returns a vector of objects of class
#' \linkS4class{Bootstrap}. If a statistic function is provided bootstrap
#' returns a vector of objects of class \linkS4class{BootstrapStatistic},
#' which in addition to the in-bag and out-of-bag samples contains the name of
#' the statistic, variable on which the statistic is computed, and the numerical
#' result of the statistic for each in-bag and out-of-bag sample.
#' @seealso \linkS4class{Bootstrap}, \linkS4class{BootstrapStatistic}
#' @param x Source data to bootstrap.
#' @param trt Treatment variable. (optional)
#' @param trt_control Value for treatment control arm. Default value is 'Control'.
#' @param FUN Function to compute statistic for each bootstrap sample. (optional)
#' @param varname Name of variable in x on which to compute FUN. If x has only
#' one column varname is not needed. If x has more than one column then either
#' varname or varcol must be specified.  
#' @param varcol Column index of x on which to compute FUN. If x has only one
#' column varcol is not needed. If x has more than one column then either varname
#' or varcol must be specified.  
#' @param arglist List of additional arguments to pass to FUN.
#' @param n_samples Number of bootstrap samples to generate.
#' @return  If FUN is NULL returns a vector of objects of class
#' \linkS4class{Bootstrap}. If FUN is non-NULL returns a vector of objects
#' of class \linkS4class{BootstrapStatistic}
#' @export
#' @examples
#' ## Generate example data frame containing response and treatment
#' N <- 20
#' x <- data.frame( runif( N ) )
#' names( x ) <- "response"
#' x$treatment <- factor( sample( c("Control","Experimental"), size = N,
#'                        prob = c(0.8,0.2), replace = TRUE ) )
#' 
#' ## Generate two bootstrap samples without regard to treatment
#' ex1 <- bootstrap( x, n_samples = 2 )
#' 
#' ## Generate two bootstrap samples stratified by treatment
#' ex2 <- bootstrap( x, trt = x$treatment, trt_control = "Control", n_samples = 2 )
#' 
#' ## For each bootstrap sample compute a statistic on the in-bag and out-of-bag data
#' ex3 <- bootstrap( x, FUN = mean, varname = "response", n_samples = 2 )
#' 
#' ## Specify a user-defined function that takes a numeric vector input and
#' ## returns a numeric result
#' sort_and_rank <- function( z, rank ){
#'   z <- sort( z )
#'   return( z[rank] )
#' }
#' 
#' ex4 <- bootstrap( x, FUN = sort_and_rank, arglist = list( rank = 1 ),
#'                   varname = "response", n_samples = 2 )
bootstrap <- function( x,
                       trt = NULL,
                       trt_control = "Control",
                       FUN = NULL,
                       varname = NULL,
                       varcol = NULL,
                       arglist = NULL,
                       n_samples = 1  ){
  
  if( !is.data.frame( x ) )
      x <- data.frame( x )
  
  if( !is.null( varname ) && ( varname %nin% names( x ) ) )
      stop( paste0( "ERROR: variable ", varname, " not found" ) )
  
  if( is.null( trt ) ){
    trt <- rep( "SINGLE_ARM", times = NROW( x ) )
    trt_control <- "SINGLE_ARM"
  }
  
  samples <- lapply( rep( "Bootstrap", n_samples ), new )
  
  for( i in 1:n_samples ){
 
    bags <- bootstrap_bags( x = x, trt = trt, trt_control = trt_control )
    
    if( is.null( FUN ) ){
      next_sample <- new( "Bootstrap", inbag = bags$inbag, oob = bags$oob )
      samples[[i]] <- next_sample
    }
    
    else{
      FUN_name <- deparse( substitute( FUN ) )
      
      if( is.null( varname ) & is.null( varcol ) &  NCOL(x) > 1 ){
        stop( "ERROR: If FUN is non-null and NCOL(x) > 1 then either varname or varcol must be non-null" )
      }else if( is.null( varname ) & is.null( varcol ) &  NCOL(x) == 1 ){
        varcol <- 1
      }
      
      if( !is.null( varname ) ){
        inbagX <- paste0( "bags$inbag$", varname )
        oobX   <- paste0( "bags$oob$", varname )
        variable <- varname
      }else if( !is.null( varcol ) ){
        inbagX <- paste0( "bags$inbag[,", varcol, "]" )
        oobX   <- paste0( "bags$oob[,", varcol, "]" )
        variable <- names(x)[varcol]
      }
     
      if( is.null(arglist) ){
        inbag_stat <- do.call( what = FUN, args = list( eval( parse( text = inbagX ) ) ) )
        oob_stat   <- do.call( what = FUN, args = list( eval( parse( text = oobX ) ) ) )
        
        next_sample <- new( "BootstrapStatistic", inbag = bags$inbag, oob = bags$oob,
                           statname = FUN_name, arglist = list(NULL), variable = variable,
                           inbag_stat = inbag_stat, oob_stat = oob_stat )
      }
      else{
        inbag_stat <- do.call( FUN, args = c( list( eval(parse( text = inbagX  ) ) ), arglist ) )
        oob_stat   <- do.call( FUN, args = c( list( eval(parse( text = oobX ) ) ), arglist ) )
        
        next_sample <- new( "BootstrapStatistic", inbag = bags$inbag, oob = bags$oob,
                           statname = FUN_name, arglist = arglist, variable = variable,
                           inbag_stat = inbag_stat, oob_stat = oob_stat )
      }
      samples[[i]] <- next_sample
    }
    rm( next_sample )
  }
  return( samples )
}

## END OF FILE
