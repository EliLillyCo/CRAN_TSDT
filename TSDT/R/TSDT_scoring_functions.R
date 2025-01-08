#!usr/bin/dev/R
#################################################################################
# FILENAME   : TSDT_scoring_functions.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 09/24/2013
# DESCRIPTION: 
#################################################################################

scoring_function_wrapper <- function( scoring_function_name, data, scoring_function_parameters = NULL ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  result <- NULL;rm( result )
  
  scoring_function_text <- paste0( 'result <- ', scoring_function_name, '( data' )

  if( !is.null( scoring_function_parameters ) )
      scoring_function_text <- paste0( scoring_function_text, ', scoring_function_parameters' )

  scoring_function_text <- paste0( scoring_function_text, ' )' )
  
  eval( parse( text = scoring_function_text ) )
                   
  return( result )
}

#' @title mean_response
#' @description Compute the mean response.
#' @details This function will compute the mean of the response variable. If
#' a value for trt_arm is provided the mean in that treatment arm only will be
#' computed (and the trt variable must also be provided), otherwise the mean
#' for all data passed to the function will be computed.
#' @seealso \link{TSDT}, \link{treatment_effect}
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control
#' parameters
#' @return The mean of the provided response variable.
#' @examples
#' N <- 50
#' 
#' data <- data.frame( continuous_response = numeric(N),
#'                    trt = character(N) )
#' 
#' data$continuous_response <- runif( min = 0, max = 20, n = N )
#' data$trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6), replace = TRUE )
#' 
#' ## Compute mean response for all data
#' mean_response( data, scoring_function_parameters = list( y_var = 'continuous_response' ) )
#' mean( data$continuous_response ) # Function return value should match this value
#' 
#' ## Compute mean response for Experimental treatment arm only
#' scoring_function_parameters <- list( y_var = 'continuous_response', trt_arm = 'Experimental' )
#' mean_response( data, scoring_function_parameters = scoring_function_parameters )
#' # Function return value should match this value
#' mean( data$continuous_response[ data$trt == 'Experimental' ] )
#' @export
mean_response <- function( data,
                           scoring_function_parameters = NULL ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  trt_arm <- NULL;rm( trt_arm )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )

  response <- get_y( data, scoring_function_parameters )
  trt <- get_trt( data, scoring_function_parameters )

  if( !exists( "trt_arm" ) ){
    if( "WEIGHTS__" %nin% names( data ) )
        return( mean( response, na.rm = TRUE ) )
    else
        return( weighted.mean( x = response, w = data$WEIGHTS__, na.rm = TRUE ) )
  }else{
    if( !any( trt == trt_arm ) )
        stop( "ERROR: trt_arm value not found" )
        
    if( !exists( "trt_arm" ) ){
      if( "WEIGHTS__" %nin% names( data ) )
          return( mean( response[ trt == trt_arm ], na.rm = TRUE ) )
      else
          return( weighted.mean( response[ trt == trt_arm ], w = data$WEIGHTS__[ trt == trt_arm ], na.rm = TRUE ) )
    }
  }
}


#' @title quantile_response
#' @description Return the specified quantile of the response distribution.
#' @details This function returns the response quantiles associated with a
#' specified percentile. The default behavior is to return the median -- i.e.
#' 50th-percentile.
#' @seealso \link{TSDT}, \link{diff_quantile_response}, \link{quantile}
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control
#' parameters
#' @return A quantile of the response variable.
#' @examples
#' ## Generate example data containing response and treatment
#' N <- 100
#' y = runif( min = 0, max = 20, n = N )
#' df <- as.data.frame( y )
#' names( df )  <- "y"
#' df$trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6),
#'                   replace = TRUE )
#' 
#' ## Default behavior is to return the median
#' quantile_response( df )
#' median( df$y ) # should match previous result from quantile_response
#' 
#' ## Get Q1 response
#' quantile_response( df, scoring_function_parameters = list( percentile = 0.25 ) )
#' quantile( df$y, 0.25 ) # should match previous result from quantile_response
#' 
#' ## Get max response
#' quantile_response( df, scoring_function_parameters = list( percentile = 1 ) )
#' max( df$y ) # should match previous result from quantile_response
#' @export
quantile_response <- function( data,
                               scoring_function_parameters = NULL ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  trt_arm  <- NULL;rm( trt_arm )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  
  if( !exists( "percentile" ) )
      percentile <- 0.50
  
  response <- get_y( data, scoring_function_parameters )
  
  trt <- get_trt( data, scoring_function_parameters )
  
  if( !exists( "trt_arm" ) ){
    return( quantile( response, probs = percentile ) )
  }
  else{
    if( !any( trt == trt_arm ) )
        stop( "ERROR: trt_arm value not found" )
    
    return( quantile( response[ trt == trt_arm ], probs = percentile ) )
  }  
}


#' @title diff_quantile_response
#' @description Return the difference across treatment arms of a specified
#' response quantile
#' @details This function returns the difference across treatment arms of the
#' response quantile associated with a specified percentile. The default behavior
#' is to return the difference in medians.
#' @seealso \link{TSDT}, \link{quantile_response}, \link{quantile}
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control parameters
#' @return A difference of response quantiles across treatment arms
#' @examples
#' ## Generate example data containing response and treatment
#' N <- 100
#' y = runif( min = 0, max = 20, n = N )
#' df <- as.data.frame( y )
#' names( df )  <- "y"
#' df$trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6),
#'                   replace = TRUE )
#' 
#' ## Default behavior is to return the median
#' diff_quantile_response( df )
#' 
#' # should match previous result from quantile_response
#' median( df$y[df$trt!='Control'] ) - median( df$y[df$trt=='Control'] )
#' 
#' ## Get Q1 response
#' diff_quantile_response( df, scoring_function_parameters = list( percentile = 0.25 ) )
#'
#' # should match previous result from quantile_response
#' quantile( df$y[df$trt!='Control'], 0.25 ) - quantile( df$y[df$trt=='Control'], 0.25 )
#' 
#' ## Get max response
#' diff_quantile_response( df, scoring_function_parameters = list( percentile = 1 ) )
#'
#' # should match previous result from quantile_response
#' max( df$y[df$trt!='Control'] ) -  max( df$y[df$trt=='Control'] )
#' @export
diff_quantile_response <- function( data,
                                    scoring_function_parameters = NULL ){
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  
  if( !exists( "trt_control" ) )
      trt_control <- 'Control'
  
  trt <- get_trt( data, scoring_function_parameters )
  
  if( !any( trt == trt_control ) )
    stop( paste0( 'ERROR: trt_control value (', trt_control, ') not found in trt variable' ) )
  
  scoring_function_parameters$trt_arm <- trt_control
  quantile_control <- quantile_response( data, scoring_function_parameters )
  
  scoring_function_parameters$trt_arm <- get_experimental_trt_arm( trt, trt_control )
  quantile_trt <- quantile_response( data, scoring_function_parameters )
  
  return( quantile_trt - quantile_control )
}


#' @title treatment_effect
#' @description Compute treatment effect as mean( treatment response ) - mean( control response )
#' @details This function will compute the treatment for the response. The
#' treatment effect is computed as the difference in means between the non-control
#' treatment arm and the control treatment arm. The user must provide the
#' treatment variable as well as the control value.
#' @seealso \link{TSDT}, \link{mean_response}
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control parameters
#' @return The difference in mean response across treatment arms.
#' @examples
#' N <- 100
#' 
#' df <- data.frame( continuous_response = numeric(N),
#'                   trt = integer(N) )
#' 
#' df$continuous_response <- runif( min = 0, max = 20, n = N )
#' df$trt <- sample( c(0,1), size = N, prob = c(0.4,0.6), replace = TRUE )
#' 
#' # Compute the treatment effect
#' treatment_effect( df, list( y_var = 'continuous_response', trt_control = 0 ) )
#' 
#' # Function return value should match this value
#' mean( df$continuous_response[df$trt == 1] ) - mean( df$continuous_response[df$trt == 0] )
#' @export
treatment_effect <- function( data,
                              scoring_function_parameters = NULL ){
    
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  
  response <- get_y( data, scoring_function_parameters )
  trt <- get_trt( data, scoring_function_parameters )
  
  if( !exists( "trt_control" ) )
      trt_control <- 'Control'
  
  if( !any( trt == trt_control ) )
      stop( paste0( 'ERROR: trt_control value (', trt_control, ') not found in trt variable' ) )
  if( "WEIGHTS__" %nin% names( data ) ){
    mean_trt_response <- mean( response[ trt != trt_control ], na.rm = TRUE )
    mean_control_response <- mean( response[ trt == trt_control ], na.rm = TRUE )
  }else{
    mean_trt_response <- weighted.mean( response[ trt != trt_control ], w = data$WEIGHTS__[ trt != trt_control ], na.rm = TRUE )
    mean_control_response <- weighted.mean( response[ trt == trt_control ], w = data$WEIGHTS__[ trt == trt_control ], na.rm = TRUE )
  }
  return( mean_trt_response - mean_control_response )  
}


## #' @title desirable_response_proportion
## #' @description Compute proportion of subjects that have a desirable (binary) response
## #' @details This function will compute the proportion of subjects for which
## #' the response indicates a desirable response. The response proportion can
## #' be computed for a single treatment arm only (if valued for trt_arm and trt
## #' are provided) or for all data passed to the function.
## #' @seealso link{TSDT}, \link{binary_transform}
## #' @param data data.frame containing response data
## #' @param scoring_function_parameters named list of scoring function control parameters
## #' @return Proportion of (binary) response values that have a desirable value.
## #' The desirable value is 1 if desirable_response = 'increasing' and the
## #' desirable value is 0 if desirable_response = 'decreasing'.
## #' @examples
## #' N <- 50
## #' 
## #' data <- data.frame( binary_response = numeric(N) )
## #' 
## #' data$binary_response <- sample( c(0,1), size = N, prob = c(0.3,0.7), replace = TRUE )
## #' 
## #' ## Compute desirable response proportion with increasing desirable response
## #' ## (i.e. larger response value is better)
## #' 
## #' desirable_response_proportion( data, list( y_var = 'binary_response', desirable_response = 'increasing' ) )
## #' mean( data$binary_response ) # Function return value should match this value
## #'
## #' 
## #' ## Compute desirable response proportion for Experimental treatment arm only
## #' ## with decreasing desirable response (i.e. smaller response value is better).
## #'
## #' data$trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6), replace = TRUE )
## #' 
## #' desirable_response_proportion( data, list( y_var = 'binary_response', desirable_response = 'decreasing', trt_control = 'Control' ) )
## #'
## #' prop.table( table( subset( data, trt != 'Control', select = 'binary_response', drop = TRUE ) ) )
## #' 
## #' @export
## desirable_response_proportion <- function( data,
##                                            scoring_function_parameters = NULL ){

##   if( !is.null( scoring_function_parameters ) )
##       unpack_args( scoring_function_parameters )

##   response <- binary_transform( get_y( data, scoring_function_parameters ) )
##   trt <- get_trt( data, scoring_function_parameters )
  
##   if( desirable_response %nin% c("increasing","decreasing") )
##       stop( "ERROR: desirable_response must be in {increasing,decreasing}" )      

##   if( !is.null( trt ) && !exists( "trt_control" ) )
##       stop( "ERROR: if a treatment variable is present in data a trt_control value must be provided in scoring_function_parameters" )
  
##   # For computing the proportion when no treatment variable is provided
##   #if( !exists( "trt_var" ) ){
##   if( is.null( trt ) ){
##     # Always count score as proportion of response values equal to 1.
##     # Do not adjust this definition for desirable_response. This will
##     # be handled by the superior_subgrous() function.

##     N <- length( response )
    
##     if( desirable_response == "increasing" ){
##       proportion <- length( subset( response, response == 1 ) )/N
##     }else{
##       proportion <- length( subset( response, response == 0 ) )/N
##     }
##   }
##   # For computing the proportion only on the treatment arm
##   else{
    
##     response <- subset( response, trt != trt_control )
##     N <- length( response )
    
##     if( desirable_response == "increasing" ){
##       response <- subset( response, response == 1 )
##     }else{ #decreasing
##       response <- subset( response, response == 0 )
##     }
##     proportion <- length( response )/N 
##   } 
##   return( proportion )
## }

#################################################################################
# Scoring functions for a survival outcome: difference in median (or other      #
# quantile) survival time.                                                      #
################################################################################# 
#' @title survival_time_quantile
#' @description Computes the quantile of a survival function.
#' @details Computes the quantile of a survival function. The user specifies the
#' percentile associated with the desired quantile in
#' scoring_function_parameters. The default is percentile = 0.50, which returns
#' the median survival. A user may also specify a value for the trt_arm parameter
#' in scoring_function_parameters to compute the survival quantile for only one
#' arm.     
#' @seealso \link{TSDT}, \link{diff_survival_time_quantile},
#' \link[survival]{Surv}, \link[survival]{coxph}, \link[survival]{survfit},
#' \link[survival]{survreg}, \link[survival]{quantile.survfit},
#' \link[survival]{predict.coxph}, \link[survival]{predict.survreg}
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control parameters
#' @return A quantile of the response survival time.
#' @examples
#' N <- 200
#' time <- runif( min = 0, max = 20, n = N )
#' event <- sample( c(0,1), size = N, prob = c(0.2,0.8), replace = TRUE )
#' trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6), replace = TRUE )
#' df <- data.frame( y = survival::Surv( time, event ), trt = trt )
#' 
#' ## Compute median survival time in Experimental treatment arm.
#' ex1 <- survival_time_quantile( data = df,
#'                                scoring_function_parameters = list( trt_var = "trt",
#'                                trt_arm = "Experimental",
#'                                percentile = 0.50 ) )
#' 
#' ## Compute Q1 survival time for all data. It is necessary here to explicitly
#' ## specify trt = NULL because a variable called trt exists in df. The default
#' ## behavior is to use this variable as the treatment variable. To override
#' ## the default behavior trt = NULL is included in scoring_function_parameters.
#' ex2 <- survival_time_quantile( data = df,
#'                                scoring_function_parameters = list( trt = NULL, percentile = 0.25 ) )
#' @export
survival_time_quantile <- function( data,
                                    scoring_function_parameters = NULL ){ 
  
  requireNamespace( "survival", quietly = TRUE )
  
  ## Create NULL placeholders to prevent NOTE in R CMD check
  trt_control <- NULL;rm( trt_control )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  
  # Cox proportional hazard
  COXPH_MODELS <- "coxph"
  
  # Non-parametric and semi-parametric models
  SURVFIT_MODELS <- c("kaplan-meier","fleming-harrington","fh2")
  
  # Parametric models
  SURVREG_MODELS <- c("weibull","exponential","gaussian","logistic",
                      "lognormal","loglogistic")

  MODEL_SET <- paste( c( COXPH_MODELS, SURVFIT_MODELS, SURVREG_MODELS ), sep = "", collapse = "," )
  
  # Validate input parameters and define default values
  if( !exists( "survival_model" ) )
      survival_model <- "kaplan-meier" # Compute Kaplan-Meier estimator by default

  if( !exists( "percentile" ) )        # Compute median survival by default
      percentile <- 0.50
  
  if( survival_model %nin% c( COXPH_MODELS, SURVFIT_MODELS, SURVREG_MODELS ) )
      stop( paste0( "ERROR: survival_model is ", survival_model, ". It must be one of {", MODEL_SET, "}." ) )
  
  # The class of response is changed from Surv to matrix when retrieved with
  # get_y(). Cast it back to Surv by wrapping the return value in Surv().
  response <- get_y( data, scoring_function_parameters )
  
  if( !is(response, "Surv") )
      response <- survival::Surv( response )
  
  trt <- get_trt( data, scoring_function_parameters )

  # If no trt_arm is provided for which arm to fit the survival model on then
  # assume the experimental arm.
  if( !is.null( trt ) && !exists("trt_arm") )
    trt_arm <- as.character( trt[ trt != trt_control ][[1]] )
  
  # Construct Surv object for the survival model formula
  if( !is.null( trt ) ){
    # It is neccessary to re-create the Surv object after subsetting
    surv_time <- eval( parse( text = paste0( 'response[ as.character( trt ) == "', trt_arm, '","time"]' ) ) )
    surv_status <- eval( parse( text = paste0( 'response[ as.character( trt ) == "', trt_arm, '","status"]' ) ) )
    response <- survival::Surv( surv_time, surv_status )
  }

  surv_formula_y <- 'response ~ '
  surv_formula_X <- '1'
  SURV_FORMULA <- as.formula( paste0( surv_formula_y, surv_formula_X ) )
  
  # Compute survival model
  if( survival_model %in% COXPH_MODELS ){
    surv <- survival::survfit( coxph( formula = SURV_FORMULA ), se.fit = FALSE, control = survival::coxph.control( iter.max = 1E4 ) )
  }else if( survival_model %in% SURVFIT_MODELS ){
    surv <- survival::survfit( formula = SURV_FORMULA, se.fit = FALSE, type = survival_model )
  }else if( survival_model %in% SURVREG_MODELS ){
    surv <- survival::survreg( formula = SURV_FORMULA, dist = survival_model, control = survival::survreg.control( iter.max = 1E3 ) )
  }
  
  # Get survival function quantile for specified percentile
  if( survival_model %in% c( COXPH_MODELS, SURVFIT_MODELS ) ){
    
    quantile__ <- quantile( surv, probs = percentile )
    
    survival_time_quantile <- NA
    
    if( !is(quantile__, "logical") ){
      
      if( is( quantile__, "numeric" ) ){
        survival_time_quantile <- quantile__
      }else{
        survival_time_quantile <- quantile__$quantile[[1]]
      }
    }# END if quantile__ != logical
  }
    
  else if( survival_model %in% c( SURVREG_MODELS ) ){
    survival_time_quantile <- predict( surv, type = "quantile", p = percentile )[[1]]
  }
  
  # Return quantile
  return( survival_time_quantile )
}

#' @title diff_survival_time_quantile 
#' @description Computes the difference in the quantile of a survival function
#' across treatment groups.
#' @details Computes the survival function quantile for the treatment and
#' control arms and returns the difference.
#' @seealso \link{TSDT}, \link{survival_time_quantile}, \link[survival]{Surv},
#' \link[survival]{coxph}, \link[survival]{survfit}, \link[survival]{survreg},
#' \link[survival]{predict.coxph}, \link[survival]{predict.survreg} 
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control parameters
#' @return A difference in a survival time quantile across treatment arms.
#' @examples
#' requireNamespace( "survival", quiet = TRUE )
#' N <- 200
#' df <- data.frame( y = survival::Surv( runif( min = 0, max = 20, n = N ),
#'                             sample( c(0,1), size = N, prob = c(0.2,0.8), replace = TRUE ) ),
#'                   trt = sample( c('Control','Experimental'), size = N,
#'                                 prob = c(0.4,0.6), replace = TRUE ) )
#' 
#' ## Compute difference in median survival time between Experimental arm and
#' ## Control arm.  It is not actually necessary to provide the value for the
#' ## time_var, trt_var, trt_control, and percentile parameters because these
#' ## values are all equal to their default values. The value are explicitly
#' ## provided here simply for clarity.
#' ex1 <- diff_survival_time_quantile( data = df,
#'                                     scoring_function_parameters = list( trt_var = "trt",
#'                                     trt_control = "Control",
#'                                     percentile = 0.50 ) )
#' 
#' ## Compute difference in Q1 survival time. In this example the default value
#' ## for all scoring function parameters are used except percentile, which here
#' ## takes the value 0.25.
#' ex2 <- diff_survival_time_quantile( data = df,
#'                                     scoring_function_parameters = list( percentile = 0.25 ) )
#' @export
diff_survival_time_quantile <- function( data,
                                         scoring_function_parameters = NULL ){
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  
  # Get treatment variable and control and experimental treatment values
  trt <- get_trt( data, scoring_function_parameters )
  
  if( is.null( trt ) )
      stop( "ERROR: trt variable not specified" )
  
  if( !exists( "trt_control" ) )
      trt_control <- 'Control'
  
  experimental_trt_value <- as.character( trt[ trt != trt_control ][[1]] )
  
  # Get survival quantile for control arm
  scoring_function_parameters$trt_arm <- trt_control
  survival_time_control <- survival_time_quantile( data, scoring_function_parameters )
  
  # Get survival quantile for experimental arm
  scoring_function_parameters$trt_arm <- experimental_trt_value
  survival_time_trt <- survival_time_quantile( data, scoring_function_parameters )
  
  difference <- survival_time_trt - survival_time_control
  
  return( difference )
}


#' @title mean_deviance_residuals
#' @description Computes the mean of the deviance residuals from a survival model
#' @details Computes the mean of the deviance residuals from a survival model.
#' The deviance residual at time t is computed as the observed number of events
#' at time t minus the expected number of events at time t (see Therneau, et. al.
#' linked below). The expected number of events is the number of events predicted
#' by the survival model. If the event under study is an undesirable event (as
#' would likely be the case in a clinical context), then a smaller value for the
#' deviance residual is desirable -- i.e. it is desirable to observe fewer events
#' than expected from the survival model. In this case the appropriate value for
#' desirable_response in TSDT is desirable_response = 'decreasing'. If the event
#' under study is desirable then the appropriate value for desirable_response is
#' desirable_response = 'increasing'. It is assumed that most survival models are
#' modeling an undesirable event. Therefore, when the user specifies
#' mean_deviance_residual or diff_mean_deviance_residual, the default value for
#' desirable_repsonse is changed to 'decreasing', unless the user explicitly
#' provides desirable_response = 'increasing'. Note this differs from all other
#' TSDT configurations, for which the default value for desirable_response is
#' desirable_response = 'increasing'.
#' @seealso \link{TSDT}, \link[survival]{Surv}, \link[survival]{coxph},
#' \link[survival]{survfit}
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control parameters
#' @return Mean of deviance residuals
#' @references
#' Therneau, T.M.,  Grambsch, P.M., and Fleming, T.R. (1990).  Martingale-based
#' residuals for survival models.  Biometrika, 77(1), 147-160.
#' \doi{10.1093/biomet/77.1.147}
#' @export
mean_deviance_residuals <- function( data,
                                     scoring_function_parameters = NULL ){

  requireNamespace( "survival", quietly = TRUE )

  ## Create NULL placeholders to prevent NOTE in R CMD check
  trt_control <- NULL;rm( trt_control )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )

  # Cox proportional hazard
  COXPH_MODELS <- "coxph"

  # Parametric models
  SURVREG_MODELS <- c( "weibull", "exponential", "gaussian", "logistic",
                       "lognormal", "loglogistic")

  MODEL_SET <- paste( c( COXPH_MODELS, SURVREG_MODELS ), sep = "", collapse = "," )
  
  # Validate input parameters and define default values
  if( !exists( "survival_model" ) )
      survival_model <- "coxph"    # Fit Cox Proportional Hazards model by default
  
  if( survival_model %nin% c( COXPH_MODELS, SURVREG_MODELS ) )
      stop( paste0( "ERROR: survival_model is ", survival_model,". To use ",
                    "mean_deviance_residuals as the scoring function survival_model must ",
                    "be one of {", MODEL_SET, "}." ) )

  covariates <- get_covariates( data, scoring_function_parameters )
  # The class of response is changed from Surv to matrix when retrieved with
  # get_y(). Cast it back to Surv by wrapping the return value in Surv().
  response <- get_y( data, scoring_function_parameters )
  if( !is(response,"Surv") )
      response <- survival::Surv( response )
  
  trt <- get_trt( data, scoring_function_parameters )

  # If no trt_arm is provided for which arm to fit the survival model on then
  # assume the experimental arm.
  if( !is.null( trt ) && !exists("trt_arm") )
    trt_arm <- as.character( trt[ trt != trt_control ][[1]] )
  
  # Construct Surv object for the survival model formula
  if( is.null( trt ) ){
    surv_formula_y <- 'response ~ '
  }else{
    surv_formula_y <- 'response[trt==trt_arm] ~ '
  }
  
  COVARS <- names( covariates )
    
  df <- as.data.frame( covariates )
  names( df ) <- COVARS

  surv_formula_X <- '1'
  
  SURV_FORMULA <- as.formula( paste0( surv_formula_y, surv_formula_X ) )

  # Compute survival model
  if( survival_model %in% COXPH_MODELS )
    surv <- survival::coxph( formula = SURV_FORMULA, control = survival::coxph.control( iter.max = 1E3 )  )
     
  else if( survival_model %in% SURVREG_MODELS )
    surv <- survival::survreg( formula = SURV_FORMULA, dist = survival_model, control = survival::survreg.control( iter.max = 1E3 ) )
  
  # Get deviance residuals
  deviance_residuals <- residuals( surv, type  = "deviance" )
  
  # Return mean of deviance residuals
  return( mean( deviance_residuals, na.rm = TRUE ) )
}

#' @title diff_mean_deviance_residuals
#' @description Computes the difference in the mean of deviance residuals
#' function across treatment groups.
#' @details The deviance residual is the observed number of
#' events at time t minus the expected number of events at time t. See
#' documentation for mean_deviance_residuals (linked below) for more details.
#' A smaller value for the deviance residual is preferred when the event under
#' study is an undesirable event -- i.e. it is preferred to observe fewer events
#' than predicted by the survival model. A two-arm TSDT model computes the mean
#' deviance residual in the treatment arm minus the mean deviance residual in the
#' control arm. The treatment arm is superior to the control arm when the mean
#' deviance residual in the treatment arm is less than the mean deviance residual
#' in the control arm. Thus, the appropriate value for desirable_response is
#' desirable_response = 'decreasing'. If the event under study is a desirable
#' event the appropriate value for desirable_response is desirable_response =
#' 'increasing'. It is assumed most survival models will model an undesirable
#' event, so the default value for desirable_response when the scoring_function
#' is diff_mean_deviance_residuals is desirable_response = 'decreasing'. Note
#' this differs from all other TSDT configurations, for which the default value
#' for desirable_response is desirable_response = 'increasing'.
#' @seealso \link{mean_deviance_residuals}, \link[survival]{Surv},
#' \link[survival]{coxph}, \link[survival]{survreg},
#' \link[survival]{residuals.coxph},
#' \link[survival]{residuals.survreg},
#' \link{TSDT}
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control parameters
#' @return Difference in mean deviance residuals across treatment arms.
#' @export
diff_mean_deviance_residuals <- function( data,
                                         scoring_function_parameters = NULL ){
  requireNamespace( "survival", quietly = TRUE )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
          
  # Get treatment variable and control and experimental treatment values
  trt <- get_trt( data, scoring_function_parameters )
  
  if( is.null( trt ) )
      stop( "ERROR: trt variable not specified" )

  if( !exists( "trt_control" ) )
    trt_control <- 'Control'
  
  experimental_trt_value <- as.character( trt[ trt != trt_control ][[1]] )

  # Get survival quantile for control arm
  scoring_function_parameters$trt_arm <- trt_control
  mean_deviance_residuals_control <- mean_deviance_residuals( data, scoring_function_parameters )

  # Get survival quantile for experimental arm
  scoring_function_parameters$trt_arm <- experimental_trt_value
  mean_deviance_residuals_trt <- mean_deviance_residuals( data, scoring_function_parameters )
  
  difference <- mean_deviance_residuals_trt - mean_deviance_residuals_control
  
  return( difference )
}

#' @title diff_restricted_mean_survival_time
#' @description Computes the difference in restricted mean survival time across
#' treatment arms.
#' @details Computes the restricted mean survival time for the treatment and
#' control arms and returns the difference.
#' @seealso \link{TSDT}, \link[survival]{Surv}, \link[survRM2]{rmst2}
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control parameters
#' @return Difference in restricted mean survival time across treatment arms.
#' @export
diff_restricted_mean_survival_time <- function( data,
                                               scoring_function_parameters = NULL ){
  
  requireNamespace( "survRM2", quietly = TRUE )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  # The class of response is changed from Surv to matrix when retrieved with
  # get_y(). Cast it back to Surv by wrapping the return value in Surv().
  response <- get_y( data, scoring_function_parameters )
  
  if( !is(response, "Surv") )
      response <- survival::Surv( response )
  
  # Get treatment variable and control and experimental treatment values
  trt <- get_trt( data, scoring_function_parameters )
  
  if( is.null( trt ) )
      stop( "ERROR: trt variable not specified" )

  if( !exists( "trt_control" ) )
    trt_control <- 'Control'
  
  experimental_trt_value <- as.character( trt[ trt != trt_control ][[1]] )

  rmst <- rmst2( time = response[,1],
                 status = response[,2],
                 arm = ifelse( trt != trt_control, 1, 0 ) )

  difference <- rmst$unadjusted.result['RMST (arm=1)-(arm=0)','Est.']
  
  return( difference )
}

#' @title hazard_ratio
#' @description Computes the hazard ratio across treatment arms using a CoxPH
#' model.
#' @seealso \link{TSDT}, \link[survival]{Surv},  \link[survival]{coxph},
#' @param data data.frame containing response data
#' @param scoring_function_parameters named list of scoring function control parameters
#' @return Hazard ratio across treatment arms.
#' @export
hazard_ratio <- function( data,
                         scoring_function_parameters = NULL ){

  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  # The class of response is changed from Surv to matrix when retrieved with
  # get_y(). Cast it back to Surv by wrapping the return value in Surv().
  response <- get_y( data, scoring_function_parameters )
  
  if( !is(response, "Surv") )
      response <- survival::Surv( response )
  
  # Get treatment variable and control and experimental treatment values
  trt <- get_trt( data, scoring_function_parameters )
  
  if( is.null( trt ) )
      stop( "ERROR: trt variable not specified" )
  
  if( !exists( "trt_control" ) )
      trt_control <- 'Control'

  ## If the coxph model successfully completes return the hazard ratio value.
  ## If the model returns a warning (e.g. if the model does not converge) return
  ## NA. If the model returns an error pass the error to the user.
  hazard_ratio__ <- tryCatch( exp( survival::coxph( response ~ trt )$coefficients[1] ),
                              warning = function(w){return(NA)},
                              error = function(e){stop(e)})
  
  return( hazard_ratio__ )
}

## END OF FILE
