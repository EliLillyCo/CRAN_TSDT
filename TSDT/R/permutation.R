#!usr/bin/dev/R
#################################################################################
# FILENAME   : permutation.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 09/09/2012
# DESCRIPTION: Permute data using one of three methods:
#                1. Permute treatment
#                2. Permute response
#                3. Permute response for one treament arm only
#################################################################################

# Permute a single column or multiple columns together
simple_permute <- function( x ){
  
  index <- sample( 1:NROW(x), size = NROW(x), replace = FALSE )

  if( is.data.frame( x ) )
      res <- x[index,]

  else
      res <-x[ index ]

  return( res )
  
}

# Permute response within a single treatment arm
permute_response_one_arm <- function( response, trt, permute_arm ){

  if( !is.null( names( response ) ) )
      response_names <- names( response )
  else
      response_names <- "response"


  tmp__ <- data.frame( 1:NROW(response)  )
  
  tmp__ <- cbind( tmp__, response, trt )
  names( tmp__ ) <- c( "rowname__", response_names, "trt" )
  
  
  # Subset on treatment arm to be permuted
  non_permute_arm__ <- subset( tmp__, as.character( trt ) != as.character( permute_arm ) )
  permute_arm__ <- subset( tmp__, as.character( trt ) == as.character( permute_arm ) )
  
  # Permute permute_arm
  permuted__ <- as.data.frame( simple_permute( subset( permute_arm__, select = response_names ) ) )
  names( permuted__ ) <- response_names
  permuted__$rowname__ <- permute_arm__$rowname__
  permuted__$trt <- permute_arm__$trt
  
  permuted__ <- permuted__[,c("rowname__",response_names,"trt")]
  
  # Combine results and resort to original row order
  tmp__ <- rbind( non_permute_arm__, permuted__ )
  tmp__ <- tmp__[ order( tmp__$rowname__ ), ]

  tmp__$rowname__ <- NULL
  tmp__$trt <- NULL
  names( tmp__ ) <- response_names
  
  return( tmp__ )
  
}

# Permute response treatment together (same as simple permute on
# data.frame of covariates)
permute_response_and_treatment <- function( response, trt ){

  if( !is.data.frame( response ) )
      response <- as.data.frame( response )
  
  i__ <- sample( 1:NROW(response), size = NROW(response), replace = FALSE )
  response__ <- response[ i__, ]
  trt__  <-  trt[ i__ ]
  
  
  return( list( response = response__, trt = trt__ ) )
  
}

#' @title permutation
#' @description Permute response, treatment, or response for one treatment arm only.
#' @details If a response variable is provided and treatment is not provided 
#' the response variable is permuted.
#'
#' If a treatment variable is provided and response is not provided 
#' the treatment variable is permuted.
#'
#' If a response variable and treatment variable and permute are provided the
#' response variable is permuted only for the treatment arm indicated by
#' permute_arm.
#'
#' If a response variable and treatment variable are provided, but permute_arm
#' @param response Response (or other) variable(s) to be permuted. This can be a
#' data.frame of multiple variables (e.g. a data.frame of covariates or a
#' multivariate response).
#' @param trt Treatment variable.
#' @param permute_arm reatment arm to permute.
#' @return If permuting response or treatment, returns vector of permuted response
#' or treatment. If permuting response and treatment, returns a list of
#' permuted response and treatment. 
#' @examples
#' N <- 20
#' x <- data.frame( 1:N )
#' names( x ) <- "response"
#' x$trt <- factor( c( rep( "Experimental", 9 ), rep( "Control", N - 9 ) ) )
#' x$time <- x$response
#' x$event <- 0:1
#' 
#' ## Permute treatment variable
#' ex1 <- x[,c("response","trt")]
#' ex1$permuted_trt <- permutation( trt = ex1$trt )
#' 
#' ## Permute response variable
#' ex2 <- x[,c("response","trt")]
#' ex2$permuted_response <- permutation( response = ex2$response )
#' 
#' ## Permute the response for treatment arm only
#' ex3 <- x[,c("response","trt")]
#' permuted3 <- permutation( response = ex3$response, trt = ex3$trt, permute_arm = "Experimental" )
#' names( permuted3 ) <- paste( "permuted_", names(permuted3), sep = "" )
#' ex3 <- cbind( ex3, permuted3 )
#' 
#' ## Permute response and treatment together
#' ex4 <- x[,c("response","trt")]
#' permutation_list4 <- permutation( response = ex4$response, trt = ex4$trt )
#' ex4$permuted_response <- permutation_list4$response
#' ex4$permuted_trt <- permutation_list4$trt
#' 
#' ## Permute a survival response for treatment arm only
#' ex5 <- x[,c("time","event","trt")]
#' permuted5 <- permutation( response = ex5[,c("time","event")], trt = ex5$trt,
#'                           permute_arm = "Experimental" )
#' names( permuted5 ) <- paste( "permuted_", names(permuted5), sep = "" )
#' ex5 <- cbind( ex5, permuted5 )
#' 
#' ## Permute a survival outcome and treatment together
#' ex6 <- x[,c("time","event","trt")]
#' permutation_list6 <- permutation( response = ex6[,c("time","event")], trt = ex6$trt )
#' ex6$permuted_time <- permutation_list6$response$time
#' ex6$permuted_event <- permutation_list6$response$event
#' @export
permutation <- function( response = NULL,
                         trt = NULL,
                         permute_arm = NULL ){
  
  if( is.null(trt) && is.null( response ) )
      stop( "ERROR: must provide response, trt, or both" )
  
  if( !is.null(permute_arm) && any( as.character( trt ) == as.character( permute_arm ) ) == FALSE )
      stop( paste0( "ERROR: treatment value ", permute_arm, " not found in trt" ) )
  
  # Permute response for one treatment arm only
  if( !is.null( trt ) && !is.null( response ) && !is.null( permute_arm ) )
      permuted <- permute_response_one_arm( response, trt, permute_arm )
  
  # Permute response and treatment together (preserve response/treatment relationship)
  else if( !is.null( trt ) && !is.null( response ) && is.null( permute_arm ) )
      permuted <- permute_response_and_treatment( response, trt )
  
  # Simple permute response
  else if( is.null( response ) )
      permuted <- simple_permute( trt )
  
  # Simple permute treatment
  else if( is.null( trt ) )
      permuted <- simple_permute( response )
  
  return( permuted ) 
}

## END OF FILE
