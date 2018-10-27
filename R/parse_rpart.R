#!usr/bin/dev/R
#################################################################################
# FILENAME   : parse_rpart.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 10/23/2012
# DESCRIPTION: Functions for constructing a tree with rpart and parsing the
#              output to extract subgroups.
#################################################################################

#' @title rpart_nodes
#' @description Extract node information from an rpart.object.
#' @details Information about nodes and splits returned in an rpart.object
#' is contained in strings printed to the console. This function parses those
#' strings and populates a data.frame.
#' @seealso \link[rpart]{rpart}, \link[rpart]{rpart.object}
#' @param tree An rpart.object returned from call to rpart().
#' @return A data.frame containing the nodes of a parsed tree.
#' @examples
#' requireNamespace( "rpart", quietly = TRUE )
#' ## Generate example data containing response, treatment, and covariates
#' N <- 50
#' continuous_response = runif( min = 0, max = 20, n = N )
#' binary_response <- sample( c('A','B'), size = N, prob = c(0.5,0.5),
#'                            replace = TRUE )
#' trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6),
#'                replace = TRUE )
#' X1 <- runif( N, min = 0, max = 1 )
#' X2 <- runif( N, min = 0, max = 1 )
#' X3 <- sample( c(0,1), size = N, prob = c(0.2,0.8), replace = TRUE )
#' X4 <- sample( c('A','B','C'), size = N, prob = c(0.6,0.3,0.1), replace = TRUE )
#' 
#' ## Fit an rpart model with continuous response (i.e. regression)
#' fit1 <- rpart::rpart( continuous_response ~ trt + X1 + X2 + X3 + X4 )
#' fit1
#' 
#' ## Parse the results into a new data.frame
#' ex1 <- rpart_nodes( fit1 )
#' ex1
#' 
#' ## Fit an rpart model with binary response (i.e. classification)
#' fit2 <- rpart::rpart( binary_response ~ trt + X1 + X2 + X3 + X4 )
#' fit2
#' @export
rpart_nodes <- function( tree ){
  
  rpart_output <- capture.output( print(tree) )
  rpart_output <- gsub( pattern = "^[ ]*", replacement = "", rpart_output )
  rpart_output <- gsub( pattern = "< ", replacement = "<", rpart_output )
  rpart_output <- gsub( pattern = "> ", replacement = ">", rpart_output )
  
  classification <- ifelse( grepl( "(yprob)", rpart_output[[3]] ), TRUE, FALSE )
  
  rpart_output <- rpart_output[ grepl( "^[0-9]+)", rpart_output ) ]
  
  ROWS <- length( rpart_output )
  
  if( classification ){
    nodes <- data.frame( NodeID = integer( ROWS ),
                         Split = character( ROWS ),
                         n = integer( ROWS ),
                         deviance = numeric( ROWS ),
                         yval = character( ROWS ),
                         ClassificationProbability = numeric( ROWS ),
                         Leaf = logical( ROWS ), stringsAsFactors = FALSE )
  }
  else{
    nodes <- data.frame( NodeID = integer( ROWS ),
                         Split = character( ROWS ),
                         n = integer( ROWS ),
                         deviance = numeric( ROWS ),
                         yval = numeric( ROWS ),
                         Leaf = logical( ROWS ), stringsAsFactors = FALSE )
  }
  
  for( i in 1:ROWS ){
    
    rpart_output[[i]] <- gsub( pattern = "[\\(\\)]", replacement = "", rpart_output[[i]] )
    
    tokens <- strsplit( rpart_output[[i]], split = " " )[[1]]
    
    tokens <- tokens[ tokens != "" ]
    
    nodes$NodeID[i] <- as.integer( tokens[1] ) 
    nodes$Split[i] <- tokens[2]
    nodes$n[i] <- as.integer( tokens[3] )
    nodes$deviance[i] <- as.numeric( tokens[4] )
    
    if( classification ){
      nodes$yval[i] <- tokens[5]
      nodes$ClassificationProbability[i] <- as.numeric( tokens[6] )
      nodes$Leaf[i] <- ifelse( length(tokens) == 8, TRUE, FALSE )
    }else{
      nodes$yval[i] <- as.numeric( tokens[5] )
      nodes$Leaf[i] <- ifelse( length(tokens) == 6, TRUE, FALSE )
    }
  }
  return( nodes ) 
}

#' @title parse_rpart
#' @description Extract splits from an rpart.object returned from a call to
#' rpart().
#' @details This function takes as its input an rpart.object returned from a
#' call to rpart.  It parses this rpart.object using rpart_nodes() and returns
#' the splits in the tree. The data returned include the NodeID of the node to
#' split, the NodeID of that node's parent, the NodeID of that nodes left child
#' and right child, the number of observations in that node, the variable
#' used in the split, the data type for the splitting variable, the logic
#' indicating which observations will go to the node's left child, the
#' value of the splitting variable at which the split ocurrs, the mean response
#' value of the node, and (optionally) the string representation of the node's
#' subgroup. A node's subgroup is defined by the sequence of splits from the
#' root to that node.
#' @seealso \link{rpart_nodes}, \link[rpart]{rpart}, \link[rpart]{rpart.object}
#' @param tree An rpart.object returned from call to rpart().
#' @param include_subgroups A logical value indicating whether or not to include
#' a string representation of the subgroups in the results. Defaults to FALSE.
#' @return A data.frame containing a parsed tree.
#' @examples
#' requireNamespace( "rpart", quietly = TRUE )
#' ## Generate example data containing response, treatment, and covariates
#' N <- 50
#' continuous_response = runif( min = 0, max = 20, n = N )
#' trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6),
#'               replace = TRUE )
#' X1 <- runif( N, min = 0, max = 1 )
#' X2 <- runif( N, min = 0, max = 1 )
#' X3 <- sample( c(0,1), size = N, prob = c(0.2,0.8), replace = TRUE )
#' X4 <- sample( c('A','B','C'), size = N, prob = c(0.6,0.3,0.1), replace = TRUE )
#' 
#' ## Fit an rpart model
#' fit <- rpart::rpart( continuous_response ~ trt + X1 + X2 + X3 + X4,
#'                      control = rpart::rpart.control( maxdepth = 3L ) )
#' fit
#' 
#' ## Parse the results into a new data.frame
#' ex1 <- parse_rpart( fit, include_subgroups = TRUE )
#' ex1
#' @export
parse_rpart <- function( tree,
                         include_subgroups = FALSE ){

  nodes <- rpart_nodes( tree )
  
  dataClasses <- attr(tree$terms, "dataClasses" )
  dataClasses <- dataClasses[ names(dataClasses) != "response" ]
  
  numeric_covars <- names( dataClasses )[ dataClasses == "numeric" ]
  factor_covars <-  names( dataClasses )[ dataClasses == "factor" ]
  
  frame <- tree$frame
  splits_df <- data.frame( nodes$NodeID )
  names( splits_df ) <- "NodeID"
  
  splits_df$NodeID <- as.numeric( splits_df$NodeID )
  splits_df$Parent <- ifelse( splits_df$NodeID > 1, floor( splits_df$NodeID/2 ), NA )
  splits_df$LeftChild <- ifelse( tree$frame$var != '<leaf>', 2*splits_df$NodeID, NA )
  splits_df$RightChild <- ifelse( tree$frame$var != '<leaf>',2*splits_df$NodeID + 1, NA )
  
  splits_df$NodeSize <- as.numeric( tree$frame$n )
  
  splits_df$SplitVariable <- ifelse( frame$var != "<leaf>", as.character(frame$var), "" )
  splits_df$SplitVariable <- as.factor( splits_df$SplitVariable )
      
  splits_df$SplitVariableType <- ifelse( splits_df$SplitVariable %in% factor_covars,"NOMINAL","CONTINUOUS" )
  splits_df$SplitVariableType <- ifelse( splits_df$SplitVariable %in% c( numeric_covars, factor_covars ),
                                        splits_df$SplitVariableType, "" )
  
  splits_df$SplitVariableType <- as.factor( splits_df$SplitVariableType )
  
  
  splits_df$LeftLogic <- NA
  splits_df$SplitValue <- NA
  
  for( i in 1:NROW( splits_df ) ){
    
    if( !is.na( splits_df$LeftChild[i] ) ){
      
      left_split <- nodes$Split[ nodes$NodeID == splits_df$LeftChild[[i]] ]
      
      splits_df$LeftLogic[i] <- ifelse( grepl( "<=", left_split ), "<=",
                                       ifelse( grepl( ">=", left_split ), ">=",
                                              ifelse( grepl( "<", left_split ), "<",
                                                     ifelse( grepl( ">", left_split ), ">","=" ))))
      
      splits_df$SplitValue[i] <- strsplit( left_split, splits_df$LeftLogic[i] )[[1]][[2]]
                                         
    }
  }
  
  splits_df$LeftLogic <- ifelse( is.na( splits_df$LeftLogic), '', splits_df$LeftLogic )
  splits_df$SplitValue <- ifelse( is.na( splits_df$SplitValue ), '', splits_df$SplitValue )
  
  splits_df$MeanResponse <- nodes$yval
  
  splits_df$SplitVariable <- gsub(  pattern = "covariates$", replacement = "", splits_df$SplitVariable, fixed = TRUE )
  
  # Re-order columns
  VARLIST <- c( "NodeID","Parent","LeftChild","RightChild","NodeSize",
                "SplitVariable","LeftLogic","SplitValue","SplitVariableType",
                "MeanResponse" )
  
  splits_df <- splits_df[,VARLIST]
  
  if( include_subgroups ){
    
    subgroups <- get_subgroup_splits( nodes, splits_df )
    splits_df$Subgroup <- factor( subgroups )
    splits_df$Subgroup <- gsub( pattern = "covariates$", replacement = "", splits_df$Subgroup, fixed = TRUE ) 
  }
  return( splits_df )
}

#' @title rpart_wrapper
#' @description A wrapper function to rpart.
#' @details This function provides a wrapper to rpart that provides a convenient
#' interface for specifying the response variable and covariates for the
#' rpart model. The user may indicate whether the tree should be pruned to
#' the size that yields the smallest cross-validation error. An rpart.object
#' is returned.
#' @seealso \link[rpart]{rpart}, \link[rpart]{rpart.object},
#' \link[survival]{Surv}
#' @param response Response variable to use in rpart model.
#' @param response_type Class of response variable.
#' @param covariates Covariates to use in rpart model.
#' @param tree_builder_parameters A named list of parameters to pass to rpart.
#' This includes all input parameters that rpart can take.
#' @param prune Logical variable indicating whether the tree shold be
#' pruned to the subtree with the smallest cross-validation error. Defaults to
#' FALSE.                        
#' @return An object of class rpart.
#' @examples
#' ## Generate example data containing response, treatment, and covariates
#' N <- 100
#' continuous_response = runif( min = 0, max = 20, n = N )
#' trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6), replace = TRUE )
#' X1 <- runif( N, min = 0, max = 1 )
#' X2 <- runif( N, min = 0, max = 1 )
#' X3 <- sample( c(0,1), size = N, prob = c(0.2,0.8), replace = TRUE )
#' X4 <- sample( c('A','B','C'), size = N, prob = c(0.6,0.3,0.1), replace = TRUE )

#' covariates <- data.frame( trt )
#' names( covariates ) <- "trt"
#' covariates$X1 <- X1
#' covariates$X2 <- X2
#' covariates$X3 <- X3
#' covariates$X4 <- X4

#' ## Fit an rpart model
#' ex1 <- rpart_wrapper( response = continuous_response, covariates = covariates )
#' ex1
#' @export
rpart_wrapper <- function( response,
                           response_type = NULL,
                           covariates = NULL,
                           tree_builder_parameters = NULL,
                           prune = FALSE ){
  requireNamespace( "rpart", quietly = TRUE )
  
  if( is.null( tree_builder_parameters ) )
      tree_builder_parameters <- list()
  
  if( is.null( response_type ) ){
      
      if( class( response ) == "Surv")
          response_type <- "survival"
      if( is.numeric( response ) )
          response_type <- "continuous"
      if( is.binary( response_type ) )
          response_type <- "binary"
    }
  
  if( is.null( covariates ) ){
    covariates <- data.frame( rep( 1, length( response ) ) )
    names( covariates ) <- "Intercept"
  }
  
  if( "WEIGHTS__" %in% names( covariates ) ){
    tree_builder_parameters$weights <- covariates$WEIGHTS__
    covariates$WEIGHTS__ <- NULL
  }
  
  # Populate missing parameters as needed
  if( "method" %nin% names( tree_builder_parameters ) ){
      
    if( is.factor( response ) )
        tree_builder_parameters$method <- "class"
      
    if( is.numeric( response ) )
        tree_builder_parameters$method <- "anova"
    
    if( response_type == "survival" )
        tree_builder_parameters$method <- "exp"
  }else{
      
    if( is.null( tree_builder_parameters$method ) &&
       tree_builder_parameters$method %nin% c("anova","poisson","class","exp") )
        stop( "ERROR: method must be one of {anova,poisson,class,exp}" )
    
    if( is.null( tree_builder_parameters$method ) && is.factor( response ) )
        tree_builder_parameters$method <- "class"
    
    if( is.null( tree_builder_parameters$method ) && is.numeric( response ) )
        tree_builder_parameters$method <- "anova"
    
    if( is.null( tree_builder_parameters$method ) && response_type == "survival" )
        tree_builder_parameters$method <- "exp"
  }
  
  if( "subset" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$subset <- NULL
  
  if( "na.action" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$na.action <- rpart::na.rpart
  
  if( "model" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$model <- FALSE
  
  if( "x" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$x <- FALSE
  
  if( "y" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$y <- TRUE
  
  if( "control" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$control <- rpart::rpart.control()
  
  if( "weights" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$weights <- NULL
  
  if( "cost" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$cost <- rep(1,NCOL(covariates))
  
  if( "maxdepth" %in% names( tree_builder_parameters ) )
      tree_builder_parameters$control$maxdepth <- tree_builder_parameters$maxdepth
  
  if( "maxcompete" %in% names( tree_builder_parameters ) )
      tree_builder_parameters$control$maxcompete <- tree_builder_parameters$maxcompete
  
  formulaY <- "response~covariates$"
  formulaX <- paste( names(covariates), sep = "", collapse = "+covariates$" )
  
  formula__ <- as.formula( paste0( formulaY, formulaX ) )
  
  #TODO: what is rpart parms?
  fit <- rpart::rpart( formula__,
                      weights = tree_builder_parameters$weights,
                      subset = tree_builder_parameters$subset,
                      na.action = tree_builder_parameters$na.action,
                      #parms = tree_builder_parameters$parms,
                      method = tree_builder_parameters$method,
                      model = tree_builder_parameters$model,
                      x = tree_builder_parameters$x,
                      y = tree_builder_parameters$y,
                      control = tree_builder_parameters$control,
                      cost = tree_builder_parameters$cost )
  
 # Prune fit according to smallest cross-validation error in cptable
  if( prune )
      fit <- rpart::prune( fit, cp = fit$cptable[ which.min(fit$cptable[,"xerror"]),"CP" ] )
  
  return( fit )
  
}

## END OF FILE

