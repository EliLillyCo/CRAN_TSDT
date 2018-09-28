#!usr/bin/dev/R
#################################################################################
# FILENAME   : parse_party.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 08/13/2013
# DESCRIPTION: Functions for parsing a tree constucted with ctree() or mob() from
#              the party package to extract subgroups.
#################################################################################

#################################################################################
# Implement a container class for trees created by party::ctree()               #
#################################################################################
requireNamespace( "party", quietly = TRUE )

#' @title CTree
#' @description CTree is a container class for trees created by \link[party]{ctree}.
#' @seealso \link[party]{ctree}, \link[party]{BinaryTree-class}
#' @slot tree An object of class \link[party]{BinaryTree-class} produced by
#' \link[party]{ctree}.
#' @slot data Training data.
#' @slot parameters Control parameters
#' @return An object of class CTree
setClass( "CTree",
         
         representation( tree = "BinaryTree",
                         data = "data.frame",
                         parameters = "list"
                        ),
         package = "TSDT" )

setMethod( f = "initialize", signature = "CTree",
          
           definition = function( .Object, tree, data, parameters ){

            .Object@tree = tree
            .Object@data = data
            .Object@parameters = parameters
           
           return( .Object )
})

setMethod( f = "show", signature = "CTree",

           definition = function( object ){
             
              base::print( object@tree )
              flush.console()
})

#################################################################################
# Implement a container class for trees created by party::mob()                 #
#################################################################################
#' @title MOB
#' @description MOB is a container class for trees created by \link[party]{mob}.
#' @seealso \link[party]{mob}, \link[party]{BinaryTree-class}
#' @slot tree An object of class \link[party]{BinaryTree-class} produced by
#' \link[party]{mob}.
#' @slot data Training data.
#' @slot parameters Control parameters
#' @return An object of class MOB
setClass( "MOB",
         
         representation( tree = "BinaryTree",
                         data = "data.frame",
                         parameters = "list" ),
         package = "TSDT" )

setMethod( f = "initialize", signature = "MOB",
          
           definition = function( .Object, tree, data, parameters ){

            .Object@tree = tree
            .Object@data = data
            .Object@parameters = parameters
           
           return( .Object )
})

setMethod( f = "show", signature = "MOB",

           definition = function( object ){
             
              base::print( object@tree )
              flush.console()
})

#' @title ctree_wrapper
#' @description A wrapper function to \link[party]{ctree}
#' @param response Response variable to use in ctree model.
#' @param covariates Covariates to use in ctree model.
#' @param tree_builder_parameters A named list of parameters to pass to
#' \link[party]{ctree}.
#' @seealso \link[party]{ctree}
#' @return An object of class \linkS4class{CTree}
#' @examples
#' requireNamespace( "party", quietly = TRUE )
#' ## From party::ctree() examples:
#' set.seed(290875)
#' airq <- subset(airquality, !is.na(Ozone))
#' 
#' ## Provide response and covariates to fit ctree
#' ex1 <- ctree_wrapper( response = airq$Ozone,
#'                       covariates = subset( airq, select = -Ozone ) )
#' 
#' ## Pass list of control parameters. Note that ctree takes a parameter called
#' ## 'controls' (with an 's'), rather than 'control' as in rpart.
#' ex2 <- ctree_wrapper( response = airq$Ozone,
#'                       covariates = subset( airq, select = -Ozone ),
#'                       tree_builder_parameters = list( controls =
#'                                              party::ctree_control( maxdepth = 2 ) ) )
#' @export
ctree_wrapper <- function( response,
                          covariates = NULL,
                          tree_builder_parameters = list() ){

  requireNamespace( "party", quietly = TRUE )

  df__ <- as.data.frame( response )
  if( !is.null( names( response ) ) ){
    names( df__ ) <- names( response )[[1]]
  }else{
    names( df__ ) <- "response"
  }
  df__ <- cbind( df__, covariates )
  
  if( "subset" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$subset <- NULL
  
  if( "weights" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$weights <- NULL
  
  if( "controls" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$controls <- party::ctree_control()
  
  if( "xtrafo" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$xtrafo <- ptrafo
  
  if( "ytrafo" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$ytrafo <- ptrafo
  
  if( "scores"  %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$scores <- NULL
  
  if( "maxdepth" %in% names( tree_builder_parameters ) )
      tree_builder_parameters$controls@tgctrl@maxdepth <- as.integer( tree_builder_parameters$maxdepth )
  
  formulaY <- "response"
  formulaX <- paste( names(covariates), sep = "", collapse = "+" )
  formula__ <- as.formula( paste0( formulaY, " ~ ", formulaX ) )
  
  fit <- party::ctree( formula__,
                      data = df__,
                      subset = tree_builder_parameters$subset,
                      weights = tree_builder_parameters$weights,
                      controls = tree_builder_parameters$controls,
                      xtrafo = tree_builder_parameters$xtrafo,
                      ytrafo = tree_builder_parameters$ytrafo,
                      scores = tree_builder_parameters$scores )


  ctree__ <- new( "CTree", tree = fit, data = df__, parameters = tree_builder_parameters )
  
  return( ctree__ )
}

#' @title mob_wrapper
#' @description Wrapper function for \link[party]{mob}.
#' @param response Response variable to use in mob model.
#' @param x Covariates passed to model in mob. mob uses fits the formula
#' y ~ x1 + ... + xk | z1 + ... + zl where the variables before the | are passed
#' to the model and the variables after the | are used for partitioning. x
#' represents the x variables. See mob help page for more information.
#' @param z Covariates used to parition the mob model. mob uses fits the formula
#' y ~ x1 + ... + xk | z1 + ... + zl where the variables before the | are passed
#' to the model and the variables after the | are used for partitioning. z
#' represents the z variables. See mob help page for more information.
#' @param covariates An alias for z.
#' @param tree_builder_parameters A named list of parameters to pass to
#' \link[party]{mob}.
#' @seealso \link[party]{mob}
#' @return An object of class \linkS4class{MOB}
#' @export
mob_wrapper <- function( response,
                         x = NULL,
                         z = NULL,
                         covariates = NULL,
                         tree_builder_parameters = list() ){
    
    requireNamespace( "party", quietly = TRUE )

    if( is.null( z ) && !is.null( covariates ) )
        z <- covariates
    
    df__ <- as.data.frame( response )
    if( !is.null( names( response ) ) ){
      names( df__ ) <- names( response )[[1]]
    }else{
      names( df__ ) <- "response"
    }
    
    if( is.null( x ) ){
      x <- data.frame( rep( 1, NROW( df__ ) ) )
      names( x ) <- "ONE"
    }
   
    df__ <- cbind( df__, x, z )

    if( "weights" %nin% names( tree_builder_parameters ) )
        tree_builder_parameters$weights <- rep(1,NROW(response) )
    
    if( "na.action" %nin% names( tree_builder_parameters ) )
        tree_builder_parameters$na.action <- na.omit
    
    if( "model" %nin% names( tree_builder_parameters ) )
        tree_builder_parameters$model <- modeltools::glinearModel

    if( "control" %nin% names( tree_builder_parameters ) )
        tree_builder_parameters$control <- mob_control()

    formulaY <- "response"

    # If no model fitting terms (i.e. x terms) provided fit model on
    # intercept-only -- i.e. mean response in each node.
    formulaX <- paste( names(x), sep = "", collapse = "+" )
    
    formulaZ <- paste( names(z), sep = "", collapse = "+" )
    
    formula__ <- as.formula( paste0( formulaY, ' ~ ', formulaX, ' | ', formulaZ ) )
    
    fit <- party::mob( formula__,
                       data = df__,
                       weights = tree_builder_parameters$weights,
                       na.action = tree_builder_parameters$na.action,
                       model = tree_builder_parameters$model,
                       control = tree_builder_parameters$control )
    
    mob__ <- new( "MOB", tree = fit, data = df__, parameters = tree_builder_parameters )

    return( mob__ ) 
}


# Write a function to extract parent nodes
party_parents <- function( tree, node_id ){

  if( node_id %nin% tree$NodeID )
      stop( "ERROR: node_id not found in tree" )
  
  PARENTS <- NULL
  
  parent_id <- tree[ tree$NodeID == node_id, "Parent" ]

  while( !is.na( parent_id ) ){
    PARENTS <- c( PARENTS, parent_id )
    parent_id <- tree[ tree$NodeID == parent_id, "Parent" ]
  }

  return( PARENTS )
}

get_party_subgroup <- function( tree, node_id ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  NodeID <- NULL;rm( NodeID )
  
  if( node_id %nin% tree$NodeID )
      stop( "ERROR: node_id not found in tree" )

  next_node <- node_id
  
  sg <- NULL
  
  PARENTS <- party_parents( tree, node_id )

  for( i in 1:length(PARENTS) ){

    split <- subset( tree, NodeID == PARENTS[i] )

    if( next_node == split$LeftChild ){
     logic <- split$LeftLogic
    }
    else if( next_node == split$RightChild ){
      logic <- complement_logic( split$LeftLogic )
    }

    sg <- paste0( sg, split$SplitVariable, logic, split$SplitValue )
    
    if( i < length(PARENTS) )
        sg <- paste0( sg, ' & ' )

    next_node <- PARENTS[i]

    rm( split )
    rm( logic )
  }

  return( sg )
}

#' @title parse_party
#' @description Parse output from ctree() and mob() functions in party package.
#' @details Collects text output from party::ctree() or party::mob(), parses the
#' splits, and populates a data.frame with the relevant data.
#' @seealso \link[party]{ctree}, \link[party]{mob}
#' @param tree An object of class BinaryTree or mob resulting from a call to the
#' ctree() or mob() function.
#' @param data data.frame containing covariates used to create tree.
#' @param include_subgroups A logical value indicating whether or not to include
#' a string representation of the subgroups in the results. Defaults to FALSE.
#' @param digits Number of digits for rounding.
#' @return A data.frame containing a parsed tree.
#' @examples
#' requireNamespace( "party", quietly = TRUE )
#' requireNamespace( "modeltools", quietly = TRUE )
#' ## From party::ctree() examples:
#' set.seed(290875)
#' ## regression
#' airq <- subset(airquality, !is.na(Ozone))
#' airct <- party::ctree(Ozone ~ ., data = airq, 
#'                controls = party::ctree_control(maxsurrogate = 3))
#' 
#' ## Parse the results into a new data.frame
#' ex1 <- parse_party( airct )
#' ex1
#' 
#' ## From party::mob() examples:
#' data("BostonHousing", package = "mlbench")
#' ## and transform variables appropriately (for a linear regression)
#' BostonHousing$lstat <- log(BostonHousing$lstat)
#' BostonHousing$rm <- BostonHousing$rm^2
#' ## as well as partitioning variables (for fluctuation testing)
#' BostonHousing$chas <- factor( BostonHousing$chas, levels = 0:1, 
#'                               labels = c("no", "yes") )
#' BostonHousing$rad <- factor(BostonHousing$rad, ordered = TRUE)
#' 
#' ## partition the linear regression model medv ~ lstat + rm
#' ## with respect to all remaining variables:
#' fmBH <- party::mob( medv ~ lstat + rm | zn + indus + chas + nox + age + 
#'              dis + rad + tax + crim + b + ptratio,
#'              control = party::mob_control(minsplit = 40), data = BostonHousing, 
#'              model = modeltools::linearModel )
#'
#' ## Parse the results into a new data.frame
#' ex2 <- parse_party( fmBH )
#' ex2
#' @export
parse_party <- function( tree,
                         data = NULL,
                         include_subgroups = FALSE,
                         digits = NULL ){

  requireNamespace( "party", quietly = TRUE )
  requireNamespace( "modeltools", quietly = TRUE )

  ## Create NULL placeholders to prevent NOTE in R CMD check
  node <- NULL;rm( node )
  NodeID <- NULL;rm( NodeID )
  
  # If the input tree is of one of the container classes extract components
  if( class( tree ) %in% c("CTree","MOB") ){
    data <- tree@data
    parameters <- tree@parameters
    tree_class <- class( tree )
    tree <- tree@tree
  }
  
  if( class( tree )[[1]] %nin% c( "BinaryTree", "mob" ) )
      stop( "ERROR: tree should be an object of class 'BinaryTree' or 'mob'" )
  
  party.output <- as.data.frame( capture.output( tree ) )
  names( party.output ) <- "capture"
  party.output$capture <- as.character( gsub( pattern = "^\\s*", replacement = "", party.output$capture ) )
  party.output$node <- grepl( pattern = "^[0-9]+)", party.output$capture )
  party.output <- subset( party.output, node == TRUE )
  party.output$NodeID <- NA
  party.output$Split <- ""
  
  party.output$SplitVariable <- NA
  party.output$SplitVariableType <- NA
  party.output$LeftLogic <- NA
  party.output$SplitValue <- NA
  
  party.output$Parent <- NA
  party.output$LeftChild <- NA
  party.output$RightChild <- NA
  party.output$NodeSize <- NA
  party.output$Depth <- NA
  party.output$Depth[1] <- 0

  for( i in 1:NROW(party.output) ){
    
    party.output$NodeID[i] <- as.integer( strsplit( party.output$capture[[i]], split = ")" )[[1]][[1]] )

    # Extract SplitVariable, SplitVariableType, LeftLogic, and SplitValue
    
    right_parenth <- gregexpr( ')', party.output$capture[[i]], fixed = TRUE )[[1]][[1]]
    semi_colon <-  gregexpr( ';', party.output$capture[[i]], fixed = TRUE )[[1]][[1]]
    asterisk <- gregexpr( '*', party.output$capture[[i]], fixed = TRUE )[[1]][[1]]
    double_equal <- gregexpr( '==', party.output$capture[[i]], fixed = TRUE )[[1]][[1]]
    
    if( asterisk == -1 && semi_colon == -1 )
        stop_ind <- .Machine$integer.max
    
    if( semi_colon > -1 )
        stop_ind <- semi_colon - 1
    
    if( asterisk > -1 )
        stop_ind <- asterisk

    # Extract sustring containing split data (variable, logic, and value)
    party.output$Split[[i]] <- substr( party.output$capture[[i]],
                                       start = right_parenth + 1,
                                       stop = stop_ind )
    
    party.output$Split[[i]] <- as.character( gsub( pattern = "^\\s*", replacement = "", party.output$Split[[i]] ) )

    # If split does not contain an asterisk it is not a terminal node and
    # therefore contains split information
    if( asterisk == -1 ){

      # Nominal splits are of the form: SplitVariable == {val1,val2,...,valn}
      if( double_equal > -1 ){
        
        party.output$SplitVariableType[[i]] <- 'NOMINAL'
        
        party.output$SplitVariable[[i]] <- strsplit( party.output$Split[[i]], split = ' == ' )[[1]][[1]]
        party.output$SplitValue[[i]] <- strsplit( party.output$Split[[i]], split = ' == ' )[[1]][[2]]

        # Remove spaces, '{', and '}' from split value
        party.output$SplitValue[[i]] <- gsub( pattern = " ", replacement = "", party.output$SplitValue[[i]], fixed = TRUE )
        party.output$SplitValue[[i]] <- gsub( pattern = "}", replacement = "", party.output$SplitValue[[i]], fixed = TRUE )
        party.output$SplitValue[[i]] <- gsub( pattern = "{", replacement = "", party.output$SplitValue[[i]], fixed = TRUE )
        party.output$SplitValue[[i]] <- gsub( pattern = "}", replacement = "", party.output$SplitValue[[i]], fixed = TRUE )

        # Enclose individual elements in quotes
        party.output$SplitValue[[i]] <- gsub( pattern = ",", replacement = "','", party.output$SplitValue[[i]], fixed = TRUE )


        # Enclose split value in R collection operator: c(...)
        party.output$SplitValue[[i]] <- paste0( "c('", party.output$SplitValue[[i]], "')" )
        
        party.output$LeftLogic[[i]] <- "%in%"
        
        
      }else{ #Continuous (or ordinal) splits are of the form: SplitVariable < SplitValue (or <=, >, >=).
        
        party.output$SplitVariableType[[i]] <- 'CONTINUOUS'
        
        party.output$SplitVariable[[i]] <- strsplit( party.output$Split[[i]], split = ' ' )[[1]][[1]]
        party.output$LeftLogic[[i]] <- strsplit( party.output$Split[[i]], split = ' ' )[[1]][[2]]
        party.output$SplitValue[[i]] <- strsplit( party.output$Split[[i]], split = ' ' )[[1]][[3]]
      }
    }
  }
  
  # Populate Parent with node ID of parent
  for( i in 1:NROW( party.output ) ){
    
    if( party.output$NodeID[[i]] > 1 ){
      
      if( party.output$NodeID[[i]] > party.output$NodeID[[i-1]] )
          
          party.output$Parent[[i]] <- party.output$NodeID[[i-1]]
    }
  }
  
  # Populate LeftChild and RightChild
  for( i in 1:NROW( party.output ) ){  
    
    ChildNodes <- ( party.output$NodeID[ party.output$Parent %in% party.output$NodeID[[i]] ] )
    
    if( length( ChildNodes ) == 2 ){
      party.output$LeftChild[[i]] <- min( ChildNodes )
      party.output$RightChild[[i]] <- max( ChildNodes ) 
    }
    rm( ChildNodes )
  }  
  
  # Remove duplicate rows
  party.output <- party.output[!duplicated( party.output$NodeID),]

  # Populate depth
  for( i in 1:NROW( party.output ) ){
    
    if( i > 1 ){
      parent_id <- party.output$Parent[[i]]
      parent <- subset( party.output, NodeID == parent_id )
      party.output$Depth[[i]] <- parent$Depth + 1
      rm( parent_id )
      rm( parent )
    }
  }

  # Replace <NA> text with empty string
  party.output$SplitVariable <- na2empty( as.character( party.output$SplitVariable ) )
  party.output$LeftLogic <- na2empty( as.character( party.output$LeftLogic ) )
  party.output$SplitValue <- na2empty( as.character( party.output$SplitValue ) )
  party.output$SplitVariableType <- na2empty( as.character( party.output$SplitVariableType ) )

  row.names( party.output ) <- NULL

  if( exists( "parameters" ) ){
    if( "maxdepth" %in% names( parameters ) ){
      party.output <- subset( party.output, party.output$Depth <= parameters$maxdepth )
    }
  }
  
  # Populate Subgroup, which is needed to compute NodeSize.
  # If include_subgroups == FALSE then remove Subgroup from output
  party.output$Subgroup <- ""
  
  if( NROW( party.output ) > 1 ){
    
    for( i in 2:NROW( party.output ) ){
        
      party.output$Subgroup[[i]] <- get_party_subgroup( party.output, party.output$NodeID[[i]] )
      
    }
  }
  
  # For each node get number of records
  if( !is.null( data ) ){
    
    party.output$NodeSize <- NA
    party.output$MeanResponse <- NA
    
    party.output$NodeSize[[1]] <- NROW( data )
    party.output$MeanResponse[[1]] <- mean( data$response, na.rm = TRUE )
    
    if( NROW( party.output ) > 1 ){
      
      for( i in 2:NROW( party.output ) ){
        
        data__ <- get_subset( data, party.output$Subgroup[[i]], digits = digits )
        
        party.output$NodeSize[[i]] <- NROW( data__ )
        party.output$MeanResponse[[i]] <- mean( data__$response, na.rm = TRUE )
        rm( data__ )
      }
    }
    
  }
  
  # Re-order columns

  COLUMNS <- c( "NodeID", "Parent", "LeftChild", "RightChild" )
  
  if( !is.null( data ) )
      COLUMNS <- c( COLUMNS, "NodeSize" )
  
  COLUMNS <- c(COLUMNS,"Depth","SplitVariable","LeftLogic","SplitValue",
               "SplitVariableType" )
  
  if( !is.null( data ) )
      COLUMNS <- c( COLUMNS, "MeanResponse" )
  
  if( include_subgroups )
      COLUMNS <- c( COLUMNS, "Subgroup" )
  
  party.output <- party.output[,COLUMNS]
  
  return( party.output )
}

## END OF FILE
