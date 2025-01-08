#!usr/bin/dev/R
#################################################################################
# FILENAME   : tree_functions.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 01/16/2013
# DESCRIPTION: Functions for calling tree builder algorithms and selecting
#              subgroups.
#################################################################################

get_subgroup_splits <- function( nodes, splits_df ){

  N_splits <- NROW( splits_df )
  
  subgroups <- rep( "", N_splits )

  if( N_splits > 1 ){

    for( i in 2:N_splits ){
    
      new_subgroup <- NULL
    
      node_parents <- sort( parents( splits_df$NodeID[[i]] ) )[-1]
      node_parents <- c( node_parents, splits_df$NodeID[[i]] )
      
      for( p in node_parents ){

        new_split <- nodes$Split[ which(nodes$NodeID == p) ]

        if( !grepl( pattern = ">", new_split ) && !grepl( pattern = "<", new_split )  ){
         new_split <- strsplit( new_split, "=" )[[1]]
         new_split_variable <- new_split[1]
         new_split_values <- gsub( pattern = ",", replacement = "','", new_split[2] )
         new_split_values <- paste( "'", new_split_values, "'", sep = "" )
         new_split <- paste( new_split_variable, "%in%c(", new_split_values,")", sep = "" )
        }
           
        new_subgroup <- c( new_subgroup, new_split )

        new_split <- NULL
        new_split_variable <- NULL
        new_splitvalues <- NULL

      }

      new_subgroup <- collapse_redundant_splits( new_subgroup )
      
      subgroups[[i]] <- paste( new_subgroup, sep = "", collapse = " & " )
    }
    
  }
 
  return( subgroups )
}

parents <- function( nodeid ){

  parent_nodes <- NULL
  
  while( nodeid > 1 ){
    
    parent_nodes <- c( parent_nodes, floor( nodeid/2 ) )

    nodeid <- floor( nodeid/2 )

  }
  return( parent_nodes )
}


superior_subgroups <- function( splits,
                                score,
                                threshold,
                                desirable_response = "increasing",
                                mean_response = NULL,
                                eps = 1E-6 ){

  ## Avoid NOTE in R CMD check about no visible binding for global variable
  superior <- NULL
  
  if( desirable_response %nin% c("increasing","decreasing") )
      stop( "ERROR: desirable_response must be in {decreasing, increasing}" )

  superior_subgroups__ <- splits
  superior_subgroups__[,"superior"] <- FALSE
  
  if( desirable_response == "increasing" ){
    superior_subgroups__[,"superior"] <- score - eps > threshold
  }else{
    superior_subgroups__[,"superior"] <- score + eps < threshold
  }

  ## If the scoring_function is not in {desirable_response_proportion} then also
  ## ensure the mean subgroup response is superior to the mean overall response
  if( !is.null( mean_response ) ){
    if( desirable_response == "increasing" ){
      superior_subgroups__[,"superior"] <- superior_subgroups__[,"superior"] & mean_response > mean_response[1]
    }else{
      superior_subgroups__[,"superior"] <- superior_subgroups__[,"superior"] & mean_response < mean_response[1]
    } 
  }

  superior_subgroups__ <- subset( superior_subgroups__, superior == TRUE )
  superior_subgroups__ [,"superior"] <- NULL

  return( superior_subgroups__ )
}

#' @title subgroup
#' @description Subset a user-provided data.frame according to the subgroup
#' specified by a node in a tree.
#' @details After the splits from an rpart.object are extracted by a call to
#' parse_rpart(), the extracted splits define a subgroup for each node. This
#' subgroup can be used to subset a user-provided data.frame. This function takes
#' as its input a data.frame of splits obtained from a call to parse_rpart(), a
#' NodeID indicating which node specifies the desired subgroup, a data.frame of
#' covariates to subset, and (optionally) the associated response data to subset.
#' If only xdata is specified by the user, the subset of xdata implied by the
#' subgroup will be returned. If xdata and ydata are provided by the user, the
#' subset of ydata will be returned (xdata is still required from the user
#' because the subsetting is computed on the covariate values even when the data
#' returned to the user are from ydata).
#' @seealso \link{parse_rpart}, \link[rpart]{rpart}, \link[rpart]{rpart.object}
#' @param splits A data.frame of splits returned from a call to parse_rpart().
#' @param node The NodeID of the node defining the desired split.
#' @param xdata The data.frame of covariates to subset according to the subgroup
#' definition.
#' @param ydata The associated vector of response values to subset according to
#' the subgroup definition. (optional)
#' @return A data.frame containing the data consistent with the specified
#' subgroup.
#' @examples
#' requireNamespace( "rpart", quietly = TRUE )
#' 
#' ## Generate example data containing response, treatment, and covariates
#' N <- 20
#' continuous_response = runif( min = 0, max = 20, n = N )
#' trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6), replace = TRUE )
#' X1 <- runif( N, min = 0, max = 1 )
#' X2 <- runif( N, min = 0, max = 1 )
#' X3 <- sample( c(0,1), size = N, prob = c(0.2,0.8), replace = TRUE )
#' X4 <- sample( c('A','B','C'), size = N, prob = c(0.6,0.3,0.1), replace = TRUE )
#' 
#' covariates <- data.frame( trt )
#' names( covariates ) <- "trt"
#' covariates$X1 <- X1
#' covariates$X2 <- X2
#' covariates$X3 <- X3
#' covariates$X4 <- X4
#' 
#' ## Fit an rpart model
#' fit <- rpart::rpart( continuous_response ~ trt + X1 + X2 + X3 + X4 )
#'
#' ## Return parsed splits with subgroups
#' splits1 <- parse_rpart( fit, include_subgroups = TRUE )
#' splits1
#' 
#' ## Subset covariate data according to split for NodeID 3
#' ex1 <- subgroup( splits = splits1, node = 3, xdata = covariates )
#' ex1
#' 
#' ## Subset response data according to split for NodeID 3
#' ex2 <- subgroup( splits = splits1, node = 3, xdata = covariates, ydata = continuous_response )
#' ex2
#' @export
subgroup <- function( splits,
                      node,
                      xdata,
                      ydata = xdata ){
  
  if( node %nin% splits$NodeID )
      stop( "ERROR: specified node not found" )
  
  subgroup <- splits$Subgroup[ splits$NodeID == node ]
  
  subgroup_data <- ydata
  
  if( subgroup != "" ){
    
    subgroup <- paste( 'xdata$', subgroup, paste = "" )
    
    subgroup <- gsub( pattern = " ", replacement = "", subgroup )
    
    # Pad logical operators with spaces. If this is not done less than
    # inequalities followed by a negative will be intepreted as assignments.
    # For example, X1 < -5 with no spaces will be interpreted as assigning
    # 5 to X1 as follows:  X1<-5 -- i.e. X1 <- 5.
    subgroup <- gsub( pattern = "<-", replacement = "< -", subgroup )
    
    subgroup <- gsub( pattern = "&", replacement = " & xdata$", subgroup )
    
    subgroup_parse_string <- paste( "subset( ydata, ", subgroup, ")", sep = "" )
    
    subgroup_data <- eval( parse( text = subgroup_parse_string ) )
  }
  return( subgroup_data )
}

tree <- function( response,
                  response_type = NULL,
                  trt = NULL,
                  trt_control = 0,
                  covariates = NULL,
                  tree_builder = "rpart",
                  tree_builder_parameters ){
  
  if( is.null( covariates ) )
      covariates <- rep( 1, length( response ) )
  
  if( tree_builder == "rpart" ){
    
    if( !is.null( trt ) && "WEIGHTS__" %in% names( covariates ) ){
      tree_builder_parameters$weights <- covariates$WEIGHTS__
      covariates$WEIGHTS__ <- NULL
    }
    
    tree__ <- rpart_wrapper( response = response,
                             response_type = response_type,
                             covariates = covariates,
                             tree_builder_parameters = tree_builder_parameters )
    
  }else if( tree_builder == "ctree" ){
    
    
    tree__ <- ctree_wrapper( response = response,
                             covariates = covariates,
                             tree_builder_parameters = tree_builder_parameters )
    
  }else if( tree_builder == "mob" ){
    
    tree__ <- mob_wrapper( response = response,
                          covariates = covariates,
                          tree_builder_parameters = tree_builder_parameters ) 
  }
  
  return( tree__ )
}

parse_tree <- function( tree, tree_builder = NULL ){
  
  if( is.null( tree_builder ) ){
    
    if( is( tree, "rpart" ) ){
      tree_builder <- "rpart"
    }else if( is( tree, "BinaryTree" ) ){
      tree_builder <- "ctree"
    }else if( is( tree, "mob" ) ){
      tree_builder <- "mob"
    }
    
  }
  
  if( tree_builder == "rpart" ){
    parsed__ <- parse_rpart( tree, include_subgroups = TRUE )
  }else if( tree_builder %in% c("ctree","mob") ){
    # Precisions in ctree limited to six decimal places
    parsed__ <- parse_party( tree, include_subgroups = TRUE, digits = 6 )
  }
  
  return( parsed__ )
}

get_competitors <- function( splits,
                             response, response_type,
                             trt = NULL, trt_control = 0,
                             covariates, tree_builder, tree_builder_parameters ){
    
  VARIABLES_TO_EXCLUDE <- splits$SplitVariable[1] #exclude primary split variable as possible competitor
  running_count_of_competitors <- 0 #count number of competitor splits obtained collected so far

  if( "rootcompete" %nin% names( tree_builder_parameters ) )
      tree_builder_parameters$rootcompete <- 0
  
  n_competitors <- min( tree_builder_parameters$rootcompete, NCOL( covariates ) - 1 )
  
  if( tree_builder %in% c("rpart") ){
    tree_builder_parameters$control$maxcompete <- 0 # trees for competitor splits do not need competitors themselves
    tree_builder_parameters$control$maxdepth <- tree_builder_parameters$competedepth
  }else if( tree_builder %in% c("ctree") ){
    tree_builder_parameters$controls@tgctrl@maxdepth <- as.integer( tree_builder_parameters$competedepth )
  }
  
  COMPETITORS <- splits[0,] #empty data.frame with same columns as in splits
  
  while( running_count_of_competitors < n_competitors ){
    
    ## Remove VARIABLES_TO_EXCLUDE from covariates and tree_builder_parameters$cost
    INDICES_TO_KEEP <- which( names( covariates ) %nin% VARIABLES_TO_EXCLUDE )
    covariates <- subset( covariates, select = names( covariates )[INDICES_TO_KEEP] )
    if( !is(covariates, 'data.frame') )
        stop( 'ERROR: covariates is not a data.frame\n' )
    if( NCOL( covariates ) == 0 )
        stop( 'ERROR: covariates has no columns\n' )
    
    
    tree_builder_parameters$cost <- tree_builder_parameters$cost[INDICES_TO_KEEP]
    rm( INDICES_TO_KEEP )

    tree_builder_parameters$maxdepth <- tree_builder_parameters$competedepth
    
    competitor_tree__ <- tree( response = response,
                               response_type = response_type,
                               trt = trt,
                               trt_control = trt_control,
                               covariates = covariates,
                               tree_builder = tree_builder,
                               tree_builder_parameters = tree_builder_parameters )
    
    competitor_splits__ <- parse_tree( competitor_tree__, tree_builder = tree_builder )

    if( NROW( competitor_splits__ ) > 1 ){

      competitor_splits__[,"NodeID"] <- -(competitor_splits__[,"NodeID"] + 100*(running_count_of_competitors+1))
      competitor_splits__[,"Parent"] <- -(competitor_splits__[,"Parent"] + 100*(running_count_of_competitors+1))
      competitor_splits__[,"LeftChild"] <- -(competitor_splits__[,"LeftChild"] + 100*(running_count_of_competitors+1))
      competitor_splits__[,"RightChild"] <- -(competitor_splits__[,"RightChild"] + 100*(running_count_of_competitors+1))
      
      COMPETITORS <- rbind( COMPETITORS, competitor_splits__ )
      
      VARIABLES_TO_EXCLUDE <- c( VARIABLES_TO_EXCLUDE, competitor_splits__$SplitVariable[1] )
      
      running_count_of_competitors <- running_count_of_competitors + 1
    }
    else{
      running_count_of_competitors <- n_competitors
    }
  }
  return( COMPETITORS )
}

## END OF FILE

