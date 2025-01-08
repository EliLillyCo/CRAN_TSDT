#!usr/bin/dev/R
#################################################################################
# FILENAME   : TSDT_helper_functions.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 08/28/2013
# DESCRIPTION: Collection of various helper functions to perform tasks needed
#              for TSDT subgroup identification tool.
#################################################################################

space2underscore <- function( x ){

  if( class( x ) %in% c( "data.frame", "matrix" ) ){
    for( j in 1:NCOL( x ) ){
      if( any( class( x[,j] ) %in% c( "character","factor" ) ) )
          x[,j] <- gsub( pattern = "\\s", replacement = "_", perl = TRUE, x[,j] )
    }
  }else if( any( class( x ) %in% c( "character","factor" ) ) ){
    
    x <- gsub( pattern = "\\s", replacement = "_", perl = TRUE, x )    
  }
  return( x )
}

get_generic_subgroup <- function( subgroup, anon = "xxxxx" ){
  
  common_regex <- paste0( "[^=", anon, "].*[0-9]*" )

  subgroup <- strsplit( subgroup, " & " )[[1]]
  
  if( all( subgroup %nin% 'Overall' ) )
    subgroup <- collapse_redundant_splits( subgroup )
  
  generic_subgroup <- NA
  
  if( length( subgroup ) > 0 ){
    
    for( i in 1:length(subgroup) ){
      
      subgroup[i] <- gsub( pattern = paste0( ">=", common_regex ),
                          replacement = paste0(">=", anon ), subgroup[i] )

      subgroup[i] <- gsub( pattern = paste0( "<=", common_regex ),
                          replacement = paste0("<=", anon ), subgroup[i] )
      
      subgroup[i] <- gsub( pattern = paste0( "<", common_regex ),
                          replacement = paste0("<", anon ), subgroup[i] )
      
      subgroup[i] <- gsub( pattern = paste0( ">", common_regex ),
                          replacement = paste0(">", anon ), subgroup[i] )
      
    }
    
    generic_subgroup <- paste( subgroup, sep = "", collapse = " & " )
    
  }
  
  return( generic_subgroup )
}


#' @title get_suggested_subgroup
#' @description Get a string definition of the suggested subgroup definition.
#' @details Subgroups are reported in an anonymized fashion -- e.g. a subgroup
#' defined on a variable X1 could be reported as X1<xxxxx, 'xxxxx' is a string
#' used to represent an exact numeric cutoff. For each anonymized subgroup,
#' the distribution of exact numeric cutpoints is retained across all
#' bootrstrapped samples. TSDT then provides a suggested cutoff got each
#' anonymized subgroup. By default, this suggested cutoff is the median of the
#' observed cutpoints. Note that this anonymization applies only to numeric
#' splitting variables. Categorical splitting variables are not anonymized.
#' @param anonymized_subgroup A string containing the the anonymized subgroup.
#' @param suggested_cutoff A string containing the suggested cutoff.
#' @param anon The anonymization string. By default this is 'xxxxx'.
#' @export
#' @examples
#' set.seed(0)
#' N <- 200
#' continuous_response = runif( min = 0, max = 20, n = N )
#' trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6), replace = TRUE )
#' X1 <- runif( N, min = 0, max = 1 )
#' X2 <- runif( N, min = 0, max = 1 )
#' X3 <- sample( c(0,1), size = N, prob = c(0.2,0.8), replace = TRUE )
#' X4 <- sample( c('A','B','C'), size = N, prob = c(0.6,0.3,0.1), replace = TRUE )
#' covariates <- data.frame( X1 )
#' covariates$X2 <- X2
#' covariates$X3 <- factor( X3 )
#' covariates$X4 <- factor( X4 )
#' 
#' ## Create a TSDT object
#' ex1 <- TSDT( response = continuous_response,
#'              trt = trt, trt_control = 'Control',
#'              covariates = covariates[,1:4],
#'              inbag_score_margin = 0,
#'              desirable_response = "increasing",
#'              oob_score_margin = 0,
#'              min_subgroup_n_control = 10,
#'              min_subgroup_n_trt = 20,
#'              maxdepth = 2,
#'              rootcompete = 2 )
#' 
#' ## Show summary statistics
#' summary( ex1 )
#'
#' ## Get the anonymized subgroup defined on X1
#' anonymized_subgroup <- as.character( ex1@superior_subgroups$Subgroup[2] )
#'
#' ## Get the suggested cutoff for this subgroup
#' suggested_cutoff <- as.character( ex1@superior_subgroups$Suggested_Cutoff[2] )
#'
#' ## Get the suggested subgroup
#' get_suggested_subgroup( anonymized_subgroup = anonymized_subgroup,
#'                         suggested_cutoff = suggested_cutoff )
#' @export
get_suggested_subgroup <- function( anonymized_subgroup, suggested_cutoff, anon = "xxxxx" ){

  anonymized_subgroup <- as.character( anonymized_subgroup )
  suggested_cutoff <- as.character( suggested_cutoff )

  # Remove anon substring from anonymized subgroups
  anonymized_subgroup <- gsub( pattern = anon, replacement = "", anonymized_subgroup )

  # Split compound subgroup definitions into individual splits
  anon_sg <- strsplit( anonymized_subgroup, split = ' & ', fixed = TRUE )[[1]]
  
  cutoffs <- gsub( pattern = '[', replacement = "", fixed = TRUE, suggested_cutoff )
  cutoffs <- gsub( pattern = ']', replacement = "", fixed = TRUE, cutoffs )

  # Split suggested cutoffs into individual cutoffs
  cutoffs <- strsplit( cutoffs, split = ";" )[[1]]
  
  if( length( anon_sg ) != length( cutoffs ) )
      stop( "ERROR: number of anonymized subgroups must equal the number of suggested cutoffs" )

  # Append suggested cutoff to each subgroup (numeric cutoffs only)
  for( i in 1:length(anon_sg) )
      if( !grepl( pattern = "%in%", anon_sg[[i]], fixed = TRUE ) )
          anon_sg[[i]] <- paste0( anon_sg[[i]], cutoffs[[i]] )

  anon_sg <- paste( anon_sg, sep = "", collapse = " & " )
  
  return( anon_sg )
}





# Get proportion of TRUE values by group
prop_true <- function( x, by ){

  result <- data.frame( Subgroup = character(0),
                        Frequency = numeric(0) )

  GROUPS <- unique( by )

  for( group in GROUPS ){
    
    freq <- ifelse( any( x[by==group] ), prop.table(table( x[by==group] ))[["TRUE"]], 0 )
    result <- rbind( result, cbind( group, freq ) )
    
  }
  
  names( result ) <- c( "Subgroup", "Frequency" )
  
  return( result )
  
}


get_experimental_trt_arm <- function( trt, trt_control ){
  # Get non-control trt
  trt_experimental <- NULL
  if( !is.null( trt ) ){
    
    TRT_VALUES <- unique( trt )

    if( trt_control %nin% TRT_VALUES )
      stop( "ERROR: trt_control value not found in trt" )
    
    if( is.factor( TRT_VALUES ) )
      TRT_VALUES <- as.character( TRT_VALUES )
    
    trt_experimental <- TRT_VALUES[ TRT_VALUES != trt_control ]
  }
  return( trt_experimental )
}

# Get scoring function and associated parameters (if any)
get_scoring_function <- function( scoring_function = NULL,
                                  data,
                                  scoring_function_parameters = NULL ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  response_type <- NULL;rm( response_type )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  
  if( !is.null( scoring_function ) && !exists( "scoring_function_name" ) )
    scoring_function_name <- deparse( substitute( scoring_function ) )
  
  # Select scoring_function if none provided
  if( is.null( scoring_function ) ){
    
    # Default scoring function for binary or continuous response, single-arm is mean_response
    if( response_type %in% c("binary","continuous") && !exists( "trt_var" ) ){
      scoring_function <- mean_response
      scoring_function_name <- "mean_response"
    }
    
    # Default scoring function for binary or continuous response, two-arm is treatment_effect
    else if( response_type %in% c("binary","continuous") && exists( "trt_var" ) ){
      scoring_function <- treatment_effect
      scoring_function_name <- "treatment_effect"
    }

    else if( response_type == "survival" ){

      if( !exists( "trt_var" ) ){
        scoring_function <- survival_time_quantile
        scoring_function_name <- "survival_time_quantile"
      }
        
      else if( exists( "trt_var" ) ){
        scoring_function <- diff_survival_time_quantile
        scoring_function_name <- "diff_survival_time_quantile"
      }
    }
  }
  
  return( list( scoring_function = scoring_function,
                scoring_function_name = scoring_function_name,
                scoring_function_parameters = scoring_function_parameters ) )
}

get_samples <- function( data, trt, trt_control, n_samples, sampling_method, inbag_proportion ){

  if( sampling_method %nin% c("bootstrap","subsample") )
    stop( "ERROR: sampling_method must be one of {bootstrap,subsample}" )
    
  if( sampling_method == "bootstrap" )
    return( bootstrap( data, trt = trt, trt_control = trt_control, n_samples = n_samples ) )

  else{

    samples <- subsample( data, trt = trt, trt_control = trt_control, n_samples = n_samples,
                          training_fraction = inbag_proportion,
                          validation_fraction = 0,
                          test_fraction = 1 - inbag_proportion )

    # NOTE: the container class Bootstrap is used regardless of whether
    # subsample or bootstrap is used for subsampling.
    samples_to_bags <- lapply( rep( "Bootstrap", n_samples ), new )

    # Map Subsample slots to Bootstrap slots
    for( i in 1:n_samples ){

      samples_to_bags[[i]]@inbag <- samples[[i]]@training
      samples_to_bags[[i]]@oob <- samples[[i]]@test
    }
    
    return( samples_to_bags )
  }
}


get_superior_subgroups <- function( splits,
                                    desirable_response,
                                    inbag_score_margin,
                                    oob_score_margin,
                                    eps,
                                    scoring_function_parameters = NULL ){
  

  ## Avoid NOTE in R CMD check about no visible binding for global variable
  NodeID <- Subgroup <- response_type <- NULL
  
  unpack_args( scoring_function_parameters )

  splits <- subset( splits, NodeID==1 | nchar(trimws(Subgroup)) > 0 )
  
  if( response_type != "survival" ){
    inbag_mean_response <- splits$Mean_Inbag_Response  
    oob_mean_response <- splits$Mean_OOB_Response
  }
  else{
    inbag_mean_response <- splits$Inbag_Score
    oob_mean_response <- splits$OOB_Score
  }

  inbag_mean_response__ <- subset( inbag_mean_response, splits[,"NodeID"]==1 | nchar(trimws(splits[,"Subgroup"])) > 0 )
  oob_mean_response__ <- subset( oob_mean_response, splits[,"NodeID"]==1 | nchar(trimws(splits[,"Subgroup"])) > 0 )
  splits__ <- subset( splits, splits[,"NodeID"]==1 | nchar(trimws(splits[,"Subgroup"])) > 0 )

  # Identify superior in-bag subgroups                              
  inbag_superior_subgroup <- superior_subgroups( splits = splits__,
                                                 mean_response = inbag_mean_response__,
                                                 score = splits__$Inbag_Score,
                                                 threshold = splits__$Inbag_Score[[1]] + abs(splits__$Inbag_Score[[1]]) * inbag_score_margin,
                                                 desirable_response = desirable_response,
                                                 eps = eps )

  splits$Superior_Inbag_Subgroup <- ifelse( splits$NodeID %in% inbag_superior_subgroup$NodeID, TRUE, FALSE )
        
  if( desirable_response == "decreasing" && inbag_score_margin > 0 )
      cat( "\nNOTE: If desirable_response == decreasing then inbag_score_margin should be negative\n\n" )
  
  else if( desirable_response == "increasing" && inbag_score_margin < 0 )
      cat( "\nNOTE: If desirable_response == increasing then inbag_score_margin should be positive\n\n" )
  
  # Identify superior out-of-bag subgroups                              
  oob_superior_subgroup <- superior_subgroups( splits = splits__,
                                               mean_response = oob_mean_response__,
                                               score = splits__$OOB_Score,
                                               threshold = splits__$OOB_Score[[1]] + abs( splits__$OOB_Score[[1]] ) * oob_score_margin,
                                               desirable_response = desirable_response,
                                               eps = eps )
  
  splits$Superior_OOB_Subgroup <- ifelse( splits$NodeID %in% oob_superior_subgroup$NodeID, TRUE, FALSE )

  # Add root node to superior subgroups to get Overall mean response, treatment effect, etc.
  splits$Superior_Inbag_Subgroup[ which( splits$NodeID == 1 ) ] <- TRUE
  splits$Superior_OOB_Subgroup[ which( splits$NodeID == 1 ) ] <- TRUE
  splits$Subgroup[ which( splits$NodeID == 1 ) ] <- "Overall"
  
  # Check external consistency of identified subgroups -- i.e. if a subgroup
  # is identified as superior in the in-bag data is it also superior in the
  # out-of-bag data?
      
  splits$External_Consistency <- ifelse( splits$Superior_Inbag_Subgroup,
                                         splits$Superior_Inbag_Subgroup == splits$Superior_OOB_Subgroup,
                                         NA )
  
  return( splits )
}

get_mean_response <- function( splits, data, trt, trt_control ){

  Mean_Response <- rep( NA, NROW( splits ) )  
    
  for( j in 1:NROW( splits ) ){
    
    node <- splits$NodeID[j]
    
    data_subgroup <- subgroup( splits, node = node, xdata = data )
    
    if( is.null( trt ) ){
      if( "WEIGHTS__" %nin% names( data ) )
          Mean_Response[j] <- mean( data_subgroup$response, na.rm = TRUE )
      else
          Mean_Response[j] <- weighted.mean( data_subgroup$response, w = data_subgroup$WEIGHTS__, na.rm = TRUE )
    } # END is.null( trt )
    else{
      if( "WEIGHTS__" %nin% names( data ) )
          Mean_Response[j] <- mean( data_subgroup$response[ data_subgroup$trt != trt_control ], na.rm = TRUE )
      else
          Mean_Response[j] <- weighted.mean( data_subgroup$response[ data_subgroup$trt != trt_control ],
                                             w = data_subgroup$WEIGHTS__[ data_subgroup$trt != trt_control ],
                                             na.rm = TRUE)
    }# END else
  }# END for loop
  
  return( Mean_Response )
}


get_score <- function( splits,
                       inbag,
                       oob,
                       scoring_function,
                       scoring_function_parameters = NULL,
                       min_subgroup_n_control,
                       min_subgroup_n_trt,
                       min_subgroup_n_oob_control,
                       min_subgroup_n_oob_trt ){


  ## Create NULL placeholders to prevent NOTE in R CMD check
  scoring_function_name <- NULL;rm( scoring_function_name )
  trt_control <- NULL;rm( trt_control )
  trt_var <- NULL;rm( trt_var )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )

  
  Inbag_Score <- rep( NA, NROW( splits ) )
  OOB_Score <- rep( NA, NROW( splits ) )
  
  for( j in 1:NROW( splits ) ){
      
    node <- splits$NodeID[j]
    
    inbag_subgroup <- subgroup( splits, node = node, xdata = inbag )
    
    oob_subgroup   <- subgroup( splits, node = node, xdata = oob )
    
    if( !exists( "trt_var" ) ){  # single-arm
        # inbag_arglist is a list of arguments to pass to the default scoring
        # functions. It is a list nested within a list. The inner list contains
        # all the parameters the scoring function needs to compute its score
        # value. The outer list wraps the inner list and names the inner list
        # arglist. This is required because the do.call function passes the
        # inner list to the scoring function and the scoring function takes
        # a list argument named arglist. Similarly for oob_arglist.
      if( NROW( inbag_subgroup ) >= min_subgroup_n_trt &&  NROW( oob_subgroup ) >= min_subgroup_n_oob_trt ){
        
        inbag_arglist <- list( scoring_function_name = scoring_function_name,
                              data = inbag_subgroup,
                              scoring_function_parameters = scoring_function_parameters )
        
        Inbag_Score[[j]] <- do.call( scoring_function_wrapper, args = inbag_arglist )
        
        oob_arglist <- list( scoring_function_name = scoring_function_name,
                            data = oob_subgroup,
                            scoring_function_parameters = scoring_function_parameters )
        
        OOB_Score[[j]] <- do.call( scoring_function_wrapper, args = oob_arglist )
      }
    }
    
    else{ # two-arm

      inbag_control_n <- NROW( subset( inbag_subgroup, inbag_subgroup[,trt_var] == trt_control ) )
      inbag_trt_n     <- NROW( subset( inbag_subgroup, inbag_subgroup[,trt_var] != trt_control ) )
      
      oob_control_n <- NROW( subset( oob_subgroup, oob_subgroup[,trt_var] == trt_control ) )
      oob_trt_n     <- NROW( subset( oob_subgroup, oob_subgroup[,trt_var] != trt_control ) )
      
      if( inbag_control_n >= min_subgroup_n_control && inbag_trt_n >= min_subgroup_n_trt &&
          oob_control_n >= min_subgroup_n_oob_control && oob_trt_n >= min_subgroup_n_oob_trt ){
      
        inbag_arglist <- list( scoring_function_name = scoring_function_name,
                               data = inbag_subgroup,
                               scoring_function_parameters = scoring_function_parameters )
      
        Inbag_Score[[j]] <- do.call( scoring_function_wrapper, args = inbag_arglist )
        
        
        oob_arglist <- list( scoring_function_name = scoring_function_name,
                             data = oob_subgroup,
                             scoring_function_parameters = scoring_function_parameters )
        
        OOB_Score[[j]] <- do.call( scoring_function_wrapper, args = oob_arglist )
        
      }
    }
  }
  
  splits$Inbag_Score <- Inbag_Score
  splits$OOB_Score <- OOB_Score
  
  return( splits )
}

populate_tsdt_samples <- function( samples,
                                   response_type,
                                   trt,
                                   trt_control,
                                   min_subgroup_n_control,
                                   min_subgroup_n_trt,
                                   min_subgroup_n_oob_control,
                                   min_subgroup_n_oob_trt,
                                   scoring_function,
                                   scoring_function_parameters = NULL,
                                   desirable_response,
                                   eps,
                                   inbag_score_margin,
                                   oob_score_margin,
                                   tree_builder,
                                   tree_builder_parameters ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  y_var <- NULL;rm( y_var )
  covariate_vars <- NULL;rm( covariate_vars )
  trt_var <- NULL;rm( trt_var )
  INVALID_TREE__ <- NULL;rm( INVALID_TREE__ )
  Superior_Inbag_Subgroup <- NULL;rm( Superior_Inbag_Subgroup )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )

  if( !is.null( tree_builder_parameters ) )
      unpack_args( tree_builder_parameters )
  
  if( exists( "trt_var" ) && is.null( trt_var ) )
      rm( trt_var )
  
  OverallInternalConsistency <- NULL
  OverallExternalConsistency <- data.frame( Subgroup = character(0),
                                            OOB_Superior_Subgroup = logical(0) )

  TSDT_SAMPLES <- lapply( rep( "TSDT_Sample", length(samples) ), new )
  
  CutpointDistribution <- new( "TSDT_CutpointDistribution", Cutpoints = NULL )
  
  for( i in 1:length(samples) ){
       
    inbag_response <- inbag_covariates <- inbag_trt <- NULL
    oob_response <- oob_covariates <- oob_trt <- NULL
    
    # For each TSDT_Sample object populate the inbag and oob slots
    
    inbag_response <- samples[[i]]@inbag[,c(y_var)]
    oob_response   <- samples[[i]]@oob[,c(y_var)]
    
    inbag_covariates <- subset( samples[[i]]@inbag, select = covariate_vars ) 
    oob_covariates   <- subset( samples[[i]]@oob, select = covariate_vars )
    
    ## If two-arm TSDT train model only on experimental arm.
    if( exists( 'trt_var' ) ){
      inbag_trt <- samples[[i]]@inbag[,c(trt_var)]
      oob_trt   <- samples[[i]]@oob[,c(trt_var)]
      
      inbag_response <- subset( inbag_response, inbag_trt != trt_control )
      oob_response   <- subset( oob_response,   oob_trt != trt_control )
      
      inbag_covariates <- subset( inbag_covariates, inbag_trt != trt_control )
      oob_covariates   <- subset( oob_covariates, oob_trt != trt_control )
      
    }else{
      inbag_trt <- NULL
      oob_trt   <- NULL
    }
    
    if( NROW( inbag_covariates ) == 0 )
        stop( "ERROR: no in-bag covariate data" )

    # Subsetting destroys the Surv object, splitting it into a matrix with
    # two columns: time and status. Re-construct the Surv object before
    # proceeding.
    if( response_type == "survival" ){
      inbag_response <- Surv( inbag_response[,"time"], inbag_response[,"status"] )
      oob_response <- Surv( oob_response[,"time"], oob_response[,"status"] )
    }
    
    # Fit tree
    tree__ <- tree( response = inbag_response,
                    response_type = response_type,
                    trt = inbag_trt,
                    trt_control = trt_control,
                    covariates = inbag_covariates,
                    tree_builder = tree_builder,
                    tree_builder_parameters = tree_builder_parameters )

    # Extract data.frame of splits (i.e. subgroups) from tree
    splits__ <- tryCatch( parse_tree( tree__, tree_builder = tree_builder ),
                          warning = function(w){print(w);flush.console()},
                          error = function(e){ print(e)
                                               flush.console()
                                               INVALID_TREE__ <<- tree__
                                               stop("NOTE: invalid tree accessible in variable INVALID_TREE__\n")})
    
    # Add competitor root node competitor splits
    if( tree_builder_parameters$rootcompete > 0 && NROW( splits__ ) > 1 ){
      
      competitors__ <- get_competitors( splits = splits__,
                                        response = inbag_response,
                                        response_type = response_type,
                                        trt = inbag_trt,
                                        trt_control = trt_control,
                                        covariates = inbag_covariates,
                                        tree_builder = tree_builder,
                                        tree_builder_parameters = tree_builder_parameters )
      
      # Give competitor splits negative NodeIDs
      if( NROW(competitors__) > 0 ){
        splits__ <- rbind( splits__, competitors__ )
      }
      
    }
    
    # Rename MeanResponse to Mean_Inbag_Response
    if( "MeanResponse" %in% names( splits__ ) ){
      splits__$Mean_Inbag_Response <- splits__$MeanResponse
      splits__$MeanResponse <- NULL
    }else{
      
      splits__$Mean_Inbag_Response <- get_mean_response( splits = splits__,
                                                         data = samples[[i]]@inbag,
                                                         trt = inbag_trt,
                                                         trt_control = trt_control )
    }
    
    # The mean response values are meaningless for a survival outcome
    # so do not include these in the output if the outcome is survival
    if( response_type == "survival" ){
      splits__$Mean_Inbag_Response <- NULL
    }else{
      # Compute Mean_OOB_Response if not a survival outcome
      splits__$Mean_OOB_Response <- get_mean_response( splits = splits__,
                                                       data = samples[[i]]@oob,
                                                       trt = oob_trt,
                                                       trt_control = trt_control )
    }

    # Use scoring function to compute in-bag and oob scores for each subgroup
    splits__ <- get_score( splits = splits__,
                           inbag = samples[[i]]@inbag,
                           oob = samples[[i]]@oob,
                           scoring_function = scoring_function,
                           scoring_function_parameters = scoring_function_parameters,
                           min_subgroup_n_control = min_subgroup_n_control,
                           min_subgroup_n_trt = min_subgroup_n_trt,
                           min_subgroup_n_oob_control = min_subgroup_n_oob_control,
                           min_subgroup_n_oob_trt = min_subgroup_n_oob_trt )

    
    # Identify superior in-bag and oob subgroups
    splits__ <- get_superior_subgroups( splits = splits__,
                                        desirable_response = desirable_response,
                                        inbag_score_margin = inbag_score_margin,
                                        oob_score_margin = oob_score_margin,
                                        eps = eps,
                                        scoring_function_parameters = scoring_function_parameters )
    
    # Collect superior subgroups for computation of internal and external
    # consistency.
    consistency_subgroups <- subset( splits__, Superior_Inbag_Subgroup == TRUE )

    # Anonymize subgroups
    if( NROW(consistency_subgroups) > 0 ){
        for( k in 1:NROW(consistency_subgroups) ){
            
            if( consistency_subgroups$Subgroup[[k]] %nin% c('Overall') )
                set_cutpoints( CutpointDistribution, consistency_subgroups$Subgroup[[k]] )
            
            consistency_subgroups$Subgroup[k] <- get_generic_subgroup( subgroup = consistency_subgroups$Subgroup[[k]] )
        }
    }
    
    OverallExternalConsistency <- rbind( OverallExternalConsistency, consistency_subgroups ) 
      
      # Check internal consistency of identified subgroups -- i.e. what are the
      # most commonly identified subgroups across in-bag bootstrapped samples?
      # For the purpose of the comparison we are not concerned with the individual
      # cutpoints for numeric splitting variables. That is, a subgroup such as
      # X1 < 5.4453 would be considered the same as X1 < 5.7656 and we construct a
      # generic cutpoint of the form X1 < xxxxx. Thus, we would consider this
      # subgroup to have appeared twice here.
      
      OverallInternalConsistency <- c( OverallInternalConsistency,
                                       unique( consistency_subgroups$Subgroup ) )
      
      rm( consistency_subgroups )
      row.names( splits__ ) <- NULL
      
      next_TSDT_sample <- new( "TSDT_Sample",
                               inbag = samples[[i]]@inbag,
                               oob   = samples[[i]]@oob,
                               subgroups  = splits__ )
      
      TSDT_SAMPLES[[i]] <- next_TSDT_sample
      
      rm( next_TSDT_sample )
    }
  
  return( list( TSDT_SAMPLES = TSDT_SAMPLES,
                OverallExternalConsistency = OverallExternalConsistency,
                OverallInternalConsistency = OverallInternalConsistency,
                CutpointDistribution = CutpointDistribution ) )
}

compute_consistency <- function( OverallExternalConsistency,
                                 OverallInternalConsistency,
                                 n_samples ){
          
  if( length( OverallInternalConsistency ) == 0 ){
    
      consistency <- data.frame( Subgroup = "None",
                                 Internal_Consistency = NA,
                                 External_Consistency = NA )
  }
  else{
    consistency <- data.frame( table( OverallInternalConsistency ) )
    names( consistency ) <- c("Subgroup","Internal_Consistency")
    consistency$Internal_Consistency <- consistency$Internal_Consistency/n_samples
    
    rownames( consistency ) <- NULL
    
    OverallExternalConsistency <- prop_true( OverallExternalConsistency$Superior_OOB_Subgroup, by = OverallExternalConsistency$Subgroup )
    
    names( OverallExternalConsistency ) <- c( "Subgroup", "External_Consistency" )
    
    consistency <- merge( x = consistency, y = OverallExternalConsistency,
                          merge.x = consistency$Subgroup,
                          merge.y = OverallExternalConsistency$Subgroup,
                          all.x = TRUE, all.y = TRUE )
    
    consistency <- consistency[ order( consistency$Internal_Consistency,
                                       consistency$External_Consistency, decreasing = TRUE ),]
    
    consistency$Internal_Consistency <- as.numeric( as.character( consistency$Internal_Consistency ) )
    consistency$External_Consistency <- as.numeric( as.character( consistency$External_Consistency ) )
    
    if( any( consistency$Internal_Consistency < 0 ) || any( consistency$Internal_Consistency > 1 ) )
      stop( "ERROR: Internal_Consistency must be in [0,1]" )

    if( any( consistency$External_Consistency < 0 ) || any( consistency$External_Consistency > 1 ) )
      stop( "ERROR: External_Consistency must be in [0,1]" )
  }

  return( consistency )
}
  
get_null_scores <- function( n_permutations,
                             response,
                             trt,
                             covariates,
                             tree_builder = tree_builder,
                             tree_builder_parameters = tree_builder_parameters,
                             permute_method,
                             permute_arm,
                             n_samples,
                             sampling_method,
                             inbag_proportion,
                             desirable_response,
                             scoring_function,
                             scoring_function_parameters,
                             inbag_score_margin,
                             oob_score_margin,
                             eps,
                             min_subgroup_n,
                             min_subgroup_n_control,
                             min_subgroup_n_trt,
                             min_subgroup_n_oob_control,
                             min_subgroup_n_oob_trt,
                             n_cpu,
                            trace ){
  
  ## Create NULL placeholders to prevent NOTE in R CMD check
  trt_control <- NULL;rm( trt_control )
  OverallExternalConsistency <- NULL;rm( OverallExternalConsistency )
  OverallInternalConsistency <- NULL;rm( OverallInternalConsistency ) 
  Subgroup <- NULL;rm( Subgroup )
  TSDT_SAMPLES <- NULL;rm( TSDT_SAMPLES )
  
  if( !is.null(  scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  
  NULL_SCORES <- NULL
  while( length( NULL_SCORES ) < n_permutations ){
    
    permuted_data <- permutation( response = response,
                                 trt = trt,
                                 permute_arm = permute_arm )

    if( !is.data.frame( permuted_data ) ){
      permuted_data <- as.data.frame( permuted_data )
      names( permuted_data ) <- "response"
      
    }
    
    if( is.null( permuted_data$trt ) )
        permuted_data$trt <- trt
    
    if( !is.null( covariates ) )
        permuted_data <- cbind( permuted_data, covariates )
   
    samples <- get_samples( data = permuted_data,
                           trt = permuted_data$trt,
                           trt_control = trt_control,
                           n_samples = n_samples,
                           sampling_method = sampling_method,
                           inbag_proportion = inbag_proportion )
    
    tsdt_samples_list <- populate_tsdt_samples( samples = samples,
                                                min_subgroup_n_control = min_subgroup_n_control,
                                                min_subgroup_n_trt = min_subgroup_n_trt,
                                                min_subgroup_n_oob_control = min_subgroup_n_oob_control,
                                                min_subgroup_n_oob_trt = min_subgroup_n_oob_trt,
                                               
                                                scoring_function = scoring_function,
                                                scoring_function_parameters = scoring_function_parameters,
                                                desirable_response = desirable_response,
                                                eps = eps,
                                               
                                                inbag_score_margin = inbag_score_margin,
                                                oob_score_margin = oob_score_margin,
                                               
                                                tree_builder = tree_builder,
                                                tree_builder_parameters = tree_builder_parameters )
    
    unpack_args( tsdt_samples_list )
    
    rm( tsdt_samples_list )
    
    consistency <- compute_consistency( OverallExternalConsistency = OverallExternalConsistency,
                                        OverallInternalConsistency = OverallInternalConsistency,
                                        n_samples = n_samples )
    
    if( all( consistency$Subgroup %nin% c('Overall') ) )
        stop( "ERROR: Overall subgroup not found" )
    
    consistency <- subset( consistency, Subgroup %nin% c('Overall') )
    
    if( NROW( consistency ) > 0 ){
      NULL_SCORES <- c( NULL_SCORES,
                       max( consistency$Internal_Consistency * consistency$External_Consistency ) )
    }else{
      NULL_SCORES <- c( NULL_SCORES, 0 )
    }

    if( trace )
        cat( paste0( "Completed permutation ", n_cpu * length(NULL_SCORES), " of ", n_cpu * n_permutations, "\n" ) )
    
    rm( samples )
    rm( TSDT_SAMPLES )
    rm( OverallExternalConsistency )
    rm( OverallInternalConsistency )
    rm( consistency )
    
  }
  
  return( NULL_SCORES )  
}

#' @title get_y
#' @description Returns the response variable in the in-bag or out-of-bag data.
#' @details If the user provides a y_var parameter in the list of
#' scoring_function_parameters this function will return the variable specified
#' by that parameter. If the user specifies a y_col parameter in the list of
#' scoring_function_parameters the function returns the column in data indexed
#' by that parameter. Lastly, if data contains a variable called 'y' that
#' variable is returned. Otherwise, NULL is returned.
#' @seealso \link{get_trt}, \link{get_covariates}
#' @param data A data.frame containing in-bag or out-of-bag data
#' @param scoring_function_parameters A list of named elements containing control
#' parameters and other data required by the scoring function
#' @return Response variable (if present) or NULL.
#' @examples
#' ## Create an example data.frame
#' df <- data.frame( y <- 1:5 )
#' names( df ) <- "y"
#' df$time <- 10:14
#' df$time2 <- 20:24
#' df$event <- sample( c(0:1), size = 5, replace = TRUE )
#' df$trt <- sample( c("Control","Treatment"), size = 5, replace = TRUE )
#' df$x1 <- runif( n = 5 )
#' df$x2 <- LETTERS[1:5]
#' 
#' ## Select the y variable by name
#' get_y( df, scoring_function_parameters = list( y_var = 'y' ) )
#' 
#' ## Select the y variable by column index
#' get_y( df, scoring_function_parameters = list( y_col = 1 ) )
#' 
#' ## The default behavior works for this example because the y variable in df
#' ## is actually called y.
#' get_y( df )
#' 
#' ## If the user's data does not contain a variable called
#' ## 'y' the default behavior will fail. In this case the user must explicitly
#' ## identify the 'y' variable via one of the two previous methods.
#' names( df )[which(names(df) == "y")] <- "response" # rename the 'y' variable to 'response'
#' 
#' get_y( df )  # now default behavior fails (i.e. returns NULL)
#' 
#' get_y( df, scoring_function_parameters = list( y_var = 'response' ) ) # this works
#' @export
get_y <- function( data,
                  scoring_function_parameters = NULL ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  y_var <- NULL;rm( y_var )
  y_col <- NULL;rm( y_col )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  
  if( exists( "y_var" ) )
      response <- data[,c(y_var)]
  
  else if( exists( "y_col" ) )
      response <- data[,y_col]
  
  else if( "y" %in% names( data ) )
      response <- data$y
  
  else
      response <- NULL

  return( response )
}

#' @title get_trt
#' @description Returns the treatment variable in the in-bag or out-of-bag data.
#' @details If the user provides a trt_var parameter in the list of
#' scoring_function_parameters this function will return the variable specified
#' by that parameter. If the user specifies a trt_col parameter in the list of
#' scoring_function_parameters the function returns the column in data indexed
#' by that parameter. Lastly, if data contains a variable called 'trt' that
#' variable is returned. Otherwise, NULL is returned.
#' @seealso \link{get_y}, \link{get_covariates}
#' @param data A data.frame containing in-bag or out-of-bag data
#' @param scoring_function_parameters A list of named elements containing control
#' parameters and other data required by the scoring function
#' @return Treatment variable (if available) or NULL.
#' @examples
#' ## Create an example data.frame
#' df <- data.frame( y <- 1:5 )
#' names( df ) <- "y"
#' df$time <- 10:14
#' df$time2 <- 20:24
#' df$event <- sample( c(0:1), size = 5, replace = TRUE )
#' df$trt <- sample( c("Control","Treatment"), size = 5, replace = TRUE )
#' df$x1 <- runif( n = 5 )
#' df$x2 <- LETTERS[1:5]
#' 
#' ## Select the trt variable by name
#' get_trt( df, scoring_function_parameters = list( trt_var = 'trt' ) )
#' 
#' ## Select the trt variable by column index
#' get_trt( df, scoring_function_parameters = list( trt_col = 5 ) )
#' 
#' ## The default behavior works for this example because the trt variable in df
#' ## is actually called trt.
#' get_trt( df )
#' 
#' ## If the user's data does not contain a variable called
#' ## 'y' the default behavior will fail. In this case the user must explicitly
#' ## identify the 'y' variable via one of the two previous methods.
#' names( df )[which(names(df) == "trt")] <- "treatment" # rename the 'trt' variable to 'treatment'
#' 
#' get_trt( df )  # now default behavior fails (i.e. returns NULL)
#' 
#' get_trt( df, scoring_function_parameters = list( trt_var = 'treatment' ) ) # this works
#' @export
get_trt <- function( data,
                     scoring_function_parameters = NULL ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  trt_var <- NULL;rm( trt_var )
  trt_col <- NULL;rm( trt_col )
  
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )

  if( exists( "trt_var" ) )
      trt <- data[,c(trt_var)]
  
  else if( exists( "trt_col" ) ){
    ## cat( "got here\n" )
    ## print( data )
    ## print( trt_col )
    ## flush.console()
    
    trt <- data[,trt_col]
  }
  
  else if( "trt" %in% names( data ) ) {
    # special case: user explicitely defined trt as NULL
    if(exists("trt")) {
      if(!is.null(trt) ) {
        trt <- data$trt
      } else trt <- NULL
    }
    else {
      trt <- data$trt #trt <- NULL
    }
  } else
      trt <- NULL
  return( trt )
}

##################################################

#' @title get_covariates
#' @description Returns the covariate variables in the in-bag or out-of-bag data.
#' @details If the user provides a covariate_vars parameter in the list of
#' scoring_function_parameters this function will return the variables specified
#' by that parameter. If the user specifies a covariate_cols parameter in the list of
#' scoring_function_parameters the function returns the columns in data indexed
#' by that parameter. Otherwise, NULL is returned.
#' @seealso \link{get_y}, \link{get_trt}
#' @param data A data.frame containing in-bag or out-of-bag data
#' @param scoring_function_parameters A list of named elements containing control
#' parameters and other data required by the scoring function
#' @return A data.frame of covariates.
#' @examples
#' ## Create an example data.frame
#' df <- data.frame( y <- 1:5 )
#' names( df ) <- "y"
#' df$time <- 10:14
#' df$time2 <- 20:24
#' df$event <- sample( c(0:1), size = 5, replace = TRUE )
#' df$trt <- sample( c("Control","Treatment"), size = 5, replace = TRUE )
#' df$x1 <- runif( n = 5 )
#' df$x2 <- LETTERS[1:5]
#' 
#' ## Select the covariate variables by name
#' get_covariates( df, scoring_function_parameters = list( covariate_vars = c("x1","x2") ) )
#' 
#' ## Select the covariate variables by column index
#' get_covariates( df, scoring_function_parameters = list( covariate_cols = c(6:7) ) )
#' @export
get_covariates <- function( data,
                            scoring_function_parameters ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  covariate_vars <- NULL;rm( covariate_vars )
  covariate_cols <- NULL;rm( covariate_cols )
    
  if( !is.null( scoring_function_parameters ) )
      unpack_args( scoring_function_parameters )
  
  if( exists( "covariate_vars" ) )
      covariates <- subset( data, select = c(covariate_vars) )
  
  else if( exists( "covariate_cols" ) ){
    varnames <- names( data )[covariate_cols]
    covariates <- as.data.frame( data[,covariate_cols] )
    names( covariates ) <- varnames
  }
  
  else
      covariates <- NULL

  return( covariates )
}

#' @title distribution
#' @description Returns the distribution of values used to compute TSDT summary
#' statistics.
#' @details This function returns the distribution of all values used to compute
#' summary statistics for superior subgroups identified by the TSDT algorithm.
#' The summary statistics returned for a TSDT object include the mean
#' subgroup size, mean response value, and median value of the scoring function.
#' These statistics reported seperately for in-bag and out-of-bag data sets, and
#' also stratified by treatment arm. This function can also provide the
#' distribution of all cutpoints for a numeric splitting variable in a subgroup
#' definition.
#' @seealso \link{TSDT}, \link[TSDT]{summary-methods}
#' @param object An object of class TSDT
#' @param statistic The desired statistic distribution
#' @param subgroup The desired subgroup
#' @param subsub A subset of the subgroup
#' @return A vector containing the observed values for the specified subgroup
#' @export
#' @examples
#' set.seed(0)
#' N <- 200
#' continuous_response = runif( min = 0, max = 20, n = N )
#' trt <- sample( c('Control','Experimental'), size = N, prob = c(0.4,0.6),
#'                replace = TRUE )
#' X1 <- runif( N, min = 0, max = 1 )
#' X2 <- runif( N, min = 0, max = 1 )
#' X3 <- sample( c(0,1), size = N, prob = c(0.2,0.8), replace = TRUE )
#' X4 <- sample( c('A','B','C'), size = N, prob = c(0.6,0.3,0.1), replace = TRUE )
#' covariates <- data.frame( X1 )
#' covariates$X2 <- X2
#' covariates$X3 <- factor( X3 )
#' covariates$X4 <- factor( X4 )
#' 
#' ## Create a TSDT object
#' ex1 <- TSDT( response = continuous_response,
#'             trt = trt, trt_control = 'Control',
#'             covariates = covariates[,1:4],
#'             inbag_score_margin = 0,
#'             desirable_response = "increasing",
#'             oob_score_margin = 0,
#'             min_subgroup_n_control = 5,
#'             min_subgroup_n_trt = 5,
#'             n_sample = 5 )
#' 
#' ## Show summary statistics
#' summary( ex1 )
#' 
#' ## Get the number of subjects in each superior in-bag subgroup
#' distribution( ex1, statistic = 'Inbag_Subgroup_Size' )
#' 
#' ## Get the vector of subgroup sample sizes for a particular subgroup
#' distribution( ex1, statistic = 'Inbag_Subgroup_Size',
#'               subgroup = 'X1<xxxxx & X1>=xxxxx' )
#' 
#' ## Get the observed cutpoints for the numeric splitting variables in a subgroup
#' distribution( ex1, statistic = 'Cutpoints', subgroup = 'X1<xxxxx & X1>=xxxxx' )
#' 
#' ## If the subgroup definition has more than one numeric splitting variable you
#' ## can retrieve the numeric cutpoints for the splitting variables individually
#' distribution( ex1, statistic = 'Cutpoints', subgroup = 'X1<xxxxx & X1>=xxxxx',
#'               subsub = 'X1<xxxxx' )
#' distribution( ex1, statistic = 'Cutpoints', subgroup = 'X1<xxxxx & X1>=xxxxx',
#'               subsub = 'X1>=xxxxx' )
#' 
#' ## Valid statistic names come from the column names in the summary output. If
#' ## you are uncertain what the possible statistic values could be, you can pass
#' ## any arbitrary string as the statistic and an error message is returned
#' ## listing valid statistic values.
#' \dontrun{
#' distribution( ex1, statistic = 'Invalid_Statistic' )
#' }
#' @export
distribution <-function( object,
                         statistic,
                         subgroup = NULL,
                         subsub = NULL ){
  
  VALID_STATISTICS <- c( 'Cutpoints', sort( names( object@distributions ) ) )
  
  STATISTICS_SET <- paste( "\t\t", VALID_STATISTICS, sep = "", collapse = "\n" )
  
  ERROR_MESSAGE <- paste0( "ERROR: ", statistic, " not a valid statistic in object@distributions.\n",
                          "\tShould be one of:\n", STATISTICS_SET )
  
  if( statistic %nin% VALID_STATISTICS )
      stop( ERROR_MESSAGE ) 

  # If user requests Cutpoints call cutpoints()
  if( statistic == 'Cutpoints' )
      return( cutpoints( object, subgroup, subsub ) )

  # For other statistics return hash
  else{
  
    eval( parse( text = paste0( 'hash <- object@distributions$', statistic ) ) )
  
    if( is.null( subgroup ) )
        return( hash )
  
    #else
    subgroup <- as.character( subgroup )
    
    return( hash[[ subgroup ]] )
  }
}

flatten_parameters <- function( parameters ){

  PARAMETERS <- list()
  
  for( p1 in 1:length(parameters) ){
    
    if( !is(parameters[[p1]], 'list') ){
      PARAMETERS <- unlist( c( PARAMETERS, parameters[p1] ) )
    }else{
      
      param__ <- parameters[[p1]]
      
      for( p2 in 1:length(param__) ){
        
        if( class( param__[[p2]] ) %in% c("character","numeric") ){
          
          if( names( param__ )[p2] %nin% names( PARAMETERS ) )
              PARAMETERS[ names( param__ )[p2] ] <- param__[p2]
          
        }else{ # any other class of parameter just convert to string representation
          string__ <- paste( capture.output( str( param__[[p2]] ) ), sep = '', collapse  = '' )
          PARAMETERS[[ names( param__ )[p2] ]] <- string__
          rm( string__ )
        }
      }
      rm( param__ )
    }
  }
  
  PARAMETERS <- as.data.frame( PARAMETERS, stringsAsFactors = FALSE )
  
  return( PARAMETERS )  
}

equal <- function( x, y, verbose = FALSE ){

  # Confirm both objects are of class TSDT
  if( !is(x, "TSDT") || !is(y,"TSDT") )
      stop( "ERROR: both objects must be of class TSDT" )

  ## Confirm equivalence of object dimensions
  if( length( x@parameters ) != length( y@parameters ) ){
    if( verbose )
        cat( 'NOTE: paramater slots do have not equal length\n' )
    return( FALSE )
  }
  
  if( length( x@samples ) != length( y@samples ) ){
    if( verbose )
        cat( 'NOTE: samples slots do not have equal length\n' )
    return( FALSE )
  }

  if( !identical( dim( x@superior_subgroups ), dim( y@superior_subgroups ) ) ){
    if( verbose )
        cat( 'NOTE: superior_subgroups slots do not have equal dimensions\n' )
    return( FALSE )
  }

  if( length( x@cutpoints@Cutpoints ) != length( y@cutpoints@Cutpoints ) ){
    if( verbose )
        cat( 'NOTE: cutpoints slots do not have equal length\n' )
    return( FALSE )
  }
  
  if( length( setdiff( names( x@distributions ), names( y@distributions ) ) ) != 0 ){
    if( verbose )
        cat( 'NOTE: distributions slots do not contain identical contents\n' )
    return( FALSE )
  }

  ## Confirm equivalence of parameter contents
  xparm <- flatten_parameters( x@parameters )
  names(xparm) <- 'x'
  xparm$parm <- row.names( xparm )
  row.names( xparm ) <- NULL

  yparm <- flatten_parameters( y@parameters )
  names(yparm) <- 'y'
  yparm$parm <- row.names( yparm )
  row.names( yparm ) <- NULL
  
  if( length( setdiff( xparm$parm, yparm$parm ) ) != 0 ){
    if( verbose )
        cat( 'NOTE: parameter slots not equal\n' )
    return( FALSE )
  }

  parms <- merge( x = xparm, y = yparm, all = TRUE )
  
  if( !identical( parms$x, parms$y ) ){
    if( verbose )
        cat( 'NOTE: parameter slots do not contain identical contents\n' )
    return( FALSE )
  }

  rm( xparm, yparm, parms )
  
  ## Confirm equivalence of samples contents
  for( i in 1:length( x@samples ) ){

    if( !identical( x@samples[[i]]@inbag, y@samples[[i]]@inbag  ) ){
      if( verbose )
          cat( paste0( 'NOTE: samples[[', i, ']]@inbag slots do not contain identical contents\n' ) )
      return( FALSE )
    }
    
    if( !identical( x@samples[[i]]@oob, y@samples[[i]]@oob ) ){
      if( verbose )
          cat( paste0( 'NOTE: samples[[', i, ']]@oob slots do not contain identical contents\n' ) )
      return( FALSE )
    }
    
    if( !identical( x@samples[[i]]@subgroups, y@samples[[i]]@subgroups ) ){
      if( verbose )
          cat( paste0( 'NOTE: samples[[', i, ']]@subgroups slots do not contain identical contents\n' ) )
      return( FALSE )
    }
  }
  rm( i )
  ## Confirm equivalence of superior_subgroups contents
  if( !identical( x@superior_subgroups, y@superior_subgroups ) ){
    if( verbose )
        cat( 'NOTE: superior_subgroups slots do not contain identical contents\n' )
    return( FALSE )
  }
  
  ## Confirm equivalence of cutpoints contents
  xcutpoints_unlist <- unlist( values( x@cutpoints@Cutpoints ) )
  ycutpoints_unlist <- unlist( values( y@cutpoints@Cutpoints ) )
  
  xcutpoints <- as.data.frame( xcutpoints_unlist )
  names(xcutpoints) <- 'x'
  xcutpoints$cutpoints <- names( xcutpoints_unlist )
  rownames( xcutpoints ) <- NULL

  ycutpoints <- as.data.frame( ycutpoints_unlist )
  names(ycutpoints) <- 'y'
  ycutpoints$cutpoints <- names( ycutpoints_unlist )
  rownames( ycutpoints ) <- NULL
  
  rm( xcutpoints_unlist, ycutpoints_unlist )
  
  if( length( setdiff( xcutpoints$cutpoints, ycutpoints$cutpoints ) ) != 0 ){
    if( verbose )
        cat( 'NOTE: cutpoints slots do not contain identical contents\n' )
    return( FALSE )
  }

  cutpoints <- merge( x = xcutpoints, y = ycutpoints, all = TRUE )

  if( !identical( cutpoints$x, cutpoints$y) ){
    if( verbose )
        cat( 'NOTE: cutpoints slots do not contain identical contents\n' )
    return( FALSE )
  } 

  rm( xcutpoints, ycutpoints, cutpoints )

  ## Confirm equivalence of distributions contents
  for( j in names( x@distributions ) ){
    
    xdistribution_unlist <- unlist( values( x@distributions[[j]] ) )
    ydistribution_unlist <- unlist( values( y@distributions[[j]] ) )

    
    xdistribution <- as.data.frame( xdistribution_unlist )
    rownames( xdistribution ) <- NULL
    
    ydistribution <- as.data.frame( ydistribution_unlist )
    rownames( ydistribution ) <- NULL

    # If distribution extracts to a matrix
    if( is( xdistribution_unlist, "matrix" ) ){

      if( length( setdiff( names( xdistribution ), names( ydistribution ) ) ) != 0 ){
        if( verbose )
            cat( 'NOTE: distributions slots do not contain identical contents 0\n' )
        return( FALSE )
      }
      
      for( d in names( xdistribution ) ){
        
        if( !identical( xdistribution[,d], ydistribution[,d] ) ){
          if( verbose )
            cat( 'NOTE: distributions slots do not contain identical contents 1\n' )
          return( FALSE )
        }
      }
      
      rm( d )
    }
    # If distribution extracts to a vector
    if( is( xdistribution_unlist,"numeric" ) ){
      
      names(xdistribution) <- 'x'
      xdistribution$distributions <- names( xdistribution_unlist )

      names(ydistribution) <- 'y'
      ydistribution$distributions <- names( ydistribution_unlist )

      if( length( setdiff( xdistribution$distributions, ydistribution$distributions ) ) != 0 ){
        if( verbose )
            cat( 'NOTE: distributions slots do not contain identical contents\n' )
        return( FALSE )
      }

      distributions <- merge( x = xdistribution, y = ydistribution, all = TRUE )

      if( !identical( distributions$x, distributions$y ) ){
        if( verbose )
            cat( 'NOTE: distributions slots do not contain identical contents\n' )
        return( FALSE )
      }
      rm( distributions ) 
    }
    
    rm( xdistribution, ydistribution ) 
  }

  rm( j )

  return(TRUE)
}


## END OF FILE
