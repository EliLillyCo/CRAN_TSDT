#!usr/bin/dev/R
#################################################################################
# FILENAME   : TSDT.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 01/24/2013
# DESCRIPTION: R implementation of TSDT subgroup identification tool.
#################################################################################

#' @title Treatment-Specific Subgroup Detection Tool
#' 
#' @description
#' Implements a method for identifying subgroups with superior response relative
#' to the overall sample.
#'  
#' @details
#' The Treatment-Specific Subgroup Detection Tool (TSDT) creates several
#' bootstrapped samples from the input data. For each of these bootstrapped
#' samples the in-bag and out-of-bag data are retained. A tree is grown on the
#' in-bag data of each bootstrapped sample using the response variable and
#' supplied covariates. Each split in the tree defines a subgroup. The overall
#' mean response for the in-bag data is computed as well as the mean response
#' within each subgroup. Additionally, a scoring function is provided. Example
#' scoring functions might be mean response, difference in mean response between
#' treatment arms (i.e. treatment effect), or a quantile of the response (e.g.
#' median), or a difference in quantiles across treatment arms. Sensible
#' defaults are provided given the data type of the response and treatment
#' variables. The user can also specify a custom scoring function. The value of
#' the scoring function is computed for the overall in-bag data and each
#' subgroup. Subgroups with mean response larger than the overall in-bag mean
#' response and a mean scoring function value larger than the overall in-bag
#' scoring function value are identified as superior subgroups. This definition
#' of a superior subgroup assumes a larger value of the response variable is
#' desirable. If a smaller value of the response is desirable then subgroups
#' with mean response and mean scoring function smaller than the overall in-bag
#' mean are superior. The same computation of overall and subgroup mean
#' response and mean scoring function are done for the out-of-bag data. This is
#' repeated for all bootstrapped samples.  Measures of internal and external
#' consistency are then computed. Internal consistency is computed for each
#' subgroup that is identified as superior in one of the in-bag samples.
#' Internal consistency for each of these subgroups is the fraction of
#' bootstrapped samples where that subgroup is identified as superior in the
#' in-bag data. External consistency is also defined only for subgroups that
#' are identified as superior in at least one of the in-bag samples. For each
#' of these subgroups, external consistency is the number of bootstrapped
#' samples where the subgroup is defined as superior in the in-bag and
#' out-of-bag data divided by the number of bootstrapped samples where the
#' subgroup is identified as superior in the in-bag data. The internal and
#' external consistency results are returned for each subgroup that identified
#' as superior in the in-bag data of at least one bootstrapped sample. A score
#' for the overall strength of each subgroup is computed as the product of the
#' internal and external consistency. Optionally, a permutation-adjusted
#' p-value for the strength of each subgroup can be computed. Based on this
#' p-value subgroups are classified as strong, moderate, weak, or not
#' confirmed. A suggested cutoff for each subgroup is also provided. This is
#' helpful because two subgroups defined on the same continuous splitting
#' variable but with different cutpoints are considered equivalent. That is,
#' one subgroup X1<0.6 and another X1<0.7 would be considered equivalent and
#' listed in the results as X1<xxxxx. (Note that X1<0.6 and X1>=0.7 would be
#' considered distinct subgroups and listed in the output as X1<xxxxx and
#' X1>=xxxxx, respectively.) So if a subgroup listed in the output as X1<xxxxx
#' could actually represent many different numeric values for xxxxx it is
#' helpful to provide a final suggestion for the cutpoint. The algorithm
#' retains all the numeric values and uses the median as the suggested cutoff.
#' The user can also request the vector of numeric cutpoints and use any
#' function of their choosing to compute a suggested cutoff.
#' @param response Response variable.
#' @param response_type Data type of response. Must be one of binary, continuous,
#' survival. If none provided it will be inferred from the data type of response.
#' (optional)
#' @param survival_model The model to use for a survival response. Defaults to
#' kaplan-meier. Other possible values are: coxph, fleming-harrington, fh2,
#' weibull, exponential, gaussian, logistic, lognormal, and loglogistic.
#' (optional)
#' @param percentile For a two-arm study this parameter specifies a test for
#' the difference in response percentile across the two treatment arms. For a
#' continuous response the default value for percentile is NULL. Instead, the
#' difference in mean response is computed by default for a continuous response.
#' If the user provides a values of percentile = 0.50 then the difference in
#' median response would be computed. For a survival outcome, the default value
#' for percentile is 0.50, which computes the difference in median survival.
#' @param tree_builder The algorithm to use for building the trees. Defaults to
#' rpart. Other possible values include ctree and mob (both from the party
#' package). (optional)
#' @param tree_builder_parameters A named list of parameters to pass to the
#' tree-builder. The default tree-builder is rpart. In this case, the parameters
#' passed here would be rpart parameters. Examples might include parameters such
#' as control, cost, weights, na.action, etc. Consult the rpart documentation (or
#' the documentation of your selected tree-builder) for a complete list. (optional)
#' @param covariates A data.frame containing the covariates.
#' @param trt Treatment variable. Only needed if there are two treatment arms.
#' (optional)
#' @param trt_control Value for treatment control arm. This parameter is relevant
#' only for two-arm data. (defaults to 0)
#' @param permute_method Indicates whether only the response variable should be
#' permuted in the computation of the p-value, or the response and treatment
#' variable should be permuted together (preserving the treatment-response
#' correlation, but eliminating the correlation with the covariates), or the
#' response variable should be permuted within one treatment arm only. The
#' parameter values for these permutation schemes are (respectively) simple,
#' permute_response_and_treatment, and permute_response_one_arm. See permute_arm
#' to specify which treatment arm is to be permuted. The default permutation
#' scheme is response_one_arm. As noted in the documentation for the permute_arm
#' parameter is to permute the non-control arm. Taken together, this implies the
#' default permutation method for p-value computation is to permute the response
#' in the non-control arm only. For one-arm data only the response is permuted.
#' (optional)
#' @param permute_arm Which treatment arm should be permuted? Defaults to the
#' experimental treatment arm -- i.e. the treatment arm not matching the value
#' provided in trt_control. For one-arm data only the response is permuted.
#' (optional)
#' @param n_samples Number of TSDT_Samples to draw.
#' @param desirable_response Direction of desirable response. Valid values are
#' 'increasing' or 'decreasing'. The default value is 'increasing'. It is
#' important to note that although the parameter is called desirable_response, it
#' actually refers to the desirable direction of scoring function values. In most
#' cases there is a positive correlation bewteen the response and scoring
#' function values -- i.e. as the response increases the scoring function also
#' increases. One instance for which this relationship between response and
#' scoring function may not hold is when mean_deviance_residuals or
#' diff_mean_deviance_residuals is used as the scoring function. See the help for
#' these scorings function for further details.
#' @param sampling_method Sampling method used to populate samples for TSDT
#' in-bag and out-of-bag data. Must be either bootstrap or subsample. Default
#' is bootstrap.
#' @param inbag_proportion The proportion of the data to use as the in-bag
#' subset when sampling_method is subsample.
#' @param scoring_function Scoring function to compute treatment effect. Links to
#' several possible scoring functions are provided in the See Also section below.
#' @param scoring_function_parameters Parameters passed to the scoring function.
#' As an example, the scoring function quantile_response takes a parameter
#' "percentile" which indicates the desired percentile of the response distribution.
#' Thus, if the median response is desired, this parameter could be set as follows:
#' scoring_function_parameters = list( percentile = 0.50 ). Most of the built-in
#' scoring functions have sensible defaults for the scoring function parameters
#' so it is not necessary to specify them explicitly in the call to TSDT. But this
#' parameter could be very useful for user-defined custom scoring functions.
#' (optional)
#' @param inbag_score_margin Required margin above overall mean for a subgroup to
#' be considered superior. If a subgroup mean must be 10\% larger than the overall
#' subgroup mean to be superior then inbag_score_margin = 0.10. If
#' desirable_response = "decreasing" then inbag_score_margin should be negative
#' or zero.
#' @param oob_score_margin Similar to inbag_score_margin but for classifying
#' out-of-bag subgroups as superior.
#' @param eps Tolerance value for floating-point precision. The default is 1E-5. (optional)
#' @param min_subgroup_n_control Minimum number of Control arm observations in an
#' in-bag subgroup. A value greater than or equal to one will be interpreted as
#' the required minimum number of observations. A value between zero and one
#' will be interpreted as a proportion of the in-bag Control observations. For a
#' bootstrapped in-bag sample the default for this parameter is 10% of the number
#' of Control observations in the overall sample. For an in-bag sample obtained
#' via subsampling the default value is the inbag_proportion times 10% of the
#' number of Control observations in the overall sample.
#' @param min_subgroup_n_trt Minimum number of Experimental arm observations in
#' an in-bag subgroup. A value greater than or equal to one will be interpreted
#' as the required minimum number of observations. A value between zero and one
#' will be interpreted as a proportion of the in-bag Experimental observations.
#' For a bootstrapped in-bag sample the default for this parameter is 10% of the
#' number of Experimental observations in the overall sample. For an in-bag
#' sample obtained via subsampling the default value is the inbag_proportion
#' times 10\% of the number of Experimental observations in the overall sample.
#' @param min_subgroup_n_oob_control Minimum number of Control arm observations
#' in an out-of-bag subgroup. A value greater than or equal to one will be
#' interpreted as the required minimum number of observations. A value between
#' zero and one will be interpreted as a proportion of the out-of-bag Control
#' observations. For a bootstrapped out-of-bag sample the default for this
#' parameter is exp(-1)*10\% of the number of Control observations in the overall
#' sample. For an out-of-bag sample obtained via subsampling the default value
#' is the inbag_proportion times (1-inbag_proportion)*10% of the number of
#' Control observations in the overall sample.
#' @param min_subgroup_n_oob_trt Minimum number of Experimental arm observations
#' in an out-of-bag subgroup. A value greater than or equal to one will be
#' interpreted as the required minimum number of observations. A value between
#' zero and one will be interpreted as a proportion of the out-of-bag
#' Experimental observations. For a bootstrapped out-of-bag sample the default
#' for this parameter is exp(-1)*10\% of the number of Experimental observations
#' in the overall sample. For an out-of-bag sample obtained via subsampling the
#' default value is the inbag_proportion times (1-inbag_proportion)*10\% of the
#' number of Experimental  observations in the overall sample.
#' @param maxdepth Maximum depth of trees.
#' @param rootcompete Number of competitor splits to retain for root node split.
#' @param competedepth Depth of competitor split trees (defaults to 1)
#' @param strength_cutpoints Cutpoints for permuted p-values to classify a
#' subgroup as Strong, Moderate, Weak, or Not Confirmed. The default cutpoints are
#' 0.10, 0.20, and 0.30 for Strong, Moderate, and Weak subgroups, respectively. (optional)
#' @param n_permutations Number of permutations to compute for adjusted p-value.
#' Defaults to zero (no p-value computation). If p-values are desired, it is
#' recommended to use at least 500 permutations.
#' @param n_cpu Number of CPUs to use. Defaults to 1.
#' @param trace Report number of permutations computed as algorithm proceeds.
#' @seealso \link{mean_response}, \link{quantile_response},
#' \link{diff_quantile_response}, \link{treatment_effect},
#' \link{survival_time_quantile}, \link{diff_survival_time_quantile},
#' \link{mean_deviance_residuals}, \link{diff_mean_deviance_residuals},
#' \link{diff_restricted_mean_survival_time}, \linkS4class{TSDT},
#' \link[rpart]{rpart}, \link[party]{ctree}, \link[party]{mob}
#' @return An object of class \linkS4class{TSDT}
#' @author Brian Denton \email{denton_brian_david@@lilly.com},
#' Chakib Battioui \email{battioui_chakib@@lilly.com},
#' Lei Shen \email{shen_lei@@lilly.com}
#' @references
#' Battioui, C., Shen, L., Ruberg, S., (2014). A Resampling-based Ensemble Tree Method to Identify Patient Subgroups with Enhanced Treatment Effect. JSM proceedings, 2014
#' 
#' Shen, L., Battioui, C., Ding, Y., (2013). Chapter "A Framework of Statistical methods for Identification of Subgroups with Differential Treatment Effects in Randomized Trials" in the book "Applied Statistics in Biomedicine and Clinical Trials Design" 
#' @examples
#' ## Create example data for constructing TSDT object
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
#'
#' ## In the following examples n_samples and n_permutations are set to small
#' ## values so the examples complete quickly. The intent here is to provide
#' ## a small functional example to demonstrate the structure of the output. In
#' ## a real-world use of TSDT these values should be at least 100 and 500,
#' ## respectively.
#' 
#' ## Single-arm TSDT
#' ex1 <- TSDT( response = continuous_response,
#'             covariates = covariates[,1:4],
#'             inbag_score_margin = 0,
#'             desirable_response = "increasing",
#'             n_samples = 5,       ## use value >= 100 in real world application
#'             n_permutations = 5,  ## use value >= 500 in real world application
#'             rootcompete = 1,
#'             maxdepth = 2 )
#' 
#' ## Two-arm TSDT
#' ex2 <- TSDT( response = continuous_response,
#'             trt = trt, trt_control = 'Control',
#'             covariates = covariates[,1:4],
#'             inbag_score_margin = 0,
#'             desirable_response = "increasing",
#'             oob_score_margin = 0,
#'             min_subgroup_n_control = 10,
#'             min_subgroup_n_trt = 20,
#'             maxdepth = 2,
#'             rootcompete = 1,
#'             n_samples = 5,      ## use value >= 100 in real world application
#'             n_permutations = 5 ) ## use value >= 500 in real world application
#' @export
TSDT <- function( response = NULL,
                  response_type = NULL,
                  survival_model = "kaplan-meier",
                  percentile = 0.50,
                  tree_builder = "rpart",
                  tree_builder_parameters = list(),
                  covariates,
                  trt = NULL,
                  trt_control = 0,
                  permute_method = NULL,
                  permute_arm = NULL,
                  n_samples = 1,
                  desirable_response = NULL,
                  sampling_method = 'bootstrap',
                  inbag_proportion = 0.5,
                  scoring_function = NULL,
                  scoring_function_parameters = list(),
                  inbag_score_margin = 0,
                  oob_score_margin = 0,
                  eps = 1E-5,
                  min_subgroup_n_control = NULL,
                  min_subgroup_n_trt = NULL,
                  min_subgroup_n_oob_control = NULL,
                  min_subgroup_n_oob_trt = NULL,
                  maxdepth = 30,
                  rootcompete = 0,
                  competedepth = 1,
                  strength_cutpoints = c(0.10,0.20,0.30),
                  n_permutations = 0,
                  n_cpu = 1,
                  trace = FALSE ){

  ## Avoid NOTE in R CMD check about no visible binding for global variable
  PID0 <- scoring_function_name <- Subgroup <- NULL
  
  if( tree_builder %nin% c('rpart','ctree','mob') )
      stop( 'ERROR: tree_builder must be one of {rpart,ctree,mob}' )

  ## validate covariate names
  invalid_names <- NULL
  for( k in names( covariates ) ){
    chars <- strsplit( split = "", k )[[1]]
    if( any( chars %nin% c( LETTERS, letters, 0:9, '.', '_' ) ) ){
      invalid_names <- c( invalid_names, k )
    }
    rm( chars )
  }
  
  if( !is.null( invalid_names ) ){
    msg__ <- paste0( 'ERROR: TSDT requires alphanumeric covariate names. The ',
                    'following covariate names are invalid:\n' )
    
    msg__ <- paste0( msg__, paste( invalid_names, sep = "", collapse = "\n" ) )
    
    stop( msg__ )
  }

  
  if( tree_builder == 'rpart'){
    requireNamespace( "rpart", quietly = TRUE )
  }else if( tree_builder %in% c('ctree','mob') ){
    requireNamespace( "party", quietly = TRUE )
  }
  
  if( n_cpu > 1 ){
    requireNamespace( "parallel", quietly = TRUE )
    RNGkind( kind = "L'Ecuyer-CMRG" )
    parallel::mc.reset.stream()
  }
  
  if( n_cpu > n_samples )
      stop( "ERROR: n_cpu must be less than or equal to n_samples" )

  # If trt provided make sure trt_control value exists in data
  if( !is.null( trt ) && !any( trt == trt_control ) )
      stop( "ERROR: trt_control value not found" )

  ## If trt is NULL make sure other trt-related variables are also NULL
  if( is.null( trt ) ){
    trt_control <- NULL
    permute_arm <- NULL
  }
  
  # Remove records with missing response
  if( any( is.na( response ) ) ){
    covariates <- subset( covariates, !is.na( response ) )
    
    if( !is.null( trt ) )
        trt <- subset( trt, !is.na( response ) )

    if( "weights" %in% names( tree_builder_parameters ) )
        tree_builder_parameters$weights <- subset( tree_builder_parameters$weights, !is.na( response ) )


    if( "cost" %in% names( tree_builder_parameters ) ){
      if( is.null( names( tree_builder_parameters$cost ) ) )
          names( tree_builder_parameters$cost ) < names( covariates )
    }else{
      tree_builder_parameters$cost <- rep( 1, NCOL( covariates ) )
      names( tree_builder_parameters$cost ) < names( covariates )
    }
    
    response <-  subset( response, !is.na( response ) )
  }
  
  # Convert spaces in text/factor covariates to underscores
  covariates <- space2underscore( covariates )
  
  # Populate response_type if none provided
  if( is.null( response_type ) ){
    
    if( is( response, "Surv" ) )
        response_type <- "survival"
    
    else if( ( is(response, "character") || is( response, "factor" ) ) && is.binary( response ) )
        response_type <- "binary"
    
    else if( is( response, "numeric" ) )
        response_type <- "continuous"
    
    else
        stop( "ERROR: response_type must be one of {binary,continuous,survival}" )
  }

  if( response_type == "binary" )
      response <- binary_transform( response )
  
  if( response_type == "survival" )
      requireNamespace( "survival", quietly = TRUE )


  ## Populate default desirable_response
  if( is.null( desirable_response ) ){
    if( !is.null( scoring_function) && substitute( scoring_function ) %in% c('mean_deviance_residuals','diff_mean_deviance_residuals') ){
      desirable_response <- 'decreasing'
      cat( "NOTE: setting desirable_response = 'decreasing'\n" )
      cat( paste0( " See help for ", substitute( scoring_function ), '\n' ) )
    }else{
      desirable_response <- 'increasing'
    }
  }
  
  scoring_function_parameters$response_type <- response_type
  
  # Construct a data.frame from provided data and populate relevant scoring
  # function parameters
  source_data <- data.frame( response )
  names( source_data ) <- "response"
  scoring_function_parameters$y_var <- "response"
      
  if( response_type == "survival" ){
    scoring_function_parameters$survival_model <- survival_model
    scoring_function_parameters$percentile <- percentile
  }
  
  if( !is.null( trt ) ){
    source_data$trt <- trt
    scoring_function_parameters$trt_var <- "trt"
    scoring_function_parameters$trt_control <- trt_control
  }
  
  scoring_function_parameters$desirable_response <- desirable_response
    
  covariates <- as.data.frame( covariates )  
  source_data <- cbind( source_data, covariates )
  COVARS <- names( covariates )
  scoring_function_parameters$covariate_vars <- COVARS

  
  if( !is.null( scoring_function ) && !exists( "scoring_function_parameters$scoring_function_name" ) )
      scoring_function_parameters$scoring_function_name <- deparse( substitute( scoring_function ) )
  
  
  if( is.null( permute_method ) )
      permute_method <- "permute_response_one_arm"
  
  else{
    if( !is.character( permute_method ) )
        permute_method <- as.character( substitute( permute_method ) )
    
    if( permute_method %nin% c( "permute_response_and_treatment",
                                "permute_response_one_arm",
                                "simple_permute" ) )
        stop( "ERROR: permute_method must be one of {permute_response_and_treatment,permute_response_one_arm,simple_permute}" )
  }
  
  if( !is.null( trt ) &&
     permute_method %in% c("permute_response_and_treatment", "permute_response_one_arm" ) ){
    
    if( is.null( permute_arm ) )
        permute_arm <- treatment_value( trt, trt_control )
    
    if( permute_arm %nin% trt )
        stop( "ERROR: permute_arm value not found in trt" )
  }

  if( rootcompete > NCOL( covariates ) - 1 )
      stop( "ERROR: rootcompete must be less than or equal to the number of covariates minus one." )

  if( sampling_method %nin% c('bootstrap','subsample') )
      stop( 'ERROR: sampling_method must be one of {bootstrap,subsample}' )

  ## Reset factor levels
  covariates <- droplevels( covariates )
  if( is( response, 'factor' ) ){
    response <- droplevels( response )
  }

  if( !is.null( trt ) && is( trt, 'factor' ) ){
    trt <- droplevels( trt )
  }
  
  ###############################################################################
  # Populate/validate control parameters                                        #
  ###############################################################################
  # Populate defualt values for min_subgroup_n_control, min_subgroup_n_trt,
  # min_subgroup_n_oob_control, min_subgroup_n_oob_trt when these parameters
  # are NULL.

  if( !is.null( trt ) ){
    n_control <- length( trt[ trt == trt_control ] )
    n_trt <- length( trt[ trt != trt_control ] )
  }else{
    n_control <- 0
    n_trt <- NROW( response )
  }
    
  if( sampling_method == 'bootstrap' ){

    # Set default bootstrapped in-bag subgroup size to 10% of in-bag control
    # and treatment arms
    if( is.null( min_subgroup_n_control ) )
        min_subgroup_n_control <- ceiling( 0.10 * n_control )
    
    if( is.null( min_subgroup_n_trt ) )
        min_subgroup_n_trt <- ceiling( 0.10 * n_trt )
    
    # Set default bootstrappped out-of-bag subgroup sizes to
    # exp(-1) * in-bag subgroup sizes for control and treatment arms 
    if( is.null( min_subgroup_n_oob_control ) )
        min_subgroup_n_oob_control <- ceiling( exp(-1) * min_subgroup_n_control )
    
    if( is.null( min_subgroup_n_oob_trt ) )
        min_subgroup_n_oob_trt <- ceiling( exp(-1) * min_subgroup_n_trt )
    
  } # END bootstrap
  
  if( sampling_method == 'subsample' ){ # BEGIN subsample
    
    # Set default subsampled in-bag subgroup sizes to 10% of in-bag sample size
    # for control and treament arms
    if( is.null( min_subgroup_n_control ) )
        min_subgroup_n_control <- ceiling( 0.10 * inbag_proportion * n_control )
    
    if( is.null( min_subgroup_n_trt ) )
        min_subgroup_n_trt <- ceiling( 0.10 * inbag_proportion * n_trt )
    
    # Set default subsampled out-of-bag subgroup sizes to 10% of out-of-bag
    # sample size for control and treament arms
    if( is.null( min_subgroup_n_oob_control ) )
        min_subgroup_n_oob_control <- ceiling( 0.10 * ( 1 - inbag_proportion ) * n_control )
    
    if( is.null( min_subgroup_n_oob_trt ) )
        min_subgroup_n_oob_trt <- ceiling( 0.10 * ( 1 - inbag_proportion ) * n_trt ) 
  }# END subsample
  
  # Convert minimum subgroup sizes expressed as proportions to absolute numbers.
  if( sampling_method == 'bootstrap' ){
    
    if( min_subgroup_n_control < 1 )
      min_subgroup_n_control <- min_subgroup_n_control * n_control

    if( min_subgroup_n_trt < 1 )
        min_subgroup_n_trt <- min_subgroup_n_trt * n_trt

    if( min_subgroup_n_oob_control < 1 )
        min_subgroup_n_oob_control <- exp(-1) * min_subgroup_n_oob_control * n_control

    if( min_subgroup_n_oob_trt < 1 )
        min_subgroup_n_oob_trt <- exp(-1) * min_subgroup_n_oob_trt * n_trt
    
  }else if( sampling_method == 'subsample' ){
    
    if( min_subgroup_n_control < 1 )
        min_subgroup_n_control <- inbag_proportion * min_subgroup_n_control * n_control

    if( min_subgroup_n_trt < 1 )
        min_subgroup_n_trt <- inbag_proportion * min_subgroup_n_trt * n_trt

    if( min_subgroup_n_oob_control < 1 )
        min_subgroup_n_oob_control <- ( 1 - inbag_proportion ) * min_subgroup_n_oob_control * n_control

    if( min_subgroup_n_oob_trt < 1 )
        min_subgroup_n_oob_trt <- ( 1 - inbag_proportion ) * min_subgroup_n_oob_trt * n_trt
  }

  # Use ceiling of computed subgroup sizes
  min_subgroup_n_control <- ceiling( min_subgroup_n_control )
  min_subgroup_n_trt <- ceiling( min_subgroup_n_trt )
  min_subgroup_n_oob_control <- ceiling( min_subgroup_n_oob_control )
  min_subgroup_n_oob_trt <- ceiling(  min_subgroup_n_oob_trt )


  if( !is.null( trt ) ){
    if( min_subgroup_n_control > length( which( trt == trt_control ) ) )
        stop( 'ERROR: min_subgroup_n_control cannot be larger than the total number of control subjects' )
    if( min_subgroup_n_oob_control > length( which( trt == trt_control ) ) )
        stop( 'ERROR: min_subgroup_n_oob_control cannot be larger than the total number of control subjects' )
    
    if( min_subgroup_n_trt > length( which( trt != trt_control ) ) )
        stop( 'ERROR: min_subgroup_n_trt cannot be larger than the total number of trt subjects' )
    if( min_subgroup_n_oob_trt > length( which( trt != trt_control ) ) )
        stop( 'ERROR: min_subgroup_n_oob_trt cannot be larger than the total number of trt subjects' )
 
  }else{
    if( min_subgroup_n_control > NROW( response ) )
        stop( 'ERROR: min_subgroup_n_control cannot be larger than the total number of subjects' )
    if( min_subgroup_n_oob_control > NROW( response ) )
        stop( 'ERROR: min_subgroup_n_oob_control cannot be larger than the total number of subjects' )

    if( min_subgroup_n_trt > NROW( response ) )
        stop( 'ERROR: min_subgroup_n_trt cannot be larger than the total number of subjects' )
    if( min_subgroup_n_oob_trt > NROW( response ) )
        stop( 'ERROR: min_subgroup_n_oob_trt cannot be larger than the total number of subjects' )
  }
            
  scoring_function_parameters$min_subgroup_n_control <- min_subgroup_n_control
  scoring_function_parameters$min_subgroup_n_trt <- min_subgroup_n_trt

  scoring_function_parameters$min_subgroup_n_oob_control <- min_subgroup_n_oob_control
  scoring_function_parameters$min_subgroup_n_oob_trt <- min_subgroup_n_oob_trt
  
  # Get scoring function and associated parameters (if any)
  scoring_function_parms <- get_scoring_function( scoring_function = scoring_function,
                                                  data = source_data,
                                                  scoring_function_parameters = scoring_function_parameters )

  # Unpack the list populated in the previous step
  unpack_args( scoring_function_parms )
  
  scoring_function_parameters$scoring_function_name <- scoring_function_name
    
  tree_builder_parameters$maxdepth <- maxdepth
  tree_builder_parameters$rootcompete <- rootcompete
  tree_builder_parameters$competedepth <- competedepth
  
  # Set parameters specific to each tree-building algorithm
  if( tree_builder == "rpart" ){
    
    if( "control" %nin% names( tree_builder_parameters ) )
        tree_builder_parameters$control <- rpart.control()
          
    tree_builder_parameters$control$maxdepth <- maxdepth

    # competitor splits are obtained via get_competitors(), not rpart's
    # internal mechanism for retaining competitor splits
    tree_builder_parameters$control$maxcompete <- 0
          
    tree_builder_parameters$control$cp <- 0 # grow largest possible tree
    tree_builder_parameters$control$xval <- 0 # do not use cross-validation
    
    tree_builder_parameters$control$minbucket <- min_subgroup_n_trt
    tree_builder_parameters$control$minsplit <- 2 * tree_builder_parameters$control$minbucket
    
    
    if( "weights" %in% names( tree_builder_parameters ) ){
      source_data$WEIGHTS__ <- tree_builder_parameters$weights
      scoring_function_parameters$covariate_vars <- c( scoring_function_parameters$covariate_vars, "WEIGHTS__" )
    }
    
  }else if( tree_builder == "ctree" ){
    
    if( "controls" %nin% names( tree_builder_parameters ) )
        tree_builder_parameters$controls <- party::ctree_control()
    
    tree_builder_parameters$controls@tgctrl@maxdepth <- as.integer( maxdepth )
    
    tree_builder_parameters$controls@gtctrl@mincriterion <- 0 # Grow largest possible tree
    tree_builder_parameters$controls@splitctrl@minbucket <- min_subgroup_n_trt
    tree_builder_parameters$controls@splitctrl@minsplit <- 2*tree_builder_parameters$controls@splitctrl@minbucket
    
  }else if( tree_builder == "mob" ){
    
    
    if( "control" %nin% names( tree_builder_parameters ) )
        tree_builder_parameters$control <- mob_control()
    
    tree_builder_parameters$control$minsplit <- 2 * min_subgroup_n_trt
    tree_builder_parameters$control$alpha <- 1 # Grow large tree (do not test for parameter instability, see mob)
    tree_builder_parameters$control$bonferroni <- FALSE
    
  }else{
    stop( "ERROR: tree_builder must be one of {rpart,ctree,mob}" )
  }
  
  # Validate inbag_proportion for subsample sampling method
  if( sampling_method == "subsample" ){
    if( inbag_proportion < 0 || inbag_proportion > 1 ){
      stop( "ERROR: inbag_proportion must be in [0,1]" )
    }
  }

  ###############################################################################
  # Convert logical inequalities for numeric variables (i.e. <, >, <=, and >=)  #
  # that are embedded within character variable values to a character           #
  # representation. For example, if a character variable X1 has the value:      #
  # "Subjects with temperature >= 99 degrees", this will be converted to        #
  # "Subjects with temperature __GE__ 99 degrees".                              #
  ###############################################################################

  for( j in 1:NCOL( source_data ) ){
    if( is( source_data[,j], 'character') ){
      source_data[,j] <- gsub( pattern = '>=', replacement = '%%__GE__%%', fixed = TRUE, source_data[,j] )
      source_data[,j] <- gsub( pattern = '<=', replacement = '%%__LE__%%', fixed = TRUE, source_data[,j] )
      source_data[,j] <- gsub( pattern = '>',  replacement = '%%__GT__%%', fixed = TRUE, source_data[,j] )
      source_data[,j] <- gsub( pattern = '<',  replacement = '%%__LT__%%', fixed = TRUE, source_data[,j] )
      source_data[,j] <- gsub( pattern = '=',  replacement = '%%__EQ__%%', fixed = TRUE, source_data[,j] )

      source_data[,j] <- ifelse( grepl( pattern = '^\\s*$', source_data[,j] ), '%%__EMPTY_STRING__%%', source_data[,j] )
    }
  }
  
  ###############################################################################
  # Generate bootstrapped or subsampled subsets for constructing TSDT in-bag    #
  # and out-of-bag data sets.                                                   #
  ###############################################################################
  
  samples <- get_samples( data = source_data,
                          trt = source_data$trt,
                          trt_control = trt_control,
                          n_samples = n_samples,
                          sampling_method = sampling_method,
                          inbag_proportion = inbag_proportion )
  
  ###############################################################################
  # Populate vector of TSDT_Sample objects                                      #
  ###############################################################################
 
  if( n_cpu == 1 ){
    
    tsdt_samples_list <- populate_tsdt_samples( samples = samples,
                                                response_type = response_type,
                                                trt = trt,
                                                trt_control = trt_control,
                                                tree_builder = tree_builder,
                                                tree_builder_parameters = tree_builder_parameters,
                                                min_subgroup_n_control = min_subgroup_n_control,
                                                min_subgroup_n_trt = min_subgroup_n_trt,
                                                min_subgroup_n_oob_control = min_subgroup_n_oob_control,
                                                min_subgroup_n_oob_trt = min_subgroup_n_oob_trt,
                                                scoring_function = scoring_function,
                                                scoring_function_parameters = scoring_function_parameters,
                                                desirable_response = desirable_response,
                                                eps = eps,
                                                inbag_score_margin = inbag_score_margin,
                                                oob_score_margin = oob_score_margin )
    
    unpack_args( tsdt_samples_list )
    rm( tsdt_samples_list )
  }
  
  else if( n_cpu > 1 ){

    requireNamespace( "parallel", quietly = TRUE )
    
    SAMPLES_PARTITION <- partition( 1:length(samples), n = n_cpu )
    
    PID0 <<- Sys.getpid()
    
    for( p in 1:length( SAMPLES_PARTITION ) ){

      parallel::mcparallel({
        
        tsdt_samples_list <- populate_tsdt_samples( samples = samples[ SAMPLES_PARTITION[[p]] ],
                                                   response_type = response_type,
                                                   trt = trt,
                                                   trt_control = trt_control,
                                                   tree_builder = tree_builder,
                                                   tree_builder_parameters = tree_builder_parameters,
                                                   min_subgroup_n_control = min_subgroup_n_control,
                                                   min_subgroup_n_trt = min_subgroup_n_trt,
                                                   min_subgroup_n_oob_control = min_subgroup_n_oob_control,
                                                   min_subgroup_n_oob_trt = min_subgroup_n_oob_trt,
                                                   scoring_function = scoring_function,
                                                   scoring_function_parameters = scoring_function_parameters,
                                                   desirable_response = desirable_response,
                                                   eps = eps,
                                                   inbag_score_margin = inbag_score_margin,
                                                   oob_score_margin = oob_score_margin )
      })

    }

    JOBS <- parallel::mccollect()
    
    TSDT_SAMPLES <- NULL
    OverallExternalConsistency <- data.frame( Subgroup = character(0),
                                            OOB_Superior_Subgroup = logical(0) )
    OverallInternalConsistency <- NULL
    CutpointDistribution <- new( "TSDT_CutpointDistribution" )

    for( job in JOBS ){
      if( length( job ) ==  4 ){
        TSDT_SAMPLES <- c( TSDT_SAMPLES, job$TSDT_SAMPLES )
        OverallExternalConsistency <- rbind( OverallExternalConsistency, job$OverallExternalConsistency )
        OverallInternalConsistency <- c( OverallInternalConsistency, job$OverallInternalConsistency )
        
        append_cutpoints( CutpointDistribution, job$CutpointDistribution )
      }
    }
  }
  
  ###############################################################################
  # Compute overall external and internal consistency of identified superior    #
  # subgroups.                                                                  #
  ###############################################################################
      
  superior_subgroups <- compute_consistency( OverallExternalConsistency = OverallExternalConsistency,
                                      OverallInternalConsistency = OverallInternalConsistency,
                                            n_samples = n_samples )

  ###############################################################################
  # For each subgroup in superior_subgroups collect across all bootstrapped     #
  # samples:                                                                    #
  #   1. NodeSize (for treatment and control)                                   #
  #   2. Mean_Inbag_Response, Mean_OOB_Response (for treatment and control)     #
  #   3. Inbag_Score, OOB_Score (for treatment and control)                     #
  ###############################################################################

  Inbag_Control_Subgroup_Size_Hash <- hash()
  Inbag_Treatment_Subgroup_Size_Hash <- hash()
  Inbag_Subgroup_Size_Hash <- hash()
  
  OOB_Control_Subgroup_Size_Hash <- hash()
  OOB_Treatment_Subgroup_Size_Hash <- hash()
  OOB_Subgroup_Size_Hash <- hash()

  Inbag_Control_Mean_Response_Hash <- hash()
  Inbag_Treatment_Mean_Response_Hash <- hash()
  Inbag_Mean_Response_Hash <- hash()

  OOB_Control_Mean_Response_Hash <- hash()
  OOB_Treatment_Mean_Response_Hash <- hash()
  OOB_Mean_Response_Hash <- hash()

  Inbag_Control_Score_Hash <- hash()
  Inbag_Treatment_Score_Hash <- hash()
  Inbag_Score_Hash <- hash()

  OOB_Control_Score_Hash <- hash()
  OOB_Treatment_Score_Hash <- hash()
  OOB_Score_Hash <- hash()
  
  # Initialize all hashes with the subgroups from the superior_subgroups data.frame
  for( i in 1:NROW( superior_subgroups ) ){


    .set( Inbag_Control_Subgroup_Size_Hash, superior_subgroups$Subgroup[[i]], NULL )
    .set( Inbag_Treatment_Subgroup_Size_Hash, superior_subgroups$Subgroup[[i]], NULL )
    .set( Inbag_Subgroup_Size_Hash, superior_subgroups$Subgroup[[i]], NULL )
    
    .set( OOB_Control_Subgroup_Size_Hash, superior_subgroups$Subgroup[[i]], NULL )
    .set( OOB_Treatment_Subgroup_Size_Hash, superior_subgroups$Subgroup[[i]], NULL )
    .set( OOB_Subgroup_Size_Hash, superior_subgroups$Subgroup[[i]], NULL )
    
    .set( Inbag_Control_Mean_Response_Hash, superior_subgroups$Subgroup[[i]], NULL )
    .set( Inbag_Treatment_Mean_Response_Hash, superior_subgroups$Subgroup[[i]], NULL )
    .set( Inbag_Mean_Response_Hash, superior_subgroups$Subgroup[[i]], NULL )
    
    .set( OOB_Control_Mean_Response_Hash, superior_subgroups$Subgroup[[i]], NULL )
    .set( OOB_Treatment_Mean_Response_Hash, superior_subgroups$Subgroup[[i]], NULL )
    .set( OOB_Mean_Response_Hash, superior_subgroups$Subgroup[[i]], NULL )

    # The score cannot be computed seperately for the two treament arms because
    # the score is frequently a difference between treatment arms
    .set( Inbag_Score_Hash, superior_subgroups$Subgroup[[i]], NULL )
    .set( OOB_Score_Hash, superior_subgroups$Subgroup[[i]], NULL )
    
  }
  
  # Iterate through all bootstrapped samples and append the node size,
  # mean response, and score values to their associated hashes

  for( i in 1:length( TSDT_SAMPLES ) ){

    tree__ <- TSDT_SAMPLES[[i]]@subgroups
    tree__ <- subset( tree__, nchar(trimws(Subgroup)) > 0 )
    

    if( NROW( tree__ ) > 1 ){
    
      for( j in 1:NROW( tree__ ) ){

        # Subset inbag and oob data according to subgroup
        if( tree__$Subgroup[[j]] == 'Overall' ){
          inbag_subgroup__ <- TSDT_SAMPLES[[i]]@inbag
          oob_subgroup__   <- TSDT_SAMPLES[[i]]@oob

        }
        else{
          inbag_subgroup__ <- subgroup( splits = tree__,
                                        node = tree__$NodeID[[j]],
                                        xdata = TSDT_SAMPLES[[i]]@inbag )
          
          
          oob_subgroup__ <- subgroup( splits = tree__,
                                      node = tree__$NodeID[[j]],
                                      xdata = TSDT_SAMPLES[[i]]@oob )
        }
        
        if( !is.null( trt ) ){
          inbag_control_subgroup__ <- subset( inbag_subgroup__, trt == trt_control )
          inbag_treatment_subgroup__ <- subset( inbag_subgroup__, trt != trt_control )
          
          oob_control_subgroup__ <- subset( oob_subgroup__, trt == trt_control )
          oob_treatment_subgroup__ <- subset( oob_subgroup__, trt != trt_control )
        }
        
        # Get generic version of subgroup definition
        generic_subgroup__ <- get_generic_subgroup( tree__$Subgroup[[j]] )

        
        # Populate hashes for each superior subgroup found in superior_subgroups results
        if( has.key( generic_subgroup__, Inbag_Subgroup_Size_Hash ) ){
          
          Inbag_Subgroup_Size_Hash[[ generic_subgroup__ ]] <- c( Inbag_Subgroup_Size_Hash[[ generic_subgroup__ ]], NROW( inbag_subgroup__ ) )
          OOB_Subgroup_Size_Hash[[ generic_subgroup__ ]] <- c( OOB_Subgroup_Size_Hash[[ generic_subgroup__ ]], NROW( oob_subgroup__ ) )

          if( response_type %nin% c( "survival" ) ){
            if( "WEIGHTS__" %nin% names( inbag_subgroup__ ) ) 
                Inbag_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( Inbag_Mean_Response_Hash[[ generic_subgroup__ ]], mean( inbag_subgroup__$response, na.rm = TRUE ) )
            else
                Inbag_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( Inbag_Mean_Response_Hash[[ generic_subgroup__ ]], weighted.mean( inbag_subgroup__$response,
                                                                                                                                        inbag_subgroup__$WEIGHTS__,
                                                                                                                                        na.rm = TRUE ) )
            if( "WEIGHTS__" %nin% names( oob_subgroup__ ) ) 
                OOB_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( OOB_Mean_Response_Hash[[ generic_subgroup__ ]], mean( oob_subgroup__$response, na.rm = TRUE ) )
            else
                OOB_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( OOB_Mean_Response_Hash[[ generic_subgroup__ ]], weighted.mean( oob_subgroup__$response,
                                                                                                                                    oob_subgroup__$WEIGHTS__,
                                                                                                                                    na.rm = TRUE ) )
          }
          
          if( !is.na( tree__$Inbag_Score[[j]] ) )
              Inbag_Score_Hash[[ generic_subgroup__ ]] <- c( Inbag_Score_Hash[[ generic_subgroup__ ]], tree__$Inbag_Score[[j]] )
          
          if( !is.na( tree__$OOB_Score[[j]] ) )
              OOB_Score_Hash[[ generic_subgroup__ ]] <- c( OOB_Score_Hash[[ generic_subgroup__ ]], tree__$OOB_Score[[j]] )


          if( !is.null( trt ) ){

            Inbag_Control_Subgroup_Size_Hash [[ generic_subgroup__ ]] <- c( Inbag_Control_Subgroup_Size_Hash [[ generic_subgroup__ ]], NROW( inbag_control_subgroup__ ) )
            Inbag_Treatment_Subgroup_Size_Hash [[ generic_subgroup__ ]] <- c( Inbag_Treatment_Subgroup_Size_Hash [[ generic_subgroup__ ]], NROW( inbag_treatment_subgroup__ ) )

            OOB_Control_Subgroup_Size_Hash [[ generic_subgroup__ ]] <- c( OOB_Control_Subgroup_Size_Hash [[ generic_subgroup__ ]], NROW( oob_control_subgroup__ ) )
            OOB_Treatment_Subgroup_Size_Hash [[ generic_subgroup__ ]] <- c( OOB_Treatment_Subgroup_Size_Hash [[ generic_subgroup__ ]], NROW( oob_treatment_subgroup__ ) )

            if( response_type %nin% c( "survival" ) ){
              if( "WEIGHTS__" %nin% names( inbag_subgroup__ ) ){
                Inbag_Control_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( Inbag_Control_Mean_Response_Hash[[ generic_subgroup__ ]], mean( inbag_control_subgroup__$response, na.rm = TRUE ) )
                Inbag_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( Inbag_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]], mean( inbag_treatment_subgroup__$response, na.rm = TRUE ) )
              }else{
                Inbag_Control_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( Inbag_Control_Mean_Response_Hash[[ generic_subgroup__ ]], weighted.mean( inbag_control_subgroup__$response,
                                                                                                                                                        inbag_control_subgroup__$WEIGHTS__,
                                                                                                                                                        na.rm = TRUE ) )
                Inbag_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( Inbag_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]], weighted.mean( inbag_treatment_subgroup__$response,
                                                                                                                                                            inbag_treatment_subgroup__$WEIGHTS__,
                                                                                                                                                            na.rm = TRUE ) )
              }
              
              if( "WEIGHTS__" %nin% names( oob_subgroup__ ) ){
                OOB_Control_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( OOB_Control_Mean_Response_Hash[[ generic_subgroup__ ]], mean( oob_control_subgroup__$response, na.rm = TRUE ) )
                OOB_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( OOB_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]], mean( oob_treatment_subgroup__$response, na.rm = TRUE ) )
              }else{
                OOB_Control_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( OOB_Control_Mean_Response_Hash[[ generic_subgroup__ ]], weighted.mean( oob_control_subgroup__$response,
                                                                                                                                                    oob_control_subgroup__$WEIGHTS__,
                                                                                                                                                    na.rm = TRUE ) )
                OOB_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]] <- c( OOB_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]], weighted.mean( oob_treatment_subgroup__$response,
                                                                                                                                                        oob_treatment_subgroup__$WEIGHTS__,
                                                                                                                                                        na.rm = TRUE ) )
              }
            }
          }
        }
        rm( generic_subgroup__ )
      } 
      rm( tree__ )
    }
  }
  
  # Create named list of distributions
  
  DISTRIBUTIONS <- list( Inbag_Subgroup_Size = Inbag_Subgroup_Size_Hash,
                         OOB_Subgroup_Size = OOB_Subgroup_Size_Hash,
                         Inbag_Mean_Response = Inbag_Mean_Response_Hash,
                         OOB_Mean_Response = OOB_Mean_Response_Hash,
                         Inbag_Score = Inbag_Score_Hash,
                         OOB_Score = OOB_Score_Hash )
  
  if( response_type %nin% c( "survival" ) ){
    DISTRIBUTIONS <- c( DISTRIBUTIONS,
                       list( Inbag_Mean_Response = Inbag_Mean_Response_Hash,
                             OOB_Mean_Response = OOB_Mean_Response_Hash ) )
  }
  
  if( !is.null( trt ) ){
    DISTRIBUTIONS <- c( DISTRIBUTIONS,
                       list( Inbag_Control_Subgroup_Size = Inbag_Control_Subgroup_Size_Hash,
                             Inbag_Treatment_Subgroup_Size = Inbag_Treatment_Subgroup_Size_Hash,
                             OOB_Control_Subgroup_Size = OOB_Control_Subgroup_Size_Hash,
                             OOB_Treatment_Subgroup_Size = OOB_Treatment_Subgroup_Size_Hash ) )
    
    if( response_type %nin% c( "survival" ) ){
      DISTRIBUTIONS <- c( DISTRIBUTIONS,
                         list( Inbag_Control_Mean_Response = Inbag_Control_Mean_Response_Hash,
                               Inbag_Treatment_Mean_Response = Inbag_Treatment_Mean_Response_Hash,
                               OOB_Control_Mean_Response = OOB_Control_Mean_Response_Hash,
                               OOB_Treatment_Mean_Response = OOB_Treatment_Mean_Response_Hash ) )
      
    }
  }
  
  # Populate additional columns in superior_subgroups
  
  superior_subgroups$Inbag_Mean_Subgroup_Size <- NA
  superior_subgroups$OOB_Mean_Subgroup_Size <- NA

  if( response_type %nin% c( "survival" ) ){
    superior_subgroups$Inbag_Mean_Response <- NA
    superior_subgroups$OOB_Mean_Response <- NA
  }
  
  superior_subgroups$Inbag_Median_Score <- NA
  superior_subgroups$OOB_Median_Score <- NA


  if( !is.null( trt ) ){

    superior_subgroups$Inbag_Control_Mean_Subgroup_Size <- NA
    superior_subgroups$Inbag_Treatment_Mean_Subgroup_Size <- NA
    superior_subgroups$OOB_Control_Mean_Subgroup_Size <- NA
    superior_subgroups$OOB_Treatment_Mean_Subgroup_Size <- NA

    if( response_type %nin% c( "survival" ) ){
      superior_subgroups$Inbag_Control_Mean_Response <- NA
      superior_subgroups$Inbag_Treatment_Mean_Response <- NA
      superior_subgroups$OOB_Control_Mean_Response <- NA
      superior_subgroups$OOB_Treatment_Mean_Response <- NA
    }
  }
  
  for( i in 1:NROW( superior_subgroups ) ){

    generic_subgroup__ <- as.character( superior_subgroups$Subgroup[[i]] )
    
    if( !is.null( Inbag_Subgroup_Size_Hash[[ generic_subgroup__ ]] ) )
        superior_subgroups$Inbag_Mean_Subgroup_Size[[i]] <- mean( Inbag_Subgroup_Size_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
    if( !is.null( OOB_Subgroup_Size_Hash[[ generic_subgroup__ ]] ) )
        superior_subgroups$OOB_Mean_Subgroup_Size[[i]] <- mean( OOB_Subgroup_Size_Hash[[ generic_subgroup__ ]], na.rm = TRUE )

    if( response_type %nin% c( "survival" ) ){
      if( !is.null( Inbag_Mean_Response_Hash[[ generic_subgroup__ ]] ) )
          superior_subgroups$Inbag_Mean_Response[[i]] <- mean( Inbag_Mean_Response_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
      if( !is.null( OOB_Mean_Response_Hash[[ generic_subgroup__ ]] ) )
          superior_subgroups$OOB_Mean_Response[[i]] <- mean( OOB_Mean_Response_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
    }

    if( !is.null( Inbag_Score_Hash[[ generic_subgroup__ ]] ) )
        superior_subgroups$Inbag_Median_Score[[i]] <- median( Inbag_Score_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
    if( !is.null( OOB_Score_Hash[[ generic_subgroup__ ]] ) )
        superior_subgroups$OOB_Median_Score[[i]] <- median( OOB_Score_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
    
    if( !is.null( trt ) ){
      
      if( !is.null( Inbag_Control_Subgroup_Size_Hash[[ generic_subgroup__ ]] ) )
          superior_subgroups$Inbag_Control_Mean_Subgroup_Size[[i]] <- mean( Inbag_Control_Subgroup_Size_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
      if( !is.null( Inbag_Treatment_Subgroup_Size_Hash[[ generic_subgroup__ ]] ) )
          superior_subgroups$Inbag_Treatment_Mean_Subgroup_Size[[i]] <- mean( Inbag_Treatment_Subgroup_Size_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
      
      if( !is.null( OOB_Control_Subgroup_Size_Hash[[ generic_subgroup__ ]] ) )
          superior_subgroups$OOB_Control_Mean_Subgroup_Size[[i]] <- mean( OOB_Control_Subgroup_Size_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
      if( !is.null( OOB_Treatment_Subgroup_Size_Hash[[ generic_subgroup__ ]] ) )
          superior_subgroups$OOB_Treatment_Mean_Subgroup_Size[[i]] <- mean( OOB_Treatment_Subgroup_Size_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
      
      if( response_type %nin% c( "survival" ) ){
        if( !is.null( Inbag_Control_Mean_Response_Hash[[ generic_subgroup__ ]] ) )
            superior_subgroups$Inbag_Control_Mean_Response[[i]] <- mean( Inbag_Control_Mean_Response_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
        if( !is.null( Inbag_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]] ) )
            superior_subgroups$Inbag_Treatment_Mean_Response[[i]] <- mean( Inbag_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]], na.rm = TRUE )

        if( !is.null( OOB_Control_Mean_Response_Hash[[ generic_subgroup__ ]] ) )
            superior_subgroups$OOB_Control_Mean_Response[[i]] <- mean( OOB_Control_Mean_Response_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
        if( !is.null( OOB_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]] ) )
            superior_subgroups$OOB_Treatment_Mean_Response[[i]] <- mean( OOB_Treatment_Mean_Response_Hash[[ generic_subgroup__ ]], na.rm = TRUE )
      }
    }
    
    rm( generic_subgroup__ )
  }
  
  superior_subgroups$score <- superior_subgroups$External_Consistency * superior_subgroups$Internal_Consistency
        
  ###############################################################################
  # Compute permutation-adjusted p-value for each subgroup in superior_subgroups#
  # data.frame.                                                                 #
  ###############################################################################

  if( n_permutations > 0 ){
       
    if( n_cpu == 1 ){
      NULL_SCORES <- get_null_scores( n_permutations = n_permutations,
                                      response = response,
                                      trt = trt,
                                      covariates = source_data,
                                      tree_builder = tree_builder,
                                      tree_builder_parameters = tree_builder_parameters,
                                      permute_method = permute_method,
                                      permute_arm = permute_arm,
                                      n_samples = n_samples,
                                      sampling_method = sampling_method,
                                      inbag_proportion = inbag_proportion,
                                      desirable_response = desirable_response,
                                      scoring_function = scoring_function,
                                      scoring_function_parameters = scoring_function_parameters,
                                      inbag_score_margin = inbag_score_margin,
                                      oob_score_margin = oob_score_margin,
                                      eps = eps,
                                      min_subgroup_n_control = min_subgroup_n_control,
                                      min_subgroup_n_trt = min_subgroup_n_trt,
                                      min_subgroup_n_oob_control = min_subgroup_n_oob_control,
                                      min_subgroup_n_oob_trt = min_subgroup_n_oob_trt,
                                      n_cpu = n_cpu,
                                      trace = trace )
    }


    else if( n_cpu > 1 ){# BEGIN parallel computing

      for( p in 1:n_cpu ){
      
        parallel::mcparallel({
          
          NULL_SCORES <- get_null_scores( n_permutations = ceiling(n_permutations/n_cpu),
                                          response = response,
                                          trt = trt,
                                          covariates = source_data,
                                          tree_builder = tree_builder,
                                          tree_builder_parameters = tree_builder_parameters,
                                          permute_method = permute_method,
                                          permute_arm = permute_arm,
                                          n_samples = n_samples,
                                          sampling_method = sampling_method,
                                          inbag_proportion = inbag_proportion,
                                          desirable_response = desirable_response,
                                          scoring_function = scoring_function,
                                          scoring_function_parameters = scoring_function_parameters,
                                          inbag_score_margin = inbag_score_margin,
                                          oob_score_margin = oob_score_margin,
                                          eps = eps,
                                          min_subgroup_n_control = min_subgroup_n_control,
                                          min_subgroup_n_trt = min_subgroup_n_trt,
                                          min_subgroup_n_oob_control = min_subgroup_n_oob_control,
                                          min_subgroup_n_oob_trt = min_subgroup_n_oob_trt,
                                          n_cpu = n_cpu,
                                          trace = trace )
        })
      }
      
      JOBS <- parallel::mccollect()
      
      NULL_SCORES <- NULL
    
      for( job in JOBS )
          NULL_SCORES <- c( NULL_SCORES, job )
      
    }# END parallel computing
    
    if( NROW(superior_subgroups) == 0 ){# BEGIN no superior subgroups found
      superior_subgroups$Adjusted_Pvalue <- NA
      superior_subgroups$Strength <- ""
    }# END no superior subgroups found
    
    else{ # BEGIN compute adjusted p-value for each superior subgroup
      superior_subgroups$Adjusted_Pvalue <- 1
      superior_subgroups$Strength <- "Not confirmed"
      
      for( i in 1:NROW(superior_subgroups) ){
 
        # For computation of Adjusted_Pvalue:
        #  numerator = number of NULL_SCORES > superior_subgroups$score + 1
        #  denominator = number of NULL_SCORES + 1
        # Then Adjusted_Pvalue = numerator/denominator)
        
        numerator__ <- 1
        if( any( NULL_SCORES > superior_subgroups$score[[i]] ) )
            numerator__ <- table( NULL_SCORES > superior_subgroups$score[[i]] )[["TRUE"]]  + 1

        ## Account for corner case where External Consistency == 0. P-value
        ## should always be 1.
        if( superior_subgroups$score[[i]] == 0 )
            numerator__ <- length( NULL_SCORES ) + 1
        
        denominator__ <- length( NULL_SCORES ) + 1
        
        superior_subgroups$Adjusted_Pvalue[[i]] <- numerator__/denominator__
        
        rm( numerator__ )
        rm( denominator__ )
  
        superior_subgroups$Strength[i] <- ifelse( superior_subgroups$Adjusted_Pvalue[i] < strength_cutpoints[3], "Weak", superior_subgroups$Strength[i] )
        superior_subgroups$Strength[i] <- ifelse( superior_subgroups$Adjusted_Pvalue[i] < strength_cutpoints[2], "Moderate", superior_subgroups$Strength[i] )
        superior_subgroups$Strength[i] <- ifelse( superior_subgroups$Adjusted_Pvalue[i] < strength_cutpoints[1], "Strong", superior_subgroups$Strength[i] )

      }# END compute adjusted p-value for each superior subgroup
    }# END at least one superior subgroup found
  }# END n_permutations > 0 (compute p-pvalue)

  
  superior_subgroups$Suggested_Cutoff <- ""
  
  for( i in 1:NROW( superior_subgroups ) ){ # BEGIN get suggested cutoff

    if( superior_subgroups$Subgroup[[i]] %nin% c('Overall') ){
      cutpoints__ <-  get_cutpoints( CutpointDistribution, as.character( superior_subgroups$Subgroup[[i]] ) )
  
      SUBSUBS <- strsplit( as.character( superior_subgroups$Subgroup[[i]] ), split = " & " )[[1]]
  
      cutoff <- "["
    
      for( subsub in SUBSUBS ){
        
        if( grepl( pattern = "%in%c(", subsub, fixed = TRUE ) ){
          start_pos <- regexpr( pattern = "%in%c(", subsub, fixed = TRUE ) + 5
          stop_pos <- nchar( subsub )
          
          cutoff <- paste0( cutoff, substr( subsub, start = start_pos, stop = stop_pos ),";" )
        }else if( grepl( pattern = "%nin%c(", subsub, fixed = TRUE ) ){
          start_pos <- regexpr( pattern = "%nin%c(", subsub, fixed = TRUE ) + 6
          stop_pos <- nchar( subsub )
          
          cutoff <- paste0( cutoff, substr( subsub, start = start_pos, stop = stop_pos ),";" )
        }else{
          
          cutpoint_dist <- get_cutpoints( CutpointDistribution,
                                          as.character( superior_subgroups$Subgroup[[i]] ),
                                          subsub )
          
          cutoff <- paste0( cutoff, median( cutpoint_dist ),";" )
        }
        
      }
      
      # Remove trailing comma and add closing right bracket
      cutoff <- substr( cutoff, start = 1, stop = nchar(cutoff) - 1 )
      cutoff <- paste0( cutoff, "]" )
      
      superior_subgroups$Suggested_Cutoff[[i]] <- cutoff
    }
  } # END get suggested cutoff

  # Sort by p-value, descending score, internal and external consistencies,
  # and subgroup. Overall subgroup should always be first.
  superior_subgroups$is_overall <- ifelse( superior_subgroups$Subgroup == 'Overall', 1, 0 )
  
  if( n_permutations > 0 )
  
      superior_subgroups <- superior_subgroups[ order( -superior_subgroups$is_overall,
                                                       superior_subgroups$Adjusted_Pvalue,
                                                       -superior_subgroups$score,
                                                       -superior_subgroups$Internal_Consistency,
                                                       -superior_subgroups$External_Consistency,
                                                       superior_subgroups$Subgroup ), ]
  else
      superior_subgroups <- superior_subgroups[ order( -superior_subgroups$is_overall,
                                                       -superior_subgroups$score,
                                                       -superior_subgroups$Internal_Consistency,
                                                       -superior_subgroups$External_Consistency,
                                                       superior_subgroups$Subgroup ), ]

      
  superior_subgroups$score <- NULL
  superior_subgroups$is_overall <- NULL
  row.names( superior_subgroups ) <- NULL
  
  # Order columns in superior_subgroups      
  TWO_ARM_COLS <-  c( "Subgroup",
                     "Internal_Consistency",
                     "External_Consistency",
                     
                     "Inbag_Mean_Subgroup_Size",
                     "Inbag_Control_Mean_Subgroup_Size",
                     "Inbag_Treatment_Mean_Subgroup_Size",
                     
                     "OOB_Mean_Subgroup_Size",
                     "OOB_Control_Mean_Subgroup_Size",
                     "OOB_Treatment_Mean_Subgroup_Size",
                     
                     "Inbag_Mean_Response",
                     "Inbag_Control_Mean_Response",
                     "Inbag_Treatment_Mean_Response",
                     
                     "OOB_Mean_Response",
                     "OOB_Control_Mean_Response",
                     "OOB_Treatment_Mean_Response" )

  if( response_type == "survival" ){

    # Remove mean response columns b/c they aren't relevant for survival outcome
    TWO_ARM_COLS <- TWO_ARM_COLS[ !grepl( pattern = "Response", TWO_ARM_COLS ) ]
      
    SURVIVAL_COLS <- c( "CoxPH_Median_Survival_OBSERVED",
                        "CoxPH_Control_Median_Survival_OBSERVED",
                        "CoxPH_Treatment_Median_Survival_OBSERVED",
                        "CoxPH_Hazard_Ratio_OBSERVED" )
    TWO_ARM_COLS <- c( TWO_ARM_COLS, SURVIVAL_COLS )
    superior_subgroups[,SURVIVAL_COLS] <- NA
  }

  if( n_permutations > 0 )
      TWO_ARM_COLS <- c( TWO_ARM_COLS, "Adjusted_Pvalue","Strength" )

  TWO_ARM_COLS <- c( TWO_ARM_COLS, "Inbag_Median_Score", "OOB_Median_Score","Suggested_Cutoff" )
      
  ONE_ARM_COLS <- TWO_ARM_COLS[ !grepl( pattern = "Control|Treatment", TWO_ARM_COLS ) ]
  
  if( is.null( trt ) )
      superior_subgroups <- superior_subgroups[,ONE_ARM_COLS]
  
  else
      superior_subgroups <- superior_subgroups[,TWO_ARM_COLS]

  if(response_type == "survival" ){

    superior_subgroups__ <- superior_subgroups # create temporary copy of superior_subgroups
    superior_subgroups__$Subgroup <- as.character( superior_subgroups__$Subgroup )
    superior_subgroups__$NodeID <- 1:NROW( superior_subgroups ) # add temp NodeID 

    df__ <- as.data.frame( response )
    names( df__ ) <- "response"
    df__$trt <- trt
    df__ <- cbind( df__, covariates )
    
    # Construct formula for Cox PH model
    formula_y <- 'response ~ '
    one_arm_formula <- as.formula( paste0( formula_y, '1' ) )
    two_arm_formula <- as.formula( paste0( formula_y, 'trt' ) )
    
    # Subset data for each suggested subgroup in superior subgroups, fit the
    # Cox PH model, and compute median survival, and hazard ratio
    for( i in 1:NROW( superior_subgroups__ ) ){
      
      if( superior_subgroups__$Subgroup[[i]] == 'Overall' )
          sg__ <- df__
      
      else{
        # Replace anonymized subgroup with suggested subgroup
        superior_subgroups__$Subgroup[[i]] <- get_suggested_subgroup( anonymized_subgroup = superior_subgroups__$Subgroup[[i]], 
                                                               suggested_cutoff = superior_subgroups__$Suggested_Cutoff[[i]] )

        # Get subset of original data associated with suggested subgroup
        sg__ <- subgroup( superior_subgroups__, node = superior_subgroups__$NodeID[[i]], xdata = df__ )
      }

      if( !is.null( trt ) ){

        two_arm_coxph__ <- coxph( two_arm_formula, data = sg__ )

        control_coxph__ <- coxph( one_arm_formula, data = subset( sg__, trt == trt_control ) )
        treatment_coxph__ <- coxph( one_arm_formula, data = subset( sg__, trt != trt_control ) )
      
        # Update values in original superior_subgroups
        superior_subgroups$CoxPH_Hazard_Ratio_OBSERVED[[i]] <- exp( coef( two_arm_coxph__ ) )

        superior_subgroups$CoxPH_Median_Survival_OBSERVED[[i]] <- quantile( survfit( two_arm_coxph__ ), probs = 0.50 )$quantile[[1]]
        superior_subgroups$CoxPH_Control_Median_Survival_OBSERVED[[i]] <- quantile( survfit( control_coxph__ ), probs = 0.50 )$quantile[[1]]
        superior_subgroups$CoxPH_Treatment_Median_Survival_OBSERVED[[i]] <- quantile( survfit( treatment_coxph__ ), probs = 0.50 )$quantile[[1]]

        rm( two_arm_coxph__ )
        rm( control_coxph__ )
        rm( treatment_coxph__ )
      }
      else{

        one_arm_coxph__ <- coxph( one_arm_formula, data = sg__ )
      
        # Update values in original superior_subgroups
        superior_subgroups$CoxPH_Hazard_Ratio_OBSERVED <- NULL

        superior_subgroups$CoxPH_Median_Survival_OBSERVED[[i]] <- quantile( survfit( one_arm_coxph__ ), probs = 0.50 )$quantile[[1]]
        superior_subgroups$CoxPH_Control_Median_Survival_OBSERVED <- NULL
        superior_subgroups$CoxPH_Treatment_Median_Survival_OBSERVED <- NULL

        rm( one_arm_coxph__ )
        
      }

      rm( sg__ ) 
    }
    
    rm( df__ )
    rm( superior_subgroups__ )
  }
      
  # Remove p-value and strength for Overall subgroup
  for( i in 1:NROW( superior_subgroups ) ){
    if( superior_subgroups$Subgroup[[i]] == 'Overall' ){
      superior_subgroups$Internal_Consistency[[i]] <- NA
      superior_subgroups$External_Consistency[[i]] <- NA
      superior_subgroups$Adjusted_Pvalue[[i]] <- NA
      superior_subgroups$Strength[[i]] <- 'N/A'
    }
  }

  ## Make superior_subgroups$Subgroup character rather than factor
  superior_subgroups$Subgroup <- as.character( superior_subgroups$Subgroup )
  
  # Add NULL_SCORES to DISTRIBUTIONS list
  if( exists( "NULL_SCORES" ) )
      DISTRIBUTIONS$Null_Scores <- hash( 'Null_Scores', NULL_SCORES )

  ###############################################################################
  # Convert translated logical inequalities back to their original values       #
  ###############################################################################
  decode_sample_inequalities <- function( pattern, replacement ){
  
    ## Convert values in TSDT_SAMPLES
    for( s in 1:n_samples ){
      
      ## Convert %%__GE__%% back to >=
      for( j in 1:NCOL( TSDT_SAMPLES[[s]]@inbag ) ){
        if( is( TSDT_SAMPLES[[s]]@inbag[,j], 'character' ) ){
          TSDT_SAMPLES[[s]]@inbag[,j] <- gsub( pattern = pattern,
                                               replacement = replacement,
                                               fixed = TRUE,
                                               TSDT_SAMPLES[[s]]@inbag[,j] )
        }
      }
      rm( j )
      
      for( j in 1:NCOL( TSDT_SAMPLES[[s]]@oob ) ){
        if( is( TSDT_SAMPLES[[s]]@oob[,j], 'character' ) ){
          TSDT_SAMPLES[[s]]@oob[,j] <- gsub( pattern = pattern,
                                             replacement = replacement,
                                             fixed = TRUE,
                                             TSDT_SAMPLES[[s]]@oob[,j] )
        }
      }
      rm( j )
      
      for( j in 1:NCOL( TSDT_SAMPLES[[s]]@subgroups ) ){
        if( is( TSDT_SAMPLES[[s]]@subgroups[,j], 'character' ) ){
          TSDT_SAMPLES[[s]]@subgroups[,j] <- gsub( pattern = pattern,
                                                   replacement = replacement,
                                                   fixed = TRUE,
                                                   TSDT_SAMPLES[[s]]@subgroups[,j] )
        }
      }
      rm( j )
      
    }## END s in 1:n_samples
  }## END decode_sample_inqualities
  
  decode_sample_inequalities( pattern = '%%__GE__%%', replacement = '>=' )
  decode_sample_inequalities( pattern = '%%__LE__%%', replacement = '<=' )
  decode_sample_inequalities( pattern = '%%__GT__%%', replacement = '>' )
  decode_sample_inequalities( pattern = '%%__LE__%%', replacement = '<' )
  decode_sample_inequalities( pattern = '%%__EQ__%%', replacement = '=' )
  decode_sample_inequalities( pattern = '%%__EMPTY_STRING__%%', replacement = '' )
  
  ## Convert values in DISTRIBUTIONS
  decode_distribution_inequalities <- function( pattern, replacement ){
    
    for( i in 1:length( DISTRIBUTIONS ) ){
      
      keys__ <- names( DISTRIBUTIONS[[i]]@.xData )
      
      INDEX <- which( grepl( pattern = pattern, keys__ ) )
      for( index in INDEX ){
        old__ <- keys__[index]
        new__ <- gsub( pattern = pattern, replacement = replacement, fixed = TRUE, old__ )
        DISTRIBUTIONS[[i]]@.xData[[ new__ ]] <- DISTRIBUTIONS[[i]]@.xData[[ old__ ]]
        eval( parse( text = paste0( 'rm( "', old__,'", envir = DISTRIBUTIONS[[i]]@.xData )' ) ) )
        rm( old__, new__ )
      }
      rm( INDEX, keys__ )
    }
    rm( i )
  }
  
  decode_distribution_inequalities( pattern = '%%__GE__%%', replacement = '>=' )
  decode_distribution_inequalities( pattern = '%%__LE__%%', replacement = '<=' )
  decode_distribution_inequalities( pattern = '%%__GT__%%', replacement = '>' )
  decode_distribution_inequalities( pattern = '%%__LE__%%', replacement = '<' )
  decode_distribution_inequalities( pattern = '%%__EQ__%%', replacement = '=' )
  decode_distribution_inequalities( pattern = '%%__EMPTY_STRING__%%', replacement = '' )

  ## Convert values in superior_subgroups
  superior_subgroups$Subgroup <- gsub( pattern = '%%__GE__%%', replacement = '>=', fixed = TRUE, superior_subgroups$Subgroup )
  superior_subgroups$Subgroup <- gsub( pattern = '%%__LE__%%', replacement = '<=', fixed = TRUE, superior_subgroups$Subgroup )
  superior_subgroups$Subgroup <- gsub( pattern = '%%__GT__%%', replacement = '>', fixed = TRUE, superior_subgroups$Subgroup )
  superior_subgroups$Subgroup <- gsub( pattern = '%%__LT__%%', replacement = '<', fixed = TRUE, superior_subgroups$Subgroup )
  superior_subgroups$Subgroup <- gsub( pattern = '%%__EQ__%%', replacement = '=', fixed = TRUE, superior_subgroups$Subgroup )
  superior_subgroups$Subgroup <- gsub( pattern = '%%__EMPTY_STRING__%%', replacement = '', fixed = TRUE, superior_subgroups$Subgroup )

  
  superior_subgroups$Suggested_Cutoff <- gsub( pattern = '%%__GE__%%', replacement = '>=', fixed = TRUE, superior_subgroups$Suggested_Cutoff )
  superior_subgroups$Suggested_Cutoff <- gsub( pattern = '%%__LE__%%', replacement = '<=', fixed = TRUE, superior_subgroups$Suggested_Cutoff )
  superior_subgroups$Suggested_Cutoff <- gsub( pattern = '%%__GT__%%', replacement = '>', fixed = TRUE, superior_subgroups$Suggested_Cutoff )
  superior_subgroups$Suggested_Cutoff <- gsub( pattern = '%%__LT__%%', replacement = '<', fixed = TRUE, superior_subgroups$Suggested_Cutoff )
  superior_subgroups$Suggested_Cutoff <- gsub( pattern = '%%__EQ__%%', replacement = '=', fixed = TRUE, superior_subgroups$Suggested_Cutoff )
  superior_subgroups$Suggested_Cutoff <- gsub( pattern = '%%__EMPTY_STRING__%%', replacement = '', fixed = TRUE, superior_subgroups$Suggested_Cutoff )
  
  ###############################################################################
  # Instantiate an object of class TSDT with the vector of TSDT_Sample objects, #
  # the data.frame containing the superior_subgroups values, and a list of the  #
  # parameters used to create this data.                                        #
  ###############################################################################
  
  TSDT1__ <- new( "TSDT",
                parameters = list( response_type = response_type,
                                   percentile = percentile,
                                   min_subgroup_n_control = min_subgroup_n_control,
                                   min_subgroup_n_trt = min_subgroup_n_trt,
                                   min_subgroup_n_oob_control = min_subgroup_n_oob_control,
                                   min_subgroup_n_oob_trt = min_subgroup_n_oob_trt,
                                   maxdepth = maxdepth,
                                   rootcompete = rootcompete,
                                   competedepth = competedepth,
                                   trt_control = trt_control,
                                   permute_method = as.character( substitute( permute_method ) ),
                                   permute_arm = as.character( substitute( permute_arm ) ),
                                   n_samples = n_samples,
                                   sampling_method = sampling_method,
                                   inbag_proportion = inbag_proportion,
                                   inbag_score_margin = inbag_score_margin,
                                   oob_score_margin = oob_score_margin,
                                   strength_cutpoints = strength_cutpoints,
                                   n_permutations = n_permutations,
                                   n_cpu = n_cpu,
                                   tree_builder = tree_builder,
                                   tree_builder_parameters = tree_builder_parameters,
                                   desirable_response = desirable_response,
                                   scoring_function = scoring_function_name,
                                   scoring_function_parameters = scoring_function_parameters,
                                   eps = eps ),
               
                samples = TSDT_SAMPLES,
                superior_subgroups = superior_subgroups,
                cutpoints = CutpointDistribution,
                distributions = DISTRIBUTIONS )
  
  return( TSDT1__ )
}
        
## END OF FILE

