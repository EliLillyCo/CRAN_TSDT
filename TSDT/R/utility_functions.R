#!usr/bin/dev/R
#################################################################################
# FILENAME   : utility_functions.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 08/09/2012
# DESCRIPTION: Collection of helpful utility functions.
#################################################################################

#' @title %nin%
#' @description Negation of %in%
#' @param a,b two data objects to test for non-equality
#' @return Logical value TRUE or FALSE
#' @examples
#' %nin%(1,2)
#' %nin%(1,1)
#' @export
`%nin%` <- function(a, b){
  return( !(a %in% b) )
}

#' @title binary transform
#' @description Converts any variable with two possible values to a \{0,1\} binary
#' variable.
#' @param x A variable with two possible values.
#' @return A vector with values in \{0,1\}.
#' @examples
#' ## Convert a variable that takes values 'A' and 'B' to 0 and 1
#' x <- sample( c('A','B'), size = 10, prob = c(0.5,0.5), replace = TRUE )
#' print(x);flush.console()
#' binary_transform( x )
#' @export
binary_transform <- function( x ){

  # If x is already a \{0,1\} variable do nothing
  if( any( x %nin% c(0,1) ) ){
    
    xnum <- as.numeric( as.factor(x) ) - 1
    
    # Confirm that xnum has only two unique values
    if( any( xnum %nin% c(0,1) ) )
        stop( "ERROR: binary variable must have at most two levels" )
    
    table__ <- table( x, xnum )
    
    cat( paste0( "NOTE: binary mappings are as follows:\t",
                row.names(table__)[1], "->", colnames(table__)[1],
                ", ", row.names(table__)[2], "->", colnames(table__)[2], "\n" ) )
    x <- xnum
  }
  
  return( x )
}

#' @title partition
#' @description Partitions a vector x into n groups of roughly equal size.
#' @param x Vector to partition.
#' @param n Number of (roughly) equally-sized groups
#' @return A list of partitions of the vector x.
#' @examples
#' x <- 1:10
#' partition( x, 3 )
#' @export
partition <- function( x, n ){

  if( floor(n) != n || n < 1 )
    stop( "ERROR: n must be a positive integer" )
  
  if( n > length(x) )
    stop( "ERROR: n must be <= length(x)" )
  
  partitions__ <- split(x, as.integer((seq_along(x) - 1) / (length(x)/n) ) )

  names( partitions__ ) <- NULL
  
  return( partitions__ )
}


#' @title folds
#' @description Partition data into k folds for k-fold cross-validation. Adds a variable
#' fold_id to the data.frame.
#' @param x data.frame to partition into k folds for k-fold cross-validation.
#' @param k Number of folds to use in cross-validation
#' @return A list of partitions of the vector x.
#' @examples
#' # Generate random example data
#' N <- 200
#' ID <- 1:N
#' continuous_response = runif( min = 0, max = 20, n = N )
#' X1 <- runif( N, min = 0, max = 1 )
#' X2 <- runif( N, min = 0, max = 1 )
#' X3 <- sample( c(0,1), size = N, prob = c(0.2,0.8), replace = TRUE )
#' X4 <- sample( c('A','B','C'), size = N, prob = c(0.6,0.3,0.1), replace = TRUE )
#'
#' df <- data.frame( ID )
#' names( df ) <- "ID"
#' df$response <- continuous_response
#' df$X1 <- X1
#' df$X2 <- X2
#' df$X3 <- factor( X3 )
#' df$X4 <- factor( X4 )
#'
#' ## Partition data into 5 folds
#' ex1 <- folds( df, k = 5 )
#'
#' ## Partition data into 10 folds
#' ex2 <- folds( df, k = 10 )
#' @export
folds <- function( x, k ){

  # create row_id so rows can be returned in original order
  x$row_id <- 1:NROW(x)

  # scramble row id
  rows <- sample( x$row_id, size = NROW(x), replace = FALSE )
  
  # partition row_id into k parts
  fold_id <- partition( rows, k )
  
  # shuffle rows of x according to scrambled row_id
  x <- x[rows,]
  
  # create fold_id
  for( i in 1:length(fold_id) )
      fold_id[[i]] <- rep( i, length( fold_id[[i]] ) )
  
  fold_id <- unlist( fold_id )
  
  x$fold_id <- fold_id
  
  # reorder rows to original order
  x <- x[order(x$row_id),]
  x$row_id <- NULL
  
  return( x ) 
}


#' @title unpack_args
#' @description Assign the elements of a named list in current environment.
#' @details This function takes a list of named entities and assigns each element
#' of the list to its name in the calling environment.
#' @param args List of entities to be assigned.
#' @seealso \link{assign}, \link{parent.frame}
#' @examples
#' ## Create a list of named elements
#' arglist <- list( one = 1, two = 2, color = "blue" )
#' 
#' ## The variables one, two, and color do not exist in the current environment
#' ls()
#' 
#' ## Unpack the elements in arglist
#' unpack_args( arglist )
#' 
#' ## Now the variables one, two, and color do exist in the current environment
#' ls()
#' one
#' @export
unpack_args <- function( args ){
  for( i in 1:length(args) ){
    #print( names(args)[i] )
    #print( args[[i]] )
    #flush.console()
    assign( names(args)[i], args[[i]], pos = parent.frame(1) )
  }
}


# Return vector of column indices of factor columns
factor_cols <- function( df ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  x <- NULL
  
  return( which( sapply( x[,1:NCOL(x)], is.factor ) ) )
}

#' @title function_parameter_names
#' @description Returns a character vector of the specified function's parameters
#' @param FUN The name of a function
#' @return A character vector of function parameter names
#' @examples
#' ## Define a function
#' example_function <- function( parm1, arg2, x, bool = FALSE ){
#'   cat( "This is an example function.\n" )
#' }
#' 
#' ## Return the function parameter names
#' function_parameter_names( example_function )
#' @export
function_parameter_names <- function( FUN ){
    
  raw_func_str <- capture.output( FUN )

  # Strip comments
  concat_func_str <- NULL
  for( i in 1:NROW( raw_func_str ) ){
    comment_start <- regexpr( pattern = "#", raw_func_str[[i]], useBytes = TRUE )[[1]]
    if( comment_start > 0 )
      raw_func_str[[i]] <- substr( raw_func_str[[i]], 1, comment_start -1 )
    concat_func_str <- paste0( concat_func_str, raw_func_str[[i]] )
  }
  
  concat_func_str <- gsub( pattern = " ", replacement = "", concat_func_str )
  args_start <- regexpr( pattern = "\\(", concat_func_str, useBytes = TRUE ) + 1
  args_stop <- regexpr( pattern = "(\\)\\{)|(\\)UseMethod)", concat_func_str, useBytes = TRUE ) - 1
  params <- strsplit( substr( concat_func_str, args_start, args_stop ), split = "," )[[1]]
  
  if( length( params ) > 0 ){
    for( i in 1:length(params) ){
    
      equal_sign <- regexpr( pattern = "=", params[[i]], useBytes = TRUE )
      if( equal_sign > 0 )
        params[[i]] <- substr( params[[i]], 1, equal_sign - 1 )
    }
  }


  else if( length( params ) == 0 ){
    params <- NULL
  }
  
  return( params )

}


#' @title na2empty
#' @description Replace all instances of NA in character variable with empty
#' string.
#' @param x A character vector.
#' @seealso \link{unfactor}
#' @return A character vector with NA values replaced with empty string.
#' @examples
#' ## Create character variable with missing values
#' ex1 <- c( 'A', NA, 'B', NA, 'C', NA )
#' ex1
#' 
#' ## Replace NAs with empty string
#' ex1 <- na2empty( ex1 )
#' ex1
#' @export
na2empty <- function( x ){

  if( !is.character( x ) )
      stop( "ERROR: x should be a character vector" )
  
  return( ifelse( !is.na( x ), x, "" ) )
  
}

complement_logic <- function( logic ){

  logic <- gsub( pattern = " ", replacement = "", logic )
  
  if( logic == "<" )
      return( ">=" )

  if( logic == "<=" )
      return( ">" )

  if( logic == ">" )
      return( "<=" )
  
  if( logic == ">=" )
      return( "<" )
  
  if( logic == "%in%" )
      return( "%nin%" )

  if( logic == "%nin%" )
      return( "%in%" )

  stop( "ERROR: logic value not recognized" )
  
}

parse_subgroup <- function( subgroup, logic ){
  
  if( logic %nin% c('<','<=','>','>=','%in%','%nin%') )
      stop( "ERROR: logic must be one of {<,<=,>,>=,%in%,%nin%}" )
  
  SplitVariable <- strsplit( subgroup, split = logic )[[1]][1]
  SplitValue    <- strsplit( subgroup, split = logic )[[1]][2]

  if( logic %in% c('<','<=','>','>=') ){
    SplitValue <- as.numeric( SplitValue )
  }

  else if( logic %in% c('%in%','%nin%') ){
    SplitValue <- gsub( pattern = "c(", replacement = "", SplitValue, fixed = TRUE )
    SplitValue <- gsub( pattern = "(", replacement = "", SplitValue, fixed = TRUE )
    SplitValue <- gsub( pattern = ")", replacement = "", SplitValue, fixed = TRUE )
  }

  return( list( SplitVariable = SplitVariable, SplitValue = SplitValue ) )
}


# If a subgroup has multiple terms that use the same split variable and
# same split logic multiple times then collapse as appropriate.
# Example:  X1<0.5 & X<0.6 should be collapsed to X1<0.5
collapse_redundant_splits <- function( subgroup ){

  requireNamespace( "hash", quietly = TRUE )
  
  LT_HASH <- hash()     # splits with <
  LE_HASH <- hash()     # splits with <=
  GT_HASH <- hash()     # splits with >
  GE_HASH <- hash()     # splits with >=
  NOMINAL_IN_HASH  <- hash() # splits with nominal values
  NOMINAL_NIN_HASH  <- hash() 

  SPLIT_VARIABLES <- NULL
  
  # Populate hashes
  for( sg in subgroup ){
    
    # Populate LE_HASH if split contains <=
    if( grepl( pattern = "<=", sg ) ){

      ## Create NULL placeholders to prevent NOTE in R CMD check
      SplitVariable <- NULL
      SplitValue <- NULL
      
      unpack_args( parse_subgroup( subgroup = sg, logic = '<=' ) )
      
      if( has.key( SplitVariable, LE_HASH ) )
          LE_HASH[[ SplitVariable ]] <- c( LE_HASH[[ SplitVariable ]], SplitValue )
      
      else
          LE_HASH[[ SplitVariable ]] <- SplitValue

      SPLIT_VARIABLES <- c( SPLIT_VARIABLES, SplitVariable )
      rm( SplitVariable )
      rm( SplitValue )
    }
    # Populate GE_HASH if split contains >=
    else if( grepl( pattern = ">=", sg ) ){

      unpack_args( parse_subgroup( subgroup = sg, logic = '>=' ) )
      
      if( has.key( SplitVariable, GE_HASH ) )
          GE_HASH[[ SplitVariable ]] <- c( GE_HASH[[ SplitVariable ]], SplitValue )
      
      else
          GE_HASH[[ SplitVariable ]] <- SplitValue

      SPLIT_VARIABLES <- c( SPLIT_VARIABLES, SplitVariable )
      rm( SplitVariable )
      rm( SplitValue )
    }
    # Populate LT_HASH if split contains <
    else if( grepl( pattern = "<", sg ) ){
      
      unpack_args( parse_subgroup( subgroup = sg, logic = '<' ) )
      
      if( has.key( SplitVariable, LT_HASH ) )
          LT_HASH[[ SplitVariable ]] <- c( LT_HASH[[ SplitVariable ]], SplitValue )
      
      else
          LT_HASH[[ SplitVariable ]] <- SplitValue
      
      SPLIT_VARIABLES <- c( SPLIT_VARIABLES, SplitVariable )
      rm( SplitVariable )
      rm( SplitValue )
    }
    # Populate GT_HASH if split contains >
    else if( grepl( pattern = ">", sg ) ){

      unpack_args( parse_subgroup( subgroup = sg, logic = '>' ) )
      
      if( has.key( SplitVariable, GT_HASH ) )
          GT_HASH[[ SplitVariable ]] <- c( GT_HASH[[ SplitVariable ]], SplitValue )
      
      else
          GT_HASH[[ SplitVariable ]] <- SplitValue
      
      SPLIT_VARIABLES <- c( SPLIT_VARIABLES, SplitVariable )
      rm( SplitVariable )
      rm( SplitValue )
    }
    # Populate NOMINAL_IN_HASH if split contains %in%
    else if( grepl( pattern = "%in%", sg ) ){
      
      unpack_args( parse_subgroup( subgroup = sg, logic = '%in%' ) )
      
      if( has.key( SplitVariable, NOMINAL_IN_HASH ) )
          NOMINAL_IN_HASH[[ SplitVariable ]] <- c( NOMINAL_IN_HASH[[ SplitVariable ]], SplitValue )
      
      else
          NOMINAL_IN_HASH[[ SplitVariable ]] <- SplitValue
      
      SPLIT_VARIABLES <- c( SPLIT_VARIABLES, SplitVariable )
      rm( SplitVariable )
      rm( SplitValue )
    }
    # Populate NOMINAL_NIN_HASH if split contains %nin%
    else if( grepl( pattern = "%nin%", sg ) ){
          
      unpack_args( parse_subgroup( subgroup = sg, logic = '%nin%' ) )
      
      if( has.key( SplitVariable, NOMINAL_NIN_HASH ) )
          NOMINAL_NIN_HASH[[ SplitVariable ]] <- c( NOMINAL_NIN_HASH[[ SplitVariable ]], SplitValue )
      
      else
          NOMINAL_NIN_HASH[[ SplitVariable ]] <- SplitValue
      
      SPLIT_VARIABLES <- c( SPLIT_VARIABLES, SplitVariable )
      rm( SplitVariable )
      rm( SplitValue )
    }
    
    else
        stop( paste0( "ERROR: valid logic operator not found in: ", sg,
                      ". logic must be one of {<,<=,>,>=,%in%,%nin%}" ) )
    
  }

  SPLIT_VARIABLES <- unique( SPLIT_VARIABLES )
  
  # Construct a non-redundant subgroup (NRS)
  NRS <- NULL

  for( k in SPLIT_VARIABLES ){

    # If a split variable has both a < and <= logic then use < as the split logic
    # if the minimum split value in LT_HASH is less than or equal to the minimum
    # split value in the LE_HASH. Otherwise, if the LE_HASH has the minimum split
    # value use <= as the split logic.
    if( k %in% keys( LT_HASH ) && k %in% keys( LE_HASH ) ){
      if( min( LT_HASH[[ k ]] ) <= min( LE_HASH[[ k ]] ) )
          NRS <- c( NRS,  paste0( k, "<", min( LT_HASH[[ k ]] ) ) )
      else
          NRS <- c( NRS,  paste0( k, "<=", min( LE_HASH[[ k ]] ) ) )
    }# end both < and <=
       
    # Split variable has < logic but not <= logic
    if( k %in% keys( LT_HASH ) && k %nin% keys( LE_HASH ) )
        NRS <- c( NRS, paste0( k, "<", min( LT_HASH[[ k ]] ) ) )
      
    # Split variable has <= logic but not < logic
    if( k %in% keys( LE_HASH ) && k %nin% keys( LT_HASH ) )
        NRS <- c( NRS, paste0( k, "<=", min( LE_HASH[[ k ]] ) ) )

    # Similar logic as above, but for > and >= logic.
    if( k %in% keys( GT_HASH ) && k %in% keys( GE_HASH ) ){
      if( max( GT_HASH[[ k ]] ) >= max( GE_HASH[[ k ]] ) )
          NRS <- c( NRS,  paste0( k, ">", max( GT_HASH[[ k ]] ) ) )
      else
          NRS <- c( NRS,  paste0( k, ">=", max( GE_HASH[[ k ]] ) ) )
    }# end both > and >=
    
    if( k %in% keys( GT_HASH ) && k %nin% keys( GE_HASH ) )
        NRS <- c( NRS, paste0( k, ">", max( GT_HASH[[ k ]] ) ) )
    
    if( k %in% keys( GE_HASH ) && k %nin% keys( GT_HASH ) )
        NRS <- c( NRS, paste0( k, ">=", max( GE_HASH[[ k ]] ) ) )
    
    # For nominal split values use the intersection of all values in the hash
    # for each variable
    else if( k %in% keys( NOMINAL_IN_HASH ) ){
    
      shortest_split_value_length <- .Machine$integer.max
      shortest_subgroup <- NULL
      
      split_value_vector <- NOMINAL_IN_HASH[[ k ]]
      
      for( val in split_value_vector ){
        
        if( nchar( val ) < shortest_split_value_length ){
          shortest_split_value_length <- nchar( val )
          shortest_subgroup <- val
        }
      }
      NRS <- c( NRS, paste0( k, "%in%c(", shortest_subgroup, ")" ) )
    }


  # For nominal split values use the union of all split values in the hash for each variable
    else if( k %in% keys( NOMINAL_NIN_HASH ) ){
    
      longest_split_value_length <- 0
      longest_subgroup <- NULL
      
      split_value_vector <- NOMINAL_NIN_HASH[[ k ]]
      
      for( val in split_value_vector ){
        
        if( nchar( val ) == longest_split_value_length ){ #Should never happen
          stop( paste0( "ERROR: inconsistent subgroup definition:  ",
                       k, "%nin%c(", longest_subgroup, ") AND ",
                       k, "%nin%c(", val, ")." ) )
          
        }
        
        if( nchar( val ) > longest_split_value_length ){
          longest_split_value_length <- nchar( val )
          longest_subgroup <- val
        }
      }
      
      NRS <- c( NRS, paste0( k, "%nin%c(", longest_subgroup, ")" ) )
    }
  }
  
  return( sort( NRS ) )
}


is.binary <- function( x ){
  return( length( unique( x ) ) == 2 )
}

treatment_value <- function( trt, trt_control ){
  TRT_VALUES <- unique( trt )
  return( TRT_VALUES[ TRT_VALUES != trt_control ] )
}

#' @title reset_factor_levels
#' @description Reset the list of levels associated with a factor variable.
#' @details After subsetting a factor variable some factor levels that were
#' previously present may be lost. This is particularly true for relatively rare
#' factor levels. This function resets the list of factor levels to include only
#' the levels currently present.
#' @param data A data.frame containing factor variables.
#' @return A data.frame with factor variable that now have reset levels.
#' @examples
#' ex1 = as.factor( c( rep('A', 3), rep('B',3), rep('C',3) ) )
#' 
#' ## The levels associated with the factor variable include the letters A, B, C
#' ex1  # Levels are A, B, C
#' 
#' ## If the last three observations are dropped the value C no longer occurs
#' ## in the data, but the list of associated factor levels still contains C.
#' ## This mismatch between the data and the list of factor levels may cause
#' ## problems, particularly for algorithms that iterate over the factor levels.
#'
#' ex1 <- ex1[1:6]
#' ex1 # Levels are still A, B, C, but the data contains only A and B
#' 
#' ## If the factor levels are reset the data and list of levels will once again
#' ## be consistent
#' ex1 <- reset_factor_levels( ex1 )
#' ex1 # Levels now contain only A and B, which is consistent with data
#' @export
reset_factor_levels <- function( data ){

  data_is_factor <- FALSE
  
  if( is.factor( data ) ){
    data_is_factor <- TRUE
    data <- as.data.frame( data )#if data is just a single factor convert to data.frame
  }          
  
  for( j in 1:NCOL( data ) ){
    if( is.factor( data[,j] ) ){
      data[,j] <- factor( data[,j] )
    }
  }

  if( data_is_factor )
      data <- as.factor( data[,1] ) # convert from data.frame back to factor
  
  return( data )
}

#' @title unfactor
#' @description Convert the factor columns of a data.frame to character or numeric.
#' @details If the levels of a factor variable in data represent numeric values
#' the variable will be converted to a numeric data type, otherwise it is
#' converted to a character data type.
#' @seealso \link{na2empty}
#' @param data A factor variable or a data.frame containing factor variables.
#' @return A vector or data.frame no longer containing any factor variables.
#' @examples
#' ## Generate example data.frame of factors with factor levels of numeric,
#' ## character and mixed data types.
#' N <- 20
#' ex1 <- data.frame( factor( sample( c(0,1,NA), size = N, prob = c(0.4,0.3,0.3),
#'                            replace = TRUE ) )  )
#' names( ex1 ) <- "num"
#' ex1$char <- factor( sample( c("Control","Experimental", NA ), size = N,
#'                     prob = c(0.4,0.3,0.3), replace = TRUE ) )
#' ex1$mixed <- factor( sample( c(10,'A',NA), size = N, prob = c(0.4,0.3,0.3),
#'                      replace = TRUE ) )
#' 
#' ## Initially the data type of all variables in ex1 is factor
#' ex1
#' class( ex1$num )   #factor
#' class( ex1$char )  #factor
#' class( ex1$mixed ) #factor
#' 
#' ## Now convert all factor variables to numeric or character
#' ex2 <- unfactor( ex1 )
#' ex2
#' 
#' ## The data types are now numeric or character
#' class( ex2$num )   # numeric
#' class( ex2$char )  # character
#' class( ex2$mixed ) # character
#' 
#' ## The <NA> notation for missing factor values that have been converted to
#' ## character can be changed to an empty string for easier reading by use of
#' ## the function na2empty().
#' ex2$char <- na2empty( ex2$char )
#' ex2$mixed <- na2empty( ex2$mixed )
#' ex2
#' @export
unfactor <- function( data ){
  
  data_is_factor <- FALSE
  
  if( is.factor( data ) ){
    data_is_factor <- TRUE
    data <- as.data.frame( data )#if data is just a single factor convert to data.frame
  }
          
  for( j in 1:NCOL( data ) ){
    if( is.factor( data[,j] ) ){

      char_var <- as.character( data[,j] )
      num_var <- suppressWarnings( as.numeric( char_var ) )

      ## A factor should be converted to a character if any of the following
      ## conditions are true:
      ## 1. all elements are missing when converted to numeric
      ## 2. any element is non-missing when converted to character, but missing
      ##    when converted to numeric
      if( all( is.na( num_var ) ) || ( any( is.na( num_var ) & !is.na( char_var ) ) ) ){
        
        data[,j] <- char_var
        if( data_is_factor )# not a data.frame
            data <- as.character( data[,1] )
      }
      else{
        data[,j] <- num_var
        if( data_is_factor )# not a data.frame
            data <- as.numeric( data[,1] )
      }

      rm( char_var )
      rm( num_var )
    }
  }
  return( data ) 
}

trim <- function( x ){
  return( gsub( "^\\s+|\\s+$", "", x ) )
}


get_rpart_subgroup <- function( tree, node_id ){

  ## Create NULL placeholders to prevent NOTE in R CMD check
  NodeID <- NULL
  
  if( node_id %nin% tree$NodeID )
      stop( "ERROR: node_id not found in tree" )

  next_node <- node_id
  
  sg <- NULL
  
  PARENTS <- parents( node_id )

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

get_subset <- function( data, subgroup, digits = NULL ){

  # If digits is non-null then round to digits
  if( !is.null( digits ) ){

    splits__ <- strsplit( split = "&", subgroup )[[1]]

    for( j in 1:length( splits__ ) ){


      if( grepl( pattern = "<|>", perl = TRUE, splits__[[j]] ) ){

        splits__[[j]] <- gsub( pattern = "<", perl = TRUE, replacement = paste0( ",", digits, ")-1E-", digits, "<"), splits__[[j]] )
        splits__[[j]] <- gsub( pattern = ">", perl = TRUE, replacement = paste0( ",", digits, ")+1E-", digits, ">"), splits__[[j]] )
        splits__[[j]] <- paste0( "round(", splits__[[j]] )

      }
    }
    subgroup <- paste( splits__, collapse = "&" )
  }
  
  subset_str <- paste0( 'subset( ', quote( data ), ', ', subgroup, ' )' )
  return( eval( parse( text = subset_str ) ) )
}


get_subgroup_response <- function( data,
                                   response = NULL, y_var = 'y',
                                   trt = NULL, trt_var = NULL, trt_value = 'Control',
                                   subgroup = NULL,
                                   complement = FALSE ){

  SUBGROUP_RESPONSE <- NULL
  
  if( is.null( subgroup ) || subgroup == "" )
      stop( "ERROR: subgroup must not be null" )

  if( subgroup == "Overall" )
      subgroup <- "TRUE"

  # If a seperate response variable is provided add it to data
  if( !is.null( response ) )
      data[,y_var] <- response

  # If a seperate treatment variable is provided add it to data
  if( !is.null( trt ) ){
    
    if( is.null( trt_var ) )
      trt_var <- 'trt'
      
    data[,trt_var] <- trt
  }

  
  if( is.null( trt_var ) ){ # Subgroup with no treatment arm specified
    
    if( !complement )
        SUBGROUP_RESPONSE <- eval( parse( text = paste0( 'subset( data, ', subgroup, ', select = y_var )' ) ) )
    if( complement )
        SUBGROUP_RESPONSE <- eval( parse( text = paste0( 'subset( data, !(', subgroup, '), select = y_var )' ) ) )
    
  }else{ # Subgroup within a specified treatment arm
    
    if( trt_value %nin% data[,trt_var] )
        stop( paste0( "ERROR: ", trt_value, " not found in ", trt_var ) )
    if( !complement )
        subset_string <- paste0( 'subset( data, ', subgroup, ' & ', trt_var, ' == ', trt_value, ', select = y_var )' )
    if( complement )
        subset_string <- paste0( 'subset( data, !(', subgroup, ') & ', trt_var, ' == ', trt_value, ', select = y_var )' )
    
    SUBGROUP_RESPONSE <- eval( parse( text = subset_string ) )
  }
  
  return( SUBGROUP_RESPONSE[[1]] )
}

## END OF FILE
