#!usr/bin/dev/R
#################################################################################
# FILENAME   : Create_Package.R
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 02/25/2013
# DESCRIPTION: Create and populate package skeleton and generate help files.
#################################################################################

options( warn = 2 )

library( devtools )
library( roxygen2 )

ProjectPath <- '/home/c065288/Projects/CRAN_TSDT/'
SRC <- paste( ProjectPath, 'R', sep = "" )

system( paste0( 'rm -rf ', ProjectPath, 'TSDT' ) )

R_FILES <- list.files( path = SRC, pattern = "*.R", full.names = TRUE )

#################################################################################
# Copy R files to temp folder and remove source(...) lines
#################################################################################

TMP <- tempdir()

cmd1 <- paste( 'cp', paste( R_FILES, collapse = ' ' ), TMP, collapse = ' ' ) 
cmd2 <- paste0( 'chmod ugo+rw ', TMP, '/*.*' )
cmd3 <- paste0( "sed --in-place '/source(/d' ", TMP, '/*.R' )

system( cmd1 )
system( cmd2 )
system( cmd3 )

R_FILES <- list.files( path = TMP, pattern = "*.R", full.names = TRUE )

#################################################################################
# Create package skeleton                                                       #
#################################################################################

package.skeleton( name = "TSDT", environment = "package:TSDT",
                  path = ProjectPath, force = TRUE,
                  code_files = R_FILES )

#################################################################################
# Copy original source code into package skeleton                               #
#################################################################################

cmd4 <- paste0( "cp ", paste( R_FILES, collapse = ' ' ), ' ', ProjectPath, "TSDT/R/" )

cmd5.1 <- paste0( "cp ", ProjectPath, 'Manually_Edited_Files/DESCRIPTION ', ProjectPath, "TSDT/" )
cmd5.2 <- paste0( "cp ", ProjectPath, 'Manually_Edited_Files/NAMESPACE ', ProjectPath, "TSDT/" )
cmd5.3 <- paste0( "cp ", ProjectPath, 'Manually_Edited_Files/CHANGELOG ', ProjectPath, "TSDT/" )
cmd5.4 <- paste0( "cp ", ProjectPath, 'Manually_Edited_Files/*.Rd ', ProjectPath, "TSDT/man" )

cmd6 <- paste0( "chmod --recursive u+rw ", ProjectPath, "TSDT" )
cmd7 <- paste0( "rm -f ", ProjectPath, "TSDT/Read-and-delete-me" )

cat( "Copying *.R files\n" )
system( cmd4 )

cat( "Copying manually edited files\n" )
system( cmd5.1 )
system( cmd5.2 )
system( cmd5.3 )
system( cmd5.4 )

cat( "Make TSDT package folder writeable\n" )
system( cmd6 )

cat( "Deleting Read-and-delete-me\n" )
system( cmd7 )

#################################################################################
# Remove unwanted help files                                                    #
#################################################################################

## HELP_FILES_TO_RETAIN <- c(
##     'binary_transform.Rd',
##     'bootstrap.Rd',
##     'ctree_wrapper.Rd',
##     'desirable_response_proportion.Rd',
##     'diff_mean_deviance_residuals.Rd',
##     'diff_quantile_response.Rd',
##     'diff_restricted_mean_survival_time.Rd',
##     'diff_survival_time_quantile.Rd',
##     'distribution.Rd',
##     'get_covariates.Rd',
##     'get_trt.Rd',
##     'get_y.Rd',
##     'mean_deviance_residuals.Rd',
##     'mean_response.Rd',
##     'mob_wrapper.Rd',
##     'parse_party.Rd',
##     'parse_rpart.Rd',
##     'permutation.Rd',
##     'quantile_response.Rd',
##     'rpart_nodes.Rd',
##     'rpart_wrapper.Rd',
##     'subgroup.Rd',
##     'subsample.Rd',
##     'survival_time_quantile.Rd',
##     'treatment_effect.Rd',
##     'TSDT.Rd',
##     'unpack_args.Rd',
##     'Bootstrap.ClassDefinition.Rd',
##     'BootstrapStatistic.ClassDefinition.Rd',
##     'Subsample.ClassDefinition.Rd',
##     'TSDT.ClassDefinition.Rd',
##     'TSDT_Sample.ClassDefinition.Rd',
##     'TSDT_CutpointDistribution.ClassDefinition.Rd'
##  )


HELP_FILES_TO_RETAIN <- NULL

RD_FILES <- list.files( path = paste0( ProjectPath, 'TSDT/man' ), pattern = '*.Rd' )

HELP_FILES_TO_REMOVE <- setdiff( RD_FILES, HELP_FILES_TO_RETAIN )

cmd8 <- paste0( 'rm -rf ', paste( ProjectPath, 'TSDT/man/', HELP_FILES_TO_REMOVE, sep = '', collapse = ' ' ) )

cat( "Remove unwanted help files\n" )
system( cmd8 )



#################################################################################
# Generate help files and build package                                         #
#################################################################################

roxygenise( paste0( ProjectPath, "TSDT" ) )

## devtools::load_all( paste0( ProjectPath, "TSDT" ) )
#devtools::document( paste0( ProjectPath, "TSDT" ) )

cmd9 <- paste0( "cp ", ProjectPath, 'Manually_Edited_Files/nin.Rd ', ProjectPath, "TSDT/man/grapes-nin-grapes.Rd" )
system( cmd9 )
    
build( paste0( ProjectPath, "TSDT" ) )

## END OF FILE
