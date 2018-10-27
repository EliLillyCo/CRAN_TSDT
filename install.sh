#!/bin/bash
#################################################################################
# FILENAME   : install.sh
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 2015-09-08
# PROJECT    : 
# DESCRIPTION: 
#################################################################################

#R_LIBS=/home/c065288/R/x86_64-unknown-linux-gnu-library
R_LIBS=/home/c065288/Projects/CRAN_TSDT/packrat/lib/x86_64-pc-linux-gnu/3.5.0
TSDT=/home/c065288/Projects/CRAN_TSDT/TSDT_*.tar.gz

export R_LIBS=$R_LIBS
R CMD INSTALL -l $R_LIBS $TSDT

# export R_LIB=/home/c065288/R/x86_64-unknown-linux-gnu-library
# export R_LIBS=/home/c065288/R/x86_64-unknown-linux-gnu-library

# library( TSDT, lib.loc = '/home/c065288/R/x86_64-unknown-linux-gnu-library' )

# library( BSID, lib.loc = '/lrlhps/apps/R/R_library_latest' )

# remove.packages( 'TSDT', lib = '/home/c065288/R/x86_64-unknown-linux-gnu-library' )

#Sys.setenv( R_LIBS = '/lrlhps/apps/R/R_library_latest' )

## END OF FILE
