#!/bin/bash
#################################################################################
# FILENAME   : lilly_clean.sh
# AUTHOR     : Brian Denton <denton_brian_david@lilly.com>
# DATE       : 2018-11-11
# PROJECT    : 
# DESCRIPTION: 
#################################################################################

TMPDIR=$(/bin/mktemp --directory --tmpdir=.)

#TMPDIR="/home/c065288/Projects/CRAN_TSDT/temp"
#mkdir --parents $TMPDIR

echo $TMPDIR

TARBALL=$(/bin/ls TSDT_*.tar.gz)

## Un-tar
/bin/tar -xvf $TARBALL --directory=$TMPDIR

## Remove Lilly userid
/bin/sed --in-place 's/c065288/Brian David Denton/g' $TMPDIR/TSDT/DESCRIPTION

## Re-tar
cd $TMPDIR
/bin/tar -cvf $TARBALL TSDT
cd /home/c065288/Projects/CRAN_TSDT

## Replace tarball in project folder
/bin/mv $TMPDIR/$TARBALL .

## Delete temp folder
/bin/rm -rf $TMPDIR

## END OF FILE
