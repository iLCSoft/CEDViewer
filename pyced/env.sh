#!/bin/bash
####################################################
# set the environment for pyced
# call with ilcsoft init script as first argument
#
# F.Gaede, DESY
####################################################

if [ "$#" -ne 1 ]
then
  echo "usage:  . env.sh _path_to_ilcsoft_init.sh"
else

ILCSOFTINITFILE=$1

if [ -f $ILCSOFTINITFILE ];
then

source $ILCSOFTINITFILE

# find directory with CED lib
CEDLIB=$(which glced)
CEDLIB=${CEDLIB%/*}
CEDLIB=${CEDLIB%/*}/lib

#echo ${CEDLIB}

# finn directory with MarlinUtil lib
MULIB=$(find $ILCSOFT/MarlinUtil -maxdepth 1 -mindepth 1 -type d)
MULIB=${MULIB}/lib

#echo $MULIB

export PYTHONPATH=${LCIO}/src/python:${ROOTSYS}/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CEDLIB:$MULIB
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$CEDLIB:$MULIB

else
  echo "usage:  . env.sh _path_to_ilcsoft_init.sh"
fi
fi

