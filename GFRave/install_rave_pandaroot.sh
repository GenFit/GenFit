#!/bin/bash

# check if necessary environment variables are set
: ${SIMPATH:?"SIMPATH is not set!"}
: ${PANDABUILD:?"PANDABUILD is not set!"}
: ${RAVEPATH:?"RAVEPATH is not set!"}

# set important variables for rave installation
export CLHEPPATH=$SIMPATH/cern/clhep 
export CLHEPLIBPATH=$SIMPATH/cern/clhep/lib

export CLHEP_VECTORPATH=$SIMPATH/cern/clhep/Vector
export CLHEP_VECTORLIBPATH=$SIMPATH/cern/clhep/lib

export CLHEP_MATRIXPATH=$SIMPATH/cern/clhep/Matrix
export CLHEP_MATRIXLIBPATH=$SIMPATH/cern/clhep/lib



cd $RAVEPATH


buildRave=false

echo "make rave:"

if [ ! -f configure ] # bootstrap rave if necessary
then
    echo "bootsrap rave"
    ./bootstrap
    buildRave=true
elif [ ! -f installTimeCheck ] # check if file installTimeCheck is there
then
    echo "build rave for the first time"
    buildRave=true
else # check for new or modified files
  echo "checking for new/modified files"
  NEWFILES=`find . -path '*\.svn*' -prune -o -newer installTimeCheck -print`
  if [ ${#NEWFILES} != 0 ]
  then
    #echo $NEWFILES
    echo "found file newer than last build -> rebuild rave"
    buildRave=true
  fi
fi


if $buildRave
then
  # configure make make install rave
  ./configure --prefix=$PANDABUILD
  make
  make install
  
  touch installTimeCheck # create (or update) installTimeCheck file
else
  echo "nothing to be done"
fi

cd -


