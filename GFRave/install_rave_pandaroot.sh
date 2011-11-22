export CLHEPPATH=$SIMPATH/cern/clhep 
export CLHEPLIBPATH=$SIMPATH/cern/clhep/lib

export CLHEP_VECTORPATH=$SIMPATH/cern/clhep/Vector
export CLHEP_VECTORLIBPATH=$SIMPATH/cern/clhep/lib

export CLHEP_MATRIXPATH=$SIMPATH/cern/clhep/Matrix
export CLHEP_MATRIXLIBPATH=$SIMPATH/cern/clhep/lib

cd $RAVEPATH
./configure --prefix=$PANDABUILD
make
make install
