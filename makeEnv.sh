#!/bin/bash
ENVFILE=env.sh
if [ -z $GENFIT ]; then
  echo "define the GENFIT env variable before executing this script"
elif [ -z $VMC ]; then
  echo "WARNING: the environment variable VMC is not set. You will not be able to compile/use the GeaneTrackRep2!"
  echo "if [ -z \$ROOTSYS ]; then" > $ENVFILE
    echo "echo \"ROOTSYS is not set. Check your ROOT installation.\"" >> $ENVFILE
  echo "else" >> $ENVFILE
    echo "export GENFIT=$GENFIT" >> $ENVFILE
    echo "export VMC=\$VMC" >> $ENVFILE
    echo "if [ \`root-config --arch\` = macosx ]; then" >>$ENVFILE
      echo "export DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:\$VMC/lib/tgt_macosx:\$GENFIT/lib" >>$ENVFILE
      echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$VMC/lib/tgt_macosx:\$GENFIT/lib" >>$ENVFILE
    echo "else" >>$ENVFILE
      echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$VMC/lib/tgt_linux:\$GENFIT/lib" >>$ENVFILE
    echo "fi" >>$ENVFILE
    echo "if [ -z \$G4LIB ]; then" >> $ENVFILE
      echo "echo \"G4LIB (needed for G4eTrackRep) is not set. Run env.sh of your Geant4 installation.\"" >> $ENVFILE
    echo "else" >> $ENVFILE
      echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$G4LIB:\$G4LIB/\$G4SYSTEM" >>$ENVFILE
    echo "fi" >> $ENVFILE
    echo "if [ -z \$CLHEP_LIB_DIR ]; then" >> $ENVFILE
      echo "echo \"CLHEP_LIB_DIR (needed for G4eTrackRep) is not set. Run env.sh of your Geant4 installation.\"" >> $ENVFILE
    echo "else" >> $ENVFILE
      echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$CLHEP_LIB_DIR" >>$ENVFILE
    echo "fi" >> $ENVFILE
  echo "fi" >> $ENVFILE
fi

