#!/bin/bash
ENVFILE=env.sh
if [ -z $GENFIT ]; then
echo "define the GENFIT env variable before executing this script"
else
echo "if [ -z \$ROOTSYS ]; then" > $ENVFILE
echo "echo \"ROOTSYS is not set. Check your ROOT installation.\"" >> $ENVFILE
echo "else" >> $ENVFILE
echo "export GENFIT=$GENFIT" >> $ENVFILE
echo "export VMC=\$GENFIT/geant3" >> $ENVFILE
echo "if [ \`root-config --arch\` = macosx ]; then" >>$ENVFILE
echo "export DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:\$VMC/lib/tgt_macosx:\$GENFIT/lib" >>$ENVFILE
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$VMC/lib/tgt_macosx:\$GENFIT/lib" >>$ENVFILE
echo "else" >>$ENVFILE
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$VMC/lib/tgt_linux:\$GENFIT/lib" >>$ENVFILE
echo "fi" >>$ENVFILE
echo "fi" >> $ENVFILE
fi