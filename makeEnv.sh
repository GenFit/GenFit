#!/bin/bash
ENVFILE=env.sh
if [ -z $GENFIT ]; then
  echo "define the GENFIT env variable before executing this script"
else
  echo "if [ -z \$ROOTSYS ]; then" > $ENVFILE
    echo "echo \"ROOTSYS is not set. Check your ROOT installation.\"" >> $ENVFILE
  echo "else" >> $ENVFILE
    echo "export GENFIT=$GENFIT" >> $ENVFILE
    echo "if [ \`root-config --arch\` = macosx ]; then" >>$ENVFILE
      echo "export DYLD_LIBRARY_PATH=\$DYLD_LIBRARY_PATH:\$GENFIT/lib" >>$ENVFILE
      echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$GENFIT/lib" >>$ENVFILE
    echo "else" >>$ENVFILE
      echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$GENFIT/lib" >>$ENVFILE
    echo "fi" >>$ENVFILE
  echo "fi" >> $ENVFILE
  
  if [ $RAVEPATH ]; then
    echo "export RAVEPATH=$RAVEPATH" >> $ENVFILE
  fi
fi

