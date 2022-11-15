#!/bin/bash

############# COMPILE NEMO ##########################
cd $NEMO_DIR

if [ "$ARCHER2" = true ] ; then
yes | cp $CODE_DIR/archer2-files/nemo/Config_gcc.pm $NEMO_DIR/ext/FCM/lib/Fcm/Config.pm
yes | cp -r $CODE_DIR/archer2-files/eORCA1 $NEMO_DIR/cfgs/
fi

printf 'y\nn\nn\ny\nn\nn\nn\nn\n' |./makenemo -r $NEMO_REF -m $NEMO_ARCH -n $NEMO_CFG -j 0
./makenemo -r $NEMO_REF -m $NEMO_ARCH -n $NEMO_CFG -j 4 clean
./makenemo -m $NEMO_ARCH -r $NEMO_REF -n $NEMO_CFG -j 8

