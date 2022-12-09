#!/bin/bash
  
# Clone code bases
source 0_set_environment.sh

svn checkout http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5 $XIOS_CLONE

git clone https://github.com/pmlmodelling/ersem.git $ERSEM_CLONE
git clone https://github.com/FABM-MEDUSA/fabm-medusa.git $MEDUSA_CLONE

git clone https://github.com/fabm-model/fabm.git $FABM_CLONE

git clone https://github.com/pmlmodelling/NEMO4.0-FABM.git $NEMO_DIR

