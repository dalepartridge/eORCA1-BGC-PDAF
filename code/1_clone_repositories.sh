#!/bin/bash
  
# Clone code bases
source 0_set_environment.sh

svn checkout http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-2.5 $XIOS_CLONE

git clone https://github.com/pmlmodelling/ersem.git $ERSEM_CLONE
git clone https://github.com/FABM-MEDUSA/fabm-medusa.git $MEDUSA_CLONE

git clone https://github.com/fabm-model/fabm.git $FABM_CLONE

git clone https://github.com/pmlmodelling/NEMO4.0-FABM.git $NEMO_DIR

# insert PDAF directory here:
if [ "$USE_PDAF" = true ] ; then
   wget -O pdaf.tar.gz 'http://pdaf.awi.de/download/index.php?id=54f25d89da3f341c748c39966fb8f073&package=PDAF_V2.0.tar.gz'
   gzip -d pdaf.tar.gz
   tar -zxvf pdaf.tar
   rm pdaf.tar
   mv PDAF_V2.0 $PDAF_CLONE
   
   git clone https://github.com/pmlmodelling/NEMO4.0-FABM-PDAF.git $NEMO_PDAF_CLONE
fi
