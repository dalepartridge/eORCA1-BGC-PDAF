#!/bin/bash

#************* CHANGE TO WORKING AREA *******************
export WORK=/work/n01/n01/dapa/NCEO/eORCA1-BGC-PDAF
#********************************************************

#Config options
export CODE_DIR=$WORK/code
export ARCHER2=true

#Load modules
source $CODE_DIR/archer2-files/ucx_env

#XIOS options
export XIOS_CLONE=$CODE_DIR/xios
export XIOS_HOME=$CODE_DIR/xios-build
export CC=cc export CXX=CC export FC=ftn export F77=ftn export F90=ftn
export XIOS_ARCH=GCC_ARCHER2

#BGC options
export ERSEM_CLONE=$CODE_DIR/ersem
export MEDUSA_CLONE=$CODE_DIR/medusa

#PDAF Options
export PDAF_CLONE=$CODE_DIR/pdaf
export PDAF_MYSRC=$CODE_DIR/pdaf-src

#FABM options
export FABM_CLONE=$CODE_DIR/fabm
export FABM_BUILD=$CODE_DIR/fabm-build
export FABM_DEBUG=$CODE_DIR/fabm-debug

#NEMO options
export FABM_HOME=$FABM_BUILD 
export NEMO_DIR=$CODE_DIR/nemo
export NEMO_ARCH=GCC_ARCHER2
export NEMO_REF=eORCA1
export NEMO_CFG=eORCA1-build
