#!/bin/bash

# Build FABM with cray compiler

mkdir -p $FABM_BUILD
cd $FABM_BUILD
cmake $FABM_CLONE/src -DFABM_HOST=nemo -DFABM_INSTITUTES='ersem;medusa' \
                    -DFABM_ERSEM_BASE=$ERSEM_CLONE -DFABM_MEDUSA_BASE=$MEDUSA_CLONE \
            -DFABM_EMBED_VERSION=ON -DCMAKE_INSTALL_PREFIX=$FABM_BUILD -DCMAKE_Fortran_COMPILER=ftn
make
make install -j4
cd $CODE_DIR
