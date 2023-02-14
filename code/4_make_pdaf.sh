source 0_set_environment.sh

if [ "$ARCHER2" = true ] ; then
   cp -v $CODE_DIR/archer2-files/pdaf/cray_gfortran.h $PDAF_CLONE/make.arch/cray_gfortran.h  
   export PDAF_ARCH=cray_gfortran
fi

cd $PDAF_CLONE/src && make clean && make
cd $CODE_DIR
