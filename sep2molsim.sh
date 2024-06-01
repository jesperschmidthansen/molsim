#!/bin/bash

if [[ $1 == "dirac" ]]; then
	SEPLIB_PATH=/net/dirac/jschmidt/software/seplib
else
	SEPLIB_PATH=/home/jschmidt/software/seplib
fi

cp -r $SEPLIB_PATH/include/* src/include/
cp -r $SEPLIB_PATH/source/* src/source/
cp -r $SEPLIB_PATH/tools/* src/tools/
cp -r $SEPLIB_PATH/tools/molconfigs/* resources/
cp -r $SEPLIB_PATH/cuda/* cuda/
