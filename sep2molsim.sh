#!/bin/bash

SEPLIB_PATH=~/software/seplib/

cp -r $SEPLIB_PATH/include/* src/include/
cp -r $SEPLIB_PATH/source/* src/source/
cp -r $SEPLIB_PATH/tools/* src/tools/
cp -r $SEPLIB_PATH/tools/molconfigs/* resources/
