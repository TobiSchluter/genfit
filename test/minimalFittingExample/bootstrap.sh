#!/bin/sh

root -q -l makeGeom.C
mkdir build
cd build
ln -s ../genfitGeom.root
cmake ..
make
./main
