#!/bin/sh

tar -xzvf mpich-4.1.2.tar.gz
MPICH_DIR=/scratch/mpich-4.1.2

MPICH_INSTALL_DIR=$MPICH_DIR/install
cd $MPICH_DIR
# from petsc 3.20 --download-mpich=yes configure  for mpich 4.1.2
./configure --prefix=$MPICH_INSTALL_DIR MAKE=/usr/bin/gmake --libdir=$MPICH_INSTALL_DIR/lib CC=gcc CFLAGS="-fPIC -Wno-lto-type-mismatch -Wno-stringop-overflow -g -O" AR=/usr/bin/ar ARFLAGS=cr CXX=g++ CXXFLAGS="-Wno-lto-type-mismatch -Wno-psabi -g -O -std=gnu++20" FFLAGS="-ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -g -O -fallow-argument-mismatch" FC=gfortran F77=gfortran FCFLAGS="-ffree-line-length-none -ffree-line-length-0 -Wno-lto-type-mismatch -g -O -fallow-argument-mismatch" --disable-shared --with-pm=hydra --disable-java --with-hwloc=embedded --with-device=ch3:nemesis --enable-g=meminit PYTHON="/usr/bin/python3" --disable-maintainer-mode --disable-dependency-tracking
make all; make install
