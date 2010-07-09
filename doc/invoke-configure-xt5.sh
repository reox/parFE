#!/bin/bash

CC=cc 
CXX=CC 
CPP="CC -E"

INST_DIR="$HOME/parfe-build"
TRILINOS_DIR="/path/to/trilinos-10.2.2"
HDF5_DIR="/path/to/hdf5-1.8"
PARMETIS_DIR="/path/to/parmetis-3.1"

./configure \
  --prefix=$INST_DIR \
  --enable-mpi \
  --with-trilinos=$TRILINOS_DIR \
  --with-hdf5=$HDF5_DIR \
  --with-parmetis=$PARMETIS_DIR \
  --with-cxxflags="-DH5_USE_16_API"
