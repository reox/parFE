#!/bin/bash
./configure \
  --with-mpi-compilers \
  --with-trilinos=$HOME/builds/trilinos-10.6.4/ \
  --with-cxxflags="-g -march=core2 -O3" \
  --with-cflags="-g -march=core2 -O3"
