#!/bin/bash
./configure \
  --with-mpi-compilers \
  --with-trilinos=/home/cflaig/builds/trilinos-10.2.2/ \
  --with-cxxflags="-DH5_USE_16_API"
