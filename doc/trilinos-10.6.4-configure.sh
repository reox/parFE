#!/bin/bash
CC=mpicc
CXX=mpicxx
CPP="mpicxx -E" 
F77=mpif77 

cmake \
           --prefix=$HOME/builds/trilinos-10.6.4 \
           -DCMAKE_INSTALL_PREFIX:PATH=$HOME/builds/trilinos-10.6.4 \
		   -DCMAKE_CXX_FLAGS:STRING="-DMPICH_IGNORE_CXX_SEEK -fPIC -march=core2 -O2 -g" \
           -DCMAKE_C_FLAGS:STRING="-DMPICH_IGNORE_CXX_SEEK -fPIC -march=core2 -O2 -g" \
           -DCMAKE_Fortran_FLAGS:STRING="-fPIC" \
           -D CMAKE_BUILD_TYPE:STRING=DEBUG \
           -D TPL_ENABLE_MPI:BOOL=ON \
           -D TPL_ENABLE_BLAS:BOOL=ON \
           -D TPL_ENABLE_LAPACK:BOOL=ON \
           -D TPL_ENABLE_ParMETIS:BOOL=ON \
           -D TPL_ENABLE_METIS:BOOL=ON \
           -D TPL_METIS_INCLUDE_DIRS:PATH=/usr/include/metis \
           -D Trilinos_ENABLE_Epetra:BOOL=ON \
           -D Trilinos_ENABLE_EpetraExt:BOOL=ON \
           -D Trilinos_ENABLE_Tpetra:BOOL=ON \
           -D Trilinos_ENABLE_Ifpack:BOOL=ON \
           -D Trilinos_ENABLE_ML:BOOL=ON \
           -D Trilinos_ENABLE_Amesos:BOOL=ON \
           -D Trilinos_ENABLE_AztecOO:BOOL=ON \
           -D Trilinos_ENABLE_Teuchos:BOOL=ON \
           -D Trilinos_ENABLE_Aztecoo-Teuchos:BOOL=ON \
           -D Trilinos_ENABLE_Teuchos-Extended:BOOL=ON \
           -D Trilinos_ENABLE_Isorropia:BOOL=ON \
           -D Trilinos_ENABLE_Isorropia-Epetraext:BOOL=ON \
           -D Trilinos_ENABLE_TESTS:BOOL=OFF \
           ../trilinos-10.6.4-Source/
