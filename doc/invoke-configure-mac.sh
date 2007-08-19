./configure \
--prefix=$PWD \
CXX=mpic++ \
CXXFLAGS=-g \
--with-cxxflags="-DLAM_BUILDING" \
--with-incdirs="-I$HOME/include -I$HOME/Trilinos/G4_MPI/include" \
--with-ldflags="-L$HOME/lib -L$HOME/Trilinos/G4_MPI/lib -L$HOME/lib -L/sw/lib" \
--with-libs="-lml -lifpack -laztecoo -lamesos \
-lepetraext -lepetra -lteuchos -framework vecLib -lSystemStubs -lf2c -lexpat"
