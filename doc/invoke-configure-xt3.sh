./configure \
--prefix=$PWD \
CC=cc  \
CXX=CC \
CPP="CC -E" \
FC=ftn \
--host=x86_64-unknown-linux-gnu \
--enable-static \
--with-ldflags="-L/apps/zlib/lib" \
--with-hdf5=/apps/hdf5 \
--with-trilinos=$HOME/Trilinos/CRAY_XT3 \
--with-parmetis=/apps/metis/parmetis-3.1_PE1.4.10 \
--with-cxxflags="-DMPICH_IGNORE_CXX_SEEK"
