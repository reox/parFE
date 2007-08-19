AC_INIT(parfe, 0.1, ETHZ)

echo "--> Configuring PARFE code..."

# This is to protect against accidentally specifying the wrong
# directory with --srcdir.  Any file in that directory will do,
# preferably one that is unlikely to be removed or renamed.

AC_CONFIG_SRCDIR(/src/Jacobi.cpp)
AC_CONFIG_AUX_DIR(config)
AM_CONFIG_HEADER(src/parfe_config.h:src/parfe_config.h.in)
AC_ARG_WITH(install,
 [AC_HELP_STRING([--with-install=INSTALL_PROGRAM],
 [Use the installation program INSTALL_PROGRAM rather the default that is provided.  For example --with-install="/path/install -p"])],
 [
   INSTALL=$withval
   INSTALL_PROGRAM=$withval
   INSTALL_SCRIPT=$withval
   INSTALL_DATA="$withval -m 644"
 ],)

                                                                            
# AM_MAINTAINER_MODE turns off maintainer-only makefile targets by
# default, and changes configure to understand a
# --enable-maintainer-mode option. --enable-maintainer-mode turns the
# maintainer-only targets back on. The maintainer-only makefile
# targets permit end users to clean automatically-generated files such
# as configure, which means they have to have autoconf and automake
# installed to repair the damage. AM_MAINTAINER_MODE makes it a bit
# harder for users to shoot themselves in the foot.

AM_MAINTAINER_MODE

# Define $build, $host, $target, etc

AC_CANONICAL_TARGET

# Use automake

#  - Required version of automake.
AM_INIT_AUTOMAKE(1.6.3)

# Specify required version of autoconf.

AC_PREREQ(2.52)

AC_PROG_F77(g95 f90 xlf90 ftn pgf90 gfortran)
# Satisfy newer automake versions
test -n "$FC" || FC=$F77
AC_SUBST(FC)

# Determine libraries needed for linking with Fortran
AC_F77_LIBRARY_LDFLAGS

ac_cv_flibs=$FLIBS
#Workaround for CRAY XT3 compiler wrappers
#remove confusing quotes
ac_cv_flibs=`echo $ac_cv_flibs | sed "s/'//g"`
# remove a -lcrtend' that sometimes appears on the XT3
ac_cv_flibs=`echo $ac_cv_flibs | sed "s/-lcrtend\'//g"`
#remove -lcrtend flags (I don't know why, but otherwise the program
#crashes)
ac_cv_flibs=`echo $ac_cv_flibs | sed "s/-lcrtend//g"`
# remove a -ldl" that appears with INTEL compilers under LINUX.
ac_cv_flibs=`echo $ac_cv_flibs | sed "s/-ldl\"//g"`
#finally set FLIBS to the new value of $ac_cv_flibs
test -n "$ac_cv_flibs" && FLIBS=$ac_cv_flibs

# ------------------------------------------------------------------------
# Check to see if MPI enabled and if any special configuration done
# ------------------------------------------------------------------------

TAC_ARG_CONFIG_MPI

# ------------------------------------------------------------------------
# Checks for programs
# ------------------------------------------------------------------------

AC_PROG_CC(cc gcc)
AC_PROG_CXX(CC g++ gcc c++ cxx)

AC_PROG_RANLIB
AC_PROG_INSTALL

# ------------------------------------------------------------------------
# Checks for extra flags
# ------------------------------------------------------------------------
# Check if --with-flags present, prepend any specs to FLAGS

TAC_ARG_WITH_FLAGS(ccflags, CCFLAGS)
TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
TAC_ARG_WITH_FLAGS(cflags, CFLAGS)
TAC_ARG_WITH_FLAGS(fflags, FFLAGS)
TAC_ARG_WITH_LIBS
TAC_ARG_WITH_FLAGS(ldflags, LDFLAGS)

# ------------------------------------------------------------------------
# Alternate archiver
# ------------------------------------------------------------------------

TAC_ARG_WITH_AR

# ------------------------------------------------------------------------
# MPI link check
# ------------------------------------------------------------------------

TAC_ARG_CHECK_MPI

# ------------------------------------------------------------------------
# Checks for header files
# ------------------------------------------------------------------------

dnl AC_CHECK_HEADERS([cstdlib stdlib.h], [break],)
dnl AC_CHECK_HEADERS([cstdio stdio.h], [break],)
dnl AC_CHECK_HEADERS([cassert assert.h], [break],)
dnl AC_CHECK_HEADERS([cmath math.h], [break],)
dnl AC_CHECK_HEADERS([iomanip iomanip.h], [break],)
dnl AC_CHECK_HEADERS([iostream iostream.h], [break],)
dnl AC_CHECK_HEADERS([string string.h], [break],)

dnl if test "X$ac_cv_header_cstdlib" != "Xyes" && test "X$ac_cv_header_stdlib_h" !=
dnl "Xyes"; then
dnl   AC_MSG_ERROR([[Cannot find cstdlib or stdlib.h.  New_Package requires at least
dnl  one or the other.]])
dnl fi

dnl if test "X$ac_cv_header_cstdio" != "Xyes" && test "X$ac_cv_header_stdio_h" != "X
dnl yes"; then
dnl   AC_MSG_ERROR([[Cannot find cstdio or stdio.h.  New_Package requires at least o
dnl ne or the other.]])
dnl fi

dnl if test "X$ac_cv_header_cassert" != "Xyes" && test "X$ac_cv_header_assert_h" !=
dnl "Xyes"; then
dnl   AC_MSG_ERROR([[Cannot find cassert or assert.h.  New_Package requires at least
dnl  one or the other.]])
dnl fi

dnl if test "X$ac_cv_header_string" != "Xyes" && test "X$ac_cv_header_string_h" != "
dnl Xyes"; then
dnl   AC_MSG_ERROR([[Cannot find string or string.h.  New_Package requires at least
dnl one or the other.]])
dnl fi

dnl if test "X$ac_cv_header_iostream" != "Xyes" && test "X$ac_cv_header_iostream_h"
dnl != "Xyes"; then
dnl   AC_MSG_ERROR([[Cannot find iostream or iostream.h.  New_Package requires at le
dnl ast one or the other.]])
dnl fi

dnl if test "X${host}" != "Xtflop"; then
dnl   if test "X$ac_cv_header_cmath" != "Xyes" && test "X$ac_cv_header_math_h" != "X
dnl yes"; then
dnl     AC_MSG_ERROR([[Cannot find cmath or math.h.  New_Package requires at least o
dnl ne or the other.]])
dnl   fi
dnl else
dnl   # TFLOP
dnl   if test "X$ac_cv_header_iomanip" != "Xyes" && test "X$ac_cv_header_iomanip_h"
dnl != "Xyes"; then
dnl     AC_MSG_ERROR([[Cannot find iomanip or iomanip.h.  New_Package requires at le
dnl ast one or the other.]])
dnl   fi
dnl fi


# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# ------------------------------------------------------------------------
# Checks for compiler characteristics - Epetra_ConfigDefs.h uses these, 
# so it makes sense to leave them in - see PORTABILITY in this directory
# ------------------------------------------------------------------------

AC_CXX_NAMESPACES
AC_CXX_HAVE_STL
AC_CXX_BOOL
AC_CXX_MUTABLE

# ---------------------------------------------------------------------
# These are pretty obscure
# ---------------------------------------------------------------------

AC_CXX_NEW_FOR_SCOPING
AC_CXX_STD_SPRINTF

# ------------------------------------------------------------------------
# Specify other directories 
# ------------------------------------------------------------------------


# enable use of --with-libdirs="-Llibdir1 -Llibdir2 ..." to prepend to LDFLAGS
TAC_ARG_WITH_LIBDIRS
# enable use of --with-incdirs="-Iincdir1 -Iincdir2 ..." to prepend to CPPFLAGS
TAC_ARG_WITH_INCDIRS

AC_CHECK_LIB([f2c], [f_exit], ,
	[#if not found, try g2c
	AC_CHECK_LIB([g2c], [f_exit],,)]
)

#workaround for fortran pgf compilers
dnl AC_CHECK_LIB([pghpf_mpi],[__hpf_exit],[FLIBS="$FLIBS -lpghpf_mpi $FLIBS"],,[$FLIBS])

# Define F77_FUNC that will be used to link with Fortran subroutines, e.g., BLAS/LAPACK - trash WORKGXX 
AC_F77_WRAPPERS

# Find math library
AC_SEARCH_LIBS(pow,[m],,AC_MSG_ERROR(Cannot find math library))
AC_SEARCH_LIBS(sqrt,[m],,AC_MSG_ERROR(Cannot find math library))

TAC_ARG_WITH_BLASLIB
TAC_ARG_WITH_LAPACKLIB

ACX_BLAS([], [AC_MSG_ERROR([[Cannot find blas library, specify a path using --with-ldflags=-I<DIR> (ex. --with-ldflags=-I/usr/path/lib) or a specific library using --with-blas=DIR/LIB (for example --with-blas=/usr/path/lib/libcxml.a)]])])

ACX_LAPACK([], [AC_MSG_ERROR([[Cannot find lapack library, specify a path using --with-ldflags=-I<DIR> (ex. --with-dlflags=-I/usr/path/lib) or a specific library using --with-lapack=DIR/LIB (ex. --with-lapack=/usr/path/lib/libcxml.a)]])])

ACX_SUPERLU([],[])

ACX_PARMETIS([], [AC_MSG_ERROR([[Cannot find ParMETIS library, specify a path using --with-parmetis=DIR (ex. --with-parmetis=/opt/ParMETIS)]])])

ACX_TRILINOS([], [AC_MSG_ERROR([[Cannot find all necessary trilinos libraries, specify a path using --with-trilinos=DIR (ex. --with-trilinos=/opt/trilinos)]])])

ACX_HDF5([],[])

# ------------------------------------------------------------------------
# Checks for linker characteristics
# ------------------------------------------------------------------------

AC_MSG_CHECKING(for mallinfo() function)
AC_TRY_LINK_FUNC(mallinfo,,mallinfook="yes",mallinfook="no")
if test "X$mallinfook" != "Xyes"; then
  AC_DEFINE(HAVE_MALLINFO,[],[mallinfo])  AC_MSG_RESULT(yes)
else
  AC_MSG_RESULT(no)
fi 

AC_CHECK_FUNC(heap_info, [AC_DEFINE([HAVE_HEAP_INFO],[],[heap_info])])


# ------------------------------------------------------------------------
# doxygen support, Copyright (C) 2004 Oren Ben-Kiki.
# ------------------------------------------------------------------------

AM_PROG_CC_C_O

DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(ON)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)
DX_INIT_DOXYGEN(parfe, doc/Doxyfile)


# ------------------------------------------------------------------------
# Perform substitutions in output files
# ------------------------------------------------------------------------

AC_SUBST(ac_aux_dir)

# ------------------------------------------------------------------------
# Output files
# ------------------------------------------------------------------------

AC_CONFIG_FILES([
		Makefile
		src/Makefile
		example/Makefile
                test/Makefile
		])

AC_OUTPUT()

echo "--> Finished Running PARFE Configure Script"
