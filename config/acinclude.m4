dnl acinclude.m4

dnl Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
dnl Free Software Foundation, Inc.
dnl This file is free software; the Free Software Foundation
dnl gives unlimited permission to copy and/or distribute it,
dnl with or without modifications, as long as this notice is preserved.

dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY, to the extent permitted by law; without
dnl even the implied warranty of MERCHANTABILITY or FITNESS FOR A
dnl PARTICULAR PURPOSE.

dnl Like AC_CONFIG_HEADER, but automatically create stamp file. -*- Autoconf -*-

dnl Copyright 1996, 1997, 2000, 2001 Free Software Foundation, Inc.

dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 2, or (at your option)
dnl any later version.

dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.

dnl You should have received a copy of the GNU General Public License
dnl along with this program; if not, write to the Free Software
dnl Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
dnl 02111-1307, USA.

AC_PREREQ([2.52])

dnl check wether non-deprecated macro AS_HELP_STRING exists
m4_ifdef([AS_HELP_STRING], , [m4_define([AS_HELP_STRING], m4_defn([AC_HELP_STRING]))])


dnl @synopsis TAC_ARG_CONFIG_MPI
dnl
dnl Test a variety of MPI options:
dnl --enable-mpi       - Turns MPI compiling mode on
dnl --with-mpi         - specify root directory of MPI
dnl --with-mpi-compilers - Turns on MPI compiling mode and sets the MPI C++
dnl                       compiler = mpicxx, mpic++ or mpiCC,
dnl                       the MPI C compiler = mpicc and 
dnl                       the MPI Fortran compiler = mpif77
dnl --with-mpi-incdir - specify include directory for MPI 
dnl --with-mpi-libs    - specify MPI libraries
dnl --with-mpi-libdir  - specify location of MPI libraries
dnl
dnl If any of these options are set, HAVE_MPI will be defined for both
dnl Autoconf and Automake, and HAVE_MPI will be defined in the
dnl generated config.h file
dnl
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl edited by Uche Mennel <umennel@student.ethz.ch> since the parfe package
dnl requires MPI. Also prevented that compiler variables (CC, CXX, F77)
dnl are overwritten if already set.
AC_DEFUN([TAC_ARG_CONFIG_MPI],
[

# AC_ARG_ENABLE(mpi,
# [AC_HELP_STRING([--enable-mpi],[MPI support])],
# [HAVE_PKG_MPI=$enableval],
# [HAVE_PKG_MPI=no]
# )

#umennel: MPI is required!
HAVE_PKG_MPI=yes

AC_ARG_WITH(mpi-compilers,
[AC_HELP_STRING([--with-mpi-compilers=PATH],
[use MPI compilers mpicc, mpif77, and mpicxx, mpic++ or mpiCC in the specified path or in the default path if no path is specified. Enables MPI])],
[
  if test X${withval} != Xno; then
    HAVE_PKG_MPI=yes
    if test X${withval} = Xyes; then
      # Check for mpicxx, if it does not exist, check for mpic++, if it does 
      # not exist, use mpiCC instead.
      AC_CHECK_PROG(MPI_TEMP_CXX, mpicxx, mpicxx, no)
      if test X${MPI_TEMP_CXX} = Xno; then
	AC_CHECK_PROG(MPI_CXX, mpic++, mpic++, mpiCC)
      else 
	MPI_CXX=${MPI_TEMP_CXX}
      fi
      MPI_CC=mpicc
      MPI_F77=mpif77
    else
      if test -f ${withval}/mpicxx; then
        MPI_CXX=${withval}/mpicxx
      elif test -f ${withval}/mpic++; then
	MPI_CXX=${withval}/mpic++
      else
        MPI_CXX=${withval}/mpiCC
      fi
      MPI_CC=${withval}/mpicc
      MPI_F77=${withval}/mpif77
    fi
  fi
]
)

AC_ARG_WITH(mpi,
[AC_HELP_STRING([--with-mpi=MPIROOT],[use MPI root directory (enables MPI)])],
[
  HAVE_PKG_MPI=yes
  MPI_DIR=${withval}
  AC_MSG_CHECKING(MPI directory)
  AC_MSG_RESULT([${MPI_DIR}])
]
)

#AC_ARG_WITH(mpi-include,
#[AC_HELP_STRING([--with-mpi-include],[Obsolete.  Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.])],
#[AC_MSG_ERROR([--with-mpi-include is an obsolte option.   Use --with-mpi-incdir=DIR instead.  Do not prefix DIR with '-I'.  For example '--with-mpi-incdir=/usr/lam_path/include'.])]
#)

AC_ARG_WITH(mpi-libs,
[AC_HELP_STRING([--with-mpi-libs="LIBS"],[MPI libraries @<:@"-lmpi"@:>@])],
[
  MPI_LIBS=${withval}
  AC_MSG_CHECKING(user-defined MPI libraries)
  AC_MSG_RESULT([${MPI_LIBS}])
]
)

AC_ARG_WITH(mpi-incdir,
[AC_HELP_STRING([--with-mpi-incdir=DIR],[MPI include directory @<:@MPIROOT/include@:>@  Do not use -I])],
[
  MPI_INC=${withval}
  AC_MSG_CHECKING(user-defined MPI includes)
  AC_MSG_RESULT([${MPI_INC}])
]
)

AC_ARG_WITH(mpi-libdir,
[AC_HELP_STRING([--with-mpi-libdir=DIR],[MPI library directory @<:@MPIROOT/lib@:>@  Do not use -L])],
[
  MPI_LIBDIR=${withval}
  AC_MSG_CHECKING(user-defined MPI library directory)
  AC_MSG_RESULT([${MPI_LIBDIR}])
]
)

AC_MSG_CHECKING(whether we are using MPI)
AC_MSG_RESULT([${HAVE_PKG_MPI}])

if test "X${HAVE_PKG_MPI}" = "Xyes"; then
   AC_DEFINE(HAVE_MPI,,[define if we want to use MPI])
fi

dnl Define Automake version of HAVE_MPI if appropriate

AM_CONDITIONAL(HAVE_MPI, [test "X${HAVE_PKG_MPI}" = "Xyes"])


dnl
dnl --------------------------------------------------------------------
dnl Check for MPI compilers (must be done *before* AC_PROG_CXX,
dnl AC_PROG_CC and AC_PROG_F77)
dnl 
dnl --------------------------------------------------------------------

if test -n "${MPI_CXX}"; then
  if test -f ${MPI_CXX}; then
    MPI_CXX_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CXX_EXISTS, ${MPI_CXX}, yes, no)
  fi

  if test "X${MPI_CXX_EXISTS}" = "Xyes"; then
#umennel: do not overwrite CXX, if already set
    test -n "$CXX" || CXX=${MPI_CXX};
  else
    echo "-----"
    echo "Cannot find MPI C++ compiler ${MPI_CXX}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a C++ compiler using CXX=<compiler>"
    echo "Do not use --with-mpi-compilers if using CXX=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI C++ compiler (${MPI_CXX}) not found.])
  fi
fi

if test -n "${MPI_CC}"; then
  if test -f ${MPI_CC}; then
    MPI_CC_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_CC_EXISTS, ${MPI_CC}, yes, no)
  fi

  if test "X${MPI_CC_EXISTS}" = "Xyes"; then
#umennel: do not overwrite CC, if already set
    test -n "$CC" || CC=${MPI_CC};
  else
    echo "-----"
    echo "Cannot find MPI C compiler ${MPI_CC}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a C compiler using CC=<compiler>"
    echo "Do not use --with-mpi-compilers if using CC=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI C compiler (${MPI_CC}) not found.])
  fi
fi

if test -n "${MPI_F77}"; then
  if test -f ${MPI_F77}; then
    MPI_F77_EXISTS=yes
  else
    AC_CHECK_PROG(MPI_F77_EXISTS, ${MPI_F77}, yes, no)
  fi

  if test "X${MPI_F77_EXISTS}" = "Xyes"; then
#umennel: do not overwrite F77, if already set
    test -n "$F77" || F77=${MPI_F77};
  else
    echo "-----"
    echo "Cannot find MPI Fortran compiler ${MPI_F77}."
    echo "Specify a path to all mpi compilers with --with-mpi-compilers=PATH"
    echo "or specify a Fortran 77 compiler using F77=<compiler>"
    echo "Do not use --with-mpi-compilers if using F77=<compiler>"
    echo "-----"
    AC_MSG_ERROR([MPI Fortran 77 compiler (${MPI_F77}) not found.])
  fi
fi
]) dnl TAC_ARG_CONFIG_MPI


dnl @synopsis TAC_ARG_WITH_FLAGS(lcase_name, UCASE_NAME)
dnl
dnl Test for --with-lcase_name="compiler/loader flags".  if defined, prepend 
dnl flags to standard UCASE_NAME definition.
dnl
dnl Use this macro to facilitate additional special flags that should be
dnl passed on to the preprocessor/compilers/loader.
dnl
dnl Example use
dnl 
dnl TAC_ARG_WITH_FLAGS(cxxflags, CXXFLAGS)
dnl 
dnl tests for --with-cxxflags and pre-pends to CXXFLAGS
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_FLAGS],
[
AC_MSG_CHECKING([whether additional [$2] flags should be added])
AC_ARG_WITH($1,
AC_HELP_STRING([--with-$1], 
[additional [$2] flags to be added: will prepend to [$2]]),
[
$2="${withval} ${$2}"
AC_MSG_RESULT([$2 = ${$2}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis TAC_ARG_WITH_LIBS
dnl
dnl Test for --with-libs="name(s)".
dnl 
dnl Prepends the specified name(s) to the list of libraries to link 
dnl with.  
dnl
dnl Example use
dnl
dnl TAC_ARG_WITH_LIBS
dnl 
dnl tests for --with-libs and pre-pends to LIBS
dnl
dnl @author Jim Willenbring <jmwille@sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LIBS],
[
AC_MSG_CHECKING([whether additional libraries are needed])
AC_ARG_WITH(libs,
AC_HELP_STRING([--with-libs], 
[List additional libraries here.  For example, --with-libs=-lsuperlu
or --with-libs=/path/libsuperlu.a]),
[
LIBS="${withval} ${LIBS}"
AC_MSG_RESULT([LIBS = ${LIBS}])
],
AC_MSG_RESULT(no)
)
]
)

dnl @synopsis TAC_ARG_WITH_AR
dnl
dnl Test for --with-ar="ar_program ar_flags".
dnl Default is "ar cru"
dnl 
dnl Generates an Automake conditional USE_ALTERNATE_AR that can be tested.  
dnl Generates the user-specified archiver command in @ALTERNATE_AR@.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_AR],
[
AC_ARG_WITH(ar,
AC_HELP_STRING([--with-ar], [override archiver command (default is "ar cru")]),
[
AC_MSG_CHECKING(user-defined archiver)
AC_MSG_RESULT([${withval}])
USE_ALTERNATE_AR=yes
ALTERNATE_AR="${withval}"
]
)

if test -n "${SPECIAL_AR}" && test "X${USE_ALTERNATE_AR}" != "Xyes";
then
  USE_ALTERNATE_AR=yes
  ALTERNATE_AR="${SPECIAL_AR}"
fi

AC_MSG_CHECKING(for special archiver command)
if test "X${USE_ALTERNATE_AR}" = "Xyes"; then
   AC_MSG_RESULT([${ALTERNATE_AR}])
   AM_CONDITIONAL(USE_ALTERNATE_AR, true)
else
   AC_MSG_RESULT([none])
   AM_CONDITIONAL(USE_ALTERNATE_AR, false)
fi
AC_SUBST(ALTERNATE_AR)
])


dnl @synopsis TAC_ARG_CHECK_MPI
dnl
dnl Check to make sure any definitions set in TAC_ARG_CONFIG_MPI
dnl are valid, set the MPI flags.  Test MPI compile using C++ compiler.
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_CHECK_MPI],
[

if test "X${HAVE_PKG_MPI}" = "Xyes"; then

  if test -n "${MPI_DIR}" && test -z "${MPI_INC}"; then
    MPI_INC="${MPI_DIR}/include"
  fi

  if test -n "${MPI_INC}"; then
    CPPFLAGS="${CPPFLAGS} -I${MPI_INC}"
  fi

  AC_LANG_CPLUSPLUS 
  AC_MSG_CHECKING(for mpi.h)
  AC_TRY_CPP([#include "mpi.h"],
    [AC_MSG_RESULT(yes)], 
    [
     AC_MSG_RESULT(no)  
     echo "-----"
     echo "Cannot link simple MPI program."
     echo "Try --with-mpi-compilers to specify MPI compilers."
     echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
     echo "to specify all the specific MPI compile options."
     echo "-----"
     AC_MSG_ERROR(MPI cannot link)
    ])

  if test -n "${MPI_DIR}" && test -z "${MPI_LIBDIR}"; then
    MPI_LIBDIR="${MPI_DIR}/lib"
  fi

  if test -n "${MPI_LIBDIR}"; then
    LDFLAGS="${LDFLAGS} -L${MPI_LIBDIR}"
  fi

  if test -z "${MPI_LIBS}" && test -n "${MPI_LIBDIR}"; then
    MPI_LIBS="-lmpi"
  fi

  if test -n "${MPI_LIBS}"; then
    LIBS="${MPI_LIBS} ${LIBS}"
  fi

#   AC_LANG_CPLUSPLUS 
#   AC_MSG_CHECKING(whether MPI will link using C++ compiler)
#   AC_TRY_LINK([#include <mpi.h>],
#   [int c; char** v; MPI_Init(&c,&v);],
#   [AC_MSG_RESULT(yes)], 
#   [AC_MSG_RESULT(no)  
#    echo "-----"
#    echo "Cannot link simple MPI program."
#    echo "Try --with-mpi-cxx to specify MPI C++ compile script."
#    echo "Or try --with-mpi-libs, --with-mpi-incdir, --with-mpi-libdir"
#    echo "to specify all the specific MPI compile options."
#    echo "-----"
#    AC_MSG_ERROR(MPI cannot link)]
#   )

fi
])

dnl @synopsis AC_CXX_NAMESPACES
dnl
dnl If the compiler can prevent names clashes using namespaces, define
dnl HAVE_NAMESPACES.
dnl
dnl @version $Id: ac_cxx_namespaces.m4,v 1.1 2003/01/08 19:18:58 jmwille Exp $
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ac_cv_cxx_namespaces,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([namespace Outer { namespace Inner { int i = 0; }}],
                [using namespace Outer::Inner; return i;],
 ac_cv_cxx_namespaces=yes, ac_cv_cxx_namespaces=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])


dnl @synopsis AC_CXX_NAMESPACES
dnl
dnl If the compiler can prevent names clashes using namespaces, define
dnl HAVE_NAMESPACES.
dnl
dnl @version $Id: ac_cxx_namespaces.m4,v 1.1 2003/01/08 19:18:58 jmwille Exp $
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ac_cv_cxx_namespaces,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([namespace Outer { namespace Inner { int i = 0; }}],
                [using namespace Outer::Inner; return i;],
 ac_cv_cxx_namespaces=yes, ac_cv_cxx_namespaces=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])

dnl @synopsis AC_CXX_HAVE_STL
dnl
dnl If the compiler supports the Standard Template Library, define HAVE_STL.
dnl
dnl @version $Id: ac_cxx_have_stl.m4,v 1.1 2003/01/08 19:18:58 jmwille Exp $
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_HAVE_STL],
[AC_CACHE_CHECK(whether the compiler supports Standard Template Library,
ac_cv_cxx_have_stl,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <list>
#include <deque>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],[list<int> x; x.push_back(5);
list<int>::iterator iter = x.begin(); if (iter != x.end()) ++iter; return 0;],
 ac_cv_cxx_have_stl=yes, ac_cv_cxx_have_stl=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_stl" = yes; then
  AC_DEFINE(HAVE_STL,,[define if the compiler supports Standard Template Library])
fi
])

dnl @synopsis AC_CXX_BOOL
dnl
dnl If the compiler recognizes bool as a separate built-in type,
dnl define HAVE_BOOL. Note that a typedef is not a separate
dnl type since you cannot overload a function such that it accepts either
dnl the basic type or the typedef.
dnl
dnl @version $Id: ac_cxx_bool.m4,v 1.1 2003/01/08 19:18:58 jmwille Exp $
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_BOOL],
[AC_CACHE_CHECK(whether the compiler recognizes bool as a built-in type,
ac_cv_cxx_bool,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
int f(int  x){return 1;}
int f(char x){return 1;}
int f(bool x){return 1;}
],[bool b = true; return f(b);],
 ac_cv_cxx_bool=yes, ac_cv_cxx_bool=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_bool" = yes; then
  AC_DEFINE(HAVE_BOOL,,[define if bool is a built-in type])
fi
])

dnl @synopsis AC_CXX_MUTABLE
dnl
dnl If the compiler allows modifying class data members flagged with
dnl the mutable keyword even in const objects (for example in the
dnl body of a const member function), define HAVE_MUTABLE.
dnl
dnl @version $Id: ac_cxx_mutable.m4,v 1.1 2003/01/08 19:18:58 jmwille Exp $
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_MUTABLE],
[AC_CACHE_CHECK(whether the compiler supports the mutable keyword,
ac_cv_cxx_mutable,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
class A { mutable int i;
          public:
          int f (int n) const { i = n; return i; }
        };
],[A a; return a.f (1);],
 ac_cv_cxx_mutable=yes, ac_cv_cxx_mutable=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_mutable" = yes; then
  AC_DEFINE(HAVE_MUTABLE,,[define if the compiler supports the mutable keyword])
fi
])

dnl @synopsis AC_CXX_NEW_FOR_SCOPING
dnl
dnl If the compiler accepts the new for scoping rules (the scope of a
dnl variable declared inside the parentheses is restricted to the
dnl for-body), define HAVE_NEW_FOR_SCOPING.
dnl
dnl @version $Id: ac_cxx_new_for_scoping.m4,v 1.1 2003/01/08 19:18:58 jmwille Exp $
dnl @author Luc Maisonobe
dnl
AC_DEFUN([AC_CXX_NEW_FOR_SCOPING],
[AC_CACHE_CHECK(whether the compiler accepts the new for scoping rules,
ac_cv_cxx_new_for_scoping,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE(,[
  int z = 0;
  for (int i = 0; i < 10; ++i)
    z = z + i;
  for (int i = 0; i < 10; ++i)
    z = z - i;
  return z;],
 ac_cv_cxx_new_for_scoping=yes, ac_cv_cxx_new_for_scoping=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_new_for_scoping" = yes; then
  AC_DEFINE(HAVE_NEW_FOR_SCOPING,,[define if the compiler accepts the new for scoping rules])
fi
])

dnl @synopsis AC_CXX_STD_SPRINTF
dnl
dnl If the compiler recognizes std::sprintf as a function for IO,
dnl define HAVE_STD_SPRINTF.  If this test fails, use sprintf with no std prefix
dnl Note that we try to compile two versions of this routine, one using cstdio and
dnl another using stdio.h.  This approach is used to eliminate the need to test which
dnl of the two header files is present.  If one or both is usable the test will return true.
dnl

AC_DEFUN([AC_CXX_STD_SPRINTF],
[AC_CACHE_CHECK([[whether the compiler recognizes std::sprintf as supported IO function]],
ac_cv_cxx_std_sprintf,
[ AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
#include <cstdio>
#include <string>
],
[
     int x = 100;
     char s[100];
     std::sprintf(s, "%d", x);
],
 ac_cv_cxx_std_sprintf1=yes, ac_cv_cxx_std_sprintf1=no)

AC_TRY_COMPILE([
#include <stdio.h>
#include <string>
],
[
     int x = 100;
     char s[100];
     std::sprintf(s, "%d", x);
],
 ac_cv_cxx_std_sprintf2=yes, ac_cv_cxx_std_sprintf2=no)

if (test "$ac_cv_cxx_std_sprintf1" = yes || test "$ac_cv_cxx_std_sprintf2" = yes); then
 ac_cv_cxx_std_sprintf=yes
else
 ac_cv_cxx_std_sprintf=no
fi
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_std_sprintf" = yes; then
  AC_DEFINE(HAVE_STD_SPRINTF,,[define if std::sprintf is supported])
fi
])

dnl @synopsis TAC_ARG_WITH_LIBDIRS
dnl
dnl Test for --with-libdirs="-Llibdir1 -Llibdir2".  if defined, 
dnl prepend "-Llibdir1 -Llibdir2" to LDFLAGS
dnl
dnl Use this macro to facilitate addition of directories to library search path.
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LIBDIRS],
[
AC_MSG_CHECKING([whether additional library search paths defined])
AC_ARG_WITH(libdirs,
AC_HELP_STRING([--with-libdirs], 
[OBSOLETE use --with-ldflags instead.  (ex. --with-ldflags="-L<DIR> -L<DIR2>")]),
[
LDFLAGS="${withval} ${LDFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis TAC_ARG_WITH_INCDIRS
dnl
dnl Test for --with-incdirs="-Iincdir1 -Iincdir2".  if defined, prepend 
dnl "-Iincdir1 -Iincdir2" to CPPFLAGS
dnl
dnl Use this macro to facilitate addition of directories to include file search path.
dnl 
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_INCDIRS],
[
AC_MSG_CHECKING([whether additional include search paths defined])
AC_ARG_WITH(incdirs,
AC_HELP_STRING([--with-incdirs], 
[additional directories containing include files: will prepend to search here for includes, use -Idir format]),
[
CPPFLAGS="${withval} ${CPPFLAGS}"
AC_MSG_RESULT([${withval}])
],
AC_MSG_RESULT(no)
)
])


dnl @synopsis TAC_ARG_WITH_BLASLIB
dnl
dnl Test for --with-blaslib="name".
dnl 
dnl Prepends the specified name to the list of files to check for BLAS
dnl routines.  
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_BLASLIB],
[
AC_ARG_WITH(blaslib,
AC_HELP_STRING([--with-blaslib], 
[name of library containing BLAS: will search lib directories for
-lname]),
[
USE_BLASLIB=yes
NEWBLASLIB=${withval}
]
)

   BLASLIBS="cxml blas complib.sgimath"

if test "X${USE_BLASLIB}" = "Xyes"; then

   BLASLIBS="${NEWBLASLIB} ${BLASLIBS}"

fi
])

dnl @synopsis TAC_ARG_WITH_LAPACKLIB
dnl
dnl Test for --with-lapacklib="name".
dnl 
dnl Prepends the specified name to the list of files to check for LAPACK
dnl routines.  
dnl
dnl @author Mike Heroux <mheroux@cs.sandia.gov>
dnl
AC_DEFUN([TAC_ARG_WITH_LAPACKLIB],
[
AC_ARG_WITH(lapacklib,
AC_HELP_STRING([--with-lapacklib], 
[name of library containing LAPACK: will search lib directories for -lname]),
[
USE_LAPACKLIB=yes
NEWLAPACKLIB=${withval}
]
)

   LAPACKLIBS="cxml lapack complib.sgimath"

if test "X${USE_LAPACKLIB}" = "Xyes"; then

   LAPACKLIBS="${NEWLAPACKLIB} ${LAPACKLIBS}"
fi
])


dnl @synopsis ACX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the BLAS
dnl linear-algebra interface (see http://www.netlib.org/blas/).
dnl On success, it sets the BLAS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with BLAS, you should link with:
dnl
dnl 	$BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl Many libraries are searched for, from ATLAS to CXML to ESSL.
dnl The user may also use --with-blas=<lib> in order to use some
dnl specific BLAS library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the BLAS library.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a BLAS
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_BLAS.
dnl
dnl This macro requires autoconf 2.50 or later.
dnl
dnl @version $Id: acx_blas.m4,v 1.2 2004/05/26 21:26:23 jmwille Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl
dnl Edited by Jim Willenbring on 5-14-2004 to check for dgemm instead of
dnl sgemm.
AC_DEFUN([ACX_BLAS], [
AC_PREREQ(2.50)
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
AC_F77_FUNC(sgemm)
AC_F77_FUNC(dgemm)

acx_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $acx_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $dgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC($dgemm, [acx_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($acx_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $acx_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC($dgemm, [acx_blas_ok=yes])
	LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $dgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[acx_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $dgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Alpha CXML library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(cxml, $dgemm, [acx_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(dxml, $dgemm, [acx_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $acx_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $dgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 acx_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(scs, $dgemm, [acx_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $dgemm,
		     [acx_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $dgemm,
		[AC_CHECK_LIB(essl, $dgemm,
			[acx_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $acx_blas_ok = no; then
	AC_CHECK_LIB(blas, $dgemm, [acx_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$acx_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        acx_blas_ok=no
        $2
fi
])dnl ACX_BLAS


dnl @synopsis ACX_LAPACK([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for a library that implements the LAPACK
dnl linear-algebra interface (see http://www.netlib.org/lapack/).
dnl On success, it sets the LAPACK_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with LAPACK, you should link with:
dnl
dnl     $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order.  BLAS_LIBS is the output variable of the ACX_BLAS
dnl macro, called automatically.  FLIBS is the output variable of the
dnl AC_F77_LIBRARY_LDFLAGS macro (called if necessary by ACX_BLAS),
dnl and is sometimes necessary in order to link with F77 libraries.
dnl Users will also need to use AC_F77_DUMMY_MAIN (see the autoconf
dnl manual), for the same reason.
dnl
dnl The user may also use --with-lapack=<lib> in order to use some
dnl specific LAPACK library <lib>.  In order to link successfully,
dnl however, be aware that you will probably need to use the same
dnl Fortran compiler (which can be set via the F77 env. var.) as
dnl was used to compile the LAPACK and BLAS libraries.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if a LAPACK
dnl library is found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if it is not found.  If ACTION-IF-FOUND is not specified,
dnl the default action will define HAVE_LAPACK.
dnl
dnl @version $Id: acx_lapack.m4,v 1.2 2004/05/14 18:00:11 jmwille Exp $
dnl @author Steven G. Johnson <stevenj@alum.mit.edu>
dnl edited by Jim Willenbring <jmwille@sandia.gov> to check for sgecon
dnl rather than cheev because by default (as of 8-13-2002) Trilinos
dnl does not build the complex portions of the lapack library.  Edited
dnl again on 5-13-2004 to check for dgecon instead of sgecon.

AC_DEFUN([ACX_LAPACK], [
AC_REQUIRE([ACX_BLAS])
acx_lapack_ok=no

AC_ARG_WITH(lapack,
        [AC_HELP_STRING([--with-lapack=<lib>], [use LAPACK library <lib>])])
case $with_lapack in
        yes | "") ;;
        no) acx_lapack_ok=disable ;;
        -* | */* | *.a | *.so | *.so.* | *.o) LAPACK_LIBS="$with_lapack" ;;
        *) LAPACK_LIBS="-l$with_lapack" ;;
esac

# Get fortran linker name of LAPACK function to check for.
AC_F77_FUNC(dgecon)

# We cannot use LAPACK if BLAS is not found
if test "x$acx_blas_ok" != xyes; then
        acx_lapack_ok=noblas
fi

# First, check LAPACK_LIBS environment variable
if test "x$LAPACK_LIBS" != x; then
        save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
        AC_MSG_CHECKING([for $dgecon in $LAPACK_LIBS])
        AC_TRY_LINK_FUNC($dgecon, [acx_lapack_ok=yes], [LAPACK_LIBS=""])
        AC_MSG_RESULT($acx_lapack_ok)
        LIBS="$save_LIBS"
        if test acx_lapack_ok = no; then
                LAPACK_LIBS=""
        fi
fi

# LAPACK linked to by default?  (is sometimes included in BLAS lib)
if test $acx_lapack_ok = no; then
        save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS $FLIBS"
        AC_CHECK_FUNC($dgecon, [acx_lapack_ok=yes])
        LIBS="$save_LIBS"
fi

# Generic LAPACK library?
for lapack in lapack lapack_rs6k; do
        if test $acx_lapack_ok = no; then
                save_LIBS="$LIBS"; LIBS="$LIBS $FLIBS"
                AC_CHECK_LIB($lapack, $dgecon,
                    [acx_lapack_ok=yes; LAPACK_LIBS="-l$lapack"], [], [$BLAS_LIBS])
                LIBS="$save_LIBS"
        fi
done

AC_SUBST(LAPACK_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_lapack_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_LAPACK,1,[Define if you have LAPACK library.]),[$1])
        :
else
        acx_lapack_ok=no
        $2
fi
])dnl ACX_LAPACK




dnl @synopsis ACX_TRILINOS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for libraries that implement the Trilinos ML
dnl and AztecOO packages.
dnl On success, it sets the TRILINOS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with Trilinos, you should link with:
dnl
dnl     $TRILINOS_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. LAPACK_LIBS is the output variable of the ACX_LAPACK
dnl macro, called automatically. BLAS_LIBS is the output variable of the 
dnl ACX_BLAS macro, called automatically. FLIBS is the output variable
dnl of the AC_F77_LIBRARY_LDFLAGS macro (called if necessary by
dnl ACX_BLAS), and is sometimes necessary in order to link with F77
dnl libraries. Users will also need to use AC_F77_DUMMY_MAIN (see the
dnl autoconf manual), for the same reason.
dnl
dnl The user may also use --with-trilinos=DIR in order to specifiy 
dnl the Trilinos install directory.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if the necessary
dnl libraries are found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if they are not found.
dnl
dnl @author Uche Mennel <umennel@student.ethz.ch>
dnl
AC_DEFUN([ACX_TRILINOS], [
AC_REQUIRE([ACX_LAPACK])
AC_REQUIRE([AC_F77_LIBRARY_LDFLAGS])
acx_trilinos_ok=yes

acx_trilinos_save_LIBS="$LIBS"
TRILINOS_LIBS="-laztecoo -ltriutils -lepetraext -lepetra -lteuchos"
EXTRA_LIBS="$PARMETIS_LIBS $SUPERLU_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS $LIBS"


AC_ARG_WITH(trilinos,
AC_HELP_STRING([--with-trilinos], 
[Specify Trilinos directory.  For example, --with-trilinos=/opt/trilinos]),
[
CPPFLAGS="$CPPFLAGS -I$withval/include"
LDFLAGS="$LDFLAGS -L$withval/lib"
],)

# Find Trilinos libraries
AC_MSG_CHECKING(for AztecOO in Trilinos)
LIBS="$TRILINOS_LIBS $EXTRA_LIBS"
AC_TRY_LINK([
class AztecOO;
],
[AztecOO solver();],
[echo "yes"], 
[acx_trilinos_ok=no; echo "no"])


AC_MSG_CHECKING(wether Ifpack was enabled in Trilinos)
LIBS="-lifpack $TRILINOS_LIBS $LIBS"
AC_TRY_LINK([
class Ifpack;
], 
[Ifpack factory();],
[AC_DEFINE(HAVE_IFPACKLIB, [], [ifpack]) TRILINOS_LIBS="-lifpack $TRILINOS_LIBS"; echo "yes"], 
[echo "no"])


AC_MSG_CHECKING(wether Amesos was enabled in Trilinos)
LIBS="-lamesos $TRILINOS_LIBS $EXTRA_LIBS"
AC_TRY_LINK([
class Amesos;
], 
[Amesos factory();],
[AC_DEFINE(HAVE_AMESOSLIB, [], [amesos]) TRILINOS_LIBS="-lamesos $TRILINOS_LIBS"; echo "yes"], 
[echo "no"])


AC_MSG_CHECKING(wether ML was enabled in Trilinos)
LIBS="-lml $TRILINOS_LIBS $EXTRA_LIBS"
AC_TRY_LINK([
namespace ML_Epetra {
class MultiLevelPreconditioner;
}
class Epetra_RowMatrix;
], 
[Epetra_RowMatrix *A; 
ML_Epetra::MultiLevelPreconditioner Prec();],
[AC_DEFINE(HAVE_MLLIB, [], [ml]) TRILINOS_LIBS="-lml $TRILINOS_LIBS"; echo "yes"], 
[echo "no"])

LIBS="$acx_trilinos_save_LIBS"

AC_SUBST(TRILINOS_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_trilinos_ok" = xyes; then
        $1
        :
else
        acx_trilinos_ok=no
        $2
fi
]) dnl ACX_TRILINOS


dnl @synopsis ACX_PARMETIS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for libraries that implement the ParMETIS
dnl graph partitioning routines.
dnl On success, it sets the PARMETIS_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with ParMETIS, you should link with:
dnl
dnl     $PARMETIS_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. LAPACK_LIBS is the output variable of the ACX_LAPACK
dnl macro, called automatically. BLAS_LIBS is the output variable of the 
dnl ACX_BLAS macro, called automatically. FLIBS is the output variable
dnl of the AC_F77_LIBRARY_LDFLAGS macro (called if necessary by
dnl ACX_BLAS), and is sometimes necessary in order to link with F77
dnl libraries. Users will also need to use AC_F77_DUMMY_MAIN (see the
dnl autoconf manual), for the same reason.
dnl
dnl The user may also use --with-parmetis=DIR in order to specifiy 
dnl the ParMETIS install directory.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if the necessary
dnl libraries are found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if they are not found.
dnl
dnl @author Uche Mennel <umennel@student.ethz.ch>
dnl
AC_DEFUN([ACX_PARMETIS], [
acx_parmetis_ok=no

acx_parmetis_save_LIBS="$LIBS"
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
PARMETIS_LIBS="-lmetis"

AC_ARG_WITH(parmetis,
AC_HELP_STRING([--with-parmetis], 
[Specify ParMETIS directory.  For example, --with-parmetis=/opt/ParMETIS]),
[
CPPFLAGS="$CPPFLAGS -I$withval"
LDFLAGS="$LDFLAGS -L$withval"
],)

AC_CHECK_LIB(parmetis,ParMETIS_V3_PartKway,[acx_parmetis_ok=yes; PARMETIS_LIBS="-lparmetis $PARMETIS_LIBS"],[],[$PARMETIS_LIBS])

LIBS="$acx_parmetis_save_LIBS"

AC_SUBST(PARMETIS_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_parmetis_ok" = xyes; then
        $1
        :
else
        acx_parmetis_ok=no
        $2
fi
]) dnl ACX_PARMETIS

dnl @synopsis ACX_HDF5([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for libraries that implement the Hdf5
dnl graph partitioning routines.
dnl On success, it sets the HDF5_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with Hdf5, you should link with:
dnl
dnl     $HDF5_LIBS $LIBS 
dnl
dnl in that order.
dnl
dnl The user may also use --with-hdf5=DIR in order to specifiy 
dnl the Hdf5 install directory.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if the necessary
dnl libraries are found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if they are not found.
dnl
dnl @author Uche Mennel <umennel@student.ethz.ch>
dnl
AC_DEFUN([ACX_HDF5], [
acx_hdf5_ok=no

acx_hdf5_save_LIBS="$LIBS"
HDF5_LIBS="-lz"

AC_ARG_WITH(hdf5,
AC_HELP_STRING([--with-hdf5], 
[Specify Hdf5 directory.  For example, --with-hdf5=/opt/Hdf5]),
[
CPPFLAGS="$CPPFLAGS -I$withval/include"
LDFLAGS="$LDFLAGS -L$withval/lib"
],)

AC_CHECK_LIB(hdf5,H5Fopen,[acx_hdf5_ok=yes; HDF5_LIBS="-lhdf5 $HDF5_LIBS"],[],[$HDF5_LIBS])

LIBS="$acx_hdf5_save_LIBS"

AC_SUBST(HDF5_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_hdf5_ok" = xyes; then
        $1
        :
else
        acx_hdf5_ok=no
        $2
fi
]) dnl ACX_HDF5


dnl @synopsis ACX_SUPERLU([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
dnl
dnl This macro looks for libraries that implement the SuperLU
dnl routines.
dnl On success, it sets the SUPERLU_LIBS output variable to
dnl hold the requisite library linkages.
dnl
dnl To link with SuperLU, you should link with:
dnl
dnl     $SUPERLU_LIBS $LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS
dnl
dnl in that order. LAPACK_LIBS is the output variable of the ACX_LAPACK
dnl macro, called automatically. BLAS_LIBS is the output variable of the 
dnl ACX_BLAS macro, called automatically. FLIBS is the output variable
dnl of the AC_F77_LIBRARY_LDFLAGS macro (called if necessary by
dnl ACX_BLAS), and is sometimes necessary in order to link with F77
dnl libraries. Users will also need to use AC_F77_DUMMY_MAIN (see the
dnl autoconf manual), for the same reason.
dnl
dnl The user may also use --with-superlu=DIR in order to specifiy 
dnl the SuperLU install directory.
dnl
dnl ACTION-IF-FOUND is a list of shell commands to run if the necessary
dnl libraries are found, and ACTION-IF-NOT-FOUND is a list of commands
dnl to run it if they are not found.
dnl
dnl @author Uche Mennel <umennel@student.ethz.ch>
dnl
AC_DEFUN([ACX_SUPERLU], [
acx_superlu_ok=no

acx_superlu_save_LIBS="$LIBS"
LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"

AC_ARG_WITH(superlu,
AC_HELP_STRING([--with-superlu], 
[Specify SuperLU directory.  For example, --with-SuperLU=/opt/SuperLU]),
[
CPPFLAGS="$CPPFLAGS -I$withval"
LDFLAGS="$LDFLAGS -L$withval"
],)

AC_CHECK_LIB(superlu,dgssv,[acx_superlu_ok=yes; SUPERLU_LIBS="-lsuperlu"],[],[])

LIBS="$acx_superlu_save_LIBS"

AC_SUBST(SUPERLU_LIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$acx_superlu_ok" = xyes; then
        $1
        :
else
        acx_superlu_ok=no
        $2
fi
]) dnl ACX_SUPERLU


# Copyright (C) 2004 Oren Ben-Kiki
# This file is distributed under the same terms as the Autoconf macro files.

# Generate automatic documentation using Doxygen. Works in concert with the
# aminclude.m4 file and a compatible doxygen configuration file. Defines the
# following public macros:
#
# DX_???_FEATURE(ON|OFF) - control the default setting fo a Doxygen feature.
# Supported features are 'DOXYGEN' itself, 'DOT' for generating graphics,
# 'HTML' for plain HTML, 'CHM' for compressed HTML help (for MS users), 'CHI'
# for generating a seperate .chi file by the .chm file, and 'MAN', 'RTF',
# 'XML', 'PDF' and 'PS' for the appropriate output formats. The environment
# variable DOXYGEN_PAPER_SIZE may be specified to override the default 'a4wide'
# paper size.
#
# By default, HTML, PDF and PS documentation is generated as this seems to be
# the most popular and portable combination. MAN pages created by Doxygen are
# usually problematic, though by picking an appropriate subset and doing some
# massaging they might be better than nothing. CHM and RTF are specific for MS
# (note that you can't generate both HTML and CHM at the same time). The XML is
# rather useless unless you apply specialized post-processing to it.
#
# The macro mainly controls the default state of the feature. The use can
# override the default by specifying --enable or --disable. The macros ensure
# that contradictory flags are not given (e.g., --enable-doxygen-html and
# --enable-doxygen-chm, --enable-doxygen-anything with --disable-doxygen, etc.)
# Finally, each feature will be automatically disabled (with a warning) if the
# required programs are missing.
#
# Once all the feature defaults have been specified, call DX_INIT_DOXYGEN with
# the following parameters: a one-word name for the project for use as a
# filename base etc., an optional configuration file name (the default is
# 'Doxyfile', the same as Doxygen's default), and an optional output directory
# name (the default is 'doc/doxygen').

## ----------##
## Defaults. ##
## ----------##

DX_ENV=""
AC_DEFUN([DX_FEATURE_doc],  ON)
AC_DEFUN([DX_FEATURE_dot],  ON)
AC_DEFUN([DX_FEATURE_man],  OFF)
AC_DEFUN([DX_FEATURE_html], ON)
AC_DEFUN([DX_FEATURE_chm],  OFF)
AC_DEFUN([DX_FEATURE_chi],  OFF)
AC_DEFUN([DX_FEATURE_rtf],  OFF)
AC_DEFUN([DX_FEATURE_xml],  OFF)
AC_DEFUN([DX_FEATURE_pdf],  ON)
AC_DEFUN([DX_FEATURE_ps],   ON)

## --------------- ##
## Private macros. ##
## --------------- ##

# DX_ENV_APPEND(VARIABLE, VALUE)
# ------------------------------
# Append VARIABLE="VALUE" to DX_ENV for invoking doxygen.
AC_DEFUN([DX_ENV_APPEND], [AC_SUBST([DX_ENV], ["$DX_ENV $1='$2'"])])

# DX_DIRNAME_EXPR
# ---------------
# Expand into a shell expression prints the directory part of a path.
AC_DEFUN([DX_DIRNAME_EXPR],
         [[expr ".$1" : '\(\.\)[^/]*$' \| "x$1" : 'x\(.*\)/[^/]*$']])

# DX_IF_FEATURE(FEATURE, IF-ON, IF-OFF)
# -------------------------------------
# Expands according to the M4 (static) status of the feature.
AC_DEFUN([DX_IF_FEATURE], [ifelse(DX_FEATURE_$1, ON, [$2], [$3])])

# DX_REQUIRE_PROG(VARIABLE, PROGRAM)
# ----------------------------------
# Require the specified program to be found for the DX_CURRENT_FEATURE to work.
AC_DEFUN([DX_REQUIRE_PROG], [
AC_PATH_TOOL([$1], [$2])
if test "$DX_FLAG_[]DX_CURRENT_FEATURE$$1" = 1; then
    AC_MSG_WARN([$2 not found - will not DX_CURRENT_DESCRIPTION])
    AC_SUBST([DX_FLAG_[]DX_CURRENT_FEATURE], 0)
fi
])

# DX_TEST_FEATURE(FEATURE)
# ------------------------
# Expand to a shell expression testing whether the feature is active.
AC_DEFUN([DX_TEST_FEATURE], [test "$DX_FLAG_$1" = 1])

# DX_CHECK_DEPEND(REQUIRED_FEATURE, REQUIRED_STATE)
# -------------------------------------------------
# Verify that a required features has the right state before trying to turn on
# the DX_CURRENT_FEATURE.
AC_DEFUN([DX_CHECK_DEPEND], [
test "$DX_FLAG_$1" = "$2" \
|| AC_MSG_ERROR([doxygen-DX_CURRENT_FEATURE ifelse([$2], 1,
                            requires, contradicts) doxygen-DX_CURRENT_FEATURE])
])

# DX_CLEAR_DEPEND(FEATURE, REQUIRED_FEATURE, REQUIRED_STATE)
# ----------------------------------------------------------
# Turn off the DX_CURRENT_FEATURE if the required feature is off.
AC_DEFUN([DX_CLEAR_DEPEND], [
test "$DX_FLAG_$1" = "$2" || AC_SUBST([DX_FLAG_[]DX_CURRENT_FEATURE], 0)
])

# DX_FEATURE_ARG(FEATURE, DESCRIPTION,
#                CHECK_DEPEND, CLEAR_DEPEND,
#                REQUIRE, DO-IF-ON, DO-IF-OFF)
# --------------------------------------------
# Parse the command-line option controlling a feature. CHECK_DEPEND is called
# if the user explicitly turns the feature on (and invokes DX_CHECK_DEPEND),
# otherwise CLEAR_DEPEND is called to turn off the default state if a required
# feature is disabled (using DX_CLEAR_DEPEND). REQUIRE performs additional
# requirement tests (DX_REQUIRE_PROG). Finally, an automake flag is set and
# DO-IF-ON or DO-IF-OFF are called according to the final state of the feature.
AC_DEFUN([DX_ARG_ABLE], [
    AC_DEFUN([DX_CURRENT_FEATURE], [$1])
    AC_DEFUN([DX_CURRENT_DESCRIPTION], [$2])
    AC_ARG_ENABLE(doxygen-$1,
                  [AS_HELP_STRING(DX_IF_FEATURE([$1], [--disable-doxygen-$1],
                                                      [--enable-doxygen-$1]),
                                  DX_IF_FEATURE([$1], [don't $2], [$2]))],
                  [
case "$enableval" in
#(
y|Y|yes|Yes|YES)
    AC_SUBST([DX_FLAG_$1], 1)
    $3
;; #(
n|N|no|No|NO)
    AC_SUBST([DX_FLAG_$1], 0)
;; #(
*)
    AC_MSG_ERROR([invalid value '$enableval' given to doxygen-$1])
;;
esac
], [
AC_SUBST([DX_FLAG_$1], [DX_IF_FEATURE([$1], 1, 0)])
$4
])
if DX_TEST_FEATURE([$1]); then
    $5
    :
fi
if DX_TEST_FEATURE([$1]); then
    AM_CONDITIONAL(DX_COND_$1, :)
    $6
    :
else
    AM_CONDITIONAL(DX_COND_$1, false)
    $7
    :
fi
])

## -------------- ##
## Public macros. ##
## -------------- ##

# DX_XXX_FEATURE(DEFAULT_STATE)
# -----------------------------
AC_DEFUN([DX_DOXYGEN_FEATURE], [AC_DEFUN([DX_FEATURE_doc],  [$1])])
AC_DEFUN([DX_MAN_FEATURE],     [AC_DEFUN([DX_FEATURE_man],  [$1])])
AC_DEFUN([DX_HTML_FEATURE],    [AC_DEFUN([DX_FEATURE_html], [$1])])
AC_DEFUN([DX_CHM_FEATURE],     [AC_DEFUN([DX_FEATURE_chm],  [$1])])
AC_DEFUN([DX_CHI_FEATURE],     [AC_DEFUN([DX_FEATURE_chi],  [$1])])
AC_DEFUN([DX_RTF_FEATURE],     [AC_DEFUN([DX_FEATURE_rtf],  [$1])])
AC_DEFUN([DX_XML_FEATURE],     [AC_DEFUN([DX_FEATURE_xml],  [$1])])
AC_DEFUN([DX_XML_FEATURE],     [AC_DEFUN([DX_FEATURE_xml],  [$1])])
AC_DEFUN([DX_PDF_FEATURE],     [AC_DEFUN([DX_FEATURE_pdf],  [$1])])
AC_DEFUN([DX_PS_FEATURE],      [AC_DEFUN([DX_FEATURE_ps],   [$1])])

# DX_INIT_DOXYGEN(PROJECT, [CONFIG-FILE], [OUTPUT-DOC-DIR])
# ---------------------------------------------------------
# PROJECT also serves as the base name for the documentation files.
# The default CONFIG-FILE is "Doxyfile" and OUTPUT-DOC-DIR is "doc/doxygen".
AC_DEFUN([DX_INIT_DOXYGEN], [

# Files:
AC_SUBST([DX_PROJECT], [$1])
AC_SUBST([DX_CONFIG], [ifelse([$2], [], Doxyfile, [$2])])
AC_SUBST([DX_DOCDIR], [ifelse([$3], [], doc/doxygen, [$3])])

# Environment variables used inside doxygen.cfg:
DX_ENV_APPEND(SRCDIR, $srcdir)
DX_ENV_APPEND(PROJECT, $DX_PROJECT)
DX_ENV_APPEND(DOCDIR, $DX_DOCDIR)
DX_ENV_APPEND(VERSION, $PACKAGE_VERSION)

# Doxygen itself:
DX_ARG_ABLE(doc, [generate any doxygen documentation],
            [],
            [],
            [DX_REQUIRE_PROG([DX_DOXYGEN], doxygen)
             DX_REQUIRE_PROG([DX_PERL], perl)],
            [DX_ENV_APPEND(PERL_PATH, $DX_PERL)])

# Dot for graphics:
DX_ARG_ABLE(dot, [generate graphics for doxygen documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [DX_REQUIRE_PROG([DX_DOT], dot)],
            [DX_ENV_APPEND(HAVE_DOT, YES)
             DX_ENV_APPEND(DOT_PATH, [`DX_DIRNAME_EXPR($DX_DOT)`])],
            [DX_ENV_APPEND(HAVE_DOT, NO)])

# Man pages generation:
DX_ARG_ABLE(man, [generate doxygen manual pages],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [],
            [DX_ENV_APPEND(GENERATE_MAN, YES)],
            [DX_ENV_APPEND(GENERATE_MAN, NO)])

# RTF file generation:
DX_ARG_ABLE(rtf, [generate doxygen RTF documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [],
            [DX_ENV_APPEND(GENERATE_RTF, YES)],
            [DX_ENV_APPEND(GENERATE_RTF, NO)])

# XML file generation:
DX_ARG_ABLE(xml, [generate doxygen XML documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [],
            [DX_ENV_APPEND(GENERATE_XML, YES)],
            [DX_ENV_APPEND(GENERATE_XML, NO)])

# (Compressed) HTML help generation:
DX_ARG_ABLE(chm, [generate doxygen compressed HTML help documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [DX_REQUIRE_PROG([DX_HHC], hhc)],
            [DX_ENV_APPEND(HHC_PATH, $DX_HHC)
             DX_ENV_APPEND(GENERATE_HTML, YES)
             DX_ENV_APPEND(GENERATE_HTMLHELP, YES)],
            [DX_ENV_APPEND(GENERATE_HTMLHELP, NO)])

# Seperate CHI file generation.
DX_ARG_ABLE(chi, [generate doxygen seperate compressed HTML help index file],
            [DX_CHECK_DEPEND(chm, 1)],
            [DX_CLEAR_DEPEND(chm, 1)],
            [],
            [DX_ENV_APPEND(GENERATE_CHI, YES)],
            [DX_ENV_APPEND(GENERATE_CHI, NO)])

# Plain HTML pages generation:
DX_ARG_ABLE(html, [generate doxygen plain HTML documentation],
            [DX_CHECK_DEPEND(doc, 1) DX_CHECK_DEPEND(chm, 0)],
            [DX_CLEAR_DEPEND(doc, 1) DX_CLEAR_DEPEND(chm, 0)],
            [],
            [DX_ENV_APPEND(GENERATE_HTML, YES)],
            [DX_TEST_FEATURE(chm) || DX_ENV_APPEND(GENERATE_HTML, NO)])

# PostScript file generation:
DX_ARG_ABLE(ps, [generate doxygen PostScript documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [DX_REQUIRE_PROG([DX_LATEX], latex)
             DX_REQUIRE_PROG([DX_MAKEINDEX], makeindex)
             DX_REQUIRE_PROG([DX_DVIPS], dvips)
             DX_REQUIRE_PROG([DX_EGREP], egrep)])

# PDF file generation:
DX_ARG_ABLE(pdf, [generate doxygen PDF documentation],
            [DX_CHECK_DEPEND(doc, 1)],
            [DX_CLEAR_DEPEND(doc, 1)],
            [DX_REQUIRE_PROG([DX_PDFLATEX], pdflatex)
             DX_REQUIRE_PROG([DX_MAKEINDEX], makeindex)
             DX_REQUIRE_PROG([DX_EGREP], egrep)])

# LaTeX generation for PS and/or PDF:
if DX_TEST_FEATURE(ps) || DX_TEST_FEATURE(pdf); then
    AM_CONDITIONAL(DX_COND_latex, :)
    DX_ENV_APPEND(GENERATE_LATEX, YES)
else
    AM_CONDITIONAL(DX_COND_latex, false)
    DX_ENV_APPEND(GENERATE_LATEX, NO)
fi

# Paper size for PS and/or PDF:
AC_ARG_VAR(DOXYGEN_PAPER_SIZE,
           [a4wide (default), a4, letter, legal or executive])
case "$DOXYGEN_PAPER_SIZE" in
#(
"")
    AC_SUBST(DOXYGEN_PAPER_SIZE, "")
;; #(
a4wide|a4|letter|legal|executive)
    DX_ENV_APPEND(PAPER_SIZE, $DOXYGEN_PAPER_SIZE)
;; #(
*)
    AC_MSG_ERROR([unknown DOXYGEN_PAPER_SIZE='$DOXYGEN_PAPER_SIZE'])
;;
esac

#For debugging:
#echo DX_FLAG_doc=$DX_FLAG_doc
#echo DX_FLAG_dot=$DX_FLAG_dot
#echo DX_FLAG_man=$DX_FLAG_man
#echo DX_FLAG_html=$DX_FLAG_html
#echo DX_FLAG_chm=$DX_FLAG_chm
#echo DX_FLAG_chi=$DX_FLAG_chi
#echo DX_FLAG_rtf=$DX_FLAG_rtf
#echo DX_FLAG_xml=$DX_FLAG_xml
#echo DX_FLAG_pdf=$DX_FLAG_pdf
#echo DX_FLAG_ps=$DX_FLAG_ps
#echo DX_ENV=$DX_ENV
])
