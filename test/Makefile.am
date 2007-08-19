EXEEXT = .exe

noinst_PROGRAMS = test_preconditioner

test_preconditioner_SOURCES = $(srcdir)/test_preconditioner.cpp

test_preconditioner_DEPENDENCIES = \
	$(top_builddir)/src/libparfe.a 

test_preconditioner_CXXFLAGS = -I$(top_srcdir)/src

test_preconditioner_LDADD = $(top_builddir)/src/libparfe.a \
	$(TRILINOS_LIBS) $(PARMETIS_LIBS) $(SUPERLU_LIBS) $(LAPACK_LIBS) $(BLAS_LIBS) $(HDF5_LIBS) $(FLIBS)

SUBDIRS = 
