EXEEXT = .exe

noinst_PROGRAMS = pfaim new_pfaim converter analyze_all

pfaim_SOURCES      = $(srcdir)/pfaim.cpp
pfaim_DEPENDENCIES = $(top_builddir)/src/libparfe.a 
pfaim_CXXFLAGS     = -I$(top_srcdir)/src
pfaim_LDADD        = $(top_builddir)/src/libparfe.a \
		     $(TRILINOS_LIBS) $(PARMETIS_LIBS) \
		     $(LAPACK_LIBS) $(BLAS_LIBS) $(HDF5_LIBS) $(FLIBS)

new_pfaim_SOURCES      = $(srcdir)/new_pfaim.cpp
new_pfaim_DEPENDENCIES = $(top_builddir)/src/libparfe.a 
new_pfaim_CXXFLAGS     = -I$(top_srcdir)/src
new_pfaim_LDADD        = $(top_builddir)/src/libparfe.a \
		         $(TRILINOS_LIBS) $(PARMETIS_LIBS) \
		         $(LAPACK_LIBS) $(BLAS_LIBS) $(HDF5_LIBS) $(FLIBS)

converter_SOURCES      = $(srcdir)/converter.cpp
converter_DEPENDENCIES = $(top_builddir)/src/libparfe.a 
converter_CXXFLAGS     = -I$(top_srcdir)/src
converter_LDADD        = $(top_builddir)/src/libparfe.a \
		         $(TRILINOS_LIBS) $(PARMETIS_LIBS) \
		         $(LAPACK_LIBS) $(BLAS_LIBS) $(HDF5_LIBS) $(FLIBS)

analyze_all_SOURCES      = $(srcdir)/analyze_all.cpp
analyze_all_DEPENDENCIES = $(top_builddir)/src/libparfe.a 
analyze_all_CXXFLAGS     = -I$(top_srcdir)/src
analyze_all_LDADD        = $(top_builddir)/src/libparfe.a \
		          $(TRILINOS_LIBS) $(PARMETIS_LIBS) \
		          $(LAPACK_LIBS) $(BLAS_LIBS) $(HDF5_LIBS) $(FLIBS)

SUBDIRS = 
