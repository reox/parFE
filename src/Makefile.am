CORE_H = \
	$(srcdir)/parfe_ConfigDefs.h \
	$(srcdir)/BoundaryCondition.h \
	$(srcdir)/ElasticityProblem.h \
	$(srcdir)/VbrElasticityProblem.h \
	$(srcdir)/EbeElasticityProblem.h \
	$(srcdir)/Jacobi.h \
	$(srcdir)/DistMesh.h \
	$(srcdir)/ElementByElementMatrix.h \
	$(srcdir)/fem.h \
	$(srcdir)/FEParameterList.hpp \
	$(srcdir)/MEDIT_MeshWriter.hpp \
	$(srcdir)/PMVIS_MeshWriter.hpp \
	$(srcdir)/MEDIT_DispWriter.hpp \
	$(srcdir)/IBT_ProblemReader.hpp \
	$(srcdir)/IBT_ProblemWriter.hpp \
	$(srcdir)/IBT_DispWriter.hpp \
	$(srcdir)/IBT_SolutionWriter.hpp \
	$(srcdir)/MATLAB_MatrixWriter.hpp \
	$(srcdir)/MATLAB_VectorWriter.hpp \
	$(srcdir)/GWriter.hpp \
	$(srcdir)/GReader.hpp \
	$(srcdir)/ElementIntegrator.h \
	$(srcdir)/IntMap.h \
	$(srcdir)/MeshWriter.h \
	$(srcdir)/VectorWriter.h \
	$(srcdir)/ProblemReader.h \
	$(srcdir)/ProblemWriter.h \
	$(srcdir)/MatrixWriter.h  \
	$(srcdir)/SolutionWriter.h  \
	$(srcdir)/Tools.h

CORE = \
	$(srcdir)/fem.cpp \
	$(srcdir)/BoundaryCondition.cpp \
	$(srcdir)/ElementByElementMatrix.cpp \
	$(srcdir)/DistMesh.cpp \
	$(srcdir)/Jacobi.cpp \
	$(srcdir)/ElasticityProblem.cpp \
	$(srcdir)/EbeElasticityProblem.cpp \
	$(srcdir)/VbrElasticityProblem.cpp \
	$(srcdir)/GWriter.cpp \
	$(srcdir)/GReader.cpp \
	$(srcdir)/ElementIntegrator.cpp \
	$(srcdir)/IntMap.cpp  \
	$(srcdir)/Tools.cpp

PARFE_LIB = libparfe.a

PARFE_H = \
	$(CORE_H)

libparfe_a_SOURCES = \
	$(CORE) 

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libparfe_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(PARFE_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(PARFE_H) 

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 






