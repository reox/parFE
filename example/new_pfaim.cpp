/*
 * ParFE: a micro-FE solver for trabecular bone modeling
 * Copyright (C) 2006 ParFE developers, see
 * http://parfe.sourceforge.net/developers.php
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  
 * 02110-1301, USA.
 */

/*
 * This code should be used to perform a single finite element analysis
 * of a trabecular bone, stored in either ASCII or HDF5 format. Three
 * preconditioners are available: "Jacobi" (for comparisons only, not
 * for production), "ml" (matrix-ready multilevel preconditioner)
 * and "matrixfree" (matrix-free multilevel preconditioner).
 *
 * \author Marzio Sala
 *
 * \date Last modified on 15-Nov-06
 */

#include "parfe_ConfigDefs.h"
#include <string>
#include <fstream>
#include <vector>
#include <Epetra_Time.h>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <time.h>
#include "DistMesh.h"
#include "MEDIT_MeshWriter.hpp"
#include "PMVIS_MeshWriter.hpp"
#include "MEDIT_DispWriter.hpp"
#include "IBT_ProblemReader.hpp"
#include "IBT_ProblemWriter.hpp"
#include "MATLAB_MatrixWriter.hpp"
#include "MATLAB_VectorWriter.hpp"
#include "IBT_DispWriter.hpp"
#include "IBT_SolutionWriter.hpp"
#include "EbeElasticityProblem.h"
#include "VbrElasticityProblem.h"
#include "Jacobi.h"
#include "Tools.h"
#include "FEParameterList.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <ml_MultiLevelPreconditioner.h>
#include <ml_MatrixFreePreconditioner.h>
#include <AztecOO.h>

using namespace ML_Epetra;
using namespace Teuchos;

// =========== //
// main driver //
// =========== //

int main( int argc, char* argv[] )
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  // define an Epetra communicator
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // get the proc ID of this process
  int MyPID = comm.MyPID();

  // get the total number of processes
  int NumProc = comm.NumProc();

  Epetra_Time timer(comm);
  double start_time = timer.WallTime();
  double temp_time = start_time;

  //Macros to print information
#define INFO(mesg) \
  if (MyPID == 0) cout << "*** INFO ***: " << mesg << endl;

#define WARN(mesg) \
  if (MyPID == 0) cout << "*** WARNING ***: " << mesg << endl;

#define ERROR(mesg, exit_code)\
  if (MyPID == 0) cerr << "*** ERROR ***: " << mesg << endl; \
  exit(exit_code);		

  CommandLineProcessor CLP;

  char help[] = "This is the latest version of the ParFE driver.\n" 
    "It can solve problems stored in ASCII or HDF5 format, perform\n"
    "load balance, build the preconditioner for matrix-free or\n"
    "matrix-ready problems, and print out the solution.\n";
  CLP.setDocString(help);

  string inputFileName = "../mesh/simple_test.mesh";
  string outputFileName = "./output";
  string precType = "ml";
  int maxIters = 1550;
  double tolerance = 1.e-5;
  bool verbose = true;
  bool printMEDIT = false;
  bool printPMVIS = false;
  int outputFrequency = 16;

  CLP.setOption("filename", &inputFileName, "Name of input file");
  CLP.setOption("output", &outputFileName, "Base name of output file");
  CLP.setOption("precond", &precType, 
                "Preconditioner to be used [Jacobi/ml/matrixfree]");
  CLP.setOption("maxiters", &maxIters, "Maximum CG iterations");
  CLP.setOption("tolerance", &tolerance, "Tolerance for CG");
  CLP.setOption("verbose", "silent", &verbose, 
                "Set verbosity level");
  CLP.setOption("printmedit", "noprintmedit", &printMEDIT, 
                "Print solution in MEDIT format");
  CLP.setOption("printpmvis", "noprintpmvis", &printPMVIS, 
                "Print partition in PMVIS format");
  CLP.setOption("frequency", &outputFrequency, 
                "Prints out residual every specified iterations");

  CLP.recogniseAllOptions(false);
  CLP.throwExceptions(false);

  if (CLP.parse(argc, argv) == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
  {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  FEParameterList fe_param(NumProc);
  fe_param.set("input", outputFileName);

  Epetra_MultiVector* BCs = 0;

  ProblemReader* pr = 0;
  if (inputFileName.find(".mesh.h5") != string::npos)
    pr = new IBT_ProblemReader<HDF5_GReader>(inputFileName, comm);
  else if (inputFileName.find(".mesh") != string::npos)
    pr = new IBT_ProblemReader<C_ASCII_GReader>(inputFileName, comm);
  else
  {
    cerr << "File name (" << inputFileName << ") not in ASCII or HDF5 format" << endl;
    exit(EXIT_FAILURE);
  }

  fe_param.Scan(pr);

  INFO("#nodes       = " << fe_param.get<int>("#nodes"));
  INFO("#dimensions  = " << fe_param.get<int>("#dimensions"));
  INFO("#elements    = " << fe_param.get<int>("#elements"));

  temp_time = timer.WallTime();

  //If these parameters are still unset, provide some default values
  fe_param.get("iteration limit", maxIters);
  fe_param.get("tolerance", tolerance);

  DistMesh* mesh = new DistMesh(fe_param.get<int>("#nodes"), fe_param.get<int>("#dimensions"), fe_param.get<int>("#elements"), fe_param.get<int>("#nodes per element"), fe_param.get<IntMap>("material ids"), comm);

  mesh->Scan(pr);
  BoundaryCondition* bcond = new BoundaryCondition(fe_param.get<int>("#dofs per node"), comm);
  bcond->Scan(pr);

  delete pr;
  
  //--------------------------end of input--------------------------------
  
  INFO("Input Time: " << timer.WallTime() - temp_time);
  temp_time = timer.WallTime();

  mesh->Redistribute(true, verbose);
  if (fe_param.get<bool>("load balancing")) {
    INFO("Time used to partition mesh: " << timer.WallTime() - temp_time);
  } else {
    WARN("No load balancing used");
  }

  //redistribute the boundary condition data according to the node map
  bcond->Redistribute(mesh->NodeMap()->NumMyElements(), mesh->NodeMap()->MyGlobalElements());
  
  if (printMEDIT)
  {
    MEDIT_MeshWriter meshw(outputFileName + ".medit.mesh", 3, comm);
    mesh->Print(&meshw);
  }

  if (fe_param.get<bool>("print for PMVIS")) 
  {
    std::string con = outputFileName + ".con";
    std::string xyz = outputFileName + ".xyz";
    std::string part = outputFileName + ".part";
    PMVIS_MeshWriter meshw(xyz, con, part, comm);
    mesh->Print(&meshw);
  }

  int memused = meminfo();
  INFO("Mem used to hold mesh information: " << memused <<  "MB");

  INFO("Total Mesh setup time: " << timer.WallTime() - start_time);

  //-------------------------beginning of assembly------------------------
  
  temp_time = timer.WallTime();

  ElasticityProblem* ep;
  if (precType == "Jacobi" || precType == "matrixfree")
  {
    fe_param.set("element by element", true);
    ep = new EbeElasticityProblem(*mesh, fe_param);
  }
  else
    ep = new VbrElasticityProblem(*mesh, fe_param);

  ep->Impose(*bcond);

  BCs = new Epetra_MultiVector(*ep->GetRHS());

  INFO("Time used to assemble matrix: " << timer.WallTime() - temp_time);

  memused = meminfo() - memused;
  INFO("Mem used to hold stiffness matrix operator: " << memused << "MB");

  Epetra_Vector* X = dynamic_cast<Epetra_Vector*>(ep->GetLHS());
  Epetra_Vector* B = dynamic_cast<Epetra_Vector*>(ep->GetRHS());
  if (X == 0 || B == 0) 
  {
    ERROR("Use Epetra_Vector in ElastictyProblem", 1);
  }

  vector<double> null_space;

  //Compute null space for ml preconditioner
  if (precType == "ml" || precType == "matrixfree")
    set_null_space(mesh, null_space);

  RefCountPtr<Epetra_CrsGraph> graph;

  if (precType == "matrixfree") 
  {
    ElementByElementMatrix* ebe = dynamic_cast<ElementByElementMatrix*>(ep->GetOperator());
    if (ebe != 0)
    {
      graph = rcp(new Epetra_CrsGraph(Copy, ebe->Map(), 0));
      mesh->MatrixGraph(*(graph.get()));
    }
  }

  //------------------------beginning of solution------------------------ 

  RefCountPtr<Epetra_Operator> Prec;
  RefCountPtr<Epetra_Vector> diagonal;

  AztecOO solver(*ep);

  solver.SetAztecOption(AZ_solver, AZ_cg);
  solver.SetAztecOption(AZ_output, outputFrequency);

  temp_time = timer.WallTime();

  if (precType == "Jacobi") 
  {
    ElementByElementMatrix* ebe = dynamic_cast<ElementByElementMatrix*>(ep->GetOperator());

    // impose the boundary conditions for nodes with fixed displacement.
    // This is done by applying the operator without the fix for BC,
    // then resetting the nodes corresponding to fixed nodes in RHS.
    ebe->Apply_NoResetBoundaries(*B, *X);
    *B = *X;
    B->Scale(-1.0);

    ebe->ResetBoundaries(B);
    X->PutScalar(0.0);

    //Use user defined jacobi preconditioner
    diagonal = rcp(new Epetra_Vector(ebe->Map()));
    ebe->ExtractDiagonalCopy(*diagonal);
    Prec = rcp(new Jacobi(*diagonal));  
    solver.SetPrecOperator(Prec.get());

    for (int i = 0; i < BCs->MyLength(); ++i)
      if ((*BCs)[0][i] != 0.0)
        (*B)[i] = (*BCs)[0][i];
  }
  else if (precType == "matrixfree") 
  {
    ElementByElementMatrix* ebe = dynamic_cast<ElementByElementMatrix*>(ep->GetOperator());

    // impose the boundary conditions for nodes with fixed displacement.
    // This is done by applying the operator without the fix for BC,
    // then resetting the nodes corresponding to fixed nodes in RHS.
    ebe->Apply_NoResetBoundaries(*B, *X);
    *B = *X;
    B->Scale(-1.0);

    ebe->ResetBoundaries(B);
    X->PutScalar(0.0);

    ParameterList MLList;
    MLList.set("prec: type", "hybrid");
    MLList.set("smoother: type", "Chebyshev");
    MLList.set("low memory", true);
    MLList.set("aggregation: type", "METIS");
    MLList.set("aggregation: nodes per aggregate", 128);
    MLList.set("eigen-analysis: max iters", 10);
    MLList.set("smoother: degree", 5);
    MLList.set("AP allocation factor", 0.3);

    SetDefaults("SA", MLList.sublist("ML list"));
    MLList.sublist("ML list").set("aggregation: damping factor", 1.333);
    MLList.sublist("ML list").set("coarse: max size", 1024);
    MLList.sublist("ML list").set("smoother: type", "symmetric Gauss-Seidel");
    MLList.sublist("ML list").set("aggregation: type", "Uncoupled-MIS");
    MLList.sublist("ML list").set("coarse: type", "Amesos-KLU");
    MLList.sublist("ML list").set("low memory usage", true);
    MLList.sublist("ML list").set("cycle applications", 10);

    if (verbose)
      MLList.set("output", 10);
    else
      MLList.set("output", 0);
    MLList.sublist("ML list").set("max levels", 10);

    Epetra_MultiVector NullSpace(Copy, graph->Map(),
                                 &null_space[0], graph->Map().NumMyPoints(), 6);

    diagonal = rcp(new Epetra_Vector(ebe->Map()));
    ebe->ExtractDiagonalCopy(*diagonal);

    try
    {
      Prec = rcp(new MatrixFreePreconditioner(*ebe, *graph, NullSpace,
                                              *diagonal, MLList));
    }
    catch (...)
    {
      ERROR("Caught generic exception in the construction" << endl \
            << "of the matrix-free multilevel preconditioner," << endl \
            << "maybe this is a memory allocation error.", 1);
    }
    graph = null;
    null_space.resize(0);
    solver.SetPrecOperator(Prec.get());

    for (int i = 0; i < BCs->MyLength(); ++i)
      if ((*BCs)[0][i] != 0.0)
        (*B)[i] = (*BCs)[0][i];
  }
  else if (precType == "ml")
  {
    // allocate the ML preconditioner
    ParameterList MLList;
    SetDefaults("SA", MLList);
    MLList.set("smoother: type (level 0)", "Chebyshev");
    MLList.set("smoother: type (level 1)", "symmetric Gauss-Seidel");
    MLList.set("coarse: max size", 1024);
    MLList.set("null space: dimension", 6);
    MLList.set("null space: type", "pre-computed");
    MLList.set("null space: vectors", &null_space[0]);
    MLList.set("aggregation: type (level 0)", "METIS");
    MLList.set("aggregation: type (level 1)", "Uncoupled");
    MLList.set("aggregation: type (level 2)", "MIS");
    MLList.set("low memory usage", true);
    if (verbose) MLList.set("output", 10);
    else         MLList.set("output", 0);

    Prec = rcp(new MultiLevelPreconditioner(*ep->GetMatrix(), MLList, true));
    solver.SetPrecOperator(Prec.get());
    null_space.resize(0);
  } 
  else 
  {
    solver.SetAztecOption(AZ_precond, AZ_none);
    WARN("No preconditioner used");
  }

  INFO("Time used to build preconditioner: " << timer.WallTime() - temp_time);
  memused = meminfo() - memused;
  INFO("Mem used to hold preconditioner: " << memused << "MB");

  X->PutScalar(0.0);

  solver.Iterate(maxIters, tolerance);

  if (precType == "matrixfree" || precType == "Jacobi")
  {
    ElementByElementMatrix* ebe = dynamic_cast<ElementByElementMatrix*>(ep->GetOperator());
    ebe->ResetRestrainedBoundaries(X);
  }

  diagonal = null;
  Prec = null;

  //------------------------------Print solution-------------------------

  temp_time = timer.WallTime();

  if (inputFileName.find(".mesh.h5") != string::npos)
  {
    IBT_SolutionWriter<HDF5_GWriter> sw(inputFileName, comm);
    ep->PrintSolution(&sw, fe_param.get<Epetra_SerialDenseMatrix>("material properties"));
  }

  if (printMEDIT)
  {
    MEDIT_DispWriter dispw(outputFileName + ".medit.bb", comm, 3);
    ep->PrintLHS(&dispw);
  }

  INFO("output time: " << timer.WallTime() -temp_time);

  comm.Barrier();
  INFO("Total time used: " << timer.WallTime() -start_time);

  // ======================= //
  // Finalize MPI and exit //
  // ----------------------- //

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  exit(EXIT_SUCCESS);
}
