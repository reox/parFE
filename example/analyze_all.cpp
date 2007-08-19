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
 * This code should be used to perform the finite element analysis
 * of a set of trabecular bones, stored in HDF5 format. The preconditioner
 * is based on matrix-free multigrid.
 *
 * \author Marzio Sala
 *
 * \date Last modified on 15-Nov-06
 */

#include "parfe_ConfigDefs.h"
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
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
#include "IBT_ProblemReader.hpp"
#include "IBT_ProblemWriter.hpp"
#include "IBT_DispWriter.hpp"
#include "IBT_SolutionWriter.hpp"
#include "EbeElasticityProblem.h"
#include "Tools.h"
#include "FEParameterList.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <ml_MatrixFreePreconditioner.h>
#include <AztecOO.h>

using namespace ML_Epetra;
using namespace Teuchos;

void Tokenize(const string& str,
              vector<string>& tokens,
              const string& delimiters = " ")
{
  // Skip delimiters at beginning.
  string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  // Find first "non-delimiter".
  string::size_type pos     = str.find_first_of(delimiters, lastPos);

  while (string::npos != pos || string::npos != lastPos)
  {
    // Found a token, add it to the vector.
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
  }
}

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

  char help[] = "This example should be used to analyze\n"
    "a set of problems stored in HDF5 format. The list of\n"
    "input files is given through the --filename option, where\n"
    "each input file is separated by commas. For example, one can\n"
    "write --filename=cube1.mesh.h5,cube2.mesh.h5\n";
  CLP.setDocString(help);

  string inputFileNames = "../mesh/simple_test.mesh";
  int maxIters = 1550;
  double tolerance = 1.e-5;
  bool verbose = false;
  int outputFrequency = 16;

  CLP.setOption("filename", &inputFileNames, "Name of input file(s)");
  CLP.setOption("maxiters", &maxIters, "Maximum CG iterations");
  CLP.setOption("tolerance", &tolerance, "Tolerance for CG");
  CLP.setOption("verbose", "silent", &verbose, 
                "Set verbosity level");
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

  vector<string> tokens;
  vector<string>::iterator it;

  Tokenize(inputFileNames, tokens, ",");

  // ----------------------------------------------------- //
  // loop over all specified input files. If an input file //
  // is not in HDF5 format, skip it.                       //
  // ----------------------------------------------------- //
  
  for (it = tokens.begin(); it != tokens.end(); ++it)
  {
    string inputFileName = *it;
    INFO("Processing file " << inputFileName);
    if (inputFileName.find(".mesh.h5") == string::npos)
    {
      INFO("File NOT in HDF5 format; skip it");
      continue;
    }

    fe_param.set("input", inputFileName);

    ProblemReader* pr = new IBT_ProblemReader<HDF5_GReader>(inputFileName, comm);

    fe_param.Scan(pr);

    temp_time = timer.WallTime();

    fe_param.get("iteration limit", maxIters);
    fe_param.get("tolerance", tolerance);

    DistMesh mesh(fe_param.get<int>("#nodes"), fe_param.get<int>("#dimensions"), fe_param.get<int>("#elements"), fe_param.get<int>("#nodes per element"), fe_param.get<IntMap>("material ids"), comm);

    mesh.Scan(pr);
    BoundaryCondition* bcond = new BoundaryCondition(fe_param.get<int>("#dofs per node"), comm);
    bcond->Scan(pr);

    delete pr;

    //--------------------------end of input--------------------------------

    INFO("Input Time: " << timer.WallTime() - temp_time);
    temp_time = timer.WallTime();

    mesh.Redistribute(true, verbose);
    bcond->Redistribute(mesh.NodeMap()->NumMyElements(), mesh.NodeMap()->MyGlobalElements());

    INFO("Time used to partition mesh: " << timer.WallTime() - temp_time);

    int memused = meminfo();
    INFO("Mem used to hold mesh information: " << memused <<  "MB");

    INFO("Total Mesh setup time: " << timer.WallTime() - start_time);

    //-------------------------beginning of assembly------------------------

    temp_time = timer.WallTime();

    EbeElasticityProblem ep(mesh, fe_param);
    fe_param.set("element by element", true);

    ep.Impose(*bcond);
    delete bcond;

    Epetra_MultiVector BCs(*ep.GetRHS());

    INFO("Time used to assemble matrix: " << timer.WallTime() - temp_time);

    memused = meminfo() - memused;
    INFO("Mem used to hold stiffness matrix operator: " << memused << "MB");

    ElementByElementMatrix* ebe = dynamic_cast<ElementByElementMatrix*>(ep.GetOperator());
    Epetra_Vector* X = dynamic_cast<Epetra_Vector*>(ep.GetLHS());
    Epetra_Vector* B = dynamic_cast<Epetra_Vector*>(ep.GetRHS());
    assert (X != 0); assert (B != 0);

    // ----------------------------------------------------------------- //
    // Build the preconditioner for matrix-free multigrid.               //
    // Impose the boundary conditions for nodes with fixed displacement. //
    // This is done by applying the operator without the fix for BC,     //
    // then resetting the nodes corresponding to fixed nodes in RHS.     //
    // ----------------------------------------------------------------- //

    temp_time = timer.WallTime();

    RefCountPtr<Epetra_Operator> Prec;
    RefCountPtr<Epetra_Vector> diagonal;
    RefCountPtr<Epetra_CrsGraph> graph;

    vector<double> null_space;
    set_null_space(&mesh, null_space);

    graph = rcp(new Epetra_CrsGraph(Copy, ebe->Map(), 0));
    mesh.MatrixGraph(*(graph.get()));

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
    MLList.set("eigen-analysis: max iters", 20);
    MLList.set("smoother: degree", 3);
    MLList.set("AP allocation factor", 0.5);
    MLList.set("aggregation: nodes per aggregate", 128);

    SetDefaults("SA", MLList.sublist("ML list"));
    MLList.sublist("ML list").set("aggregation: damping factor", 1.333);
    MLList.sublist("ML list").set("coarse: max size", 1024);
    MLList.sublist("ML list").set("smoother: type", "symmetric Gauss-Seidel");
    MLList.sublist("ML list").set("aggregation: type", "Uncoupled-MIS");
    MLList.sublist("ML list").set("coarse: type", "Amesos-KLU");
    MLList.sublist("ML list").set("low memory usage", true);
    MLList.sublist("ML list").set("cycle applications", 10);
    MLList.sublist("ML list").set("max levels", 10);

    if (verbose)
    {
      MLList.set("output", 10);
      MLList.sublist("ML list").set("output", 10);
    }
    else
    {
      MLList.set("output", 0);
      MLList.sublist("ML list").set("output", 0);
    }

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

    for (int i = 0; i < BCs.MyLength(); ++i)
      if (BCs[0][i] != 0.0)
        (*B)[i] = BCs[0][i];

    INFO("Time used to build preconditioner: " << timer.WallTime() - temp_time);
    memused = meminfo() - memused;
    INFO("Mem used to hold preconditioner: " << memused << "MB");

    // ----------------------------------------------------------------- //
    // setup AztecOO's conjugate gradient, using zero starting solution. //
    // Output level is as specified via command line option.             //
    // After convergence, destroy objects and reset the boundaries       //
    // related to restrained nodes.                                      //
    // ----------------------------------------------------------------- //
    
    AztecOO solver(ep);

    solver.SetAztecOption(AZ_solver, AZ_cg);
    solver.SetAztecOption(AZ_output, outputFrequency);
    solver.SetPrecOperator(Prec.get());

    X->PutScalar(0.0);

    solver.Iterate(maxIters, tolerance);

    ebe->ResetRestrainedBoundaries(X);

    diagonal = null;
    Prec = null;

    // --------------------------- //
    // print solution in HDF5 file //
    // --------------------------- //

    temp_time = timer.WallTime();

    IBT_SolutionWriter<HDF5_GWriter> sw(inputFileName, comm);
    ep.PrintSolution(&sw, fe_param.get<Epetra_SerialDenseMatrix>("material properties"));

    INFO("output time: " << timer.WallTime() -temp_time);

    comm.Barrier();
    INFO("Total time used for " << inputFileName << ": " << timer.WallTime() -start_time);
  }

  // ----------------------- //
  // Finalize MPI and exit //
  // ----------------------- //

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  exit(EXIT_SUCCESS);
}
