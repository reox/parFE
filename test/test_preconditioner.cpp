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
 * This test aims to compare the "classical" ML preconditioner with 
 * the matrix-free version. A command line parameter should specify the
 * input file (typically contained in parfe/mesh) and the preconditioner
 * type. Note that more convenient settings for both preconditioners could be
 * found; however what we want to do here is only to compare them and be
 * sure that we are solving the same problem.
 *
 * \author Marzio Sala
 *
 * \date November 2006
 */

//include variables set by configure
#include "parfe_ConfigDefs.h"
#include <string>
#include <fstream>
#include <vector>
#ifdef HAVE_MPI
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_Vector.h>
#include "DistMesh.h"
#include <AztecOO.h>
#include "Tools.h"
#include "IBT_ProblemReader.hpp"
#include "EbeElasticityProblem.h"
#include "VbrElasticityProblem.h"
#include "MEDIT_DispWriter.hpp"
#include "MEDIT_MeshWriter.hpp"
#include <unistd.h>
#include "FEParameterList.hpp"
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <ml_MultiLevelPreconditioner.h>
#include <ml_MatrixFreePreconditioner.h>

using namespace ML_Epetra;
using namespace Teuchos;

inline void EXIT(const int exit_code)
{
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(exit_code);
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

  Teuchos::CommandLineProcessor CLP;

  CLP.setDocString("This test compares the multilevel preconditioners");
  string fileName = "../mesh/simple_test.mesh";
  string precType = "ml";
  int maxIters = 1550;
  double tolerance = 1.e-8;

  CLP.setOption("filename", &fileName, "name of input file");

  CLP.recogniseAllOptions(true);
  CLP.throwExceptions(false);

  CLP.parse(argc, argv);

  FEParameterList fe_param(comm.NumProc());
  //set default values
  fe_param.set("#processors", comm.NumProc());
  fe_param.set("load balancing", true);
  fe_param.set("hdf5 I/O", true);
  fe_param.set("solve", true);
  fe_param.set("print for MEDIT", false);
  fe_param.set("print for PMVIS", false);
  fe_param.set("print matrix", false);
  fe_param.set("print lhs", false);
  fe_param.set("print rhs", false);
  fe_param.set("convert to ascii", false);
  fe_param.set("convert to hdf5", false);
  fe_param.set("verbose", true);
  fe_param.set("element by element", false);
  fe_param.set("input", fileName);

#if 0
  if (precType == "ml")
  {
    fe_param.set("preconditioner", "ml");
    fe_param.set("element by element", false);
  }
  else
  {
    fe_param.set("preconditioner", "matrixfree");
    fe_param.set("element by element", true);
  }
#endif

  Epetra_MultiVector* BCs = 0;

  ProblemReader* pr = 0;
  if (fileName.find(".mesh.h5") != string::npos)
    pr = new IBT_ProblemReader<HDF5_GReader>(fileName, comm);
  else if (fileName.find(".mesh") != string::npos)
    pr = new IBT_ProblemReader<C_ASCII_GReader>(fileName, comm);
  else
  {
    cerr << "File name (" << fileName << ") not in ASCII or HDF5 format" << endl;
    EXIT(EXIT_FAILURE);
  }

  fe_param.Scan(pr);

  //If these parameters are still unset, provide some default values
  fe_param.get("iteration limit", 1000);
  fe_param.get("tolerance", 1e-6);

  DistMesh* mesh = new DistMesh(fe_param.get<int>("#nodes"), fe_param.get<int>("#dimensions"), fe_param.get<int>("#elements"), fe_param.get<int>("#nodes per element"), fe_param.get<IntMap>("material ids"), comm);

  mesh->Scan(pr);
  BoundaryCondition* bcond = new BoundaryCondition(fe_param.get<int>("#dofs per node"), comm);
  bcond->Scan(pr);

  delete pr;

  mesh->Redistribute(fe_param.get<bool>("load balancing"), fe_param.get<bool>("verbose"));

  bcond->Redistribute(mesh->NodeMap()->NumMyElements(), mesh->NodeMap()->MyGlobalElements());

  //------------------------- null space ------------------------
  
  int nullSpaceDimension = 6;
  vector<double> null_space;

  set_null_space(mesh, null_space);

  // -------------------- //
  // Matrix-free solution //
  // -------------------- //

  EbeElasticityProblem ep_free(*mesh, fe_param);

  ep_free.Impose(*bcond);

  BCs = new Epetra_MultiVector(*ep_free.GetRHS());

  RefCountPtr<Epetra_CrsGraph> graph;

  ElementByElementMatrix* ebe = dynamic_cast<ElementByElementMatrix*>(ep_free.GetOperator());
  graph = rcp(new Epetra_CrsGraph(Copy, ebe->Map(), 0));
  mesh->MatrixGraph(*(graph.get()));

  Epetra_Vector* X_free = dynamic_cast<Epetra_Vector*>(ep_free.GetLHS());
  Epetra_Vector* B_free = dynamic_cast<Epetra_Vector*>(ep_free.GetRHS());

  RefCountPtr<Epetra_Operator> prec_free;
  RefCountPtr<Epetra_Vector> diagonal;

  AztecOO solver_free(ep_free);

  solver_free.SetAztecOption(AZ_solver, AZ_cg);
  solver_free.SetAztecOption(AZ_output, 16);

  ebe->Apply_NoResetBoundaries(*B_free, *X_free);
  *B_free = *X_free;
  B_free->Scale(-1.0);

  ebe->ResetBoundaries(B_free);
  X_free->PutScalar(0.0);

  ParameterList MLList_free;
  MLList_free.set("prec: type", "hybrid");
  MLList_free.set("smoother: type", "Chebyshev");
  MLList_free.set("smoother: degree", 3);
  MLList_free.set("low memory", true);
  MLList_free.set("aggregation: type", "Uncoupled");
  MLList_free.set("eigen-analysis: max iters", 10);

  MLList_free.sublist("ML list").set("max levels", 10);

  SetDefaults("SA", MLList_free.sublist("ML list"));
  MLList_free.sublist("ML list").set("aggregation: damping factor", 0.0);
  MLList_free.sublist("ML list").set("coarse: max size", 64);
  MLList_free.sublist("ML list").set("smoother: type", "symmetric Gauss-Seidel");
  MLList_free.sublist("ML list").set("aggregation: type", "Uncoupled");
  MLList_free.sublist("ML list").set("coarse: type", "Amesos-KLU");
  MLList_free.sublist("ML list").set("low memory usage", true);

  MLList_free.set("output", 0);
  MLList_free.sublist("ML list").set("output", 0);

  MLList_free.set("AP allocation factor", 0.5);
  MLList_free.set("aggregation: nodes per aggregate", 27);
  MLList_free.sublist("ML list").set("cycle applications", 1);

  Epetra_MultiVector NullSpace(Copy, graph->Map(), &null_space[0], 
                               graph->Map().NumMyPoints(), nullSpaceDimension);

  diagonal = rcp(new Epetra_Vector(ebe->Map()));
  ebe->ExtractDiagonalCopy(*diagonal);

  try
  {
    prec_free = rcp(new MatrixFreePreconditioner(*ebe, *graph, NullSpace,
                                                 *diagonal, MLList_free));
  }
  catch (...)
  {
    cerr << "Caught generic exception in the construction" << endl 
      << "of the matrix-free multilevel preconditioner," << endl
      << "maybe this is a memory allocation error." << endl;
    exit(EXIT_FAILURE);
  }
  solver_free.SetPrecOperator(prec_free.get());

  for (int i = 0; i < BCs->MyLength(); ++i)
    if ((*BCs)[0][i] != 0.0)
      (*B_free)[i] = (*BCs)[0][i];

  X_free->PutScalar(0.0);

  solver_free.Iterate(maxIters, tolerance);

  ebe->ResetRestrainedBoundaries(X_free);

  prec_free = null;

  // --------------------- //
  // matrix-ready solution //
  // --------------------- //

  VbrElasticityProblem ep_ready(*mesh, fe_param);

  ep_ready.Impose(*bcond);

  Epetra_Vector* X_ready = dynamic_cast<Epetra_Vector*>(ep_ready.GetLHS());
  Epetra_Vector* B_ready = dynamic_cast<Epetra_Vector*>(ep_ready.GetRHS());

  RefCountPtr<Epetra_Operator> prec_ready;

  ParameterList MLList_ready;
  SetDefaults("SA", MLList_ready);
  MLList_ready.set("smoother: type (level 0)", "Chebyshev");
  MLList_ready.set("smoother: sweeps (level 0)", 3);
  MLList_ready.set("smoother: type (level 1)", "symmetric Gauss-Seidel");
  MLList_ready.set("aggregation: damping factor", 0.0);
  MLList_ready.set("coarse: max size", 64);
  MLList_ready.set("null space: dimension", nullSpaceDimension);
  MLList_ready.set("null space: type", "pre-computed");
  MLList_ready.set("null space: vectors", &null_space[0]);
  MLList_ready.set("aggregation: type", "Uncoupled");
  MLList_ready.set("aggregation: nodes per aggregate", 27);
  MLList_ready.set("low memory usage", true);
  MLList_ready.set("output", 0);
  prec_ready = rcp(new MultiLevelPreconditioner(*ep_ready.GetMatrix(), MLList_ready, true));

  X_ready->PutScalar(0.0);

  AztecOO solver_ready(ep_ready);

  solver_ready.SetPrecOperator(prec_ready.get());

  solver_ready.SetAztecOption(AZ_solver, AZ_cg);
  solver_ready.SetAztecOption(AZ_output, 16);
  solver_ready.Iterate(maxIters, tolerance);

  prec_ready = null;

  // ----------------- //
  // compare solutions //
  // ----------------- //
  
  double norm2_diff, norm2_ready, norm2_scaled;
  X_ready->Norm2(&norm2_ready);
  X_ready->Update(1.0, *X_free, -1.0);
  X_ready->Norm2(&norm2_diff);
  norm2_scaled = norm2_diff / norm2_ready;

  if (comm.MyPID() == 0)
  {
    cout << "||X_ready||_2                           = " << norm2_ready << endl;
    cout << "||X_ready - X_ free||_2                 = " << norm2_diff << endl;
    cout << "||X_ready - X_ free||_2 / ||X_ready||_2 = " << norm2_scaled << endl;
  }

  if (norm2_scaled > 1e-4)
  {
    if (comm.MyPID() == 0)
      cout << "*** TEST FAILED! ***" << endl;
    exit(EXIT_FAILURE);
  }
  else
  {
    if (comm.MyPID() == 0)
      cout <<  "*** TEST PASSED! ***" << endl;
  }

#if 0
  // uncomment to print out using MEDIT
  MEDIT_MeshWriter meshw("output.mesh", 3, comm);
  mesh->Print(&meshw);

  MEDIT_DispWriter dispw("output.bb", comm, 3);
  ep_ready.PrintLHS(&dispw);
#endif

  // free what is left
  delete mesh;
  graph = null;
  diagonal = null;
  null_space.resize(0);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  exit(EXIT_SUCCESS);
}
