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
#include <Epetra_IntVector.h>
#include "DistMesh.h"
#include "BoundaryCondition.h"
#include "MEDIT_MeshWriter.hpp"
#include "PMVIS_MeshWriter.hpp"
#include "MEDIT_DispWriter.hpp"
#include "IBT_ProblemReader.hpp"
#include "IBT_ProblemWriter.hpp"
#include "IBT_DispWriter.hpp"
#include "IBT_SolutionWriter.hpp"
#include "Tools.h"
#include "FEParameterList.hpp"
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RefCountPtr.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

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

  Teuchos::CommandLineProcessor CLP;

  char help[] = "This executable should be used to convert mesh\n"
    "to/from ASCII and HDF5 format. It reads the mesh in the given\n"
    "format (say, HDF5), and creates the equivalent mesh in the\n"
    "other format (in this case, ASCII)";
  CLP.setDocString(help);

  string fileName = "../mesh/simple_test.mesh";
  string output = "../mesh/simple_test";

  CLP.setOption("filename", &fileName, "Name of input file");
  CLP.setOption("output", &output, "Base name of output file");

  CLP.recogniseAllOptions(false);
  CLP.throwExceptions(false);

  if (CLP.parse(argc, argv) == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED)
  {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    exit(EXIT_SUCCESS);
  }

  string what = "not-set";

  FEParameterList fe_param(NumProc);

  ProblemReader* pr = 0;
  if (fileName.find(".mesh.h5") != string::npos)
  {
    what = "toASCII";;
    pr = new IBT_ProblemReader<HDF5_GReader>(fileName, comm);
  }
  else if (fileName.find(".mesh") != string::npos)
  {
    what = "toHDF5";;
    pr = new IBT_ProblemReader<C_ASCII_GReader>(fileName, comm);
  }
  else
  {
    cerr << "File name (" << fileName << ") not in ASCII or HDF5 format" << endl;
    exit(EXIT_FAILURE);
  }

  fe_param.Scan(pr);

  if (MyPID == 0)
  {
    cout << "#nodes       = " << fe_param.get<int>("#nodes") << endl;
    cout << "#dimensions  = " << fe_param.get<int>("#dimensions") << endl;
    cout << "#elements    = " << fe_param.get<int>("#elements") << endl;
  }

  DistMesh mesh(fe_param.get<int>("#nodes"), fe_param.get<int>("#dimensions"), fe_param.get<int>("#elements"), fe_param.get<int>("#nodes per element"), fe_param.get<IntMap>("material ids"), comm);

  mesh.Scan(pr);
  BoundaryCondition bcond(fe_param.get<int>("#dofs per node"), comm);
  bcond.Scan(pr);

  delete pr;
  
  //--------------------------end of input--------------------------------
  
  mesh.Redistribute(true, true);

  //redistribute the boundary condition data according to the node map
  bcond.Redistribute(mesh.NodeMap()->NumMyElements(), 
                     mesh.NodeMap()->MyGlobalElements());
  
  ProblemWriter* pw = 0;
  if (what == "toHDF5")
    pw = new IBT_ProblemWriter<HDF5_GWriter>(output + ".mesh.h5", comm);
  else 
    pw = new IBT_ProblemWriter<C_ASCII_GWriter>(output + ".mesh", comm);

  fe_param.Print(pw);
  mesh.Print(pw);
  bcond.Print(pw);

  delete pw;

#ifdef HAVE_MPI
  MPI_Finalize();
#endif
  exit(EXIT_SUCCESS);
}
