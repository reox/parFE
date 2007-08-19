/*
 * ParFE: a micro-FE solver for trabecular bone modeling
 * Copyright (C) 2006, Uche Mennel and Marzio Sala
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

#include "ElasticityProblem.h"
#include <Epetra_Vector.h>
#include "DistMesh.h"
#include "BoundaryCondition.h"
#include "ElementByElementMatrix.h"
#include <Epetra_MpiComm.h>
#include <Epetra_FEVbrMatrix.h>
#include <Epetra_SerialDenseMatrix.h>
#include "MatrixWriter.h"
#include "VectorWriter.h"
#include "SolutionWriter.h"
#include "fem.h"
#include <Teuchos_ParameterList.hpp>

//! A class to combine parameters, mesh and boundary data into a Epetra_LinearProblem

class EbeElasticityProblem : public ElasticityProblem
{
 public:
  //! ElasticityProblem constructor

  /*! The constructor assembles the global stiffness, imposes the boundary conditions and provides a
      complete linear problem Ax = b to be solved.
      After the constructor was called, the input datastructures can be deleted to save memory.

      \param mesh
      (In) Object containing mesh datastructures.

      \param param
      (In) A  Teuchos::ParameterList containing the relevant parameters

      \return Pointer the created ElasticityProblem.
  */
  EbeElasticityProblem(DistMesh& mesh, Teuchos::ParameterList& param);

  ~EbeElasticityProblem();

  //! Impose boundary conditions
  int Impose(BoundaryCondition& bcond);
  
  //! Remove all boundary conditions from the system
  int Restore();

 private:
  int Assemble();

};


