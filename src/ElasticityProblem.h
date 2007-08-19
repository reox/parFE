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

#ifndef _ELASTICITY_PROBLEM_H_
#define _ELASTICITY_PROBLEM_H_
#include <Epetra_LinearProblem.h>
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

class ElasticityProblem : public Epetra_LinearProblem
{
 public:
  //! ElasticityProblem constructor

  /*! The constructor assembles the global stiffness, imposes the boundary conditions and provides a
      complete linear problem Ax = b to be solved.
      After the constructor was called, the input datastructures can be deleted to save memory.

      \param mesh
      (In) Object containing mesh datastructures.
      
      \param bc
      (In) Object containing boundary condition datatrsuctures.

      \param param
      (In) A  Teuchos::ParameterList containing the relevant parameters

      \param comm
      (In) A Epetra_Comm object to be associated with this ElasticityProblem

      \return Pointer the created ElasticityProblem.
  */
  ElasticityProblem(DistMesh& mesh, Teuchos::ParameterList& param);

  //! ElasticityProblem destructor
  virtual ~ElasticityProblem();

  //! Print left hand side to a file
  /*!
    \param vw
    (In) Pointer to the VectorWriter used to print the lhs vector to a file.
  */
  void PrintLHS(VectorWriter* vw);
  
  //! Print right hand side to a file
  /*!
    \param vw
    (In) Pointer to the VectorWriter used to print the rhs vector to a file.
  */
  void PrintRHS(VectorWriter* vw);

  //! Print the matrix into a file
   /*!
    \param mw
    (In) Pointer to the MatrixWriter used to print the matrix to a file.
  */
  void PrintMatrix(MatrixWriter* mw);

  //! Impose boundary condition
  virtual int Impose(BoundaryCondition& bcond) = 0;

  //! Remove all boundary conditions from the system
  virtual int Restore() = 0;

  //! Reference to the object containing mesh datastructures.
  DistMesh& mesh;

  //! Number of DOFs per node
  const int NumDofsPerNode;
  //! Number of spatial dimensions
  const int Dimension;
  //! Global number of nodes
  const int NumNodes;
  //! Number of nodes per element
  const int NumNodesPerElement;
  //! Number of DOFs per element
  const int NumDofsPerElement;

  const Epetra_SerialDenseMatrix MaterialProperties;

  //! Weight vector to cancel out penalty terms from the residual during an iterative solution
  Epetra_Vector* ResidualWeightVector;

  void PrintSolution(SolutionWriter* sw, Epetra_SerialDenseMatrix& MaterialProps);

 
 protected:
 virtual int Assemble() = 0;
  int SetBoundaryConditions(BoundaryCondition& bcond, Epetra_Vector* diagonal, Epetra_Vector* rhs, double penalty=1.0);
  Epetra_MpiComm comm;
  ElementIntegrator ElementMatrix;
    
};


#endif
