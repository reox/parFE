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

#ifndef _ELEMENT_BY_ELEMENT_MATRIX_H_
#define _ELEMENT_BY_ELEMENT_MATRIX_H_
#include <iostream>
#include <vector>

// \todo: find better solution
#define EBE_PENALTY 1.0e20

#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_Export.h>
#include <Epetra_Import.h>
#include <Epetra_Operator.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_MpiComm.h>
#include "ElementIntegrator.h"

//! A class to perform an EBE matrix-vector product

/*! The ElementByElementMatrix is used to multiply a vector by the global stiffness matrix, whithout ever
    assembling the matrix.
    Instead, chunks of the vector are gathered and multiplied by the corresponding element stiffness matrix.
    Afterwards, the resulting products are scattered and summed into the result vector.
 */

class ElementByElementMatrix : public Epetra_Operator, public Epetra_SrcDistObject
{
 public:

  //! ElementByElementMatrix constructor
  
   /*!
    \param a_dof_map
    (In) A Epetra_BlockMap indicating how the DOFs are distributed

    \param element_matrix
    (In) A ElementIntergrator object that is able to compute any needed element stiffness matrices.
    Normally, the element stiffness matrix has to be computed only once.
    
    \param elem2node
    (In) This Epetra_IntVector contains the element-to-node table.

     \param mat_ids
    (In) This Epetra_IntVector contains the material IDs.

    \return Pointer to the ElementByElementMatrix object.

  */
  ElementByElementMatrix(Epetra_BlockMap& a_dof_map, ElementIntegrator& element_matrix, Epetra_IntVector* elem2node, Epetra_IntVector* mat_ids);
  
  //! ElementByElementMatrix destructor
  ~ElementByElementMatrix();
  
  //! Return -1. Implemented only for convenience.
  int SetUseTranspose(bool UseTranspose);
  
  //! Returns the result of the ElementByElementMatrix multiplied with a Epetra_MultiVector X in Y.
  /*!
    \param X
    (In) A Epetra_MultiVector of dimension NumVectors to be multiplied by the ElementByElementMatrix.
    
    \param Y
    (Out) A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  //! Applies the operator but does not impose any Dirichlet boundary condition.
  int Apply_NoResetBoundaries(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  //! Impose the boundary conditions to vector x and store the result in y.
  
  /*!
    \param y
    (Out) Pointing to resulting vector.

    \param x
    (In) Pointer to the original vector.
  */
  void SetBoundaries(Epetra_Vector* y, const Epetra_Vector* x) const;
  
  //! Reset the values on the boundaries.
  void ResetBoundaries(Epetra_Vector* y) const;

  //! Reset the values on the boundaries, for restrained nodes only.
  void ResetRestrainedBoundaries(Epetra_Vector* y) const;

  //! Reset the values on the boundaries, for fixed nodes only.
  void ResetFixedBoundaries(Epetra_Vector* y) const;

  //! Returns -1. Only implemented for convenience.
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  //! Returns -1.0. Only implemented for convenience.
  double NormInf() const;
  
  //! Returns a character string describing the operator. 
  const char* Label() const;
  
  //! Returns false
  bool UseTranspose() const;
  
  //! Returns false
  bool HasNormInf() const;
  
  //! Returns the Epetra_Comm object associated with this operator.
  const Epetra_Comm& Comm() const;
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map& OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator. 
  const Epetra_Map& OperatorRangeMap() const;
  
  //! Returns the Epetra_BlockMap object associated with this operator. 
  const Epetra_BlockMap& Map() const;

  //! Returns a copy of the main diagonal in a user-provided vector.
  /*! 
  \param Diagonal
         (Out) Extracted main diagonal.

  \return Integer error code, set to 0 if successful.
  */
   int ExtractDiagonalCopy(Epetra_Vector & Diagonal) const ;

  //! Replaces diagonal values of the with those in the user-provided vector.
  /*! This routine is meant to allow replacement of existing diagonal values.
      The Epetra_Map associated with the input Epetra_Vector must be compatible
with the RowMap of the matrix.

      \param Diagonal
      (In) New values to be placed in the main diagonal.

      \return Integer error code, set to 0 if successful.
  */
  int ReplaceDiagonalValues(Epetra_Vector& Diagonal) const;


  int PrepareApply(Epetra_Vector* fixed_dofs, Epetra_IntVector* restrained_dofs);

  //! Penalty value used to impose the boundary conditions
  const double Penalty;

 private:
  int BlockMap2PointMap(const Epetra_BlockMap& BlockMap, Epetra_Map * & PointMap) const;
  void gather( const Epetra_Vector* p, int matid ) const;
  void gather2( const Epetra_Vector* p, int matid ) const;
  void scatter( Epetra_Vector* u, int matid ) const;
  
  Epetra_BlockMap dof_map;
  Epetra_Map* domain_map;
  Epetra_Map* range_map;
  Epetra_Vector* diagonal;
  Epetra_Vector* pmul_vec;
  const Epetra_MpiComm comm;
  
  int ndof;
  int nelem;

  Epetra_Import* importer;
  double* g_pmul;
  double* u_tmp;
  ElementIntegrator ElementMatrix;

  //datastructures for precalculation of indices
  std::vector<int> fixed_domain_inds;
  std::vector<int> restrained_domain_inds;
  std::vector<int> scatter_inds;
  std::vector<int> gather_inds;
  int* inds_ptr;
  int* N;
};


#endif
