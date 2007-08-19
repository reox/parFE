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

#include <iostream>
#include <vector>

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_MpiComm.h"


//! A Jacobi preconditioner class

/*! The Jacobi class reperesents the Jacobi preconditioner.
    The diagonal of the matrix is supplied with the constructor.
    When applied to a vector, the Jacobi preconditioner performs a simple element-wise multiplication
    with the reciprocal elements of the diagonal
*/  

class Jacobi : public Epetra_Operator
{
  public:

  //! Jacobi constructor

  /*!
    \param diagonal
    (In) The diagonal of a matrix.
    
    \return Pointer to a Jacobi object.
  */
  Jacobi(Epetra_Vector& diagonal);
  
  //! Jacobi destructor
  ~Jacobi();
  
  //! Does nothing. Only implemented for convenience.
  int SetUseTranspose(bool UseTranspose);
  
  //! Returns the result of the Jacobi preconditioner applied to a Epetra_MultiVector X in Y.
  /*!
    \param X
    (In) A Epetra_MultiVector of dimension NumVectors to apply to the Jacobi preconditioner.
    
    \param Y
    (Out) A Epetra_MultiVector of dimension NumVectors containing the result.

    \return Integer error code, set to 0 if successful.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  //! Returns -1. Only implemented for convenience.
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  //! Return -1. Only implemented for convenience.
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
  
 private:
  Epetra_Vector& diagonal;
  Epetra_Map* range_domain_map;
  int BlockMap2PointMap(const Epetra_BlockMap* BlockMap, Epetra_Map * & PointMap) const;

};
