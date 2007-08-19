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

#include "Jacobi.h"


Jacobi::Jacobi(Epetra_Vector& a_diagonal)
  : diagonal(a_diagonal)
{
  BlockMap2PointMap(&diagonal.Map(), range_domain_map);
}


Jacobi::~Jacobi()
{
  delete range_domain_map;
}


int Jacobi::SetUseTranspose(bool UseTranspose)
{
  //EPETRA_CHK_ERR(-1);
  return -1;
}


int Jacobi::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  assert (X.NumVectors() == Y.NumVectors());
  Y.Multiply(1.0, X, diagonal, 0.0);

  return(0);
}

 int Jacobi::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  assert (X.NumVectors() == Y.NumVectors());
  Y.ReciprocalMultiply(1.0, diagonal, X, 0.0);

  return(0);
}


double Jacobi::NormInf() const
{
  return(-1.0);
}

const char* Jacobi::Label() const
{
  return("Jacobi Preconditioner");
}

bool Jacobi::UseTranspose() const
{
  return(false);
}

bool Jacobi::HasNormInf() const
{
  return(false);
}

const Epetra_Comm& Jacobi::Comm() const
{
  return dynamic_cast<const Epetra_MpiComm&>(diagonal.Comm());
}

const Epetra_Map& Jacobi::OperatorDomainMap() const
{
  return *range_domain_map;
}

const Epetra_Map& Jacobi::OperatorRangeMap() const
{
  return *range_domain_map;
}

const Epetra_BlockMap& Jacobi::Map() const 
{
  return diagonal.Map();
}

int Jacobi::BlockMap2PointMap(const Epetra_BlockMap* BlockMap, Epetra_Map * & PointMap) const
{
  // Generate an Epetra_Map that has the same number and distribution of points
  // as the input Epetra_BlockMap object.  The global IDs for the output PointMap
  // are computed by using the MaxElementSize of the BlockMap.  For variable block
  // sizes this will create gaps in the GID space, but that is OK for Epetra_Maps.

  int MaxElementSize = BlockMap->MaxElementSize();
  int PtNumMyElements = BlockMap->NumMyPoints();
  int * PtMyGlobalElements = 0;
  if (PtNumMyElements>0) PtMyGlobalElements = new int[PtNumMyElements];

  int NumMyElements = BlockMap->NumMyElements();

  int curID = 0;
  for (int i=0; i<NumMyElements; i++) {
    int StartID = BlockMap->GID(i)*MaxElementSize;
    int ElementSize = BlockMap->ElementSize(i);
    for (int j=0; j<ElementSize; j++) PtMyGlobalElements[curID++] = StartID+j;
  }
  assert(curID==PtNumMyElements); // Sanity test

  PointMap = new Epetra_Map(-1, PtNumMyElements, PtMyGlobalElements, BlockMap->IndexBase(), BlockMap->Comm());

  if (PtNumMyElements>0) delete [] PtMyGlobalElements;

  if (!BlockMap->PointSameAs(*PointMap)) {EPETRA_CHK_ERR(-1);} // Maps not compatible
  return(0);
}
