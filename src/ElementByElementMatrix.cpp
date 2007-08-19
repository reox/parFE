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
#include <set>
#include <algorithm>
#include <Epetra_Map.h>
#include <Epetra_Vector.h>
#include "ElementByElementMatrix.h"

ElementByElementMatrix::ElementByElementMatrix(Epetra_BlockMap& a_dof_map, 
					       ElementIntegrator& elem_mat, 
					       Epetra_IntVector* element_nodes,
					       Epetra_IntVector* mat_ids)
  : dof_map(a_dof_map), ElementMatrix(elem_mat), comm(dynamic_cast<const Epetra_MpiComm&>(a_dof_map.Comm())), Penalty(EBE_PENALTY)
{
  //build point map from blockmap
  BlockMap2PointMap(dof_map, domain_map);
  range_map = domain_map;

  //build importer
  //collect all nodes
  std::set<int> all_nodes_set;
  std::vector<int> all_nodes_vec;
  int node;
  for (int i = 0; i < element_nodes->MyLength(); ++i) {
    node = element_nodes->operator[](i);
    //node is only inserted if it doesn't exist in the set
    all_nodes_set.insert( node );
  }
  all_nodes_vec.insert(all_nodes_vec.begin(), all_nodes_set.begin(), all_nodes_set.end());
  //all_node_map is the map comprising my nodes plus ghost nodes
  Epetra_BlockMap all_node_map(-1, all_nodes_vec.size(), &all_nodes_vec[0], dof_map.ElementSize(), 0, comm);
  Epetra_Map* all_dof_map;
  BlockMap2PointMap(all_node_map, all_dof_map);
  importer = new Epetra_Import(*all_dof_map, *domain_map);
  pmul_vec = new Epetra_Vector(importer->TargetMap());

  int nodof = dof_map.ElementSize();
  ndof = element_nodes->Map().MaxElementSize()*nodof;
  nelem =  element_nodes->Map().NumMyElements();
  g_pmul = new double[ndof * nelem];
  u_tmp = new double[ndof * nelem];
  diagonal = new Epetra_Vector(dof_map);

  //precalculation of time-critical datastructures
  int element, ldof;

  //local indices mapping from owned element dofs plus interface dofs to domain map
  int ptr = 0;
  int matid=0;
  int num_elems=0;
  N = new int[ElementMatrix.Size()];
  inds_ptr = new int[ ElementMatrix.Size()+1];
  for (matid =0; matid < ElementMatrix.Size(); ++matid) {
    inds_ptr[matid] = ptr;
    for (int i = 0 ; i < element_nodes->MyLength(); ++i) {
      int element, ldof;
      element_nodes->Map().FindLocalElementID(i, element, ldof);
      if (mat_ids->operator[](element) == matid) {
	if (ldof==0) ++num_elems;
	int loc = all_dof_map->LID(element_nodes->operator[](i)*nodof);
	for (int j=0; j< nodof; ++j) {
	  gather_inds.push_back(loc + j);
	  ++ptr;
	}
      }
    }
    N[matid] = num_elems;
    num_elems = 0;
  }
  inds_ptr[matid] = ptr;


  //local indices mapping from my element dofs to domain map
  for (matid =0; matid < ElementMatrix.Size(); ++matid) {
    for (int i = 0 ; i < element_nodes->MyLength(); ++i) {
      int element, ldof;
      element_nodes->Map().FindLocalElementID(i, element, ldof);
      if (mat_ids->operator[](element) == matid) {
	int loc = domain_map->LID(element_nodes->operator[](i)*nodof);
	if (loc > -1) {
	  for (int j=0; j< nodof; ++j)
	    scatter_inds.push_back(loc + j);
	} else {
	  for (int j=0; j< nodof; ++j)
	    scatter_inds.push_back(loc);
	}
      }
    }
  }

  //Build diagonal (no communication needed)
  for (int i=0; i<element_nodes->MyLength(); ++i) {
    int gdof = element_nodes->operator[](i);
    if (dof_map.MyGID(gdof)) {
      int element, ldof;
      element_nodes->Map().FindLocalElementID(i, element, ldof);
      ldof = ldof*nodof;
      for (int j=0; j< nodof; ++j) {
	double value = ElementMatrix[mat_ids->operator[](element)](ldof,ldof);
	diagonal->SumIntoGlobalValue(gdof, j, 0, value);
      }
    }
  }

}

int ElementByElementMatrix::PrepareApply(Epetra_Vector* fixed_dofs, Epetra_IntVector* restrained_dofs)
{
  fixed_domain_inds.resize(0);
  //precalculation of time-critical datastructures
  int element, ldof;
  int nodof = dof_map.ElementSize(); 
  //local indices mapping from fixed dofs to domain map
  for (int i = 0 ; i < fixed_dofs->MyLength(); ++i) {
    if (fixed_dofs->operator[](i) != 0.0) {
      fixed_dofs->Map().FindLocalElementID(i, element, ldof);
      fixed_domain_inds.push_back( domain_map->LID(fixed_dofs->Map().GID(element)*nodof + ldof) );
    }
  }
  restrained_domain_inds.resize(0);
  //local indices mapping from restrained dofs to domain map
  for (int i = 0 ; i < restrained_dofs->MyLength(); ++i) {
    if (restrained_dofs->operator[](i) == 0) {
      restrained_dofs->Map().FindLocalElementID(i, element, ldof);
      restrained_domain_inds.push_back( domain_map->LID(restrained_dofs->Map().GID(element)*nodof + ldof) );
    }
  }

  return 0;
}


ElementByElementMatrix::~ElementByElementMatrix()
{
  delete importer;
  delete pmul_vec;
  delete[] g_pmul;
  delete[] u_tmp;
  delete diagonal;
  delete[] inds_ptr;
  delete[] N;
}

int ElementByElementMatrix::SetUseTranspose(bool UseTranspose)
{
  //EPETRA_CHK_ERR(-1);
  return -1;
}

//! Does not impose any BC before and after the matrix-vector product. 
int ElementByElementMatrix::Apply_NoResetBoundaries(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  const Epetra_Vector* X_vector= dynamic_cast<const Epetra_Vector*>(&X);
  Epetra_Vector* Y_vector = dynamic_cast<Epetra_Vector*>(&Y);

  assert (X_vector != 0 && Y_vector != 0);

  pmul_vec->Import(*X_vector, *importer, Insert);

  Y_vector->PutScalar(0.0);

  for (int i=0; i<ElementMatrix.Size(); ++i) {
    //Use creation in View mode, so hope no heap allocation is done here
    Epetra_SerialDenseMatrix PMatrix(View, g_pmul, ElementMatrix[i].LDA(), ElementMatrix[i].M(), N[i]);
    Epetra_SerialDenseMatrix UMatrix(View, u_tmp, ElementMatrix[i].LDA(), ElementMatrix[i].M(), N[i]);
    gather2(pmul_vec, i);
    UMatrix.Multiply('N', 'N', 1.0, ElementMatrix[i], PMatrix, 0.0);
    scatter(Y_vector, i);
  }

  return(0);
}


void ElementByElementMatrix::SetBoundaries(Epetra_Vector* y, const Epetra_Vector* x) const
{
  std::vector<int>::const_iterator it;
  for (it = fixed_domain_inds.begin(); it != fixed_domain_inds.end(); ++it) {
    y->operator[](*it) = x->operator[](*it);
  }

  for (it = restrained_domain_inds.begin(); it != restrained_domain_inds.end(); ++it) {
    y->operator[](*it) = 0.0;
  }

}

//! Sets all boundary nodes (both fixed and restrained) to zero.
void ElementByElementMatrix::ResetBoundaries(Epetra_Vector* y) const
{
  std::vector<int>::const_iterator it;
  for (it = fixed_domain_inds.begin(); it != fixed_domain_inds.end(); ++it) {
    y->operator[](*it) = 0.0;
  }

  for (it = restrained_domain_inds.begin(); it != restrained_domain_inds.end(); ++it) {
    y->operator[](*it) = 0.0;
  }
}

void ElementByElementMatrix::ResetRestrainedBoundaries(Epetra_Vector* y) const
{
  std::vector<int>::const_iterator it;

  for (it = restrained_domain_inds.begin(); it != restrained_domain_inds.end(); ++it) {
    y->operator[](*it) = 0.0;
  }
}

void ElementByElementMatrix::ResetFixedBoundaries(Epetra_Vector* y) const
{
  std::vector<int>::const_iterator it;
  for (it = fixed_domain_inds.begin(); it != fixed_domain_inds.end(); ++it) {
    y->operator[](*it) = 0.0;
  }
}

int ElementByElementMatrix::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  //EPETRA_CHK_ERR(-1);
  return -1;
}

double ElementByElementMatrix::NormInf() const
{
  return(-1.0);
}

const char* ElementByElementMatrix::Label() const
{
  return("ElementByElementMatrix");
}

bool ElementByElementMatrix::UseTranspose() const
{
  return(false);
}

bool ElementByElementMatrix::HasNormInf() const
{
  return(false);
}

const Epetra_Comm& ElementByElementMatrix::Comm() const
{
  return(comm);
}

const Epetra_Map& ElementByElementMatrix::OperatorDomainMap() const
{
  return(*domain_map);
}

const Epetra_Map & ElementByElementMatrix::OperatorRangeMap() const
{
  return(*range_map);
}

const Epetra_BlockMap& ElementByElementMatrix::Map() const 
{
  return(dof_map);
}


void ElementByElementMatrix::gather( const Epetra_Vector* p, int matid ) const {

  int i=0;
  std::vector<int>::const_iterator it;
  for (it = gather_inds.begin()+inds_ptr[matid]; it != gather_inds.begin()+inds_ptr[matid+1]; ++it)
    g_pmul[i++] = pmul_vec->operator[](*it);
}

void ElementByElementMatrix::gather2( const Epetra_Vector* p, int matid ) const {

  const double* p_ptr = p->Values();
  int i=0;
  std::vector<int>::const_iterator it;
  for (it = gather_inds.begin()+inds_ptr[matid]; it != gather_inds.begin()+inds_ptr[matid+1]; ++it)
    g_pmul[i++] = p_ptr[*it]; // FIXME: it was pmul_vec
}


void ElementByElementMatrix::scatter( Epetra_Vector* u, int matid ) const  {

  double* u_ptr = u->Values();

  int i=0;
  std::vector<int>::const_iterator it;
  for (it = scatter_inds.begin()+inds_ptr[matid]; it != scatter_inds.begin()+inds_ptr[matid+1] ; ++it) {
    if (*it > -1)
      u_ptr[*it] += u_tmp[i];
    ++i;
  }
  
  //No communication necessary!
} 


int ElementByElementMatrix::ExtractDiagonalCopy(Epetra_Vector& Diagonal) const
{  
  //Copy by assignment
  Diagonal = *diagonal;
  return 0;
}


int ElementByElementMatrix::ReplaceDiagonalValues(Epetra_Vector& Diagonal) const
{
  //Copy by assignment
  *diagonal = Diagonal;
  return 0;
}


int ElementByElementMatrix::BlockMap2PointMap(const Epetra_BlockMap& BlockMap, Epetra_Map * & PointMap) const
{
  // Generate an Epetra_Map that has the same number and distribution of points
  // as the input Epetra_BlockMap object.  The global IDs for the output PointMap
  // are computed by using the MaxElementSize of the BlockMap.  For variable block
  // sizes this will create gaps in the GID space, but that is OK for Epetra_Maps.

  int MaxElementSize = BlockMap.MaxElementSize();
  int PtNumMyElements = BlockMap.NumMyPoints();
  int * PtMyGlobalElements = 0;
  if (PtNumMyElements>0) PtMyGlobalElements = new int[PtNumMyElements];

  int NumMyElements = BlockMap.NumMyElements();

  int curID = 0;
  for (int i=0; i<NumMyElements; i++) {
    int StartID = BlockMap.GID(i)*MaxElementSize;
    int ElementSize = BlockMap.ElementSize(i);
    for (int j=0; j<ElementSize; j++) PtMyGlobalElements[curID++] = StartID+j;
  }
  assert(curID==PtNumMyElements); // Sanity test

  PointMap = new Epetra_Map(-1, PtNumMyElements, PtMyGlobalElements, BlockMap.IndexBase(), BlockMap.Comm());

  if (PtNumMyElements>0) delete [] PtMyGlobalElements;

  if (!BlockMap.PointSameAs(*PointMap)) {EPETRA_CHK_ERR(-1);} // Maps not compatible
  return(0);
}


int ElementByElementMatrix::Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{

  Y.PutScalar(0.0);

  const Epetra_Vector* X_vector = dynamic_cast<const Epetra_Vector*>(&X);
  Epetra_Vector* Y_vector = dynamic_cast<Epetra_Vector*>(&Y);

  if (X_vector != 0 && Y_vector != 0)
  {
    // in this case they are Epetra_Vector's
    // FIXME: one additional allocation here.
    Epetra_Vector X2_vector(*X_vector);
    ResetBoundaries(&X2_vector);
    pmul_vec->Import(X2_vector, *importer, Insert);

    for (int i=0; i<ElementMatrix.Size(); ++i) {
      //Use creation in View mode, so hope no heap allocation is done here
      Epetra_SerialDenseMatrix PMatrix(View, g_pmul, ElementMatrix[i].LDA(), ElementMatrix[i].M(), N[i]);
      Epetra_SerialDenseMatrix UMatrix(View, u_tmp, ElementMatrix[i].LDA(), ElementMatrix[i].M(), N[i]);
      gather2(pmul_vec, i);
      UMatrix.Multiply('N', 'N', 1.0, ElementMatrix[i], PMatrix, 0.0);
      scatter(Y_vector, i);
    }

    // simply passes the values of X_vector to U, for all matrix rows corresponding
    // to any BC. The resulting matrix is no longer singular.
    SetBoundaries(Y_vector, X_vector);
  }
  else
  {
    assert (X.NumVectors() == Y.NumVectors());

    // Generic case, they are Epetra_MultiVector's.
    //
    // It appears that the operator() does not work properly,
    // at least using INTEL compilers, so  I create a new vector
    // with View properties, and use it. Ideally, only one
    // matrix-vector product will be performed, instead of applying the
    // operator one vector at-a-time.

    for (int v = 0; v < X.NumVectors(); ++v)
    {
      const Epetra_Vector X_vector(View, X.Map(), X[v]); // cheap
      Epetra_Vector X2_vector(Copy, X.Map(), X[v]); // more expensive
      ResetBoundaries(&X2_vector);
      pmul_vec->Import(X2_vector, *importer, Insert);
      Epetra_Vector Y_vector(View, X.Map(), Y[v]);

      for (int i=0; i<ElementMatrix.Size(); ++i) {
        //Use creation in View mode, so hope no heap allocation is done here
        Epetra_SerialDenseMatrix PMatrix(View, g_pmul, ElementMatrix[i].LDA(), ElementMatrix[i].M(), N[i]);
        Epetra_SerialDenseMatrix UMatrix(View, u_tmp, ElementMatrix[i].LDA(), ElementMatrix[i].M(), N[i]);
        gather2(pmul_vec, i);
        UMatrix.Multiply('N', 'N', 1.0, ElementMatrix[i], PMatrix, 0.0);
        scatter(&Y_vector, i);
      }

      // simply passes the values of P to U, for all matrix rows corresponding
      // to any BC. The resulting matrix is no longer singular.
      SetBoundaries(&Y_vector, &X_vector);
    }
  }

#if 0
  U->PutScalar(0.0);

  //import my nodes and ghost nodes
  pmul_vec->Import(*P, *importer, Insert);
  for (int i=0; i<ElementMatrix.Size(); ++i) {
    //Use creation in View mode, so hope no heap allocation is done here
    Epetra_SerialDenseMatrix PMatrix(View, g_pmul, ElementMatrix[i].LDA(), ElementMatrix[i].M(), N[i]);
    Epetra_SerialDenseMatrix UMatrix(View, u_tmp, ElementMatrix[i].LDA(), ElementMatrix[i].M(), N[i]);
    gather(P, i);
    UMatrix.Multiply('N', 'N', 1.0, ElementMatrix[i], PMatrix, 0.0);
    scatter(U, i);
  }
  
  SetBoundaries(U, P);
#endif
  return(0);
}
