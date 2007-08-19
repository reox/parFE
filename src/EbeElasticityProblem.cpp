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

#include <Epetra_FECrsGraph.h>
#include <Epetra_MultiVector.h>
#include "EbeElasticityProblem.h"


EbeElasticityProblem::EbeElasticityProblem(DistMesh& a_mesh, 
				     Teuchos::ParameterList& fe_param)
  : ElasticityProblem(a_mesh, fe_param)
{  
  int ref_element = mesh.ReferenceElement();
  //Integrate elements
  //Collect coordinates
  const Epetra_BlockMap& EMap = mesh.ElementNodes()->Map();
  const Epetra_BlockMap& NMap = mesh.Coordinates()->Map();
  //ref_element = EMap.LID(ref_element);
  double* coords = new double[Dimension*NumNodesPerElement];
  int* nodes = mesh.ElementNodes()->Values()+EMap.FirstPointInElement(ref_element);
  for (int i=0; i < NumNodesPerElement; ++i)
    for (int j=0; j < Dimension; ++j)
      coords[i*Dimension+j] = mesh.Coordinates()->operator[](NMap.FirstPointInElement(NMap.LID(nodes[i]))+j);
  ElementMatrix(coords, MaterialProperties.A(), MaterialProperties.N());

  Epetra_BlockMap row_map(-1, mesh.NodeMap()->NumMyElements(), mesh.NodeMap()->MyGlobalElements(), NumDofsPerNode, 0, comm);
  SetOperator(new ElementByElementMatrix(row_map, ElementMatrix, mesh.ElementNodes(), mesh.MaterialIDs()));
  Epetra_Vector* X = new Epetra_Vector(row_map);
  Epetra_Vector* B = new Epetra_Vector(row_map);

  SetLHS(X);
  SetRHS(B);
}

EbeElasticityProblem::~EbeElasticityProblem()
{
  delete GetLHS();
  delete GetRHS();
  delete GetOperator();
}

int EbeElasticityProblem::Impose(BoundaryCondition& bcond)
{


  ElementByElementMatrix* ebe_mat = dynamic_cast<ElementByElementMatrix*>(GetOperator());
  Epetra_Vector* X = dynamic_cast<Epetra_Vector*>(GetLHS());
  Epetra_Vector* B = dynamic_cast<Epetra_Vector*>(GetRHS());

  //X->ReplaceMap(ebe_mat->Map());
  // the start vector is all zeros except at the eqns with fixed displacements
  // the reason is that we use a weighted residual stopping criterion
  // if the vector was all zeros, the first residual would be zero as well!
  Epetra_Export fn_exporter(bcond.FixedNodes()->Map(), ebe_mat->Map());
  X->Export(*bcond.FixedNodes(), fn_exporter, Insert);
  
  ebe_mat->PrepareApply(bcond.FixedNodes(), bcond.RestrainedNodes());
  //B->ReplaceMap(ebe_mat->Map());
  // apply loaded nodes on RHS
  Epetra_Export loads_exporter(bcond.LoadedNodes()->Map(), ebe_mat->Map());
  B->Export(*bcond.LoadedNodes(), loads_exporter, Insert);
  
  Epetra_Vector diagonal(ebe_mat->Map());
  ebe_mat->ExtractDiagonalCopy(diagonal);
  
  SetBoundaryConditions(bcond, &diagonal, B, 1.0); // ebe_mat->Penalty);
  
  ResidualWeightVector = new Epetra_Vector(ebe_mat->Map());
  ResidualWeightVector->PutScalar(1.0);
  SetBoundaryConditions(bcond, ResidualWeightVector, ResidualWeightVector, 0.0);

  ebe_mat->ReplaceDiagonalValues(diagonal);
  
  X->ReplaceMap(ebe_mat->OperatorDomainMap());
  B->ReplaceMap(ebe_mat->OperatorRangeMap());
  
  return 0;
}


int EbeElasticityProblem::Assemble()
{
  //Nothing to do here
  return 0;
}


int EbeElasticityProblem::Restore()
{
  //Clear boundary conditions
  //Create Null-Map
  Epetra_Map NullMap(0,0,comm);
  Epetra_Vector fixed_nodes(NullMap);
  Epetra_IntVector restrained_nodes(NullMap);
  ElementByElementMatrix* ebe_mat = dynamic_cast<ElementByElementMatrix*>(GetOperator());
  ebe_mat->PrepareApply(&fixed_nodes, &restrained_nodes);

  return 0;
}
