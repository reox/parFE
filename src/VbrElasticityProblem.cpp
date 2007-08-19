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
#include "VbrElasticityProblem.h"


VbrElasticityProblem::VbrElasticityProblem(DistMesh& a_mesh, Teuchos::ParameterList& fe_param)
  : ElasticityProblem(a_mesh, fe_param)
{  
  //Create block map with degrees of freedoms
  Epetra_BlockMap row_map(-1, mesh.NodeMap()->NumMyElements(), mesh.NodeMap()->MyGlobalElements(), NumDofsPerNode, 0, comm);

  Epetra_FECrsGraph graph(Copy, row_map, 64);
  mesh.MatrixGraph(graph);
  
  Epetra_VbrMatrix* block_mat = new Epetra_FEVbrMatrix(Copy, graph);

  Epetra_Vector* X = new Epetra_Vector(row_map);
  Epetra_Vector* B = new Epetra_Vector(row_map);
  SetOperator(block_mat);
  SetLHS(X);
  SetRHS(B);
  Assemble();
}

VbrElasticityProblem::~VbrElasticityProblem()
{
  delete GetLHS();
  delete GetRHS();
  delete GetOperator();
}


int VbrElasticityProblem::Assemble()
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

  Epetra_VbrMatrix* block_mat = dynamic_cast<Epetra_VbrMatrix*>(GetOperator());

  Epetra_IntVector* num = mesh.ElementNodes();
  Epetra_IntVector* mat_ids = mesh.MaterialIDs();
  int* local_nodes = new int[NumNodesPerElement];
  for (int i=0; i<NumNodesPerElement; ++i)
    local_nodes[i] = i;
  int* global_nodes = num->Values();

  //the assembly loop begins
  int iel=0;
  for (int i=0; i< num->MyLength(); i+=NumNodesPerElement) {
    int* it = global_nodes+i;
    for (int j = 0; j<NumNodesPerElement; ++j) {
      if (block_mat->MyGRID(it[j])) {
	block_mat->BeginSumIntoGlobalValues(it[j], NumNodesPerElement, it);
	for (int k=0; k< NumNodesPerElement; ++k) {
	  block_mat->SubmitBlockEntry (ElementMatrix[mat_ids->operator[](iel)].A() + (local_nodes[j] + local_nodes[k]*NumDofsPerElement)*NumDofsPerNode, NumDofsPerElement, NumDofsPerNode, NumDofsPerNode);
	}
	block_mat->EndSubmitEntries();
      }
    }
    ++iel;
  }
  
  block_mat->FillComplete();
  block_mat->OptimizeStorage();

  delete[] local_nodes;

  return 0;
}



int VbrElasticityProblem::Impose(BoundaryCondition& bcond)
{

  Epetra_VbrMatrix* block_mat = dynamic_cast<Epetra_VbrMatrix*>(GetOperator());
  Epetra_Vector* X = dynamic_cast<Epetra_Vector*>(GetLHS());
  Epetra_Vector* B = dynamic_cast<Epetra_Vector*>(GetRHS());

  //B->ReplaceMap(block_mat->Map());
  
  // apply loaded nodes on RHS
  Epetra_Export loads_exporter(bcond.LoadedNodes()->Map(), block_mat->Map());
  B->Export(*bcond.LoadedNodes(), loads_exporter, Insert);
  
  
  //apply fixed nodes on Matrix
  int* gids = bcond.FixedNodes()->Map().MyGlobalElements();
  double* values = bcond.FixedNodes()->Values();
  int num_mygids = bcond.FixedNodes()->Map().NumMyElements();
  int num_myvalues = bcond.FixedNodes()->MyLength();


  Epetra_Vector corr(block_mat->Map());
  Epetra_Vector fixed_nodes(block_mat->Map());
  num_mygids = bcond.FixedNodes()->Map().NumMyElements();
  num_myvalues = bcond.FixedNodes()->MyLength();
  gids = bcond.FixedNodes()->Map().MyGlobalElements();
  values = bcond.FixedNodes()->Values();
  
  for (int i=0; i<num_mygids; ++i) {
    for (int j=0; j<NumDofsPerNode; ++j) {
      if (values[j] != 0.0) {
	fixed_nodes.ReplaceGlobalValue(gids[i], j, 0, values[j]);
      }
    }
    values += NumDofsPerNode;
  }
  block_mat->Multiply1(false, fixed_nodes, corr);
  B->Update(-1, corr, 1);
  

  //Build vector to filter out BC
  Epetra_Vector bc_filter(block_mat->Map());
  bc_filter.PutScalar(1.0);

  
  num_mygids = bcond.RestrainedNodes()->Map().NumMyElements();
  num_myvalues = bcond.RestrainedNodes()->MyLength();
  gids = bcond.RestrainedNodes()->Map().MyGlobalElements();
  int* ivalues = bcond.RestrainedNodes()->Values();

  for (int i=0; i<num_mygids; ++i) {
    for (int j=0; j<NumDofsPerNode; ++j) {
      if (ivalues[j] == 0) {
	bc_filter.ReplaceGlobalValue(gids[i], j, 0, 0.0);
      }
    }
    ivalues += NumDofsPerNode;
  }


  num_mygids = bcond.FixedNodes()->Map().NumMyElements();
  num_myvalues = bcond.FixedNodes()->MyLength();
  gids = bcond.FixedNodes()->Map().MyGlobalElements();
  values = bcond.FixedNodes()->Values();

  for (int i=0; i<num_mygids; ++i) {
    for (int j=0; j<NumDofsPerNode; ++j) {
      if (values[j] != 0) {
	bc_filter.ReplaceGlobalValue(gids[i], j, 0, 0.0);
      }
    }
    values += NumDofsPerNode;
  }


  //filter out bc by applying leftscale and rightscale to the Matrix
  block_mat->LeftScale(bc_filter);
  block_mat->RightScale(bc_filter);

  Epetra_Vector diagonal(block_mat->Map());
  block_mat->ExtractDiagonalCopy(diagonal);
    
  SetBoundaryConditions(bcond, &diagonal, B);

  block_mat->ReplaceDiagonalValues(diagonal);

  X->ReplaceMap(block_mat->OperatorDomainMap());
  B->ReplaceMap(block_mat->OperatorRangeMap());

  return 0;
}


int VbrElasticityProblem::Restore()
{
  Epetra_VbrMatrix* block_mat = dynamic_cast<Epetra_VbrMatrix*>(GetOperator());
  
  block_mat->PutScalar(0.0);
  
  //Reassemble matrix, but this time Epetra_VbrMatrix is transformed to
  //local space
  
  Epetra_IntVector* num = mesh.ElementNodes();
  Epetra_IntVector* mat_ids = mesh.MaterialIDs();
  int* local_nodes = new int[NumNodesPerElement];
  for (int i=0; i<NumNodesPerElement; ++i)
    local_nodes[i] = i;
  int* global_nodes = num->Values();
  int* local_column_indices = new int[NumNodesPerElement];
  //the assembly loop begins
  int iel=0;
  for (int i=0; i< num->MyLength(); i+=NumNodesPerElement) {
    int* it = global_nodes+i;
    for (int l = 0; l<NumNodesPerElement; ++l)
      local_column_indices[l] = block_mat->LCID(it[l]);
    for (int j = 0; j<NumNodesPerElement; ++j) {
      if (block_mat->MyGRID(it[j])) {
	block_mat->BeginSumIntoMyValues(block_mat->LRID(it[j]), NumNodesPerElement, local_column_indices);
	for (int k=0; k< NumNodesPerElement; ++k) {
	  block_mat->SubmitBlockEntry (ElementMatrix[mat_ids->operator[](iel)].A() + (local_nodes[j] + local_nodes[k]*NumDofsPerElement)*NumDofsPerNode, NumDofsPerElement, NumDofsPerNode, NumDofsPerNode);
	}
	block_mat->EndSubmitEntries();
      }
    }
    ++iel;
  }
  return 0;
}
