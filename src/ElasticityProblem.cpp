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
#include "ElasticityProblem.h"
#include <set>

ElasticityProblem::ElasticityProblem(DistMesh& a_mesh, 
				     Teuchos::ParameterList& fe_param)
  : mesh(a_mesh),
    comm(dynamic_cast<const Epetra_MpiComm&>(mesh.NodeMap()->Comm())),
    NumDofsPerNode(fe_param.get<int>("#dofs per node")),
    Dimension(fe_param.get<int>("#dimensions")),
    NumNodes(fe_param.get<int>("#nodes")),
    NumNodesPerElement(fe_param.get<int>("#nodes per element")),
    NumDofsPerElement(NumDofsPerNode*NumNodesPerElement),
    MaterialProperties(fe_param.get<Epetra_SerialDenseMatrix>("material properties")),
    ElementMatrix(NumDofsPerNode,
		  Dimension,
		  NumNodesPerElement,
		  fe_param.get<int>("#integration points"),
		  MaterialProperties.LDA(),
		  fe_param.get<int>("size of stress-strain matrix")),
    //E(fe_param.get<double>("Young's modulus")),
    //V(fe_param.get<double>("Poisson's ratio")),
    ResidualWeightVector(0),
    Epetra_LinearProblem()
{  

}


ElasticityProblem::~ElasticityProblem()
{
  delete ResidualWeightVector;
}


int ElasticityProblem::SetBoundaryConditions(BoundaryCondition& bcond, Epetra_Vector* diagonal, Epetra_Vector* rhs, double penalty)
{

  // apply restrained nodes on vec
  int* gids = bcond.RestrainedNodes()->Map().MyGlobalElements();
  int* ivalues = bcond.RestrainedNodes()->Values();
  int num_mygids = bcond.RestrainedNodes()->Map().NumMyElements();
  int num_myvalues = bcond.RestrainedNodes()->MyLength();
  
  for (int i=0; i < num_mygids; ++i) {
    for (int j=0; j < NumDofsPerNode; ++j) {
      if (ivalues[j] == 0) {
	//1.0 for diagonal, 0.0 for rhs
	diagonal->ReplaceGlobalValue(gids[i], j, 0, 1.0);
	rhs->ReplaceGlobalValue(gids[i], j, 0, 0.0);
      }
    }
    ivalues+=NumDofsPerNode;
  }

  //apply fixed nodes on vec
  gids = bcond.FixedNodes()->Map().MyGlobalElements();
  double* values = bcond.FixedNodes()->Values();
  num_mygids = bcond.FixedNodes()->Map().NumMyElements();
  num_myvalues = bcond.FixedNodes()->MyLength();

  for (int i=0; i < num_mygids; ++i) {
    for (int j=0; j < NumDofsPerNode; ++j) {
      if (values[j] != 0) {
	//1.0 for diag, values for RHS
	diagonal->ReplaceGlobalValue(gids[i], j, 0, penalty);
	rhs->ReplaceGlobalValue(gids[i], j, 0, penalty*values[j]);
      }
    }
    values+=NumDofsPerNode;
  }

  return 0;
}



void ElasticityProblem::PrintLHS(VectorWriter* vw) {
    Epetra_Vector& vec = *dynamic_cast<Epetra_Vector*>(GetLHS());
    //std::string s("ElementByElementMatrix");
    vec.ReplaceMap(dynamic_cast<Epetra_SrcDistObject*>(GetOperator())->Map());
    vw->PrintVector(vec);
}

void ElasticityProblem::PrintRHS(VectorWriter* vw) {
    Epetra_Vector& vec = *dynamic_cast<Epetra_Vector*>(GetRHS());
    //std::string s("ElementByElementMatrix");
    vec.ReplaceMap(dynamic_cast<Epetra_SrcDistObject*>(GetOperator())->Map());
    vw->PrintVector(vec);
}
  

void ElasticityProblem::PrintMatrix(MatrixWriter* mw) {
    if (GetOperator()->Label() == std::string("Epetra_VbrMatrix"))
	mw->PrintRows(dynamic_cast<Epetra_VbrMatrix*>(GetOperator()));
}


void ElasticityProblem::PrintSolution(SolutionWriter* sw, Epetra_SerialDenseMatrix& MaterialProps) {

    int NumGaussPoints=1;
    int SSMatrixSize = 6;
    //compute forces by applying a matrix-vector product 
    Epetra_Vector& disp = *dynamic_cast<Epetra_Vector*>(GetLHS());
    Epetra_Vector& force = *dynamic_cast<Epetra_Vector*>(GetRHS());
    disp.ReplaceMap(GetOperator()->OperatorDomainMap());
    force.ReplaceMap(GetOperator()->OperatorRangeMap());
    //Reassemble matrix
    Restore();
    GetOperator()->Apply(disp, force);
    disp.ReplaceMap(dynamic_cast<Epetra_SrcDistObject*>(GetOperator())->Map());
    force.ReplaceMap(dynamic_cast<Epetra_SrcDistObject*>(GetOperator())->Map());

    sw->PrintDisplacements(disp);
    sw->PrintForces(force);

    //Construct vector with coordinates of element nodes
    Epetra_BlockMap coordinate_map(mesh.ElementMap()->NumGlobalElements(),
				   mesh.ElementMap()->NumMyElements(),
				   mesh.ElementMap()->MyGlobalElements(),
				   mesh.NodeMap()->ElementSize()*mesh.ElementMap()->ElementSize(),
				   0, comm);

   
    Epetra_BlockMap displacement_map(mesh.ElementMap()->NumGlobalElements(),
				     mesh.ElementMap()->NumMyElements(),
				     mesh.ElementMap()->MyGlobalElements(),
				     disp.Map().ElementSize()*mesh.ElementMap()->ElementSize(),		     
				     0, comm);
    
    //build importer
    //collect all nodes
    std::set<int> all_nodes_set;
    std::vector<int> all_nodes_vec;
    int node;
    for (int i = 0; i < mesh.ElementNodes()->MyLength(); ++i) {
	node = mesh.ElementNodes()->operator[](i);
	//node is only inserted if it doesn't exist in the set
	all_nodes_set.insert( node );
    }
    all_nodes_vec.insert(all_nodes_vec.begin(), all_nodes_set.begin(), all_nodes_set.end());
    
    Epetra_BlockMap all_coord_map(-1, all_nodes_vec.size(), &all_nodes_vec[0], mesh.NodeMap()->ElementSize(), 0, comm);
    Epetra_BlockMap all_disp_map(-1, all_nodes_vec.size(), &all_nodes_vec[0], disp.Map().ElementSize(), 0, comm);
    
    Epetra_Vector all_disp(all_disp_map);
    Epetra_Vector all_coord(all_coord_map);

    Epetra_Vector coordinates(coordinate_map);
    Epetra_Vector displacements(displacement_map);
    
    Epetra_Import disp_importer(all_disp_map, disp.Map());
    all_disp.Import(disp, disp_importer, Insert);

    Epetra_Import coord_importer(all_coord_map, *mesh.NodeMap());
    all_coord.Import(*mesh.Coordinates(), coord_importer, Insert);
    
    int k_d=0;
    int k_c=0;
    for (int i = 0 ; i < mesh.ElementNodes()->MyLength(); ++i) {
      int loc_d = all_disp_map.LID(mesh.ElementNodes()->operator[](i));
      int p_d = all_disp_map.FirstPointInElement(loc_d);
      int loc_c = all_coord_map.LID(mesh.ElementNodes()->operator[](i));
      int p_c = all_coord_map.FirstPointInElement(loc_c);
      for (int j=0; j< all_disp_map.ElementSize(); ++j) {
	  displacements[k_d++] = all_disp[p_d+j];
      }
      for (int j=0; j< all_coord_map.ElementSize(); ++j) {
	  coordinates[k_c++] = all_coord[p_c+j];
      }
    }
    
    Epetra_BlockMap strain_map(mesh.ElementMap()->NumGlobalElements(), 
		    mesh.ElementMap()->NumMyElements(),
		    mesh.ElementMap()->MyGlobalElements(),
		    (SSMatrixSize +2),
		    0, comm);

    Epetra_BlockMap stress_map(mesh.ElementMap()->NumGlobalElements(), 
		    mesh.ElementMap()->NumMyElements(),
		    mesh.ElementMap()->MyGlobalElements(),
		    (SSMatrixSize +1),
		    0, comm);

    //Use multivector here because of multiple gauss points
    Epetra_MultiVector strains(strain_map, NumGaussPoints);
    Epetra_MultiVector stresses(stress_map, NumGaussPoints);

    Epetra_BlockMap foomap( mesh.ElementMap()->NumGlobalElements(), 
		mesh.ElementMap()->NumMyElements(),
		mesh.ElementMap()->MyGlobalElements(),
		NumGaussPoints,
		0, comm);

    Epetra_Vector theta(foomap);
    Epetra_Vector sigma(foomap);
    
    //collect strains and stresses of each element

    double* strainbuf = new double[(SSMatrixSize +1) * NumGaussPoints ];
    double* stressbuf = new double[(SSMatrixSize +1) * NumGaussPoints ];

    for (int i=0; i<mesh.ElementMap()->NumMyElements(); ++i) {
	Element_Stress(MaterialProps[mesh.MaterialIDs()->operator[](i)],
		       MaterialProps.LDA(), 
		       NumNodesPerElement, NumDofsPerElement, Dimension, NumGaussPoints, SSMatrixSize, 
                       &coordinates[i*Dimension*NumNodesPerElement],
                       &displacements[i*NumDofsPerElement],
		       strainbuf,
		       stressbuf,
		       &sigma[foomap.FirstPointInElement(i)],
		       &theta[foomap.FirstPointInElement(i)]);

	
	for (int j=0; j<NumGaussPoints; ++j) {
	    memcpy(strains[j] + strains.Map().FirstPointInElement(i), strainbuf+j*NumGaussPoints, (SSMatrixSize +1) *  sizeof(double));
	    memcpy(stresses[j] + stresses.Map().FirstPointInElement(i), stressbuf+j*NumGaussPoints, (SSMatrixSize +1) *  sizeof(double));
	    //compute effective strain
	    double* sed = strains[j] + strains.Map().FirstPointInElement(i)+SSMatrixSize;
	    double* eff = sed+1;
	    double E_tissue_FE = MaterialProps[mesh.MaterialIDs()->operator[](i)][0];
	    *eff = std::sqrt(2.* (*sed)/E_tissue_FE);
	}
    }

    delete[] strainbuf;
    delete[] stressbuf;

    sw->PrintStrains(strains);
    sw->PrintStresses(stresses);

    
}
