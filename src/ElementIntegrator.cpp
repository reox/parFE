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

#include "ElementIntegrator.h"
#include "fem.h"
#include <Epetra_MpiComm.h>

ElementIntegrator::ElementIntegrator(int dofs_per_node,
				     int dimension,
				     int num_nodes_per_element,
				     int num_integration_points,
				     int num_material_props,
				     int ss_matrix_size)
  :NumDofsPerNode(dofs_per_node),
   Dimension(dimension),
   NumNodesPerElement(num_nodes_per_element),
   NumDofsPerElement(NumDofsPerNode*NumNodesPerElement),
   NumIntegrationPoints(num_integration_points),
   NumMaterialProps(num_material_props),
   SSMatrixSize(ss_matrix_size)
{

}


const Epetra_SerialDenseMatrix* ElementIntegrator::operator()(double* coord, double* matprops, int num_mattypes)
{

  double* ref_coord = new double[Dimension];

  for (int j=0; j<Dimension; ++j) {
    ref_coord[j] = coord[j];
  }

  for(int i=0; i<NumNodesPerElement; ++i) {
    for (int j=0; j<Dimension; ++j) {
      coord[i*Dimension+j] -= ref_coord[j];
    }
  }

  if (NumNodesPerElement != 8 || NumDofsPerElement != 24 || Dimension != 3) {
    cerr << "The current version of parfe can be used ONLY with:" << endl;
    cerr << "- NumNodesPerElement = 8" << endl;
    cerr << "- NumDofsPerElement = 24" << endl;
    cerr << "- Dimension = 3" << endl;
    exit(EXIT_FAILURE);
  }
      
  km.resize(num_mattypes);
  for (int i=0; i<num_mattypes; ++i) {
    //Macro to call fortran routine
    km[i].Shape(NumDofsPerElement, NumDofsPerElement);
    Stiffness_Matrix(&matprops[i*NumMaterialProps], NumMaterialProps, 
                     NumNodesPerElement, NumDofsPerElement,
                     Dimension, NumIntegrationPoints, 
                     SSMatrixSize, coord, km[i].A());
  }

  delete[] ref_coord;
  
  return &km[0];
}
