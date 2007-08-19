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

#ifndef _ELEMENT_INTEGRATOR_H_
#define _ELEMENT_INTEGRATOR_H_
#include <Epetra_IntVector.h>
#include <Epetra_Vector.h>
#include <Epetra_SerialDenseMatrix.h>
#include <vector>
#include <Teuchos_ParameterList.hpp>
#include "IntMap.h"

//! A class to perform element stiffness matrix integrations

/*! This class provides functions to integrate an arbitrary locally stored element
    by using the desired material parameters.
    For the case the same element stiffness matrix is required several time, the last
    computed element stiffness matrices are stored and can be quickly accessed.
*/

class ElementIntegrator {
 public:
  
  //! ElementIntegrator constructor
  
  /*! 
    \param DofsPerNode
    (In) Number of nodal degrees of freedom.

    \param Dimension
    (In) Number of spatial dimensions.

    \param NumNodesPerElement
    (In) Number of nodes per element.

    \param NumIntegrationPoints
    (In) Number of gaussian integration points.

    \param NumMaterialProperties
    (In) Number of different material parameters.

    \param SSMatrixSize
    (In) Size of stress-strain matrix.

    \return Pointer to the ElementIntegrator object.
  */
  ElementIntegrator(int DofsPerNode,
		    int Dimension,
		    int NumNodesPerElement,
		    int NumIntergrationPoints,
		    int NumMaterialProperties,
		    int SSMatrixSize);
  
  //! Computes the element stiffness matrix by integration.

  /*!
    \param lid
    (In) Local ID of the element to be integrated.

    \param matid
    (In) ID of the material parameters to be used for each integration.
    If matid == -1, the element is integrated with every material.

    \return Pointer to the computed element stiffness matrix array.
  */
  const Epetra_SerialDenseMatrix* operator()(double* coords, double* matprops, int num_mattypes = 1);
  
  //! Returns the last computed element stiffness matrix which has been intergrated by using the specified material.
  
  /*!
    \param matid
    (In) The material index an element has been intergrated with

    \return the element last integrated with that material
  */
  inline const Epetra_SerialDenseMatrix& operator[](int matid) const {
    return km[matid];
  }

  //! Return the size of the matrices array, it is usually equal to the number of material parameters.
  inline int Size() const {
    return km.size();
  }

 private:
  const int NumDofsPerNode;
  const int Dimension;
  const int NumNodesPerElement;
  const int NumDofsPerElement;
  const int NumIntegrationPoints;
  const int SSMatrixSize;
  const int NumMaterialProps;
  std::vector<Epetra_SerialDenseMatrix> km;

};


#endif
