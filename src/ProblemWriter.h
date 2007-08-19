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

#ifndef _PROBLEMWRITER_H_
#define _PROBLEMWRITER_H_

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
  
//! A pure virtual class for writing the data of an elasticity problem into a file.

/*! The VectoWriter class is a pure virtual class (specifies interface only) that enables writing data
    of an elasticity problem to file by using a particular file format.
*/

class ProblemWriter {
 public:
  virtual ~ProblemWriter() {}

  //! Prints the parameters relevant for an elasticity problem into a file.
  /*! The paramters are extracted from a Teuchos::ParameterList.
  \param plist
  (In) A Teuchos::ParameterList that contains all the relevant parameters

  \return Number of records written.
  */
  virtual int PrintParameters(Teuchos::ParameterList plist)=0;
  
  //! Prints the vector containing the coordinates of the nodes of a mesh
   /*!
  \param coordinates
  (In) An Epetra_Vector that contains the coordinates

  \return Number of records written.
  */
  virtual int PrintCoordinates(const Epetra_Vector& coordinates)=0;

  //! Prints the vector forming the connectivity information of the nodes of a mesh
   /*!
  \param elements
  (In) An Epetra_IntVector that contains the node numbers of each element

  \return Number of records written.
  */
  virtual int PrintElements(const Epetra_IntVector& elements)=0;

   //! Print the material ids corresponding to the elements
  /*!
    \param matids
    (In) An Epetra_IntVector that contains the material ids associated with the elements.
    
    \return Number of records read.
  */
  virtual int PrintMaterialIDs(const Epetra_IntVector& matids)=0;
  
  //! Prints the values of the load vector
   /*! The load vector is used to build the right hand side of the system
  \param loaded_nodes
  (In) An Epetra_Vector that values of the load vector

  \return Number of records written.
  */
  virtual int PrintLoadedNodes(const Epetra_Vector& loaded_nodes)=0;
  
  //! Prints the vector containing the nodes that have restrained DOFs assigned.
   /*! A restrained DOF reflects a direction for which the displacement will always remain zero.
  \param restrained_nodes
  (In) An Epetra_IntVector that contains for each node the information about the DOFs restrained

  \return Number of records written.
  */
  virtual int PrintRestrainedNodes(const Epetra_IntVector& restrained_nodes)=0;
  
  //! Prints the vector containing the nodes that have fixed DOFs assigned.
   /*! A fixed DOF reflects a direction for which the displacement will always remain constant, but not zero.
  \param fixed_nodes
  (In) An Epetra_Vector that contains for each node the information about the DOFs fixed.

  \return Number of records written.
  */
  virtual int PrintFixedNodes(const Epetra_Vector& fixed_nodes)=0;

};

#endif
