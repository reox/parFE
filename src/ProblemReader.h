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

#ifndef _PROBLEMREADER_H_
#define _PROBLEMREADER_H_  

#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>

//! A pure virtual class for reding the data of an elasticity problem from a file.

/*! The ProblemReader class is a pure virtual class (specifies interface only) that enables reading data
    of an elasticity problem from a file by using a particular file format.
*/

class ProblemReader {
 public:

   virtual ~ProblemReader() {}

  //! Reads the parameters relevant for an elasticity problem from a file.
  /*! The paramters are extracted from the file entered to the Teuchos::ParameterList.
    \param param
    (Out) A Teuchos::ParameterList that will contains all the relevant parameters.
    
    \return Number of records read.
  */
  virtual int ScanParameters(Teuchos::ParameterList& param)=0;

  //! Reads the coordinates of the nodes of a mesh and stores them in an Epetra_Vector
  /*!
    \param coordinates
    (Out) An Epetra_Vector that will contain the coordinates.
    
    \return Number of records read.
  */
  virtual int ScanCoordinates(Epetra_Vector& coordinates)=0;

  //! Reads the connectivity information of the nodes of a mesh and stores them in an Epetra_IntVector.
  /*!
    \param elements
    (Out) An Epetra_IntVector that will contains the node numbers of each element.
    
    \return Number of records read.
  */
  virtual int ScanElements(Epetra_IntVector& elements)=0;

  //! Reads the material ids corresponding to the elements and stores them in an Epetra_IntVector.
  /*!
    \param matids
    (Out) An Epetra_IntVector that contains the material id associated with the elements.
    
    \return Number of records read.
  */
  virtual int ScanMaterialIDs(Epetra_IntVector& matids)=0;

  //! Reads the values of the load vector and stores them in an Epetra_Vector.
  /*! The load vector is used to build the right hand side of the system.
    \param loaded_nodes
    (Out) An Epetra_Vector pointer that will point to the values of the load vector.
     
    \param num_nodal_dofs
    (In) Number of DOFs per node.
    
    \return Number of records read.
  */
  virtual int ScanLoadedNodes(Epetra_Vector*& loaded_nodes, int num_nodal_dofs)=0;

  //! Reads the nodes that have restrained DOFs assigned and stores them in an Epetra_IntVector.
  /*! A restrained DOF reflects a direction for which the displacement will always remain zero.
    \param restrained_nodes
    (Out) An Epetra_IntVector pointer that will point to the information about the DOFs restrained for each node.

    \param num_nodal_dofs
    (In) Number of DOFs per node.
    
    \return Number of records read.
  */
  virtual int ScanRestrainedNodes(Epetra_IntVector*& restrained_nodes, int num_nodal_dofs)=0;

  //! Reads the nodes that have fixed DOFs assigned and stores them in an Epetra_Vector.
  /*! A fixed DOF reflects a direction for which the displacement will always remain constant, but not zero.
    \param fixed_nodes
    (Out) An Epetra_Vector pointer that will point to the information about the DOFs fixed for each node.
     
    \param num_nodal_dofs
    (In) Number of DOFs per node.
    
    \return Number of records read.
  */
  virtual int ScanFixedNodes(Epetra_Vector*& fixed_nodes, int num_nodal_dofs)=0;
};

#endif
