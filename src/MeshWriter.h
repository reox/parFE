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

#ifndef _MESHWRITER_H_
#define _MESHWRITER_H_

#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>

//! A pure virtual class for writing the mesh data into a file.

/*! The MeshWriter class is a pure virtual class (specifies interface only) that enables writing
  mesh data into a file by using a particular file format.
*/  

class MeshWriter {
 public:
  
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
};


#endif
