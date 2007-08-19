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

#ifndef _PMVIS_MESHWRITER_HPP_
#define _PMVIS_MESHWRITER_HPP_

#include "MeshWriter.h"
#include "GWriter.hpp"
#include <Epetra_MpiComm.h>

//! A class to write the mesh into PMVIS input files.

/*! PMVIS (Partitioned Mesh Visualizer) is a program to view partitioned, unstructured meshes,
    consisting of triangular, quadrilateral, tetrahedral, or hexahedral elements.
    It is available at http://www-users.cs.umn.edu/~oztekin/pmvis/. 
    PMVIS_MeshWriter is a class to output a distributed mesh into files that can be read by PMVIS.
    It creates 3 file: a coordintes file, a connectivity file and a partition file. The filenames
    are specified at construction of a PMVIS_MeshWriter object.
    To invoke PMVIS to visualize the partitioned mesh, enter
    \verbatim > pmvis -n <coordintes file> -c <connectivity file> -p <partition file> -g hex -o 1 \endverbatim
*/
class PMVIS_MeshWriter : public MeshWriter
{
 public:
  //! PMVIS_MeshWriter constructor

  /*!
    \param sxyz
    (In) Path name of the coordinates file.
    
    \param scon
    (In) Path name of the connectivity file.

    \param spart
    (In) Path name of the partition file.

    \param comm
    (In) The Epetra_MpiComm object associated with this PMVIS_MeshWriter.

    \return A pointer to the cretated PMVIS_MeshWriter object.
  */
  PMVIS_MeshWriter(const std::string& sxyz, const std::string& scon, const std::string& spart, Epetra_MpiComm& comm);

  //! Prints the vector containing the coordinates of the nodes of a mesh
   
  /*! The coordinates are written into the coordinate file.
  \param coordinates
  (In) An Epetra_Vector that contains the coordinates

  \return Number of records written.
  */
  int PrintCoordinates(const Epetra_Vector& coordinates);

  //! Prints the vector forming the connectivity information of a mesh
  /*! This generates two files: the connectivity file and the partition file.
      The partition file is generated with repsect to the element distribution.
      
    \param elements
    (In) An Epetra_IntVector that contains the node numbers of each element
    
    \return Number of records written.
  */
  int PrintElements(const Epetra_IntVector& elements);
	
 private:
  Epetra_MpiComm& comm;
  C_ASCII_GWriter fwriter_con;
  C_ASCII_GWriter fwriter_xyz;
  C_ASCII_GWriter fwriter_part;
};


PMVIS_MeshWriter::PMVIS_MeshWriter(const std::string& sxyz, const std::string& scon, const std::string& spart, Epetra_MpiComm& c)
  : comm(c), fwriter_con(scon, comm.GetMpiComm()), fwriter_xyz(sxyz, comm.GetMpiComm()), fwriter_part(spart, comm.GetMpiComm())
{
}


int PMVIS_MeshWriter::PrintCoordinates(const Epetra_Vector& vec) {
  
  //build a linear map
  Epetra_BlockMap linear_map(vec.Map().MaxAllGID()+1, vec.Map().ElementSize(), 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);
  
  fwriter_xyz.Write("Coordinates", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());

  return 0;
}

int PMVIS_MeshWriter::PrintElements(const Epetra_IntVector& vec) {

  //build a linear map
  Epetra_BlockMap linear_map(vec.Map().MaxAllGID()+1, vec.Map().ElementSize(), 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);
  
  Epetra_IntVector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);
  
  int my_size = linear_map.NumMyElements();
  int global_size = linear_map.NumGlobalElements();
  int* PIDList = new int[my_size];
  vec.Map().RemoteIDList (my_size, linear_map.MyGlobalElements(), PIDList, NULL);
  for (int i=0; i<my_size; ++i)
    ++PIDList[i];

  for (int i=0; i<linear_vec.MyLength(); ++i)
    ++linear_vec[i];

  fwriter_con.Write("Elements", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());
 
  fwriter_part.Write("Partition", PIDList, linear_map.NumGlobalElements(), linear_map.NumMyElements(), 1, linear_map.MinMyGID());

  delete[] PIDList;
  return 0;
}

#endif
