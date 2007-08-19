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

#ifndef _MEDIT_MESHWRITER_HPP_
#define _MEDIT_MESHWRITER_HPP_

#include "MeshWriter.h"
#include "GWriter.hpp"
#include <Epetra_Export.h>
#include <Epetra_MpiComm.h>

//! A class to write the mesh into MEDIT input files.

/*! MEDIT is a program for scientific visualization. It is designed to visualize results of computations
    involving 2D or 3D meshes.
    It is available at http://www.ann.jussieu.fr/~frey/logiciels/medit.html.
    MEDIT_MeshWriter is a class to write the mesh into a MEDIT mesh file. The mesh file must end with '.mesh'.
    The corresponding solution file must have the same name execept for the ending. It must be '.bb' in this case.
    Note that the solution file is not created. There exist the class MEDIT_DispWriter for just that purpose.
    However, it is possible to view the mesh without any solution given.
    To invoke MEDIT to visualize the mesh, enter
    \verbatim > medit <mesh file> \endverbatim
    If there exist a solution file, it will automatically read and included in the visualization.
*/
class MEDIT_MeshWriter : public MeshWriter
{
 public:
  //! MEDIT_MeshWriter constructor

  /*! Constructs a new MEDIT_MeshWriter and writes the file header

    \param meshfile
    (In) Path name of the mesh file. Note that is must end with '.mesh'.

    \param comm
    (In) A Epetra_MpiComm that is associated with this MEDIT_MeshWriter.

    \return A pointer to the cretated MEDIT_MeshWriter object.
  */
  MEDIT_MeshWriter(const std::string& meshfile, int, Epetra_MpiComm& comm);
  
  //! MEDIT_MeshWriter destructor
  
  /*! Writes the file footer before closing
   */
  ~MEDIT_MeshWriter();

  //! Prints the vector containing the coordinates of the nodes of a mesh into the mesh file.
   
  /*!
  \param coordinates
  (In) An Epetra_Vector that contains the coordinates

  \return Number of records written.
  */
  int PrintCoordinates(const Epetra_Vector& coordinates);

  //! Prints the vector forming the connectivity information of a mesh into the mesh file.
  /*!      
    \param elements
    (In) An Epetra_IntVector that contains the node numbers of each element
    
    \return Number of records written.
  */
  int PrintElements(const Epetra_IntVector& elements);
	
 private:
  Epetra_MpiComm& comm;
  C_ASCII_GWriter fwriter;
};


MEDIT_MeshWriter::MEDIT_MeshWriter(const std::string& s, int dim, Epetra_MpiComm& c)
  : comm(c), fwriter(s, comm.GetMpiComm())
{
  fwriter.Write("Header", "MeshVersionFormatted 1\nDimension");
  fwriter.Write("Dimension", dim);
}


MEDIT_MeshWriter::~MEDIT_MeshWriter()
{
  //Write the footer
  fwriter.Write("Footer", "End");
}

int MEDIT_MeshWriter::PrintCoordinates(const Epetra_Vector& vec) {
  
  fwriter.Write("Vertex type", "Vertices");
  fwriter.Write("Vertices size", vec.Map().MaxAllGID()+1);

  //build a linear map
  Epetra_BlockMap linear_map(vec.Map().MaxAllGID()+1, vec.Map().ElementSize(), 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);
  int my_size = linear_map.NumMyElements();
  int global_size = linear_map.NumGlobalElements();
  int* PIDList = new int[my_size];
  vec.Map().RemoteIDList (my_size, linear_map.MyGlobalElements(), PIDList, NULL);
  
  for (int i=0; i<my_size; ++i)
    ++PIDList[i];
  
  fwriter.Write("Coordinates", "xyz", linear_vec.Values(), "PID", PIDList, linear_map.NumGlobalElements(), my_size, linear_map.ElementSize(), 1, linear_map.MinMyGID());

  delete[] PIDList;
  return 0;
}

int MEDIT_MeshWriter::PrintElements(const Epetra_IntVector& vec) {

  //writer header
  fwriter.Write("Element type", "Hexahedra");
  fwriter.Write("Elements size", vec.Map().MaxAllGID()+1);
  
  //build a linear map
  Epetra_BlockMap linear_map(vec.Map().MaxAllGID()+1, vec.Map().ElementSize(), 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);
  
  Epetra_IntVector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);
  
  std::vector<int> material(linear_map.NumMyElements(),1);
  
  for (int i=0; i<linear_vec.MyLength(); ++i)
    ++linear_vec[i];

  fwriter.Write("Elements", "Node numbers", linear_vec.Values(), "Material", &material[0], linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), 1, linear_map.MinMyGID());

  return 0;
}

#endif
