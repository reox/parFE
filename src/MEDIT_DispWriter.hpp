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

#ifndef _MEDIT_DISPWRITER_HPP_
#define _MEDIT_DISPWRITER_HPP_

#include "VectorWriter.h"
#include "GWriter.hpp"
#include <Epetra_MpiComm.h>

//! A class to write the displacements into MEDIT input files.

/*! MEDIT is a program for scientific visualization. It is designed to visualize results of computations
    involving 2D or 3D meshes.
    It is available at http://www.ann.jussieu.fr/~frey/logiciels/medit.html.
    MEDIT_DispWriter is a class to write the displacements into a MEDIT solution file. The solution file must end with '.bb'.
    The corresponding mesh file must have the same name execept for the ending. It must be '.mesh' in this case.
    Note that is not possible to run MEDIT without any mesh file given.
    To invoke MEDIT to visualize the solution, enter
    \verbatim > medit <mesh file> \endverbatim
    The solution file will be automatically read and included in the visualization.
*/

class MEDIT_DispWriter : public VectorWriter
{
 public:
  
  //! MEDIT_DispWriter constructor

  /*!
    \param solutionfile
    (In) Path name of the solution file. Note that is must end with '.bb'.

    \param comm
    (In) A Epetra_MpiComm that is associated with this MEDIT_DispWriter.

    \param dim
    (In) Number of solution values per node.

    \return A pointer to the cretated MEDIT_DispWriter object.
  */
  MEDIT_DispWriter(const std::string& solutionfile, Epetra_MpiComm& comm, int dim);
  
  //! Writes a Epetra_Vector containing the displacements into a file.
  
  /*!
     \param displacements
     (In) An Epetra_Vector object containing the displacements to be written to file.

     \return Number of records written.
  */
  int PrintVector(const Epetra_Vector& displacements);

 private:
  Epetra_MpiComm& comm;
  C_ASCII_GWriter fwriter;
  int dim_;

};


MEDIT_DispWriter::MEDIT_DispWriter(const std::string& s, Epetra_MpiComm& c, int dim)
    : comm(c), fwriter(s, comm.GetMpiComm()), dim_(dim)
{}


int MEDIT_DispWriter::PrintVector(const Epetra_Vector& vec) {
  
  int dofs_per_node = vec.Map().ElementSize();
  int num_nodes = vec.Map().MaxAllGID()+1;
  //write header
  int header[4];
  header[0] = dim_;
  header[1] = dofs_per_node;
  header[2] = num_nodes;
  header[3] = 2;

  int i = (comm.MyPID()==0);
  fwriter.Write("Header", header, 1, 1, 4, 0);

  //Print the displacements
  Epetra_BlockMap linear_map(num_nodes, dofs_per_node, 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);

  fwriter.Write("Displacements", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());

  return 0;
}

#endif
