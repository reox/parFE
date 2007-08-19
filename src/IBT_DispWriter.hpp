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

#ifndef _IBT_DISPWRITER_HPP_
#define _IBT_DISPWRITER_HPP_

#include "VectorWriter.h"
#include "GWriter.hpp"
#include <Epetra_MpiComm.h>

//! A class to write the solution of an elasticity problem file of the institute of biotechnical engineering at the ETH zurich.

/*! The IBT_ProblemReader class is designed to write the solution, i.e. the displacements of an elasticity problem
    provided by the institute of biotechnical engineering at the ETH zurich.

    \param Writer_T
    (In) Perform I/O through this type. Chose one of C_ASCII_GWriter, CPP_ASCII_GWriter or HDF5_GWriter.
*/

template <typename Writer_T>
class IBT_DispWriter : public VectorWriter
{
 public:
  
  //! IBT_DispWriter constructor

   /*!
    \param filename
    (In) Path name of the solution file.

    \param comm
    (In) A Epetra_MpiComm that is associated with this IBT_DispWriter.

    \return A pointer to the cretated IBT_DispWriter object.
  */
  IBT_DispWriter(const std::string& filename, Epetra_MpiComm& comm);

  //! Writes a Epetra_Vector containing the displacements into a file.
  
  /*!
     \param displacements
     (In) An Epetra_Vector object containing the displacements to be written to file.

     \return Number of records written.
  */
  int PrintVector(const Epetra_Vector& displacements);
	
 private:
  Epetra_MpiComm& comm;
  Writer_T fwriter;

};


template<typename Writer_T> IBT_DispWriter<Writer_T>::IBT_DispWriter(const std::string& s, Epetra_MpiComm& c)
  : comm(c), fwriter(s, comm.GetMpiComm())
{
}


//Wrapper function
template<typename Writer_T> int IBT_DispWriter<Writer_T>::PrintVector(const Epetra_Vector& vec) {
  
  fwriter.Select("/Solution");
  int dofs_per_node = vec.Map().ElementSize();
  int num_nodes = vec.Map().MaxAllGID()+1;

  //Print the displacements
  Epetra_BlockMap linear_map(num_nodes, dofs_per_node, 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);
  fwriter.Write("Nodal displacements", "Node", linear_map.MyGlobalElements(), "Ux Uy Uz", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), 1, linear_map.ElementSize(), linear_map.MinMyGID());

  return 0;
}

#endif
