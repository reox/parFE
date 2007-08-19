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

#ifndef _MATLAB_VECTORWRITER_HPP_
#define _MATLAB_VECTORWRITER_HPP_

#include "VectorWriter.h"
#include "GWriter.hpp"
#include <Epetra_MpiComm.h>


//! A class to write a vector into a MATLAB input file

/*! MATLAB is a numerical computing environment and programming language created by The MathWorks.
    MATLAB_VectorWriter is a class to output a distributed vector into a file that can be read by MATLAB.
    If using HDF5 I/O, enter the MATLAB command 'hdf5read' to read the vector:
    \verbatim >> data = hdf5read(filename, datasetname);\endverbatim
    If using just ASCII I/O, enter the command 'load'. Note that an ASCII file must not end with '.mat'.
    \verbatim >> data = load(filename);\endverbatim

    \param Writer_T
    (In) Perform I/O through this type. Chose one of C_ASCII_GWriter, CPP_ASCII_GWriter or HDF5_GWriter.
*/
template <typename Writer_T>
class MATLAB_VectorWriter : public VectorWriter
{
 public:
  
  //! MATLAB_VectorWriter constructor

  /*!
    \param filename
    (In) Path name of the file. Note that it must not end with '.mat' in case of ASCII I/O

    \param comm
    (In) The Epetra_MpiComm object associated with this MATLAB_VectorWriter.

     \return A pointer to the created MATLAB_VectorWriter object.
  */
  MATLAB_VectorWriter(const std::string& filename, Epetra_MpiComm& comm);
 
  //! Writes a Epetra_Vector containing the data into a file.
  
  /*!
     \param data
     (In) An Epetra_Vector object containing the data to be written to file.

     \return Number of records written.
  */
  int PrintVector(const Epetra_Vector& data);
	
 private:
  Epetra_MpiComm& comm;
  Writer_T fwriter;

};


template<typename Writer_T> MATLAB_VectorWriter<Writer_T>::MATLAB_VectorWriter(const std::string& s, Epetra_MpiComm& c)
  : comm(c), fwriter(s, comm.GetMpiComm())
{
}


template<typename Writer_T> int MATLAB_VectorWriter<Writer_T>::PrintVector(const Epetra_Vector& vec) {
  
  int elem_size = vec.Map().ElementSize();
  int num_nodes = vec.Map().MaxAllGID()+1;

  //Print the displacements
  Epetra_BlockMap linear_map(num_nodes, elem_size, 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);

  fwriter.Write("Vector", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());

  return 0;
}

#endif
