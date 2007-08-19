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

#ifndef _MATLAB_MATRIXWRITER_HPP_
#define _MATLAB_MATRIXWRITER_HPP_

#include "MatrixWriter.h"
#include "GWriter.hpp"
#include <Epetra_MpiComm.h>

//! A class to write a matrix into a MATLAB input file using sparse format.

/*! MATLAB is a numerical computing environment and programming language created by The MathWorks.
    MATLAB_MatrixWriter is a class to output a distributed matrix into a file that can be read by MATLAB.
    The matrix is stored using the MATLAB sparse format.
    If using HDF5 I/O, enter the MATLAB command 'hdf5read' to read the matrix:
    \verbatim >> data = hdf5read(filename, datasetname);\endverbatim
    If using just ASCII I/O, enter the command 'load'. Note that an ASCII file must not end with '.mat'.
    \verbatim >> data = load(filename);\endverbatim

    \param Writer_T
    (In) Perform I/O through this type. Chose one of C_ASCII_GWriter, CPP_ASCII_GWriter or HDF5_GWriter.
*/
template <typename Writer_T>
class MATLAB_MatrixWriter : public MatrixWriter
{
 public:
   
  //! MATLAB_VectorWriter constructor

  /*!
    \param filename
    (In) Path name of the file. Note that it must not end with '.mat' in case of ASCII I/O

    \param comm
    (In) The Epetra_MpiComm object associated with this MATLAB_MatrixWriter.

     \return A pointer to the created MATLAB_MatrixWriter object.
  */
  MATLAB_MatrixWriter(const std::string& filename, Epetra_MpiComm& comm);
  
  //! Writes a Epetra_VbrMatrix containing the data into a file.
  
  /*!
     \param matrix
     (In) An Epetra_VbrMatrix object containing the data to be written to file.
     The matrix will be written in MATLAB sparse format.
     
     \return Number of records written.
  */
  int PrintRows(const Epetra_VbrMatrix* matrix);
 
 private:
  Epetra_MpiComm& comm;
  Writer_T fwriter;

};


template<typename Writer_T> MATLAB_MatrixWriter<Writer_T>::MATLAB_MatrixWriter(const std::string& s, Epetra_MpiComm& c)
  : comm(c), fwriter(s, comm.GetMpiComm())
{
}

template<typename Writer_T> int MATLAB_MatrixWriter<Writer_T>::PrintRows(const  Epetra_VbrMatrix* mat)
{  
  int NumMyRows = mat->NumMyRows();
  int MaxNumRows;
  int my_offset;
  int num_entries;
  int NumGlobalEntries;
  comm.MaxAll(&NumMyRows, &MaxNumRows, 1);
  double* vals;
  int* inds;
  int* rows;
  for (int i=0; i < MaxNumRows; ++i) {
    if (i < NumMyRows) {
      mat->NumMyRowEntries(i, num_entries);
      vals = new double[num_entries];
      inds = new int[num_entries];
      rows = new int[num_entries];
      mat->ExtractMyRowCopy (i, num_entries, num_entries, vals, inds);
      for (int j=0; j < num_entries; ++j) {
	int element_id, offset;
	mat->ColMap().FindLocalElementID(inds[j], element_id, offset);
	inds[j] = mat->GCID(element_id)*mat->ColMap().ElementSize()+offset+1;
	rows[j] = mat->RowMatrixRowMap().GID(i)+1;
      }
    } else {
      num_entries = 0;
      vals=0;
      inds=0;
      rows=0;
    }
    //Write the row
    comm.ScanSum(&num_entries, &my_offset, 1);
    comm.SumAll(&num_entries, &NumGlobalEntries, 1);
    my_offset -= num_entries;
      
    fwriter.Write("Matrix", "Row indices", rows, "Column indices", inds, "Values", vals, NumGlobalEntries, num_entries, 1, 1, 1, my_offset);
    delete[] vals;
    delete[] inds;
    delete[] rows;
    my_offset += num_entries;
  }

  return 0;
}

//Specialization for HDF5
template<> int MATLAB_MatrixWriter<HDF5_GWriter>::PrintRows(const  Epetra_VbrMatrix* mat)
{  
  int my_offset;
  int my_size = mat->NumMyNonzeros();
  int max_size;
  int num_entries;
  comm.SumAll(&my_size, &max_size, 1);
  comm.ScanSum(&my_size, &my_offset, 1);
  my_offset -= my_size;
  double* vals;
  int* inds;
  int* rows;
  for (int i=0; i < max_size; ++i) {
    if (i < my_size ) {
      mat->NumMyRowEntries(i, num_entries);
      vals = new double[num_entries];
      inds = new int[num_entries];
      rows = new int[num_entries];
      mat->ExtractMyRowCopy (i, num_entries, num_entries, vals, inds);
      for (int j=0; j < num_entries; ++j) {
	int element_id, offset;
	mat->ColMap().FindLocalElementID(inds[j], element_id, offset);
	inds[j] = mat->GCID(element_id)*mat->ColMap().ElementSize()+offset+1;
	rows[j] = mat->RowMatrixRowMap().GID(i)+1;
      }
    } else {
      num_entries = 0;
      vals=0;
      inds=0;
      rows=0;
    }
    //Write the row
    fwriter.Write("Matrix", "Row indices", rows, "Column indices", inds, "Values", vals, mat->NumGlobalNonzeros(), num_entries, 1, 1, 1, my_offset);
    delete[] vals;
    delete[] inds;
    delete[] rows;
    my_offset += num_entries;
  }

  return 0;
}



#endif
