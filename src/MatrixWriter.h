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

#ifndef _MATRIXWRITER_H_
#define _MATRIXWRITER_H_

#include <Epetra_VbrMatrix.h>

//! A pure virtual class for writing an Epetra_VbrMatrix into a file.

/*! The MatrixWriter class is a pure virtual class (specifies interface only) that enables writing
  matrix data into a file by using a particular file format.
*/  

class MatrixWriter {
 public:

  //! Prints the matrix row-wise 
  /*!
  \param mat
  (In) A pointer to an Epetra_VbrMatrix

  \return Number of records written.
  */
  virtual int PrintRows(const Epetra_VbrMatrix* mat)=0;
};

#endif
