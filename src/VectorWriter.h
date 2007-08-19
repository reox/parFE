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

#ifndef _VECTORWRITER_H_
#define _VECTORWRITER_H_

#include <Epetra_Vector.h>

//! A pure virtual class for writing objects of type Epetra_Vector into a file

/*! The VectoWriter class is a pure virtual class (specifies interface only) that enables writing an Epetra_Vector
    object to file by using a particular file format.
*/

class VectorWriter {
 public:
  //! Writes the entire Epetra_Vector into a file.
  
  /*!
     \param vec
     (In) An Epetra_Vector object to be written to file.

     \return Number of records written.
  */
  virtual int PrintVector(const Epetra_Vector& vec)=0;

  virtual ~VectorWriter() {}
};

#endif
