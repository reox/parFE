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

#ifndef _FEPARAMETERLIST_HPP_
#define _FEPARAMETERLIST_HPP_

#include <Teuchos_ParameterList.hpp>
#include "ProblemReader.h"
#include "ProblemWriter.h"

//! A Class that implement an I/O capable Teuchos::ParameterList 

class FEParameterList : public Teuchos::ParameterList
{
 public:
   FEParameterList(const int NumProc)
   {
     //set default values
     set("#processors", NumProc);
     set("preconditioner", "ml");
     set("load balancing", true);
     set("hdf5 I/O", true);
     set("solve", true);
     set("print for MEDIT", false);
     set("print for PMVIS", false);
     set("print matrix", false);
     set("print lhs", false);
     set("print rhs", false);
     set("convert to ascii", false);
     set("convert to hdf5", false);
     set("verbose", false);
     set("memory usage", std::string("normal"));
     set("element by element", false);
   }
  
  //! Scan a file for relevant parameters by using a ProblemReader
  void Scan(ProblemReader* pr) {
    pr->ScanParameters(*this);
  }

  //! Print relevant parameters into a file by using a ProblemWriter
  void Print(ProblemWriter* pw) {
    pw->PrintParameters(*this);
  }

};

#endif
