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

#ifndef _BOUNDARY_CONDITON_H_
#define _BOUNDARY_CONDITON_H_
#include <iostream>
#include "ProblemReader.h"
#include "ProblemWriter.h"
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Export.h>
#include "Epetra_MpiComm.h"



class BoundaryCondition
{

 public:
  BoundaryCondition(int nodof, Epetra_Comm&);
  ~BoundaryCondition();

  Epetra_Vector* FixedNodes();
  Epetra_IntVector* RestrainedNodes();
  Epetra_Vector* LoadedNodes();
  //const Epetra_BlockMap& Map();

  int NumNodalDofs();
    
  void Scan(ProblemReader* pr);
  void Print(ProblemWriter* pw);
 
  int Redistribute(int num_myglobals, int* my_globals);
  int RedistributeFixedNodes(int num_myglobals, int* my_globals);
  int RedistributeRestrainedNodes(int num_myglobals, int* my_globals);
  int RedistributeLoads(int num_myglobals, int* my_globals);

 private:
  int num_nodal_dofs;
  Epetra_Vector* fixed_nodes;
  Epetra_IntVector* restrained_nodes;
  Epetra_Vector* loaded_nodes;

  Epetra_BlockMap* node_map;
  Epetra_MpiComm& comm;
};


#endif
