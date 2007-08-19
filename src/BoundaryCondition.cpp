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

#include "BoundaryCondition.h"
#include <vector>
#include <set>
#include <map>


BoundaryCondition::BoundaryCondition(int nodof, Epetra_Comm& c)
  : num_nodal_dofs(nodof), comm(dynamic_cast<Epetra_MpiComm&>(c))
{  
}


BoundaryCondition::~BoundaryCondition()
{
  delete fixed_nodes;
  delete restrained_nodes;
  delete loaded_nodes;
}


Epetra_Vector* BoundaryCondition::FixedNodes()
{
  return fixed_nodes;
}


Epetra_IntVector* BoundaryCondition::RestrainedNodes()
{
  return restrained_nodes;
}


Epetra_Vector* BoundaryCondition::LoadedNodes()
{
  return loaded_nodes;
}

int BoundaryCondition::NumNodalDofs()
{
  return num_nodal_dofs;
}

//const Epetra_BlockMap& BoundaryCondition::Map()
//{
//  return *node_map;
//}

void BoundaryCondition::Scan(ProblemReader* pr) {
    pr->ScanRestrainedNodes(restrained_nodes, num_nodal_dofs);
    pr->ScanLoadedNodes(loaded_nodes, num_nodal_dofs);
    pr->ScanFixedNodes(fixed_nodes, num_nodal_dofs);
}

void BoundaryCondition::Print(ProblemWriter* pw) {
    pw->PrintRestrainedNodes(*restrained_nodes);
    pw->PrintLoadedNodes(*loaded_nodes);
    pw->PrintFixedNodes(*fixed_nodes);
}

int BoundaryCondition::Redistribute(int num_myglobals, int* my_globals)
{
  RedistributeFixedNodes(num_myglobals, my_globals);
  RedistributeRestrainedNodes(num_myglobals, my_globals);
  RedistributeLoads(num_myglobals, my_globals);
  return 0;
}

int BoundaryCondition::RedistributeFixedNodes(int num_myglobals, int* my_globals)
{
  std::vector<int> node_vec;
  int* PIDList = new int[num_myglobals];
  fixed_nodes->Map().RemoteIDList(num_myglobals, my_globals, PIDList, NULL);
  for (int i=0; i<num_myglobals; ++i)
    if (PIDList[i] >= 0)
      node_vec.push_back(my_globals[i]);

  Epetra_BlockMap new_fixed_map(-1, node_vec.size(), &node_vec[0], num_nodal_dofs, 0, comm);
  Epetra_Export exporter(fixed_nodes->Map(), new_fixed_map);
  Epetra_Vector* new_fixed_nodes = new Epetra_Vector(new_fixed_map);
  new_fixed_nodes->Export(*fixed_nodes, exporter, Insert);
  delete fixed_nodes;
  fixed_nodes = new_fixed_nodes;

  delete[] PIDList;
  return 0;
}


int BoundaryCondition::RedistributeRestrainedNodes(int num_myglobals, int* my_globals)
{
  std::vector<int> node_vec;
  int* PIDList = new int[num_myglobals];
  restrained_nodes->Map().RemoteIDList(num_myglobals, my_globals, PIDList, NULL);

  for (int i=0; i<num_myglobals; ++i)
    if (PIDList[i] >= 0)
      node_vec.push_back(my_globals[i]);

  Epetra_BlockMap new_restrained_map(-1, node_vec.size(), &node_vec[0], num_nodal_dofs, 0, comm);
  Epetra_Export exporter(restrained_nodes->Map(), new_restrained_map);
  Epetra_IntVector* new_restrained_nodes = new Epetra_IntVector(new_restrained_map);
  new_restrained_nodes->Export(*restrained_nodes, exporter, Insert);
  delete restrained_nodes;
  restrained_nodes = new_restrained_nodes;

  delete[] PIDList;
  return 0;
}


int BoundaryCondition::RedistributeLoads(int num_myglobals, int* my_globals)
{
  std::vector<int> node_vec;
  int* PIDList = new int[num_myglobals];
  loaded_nodes->Map().RemoteIDList(num_myglobals, my_globals, PIDList, NULL);

  for (int i=0; i<num_myglobals; ++i)
    if (PIDList[i] >= 0)
      node_vec.push_back(my_globals[i]);

  Epetra_BlockMap new_loaded_map(-1, node_vec.size(), &node_vec[0], num_nodal_dofs, 0, comm);
  Epetra_Export exporter(loaded_nodes->Map(), new_loaded_map);
  Epetra_Vector* new_loaded_nodes = new Epetra_Vector(new_loaded_map);
  new_loaded_nodes->Export(*loaded_nodes, exporter, Insert);
  delete loaded_nodes;
  loaded_nodes = new_loaded_nodes;
  
  delete[] PIDList;
  return 0;
}
