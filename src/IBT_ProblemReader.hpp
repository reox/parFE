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

#ifndef _IBT_PROBLEMREADER_HPP_
#define _IBT_PROBLEMREADER_HPP_

#include "ProblemReader.h"
#include "GReader.hpp"
#include <Epetra_MpiComm.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_Map.h>
#include <vector>
#include <set>
#include "IntMap.h"

//! A class to read an elasticity problem file of the Institute for Biomechanics at ETH Zurich.

/*! The IBT_ProblemReader class is designed to read micro-CT scanner data, provided by the
    institute of biotechnical engineering at the ETH zurich.
    The data file contains the elasticity problem to be solved.

    \param Reader_T
    (In) Perform I/O through this type. Chose one of C_ASCII_GWriter, CPP_ASCII_GWriter or HDF5_GWriter.
*/

template <typename Reader_T>
class IBT_ProblemReader : public ProblemReader
{
 public:
  //! ProblemReader constructor

  /*!
    \param filename
    (In) Path name of the input file.

    \param comm
    (In) The Epetra_MpiComm object associated with this IBT_ProblemReader.

    \return A pointer to the created IBT_ProblemReader object.
    
    \throw std::string
    error string
  */
  IBT_ProblemReader(const std::string& filename, Epetra_MpiComm& comm);

  //! Reads the elasticity problem parameters from the input file.
  /*! The paramters are extracted from the file and entered into a Teuchos::ParameterList. 
    \param param
    (Out) A Teuchos::ParameterList that will contains all the relevant parameters.
    
    \return Number of records read.
  */
  int ScanParameters(Teuchos::ParameterList& param);
  
   //! Reads the coordinates and stores them in an Epetra_Vector
  /*!
    \param coordinates
    (Out) An Epetra_Vector that will contain the coordinates.
    
    \return Number of records read.
  */
  int ScanCoordinates(Epetra_Vector& coordinates);
  
  //! Reads the elements and stores them in an Epetra_IntVector.
  /*!
    \param elements
    (Out) An Epetra_IntVector that will contains the node numbers of each element.
    
    \return Number of records read.
  */
  int ScanElements(Epetra_IntVector& elements);

  //! Reads the material ids corresponding to the elements and stores them in an Epetra_IntVector.
  /*!
    \param matids
    (Out) An Epetra_IntVector that contains the material id associated with the elements.
    
    \return Number of records read.
  */
  int ScanMaterialIDs(Epetra_IntVector& matids);

  //! Reads the values of the load vector and stores them in an Epetra_Vector.
  /*! The load vector is used to build the right hand side of the system.
    \param loaded_nodes
    (Out) An Epetra_Vector pointer that will point to the values of the load vector.
     
    \param num_nodal_dofs
    (In) Number of DOFs per node.
    
    \return Number of records read.
  */
  int ScanLoadedNodes(Epetra_Vector*& loaded_nodes, int num_nodal_dofs);
  
  //! Reads the nodes that have restrained DOFs assigned and stores them in an Epetra_IntVector.
  /*! A restrained DOF reflects a direction for which the displacement will always remain zero.
    \param restrained_nodes
    (Out) An Epetra_IntVector pointer that will point to the information about the DOFs restrained for each node.

    \param num_nodal_dofs
    (In) Number of DOFs per node.
    
    \return Number of records read.
  */
  int ScanRestrainedNodes(Epetra_IntVector*& restrained_nodes, int num_nodal_dofs);
  
  //! Reads the nodes that have fixed DOFs assigned and stores them in an Epetra_Vector.
  /*! A fixed DOF reflects a direction for which the displacement will always remain constant, but not zero.
    \param fixed_nodes
    (Out) An Epetra_Vector pointer that will point to the information about the DOFs fixed for each node.
     
    \param num_nodal_dofs
    (In) Number of DOFs per node.
    
    \return Number of records read.
  */
  int ScanFixedNodes(Epetra_Vector*& fixed_nodes, int num_nodal_dofs);

 private:
  Epetra_MpiComm& comm;
  Reader_T freader;

};


template<typename Reader_T> IBT_ProblemReader<Reader_T>::IBT_ProblemReader(const std::string& s, Epetra_MpiComm& c)
  : comm(c), freader(s, comm.GetMpiComm())
{
  if (!freader) throw "Unable to open file "+s;
}

template<typename Reader_T> int IBT_ProblemReader<Reader_T>::ScanParameters(Teuchos::ParameterList& param) {
  int nels,nn,nip,nodof,nod,nst,ndim,nprops,ntypes,limit,ndof;
  double aa,bb,cc,tol;
  std::string element;
  
  freader.Select("/Parameters");
  freader.Skip('#');
  freader.Read("element type", element);
  freader.Skip('#');
  freader.Read("#elements", nels);
  freader.Skip('#');
  freader.Read("#nodes", nn);
  freader.Skip('#');
  freader.Read("#integration points", nip);
  freader.Skip('#');
  freader.Read("#dofs per node", nodof);
  freader.Skip('#');
  freader.Read("#nodes per element", nod);
  freader.Skip('#');
  freader.Read("size of stress-strain matrix", nst);
  freader.Skip('#');
  freader.Read("#dimensions", ndim);
  freader.Skip('#');
  element = element.substr(1, element.length()-2);
  freader.Skip('#');
  freader.Read("aa", aa);
  freader.Read("bb", bb);
  freader.Read("cc", cc);

  freader.Skip('#');
  freader.Read("#material properties", nprops);
  freader.Read("#material types", ntypes);
  freader.Skip('#');  
  
  Epetra_SerialDenseMatrix MaterialProps(nprops, ntypes);
  IntMap MaterialIDs;
  int* material_ids = new int[ntypes];
  material_ids[0] = 0;
  if (ntypes > 1) {
    freader.Read("materials", "ids", material_ids, "properties", MaterialProps.A(), ntypes, 1, nprops);
  } else {
    freader.Read("materials", MaterialProps.A(), ntypes, nprops);
  }
  for (int i=0; i<ntypes; ++i)
    MaterialIDs[material_ids[i]] = i;

  freader.Skip('#');
  freader.Read("tolerance", tol);
  freader.Read("iteration limit", limit);

  ndof = nod*nodof;

  param.set("material ids", MaterialIDs);
  param.set("material properties", MaterialProps);
  //use get to prevent values to be overridden
  param.get("element type", element);
  param.get("#elements", nels);
  param.get("#nodes", nn);
  param.get("#nodes per element", nod);
  param.get("#dofs per node", nodof);
  param.get("#dofs per element", ndof);
  param.get("#integration points", nip);
  param.get("#dimensions", ndim);
  param.get("size of stress-strain matrix", nst);
  param.get("iteration limit", limit);
  param.get("tolerance", tol);
  param.get("#material properties", nprops);
  param.get("#material types", ntypes);
  //param.get("Young's modulus", E);
  //param.get("Poisson's ratio", v);
  param.get("aa", aa);
  param.get("bb", bb);
  param.get("cc", cc);

  delete[] material_ids;
  
  return 0;
}

template<typename Reader_T> int IBT_ProblemReader<Reader_T>::ScanCoordinates(Epetra_Vector& coordinates) {
  //Scan nodes sequentially
  freader.Skip('#');  
  freader.Select("/Mesh");

  const Epetra_BlockMap& coordmap = coordinates.Map();
  if (!coordmap.LinearMap())
    return -1;
  return freader.Read("Coordinates", coordinates.Values(), coordmap.NumGlobalElements(), coordmap.NumMyElements(), coordmap.ElementSize(), coordmap.MinMyGID());
}


template<typename Reader_T> int IBT_ProblemReader<Reader_T>::ScanElements(Epetra_IntVector& elements)
{
  //Scan elements sequentially
  freader.Skip('#');
  freader.Select("/Mesh");

  const Epetra_BlockMap& elementmap = elements.Map();
  if (!elementmap.LinearMap())
    return -1;
  int res = freader.Read("Elements", elements.Values(), elementmap.NumGlobalElements(), elementmap.NumMyElements(), elementmap.ElementSize(), elementmap.MinMyGID());
  for (int i=0; i<elements.MyLength(); ++i)
    --elements[i];
  return res;
}

template<typename Reader_T> int IBT_ProblemReader<Reader_T>::ScanMaterialIDs(Epetra_IntVector& matids)
{
  //Scan nodes sequentially
  freader.Skip('#');  
  freader.Select("/Mesh");

  const Epetra_BlockMap& matidmap = matids.Map();
  if (!matidmap.LinearMap())
    return -1;
  return freader.Read("Material IDs", matids.Values(), matidmap.NumGlobalElements(), matidmap.NumMyElements(), matidmap.ElementSize(), matidmap.MinMyGID());
}


template<typename Reader_T> int IBT_ProblemReader<Reader_T>::ScanFixedNodes(Epetra_Vector*& fixed_nodes, int num_nodal_dofs)
{
  int num_fn;
  std::set<int> node_set;
  std::vector<int> node_vec;
  freader.Skip('#');
  freader.Select("/Boundary conditions");
  freader.Read("Fixed nodes size", num_fn);
  int* nodes = new int[num_fn];
  int* senses = new int[num_fn];
  double* values = new double[num_fn];

  Epetra_Map dummy_map(num_fn, 0, comm);
  int my_len = dummy_map.NumMyElements();
  freader.Read("Fixed nodes", "Node number", nodes, "Sense", senses, "Value", values, num_fn, my_len, 1, 1, 1, dummy_map.MinMyGID());
  for (int i=0; i<my_len; ++i) {
    --nodes[i];
    --senses[i];
  }

  //find unique nodes
  node_set.insert(nodes, nodes+my_len);
  node_vec.insert(node_vec.end(), node_set.begin(), node_set.end());

  Epetra_BlockMap fixed_map(-1, node_vec.size(), &node_vec[0], num_nodal_dofs, 0, comm);
  fixed_nodes = new Epetra_Vector(fixed_map);
  double* it = fixed_nodes->Values();
  int* first_points = fixed_map.FirstPointInElementList();
  
  for (int i=0; i<my_len; ++i)
    it[first_points[fixed_map.LID(nodes[i])] + senses[i]] = values[i];
  
  delete[] nodes;
  delete[] senses;
  delete[] values;

  return fixed_map.NumMyElements();
}


template<typename Reader_T> int IBT_ProblemReader<Reader_T>::ScanRestrainedNodes(Epetra_IntVector*& restrained_nodes, int num_nodal_dofs)
{
  int num_rn;
  std::vector<int> node_vec;
  freader.Skip('#');
  freader.Select("/Boundary conditions");
  freader.Read("Restrained nodes size", num_rn);
  int* nodes = new int[num_rn];
  int* freedoms = new int[num_rn*num_nodal_dofs];

  Epetra_Map dummy_map(num_rn, 0, comm);
  int my_len = dummy_map.NumMyElements();
  freader.Read("Restrained nodes", "Node number", nodes, "Nodal freedom", freedoms, num_rn, my_len, 1, num_nodal_dofs, dummy_map.MinMyGID());
  for (int i=0; i<my_len; ++i) {
    --nodes[i];
  }

  Epetra_BlockMap restrained_map(-1, my_len, nodes, num_nodal_dofs, 0, comm);
  restrained_nodes = new Epetra_IntVector(Copy, restrained_map, freedoms);

  delete[] nodes;
  delete[] freedoms;

  return restrained_map.NumMyElements();

}


template<typename Reader_T> int IBT_ProblemReader<Reader_T>::ScanLoadedNodes(Epetra_Vector*& loaded_nodes, int num_nodal_dofs)
{
  int num_ln;
  std::vector<int> node_vec;
  freader.Skip('#');
  freader.Select("/Boundary conditions");
  freader.Read("Loaded nodes size", num_ln);
  int* nodes = new int[num_ln];
  double* loads = new double[num_ln*num_nodal_dofs];

  Epetra_Map dummy_map(num_ln, 0, comm);
  int my_len = dummy_map.NumMyElements();
  freader.Read("Loaded nodes", "Node number", nodes, "Loads", loads, num_ln, my_len, 1, 3, dummy_map.MinMyGID());
  for (int i=0; i<my_len; ++i) {
    --nodes[i];
  }

  Epetra_BlockMap loaded_map(-1, my_len, nodes, num_nodal_dofs, 0, comm);
  loaded_nodes = new Epetra_Vector(Copy, loaded_map, loads);

  delete[] nodes;
  delete[] loads;


  return loaded_map.NumMyElements();

}

#endif
