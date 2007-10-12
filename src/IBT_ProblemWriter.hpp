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

#ifndef _IBT_PROBLEMWRITER_HPP_
#define _IBT_PROBLEMWRITER_HPP_

#include "ProblemWriter.h"
#include "GWriter.hpp"
#include <Epetra_MpiComm.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_BlockMap.h>
#include <string>
#include "IntMap.h"

//! A class to write an elasticity problem file of the Institute for Biomechanics at ETH Zurich.

/*! The IBT_ProblemReader class is designed to write micro-CT scanner data, provided by the
    institute of biotechnical engineering at the ETH zurich.
    The data file contains the elasticity problem to be solved.

    \param Writer_T
    (In) Perform I/O through this type. Chose one of C_ASCII_GWriter, CPP_ASCII_GWriter or HDF5_GWriter.
*/

template <typename Writer_T>
class IBT_ProblemWriter : public virtual ProblemWriter
{
 public:
   //! ProblemWriter constructor

  /*!
    \param filename
    (In) Path name of the input file.

    \param comm
    (In) The Epetra_MpiComm object associated with this IBT_ProblemReader.

    \return A pointer to the created IBT_ProblemReader object.
  */
  IBT_ProblemWriter(const std::string& filename, Epetra_MpiComm& comm);
  
    //! Prints the parameters relevant for an elasticity problem into a file.
  /*! The paramters are extracted from a Teuchos::ParameterList.
  \param plist
  (In) A Teuchos::ParameterList that contains all the relevant parameters

  \return Number of records written.
  */
  int PrintParameters(Teuchos::ParameterList param);

  //! Prints the vector containing the coordinates of the nodes of a mesh
   /*!
  \param coordinates
  (In) An Epetra_Vector that contains the coordinates

  \return Number of records written.
  */
  int PrintCoordinates(const Epetra_Vector& coordinates);

  //! Prints the vector forming the connectivity information of the nodes of a mesh
   /*!
  \param elements
  (In) An Epetra_IntVector that contains the node numbers of each element

  \return Number of records written.
  */
  int PrintElements(const Epetra_IntVector& elements);
 
  //! Print the material ids corresponding to the elements
  /*!
    \param matids
    (In) An Epetra_IntVector that contains the material ids associated with the elements.
    
    \return Number of records read.
  */
  int PrintMaterialIDs(const Epetra_IntVector& matids);
  
  //! Prints the values of the load vector
   /*! The load vector is used to build the right hand side of the system
  \param loaded_nodes
  (In) An Epetra_Vector that values of the load vector

  \return Number of records written.
  */
  int PrintLoadedNodes(const Epetra_Vector& loaded_nodes);
  
  //! Prints the vector containing the nodes that have restrained DOFs assigned.
   /*! A restrained DOF reflects a direction for which the displacement will always remain zero.
  \param restrained_nodes
  (In) An Epetra_IntVector that contains for each node the information about the DOFs restrained

  \return Number of records written.
  */
  int PrintRestrainedNodes(const Epetra_IntVector& restrained_nodes);
  
  //! Prints the vector containing the nodes that have fixed DOFs assigned.
   /*! A fixed DOF reflects a direction for which the displacement will always remain constant, but not zero.
  \param fixed_nodes
  (In) An Epetra_Vector that contains for each node the information about the DOFs fixed.

  \return Number of records written.
  */
  int PrintFixedNodes(const Epetra_Vector& fixed_nodes);

 private:
  Epetra_MpiComm& comm;
  Writer_T fwriter;

};

template<typename Writer_T> IBT_ProblemWriter<Writer_T>::IBT_ProblemWriter(const std::string& s, Epetra_MpiComm& c)
  : comm(c), fwriter(s, comm.GetMpiComm())
{
}



template<typename Writer_T> int IBT_ProblemWriter<Writer_T>::PrintParameters(Teuchos::ParameterList param)
{
  fwriter.Select("/Parameters");
  fwriter.Write("element type", "'"+param.get("element type","")+"'");
  fwriter.Write("#elements", param.get("#elements",0));
  fwriter.Write("#nodes", param.get("#nodes",0));
  fwriter.Write("#integration points", param.get("#integration points",0));
  fwriter.Write("#dofs per node",param.get("#dofs per node",0));
  fwriter.Write("#nodes per element", param.get("#nodes per element",0));
  fwriter.Write("size of stress-strain matrix", param.get("size of stress-strain matrix",0));
  fwriter.Write("#dimensions", param.get("#dimensions",0));

  fwriter.Write("aa", param.get("aa",0.0));
  fwriter.Write("bb", param.get("bb",0.0));
  fwriter.Write("cc", param.get("cc",0.0));

  int nprops = param.get("#material properties",0);
  int ntypes = param.get("#material types",0);
  fwriter.Write("#material properties", nprops);
  fwriter.Write("#material types", ntypes);
  
  //fwriter.Write("Young's modulus", param.get("Young's modulus",0.0));
  //fwriter.Write("Poisson's ratio", param.get("Poisson's ratio",0.0));
  
  Epetra_SerialDenseMatrix material_props;
  IntMap matid_map;
  matid_map = param.get("material ids", matid_map);
  material_props = param.get("material properties", material_props);
  std::vector<int> material_ids(matid_map.size());
  for (IntMap::const_iterator it = matid_map.begin(); it !=  matid_map.end(); ++it)
    material_ids[(*it).second] = (*it).first;
  if (matid_map.size() > 1)
    fwriter.Write("materials", "ids", &material_ids[0], "properties", material_props.A(), ntypes, ntypes, 1, nprops, 0);
  else
    fwriter.Write("materials", material_props.A(), ntypes, ntypes, nprops, 0);
  fwriter.Write("tolerance", param.get("tolerance",0.0));
  fwriter.Write("iteration limit", param.get("iteration limit",0));

  return 0;
}


template<typename Writer_T> int IBT_ProblemWriter<Writer_T>::PrintCoordinates(const Epetra_Vector& coordinates) {
  //Print nodes
  fwriter.Select("/Mesh");

  const Epetra_BlockMap& coordmap = coordinates.Map();
  //build a linear map
  Epetra_BlockMap linear_map(coordmap.MaxAllGID()+1, coordmap.ElementSize(), 0, comm);

  Epetra_Export linear_exporter(coordmap, linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(coordinates, linear_exporter, Insert);

  return fwriter.Write("Coordinates", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());
}


template<typename Writer_T> int IBT_ProblemWriter<Writer_T>::PrintElements(const Epetra_IntVector& elements)
{

  //Print elements
  fwriter.Select("/Mesh");  
  const Epetra_BlockMap& elementmap = elements.Map();  
  //build a linear map
  Epetra_BlockMap linear_map(elementmap.MaxAllGID()+1, elementmap.ElementSize(), 0, comm);
  Epetra_Export linear_exporter(elementmap, linear_map);
  
  Epetra_IntVector linear_vec(linear_map);
  linear_vec.Export(elements, linear_exporter, Insert);
  
  for (int i=0; i<linear_vec.MyLength(); ++i)
    ++linear_vec[i];

  return fwriter.Write("Elements", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());
}

template<typename Writer_T> int IBT_ProblemWriter<Writer_T>::PrintMaterialIDs(const Epetra_IntVector& matids)
{

  //Print elements
  fwriter.Select("/Mesh");  
  const Epetra_BlockMap& matidmap = matids.Map();  
  //build a linear map
  Epetra_BlockMap linear_map(matidmap.MaxAllGID()+1, matidmap.ElementSize(), 0, comm);
  Epetra_Export linear_exporter(matidmap, linear_map);
  
  Epetra_IntVector linear_vec(linear_map);
  linear_vec.Export(matids, linear_exporter, Insert);

  return fwriter.Write("Material IDs", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());
}

template<typename Writer_T> int IBT_ProblemWriter<Writer_T>::PrintFixedNodes(const Epetra_Vector& fixed_nodes)
{
  //Print fixed nodes
  fwriter.Select("/Boundary conditions");  
  const Epetra_BlockMap& fixedmap = fixed_nodes.Map();
  
  //build a linear map
  Epetra_BlockMap linear_map(fixedmap.MaxAllGID()+1, fixedmap.ElementSize(), 0, comm);
  Epetra_Export linear_exporter(fixedmap, linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(fixed_nodes, linear_exporter, Insert);

  double* it = linear_vec.Values();
  //filter out non-zeros
  std::vector<int> fn_vec;
  std::vector<int> fsns_vec;
  std::vector<double> fval_vec;
  for (int i=0; i<linear_map.NumMyElements(); ++i) {
    for (int j=0; j<linear_map.ElementSize(); ++j) {
      double value = it[j];
      if (value != 0.0) {
	fn_vec.push_back(linear_map.GID(i)+1);
	fsns_vec.push_back(j+1);
	fval_vec.push_back(value);
      }
    }
    it+=linear_map.ElementSize();
  }
  int my_size=fn_vec.size();
  int global_size;
  int my_offset;
  comm.SumAll(&my_size, &global_size, 1);
  comm.ScanSum(&my_size, &my_offset, 1);
  my_offset -= my_size;
  fwriter.Write("Fixed nodes size", global_size);
  return fwriter.Write("Fixed nodes", "Node number", &fn_vec[0], "Sense", &fsns_vec[0], "Value", &fval_vec[0], global_size, my_size, 1,1,1, my_offset);
}


template<typename Writer_T> int IBT_ProblemWriter<Writer_T>::PrintRestrainedNodes(const Epetra_IntVector& restrained_nodes)
{

  //Print restrained nodes
  fwriter.Select("/Boundary conditions");  
  const Epetra_BlockMap& restrainedmap = restrained_nodes.Map();  
  //build a linear map

  Epetra_BlockMap linear_map(restrainedmap.MaxAllGID()+1, restrainedmap.ElementSize(), 0, comm);
  Epetra_Export linear_exporter(restrainedmap, linear_map);

  Epetra_IntVector linear_vec(linear_map);
  linear_vec.PutValue(1);
  linear_vec.Export(restrained_nodes, linear_exporter, Insert);

  int* it = linear_vec.Values();
  
  //filter out non-ones
  std::vector<int> rn_vec;
  std::vector<int> rval_vec;
  for (int i=0; i<linear_map.NumMyElements(); ++i) {
    bool one=true;
    for (int j=0; j<linear_map.ElementSize(); ++j) {
      if (it[j] == 0)
	one = false;
    }
    if (!one) {
      rn_vec.push_back(linear_map.GID(i)+1);
      for (int j=0; j<linear_map.ElementSize(); ++j)
	rval_vec.push_back(it[j]);
    }
    it+=linear_map.ElementSize();
  }
  int element_size = linear_map.ElementSize();
  int my_size=rn_vec.size();
  int global_size;
  int my_offset;
  comm.SumAll(&my_size, &global_size, 1);
  comm.ScanSum(&my_size, &my_offset, 1);
  my_offset -= my_size;
  fwriter.Write("Restrained nodes size", global_size);
  //write into compound data type
  return fwriter.Write("Restrained nodes", "Node number", &rn_vec[0], "Nodal freedom", &rval_vec[0], global_size, my_size, 1, element_size, my_offset);
}


template<typename Writer_T> int IBT_ProblemWriter<Writer_T>::PrintLoadedNodes(const Epetra_Vector& loaded_nodes)
{
  //Print restrained nodes
  fwriter.Select("/Boundary conditions");  
  const Epetra_BlockMap& loadmap = loaded_nodes.Map();  
  //build a linear map
  Epetra_BlockMap linear_map(loadmap.MaxAllGID()+1, loadmap.ElementSize(), 0, comm);
  Epetra_Export linear_exporter(loadmap, linear_map);
  
  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(loaded_nodes, linear_exporter, Insert);
  
  double* it = linear_vec.Values();
  
  //filter out non-zeros
  std::vector<int> ln_vec;
  std::vector<double> lval_vec;
  for (int i=0; i<linear_map.NumMyElements(); ++i) {
    bool zero = true;
    for (int j=0; j<linear_map.ElementSize(); ++j) {
      if (it[j] != 0.0)
	zero = false;
    }
 
    if (!zero) {
      ln_vec.push_back(linear_map.GID(i)+1);
      for (int j=0; j<linear_map.ElementSize(); ++j)
	lval_vec.push_back(it[j]);
    }
    
    it+=linear_map.ElementSize();
  }
  int element_size = linear_map.ElementSize();
  int my_size=ln_vec.size();
  int global_size;
  int my_offset;
  comm.SumAll(&my_size, &global_size, 1);
  comm.ScanSum(&my_size, &my_offset, 1);
  my_offset -= my_size;
  fwriter.Write("Loaded nodes size", global_size);
  //write into compound data type
  return fwriter.Write("Loaded nodes", "Node number", &ln_vec[0], "Loads", &lval_vec[0], global_size, my_size, 1, element_size, my_offset);
}

//TODO: Node and element sets
#endif
