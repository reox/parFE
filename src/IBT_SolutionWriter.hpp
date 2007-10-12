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

#ifndef _IBT_SOLUTIONWRITER_HPP_
#define _IBT_SOLUTIONWRITER_HPP_

#include "VectorWriter.h"
#include "GWriter.hpp"
#include "SolutionWriter.h"
#include <Epetra_MpiComm.h>
#include <Epetra_BlockMap.h>
#include <Epetra_Export.h>

//! A class to write the solution of an elasticity problem file of the Institute for Biomechanics at ETH Zurich.

/*! The IBT_SolutionWriter class is designed to write the solution, i.e. the displacements of an elasticity problem
    provided by the Institute for Biomechanics at ETH Zurich.

    \param Writer_T
    (In) Perform I/O through this type. Chose one of C_ASCII_GWriter, CPP_ASCII_GWriter or HDF5_GWriter.
*/

template <typename Writer_T>
class IBT_SolutionWriter : public SolutionWriter
{
 public:
  
  //! IBT_SolutionWriter constructor

   /*!
    \param filename
    (In) Path name of the solution file.

    \param comm
    (In) A Epetra_MpiComm that is associated with this IBT_SolutionWriter.

    \return A pointer to the cretated IBT_SolutionWriter object.
  */
  IBT_SolutionWriter(const std::string& filename, Epetra_MpiComm& comm);

  //! Writes a Epetra_Vector containing the displacements into a file.
  
  /*!
     \param displacements
     (In) An Epetra_Vector object containing the displacements to be written to file.

     \return Number of records written.
  */
  int PrintDisplacements(const Epetra_Vector& displacements);

  int PrintForces(const Epetra_Vector& displacements);
  int PrintStrains(const Epetra_MultiVector& strains);
  int PrintStresses(const Epetra_MultiVector& stresses); 
  
	
 private:
  Epetra_MpiComm& comm;
  Writer_T fwriter;

};


template<typename Writer_T> IBT_SolutionWriter<Writer_T>::IBT_SolutionWriter(const std::string& s, Epetra_MpiComm& c)
  : comm(c), fwriter(s, comm.GetMpiComm())
{
}


//Wrapper function
template<typename Writer_T> int IBT_SolutionWriter<Writer_T>::PrintDisplacements(const Epetra_Vector& vec) {
  
  fwriter.Select("/Solution");
  int dofs_per_node = vec.Map().ElementSize();
  int num_nodes = vec.Map().MaxAllGID()+1;
  
  //Print the displacements
  Epetra_BlockMap linear_map(num_nodes, dofs_per_node, 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);
  int* nnumbers = new int[linear_map.NumMyElements()];
  int* it = linear_map.MyGlobalElements();
  for (int i =0; i<linear_map.NumMyElements(); ++i) {
      nnumbers[i] = (*it)+1;
      ++it;
  }
  fwriter.Write("Nodal displacements", "Node", nnumbers, "Ux Uy Uz", linear_vec.Values(), linear_map.NumGlobalElements(), 1, linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());

  delete[] nnumbers;
  return 0;
}

template<typename Writer_T> int IBT_SolutionWriter<Writer_T>::PrintForces(const Epetra_Vector& vec) {
  
  fwriter.Select("/Solution");
  int dofs_per_node = vec.Map().ElementSize();
  int num_nodes = vec.Map().MaxAllGID()+1;

  
  //Print the forces
  Epetra_BlockMap linear_map(num_nodes, dofs_per_node, 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);

  int* nnumbers = new int[linear_map.NumMyElements()];
  int* it = linear_map.MyGlobalElements();
  for (int i =0; i<linear_map.NumMyElements(); ++i) {
      nnumbers[i] = (*it)+1;
      it++;
  }
  fwriter.Write("Nodal forces", "Node", nnumbers, "Fx Fy Fz", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), 1, linear_map.ElementSize(), linear_map.MinMyGID());

  delete[] nnumbers;
  return 0;
}

template<typename Writer_T> int IBT_SolutionWriter<Writer_T>::PrintStrains(const Epetra_MultiVector& vec) {
  
  fwriter.Select("/Solution");
  int num_strains = vec.Map().ElementSize();
  int num_elements = vec.Map().MaxAllGID()+1;

  //Print the strains
  Epetra_BlockMap linear_map(num_elements, num_strains, 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);

  int* enumbers = new int[linear_map.NumMyElements()];
  int* it = linear_map.MyGlobalElements();
  for (int i =0; i<linear_map.NumMyElements(); ++i) {
      enumbers[i] = (*it)+1;
      it++;
  }
  int *gausspts = new int[linear_map.NumMyElements()];
  for (int i=0; i<vec.NumVectors(); ++i) {

      for (int j=0; j<linear_map.NumMyElements(); ++j)
	  gausspts[j] = i+1;

      fwriter.Write("Element strains", "Element", enumbers, "GaussPt", gausspts, "e11 e22 e33 e12 e23 e31 SED", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), 1, 1, linear_map.ElementSize(), linear_map.MinMyGID()+linear_map.NumGlobalElements()*i);
  }

  delete[] enumbers; 
  delete[] gausspts;
  return 0;
}

template<typename Writer_T> int IBT_SolutionWriter<Writer_T>::PrintStresses(const Epetra_MultiVector& vec) {

    fwriter.Select("/Solution");
    int num_stresses = vec.Map().ElementSize();
    int num_elements = vec.Map().MaxAllGID()+1;
    
    //Print the stresses
    Epetra_BlockMap linear_map(num_elements, num_stresses, 0, comm);
    Epetra_Export linear_exporter(vec.Map(), linear_map);
    
    Epetra_Vector linear_vec(linear_map);
    linear_vec.Export(vec, linear_exporter, Insert);
    
    int* enumbers = new int[linear_map.NumMyElements()];
    int* it = linear_map.MyGlobalElements();
    for (int i =0; i<linear_map.NumMyElements(); ++i) {
	enumbers[i] = (*it)+1;
	it++;
    }
    int *gausspts = new int[linear_map.NumMyElements()];
    for (int i=0; i<vec.NumVectors(); ++i) {
	
	for (int j=0; j<linear_map.NumMyElements(); ++j)
	    gausspts[j] = i+1;
	
	fwriter.Write("Element stresses", "Element", enumbers, "GaussPt", gausspts, "s11 s22 s33 s12 s23 s31 vonMises", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), 1, 1, linear_map.ElementSize(), linear_map.MinMyGID()+linear_map.NumGlobalElements()*i);
    }
    
    delete[] enumbers; 
    delete[] gausspts;
    return 0;
}  

//Specialisations for HDF5:
//Don't print node numbers since they are obvious
template<> int IBT_SolutionWriter<HDF5_GWriter>::PrintDisplacements(const Epetra_Vector& vec) {
  
  fwriter.Select("/Solution");
  int dofs_per_node = vec.Map().ElementSize();
  int num_nodes = vec.Map().MaxAllGID()+1;

  //Print the displacements
  Epetra_BlockMap linear_map(num_nodes, dofs_per_node, 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);

  fwriter.Write("Nodal displacements", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());

  return 0;
}

template<> int IBT_SolutionWriter<HDF5_GWriter>::PrintForces(const Epetra_Vector& vec) {
  
  fwriter.Select("/Solution");
  int dofs_per_node = vec.Map().ElementSize();
  int num_nodes = vec.Map().MaxAllGID()+1;

  
  //Print the forces
  Epetra_BlockMap linear_map(num_nodes, dofs_per_node, 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);

  fwriter.Write("Nodal forces", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID());

  return 0;
}

template<> int IBT_SolutionWriter<HDF5_GWriter>::PrintStrains(const Epetra_MultiVector& vec) {
  
  fwriter.Select("/Solution");
  int num_gausspts = vec.NumVectors();
  fwriter.Write("#Gauss points", num_gausspts);

  int num_strains = vec.Map().ElementSize();
  int num_elements = vec.Map().MaxAllGID()+1;

  //Print the strains
  Epetra_BlockMap linear_map(num_elements, num_strains, 0, comm);
  Epetra_Export linear_exporter(vec.Map(), linear_map);

  Epetra_Vector linear_vec(linear_map);
  linear_vec.Export(vec, linear_exporter, Insert);

  for (int i=0; i<num_gausspts; ++i) {

      fwriter.Write("Element strain", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID()+linear_map.NumGlobalElements()*i);
  }

  return 0;
}

//specialization for HDF5
template<> int IBT_SolutionWriter<HDF5_GWriter>::PrintStresses(const Epetra_MultiVector& vec) {

    fwriter.Select("/Solution");
    int num_gausspts = vec.NumVectors();
    fwriter.Write("#Gauss points", num_gausspts);
    int num_stresses = vec.Map().ElementSize();
    int num_elements = vec.Map().MaxAllGID()+1;
    
    //Print the stresses
    Epetra_BlockMap linear_map(num_elements, num_stresses, 0, comm);
    Epetra_Export linear_exporter(vec.Map(), linear_map);
    
    Epetra_Vector linear_vec(linear_map);
    linear_vec.Export(vec, linear_exporter, Insert);
    
    for (int i=0; i<num_gausspts; ++i) {
	
	fwriter.Write("Element stress", linear_vec.Values(), linear_map.NumGlobalElements(), linear_map.NumMyElements(), linear_map.ElementSize(), linear_map.MinMyGID()+linear_map.NumGlobalElements()*i);
    }
    
    return 0;
}  


#endif
