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

#ifndef _DISTMESH_H_
#define _DISTMESH_H_

#include "parfe_ConfigDefs.h"

#ifdef HAVE_MPI 
#include "mpi.h"
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Epetra_BlockMap.h>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_CrsGraph.h>
#include <vector>
#include <string>
#include "MeshWriter.h"
#include "ProblemReader.h"
#include "ProblemWriter.h"

//! A class to store and repartition a FE mesh.

/*! The DistMesh class comprises two fundamental datastructures:
    On one hand the element-to-node table, realized as Epetra_IntVector,
    on the other hand the coordinate table implemented as Epetra_Vector.
    The element-to-node table is composed of blocks of global node numbers making up the elements.
    Thus, the block size equals the number of nodes per element.
    The coordinate table is composed of blocks of d real values, where d denotes the number of dimensions (usually 3).
    Both block data structures are achieved by constructing the vector with Epetra_BlockMaps.
    Note that the mesh is distributed in a way that the nodes are always uniquely assigned to
    the processes, whereas some of the elements may be "shared" by multiple processes, i.e. they are
    replicated in all the processes refering to the adjacent subdomains.
    Neighborhood nodes of a mesh partition are implicitely stored in the element-to-node table:
    A node is a neighbor node, if it appears in the local element-to-node table, but not in the local node map,
    i.e. the map determining the distribution of the coordinates.
    In addition, each element has a material ID assigned. The material IDs are stored in a separate Epetra_IntVector
    that is distributed in exactly the same way than the element-to-node table.
*/

class DistMesh {

 public:
  //! DistMesh constructor
  
  /*!
    \param num_v
    (In) Global number of nodes.
    
    \param dim
    (In) Number of spatial dimensions (usually 3)

    \param num_e
    (In) Global number of elements.

    \param v_per_e
    (In) Number of nodes per elements.
 
    \param  MaterialID_map
    (In) Maps the given material IDs to the internal ones.

    \para comm
    (In) A Epetra_MpiComm object associated with the DistMesh object.

    \return A pointer to the created DistMesh object.
  */
  DistMesh(int num_v, int dim, int num_e, int v_per_e, std::map<int,int>& MaterialID_map, Epetra_MpiComm& comm);
  
  //! DistMesh destructor
  virtual ~DistMesh();
  
  //! Global number of nodes
  const int NumNodes;

  //! Number of spatial dimensions
  const int Dimension;

  //! Number of global elements
  const int NumElements;

  //! Number of nodes per element
  const int NumNodesPerElement;

  //! Print the mesh into a seperate mesh file
  
  /*!
    \param mw
    A pointer to the MeshWriter object used to ouput the mesh.
  */
  void Print(MeshWriter* mw);
  
  //! Print the mesh into a problem file

  /*!
    \param pw
    (In) A pointer to the ProblemWriter object used to ouput the mesh.
  */
  void Print(ProblemWriter* pw);

  //! Scan the problem file for mesh datastructures

  /*!
    \param pr
    (In) A pointer to the ProblemReader object used to read in the mesh.
  */
  void Scan(ProblemReader* pr);

  //! Redistributes nodes and elements to obtain a load balanced system

  /*!
    \param repart
    (In) Enable (true) or disable (false) load balancing.

    \param debug
    (In) Enable (true) or disable (false) detailed console output.

    \return 0 on success.
  */
  int Redistribute(bool repart=true, bool debug=false);

  //! Redistributes the nodes according to user defined map.

  /*!
    \param num_myglobals
    (In) Number of locally stored nodes

    \param my_globals
    (In) Global node numbers of locally stored nodes.

    \return 0 on success.
  */
  int RedistributeNodes(int num_myglobals, int* my_globals);

  //! Redistributes the elements according to user defined map.

  /*!
    \param num_myglobals
    Number of locally stored elements

    \param my_globals
    Global element numbers of locally stored nodes.

    \return 0 on success.
  */
  int RedistributeElements(int num_myglobals, int* my_globals);

  //! Return the node map, i.e. the map determining the distribution of the coordinates.
  const Epetra_BlockMap* NodeMap();
  
  //! Return the element map, i.e. the map determining the distribution of the elements.
  const Epetra_BlockMap* ElementMap();

  //! Return the element-to-node table
  Epetra_IntVector* ElementNodes();

  //! Return the material id table
  Epetra_IntVector* MaterialIDs();
  
  //! Return the coordinate table
  Epetra_Vector* Coordinates();

  //! Return the local index of a reference element, i.e. an element having all the coordinates of the element nodes locally available.
  int ReferenceElement();

  //! Return the Epetra_MpiComm object associated with the DistMesh object.
  Epetra_MpiComm& Comm();

  //! Build the vertex graph

  /*! This function works only for hexahedral elements number by convention of the book:
      Ian M. Smith, Vaughan Griffiths,  Programming the Finite Element Method

      \param graph
      (Out) A Epetra_CrsGraph that will contain the vertex graph

      \return 0 on success.
  */
  int VertexGraph(Epetra_CrsGraph& graph);
  
  // Not yet implemented ...
  int ElementGraph(Epetra_CrsGraph& graph);

  //! Build the matrix graph

  /*!
  \param graph
  (Out) A Epetra_CrsGraph that will contain the matrix graph

  \return 0 on success.
  */
  int MatrixGraph(Epetra_CrsGraph& graph);

 private:
  int Repartition(const Epetra_CrsGraph&, int&, int*&, bool=false);
  int Repartition(const Epetra_Vector&, int&, int*&, bool=false);
  int Repartition(const Epetra_BlockMap&, int* adjacy, int* xadj, int& num_myglobals, int*& my_globals, bool debug);
  //int Repartition(int&, int*&, bool=false);
  int ComputeElementEnvelopeOfElements();
  int ComputeElementEnvelopeOfNodes(const Epetra_BlockMap&, int&, int*&);
  Epetra_IntVector* element_nodes;
  Epetra_IntVector* mat_ids;
  Epetra_Vector* coordinates;
  Epetra_MpiComm& comm;
  std::map<int,int> MaterialIDMap;
  
};

#endif
