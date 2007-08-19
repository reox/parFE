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

#include "DistMesh.h"
#include <Epetra_FECrsGraph.h>
#include <Epetra_Import.h>
#include "parmetis.h"
#include <set>
#include <map>
//#include <pat_api.h>

//Macros to print trace information
//#define CALL(mesg)\
//     std::cerr << "Proc " << comm.MyPID() << ": "; \
//     std::cerr << "Calling " << mesg << "..." << std::endl;

//#define RET(mesg, ret_val)\
//     std::cerr << "Proc " << comm.MyPID() << ": "; \
//     std::cerr << mesg << " returned " << ret_val << "." << std::endl;

#define CALL(mesg) ;
#define RET(mesg, ret_val) ;

DistMesh::DistMesh(int num_nodes, int dimension, int num_elements, int nodes_per_element, std::map<int,int>& mat_id_map, Epetra_MpiComm& c)
  : NumNodes(num_nodes),
    Dimension(dimension),
    NumElements(num_elements),
    NumNodesPerElement(nodes_per_element),
    MaterialIDMap(mat_id_map),
    comm(c),
    coordinates(NULL),
    element_nodes(NULL),
    mat_ids(NULL)
{}


DistMesh::~DistMesh()
{
  delete coordinates;
  delete element_nodes;
  delete mat_ids;
}


int DistMesh::VertexGraph(Epetra_CrsGraph& graph) {
  //PAT_region_begin (21, "VertexGraph");
  
  //PAT_region_begin (31, "Init distro in VertexGraph");
  int* my_elements;
  int num_myelements;
  CALL("DistMesh::ComputeElementEnvelopeOfNodes");
  int ret = ComputeElementEnvelopeOfNodes(graph.Map(), num_myelements, my_elements);
  RET("DistMesh::ComputeElementEnvelopeOfNodes", ret);
  //RedistributeElements(num_myelements, my_elements);
  Epetra_BlockMap elementmap(-1, num_myelements, my_elements, NumNodesPerElement, 0, comm);
  delete[] my_elements;
  Epetra_IntVector* new_elements = new Epetra_IntVector(elementmap);
  Epetra_Import elements_importer(elementmap, element_nodes->Map());
  ret = new_elements->Import(*element_nodes, elements_importer, Insert);
  //PAT_region_end(31);

  //PAT_region_begin(32, "build graph in VertexGraph");
  const Epetra_BlockMap& vtx_map = graph.Map();
  int neighbors[3]; //each vertex has three neighbors
  
  int* global_nodes = new_elements->Values();
  for (int ie = 0 ; ie < new_elements->MyLength(); ie+=NumNodesPerElement)
  {
    int* k = global_nodes+ie;
    for (int i = 0; i < NumNodesPerElement; ++i) {
      int j = k[i];
      if (vtx_map.MyGID(j)) {
	// compute neighbors
	neighbors[0] = k[i ^ 1];
	neighbors[1] = k[i ^ 3];
	neighbors[2] = k[i ^ 4];
	graph.InsertGlobalIndices(j, 3, neighbors);
      }
    }
  }
  //PAT_region_end(32);
  //PAT_region_begin(33, "FillComplete() in VertexGraph");
  // Transform to local indices
  ret = graph.FillComplete();
  //PAT_region_end(33);
  // Commented out because of bad performance
  // Arrange elements to build contigous memory blocks
  //if (graph.OptimizeStorage())
  //return -1;
  delete new_elements; 
  //PAT_region_end(21);
  return 0;
}

int DistMesh::MatrixGraph(Epetra_CrsGraph& graph) {
  //PAT_region_begin (20, "MatrixGraph");
  const Epetra_BlockMap& matrix_map = graph.Map();
   
  int* global_nodes = element_nodes->Values();
  for (int ie = 0 ; ie < element_nodes->MyLength(); ie+=NumNodesPerElement)
  {
    int* k = global_nodes+ie;
    for (int i = 0; i < NumNodesPerElement; ++i) {
      int j = k[i];
      if (matrix_map.MyGID(j)) {
	graph.InsertGlobalIndices(j, NumNodesPerElement, k);
      }
    }
  }
  
  // Transform to local indices
  if (graph.FillComplete() != 0)
    return -1;

  // Arrange elements to build contigous memory blocks
  //if (graph.OptimizeStorage() != 0)
  //return -1;
  //PAT_region_end(20);
  return 0;
}


int DistMesh::ComputeElementEnvelopeOfElements()
{
  Epetra_Map elem_map(NumElements, 0, comm);
  Epetra_CrsGraph graph(Copy, elem_map, 32);
  ElementGraph(graph);

  //Get set of elements in the graph
  std::set<int> my_element_set;
  std::vector<int> my_elements;
  for ( int i=0; i<graph.NumMyRows(); ++i ) {
    int len;
    int* indices;
    graph.ExtractMyRowView(i, len, indices);
    my_element_set.insert(indices, indices+len);
  }
  //my_elements.insert(my_elements.end(), my_element_set.begin(), my_element_set.end());
  //TODO: store set of my element nodes
  //TODO: Redistribute Nodes, such that they are unique
  //1. Find ghost elements and the pids of the owners
  //2. Find shared nodes: Intersection of my element nodes and ghost element nodes
  //3. From each sharer, request random value for a specified node (with CreateFromRecvs)
  //4. Choose highest value. If I have highest value, keep the node, else discard.
  //5. ComputeElementEnvelopeOfNodes();
  return RedistributeElements(my_elements.size(), &my_elements[0]);
}


int DistMesh::ComputeElementEnvelopeOfNodes(const Epetra_BlockMap& vtx_map, int& num_myelements, int*& my_elements )
{
    //PAT_region_begin (11, "ComputeElementEnvelopeOfNodes");
  //Redistribute the elements
  //int max = element_nodes->MaxValue();
  //int min = element_nodes->MinValue();
  //int len = max - min + 1;
  //get node gids in the element-to-node table
  std::map<int, int> pid_map;
  //initialize gid_map
  for (int i=0; i<element_nodes->MyLength();++i)
    pid_map[element_nodes->operator[](i)] = -1;

  int len = pid_map.size();
  int* pids = new int[len];
  int* gids = new int[len];
  int* int_it = gids;
  for (std::map<int, int>::iterator it = pid_map.begin(); it != pid_map.end(); ++it)
    *int_it++ = it->first;
  int ret = vtx_map.RemoteIDList(len, gids, pids, 0);
  int_it = pids;
  for (std::map<int, int>::iterator it = pid_map.begin(); it != pid_map.end(); ++it)
     it->second = *int_it++;

  delete[] pids;
  delete[] gids;

  int* my_lids = element_nodes->Map().PointToElementList();
  std::set<int> my_element_set;
  std::set<int> export_PID_set;
  std::vector<int> export_PIDs;
  std::vector<int> export_objs;
  int num_export_IDs = 0;
  int num_remote_IDs = 0;

  len = element_nodes->MyLength();
  int* values =  element_nodes->Values();
  int count = 0;
  //let's examine the PIDList
  for (int i=0; i<len; i+=NumNodesPerElement) {
    //get actual element
    int gid = element_nodes->Map().GID(my_lids[i]);
    for (int j=0; j<NumNodesPerElement; ++j) {

      if (pid_map[values[i+j]] == comm.MyPID()) ++count;
      else export_PID_set.insert(pid_map[values[i+j]]);
    }

    if (count > 0) {     //if some nodes belong to me, I keep this element
      my_element_set.insert(gid);
    }
    if (export_PID_set.size() > 0) { //this is a remote or shared element
      export_PIDs.insert(export_PIDs.end(), export_PID_set.begin(), export_PID_set.end());
      export_objs.insert(export_objs.end(), export_PID_set.size(), gid);
    }
    count = 0;
    export_PID_set.clear();
  }
  
  num_export_IDs = export_PIDs.size();
  Epetra_Distributor* part_dist = comm.CreateDistributor();

  ret = part_dist->CreateFromSends(num_export_IDs, &(*(export_PIDs.begin())), true, num_remote_IDs);
  int len_import_objs = num_remote_IDs*sizeof(int);
  int* import_objs = new int[len_import_objs];
  ret = part_dist->Do( reinterpret_cast<char*>(&(*export_objs.begin())), static_cast<int>(sizeof(int)), len_import_objs, reinterpret_cast<char*&>(import_objs) );
  //add imported nodes to the new map
  int num_import_objs = 0;
  for (int i=0; i<num_remote_IDs; ++i) {
    my_element_set.insert(import_objs[num_import_objs]);
    ++num_import_objs;
  }

  delete part_dist;
  delete[] import_objs;
  num_myelements =  my_element_set.size();
  my_elements = new int[num_myelements];
  int* my_elements_it = my_elements;
  for (std::set<int>::iterator it = my_element_set.begin(); it != my_element_set.end(); ++my_elements_it)
    *my_elements_it = *it++;
  //PAT_region_end(11);
  return num_myelements;
}


int DistMesh::ElementGraph(Epetra_CrsGraph& graph) {
  const Epetra_BlockMap& elem_map = graph.Map();

  //Build ParMETIS datastructures
  // build elmdist array
  int* elmdst = new int[comm.NumProc()+1];
  int myMaxGid = elem_map.MaxMyGID() + 1;
  comm.GatherAll (&myMaxGid, elmdst+1, 1 );
  elmdst[0] = 0;
  
  //build eptr array
  int num_elems = ElementMap()->NumMyElements();
  int* eptr = new int[num_elems+1];
  memcpy(eptr, ElementMap()->FirstPointInElementList(),num_elems*sizeof(int));
  eptr[num_elems] = element_nodes->MyLength();
  
  //build eind array
  int* eind = element_nodes->Values();
  int numflag = 0;
  //common nodes is equal to 1 because every common node results in additional
  //communication
  int ncommonnodes = 1;
  
  MPI_Comm mpi_comm = comm.Comm();
  
  //output
  int* xadj;
  int* adjncy;

  //call
  ParMETIS_V3_Mesh2Dual(elmdst, eptr, eind, &numflag, &ncommonnodes, &xadj, &adjncy, &mpi_comm);

  //cleanup
  delete[] eptr;
  delete[] elmdst;

  //fill graph
  for (int i=0; i<num_elems+1; ++i) {
    int num_neighbors = xadj[i+1] - xadj[i];
    graph.InsertGlobalIndices(ElementMap()->GID(i), num_neighbors, adjncy+xadj[i]);
  }

  // Transform to local indices
  if (graph.FillComplete())
    return -1;

  // Arrange elements to build contigous memory blocks
  if (graph.OptimizeStorage())
    return -1;
 
 //clean up
  delete[] xadj;
  delete[] adjncy;

  return 0;
}


int DistMesh::Redistribute(bool repart, bool debug)
{
#if 1
  // Do it the vertex based way  
  if (repart) {
    int* my_globals;
    int num_myglobals;
    int ret = Repartition(*coordinates, num_myglobals, my_globals, debug);
    Epetra_Map vtx_map(-1, num_myglobals, my_globals, 0, comm);
    Epetra_FECrsGraph graph(Copy, vtx_map, 64);
    ret = VertexGraph(graph);
    delete[] my_globals;
    ret = Repartition(graph, num_myglobals, my_globals, debug);
    //Repartition(num_myglobals, my_globals, debug);
    ret = RedistributeNodes(num_myglobals, my_globals);
    delete[] my_globals;
  }
  int* my_elements;
  int num_myelements;
  CALL("DistMesh::ComputeElementEnvelopeOfNodes");
  int ret = ComputeElementEnvelopeOfNodes(coordinates->Map(), num_myelements, my_elements);
  RET("DistMesh::ComputeElementEnvelopeOfNodes", ret);
  int res = RedistributeElements(num_myelements, my_elements);
  delete[] my_elements;
  return res;
#else
  // Do it the element based way
  if (repart) {
    Epetra_Map elem_map(NumElements, 0, comm);
    Epetra_CrsGraph graph = new Epetra_CrsGraph(Copy, elem_map, 0);
    ElementGraph(*graph);
    int* my_globals;
    int num_myglobals;
    Repartition(graph, num_myglobals, my_globals);
    RedistributeElements(num_myglobals, my_globals);
    delete[] my_globals;
  }
  ComputeElementEnvelopeOfElements();
#endif
}


int DistMesh::RedistributeElements(int num_myglobals, int* my_globals)
{
  //PAT_region_begin(10, "RedistributeElements");
  //TODO: Check:Works only if original element map has unique gids!!
  Epetra_BlockMap elementmap(-1, num_myglobals, my_globals, NumNodesPerElement, 0, comm); 
  //now we have the new element-map. Remains exporting the nodes.
  Epetra_IntVector* new_elements  = new Epetra_IntVector(elementmap);
  Epetra_Import elements_importer(elementmap, element_nodes->Map());
  int ret = new_elements->Import(*element_nodes, elements_importer, Insert);
  //done. Reset pointers to new Vectors
  delete element_nodes;
  element_nodes = new_elements;

  //TODO: Check:Works only if original material id map has unique gids!!
  Epetra_Map matidmap(-1, num_myglobals, my_globals, 0, comm); 
  //now we have the new material id map. Remains exporting ids
  Epetra_IntVector* new_mat_ids  = new Epetra_IntVector(matidmap);
  Epetra_Import matid_importer(matidmap, mat_ids->Map());
  ret = new_mat_ids->Import(*mat_ids, matid_importer, Insert);
  //done. Reset pointers to new Vectors
  delete mat_ids;
  mat_ids = new_mat_ids;
  //PAT_region_end(10);
  return 0;
}

int DistMesh::RedistributeNodes(int num_myglobals, int* my_globals)
{
  //PAT_region_begin (9, "RedistributeNodes");
  //TODO: Check:Works only if original node map has unique gids!!
  //Redistribute the coordinates 
  Epetra_BlockMap coordmap(NumNodes, num_myglobals, my_globals, Dimension, 0, comm);
  Epetra_Vector* new_coords = new Epetra_Vector(coordmap);
  Epetra_Import coords_importer(coordmap, coordinates->Map());
  int ret = new_coords->Import(*coordinates, coords_importer, Insert);
  //done. Reset pointers to new Vectors
  delete coordinates;
  coordinates = new_coords;
  //PAT_region_end(9);
  return 0;
  
}


int DistMesh::Repartition(const Epetra_CrsGraph& graph, int& num_myglobals, int*& my_globals, bool debug)
//int DistMesh::Repartition(int& num_myglobals, int*& my_globals, bool debug)
{
  // //copy the graph map
  const Epetra_BlockMap& graph_map = graph.Map();
  //build datastructures for parmetis
  int num_indices;
  int* xadj = new int[graph.NumMyRows()+1];
  int adjacy_len = graph.NumMyNonzeros();
  int* adjacy = new int[adjacy_len];
  int ind_sum = 0;
  xadj[0] = 0;
  for (int i=0; i<graph.NumMyRows(); ++i) {
    graph.ExtractGlobalRowCopy(graph.GRID(i), adjacy_len, num_indices, adjacy + ind_sum);
    ind_sum += num_indices;
    adjacy_len -= num_indices;
    xadj[i+1] = ind_sum;
  }

  // build vtxdist array
  //int* vtxdst = new int[comm.NumProc()+1];
  //int myMaxGid = graph_map.MaxMyGID()+ 1;
  //comm.GatherAll (&myMaxGid, vtxdst+1, 1 );
  //vtxdst[0] = 0;

  return Repartition(graph_map, adjacy, xadj, num_myglobals, my_globals, debug);
}


int DistMesh::Repartition(const Epetra_BlockMap& graph_map, int* adjacy, int* xadj, int& num_myglobals, int*& my_globals, bool debug)
{

  //PAT_region_begin (4, "Repartition");
  int zero = 0;
  int ncon = 1;
  int nparts = comm.NumProc();
  int edgecut;
  int options[3] = {1, 0, 7};
  if (debug) {
      //options[0] = 1;
      options[1] = 7;
  }
  int* part = new int[graph_map.NumMyElements()];
  float* tpwgts = new float[nparts];
  for (int i=0; i<nparts; ++i)
      tpwgts[i] = 1.0/nparts;
  float ubvec = 1.05;

  int* vtxdst = new int[comm.NumProc()+1];
  int NumMyElements = graph_map.NumMyElements();
  int myMaxGid;
  int ret = comm.ScanSum(&NumMyElements, &myMaxGid, 1);
  ret = comm.GatherAll(&myMaxGid, vtxdst+1, 1 );
  vtxdst[0] = 0;

//   for(int i=0; i<graph_map.NumMyElements(); ++i) {
//       for (int j=xadj[i]; j<xadj[i+1]; ++j) { 
// 	  int lid = graph_map.LID(adjacy[j]);
// 	  int pid = comm.MyPID();
// 	  if (lid < 0) {
// 	      graph_map.RemoteIDList(1, &adjacy[j], &pid, &lid);
// 	  }
// 	  adjacy[j] = vtxdst[pid]+lid;
//       }
//       //sort(adjacy+xadj[i], adjacy+xadj[i+1]);
//   }



//    for (int j=0; j<xadj[graph_map.NumMyElements()]; ++j) {
//         int lid = graph_map.LID(adjacy[j]);
//         int pid = comm.MyPID();
//         if (lid < 0) {
//             graph_map.RemoteIDList(1, &adjacy[j], &pid, &lid);
//         }
//         cout <<  "pid: " << pid  << ", lid: " << lid << endl;
//         adjacy[j] = vtxdst[pid]+lid;
//    }

  std::map<int, std::pair<int,int> >* pid_map = new std::map<int, std::pair<int,int> >;
   //initialize gid_map
  std::pair<int, int> a_pair(-1,-1);
   for (int i=0; i<xadj[graph_map.NumMyElements()];++i)
       pid_map->operator[](adjacy[i]) = a_pair;

   int len = pid_map->size();
   int* pids = new int[len];
   int* gids = new int[len];
   int* lids = new int[len];
   int* int_it = gids;
   for (std::map<int, std::pair<int,int> >::iterator it = pid_map->begin(); it != pid_map->end(); ++it)
       *int_it++ = it->first;

   ret = graph_map.RemoteIDList(len, gids, pids, lids);

   int* pids_it = pids;
   int* lids_it = lids;
   for (std::map<int, std::pair<int,int> >::iterator it = pid_map->begin(); it != pid_map->end(); ++it) {
       it->second.first = *pids_it++;
       it->second.second = *lids_it++;
   }
   
   delete[] pids;
   delete[] gids;
   delete[] lids;

   for (int j=0; j<xadj[graph_map.NumMyElements()]; ++j) {
       adjacy[j] = vtxdst[pid_map->operator[](adjacy[j]).first]+pid_map->operator[](adjacy[j]).second;
   }

   delete pid_map;

  MPI_Comm mpi_comm = comm.Comm();

  CALL("ParMETIS_V3_PartKway");
  ParMETIS_V3_PartKway(vtxdst, xadj, adjacy, NULL, NULL,
	    &zero, &zero, &ncon, &nparts, tpwgts, &ubvec,
	    options, &edgecut, part, &mpi_comm );
  RET("ParMETIS_V3_PartKway", 0);
  //free graph memory
  delete[] xadj;
  delete[] adjacy;
  delete[] vtxdst;
  delete[] tpwgts;

  //PAT_region_begin (6, "Redistribution in Repartition");
  //Redistribution using Epetra_Distributor
  std::vector<int> export_PIDs;
  std::vector<int> export_objs;
  int num_export_IDs = 0;
  int num_remote_IDs;
  int count = 0;
  for (int i=0; i<graph_map.NumMyElements(); ++i) {
    if (part[i] != comm.MyPID()) {
      export_PIDs.push_back(part[i]);
      export_objs.push_back(graph_map.GID(i));
    } else {
      ++count;
    }
  }

  num_export_IDs = export_PIDs.size();

  Epetra_Distributor* part_dist = comm.CreateDistributor();
  ret = part_dist->CreateFromSends(num_export_IDs,&(*(export_PIDs.begin())), true, num_remote_IDs);

  int len_import_objs = num_remote_IDs*sizeof(int);
  int* import_objs = new int[len_import_objs];
  
  ret = part_dist->Do( reinterpret_cast<char*>(&(*export_objs.begin())), static_cast<int>(sizeof(int)), len_import_objs, reinterpret_cast<char*&>(import_objs) );
  //PAT_region_end (6);

  num_myglobals = count+num_remote_IDs;
  my_globals = new int[num_myglobals];
  int* it = my_globals;
  //add kept nodes to the new map
  for (int i=0; i<graph_map.NumMyElements(); ++i) {
    if (part[i] == comm.MyPID()) {
      *it++ = graph_map.GID(i);
    }
  }

  //add imported nodes to the new map
  int num_import_objs = 0;
  for (int i=0; i<num_remote_IDs; ++i) {
    *it++ = import_objs[num_import_objs];
    ++num_import_objs;
  }

  //sort array in case we have to reread the mesh file (there the indices are sorted)
  //sort(my_globals, it);
  
  // clean up datastructures
  delete[] part;
  delete[] import_objs;
  delete part_dist;
  //PAT_region_end(4);
  return num_myglobals;
}


int DistMesh::Repartition(const Epetra_Vector& coord_vec, int& num_myglobals, int*& my_globals, bool debug)
{
  //PAT_region_begin (4, "Repartition");
  // build vtxdist array
  const Epetra_BlockMap& coord_map = coord_vec.Map();
  int* vtxdst = new int[comm.NumProc()+1];
  int myMaxGid = coord_map.MaxMyGID()+ 1;
  int ret = comm.GatherAll (&myMaxGid, vtxdst+1, 1 );
  vtxdst[0] = 0;
  int ndims = Dimension;
  float* xyz = new float[coord_vec.MyLength()];
  for (int i=0; i< coord_vec.MyLength(); ++i)
    xyz[i] = float(coord_vec.operator[](i));
  int* part = new int[coord_map.NumMyElements()];
  MPI_Comm mpi_comm = comm.Comm();
  CALL("ParMETIS_V3_PartGeom");
  ParMETIS_V3_PartGeom (vtxdst, &ndims, xyz, part, &mpi_comm);
  RET("ParMETIS_V3_PartGeom", 0);
  delete[] vtxdst;
  delete[] xyz;
  
  //Redistribution using Epetra_Distributor
  std::vector<int> export_PIDs;
  std::vector<int> export_objs;
  int num_export_IDs = 0;
  int num_remote_IDs;
  int count = 0;
  for (int i=0; i<coord_map.NumMyElements(); ++i) {
    if (part[i] != comm.MyPID()) {
      export_PIDs.push_back(part[i]);
      export_objs.push_back(coord_map.GID(i));
    } else {
      ++count;
    }
  }

  num_export_IDs = export_PIDs.size();

  Epetra_Distributor* part_dist = comm.CreateDistributor();
  ret = part_dist->CreateFromSends(num_export_IDs,&(*(export_PIDs.begin())), true, num_remote_IDs);
  int len_import_objs = num_remote_IDs*sizeof(int);
  int* import_objs = new int[len_import_objs];
  ret = part_dist->Do( reinterpret_cast<char*>(&(*export_objs.begin())), static_cast<int>(sizeof(int)), len_import_objs, reinterpret_cast<char*&>(import_objs) );
  //PAT_region_end (6);
  
  num_myglobals = count+num_remote_IDs;
  my_globals = new int[num_myglobals];
  int* it = my_globals;
  //add kept nodes to the new map
  for (int i=0; i<coord_map.NumMyElements(); ++i) {
    if (part[i] == comm.MyPID()) {
      *it++ = coord_map.GID(i);
    }
  }

  //add imported nodes to the new map
  int num_import_objs = 0;
  for (int i=0; i<num_remote_IDs; ++i) {
    *it++ = import_objs[num_import_objs];
    ++num_import_objs;
  }

  //sort array in case we have to reread the mesh file (there the indices are sorted)
  //sort(my_globals, it);
  
  // clean up datastructures
  delete[] part;
  delete[] import_objs;
  delete part_dist;
  //PAT_region_end(4);
  return num_myglobals;
}


const Epetra_BlockMap* DistMesh::NodeMap()
{
  return &coordinates->Map();
}


const Epetra_BlockMap* DistMesh::ElementMap()
{
  return &element_nodes->Map();
}


Epetra_IntVector* DistMesh::ElementNodes()
{
  return element_nodes;
}


Epetra_IntVector* DistMesh::MaterialIDs()
{
  return mat_ids;
}


Epetra_Vector* DistMesh::Coordinates()
{
  return coordinates;
}

int DistMesh::ReferenceElement()
{
  int* values = element_nodes->Values();
  int len = element_nodes->MyLength();
  int mine = true;
  for (int i=0; i<len; i+=NumNodesPerElement) {
    int* elem_values = values+i;
    for (int j=0; j<NumNodesPerElement; ++j) {
      if (!coordinates->Map().MyGID(elem_values[j]))
	mine = false;
    }
    if (mine) 
      return element_nodes->Map().PointToElementList()[i];
    mine = true;
  }
  return -1;
}


void DistMesh::Print(MeshWriter* mw) {
    mw->PrintCoordinates(*coordinates);
    mw->PrintElements(*element_nodes);
}

void DistMesh::Print(ProblemWriter* pw) {
    pw->PrintCoordinates(*coordinates);
    pw->PrintElements(*element_nodes);
    if (MaterialIDMap.size() > 1) {
      std::vector<int> material_ids(MaterialIDMap.size());
      for (std::map<int, int>::const_iterator it = MaterialIDMap.begin(); it !=  MaterialIDMap.end(); ++it)
	material_ids[(*it).second] = (*it).first;
      
      for (int i=0; i<mat_ids->MyLength(); ++i)
	mat_ids->operator[](i) = material_ids[mat_ids->operator[](i)];
      pw->PrintMaterialIDs(*mat_ids);
    }
}

void DistMesh::Scan(ProblemReader* pr) {
    Epetra_BlockMap coordmap(NumNodes, Dimension, 0, comm);
    coordinates = new Epetra_Vector(coordmap);
    pr->ScanCoordinates(*coordinates);
    
    Epetra_BlockMap elemmap(NumElements, NumNodesPerElement, 0, comm);
    element_nodes = new Epetra_IntVector(elemmap);
    pr->ScanElements(*element_nodes);

    Epetra_Map matidmap(NumElements, 0, comm);
    mat_ids = new Epetra_IntVector(matidmap);
    if (MaterialIDMap.size() > 1) {
      pr->ScanMaterialIDs(*mat_ids);
      for (int i=0; i<mat_ids->MyLength(); ++i)
	mat_ids->operator[](i) = MaterialIDMap[mat_ids->operator[](i)];
    }
}

  
