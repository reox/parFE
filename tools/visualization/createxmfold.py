#!/usr/bin/python

import h5py
import sys

if  len(sys.argv) != 2:
  print "usage: "+sys.argv[0]+" filename"
  sys.exit(0)

filename = sys.argv[1]
try:
  f = h5py.File(filename,'r+')
except:
  print "usgage: "+sys.argv[0]+" filename"
  sys.exit(2)

#get numbers
nr_nodes = f['Parameters']['#nodes'].value
nr_elements = f['Parameters']['#elements'].value

print "nodes: ", nr_nodes, " elements: ", nr_elements

# xml
sStart = "<?xml version=\"1.0\" ?>\n" \
         "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" \
         "<Xdmf Version=\"2.0\">\n" \
         " <Domain>\n" \
         "   <Grid Name=\"mesh1\" GridType=\"Uniform\">\n" 

sMesh =  "     <Topology TopologyType=\"Hexahedron\" NumberOfElements=\""+repr(nr_elements)+"\" BaseOffset=\"1\" >\n" \
         "       <DataItem Dimensions=\""+repr(nr_elements)+" 8\" NumberType=\"Int\" Format=\"HDF\">\n" \
         "         "+filename+":/Mesh/Elements\n" \
         "       </DataItem>\n" \
         "     </Topology>\n" \
         "     <Geometry GeometryType=\"XYZ\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_nodes)+" 3\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n" \
         "         "+filename+":/Mesh/Coordinates\n" \
         "       </DataItem>\n" \
         "     </Geometry>\n"

sDisp =  "     <Attribute Name=\"Displacement\" AttributeType=\"Vector\" Center=\"Node\">\n" \
         "       <DataItem Dimensions=\""+repr(nr_nodes)+" 3\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "         "+filename+":/Solution/Nodal displacements\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

sSED =   "     <Attribute Name=\"SED\" AttributeType=\"Scalar\" Center=\"Cell\">\n" \
         "       <DataItem ItemType=\"HyperSlab\" Dimensions=\""+repr(nr_elements)+" 1\" Type=\"HyperSlab\">\n" \
         "         <DataItem Dimensions=\"3 2\" Format=\"XML\">\n" \
         "           0 6      <!-- START  : from the small dimension the 6th element (7th column) --> \n" \
         "           1 1      <!-- STRIDE : from the 2nd dim all elements -->\n" \
         "           "+repr(nr_elements)+" 1  <!-- COUNT  : -->\n" \
         "         </DataItem>\n" \
         "         <DataItem Dimensions=\""+repr(nr_elements)+" 8\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "           "+filename+":/Solution/Element strain\n" \
         "         </DataItem>\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

sVM =    "     <Attribute Name=\"S_vonMises\" AttributeType=\"Scalar\" Center=\"Cell\">\n" \
         "       <DataItem ItemType=\"HyperSlab\" Dimensions=\""+repr(nr_elements)+" 1\" Type=\"HyperSlab\">\n" \
         "         <DataItem Dimensions=\"3 2\" Format=\"XML\">\n" \
         "           0 6\n" \
         "           1 1\n" \
         "           "+repr(nr_elements)+" 1\n" \
         "         </DataItem>\n" \
         "         <DataItem Dimensions=\""+repr(nr_elements)+" 7\" Format=\"HDF\" NumberType=\"Float\" Precision=\"8\" >\n" \
         "           "+filename+":/Solution/Element stress\n" \
         "         </DataItem>\n" \
         "       </DataItem>\n" \
         "     </Attribute>\n"

sEnd =   "   </Grid>\n" \
         " </Domain>\n" \
         "</Xdmf>\n" 
#write file
outfilename = filename.replace("mesh.h5","xmf")
outfile = open(outfilename, 'w')
outfile.write("%s%s%s%s%s%s" % (sStart, sMesh, sDisp, sSED, sVM, sEnd))
outfile.close()
