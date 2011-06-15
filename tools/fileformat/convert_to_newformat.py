#!/usr/bin/python2

import h5py
import sys

if  len(sys.argv) != 2:
  print "usage: "+sys.argv[0]+" filename"
  sys.exit(0)
try:
  f = h5py.File(sys.argv[1],'r+')
except:
  print "usgage: "+sys.argv[0]+" filename"
  sys.exit(2)

try:
  f.id.move('/Boundary conditions' , '/Boundary_conditions')
except:
  print "/Boundary conditions does not exist"

changing_sets = {
  '/Boundary_conditions/Fixed nodes'             : '/Boundary_conditions/fixed_nodes',
  '/Boundary_conditions/Fixed nodes size'        : '/Boundary_conditions/fixed_nodes_size',
  '/Boundary_conditions/Loaded nodes'            : '/Boundary_conditions/loaded_nodes',
  '/Boundary_conditions/Loaded nodes size'       : '/Boundary_conditions/loaded_nodes_size',
  '/Boundary_conditions/Restrained nodes'        : '/Boundary_conditions/restrained_nodes',
  '/Boundary_conditions/Restrained nodes size'   : '/Boundary_conditions/restrained_nodes_size',
  '/Parameters/#dimensions'                      : '/Parameters/nr_dimensions',
  '/Parameters/#dofs per node'                   : '/Parameters/nr_dofs_per_node',
  '/Parameters/#elements'                        : '/Parameters/nr_elements',
  '/Parameters/#integration points'              : '/Parameters/nr_integration_points',
  '/Parameters/#material properties'             : '/Parameters/nr_material_properties',
  '/Parameters/#material types'                  : '/Parameters/nr_material_types',
  '/Parameters/#nodes'                           : '/Parameters/nr_nodes',
  '/Parameters/#nodes per element'               : '/Parameters/nr_nodes_per_element',
  '/Parameters/element type'                     : '/Parameters/element_type',
  '/Parameters/iteration limit'                  : '/Parameters/iteration_limit',
  '/Parameters/size of stress-strain matrix'     : '/Parameters/size_of_stress_strain_matrix',
  '/Solution/#Gauss points'                      : '/Solution/nr_gauss_points',
  '/Solution/Element stress'                     : '/Solution/element_stress',
  '/Solution/Element strain'                     : '/Solution/element_strain',
  '/Solution/Nodal displacements'                : '/Solution/nodal_displacements',
  '/Solution/Nodal forces'                       : '/Solution/nodal_forces',
  }
for key in changing_sets.keys():
  try:
    f.id.move(key, changing_sets[key])
  except:
    print key+" does not exist"
f.close()
