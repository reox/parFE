#ifndef TOOLS_H
#define TOOLS_H
#endif

#include <vector>
class DistMesh;

//Helper function to compute memory usage
unsigned meminfo();

void set_null_space(DistMesh* mesh, vector<double>& null_space);
