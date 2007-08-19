#include "parfe_ConfigDefs.h"
#include "Tools.h"
#if defined(HAVE_MALLINFO)
#include <malloc.h>
#elif defined(HAVE_HEAP_INFO)
#include <catamount/catmalloc.h>
#endif
#include "DistMesh.h"

unsigned meminfo()
{
#if defined(HAVE_MALLINFO)
    long int fragments, total_free, largest_free, total_used;
    /* int percent; */

    /* use system call to get memory used information */

    struct mallinfo M = mallinfo();
    fragments = M.ordblks + M.smblks + M.hblks;
    total_free = M.fsmblks + M.fordblks;
    total_used = M.hblkhd + M.usmblks + M.uordblks;
    /*  total_free = ml_total_mem - total_used;  */
    largest_free = -1024;

    /* convert to Mbytes */

    return( (unsigned)(total_used/(1024*1024)) );

#elif defined(HAVE_HEAP_INFO)
    size_t fragments;
    unsigned long total_free;
    unsigned long largest_free;
    unsigned long total_used;

    heap_info(&fragments, &total_free, &largest_free, &total_used);
   
    return( (unsigned)(total_used/(1024*1024)) );

#else
    return(0);
#endif
}

// ============================================================================ 
void set_null_space(DistMesh* mesh, vector<double>& null_space)
{
  int num_rows =  mesh->NodeMap()->NumMyPoints();
  null_space.resize(6 * num_rows);

  double* row_it = &null_space[0];

  //clear vector
  for (double* it = &null_space[0]; it < &null_space[0] + 6 * num_rows; ++it)
    *it = 0;

  //first coordinate is constant
  for (double* it = row_it; it < row_it + num_rows; it += 3)
    *it = 1;

  //second coordinate is constant
  row_it += num_rows;
  for (double* it = row_it+1;it < row_it+num_rows; it+=3)
    *it = 1;

  //third coordinate is constant
  row_it += num_rows;
  for (double* it = row_it+2;it < row_it+num_rows; it+=3)
    *it = 1;

  // u = y and v = -x
  double* values = mesh->Coordinates()->Values();
  row_it += num_rows;
  for (double* it = row_it;it < row_it+num_rows; it+=3) {
    *it = *(values+1);
    *(it+1) = -*values;
    values+=3;
  }

  // v = z and w = -y
  values = mesh->Coordinates()->Values();
  row_it += num_rows;
  for (double* it = row_it;it < row_it+num_rows; it+=3) {
    *(it+1) = *(values+2);
    *(it+2) = -*(values+1);
    values+=3;
  }

  // u = z and w = -x
  values = mesh->Coordinates()->Values();
  row_it += num_rows;
  for (double* it = row_it;it < row_it+num_rows; it+=3) {
    *it = *(values+2);
    *(it+2) = -*values;
    values+=3;
  }
}
