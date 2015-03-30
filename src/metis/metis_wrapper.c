#include "metislib.h"

#define FRENAME(name, dargs, cargs, name1, name2, name3, name4)   \
  int name1 dargs { return name cargs; }                          \
  int name2 dargs { return name cargs; }                          \
  int name3 dargs { return name cargs; }                          \
  int name4 dargs { return name cargs; }

int METIS_PartGraphRecursive_wrapper(int *nvtxs, int *xadj, int *adjncy, 
	int *vwgt, int *adjwgt, int *nparts, int *objval, int *part)
{
   int ncon;
   ncon = 1;
   METIS_PartGraphRecursive(nvtxs, &ncon, xadj, adjncy, vwgt, 
          NULL, adjwgt, nparts, NULL, NULL, NULL, objval, part);
   return 0;
}

FRENAME(
    METIS_PartGraphRecursive_wrapper, 
    (int *nvtxs, int *xadj, int *adjncy, int *vwgt,
     int *adjwgt, int *nparts, int *edgecut, int *part),
    (nvtxs, xadj, adjncy, vwgt, adjwgt, nparts,
     edgecut, part),
    METIS_PARTGRAPHRECURSIVE_WRAPPER, 
    metis_partgraphrecursive_wrapper, 
    metis_partgraphrecursive_wrapper_, 
    metis_partgraphrecursive_wrapper__
) 
