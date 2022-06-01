/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * frename.c
 * 
 * This file contains some renaming routines to deal with different Fortran compilers
 *
 * Started 9/15/97
 * George
 *
 * $Id: frename.c,v 1.1 1998/11/27 17:59:14 karypis Exp $
 *
 */

#include <metis.h>


void METIS_PARTGRAPHRECURSIVE(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphrecursive(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{ 
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part); 
}
void metis_partgraphrecursive_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphrecursive__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_PartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}


void METIS_WPARTGRAPHRECURSIVE(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphrecursive(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphrecursive_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphrecursive__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_WPartGraphRecursive(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}



void METIS_PARTGRAPHKWAY(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphkway(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphkway_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_partgraphkway__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_PartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}



void METIS_WPARTGRAPHKWAY(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphkway(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphkway_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}
void metis_wpartgraphkway__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_WPartGraphKway(nvtxs, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, tpwgts, options, edgecut, part);
}



void METIS_EDGEND(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_edgend(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_edgend_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_edgend__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_EdgeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}



void METIS_NODEND(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_nodend(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_nodend_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}
void metis_nodend__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeND(nvtxs, xadj, adjncy, numflag, options, perm, iperm);
}



void METIS_NODEWND(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
void metis_nodewnd(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
void metis_nodewnd_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}
void metis_nodewnd__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, my_int *numflag, my_int *options, idxtype *perm, idxtype *iperm)
{
  METIS_NodeWND(nvtxs, xadj, adjncy, vwgt, numflag, options, perm, iperm);
}



void METIS_PARTMESHNODAL(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, my_int *nparts, my_int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshnodal(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, my_int *nparts, my_int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshnodal_(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, my_int *nparts, my_int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshnodal__(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, my_int *nparts, my_int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshNodal(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}


void METIS_PARTMESHDUAL(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, my_int *nparts, my_int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshdual(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, my_int *nparts, my_int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshdual_(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, my_int *nparts, my_int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}
void metis_partmeshdual__(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, my_int *nparts, my_int *edgecut, idxtype *epart, idxtype *npart)
{
  METIS_PartMeshDual(ne, nn, elmnts, etype, numflag, nparts, edgecut, epart, npart);
}


void METIS_MESHTONODAL(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtonodal(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtonodal_(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtonodal__(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToNodal(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}


void METIS_MESHTODUAL(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtodual(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtodual_(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}
void metis_meshtodual__(my_int *ne, my_int *nn, idxtype *elmnts, my_int *etype, my_int *numflag, idxtype *dxadj, idxtype *dadjncy)
{
  METIS_MeshToDual(ne, nn, elmnts, etype, numflag, dxadj, dadjncy);
}


void METIS_ESTIMATEMEMORY(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *optype, my_int *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}
void metis_estimatememory(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *optype, my_int *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}
void metis_estimatememory_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *optype, my_int *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}
void metis_estimatememory__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, my_int *numflag, my_int *optype, my_int *nbytes)
{
  METIS_EstimateMemory(nvtxs, xadj, adjncy, numflag, optype, nbytes);
}



void METIS_MCPARTGRAPHRECURSIVE(my_int *nvtxs, my_int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_mcpartgraphrecursive(my_int *nvtxs, my_int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_mcpartgraphrecursive_(my_int *nvtxs, my_int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}
void metis_mcpartgraphrecursive__(my_int *nvtxs, my_int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_mCPartGraphRecursive(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, options, edgecut, part);
}


void METIS_MCPARTGRAPHKWAY(my_int *nvtxs, my_int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *rubvec, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}
void metis_mcpartgraphkway(my_int *nvtxs, my_int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *rubvec, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}
void metis_mcpartgraphkway_(my_int *nvtxs, my_int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *rubvec, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}
void metis_mcpartgraphkway__(my_int *nvtxs, my_int *ncon, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt, my_int *wgtflag, my_int *numflag, my_int *nparts, float *rubvec, my_int *options, my_int *edgecut, idxtype *part)
{
  METIS_mCPartGraphKway(nvtxs, ncon, xadj, adjncy, vwgt, adjwgt, wgtflag, numflag, nparts, rubvec, options, edgecut, part);
}


void METIS_PARTGRAPHVKWAY(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}
void metis_partgraphvkway(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}
void metis_partgraphvkway_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}
void metis_partgraphvkway__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, my_int *wgtflag, my_int *numflag, my_int *nparts, my_int *options, my_int *volume, idxtype *part)
{
  METIS_PartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, options, volume, part);
}

void METIS_WPARTGRAPHVKWAY(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}
void metis_wpartgraphvkway(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}
void metis_wpartgraphvkway_(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}
void metis_wpartgraphvkway__(my_int *nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *vsize, my_int *wgtflag, my_int *numflag, my_int *nparts, float *tpwgts, my_int *options, my_int *volume, idxtype *part)
{
  METIS_WPartGraphVKway(nvtxs, xadj, adjncy, vwgt, vsize, wgtflag, numflag, nparts, tpwgts, options, volume, part);
}



