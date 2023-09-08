/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h,v 1.1 1998/11/27 17:59:28 karypis Exp $
 *
 */

/* balance.c */
void Balance2Way(CtrlType *, GraphType *, my_int *, float);
void Bnd2WayBalance(CtrlType *, GraphType *, my_int *);
void General2WayBalance(CtrlType *, GraphType *, my_int *);

/* bucketsort.c */
void BucketSortKeysInc(my_int, my_int, idxtype *, idxtype *, idxtype *);

/* ccgraph.c */
void CreateCoarseGraph(CtrlType *, GraphType *, my_int, idxtype *, idxtype *);
void CreateCoarseGraphNoMask(CtrlType *, GraphType *, my_int, idxtype *, idxtype *);
void CreateCoarseGraph_NVW(CtrlType *, GraphType *, my_int, idxtype *, idxtype *);
GraphType *SetUpCoarseGraph(GraphType *, my_int, my_int);
void ReAdjustMemory(GraphType *, GraphType *, my_int);

/* coarsen.c */
GraphType *Coarsen2Way(CtrlType *, GraphType *);

/* compress.c */
void CompressGraph(CtrlType *, GraphType *, my_int, idxtype *, idxtype *, idxtype *, idxtype *);
void PruneGraph(CtrlType *, GraphType *, my_int, idxtype *, idxtype *, idxtype *, float);

/* debug.c */
my_int ComputeCut(GraphType *, idxtype *);
my_int CheckBnd(GraphType *);
my_int CheckBnd2(GraphType *);
my_int CheckNodeBnd(GraphType *, my_int);
my_int CheckRInfo(RInfoType *);
my_int CheckNodePartitionParams(GraphType *);
my_int IsSeparable(GraphType *);

/* estmem.c */
void METIS_EstimateMemory(my_int *, idxtype *, idxtype *, my_int *, my_int *, my_int *);
void EstimateCFraction(my_int, idxtype *, idxtype *, float *, float *);
my_int ComputeCoarseGraphSize(my_int, idxtype *, idxtype *, my_int, idxtype *, idxtype *, idxtype *);

/* fm.c */
void FM_2WayEdgeRefine(CtrlType *, GraphType *, my_int *, my_int);

/* fortran.c */
void Change2CNumbering(my_int, idxtype *, idxtype *);
void Change2FNumbering(my_int, idxtype *, idxtype *, idxtype *);
void Change2FNumbering2(my_int, idxtype *, idxtype *);
void Change2FNumberingOrder(my_int, idxtype *, idxtype *, idxtype *, idxtype *);
void ChangeMesh2CNumbering(my_int, idxtype *);
void ChangeMesh2FNumbering(my_int, idxtype *, my_int, idxtype *, idxtype *);
void ChangeMesh2FNumbering2(my_int, idxtype *, my_int, my_int, idxtype *, idxtype *);

/* frename.c */
void METIS_PARTGRAPHRECURSIVE(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void metis_partgraphrecursive(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void metis_partgraphrecursive_(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void metis_partgraphrecursive__(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void METIS_WPARTGRAPHRECURSIVE(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
void metis_wpartgraphrecursive(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
void metis_wpartgraphrecursive_(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
void metis_wpartgraphrecursive__(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
void METIS_PARTGRAPHKWAY(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void metis_partgraphkway(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void metis_partgraphkway_(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void metis_partgraphkway__(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void METIS_WPARTGRAPHKWAY(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
void metis_wpartgraphkway(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
void metis_wpartgraphkway_(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
void metis_wpartgraphkway__(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
void METIS_EDGEND(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void metis_edgend(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void metis_edgend_(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void metis_edgend__(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void METIS_NODEND(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void metis_nodend(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void metis_nodend_(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void metis_nodend__(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void METIS_NODEWND(my_int *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void metis_nodewnd(my_int *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void metis_nodewnd_(my_int *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void metis_nodewnd__(my_int *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void METIS_PARTMESHNODAL(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);
void metis_partmeshnodal(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);
void metis_partmeshnodal_(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);
void metis_partmeshnodal__(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);
void METIS_PARTMESHDUAL(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);
void metis_partmeshdual(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);
void metis_partmeshdual_(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);
void metis_partmeshdual__(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);
void METIS_MESHTONODAL(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void metis_meshtonodal(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void metis_meshtonodal_(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void metis_meshtonodal__(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void METIS_MESHTODUAL(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void metis_meshtodual(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void metis_meshtodual_(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void metis_meshtodual__(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void METIS_ESTIMATEMEMORY(my_int *, idxtype *, idxtype *, my_int *, my_int *, my_int *);
void metis_estimatememory(my_int *, idxtype *, idxtype *, my_int *, my_int *, my_int *);
void metis_estimatememory_(my_int *, idxtype *, idxtype *, my_int *, my_int *, my_int *);
void metis_estimatememory__(my_int *, idxtype *, idxtype *, my_int *, my_int *, my_int *);
void METIS_MCPARTGRAPHRECURSIVE(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void metis_mcpartgraphrecursive(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void metis_mcpartgraphrecursive_(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void metis_mcpartgraphrecursive__(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void METIS_MCPARTGRAPHKWAY(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
void metis_mcpartgraphkway(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
void metis_mcpartgraphkway_(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
void metis_mcpartgraphkway__(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
void METIS_PARTGRAPHVKWAY(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void metis_partgraphvkway(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void metis_partgraphvkway_(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void metis_partgraphvkway__(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void METIS_WPARTGRAPHVKWAY(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
void metis_wpartgraphvkway(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
void metis_wpartgraphvkway_(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
void metis_wpartgraphvkway__(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);

/* graph.c */
void SetUpGraph(GraphType *, my_int, my_int, my_int, idxtype *, idxtype *, idxtype *, idxtype *, my_int);
void SetUpGraphKway(GraphType *, my_int, idxtype *, idxtype *);
void SetUpGraph2(GraphType *, my_int, my_int, idxtype *, idxtype *, float *, idxtype *);
void VolSetUpGraph(GraphType *, my_int, my_int, my_int, idxtype *, idxtype *, idxtype *, idxtype *, my_int);
void RandomizeGraph(GraphType *);
my_int IsConnectedSubdomain(CtrlType *, GraphType *, my_int, my_int);
my_int IsConnected(CtrlType *, GraphType *, my_int);
my_int IsConnected2(GraphType *, my_int);
my_int FindComponents(CtrlType *, GraphType *, idxtype *, idxtype *);

/* initpart.c */
void Init2WayPartition(CtrlType *, GraphType *, my_int *, float);
void InitSeparator(CtrlType *, GraphType *, float);
void GrowBisection(CtrlType *, GraphType *, my_int *, float);
void GrowBisectionNode(CtrlType *, GraphType *, float);
void RandomBisection(CtrlType *, GraphType *, my_int *, float);

/* kmetis.c */
void METIS_PartGraphKway(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void METIS_WPartGraphKway(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
my_int MlevelKWayPartitioning(CtrlType *, GraphType *, my_int, idxtype *, float *, float);

/* kvmetis.c */
void METIS_PartGraphVKway(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void METIS_WPartGraphVKway(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
my_int MlevelVolKWayPartitioning(CtrlType *, GraphType *, my_int, idxtype *, float *, float);

/* kwayfm.c */
void Random_KWayEdgeRefine(CtrlType *, GraphType *, my_int, float *, float, my_int, my_int);
void Greedy_KWayEdgeRefine(CtrlType *, GraphType *, my_int, float *, float, my_int);
void Greedy_KWayEdgeBalance(CtrlType *, GraphType *, my_int, float *, float, my_int);

/* kwayrefine.c */
void RefineKWay(CtrlType *, GraphType *, GraphType *, my_int, float *, float);
void AllocateKWayPartitionMemory(CtrlType *, GraphType *, my_int);
void ComputeKWayPartitionParams(CtrlType *, GraphType *, my_int);
void ProjectKWayPartition(CtrlType *, GraphType *, my_int);
my_int IsBalanced(idxtype *, my_int, float *, float);
void ComputeKWayBoundary(CtrlType *, GraphType *, my_int);
void ComputeKWayBalanceBoundary(CtrlType *, GraphType *, my_int);

/* kwayvolfm.c */
void Random_KWayVolRefine(CtrlType *, GraphType *, my_int, float *, float, my_int, my_int);
void Random_KWayVolRefineMConn(CtrlType *, GraphType *, my_int, float *, float, my_int, my_int);
void Greedy_KWayVolBalance(CtrlType *, GraphType *, my_int, float *, float, my_int);
void Greedy_KWayVolBalanceMConn(CtrlType *, GraphType *, my_int, float *, float, my_int);
void KWayVolUpdate(CtrlType *, GraphType *, my_int, my_int, my_int, idxtype *, idxtype *, idxtype *);
void ComputeKWayVolume(GraphType *, my_int, idxtype *, idxtype *, idxtype *);
my_int ComputeVolume(GraphType *, idxtype *);
void CheckVolKWayPartitionParams(CtrlType *, GraphType *, my_int);
void ComputeVolSubDomainGraph(GraphType *, my_int, idxtype *, idxtype *);
void EliminateVolSubDomainEdges(CtrlType *, GraphType *, my_int, float *);
void EliminateVolComponents(CtrlType *, GraphType *, my_int, float *, float);

/* kwayvolrefine.c */
void RefineVolKWay(CtrlType *, GraphType *, GraphType *, my_int, float *, float);
void AllocateVolKWayPartitionMemory(CtrlType *, GraphType *, my_int);
void ComputeVolKWayPartitionParams(CtrlType *, GraphType *, my_int);
void ComputeKWayVolGains(CtrlType *, GraphType *, my_int);
void ProjectVolKWayPartition(CtrlType *, GraphType *, my_int);
void ComputeVolKWayBoundary(CtrlType *, GraphType *, my_int);
void ComputeVolKWayBalanceBoundary(CtrlType *, GraphType *, my_int);

/* match.c */
void Match_RM(CtrlType *, GraphType *);
void Match_RM_NVW(CtrlType *, GraphType *);
void Match_HEM(CtrlType *, GraphType *);
void Match_SHEM(CtrlType *, GraphType *);

/* mbalance.c */
void MocBalance2Way(CtrlType *, GraphType *, float *, float);
void MocGeneral2WayBalance(CtrlType *, GraphType *, float *, float);

/* mbalance2.c */
void MocBalance2Way2(CtrlType *, GraphType *, float *, float *);
void MocGeneral2WayBalance2(CtrlType *, GraphType *, float *, float *);
void SelectQueue3(my_int, float *, float *, my_int *, my_int *, PQueueType [MAXNCON][2], float *);

/* mcoarsen.c */
GraphType *MCCoarsen2Way(CtrlType *, GraphType *);

/* memory.c */
void AllocateWorkSpace(CtrlType *, GraphType *, my_int);
void FreeWorkSpace(CtrlType *, GraphType *);
my_int WspaceAvail(CtrlType *);
idxtype *idxwspacemalloc(CtrlType *, my_int);
void idxwspacefree(CtrlType *, my_int);
float *fwspacemalloc(CtrlType *, my_int);
void fwspacefree(CtrlType *, my_int);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeGraph(GraphType *);

/* mesh.c */
void METIS_MeshToDual(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void METIS_MeshToNodal(my_int *, my_int *, idxtype *, my_int *, my_int *, idxtype *, idxtype *);
void GENDUALMETIS(my_int, my_int, my_int, idxtype *, idxtype *, idxtype *adjncy);
void TRINODALMETIS(my_int, my_int, idxtype *, idxtype *, idxtype *adjncy);
void TETNODALMETIS(my_int, my_int, idxtype *, idxtype *, idxtype *adjncy);
void HEXNODALMETIS(my_int, my_int, idxtype *, idxtype *, idxtype *adjncy);
void QUADNODALMETIS(my_int, my_int, idxtype *, idxtype *, idxtype *adjncy);

/* meshpart.c */
void METIS_PartMeshNodal(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);
void METIS_PartMeshDual(my_int *, my_int *, idxtype *, my_int *, my_int *, my_int *, my_int *, idxtype *, idxtype *);

/* mfm.c */
void MocFM_2WayEdgeRefine(CtrlType *, GraphType *, float *, my_int);
void SelectQueue(my_int, float *, float *, my_int *, my_int *, PQueueType [MAXNCON][2]);
my_int BetterBalance(my_int, float *, float *, float *);
float Compute2WayHLoadImbalance(my_int, float *, float *);
void Compute2WayHLoadImbalanceVec(my_int, float *, float *, float *);

/* mfm2.c */
void MocFM_2WayEdgeRefine2(CtrlType *, GraphType *, float *, float *, my_int);
void SelectQueue2(my_int, float *, float *, my_int *, my_int *, PQueueType [MAXNCON][2], float *);
my_int IsBetter2wayBalance(my_int, float *, float *, float *);

/* mincover.o */
void MinCover(idxtype *, idxtype *, my_int, my_int, idxtype *, my_int *);
my_int MinCover_Augment(idxtype *, idxtype *, my_int, idxtype *, idxtype *, idxtype *, my_int);
void MinCover_Decompose(idxtype *, idxtype *, my_int, my_int, idxtype *, idxtype *, my_int *);
void MinCover_ColDFS(idxtype *, idxtype *, my_int, idxtype *, idxtype *, my_int);
void MinCover_RowDFS(idxtype *, idxtype *, my_int, idxtype *, idxtype *, my_int);

/* minitpart.c */
void MocInit2WayPartition(CtrlType *, GraphType *, float *, float);
void MocGrowBisection(CtrlType *, GraphType *, float *, float);
void MocRandomBisection(CtrlType *, GraphType *, float *, float);
void MocInit2WayBalance(CtrlType *, GraphType *, float *);
my_int SelectQueueoneWay(my_int, float *, float *, my_int, PQueueType [MAXNCON][2]);

/* minitpart2.c */
void MocInit2WayPartition2(CtrlType *, GraphType *, float *, float *);
void MocGrowBisection2(CtrlType *, GraphType *, float *, float *);
void MocGrowBisectionNew2(CtrlType *, GraphType *, float *, float *);
void MocInit2WayBalance2(CtrlType *, GraphType *, float *, float *);
my_int SelectQueueOneWay2(my_int, float *, PQueueType [MAXNCON][2], float *);

/* mkmetis.c */
void METIS_mCPartGraphKway(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
my_int MCMlevelKWayPartitioning(CtrlType *, GraphType *, my_int, idxtype *, float *);

/* mkwayfmh.c */
void MCRandom_KWayEdgeRefineHorizontal(CtrlType *, GraphType *, my_int, float *, my_int);
void MCGreedy_KWayEdgeBalanceHorizontal(CtrlType *, GraphType *, my_int, float *, my_int);
my_int AreAllHVwgtsBelow(my_int, float, float *, float, float *, float *);
my_int AreAllHVwgtsAbove(my_int, float, float *, float, float *, float *);
void ComputeHKWayLoadImbalance(my_int, my_int, float *, float *);
my_int MocIsHBalanced(my_int, my_int, float *, float *);
my_int IsHBalanceBetterFT(my_int, my_int, float *, float *, float *, float *);
my_int IsHBalanceBetterTT(my_int, my_int, float *, float *, float *, float *);

/* mkwayrefine.c */
void MocRefineKWayHorizontal(CtrlType *, GraphType *, GraphType *, my_int, float *);
void MocAllocateKWayPartitionMemory(CtrlType *, GraphType *, my_int);
void MocComputeKWayPartitionParams(CtrlType *, GraphType *, my_int);
void MocProjectKWayPartition(CtrlType *, GraphType *, my_int);
void MocComputeKWayBalanceBoundary(CtrlType *, GraphType *, my_int);

/* mmatch.c */
void MCMatch_RM(CtrlType *, GraphType *);
void MCMatch_HEM(CtrlType *, GraphType *);
void MCMatch_SHEM(CtrlType *, GraphType *);
void MCMatch_SHEBM(CtrlType *, GraphType *, my_int);
void MCMatch_SBHEM(CtrlType *, GraphType *, my_int);
float BetterVBalance(my_int, my_int, float *, float *, float *);
my_int AreAllVwgtsBelowFast(my_int, float *, float *, float);

/* mmd.c */
void genmmd(my_int, idxtype *, idxtype *, idxtype *, idxtype *, my_int , idxtype *, idxtype *, idxtype *, idxtype *, my_int, my_int *);
void mmdelm(my_int, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, my_int, my_int);
my_int  mmdint(my_int, idxtype *xadj, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *);
void mmdnum(my_int, idxtype *, idxtype *, idxtype *);
void mmdupd(my_int, my_int, idxtype *, idxtype *, my_int, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, my_int, my_int *tag);

/* mpmetis.c */
void METIS_mCPartGraphRecursive(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *);
void METIS_mCHPartGraphRecursive(my_int *, my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *);
void METIS_mCPartGraphRecursiveInternal(my_int *, my_int *, idxtype *, idxtype *, float *, idxtype *, my_int *, my_int *, my_int *, idxtype *);
void METIS_mCHPartGraphRecursiveInternal(my_int *, my_int *, idxtype *, idxtype *, float *, idxtype *, my_int *, float *, my_int *, my_int *, idxtype *);
my_int MCMlevelRecursiveBisection(CtrlType *, GraphType *, my_int, idxtype *, float, my_int);
my_int MCHMlevelRecursiveBisection(CtrlType *, GraphType *, my_int, idxtype *, float *, my_int);
void MCMlevelEdgeBisection(CtrlType *, GraphType *, float *, float);
void MCHMlevelEdgeBisection(CtrlType *, GraphType *, float *, float *);

/* mrefine.c */
void MocRefine2Way(CtrlType *, GraphType *, GraphType *, float *, float);
void MocAllocate2WayPartitionMemory(CtrlType *, GraphType *);
void MocCompute2WayPartitionParams(CtrlType *, GraphType *);
void MocProject2WayPartition(CtrlType *, GraphType *);

/* mrefine2.c */
void MocRefine2Way2(CtrlType *, GraphType *, GraphType *, float *, float *);

/* mutil.c */
my_int AreAllVwgtsBelow(my_int, float, float *, float, float *, float);
my_int AreAnyVwgtsBelow(my_int, float, float *, float, float *, float);
my_int AreAllVwgtsAbove(my_int, float, float *, float, float *, float);
float ComputeLoadImbalance(my_int, my_int, float *, float *);
my_int AreAllBelow(my_int, float *, float *);

/* myqsort.c */
void iidxsort(my_int, idxtype *);
void iintsort(my_int, my_int *);
void ikeysort(my_int, KeyValueType *);
void ikeyvalsort(my_int, KeyValueType *);

/* ometis.c */
void METIS_EdgeND(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void METIS_NodeND(my_int *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void METIS_NodeWND(my_int *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, idxtype *, idxtype *); 
void MlevelNestedDissection(CtrlType *, GraphType *, idxtype *, float, my_int);
void MlevelNestedDissectionCC(CtrlType *, GraphType *, idxtype *, float, my_int);
void MlevelNodeBisectionMultiple(CtrlType *, GraphType *, my_int *, float);
void MlevelNodeBisection(CtrlType *, GraphType *, my_int *, float);
void SplitGraphOrder(CtrlType *, GraphType *, GraphType *, GraphType *);
void MMDOrder(CtrlType *, GraphType *, idxtype *, my_int);
my_int SplitGraphOrderCC(CtrlType *, GraphType *, GraphType *, my_int, idxtype *, idxtype *);

/* parmetis.c */
void METIS_PartGraphKway2(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void METIS_WPartGraphKway2(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
void METIS_NodeNDP(my_int, idxtype *, idxtype *, my_int, my_int *, idxtype *, idxtype *, idxtype *);
void MlevelNestedDissectionP(CtrlType *, GraphType *, idxtype *, my_int, my_int, my_int, idxtype *);
void METIS_NodeComputeSeparator(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, idxtype *); 
void METIS_EdgeComputeSeparator(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, idxtype *); 

/* pmetis.c */
void METIS_PartGraphRecursive(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, my_int *, my_int *, idxtype *); 
void METIS_WPartGraphRecursive(my_int *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, my_int *, my_int *, float *, my_int *, my_int *, idxtype *); 
my_int MlevelRecursiveBisection(CtrlType *, GraphType *, my_int, idxtype *, float *, float, my_int);
void MlevelEdgeBisection(CtrlType *, GraphType *, my_int *, float);
void SplitGraphPart(CtrlType *, GraphType *, GraphType *, GraphType *);
void SetUpSplitGraph(GraphType *, GraphType *, my_int, my_int);

/* pqueue.c */
void PQueueInit(CtrlType *ctrl, PQueueType *, my_int, my_int);
void PQueueReset(PQueueType *);
void PQueueFree(CtrlType *ctrl, PQueueType *);
my_int PQueueGetSize(PQueueType *);
my_int PQueueInsert(PQueueType *, my_int, my_int);
my_int PQueueDelete(PQueueType *, my_int, my_int);
my_int PQueueUpdate(PQueueType *, my_int, my_int, my_int);
void PQueueUpdateUp(PQueueType *, my_int, my_int, my_int);
my_int PQueueGetMax(PQueueType *);
my_int PQueueSeeMax(PQueueType *);
my_int PQueueGetKey(PQueueType *);
my_int CheckHeap(PQueueType *);

/* refine.c */
void Refine2Way(CtrlType *, GraphType *, GraphType *, my_int *, float ubfactor);
void Allocate2WayPartitionMemory(CtrlType *, GraphType *);
void Compute2WayPartitionParams(CtrlType *, GraphType *);
void Project2WayPartition(CtrlType *, GraphType *);

/* separator.c */
void ConstructSeparator(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator0(CtrlType *, GraphType *, float);
void ConstructMinCoverSeparator(CtrlType *, GraphType *, float);

/* sfm.c */
void FM_2WayNodeRefine(CtrlType *, GraphType *, float, my_int);
void FM_2WayNodeRefineEqWgt(CtrlType *, GraphType *, my_int);
void FM_2WayNodeRefine_OneSided(CtrlType *, GraphType *, float, my_int);
void FM_2WayNodeBalance(CtrlType *, GraphType *, float);
my_int ComputeMaxNodeGain(my_int, idxtype *, idxtype *, idxtype *);

/* srefine.c */
void Refine2WayNode(CtrlType *, GraphType *, GraphType *, float);
void Allocate2WayNodePartitionMemory(CtrlType *, GraphType *);
void Compute2WayNodePartitionParams(CtrlType *, GraphType *);
void Project2WayNodePartition(CtrlType *, GraphType *);

/* stat.c */
void ComputePartitionInfo(GraphType *, my_int, idxtype *);
void ComputePartitionInfoBipartite(GraphType *, my_int, idxtype *);
void ComputePartitionBalance(GraphType *, my_int, idxtype *, float *);
float ComputeElementBalance(my_int, my_int, idxtype *);

/* subdomains.c */
void Random_KWayEdgeRefineMConn(CtrlType *, GraphType *, my_int, float *, float, my_int, my_int);
void Greedy_KWayEdgeBalanceMConn(CtrlType *, GraphType *, my_int, float *, float, my_int);
void PrintSubDomainGraph(GraphType *, my_int, idxtype *);
void ComputeSubDomainGraph(GraphType *, my_int, idxtype *, idxtype *);
void EliminateSubDomainEdges(CtrlType *, GraphType *, my_int, float *);
void MoveGroupMConn(CtrlType *, GraphType *, idxtype *, idxtype *, my_int, my_int, my_int, idxtype *);
void EliminateComponents(CtrlType *, GraphType *, my_int, float *, float);
void MoveGroup(CtrlType *, GraphType *, my_int, my_int, my_int, idxtype *, idxtype *);

/* timing.c */
void InitTimers(CtrlType *);
void PrintTimers(CtrlType *);
double seconds(void);

/* util.c */
void errexit(char *,...);
#ifndef DMALLOC
my_int *imalloc(my_int, char *);
idxtype *idxmalloc(my_int, char *);
float *fmalloc(my_int, char *);
my_int *ismalloc(my_int, my_int, char *);
idxtype *idxsmalloc(my_int, idxtype, char *);
void *GKmalloc(my_int, char *);
#endif
/*void GKfree(void **,...); */
my_int *iset(my_int n, my_int val, my_int *x);
idxtype *idxset(my_int n, idxtype val, idxtype *x);
float *sset(my_int n, float val, float *x);
my_int iamax(my_int, my_int *);
my_int idxamax(my_int, idxtype *);
my_int idxamax_strd(my_int, idxtype *, my_int);
my_int samax(my_int, float *);
my_int samax2(my_int, float *);
my_int idxamin(my_int, idxtype *);
my_int samin(my_int, float *);
my_int idxsum(my_int, idxtype *);
my_int idxsum_strd(my_int, idxtype *, my_int);
void idxadd(my_int, idxtype *, idxtype *);
my_int charsum(my_int, char *);
my_int isum(my_int, my_int *);
float ssum(my_int, float *);
float ssum_strd(my_int n, float *x, my_int);
void sscale(my_int n, float, float *x);
float snorm2(my_int, float *);
float sdot(my_int n, float *, float *);
void saxpy(my_int, float, float *, my_int, float *, my_int);
void RandomPermute(my_int, idxtype *, my_int);
double drand48();
void srand48(long);
my_int ispow2(my_int);
void InitRandom(my_int);
my_int log2(my_int);










/***************************************************************
* Programs Directory
****************************************************************/

/* io.c */
void ReadGraph(GraphType *, char *, my_int *);
void WritePartition(char *, idxtype *, my_int, my_int);
void WriteMeshPartition(char *, my_int, my_int, idxtype *, my_int, idxtype *);
void WritePermutation(char *, idxtype *, my_int);
my_int CheckGraph(GraphType *);
idxtype *ReadMesh(char *, my_int *, my_int *, my_int *);
void WriteGraph(char *, my_int, idxtype *, idxtype *);

/* smbfactor.c */
void ComputeFillIn(GraphType *, idxtype *);
idxtype ComputeFillIn2(GraphType *, idxtype *);
my_int smbfct(my_int, idxtype *, idxtype *, idxtype *, idxtype *, idxtype *, my_int *, idxtype *, idxtype *, my_int *);


/***************************************************************
* Test Directory
****************************************************************/
void Test_PartGraph(my_int, idxtype *, idxtype *);
my_int VerifyPart(my_int, idxtype *, idxtype *, idxtype *, idxtype *, my_int, my_int, idxtype *);
my_int VerifyWPart(my_int, idxtype *, idxtype *, idxtype *, idxtype *, my_int, float *, my_int, idxtype *);
void Test_PartGraphV(my_int, idxtype *, idxtype *);
my_int VerifyPartV(my_int, idxtype *, idxtype *, idxtype *, idxtype *, my_int, my_int, idxtype *);
my_int VerifyWPartV(my_int, idxtype *, idxtype *, idxtype *, idxtype *, my_int, float *, my_int, idxtype *);
void Test_PartGraphmC(my_int, idxtype *, idxtype *);
my_int VerifyPartmC(my_int, my_int, idxtype *, idxtype *, idxtype *, idxtype *, my_int, float *, my_int, idxtype *);
void Test_ND(my_int, idxtype *, idxtype *);
my_int VerifyND(my_int, idxtype *, idxtype *);

