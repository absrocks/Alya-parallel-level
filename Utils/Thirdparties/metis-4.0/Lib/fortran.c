/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * fortran.c
 *
 * This file contains code for the fortran to C interface
 *
 * Started 8/19/97
 * George
 *
 * $Id: fortran.c,v 1.1 1998/11/27 17:59:14 karypis Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function changes the numbering to start from 0 instead of 1
**************************************************************************/
void Change2CNumbering(my_int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  my_int i, nedges;

  for (i=0; i<=nvtxs; i++)
    xadj[i]--;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]--;
}

/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void Change2FNumbering(my_int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *vector)
{
  my_int i, nedges;

  for (i=0; i<nvtxs; i++)
    vector[i]++;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;
}

/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void Change2FNumbering2(my_int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  my_int i, nedges;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;
}



/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void Change2FNumberingOrder(my_int nvtxs, idxtype *xadj, idxtype *adjncy, idxtype *v1, idxtype *v2)
{
  my_int i, nedges;

  for (i=0; i<nvtxs; i++) {
    v1[i]++;
    v2[i]++;
  }

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;

}



/*************************************************************************
* This function changes the numbering to start from 0 instead of 1
**************************************************************************/
void ChangeMesh2CNumbering(my_int n, idxtype *mesh)
{
  my_int i;

  for (i=0; i<n; i++)
    mesh[i]--;

}


/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void ChangeMesh2FNumbering(my_int n, idxtype *mesh, my_int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  my_int i, nedges;

  for (i=0; i<n; i++)
    mesh[i]++;

  nedges = xadj[nvtxs];
  for (i=0; i<nedges; i++)
    adjncy[i]++;

  for (i=0; i<=nvtxs; i++)
    xadj[i]++;

}


/*************************************************************************
* This function changes the numbering to start from 1 instead of 0
**************************************************************************/
void ChangeMesh2FNumbering2(my_int n, idxtype *mesh, my_int ne, my_int nn, idxtype *epart, idxtype *npart)
{
  my_int i, nedges;

  for (i=0; i<n; i++)
    mesh[i]++;

  for (i=0; i<ne; i++)
    epart[i]++;

  for (i=0; i<nn; i++)
    npart[i]++;

}

