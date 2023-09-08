#include "outside_ple.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//
//   elt_type    <-- type of element
//   uvw[]       <-- parametric coordinates
//   shapef[]    --> barycenter's coordinates
//   deriv [][]  --> derivative of shape function
//
static void
_compute_shapef_3d(syr_cfd_element_t elt_type,
                   const double   uvw[3],
                   double         shapef[8],
                   double         deriv[8][3])

{
  assert(elt_type == FVM_CELL_HEXA);

  memset(shapef, 0.0, 8*sizeof(double));

  shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
  shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
  shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
  shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
  shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * uvw[2];
  shapef[5] = uvw[0] * (1.0 - uvw[1]) * uvw[2];
  shapef[6] = uvw[0] * uvw[1] * uvw[2];
  shapef[7] = (1.0 - uvw[0]) * uvw[1] * uvw[2];

  if (deriv != NULL)
  {
      deriv[0][0] = -(1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[1]) * (1.0 - uvw[2]);
      deriv[1][1] = -uvw[0] * (1.0 - uvw[2]);
      deriv[1][2] = -uvw[0] * (1.0 - uvw[1]);
      deriv[2][0] =  uvw[1] * (1.0 - uvw[2]);
      deriv[2][1] =  uvw[0] * (1.0 - uvw[2]);
      deriv[2][2] = -uvw[0] * uvw[1];
      deriv[3][0] = -uvw[1] * (1.0 - uvw[2]);
      deriv[3][1] =  (1.0 - uvw[0]) * (1.0 - uvw[2]);
      deriv[3][2] = -(1.0 - uvw[0]) * uvw[1];
      deriv[4][0] = -(1.0 - uvw[1]) * uvw[2];
      deriv[4][1] = -(1.0 - uvw[0]) * uvw[2];
      deriv[4][2] =  (1.0 - uvw[0]) * (1.0 - uvw[1]);
      deriv[5][0] =  (1.0 - uvw[1]) * uvw[2];
      deriv[5][1] = -uvw[0] * uvw[2];
      deriv[5][2] =  uvw[0] * (1.0 - uvw[1]);
      deriv[6][0] =  uvw[1] * uvw[2];
      deriv[6][1] =  uvw[0] * uvw[2];
      deriv[6][2] =  uvw[0] * uvw[1];
      deriv[7][0] = -uvw[1] * uvw[2];
      deriv[7][1] =  (1.0 - uvw[0]) * uvw[2];
      deriv[7][2] =  (1.0 - uvw[0]) * uvw[1];
  }

} // _compute_shapef_3d

//-----------------------------------------------------------------------||---//
static int
_inverse_3x3(double  m[3][3],
             double  b[3],
             double  x[3])
{
  double det, det_inv, x0, x1, x2;
  double _epsilon_denom = 1.e-28;

  det =   m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - m[1][0]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + m[2][0]*(m[0][1]*m[1][2] - m[1][1]*m[0][2]);

  if(fabs(det) < _epsilon_denom)
  {
    return 1;
  }
  else
    det_inv = 1./det;

  // Use local variables to ensure no aliasing
  x0 = (  b[0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - b[1]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + b[2]*(m[0][1]*m[1][2] - m[1][1]*m[0][2])) * det_inv;

  x1 = (  m[0][0]*(b[1]*m[2][2] - b[2]*m[1][2])
        - m[1][0]*(b[0]*m[2][2] - b[2]*m[0][2])
        + m[2][0]*(b[0]*m[1][2] - b[1]*m[0][2])) * det_inv;

  x2 = (  m[0][0]*(m[1][1]*b[2] - m[2][1]*b[1])
        - m[1][0]*(m[0][1]*b[2] - m[2][1]*b[0])
        + m[2][0]*(m[0][1]*b[1] - m[1][1]*b[0])) * det_inv;

  // Copy local variables to output

  x[0] = x0; x[1] = x1; x[2] = x2;

  return 0;
}

/*----------------------------------------------------------------------------
 * Calculation of natural coordinates of the element
 *
 *    elt_type            <-- type of element
 *    point_coords        <-- point coordinates
 *    vertex_coords[]     <-- pointer to element vertex coordinates
 *    tolerance           <-- location tolerance factor
 *    uvw[]               --> parametric coordinates of point in element
 *
 * returns.
 *    uvw[]
 *----------------------------------------------------------------------------*/
static int
_compute_uvw(syr_cfd_element_t  elt_type,
             const ple_coord_t  point_coords[],
             double             vertex_coords[8][3],
             double             tolerance,
             double             uvw[3])

{
  int i, j, n_elt_vertices, iter;
  int max_iter = 20;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  n_elt_vertices = syr_cfd_mesh_n_vertices_element[elt_type];

  assert(elt_type == FVM_CELL_HEXA); // || elt_type == FVM_CELL_PRISM || elt_type == FVM_CELL_PYRAM);

  // Use Newton-method to determine parametric coordinates and shape function
  for (i = 0; i < 3; i++) uvw[i] = 0.5;

  for (iter = 0; iter < max_iter; iter++)
  {
    _compute_shapef_3d(elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) a[i][j] = 0.0;

    for (i = 0; i < n_elt_vertices; i++)
    {
      b[0] += (shapef[i] * vertex_coords[i][0]);
      b[1] += (shapef[i] * vertex_coords[i][1]);
      b[2] += (shapef[i] * vertex_coords[i][2]);

      for (j = 0; j < 3; j++)
      {
        a[0][j] -= (dw[i][j] * vertex_coords[i][0]);
        a[1][j] -= (dw[i][j] * vertex_coords[i][1]);
        a[2][j] -= (dw[i][j] * vertex_coords[i][2]);
      }

    }

    if (_inverse_3x3(a, b, x)) return 0;

    dist = 0.0;

    for (i = 0; i < 3; i++)
    {
      dist   += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dist <= (tolerance * tolerance)) return 1;
  } // for

  return 0;
}

/*----------------------------------------------------------------------------
 * Calculation of shape functions iteratively
 *
 * This subroutine is executed for all the elements and points which are
 * potentially associated with the coupling.
 * Shape functions are used for localization and for the interpoaltion
 *
 *    elt_type             <-- type of element
 *    element_vertex_num[] <-- id (index) element vertices
 *    vertex_coords[]      <-- array of coordinates
 *    point[]              <-- coordinates point
 *    tolerance            <-- tolerance
 *    shapef               --> shape functions (weigth for the calculation)
 * returns.
 *    shape functions
 *----------------------------------------------------------------------------*/
void
_locate_in_cell_3d(
                   syr_cfd_element_t   elt_type,
                   const ple_lnum_t    element_vertex_num[],
                   const ple_coord_t   vertex_coords[],
                   const ple_coord_t   point[],
                   double              tolerance,
                   double*             shapef
                   )
{
  int j, k, n_vertices;
  ple_lnum_t coord_idx, vertex_id;

  double uvw[3], dist, max_dist;
  double  _vertex_coords[8][3];

  n_vertices = syr_cfd_mesh_n_vertices_element[elt_type];

  // Initialize local element coordinates copy
  for (vertex_id = 0; vertex_id < n_vertices; vertex_id++)
  {
    coord_idx = element_vertex_num[vertex_id] - 1; //<-NOTE: fortran numeration
    for (j = 0; j < 3; j++)  _vertex_coords[vertex_id][j] = vertex_coords[(coord_idx * 3) + j];
  }

  memset(shapef, 0.0, 8*sizeof(double));
  if (elt_type == SYR_CFD_TETRA)
  {
    // Shape functions may be computed directly with tetrahedra
    printf("LOCATE_IN_CELL_3D: ERROR WE HAVE A TETRA!!!");
    exit(0);
  }
  else
  {
    // For cell shapes other than tetrahedra, find shape functions iteratively
    if(_compute_uvw(elt_type, point, _vertex_coords, tolerance, uvw) )
    {
      max_dist = -1.0;

      // For hexahedra, no need to compute shape functions, as the 3 parametric coordinates are simpler to use
      if(elt_type == FVM_CELL_HEXA)
      {
        _compute_shapef_3d(elt_type, uvw, shapef, NULL);

        for (j = 0; j < 3; j++)
        {
          dist = 2.0 * fabs(uvw[j] - 0.5);
          if(max_dist < dist) max_dist = dist;
        }
      }
      else
      {
         // For pyramids ands prisms, we need to compute shape functions
         printf("LOCATE_IN_CELL_3D: PYRAMIDS AND PRISMS NOT IMPLEMENTED");
         exit(0);
      }

    }
    else
    {
      printf("LOCATE_IN_CELL_3D: ERROR calculation natural coordinates HEXA:");
      printf("       <u,v,w>: %f, %f, %f \n", uvw[0], uvw[1], uvw[2]);
      exit(0);
    }

  }
} // _locate_in_cell_3d

//-----------------------------------------------------------------------||---//
static void
__print_cell__(
                   ple_lnum_t         __elt_num,
                   syr_cfd_element_t   elt_type,
                   const ple_lnum_t    element_vertex_num[],
                   const ple_coord_t   vertex_coords[]
                   )
{
  int j, k, n_vertices;
  ple_lnum_t vertex_id;

  n_vertices = syr_cfd_mesh_n_vertices_element[elt_type];

  printf(" %d) ", __elt_num );
  for (vertex_id = 0; vertex_id < n_vertices; vertex_id++) printf(" %d ", element_vertex_num[vertex_id] );
  printf(" \n");

}
/*----------------------------------------------------------------------------
 * Main subroutine for hexahedron
 *
 * This subroutine locates a point inside an hexaedron
 *
 * returns.
 *    __location
 *    __distance
 *----------------------------------------------------------------------------*/
static void
_locate_in_hexa(ple_lnum_t         __elt_num,               //  input: for (i = 0; i < mesh->n_elements; i++)
                const ple_lnum_t   __element_vertex_num[],  //  input: vertex_num + i*stride
                const ple_coord_t  __vertex_coords[],       //  input: vertex_coords
                const ple_coord_t  __point_coords[],        //  input: vertex_coords_j ??
                ple_lnum_t         __n_points_in_extents,   //  input: _query_octree -> n_points_in_extents
                const ple_lnum_t   __points_in_extents[],   //  input: _query_octree -> points_in_extents
                double             __tolerance,             //  input: 1e-8
                ple_lnum_t         __location[],            // output:
                float              __distance[])            // output:
{
  double vol6;
  double dist, max_dist;
  int i, j, k;

/*
  double coord[3];
  for (i = 0; i < __n_points_in_extents; i++)
  {
    j = __points_in_extents[i];
    for(k=0; k<3; k++) coord[k] = __point_coords[i*3 + k];

    printf("%d) ", i);
    for(k=0; k<3; k++) printf("%f ", coord[k] );
    printf("\n ");
  }
// */

//  if(__n_points_in_extents) printf("%d) n_points_in_extents:%d \n ", __elt_num, __n_points_in_extents);

/*
  __print_cell__(
                     __elt_num,
                     FVM_CELL_HEXA,
                     __element_vertex_num,
                     __vertex_coords
                 );
*/

  for (k = 0; k < __n_points_in_extents; k++)
  {
    i = __points_in_extents[k];

    double coord[3];
    for(j=0; j<3; j++) coord[j] = __point_coords[i*3 + j];

    double shapef[8];
    //printf("<element>: %d\n", __elt_num);
    _locate_in_cell_3d(
                        FVM_CELL_HEXA,
                        __element_vertex_num,
                        __vertex_coords,
                        coord,
                        __tolerance,
                        shapef
                      );
    //
    // Criterion to find if the point is inside or outside the element
    //
    max_dist = -1.0;
    for (j = 0; j < 4; j++)
    {
      dist = 2.0 * PLE_ABS(shapef[j] - 0.5);
      if (max_dist < dist) max_dist = dist;
    }
    //
    // Localization (We only localize if the point is inside the element)
    //
    if( (max_dist > -0.5 && max_dist < (1. + 2.*__tolerance)) && (max_dist < __distance[i] || __distance[i] < 0))
    {
      __location[i] = __elt_num;
      __distance[i] = max_dist;
    }

  } // for

}  // _locate_in_hexa

