#include "cs_alya.h" // fvm_element_t 
#include "build_octree.h"
double __ERROR__ = 1.e-28; 

//---------------------------------------------------------------| TETRAS |---//
void
__tetra__(const double* coords, const int* vertices, const double* point, double* shapef)
{
  double _epsilon_denom = __ERROR__; 
  double isop_0, isop_1, isop_2;
  double t00, t10, t20, t01, t02, t03, t11, t12, t13, t21, t22, t23;
  double v01[3], v02[3], v03[3]; //, shapef[4];
  double vol6;
  int i, j, k;


  ple_coord_t  tetra_coords[4][3];
  for(i = 0; i < 4; i++)
  {
    ple_lnum_t   coord_idx = vertices[i] - 1; //<-NOTE: fortran numeration 
    for(j = 0; j < 3; j++) tetra_coords[i][j] = coords[(coord_idx * 3) + j];
  }

  // vol[i] = (<x_j,y_j,z_j> - <x_0,y_0,z_0>)_i = <x_j-x_0, y_j-y_0, z_j-z_0>_i = R_j - R_0 = R_0j 
  for(i = 0; i < 3; i++)
  {
    v01[i] = tetra_coords[1][i] - tetra_coords[0][i];
    v02[i] = tetra_coords[2][i] - tetra_coords[0][i];
    v03[i] = tetra_coords[3][i] - tetra_coords[0][i];
  }

  // 6*Vol = J = A.BxC = R_01 . R_02 x R_03
  vol6 = fabs(  v01[0] * (v02[1]*v03[2] - v02[2]*v03[1])
              - v02[0] * (v01[1]*v03[2] - v01[2]*v03[1])
              + v03[0] * (v01[1]*v02[2] - v01[2]*v02[1]));
  if(vol6 < _epsilon_denom) return;


  // Tetrahedral coordinates 
  t00  = - tetra_coords[0][0] + point[0];
  t10  = - tetra_coords[0][1] + point[1];
  t20  = - tetra_coords[0][2] + point[2];

  t01  = - tetra_coords[0][0] + tetra_coords[1][0];
  t02  = - tetra_coords[0][0] + tetra_coords[2][0];
  t03  = - tetra_coords[0][0] + tetra_coords[3][0];

  t11  = - tetra_coords[0][1] + tetra_coords[1][1];
  t12  = - tetra_coords[0][1] + tetra_coords[2][1];
  t13  = - tetra_coords[0][1] + tetra_coords[3][1];

  t21  = - tetra_coords[0][2] + tetra_coords[1][2];
  t22  = - tetra_coords[0][2] + tetra_coords[2][2];
  t23  = - tetra_coords[0][2] + tetra_coords[3][2];

  isop_0 = (  t00 * (t12*t23 - t13*t22)
            - t10 * (t02*t23 - t22*t03)
            + t20 * (t02*t13 - t12*t03)) / vol6;
  isop_1 = (- t00 * (t11*t23 - t13*t21)
            + t10 * (t01*t23 - t21*t03)
            - t20 * (t01*t13 - t03*t11)) / vol6;
  isop_2 = (  t00 * (t11*t22 - t21*t12)
            - t10 * (t01*t22 - t21*t02)
            + t20 * (t01*t12 - t11*t02)) / vol6;


  memset(shapef, 0.0, 4*sizeof(double));

  shapef[0] = 1.0 - isop_0 - isop_1 - isop_2;
  shapef[1] =       isop_0;
  shapef[2] =                isop_1;
  shapef[3] =                         isop_2;

}
//-----------------------------------------------------------------------||---//

//-----------------------------------------------------------------------||---//
static int
_inverse_3x3(double  m[3][3],
             double  b[3],
             double  x[3])
{
  double det, det_inv, x0, x1, x2;
  double _epsilon_denom = __ERROR__;

  det =   m[0][0]*(m[1][1]*m[2][2] - m[2][1]*m[1][2])
        - m[1][0]*(m[0][1]*m[2][2] - m[2][1]*m[0][2])
        + m[2][0]*(m[0][1]*m[1][2] - m[1][1]*m[0][2]);

  if(fabs(det) < _epsilon_denom) 
  {
    return 1;
/*
  if (CS_ABS(det) < _epsilon_denom)
    return 1;
*/
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
//-----------------------------------------------------------------------||---//

//-----------------------------------------------------------------------||---//
//
//   elt_type    <-- type of element
//   uvw[]       <-- parametric coordinates
//   shapef[]    --> barycenter's coordinates
//   deriv [][]  --> derivative of shape function
//
static void
_compute_shapef_3d(fvm_element_t  elt_type,
                   const double   uvw[3],
                   double         shapef[8],
                   double         deriv[8][3])

{
  memset(shapef, 0.0, 8*sizeof(double));

  switch (elt_type) 
  {
  case FVM_CELL_HEXA:
    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * uvw[2];
    shapef[5] = uvw[0] * (1.0 - uvw[1]) * uvw[2];
    shapef[6] = uvw[0] * uvw[1] * uvw[2];
    shapef[7] = (1.0 - uvw[0]) * uvw[1] * uvw[2];

    if (deriv != NULL) {
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
    break;

  case FVM_CELL_PRISM:
    shapef[0] = (1.0 - uvw[0] - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[2]);
    shapef[2] = uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0] - uvw[1]) * uvw[2];
    shapef[4] = uvw[0] * uvw[2];
    shapef[5] = uvw[1] * uvw[2];

    if (deriv != NULL) {
      deriv[0][0] = -(1.0 - uvw[2]);
      deriv[0][1] = -(1.0 - uvw[2]);
      deriv[0][2] = -(1.0 - uvw[0] - uvw[1]);
      deriv[1][0] =  (1.0 - uvw[2]);
      deriv[1][1] =  0.0;
      deriv[1][2] = -uvw[0];
      deriv[2][0] =  0.0;
      deriv[2][1] =  (1.0 - uvw[2]);
      deriv[2][2] = -uvw[1];
      deriv[3][0] = -uvw[2];
      deriv[3][1] = -uvw[2];
      deriv[3][2] =  (1.0 - uvw[0] - uvw[1]);
      deriv[4][0] =  uvw[2];
      deriv[4][1] =  0.0;
      deriv[4][2] =  uvw[0];
      deriv[5][0] =  0.0;
      deriv[5][1] =  uvw[2];
      deriv[5][2] =  uvw[1];
    }
    break;

  case FVM_CELL_PYRAM:
    shapef[0] = (1.0 - uvw[0]) * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[1] = uvw[0] * (1.0 - uvw[1]) * (1.0 - uvw[2]);
    shapef[2] = uvw[0] * uvw[1] * (1.0 - uvw[2]);
    shapef[3] = (1.0 - uvw[0]) * uvw[1] * (1.0 - uvw[2]);
    shapef[4] = uvw[2];

    if (deriv != NULL) {
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
      deriv[4][0] =  0.0;
      deriv[4][1] =  0.0;
      deriv[4][2] =  1.0;
    }
    break;

  default:
    printf("ERROR: FVM_CELL NOT FOUND!!\n\n");
    exit(0); 
  //bft_error(__FILE__, __LINE__, 0,
  //          _("_compute_shapef: unhandled element type %s\n"),
  //          fvm_element_type_name[elt_type]);
  }

}

//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
//
//    elt_type            <-- type of element
//    point_coords        <-- point coordinates
//    vertex_coords[]     <-- pointer to element vertex coordinates
//     tolerance           <-- location tolerance factor
//     uvw[]               --> parametric coordinates of point in element
//
static int
_compute_uvw(fvm_element_t      elt_type,
             const ple_coord_t  point_coords[],
             double             vertex_coords[8][3],
             double             tolerance,
             double             uvw[3])

{
  int i, j, n_elt_vertices, iter;
  int max_iter = 20;
  double dist;
  double a[3][3], b[3], x[3], shapef[8], dw[8][3];

  n_elt_vertices = fvm_nodal_n_vertices_element[elt_type];

  assert(   elt_type == FVM_CELL_HEXA
         || elt_type == FVM_CELL_PRISM
         || elt_type == FVM_CELL_PYRAM);

  // Use Newton-method to determine parametric coordinates and shape function

  for (i = 0; i < 3; i++)
    uvw[i] = 0.5;

  for (iter = 0; iter < max_iter; iter++) {

    _compute_shapef_3d(elt_type, uvw, shapef, dw);

    b[0] = - point_coords[0];
    b[1] = - point_coords[1];
    b[2] = - point_coords[2];

    for (i = 0; i < 3; i++) {
      for (j = 0; j < 3; j++)
        a[i][j] = 0.0;
    }

    for (i = 0; i < n_elt_vertices; i++) {

      b[0] += (shapef[i] * vertex_coords[i][0]);
      b[1] += (shapef[i] * vertex_coords[i][1]);
      b[2] += (shapef[i] * vertex_coords[i][2]);

      for (j = 0; j < 3; j++) {
        a[0][j] -= (dw[i][j] * vertex_coords[i][0]);
        a[1][j] -= (dw[i][j] * vertex_coords[i][1]);
        a[2][j] -= (dw[i][j] * vertex_coords[i][2]);
      }

    }

    if (_inverse_3x3(a, b, x))
      return 0;

    dist = 0.0;

    for (i = 0; i < 3; i++) {
      dist += x[i] * x[i];
      uvw[i] += x[i];
    }

    if (dist <= (tolerance * tolerance))
      return 1;

  }

  return 0;
}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
static void
_locate_in_cell_3d(
//                   ple_lnum_t          elt_num,
                   fvm_element_t       elt_type,
                   const ple_lnum_t    element_vertex_num[],
//                   const ple_lnum_t   *parent_vertex_num,
                   const ple_coord_t   vertex_coords[],
//                   const ple_coord_t   point_coords[],
                   const ple_coord_t   point[], 
//                   ple_lnum_t          n_points_in_extents,
//                   const ple_lnum_t    points_in_extents[],
                   double              tolerance, //,
                   //ple_lnum_t          location[],
                   //float               distance[]
  double *shapef
                   )
{
  int j, k, n_vertices;
  ple_lnum_t coord_idx, vertex_id;

  double uvw[3], dist, max_dist;
  double  _vertex_coords[8][3];

  n_vertices = fvm_nodal_n_vertices_element[elt_type];

  // Initialize local element coordinates copy
  for (vertex_id = 0; vertex_id < n_vertices; vertex_id++) 
  {
    coord_idx = element_vertex_num[vertex_id] - 1; //<-NOTE: fortran numeration 
    for (j = 0; j < 3; j++)  _vertex_coords[vertex_id][j] = vertex_coords[(coord_idx * 3) + j];
  }

  memset(shapef, 0.0, 8*sizeof(double));
  // Shape functions may be computed directly with tetrahedra 
  if (elt_type == FVM_CELL_TETRA)
  {
    __tetra__(vertex_coords, element_vertex_num, point, shapef);
  }
  else 
  {// For cell shapes other than tetrahedra, find shape functions iteratively 
      if(
        _compute_uvw(elt_type,
                       point, 
                       _vertex_coords,
                       tolerance,
                       uvw) 
         ) 
      {
          _compute_shapef_3d(elt_type, uvw, shapef, NULL);

/*
        max_dist = -1.0;

        // For hexahedra, no need to compute shape functions, as
        // the 3 parametric coordinates are simpler to use 
        if (elt_type == FVM_CELL_HEXA) 
        {
          _compute_shapef_3d(elt_type, uvw, shapef, NULL);

          for (j = 0; j < 3; j++)
          {
            dist = 2.*fabs(uvw[j] - 0.5);
            if(max_dist < dist) max_dist = dist;
          }

        }
        else 
        {// For pyramids ands prisms, we need to compute shape functions
          _compute_shapef_3d(elt_type, uvw, shapef, NULL);

          for (j = 0; j < n_vertices; j++)
          {
            dist = 2.*fabs(shapef[j] - 0.5);
            if(max_dist < dist) max_dist = dist;
          }

        }
*/
      }
      else
      {
        printf("ERROR!:");
        printf("       <u,v,w>: %f, %f, %f \n", uvw[0], uvw[1], uvw[2]); 
        exit(0); 
      }

  }

}
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
//-----------------------------------------------------------------------||---//


//-----------------------------------------------------------------------||---//
void 
get_shape_func
             (
//                   ple_lnum_t          elt_num,
                   fvm_element_t       elt_type,
                   const ple_lnum_t    element_vertex_num[],
//                   const ple_lnum_t   *parent_vertex_num,
                   const ple_coord_t   vertex_coords[],
                   const ple_coord_t   point_coords[],
//                   ple_lnum_t          n_points_in_extents,
//                   const ple_lnum_t    points_in_extents[],
                   double              tol, 
                   double             *shapef
             ) 
{
 _locate_in_cell_3d(
//                    elt_num,
                    elt_type,
                    element_vertex_num,
//                    parent_vertex_num,
                    vertex_coords,
                    point_coords,
//                    n_points_in_extents,
//                    points_in_extents,
                    tol, 
                    shapef
                   ); 
} 

//-----------------------------------------------------------------------||---//
