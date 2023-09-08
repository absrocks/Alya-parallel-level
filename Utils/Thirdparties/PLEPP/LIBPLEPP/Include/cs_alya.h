//#ifndef CS_ALYA_H
//#define HAVE_CONFIG_H
//#endif
 
#include "cs_config.h"
#include "cs_defs.h"
#include "cs_coupling.h"        // base/cs_coupling.cs_coupling_point_in_mesh
#include "fvm_defs.h"
#include "fvm_nodal.h"
#include "fvm_point_location.h" // fvm/fvm_point_location.fvm_point_location_nodal
#include "fvm_nodal_priv.h"
//#include "cs_coupling.h"        // base/cs_coupling.cs_coupling_point_in_mesh

/*
#include "code_saturne/cs_config.h"
#include "code_saturne/cs_defs.h"
#include "code_saturne/fvm_defs.h"
#include "code_saturne/fvm_nodal.h"
#include "code_saturne/fvm_point_location.h" // fvm/fvm_point_location.fvm_point_location_nodal
#include "code_saturne/cs_coupling.h"        // base/cs_coupling.cs_coupling_point_in_mesh
#include "code_saturne/fvm_nodal_priv.h"
//#include "code_saturne/fvm_nodal_append.h"
//#include "code_saturne/"
*/
#ifndef CS_ALYA_H
#define CS_ALYA_H 
//=======================================================================||===//
typedef enum 
{
  TRI03 = 10,
  QUA04 = 12,
  TET04 = 30,
  PYR05 = 32,
  PEN06 = 34,
  HEX08 = 37
} alya_element_t;


#ifdef __cplusplus
extern "C" {  // c/fortran visibility!!
#endif

void
get_shape_func(
                   //ple_lnum_t          elt_num,
                   fvm_element_t       elt_type,
                   const ple_lnum_t    element_vertex_num[],
                   //const ple_lnum_t   *parent_vertex_num,
                   const ple_coord_t   vertex_coords[],
                   const ple_coord_t   point_coords[],
                   //ple_lnum_t          n_points_in_extents,
                   //const ple_lnum_t    points_in_extents[],
                   double              tol, 
                   double* shapef
             ); 

#ifdef __cplusplus
}
#endif


#endif // COMMDOM_H
//=======================================================================||===//
//=======================================================================||===//

/*

  fvm_element_t  cell_type = FVM_CELL_HEXA;

  for(int type_id = 0; type_id < FVM_N_ELEMENT_TYPES; type_id++) 
  {
    n_elements_type[type_id] = 0;
    sections[type_id] = NULL;
  }
*/

/*
  SATURNE: ./src/fvm/fvm_defs.h

typedef enum 
{
  FVM_EDGE,               // Edge 
  FVM_FACE_TRIA,          // Triangle 
  FVM_FACE_QUAD,          // Quadrangle 
  FVM_FACE_POLY,          // Simple Polygon 
  FVM_CELL_TETRA,         // Tetrahedron 
  FVM_CELL_PYRAM,         // Pyramid 
  FVM_CELL_PRISM,         // Prism (pentahedron) 
  FVM_CELL_HEXA,          // Hexahedron (brick) 
  FVM_CELL_POLY,          // Simple Polyhedron (convex or quasi-convex)
  FVM_N_ELEMENT_TYPES     // Number of element types
} fvm_element_t;
*/


/*
def Alya_cell_types():
    # col=1: kernel/domain/elmtyp.f90        <<==
    # col=2: kernel/elsest/elsest_geogid.f90 <<=
    alya_types = []
    alya_types.append( ("TRI03", 10) )
    alya_types.append( ("QUA04", 12) )
    alya_types.append( ("TET04", 30) )
    alya_types.append( ("PYR05", 32) )
    alya_types.append( ("PEN06", 34) )
    alya_types.append( ("HEX08", 37) )

    return alya_types

    #self.salome2vtk = {}
    #for salome_type, vtk_type in zip(self.salome_types, self.VTK.vtk_types):
    #  self.salome2vtk[salome_type] = vtk_type
*/

/*
def Vtk_cell_types():
    vtk_types = []    #(vtk_type, vtk_id, vtkCell) <= linear cell types found in VTK 
    vtk_types.append( ("VTK_TRIANGLE",     5, vtk.vtkTriangle()   ) )
    vtk_types.append( ("VTK_QUAD",         9, vtk.vtkQuad()       ) )
    vtk_types.append( ("VTK_TETRA",       10, vtk.vtkTetra()      ) )
    vtk_types.append( ("VTK_PYRAMID",     14, vtk.vtkPyramid() ) )
    vtk_types.append( ("VTK_WEDGE",       13, vtk.vtkWedge()       ) )
    vtk_types.append( ("VTK_HEXAHEDRON",  12, vtk.vtkHexahedron() ) )

    return vtk_types
*/

