/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief External triangle-mesh data importer.

#include "triangle_mesh_importer.h"

#include <cstdio>
#include <cmath>
#include <limits>
#include <vector>

#include <nv/index/itriangle_mesh_subset.h>

#include "forwarding_logger.h"

namespace nv {
namespace index_common {

static bool dbg_calc_trimesh_bbox = false;

/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Moeller                             */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/

#define X 0
#define Y 1
#define Z 2

#define CROSS(dest,v1,v2) \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0]; 

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2) \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2]; 

#define FINDMINMAX(x0,x1,x2,min,max) \
  min = max = x0;   \
  if(x1<min) min=x1;\
  if(x1>max) max=x1;\
  if(x2<min) min=x2;\
  if(x2>max) max=x2;


/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)             \
    p0 = a*v0[Y] - b*v0[Z];                    \
    p2 = a*v2[Y] - b*v2[Z];                    \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fa, fb)              \
    p0 = a*v0[Y] - b*v0[Z];                    \
    p1 = a*v1[Y] - b*v1[Z];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];   \
    if(min>rad || max<-rad) return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)             \
    p0 = -a*v0[X] + b*v0[Z];                   \
    p2 = -a*v2[X] + b*v2[Z];                       \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fa, fb)              \
    p0 = -a*v0[X] + b*v0[Z];                   \
    p1 = -a*v1[X] + b*v1[Z];                       \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];   \
    if(min>rad || max<-rad) return 0;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)             \
    p1 = a*v1[X] - b*v1[Y];                    \
    p2 = a*v2[X] - b*v2[Y];                    \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
    if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fa, fb)              \
    p0 = a*v0[X] - b*v0[Y];                \
    p1 = a*v1[X] - b*v1[Y];                    \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
    if(min>rad || max<-rad) return 0;

// -- note: the load bbox can actually be different from the object bbox, should
// be ok if there are two boxes
//----------------------------------------------------------------------
Triangle_mesh_importer::Triangle_mesh_importer(const std::string& filename)
{
    m_filename = filename;

    m_configuration = std::string() +
        "importer=raw\n" +
        "input_file=" + m_filename + "\n";
}

//----------------------------------------------------------------------
// constructor
Triangle_mesh_importer::Triangle_mesh_importer()
{
}

//----------------------------------------------------------------------
// destructor
Triangle_mesh_importer::~Triangle_mesh_importer()
{
}

//----------------------------------------------------------------------
// execute importing data
//
// file format (todo):
//
// entire box <bounding box, 6x float32>
// total number of triangle_meshs <uint32>
//   radius <float32>
//   [triangle_mesh bbox <6x float32>
//    number of points in triangle_mesh <uint32>
//    ... 3d float32 vectors (x no. points)]
//   [...]
//

namespace
{

static inline mi::Sint32 planeBoxOverlap(mi::Float32 normal[3], mi::Float32 vert[3], mi::Float32 maxbox[3])    // -NJMP-
{
  mi::Sint32 q;
  mi::Float32 vmin[3],vmax[3],v;
  for(q=X;q<=Z;q++)
  {
    v=vert[q];                  // -NJMP-
    if(normal[q]>0.0f)
    {
        vmin[q]=-maxbox[q] - v; // -NJMP-
        vmax[q]= maxbox[q] - v; // -NJMP-
    }
    else
    {
        vmin[q]= maxbox[q] - v; // -NJMP-
        vmax[q]=-maxbox[q] - v; // -NJMP-
    }
  }
  if(DOT(normal,vmin)>0.0f) return 0;   // -NJMP-
  if(DOT(normal,vmax)>=0.0f) return 1;  // -NJMP-
  
  return 0;
}

static inline mi::Sint32 triBoxOverlap(mi::Float32 boxcenter[3],mi::Float32 boxhalfsize[3],
    mi::Float32 trivert_a[3],
    mi::Float32 trivert_b[3],
    mi::Float32 trivert_c[3])
{

  /*    use separating axis theorem to test overlap between triangle and box */
  /*    need to test for overlap in these directions: */
  /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
  /*       we do not even need to test these) */
  /*    2) normal of the triangle */
  /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
  /*       this gives 3x3=9 more tests */
   mi::Float32 v0[3],v1[3],v2[3];
   mi::Float32 min,max,p0,p1,p2,rad,fex,fey,fez;      // -NJMP- "d" local variable removed
   mi::Float32 normal[3],e0[3],e1[3],e2[3];

   /* This is the fastest branch on Sun */
   /* move everything so that the boxcenter is in (0,0,0) */
   SUB(v0,trivert_a,boxcenter);
   SUB(v1,trivert_b,boxcenter);
   SUB(v2,trivert_c,boxcenter);

   /* compute triangle edges */
   SUB(e0,v1,v0);      /* tri edge 0 */
   SUB(e1,v2,v1);      /* tri edge 1 */
   SUB(e2,v0,v2);      /* tri edge 2 */

   /* Bullet 3:  */
   /*  test the 9 tests first (this was faster) */
   fex = fabsf(e0[X]);
   fey = fabsf(e0[Y]);
   fez = fabsf(e0[Z]);
   AXISTEST_X01(e0[Z], e0[Y], fez, fey);
   AXISTEST_Y02(e0[Z], e0[X], fez, fex);
   AXISTEST_Z12(e0[Y], e0[X], fey, fex);

   fex = fabsf(e1[X]);
   fey = fabsf(e1[Y]);
   fez = fabsf(e1[Z]);
   AXISTEST_X01(e1[Z], e1[Y], fez, fey);
   AXISTEST_Y02(e1[Z], e1[X], fez, fex);
   AXISTEST_Z0(e1[Y], e1[X], fey, fex);

   fex = fabsf(e2[X]);
   fey = fabsf(e2[Y]);
   fez = fabsf(e2[Z]);
   AXISTEST_X2(e2[Z], e2[Y], fez, fey);
   AXISTEST_Y1(e2[Z], e2[X], fez, fex);
   AXISTEST_Z12(e2[Y], e2[X], fey, fex);

   /* Bullet 1: */
   /*  first test overlap in the {x,y,z}-directions */
   /*  find min, max of the triangle each direction, and test for overlap in */
   /*  that direction -- this is equivalent to testing a minimal AABB around */
   /*  the triangle against the AABB */

   /* test in X-direction */
   FINDMINMAX(v0[X],v1[X],v2[X],min,max);
   if(min>boxhalfsize[X] || max<-boxhalfsize[X]) return 0;

   /* test in Y-direction */
   FINDMINMAX(v0[Y],v1[Y],v2[Y],min,max);
   if(min>boxhalfsize[Y] || max<-boxhalfsize[Y]) return 0;

   /* test in Z-direction */
   FINDMINMAX(v0[Z],v1[Z],v2[Z],min,max);
   if(min>boxhalfsize[Z] || max<-boxhalfsize[Z]) return 0;

   /* Bullet 2: */
   /*  test if the box intersects the plane of the triangle */
   /*  compute plane equation of triangle: normal*x+d=0 */
   CROSS(normal,e0,e1);
   // -NJMP- (line removed here)
   if(!planeBoxOverlap(normal,v0,boxhalfsize)) return 0;    // -NJMP-

   return 1;   /* box and triangle overlaps */
}

static inline bool read_file_index_item(
    mi::Uint64& d,
    FILE*       s,
    bool        is64bit)
{
    if (is64bit) {
        return (fread(&d, sizeof(mi::Uint64), 1, s) == 1);
    }
    else {
        mi::Uint32 d32;
        if (fread(&d32, sizeof(mi::Uint32), 1, s) == 1) {
            d = d32;
            return true;
        }
        else {
            return false;
        }
    }
}

static inline bool read_file_index_items(
    mi::Uint64* d,
    mi::Uint64  n,
    FILE*       s,
    bool        is64bit)
{
    if (is64bit) {
        return (fread(d, sizeof(mi::Uint64), n, s) != n);
    }
    else {
        if (n > (std::numeric_limits<mi::Uint32>::max)()) {
            return false;
        }
        mi::Uint32* d32 = new mi::Uint32[n];
        if (fread(d32, sizeof(mi::Uint32), n, s) == n) {
            for (mi::Uint64 i = 0; i < n; ++i) {
                d[i] = d32[i];
            }
            delete[] d32;
            return true;
        }
        else {
            delete[] d32;
            return false;
        }
    }
}

} // namespace

mi::Size Triangle_mesh_importer::estimate(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    DEBUG_LOG << "The estimate still needs to be given!";
    return 8; // TODO!!
}

nv::index::IDistributed_data_subset* Triangle_mesh_importer::create(
    const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    VERBOSE_LOG << "invoked triangle mesh importer: " << mi::math::Bbox<mi::Sint32, 3>(bounding_box);

    // open and get the file magic
    FILE * p_file = fopen(m_filename.c_str(), "rb");
    if(p_file == 0)
    {
        ERROR_LOG << "Cannot open triangle-mesh data file: " << m_filename;
        return NULL;
    }

    // read the magic "ssvt"
    std::string magic_str;
    for(mi::Sint32 i = 0; (i < 4) && (!feof(p_file)); ++i)
    {
        mi::Sint32 const c = fgetc(p_file);
        if(c == -1)
        {
            ERROR_LOG << "Unexpected EOF when loading a triangle-mesh data file: " << m_filename;
            fclose(p_file);
            return NULL;
        }
        magic_str += static_cast< char >(c);
    }
    if(magic_str != "ssvt")
    {
        fclose(p_file);
        ERROR_LOG << m_filename << " is not a binary format triangle-mesh data file.";
        return NULL;
    }

    bool is_64bit_trimesh = false;

    // read version number
    mi::Sint32 version;
    if (fread(&version, sizeof(int), 1, p_file) != 1) 
    {
        fclose(p_file);
        ERROR_LOG << "Failed to read the triangle-mesh data version number [" << m_filename << "]";
        return NULL;
    }

    switch (version) {
        case 1:
            is_64bit_trimesh = false;
            break;
        case 2:
            is_64bit_trimesh = true;
            break;
        default:
            fclose(p_file);
            ERROR_LOG << "triangle-mesh data invalid version number [" << m_filename << "]";
            return NULL;
    }

    // read colormap mode
    mi::Sint32 has_vertex_colormap_indices;
    if (fread(&has_vertex_colormap_indices, sizeof(mi::Uint32), 1, p_file) != 1)
    {
        fclose(p_file);
        ERROR_LOG << "Failed to read the triangle-mesh data vertex colormap mode [" << m_filename << "]";
        return NULL;
    }

    // read bounding box of entire dataset
    mi::math::Bbox_struct<mi::Float32, 3> glob_bbox;
    if (fread(&glob_bbox.min.x, sizeof(mi::Float32), 6, p_file) != 6)
    {
        fclose(p_file);
        ERROR_LOG << "Failed to read the triangle-mesh data global bbox [" << m_filename << "]";
        return NULL;
    }

    // read positions
    mi::Uint64 nvectors;
    if (!read_file_index_item(nvectors, p_file, is_64bit_trimesh))
    {
        fclose(p_file);
        ERROR_LOG << "Failed to read the triangle-mesh data no. vectors [" << m_filename << "]";
        return NULL;
    }
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > vectors;
    vectors.resize(nvectors);
    if (nvectors > 0)
    {
        if (fread(&vectors[0], sizeof(mi::Float32)*3, nvectors, p_file) != nvectors)
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read triangle-mesh data vectors [" << m_filename << "]";
            return NULL;
        }
    }

    // read normals
    mi::Uint64 nnormals;
    if (!read_file_index_item(nnormals, p_file, is_64bit_trimesh))
    {
        fclose(p_file);
        ERROR_LOG << "Failed to read the triangle-mesh data no. normals [" << m_filename << "]";
        return NULL;
    }
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > normals;
    normals.resize(nnormals);
    if (nnormals > 0)
    {
        if (fread(&normals[0], sizeof(mi::Float32)*3, nnormals, p_file) != nnormals)
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read triangle-mesh data normals [" << m_filename << "]";
            return NULL;
        }
    }

    // read 2d-texture coordinates
    mi::Uint64 ntexcoords;
    if (!read_file_index_item(ntexcoords, p_file, is_64bit_trimesh))
    {
        fclose(p_file);
        ERROR_LOG << "Failed to read the triangle-mesh data no. texture coordinates [" << m_filename << "]";
        return NULL;
    }
    std::vector<mi::math::Vector_struct<mi::Float32, 2> > tex_coords;
    tex_coords.resize(ntexcoords);
    if (ntexcoords > 0)
    {
        if (fread(&tex_coords[0], sizeof(mi::Float32)*2, ntexcoords, p_file) != ntexcoords)
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read triangle-mesh data texture coordinates [" << m_filename << "]";
            return NULL;
        }
    }

    // read vertex colors
    mi::Uint64 nvc;
    if (!read_file_index_item(nvc, p_file, is_64bit_trimesh))
    {
        fclose(p_file);
        ERROR_LOG << "Failed to read the triangle-mesh data no. vertex colors [" << m_filename << "]";
        return NULL;
    }
    std::vector<mi::math::Color_struct> vertex_colors;
    vertex_colors.resize(nvc);
    if (nvc > 0)
    {
        if (fread(&vertex_colors[0], sizeof(mi::Float32)*4, nvc, p_file) != nvc)
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read triangle-mesh data vertex colors [" << m_filename << "]";
            return NULL;
        }
    }
    // read number of triangles contained in dataset
    mi::Uint64 ntriangles;
    if (!read_file_index_item(ntriangles, p_file, is_64bit_trimesh))
    {
        fclose(p_file);
        ERROR_LOG << "Failed to read the triangle-mesh data no. triangles [" << m_filename << "]";
        return NULL;
    }
    std::vector<mi::Uint64> pos_indices;
    pos_indices.resize(ntriangles * 3);
    if (ntriangles > 0)
    {
        if (!read_file_index_items(&pos_indices[0], ntriangles * 3, p_file, is_64bit_trimesh))
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read triangle-mesh data position indices [" << m_filename << "]";
            return NULL;
        }
    }

    std::vector<mi::Uint64> n_indices;
    if (nnormals)
    {
        n_indices.resize(ntriangles * 3);
        if (!read_file_index_items(&n_indices[0], ntriangles * 3, p_file, is_64bit_trimesh))
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read triangle-mesh data normal indices [" << m_filename << "]";
            return NULL;
        }
    }

    std::vector<mi::Uint64> t_indices;
    if (ntexcoords) {
        t_indices.resize(ntriangles * 3);
        if (!read_file_index_items(&t_indices[0], ntriangles * 3, p_file, is_64bit_trimesh))
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read triangle-mesh data texture coordinate indices [" << m_filename << "]";
            return NULL;
        }
    }

    std::vector<mi::Uint64> vc_indices;
    if (nvc)
    {
        vc_indices.resize(ntriangles * 3);
        if (!read_file_index_items(&vc_indices[0], ntriangles * 3, p_file, is_64bit_trimesh))
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read triangle-mesh data vertex color indices [" << m_filename << "]";
            return NULL;
        }
    }

    std::vector<mi::Uint64> vcmap_indices;
    if (has_vertex_colormap_indices)
    {
        vcmap_indices.resize(ntriangles * 3);
        if (!read_file_index_items(&vcmap_indices[0], ntriangles * 3, p_file, is_64bit_trimesh))
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read triangle-mesh data vertex colormap indices [" << m_filename << "]";
            return NULL;
        }
    }
    
    // read materials
    mi::Uint32 nmaterials;
    if (fread(&nmaterials, sizeof(mi::Uint32), 1, p_file) != 1)
    {
        fclose(p_file);
        ERROR_LOG << "Failed to read a triangle-mesh data no. materials [" << m_filename << "]";
        return NULL;
    }
    std::vector<mi::Uint32> materials;
    materials.resize(nmaterials);
    if (nmaterials > 0)
    {
        if (fread(&materials[0], sizeof(mi::Uint32), nmaterials, p_file) != nmaterials)
        {
            fclose(p_file);
            ERROR_LOG << "Failed to read the triangle-mesh data materials [" << m_filename << "]";
            return NULL;
        }
    }
    
    // read optional triangle flags
    mi::Uint32 ntriangle_flags;
    std::vector<mi::Uint32> triangle_flags;

    if (fread(&ntriangle_flags, sizeof(mi::Uint32), 1, p_file) == 1)
    {
        if (ntriangle_flags > 0)
        {
            triangle_flags.resize(ntriangle_flags);
            if (fread(&triangle_flags[0], sizeof(mi::Uint32), ntriangle_flags, p_file) != ntriangle_flags)
            {
                fclose(p_file);
                ERROR_LOG << "Failed to read the triangle-mesh data flags [" << m_filename << "]";
                return NULL;
            }
        }
    }

    fclose(p_file);
    
    mi::Float32 boxcenter[3], boxhalfsize[3];
    boxcenter[0]   = (bounding_box.max.x + bounding_box.min.x) / 2.f;
    boxcenter[1]   = (bounding_box.max.y + bounding_box.min.y) / 2.f;
    boxcenter[2]   = (bounding_box.max.z + bounding_box.min.z) / 2.f;
    boxhalfsize[0] = (bounding_box.max.x - bounding_box.min.x) / 2.f;
    boxhalfsize[1] = (bounding_box.max.y - bounding_box.min.y) / 2.f;
    boxhalfsize[2] = (bounding_box.max.z - bounding_box.min.z) / 2.f;

    typedef mi::math::Vector_struct<mi::Float32, 2> Vec2f;
    typedef mi::math::Vector_struct<mi::Float32, 3> Vec3f;
    typedef mi::math::Color_struct                  Colorf;
    typedef mi::math::Bbox_struct<mi::Float32, 3>   Bbox3f;

    static const mi::Uint32 invalid_remap_idx = ~0u;

    struct SR_tm_data
    {
        std::vector<Vec3f>      sr_vectors;
        std::vector<Vec3f>      sr_normals;
        std::vector<Vec2f>      sr_tex_coords;
        std::vector<Colorf>     sr_vertex_colors;
        std::vector<mi::Uint32> sr_pos_indices;
        std::vector<mi::Uint32> sr_normal_indices;
        std::vector<mi::Uint32> sr_tex_coord_indices;
        std::vector<mi::Uint32> sr_vertex_color_indices;
        std::vector<mi::Uint32> sr_vertex_colormap_indices;
        std::vector<mi::Uint16> sr_materials;
        std::vector<nv::index::ITriangle_mesh_subset::Triflags> sr_triangle_flags;

        // debug global bbox
        Bbox3f                  global_bbox;

        // global 64bit to 32bit per-subregion index translation
        std::vector<mi::Uint32> remap_vectors;
        std::vector<mi::Uint32> remap_normals;
        std::vector<mi::Uint32> remap_tex_coords;
        std::vector<mi::Uint32> remap_vertex_colors;

        // 32bit per-subregion to global 64bit index translation
        std::vector<mi::Uint64> remap_triangles;

        SR_tm_data(
            const mi::Uint64 nvectors,
            const mi::Uint64 nnormals,
            const mi::Uint64 ntexcoords,
            const mi::Uint64 nvc)
        {
            remap_vectors.resize(nvectors,      invalid_remap_idx);
            remap_normals.resize(nnormals,      invalid_remap_idx);
            remap_tex_coords.resize(ntexcoords, invalid_remap_idx);
            remap_vertex_colors.resize(nvc,     invalid_remap_idx);

            global_bbox.min = mi::math::Vector<mi::Float32, 3>( (std::numeric_limits<mi::Float32>::max)());
            global_bbox.max = mi::math::Vector<mi::Float32, 3>(-(std::numeric_limits<mi::Float32>::max)());
        }

        bool add_position(
            const mi::Uint64  g_pidx,
            const Vec3f*const positions)
        {
            mi::Uint32       l_pidx = remap_vectors[g_pidx];
            if (l_pidx == invalid_remap_idx) {
                if (sr_vectors.size() == (std::numeric_limits<mi::Uint32>::max)()) {
                    ERROR_LOG << "Subregion limit of 32-bit triangle-mesh index exceeded (positions): "
                              << "lower the per-subregion triangle count!";
                    return false;
                }
                l_pidx = static_cast<mi::Uint32>(sr_vectors.size());
                const Vec3f& v = positions[g_pidx];

                sr_vectors.push_back(v);
                remap_vectors[g_pidx] = l_pidx;
            }
            sr_pos_indices.push_back(l_pidx);
            return true;
        }
        void calc_bbox(
            const mi::Uint64  g_pidx,
            const Vec3f*const positions)
        {
            const Vec3f& v = positions[g_pidx];

            global_bbox.min.x = mi::math::min(v.x, global_bbox.min.x);
            global_bbox.min.y = mi::math::min(v.y, global_bbox.min.y);
            global_bbox.min.z = mi::math::min(v.z, global_bbox.min.z);
            global_bbox.max.x = mi::math::max(v.x, global_bbox.max.x);
            global_bbox.max.y = mi::math::max(v.y, global_bbox.max.y);
            global_bbox.max.z = mi::math::max(v.z, global_bbox.max.z);
        }
        bool add_position(
            const mi::Uint64  g_pidx,
            const Vec3f*const positions,
                  Bbox3f&     bbox)
        {
            mi::Uint32       l_pidx = remap_vectors[g_pidx];
            if (l_pidx == invalid_remap_idx) {
                if (sr_vectors.size() == (std::numeric_limits<mi::Uint32>::max)()) {
                    ERROR_LOG << "Subregion limit of 32-bit triangle-mesh index exceeded (positions): "
                              << "lower the per-subregion triangle count!";
                    return false;
                }
                l_pidx = static_cast<mi::Uint32>(sr_vectors.size());
                const Vec3f& v = positions[g_pidx];

                bbox.min.x = mi::math::min(v.x, bbox.min.x);
                bbox.min.y = mi::math::min(v.y, bbox.min.y);
                bbox.min.z = mi::math::min(v.z, bbox.min.z);
                bbox.max.x = mi::math::max(v.x, bbox.max.x);
                bbox.max.y = mi::math::max(v.y, bbox.max.y);
                bbox.max.z = mi::math::max(v.z, bbox.max.z);

                sr_vectors.push_back(v);
                remap_vectors[g_pidx] = l_pidx;
            }
            sr_pos_indices.push_back(l_pidx);
            return true;
        }
        bool add_normal(
            const mi::Uint64  g_nidx,
            const Vec3f*const normals)
        {
            mi::Uint32       l_nidx = remap_normals[g_nidx];
            if (l_nidx == invalid_remap_idx) {
                if (sr_normals.size() == (std::numeric_limits<mi::Uint32>::max)()) {
                    ERROR_LOG << "Subregion limit of 32-bit triangle-mesh index exceeded (normals): "
                              << "lower the per-subregion triangle count!";
                    return false;
                }
                l_nidx = static_cast<mi::Uint32>(sr_normals.size());
                sr_normals.push_back(normals[g_nidx]);
                remap_normals[g_nidx] = l_nidx;
            }
            sr_normal_indices.push_back(l_nidx);
            return true;
        }
        bool add_tex_coord(
            const mi::Uint64  g_tidx,
            const Vec2f*const tex_coords)
        {
            mi::Uint32       l_tidx = remap_tex_coords[g_tidx];
            if (l_tidx == invalid_remap_idx) {
                if (sr_tex_coords.size() == (std::numeric_limits<mi::Uint32>::max)()) {
                    ERROR_LOG << "Subregion limit of 32-bit triangle-mesh index exceeded (texture coordinates): "
                              << "lower the per-subregion triangle count!";
                    return false;
                }
                l_tidx = static_cast<mi::Uint32>(sr_tex_coords.size());
                sr_tex_coords.push_back(tex_coords[g_tidx]);
                remap_tex_coords[g_tidx] = l_tidx;
            }
            sr_tex_coord_indices.push_back(l_tidx);
            return true;
        }
        bool add_vertex_color(
            const mi::Uint64   g_vcidx,
            const Colorf*const vertex_colors)
        {
            mi::Uint32       l_vcidx = remap_vertex_colors[g_vcidx];
            if (l_vcidx == invalid_remap_idx) {
                if (sr_vertex_colors.size() == (std::numeric_limits<mi::Uint32>::max)()) {
                    ERROR_LOG << "Subregion limit of 32-bit triangle-mesh index exceeded (vertex colors): "
                              << "lower the per-subregion triangle count!";
                    return false;
                }
                l_vcidx = static_cast<mi::Uint32>(sr_vertex_colors.size());
                sr_vertex_colors.push_back(vertex_colors[g_vcidx]);
                remap_vertex_colors[g_vcidx] = l_vcidx;
            }
            sr_vertex_color_indices.push_back(l_vcidx);
            return true;
        }
        bool add_vertex_cmap(
            const mi::Uint64   g_vcmap_idx)
        {
            if (g_vcmap_idx > (std::numeric_limits<mi::Uint32>::max)()) {
                ERROR_LOG << "Subregion limit of 32-bit triangle-mesh index exceeded (vertex colormap): "
                          << "lower the per-subregion triangle count!";
                return false;
            }
            // no remapping, colormap index
            sr_vertex_colormap_indices.push_back(static_cast<mi::Uint32>(g_vcmap_idx));
            return true;
        }
        void add_material(
            const mi::Uint32 m)
        {
            sr_materials.push_back(static_cast<mi::Uint16>(m));
        }
        void add_triangle_flags(
            const nv::index::ITriangle_mesh_subset::Triflags f)
        {
            sr_triangle_flags.push_back(f);
        }

        void add_remap_triangle(
            const mi::Uint64 t)
        {
            remap_triangles.push_back(t);
        }
    };

    SR_tm_data sr_tmesh(nvectors, nnormals, ntexcoords, nvc);

    // put only triangles into mesh which are crossing/inside voxel
    Bbox3f bbox;
    bbox.min = mi::math::Vector<mi::Float32, 3>( (std::numeric_limits<mi::Float32>::max)());
    bbox.max = mi::math::Vector<mi::Float32, 3>(-(std::numeric_limits<mi::Float32>::max)());

    if (dbg_calc_trimesh_bbox) {
        // debug, calculate global bbox
        for (mi::Uint64 i=0; i < ntriangles; i++)
        {
            for (mi::Uint32 k=0; k < 3; k++) {
                const mi::Uint64 g_pidx = pos_indices[i*3+k];
                sr_tmesh.calc_bbox(g_pidx, &vectors.front());
            }
        }
    }

    std::vector<unsigned char> added;
    added.resize(ntriangles, 0);
    for (mi::Uint64 i=0; i < ntriangles; i++)
    {
        if (triBoxOverlap(boxcenter, boxhalfsize,
            &vectors[pos_indices[i*3+0]].x,
            &vectors[pos_indices[i*3+1]].x,
            &vectors[pos_indices[i*3+2]].x))
        {
            // add triangle
            for (mi::Uint32 k=0; k < 3; k++) {
                const mi::Uint64 g_pidx = pos_indices[i*3+k];
                sr_tmesh.add_position(g_pidx, &vectors.front(), bbox);
                if (nnormals) {
                    const mi::Uint64 g_nidx = n_indices[i*3+k];
                    sr_tmesh.add_normal(g_nidx, &normals.front());
                }
                if (ntexcoords) {
                    const mi::Uint64 g_tidx = t_indices[i*3+k];
                    sr_tmesh.add_tex_coord(g_tidx, &tex_coords.front());
                }
                if (nvc) {
                    const mi::Uint64 g_vcidx = vc_indices[i*3+k];
                    sr_tmesh.add_vertex_color(g_vcidx, &vertex_colors.front());
                }
                if (has_vertex_colormap_indices) {
                    const mi::Uint64 g_vcmap_idx = vcmap_indices[i*3+k];
                    sr_tmesh.add_vertex_cmap(g_vcmap_idx);
                }
            }

            if (!(materials.empty())) {
                sr_tmesh.add_material(i < materials.size() ? materials[i] : materials[0]);
            }
                
            if (!triangle_flags.empty()) {
                sr_tmesh.add_triangle_flags((nv::index::ITriangle_mesh_subset::Triflags) triangle_flags[i]);
            }
                
            sr_tmesh.add_remap_triangle(i);
            added[i] = 1;
        }
    }

    if (nnormals == 0)
    {
        // need to compute normals: add 1-ring around border
        std::vector<mi::Uint64> one_ring;

        for (mi::Uint64 i=0; i < ntriangles; i++)
        {
            if (sr_tmesh.remap_vectors[pos_indices[i*3+0]] != invalid_remap_idx ||
                sr_tmesh.remap_vectors[pos_indices[i*3+1]] != invalid_remap_idx || // a vertex of this triangle was remapped,
                sr_tmesh.remap_vectors[pos_indices[i*3+2]] != invalid_remap_idx)
            {   // i.e. a triangle sharing this vertex has been added
                // check whether this triangle was added already
                if (!added[i])
                    one_ring.push_back(i);
            }
        }
        for (mi::Uint64 j=0; j < one_ring.size(); j++)
        {
            mi::Uint64 i = one_ring[j];
            // add triangle
            for (mi::Uint32 k=0; k < 3; k++) {
                const mi::Uint64 g_pidx = pos_indices[i*3+k];
                sr_tmesh.add_position(g_pidx, &vectors.front(), bbox);
                if (nnormals) {
                    const mi::Uint64 g_nidx = n_indices[i*3+k];
                    sr_tmesh.add_normal(g_nidx, &normals.front());
                }
                if (ntexcoords) {
                    const mi::Uint64 g_tidx = t_indices[i*3+k];
                    sr_tmesh.add_tex_coord(g_tidx, &tex_coords.front());
                }
                if (nvc) {
                    const mi::Uint64 g_vcidx = vc_indices[i*3+k];
                    sr_tmesh.add_vertex_color(g_vcidx, &vertex_colors.front());
                }
                if (has_vertex_colormap_indices) {
                    const mi::Uint64 g_vcmap_idx = vcmap_indices[i*3+k];
                    sr_tmesh.add_vertex_cmap(g_vcmap_idx);
                }
            }

            if (!(materials.empty())) {
                sr_tmesh.add_material(i < materials.size() ? materials[i] : materials[0]);
            }
                
            if (!triangle_flags.empty()) {
                sr_tmesh.add_triangle_flags((nv::index::ITriangle_mesh_subset::Triflags) triangle_flags[i]);
            }
                
            sr_tmesh.add_remap_triangle(i);
        }
    }

    bool single_material = !sr_tmesh.sr_materials.empty();
    for (size_t i=1; i < sr_tmesh.sr_materials.size(); i++)
        if (sr_tmesh.sr_materials[i] != sr_tmesh.sr_materials[0])
            single_material = false;
    if (single_material)
        sr_tmesh.sr_materials.resize(1);
        
    if (dbg_calc_trimesh_bbox)
    {
        INFO_LOG << "triangle-mesh import job: loaded triangle mesh with global bbox " << sr_tmesh.global_bbox;
    }

    // ---------------------------------------------------------------------------------------------------
    // Create submesh (must be derived from ITriangle_mesh_subset)
    // and returns mesh to the caller. The caller takes ownership!
    //
    mi::base::Handle<nv::index::ITriangle_mesh_subset> submesh(
        factory->create<nv::index::ITriangle_mesh_subset>());
    if(submesh.is_valid_interface() && !sr_tmesh.sr_vectors.empty())
    {
        const bool valid_trianglemesh = submesh->initialize(
            bbox,
            sr_tmesh.sr_vectors.empty()                 ? 0 : &sr_tmesh.sr_vectors[0],                    sr_tmesh.sr_vectors.size(),
            sr_tmesh.sr_pos_indices.empty()             ? 0 : &sr_tmesh.sr_pos_indices[0],                
            sr_tmesh.sr_pos_indices.size(),

            sr_tmesh.remap_triangles.empty()            ? 0 : &sr_tmesh.remap_triangles[0],               sr_tmesh.remap_triangles.size(),
            
            sr_tmesh.sr_normals.empty()                 ? 0 : &sr_tmesh.sr_normals[0],                    sr_tmesh.sr_normals.size(),
            sr_tmesh.sr_tex_coords.empty()              ? 0 : &sr_tmesh.sr_tex_coords[0],                 sr_tmesh.sr_tex_coords.size(),
            sr_tmesh.sr_vertex_colors.empty()           ? 0 : &sr_tmesh.sr_vertex_colors[0],              sr_tmesh.sr_vertex_colors.size(),
            
            sr_tmesh.sr_normal_indices.empty()          ? 0 : &sr_tmesh.sr_normal_indices[0],             
            sr_tmesh.sr_tex_coord_indices.empty()       ? 0 : &sr_tmesh.sr_tex_coord_indices[0],          
            sr_tmesh.sr_vertex_color_indices.empty()    ? 0 : &sr_tmesh.sr_vertex_color_indices[0],       
            sr_tmesh.sr_vertex_colormap_indices.empty() ? 0 : &sr_tmesh.sr_vertex_colormap_indices[0],
            
            sr_tmesh.sr_materials.empty()               ? 0 : &sr_tmesh.sr_materials[0],                  sr_tmesh.sr_materials.size(),
            sr_tmesh.sr_triangle_flags.empty()          ? 0 : &sr_tmesh.sr_triangle_flags[0],             sr_tmesh.sr_triangle_flags.size()
            );

        if(valid_trianglemesh)
        {
            submesh->retain();  // since the handle will be out of scope.
            return submesh.get();
        }
    }

    return NULL;
}

//----------------------------------------------------------------------
mi::base::Uuid Triangle_mesh_importer::subset_id() const
{
    return nv::index::ITriangle_mesh_subset::IID();
}

//----------------------------------------------------------------------
const char* Triangle_mesh_importer::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
void Triangle_mesh_importer::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    mi::Uint32 nb_elements = mi::Uint32(m_filename.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_filename.c_str()), nb_elements);
}

//----------------------------------------------------------------------
void Triangle_mesh_importer::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_filename.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_filename[0]), nb_elements);
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
