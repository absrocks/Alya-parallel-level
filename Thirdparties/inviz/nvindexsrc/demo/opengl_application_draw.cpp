/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "opengl_application_draw.h"

#include <cassert>
#include <vector>
#include <string>
#include <sstream>

#include <GL/glew.h>
#include <GL/gl.h>

#include "common/forwarding_logger.h"
#include "common/string_dict.h"

#include "nvindex_appdata.h"
#include "opengl_application_buffer.h"
#include "opengl_drawing_utilities.h"
#include "polygonal_cylinder_gen.h"

//----------------------------------------------------------------------

/// vector 3 float
typedef mi::math::Vector< mi::Float32, 3 > Float32_3;
/// vector 3 Sint32
typedef mi::math::Vector< mi::Sint32, 3 >  Sint32_3;
/// vector 4 Uint8 (for color)
typedef mi::math::Vector< mi::Uint8, 4>    Uint8_4;

/// application sample triangle mesh
class AppGLTriMesh
{
public:
    /// default constructor
    AppGLTriMesh()
    {
        m_vertex_list.clear();
        m_face_list.clear();
    }
    /// destructor
    virtual ~AppGLTriMesh()
    {
        m_vertex_list.clear();
        m_face_list.clear();
    }

    /// get string representation of this object.
    /// \return trimesh information
    std::string to_string() const {
        std::stringstream sstr;
        sstr << "# vertices: " << m_vertex_list.size()
             << ", # faces: "  << m_face_list.size();
        return sstr.str();
    }

public:
    // draw

    /// point mode
    void draw_points();
    /// wireframe mode
    void draw_wireframe();
    /// solid base color mode
    void draw_solid_color();
    /// solid flat mode
    void draw_solid_flat();


public:
    /// triangle vertices
    std::vector< Float32_3 > m_vertex_list;
    /// triangle face index
    std::vector< Sint32_3  > m_face_list;

private:
    // use default copy constructor.
    // AppGLTriMesh(AppGLTriMesh const &);
    // use default operator=.
    // AppGLTriMesh const & operator=(AppGLTriMesh const &);
};

//----------------------------------------------------------------------
void AppGLTriMesh::draw_points()
{
    DEBUG_LOG << "NIN: draw_points";
}

//----------------------------------------------------------------------
void AppGLTriMesh::draw_wireframe()
{
    DEBUG_LOG << "NIN: draw_wireframe";
}

//----------------------------------------------------------------------
void AppGLTriMesh::draw_solid_color()
{
    glBegin(GL_TRIANGLES);
    {
        mi::Uint32 fsize = m_face_list.size();
        for(mi::Uint32 fi = 0; fi < fsize; ++fi){
            glColor4f(0.89f, 0.2f, 0.13f, 1.f);

            assert(fi < m_face_list.size());
            Sint32_3 const fidx_list = m_face_list[fi];
            for(mi::Sint32 i = 0; i < 3; ++i){
                assert(fidx_list[i] < static_cast< mi::Sint32 >(m_vertex_list.size()));
                Float32_3 pos = m_vertex_list[fidx_list[i]];
                glVertex3f(pos[0], pos[1], pos[2]);
            }
        }
    }
    glEnd();
}

//----------------------------------------------------------------------
void AppGLTriMesh::draw_solid_flat()
{
    DEBUG_LOG << "NIN: draw_solid_flat";
}

//----------------------------------------------------------------------
/// application primitive correction with OpenGL draw().
///
/// singleton
class AppGLPrimitive
{
public:
    /// draw mode
    enum base_mode_bitmap_e {
        /// invalid mode
        No_mode        = 0x0000,
        /// points
        Points         = 0x0002,
        /// wireframe
        Wireframe      = 0x0004,
        /// solid base color
        Solid_color    = 0x0010,
        /// solid flat shading
        Solid_flat     = 0x0020,
        /// solid Gouraud shading
        Solid_Gouraud  = 0x0040,
        /// solid texture
        Solid_texture  = 0x0080,
        /// bounding box
        Bounding_box   = 0x0100,
        /// Pick mode for OpenGL selection
        Pick_mode      = 0x0200,
        /// Custom mode starts
        Custom_bottom  = 0x1000
    };


public:
    /// get the instance
    /// \return colormap manager singleton instance
    static AppGLPrimitive * instance()
    {
        if(G_p_app_primitive == 0){
            G_p_app_primitive = new AppGLPrimitive();
        }
        return G_p_app_primitive;
    }

    /// delete the singleton. (For unit test purpose. I don't want to
    /// confuse a memory checker.)
    static void delete_instance()
    {
        if(G_p_app_primitive != 0){
            delete G_p_app_primitive;
            G_p_app_primitive = 0;
        }
    }

private:
    // singleton instance
    static AppGLPrimitive * G_p_app_primitive;

public:
    /// destructor
    virtual ~AppGLPrimitive()
    {
        m_app_trimesh_vec.clear();
    }

public:
    //------------------------------------------------------------
    // primitive creation
    //------------------------------------------------------------

    /// create a primitives (a triangle) for test
    void create_primitve_test_triangle();


public:
    /// get number of meshes
    /// \return number of meshes
    mi::Sint32 get_nb_mesh() const
    {
        return static_cast< mi::Sint32 >(m_app_trimesh_vec.size());
    }

    /// peek a mesh at index
    /// \param[in] idx index of the mesh
    /// \return pointer to the mesh of idx
    AppGLTriMesh * peek_mesh(mi::Sint32 idx)
    {
        assert((0 <= idx) && (idx < this->get_nb_mesh()));
        return &(m_app_trimesh_vec.at(idx));
    }

    /// append a mesh
    /// \param[in] tmesh a triangle mesh. Note: this is copied.
    void append_mesh(AppGLTriMesh const & tmesh)
    {
        m_app_trimesh_vec.push_back(tmesh);
    }

    /// clear all the mesh
    void clear()
    {
        m_app_trimesh_vec.clear();
    }

    /// has a mesh?
    /// \return true when there are some meshes. false after clear().
    bool has_mesh() const
    {
        return !(m_app_trimesh_vec.empty());
    }

    /// get string representation of this object.
    /// \return primitive information
    std::string to_string() const
    {
        std::stringstream sstr;
        sstr << "# meshes:   " << m_app_trimesh_vec.size() << "\n";

        mi::Uint32 psize = m_app_trimesh_vec.size();
        for(mi::Uint32 i = 0; i < psize; ++i)
        {
            sstr << m_app_trimesh_vec.at(i).to_string();
            if((i+1) << psize){
                sstr << "\n";
            }
        }
        return sstr.str();
    }

public:
    // draw related functions

    /// set opengl draw mode
    /// \param[in] draw_mode opengl draw mode (not all are supported)
    void set_draw_mode(mi::Uint32 draw_mode)
    {
        m_draw_mode = draw_mode;
    }
    /// get opengl draw mode
    /// \return OpenGL draw mode
    mi::Uint32 get_draw_mode() const
    {
        return m_draw_mode;
    }

    /// draw all the mesh
    void draw_all_mesh();

private:
    /// set light enable
    /// \param[in] is_light_on light status to be set. true when light
    /// is on.
    void set_light_enable(GLboolean const is_light_on)
    {
        if(is_light_on == GL_TRUE){
            glEnable(GL_LIGHTING);
        }else{
            glDisable(GL_LIGHTING);
        }
    }

    /// draw a mesh
    void draw_a_mesh(mi::Uint32 draw_mode, AppGLTriMesh * p_tmesh);

private:
    /// application triangle meshes
    std::vector< AppGLTriMesh > m_app_trimesh_vec;
    /// draw mode
    mi::Uint32 m_draw_mode;

private:
    /// default constructor
    AppGLPrimitive()
        :
        m_draw_mode(Solid_color)
    {
        m_app_trimesh_vec.clear();
    }

private:
    /// copy constructor. prohibit until proved useful.
    AppGLPrimitive(AppGLPrimitive const &);
    /// operator=. prohibit until proved useful.
    AppGLPrimitive const & operator=(AppGLPrimitive const &);
};

/// singleton instance
AppGLPrimitive * AppGLPrimitive::G_p_app_primitive = 0;

//----------------------------------------------------------------------
void AppGLPrimitive::create_primitve_test_triangle()
{
    this->clear();

    // here is a simple one triangle.
    AppGLTriMesh sample_mesh;
    sample_mesh.m_vertex_list.push_back(Float32_3( 1.0f, 0.0f,    0.0f));
    sample_mesh.m_vertex_list.push_back(Float32_3( 1.0f, 1.7321f, 0.0f));
    sample_mesh.m_vertex_list.push_back(Float32_3(-1.0f, 0.0f,    0.0f));

    mi::Float32 scale_factor = 100.0f;
    for(std::vector< Float32_3 >::iterator vi = sample_mesh.m_vertex_list.begin();
        vi != sample_mesh.m_vertex_list.end(); ++vi)
    {
        (*vi) = scale_factor * (*vi);
    }

    sample_mesh.m_face_list.push_back(Sint32_3(0, 1, 2));

    DEBUG_LOG << this->to_string();

    this->append_mesh(sample_mesh);
}

//----------------------------------------------------------------------
void AppGLPrimitive::draw_all_mesh()
{
    // assume all meshes are drawn with the same mode.

    // push lighting/shading status
    const GLboolean is_light_on = glIsEnabled(GL_LIGHTING);
    GLint shading_model;
    glGetIntegerv(GL_SHADE_MODEL, &shading_model);

    // set light
    // NIN FIXME

    // push materials
    GLfloat push_diffuse_color_front[4];
    GLfloat push_diffuse_color_back[4];
    glGetMaterialfv(GL_FRONT, GL_DIFFUSE, push_diffuse_color_front);
    glGetMaterialfv(GL_BACK,  GL_DIFFUSE, push_diffuse_color_back);

    // set materials
    GLfloat diffuse_color[4] = { 1.0f, 0.0f, 0.0f, 1.0f, };
    glMaterialfv(GL_FRONT, GL_DIFFUSE, diffuse_color);
    glMaterialfv(GL_BACK,  GL_DIFFUSE, diffuse_color);

    {
        // draw meshes
        for(mi::Sint32 i = 0; i < this->get_nb_mesh(); ++i){
            this->draw_a_mesh(this->get_draw_mode(), this->peek_mesh(i));
        }
    }

    // pop materials
    glMaterialfv(GL_FRONT, GL_DIFFUSE, push_diffuse_color_front);
    glMaterialfv(GL_BACK,  GL_DIFFUSE, push_diffuse_color_back);

    // pop lighting/shading status
    this->set_light_enable(is_light_on);
    glShadeModel(shading_model);

    glFlush();
}

//----------------------------------------------------------------------
void AppGLPrimitive::draw_a_mesh(mi::Uint32 draw_mode, AppGLTriMesh * p_tmesh)
{
    if((draw_mode & Points) != 0){
        p_tmesh->draw_points();
    }

    if((draw_mode & Wireframe) != 0){
        p_tmesh->draw_wireframe();
    }

    if((draw_mode & Solid_color) != 0){
        glDisable(GL_LIGHTING); // always off for this mode
        p_tmesh->draw_solid_color();
    }

    if((draw_mode & Solid_flat) != 0){
        this->set_light_enable(GL_TRUE);
        p_tmesh->draw_solid_flat();
    }
}

//======================================================================
//----------------------------------------------------------------------
void gl_application_initialize(Opengl_application_buffer * p_app_buf,
                               const mi::math::Vector<mi::Sint32,2> & window_resolution)
{
    assert(p_app_buf != 0);
    const mi::math::Vector_struct< mi::Sint32, 2 > buff_size = window_resolution;
    p_app_buf->resize_buffer(buff_size);

    mi::Sint32 depth_bits = 0;
    glGetIntegerv(GL_DEPTH_BITS, &depth_bits);
    p_app_buf->set_z_buffer_precision(depth_bits);

    AppGLPrimitive::instance()->clear();
}

//----------------------------------------------------------------------
void gl_application_draw_object(
    Opengl_application_buffer*  p_app_buf)
{
    AppGLPrimitive::instance()->draw_all_mesh();

    assert(p_app_buf->is_buffer_allocated());
    mi::math::Vector< mi::Sint32, 2 > const res = p_app_buf->get_resolution();
    GLsizei const width  = res.x;
    GLsizei const height = res.y;

    // read the front buffer
    //glReadBuffer(GL_FRONT);

    // get z-buffer
    GLuint* p_z_buffer = p_app_buf->get_z_buffer_ptr();
    assert(p_z_buffer != 0);
    glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_UNSIGNED_INT, p_z_buffer);
    
/*
        // dump the read buffers: Note: slow. But only in debug build.
        mi::Sint32 const dump_buffer_per_frame = get_sint32(
            nv::index_common::Index_AppData::instance()->peek_app_proj()->
            get("app::experimental::dump_opengl_z_rgba_buffer_per_frame", "0"));
        if(dump_buffer_per_frame > 0){
            static mi::Sint32 frame_num = 0;
            ++frame_num;
            if((frame_num % dump_buffer_per_frame) == 0){
                std::stringstream sstr;
                sstr << "frame." << frame_num;
                bool ret = p_app_buf->debug_write_buffer_to_file(sstr.str());
                assert(ret);
                INFO_LOG << "frame " << frame_num << " is written with ret = " << ret;
            }
        }
    }
*/
}

//======================================================================

/// convert polygonal cylinder gen object to AppGLTriMesh.
/// \param[in] agt_mesh an AppGLTriMesh
/// \param[in] pc_gen polygonal cylinder generator. The cylinder must
/// have been generated.
static void convert_polygonal_cylinder_gen_to_appgltrimesh(
    AppGLTriMesh & agt_mesh,  Polygonal_cylinder_gen const & pc_gen)
{
    assert(pc_gen.has_cylinder());

    // set generated vertices
    Polygonal_cylinder_gen::Float32_3_vec const & v_vec = pc_gen.get_vertex_vec_ref();
    for(Polygonal_cylinder_gen::Float32_3_vec::const_iterator vi = v_vec.begin();
        vi != v_vec.end(); ++vi)
    {
        agt_mesh.m_vertex_list.push_back(*vi);
    }

    // and set generated mesh topology
    Polygonal_cylinder_gen::Sint32_3_vec const & f_vec = pc_gen.get_side_face_ref();
    for(Polygonal_cylinder_gen::Sint32_3_vec::const_iterator fi = f_vec.begin();
        fi != f_vec.end(); ++fi)
    {
        agt_mesh.m_face_list.push_back(*fi);
    }
}

//----------------------------------------------------------------------
/// generate a simple well. Just a cylinder.
///
/// \param[in] center_begin_point well center begin point (well starting point)
/// \param[in] roi_xyz region of interest. The wells are created in this region.
/// \param[in] center_begin_point well's center point.
/// \param[in] n_gon  n-polygon's n. well cylinder is n-gon shape.
/// \param[in] radius well's radius
/// \param[in] seg_count number of segment of the well.
static void generate_well_simple(
    mi::math::Bbox< mi::Float32, 3 > const & roi_xyz,
    mi::math::Vector< mi::Float32, 3 > const & center_begin_point,
    mi::Sint32  n_gon,
    mi::Float32 radius,
    mi::Sint32  seg_count)
{
    Polygonal_cylinder_gen pc_gen;

    pc_gen.clear();
    pc_gen.set_n_gon(n_gon);
    bool const is_gen_segment_tris = false;
    pc_gen.set_generate_segment_tris(is_gen_segment_tris);

    // not goes to the bottom (one less)
    mi::math::Vector< mi::Float32, 3 > bbox_size = roi_xyz.max - roi_xyz.min;
    mi::Float32 seg_size = bbox_size.z / static_cast< mi::Float32 >(seg_count);
    assert(seg_size > 0.0f);

    // generate simple well
    mi::math::Vector< mi::Float32, 3 > center_point = center_begin_point;
    for(mi::Sint32 i = 0; i < seg_count; ++i){
        pc_gen.append_center_point(center_point, radius);
        center_point.z += seg_size;
    }
    pc_gen.gen_cylinder();

    // create a mesh and add
    AppGLTriMesh well_mesh;
    convert_polygonal_cylinder_gen_to_appgltrimesh(well_mesh, pc_gen);
    AppGLPrimitive::instance()->append_mesh(well_mesh);

    DEBUG_LOG << "generated a simple well at " << center_begin_point;
}

//----------------------------------------------------------------------
/// generate a well with joints.
///
/// \param[in] center_begin_point well center begin point (well starting point)
/// \param[in] roi_xyz region of interest. The wells are created in this region.
/// \param[in] center_begin_point well's center point.
/// \param[in] n_gon  n-polygon's n. well cylinder is n-gon shape.
/// \param[in] radius well's radius
/// \param[in] seg_count number of segment of the well.
/// \param[in] joint_int_count joint interval count
/// \param[in] joint_height_ratio joint relative height agaist to the
/// segment length (must be in (0.0, 1.0))
/// \param[in] joint_radius joint radius
/// \param[in] zig_degree zigzag degree
static void generate_well_joint(
    mi::math::Bbox< mi::Float32, 3 > const & roi_xyz,
    mi::math::Vector< mi::Float32, 3 > const & center_begin_point,
    mi::Sint32  n_gon,
    mi::Float32 radius,
    mi::Sint32  seg_count,
    mi::Sint32  joint_int_count,
    mi::Float32 joint_height_ratio,
    mi::Float32 joint_radius,
    mi::Float32 zig_degree)
{
    Polygonal_cylinder_gen pc_gen;

    pc_gen.clear();
    pc_gen.set_n_gon(n_gon);
    bool const is_gen_segment_tris = false;
    pc_gen.set_generate_segment_tris(is_gen_segment_tris);

    // not goes to the bottom (one less)
    mi::math::Vector< mi::Float32, 3 > bbox_size = roi_xyz.max - roi_xyz.min;
    mi::Float32 const seg_size = bbox_size.z / static_cast< mi::Float32 >(seg_count);
    assert(seg_size > 0.0f);

    if(joint_int_count >= seg_count){
        WARN_LOG << "joint_int_count(" << joint_int_count
                 << ") >= seg_count(" << seg_count << "), no joint will be created.";
    }
    if(joint_height_ratio <= 0.0){
        ERROR_LOG << "joint_height_ratio(" << joint_height_ratio
                  << ") <= 0.0, cannot create joints. Reset to 0.1";
        joint_height_ratio = 0.1f;
    }
    if(joint_height_ratio >= 1.0){
        ERROR_LOG << "joint_height_ratio(" << joint_height_ratio
                  << ") >= 1.0, cannot create joints. Reset to 0.1";
        joint_height_ratio = 0.1f;
    }
    assert(0.0f < joint_height_ratio);
    assert(joint_height_ratio < 1.0f);

    // generate well + joint (zigzag)
    mi::math::Vector< mi::Float32, 3 > center_point = center_begin_point;
    bool is_zig = false;
    for(mi::Sint32 i = 0; i < seg_count; ++i){
        if((i > 0) && ((i % joint_int_count) == 0)){
            // create a joint
            mi::Float32 const joint_shift = seg_size * joint_height_ratio;

            // joint up
            mi::math::Vector< mi::Float32, 3 > joint_up = center_point;
            joint_up.z -= joint_shift;

            // joint down
            mi::math::Vector< mi::Float32, 3 > joint_down = center_point;
            joint_down.z += joint_shift;

            pc_gen.append_center_point(joint_up, radius);
            pc_gen.append_center_point(center_point, joint_radius);
            pc_gen.append_center_point(joint_down, radius);

            is_zig = !is_zig;
        }
        else{
            // normal segment
            pc_gen.append_center_point(center_point, radius);
        }

        center_point.z += seg_size;
        if(is_zig){
            // zig
            center_point.x += zig_degree;
            center_point.y += zig_degree;
        }
        else{
            // zag
            center_point.x -= zig_degree;
            center_point.y -= zig_degree;
        }
    }
    pc_gen.gen_cylinder();

    // create a mesh and add
    AppGLTriMesh well_mesh;
    convert_polygonal_cylinder_gen_to_appgltrimesh(well_mesh, pc_gen);
    AppGLPrimitive::instance()->append_mesh(well_mesh);

    DEBUG_LOG << "generated a joint/zigzag well at " << center_begin_point;
}


//----------------------------------------------------------------------
void gl_application_example_well_creation(
    mi::math::Bbox< mi::Float32, 3 > const & roi_xyz,
    nv::index_common::String_dict const & well_opt)
{
    if(AppGLPrimitive::instance()->has_mesh()){
        // if we have already well, do not create them.
        return;
    }
    INFO_LOG << "creating well objects in bbox " << roi_xyz;

    // get parameters from well_opt
    std::string const well_type_str =
        well_opt.get("app::experimental::example_well::type", "simple");
    mi::Sint32 const n_gon =
        nv::index_common::get_sint32(well_opt.get("app::experimental::example_well::n_gon", "6"));
    mi::Sint32  const seg_count =
        nv::index_common::get_sint32(well_opt.get("app::experimental::example_well::segment_count", "10"));
    mi::Float32 const base_r =
        nv::index_common::get_float32(well_opt.get("app::experimental::example_well::base_radius",  "1.0"));
    mi::Sint32 const i_count =
        nv::index_common::get_sint32(well_opt.get("app::experimental::example_well::i_count", "1"));
    mi::Sint32 const j_count =
        nv::index_common::get_sint32(well_opt.get("app::experimental::example_well::j_count", "1"));

    DEBUG_LOG << "i_count: " << i_count << ", j_count: " << j_count;
    assert(n_gon      >= 3);
    assert(seg_count >= 2);
    assert(i_count   >  0);
    assert(j_count   >  0);

    mi::math::Vector< mi::Float32, 3 > bbox_size = roi_xyz.max - roi_xyz.min;
    mi::Float32 const i_delta = bbox_size.x / static_cast< mi::Float32 >(i_count);
    mi::Float32 const j_delta = bbox_size.y / static_cast< mi::Float32 >(j_count);

    mi::math::Vector< mi::Float32, 3 > center_begin_point = roi_xyz.min;

    if(well_type_str == "simple"){
        for(mi::Sint32 i = 0; i < i_count; ++i){
            center_begin_point.x = roi_xyz.min.x + static_cast< mi::Float32 >(i) * i_delta;
            for(mi::Sint32 j = 0; j < j_count; ++j){
                center_begin_point.y = roi_xyz.min.y + static_cast< mi::Float32 >(j) * j_delta;
                center_begin_point.z = roi_xyz.min.z;
                generate_well_simple(roi_xyz, center_begin_point, n_gon, base_r, seg_count);
            }
        }
    }
    else if(well_type_str == "joint"){
        mi::Sint32 const joint_int_count =
            nv::index_common::get_sint32(well_opt.get("app::experimental::example_well::joint_interval", "7"));
        mi::Float32 const joint_height_ratio =
            nv::index_common::get_float32(well_opt.get("app::experimental::example_well::joint_h_ratio", "0.1"));
        mi::Float32 const joint_radius =
            nv::index_common::get_float32(well_opt.get("app::experimental::example_well::joint_radius", "2.0"));
        mi::Float32 const zig_degree =
            nv::index_common::get_float32(well_opt.get("app::experimental::example_well::zig_degree", "0.0"));

        for(mi::Sint32 i = 0; i < i_count; ++i){
            center_begin_point.x = roi_xyz.min.x + static_cast< mi::Float32 >(i) * i_delta;
            for(mi::Sint32 j = 0; j < j_count; ++j){
                center_begin_point.y = roi_xyz.min.y + static_cast< mi::Float32 >(j) * j_delta;
                center_begin_point.z = roi_xyz.min.z;
                generate_well_joint(roi_xyz, center_begin_point, n_gon, base_r, seg_count,
                                    joint_int_count, joint_height_ratio, joint_radius,
                                    zig_degree);
            }
        }
    }
    else{
        ERROR_LOG << "unknown well type: " << well_type_str << ", no well creation.";
    }
}

//----------------------------------------------------------------------
void gl_application_shutdown()
{
    AppGLPrimitive::delete_instance();
}


//----------------------------------------------------------------------
