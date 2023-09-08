/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "nvindex_appdata.h"
#include <nv/index/iregular_volume.h>
#include <nv/index/iregular_heightfield.h>

#include <cassert>
#include <fstream>

#include "nvindex_library_accessor.h"

#include "common/forwarding_logger.h"
#include "../alya/alyaglobal.h"
/// singleton instance
Nvindex_AppData * Nvindex_AppData::G_p_nvindex_appdata = 0;

//----------------------------------------------------------------------
Nvindex_AppData::~Nvindex_AppData()
{
    if(m_p_nvindex_interface != 0)
    {
        delete m_p_nvindex_interface;
        m_p_nvindex_interface = 0;
    }

    m_user_interaction_vec.clear();
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_nvindex_library_accessor(
    Nvindex_library_accessor * p_nvindex_interface)
{
    assert(p_nvindex_interface != 0);
    assert(m_p_nvindex_interface == 0); // check double initialization

    m_p_nvindex_interface = p_nvindex_interface;
}

//----------------------------------------------------------------------
Nvindex_library_accessor * Nvindex_AppData::get_nvindex_library_accessor() const
{
    return m_p_nvindex_interface;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_any_video_stream_enabled(bool is_enabled)
{
    m_is_any_video_stream_enabled = is_enabled;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_any_video_stream_enabled() const
{
    return m_is_any_video_stream_enabled;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_show_bounding_boxes_overlay(bool show)
{
    m_show_bounding_boxes_overlay = show;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_show_bounding_boxes_overlay() const
{
    return m_show_bounding_boxes_overlay;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_show_region_of_interest_overlay(bool show)
{
    m_show_region_of_interest_overlay = show;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_show_region_of_interest_overlay() const
{
    return m_show_region_of_interest_overlay;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_show_horizontal_spans_overlay(bool show)
{
    m_show_horizontal_spans_overlay = show;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_show_horizontal_spans_overlay() const
{
    return m_show_horizontal_spans_overlay;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_show_coordinate_system_overlay(bool show)
{
    m_show_coordinate_system_overlay = show;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_show_coordinate_system_overlay() const
{
    return m_show_coordinate_system_overlay;
}

//----------------------------------------------------------------------
nv::index_common::String_dict * Nvindex_AppData::peek_app_proj()
{
    return &m_app_dict;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_app_run(bool is_run)
{
    m_is_app_run = is_run;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_app_run() const
{
    return m_is_app_run;
}

//----------------------------------------------------------------------
mi::Sint32 Nvindex_AppData::get_unique_identifier()
{
    mi::base::Lock::Block block(&(this->m_lock_unique_identifier));
    mi::Sint32 const retval = m_unique_id;
    {
        ++m_unique_id;          // critical section
    }
    return retval;
}

//----------------------------------------------------------------------
Colormap_manager* Nvindex_AppData::get_colormap_manager()
{
    return &m_colormap_manager;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_current_colormap_index(mi::Uint32 current_colormap_index)
{
    m_current_colormap_index = current_colormap_index;
}

//----------------------------------------------------------------------
mi::Uint32 Nvindex_AppData::get_current_colormap_index() const
{
    return m_current_colormap_index;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_render_volume(bool is_render_volume)
{
    m_is_render_volumes = is_render_volume;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_render_volume() const
{
    return m_is_render_volumes;
}

//----------------------------------------------------------------------
Cumulative_stat * Nvindex_AppData::peek_cumulative_stat()
{
    return &m_cumulative_fps_stat;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_enabled_recent_n_stat() const
{
    return m_is_enabled_recent_n_stat;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_enabled_recent_n_stat(bool is_enabled_stat)
{
    m_is_enabled_recent_n_stat = is_enabled_stat;
}

//----------------------------------------------------------------------
mi::Float32 Nvindex_AppData::get_stat_recent_n_max_fps() const
{
    return m_stat_recent_n_max_fps;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_stat_recent_n_max_fps(mi::Float32 recent_n_max_fps)
{
    m_stat_recent_n_max_fps = recent_n_max_fps;
}

//----------------------------------------------------------------------
mi::Float32 Nvindex_AppData::get_stat_recent_n_ave_fps() const
{
    return m_stat_recent_n_ave_fps;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_stat_recent_n_ave_fps(mi::Float32 recent_n_ave_fps)
{
    m_stat_recent_n_ave_fps = recent_n_ave_fps;
}

//----------------------------------------------------------------------
Host_info * Nvindex_AppData::peek_host_info()
{
    return &m_host_info;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_request_canvas_resolution(const mi::Sint32_2 & requesting_canvas_resolution)
{
    if((requesting_canvas_resolution.x <= 0) || (requesting_canvas_resolution.y <= 0)){
        ERROR_LOG << "Can't request negative canvas size: " << requesting_canvas_resolution << ", ignored.";
        return;
    }
    mi::base::Lock::Block block(&m_request_canvas_resolution_lock);    
    m_requested_canvas_size = requesting_canvas_resolution;
}

//----------------------------------------------------------------------
mi::Sint32_2 Nvindex_AppData::get_request_canvas_resolution() const
{
    mi::base::Lock::Block block(&m_request_canvas_resolution_lock);    
    return m_requested_canvas_size;    
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_quit_after_frame_number(mi::Sint32 quit_after_frame_number)
{
    m_quit_after_frame_number = quit_after_frame_number;
}

//----------------------------------------------------------------------
mi::Sint32 Nvindex_AppData::get_quit_after_frame_number() const
{
    return m_quit_after_frame_number;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_use_opengl() const
{
#ifdef USE_OPENGL
    return true;
#else
    return false;
#endif
}

//----------------------------------------------------------------------
OpenGL_AppData* Nvindex_AppData::get_opengl_appdata()
{
    return &m_opengl_appdata;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_stereo_mode(bool is_on)
{
    m_is_stereo_on = is_on;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_stereo_mode() const
{
    return m_is_stereo_on;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_storage_camera_index(mi::Sint32 cam_idx)
{
    m_selected_storage_camera_idx = cam_idx;
}

//----------------------------------------------------------------------
mi::Sint32 Nvindex_AppData::get_storage_camera_index() const
{
    return m_selected_storage_camera_idx;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_snapshot_filename(std::string const & snapshot_fname)
{
    if(snapshot_fname.empty()){
        ERROR_LOG << "snapshot filename is empty. ignored.";
        return;
    }
    m_snapshot_fname = snapshot_fname;
}

//----------------------------------------------------------------------
std::string Nvindex_AppData::get_snapshot_filename(mi::Sint32 idx) const
{
    mi::Sint32 const BUFSIZE = 1024;
    char buf[BUFSIZE];
    snprintf(buf, BUFSIZE, m_snapshot_fname.c_str(), idx);
    std::string const fname(buf);

    return fname;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_snapshot_index(mi::Sint32 idx)
{
    m_snapshot_sequence_idx = idx;
}

//----------------------------------------------------------------------
void Nvindex_AppData::inc_snapshot_index()
{
    ++m_snapshot_sequence_idx;
}

//----------------------------------------------------------------------
mi::Sint32 Nvindex_AppData::get_snapshot_index() const
{
    return m_snapshot_sequence_idx;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_snapshot_sequence_mode(bool is_sequence_mode)
{
    m_snapshot_sequence_index_mode = is_sequence_mode;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_snapshot_sequence_mode() const
{
    return m_snapshot_sequence_index_mode;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_snapshot_on(bool is_shot)
{
    m_is_snapshot_on = is_shot;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_snapshot_on() const
{
    return m_is_snapshot_on;
}

//----------------------------------------------------------------------
std::string Nvindex_AppData::prepare_last_snapshot() const
{
    mi::Sint32 sidx = get_snapshot_index();

    if (sidx <= 0)
    {
        ERROR_LOG << "No snapshot file available yet";
        return "";
    }

    std::string last = get_snapshot_filename(sidx - 1);
    std::string dest;

#ifdef LINUX
    // Try to convert the PPM to PNG using the 'convert' command
    dest = "index_snapshot.png";

    // None of the filenames are user-defined, so this call should be safe
    if (std::system(("convert " + last + " " + dest + " 2> /dev/null").c_str()) == 0)
        return dest;
#endif

    // Otherwise just copy the PPM file and return it
    dest = "index_snapshot.ppm";
    std::ifstream ifs(last.c_str(), std::ios::binary);
    if (ifs.fail())
    {
        ERROR_LOG << "Could not open file '" << last << "' for reading.";
        return "";
    }

    std::ofstream ofs(dest.c_str(), std::ios::binary | std::ios::trunc);
    if (ofs.fail())
    {
        ERROR_LOG << "Could not open file '" << dest << "' for writing.";
        return "";
    }

    ofs << ifs.rdbuf();

    return dest;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_animation_mode(bool enable)
{
    m_animation_mode = enable;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_animation_mode() const
{
    return m_animation_mode;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_camera_animator_data(Camera_animator const & animator)
{
    m_camera_animator = animator;
}

//----------------------------------------------------------------------
Camera_animator Nvindex_AppData::get_camera_animator_data() const
{
    return m_camera_animator;
}

//----------------------------------------------------------------------
Camera_animator * Nvindex_AppData::peek_camera_animator_data()
{
    return &m_camera_animator;
}

//----------------------------------------------------------------------
void Nvindex_AppData::create_user_interaction()
{
    m_user_interaction_vec.push_back(User_interaction());    
}

//----------------------------------------------------------------------
mi::Size Nvindex_AppData::get_nb_user_interaction() const
{
    return m_user_interaction_vec.size();
}

//----------------------------------------------------------------------
User_interaction* Nvindex_AppData::get_user_interaction(mi::Size ui_idx)
{
    if (ui_idx >= m_user_interaction_vec.size())
    {
        ERROR_LOG << "invalid ui_idx: " << ui_idx;
        return 0;
    }

    return &(m_user_interaction_vec[ui_idx]);
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_pick_operation_enabled(bool is_enable)
{
    m_enable_pick_operation = is_enable;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_pick_operation_enabled() const
{
    return m_enable_pick_operation;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_volume_tag_vec(const std::vector<mi::neuraylib::Tag> & volume_tag_vec)
{
    mi::base::Lock::Block block(&m_scene_hierarchy_access_lock);    
    m_volume_tag_vec = volume_tag_vec;
}

//----------------------------------------------------------------------
std::vector<mi::neuraylib::Tag> Nvindex_AppData::get_volume_tag_vec() const
{
    std::vector<mi::neuraylib::Tag> tmp_tag_vec;
    {
        mi::base::Lock::Block block(&m_scene_hierarchy_access_lock);
        tmp_tag_vec = m_volume_tag_vec;
    }
    return tmp_tag_vec;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_sparse_volume_tag_vec(const std::vector<mi::neuraylib::Tag> & svol_tag_vec)
{
    mi::base::Lock::Block block(&m_scene_hierarchy_access_lock);    
    m_sparse_volume_tag_vec = svol_tag_vec;
}

//----------------------------------------------------------------------
std::vector<mi::neuraylib::Tag> Nvindex_AppData::get_sparse_volume_tag_vec() const
{
    std::vector<mi::neuraylib::Tag> tmp_tag_vec;
    {
        mi::base::Lock::Block block(&m_scene_hierarchy_access_lock);
        tmp_tag_vec = m_sparse_volume_tag_vec;
    }
    return tmp_tag_vec;
}

//----------------------------------------------------------------------
void Nvindex_AppData::get_steering_names(
    std::vector< std::string > &              result_name_vec) const
{
    result_name_vec.clear();
    result_name_vec.push_back("Depolarize/Reset");
    result_name_vec.push_back("Centers 1");
    result_name_vec.push_back("Centers 2");
    result_name_vec.push_back("Centers 3");
    result_name_vec.push_back("Centers 4");
}

/// set steering technique name vector
void Nvindex_AppData::set_steering(
    mi::Sint32          method_id,
    const std::string&  method_name)
{
  WARN_LOG << "Please invoke your specific steering technique (id: "
	   << method_id << ", name: " << method_name << ") and apply it to your heart/AYLA simulation from here: \n"
	   << "(file " << __FILE__ << ", line: " << __LINE__ << ")";
}

//----------------------------------------------------------------------
void Nvindex_AppData::get_volume_name_vec(
    mi::base::Handle<mi::neuraylib::IScope> & scope,
    std::vector< std::string > &              result_volume_name_vec) const
{
    assert(scope.is_valid_interface());

    result_volume_name_vec.clear();

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        const std::vector<mi::neuraylib::Tag> volume_tag_vec = this->get_volume_tag_vec();
        const mi::Uint32 nb_volumes = volume_tag_vec.size();

        for(mi::Uint32 i=0; i<nb_volumes; ++i){
            mi::base::Handle<const nv::index::IRegular_volume> volume_scene_element(
                dice_transaction->access<const nv::index::IRegular_volume>(volume_tag_vec.at(i)));
            assert(volume_scene_element.is_valid_interface());
        
            result_volume_name_vec.push_back(volume_scene_element->get_name());
        }
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
mi::neuraylib::Tag Nvindex_AppData::get_volume_tag_from_idx(
    mi::Sint32 volume_index) const
{
    // invalid volume index.
    if(volume_index < 0){
        WARN_LOG << "get_volume_tag_from_idx: specified seicmic volume index ("
                 << volume_index << ") is out of range. return NULL tag.";
        return mi::neuraylib::NULL_TAG;
    }

    mi::Uint32 ui_volume_index = static_cast< mi::Uint32 >(volume_index);
    mi::neuraylib::Tag volume_tag = mi::neuraylib::NULL_TAG;

    const std::vector<mi::neuraylib::Tag> volume_tag_vec = this->get_volume_tag_vec();
    if(ui_volume_index < volume_tag_vec.size()){
        volume_tag = volume_tag_vec.at(ui_volume_index);
    }
    else{
        WARN_LOG << "get_volume_tag_from_idx: specified volume index ("
                 << ui_volume_index << ") is out of range. return NULL tag.";
    }

    return volume_tag;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_heightfield_tag_vec(const std::vector<mi::neuraylib::Tag> & heightfield_tag_vec)
{
    mi::base::Lock::Block block(&m_scene_hierarchy_access_lock);    
    m_heightfield_tag_vec = heightfield_tag_vec;
}

//----------------------------------------------------------------------
std::vector<mi::neuraylib::Tag> Nvindex_AppData::get_heightfield_tag_vec() const
{
    std::vector<mi::neuraylib::Tag> tmp_tag_vec;
    {
        mi::base::Lock::Block block(&m_scene_hierarchy_access_lock);    
        tmp_tag_vec = m_heightfield_tag_vec;
    }
    return tmp_tag_vec;
}

//----------------------------------------------------------------------
void Nvindex_AppData::get_heightfield_name_vec(
    mi::base::Handle<mi::neuraylib::IScope> & scope,
    std::vector< std::string > &              result_heightfield_name_vec) const
{
    assert(scope.is_valid_interface());

    result_heightfield_name_vec.clear();

    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        const std::vector<mi::neuraylib::Tag> hf_tag_vec = this->get_heightfield_tag_vec();
        const mi::Uint32 nb_heightfield = hf_tag_vec.size();
        for(mi::Uint32 i=0; i<nb_heightfield; ++i)
        {
            mi::base::Handle<const nv::index::IRegular_heightfield> heightfield_scene_element(
                dice_transaction->access<const nv::index::IRegular_heightfield>(hf_tag_vec.at(i)));
            assert(heightfield_scene_element.is_valid_interface());
            assert(heightfield_scene_element->get_name() != 0);

            result_heightfield_name_vec.push_back(std::string(heightfield_scene_element->get_name()));
        }
    }
    dice_transaction->commit();
}

//----------------------------------------------------------------------
mi::neuraylib::Tag Nvindex_AppData::get_heightfield_tag_from_idx(
    mi::Sint32 heightfield_idx) const
{
    // invalid heightfield index.
    if(heightfield_idx < 0){
        WARN_LOG << "get_heightfield_tag_from_idx: specified heightfield index ("
                 << heightfield_idx << ") is out of range. return NULL tag.";
        return mi::neuraylib::NULL_TAG;
    }

    mi::Uint32 ui_heightfield_idx = static_cast< mi::Uint32 >(heightfield_idx);
    mi::neuraylib::Tag heightfield_tag = mi::neuraylib::NULL_TAG;

    const std::vector<mi::neuraylib::Tag> hf_tag_vec = this->get_heightfield_tag_vec();
    if(ui_heightfield_idx < hf_tag_vec.size()){
        heightfield_tag = hf_tag_vec.at(ui_heightfield_idx);
    }
    else{
        WARN_LOG << "get_heightfield_tag_from_idx: specified heightfield index ("
                 << ui_heightfield_idx << ") is out of range. return NULL tag.";
    }

    return heightfield_tag;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_rtc_param_buffer_tag_vec(const std::vector<mi::neuraylib::Tag> & rtc_buffer_tag_vec)
{
    mi::base::Lock::Block block(&m_scene_hierarchy_access_lock);    
    m_rtc_parameter_buffer_tag_vec = rtc_buffer_tag_vec;
}

//----------------------------------------------------------------------
std::vector<mi::neuraylib::Tag> Nvindex_AppData::get_rtc_param_buffer_tag_vec() const
{
    std::vector<mi::neuraylib::Tag> tmp_tag_vec;
    {
        mi::base::Lock::Block block(&m_scene_hierarchy_access_lock);    
        tmp_tag_vec = m_rtc_parameter_buffer_tag_vec;
    }
    return tmp_tag_vec;
}
//----------------------------------------------------------------------
Performance_logger_app_state & Nvindex_AppData::peek_log_state()
{
    return m_log_state;
}

//----------------------------------------------------------------------
Performance_logger_set & Nvindex_AppData::peek_logger_set()
{
    return m_perf_logger_set;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_remaining_performance_measuring_frames(mi::Sint32 nb_frames)
{
    m_remaining_performance_measuring_frames = nb_frames;
}

//----------------------------------------------------------------------
mi::Sint32 Nvindex_AppData::get_remaining_performance_measuring_frames() const
{
    return m_remaining_performance_measuring_frames;
}

//----------------------------------------------------------------------
void Nvindex_AppData::decrement_remaining_performance_measuring_frames()
{
    if(m_remaining_performance_measuring_frames > 0){
        --m_remaining_performance_measuring_frames;
    }
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_enable_opengl_z_buffer_for_rendering(bool is_z_enabled_rendering)
{
    m_is_z_enabled_rendering = is_z_enabled_rendering;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_enable_opengl_z_buffer_for_rendering() const
{
    return m_is_z_enabled_rendering;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_enable_multi_view_mode(bool enable)
{
    m_is_multi_view_mode_on = enable;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_enable_multi_view_mode() const
{
    return m_is_multi_view_mode_on;
}

//----------------------------------------------------------------------
void Nvindex_AppData::set_enable_multi_view_advisory(bool is_advisory)
{
    m_is_multi_view_advisory = is_advisory;
}

//----------------------------------------------------------------------
bool Nvindex_AppData::is_enable_multi_view_advisory() const
{
    return m_is_multi_view_advisory;
}

//----------------------------------------------------------------------
Nvindex_AppData::Nvindex_AppData()
    :
    m_show_statistics_overlay(false),
    m_statistics_overlay_large(false),
    m_show_rendering_statistics_overlay(true),
    m_show_compositing_statistics_overlay(true),
    m_show_all_performance_statistics_overlay(true),
    m_render_always(true),
    m_last_nb_video_connections(-10),
    m_application_compute_time(0.f),
    m_application_compute_memory(0),
    m_immediate_final_parallel_compositing(true),
    m_p_nvindex_interface(0),
    m_is_any_video_stream_enabled(false),
    m_show_bounding_boxes_overlay(false),
    m_show_region_of_interest_overlay(false),
    m_show_horizontal_spans_overlay(false),
    m_show_coordinate_system_overlay(false),
    m_is_app_run(true),
    m_unique_id(0),
    m_current_colormap_index(0),
    m_is_render_volumes(true),
    m_is_enabled_recent_n_stat(false),
    m_stat_recent_n_max_fps(-1.0f),
    m_stat_recent_n_ave_fps(-1.0f),
    m_requested_canvas_size(1024, 1024),
    m_quit_after_frame_number(-1),
    m_animation_mode(false),
    m_is_stereo_on(false),
    m_selected_storage_camera_idx(0),
    m_snapshot_fname("snap_%03d.ppm"),
    m_snapshot_sequence_index_mode(true),
    m_snapshot_sequence_idx(0),
    m_is_snapshot_on(false),
    m_temp_is_encoder_snapshot_on(false),
    m_camera_animator(),
    m_enable_pick_operation(true),
    m_remaining_performance_measuring_frames(0),
    m_is_z_enabled_rendering(false),
    m_is_multi_view_mode_on(false),
    m_is_multi_view_advisory(false)
{
    // empty
}

//----------------------------------------------------------------------
