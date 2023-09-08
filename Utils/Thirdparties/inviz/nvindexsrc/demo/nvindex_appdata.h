/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief IndeX viewer application data.

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_INDEX_APPDATA_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_INDEX_APPDATA_H

#include <mi/dice.h>
#include <string>

#include <nv/index/iindex.h>
#include <nv/index/iperformance_values.h>

#include "common/string_dict.h"

#include "camera_animator.h"
#include "colormap_manager.h"
#include "host_info.h"
#include "opengl_appdata.h"
#include "performance_logger_app_state.h"
#include "performance_logger_set.h"
#include "user_interaction.h"
#include "utilities.h"


//----------------------------------------------------------------------
// forward declaration
class Nvindex_library_accessor;
//----------------------------------------------------------------------

/// NVIDIA IndeX viewer application data.
///
/// This might remove all the global variables, but, currently we have
/// necessary ones here.
///
/// singleton
class Nvindex_AppData
{
public:
    /// get the instance
    /// \return geostream application data singleton instance
    static Nvindex_AppData * instance()
    {
        if(G_p_nvindex_appdata == 0){
            G_p_nvindex_appdata = new Nvindex_AppData();
        }
        return G_p_nvindex_appdata;
    }

    /// delete the singleton. (For unit test purpose. I don't want to
    /// confuse a memory checker.)
    static void delete_instance()
    {
        if(G_p_nvindex_appdata != 0){
            delete G_p_nvindex_appdata;
            G_p_nvindex_appdata = 0;
        }
    }

private:
    // singleton instance
    static Nvindex_AppData * G_p_nvindex_appdata;

public:
    /// destructor
    virtual ~Nvindex_AppData();

    /// set new-ed library accesor
    /// \param[in] p_nvindex_interface nvindex library access interface
    void set_nvindex_library_accessor(
        Nvindex_library_accessor * p_nvindex_interface);

    /// get current nvindex library accesor
    /// \return nvindex accesor
    Nvindex_library_accessor * get_nvindex_library_accessor() const;

    /// set any (one of) video streaming enabled status
    /// \param[in] is_enabled true when any of the video streaming is enabled
    void set_any_video_stream_enabled(bool is_enabled);
    /// get any of video streaming enabled status
    /// \return true when any of the video streaming is enabled
    bool is_any_video_stream_enabled() const;

    /// Returns whether the bounding boxes of scene elements should be shown.
    bool is_show_bounding_boxes_overlay() const;
    /// Sets whether the bounding boxes of scene elements should be shown.
    void set_show_bounding_boxes_overlay(bool show);

    /// Returns whether the region of interest should be shown.
    bool is_show_region_of_interest_overlay() const;
    /// Sets whether the region of interest should be shown.
    void set_show_region_of_interest_overlay(bool show);

    /// Returns whether horizontal spans should be visualized as rectangles.
    bool is_show_horizontal_spans_overlay() const;
    /// Controls whether horizontal spans should be visualized as rectangles.
    void set_show_horizontal_spans_overlay(bool show);

    /// Returns whether the axes of the coordinate system should be visualized.
    bool is_show_coordinate_system_overlay() const;
    /// Controls whether the axes of the coordinate system should be visualized.
    void set_show_coordinate_system_overlay(bool show);

    /// peek the application option dictionary.
    /// Note: this method is peeking inside the object. Please mind
    /// the dangerous side effect.
    /// \return application dictionary's reference pointer.
    nv::index_common::String_dict * peek_app_proj();

    /// set application run
    /// \param[in] is_run true when running
    void set_app_run(bool is_run);
    /// get run status of application
    /// \return true when application is run
    bool is_app_run() const;

    /// get unique identifier.
    ///
    /// This conunt up the identifier number, so every time you call,
    /// you got an identifier.
    ///
    /// \return unique ID
    mi::Sint32 get_unique_identifier();

    /// get colormap manager
    /// \return reference to the colormap manager
    Colormap_manager* get_colormap_manager();

    /// set current colormap index
    /// \param[in] current_colormap_index current colormap index
    void set_current_colormap_index(mi::Uint32 current_colormap_index);

    /// get current colormap index
    /// \return current_colormap_index
    mi::Uint32 get_current_colormap_index() const;

    /// set rendering availability on volume
    /// \param[in] is_render_volume true when render the volume
    void set_render_volume(bool is_render_volume);

    /// get rendering state of the volume
    /// \return true when rendering the volume
    bool is_render_volume() const;

    /// peek recent n frame statistics buffer
    // Ring_buffer * peek_recent_n_fps_buffer();

    /// peek cumulative statistics
    /// \return reference to the cumulative statistics object. Note:
    /// returning to the member reference
    Cumulative_stat * peek_cumulative_stat();

    /// enabled recent N statistics?
    /// \return true when recent N statistics is on
    bool is_enabled_recent_n_stat() const;
    /// enabled recent N statistics?
    /// \param[in] is_enabled_stat true when recent N statistics is on
    void set_enabled_recent_n_stat(bool is_enabled_stat);
    /// get recent N max fps
    /// \return recent N max fps, when no data -1.
    mi::Float32 get_stat_recent_n_max_fps() const;
    /// set recent N max fps
    /// \param[in] recent_n_max_fps recent N max fps. -1 for no data.
    void set_stat_recent_n_max_fps(mi::Float32 recent_n_max_fps);
    /// get recent N average fps
    /// \return recent N average fps, when no data -1.
    mi::Float32 get_stat_recent_n_ave_fps() const;
    /// set recent N average fps
    /// \param[in] recent_n_ave_fps recent N average fps. -1 for no data.
    void set_stat_recent_n_ave_fps(mi::Float32 recent_n_ave_fps);

    /// Peek the host information (Note: this can change the internal state.)
    Host_info * peek_host_info();

    /// requested canvas resolution
    void set_request_canvas_resolution(const mi::Sint32_2 & requesting_canvas_resolution);
    /// requested canvas resolution
    mi::Sint32_2 get_request_canvas_resolution() const;

    /// set quit after frame number
    void set_quit_after_frame_number(mi::Sint32 quit_after_frame_number);
    /// get quit after frame number
    mi::Sint32 get_quit_after_frame_number() const;

    /// is OpenGL enable application?
    bool is_use_opengl() const;
    /// Get OpenGL application data instance
    OpenGL_AppData* get_opengl_appdata();

public:
    /// \addgroup camera_related_application_state camera_related_application_state
    /// FIXME move to User_Interaction
    /// \brief Unit test for geostream_viewer.
    /// @{

    /// set stero mode
    /// FIXME move to User_Interaction::Multiple_camera
    /// \param[in]  is_on true when on
    void set_stereo_mode(bool is_on);
    /// is stereo mode on?
    /// FIXME move to User_Interaction
    /// \return true when on
    bool is_stereo_mode() const;

    /// set application selected storage camera index
    /// FIXME move to User_Interaction::Multiple_camera
    /// \param[in]  cam_idx selected storage camera index
    void set_storage_camera_index(mi::Sint32 cam_idx);
    /// set application selected storage camera index
    /// FIXME move to User_Interaction
    /// \return selected storage camera index
    mi::Sint32 get_storage_camera_index() const;

    /// set snapshot filename.
    /// The name will get one integer digits with printf
    /// form. (ex. "snap_%03d.ppm")
    /// \param[in]  snapshot_fname snapshot filename
    void set_snapshot_filename(std::string const & snapshot_fname);
    /// set application selected storage camera index
    /// \param[in]  idx snapshot index
    /// \return selected storage camera index
    std::string get_snapshot_filename(mi::Sint32 idx) const;

    /// set snapshot index (index sequence mode)
    /// \param[in] idx shapshot index
    void set_snapshot_index(mi::Sint32 idx);
    /// increment snapshot index
    void inc_snapshot_index();
    /// get snapshot index
    mi::Sint32 get_snapshot_index() const;

    /// set snapshot sequence mode.
    ///   snapshot has two modes.
    ///    - sequence mode (index starts with 0 and continue to 1, 2...)
    ///    - frame number mode (index is equal to rendering frame
    ///      sequence). The frame number should be given by
    ///      application
    ///
    /// \param[in] is_sequence_mode true when sequence mode
    void set_snapshot_sequence_mode(bool is_sequence_mode);
    /// get snapshot sequence mode
    bool is_snapshot_sequence_mode() const;

    // Moves the last snapshot image in place for download, returning the URL to use.
    std::string prepare_last_snapshot() const;

    /// set snapshot on next frame
    /// \param[in] is_shot when true, export snapshot file in the next
    /// chance
    void set_snapshot_on(bool is_shot);
    /// is the snap shot on state?
    /// \return true when the snapshot on
    bool is_snapshot_on() const;

    /// Set animation mode.
    /// In animation mode a project file snippet will be loaded for each frame and the rendering
    /// result stored as a snapshot.
    void set_animation_mode(bool enable);
    /// Return whether animation mdoe is enabled.
    bool is_animation_mode() const;

    /// set image file canvas snapshot on next frame.
    /// The different between snapshot and image file snapshot is
    ///  - snapshot: dump the raw current span buffer
    ///  - image file snapshot: render with specified resolution into
    ///    image file canvas and dump the image file canvas
    /// In short, image file snapshot can change the output resolution
    /// on the fly.
    ///
    /// \param[in] is_shot when true, trigger image file canvas
    /// snapshot file in the next chance
    // void set_image_file_snapshot_on(bool is_shot);
    /// is the image file snap shot on state?
    /// \return true when the snapshot on
    // bool is_image_file_snapshot_on() const;

    /// set camera animator data
    /// \param[in] animator camera animator data
    void set_camera_animator_data(Camera_animator const & animator);
    /// get the camera animator data
    /// \return camera animator data by copy
    Camera_animator get_camera_animator_data() const;
    /// peek the camera animator data
    /// \return camera animator data reference (use with care)
    Camera_animator * peek_camera_animator_data();

    /// @}
    /// camera_related_application_state

    /// \addgroup user interaction
    /// \brief User interaction object access methods
    /// @{

    /// Append an user interaction
    ///
    /// Create and append a user interaction at end of the vector.
    void create_user_interaction();

    /// Get the number of user interaction
    ///
    /// \return get current number of user interaction
    mi::Size get_nb_user_interaction() const;

    /// Get a user interaction
    ///
    /// \param[in] ui_idx user interaction index
    /// \return an user interaction. 0 when error (e.g., invalid index)
    ///
    /// Note: currently only ui_idx == 0 is valid.
    User_interaction* get_user_interaction(mi::Size ui_idx);

    /// @}

    /// Enable/disable pick operation
    /// FIXME move to User_Interaction
    void set_pick_operation_enabled(bool is_enable);
    /// Get enable/disable pick operation
    /// FIXME move to User_Interaction
    bool is_pick_operation_enabled() const;

public:
    /// \addgroup scene_status_related_application_state scene_status_related_application_state
    /// \brief Current scene elements snapshot in the scene
    /// @{

    /// set current volume tags
    void set_volume_tag_vec(const std::vector<mi::neuraylib::Tag> & volume_tag_vec);
    /// get current volume tags
    std::vector<mi::neuraylib::Tag> get_volume_tag_vec() const;

    /// set current sparse volume tags
    void set_sparse_volume_tag_vec(const std::vector<mi::neuraylib::Tag> & svol_tag_vec);
    /// get current sparse volume tags
    std::vector<mi::neuraylib::Tag> get_sparse_volume_tag_vec() const;

    /// get volume name vector from stored volume tags
    ///
    /// \param[in]  scope         dice scope
    /// \param[out] result_volume_name_vec (output) all volume name vector
    void get_volume_name_vec(
        mi::base::Handle<mi::neuraylib::IScope> & scope,
        std::vector< std::string > &              result_volume_name_vec) const;

    /// get steering technique name vector
    ///
    /// \param[out] result_steering_name_vec (output) all steering technique name vector
    void get_steering_names(
        std::vector< std::string > &              result_steering_name_vec) const;

    /// set steering technique name vector
    void set_steering(
        mi::Sint32          method_id,
        const std::string&  method_name);

    /// get volume tag from volume index in the get_volume_tag_vec()'s vector
    ///
    /// \param[in] volume_index volume index.
    /// \return volume tag. NULL_TAG when the volume_index is invalid.
    mi::neuraylib::Tag get_volume_tag_from_idx(mi::Sint32 volume_index) const;

    /// set current heightfield tags
    void set_heightfield_tag_vec(const std::vector<mi::neuraylib::Tag> & heightfield_tag_vec);
    /// get current heightfield tags
    std::vector<mi::neuraylib::Tag> get_heightfield_tag_vec() const;

    /// get heightfield name vector from stored volume tags
    ///
    /// \param[in]  scope        dice scope
    /// \param[out] result_heightfield_name_vec (output) all heightfield name vector
    void get_heightfield_name_vec(
        mi::base::Handle<mi::neuraylib::IScope> & scope,
        std::vector< std::string > &              result_heightfield_name_vec) const;

    /// get heightfield tag from heightfield index in the get_heightfield_tag_vec()'s vector
    ///
    /// \param[in] heightfield_idx heightfield index.
    /// \return heightfield tag. NULL_TAG when the heightfield_index is invalid.
    mi::neuraylib::Tag get_heightfield_tag_from_idx(
        mi::Sint32 heightfield_idx) const;

    /// set current rtc parameter buffer tags
    void set_rtc_param_buffer_tag_vec(const std::vector<mi::neuraylib::Tag> & rtc_buffer_tag_vec);
    /// get current rtc parameter buffer tags
    std::vector<mi::neuraylib::Tag> get_rtc_param_buffer_tag_vec() const;

    /// @}
    /// scene_status_related_application_state

public:
    /// \addgroup application_logging_state application_logging_state 
    /// \brief Application logging out state
    /// @{

    /// Peek application logging state.
    /// peek (instead of get) returns reference to the member, this
    /// means you can change the member state.
    Performance_logger_app_state & peek_log_state();

    /// Peek performance logger
    /// peek (instead of get) returns reference to the member, this
    /// means you can change the member state.
    Performance_logger_set & peek_logger_set();
    
    /// Set performance measuring frames for performance test playback
    /// \param[in] nb_frames 
    void set_remaining_performance_measuring_frames(mi::Sint32 nb_frames);
    /// Get the current performance measuring frames
    /// \return current remaining frames
    mi::Sint32 get_remaining_performance_measuring_frames() const;
    /// decrements the remaining performance measure frames
    void decrement_remaining_performance_measuring_frames();

    /// @}
    /// application_logging_state 

    /// enable/disable multi-view
    void set_enable_multi_view_mode(bool enable);
    /// get multi-view state
    bool is_enable_multi_view_mode() const;
    /// multi-view: set advisory
    void set_enable_multi_view_advisory(bool is_advisory);
    /// multi-view: get advisory
    bool is_enable_multi_view_advisory() const;


public:
    // experimental

    /// experimental feature: Enable OpenGL primitives
    /// \param[in] is_z_enabled_rendering enabled z-buffer rendering when true.
    void set_enable_opengl_z_buffer_for_rendering(bool is_z_enabled_rendering);
    /// experimental feature: Enable OpenGL primitives
    /// \return true when enabled OpenGL Z-buffer rendering
    bool is_enable_opengl_z_buffer_for_rendering() const;
public:
    //TODO: These are the former glboals from geostream_viewer.cpp.
    // Before creating setters and getters for them, it should be checked whether they actually
    // should be stored here, or if there is a more localized place for them.
    bool                                             m_show_statistics_overlay;
    bool                                             m_statistics_overlay_large;
    bool                                             m_show_rendering_statistics_overlay;
    bool                                             m_show_compositing_statistics_overlay;
    bool                                             m_show_all_performance_statistics_overlay;

    mi::math::Color                                  m_background_color;

    /// Render always, even if no video stream is connected
    bool                                             m_render_always;
    /// Use negative value to render initial frames before waiting for connection
    mi::Sint32                                       m_last_nb_video_connections;

    mi::base::Handle<nv::index::IPerformance_values> m_performance_values;
    mi::base::Lock                                   m_performance_values_lock;
    mi::Float32                                      m_application_compute_time;
    mi::Uint64                                       m_application_compute_memory;

    bool                                             m_immediate_final_parallel_compositing;

    mi::base::Lock                                   m_scene_edit_lock;
    mi::base::Lock                                   m_scene_update_lock;

    /// Initial region of interest (as specified in the project file)
    mi::math::Bbox<mi::Float32, 3>                   m_initial_roi;

private:
    /// dice library interface
    Nvindex_library_accessor* m_p_nvindex_interface;

    /// any of the video streaming enabled
    bool m_is_any_video_stream_enabled;

    bool m_show_bounding_boxes_overlay;
    bool m_show_region_of_interest_overlay;
    bool m_show_horizontal_spans_overlay;
    bool m_show_coordinate_system_overlay;

    /// application dict option
    nv::index_common::String_dict m_app_dict;
    /// application running status
    bool m_is_app_run;
    /// for an application unique number
    mi::Sint32 m_unique_id;
    /// lock for m_unique_id.
    mi::base::Lock m_lock_unique_identifier;

    /// Colormap manager instance
    Colormap_manager m_colormap_manager;
    /// current colormap index
    mi::Uint32 m_current_colormap_index;

    /// now rendering volume
    bool m_is_render_volumes;
    // recent N frame statistics buffer
    // Ring_buffer m_recent_n_fps;
    /// cumulative statistics
    Cumulative_stat m_cumulative_fps_stat;

    /// recent N frame statistics on?
    bool m_is_enabled_recent_n_stat;
    /// statistics recent N max fps
    mi::Float32 m_stat_recent_n_max_fps;
    /// statistics recent N average fps
    mi::Float32 m_stat_recent_n_ave_fps;

    /// requested canvas size
    mi::Sint32_2 m_requested_canvas_size;

    /// Quit after number of frame: a user option
    mi::Sint32 m_quit_after_frame_number;

    /// Animation mode
    bool m_animation_mode;

    /// OpenGL application data
    OpenGL_AppData m_opengl_appdata;

    //----------------------------------------------------------------------
    // host information
    //----------------------------------------------------------------------
    Host_info m_host_info;

    //----------------------------------------------------------------------

    /// lock for requesting resolution change
    mutable mi::base::Lock m_request_canvas_resolution_lock;

    //----------------------------------------------------------------------
    // camera related
    //----------------------------------------------------------------------
    /// stereo mode
    bool m_is_stereo_on;
    /// GUI selected camera index
    mi::Sint32 m_selected_storage_camera_idx;
    /// snapshot filename ("snap_%03d.ppm")
    std::string m_snapshot_fname;
    /// snapshot index mode (true when sequence mode)
    bool m_snapshot_sequence_index_mode;
    /// snapshot index (sequence)
    mi::Sint32 m_snapshot_sequence_idx;
    /// snapshot on state
    bool m_is_snapshot_on;
    /// temp for bug 130052
    bool m_temp_is_encoder_snapshot_on;
    /// camera animation data
    Camera_animator m_camera_animator;

    //----------------------------------------------------------------------
    // User interaction
    //----------------------------------------------------------------------
    /// User interaction object vector
    std::vector<User_interaction> m_user_interaction_vec;
    /// pick operation 
    bool                          m_enable_pick_operation;

    //----------------------------------------------------------------------
    // scene status
    //----------------------------------------------------------------------
    /// lock for scene status access
    mutable mi::base::Lock m_scene_hierarchy_access_lock;
    /// current volume tags (checked at each frame)
    std::vector<mi::neuraylib::Tag> m_volume_tag_vec;
    /// current sparse volume tags (checked at each frame)
    std::vector<mi::neuraylib::Tag> m_sparse_volume_tag_vec;
    /// current heightfield tags (checked at each frame)
    std::vector<mi::neuraylib::Tag> m_heightfield_tag_vec;
    /// current rtc parameter buffer tags (checked at each frame)
    std::vector<mi::neuraylib::Tag> m_rtc_parameter_buffer_tag_vec;

    //----------------------------------------------------------------------
    // application logging state
    //----------------------------------------------------------------------
    // application state for logging 
    Performance_logger_app_state m_log_state;
    // performance logger set
    Performance_logger_set m_perf_logger_set;
    // performance measuring frames
    mi::Sint32 m_remaining_performance_measuring_frames;

    //----------------------------------------------------------------------
    // experimental features
    //----------------------------------------------------------------------
    /// OpenGL primitive Z buffer rendering
    bool m_is_z_enabled_rendering;
    /// multi-view mode
    bool m_is_multi_view_mode_on;
    /// is multi-view advisory mode
    bool m_is_multi_view_advisory;
    
private:
    /// default constructor
    Nvindex_AppData();

private:
    /// copy constructor. prohibit until proved useful.
    Nvindex_AppData(Nvindex_AppData const &);
    /// operator=. prohibit until proved useful.
    Nvindex_AppData const & operator=(Nvindex_AppData const &);
};

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_INDEX_APPDATA_H
