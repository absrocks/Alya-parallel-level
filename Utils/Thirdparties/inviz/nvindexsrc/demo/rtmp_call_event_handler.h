/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTMP_CALL_EVENT_HANDLER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTMP_CALL_EVENT_HANDLER_H

#include <mi/dice.h>

#include "common/common_utility.h"
#include "common/forwarding_logger.h"
#include "common/string_dict.h"
#include "common/type_conversion_utility.h"


#include "colormap_manager.h"
#include "colormap_util.h"
#include "command_handler_base.h"
#include "compute_utility.h"
#include "heightfield_data_export.h"
#include "multiple_camera.h"
#include "nvindex_appdata.h"
#include "nvindex_rendering_context.h"
#include "scene_utility.h"
#include "utilities.h"
#include "volume_data_export.h"
#include "volume_data_filter.h"

// An RTMP call event handler
class Call_event_handler : public mi::base::Interface_implement<mi::rtmp::ICall_event_handler> //-V553 PVS
{
protected:
    mi::Uint32                                        m_viewing_scenario;
    Nvindex_rendering_context&                        m_irc;
    mi::base::Handle<nv::index::IIndex>               m_iindex_if;
    mi::base::Handle<nv::index::IIndex_session>       m_iindex_session;
    mi::base::Handle<nv::index::IIndex_rendering>     m_iindex_rendering;
    mi::neuraylib::Tag                                m_session_tag;
    mi::base::Handle<nv::index::IProgress_callback>   m_progress_callback;
    mi::base::Handle<Frame_info_callbacks>            m_frame_info_callbacks;
    // GUI state for event handler method communication (some set
    // values and the other read it thorough this object.)
    nv::index_common::String_dict                     m_gui_state_dict;

    mi::base::Handle<nv::index::IIndex_scene_query>   m_iindex_scene_query;

    std::vector<ICall_event_command_handler*>         m_command_handlers;

    /// handle RTMP call arguments: get
    struct Arguments_get_wrapper
    {
        Arguments_get_wrapper(mi::base::Handle<const mi::IMap> imap)
          : m_imap(imap)
        {
        }

        const char* get_value(const char* key)
        {
            assert(key != 0);
            if (m_string_map.find(key) != m_string_map.end())
            {
                return m_string_map[key].c_str();
            }

            mi::base::Handle<const mi::IString> tmp(m_imap->get_value<mi::IString>(key));
            if (tmp == 0)
            {
                return 0;
            }

            m_string_map[key] = tmp->get_c_str();
            return m_string_map[key].c_str();
        }

        /// get the contents as a string
        std::string to_string()
        {
            std::stringstream sstr;
            for(std::map< std::string, std::string >::const_iterator si = m_string_map.begin();
                si != m_string_map.end(); ++si)
            {
                sstr << si->first << ": " << si->second << "\n";
            }
            sstr << "std::map_size: " << m_string_map.size() << "\n";

            const mi::Size len = m_imap->get_length();
            for(mi::Size i = 0; i < len; ++i){
                sstr << "IMap: " << m_imap->get_key(i) << ": " << this ->get_value(m_imap->get_key(i)) << "\n";
            }
            sstr << "IMap_size: " << len;                 

            return sstr.str();
        }

        std::map< std::string, std::string > m_string_map;            
        mi::base::Handle<const mi::IMap>     m_imap;
    };

    /// handle RTMP call arguments: set
    struct Arguments_set_wrapper
    {
        /// argument set wrapper
        ///
        /// \param[in] response_args response argument
        /// \param[in] iindex_if nvindex library interface
        Arguments_set_wrapper(
            mi::IData** response_args,
            mi::base::Handle<nv::index::IIndex>& iindex_if)
          : m_iindex_if(iindex_if)
        {
            assert(m_iindex_if.is_valid_interface());
            mi::base::Handle<mi::neuraylib::IFactory> factory(
                m_iindex_if->get_api_component<mi::neuraylib::IFactory>());
            if (factory.is_valid_interface())
            {
                m_imap = factory->create<mi::IMap>("Map<Interface>");
                if (m_imap)
                    *response_args = m_imap;
            }
            else
            {
                ERROR_LOG << "fail to create Arguments_set_wrapper.";
            }
        }

        /// set value
        ///
        /// \param[in] key   the key
        /// \param[in] value the value to set
        mi::Sint32 set_value(const char* key, const char* value)
        {
            assert(m_iindex_if.is_valid_interface() != 0);
            mi::base::Handle<mi::neuraylib::IFactory> factory(m_iindex_if->get_api_component<mi::neuraylib::IFactory>());

            if (!factory.is_valid_interface())
            {
                ERROR_LOG << "factory is invalid in Arguments_set_wrapper.";
                return -1;
            }

            if (m_imap && key && value)
            {
                mi::base::Handle<mi::IString> str(factory->create<mi::IString>());
                if (!str.is_valid_interface())
                {
                    return -2;
                }
                str->set_c_str(value);
                if (m_imap->has_key(key))
                {
                    m_imap->set_value(key, str.get());
                }
                else
                {
                    m_imap->insert(key, str.get());
                }
            }

            return 0;
        }

        mi::IMap*                           m_imap;
        mi::base::Handle<nv::index::IIndex> m_iindex_if;
    };

    mi::base::Handle<mi::neuraylib::IScope> get_scope() const
    {
        return m_irc.get_current_scope();
    }

public:
    /// constructor
    Call_event_handler(
        mi::Uint32                                       viewing_scenario,
        Nvindex_rendering_context&                       irc);

    /// destructor
    virtual ~Call_event_handler();

    /// Handle the remote command call
    virtual bool handle( //-V553 PVS
        mi::rtmp::IConnection*      /*connection*/,
        const char*                 /*procedure_name*/,
        const mi::IData*            /*command_args*/,
        const mi::IData*            user_args,
        mi::IData**                 response_args);

    /// Play back the recorded commands from a text file, starting at the given line
    mi::Sint32 play_recording(const std::string& filename, mi::Sint32 start_line);

    /// For demonstrating multi-view
    void render_viewing_scenario(
        mi::Uint32                                                id,
        const mi::base::Handle<mi::neuraylib::IDice_transaction>& dice_transaction);

private:
    //======================================================================
    // utility to accessing arguments
    //======================================================================
    /// get std::string value from userarguments
    ///
    /// \param[in]  user_arguments user argument wrapper
    /// \param[in]  arg_key user argument key
    /// \return string value corresponging to arg_key in arguments
    static std::string get_string_userarg(Arguments_get_wrapper& user_arguments,
                                          const std::string& arg_key)
    {
        assert(!arg_key.empty());

        const char * p_str = user_arguments.get_value(arg_key.c_str());
        if (p_str == 0)
        {
            WARN_LOG << "Could not find the [" << arg_key 
                     << "]. Please check the project/record files. Return empty std::string.";
            return std::string();
        }
        return std::string(p_str);
    }

    /// get mi::Float32 value from userarguments
    ///
    /// \param[in]  user_arguments user argument wrapper
    /// \param[in]  arg_key user argument key
    /// \return mi::Float32 value corresponging to arg_key in arguments
    static mi::Float32 get_float32_userarg(Arguments_get_wrapper& user_arguments,
                                           const std::string& arg_key)
    {
        assert(!arg_key.empty());

        const char * p_float_str = user_arguments.get_value(arg_key.c_str());
        if (p_float_str == 0)
        {
            WARN_LOG << "Could not find the [" << arg_key 
                     << "]. Please check the project/record files. Return Float32 0.0.";
            return 0.0;
        }
        return nv::index_common::get_float32(p_float_str);
    }

    /// get Sint32 value from userarguments
    ///
    /// \param[in]  user_arguments user argument wrapper
    /// \param[in]  arg_key user argument key
    /// \return Sint32 value corresponging to arg_key in arguments
    static mi::Sint32 get_sint32_userarg(Arguments_get_wrapper& user_arguments,
                                         const std::string& arg_key)
    {
        assert(!arg_key.empty());

        const char * p_sint_str = user_arguments.get_value(arg_key.c_str());
        if (p_sint_str == 0)
        {
            WARN_LOG << "Could not find the [" << arg_key 
                     << "]. Please check the project/record files. Return Sint32 0.";
            return 0;
        }
        return nv::index_common::get_sint32(p_sint_str);
    }

    /// get Uint32 value from userarguments
    ///
    /// \param[in]  user_arguments user argument wrapper
    /// \param[in]  arg_key user argument key
    /// \return Uint32 value corresponging to arg_key in arguments
    static mi::Uint32 get_uint32_userarg(Arguments_get_wrapper& user_arguments,
                                         const std::string& arg_key)
    {
        assert(!arg_key.empty());
        
        const char * p_uint_str = user_arguments.get_value(arg_key.c_str());
        if (p_uint_str == 0)
        {
            WARN_LOG << "Could not find the [" << arg_key 
                     << "]. Please check the project/record files. Return Uint32 0.";
            return 0;
        }
        return nv::index_common::get_uint32(p_uint_str);
    }

    /// get bool ("on" or "off") value from userarguments
    ///
    /// \param[in]  user_arguments user argument wrapper
    /// \param[in]  arg_key user argument key
    /// \return int value corresponging to arg_key in arguments
    static bool get_bool_userarg(Arguments_get_wrapper& user_arguments,
                                 const std::string& arg_key)
    {
        assert(!arg_key.empty());

        const char * p_bool_str = user_arguments.get_value(arg_key.c_str());
        if (p_bool_str == 0)
        {
            WARN_LOG << "Could not find the [" << arg_key 
                     << "]. Please check the project/record files. Return bool false.";
            return false;
        }
        const std::string bool_str(p_bool_str);
        if (bool_str == "on")
        {
            return true;
        }
        else if (bool_str == "off")
        {
            return false;            
        }

        WARN_LOG << "Unexpected bool userarg [" << bool_str << "], return bool false;" ;
        return false;
    }

    //======================================================================
    // process a user command (command line command)
    //======================================================================

    /// user command: edit_volume [volume_tag] [amplitude_value]
    /// See the com_usage_msg in the method for the details.
    ///
    /// \param[in] com command string from the user
    /// \return result message
    std::string user_command_edit_volume(const std::string& com);

    /// user command: copy_volume [src_volume_tag] [dst_group_tag] [x] [y] [z]
    ///
    /// - src_volume_tag: copy source volume tag. Must be a volume
    /// - dst_group_tag:  copy destination group tag. Must be a Static_group_node.
    /// - x y z:          translation (x, y, z)
    ///
    /// \param[in] com command string from the user
    /// \return result message
    std::string user_command_copy_volume(const std::string& com);

    /// user command: filter_volume <volume_index>
    ///
    /// \param[in] com command string from the user
    /// \return result message
    std::string user_command_filter_volume(const std::string& com);

    /// user command: distributed_job <volume_index>
    /// distributed job example template test.
    /// \param[in] com command string from the user
    /// \return result message
    std::string user_command_template_distributed_job(const std::string& com);

    /// user command: save_colormap
    ///
    /// \param[in] com command string from the user
    /// \return result message
    std::string user_command_save_colormap(const std::string& com);

    //======================================================================
    // camera related commands
    //======================================================================
    /// predefined view parameter set
    void camera_predefined_view_set_option(Arguments_get_wrapper& user_args);
    /// predefined view parameter get
    void camera_predefined_view_get_option(Arguments_set_wrapper& response_arguments);
    /// predefined view action
    void camera_predefined_view_action();

    /// icamera camera set option to the application
    void camera_set_option(Arguments_get_wrapper& user_args);
    /// icamera camera option. GUI gets options
    void camera_get_option(Arguments_set_wrapper& response_arguments);
    /// icamera set camera parameters. Used by recording. 
    void camera_set_parameter(Arguments_get_wrapper& args);
    /// icamera set current view to camera x
    void camera_set_view_action();
    /// icamera restore camera x to current view
    void camera_restore_view_action();
    /// icamera print parameter action
    void camera_print_param_action();

    //======================================================================
    // colormap handling: user defined colormap implementation
    //======================================================================
    /// get current user colormap. all the colormap values
    /// \param[out] color_table current colormap table
    void get_user_current_colormap_value(std::vector<mi::math::Color_struct>& color_table);
    /// get current user colormap
    /// \param[out] response_arguments return value arguments to the flash player
    void get_user_current_colormap(Arguments_set_wrapper& response_arguments);
    /// get current user colormap all
    /// \param[out] response_arguments return value arguments to the flash player
    void get_user_current_colormap_all(Arguments_set_wrapper& response_arguments);
    /// set user colormap
    /// \param[out] user_arguments user arguments from the flash player
    void set_user_colormap(Arguments_get_wrapper& user_arguments);

    //======================================================================
    // export heightfield data related methods
    //======================================================================
    /// export a heightfield data to a file.
    /// \param[in] user_arguments user command
    void export_heightfield_by_gui(Arguments_get_wrapper& user_arguments);
    /// get export heightfield candidates name list
    /// \param[in] user_arguments user command
    void get_export_heightfield_name_list(Arguments_set_wrapper& response_arguments);
    /// set exporting heightfield GUI index change
    /// \param[in] user_arguments user command
    void set_export_heightfield_index_change(Arguments_get_wrapper& user_arguments);

    //======================================================================
    // export volume data related methods
    //======================================================================
    /// export a volume data to a file.
    /// \param[in] user_arguments user command
    void export_volume_by_gui(Arguments_get_wrapper& user_arguments);
    /// get export volume candidates name list
    /// \param[in] user_arguments user command
    void get_export_volume_name_list(Arguments_set_wrapper& response_arguments);

    //======================================================================
    // attribute generation/compute volume related methods
    //======================================================================
    /// copy the volume under the same static group. Called by GUI.
    /// \param[in] user_arguments user arguments from the flash GUI
    void attrgen_copy_volume_by_gui(Arguments_get_wrapper& user_args);
    /// filter the volume. Called by GUI.
    /// \param[in] user_arguments user arguments from the flash GUI
    void attrgen_filter_volume_by_gui(Arguments_get_wrapper& user_args);
    /// edit the volume. Called by GUI.
    /// \param[in] user_arguments user arguments from the flash GUI
    void attrgen_edit_volume_by_gui(Arguments_get_wrapper& user_args);
    /// get attribute generation GUI related information from the server.
    /// \param[in] p_response_args responce argument to GUI
    /// send attribute generation related information to the flash GUI
    void get_attrgen_gui_section_info(Arguments_set_wrapper& p_response_args);

    //======================================================================
    // RTC kernel program parameter buffer manipulation
    //======================================================================
    /// edit RTC parameter buffer. Called by GUI.
    /// \param[in] user_arguments user arguments from the flash GUI
    void rtc_edit_param_buffer_by_gui(Arguments_get_wrapper& user_args);

    //======================================================================
    // Performance monitor logging state
    //======================================================================
    /// get performance logging state of application (GUI will get the value)
    /// \param[in] p_response_args responce argument to GUI
    void get_performance_logging_state_by_gui(Arguments_set_wrapper& p_response_args);

    //======================================================================
    // Recording for regression test
    //======================================================================
    // Start recording all received commands to the given file
    static void start_recording(const std::string& filename);

    // Stop recording and close the recording file
    static void stop_recording();

    // Write the command with its parameters to the recording file
    static void record_received_command(const std::string& command, mi::base::Handle<const mi::IMap> imap);

    static bool           s_is_recording;
    static mi::base::Lock s_recording_lock;
    static std::ofstream  s_recording_file;
};

class Scene_setup_call_event_handler : public Call_event_handler
{
public:
    Scene_setup_call_event_handler(
        Nvindex_rendering_context& irc,
        bool                       admin = true);

    // Handle the remote command call
    virtual bool handle(
        mi::rtmp::IConnection*      /*connection*/,
        const char*                 /*procedure_name*/,
        const mi::IData*            /*command_args*/,
        const mi::IData*            user_args,
        mi::IData**                 response_args);

private:
    bool                     m_admin;
    std::vector<std::string> m_scene_names;
    std::vector<std::string> m_scene_files;
    std::vector<mi::Sint32>  m_scene_cluster_size_min;

    static mi::base::Lock s_scene_setup_lock;
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTMP_CALL_EVENT_HANDLER_H
