/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief NVIDIA IndeX DiCE Bridge Server

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_BRIDGE_SERVER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_BRIDGE_SERVER_H

#include <mi/base/lock.h>

#include <nv/index/iperformance_values.h>
#include <nv/index/icolormap.h>

#include "common/forwarding_logger.h"
#include "common/scene_parameter_transfer.h"
#include "common/statistics_transfer.h"
#include "common/string_dict.h"

#include "bridge_video_stream.h"
#include "nvindex_appdata.h"
#include "nvindex_rendering_context.h"
#include "examiner_manipulator.h"

#include <iostream>

//----------------------------------------------------------------------
// An application session handler that accepts all sessions that provide a given security token.
class Application_session_handler
  : public mi::base::Interface_implement<mi::bridge::IApplication_session_handler>
{
public:
    Application_session_handler( const std::string& security_token)
      : m_security_token( security_token) { }

    bool on_session_connect( mi::bridge::IServer_session* session)
    {
        const char* security_token = session->get_security_token();
        return security_token && m_security_token == security_token;

    }
private:
    std::string m_security_token;
};

//----------------------------------------------------------------------
// The server-side version of a Bridge job that accesses a database element identified by the tag,
// follows the reference, and calls get_data() called on that element. It returns the sum of the
// return values of that call and the get_data() calls of the serializables in the accessed database
// element and this job.
class Bridge_job_serverside : public mi::bridge::Server_job<0x4e368363,0x49f6,0x4451,0x9f,0x10,0xc,0xc8,0xc3,0x3e,0xba,0x88>
{
public:
    /// constructor
    ///
    /// \param[in] irc index rendering context
    /// \param[in] is_out_command when true, communication command will output to the stderr.
    Bridge_job_serverside(Nvindex_rendering_context* irc, bool is_out_command)
      : 
        m_irc(irc),
        m_scene_parameter_transfer(new nv::index_common::Scene_parameter_transfer),
        m_statistics_transfer(new nv::index_common::Statistics_transfer),
        m_is_out_command(is_out_command),  // for debug
        m_client_command(),
        m_client_argument()
    {
        // empty
    }

    virtual void deserialize( mi::neuraylib::IDeserializer* deserializer);

    void execute(
        mi::bridge::IServer_transaction*    server_transaction,
        mi::neuraylib::ISerializer*         serializer,
        mi::bridge::IServer_job_info*);

    

private: 
    /// Is output mode
    bool is_out_command() const;

    /// Receive client command and arg
    void receive_command(
        mi::neuraylib::IDeserializer* deserializer,
        std::string& com_str,
        std::string& arg_str);

private:
    /// command: client::get_camera_parameter
    /// Send back the performance value to the bridge client
    void serialize_performance_value(mi::neuraylib::ISerializer* serializer);
    /// Send back the result string to the bridge client
    void serialize_string_to_client(mi::neuraylib::ISerializer* serializer, const std::string& str);

private:
    /// command: client::set_camera_parameter
    /// Store camera parameter to internal string
    void store_camera_parameter_to_string(mi::bridge::IServer_transaction* server_transaction);
    /// Set IndeX camera parameter from camera parameter string 
    void set_index_camera_parameter_from_string(
        mi::bridge::IServer_transaction* server_transaction,
        const std::string& camera_param);
private:
    void store_animation_parameter_to_string();

private:
    // command: client::get_colormap_data
    /// Store the current colormap contents to this object
    ///
    /// \param[in] server_transaction
    /// \return true if the colomap gets successfully 
    bool store_colormap_contents(mi::bridge::IServer_transaction* server_transaction);
    /// Send colormap contents to the client 
    void send_colormap_contents_to_client(mi::neuraylib::ISerializer* serializer);

private:
    // command: client::set_colormap_data
    /// Receive the colormap buffer contents from the client
    void receive_colormap_buffer_contents(mi::neuraylib::IDeserializer* deserializer);
    /// Set the received colormap buffer contents to IndeX
    void set_index_colormap(mi::bridge::IServer_transaction* server_transaction, 
                            const std::string& client_command_arg);

private:
    // command: client::get_colormap_list_size
    /// get colormap list size as a string
    std::string get_colormap_list_size();

private:
    // command: client::get_scene_bounding_box, client::get_scene_region_of_interest,
    //          client::set_scene_bounding_box

    /// Get scene bounding box as a string
    ///
    /// \param[in] session_tag session tag
    /// \param[in] server_transaction server transaction to access to the scene
    /// \return scene bounding box information as a string
    std::string get_index_scene_bounding_box(
        mi::neuraylib::Tag               session_tag,        
        mi::bridge::IServer_transaction* server_transaction);

    /// Get scene region of interest
    ///
    /// \param[in] session_tag session tag
    /// \param[in] server_transaction server transaction to access to the scene
    /// \return scene region_of_interest information as a string
    std::string get_index_scene_region_of_interest(
        mi::neuraylib::Tag               session_tag,
        mi::bridge::IServer_transaction* server_transaction);


    /// set scene region of interest
    ///
    /// \param[in] session_tag session tag
    /// \param[in] client_command_arg client argument string, scene region of interest.
    /// \param[in] server_transaction server transaction to access to the scene
    /// \return result status
    std::string set_index_scene_region_of_interest(
        mi::neuraylib::Tag               session_tag,
        const std::string&               client_command_arg, 
        mi::bridge::IServer_transaction* server_transaction);

private: 
    Nvindex_rendering_context*  m_irc;

    // local objects
    mi::base::Handle<nv::index_common::Scene_parameter_transfer>  m_scene_parameter_transfer;
    mi::base::Handle<nv::index_common::Statistics_transfer>       m_statistics_transfer;

    /// colormap local storage for sending to the client
    std::vector<mi::math::Color_struct>                           m_colormap_transfer;

    /// output the communication command and data
    bool m_is_out_command;
    /// command that client send to here
    std::string m_client_command;
    /// string that client send to here
    std::string m_client_argument;
    /// string that the result of request
    std::string m_result_str;
};


class Bridge_job_serverside_factory : public mi::base::Interface_implement<mi::neuraylib::IUser_class_factory>
{
public:

    Bridge_job_serverside_factory(Nvindex_rendering_context* irc, bool is_output_command)
        : 
        m_irc(irc),
        m_is_output_command(is_output_command)
    {
        // empty
    }

    virtual mi::base::IInterface* create(
        mi::neuraylib::ITransaction*    transaction, 
        mi::Uint32                      argc,
        const mi::base::IInterface*     argv[])
    {
        return new Bridge_job_serverside(m_irc, m_is_output_command);
    }

private: 
    Nvindex_rendering_context*  m_irc;
    /// True when output the each command
    bool m_is_output_command;
};


//----------------------------------------------------------------------
// The server-side Bridge initializtion job

class Bridge_init_job_serverside : public mi::bridge::Server_job<0x1ed5bcf,0xf720,0x42cf,0xb5,0xde,0x72,0x8d,0x2f,0x1d,0x95,0xa8>
{
public:
    Bridge_init_job_serverside(Nvindex_rendering_context *irc)
        : 
        m_video_context_id(0),
        m_video_codec(),
        m_video_frame_rate(0),
        m_video_bit_rate(0),
        m_video_max_pending_frames(0),
        m_irc(irc)
    {
        // empty
    }
    
    virtual void deserialize( mi::neuraylib::IDeserializer* deserializer)
    {
        deserializer->read( &m_video_context_id);
        
        mi::Uint32 nb_elements = 0;
        deserializer->read(&nb_elements, 1);
        m_video_codec.resize(nb_elements);
        deserializer->read(reinterpret_cast<mi::Uint8*>(&m_video_codec[0]), nb_elements);

        deserializer->read( &m_video_frame_rate);
        deserializer->read( &m_video_bit_rate);
        deserializer->read( &m_video_max_pending_frames);
        
    }

    void execute(
        mi::bridge::IServer_transaction* server_transaction,
        mi::neuraylib::ISerializer* serializer,
        mi::bridge::IServer_job_info*)
    {
        mi::Uint64 result = 0; // success
        
        mi::bridge::IServer_session *session = server_transaction->get_session();
        
        Bridge_video_stream *video_source = new Bridge_video_stream(session, m_video_context_id, m_irc);
        m_irc->m_bridge_video_source = video_source;
        
        video_source->set_stream_format(m_video_codec);
        video_source->set_stream_bitrate(m_video_bit_rate);
        video_source->set_stream_fps(m_video_frame_rate);
        video_source->set_max_pending_frames(m_video_frame_rate);
        
        const mi::math::Vector<mi::Sint32, 2> main_window_resolution =
            Nvindex_AppData::instance()->get_user_interaction(0)->
            get_examiner()->get_main_window_resolution();

        serializer->write( &result);
        serializer->write( &main_window_resolution.x);
        serializer->write( &main_window_resolution.y);
    }

private:
    mi::Sint32                  m_video_context_id;
    std::string                 m_video_codec;
    mi::Uint32                  m_video_frame_rate;
    mi::Uint32                  m_video_bit_rate;
    mi::Uint32                  m_video_max_pending_frames;
    
    Nvindex_rendering_context*  m_irc;
};

class Bridge_init_job_serverside_factory : public mi::base::Interface_implement<mi::neuraylib::IUser_class_factory>
{
public:

    Bridge_init_job_serverside_factory(Nvindex_rendering_context *irc)
        : m_irc(irc)
    {}

    virtual mi::base::IInterface* create(mi::neuraylib::ITransaction* transaction, 
        mi::Uint32 argc, const mi::base::IInterface* argv[])
    {
        return new Bridge_init_job_serverside(m_irc);
    }

private: 
    Nvindex_rendering_context*  m_irc;
};

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_BRIDGE_SERVER_H
