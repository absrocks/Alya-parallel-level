/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief NVIDIA IndeX DiCE Bridge server side implementation

#include "bridge_server.h"

#include "colormap_util.h"

#include "common/clock_pulse_generator.h"
#include "common/string_dict.h"
#include "common/type_conversion_utility.h"

#include "scene_utility.h"

#include <iterator>

//----------------------------------------------------------------------
void Bridge_job_serverside::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    receive_command(deserializer, m_client_command, m_client_argument);
    m_scene_parameter_transfer->deserialize(deserializer);

    // additional data from the client
    if (m_client_command == "client::set_colormap_data")
    {
        receive_colormap_buffer_contents(deserializer);
    }

    mi::base::Handle<nv::index::IClock_pulse_generator> iclock(
        m_irc->m_iindex_session->get_clock());
    if (iclock)
    {
        mi::base::Handle<nv::index_common::Clock_pulse_generator> clock(
              iclock->get_interface<nv::index_common::Clock_pulse_generator>());
        clock->set_tick(m_scene_parameter_transfer->get_animation().m_t_current);
        
        if(! (m_scene_parameter_transfer->get_animation().m_t_start==0.
                && m_scene_parameter_transfer->get_animation().m_t_end==0) )
        {
            clock->set_start(m_scene_parameter_transfer->get_animation().m_t_start);
            clock->set_end(m_scene_parameter_transfer->get_animation().m_t_end);
    
            /*INFO_LOG << "Set clock: t="
                     << m_scene_parameter_transfer->get_animation().m_t_current << " in ["
                     << m_scene_parameter_transfer->get_animation().m_t_start << ","
                     << m_scene_parameter_transfer->get_animation().m_t_end << "].";*/
        }

        clock->use_external_setter(true);
    }
}


//----------------------------------------------------------------------
void Bridge_job_serverside::execute(
    mi::bridge::IServer_transaction*    server_transaction,
    mi::neuraylib::ISerializer*         serializer,
    mi::bridge::IServer_job_info*)
{
    if (m_client_command == "client::get_camera_parameter")
    {
        // get the camera parameters and store to the result string
        store_camera_parameter_to_string(server_transaction);
        // send back the result string
        serialize_string_to_client(serializer, m_result_str);
    }
    else if (m_client_command == "client::set_camera_parameter")
    {
        // set camera parameter from client to IndeX camera
        set_index_camera_parameter_from_string(server_transaction, m_client_argument);
        // send back the result string
        serialize_string_to_client(serializer, m_result_str);        
    }
    else if (m_client_command == "client::get_animation_parameter")
    {
        // get the camera parameters and store to the result string
        store_animation_parameter_to_string();
        // send back the result string
        serialize_string_to_client(serializer, m_result_str);
    }
    else if (m_client_command == "client::get_colormap_data")
    {
        // set camera parameter from client to IndeX camera
        const bool is_colormap_ready = store_colormap_contents(server_transaction);

        // send back the result string
        serialize_string_to_client(serializer, m_result_str);        
        if (is_colormap_ready)
        {
            // send back the colormap contents to the client
            send_colormap_contents_to_client(serializer);
        }
    }
    else if (m_client_command == "client::set_colormap_data")
    {
        // has got the colormap from the client, set it to IndeX
        set_index_colormap(server_transaction, m_client_argument);
        // send back the result string
        serialize_string_to_client(serializer, m_result_str);        
    }
    else if (m_client_command == "client::get_colormap_list_size")
    {
        m_result_str = get_colormap_list_size();
        // send back the result string
        serialize_string_to_client(serializer, m_result_str);                
    }
    else if (m_client_command == "client::get_scene_bounding_box")
    {
        m_result_str = get_index_scene_bounding_box(m_irc->m_session_tag, server_transaction);
        // send back the result string
        serialize_string_to_client(serializer, m_result_str);                
    }
    else if (m_client_command == "client::get_scene_region_of_interest")
    {
        m_result_str = get_index_scene_region_of_interest(m_irc->m_session_tag, server_transaction);
        // send back the result string
        serialize_string_to_client(serializer, m_result_str);                
    }
    else if (m_client_command == "client::set_scene_region_of_interest")
    {
        m_result_str = set_index_scene_region_of_interest(
            m_irc->m_session_tag, m_client_argument, server_transaction);
        // send back the result string
        serialize_string_to_client(serializer, m_result_str);                
    }
    else
    {
        std::cerr << "execute: Unknown command: " << m_client_command << std::endl;
        // no send back the result string.
    }
    
    // currently always send back the performance value regardless the
    // client command.
    serialize_performance_value(serializer);
}

//----------------------------------------------------------------------
bool Bridge_job_serverside::is_out_command() const
{
    return m_is_out_command;
}

//----------------------------------------------------------------------
void Bridge_job_serverside::receive_command(
    mi::neuraylib::IDeserializer* deserializer,
    std::string& com_str,
    std::string& arg_str)
{
    // receive command
    {
        mi::Uint32 nb_elements = 0;
        deserializer->read(&nb_elements, 1);
        assert(nb_elements > 0);
        com_str.resize(nb_elements);
        deserializer->read(reinterpret_cast<mi::Uint8*>(&(com_str[0])), nb_elements);
    }

    // receive arg
    {
        mi::Uint32 nb_elements = 0;
        deserializer->read(&nb_elements, 1);
        assert(nb_elements > 0);
        arg_str.resize(nb_elements);
        deserializer->read(reinterpret_cast<mi::Uint8*>(&(arg_str[0])), nb_elements);
    }

    if (is_out_command())
    {
        std::cerr << "command: " << com_str << std::endl;
        std::cerr << "arg: "     << arg_str << std::endl;
    }
}

//----------------------------------------------------------------------
void Bridge_job_serverside::serialize_performance_value(mi::neuraylib::ISerializer* serializer)
{
    // NVIDIA IndeX Application Data is still a singleton.
    // Lock and grab performance values from there:
    //
    mi::base::Lock::Block block(&Nvindex_AppData::instance()->m_performance_values_lock);
    if(Nvindex_AppData::instance()->m_performance_values)
    {
        nv::index::IPerformance_values* perf = Nvindex_AppData::instance()->m_performance_values.get();
        const mi::Float32 fps               = perf->get_time("frames_per_second");
        const mi::Float32 total_rendering   = perf->get_time("time_total_rendering");
        const mi::Float32 total_compositing = perf->get_time("time_total_compositing");
                
        const mi::Float32 rendering_only    = perf->get_time("time_rendering_only");

        const mi::Float32 gpu_upload        = perf->get_time("time_gpu_upload");
        if (false && gpu_upload)
        {
            // LET US KNOW!
            WARN_LOG << "Time for data uploading to the GPU: " << gpu_upload;
        }
        nv::index_common::Statistics_transfer::Statistics_values values;
        values.m_fps                    = fps;
        values.m_total_rendering_time   = total_rendering;
        values.m_total_compositing_time = total_compositing;

        values.m_rendering_time_only    = rendering_only;
        values.m_GPU_uploading_time     = gpu_upload;

        m_statistics_transfer->set_values(values);
    }

    // Something has to be serialized as the client expects data .... for now!
    m_statistics_transfer->serialize(serializer);
}

//----------------------------------------------------------------------
void Bridge_job_serverside::serialize_string_to_client(
    mi::neuraylib::ISerializer* serializer,
    const std::string& str)
{
    if (is_out_command())
    {
        std::cerr << "serialize_string: " << str << std::endl;
    }

    mi::Uint32 nb_elements = static_cast<mi::Uint32>(str.size());
    serializer->write(&nb_elements, 1);
    assert(nb_elements > 0);
    serializer->write(reinterpret_cast<const mi::Uint8*>(&(str[0])), nb_elements);
}

//----------------------------------------------------------------------
void Bridge_job_serverside::store_animation_parameter_to_string()
{
    mi::base::Handle<nv::index::IClock_pulse_generator> iclock(
        m_irc->m_iindex_session->get_clock());
    if (iclock)
    {
        mi::base::Handle<nv::index_common::Clock_pulse_generator> clock(
              iclock->get_interface<nv::index_common::Clock_pulse_generator>());
        mi::Float64 start = clock->get_start();
        mi::Float64 end   = clock->get_end();

        nv::index_common::String_dict anim_param_dict;
        anim_param_dict.insert("result", nv::index_common::to_string(true)); // return as success
        anim_param_dict.insert("start",  nv::index_common::to_string(start));
        anim_param_dict.insert("end",    nv::index_common::to_string(end));

        std::stringstream sstr;
        anim_param_dict.write(sstr, "");

        m_result_str = sstr.str();
    }
    else
    {
        INFO_LOG << "No animation parameters to store. Please check app::clock_pulse::interval "
                 << "entry exists in your project file if you needed animation.";
        m_result_str = 
            "result = false\n"  // This tells to the client no animation data set at the server.
            "start = 0\n"
            "end = 1\n";
    }
}

//----------------------------------------------------------------------
void Bridge_job_serverside::store_camera_parameter_to_string(
    mi::bridge::IServer_transaction* server_transaction)
{
    assert(server_transaction != 0);
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        server_transaction->get_database_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc->m_session_tag));
        assert(session.is_valid_interface());

        mi::base::Handle<nv::index::IScene> scene(dice_transaction->edit<nv::index::IScene>(session->get_scene()));
        assert(scene.is_valid_interface());

        const mi::neuraylib::Tag& cam_tag = scene->get_camera();
        assert(cam_tag.is_valid());
        mi::base::Handle<const nv::index::ICamera> cam(dice_transaction->access<const nv::index::ICamera>(cam_tag));
        assert(cam.is_valid_interface());        

        // get the camera parameter
        Camera_parameter camparam;
        camparam.copy_camera_parameter_to_this(cam.get());

        const std::string string_dict_prefix = "";
        nv::index_common::String_dict cam_param_dict;
        camparam.get_parameter_to_string_dict(string_dict_prefix, cam_param_dict);

        std::stringstream sstr;
        cam_param_dict.write(sstr, "");

        m_result_str = sstr.str();
    }
    // dice_transaction->commit(); This crashes the client!
}

//----------------------------------------------------------------------
void Bridge_job_serverside::set_index_camera_parameter_from_string(
    mi::bridge::IServer_transaction* server_transaction,
    const std::string& camera_param)
{
    assert(server_transaction != 0);

    nv::index_common::String_dict cam_param_dict;

    std::stringstream sstr(camera_param);
    cam_param_dict.read(sstr);

    assert(cam_param_dict.is_defined("from"));
    assert(cam_param_dict.is_defined("dir"));
    assert(cam_param_dict.is_defined("up"));

    const mi::math::Vector<mi::Float32, 3> from = nv::index_common::get_vec_float32_3(cam_param_dict.get("from", "0 0  0"));
    const mi::math::Vector<mi::Float32, 3> dir  = nv::index_common::get_vec_float32_3(cam_param_dict.get("dir",  "0 0 -1"));
    const mi::math::Vector<mi::Float32, 3> up   = nv::index_common::get_vec_float32_3(cam_param_dict.get("up",   "0 1  0"));

    assert(server_transaction != 0);
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        server_transaction->get_database_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc->m_session_tag));
        assert(session.is_valid_interface());

        mi::base::Handle<nv::index::IScene> scene(dice_transaction->edit<nv::index::IScene>(session->get_scene()));
        assert(session.is_valid_interface());

        const mi::neuraylib::Tag& cam_tag = scene->get_camera();
        assert(cam_tag.is_valid());

        mi::base::Handle<nv::index::ICamera> cam(dice_transaction->edit<nv::index::ICamera>(cam_tag));
        assert(cam.is_valid_interface());

        cam->set(from, dir, up);
    }
    // dice_transaction->commit(); This crashes the client!

    m_result_str = "result = true\n"; // successfully finished
}


//----------------------------------------------------------------------
bool Bridge_job_serverside::store_colormap_contents(mi::bridge::IServer_transaction* server_transaction)
{
    nv::index_common::String_dict result_dict;
    m_result_str = "result = false";
    result_dict.insert("result", "true");    

    const std::string mn = "Bridge_job_serverside::store_colormap:"; // method name

    // For demo, always the first volume
    const mi::Sint32 volume_idx = 0;
    const mi::neuraylib::Tag volume_tag = Nvindex_AppData::instance()->get_volume_tag_from_idx(volume_idx);
    if (!volume_tag.is_valid())
    {
        ERROR_LOG << mn << "no volume in the scene found.";
        return false;
    }

    ERROR_LOG << mn << "Currently this is the first volume not requested colormap index as set command. FIXME";

    assert(server_transaction != 0);
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        server_transaction->get_database_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        mi::base::Handle<const nv::index::IRegular_volume> volume(
            dice_transaction->access<const nv::index::IRegular_volume>(volume_tag));
        assert(volume.is_valid_interface());
        result_dict.insert("volume_tag", nv::index_common::to_string(volume_tag.id));

        // INFO_LOG << "DEBUG: found a volume: " << volume_tag << " for colormap.";
        const mi::neuraylib::Tag colormap_tag = volume->assigned_colormap();
        if (!colormap_tag.is_valid())
        {
            ERROR_LOG << mn << "no colormap is assigned for volume: " << volume_tag;
            // no need to commit/abort this transaction
            return false;
        }
        // INFO_LOG << "DEBUG: found a colormap tag: " << colormap_tag;

        // keep the colormap local in this object for later transfer
        mi::base::Handle<nv::index::IColormap> colormap(
            dice_transaction->edit<nv::index::IColormap>(colormap_tag));
        assert(colormap.is_valid_interface());

        // store the colormap to local
        const mi::Size nb_colormap_entry = colormap->get_number_of_entries();
        assert(nb_colormap_entry > 0);
        m_colormap_transfer.clear();
        for (mi::Size i = 0; i < nb_colormap_entry; ++i)
        {
            m_colormap_transfer.push_back(colormap->get_colormap()[i]);
            // std::cerr << "Send Colormap to client[" << i << "]: " << m_colormap_transfer[i].r
            //           << " " << m_colormap_transfer[i].g
            //           << " " << m_colormap_transfer[i].b            
            //           << " " << m_colormap_transfer[i].a
            //           << std::endl;
        }
        result_dict.insert("colomap_entry_size", nv::index_common::to_string(nb_colormap_entry));

        mi::Float32 domain_min = 0.0f;
        mi::Float32 domain_max = 0.0f;
        colormap->get_domain(domain_min, domain_max);
        nv::index::IColormap::Domain_boundary_mode dmode = colormap->get_domain_boundary_mode();

        mi::Uint32 trans_min = 0;
        mi::Uint32 trans_max = 0;
        colormap->get_non_transparent_range(trans_min, trans_max);

        result_dict.insert("colomap_domain_min",  nv::index_common::to_string(domain_min));
        result_dict.insert("colomap_domain_max",  nv::index_common::to_string(domain_max));
        result_dict.insert("colomap_domain_mode", nv::index_common::to_string(static_cast<mi::Sint32>(dmode)));
        result_dict.insert("colomap_transparent_min", nv::index_common::to_string(trans_min));
        result_dict.insert("colomap_transparent_max", nv::index_common::to_string(trans_max));
    }
    // don't commit the transaction

    // set the result string
    std::stringstream sstr;
    result_dict.write(sstr, "");
    m_result_str = sstr.str();

    return true;
}

//----------------------------------------------------------------------
void Bridge_job_serverside::send_colormap_contents_to_client(
    mi::neuraylib::ISerializer* serializer)
{
    assert(serializer != 0);
    assert(!(m_colormap_transfer.empty()));    

    const mi::Uint32 nb_elements = static_cast<mi::Uint32>(m_colormap_transfer.size());
    assert(nb_elements > 0);
    serializer->write(&nb_elements, 1);
    serializer->write((&(m_colormap_transfer[0].r)), nb_elements * 4);
}

//----------------------------------------------------------------------
void Bridge_job_serverside::receive_colormap_buffer_contents(mi::neuraylib::IDeserializer* deserializer)
{
    m_colormap_transfer.clear();

    mi::Uint32 nb_entries = 0;
    deserializer->read(&nb_entries, 1);
    m_colormap_transfer.resize(nb_entries);
    assert(nb_entries > 0);
    deserializer->read(&m_colormap_transfer[0].r, nb_entries * 4);

    // for (mi::Size i = 0; i < m_colormap_transfer.size(); ++i)
    // {
    //     std::cerr << "Recv Colormap from client[" << i << "]: " << m_colormap_transfer[i].r
    //               << " " << m_colormap_transfer[i].g
    //               << " " << m_colormap_transfer[i].b            
    //               << " " << m_colormap_transfer[i].a
    //               << std::endl;
    // }
}

//----------------------------------------------------------------------
void Bridge_job_serverside::set_index_colormap(
    mi::bridge::IServer_transaction* server_transaction,
    const std::string& client_command_arg)
{
    const std::string mn = "Bridge_job_serverside::set_index_colormap:"; // method name

    nv::index_common::String_dict param = nv::index_common::get_string_dict_from_string(client_command_arg);
    // Necessary defined keys
    char const * const p_def_param_key[] = {
        "colormap_index", 
        "colormap_domain_min",
        "colormap_domain_max",
        "colormap_domain_mode",    
        // "colormap_transparent_min",
        // "colormap_transparent_max",
        0,
    };
    std::vector< std::string > undef_list;
    if (!is_all_keys_defined(param, p_def_param_key, &undef_list))
    {
        std::stringstream sstr;
        std::copy(undef_list.begin(), undef_list.end(), std::ostream_iterator< std::string >(sstr, " "));
        ERROR_LOG << mn << "Undefined camera key entries (use default):\n" << sstr.str();
    }

    mi::Sint32 colormap_index = nv::index_common::get_sint32(param.get("colormap_index", "0"));
    assert(colormap_index >= 0);
    assert(colormap_index < static_cast<mi::Sint32>(get_number_of_colormap()));
    mi::Float32 domain_min  = nv::index_common::get_float32(param.get("colormap_domain_min", "0.0"));
    mi::Float32 domain_max  = nv::index_common::get_float32(param.get("colormap_domain_max", "1.0"));
    mi::Sint32  domain_mode = nv::index_common::get_sint32(param.get("colormap_domain_mode", "0"));
    assert((domain_mode == 0) || (domain_mode == 1));

    mi::neuraylib::Tag colormap_tag = get_colormap_tag(static_cast<mi::Uint32>(colormap_index));
    if (!colormap_tag.is_valid())
    {
        ERROR_LOG << mn << "no colormap of index: " << colormap_index;
        m_result_str = "result = false\n"; // failed
        return;
    }

    assert(server_transaction != 0);
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        server_transaction->get_database_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        // keep the colormap local in this object for later transfer
        mi::base::Handle<nv::index::IColormap> colormap(
            dice_transaction->edit<nv::index::IColormap>(colormap_tag));
        assert(colormap.is_valid_interface());

        // store the colormap to local
        assert(!(m_colormap_transfer.empty()));
        colormap->set_colormap(&(m_colormap_transfer[0]), m_colormap_transfer.size());

        colormap->set_domain(domain_min, domain_max);
        nv::index::IColormap::Domain_boundary_mode dmode =
            static_cast<nv::index::IColormap::Domain_boundary_mode>(domain_mode);
        colormap->set_domain_boundary_mode(dmode);
    }
    // don't commit the transaction

    m_result_str = "result = true\n"; // successfully finished
}

//----------------------------------------------------------------------
std::string Bridge_job_serverside::get_colormap_list_size()
{
    mi::Uint32 nb_colormap = get_number_of_colormap();
    std::stringstream sstr;
    sstr << "result = true\n"
         << "colormap_list_size = " << nb_colormap << "\n";

    return sstr.str();
}

//----------------------------------------------------------------------
/// get global scope from a given dice transaction
///
/// \param[in]  current_dice_transaction current dice transaction
/// \param[out] global_scope             (out) global scope
static void get_global_scope(
    const mi::neuraylib::IDice_transaction*  current_dice_transaction,
    mi::base::Handle<mi::neuraylib::IScope>& global_scope)
{
    assert(current_dice_transaction != 0);
    mi::base::Handle<mi::neuraylib::IScope> cur_scope(current_dice_transaction->get_scope());
    while (true)
    {
        mi::base::Handle<mi::neuraylib::IScope> parent_scope(cur_scope->get_parent());
        if (!parent_scope.is_valid_interface())
        {
            global_scope.swap(cur_scope);
            break;
        }
        cur_scope.swap(parent_scope);
        assert(cur_scope.is_valid_interface());
    }
    assert(global_scope.is_valid_interface());
}

//----------------------------------------------------------------------
std::string Bridge_job_serverside::get_index_scene_bounding_box(
    mi::neuraylib::Tag               session_tag,
    mi::bridge::IServer_transaction* server_transaction)
{
    assert(server_transaction != 0);
    mi::base::Handle<mi::neuraylib::IDice_transaction> current_dice_transaction(
        server_transaction->get_database_transaction<mi::neuraylib::IDice_transaction>());
    assert(current_dice_transaction.is_valid_interface());
    std::stringstream sstr;
    {
        mi::base::Handle<mi::neuraylib::IScope> global_scope;
        get_global_scope(current_dice_transaction.get(), global_scope);

        mi::base::Handle<mi::neuraylib::IDice_transaction> global_scope_dice_transaction(
            global_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(global_scope_dice_transaction.is_valid_interface());
        {
            const mi::math::Bbox<mi::Float32, 3> scene_bbox_f32 = get_scene_bounding_box(session_tag, global_scope_dice_transaction.get());
            // INFO_LOG << "scene_bbox: " << tmp_scene_bbox_f32 << ", " << scene_bbox_f32;

            assert(!(scene_bbox_f32.empty()));
            const mi::math::Bbox<mi::Sint32, 3>  scene_bbox_s32 =
                nv::index_common::convert_bbox_type<mi::Sint32, mi::Float32, 3>(scene_bbox_f32);
            assert(!(scene_bbox_s32.empty()));

            sstr << "result = true\n"
                 << "scene_bounding_box = "
                 << scene_bbox_s32.min.x << " " << scene_bbox_s32.min.y << " " << scene_bbox_s32.min.z << " " 
                 << scene_bbox_s32.max.x << " " << scene_bbox_s32.max.y << " " << scene_bbox_s32.max.z << "\n";
        }
        global_scope_dice_transaction->commit();
    }

    return sstr.str();
}

//----------------------------------------------------------------------
std::string Bridge_job_serverside::get_index_scene_region_of_interest(
    mi::neuraylib::Tag               session_tag,
    mi::bridge::IServer_transaction* server_transaction)
{
    assert(server_transaction != 0);
    mi::base::Handle<mi::neuraylib::IDice_transaction> current_dice_transaction(
        server_transaction->get_database_transaction<mi::neuraylib::IDice_transaction>());
    assert(current_dice_transaction.is_valid_interface());
    std::stringstream sstr;
    {
        mi::base::Handle<mi::neuraylib::IScope> global_scope;
        get_global_scope(current_dice_transaction.get(), global_scope);

        mi::base::Handle<mi::neuraylib::IDice_transaction> global_scope_dice_transaction(
            global_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(global_scope_dice_transaction.is_valid_interface());
        {
            mi::base::Handle<const nv::index::ISession> session(
                global_scope_dice_transaction->access<nv::index::ISession>(session_tag));
            assert(session.is_valid_interface());

            mi::base::Handle<const nv::index::IScene> scene(
                global_scope_dice_transaction->access<const nv::index::IScene>(session->get_scene()));
            assert(scene.is_valid_interface());

            const mi::math::Bbox<mi::Float32, 3> roi_bbox_f32 = scene->get_clipped_bounding_box();
            assert(!(roi_bbox_f32.empty()));

            const mi::math::Bbox<mi::Sint32, 3>  roi_bbox_s32 =
                nv::index_common::convert_bbox_type<mi::Sint32, mi::Float32, 3>(roi_bbox_f32);
            assert(!(roi_bbox_s32.empty()));
            
            sstr << "result = true\n"
                 << "scene_region_of_interest = "
                 << roi_bbox_s32.min.x << " " << roi_bbox_s32.min.y << " " << roi_bbox_s32.min.z << " " 
                 << roi_bbox_s32.max.x << " " << roi_bbox_s32.max.y << " " << roi_bbox_s32.max.z << "\n";
        }
        global_scope_dice_transaction->commit();
    }

    return sstr.str();
}

//----------------------------------------------------------------------
std::string Bridge_job_serverside::set_index_scene_region_of_interest(
    mi::neuraylib::Tag               session_tag,
    const std::string&               client_command_arg,
    mi::bridge::IServer_transaction* server_transaction)
{
    assert(server_transaction != 0);
    mi::base::Handle<mi::neuraylib::IDice_transaction> current_dice_transaction(
        server_transaction->get_database_transaction<mi::neuraylib::IDice_transaction>());
    assert(current_dice_transaction.is_valid_interface());
    {
        mi::base::Handle<mi::neuraylib::IScope> global_scope;
        get_global_scope(current_dice_transaction.get(), global_scope);

        mi::base::Handle<mi::neuraylib::IDice_transaction> global_scope_dice_transaction(
            global_scope->create_transaction<mi::neuraylib::IDice_transaction>());
        assert(global_scope_dice_transaction.is_valid_interface());
        {
            mi::base::Handle<const nv::index::ISession> session(
                global_scope_dice_transaction->access<nv::index::ISession>(session_tag));
            assert(session.is_valid_interface());

            mi::base::Handle<nv::index::IScene> scene(
                global_scope_dice_transaction->edit<nv::index::IScene>(session->get_scene()));
            assert(scene.is_valid_interface());

            nv::index_common::String_dict param = nv::index_common::get_string_dict_from_string(client_command_arg);            
            assert(param.is_defined("scene_region_of_interest"));
            const mi::math::Bbox<mi::Float32, 3>  roi_bbox_f32 =
                nv::index_common::get_bbox_float32_3(param.get("scene_region_of_interest"));
            assert(!(roi_bbox_f32.empty()));
            scene->set_clipped_bounding_box(roi_bbox_f32);
        }
        global_scope_dice_transaction->commit();
    }

    const std::string ret("result = true\n");
    return ret;
}
//----------------------------------------------------------------------

