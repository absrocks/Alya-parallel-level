/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger global/per_host csv output implementation

#include "performance_logger_per_host_csv.h"

#include "utilities.h"

//----------------------------------------------------------------------
Performance_logger_per_host_csv::Performance_logger_per_host_csv()
    :
    Performance_logger_base(),
    m_is_system_global_logging_mode(true)
{
    // always output the frame_id and host_id
    m_per_host_typename_entry.push_back("frame_id");
    m_per_host_typename_entry.push_back("host_id");
}

//----------------------------------------------------------------------
Performance_logger_per_host_csv::~Performance_logger_per_host_csv()
{
    // empty
}

//----------------------------------------------------------------------
void Performance_logger_per_host_csv::set_system_global_logging(bool is_system_global)
{
    m_is_system_global_logging_mode = is_system_global;
}

//----------------------------------------------------------------------
bool Performance_logger_per_host_csv::is_system_global_logging() const
{
    return m_is_system_global_logging_mode;
}

//----------------------------------------------------------------------
void Performance_logger_per_host_csv::append_header(
    nv::index_common::String_dict & app_proj)
{
    append_caption_entry(app_proj);
    append_header_entry();
}

//----------------------------------------------------------------------
void Performance_logger_per_host_csv::append_performance_log(
    nv::index::IPerformance_values* performance_value,
    mi::Sint32                      frame_num)
{
    assert(performance_value != 0);
    if(!this->good()){
        // no stream opened
        return;
    }
    assert(this->get_os() != 0);
    std::ostream & os = *(this->get_os());

    if(is_system_global_logging()){
        // system global logging
        const mi::Sint32 host_id = 0;
        append_perf_log_one_host(os, performance_value, frame_num, host_id);
    }
    else{
        // per-host logging
        const mi::Size host_count = performance_value->get_nb_host_ids();
        if(host_count == 0){
            os << "No host ids. All hosts are idle and no per-host performance values. Empty scene?";
        }
        else{
            // retrieve host id array
            mi::Uint32 * host_id_array = new mi::Uint32[host_count];
            {
                performance_value->get_host_id_array(host_id_array);
                for(mi::Size i = 0; i < host_count; ++i){
                    append_perf_log_one_host(os, performance_value, frame_num, host_id_array[i]);
                }
            }
            delete [] host_id_array;
        }
    }

    if(this->is_flush_each_write()){
        os.flush();
    }
}


//----------------------------------------------------------------------
void Performance_logger_per_host_csv::append_caption_entry(
    nv::index_common::String_dict & app_dict)
{
    if(!this->good()){
        return;
    }
    std::ostream & os = *(this->get_os());

    // Additional information
    os << "# Performance_logger_name: Performance_logger_per_host_csv "
       << (is_system_global_logging() ? "global_system" : "per_host") << "\n";
    os << "# Title: "
       << app_dict.get("app::performance::title", "(none)") << "\n";
    os << "# Description: "
       << app_dict.get("app::performance::description", "(none)") << "\n";

    // record versions
    os << "# dice_version: "   << app_dict.get("app::info::dice_version",   "not available.") << "\n";
    os << "# index_revision: " << app_dict.get("app::info::index_revision", "not available.") << "\n";
    os << "# driver_version: " << app_dict.get("app::info::driver_version", "not available.") << "\n";
    os << "# cuda_runtime_version: " << app_dict.get("app::info::cuda_version",  "not available.") << "\n";

    // span renderer
    os << "# Span_renderer: "  << app_dict.get("span_renderer_name", "(none).") << "\n";

    // recording time
    os << "# Recordtime: " << current_system_calender_str() << "\n";

    // host information
    os << "# Localhost: " << app_dict.get("app::performance::local_host_name", "localhost") << "\n"
       << "# Hosts: "     << app_dict.get("app::performance::host_list") << "\n";

    os << "# build_environment: " << get_app_build_info() << "\n";
    
    // video encoder
    os << "# RTMP Video_encoder: " << app_dict.get("dice::rtmp_video_streaming::video_codec", "(none)") << "\n";

    // Screen resolution
    os << "# Screen_resolution: " << app_dict.get("index::canvas_resolution") << "\n";
}

//----------------------------------------------------------------------
void Performance_logger_per_host_csv::append_header_entry()
{
    if(!this->good()){
        return;
    }
    std::ostream & os = *(this->get_os());

    // The header line is also a special line
    os << "# Performance_item: ";

    // per_host/system related performance values
    for(std::vector< std::string >::const_iterator ei = m_per_host_typename_entry.begin();
        ei != m_per_host_typename_entry.end(); ++ei)
    {
        os << (*ei) << ",";
    }

    // IPerformance_values related performance values
    const std::vector<std::string> & perf_entry_vec   = peek_performance_logging_entry();
    for(std::vector< std::string >::const_iterator ei = perf_entry_vec.begin();
        ei != perf_entry_vec.end(); ++ei)
    {
        os << (*ei) << ",";
    }
    os << "\n";
}

//----------------------------------------------------------------------
void Performance_logger_per_host_csv::append_per_host_proper_perf_log(
    std::ostream & os, 
    mi::Sint32 frame_num,
    mi::Sint32 host_id)
{
    assert(os.good());
    os << frame_num << "," << host_id << ",";    
}

//----------------------------------------------------------------------
void Performance_logger_per_host_csv::append_perf_log_one_host(
    std::ostream & os, 
    nv::index::IPerformance_values* performance_value,
    mi::Sint32 frame_num,
    mi::Sint32 host_id)
{
    assert(performance_value != 0);

    // special entries
    append_per_host_proper_perf_log(os, frame_num, host_id);

    // IPerformance_values
    const std::vector< std::string > & iperf_typename_vec = peek_performance_logging_entry();

    for(std::vector< std::string >::const_iterator ei = iperf_typename_vec.begin();
        ei != iperf_typename_vec.end(); ++ei)
    {
        const std::string entry_name(*ei);
        const Performance_data_typename_e dtype = Performance_logger_base::get_data_type(entry_name);
        if(dtype == PDT_Uint64){
            os << performance_value->get(entry_name.c_str(), host_id) << ",";
        }
        else if(dtype == PDT_Float32){
            os << performance_value->get_time(entry_name.c_str(), host_id) << ",";
        }
        else{
            ERROR_LOG << "Performance_logger_per_host_csv::append_performance_value_perf_log: "
                      << "Unsupported data type [" << dtype << "] of [" << entry_name << "]";
        }
    }
    os << "\n";
}

//----------------------------------------------------------------------
