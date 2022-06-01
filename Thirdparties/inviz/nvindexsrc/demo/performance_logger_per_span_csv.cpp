/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger per_span csv output implementation

#include "performance_logger_per_span_csv.h"

#include "utilities.h"

//----------------------------------------------------------------------
Performance_logger_per_span_csv::Performance_logger_per_span_csv()
    :
    Performance_logger_base(),
    m_is_system_global_logging_mode(true)
{
    // always output the frame_id and host_id
    m_per_host_typename_entry.push_back("frame_id");
    m_per_host_typename_entry.push_back("span_id");
    m_per_host_typename_entry.push_back("cluster_id");
}

//----------------------------------------------------------------------
Performance_logger_per_span_csv::~Performance_logger_per_span_csv()
{
    // empty
}

//----------------------------------------------------------------------
void Performance_logger_per_span_csv::append_header(
    nv::index_common::String_dict & app_proj)
{
    append_caption_entry(app_proj);
    append_header_entry();
}

//----------------------------------------------------------------------
void Performance_logger_per_span_csv::append_performance_log(
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

    const mi::Uint32 system_host_id = 0; // for query the number of spans
    const mi::Uint64 nb_spans = performance_value->get("nb_horizontal_spans", system_host_id);
    if(nb_spans == 0){
        os << "No spans. Empty scene?";
    }
    else{
        const mi::Uint32 nb_spans_u32 = static_cast<mi::Uint32>(nb_spans);
        for(mi::Uint32 i = 0; i < nb_spans_u32; ++i){
            append_perf_log_one_span(os, performance_value, frame_num, i);
        }
    }

    if(this->is_flush_each_write()){
        os.flush();
    }
}


//----------------------------------------------------------------------
void Performance_logger_per_span_csv::append_caption_entry(
    nv::index_common::String_dict & app_dict)
{
    if(!this->good()){
        return;
    }
    std::ostream & os = *(this->get_os());

    // Additional information
    os << "# Performance_logger_name: Performance_logger_per_span_csv\n";
    os << "# Title: "
       << app_dict.get("app::performance::title", "(none)") << "\n";
    os << "# Description: "
       << app_dict.get("app::performance::description", "(none)") << "\n";

    // record library version
    os << "# Version: "
       << app_dict.get("app::info::dice_version", "no version available.") << "\n";

    // span renderer
    os << "# Span_renderer: "
       << app_dict.get("span_renderer_name", "(none).") << "\n";

    // recording time
    os << "# Recordtime: " << current_system_calender_str() << "\n";

    // hosts
    os << "# Localhost: " << app_dict.get("app::performance::local_host_name", "localhost") << "\n"
       << "# Hosts: " << app_dict.get("app::performance::host_list") << "\n";

    os << "# build_environment: " << get_app_build_info() << "\n";
    
    // RTMP video encoder
    os << "# RTMP_video_encoder: "
       << app_dict.get("dice::rtmp_video_streaming::video_codec", "(none)") << "\n";

    // Screen resolution
    os << "# Screen_resolution: "
      << app_dict.get("index::canvas_resolution") << "\n";
}

//----------------------------------------------------------------------
void Performance_logger_per_span_csv::append_header_entry()
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
void Performance_logger_per_span_csv::append_per_span_proper_perf_log(
    std::ostream & os, 
    mi::Sint32 frame_num,
    mi::Uint32 span_id,
    mi::Uint32 cluster_id)
{
    assert(os.good());
    os << frame_num << "," << span_id << "," << cluster_id << ",";    
}

//----------------------------------------------------------------------
void Performance_logger_per_span_csv::append_perf_log_one_span(
    std::ostream & os, 
    nv::index::IPerformance_values* performance_value,
    mi::Sint32 frame_num,
    mi::Uint32 span_idx)
{
    assert(performance_value != 0);

    // per-span logging
    mi::base::Handle<nv::index::IPer_span_statistics> per_span_stat(
        performance_value->get_per_span_statistics(span_idx));

    if(per_span_stat == 0){
        WARN_LOG << "Performance_logger_per_span_csv::append_per_span_proper_perf_log: invalid span_idx ["
                 << span_idx << "], skip it.";
        return;
    }

    // export extra info: span_id, cluster_id
    // Note: we have now span_idx -> span_id. This holds in the current implementation.
    mi::Uint32 const span_id    = per_span_stat->get_span_id();
    mi::Uint32 const cluster_id = per_span_stat->get_cluster_node_id();

    // output special entries
    append_per_span_proper_perf_log(os, frame_num, span_id, cluster_id);

    const std::vector< std::string > & iperf_typename_vec = peek_performance_logging_entry();

    for(std::vector< std::string >::const_iterator ei = iperf_typename_vec.begin();
        ei != iperf_typename_vec.end(); ++ei)
    {
        std::string const entry_name(*ei);
        const Performance_data_typename_e dtype = Performance_logger_base::get_data_type(entry_name);
        if(dtype == PDT_Uint64){
            os << per_span_stat->get(ei->c_str())      << ",";
        }
        else if(dtype == PDT_Float32){
            os << per_span_stat->get_time(ei->c_str()) << ",";
        }
        else{
            ERROR_LOG << "Performance_logger_per_span_csv::append_perf_log_one_host: "
                      << "Unsupported data type [" << dtype << "] of [" << entry_name << "]";
        }
    }
    os << "\n";
}

//----------------------------------------------------------------------
