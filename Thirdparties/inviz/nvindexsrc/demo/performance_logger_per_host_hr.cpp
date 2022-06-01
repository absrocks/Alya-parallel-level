/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger global/per_host csv output implementation

#include "performance_logger_per_host_hr.h"

#include "utilities.h"

//----------------------------------------------------------------------
Performance_logger_per_host_hr::Performance_logger_per_host_hr()
    :
    Performance_logger_base(),
    m_is_system_global_logging_mode(true),
    m_out_type(HROF_Simple)
{
    // always output the frame_id and host_id
    m_per_host_typename_entry.push_back("frame_id");
    m_per_host_typename_entry.push_back("host_id");
}

//----------------------------------------------------------------------
Performance_logger_per_host_hr::~Performance_logger_per_host_hr()
{
    // empty
}

//----------------------------------------------------------------------
void Performance_logger_per_host_hr::set_system_global_logging(bool is_system_global)
{
    m_is_system_global_logging_mode = is_system_global;
}

//----------------------------------------------------------------------
bool Performance_logger_per_host_hr::is_system_global_logging() const
{
    return m_is_system_global_logging_mode;
}

//----------------------------------------------------------------------
void Performance_logger_per_host_hr::set_output_format(const nv::index_common::String_dict & opt)
{
    const std::string kind_key = (is_system_global_logging() ? "global" : "per_host");
    const std::string out_type = opt.get("app::performance::" + kind_key + "::human_readable_type", "simple");

    if(out_type == "simple"){
        m_out_type = HROF_Simple;
    }
    else if(out_type == "transfer"){
        m_out_type = HROF_Transfer;        
    }
    else{
        WARN_LOG << "Performance_logger_per_host_hr::set_output_format: unknown output type [" 
                 << out_type <<"], set to [simple].";
        m_out_type = HROF_Simple;        
    }
}

//----------------------------------------------------------------------
void Performance_logger_per_host_hr::append_header(
    nv::index_common::String_dict & app_proj)
{
    append_caption_entry(app_proj);
    append_header_entry();
}

//----------------------------------------------------------------------
void Performance_logger_per_host_hr::append_performance_log(
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
void Performance_logger_per_host_hr::append_caption_entry(
    nv::index_common::String_dict & app_dict)
{
    if(!this->good()){
        return;
    }
    std::ostream & os = *(this->get_os());

    // show the logger name
    os << "# Performance_logger_name: Performance_logger_per_host_hr "
       << (is_system_global_logging() ? "global_system" : "per_host") << "\n";
}

//----------------------------------------------------------------------
void Performance_logger_per_host_hr::append_header_entry()
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
void Performance_logger_per_host_hr::append_per_host_proper_perf_log(
    std::ostream & os, 
    mi::Sint32 frame_num,
    mi::Sint32 host_id)
{
    assert(os.good());
    os << frame_num << "," << host_id << ",";    
}

//----------------------------------------------------------------------
void Performance_logger_per_host_hr::append_perf_log_one_host(
    std::ostream & os, 
    nv::index::IPerformance_values* performance_value,
    mi::Sint32 frame_num,
    mi::Sint32 host_id)
{
    if(m_out_type == HROF_Simple){
        append_perf_log_one_host_simple(os, performance_value);
    }
    else if(m_out_type == HROF_Transfer){
        append_perf_log_one_host_transfer(os, performance_value);
    }
    else{
        // internal error
        assert(false);
    }
}

//----------------------------------------------------------------------
void Performance_logger_per_host_hr::append_perf_log_one_host_simple(
    std::ostream & os, 
    nv::index::IPerformance_values* performance_value)
{
    assert(performance_value != 0);

    const mi::Float32 time_per_frame   = performance_value->get_time("time_complete_frame")    / 1000.f;
    const mi::Float32 rendering_time   = performance_value->get_time("time_total_rendering")   / 1000.f;
    const mi::Float32 compositing_time = performance_value->get_time("time_total_compositing") / 1000.f;
    const mi::Float32 frames_per_sec   = performance_value->get_time("frames_per_second");

    os  << "rendering: " << rendering_time << "\t" << "compositing: " << compositing_time << "\t"
        << "entire_frame: " << time_per_frame << "\t" << "fps: " << frames_per_sec << "\n";
}

//----------------------------------------------------------------------
void Performance_logger_per_host_hr::append_perf_log_one_host_transfer(
    std::ostream & os, 
    nv::index::IPerformance_values* performance_value)
{
    const mi::Float32 compositing_time = performance_value->get_time("time_total_compositing") / 1000.f;
    const mi::Float32 frames_per_sec   = performance_value->get_time("frames_per_second");

    os << "Performance values:\n"
       << "Frames per second: " << frames_per_sec << " fps\n"
       << "Rendering only :   " << performance_value->get_time("time_rendering_only") << " ms\n"
       << "Rendering volume:  " << performance_value->get_time("time_rendering_volume") << " ms\n"
       << "Rendering:         " << performance_value->get_time("time_rendering") << " ms\n"
       << "Total rendering:   " << performance_value->get_time("time_total_rendering") << " ms\n"
       << "GPU download:      " << performance_value->get_time("time_gpu_download") << " ms\n"
       << "GPU upload:        " << performance_value->get_time("time_gpu_upload") << " ms\n"
       << "Total compositing: " << performance_value->get_time("time_total_compositing") << " ms\n"
       << "Final compositing: " << performance_value->get_time("time_total_final_compositing") << " ms\n"
       << "Frame setup:       " << performance_value->get_time("time_frame_setup") << " ms\n"
       << "Frame finish:      " << performance_value->get_time("time_frame_finish") << " ms\n"
       << "Complete frame:    " << performance_value->get_time("time_complete_frame") << " ms\n"
       << "\n"

       << "Uncompressed span images: "
       << (mi::Float32(performance_value->get("size_span_transfer"))/(1024.f*1024.f)) << " MB \t ("
       << ((mi::Float32(performance_value->get("size_span_transfer"))/(1024.f*1024.f)) / compositing_time) << " MB/s)\n"

       << "Compressed span images:   "
       << (mi::Float32(performance_value->get("size_span_transfer_compressed"))/(1024.f*1024.f)) << " MB \t ("
       << ((mi::Float32(performance_value->get("size_span_transfer_compressed"))/(1024.f*1024.f) ) / compositing_time) << " MB/s)\n"

       << "-> Compression ratio:     "
       << ( (mi::Float32(performance_value->get("size_span_transfer_compressed"))/(1024.f*1024.f)) /
            (mi::Float32(performance_value->get("size_span_transfer"))/(1024.f*1024.f) ) ) << "\n"
       << "\n"

       << "Uncompressed tile images: "
       << (mi::Float32(performance_value->get("size_intermediate_image_transfer"))/(1024.f*1024.f)) << " MB \t ("
       << ((mi::Float32(performance_value->get("size_intermediate_image_transfer"))/(1024.f*1024.f)) / compositing_time) << " MB/s)\n"

       << "Compressed tile images:   "
       << (mi::Float32(performance_value->get("size_intermediate_image_transfer_compressed"))/(1024.f*1024.f)) << " MB \t ("
       << ( (mi::Float32(performance_value->get("size_intermediate_image_transfer_compressed"))/(1024.f*1024.f) ) / compositing_time) << " MB/s)\n"

       << "-> Compression ratio:     "
       << ( (mi::Float32(performance_value->get("size_intermediate_image_transfer_compressed"))/(1024.f*1024.f)) /
            (mi::Float32(performance_value->get("size_intermediate_image_transfer"))/(1024.f*1024.f) ) ) << "\n"
       << "\n";

    const mi::Float32 image_compositing_time = performance_value->get_time("time_image_compositing");
    os << "Image compositing time:   " << image_compositing_time << " ms\n"

       << "Composited tile images:   "
       << (mi::Float32(performance_value->get("size_intermediate_image_composited"))/(1024.f*1024.f)) << " MB \t ("
       << ( (mi::Float32(performance_value->get("size_intermediate_image_composited"))/(1024.f*1024.f) ) / compositing_time) << " MB/s)\n\n";
}

//----------------------------------------------------------------------
