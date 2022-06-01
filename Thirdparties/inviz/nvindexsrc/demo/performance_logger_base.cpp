/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger base class for global/per_host/per_span logger

#include "performance_logger_base.h"

//----------------------------------------------------------------------
std::vector<std::string> Performance_logger_base::get_per_host_performance_value_entry_from_performance_value(
    nv::index::IPerformance_values* performance_value)
{
    assert(performance_value != 0);
    std::vector<std::string> perf_type_vec;
    const mi::Uint32 entry_count = performance_value->get_nb_type_names();
    for(mi::Uint32 i = 0; i < entry_count; ++i){
        const std::string entry_name = performance_value->get_type_name(i);
        perf_type_vec.push_back(entry_name);
    }
    return perf_type_vec;
}

//----------------------------------------------------------------------
std::vector<std::string> Performance_logger_base::get_per_span_performance_value_entry_from_performance_value(
    nv::index::IPerformance_values* performance_value)
{
    assert(performance_value != 0);
    const std::string mn = "Performance_logger_base::get_per_span_performance_value_entry_from_performance_value: ";
    std::vector<std::string> perf_type_vec;

    const mi::Uint32 system_host_id = 0; // for query the number of spans
    const mi::Uint64 nb_spans = performance_value->get("nb_horizontal_spans", system_host_id);
    if(nb_spans == 0){
        ERROR_LOG << mn << "No spans. Fail to get the span statistics key list.";
        return perf_type_vec;
    }

    // nb_spans > 0, we can access to 0.
    mi::Uint32 span_id = 0; 
    mi::base::Handle<nv::index::IPer_span_statistics> per_span_stat(
        performance_value->get_per_span_statistics(span_id));
    if(!(per_span_stat.is_valid_interface())){
        ERROR_LOG << mn << "invalid span_id: " << span_id << "], Fail to get the span statistics key list.";
        return perf_type_vec;
    }

    const mi::Uint32 entry_count = per_span_stat->get_nb_type_names();
    for(mi::Uint32 i = 0; i < entry_count; ++i){
        const std::string entry_name = per_span_stat->get_type_name(i);
        perf_type_vec.push_back(entry_name);
    }
    return perf_type_vec;
}

//----------------------------------------------------------------------

