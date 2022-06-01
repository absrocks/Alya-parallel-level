/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger per_span csv output implementation

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_PER_SPAN_CSV_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_PER_SPAN_CSV_H

#include <nv/index/iperformance_values.h>

#include <vector>
#include <string>
#include <algorithm>

#include "common/string_dict.h"
#include "performance_logger_base.h"

//======================================================================
/// performance log file information.
/// - Output kind:   per_span
/// - Output format: csv 
class Performance_logger_per_span_csv : public Performance_logger_base
{
public:
    /// constructor
    Performance_logger_per_span_csv();
    /// destructor
    virtual ~Performance_logger_per_span_csv();

public:
    // Implemented Performance_logger_base

    virtual void append_header(
        nv::index_common::String_dict & app_proj);

    virtual void append_performance_log(
        nv::index::IPerformance_values* performance_value,
        mi::Sint32                      frame_num);

private:
    /// append caption entry.
    ///
    /// Caption: build number, version, record time, title, etc.
    ///
    /// \param[in] performance_value perfromance values
    /// \param[in] frame current number of frame
    void append_caption_entry(nv::index_common::String_dict & app_dict);

    /// append header entries
    ///
    /// Header: Items what performance values we output
    void append_header_entry();

    /// append per_host/system performance logging
    ///
    /// \param[in] os         output stream
    /// \param[in] frame_num  frame number
    /// \param[in] span_id    span id
    /// \param[in] cluster_id cluster id
    void append_per_span_proper_perf_log(std::ostream & os,
                                         mi::Sint32 frame_num,
                                         mi::Uint32 span_id,
                                         mi::Uint32 cluster_id);

    /// append one host logging
    ///
    /// \param[in] os output stream
    /// \param[in] performance_value an IPerformance_values
    /// \param[in] frame_num frame number
    /// \param[in] host_id   host id (when system global, it is 0)
    void append_perf_log_one_span(std::ostream & os, 
                                  nv::index::IPerformance_values* performance_value,
                                  mi::Sint32 frame_num,
                                  mi::Uint32 span_id);

private:
    /// system global logging mode
    bool m_is_system_global_logging_mode;
    /// performance type name entries for per_host
    std::vector<std::string> m_per_host_typename_entry;

private:
    /// copy constructor. prohibit until proved useful.
    Performance_logger_per_span_csv(Performance_logger_per_span_csv const &);
    /// operator=. prohibit until proved useful.
    Performance_logger_per_span_csv const & operator=(Performance_logger_per_span_csv const &);
};

//======================================================================
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_PER_SPAN_CSV_H
