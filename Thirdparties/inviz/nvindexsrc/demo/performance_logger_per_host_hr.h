/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger global/per_host human readable output implementation

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_PER_HOST_HR_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_PER_HOST_HR_H

#include <nv/index/iperformance_values.h>

#include <vector>
#include <string>
#include <algorithm>

#include "common/string_dict.h"
#include "performance_logger_base.h"

//======================================================================
/// performance log file information.
/// - Output kind:   per_host and system global (default: system global)
/// - Output format: human readable format
///
/// Project entries:
/// - app::performance::global::human_readable_type   ... {simple, transfer}
/// - app::performance::per_host::human_readable_type ... {simple, transfer}
class Performance_logger_per_host_hr : public Performance_logger_base
{
public:
    /// constructor
    Performance_logger_per_host_hr();
    /// destructor
    virtual ~Performance_logger_per_host_hr();

    /// set system global logging mode
    /// 
    /// \param[in] is_system_global system global logging when true
    void set_system_global_logging(bool is_system_global);

    /// is system global logging mode
    /// 
    /// \return true when system global logging mode
    bool is_system_global_logging() const;

    /// initialize format depends on the current option state
    ///
    /// \param[in] opt an extra option for this logger
    void set_output_format(const nv::index_common::String_dict & opt);

public:
    // Implemented Performance_logger_base

    virtual void append_header(
        nv::index_common::String_dict & app_proj);

    virtual void append_performance_log(
        nv::index::IPerformance_values* performance_value,
        mi::Sint32                      frame_num);

private:
    /// human readable output format type
    enum Human_readable_output_format_e {
        /// simple format
        HROF_Simple, 
        /// transfer data format
        HROF_Transfer, 
        /// sentinel
        HRDF_Count
    };

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
    /// \param[in] os output stream
    /// \param[in] frame_num frame number
    /// \param[in] host_id   host id (when system global, it is 0)
    void append_per_host_proper_perf_log(std::ostream & os, mi::Sint32 frame_num, mi::Sint32 host_id);

    /// append one host logging
    ///
    /// \param[in] os output stream
    /// \param[in] performance_value an IPerformance_values
    /// \param[in] frame_num frame number
    /// \param[in] host_id   host id (when system global, it is 0)
    void append_perf_log_one_host(std::ostream & os, 
                                  nv::index::IPerformance_values* performance_value,
                                  mi::Sint32 frame_num,
                                  mi::Sint32 host_id);

    /// append one host logging: simple
    ///
    /// \param[in] os output stream
    /// \param[in] performance_value an IPerformance_values
    void append_perf_log_one_host_simple(std::ostream & os, 
                                         nv::index::IPerformance_values* performance_value);

    /// append one host logging: transfer
    ///
    /// \param[in] os output stream
    /// \param[in] performance_value an IPerformance_values
    void append_perf_log_one_host_transfer(std::ostream & os, 
                                           nv::index::IPerformance_values* performance_value);

private:
    /// system global logging mode
    bool m_is_system_global_logging_mode;
    /// performance type name entries for per_host
    std::vector<std::string> m_per_host_typename_entry;
    /// output type
    mi::Uint32 m_out_type;

private:
    /// copy constructor. prohibit until proved useful.
    Performance_logger_per_host_hr(Performance_logger_per_host_hr const &);
    /// operator=. prohibit until proved useful.
    Performance_logger_per_host_hr const & operator=(Performance_logger_per_host_hr const &);
};

//======================================================================
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_PER_HOST_HR_H
