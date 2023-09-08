/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger set holder

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_SET_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_SET_H

#include <mi/dice.h>

#include <vector>
#include <string>

#include "common/string_dict.h"

class Performance_logger_base;

//======================================================================
/// Performance logger set managers.
///
/// This contains performance logger set. This depends on:
/// - What performance logger is used (for csv output, for QA output).
/// - Request the type of the logger (for system, per_host, per_span).
class Performance_logger_set
{
public:
    /// default constructor
    Performance_logger_set();

    /// destructor
    ~Performance_logger_set();

    /// initialize loggers.
    ///
    /// \param[in] logger_type_name logger type ("csv", "human_readable")
    /// \param[in] opt    extra option (depends on the implementation)
    /// \return true when the type name has been accepted
    bool initialize_logger(const std::string & logger_type_name,
                           const nv::index_common::String_dict & opt);

    /// get current logger type name
    /// \return logger type name (or "" when undefined)
    std::string get_logger_type_name() const;

    /// clear loggers
    void clear();

    /// the logger set state is valid?
    bool is_valid() const;

    /// access to the performance logger
    Performance_logger_base * peek_logger(mi::Uint32 logger_type);

private:
    /// create csv format performance logger set.
    void create_csv_performance_logger(const nv::index_common::String_dict & opt);
    /// create human readable performance logger set.
    void create_human_readable_performance_logger(const nv::index_common::String_dict & opt);
    /// create QA performance logger set.
    void create_qa_performance_logger(const nv::index_common::String_dict & opt);

private:
    /// system, per-host, per_span.
    const static mi::Sint32 Logger_count = 3;
    /// the real logger containers for three types (system, per-host, per_span)
    std::vector<Performance_logger_base *> m_perf_logger_ref_vec;
    /// current logger type name
    std::string m_logger_type_name;

private:
    /// copy constructor. prohibit until proved useful.
    Performance_logger_set(Performance_logger_set const &);
    /// operator=. prohibit until proved useful.
    Performance_logger_set const & operator=(Performance_logger_set const &);
};

//======================================================================
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_SET_H
