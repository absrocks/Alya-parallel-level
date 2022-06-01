/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger base class for global/per_host/per_span logger

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_BASE_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_BASE_H

#include "filename_ostream.h"

#include <vector>
#include <string>

#include "common/string_dict.h"

/// Performance logger base class
class Performance_logger_base : public Filename_ostream
{
public:
    /// performance data type
    enum Performance_data_typename_e
    {
        /// Uint64 data
        PDT_Uint64,
        /// Float32 data
        PDT_Float32,
        /// sentinel
        PDT_COUNT
    };

    /// output performance log type
    enum Logging_kind_e {
        /// system total (global) performance
        LK_System_global,
        /// per host performance
        LK_Per_host,
        /// per span performance
        LK_Per_span,
        /// sentinel
        LK_COUNT
    };

public:
    /// constructor
    Performance_logger_base()
        :
        Filename_ostream(),
        m_is_flush(true)
    {
        // empty
    }

    /// destructor
    virtual ~Performance_logger_base()
    {
        // empty
    }

    /// set flush stream each write
    ///
    /// \param[in] is_flush is flush stream each write?
    virtual void set_flush_each_write(bool is_flush)
    {
        m_is_flush = is_flush;
    }
    /// is flush stream each write state
    ///
    /// \return true when stream is flushed each time.
    virtual bool is_flush_each_write() const
    {
        return m_is_flush;
    }

    /// set performance logging typename entries
    ///
    /// \param[in] typeneme_vec performance typename vector.
    virtual void set_performance_logging_entry(const std::vector<std::string> & typename_vec)
    {
        m_perf_typename_entry = typename_vec;
    }

    /// peek current performance logging typename entries. Not this is
    /// peeking the member.
    ///
    /// \return performance typename vector.
    virtual const std::vector<std::string> & peek_performance_logging_entry() const
    {
        return m_perf_typename_entry;
    }

    /// append performance logging header
    ///
    /// \param[in] performance_value IndeX performance value map object
    /// \param[in] app_proj          IndeX application project
    virtual void append_header(
        nv::index_common::String_dict & app_proj) = 0;

    /// append performance log to the current stream
    ///
    /// \param[in] performance_value IndeX performance value map object
    /// \param[in] frame_num         current frame number
    virtual void append_performance_log(
        nv::index::IPerformance_values* performance_value,
        mi::Sint32                      frame_num) = 0;

public:
    /// get all per_host performance value entries from the performance value
    ///
    /// \param[in] performance_value performance value
    /// \return available all per_host/system performance typename vector
    static std::vector<std::string> get_per_host_performance_value_entry_from_performance_value(
        nv::index::IPerformance_values* performance_value);

    /// get all per_span performance value entries from the performance value
    ///
    /// \param[in] performance_value performance value
    /// \return available all span performance typename vector
    static std::vector<std::string> get_per_span_performance_value_entry_from_performance_value(
        nv::index::IPerformance_values* performance_value);

    /// enum type to string
    static std::string get_log_type_string(mi::Uint32 log_type)
    {
        if(log_type == LK_System_global){
            return std::string("global");
        }
        else if(log_type == LK_Per_host){
            return std::string("per_host");
        }
        else if(log_type == LK_Per_span){
            return std::string("per_span");
        }
        ERROR_LOG << "Performance_log_app_state::get_log_type_string: "
                  << "failed to convert Logging_kind_e to string.";
        return std::string("unknown_log_type");            
    }

    /// string to enum type
    static mi::Uint32 get_log_type_enum(const std::string & log_type_name)
    {
        if(log_type_name == "global"){
            return LK_System_global;
        }
        else if(log_type_name == "per_host"){
            return LK_Per_host;
        }
        else if(log_type_name == "per_span"){
            return LK_Per_host;
        }
        ERROR_LOG << "Performance_log_app_state::get_log_type_enum: "
                  << "failed to convert log_type_name [" << log_type_name << "] to enum.";
        return LK_COUNT;
    }


    /// name to data type lookup
    ///
    /// \param[in] entry_name entry name
    /// \return entry name's data type
    inline static Performance_data_typename_e get_data_type(const std::string & entry_name)
    {
        if((entry_name.find("time_") == 0) || (entry_name == "frames_per_second")){
            return PDT_Float32;    // access by get_time()
        }
        return PDT_Uint64;         // access by get()
    }

private:
    /// each time flush the stream?
    /// The derived class should reflect this flag.
    bool m_is_flush;
    /// performance type name entries
    std::vector<std::string> m_perf_typename_entry;

private:
    /// copy constructor. prohibit until proved useful.
    Performance_logger_base(Performance_logger_base const &);
    /// operator=. prohibit until proved useful.
    Performance_logger_base const & operator=(Performance_logger_base const &);
};

//======================================================================
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_BASE_H
