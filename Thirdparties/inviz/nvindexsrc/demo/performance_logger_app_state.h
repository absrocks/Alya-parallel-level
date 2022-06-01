/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger application state

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_APP_STATE_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_APP_STATE_H

#include <mi/dice.h>

#include <string>

#include "performance_logger_base.h"

//======================================================================
/// Performance logging application state.
///
/// This keeps application state what kind of performance information
/// to be logged.
class Performance_logger_app_state
{
public:
    /// default constructor
    Performance_logger_app_state();
    /// destructor
    ~Performance_logger_app_state();

public:
    /// set monitoring on/off
    ///
    /// Monitoring is IndeX performance monitoring. Performance value
    /// output or not is independent from monitoring, but without
    /// monitoring, output (logging) has no meaning. This is a local
    /// copy for logging, the real value is in the IndeX's
    /// IConfig_setting.
    ///
    /// \param[in] is_on is monitoring is currently on?
    void set_monitoring_on(bool is_on)
    {
        m_monitoring_on = is_on;
    }

    /// query monitoring state on/off
    ///
    /// \return true when monitoring is on
    bool is_monitoring_on() const
    {
        return m_monitoring_on;
    }

    /// set logging state on/off
    ///
    /// \param[in] log_type log type, which logging is on?
    /// \param[in] is_on    on when true
    void set_logging_on(mi::Uint32 log_type, bool is_on);

    /// query logging state on/off
    ///
    /// \param[in] log_type log type, which logging is on?
    /// \return true when logging is on
    bool is_logging_on(mi::Uint32 log_type) const;

    /// set logging filename
    ///
    /// \param[in] log_type log type, which logging is on?
    /// \param[in] fname    log file filename
    void set_logging_filename(mi::Uint32 log_type, const std::string & fname);

    /// query logging filename
    ///
    /// \param[in] log_type log type, which logging is on?
    /// \return filename of log_type
    std::string get_logging_filename(mi::Uint32 log_type) const;

private:
    /// IndeX performance monitoring state cache.
    bool m_monitoring_on;
    /// output logging state for system global, per_host, per_span.
    bool m_logging_on[Performance_logger_base::LK_COUNT];
    /// output logging filename.
    std::vector<std::string> m_log_filename_vec;

private:
    /// copy constructor. prohibit until proved useful.
    Performance_logger_app_state(Performance_logger_app_state const &);
    /// operator=. prohibit until proved useful.
    Performance_logger_app_state const & operator=(Performance_logger_app_state const &);
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOGGER_APP_STATE_H
