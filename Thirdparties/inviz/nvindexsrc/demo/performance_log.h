/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief NVIDIA IndeX demo viewer performance logging functions

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOG_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOG_H

#include <mi/dice.h>

#include <nv/index/iperformance_values.h>
#include <nv/index/iconfig_settings.h>

#include "common/string_dict.h"

class Nvindex_rendering_context;

//----------------------------------------------------------------------
/// Update recent n statistics
///
/// \param[in] performance_values a IPerformance_values
void update_recent_n_statistics(
    nv::index::IPerformance_values * performance_values);

//----------------------------------------------------------------------
/// Set monitor performance value.
/// Please use this to synchronize IndeX state and application state.
///
/// \param[in] edit_config_settings editable valid config settings
/// \param[in] is_monitor_on true when performance monitoring is on
///            (actually output is depends on the performance logger
///            instance)
void set_monitor_performance_values(
    mi::base::Handle<nv::index::IConfig_settings> & edit_config_settings,
    bool is_monitor_on);

//----------------------------------------------------------------------
/// set performance monitoring state by project settings
///
/// \param[in] irc              IndeX rendering context
/// \param[in] opt              option for performance monitoring setting
/// \param[in] dice_transaction dice transaction
void set_performance_monitoring_by_project(
    Nvindex_rendering_context& irc,
    nv::index_common::String_dict & opt,
    mi::neuraylib::IDice_transaction * dice_transaction);

//----------------------------------------------------------------------
/// monitor all performance data depends on the performance logging state
///
/// \param[in] performance_values an IPerformance_value
/// \param[in] frame_num          current frame number
void log_performance_information(
    nv::index::IPerformance_values * performance_value,
    mi::Sint32 frame_num);

//----------------------------------------------------------------------
/// set up performance logger by a command line string
///
/// Command example:
///  - "perf_logger csv"
///     set cvs performance logger
///  - "perf_logger human_readable simple"
///     set human_readable performance logger with option "simple"
///  - "perf_logger human_readable transfer"
///     set human_readable performance logger with option "transfer"
///
/// \param[in]  command a command string. The first token should be perf_logger.
/// \param[out] ret_mes return message
/// \return true when success
bool setup_performance_logger_by_command_str(
    const std::string & command,
    std::string & ret_mes);

//----------------------------------------------------------------------
/// log performance via command line GUI
///     for performance regression test
bool log_performance_by_command_str(
    const std::string & command,
    std::string & ret_mes);

//----------------------------------------------------------------------
/// turn off all loggers
void turn_off_logging();

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_PERFORMANCE_LOG_H
