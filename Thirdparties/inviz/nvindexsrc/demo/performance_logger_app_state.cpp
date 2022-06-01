/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger application state

#include "performance_logger_app_state.h"

//----------------------------------------------------------------------
// lock for the logging state
static mi::base::Lock Performance_logger_logging_state_lock;

//----------------------------------------------------------------------
Performance_logger_app_state::Performance_logger_app_state()
    :
    m_monitoring_on(false)
{
    for(mi::Uint32 i = 0; i < Performance_logger_base::LK_COUNT; ++i){
        m_logging_on[i] = false;
        m_log_filename_vec.push_back("/dev/null");
    }
}
//----------------------------------------------------------------------
Performance_logger_app_state::~Performance_logger_app_state()
{
    // empty
}

//----------------------------------------------------------------------
void Performance_logger_app_state::set_logging_on(mi::Uint32 log_type, bool is_on)
{
    // lock for the state change (turn on/off)
    mi::base::Lock::Block block(&Performance_logger_logging_state_lock);

    if(log_type >= Performance_logger_base::LK_COUNT){
        ERROR_LOG << "Performance_logger_app_state::set_logging_on: invalid log_type ["
                  << log_type << "], ignored.";
        return;
    }
    m_logging_on[log_type] = is_on;
}

//----------------------------------------------------------------------
bool Performance_logger_app_state::is_logging_on(mi::Uint32 log_type) const
{
    if(log_type >= Performance_logger_base::LK_COUNT){
        ERROR_LOG << "Performance_logger_app_state::is_logging_on: invalid log_type ["
                  << log_type << "], ignored.";
        return false;
    }
    return m_logging_on[log_type];
}

//----------------------------------------------------------------------
void Performance_logger_app_state::set_logging_filename(mi::Uint32 log_type, const std::string & fname)
{
    if(log_type >= Performance_logger_base::LK_COUNT){
        ERROR_LOG << "Performance_logger_app_state::set_logging_filename: invalid log_type ["
                  << log_type << "], ignored.";
        return;
    }
    assert(m_log_filename_vec.size() == Performance_logger_base::LK_COUNT);
    m_log_filename_vec[log_type] = fname;
}

//----------------------------------------------------------------------
std::string Performance_logger_app_state::get_logging_filename(mi::Uint32 log_type) const
{
    if(log_type >= Performance_logger_base::LK_COUNT){
        ERROR_LOG << "Performance_logger_app_state::is_logging_on: invalid log_type ["
                  << log_type << "], ignored.";
        return std::string("/dev/null");
    }
    assert(m_log_filename_vec.size() == Performance_logger_base::LK_COUNT);
    return m_log_filename_vec[log_type];
}
//----------------------------------------------------------------------
