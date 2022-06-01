/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief performance logger set holder

#include "performance_logger_set.h"

#include "common/forwarding_logger.h"

// Performance logger implementation
#include "performance_logger_per_host_csv.h"
#include "performance_logger_per_host_hr.h"
#include "performance_logger_per_host_qa.h"
#include "performance_logger_per_span_csv.h"

//----------------------------------------------------------------------
Performance_logger_set::Performance_logger_set()
{
    m_logger_type_name = "";
}

//----------------------------------------------------------------------
Performance_logger_set::~Performance_logger_set()
{
    this->clear();
}

//----------------------------------------------------------------------
bool Performance_logger_set::initialize_logger(const std::string & logger_type_name,
                                               const nv::index_common::String_dict & opt)
{
    if(logger_type_name == "csv"){
        create_csv_performance_logger(opt);
        m_logger_type_name = logger_type_name;
        return true;
    }
    else if(logger_type_name == "human_readable"){
        create_human_readable_performance_logger(opt);
        m_logger_type_name = logger_type_name;
        return true;
    }
    else if(logger_type_name == "qa"){
        create_qa_performance_logger(opt);
        m_logger_type_name = "qa";
        return true;
    }
    ERROR_LOG << "No such performance logger type [" << logger_type_name
              << "], use default performance logger.";
    create_csv_performance_logger(opt);
    m_logger_type_name = "csv";

    return false;
}

//----------------------------------------------------------------------
std::string Performance_logger_set::get_logger_type_name() const
{
    return m_logger_type_name;
}

//----------------------------------------------------------------------
void Performance_logger_set::clear()
{
    if(!(m_perf_logger_ref_vec.empty())){
        for(mi::Sint32 i = 0; i < Logger_count; ++i){
            Performance_logger_base * p_logbase = m_perf_logger_ref_vec.at(i);
            if(p_logbase != 0){
                delete p_logbase;
                m_perf_logger_ref_vec.at(i) = 0;            
            }
        }
        m_perf_logger_ref_vec.clear();
    }
    m_logger_type_name = "";
}

//----------------------------------------------------------------------
bool Performance_logger_set::is_valid() const
{
    bool is_valid = true;
    const mi::Size logger_count = m_perf_logger_ref_vec.size();
    for(mi::Size i = 0; i < logger_count; ++i){
        if(m_perf_logger_ref_vec.at(i) == 0){
            ERROR_LOG << "Performance logger [" << i << "] is invalid.";
            is_valid = false;
            break;
        }
    }
    return is_valid;
}

//----------------------------------------------------------------------
Performance_logger_base * Performance_logger_set::peek_logger(mi::Uint32 logger_type)
{
    if(logger_type >= Performance_logger_base::LK_COUNT){
        ERROR_LOG << "invalid logger_type [" << logger_type << "]";
        return 0;
    }
    if(m_perf_logger_ref_vec.empty()){
        return 0;
    }

    return m_perf_logger_ref_vec.at(logger_type);
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
void Performance_logger_set::create_csv_performance_logger(
    const nv::index_common::String_dict & opt)
{
    this->clear();

    // added pointer area (but 0)
    m_perf_logger_ref_vec.resize(Logger_count, 0);


    // system global performance logger
    {
        Performance_logger_per_host_csv * system_logger = 
            new Performance_logger_per_host_csv();
        const bool is_system_global = true;
        system_logger->set_system_global_logging(is_system_global);
        m_perf_logger_ref_vec.at(Performance_logger_base::LK_System_global) = system_logger;
    }

    // per-host performance logger
    {
        Performance_logger_per_host_csv * per_host_logger = 
            new Performance_logger_per_host_csv();
        const bool is_system_global = false;
        per_host_logger->set_system_global_logging(is_system_global);
        m_perf_logger_ref_vec.at(Performance_logger_base::LK_Per_host) = per_host_logger;
    }

    // per-span performance logger
    {
        Performance_logger_per_span_csv * per_span_logger =
            new Performance_logger_per_span_csv();
        m_perf_logger_ref_vec.at(Performance_logger_base::LK_Per_span) = per_span_logger;
    }
}

//----------------------------------------------------------------------
void Performance_logger_set::create_human_readable_performance_logger(
    const nv::index_common::String_dict & opt)
{
    this->clear();

    // added pointer area (but 0)
    m_perf_logger_ref_vec.resize(Logger_count, 0);

    // system global performance logger
    {
        Performance_logger_per_host_hr * system_logger = 
            new Performance_logger_per_host_hr();
        const bool is_system_global = true;
        system_logger->set_system_global_logging(is_system_global);
        system_logger->set_output_format(opt);
        m_perf_logger_ref_vec.at(Performance_logger_base::LK_System_global) = system_logger;
    }

    // per-host performance logger
    {
        Performance_logger_per_host_hr * per_host_logger = 
            new Performance_logger_per_host_hr();
        const bool is_system_global = false;
        per_host_logger->set_system_global_logging(is_system_global);
        per_host_logger->set_output_format(opt);
        m_perf_logger_ref_vec.at(Performance_logger_base::LK_Per_host) = per_host_logger;
    }

    // per-span performance logger
    {
        Performance_logger_per_span_csv * per_span_logger =
            new Performance_logger_per_span_csv();
        m_perf_logger_ref_vec.at(Performance_logger_base::LK_Per_span) = per_span_logger;
    }
}

//----------------------------------------------------------------------
void Performance_logger_set::create_qa_performance_logger(
    const nv::index_common::String_dict & opt)
{
    this->clear();

    // added pointer area (but 0)
    m_perf_logger_ref_vec.resize(Logger_count, 0);

    // system global performance logger
    {
        Performance_logger_per_host_qa * system_logger = 
            new Performance_logger_per_host_qa();
        const bool is_system_global = true;
        system_logger->set_system_global_logging(is_system_global);
        m_perf_logger_ref_vec.at(Performance_logger_base::LK_System_global) = system_logger;
    }

    // per-host performance logger
    {
        Performance_logger_per_host_hr * per_host_logger = 
            new Performance_logger_per_host_hr();
        const bool is_system_global = false;
        per_host_logger->set_system_global_logging(is_system_global);
        per_host_logger->set_output_format(opt);
        m_perf_logger_ref_vec.at(Performance_logger_base::LK_Per_host) = per_host_logger;
    }

    // per-span performance logger
    {
        Performance_logger_per_span_csv * per_span_logger =
            new Performance_logger_per_span_csv();
        m_perf_logger_ref_vec.at(Performance_logger_base::LK_Per_span) = per_span_logger;
    }
}

//======================================================================
