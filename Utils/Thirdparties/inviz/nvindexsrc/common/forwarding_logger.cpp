/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "forwarding_logger.h"

#include <cassert>
#include <cstdio>
#include <iostream>

#include "common_utility.h"

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
/// singleton instance
Forwarding_logger_factory * Forwarding_logger_factory::G_p_forwarding_logger_factory = 0;

//----------------------------------------------------------------------
Forwarding_logger::Forwarding_logger()
    :
    m_os(),
    m_level(mi::base::MESSAGE_SEVERITY_INFO),
    m_forwarding_logger()
{
    // The forwarding logger pointer must be captured by a handle,
    // otherwise leaking some memory.
    if(Forwarding_logger_factory::instance()->is_enabled()){
        mi::base::Handle< mi::base::ILogger > forwarding_logger(
            Forwarding_logger_factory::instance()->get_forwarding_logger());
        m_forwarding_logger.swap(forwarding_logger);
        assert(m_forwarding_logger.is_valid_interface());
    }
    else{
        // no dice available case. Show the local logging
    }

    std::string const hd = Forwarding_logger_factory::instance()->get_message_header();
    if(!(hd.empty())){
        m_os << hd << " ";
    }
}

//----------------------------------------------------------------------
Forwarding_logger::~Forwarding_logger()
{
    m_os << std::endl;

    // Always show the message if it has the level VERBOSE or lower
    if(m_level <= mi::base::MESSAGE_SEVERITY_VERBOSE)
    {
        if(m_forwarding_logger.is_valid_interface())
        {
            m_forwarding_logger->message(m_level, "APP:MAIN", m_os.str().c_str());
        }
        else
        {
            // No forwarding logger available (yet)
            if(Forwarding_logger_factory::instance()->get_fallback_log_severity() >= 
               static_cast<mi::Uint32>(m_level))
            {
                std::string level = level_to_string(m_level);
                while (level.size() < 5)
                {
                    level += ' ';
                }
                std::cout << "        APP    init " << level << ": " << m_os.str();
            }
        }
    }
#ifdef DEBUG
    else
    {
        // Show message with level DEBUG only in debug build
        const std::string& str = m_os.str(); // for cpp check.

        if(m_forwarding_logger.is_valid_interface())
        {
            m_forwarding_logger->message(m_level, "APP:MAIN", str.c_str());
        }
        else
        {
            // No forwarding logger available
            if(Forwarding_logger_factory::instance()->get_fallback_log_severity() >=
               static_cast<mi::Uint32>(m_level))
            {
                std::string level = level_to_string(m_level);
                while (level.size() < 5)
                {
                    level += ' ';
                }
                std::cout << "        APP    init " << level << ": " << m_os.str();
            }
        }
    }
#endif
    // m_forwarding_logger is destructed and ref-count is down here.
    m_forwarding_logger = 0;
}

//----------------------------------------------------------------------
std::ostringstream& Forwarding_logger::get_message(mi::base::Message_severity level)
{
    m_level = level;
    return m_os;
}

//----------------------------------------------------------------------
std::string Forwarding_logger::level_to_string(mi::base::Message_severity level)
{
    static const char* const g_buffer[] = {
        "fatal",
        "error",
        "warn",
        "info",
        "verbose",
        "debug",
        0,
    };

    assert(0 <= level);
    assert(level <= (mi::base::MESSAGE_SEVERITY_DEBUG));

    return g_buffer[level];
}

//----------------------------------------------------------------------
//======================================================================
//----------------------------------------------------------------------
Forwarding_logger_factory::Forwarding_logger_factory()
    :
    m_iindex_if(),
    m_fallback_severity_level(mi::base::MESSAGE_SEVERITY_INFO),
    m_header_str("")
{
    // empty
}

//----------------------------------------------------------------------
Forwarding_logger_factory::~Forwarding_logger_factory()
{
    this->shutdown();
}

//----------------------------------------------------------------------
void Forwarding_logger_factory::initialize(mi::base::Handle<nv::index::IIndex>& iindex_if)
{
    // This must be initialize once before the first use the
    // forwarding logger.
    assert(iindex_if.is_valid_interface());
    m_iindex_if = iindex_if;
}

//----------------------------------------------------------------------
void Forwarding_logger_factory::set_message_header(std::string const & header_str)
{
    m_header_str = header_str;
}

//----------------------------------------------------------------------
std::string Forwarding_logger_factory::get_message_header() const
{
    return m_header_str;
}

//----------------------------------------------------------------------
bool Forwarding_logger_factory::shutdown()
{
    if(!this->is_enabled()){
        return false;
    }

    m_iindex_if = 0;            // unref
    return true;
}

//----------------------------------------------------------------------
bool Forwarding_logger_factory::is_enabled() const
{
    return m_iindex_if.is_valid_interface();
}

//----------------------------------------------------------------------
mi::base::ILogger * Forwarding_logger_factory::get_forwarding_logger() const
{
    if(!m_iindex_if.is_valid_interface()){
        fprintf(stdout, "error: forwarding logger factory has not been initialized. May cause a crash.\n");
        return 0;
    }
    return m_iindex_if->get_forwarding_logger();
}

//----------------------------------------------------------------------
bool Forwarding_logger_factory::set_fallback_log_severity(mi::Uint32 fb_level)
{
    if((fb_level < mi::base::MESSAGE_SEVERITY_FATAL) ||
       (fb_level > mi::base::MESSAGE_SEVERITY_DEBUG))
    {
        return false;
    }

    m_fallback_severity_level = fb_level;
    return true;
}

//----------------------------------------------------------------------
mi::Uint32 Forwarding_logger_factory::get_fallback_log_severity() const
{
    return m_fallback_severity_level;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common

