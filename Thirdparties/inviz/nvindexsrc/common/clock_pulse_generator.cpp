/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "clock_pulse_generator.h"

#include "common_utility.h"

namespace nv {
namespace index_common
{

Clock_pulse_generator::Clock_pulse_generator(
    mi::Float64 start_time,
    mi::Float64 end_time,
    mi::Float64 delta) 
  : m_is_running(false),
    m_elapsed_time(0.f),
    m_start_time(0.f),
    m_playback(FORWARD),
    m_repeat(true),
    m_delta(delta),
    m_use_external_setter(false),
    m_t(0.f),
    m_t_start(start_time),
    m_t_end(end_time)
{
}

/// --------------------------------------------------------------------------

mi::Float64 Clock_pulse_generator::get_tick()
{
    if(!m_use_external_setter)
    {
        const mi::Float64 elapsed  = this->elapsed();
        const mi::Float64 duration = m_t_end - m_t_start;
        if(duration==0.)
        {
            return 0.;
        }

        // FORWARD
        if(m_playback == FORWARD)
        {
            m_t = mi::math::fmod(elapsed, duration);
        }
        else // REWIND
        {
            m_t = duration - mi::math::fmod(elapsed, duration);
        }
        
        // REPEAT
        if(!m_repeat && elapsed>=duration)
        {
            this->stop();
        }
    }

    return m_t;
}

/// --------------------------------------------------------------------------

void Clock_pulse_generator::restart()
{
    m_is_running = true;
    m_elapsed_time = 0;
    m_start_time = nv::index_common::get_time();
}

void Clock_pulse_generator::resume()
{
    if(!m_is_running)
    {
        m_is_running = true;
        m_start_time = nv::index_common::get_time();
    }
}

void Clock_pulse_generator::stop()
{
    if(m_is_running)
    {
        m_is_running = false;
        m_elapsed_time += nv::index_common::get_time() - m_start_time;
    }
}

void Clock_pulse_generator::reset(mi::Float64 time)
{
    m_elapsed_time = time;
    if(m_is_running)
    {
        m_start_time = nv::index_common::get_time();
    }
}

mi::Float64 Clock_pulse_generator::elapsed() const
{
    mi::Float64 time = m_elapsed_time;

    if( m_is_running )
    {
        const mi::Float64 passed_time = (nv::index_common::get_time() - m_start_time);
        time += passed_time;
    }
    
    // can never fully trust timers
    // just make sure that elapsed times are at least never negative
    return (time > 0.f) ? time : 0.f;
}

} // namespace index_common
} // namespace nv
