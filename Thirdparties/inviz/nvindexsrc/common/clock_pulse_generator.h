/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Defines the some basic clock pulse generator for time-dependent data vis.

#ifndef NVIDIA_INDEX_BIN_COMMON_CLOCK_PULSE_GENERATOR_H
#define NVIDIA_INDEX_BIN_COMMON_CLOCK_PULSE_GENERATOR_H

#include <mi/dice.h>
#include <mi/base/interface_declare.h>
#include <mi/base/types.h>

#include <nv/index/itime_mapping.h>

namespace nv
{
namespace index_common
{

/// Application/user-based affinity information for host assignments base class.
///
class Clock_pulse_generator_base :
    public mi::base::Interface_declare<0x203a51f4,0x9e0a,0x4456,0xb0,0xee,0x61,0xab,0xf9,0x54,0xfd,0x20,
        nv::index::IClock_pulse_generator>
{
};


/// The interface class \c IClock_pulse_generator generates the time for playing animations.
///
/// The user has the ability to specify and implement a custom clock pulse generator to steer the
/// time-dependent visualization of large-scale data.
///
/// Application/user-based affinity information for host assignments class.
class Clock_pulse_generator : public mi::base::Interface_implement<Clock_pulse_generator_base>
{
public:
    enum Playback
    {
        FORWARD = 0,
        REWIND  = 1
    };
	Clock_pulse_generator(mi::Float64 start = 0., mi::Float64 end = 0., mi::Float64 delta = 0.);

	/// --------------------------------------------------------------------------
    void restart();
    void resume();
    void start() { resume(); }
    void stop();
    void reset(mi::Float64 time = 0.);
    mi::Float64 elapsed() const;
    
    /// --------------------------------------------------------------------------
    bool is_running() const { return m_is_running; }

    void set_playback(Playback playback) { m_playback = playback; }
    Playback get_playback() const        { return m_playback;     }

    void set_repeat(bool repeat)         { m_repeat = repeat;     }
    bool get_repeat() const              { return m_repeat;       }

    void set_delta(mi::Float64 delta)    { m_delta = delta;       }
    mi::Float64 get_delta() const        { return m_delta;        }

    void set_tick(mi::Float64 tock)      { m_t = tock;            }
    void set_start(mi::Float64 start)    { m_t_start = start;     }
    void set_end(mi::Float64 end)        { m_t_end = end;         }

    void use_external_setter(bool s)     { m_use_external_setter = s; }

    /// --------------------------------------------------------------------------
    /// Implementaion of IClock_pulse_generator.
    virtual mi::Float64 get_tick();
    virtual mi::Float64 get_start() const { return m_t_start; }
    virtual mi::Float64 get_end() const   { return m_t_end;   }

private:
    bool        	m_is_running;       ///< true if running
    mi::Float64     m_elapsed_time;     ///< cumulative time
    mi::Float64     m_start_time;       ///< start time

    Playback        m_playback;
    bool            m_repeat;
    mi::Float64     m_delta;

    bool            m_use_external_setter;

    mi::Float64     m_t;
	mi::Float64     m_t_start;
	mi::Float64     m_t_end;
};


} // namespace index_common
} // namespace nv

#endif // NVIDIA_INDEX_BIN_COMMON_CLOCK_PULSE_GENERATOR_H