/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief

#ifndef NVIDIA_INDEX_BIN_FRAME_INFO_CALLBACKS_H
#define NVIDIA_INDEX_BIN_FRAME_INFO_CALLBACKS_H

#include <vector>

#include <mi/base/lock.h>

#include <nv/index/iindex.h>

#include "common/forwarding_logger.h"


class Frame_info_callbacks : public mi::base::Interface_implement<nv::index::IFrame_info_callbacks>
{
public:
    struct Dynamic_alloc_event
    {
        mi::Uint32  m_host_id;
        mi::Uint32  m_device_id;
        mi::Size    m_memory_allocation_size;
        mi::Size    m_memory_available;
        mi::Size    m_memory_freed_up;
    };
    struct Device_reset_event
    {
        mi::Uint32  m_host_id;
        mi::Uint32  m_device_id;
    };

public:
    Frame_info_callbacks()
    {
    }
    
    ///
    virtual void set_frame_identifier(const nv::index::IFrame_identifier& frame_id)
    {
        mi::base::Lock::Block block(&m_callback_lock);
        m_frame_id = frame_id;
    }

    const nv::index::IFrame_identifier& get_frame_identifier() const
    {
        mi::base::Lock::Block block(&m_callback_lock);
        return m_frame_id;
    }

    // dynamic memory eviction events /////////////////////////////////////////////////////////////
    virtual void report_dynamic_memory_eviction(
        mi::Uint32  host_id,
        mi::Uint32  device_id,
        mi::Size    memory_allocation_size,
        mi::Size    memory_available,
        mi::Size    memory_freed_up)
    {
        mi::base::Lock::Block block(&m_callback_lock);

        Dynamic_alloc_event e;
        e.m_host_id                 = host_id;
        e.m_device_id               = device_id;
        e.m_memory_allocation_size  = memory_allocation_size;
        e.m_memory_available        = memory_available;
        e.m_memory_freed_up         = memory_freed_up;

        m_dynamic_alloc_events.push_back(e);
    }

    std::vector<Dynamic_alloc_event> get_dynamic_allocation_events() const
    {
        mi::base::Lock::Block block(&m_callback_lock);

        return m_dynamic_alloc_events;
    }

    void clear_dynamic_allocation_events()
    {
        mi::base::Lock::Block block(&m_callback_lock);

        m_dynamic_alloc_events.clear();
    }

    // device reset events ////////////////////////////////////////////////////////////////////////
    virtual void report_device_memory_reset(
        mi::Uint32  host_id,
        mi::Uint32  device_id)
    {
        mi::base::Lock::Block block(&m_callback_lock);

        Device_reset_event e;
        e.m_host_id                 = host_id;
        e.m_device_id               = device_id;

        m_device_reset_events.push_back(e);
    }

    std::vector<Device_reset_event> get_device_reset_events() const
    {
        mi::base::Lock::Block block(&m_callback_lock);

        return m_device_reset_events;
    }

    void clear_device_reset_events()
    {
        mi::base::Lock::Block block(&m_callback_lock);

        m_device_reset_events.clear();
    }

    // workload balancing events //////////////////////////////////////////////////////////////////
    virtual void report_workload_balancing_operations(
        nv::index::IBalancing_operations* balancing_ops)
    {
        mi::base::Lock::Block block(&m_callback_lock);

        m_balancing_ops.push_back(mi::base::Handle<nv::index::IBalancing_operations>(balancing_ops));
    }

    std::vector<mi::base::Handle<nv::index::IBalancing_operations> > get_workload_balancing_operations() const
    {
        mi::base::Lock::Block block(&m_callback_lock);

        return m_balancing_ops;
    }

    void clear_workload_balancing_operations()
    {
        mi::base::Lock::Block block(&m_callback_lock);

        m_balancing_ops.clear();
    }

private:
    // callback lock
    mutable mi::base::Lock              m_callback_lock;
    nv::index::IFrame_identifier        m_frame_id;
    std::vector<Dynamic_alloc_event>    m_dynamic_alloc_events;
    std::vector<Device_reset_event>     m_device_reset_events;

    std::vector<mi::base::Handle<nv::index::IBalancing_operations> > m_balancing_ops;
};


#endif // NVIDIA_INDEX_BIN_FRAME_INFO_CALLBACKS_H
