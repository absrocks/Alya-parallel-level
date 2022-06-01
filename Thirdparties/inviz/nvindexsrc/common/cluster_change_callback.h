/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief  Interface calls for implementing progress callbacks.

#ifndef MI_NVIDIA_INDEX_COMMON_CLUSTER_CHANGE_CALLBACK_H
#define MI_NVIDIA_INDEX_COMMON_CLUSTER_CHANGE_CALLBACK_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/icluster_change_callback.h>

#include "forwarding_logger.h"

namespace nv {
namespace index_common {

/// The class allows implementing user-defined callbacks
/// issued whenever the cluster topology has changed.
///
class Cluster_change_callback : public mi::base::Interface_implement<nv::index::ICluster_change_callback>
{
public:
    Cluster_change_callback() {}
    virtual ~Cluster_change_callback() {}

    /// The updates to the cluster topology will be reported to the callback.
    ///
    /// \param[in] machine_id   The unique id of the machine that left or joined
    ///                         the cluster.
    ///
    /// \param[in] machine_name The name of the machine that left or joined
    ///                         the cluster.
    ///
    /// \param[in] change       The flag that notifies a joining or leaving
    ///                         event.
    ///
    virtual void cluster_change(
        mi::Uint32  machine_id,
        const char* machine_name,
        bool        change)
    {
        if(change)
            VERBOSE_LOG << "User-defined callback reports: host " << machine_name << " (id " << machine_id << ") joined the cluster.";
        else
            VERBOSE_LOG << "User-defined callback reports: host " << machine_name << " (id " << machine_id << ") left the cluster.";
    }
};


}} // namespace nv::index_common

#endif // MI_NVIDIA_INDEX_COMMON_CLUSTER_CHANGE_CALLBACK_H
