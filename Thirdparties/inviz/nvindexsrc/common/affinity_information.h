/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Defines the affinity of spatial areas to machines/GPUs in the cluster.

#ifndef NVIDIA_INDEX_BIN_COMMON_AFFINITY_H
#define NVIDIA_INDEX_BIN_COMMON_AFFINITY_H

#include <vector>

#include <mi/dice.h>
#include <mi/base/interface_declare.h>
#include <mi/math/bbox.h>
#include <mi/neuraylib/iserializer.h>

#include <nv/index/iaffinity_information.h>

namespace nv
{
namespace index_common
{

/// Simple affinity information.
struct Affinity
{
    Affinity()
      : m_host_id(~0u),
        m_gpu_id(~0u)
    {}

    Affinity(
        const mi::math::Bbox<mi::Float32, 3>& bbox,
        mi::Uint32                            host_id,
        mi::Uint32                            gpu_id)
      : m_bbox(bbox),
        m_host_id(host_id),
        m_gpu_id(gpu_id)
    {}

    mi::math::Bbox<mi::Float32, 3> m_bbox;
    mi::Uint32                     m_host_id;
    mi::Uint32                     m_gpu_id;
};

/// Application/user-based affinity information for host assignments base class.
///
class Affinity_information_base :
    public mi::base::Interface_declare<0x3c481110,0xf4d0,0x464e,0x80,0x62,0xd7,0xf4,0xb5,0x31,0x23,0xf4,
        nv::index::IAffinity_information>
{
};

/// Application/user-based affinity information for host assignments class.
class Affinity_information : public mi::base::Interface_implement<Affinity_information_base>
{
public:
    Affinity_information();

    // -------------------------------------------------------------------------------------------
    void add_affinity(
        const mi::math::Bbox<mi::Float32, 3>& bbox,
        mi::Uint32                            host_id,
        mi::Uint32                            gpu_id);

    // -------------------------------------------------------------------------------------------
    virtual bool get_affinity(
        const mi::math::Bbox_struct<mi::Float32, 3>& subregion,
        mi::Uint32&                                  host_id,
        mi::Uint32&                                  gpu_id) const;

    // -------------------------------------------------------------------------------------------
    virtual mi::base::Uuid get_class_id() const { return IID(); }
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    // Affinity information consists of the bounding box that maps to a host and a GPU id.
    std::vector<Affinity>   m_affinity_information;
};


/// Application/user-based affinity information for host assignments base class.
///
class Domain_specific_subdivision_base :
    public mi::base::Interface_declare<0x3c482120,0xf4d0,0x124e,0x77,0x68,0xd1,0xe4,0xb5,0x31,0x23,0xf2,
        nv::index::IDomain_specific_subdivision>
{
};

/// Application/user-based affinity information for host assignments class.
class Domain_specific_subdivision : public mi::base::Interface_implement<Domain_specific_subdivision_base>
{
public:
    Domain_specific_subdivision();

    // -------------------------------------------------------------------------------------------
    void add(
        const mi::math::Bbox<mi::Float32, 3>& bbox,
        mi::Uint32                            host_id = ~0u,
        mi::Uint32                            gpu_id  = ~0u);

    // -------------------------------------------------------------------------------------------
    virtual bool get_affinity(
        const mi::math::Bbox_struct<mi::Float32, 3>& subregion,
        mi::Uint32&                                  host_id,
        mi::Uint32&                                  gpu_id) const;

    // -------------------------------------------------------------------------------------------
    virtual mi::Uint32 get_nb_subregions() const { return m_spatial_subdivision.size(); }
    virtual mi::math::Bbox_struct<mi::Float32, 3> get_subregion(mi::Uint32 index) const;

    // -------------------------------------------------------------------------------------------
    virtual mi::base::Uuid get_class_id() const { return IID(); }
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    // Affinity information consists of the bounding box that maps to a host and a GPU id.
    std::vector<Affinity>   m_spatial_subdivision;
};

} // namespace index_common
} // namespace nv

#endif // NVIDIA_INDEX_BIN_COMMON_AFFINITY_H
