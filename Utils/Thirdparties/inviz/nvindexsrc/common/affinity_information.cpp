/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Defines the affinity of spatial areas to machines/GPUs in the cluster.

#include "affinity_information.h"

namespace nv
{
namespace index_common
{

// -------------------------------------------------------------------------------------------
Affinity_information::Affinity_information()
{
}

// -------------------------------------------------------------------------------------------
void Affinity_information::add_affinity(
    const mi::math::Bbox<mi::Float32, 3>& bbox,
    mi::Uint32                           host_id,
    mi::Uint32                           gpu_id)
{
    const Affinity affinity(bbox, host_id, gpu_id);
    m_affinity_information.push_back(affinity);
}

// -------------------------------------------------------------------------------------------
bool Affinity_information::get_affinity(
    const mi::math::Bbox_struct<mi::Float32, 3>& subregion_st,
    mi::Uint32&                                 host_id,
    mi::Uint32&                                 gpu_id) const
{
    const mi::math::Bbox<mi::Float32, 3> subregion(subregion_st);

    const mi::Uint32 nb_elements = m_affinity_information.size();
    for (mi::Size i=0; i < nb_elements; ++i)
    {
        const Affinity& affinity = m_affinity_information[i];
        if(affinity.m_bbox.contains(subregion.min) && affinity.m_bbox.contains(subregion.max) )
        {
            host_id = affinity.m_host_id;
            gpu_id  = affinity.m_gpu_id;
            return true;
        }
    }

    return false;
}

// -------------------------------------------------------------------------------------------
void Affinity_information::serialize(mi::neuraylib::ISerializer* serializer) const
{
    const mi::Uint32 nb_elements = m_affinity_information.size();
    serializer->write(&nb_elements, 1);
    for (mi::Size i=0; i < nb_elements; ++i)
    {
        const Affinity& affinity = m_affinity_information[i];
        serializer->write(&affinity.m_bbox.min.x, 6);
        serializer->write(&affinity.m_host_id,    1);
        serializer->write(&affinity.m_gpu_id,     1);
    }
}
void Affinity_information::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    for (mi::Size i=0; i < nb_elements; ++i)
    {
        Affinity affinity;
        deserializer->read(&affinity.m_bbox.min.x, 6);
        deserializer->read(&affinity.m_host_id,    1);
        deserializer->read(&affinity.m_gpu_id,     1);
        m_affinity_information.push_back(affinity);
    }
}

// ==================================================================================================

Domain_specific_subdivision::Domain_specific_subdivision()
{
}

// -------------------------------------------------------------------------------------------
void Domain_specific_subdivision::add(
    const mi::math::Bbox<mi::Float32, 3>& bbox,
    mi::Uint32                            host_id,
    mi::Uint32                            gpu_id)
{
    const Affinity affinity(bbox, host_id, gpu_id);
    m_spatial_subdivision.push_back(affinity);
}

// -------------------------------------------------------------------------------------------
bool Domain_specific_subdivision::get_affinity(
    const mi::math::Bbox_struct<mi::Float32, 3>& subregion_st,
    mi::Uint32&                                  host_id,
    mi::Uint32&                                  gpu_id) const
{
    const mi::math::Bbox<mi::Float32, 3> subregion(subregion_st);

    const mi::Uint32 nb_elements = m_spatial_subdivision.size();
    for (mi::Size i=0; i < nb_elements; ++i)
    {
        const Affinity& affinity = m_spatial_subdivision[i];
        if(affinity.m_bbox.contains(subregion.min) && affinity.m_bbox.contains(subregion.max) )
        {
            host_id = affinity.m_host_id;
            gpu_id  = affinity.m_gpu_id;
            return true;
        }
    }

    return false;
}

// -------------------------------------------------------------------------------------------
void Domain_specific_subdivision::serialize(mi::neuraylib::ISerializer* serializer) const
{
    const mi::Uint32 nb_elements = m_spatial_subdivision.size();
    serializer->write(&nb_elements, 1);
    for (mi::Size i=0; i < nb_elements; ++i)
    {
        const Affinity& affinity = m_spatial_subdivision[i];
        serializer->write(&affinity.m_bbox.min.x, 6);
        serializer->write(&affinity.m_host_id,    1);
        serializer->write(&affinity.m_gpu_id,     1);
    }
}
void Domain_specific_subdivision::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    for (mi::Size i=0; i < nb_elements; ++i)
    {
        Affinity affinity;
        deserializer->read(&affinity.m_bbox.min.x, 6);
        deserializer->read(&affinity.m_host_id,    1);
        deserializer->read(&affinity.m_gpu_id,     1);
        m_spatial_subdivision.push_back(affinity);
    }
}

mi::math::Bbox_struct<mi::Float32, 3> Domain_specific_subdivision::get_subregion(mi::Uint32 index) const
{
    const Affinity& spatial_subdivision = m_spatial_subdivision[index];
    return spatial_subdivision.m_bbox;
}

} // namespace index_common
} // namespace nv
