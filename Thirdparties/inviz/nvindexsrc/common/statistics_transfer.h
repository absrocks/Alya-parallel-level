/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Transfering statistics values over a serializer

#ifndef NVIDIA_INDEX_BIN_COMMON_STATISTICS_TRANSFER_H
#define NVIDIA_INDEX_BIN_COMMON_STATISTICS_TRANSFER_H

#include <mi/neuraylib/dice.h>
#include <mi/neuraylib/iserializer.h>

namespace nv {
namespace index_common {

// Statistics_transfer
class Statistics_transfer : public mi::base::Interface_implement<mi::neuraylib::ISerializable>
{
public:
    struct Statistics_values // The statistics that are supported for transfer.
    {
        mi::Float32 m_fps;
        mi::Float32 m_total_rendering_time;
        mi::Float32 m_total_compositing_time;

        mi::Float32 m_rendering_time_only;
        mi::Float32 m_GPU_uploading_time;
    };
    Statistics_transfer() { }

    const Statistics_values& get_values() const      { return m_statistics_values;   }
    void set_values(const Statistics_values& values) { m_statistics_values = values; }

    // Class' uuid
    typedef mi::base::Uuid_t<0x3984fdfe,0xd38e,0x4f53,0x90,0x3f,0x71,0x0f,0x3b,0xbc,0x80,0x3e> IID;
    virtual mi::base::Uuid get_class_id() const { return IID(); }

    // Implement serialization:
    // CAUTION: PLEASE TAKE CARE THAT THE DE/SERILIZER IS UPDATED EACH TIME A PARAMTER IS ADDED TO THE STATISTICS VALUES!
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const
    {
        serializer->write(&m_statistics_values.m_fps, 5);

    } 
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer)
    {
        deserializer->read(&m_statistics_values.m_fps, 5);
    }

private:
    Statistics_values   m_statistics_values;
};

//----------------------------------------------------------------------

} // namespace index_common
} // namespace nv

#endif // NVIDIA_INDEX_BIN_COMMON_STATISTICS_TRANSFER_H
