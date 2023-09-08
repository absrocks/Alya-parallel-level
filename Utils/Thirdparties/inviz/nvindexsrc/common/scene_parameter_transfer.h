/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Transfering scene parameter over a serializer

#ifndef NVIDIA_INDEX_BIN_COMMON_SCENE_PARAMETER_TRANSFER_H
#define NVIDIA_INDEX_BIN_COMMON_SCENE_PARAMETER_TRANSFER_H

#include <mi/neuraylib/dice.h>
#include <mi/neuraylib/iserializer.h>

namespace nv {
namespace index_common {

// Statistics_transfer
class Scene_parameter_transfer : public mi::base::Interface_implement<mi::neuraylib::ISerializable>
{
public:
    struct Animation_parameter // The camera parameter that are supported for transfer.
    {
        mi::Float32 m_t_start;
        mi::Float32 m_t_end;
        mi::Float32 m_t_current;
    };

    Scene_parameter_transfer() { }

    const Animation_parameter& get_animation() const         { return m_animation_parameter;       }
    void set_animation(const Animation_parameter& parameter) { m_animation_parameter = parameter;  }

    // Class' uuid
    typedef mi::base::Uuid_t<0x3cc8c59d,0x89f8,0x4be8,0xbf,0xea,0x17,0x77,0x3e,0xb4,0xac,0x70> IID;
    virtual mi::base::Uuid get_class_id() const { return IID(); }

    // Implement serialization:
    // CAUTION: PLEASE TAKE CARE THAT THE DE/SERILIZER IS UPDATED EACH TIME A PARAMTER IS ADDED TO THE STATISTICS VALUES!
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const
    {
        serializer->write(&m_animation_parameter.m_t_start,    3);
    }

    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer)
    {
        deserializer->read(&m_animation_parameter.m_t_start,   3);
    }

private:
    Animation_parameter     m_animation_parameter;
};

//----------------------------------------------------------------------

} // namespace index_common
} // namespace nv

#endif // NVIDIA_INDEX_BIN_COMMON_SCENE_PARAMETER_TRANSFER_H
