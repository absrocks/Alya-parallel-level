/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Discrete normal encoder
#ifndef NVIDIA_INDEX_BIN_COMMON_NORMAL_ENCODER_H
#define NVIDIA_INDEX_BIN_COMMON_NORMAL_ENCODER_H

#include <mi/dice.h>

namespace nv {
namespace index_common {

/// Discrete normal encoder
/// The patch normal representation will be upgraded to Vector<Float32, 3>.
/// \deprecated
class Normal_encoder
{
public:
    /// get the instance
    /// \return normal encoder singleton instance
    static Normal_encoder* instance()
    {
        if (G_p_normal_encoder == 0)
        {
            G_p_normal_encoder = new Normal_encoder();
        }
        return G_p_normal_encoder;
    }

    /// delete the singleton. (For unit test purpose. I don't want to
    /// confuse a memory checker.)
    static void delete_instance()
    {
        if (G_p_normal_encoder != 0)
        {
            delete G_p_normal_encoder;
            G_p_normal_encoder = 0;
        }
    }

public:
    mi::math::Vector<mi::Float32, 3> get_normal(mi::Uint32 v) const
    {
        const mi::Uint32 theta = v >> 16;
        const mi::Uint32 phi   = v & 0xffff;

        const mi::math::Vector<mi::Float32, 3> normal(
            sinTheta[theta] * cosPhi[phi],
            sinTheta[theta] * sinPhi[phi],
            cosTheta[theta]);

        return normal;
    }

    /// convert normal vector to encoded normal
    static mi::Uint32 get_encoded_normal(const mi::math::Vector<mi::Float32, 3> & normal_vec)
    {
        const mi::Uint32 theta   = static_cast<mi::Uint32>(4095.0f * acos(normal_vec.z) / M_PI);
        const mi::Uint32 phi     = static_cast<mi::Uint32>(4095.0f * (atan2(normal_vec.y, normal_vec.x) + M_PI) / (2 * M_PI));
        const mi::Uint32 enc_val = theta << 16 | phi;

        return enc_val;
    }

private:
    /// constructor
    Normal_encoder()
    {
        static mi::Float64 pi = 3.14159265358979323846;
        for (mi::Sint32 i = 0; i < 4096; i++)
        {
            cosTheta[i] = static_cast<mi::Float32>(cos(i * pi / 4095.0));
            cosPhi[i]   = static_cast<mi::Float32>(cos(-pi + i * (2.0 * pi) / 4095.0));
            sinTheta[i] = static_cast<mi::Float32>(sin(i * pi / 4095.0));
            sinPhi[i]   = static_cast<mi::Float32>(sin(-pi + i * (2.0 * pi) / 4095.0));
        }
    }

private:
    // singleton instance
    static Normal_encoder * G_p_normal_encoder;

private:
    mi::Float32   sinPhi[4096];
    mi::Float32   cosPhi[4096];
    mi::Float32   sinTheta[4096];
    mi::Float32   cosTheta[4096];
};

} // namespace index_common
} // namespace nv
#endif // #ifndef NVIDIA_INDEX_BIN_COMMON_NORMAL_ENCODER_H
