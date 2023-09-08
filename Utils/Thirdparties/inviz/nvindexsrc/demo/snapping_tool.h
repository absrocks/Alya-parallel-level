/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief snapping tool example based on a picking position
#ifndef NVIDIA_INDEX_EXTERNAL_SNAPPING_TOOL_H
#define NVIDIA_INDEX_EXTERNAL_SNAPPING_TOOL_H

#include <mi/dice.h>
#include <mi/base/interface_implement.h>
#include <mi/base/handle.h>

#include <nv/index/isnapping_algorithm.h>
#include <nv/index/iregular_volume.h>

#include "common/forwarding_logger.h"

#include "utilities.h"


/// A simple user-defined snapping tool that illustrates the use of the snapping algorithm
/// for manual picking
class Snapping_tool : public mi::base::Interface_implement<nv::index::ISnapping_algorithm>
{
public:
    /// constructor.
    Snapping_tool(mi::Uint32 k_range_min = 10,
                  mi::Uint32 k_range_max = 10)
        : m_k_range_min(k_range_min),
          m_k_range_max(k_range_max)
    {
        // empty
    }

    /// destructor
    virtual ~Snapping_tool()
    {
        // empty
    }

    /// Adjusting a manual pick applied to a heightfield
    /// \param[in] pick The original manual pick position in IJK space resulting from,
    ///                 for instance, the ray-intersection with a slice.
    /// \param[in] volume_data_access_techniques The volume data
    ///                 access provides means to access the volume amplitude values
    ///                 required for snapping the pick accorging to the volume
    ///                 amplitude values below and above the intersection.
    /// \param[in] dice_transaction The dice transaction.
    /// \param[in] result The adjusted pick position in IJK space of the slice.
    /// \return Tells if the returned pick position is valid.
    virtual bool adjust_pick(
        const mi::math::Vector_struct<mi::Float32, 3>&  ijk_pick_position,
        nv::index::IRegular_volume_data_access*         volume_data_access,
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::math::Vector_struct<mi::Float32, 3>&        result) const;

private:
    mi::Uint32 m_k_range_min;
    mi::Uint32 m_k_range_max;
};


#endif // NVIDIA_INDEX_EXTERNAL_SNAPPING_TOOL_H
