/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "heightfield_data_retrieval_appjob.h"

#include "common/type_conversion_utility.h"

//----------------------------------------------------------------------
Heightfield_data_retrieval_appjob::Heightfield_data_retrieval_appjob(
    const mi::neuraylib::Tag&                                                   heightfield_tag,
    const std::vector<mi::math::Bbox<mi::Sint32, 2> >&                          bounds,
    const mi::base::Handle<const nv::index::IDistributed_data_access_factory>&  data_access_factory)
    :
    m_heightfield_tag(heightfield_tag),
    m_bounds(bounds),
    m_data_access_factory(data_access_factory)
{
    m_accessed_heightfield.resize(m_bounds.size());
}

//----------------------------------------------------------------------
void Heightfield_data_retrieval_appjob::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::Size nb_bounds = m_bounds.size();
    for (mi::Size i = 0; i < nb_bounds; ++i)
    {
        const mi::math::Bbox<mi::Uint32, 2> query_bounds = 
            nv::index_common::convert_bbox_type<mi::Uint32, mi::Sint32, 2>(m_bounds[i]);
        INFO_LOG << "Fragment " << index << " (of " << count << ") queries heightfield data contained in: "
                 << query_bounds;
        m_accessed_heightfield[i] = m_data_access_factory->create_regular_heightfield_data_access(m_heightfield_tag);
        m_accessed_heightfield[i]->access(query_bounds, dice_transaction);
    }
}

//----------------------------------------------------------------------
