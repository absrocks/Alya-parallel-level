/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "volume_data_retrieval_appjob.h"

#include "common/type_conversion_utility.h"

//----------------------------------------------------------------------
Volume_data_retrieval_appjob::Volume_data_retrieval_appjob(
    const mi::neuraylib::Tag & volume_tag,
    const std::vector<mi::math::Bbox<mi::Sint32, 3> >& bbox_vec,
    const mi::base::Handle<const nv::index::IDistributed_data_access_factory>& data_access_factory)
    :
    m_volume_tag(volume_tag),
    m_bounds_vec(bbox_vec),
    m_data_access_factory(data_access_factory)
{
    m_accessed_volume_vec.resize(m_bounds_vec.size());
}

//----------------------------------------------------------------------
void Volume_data_retrieval_appjob::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    const mi::math::Bbox_struct<mi::Uint32, 3> query_bounds = 
        nv::index_common::convert_bbox_type<mi::Uint32, mi::Sint32, 3>(m_bounds_vec[index]);
    INFO_LOG << "Fragment " << index << "(of " << count
             << ") queries volume data contained in: " << query_bounds;

    m_accessed_volume_vec[index] =
        m_data_access_factory->create_regular_volume_data_access(m_volume_tag);
    m_accessed_volume_vec[index]->access(query_bounds, dice_transaction);
}

//----------------------------------------------------------------------
