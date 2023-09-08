/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief volume data retrieval job for an application

#ifndef GEOSPATIAL_STREAM_VIEWER_VOLUME_DATA_RETRIEVAL_APPJOB_H
#define GEOSPATIAL_STREAM_VIEWER_VOLUME_DATA_RETRIEVAL_APPJOB_H

#include <mi/math/bbox.h>
#include <mi/dice.h>
#include <mi/base/handle.h>
#include <mi/base/interface_declare.h>
#include <nv/index/idistributed_data_locality.h>
#include <nv/index/idistributed_data_access.h>

#include <cassert>
#include <vector>

#include "common/forwarding_logger.h"

//======================================================================

/// Retrieve volume data from cluster nodes for an application.
///
/// This job itself runs only local.
class Volume_data_retrieval_appjob :
        public mi::neuraylib::Fragmented_job<0x481ab536,0xfc47,0x4be0,0xa1,0x94,0xa2,0x35,0x85,0x03,0x04,0x5d>
{
public:
    /// constructor
    ///
    /// \param[in] volume_tag volume data to be retrieved
    /// \param[in] bbox_vec bounding box vector to be retrieved. Each
    /// bounding box should only contains one brick.
    /// \param[in] data_access_factory data access factory. Each
    /// bounding box data is retrieved by jobs created by this
    /// factory.
    Volume_data_retrieval_appjob(
        const mi::neuraylib::Tag&                                                  volume_tag,
        const std::vector<mi::math::Bbox<mi::Sint32, 3> >&                         bbox_vec,
        const mi::base::Handle<const nv::index::IDistributed_data_access_factory>& data_access_factory);

    /// default constructor for serialization only
    Volume_data_retrieval_appjob()
        :
        m_volume_tag(mi::neuraylib::NULL_TAG),
        m_bounds_vec(),
        m_data_access_factory(),
        m_accessed_volume_vec()
    {
        // do nothing, just for serialization
    }

    /// destructor
    virtual ~Volume_data_retrieval_appjob() { /* do nothing as handles are used */ }

    /// Access the queried volume data ...
    const mi::base::Handle<nv::index::IRegular_volume_data_access>& get_seismic_data(
        mi::Uint32 index) const
    {
        return m_accessed_volume_vec[index];
    }

    /// Query job that executes a bunch of sub cubes and renderes their image
    /// into intermediate rendering results
    /// \param dice_transaction The transaction the job is executing in.
    /// \param index The index identifies for which fragment the execute function is called.
    /// \param count The number of fragments in which the job is split.
    /// \param context
    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

private:
    mi::neuraylib::Tag                                                     m_volume_tag;
    std::vector<mi::math::Bbox<mi::Sint32, 3> >                            m_bounds_vec;
    mi::base::Handle<const nv::index::IDistributed_data_access_factory>    m_data_access_factory;

    // Resulting volume data, do not serialize!
    std::vector<mi::base::Handle<nv::index::IRegular_volume_data_access> > m_accessed_volume_vec;
};

//----------------------------------------------------------------------

#endif // GEOSPATIAL_STREAM_VIEWER_VOLUME_DATA_RETRIEVAL_APPJOB_H
