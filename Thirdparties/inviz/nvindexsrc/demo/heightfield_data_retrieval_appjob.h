/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Retrieve heightfield data from cluster nodes

#ifndef GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_DATA_RETRIEVAL_APPJOB_H
#define GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_DATA_RETRIEVAL_APPJOB_H

#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_data_access.h>

#include <cassert>
#include <vector>

#include "common/forwarding_logger.h"

//----------------------------------------------------------------------
/// Retrieve heightfield data from cluster nodes.
class Heightfield_data_retrieval_appjob : 
    public mi::neuraylib::Fragmented_job<0xec9d7a1c,0xf03d,0x44de,0x8d,0x7e,0x80,0xd1,0xfb,0x02,0x93,0xf5>
{
public:
    /// constructor
    Heightfield_data_retrieval_appjob(
        const mi::neuraylib::Tag&                                                   heightfield_tag,
        const std::vector<mi::math::Bbox<mi::Sint32, 2> >&                          bounds,
        const mi::base::Handle<const nv::index::IDistributed_data_access_factory>&  data_access_factory);

    /// default constructor
    Heightfield_data_retrieval_appjob()
        :
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_bounds(),
        m_data_access_factory(),
        m_accessed_heightfield()
    {
        // do nothing, just for serialization
    }

    /// destructor
    virtual ~Heightfield_data_retrieval_appjob()
    {
        // do nothing as handles are used
    }

    /// Access the queried heightfield data ...
    const mi::base::Handle<nv::index::IRegular_heightfield_data_access>& get_heightfield_data(
        mi::Uint32 index) const
    {
        return m_accessed_heightfield[index];
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
    mi::neuraylib::Tag                                                  m_heightfield_tag;
    std::vector<mi::math::Bbox<mi::Sint32, 2> >                         m_bounds;
    mi::base::Handle<const nv::index::IDistributed_data_access_factory> m_data_access_factory;

    // Resulting heightfield data, do not serialize!
    std::vector<mi::base::Handle<nv::index::IRegular_heightfield_data_access> > m_accessed_heightfield;
};

//----------------------------------------------------------------------
#endif // GEOSPATIAL_STREAM_VIEWER_HEIGHTFIELD_DATA_RETRIEVAL_APPJOB_H
