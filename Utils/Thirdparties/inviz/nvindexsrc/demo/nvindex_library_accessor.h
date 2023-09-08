/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief nvindex library access utilities

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_NVINDEX_LIBRARY_ACCESSOR_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_NVINDEX_LIBRARY_ACCESSOR_H

#include <sstream>
#include <string>
#include <stdio.h>
#include <iostream>

#include <mi/dice.h>
#include <mi/base/handle.h>
#include <mi/base/ilogger.h>

#include <nv/index/iindex.h>

namespace nv {
namespace index_common {
class String_dict;
}} // namespace
class Nvindex_rendering_context;

/// Creating and accessing an instance of IndeX library interface
class Nvindex_library_accessor
{
public:
    /// constructor
    Nvindex_library_accessor();

    /// destructor
    virtual ~Nvindex_library_accessor();

    /// is initialized dice library
    /// \return true when it has been initialized.
    bool is_initialized() const
    {
        return m_nvindex_interface.is_valid_interface();
    }

    /// initialize library accessor
    ///
    /// \param[in] app_project application project options
    /// \param[in] irc         IndeX rendering context
    bool initialize(
        const nv::index_common::String_dict& app_project,
        Nvindex_rendering_context&           irc);

    /// get nvindex library interface
    /// \return nvindex library interface
    mi::base::Handle<nv::index::IIndex>& get_interface()
    {
        return m_nvindex_interface;
    }

private:
    /// access IndeX interface. loading nvindex library.
    /// \param[in] nvindex_library_variant nvindex library name
    /// \return IIndex interface
    nv::index::IIndex* access_index_interface(
        const std::string& nvindex_library_variant);

    /// clean up when shutdown the interface. but not shutdown
    /// nvindex_interface itself. This is for the failure case.
    void clear();

private:
    /// nvindex interface
    mi::base::Handle<nv::index::IIndex> m_nvindex_interface;

private:
    /// copy constructor. prohibit until proved useful.
    Nvindex_library_accessor(Nvindex_library_accessor const &);
    /// operator=. prohibit until proved useful.
    Nvindex_library_accessor const & operator=(Nvindex_library_accessor const &);
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_NVINDEX_LIBRARY_ACCESSOR_H

