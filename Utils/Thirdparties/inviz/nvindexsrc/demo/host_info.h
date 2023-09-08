/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief IndeX connected host information

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HOST_INFO_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HOST_INFO_H

#include <mi/dice.h>

#include <algorithm>
#include <map>
#include <string>
#include <vector>

/// NVIDIA IndeX connected host information
class Host_info
{
public:
    // comparison object
    struct id_comp
    {
        bool operator()(const std::pair<mi::Uint32, std::string> & rhs,
                        const std::pair<mi::Uint32, std::string> & lhs)
        {
            return rhs.first < lhs.first;
        }
    };

public:
    /// default constructor
    Host_info();
    /// destructor
    ~Host_info();

    /// set local host id
    /// \param[in] local_host_id current local host id
    void set_local_host_id(mi::Uint32 local_host_id);
    /// get local host id
    /// \return current local host id
    mi::Uint32 get_local_host_id() const;
    /// get local host name
    /// \return local host name
    std::string get_local_host_name() const;
    /// set host id to host name map
    /// \param[in] hi_hn_map a host id to host name map
    void set_host_id_host_name_map(const std::map<mi::Uint32, std::string> & hi_hn_map);
    /// get host id to host name map
    /// \return a host id to host name map
    std::map<mi::Uint32, std::string> get_host_id_host_name_map() const;
    /// get current host map size
    /// \return size of the map 
    mi::Size size() const;
    /// get sorted host_name vec
    /// \return host_name vector sorted by host id
    std::vector<std::string> get_sorted_host_name_vec() const;

private:
    /// updated sorted (by host id) host vector
    void update_sorted_host_vector();

private:
    // local host ID
    mi::Uint32 m_local_host_id;
    // map: host ID -> host_name 
    std::map<mi::Uint32, std::string> m_host_id_host_name_map;
    // id sorted host list
    std::vector<std::string> m_sorted_host_vec;

private:
    /// copy constructor. prohibit until proved useful.
    Host_info(Host_info const &);
    /// operator=. prohibit until proved useful.
    Host_info const & operator=(Host_info const &);
};

//======================================================================

/// forward declaration
class Nvindex_rendering_context;

//----------------------------------------------------------------------
/// update the current host map information
///
/// \param[in]     irc       IndeX rendering context
/// \param[in,out] hostinfo  host map information (This will be updated.)
void update_host_map(Nvindex_rendering_context& irc,
                     Host_info & hostinfo);

//----------------------------------------------------------------------
#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_HOST_INFO_H
