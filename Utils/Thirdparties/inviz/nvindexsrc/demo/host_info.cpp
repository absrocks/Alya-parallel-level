/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief IndeX connected host information

#include "host_info.h"

#include "nvindex_rendering_context.h"
#include "nvindex_appdata.h"

#include <cassert>
#include <iterator>
#include <sstream>
#include <vector>

//----------------------------------------------------------------------
Host_info::Host_info() :
    m_local_host_id(0)
{
    // empty
}
//----------------------------------------------------------------------
Host_info::~Host_info()
{
    // empty
}

//----------------------------------------------------------------------
void Host_info::set_local_host_id(mi::Uint32 local_host_id)
{
    m_local_host_id = local_host_id;
}

//----------------------------------------------------------------------
mi::Uint32 Host_info::get_local_host_id() const
{
    return m_local_host_id;
}

//----------------------------------------------------------------------
std::string Host_info::get_local_host_name() const
{
    std::map<mi::Uint32, std::string>::const_iterator hi = m_host_id_host_name_map.find(this->get_local_host_id());
    if(hi != m_host_id_host_name_map.end()){
        return hi->second;
    }

    WARN_LOG << "Cannot retrieve local host name.";
    return "localhost";
}

//----------------------------------------------------------------------
void Host_info::set_host_id_host_name_map(const std::map<mi::Uint32, std::string> & hi_hn_map)
{
    m_host_id_host_name_map = hi_hn_map;
    update_sorted_host_vector();
}

//----------------------------------------------------------------------
std::map<mi::Uint32, std::string> Host_info::get_host_id_host_name_map() const
{
    return m_host_id_host_name_map;
}

//----------------------------------------------------------------------
mi::Size Host_info::size() const
{
    return m_host_id_host_name_map.size();
}

//----------------------------------------------------------------------
std::vector<std::string> Host_info::get_sorted_host_name_vec() const
{
    return m_sorted_host_vec;
}

//----------------------------------------------------------------------
void Host_info::update_sorted_host_vector()
{
    m_sorted_host_vec.clear();

    const mi::Size len = m_host_id_host_name_map.size();
    if(len == 0){
        return;             // nothing, no need to sort
    }
            
    std::vector<std::pair<mi::Uint32, std::string> > sort_buf;
    for(std::map<mi::Uint32, std::string>::const_iterator hi = m_host_id_host_name_map.begin();
        hi != m_host_id_host_name_map.end();
        ++hi)
    {
        sort_buf.push_back(*hi);
    }
    assert(sort_buf.size() == len);
    std::sort(sort_buf.begin(), sort_buf.end(), id_comp());
    for(std::vector<std::pair<mi::Uint32, std::string> >::const_iterator hi = sort_buf.begin();
        hi != sort_buf.end();
        ++hi)
    {
        m_sorted_host_vec.push_back(hi->second);
    }
}

//======================================================================

void update_host_map(Nvindex_rendering_context& irc,
                     Host_info & hostinfo)
{
    // At this point, we don't know the number of active hosts.
    // Therefore, we use DiCE joined number of hosts until updated.
    mi::Uint32 nb_hosts = irc.m_icluster_configuration->get_number_of_hosts();

    if(nb_hosts == hostinfo.size()){
        return;
    }

    // Retrieve the local host id.
    const mi::Uint32 local_host_id = irc.m_icluster_configuration->get_local_host_id();

    std::map<mi::Uint32, std::string> host_id_host_name_map;

    // Retrieve host_names of all connected hosts.
    for (mi::Uint32 i = 0; i < nb_hosts; ++i)
    {
        const mi::Uint32 host_id = irc.m_icluster_configuration->get_host_index(i);

        std::ostringstream os;
        os << "Host " << host_id;
        std::string host_name = os.str();
        
        if (host_id == local_host_id)
        {
            // Add our local host_name as stored in the DiCE config
            mi::base::Handle<mi::neuraylib::IGeneral_configuration> igeneral_configuration(
                irc.m_iindex_if->get_api_component<mi::neuraylib::IGeneral_configuration>());
            assert(igeneral_configuration.is_valid_interface());
            
            mi::base::Handle<const mi::neuraylib::IHost_properties> host_prop(
                igeneral_configuration->get_host_properties());
            assert(host_prop.is_valid_interface());
            mi::base::Handle<const mi::IString > hoststr(host_prop->get_property("host_name"));
            if (hoststr.is_valid_interface())
            {
                host_name = hoststr->get_c_str();
            }
            else
            {
                host_name = "unknown_localhost"; // no network mode
            }
        }
        else
        {
            // Get remote host_name
            const char* s = irc.m_icluster_configuration->get_host_name(host_id);
            if (s == 0){
                INFO_LOG << "Cannot retrieve the remote host_name for host id " << host_name << ", "
                         << "setting the name to '" << host_name << "'";
            }
            else{
                host_name = s;
            }
        }
        host_id_host_name_map[host_id] = host_name;
    }

    // Set updated host id -> host name map
    hostinfo.set_host_id_host_name_map(host_id_host_name_map);
    // Set local host id
    hostinfo.set_local_host_id(local_host_id);

    // updated the app_proj host list
    std::vector<std::string> hostvec = hostinfo.get_sorted_host_name_vec();
    std::stringstream sstr;
    std::copy(hostvec.begin(), hostvec.end(), std::ostream_iterator<std::string>(sstr, " "));

    Nvindex_AppData::instance()->peek_app_proj()->insert("app::performance::host_list",
                                                         sstr.str());
    Nvindex_AppData::instance()->peek_app_proj()->insert("app::performance::local_host_name", 
                                                         hostinfo.get_local_host_name());

}

//----------------------------------------------------------------------
