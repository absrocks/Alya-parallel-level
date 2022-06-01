/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief IndeX and dice configuration utiliies

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_CONFIG_UTILILITY_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_CONFIG_UTILILITY_H

#include <mi/dice.h>
#include <nv/index/iindex.h>
#include <nv/index/iline_set.h>

#include <string>
#include <vector>
#include <cassert>

#include "common/forwarding_logger.h"
#include "common/receiving_logger.h"
#include "common/string_dict.h"

//----------------------------------------------------------------------
/// configure ILog module
///
/// \param[in] iindex_if IIndex API object
/// \param[in] app_project          application project setting dict
inline void configure_log_module(mi::base::Handle<nv::index::IIndex>& iindex_if,
                                 const nv::index_common::String_dict& app_project)
{
    assert(iindex_if.is_valid_interface());

    // Initialize the forwarding logger
    nv::index_common::Forwarding_logger_factory::instance()->initialize(iindex_if);

    mi::base::Handle<mi::neuraylib::ILogging_configuration> logging_configuration(
        iindex_if->get_api_component<mi::neuraylib::ILogging_configuration>());
    assert(logging_configuration.is_valid_interface());

    // Logging format
    mi::Uint32 log_prefix =
        mi::neuraylib::LOG_PREFIX_MODULE   +
        mi::neuraylib::LOG_PREFIX_CATEGORY +
        mi::neuraylib::LOG_PREFIX_SEVERITY;

    // Print host.thread only when network mode is active
    if (app_project.get("dice::network::mode", "OFF") != "OFF")
    {
        log_prefix += mi::neuraylib::LOG_PREFIX_HOST_THREAD;
    }

    // Timestamp logging
    if (nv::index_common::get_bool(app_project.get("dice::log_timestamp", "no")))
    {
        log_prefix += mi::neuraylib::LOG_PREFIX_TIME;
    }

    // Host name logging
    if (nv::index_common::get_bool(app_project.get("dice::log_hostname", "no")))
    {
        log_prefix += mi::neuraylib::LOG_PREFIX_HOST_NAME;
    }

    logging_configuration->set_log_prefix(log_prefix);

    // Set the log level (severity)
    if(app_project.is_defined("dice::verbose"))
    {
        const mi::Uint32 lv_uint32 = nv::index_common::get_uint32(app_project.get("dice::verbose"));

        INFO_LOG << "Setting log level to " << lv_uint32 << " ("
                 << nv::index_common::Forwarding_logger::level_to_string(mi::base::Message_severity(lv_uint32))
                 << ")";

        logging_configuration->set_log_level(mi::base::Message_severity(lv_uint32));

        // Also set log level for the category used by IndeX
        logging_configuration->set_log_level_by_category("MAIN", mi::base::Message_severity(lv_uint32));

        nv::index_common::Forwarding_logger_factory::instance()->set_fallback_log_severity(lv_uint32);
    }

    // Install the receiving logger.
    // This should be done after setting the logging configuration, so that delayed log messages are
    // also formated correctly, since they are emmitted when set_receiving_logger() is called.
    mi::base::Handle<mi::base::ILogger> receiving_logger(new nv::index_common::Receiving_logger());
    logging_configuration->set_receiving_logger(receiving_logger.get());


    // Handle 'dice::log_locally'
    mi::base::Handle<mi::neuraylib::IDebug_configuration> idebug_configuration(
        iindex_if->get_api_component<mi::neuraylib::IDebug_configuration>());
    assert(idebug_configuration.is_valid_interface());
    {
        const bool is_dice_log_locally = nv::index_common::get_bool(app_project.get("dice::log_locally", "0"));
        std::stringstream sstr;
        sstr << "log_locally=" << (is_dice_log_locally ? "1" : "0");
        idebug_configuration->set_option(sstr.str().c_str());
        INFO_LOG << "Local logging " << (is_dice_log_locally ? "enabled" : "disabled");
    }
}

//----------------------------------------------------------------------
/// string pair vector comparison function with the first of the pair.
/// For the get_dice_idebug_configuration_vector().
class String_pair_vector_comp_with_first
{
public:
    /// return the component ordering
    int operator()(std::pair< std::string, std::string > const arg0,
                   std::pair< std::string, std::string > const arg1)
    {
        return (arg0.first < arg1.first);
    }
};

//----------------------------------------------------------------------
/// get debug configuration vector
///
/// Filter with the key prefix (e.g., "dice::debug_configuration::",
/// "index::debug_configuration::") and return a sorted
/// vector according to the key string.
///
/// \param[in] opt_dict option string dictionary
/// \param[in] delete_prefix remove the prefix from the result
/// \return sorted vector of {key,value} vector. The values are strings.
inline std::vector< std::pair< std::string, std::string > >
get_configuration_vector(
    const std::string& config_key_prefix,
    const nv::index_common::String_dict& opt_dict,
    bool  delete_prefix = false)
{
    nv::index_common::String_dict extra_opt_dict;
    nv::index_common::string_dict_key_prefix_filter(opt_dict, config_key_prefix, extra_opt_dict, delete_prefix);
    std::vector< std::pair< std::string, std::string > > retvec;
    for (nv::index_common::String_dict::const_iterator ii = extra_opt_dict.begin(); 
        ii != extra_opt_dict.end();
        ++ii)
    {
        retvec.push_back(std::make_pair(ii->first, ii->second));
    }
    std::sort(retvec.begin(), retvec.end(), String_pair_vector_comp_with_first());

    return retvec;
}

//----------------------------------------------------------------------
/// set idebug configuration option helper function
///
/// \param[in] idebug_configuration a valid idebug configuration
/// \param[in] opt    dice debug option to be set
/// \param[in] mes    message output (opt prefix)
/// \param[in,out] set_counter count up this when set happens.
/// \return false when error. true when no error (includes ignored)
inline bool set_idebug_configuration_option_sub(
    mi::base::Handle<mi::neuraylib::IDebug_configuration> & idebug_configuration,
    std::string const & opt,
    std::string const & mes,
    mi::Sint32 & set_counter)
{
    assert(idebug_configuration.is_valid_interface());
    INFO_LOG << mes;
    if(idebug_configuration->set_option(opt.c_str()) == 0){
        ++set_counter;
        return true;
    }

    ERROR_LOG << mes << " failed";
    return false;
}

//----------------------------------------------------------------------
/// set dice network configuration
///
/// This configuration settings are actually not officially supported,
/// therefore no get method available.
///
/// Accepted parameter keys (see details in the project file documentation)
///
/// - dice::network::max_bandwidth
/// - dice::network::unicast_nak_interval
/// - dice::network::bandwidth_increment
/// - dice::network::bandwidth_decrement
/// - dice::network::retransmission_interval
/// - dice::network::additional_unicast_sockets
/// - dice::network::send_elements_only_to_owners
/// - dice::network::disk_cache_path
/// - dice::network::retention
///
/// \param[in] opt_dict dice network configuration options
/// \return number of set options, -1 when error.
inline mi::Sint32 set_config_dice_network_setting(
    mi::base::Handle<mi::neuraylib::IDebug_configuration> & idebug_configuration,
    nv::index_common::String_dict const & opt_dict)
{
    assert(idebug_configuration.is_valid_interface());

    mi::Sint32 set_count = 0;
    bool is_ok = true;

    // dice network configuration.
    // These parameters have not yet officially supported by DiCE.

    // The maximum bandwidth leveraged by DiCE for packet transfers
    // Default is 1000000000 bit, i.e., 1 Gigabit.
    const std::string max_bandwidth = opt_dict.get("dice::network::max_bandwidth", "");
    if(!max_bandwidth.empty())
    {
        std::string const opt = "max_bandwidth=" + max_bandwidth;
        std::string const mes = "Setting dice::network::" + opt;
        is_ok = is_ok && set_idebug_configuration_option_sub(idebug_configuration, opt, mes, set_count);
    }

    // The unicast_nak_interval controls how fast missing packets are
    // reported in case of unicast transmissions.
    const std::string unicast_nak_interval = opt_dict.get("dice::network::unicast_nak_interval", "");
    if(!unicast_nak_interval.empty())
    {
        std::string const opt = "unicast_nak_interval=" + unicast_nak_interval;
        std::string const mes = "Setting dice::network::" + opt;
        is_ok = is_ok && set_idebug_configuration_option_sub(idebug_configuration, opt, mes, set_count);
    }

    // The bandwidth increment controls how fast the bandwidth increases
    // if a packet was confirmed as going through without loss (in bit/s).
    const std::string bandwidth_increment = opt_dict.get("dice::network::bandwidth_increment", "");
    if(!bandwidth_increment.empty())
    {
        std::string const opt = "bandwidth_increment=" + bandwidth_increment;
        std::string const mes = "Setting dice::network::" + opt;
        is_ok = is_ok && set_idebug_configuration_option_sub(idebug_configuration, opt, mes, set_count);
    }

    // The bandwidth decrement controls how much the bandwidth is cut down
    // in case of a lost packet. The factor is multiplied with the bandwidth.
    const std::string bandwidth_decrement = opt_dict.get("dice::network::bandwidth_decrement", "");
    if(!bandwidth_decrement.empty())
    {
        std::string const opt = "bandwidth_decrement=" + bandwidth_decrement;
        std::string const mes = "Setting dice::network::" + opt;
        is_ok = is_ok && set_idebug_configuration_option_sub(idebug_configuration, opt, mes, set_count);
    }

    const std::string retransmission_interval = opt_dict.get("dice::network::retransmission_interval", "");
    if(!retransmission_interval.empty())
    {
        std::string const opt = "retransmission_interval=" + retransmission_interval;
        std::string const mes = "Setting dice::network::" + opt;
        is_ok = is_ok && set_idebug_configuration_option_sub(idebug_configuration, opt, mes, set_count);
    }

    const mi::Uint32 additional_unicast_sockets =
        nv::index_common::get_uint32(opt_dict.get("dice::network::additional_unicast_sockets", "0"));
    if(additional_unicast_sockets > 0)
    {
        std::stringstream sstr;
        sstr << "additional_unicast_sockets=" << additional_unicast_sockets;
        std::string const opt = sstr.str();
        std::string const mes = "Setting dice::network::" + opt;
        is_ok = is_ok && set_idebug_configuration_option_sub(idebug_configuration, opt, mes, set_count);
    }

    const bool send_elements_only_to_owners = 
        nv::index_common::get_bool(opt_dict.get("dice::network::send_elements_only_to_owners", "0"));
    if(send_elements_only_to_owners)
    {
        std::stringstream sstr;
        sstr << "send_elements_only_to_owners=1";
        std::string const opt = sstr.str();
        std::string const mes = "Setting dice::network::" + opt;
        is_ok = is_ok && set_idebug_configuration_option_sub(idebug_configuration, opt, mes, set_count);
    }

    const std::string disk_cache_path = opt_dict.get("dice::network::disk_cache_path", "\".\"");
    if(disk_cache_path != "\".\"")
    {
        std::stringstream sstr;
        sstr << "disk_cache_path=" << disk_cache_path;
        std::string const opt = sstr.str();
        std::string const mes = "Setting dice::network::" + opt;
        is_ok = is_ok && set_idebug_configuration_option_sub(idebug_configuration, opt, mes, set_count);
    }

    mi::Sint32 const dice_retention = nv::index_common::get_sint32(opt_dict.get("dice::network::retention", "-1"));
    if ((dice_retention >= 0 && dice_retention <= 1) || (dice_retention > 6000)) {
        WARN_LOG << "DiCE retention parameter has an extreme value [" << dice_retention
                 << "]. Accepted as it is, but this affects performance.";
    }
    if(dice_retention >= 0)
    {
        std::stringstream sstr;
        sstr << "retention=" << dice_retention;
        std::string const opt = sstr.str();
        std::string const mes = "Setting dice::network::" + opt;
        is_ok = is_ok && set_idebug_configuration_option_sub(idebug_configuration, opt, mes, set_count);
    }

    if(is_ok){
        return set_count;
    }
    // not ok now
    return -1;
}

//----------------------------------------------------------------------
#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_CONFIG_UTILILITY_H
