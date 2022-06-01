/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "command_handler_base.h"

#include <cassert>

//////////////////// Call_event_arguments ////////////////////

Call_event_arguments::Call_event_arguments(
    const mi::base::Handle<const mi::IMap>& args)
  : m_input(args)
{
}

std::string Call_event_arguments::get(const std::string& key)
{
    std::map<std::string, std::string>::const_iterator it = m_input_cache.find(key);
    if (it != m_input_cache.end())
    {
        return it->second;
    }

    mi::base::Handle<const mi::IString> value(m_input->get_value<mi::IString>(key.c_str()));
    if (value.is_valid_interface())
    {
        std::string s = value->get_c_str();
        m_input_cache[key] = s;
        return s;
    }
    else
    {
        WARN_LOG << "Call_event_handler: value for key '" << key << "' does not exist "
                 << "or is is not of string type";
        return "";
    }
}

bool Call_event_arguments::has(const std::string& key)
{
    std::map<std::string, std::string>::const_iterator it = m_input_cache.find(key);
    if (it != m_input_cache.end())
    {
        return true;
    }

    mi::base::Handle<const mi::IString> value(m_input->get_value<mi::IString>(key.c_str()));
    return  value.is_valid_interface();
}

void Call_event_arguments::set(
    const std::string& key,
    const std::string& value)
{
    m_output_cache[key] = value;
}

mi::IMap* Call_event_arguments::generate_response(
    const mi::base::Handle<nv::index::IIndex>& iindex_if) const
{
    assert(iindex_if.is_valid_interface());

    mi::base::Handle<mi::neuraylib::IFactory> factory(
        iindex_if->get_api_component<mi::neuraylib::IFactory>());
    assert(factory.is_valid_interface());

    mi::IMap* map = factory->create<mi::IMap>("Map<Interface>");
    assert(map != 0);

    std::map<std::string, std::string>::const_iterator it;
    for (it = m_output_cache.begin(); it != m_output_cache.end(); ++it)
    {
        const char* key = it->first.c_str();
        const char* value = it->second.c_str();
        mi::base::Handle<mi::IString> str(factory->create<mi::IString>());
        assert(str.is_valid_interface());
        str->set_c_str(value);
        if (map->has_key(key))
            map->set_value(key, str.get());
        else
            map->insert(key, str.get());
    }

    return map;
}

void Call_event_arguments::clear_response()
{
    m_output_cache.clear();
}

void Call_event_arguments::pass_through(const std::string& key)
{
    if (has(key))
        set(key, get(key));
}

//////////////////// ICall_event_command_handler ////////////////////

ICall_event_command_handler::ICall_event_command_handler(
    Nvindex_rendering_context& irc,
    Nvindex_AppData*           appdata)
  : m_irc(irc),
    m_appdata(appdata)
{
}

ICall_event_command_handler::~ICall_event_command_handler()
{
}
