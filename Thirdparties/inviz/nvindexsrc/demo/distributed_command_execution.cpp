/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "distributed_command_execution.h"

#include <iostream>

#include <nv/index/isession.h>

#include "common/common_utility.h"
#include "common/forwarding_logger.h"

Distributed_command_execution::Distributed_command_execution(
    const std::string& command,
    mi::neuraylib::Tag session_tag)
  : m_command(command),
    m_session_tag(session_tag)
{
}

Distributed_command_execution::Distributed_command_execution()
{
}

mi::Size Distributed_command_execution::get_nb_of_fragments() const
{
    // Together with scheduling mode ONCE_PER_HOST, 0 means run exactly one fragment per remote host
    return 0;
}

mi::neuraylib::IFragmented_job::Scheduling_mode Distributed_command_execution::get_scheduling_mode() const
{
    return mi::neuraylib::IFragmented_job::ONCE_PER_HOST;
}

void Distributed_command_execution::execute_fragment(
    mi::neuraylib::IDice_transaction*               dice_transaction,
    mi::Size                                        index,
    mi::Size                                        count,
    const mi::neuraylib::IJob_execution_context*    context)
{
    run_command(dice_transaction);
}

void Distributed_command_execution::execute_fragment_remote(
    mi::neuraylib::ISerializer*                   serializer,
    mi::neuraylib::IDice_transaction*             dice_transaction,
    mi::Size                                      index,
    mi::Size                                      count,
    const mi::neuraylib::IJob_execution_context*  context)
{
    run_command(dice_transaction);
}

void Distributed_command_execution::serialize(mi::neuraylib::ISerializer *serializer) const
{
    mi::Size nb_elements = m_command.size();
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(&m_command[0]), nb_elements);

    serializer->write(&m_session_tag.id, 1);
}

void Distributed_command_execution::deserialize(mi::neuraylib::IDeserializer* deserializer)
{
    mi::Size nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_command.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_command[0]), nb_elements);

    deserializer->read(&m_session_tag.id, 1);
}

void Distributed_command_execution::run_command(mi::neuraylib::IDice_transaction* dice_transaction)
{
    INFO_LOG << "running command '" << m_command << "' on host '" << nv::index_common::get_host_name() << "':";

    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<const nv::index::ISession>(m_session_tag));

    {
        INFO_LOG << "Unknown command '" << m_command << "'";
    }
}
