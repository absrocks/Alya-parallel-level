/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
///

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_DISTRIBUTED_COMMAND_EXECUTION_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_DISTRIBUTED_COMMAND_EXECUTION_H

#include <string>

#include <mi/neuraylib/dice.h>

// Passes a string command to all hosts in the cluster and executes it there.
class Distributed_command_execution :
    public mi::neuraylib::Fragmented_job<0x0b7e59a2,0x6308,0x4ee4,0xbb,0x0c,0xce,0x6b,0x7c,0xdb,0x71,0xb8>
{
public:
    Distributed_command_execution(
        const std::string& command,
        mi::neuraylib::Tag session_tag);

    Distributed_command_execution();

    mi::Size get_nb_of_fragments() const;

    virtual mi::neuraylib::IFragmented_job::Scheduling_mode get_scheduling_mode() const;

    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

   virtual void execute_fragment_remote(
        mi::neuraylib::ISerializer*                   serializer,
        mi::neuraylib::IDice_transaction*             dice_transaction,
        mi::Size                                      index,
        mi::Size                                      count,
        const mi::neuraylib::IJob_execution_context*  context);

    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    void run_command(mi::neuraylib::IDice_transaction* dice_transaction);

    std::string                    m_command;
    mi::neuraylib::Tag             m_session_tag;
};

#endif // NVIDIA_INDEX_DISTRIBUTED_COMMAND_EXECUTION_H
