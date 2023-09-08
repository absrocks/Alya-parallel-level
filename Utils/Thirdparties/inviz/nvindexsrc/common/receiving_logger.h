/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Implementation example of a receiving message logger

#ifndef NVIDIA_INDEX_BIN_COMMON_RECEIVING_LOGGER_H
#define NVIDIA_INDEX_BIN_COMMON_RECEIVING_LOGGER_H

#include <sstream>
#include <string>

#include <mi/dice.h>
#include <mi/base/ilogger.h>

namespace nv {
namespace index_common {
//----------------------------------------------------------------------
/// A simple thread-safe receiving logger
class Receiving_logger : public mi::base::Interface_implement<mi::base::ILogger>
{
public:
    /// Constructor
    Receiving_logger();

    /// Destructor. Message out.
    virtual ~Receiving_logger();

    /// Get message stream
    /// \param[in] level message severity level
    std::ostringstream& get_message(mi::Uint32 level);

    virtual void message(
        mi::base::Message_severity level,
        const char*                category,
        const char*                message);

private:
    /// Output stream
    std::ostringstream      m_os;
    /// Message severity level
    mi::Uint32              m_level;
    /// Use color output on terminal
    bool                    m_color;

    Receiving_logger(const Receiving_logger&);               // forbid this
    Receiving_logger& operator=(const Receiving_logger&);    // forbid this
};
//----------------------------------------------------------------------

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_RECEIVING_LOGGER_H
