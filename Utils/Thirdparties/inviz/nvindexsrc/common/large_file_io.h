/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief File abstraction for transparent handling of large (>4GiB) on
///        Linux and Windows platforms.

#ifndef NVIDIA_INDEX_BIN_COMMON_LARGE_FILE_IO_H
#define NVIDIA_INDEX_BIN_COMMON_LARGE_FILE_IO_H

#include <ios>
#include <cstddef>
#include <string>

#include <mi/base/types.h>

namespace nv {
namespace index_common {
namespace io {

/// File abstraction for a large (>4GiB) file.
class File
{
public:
    File();
    explicit File(mi::Uint32 gzip_compression_level);
    explicit File(const std::string&       file_path,
                  std::ios_base::openmode  open_mode = std::ios_base::in | std::ios_base::out);
    virtual ~File();

    bool                open(const std::string&       file_path,
                             std::ios_base::openmode  open_mode = std::ios_base::in | std::ios_base::out);
    bool                is_open() const;
    void                close();

    mi::Uint64          read(void*          output_buffer,
                             mi::Uint64     start_position,
                             mi::Uint64     num_bytes_to_read);
    mi::Uint64          write(const void*   input_buffer,
                              mi::Uint64    start_position,
                              mi::Uint64    num_bytes_to_write);

    bool                flush_buffers() const;

    mi::Uint64          size() const;
    const std::string&  file_path() const;

    operator            bool() const;

public:
    class File_impl;

protected:
    File_impl*          m_file_impl;

private:
    // declared, never defined
    File(const File&) /*= delete*/;
    File& operator=(const File&) /*= delete*/;
}; // class File

} // namespace io
}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_LARGE_FILE_IO_H
