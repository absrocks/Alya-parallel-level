/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "large_file_io.h"

#include <cassert>
#include <fstream>

#include "forwarding_logger.h"

#ifdef _WIN32
#   ifndef WIN32_LEAN_AND_MEAN
#       define WIN32_LEAN_AND_MEAN 1
#   endif
#   include <windows.h>
#else // _WIN32
#   include <sys/types.h>
#   include <sys/stat.h>
#   include <errno.h>
#   include <fcntl.h>
#   include <unistd.h>
#endif // _WIN32

#include <zlib.h>

namespace nv {
namespace index_common {
namespace io {
/// namespace util for io
namespace util {
/// check file exists
bool file_exists(const std::string& file_path)
{
    using namespace std;
    ifstream f(file_path.c_str(), ios_base::in);

    return f.is_open();
}

} // namespace util

#ifdef _WIN32

class File::File_impl
{
public:
    File_impl()
    {
        m_file_handle = INVALID_HANDLE_VALUE;
        m_position    = 0;
        m_file_size   = 0;
        m_open_mode   = static_cast<std::ios_base::openmode>(0);
    }

    virtual ~File_impl()
    {
        close();
    }

    virtual bool open(const std::string&       file_path,
              std::ios_base::openmode  open_mode)
    {
        // translate open mode to access mode
        DWORD desired_access           = 0;
        DWORD read_write_buffer_access = 0;

        if (open_mode & std::ios_base::in) {
            desired_access           |= GENERIC_READ;
            read_write_buffer_access  = PAGE_READWRITE;
        }
        if (open_mode & std::ios_base::out) {
            desired_access           |= GENERIC_WRITE;
            read_write_buffer_access  = PAGE_READWRITE;
        }

        // share mode
        DWORD share_mode = FILE_SHARE_READ;

        // translate open mode to creation modes
        DWORD creation_disposition = 0;

        if (util::file_exists(file_path)) {
            if (    ((open_mode & std::ios_base::out)
                  || (open_mode & std::ios_base::in))
                && !(open_mode & std::ios_base::trunc)) {
                creation_disposition = OPEN_ALWAYS;
            }
            else if (   (open_mode & std::ios_base::out)
                     && (open_mode & std::ios_base::trunc)) {
                creation_disposition = TRUNCATE_EXISTING;
            }
            else {
                DEBUG_LOG << "File::open(): "
                          << "illegal open mode 0x" << std::hex << open_mode << "on existing file "
                          <<  "'" << file_path << "'";
                return false;
            }
        }
        else {
            if (      (open_mode & std::ios_base::out)
                  //&& !(open_mode & std::ios_base::in)
                  && !(open_mode & std::ios_base::trunc)) {
                creation_disposition = CREATE_NEW;
            }
            else {
                //DEBUG_LOG << "File::open(): "
                //          << "illegal open mode 0x" << std::hex << open_mode << " on non existing file "
                //          <<  "'" << file_path << "'";
                return false;
            }
        }

        // file attributes
        DWORD flags_and_attributes = FILE_ATTRIBUTE_NORMAL;// | FILE_FLAG_SEQUENTIAL_SCAN;

        // use shared pointer to manage the close in case we miss it somehow
        m_file_handle = CreateFile(file_path.c_str(),
                                   desired_access,
                                   share_mode,
                                   0,
                                   creation_disposition,
                                   flags_and_attributes,
                                   0);

        if (m_file_handle == INVALID_HANDLE_VALUE) {
            DEBUG_LOG << "File::open(): "
                      << "error creating/opening file: "
                      <<  "'" << file_path << "'";
            return false;
        }

        if (   (open_mode & std::ios_base::ate)
            || (open_mode & std::ios_base::app)) {

            m_position = actual_file_size();
        }

        m_file_path  = file_path;
        m_open_mode  = open_mode;
        m_file_size  = actual_file_size();

        return true;
    }
    
    virtual bool is_open() const
    {
        if (m_file_handle == INVALID_HANDLE_VALUE) {
            return false;
        }
        else {
            return true;
        }
    }

    virtual void close()
    {
        if (is_open()) {
            CloseHandle(m_file_handle);
        }

        m_file_handle = INVALID_HANDLE_VALUE;
        m_position    = 0;
        m_file_size   = 0;
        m_open_mode   = static_cast<std::ios_base::openmode>(0);
        m_file_path.clear();
    }

    virtual mi::Uint64 read(
        void*          output_buffer,
        mi::Uint64     start_position,
        mi::Uint64     num_bytes_to_read)
    {
        if (!is_open()) {
            DEBUG_LOG << "File::read(): "
                      << "read access on invalid file handle.";
            return 0;
        }

        using namespace mi;

        Uint8*      output_byte_buffer  = reinterpret_cast<Uint8*>(output_buffer);
        Uint64      bytes_read          = 0;

        if (!set_file_pointer(start_position)) {
            DEBUG_LOG << "File::read(): "
                      << "unable to set file pointer to current position.";
            return 0;
        }
        m_position = start_position;

        // TODO loop for read requests larger than 4GiB!
        DWORD file_bytes_read = 0;

        if (ReadFile(m_file_handle, output_byte_buffer, static_cast<DWORD>(num_bytes_to_read), &file_bytes_read, 0) == 0) {
            DEBUG_LOG << "File::read(): "
                      << "error reading from file " << m_file_path;
            return 0;
        }

        if (file_bytes_read == 0) {
            // eof
            return 0;
        }

        if (file_bytes_read <= num_bytes_to_read) {
            m_position  += file_bytes_read;
            bytes_read   = file_bytes_read;
        }
        else {
            DEBUG_LOG << "File::read(): "
                      << "unknown error reading from file " << m_file_path;
            return 0;
        }

        assert(bytes_read > 0);

        return bytes_read;
    }

    virtual mi::Uint64 write(
        const void*   input_buffer,
        mi::Uint64    start_position,
        mi::Uint64    num_bytes_to_write)
    {
        using namespace mi;

        if (!is_open()) {
            DEBUG_LOG << "File::write(): "
                      << "write access on invalid file handle.";
            return 0;
        }

        if (m_open_mode & std::ios_base::app) {
            m_position = m_file_size;
        }

        const Uint8*    input_byte_buffer   = reinterpret_cast<const Uint8*>(input_buffer);
        Uint64          bytes_written       = 0;

        m_position = start_position;
        if (!set_file_pointer(m_position)) {
            return 0;
        }

        // TODO loop for write requests lager than 4GiB!
        DWORD   file_bytes_written  = 0;

        if (WriteFile(m_file_handle,
                      input_byte_buffer,
                      static_cast<DWORD>(num_bytes_to_write),
                      &file_bytes_written,
                      0) == 0) {
            DEBUG_LOG << "File::write(): "
                      << "error writing to file " << m_file_path;
            return 0;
        }

        if (file_bytes_written <= num_bytes_to_write) {
            m_position      += file_bytes_written;
            bytes_written    = file_bytes_written;
        }
        else {
            DEBUG_LOG << "File::write(): "
                      << "unknown error writing to file " << m_file_path;
        }

        return bytes_written;
    }

    virtual bool flush_buffers() const
    {
        if (!is_open()) {
            DEBUG_LOG << "File::flush_buffers(): "
                      << "write access on invalid file handle.";
            return false;
        }
        else {
            return FlushFileBuffers(m_file_handle) != FALSE ? true : false;
        }
    }

    virtual mi::Uint64 size() const
    {
        return m_file_size;
    }

    virtual const std::string& file_path() const
    {
        return m_file_path;
    }

    virtual mi::Uint64 actual_file_size() const
    {
        assert(m_file_handle != INVALID_HANDLE_VALUE);

        LARGE_INTEGER   cur_size_li;

        if (GetFileSizeEx(m_file_handle, &cur_size_li) == 0) {
            DEBUG_LOG << "File::actual_file_size(): "
                      << "error retrieving current file size: " << m_file_path;
            return 0;
        }

        return static_cast<mi::Uint64>(cur_size_li.QuadPart);
    }

    virtual bool set_file_pointer(mi::Uint64 new_pos)
    {
        assert(m_file_handle != INVALID_HANDLE_VALUE);

        LARGE_INTEGER   position_li;
        position_li.QuadPart = new_pos;

        if (   SetFilePointer(m_file_handle, position_li.LowPart, &position_li.HighPart, FILE_BEGIN) == INVALID_SET_FILE_POINTER
            && GetLastError() != NO_ERROR) {
            DEBUG_LOG << "File::set_file_pointer(): "
                      << "error setting file pointer to position "
                      << std::hex << new_pos
                      << " on file '" << m_file_path << "'";
            return false;
        }
        return true;
    }

protected:
    HANDLE                  m_file_handle;

    mi::Uint64              m_position;

    std::string             m_file_path;
    mi::Uint64              m_file_size;

    std::ios_base::openmode m_open_mode;

}; // class File::File_impl

#else // _WIN32

/// Implementation of internal class of File class
class File::File_impl
{
public:
    File_impl()
    {
        m_file_handle = -1;
        m_position    = 0;
        m_file_size   = 0;
        m_open_mode   = static_cast<std::ios_base::openmode>(0);
    }

    virtual ~File_impl()
    {
    }

    virtual bool open(const std::string&       file_path,
              std::ios_base::openmode  open_mode)
    {
        using namespace mi;

        mi::Sint32    open_flags  = O_LARGEFILE; // yes, we mainly go through this pain for large files
        mode_t create_mode = S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH;

        if (   (open_mode & std::ios_base::in)
            && (open_mode & std::ios_base::out)) {
            open_flags |= O_RDWR;
        }
        else if (open_mode & std::ios_base::in) {
            open_flags |= O_RDONLY;
        }
        else if (open_mode & std::ios_base::out) {
            open_flags |= O_WRONLY;
        }
        else {
            DEBUG_LOG << "File::open(): "
                      << "illegal open mode 0x" << std::hex << open_mode
                      <<  "'" << file_path << "'";
            return false;
        }

        if (util::file_exists(file_path)) {
            if (   ((open_mode & std::ios_base::out)
                ||  (open_mode & std::ios_base::in))
                && !(open_mode & std::ios_base::trunc)) {
                // everything ok
            }
            else if (   (open_mode & std::ios_base::out)
                     && (open_mode & std::ios_base::trunc)) {
                open_flags |= O_TRUNC;
            }
            else {
                DEBUG_LOG << "File::open(): "
                          << "illegal open mode 0x" << std::hex << open_mode << "on existing file "
                          <<  "'" << file_path << "'";
                return false;
            }
        }
        else {
            if (      (open_mode & std::ios_base::out)
                  //&& !(open_mode & std::ios_base::in)
                  && !(open_mode & std::ios_base::trunc)) {
                open_flags |= O_CREAT;
            }
            else {
                //DEBUG_LOG << "File::open(): "
                //          << "illegal open mode 0x" << std::hex << open_mode << " on non existing file "
                //          <<  "'" << file_path << "'";
                return false;
            }
        }

        m_file_handle = ::open64(file_path.c_str(), open_flags, create_mode);

        if (m_file_handle < 0) {
            std::string ret_error;
            switch (errno) {
                case EACCES:    ret_error.assign("requested access to the file is not allowed (insufficient permissions?)"); break;
                case EEXIST:    ret_error.assign("pathname already exists and O_CREAT and O_EXCL were used"); break;
                case EFAULT:    ret_error.assign("pathname points outside your accessible address space"); break;
                case EISDIR:    ret_error.assign("pathname refers to a directory and the access requested involved writing"); break;
                default:        ret_error.assign("unknown error"); break;
            }
            DEBUG_LOG << "File::open(): "
                      << "error opening file "
                      << "(" << ret_error << ")"
                      << " '" << file_path << "'";
            return false;
        }


        if (   (open_mode & std::ios_base::ate)
            || (open_mode & std::ios_base::app)) {

            m_position = actual_file_size();
        }

        m_file_path  = file_path;
        m_open_mode  = open_mode;
        m_file_size  = actual_file_size();

        return true;
    }

    virtual bool is_open() const
    {
        return m_file_handle > -1;
    }

    virtual void close()
    {
        if (is_open()) {
            if (0 != ::close(m_file_handle)) {
                std::string ret_error;
                switch (errno) {
                    case EBADF: ret_error.assign("invalid file descriptor"); break;
                    case EINTR: ret_error.assign("close interrupted by signal"); break;
                    case EIO:   ret_error.assign("I/O error"); break;
                    default:    ret_error.assign("unknown error"); break;
                }
                DEBUG_LOG << "File::close(): "
                           << "error closing file "
                           << "(" << ret_error << ")"
                           << " '" << m_file_path << "'";
            }
        }

        m_file_handle = -1;
        m_position    = 0;
        m_file_size   = 0;
        m_open_mode   = static_cast<std::ios_base::openmode>(0);
        m_file_path.clear();
    }

    virtual mi::Uint64 read(void*          output_buffer,
                            mi::Uint64     start_position,
                            mi::Uint64     num_bytes_to_read)
    {
        if (!is_open()) {
            DEBUG_LOG << "File::read(): "
                      << "read access on invalid file handle.";
            return 0;
        }

        using namespace mi;

        Uint8*      output_byte_buffer  = reinterpret_cast<Uint8*>(output_buffer);
        Uint64      bytes_read          = 0;

        m_position = start_position;

        if (num_bytes_to_read <= 0) {
            return 0;
        }

        while (bytes_read < num_bytes_to_read)
        {
            const Uint64 file_bytes_to_read = num_bytes_to_read - bytes_read;
            const ssize_t file_bytes_read = ::pread64(m_file_handle, output_byte_buffer, file_bytes_to_read, m_position);

            if (file_bytes_read < 0) {
                DEBUG_LOG << "File::read(): "
                          << "error reading from file " << m_file_path;
                return 0;
            }
            else if (file_bytes_read == 0) {
                break; // end of file
            }
            else {
                m_position  += file_bytes_read;
                bytes_read  += file_bytes_read;
            }
        }

        return bytes_read;
    }

    virtual mi::Uint64 write(const void*   input_buffer,
                             mi::Uint64    start_position,
                             mi::Uint64    num_bytes_to_write)
    {
        using namespace mi;

        if (!is_open()) {
            DEBUG_LOG << "File::write(): "
                      << "write access on invalid file handle.";
            return 0;
        }

        if (m_open_mode & std::ios_base::app) {
            m_position = m_file_size;
        }

        const Uint8*    input_byte_buffer   = reinterpret_cast<const Uint8*>(input_buffer);
        Uint64          bytes_written       = 0;

        m_position = start_position;

        while (bytes_written < num_bytes_to_write)
        {
            const Uint64 file_bytes_to_write = num_bytes_to_write - bytes_written;
            const ssize_t file_bytes_written = ::pwrite64(m_file_handle, input_byte_buffer, file_bytes_to_write, m_position);

            if (file_bytes_written < 0) {
                DEBUG_LOG << "File::write(): "
                          << "error writing to file " << m_file_path;
                return 0;
            }
            else if (file_bytes_written == 0) {
                break; // nothing was written
            }
            else {
                m_position    += file_bytes_written;
                bytes_written += file_bytes_written;
            }
        }

        return bytes_written;
    }

    virtual bool flush_buffers() const
    {
        if (!is_open()) {
            DEBUG_LOG << "File::flush_buffers(): "
                      << "write access on invalid file handle.";
            return false;
        }
        else {
            return true;
        }
    }

    virtual mi::Uint64 size() const
    {
        return m_file_size;
    }

    virtual const std::string& file_path() const
    {
        return m_file_path;
    }

    virtual mi::Uint64 actual_file_size() const
    {
        assert(is_open());

        off64_t new_pos = ::lseek64(m_file_handle, 0, SEEK_END);

        if (new_pos < 0) {
            std::string ret_error;
            switch (errno) {
                case EBADF:     ret_error.assign("invalid file descriptor"); break;
                case EINVAL:    ret_error.assign("invalid direction specified on lseek"); break;
                case EOVERFLOW: ret_error.assign("overflow on returned offset"); break;
                case ESPIPE:    ret_error.assign("file descripor is associated with a pipe, socket, or FIFO"); break;
                default:        ret_error.assign("unknown error"); break;
            }
            DEBUG_LOG << "File::actual_file_size(): "
                       << "error retrieving file size "
                       << "(" << ret_error << ")"
                       << " on file '" << m_file_path << "'";
            return 0;
        }

        return new_pos;
    }

protected:
    mi::Sint32              m_file_handle;

    mi::Uint64              m_position;

    std::string             m_file_path;
    mi::Uint64              m_file_size;

    std::ios_base::openmode m_open_mode;

}; // class File::File_impl

#endif // _WIN32

/// Gzip compressed file implementation.
class File_gzip_impl : public File::File_impl
{
public:
    File_gzip_impl(mi::Uint32 compression_level)
      : File_impl(),
        m_compression_level(compression_level),
        m_gz_handle(Z_NULL)
    {
    }

    virtual ~File_gzip_impl()
    {
    }

    virtual bool open(const std::string&       file_path,
                      std::ios_base::openmode  open_mode)
    {
        std::string open_flags;
        if (open_mode == std::ios_base::in)
        {
            open_flags = "rb";
        }
        else if (open_mode == std::ios_base::out)
        {
            open_flags = "wb";
        }
        else
        {
            ERROR_LOG << "File::open(): "
                      << "illegal open mode 0x" << std::hex << open_mode
                      << " for use with gzip on file "
                      <<  "'" << file_path << "'";
            return false;
        }

        std::ostringstream os;
        os << open_flags << m_compression_level;

        m_gz_handle = gzopen(file_path.c_str(), os.str().c_str());

        if (m_gz_handle == Z_NULL)
            return false;

        m_file_path  = file_path;
        m_open_mode  = open_mode;

        return true;
    }

    virtual bool is_open() const
    {
        return m_gz_handle != Z_NULL;
    }

    virtual void close()
    {
        if (is_open()) {
            gzclose(m_gz_handle);
        }

        m_gz_handle = Z_NULL;
        m_position    = 0;
        m_file_size   = 0;
        m_open_mode   = static_cast<std::ios_base::openmode>(0);
        m_file_path.clear();
    }

    virtual mi::Uint64 read(void*          output_buffer,
                            mi::Uint64     start_position,
                            mi::Uint64     num_bytes_to_read)
    {
        if (!is_open()) {
            DEBUG_LOG << "File::read(): "
                      << "read access on invalid file handle.";
            return 0;
        }

        if (start_position != 0)
        {
            ERROR_LOG << "Can't seek in gzip file";
            return 0;
        }

        using namespace mi;

        Uint8*      output_byte_buffer  = reinterpret_cast<Uint8*>(output_buffer);
        Uint64      bytes_read          = 0;

        if (num_bytes_to_read <= 0) {
            return 0;
        }

        int         file_bytes_read = 0;
        unsigned    file_bytes_to_read = static_cast<unsigned>(num_bytes_to_read);
        file_bytes_read = gzread(m_gz_handle, output_byte_buffer, file_bytes_to_read);

        if (file_bytes_read <= 0) {
            DEBUG_LOG << "File::read(): "
                      << "error reading from file " << m_file_path;
            return 0;
        }

        if (file_bytes_read == 0) {
            // eof
            return ~0ULL;
        }
        if (static_cast<Uint64>(file_bytes_read) <= num_bytes_to_read) {
            m_position  += file_bytes_read;
            bytes_read   = file_bytes_read;
        }
        else {
            DEBUG_LOG << "File::read(): "
                      << "unknown error reading from file " << m_file_path;
            return 0;
        }
        assert(bytes_read > 0);

        return bytes_read;
    }

    virtual mi::Uint64 write(const void*   input_buffer,
                     mi::Uint64    start_position,
                     mi::Uint64    num_bytes_to_write)
    {
        using namespace mi;

        if (!is_open()) {
            DEBUG_LOG << "File::write(): "
                      << "write access on invalid file handle.";
            return 0;
        }

        if (start_position != 0)
        {
            ERROR_LOG << "Can't seek in gzip file";
            return 0;
        }

        const Uint8*    input_byte_buffer   = reinterpret_cast<const Uint8*>(input_buffer);
        Uint64          bytes_written       = 0;

        int         file_bytes_written  = 0;
        unsigned    file_bytes_to_write = static_cast<unsigned>(num_bytes_to_write);
        file_bytes_written = gzwrite(m_gz_handle, input_byte_buffer, file_bytes_to_write);

        if (file_bytes_written <= 0) {
            DEBUG_LOG << "File::write(): "
                      << "error writing to file " << m_file_path;
            return 0;
        }

        if (static_cast<Uint64>(file_bytes_written) <= num_bytes_to_write) {
            m_position    += file_bytes_written;
            bytes_written  = file_bytes_written;
        }
        else {
            DEBUG_LOG << "File::write(): "
                      << "unknown error writing to file " << m_file_path;
            return 0;
        }

        return bytes_written;
    }

    virtual bool flush_buffers() const
    {
        if (!is_open()) {
            DEBUG_LOG << "File::flush_buffers(): "
                      << "write access on invalid file handle.";
            return false;
        }
        else {
            gzflush(m_gz_handle, Z_SYNC_FLUSH);
            return true;
        }
    }

    virtual mi::Uint64 size() const
    {
        ERROR_LOG << "size() no implemented for gzip";
        return m_file_size;
    }

    virtual const std::string& file_path() const
    {
        return m_file_path;
    }

    virtual mi::Uint64 actual_file_size() const
    {
        ERROR_LOG << "actual_file_size() no implemented for gzip";
        return 0;
    }

protected:
    mi::Uint32              m_compression_level;
    gzFile                  m_gz_handle;

}; // class File::File_gzip_impl

File::File()
{
    m_file_impl = new File_impl();
}

File::File(mi::Uint32 gzip_compression_level)
{
    if (gzip_compression_level == 0)
    {
        m_file_impl = new File_impl();
        return;
    }

    m_file_impl = new File_gzip_impl(gzip_compression_level);
}

File::File(const std::string&       file_path,
           std::ios_base::openmode  open_mode)
{
    m_file_impl = new File_impl();

    open(file_path, open_mode);
}

File::~File()
{
    assert(0 != m_file_impl);

    close();
    delete m_file_impl;
}

bool File::open(const std::string&       file_path,
                std::ios_base::openmode  open_mode)
{
    assert(0 != m_file_impl);
    return m_file_impl->open(file_path, open_mode);
} 

bool File::is_open() const
{
    assert(0 != m_file_impl);
    return m_file_impl->is_open();
}

void File::close()
{
    assert(0 != m_file_impl);
    return m_file_impl->close();
}

mi::Uint64 File::read(void*          output_buffer,
                      mi::Uint64     start_position,
                      mi::Uint64     num_bytes_to_read)
{
    assert(0 != m_file_impl);
    return m_file_impl->read(output_buffer, start_position, num_bytes_to_read);
}

mi::Uint64 File::write(const void*   input_buffer,
                       mi::Uint64    start_position,
                       mi::Uint64    num_bytes_to_write)
{
    assert(0 != m_file_impl);
    return m_file_impl->write(input_buffer, start_position, num_bytes_to_write);
}

bool File::flush_buffers() const
{
    assert(0 != m_file_impl);
    return m_file_impl->flush_buffers();
}

mi::Uint64 File::size() const
{
    assert(0 != m_file_impl);
    return m_file_impl->size();
}

const std::string& File::file_path() const
{
    assert(0 != m_file_impl);
    return m_file_impl->file_path();
}

File::operator bool() const
{
    return is_open();
}

} // namespace io
}} // namespace nv::index_common
