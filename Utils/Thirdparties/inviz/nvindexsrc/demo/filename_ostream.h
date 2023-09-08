/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief a simple stream associated with a filename

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_FILENAME_STREAM_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_FILENAME_STREAM_H

#include <cassert>
#include <iostream>
#include <fstream>
#include <string>

#include <sys/stat.h>

#include "common/forwarding_logger.h"

//======================================================================
/// An output stream associated with a filename.
/// A common usage is for a logger.
class Filename_ostream
{
public:
    /// Check if the path file exist. os.path.isfile() in python.
    /// 
    /// \param[in] path path string
    /// \return true when the path exists.
    static bool isfile(const std::string & path)
    {
        struct stat fileinfo;
        const mi::Sint32 ret_val = stat(path.c_str(), &fileinfo);
        if(ret_val == 0){
            return true;
        }
        return false;
    }


public:
    /// constructor
    Filename_ostream()
        :
        m_filename("/dev/null"),
        m_os(0),
        m_is_open(false),
        m_is_stream_file_on_disk(false)
    {
        // empty
    }

    /// constructor
    /// \param[in] fname output filename
    explicit Filename_ostream(const std::string & fname)
        :
        m_filename(""),
        m_os(0),
        m_is_open(false),
        m_is_stream_file_on_disk(false)
    {
        this->set_filename(fname);
        this->open_stream();
    }

    /// destructor
    virtual ~Filename_ostream()
    {
        this->clear();
    }

    /// clear the state. The stream will be closed if it is opened. 
    virtual void clear()
    {
        this->close_stream();
        m_filename = "/dev/null";
        m_is_stream_file_on_disk = true;
    }

    /// set file name. If there is open stream, the opened file will be closed.
    ///
    /// Two special filename.
    /// - stdout:    ... standard output
    /// - /dev/null  ... no output open
    /// 
    /// \param[in] fname performance value output file name
    virtual void set_filename(const std::string& fname)
    {
        this->clear();
        m_filename = fname;
        // stdout: is not a file on a disk
        m_is_stream_file_on_disk = (m_filename == "stdout:") ? false : true;
    }

    /// get performance log file names
    ///
    /// \return performance value output file name
    virtual std::string get_filename() const
    {
        return m_filename;
    }

    /// open the output stream 
    ///
    /// \return true when open succeeded
    virtual bool open_stream()
    {
        if(m_filename == "stdout:"){
            m_os = &(std::cout);
            m_is_open = true;
            return true;
        }

        m_os = new std::ofstream(m_filename.c_str(), std::ios_base::out);
        if(!(m_os->good())){
            ERROR_LOG << "Fail to open a output file [" << m_filename << "]";
            return false;
        }
        m_is_open = true;
        return true;
    }

    /// reopen the output stream. Open the file with the append mode.
    ///
    /// \return true when open succeeded
    virtual bool reopen_stream()
    {
        if(m_filename == "stdout:"){
            m_os = &(std::cout);
            m_is_open = true;
            return true;
        }

        // we can open "/dev/null"
        m_os = new std::ofstream(m_filename.c_str(), std::ios_base::out|std::ios_base::app);
        if(!(m_os->good())){
            ERROR_LOG << "Fail to reopen a output file [" << m_filename << "]";
            return false;
        }
        m_is_open = true;
        return true;
    }

    /// close the stream
    ///
    /// \return true when open succeeded
    virtual bool close_stream()
    {
        // handle special files
        if(m_filename == "stdout:"){
            m_os = 0;
            m_is_open = false;
            return true;
        }

        // not yet opened
        if(m_os == 0){
            assert(m_is_open == false);
            return false;       // fail
        }

        std::ofstream * ofs = dynamic_cast<std::ofstream *>(m_os);
        assert(ofs != 0);
        if(ofs->is_open()){
            ofs->close();

            delete m_os;
            m_os = 0;
            m_is_open = false;
            return true;
        }

        // fail to close
        m_is_open = false;
        return false;
    }

    /// state good()
    ///
    /// \return true when we can use get_os() to write.
    virtual bool good() const
    {
        if((m_os != 0) && (m_os->good()) && (m_is_open == true)){
            return true;
        }
        return false;
    }

    /// get the current stream reference
    ///
    /// \return current stream reference, may 0
    virtual std::ostream * get_os()
    {
        assert(this->good());
        return m_os;
    }

    /// is the stream a file on disk?
    /// \return true when the stream is a file on disk
    virtual bool is_stream_file_on_disk() const 
    {
        return m_is_stream_file_on_disk;
    }


private:
    /// stream file name
    std::string m_filename;
    /// reference to the current stream
    std::ostream * m_os;
    /// is open?
    bool m_is_open;
    /// is this stream a file on a disk?
    bool m_is_stream_file_on_disk;
    
private:
    /// copy constructor. prohibit until proved useful.
    Filename_ostream(Filename_ostream const &);
    /// operator=. prohibit until proved useful.
    Filename_ostream const & operator=(Filename_ostream const &);
};
//======================================================================


#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_FILENAME_STREAM_H
