/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Implementation of a forwarding message logger

#ifndef NVIDIA_INDEX_BIN_COMMON_FORWARDING_LOGGER_H
#define NVIDIA_INDEX_BIN_COMMON_FORWARDING_LOGGER_H

#include <sstream>
#include <string>

#include <mi/dice.h>
#include <mi/base/ilogger.h>
#include <mi/base/interface_implement.h>
#include <nv/index/iindex.h>

namespace nv {
namespace index_common {

/// Implementation of a forwarding message logger.
///
/// This is a container of forwarding logger from the dice. If the
/// dice is valid (started) and the IndeX interface is registered
/// in the Forwarding_logger_factory, we use the dice forwarding
/// logger. However, is dice is not present, just use a fall back
/// logger which outputs the message to stdout.
class Forwarding_logger
{
public:
    /// Constructor
    Forwarding_logger();

    /// Destructor. Message out.
    virtual ~Forwarding_logger();

    /// Get message stream
    /// \param[in] level message severity level
    /// \return output stream with the level
    std::ostringstream& get_message(mi::base::Message_severity level);

    /// Level's string representation
    ///
    /// \param[in] level severity level
    /// \return level string
    static std::string level_to_string(mi::base::Message_severity level);

private:
    /// output stream
    std::ostringstream m_os;
    /// Message severity level
    mi::base::Message_severity m_level;
    /// forwarding logger
    mi::base::Handle< mi::base::ILogger > m_forwarding_logger;

private:
    Forwarding_logger(const Forwarding_logger&);               // forbid this
    Forwarding_logger& operator=(const Forwarding_logger&);    // forbid this
};


/// Forwarding_logger factory singleton.
class Forwarding_logger_factory
{
public:
    /// get the instance
    /// \return forwarding logger singleton instance
    static Forwarding_logger_factory * instance()
    {
        if(G_p_forwarding_logger_factory == 0){
            G_p_forwarding_logger_factory = new Forwarding_logger_factory();
        }
        return G_p_forwarding_logger_factory;
    }

    /// delete the singleton. (For unit test purpose. I don't want to
    /// confuse a memory checker.)
    static void delete_instance()
    {
        if(G_p_forwarding_logger_factory != 0){
            delete G_p_forwarding_logger_factory;
            G_p_forwarding_logger_factory = 0;
        }
    }

private:
    // singleton instance
    static Forwarding_logger_factory * G_p_forwarding_logger_factory;

public:
    /// constructor
    Forwarding_logger_factory();
    /// destructor
    ~Forwarding_logger_factory();

    /// initialize the factory
    /// initialize the factory with IndeX interface. This should be
    /// only once done.
    /// \param[in] iindex_if a IIndex interface
    void initialize(mi::base::Handle<nv::index::IIndex>& iindex_if);

    /// set a message header.
    /// If the message header is not empty, it is shown up in the log.
    /// \param[in] header_str message header string to be set.
    void set_message_header(std::string const & header_str);

    /// get the current header.
    /// \return current message header string.
    std::string get_message_header() const;

    /// shutdown the factory
    /// \return true when succeeded.
    bool shutdown();

    /// is enabled the factory
    /// \return true when this factory is ready to use
    bool is_enabled() const;

    /// set fall back (= dice is not available) log severity
    /// \param[in] fb_level message output severity level
    /// \return true when succeeded
    bool set_fallback_log_severity(mi::Uint32 fb_level);

    /// get fall back log severity
    /// \return current level message output severity level
    mi::Uint32 get_fallback_log_severity() const;

    /// get a forwarding logger.
    ///
    /// This need to pass to a Handle, otherwise memory leak. Since
    /// the pointer is retained.
    /// \return a retained logger
    mi::base::ILogger * get_forwarding_logger() const;

private:
    /// IndeX interface to access to the forwarding logger
    mi::base::Handle<nv::index::IIndex> m_iindex_if;
    /// fall back message severity level when IndeX interface is not available.
    mi::Uint32 m_fallback_severity_level;
    /// message header string
    std::string m_header_str;

private:
    Forwarding_logger_factory(const Forwarding_logger_factory&);               // forbid this
    Forwarding_logger_factory& operator=(const Forwarding_logger_factory&);    // forbid this
};

}} // namespace nv::index_common

//
// Overload the stream operator for some common types.
//

// To be able to call the operators from any other namespace, they need to be defined in the same
// namespace where their argument type is defined (argument-dependend name lookup / Koenig lookup).
namespace mi {
namespace math {

/// overloaded ostream for Vector<*,2>
/// \param[in] str ostream
/// \param[in] vec output vector.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Vector<T, 2>& vec)
{
    return (str << "[" << vec.x << " " << vec.y << "]");
}

/// overloaded ostream for Vector<*,3>
/// \param[in] str ostream
/// \param[in] vec output vector.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Vector<T, 3>& vec)
{
    return (str << "[" << vec.x << " " << vec.y << " " << vec.z << "]");
}

/// overloaded ostream for Vector<*,4>
/// \param[in] str ostream
/// \param[in] vec output vector.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Vector<T, 4>& vec)
{
    return (str << "[" << vec.x << " " << vec.y << " " << vec.z << " " << vec.w << "]");
}

/// overloaded ostream for color
/// \param[in] str ostream
/// \param[in] col output color.
/// \return ostream
inline std::ostream& operator<< (std::ostream& str, const mi::math::Color& col)
{
    return (str << "[" << col.r << " " << col.g << " " << col.b << " " << col.a << "]");
}

/// overloaded ostream for Bbox<*,2>
/// \param[in] str ostream
/// \param[in] bbox output bounding box.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Bbox<T, 2>& bbox)
{
    return (str << "["
            <<  bbox.min.x << " " <<  bbox.min.y
            << "; "
            <<  bbox.max.x << " " <<  bbox.max.y
            << "]");
}

/// overloaded ostream for Bbox<*,3>
/// \param[in] str ostream
/// \param[in] bbox output bounding box.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Bbox<T, 3>& bbox)
{
    return (str << "["
            <<  bbox.min.x << " " <<  bbox.min.y << " " <<  bbox.min.z << "; "
            <<  bbox.max.x << " " <<  bbox.max.y << " " <<  bbox.max.z << "]");
}

/// overloaded ostream for Matrix<*,3,3>
/// \param[in] str ostream
/// \param[in] mat output a 3x3 matrix.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Matrix<T, 3, 3>& mat)
{
    return (str
            << "\t[ " << mat.xx << "\t " << mat.yx << "\t " << mat.zx << "\n"
            << "\t  " << mat.xy << "\t " << mat.yy << "\t " << mat.zy << "\n"
            << "\t  " << mat.xz << "\t " << mat.yz << "\t " << mat.zz << "\t ]");
}

/// overloaded ostream for Matrix<*,4,4>
/// \param[in] str ostream
/// \param[in] mat output a 4x4 matrix.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Matrix<T, 4, 4>& mat)
{
    return (str
            << "\t[ " << mat.xx << "\t " << mat.yx << "\t " << mat.zx << "\t " << mat.wx << "\n"
            << "\t  " << mat.xy << "\t " << mat.yy << "\t " << mat.zy << "\t " << mat.wy << "\n"
            << "\t  " << mat.xz << "\t " << mat.yz << "\t " << mat.zz << "\t " << mat.wz << "\n"
            << "\t  " << mat.xw << "\t " << mat.yw << "\t " << mat.zw << "\t " << mat.ww << "\t ]");
}

/// overloaded ostream for Vector_struct<*,2>
/// \param[in] str ostream
/// \param[in] vec output vector.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Vector_struct<T, 2>& vec)
{
    return (str << mi::math::Vector<T, 2>(vec));
}

/// overloaded ostream for Vector_struct<*,3>
/// \param[in] str ostream
/// \param[in] vec output vector.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Vector_struct<T, 3>& vec)
{
    return (str << mi::math::Vector<T, 3>(vec));
}

/// overloaded ostream for Vector_struct<*,4>
/// \param[in] str ostream
/// \param[in] vec output vector.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Vector_struct<T, 4>& vec)
{
    return (str << mi::math::Vector<T, 4>(vec));
}

/// overloaded ostream for Color_struct
/// \param[in] str ostream
/// \param[in] col output color.
/// \return ostream
inline std::ostream& operator<< (std::ostream& str, const mi::math::Color_struct& col)
{
    return (str << mi::math::Color(col));
}

/// overloaded ostream for Bbox_struct<*,2>
/// \param[in] str ostream
/// \param[in] bbox output bounding box.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Bbox_struct<T, 2>& bbox)
{
    return (str << mi::math::Bbox<T, 2>(bbox));
}

/// overloaded ostream for Bbox_struct<*,3>
/// \param[in] str ostream
/// \param[in] bbox output bounding box.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Bbox_struct<T, 3>& bbox)
{
    return (str << mi::math::Bbox<T, 3>(bbox));
}

/// overloaded ostream for Matrix_struct<*,3,3>
/// \param[in] str ostream
/// \param[in] mat output a 3x3 matrix.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Matrix_struct<T, 3, 3>& mat)
{
    return (str << mi::math::Matrix<T, 3, 3>(mat));
}

/// overloaded ostream for Matrix_struct<*,4,4>
/// \param[in] str ostream
/// \param[in] mat output a 4x4 matrix.
/// \return ostream
template<typename T>
inline std::ostream& operator<< (std::ostream& str, const mi::math::Matrix_struct<T, 4, 4>& mat)
{
    return (str << mi::math::Matrix<T, 4, 4>(mat));
}

} // namespace mi::math

namespace neuraylib {

/// overloaded ostream for Tag
/// \param[in] str ostream
/// \param[in] tag tag to output.
/// \return ostream
inline std::ostream& operator<< (std::ostream& str, const mi::neuraylib::Tag& tag)
{
    return (str << tag.id);
}

} // namespace mi::neuraylib
} // namespace mi

/// Log output stream for debug messages, will only be used in debug builds
#ifdef DEBUG
#define DEBUG_LOG nv::index_common::Forwarding_logger().get_message(mi::base::MESSAGE_SEVERITY_DEBUG)
#else
// Do not only disable output but also prevent evaluation of stream arguments in release mode
#define DEBUG_LOG while (false) nv::index_common::Forwarding_logger().get_message(mi::base::MESSAGE_SEVERITY_DEBUG)
#endif // DEBUG

/// Log output stream for verbose information
#define VERBOSE_LOG nv::index_common::Forwarding_logger().get_message(mi::base::MESSAGE_SEVERITY_VERBOSE)
/// Log output stream for general information
#define INFO_LOG  nv::index_common::Forwarding_logger().get_message(mi::base::MESSAGE_SEVERITY_INFO)
/// Log output stream for warnings
#define WARN_LOG  nv::index_common::Forwarding_logger().get_message(mi::base::MESSAGE_SEVERITY_WARNING)
/// Log output stream for errors
#define ERROR_LOG nv::index_common::Forwarding_logger().get_message(mi::base::MESSAGE_SEVERITY_ERROR)

#endif // NVIDIA_INDEX_BIN_COMMON_FORWARDING_LOGGER_H
