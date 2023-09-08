/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Heightfield data type (file type, temporary)

#ifndef NVIDIA_INDEX_BIN_COMMON_RAW_HEIGHTFIELD_DATA_TYPE_H
#define NVIDIA_INDEX_BIN_COMMON_RAW_HEIGHTFIELD_DATA_TYPE_H

namespace nv {
namespace index_common {

/// FIXME: replace with Raw_heightfield_file_type inside the class. FILE_TYPE -> RHDI_NONE, ...
/// \deprecated
enum File_type
{
    /// obsolete encoded normal file
    FILE_TYPE_SSV,
    /// a demo
    FILE_TYPE_DEMO,
    /// (sentinel)
    FILE_TYPE_COUNT
};

/// Raw heightfield file format type
enum Raw_heightfield_file_format_type
{
    /// Reserved as none
    RHFFT_NONE = 0,
    /// IndeX raw heightfield file format using -1 as hole marker and only supporting positive
    /// height values
    RHFFT_INDEX_RAW_v1,
    /// IndeX raw heightfield file format using NaN as hole marker and therefore also supporing
    /// negative height values
    RHFFT_INDEX_RAW_v2,
    /// (sentinel)
    RHFFT_COUNT,
};

} // namespace index_common
} // namespace nv
#endif // NVIDIA_INDEX_BIN_COMMON_RAW_HEIGHTFIELD_DATA_TYPE_H
