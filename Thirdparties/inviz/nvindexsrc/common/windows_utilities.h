/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Utility functions for Windows builds.

#ifndef NVIDIA_INDEX_BIN_COMMON_WINDOWS_UTILITIES_H
#define NVIDIA_INDEX_BIN_COMMON_WINDOWS_UTILITIES_H

#include <string>

#include <mi/base/types.h>

namespace nv {
namespace index_common {

// some Unix function replacements
#if _WIN32

mi::Sint32 unsetenv(const char* env_name);

#endif // _WIN32

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_WINDOWS_UTILITIES_H
