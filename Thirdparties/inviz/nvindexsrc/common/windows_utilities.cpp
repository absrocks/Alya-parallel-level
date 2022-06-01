/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Utility functions for Windows builds.

#include "windows_utilities.h"

#include <cassert>
#include <stdlib.h>

#if _WIN32
#   ifndef WIN32_LEAN_AND_MEAN
#       define WIN32_LEAN_AND_MEAN 1
#   endif
#   include <windows.h>
#   include <winbase.h>
#endif // _WIN32

namespace nv {
namespace index_common {

#if _WIN32

mi::Sint32 unsetenv(const char* env_name)
{
    if (SetEnvironmentVariable(env_name, NULL) != 0) {
        return 0;
    }
    else {
        return -1;
    }
}

#endif // _WIN32

}} // namespace nv::index_common
