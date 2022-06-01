/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Rendering kernel program parameter buffer manipulation example implementation

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTC_PARAMETER_BUFFER_MANIP_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTC_PARAMETER_BUFFER_MANIP_H

#include <mi/dice.h>
#include <mi/base/types.h>

#include <nv/index/irendering_kernel_programs.h>

#include <cassert>

#include "common/forwarding_logger.h"


/// IndeX rendering kernel parameter buffer manipulation functionality.
///
/// Singleton
///
class RTC_parameter_buffer_manip
{
public:
    /// set heightfield workflow type
    void add_parameter_buffer_data(
            mi::base::Handle<nv::index::IRendering_kernel_program_parameters>& param_buffer,
            mi::Float32                                                        param_value = 0.0f);

    void edit_parameter_buffer_instance(
            const mi::neuraylib::Tag&         param_buffer_tag,
            const mi::Float32                 param_buffer_value,
            mi::neuraylib::IDice_transaction* dice_transaction);

public:
    /// Access workflow instance
    static RTC_parameter_buffer_manip* instance();
    /// delete the heightfield workflow instance
    static void delete_instance();

private:
    /// Singleton instance
    static RTC_parameter_buffer_manip* g_rtc_parameter_buffer_manip;

    RTC_parameter_buffer_manip();

private:
};

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RTC_PARAMETER_BUFFER_MANIP_H
