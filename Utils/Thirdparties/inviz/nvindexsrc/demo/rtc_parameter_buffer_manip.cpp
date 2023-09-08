/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "rtc_parameter_buffer_manip.h"

RTC_parameter_buffer_manip* RTC_parameter_buffer_manip::g_rtc_parameter_buffer_manip = 0;

RTC_parameter_buffer_manip::RTC_parameter_buffer_manip()
{
}

/// set heightfield workflow type
void RTC_parameter_buffer_manip::add_parameter_buffer_data(
        mi::base::Handle<nv::index::IRendering_kernel_program_parameters>& param_buffer,
        mi::Float32                                                        param_value)
{
    // testing code, this needs to match what is used in the RTC program
    typedef mi::math::Vector<mi::Float32, 3>    Vec3f;
    typedef mi::math::Vector<mi::Float32, 4>    Vec4f;

    struct Plane_data_struct
    {
        Vec4f   plane_equation;
    };

    Plane_data_struct planes[2];

    Vec3f clip_plane_normal = Vec3f(1.0f, 1.0f, 1.0f);
    clip_plane_normal.normalize();

    planes[0].plane_equation = Vec4f(clip_plane_normal,   0.0f - param_value);
    planes[1].plane_equation = Vec4f(-clip_plane_normal, 50.0f + param_value);

    param_buffer->set_buffer_data(0, planes,     sizeof(Plane_data_struct) * 2);
    param_buffer->set_buffer_data(3, planes,     sizeof(Plane_data_struct));
    param_buffer->set_buffer_data(7, planes + 1, sizeof(Plane_data_struct));
}

void RTC_parameter_buffer_manip::edit_parameter_buffer_instance(
        const mi::neuraylib::Tag&         param_buffer_tag,
        const mi::Float32                 param_buffer_value,
        mi::neuraylib::IDice_transaction* dice_transaction)
{
    if (!param_buffer_tag.is_valid())
    {
        return;
    }

    mi::base::Handle<nv::index::IRendering_kernel_program_parameters> rtc_param_buffer(
        dice_transaction->edit<nv::index::IRendering_kernel_program_parameters>(param_buffer_tag));
    
    if (rtc_param_buffer.is_valid_interface())
    {
        add_parameter_buffer_data(rtc_param_buffer, param_buffer_value);
    }
}

RTC_parameter_buffer_manip* RTC_parameter_buffer_manip::instance()
{
    if (g_rtc_parameter_buffer_manip == 0)
    {
        g_rtc_parameter_buffer_manip = new RTC_parameter_buffer_manip();
    }
    return g_rtc_parameter_buffer_manip;
}

void RTC_parameter_buffer_manip::delete_instance()
{
    if (g_rtc_parameter_buffer_manip != 0)
    {
        delete g_rtc_parameter_buffer_manip;
        g_rtc_parameter_buffer_manip = 0;
    }
}

