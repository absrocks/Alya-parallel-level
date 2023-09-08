/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief progress callback application implementation

#ifndef NVIDIA_INDEX_BIN_PROGRESS_CALLBACK_H
#define NVIDIA_INDEX_BIN_PROGRESS_CALLBACK_H

#include <mi/math/vector.h>
#include <mi/math/matrix.h>
#include <mi/base/types.h>

#include <nv/index/iprogress_callback.h>

#include "common/forwarding_logger.h"


class Progress_callback : public mi::base::Interface_implement<nv::index::IProgress_callback>
{
public:
    Progress_callback()
      : m_percentage(0.f)
    {
    }
    
    virtual void set_progress(mi::Float32 progress)
    {
        m_percentage = progress;
    }

    virtual mi::Float32 get_progress() const 
    {
        return m_percentage;
    }
    
private:
    mi::Float32 m_percentage;
};


#endif // NVIDIA_INDEX_BIN_PROGRESS_CALLBACK_H
