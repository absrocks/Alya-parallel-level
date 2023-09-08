/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 ******************************************************************************/
/// \file

#ifndef MI_XWIN_SCOPED_ERROR_HANDLER
#define MI_XWIN_SCOPED_ERROR_HANDLER

#include <X11/Xlib.h>

namespace MI {
namespace GPU {

/// Override the X error handler for local scope
class Scoped_error_handler
{
  public:

    Scoped_error_handler();
    ~Scoped_error_handler();

  private:
    int (* m_old_error_handler) (Display*, XErrorEvent*);
};

}}

#endif
