/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief opengl offscreen context

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_OFFSCREEN_CONTEXT_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_OFFSCREEN_CONTEXT_H

#include <GL/glew.h>

#if defined (LINUX)
#  include <GL/glx.h>
#  include "xwin_fastxopendisplay.h"
#  include "NVCtrlLib.h"
#elif defined (WIN_NT)
#  include <GL/wglew.h>
#endif

#ifndef WIN_NT
#  include <X11/Xlib.h>
#  include "xwin_scoped_error_handler.h"
#endif

#include "opengl_drawing_utilities.h" // gl_initialize_offscreen_context()


//======================================================================
#if defined (LINUX)

/// OpengGL context with Pbuffer.
class Offscreen_context
{
public:
    /// Constructor.
    explicit Offscreen_context(Display* dpy, mi::Sint32 width, mi::Sint32 height);

    /// resize the render target
    void resize(mi::Sint32 width, mi::Sint32 height);

    /// Destructor.
    ~Offscreen_context();

    /// Make context current.
    bool  make_current();

private:
    Display*    m_display;
    GLXContext  m_context;
    GLXPbuffer	m_pbuffer;

    GLuint  m_fbo;
    GLuint  m_color;
    GLuint  m_depth;

    // current off screen window size
    mi::Sint32 m_width;
    mi::Sint32 m_height;

private:
    Offscreen_context(Offscreen_context const&);
    Offscreen_context const& operator=(Offscreen_context const&);
};

//======================================================================
////////////////////////////
#elif defined (WIN_NT)  ///
//////////////////////////
//======================================================================

class Offscreen_context
{
  public:
    Offscreen_context(mi::Sint32 gpu_index, mi::Sint32 width, mi::Sint32 height);
    ~Offscreen_context();

    void resize(mi::Sint32 width, mi::Sint32 height);

    void make_current();

  private:
    HDC        m_hdc;          // Window device context
    HGLRC      m_hrc;          // GL render context of window
    mi::Sint32 m_pix_format;   // Pixel format

    GLuint  m_fbo;
    GLuint  m_color;
    GLuint  m_depth;

    // current off screen window size
    mi::Sint32 m_width;
    mi::Sint32 m_height;

    static PFNWGLCREATEAFFINITYDCNVPROC        mi_wglCreateAffinityDCNV;
    static PFNWGLDELETEDCNVPROC                mi_wglDeleteDCNV;
    static PFNWGLENUMGPUDEVICESNVPROC          mi_wglEnumGpuDevicesNV;
    static PFNWGLENUMGPUSFROMAFFINITYDCNVPROC  mi_wglEnumGpusFromAffinityDCNV;
    static PFNWGLENUMGPUSNVPROC                mi_wglEnumGpusNV;
};

//======================================================================
#else
#  error "Unknown environment"
#endif
//======================================================================

#endif // NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_OPENGL_OFFSCREEN_CONTEXT_H
