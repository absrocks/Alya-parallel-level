/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "opengl_offscreen_context.h"

#include "examiner_manipulator.h"
#include "nvindex_appdata.h"
#include "opengl_appdata.h"

#include "common/forwarding_logger.h"

#if defined (LINUX)

Offscreen_context::Offscreen_context(Display *dpy, int width, int height)
:   m_display(dpy),
    m_context(0),
    m_pbuffer(0),
    m_fbo(0),
    m_color(0),
    m_depth(0),
    m_width(width),
    m_height(height)
{
    if (!m_display) {
        return;
    }

    if((width <= 0) || (height <= 0)){
        ERROR_LOG << "Illegal window size [" << width << "," << height << "]";
        return;
    }

    // Check for GLX.
    mi::Sint32 errorBase;
    mi::Sint32 eventBase;
    if(!glXQueryExtension(m_display, &errorBase, &eventBase)) {
        ERROR_LOG << "No glX extension detected.";
        return;
    }

    mi::Sint32 attrib[] =
    {
        GLX_RENDER_TYPE,   GLX_RGBA_BIT,
        GLX_DRAWABLE_TYPE, GLX_PBUFFER_BIT,
        GLX_CONFIG_CAVEAT, GLX_NONE,
        GLX_DOUBLEBUFFER,  False,
        GLX_RED_SIZE,      8,
        GLX_GREEN_SIZE,    8,
        GLX_BLUE_SIZE,     8,
        GLX_ALPHA_SIZE,    8,
        GLX_DEPTH_SIZE,    24,
        GLX_STENCIL_SIZE,  8,
        None
    };

    mi::Sint32 nitems = 0;
    GLXFBConfig* fbconfig = glXChooseFBConfig(
        m_display, DefaultScreen(m_display), attrib, &nitems);
    if (!fbconfig || nitems == 0) {
        return;
    }

    const mi::Sint32 pbuf_attr_list[] = {
        GLX_PBUFFER_WIDTH, 1,  // just dummy values
        GLX_PBUFFER_HEIGHT, 1, // just dummy values
        GLX_LARGEST_PBUFFER, False,
        GLX_PRESERVED_CONTENTS, False,
        None
    };

    m_pbuffer = glXCreatePbuffer(m_display, *fbconfig, pbuf_attr_list);
    if (!m_pbuffer) {
	if (fbconfig) XFree(fbconfig);
	return;
    }

    m_context = glXCreateNewContext(
        m_display, *fbconfig, GLX_RGBA_TYPE,
        0, // no share
        GL_TRUE); // direct

    if (!m_context) {
	if (fbconfig) XFree(fbconfig);
        return;
    }

    if (glXMakeContextCurrent(m_display, m_pbuffer, m_pbuffer, m_context) == False) {
	if (fbconfig) XFree(fbconfig);
        return;
    }

    if (fbconfig != 0)
    {
        XFree(fbconfig);
    }

    gl_initialize_glew();

    //
    // An affinity context does not have a default framebuffer, so we need to make one
    //
    glGenFramebuffers(1, &m_fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);

    glGenRenderbuffers(1, &m_color);
    glBindRenderbuffer(GL_RENDERBUFFER, m_color);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);

    glGenRenderbuffers(1, &m_depth);
    glBindRenderbuffer(GL_RENDERBUFFER, m_depth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);

    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, m_color);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT , GL_RENDERBUFFER, m_depth);
}


Offscreen_context::~Offscreen_context()
{
    if (m_fbo)      glDeleteFramebuffers(1, &m_fbo);
    if (m_color)    glDeleteRenderbuffers(1, &m_color);
    if (m_depth)    glDeleteRenderbuffers(1, &m_depth);

    if (m_context) {
        glXMakeContextCurrent(m_display, None, None, 0);
        glXDestroyContext(m_display, m_context);
    }
    if (m_pbuffer) {
        glXDestroyPbuffer(m_display, m_pbuffer);
    }
}

void Offscreen_context::resize(mi::Sint32 width, mi::Sint32 height)
{
    if((width <= 0) || (height <= 0)){
        ERROR_LOG << "Offscreen_context::resize: Illegal window size [" << width << "," << height << "]";
        return;
    }

    if((m_width == width) && (m_height == height)){
        return;                 // same size, no update needed
    }

    // updated the current size
    m_width  = width;
    m_height = height;

    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);

    glBindRenderbuffer(GL_RENDERBUFFER, m_color);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);

    glBindRenderbuffer(GL_RENDERBUFFER, m_depth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);

    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, m_color);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT , GL_RENDERBUFFER, m_depth);
}

bool Offscreen_context::make_current()
{
    if (m_display && m_pbuffer) {
	if (glXMakeContextCurrent(m_display, m_pbuffer, m_pbuffer, m_context) == False) {
	    return false;
	}
	return true;
    }
    return false;
}


static inline mi::Sint32 GetNvXScreen(Display *dpy)
{
    mi::Sint32 defaultScreen, screen;
    defaultScreen = DefaultScreen(dpy);

    if (XNVCTRLIsNvScreen(dpy, defaultScreen)) {
        return defaultScreen;
    }

    for (screen = 0; screen < ScreenCount(dpy); screen++) {
        if (XNVCTRLIsNvScreen(dpy, screen)) {
            return screen;
        }
    }
    return -1;
}

/// initialize offscreen context for linux
void gl_initialize_offscreen_context_linux()
{
    MI::GPU::Scoped_error_handler x_error_handler;

    // we assume there is a display 0.
    bool use_fast_open_display = true;
    Display* dpy = fast_XOpenDisplay(":0");
    if (!dpy) {
        ERROR_LOG << "The connection to display 0 was not possible, trying again ....";
        dpy = XOpenDisplay(":0");
        if (!dpy) {
            ERROR_LOG << "The connection to display 0 was not possible, aborting.";
            ERROR_LOG << "Note: exit status of this process is 0 as success.";
            Nvindex_AppData::instance()->delete_instance();
            exit(0);
        }
        use_fast_open_display = false;
    }

    mi::Sint32 first_nvscreen = GetNvXScreen(dpy);
    if (first_nvscreen < 0) {
        ERROR_LOG << "No nv-screens found in display 0.";
        XCloseDisplay(dpy);
        Nvindex_AppData::instance()->delete_instance();
        exit(1);
    }
    XCloseDisplay(dpy);

    std::ostringstream oss;
    oss << ":0." << first_nvscreen;

    if(use_fast_open_display)
        dpy = fast_XOpenDisplay(oss.str().c_str());
    else
        dpy = XOpenDisplay(oss.str().c_str());

    if (!dpy) {
        ERROR_LOG << "Failed to open display " << oss.str();
        Nvindex_AppData::instance()->delete_instance();
        exit(1);
    }

    const mi::math::Vector<mi::Sint32,2> main_window_resolution = 
        Nvindex_AppData::instance()->get_user_interaction(0)->
        get_examiner()->get_main_window_resolution();
    Offscreen_context* context = new Offscreen_context(
        dpy, main_window_resolution.x, main_window_resolution.y);
    // to delete this pointer later, keep it in OpenGL_Appdata
    Nvindex_AppData::instance()->get_opengl_appdata()->set_offscreen_context(context);

    context->make_current();
}

//======================================================================
////////////////////////////
#elif defined (WIN_NT)  ///
//////////////////////////
//======================================================================

static void print_last_error()
{
    LPVOID lpMsgBuf;
    FormatMessage(
	FORMAT_MESSAGE_ALLOCATE_BUFFER |
	FORMAT_MESSAGE_FROM_SYSTEM |
	FORMAT_MESSAGE_IGNORE_INSERTS,
	NULL,
	GetLastError(),
	0, // Default language
	(LPTSTR) &lpMsgBuf,
	0,
	NULL
	);
    ERROR_LOG << (LPCTSTR)lpMsgBuf;
}


class Mini_dc_rc
{
  public:
    Mini_dc_rc();
    ~Mini_dc_rc();

  private:
    HINSTANCE	m_hinstance;	// handle to the application
    HWND	m_hwnd;		// window handle
    HDC         m_hdc;          // Window device context
    HGLRC       m_hrc;          // GL render context of window
    mi::Sint32  m_pixelformat;  // index of pixel format

    static const char* m_wnd_class_name;
    static LONG WINAPI main_win_proc(
	HWND    hwnd,		// window handle
	UINT    umsg,		// message id
	WPARAM  wparam,		// wparam
	LPARAM  lparam);	// lparam
};


//
// minimal window event handling routine: Ignore all events.
//
LONG WINAPI Mini_dc_rc::main_win_proc(
    HWND    hwnd,		// window handle
    UINT    umsg,		// message id
    WPARAM  wparam,		// wparam
    LPARAM  lparam)		// lparam
{
    LONG retval = 0;
    switch (umsg)
    {
    case WM_CREATE:
    case WM_CLOSE:
    case WM_DESTROY:
    case WM_PAINT:
    case WM_MOUSEMOVE:
    case WM_LBUTTONUP:
    case WM_SIZE: break;

    default:
	retval = (LONG) DefWindowProc(hwnd, umsg, wparam, lparam);
    }
    return retval;
}


//
// the window class name, needed to register the window class
//
const char* Mini_dc_rc::m_wnd_class_name = "Mini_dc_rc";


//
// constructor, creates window
//
Mini_dc_rc::Mini_dc_rc()
  : m_hinstance(GetModuleHandle(NULL)),
    m_hwnd(0),
    m_hdc(0),
    m_hrc(0),
    m_pixelformat(0)
{
    WNDCLASS    wc;
    memset(&wc, 0, sizeof(wc));

    // Create and register the window class only if it does not exist
    if (!GetClassInfo(m_hinstance, m_wnd_class_name, &wc)) {
	// Application specific
	wc.style	    = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
	wc.lpfnWndProc	    = (WNDPROC) main_win_proc;
	wc.cbClsExtra	    = 0;
	wc.cbWndExtra	    = 0;
	wc.hInstance	    = m_hinstance;
	wc.hIcon	    = LoadIcon(NULL, IDI_WINLOGO);
	wc.hCursor	    = LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground    = NULL;
	wc.lpszMenuName	    = NULL;
	wc.lpszClassName    = m_wnd_class_name;
	if (!RegisterClass(&wc)) {
	    return;
	}
    }

    // Window specific
    DWORD dwExStyle = WS_EX_TOPMOST | WS_EX_APPWINDOW | WS_EX_WINDOWEDGE;
    DWORD dwStyle   = WS_OVERLAPPEDWINDOW;

    RECT win_rect;
    win_rect.left   = 0;
    win_rect.right  = 1;
    win_rect.top    = 0;
    win_rect.bottom = 1;
    // Adjust Window To True Requested Size
    AdjustWindowRectEx(&win_rect, dwStyle, FALSE, dwExStyle);

    // Create The Window
    m_hwnd = CreateWindowEx(
	dwExStyle,			// Extended Style For The Window
	m_wnd_class_name,		// Class Name
	m_wnd_class_name,		// Window Title
	dwStyle |			// Defined Window Style
	WS_CLIPSIBLINGS |		// Required Window Style
	WS_CLIPCHILDREN,		// Required Window Style
	0, 0,				// Window Position
	1,                              // Window Width
        1,                              // Window Height
	NULL,				// No Parent Window
	NULL,				// No Menu
	m_hinstance,			// Instance
	NULL);	 			// Don't Pass Anything To WM_CREATE


    // get device context
    if (!(m_hdc = GetDC(m_hwnd))) {
	return;
    }

    // Description of the pixel
    PIXELFORMATDESCRIPTOR pfd;
    memset(&pfd, 0, sizeof(pfd));
    pfd.nSize        = sizeof(pfd);
    pfd.nVersion     = 1;
    pfd.dwFlags      = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL;
    pfd.iPixelType   = PFD_TYPE_RGBA;
    pfd.cColorBits   = 0;
    pfd.cAlphaBits   = 0;
    pfd.cStencilBits = 0;
    pfd.cDepthBits   = 0;
    pfd.iLayerType   = PFD_MAIN_PLANE;


    // This is to find the first standard pixel format
    m_pixelformat = ChoosePixelFormat(m_hdc, &pfd);
    if (!m_pixelformat) {
	return;
    }

    // Are we able to set the pixel format ?
    if (!SetPixelFormat(m_hdc, m_pixelformat, &pfd)) {
	return;
    }

    // Are we able to get a rendering context?
    if (!(m_hrc = wglCreateContext(m_hdc))) {
	return;
    }

    // Try to activate the rendering context
    if (!wglMakeCurrent(m_hdc, m_hrc)) {
	return;
    }

}


//
// destructor, destroys window and unregisters window class
Mini_dc_rc::~Mini_dc_rc()
{
    if (m_hrc) {
	wglMakeCurrent(m_hdc, NULL);
        wglDeleteContext(m_hrc);
	m_hrc = 0;
    }

    if (m_hdc) {
	ReleaseDC(m_hwnd, m_hdc);
        m_hdc = 0;
    }

    if (m_hwnd) {
	DestroyWindow(m_hwnd);
	m_hwnd = 0;
    }

    if (m_hinstance) {
	UnregisterClass(m_wnd_class_name, m_hinstance);
	m_hinstance = 0;
    }
}

PFNWGLCREATEAFFINITYDCNVPROC        Offscreen_context::mi_wglCreateAffinityDCNV          = NULL;
PFNWGLDELETEDCNVPROC                Offscreen_context::mi_wglDeleteDCNV                  = NULL;
PFNWGLENUMGPUDEVICESNVPROC          Offscreen_context::mi_wglEnumGpuDevicesNV            = NULL;
PFNWGLENUMGPUSFROMAFFINITYDCNVPROC  Offscreen_context::mi_wglEnumGpusFromAffinityDCNV    = NULL;
PFNWGLENUMGPUSNVPROC                Offscreen_context::mi_wglEnumGpusNV                  = NULL;

Offscreen_context::Offscreen_context(mi::Sint32 gpu_index, mi::Sint32 width, mi::Sint32 height)
{
    m_hdc	    = 0;
    m_hrc	    = 0;
    m_pix_format    = 0;
    m_fbo           = 0;
    m_color         = 0;
    m_depth         = 0;
    m_width         = width;
    m_height        = height;

    {
        Mini_dc_rc mini_dc_rc;
        mi_wglCreateAffinityDCNV       = (PFNWGLCREATEAFFINITYDCNVPROC)
            wglGetProcAddress((LPCSTR)"wglCreateAffinityDCNV");
        mi_wglDeleteDCNV               = (PFNWGLDELETEDCNVPROC)
            wglGetProcAddress((LPCSTR)"wglDeleteDCNV");
        mi_wglEnumGpuDevicesNV         = (PFNWGLENUMGPUDEVICESNVPROC)
            wglGetProcAddress((LPCSTR)"wglEnumGpuDevicesNV");
        mi_wglEnumGpusFromAffinityDCNV = (PFNWGLENUMGPUSFROMAFFINITYDCNVPROC)
            wglGetProcAddress((LPCSTR)"wglEnumGpusFromAffinityDCNV");
        mi_wglEnumGpusNV               = (PFNWGLENUMGPUSNVPROC)
            wglGetProcAddress((LPCSTR)"wglEnumGpusNV");
    }

    // Try to setup an affinity dc
    if (mi_wglEnumGpusNV != NULL) {
	HGPUNV     gpu_handle=0;
        GPU_DEVICE gpu_device;
        memset(&gpu_device, 0, sizeof gpu_device);
        gpu_device.cb = sizeof gpu_device;

        if (!mi_wglEnumGpusNV(gpu_index, &gpu_handle)) {
            ERROR_LOG << "wglEnumGpusNV failed";
            return;
        }
        if (!mi_wglEnumGpuDevicesNV(gpu_handle, 0, &gpu_device)) {
            ERROR_LOG << "wglEnumGpuDevicesNV failed";
            return;
        }

        HGPUNV gpu_mask[2] = {gpu_handle, 0};
	m_hdc = mi_wglCreateAffinityDCNV(gpu_mask);
    }
    else {
        ERROR_LOG << "wglEnumGpusNV " << gpu_index << " is missing (it is available on Quadro cards).";
        return;
    }

    if (!m_hdc) {
        ERROR_LOG << "Creating an affinity OpenGL context on GPU " << gpu_index << " failed.";
        return;
    }

    // Description of the pixel
    PIXELFORMATDESCRIPTOR pfd;
    memset(&pfd, 0, sizeof(pfd));
    pfd.nSize        = sizeof(pfd);
    pfd.nVersion     = 1;
    pfd.dwFlags      = PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL;
    pfd.iPixelType   = PFD_TYPE_RGBA;
    pfd.cColorBits   = 0;
    pfd.cAlphaBits   = 0;
    pfd.cStencilBits = 0;
    pfd.cDepthBits   = 0;
    pfd.iLayerType   = PFD_MAIN_PLANE;


    // This is to find the first standard pixel format
    m_pix_format = ChoosePixelFormat(m_hdc, &pfd);
    DescribePixelFormat(m_hdc, m_pix_format, sizeof(pfd), &pfd);

    if (!SetPixelFormat(m_hdc, m_pix_format, &pfd)) {
	ERROR_LOG <<  "Can't set the PixelFormat";
	print_last_error();
	return;
    }

    if (!(m_hrc = wglCreateContext(m_hdc))) {
	ERROR_LOG <<  "Can't create a GL rendering context";
	print_last_error();
	return;
    }

    if (!wglMakeCurrent(m_hdc, m_hrc)) {
	ERROR_LOG <<  "Can't activate the GL rendering context";
	print_last_error();
	return;
    }

    gl_initialize_glew();

    //
    // An affinity context does not have a default framebuffer, so we need to make one
    //
    glGenFramebuffers(1, &m_fbo);
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);

    glGenRenderbuffers(1, &m_color);
    glBindRenderbuffer(GL_RENDERBUFFER, m_color);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);

    glGenRenderbuffers(1, &m_depth);
    glBindRenderbuffer(GL_RENDERBUFFER, m_depth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);

    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, m_color);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT , GL_RENDERBUFFER, m_depth);
}


Offscreen_context::~Offscreen_context()
{
    if (m_fbo)      glDeleteFramebuffers(1, &m_fbo);
    if (m_color)    glDeleteRenderbuffers(1, &m_color);
    if (m_depth)    glDeleteRenderbuffers(1, &m_depth);

    if (m_hrc && !wglDeleteContext(m_hrc)) {
	ERROR_LOG << "Failed to release rendering context.";
	print_last_error();
    }
    m_hrc = NULL;

    if (m_hdc && !mi_wglDeleteDCNV(m_hdc)) {
        ERROR_LOG << "Failed to release affinity device context.";
        print_last_error();
    }
    m_hdc = NULL;
}

void Offscreen_context::resize(mi::Sint32 width, mi::Sint32 height)
{
    if((width <= 0) || (height <= 0)){
        ERROR_LOG << "Offscreen_context::resize: Illegal window size [" << width << "," << height << "]";
        return;
    }

    if((m_width == width) && (m_height == height)){
        return;                 // same size, no update needed
    }

    // updated the current size
    m_width  = width;
    m_height = height;

    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);

    glBindRenderbuffer(GL_RENDERBUFFER, m_color);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, width, height);

    glBindRenderbuffer(GL_RENDERBUFFER, m_depth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, height);

    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, m_color);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT , GL_RENDERBUFFER, m_depth);
}

// Make this render context current for the thread (like OpenGL makeCurrent)
void  Offscreen_context::make_current()
{
    if (m_hdc && m_hrc) {
        if (!wglMakeCurrent(m_hdc, m_hrc)) {
	    ERROR_LOG << "Can't activate the GL rendering context";
	    print_last_error();
        }
    }
}


/// initialize offscreen context for windows
void gl_initialize_offscreen_context_win_nt()
{
    const mi::math::Vector<mi::Sint32,2> main_window_resolution = 
        Nvindex_AppData::instance()->get_user_interaction(0)->
        get_examiner()->get_main_window_resolution();

    Offscreen_context* context = new Offscreen_context(
        0, main_window_resolution.x, main_window_resolution.y);

    // to delete this pointer later, keep it in OpenGL_Appdata
    Nvindex_AppData::instance()->get_opengl_appdata()->set_offscreen_context(context);

    context->make_current();
}

#else
# error "unknown environment"
#endif  // defined (LINUX)

//----------------------------------------------------------------------
// initialize offscreen context
void gl_initialize_offscreen_context()
{
#if defined (LINUX)
    gl_initialize_offscreen_context_linux();
#elif defined (WIN_NT)
    gl_initialize_offscreen_context_win_nt();
#else
    ERROR_LOG << "No OpenGL offscreen canvas initialized. Shutting down the application, now.";
    Nvindex_AppData::instance()->delete_instance();
    exit(1);
#endif  // defined (LINUX)
}

//----------------------------------------------------------------------

