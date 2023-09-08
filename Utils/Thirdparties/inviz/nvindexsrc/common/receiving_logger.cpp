/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "receiving_logger.h"

#include <iostream>

#ifdef LINUX
#include <execinfo.h>
#include <signal.h>
#include <unistd.h>
#endif // LINUX

// color cout support could be put in a separate header
#ifdef WIN_NT
#include <io.h>

#ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN 1
#endif
#include <windows.h>

namespace util {

class Default_console
{
public:
    static const WORD bg_mask = BACKGROUND_BLUE | BACKGROUND_GREEN | BACKGROUND_RED | BACKGROUND_INTENSITY;
    static const WORD fg_mask = FOREGROUND_BLUE | FOREGROUND_GREEN | FOREGROUND_RED | FOREGROUND_INTENSITY;

    static const WORD fg_black          = 0;
    static const WORD fg_low_red        = FOREGROUND_RED;
    static const WORD fg_low_green      = FOREGROUND_GREEN;
    static const WORD fg_low_blue       = FOREGROUND_BLUE;
    static const WORD fg_low_cyan       = fg_low_green   | fg_low_blue;
    static const WORD fg_low_magenta    = fg_low_red     | fg_low_blue;
    static const WORD fg_low_yellow     = fg_low_red     | fg_low_green;
    static const WORD fg_low_white      = fg_low_red     | fg_low_green | fg_low_blue;
    static const WORD fg_gray           = fg_black       | FOREGROUND_INTENSITY; 
    static const WORD fg_hi_white       = fg_low_white   | FOREGROUND_INTENSITY; 
    static const WORD fg_hi_blue        = fg_low_blue    | FOREGROUND_INTENSITY; 
    static const WORD fg_hi_green       = fg_low_green   | FOREGROUND_INTENSITY; 
    static const WORD fg_hi_red         = fg_low_red     | FOREGROUND_INTENSITY; 
    static const WORD fg_hi_cyan        = fg_low_cyan    | FOREGROUND_INTENSITY; 
    static const WORD fg_hi_magenta     = fg_low_magenta | FOREGROUND_INTENSITY; 
    static const WORD fg_hi_yellow      = fg_low_yellow  | FOREGROUND_INTENSITY;

    static const WORD bg_black          = 0;
    static const WORD bg_low_red        = BACKGROUND_RED;
    static const WORD bg_low_green      = BACKGROUND_GREEN;
    static const WORD bg_low_blue       = BACKGROUND_BLUE;
    static const WORD bg_low_cyan       = bg_low_green   | bg_low_blue;
    static const WORD bg_low_magenta    = bg_low_red     | bg_low_blue;
    static const WORD bg_low_yellow     = bg_low_red     | bg_low_green;
    static const WORD bg_low_white      = bg_low_red     | bg_low_green | bg_low_blue;
    static const WORD bg_gray           = bg_black       | BACKGROUND_INTENSITY;
    static const WORD bg_hi_white       = bg_low_white   | BACKGROUND_INTENSITY;
    static const WORD bg_hi_blue        = bg_low_blue    | BACKGROUND_INTENSITY;
    static const WORD bg_hi_green       = bg_low_green   | BACKGROUND_INTENSITY;
    static const WORD bg_hi_red         = bg_low_red     | BACKGROUND_INTENSITY;
    static const WORD bg_hi_cyan        = bg_low_cyan    | BACKGROUND_INTENSITY;
    static const WORD bg_hi_magenta     = bg_low_magenta | BACKGROUND_INTENSITY;
    static const WORD bg_hi_yellow      = bg_low_yellow  | BACKGROUND_INTENSITY;
    

private:
    HANDLE                      m_con_handle;
    DWORD                       m_con_size;
    DWORD                       m_chars_written; 
    CONSOLE_SCREEN_BUFFER_INFO  m_info; 
    CONSOLE_SCREEN_BUFFER_INFO  m_def_info; 

public:
    Default_console()
      : m_con_handle(GetStdHandle(STD_OUTPUT_HANDLE))
    { 
        GetConsoleScreenBufferInfo(m_con_handle, &m_def_info);
    }
public:
    void
    clear()
    {
        COORD screen_coord = {0, 0};
            
        get_info(); 
        FillConsoleOutputCharacter(m_con_handle, ' ', m_con_size, screen_coord, &m_chars_written); 
        get_info(); 
        FillConsoleOutputAttribute(m_con_handle, m_info.wAttributes, m_con_size, screen_coord, &m_chars_written); 
        SetConsoleCursorPosition(m_con_handle, screen_coord); 
    }
    void
    set_color(const WORD rgbi, const WORD mask)
    {
        get_info();
        m_info.wAttributes &= mask; 
        m_info.wAttributes |= rgbi; 
        SetConsoleTextAttribute(m_con_handle, m_info.wAttributes);
    }
    void
    reset_color()
    {
        SetConsoleTextAttribute(m_con_handle, m_def_info.wAttributes);
    }
private:
    void
    get_info()
    {
        GetConsoleScreenBufferInfo(m_con_handle, &m_info);
        m_con_size = m_info.dwSize.X * m_info.dwSize.Y; 
    }
}; // class Default_console

static Default_console def_console;
    
inline std::ostream& clr        (std::ostream& os) { os.flush(); def_console.clear();       return os; };
inline std::ostream& reset_color(std::ostream& os) { os.flush(); def_console.reset_color(); return os; };

inline std::ostream& fg_red         (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_red,           Default_console::bg_mask); return os; }
inline std::ostream& fg_green       (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_green,         Default_console::bg_mask); return os; }
inline std::ostream& fg_blue        (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_blue,          Default_console::bg_mask); return os; }
inline std::ostream& fg_white       (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_white,         Default_console::bg_mask); return os; }
inline std::ostream& fg_cyan        (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_cyan,          Default_console::bg_mask); return os; }
inline std::ostream& fg_magenta     (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_magenta,       Default_console::bg_mask); return os; }
inline std::ostream& fg_yellow      (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_yellow,        Default_console::bg_mask); return os; }
inline std::ostream& fg_dk_red      (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_red,          Default_console::bg_mask); return os; }
inline std::ostream& fg_dk_green    (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_green,        Default_console::bg_mask); return os; }
inline std::ostream& fg_dk_blue     (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_blue,         Default_console::bg_mask); return os; }
inline std::ostream& fg_dk_white    (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_white,        Default_console::bg_mask); return os; }
inline std::ostream& fg_dk_cyan     (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_cyan,         Default_console::bg_mask); return os; }
inline std::ostream& fg_dk_magenta  (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_magenta,      Default_console::bg_mask); return os; }
inline std::ostream& fg_dk_yellow   (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_yellow,       Default_console::bg_mask); return os; }
inline std::ostream& fg_black       (std::ostream& os) { os.flush(); def_console.set_color(Default_console::fg_black,            Default_console::bg_mask); return os; }
inline std::ostream& fg_gray        (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_gray,             Default_console::bg_mask); return os; }
inline std::ostream& bg_red         (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_red,           Default_console::fg_mask); return os; }
inline std::ostream& bg_green       (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_green,         Default_console::fg_mask); return os; }
inline std::ostream& bg_blue        (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_blue,          Default_console::fg_mask); return os; }
inline std::ostream& bg_white       (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_white,         Default_console::fg_mask); return os; }
inline std::ostream& bg_cyan        (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_cyan,          Default_console::fg_mask); return os; } 
inline std::ostream& bg_magenta     (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_magenta,       Default_console::fg_mask); return os; }
inline std::ostream& bg_yellow      (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_yellow,        Default_console::fg_mask); return os; }
inline std::ostream& bg_dk_red      (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_red,          Default_console::fg_mask); return os; }
inline std::ostream& bg_dk_green    (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_green,        Default_console::fg_mask); return os; }
inline std::ostream& bg_dk_blue     (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_blue,         Default_console::fg_mask); return os; }
inline std::ostream& bg_dk_white    (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_white,        Default_console::fg_mask); return os; }
inline std::ostream& bg_dk_cyan     (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_cyan,         Default_console::fg_mask); return os; } 
inline std::ostream& bg_dk_magenta  (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_magenta,      Default_console::fg_mask); return os; }
inline std::ostream& bg_dk_yellow   (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_yellow,       Default_console::fg_mask); return os; }
inline std::ostream& bg_black       (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_black,            Default_console::fg_mask); return os; }
inline std::ostream& bg_gray        (std::ostream& os) { os.flush(); def_console.set_color(Default_console::bg_gray,             Default_console::fg_mask); return os; }

inline std::wostream& clr        (std::wostream& os) { os.flush(); def_console.clear();       return os; };
inline std::wostream& reset_color(std::wostream& os) { os.flush(); def_console.reset_color(); return os; };

inline std::wostream& fg_red         (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_red,           Default_console::bg_mask); return os; }
inline std::wostream& fg_green       (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_green,         Default_console::bg_mask); return os; }
inline std::wostream& fg_blue        (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_blue,          Default_console::bg_mask); return os; }
inline std::wostream& fg_white       (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_white,         Default_console::bg_mask); return os; }
inline std::wostream& fg_cyan        (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_cyan,          Default_console::bg_mask); return os; }
inline std::wostream& fg_magenta     (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_magenta,       Default_console::bg_mask); return os; }
inline std::wostream& fg_yellow      (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_hi_yellow,        Default_console::bg_mask); return os; }
inline std::wostream& fg_dk_red      (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_red,          Default_console::bg_mask); return os; }
inline std::wostream& fg_dk_green    (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_green,        Default_console::bg_mask); return os; }
inline std::wostream& fg_dk_blue     (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_blue,         Default_console::bg_mask); return os; }
inline std::wostream& fg_dk_white    (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_white,        Default_console::bg_mask); return os; }
inline std::wostream& fg_dk_cyan     (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_cyan,         Default_console::bg_mask); return os; }
inline std::wostream& fg_dk_magenta  (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_magenta,      Default_console::bg_mask); return os; }
inline std::wostream& fg_dk_yellow   (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_low_yellow,       Default_console::bg_mask); return os; }
inline std::wostream& fg_black       (std::wostream& os) { os.flush(); def_console.set_color(Default_console::fg_black,            Default_console::bg_mask); return os; }
inline std::wostream& fg_gray        (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_gray,             Default_console::bg_mask); return os; }
inline std::wostream& bg_red         (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_red,           Default_console::fg_mask); return os; }
inline std::wostream& bg_green       (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_green,         Default_console::fg_mask); return os; }
inline std::wostream& bg_blue        (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_blue,          Default_console::fg_mask); return os; }
inline std::wostream& bg_white       (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_white,         Default_console::fg_mask); return os; }
inline std::wostream& bg_cyan        (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_cyan,          Default_console::fg_mask); return os; } 
inline std::wostream& bg_magenta     (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_magenta,       Default_console::fg_mask); return os; }
inline std::wostream& bg_yellow      (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_hi_yellow,        Default_console::fg_mask); return os; }
inline std::wostream& bg_dk_red      (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_red,          Default_console::fg_mask); return os; }
inline std::wostream& bg_dk_green    (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_green,        Default_console::fg_mask); return os; }
inline std::wostream& bg_dk_blue     (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_blue,         Default_console::fg_mask); return os; }
inline std::wostream& bg_dk_white    (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_white,        Default_console::fg_mask); return os; }
inline std::wostream& bg_dk_cyan     (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_cyan,         Default_console::fg_mask); return os; } 
inline std::wostream& bg_dk_magenta  (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_magenta,      Default_console::fg_mask); return os; }
inline std::wostream& bg_dk_yellow   (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_low_yellow,       Default_console::fg_mask); return os; }
inline std::wostream& bg_black       (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_black,            Default_console::fg_mask); return os; }
inline std::wostream& bg_gray        (std::wostream& os) { os.flush(); def_console.set_color(Default_console::bg_gray,             Default_console::fg_mask); return os; }

} // namespace util

#endif // WIN_NT

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
Receiving_logger::Receiving_logger()
  : m_os(),
    m_level(mi::base::MESSAGE_SEVERITY_INFO),
    m_color(false)
{
    // Only enable color output if we are connected to a terminal and not if stdout is redirected to
    // a file or connected to a pipe
#ifdef LINUX
    if (isatty(fileno(stdout)) == 1)
    {
        // Check terminal type
        std::string term;
        const char* tmp = getenv("TERM");
        if (tmp)
            term = tmp;
        // Allow color output only on supported terminal types
        if (term == "xterm" || term == "xterm-color" || term == "xterm-256color" ||
            term == "screen" || term == "vt100")
        {
            m_color = true;
        }
    }
#endif // LINUX
#ifdef WIN_NT
    if (_isatty(_fileno(stdout)) != 0)
    {
        m_color = true;
    }
#endif // WIN_NT
}

#ifdef LINUX
#ifndef DEBUG

namespace {

static int s_running_in_debugger = -1;

static void sigtrap_handler(int signum)
{
    s_running_in_debugger = 0; // we got the signal, so no debugger
    signal(SIGTRAP, SIG_DFL);  // restore default handler
}

// Will trigger a breakpoint, but only when the application is running in a debugger.
void debug_break()
{
    if (s_running_in_debugger == -1)
    {
        s_running_in_debugger = 1; // assume we are running in a debugger

        signal(SIGTRAP, sigtrap_handler); // install our own handler

        // Now trigger the signal. It will either be called by our handler or by a debugger, if
        // there is one running
        raise(SIGTRAP);
    }
}

} // namespace

#endif // DEBUG
#endif // LINUX

//----------------------------------------------------------------------
Receiving_logger::~Receiving_logger()
{
#ifndef DEBUG
    if (m_level >= mi::base::MESSAGE_SEVERITY_DEBUG) // No debug output when optimized build
        return;
#endif

    if (m_os.str().empty() && m_level > mi::base::MESSAGE_SEVERITY_FATAL)
        return;                 // no message

    // Show the message
    if (m_color && m_level <= mi::base::MESSAGE_SEVERITY_WARNING)
    {
#ifdef LINUX
        // See http://en.wikipedia.org/wiki/ANSI_escape_code#Colors
        std::string color;
        if (m_level == mi::base::MESSAGE_SEVERITY_WARNING)
            color = "\033[31m";      // red text
        else
            color = "\033[1;37;41m"; // bright white text on red background

        const std::string reset_color = "\033[0m";
        std::cout << color << m_os.str() << reset_color << std::endl;
#endif // LINUX
#ifdef WIN_NT
        if (m_level == mi::base::MESSAGE_SEVERITY_WARNING)
            std::cout << util::bg_black << util::fg_yellow
                      << m_os.str()
                      << util::reset_color << std::endl;
        else
            std::cout << util::bg_red << util::fg_white
                      << m_os.str()
                      << util::reset_color << std::endl;
#endif // WIN_NT
    }
    else
    {
        std::cout << m_os.str() << std::endl;
    }

    if (m_level == mi::base::MESSAGE_SEVERITY_FATAL)
    {
#ifdef LINUX
        // Print a simple backtrace. To make it more useful, use c++filt to demangle the symbols and
        // addr2line to get line numbers.
        std::cout << "Backtrace:" << std::endl;
        void* array[255];
        size_t size = backtrace(array, 255);
        char** buffer = backtrace_symbols(array, size);
        if (buffer)
        {
            for (size_t i = 0; i < size; ++i)
            {
                std::cout << "  (" << (i + 1) << ") " << buffer[i] << std::endl;
            }
            free(buffer);
        }
#endif // LINUX

        // When a log message of severity 'fatal' is received, then DiCE will terminate the
        // application after printing the message. To find out what caused the error, we call
        // abort() here so that a proper stack trace is created when running with a debugger.

#ifdef DEBUG
        std::cout << "Receiving_logger got fatal error, aborting." << std::endl;
        abort(); // to get backtrace in a debugger
#else
#ifdef LINUX
        debug_break();
#endif // LINUX
#endif // DEBUG
    }
}

//----------------------------------------------------------------------
std::ostringstream& Receiving_logger::get_message(mi::Uint32 level)
{
    m_level = level;
    return m_os;
}

//----------------------------------------------------------------------
void Receiving_logger::message(
    mi::base::Message_severity level,
    const char*                category,
    const char*                message)
{
    Receiving_logger().get_message(level) << message;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
