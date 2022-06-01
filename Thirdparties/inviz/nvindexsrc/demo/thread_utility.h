/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief pthread helper. If possible, better to use boost or c++11 thread.

#ifndef NVIDIA_INDEX_THREAD_UTILITY_H
#define NVIDIA_INDEX_THREAD_UTILITY_H

#include <mi/dice.h>

#ifndef MI_PLATFORM_WINDOWS
#include <pthread.h>
#else
#include <windows.h>
#include <process.h>
#endif

// A simple thread class wrapping the Windows API or POSIX threads.
class Application_thread_base
{
public:
    // Constructor
    Application_thread_base();

    // Destructor
    virtual ~Application_thread_base();

    // Starts this thread.
    //
    // \return   \c true for success or \c false for failure.
    bool start();

    // Waits for this thread to finish.
    void join();

protected:
    // The actual function for doing work on this thread. This function must be
    // implemented by any subclass.
    virtual void run() = 0;

private:
    // Lock for #m_was_started.
    mi::base::Lock m_thread_lock;

    // Flag that indicates whether the thread was already started.
    bool m_was_started;

#ifndef MI_PLATFORM_WINDOWS
    // The thread ID.
    pthread_t m_tid;

    // Helper function for starting the thread.
    static void* do_run( void* thread);
#else
    // The thread ID
    unsigned int m_tid;

    // The thread handle.
    HANDLE m_handle;

    // Helper function for starting the thread.
    static unsigned __stdcall do_run( void* thread);
#endif
};

#endif // #ifndef NVIDIA_INDEX_THREAD_UTILITY_H
