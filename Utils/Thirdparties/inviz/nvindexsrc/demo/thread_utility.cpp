/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief pthread helper. If possible, better to use boost or c++11 thread.

#include "thread_utility.h"

#ifndef MI_PLATFORM_WINDOWS

Application_thread_base::Application_thread_base()
    :
    m_was_started(false)
{
    // empty
}

Application_thread_base::~Application_thread_base()
{
    // empty
}

bool Application_thread_base::start()
{
    mi::base::Lock::Block block( &m_thread_lock);

    if (m_was_started)
    {
        return false;
    }
    m_was_started = true;

    pthread_create( &m_tid, 0, do_run, this);

    return true;
}

void Application_thread_base::join()
{
    mi::base::Lock::Block block( &m_thread_lock);
    if (m_was_started)
    {
        pthread_join( m_tid, 0);
    }
}

void* Application_thread_base::do_run( void* thread)
{
    ((Application_thread_base *) thread)->run();
    return 0;
}

#else // MI_PLATFORM_WINDOWS

Application_thread_base::Application_thread_base() : m_was_started( false)
{
    m_handle = (HANDLE) _beginthreadex( 0, 0, &Application_thread_base::do_run, this, CREATE_SUSPENDED, &m_tid);
}

Application_thread_base::~Application_thread_base()
{
    CloseHandle( m_handle);
}

bool Application_thread_base::start()
{
    mi::base::Lock::Block block( &m_thread_lock);
    if (m_was_started)
    {
        return false;
    }
    m_was_started = true;

    DWORD result = ResumeThread( m_handle);
    return result != -1;
}

void Application_thread_base::join()
{
    mi::base::Lock::Block block( &m_thread_lock);
    if (m_was_started)
    {
        WaitForSingleObject( m_handle, INFINITE);
    }
}

unsigned __stdcall Application_thread_base::do_run( void *thread)
{
    ((Application_thread_base *) thread)->run();
    _endthreadex( 0);
    return 0;
}

#endif // MI_PLATFORM_WINDOWS
