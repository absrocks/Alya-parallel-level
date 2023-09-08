/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "ring_buffer.h"

#include <sstream>
#include <iostream>

//----------------------------------------------------------------------
Ring_buffer::Ring_buffer()
    :
    m_data_buf(),
    m_buf_size(0),
    m_begin_idx(0),
    m_end_idx(0),
    m_data_size(0)
{
    this->resize_buffer(1024);
}

//----------------------------------------------------------------------
Ring_buffer::Ring_buffer(size_t buf_size)
    :
    m_data_buf(),
    m_buf_size(0),
    m_begin_idx(0),
    m_end_idx(0),
    m_data_size(0)
{
    this->resize_buffer(buf_size);
}

//----------------------------------------------------------------------
Ring_buffer::~Ring_buffer()
{
    m_data_buf.clear();
    m_buf_size  = 0;
    m_begin_idx = 0;
    m_end_idx   = 0;
    m_data_size = 0;
}

//----------------------------------------------------------------------
void Ring_buffer::resize_buffer(size_t buf_size)
{
    assert(buf_size > 0);

    m_data_buf.resize(buf_size);
    m_buf_size  = buf_size;
    m_begin_idx = 0;
    m_end_idx   = 0;
    m_data_size = 0;

    assert(m_buf_size == m_data_buf.size());
}

//----------------------------------------------------------------------
void Ring_buffer::push_back(value_type const & dat)
{
    if(this->full()){
        // remove the begin and push back at the end
        // std::cout << "full" << std::endl;
        this->pop_front();
    }
    if(this->empty()){
        // if empty(), we just insert the data at the current end
        // position. No end move, but the size is increase.
        ++m_data_size;
        assert(m_end_idx < m_data_buf.size());
        m_data_buf[m_end_idx] = dat;
    }
    else{
        this->inc_end();
        assert(m_end_idx < m_data_buf.size());
        m_data_buf[m_end_idx] = dat;
    }
}

//----------------------------------------------------------------------
void Ring_buffer::pop_front()
{
    this->inc_begin();
}

//----------------------------------------------------------------------
Ring_buffer::value_type const & Ring_buffer::front() const
{
    assert(!this->empty());
    assert(m_begin_idx < m_data_buf.size());

    return m_data_buf[m_begin_idx];
}

//----------------------------------------------------------------------
Ring_buffer::value_type const & Ring_buffer::back() const
{
    assert(!this->empty());
    assert(m_end_idx < m_data_buf.size());

    return m_data_buf[m_end_idx];
}

//----------------------------------------------------------------------
size_t Ring_buffer::size() const
{
    return m_data_size;
}

//----------------------------------------------------------------------
size_t Ring_buffer::capacity() const
{
    return m_buf_size;
}

//----------------------------------------------------------------------
bool Ring_buffer::empty() const
{
    if(m_data_size == 0){
        assert(m_begin_idx == m_end_idx);
        return true;
    }
    assert(m_data_size > 0);
    return false;
}

//----------------------------------------------------------------------
bool Ring_buffer::full() const
{
    if(m_data_size == m_buf_size){
        return true;
    }
    return false;
}

//----------------------------------------------------------------------
std::string Ring_buffer::to_string() const
{
    assert(m_data_buf.size() == m_buf_size);

    std::stringstream sstr;
    sstr << "buffsize: " << m_data_buf.size()
         << ", internal index [" << m_begin_idx << "," << m_end_idx
         << "], datasize: " << m_data_size;
    if(this->empty()){
        sstr << ", empty.";
    }
    if(this->full()){
        sstr << ", full.";
    }

    return sstr.str();
}

//----------------------------------------------------------------------
Ring_buffer::iterator Ring_buffer::begin()
{
    return Ring_buffer::iterator(this, true);
}

//----------------------------------------------------------------------
Ring_buffer::iterator Ring_buffer::end()
{
    return Ring_buffer::iterator(this, false);
}

//----------------------------------------------------------------------
void Ring_buffer::inc_begin()
{
    assert(!this->empty());
    assert(m_data_size > 0);

    ++m_begin_idx;
    --m_data_size;

    // if index is over the buffer, back to 0
    if(m_begin_idx == m_buf_size){
        m_begin_idx = 0;
    }

    assert(m_begin_idx < m_buf_size);
}

//----------------------------------------------------------------------
void Ring_buffer::inc_end()
{
    assert(!this->full());

    ++m_end_idx;
    ++m_data_size;
    assert(m_data_size <= m_buf_size);

    // if index is over the buffer, back to 0
    if(m_end_idx == m_buf_size){
        m_end_idx = 0;
    }

    assert(m_end_idx < m_buf_size);
}

//----------------------------------------------------------------------
