/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief a simple ring buffer

#ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RING_BUFFER_H
#define NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RING_BUFFER_H

#include <mi/dice.h>

#include <vector>
#include <string>
#include <cassert>
#include <iostream>

/// a fixed size ring buffer
class Ring_buffer {
public:
    /// data type
    typedef mi::Float64 value_type;

public:
    /// default constructor. default buffer size = 1024.
    Ring_buffer();
    /// constructor
    /// \param[in] buf_size initial buffer size
    explicit Ring_buffer(size_t buf_size);
    /// destructor
    ~Ring_buffer();

    /// set the buffer size. The data is lost.
    /// \param[in] buf_size buffer size
    void resize_buffer(size_t buf_size);

    /// push back a data at the end
    void push_back(value_type const & dat);

    /// pop front the data
    void pop_front();

    /// get front. access to the 0th data.
    value_type const & front() const;

    /// get back. access to the n-1th data.
    value_type const & back() const;

    /// get valid data size
    size_t size() const;

    /// get the capacity (buffer size)
    size_t capacity() const;

    /// is empty the data?
    bool empty() const;

    /// is full the buffer?
    bool full() const;

    /// to string
    std::string to_string() const;

public:
    /// iterator
    class iterator
    {
    public:
        friend class Ring_buffer;

        /// constructor
        /// \param[in] p_rb     a ring buffer to iterate the elements
        /// \param[in] is_start iterator start when true, otherwise end.
        iterator(Ring_buffer * p_rb, bool is_start)
            :
            m_p_ringbuffer_ref(p_rb),
            m_i_begin_idx(p_rb->m_begin_idx),
            m_i_end_idx(  p_rb->m_end_idx),
            m_i_data_size(p_rb->m_data_size),
            m_i_idx(      p_rb->m_begin_idx)
        {
            // empty or !is_start, index is invalid.
            if((m_p_ringbuffer_ref->empty()) || (!is_start)){
                m_i_idx = ~0;   // invalid
            }
        }
        /// destructor
        ~iterator(){
            m_p_ringbuffer_ref = 0;
            m_i_idx = ~0;   // invalid
        }
        /// advance by one
        iterator& operator++(){
            if(m_i_idx == size_t(~0)){
                return (*this); // once invalid, stay it
            }

            if(m_i_idx == m_i_end_idx){
                m_i_idx = ~0;   // invalid
                return (*this);
            }

            ++m_i_idx;
            // std::cout << "m_i_idx: " << m_i_idx << "/" << m_i_data_size << " "
            //           << m_i_data_size << ", end: " << m_i_end_idx << std::endl;
            if(m_i_idx >= m_i_data_size){
                m_i_idx = 0;    // wrap around
            }
            return (*this);
        }

        /// equal?
        bool operator==(const iterator & ii) const {
            assert(m_p_ringbuffer_ref == ii.m_p_ringbuffer_ref);
            assert(m_i_begin_idx == ii.m_i_begin_idx);
            assert(m_i_end_idx   == ii.m_i_end_idx);
            assert(m_i_data_size == ii.m_i_data_size);

            return m_i_idx == ii.m_i_idx;
        }
        /// not equal?
        bool operator!=(const iterator& ii) const {
            assert(m_p_ringbuffer_ref == ii.m_p_ringbuffer_ref);
            assert(m_i_begin_idx == ii.m_i_begin_idx);
            assert(m_i_end_idx   == ii.m_i_end_idx);
            assert(m_i_data_size == ii.m_i_data_size);

            return m_i_idx != ii.m_i_idx;
        }

        /// dereference
        Ring_buffer::value_type operator*() const {
            // std::cout << "m_i_idx: " << m_i_idx << std::endl;
            assert(m_i_idx < m_p_ringbuffer_ref->m_data_buf.size());

            return m_p_ringbuffer_ref->m_data_buf[m_i_idx];
        }

    private:
        /// default constructor, prohibit until proben useful
        iterator();

    private:

        /// refeence to the ring buffer
        Ring_buffer * m_p_ringbuffer_ref;
        /// begin index const copy
        size_t const m_i_begin_idx;
        /// end index const copy
        size_t const m_i_end_idx;
        /// valid data size const copy
        size_t const m_i_data_size;
        /// current iterator index
        size_t m_i_idx;
    };

public:
    /// iterator begin
    iterator begin();

    /// iterator end
    iterator end();

private:
    /// increment the begin when not empty
    void inc_begin();
    /// increment the end when not full
    void inc_end();

private:
    /// data buffer
    std::vector< value_type > m_data_buf;
    /// buffer size
    size_t m_buf_size;
    /// begin index
    size_t m_begin_idx;
    /// end index
    size_t m_end_idx;
    /// valid data size
    size_t m_data_size;

private:
    /// copy constructor, disabled unless proven useful
    Ring_buffer(const Ring_buffer& _rhs);
    /// operator=, disabled unless proven useful
    const Ring_buffer& operator=(const Ring_buffer& _rhs);
};

#endif // #ifndef NVIDIA_INDEX_BIN_GEOSPATIAL_STREAM_VIEWER_RING_BUFFER_H
