/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief distributed seed point evaluation computation

#ifndef NVIDIA_INDEX_DISTRIBUTED_SEED_POINT_EVALUATION_H
#define NVIDIA_INDEX_DISTRIBUTED_SEED_POINT_EVALUATION_H

#include <mi/base/types.h>
#include <mi/dice.h>
#include <mi/base/interface_declare.h>

#include <nv/index/idistributed_compute_algorithm.h>
#include <nv/index/idistributed_data_locality.h>

#include <nv/index/iregular_volume.h>
#include <nv/index/idistributed_data_access.h>
#include <nv/index/idistributed_data_locality.h>

#include <vector>


///
class Flow_grid : public mi::neuraylib::Element<0x9ecd5093,0xe5f1,0x43bf,0x8e,0x5c,0xf0,0x13,0x69,0xf2,0x74,0x1e>
{
public:
    ///
    Flow_grid(const mi::math::Bbox_struct<mi::Uint32, 2>& range);
    Flow_grid()
        :
        m_range(),
        m_x_range(0),
        m_y_range(0),
        m_flow_grid(),
        m_max_value(0.0f),
        m_min_value(0.0f)
    {
        // for serialization only
    }

    virtual ~Flow_grid();

    /// ACCESSING
    const mi::math::Bbox_struct<mi::Uint32, 2>& get_range() const { return m_range; }
    const std::vector<mi::Float32>& get_flow_grid() const { return m_flow_grid; }
    mi::Float32 get_sample_value(
        mi::Uint32 i_value,
        mi::Uint32 j_value) const;
    mi::Float32 get_sample(
        mi::Uint32 i_value,
        mi::Uint32 j_value) const;

    void get_range_with_border(mi::Uint32& flow_x_range, mi::Uint32& flow_y_range) const
    {
        flow_x_range = m_x_range;
        flow_y_range = m_y_range;
    }

    void get_min_max_heights(
        mi::math::Vector_struct<mi::Float32, 2>& min_max_heights) const
    {
        min_max_heights.x = m_min_value;
        min_max_heights.y = m_max_value;
    }

    /// EDITING
    void set_sample_value(mi::Uint32 i_value, mi::Uint32 j_value, mi::Float32 k_value) const;
    void set_border_sample(const mi::math::Vector_struct<mi::Float32, 4>& border) const;

    /// -------------------------------------------------------------------------------------------
    /// Serialize the class to the given serializer.
    /// \param serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;

    /// Deserialize the class from the given deserializer.
    /// \param serializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

    /// The copy needs to create and return a newly allocated object of the same
    /// class and has to copy all the fields which should be the same in the copy.
    /// \return The new copy of the heightfield database element
    virtual mi::neuraylib::IElement* copy() const;

    /// Give back a human readable representation of the class name
    virtual const char* get_class_name() const  { return "Flow_grid"; }

    ///
    virtual void get_references(mi::neuraylib::ITag_set* result) const {  }

private:
    mi::math::Bbox_struct<mi::Uint32, 2>    m_range;
    mi::Uint32                              m_x_range;
    mi::Uint32                              m_y_range;

    mutable std::vector<mi::Float32>        m_flow_grid;
    mutable mi::Float32                     m_max_value;
    mutable mi::Float32                     m_min_value;
};

class ICompute_height_value
{
public:
    virtual void compute_height_value(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        const mi::math::Vector_struct<mi::Float32, 4>&  input_seed,
        mi::math::Vector_struct<mi::Float32, 4>&        output_seed) = 0;
};

/// Proof of concept computing algorithm .. .
class Distributed_seed_point_evaluation :
    public nv::index::Distributed_compute_algorithm<0xfef5cfb3,0xbb7a, 0x420c,0xad,0x39,0xdf,0xc1,0x0f,0x4e,0xe4,0x2d>
{
public:
    Distributed_seed_point_evaluation(
        const mi::neuraylib::Tag&                                       session_tag,
        const mi::neuraylib::Tag&                                       heightfield_tag,
        const mi::neuraylib::Tag&                                       volume_tag,
        const std::vector<mi::math::Vector_struct<mi::Float32, 3> >&    seed_points,
        const std::vector<mi::Uint32>&                                  cluster_hosts,
        const mi::math::Vector_struct<mi::Uint32, 2>&                   heightfield_size,
        bool                                                            eight_way_seedings);
    Distributed_seed_point_evaluation()
        :
        m_session_tag(mi::neuraylib::NULL_TAG),
        m_heightfield_tag(mi::neuraylib::NULL_TAG),
        m_volume_tag(mi::neuraylib::NULL_TAG),
        m_seed_points(),
        m_cluster_hosts(),
        m_update_heightfield(0),
        m_dangling_points(),
        m_border_points(),
        m_samples(),
        m_dangling_samples(),
        m_border_samples(),
        m_flow_grid_tags(),
        m_volume_brick_tags(),
        m_patch_tag_vectors(),
        m_heightfield_size(),
        m_eight_way_seedings(0)
    {
        // for serialization only
    }
    virtual ~Distributed_seed_point_evaluation();

public:
    virtual mi::Size get_nb_of_fragments() const { return m_cluster_hosts.size(); }

    virtual void get_updated_bounding_box(mi::math::Bbox_struct<mi::Float32, 3>& bbox) const
    {
        // empty
    }
    
    const std::vector<std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> > >& get_patch_tag_lists() const
    {
        return m_patch_tag_vectors;
    }

public:
    ///
    void update_seed_points();
    mi::Uint32 get_nb_seed_points() const;

    ///
    void clean_up(
        mi::neuraylib::IDice_transaction*   dice_transaction);

    /// Enable heightfield elevation value updates.
    void enable_elevation_updates() { m_update_heightfield = 1; }
    void disable_elevation_updates() { m_update_heightfield = 0; }

public:
    // Implemented fragmented jobs method

    virtual mi::neuraylib::IFragmented_job::Scheduling_mode get_scheduling_mode() const
    {
        return mi::neuraylib::IFragmented_job::USER_DEFINED;
    }

    virtual void assign_fragments_to_hosts(
        mi::Uint32* slots,
        mi::Size    nr_slots);

    virtual void execute_fragment(
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

    virtual void execute_fragment_remote(
		mi::neuraylib::ISerializer*                     serializer,
        mi::neuraylib::IDice_transaction*               dice_transaction,
        mi::Size                                        index,
        mi::Size                                        count,
        const mi::neuraylib::IJob_execution_context*    context);

    virtual void receive_remote_result(
        mi::neuraylib::IDeserializer*       deserializer,
        mi::neuraylib::IDice_transaction*   dice_transaction,
        mi::Size                            index,
        mi::Size                            count);

    virtual void serialize(mi::neuraylib::ISerializer *serializer) const;

    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);

private:
    void evaluate_seed_flow(
        mi::neuraylib::IDice_transaction*                       dice_transaction,
        nv::index::IRegular_heightfield_data_locality*          heightfield_locality,
        nv::index::IRegular_volume_data_locality*               volume_locality,
        mi::Uint32                                              host_id,
        std::vector<mi::neuraylib::Tag>&                        flow_grid_tags,
        std::vector<std::vector<mi::neuraylib::Tag> >&          volume_brick_tags,
        std::vector<mi::math::Vector_struct<mi::Float32, 4> >&  result_samples,
        std::vector<mi::math::Vector_struct<mi::Float32, 4> >&  dangling_samples,
        std::vector<mi::math::Vector_struct<mi::Float32, 4> >&  border_samples);

    void evaluate_seed_flow(
        mi::neuraylib::IDice_transaction*                       dice_transaction,
        const mi::math::Bbox_struct<mi::Uint32, 2>&             patch,
        const Flow_grid*                                        flow_grid,
        ICompute_height_value*                                  compute_height_value,
        std::vector<mi::math::Vector_struct<mi::Float32, 4> >&  result_samples,
        std::vector<mi::math::Vector_struct<mi::Float32, 4> >&  dangling_samples) const;

private:
    mi::neuraylib::Tag                                      m_session_tag;
    mi::neuraylib::Tag                                      m_heightfield_tag;
    mi::neuraylib::Tag                                      m_volume_tag;

    std::vector<mi::math::Vector_struct<mi::Float32, 4> >   m_seed_points;

    // First implemented the user-defined scheduling.
    std::vector<mi::Uint32>                                 m_cluster_hosts;
    mi::Uint32                                              m_update_heightfield;

    std::vector<mi::math::Vector_struct<mi::Float32, 4> >   m_dangling_points;
    std::vector<mi::math::Vector_struct<mi::Float32, 4> >   m_border_points;

    std::vector<std::vector<mi::math::Vector_struct<mi::Float32, 4> > > m_samples;
    std::vector<std::vector<mi::math::Vector_struct<mi::Float32, 4> > > m_dangling_samples;
    std::vector<std::vector<mi::math::Vector_struct<mi::Float32, 4> > > m_border_samples;
    std::vector<std::vector<mi::neuraylib::Tag> >                       m_flow_grid_tags;

    std::vector<std::vector<std::vector<mi::neuraylib::Tag> > >         m_volume_brick_tags;
    std::vector<std::vector<std::pair<mi::math::Bbox_struct<mi::Uint32, 2>, mi::neuraylib::Tag> > > m_patch_tag_vectors;

    mi::math::Vector_struct<mi::Uint32, 2>                              m_heightfield_size;
    mi::Uint32                                                          m_eight_way_seedings;
};


#endif // NVIDIA_INDEX_DISTRIBUTED_SEED_POINT_EVALUATION_H
