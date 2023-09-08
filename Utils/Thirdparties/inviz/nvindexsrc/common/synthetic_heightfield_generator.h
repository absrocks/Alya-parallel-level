/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external heightfield data generator example implementation

#ifndef NVIDIA_INDEX_SYNTHETIC_HEIGHTFIELD_GENERATOR_H
#define NVIDIA_INDEX_SYNTHETIC_HEIGHTFIELD_GENERATOR_H

#include <string>

#include <mi/dice.h>

#include <nv/index/idistributed_data_import_callback.h>

namespace nv {
namespace index_common {

/// Configuration of synthetic heightfield data generation declaration.
///
/// The NVIDIA IndeX reference viewer uses the following configuration
/// in its project file to generate a synthetic heightfield.
///
/// <pre>
/// #
/// # synthetic heightfield (type: ij)
/// #
/// app::heightfield::synthetic_ij::name = synthetic_ij
/// # heightfield whole size
/// app::heightfield::synthetic_ij::size = 200 100
/// app::heightfield::synthetic_ij::clip = 0 0 0 200 200 200
/// app::heightfield::synthetic_ij::range = 0.0 250.0
/// app::heightfield::synthetic_ij::importer = synthetic
/// app::heightfield::synthetic_ij::color = 0.4 0.4 0.7 1.0
/// # synthetic_type = {ij, i , j, 0}
/// app::heightfield::synthetic_ij::synthetic_type = ij
/// </pre>
///
/// The synthetic heightfield patch generation allows the following techniques:
/// - Constant height surface.
/// - Sawtooth like surface.
/// - Perlin noise based turbulence (default).
///
/// The following generation technique generates the distributed and independent patch data.
/// 
class Synthetic_heightfield_generator :
        public nv::index::Distributed_discrete_data_import_callback<0x938d6269,0xb929,0x43a6,0xab,0x79,0x9a,0xbd,0x68,0x6e,0xc8,0x3f>
{
public:
    /// The default constructor is merely used for serialization/deserialization.
    Synthetic_heightfield_generator();

    /// The constructor requires the following parameter to generate the patches of an 
    /// synthetic heightfield:
    ///
    /// \param[in] technique            The identifier used to choose the generation technique.
    ///
    /// \param[in] technique_param      The parameters of generation technique. (e.g., noise parameter)
    ///
    /// \param[in] heightfield_size     The size of the heightfield.
    ///
    Synthetic_heightfield_generator(
        const std::string&                              technique,
        const std::string&                              technique_param,
        const mi::math::Vector_struct<mi::Uint32, 2>&   heightfield_size);

    /// The destructor will be called by the NVIDIA IndeX library.
    virtual ~Synthetic_heightfield_generator();

    // -------------------------------------------------------------------------
    // Implements the methods of the importer callback interface \c IDistributed_data_import_callback to
    // estimate memory consumption and to create a volume brick.
    virtual mi::Size estimate(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;
    
    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    virtual mi::base::Uuid subset_id() const;

    /// -------------------------------------------------------------------------------------------
    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    /// -------------------------------------------------------------------------------------------
    /// Serialize the class to the given serializer.
    virtual void serialize(mi::neuraylib::ISerializer * serializer) const;

    /// Deserialize the class from the given deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer * deserializer);

private:
    /// Generate height with a function without parameter
    void generate_height_no_parameter(const mi::math::Bbox<mi::Sint64, 2> & patch_raw_bbox,
                                      mi::Float32 * const p_heightfield_data,
                                      const mi::Sint32 syn_type) const;

    /// Generate height with a function of HSM_0 (with the height parameter)
    void generate_height_0(const mi::math::Bbox<mi::Sint64, 2> & patch_raw_bbox,
                           mi::Float32 * const p_heightfield_data) const;

    /// Generate height with a function with parameter
    void generate_height_perlin_noise(const mi::math::Bbox<mi::Sint64, 2> & patch_raw_bbox,
                                      mi::Float32 * const p_heightfield_data) const;

    // Name of the generation technique.
    std::string                             m_technique;
    // Parameter that steer the generation technique (e.g., noise parameter)
    std::string                             m_technique_param;
    // Size of the height field
    mi::math::Vector<mi::Uint32, 2>         m_size;

    // Configuration information for the session export
    std::string                             m_configuration;
};

//----------------------------------------------------------------------
}} // namespace nv::index_common

#endif // NVIDIA_INDEX_SYNTHETIC_HEIGHTFIELD_GENERATOR_H
