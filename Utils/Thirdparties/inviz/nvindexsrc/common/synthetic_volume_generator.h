/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief Volume data generator.

#ifndef NVIDIA_INDEX_BIN_COMMON_SYNTHETIC_VOLUME_DATA_GENERATOR_H
#define NVIDIA_INDEX_BIN_COMMON_SYNTHETIC_VOLUME_DATA_GENERATOR_H

#include <string>

#include <mi/dice.h>
#include <mi/base/interface_declare.h>
#include <mi/base/interface_implement.h>

#include <nv/index/idistributed_data_import_callback.h>

#include "encode_voxel.h"


namespace nv {
namespace index_common {


/// The synthetic volume data generation enables the display of large-scale volumes,
/// the benchmarking and testing of the large-scale, distributed volume rendering. 
/// The synthetic volume generation thus enables illustrating NVIDIA IndeX features
/// without the need of proprietary datasets and data importers.
/// The synthetic data generation allows to scale the amount of volume data to any
/// appropriate size.
class Synthetic_volume_generator :
    public nv::index::Distributed_discrete_data_import_callback<0x46cfd7d8,0x35d9,0x49e1,0x8d,0x8d,0xb2,0x5c,0x01,0xea,0xfa,0x83>
{
public:
    /// Available set of synthesis methods.
    enum Synthesis_method
    {
        /// default synthesis method
        DEFAULT = 0,
        /// ijk encoding for i and j
        IJ = 1,
        /// ijk encoding for j and k
        JK = 2,
        /// ijk encoding for k and i
        KI = 3,
        /// ijk encoding for i
        I = 4,
        /// ijk encoding for j
        J = 5,
        /// ijk encoding for k
        K = 6,
        /// sphere type 0
        SPHERE_0 = 7,
        /// for a transparent test
        TRANSPARENT_TEST = 8,
        /// Perlin noise (need parameters)
        PERLIN_NOISE = 9,
        /// fill the volume with zeroes
        ZERO = 10,
        /// sentinel
        COUNT
    };

    /// Voxel format
    enum Voxel_format
    {
        /// 8 bits unsigned (colormap index)
        VOXEL_FORMAT_UINT_8  = 0,
        /// 32 bits Packed RGBA (raw)
        VOXEL_FORMAT_RGBA_8  = 1,
        /// 8 bits unsigned (colormap index)
        VOXEL_FORMAT_UINT_16  = 2,
        /// 32 float normalized values
        VOXEL_FORMAT_FLOAT_32 = 3
    };
    
    /// Synthetic volume generator constructor
    ///
    /// \param[in] bounding_box      whole volume bounding box
    /// \param[in] voxel_format      voxel format type
    /// \param[in] synthesis_method  synthesis method name
    /// \param[in] parameters        synthesis method parameters in a string format
    ///
    Synthetic_volume_generator(
        const mi::math::Bbox_struct<mi::Uint32, 3>& bounding_box,
        Voxel_format                                voxel_format,
        const std::string&                          synthesis_method,
        const std::string&                          parameters);

    /// The default constructor is required for serialization only.
    Synthetic_volume_generator();
    
    // -------------------------------------------------------------------------
    // Implements the methods of the importer callback interface \c IDistributed_data_import_callback to
    // estimate memory consuptions and to create a volume brick.
    virtual mi::Size estimate(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bounding_box,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;
    
    virtual nv::index::IDistributed_data_subset* create(
        const mi::math::Bbox_struct<mi::Sint32, 3>&     bbox,
        nv::index::IData_subset_factory*                factory,
        mi::neuraylib::IDice_transaction*               dice_transaction) const;

    virtual mi::base::Uuid subset_id() const;

    // -------------------------------------------------------------------------------------------
    /// Returns configuration settings that may be used by the library for the session export
    /// mechanism provided by ISession::export_session().
    virtual const char* get_configuration() const;

    // -------------------------------------------------------------------------------------------
    /// Serialize the class to the given serializer.
    /// \param serializer Write to this serializer.
    virtual void serialize(mi::neuraylib::ISerializer* serializer) const;

    /// Deserialize the class from the given deserializer.
    /// \param deserializer Read from this deserializer.
    virtual void deserialize(mi::neuraylib::IDeserializer* deserializer);
    
private:
    // -------------------------------------------------------------------------
    /// Generate synthetic voxel data based on the Perlin noise function.
    ///
    /// The parameter names are based on the paper, 
    ///
    ///    Ken Perlin, An Image Synthesizer, SIGGRAPH Volume 19, No. 3, 1985
    ///
    /// Here is a simple explanation only, for details, please refer the
    /// original paper.
    ///
    /// The parameter format is comma separated "key=value", e.g., 
    /// "cube_unit=512,time=0.0".
    ///
    /// Accepted parameter keys
    /// - cube_unit: integer. The cube size of the noise
    ///   generator. This defines repeating period. If you have 512,
    ///   then you will see the discontinuity of the noise at the 512
    ///   position.
    /// - time: float.  Time variant parameter of Perlin noise.
    ///   If it does not depend on time, set to 0.0f.
    /// - terms: integer. How many terms will be used. Recommended 1 to 5,
    ///   the large number needs a heavier computation.
    /// - turbulence_weight: Vector<Float32, 4>.  turbulence flow
    ///   weight. Useful for Marble solid texture. try "turbulence_weight=0 0 0 0"
    ///   "turbulence_weight=1 1 1 1"
    /// - abs_noise: boolean. Take absolute value of the noise. try 0 first.
    /// - ridged: boolean.  Add ridged effect. try 1 when you want to see
    ///   some discontinuity in your noise.
    ///
    /// Example.
    ///   "cube_unit=512,time=0,terms=2,turbulence_weight=1 1 1 1,abs_noise=0,ridged=0"
    ///
    /// \param[in] brick_data               The voxel data contained in a brick.
    /// \param[in] bounding_box             The brick's bounding box
    ///
    template<typename T>
    void generate_perlin_noise_volume(
        T*                                      brick_data,
        const mi::math::Bbox< mi::Sint64, 3 >&  raw_bbox,
        mi::neuraylib::IDice_transaction*       dice_transaction) const
    {    
        // get the parameters for Perlin noise
        String_dict perlin_param_dict;
        std::string err_mes;
        bool is_success = get_string_dict_from_inline_parameter_string(
            m_parameter_str, perlin_param_dict, err_mes);

        mi::Sint64 cube_unit = 512;
        mi::Float32 time = 0.0f;
        mi::Uint32 terms = 5;
        mi::math::Vector< mi::Float32, 4 > turbulence_weight(1.0f, 1.0f, 1.0f, 1.0f);
        bool abs_noise = false;
        bool ridged = false;

        if (is_success)
        {
            cube_unit = get_sint64(perlin_param_dict.get("cube_unit", "512"));
            time  = get_float32(perlin_param_dict.get("time", "0.0"));
            terms = get_uint32(perlin_param_dict.get("terms", "5"));
            turbulence_weight = get_vec_float32_4(perlin_param_dict.get("turbulence_weight", "1 1 1 1"));
            abs_noise = get_bool(perlin_param_dict.get("abs_noise", "0"));
            ridged    = get_bool(perlin_param_dict.get("ridged", "0"));
        }
        else
        {
            ERROR_LOG << "Failed to get Perlin noise parameters from [" << m_parameter_str << "] Use default.";
        }

        INFO_LOG << "Creating synthetic volume based on the Perlin noise function and using the following parameters: "
                  << "\n  cube unit [" << cube_unit
                  << "], time [" << time
                  << "], terms [" << terms 
                  << "], turbulence weight [" << turbulence_weight
                  << "], absolute noise [" << abs_noise
                  << "], ridged [" << ridged << "]";

        INFO_LOG << "Creating synthetic volume for " << raw_bbox;

        Parallel_perlin_noise<T> generate_perlin_noise(brick_data, raw_bbox, cube_unit, time, terms, turbulence_weight, abs_noise, ridged);
        dice_transaction->execute_fragmented(&generate_perlin_noise, (raw_bbox.max[0]-raw_bbox.min[0]));
    }


    /// Generate a synthetic volume based on the given types.
    ///
    /// \param[in] brick_data               The voxel data contained in a brick.
    /// \param[in] bounding_box             The brick's bounding box
    ///
    void generate_synthetic_data(
        mi::Uint8*                              brick_data,
        const mi::math::Bbox<mi::Sint64, 3>&    bounding_box) const;
        
    void generate_synthetic_data(
        mi::math::Vector_struct<mi::Uint8, 4>*  brick_data,
        const mi::math::Bbox<mi::Sint64, 3>&    bounding_box) const;
        
    void generate_synthetic_data(
        mi::Float32*                            brick_data,
        const mi::math::Bbox<mi::Sint64, 3>&    bounding_box,
        const mi::math::Vector<mi::Float32, 2>& range,
        bool                                    normalize) const;

private:
    // Whole volume bbox in datasets local space.
    mi::math::Bbox<mi::Uint32, 3>           m_whole_bbox;
    // Associated session.
    mi::neuraylib::Tag                      m_session_tag;
    // Synthetic volume generation technique.
    mi::Uint32                              m_synthetic_type;
    // Volume voxel format.
    Voxel_format                            m_voxel_format;
    // Additional importer parameter.
    std::string                             m_parameter_str;
    // configuration information for the session exporter
    std::string                             m_configuration;
};

}} // namespace nv::index_common

#endif // NVIDIA_INDEX_BIN_COMMON_SYNTHETIC_VOLUME_DATA_GENERATOR_H
