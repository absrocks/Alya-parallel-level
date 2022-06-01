#*****************************************************************************
# Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
#*****************************************************************************

SRC := \
        affinity_information.cpp                        \
        clock_pulse_generator.cpp                       \
        color_space.cpp                                 \
        colormap_io.cpp                                 \
        colormap_ppm.cpp                                \
        distributed_compute_techniques.cpp              \
        distributed_heightfield_elevation_change.cpp    \
        distributed_voxel_value_change.cpp              \
        forwarding_logger.cpp                           \
        irregular_volume_importer.cpp                   \
        irregular_volume_importer_vtk.cpp               \
        large_file_io.cpp                               \
        multi_attribute_scaling_sequence_importer.cpp   \
        normal_encoder.cpp                              \
        pipe.cpp                                        \
        ppm_heightfield_importer.cpp                    \
        ppm_io.cpp                                      \
        raw_heightfield_data_importer.cpp               \
        raw_heightfield_data_importer_ssv.cpp           \
        raw_volume_data_importer.cpp                    \
        raw_volume_data_sequence_importer.cpp           \
        receiving_logger.cpp                            \
        repeat_raw_volume_data_importer.cpp             \
        reservoir_grid_importer.cpp                     \
        rgba_raw_volume_data_importer.cpp               \
		scene_logger.cpp                                \
        scenegraph_utility.cpp                          \
        short_value_raw_volume_data_importer.cpp        \
        single_value_raw_volume_data_importer.cpp       \
        spline.cpp                                      \
        string_dict.cpp                                 \
        synthetic_heightfield_generator.cpp             \
        synthetic_reservoir_patch_generator.cpp         \
        synthetic_volume_generator.cpp                  \
        transferfunction_colormap.cpp                   \
        triangle_mesh_importer.cpp                      \
        volume_replica_generation.cpp                   \
        windows_utilities.cpp                           \

SRC_OPENVDB :=                                       \
        sparse_volume_importer.cpp                   \

# Local Variables: ***
# mode: Makefile ***
# End: ***
