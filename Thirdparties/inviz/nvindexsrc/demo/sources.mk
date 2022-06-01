#*****************************************************************************
# Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
#*****************************************************************************

SRC := \
        autotracking_workflow_operation.cpp             \
        bridge_video_stream.cpp                         \
        bridge_server.cpp                               \
        camera_animator.cpp                             \
        camera_utility.cpp                              \
        colormap_manager.cpp                            \
        colormap_util.cpp                               \
        command_handler_base.cpp                        \
        command_handler_colormap.cpp                    \
        command_handler_scene_element.cpp               \
        command_handler_settings.cpp                    \
        compute_utility.cpp                             \
        distributed_amplitude_map_job.cpp               \
        distributed_command_execution.cpp               \
        distributed_heightfield_elevation_delete.cpp    \
        distributed_heightfield_seeding_updates.cpp     \
        distributed_seed_point_evaluation.cpp           \
        distributed_volume_bordermap_generation.cpp     \
        distributed_volume_filter.cpp                   \
        examiner_manipulator.cpp                        \
        geostream_io.cpp                                \
        geostream_viewer.cpp                            \
        geostream_viewer_main.cpp                       \
        heightfield_data_export.cpp                     \
        heightfield_data_retrieval_appjob.cpp           \
        heightfield_seed_manipuration.cpp               \
        heightfield_workflow_functionality.cpp          \
        host_info.cpp                                   \
        html5_video_stream.cpp                          \
        http_handler.cpp                                \
        image_tile.cpp                                  \
        irregular_volume_data_editing.cpp               \
        jsoncpp.cpp                                     \
        multiple_camera.cpp                             \
        nvindex_appdata.cpp                             \
        nvindex_library_accessor.cpp                    \
        opengl_appdata.cpp                              \
        opengl_application_buffer.cpp                   \
        performance_log.cpp                             \
        performance_logger_app_state.cpp                \
        performance_logger_base.cpp                     \
        performance_logger_per_host_csv.cpp             \
        performance_logger_per_host_hr.cpp              \
        performance_logger_per_host_qa.cpp              \
        performance_logger_per_span_csv.cpp             \
        performance_logger_set.cpp                      \
        pick_tool.cpp                                   \
        polygonal_cylinder_gen.cpp                      \
        ring_buffer.cpp                                 \
        rtc_parameter_buffer_manip.cpp                  \
        rtmp_call_event_handler.cpp                     \
        rtmp_connect_event_handler.cpp                  \
        rtmp_handler.cpp                                \
        scene_element_importer.cpp                      \
        scene_element_importer_attributes.cpp           \
        scene_element_importer_massive_shapes.cpp       \
        scene_element_importer_raster_shapes.cpp        \
        scene_element_importer_simple_shapes.cpp        \
        scene_utility.cpp                               \
        shapes_utility.cpp                              \
        simple_gridding_operation.cpp                   \
        snapping_tool.cpp                               \
        span_renderer_no_gl.cpp                         \
        stereo_camera.cpp                               \
        thread_utility.cpp                              \
        user_interaction.cpp                            \
        utilities.cpp                                   \
        volume_bordermap_database_element.cpp           \
        volume_brick_element.cpp                        \
        volume_data_export.cpp                          \
        volume_data_filter.cpp                          \
        volume_data_retrieval_appjob.cpp                \
        volume_filter_compute_task.cpp                  \
        websocket_utility.cpp                           \


# OpenGL specific files
SRC_GL  :=                                      \
        glew.c                                  \
        opengl_application_draw.cpp             \
        opengl_drawing_utilities.cpp            \
        opengl_offscreen_context.cpp            \
        span_renderer_gl.cpp                    \

# OpenGL specific files only for Linux
SRC_GL_LINUX :=                                 \
        NVCtrl.c                                \
        transport.c                             \
        xwin_fastxopendisplay.cpp               \
        xwin_scoped_error_handler.cpp           \

SRC_GL_WINDOWS :=                               \

# Files for building without OpenGL support
SRC_NO_GL       :=                              \
        opengl_application_draw_no_gl.cpp       \
        opengl_drawing_utilities_no_gl.cpp      \

# Files for MPI IPC Compute
SRC_MPI_IPC_COMPUTE :=                          \
        cuda_ipc_volume_editing.cpp             \
        distributed_volume_application.cpp      \

SRC_MPI_IPC_COMPUTE_CUDA :=                     \
        cuda_ipc_volume_editing_kernel.cu       \


# Local Variables: ***
# mode: Makefile ***
# End: ***
