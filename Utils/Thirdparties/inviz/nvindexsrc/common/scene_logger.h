/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#ifndef NVIDIA_INDEX_BIN_COMMON_SCENE_LOGGER_H
#define NVIDIA_INDEX_BIN_COMMON_SCENE_LOGGER_H

#include <stack>
#include <vector>

#include <nv/index/iscene_visitor.h>

namespace nv {
namespace index_common {

/// Implemented a simple visitor that logs the scene contents.
/// Should better be used outside of the library to illustrate how
/// to traverse the scene graph and visualize dependencies.
class Scene_logger : public nv::index::IScene_visitor
{
public:
    Scene_logger(bool forward_evaluation_only = true);
    virtual  ~Scene_logger() {}

    /// Evaluate each and every scene graph node (tag) separately.
    virtual void evaluate(
        mi::neuraylib::Tag                                  scene_element,
        const mi::base::Uuid&                               uuid,
        nv::index::IScene_visitor::Scene_evaluation_mode    evaluation_mode,
        mi::neuraylib::IDice_transaction*                   dice_transaction);

    const std::vector<mi::neuraylib::Tag>& get_color_maps() const;

private:
    bool                                m_forward_evaluation_only;
    mi::Uint32                          m_indentation;

    std::stack<mi::Uint32>              m_nb_scene_nodes;

    std::vector<mi::neuraylib::Tag>     m_color_maps;
};

}} // namespace index_common / nv

#endif // NVIDIA_INDEX_BIN_COMMON_SCENE_LOGGER_H
