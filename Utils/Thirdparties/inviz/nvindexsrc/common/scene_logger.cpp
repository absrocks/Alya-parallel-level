/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include <string>
#include <assert.h>

#include "scene_logger.h"

#include <nv/index/icolormap.h>
#include <nv/index/iscene.h>
#include <nv/index/iscene_element.h>
#include <nv/index/iscene_group.h>

#include "forwarding_logger.h"

namespace nv {

namespace index_common {

Scene_logger::Scene_logger(bool forward_evaluation_only) : m_forward_evaluation_only(forward_evaluation_only)
{
    m_indentation = 0;
}

void Scene_logger::evaluate(
    mi::neuraylib::Tag                                  element,
    const mi::base::Uuid&                               uuid,
    nv::index::IScene_visitor::Scene_evaluation_mode    evaluation_mode,
    mi::neuraylib::IDice_transaction*                   dice_transaction)
{
    const bool is_scene_group = (uuid == nv::index::IStatic_scene_group::IID()
        || uuid == nv::index::ITransformed_scene_group::IID()
        || uuid == nv::index::IShape_scene_group::IID()
        || uuid == nv::index::IScene::IID());

    // Logging the scene element preceeded by the indentation.
    mi::base::Handle<const nv::index::IScene_element> scene_element(
        dice_transaction->access<const nv::index::IScene_element>(element));
    assert(scene_element);

    if (evaluation_mode == IScene_visitor::DEPTH_FIRST_FORWARD_EVALUATION)
    {
        m_nb_scene_nodes.top()++;
        const std::string indents(m_indentation + 1, '\t');
        INFO_LOG << indents << "Level " << m_nb_scene_nodes.size() << ", element " << m_nb_scene_nodes.top() << ": " << scene_element->get_class_name() << " (tag id: " << element.id << "), (forward evaluation).";

        if (is_scene_group)
        {
            m_indentation++;
            m_nb_scene_nodes.push(0);
        }
        if (uuid==nv::index::IColormap::IID())
        {
            m_color_maps.push_back(element);
        }
    }
    else if (evaluation_mode == IScene_visitor::DEPTH_FIRST_BACKWARD_EVALUATION)
    {
        if (is_scene_group)
        {
            m_indentation--;
            m_nb_scene_nodes.pop();
        }

        const std::string indents(m_indentation + 1, '\t');
        if (!m_forward_evaluation_only)
        {
            INFO_LOG << indents << "Level " << m_nb_scene_nodes.size() << ", element " << m_nb_scene_nodes.top() << ": " << scene_element->get_class_name() << " (tag id: " << element.id << "), (backward evaluation).";
        }

        if (is_scene_group)
        {
            m_nb_scene_nodes.top()--;
        }
    }
    else if (evaluation_mode == IScene_visitor::START_EVALUATION)
    {
        const std::string indents(m_indentation + 1, '\t');
        INFO_LOG << indents << "Base Level: " << scene_element->get_class_name() << " (tag id: " << element.id << "), (started scene evaluation/traversal).";

        m_indentation++;
        m_nb_scene_nodes.push(0);
    }
}

const std::vector<mi::neuraylib::Tag>& Scene_logger::get_color_maps() const
{
    return m_color_maps;
}


}} // namespace index_common / namespace nv
