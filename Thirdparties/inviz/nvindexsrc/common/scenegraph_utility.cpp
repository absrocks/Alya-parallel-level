/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "scenegraph_utility.h"

#include <nv/index/isession.h>
#include <nv/index/iscene.h>
#include <nv/index/iscene_group.h>
#include <nv/index/iregular_volume.h>

#include "common_utility.h"
#include "forwarding_logger.h"

#include <cassert>

namespace nv {
namespace index_common {

//----------------------------------------------------------------------
mi::neuraylib::Tag get_scene_group(
    const mi::neuraylib::Tag&         session_tag,
    const std::string&                group_name,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    if(dice_transaction == 0)
    {
        ERROR_LOG << "Invalid DiCE transaction (file: " << __FILE__ ", line: " << __LINE__ << ").";
        return mi::neuraylib::NULL_TAG;
    }

    if(group_name.empty())
    {
        ERROR_LOG << "No scenegraph name defined (file: " << __FILE__ ", line: " << __LINE__ << ").";
        return mi::neuraylib::NULL_TAG;
    }
    
    mi::neuraylib::Tag group_tag = dice_transaction->name_to_tag(group_name.c_str());
    if(!group_tag.is_valid())
    {
        INFO_LOG << "No scene group stored with the name " << group_name << ". Creating a new scene group.";
        group_tag = create_scene_group(session_tag, group_name, dice_transaction);
    }
    
    mi::base::Handle<const nv::index::ITransformed_scene_group> group(
        dice_transaction->access<const nv::index::ITransformed_scene_group>(group_tag));
    if(!group.is_valid_interface())
    {
        INFO_LOG << "The database element stored with the name " << group_name << " is not a valid scene group.";
        return mi::neuraylib::NULL_TAG;
    }

    return group_tag;
}

mi::neuraylib::Tag create_scene_group(
    const mi::neuraylib::Tag&         session_tag,
    const std::string&                group_name,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    if(dice_transaction == 0)
    {
        ERROR_LOG << "Invalid DiCE transaction (file: " << __FILE__ ", line: " << __LINE__ << ").";
        return mi::neuraylib::NULL_TAG;
    }

    if(!session_tag.is_valid())
    {
        ERROR_LOG << "Invalid session tag (file: " << __FILE__ ", line: " << __LINE__ << ").";
        return mi::neuraylib::NULL_TAG;
    }

    mi::base::Handle<const nv::index::ISession> session(
        dice_transaction->access<nv::index::ISession>(session_tag));
    assert(session.is_valid_interface());

    mi::base::Handle<nv::index::IScene> scene(
        dice_transaction->edit<nv::index::IScene>(session->get_scene()));
    assert(scene.is_valid_interface());

    mi::base::Handle<nv::index::ITransformed_scene_group> transformed(
        scene->create_scene_group<nv::index::ITransformed_scene_group>());
        
    mi::neuraylib::Tag group_tag = dice_transaction->store_for_reference_counting(
            transformed.get(),
            mi::neuraylib::NULL_TAG,
            group_name.c_str());
    scene->append(group_tag, dice_transaction);
    
    return group_tag;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
