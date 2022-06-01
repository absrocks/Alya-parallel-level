/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/

#include "command_handler_scene_element.h"

#ifdef LINUX
#include <sys/types.h>
#include <dirent.h>
#endif

#include <set>
#include <sstream>
#include <stack>

#include <nv/index/irendering_kernel_programs.h>
#include "scene_element_importer.h"
#include "scene_utility.h"

using namespace nv::index_common;

namespace {

std::string quote(const std::string& s)
{
    return '"' + escape_JSON(s) + '"';
}

std::string get_name(
    mi::base::Handle<const nv::index::IScene_element>& scene_element,
    mi::neuraylib::Tag                                 scene_element_tag,
    mi::neuraylib::IDice_transaction*                  dice_transaction)
{
    std::string name;
    const std::string class_name = scene_element->get_class_name();
    const char* tag_name = dice_transaction->tag_to_name(scene_element_tag);
    if (tag_name != 0)
    {
        name = tag_name;

        // Skip scene elements used internally by the viewer
        if (name.compare(0, 5, "__nv_") == 0)
            return "";
    }

    if (name.empty())
    {
        // No name assigned: use the class name
        name = class_name;

        if (name == "Interpreter_scene")
            name = "Scene";
    }

    return name;
}

bool is_internal_scene_element(
    mi::neuraylib::Tag                scene_element_tag,
    mi::neuraylib::IDice_transaction* dice_transaction)
{
    const char* tag_name = dice_transaction->tag_to_name(scene_element_tag);
    if (tag_name != 0)
    {
        if (std::string(tag_name).compare(0, 5, "__nv_") == 0)
            return true;
    }
    return false;
}

// Traverses the scene hierarchy and creates a tree and a list representation of all elements in
// JSON format.
class Scene_traverse_list_JSON : public Scene_traverse_strategy_if
{
public:
    void set_filter(const std::vector<std::string>& filter)
    {
        m_filter = filter;
    }

    void run_before_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction)
    {
        assert(scene_element_tag.is_valid());
        assert(dice_transaction != 0);

        mi::base::Handle<const nv::index::IScene_element> scene_element(
            dice_transaction->access<const nv::index::IScene_element>(scene_element_tag));

        const std::string name = get_name(scene_element, scene_element_tag, dice_transaction);
        if (name.empty())
            return;

        std::ostringstream label;
        label << name << " (" << scene_element_tag << ")";

        // Add a marker to the label if the scene element is localized to the current scope
        mi::base::Handle<mi::neuraylib::IScope> scope(dice_transaction->get_scope());
        const mi::Sint32 scope_privacy_level = static_cast<mi::Sint32>(scope->get_privacy_level());
        const mi::Sint32 elem_privacy_level  = dice_transaction->get_privacy_level(scene_element_tag);
        if (elem_privacy_level > 0 && elem_privacy_level == scope_privacy_level)
        {
            label << " [local]";
        }

        if (!m_filter.empty())
        {
            const std::string class_name = scene_element->get_class_name();
            bool found = false;
            for (size_t i=0; i < m_filter.size(); ++i)
            {
                if (m_filter[i] == class_name)
                {
                    found = true;
                    break;
                }
            }

            if (found)
            {
                if (m_os_list.tellp() > 0)
                    m_os_list << ",\n";
                m_os_list << "{"
                          << quote("label") << ":" << quote(label.str()) << ","
                          << quote("tag") << ":" << scene_element_tag
                          << "}";
            }

            return;
        }

        if (m_parent_stack.empty())
            m_current_parent_tag = mi::neuraylib::Tag();
        else
            m_current_parent_tag = m_parent_stack.top();

        mi::base::Handle<const nv::index::IScene_group> scene_group(
            scene_element->get_interface<nv::index::IScene_group>());

        if (m_os_list.tellp() > 0)
            m_os_list << ",\n";

        m_os_list << "{"
                  << quote("label") << ":" << quote(label.str()) << ","
                  << quote("tag") << ":" << scene_element_tag << ","
                  << quote("parent_tag") << ":" << m_current_parent_tag
                  << "}";

        std::ostringstream os;
        if (!m_child_count.empty() && m_child_count.top() > 0 && m_os_tree.tellp() > 0)
            m_os_tree << ",";
        if (m_os_tree.tellp() > 0)
            m_os_tree << "\n";

        m_os_tree << "{"
                  << quote("label") << ":" << quote(label.str()) << ","
                  << quote("tag") << ":" << scene_element_tag << ","
                  << quote("parent_tag") << ":" << m_current_parent_tag;

        if (scene_group.is_valid_interface())
            m_os_tree << "," << quote("children") << ": [";
        else
            m_os_tree << "}";

        if (!m_child_count.empty())
            m_child_count.top()++;

        if (scene_group.is_valid_interface())
        {
            m_parent_stack.push(scene_element_tag);
            m_child_count.push(0);
        }
    }

    void run_after_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction)
    {
        if (!m_filter.empty())
            return;

        mi::base::Handle<const nv::index::IScene_element> scene_element(
            dice_transaction->access<const nv::index::IScene_element>(scene_element_tag));

        const std::string name = get_name(scene_element, scene_element_tag, dice_transaction);
        if (name.empty())
            return;

        mi::base::Handle<const nv::index::IScene_group> scene_group(
            scene_element->get_interface<nv::index::IScene_group>());
        if (scene_group.is_valid_interface())
        {
            m_parent_stack.pop();
            m_child_count.pop();
            m_os_tree << "\n]}";
        }
    }

    std::string get_list() const
    {
        return "[\n" + m_os_list.str() + "]\n";
    }

    std::string get_tree() const
    {
        return m_os_tree.str();
    }

private:
    std::vector<std::string>       m_filter;
    std::stack<mi::neuraylib::Tag> m_parent_stack;
    std::stack<std::string>        m_children;
    std::stack<int>                m_child_count;

    std::ostringstream             m_os_list;
    std::ostringstream             m_os_tree;

    mi::neuraylib::Tag             m_current_parent_tag;
};

// Traverses the scene hierarchy to ensure that the given identifier is uniuqe
class Scene_traverse_unique_identifier : public Scene_traverse_strategy_if
{
public:
    Scene_traverse_unique_identifier(const std::string& id)
      : m_id(id)
    {
    }

    void run_before_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction)
    {
        mi::base::Handle<const nv::index::IScene_element> scene_element(
            dice_transaction->access<const nv::index::IScene_element>(scene_element_tag));

        m_ids.insert(get_name(scene_element, scene_element_tag, dice_transaction));
    }

    std::string get() const
    {
        std::string id = m_id;
        mi::Sint32 suffix = 1;
        while (true)
        {
            if (m_ids.find(id) == m_ids.end())
                return id;

            std::ostringstream os;
            os << m_id << "_" << suffix;
            id = os.str();
            suffix++;
        }
    }

private:
    const std::string     m_id;
    std::set<std::string> m_ids;
};

// Traverses the scene hierarchy and removes all instances of the given scene element
class Scene_traverse_remove : public Scene_traverse_strategy_if
{
public:
    Scene_traverse_remove(mi::neuraylib::Tag remove_tag)
      : m_remove_tag(remove_tag)
    {
    }

    void run_before_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction)
    {
        if (is_internal_scene_element(scene_element_tag, dice_transaction))
            return;

        if (scene_element_tag == m_remove_tag)
        {
            m_containing_groups.push_back(m_parent_stack.top());
        }

        mi::base::Handle<const nv::index::IScene_group> scene_group(
            dice_transaction->access<nv::index::IScene_group>(scene_element_tag));
        if (scene_group.is_valid_interface())
        {
            m_parent_stack.push(scene_element_tag);
        }
    }

    void run_after_recursion(
        mi::neuraylib::Tag                scene_element_tag,
        mi::neuraylib::IDice_transaction* dice_transaction)
    {
        mi::base::Handle<const nv::index::IScene_group> scene_group(
            dice_transaction->access<nv::index::IScene_group>(scene_element_tag));
        if (scene_group.is_valid_interface())
        {
            m_parent_stack.pop();
        }

        if (m_parent_stack.empty())
        {
            for (size_t i=0; i < m_containing_groups.size(); ++i)
            {
                mi::base::Handle<nv::index::IScene_group> scene_group(
                    dice_transaction->edit<nv::index::IScene_group>(m_containing_groups[i]));
                if (scene_group)
                    scene_group->remove(m_remove_tag, dice_transaction);
            }
        }
    }

private:
    mi::neuraylib::Tag              m_remove_tag;
    std::stack<mi::neuraylib::Tag>  m_parent_stack;

    std::vector<mi::neuraylib::Tag> m_containing_groups;
};

//
// Handling of dataset description files
//

std::string get_dataset_description_file_suffix(const std::string& type)
{
    std::string suffix;
    if (type == "volume")
        suffix = "vol";
    else if (type == "irregular_volume")
        suffix = "ivol";
    else if (type == "horizon" || type == "heightfield")
        suffix = "hor";
    else if (type == "triangle_mesh")
        suffix = "tri";
    else
        return "";

    return "." + suffix + ".prj";
}

void list_regular_files(
    const std::string&        directory,
    const std::string&        suffix,
    std::vector<std::string>& files)
{
#ifdef LINUX
    DIR* dir = opendir(directory.c_str());
    if (dir)
    {
        dirent* entity;
        while ((entity = readdir(dir)) != 0)
        {
            // Only use regular files, but no symlinks (for security)
            if (entity->d_type & DT_REG)
            {
                // Only use files with the correct suffix
                const std::string name = entity->d_name;
                if (name.size() > suffix.size() &&
                    name.compare(name.size() - suffix.size(), suffix.size(), suffix) == 0)
                {
                    files.push_back(name);
                }
            }
        }
    }
    else
    {
        ERROR_LOG << "Could not open the directory '" << directory << "'";
    }
#else
        ERROR_LOG << "list_regular_files is not implemented for the current platform!";
#endif // LINUX
}

bool validate_dataset_description_file(
    const std::string& directory,
    const std::string& description_file,
    const std::string& type)
{
    std::vector<std::string> files;
    list_regular_files(directory, get_dataset_description_file_suffix(type), files);

    for (size_t i=0; i < files.size(); ++i)
    {
        if (files[i] == description_file)
            return true;
    }

    return false;
}

} // namespace

bool Command_handler_scene_element::handle(
    const std::string&    cmd,
    Call_event_arguments& args)
{
    mi::base::Handle<mi::neuraylib::IScope> scope = m_irc.get_current_scope();

    if (cmd == "get_scene_graph")
    {
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            scope->create_transaction<mi::neuraylib::IDice_transaction>());

        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<const nv::index::ISession>(m_irc.m_session_tag));

        // Retrieve a list of all scene elements
        Scene_traverse_list_JSON list_scene;

        if (args.has("filter"))
        {
            std::vector<std::string> filter;
            std::stringstream s(args.get("filter"));
            std::string item;
            while (std::getline(s, item, ',')) {
                filter.push_back(item);
            }
            list_scene.set_filter(filter);
        }

        traverse_scene_hierarchy_by_tag(session->get_scene(), &list_scene, dice_transaction.get());

        args.set("cmd", cmd);

        if (args.has("filter"))
            args.set("elements", list_scene.get_list());
        else
            args.set("tree", list_scene.get_tree());

        mi::base::Handle<mi::neuraylib::IScope> scope(dice_transaction->get_scope());

        std::ostringstream os;
        os << static_cast<int>(scope->get_privacy_level());
        args.set("scope_privacy_level", os.str().c_str());

        args.pass_through("select_tag");
        args.pass_through("update_scene_element");

        dice_transaction->commit();
        return true;
    }
    else if (cmd == "get_scene_element")
    {
        //TODO: The session exporter should be modified so that the scene update lock is not needed
        mi::base::Lock::Block block_update(&Nvindex_AppData::instance()->m_scene_update_lock);

        mi::neuraylib::Tag tag;
        std::istringstream is(args.get("tag"));
        is >> tag.id;

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            scope->create_transaction<mi::neuraylib::IDice_transaction>());

        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc.m_session_tag));

        mi::base::Handle<mi::neuraylib::IFactory> factory(
            m_irc.m_iindex_if->get_api_component<mi::neuraylib::IFactory>());
        mi::base::Handle<mi::IString> str(factory->create<mi::IString>("String"));

        session->export_scene_element(
            tag,
            nv::index::ISession::EXPORT_FORMAT_JSON |
            nv::index::ISession::EXPORT_VERBOSE  |
            nv::index::ISession::EXPORT_DEBUG  |
            nv::index::ISession::EXPORT_USER_DATA |
            nv::index::ISession::EXPORT_HINTS,
            str.get(),
            dice_transaction.get());

        args.set("cmd", cmd);
        args.set("tag", args.get("tag"));
        args.set("element", str->get_c_str());
        args.pass_through("previous_tag");

        dice_transaction->commit();
        return true;
    }
    else if (cmd == "edit_scene_element")
    {
        // First aquire the scene edit and update locks
        //TODO: The session exporter should be modified so that the scene update lock is not needed
        mi::base::Lock::Block block_edit(&Nvindex_AppData::instance()->m_scene_edit_lock);
        mi::base::Lock::Block block_update(&Nvindex_AppData::instance()->m_scene_update_lock);

        mi::neuraylib::Tag tag;
        std::istringstream is(args.get("tag"));
        is >> tag.id;

        args.set("cmd", cmd);
        args.set("element_tag", args.get("tag"));

        args.pass_through("reload_scene_graph");

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            scope->create_transaction<mi::neuraylib::IDice_transaction>());

        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc.m_session_tag));

        const std::string key = args.get("property");

        // Retrieve the scene element name (unique identifier)
        std::string elem_name;
        if (!args.get("id").empty())
            elem_name = args.get("id");

        if (elem_name.empty())
        {
            elem_name = "__unknown__";
        }
        else
        {
            // Query the tag assigned to the identifier and use that one instead of the passed tag.
            // This shouldn't make any difference during interactive use, but is required for
            // playing back regression test recordings, since the element tags may change between
            // multiple runs of the application.
            mi::neuraylib::Tag id_tag = dice_transaction->name_to_tag(elem_name.c_str());
            if (id_tag)
            {
                tag = id_tag;
            }
        }

        // Special handling for editing the root scene group
        if (tag == session->get_scene())
        {
            const std::string val = args.get("value");

            {
                mi::base::Handle<nv::index::IScene> scene(
                    dice_transaction->edit<nv::index::IScene>(session->get_scene()));

                if (key == "region_of_interest")
                {
                    scene->set_clipped_bounding_box(get_bbox_float32_3(val));
                }
                else if (key == "enabled")
                {
                    scene->set_enabled(get_bool(val));
                }
                else if (key == "children")
                {
                    scene->clear(dice_transaction.get()); // remove all existing children
                    std::istringstream children_list(val);
                    while (!children_list.eof())
                    {
                        mi::neuraylib::Tag tag;
                        children_list >> tag.id;
                        if (tag.is_valid())
                            scene->append(tag, dice_transaction.get());
                    }
                }
            }
            dice_transaction->commit();
            return true;
        }

        mi::base::Handle<mi::neuraylib::IFactory> factory(
            m_irc.m_iindex_if->get_api_component<mi::neuraylib::IFactory>());
        mi::base::Handle<mi::IString> str(factory->create<mi::IString>("String"));

        // Export the single scene element specified by the tag
        session->export_scene_element(
            tag,
            nv::index::ISession::EXPORT_FORMAT_PRJ | nv::index::ISession::EXPORT_USER_DATA,
            str.get(),
            dice_transaction.get());

        std::ostringstream changes;

        // Add the requested extra property before the other changes
        if (args.has("extra_prop_key") && args.has("extra_prop_value"))
        {
            changes << args.get("extra_prop_key") << " = " << args.get("extra_prop_value") << "\n";
        }

        if (args.has("keys"))
        {
            const mi::Sint32 size = get_sint32(args.get("size"));

            std::stringstream os_size(args.get("size_props"));
            std::string size_prop;
            while (std::getline(os_size, size_prop, ' '))
            {
                changes << size_prop << " = " << size << "\n";
            }

            std::stringstream os_keys(args.get("keys"));
            std::string key;
            mi::Sint32 i = 0;
            while (std::getline(os_keys, key, ' '))
            {
                for (mi::Sint32 j = 0; j < size; ++j)
                {
                    std::ostringstream input_key;
                    input_key << "value" << i << "_" << j;
                    changes << key << j << " = " << args.get(input_key.str()) << "\n";
                }
                i++;
            }
        }
        else
        {
            changes << key << " = " << args.get("value") << "\n";
        }

        // Construct the scene element as the exported state plus the newly set property, which
        // will overwrite the existing one (in .prj format)
        const std::string prj = str->get_c_str() + changes.str();

        // Read the .prj string into a String_dict
        String_dict dict;
        std::istringstream prj_is(prj);
        dict.read(prj_is);

        // When an expensive property was changed explicitly then allow reloading of data
        if (args.has("expensive"))
        {
            dict.insert("allow_reloading", "yes");
        }

        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<nv::index::IScene>(session->get_scene()));

        // Dummy data, only required for handling references
        String_dict& scene_dict = dict;
        std::map<std::string, mi::neuraylib::Tag> name_to_tag;
        std::vector<std::string> referencing_elements;

        // Run the scene importer, passing the current scene element tag, so that it will be
        // overwritten instead of creating a new element
        Scene_element_importer importer(
            dice_transaction, scene.get(), m_irc.m_session_tag, name_to_tag, referencing_elements, scene_dict);
        bool is_group;
        importer.import(elem_name, dict.get("type"), dict, is_group, tag);

        dice_transaction->commit();

        return true;
    }
    else if (cmd == "localize_scene_element")
    {
        // Localize the scene element in the current scope
        mi::base::Lock::Block block_edit(&Nvindex_AppData::instance()->m_scene_edit_lock);
        mi::base::Lock::Block block_update(&Nvindex_AppData::instance()->m_scene_update_lock);

        args.set("cmd", cmd);

        args.pass_through("element_tag");
        args.pass_through("reload_scene_graph");

        mi::neuraylib::Tag tag;
        tag.id = get_uint32(args.get("element_tag"));

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            scope->create_transaction<mi::neuraylib::IDice_transaction>());

        dice_transaction->localize(tag, mi::neuraylib::IDice_transaction::LOCAL_SCOPE);

        dice_transaction->commit();

        return true;
    }
    else if (cmd == "move_scene_element")
    {
        // Move scene element from one scene group to another
        mi::base::Lock::Block block_edit(&Nvindex_AppData::instance()->m_scene_edit_lock);
        mi::base::Lock::Block block_update(&Nvindex_AppData::instance()->m_scene_update_lock);

        args.set("cmd", cmd);

        args.pass_through("element_tag");
        args.pass_through("reload_scene_graph");

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            scope->create_transaction<mi::neuraylib::IDice_transaction>());
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc.m_session_tag));

        std::string error;

        // Edit the source and destination group. When the scene element was moved inside its parent
        // group then there will only be one group parameter.
        static const char* NUMS[] = { "1", "2" };
        for (mi::Sint32 i = 0; i < 2; ++i)
        {
            const std::string key = "group" + std::string(NUMS[i]) + "_tag";
            if (!args.has(key))
                continue;

            mi::neuraylib::Tag group_tag;
            std::istringstream is(args.get(key));
            is >> group_tag.id;

            mi::base::Handle<nv::index::IScene_group> group(
                dice_transaction->edit<nv::index::IScene_group>(group_tag));
            if (group)
            {
                // Set the updated list of children
                group->clear(dice_transaction.get());

                std::stringstream s(args.get("group" + std::string(NUMS[i]) + "_children"));
                std::string child;
                while (std::getline(s, child, ' '))
                {
                    mi::neuraylib::Tag child_tag;
                    child_tag.id = get_uint32(child);
                    if (!group->append(child_tag, dice_transaction.get()))
                    {
                        std::ostringstream os;
                        os << "Could not append scene element " << child_tag << " to scene group " << group_tag << ". "
                           << "This most likely failed because you were trying to add a massive scene element to a "
                           << "Transformed_scene_group instead of a Static_scene_group.";
                        error = os.str();
                        break;
                    }
                }
            }
            else
            {
                std::ostringstream os;
                os << "Could not edit scene group " << group_tag;
                        error = os.str();
                break;
            }
        }

        if (error.empty())
        {
            dice_transaction->commit();
        }
        else
        {
            ERROR_LOG << error;
            args.set("error", error);
            dice_transaction->abort();
        }

        return true;
    }
    else if (cmd == "add_scene_element")
    {
        mi::base::Lock::Block block_edit(&Nvindex_AppData::instance()->m_scene_edit_lock);
        mi::base::Lock::Block block_update(&Nvindex_AppData::instance()->m_scene_update_lock);

        mi::neuraylib::Tag group_tag;
        group_tag.id = get_uint32(args.get("parent_group"));

        // Create the scene element either in the current local or in global scope, depending on
        // whether the parent scene group is localized or not.
        // This prevents an invalid tag access when the scene element is created in the local scope
        // added to a scene group that is not localized, and the element's tag is unknown when
        // switching to global scope.
        mi::neuraylib::IScope* create_scope;
        {
            // Retrieve the privacy level of the parent scene group
            mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
                scope->create_transaction<mi::neuraylib::IDice_transaction>());
            const mi::Sint32 group_privacy_level = dice_transaction->get_privacy_level(group_tag);
            dice_transaction->commit();

            if (group_privacy_level > 0)
            {
                // Parent scene group is localized, so create new scene element in local scope
                create_scope = scope.get();
            }
            else
            {
                // Parent scene group is in global scope, so also create new scene element there
                create_scope = m_irc.m_global_scope.get();
            }
        }

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            create_scope->create_transaction<mi::neuraylib::IDice_transaction>());

        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc.m_session_tag));
        mi::base::Handle<const nv::index::IScene> scene(
            dice_transaction->access<nv::index::IScene>(session->get_scene()));

        mi::neuraylib::Tag elem_tag;
        if (args.has("reference"))
            elem_tag.id = get_uint32(args.get("reference"));

        if (!elem_tag)
        {
            std::string type = args.get("type");

            Scene_traverse_unique_identifier unique(args.get("name"));
            traverse_scene_hierarchy_by_tag(session->get_scene(), &unique, dice_transaction.get());
            const std::string elem_name = unique.get();

            String_dict dict;

            // Add dataset description entry to dictionary if specified (for massive scene elements)
            if (!args.get("dataset_description").empty())
            {
                const std::string description_file = args.get("dataset_description");
                // Do not blindly trust the user input, but check that the specified file is in the
                // list of valid description files
                const std::string dir = Nvindex_AppData::instance()->peek_app_proj()->get(
                    "app::dataset_description_directory");
                if (!dir.empty())
                {
                    if (validate_dataset_description_file(dir, description_file, type))
                        dict.insert("description_file", dir + "/" + description_file);
                    else
                        ERROR_LOG << "Received invalid dataset description file '" << description_file << "'";
                }
            }

            // Use nicer default material than the "no material defined" pink
            if (type == "phong_gl")
            {
                dict.insert("ambient",   "0.3 0.3 0.3");
                dict.insert("diffuse",   "0.4 0.4 0.4");
                dict.insert("specular",  "0.4 0.4 0.4");
                dict.insert("shininess", "100");
            }

            // Dummy data, only required for handling references
            String_dict& scene_dict = dict;
            std::map<std::string, mi::neuraylib::Tag> name_to_tag;
            std::vector<std::string> referencing_elements;

            // Run the scene importer
            Scene_element_importer importer(
                dice_transaction, scene.get(), m_irc.m_session_tag, name_to_tag, referencing_elements, scene_dict);
            bool is_group;
            elem_tag = importer.import(elem_name, type, dict, is_group, mi::neuraylib::Tag());

            if (!elem_tag)
            {
                ERROR_LOG << "Failed to create scene element of type '" << type << "'";
                args.set("error", "Failed to create scene element of type '" + type + "'.");
            }
        }

        if (elem_tag == group_tag)
        {
            ERROR_LOG << "Trying to append group with tag " << group_tag << " to itself, ignored.";
            elem_tag = mi::neuraylib::Tag();
        }

        if (elem_tag)
        {
            bool success = false;
            mi::base::Handle<nv::index::IScene_group> group(
                dice_transaction->edit<nv::index::IScene_group>(group_tag));
            if (group)
            {
                group->clear(dice_transaction.get());

                std::stringstream s(args.get("child_order"));
                std::string child;
                while (std::getline(s, child, ' '))
                {
                    mi::neuraylib::Tag child_tag;
                    child_tag.id = get_uint32(child);
                    if (child_tag.is_valid())
                        group->append(child_tag, dice_transaction.get());
                    else
                        // Null tag, insert the new element
                        success = group->append(elem_tag, dice_transaction.get());
                }
            }

            if (!success)
            {
                ERROR_LOG << "Could not append to group " << group_tag;
                args.set("error", "Could not append the new scene element to the scene group.");
            }
        }

        dice_transaction->commit();

        std::ostringstream os;

        if (args.has("reference"))
            os << group_tag; // select parent scene group
        else
            os << elem_tag;  // select new scene element
        args.set("cmd", cmd);
        args.set("element_tag", os.str());

        return true;
    }
    else if (cmd == "remove_scene_element")
    {
        mi::base::Lock::Block block_edit(&Nvindex_AppData::instance()->m_scene_edit_lock);
        mi::base::Lock::Block block_update(&Nvindex_AppData::instance()->m_scene_update_lock);

        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            scope->create_transaction<mi::neuraylib::IDice_transaction>());
        mi::base::Handle<const nv::index::ISession> session(
            dice_transaction->access<nv::index::ISession>(m_irc.m_session_tag));

        mi::neuraylib::Tag remove_tag;
        remove_tag.id = get_uint32(args.get("tag"));

        Scene_traverse_remove remover(remove_tag);
        traverse_scene_hierarchy_by_tag(session->get_scene(), &remover, dice_transaction.get());

        dice_transaction->commit();

        args.set("cmd", cmd);

        return true;
    }
    else if (cmd == "get_scene_element_types")
    {
        args.set("cmd", cmd);

        std::multimap<std::string, std::string> types;
        Scene_element_importer::get_supported_types(types);

        std::ostringstream os;
        os << "[";
        std::multimap<std::string, std::string>::const_iterator it = types.begin();
        while (it != types.end())
        {
            if (it != types.begin())
                os << ",";
            os << "{" << quote("label") << ":" << quote(it->first) << ", "
               << quote("children") << ":[\n";
            std::multimap<std::string, std::string>::const_iterator current = it;
            while (it != types.end() && it->first == current->first)
            {
                if (it != current)
                    os << ",";
                os << "{" << quote("label") << ":" << quote(it->second) << ", "
                   << quote("parent") << ":" << quote(it->first) << "}";
                ++it;
            }
            os << "]}";
        }
        os << "]";

        args.set("types", os.str());
        return true;
    }
    else if (cmd == "get_dataset_descriptions")
    {
        args.set("cmd", cmd);

        const std::string type = args.get("type");
        const std::string suffix = get_dataset_description_file_suffix(type);
        if (suffix.empty())
        {
            ERROR_LOG << "Unsupported type for get_dataset_descriptions: '" << type << "'";
            return true;
        }

        // Look in this directory for description files (or current director if not specified)
        const std::string dir = Nvindex_AppData::instance()->peek_app_proj()->get("app::dataset_description_directory");
        std::vector<std::string> files;
        if (!dir.empty())
            list_regular_files(dir, suffix, files);

        std::string s;
        for (size_t i=0; i < files.size(); ++i)
        {
            if (i > 0)
                s += "/"; // Slash used as entry separator
            s += files[i];
        }
        args.set("files", s);
        return true;
    }
    else if (cmd == "get_is_scene_initialized")
    {
        args.set("cmd", cmd);
        args.set("initialized", "true");
        args.set("scene_setup_done", "true");
        return true;
    }
    else if (cmd == "update_source_code")
    {
        // lock
        mi::base::Lock::Block block_edit(&Nvindex_AppData::instance()->m_scene_edit_lock);
        mi::base::Lock::Block block_update(&Nvindex_AppData::instance()->m_scene_update_lock);

        // must have: tag, source_code
        assert(args.has("tag"));
        assert(args.has("source_code"));

        // get tag
        mi::neuraylib::Tag tag;
        std::istringstream is(args.get("tag"));
        is >> tag.id;
        assert(tag.is_valid());
        // INFO_LOG << "update_source_code: tag: " << tag;

        // results
        args.set("cmd", cmd);
        args.set("element_tag", args.get("tag"));
        args.set("exit_value", "ok"); // Depends on the compile result.
        args.set("error_message", "");

        // get source code
        std::string source_code(args.get("source_code"));
        if (source_code.empty())
        {
            ERROR_LOG << "Empty source code. ignored.";
            // Set the error status.
            args.set("exit_value", "error");
            args.set("error_message", "Empty source code.");
            
            return true;        // command has been accepted
        }

        // set the code to the attribute
        mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
            scope->create_transaction<mi::neuraylib::IDice_transaction>());
        {
            mi::base::Handle<nv::index::IRendering_kernel_program> rkc(
                dice_transaction->edit<nv::index::IRendering_kernel_program>(tag));
            assert(rkc.is_valid_interface());
            rkc->set_program_source(source_code.c_str());
        }
        dice_transaction->commit();

        return true;
    }

    return false;
}
