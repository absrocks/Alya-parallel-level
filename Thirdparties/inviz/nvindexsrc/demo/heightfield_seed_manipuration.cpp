/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file

#include "heightfield_seed_manipuration.h"

#include <nv/index/isession.h>
#include <nv/index/iscene.h>
#include <nv/index/iregular_heightfield.h>

#include <cassert>
#include <cstdio>
#include <ctime>
#include <vector>
#include <map>

#include "common/forwarding_logger.h"
#include "common/string_dict.h"

#include "nvindex_appdata.h"



//----------------------------------------------------------------------
void import_seed_lines(
    mi::base::Handle<mi::neuraylib::IScope>&    scope,
    const mi::neuraylib::Tag&                   session_tag,
    mi::Uint32                                  heightfield_id,
    const std::string&                          seed_lines_file_name)
{
    if(seed_lines_file_name.size()==0)
    {
        INFO_LOG << "No seed lines file";
        return;
    }

    FILE* seed_lines_file = fopen(seed_lines_file_name.c_str(), "r");
    if(!seed_lines_file)
    {
        INFO_LOG << "Seed lines file '" << seed_lines_file_name.c_str() << "' cannot be opened for reading.";
        return;
    }

    std::map<mi::Float32, std::vector<mi::math::Vector_struct<mi::Float32, 3> > > seed_lines;
    std::vector<mi::math::Vector_struct<mi::Float32, 3> > seed_points;

    char line[80];
    while(fgets(line, sizeof(line), seed_lines_file))
    {
        const std::string line_in = std::string(&line[0]);
        if(line_in[0]=='m')
            continue;

        std::string::size_type pos_0 = line_in.find(" ");
        const std::string val_0_str = line_in.substr(0, pos_0);

        const std::string rest = line_in.substr(pos_0 + 1);
        std::string::size_type pos_1 = rest.find(" ");

        const std::string val_1_str = rest.substr(0, pos_1);
        const std::string val_2_str = rest.substr(pos_1+1);

        mi::math::Vector_struct<mi::Float32, 3> vector;
        vector.x = nv::index_common::get_float32(val_0_str);
        vector.y = nv::index_common::get_float32(val_1_str);
        vector.z = nv::index_common::get_float32(val_2_str);

        seed_points.push_back(vector);

        std::map<mi::Float32, std::vector<mi::math::Vector_struct<mi::Float32, 3> > >::iterator itr_find = seed_lines.find(vector.y);
        if(itr_find==seed_lines.end())
        {
            std::vector<mi::math::Vector_struct<mi::Float32, 3> > seed_line(0);
            seed_line.push_back(vector);
            seed_lines[vector.y] = seed_line;
        }
        else
        {
            std::vector<mi::math::Vector_struct<mi::Float32, 3> >& seed_line = itr_find->second;
            seed_line.push_back(vector);
        }
    }
    fclose(seed_lines_file);

    // get dice transaction
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        const std::vector<mi::neuraylib::Tag> hf_tag_vec = Nvindex_AppData::instance()->get_heightfield_tag_vec();
        assert(heightfield_id < hf_tag_vec.size());

        const mi::neuraylib::Tag& heightfield_tag = hf_tag_vec.at(heightfield_id);
        {   // Editing the heightfield dataset
            mi::base::Handle<nv::index::IRegular_heightfield> heightfield(
                dice_transaction->edit<nv::index::IRegular_heightfield>(heightfield_tag));
            assert(heightfield.is_valid_interface());
            if (!(heightfield.is_valid_interface()))
            {
                ERROR_LOG << "No heightfield found with tag " << heightfield_tag.id;
                return;
            }

            std::map<mi::Float32, std::vector<mi::math::Vector_struct<mi::Float32, 3> > >::iterator itr = seed_lines.begin();
            mi::Sint32 nb_seed_lines = 0;
            mi::Sint32 nb_seed_line_vertices = 0;
            for(; itr!=seed_lines.end(); ++itr)
            {
                std::vector<mi::math::Vector_struct<mi::Float32, 3> >& seed_line = itr->second;
                const mi::Uint32 nb_seed_line_values = seed_line.size();
                for(mi::Uint32 i=0; i<nb_seed_line_values; ++i)
                {
                    INFO_LOG << "Seed line value [" << seed_line[i].y << "]: "
                             << seed_line[i].x << "," << seed_line[i].y << "," << seed_line[i].z;
                }
                heightfield->add_seed_line(&seed_line[0], nb_seed_line_values);
                nb_seed_lines++;
                nb_seed_line_vertices += nb_seed_line_values;
            }

            // Adding seed points
            const mi::Uint32 nb_seed_points = seed_points.size();
            heightfield->add_seed_points(&seed_points[0], nb_seed_points);
            
            INFO_LOG << "Added " << nb_seed_lines << " seed lines with " << nb_seed_line_vertices
                     << " vertices and " << nb_seed_points << " seed points.";
        }
    }
    // commit changes
    dice_transaction->commit();
}

//----------------------------------------------------------------------
void remove_arbitrary_seed_line(
    mi::base::Handle<mi::neuraylib::IScope>&    scope,
    const mi::neuraylib::Tag&                   session_tag,
    mi::Uint32                                  heightfield_id)
{
    // DiCE transaction
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        const std::vector<mi::neuraylib::Tag> hf_tag_vec = Nvindex_AppData::instance()->get_heightfield_tag_vec();
        assert(heightfield_id < hf_tag_vec.size());
        
        // Heightfield scene element tag
        const mi::neuraylib::Tag& heightfield_tag = hf_tag_vec.at(heightfield_id);
        {   // Heightfield scene element for editing
            mi::base::Handle<nv::index::IRegular_heightfield> heightfield(
                dice_transaction->edit<nv::index::IRegular_heightfield>(heightfield_tag));
            assert(heightfield.is_valid_interface());
            
            // Heightfield name and number of seed lines
            const char* heightfield_name = (heightfield->get_name() != NULL ? heightfield->get_name() : "<unknown>");
            const mi::Uint32 nb_seed_lines = heightfield->get_nb_seed_lines();
            INFO_LOG << "Heightfield dataset " << heightfield_name << " has " << nb_seed_lines << " seed lines";
            
            // Access, print, and remove a random seed line.
            if(nb_seed_lines>0)
            {
                // Compute a random number in the range [0, nb_seed_lines-1] that indexes a seed line
                srand((int)(time(NULL)));
                mi::Uint32 random_line_index = rand() % nb_seed_lines;
                
                // Access the seed line that corresponds to the index/random number
                mi::Uint32 nb_seed_line_elements = 0;
                
                // Access ...
                const mi::math::Vector_struct<mi::Float32, 3>* seed_line = heightfield->get_seed_line(random_line_index, nb_seed_line_elements);
                
                // ... print ...
                INFO_LOG << "Seed line with index " << random_line_index << " has " << nb_seed_line_elements << " elements.";
                for(mi::Uint32 i=0; i<nb_seed_line_elements; ++i)
                {
                    const mi::math::Vector_struct<mi::Float32, 3> vertex = seed_line[i];
                    INFO_LOG << i << ":\t " << vertex;
                }
                delete[] seed_line;
                
                // ... remove
                INFO_LOG << "Removing line (index: " << random_line_index << ") from heightfield dataset " << heightfield_name;
                heightfield->remove_seed_line(random_line_index);
            }
        }
    }
    // commit changes
    dice_transaction->commit();
}
    
//----------------------------------------------------------------------
void remove_arbitrary_seed_point(
    mi::base::Handle<mi::neuraylib::IScope>&    scope,
    const mi::neuraylib::Tag&                   session_tag,
    mi::Uint32                                  heightfield_id)
{
    // DiCE transaction
    mi::base::Handle<mi::neuraylib::IDice_transaction> dice_transaction(
        scope->create_transaction<mi::neuraylib::IDice_transaction>());
    assert(dice_transaction.is_valid_interface());
    {
        const std::vector<mi::neuraylib::Tag> hf_tag_vec = Nvindex_AppData::instance()->get_heightfield_tag_vec();
        assert(heightfield_id < hf_tag_vec.size());

        // Heightfield scene element tag
        const mi::neuraylib::Tag& heightfield_tag = hf_tag_vec.at(heightfield_id);
        
        {   // Heightfield scene element for editing
            mi::base::Handle<nv::index::IRegular_heightfield> heightfield(
                dice_transaction->edit<nv::index::IRegular_heightfield>(heightfield_tag));
            assert(heightfield.is_valid_interface());
            
            // Heightfield name and number of seed points
            const char* heightfield_name = (heightfield->get_name() != NULL ? heightfield->get_name() : "<unknown>");
            const mi::Uint32 nb_seed_points = heightfield->get_nb_seed_points();
            INFO_LOG << "Heightfield dataset " << heightfield_name << " has " << nb_seed_points << " seed points";
            
            // Access, print, and remove a random seed point.
            if(nb_seed_points>0)
            {
                // Compute a random number in the range [0, nb_seed_points-1] that indexes a seed point
                srand((int)(time(NULL)));
                mi::Uint32 random_line_index = rand() % nb_seed_points;
                
                // Access the seed point that corresponds to the index/random number ...
                mi::math::Vector_struct<mi::Float32, 3> seed_point;
                heightfield->get_seed_point(random_line_index, seed_point);
                
                // ... print and remove
                INFO_LOG << "Removing point " << seed_point << " (index: " << random_line_index << ") from heightfield dataset " << heightfield_name;
                heightfield->remove_seed_point(random_line_index);
            }
        }
    }        
    // commit changes
    dice_transaction->commit();
}

//----------------------------------------------------------------------

