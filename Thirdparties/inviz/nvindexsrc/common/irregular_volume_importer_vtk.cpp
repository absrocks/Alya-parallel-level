/******************************************************************************
 * Copyright 1986, 2017 NVIDIA Corporation. All rights reserved.
 *****************************************************************************/
/// \file
/// \brief external irregular volume data importer job

#include "irregular_volume_importer_vtk.h"

#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cassert>
#include <limits>
#include <vector>
#include <map>
#include <stack>
#include <nv/index/iirregular_volume_subset.h>

#include "forwarding_logger.h"

namespace nv {
namespace index_common {
    
// split a long string with spaces into separated substring without spaces    
std::vector<std::string> string_split(std::string s)
{
    
    std::vector<std::string> output;
    std::string substring;
    std::istringstream iss(s);    
    
    while (!iss.eof())      
    {
        iss >> substring;
        output.push_back(substring);
    }
    
    return output;
}
 
// Face definitions for a VTK_HEXADRON
const mi::Uint32 FACE_UNDEFINED = ~0;
const mi::Uint32 FACE_OPPOSITE[6] = {1, 0, 3, 2, 5, 4};
const mi::Uint32 FACE_IDS[6][4] =
{
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {1, 2, 6, 5},
    {0, 3, 7, 4},
    {0, 1, 5, 4},
    {3, 2, 6, 7}
};

// This struct represent a single face of an Hexahedron
struct Hexahedron_face
{
    mi::Uint32 vids[4];
    
    Hexahedron_face(mi::Uint32 a, mi::Uint32 b, mi::Uint32 c, mi::Uint32 d)
    {
        // Add this code instead if it's needed faces ids to be sorted by lower vertice id.
        // mi::Uint32 tids[4] = {a, b, c, d};
        // mi::Uint32 offset = 0;
        // mi::Uint32 min_id = a;
        
        // for(mi::Uint32 i=1; i<4; i++)
        // {
            // if(tids[i] < min_id)
            // {
                // min_id = tids[i];
                // offset = i;
            // }
        // }
        
        // for(mi::Uint32 i=0; i<4; i++)
            // vids[i] = tids[(i+offset)%4];
        
        vids[0] = a, vids[1] = b, vids[2] = c, vids[3] = d;
    }
};

// Hexahedron face compare struct. Used to sort faces based its 4 vertice indices
struct Hexahedron_face_compare 
{
    bool operator() (const Hexahedron_face& lhs, const Hexahedron_face& rhs) const
    {
        if(lhs.vids[0] < rhs.vids[0])
            return true;
        else if (lhs.vids[0] > rhs.vids[0])
            return false;
        else
        {
            if(lhs.vids[1] < rhs.vids[1])
                return true;
            else if (lhs.vids[1] > rhs.vids[1])
                return false;
            else
            {
                if(lhs.vids[2] < rhs.vids[2])
                    return true;
                else if (lhs.vids[2] > rhs.vids[2])
                    return false;
                else
                {
                    if(lhs.vids[3] < rhs.vids[3])
                        return true;
                    else
                        return false;
                }
            }
        }
    }
};

// This function is used to calculate the cells(hexahedrons) orientation to be 'normal oriented' or 'inverse oriented'.
// This information is used later to decompose the hexahedrons into tetrahedrons and decide to use 
// a 'normal' or 'inverse' decomposition in order to get continuity between hexahedron faces. 
bool set_cell_orientation(
    std::vector<mi::Sint8>&                         cells_orientation, 
    const std::vector<std::vector<mi::Uint32> >&    hexahedron_adjacent_list)
{
    std::pair <std::multimap<Hexahedron_face,mi::Uint32>::const_iterator, std::multimap<Hexahedron_face,mi::Uint32>::const_iterator> ret;

    // invalid input data test
    if(cells_orientation.empty() || cells_orientation.size() != hexahedron_adjacent_list.size())
    {
        ERROR_LOG << "VTK IRV Importer| Cell orientation: Invalid input data.";
        return false;
    }
    
    std::stack<mi::Uint32> cells_stack;

    const mi::Uint32 total_cells = cells_orientation.size();
    mi::Uint32 calculated_faces = 0u;
    mi::Uint32 pivot = 0u;
    
    while(calculated_faces < total_cells)
    {
        // look for next sparse connected cells group
        for(mi::Uint32 i=pivot; i<total_cells; i++)
        {
            // define first cell to be normal oriented
            if(cells_orientation[i] == 0)
            {
                cells_orientation[i] = 1;
                cells_stack.push(i);
                pivot = i+1;
                break;
            }
        }
        
        while(!cells_stack.empty())
        {
            // get next cell to check its neighboor cells
            const mi::Uint32 cur_cell_id = cells_stack.top();
            cells_stack.pop();
            calculated_faces++;
            
            const mi::Sint8 cur_orientation = cells_orientation[cur_cell_id];
            const std::vector<mi::Uint32> &cur_adj_list = hexahedron_adjacent_list[cur_cell_id];
            
            for(mi::Uint32 k=0; k<cur_adj_list.size(); k++)
            {
                const mi::Uint32 next_cell_id = cur_adj_list[k];
                if(cells_orientation[next_cell_id] == 0)
                {
                    cells_orientation[next_cell_id] = -cur_orientation;
                    cells_stack.push(next_cell_id);
                }
                // this can only happen if the mesh is not concave or presents holes.
                else if(cells_orientation[next_cell_id] == cur_orientation)
                {
                    ERROR_LOG << "VTK IRV Importer| Cell orientation: Adjacent cell conflict.";
                    return false;
                }
            }
        }
    }
    
    return true;
}

//----------------------------------------------------------------------
Irregular_volume_importer_vtk::IVOL_ts_mesh::IVOL_ts_mesh()
  : bbox(mi::math::Vector<mi::Float32, 3>( (std::numeric_limits<float>::max)()),
         mi::math::Vector<mi::Float32, 3>(-(std::numeric_limits<float>::max)()))
  , attributes_value_range( (std::numeric_limits<float>::max)(),
                           -(std::numeric_limits<float>::max)())
{
}

mi::base::Lock                          Irregular_volume_importer_vtk::m_cached_input_lock;
std::string                             Irregular_volume_importer_vtk::m_cached_input_file;
Irregular_volume_importer_vtk::IVOL_ts_mesh Irregular_volume_importer_vtk::m_cached_input_ivol_mesh;

//----------------------------------------------------------------------
// constructor
Irregular_volume_importer_vtk::Irregular_volume_importer_vtk(const std::string& filename, bool use_cache)
{
    m_use_cache = use_cache;
    
    m_filename = filename;

    m_configuration = std::string() +
        "importer=vtk\n" +
        "input_file=" + m_filename + "\n";
}

//----------------------------------------------------------------------
// constructor
Irregular_volume_importer_vtk::Irregular_volume_importer_vtk()
{
    m_use_cache = false;
}

//----------------------------------------------------------------------
// destructor
Irregular_volume_importer_vtk::~Irregular_volume_importer_vtk()
{
}

//----------------------------------------------------------------------
mi::Size Irregular_volume_importer_vtk::estimate(
    const mi::math::Bbox_struct<mi::Float32, 3>&    bounding_box,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{
    return 0;
}

//----------------------------------------------------------------------
nv::index::IDistributed_data_subset* Irregular_volume_importer_vtk::create(
    const mi::math::Bbox_struct<mi::Float32, 3>&    bounding_box,
    nv::index::IData_subset_factory*                factory,
    mi::neuraylib::IDice_transaction*               dice_transaction) const
{    
    INFO_LOG << "VTK IRV Importer| VTK Irregular volume importer loads '" << m_filename << "', bounds: " << bounding_box;

    const std::string file_extension =
          m_filename.find_last_of(".") != std::string::npos
        ? m_filename.substr(m_filename.find_last_of(".") + 1)
        : "";
      
    if( file_extension == "" || !(file_extension == "vtk") )
    {
        ERROR_LOG << "VTK IRV Importer: Unknown file type: '" << m_filename << "'.";
        return NULL;
    }
    
    mi::base::Handle<nv::index::IIrregular_volume_subset>
        irregular_volume_subset(factory->create<nv::index::IIrregular_volume_subset>());

    if (!irregular_volume_subset.is_valid_interface()) 
    {
        ERROR_LOG << "VTK IRV Importer: Cannot create an irregular volume subset.";
        return NULL;
    }
    
    std::ostringstream cache_file_name_str;
    std::ifstream      ifile_cache;

    if (m_use_cache) 
    {
        cache_file_name_str << m_filename << "_"
                            << bounding_box.min.x << "." << bounding_box.max.x << "x"
                            << bounding_box.min.y << "." << bounding_box.max.y << "x"
                            << bounding_box.min.z << "." << bounding_box.max.z
                            << ".ivc";
                            
        ifile_cache.open(cache_file_name_str.str().c_str(), std::ios_base::in | std::ios_base::binary);
    }
    
    // try read irregular volume subset from cache file
    if (m_use_cache && ifile_cache) 
    {
        if(!load_irregular_volume_from_cache(ifile_cache, irregular_volume_subset.get()))
            return NULL;
    }
    else 
    {
        const IVOL_ts_mesh* ivol_mesh = 0;

        {
            mi::base::Lock::Block block(&m_cached_input_lock);

            if (m_cached_input_file != m_filename) 
            {
                std::ifstream ifile(m_filename.c_str(), std::ios_base::in);

                if (!ifile.is_open()) 
                {
                    ERROR_LOG << "VTK IRV Importer: Error opening input file '" << m_filename << "'.";
                    return NULL;
                }
                
                std::string line;
                
                // read header
                getline(ifile, line);
                INFO_LOG << "VTK IRV Importer| Header: " << line;
                
                // read title
                getline(ifile, line);
                INFO_LOG << "VTK IRV Importer| Title: " << line;
                
                // read data type
                getline(ifile, line);
                INFO_LOG << "VTK IRV Importer| Data type: " << "'" << line << "'";
                if(line.compare(0, 5, "ASCII") != 0)
                {
                    ERROR_LOG << "VTK IRV Importer: Datatype unsupported, please use ASCII file";
                    return NULL;
                }

                // Dataset description
                getline(ifile, line);
                INFO_LOG << "VTK IRV Importer| Dataset description: " << line;
                if(line.compare(0, 7, "DATASET") != 0)
                {
                    ERROR_LOG << "VTK IRV Importer: Can't find dataset description";
                    return NULL;
                }
                else
                {
                    std::vector<std::string> line_split = string_split(line);
                    if(line_split.size() != 2 || line_split[1].compare("UNSTRUCTURED_GRID") != 0)
                    {
                        ERROR_LOG << "VTK IRV Importer: DataSet unsupported, please use UNSTRUCTURED_GRID file";
                        return NULL;
                    }
                }
                
                // Points line
                mi::Uint32 nb_global_vertices = 0u;
                mi::Uint32 nb_global_tetrahedrons = 0u;
                
                getline(ifile, line);
                {
                    std::vector<std::string> line_split = string_split(line);
                    if(line_split.size() != 3 || line_split[0].compare("POINTS") != 0)
                    {
                        ERROR_LOG << "VTK IRV Importer: Error in POINTS line";
                        return NULL;
                    }
                    
                    nb_global_vertices = std::strtoul(line_split[1].c_str(), NULL, 0);
                    INFO_LOG << "VTK IRV Importer| Mesh points: " << nb_global_vertices;
                }
                    
                // reset the cached mesh
                m_cached_input_ivol_mesh = IVOL_ts_mesh();

                // read global mesh vertices
                m_cached_input_ivol_mesh.vertices.resize(nb_global_vertices);
                mi::math::Vector<mi::Float32, 3> *cur_vert   = &m_cached_input_ivol_mesh.vertices[0];
                
                mi::math::Bbox<mi::Float32, 3> &bbox = m_cached_input_ivol_mesh.bbox;
                bbox.clear();
                            
                for(mi::Uint32 i=0; i<nb_global_vertices; i++, cur_vert++)
                {
                    getline(ifile, line);
                    std::istringstream iss(line);
                    
                    iss >> cur_vert->x;
                    iss >> cur_vert->y;
                    iss >> cur_vert->z;
                    
                    bbox.insert(*cur_vert);
                }
                
                INFO_LOG << "VTK IRV Importer| Point's Bounding Box: " << m_cached_input_ivol_mesh.bbox << ", " << bbox;

                // read empty line
                getline(ifile, line);
                
                // read cells info
                mi::Uint32 nb_cells = 0u;
                mi::Uint32 nb_to_read = 0u;
                mi::Uint32 points_per_cell = 0u; // 4: tetrahedral, 8: voxel or hexahedral
                bool is_hexa = false; 
                
                // data point data
                mi::Uint32 nb_data_points = 0u;
               
                while(ifile)
                {
                    getline(ifile, line);
                    
                    if(line.compare(0, 5, "CELLS")==0)
                    {
                        {
                            std::vector<std::string> line_split = string_split(line);
                            if(line_split.size() != 3 || line_split[0].compare("CELLS")!=0)
                            {
                                ERROR_LOG << "VTK IRV Importer: Error in CELLS line";
                                return NULL;
                            }
                            
                            nb_cells = std::strtoul(line_split[1].c_str(), NULL, 0);
                            nb_to_read = std::strtoul(line_split[2].c_str(), NULL, 0);
                            points_per_cell = nb_to_read/nb_cells - 1;
                            is_hexa = points_per_cell == 8;

                            nb_global_tetrahedrons = nb_cells;
                            if(is_hexa)
                                nb_global_tetrahedrons *= 5u;
                            
                            INFO_LOG << "VTK IRV Importer| Number of tetrahedrons: " << nb_global_tetrahedrons;
                            
                            // Only VTK_TETRA, VTK_VOXEL or VTK_HEXADRON are supported currently
                            if(points_per_cell != 4 && points_per_cell != 8)
                            {
                                ERROR_LOG << "VTK IRV Importer: CELLS type not suppoerted, please use VTK_TETRA or VTK_VOXEL or VTK_HEXADRON";
                                return NULL;
                               
                            }

                            m_cached_input_ivol_mesh.tetrahedron_indices.resize(nb_global_tetrahedrons);
                            mi::math::Vector<mi::Uint32, 4> *tetrahedrons = &m_cached_input_ivol_mesh.tetrahedron_indices[0];
                            
                            // it's hexahedron or voxel?
                            if(is_hexa)
                            {
                                std::vector<std::vector<mi::Uint32> > hexahedrons(nb_cells);
                                
                                std::vector<std::vector<mi::Uint32> > hexahedron_adjacent_list(nb_cells);
                                std::map<Hexahedron_face, mi::Uint32, Hexahedron_face_compare> face_to_cell_unique;
                                
                                // read hexahedrons and build face to cell map
                                for(mi::Uint32 i=0; i<nb_cells; i++)
                                {
                                    hexahedrons[i].resize(points_per_cell);
                                    
                                    getline(ifile, line);

                                    // parse cell vertex indices
                                    std::istringstream iss(line); 
                                    mi::Uint32 idx = 0u;

                                    mi::Uint32 dummy_value;
                                    iss >> dummy_value; // skip first value (nb of vertices) which is not needed
                                    while (!iss.eof())      
                                        iss >> hexahedrons[i][idx++];
                                    
                                    //update face to cell map
                                    for(mi::Uint32 j=0; j<6; j++)
                                    {
                                        Hexahedron_face face(
                                            hexahedrons[i][FACE_IDS[j][0]],
                                            hexahedrons[i][FACE_IDS[j][1]],
                                            hexahedrons[i][FACE_IDS[j][2]],
                                            hexahedrons[i][FACE_IDS[j][3]]);
                                        
                                        std::pair<std::map<Hexahedron_face, mi::Uint32, Hexahedron_face_compare>::iterator, bool> ret = 
                                            face_to_cell_unique.insert(std::pair<Hexahedron_face, mi::Uint32>(face, i));
                                            
                                        // face is already insert, update hexahedron adjacent list
                                        if(ret.second == false)
                                        {
                                            hexahedron_adjacent_list[ret.first->second].push_back(i);
                                            hexahedron_adjacent_list[i].push_back(ret.first->second);
                                        }
                                        
                                    }
                                }
                                
                                // build tetrahedrons
                                std::vector<mi::Sint8> cells_orientation(nb_cells, 0);
                                
                                if(!set_cell_orientation(cells_orientation, hexahedron_adjacent_list))
                                {
                                    ERROR_LOG << "VTK IRV Importer: Error decomposing hexahedrons";
                                    return NULL;
                                }
                                                                
                                for(mi::Uint32 i=0; i<nb_cells; i++)
                                {
                                    if(cells_orientation[i] > 0) //normall decompostion
                                    {
                                        // tetrahedron(0, 3, 1, 4)
                                        tetrahedrons->x = hexahedrons[i][0];
                                        tetrahedrons->y = hexahedrons[i][3];
                                        tetrahedrons->z = hexahedrons[i][1];
                                        tetrahedrons->w = hexahedrons[i][4];
                                        tetrahedrons++;

                                        // tetrahedron(1, 6, 4, 3) // internal tetrahedron
                                        tetrahedrons->x = hexahedrons[i][1];
                                        tetrahedrons->y = hexahedrons[i][6];
                                        tetrahedrons->z = hexahedrons[i][4];
                                        tetrahedrons->w = hexahedrons[i][3];
                                        tetrahedrons++;
                                        
                                        // tetrahedron(2, 6, 1, 3)
                                        tetrahedrons->x = hexahedrons[i][2];
                                        tetrahedrons->y = hexahedrons[i][6];
                                        tetrahedrons->z = hexahedrons[i][1];
                                        tetrahedrons->w = hexahedrons[i][3];
                                        tetrahedrons++;

                                        // tetrahedron(5, 1, 6, 4)
                                        tetrahedrons->x = hexahedrons[i][5];
                                        tetrahedrons->y = hexahedrons[i][1];
                                        tetrahedrons->z = hexahedrons[i][6];
                                        tetrahedrons->w = hexahedrons[i][4];
                                        tetrahedrons++;
                                        
                                        // tetrahedron(7, 3, 4, 6)
                                        tetrahedrons->x = hexahedrons[i][7];
                                        tetrahedrons->y = hexahedrons[i][3];
                                        tetrahedrons->z = hexahedrons[i][4];
                                        tetrahedrons->w = hexahedrons[i][6];
                                        tetrahedrons++;
                                    }
                                    else // inverted decomposition
                                    {
                                        // tetrahedron(1, 2, 5, 0)
                                        tetrahedrons->x = hexahedrons[i][1];
                                        tetrahedrons->y = hexahedrons[i][2];
                                        tetrahedrons->z = hexahedrons[i][5];
                                        tetrahedrons->w = hexahedrons[i][0];
                                        tetrahedrons++;

                                        // tetrahedron(5, 7, 0, 2) // internal tetrahedron
                                        tetrahedrons->x = hexahedrons[i][5];
                                        tetrahedrons->y = hexahedrons[i][7];
                                        tetrahedrons->z = hexahedrons[i][0];
                                        tetrahedrons->w = hexahedrons[i][2];
                                        tetrahedrons++;
                                        
                                        // tetrahedron(6, 7, 5, 2)
                                        tetrahedrons->x = hexahedrons[i][6];
                                        tetrahedrons->y = hexahedrons[i][7];
                                        tetrahedrons->z = hexahedrons[i][5];
                                        tetrahedrons->w = hexahedrons[i][2];
                                        tetrahedrons++;

                                        // tetrahedron(4, 5, 7, 0)
                                        tetrahedrons->x = hexahedrons[i][4];
                                        tetrahedrons->y = hexahedrons[i][5];
                                        tetrahedrons->z = hexahedrons[i][7];
                                        tetrahedrons->w = hexahedrons[i][0];
                                        tetrahedrons++;
                                        
                                        // tetrahedron(3, 2, 0, 7)
                                        tetrahedrons->x = hexahedrons[i][3];
                                        tetrahedrons->y = hexahedrons[i][2];
                                        tetrahedrons->z = hexahedrons[i][0];
                                        tetrahedrons->w = hexahedrons[i][7];
                                        tetrahedrons++;
                                    }
                                }
                            }
                            // it's tetrahedron
                            else
                            {
                                std::vector<mi::Uint32> cell_vertex_ids(5);

                                for(mi::Uint32 i=0; i<nb_cells; i++)
                                {
                                    getline(ifile, line);

                                    // parse cell vertex indices
                                    std::istringstream iss(line); 
                                    mi::Uint32 idx = 0u;

                                    while (!iss.eof())      
                                        iss >> cell_vertex_ids[idx++];
                                    
                                    tetrahedrons->x = cell_vertex_ids[1];
                                    tetrahedrons->y = cell_vertex_ids[2];
                                    tetrahedrons->z = cell_vertex_ids[3];
                                    tetrahedrons->w = cell_vertex_ids[4];
                                    tetrahedrons++;
                                    
                                }
                            }

                        }
                    }
                    else if (line.compare(0, 9, "CELL_DATA")==0)
                    {
                        // parse cell data
                    }
                    else if (line.compare(0, 10, "POINT_DATA")==0)
                    {
                        // parse point data
                        std::vector<std::string> line_split = string_split(line);
                        if(line_split.size() != 2)
                        {
                            ERROR_LOG << "VTK IRV Importer: Error in POINT_DATA line";
                            return NULL;
                        }

                        // number of data points
                        nb_data_points = std::strtoul(line_split[1].c_str(), NULL, 0);
                        INFO_LOG << "VTK IRV Importer| Number of data points: " << nb_data_points;
                        
                        // last section...reads until the end
                        while(ifile)
                        {
                            getline(ifile, line);
                            //scalar or field ??
                            
                            if(line.compare(0, 7, "SCALARS") == 0)
                            {
                                std::vector<std::string> line_split = string_split(line);
                                if(line_split.empty() || line_split.size() > 4)
                                {
                                    ERROR_LOG << "VTK IRV Importer: Error in SCALARS line";
                                    return NULL;
                                }

                                INFO_LOG << "VTK IRV Importer| SCALAR name: " << line_split[1];
                                
                                // read global attributes
                                m_cached_input_ivol_mesh.attributes.resize(nb_data_points);
                                mi::Float32 *cur_attrib   = &m_cached_input_ivol_mesh.attributes[0];
                                mi::math::Vector<mi::Float32, 2> &attrib_range =  m_cached_input_ivol_mesh.attributes_value_range;
                                
                                for(mi::Uint32 i=0; i<nb_data_points; i++, cur_attrib++)
                                {
                                    // jump line if 'LOOKUP_TABLE'
                                    getline(ifile, line);
                                    if(i==0 && line.compare(0, 12, "LOOKUP_TABLE") == 0)
                                        getline(ifile, line);
                                    
                                    *cur_attrib = static_cast<mi::Float32>(atof(line.c_str()));
                                    
                                    if(i==0)
                                    {
                                        attrib_range.x = attrib_range.y = *cur_attrib;
                                    }
                                    else
                                    {
                                        if(*cur_attrib < attrib_range.x)
                                            attrib_range.x = *cur_attrib;
                                        else if(*cur_attrib > attrib_range.y)
                                            attrib_range.y = *cur_attrib;
                                    }
                                }
                                
                                INFO_LOG << "VTK IRV Importer| Scalar range: " << attrib_range;
                                
                            }
                            else if(line.compare(0, 5, "FIELD") == 0)
                            {
                                ERROR_LOG << "VTK IRV Importer: FIELD data points not supported, please use SCALARS instead.";
                                return NULL;
                            }
                        }
                    }
                    else if (line.compare(0, 10, "CELL_TYPES")==0)
                    {
                        // we can ignore this section for now, since we expect vtk files with homogeneous
                        // cell type of the following types: VTK_TETRA, VTK_VOXEL or VTK_HEXADRON and for those
                        // the correct IndeX internal representation can be deduced by only looking at the
                        // CELLS lines.
                        continue;

                        // This code should handle other cases not mentioned in the previous statement

                        // // parse cell type block
                        // std::vector<std::string> line_split = string_split(line);
                        
                        // if(line_split.empty() || line_split.size() < 2)
                        // {
                            // ERROR_LOG << "VTK IRV Importer: Error in CELL_TYPES line";
                            // return NULL;
                        // }
                        
                        // mi::Uint32 nb_cell_types = static_cast<mi::Uint32>(atoi(line_split[1].c_str()));
                        // if(nb_cell_types < nb_cells)
                        // {
                            // ERROR_LOG << "VTK IRV Importer: Error in CELL_TYPES line";
                            // return NULL;
                        // }
                        
                        
                        // for(mi::Uint32 cell_index=0; cell_index < nb_cells;)
                        // {
                            // getline(ifile, line);
                            // std::istringstream iss(line);    
                            
                            // while (!iss.eof())      
                            // {
                                // mi::Uint32 cell_type;
                                // iss >> cell_type;
                                
                                // cell_index++;
                                
                                // if(cell_index >= nb_cells)
                                    // break;
                                
                            // }
                        // }
                    }
                    
                }
                ifile.close();
                
            }

            m_cached_input_file = m_filename;
            ivol_mesh           = &m_cached_input_ivol_mesh; 
        }
        
        // read tetrahedrons and collect only those that intersects the subset bounding box
        mi::math::Bbox<mi::Float32, 3>      subset_bbox(static_cast<mi::Float32>(bounding_box.min.x),
                                                        static_cast<mi::Float32>(bounding_box.min.y),
                                                        static_cast<mi::Float32>(bounding_box.min.z),
                                                        static_cast<mi::Float32>(bounding_box.max.x),
                                                        static_cast<mi::Float32>(bounding_box.max.y),
                                                        static_cast<mi::Float32>(bounding_box.max.z));

        using mi::math::Vector;
        using mi::math::Vector_struct;
                                                        
        typedef Vector<mi::Uint32, 4>       Vec4ui;
        typedef Vector<mi::Float32, 3>      Vec3f;
                                                        
        mi::Uint32                          nb_subset_vertices     = 0u;
        mi::Uint32                          nb_subset_tetrahedrons = 0u;

        std::vector<Vec4ui>                 subset_tetrahedrons;
        std::vector<Vec3f>                  subset_vertices;
        std::vector<mi::Uint8>              subset_attributes;
        std::map<mi::Uint32, mi::Uint32>    global_to_local_vtx_idx_map; // TODO: change to unordered_map when available
        
        // generate the irregular volume mesh subset (all tetrahedra contained in the subset-bbox)
        // * also generate the global max edge length required by the renderer
        
        const std::vector<Vec4ui>&      global_tet_indices = ivol_mesh->tetrahedron_indices;
        const std::vector<Vec3f>&       global_vertices    = ivol_mesh->vertices;
        const std::vector<mi::Float32>& global_attribs     = ivol_mesh->attributes;

        mi::Float32 max_edge_len_sqr = -(std::numeric_limits<mi::Float32>::max)();
        for (mi::Size t = 0u; t < ivol_mesh->tetrahedron_indices.size(); ++t) 
        {
            const Vec4ui&                   tet_vtx_indices = global_tet_indices[t];
            mi::math::Bbox<mi::Float32, 3>  tet_bbox;

            const Vec3f a = Vec3f(global_vertices[tet_vtx_indices.x]);
            const Vec3f b = Vec3f(global_vertices[tet_vtx_indices.y]);
            const Vec3f c = Vec3f(global_vertices[tet_vtx_indices.z]);
            const Vec3f d = Vec3f(global_vertices[tet_vtx_indices.w]);

            // TODO: use triangle-box intersection for all four faces, much more precise
            tet_bbox.clear();
            tet_bbox.insert(a);
            tet_bbox.insert(b);
            tet_bbox.insert(c);
            tet_bbox.insert(d);
            
            if (tet_bbox.intersects(subset_bbox)) {
                subset_tetrahedrons.push_back(tet_vtx_indices);
            }

            // calculate max edge length
            const Vec3f e0 = b - a;
            const Vec3f e1 = c - a;
            const Vec3f e2 = d - a;
            const Vec3f e3 = c - b;
            const Vec3f e4 = d - b;
            const Vec3f e5 = d - c;

            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e0));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e1));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e2));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e3));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e4));
            max_edge_len_sqr = mi::math::max(max_edge_len_sqr, mi::math::square_length(e5));
        }
        
        // build subset vertex list and remap indices
        nb_subset_tetrahedrons = subset_tetrahedrons.size();

        for (mi::Uint32 t = 0u; t < nb_subset_tetrahedrons; ++t) 
        {
            Vec4ui& tet_vtx_indices = subset_tetrahedrons[t];
            for(mi::Uint32 j = 0; j < 4u; ++j) \
            {
                mi::Uint32& vtx_index = tet_vtx_indices[j];

                std::map<mi::Uint32, mi::Uint32>::const_iterator kt = global_to_local_vtx_idx_map.find(vtx_index);
                if(kt != global_to_local_vtx_idx_map.end()) {
                    vtx_index = kt->second;
                }
                else {
                    const mi::Uint32 new_vtx_idx = subset_vertices.size();
                    global_to_local_vtx_idx_map[vtx_index] = new_vtx_idx;
                    
                    subset_vertices.push_back(global_vertices[vtx_index]);

                    mi::Float32 v_attrib_f = global_attribs[vtx_index];

                    // normalize vertex attributes (disregard negative values here)
                    v_attrib_f =   (v_attrib_f - ivol_mesh->attributes_value_range.x)
                                 / (ivol_mesh->attributes_value_range.y - ivol_mesh->attributes_value_range.x);

                    subset_attributes.push_back(static_cast<mi::Uint8>(v_attrib_f * 255.0f));

                    vtx_index = new_vtx_idx;
                }
            }
        }
        
        nb_subset_vertices                 = subset_vertices.size();
        mi::Uint32 nb_cells                = nb_subset_tetrahedrons;
        mi::Uint32 nb_cell_face_indices    = nb_cells * 4u;

        mi::Uint32 nb_faces                = nb_cell_face_indices;
        mi::Uint32 nb_face_vtx_indices     = nb_faces * 3u;

        nv::index::IIrregular_volume_subset::Mesh_parameters mesh_params;

        // general mesh geometry and topology info
        mesh_params.nb_vertices            = nb_subset_vertices;
        mesh_params.nb_face_vtx_indices    = nb_face_vtx_indices;
        mesh_params.nb_faces               = nb_faces;
        mesh_params.nb_cell_face_indices   = nb_cell_face_indices;
        mesh_params.nb_cells               = nb_cells;
        
        INFO_LOG << "VTK IRV Importer| Vertices   : " << nb_subset_vertices;
        INFO_LOG << "VTK IRV Importer| Cells      : " << nb_cells;
        
        // required mesh information for the renderer
        mesh_params.global_max_edge_length = mi::math::sqrt(max_edge_len_sqr);

        nv::index::IIrregular_volume_subset::Mesh_storage mesh_storage;
        if (!irregular_volume_subset->generate_mesh_storage(mesh_params, mesh_storage)) {
            ERROR_LOG << "VTK IRV Importer| Irregular volume importer (.vtk file): "
                      << "unable to generate irregular volume mesh storage, "
                      << "file: '" << m_filename << "'";
            return NULL;
        }
            
        nv::index::IIrregular_volume_subset::Attribute_parameters attrib_params;
        attrib_params.affiliation        = nv::index::IIrregular_volume_subset::ATTRIB_AFFIL_PER_VERTEX;
        attrib_params.type               = nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT8;
        attrib_params.nb_attrib_values   = nb_subset_vertices;

        nv::index::IIrregular_volume_subset::Attribute_storage attribute_storage;
        if (!irregular_volume_subset->generate_attribute_storage(0u, attrib_params, attribute_storage)) {
            ERROR_LOG << "VTK IRV Importer| Irregular volume importer (.vtk file): "
                      << "unable to generate irregular volume attribute storage, "
                      << "file: '" << m_filename << "'";
            return NULL;
        }
                    
        mi::Uint8* subset_attrib_values = reinterpret_cast<mi::Uint8*>(attribute_storage.attrib_values);
        
        // copy vertices and attributes
        for (mi::Uint32 v = 0u; v < nb_subset_vertices; ++v) {
            mesh_storage.vertices[v] = subset_vertices[v];
            subset_attrib_values[v]  = subset_attributes[v];
        }
        
        mi::Uint32 next_vidx = 0u;
        mi::Uint32 next_fidx = 0u;

        // generate cells, cell's face indices, faces, face's vertex indices 
        for (mi::Uint32 t = 0u; t < nb_subset_tetrahedrons; ++t) {
            const mi::Uint32 a   = subset_tetrahedrons[t].x;
            const mi::Uint32 b   = subset_tetrahedrons[t].y;
            const mi::Uint32 c   = subset_tetrahedrons[t].z;
            const mi::Uint32 d   = subset_tetrahedrons[t].w;

            const Vec3f& av      = subset_vertices[a];
            const Vec3f& bv      = subset_vertices[b];
            const Vec3f& cv      = subset_vertices[c];
            const Vec3f& dv      = subset_vertices[d];

            const Vec3f centroid = (av + bv + cv + dv) * 0.25f;
            
            // TODO: when available, use C++11 lambda
            // * adds a face to the ivol mesh storage
            // * tries to orient faces to have correct CCW vertex ordering
#define IVOL_ADD_TET_FACE(i0, i1, i2, p0, p1, p2)                               \
            do {                                                                \
                const Vec3f e1 = ((p1) - (p0));                                 \
                const Vec3f e2 = ((p2) - (p0));                                 \
                                                                                \
                /* face plane */                                                \
                Vec3f       n = cross(e1, e2);                                  \
                n.normalize();                                                  \
                const mi::Float32 dst = -(dot(n, (p0)));                        \
                                                                                \
                /* tetrahedron centroid distance to face plane */               \
                const mi::Float32 cd = dot(n, centroid) + dst;                  \
                                                                                \
                if (cd < 0.0f) {                                                \
                    /* correct ordering */                                      \
                    mesh_storage.face_vtx_indices[next_vidx + 0u] = (i0);       \
                    mesh_storage.face_vtx_indices[next_vidx + 1u] = (i1);       \
                    mesh_storage.face_vtx_indices[next_vidx + 2u] = (i2);       \
                }                                                               \
                else {                                                          \
                    /* invert ordering */                                       \
                    mesh_storage.face_vtx_indices[next_vidx + 0u] = (i0);       \
                    mesh_storage.face_vtx_indices[next_vidx + 1u] = (i2);       \
                    mesh_storage.face_vtx_indices[next_vidx + 2u] = (i1);       \
                }                                                               \
                                                                                \
                mesh_storage.faces[next_fidx].nb_vertices        = 3u;          \
                mesh_storage.faces[next_fidx].start_vertex_index = next_vidx;   \
                                                                                \
                next_fidx += 1;                                                 \
                next_vidx += 3;                                                 \
            } while (0)

            // cell's face indices
            // * this might seem redundant here, but consider datasets with shared faces among cells
            mesh_storage.cell_face_indices[next_fidx + 0] = next_fidx + 0;
            mesh_storage.cell_face_indices[next_fidx + 1] = next_fidx + 1;
            mesh_storage.cell_face_indices[next_fidx + 2] = next_fidx + 2;
            mesh_storage.cell_face_indices[next_fidx + 3] = next_fidx + 3;

            // fill cell data
            mesh_storage.cells[t].nb_faces         = 4u;
            mesh_storage.cells[t].start_face_index = next_fidx;

            // needs to happen down here as it increases next_fidx and next_vidx
            IVOL_ADD_TET_FACE(a, c, d, av, cv, dv);
            IVOL_ADD_TET_FACE(a, b, c, av, bv, cv);
            IVOL_ADD_TET_FACE(a, d, b, av, dv, bv);
            IVOL_ADD_TET_FACE(b, d, c, bv, dv, cv);

#undef IVOL_ADD_TET_FACE

            // TODO: determine max edge length
            // * needs to be a global measure
            //{
            //    ivol_mesh->m_max_edge_length = max(length(bv - av), ivol_mesh->m_max_edge_length);
            //    ivol_mesh->m_max_edge_length = max(length(cv - av), ivol_mesh->m_max_edge_length);
            //    ivol_mesh->m_max_edge_length = max(length(dv - av), ivol_mesh->m_max_edge_length);

            //    ivol_mesh->m_max_edge_length = max(length(cv - bv), ivol_mesh->m_max_edge_length);
            //    ivol_mesh->m_max_edge_length = max(length(dv - bv), ivol_mesh->m_max_edge_length);

            //    ivol_mesh->m_max_edge_length = max(length(dv - cv), ivol_mesh->m_max_edge_length);
            //}
        }
            
        // mi::Size subregion_storage = 
            // sizeof(mesh_storage.vertices[0])*mesh_params.nb_vertices +
            // sizeof(mesh_storage.face_vtx_indices[0])*mesh_params.nb_face_vtx_indices +
            // sizeof(mesh_storage.faces[0])*mesh_params.nb_faces +
            // sizeof(mesh_storage.cell_face_indices[0])*mesh_params.nb_cell_face_indices +
            // sizeof(mesh_storage.cells[0])*mesh_params.nb_cells;
            
            
        // INFO_LOG << "VTK IRV Importer| Subregion storage: " << subregion_storage << " [bytes]";
        
        // save irregular volume subset to cache file
        if (m_use_cache) 
        {
            std::ofstream ofile_cache(cache_file_name_str.str().c_str(), std::ios_base::out | std::ios_base::binary);
            if (ofile_cache) 
                save_irregular_volume_to_cache(ofile_cache, mesh_storage, irregular_volume_subset.get());
        }
        
    
    }
                    
    irregular_volume_subset->retain();
    return irregular_volume_subset.get();
}

//----------------------------------------------------------------------
const char* Irregular_volume_importer_vtk::get_configuration() const
{
    return m_configuration.c_str();
}

//----------------------------------------------------------------------
void Irregular_volume_importer_vtk::serialize(
    mi::neuraylib::ISerializer * serializer) const
{
    serializer->write(&m_use_cache, 1);
    
    mi::Uint32 nb_elements = mi::Uint32(m_filename.size());
    serializer->write(&nb_elements, 1);
    serializer->write(reinterpret_cast<const mi::Uint8*>(m_filename.c_str()), nb_elements);
}

//----------------------------------------------------------------------
void Irregular_volume_importer_vtk::deserialize(
    mi::neuraylib::IDeserializer * deserializer)
{
    deserializer->read(&m_use_cache, 1);
    
    mi::Uint32 nb_elements = 0;
    deserializer->read(&nb_elements, 1);
    m_filename.resize(nb_elements);
    deserializer->read(reinterpret_cast<mi::Uint8*>(&m_filename[0]), nb_elements);
}

//----------------------------------------------------------------------
bool Irregular_volume_importer_vtk::load_irregular_volume_from_cache(
    std::ifstream&                                  ifile,
    nv::index::IIrregular_volume_subset*            irregular_volume_subset) const
{
    // read mesh storage params
    nv::index::IIrregular_volume_subset::Mesh_parameters mesh_params;
    ifile.read(reinterpret_cast<char*>(&mesh_params), 
        sizeof(nv::index::IIrregular_volume_subset::Mesh_parameters));
            
    // create and mesh storage arrays
    nv::index::IIrregular_volume_subset::Mesh_storage mesh_storage;
    if(!irregular_volume_subset->generate_mesh_storage(mesh_params, mesh_storage))
        return false;
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.vertices),
        sizeof(mi::math::Vector_struct<mi::Float32, 3>)*mesh_params.nb_vertices);
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.face_vtx_indices),
        sizeof(mi::Uint32)*mesh_params.nb_face_vtx_indices);
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.faces),
        sizeof(nv::index::IIrregular_volume_subset::Face)*mesh_params.nb_faces);
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.cell_face_indices),
        sizeof(mi::Uint32)*mesh_params.nb_cell_face_indices);
        
    ifile.read(reinterpret_cast<char*>(mesh_storage.cells),
        sizeof(nv::index::IIrregular_volume_subset::Cell)*mesh_params.nb_cells);
    
    // read attribute parameters
    mi::Uint32 nb_attributes = 0u;
    ifile.read(reinterpret_cast<char*>(&nb_attributes), sizeof(mi::Uint32));

    for(mi::Uint32 i=0; i<nb_attributes; ++i)
    {
        mi::Uint32 attribute_index = 0u;
        ifile.read(reinterpret_cast<char*>(&attribute_index), sizeof(mi::Uint32));
        
        nv::index::IIrregular_volume_subset::Attribute_parameters attrib_params;
        ifile.read(reinterpret_cast<char*>(&attrib_params), 
            sizeof(nv::index::IIrregular_volume_subset::Attribute_parameters));

        // create and read attribute storage
        nv::index::IIrregular_volume_subset::Attribute_storage attribute_storage;
        if(!irregular_volume_subset->generate_attribute_storage(attribute_index, attrib_params, attribute_storage))
            return false;
            
        switch(attrib_params.type)
        {
            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT8:
            
                ifile.read(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Uint8)*attrib_params.nb_attrib_values);
                break;

            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT16:

                ifile.read(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Uint16)*attrib_params.nb_attrib_values);
                break;
                
            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_FLOAT32:
            
                ifile.read(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Float32)*attrib_params.nb_attrib_values);
                break;
        }
    }
    
    return true;
}

//----------------------------------------------------------------------
bool Irregular_volume_importer_vtk::save_irregular_volume_to_cache(
    std::ofstream                                           &ofile,
    const nv::index::IIrregular_volume_subset::Mesh_storage &mesh_storage,
    const nv::index::IIrregular_volume_subset*              irregular_volume_subset) const
{
    nv::index::IIrregular_volume_subset::Mesh_parameters mesh_params = 
        irregular_volume_subset->get_mesh_parameters();

    // write mesh storage parameters
    ofile.write(reinterpret_cast<char*>(&mesh_params), 
        sizeof(nv::index::IIrregular_volume_subset::Mesh_parameters));
        
    // write mesh storage
    ofile.write(reinterpret_cast<char*>(mesh_storage.vertices),
        sizeof(mi::math::Vector_struct<mi::Float32, 3>)*mesh_params.nb_vertices);
        
    ofile.write(reinterpret_cast<char*>(mesh_storage.face_vtx_indices),
        sizeof(mi::Uint32)*mesh_params.nb_face_vtx_indices);
        
    ofile.write(reinterpret_cast<char*>(mesh_storage.faces),
        sizeof(nv::index::IIrregular_volume_subset::Face)*mesh_params.nb_faces);
        
    ofile.write(reinterpret_cast<char*>(mesh_storage.cell_face_indices),
        sizeof(mi::Uint32)*mesh_params.nb_cell_face_indices);
        
    ofile.write(reinterpret_cast<char*>(mesh_storage.cells),
        sizeof(nv::index::IIrregular_volume_subset::Cell)*mesh_params.nb_cells);
        
    // write attribute storage
    mi::Uint32 nb_attributes = irregular_volume_subset->get_nb_attributes();
    for(mi::Uint32 i=0; i<nb_attributes; ++i)
    {
        // attribute index
        ofile.write(reinterpret_cast<char*>(&i), sizeof(mi::Uint32));
        
        //attribute params
        nv::index::IIrregular_volume_subset::Attribute_parameters attrib_params;
        irregular_volume_subset->get_attribute_parameters(i, attrib_params);
        ofile.write(reinterpret_cast<char*>(&attrib_params), 
            sizeof(nv::index::IIrregular_volume_subset::Attribute_parameters));
            
        nv::index::IIrregular_volume_subset::Attribute_storage attribute_storage;
        irregular_volume_subset->get_attribute(i, attribute_storage);
        
        switch(attrib_params.type)
        {
            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT8:
            
                ofile.write(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Uint8)*attrib_params.nb_attrib_values);
                break;
                
            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_UINT16:
            
                ofile.write(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Uint16)*attrib_params.nb_attrib_values);
                break;

            case nv::index::IIrregular_volume_subset::ATTRIB_TYPE_FLOAT32:
            
                ofile.write(reinterpret_cast<char*>(attribute_storage.attrib_values),
                    sizeof(mi::Float32)*attrib_params.nb_attrib_values);
                break;
        }
    }
    return true;
}

//----------------------------------------------------------------------
}} // namespace nv::index_common
 
