#include <vector>
#include "registry.h"

REGISTRY_SIM::REGISTRY_SIM()
{
  boxlist.resize(0);
}

REGISTRY_SIM::~REGISTRY_SIM()
{
    
}

int REGISTRY_SIM::bbox_exist(float minx,float miny,float minz,float maxx,float maxy,float maxz)
{
  for(int i=0;i<boxlist.size();i++)
    {
      if( boxlist[i].minx == minx && boxlist[i].maxx == maxx && boxlist[i].miny == miny && boxlist[i].maxy == maxy && boxlist[i].minz == minz && boxlist[i].maxz == maxz)
	{
	  return 1;
	}
    }
  return -1;
}

int REGISTRY_SIM::add_bbox(float minx,float miny,float minz,float maxx,float maxy,float maxz)
{
  BBOX_SIM newbox;
  newbox.minx = minx;
  newbox.miny = miny;
  newbox.minz = minz;
  newbox.maxx = maxx;
  newbox.maxy = maxy;
  newbox.maxz = maxz;
  newbox.boxid = boxlist.size();
  boxlist.push_back(newbox);
  return newbox.boxid;
}

std::map<unsigned int,unsigned int>::const_iterator REGISTRY_SIM::localfind(int boxid,unsigned int vertex)
{
  return boxlist[boxid].global_to_local_vtx_idx_map.find(vertex);
}

std::map<unsigned int,unsigned int>::const_iterator REGISTRY_SIM::localend(int boxid)
{
  return boxlist[boxid].global_to_local_vtx_idx_map.end();
}

void REGISTRY_SIM::localset(int boxid,unsigned int key,unsigned int value)
{
  boxlist[boxid].global_to_local_vtx_idx_map[key] = value;
  return;
}
void REGISTRY_SIM::globalset(int boxid,unsigned int key,unsigned int value)
{
  boxlist[boxid].local_to_global_vtx_idx_map[key] = value;
  return;
}

std::map<unsigned int,unsigned int>::const_iterator REGISTRY_SIM::globalfind(int boxid,unsigned int vertex)
{
  return boxlist[boxid].local_to_global_vtx_idx_map.find(vertex);
}

std::map<unsigned int,unsigned int>::const_iterator REGISTRY_SIM::globalend(int boxid)
{
  return boxlist[boxid].local_to_global_vtx_idx_map.end();
}
