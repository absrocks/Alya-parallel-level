#ifndef __SIMregistry__
#define __SIMregistry__

#include <map>
#include <vector>

class BBOX_SIM {
 public:
  std::map<unsigned int,unsigned int>    global_to_local_vtx_idx_map;
  std::map<unsigned int,unsigned int> local_to_global_vtx_idx_map;
  int *global_numbering;
  double *pt_to_scalardata;
  float minx,miny,minz;
  float maxx,maxy,maxz;
  int boxid;
};

class REGISTRY_SIM {
 
  std::vector<BBOX_SIM> boxlist;
  
 public:
  int bbox_exist(float minx,float miny,float minz,float maxx,float maxy,float maxz);
  REGISTRY_SIM();
  ~REGISTRY_SIM();
  int getboxcount();
  int add_bbox(float minx,float miny,float minz,float maxx,float maxy,float maxz);
  std::map<unsigned int,unsigned int>::const_iterator localfind(int boxid,unsigned int vertex);
  std::map<unsigned int,unsigned int>::const_iterator localend(int boxid);
  void localset(int boxid,unsigned int key,unsigned int value);
  void globalset(int boxid,unsigned int key,unsigned int value);
  std::map<unsigned int,unsigned int>::const_iterator globalfind(int boxid,unsigned int vertex);
  std::map<unsigned int,unsigned int>::const_iterator globalend(int boxid);

};

#endif
