#ifndef STRING_F2C_H
#define STRING_F2C_H 

#include <vector>
#include <string>
#include <cstring>
#include <stdlib.h> 
#include <iostream>
#include <algorithm>

class string_f2c  
{
private:
  std::vector<std::string>              argvs;
public:

  std::string                      argv_found;

  void
  set_argvs( std::string dummy )
  {
    argvs.push_back( dummy );
  }


  void
  get_argvs( std::string& dummy )
  {
    dummy = argv_found;
  }


  void
  analyse_argvs( std::string token )
  {
    argv_found = "";

    std::string s;

    int argc = argvs.size();
//std::cout<<" argc:'"<< argc <<"' \n";

    int numarg = -1;
    while(++numarg < argc)
    {
      s = argvs[numarg];
//std::cout<<" s:'"<< s <<"' \n";
      if( strncmp(s.c_str(), token.c_str(), token.size()) == 0 )
      {
        if(++numarg < argc)
        {
          s   = argvs[numarg];
//std::cout<<" s:'"<< s <<"' \n";

          int len = s.size();
          char* log_name = (char *)malloc(sizeof(char)*(len+1));
          strncpy(log_name, s.c_str(), len);
          log_name[len] = '\0';
          argv_found = s;
        }
      }
    }

//std::cout<<" argv_found:'"<< argv_found <<"' \n";
  }


  void
  split_string( std::string& cname ) 
  {
    const std::string     to_split( cname );
    const std::vector<std::string>  words = __split_string__(to_split, "/");
    if(words.size()>1)    cname = words[words.size()-1];
  }


//-----------------------------------------------------------------------||---//
  const char*
  f2c(const char* fchar, int fchar_length)
  {
    int i = fchar_length-1;
    while( ((fchar[i]==' ')||(fchar[i]==0)) && (i>=0) ) i--;
    i+=1;

    std::string  cstr;
    for(int j=0; j<i; j++) cstr.push_back(fchar[j]);

    return cstr.c_str();
  }


  std::vector<std::string>
  __split_string__(const std::string& s, const std::string& delim, const bool keep_empty=true)
  {
    std::vector<std::string> result;
    if(delim.empty())
    {
        result.push_back(s);
        return result;
    }
    std::string::const_iterator substart = s.begin(), subend;

    while(true)
    {
        subend = search(substart, s.end(), delim.begin(), delim.end());
        std::string temp(substart, subend);
        if(keep_empty || !temp.empty())
        {
            result.push_back(temp);
        }

        if(subend == s.end())
        {
            break;
        }

        substart = subend + delim.size();
    }
    return result;
  }


  void
  __strncmp__(const char* x, const char* y, int* size, int* ok)
  {
    ok[0] = -1;
    if(strncmp(x, y, size[0]) == 0) ok[0] = 1;
  }


  void
  __error__(bool ok, std::string message)
  {
    if(ok)
    {
      std::cout<<"\n\nERROR: "<< message  <<"!! \n\n";
      exit(1);
    }
  }


}; // class 

#endif
//-----------------------------------------------------------------------||---//


 
static string_f2c* __ptr_class = NULL;

extern "C"
{
  void
  mui_create_()
  {
    __ptr_class = new string_f2c(); 
  }

  void
  mui_delete_()
  {
    delete __ptr_class;
  }


  void
  mui_set_argvs_(const char* ftype, int* ntype)
  {
    __ptr_class->set_argvs( __ptr_class->f2c(ftype, ntype[0]) ); 
  }


  void 
  mui_get_argvs_(char* ftype)
  {
    if(__ptr_class->argv_found.size()>0) strcpy(ftype, __ptr_class->argv_found.c_str());
  }


  void 
  mui_analyse_argvs_(const char* ftype, int* ntype)
  {
    __ptr_class->analyse_argvs( __ptr_class->f2c(ftype, ntype[0]) );
  }
 


}


//-----------------------------------------------------------------------||---//
