#include <cstdlib>    // exit() 
#include <iostream>   // cout  
#include "petsc.h"

#ifndef PETSC_HEADER
#define PETSC_HEADER


#define GIGA           1e9 
#define TEMPLATE       template<class T>
#define DEBUG_INDICES  0 

namespace Temple{


TEMPLATE
class Solver
{
  public:   
    Solver();   
   ~Solver(); 
   
    void init(int, char**); 
    void end(); 
    
    void init_matrix(int, int, int, int*); 
    void set_matrix(int*, int*, T*);
    void set_rhs(T*);  
    void solver(); 
    void get_x(T*); 

  protected: 
   
  private:   
    void __set_vector(T*, Vec);    

    KSP  __sles;
    Mat  __A;
    Vec  __b; 
    Vec  __x;
       
    int __nz;
    int __m;  
    int ii; 
};


#include "Detail/petsc_v1.inl"
}

#undef TEMPLATE

#endif
