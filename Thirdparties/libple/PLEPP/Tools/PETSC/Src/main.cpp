#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <blitz/array.h>
#include "mmio.hpp"
#include "petsc_v1.hpp"
   

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE *f;
    int M, N, nz;   
    int i, *I, *J;
    double *val;

    if(argc < 2)
    {
      fprintf(stderr, "Usage: %s [martix-market-filename]\n", argv[0]);
      exit(1);
    }
    else    
    { 
      if((f = fopen(argv[1], "r")) == NULL) 
      exit(1);
    }

    if(mm_read_banner(f, &matcode) != 0)
    {
      printf("Could not process Matrix Market banner.\n");
      exit(1);
    }

    if( mm_is_complex(matcode)&&mm_is_matrix(matcode)&&mm_is_sparse(matcode) )
    {
      printf("Sorry, this application does not support ");
      printf("Market Market type: [%s]\n", mm_typecode_to_str(matcode));
      exit(1);
    }

    if ((ret_code=mm_read_mtx_crd_size(f, &M, &N, &nz))!=0) exit(1);

//    nnz = (int *) malloc(M * sizeof(double));
    I   = (int *) malloc(nz * sizeof(int));
    J   = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));


    blitz::Array<int, 1> nnz(M); 
    for (i=0; i<nz; i++)
    {
      fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
      I[i]--;  
      J[i]--;
      
      nnz( I[i] ) += 1; 
/*      
      std::cout<<   I[i] <<" "; 
      std::cout<<   J[i] <<" "; 
      std::cout<< val[i] <<" "; 
      std::cout<<" \n"; 
**/
    }
    if (f !=stdin) fclose(f);

blitz::Array<double,1>  RHS( M );
blitz::Array<double,1>   Xi( M );
RHS = 1; 
Xi  = 0; 

double*  xiv =  Xi.data();
double* rhsv = RHS.data();


//-------------------------------------------------------------| SOLVER |---//   
Temple::Solver<double> S;
S.init(argc, argv);
S.init_matrix(M, N, nz, nnz.data());

S.set_matrix(I, J, val);
S.set_rhs(rhsv);

S.solver();
S.get_x(xiv);
S.end(); 

std::cout<< Xi; 









  std::cout<<"ok\n";
  return 0;
}

/*
*/




