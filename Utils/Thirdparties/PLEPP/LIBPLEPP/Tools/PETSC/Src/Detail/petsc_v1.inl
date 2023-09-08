

TEMPLATE
Solver<T>::Solver()
{
  ii=0; 
  std::cout<<"\n +Solver ";
}   


TEMPLATE
Solver<T>::~Solver()
{
  std::cout<<"\n +Solver end \n";
}   


TEMPLATE
void Solver<T>::init(int argc, char* argv[])
{
  PetscInitialize( &argc, &argv, 0, 0 );
  PetscPrintf( PETSC_COMM_WORLD, "init \n" );
}   


TEMPLATE
void Solver<T>::end()
{
  KSPDestroy( &__sles );
  
  MatDestroy( &__A ); 
  
  VecDestroy( &__b ); 
  VecDestroy( &__x );
    
  PetscFinalize();
}   


TEMPLATE
void Solver<T>::init_matrix(int m, int n, int nz, int* nnz)
{
/*
  if( m!=nz || n!=nz ) 
  {
    std::cout<<" |\n |_ERROR: \'set_matrix\' m != nz || n !=nz  !!\n\n";  
    exit(1); 
  }
*/
  std::cout<<" |_"<< ii++ <<") init_matrix \n"; 
  
  MatCreateSeqAIJ(PETSC_COMM_SELF, m, n, 0, nnz, &__A);
  MatSetOption(__A, MAT_IGNORE_ZERO_ENTRIES, PETSC_TRUE);

  VecCreate( PETSC_COMM_WORLD, &__x );
  VecSetSizes( __x, PETSC_DECIDE, m);
  VecSetFromOptions( __x );

  VecCreate( PETSC_COMM_WORLD, &__b );
  VecSetSizes( __b, PETSC_DECIDE, m);
  VecSetFromOptions( __b );

  __nz = nz; 
  __m  = m; 
}


TEMPLATE
void Solver<T>::set_matrix(int* idxm, int* idxn, T* vals)
{
/*
  MatSetValues(__A, __nz, idxm, __nz, idxn, vals, INSERT_VALUES); 
*/

  for(int i=0; i<__nz; i++)
  {
    int    m = idxm[i];  
    int    n = idxn[i];  
    double v = vals[i]; 

    MatSetValue(__A, m, n, v, INSERT_VALUES); 
  }

  MatAssemblyBegin(__A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(  __A, MAT_FINAL_ASSEMBLY);

  //MatView(__A, PETSC_VIEWER_STDOUT_WORLD);

  std::cout<<"  | |_set_matrix \n"; 
}


TEMPLATE
void Solver<T>::set_rhs(T* vals)
{
  __set_vector(vals, __b); 

  std::cout<<"  | |_set_rhs \n"; 
}


TEMPLATE
void Solver<T>::solver()
{
  int its; 
  
  KSPCreate(PETSC_COMM_WORLD, &__sles); 
  KSPSetOperators(__sles, __A, __A, DIFFERENT_NONZERO_PATTERN);
  KSPSetFromOptions(__sles);

  KSPSolve(__sles, __b, __x);
  KSPGetIterationNumber(__sles, &its);

  PetscPrintf(PETSC_COMM_WORLD, "Solution in %d iterations is:\n", its); 
  
  //VecView(__x, PETSC_VIEWER_STDOUT_WORLD);

  std::cout<<"  | |_set_solver \n"; 
} 



TEMPLATE
void Solver<T>::get_x(T* vals)
{
/*
  double* aux[1]; 
  aux[0] = vals; 
  VecGetArray(__x, aux);
*/
  int* ids = new int[__m]; 
  for(int i=0; i<__m; i++) ids[i] = i; 

  VecGetValues(__x, __m, ids, vals); 

  delete[] ids; 
}


//-----------------------------------------------------------------| AUXs |---//

TEMPLATE
void Solver<T>::__set_vector(T* vals, Vec V)
{
  for(int i=0; i<__m; i++) VecSetValue(V, i, vals[i], INSERT_VALUES);

  VecAssemblyBegin(V);
  VecAssemblyEnd(V);

//  VecView(V, PETSC_VIEWER_STDOUT_WORLD);
}

