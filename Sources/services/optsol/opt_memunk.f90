
subroutine opt_memunk()

  !-----------------------------------------------------------------------
  !    
  ! This routine allocates memory for all the adjoint unknowns of the problem
  ! When using Parall, Master allocates minimum memory
  !
  !----------------------------------------------------------------------
  use def_parame
  use def_master 
  use def_domain
  use def_solver
  use def_optsol
  use mod_memory
  implicit none
  integer(4) :: istat

  ! Adjoint Algebraic solver arrays
  call memory_alloca(memma,'DAMATR','opt_memunk',damatr,nzmat)
  call memory_alloca(memma,'DRHSID','opt_memunk',drhsid,nzrhs)
  call memory_alloca(memma,'AUNKNO','opt_memunk',aunkno,nzrhs)
  
  ! Adjoint Algebraic complex solver arrays
  call memory_alloca(memma,'DAMATX','opt_memunk',damatx,nzmax)
  call memory_alloca(memma,'DRHSIX','opt_memunk',drhsix,nzrhx)
  call memory_alloca(memma,'AUNKNX','opt_memunk',aunknx,nzrhx)

  ! Right-hand side of adjoint problem
  call memory_alloca(memma,'DCOSTR','opt_memunk',dcostr,nzrhs)
  call memory_alloca(memma,'DCOSTX','opt_memunk',dcostx,nzrhx)

  ! Reduced gradient related arrays
  call memory_alloca(memma,'DIFFJ','opt_memunk',diffj,kfl_ndvars_opt)
  call memory_alloca(memma,'DIFFJ_PREV','opt_memunk',diffj_prev,kfl_ndvars_opt)
  call memory_alloca(memma,'DIFFJ_SHOT','opt_memunk',diffj_shot,kfl_ndvars_opt)
  call memory_alloca(memma,'DIFFJ_ISINSIDE','opt_memunk',diffj_isInside,kfl_ndvars_opt)
  call memory_alloca(memma,'DIFFJ_ILLUM','opt_memunk',diffj_illum,kfl_ndvars_opt)

  ! Descent direction related arrays
  call memory_alloca(memma,'DESCDIR','opt_memunk',descdir,kfl_ndvars_opt)
  call memory_alloca(memma,'DESCDIR_PREV','opt_memunk',descdir_prev,kfl_ndvars_opt)

  ! Design variables related arrays
  call memory_alloca(memma,'DESIGN_VARS','opt_memunk',design_vars,kfl_ndvars_opt)
  call memory_alloca(memma,'DESIGN_VARS_REF','opt_memunk',design_vars_ref,kfl_ndvars_opt)
  call memory_alloca(memma,'DESIGN_VARS_PREV','opt_memunk',design_vars_prev,kfl_ndvars_opt)
  call memory_alloca(memma,'DESIGN_VARS_TMP','opt_memunk',design_vars_tmp,kfl_ndvars_opt)

  ! 
  call memory_alloca(memma,'MODULO_DCOST','opt_memunk',modulo_dcost,kfl_ndvars_opt)
  
end subroutine opt_memunk
