subroutine hlm_matrix()

  !----------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_matrix.f90
  ! NAME 
  !    hlm_matrix
  ! DESCRIPTION
  !    This routine constructs the system matrix and the system right-hand-side.
  !    Also, this routine imposes Dirichlet boundary conditions on the system equations.
  ! USES
  !    hlm_dirbcs
  !    hlm_elmope
  ! USED BY
  !    hlm_solite
  !----------------------------------------------------------------------------

  use def_parame
  use def_master
  use def_helmoz
  use def_domain
  use def_kermod, only       :  kfl_ndvars_opt
  use mod_communications
  implicit none

  real(rp)    :: cpu_syma1,cpu_syma2,cpu_syma
  real(rp)    :: cpu_bcs1,cpu_bcs2,cpu_bcs

  call cputim(cpu_syma1)

  if(kfl_servi(ID_OPTSOL)==1) then
   
     if(kfl_first_opt==1)then
        if (INOTMASTER) then
           call hlm_elmopefirst() 
        end if
        call pararr('MIN',0_ip,kfl_ndvars_opt,design_vars)
        call parari('MAX',0_ip,kfl_ndvars_opt,diffj_isInside)
     end if


     if (INOTMASTER) then
        call hlm_elmopedesig() 
     end if
     


     call pararr('MAX',0_ip,kfl_ndvars_opt,diffj_illum)
  
  else

     if (INOTMASTER) then
        call hlm_elmope() 
     end if

  end if


  call cputim(cpu_syma2)
  call PAR_BARRIER() 
  !if (INOTSLAVE) then
  !   write(*,*) 'Anisotropy level:', aniso_hlm
  !   cpu_syma = cpu_syma2 - cpu_syma1
  !   write(*,*) 'Time to construct the system matrix:',cpu_syma
  !endif
     
  call cputim(cpu_bcs1)
  if (INOTMASTER) then
     call hlm_dirbcs()
  end if
  call PAR_BARRIER() 
  call cputim(cpu_bcs2)
  !if (INOTSLAVE) then
  !   cpu_bcs = cpu_bcs2 - cpu_bcs1
  !   write(*,*) 'Time to impose boundary conditions:',cpu_bcs
  !endif

end subroutine hlm_matrix
