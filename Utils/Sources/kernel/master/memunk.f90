
#ifdef NINJA
module pinnedmemorygpu
  integer*8 :: c_amatr,c_rhsid,c_unkno  
end module pinnedmemorygpu
#endif

subroutine memunk(itask)

  !-----------------------------------------------------------------------
  !****f* master/memunk
  ! NAME
  !    Turnon
  ! DESCRIPTION
  !    This routine allocates memory for all the unknowns of the problem
  !    When using Parall, Master allocates minimum memory
  ! USES
  ! USED BY
  !    Alya
  !***
  !-----------------------------------------------------------------------
#ifdef NINJA
  use pinnedmemorygpu
#endif
  
  use def_parame
  use def_master 
  use def_kermod
  use def_domain
  use def_solver
  use mod_memory
#ifdef NINJA
  use iso_c_binding
#endif
  implicit none
  integer(ip), intent(in) :: itask

#ifdef NINJA
  integer*8 :: d_temp,sz
  type(c_ptr) :: cptr
#endif
  !
  ! Deallocate before allocating
  !
  if( itask == 2 ) then
  end if
  !
  ! Algebraic solver
  !

#ifdef NINJA
  sz = nzmat
  sz = sz * 8

!  call gpumallochost(d_temp,sz)
!  call int2cptr(d_temp,cptr)
!  call c_f_pointer(cptr,amatr,[nzmat])
!  c_amatr = d_temp
  call memory_alloca(memma,'AMATR','memunk',amatr,nzmat)

  sz = nzrhs
  sz = sz * 8

!  call gpumallochost(d_temp,sz)
!  call int2cptr(d_temp,cptr)
!  call c_f_pointer(cptr,rhsid,[nzrhs])
!  c_rhsid = d_temp
!  call gpumallochost(d_temp,sz)
!  call int2cptr(d_temp,cptr)
!  call c_f_pointer(cptr,unkno,[nzrhs])
!  c_unkno = d_temp
  call memory_alloca(memma,'RHSID','memunk',rhsid,nzrhs)
  call memory_alloca(memma,'UNKNO','memunk',unkno,nzrhs)

#else
  
  call memory_alloca(memma,'AMATR','memunk',amatr,nzmat)
  call memory_alloca(memma,'RHSID','memunk',rhsid,nzrhs)
  call memory_alloca(memma,'UNKNO','memunk',unkno,nzrhs)

#endif

  call memory_alloca(memma,'BMATR','memunk',bmatr,nzmbt)
  call memory_alloca(memma,'PMATR','memunk',pmatr,nzpre)
  call memory_alloca(memma,'ERRES','memunk',erres,nzerr)
  !
  ! Eigen solver
  !
  call memory_alloca(memma,'EIGEN','memunk',eigen,neige)
  call memory_alloca(memma,'EIGVA','memunk',eigva,neiva)
  !
  ! Algebraic complex solver 
  !
  call memory_alloca(memma,'AMATX','memunk',amatx,nzmax)
  call memory_alloca(memma,'RHSIX','memunk',rhsix,nzrhx)
  call memory_alloca(memma,'UNKNX','memunk',unknx,nzrhx)
  call memory_alloca(memma,'PMATX','memunk',pmatx,nzprx)
  !
  ! lumpped mass matrix - for use in dual time step preconditioner
  !
   if( INOTMASTER ) call memgeo(63_ip) ! not sure if this is the optimal place
   
end subroutine memunk
