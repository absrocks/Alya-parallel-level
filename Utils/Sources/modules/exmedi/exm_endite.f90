subroutine exm_endite(itask)
  !-----------------------------------------------------------------------
  !****f* Exmedi/exm_endite
  ! NAME 
  !    exm_endite
  ! DESCRIPTION
  !    This routine checks convergence and updates unknowns at:
  !    - itask=1 The end of an internal iteration
  !    - itask=2 The end of the internal loop iteration
  ! USES
  !    exm_cvgunk
  ! USED BY
  !    exm_doiter
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_parame
  use      def_exmedi
  use mod_messages, only : livinf


  implicit none
  integer(ip) :: itask, ipoin
  character(300)   :: messa
  integer(ip)      :: maxiter

  select case(itask)

  case(one)
     !
     !  Compute convergence residual of the internal iteration (that is,
     !  || u(n,i,j) - u(n,i,j-1)|| / ||u(n,i,j)||) and update unknowns:
     !  u(n,i,j-1) <-- u(n,i,j)
     !
     call exm_cvgunk(  one)
     call exm_updunk(ITERK_EQ_UNKNO)     !  u(,ITER_K) <-- unkno  

     maxiter= solve_sol(1) % miter
     messa = &
          ' (SUBIT: '//trim(intost(itinn(modul)))//'/'//trim(intost(miinn_exm))//' IT: '//trim(intost(last_iters_exm))//'/'//trim(intost(maxiter))//')'
     call livinf(-3_ip,messa,one)
     call livinf(56_ip,' ',modul)


  case(two)
     !
     !  Compute convergence residual of the external iteration (that is,
     !  || u(n,i,*) - u(n,i-1,*)|| / ||u(n,i,*)||) and update unknowns:
     !  u(n,i-1,*) <-- u(n,i,*)
     !
     call livinf(16_ip,' ',itinn(modul))
     call exm_cvgunk( two)
     call exm_updunk(ITERAUX_EQ_ITERK) !  u(,ITER_AUX) <-- u(,ITER_K)
     call exm_upcell(ITERAUX_EQ_ITERK) !  u(,ITER_AUX) <-- u(,ITER_K)

  end select

end subroutine exm_endite
