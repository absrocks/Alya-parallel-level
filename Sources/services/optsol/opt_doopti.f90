subroutine opt_doopti
  !------------------------------------------------------------------------
  !****f* Optsol/opt_doopti
  ! NAME
  !    opt_doopti
  ! DESCRIPTION
  ! OUTPUT
  ! USED BY
  !    Optsol
  !***
  !------------------------------------------------------------------------
  use      def_parame
  use      def_master
  use      def_domain
  use      def_optsol
  use      def_inpout
  use      def_solver
  implicit none

  integer(ip) :: indvars
  character(len=8) :: bufftmp,bufftmp2,bufftmp3,bufftmp4
  real(rp) :: normdiffj,normdiffjprev,normdiffjinv

  bufftmp='(I4.0)'
  bufftmp2='(I4.0)'

  if(kfl_goopt==1) then

     call opt_descdir()

     call opt_upddes()

     !---------OPTIMIZATION REPORT---------------------------------------
     !
     !---------(MUST BE INTEGRATED INTO ALYA's MESSAGE SYSTEM)-----------
     !
     !-------------------------------------------------------------------

     !if(kfl_curlin_opt==1) then
     !   normdiffj=0.0_rp
     !   do indvars=1,kfl_ndvars_opt
     !      normdiffj=normdiffj+diffj(indvars)**2
     !   end do
     !   normdiffj=sqrt(normdiffj)
     !   normdiffjinv = 1.0_rp / normdiffj
     !end if

     write(bufftmp3,bufftmp) kfl_curstp_opt
     write(bufftmp4,bufftmp) kfl_curlin_opt

     if(kfl_paral==-1 .or. kfl_paral==1) then 
         !normdiffj=0.0_rp
         !normdiffjprev=0.0_rp
         !do indvars=1,kfl_ndvars_opt
         !   normdiffj=normdiffj+diffj(indvars)**2
         !   normdiffjprev=normdiffjprev+diffj_prev(indvars)**2
         !end do
         !normdiffj=sqrt(normdiffj)
         !normdiffjprev=sqrt(normdiffjprev)
         write(*,*) 'DBG(costf)'//trim(bufftmp3)//' '//trim(bufftmp4),   costf
         write(*,*) 'DBG(costfp)'//trim(bufftmp3)//' '//trim(bufftmp4),   costf_prev
         write(*,*) 'DBG(normd)'//trim(bufftmp3)//' '//trim(bufftmp4),   sqrt(dot_product(design_vars,design_vars))
         write(*,*) 'DBG(normdt)'//trim(bufftmp3)//' '//trim(bufftmp4),  sqrt(dot_product(design_vars_tmp,design_vars_tmp)) 
         write(*,*) 'DBG(normdp)'//trim(bufftmp3)//' '//trim(bufftmp4),  sqrt(dot_product(design_vars_prev,design_vars_prev)) 
         write(*,*) 'DBG(ngrad)'//trim(bufftmp3)//' '//trim(bufftmp4),   sqrt(dot_product(diffj,diffj))
         write(*,*) 'DBG(ngradp)'//trim(bufftmp3)//' '//trim(bufftmp4),  sqrt(dot_product(diffj_prev,diffj_prev))
         write(*,*) 'DBG(step)'//trim(bufftmp3)//' '//trim(bufftmp4),   stepj
     end if

     modul=0
     call moddef(9_ip)

  end if

  call Parall(20_ip) !MPI_Barrier

end subroutine opt_doopti
