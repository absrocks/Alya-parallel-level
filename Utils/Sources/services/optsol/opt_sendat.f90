subroutine opt_sendat(order)
  !-----------------------------------------------------------------------
  !****f* Optsol/opt_sendat
  ! NAME
  !    opt_sendat
  ! DESCRIPTION
  !    This routine exchange data related to the number of design variables
  !    and numerical parameters
  ! USES
  ! USED BY
  !    opt_turnon
  !***
  !-----------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_optsol
  use def_inpout
  use mod_memchk
  use mod_opebcs
  implicit none
  integer(ip), intent(in) :: order
  integer(ip)             :: ji,jr,ki,kr,ir,kfl_ptask_old,dummi
  integer(4)              :: istat
  real(rp)                :: extmp1, extmp2, extmp3

  select case (order)


  case(-1_ip) 

     !------------------------------------------------------------------- 
     !
     ! Exchange kfl_ndvars_opt read in opt_reapro
     ! to be used in master/memunk
     !
     !------------------------------------------------------------------- 

     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1
     do parii=1,2 
        npari=0
        nparr=0
        nparc=0
        call iexcha(kfl_ndvars_opt)
        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','opt_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','opt_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(2_ip)
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','opt_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','opt_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','opt_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','opt_sendat',0_ip)     

  case(1_ip)    

     !------------------------------------------------------------------- 
     !
     ! Exchange other data read in opt_reapro
     !
     !------------------------------------------------------------------- 

     kfl_ptask_old= kfl_ptask
     kfl_ptask    = 1

     do parii=1,2 
        npari=0
        nparr=0
        nparc=0

        extmp1=kfl_conver_opt 
        extmp2=kfl_scale_opt 
        extmp3=kfl_steplength_opt 

        call rexcha(extmp1)
        call rexcha(extmp2)
        call rexcha(extmp3)


        call iexcha(kfl_method_opt)
        call iexcha(kfl_maxstp_opt)
        call iexcha(kfl_maxlin_opt)
        call iexcha(kfl_maxbar_opt)
        call iexcha(kfl_maxres_opt)
        call iexcha(kfl_first_opt)
        call iexcha(kfl_ndvars_opt)
        call iexcha(kfl_scheme_opt)


        !
        ! Allocate memory for the first pass
        !
        if(parii==1) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parin','opt_sendat',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,ID_PARALL),'parre','opt_sendat',parre)
           if(kfl_paral>=1.or.kfl_ptask==2) call Parall(2_ip)
     
        end if
     end do

     if(kfl_paral==0.and.kfl_ptask/=2) call Parall(2_ip)

     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parin','opt_sendat',parin)
     deallocate(parin,stat=istat)
     if(istat/=0) call memerr(two,'parin','opt_sendat',0_ip)
     call memchk(two,istat,mem_servi(1:2,ID_PARALL),'parre','opt_sendat',parre)
     deallocate(parre,stat=istat)
     if(istat/=0) call memerr(two,'parre','opt_sendat',0_ip)     

  end select

  npari=0
  nparr=0
  nparc=0

end subroutine opt_sendat
