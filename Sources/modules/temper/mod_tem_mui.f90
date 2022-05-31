!-------------------------------------------------------------------------||---!
!
!< 2016AUG25  
!
!-------------------------------------------------------------------------||---!
module mod_tem_mui 
  use def_parame,           only: ip, rp
  use def_master,           only: IMASTER, INOTMASTER, ISEQUEN
  use def_domain,           only: npoin, nboun
  use def_kintyp,           only: soltyp
  use def_master,           only: momod, modul
  use def_domain,           only: coord, mnode, ndime, npoin
  use mod_mui_driver,       only: MUI_SENDRECV
  use mod_mui,              only: MUI_CPLNG
  use mod_mui,              only: mui_exchange  
#ifdef MUI  
#endif
  implicit none
  type(soltyp), pointer :: solve(:)
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!
  private
  !
  public:: tem_mui_plugin 
  !
  !-----------------------------------------------------------------------||---!
  !
  !-----------------------------------------------------------------------||---!

  !=============================================================| contains |===!
contains
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!
!-------------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  subroutine tem_mui_plugin()
    call tem_mui_plugin01()
  end subroutine
  !-----------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !   +
  !   |_Alya                                       
  !     |_call Turnon()                            
  !     |_call Iniunk()                             
  !     |_time: do while
  !       |_call Timste()                         
  !       |_reset: do 
  !         |_call Begste()              TEMPE-->, HEATF<--0  
  !           |_block: do while                          
  !             |_coupling: do while                    
  !               |_call Begzon()        TEMPE-->, HEATF<--0  [sendrecv]
  !               |_modules: do while                           
  !                 |_call Doiter()                
  !                 |_call Concou()                
  !               |_call Endzon()        TEMPE<--, HEATF-->                  
  !                                                                           
  !             |_call Conblk()                             
  !       |_call Endste()                                     
  !   |_call Turnof()                    
  !
  !-----------------------------------------------------------------------||---!
  subroutine tem_mui_plugin01()                                             !  
  use def_temper,          only: kfl_regim_tem, bvess_tem, bvnat_tem, kfl_fixno_tem
  use def_temper,          only: kfl_plepp_tem
  use def_master,          only: therm, title
  !
#ifdef MUI
#endif 
  !
  implicit none
  real(rp)     :: relax_temp = 1.0_rp 
  character(4) :: nameij = '-->'
  character(4) :: nameji = '<--' 
  integer(ip)  :: n_recv
  !
  real(rp), pointer       :: aux(:,:) => null()
  ! 
  !-----------------------------------------------------------------------||---!
  n_recv = 0
  if(INOTMASTER) n_recv = npoin

  allocate( aux(ndime,n_recv) )
  aux = huge(0_ip)
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  if( any(MUI_SENDRECV) ) then
    !
#ifdef MUI  
    !---------------------------------------------------------------------||---!
    !---------------------------------------------------------------------||---!
    if(MUI_CPLNG%current_code==MUI_CPLNG%code_i) then
      if( MUI_SENDRECV(7) ) then ! -ENDZON+
        !-------------------------------------------------------| HEATF--> |---!
        if(inotmaster) then 
          !
!          MUI_CPLNG%var_ij(1,1:npoin) = 0.0
!          call mui_alya_cht_node_flux( MUI_CPLNG%var_ij(1,1:npoin) )
!          MUI_CPLNG%var_ij(1,1:npoin) = -MUI_CPLNG%var_ij(1,1:npoin)
          !
        endif 
        !-----------------------------------------------------------------||---!
        !
        call mui_exchange( MUI_CPLNG, &
                           nameji="YYY", tji=0.0_rp, varji=aux(2_ip,1:n_recv), &
                           nameij="YYY", tij=0.0_rp, varij=aux(1_ip,1:n_recv)  &
                         )
        !
        !-------------------------------------------------------| TEMPE<-- |---!
        if(inotmaster) then 
!          call temp2enth(  MUI_CPLNG%var_ji(1,1:npoin) )
!          call mui_dynamic_set_vals( MUI_CPLNG%var_ji(1,1:npoin),     therm(  1:npoin,1), relax_op=relax_temp ) ! T(n)
        endif
        !-----------------------------------------------------------------||---!
      endif
    endif
    !---------------------------------------------------------------------||---!
    if(MUI_CPLNG%current_code==MUI_CPLNG%code_j) then
      if( MUI_SENDRECV(7) ) then ! +BEGZON-
        !-------------------------------------------------------| TEMPE--> |---!
        if(inotmaster) then 
          if (kfl_regim_tem == 4) call tem_clippi()
!          MUI_CPLNG%var_ij(1,1:npoin) = therm(1:npoin,1)
        endif
        !-----------------------------------------------------------------||---!
        !
        call mui_exchange( MUI_CPLNG, &
                           nameji="YYY", tji=0.0_rp, varji=aux(2_ip,1:n_recv), &
                           nameij="YYY", tij=0.0_rp, varij=aux(1_ip,1:n_recv) &
                         )
        !
        !-------------------------------------------------------| HEATF<-- |---!
        if(inotmaster) then 
!          call mui_alya_cht_nodes2bvnat( MUI_CPLNG%var_ji(1,1:npoin), bvnat_tem(3,1:nboun,1) )
!          MUI_CPLNG%var_ji(1,1:npoin) = 0.0_rp                              ! <--- reset values for the next step...
        endif 
        !-----------------------------------------------------------------||---!
      endif
    endif    
    !---------------------------------------------------------------------||---!
    ! 
!    print *, "[tem_mui_plugin01] ",  "'"//trim(title)//"."//trim(nameji)//"'"
    !
#endif
    !
    deallocate( aux )
    !
  endif
  !-----------------------------------------------------------------------||---!
  !                                                                            !
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



!-------------------------------------------------------------------------||---!
!----------------------------------------------------------------| PRIVATE |---!
!-------------------------------------------------------------------------||---!


  !-----------------------------------------------------------------------||---!
  !> @author Dani mira 
  !> @date    
  !> @brief   
  !> @details 
  !
  !-----------------------------------------------------------------------||---!
  subroutine temp2enth(temperature) !< ok  
  use def_temper,           only: kfl_regim_tem, kfl_fixno_tem
  use def_master,           only: sphec  
  use def_domain,           only: npoin, nboun 
  use mod_physics,          only: physics_T_2_HCp
  implicit none 
  real(rp), intent(out) :: temperature(npoin) 
  real(rp) :: dummr = 0.0_rp, cploc(6,2) = 0.0_rp, hnew=0.0_rp  
  integer(ip) :: ipoin  
 
  if(kfl_regim_tem==4) then
     do ipoin=1,npoin
        if(kfl_fixno_tem(1,ipoin)==1) then
           cploc(1:6,1:2) = sphec(ipoin,1:6,1:2)
           call physics_T_2_HCp(temperature(ipoin),cploc,hnew,dummr)
           temperature(ipoin) = hnew 
        endif 
     enddo 
  end if
  !-----------------------------------------------------------------------||---!
  end subroutine
  !-----------------------------------------------------------------------||---!



  !-----------------------------------------------------------------------||---!
  subroutine mui_alya_cht_node_flux( prop ) !< ok 
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_memchk
  use mod_postpr
  use mod_gradie
  implicit none
  real(rp), intent(inout) :: prop(npoin) 
  integer(ip)             :: ipoin,ibopo,idime
  integer(4)              :: istat
  real(rp), allocatable   :: gradt(:,:)
  !
  ! Allocate memory
  allocate( gradt(ndime,npoin), stat=istat)
  call memchk(zero,istat,mem_modul(1:2,modul),'GRADT','tem_bouflux',gradt)
  !
  ! Compute temperature gradients
  call tem_heatfl( gradt )

  ! Compute heat flux
  do ipoin = 1,npoin
    prop(ipoin) = 0.0_rp 
    ibopo = lpoty(ipoin)
    if(ibopo >= 1) then
        do idime = 1,ndime
          prop(ipoin) = prop(ipoin) + gradt(idime,ipoin) * exnor(idime,1,ibopo)
        end do
    endif 
  end do
  !
  ! Deallocate memory
  !
  call memchk(two,istat,mem_modul(1:2,modul),'GRADT','tem_outhfl',gradt)
  deallocate(gradt,stat=istat)
  if(istat/=0) call memerr(two,'GRADT','tem_outhfl',0_ip)
  !-----------------------------------------------------------------------||---!
  end subroutine mui_alya_cht_node_flux
  !-----------------------------------------------------------------------||---!

  !=============================================================| contains |===!
end module mod_tem_mui 
!==============================================================================!
!==============================================================================!
