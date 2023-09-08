subroutine tem_outvar(ivari)
  !------------------------------------------------------------------------
  !****f* Temper/tem_output
  ! NAME 
  !    tem_output
  ! DESCRIPTION
  !    Output a postprocess variable
  ! USES
  !    postpr
  !    memgen
  ! USED BY
  !    tem_output
  !***
  !------------------------------------------------------------------------
  use def_parame
  use def_master
  use def_domain
  use def_temper
  use mod_gradie
  use mod_lodi_tem
  use mod_communications,  only : PAR_INTERFACE_NODE_EXCHANGE
  use mod_finite_volume,   only : finite_volume_element_to_nodes
  use mod_ADR,             only : ADR_manufactured_nodal_error
  use mod_projec,          only : projec_elements_to_nodes
  use mod_tem_commdom,     only : commdom_alya_cht_node_flux
  use mod_commdom_driver,  only : commdom_driver_get_total_flux
  use mod_commdom_driver,  only : commdom_driver_get_residual 
  use mod_commdom_driver,  only : commdom_driver_n_fixno
#ifdef COMMDOM
  use mod_commdom_dynamic, only: commdom_dynamic_outvar
#endif
 
  implicit none
  integer(ip), intent(in) :: ivari
  integer(ip)             :: ipoin,iline,jpoin,icont,idime
  real(rp)                :: dummr,rutim
  real(rp), pointer       :: aux(:) 
  real(rp), pointer       :: tempeAux(:) 
  real(rp), pointer       :: grad_tempe(:,:) 

  rutim = cutim

  select case (ivari)  

  case(0_ip)
     !
     ! Do nothing
     !
     return

  case(1_ip)
     !
     ! Temperature
     !
     if( kfl_discr_tem == NODAL_SCHEME ) then
        gesca => tempe(:,1) 
     else        
        if( INOTMASTER ) then
          call memgen(zero,npoin,zero)
         !  call finite_volume_element_to_nodes(meshe(ndivi),1_ip,tempe,gesca)
          tempeAux => tempe(:,1)
          call projec_elements_to_nodes(tempeAux,gesca)
        end if
     end if

  case(2_ip)
     !
     ! Heat flux
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call tem_outhfl()
     end if

  case(3_ip)
     !
     ! Tesgs
     !
     !if( INOTMASTER ) ger3p => ADR_tem % sgs 
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call projec_elements_to_nodes(ADR_tem % sgs,gesca)
     end if

  case(4_ip)
     !
     ! Error w/r manufactured solution
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        call ADR_manufactured_nodal_error(ADR_tem,cutim,tempe,gesca)
     end if

  case(5_ip)
     !
     ! Average temperature
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  avtem_tem(ipoin)/dummr
              avtem_tem(ipoin)=0.0_rp
           end do
        end if
     end if

  case(6_ip)
     !
     ! Velocity
     !
     if( kfl_advec_tem /= 0 ) then
        if( INOTMASTER ) then        
           if( kfl_advec_tem == 1 ) then
              gevec => veloc(:,:,1)
           else if( kfl_advec_tem > 1 .and. kfl_paral /= 0 ) then
              call memgen(zero,ndime,npoin)
              call tem_velfun(npoin,coord,gevec)
           end if
        end if
     else
        return
     end if

  case(7_ip)
     !
     ! Turbulent viscosity
     !
     if( INOTMASTER ) then
        if(size(turmu,1)>1) then
           gesca => turmu
        end if
     end if

  case(8_ip)
     !
     ! Residual
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = therm(ipoin,1)-teold_tem(ipoin)
        end do
     end if
     rutim = real(ittot_tem,rp)

  case(9_ip)
     !
     ! GROUPS FOR DEFLATED CG
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,0_ip)
        do ipoin=1,npoin
           gesca(ipoin)=real(solve(1)%lgrou(ipoin),rp)
        end do
     end if

  case(10_ip)
     !
     ! Projection
     !
     gesca => ADR_tem % proje1

  case(11_ip)
     !
     ! Limiter
     !
     if( kfl_limit_tem /= 0 ) then
       call runend('TEM_OUTVAR: LIMITER NOT CODED') 
     end if

  case(12_ip)
     !
     ! LINTE: Linelets of preconditioner 
     !
     if( INOTMASTER ) then
        icont=0
        do ipoin=1,npoin
           rhsid(ipoin)=0.0_rp
        end do
        do iline=1,solve(1)%nline
           icont=icont+1
           do ipoin=solve(1)%lline(iline),solve(1)%lline(iline+1)-1
              jpoin=solve(1)%lrenup(ipoin)
              rhsid(jpoin)=real(icont,rp)
           end do
        end do
        gesca => rhsid
     end if

  case(13_ip)
     !
     ! TEST
     !
     if( INOTMASTER ) then
        call memgen(0_ip,ndime,npoin)
        call memgen(0_ip,npoin,0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = 2.0_rp*coord(1,ipoin)+3.0_rp*coord(2,ipoin)
        end do
        call grasca(gesca,gevec)
        call memgen(2_ip,npoin,0_ip)
     end if

  case(14_ip)
     !
     ! WAT_VAPOR: Water vapor
     !
     if( INOTMASTER ) then
        call memgen(0_ip,npoin,0_ip)
        call tem_poswat()
     end if

  case(15_ip)
     !
     ! KFL_FIXNO_TEM
     !
     if( INOTMASTER ) then
        call memgen(zero,npoin,zero)
        do ipoin = 1,npoin
           gesca(ipoin) = real(kfl_fixno_tem(1,ipoin),rp)
        end do
     end if

  case(18_ip)
     !
     ! Average tempe*tempe
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  avte2_tem(ipoin)/dummr
              avte2_tem(ipoin)=0.0_rp
           end do
        end if
     end if

  case(19_ip)
     !
     ! Average veloc*tempe
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime=1, ndime
                 gevec(idime,ipoin) =  avtev_tem(idime,ipoin)/dummr
                 avtev_tem(idime,ipoin)=0.0_rp
              end do
           end do
        end if
     end if
  case(20_ip)
     !
     ! Average temperature
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  avden_tem(ipoin)/dummr
              avden_tem(ipoin)=0.0_rp
           end do
        end if
     end if
  case(21_ip)
     !
     ! Favre average velocity rho*veloc
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,ndime,npoin)
           do ipoin = 1,npoin
              do idime=1, ndime
                 gevec(idime,ipoin) =  fvvel_tem(idime,ipoin)/dummr
                 fvvel_tem(idime,ipoin)=0.0_rp
              end do
           end do
        end if
     end if

  case(22_ip)
     !
     ! Heat flux computed from matrix RHS
     ! 
     if( INOTMASTER ) then
        call memgen(zero,npoin,0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = solve_sol(1) % reaction(1,ipoin)
        end do
        !call PAR_INTERFACE_NODE_EXCHANGE(gesca,'SUM','IN MY CODE')
     end if 

  case(23_ip)
     !
     ! grad(T)
     !
     if( INOTMASTER ) then
        call memgen(zero,ndime,npoin)
        if( kfl_discr_tem == NODAL_SCHEME ) then
           call grasca(tempe,gevec)
        else
           allocate(grad_tempe(ndime,nelem) )
           !select case( my_grad ) 
           !case ( fv_grad_method_gauss )
              call tem_gradient_gauss(tempe,grad_tempe)
           !case (fv_grad_method_ls)
           !   call tem_gradient_ls(tempe,grad_tempe)
           !case default
            !  call tem_gradient_gauss(tempe,grad_tempe) ! is the one more robust
           !end select
           call projec_elements_to_nodes(grad_tempe,gevec)
           deallocate(grad_tempe)
        end if
     end if

  case(24_ip)
     if( INOTMASTER ) then
        !call memgen(zero,npoin,zero)
        !call ker_proper('DENSI','NPOIN',dummi,dummi,gesca)
        !gesca(1:npoin) = prthe(1)/(gasco*tempe(1:npoin,1))
        !CHARAs%detem => gesca 

       !CHARAs%gamme = gasco
        call lodi_tem_allocate( CHARAs )
        call lodi_tem_get_characteristics( CHARAs )
        !
        call memgen(zero,ndime,npoin) 
        do idime = 1,ndime
          do ipoin = 1,npoin
            gevec(idime,ipoin) = CHARAs%chrc(CHARAs%idofn,ipoin,idime) ! chrc_tem(ndofn_tem,npoin,ndime) 
          enddo
        enddo
        !
        call lodi_tem_deallocate( CHARAs )
     end if

  case(25_ip)
     !
     ! Heat flux computed from matrix RHS
     ! 
     if( INOTMASTER ) then
        call memgen(zero,npoin,0_ip)
        do ipoin = 1,npoin
           gesca(ipoin) = solve_sol(1) % bvnat(1,ipoin) 
        end do
        call PAR_INTERFACE_NODE_EXCHANGE(gesca,'SUM','IN MY CODE')
     end if 

  case(26_ip)
     !
     ! Ethalpy
     !
     gesca => therm(:,1)


  case(27_ip)
     !
     ! 'REACT'
     !
     if( INOTMASTER ) then
       if(.not.associated(solve_sol(1)%lpoin_reaction) )  call runend('ERROR: -->POSTPROCESS REACT<-- ')
       call memgen(zero,npoin,0_ip)
       gesca(1:npoin) = -1.0_rp
       where( solve_sol(1) % lpoin_reaction(1:npoin) ) gesca(1:npoin) = 1.0_rp
     endif 

  case(28_ip)
     !
     ! 'TFLUX' Heat flux interpolated at the nodes computed from Fourier law
     !
     if( INOTMASTER ) then
       call memgen(zero,npoin,zero)
       call tem_outhfl()
       nullify(aux)
       allocate( aux(npoin) )   
         aux(1:npoin) = gesca(1:npoin) 
       gesca(1:npoin) = 0.0_rp  
       call commdom_driver_get_total_flux(  aux(1:npoin), gesca(1:npoin) ) 
       deallocate( aux )
     endif


  case(29_ip )
    ! 'RESID' 
    call commdom_driver_get_residual() 

  case(30_ip)
     !
     ! Projection
     !
     if( associated(ADR_tem % proje2) ) then
        gesca => ADR_tem % proje2
     else
        return
     end if

  case(31_ip)
     !
     ! 'RESHE' Heat flux interpolated at the nodes computed from Residuals
     !
     if( INOTMASTER ) then
       call memgen(zero,npoin,zero)
       gesca(1:npoin) =  solve(1)%reaction(1,1:npoin)
     endif

  case(32_ip)
     !
     ! Average heat flux interpolated at the nodes computed from Residuals
     !
     if( rutim > avtim_tem ) then
        dummr = rutim - avtim_tem
        if( INOTMASTER ) then
           call memgen(zero,npoin,zero)
           do ipoin = 1,npoin
              gesca(ipoin) =  avres_tem(ipoin)/dummr
              avres_tem(ipoin)=0.0_rp
           end do
        end if
     end if


  case(33_ip) !< 2016Feb24   
     ! NFIXN 
     if( INOTMASTER ) then
       call memgen(0_ip,npoin,0_ip)
       gesca(1:npoin) = 0.0_rp
       call commdom_driver_n_fixno( gesca ) 
     endif 

  case(34_ip) !< 2016Feb24   
     ! TOUCH 
     if( INOTMASTER ) then
#ifdef COMMDOM
       call commdom_dynamic_outvar( )
#else 
       call memgen(0_ip,npoin,0_ip)
       gesca(1:npoin) = -69.69_rp
#endif
     endif 


  end select


  call outvar(&
       ivari,&
       ittim,rutim,postp(1) % wopos(1,ivari))

end subroutine tem_outvar
