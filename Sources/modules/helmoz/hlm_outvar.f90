subroutine hlm_outvar(ivari)

  !-----------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_outvar.f90
  ! NAME 
  !    hlm_outvar
  ! DESCRIPTION
  !    This routine outputs postprocess variables.
  !    Also, master gets values of FE-computed potentials that 
  !    it needs for computation of field vectors. 
  ! USES
  !    postpr
  !    memgen
  !    hlm_postpr
  ! USED BY
  !    hlm_output
  !-----------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain
  use def_helmoz
  use def_kermod
  use mod_outvar,         only : outvar

  implicit none

  integer(ip), intent(in) :: ivari

  !integer(ip)             :: iar3p,ibopo,ipoin
  integer(ip)             :: iascatmp,iavectmp,iar3ptmp
  integer(ip)             :: ibopo,ipoin,jttim,iedge
  integer(ip)             :: jpoin
  real(rp)                :: dummr,rutim,normgradinv 
  character(5)            :: wopos(3)

  iascatmp = 0_ip
  iavectmp = 0_ip
  iar3ptmp = 0_ip
  !ibopo = 0_ip
  rutim = cutim


  !Define postprocess variables
  select case (ivari)  

  case(0_ip)

     return

  case(1_ip)     

     !kfl_cntrl = 1_ip
     !Secondary magnetic vector potential  
     if (INOTMASTER) then
        call memgen(6_ip,ndime,nsite_hlm*nmlsi_hlm)
        gevex => smgvp_hlm
     endif

  case(2_ip)       

     !kfl_cntrl = 1_ip
     !Secondary electric scalar potential  
     if (INOTMASTER) then       
        call memgen(6_ip,nsite_hlm*nmlsi_hlm,0_ip)
        gescx => selsp_hlm
     endif

  case(5_ip)

     if(INOTMASTER) then
        call memgen(0_ip,npoin,0_ip)
        do ipoin=1,npoin
           diffj_hlm(ipoin)=diffj(lninv_loc(ipoin))
           !design_hlm(ipoin)=design_vars(lninv_loc(ipoin))
        end do
        gesca => diffj_hlm
     end if

  case(6_ip)

     if(INOTMASTER) then
        call memgen(0_ip,npoin,0_ip)
        do ipoin=1,npoin
           !diffj_hlm(ipoin)=diffj_scale(lninv_loc(ipoin))
           design_hlm(ipoin)=design_vars(lninv_loc(ipoin))
           !design_hlm(ipoin)=design_vars_tmp(lninv_loc(ipoin))
        end do
        gesca => design_hlm
     end if

  case(7_ip)

     normgradinv = 1.0_rp / sqrt(dot_product(diffj,diffj))
     diffj = diffj * normgradinv
     if(INOTMASTER) then
        do ipoin=1,npoin
           diffj(lninv_loc(ipoin)) = diffj(lninv_loc(ipoin)) * diffj_illum(lninv_loc(ipoin))
           !      !diffj_hlm(ipoin)=diffj_scale(lninv_loc(ipoin))
           !      !diffj_hlm(ipoin)=(diffj(lninv_loc(ipoin)) / diffj_illum(lninv_loc(ipoin))) * (1.0_rp/(exp(design_vars(lninv_loc(ipoin)))))
           !      !diffj_hlm(ipoin)=exp( 1.0e+10_rp * diffj(lninv_loc(ipoin)) / diffj_illum(lninv_loc(ipoin)) )
           !      !diffj_hlm(ipoin)=exp(-diffj(lninv_loc(ipoin))) * exp(design_vars(lninv_loc(ipoin)))
        end do
     end if
     if(INOTMASTER) then
        call memgen(0_ip,npoin,0_ip)
        do ipoin=1,npoin
           diffj_hlm(ipoin)=exp(-diffj(lninv_loc(ipoin))) * exp(design_vars(lninv_loc(ipoin)))
        end do

        gesca => diffj_hlm
     end if

     !      if(INOTMASTER) then
     !         do ipoin=1,npoin
     !            diffj_hlm(ipoin)= diffj(lninv_loc(ipoin))
     !         end do
     !
     !
     !         gesca => diffj_hlm
     !      end if

  case(8_ip)

     if(INOTMASTER) then
        call memgen(0_ip,npoin,0_ip)
        do ipoin=1,npoin
           diffj_hlm(ipoin)=descdir(lninv_loc(ipoin))
           !design_hlm(ipoin)=design_vars(lninv_loc(ipoin))
        end do

        gesca => diffj_hlm
     end if

  case ( 9_ip )
     !
     ! Boundary condition
     !
     if( INOTMASTER ) then
        call memgen(0_ip,npoin,0_ip)
        do iedge = 1,meshe(ndivi) % nedge
           ipoin = lninv_loc(meshe(ndivi) % edge_to_node(1,iedge))
           jpoin = lninv_loc(meshe(ndivi) % edge_to_node(2,iedge))
if(ipoin==3095.and.jpoin==2956) print*,'IEDGE 1=',iedge,kfl_fixno_hlm(1,iedge)
if(jpoin==3095.and.ipoin==2956) print*,'IEDGE 2=',iedge,kfl_fixno_hlm(1,iedge),kfl_coded(1:2,iedge)
           gesca(ipoin) = max(gesca(ipoin),real(kfl_fixno_hlm(1,iedge),rp))
           gesca(jpoin) = max(gesca(jpoin),real(kfl_fixno_hlm(1,iedge),rp))
        end do
     end if
     
  end select
  !
  ! Output postprocess variables
  !
  call outvar(&
       ivari,&
       ittim,cutim,postp(1) % wopos(:,ivari))

  !ittim,rutim,postp(1) % wopos(:,ivari))

  !call outvax(&
  !     ivari,ittim,rutim,wopos(:))
  !     !ivari,ittim,rutim,postp(1) % wopos(1,ivari))

  !Master gets values of FE-computed potentials that it needs for computation of field vectors
  !if (kfl_cntrl /= 0_ip .and. IMASTER) call hlm_postpr(ivari)

  !kfl_cntrl = 0_ip
  !cntrl_hlm = 1_ip

  !Deallocate postprocess variables
  select case (ivari)  
  case(0_ip)
     return
  case(1_ip)     
     !Secondary magnetic vector potential  
     if (INOTMASTER) then
        call memgen(7_ip,ndime,nsite_hlm*nmlsi_hlm)
     endif

  case(2_ip)       
     !kfl_cntrl = 1_ip
     !Secondary electric scalar potential  
     if (INOTMASTER) then       
        call memgen(7_ip,nsite_hlm*nmlsi_hlm,0_ip)
     endif

  endselect

end subroutine hlm_outvar
