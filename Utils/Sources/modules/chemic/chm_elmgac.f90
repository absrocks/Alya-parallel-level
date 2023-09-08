subroutine chm_elmgac(&
     itask,ielem,pnode,lnods,elden,elcod,elcon,elco2,elvel,eltem,elDik,elmas,&
     elmol,elmut,elkey,eleps,elrrt,elsen)
  !------------------------------------------------------------------------
  !****f* partis/chm_elmgac
  ! NAME 
  !    chm_elmgat
  ! DESCRIPTION
  !    Gather operations for the combustion models
  ! USES
  ! USED BY
  !    chm_elmope
  !***
  !------------------------------------------------------------------------ 
  use def_kintyp, only     :  ip,rp 
  use def_domain, only     :  ndime,mnode,npoin,coord
  use def_master, only     :  conce,tempe,wmean,speci,visck,turmu,untur,ID_CHEMIC,ID_NASTAL,&
                              kfl_coupl,kfl_paral,condk,sphek,lescl,massk
  use def_chemic, only     :  kfl_timei_chm,kfl_advec_chm,&
       &                      kfl_tiacc_chm,kfl_tisch_chm,&
       &                      nclas_chm,nspec_chm,iclas_chm,&
       &                      proje_chm,kfl_stabi_chm,kfl_advec_chm,&
       &                      veloc_chm,vterm_chm,tmrat_chm,densi_chm,&
       &                      lawvt_chm,kfl_sourc_chm,lawte_chm, &
       &                      kfl_stagg_chm, &
       &                      denma_chm,temma_chm,sourc_chm,lawdi_chm,&
       &                      diffu_chm,kfl_diffu_chm,kfl_tfles_chm,tfles_chm,&
       &                      kfl_model_chm,kfl_cotur_chm,flsen_chm, &
       &                      ADR_chm
  use def_kermod, only     :  gasco
  use mod_ADR,    only     :  BDF
  use mod_ker_proper

  implicit none 
  integer(ip), intent(in)  :: itask,pnode,ielem
  integer(ip), intent(in)  :: lnods(pnode)
  real(rp),    intent(in)  :: elden(pnode)
  real(rp),    intent(out) :: elcod(ndime,pnode)
  real(rp),    intent(out) :: elcon(pnode,nspec_chm,*)
  real(rp),    intent(out) :: elmut(pnode)                          ! Turbulent viscosity
  real(rp),    intent(out) :: elrrt(pnode)                          ! {w_c * c} for CFI model
  real(rp),    intent(out) :: elco2(pnode,nspec_chm)
  real(rp),    intent(out) :: elvel(ndime,pnode)
  real(rp),    intent(out) :: eltem(pnode)
  real(rp),    intent(out) :: elDik(pnode,nspec_chm)                ! Species diffusion coefficient
  real(rp),    intent(out) :: elmas(pnode,nspec_chm)                ! Mass source terms
  real(rp),    intent(out) :: elmol(pnode)                          ! Average molar mass 
  real(rp),    intent(out) :: elkey(pnode),eleps(pnode)
  real(rp),    intent(out) :: elsen(pnode)
  integer(ip)              :: inode,ipoin,idime,itime,iclas,ispec,dummi
  real(rp)                 :: elvis(pnode),dummr(3,3)
  real(rp)                 :: elhco(pnode),elsph(pnode)
  real(rp), pointer        :: prope_tmp(:)


  if (itask >= 1_ip) then  ! Everybody calls this, elcon(1), elcod, elmas and eldik
     !
     ! Current concentration and coordinates
     !
     do inode=1,pnode
        ipoin=lnods(inode)
        do iclas=1,nspec_chm
           elcon(inode,iclas,1) = conce(ipoin,iclas,1)
        end do
        do idime=1,ndime
           elcod(idime,inode)   = coord(idime,ipoin)
        end do
     end do

     !
     ! Turbulent viscosity if needed
     !
     if (kfl_tfles_chm >= 1_ip) then
        elsen = 0.0_rp
        do inode=1,pnode
           ipoin=lnods(inode)
           elsen(inode) = flsen_chm(ipoin)
        end do
        
        if ( kfl_coupl(ID_CHEMIC,ID_NASTAL) > 0_ip ) then    !! bring the turbulent viscosity to chemic
           elmut = 0.0_rp
           do inode=1,pnode
              ipoin=lnods(inode)
              elmut(inode) = turmu(ipoin)
           end do
        end if        
     end if

     !
     ! TKE and EPS for turbulence in CFI model
     !
     if (kfl_model_chm == 5_ip) then
        elkey = 0.0_rp
        eleps = 0.0_rp
        elrrt = 0.0_rp          ! Transport of reaction rate by turbulent fluctuations {w_c * c}
        do inode=1,pnode
           ipoin=lnods(inode)
           elrrt(inode) = lescl(ipoin)
        end do
        if (kfl_cotur_chm == 1_ip ) then  !! K-EPSILON MODEL
           elmut = 0.0_rp
           do inode=1,pnode
              ipoin=lnods(inode)
              elkey(inode) = untur(1,ipoin,1)
              eleps(inode) = untur(2,ipoin,1)
              elmut(inode) = turmu(ipoin)
           end do
        elseif (kfl_cotur_chm == 2_ip ) then   !! K-OMEGA MODEL
           elmut = 0.0_rp
           do inode=1,pnode
              ipoin=lnods(inode)
              elkey(inode) = untur(1,ipoin,1)
              eleps(inode) = 0.09_rp * untur(1,ipoin,1) * untur(2,ipoin,1)
              elmut(inode) = turmu(ipoin)
           end do
        endif       

     end if

     !
     ! Mass source terms coefficients
     !
     if (kfl_stagg_chm==0 .or. itask==1_ip) then ! Complete step with sources
        select case (kfl_sourc_chm)
        case (0_ip)
           do ispec = 1, nspec_chm
              do inode = 1,pnode
                 elmas(inode,ispec) = 0.0_rp
              enddo
           end do
        case (1_ip)
           do ispec = 1, nspec_chm
              do inode = 1,pnode
                 ipoin = lnods(inode)
                 elmas(inode,ispec) = sourc_chm
              enddo
           end do
        case (4_ip)
           do ispec = 1, nspec_chm
              do inode = 1,pnode
                 ipoin = lnods(inode)
                 elmas(inode,ispec) = massk(ipoin,ispec)
              enddo
           end do
        case default
           call runend('CHEMIC: Wrong source option for combustion')
        end select
     else !Staggered step, mass sources are solved later
        elmas=0.0_rp
     endif
     !
     ! Diffusion coefficients
     !
     do ispec = 1, nspec_chm
        select case (lawdi_chm(1,ispec))

        case (0_ip)
           ! 
           ! Constant diffusion
           !
           do inode = 1,pnode
              ipoin = lnods(inode)
              elDik(inode,ispec) = diffu_chm(1,ispec)
           enddo

        case (1_ip)
           ! 
           ! Each species has its own diffusion coeff
           !
           do inode = 1,pnode
              ipoin = lnods(inode)
              elDik(inode,ispec) = visck(ipoin,ispec)/(elden(inode)*speci(ispec)%lewis*speci(ispec)%prand) 
           enddo

        case (2_ip) 
           !
           ! Molecular model
           !
           call runend("CHEMIC: Molecular diffusion coefficient not coded")

        case (3_ip) ! Uniform Lewis and Prandtl numbers
           call ker_proper('VISCO','PNODE',dummi,ielem,elvis,pnode,dummi)
           do inode = 1,pnode
              elDik(inode,ispec) = elvis(inode)/(elden(inode)* diffu_chm(1,1)* diffu_chm(2,1)) !!!! FER: Should I use viscosity from previous step or this one????
           enddo

        case (4_ip)
           !
           ! Watervapor
           !
           do inode = 1,pnode
              ipoin = lnods(inode)
              !elDik(inode,ispec) = 21.2e-6_rp*(1.0_rp+0.0071_rp*tempe(ipoin,1))
              elDik(inode,ispec) = 2.338e-5_rp ! old value: 21.2e-6_rp
           end do     

        case (5_ip)
           ! 
           ! Individual diffusion coefficients + effects of turbulence (RANS/LES)
           !
           if (kfl_model_chm == 5) then
              call ker_proper('CONDU','PNODE',dummi,ielem,elhco,pnode,dummi)
              call ker_proper('SPHEA','PNODE',dummi,ielem,elsph,pnode,dummi)
              do inode = 1,pnode    !! for CFI model elDik is lambda/cp and only used for gpgrd
                 elDik(inode,ispec) = elhco(inode) / elsph(inode)
              end do
           else
              nullify ( prope_tmp )
              allocate( prope_tmp(npoin) )
  
              prope_tmp = 0.0_rp  !! Laminar flow

              if (kfl_cotur_chm > 0_ip .or. kfl_coupl(ID_CHEMIC,ID_NASTAl) > 0_ip) then  !! RANS / COMPRESSIBLE
                 prope_tmp = turmu
              else if (kfl_cotur_chm < 0_ip ) then           !! LES
                 call ker_proper('TURBU','NPOIN',dummi,dummi,prope_tmp)
                 do inode=1,pnode
                    ipoin = lnods(inode)
                    prope_tmp(ipoin) = prope_tmp(ipoin)*elden(inode)
                 end do 
              endif

              do inode = 1,pnode
                 ipoin = lnods(inode)
                 elDik(inode,ispec) = visck(ipoin,ispec)/(elden(inode)*speci(ispec)%lewis*speci(ispec)%prand) + & ! Resolved scale contribution
                                      prope_tmp(ipoin) / (diffu_chm(1,1) * elden(inode))                          ! Turbulence contribution
              end do

              deallocate(prope_tmp) 
           end if

        end select
     end do

  endif

  if (itask >= 2) then ! Assembly and critical time call this
     do inode=1,pnode
        ipoin=lnods(inode)
        do iclas=1,nspec_chm
           elco2(inode,iclas)=conce(ipoin,iclas,2)
        end do
     end do
     !
     ! Advection
     !
     if( kfl_advec_chm /= 0 ) then
        do inode = 1,pnode
           ipoin = lnods(inode)
           do idime = 1,ndime
              elvel(idime,inode) = veloc_chm(idime,ipoin)
           end do
        end do
     else
        do inode = 1,pnode
           do idime = 1,ndime
              elvel(idime,inode) = 0.0_rp
           end do
        end do
     end if
  endif

  if (itask >= 3_ip) then ! Only assembly calls this
     !
     ! Time integration
     !
     if( ADR_chm(1) % kfl_time_integration == 1 ) then
        do iclas = 1,nclas_chm
           do inode = 1,pnode
              ipoin = lnods(inode)
              elcon(inode,iclas,2) = conce(ipoin,iclas,3)
           end do
        end do
        if( ADR_chm(1) % kfl_time_scheme == BDF ) then
           do iclas = 1,nclas_chm
              do itime = 3,ADR_chm(iclas) % kfl_time_order + 1
                 do inode = 1,pnode
                    ipoin = lnods(inode)
                    elcon(inode,iclas,itime) = conce(ipoin,iclas,itime+1)
                 end do
              end do
           end do
        end if
     end if
     !
     ! Average molar mass
     !
     if (kfl_diffu_chm==2) then
       do inode = 1,pnode
          ipoin = lnods(inode)
          elmol(inode) = wmean(ipoin,1)
       enddo
     endif
  endif

  if (itask == 1_ip .or. itask == 2_ip ) then ! Heat source term also cares about temperature
     !
     ! Temperature
     !
     select case (lawte_chm)
     case(0)
        do inode = 1,pnode
           ipoin = lnods(inode)
           eltem(inode) = temma_chm(1)
        end do
     case(-2)
        do inode = 1,pnode
           ipoin = lnods(inode)
           eltem(inode) = tempe(ipoin,1)
        end do
     end select
  endif
  


end subroutine chm_elmgac



!!$  !
!!$  ! Projection
!!$  !
!!$  if( kfl_stabi_chm >= 1 ) then
!!$    do inode = 1,pnode
!!$       ipoin = lnods(inode)
!!$       elpro(inode) = 0.0_rp !proje_chm(ipoin,iclas_chm)
!!$    end do
!!$  end if




