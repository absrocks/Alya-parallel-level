subroutine chm_omegak(inisp,endsp)
  !-----------------------------------------------------------------------
  !****f* partis/chm_omegak
  ! NAME 
  !    chm_omegak
  ! DESCRIPTION
  !    Compute mass source of each species
  ! USES
  ! USED BY
  !    chm_begite
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only      :  ip,rp
  use def_domain, only      :  ndime,npoin,vmass,nboun,nnode,ltypb,lnodb
  use def_master, only      :  tempe,conce,inotmaster,wmean,speci,turmu,&
                               ID_CHEMIC,ID_NASTAL,kfl_coupl,kfl_paral,massk
  use mod_communications, only: PAR_MAX
  use def_chemic
  use mod_ker_proper 
  use def_kermod

  implicit none
  integer(ip),intent(in)    :: inisp,endsp
  integer(ip)               :: ireac,ispec,ipoin,icoef
  integer(ip)               :: ifuel,ioxo2,ioxn2           ! Index for fuel and oxidizer 
  integer(ip)               :: pblty,pnodb,inodb,iboun 
  real(rp)                  :: Q(nreac_chm)  ! Progress rate of each reaction
  real(rp)                  :: kforw,kback,xconF,xconB,T,kpres,chape, Fcent, Ntroe,Ctroe,Ftroe
  real(rp)                  :: conce_temp(nspec_chm)
  real(rp)                  :: difnu(nreac_chm,nspec_chm)  ! Stoichiometric coefficients: backward-forward
  real(rp)                  :: fuel_air_ratio              ! Fuel / Air for DTFLES
  real(rp)                  :: flame_speed,flame_thick     ! Flame speed and flame thickness
  real(rp)                  :: sgs_wrink                   ! Subgrid scale wrinkling factor for DTFLES
  real(rp)                  :: tfles_factor                ! Dynamic thickened flame factor
  real(rp)                  :: rate_max(nspec_chm)         ! Maximum reaction rate in the domain
  real(rp)                  :: beta                        ! Model constant controlling flame transition (DTFLES)
  integer(ip)               :: dummi
  real(rp)                  :: dummr(3)
  real(rp), pointer         :: prope_tmp(:)

  !
  ! Init
  !
  if (INOTMASTER) then
     do ispec = 1,nspec_chm
        do ireac = 1,nreac_chm
           difnu(ireac,ispec) = 0.0_rp
        enddo
     enddo
     do ireac = 1,nreac_chm
        do icoef=6, 5+lreac_chm(ireac)%l(1)          !Reactants
           ispec = lreac_chm(ireac)%l(  icoef )
           difnu(ireac,ispec) = difnu(ireac,ispec)-stoic_chm(ispec,ireac,1)
        enddo
        do icoef=6 +lreac_chm(ireac)%l(1),5 +lreac_chm(ireac)%l(1)+lreac_chm(ireac)%l(2) !Products
           ispec = lreac_chm(ireac)%l(icoef )
           difnu(ireac,ispec) = difnu(ireac,ispec)+stoic_chm(ispec,ireac,2)
        enddo
     enddo

     nullify ( prope_tmp )
     allocate( prope_tmp(npoin) )
     
     call ker_proper('DENSI','NPOIN',dummi,dummi,prope_tmp)

     do ipoin = 1,npoin
        select case (lawte_chm)
        case(0)
           T=temma_chm(1)
        case(-2)
           T=tempe_chm(ipoin)
        end select

        if (kfl_temli_chm /= 0 .and. T > temli_chm) T = temli_chm

        conce_temp=conce(ipoin,:,1)

        call chm_heat_of_reactions(T,prope_tmp(ipoin),conce_temp,Q)

        rate_max = 0.0_rp
        do ispec= inisp,endsp
           massk(ipoin,ispec) = 0.0_rp
           do ireac = 1,nreac_chm
              massk(ipoin,ispec) = massk(ipoin,ispec) + speci(ispec)%weigh * difnu(ireac,ispec) * Q(ireac)
              rate_max(ispec) = max(rate_max(ispec),abs(massk(ipoin,ispec)))
           enddo
        enddo        
     enddo

     if (kfl_wallc_chm == 1_ip ) then
       !
       ! Imposing reaction source term to zero at walls
       !
       boundaries: do iboun = 1,nboun
          pblty = ltypb(iboun)
          pnodb = nnode(pblty)
          do inodb = 1,pnodb
             ipoin = lnodb(inodb,iboun)
             massk(ipoin,inisp:endsp) = 0.0_rp
          end do
       end do boundaries
     endif

     if (kfl_tfles_chm >= 1) then
        !
        ! Correction of the reaction rate for the DTFLES
        !
        if (kfl_coupl(ID_CHEMIC,ID_NASTAl) > 0_ip) then
           prope_tmp = turmu
        else
           call ker_proper('TURBU','NPOIN',dummi,dummi,prope_tmp) 
        endif
        ifuel = 1  ! Fuel
        ioxo2 = 2  ! O2
        ioxn2 = 5  ! N2
        beta  = 5.0_rp

        !
        ! Maximum reaction rate in the domain for flame front sensor
        !
        call PAR_MAX(endsp,rate_max,'IN MY CODE')
        flsen_chm = 0.0_rp
 
        do ipoin = 1,npoin

           !
           ! Fuel air ratio 
           !
           fuel_air_ratio = conce(ipoin,ifuel,1) / ( conce(ipoin,ioxo2,1) + conce(ipoin,ioxn2,1))

           !
           ! Flame properties dependent on equivalence ratio
           !
           call chm_flame_prop(fuel_air_ratio,flame_speed,flame_thick) 

           !
           ! Flame sensor OMEGA
           !
           if (massk(ipoin,1_ip) /= 0.0_rp) flsen_chm(ipoin) = tanh(beta*abs(massk(ipoin,1_ip))/rate_max(1_ip))

           do ispec = inisp,endsp

              !
              ! Efficiency function and dynamic thickening factor
              !
              call chm_tfles_sgs(flame_speed,flame_thick,&
                           prope_tmp,vmass(ipoin),sgs_wrink,&
                           tfles_factor,flsen_chm(ipoin))
         
              !
              ! Correction of reaction rate
              !
              massk(ipoin,ispec) = massk(ipoin,ispec) * sgs_wrink / tfles_factor
           enddo
        enddo
     endif
     deallocate( prope_tmp )

  else
     if (kfl_tfles_chm >= 1) call PAR_MAX(endsp,rate_max,'IN MY CODE')  ! Master also has to participate in MPI_ALLREDUCE
  endif
end subroutine chm_omegak

