subroutine hlm_reaphy()

  !------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_reaphy.f90
  ! NAME 
  !    hlm_reaphy
  ! DESCRIPTION
  !    This routine reads the physical problem. 
  ! USES
  ! USED BY
  !    hlm_turnon
  !------------------------------------------------------------------------

  use def_parame
  use def_inpout
  use def_master
  use def_helmoz
  use def_domain
  use mod_ecoute, only :  ecoute

  implicit none

  integer(ip) :: imate,isite,ii,jj,ishot

  if (INOTSLAVE) then
     !Read and write files
     lispa = 0
     lisda = lun_pdata_hlm       !Reading file
     lisre = lun_outpu_hlm       !Writing file

     !Initializations (default values)
     nequs_hlm = 1_ip            !Number of equations (number of unknown parameters: Asx, Asy, Asz, Psis)
     ppout_hlm = 2_ip            !Output secondary fields
     ncond_hlm = 3_ip            !By default we have 3 variables in the conductivity tensor
     frequ_hlm = 1.0_rp          !f = Frequency
     anguf_hlm = 1.0_rp          !w = Angular frequency
     nsite_hlm = 50_ip           !Number of sites
     nmlsi_hlm = 30_ip           !Number of closest mesh nodes to a site needed for the MLSI
     airpl_hlm = 1.0e30_rp       !Off by default, MLSI without limits
     nshot_hlm = 0_ip
     kfl_edges_hlm = 0           ! Edge element formulation

     !Allocate memory for the properties
     call hlm_memphy(1_ip)

     !Reach the section
     call ecoute('hlm_reaphy')
     do while (words(1) /= 'PHYSI')
        call ecoute('hlm_reaphy')
     enddo

     !Begin to read data
     do while (words(1) /= 'ENDPH')
        call ecoute('hlm_reaphy')

        if (words(1) == 'PROBL') then                        !Problem definition data

           call ecoute('hlm_reaphy')

           do while (words(1) /= 'ENDPR')

              if( words(1) == 'EDGEE') then
                 !
                 ! Edge element formulation
                 !
                 if( option('EDGEE') ) then
                    kfl_edges_hlm = 1
                 else
                    kfl_edges_hlm = 0
                 end if
                 
              else if (words(1) == 'SOURC') then      
                 !
                 ! Sources
                 !
                 if(exists('CSEMC')) emmet_hlm = 1_ip      !Current loop as a source
                 if(exists('CSEMH')) emmet_hlm = 2_ip      !CSEM with HED
                 if(exists('CSEMV')) emmet_hlm = 3_ip      !CSEM with VED
                 if(exists('MTELX')) emmet_hlm = 4_ip      !Magnetotellurics (natural-source) method (MT) in one of orthogonal directions
                 if(exists('MTELY')) emmet_hlm = 5_ip      !Magnetotellurics (natural-source) method (MT) in the other orthogonal direction
                 if(exists('INTER')) ppcod_hlm = 1_ip      !Read potential values in the z-r plane and interpolate on nodes
                 if(exists('DIREC')) ppcod_hlm = 2_ip      !Read potential values directly on nodes
                 if(exists('MYSEL')) ppcod_hlm = 3_ip      !Generate them here

              else if (words(1) == 'SHOTS') then
                 nshot_hlm = getint('SHOTS',50_ip,'#Number of shots')
                 if (nshot_hlm /= 0_ip) then 
                    call hlm_memphy(5_ip)     !Allocate memory for the Cartesian coordinates, current and length of the shots
                 end if

              else if (words(1) == 'NUMBE') then           
                 !
                 ! Number of equations 
                 !
                 nequs_hlm = getint('NUMBE',4_ip,'#Number of equations')

              else if (words(1) == 'FREQU') then       
                 !    
                 ! f = Frequency
                 !
                 frequ_hlm = getrea('FREQU',1.0_rp,'#Frequency')
                 anguf_hlm = 2 * frequ_hlm * 3.141592653589793E0

              else if (words(1) == 'ANGUL') then           
                 !    
                 ! w = Angular frequency
                 !    
                 anguf_hlm = getrea('ANGUL',1.0_rp,'#Angular frequency')

              else if (words(1) == 'SITES') then           
                 !    
                 ! Number of sites
                 !    
                 nsite_hlm = getint('SITES',50_ip,'#Number of sites')
                 if (nsite_hlm /= 0_ip) then 
                    call hlm_memphy(3_ip)                  !Allocate memory for the Cartesian coordinates of the sites
                 else 
                    nsite_hlm = npoin                      !Calculate values of field intensities in all mesh nodes
                    site_hlm => coord
                 endif

              else if (words(1) == 'NMLSI') then           
                 !    
                 ! Number of closest mesh nodes to a site needed for the MLSI
                 !
                 nmlsi_hlm = getint('NMLSI',30_ip,'#Number of closest mesh nodes to a site needed for the MLSI')
                 call hlm_memphy(4_ip)                     !Allocate memory for clnod1_hlm, clnod2_hlm and clcoor_hlm

              else if (words(1) == 'AIRPL') then           
                 !    
                 ! Z-coordinate of air-plane that separates MLSI regions
                 !
                 airpl_hlm = getrea('AIRPL',0.0_rp,'#Z-coordinate of airplane')

              else if (words(1) == 'OUTPU') then           !What to output

                 if(exists('PRIMA')) ppout_hlm = 1_ip      !Primary
                 if(exists('SECON')) ppout_hlm = 2_ip      !Secondary
                 if(exists('TOTAL')) ppout_hlm = 3_ip      !Total

              end if
              call ecoute('hlm_reaphy')
           enddo

        else if (words(1) == 'PROPE') then                   !Properties

           call ecoute('hlm_reaphy')

           do while (words(1) /= 'ENDPR')
              if (words(1) == 'MAGNE') then                !mu = Magnetic permeability  
                 call ecoute('hlm_reaphy')
                 do while (words(1) /= 'ENDMA')
                    imate = int(param(1))
                    if (imate < 1 .or. imate > nmate) call runend('HLM_REAPHY: WRONG MAGNETIC PERMEABILITY')
                    perma_hlm(imate) = param(2)
                    call ecoute('hlm_reaphy')
                 enddo
              else if (words(1) == 'DIELE') then           !eps = Dielectric permittivity
                 call ecoute('hlm_reaphy')
                 do while (words(1) /= 'ENDDI')
                    imate = int(param(1))
                    if (imate < 1 .or. imate > nmate) call runend('HLM_REAPHY: WRONG DIELECTRIC PERMITTIVITY')
                    epsil_hlm(imate) = param(2)
                    call ecoute('hlm_reaphy')
                 enddo
              else if (words(1) == 'ELECT') then           !sig = Electric conductivity
                 call ecoute('hlm_reaphy')
                 do while (words(1) /= 'ENDEL')
                    imate = int(param(1))
                    if (imate < 1 .or. imate > nmate) call runend('HLM_REAPHY: WRONG ELECTRIC CONDUCTIVITY')
                    do ii=1,ncond_hlm
                       sigma_hlm(ii,imate) = param(ii+1)
                    enddo
                    call ecoute('hlm_reaphy')
                 enddo
                 aniso_hlm = 0_ip                         !Isotropic by default
                 do ii=1,nmate
                    if (aniso_hlm == 0_ip .and. sigma_hlm(1,ii) /= sigma_hlm(3,ii)) aniso_hlm = 1_ip    ! TIV: Gx = Gy /= Gz
                    if (aniso_hlm <= 1_ip .and. sigma_hlm(1,ii) /= sigma_hlm(2,ii)) aniso_hlm = 2_ip    ! TIH: Gx /= Gy = Gz
                    if (((sigma_hlm(1,ii) /= sigma_hlm(2,ii)) .and. sigma_hlm(2,ii) /= sigma_hlm(3,ii))) aniso_hlm = 3_ip     ! Triaxial: Gx /= Gy /= Gz
                 enddo
              else if (words(1) == 'BACKG') then           !bcksig = Background electric conductivity
                 call ecoute('hlm_reaphy')
                 do while (words(1) /= 'ENDBA')
                    do ii=1,ncond_hlm
                       bckco_hlm(ii) = param(ii)
                    enddo
                    call ecoute('hlm_reaphy')
                 enddo
              endif
              call ecoute('hlm_reaphy')
           enddo
	else if (words(1) == 'SHOTS') then               !Shots
    !print *,'nshot_hlm=',nshot_hlm
           call ecoute('hlm_reaphy')
           do while (words(1) /= 'ENDSH')
              ishot = int(param(1))
              if (ishot < 1 .or. ishot > nshot_hlm) call runend('HLM_REAPHY: WRONG SHOT')

              xoffsv_hlm(ishot)=param(2)  !X-offset of the shot 
              yoffsv_hlm(ishot)=param(3)  !Y-offset of the shot
              zoffsv_hlm(ishot)=param(4)  !Z-offset of the shot
              elcurv_hlm(ishot)=param(5)  !Dipole current
              lengthv_hlm(ishot)=param(6) !Dipole length


              !print *,'xoffsv_hlm(',ishot,')=', xoffsv_hlm(ishot) 
              !print *,'yoffsv_hlm(',ishot,')=', yoffsv_hlm(ishot) 
              !print *,'zoffsv_hlm(',ishot,')=', zoffsv_hlm(ishot) 
              !print *,'elcurv_hlm(',ishot,')=', elcurv_hlm(ishot) 
              !print *,'lengthv_hlm(',ishot,')=', lengthv_hlm(ishot) 

              call ecoute('hlm_reaphy')
           enddo
        else if (words(1) == 'SITES') then                   !Sites
           call ecoute('hlm_reaphy')
           do while (words(1) /= 'ENDSI')
              isite = int(param(1))
              !print *,'reading site ',isite
              if (isite < 1 .or. isite > nsite_hlm) call runend('HLM_REAPHY: WRONG SITE')
              site_hlm(1,isite) = param(2)
              site_hlm(2,isite) = param(3)
              site_hlm(3,isite) = param(4)
              call ecoute('hlm_reaphy')
           enddo
        else if (words(1) == 'OBSER') then                   !Observations for all shots in each site
           call ecoute('hlm_reaphy')
           if (nshot_hlm /= 0_ip .and. nsite_hlm/=0) then 
              call hlm_memphy(6_ip)            !Allocate memory for observations
           else
              call runend('HLM_REAPHY: SHOTS_NUMBER AND SITES_NUMBER MUST BE DEFINED (OBSERVATIONS)')
           end if

           do while (words(1) /= 'ENDOB')
              ishot = int(param(1))
              isite = int(param(2))
              !print *,'reading site ',isite
              if (ishot < 1 .or. ishot > nshot_hlm) call runend('HLM_REAPHY: WRONG SHOT (OBSERVATIONS)')
              if (isite < 1 .or. isite > nsite_hlm) call runend('HLM_REAPHY: WRONG SITE (OBSERVATIONS)')
              smgvpX_obs(ishot,isite) = cmplx(real(param(3)),real(param(4)),kind=rp)
              call ecoute('hlm_reaphy')
           enddo
        else if (words(1) == 'INCID') then      !Incidence matrix: sites associated to each shot

           call hlm_memphy(7_ip)            !Allocate memory for incidence matrix

           if(exists('ALLTO')) then
              incidence_obs(:,:)=1_ip
           else
              incidence_obs(:,:)=0_ip
              call ecoute('hlm_reaphy')
              do while (words(1) /= 'ENDIN')
                 ishot = int(param(1))
                 if (ishot < 1 .or. ishot > nshot_hlm) call runend('HLM_REAPHY: WRONG SHOT (INCIDENCE)')
                 ii = 2_ip
                 do while(param(ii)/=0.0_rp)
                    isite = int(param(ii))
                    if (isite < 1 .or. isite > nsite_hlm) call runend('HLM_REAPHY: WRONG SITE (INCIDENCE)')
                    incidence_obs(ishot,isite) = 1_ip
                    ii=ii+1_ip
                 enddo
                 if (ii==2_ip) call runend('HLM_REAPHY: NO SITES ASSOCIATED TO A SHOT (INCIDENCE)')
                 call ecoute('hlm_reaphy')
              enddo
           end if
        endif
     enddo

     return ! Vladimir

     if (emmet_hlm <= 3 .and. ppcod_hlm <= 2) then
        call hlm_reappo()       !Read in values of primary vector potential
     endif
     !    call hlm_memphy(2_ip)       !For testing Primary Potentials

     do imate = 1,nmate
        do ii=1,ncond_hlm
           dsigma_hlm(ii,imate) = sigma_hlm(ii,imate) - bckco_hlm(ii)       !dsig = Electric conductivity - Background electric conductivity
        enddo
     enddo

  endif

end subroutine hlm_reaphy

