subroutine hlm_parall(order)

  !---------------------------------------------------------------------------
  ! Sources/modules/helmoz/hlm_parall.f90
  ! NAME
  !    hlm_parall
  ! DESCRIPTION
  !    This routine exchanges module specific data for a parallel execution.
  ! USES
  ! USED BY
  !    hlm_turnon
  !----------------------------------------------------------------------------

  use def_parame
  use def_master
  use def_domain
  use def_solver
  use def_helmoz
  use mod_memchk
  use mod_opebcs

  implicit none

  integer(ip), intent(in) :: order

  integer(ip)             :: ii,jj,ji,ki,kfl_ptask_old,isour
  integer(ip)             :: ireac,icoef,dummi
  integer(4)              :: istat



  if (ISEQUEN) return

  select case (order)
  case (1_ip)     
     !Exchange data read from 'hlm_reaphy', 'hlm_reanut' and 'hlm_reaous' files
     kfl_ptask_old = kfl_ptask
     kfl_ptask     = 1_ip
     call Parall(29_ip)
     do parii = 1_ip,2_ip 
        npari = 0_ip
        nparr = 0_ip
        nparc = 0_ip
        nparx = 0_ip
        !Exchange of 'hlm_reaphy' variables 
        call iexcha(emmet_hlm)       !Natural (MT) or controlled (CSEM) source method
        call iexcha(ppcod_hlm)       !Way to read the potentials
        call iexcha(ppout_hlm)       !Way to output the potentials
        call iexcha(nequs_hlm)       !Number of equations
        call iexcha(ncond_hlm)       !Number of elements in conductivity tensor
        call rexcha(frequ_hlm)       !f = Frequency
        call rexcha(anguf_hlm)       !w = Angular frequency 2*pi*f
        call rexcha(xoffs_hlm)       !Source offsets
        call rexcha(yoffs_hlm)
        call rexcha(zoffs_hlm)
        call rexcha(elcur_hlm)       ! I = Electric current
        call rexcha(length_hlm)      ! dl = Dipole length
        call rexcha(airpl_hlm)       ! Air-plane that separates MLSI regions

        call iexcha(nshot_hlm)       ! Number of shots (transmiters)
        call iexcha(kfl_edges_hlm)   ! Edge elements

        call iexcha(nsite_hlm)       ! Number of sites for sources (receivers)



        !Exchange of 'hlm_reanut' variables 
        solve_sol => solve
        !if(kfl_servi(ID_OPTSOL)==1) then
        solad_sol => solad
        !end if
        call soldef(1_ip)
        call iexcha(solve_sol(1)%ndofn)
        !if(kfl_servi(ID_OPTSOL)==1) then
        call iexcha(solad_sol(1)%ndofn)
        !end if
        !Exchange data read from 'hlm_reaous'
        call posdef(1_ip,dummi)
        !Allocate memory for the first pass
        if (parii == 1_ip) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'parin','hlm_parall',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'parre','hlm_parall',parre)
           allocate(parcx(nparx),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'parcx','hlm_parall',parcx)
           if (ISLAVE .or. kfl_ptask == 2_ip) call Parall(2_ip)
        endif
     enddo

     if(IMASTER .and. kfl_ptask /= 2_ip) call Parall(2_ip)
     call memchk(two,istat,mem_servi(1:2,servi),'parin','hlm_parall',parin)
     deallocate(parin,stat=istat)
     if (istat /= 0_ip) call memerr(two,'parin','hlm_parall',0_ip)
     call memchk(two,istat,mem_servi(1:2,servi),'parre','hlm_parall',parre)
     deallocate(parre,stat=istat)
     if (istat /= 0_ip) call memerr(two,'parre','hlm_parall',0_ip)     
     call memchk(two,istat,mem_servi(1:2,servi),'parcx','hlm_parall',parcx)
     deallocate(parcx,stat=istat)
     if (istat /= 0_ip) call memerr(two,'parcx','hlm_parall',0_ip) 

     !Exchange of 'hlm_reaphy' variables which depend on what is read from 'hlm_reaphy' file
     do parii = 1_ip,2_ip 
        npari = 0_ip
        nparr = 0_ip
        nparc = 0_ip
        nparx = 0_ip
        if (emmet_hlm <= 3_ip) then
           call iexcha(nz_hlm)          !Number of points in the z direction, nz
           call iexcha(nr_hlm)          !Number of points in the r direction, nr                               
        endif
        !Allocate memory for the first pass
        if (parii == 1_ip) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'parin','hlm_parall',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'parre','hlm_parall',parre)
           allocate(parcx(nparx),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'parcx','hlm_parall',parcx)
           if (ISLAVE .or. kfl_ptask == 2_ip) call Parall(2_ip)
        endif
     enddo

     if(IMASTER .and. kfl_ptask /= 2_ip) call Parall(2_ip)
     call memchk(two,istat,mem_servi(1:2,servi),'parin','hlm_parall',parin)
     deallocate(parin,stat=istat)
     if (istat /= 0_ip) call memerr(two,'parin','hlm_parall',0_ip)
     call memchk(two,istat,mem_servi(1:2,servi),'parre','hlm_parall',parre)
     deallocate(parre,stat=istat)
     if (istat /= 0_ip) call memerr(two,'parre','hlm_parall',0_ip)     
     call memchk(two,istat,mem_servi(1:2,servi),'parcx','hlm_parall',parcx)
     deallocate(parcx,stat=istat)
     if (istat /= 0_ip) call memerr(two,'parcx','hlm_parall',0_ip)   

     !Allocate memory for allocatable arrays 
     call Parall(30_ip)
     if (ISLAVE) then
        call hlm_memphy(1_ip)           !Allocate memory for material properties
        if (emmet_hlm <= 3_ip) then
           call hlm_memphy(2_ip)       !Allocate memory for primary vector potential
        endif

        call hlm_memphy(3_ip)           !Allocate memory for sites
        call hlm_memphy(5_ip)           !Allocate memory for shots

        call hlm_memphy(6_ip) ! Allocate memory for cost function eval

        call hlm_memphy(7_ip) ! Allocate memory for incidence matrix

     endif
     !call hlm_memphy(6_ip) ! Allocate memory for cost function eval

     kfl_ptask = kfl_ptask_old

     !Exchange physical properties whose dimensions depend on what is read from 'hlm_reaphy' file
     do parii = 1_ip,2_ip 
        npari = 0_ip
        nparr = 0_ip
        nparc = 0_ip
        nparx = 0_ip
        !Exchange of 'hlm_reaphy' variables whose dimensions depend on what is read in 'hlm_reaphy' file
        do ii=1,ncond_hlm
           call rexcha(bckco_hlm(ii))       !Background conductivity
        enddo

        do ii=1,nshot_hlm
           call rexcha(xoffsv_hlm(ii))     !x coordinate of shot
           call rexcha(yoffsv_hlm(ii))     !y coordinate of shot
           call rexcha(zoffsv_hlm(ii))     !z coordinate of shot
           call rexcha(elcurv_hlm(ii))     !Dipole current
           call rexcha(lengthv_hlm(ii))    !Dipole length
        enddo


        do ii=1,nsite_hlm
           call rexcha(site_hlm(1,ii))     !x coordinate of site
           call rexcha(site_hlm(2,ii))     !y coordinate of site
           call rexcha(site_hlm(3,ii))     !z coordinate of site
        enddo


        do ii=1,(nsite_hlm*1_ip)
           call rexcha(clsite_hlm(1,ii))   !x coordinate of closest node to site
           call rexcha(clsite_hlm(2,ii))   !y coordinate of closest node to site
           call rexcha(clsite_hlm(3,ii))   !z coordinate of closest node to site

           call iexcha(clsite1_hlm(ii))    !Global node number of closest node to site
           call rexcha(clsite2_hlm(ii))                 
        enddo

        if(kfl_servi(ID_OPTSOL)==1) then
           do ii=1,nshot_hlm
              do jj=1,nsite_hlm
                 call xexcha(smgvpX_obs(ii,jj))    
              enddo
           enddo
        end if

        do ii=1,nshot_hlm
           do jj=1,nsite_hlm
              call iexcha(incidence_obs(ii,jj))   ! 1 or 0
           enddo
        enddo

        do ji = 1,nmate
           call rexcha(perma_hlm(ji))
           call rexcha(epsil_hlm(ji))
           do ii=1,ncond_hlm
              call rexcha(sigma_hlm(ii,ji))
              call rexcha(dsigma_hlm(ii,ji))
           enddo
        enddo
        if (emmet_hlm <= 3_ip) then
           if (ppcod_hlm == 1) then
              do ii = 1,nz_hlm
                 call rexcha(z_hlm(ii))
              enddo
              do ii = 1,nr_hlm
                 call rexcha(r_hlm(ii))
              enddo
              do jj = 1,nr_hlm
                 do ii = 1,nz_hlm
                    call xexcha(pvepo_hlm(ii,jj))
                 enddo
              enddo
           endif
           if (ppcod_hlm == 2) then
              do jj = 1,npoin
                 do ii = 1,ndime
                    call xexcha(pmgvp_hlm(ii,jj))
                 enddo
                 call xexcha(pelsp_hlm(jj))
              enddo
           endif
        endif

        if (parii == 1_ip) then
           allocate(parin(npari),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'parin','hlm_parall',parin)
           allocate(parre(nparr),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'parre','hlm_parall',parre)
           allocate(parcx(nparx),stat=istat)
           call memchk(zero,istat,mem_servi(1:2,servi),'parcx','hlm_parall',parcx)
           if (ISLAVE .or. kfl_ptask == 2_ip) call Parall(2_ip)
        endif
     enddo
     !write (*,*) 'I am ', kfl_paral, emmet_hlm, ppcod_hlm, anguf_hlm, bckco_hlm(1), bckco_hlm(2), bckco_hlm(3), xoffs_hlm, yoffs_hlm, zoffs_hlm, elcur_hlm, length_hlm

     if (IMASTER .and. kfl_ptask /= 2_ip) call Parall(2_ip)
     call memchk(two,istat,mem_servi(1:2,servi),'parin','hlm_parall',parin)
     deallocate(parin,stat=istat)
     if ( istat /= 0_ip ) call memerr(two,'parin','hlm_parall',0_ip)
     call memchk(two,istat,mem_servi(1:2,servi),'parre','hlm_parall',parre)
     deallocate(parre,stat=istat)
     if ( istat /= 0_ip ) call memerr(two,'parre','hlm_parall',0_ip)     
     call memchk(two,istat,mem_servi(1:2,servi),'parcx','hlm_parall',parcx)
     deallocate(parcx,stat=istat)
     if ( istat /= 0_ip ) call memerr(two,'parcx','hlm_parall',0_ip)   

     !Exchange boundary conditions
     call spnbcs(tncod_hlm)
  end select


end subroutine hlm_parall
