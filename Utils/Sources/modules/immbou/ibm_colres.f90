subroutine ibm_colres(ittib)
  !-----------------------------------------------------------------------
  !****f* ibm_colres/ibm_colres
  ! NAME
  !    ibm_colres
  ! DESCRIPTION
  !    This routines computes collision
  !    In ibm_coldet, we determined the minimum time step during which we have a collision
  !    Now, the particles that have collision can influence that particles that
  !    in ibm_coldet were not supposed to collision.
  !    To check this, we compute the influence radius of the particle that
  !    has collisioned. This will affect the particles that are located inside this 
  !    influence radius. Influence radius: mbbox
  !    We have computed in ibm_newmak:
  !    - new position
  !    - new velocity
  !    - new acceleration 
  ! USED BY
  !***
  !-----------------------------------------------------------------------
  use def_kintyp
  use def_master
  use def_domain
  use def_immbou
  use mod_memchk
  use mod_messages, only : livinf
  implicit none
  integer(ip), intent(in) :: ittib
  integer(ip)             :: iimbo,jimbo,iwaib,idime,ifind
  integer(ip)             :: icond,icoll,icolw,kk,ii
  integer(ip)             :: grows,icont,itera
  integer(ip)             :: group(nimbo),igrou,jgrou,ngrou,ninvi
  integer(ip)             :: ibodi,jbodi,tbodi,nbodi(nimbo)
  integer(ip)             :: iwall,twall,nwall(nimbo)
  integer(ip)             :: lcoll(nimbo)
  integer(4)              :: istat

  real(rp)                :: dt,toler,vrel
  real(rp)                :: dista,coori(3),coorj(3),norma(3)
  real(rp)                :: boboi(3,2),boboj(3,2)
  real(rp)                :: cotim(nimbo + nwaib)
  real(rp)                :: ri(3),rj(3),wtoti(3),wtotj(3)
  real(rp), target        :: dummr(2)
  real(rp), pointer       :: lwall(:)

  type cotyp
     integer(ip)          :: igrou
     integer(ip)          :: iimbo
     integer(ip)          :: jimbo
     real(rp)             :: norma(3)
     real(rp)             :: coori(3)
     real(rp)             :: coorj(3)
  end type cotyp
  type(cotyp), pointer    :: lbodi(:)

  if( kfl_colli_ibm /= 0 ) then     
     !
     ! Initialization
     ! 
     dt    = cutim - cutim_ibm - dtime_ibm
     toler = 2.5e-2_rp
     icoll = 0
     icolw = 0
     do iimbo = 1,nimbo 
        cotim(iimbo) = imbou(iimbo) % cotim
        lcoll(iimbo) = 0_ip
     end do
     do iwaib = 1,nwaib  
        cotim(nimbo + iwaib) = twall_ibm(iwaib) % cotim
     end do

     !----------------------------------------
     !
     ! Find groups of interconnected particles
     !
     !----------------------------------------
     grows           = 1_ip
     ngrou           = 1_ip
     ninvi           = 0_ip
     do iimbo=1,nimbo
        group(iimbo) = 0_ip        
     end do
     group(1)        = 1_ip
     nbodi(1)        = 0_ip
     nwall(1)        = 0_ip
     do while (grows == 1_ip)
        grows = 0_ip     
        do iimbo=1,nimbo
           !
           ! Check new contacts for the current group
           !
           if ( group(iimbo) == ngrou ) then
              if ( imbou(iimbo) % cotim <= cutim_ibm + dtime_ibm + 0.1_rp*dtime_ibm ) then
                 group(iimbo) = -ngrou
                 do idime = 1,3
                    boboi(idime,1) = imbou(iimbo) % bobox(idime,1) - abs(imbou(iimbo)%maxdi*toler)
                    boboi(idime,2) = imbou(iimbo) % bobox(idime,2) + abs(imbou(iimbo)%maxdi*toler)
                 end do
                 !
                 ! IIMBO vs. JIMBO PARTICLES
                 !
                 do jimbo=1,nimbo
                    if ( group(jimbo) == 0_ip .or. group(jimbo) == ngrou ) then
                       if ( imbou(jimbo) % cotim <= cutim_ibm + dtime_ibm + 0.1_rp*dtime_ibm ) then
                          do idime = 1,ndime
                             boboj(idime,1) = imbou(jimbo) % bobox(idime,1) - abs(imbou(jimbo)%maxdi*toler)
                             boboj(idime,2) = imbou(jimbo) % bobox(idime,2) + abs(imbou(jimbo)%maxdi*toler)
                          end do
                          icond=0_ip
                          !
                          ! Determine if the bounding boxes are intersected
                          !
                          do idime=1,ndime
                             if ( boboi(idime,1) < boboj(idime,2) .and. boboi(idime,2) > boboj(idime,1) ) then
                                icond=icond+1_ip
                             end if
                          end do
                          if ( icond==ndime ) then            
                             !
                             ! Check only the particles that are not checked before
                             !                             
                             call ibm_dpapar(iimbo,jimbo,coori,coorj,norma,dista)
                             if ( dista < toler*min(imbou(iimbo)%maxdi,imbou(jimbo)%maxdi) ) then
                                grows        = 1_ip
                                group(jimbo) = ngrou
                                nbodi(ngrou) = nbodi(ngrou)  + 1_ip
                             end if
                          end if
                       end if
                    end if
                 end do
                 !
                 ! IIMBO PARTCILE vs. IWAIB WALL
                 !
                 do iwaib = 1,nwaib
                    if ( twall_ibm(iwaib) % cotim <= cutim_ibm + dtime_ibm + 0.1_rp*dtime_ibm ) then                       
                       icond=0_ip
                       !
                       ! Determine if the bounding boxes are intersected
                       !
                       do idime=1,ndime           
                          if ( boboi(idime,1) < twall_ibm(iwaib) % bobox(idime,2) &
                               .and. boboi(idime,2) > twall_ibm(iwaib) % bobox(idime,1) ) then
                             icond=icond+1_ip
                          end if
                       end do
                       if( icond == ndime ) then
                          call ibm_boxwal(boboi,iwaib,ifind)                    
                          if ( ifind == 1 ) then
                    
                             !
                             ! Determine the minimum distance betwwen the particle and the wall
                             ! 
                             call ibm_dpawal(iimbo,iwaib,coori,coorj,norma,dista) 
                             if ( dista < toler*imbou(iimbo)%maxdi ) then
                                grows        = 1_ip
                                nwall(ngrou) = nwall(ngrou) + 1_ip
                             end if
                          end if
                       end if
                    end if
                 end do

              end if
           end if
        end do
        !
        ! Find a new group
        !     
        if (grows == 0) then
           do iimbo = 1,nimbo              
              if (abs(group(iimbo)) == ngrou) then
                 group(iimbo) = abs(group(iimbo))
              end if
           end do
           if ( nbodi(ngrou) == 0 )  then
              if ( nwall(ngrou) /= 0 ) then
                 do iimbo = 1,nimbo
                    if (group(iimbo) == ngrou) then
                       group(iimbo) = nimbo+1_ip
                    end  if
                 end do
                 ninvi        = ninvi + nwall(ngrou)
                 nwall(ngrou) = 0_ip 
              else
                 do iimbo = 1,nimbo
                    if (group(iimbo) == ngrou) then
                       group(iimbo) = -nimbo-1_ip
                    end  if
                 end do
              end if
              ngrou = ngrou - 1_ip
           end if
           iimbo   = 0_ip
           do while (iimbo < nimbo .and. grows == 0)
              iimbo = iimbo + 1_ip
              if (group(iimbo) == 0_ip) then
                 grows        = 1_ip
                 ngrou        = ngrou + 1_ip
                 group(iimbo) = ngrou
                 nbodi(ngrou) = 0_ip
                 nwall(ngrou) = 0_ip
              end if
           end do
        end if
     end do

     !----------------------------------------------
     !
     ! Obtain contact information
     !     
     !----------------------------------------------

     tbodi = 0_ip
     twall = 0_ip
     do igrou = 1,ngrou
        tbodi = tbodi + nbodi(igrou)
        twall = twall + nwall(igrou)  
     end do
     if ( tbodi > 0 ) allocate( lbodi(tbodi) )
     if ( twall+ninvi > 0 ) allocate( lwall(twall*8 + ninvi*8) )
     iwall = 0_ip     
     ibodi = 0_ip
     do igrou = 1,ngrou
        !
        ! IIMBO PARTCILE vs. IWAIB WALL
        !
        do iimbo=1,nimbo
           if (group(iimbo) == igrou) then
              do idime = 1,3
                 boboi(idime,1) = imbou(iimbo) % bobox(idime,1) - abs(imbou(iimbo)%maxdi*toler)
                 boboi(idime,2) = imbou(iimbo) % bobox(idime,2) + abs(imbou(iimbo)%maxdi*toler)
              end do
              do iwaib = 1,nwaib
                 icond=0_ip
                 !
                 ! Determine if the bounding boxes are intersected
                 !
                 do idime=1,ndime           
                    if ( boboi(idime,1) < twall_ibm(iwaib) % bobox(idime,2) .and. & 
                         boboi(idime,2) > twall_ibm(iwaib) % bobox(idime,1) ) then
                       icond=icond+1_ip
                    end if
                 end do                    
                 if( icond == ndime ) then
                    
                    call ibm_boxwal(boboi,iwaib,ifind)                    
                    if ( ifind == 1 ) then

                       !
                       ! Determine the minimum distance betwwen the particle and the wall
                       !  
                       call ibm_dpawal(iimbo,iwaib,coori,coorj,norma,dista) 
                       if ( dista < toler*imbou(iimbo)%maxdi ) then                          
                          iwall = iwall + 1_ip

                          lwall( (iwall-1)*8 + 1 ) = real(igrou)
                          lwall( (iwall-1)*8 + 2 ) = real(iimbo)
                          lwall( (iwall-1)*8 + 3 ) = norma(1)
                          lwall( (iwall-1)*8 + 4 ) = norma(2)
                          lwall( (iwall-1)*8 + 5 ) = norma(ndime)
                          lwall( (iwall-1)*8 + 6 ) = coori(1)
                          lwall( (iwall-1)*8 + 7 ) = coori(2)
                          lwall( (iwall-1)*8 + 8 ) = coori(ndime)

                       end if
                    end if
                 end if
              end do
              !
              ! IIMBO vs. JIMBO PARTICLES
              !
              do jimbo=iimbo+1,nimbo                 
                 if (group(jimbo) == igrou) then
                    do idime = 1,3
                       boboj(idime,1) = imbou(jimbo) % bobox(idime,1) - abs(imbou(jimbo)%maxdi*toler)
                       boboj(idime,2) = imbou(jimbo) % bobox(idime,2) + abs(imbou(jimbo)%maxdi*toler)
                    end do
                    icond=0_ip
                    !
                    ! Determine if the bounding boxes are intersected
                    !
                    do idime=1,ndime
                       if ( boboi(idime,1) < boboj(idime,2) .and. boboi(idime,2) > boboj(idime,1) ) then
                          icond=icond+1_ip
                       end if
                    end do
                    if ( icond==ndime ) then      
                       call ibm_dpapar(iimbo,jimbo,coori,coorj,norma,dista)
                       if ( dista < toler*min(imbou(iimbo)%maxdi,imbou(jimbo)%maxdi) ) then

                          ibodi = ibodi + 1_ip
                          lbodi(ibodi)%igrou = igrou
                          lbodi(ibodi)%iimbo = iimbo
                          lbodi(ibodi)%jimbo = jimbo
                          do idime = 1,ndime
                             lbodi(ibodi)%norma(idime) = norma(idime)
                             lbodi(ibodi)%coori(idime) = coori(idime)
                             lbodi(ibodi)%coorj(idime) = coorj(idime)
                          end do
                       end if
                    end if
                 end if
              end do
           end if
        end do
     end do

     !----------------------------------------------
     !
     ! Obtain invisible contact information
     !     
     !----------------------------------------------
     iwall = twall
     do iimbo=1,nimbo
        if (group(iimbo) == nimbo+1_ip) then
           do idime = 1,3
              boboi(idime,1) = imbou(iimbo) % bobox(idime,1) - abs(imbou(iimbo)%maxdi*toler)
              boboi(idime,2) = imbou(iimbo) % bobox(idime,2) + abs(imbou(iimbo)%maxdi*toler)
           end do
           do iwaib = 1,nwaib
              if ( twall_ibm(iwaib) % cotim <= cutim_ibm + dtime_ibm + 0.1_rp*dtime_ibm ) then                       

                 icond=0_ip
                 !
                 ! Determine if the bounding boxes are intersected
                 !
                 do idime=1,ndime           
                    if ( boboi(idime,1) < twall_ibm(iwaib) % bobox(idime,2) .and. & 
                         boboi(idime,2) > twall_ibm(iwaib) % bobox(idime,1) ) then
                       icond=icond+1_ip
                    end if
                 end do                    
                 if( icond == ndime ) then
                    call ibm_boxwal(boboi,iwaib,ifind)                    
                    if ( ifind == 1 ) then
                       !
                       ! Determine the minimum distance betwwen the particle and the wall
                       !  
                       call ibm_dpawal(iimbo,iwaib,coori,coorj,norma,dista)     
                                   
                       if ( dista < toler*imbou(iimbo)%maxdi ) then      
                          iwall = iwall + 1_ip
                          lwall( (iwall-1)*8 + 1 ) = real(nimbo+1_ip)
                          lwall( (iwall-1)*8 + 2 ) = real(iimbo)
                          lwall( (iwall-1)*8 + 3 ) = norma(1)
                          lwall( (iwall-1)*8 + 4 ) = norma(2)
                          lwall( (iwall-1)*8 + 5 ) = norma(ndime)
                          lwall( (iwall-1)*8 + 6 ) = coori(1)
                          lwall( (iwall-1)*8 + 7 ) = coori(2)
                          lwall( (iwall-1)*8 + 8 ) = coori(ndime)
                       end if                       
                    end if
                 end if
              end if
           end do
        end if
     end do

     !---------------------------------------------------------------------------------------------
     !
     ! Parall
     ! - Each slaves have a list of walls
     !              
     !---------------------------------------------------------------------------------------------

     twall = twall + ninvi
     
     if ( IPARALL ) then
        kk = twall      
        call parari('SUM',0_ip,1_ip,kk)
        
        if( kk > 0 ) then
           
           allocate( parig(npart+1) , stat = istat )    
           parin => parig
           call parari('AGA',0_ip,1_ip,twall)
           
           kk = 0
           do ii = 1,npart+1
              parig(ii) = parig(ii) * 8
              kk        = kk + parig(ii)
           end do
           
           if( kk > 0 ) then
              allocate( pari1(npart+1) , stat = istat )    
              pari1(1) = 0
              do ii = 2,npart+1
                 pari1(ii) = pari1(ii-1) + parig(ii-1)
              end do
              allocate( parre(max(1_ip,kk)) , stat = istat )
              npasr =  twall * 8
              if( npasr == 0 ) then
                 parrs => dummr
              else
                 parrs => lwall
              end if
              
              call Parall(47_ip)
              
              twall = kk / 8
              deallocate( lwall , stat = istat )
              allocate( lwall(kk) , stat = istat )
              do ii = 1,kk
                 lwall(ii) = parre(ii)
              end do
              deallocate( parre )
              deallocate( pari1 )
              
           end if
           deallocate( parig )
        end if
     end if
     !---------------------------------------------------------------------------------------------
     !
     ! Solves the contacts for each pair of particles until all the contact in the group are solved
     !              
     !---------------------------------------------------------------------------------------------


     jbodi = 0_ip
     do igrou = 1,ngrou
        icont = 1_ip
        itera = 0_ip
        do while (icont == 1 .and. itera < 1000)
           itera = itera + 1_ip
           icont = 0_ip
           do iwall = 1,twall
              jgrou = int(lwall( (iwall-1_ip)*8 + 1 ))
              if (igrou == jgrou) then
                 iimbo = int(lwall( (iwall-1)*8 + 2 ))
                 do idime = 1,ndime
                    norma(idime) = lwall( (iwall-1)*8 + 2 + idime )
                    coori(idime) = lwall( (iwall-1)*8 + 5 + idime )
                 end do
                 ri(3) =  0.0_rp
                 do idime = 1,ndime
                    ri(idime) = coori(idime) - imbou(iimbo)%posil(idime,1)
                 end do
                 call vecpro(imbou(iimbo) % veloa,ri,wtoti,3_ip)
                 vrel = 0.0_rp
                 do idime = 1,3
                    ! Calculate the relative velocity between the particles
                    vrel = vrel + norma(idime)*(imbou(iimbo) % velol(idime,1) + wtoti(idime))
                 end do
                 if ( vrel < 0.0_rp ) then            
       
                    icont = 1_ip
                    if( INOTSLAVE ) then
                       write(momod(modul) % lun_outpu,'(a)')          ' *****************************************************************'
                       write(momod(modul) % lun_outpu,'(a)')          '             A COLLISION HAS BEEN FOUNDED (Particle-Wall)'
                       write(momod(modul) % lun_outpu,'(a)')          ' *****************************************************************'
                       write(momod(modul) % lun_outpu,'(a,i4,1x,i4)') '   PARTICLES INVOLVED:              ',iimbo, iwaib
                       write(momod(modul) % lun_outpu,'(a,e12.6)')    '   TIME OF COLLISION:               ',cutim_ibm + dtime_ibm
                       !write(momod(modul) % lun_outpu,'(a,e12.6)')    '   DISTANCE BEETWEEN THE PARTICLES: ',dista              
                       write(momod(modul) % lun_outpu,'(a,e12.6)')    '   TOLERANCY:                       ',toler*imbou(iimbo)%maxdi
                       write(momod(modul) % lun_outpu,'(a,e12.6)')    '   RELATIVE VELOCITY:               ',vrel
                       write(momod(modul) % lun_outpu,'(a)')          ' *****************************************************************'
                       flush(momod(modul) % lun_outpu)
                    end if
                    !
                    ! Calculate the impulse and the velocities of the particles after the collision
                    !
                    call ibm_ipawal(iimbo,coori,norma)
                    lcoll(iimbo) = 1_ip
                 end if
              end if
           end do
           jbodi = 0_ip
           do jgrou = 1,igrou-1
              jbodi = jbodi + nbodi(jgrou)
           end do
           do ibodi = 1,nbodi(igrou)
              iimbo = lbodi(ibodi+jbodi) % iimbo
              jimbo = lbodi(ibodi+jbodi) % jimbo
              do idime = 1,ndime
                 norma(idime) = lbodi(ibodi+jbodi)%norma(idime)
                 coori(idime) = lbodi(ibodi+jbodi)%coori(idime)
                 coorj(idime) = lbodi(ibodi+jbodi)%coorj(idime)
              end do
              ri(3) =  0.0_rp
              rj(3) =  0.0_rp
              do idime = 1,ndime
                 ri(idime) = coori(idime) - imbou(iimbo)%posil(idime,1) 
                 rj(idime) = coorj(idime) - imbou(jimbo)%posil(idime,1)
              end do
              call vecpro(imbou(iimbo)%veloa,ri,wtoti,3_ip)
              call vecpro(imbou(jimbo)%veloa,rj,wtotj,3_ip)
              vrel = 0.0_rp
              do idime = 1,3
                 ! Calculate the relative velocity between the particles
                 vrel = vrel + norma(idime)*(imbou(iimbo)%velol(idime,1) + wtoti(idime) - imbou(jimbo)%velol(idime,1) - wtotj(idime))
              end do
              if ( vrel < 0.0_rp ) then
                 icont = 1_ip
                 if( INOTSLAVE ) then
                    write(momod(modul) % lun_outpu,'(a)')          ' *****************************************************************'
                    write(momod(modul) % lun_outpu,'(a)')          '            A COLLISION HAS BEEN FOUNDED (Particle-Particle)'
                    write(momod(modul) % lun_outpu,'(a)')          ' *****************************************************************'
                    write(momod(modul) % lun_outpu,'(a,i4,1x,i4)') '   PARTICLES INVOLVED:              ',iimbo, jimbo
                    write(momod(modul) % lun_outpu,'(a,e12.6)')    '   TIME OF COLLISION:               ',cutim_ibm + dtime_ibm
                    !write(momod(modul) % lun_outpu,'(a,e12.6)')    '   DISTANCE BEETWEEN THE PARTICLES: ',dista              
                    write(momod(modul) % lun_outpu,'(a,e12.6)')    '   TOLERANCY:                       ',toler* &
                         min(imbou(iimbo)%maxdi,imbou(jimbo)%maxdi)
                    write(momod(modul) % lun_outpu,'(a,e12.6)')    '   RELATIVE VELOCITY:               ',vrel
                    write(momod(modul) % lun_outpu,'(a)')          ' *****************************************************************'
                    flush(momod(modul) % lun_outpu)
                 end if
                 !
                 ! Calculate the impulse and the velocities of the particles after the collision
                 !
                 call ibm_ipapar(iimbo,jimbo,coori,coorj,norma)
                 lcoll(iimbo) = 1_ip
                 lcoll(jimbo) = 1_ip
              end if
           end do
        end do
     end do
     !---------------------------------------------------------------------------------------------
     !
     ! Solves contacts involves with the invisible groups
     !              
     !---------------------------------------------------------------------------------------------

     do iimbo = 1,nimbo 
        icont = 1_ip
        itera = 0_ip
        do while (icont == 1 .and. itera < 1000)
           itera = itera + 1_ip
           icont = 0_ip        
           do iwall = 1,twall
              igrou = int(lwall( (iwall-1_ip)*8_ip + 1_ip ))
              jimbo = int(lwall( (iwall-1_ip)*8_ip + 2_ip ))
              if (igrou > nimbo .and. iimbo == jimbo) then
                 do idime = 1,ndime
                    norma(idime) = lwall( (iwall-1)*8 + 2 + idime )
                    coori(idime) = lwall( (iwall-1)*8 + 5 + idime )
                 end do
                 ri(3) =  0.0_rp
                 do idime = 1,ndime
                    ri(idime) = coori(idime) - imbou(iimbo)%posil(idime,1)
                 end do
                 call vecpro(imbou(iimbo) % veloa,ri,wtoti,3_ip)
                 vrel = 0.0_rp
                 do idime = 1,3
                    ! Calculate the relative velocity between the particles
                    vrel = vrel + norma(idime)*(imbou(iimbo) % velol(idime,1) + wtoti(idime))
                 end do                 
                 if ( vrel < 0.0_rp ) then
                    icont = 1_ip
                    if( INOTSLAVE ) then                       
                       write(momod(modul) % lun_outpu,'(a)')          ' *****************************************************************'
                       write(momod(modul) % lun_outpu,'(a)')          '             A COLLISION HAS BEEN FOUNDED (Particle-Wall)'
                       write(momod(modul) % lun_outpu,'(a)')          ' *****************************************************************'
                       write(momod(modul) % lun_outpu,'(a,i4,1x,i4)') '   PARTICLES INVOLVED:              ',iimbo, iwaib
                       write(momod(modul) % lun_outpu,'(a,e12.6)')    '   TIME OF COLLISION:               ',cutim_ibm + dtime_ibm
                       !write(momod(modul) % lun_outpu,'(a,e12.6)')    '   DISTANCE BEETWEEN THE PARTICLES: ',dista              
                       write(momod(modul) % lun_outpu,'(a,e12.6)')    '   TOLERANCY:                       ',toler*imbou(iimbo)%maxdi
                       write(momod(modul) % lun_outpu,'(a,e12.6)')    '   RELATIVE VELOCITY:               ',vrel
                       write(momod(modul) % lun_outpu,'(a)')          ' *****************************************************************'
                       flush(momod(modul) % lun_outpu)
                    end if
                    !
                    ! Calculate the impulse and the velocities of the particles after the collision
                    !
                    call ibm_ipawal(iimbo,coori,norma)
                    lcoll(iimbo) = 1_ip
                 end if
              end if
           end do
        end do
     end do
     !if (twall > 0) pause

     if ( tbodi > 0 ) deallocate( lbodi )
     if ( twall > 0 ) deallocate( lwall )
     !
     ! Count number of collisions
     !
     if( icoll /= 0 ) then
        ioutp(1) = icoll
        call livinf(99_ip,'COLLISION BETWEEN PARTICLES AND PARTICLES:',0_ip)
     end if
     call parari('SUM',0_ip,1_ip,icolw)
     if( icolw /= 0 ) then
        ioutp(1) = icolw
        call livinf(99_ip,'COLLISION BETWEEN PARTICLES AND WALLS:',0_ip)
     end if

     !-----------------------------------------------------------------------------------------
     !
     ! Force to check the time of collision again for the particles whose velocity have changed
     !
     !-----------------------------------------------------------------------------------------
     do iimbo = 1 , nimbo
        if ( lcoll(iimbo) == 1 ) then
           call ibm_movbox(iimbo,1,dt,boboi)
           cotim(iimbo) = -1.0e10_rp
           !
           ! IIMBO vs. JIMBO PARTICLES
           !
           do jimbo = 1,nimbo
              if( iimbo /= jimbo ) then
                 call ibm_movbox(jimbo,1,dt,boboj)
                 icond=0_ip
                 do idime=1,ndime
                    if ( boboj(idime,1) < boboi(idime,2) .and. boboj(idime,2) > boboi(idime,1) ) then
                       icond = icond+1_ip
                    end if
                 end do
                 if (icond == ndime) cotim(jimbo) = -1.0e10_rp

              end if
           end do
           !
           ! IIMBO PARTICLE vs. IWAIB WALL
           !
           do iwaib = 1,nwaib    
              icond=0_ip
              do idime=1,ndime
                 if ( twall_ibm(iwaib) % bobox(idime,1) < boboi(idime,2) .and. twall_ibm(iwaib) % bobox(idime,2) > boboi(idime,1) ) then
                    icond = icond+1_ip
                 end if
              end do
              if (icond == ndime) then 
                 call ibm_boxwal(boboi,iwaib,ifind)                    
                 if ( ifind == 1 ) cotim(nimbo + iwaib) = -1.0e10_rp
              end if
           end do
        end if
     end do
     !
     ! Parall: Force to check the time of collision again 
     !
     call pararr('MIN',0_ip,nimbo,cotim)
     do iimbo = 1,nimbo  
        imbou(iimbo) % cotim = cotim(iimbo)
     end do
     do iwaib = 1,nwaib  
        twall_ibm(iwaib) % cotim = cotim(nimbo + iwaib)
     end do
  end if


end subroutine ibm_colres



!-----------------------------------------------------------------------
! NAME
!    ibm_ipapar
! DESCRIPTION
!    This routines calculate the impulse of the collision and the velocities 
!    between two particles
! USED BY
!    
!-----------------------------------------------------------------------
subroutine ibm_ipapar(iimbo,jimbo,coori,coorj,norma)
  use def_kintyp, only     :  ip,rp
  use def_master, only     :  imbou
  use def_domain, only     :  ndime,nimbo,nnode
  use def_immbou, only     :  kfl_rotib_ibm,dtime_ibm
  implicit none
  integer(ip), intent(in)    :: iimbo,jimbo
  integer(ip)                :: idime

  real(rp),    intent(in)    :: coori(3),coorj(3),norma(3)
  real(rp)                   :: impul(3)
  real(rp)                   :: Io(3,3),Item(3,3),I(3,3),invIi(3,3),invIj(3,3)
  real(rp)                   :: ri(3),rj(3),vti(3),vtj(3),vtem(3),resi(3),resj(3)
  real(rp)                   :: dt
  real(rp)                   :: C,numer,denom,immag,deter
  real(rp),    pointer       :: xi(:,:),ai(:,:),vi(:,:),wi(:,:),ITi(:),Roti(:,:) ! Linear and angular quantities of iimbo
  real(rp),    pointer       :: xj(:,:),aj(:,:),vj(:,:),wj(:,:),ITj(:),Rotj(:,:) ! Linear and angular quantities of jimbo
  real(rp)                   :: Et,v2i,w2i,v2j,w2j
  !
  ! Initialization
  ! Restitution coefficient. When C=1.0, perfect elasticity, if C=0.0, perfect inelasticity
  !
  C     = 0.8_rp  
  dt    = dtime_ibm
  do idime = 1,3
     vti(idime) = 0.0_rp
     vtj(idime) = 0.0_rp
  end do
  xi   => imbou(iimbo)%posil
  vi   => imbou(iimbo)%velol
  ai   => imbou(iimbo)%accel
  wi   => imbou(iimbo)%veloa
  ITi  => imbou(iimbo)%momin
  Roti => imbou(iimbo)%rotac

  xj   => imbou(jimbo)%posil
  vj   => imbou(jimbo)%velol
  aj   => imbou(jimbo)%accel
  wj   => imbou(jimbo)%veloa
  ITj  => imbou(jimbo)%momin
  Rotj => imbou(jimbo)%rotac
  !
  ! Coordinates of contact point respect to the center of mass of the particles 
  !
  ri(3) = 0.0_rp
  rj(3) = 0.0_rp
  resi(3) = 0.0_rp
  resj(3) = 0.0_rp
  do idime=1,ndime
     ri(idime)   = coori(idime) - imbou(iimbo)%posil(idime,1)
     rj(idime)   = coorj(idime) - imbou(jimbo)%posil(idime,1)
     resi(idime) = 0.0_rp
     resj(idime) = 0.0_rp
  end do

  !------------------------------------------------------------------
  !
  ! Calculate the impulse of the collision
  !
  !------------------------------------------------------------------
  if( kfl_rotib_ibm /= 0 ) then
     if (ndime == 3) then

        ! Compute the actual inertia tensor for iimbo
        Io(1,1)=ITi(1); Io(1,2)=ITi(4); Io(1,3)=ITi(5)
        Io(2,1)=ITi(4); Io(2,2)=ITi(2); Io(2,3)=ITi(6)
        Io(3,1)=ITi(5); Io(3,2)=ITi(6); Io(3,3)=ITi(3)        
        call mbmab0(Item,Roti,Io ,3,3,3)
        call mbmabt(I   ,Item,Roti,3,3,3)        
        ! Compute the actual inverse inertia tensor for iimbo
        call invmtx(I   ,invIi,deter,3)
        !
        ! resi = (invIj � (rj x norma)) x rj   (for jimbo)
        !
        call vecpro(ri  ,norma,resi,3)
        call mbvab0(vtem,invIi,resi,3,3) 
        call vecpro(vtem,ri   ,resi,3) 
        !        
        ! Compute the actual inertia tensor for jimbo
        !
        Io(1,1)=ITj(1); Io(1,2)=ITj(4); Io(1,3)=ITj(5)
        Io(2,1)=ITj(4); Io(2,2)=ITj(2); Io(2,3)=ITj(6)
        Io(3,1)=ITj(5); Io(3,2)=ITj(6); Io(3,3)=ITj(3)        
        call mbmab0(Item,Rotj,Io  ,3,3,3)
        call mbmabt(I   ,Item,Rotj,3,3,3)
        !
        ! Compute the actual inverse inertia tensor for jimbo
        !
        call invmtx(I,invIj,deter,3)
        !               
        ! resj = (invIj � (rj x norma)) x rj   (for jimbo)
        !
        call vecpro(rj,  norma,resj,3) 
        call mbvab0(vtem,invIj,resj,3,3) 
        call vecpro(vtem,rj   ,resj,3) 

     elseif (ndime == 2) then
        !
        ! resi = (invIj � (rj x norma)) x rj   (for jimbo)
        !
        call vecpro(ri  ,norma,vtem,3)
        call vecpro(vtem,ri   ,resi,3)
        !
        ! resj = (invIj � (rj x norma)) x rj   (for jimbo)
        !
        call vecpro(rj  ,norma,vtem,3)
        call vecpro(vtem,rj   ,resj,3)        
        do idime=1,3     
           ! Multiply by the inertia, in 2D is a scalar 
           resi(idime) = (1.0_rp/ITi(1)) * resi(idime)
           resj(idime) = (1.0_rp/ITj(1)) * resj(idime)
        end do
     end if
     call vecpro(wi,ri,vti,3)
     call vecpro(wj,rj,vtj,3)
  end if
  !
  ! Calculate the total velocity before the collision
  !  
  numer = 0.0_rp
  denom = 0.0_rp
  do idime=1,3
     vti(idime) = vi(idime,1) + vti(idime)
     vtj(idime) = vj(idime,1) + vtj(idime)
     denom = denom + norma(idime) * (resi(idime) + resj(idime))
     numer = numer + norma(idime) * ( vti(idime) -  vtj(idime))
  end do

  numer = (-1.0_rp-C)*numer
  denom = denom + (1.0_rp/imbou(iimbo)%massa) +  (1.0_rp/imbou(jimbo)%massa)
  !
  ! Magnitude of the impulse
  !
  immag = numer * (1.0_rp/denom)
  !
  ! Impulse vector
  !
  do idime=1,3
     impul(idime) = immag * norma(idime)
  end do
  !------------------------------------------------------------------
  !
  ! Calculate the new velocities of the particles after the collision
  !
  !------------------------------------------------------------------
  resi(3) = 0.0_rp
  resj(3) = 0.0_rp
  do idime=1,ndime
     resi(idime) = 0.0_rp
     resj(idime) = 0.0_rp
  end do
  if (ndime == 3) then
     ! resi = invIi � (ri x impul)  (for iimbo)
     call vecpro(ri,impul,vtem,3)      
     call mbvab0(resi,invIi,vtem,3,3) 

     ! resj = invIj � (rj x impul)   (for jimbo)
     call vecpro(rj,impul,vtem,3) 
     call mbvab0(resj,invIj,vtem,3,3)   
  elseif (ndime == 2) then
     ! resi = invIi � (ri x impul)  (for iimbo)
     call vecpro(ri,impul,resi,3)      

     ! resj = invIj � (rj x impul)   (for jimbo)
     call vecpro(rj,impul,resj,3) 
     do idime=1,3     
        ! Multiply by the inertia, in 2D is a scalar 
        resi(idime) = (1.0_rp/ITi(1)) * resi(idime)
        resj(idime) = (1.0_rp/ITj(1)) * resj(idime)
     end do
  end if
  !
  ! Angular energy in 2D: Er = (1/2) * I * w^2
  ! Linear energy:        El = (1/2) * m * v^2
  ! Toral  energy:        Et = Er + Et
  ! Fuente: http://en.wikipedia.org/wiki/Rotational_energy
  !         http://en.wikipedia.org/wiki/Kinetic_energy#Rotation_in_systems
  !
  !v2i = 0.0_rp;  w2i = 0.0_rp;   v2j = 0.0_rp;  w2j = 0.0_rp
  !do idime = 1,3
  !   v2i = v2i + vi(idime,1) * vi(idime,1)
  !   w2i = w2i + wi(idime,1) * wi(idime,1)
  !   v2j = v2j + vj(idime,1) * vj(idime,1)
  !   w2j = w2j + wj(idime,1) * wj(idime,1)
  !end do
  !Et = (0.5_rp)*imbou(iimbo)%massa*v2i + (0.5_rp)*imbou(jimbo)%massa*v2j + (0.5_rp)*ITi(1)*w2i + (0.5_rp)*ITj(1)*w2j

  do idime = 1,3
     vi(idime,1) = vi(idime,1) + impul(idime) * (1.0_rp/imbou(iimbo)%massa)
     wi(idime,1) = wi(idime,1) + resi(idime)
     vj(idime,1) = vj(idime,1) - impul(idime) * (1.0_rp/imbou(jimbo)%massa)
     wj(idime,1) = wj(idime,1) - resj(idime)
  end do
  !v2i = 0.0_rp;  w2i = 0.0_rp;   v2j = 0.0_rp;  w2j = 0.0_rp
  !do idime = 1,3
  !   v2i = v2i + vi(idime,1) * vi(idime,1)
  !   w2i = w2i + wi(idime,1) * wi(idime,1)
  !   v2j = v2j + vj(idime,1) * vj(idime,1)
  !   w2j = w2j + wj(idime,1) * wj(idime,1)
  !end do
  !Et = (0.5_rp)*imbou(iimbo)%massa*v2i + (0.5_rp)*imbou(jimbo)%massa*v2j + (0.5_rp)*ITi(1)*w2i + (0.5_rp)*ITj(1)*w2j
end subroutine ibm_ipapar


subroutine ibm_ipawal(iimbo,coori,norma)

  !-----------------------------------------------------------------------
  ! NAME
  !    ibm_ipawal
  ! DESCRIPTION
  !    This routines calculate the impulse of the collision and the velocities 
  !    between a particles and a wall
  ! USED BY
  !    
  !-----------------------------------------------------------------------

  use def_kintyp, only     :  ip,rp
  use def_master, only     :  imbou
  use def_domain, only     :  ndime,nimbo,nnode
  use def_immbou, only     :  kfl_rotib_ibm,dtime_ibm
  implicit none
  integer(ip), intent(in)    :: iimbo
  integer(ip)                :: idime

  real(rp),    intent(in)    :: coori(3),norma(3)
  real(rp)                   :: impul(3)
  real(rp)                   :: Io(3,3),Item(3,3),I(3,3),invIi(3,3)
  real(rp)                   :: ri(3),vti(3),vtem(3),resi(3)
  real(rp)                   :: dt
  real(rp)                   :: C,numer,denom,immag,deter
  real(rp),    pointer       :: xi(:,:),ai(:,:),vi(:,:),wi(:,:),ITi(:),Roti(:,:) ! Linear and angular quantities of iimbo
  !
  ! Initialization
  ! Restitution coefficient. When C=1.0, perfect elasticity, if C=0.0, perfect inelasticity
  !
  C     = 0.8_rp  
  dt    = dtime_ibm
  do idime = 1,3
     vti(idime) = 0.0_rp
  end do
  xi   => imbou(iimbo)%posil
  vi   => imbou(iimbo)%velol
  ai   => imbou(iimbo)%accel
  wi   => imbou(iimbo)%veloa
  ITi  => imbou(iimbo)%momin
  Roti => imbou(iimbo)%rotac
  !
  ! Coordinates of contact point respect to the center of mass of the particles 
  !
  ri(3) = 0.0_rp
  resi(3) = 0.0_rp
  do idime=1,ndime
     ri(idime) = coori(idime) - imbou(iimbo)%posil(idime,1)
     resi(idime) = 0.0_rp
  end do
  !------------------------------------------------------------------
  !
  ! Calculate the impulse of the collision
  !
  !------------------------------------------------------------------
  if( kfl_rotib_ibm /= 0 ) then
     if (ndime == 3) then

        ! Compute the actual inertia tensor for iimbo
        Io(1,1)=ITi(1); Io(1,2)=ITi(4); Io(1,3)=ITi(5)
        Io(2,1)=ITi(4); Io(2,2)=ITi(2); Io(2,3)=ITi(6)
        Io(3,1)=ITi(5); Io(3,2)=ITi(6); Io(3,3)=ITi(3)        
        call mbmab0(Item,Roti,Io,3,3,3)
        call mbmabt(I,Item,Roti,3,3,3)        
        ! Compute the actual inverse inertia tensor for iimbo
        call invmtx(I,invIi,deter,3)
        !
        ! resi = (invIj � (rj x norma)) x rj   (for jimbo)
        !
        call vecpro(ri,norma,resi,3)
        call mbvab0(vtem,invIi,resi,3,3) 
        call vecpro(vtem,ri,resi,3)                                 
     elseif (ndime == 2) then        
        ! resi = (invIj � (rj x norma)) x rj   (for jimbo)
        call vecpro(ri,norma,vtem,3)      
        call vecpro(vtem,ri,resi,3)                 
        do idime=1,3     
           ! Multiply by the inertia, in 2D is a scalar 
           resi(idime) = (1.0_rp/ITi(1)) * resi(idime)
        end do
     end if
     call vecpro(wi,ri,vti,3)
  end if
  numer = 0.0_rp
  denom = 0.0_rp
  immag = 0.0_rp
  !
  ! Calculate the total velocity before the collision
  !  
  do idime=1,3     
     vti(idime) = vi(idime,1) + vti(idime)
     denom = denom + norma(idime) * resi(idime)
     numer = numer + norma(idime) * vti(idime)
  end do
  denom = denom + (1.0_rp/imbou(iimbo)%massa)
  numer = (-1.0_rp-C)*numer
  !
  ! Magnitude of the impulse
  !
  immag = numer * (1.0_rp/denom)
  !
  ! Impulse vector
  !
  do idime=1,3
     impul(idime) = immag * norma(idime)
  end do
  !------------------------------------------------------------------
  !
  ! Calculate the total velocity before the collision
  !
  !------------------------------------------------------------------
  if (ndime == 3) then
     ! resi = invIi � (ri x impul)  (for iimbo)
     call vecpro(ri,impul,vtem,3)
     call mbvab0(resi,invIi,vtem,3,3)

  elseif (ndime == 2) then
     ! resi = invIi � (ri x impul)  (for iimbo)
     call vecpro(ri,impul,resi,3)
     do idime=1,3
        ! Multiply by the inertia, in 2D is a scalar 
        resi(idime) = (1.0_rp/ITi(1)) * resi(idime)
     end do
  end if
  do idime = 1,3
     vi(idime,1) = vi(idime,1) + impul(idime) * (1.0_rp/imbou(iimbo)%massa)
     wi(idime,1) = wi(idime,1) + resi(idime)
  end do

end subroutine ibm_ipawal
