module mod_dod_extens

  use def_kintyp, only : ip,rp,lg
  use def_parame, only : pi
  use def_elmtyp, only : BOEXT
  use mod_memory

contains

  subroutine dod_extens_crea2d(&       
       ipoin,boun1,boun2,mnodb,lnodb,pbopo_dod,lbopo_dod,lpobo_dod,coord,&
       nnext,lnext,lpoex_bo1,lpoex_bo2,ntext,ipex1,ipex2)
    implicit none  
    integer(ip), intent(in)    :: ipoin                        !< Fringe node from which we extend
    integer(ip), intent(in)    :: boun1                        !< Boundary 1
    integer(ip), intent(in)    :: boun2                        !< Boundary 2
    integer(ip), intent(in)    :: mnodb                        !< Boundary 2
    integer(ip), intent(in)    :: lnodb(mnodb,*)               !< Boundary connectivity
    integer(ip), intent(in)    :: pbopo_dod(*)                 !< Node-boundary graph
    integer(ip), intent(in)    :: lbopo_dod(*)                 !< Node-boundary graph
    integer(ip), intent(in)    :: lpobo_dod(*)                 !< Boundary node type (-1: fringe, 1: real boundary)
    real(rp),    intent(in)    :: coord(2,*)                   !< Node coordinates
    integer(ip), intent(inout) :: nnext                        !< Number of candidates
    integer(ip), intent(inout) :: lnext(nnext)                 !< Original list of candidates
    integer(ip), intent(inout) :: lpoex_bo1(nnext)             !< Candidate of boundary boun1
    integer(ip), intent(inout) :: lpoex_bo2(nnext)             !< Candidate of boundary boun2
    integer(ip), intent(inout) :: ntext                        !< Number of nodes (1 or 2)
    integer(ip), intent(inout) :: ipex1                        !< Extension node 1
    integer(ip), intent(inout) :: ipex2                        !< Extension node 2

    integer(ip)                :: kfl_apert,tmp,perm
    integer(ip)                :: kfl_esqui,kfl_boxxx
    integer(ip)                :: boun3,boun4,cont
    integer(ip)                :: inod1,inod2,kk,ibopo,idime
    integer(ip)                :: inext,jnext,inexb1,inexb2
    integer(ip)                :: pexti,pextj,pextk,jpoin
    integer(ip)                :: kpoin,bound
    real(rp)                   :: qaux,q,q_old,q_res,q1,q2
    real(rp)                   :: norma,norme,cosa,cosb1
    real(rp)                   :: cosb2,nord1,nord2
    real(rp)                   :: norb1(2),norb2(2),normp(2)
    real(rp)                   :: vecex(2),vecex_bo1(2),vecex_bo2(2)
    real(rp)                   :: comp1(2),tang1(2),tang2(2),dirc3(2),dirc4(2)
    real(rp)                   :: coseno_aux1,coseno_aux2,coseno_aux3
    real(rp)                   :: coseno_aux4,coseno_aux5,coseno_aux6
    real(rp)                   :: cosen1_nei,cosen2_nei
    real(rp)                   :: angle_max,cos_angle_max,typical_distance
    real(rp)                   :: toler_distance
    real(rp)                   :: coord_ipoin(2)
    real(rp)                   :: coord_ipex1(2)
    real(rp)                   :: coord_ipex2(2)
    real(rp)                   :: coord_ipex3(2)
    integer(ip)                :: cont1,cont2,inodb,pnodb
    integer(ip), pointer       :: indic(:)
    integer(ip), pointer       :: lnext_tmp(:)
    real(rp),    pointer       :: quali_bo1(:)
    real(rp),    pointer       :: quali_bo2(:)
    !
    ! Corners detection
    !
    kfl_esqui = 0
    angle_max = 89.0_rp * 1.74532925199e-2_rp ! Conversion to radians
    cos_angle_max = cos(angle_max)            ! Max cos to detect a corner
    !
    ! Typical distance (used to normalize things)
    !
    typical_distance =  sqrt( (coord(1,lnodb(1,boun1))-coord(1,lnodb(2,boun1)))**2 &
         &                  + (coord(2,lnodb(1,boun1))-coord(2,lnodb(2,boun1)))**2 )
    toler_distance   =  1.0e-03_rp * typical_distance
    !
    ! Save IPOIN coordinates
    !
    do idime = 1,2
       coord_ipoin(idime) = coord(idime,ipoin)
    end do
    !
    ! Allocate temporary memory
    !
    allocate( lnext_tmp(2*nnext) )
    allocate( indic(nnext)       )
    allocate( quali_bo1(nnext)   )
    allocate( quali_bo2(nnext)   )
    do inext = 1,nnext
       indic(inext) = 0
    end do

    kfl_boxxx = 0
    pnodb     = 2_ip
    ipex1     = 0_ip
    ipex2     = 0_ip
    !
    ! Check number of extension boundaries we have
    ! BOUND = 0 ....... both boundaries
    !       = BOUN1 ... boundary 1
    !       = BOUN2 ... boundary 2
    !
    bound = 0
    cont1 = 0
    cont2 = 0
    do inodb = 1,pnodb
       if( lpobo_dod(lnodb(inodb,boun1)) == -1 ) cont1 = cont2 + 1  
       if( lpobo_dod(lnodb(inodb,boun2)) == -1 ) cont2 = cont2 + 1  
    end do
    if( cont1 + cont2 == 2*pnodb ) then
       bound = 0
    else if( cont1 == pnodb ) then
       bound = boun1
    else if( cont2 == pnodb ) then
       bound = boun2
    else
       write(*,*) 'a=',lnodb(1,boun1),lnodb(2,boun1)
       write(*,*) 'b=',lpobo_dod(lnodb(1,boun1)),lpobo_dod(lnodb(2,boun1))
       write(*,*) 'c=',lnodb(1,boun2),lnodb(2,boun2)
       write(*,*) 'd=',lpobo_dod(lnodb(1,boun2)),lpobo_dod(lnodb(2,boun2))
       call runend('DOD_EXTENS: WRONG FRINGE NODE')
    end if

    if( bound == 0 ) then

       !-----------------------------------------------------------------     
       !
       ! I have two boundaries to extend from 
       !
       ! Normalized norms: NORB1(2), NORB2(2) of NBOU1 and NBOU2
       !
       !-----------------------------------------------------------------

       norb1(1) = -( coord(2,lnodb(1,boun1))-coord(2,lnodb(2,boun1)) ) ! n1(x)
       norb1(2) =  ( coord(1,lnodb(1,boun1))-coord(1,lnodb(2,boun1)) ) ! n1(y)
       norb2(1) = -( coord(2,lnodb(1,boun2))-coord(2,lnodb(2,boun2)) ) ! n2(x)
       norb2(2) =  ( coord(1,lnodb(1,boun2))-coord(1,lnodb(2,boun2)) ) ! n2(y)
       norme    =  1.0_rp / sqrt(norb1(1)*norb1(1)+norb1(2)*norb1(2))
       norb1(1) =  norb1(1) * norme
       norb1(2) =  norb1(2) * norme
       norme    =  1.0_rp / sqrt(norb2(1)*norb2(1)+norb2(2)*norb2(2))
       norb2(1) =  norb2(1) * norme
       norb2(2) =  norb2(2) * norme
       !
       ! Check open/close corners
       !
       kfl_apert = -1_ip 

       do inodb = 1,pnodb
          !
          ! Normalized tangent vectors TANG1(2), TANG2(2) of NBOU1 and NBOU2
          !
          if( lnodb(inodb,boun1) /= ipoin ) then
             inod1     = inodb
             tang1(1)  = coord(1,lnodb(inodb,boun1)) - coord(1,ipoin)
             tang1(2)  = coord(2,lnodb(inodb,boun1)) - coord(2,ipoin)
             nord1     = 1.0_rp / sqrt(tang1(1) * tang1(1) + tang1(2) * tang1(2))
             tang1(1)  = tang1(1) * nord1 
             tang1(2)  = tang1(2) * nord1
          end if

          if( lnodb(inodb,boun2) /= ipoin ) then
             inod2     = inodb
             tang2(1)  = coord(1,lnodb(inodb,boun2)) - coord(1,ipoin)
             tang2(2)  = coord(2,lnodb(inodb,boun2)) - coord(2,ipoin)
             nord2     = 1.0_rp / sqrt(tang2(1) * tang2(1) + tang2(2) * tang2(2))
             tang2(1)  = tang2(1) * nord2 
             tang2(2)  = tang2(2) * nord2
          end if
       end do
       !              TANG2
       !                _      x
       !     TANG1      /|    /
       !      <--      /    /  BOUN2
       !                  /
       !   x-------------o
       !     BOUN1        IPOIN
       !
       coseno_aux3 = tang1(1) * tang2(1) + tang1(2) * tang2(2)
       !
       ! Check if we have a right angle => KFL_ESQUI = 1
       !
       if( abs(coseno_aux3) < cos_angle_max ) kfl_esqui = 1
       !
       ! Average normal NORMP = 1/2 ( NORB1 + NORB2 ) 
       !
       normp(1) = norb1(1) + norb2(1)
       normp(2) = norb1(2) + norb2(2)
       norma    = 1.0_rp / sqrt( normp(1) * normp(1) + normp(2)*normp(2) )
       normp(1) = normp(1) * norma
       normp(2) = normp(2) * norma
       !
       ! COSENO_AUX3, COSENO_AUX4: scalar products between averaged normal and tangents
       !
       coseno_aux3 = tang1(1) * normp(1) + tang1(2) * normp(2)
       coseno_aux4 = tang2(1) * normp(1) + tang2(2) * normp(2)
       !
       ! Case KFL_ESQUI = 1: what kind of corner?
       !
       !  KFL_APERT = 0         KFL_APERT = 1
       !      x-------
       !      |
       !      |                   x-----o
       ! x----o                         |   
       ! |                  subdomain   |
       ! |    subdomain                 x
       !
       if( kfl_esqui == 1 ) then
          if( coseno_aux3 > 0.0_rp ) then
             kfl_apert = 0                      ! se ven los ojos
          else 
             kfl_apert = 1                      ! se ven el culo
          end if
       end if
       !
       ! Loop over candidates: VECEX = VEC(IPOIN,CANDIDATE)
       !
       inexb1 = 0
       inexb2 = 0
       do inext = 1,nnext
          pexti    = lnext(inext)          
          vecex(1) = coord(1,pexti) - coord(1,ipoin)
          vecex(2) = coord(2,pexti) - coord(2,ipoin)
          norme    = sqrt( vecex(1) * vecex(1) + vecex(2) * vecex(2) )

          if( norme < toler_distance ) then
             !
             ! Candidate is on me! Remove from list
             !
             cosb1 = 0.0_rp
             cosb2 = 0.0_rp
          else
             !
             ! Compute angles between Tangents and VEC(IPOIN,CANDIDATE)
             !
             norme       = 1.0_rp / norme
             vecex(1)    = vecex(1) * norme
             vecex(2)    = vecex(2) * norme
             coseno_aux1 = vecex(1) * tang1(1) + vecex(2) * tang1(2)
             coseno_aux2 = vecex(1) * tang2(1) + vecex(2) * tang2(2)
             cosa        = vecex(1) * normp(1) + vecex(2) * normp(2)
             cosb1       = vecex(1) * norb1(1) + vecex(2) * norb1(2)
             cosb2       = vecex(1) * norb2(1) + vecex(2) * norb2(2)
             !
             ! Take only candidates located between -60 and 60 degrees: (x)
             ! Do not consider others (-)
             !
             !          x
             !       x         -
             !     \     /    
             ! -    \   /  -      -
             !   -   \ /
             ! o------o------o
             !
             if( cosb1 > 0.4_rp ) then
                inexb1 = inexb1 + 1            ! For BOUN1
                lpoex_bo1(inexb1) = pexti
             end if
             if( cosb2 > 0.4_rp ) then
                inexb2 = inexb2 + 1            ! For BOUN2
                lpoex_bo2(inexb2) = pexti
             end if
          end if
       end do
       !
       ! Look for neighboring boundary BOUN3 to BOUN1
       !
       !               jpoin    BOUN3
       !                   x-----------x
       !                   |
       !                   | BOUN1
       !                   |
       !        x----------o   
       !        |  BOUN2    ipoin
       !  BOUN4 |    
       !        |       subdomain
       !        x
       !      kpoin
       !
       jpoin = lnodb(1,boun1)
       do ibopo = pbopo_dod(jpoin),pbopo_dod(jpoin+1)-1
          cont = 0
          do inodb = 1,pnodb
             if( lnodb(inodb,lbopo_dod(ibopo)) /= ipoin ) then 
                cont = cont + 1
             end if
          end do
          if( cont == pnodb ) boun3 = lbopo_dod(ibopo)
       end do
       !
       ! DIRC3 = VEC(BOUN3)
       !
       dirc3(1) = coord(1,lnodb(1,boun3)) - coord(1,lnodb(2,boun3))
       dirc3(2) = coord(2,lnodb(1,boun3)) - coord(2,lnodb(2,boun3))
       norme    = 1.0_rp / sqrt( dirc3(1) * dirc3(1) + dirc3(2) * dirc3(2))
       dirc3(1) = dirc3(1) * norme 
       dirc3(2) = dirc3(2) * norme

       kpoin = lnodb(2,boun2)
       do ibopo = pbopo_dod(kpoin),pbopo_dod(kpoin+1)-1
          cont = 0
          do inodb = 1,pnodb
             if( lnodb(inodb,lbopo_dod(ibopo)) /= ipoin ) then !no es este
                cont = cont + 1
             end if
          end do
          if( cont == pnodb ) boun4 = lbopo_dod(ibopo)
       end do
       !
       ! DIRC4 = VEC(BOUN4)
       !
       dirc4(1)  = coord(1,lnodb(2,boun4)) - coord(1,lnodb(1,boun4))
       dirc4(2)  = coord(2,lnodb(2,boun4)) - coord(2,lnodb(1,boun4))
       norme     = 1.0_rp / sqrt( dirc4(1) * dirc4(1) + dirc4(2) * dirc4(2) )
       dirc4(1)  = dirc4(1) * norme 
       dirc4(2)  = dirc4(2) * norme
       !
       ! Angles between TANG1.DIRC3 and TANG2.DIRC4
       !
       cosen1_nei = tang1(1) * dirc3(1) + tang1(2) * dirc3(2)
       cosen2_nei = tang2(1) * dirc4(1) + tang2(2) * dirc4(2)

       if( abs(cosen1_nei) < 1.0e-7_rp .or. abs(cosen2_nei) < 1.0e-7_rp ) kfl_boxxx = 1
       !
       ! Candidates for BOUN1
       !
       do inext = 1,inexb1
          pexti = lpoex_bo1(inext)
          pextj = lnodb(inod1,boun1)
          do idime = 1,2
             coord_ipex1(idime) = coord(idime,pextj)
             coord_ipex2(idime) = coord(idime,pexti)
          end do
          call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pextj,pexti,q,perm)
          quali_bo1(inext) = q
       end do
       !
       ! Order quality for BOUN1 triangles
       !
       do inext = 1,inexb1-1
          do jnext = inext+1,inexb1
             if( quali_bo1(inext) > quali_bo1(jnext) ) then
                tmp              = lpoex_bo1(inext)
                lpoex_bo1(inext) = lpoex_bo1(jnext)
                lpoex_bo1(jnext) = tmp
                qaux             = quali_bo1(inext)
                quali_bo1(inext) = quali_bo1(jnext)
                quali_bo1(jnext) = qaux
             end if
          end do
       end do
       !
       ! Candidates for BOUN2
       !       
       do inext = 1,inexb2
          pexti = lpoex_bo2(inext)
          pextj = lnodb(inod2,boun2)
          coord_ipex1(1) = coord(1,pextj)
          coord_ipex1(2) = coord(2,pextj)
          coord_ipex2(1) = coord(1,pexti)
          coord_ipex2(2) = coord(2,pexti)
          call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pextj,pexti,q,perm)
          quali_bo2(inext) = q
       end do
       !
       ! Order quality for BOUN2 triangles
       !
       do inext = 1,inexb2-1
          do jnext = inext+1,inexb2
             if( quali_bo2(inext) > quali_bo2(jnext) ) then
                tmp              = lpoex_bo2(inext)
                lpoex_bo2(inext) = lpoex_bo2(jnext)
                lpoex_bo2(jnext) = tmp
                qaux             = quali_bo2(inext)
                quali_bo2(inext) = quali_bo2(jnext)
                quali_bo2(jnext) = qaux
             end if
          end do
       end do
       !
       ! Choose best combination
       !
       q_old = 1000.0_rp
       cont  = 0
       do inext = 1,inexb1
          do jnext = 1,inexb2

             if( lpoex_bo1(inext) /= lpoex_bo2(jnext) ) then
                !
                ! Both extensions have chosen a different node
                !
                q1             = 0.0_rp
                q2             = 0.0_rp
                pexti          = lpoex_bo1(inext)
                pextj          = lpoex_bo2(jnext)
                coord_ipex1(1) = coord(1,pexti)
                coord_ipex2(1) = coord(1,pextj)
                coord_ipex1(2) = coord(2,pexti)
                coord_ipex2(2) = coord(2,pextj)

                vecex_bo1(1)   = coord_ipex1(1)-coord_ipoin(1)
                vecex_bo1(2)   = coord_ipex1(2)-coord_ipoin(2)
                norme          = 1.0_rp / sqrt( vecex_bo1(1) * vecex_bo1(1) + vecex_bo1(2) * vecex_bo1(2) )
                vecex_bo1(1)   = vecex_bo1(1) * norme
                vecex_bo1(2)   = vecex_bo1(2) * norme

                vecex_bo2(1)   = coord_ipex2(1)-coord_ipoin(1)
                vecex_bo2(2)   = coord_ipex2(2)-coord_ipoin(2)
                norme          = 1.0_rp / sqrt( vecex_bo2(1) * vecex_bo2(1) + vecex_bo2(2) * vecex_bo2(2) )
                vecex_bo2(1)   = vecex_bo2(1) * norme
                vecex_bo2(2)   = vecex_bo2(2) * norme

                coseno_aux4    = vecex_bo1(1) * tang1(1) + vecex_bo1(2) * tang1(2)
                coseno_aux5    = vecex_bo2(1) * tang2(1) + vecex_bo2(2) * tang2(2)
                coseno_aux6    = vecex_bo2(1) * tang1(1) + vecex_bo2(2) * tang1(2)

                if( coseno_aux6 > coseno_aux4 ) q2 = 2001.0_rp

                if( kfl_esqui == 1 .and. kfl_apert == 0 ) then             
                   if( coseno_aux4 < 0.0_rp ) q1 = 2001.0_rp                
                   if( coseno_aux5 < 0.0_rp ) q2 = 2001.0_rp
                   if( coseno_aux4 > 0.0_rp .and. coseno_aux5 > 0.0_rp .and. coseno_aux6 > coseno_aux4 ) q2 = 2001.0_rp
                end if
                if( kfl_boxxx == 1 .and. kfl_apert /= 1 ) then 
                   if( coseno_aux4 < 0.5_rp ) q1 = 2001.0_rp
                   if( coseno_aux5 < 0.5_rp ) q2 = 2001.0_rp
                end if

                call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pexti,pextj,q,perm)
                q_res = max(quali_bo1(inext),quali_bo2(jnext),q,q1,q2)

             else

                q_res = max(quali_bo1(inext),quali_bo2(jnext))

             end if

             if( q_res <= q_old ) then
                cont  = cont +  1
                q_old = q_res
                if( lpoex_bo1(inext) /= lpoex_bo2(jnext) ) then
                   ipex1 = lpoex_bo1(inext)
                   ipex2 = lpoex_bo2(jnext)
                   ntext = 3
                else
                   ntext = 2
                   ipex1 = lpoex_bo1(inext)
                   ipex2 = 0_ip
                end if
             end if

          end do
       end do
       !
       ! No good extension elements have been found
       !
       if( cont == 0 ) then

          do inext = 1, inexb1
             lnext_tmp(inext) = lpoex_bo1(inext)
          end do
          kk = inexb1
          do inext = 1,inexb2
             kk = kk + 1
             lnext_tmp(kk) = lpoex_bo2(inext)
          end do

          q_old = 1000.0_rp          

          do inext = 1,kk
             pexti          = lnext_tmp(inext)
             pextj          = lnodb(1,boun1)
             pextk          = lnodb(2,boun2)
             coord_ipex1(1) = coord(1,pexti)
             coord_ipex1(2) = coord(2,pexti)
             coord_ipex2(1) = coord(1,pextj)
             coord_ipex2(2) = coord(2,pextj)
             coord_ipex3(1) = coord(1,pextk)
             coord_ipex3(2) = coord(2,pextk)

             vecex_bo1(1)   = coord_ipex1(1)-coord_ipoin(1)
             vecex_bo1(2)   = coord_ipex1(2)-coord_ipoin(2)
             norme          = 1.0_rp / sqrt( vecex_bo1(1) * vecex_bo1(1) + vecex_bo1(2) * vecex_bo1(2) )
             vecex_bo1(1)   = vecex_bo1(1) * norme
             vecex_bo1(2)   = vecex_bo1(2) * norme

             call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pexti,pextj,q1,perm)
             call dod_extens_qual2d(coord_ipoin,coord_ipex3,coord_ipex1,pextk,pexti,q2,perm)
             q_res = max(q1,q2)

             if( kfl_esqui == 1 .and. kfl_apert == 0 ) then
                coseno_aux4 = vecex_bo1(1) * tang1(1) + vecex_bo1(2) * tang1(2)
                coseno_aux5 = vecex_bo1(1) * tang2(1) + vecex_bo1(2) * tang2(2)
             end if
             if( kfl_boxxx == 1 .and. kfl_apert /= 1 ) then  
                if( coseno_aux4 < 0.5_rp ) quali_bo1(inext) = 2001.0_rp
                if( coseno_aux5 < 0.5_rp ) quali_bo2(jnext) = 2001.0_rp
             end if

             if( coseno_aux4 < 0.0_rp .or. coseno_aux5 < 0.0_rp ) q_res = 2001.0_rp                

             if( q_res <= q_old ) then
                cont  = cont + 1
                q_old = q_res
                ntext = 2
                ipex1 = lnext_tmp(inext)
                ipex2 = 0                
             end if
          end do

       end if

       if( cont == 0 ) then
          print*,ipoin,lnext(1:nnext)     
          call runend('DOD_NORM2D: I RUN OUT OF IDEA TO CREATE EXTENSIONS...')
       end if

    else

       !-----------------------------------------------------------------     
       !
       ! I have only one boundary to extend from: BOUND
       !
       !-----------------------------------------------------------------

       comp1(1) =  coord(1,lnodb(1,bound))-coord(1,lnodb(2,bound))
       comp1(2) =  coord(2,lnodb(1,bound))-coord(2,lnodb(2,bound))
       norb1(1) = -comp1(2)
       norb1(2) =  comp1(1)
       do inodb = 1,pnodb
          if( lnodb(inodb,bound) /= ipoin ) then
             inod1     = inodb
             tang1(1)  = coord(1,lnodb(inodb,bound)) - coord(1,ipoin)
             tang1(2)  = coord(2,lnodb(inodb,bound)) - coord(2,ipoin)
             nord1     = 1.0_rp / sqrt( tang1(1) * tang1(1) + tang1(2) * tang1(2) )
             tang1(1)  = tang1(1) * nord1 
             tang1(2)  = tang1(2) * nord1
          end if
       end do
       inexb1 = 0
       inexb2 = 0
       do inext = 1,nnext

          pexti    = lnext(inext)
          vecex(1) = coord(1,pexti)-coord_ipoin(1)
          vecex(2) = coord(2,pexti)-coord_ipoin(2)
          norme    = sqrt( vecex(1) * vecex(1) + vecex(2) * vecex(2) )

          if( norme < toler_distance ) then
             cosb1 = 0.0_rp
             cosb2 = 0.0_rp
          else
             vecex(1)    = vecex(1) / norme
             vecex(2)    = vecex(2) / norme
             coseno_aux1 = vecex(1) * tang1(1) + vecex(2) * tang1(2)

             cosb1 = vecex(1)*norb1(1) + vecex(2)*norb1(2)
             if( cosb1 < 0.4_rp ) then
                inexb1 = inexb1 + 1
                lpoex_bo1(inexb1) = pexti
             end if
          end if
       end do

       do inext = 1,inexb1
          pexti = lpoex_bo1(inext)
          pextj = lnodb(inod1,bound)
          coord_ipex1(1) = coord(1,pextj)
          coord_ipex1(2) = coord(2,pextj)
          coord_ipex2(1) = coord(1,pexti)
          coord_ipex2(2) = coord(2,pexti)
          call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pextj,pexti,q,perm)
          quali_bo1(inext) = q
       end do

       do inext = 1,inexb1-1
          do jnext = inext+1,inexb1
             if( quali_bo1(inext) > quali_bo1(jnext) ) then
                tmp              = lpoex_bo1(inext)
                lpoex_bo1(inext) = lpoex_bo1(jnext)
                lpoex_bo1(jnext) = tmp
                qaux             = quali_bo1(inext)
                quali_bo1(inext) = quali_bo1(jnext)
                quali_bo1(jnext) = qaux
             end if
          end do
       end do

       ipex1 = lpoex_bo1(1)
       ipex2 = 0

    end if

    deallocate( quali_bo2 )
    deallocate( quali_bo1 )
    deallocate( indic     )
    deallocate( lnext_tmp )

  end subroutine dod_extens_crea2d

  subroutine dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,ipex1,ipex2,quali,iperm)
    implicit none
    real(rp),    intent(in)     :: coord_ipoin(2)
    real(rp),    intent(in)     :: coord_ipex1(2)
    real(rp),    intent(in)     :: coord_ipex2(2)
    integer(ip), intent(inout)  :: ipex1
    integer(ip), intent(inout)  :: ipex2
    integer(ip), intent(inout)  :: iperm
    real(rp),    intent(out)    :: quali
    integer(ip)                 :: ipex1_tmp,idime
    real(rp)                    :: perim,side1,side2,side3,hmaxi,alpha,volum
    real(rp)                    :: coora(2),coorb(2)

    iperm = 0
    alpha = 0.28867513459_rp ! sqrt(3)/6
    side1 = 0
    side2 = 0
    side3 = 0

    do idime = 1,2
       side1 = side1 + ( coord_ipex1(idime) - coord_ipoin(idime) ) * ( coord_ipex1(idime) - coord_ipoin(idime) ) 
       side2 = side2 + ( coord_ipex2(idime) - coord_ipoin(idime) ) * ( coord_ipex2(idime) - coord_ipoin(idime) ) 
       side3 = side3 + ( coord_ipex2(idime) - coord_ipex1(idime) ) * ( coord_ipex2(idime) - coord_ipex1(idime) ) 
    end do
    side1 = sqrt(side1)
    side2 = sqrt(side2)
    side3 = sqrt(side3)
    hmaxi = max(side1,side2,side3)
    perim = 0.5_rp * ( side1 + side2 + side3 )

    coora(1) = coord_ipex1(1) - coord_ipoin(1)  
    coora(2) = coord_ipex1(2) - coord_ipoin(2)  
    coorb(1) = coord_ipex2(1) - coord_ipoin(1)  
    coorb(2) = coord_ipex2(2) - coord_ipoin(2)  

    volum = 0.5_rp * ( coora(1) * coorb(2) - coorb(1) * coora(2) )
    !
    ! Quality
    !
    if ( volum == 0.0_rp ) then

       quali = 100000_rp

    else if ( volum < 0.0_rp ) then

       ipex1_tmp = ipex1
       ipex1     = ipex2
       ipex2     = ipex1_tmp
       iperm      = 1  

    end if

    if( volum /= 0.0_rp ) quali = alpha * hmaxi * perim / abs(volum)

  end subroutine dod_extens_qual2d

  subroutine dod_extension_elements_2d(&       
       ipoin,boun1,boun2,mnodb,lnodb,lboch,pbopo,lbopo,coord,&
       nnext,lnext,lpoex_bo1,lpoex_bo2,ntext,ipex1,ipex2)
    implicit none  
    integer(ip), intent(in)    :: ipoin                        !< Fringe node from which we extend
    integer(ip), intent(in)    :: boun1                        !< Boundary 1
    integer(ip), intent(in)    :: boun2                        !< Boundary 2
    integer(ip), intent(in)    :: mnodb                        !< Boundary 2
    integer(ip), intent(in)    :: lnodb(mnodb,*)               !< Boundary connectivity
    integer(ip), intent(in)    :: lboch(*)                     !< Boundary type
    integer(ip), intent(in)    :: pbopo(*)                     !< Node-boundary graph
    integer(ip), intent(in)    :: lbopo(*)                     !< Node-boundary graph
    real(rp),    intent(in)    :: coord(2,*)                   !< Node coordinates
    integer(ip), intent(inout) :: nnext                        !< Number of candidates
    integer(ip), intent(inout) :: lnext(nnext)                 !< Original list of candidates
    integer(ip), intent(inout) :: lpoex_bo1(nnext)             !< Candidate of boundary boun1
    integer(ip), intent(inout) :: lpoex_bo2(nnext)             !< Candidate of boundary boun2
    integer(ip), intent(inout) :: ntext                        !< Number of nodes (1 or 2)
    integer(ip), intent(inout) :: ipex1                        !< Extension node 1
    integer(ip), intent(inout) :: ipex2                        !< Extension node 2

    integer(ip)                :: kfl_apert,tmp,perm
    integer(ip)                :: kfl_esqui,kfl_boxxx
    integer(ip)                :: boun3,boun4,cont
    integer(ip)                :: inod1,inod2,kk,ibopo,idime
    integer(ip)                :: inext,jnext,inexb1,inexb2
    integer(ip)                :: pexti,pextj,pextk,jpoin
    integer(ip)                :: kpoin,bound
    integer(ip)                :: ipoi1_boun1,ipoi2_boun1
    integer(ip)                :: ipoi1_boun2,ipoi2_boun2 
    real(rp)                   :: qaux,q,q_old,q_res,q1,q2
    real(rp)                   :: norma,norme,cosa,cosb1
    real(rp)                   :: cosb2,nord1,nord2
    real(rp)                   :: norb1(2),norb2(2),normp(2)
    real(rp)                   :: vecex(2),vecex_bo1(2),vecex_bo2(2)
    real(rp)                   :: comp1(2),tang1(2),tang2(2),dirc3(2),dirc4(2)
    real(rp)                   :: coseno_aux1,coseno_aux2,coseno_aux3
    real(rp)                   :: coseno_aux4,coseno_aux5,coseno_aux6
    real(rp)                   :: cosen1_nei,cosen2_nei
    real(rp)                   :: angle_max,cos_angle_max,typical_distance
    real(rp)                   :: toler_distance
    real(rp)                   :: coord_ipoin(2)
    real(rp)                   :: coord_ipex1(2)
    real(rp)                   :: coord_ipex2(2)
    real(rp)                   :: coord_ipex3(2)
    integer(ip)                :: cont1,cont2,inodb,pnodb
    integer(ip), pointer       :: indic(:)
    integer(ip), pointer       :: lnext_tmp(:)
    real(rp),    pointer       :: quali_bo1(:)
    real(rp),    pointer       :: quali_bo2(:)
    !
    ! Corners detection
    ! 
    kfl_esqui = 0
    angle_max = 89.0_rp * 1.74532925199e-2_rp ! Conversion to radians
    cos_angle_max = cos(angle_max)            ! Max cos to detect a corner
    !
    ! Typical distance (used to normalize things)
    !
    typical_distance =  sqrt( (coord(1,lnodb(1,boun1))-coord(1,lnodb(2,boun1)))**2 &
         &                  + (coord(2,lnodb(1,boun1))-coord(2,lnodb(2,boun1)))**2 )
    toler_distance   =  1.0e-03_rp * typical_distance
    !
    ! Save IPOIN coordinates
    !
    do idime = 1,2
       coord_ipoin(idime) = coord(idime,ipoin)
    end do
    !
    ! Allocate temporary memory
    !
    allocate( lnext_tmp(2*nnext) )
    allocate( indic(nnext)       )
    allocate( quali_bo1(nnext)   )
    allocate( quali_bo2(nnext)   )
    do inext = 1,nnext
       indic(inext) = 0
    end do

    kfl_boxxx = 0
    pnodb     = 2_ip
    ipex1     = 0_ip
    ipex2     = 0_ip
    !
    ! Check number of extension boundaries we have
    ! BOUND = 0 ....... both boundaries
    !       = BOUN1 ... boundary 1
    !       = BOUN2 ... boundary 2
    !
    bound = 0
    cont1 = 0
    cont2 = 0
 
    if( lboch(boun1) == BOEXT .and. lboch(boun2) == BOEXT ) then
       bound = 0
    else if( lboch(boun1) == BOEXT ) then
       bound = boun1
    else
       bound = boun2
    end if


    if( bound == 0 ) then

       !-----------------------------------------------------------------     
       !
       ! I have two boundaries to extend from 
       !
       ! Normalized norms: NORB1(2), NORB2(2) of NBOU1 and NBOU2
       !
       !-----------------------------------------------------------------

       ipoi1_boun1 = lnodb(1,boun1)
       ipoi2_boun1 = lnodb(2,boun1)
       ipoi1_boun2 = lnodb(1,boun2)
       ipoi2_boun2 = lnodb(2,boun2)

       norb1(1) = -( coord(2,ipoi1_boun1) - coord(2,ipoi2_boun1) ) ! n1(x)
       norb1(2) =  ( coord(1,ipoi1_boun1) - coord(1,ipoi2_boun1) ) ! n1(y)
       norb2(1) = -( coord(2,ipoi1_boun2) - coord(2,ipoi2_boun2) ) ! n2(x)
       norb2(2) =  ( coord(1,ipoi1_boun2) - coord(1,ipoi2_boun2) ) ! n2(y)
       norme    =  1.0_rp / sqrt(norb1(1)*norb1(1)+norb1(2)*norb1(2))
       norb1(1) =  norb1(1) * norme
       norb1(2) =  norb1(2) * norme
       norme    =  1.0_rp / sqrt(norb2(1)*norb2(1)+norb2(2)*norb2(2))
       norb2(1) =  norb2(1) * norme
       norb2(2) =  norb2(2) * norme
       !
       ! Check open/close corners
       !
       kfl_apert = -1_ip 

       do inodb = 1,pnodb
          !
          ! Normalized tangent vectors TANG1(2), TANG2(2) of NBOU1 and NBOU2
          !
          if( lnodb(inodb,boun1) /= ipoin ) then
             inod1     = inodb
             tang1(1)  = coord(1,lnodb(inodb,boun1)) - coord(1,ipoin)
             tang1(2)  = coord(2,lnodb(inodb,boun1)) - coord(2,ipoin)
             nord1     = 1.0_rp / sqrt(tang1(1) * tang1(1) + tang1(2) * tang1(2))
             tang1(1)  = tang1(1) * nord1 
             tang1(2)  = tang1(2) * nord1
          end if

          if( lnodb(inodb,boun2) /= ipoin ) then
             inod2     = inodb
             tang2(1)  = coord(1,lnodb(inodb,boun2)) - coord(1,ipoin)
             tang2(2)  = coord(2,lnodb(inodb,boun2)) - coord(2,ipoin)
             nord2     = 1.0_rp / sqrt(tang2(1) * tang2(1) + tang2(2) * tang2(2))
             tang2(1)  = tang2(1) * nord2 
             tang2(2)  = tang2(2) * nord2
          end if
       end do
       !              TANG2
       !                _      x
       !     TANG1      /|    /
       !      <--      /    /  BOUN2
       !                  /
       !   x-------------o
       !     BOUN1        IPOIN
       !
       coseno_aux3 = tang1(1) * tang2(1) + tang1(2) * tang2(2)
       !
       ! Check if we have a right angle => KFL_ESQUI = 1
       !
       if( abs(coseno_aux3) < cos_angle_max ) kfl_esqui = 1
       !
       ! Average normal NORMP = 1/2 ( NORB1 + NORB2 ) 
       !
       normp(1) = norb1(1) + norb2(1)
       normp(2) = norb1(2) + norb2(2)
       norma    = 1.0_rp / sqrt( normp(1) * normp(1) + normp(2)*normp(2) )
       normp(1) = normp(1) * norma
       normp(2) = normp(2) * norma
       !
       ! COSENO_AUX3, COSENO_AUX4: scalar products between averaged normal and tangents
       !
       coseno_aux3 = tang1(1) * normp(1) + tang1(2) * normp(2)
       coseno_aux4 = tang2(1) * normp(1) + tang2(2) * normp(2)
       !
       ! Case KFL_ESQUI = 1: what kind of corner?
       !
       !  KFL_APERT = 0         KFL_APERT = 1
       !      x-------
       !      |
       !      |                   x-----o
       ! x----o                         |   
       ! |                  subdomain   |
       ! |    subdomain                 x
       !
       if( kfl_esqui == 1 ) then
          if( coseno_aux3 > 0.0_rp ) then
             kfl_apert = 0                      ! se ven los ojos
          else 
             kfl_apert = 1                      ! se ven el culo
          end if
       end if
       !
       ! Loop over candidates: VECEX = VEC(IPOIN,CANDIDATE)
       !
       inexb1 = 0
       inexb2 = 0
       do inext = 1,nnext
          pexti    = lnext(inext)          
          vecex(1) = coord(1,pexti) - coord(1,ipoin)
          vecex(2) = coord(2,pexti) - coord(2,ipoin)
          norme    = sqrt( vecex(1) * vecex(1) + vecex(2) * vecex(2) )

          if( norme < toler_distance ) then
             !
             ! Candidate is on me! Remove from list
             !
             cosb1 = 0.0_rp
             cosb2 = 0.0_rp
          else
             !
             ! Compute angles between Tangents and VEC(IPOIN,CANDIDATE)
             !
             norme       = 1.0_rp / norme
             vecex(1)    = vecex(1) * norme
             vecex(2)    = vecex(2) * norme
             coseno_aux1 = vecex(1) * tang1(1) + vecex(2) * tang1(2)
             coseno_aux2 = vecex(1) * tang2(1) + vecex(2) * tang2(2)
             cosa        = vecex(1) * normp(1) + vecex(2) * normp(2)
             cosb1       = vecex(1) * norb1(1) + vecex(2) * norb1(2)
             cosb2       = vecex(1) * norb2(1) + vecex(2) * norb2(2)
             !
             ! Take only candidates located between -60 and 60 degrees: (x)
             ! Do not consider others (-)
             !
             !          x
             !       x         -
             !     \     /    
             ! -    \   /  -      -
             !   -   \ /
             ! o------o------o
             !
             if( cosb1 > 0.4_rp ) then
                inexb1 = inexb1 + 1            ! For BOUN1
                lpoex_bo1(inexb1) = pexti
             end if
             if( cosb2 > 0.4_rp ) then
                inexb2 = inexb2 + 1            ! For BOUN2
                lpoex_bo2(inexb2) = pexti
             end if
          end if
       end do
       !
       ! Look for neighboring boundary BOUN3 to BOUN1
       !
       !               jpoin    BOUN3
       !                   x-----------x
       !                   |
       !                   | BOUN1
       !                   |
       !        x----------o   
       !        |  BOUN2    ipoin
       !  BOUN4 |    
       !        |       subdomain
       !        x
       !      kpoin
       !
       jpoin = lnodb(1,boun1)
       do ibopo = pbopo(jpoin),pbopo(jpoin+1)-1
          cont = 0
          do inodb = 1,pnodb
             if( lnodb(inodb,lbopo(ibopo)) /= ipoin ) then 
                cont = cont + 1
             end if
          end do
          if( cont == pnodb ) boun3 = lbopo(ibopo)
       end do
       !
       ! DIRC3 = VEC(BOUN3)
       !
       dirc3(1) = coord(1,lnodb(1,boun3)) - coord(1,lnodb(2,boun3))
       dirc3(2) = coord(2,lnodb(1,boun3)) - coord(2,lnodb(2,boun3))
       norme    = 1.0_rp / sqrt( dirc3(1) * dirc3(1) + dirc3(2) * dirc3(2))
       dirc3(1) = dirc3(1) * norme 
       dirc3(2) = dirc3(2) * norme

       kpoin = lnodb(2,boun2)
       do ibopo = pbopo(kpoin),pbopo(kpoin+1)-1
          cont = 0
          do inodb = 1,pnodb
             if( lnodb(inodb,lbopo(ibopo)) /= ipoin ) then !no es este
                cont = cont + 1
             end if
          end do
          if( cont == pnodb ) boun4 = lbopo(ibopo)
       end do
       !
       ! DIRC4 = VEC(BOUN4)
       !
       dirc4(1)  = coord(1,lnodb(2,boun4)) - coord(1,lnodb(1,boun4))
       dirc4(2)  = coord(2,lnodb(2,boun4)) - coord(2,lnodb(1,boun4))
       norme     = 1.0_rp / sqrt( dirc4(1) * dirc4(1) + dirc4(2) * dirc4(2) )
       dirc4(1)  = dirc4(1) * norme 
       dirc4(2)  = dirc4(2) * norme
       !
       ! Angles between TANG1.DIRC3 and TANG2.DIRC4
       !
       cosen1_nei = tang1(1) * dirc3(1) + tang1(2) * dirc3(2)
       cosen2_nei = tang2(1) * dirc4(1) + tang2(2) * dirc4(2)

       if( abs(cosen1_nei) < 1.0e-7_rp .or. abs(cosen2_nei) < 1.0e-7_rp ) kfl_boxxx = 1
       !
       ! Candidates for BOUN1
       !
       do inext = 1,inexb1
          pexti = lpoex_bo1(inext)
          pextj = lnodb(inod1,boun1)
          do idime = 1,2
             coord_ipex1(idime) = coord(idime,pextj)
             coord_ipex2(idime) = coord(idime,pexti)
          end do
          call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pextj,pexti,q,perm)
          quali_bo1(inext) = q
       end do
       !
       ! Order quality for BOUN1 triangles
       !
       do inext = 1,inexb1-1
          do jnext = inext+1,inexb1
             if( quali_bo1(inext) > quali_bo1(jnext) ) then
                tmp              = lpoex_bo1(inext)
                lpoex_bo1(inext) = lpoex_bo1(jnext)
                lpoex_bo1(jnext) = tmp
                qaux             = quali_bo1(inext)
                quali_bo1(inext) = quali_bo1(jnext)
                quali_bo1(jnext) = qaux
             end if
          end do
       end do
       !
       ! Candidates for BOUN2
       !       
       do inext = 1,inexb2
          pexti = lpoex_bo2(inext)
          pextj = lnodb(inod2,boun2)
          coord_ipex1(1) = coord(1,pextj)
          coord_ipex1(2) = coord(2,pextj)
          coord_ipex2(1) = coord(1,pexti)
          coord_ipex2(2) = coord(2,pexti)
          call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pextj,pexti,q,perm)
          quali_bo2(inext) = q
       end do
       !
       ! Order quality for BOUN2 triangles
       !
       do inext = 1,inexb2-1
          do jnext = inext+1,inexb2
             if( quali_bo2(inext) > quali_bo2(jnext) ) then
                tmp              = lpoex_bo2(inext)
                lpoex_bo2(inext) = lpoex_bo2(jnext)
                lpoex_bo2(jnext) = tmp
                qaux             = quali_bo2(inext)
                quali_bo2(inext) = quali_bo2(jnext)
                quali_bo2(jnext) = qaux
             end if
          end do
       end do
       !
       ! Choose best combination
       !
       q_old = 1000.0_rp
       cont  = 0
       do inext = 1,inexb1
          do jnext = 1,inexb2

             if( lpoex_bo1(inext) /= lpoex_bo2(jnext) ) then
                !
                ! Both extensions have chosen a different node
                !
                q1             = 0.0_rp
                q2             = 0.0_rp
                pexti          = lpoex_bo1(inext)
                pextj          = lpoex_bo2(jnext)
                coord_ipex1(1) = coord(1,pexti)
                coord_ipex2(1) = coord(1,pextj)
                coord_ipex1(2) = coord(2,pexti)
                coord_ipex2(2) = coord(2,pextj)

                vecex_bo1(1)   = coord_ipex1(1)-coord_ipoin(1)
                vecex_bo1(2)   = coord_ipex1(2)-coord_ipoin(2)
                norme          = 1.0_rp / sqrt( vecex_bo1(1) * vecex_bo1(1) + vecex_bo1(2) * vecex_bo1(2) )
                vecex_bo1(1)   = vecex_bo1(1) * norme
                vecex_bo1(2)   = vecex_bo1(2) * norme

                vecex_bo2(1)   = coord_ipex2(1)-coord_ipoin(1)
                vecex_bo2(2)   = coord_ipex2(2)-coord_ipoin(2)
                norme          = 1.0_rp / sqrt( vecex_bo2(1) * vecex_bo2(1) + vecex_bo2(2) * vecex_bo2(2) )
                vecex_bo2(1)   = vecex_bo2(1) * norme
                vecex_bo2(2)   = vecex_bo2(2) * norme

                coseno_aux4    = vecex_bo1(1) * tang1(1) + vecex_bo1(2) * tang1(2)
                coseno_aux5    = vecex_bo2(1) * tang2(1) + vecex_bo2(2) * tang2(2)
                coseno_aux6    = vecex_bo2(1) * tang1(1) + vecex_bo2(2) * tang1(2)

                if( coseno_aux6 > coseno_aux4 ) q2 = 2001.0_rp

                if( kfl_esqui == 1 .and. kfl_apert == 0 ) then             
                   if( coseno_aux4 < 0.0_rp ) q1 = 2001.0_rp                
                   if( coseno_aux5 < 0.0_rp ) q2 = 2001.0_rp
                   if( coseno_aux4 > 0.0_rp .and. coseno_aux5 > 0.0_rp .and. coseno_aux6 > coseno_aux4 ) q2 = 2001.0_rp
                end if
                if( kfl_boxxx == 1 .and. kfl_apert /= 1 ) then 
                   if( coseno_aux4 < 0.5_rp ) q1 = 2001.0_rp
                   if( coseno_aux5 < 0.5_rp ) q2 = 2001.0_rp
                end if

                call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pexti,pextj,q,perm)
                q_res = max(quali_bo1(inext),quali_bo2(jnext),q,q1,q2)

             else

                q_res = max(quali_bo1(inext),quali_bo2(jnext))

             end if

             if( q_res <= q_old ) then
                cont  = cont +  1
                q_old = q_res
                if( lpoex_bo1(inext) /= lpoex_bo2(jnext) ) then
                   ipex1 = lpoex_bo1(inext)
                   ipex2 = lpoex_bo2(jnext)
                   ntext = 3
                else
                   ntext = 2
                   ipex1 = lpoex_bo1(inext)
                   ipex2 = 0_ip
                end if
             end if

          end do
       end do
       !
       ! No good extension elements have been found
       !
       if( cont == 0 ) then

          do inext = 1, inexb1
             lnext_tmp(inext) = lpoex_bo1(inext)
          end do
          kk = inexb1
          do inext = 1,inexb2
             kk = kk + 1
             lnext_tmp(kk) = lpoex_bo2(inext)
          end do

          q_old = 1000.0_rp          

          do inext = 1,kk
             pexti          = lnext_tmp(inext)
             pextj          = lnodb(1,boun1)
             pextk          = lnodb(2,boun2)
             coord_ipex1(1) = coord(1,pexti)
             coord_ipex1(2) = coord(2,pexti)
             coord_ipex2(1) = coord(1,pextj)
             coord_ipex2(2) = coord(2,pextj)
             coord_ipex3(1) = coord(1,pextk)
             coord_ipex3(2) = coord(2,pextk)

             vecex_bo1(1)   = coord_ipex1(1)-coord_ipoin(1)
             vecex_bo1(2)   = coord_ipex1(2)-coord_ipoin(2)
             norme          = 1.0_rp / sqrt( vecex_bo1(1) * vecex_bo1(1) + vecex_bo1(2) * vecex_bo1(2) )
             vecex_bo1(1)   = vecex_bo1(1) * norme
             vecex_bo1(2)   = vecex_bo1(2) * norme

             call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pexti,pextj,q1,perm)
             call dod_extens_qual2d(coord_ipoin,coord_ipex3,coord_ipex1,pextk,pexti,q2,perm)
             q_res = max(q1,q2)

             if( kfl_esqui == 1 .and. kfl_apert == 0 ) then
                coseno_aux4 = vecex_bo1(1) * tang1(1) + vecex_bo1(2) * tang1(2)
                coseno_aux5 = vecex_bo1(1) * tang2(1) + vecex_bo1(2) * tang2(2)
             end if
             if( kfl_boxxx == 1 .and. kfl_apert /= 1 ) then  
                if( coseno_aux4 < 0.5_rp ) quali_bo1(inext) = 2001.0_rp
                if( coseno_aux5 < 0.5_rp ) quali_bo2(jnext) = 2001.0_rp
             end if

             if( coseno_aux4 < 0.0_rp .or. coseno_aux5 < 0.0_rp ) q_res = 2001.0_rp                

             if( q_res <= q_old ) then
                cont  = cont + 1
                q_old = q_res
                ntext = 2
                ipex1 = lnext_tmp(inext)
                ipex2 = 0                
             end if
          end do

       end if

       if( cont == 0 ) then
          print*,ipoin,lnext(1:nnext)     
          call runend('DOD_NORM2D: I RUN OUT OF IDEA TO CREATE EXTENSIONS...')
       end if

    else

       !-----------------------------------------------------------------     
       !
       ! I have only one boundary to extend from: BOUND
       !
       !-----------------------------------------------------------------

       comp1(1) =  coord(1,lnodb(1,bound))-coord(1,lnodb(2,bound))
       comp1(2) =  coord(2,lnodb(1,bound))-coord(2,lnodb(2,bound))
       norb1(1) = -comp1(2)
       norb1(2) =  comp1(1)
       do inodb = 1,pnodb
          if( lnodb(inodb,bound) /= ipoin ) then
             inod1     = inodb
             tang1(1)  = coord(1,lnodb(inodb,bound)) - coord(1,ipoin)
             tang1(2)  = coord(2,lnodb(inodb,bound)) - coord(2,ipoin)
             nord1     = 1.0_rp / sqrt( tang1(1) * tang1(1) + tang1(2) * tang1(2) )
             tang1(1)  = tang1(1) * nord1 
             tang1(2)  = tang1(2) * nord1
          end if
       end do
       inexb1 = 0
       inexb2 = 0
       do inext = 1,nnext

          pexti    = lnext(inext)
          vecex(1) = coord(1,pexti)-coord_ipoin(1)
          vecex(2) = coord(2,pexti)-coord_ipoin(2)
          norme    = sqrt( vecex(1) * vecex(1) + vecex(2) * vecex(2) )

          if( norme < toler_distance ) then
             cosb1 = 0.0_rp
             cosb2 = 0.0_rp
          else
             vecex(1)    = vecex(1) / norme
             vecex(2)    = vecex(2) / norme
             coseno_aux1 = vecex(1) * tang1(1) + vecex(2) * tang1(2)

             cosb1 = vecex(1)*norb1(1) + vecex(2)*norb1(2)
             if( cosb1 < 0.4_rp ) then
                inexb1 = inexb1 + 1
                lpoex_bo1(inexb1) = pexti
             end if
          end if
       end do

       do inext = 1,inexb1
          pexti = lpoex_bo1(inext)
          pextj = lnodb(inod1,bound)
          coord_ipex1(1) = coord(1,pextj)
          coord_ipex1(2) = coord(2,pextj)
          coord_ipex2(1) = coord(1,pexti)
          coord_ipex2(2) = coord(2,pexti)
          call dod_extens_qual2d(coord_ipoin,coord_ipex1,coord_ipex2,pextj,pexti,q,perm)
          quali_bo1(inext) = q
       end do

       do inext = 1,inexb1-1
          do jnext = inext+1,inexb1
            if( quali_bo1(inext) > quali_bo1(jnext) ) then
                tmp              = lpoex_bo1(inext)
                lpoex_bo1(inext) = lpoex_bo1(jnext)
                lpoex_bo1(jnext) = tmp
                qaux             = quali_bo1(inext)
                quali_bo1(inext) = quali_bo1(jnext)
                quali_bo1(jnext) = qaux
             end if
          end do
       end do

       ipex1 = lpoex_bo1(1)
       ipex2 = 0

    end if

    deallocate( quali_bo2 )
    deallocate( quali_bo1 )
    deallocate( indic     )
    deallocate( lnext_tmp )

  end subroutine dod_extension_elements_2d

  subroutine dod_extension_elements_3d(itask,epoi1,epoi2,epoi3,normal_x,normal_y,normal_z,orien,&
       nnext,lnext,isubd,jsubd,iboun,tboun )
    use def_parame
    use def_elmtyp
    use def_master
    use def_domain
    use def_dodeme
    implicit none
    integer(ip),intent(in)    :: itask,epoi1,epoi2,epoi3,isubd,jsubd,iboun,tboun
    integer(ip),intent(inout) :: nnext,lnext(nnext)
    real(rp) ,intent(in   )   :: orien(3)
    real(rp) ,intent(inout)   :: normal_x, normal_y, normal_z
    integer(ip)               :: idime,pexti,cont,lnext_tmp(nnext),inext,jnext,kk,tmp,mnext,sign
    integer(ip)               :: jelem,inode,ipoin,ipoin_subd,jpoin,kpoin,kpoin_subd
    integer(ip)               :: ii,pnode,minim,flag,izdom,epoi2_aux,pexti_aux
    integer(4)                :: istat
    real(rp)                  :: norma,norme,cosa,kappa,asrad
    real(rp)                  :: coor1(ndime),coor2(ndime),coor3(ndime),coocg(ndime)
    real(rp)                  :: vec12(ndime),vec13(ndime),vecex(ndime),q,qaux
    real(rp), allocatable     :: qual(:) , qual_aux(:)
    integer(ip), allocatable  :: lnext_aux(:),lnext_sec(:) 
    real(rp)                  :: coori(3),xieta(3)
    select case(itask) 

    case(1)
       do idime = 1,ndime
          coor1(idime)  = coord(idime,epoi1)
       end do
       do idime = 1,ndime
          coor2(idime)  = coord(idime,epoi2)
       end do
       do idime = 1,ndime
          coor3(idime)  = coord(idime,epoi3)
       end do
       do idime = 1,ndime
          vec12(idime) = coor2(idime)-coor1(idime)
       end do
       do idime = 1,ndime
          vec13(idime) = coor3(idime)-coor1(idime)
       end do
       normal_x = vec12(2)*vec13(3) - vec12(3)*vec13(2)
       normal_y = vec12(3)*vec13(1) - vec12(1)*vec13(3)
       normal_z = vec12(1)*vec13(2) - vec12(2)*vec13(1)

       norma = sqrt(normal_x*normal_x + normal_y*normal_y + normal_z*normal_z)

       normal_x = normal_x/norma
       normal_y = normal_y/norma
       normal_z = normal_z/norma

    case(2)
       allocate(qual(nnext),stat=istat)
       allocate(lnext_sec(nnext),stat=istat)

       do inext =1,nnext
          lnext_tmp(inext)=0
       end do
       !
       !Otro criterio para el parabrisas...
       !
       do idime = 1,ndime
          coor1(idime)  = coord(idime,epoi1)
       end do
       do idime = 1,ndime
          coor2(idime)  = coord(idime,epoi2)
       end do
       do idime = 1,ndime
          coor3(idime)  = coord(idime,epoi3)
       end do
       do idime = 1,ndime
          coocg(idime)  = (coor1(idime)+coor2(idime)+coor3(idime))/3.0_rp
       end do

       do inext = 1,nnext
          pexti             = lnext(inext)
          lnext_sec(inext)  = lnext(inext)
          if(pexti/=0)then
             do idime = 1,ndime
                vecex(idime) = coord(idime,pexti)-coocg(idime)
             end do
             norme = 0
             do idime = 1,ndime
                norme    = norme + vecex(idime)*vecex(idime)
             end do
             norme = sqrt(norme)
             do idime = 1,ndime
                vecex(idime)  = vecex(idime)/norme
             end do
             cosa = vecex(1)*normal_x + vecex(2)*normal_y + vecex(3)*normal_z
             call dod_qualit(epoi1,epoi2,epoi3,pexti,q,sign)!OJOOOOOO
             epoi2_aux = epoi2
             pexti_aux = pexti
             !call dod_extens_qual3d(2_ip,epoi1,epoi2_aux,epoi3,pexti_aux,kappa,asrad,q,sign,orien)

!print*,epoi1,epoi2,epoi3,pexti,q,cosa
             If(cosa<0.2_rp .and. cosa>0.0_rp)then
                if(q>7.0_rp)then
                   lnext(inext)=0
                else
                   qual(inext)= q
                end if
             else if (q>40.0_rp .or. cosa<0.0_rp)then
                lnext(inext)=0
             else
                qual(inext)= q
             end if
          end if
       end do
       do inext=1,nnext-1
          do jnext=inext+1,nnext
             if(qual(inext)>qual(jnext))then
                tmp = lnext(inext)
                lnext(inext)=lnext(jnext)
                lnext(jnext)=tmp
                qaux = qual(inext)
                qual(inext)= qual(jnext)
                qual(jnext)= qaux
             end if
          end do
       end do
       cont  = 0
       kk = 0
       do inext =1,nnext
          if(lnext(inext)/=0)then
             cont = cont + 1
             lnext_tmp(cont)= lnext(inext)
             lnext(inext) = 0
          else 
             kk = kk + 1
          end if
       end do
!!! nnext = nnext-kk
       cont = 0
       do inext = 1,nnext
          lnext(inext) = lnext_tmp(inext)
          if(lnext(inext) == 0) cont = cont + 1
       end do
       if(cont==nnext)then
          print*,epoi1,'ERROR:DOD_NORMAL:SIN LNEXT sigo probando...',lnext_sec
          !stop
          mnext = 200 ! amplio lnext!!!!
          allocate(qual_aux(mnext),stat=istat)
          allocate(lnext_aux(mnext),stat=istat)
          do inext= 1,mnext
             qual_aux(inext)  = 0.0_rp
             lnext_aux(inext) = 0_ip
          end do
          ii = 0
          do inext=1,nnext
             if(lnext_sec(inext)/=0)then
                jpoin          = lnext_sec(inext)
                !! do izdom =ppopo(jpoin), ppopo(jpoin+1)-1
                !!    kpoin   = lpopo(izdom) !!!!local!!!!!???
                !!       if(subdomain(jsubd)% lsubd_npoin(kpoin_subd) == 0)then  !no soy punto de extension
                if(ipres_dod(isubd) == 1 )then
                   !if(lpobo_dod(epoi1)==1.and. lpobo_dod(kpoin)==1 ) ipoin=kpoin
                   print*,'TIPO LEGO:'
                   print*,'TENGO QUE COMPROBAR QUE ESCOJO A UN PUNTO DE CONTORNO REAL SI SOY CONTORNO REAL'
                   stop
                else
                   ipoin=jpoin
                end if

                do idime = 1,ndime
                   vecex(idime) = coord(idime,ipoin)-coocg(idime)
                end do

                norme = 0
                do idime = 1,ndime
                   norme    = norme + vecex(idime)*vecex(idime)
                end do

                norme = sqrt(norme)
                do idime = 1,ndime
                   vecex(idime)  = vecex(idime)/norme
                end do

                cosa = vecex(1)*normal_x + vecex(2)*normal_y + vecex(3)*normal_z
                call dod_qualit(epoi1,epoi2,epoi3,ipoin,q,sign)
                if(cosa>0.15_rp.and. q<40_rp )then
                   ii = ii + 1
                   qual_aux(ii)= q
                   lnext_aux(ii)= ipoin                            
                end if

             end if
             !!       end do
             !!          end if
          end do

          if(ii==0)then
             print*,epoi1,'ERROR:DOD_EXTENSION_ELEMENTS_3D:SIN LNEXT','contorno',epoi1,epoi2,epoi3
             print*,'la lista de la que partia',lnext_sec(1:nnext)
             stop
          else if(ii==1)then
             lnext(ii)=lnext_aux(ii)
          else
             cont = 0
             do inext = 1,ii-1
                do jnext = inext+1,ii
                   if(lnext_aux(inext) == lnext_aux(jnext)  )then
                      lnext_aux(jnext) = 0
                      qual_aux(jnext) = 1000.0_rp
                      cont = cont + 1 
                   end if
                end do
             end do

             do inext=1,ii-1
                do jnext=inext+1,ii
                   if(qual_aux(inext)>qual_aux(jnext))then
                      tmp = lnext_aux(inext)
                      lnext_aux(inext)=lnext_aux(jnext)
                      lnext_aux(jnext)=tmp
                      qaux = qual_aux(inext)
                      qual_aux(inext)= qual_aux(jnext)
                      qual_aux(jnext)= qaux
                   end if
                end do
             end do


             minim = min(nnext,ii-cont)
             inext = 1
             jnext = 0
             if(ii>nnext)then
                cont = 0
                do inext = 1,nnext
                   if(lnext_aux(inext)/=0_ip)then
                      cont = cont + 1
                      lnext(cont)=lnext_aux(inext)
                   end if
                end do

             else if(ii<=nnext)then
                cont = 0
                do inext = 1,ii
                   if(lnext_aux(inext)/=0_ip)then
                      cont = cont + 1
                      lnext(cont)=lnext_aux(inext)
                   end if
                end do

             end if

             deallocate(qual_aux,stat=istat)
             deallocate(lnext_aux,stat=istat)

          end if
       end if
       deallocate(qual,stat=istat)
       deallocate(lnext_sec,stat=istat)

    end select

  end subroutine dod_extension_elements_3d

  subroutine dod_extens_qual3d(itask,ipoi1,ipoi2,ipoi3,ipoi4,kappaS,asrad,gamma,sign,norma)
    use def_dodeme
    use def_kintyp
    use def_parame
    use def_master
    use def_domain
    implicit none  
    integer(ip),intent(in)      :: itask
    real(rp),intent(in)         :: norma(3)
    integer(ip),intent(in)      :: ipoi1,ipoi3
    integer(ip),intent(inout)   :: ipoi2,ipoi4
    integer(ip),intent(out)     :: sign
    real(rp)   ,intent(out)     :: kappaS,asrad,gamma
    integer(ip)                 :: idime,inode,jnode,tetra,ii,jj,pnode,pelty
    integer(ip)                 :: lnods_tmp(mnode),jpoin_tmp
    real(rp)                    :: invW(ndime,3),tinvW(3,ndime),invWt(3,ndime)
    real(rp)                    :: A(ndime,3),S(ndime,3),invS(ndime,3)
    real(rp)                    :: normS,norminvS,StS(ndime,3)
    real(rp)                    :: invSt(ndime,3)
    real(rp)                    :: aux(ndime,3),aux2(ndime,3)
    real(rp)                    :: detinvWt,detS,detSt
    real(rp)                    :: elcod(ndime,mnode),volum
    real(rp)                    :: xjaci(ndime,ndime),xjacm(ndime,ndime) ,detjm
    real(rp)                    :: cartd(ndime,mnode) ,sumsi,cosa
    real(rp)                    :: side1(ndime),side2(ndime),side3(ndime)
    real(rp)                    :: side4(ndime),side5(ndime),side6(ndime)
    real(rp)                    :: leng1,leng2,leng3,leng4,leng5,leng6,hmax
    real(rp)                    :: surf1,surf2,surf3,surf4,surf,alpha,increm
    real(rp)                    :: orien(3),coocg(3),coor1(3),coor2(3),coor3(3)

    gamma  = 0.0_rp
    asrad  = 0.0_rp
    kappaS = 0.0_rp
    sign   = 1_ip
    increm = epsilon(1.0_rp)


    lnods_tmp(1)=ipoi1
    lnods_tmp(2)=ipoi2
    lnods_tmp(3)=ipoi3
    lnods_tmp(4)=ipoi4

    if(itask==1)then
       alpha = sqrt(6.0_rp)/12.0_rp
       do idime=1,ndime
          side1(idime) = coord(idime,lnods_tmp(1))-coord(idime,lnods_tmp(2))  
          side2(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(2)) 
          side3(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(1)) 
          side4(idime) = coord(idime,lnods_tmp(3))-coord(idime,lnods_tmp(1)) 
          side5(idime) = coord(idime,lnods_tmp(3))-coord(idime,lnods_tmp(2)) 
          side6(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(3)) 
       end do


       !
       ! hmax
       !
       leng1=0.0_rp
       leng2=0.0_rp
       leng3=0.0_rp
       leng4=0.0_rp
       leng5=0.0_rp
       leng6=0.0_rp

       do idime=1,ndime 
          leng1 = leng1 + side1(idime)*side1(idime) 
          leng2 = leng2 + side2(idime)*side2(idime)
          leng3 = leng3 + side3(idime)*side3(idime)
          leng4 = leng4 + side4(idime)*side4(idime) 
          leng5 = leng5 + side5(idime)*side5(idime) 
          leng6 = leng6 + side6(idime)*side6(idime)
       end do


       leng1 = sqrt(leng1)
       leng2 = sqrt(leng2)
       leng3 = sqrt(leng3)
       leng4 = sqrt(leng4)
       leng5 = sqrt(leng5)
       leng6 = sqrt(leng6)
       hmax = max(leng1,leng2,leng3,leng4,leng5,leng6)
       !
       ! total surface
       !
       surf1 = ( side1(1)*side2(2) - side1(2)*side2(1) )*( side1(1)*side2(2) - side1(2)*side2(1) ) +&
            ( side1(3)*side2(1) - side1(1)*side2(3) )*( side1(3)*side2(1) - side1(1)*side2(3) ) +&
            ( side1(2)*side2(3) - side1(3)*side2(2) )*( side1(2)*side2(3) - side1(3)*side2(2) )
       surf1 =  sqrt(surf1)*0.5_rp 


       surf2 = ( side3(1)*side4(2) - side4(1)*side3(2) )*( side3(1)*side4(2) - side4(1)*side3(2) ) +&
            ( side3(3)*side4(1) - side3(1)*side4(3) )*( side3(3)*side4(1) - side3(1)*side4(3) ) +&
            ( side3(2)*side4(3) - side3(3)*side4(2) )*( side3(2)*side4(3) - side3(3)*side4(2) )
       surf2 =  sqrt(surf2)*0.5_rp


       surf3 = ( side1(1)*side5(2) - side5(1)*side1(2) )*( side1(1)*side5(2) - side5(1)*side1(2) ) +&
            ( side1(3)*side5(1) - side1(1)*side5(3) )*( side1(3)*side5(1) - side1(1)*side5(3) ) +&
            ( side1(2)*side5(3) - side1(3)*side5(2) )*( side1(2)*side5(3) - side1(3)*side5(2) )
       surf3 =  sqrt(surf3)*0.5_rp

       surf4 = ( side6(1)*side2(2) - side2(1)*side6(2) )*( side6(1)*side2(2) - side2(1)*side6(2) ) +&
            ( side6(3)*side2(1) - side6(1)*side2(3) )*( side6(3)*side2(1) - side6(1)*side2(3) ) +&
            ( side6(2)*side2(3) - side6(3)*side2(2) )*( side6(2)*side2(3) - side6(3)*side2(2) )
       surf4 =  sqrt(surf4)*0.5_rp
       surf = surf1 + surf2 + surf3 + surf4

       !
       ! volume
       !


       volum = (  side1(1)*side3(2)*side4(3) + side1(2)*side3(3)*side4(1) + side1(3)*side3(1)*side4(2) &
            &   - side1(3)*side3(2)*side4(1) - side1(2)*side3(1)*side4(3) - side1(1)*side3(3)*side4(2) ) /6.0_rp

       if (volum <= 0.0_rp) sign=-1_ip

       !
       ! Quality
       !

       gamma = alpha * hmax * surf /(3.0_rp*abs(volum)+increm )


    else if(itask ==2) then
       detinvWt = 0.0_rp
       detS = 0.0_rp
       pnode = 4
       invW(1,1) = 1.0_rp
       invW(1,2) = -1.0_rp/3.0_rp*sqrt(3.0_rp)             !!!-5.77350269189626e-01!!!
       invW(1,3) = -1.0_rp/6.0_rp*sqrt(6.0_rp)             !!-4.08248290463863e-01 
       invW(2,1) = 0.0_rp
       invW(2,2) = 2.0_rp/3.0_rp*sqrt(3.0_rp)              !!!1.15470053837925e+00 
       invW(2,3) = -1.0_rp/6.0_rp*sqrt(6.0_rp)             !!!-4.08248290463863e-01
       invW(3,1) = 0.0_rp
       invW(3,2) = 0.0_rp
       invW(3,3) = 1.0_rp/2.0_rp*sqrt(6.0_rp)              !!!!1.22474487139159e+00

       normS = 0.0_rp
       norminvS=0.0_rp
       do inode=1,pnode-1
          do idime=1,ndime
             A(idime,inode) = 0.0_rp
             S(idime,inode) = 0.0_rp
             StS (idime,inode) = 0.0_rp
             invS(idime,inode) = 0.0_rp
          end do
       end do
       do inode=1,pnode-1
          do idime = 1,ndime
             tinvW (idime,inode)= invW(inode,idime)
          end do
       end do
       call invmtx(tinvW,invWt,detinvWt,ndime)
       ii = 0
       do inode=2,pnode
          ii = ii + 1
          do idime=1,ndime
             A(idime,ii) = coord(idime,lnods_tmp(inode)) - coord(idime,lnods_tmp(1))
          end do
       end do
       do inode =1,pnode-1
          do idime=1,ndime
             do jj = 1,pnode-1
                S(idime,inode) = S(idime,inode) + A(idime,jj) * invW(jj,inode)
             end do
          end do
       end do
       do inode=1,pnode-1
          do idime=1,ndime
             normS = normS + S(idime,inode) * S(idime,inode)  
             do jj = 1,pnode-1
                StS (idime,inode) = StS(idime,inode) + S(jj,idime)* S(jj,inode)
             end do
          end do
       end do
       normS = sqrt(normS)
       do inode=1,pnode-1
          do idime=1,ndime
             aux(idime,inode) = S(idime,inode)
          end do
       end do
       call invmtx(aux,aux2,detS,ndime)
       do inode=1,pnode-1
          do idime=1,ndime
             invS(idime,inode) = aux2(idime,inode) 
          end do
       end do

       do inode=1,pnode-1
          do idime=1,ndime
             norminvS = norminvS + invS(idime,inode) * invS(idime,inode)  
          end do
       end do
       norminvS = sqrt(norminvS)
       kappaS   = normS* norminvS
       kappaS   = kappaS/3.0_rp !!!!para que el eleme ideal tenga kappa 1
       do idime=1,ndime
          elcod(idime,1) = coord(idime,lnods_tmp(1))
          elcod(idime,2) = coord(idime,lnods_tmp(2))
          elcod(idime,3) = coord(idime,lnods_tmp(3))
          elcod(idime,4) = coord(idime,lnods_tmp(4))
       end do
       volum = 0.0_rp
       sumsi = 0.0_rp
       pelty = 30
       call jacobi(&
            ndime,pnode,elcod,elmar(pelty)%dercg,&
            xjacm,xjaci,cartd,detjm)
       volum = volum + elmar(pelty)%weicg*detjm

       do idime=1,ndime
          side1(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(1))
          side2(idime) = coord(idime,lnods_tmp(3))-coord(idime,lnods_tmp(1))
          side3(idime) = coord(idime,lnods_tmp(2))-coord(idime,lnods_tmp(1))
          side4(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(2))
          side5(idime) = coord(idime,lnods_tmp(3))-coord(idime,lnods_tmp(2))
          side6(idime) = coord(idime,lnods_tmp(4))-coord(idime,lnods_tmp(3)) 
       end do

       do idime=1,ndime
          sumsi = sumsi + side1(idime)*side1(idime) + side2(idime)*side2(idime)+ &
               side3(idime)*side3(idime) + side4(idime)*side4(idime)+ &
               side5(idime)*side5(idime) + side6(idime)*side6(idime)
       end do
       asrad = sumsi / (volum+increm)

    else if(itask==3)then

       do idime=1,ndime
          elcod(idime,1) = coord(idime,lnods_tmp(1))
          elcod(idime,2) = coord(idime,lnods_tmp(2))
          elcod(idime,3) = coord(idime,lnods_tmp(3))
          elcod(idime,4) = coord(idime,lnods_tmp(4))
       end do
       pnode=4
       pelty=30
       volum = 0.0_rp
       call jacobi(&
            ndime,pnode,elcod,elmar(pelty)%dercg,&
            xjacm,xjaci,cartd,detjm)
       volum = volum + elmar(pelty)%weicg*detjm

       if(volum < 0.0_rp )then
          jpoin_tmp = ipoi2
          ipoi2     = ipoi4
          ipoi4     = jpoin_tmp
       end if
    else if(itask ==4)then
  ! ipoi1: interface's point
  ! ipoi2: boundary point conected with ipoi1
  ! ipoi3,ipoi4 : the two candidates to form the extension

       sign = 0_ip
       do idime = 1,ndime
          coor1(idime)  = coord(idime,ipoi2)
       end do
       do idime = 1,ndime
          coor2(idime)  = coord(idime,ipoi3)
       end do
       do idime = 1,ndime
          coor3(idime)  = coord(idime,ipoi4)
       end do
       do idime = 1,ndime
          coocg(idime)  = (coor1(idime)+coor2(idime)+coor3(idime))/3.0_rp
       end do
       do idime=1,ndime
          orien(idime) = coocg(idime)-coord(idime,ipoi1)
       end do
       normS = 0.0_rp
       do idime=1,ndime
          normS = normS + orien(idime)*orien(idime) 
       end do
       normS =  sqrt(normS)
       do idime=1,ndime
          orien(idime) = orien(idime)/normS
       end do
       cosa = 0.0_rp
       do idime=1,ndime
          cosa = cosa + norma(idime)*orien(idime) 
       end do
       if(cosa>0.2_rp)sign = 1_ip
       
    end if
    
  end subroutine dod_extens_qual3d
  subroutine dod_extension_correc(q,npoex,nelem_aux,tboun,nnext,ipoin,&
       pobok,poboq,pobor,pobos,poex1,poexk,&
       poexq,poexr,qboun,kboun,rboun,sboun,lpoex,lext2)

  use def_kintyp
  use def_dodeme
  use mod_memchk
  implicit none

  integer(ip),intent(in)        ::npoex,nelem_aux,tboun,nnext,ipoin
  integer(ip),intent(inout)     ::poex1,poexq,poexr,poexk
  integer(ip),intent(inout)     ::pobok,poboq,pobor,pobos
  integer(ip),intent(in)        ::kboun,qboun,rboun,sboun
  integer(ip),intent(in)        ::lpoex(nnext,tboun)
  integer(ip),intent(inout)     ::lext2(tboun)
  real(rp),intent(in)           ::q
  integer(ip)                   ::inext,knext,qnext,rnext,ipoex,kpoex,qpoex,rpoex,cont,ii
  integer(ip)                   ::poex2,poex3,poex4, para(tboun),kk,jj,sign
  real(rp)                      ::q_old,q1,q2,q3,q4,q5,k,flag,dummi(3),dumm1,dumm2
  k = -10000
  para(1) = -1
  do kk=2,tboun
     para(kk)=0
  end do
  do kk=1,tboun-1
     do jj=kboun,tboun
        if(lext2(kk)==lext2(jj).and. para(jj)==0)para(jj)=-kk    
     end do
  end do
  if(nelem_aux==3 ) then
     !if(npoex==2)then
        !print*,'numero de tetras-aux=3 pero numero de puntos de exts=2',poex1,poexk
        !print*,'contornos de cambio: kboun=',kboun, 'qboun',qboun
     !end if
     q_old=q
     do inext = 1,nnext
        ipoex = lpoex(inext,1)
        if(ipoex/=0)then
           do knext = 1,nnext
              kpoex = lpoex(knext,kboun)
              if(kpoex/=0)then
                 !if(kpoex==ipoex)npoex=npoex-1
                 do qnext = 1,nnext
                    qpoex = lpoex(qnext,qboun)
                    if(qpoex/=0)then
                       !if(qpoex==kpoex .or. qpoex==ipoex)npoex=npoex-1
                       !call dod_qualit(ipoin,ipoex,pobok,kpoex,q1)!!!el que toca cambando el iboun=1
                       call dod_extens_qual3d(2_ip,ipoin,kpoex,ipoex,pobok,&  
                            dumm1,dumm2,q1,sign,dummi)!!!el que toca cambando el iboun=1
                       !call dod_qualit(ipoin,kpoex,poboq,qpoex,q2)            !!!el que afecta a 1
                        call dod_extens_qual3d(2_ip,ipoin,qpoex,kpoex,poboq,&  
                            dumm1,dumm2,q2,sign,dummi) !!!el que afecta a 1
                       !call dod_qualit(ipoin,ipoex,pobor,qpoex,q3)            !!!el que afecta a 1
                        call dod_extens_qual3d(2_ip,ipoin,pobor,ipoex,qpoex,&  
                            dumm1,dumm2,q3,sign,dummi) !!!el que afecta a 1
                       if(npoex==3)call dod_extens_qual3d(2_ip,ipoin,ipoex,kpoex,qpoex,&   !!!el catalan
                             dumm1,dumm2,q4,sign,dummi) 
                      if(q1<50 .and.q2<50 .and. q3<50 .and.q4<50)then
                          k=max(q1,q2,q3,q4)
                          if(k<q_old) then
                             q_old = k
                             lext2(1)    = ipoex
                             lext2(kboun)= kpoex
                             lext2(qboun)= qpoex
                          end if
                       end if
                    end if
                 end do
              end if
           end do
        end if
     end do
     if(k==-10000)print*,'NO HE ENCONTRADO ALTERNATIVA',ipoin
     do kk=1,tboun
        if(para(kk)==para(1).and.kk/=kpoex.and.kk/=qpoex.and. kk<kboun   )then
           lext2(kk)= lext2(1)
        end if
        if(para(kk)==para(kboun).and.kk/=ipoex.and.kk/=qpoex.and. kk<qboun)then
           lext2(kk)= lext2(kboun)
        end if
        if(para(kk)==para(qboun).and.kk/=ipoex.and.kk/=kpoex .and. kk>kboun)then
           lext2(kk)= lext2(qboun)
        end if
     end do


  else if(nelem_aux==4) then
     if(npoex<4)then
        if(para(1)==para(qboun))flag=1
        if(para(kboun)==para(rboun))flag=2
        if(para(1)==para(qboun) .and. para(kboun)==para(rboun))flag=3
     end if
     q_old=q
     do inext = 1,nnext
        ipoex = lpoex(inext,1)
        if(ipoex/=0)then
           do knext = 1,nnext
              kpoex = lpoex(knext,kboun)
              if(kpoex/=0 )then 
                 do qnext = 1,nnext
                    qpoex = lpoex(qnext,qboun)
                    if(qpoex/=0)then
                       do rnext=1,nnext
                          rpoex=lpoex(rnext,qboun)
                          if(rpoex/=0)then
                             call dod_extens_qual3d(2_ip,ipoin,kpoex, ipoex ,pobok , dumm1,dumm2,q1,sign,dummi)
                             call dod_extens_qual3d(2_ip,ipoin,qpoex, kpoex ,poboq , dumm1,dumm2,q2,sign,dummi)
                             call dod_extens_qual3d(2_ip,ipoin,rpoex, qpoex ,pobor , dumm1,dumm2,q3,sign,dummi)
                             call dod_extens_qual3d(2_ip,ipoin,pobos, ipoex, rpoex , dumm1,dumm2,q4,sign,dummi)
                             if(npoex==3)then
                                if(ipoex==kpoex)then
                                   poex2=ipoex
                                   poex3=qpoex
                                   poex4=rpoex
                                else if(ipoex==qpoex .or. kpoex==qpoex)then
                                   poex2=ipoex
                                   poex3=kpoex
                                   poex4=rpoex
                                else if(ipoex==rpoex .or. kpoex==rpoex .or. qpoex==rpoex)then
                                   poex2=ipoex
                                   poex3=kpoex
                                   poex4=qpoex
                                end if
                                call dod_extens_qual3d(2_ip,ipoin,poex2,poex3,poex4,dumm1,dumm2,q5,sign,dummi)!!!el catalan
                                if(q1<50 .and.q2<50 .and. q3<50 .and.q4<50 .and.q5<50)then
                                   k=max(q1,q2,q3,q4,q5)
                                   if(k<q_old) then
                                      q_old = k
                                      lext2(1)    = ipoex
                                      lext2(kboun)= kpoex
                                      lext2(qboun)= qpoex
                                      lext2(rboun)= rpoex
                                   end if
                                end if
                             else
                                if(q1<50 .and.q2<50 .and. q3<50 .and.q4<50)then
                                   k=max(q1,q2,q3,q4)
                                   if(k<q_old) then
                                      q_old = k
                                      lext2(1)    = ipoex
                                      lext2(kboun)= kpoex
                                      lext2(qboun)= qpoex
                                      lext2(rboun)= rpoex
                                   end if
                                end if
                             end if
                          end if
                       end do
                    end if
                 end do
              end if
           end do
        end if
     end do
     do kk=1,tboun
        if(para(kk)==para(1).and.kk/=kpoex.and.kk/=qpoex.and.kk/=rpoex .and. kk<kboun   )then
           lext2(kk)= lext2(1)
        end if
        if(para(kk)==para(kboun).and.kk/=ipoex.and.kk/=qpoex.and.kk/=rpoex .and. kk<qboun)then
           lext2(kk)= lext2(kboun)
        end if
        if(para(kk)==para(qboun).and.kk/=ipoex.and.kk/=kpoex.and.kk/=rpoex .and. kk>kboun)then
           lext2(kk)= lext2(qboun)
        end if
        if(para(kk)==para(rboun).and.kk/=1.and.kk/=kboun.and.kk/=qboun .and. kk>qboun)then
           lext2(kk)= lext2(rboun)
        end if
     end do
  end if

  end subroutine dod_extension_correc



end module mod_dod_extens
