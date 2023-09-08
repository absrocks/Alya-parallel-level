subroutine sld_bitcul(&
     itask,ielem,pnode,lnods,pgaus,pmate,gpsha,gpvol,gpgdi,&
     gppio,nfibe,gpdet)
  !-----------------------------------------------------------------------
  !****f* Solidz/sld_bitcul
  ! NAME
  !    sld_bitcul
  ! DESCRIPTION
  !    Projection to nodes
  ! INPUT
  !    GPGDI ....... Deformation tensor ...................... F = grad(phi)
  !    GPPIO ....... 1st Piola-Kirchhoff stress tensor ....... P = F.S
  !    GPDET ....... Deformation jacobian determinant ........ J = |F|
  ! OUTPUT
  !    CAUST_SLD ... Cauchy Stress tensor
  !    GREEN_SLD ... Strain tensor
  !    LEPSI_SLD ... Logarithmic strain tensor
  !    FIBEG_SLD ... Fibers on element
  !    FIBDE_SLD ... Fiber's vector output
  ! USES
  ! USED BY
  !    sld_elmope
  !***
  !-----------------------------------------------------------------------
  use def_kintyp, only       :  ip,rp
  !use def_elmtyp, only       :  ELEXT
  use def_elmtyp
  use def_parame, only       :  pi                                                ! ( *AQU* )
  use def_master, only       :  gpfib
  use def_domain, only       :  ndime,lelch,ltype,coord,vmass,elmar,npoin,&
                                mnode
  use def_solidz, only       :  nvoig_sld
  use def_solidz, only       :  nvgij_sld,caust_sld,green_sld,caunp_sld
  use def_solidz, only       :  lepsi_sld,fibde_sld,fibeg_sld
  use def_solidz, only       :  nopio_sld,cause_sld
  use def_solidz, only       :  kfl_foten_sld,kfl_fiber_sld
  use def_solidz, only       :  caulo_sld,parco_sld                               ! ( *AQU* )
  use def_solidz, only       :  grlst_sld
  use def_solidz, only       :  lepse_sld, epsel_sld
  use def_solidz, only       :  epsee_sld, svegm_sld, kfl_plast_sld

  implicit none

  integer(ip), intent(in)    :: itask,ielem,pgaus,pmate,pnode,lnods(pnode)
  real(rp),    intent(in)    :: gpsha(pnode,pgaus)
  real(rp),    intent(in)    :: gpgdi(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: gpvol(pgaus),gpdet(pgaus)
  real(rp),    intent(in)    :: gppio(ndime,ndime,pgaus)
  real(rp),    intent(in)    :: nfibe(ndime,pgaus)
  real(rp)                   :: gplep(ndime,ndime),gpgre(ndime,ndime)
  real(rp)                   :: gpcrt(ndime,ndime),eval(ndime),evec(ndime,ndime)
  integer(ip)                :: igaus,idime,jdime,kdime,ipoin,inode,ivoig,i,j,k
  integer(ip)                :: nrank,knode,pelty
  real(rp)                   :: sigma(ndime,ndime),gpinv,xmean,gpcrt3(3,3)
  real(rp)                   :: siglo(ndime,ndime),romat(ndime,ndime),theta        ! ( *AQU* )
  real(rp)                   :: eval3(3),evec3(3,3)
  real(rp)                   :: bidon
  real(rp)                   :: gpidg(ndime,ndime,pgaus)
  real(rp)                   :: elcod(ndime,mnode)
  real(rp)                   :: xjaci(9),xjacm(9),xfact
  real(rp)                   :: elunk(nvoig_sld,mnode)
  real(rp)                   :: detjm,gpcar(ndime,pnode)

  !
  ! Extension elements
  !
  if( lelch(ielem) == ELEXT ) then
     knode = 1
  else
     knode = pnode
  end if

  if( itask == 3 ) then

     !-------------------------------------------------------------------
     !
     ! Postprocess
     ! CAUST_SLD: Cauchy Stress tensor
     ! GREEN_SLD: Strain tensor
     ! LEPSI_SLD: Logarithmic strain tensor
     ! EPSEL_SLD: elastic strains
     !
     !-------------------------------------------------------------------

     grlst_sld(ielem,:) = 0.0_rp

     GAUSS_POINTS_LOOP: do igaus = 1,pgaus
        !
        ! Calculate sigma (local, for this G.P): sigma = (1/j)*F*P^T
        !
        sigma = 0.0_rp
        gpinv = 1.0_rp / gpdet(igaus)

        do jdime = 1,ndime
           do idime = 1,ndime
              do kdime = 1,ndime
                 sigma(idime,jdime) = sigma(idime,jdime) &
                      + gpinv * gpgdi(idime,kdime,igaus) * gppio(jdime,kdime,igaus)
              end do
           end do
        end do

        !
        ! Calculate E (Green--Lagrange strain tensor)
        !
        ! E_IJ = 0.5*(F_kI * F_kJ - delta_IJ)
        !
        gpgre = 0.0_rp
        do jdime = 1,ndime
           do idime = 1,ndime
              do kdime = 1,ndime
                 gpgre(idime,jdime) = gpgre(idime,jdime) &
                      + 0.5_rp*( gpgdi(kdime,idime,igaus) * gpgdi(kdime,jdime,igaus) )
              end do
           end do
        end do
        forall(idime=1:ndime) &
             gpgre(idime,idime)=gpgre(idime,idime)-0.5_rp
        !
        ! Save green strain
        !
        epsel_sld(ielem)%a(1:ndime,1:ndime,igaus) = gpgre(1:ndime,1:ndime)
        !
        ! Elemental green lagrange tensor
        !
        if( ndime.eq.3_ip )then !3D CASE
           ! <GGU> Esto no es voigt!!!!
          grlst_sld(ielem,1) = grlst_sld(ielem,1) + gpgre(1,1)*gpvol(igaus)
          grlst_sld(ielem,2) = grlst_sld(ielem,2) + gpgre(1,2)*gpvol(igaus)
          grlst_sld(ielem,3) = grlst_sld(ielem,3) + gpgre(1,3)*gpvol(igaus)
          grlst_sld(ielem,4) = grlst_sld(ielem,4) + gpgre(2,1)*gpvol(igaus)
          grlst_sld(ielem,5) = grlst_sld(ielem,5) + gpgre(2,2)*gpvol(igaus)
          grlst_sld(ielem,6) = grlst_sld(ielem,6) + gpgre(2,3)*gpvol(igaus)
          grlst_sld(ielem,7) = grlst_sld(ielem,7) + gpgre(3,1)*gpvol(igaus)
          grlst_sld(ielem,8) = grlst_sld(ielem,8) + gpgre(3,2)*gpvol(igaus)
          grlst_sld(ielem,9) = grlst_sld(ielem,9) + gpgre(3,3)*gpvol(igaus)
       elseif( ndime.eq.2_ip )then !2D CASE
          ! <GGU> Esto no es voigt!!!!
          grlst_sld(ielem,1) = grlst_sld(ielem,1) + gpgre(1,1)*gpvol(igaus)
          grlst_sld(ielem,2) = grlst_sld(ielem,2) + gpgre(1,2)*gpvol(igaus)
          grlst_sld(ielem,3) = grlst_sld(ielem,3) + gpgre(2,1)*gpvol(igaus)
          grlst_sld(ielem,4) = grlst_sld(ielem,4) + gpgre(2,2)*gpvol(igaus)
        else
          call runend('SLD_BITCUL: no invariant computed for this dimension')
        endif
        !
        ! Calculate LE (Logarithmic strain tensor)
        !
        ! LE_ij = log(V_ij) , where F_iK = V_ij * R_jK <= polar decomposition
        !
        gplep = 0.0_rp
        gpcrt = matmul(gpgdi(:,:,igaus),transpose(gpgdi(:,:,igaus)))

        if( ndime < 3_ip ) then
           gpcrt3 = 0.0_rp
           do i = 1,3
              gpcrt3(i,i) = 1.0_rp
           end do
           do idime = 1,ndime
              do jdime = 1,ndime
                 gpcrt3(idime,jdime) = gpcrt(idime,jdime)
              end do
           end do
           call spcdec(gpcrt3,eval3,evec3,nrank,1_ip,'SLD_BITCUL, NDIME=2, COMPUTING LOG STRESS TENSOR.')

           do k = ndime+1,3
              idime = 0
              do i = 1,3
                 if(abs(evec3(k,i)-1.0_rp) > 1.0e-12_rp) then
                    idime = idime + 1
                    eval(idime) = eval3(i)
                    do jdime = 1,ndime
                       evec(jdime,idime) = evec3(jdime,i)
                    end do
                 endif
              end do
           end do

        else

           call spcdec(gpcrt,eval,evec,nrank,1_ip,'SLD_BITCUL, NDIME=3, COMPUTING LOG STRESS TENSOR.')

        endif

        do idime = 1,ndime
           do jdime = 1,ndime
              do kdime = 1,ndime
                 gplep(idime,jdime) = gplep(idime,jdime) + log(sqrt(eval(kdime)))&
                      *evec(idime,kdime)*evec(jdime,kdime)
              end do
           end do
        end do
        !
        ! Save log strain
        !
        lepse_sld(ielem)%a(1:ndime,1:ndime,igaus) = gplep(1:ndime,1:ndime)

        if( kfl_foten_sld > 0 ) then
           !
           ! Transpose in Voigt notation for global storage
           !
           do inode = 1,knode

              ipoin = lnods(inode)
              xmean = gpsha(inode,igaus) * gpvol(igaus)

#ifdef NO_COLORING
              !$OMP CRITICAL (crit_caust_green_lepsi)
#endif
              do ivoig = 1,nvoig_sld
                 jdime                  = nvgij_sld(ivoig,1)
                 kdime                  = nvgij_sld(ivoig,2)
                 caust_sld(ivoig,ipoin) = caust_sld(ivoig,ipoin) + xmean * sigma(jdime,kdime)
                 green_sld(ivoig,ipoin) = green_sld(ivoig,ipoin) + xmean * gpgre(jdime,kdime)
                 lepsi_sld(ivoig,ipoin) = lepsi_sld(ivoig,ipoin) + xmean * gplep(jdime,kdime)
              end do

#ifdef NO_COLORING
              !$OMP END CRITICAL (crit_caust_green_lepsi)
#endif
           end do
        end if

     end do GAUSS_POINTS_LOOP

  else if( itask == 4 ) then

     !-------------------------------------------------------------------
     !
     ! FIBEG_SLD: Fibers on element
     !
     !-------------------------------------------------------------------

     do igaus = 1,pgaus
        do idime = 1,ndime
           fibeg_sld(ielem)%a(idime,igaus,1) = nfibe(idime,igaus)
        end do
     end do

  else if( itask == 5 ) then

     !-------------------------------------------------------------------
     !
     ! FIBDE_SLD: Fiber's vector output
     !
     !-------------------------------------------------------------------

     if( kfl_fiber_sld > 0 ) then

        if( lelch(ielem) == ELEXT ) then
           knode = 1
        else
           knode = pnode
        end if

        do igaus = 1,pgaus
           do inode = 1,knode

              ipoin = lnods(inode)
              xmean = gpsha(inode,igaus) * gpvol(igaus)

#ifdef NO_COLORING
              !$OMP CRITICAL (crit_fibde)
#endif
              do idime = 1,ndime
                 !fibde_sld(idime,ipoin) = fibde_sld(idime,ipoin) +  xmean * nfibe(idime,igaus)
                 fibde_sld(idime,ipoin) = fibde_sld(idime,ipoin) +  xmean * gpfib(idime,igaus,ielem)
              end do
#ifdef NO_COLORING
              !$OMP END CRITICAL (crit_fibde)
#endif
           end do
        end do
     end if

  else if( itask == 7 ) then

     if( lelch(ielem) == ELEXT ) then
        knode = 1
     else
        knode = pnode
     end if

     do igaus = 1,pgaus

        do inode = 1,knode

           ipoin = lnods(inode)
           xmean = gpsha(inode,igaus) * gpvol(igaus)

#ifdef NO_COLORING
           !$OMP CRITICAL (crit_fpiol)
#endif
           do idime = 1,ndime
              do jdime = 1,ndime
                 nopio_sld(jdime+(idime-1)*ndime,ipoin) = nopio_sld(jdime+(idime-1)*ndime,ipoin) &
                      + xmean * gppio(idime,jdime,igaus)
              end do
           end do
#ifdef NO_COLORING
           !$OMP END CRITICAL (crit_fpiol)
#endif
        end do

     end do

  else if( itask == 10 ) then

     !-------------------------------------------------------------------
     !
     ! Postprocess
     ! CAUSE_SLD: Cauchy Stress tensor
     !
     !-------------------------------------------------------------------

     do igaus = 1,pgaus
        !
        ! Calculate sigma (local, for this G.P): sigma = (1/j)*F*P^T
        !
        sigma = 0.0_rp
        gpinv = 1.0_rp / gpdet(igaus)
        do jdime = 1,ndime
           do idime = 1,ndime
              do kdime = 1,ndime
                 sigma(idime,jdime) = sigma(idime,jdime) &
                      + gpinv * gpgdi(idime,kdime,igaus) * gppio(jdime,kdime,igaus)
              end do
           end do
        end do

        do inode = 1,knode
           ipoin = lnods(inode)
           xmean = gpsha(inode,igaus) * gpvol(igaus)

#ifdef NO_COLORING
           !$OMP CRITICAL (crit_cause)
#endif
           do ivoig = 1,ndime
              jdime = nvgij_sld(ivoig,1)
              kdime = nvgij_sld(ivoig,2)
              cause_sld(ielem)%a(ivoig,igaus,1) = sigma(jdime,kdime)
              !cause_sld(ivoig,ipoin) = cause_sld(ivoig,ipoin) + xmean * sigma(jdime,kdime)
           end do
#ifdef NO_COLORING
           !$OMP END CRITICAL (crit_cause)
#endif
        end do

     end do

  else if( itask == 11 ) then

    !--------------------------------------------------------------------
    !
    ! Postprocess
    ! CAUNP_SLD: Cauchy Stress tensor in the global coordinate system             ( *AQU* )
    !            Extrapolate from gauss points to nodes within each element
    !-------------------------------------------------------------------

    ! Recovered stresses (only for hexahedral elements)
    pelty = ltype(ielem)
    if (pelty==37_ip) then

       do igaus = 1,pgaus
            !
            ! Calculate Cauchy from PK
            !    sigma = (1/j)*F*P^T
            !
            sigma = 0.0_rp
            gpinv = 1.0_rp / gpdet(igaus)
            do jdime = 1,ndime
                do idime = 1,ndime
                    do kdime = 1,ndime
                        sigma(idime,jdime) = sigma(idime,jdime) &
                            + gpinv * gpgdi(idime,kdime,igaus) * gppio(jdime,kdime,igaus)
                    end do
                end do
            end do

            !
            ! Mapping into the Voigt notation
            !
            do inode = 1,knode

                do ivoig = 1,nvoig_sld
                   jdime = nvgij_sld(ivoig,1)
                   kdime = nvgij_sld(ivoig,2)
                   caunp_sld(ielem)%a(ivoig,inode,1) = caunp_sld(ielem)%a(ivoig,inode,1) + gpsha(inode,igaus)*sigma(jdime,kdime)
                end do

            end do
       end do

    end if

  else if( itask == 12 ) then

    !--------------------------------------------------------------------
    !
    ! Postprocess
    ! CAULO_SLD: Cauchy Stress tensor in the local coordinate system             ( *AQU* )
    !
    !-------------------------------------------------------------------

    do igaus = 1,pgaus
        !
        ! Calculate Cauchy from PK
        !    sigma = (1/j)*F*P^T
        !
        sigma = 0.0_rp
        gpinv = 1.0_rp / gpdet(igaus)
        do jdime = 1,ndime
            do idime = 1,ndime
                do kdime = 1,ndime
                    sigma(idime,jdime) = sigma(idime,jdime) &
                        + gpinv * gpgdi(idime,kdime,igaus) * gppio(jdime,kdime,igaus)
                end do
            end do
        end do

        !
        ! Rotate sigma from global to material coordinate system
        !    sigma'_{ij} = G_{ip}*sigma_{pq}*G_{jq] -> sigma' = G*sigma*G'
        !    where the rotation matrix (G_{ij}) is respect to global axis Z = (0, 0 , 1)
        theta = parco_sld(10,pmate)*(pi/180)
        romat = 0.0_rp
        romat(1, 1) =  cos(theta)
        romat(1, 2) =  sin(theta)
        romat(2, 1) = -sin(theta)
        romat(2, 2) =  cos(theta)
        romat(3, 3) =  1.0_rp
        call sld_rotten(1_ip, sigma(:, :), romat(:, :), siglo(:,:), ndime)
        siglo = 0.0_rp

        !
        ! Mapping into the Voigt notation
        !
        do inode = 1,knode
           ipoin = lnods(inode)
           do ivoig = 1,nvoig_sld
               jdime = nvgij_sld(ivoig,1)
               kdime = nvgij_sld(ivoig,2)
               caulo_sld(ielem)%a(ivoig,inode,1) = caulo_sld(ielem)%a(ivoig,inode,1) + gpsha(inode,igaus)*siglo(jdime,kdime)
           end do
        end do

    end do

 else if ( itask == 13 ) then

     !-------------------------------------------------------------------
     !
     ! Postprocess
     ! GREEN_SLD: Infinitesimal strain tensor
     !
     !-------------------------------------------------------------------

     do igaus = 1,pgaus
        !
        ! Calculate Infinitesimal strains
        !
        ! e_IJ = 0.5*(F_IJ * F_JI) - delta_IJ
        !
        gpgre = 0.0_rp
        do i = 1, ndime
          do j = 1, ndime
            gpgre(i,j) = gpgre(i,j) + 0.5_rp*(gpgdi(i,j,igaus) + gpgdi(j,i,igaus))
          end do
        end do

        do i = 1, ndime
          gpgre(i,i) = gpgre(i,i) - 1.0_rp
        end do
        !
        ! Transpose in Voigt notation for global storage
        !
        do inode = 1, knode
          ipoin = lnods(inode)
          do ivoig = 1, nvoig_sld
            j                      = nvgij_sld(ivoig,1)
            k                      = nvgij_sld(ivoig,2)
#ifdef NO_COLORING
            !$OMP ATOMIC
#endif
            green_sld(ivoig,ipoin) = green_sld(ivoig,ipoin) + gpsha(inode,igaus)*gpvol(igaus)*gpgre(j,k)
          end do
        end do
    end do

 else if ( itask == 14 ) then

     !-------------------------------------------------------------------
     !
     ! Postprocess
     ! SIGMS_SLD: Infinitesimal stress tensor
     !
     !-------------------------------------------------------------------

     do igaus = 1,pgaus
        !
        ! Calculate infinitesimal stresses
        !
        !   sigma = inv(F)*P
        !
        call invmtx(gpgdi(:,:,igaus),gpidg(:,:,igaus),bidon,ndime)
        sigma = 0.0_rp
        do jdime = 1,ndime
          do idime = 1,ndime
            do kdime = 1,ndime
              sigma(idime,jdime) = sigma(idime,jdime) +  gpidg(idime,kdime,igaus)  &
              &                    * gppio(jdime,kdime,igaus)
            end do
          end do
        end do
        !
        ! Transpose in Voigt notation for global storage
        !
        do inode = 1, knode
          ipoin = lnods(inode)
          do ivoig = 1, nvoig_sld
            j                      = nvgij_sld(ivoig,1)
            k                      = nvgij_sld(ivoig,2)
#ifdef NO_COLORING
            !$OMP ATOMIC
#endif
            caust_sld(ivoig,ipoin) = caust_sld(ivoig,ipoin) + gpsha(inode,igaus)*gpvol(igaus)*sigma(j,k)
          end do
        end do
    end do

 else if ( itask == 15 ) then

     !-------------------------------------------------------------------
     !
     ! Postprocess
     ! EPSEL_SLD: Elastic strains
     !
     !-------------------------------------------------------------------

     if ( kfl_plast_sld == 1_ip ) then
       !
       ! project from gauss points to nodes
       !
       pelty = ltype(ielem)
       epsee_sld = 0.0_rp

       if( pelty > 0 ) then
        ! pnode = nnode(pelty)
        ! pgaus = ngaus(pelty)

         do inode = 1,pnode
           ipoin = lnods(inode)
           elcod(1:ndime,inode) = coord(1:ndime,ipoin)
           elunk(:,inode) = 0.0_rp
         end do

         gauss_points: do igaus = 1,pgaus
           call elmder(&
                pnode,ndime,elmar(pelty)%deriv(1,1,igaus),&
                elcod,gpcar,detjm,xjacm,xjaci)
           !gpvol = elmar(pelty) % weigp(igaus) * detjm

           do inode = 1, pnode
             xfact        = gpvol(igaus) * elmar(pelty) % shape(inode,igaus)
             do ivoig = 1, nvoig_sld
               elunk(ivoig,inode) = elunk(ivoig,inode) + xfact * svegm_sld(ielem) % a(8+ivoig-1,igaus,1)
             end do
           end do

         end do gauss_points

         do inode = 1,pnode
           ipoin = lnods(inode)
           do ivoig = 1, nvoig_sld
             epsee_sld(ivoig,ipoin) = epsee_sld(ivoig,ipoin) + elunk(ivoig,inode)
           end do
         end do

       end if
    end if
end if

end subroutine sld_bitcul

