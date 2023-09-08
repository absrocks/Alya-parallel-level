subroutine bsymai(itask,nbnodes,nbvar,an,ja,ia,xx,yy)
  !-----------------------------------------------------------------------
  ! Objective:  Multiply a symmetric matrix stored in BCSR by a vector.
  !                                yy = [A] xx
  !
  !             Only the LOWER triangle is stored. 
  !             JA elements are stored in increasing column order
  !-----------------------------------------------------------------------
  use def_kintyp, only             :  ip,rp
  use def_master, only             :  INOTMASTER,IPARALL,NPOIN_TYPE,npoi1
  implicit none
  integer(ip), intent(in)          :: itask,nbnodes,nbvar
  real(rp),    intent(in)          :: an(nbvar,nbvar,*)
  integer(ip), intent(in)          :: ja(*),ia(*)
  real(rp),    intent(in)          :: xx(nbvar,nbnodes)
  real(rp),    intent(out), target :: yy(nbvar,nbnodes)
  integer(ip)                      :: ii,jj,kk,ll,rc,dia_i
  real(rp)                         :: raux,raux2

  call runend('BSYMAI: ASYNCHRONOUS COMMUNCATIONS IMPOSSIBLE WITH UNSYMMETRIC ASSEMBLY')

  if( INOTMASTER ) then

     !-------------------------------------------------------------------
     !
     ! Boundary nodes
     !
     !-------------------------------------------------------------------

     if( nbvar > 1 ) then

        do ii = npoi1+1,nbnodes
           do kk= 1, nbvar
              yy(kk,ii) = 0.0_rp
           end do
        end do

        do ii = npoi1+1,nbnodes
           !
           ! Lower triangle + Upper triangle (blocks are transposed)
           !
           do jj = ia(ii), ia(ii+1)-2
              rc = ja(jj)

              do ll= 1, nbvar
                 raux  = 0.0_rp
                 raux2 = xx(ll,rc)
                 do kk= 1, nbvar
                    yy(kk,ii) = yy(kk,ii) + an(kk,ll,jj) * raux2
                    raux = raux + an(kk,ll,jj) * xx(kk,ii)
                 end do
                 yy(ll,rc) = yy(ll,rc) + raux
              end do
           end do
           !
           ! Diagonal block
           !
           dia_i = ia(ii+1)-1
           do ll= 1, nbvar
              raux = xx(ll,ii)
              do kk= ll, nbvar
                 yy(kk,ii) = yy(kk,ii) + an(kk,ll,dia_i) * raux
              end do

              raux = 0.0_rp
              do kk= ll+1, nbvar
                 raux = raux + an(kk,ll,dia_i) * xx(kk,ii)
              end do

              yy(ll,ii) = yy(ll,ii) + raux
           end do
        end do

     else if( nbvar == 1 ) then
        !
        ! Lower triangle + diag
        !
        do ii= npoi1+1,nbnodes
           raux  = 0.0_rp
           raux2 = xx(1,ii)
           do jj= ia(ii), ia(ii+1)-2
              raux         = raux + an(1,1,jj) * xx(1,ja(jj))
              yy(1,ja(jj)) = yy(1,ja(jj)) + an(1,1,jj) * raux2
           enddo
           kk       = ia(ii+1)-1
           yy(1,ii) = raux + an(1,1,kk) * xx(1,ja(kk))
        enddo

     else
        call runend('bcsrmvsym_d: Wrong NBVAR parameter')
     end if

     !-------------------------------------------------------------------
     !
     ! Modify YY due do periodicity and Parall service
     !
     !-------------------------------------------------------------------

     if( itask == 1 .and. IPARALL ) then
        call pararr('SLA',NPOIN_TYPE,nbnodes*nbvar,yy)
     end if

     !-------------------------------------------------------------------
     !
     ! Interior nodes
     !
     !-------------------------------------------------------------------

     if( nbvar > 1 ) then

        do ii= 1,npoi1
           do kk= 1, nbvar
              yy(kk,ii) = 0.0_rp
           end do
        end do

        do ii= 1,npoi1
           !
           ! Lower triangle + Upper triangle (blocks are transposed)
           !
           do jj= ia(ii), ia(ii+1)-2
              rc = ja(jj)

              do ll= 1, nbvar
                 raux  = 0.0_rp
                 raux2 = xx(ll,rc)
                 do kk= 1, nbvar
                    yy(kk,ii) = yy(kk,ii) + an(kk,ll,jj) * raux2
                    raux = raux + an(kk,ll,jj) * xx(kk,ii)
                 end do
                 yy(ll,rc) = yy(ll,rc) + raux
              end do
           end do
           !
           ! Diagonal block
           !
           dia_i = ia(ii+1)-1
           do ll= 1, nbvar
              raux = xx(ll,ii)
              do kk= ll, nbvar
                 yy(kk,ii) = yy(kk,ii) + an(kk,ll,dia_i) * raux
              end do

              raux = 0.0_rp
              do kk= ll+1, nbvar
                 raux = raux + an(kk,ll,dia_i) * xx(kk,ii)
              end do

              yy(ll,ii) = yy(ll,ii) + raux
           end do
        end do

     else if( nbvar == 1 ) then
        !
        ! Lower triangle + diag
        !
        do ii= 1,npoi1
           raux  = 0.0_rp
           raux2 = xx(1,ii)
           do jj= ia(ii), ia(ii+1)-2
              raux         = raux + an(1,1,jj) * xx(1,ja(jj))
              yy(1,ja(jj)) = yy(1,ja(jj)) + an(1,1,jj) * raux2
           enddo
           kk       = ia(ii+1)-1
           yy(1,ii) = raux + an(1,1,kk) * xx(1,ja(kk))
        enddo

     else
        call runend('bcsrmvsym_d: Wrong NBVAR parameter')
     end if

     !-------------------------------------------------------------------
     !
     ! Wait all and sum up contributions of boundary nodes
     !
     !-------------------------------------------------------------------

     if( itask == 1 .and. IPARALL ) then
        call pararr('SLA',NPOIN_TYPE,nbnodes*nbvar,yy)
     end if

  end if

end subroutine bsymai
