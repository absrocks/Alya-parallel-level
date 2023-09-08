subroutine nsa_raydam(lnods,pnode,pgaus,ivert,gpsha,elcod,elumo,elden,taray)
!-----------------------------------------------------------------------
!
! Rayleigh dumping 
!
!
!-----------------------------------------------------------------------------------------------------------
  use      def_kintyp
  use      def_parame, only : pi
  use      def_master
  use      def_domain, only : coord,ndime,mnode,mgaus
  use      def_nastal
  implicit none
  integer(ip)    :: pnode,pgaus,lnods(pnode)
  integer(ip)    :: inode,jnode,igaus,idime,ipoin,ivert
  real(rp)       :: &
       gpsha(pnode,pgaus),elcod(ndime,mnode),elumo(ndime,mnode),elden(mnode),taray(ndofn_nsa,ndofn_nsa,mgaus)
  real(rp)       :: &
       taden(ndime),tavel(ndime),tarau(ndofn_nsa,ndofn_nsa),aauxi,bauxi,cauxi,dauxi,umodu,xshap,coaux

  taray = 0.0_rp

  do inode=1,pnode     
     ipoin= lnods(inode)

     !
     ! Compute taden and tavel, rayleigh taus for pi and velocity
     !
     
     do idime=1,ndime
        tavel(idime)= 0.0_rp
        taden(idime)= 0.0_rp
        if (kranr_nsa(idime,1) == 1) then                ! lower

           ! to be done

        else if (kranr_nsa(idime,2) == 1) then           ! greater

           if (kranr_nsa(idime,4) == 1) then             !   gravitational waves
              if (coord(ipoin,idime) > vranr_nsa(idime,2)) then
                 coaux= (coord(ipoin,idime) - vranr_nsa(idime,2)) / &
                      (vranr_nsa(idime,3)- vranr_nsa(idime,2))  
                 if (coaux <= 0.5_rp) then
                    tavel(idime)= frayl_nsa(idime,1) * (1.0_rp - cos(coaux * pi) )
                 else
                    tavel(idime)= frayl_nsa(idime,1) * (1.0_rp + (coaux - 0-5)*pi )     ! esto de multiplicar por pi esta bien?????
                 end if
              end if
           end if

           if (kranr_nsa(idime,3) == 1) then        !   acoustic waves
              coaux= (vranr_nsa(idime,3) - coord(ipoin,idime)) / frayl_nsa(idime,2)                   
              coaux= tanh(coaux)
              tavel(idime) = tavel(idime) + (1.0_rp - coaux) / coaux 
           end if

        else if (kranr_nsa(idime,2) == 2) then           ! limits

           if (coord(ipoin,idime) > vranr_nsa(idime,3) ) then
              coaux= (vranr_nsa(idime,2) - coord(ipoin,idime)) / frayl_nsa(idime,2) 
              coaux= tanh(coaux)
           else
              coaux= (coord(ipoin,idime) - vranr_nsa(idime,1)) / frayl_nsa(idime,2) 
              coaux= tanh(coaux)
           end if

           tavel(idime) = tavel(idime) + (1.0_rp - coaux) / coaux 
           
           taden(idime) = tavel(idime)

        end if
     end do
     
     !
     ! Compute tarau, auxiliar tau for conservartive
     !     
     umodu= elumo(1,inode)*elumo(1,inode) + elumo(2,inode)*elumo(2,inode) 
     if (ndime == 3) umodu= umodu + elumo(ndime,inode)*elumo(ndime,inode)
     
     cauxi= umodu - 2.0_rp*(elumo(ivert,inode)*elumo(ivert,inode)+        grnor_nsa*elden(inode)*elden(inode)*elcod(ivert,inode))
     dauxi=                 umodu                                - 2.0_rp*grnor_nsa*elden(inode)*elden(inode)*elcod(ivert,inode)
     
     aauxi= elumo(1,inode)*elumo(1,inode) - 0.5_rp * cauxi     
     bauxi= aauxi * elden(inode)
     
     
     ! ver lo del ivert, que debe haber un tarau(ivert,...)
     
!     tarau(1,1) = (elden(inode) * elumo(1,inode) * ( taden(1) - tavel(1) ) )/ aauxi
!     tarau(1,2) = (elden(inode) * elumo(2,inode) * ( taden(2) - tavel ) )/ aauxi
     
!     tarau(ndime+1,ndime+1) = ( tavel * umodu - 0.5_rp * taden * dauxi ) / aauxi
     
     !     tarau(2,1) = ...
     !     tarau(2,2) = ...
     !     tarau(2,3) = ...
     
     
     do igaus=1,pgaus
        xshap= gpsha(inode,igaus)
        taray(      1,      1,igaus) = taray(      1,      1,igaus) + xshap*tarau(      1,      1)
        taray(      2,      2,igaus) = taray(      2,      2,igaus) + xshap*tarau(      2,      2)
        taray(ndime+1,ndime+1,igaus) = taray(ndime+1,ndime+1,igaus) + xshap*tarau(ndime+1,ndime+1)
        !     taray = ...
     end do
     
     if (ndime == 3) then
!        tarau(1,3) = (elden(inode) * elumo(3,inode) * ( taden - tavel ) )/ aauxi
        do igaus=1,pgaus
           xshap= gpsha(inode,igaus)
           taray(      1,      3,igaus) = taray(      1,      3,igaus) + xshap*tarau(      1,      3)
           !     taray = ...
        end do
     end if
     
     
  end do



end subroutine nsa_raydam
