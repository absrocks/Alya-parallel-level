!------------------------------------------------------------------------
!> @addtogroup Nastal
!> @{
!> @file    nsa_monyos.f90
!> @date    01/08/2012
!> @author  Mariano Vazquez
!> @brief   Computes the VMS subscale using a non-diagonal Tau matrix
!> @details Computes the VMS subscale xsube(idofn,igaus,1) on the Gauss points\n
!!          as a sum of taure*cosma_nsa+tauim*sinma_nsa\n
!!          over a number of Fourier frequencies, from 0 to mfreq for each space dimension.!!\n
!> @}
!------------------------------------------------------------------------
subroutine nsa_monyos(&
     pelty,pgaus,weigh,ielem,igaus,xresi,xsube,xortp,hleng,xunkn,gunkn,gpres,gvisc,dvelo,velmo,&
     xvelo,xvisc,sound,xlade,xldve,xdtix,taudi,qufac,conme,difme,xtime,mfreq)
  !-----------------------------------------------------------------------
  !***** nastal/nsa_monyos
  ! NAME 
  !    nsa_monyos
  ! DESCRIPTION
  ! USES
  !    nsa_...
  ! USED BY
  !    nsa_...
  !***
  !-----------------------------------------------------------------------
  use      def_master
  use      def_elmtyp
  use      def_domain
  use      def_parame
  use      def_nastal
  implicit none

  integer(ip)  :: &
       ndofc,ndofe,idime,jdime,igaus,jgaus,idofn,jdofn,ifreq,jfreq,kfreq,itisu,ntisu,ielem,&
       pelty,pgaus,frequ(3),mfreq

  real(rp)    :: &
       taudi(ndofn_nsa),xsube(ndofn_nsa,mgaus,3),xortp(ndofn_nsa,mgaus),dsube(ndofn_nsa), xsaux(ndofn_nsa,mgaus),&
       dundt(ndofn_nsa),spare(ndofn_nsa),emean(6+5*ndime+4*ntens),&
       xunkn(ndofn_nsa,mgaus,3), gunkn(ndofn_nsa,ndime,mgaus),&
       gpres(ndime,mgaus,ncomp_nsa),xresi(ndofn_nsa,mgaus),xdesi(ndofn_nsa,mgaus),&
       xvelo(ndime,mgaus,ncomp_nsa),gvisc(ndime,mgaus,ncomp_nsa),xlade(mgaus,ncomp_nsa),xldve(ndime,mgaus,ncomp_nsa),&
       vitot,ditot,sound(mgaus),velmo(mgaus),xvisc(mgaus),dvelo,qufac,tauen,taupa,tauco,hleng(ndime),&
       conme(ndofn_nsa,ndofn_nsa,ndime),difme(ndofn_nsa,ndofn_nsa,ndime,ndime)

  real(rp)     :: &
       tauxi(ndofn_nsa,ndofn_nsa),xdtix(ndofn_nsa,mgaus,2),xauxi,xsute,xunte,xtime(ndofn_nsa,mgaus),&
       vesta(ndofn_nsa,ndime),vista(ndofn_nsa,ndime),rasta(ndofn_nsa),absve(ndofn_nsa),absvi(ndofn_nsa), &
       gseig(ndime),xggde,xggvd,tiint,tiine,tiinc,deltat(ndofn_nsa), &
       taure(ndofn_nsa,ndofn_nsa),tauim(ndofn_nsa,ndofn_nsa), &
       adtre(ndofn_nsa,ndofn_nsa),adtim(ndofn_nsa,ndofn_nsa),detre,detim, &
       detmo,weigh(mgaus)
  real(rp)     :: unacosa, otracosa, supar
  complex(rp)  :: lfo1d_z(ndofn_nsa-1,ndofn_nsa-1),lfi1d_z(ndofn_nsa-1,ndofn_nsa-1),detma_z
  complex(rp)  :: lfour_z(ndofn_nsa  ,ndofn_nsa  ),lfinv_z(ndofn_nsa  ,ndofn_nsa  )


  ndofc= ndime+1
  ndofe= ndime+2

  frequ = 0
  xsute= 0.0_rp
  xunte= 0.0_rp


  do idofn=1,ndofn_nsa
     xsube(idofn,igaus,1) = 0.0_rp
     deltat(idofn) = 1.0_rp / dtinv_nsa        
     do jgaus=1,pgaus
!        xsute = xsube(idofn,jgaus,2) / xdtix(idofn,jgaus,1)       
!        xsute = xsube(idofn,jgaus,2) * dtinv_nsa       
        xdesi(idofn,jgaus) = xsute - xunte + xresi(idofn,jgaus) !+ xortp(idofn,jgaus)
        xsube(idofn,igaus,1) = xsube(idofn,igaus,1) + deltat(idofn) * xdesi(idofn,jgaus) * weigh(jgaus) ! zero frequency term

     end do
     ! used in shock-capturing
     taudi(idofn) = xdtix(idofn,igaus,1)
  end do

!!$if (ielem == 1) then
!!$print*
!!$print*, 'xdesi xsute xunte xresi', xdesi(1,igaus), xsute, xunte, xresi(1,igaus)
!!$end if
!!$if (ielem == 1) then
!!$print*
!!$print*, 'xsube 0', xsube(1,igaus,1)
!!$end if


  !! 1 - Dimentional case
  if(ndime == 1) then
     do ifreq=1,mfreq
        frequ(1)=ifreq
        lfour_z = 0.0_rp
        do idofn=1,ndofn_nsa
           lfour_z(idofn,idofn) = cmplx(1.0_rp/deltat(idofn) , 0.0_rp)
           do jdofn=1,ndofn_nsa        
              do idime=1,ndime
                 lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) + cmplx(0.0_rp , &
                      conme(idofn,jdofn,idime) * real(frequ(idime)) * 2.0_rp * pi / hnatu(pelty))                        
                 do jdime=1,ndime
                    lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) &
                         + cmplx(difme(idofn,jdofn,idime,jdime) &
                         * real(frequ(idime)) * real(frequ(jdime)) * 4.0_rp * pi * pi / hnatu(pelty) / hnatu(pelty), 0.0_rp)
                 end do
              end do
           end do
        end do
        
        call nsa_invmtz(lfour_z,lfinv_z,detma_z,ndofn_nsa)
        
        dsube= 0.0_rp
        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              taure(idofn,jdofn) = real(lfinv_z(idofn,jdofn))
              tauim(idofn,jdofn) = aimag(lfinv_z(idofn,jdofn))
              do jgaus=1,pgaus
                 dsube(idofn) = dsube(idofn) &
                      + (taure(idofn,jdofn) * cosma_nsa(jgaus,igaus, frmax_nsa+1+frequ(1), &
                      frmax_nsa+1+frequ(2), frmax_nsa+1+frequ(3),pelty) &
                      + tauim(idofn,jdofn) * sinma_nsa(jgaus,igaus, frmax_nsa+1+frequ(1), &
                      frmax_nsa+1+frequ(2), frmax_nsa+1+frequ(3),pelty)) &
                      * xdesi(jdofn,jgaus) * weigh(jgaus)
              end do
           end do
           xsube(idofn,igaus,1) = xsube(idofn,igaus,1) + 2.0_rp * dsube(idofn)
        end do
     end do
     
  end if
  
  
  
  !! 2 - Dimentional case
  if(ndime == 2) then
     do ifreq=1,mfreq
        frequ(1)=ifreq
        do jfreq=-mfreq,mfreq
           frequ(2)=jfreq
           lfour_z = 0.0_rp
           do idofn=1,ndofn_nsa
              lfour_z(idofn,idofn) = cmplx(1.0_rp/deltat(idofn) , 0.0_rp)
              do jdofn=1,ndofn_nsa        
                 do idime=1,ndime
                    lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) + cmplx(0.0_rp , &
                         conme(idofn,jdofn,idime) * real(frequ(idime)) * 2.0_rp * pi / hnatu(pelty))                        
                    do jdime=1,ndime
                       lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) &
                            + cmplx(difme(idofn,jdofn,idime,jdime) &
                            * real(frequ(idime)) * real(frequ(jdime)) * 4.0_rp * pi * pi  / hnatu(pelty) / hnatu(pelty), 0.0_rp)
                    end do
                 end do
              end do
           end do
           
           call nsa_invmtz(lfour_z,lfinv_z,detma_z,ndofn_nsa)
           
           dsube= 0.0_rp
           do idofn=1,ndofn_nsa
              do jdofn=1,ndofn_nsa
                 taure(idofn,jdofn) = real(lfinv_z(idofn,jdofn))
                 tauim(idofn,jdofn) = aimag(lfinv_z(idofn,jdofn))
                 do jgaus=1,pgaus
                    dsube(idofn) = dsube(idofn) &
                         + (taure(idofn,jdofn) * cosma_nsa(jgaus,igaus, frmax_nsa+1+frequ(1), &
                         frmax_nsa+1+frequ(2), frmax_nsa+1+frequ(3),pelty) &
                         + tauim(idofn,jdofn) * sinma_nsa(jgaus,igaus, frmax_nsa+1+frequ(1), &
                         frmax_nsa+1+frequ(2), frmax_nsa+1+frequ(3),pelty)) &
                         * xdesi(jdofn,jgaus) * weigh(jgaus)
                 end do
              end do
              xsube(idofn,igaus,1) = xsube(idofn,igaus,1) + 2.0_rp * dsube(idofn)
           end do

!!$if (ielem == 1) then
!!$print*
!!$print*, 'dsube,xsube', frequ(1), frequ(2), dsube(1), xsube(1,igaus,1)
!!$end if

        end do
     end do
     
     
     
     frequ(1)=0.0_rp
     do jfreq=1,mfreq
        frequ(2)=jfreq
        lfour_z = 0.0_rp
        do idofn=1,ndofn_nsa
           lfour_z(idofn,idofn) = cmplx(1.0_rp/deltat(idofn) , 0.0_rp)
           do jdofn=1,ndofn_nsa        
              do idime=1,ndime
                 lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) + cmplx(0.0_rp , &
                      conme(idofn,jdofn,idime) * real(frequ(idime)) * 2.0_rp * pi/ hnatu(pelty))                        
                 do jdime=1,ndime
                    lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) &
                         + cmplx(difme(idofn,jdofn,idime,jdime) &
                         * real(frequ(idime)) * real(frequ(jdime)) * 4.0_rp * pi * pi / hnatu(pelty)/ hnatu(pelty), 0.0_rp)
                 end do
              end do
           end do
        end do
        
        call nsa_invmtz(lfour_z,lfinv_z,detma_z,ndofn_nsa)
        
        dsube= 0.0_rp
        do idofn=1,ndofn_nsa
           do jdofn=1,ndofn_nsa
              taure(idofn,jdofn) = real(lfinv_z(idofn,jdofn))
              tauim(idofn,jdofn) = aimag(lfinv_z(idofn,jdofn))
              do jgaus=1,pgaus
                 dsube(idofn) = dsube(idofn) &
                      + (taure(idofn,jdofn) * cosma_nsa(jgaus,igaus, frmax_nsa+1+frequ(1), &
                      frmax_nsa+1+frequ(2), frmax_nsa+1+frequ(3),pelty) &
                      + tauim(idofn,jdofn) * sinma_nsa(jgaus,igaus, frmax_nsa+1+frequ(1), &
                      frmax_nsa+1+frequ(2), frmax_nsa+1+frequ(3),pelty)) &
                      * xdesi(jdofn,jgaus) * weigh(jgaus)
              end do
           end do
           xsube(idofn,igaus,1) = xsube(idofn,igaus,1) + 2.0_rp * dsube(idofn)
        end do

!!$if (ielem == 1) then
!!$print*
!!$print*, 'dsube,xsube', frequ(1), frequ(2), dsube(1), xsube(1,igaus,1)
!!$end if

     end do
     
  end if

  supar= 4.0_rp
  if (pelty == TRI03) supar= 2.0_rp
  do idofn=1,ndofn_nsa
     xsube(idofn,igaus,1) = xsube(idofn,igaus,1) / supar
  end do

end subroutine nsa_monyos
