subroutine nsa_monyos_integracion(&
     pelty,pgaus,weigh,ielem,igaus,xresi,xsube,hleng,xunkn,gunkn,gpres,gvisc,dvelo,velmo,&
     xvelo,xvisc,sound,xlade,xldve,xdtix,taudi,qufac,conme,difme,xtime,nodco)
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
  use      def_domain
  use      def_parame
  use      def_nastal
  implicit none

  integer(ip)  :: &
       ndofc,ndofe,idime,jdime,igaus,jgaus,idofn,jdofn,ifreq,jfreq,kfreq,itisu,ntisu,ielem,&
       pelty,pgaus,frequ(3)

  real(rp)    :: &
       taudi(ndofn_nsa),xsube(ndofn_nsa,mgaus,3),dsube(ndofn_nsa), xsaux(ndofn_nsa,mgaus),&
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
  real(rp)     :: unacosa, otracosa
  complex(rp)  :: lfo1d_z(ndofn_nsa-1,ndofn_nsa-1),lfi1d_z(ndofn_nsa-1,ndofn_nsa-1),detma_z
  complex(rp)  :: lfour_z(ndofn_nsa,ndofn_nsa),lfinv_z(ndofn_nsa,ndofn_nsa)

  integer(ip) :: imade,idegr,jdegr
  real(rp)    :: intco(ndofn_nsa), intsi(ndofn_nsa), nodco(5,5,5,ndofn_nsa)
  complex(rp) :: integ(ndofn_nsa)



  ndofc= ndime+1
  ndofe= ndime+2

  frequ = 0
  xsute= 0.0_rp
  xunte= 0.0_rp


  do idofn=1,ndofn_nsa
     xsube(idofn,igaus,1) = 0.0_rp
     deltat(idofn) = 1.0_rp / dtinv        
     do jgaus=1,pgaus
!!        xsute = xsube(idofn,jgaus,2) / xdtix(idofn,jgaus,2)       
        xdesi(idofn,jgaus) = xsute + xresi(idofn,jgaus) - xunte
        xsube(idofn,igaus,1) = xsube(idofn,igaus,1) + deltat(idofn) * xdesi(idofn,jgaus) * weigh(jgaus)
     end do
     ! used in shock-capturing
     taudi(idofn) = xdtix(idofn,igaus,1)
  end do


!if (ielem .eq. 1 .and. igaus .eq. 1) then
!write(6,*)
!write(6,*) 'xsube 0:', xsube(1:4,igaus,1)
!end if



  !! 1 - Dimentional case
  if(ndime == 1) then
     do ifreq=1,mfreq_nsa
        frequ(1)=ifreq

     end do
     
  end if
  
!write(6,*)
!write(6,*) 'pelty:', pelty
!stop

  
  !! 2 - Dimentional case
  if(ndime == 2) then
     do ifreq=1,mfreq_nsa
        frequ(1)=ifreq
        do jfreq=-mfreq_nsa,mfreq_nsa
           frequ(2)=jfreq

!if (ielem .eq. 1 .and. igaus .eq. 1) then
!write(6,*)
!write(6,*) 'frequ(1) :', frequ(1)
!write(6,*) 'frequ(2) :', frequ(2)
!write(6,*) 'geoco_nsa :', geoco_nsa(1:9,igaus,frequ(1)+mfreq_nsa+1,frequ(2)+mfreq_nsa+1,frequ(3)+mfreq_nsa+1,pelty)
!end if


           ! compute lfinv_z (inverse of the Fourier transform subscale operator)
           ! compute integ (integral terms in the subscale expression)
           lfour_z = 0.0_rp
           do idofn=1,ndofn_nsa
              lfour_z(idofn,idofn) = cmplx(1.0_rp/deltat(idofn) , 0.0_rp)
              do jdofn=1,ndofn_nsa        
                 do idime=1,ndime
                    lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) + cmplx(0.0_rp , &
                         conme(idofn,jdofn,idime) * real(frequ(idime)) * pi)                        
                    do jdime=1,ndime
                       lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) &
                            + cmplx(difme(idofn,jdofn,idime,jdime) &
                            * real(frequ(idime)) * real(frequ(jdime)) * pi * pi , 0.0_rp)
                    end do
                 end do
              end do

              integ(idofn) = 0.0_rp
              do idegr=1,3 ! 5 for the 3D case
!                 do jdegr=1,3 ! 5 for the 3D case
!                    imade = (idegr-1) * 3.0_rp + jdegr ! 5.0_rp for the 3D case
!                    integ(idofn) = integ(idofn) + nodco(idegr,jdegr,1,idofn) * &
!                         geoco_nsa(imade,igaus,mfreq_nsa+1+frequ(1),mfreq_nsa+1+frequ(2),mfreq_nsa+1+frequ(3),pelty)
!                 end do
              end do
              intco(idofn) = real(integ(idofn))
              intsi(idofn) = aimag(integ(idofn))
           end do


!if (ielem .eq. 1 .and. igaus .eq. 1) then
!write(6,*)
!write(6,*) 'integ :', integ(1:4)
!end if



           
           call nsa_invmtz(lfour_z,lfinv_z,detma_z,ndofn_nsa)
           
           ! Compute the subscale
           dsube = 0.0_rp !!!!!!!!!! OJOOOO : ho podria posar després del do idofn=..., com dsube(idofn)=0
           do idofn=1,ndofn_nsa
              do jdofn=1,ndofn_nsa
                 taure(idofn,jdofn) = real(lfinv_z(idofn,jdofn))
                 tauim(idofn,jdofn) = aimag(lfinv_z(idofn,jdofn))

                 dsube(idofn) = dsube(idofn) + taure(idofn,jdofn) * intco(jdofn) + tauim(idofn,jdofn) * intsi(jdofn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              do jgaus=1,pgaus
!                 dsube(idofn) = dsube(idofn) &
!                      + (taure(idofn,jdofn) * cosma_nsa(jgaus,igaus, mfreq_nsa+1+frequ(1), &
!                      mfreq_nsa+1+frequ(2), mfreq_nsa+1+frequ(3),pelty) &
!                      + tauim(idofn,jdofn) * sinma_nsa(jgaus,igaus, mfreq_nsa+1+frequ(1), &
!                      mfreq_nsa+1+frequ(2), mfreq_nsa+1+frequ(3),pelty)) &
!                      * weigh(jgaus) * xsube(jdofn,jgaus,2) / xdtix(jdofn,jgaus,2)
!              end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!...................................
!                 do jgaus=1,pgaus
!                    dsube(idofn) = dsube(idofn) &
!                         + (taure(idofn,jdofn) * cosma_nsa(jgaus,igaus, mfreq_nsa+1+frequ(1), &
!                         mfreq_nsa+1+frequ(2), mfreq_nsa+1+frequ(3),pelty) &
!                         + tauim(idofn,jdofn) * sinma_nsa(jgaus,igaus, mfreq_nsa+1+frequ(1), &
!                         mfreq_nsa+1+frequ(2), mfreq_nsa+1+frequ(3),pelty)) &
!                         * xdesi(jdofn,jgaus) * weigh(jgaus)
!                 end do
!....................................
              end do
              xsube(idofn,igaus,1) = xsube(idofn,igaus,1) + 2.0_rp * dsube(idofn)
           end do


!if (ielem .eq. 1 .and. igaus .eq. 1) then
!write(6,*)
!write(6,*) 'dsube :', dsube(1:4)
!end if


!...................................
        end do
     end do
     
     
     
     frequ(1) = 0
     do jfreq=1,mfreq_nsa
        frequ(2)=jfreq

!if (ielem .eq. 1 .and. igaus .eq. 1) then
!write(6,*)
!write(6,*) 'frequ(1) :', frequ(1)
!write(6,*) 'frequ(2) :', frequ(2)
!write(6,*) 'geoco_nsa :', geoco_nsa(1:9,igaus,frequ(1)+mfreq_nsa+1,frequ(2)+mfreq_nsa+1,frequ(3)+mfreq_nsa+1,pelty)
!end if

        ! compute lfinv_z (inverse of the Fourier transform subscale operator)
           ! compute integ (integral terms in the subscale expression)
            lfour_z = 0.0_rp
            do idofn=1,ndofn_nsa
              lfour_z(idofn,idofn) = cmplx(1.0_rp/deltat(idofn) , 0.0_rp)
              do jdofn=1,ndofn_nsa        
                 do idime=1,ndime
                    lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) + cmplx(0.0_rp , &
                         conme(idofn,jdofn,idime) * real(frequ(idime)) * pi)                        
                    do jdime=1,ndime
                       lfour_z(idofn,jdofn) = lfour_z(idofn,jdofn) &
                            + cmplx(difme(idofn,jdofn,idime,jdime) &
                            * real(frequ(idime)) * real(frequ(jdime)) * pi * pi , 0.0_rp)
                    end do
                 end do
              end do

              integ(idofn) = 0.0_rp
              do idegr=1,3 ! 5 for the 3D case
                 do jdegr=1,3 ! 5 for the 3D case
!                    imade = (idegr-1) * 3.0_rp + jdegr ! 5.0_rp for the 3D case
!                    integ(idofn) = integ(idofn) + nodco(idegr,jdegr,1,idofn) * &
!                         geoco_nsa(imade,igaus,mfreq_nsa+1+frequ(1),mfreq_nsa+1+frequ(2),mfreq_nsa+1+frequ(3),pelty)
                 end do
              end do
              intco(idofn) = real(integ(idofn))
              intsi(idofn) = aimag(integ(idofn))

           end do


!if (ielem .eq. 1 .and. igaus .eq. 1) then
!write(6,*)
!write(6,*) 'integ :', integ(1:4)
!end if


           
           call nsa_invmtz(lfour_z,lfinv_z,detma_z,ndofn_nsa)
           
           ! Compute the subscale
           dsube= 0.0_rp
           do idofn=1,ndofn_nsa
              do jdofn=1,ndofn_nsa
                 taure(idofn,jdofn) = real(lfinv_z(idofn,jdofn))
                 tauim(idofn,jdofn) = aimag(lfinv_z(idofn,jdofn))

                 dsube(idofn) = dsube(idofn) + taure(idofn,jdofn) * intco(jdofn) + tauim(idofn,jdofn) * intsi(jdofn)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              do jgaus=1,pgaus
!                 dsube(idofn) = dsube(idofn) &
!                      + (taure(idofn,jdofn) * cosma_nsa(jgaus,igaus, mfreq_nsa+1+frequ(1), &
!                      mfreq_nsa+1+frequ(2), mfreq_nsa+1+frequ(3),pelty) &
!                      + tauim(idofn,jdofn) * sinma_nsa(jgaus,igaus, mfreq_nsa+1+frequ(1), &
!                      mfreq_nsa+1+frequ(2), mfreq_nsa+1+frequ(3),pelty)) &
!                      * weigh(jgaus) * xsube(jdofn,jgaus,2) / xdtix(jdofn,jgaus,2)
!              end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!...................................
!                 do jgaus=1,pgaus
!                    dsube(idofn) = dsube(idofn) &
!                         + (taure(idofn,jdofn) * cosma_nsa(jgaus,igaus, mfreq_nsa+1+frequ(1), &
!                         mfreq_nsa+1+frequ(2), mfreq_nsa+1+frequ(3),pelty) &
!                         + tauim(idofn,jdofn) * sinma_nsa(jgaus,igaus, mfreq_nsa+1+frequ(1), &
!                         mfreq_nsa+1+frequ(2), mfreq_nsa+1+frequ(3),pelty)) &
!                         * xdesi(jdofn,jgaus) * weigh(jgaus)
!                 end do
!....................................
              end do
              xsube(idofn,igaus,1) = xsube(idofn,igaus,1) + 2.0_rp * dsube(idofn)
           end do


!if (ielem .eq. 1 .and. igaus .eq. 1) then
!write(6,*)
!write(6,*) 'dsube :', dsube(1:4)
!end if


!...................................................

     end do
     
  end if


!if (ielem .eq. 1 .and. igaus .eq. 1) then
!write(6,*)
!write(6,*) 'xsube :', xsube(1:4,igaus,1)
!end if




  
  do idofn=1,ndofn_nsa
     xsube(idofn,igaus,1) = xsube(idofn,igaus,1) / 4.0_rp
  end do

if (ielem .eq. 1 .and. igaus .eq. 1) then
write(6,*)
write(6,*) 'xsube FINAL:', xsube(1:4,igaus,1)
end if
!stop

end subroutine nsa_monyos_integracion
